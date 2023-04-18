use futures::{TryStreamExt, future};
use noodles::{
    bcf::{self, header::StringMaps},
    vcf::Header, csi::Index
};
use polars::prelude::*;
use std::{
    error::Error,
    iter::repeat,
    path::PathBuf,
};
use tokio::{
    fs::File,
};
use itertools::Itertools;
use noodles::csi;
use rayon::prelude::*;
use clap::{
    Command, 
    arg, 
    value_parser
};


fn cli() -> Command {
    Command::new("infochallenge")
            .args(&[
                arg!(
                    -i --in <FILE> "Path to bcf file with index"
                ).required(true)
                 .value_parser(value_parser!(PathBuf)),
            ])
}

async fn parse_header(in_path: &PathBuf) -> Result<(usize, StringMaps, Header), Box<dyn Error>> {

    let mut bcf_r = File::open(in_path).await.map(bcf::AsyncReader::new)?;

    bcf_r.read_file_format().await?;
    let r_header = bcf_r.read_header().await?;
    let header = r_header.parse()?;
    let stringmaps = r_header.parse()?;

    let mut first_record = bcf::Record::default();
    bcf_r.read_record(&mut first_record).await.unwrap();
    let num_samples = first_record.genotypes().len();

    Ok((num_samples, stringmaps, header))
}

async fn collect_snps(chrom: &str, stringmaps: &StringMaps, num_samples: usize, index: &Index) -> Result<Vec<Vec<Option<i8>>>, Box<dyn Error>> {
    let region = format!("{chrom}").parse().expect("failed to parse region");

    let mut bcf_r = File::open("/gpfs0/scratch/mvc002/info_challenge.bcf").await.map(bcf::AsyncReader::new).expect("failed to open file");

    bcf_r.read_file_format().await.expect("failed to read format");
    let _r_header = bcf_r.read_header().await.expect("failed to read header");

    let records = bcf_r.query(stringmaps.contigs(), index, &region)
                        .expect("failed to query index");
    

    
    //preallocate vectors to store data.
    let mut samples: Vec<Vec<Option<i8>>> = (0..num_samples).map(|_| {
        Vec::with_capacity(35e5 as usize) // a random guess from looking at the data
    }).collect();

    let fut = records.try_for_each(|record| {
        //let record = record.expect("failed to read record");
        let gt = record.genotypes().as_ref();
        //3 header values then genotypes are pairs

        //get only the GT array, drop other data
        let record_vals = gt[3..(num_samples * 2 + 3)].chunks_exact(2)
                                .map(|chunk| (chunk[0] + chunk[1]) as i8) //combine both alelles into one number
                                .enumerate();
        //optimizes out bounds checking
        assert!(samples.len() == record_vals.len());
        record_vals.for_each(|(i, read_value)| {

                //map allels to helpfull numbers otherwise replace with None
                let value = match read_value { 
                    0 => None,
                    4 => Some(0_i8),
                    6 => Some(1_i8),
                    8 => Some(2_i8),
                    _ => None
                };
                samples[i].push(value);
            });

        future::ready(Ok(()))
    });
    
    fut.await.expect("failed to read all records");

    Ok(samples)
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    let matches = cli().get_matches();
    let in_path = matches.get_one::<PathBuf>("in").expect("error parsing input path");
    let mut index_path = in_path.clone();
    index_path.set_extension("bcf.csi");

    let rt_handle = tokio::runtime::Handle::current();


    let mut index_r = std::fs::File::open(index_path).map(csi::Reader::new)?;
    let index = index_r.read_index()?;

    let num_samples;
    let stringmaps: StringMaps;
    let header: Header;

    let t = parse_header(in_path);

    (num_samples, stringmaps, header) = t.await?;

    let names_top = "\t".to_string() + header.sample_names()
                                            .iter()
                                            .join("\t")
                                            .as_str();

    println!("{names_top}");

    let data = (0..).map_while(|i|{

        stringmaps.contigs().get_index(i)

                        }).collect::<Vec<&str>>()
                        .into_par_iter()
                        .map(|chrom| {

        let fut_samples = collect_snps(chrom, &stringmaps, num_samples, &index);
        let samples;
        {
            samples = rt_handle.block_on(fut_samples).expect("falied to collect snps");
        }

        let data = samples.into_par_iter()
                .map(|sample| {
                    let sample = sample.into_iter();
                    ChunkedArray::<Int8Type>::from_iter(sample)
                })
                .collect::<Vec<ChunkedArray<Int8Type>>>();
        
        data
    }).reduce(
        || (0..num_samples).map(|_| ChunkedArray::<Int8Type>::default())
                            .collect::<Vec<ChunkedArray<Int8Type>>>(),
        |a_s, b_s| {
            a_s.into_iter()
             .zip(b_s)
             .map(|(mut a, b)| {
                a.append(&b);
                a
             }).collect::<Vec<ChunkedArray<Int8Type>>>()
        }
    );             

    let data = Arc::new(data);

    //preinitiallize matrix with 0.0 
    let mut matrix = (0_usize..num_samples).map(|_| repeat(0.0).take(num_samples).collect())
                                 .collect::<Vec<Vec<f64>>>();
    
    //returns all combinations of all numbers in the sequence as tuples
    (0_usize..num_samples).tuple_combinations::<(usize, usize)>()
                        .collect::<Vec<(usize, usize)>>()
                        .into_par_iter()
                        .map(|(i, j)| {

                            let mask = &data[i].is_not_null() & &data[j].is_not_null();
                            let data1 = &data[i].filter(&mask).unwrap();
                            let data2 = &data[j].filter(&mask).unwrap();
                            let num_compares = (data1.len() * 2) as f64;

                            let diff: f64 = (data1 - data2).abs().cast(&DataType::Float64).unwrap().sum().unwrap();
                            diff / num_compares

                        }).collect::<Vec<f64>>()
                        .into_iter()
                        .zip((0_usize..num_samples).tuple_combinations::<(usize, usize)>())
                        .for_each(|(distance, (i, j))| {
                            //store values in the matrix in both directions
                            matrix[j][i] = distance;
                            matrix[i][j] = distance;
                        });

    matrix.iter()
          .zip(header.sample_names().iter())
          .for_each(|(line, samp_name)| {
            let line = line.iter()
                           .join("\t");
            println!("{samp_name}\t{line}");
        });

    Ok(())
}
