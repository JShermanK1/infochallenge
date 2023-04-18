use noodles::{
    bcf::{self, header::StringMaps},
    vcf::Header
};
use polars::prelude::*;
use std::{
    error::Error,
    fs::File,
    iter::repeat,
    path::PathBuf,
};
use itertools::Itertools;
use noodles::csi;
use rayon::prelude::*;
use clap::{
    Command, 
    arg, 
    value_parser
};
use sync_file::SyncFile;



fn cli() -> Command {
    Command::new("infochallenge")
            .args(&[
                arg!(
                    -i --in <FILE> "Path to bcf file with index"
                ).required(true)
                 .value_parser(value_parser!(PathBuf)),
            ])
}


fn main() -> Result<(), Box<dyn Error>> {
    let matches = cli().get_matches();
    let in_path = matches.get_one::<PathBuf>("in").expect("error parsing input path");
    let mut index_path = in_path.clone();
    index_path.set_extension("bcf.csi");


    let mut index_r = File::open(index_path).map(csi::Reader::new)?;
    let index = index_r.read_index()?;

    let f = SyncFile::open(in_path)?;

    let num_samples;
    let stringmaps: StringMaps;
    let header: Header;

    {
        let mut bcf_r = bcf::Reader::new(f.clone());

        bcf_r.read_file_format()?;
        let r_header = bcf_r.read_header()?;
        header = r_header.parse()?;
        stringmaps = r_header.parse()?;

        let mut records = bcf_r.records().peekable();

        let first_record = records.peek().unwrap().clone();
        num_samples = first_record.as_ref().unwrap().genotypes().len();
    }

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

        let region = format!("{chrom}").parse().expect("failed to parse region");

        let mut bcf_r = bcf::Reader::new(f.clone());

        bcf_r.read_file_format().expect("failed to read format");
        let _r_header = bcf_r.read_header().expect("failed to read header");

        let records = bcf_r.query(stringmaps.contigs(), &index, &region)
                           .expect("failed to query index");
        

        
        //preallocate vectors to store data.
        let mut samples: Vec<Vec<Option<i8>>> = (0..num_samples).map(|_| {
            Vec::with_capacity(35e5 as usize) // a random guess from looking at the data
        }).collect();

        records.for_each(|record| {
            let record = record.expect("failed to read record");
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
        });


        let data = samples.into_iter()
                .map(|sample| {
                    let chunked = sample.into_iter().collect();
                    chunked
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
                        .for_each(|(avg, (i, j))| {
                            //store values in the matrix in both directions
                            matrix[j][i] = avg;
                            matrix[i][j] = avg;

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
