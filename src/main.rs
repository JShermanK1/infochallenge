use noodles::bcf::header::{StringMaps};
use noodles::bcf::{self};
use noodles::vcf::Header;
use polars::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::iter::repeat;
use itertools::Itertools;
use noodles::csi::{self};
use rayon::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {
    
    let f = File::open("/gpfs0/scratch/mvc002/info_challenge.bcf")?;
    let mut bcf_r = bcf::Reader::new(f);

    bcf_r.read_file_format()?;
    let r_header = bcf_r.read_header()?;
    let header: Header = r_header.parse()?;
    let _stringmaps: StringMaps = r_header.parse()?;

    let mut records = bcf_r.records().peekable();
   

    let mut index_r = File::open("/gpfs0/scratch/mvc002/info_challenge.bcf.csi").map(csi::Reader::new)?;
    let index = index_r.read_index()?;
    let num_records = index.reference_sequences().par_iter().map(|refseq| {
        let meta_o =  refseq.metadata();
        if let Some(meta) = meta_o {
            meta.mapped_record_count()
        } else {
            0
        }
    }).reduce(|| 0,|a,b| a + b );

    println!("{:?}", num_records);


    
    let mut samples: Vec<Vec<Option<i8>>> = (0..29).map(|_| {
        Vec::with_capacity(num_records as usize)
    }).collect();

    let first_record = records.peek().unwrap().clone();
    let num_samples = first_record.as_ref().unwrap().genotypes().len();
    

    records.for_each(|record| {
        let record = record.expect("failed to read record");
        let gt = record.genotypes().as_ref();
        //3 header values then genotypes are pairs

        //get only the GT array, drop other data
        let record_vals = gt[3..(num_samples * 2 + 3)].chunks_exact(2)
                                .map(|chunk| (chunk[0] + chunk[1]) as i8) //combine both alelles into one number
                                .enumerate();
        assert!(samples.len() == record_vals.len());
        record_vals.for_each(|(i, read_value)| {

                let value = match read_value { //
                    0 => None,
                    4 => Some(0_i8),
                    6 => Some(1_i8),
                    8 => Some(2_i8),
                    _ => None
                };
                samples[i].push(value);
          });
    });

    /*samples.iter().for_each(|v| {
        println!("{:#?}", v);
    });*/
    println!("finished reading");              

    let data = samples.into_par_iter()
                .map(|sample| {
                    let sample = sample.into_iter();
                    ChunkedArray::<Int8Type>::from_iter(sample)
                })
                .collect::<Vec<ChunkedArray<Int8Type>>>();
    let data = Arc::new(data);

    println!("{}", data[0].len());
    let mut matrix = (0_usize..num_samples).map(|_| repeat(0.0).take(num_samples).collect())
                                 .collect::<Vec<Vec<f64>>>();
    
    (0_usize..num_samples).combinations(2)
                        .collect::<Vec<Vec<usize>>>()
                        .into_par_iter()
                        .map(|indicies| {
                            let mask = &data[indicies[0]].is_not_null() & &data[indicies[1]].is_not_null();
                            let data1 = &data[indicies[0]].filter(&mask).unwrap();
                            let data2 = &data[indicies[1]].filter(&mask).unwrap();
                            let num_compares = (data1.len() * 2) as f64;
                            //println!("{}", data1.filter(&data1.equal(4).unwrap()).unwrap().len());
                            let diff: f64 = (data1 - data2).abs().cast(&DataType::Float64).unwrap().sum().unwrap();
                            diff / num_compares
                        }).collect::<Vec<f64>>()
                        .into_iter()
                        .zip((0_usize..num_samples).combinations(2).collect::<Vec<Vec<usize>>>())
                        .for_each(|(distance, indecies): (f64, Vec<usize>)| {
                            matrix[indecies[1]][indecies[0]] = distance;
                        });
    matrix.iter().for_each(|line| println!("{:?}", line));

    Ok(())
}
