use noodles::{bcf::{header::{string_maps::StringMap, StringMaps}, self},
              vcf::header::Header};
use std::io::{BufReader,
             };
use std::fs::File;
use std::error::Error;
use rayon::prelude::*;
use itertools::Itertools;
use peroxide::prelude::*;

fn main() -> Result<(), Box<dyn Error>> {
    let ir = File::open("/gpfs0/scratch/mvc002/info_challenge.bcf")?;
    let bf = BufReader::with_capacity(100 * 1024_usize, ir);
    let mut reader = bcf::Reader::new(bf);
    reader.read_file_format()?;

    let raw_header = reader.read_header()?;
    let header: Header = raw_header.parse()?;
    let string_maps: StringMaps = raw_header.parse()?;

    let combos = (0..29).combinations(2);
    println!("{:?}", combos);

    /*
    let mut record = bcf::Record::default();
    let mut i = 800;
        reader.read_record(&mut record)?;
        let vcf_record = record.genotypes()
                               .try_into_vcf_record_genotypes(&header, &string_maps.strings())?;
        //println!("{vcf_record}");
        //println!(""); */
    /*reader.records().par_bridge().map(|record| {
        let mut bin_geno = record.expect("error in data file")
                                 .genotypes()
                                 .to_owned();
        let rec_num = bin_geno.len() * 2 + 3;
        let buff = bin_geno.as_mut();
        buff.truncate(rec_num);
        let red: Vec<f64> = buff.split_off(3)
                .into_iter()
                .chunks(2)
                .into_iter()
                .map(|chunk| chunk.reduce(|acc, e| acc + e).unwrap() as f64)
                .collect();
        let matrix = 
        println!("{:?}", red);
    });*/
        /*
        i -= 1;
    } */
    //println!("{:?}", string_maps);
    Ok(())
}
