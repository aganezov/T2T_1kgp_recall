use std::str;
use std::iter::FromIterator;
use std::fs::File;
use std::env::args;

use flate2::write::GzEncoder;
use flate2::Compression;

use ahash::AHashMap;

use fastq::{parse_path, Record};

fn split_fastq_into_lanes(path: &str, prefix: &str) {
    let mut file_map: AHashMap<String, GzEncoder<File>> = AHashMap::new();

    parse_path(Some(path), |parser| {
        parser.each(|record| {
            let header = str::from_utf8(record.head()).unwrap();
            let info = Vec::from_iter(header.split(":"));
            let file_name = String::from(prefix) + "." + info[2] +  ".00" + info[3] + ".fastq.gz";

            let e = file_map.entry(file_name.clone()).or_insert(
                GzEncoder::new(File::create(&file_name).unwrap(), Compression::default())
            );
            record.write(e).expect("Problem with writing fastq record");
            true
        }).expect("Invalid fastq file");
    }).expect("Invalid compression");

    for (_, fhand) in file_map.into_iter() {
        fhand.finish().unwrap();
    }
}

fn main() {
    let filename = args().nth(1);

    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => return,
        Some(name) => name
    };

    let prefix = args().nth(2);
    let prefix = match prefix {
        None => return,
        Some(prefix) => prefix
    };

    split_fastq_into_lanes(path, &prefix);
}
