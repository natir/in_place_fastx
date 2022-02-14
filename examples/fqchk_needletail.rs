/* crate use */
use clap::Parser;

use needletail::FastxReader;

mod fqchk_utils;
use fqchk_utils::{Command, Data, PosInfo};

fn worker(record: needletail::parser::SequenceRecord, data: &mut Data) {
    let quality = record.qual().unwrap();

    // Update min and max length
    if data.min_len > quality.len() as u64 {
        data.min_len = quality.len() as u64;
    }

    if data.max_len < quality.len() as u64 {
        data.max_len = quality.len() as u64;
    }

    data.seq_count += 1;
    data.average_len = (1.0 / data.seq_count as f64) * record.seq().len() as f64
        + (1.0 - (1.0 / data.seq_count as f64)) * data.average_len as f64;

    // Update count table length
    if quality.len() > data.pos_infos.len() {
        data.pos_infos.extend(vec![
            PosInfo::new(data.threshold);
            quality.len() - data.pos_infos.len()
        ]);
    }

    for (indice, (nuc, qual)) in record.seq().iter().zip(quality).enumerate() {
        data.pos_infos[indice].quality_counts[(qual - data.phred_offset) as usize] += 1;
        data.pos_infos[indice].bases_count[(nuc >> 1 & 0b111) as usize] += 1;
    }
}

fn main() -> in_place_fastx::error::Result<()> {
    let params = Command::parse();

    let mut data = Data::new(params.qual_t);
    for input in params.inputs {
        let mut reader = needletail::parser::FastqReader::with_capacity(
            std::fs::File::open(input).unwrap(),
            params.blocksize as usize,
        );

        while let Some(Ok(record)) = reader.next() {
            worker(record, &mut data);
        }
    }

    println!("{}", data);
    Ok(())
}
