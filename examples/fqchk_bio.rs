/* crate use */
use clap::Parser;

mod fqchk_utils;
use fqchk_utils::{Command, Data, PosInfo};

fn worker(record: bio::io::fastq::Record, data: &mut Data) {
    // Update min and max length
    if data.min_len > record.qual().len() as u64 {
        data.min_len = record.qual().len() as u64;
    }

    if data.max_len < record.qual().len() as u64 {
        data.max_len = record.qual().len() as u64;
    }

    data.seq_count += 1;
    data.average_len = (1.0 / data.seq_count as f64) * record.seq().len() as f64
        + (1.0 - (1.0 / data.seq_count as f64)) * data.average_len as f64;

    // Update count table length
    if record.qual().len() > data.pos_infos.len() {
        data.pos_infos.extend(vec![
            PosInfo::new(data.threshold);
            record.qual().len() - data.pos_infos.len()
        ]);
    }

    for (indice, (nuc, qual)) in record.seq().iter().zip(record.qual()).enumerate() {
        data.pos_infos[indice].quality_counts[(qual - data.phred_offset) as usize] += 1;
        data.pos_infos[indice].bases_count[(nuc >> 1 & 0b111) as usize] += 1;
    }
}

fn main() -> in_place_fastx::error::Result<()> {
    let params = Command::parse();

    let mut data = Data::new(params.qual_t);
    for input in params.inputs {
        let mut records = bio::io::fastq::Reader::from_bufread(std::io::BufReader::with_capacity(
            params.blocksize as usize,
            std::fs::File::open(input).unwrap(),
        ))
        .records();

        while let Some(Ok(record)) = records.next() {
            worker(record, &mut data);
        }
    }

    println!("{}", data);
    Ok(())
}
