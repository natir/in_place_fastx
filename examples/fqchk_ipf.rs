/* crate use */
use clap::Parser;

mod fqchk_utils;
use fqchk_utils::{Command, Data, PosInfo};

in_place_fastx::fastq_sequential!(
    Sequential,
    Data,
    |record: in_place_fastx::block::Record, data: &mut Data| {
        // Update min and max length
        if data.min_len > record.quality.len() as u64 {
            data.min_len = record.quality.len() as u64;
        }

        if data.max_len < record.quality.len() as u64 {
            data.max_len = record.quality.len() as u64;
        }

        data.seq_count += 1;
        data.average_len = (1.0 / data.seq_count as f64) * record.sequence.len() as f64
            + (1.0 - (1.0 / data.seq_count as f64)) * data.average_len as f64;

        // Update count table length
        if record.quality.len() > data.pos_infos.len() {
            data.pos_infos.extend(vec![
                PosInfo::new(data.threshold);
                record.quality.len() - data.pos_infos.len()
            ]);
        }

        for (indice, (nuc, qual)) in record.sequence.iter().zip(record.quality).enumerate() {
            data.pos_infos[indice].quality_counts[(qual - data.phred_offset) as usize] += 1;
            data.pos_infos[indice].bases_count[(nuc >> 1 & 0b111) as usize] += 1;
        }
    }
);

fn main() -> in_place_fastx::error::Result<()> {
    let params = Command::parse();

    let mut data = Data::new(params.qual_t);
    let mut parser = Sequential::new();
    for input in params.inputs {
        parser.with_blocksize(params.blocksize, input, &mut data)?;
    }

    println!("{}", data);
    Ok(())
}
