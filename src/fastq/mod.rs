pub mod block;

pub type Record<'a> = (&'a [u8], &'a [u8], &'a [u8], &'a [u8]);

use crate::error;

pub struct Parser {
    block_producer: block::Producer,
}

pub type RecordWorker<T> = fn(Record, &mut T) -> error::Result<()>;
pub type BlockWorker<T> = fn(block::Block, &mut T) -> error::Result<()>;

impl Parser {
    pub fn new<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        Ok(Parser {
            block_producer: block::Producer::new(path)?,
        })
    }

    pub fn by_record<T>(&mut self, worker: RecordWorker<T>, data: &mut T) -> error::Result<()> {
        while let Some(current_block) = self.block_producer.next_block()? {
            let mut reader = block::Reader::new(current_block);

            while let Some(record) = reader.next()? {
                worker(record, data)?;
            }
        }

        Ok(())
    }

    pub fn by_block<T>(&mut self, worker: BlockWorker<T>, data: &mut T) -> error::Result<()> {
        while let Some(current_block) = self.block_producer.next_block()? {
            worker(current_block, data)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;

    use rand::Rng;
    use rand::SeedableRng;

    fn generate_fastq(seed: u64, nb_seq: usize, length: usize) -> tempfile::NamedTempFile {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

        let mut file = tempfile::NamedTempFile::new().unwrap();

        let dna = [b'A', b'C', b'T', b'G'];
        let qual = (0..94).collect::<Vec<u8>>();

        for i in 0..nb_seq {
            let dna_seq = (0..length)
                .map(|_| dna[rng.gen_range(0..4)] as char)
                .collect::<String>();
            let qual_seq = (0..length)
                .map(|_| (qual[rng.gen_range(0..94)] + 33) as char)
                .collect::<String>();

            writeln!(file, "@{}\n{}\n+{}\n{}", i, dna_seq, i, qual_seq).unwrap();
        }

        file
    }

    #[test]
    fn record_counter() {
        fn counter(_: Record, count: &mut u64) -> error::Result<()> {
            *count += 1;

            Ok(())
        }

        let mut count = 0;
        let mut parser = Parser::new(generate_fastq(42, 1_000, 150)).unwrap();

        parser.by_record(counter, &mut count).unwrap();

        assert_eq!(1_000, count);
    }

    #[test]
    fn record_counter_block() {
        fn counter(block: block::Block, count: &mut u64) -> error::Result<()> {
            let mut reader = block::Reader::new(block);

            while reader.next()?.is_some() {
                *count += 1;
            }

            Ok(())
        }

        let mut count = 0;
        let mut parser = Parser::new(generate_fastq(42, 1_000, 150)).unwrap();

        parser.by_block(counter, &mut count).unwrap();

        assert_eq!(1_000, count);
    }
}
