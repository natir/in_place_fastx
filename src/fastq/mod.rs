use rayon::iter::ParallelIterator;

use rayon::iter::ParallelBridge;

pub mod block;

pub type Record<'a> = (&'a [u8], &'a [u8], &'a [u8], &'a [u8]);

use crate::error;

pub type RecordWorker<T> = fn(Record, &mut T) -> error::Result<()>;
pub type BlockWorker<T> = fn(block::Block, &mut T) -> error::Result<()>;

pub trait Parser {
    fn file<P>(&mut self, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        self.file_with_blocksize(8192, path)
    }

    fn file_with_blocksize<P>(&mut self, blocksize: u64, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        let mut producer = block::Producer::with_blocksize(blocksize, path)?;

        while let Some(block) = producer.next_block()? {
            self.block(block)?
        }

        Ok(())
    }

    fn block(&mut self, block: block::Block) -> error::Result<()> {
        let mut reader = block::Reader::new(block);

        while let Some(record) = reader.next_record()? {
            self.record(record)
        }

        Ok(())
    }

    fn record(&mut self, _record: Record) {}

    fn multithread_by_block<P, T>(
        &self,
        path: P,
        data: &T,
        worker: fn(Record, &T),
    ) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
        T: std::marker::Send + std::marker::Sync,
    {
        self.multithread_by_block_with_blocksize(8092, path, data, worker)
    }

    fn multithread_by_block_with_blocksize<P, T>(
        &self,
        blocksize: u64,
        path: P,
        data: &T,
        worker: fn(Record, &T),
    ) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
        T: std::marker::Send + std::marker::Sync,
    {
        let producer = block::Producer::with_blocksize(blocksize, path)?;

        match producer
            .par_bridge()
            .map(|block| match block {
                Ok(block) => {
                    let mut reader = block::Reader::new(block);
                    while let Some(record) = reader.next_record()? {
                        worker(record, data);
                    }
                    Ok(())
                }
                Err(e) => Err(e),
            })
            .find_any(|x| x.is_err())
        {
            Some(e) => e,
            None => Ok(()),
        }
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
        struct Counter {
            pub count: u64,
        }

        impl Parser for Counter {
            fn record(&mut self, _record: Record) {
                self.count += 1;
            }
        }

        let mut parser = Counter { count: 0 };
        parser.file(generate_fastq(42, 1_000, 150)).unwrap();

        assert_eq!(1_000, parser.count);
    }

    #[test]
    fn record_counter_parallel() {
        struct Counter {
            pub count: std::sync::atomic::AtomicU64,
        }

        fn worker(_record: Record, data: &std::sync::atomic::AtomicU64) {
            data.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        }

        impl Parser for Counter {}

        let parser = Counter {
            count: std::sync::atomic::AtomicU64::new(0),
        };

        parser
            .multithread_by_block(generate_fastq(42, 1_000, 150), &parser.count, worker)
            .unwrap();

        assert_eq!(1_000, parser.count.into_inner());
    }

    type BaseCount<T> = [T; 4];
    trait AbsBaseCount {
        fn new() -> Self;
    }

    impl AbsBaseCount for BaseCount<u64> {
        fn new() -> Self {
            [0, 0, 0, 0]
        }
    }

    impl AbsBaseCount for BaseCount<std::sync::atomic::AtomicU64> {
        fn new() -> Self {
            [
                std::sync::atomic::AtomicU64::new(0),
                std::sync::atomic::AtomicU64::new(0),
                std::sync::atomic::AtomicU64::new(0),
                std::sync::atomic::AtomicU64::new(0),
            ]
        }
    }

    #[test]
    fn base_count() {
        struct Counter {
            pub bases: BaseCount<u64>,
        }

        impl Parser for Counter {
            fn record(&mut self, record: Record) {
                for nuc in record.1 {
                    self.bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        }

        let mut parser = Counter {
            bases: BaseCount::new(),
        };
        parser.file(generate_fastq(42, 1_000, 150)).unwrap();

        assert_eq!([37301, 37496, 37624, 37579], parser.bases);
    }

    #[test]
    fn base_count_parallel() {
        struct Counter {
            pub bases: BaseCount<std::sync::atomic::AtomicU64>,
        }

        fn worker(record: Record, data: &BaseCount<std::sync::atomic::AtomicU64>) {
            for nuc in record.1 {
                data[(nuc >> 1 & 0b11) as usize].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        }

        impl Parser for Counter {}

        let parser = Counter {
            bases: BaseCount::new(),
        };

        parser
            .multithread_by_block(generate_fastq(42, 1_000, 150), &parser.bases, worker)
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], unsafe {
            std::mem::transmute::<BaseCount<std::sync::atomic::AtomicU64>, BaseCount<u64>>(
                parser.bases,
            )
        });
    }
}
