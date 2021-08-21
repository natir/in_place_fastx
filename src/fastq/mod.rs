/* crate use */
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

/* project use */
use crate::error;

/* mod declaration */
pub mod block;

/* type declaratino*/
/// Record store a fastq record in public field
pub struct Record<'a> {
    pub comment: &'a [u8],
    pub sequence: &'a [u8],
    pub plus: &'a [u8],
    pub quality: &'a [u8],
}

/// Block reperesent a section of file memory mapped in file it's almost a &[u8]
pub struct Block {
    mem: memmap::Mmap,
    end: usize,
}

impl Block {
    pub fn new(end: usize, mem: memmap::Mmap) -> Self {
        Self { mem, end }
    }

    pub fn data(&self) -> &[u8] {
        &self.mem[..self.end]
    }

    pub fn len(&self) -> usize {
        self.mem[..self.end].len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Trait allow parsing of fastq
///
/// Blank implementation, do nothing you should reimplement record.
///
/// Reading is perform by block (by default of 8192 bytes). Parser get a block, this block is resize to remove incomplete record
pub trait Parser {
    /// Open a file and run record function on each function blocksize is 8192
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

    fn block(&mut self, block: Block) -> error::Result<()> {
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
    fn block() {
        let file = generate_fastq(42, 1_000, 50);

        let data = unsafe {
            memmap::MmapOptions::new()
                .offset(0)
                .len(500)
                .map(file.as_file())
                .unwrap()
        };

        let block = Block::new(200, data);

        assert_eq!(
            block.data(),
            &[
                64, 48, 10, 84, 84, 65, 71, 65, 84, 84, 65, 84, 65, 71, 84, 65, 67, 71, 71, 84, 65,
                84, 65, 71, 84, 71, 71, 84, 84, 65, 67, 84, 65, 84, 71, 84, 65, 71, 67, 67, 84, 65,
                65, 71, 84, 71, 71, 67, 71, 67, 67, 67, 71, 10, 43, 48, 10, 47, 71, 70, 88, 46, 85,
                77, 44, 112, 49, 123, 45, 124, 116, 78, 95, 45, 70, 77, 111, 38, 40, 45, 62, 43,
                101, 76, 57, 87, 71, 115, 79, 34, 33, 95, 45, 81, 53, 97, 118, 75, 59, 39, 121, 63,
                51, 33, 55, 113, 74, 10, 64, 49, 10, 84, 65, 84, 65, 65, 84, 67, 67, 71, 71, 65,
                67, 71, 71, 67, 65, 84, 71, 67, 71, 67, 65, 71, 71, 67, 65, 84, 71, 67, 67, 84, 65,
                84, 65, 84, 84, 67, 84, 65, 84, 71, 65, 67, 65, 71, 67, 65, 71, 71, 65, 10, 43, 49,
                10, 76, 83, 104, 99, 98, 69, 37, 69, 42, 113, 116, 37, 85, 109, 74, 45, 86, 112,
                102, 93, 61, 76, 100, 102, 80, 95, 58, 119, 83, 73, 40, 94, 85, 51, 60
            ]
        );

        assert_eq!(block.len(), 200);
        assert!(!block.is_empty());
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
                for nuc in record.sequence {
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
            for nuc in record.sequence {
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
