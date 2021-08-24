//! Struct that extract part of file (called block), each block is read in parallel

/* crate use */
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

/* project use */
use crate::error;
use crate::fastq;

/// Trait allow parallel parsing of block.
///
/// Reading is perform by block. Parser map a block of file in memory, this block is resize to remove incomplete record.
/// Block are read with a parallel iterator, each block is process in paralelle with function worker.
/// data should be able to be share between process.
pub trait SharedState {
    /// Parse file indicate by path with default blocksize [crate::DEFAULT_BLOCKSIZE]
    fn parse<P, T>(&self, path: P, data: &T, worker: fn(fastq::Record, &T)) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
        T: std::marker::Send + std::marker::Sync,
    {
        self.with_blocksize(8192, path, data, worker)
    }

    /// Parse file indicate by path with selected blocksize
    fn with_blocksize<P, T>(
        &self,
        blocksize: u64,
        path: P,
        data: &T,
        worker: fn(fastq::Record, &T),
    ) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
        T: std::marker::Send + std::marker::Sync,
    {
        let producer = fastq::block::Producer::with_blocksize(blocksize, path)?;

        match producer
            .par_bridge()
            .map(|block| {
                let mut reader = fastq::block::Reader::new(block?);
                while let Some(record) = reader.next_record()? {
                    worker(record, data);
                }
                Ok(())
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

    use crate::fastq::parser::tests::AbsBaseCount;

    impl crate::fastq::parser::tests::AbsBaseCount
        for crate::fastq::parser::tests::BaseCount<std::sync::atomic::AtomicU64>
    {
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
    fn record_count() {
        struct Counter {}

        impl SharedState for Counter {}

        fn worker(_record: crate::fastq::Record, data: &std::sync::atomic::AtomicU64) {
            data.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        }

        let parser = Counter {};
        let data = std::sync::atomic::AtomicU64::new(0);

        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150), &data, worker)
            .unwrap();

        assert_eq!(1000, data.into_inner());
    }

    #[test]
    fn base_count() {
        struct Counter {
            pub bases: crate::fastq::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        }

        impl SharedState for Counter {}

        fn worker(
            record: crate::fastq::Record,
            data: &crate::fastq::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        ) {
            for nuc in record.sequence {
                data[(nuc >> 1 & 0b11) as usize].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        }

        let parser = Counter {
            bases: crate::fastq::parser::tests::BaseCount::new(),
        };

        parser
            .with_blocksize(
                8192,
                crate::tests::generate_fastq(42, 1_000, 150),
                &parser.bases,
                worker,
            )
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], unsafe {
            std::mem::transmute::<[std::sync::atomic::AtomicU64; 4], [u64; 4]>(parser.bases)
        });
    }
}
