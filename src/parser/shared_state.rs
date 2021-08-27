//! Struct that extract part of file (called block), each block is read in parallel

/* crate use */
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;

/* project use */
use crate::block;
use crate::error;

/// Trait allow parallel parsing of block.
///
/// Reading is perform by block. Parser map a block of file in memory, this block is resize to remove incomplete record.
/// Block are read with a parallel iterator, each block is process in paralelle with function worker.
/// data should be able to be share between process.
pub trait SharedState<B, R>
where
    B: block::Producer + Iterator<Item = error::Result<block::Block>> + Send,
    R: block::Reader,
{
    /// Parse file indicate by path with default blocksize [crate::DEFAULT_BLOCKSIZE]
    fn parse<P, T>(&self, path: P, data: &T, worker: fn(block::Record, &T)) -> error::Result<()>
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
        worker: fn(block::Record, &T),
    ) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
        T: std::marker::Send + std::marker::Sync,
    {
        let producer = B::with_blocksize(blocksize, path)?;

        match producer
            .par_bridge()
            .map(|block| {
                let mut reader = R::new(block?);
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

    use crate::fasta;
    use crate::fastq;
    use crate::parser::tests::AbsBaseCount;

    impl crate::parser::tests::AbsBaseCount
        for crate::parser::tests::BaseCount<std::sync::atomic::AtomicU64>
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
    fn record_count_fasta() {
        struct Counter {}

        impl<'a> SharedState<fasta::Producer, fasta::Reader> for Counter {}

        fn worker(_record: block::Record, data: &std::sync::atomic::AtomicU64) {
            data.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        }

        let parser = Counter {};
        let data = std::sync::atomic::AtomicU64::new(0);

        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150), &data, worker)
            .unwrap();

        assert_eq!(1000, data.into_inner());
    }

    #[test]
    fn base_count_fasta() {
        struct Counter {
            pub bases: crate::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        }

        impl<'a> SharedState<fasta::Producer, fasta::Reader> for Counter {}

        fn worker(
            record: crate::block::Record,
            data: &crate::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        ) {
            for nuc in record.sequence {
                data[(nuc >> 1 & 0b11) as usize].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        }

        let parser = Counter {
            bases: crate::parser::tests::BaseCount::new(),
        };

        parser
            .with_blocksize(
                8192,
                crate::tests::generate_fasta(42, 1_000, 150),
                &parser.bases,
                worker,
            )
            .unwrap();

        assert_eq!([37378, 37548, 37548, 37526], unsafe {
            std::mem::transmute::<[std::sync::atomic::AtomicU64; 4], [u64; 4]>(parser.bases)
        });
    }

    #[test]
    fn record_count_fastq() {
        struct Counter {}

        impl<'a> SharedState<fastq::Producer, fastq::Reader> for Counter {}

        fn worker(_record: block::Record, data: &std::sync::atomic::AtomicU64) {
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
    fn base_count_fastq() {
        struct Counter {
            pub bases: crate::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        }

        impl<'a> SharedState<fastq::Producer, fastq::Reader> for Counter {}

        fn worker(
            record: block::Record,
            data: &crate::parser::tests::BaseCount<std::sync::atomic::AtomicU64>,
        ) {
            for nuc in record.sequence {
                data[(nuc >> 1 & 0b11) as usize].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        }

        let parser = Counter {
            bases: crate::parser::tests::BaseCount::new(),
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
