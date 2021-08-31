//! Struct that extract part of file (called block), each block is read in parallel

#[macro_export(local_inner_macros)]
macro_rules! impl_sharedstate {
    ($name:ident, $producer:expr, $reader:expr, $data_type:ty, $record:expr) => {
        pub struct $name {}

        impl $name {
            pub fn new() -> Self {
                Self {}
            }

            pub fn parse<P>(&mut self, path: P, data: &$data_type) -> $crate::error::Result<()>
            where
                P: AsRef<std::path::Path>,
            {
                self.with_blocksize($crate::DEFAULT_BLOCKSIZE, path, data)
            }

            fn with_blocksize<P>(
                &self,
                blocksize: u64,
                path: P,
                data: &$data_type,
            ) -> $crate::error::Result<()>
            where
                P: AsRef<std::path::Path>,
            {
                let producer = $producer(blocksize, path)?;

                match producer
                    .par_bridge()
                    .map(|block| {
                        let mut reader = $reader(block?);
                        while let Some(record) = reader.next_record()? {
                            $record(record, data);
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
    };
}

#[macro_export(local_inner_macros)]
macro_rules! fasta_sharedstate {
    ($name:ident, $data_type:ty, $record:expr) => {
        impl_sharedstate!(
            $name,
            $crate::fasta::Producer::with_blocksize,
            $crate::fasta::Reader::new,
            $data_type,
            $record
        );
    };
}

#[macro_export(local_inner_macros)]
macro_rules! fastq_sharedstate {
    ($name:ident, $data_type:ty, $record:expr) => {
        impl_sharedstate!(
            $name,
            $crate::fastq::Producer::with_blocksize,
            $crate::fastq::Reader::new,
            $data_type,
            $record
        );
    };
}

#[cfg(test)]
mod tests {
    /* crate use */
    use rayon::iter::ParallelBridge;
    use rayon::iter::ParallelIterator;

    /* project use */
    use crate::block;

    #[test]
    fn record_count_fasta() {
        fasta_sharedstate!(
            FastaRecordCount,
            std::sync::atomic::AtomicU64,
            |_record: block::Record, counter: &std::sync::atomic::AtomicU64| {
                counter.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        );

        let counter = std::sync::atomic::AtomicU64::new(0);

        let mut parser = FastaRecordCount::new();

        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150), &counter)
            .unwrap();

        assert_eq!(1000, counter.into_inner());
    }

    #[test]
    fn base_count_fasta() {
        fasta_sharedstate!(
            FastaRecordCount,
            [std::sync::atomic::AtomicU64; 4],
            |record: block::Record, counter: &[std::sync::atomic::AtomicU64; 4]| {
                for nuc in record.sequence {
                    counter[(nuc >> 1 & 0b11) as usize]
                        .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                }
            }
        );

        let counter = [
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
        ];

        let mut parser = FastaRecordCount::new();

        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150), &counter)
            .unwrap();

        assert_eq!([37378, 37548, 37548, 37526], unsafe {
            std::mem::transmute::<[std::sync::atomic::AtomicU64; 4], [u64; 4]>(counter)
        });
    }

    #[test]
    fn record_count_fastq() {
        fastq_sharedstate!(
            FastqRecordCount,
            std::sync::atomic::AtomicU64,
            |_record: block::Record, counter: &std::sync::atomic::AtomicU64| {
                counter.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
        );

        let counter = std::sync::atomic::AtomicU64::new(0);

        let mut parser = FastqRecordCount::new();

        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150), &counter)
            .unwrap();

        assert_eq!(1000, counter.into_inner());
    }

    #[test]
    fn base_count_fastq() {
        fastq_sharedstate!(
            FastqRecordCount,
            [std::sync::atomic::AtomicU64; 4],
            |record: block::Record, counter: &[std::sync::atomic::AtomicU64; 4]| {
                for nuc in record.sequence {
                    counter[(nuc >> 1 & 0b11) as usize]
                        .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                }
            }
        );

        let counter = [
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
            std::sync::atomic::AtomicU64::new(0),
        ];

        let mut parser = FastqRecordCount::new();

        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150), &counter)
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], unsafe {
            std::mem::transmute::<[std::sync::atomic::AtomicU64; 4], [u64; 4]>(counter)
        });
    }
}
