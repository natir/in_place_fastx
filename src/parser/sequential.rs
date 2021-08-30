#[macro_export(local_inner_macros)]
macro_rules! impl_sequential {
    ($name:ident, $producer:expr, $reader:expr,  $data_type:ty, $record:expr) => {
        pub struct $name {}

        impl $name {
            pub fn new() -> Self {
                Self {}
            }

            pub fn parse<P>(&mut self, path: P, data: &mut $data_type) -> $crate::error::Result<()>
            where
                P: AsRef<std::path::Path>,
            {
                self.with_blocksize($crate::DEFAULT_BLOCKSIZE, path, data)
            }

            pub fn with_blocksize<P>(
                &mut self,
                blocksize: u64,
                path: P,
                data: &mut $data_type,
            ) -> $crate::error::Result<()>
            where
                P: AsRef<std::path::Path>,
            {
                let mut producer = $producer(blocksize, path)?;

                while let Some(block) = producer.next_block()? {
                    self.block(block, data)?
                }

                Ok(())
            }

            fn block(
                &mut self,
                block: $crate::block::Block,
                data: &mut $data_type,
            ) -> $crate::error::Result<()> {
                let mut reader = $reader(block);

                while let Some(record) = reader.next_record()? {
                    self.record(record, data);
                }

                Ok(())
            }

            fn record(&self, record: $crate::block::Record, data: &mut $data_type) -> () {
                $record(record, data);
            }
        }
    };
}

#[macro_export(local_inner_macros)]
macro_rules! fasta_sequential {
    ($name:ident, $data_type:ty, $record:expr) => {
        impl_sequential!(
            $name,
            $crate::fasta::Producer::with_blocksize,
            $crate::fasta::Reader::new,
            $data_type,
            $record
        );
    }
}

#[macro_export(local_inner_macros)]
macro_rules! fastq_sequential {
    ($name:ident, $data_type:ty, $record:expr) => {
        impl_sequential!(
            $name,
            $crate::fastq::Producer::with_blocksize,
            $crate::fastq::Reader::new,
            $data_type,
            $record
        );
    }
}

#[cfg(test)]
mod tests {
    use crate::block;
    use crate::error;
    use crate::fasta;
    use crate::fastq;

    #[test]
    fn record_count_fasta() {
        impl_sequential!(
            FastaRecordCount,
            fasta::Producer::with_blocksize,
            fasta::Reader::new,
            u64,
            |_record: block::Record, counter: &mut u64| {
                *counter += 1;
            }
        );

        let mut counter = 0;

        let mut parser = FastaRecordCount::new();

        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150), &mut counter)
            .unwrap();

        assert_eq!(1_000, counter);
    }

    #[test]
    fn base_count_fasta() {
        impl_sequential!(
            FastaNucCount,
            fasta::Producer::with_blocksize,
            fasta::Reader::new,
            [u64; 4],
            |record: block::Record, bases: &mut [u64; 4]| {
                for nuc in record.sequence {
                    bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        );

        let mut bases = [0; 4];

        let mut parser = FastaNucCount::new();

        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150), &mut bases)
            .unwrap();

        assert_eq!([37378, 37548, 37548, 37526], bases);
    }

    #[test]
    fn record_count_fastq() {
        impl_sequential!(
            FastqRecordCount,
            fastq::Producer::with_blocksize,
            fastq::Reader::new,
            u64,
            |_record: block::Record, counter: &mut u64| {
                *counter += 1;
            }
        );

        let mut counter = 0;

        let mut parser = FastqRecordCount::new();

        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150), &mut counter)
            .unwrap();

        assert_eq!(1_000, counter);
    }

    #[test]
    fn base_count_fastq() {
        impl_sequential!(
            FastqNucCount,
            fastq::Producer::with_blocksize,
            fastq::Reader::new,
            [u64; 4],
            |record: block::Record, bases: &mut [u64; 4]| {
                for nuc in record.sequence {
                    bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        );

        let mut bases = [0; 4];

        let mut parser = FastqNucCount::new();

        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150), &mut bases)
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], bases);
    }
}
