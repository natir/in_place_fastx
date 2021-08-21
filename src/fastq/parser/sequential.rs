/* project use */
use crate::error;
use crate::fastq;

/// Trait allow sequential parsing of fastq
///
/// Blank implementation, do nothing you should reimplement record.
///
/// Reading is perform by block (by default of 8192 bytes). Parser get a block, this block is resize to remove incomplete record
pub trait Sequential {
    /// Open a file and run record function on each function blocksize is 8192
    fn parse<P>(&mut self, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        self.with_blocksize(8192, path)
    }

    fn with_blocksize<P>(&mut self, blocksize: u64, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        let mut producer = fastq::block::Producer::with_blocksize(blocksize, path)?;

        while let Some(block) = producer.next_block()? {
            self.block(block)?
        }

        Ok(())
    }

    fn block(&mut self, block: fastq::Block) -> error::Result<()> {
        let mut reader = fastq::block::Reader::new(block);

        while let Some(record) = reader.next_record()? {
            self.record(record)
        }

        Ok(())
    }

    fn record(&mut self, _record: fastq::Record);
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::fastq::parser::tests::AbsBaseCount;

    #[test]
    fn record_count() {
        struct Counter {
            pub count: u64,
        }

        impl Sequential for Counter {
            fn record(&mut self, _record: fastq::Record) {
                self.count += 1;
            }
        }

        let mut parser = Counter { count: 0 };
        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150))
            .unwrap();

        assert_eq!(1_000, parser.count);
    }

    impl crate::fastq::parser::tests::AbsBaseCount for crate::fastq::parser::tests::BaseCount<u64> {
        fn new() -> Self {
            [0, 0, 0, 0]
        }
    }

    #[test]
    fn base_count() {
        struct Counter {
            pub bases: crate::fastq::parser::tests::BaseCount<u64>,
        }

        impl Sequential for Counter {
            fn record(&mut self, record: fastq::Record) {
                for nuc in record.sequence {
                    self.bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        }

        let mut parser = Counter {
            bases: crate::fastq::parser::tests::BaseCount::new(),
        };
        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150))
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], parser.bases);
    }
}
