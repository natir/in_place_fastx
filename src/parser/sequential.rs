/* project use */
use crate::block;
use crate::error;

/// Trait allow sequential parsing of fastq
///
/// Reading is perform by block. Parser map a block of file in memory, this block is resize to remove incomplete record.
/// For each block record position is extract and `record` function is call on it.
pub trait Sequential<B, R>
where
    B: block::Producer,
    R: block::Reader,
{
    /// Parse file indicate by path with default blocksize [crate::DEFAULT_BLOCKSIZE]
    fn parse<P>(&mut self, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        self.with_blocksize(crate::DEFAULT_BLOCKSIZE, path)
    }

    /// Parse file indicate by path with selected blocksize
    fn with_blocksize<P>(&mut self, blocksize: u64, path: P) -> error::Result<()>
    where
        P: AsRef<std::path::Path>,
    {
        let mut producer = B::with_blocksize(blocksize, path)?;

        while let Some(block) = producer.next_block()? {
            self.block(block)?
        }

        Ok(())
    }

    /// Method call to parse a block
    fn block(&mut self, block: crate::block::Block) -> error::Result<()> {
        let mut reader = R::new(block);

        while let Some(record) = reader.next_record()? {
            self.record(record)
        }

        Ok(())
    }

    /// Method call to parse a record
    fn record(&mut self, _record: crate::block::Record);
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::fasta;
    use crate::fastq;
    use crate::parser::tests::AbsBaseCount;

    impl crate::parser::tests::AbsBaseCount for crate::parser::tests::BaseCount<u64> {
        fn new() -> Self {
            [0, 0, 0, 0]
        }
    }

    #[test]
    fn record_count_fasta() {
        struct Counter {
            pub count: u64,
        }

        impl<'a> Sequential<fasta::Producer, fasta::Reader> for Counter {
            fn record(&mut self, _record: block::Record) {
                self.count += 1;
            }
        }

        let mut parser = Counter { count: 0 };
        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150))
            .unwrap();

        assert_eq!(1_000, parser.count);
    }

    #[test]
    fn base_count_fasta() {
        struct Counter {
            pub bases: crate::parser::tests::BaseCount<u64>,
        }

        impl<'a> Sequential<fasta::Producer, fasta::Reader> for Counter {
            fn record(&mut self, record: block::Record) {
                for nuc in record.sequence {
                    self.bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        }

        let mut parser = Counter {
            bases: crate::parser::tests::BaseCount::new(),
        };
        parser
            .parse(crate::tests::generate_fasta(42, 1_000, 150))
            .unwrap();

        assert_eq!([37378, 37548, 37548, 37526], parser.bases);
    }

    #[test]
    fn record_count_fastq() {
        struct Counter {
            pub count: u64,
        }

        impl<'a> Sequential<fastq::Producer, fastq::Reader> for Counter {
            fn record(&mut self, _record: block::Record) {
                self.count += 1;
            }
        }

        let mut parser = Counter { count: 0 };
        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150))
            .unwrap();

        assert_eq!(1_000, parser.count);
    }

    #[test]
    fn base_count_fastq() {
        struct Counter {
            pub bases: crate::parser::tests::BaseCount<u64>,
        }

        impl<'a> Sequential<fastq::Producer, fastq::Reader> for Counter {
            fn record(&mut self, record: block::Record) {
                for nuc in record.sequence {
                    self.bases[(nuc >> 1 & 0b11) as usize] += 1;
                }
            }
        }

        let mut parser = Counter {
            bases: crate::parser::tests::BaseCount::new(),
        };
        parser
            .parse(crate::tests::generate_fastq(42, 1_000, 150))
            .unwrap();

        assert_eq!([37301, 37496, 37624, 37579], parser.bases);
    }
}
