//! Struct that extract part of file (called block) and read it as fasta file.

/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::block::AbcProducer;
use crate::block::Block;
use crate::error;

/// Struct that produce a [Block](super::Block) of file, this block contains complete record.
pub struct Producer {
    offset: u64,
    blocksize: u64,
    file: std::fs::File,
    file_length: u64,
}

impl Producer {
    /// Build a [Block](super::Block) producer, default [Block](super::Block) size is 2^16 bytes.
    pub fn new<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        Producer::with_blocksize(crate::DEFAULT_BLOCKSIZE, path)
    }

    /// Build a [Block](super::Block) producer, with a specific [Block](super::Block) size warning this block size must be larger than two records.
    pub fn with_blocksize<P>(blocksize: u64, path: P) -> Result<Self, error::Error>
    where
        P: AsRef<std::path::Path>,
    {
        Ok(Producer {
            offset: 0,
            blocksize: Self::fix_blocksize::<P>(&path, blocksize)?,
            file_length: Self::filesize::<P>(&path)?,
            file: std::fs::File::open(path).map_err(|source| error::Error::OpenFile { source })?,
        })
    }
}

impl AbcProducer for Producer {
    /// Search the begin of the partial record at the end of [Block](super::Block)
    fn correct_block_size(block: &[u8]) -> error::Result<u64> {
        let mut end = block.len();

        for _ in 0..2 {
            end = block[..end]
                .rfind_byte(b'\n')
                .ok_or(error::Error::NoNewLineInBlock)?;

            if end + 1 < block.len() && block[end + 1] == b'>' {
                return Ok((end + 1) as u64);
            }
        }

        Err(error::Error::NotAFastaFile)
    }

    /// Get current value of offset
    fn offset(&self) -> u64 {
        self.offset
    }

    /// Get file length
    fn file_length(&self) -> u64 {
        self.file_length
    }

    /// Get file
    fn file(&self) -> &std::fs::File {
        &self.file
    }

    /// Get blocksize
    fn blocksize(&self) -> u64 {
        self.blocksize
    }

    /// Set value of offset
    fn set_offset(&mut self, value: u64) {
        self.offset = value;
    }
}

impl Iterator for Producer {
    type Item = error::Result<Block>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_block() {
            Ok(Some(block)) => Some(Ok(block)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Struct that read [Block](Block) and produce [Record](super::Record)
pub struct Reader {
    offset: usize,
    block: Block,
}

impl Reader {
    /// Create a new [Block](Block) reader from [Block](Block) get in parameter
    pub fn new(block: Block) -> Self {
        Reader { offset: 0, block }
    }

    /// Search next end of line
    fn get_line(&self) -> error::Result<std::ops::Range<usize>> {
        let next = self.block.data()[self.offset..]
            .find_byte(b'\n')
            .ok_or(error::Error::PartialRecord)?;
        let range = self.offset..self.offset + next;

        Ok(range)
    }

    /// Produce [Record](super::Record) until block is empty
    pub fn next_record(&mut self) -> error::Result<Option<super::Record<'_>>> {
        if self.offset == self.block.len() {
            Ok(None)
        } else {
            let comment = &self.block.data()[self.get_line()?];
            self.offset += comment.len() + 1;

            let sequence = &self.block.data()[self.get_line()?];
            self.offset += sequence.len() + 1;

            Ok(Some(super::Record { comment, sequence }))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod producer {
        use super::*;

        use std::io::Write;

        #[test]
        fn new() {
            let mut tmp = Producer::new(crate::tests::generate_fasta(42, 100, 150)).unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 15490);

            assert!(tmp.next_block().unwrap().is_none());
        }

        #[test]
        fn with_blocksize() {
            let mut tmp =
                Producer::with_blocksize(463, crate::tests::generate_fasta(42, 1_000, 150))
                    .unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 462);

            assert_eq!(
                String::from_utf8(block.data().to_vec()).unwrap(),
                ">0
TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA
>1
GATACGTTTGGGGCAACCCGTAGCACGACCGGCTATGTGTTTTCTTGGACATAGTTTCGTCCACGATATATACAAGGACGCTTGGGAATAGGGCAGCGGAGTTATCGTGTACCTCCTAGCTTTTAGTTGTGCTTTAACAGTGTAACATTG
>2
GGACGCTATTACTCGCCGGTGAGGCGGTCTTCCTTGACTATACCGATCGTGGAGTTCATGCGCGCGGATCCCTCAGCGTTCTCGGGAAGCGCGAACAGAGCGTCCCCTTATACTAATTCCACGCAATGTACTCGCTTACGATTGCAATTT
".to_string()
            );
        }

        #[test]
        fn with_blocksize_buffer_larger_file() {
            let mut tmp =
                Producer::with_blocksize(8092, crate::tests::generate_fasta(44, 2, 150)).unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 308);
        }

        #[test]
        fn get_all_block() {
            let mut tmp = Producer::new(crate::tests::generate_fasta(42, 1_000, 150)).unwrap();

            let mut block_length = Vec::new();
            while let Ok(Some(block)) = tmp.next_block() {
                block_length.push(block.len());
            }

            assert_eq!(block_length, vec![65410, 65520, 24960]);
        }

        #[test]
        fn check_block() {
            let mut tmp =
                Producer::with_blocksize(400, crate::tests::generate_fasta(42, 5, 150)).unwrap();

            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok(">0
TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA
>1
GATACGTTTGGGGCAACCCGTAGCACGACCGGCTATGTGTTTTCTTGGACATAGTTTCGTCCACGATATATACAAGGACGCTTGGGAATAGGGCAGCGGAGTTATCGTGTACCTCCTAGCTTTTAGTTGTGCTTTAACAGTGTAACATTG
".to_string())
            );
            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok(">2
GGACGCTATTACTCGCCGGTGAGGCGGTCTTCCTTGACTATACCGATCGTGGAGTTCATGCGCGCGGATCCCTCAGCGTTCTCGGGAAGCGCGAACAGAGCGTCCCCTTATACTAATTCCACGCAATGTACTCGCTTACGATTGCAATTT
>3
GCAAATGAGGACCATCGTCCCTTCATATCGTCGATAAGGAGCTTGATCCTGAATGTCCCTCAATCCGCGGCATGGCTAAGTACCACCGTGGATGTAAATTTTTCAGTCGTCTCTTCATACTGTTCCTGTACTGTCAGGGATGCTCCCTTT
".to_string())
            );
            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok(">4
CACAGAGCTCGTATAATCAGTAAACGCCACGGTCCTTTCTCTGTTAACCGCTATGCTAGAGTTCGACGGATTGCGAACTGTTTATAAAGGTATTATTGGTGGAAGATCGACGCAGTTGGTGCCGCAGGAACCGGTCAACTTAATGCTGAG
".to_string())
            );
            assert!(tmp.next_block().is_ok());
            assert!(tmp.next_block().unwrap().is_none());
        }

        #[test]
        fn not_a_fasta() {
            let mut file = tempfile::NamedTempFile::new().unwrap();

            file.write(
                b"Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Vivamus ut nulla eget diam eleifend bibendum.
Praesent porta sapien id tortor hendrerit, a hendrerit dolor commodo. Donec sed elit enim.",
            )
            .unwrap();

            let mut producer = Producer::with_blocksize(150, file.path()).unwrap();
            assert!(producer.next_block().is_err());

            {
                let mut rewrite = file.reopen().unwrap();
                rewrite
                    .write(
                        b"+FAILLED FILE
+3
+TTGGGCATGAGGTTCA
@3ueauie
+~vGLKg+n!*iJ\\K
@iuiea
",
                    )
                    .unwrap();
            }

            let mut producer = Producer::with_blocksize(82, file.path()).unwrap();

            assert!(producer.next_block().is_err());

            let mut producer = Producer::with_blocksize(82, file).unwrap();
            assert!(producer.next().is_some());
            assert!(producer.next().unwrap().is_err());
        }
    }

    mod reader {
        use super::*;

        #[test]
        fn iterate_over_seq() {
            let mut producer =
                Producer::with_blocksize(500, crate::tests::generate_fasta(42, 5, 150)).unwrap();

            let mut comments = Vec::new();
            let mut seqs = Vec::new();

            while let Ok(Some(block)) = producer.next_block() {
                let mut reader = Reader::new(block);

                while let Ok(Some(record)) = reader.next_record() {
                    comments.push(String::from_utf8(record.comment.to_vec()).unwrap());
                    seqs.push(String::from_utf8(record.sequence.to_vec()).unwrap());
                }
            }

            assert_eq!(
                comments,
                vec![
                    ">0".to_string(),
                    ">1".to_string(),
                    ">2".to_string(),
                    ">3".to_string(),
                    ">4".to_string()
                ]
            );
            assert_eq!(
                seqs,
                vec![
"TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA".to_string(), "GATACGTTTGGGGCAACCCGTAGCACGACCGGCTATGTGTTTTCTTGGACATAGTTTCGTCCACGATATATACAAGGACGCTTGGGAATAGGGCAGCGGAGTTATCGTGTACCTCCTAGCTTTTAGTTGTGCTTTAACAGTGTAACATTG".to_string(), "GGACGCTATTACTCGCCGGTGAGGCGGTCTTCCTTGACTATACCGATCGTGGAGTTCATGCGCGCGGATCCCTCAGCGTTCTCGGGAAGCGCGAACAGAGCGTCCCCTTATACTAATTCCACGCAATGTACTCGCTTACGATTGCAATTT".to_string(), "GCAAATGAGGACCATCGTCCCTTCATATCGTCGATAAGGAGCTTGATCCTGAATGTCCCTCAATCCGCGGCATGGCTAAGTACCACCGTGGATGTAAATTTTTCAGTCGTCTCTTCATACTGTTCCTGTACTGTCAGGGATGCTCCCTTT".to_string(), "CACAGAGCTCGTATAATCAGTAAACGCCACGGTCCTTTCTCTGTTAACCGCTATGCTAGAGTTCGACGGATTGCGAACTGTTTATAAAGGTATTATTGGTGGAAGATCGACGCAGTTGGTGCCGCAGGAACCGGTCAACTTAATGCTGAG".to_string()]
            );
        }
    }
}
