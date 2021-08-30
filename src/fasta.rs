//! Struct that extract part of file (called block) and read it as fasta file.

// #![feature(trace_macros)]

// trace_macros!(true);

/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::block;
use crate::error;

impl_producer!(Producer, |block: &[u8]| {
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
});

impl_reader!(Reader, |block: &'a block::Block, offset: &mut usize| {
    if *offset == block.len() {
        Ok(None)
    } else {
        let comment = &block.data()[Self::get_line(block, offset)?];
        *offset += comment.len() + 1;

        let sequence = &block.data()[Self::get_line(block, offset)?];
        *offset += sequence.len() + 1;

        let plus = &block.data()[*offset..*offset];
        let quality = &block.data()[*offset..*offset];

        Ok(Some(crate::block::Record {
            comment,
            sequence,
            plus,
            quality,
        }))
    }
});

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
