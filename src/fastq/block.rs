/* crate use */
use bstr::ByteSlice;

/* project use */
use crate::error;

pub struct Producer {
    offset: u64,
    blocksize: u64,
    file: std::fs::File,
    file_length: u64,
}

impl Producer {
    pub fn new<P>(path: P) -> error::Result<Self>
    where
        P: AsRef<std::path::Path>,
    {
        Producer::with_blocksize(8192, path)
    }

    pub fn with_blocksize<P>(mut blocksize: u64, path: P) -> Result<Self, error::Error>
    where
        P: AsRef<std::path::Path>,
    {
        let file_length = path
            .as_ref()
            .metadata()
            .map_err(|source| error::Error::MetaDataFile { source })?
            .len();

        blocksize = file_length.min(blocksize);

        Ok(Producer {
            offset: 0,
            blocksize,
            file_length,
            file: std::fs::File::open(path).map_err(|source| error::Error::OpenFile { source })?,
        })
    }

    pub fn next_block(&mut self) -> error::Result<Option<super::Block>> {
        if self.offset == self.file_length {
            Ok(None)
        } else if self.offset + self.blocksize >= self.file_length {
            let block = unsafe {
                memmap::MmapOptions::new()
                    .offset(self.offset)
                    .len((self.file_length - self.offset) as usize)
                    .map(&self.file)
                    .map_err(|source| error::Error::MapFile { source })?
            };

            self.offset = self.file_length;

            Ok(Some(super::Block::new(block.len(), block)))
        } else {
            let tmp = unsafe {
                memmap::MmapOptions::new()
                    .offset(self.offset)
                    .len(self.blocksize as usize)
                    .map(&self.file)
                    .map_err(|source| error::Error::MapFile { source })?
            };

            let blocksize = Producer::correct_block_size(&tmp)?;
            self.offset += blocksize;

            Ok(Some(super::Block::new(blocksize as usize, tmp)))
        }
    }

    fn correct_block_size(block: &[u8]) -> error::Result<u64> {
        let mut end = block.len();
        let mut seen_plus = false;

        for _ in 0..8 {
            end = block[..end]
                .rfind_byte(b'\n')
                .ok_or(error::Error::NoNewLineInBlock)?;

            if end + 1 < block.len() {
                seen_plus = if !seen_plus {
                    block[end + 1] == b'+'
                } else {
                    seen_plus
                };

                if seen_plus && block[end + 1] == b'@' {
                    return Ok((end + 1) as u64);
                }
            }
        }

        Err(error::Error::NotAFastqFile)
    }
}

impl Iterator for Producer {
    type Item = error::Result<super::Block>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_block() {
            Ok(Some(block)) => Some(Ok(block)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub struct Reader {
    offset: usize,
    block: super::Block,
}

impl Reader {
    pub fn new(block: super::Block) -> Self {
        Reader { offset: 0, block }
    }

    fn get_line(&self) -> error::Result<std::ops::Range<usize>> {
        let next = self.block.data()[self.offset..]
            .find_byte(b'\n')
            .ok_or(error::Error::PartialRecord)?;
        let range = self.offset..self.offset + next;

        Ok(range)
    }

    pub fn next_record(&mut self) -> error::Result<Option<super::Record<'_>>> {
        if self.offset == self.block.len() {
            Ok(None)
        } else {
            let comment = &self.block.data()[self.get_line()?];
            self.offset += comment.len() + 1;

            let sequence = &self.block.data()[self.get_line()?];
            self.offset += sequence.len() + 1;

            let plus = &self.block.data()[self.get_line()?];
            self.offset += plus.len() + 1;

            let quality = &self.block.data()[self.get_line()?];
            self.offset += quality.len() + 1;

            Ok(Some(super::Record {
                comment,
                sequence,
                plus,
                quality,
            }))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod producer {
        use super::*;

        #[test]
        fn new() {
            let mut tmp = Producer::new(crate::tests::generate_fastq(42, 1_000, 150)).unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 7730);
        }

        #[test]
        fn with_blocksize() {
            let mut tmp = Producer::with_blocksize(463, crate::tests::generate_fastq(42, 1_000, 150)).unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 308);

            assert_eq!(
                String::from_utf8(block.data().to_vec()).unwrap(),
                "@0
TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA
+0
^U3<L0PV{cnrl:8N`!:=mF8M0!0Ez/d{4j$=9f5rLeAQ-H.ptT3w6aMy8Z6O-dZ}2`UX=YJ-Etg`s&B%~F!kR7S8]@lTI<2-\\';v0}hU.(T*0VHGx,>Gze)*5rFv}k@RllOE2K)\"DQJvO)bl?(dDhh
".to_string()
            );
        }

        #[test]
        fn with_blocksize_buffer_larger_file() {
            let mut tmp = Producer::with_blocksize(8092, crate::tests::generate_fastq(44, 2, 150)).unwrap();

            let block = tmp.next_block().unwrap().unwrap();

            assert_eq!(block.len(), 616);
        }

        #[test]
        fn get_all_block() {
            let mut tmp = Producer::new(crate::tests::generate_fastq(42, 1_000, 150)).unwrap();

            let mut block_length = Vec::new();
            while let Ok(Some(block)) = tmp.next_block() {
                block_length.push(block.len());
            }

            assert_eq!(
                block_length,
                vec![
                    7730, 7750, 7750, 7750, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800,
                    7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800,
                    7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800, 7800,
                    7800
                ]
            );
        }

        #[test]
        fn check_block() {
            let mut tmp = Producer::with_blocksize(800, crate::tests::generate_fastq(42, 5, 150)).unwrap();

            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok("@0
TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA
+0
^U3<L0PV{cnrl:8N`!:=mF8M0!0Ez/d{4j$=9f5rLeAQ-H.ptT3w6aMy8Z6O-dZ}2`UX=YJ-Etg`s&B%~F!kR7S8]@lTI<2-\\';v0}hU.(T*0VHGx,>Gze)*5rFv}k@RllOE2K)\"DQJvO)bl?(dDhh
@1
AGTTATCGTGTACCTCCTAGCTTTTAGTTGTGCTTTAACAGTGTAACATTGGGACGCTATTACTCGCCGGTGAGGCGGTCTTCCTTGACTATACCGATCGTGGAGTTCATGCGCGCGGATCCCTCAGCGTTCTCGGGAAGCGCGAACAGA
+1
iCW?:KL~15\\E|MNRKY)S$?~~Ub}d)dY2LX:e@b^'<<$$e56W0fdV,<Y>Yd(J<5p6xt)z+OxuPXv?/_yH8z^%Sks1*nxm$<7*YdkvNPf:>YW=$uxZ)}[v/DlZm&EW(s(cMelx\"iEV3Hp]cz3%_T@\\Ms
".to_string())
            );
            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok("@2
AATGTCCCTCAATCCGCGGCATGGCTAAGTACCACCGTGGATGTAAATTTTTCAGTCGTCTCTTCATACTGTTCCTGTACTGTCAGGGATGCTCCCTTTCACAGAGCTCGTATAATCAGTAAACGCCACGGTCCTTTCTCTGTTAACCGC
+2
Ouf)Y|l;S1!tk[U9n2(NK=#Tmg,t+CSsUMaPs7{+V'~On{hc1NR}aY^YbYlg[}Fcq1K_$v1HG\"tRBj`||g>\\)2pU_QrnWO{c@;lw8B0+urH~$#K>:iSa-I-C#gDJ(9UUFubOeRHsDX3Ko`?T--iL+j
@3
TTGGGCATGAGGTTCACCGAAGGTGGCAGATATGCGCCATAAATTGACCAGGTTGTATCCAGCATTGGAAGAACGCACCCGGGGGGAGCACAGATCCTAGCAGTACACGCTCTGGGTCCTCTACGTCTTCGGAGTCTCTAGCTTGCCTTA
+3
:~vGLKg+n!*iJ\\.*wfxK)5Qmh%<:f^$nql7OB$}M/d.F,5[=>ZW*#f=0>Ao(@~lEHbSG1%,b_Uy2!zL%2GMB0O.t[#UcQ[]ufFZJ!K<kLgDNQlx)s8+75E^[-\"!1l[i<S#G\"B]xZ5?as*@8Laq`{@r
".to_string())
            );
            assert_eq!(
                String::from_utf8(tmp.next_block().unwrap().unwrap().data().to_vec()),
                Ok("@4
TCTATAGCTTGTCATGCCTTTCGATTGAGGGCGTCACCAAGCGAATTACTCGCTGATCCGTTCCCCGCCAATTCTGAGACTCCATAATCCTATCTGTGTCCCTAGGTGCCGTGTTCCGGTCGTGAGTTCGGCCCTTGCCTAAAGTTAATG
+4
myS=C|jEWnl,aC\\7!jv9[!vh/PAK}_H&<.o]qf|y@4L:?ssLg3N!v7/N5RyPHn=5%Fyh(4-Z:<6wf]^#t~0:i(X\\l-7]9olH9WLV~`L~JQ7ye7B1RSi2N$PuHwjj\\pb}J\\R~pe?j+X>R#p@MyqBBe*
".to_string())
            );
            assert!(tmp.next_block().is_ok());
            assert!(tmp.next_block().unwrap().is_none());
        }
    }

    mod reader {
        use super::*;

        #[test]
        fn iterate_over_seq() {
            let mut producer = Producer::with_blocksize(500, crate::tests::generate_fastq(42, 5, 150)).unwrap();

            let mut comments = Vec::new();
            let mut seqs = Vec::new();
            let mut pluss = Vec::new();
            let mut quals = Vec::new();

            while let Ok(Some(block)) = producer.next_block() {
                let mut reader = Reader::new(block);

                while let Ok(Some(record)) = reader.next_record() {
                    comments.push(String::from_utf8(record.comment.to_vec()).unwrap());
                    seqs.push(String::from_utf8(record.sequence.to_vec()).unwrap());
                    pluss.push(String::from_utf8(record.plus.to_vec()).unwrap());
                    quals.push(String::from_utf8(record.quality.to_vec()).unwrap());
                }
            }

            assert_eq!(
                comments,
                vec![
                    "@0".to_string(),
                    "@1".to_string(),
                    "@2".to_string(),
                    "@3".to_string(),
                    "@4".to_string()
                ]
            );
            assert_eq!(
                seqs,
                vec![
		    "TTAGATTATAGTACGGTATAGTGGTTACTATGTAGCCTAAGTGGCGCCCGTTGTAGAGGAATCCACTTATATAACACAGGTATAATCCGGACGGCATGCGCAGGCATGCCTATATTCTATGACAGCAGGATTATGGAAGATGGTGCTCTA".to_string(),
		    "AGTTATCGTGTACCTCCTAGCTTTTAGTTGTGCTTTAACAGTGTAACATTGGGACGCTATTACTCGCCGGTGAGGCGGTCTTCCTTGACTATACCGATCGTGGAGTTCATGCGCGCGGATCCCTCAGCGTTCTCGGGAAGCGCGAACAGA".to_string(),
		    "AATGTCCCTCAATCCGCGGCATGGCTAAGTACCACCGTGGATGTAAATTTTTCAGTCGTCTCTTCATACTGTTCCTGTACTGTCAGGGATGCTCCCTTTCACAGAGCTCGTATAATCAGTAAACGCCACGGTCCTTTCTCTGTTAACCGC".to_string(), "TTGGGCATGAGGTTCACCGAAGGTGGCAGATATGCGCCATAAATTGACCAGGTTGTATCCAGCATTGGAAGAACGCACCCGGGGGGAGCACAGATCCTAGCAGTACACGCTCTGGGTCCTCTACGTCTTCGGAGTCTCTAGCTTGCCTTA".to_string(), "TCTATAGCTTGTCATGCCTTTCGATTGAGGGCGTCACCAAGCGAATTACTCGCTGATCCGTTCCCCGCCAATTCTGAGACTCCATAATCCTATCTGTGTCCCTAGGTGCCGTGTTCCGGTCGTGAGTTCGGCCCTTGCCTAAAGTTAATG".to_string()     ]
            );
            assert_eq!(
                pluss,
                vec![
                    "+0".to_string(),
                    "+1".to_string(),
                    "+2".to_string(),
                    "+3".to_string(),
                    "+4".to_string()
                ]
            );
            assert_eq!(
                quals,
                vec![
		    "^U3<L0PV{cnrl:8N`!:=mF8M0!0Ez/d{4j$=9f5rLeAQ-H.ptT3w6aMy8Z6O-dZ}2`UX=YJ-Etg`s&B%~F!kR7S8]@lTI<2-\\';v0}hU.(T*0VHGx,>Gze)*5rFv}k@RllOE2K)\"DQJvO)bl?(dDhh".to_string(),
		    "iCW?:KL~15\\E|MNRKY)S$?~~Ub}d)dY2LX:e@b^'<<$$e56W0fdV,<Y>Yd(J<5p6xt)z+OxuPXv?/_yH8z^%Sks1*nxm$<7*YdkvNPf:>YW=$uxZ)}[v/DlZm&EW(s(cMelx\"iEV3Hp]cz3%_T@\\Ms".to_string(),
		    "Ouf)Y|l;S1!tk[U9n2(NK=#Tmg,t+CSsUMaPs7{+V'~On{hc1NR}aY^YbYlg[}Fcq1K_$v1HG\"tRBj`||g>\\)2pU_QrnWO{c@;lw8B0+urH~$#K>:iSa-I-C#gDJ(9UUFubOeRHsDX3Ko`?T--iL+j".to_string(),
		    ":~vGLKg+n!*iJ\\.*wfxK)5Qmh%<:f^$nql7OB$}M/d.F,5[=>ZW*#f=0>Ao(@~lEHbSG1%,b_Uy2!zL%2GMB0O.t[#UcQ[]ufFZJ!K<kLgDNQlx)s8+75E^[-\"!1l[i<S#G\"B]xZ5?as*@8Laq`{@r".to_string(),
		    "myS=C|jEWnl,aC\\7!jv9[!vh/PAK}_H&<.o]qf|y@4L:?ssLg3N!v7/N5RyPHn=5%Fyh(4-Z:<6wf]^#t~0:i(X\\l-7]9olH9WLV~`L~JQ7ye7B1RSi2N$PuHwjj\\pb}J\\R~pe?j+X>R#p@MyqBBe*".to_string(),
		]
            );
        }
    }
}
