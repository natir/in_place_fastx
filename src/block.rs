//! Struct that extract part of file (called block) and read it as fastx file.
/* crate use */

/* project use */

/// Block reperesent a section of file memory mapped in file
#[derive(Debug)]
pub struct Block {
    mem: memmap::Mmap,
    end: usize,
}

impl Block {
    /// Create a new Block
    pub fn new(end: usize, mem: memmap::Mmap) -> Self {
        Self { mem, end }
    }

    /// Acces to data owned by block
    pub fn data(&self) -> &[u8] {
        &self.mem[..self.end]
    }

    /// Get length of block
    pub fn len(&self) -> usize {
        self.mem[..self.end].len()
    }

    /// Return true if the block is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct Record<'a> {
    pub comment: &'a [u8],
    pub sequence: &'a [u8],
    pub plus: &'a [u8],
    pub quality: &'a [u8],
}

#[macro_export(local_inner_macros)]
macro_rules! impl_producer {
    ($name:ident, $correct_block_size:expr) => {
        pub struct $name {
            offset: u64,
            blocksize: u64,
            file: std::fs::File,
            file_length: u64,
        }

        impl $name {
            /// Create a new Block producer
            #[inline(always)]
            pub fn new<P>(path: P) -> error::Result<Self>
            where
                P: AsRef<std::path::Path>,
                Self: Sized,
            {
                Self::with_blocksize(crate::DEFAULT_BLOCKSIZE, path)
            }

            pub fn with_blocksize<P>(blocksize: u64, path: P) -> error::Result<Self>
            where
                P: AsRef<std::path::Path>,
            {
                Ok(Self {
                    offset: 0,
                    blocksize: Self::fix_blocksize::<P>(&path, blocksize)?,
                    file_length: Self::filesize::<P>(&path)?,
                    file: std::fs::File::open(path)
                        .map_err(|source| error::Error::OpenFile { source })?,
                })
            }

            pub fn next_block(&mut self) -> error::Result<Option<block::Block>> {
                if self.offset() == self.file_length() {
                    Ok(None)
                } else if self.offset() + self.blocksize() >= self.file_length() {
                    let block = unsafe {
                        memmap::MmapOptions::new()
                            .offset(self.offset())
                            .len((self.file_length() - self.offset()) as usize)
                            .map(self.file())
                            .map_err(|source| error::Error::MapFile { source })?
                    };

                    self.set_offset(self.file_length());

                    Ok(Some(block::Block::new(block.len(), block)))
                } else {
                    let block = unsafe {
                        memmap::MmapOptions::new()
                            .offset(self.offset())
                            .len(self.blocksize() as usize)
                            .map(self.file())
                            .map_err(|source| error::Error::MapFile { source })?
                    };

                    let blocksize = Self::correct_block_size(&block)?;
                    self.set_offset(self.offset() + blocksize);
                    Ok(Some(block::Block::new(blocksize as usize, block)))
                }
            }

            /// Get file size
            pub fn filesize<P>(path: &P) -> error::Result<u64>
            where
                P: AsRef<std::path::Path>,
            {
                Ok(path
                    .as_ref()
                    .metadata()
                    .map_err(|source| error::Error::MetaDataFile { source })?
                    .len())
            }

            /// Fix blocksize
            pub fn fix_blocksize<P>(path: &P, blocksize: u64) -> error::Result<u64>
            where
                P: AsRef<std::path::Path>,
                Self: Sized,
            {
                Ok(Self::filesize::<P>(path)?.min(blocksize))
            }

            /// Search the begin of the partial record at the end of [Block](Block)
            #[inline(always)]
            pub fn correct_block_size(block: &[u8]) -> error::Result<u64> {
                $correct_block_size(block)
            }

            /// Get current value of offset
            pub fn offset(&self) -> u64 {
                self.offset
            }

            /// Get file length
            pub fn file_length(&self) -> u64 {
                self.file_length
            }

            /// Get file
            pub fn file(&self) -> &std::fs::File {
                &self.file
            }

            /// Get blocksize
            pub fn blocksize(&self) -> u64 {
                self.blocksize
            }

            /// Set value of offset
            pub fn set_offset(&mut self, value: u64) {
                self.offset = value;
            }
        }

        impl Iterator for $name {
            type Item = error::Result<block::Block>;

            fn next(&mut self) -> Option<Self::Item> {
                match self.next_block() {
                    Ok(Some(block)) => Some(Ok(block)),
                    Ok(None) => None,
                    Err(e) => Some(Err(e)),
                }
            }
        }
    };
}

#[macro_export(local_inner_macros)]
macro_rules! impl_reader {
    ($name:ident, $next_record:expr) => {
        pub struct $name {
            offset: usize,
            block: block::Block,
        }

        impl $name {
            pub fn new(block: block::Block) -> Self {
                Reader { offset: 0, block }
            }

            #[inline(always)]
            pub fn next_record<'a>(&'a mut self) -> error::Result<Option<block::Record<'a>>> {
                $next_record(&mut self.block, &mut self.offset)
            }

            pub fn get_line(
                block: &block::Block,
                offset: &usize,
            ) -> error::Result<std::ops::Range<usize>> {
                let next = block.data()[*offset..]
                    .find_byte(b'\n')
                    .ok_or(error::Error::PartialRecord)?;
                let range = *offset..*offset + next;

                Ok(range)
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block() {
        let file = crate::tests::generate_fastq(42, 1_000, 50);

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
}
