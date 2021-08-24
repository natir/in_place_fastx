/* mod declaration */
pub mod block;
pub mod parser;

/* type declaratino*/
/// Record store a fastq record in public field
pub struct Record<'a> {
    pub comment: &'a [u8],
    pub sequence: &'a [u8],
    pub plus: &'a [u8],
    pub quality: &'a [u8],
}

/// Block reperesent a section of file memory mapped in file
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
