/* mod declaration */
pub mod block;
pub mod parser;

/* type declaration */
/// Record store a fastq record, all field is public
pub struct Record<'a> {
    pub comment: &'a [u8],
    pub sequence: &'a [u8],
    pub plus: &'a [u8],
    pub quality: &'a [u8],
}
