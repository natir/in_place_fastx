/* mod declaration */
pub mod block;

/* type declaration */
/// Record store a fasta record all field is public
pub struct Record<'a> {
    pub comment: &'a [u8],
    pub sequence: &'a [u8],
}
