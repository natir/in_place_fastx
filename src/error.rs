#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("in_place_fastx failled to read file metadata {source}")]
    MetaDataFile { source: std::io::Error },

    #[error("in_place_fastx can't open file {source}")]
    OpenFile { source: std::io::Error },

    #[error("in_place_fastx can't map file on memory {source}")]
    MapFile { source: std::io::Error },

    #[error("in_place_fastx didn't find new line in block increase block size")]
    NoNewLineInBlock,

    #[error("Input file seems not be a fastq file")]
    NotAFastqFile,

    #[error("Input file seems not be a fasta file")]
    NotAFastaFile,

    #[error("in_place_fastx found a partial record")]
    PartialRecord,
}

pub type Result<T> = std::result::Result<T, Error>;
