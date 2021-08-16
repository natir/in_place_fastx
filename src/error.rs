#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("FastMap failled to read file metadata {source}")]
    MetaDataFile { source: std::io::Error },

    #[error("FastMap can't open file {source}")]
    OpenFile { source: std::io::Error },

    #[error("FastMap can't map file on memory {source}")]
    MapFile { source: std::io::Error },

    #[error("FastMap didn't find new line in block increase block size")]
    NoNewLineInBlock,

    #[error("Input file seems not be a fastq file")]
    NotAFastqFile,

    #[error("FastMap found a partial record")]
    PartialRecord,
}

pub type Result<T> = std::result::Result<T, Error>;
