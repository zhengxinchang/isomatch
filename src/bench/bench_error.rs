use thiserror::Error;

#[derive(Error, Debug)]
pub enum BenchError {
    #[error("Malformed GTF: {reason}")]
    MalformedGTF { reason: String },

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    ParseIntError(#[from] std::num::ParseIntError),

    #[error(transparent)]
    LibGtf(#[from] libgtf::error::Error),

    #[error(transparent)]
    Tool(#[from] crate::tools::tools_error::ToolError),

    #[error("Can not open merged GTF: {path}")]
    OpenGTF {
        path: std::path::PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("Can not read merged GTF near line {line_no}")]
    ReadGTF {
        line_no: usize,
        #[source]
        source: Box<BenchError>,
    },
}
