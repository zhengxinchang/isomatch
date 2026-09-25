pub(crate) mod attributes_index;
pub(crate) mod builder;
pub(crate) mod format;
pub(crate) mod junction_pool;
pub(crate) mod reader;
pub(crate) mod splice_site_pool;
pub(crate) mod string_pool;
pub(crate) mod tx;

pub const ISOMX_VERSION: u32 = 4;
pub const ISOM_GTF_SCHEMA: &str = "1";
pub const ISOMS_VERSION: u32 = 3;

pub(crate) use attributes_index::AttrIndexBuilder;
pub use attributes_index::{AttrIndexHeader, AttrIndexReader};
pub(crate) use builder::IndexBuilder;
pub(crate) use format::ChromBlockBuilder;
pub use format::{ChromDirectoryEntry, IndexHeader};
pub use junction_pool::{JunctionPool, JunctionSpan};
pub use reader::{ChromBlockReader, IndexReader};
pub use splice_site_pool::{SpliceSitePair, SpliceSitePool, SpliceSiteSpan};
pub use string_pool::{StringPool, StringSpan};
pub use tx::{TxBase, TxBaseFlags, TxBoundary};
