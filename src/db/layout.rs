use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct DbMeta {
    pub version: u32,
    pub dbtype: String,
    pub n_seqs: u64,
    pub total_bases: u64,
    pub max_seq_len: u32,
    pub encoding: String,
    pub endian: String,
}
