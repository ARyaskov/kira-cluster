use std::path::Path;

use crate::error::Result;
use crate::index::manifest::Manifest;
use crate::index::search::{Hit, SearchOpts, search_query_sequence};
use crate::index::tune::seed_policy::{SeedPolicy, load_seed_policy};

pub mod bloom;
pub mod build;
pub mod codec;
pub mod df;
pub mod index_info;
pub mod learned;
pub mod manifest;
pub mod scoring;
pub mod search;
pub mod segment;
pub mod sketch;
pub mod tier;
pub mod tune;
pub mod update;

pub struct IndexHandle {
    pub manifest: Manifest,
    pub segments: Vec<segment::SegmentHandle>,
    pub seed_policy: Option<SeedPolicy>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResolvedSeqId<'a> {
    pub global_id: u64,
    pub segment_id: u32,
    pub local_id: u32,
    pub name: &'a str,
}

impl IndexHandle {
    pub fn open(indexdir: &Path) -> Result<Self> {
        let manifest = manifest::load_manifest(indexdir)?;
        let mut segments = Vec::with_capacity(manifest.segments.len());
        for seg in &manifest.segments {
            segments.push(segment::SegmentHandle::open(indexdir, seg.segment_id)?);
        }
        let seed_policy = load_seed_policy(indexdir)?;
        Ok(Self {
            manifest,
            segments,
            seed_policy,
        })
    }

    pub fn search(&self, query: &[u8], opts: &SearchOpts) -> Vec<Hit> {
        let mut hits = Vec::new();
        search_query_sequence(&self.segments, query, opts, &mut hits);
        hits
    }

    pub fn resolve_global_seq_id(&self, global_id: u64) -> Option<ResolvedSeqId<'_>> {
        self.segments
            .iter()
            .find_map(|seg| seg.resolve_global_seq_id(global_id))
            .map(|seq| ResolvedSeqId {
                global_id: seq.global_id,
                segment_id: seq.segment_id,
                local_id: seq.local_id,
                name: seq.name,
            })
    }
}
