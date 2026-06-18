use std::path::Path;

use crate::error::Result;
use crate::index::manifest::load_manifest;
use crate::index::segment::{SegmentHandle, segment_dir};

pub fn print_index_info(indexdir: &Path) -> Result<()> {
    let manifest = load_manifest(indexdir)?;
    println!("index.version={}", manifest.version);
    println!("index.segments={}", manifest.segments.len());

    let mut total = 0u64;
    for seg in manifest.segments {
        let dir = segment_dir(indexdir, seg.segment_id);
        let handle = SegmentHandle::open(indexdir, seg.segment_id)?;
        let mut seg_size = 0u64;
        for name in [
            "keys.u64",
            "df.u32",
            "offsets_hot.u64",
            "offsets_cold.u64",
            "postings_hot.bin",
            "postings_cold.bin",
            "bloom.bin",
            "postings.u32",
            "offsets.u64",
        ] {
            let p = dir.join(name);
            if let Ok(meta) = std::fs::metadata(&p) {
                seg_size += meta.len();
            }
        }
        total += seg_size;

        println!("segment.id={}", seg.segment_id);
        println!("segment.keys={}", handle.keys().len());
        println!("segment.codec_hot={}", handle.codec_hot_name());
        println!("segment.codec_cold={}", handle.codec_cold_name());
        println!("segment.bytes={}", seg_size);
    }

    println!("index.bytes_total={}", total);
    Ok(())
}
