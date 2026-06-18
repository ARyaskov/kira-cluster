#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SketchSeed {
    pub hash: u64,
    pub seq_id: u32,
}

pub fn extract_sketch_seeds(_seq: &[u8], _seq_id: u32, _n: usize) -> Vec<SketchSeed> {
    Vec::new()
}
