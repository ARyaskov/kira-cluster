#[derive(Debug, Clone, Copy)]
pub struct PlaPoint {
    pub key: u64,
    pub pos: u32,
}

pub fn sample_points(keys: &[u64], stride: usize) -> Vec<PlaPoint> {
    let step = stride.max(1);
    let mut out = Vec::new();
    for i in (0..keys.len()).step_by(step) {
        out.push(PlaPoint {
            key: keys[i],
            pos: i as u32,
        });
    }
    if let Some((&last, idx)) = keys.last().map(|v| (v, keys.len() - 1)) {
        if out.last().map(|p| p.key) != Some(last) {
            out.push(PlaPoint {
                key: last,
                pos: idx as u32,
            });
        }
    }
    out
}
