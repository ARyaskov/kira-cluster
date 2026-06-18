pub fn idf(n_seqs: u64, df: u32) -> f32 {
    (((n_seqs as f32) + 1.0) / ((df as f32) + 1.0)).ln()
}
