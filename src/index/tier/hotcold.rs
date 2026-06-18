#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HotMode {
    Skip,
    Roaring,
    Bp128,
}

impl HotMode {
    pub fn as_str(self) -> &'static str {
        match self {
            HotMode::Skip => "skip",
            HotMode::Roaring => "roaring",
            HotMode::Bp128 => "bp128",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TierAssignment {
    Hot,
    Cold,
}

pub fn assign_tiers(dfs: &[u32], hot_df_threshold: u32) -> Vec<TierAssignment> {
    dfs.iter()
        .map(|&df| {
            if df >= hot_df_threshold {
                TierAssignment::Hot
            } else {
                TierAssignment::Cold
            }
        })
        .collect()
}
