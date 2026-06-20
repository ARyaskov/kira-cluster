#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CascadeLevel {
    Fast,
    Medium,
    Sensitive,
    Ultra,
}

impl CascadeLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            CascadeLevel::Fast => "fast",
            CascadeLevel::Medium => "medium",
            CascadeLevel::Sensitive => "sensitive",
            CascadeLevel::Ultra => "ultra",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CascadeThresholds {
    /// Loosened identity used by the positional prefilter (alignment can recover
    /// identity, so the prefilter is permissive).
    pub prefilter_identity: f32,
    pub final_identity: f32,
    pub final_coverage: f32,
    pub run_gapped: bool,
}

pub fn thresholds_for(
    level: CascadeLevel,
    min_identity: f32,
    min_coverage: f32,
    sensitivity: f32,
) -> CascadeThresholds {
    let s = sensitivity.clamp(0.0, 10.0);
    let loose = 0.04 * s;

    // Note: the coverage prefilter uses `final_coverage` directly (a length-based
    // necessary condition), so there is no separate loosened prefilter coverage.
    match level {
        CascadeLevel::Fast => CascadeThresholds {
            prefilter_identity: (min_identity - loose).clamp(0.0, 1.0),
            final_identity: min_identity,
            final_coverage: min_coverage,
            run_gapped: false,
        },
        CascadeLevel::Medium => CascadeThresholds {
            prefilter_identity: (min_identity - (0.15 + loose)).clamp(0.0, 1.0),
            final_identity: min_identity,
            final_coverage: min_coverage,
            run_gapped: true,
        },
        CascadeLevel::Sensitive => CascadeThresholds {
            prefilter_identity: (min_identity - (0.30 + loose)).clamp(0.0, 1.0),
            final_identity: min_identity,
            final_coverage: min_coverage,
            run_gapped: true,
        },
        CascadeLevel::Ultra => CascadeThresholds {
            prefilter_identity: (min_identity - (0.40 + loose)).clamp(0.0, 1.0),
            final_identity: min_identity,
            final_coverage: (min_coverage * 0.9).clamp(0.0, 1.0),
            run_gapped: true,
        },
    }
}
