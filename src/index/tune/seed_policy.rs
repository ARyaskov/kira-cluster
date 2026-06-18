use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;
use crate::io::mmap::map_readonly;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SeedPolicyMode {
    FixedK,
    VariableK,
    Syncmer,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedPolicy {
    pub mode: SeedPolicyMode,
    pub k_values: Vec<u32>,
    pub m: u32,
    pub rules: SeedPolicyRules,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeedPolicyRules {
    pub short_query_use_k: u32,
    pub long_query_use_k: u32,
}

impl SeedPolicy {
    pub fn from_mode_str(mode: &str) -> Option<Self> {
        let mut policy = choose_default_policy();
        policy.mode = match mode {
            "fixed_k" => SeedPolicyMode::FixedK,
            "variable_k" => SeedPolicyMode::VariableK,
            "syncmer" => SeedPolicyMode::Syncmer,
            _ => return None,
        };
        Some(policy)
    }
}

pub fn choose_default_policy() -> SeedPolicy {
    SeedPolicy {
        mode: SeedPolicyMode::FixedK,
        k_values: vec![13, 14, 15],
        m: 20,
        rules: SeedPolicyRules {
            short_query_use_k: 13,
            long_query_use_k: 15,
        },
    }
}

pub fn write_seed_policy(indexdir: &Path, policy: &SeedPolicy) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(policy)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize seed policy: {e}")))?;
    write_atomic(&indexdir.join("seed_policy.json"), &bytes)
}

pub fn load_seed_policy(indexdir: &Path) -> Result<Option<SeedPolicy>> {
    let path = indexdir.join("seed_policy.json");
    if !path.exists() {
        return Ok(None);
    }
    let mm = map_readonly(&path)?;
    let policy: SeedPolicy = serde_json::from_slice(&mm)
        .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse {}: {e}", path.display())))?;
    Ok(Some(policy))
}

pub fn choose_k_for_query(policy: &SeedPolicy, default_k: usize, query: &[u8]) -> usize {
    match policy.mode {
        SeedPolicyMode::FixedK => default_k,
        SeedPolicyMode::VariableK => {
            if query.len() < 80 {
                policy.rules.short_query_use_k as usize
            } else {
                policy.rules.long_query_use_k as usize
            }
        }
        SeedPolicyMode::Syncmer => default_k,
    }
}
