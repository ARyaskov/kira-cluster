use std::io::Write;
use std::path::Path;

use crate::cluster::SeqDb;
use crate::cluster::greedy::ClusterResult;
use crate::error::{AppError, Result};
use crate::io::atomic::write_atomic_with;

/// Number of residues per line in the representative FASTA output.
const FASTA_WRAP: usize = 60;

pub fn write_outputs(result_dir: &Path, db: &SeqDb, result: &ClusterResult) -> Result<()> {
    let mut by_cluster: Vec<Vec<u32>> = vec![Vec::new(); result.representatives.len()];
    for (seq_id, &cid) in result.assignments.iter().enumerate() {
        by_cluster[cid as usize].push(seq_id as u32);
    }
    for members in &mut by_cluster {
        members.sort_unstable();
    }

    write_atomic_with(&result_dir.join("clusters.tsv"), |w| {
        let mut line = String::new();
        for (cid, members) in by_cluster.iter().enumerate() {
            let rep = result.representatives[cid];
            for &member in members {
                line.clear();
                use std::fmt::Write as _;
                let _ = writeln!(line, "{cid}\t{rep}\t{member}");
                w.write_all(line.as_bytes())
                    .map_err(|e| AppError::io("write clusters.tsv", e))?;
            }
        }
        Ok(())
    })?;

    write_atomic_with(&result_dir.join("rep_seqs.fasta"), |w| {
        for &rep_id in &result.representatives {
            w.write_all(b">")
                .and_then(|_| w.write_all(db.name(rep_id).as_bytes()))
                .and_then(|_| w.write_all(b"\n"))
                .map_err(|e| AppError::io("write rep_seqs.fasta", e))?;
            let seq = db.seq(rep_id);
            for chunk in seq.chunks(FASTA_WRAP) {
                w.write_all(chunk)
                    .and_then(|_| w.write_all(b"\n"))
                    .map_err(|e| AppError::io("write rep_seqs.fasta", e))?;
            }
        }
        Ok(())
    })?;

    Ok(())
}
