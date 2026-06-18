use std::fs;
use std::path::Path;

use crate::cluster::SeqDb;
use crate::cluster::greedy::ClusterResult;
use crate::error::{AppError, Result};
use crate::io::atomic::write_atomic;

pub fn write_outputs(result_dir: &Path, db: &SeqDb, result: &ClusterResult) -> Result<()> {
    fs::create_dir_all(result_dir)
        .map_err(|e| AppError::io(format!("create dir {}", result_dir.display()), e))?;

    let mut by_cluster: Vec<Vec<u32>> = vec![Vec::new(); result.representatives.len()];
    for (seq_id, &cid) in result.assignments.iter().enumerate() {
        by_cluster[cid as usize].push(seq_id as u32);
    }
    for members in &mut by_cluster {
        members.sort_unstable();
    }

    let mut tsv = Vec::new();
    for (cid, members) in by_cluster.iter().enumerate() {
        let rep = result.representatives[cid];
        for &member in members {
            tsv.extend_from_slice(format!("{cid}\t{rep}\t{member}\n").as_bytes());
        }
    }
    write_atomic(&result_dir.join("clusters.tsv"), &tsv)?;

    let mut fasta = Vec::new();
    for (cid, &rep_id) in result.representatives.iter().enumerate() {
        let _ = cid;
        fasta.extend_from_slice(b">");
        fasta.extend_from_slice(db.name(rep_id).as_bytes());
        fasta.extend_from_slice(b"\n");
        fasta.extend_from_slice(db.seq(rep_id));
        fasta.extend_from_slice(b"\n");
    }
    write_atomic(&result_dir.join("rep_seqs.fasta"), &fasta)?;

    Ok(())
}
