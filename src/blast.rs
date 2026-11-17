use std::error::Error;
use std::process::Command;

/*
Gaurav Sablok
gsablok@proton.me
*/

pub fn blastquasi(fasta: &str) -> Result<(), Box<dyn Error>> {
    let _ = Command::new("blastn")
        .arg("-query")
        .arg(fasta)
        .arg("-db")
        .arg("nt")
        .arg("remote")
        .arg("-max_target_seqs")
        .arg("3")
        .arg("perc_identity")
        .arg("97")
        .arg("outfmt")
        .arg("6 qseqid sseqid sscinames pident length evalue")
        .arg("-out")
        .arg("species_classification.tsv")
        .output()
        .expect("failed");

    Ok(())
}
