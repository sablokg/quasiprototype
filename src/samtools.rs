use std::error::Error;
use std::process::Command;

/*
Gaurav Sablok
gsablok@proton.me
*/

pub fn samalignment(pathfasta: &str, pathfastq: &str) -> Result<String, Box<dyn Error>> {
    let _ = Command::new("samtools")
        .arg("faidx")
        .arg(pathfasta)
        .output()
        .expect("command failed");
    let _ = Command::new("minimap2")
        .arg("-ax")
        .arg("map-ont")
        .arg(pathfasta)
        .arg(pathfastq)
        .output()
        .expect("command failed");
    let _ = Command::new("samtools")
        .arg("sort")
        .arg("-o")
        .arg("aligned.bam")
        .output()
        .expect("command failed");
    let _ = Command::new("samtools").arg("index").arg("aligned.bam");

    let _ = Command::new("bcftools")
        .arg("mpileup")
        .arg("-Ou")
        .arg("-f")
        .arg(pathfasta)
        .arg("--max-depth")
        .arg("5000")
        .arg("aligned.bam")
        .arg("|")
        .arg("bcftools")
        .arg("call")
        .arg("-mv")
        .arg("-Oz")
        .arg("-o")
        .arg("variant.vcf.gz")
        .output()
        .expect("failed");

    let _ = Command::new("bftools").arg("index").arg("variant.vcf.gz");
    Ok("Sam alignment have been parsed".to_string())
}
