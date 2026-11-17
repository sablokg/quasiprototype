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

    let _ = Command::new
    Ok("Sam alignment have been parsed".to_string())
}


/*
# 2. Call polymorphic sites across the mixed sample
# Note: --ploidy 2 (or higher) ensures mixed alleles are captured in the VCF
bcftools mpileup -Ou -f $REF --max-depth 5000 aligned.bam | \
    bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz
 */
