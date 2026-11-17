mod args;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
use figlet_rs::FIGfont;
mod allele_linkage;
mod blast;
mod samtools;
use allele_linkage::build_allele_graph;
use allele_linkage::export_haplotypes_to_fasta;
use allele_linkage::load_reads_from_bam;
use blast::blastquasi;
use samtools::samalignment;

/*
Gaurav Sablok
gsablok@proton.me
*/

fn main() {
    let fontgenerate = FIGfont::standard().unwrap();
    let repgenerate = fontgenerate.convert("quasi");
    println!("{}", repgenerate.unwrap());
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::Quasispecies {
            fastafile,
            fastqfile,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let _ = samalignment(fastafile, fastqfile);
                let bam_path = "aligned.bam";
                let vcf_path = "variants.vcf.gz";
                let ref_fa = fastafile;
                let out_fa = "phased_haplotypes.fa";
                let reads = load_reads_from_bam(bam_path, vcf_path, 20, 20).expect("failed");
                let graph = build_allele_graph(&reads, 0.2, 5);
                export_haplotypes_to_fasta(&graph, ref_fa, out_fa).expect("failed");
                blastquasi(out_fa).expect("failed");
            });
        }
    }
}
