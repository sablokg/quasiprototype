mod args;
use crate::args::CommandParse;
use crate::args::Commands;
use clap::Parser;
use figlet_rs::FIGfont;
mod allele_linkage;
mod blast;
mod samtools;

/*
Gaurav Sablok
gsablok@proton.me
*/

fn main() {
    let fontgenerate = FIGfont::standard().unwrap();
    let repgenerate = fontgenerate.convert("quasiprototype");
    println!("{}", repgenerate.unwrap());
    let argparse = CommandParse::parse();
    match &argparse.command {
        Commands::Classifier {
            filepathinput,
            qualityinput,
            genotypeinput,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command_logistic =
                    logisticclassification(filepathinput, qualityinput, genotypeinput).unwrap();
                let command_knn =
                    knnclassification(filepathinput, qualityinput, genotypeinput).unwrap();
                let command_random =
                    randomforestclassification(filepathinput, qualityinput, genotypeinput).unwrap();
                println!(
                    "The commands have been finished:{}\n{}\n{}\n",
                    command_logistic, command_knn, command_random
                );
            });
        }
    }
}
