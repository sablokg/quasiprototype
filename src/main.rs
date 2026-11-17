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
    match &argparse.command {}
}
