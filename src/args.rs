use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "Quasispecies",
    version = "1.0",
    about = "quasiprototype- prototyping quasispecies.
       ************************************************
       Author Gaurav Sablok,
       Email: gsablok@proton.me
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Quasispecies Prototyping
    Quasispecies {
        /// fastafile
        fastafile: String,
        /// fastq file
        fastqfile: String,
        ///
        thread: String,
    },
}
