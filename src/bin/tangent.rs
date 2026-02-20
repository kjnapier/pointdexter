use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use clap::Parser;

pub fn main() ->  Result<(), Box<dyn std::error::Error>> {
    // We need to make this more flexible.
    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    Ok(())

    

}