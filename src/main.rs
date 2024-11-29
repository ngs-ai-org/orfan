use clap::{arg, command, Arg, Command, ArgAction, value_parser};

mod orf_seeker;
mod utils;
mod common;
mod transcript_mutator;
mod protein_fragmentor;


extern crate pretty_env_logger;

fn main()
{
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    pretty_env_logger::init_timed();
    
    log::trace!("Running Options");

    let matches = Command::new("ORFan")
        .version("0.0.1")
        .subcommand_required(true)
        .propagate_version(true)
        .about(
        "
        Open Reading Frame analysis (ORFan) is a software that allows to perform various operations. 
        It can be used to convert transcript sequences to ORFs. It can be used to mutate transcript sequences
        using a vcf file containing mutations overlapping exons. You can use it to perform in-silico fragmentation.
        You can also perform approximate peptide mass calculation for your peptides")
        .subcommand(
            Command::new("transcript-mutator")
                .about("Mutated transcript sequences")
                .arg(
                    Arg::new("vcf_file")
                        .short('f')
                        .long("vcf-file")
                        .help("input vcf file containing mutations")
                        .required(true)
                )
                .arg(
                    Arg::new("gtf_file")
                        .short('g')
                        .long("gtf-file")
                        .help("path to gtf annotation file")
                        .required(true)
                )
                .arg(
                    Arg::new("reference_genome")
                        .short('r')
                        .long("reference-genome")
                        .help("reference genome in fasta format")
                        .required(true)
                )
                .arg(
                    Arg::new("output_file")
                        .short('o')
                        .long("output-file")
                        .help("Path to output file")
                        .required(true)
                )
                .arg(
                    Arg::new("strand_aware_off")
                        .short('s')
                        .long("strand-aware-off")
                        .help("Whether to turn OFF strand aware sequence report. NOT RECOMMENDED TO BE TURNED OFF if the transcripts are going to be used to generate ORFs")
                        .required(false)
                        .action(ArgAction::SetFalse)
                )
        )
        .subcommand(
            Command::new("orf-seeker")
            .about("Reads a fasta file containing transcripts and creates all possible ORFs")
            .arg(
                Arg::new("input_fasta_file")
                    .short('i')
                    .long("input-file")
                    .help("Input fasta file")
                    .required(true)
            )
            .arg(
                Arg::new("normal_fasta_file")
                    .short('n')
                    .long("normal-fasta-file")
                    .help("Fasta file containing transcript sequences as obtained by seqtk")
                    .required(true)
            )
            .arg(
                Arg::new("output_fasta_file")
                    .short('o')
                    .long("output-file")
                    .help("Output fasta file")
                    .required(true)
            )
            .arg(
                Arg::new("mutated_only")
                    .short('m')
                    .long("mutated-only")
                    .help("Whether to only write ORFs containing mutation")
                    .required(false)
                    .action(ArgAction::SetTrue)
            )
            .arg(
                Arg::new("min_orf_len")
                    .short('l')
                    .long("min-len")
                    .help("Minimum length of the ORF")
                    .required(false)
                    .value_parser(clap::value_parser!(usize))
                    .default_value("6")
            )
            .arg(
                Arg::new("all_codons")
                    .short('a')
                    .long("all-codons")
                    .help("Whether to print ORFs starting with any codon")
                    .required(false)
                    .action(ArgAction::SetTrue)
            )
            .arg(
                Arg::new("write_missing")
                    .short('w')
                    .long("write-missing")
                    .help("Whether to write header line of missing ORFs in output")
                    .required(false)
                    .action(ArgAction::SetTrue)
            )
        )
        .subcommand(
            Command::new("fragment")
                .about("Fragments peptides")
                .arg(
                    Arg::new("fasta_file")
                        .short('f')
                        .long("fasta-file")
                        .help("Path to the input fasta file")
                        .required(true)
                )
                .arg(
                    Arg::new("output_file")
                        .short('o')
                        .long("output-file")
                        .help("Output file to write results")
                        .required(true)
                )
                .arg(
                    Arg::new("enzyme")
                        .short('e')
                        .long("enzyme")
                        .help("Enzyme schema to use for in silico fragmentation")
                        .required(false)
                        .default_value("trypsin")
                )
                .arg(
                    Arg::new("missed_cleavages")
                        .short('c')
                        .long("missed-cleavages")
                        .help("Simulate n number of missed cleavages")
                        .required(false)
                        .value_parser(value_parser!(usize))
                        .default_value("0")
                )
        )
        .get_matches();
    
    log::trace!("Running matching of subcommands");
    match matches.subcommand() {
        Some (("transcript-mutator", sub_matches)) =>
        {
            log::trace!("transcript-mutator subcommand.");
            let vcf_file: String = sub_matches.get_one::<String>("vcf_file").expect("Cannot read vcf_file parameter").to_string();
            let gtf_file: String = sub_matches.get_one::<String>("gtf_file").expect("Cannot read gtf_file parameter").to_string();
            let reference_genome: String = sub_matches.get_one::<String>("reference_genome").expect("Cannot read reference genome parameter").to_string();
            let output_file: String = sub_matches.get_one::<String>("output_file").expect("Cannot read output file parameter").to_string();
            let strand_aware_off: bool = sub_matches.get_flag("strand_aware_off");
            
            log::trace!("Running transcript_mutator");
            crate::transcript_mutator::transcript_mutator::transcript_mutator_main(vcf_file.as_str(), gtf_file.as_str(), reference_genome.as_str(), strand_aware_off, output_file.as_str())
        }
        Some (("orf-seeker", sub_matches)) => 
        {
            log::trace!("orf-seeker subcommand.");
            let input_fasta_file: String = sub_matches.get_one::<String>("input_fasta_file").expect("Cannot read parameter").to_string();
            let normal_fasta_file: String = sub_matches.get_one::<String>("normal_fasta_file").expect("Cannot read normal fasta file").to_string();
            let output_fasta_file: String = sub_matches.get_one::<String>("output_fasta_file").expect("Cannot read output fasta file").to_string();
            let mutated_only: bool = sub_matches.get_flag("mutated_only");
            let all_codons: bool = sub_matches.get_flag("all_codons");
            let write_missing: bool = sub_matches.get_flag("write_missing");
            let min_len: usize = *(sub_matches.get_one::<usize>("min_orf_len").expect("Cannot read minimum orf len parameter"));

            log::trace!("Running orf_seeker_main");
            crate::orf_seeker::orf_seeker_main(input_fasta_file.as_str(), normal_fasta_file.as_str(), output_fasta_file.as_str(), min_len, write_missing, mutated_only, all_codons);
        }
        Some(("fragment", sub_matches)) => 
        {
            log::trace!("fragment subcommand.");
            let fasta_file: String = sub_matches.get_one::<String>("fasta_file").expect("Cannot read parameter").to_string();
            let output_file: String = sub_matches.get_one::<String>("output_file").expect("Cannot read specified output_file").to_string();
            let enzyme: String = sub_matches.get_one::<String>("enzyme").expect("Cannot read specified enzyme").to_string();
            let missed_cleavage: usize = *(sub_matches.get_one::<usize>("missed_cleavages").expect("cannot read missed cleavages"));

            log::trace!("Running protein_fragmentation_main");
            crate::protein_fragmentor::protein_fragmentation_main(fasta_file, enzyme.to_lowercase(), output_file.to_string(), missed_cleavage);
        }
        _ => {
            log::trace!("No matches for subcommands");
            panic!("Hello World");
        }
    }
}





// This is required for unit-testing. I've added it at the end of the main on purpose
#[cfg(test)]
mod test;

