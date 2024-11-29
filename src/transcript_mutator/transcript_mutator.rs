use std::{env, str};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead,BufReader};
use std::path::Path;
use chrono::Local;
use rustc_hash::FxHashMap;
use core::panic;
use faimm::IndexedFasta;
use rust_htslib::bcf::{Reader, Read, Record};

extern crate pretty_env_logger;




pub struct MutatedTranscript 
{
    pub mutation_pos_on_sequence: usize,
    pub mutated_sequence: String
}



pub fn check_exon_ordering(exon_vector: &Vec < crate::common::gtf_parsing::ExonStruct >) -> bool
{
    let vec_len: usize = exon_vector.len();
    if vec_len == 1
    {
        return true;
    }
    else if exon_vector[0].start > exon_vector[1].start
    {
        return false; 
    }
    else
    {
        return true;   
    }
}


///
/// Build transcripts is used to reconstruct the transcript sequences
/// from a hash containing for each transcript the exon information.
/// 
/// # Arguments 
/// * transcript_exon_hash - HashMap where the key is the transcript id and values a vector of exon information
/// * fasta_file           - the reference genome fasta file
/// 
/// # Examples
/// 
/// ```
/// let gtf_test_file: &str = "./data/test.gtf";
/// let gtf_fasta_test_file: &str = "./data/test.gtf.fa";
/// let reference_genome = "/gpfs/omics/resources/Public/Genomes/hg37_p13/GRCh37.p13.genome.fa";
/// let gtf_hash = partial_read_gtf_file(&gtf_test_file);
/// let transcript_hash = build_transcripts(&gtf_hash, reference_genome);
/// ```
pub fn build_transcripts(transcript_exon_hash: &mut FxHashMap < String, Vec < crate::common::gtf_parsing::ExonStruct >>, fasta_file: &str, transcript_ids:  &FxHashMap < String, Vec < crate::common::bcf_parsing::MutBedStruct >> ) -> FxHashMap< String, String > 
{
    let mut trx_seq_hash: FxHashMap<String,String> = FxHashMap::default();

    let mut transcript_count = 0;

    for (trx_id, exon_vector) in transcript_exon_hash.clone().iter()
    {
        if ! transcript_ids.contains_key(trx_id) { continue; }
        if transcript_count % 1000 == 0 { log::info!("Processed {} transcripts", transcript_count); }
        transcript_count += 1;


        let mut transcript_sequence: String = String::new();
        
        
        let mut exons_tmp: Vec < crate::common::gtf_parsing::ExonStruct > = exon_vector.clone();


        //exons.sort_by(|a, b| a.start.cmp(&b.start));
        
        /**
        ** The different steps to deal with hg38 fuckeries of gtf
        *! Transcripts on the positive strand should be in correct orientation
        *! IF seq_starts go from bigger to smaller then there is a problem with the gtf file

        *? Is transcript on negative strand? 
        **     * If yes, then: 
        **         * Create new vector of exons with proper seq_starts in the correct orientation
        **         * Modify the exon_hash value for this transcript with new exon vector. 
        *? Is transcript on positive strand? 
        *! Make sure that seq_starts and order of exon genomic positions is correct. IF NOT then there
        *! is a problem with the GTF file and ORFan should not have to deal with that!!!!
        */
        
        let mut new_exon_vec: Vec < crate::common::gtf_parsing::ExonStruct > = Vec::new();
        if exon_vector[0].strand.as_str() == "-"
        {
            if ! check_exon_ordering(&exon_vector)
            {
                exons_tmp.sort_by(|a,b| a.start.cmp(&b.start));
                
                let mut seq_start_tmp: usize = 0;
                let mut seq_end_tmp: usize = 0;
                for exon in exons_tmp
                {
                    let new_exon: crate::common::gtf_parsing::ExonStruct = crate::common::gtf_parsing::ExonStruct::new(exon.chrom, exon.start, exon.end, seq_start_tmp, (seq_start_tmp + exon.exon_len -1), exon.strand, exon.exon_number, exon.ensembl_id, exon.gene_name);
                    seq_start_tmp = seq_start_tmp + exon.exon_len;
                    new_exon_vec.push(new_exon);
                }
                *transcript_exon_hash.get_mut(trx_id).unwrap() = new_exon_vec.clone();
            }
            else
            {
                new_exon_vec = exon_vector.clone();
            }
        }
        else
        {
            new_exon_vec = exon_vector.clone();    
        }

        for exon in new_exon_vec.iter()
        {
            let tmp: String = extract_sequence_from_fasta(fasta_file, exon.chrom.as_str(), exon.start, exon.end);
            transcript_sequence += &tmp;
        }
        
        trx_seq_hash.insert(trx_id.clone(),transcript_sequence);

    }
    return trx_seq_hash;
}



///Function that reads an indexed fasta (samtools faidx) and queries 
/// for a specific region.
/// 
/// # Arguments: 
/// 
/// * `fasta_file` - path to reference genome fasta file (the index must be in the same directory)
/// * `chrom`      - chromosome id (same as in the fasta file)
/// * `start`      - start position of the query
/// * `end`        - end position of the query
/// 
/// # Example 
/// 
/// ```
/// 
/// let fasta_file: &str = "test_file.fasta";
/// let chrom: &str      = "chr1";
/// let start: usize     = 1000;
/// let end: usize       = 1500;
/// 
/// extract_sequence_from_fasta(fasta_file, chrom, start, end);
/// ```
pub fn extract_sequence_from_fasta (fasta_file: &str, chrom:&str, start:usize, end:usize) -> String 
{
    let fa = IndexedFasta::from_file(fasta_file).expect("Error opening fasta file. Index not found");
    let chr_index = fa.fai().tid(format!("chr{}", chrom).as_str()).expect("Cannot find chr in index");
    let v = fa.view(chr_index, start, end).expect("Cannot find region");
    return v.to_string();
}


/// Function that takes a nucleotide sequence and creates chunks of length chunk_size.
/// 
/// # Arguments 
/// 
/// * `sequence_to_chunk` - Nucleotide sequence to be chunked
/// * `chunk_size`        - chunk size 
/// 
/// # Example 
/// 
/// ```
/// let sequence_to_chunk: String = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
/// let chunk_size: usize         = 10;
/// chunk_fasta_sequence(&sequence_to_chunk, chunk_size);
/// ```
pub fn chunk_fasta_sequence(sequence_to_chunk: &String, chunk_size:usize) -> Vec < &str >
{
    let subs: Vec < &str > = sequence_to_chunk.as_bytes()
        .chunks(chunk_size)
        .map(str::from_utf8)
        .collect::<Result<Vec<&str>,_>>()
        .unwrap();
    return subs;
}


/// Write fasta entries to file
/// 
/// # Arguments 
/// 
/// * `output_buffer`       - Output buffer where fasta entries are to be written
/// * `header_string`       - fasta header line
/// * `transcript_sequence` - fasta sequence.
pub fn write_fasta_entry(mut output_buffer: &File, header_string: &String, transcript_sequence: &String)
{
    write!(output_buffer, "{}\n", header_string).expect("Cannot write");
    
    let subs: Vec < &str > = chunk_fasta_sequence(transcript_sequence, 70);
    
    for i in subs
    {
        write!(output_buffer, "{}\n",i).expect("Cannot write to file");
    }
}



/// Function that filters vcf_hash for variants 
/// overlapping with exonic regions 
/// 
/// # Arguments 
/// vcf_hash - a hash map containing as key contigs and value a vector of variants
/// gtf_hash - a hash map containing as key the transcript ids and value a vector of exons.
pub fn filter_vcf_hash(vcf_hash: &FxHashMap < String, Vec < crate::common::bcf_parsing::VcfEntry>>, gtf_hash: &FxHashMap < String, Vec < crate::common::gtf_parsing::ExonStruct>>) -> FxHashMap < String, Vec < crate::common::bcf_parsing::MutBedStruct>>
{
    let mut tmp_map: FxHashMap < String, Vec < crate::common::bcf_parsing::MutBedStruct>> = FxHashMap::default();

    for (trx_id, exon_vector) in gtf_hash.iter()
    {
        
        if !vcf_hash.contains_key(&exon_vector[0].chrom) { continue;}
        let chrom_var: Vec < crate::common::bcf_parsing::VcfEntry > = vcf_hash[&exon_vector[0].chrom].clone();
        for var in chrom_var
        {
            for exon in exon_vector
            {   
                if var.chrom != exon.chrom { break;}
                
                log::debug!("{} {}:{}-{} -> {} {}:{}-{};{}", var.var_id,var.chrom, (var.start -1), var.start, trx_id, exon.chrom, exon.start, exon.end, exon.exon_number);
                if ! rust_tools::utils::genomic_overlap::is_overlap(exon.start, exon.end, (var.start - 1) as usize, var.start as usize) { continue; }
                let mut mutation_type: String = String::new();
                if var.is_snv()
                {
                    log::trace!("{} is a SNV", var.var_id);
                    mutation_type = "SNV".to_string();
                }
                else if var.is_del()
                {
                    log::trace!("{} is a DEL", var.var_id);
                    mutation_type = "DEL".to_string();
                }
                else if var.is_ins()
                {
                    log::trace!("{} is a INS", var.var_id);
                    mutation_type = "INS".to_string();
                }
                else if var.is_mnv()
                {
                    log::trace!("{} is a MNV", var.var_id);
                    mutation_type = "MNV".to_string();
                }
                else
                {
                    log::warn!("Unknown mutation type!");
                    eprintln!("WARNING: Unknown mutation type.");
                }
                let record_struct: crate::common::bcf_parsing::MutBedStruct = crate::common::bcf_parsing::MutBedStruct{
                    chrom: var.chrom.clone(), 
                    start: var.start -1,
                    end: var.start, 
                    mutation_id: var.var_id.clone(),
                    ref_allele: var.ref_allele.clone(),
                    alt_allele: var.alt_allele.clone(),
                    variant: "None".to_string(),
                    mutation_type: mutation_type.clone(),
                    trans_mutation_type: "None".to_string()
                };
                //println!("{:?}", record_struct);
                if ! tmp_map.contains_key(trx_id)
                {
                    let mut tmp_vec: Vec < crate::common::bcf_parsing::MutBedStruct > = Vec::new();
                    tmp_vec.push(record_struct.clone());
                    tmp_map.insert(trx_id.clone(), tmp_vec);
                }
                else 
                {
                    let mut tmp_vec: Vec < crate::common::bcf_parsing::MutBedStruct > = tmp_map[trx_id].clone();
                    tmp_vec.push(record_struct.clone());
                    *tmp_map.get_mut(trx_id).unwrap() = tmp_vec; 
                }
                break;
            }
        }
    }
    return tmp_map;

}


pub fn transcript_mutate_snv(trx_id: &str, mutation: &crate::common::bcf_parsing::MutBedStruct, trx_seq: &String, trx_exons: &Vec < crate::common::gtf_parsing::ExonStruct>) -> MutatedTranscript
{
    let mut mutation_pos_on_sequence: usize = 0;
    let mut mutated_transcript: String = String::new();

    for exon in trx_exons
    {
        if rust_tools::utils::genomic_overlap::is_overlap(exon.start, exon.end, mutation.start as usize, mutation.end as usize)
        {
            log::debug!("{}:{:?} -> {:?}",trx_id,exon, mutation);
            mutation_pos_on_sequence = (mutation.end as usize - exon.start - 1) + exon.seq_start;
            
            mutated_transcript= crate::transcript_mutator::sequence::point_mutation_mutator(&trx_seq, mutation_pos_on_sequence, mutation.alt_allele.as_str(), mutation.ref_allele.as_str());
            if exon.strand == "-".to_string() 
            {
                mutated_transcript = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&mutated_transcript);
                mutation_pos_on_sequence = &mutated_transcript.len() - mutation_pos_on_sequence -1 ;
            }
            break;
            //println!("{}\n{}", &trx_seq, mutated_transcript);
            //return MutatedTranscript{mutation_pos_on_sequence: mutation_pos_on_sequence, mutated_sequence:mutated_transcript};
        }
    }
    return MutatedTranscript{mutation_pos_on_sequence: mutation_pos_on_sequence, mutated_sequence:mutated_transcript};   
}

pub fn transcript_mutate_mnv(trx_id: &str, mutation: &crate::common::bcf_parsing::MutBedStruct, trx_seq: &String, trx_exons: &Vec < crate::common::gtf_parsing::ExonStruct>) -> MutatedTranscript
{
    let mut mutation_pos_on_sequence: usize = 0;
    let mut mutated_transcript: String = String::new();

    for exon in trx_exons
    {
        if rust_tools::utils::genomic_overlap::is_overlap(exon.start, exon.end, mutation.start as usize, mutation.end as usize)
        {
            log::debug!("{}:{:?} -> {:?}",trx_id,exon, mutation);
            mutation_pos_on_sequence = (mutation.end as usize - exon.start - 1) + exon.seq_start;
            
            mutated_transcript = crate::transcript_mutator::sequence::multi_mutation_mutator(&trx_seq, mutation_pos_on_sequence, mutation.alt_allele.as_str(), mutation.ref_allele.as_str());
            if exon.strand == "-".to_string() 
            {
                mutated_transcript = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&mutated_transcript);
                mutation_pos_on_sequence = &mutated_transcript.len() - mutation_pos_on_sequence -1 ;
            }
            //println!("{}\n{}", &trx_seq, mutated_transcript);
            break;
        }
    }
    return MutatedTranscript{mutation_pos_on_sequence: mutation_pos_on_sequence, mutated_sequence:mutated_transcript};
}

pub fn transcript_mutate_del(trx_id: &str, mutation: &crate::common::bcf_parsing::MutBedStruct, trx_seq: &String, trx_exons: &Vec < crate::common::gtf_parsing::ExonStruct>) -> MutatedTranscript
{
    let mut mutation_pos_on_sequence: usize = 0;
    let mut mutated_transcript: String = String::new();

    for exon in trx_exons
    {
        if rust_tools::utils::genomic_overlap::is_overlap(exon.start, exon.end, mutation.start as usize, mutation.end as usize)
        {
           
            mutation_pos_on_sequence= (mutation.end as usize - exon.start - 1) + exon.seq_start;
            if mutation_pos_on_sequence == exon.seq_end
            {
                log::debug!("{} on trx_id {} exceeds transcript's sequence boundaries", mutation.mutation_id, trx_id);
                return MutatedTranscript{mutation_pos_on_sequence: 0, mutated_sequence:"NA".to_string()};
            }
            log::debug!("{}:{:?} -> {:?} {}",trx_id,exon, mutation, mutation_pos_on_sequence);
            mutated_transcript = crate::transcript_mutator::sequence::deletion_mutator(&trx_seq, mutation_pos_on_sequence, mutation.ref_allele.as_str());

            if exon.strand == "-".to_string()
            {
                mutated_transcript = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&mutated_transcript);
                mutation_pos_on_sequence = &mutated_transcript.len() - mutation_pos_on_sequence -1 ;
            }
            //println!("{}\n{}", &trx_seq, mutated_transcript);
            break;
        }
    }
    return MutatedTranscript{mutation_pos_on_sequence: mutation_pos_on_sequence, mutated_sequence:mutated_transcript};
}


pub fn transcript_mutate_ins(trx_id: &str, mutation: &crate::common::bcf_parsing::MutBedStruct, trx_seq: &String, trx_exons: &Vec < crate::common::gtf_parsing::ExonStruct>) -> MutatedTranscript
{
    let mut mutation_pos_on_sequence: usize = 0;
    let mut mutated_transcript: String = String::new();

    for exon in trx_exons
    {
        if rust_tools::utils::genomic_overlap::is_overlap(exon.start, exon.end, mutation.start as usize, mutation.end as usize)
        {
            log::debug!("{}:{:?} -> {:?}",trx_id,exon, mutation);
            mutation_pos_on_sequence = (mutation.end as usize - exon.start - 1) + exon.seq_start;
            mutated_transcript = crate::transcript_mutator::sequence::insertion_mutator(&trx_seq, mutation_pos_on_sequence, mutation.alt_allele.as_str(), mutation.ref_allele.as_str());

            if exon.strand == "-".to_string()
            {
                mutated_transcript = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&mutated_transcript);
                mutation_pos_on_sequence = &mutated_transcript.len() - mutation_pos_on_sequence -1 ;
            }
            break;
        }
    }
    return MutatedTranscript{mutation_pos_on_sequence: mutation_pos_on_sequence, mutated_sequence:mutated_transcript};
}
  


/// Main function that loops over the retained transcripts and generates their mutated fasta sequence
/// based on the mutations provided in the bed or vcf file (will later be implemented).
///
/// Header string contains a variety of key words that here need a bit more details
///     ** R   => Genomic region of the transcript followed by the strand in parenthesis
///     ** ENS => EnsemblID of gene 
///     ** GN  => Gene name 
///     ** MPG => Position of the mutation on the genome 
///     ** MPR => Position of the mutation on the transcript sequence (similar to position on read)
///     ** MID => Mutation ID
///     ** ALT => Alternative allele 
///     ** REF => Reference allele 
///     ** PV  => Protein variant (if known, then provided as AAposAA), else None
///     ** MTG => mutation type at genomic level (SNV, MNV,INDEL)
///     ** MTP => Mutation type at protein level (Stop_gained, Stop_loss, etc.)
///
pub fn transcript_mutator_main(vcf_file: &str, gtf_file: &str, fasta_file: &str,strand_aware: bool, output_file: &str)
{
    

    //let transcript_ids: FxHashMap < String, Vec < MutBedStruct >> = read_mut_bed_file(bed_file);
    log::info!("  * Reading [{}]", vcf_file);

    //let mut transcript_ids: FxHashMap < String, Vec < crate::common::bcf_parsing::MutBedStruct >> = crate::common::bcf_parsing::read_mut_vcf_file(vcf_file);
    let mut vcf_hash: FxHashMap < String, Vec < crate::common::bcf_parsing::VcfEntry>> = crate::common::bcf_parsing::read_vcf_file(vcf_file);
    
    

    log::info!("  * Reading [{}]", gtf_file);
    
    let mut exon_hash = crate::common::gtf_parsing::read_gtf_file(gtf_file);
    
    log::info!("Filtering for variants overlaping known exons");
    let transcript_ids: FxHashMap < String, Vec <crate::common::bcf_parsing::MutBedStruct>> = filter_vcf_hash(&vcf_hash, &exon_hash);
    //println!("{:?}", transcript_ids);

    log::info!("  * Building transcripts");
    let transcript_seq = build_transcripts(&mut exon_hash, fasta_file, &transcript_ids);

    log::info!("  * Processing mutations");

    
    //println!("{:?}", &transcript_ids["ENST00000518987"]);
    
    // Open file to write 
    let output_buffer: File = File::create(output_file).expect("ERROR: cannot open file to write");

    for (transcript, mutation_vector) in transcript_ids.iter()
    {   
        let trx_seq: String = transcript_seq[transcript.as_str()].clone();
        let trx_exons: Vec < crate::common::gtf_parsing::ExonStruct > = exon_hash[transcript.as_str()].clone();

        for mutation in mutation_vector 
        {

            let mutated_transcript = match mutation.mutation_type.as_str()
            {
                "SNV" => transcript_mutate_snv(transcript.as_str(), mutation, &trx_seq, &trx_exons),
                "MNV" => transcript_mutate_mnv(transcript.as_str(), mutation, &trx_seq, &trx_exons),
                "DEL" => transcript_mutate_del(transcript.as_str(), mutation, &trx_seq, &trx_exons),
                "INS" => transcript_mutate_ins(transcript.as_str(), mutation, &trx_seq, &trx_exons),
                _     => MutatedTranscript{mutation_pos_on_sequence: 0, mutated_sequence:"NA".to_string()},
            };
            
            //println!("{}->{}", transcript, mutation.mutation_id);
            //let mutated_transcript: MutatedTranscript = mutate_transcript(transcript.as_str(), mutation, &transcript_seq, &exon_hash);
            //println!("{}\n{} -> {}",transcript.as_str(), mutated_transcript.mutation_pos_on_sequence, mutated_transcript.mutated_sequence);
            if mutated_transcript.mutated_sequence == "NA".to_string()
            {
                continue;
            }
            let header_string: String = format!(">{}_{}_{}_{}_{}",transcript.to_string(), mutation.mutation_id, mutated_transcript.mutation_pos_on_sequence + 1, mutation.mutation_type, mutation.trans_mutation_type); // We add 1 to mutation_pos_on_sequence so that it is 1 based
            write_fasta_entry(&output_buffer, &header_string, &mutated_transcript.mutated_sequence);
        }
    }

}


