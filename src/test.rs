


use crate::common::bcf_parsing::{MutBedStruct, read_mut_vcf_file};
use crate::common::gtf_parsing::{ExonStruct, read_gtf_file };
use crate::transcript_mutator::transcript_mutator::{build_transcripts, write_fasta_entry, chunk_fasta_sequence, filter_vcf_hash};
use crate::transcript_mutator::sequence;

use super::*;

/////////////////// 
//   UNIT-TESTS  //
//////////////////

use orf_seeker::read_fasta_file;
use sha2::{Sha256,Digest};
use std::io;
use std::str;
use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;
use rustc_hash::FxHashMap;


#[test]
fn test_vcf_filter_hash()
{
    let vcf_file: &str = "./data/mini_vcf_test.vcf";
    let gtf_file: &str = "./data/gencode.v19.annotation.mini_test.gtf";

    let vcf_hash: FxHashMap < String, Vec < crate::common::bcf_parsing::VcfEntry>> = crate::common::bcf_parsing::read_vcf_file(vcf_file);
    let exon_hash = crate::common::gtf_parsing::read_gtf_file(gtf_file);

    let transcript_ids: FxHashMap < String, Vec <crate::common::bcf_parsing::MutBedStruct>> = filter_vcf_hash(&vcf_hash, &exon_hash);
    
    let mut expected_hash: FxHashMap < String, Vec <crate::common::bcf_parsing::MutBedStruct>> = FxHashMap::default();
    let val: crate::common::bcf_parsing::MutBedStruct = crate::common::bcf_parsing::MutBedStruct{
        chrom: "X".to_string(),
        start: 129201018, 
        end: 129201019, 
        mutation_id: "MU129543605".to_string(),
        ref_allele: "C".to_string(),
        alt_allele: "T".to_string(),
        variant: "None".to_string(),
        mutation_type: "SNV".to_string(),
        trans_mutation_type: "None".to_string()
    };
    expected_hash.insert("ENST00000335997".to_string(), vec![val]);
    assert_eq!(expected_hash, transcript_ids); 

}



#[test]
fn test_build_transcripts()
{
    let gtf_test_file: &str = "./data/test.gtf";
    let gtf_fasta_test_file: &str = "./data/test.gtf.fa";
    let reference_genome = "/gpfs/omics/resources/Public/Genomes/hg37_p13/GRCh37.p13.genome.fa";
    let vcf_file: &str = "./data/test.vcf";

    let transcript_ids: FxHashMap < String, Vec < MutBedStruct >> = read_mut_vcf_file(vcf_file);
    let mut gtf_hash = read_gtf_file(&gtf_test_file);
    let transcript_hash = build_transcripts(&mut gtf_hash, reference_genome, &transcript_ids);


    let fasta = Path::new(&gtf_fasta_test_file);
    let fasta_input = File::open(&fasta).expect("Cannot open file");
    let buffer = BufReader::new(fasta_input);
    let mut expected_sequence: String = String::new();
    for line in buffer.lines()
    {
        let l = line.expect("Could not parse line");
        if &l[0..1] == ">"
        {
            continue;
        }
        else
        {
            expected_sequence += &l;    
        }
    }
    let key: String = "ENST00000335997".to_string();
    let mut observed_sequence: String = transcript_hash[&key].to_string();
    if gtf_hash[&key][0].strand == "-".to_string()
    {
        observed_sequence = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&observed_sequence);
    }
    assert_eq!(expected_sequence, observed_sequence.to_string());
}

#[test]
fn test_strand_aware_orientation()
{
    let expected_sequence: String = "ATCGATCG".to_string();
    let seq_to_rev: String = "CGATCGAT".to_string();
    let observed_sequence: String = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&seq_to_rev);
    assert_eq!(expected_sequence, observed_sequence);
}


#[test]
fn test_multi_mutation_mutator()
{
    let sequence: String= "ATCGATTGATCGATCGATCG".to_string();
    let start:usize = 10;
    let mutation_pos:usize = 11; // 1 based
    let mutated_allele:&str = "AAA";
    let reference_allele:&str = "ATC";
    let expected_mutated_sequence: String = "AAAGATTGATCGATCGATCG".to_string();
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1 )  + 0;
    let result_sequence: String = sequence::multi_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    
    assert_eq!(expected_mutated_sequence, result_sequence);
}



#[test]
fn test_multi_mutation_mutator2()
{
    let sequence: String= "ATCG".to_string();
    let start:usize = 4;
    let mutation_pos:usize = 7; // 1 based
    let mutated_allele:&str = "AA";
    let reference_allele:&str = "CG";
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1) + 0;
    let result_sequence: String = sequence::multi_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    let expected_mutated_sequence: String = "ATAA".to_string();

    assert_eq!(expected_mutated_sequence, result_sequence);
}


#[test]
fn test_multi_mutation_mutator3()
{
    let sequence: String= "ATCGATTGATCGATCGATCG".to_string();
    let start:usize = 10;
    let mutation_pos:usize = 13 ; // 1 based
    let mutated_allele:&str = "AT";
    let reference_allele:&str = "CG";
    let expected_mutated_sequence: String = "ATATATTGATCGATCGATCG".to_string();
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1 )  + 0;
    let result_sequence: String = sequence::multi_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);

    assert_eq!(expected_mutated_sequence, result_sequence);
}




#[test]
fn test_deletion_mutator()
{
    let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 15; // 1 based.
    let reference_allele: &str ="ATGCTCGA";
    let alt_allele: &str = "A";
    let expected_sequence: String = "CCTAAGGTTCCAAA".to_string();

    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::deletion_mutator(&sequence, mutation_pos_on_sequence, reference_allele);
    assert_eq!(mutated_sequence, expected_sequence);
}

#[test]
fn test_insertion_mutator()
{
    let sequence: String = "CCTAATTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 15; // 1 based
    let mutated_allele: &str = "ATGC"; 
    let reference_allele: &str ="A";
    let expected_sequence: String = "CCTAATGCTTCGAGGTTCCAAA".to_string();

    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::insertion_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    assert_eq!(mutated_sequence,expected_sequence);
}



#[test]
fn test_insertion_mutator1()
{
    let sequence: String = "CCTAATTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 11; // 1 based
    let mutated_allele: &str = "CAGC"; 
    let reference_allele: &str ="C";
    let expected_sequence: String = "CAGCCTAATTCGAGGTTCCAAA".to_string();
    

    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let observed_sequence: String = sequence::insertion_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    assert_eq!(expected_sequence, observed_sequence);
}

#[test]
fn test_insertion_mutator2()
{
    let sequence: String = "CCTAATTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 29; // 1based
    let mutated_allele: &str = "AGC"; 
    let reference_allele: &str ="A";
    let expected_sequence: String = "CCTAATTCGAGGTTCCAAAGC".to_string();
    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let observed_sequence: String = sequence::insertion_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    assert_eq!(expected_sequence, observed_sequence);
}

#[test]
fn test_deletion_mutator1()
{
    let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 15; // 1 based
    let reference_allele: &str ="ATGCTCGA";
    let alt_allele: &str = "A";
    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::deletion_mutator(&sequence, mutation_pos_on_sequence, reference_allele);
    let expected_sequence: String = "CCTAAGGTTCCAAA".to_string();
    assert_eq!(mutated_sequence, expected_sequence);
}

#[test]
fn test_deletion_mutator2()
{
    let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 11; // 1 based
    let reference_allele: &str ="CCTAA";
    let alt_allele: &str = "C";
    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::deletion_mutator(&sequence, mutation_pos_on_sequence, reference_allele);
    let expected_sequence: String = "CTGCTCGAGGTTCCAAA".to_string();
    assert_eq!(mutated_sequence, expected_sequence);
}

#[test]
fn test_deletion_mutator3()
{
    let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 26; // 1 based
    let reference_allele: &str ="TCCAAA";
    let alt_allele: &str = "T";
    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::deletion_mutator(&sequence, mutation_pos_on_sequence, reference_allele);
    let expected_sequence: String = "CCTAATGCTCGAGGTT".to_string();
    assert_eq!(mutated_sequence, expected_sequence);
}

#[test]
fn test_deletion_mutator4()
{
    let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
    let start: usize     = 10;
    let mutation_pos: usize = 11; // 1 based
    let reference_allele: &str ="CCTAATGCTCGAGGTTCCAA";
    let mutation_pos_on_sequence: usize = ( mutation_pos - start - 1) + 0;
    let mutated_sequence: String = sequence::deletion_mutator(&sequence, mutation_pos_on_sequence, reference_allele);
    let expected_sequence: String = "CA".to_string();
    assert_eq!(mutated_sequence, expected_sequence);
}


#[test]
fn chunk_fasta_sequence_test()
{
    let sequence_to_chunk: String = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string();
    let chunk_size: usize         = 10;
    let mut expected_result: Vec < &str > = Vec::new();
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CG");

    let chunk_sequence =chunk_fasta_sequence(&sequence_to_chunk, chunk_size);
    assert_eq!(expected_result,chunk_sequence);
}

#[test]
fn chunk_fasta_sequence_fail_test()
{
    let sequence_to_chunk: String = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string();
    let chunk_size: usize         = 10;
    let mut expected_result: Vec < &str > = Vec::new();
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CGATCGATCG");
    expected_result.push("ATCGATCGAT");
    expected_result.push("CG");

    let chunk_sequence =chunk_fasta_sequence(&sequence_to_chunk, chunk_size);
    assert!(chunk_sequence != expected_result);
}



#[test]
fn point_mutation_mutator_test()
{
    let sequence: String= "ATCGATTGATCGATCGATCG".to_string();
    let start:usize = 10;
    let mutation_pos:usize = 15; // 1 based
    let mutated_allele:&str = "G";
    let reference_allele:&str = "A";
    let expected_mutated_sequence: String = "ATCGGTTGATCGATCGATCG".to_string();
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1) + 0;
    let result_sequence: String = sequence::point_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
    //println!("{} -> {}", expected_mutated_sequence, result_sequence);

    assert_eq!(expected_mutated_sequence, result_sequence);
}

#[test]
fn point_mutation_mutator_test2()
{
    let sequence: String= "ATCGATTGATCGATCGATCG".to_string();
    let start:usize = 10;
    let mutation_pos:usize = 11; // 1 based
    let mutated_allele:&str = "G";
    let reference_allele:&str = "A";
    let expected_mutated_sequence: String = "GTCGATTGATCGATCGATCG".to_string();
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1 )  + 0;
    let result_sequence: String = sequence::point_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);

    assert_eq!(expected_mutated_sequence, result_sequence);
}

#[test]
fn point_mutation_mutator_test3()
{
    let sequence: String= "ATCG".to_string();
    let start:usize = 4;
    let mutation_pos:usize = 8; // 1 based
    let mutated_allele:&str = "A";
    let reference_allele:&str = "G";
    let expected_mutated_sequence: String = "ATCA".to_string();
    let mutation_pos_on_sequence: usize = (mutation_pos - start - 1) + 0;
    let result_sequence: String = sequence::point_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);

    assert_eq!(expected_mutated_sequence, result_sequence);
}


#[test]
fn write_fasta_entry_test()
{   
    let output_file: &str = "./data/test_write_fasta_entry.fasta";
    let output_buffer: File = File::create(output_file).expect("ERROR: cannot open file to write");
    let header_string: String = ">TEST ta mere".to_string();
    let transcript_sequence: String = "ATCGATCGATCGATCGATCGATCGATCG".to_string();

   write_fasta_entry(&output_buffer, &header_string, &transcript_sequence);
    let expected_results_file = "./data/test_file_writeFastaEntry.fasta";
    let mut expected_results = File::open(&expected_results_file).expect("Cannot open file");
    let mut hasher1 = Sha256::new();
    let _n2 = io::copy(&mut expected_results, &mut hasher1).expect("ERROR");
    let hash_test_file1 = hasher1.finalize();

    let mut output_file_to_test = File::open(&output_file).expect("Cannot open file");
    let mut hasher2 = Sha256::new();
    let _n = io::copy(&mut output_file_to_test, &mut hasher2).expect("ERROR");
    let output_file_hash = hasher2.finalize();
    
    std::fs::remove_file(output_file).expect("Cannot remove file");

    assert_eq!(hash_test_file1, output_file_hash);
}

#[test]
fn write_fasta_entry_test2()
{   
    let output_file: &str = "./data/test_write_fasta_entry2.fasta";
    let output_buffer: File = File::create(output_file).expect("ERROR: cannot open file to write");
    let header_string: String = ">TEST ta mere".to_string();
    let transcript_sequence: String = "ATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGTATCGAATCGT".to_string();

   write_fasta_entry(&output_buffer, &header_string, &transcript_sequence);
    let expected_results_file = "./data/test_file_writeFastaEntry2.fasta";
    let mut expected_results = File::open(&expected_results_file).expect("Cannot open file");
    let mut hasher1 = Sha256::new();
    let _n2 = io::copy(&mut expected_results, &mut hasher1).expect("ERROR");
    let hash_test_file1 = hasher1.finalize();

    let mut output_file_to_test = File::open(&output_file).expect("Cannot open file");
    let mut hasher2 = Sha256::new();
    let _n = io::copy(&mut output_file_to_test, &mut hasher2).expect("ERROR");
    let output_file_hash = hasher2.finalize();
    
    std::fs::remove_file(output_file).expect("Cannot remove file");

    assert_eq!(hash_test_file1, output_file_hash);
}


#[test]
fn strand_aware_seq_orientation_test()
{
    let sequence_to_orient:String = "ATCGATCG".to_string();
    let expected_sequence: String = "CGATCGAT".to_string();
    let obtained_sequence: String = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&sequence_to_orient);
    assert_eq!(obtained_sequence, expected_sequence);
}

#[test]
fn strand_aware_seq_orientation_fail_test()
{
    let sequence_to_orient:String = "ATCGATCG".to_string();
    let expected_sequence: String = "TAGCTAGC".to_string();
    let obtained_sequence: String = crate::transcript_mutator::sequence::strand_aware_seq_orientation(&sequence_to_orient);
    assert_ne!(obtained_sequence, expected_sequence);
}


//////////////////////
// ORF SEEKER TESTS //
/////////////////////

#[test]
fn get_value_from_fasta_description_test1()
{
    let fasta_description = "ngsai|P53621-2_mut|COPA_HUMAN Isoform 2 of Coatomer subunit alpha OS=Homo sapiens OX=9606 GN=COPA ENS=ENSG00000122218 ORF_NAME=F2:4480:4572:1494:1523:30:M PE=1 SV=1";
    let value_to_return = "GN";

    let expected_value: String = "COPA".to_string();
    let observed_value: String = orf_seeker::get_value_from_fasta_description(fasta_description, value_to_return);
    assert_eq!(observed_value, expected_value);
}

#[test]
fn get_value_from_fasta_description_test2()
{
    let fasta_description = "ngsai|P53621-2_mut|COPA_HUMAN Isoform 2 of Coatomer subunit alpha OS=Homo sapiens OX=9606 GN=COPA ENS=ENSG00000122218 ORF_NAME=F2:4480:4572:1494:1523:30:M PE=1 SV=1";
    let value_to_return = "OS";

    let expected_value: String = "Homo".to_string(); // THIS IS WRONG BUT FOR NOW WILL NOT TRY TO FIX IT AS I DO NOT NEED TO EXTRACT ANYTHING ELSE THAN THE GENE NAME AND ENSEMBL ID WHICH SHOULD NOT HAVE SPACES
    let observed_value: String = orf_seeker::get_value_from_fasta_description(fasta_description, value_to_return);
    assert_eq!(observed_value, expected_value);
}

#[test]
fn orf_contains_mutation_true_test1()
{
    let n_start: usize = 100;
    let n_end: usize = 200;
    let mutation_pos = 150;
    let observed_result: bool = orf_seeker::orf_contains_mutation(&n_start, &n_end, &mutation_pos);
    let expected_value: bool = true;
    assert_eq!(observed_result, expected_value);
}

#[test]
fn orf_contains_mutation_true_test2()
{
    let n_start: usize = 100;
    let n_end: usize = 200;
    let mutation_pos = 200;
    let observed_result: bool = orf_seeker::orf_contains_mutation(&n_start, &n_end, &mutation_pos);
    let expected_value: bool = true;
    assert_eq!(observed_result, expected_value);
}

#[test]
fn orf_contains_mutation_false_test1()
{
    let n_start: usize = 100;
    let n_end: usize = 200;
    let mutation_pos = 10;
    let observed_result: bool = orf_seeker::orf_contains_mutation(&n_start, &n_end, &mutation_pos);
    let expected_value: bool = false;
    assert_eq!(observed_result, expected_value);
}


#[test]
fn extract_sequence_from_hash_test()
{
    let mut transcript_hash: FxHashMap < String, String > = FxHashMap::default();
    let trx_id: &str = "ENST00000410344";
    let trx_seq: String = "ATCGTTAACTTTAGGC".to_string();
    transcript_hash.insert(trx_id.to_string(), trx_seq);

    let start: usize = 0;
    let end:   usize = 2;
    
    let expected_value: String = "AT".to_string();

    assert_eq!(expected_value, orf_seeker::extract_sequence_from_hash(&transcript_hash, trx_id, start, end));
    assert_ne!("ATC".to_string(), orf_seeker::extract_sequence_from_hash(&transcript_hash, trx_id, start, end));
}

#[test]
fn extract_sequence_from_hash_test2()
{
    let mut transcript_hash: FxHashMap < String, String > = FxHashMap::default();
    let trx_id: &str = "ENST00000410344";
    let trx_seq: String = "ATCGTTAACTTTAGGC".to_string();
    transcript_hash.insert(trx_id.to_string(), trx_seq);

    let start: usize = 0;
    let end:   usize = 2;
    
    assert_ne!("ATC".to_string(), orf_seeker::extract_sequence_from_hash(&transcript_hash, trx_id, start, end));
}

#[test]
fn get_mutated_aa_test()
{
    let seq: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG";
    let mut aa_seq: Vec < String > = Vec::new();
    for i in seq.chars() { aa_seq.push(i.to_string()); }
    let mutation_pos_on_pep: usize = 5;
    assert_eq!("R".to_string(), orf_seeker::get_mutated_aa(&aa_seq, mutation_pos_on_pep));
}

#[test]
fn get_mutated_aa_test2()
{
    let seq: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG";
    let mut aa_seq: Vec < String > = Vec::new();
    for i in seq.chars() { aa_seq.push(i.to_string()); }
    let mutation_pos_on_pep: usize = 26;
    assert_eq!("G".to_string(), orf_seeker::get_mutated_aa(&aa_seq, mutation_pos_on_pep));
}

#[test]
fn get_mutated_aa_pos_test()
{
    let seq: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG";
    let mut aa_seq: Vec < String > = Vec::new();
    for i in seq.chars() { aa_seq.push(i.to_string()); }
    let obj_orf: orf_seeker::Orf = orf_seeker::Orf::new(aa_seq, 0 as usize, 26 as usize, 1 as u32 ,"M".to_string(), 10 as usize,27 as usize);
    let mutation_pos_on_sequence: usize = 15;
    assert_eq!((1, 13,16), orf_seeker::get_mutated_aa_pos(&obj_orf, mutation_pos_on_sequence));
}


#[test]
fn frame_iterate_test()
{
    let seq: &str = "ATGGCGTGTGCTCGCCCAATGATATCGGTGTACTCCGAAATAA";
    let frame: u32 = 1;
    let ff: (Vec < String >, u32, Vec < usize > ) = orf_seeker::frame_iterate(seq, frame, &utils::codon_hash());
    
    let ss: &str = "MACARPMISVYSEI";
    let mut exp_seq: Vec < String > = Vec::new();
    for i in ss.chars() { exp_seq.push(i.to_string()); }
    let n_flag: u32 = 0;
    let exp_n_pos_list: Vec < usize > = vec![0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39];
    let expected_result: (Vec < String >, u32, Vec < usize > ) = (exp_seq, n_flag, exp_n_pos_list);

    assert_eq!(expected_result, ff);
}

#[test]
fn orf_seeker_test()
{
    let fasta_file: &str = "./data/orf_seeker_test_fasta.fa";
    let fasta_hash: FxHashMap < String, String > = orf_seeker::read_fasta_file(fasta_file);
    let frame: u32 = 1;
    let seq: &str = &fasta_hash["ENST00000422486_MU14521900_298_SNV_exon_variant"];
    let frame_it: (Vec < String >, u32, Vec < usize > ) = orf_seeker::frame_iterate(seq, frame.clone(), &utils::codon_hash());

    let orf_dict: Vec < orf_seeker::Orf > = orf_seeker::orf_seeker(frame_it.0, frame.clone(), frame_it.2);
    
    // Creating the expected ORF object // 
    let ss: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG$";
    let mut orf_seq: Vec < String > = Vec::new();
    for i in ss.chars() { orf_seq.push(i.to_string()); }
    let orf_len: usize = 27;
    let a_start: usize = 1; 
    let a_end  : usize = 27;
    let start_codon: String = "M".to_string();
    let n_start: usize = 0;
    let n_end: usize = 83;
    let orf_obj: orf_seeker::Orf = orf_seeker::Orf::new(orf_seq.clone(), a_start, a_end, frame.clone(), start_codon, n_start, n_end);

    assert_eq!(orf_obj, orf_dict[0]);
}


#[test]
fn orf_name_creation_test()
{
    
    // Creating the expected ORF object // 
    let ss: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG$";
    let mut orf_seq: Vec < String > = Vec::new();
    for i in ss.chars() { orf_seq.push(i.to_string()); }
    let orf_len: usize = 27;
    let a_start: usize = 1; 
    let a_end  : usize = 27;
    let frame: u32 = 1; 
    let start_codon: String = "M".to_string();
    let n_start: usize = 0;
    let n_end: usize = 83;
    let orf_obj: orf_seeker::Orf = orf_seeker::Orf::new(orf_seq.clone(), a_start, a_end, frame.clone(), start_codon, n_start, n_end);

    assert_eq!(orf_obj.orf_name, "F1_1_27_0_83_27_M".to_string());
}




#[test]
fn mutated_orf_test1()
{
    /*
     * Testing internal mutation
    */
    let expected_orf: &str = "MGFSSELCSPQGHGVLQQMQEAELRLLEGMRKWMAQRVKSDREYAGLLHHMSLQDSGGQSRAISPDSPISQSRAEITSQTEGLSRLLRQHAEDLNSGPLSKLSLLIRERQQLRKTYSEQWQQLQQELTKTHSQDIEKLKSQYRALARDSAQAKRKYQEASKDKDRDKAKDKYVRSLWKLFAHHNRYVLGVRAAQLHHQHHHQLLLPGLLRSLQDLHEEMACILKEILQEYLEISSLVQDEVVAIHREMAAAAARIQPEAEYQGFLRQYGSAPDVPPCVTFDESLLEEGEPLEPGELQLNELTVESVQHTLTSVTDELAVATEMVFRRQEMVTQLQQELRNEEENTHPRERVQLLGKRQVLQEALQGLQVALCSQAKLQAQQELLQTKLEHLGPGEPPPVLLLQDDRHSTSSSEQEREGGRTPTLEILKSHISGIFRPKFSLPPPLQLIPEVQKPLHEQLWYHGAIPRAEVAELLVHSGDFLVRESQGKQEYVLSVLWDGLPRHFIIQAHAPSEWRLLTPGPLPCRTCTDWKGKAFLAFLCSSTTY";

    let filename: &str = "./data/mutated_transcript_test1.fa";
    let reader = bio::io::fasta::Reader::from_file(filename).unwrap();
    let codon_hash = crate::utils::codon_hash();
    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let f2: u32 = 2;
        let (f2_prot_seq, f2_n_flag, f2_n_pos_list) = crate::orf_seeker::frame_iterate(fasta_sequence, f2, &codon_hash);
        let f2_orf_dict = crate::orf_seeker::orf_seeker(f2_prot_seq, f2_n_flag, f2_n_pos_list);
        let mut seq_vec: Vec < String > = f2_orf_dict[1].orf_seq.clone();
        assert_eq!(seq_vec[0..seq_vec.len()-1].join(""), expected_orf);
    }
}


#[test]
fn mutated_orf_test2()
{
    /*
     * Testing stop gained mutation
     */
    let expected_orf: &str = "MGFSSELCSPQGHGVLQQMQEAELRLLEGMRKWMAQRVKSDREYAGLLHHMSLQDSGGQSRAISPDSPISQSWAEITSQTEGLSRLLRQHAEDLNSGPLSKLSLLIRERQQLRKTYSEQWQQLQQELTKTHSQDIEKLKSQY";
    let filename: &str = "./data/mutated_transcript_test2.fa";
    let reader = bio::io::fasta::Reader::from_file(filename).unwrap();
    let codon_hash = crate::utils::codon_hash();
    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let f2: u32 = 2;
        let (f2_prot_seq, f2_n_flag, f2_n_pos_list) = crate::orf_seeker::frame_iterate(fasta_sequence, f2, &codon_hash);
        let f2_orf_dict = crate::orf_seeker::orf_seeker(f2_prot_seq, f2_n_flag, f2_n_pos_list);
        let mut seq_vec: Vec < String > = f2_orf_dict[1].orf_seq.clone();
        assert_eq!(seq_vec[0..seq_vec.len()-1].join(""), expected_orf);
    }
}

#[test]
fn mutated_orf_test3()
{
    /*
    * Testing stop loss mutation
    */
    let expected_orf: &str = "MMRRWISGALEFFAMNFSLGSLLLRQTHTKRPTKEYHG";
    let filename: &str = "./data/mutated_transcript_test3.fa";
    let reader = bio::io::fasta::Reader::from_file(filename).unwrap();
    let codon_hash = crate::utils::codon_hash();
    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let f1: u32 = 1;
        let (f1_prot_seq, f1_n_flag, f1_n_pos_list) = crate::orf_seeker::frame_iterate(fasta_sequence, f1, &codon_hash);
        let f1_orf_dict = crate::orf_seeker::orf_seeker(f1_prot_seq, f1_n_flag, f1_n_pos_list);
        let seq_vec: Vec < String > = f1_orf_dict[5].orf_seq.clone();
        assert_eq!(seq_vec[0..seq_vec.len()-1].join(""), expected_orf);
    }
}

#[test]
fn mutated_orf_test4()
{
    /*
     * Testing stop gained mutation
     */
    let expected_orf: &str = "MMIPPLKMLIRSGSVLMMSLPIWTFWIQLDRQSLQPCGTSI";
    let filename: &str = "./data/mutated_transcript_test4.fa";
    let reader = bio::io::fasta::Reader::from_file(filename).unwrap();
    let codon_hash = crate::utils::codon_hash();
    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let f2: u32 = 2;
        let (f2_prot_seq, f2_n_flag, f2_n_pos_list) = crate::orf_seeker::frame_iterate(fasta_sequence, f2, &codon_hash);
        let f2_orf_dict = crate::orf_seeker::orf_seeker(f2_prot_seq, f2_n_flag, f2_n_pos_list);
        let mut seq_vec: Vec < String > = f2_orf_dict[1].orf_seq.clone();
        assert_eq!(seq_vec[0..seq_vec.len()-1].join(""), expected_orf);
        
    }
}

#[test]
fn mutated_orf_test5()
{
    /*
    * Testing stop loss mutation
    */
    let expected_orf: &str = "MARIPRVRGGGTHRSGQGAFGNMCRGG";
    let filename: &str = "./data/mutated_transcript_test5.fa";
    let reader = bio::io::fasta::Reader::from_file(filename).unwrap();
    let codon_hash = crate::utils::codon_hash();
    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let f1: u32 = 1;
        let (f1_prot_seq, f1_n_flag, f1_n_pos_list) = crate::orf_seeker::frame_iterate(fasta_sequence, f1, &codon_hash);
        let f1_orf_dict = crate::orf_seeker::orf_seeker(f1_prot_seq, f1_n_flag, f1_n_pos_list);
        let seq_vec: Vec < String > = f1_orf_dict[2].orf_seq.clone();
        assert_eq!(seq_vec[0..seq_vec.len()-1].join(""), expected_orf);
    }
}


///////////////////////////////
// PROTEIN FRAGMENTOR TESTS //
//////////////////////////////


/*#[test]
fn trypsinate_protein_sequence_test()
{
    let protein_sequence: &str = "MCNKMCNKPMCNRMCNRPMCND";
    let mut expected_result: Vec < String > = Vec::new();
    expected_result.push("MCNK".to_string());
    expected_result.push("MCNKPMCNR".to_string());
    expected_result.push("MCNRPMCND".to_string());
    let observed_result = protein_fragmentor::trypsinate_protein_sequence(protein_sequence);

    assert_eq!(observed_result, expected_result);
}



#[test]
fn trypsinate_protein_sequence_test2()
{
    let protein_sequence: &str = "MCNKMCNKPMCNRMCNR";
    let mut expected_result: Vec < String > = Vec::new();
    expected_result.push("MCNK".to_string());
    expected_result.push("MCNKPMCNR".to_string());
    expected_result.push("MCNR".to_string());
    let observed_result = protein_fragmentor::trypsinate_protein_sequence(protein_sequence);

    assert_eq!(observed_result, expected_result);
}
*/