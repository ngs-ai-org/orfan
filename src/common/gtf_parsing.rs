use std::str;
use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;

use rustc_hash::FxHashMap;
use core::panic;
use faimm::IndexedFasta;


extern crate pretty_env_logger;


/// ExonStruct contains information regarding an exon like the chrom, start, end and exon_number
/// These are used for later extracting the sequence corresponding to each exon for building 
/// the transcript sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ExonStruct
{
    pub chrom : String, 
    pub start : usize, 
    pub end   : usize, 
    pub exon_len: usize,
    pub seq_start: usize, 
    pub seq_end: usize,
    pub strand: String,
    pub exon_number : u32,
    pub ensembl_id: String, 
    pub gene_name: String
}

impl ExonStruct 
{
    pub fn new(chrom: String, start: usize, end: usize,seq_start:usize, seq_end:usize, strand:String, exon_number: u32, ensembl_id: String, gene_name: String) -> ExonStruct
    {

        ExonStruct { 
            chrom:chrom, 
            start: start, 
            end: end, 
            exon_len: end - start,
            seq_start: seq_start, 
            seq_end: seq_end, 
            strand: strand, 
            exon_number: exon_number, 
            ensembl_id: ensembl_id, 
            gene_name: gene_name 
        }
    }
}

/// When reading the 9th column of a GTF file you have to 
/// remove spaces, semi-colomns and double-quotes. This function 
/// deals with that to avoid having to do it multiple times
pub fn clean_info_gtf_entry(entry: &str) -> String
{
    let tmp: Vec < &str > = entry.split(".").collect();
    let mut tmp_entry: String = tmp[0].to_string();
    tmp_entry.retain(|c| c != '"');
    tmp_entry.retain(|c| c != ';');
    return tmp_entry;
}


pub fn get_attributes(info_vec: &Vec < &str >) -> FxHashMap<String, String>
{
    let mut attributes_hash: FxHashMap<String,String> = FxHashMap::default();
    
    let mut i: usize = 0;
    while i < info_vec.len()
    {
        attributes_hash.insert(clean_info_gtf_entry(info_vec[i]), clean_info_gtf_entry(info_vec[i+1]));
        i += 2;
    }

    return attributes_hash;
}


/// Function that reads a gtf file exon entries;
pub fn read_gtf_file(filename: &str) -> FxHashMap < String, Vec < ExonStruct >>
{
    log::info!("Reading gtf_file");
    let mut linecount = 0;

    let my_file = Path::new(&filename);
    let input   = File::open(&my_file).expect("Couldn't open file");
    let buffer  = BufReader::new(input);
    
    let mut gtf_hash: FxHashMap < String, Vec < ExonStruct >> = FxHashMap::default();
    for line in buffer.lines()
    {
        let l = line.expect("Couldn't parse line");
        let e: Vec < &str > = l.split("\t").collect();

        if &e[0][..2] == "##"
        {
            continue;
        }
        if e[2] != "exon"
        {
            continue;
        }
        
        if linecount % 100000 == 0 { log::info!("Read {}", linecount); }
        linecount += 1;

        let mut chrom: String = e[0].to_string();
        chrom = chrom.replace("chr","").to_string();
        let start: usize = e[3].parse::<usize>().unwrap();
        let end: usize = e[4].parse::<usize>().unwrap();
        let strand: String = e[6].to_string();
        let exon_len: usize = end - start + 1;

        let info: Vec < &str > =  e[8].split(" ").collect();
        if info.len() % 2 != 0
        {
            panic!("The info column index 8 has unmatched elements!. Please check your GTF/GFF");
        }
        // Get all attributes. 
        let attributes: FxHashMap<String, String> = get_attributes(&info);
        // Get transcript id from attributes
        let transcript: String = attributes["transcript_id"].to_string();
        
        // ensemblID //
        let ensembl_id: String = attributes["gene_id"].to_string();

        // gene_name // 
        let gene_name: String = attributes["gene_name"].to_string();
        
        let mut exon_number: String = attributes["exon_number"].to_string();
        exon_number.retain(|c| c != ';');
        
        if !gtf_hash.contains_key(&transcript){
            let mut transcript_arr: Vec < ExonStruct > = Vec::new();
            let seq_start: usize = 0;
            let seq_end: usize = exon_len -1;
            let exon_struct: ExonStruct = ExonStruct::new(chrom, start-1, end,seq_start, seq_end, strand, exon_number.parse::<u32>().unwrap(), ensembl_id,gene_name);
            
            transcript_arr.push(exon_struct.clone());
            gtf_hash.insert(transcript, transcript_arr);

        }
        else
        {
            let mut transcript_arr: Vec < ExonStruct > = gtf_hash[&transcript].clone();
            let previous_exon: ExonStruct = transcript_arr[transcript_arr.len() -1].clone();
            let seq_start: usize = previous_exon.seq_end + 1;
            let seq_end: usize = seq_start + exon_len -1;
            let exon_struct: ExonStruct = ExonStruct::new(chrom, start-1, end,seq_start, seq_end, strand, exon_number.parse::<u32>().unwrap(), ensembl_id,gene_name);
            transcript_arr.push(exon_struct.clone());
            *gtf_hash.get_mut(&transcript).unwrap() = transcript_arr;
        }
    }
    return gtf_hash;
}



