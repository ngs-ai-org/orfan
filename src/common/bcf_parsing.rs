use std::str;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead,BufReader};
use std::path::Path;

use rustc_hash::FxHashMap;
use core::panic;
use faimm::IndexedFasta;
use rust_htslib::bcf::{Reader, Read, Record};

extern crate pretty_env_logger;


/////////////
// Struct //
////////////

/// MutBedStruct contains information regarding the mutations obtained from ICGC in bed-like format formated using 
/// /gpfs/omics/development/NGS-AI/nlykosk/proteomics/scripts/mutProtAnnotationCreation/fetchICGCData.py
/// script.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct MutBedStruct
{
    pub chrom : String, 
    pub start : i64,
    pub end   : i64,
    pub mutation_id  : String,
    pub ref_allele : String, 
    pub alt_allele : String,
    pub variant : String,
    pub mutation_type : String, 
    pub trans_mutation_type: String
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct VcfEntry 
{
    pub chrom: String, 
    pub start: i64, 
    pub ref_allele: String, 
    pub alt_allele: String, 
    pub var_id: String, 
}

impl VcfEntry
{
    /// Function that creates a new VcfEntry
    pub fn new(chrom: String, start: i64, ref_allele: String, alt_allele: String, var_id: String) -> VcfEntry
    {
        VcfEntry
        {
            chrom      : chrom, 
            start      : start, 
            ref_allele : ref_allele, 
            alt_allele : alt_allele,
            var_id     : var_id
        }
    }

    /// Function that assesses whether VcfEntry is an SNV
    /// or not. Returns true if VcfEntry contains a SNV or 
    /// false if not.
    pub fn is_snv(&self) -> bool
    {
        if self.ref_allele.len() == 1 && self.alt_allele.len() == 1
        {
            return true;
        }
        return false;
    }
    
    /// Function that returns true if VcfEntry contains an INDEL
    /// or false if not
    pub fn is_indel(&self) -> bool
    {
        if self.ref_allele.len() != self.alt_allele.len()
        {
            return true;
        }
        return false;
    }
    
    /// Function that returns true if VcfEntry is an insertion
    /// or false if not
    pub fn is_ins(&self) -> bool
    {
        if (self.ref_allele.len() != 1 || self.alt_allele.len() != 1) && self.alt_allele.len() > self.ref_allele.len()
        {
            return true;
        }
        return false;
    }
    
    /// Function that returns true if it is a deletion
    /// else returns false
    pub fn is_del(&self) -> bool 
    {
        if (self.ref_allele.len() != 1 || self.alt_allele.len() != 1) && self.ref_allele.len() > self.alt_allele.len()
        {
            return true;
        }
        return false;
    }

    pub fn is_mnv(&self) -> bool
    {
        if (self.ref_allele.len() == self.alt_allele.len()) && self.ref_allele.len() != 1
        {
            return true;
        }
        return false;
    }
}


pub fn read_vcf_file(filename: &str) -> FxHashMap < String, Vec <VcfEntry>>
{
    let mut vcf_entry_vector: Vec < VcfEntry > = Vec::new();
    let mut vcf_hash: FxHashMap<String, Vec < VcfEntry >> = FxHashMap::default();

    let mut linecount:usize = 0;

    let mut bcf = Reader::from_path(filename).expect("Cannot open vcf file!");

    let mut tmp_contig: String = String::new();

    for records in bcf.records()
    {
        if linecount % 10000 == 0 { log::info!("Read {}", linecount); }
        linecount += 1;

        let mut record = records.expect("Fail to read record");
        let mut id: String = String::from_utf8(record.id()).unwrap();
        let t = record.header().rid2name(record.rid().unwrap()).unwrap();
        let contig = String::from_utf8(t.to_vec()).unwrap();

        let pos: i64 = record.pos() + 1;
        let alleles = record.alleles();
        let ref_allele: &str = str::from_utf8(alleles[0]).expect("Could not retrieve reference allele from alleles");
        let alt_allele: &str = str::from_utf8(alleles[1]).expect("Could not retrieve alternative allele from alleles");

        if id == "."
        {
            id = contig.to_string() + "_" + &pos.to_string() + "_" + ref_allele + "_" + alt_allele;
        }
        let entry: VcfEntry = VcfEntry::new(contig.clone(), pos, ref_allele.to_string(), alt_allele.to_string(), id.clone());
        

        if tmp_contig == ""
        {
            tmp_contig = contig;
            vcf_entry_vector.push(entry);
        }
        else if tmp_contig != contig
        {
            vcf_hash.insert(tmp_contig, vcf_entry_vector);
            tmp_contig = contig.clone(); 
            vcf_entry_vector = Vec::new();
            let entry: VcfEntry = VcfEntry::new(contig.clone(), pos, ref_allele.to_string(), alt_allele.to_string(), id.clone());
            vcf_entry_vector.push(entry);
        }
        else 
        {
            vcf_entry_vector.push(entry);
            continue;
        }
    }
    vcf_hash.insert(tmp_contig, vcf_entry_vector);


    return vcf_hash;
}


/// Function that reads a vcf file containing mutations and creates a 
/// FxHashMap
pub fn read_mut_vcf_file(filename: &str) -> FxHashMap < String, Vec < MutBedStruct >>
{ 

    let mut linecount = 0;

    let mut tmp_map: FxHashMap<String, Vec < MutBedStruct >> = FxHashMap::default();

    let mut bcf = Reader::from_path(filename).expect("Cannot open vcf file");
    for (i, record_result) in bcf.records().enumerate()
    {
        if linecount % 10000 == 0 { println!("Read {}", linecount); }
        linecount += 1;

        let mut record = record_result.expect("Fail to read record");
        
        let id = String::from_utf8(record.id()).unwrap();
        
        let tmp = record.header().rid2name(record.rid().unwrap()).unwrap();
        let contig = String::from_utf8(tmp.to_vec()).unwrap();
        println!("{}", contig);
        let pos = record.pos() + 1; // The + 1 is required here because .pos() return 0-based coordinates.

        let alleles = record.alleles();
        let ref_allele = str::from_utf8(alleles[0]).expect("Could not convert allele");
        let alt_allele = str::from_utf8(alleles[1]).expect("Could not convert allele");

        let gmtype1 = record.info(b"GMTYPE").string().unwrap().unwrap();
        let gmtype = str::from_utf8(gmtype1[0]).expect("ERROR");
        
        let aamtype1 = record.info(b"AAMTYPE").string().unwrap().unwrap();
        let aamtype = str::from_utf8(aamtype1[0]).expect("ERROR");

        let tmtype1 = record.info(b"TMTYPE").string().unwrap().unwrap();
        let tmtype = str::from_utf8(tmtype1[0]).expect("ERROR");

        let funimpact1 = record.info(b"FUNIMPACT").string().unwrap().unwrap();
        let funimpact = str::from_utf8(funimpact1[0]).expect("ERROR");
        
        let trx_id = record.info(b"TRXID").string().unwrap().unwrap();
        let transcript_id: String = str::from_utf8(trx_id[0]).expect("ERROR").to_string();

        
        let record_struct: MutBedStruct = MutBedStruct{
            chrom: contig,
            start: pos -1,
            end: pos, 
            mutation_id: id,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant: aamtype.to_string(), 
            mutation_type: gmtype.to_string(),
            trans_mutation_type: tmtype.to_string()
        };
        if !tmp_map.contains_key(&transcript_id)
        {
            let mut tmp_vec: Vec < MutBedStruct > = Vec::new();
            tmp_vec.push(record_struct.clone());
            tmp_map.insert(transcript_id.clone(), tmp_vec);
        }
        else
        {   
            let mut tmp_vec: Vec < MutBedStruct > = tmp_map[&transcript_id].clone();
            tmp_vec.push(record_struct.clone());
            *tmp_map.get_mut(&transcript_id).unwrap() = tmp_vec; 
        }

    }

    return tmp_map;
}


#[test]
fn test_vcf_entry_struct()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "C";
    let alt_allele: &str = "T";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom.clone(), pos.clone(), ref_allele.to_string(), alt_allele.to_string(), id.clone());
    let expected: VcfEntry = VcfEntry{
        chrom: chrom.clone(), 
        start: pos, 
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
        var_id: id.clone()
    };
    assert_eq!(entry, expected);
}

#[test]
fn test_vcf_entry_struct_snv()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "C";
    let alt_allele: &str = "T";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom, pos, ref_allele.to_string(), alt_allele.to_string(), id);
    assert_eq!(entry.is_snv(), true);
    
}

#[test]
fn test_vcf_entry_struct_mnv()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "CCC";
    let alt_allele: &str = "TTT";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom, pos, ref_allele.to_string(), alt_allele.to_string(), id);
    assert_eq!(entry.is_mnv(), true);
}


#[test]
fn test_vcf_entry_struct_del()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "CCT";
    let alt_allele: &str = "C";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom, pos, ref_allele.to_string(), alt_allele.to_string(), id);
    assert_eq!(entry.is_del(), true);
}


#[test]
fn test_vcf_entry_struct_ins()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "C";
    let alt_allele: &str = "CTTC";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom, pos, ref_allele.to_string(), alt_allele.to_string(), id);
    assert_eq!(entry.is_ins(), true);
}


#[test]
fn test_vcf_entry_struct_indel()
{
    let chrom: String = "chr1".to_string();
    let pos: i64 = 1000;
    let ref_allele: &str = "C";
    let alt_allele: &str = "CTTC";
    let id: String = "ta_mere".to_string();

    let entry: VcfEntry = VcfEntry::new(chrom, pos, ref_allele.to_string(), alt_allele.to_string(), id);
    assert_eq!(entry.is_indel(), true);
}


#[test]
fn test_read_vcf_file()
{
    let test_vcf_file: &str = "./data/mini_vcf_test.vcf";
    // Create expected entries and populate expected hash 
    let entry_one: VcfEntry = VcfEntry::new("1".to_string(), 129200723 as i64, "G".to_string(), "A".to_string(), "MU63781225".to_string());
    let entry_two: VcfEntry = VcfEntry::new("2".to_string(), 129200780 as i64, "G".to_string(), "T".to_string(), "MU85781889".to_string());
    let entry_three: VcfEntry = VcfEntry::new("3".to_string(), 129200915 as i64, "C".to_string(), "T".to_string(), "MU2245115".to_string());
    let entry_four: VcfEntry = VcfEntry::new("4".to_string(), 129200959 as i64, "G".to_string(), "A".to_string(), "MU4625800".to_string());
    let entry_x: VcfEntry = VcfEntry::new("X".to_string(), 129201019 as i64, "C".to_string(), "T".to_string(), "MU129543605".to_string());
    let mut expected_hash: FxHashMap<String, Vec < VcfEntry >> = FxHashMap::default();
    expected_hash.insert("1".to_string(), vec![entry_one]);
    expected_hash.insert("2".to_string(), vec![entry_two]);
    expected_hash.insert("3".to_string(), vec![entry_three]);
    expected_hash.insert("4".to_string(), vec![entry_four]);
    expected_hash.insert("X".to_string(), vec![entry_x]);

    let observed_hash: FxHashMap<String, Vec < VcfEntry >> = read_vcf_file(test_vcf_file);
    assert_eq!(expected_hash, observed_hash);
    
}