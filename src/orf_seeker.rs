use std::str;
use std::fs::File;
use std::io::Write;
use rustc_hash::{FxHashMap, FxHasher};
use std::collections::HashSet;
use bio::io::{self, fasta};
use std::io::{BufRead,BufReader};
use std::path::Path;
use crate::utils::codon_hash;


/////////////
// Struct //
////////////

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Orf {
    pub orf_seq: Vec <String>,
    pub orf_len: usize,
    pub a_start: usize,
    pub a_end: usize,
    pub frame: u32,
    pub start_codon: String, 
    //pub orf_length: usize
    pub orf_name: String,
    pub n_start: usize, 
    pub n_end: usize, 

    //pub orf_id: String
}

impl Orf {

    pub fn new(orf_seq: Vec < String >, a_start:usize, a_end:usize, frame:u32, start_codon:String, n_start:usize, n_end:usize) -> Orf 
    {
        let orf_length: usize = orf_seq.len() - 1 as usize; // -1 as the last codon is the stop codon;
        let orf_name: String = format!("F{}_{}_{}_{}_{}_{}_{}", frame, a_start, a_end, n_start, n_end, orf_length, start_codon);
        Orf {
            orf_seq: orf_seq, 
            orf_len: orf_length,
            a_start: a_start,
            a_end: a_end,
            frame: frame, 
            orf_name: orf_name,
            start_codon: start_codon, 
            n_start: n_start, 
            n_end: n_end,
        }
    }
}

pub struct ProtGeneStruct 
{
    pub protein_id: String, 
    pub protein_name: String,
}


//////////////////////
// Public functions //
/////////////////////


/// Function that reads a sequence and returns
/// a vector containing the translated seq from 
/// nucleotides to protein, the NA flag and the 
/// nucleotide position list. 
pub fn frame_iterate<'a>(seq:&str,frame:u32, codon_hash: &FxHashMap < &str, &str>) -> (Vec <String>, u32, Vec <usize>)
{
    let seq_length: usize = seq.len();
    let mut pos = (frame - 1) as usize;
    let mut prot_seq: Vec < String > = Vec::new();
    let mut n_flag: u32= 0;
    let mut n_pos_list: Vec < usize > = Vec::new();

    while pos < seq_length - 2 
    {
        n_pos_list.push(pos);
        let codon = &seq[pos..pos+3];
        if codon_hash.contains_key(codon)
        {
            prot_seq.push(codon_hash[codon].to_string());
            //println!("{} => {}",codon,codon_hash[codon]);
        }
        else
        {
            //println!("sequence has N's");
            n_flag = 1;
            break;
        }
        pos += 3
    }
    return (prot_seq, n_flag, n_pos_list);
}


/// Refactored function from tama_orf_seeker.py
pub fn orf_seeker(prot_seq: Vec < String >, frame: u32, n_pos_list: Vec < usize >) -> Vec < Orf > 
{
    let mut orf_seq: Vec < String > = Vec::new();
    let mut orf_start: usize = 0;
    //let mut orf_dict: FxHashMap< String, Orf > = FxHashMap::default();
    let mut orf_dict: Vec < Orf > = Vec::new();
    for (i,amino_acid) in prot_seq.clone().iter().enumerate()
    {
        if amino_acid == "$"
        {
            if i > 0
            {
                let n_start = n_pos_list[0];
                let n_end = n_pos_list[i] + 2;

                orf_start = orf_start + 1;
                let orf_end = i;
                orf_seq.push("$".to_string());
                let orf_obj: Orf = Orf::new(orf_seq.clone(), orf_start, orf_end, frame.clone(), orf_seq[0].clone(), n_start, n_end);
                orf_dict.push(orf_obj);
                break;
            }
            else
            {
                break;
            }
        }
        orf_seq.push(amino_acid.clone().to_string());
    }

    let mut orf_seq: Vec < String > = Vec::new();
    let mut orf_start: usize = 0;

    for (i, amino_acid) in prot_seq.clone().iter().enumerate()
    {
        if orf_seq.len() == 0
        {
            if amino_acid == "M"
            {
                orf_start = i;
                orf_seq.push(amino_acid.clone().to_string());
            }
            continue;
        }

        if amino_acid == "$"
        {
            if orf_seq.len() > 0
            {
                let n_start = n_pos_list[orf_start];
                let n_end = n_pos_list[i] + 2;
                
                orf_start += 1;
                let orf_end = i;
                orf_seq.push("$".to_string());
                let orf_obj: Orf = Orf::new(orf_seq.clone(), orf_start, orf_end, frame.clone(), orf_seq[0].clone(), n_start, n_end);
                orf_dict.push(orf_obj);
                orf_seq = Vec::new();
                orf_start = i + 1;
            }
            continue;
        }
        orf_seq.push(amino_acid.clone().to_string());
    }

    return orf_dict;
}



/// Function that allows to retrieve key,value from fasta header.
/// This functions takes a fasta header/description as defined from bio and
/// a key value (needs to be present in the header) and returns the value of the 
/// key. 
/// 
/// # Arguments 
/// 
/// * fasta_description - the fasta header 
/// * value             - the key for which we want to return the value of
/// 
/// # Examples 
/// 
/// let fasta_description = "ENSG90909 GN=TA_MERE ENS=EE PE=1 SV=1";
/// get_value_from_fasta_description(fasta_description, "GN");
///
pub fn get_value_from_fasta_description(fasta_description: &str, value: &str) -> String
{
    let tmp: Vec < &str > = fasta_description.split(" ").collect();
    let mut to_return_value: String = String::new();
    for i in tmp.iter()
    {
        let tmp_i: Vec < &str > = i.split("=").collect();
        if tmp_i[0] == value
        {
            to_return_value = tmp_i[1].to_string();
            break;           
        }
    }
    return to_return_value;
}



/// Function that checks whether ORF contains mutation by using the coordinates of the
/// ORF and the coordinate of the mutation
/// 
/// # Arguments 
/// 
/// * n_start - nucleotide sequence start of the ORF 
/// * n_end   - nucleotide sequence end of the ORF 
/// * m_pos   - nucleotide sequence position of the mutation
/// 
/// # Example 
/// let n_start: usize = 100;
/// let n_end: usize = 200;
/// let mut_pos: usize = 150;
/// assert_eq!(true, orf_contains_mutation(&n_start, &n_end, &mut_pos));
pub fn orf_contains_mutation(n_start: &usize, n_end: &usize, mut_pos: &usize) -> bool
{
    if n_start <= mut_pos && n_end >= mut_pos 
    {
        return true;
    }
    else
    {
        return false;
    }
}


/// Function that reads a fasta file into a HashMap
/// where the key is the fasta header and the value the
/// sequence
/// 
/// # Arguments 
/// 
/// * fasta_file  - the path to the fasta file to read 
/// 
/// # Example
/// 
/// let fasta_file: &str = "./test/fasta_file.fasta";
/// read_fasta_file(fasta_file);
pub fn read_fasta_file(fasta_file: &str) -> FxHashMap < String, String >
{
    let mut fasta_hash: FxHashMap < String, String > = FxHashMap::default();
    let reader = fasta::Reader::from_file(fasta_file).unwrap();

    for results in reader.records()
    {
        let record = &results.unwrap();
        let fasta_header: &str = record.id();
        let fasta_sequence: &str = str::from_utf8(record.seq()).unwrap();
        fasta_hash.insert(fasta_header.to_string(), fasta_sequence.to_string());
    }
    return fasta_hash;
}



/// Function that reads a HashMap containing as key transcript_ids and as value
/// the corresponding sequence of the transcript and returns the sequence of 
/// interest. 
/// 
/// # Arguments
/// 
/// * transcripthash  - HashMap key: transcript_id; value -> transcript sequence
/// * transcript_id   - transcript id
/// * start           - start position (0-based) of the sequence to be extracted
/// * end             - end position of the sequence to be extracted
pub fn extract_sequence_from_hash(transcript_hash: &FxHashMap<String, String>, transcript_id: &str, start: usize, end:usize) -> String
{
    let transcript_seq: &String = &transcript_hash[transcript_id];
    let seq: String = transcript_seq[start..end].to_string();

    return seq;
}



/// Function that returns the amino acid at a given position
/// 
/// # Arguments 
/// 
/// * seq             - amino acid sequence 
/// * mutation_on_pep - position of the amino acid to be extracted (0-based)
/// 
/// # Examples 
/// 
/// let aa_seq: Vec < String > = ["T".to_string(), "A".to_string(), "M",to_string(),"E".to_string(), "R".to_string(), "E".to_string()];
/// let mutation_on_pep: usize = 2;
/// assert_eq!(get_mutated_aa(&aa_seq, mutation_on_pep), "M".to_string());
pub fn get_mutated_aa(aa_seq: &Vec < String >, mutation_on_pep: usize) -> String 
{
    return aa_seq[mutation_on_pep].to_string();
}



/// Function that extracts from an Orf object the start, end and the position of a
/// given mutation on a amino acid sequence.
/// 
/// # Arguments
/// 
/// * obj_orf                  - Orf object containing information regarding an ORF
/// * mutation_pos_on_sequence - the position of the mutation on the nuc sequence
pub fn get_mutated_aa_pos(obj_orf: &Orf, mutation_pos_on_sequence: usize) -> (usize, usize, usize)
{
    let mut mut_pos: usize = 0;
    let mut i:usize = obj_orf.n_start;

    while i < obj_orf.n_end 
    {
        if i <= mutation_pos_on_sequence && i+2 >= mutation_pos_on_sequence 
        {
             break // 0 based mut_pos!!!
        }
        i += 3;
        mut_pos += 1;
    }
    return (mut_pos, i, i+3);
}


/// This function creates the extension of a ORF id 
/// with the AAposAA nomenclature
/// 
/// # Arguments 
/// 
/// * obj_orf  - Orf struct containing ORF information
/// * mutation_pos_on_sequence  - the mutation position on the nuc sequence
/// * transcript_hash           - HashMap where keys are transcript_ids and values their corresponding sequence
/// * transcript_id             - a transcript_id (e.g. ENST0000014003)
pub fn update_with_mutation_pos(obj_orf: &Orf, mutation_pos_on_sequence: usize, transcript_hash: &FxHashMap<String, String > , transcript_id: &str) -> String
{
    let codon_hash: FxHashMap < &str, &str > = crate::utils::codon_hash();
    let mut i:usize = obj_orf.n_start;
    
    let mut extension: String = String::new();
    

    let (mut_pos, n_start, n_end) = get_mutated_aa_pos(&obj_orf, mutation_pos_on_sequence.clone()); 
    let mutated_aa = get_mutated_aa(&obj_orf.orf_seq, mut_pos);
    let codon = extract_sequence_from_hash(&transcript_hash, transcript_id, n_start, n_end);
    if mut_pos == 0
    {
        if &codon_hash[codon.as_str()] != &mutated_aa
        {
            extension = format!("_{}{}{}", &codon_hash[codon.as_str()], mut_pos+1, mutated_aa);
        }
        else {
            extension = "".to_string();
        }
    }
    else if mut_pos == obj_orf.orf_len
    {
        //extension = format!("_{}_{}", n_count, &obj_orf.orf_seq[n_count-1]);
        if &codon_hash[codon.as_str()] != &mutated_aa
        {
            extension = format!("_{}{}{}", &codon_hash[codon.as_str()], mut_pos+1, mutated_aa);
        }
        else {
            extension = "".to_string();
        }
    }
    else
    {
        //extension = format!("_{}_{}", n_count+1, &obj_orf.orf_seq[n_count]);
        if &codon_hash[codon.as_str()] != &mutated_aa
        {
            extension = format!("_{}{}{}", &codon_hash[codon.as_str()], mut_pos+1, mutated_aa);
        }
        else {
            extension = "".to_string();
        }
    }
    
    return extension;
}


/// This is the main function to run the Translation of nucleotide sequences to ORFs
/// 
/// # Arguments 
/// 
/// * fasta_file     - the path to the fasta file containing the transcript sequences to be translated
/// * output_file    - the path to the output file where to write the discovered ORFs
/// * prot_gene_hash - hashmap containing as key gene_name and values a struct with protein_id and protein_name
/// * min_len        - minimum length of the ORF to be considered. This is turned OFF if mutated_only is true
/// * write_missing  - whether to write fasta header for transcripts where ORF could not be determined because they contained NAs
/// * mutated_only   - whether to write all discovered ORFs or only the ones baring a mutation as specified by the fasta header line.
/// 
///  
pub fn orf_seeker_main(mutated_fasta_file: &str, normal_fasta_file: &str,  output_file: &str, min_len: usize, write_missing: bool, mutated_only: bool, all_codons: bool)
{
    
    let transcript_hash = read_fasta_file(&normal_fasta_file);

    // Open file to write 
    let mut output_buffer: File = File::create(output_file).expect("ERROR: cannot open file to write");
    
    let reader = fasta::Reader::from_file(mutated_fasta_file).unwrap();

    let codon_hash = codon_hash();

    let mut tmp_set: HashSet < String > = HashSet::new();

    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_header: &str = results_data.id();
        let tmp: Vec < &str > = fasta_header.split("_").collect();
        let mutation_pos_on_sequence: usize = tmp[2].parse::<usize>().unwrap() -1; // We need to substract 1 here as the mutation_pos_on_sequence in the 
        // fasta header is 1 based but here it needs to be 0 based to be able to be properly compared with n_start and n_end which them are 0 based.

        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();

        let f1: u32 = 1;
        let f2: u32 = 2;
        let f3: u32 = 3;

        let (f1_prot_seq, f1_n_flag, f1_n_pos_list) = frame_iterate(fasta_sequence, f1, &codon_hash);
        let (f2_prot_seq, f2_n_flag, f2_n_pos_list) = frame_iterate(fasta_sequence, f2, &codon_hash);
        let (f3_prot_seq, f3_n_flag, f3_n_pos_list) = frame_iterate(fasta_sequence, f3, &codon_hash);
        
        let all_n_flag: u32 = f1_n_flag + f2_n_flag + f3_n_flag; 
        if all_n_flag > 0
        {   
            if ! write_missing
            {
                continue;
            }
            else 
            {
                let outline: String = ">".to_string() + fasta_header + ":missing_nucleotides";
                write!(output_buffer, "{}\n", outline).expect("Cannot write to file");  
            }            
            continue;
        }

        let mut f1_orf_dict = orf_seeker(f1_prot_seq, f1, f1_n_pos_list);
        let f2_orf_dict = orf_seeker(f2_prot_seq, f2, f2_n_pos_list);
        let f3_orf_dict = orf_seeker(f3_prot_seq, f3, f3_n_pos_list);
        
        f1_orf_dict.extend(f2_orf_dict);
        f1_orf_dict.extend(f3_orf_dict);
        

        for obj_orf in f1_orf_dict.into_iter()
        {
            if !all_codons && obj_orf.start_codon != "M".to_string() { continue; }

            let mut key = fasta_header.to_string() + "_" + &obj_orf.orf_name.as_str();
            if mutated_only
            {
                
                if obj_orf.orf_seq.len() < min_len { continue; } // To filter length of ORFS that are outputed

                if orf_contains_mutation(&obj_orf.n_start, &obj_orf.n_end, &mutation_pos_on_sequence)
                {
                    let t: Vec < &str > = fasta_header.split("_").collect();
                    
                    let transcript_id: &str = t[0];
                    let mutation_type: &str = t[3];
                    if mutation_type == "INS"
                    {
                        let (mutation_pos_aa, _, _) = get_mutated_aa_pos(&obj_orf, mutation_pos_on_sequence);
                        let mutated_aa: String = get_mutated_aa(&obj_orf.orf_seq, mutation_pos_aa);
                        key += format!("_INS{}{}", mutation_pos_aa, mutated_aa).as_str();
                    }
                    else if mutation_type == "DEL"
                    {
                        let (mutation_pos_aa, _, _) = get_mutated_aa_pos(&obj_orf, mutation_pos_on_sequence);
                        let mutated_aa: String = get_mutated_aa(&obj_orf.orf_seq, mutation_pos_aa);
                        key += format!("_DEL{}{}", mutation_pos_aa, mutated_aa).as_str();
                    }
                    else
                    {
                        let extension: String = update_with_mutation_pos(&obj_orf, mutation_pos_on_sequence, &transcript_hash, transcript_id); 
                        if extension == "".to_string()
                        {
                            continue;
                        }
                        else {
                            key += extension.as_str();
                        }
                            
                    }
                    
                    //key += update_with_mutation_pos(&obj_orf, mutation_pos_on_sequence.clone()).as_str();
                    
                    
                    let outline: String = format!(">{}", key);

                    if tmp_set.contains(&outline) { continue; }
                    tmp_set.insert(outline.clone());
                    write!(output_buffer, "{}\n", outline).expect("Cannot write to file");
                    let orf_seq = &obj_orf.orf_seq[..obj_orf.orf_len].join("");
                    write!(output_buffer,"{}\n", orf_seq).expect("Cannot write to file");
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if obj_orf.orf_seq.len() < min_len { continue; } // To filter length of ORFS that are outputed


                let outline: String = format!(">{}", key);

                if tmp_set.contains(&outline) { continue; }
                tmp_set.insert(outline.clone());
                write!(output_buffer, "{}\n", outline).expect("Cannot write to file");
                let orf_seq = obj_orf.orf_seq.join("");
                write!(output_buffer,"{}\n", orf_seq).expect("Cannot write to file");
            }
        }
    }

}

