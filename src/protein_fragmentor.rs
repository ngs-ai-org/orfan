use std::str;
use std::fs::File;
use std::io::Write;
use rustc_hash::{FxHashMap, FxHasher};
use bio::io::{self, fasta};
use std::io::{BufRead,BufReader};
use std::path::Path;



pub fn get_cut_sites(protein_sequence: &str) -> Vec < usize >
{
    let mut cut_sites : Vec < usize > = vec![0];
    
    
    for (i,j) in protein_sequence.chars().into_iter().enumerate()
    {
        if j == 'K'
        {
            if i == protein_sequence.len() - 1
            {
                cut_sites.push(i+1);    
            }
            else if i+2 > protein_sequence.len()
            {
                if &protein_sequence[i+1..] != "P"
                {
                    cut_sites.push(i+1);
                }
            }
            else 
            {
                if &protein_sequence[i+1..i+2] != "P"
                {
                    cut_sites.push(i+1);
                }
            }   
        }
        else if j == 'R'
        {
            if i == protein_sequence.len() - 1
            {
                cut_sites.push(i+1);    
            }
            else if i+2 > protein_sequence.len()
            {
                if &protein_sequence[i+1..] != "P"
                {
                    cut_sites.push(i+1);
                }
            }
            else 
            {
                if &protein_sequence[i+1..i+2] != "P"
                {
                    cut_sites.push(i+1);
                }
            }
            
        }
    }
    if cut_sites[cut_sites.len()-1] != protein_sequence.len()
    {
        cut_sites.push(protein_sequence.len());
    }   
    return cut_sites;
}

/// Function that trypsinates a protein sequence
/// The function follows the standard cutting of trypsin which
/// cuts after a K or R but never when K or R are followed by a P.
/// In other words, if we take the following protein sequence 
/// MCNKMCNKPMCNRMCNR, trypsin should fragment it into three:
/// MCNK MCNKPMCNR and MCNR.
/// 
/// # Arguments 
/// 
/// * protein_seq - a protein sequence to trypsinate
/// 
/// # Examples 
/// 
/// let protein_sequence: &str = "MCNKMCNKPMCNRMCNR";
/// let observed_result = trypsinate_protein_sequence(protein_sequence);
/// println!("{:?}", observed_result);
/// let protein_sequence: &str = "MCNKMCNKPMCNRMCNRPMCND";
/// let observed_result = trypsinate_protein_sequence(protein_sequence);
/// println!("{:?}", observed_result);
pub fn trypsinate(protein_sequence: &str, missed_cleavage: usize) -> Vec < String >
{    
    let mut peptides: Vec < String > = Vec::new();
    let cut_sites: Vec<usize> = get_cut_sites(protein_sequence);

    
    if cut_sites.len() > 2
    {
        if missed_cleavage == 0
        {
            for i in 0..cut_sites.len()-1
            {
                if i == cut_sites.len()-1
                {
                    println!("{}:{}", cut_sites[i], cut_sites[i-1]);
                    //peptides.push(protein_sequence[*j..].to_string());
                }
                else
                {
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+1]].to_string());
                }
            }
        }
        else if missed_cleavage == 1
        {
            for i in 0..cut_sites.len()-2
            {
                if i == cut_sites.len()-2
                {
                    println!("{}:{}", cut_sites[i], cut_sites[i-1]);
                    //peptides.push(protein_sequence[*j..].to_string());
                }
                else
                {
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+1]].to_string());
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+2]].to_string());
                }
            }
            peptides.push(protein_sequence[cut_sites[cut_sites.len()-2]..cut_sites[cut_sites.len()-1]].to_string());
        }
        else if missed_cleavage == 2
        {
            for i in 0..cut_sites.len()-3
            {
                if i == cut_sites.len()-2
                {
                    println!("{}:{}", cut_sites[i], cut_sites[i-1]);
                    //peptides.push(protein_sequence[*j..].to_string());
                }
                else
                {
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+1]].to_string());
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+2]].to_string());
                    peptides.push(protein_sequence[cut_sites[i]..cut_sites[i+3]].to_string());
                }
            }
            peptides.push(protein_sequence[cut_sites[cut_sites.len()-3]..cut_sites[cut_sites.len()-2]].to_string());
            peptides.push(protein_sequence[cut_sites[cut_sites.len()-3]..cut_sites[cut_sites.len()-1]].to_string());
            peptides.push(protein_sequence[cut_sites[cut_sites.len()-2]..cut_sites[cut_sites.len()-1]].to_string());
        }
        
        
    }
    return peptides;
}






pub fn write_to_fasta(mut output_buffer: &File, header: String, trypsinated_peptides: Vec < String >)
{
    for (i,pep) in trypsinated_peptides.iter().enumerate()
{
        let new_header = header.to_string() + " FRAG=" + i.to_string().as_str();
        write!(output_buffer,">{}\n", new_header).expect("Cannot write to file");
        write!(output_buffer, "{}\n", pep).expect("Cannot write to file");
    }
}


pub fn protein_fragmentation_main(
    fasta_file: String,
    enzyme: String,
    output_file: String,
    missed_cleavage: usize
    ) 
{
    let reader = fasta::Reader::from_file(fasta_file).unwrap();

    // Open file to write
    let output_buffer: File = File::create(output_file).expect("Cannot open file");

    for result in reader.records()
    {
        let results_data = &result.unwrap();
        let fasta_header: &str = results_data.id();
        let fasta_sequence: &str = str::from_utf8(results_data.seq()).unwrap();
        let mut trypsinated_peptides: Vec < String > = Vec::new();
        match enzyme.as_str() 
        {
            "trypsin" => trypsinated_peptides = trypsinate(fasta_sequence, missed_cleavage),
            "tamere"  => eprintln!("Not yet implemented"),
            _         => panic!("WARNING: the enzyme you specified either does not exist or is not supported by the tool")
        }
        
        let header: String = fasta_header.to_string();
        write_to_fasta(&output_buffer, header, trypsinated_peptides);
    }
}