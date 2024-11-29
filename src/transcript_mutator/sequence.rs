use core::panic;
use std::str;


/// Function that compares to string slices
/// Returns true if identical
/// Returns false if different
/// 
/// # Arguments: 
/// * x - slice x 
/// * y - slice y 
/// 
/// # Examples
/// 
/// ```
/// let x: &str = "ABCD";
/// let y: &str = "ABCD";
/// assert_eq!(true, compare_seq_slices(x,y));
/// let xx: &str = "AAAA";
/// let yy: &str = "DDDD";
/// assert_ne!(xx, yy);
/// ```
pub fn compare_seq_slices(x: &str, y: &str) -> bool
{
    if x == y
    {
        return true
    }
    return false;
}


/// Function that takes a nucleotide sequence and mutates at the specific point mutation 
/// 
/// # Arguments
/// 
/// * `sequence`         - Nucleotide sequence 
/// * `start`            - Start position of the sequence in the genome 
/// * `mutation_pos`     - Position of the mutation in the genome
/// * `mutated_allele`   - The mutated allele at position
/// * `reference_allele` - The reference allele at the position
/// 
/// # Example 
/// 
/// ```
/// let sequence: String = "CCTAATTCGAGGTTCCAAA".to_string();
/// let start: usize     = 10;
/// let mutation_pos: usize = 15;
/// let mutated_allele: &str = "A";
/// let reference_allele: &str ="T";
/// 
/// point_mutation_mutator(sequence, start, mutation_pos, mutated_allele, reference_allele);
/// 
/// ```
pub fn point_mutation_mutator(sequence: &String,  mutation_pos_on_sequence:usize, mutated_allele:&str, reference_allele:&str) -> String
{
    let mut mutated_sequence: String = sequence.to_string();
    if mutation_pos_on_sequence == 0
    {
        
        if &sequence[..mutation_pos_on_sequence + 1] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele {} does not match {}!!!!", reference_allele,&sequence[..mutation_pos_on_sequence +1]);
        }
        mutated_sequence.replace_range((mutation_pos_on_sequence)..mutation_pos_on_sequence+1, &mutated_allele);
    }
    else if mutation_pos_on_sequence == sequence.len()-1
    {
        if &sequence[mutation_pos_on_sequence..] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele {} does not match {}!!!!", reference_allele,&sequence[mutation_pos_on_sequence..]);
        }
        mutated_sequence.replace_range(mutation_pos_on_sequence.., &mutated_allele);    
    }
    else
    {
        if &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+1] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele {} does not match {}!!!!", reference_allele, &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+1]);
        }
        mutated_sequence.replace_range(mutation_pos_on_sequence..mutation_pos_on_sequence+1, &mutated_allele);    
    }
    return mutated_sequence;
}



/// Function that takes a nucleotide sequence and mutates at the specific point mutation 
/// 
/// # Arguments
/// 
/// * `sequence`         - Nucleotide sequence 
/// * `start`            - Start position of the sequence in the genome 
/// * `mutation_pos_on_sequence`     - Position of the mutation in the genome
/// * `mutated_allele`   - The mutated allele at position
/// * `reference_allele` - The reference allele at the position
/// 
/// # Example 
/// 
/// ```
/// let sequence: String= "ATCGATTGATCGATCGATCG".to_string();
/// let start:usize = 10;
/// let mutation_pos:usize = 15; // 1 based
/// let mutated_allele:&str = "GG";
/// let reference_allele:&str = "AT";
/// let expected_mutated_sequence: String = "ATCGGTTGATCGATCGATCG".to_string();
/// let mutation_pos_on_sequence: usize = (mutation_pos - start - 1) + 0;
/// 
/// multi_mutation_mutator(&sequence, mutation_pos_on_sequence, mutated_allele, reference_allele);
/// 
/// ```
pub fn multi_mutation_mutator(sequence: &String,  mutation_pos_on_sequence:usize, mutated_allele:&str, reference_allele:&str) -> String
{
    let mut mutated_sequence: String = sequence.to_string();
    let mutation_len: usize = reference_allele.len();
    if mutation_pos_on_sequence == 0
    {
        
        if &sequence[..mutation_pos_on_sequence + mutation_len] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele does not match!!!!");
        }
        mutated_sequence.replace_range((mutation_pos_on_sequence)..mutation_pos_on_sequence+mutation_len, &mutated_allele);
    }
    else if mutation_pos_on_sequence == sequence.len()-1
    {
        if &sequence[mutation_pos_on_sequence..] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele does not match!!!!");
        }
        mutated_sequence.replace_range(mutation_pos_on_sequence.., &mutated_allele);    
    }
    else
    {
        if &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+mutation_len] != reference_allele
        {
            panic!("There is something very wrong here! Reference allele does not match!!!!");
        }
        mutated_sequence.replace_range(mutation_pos_on_sequence..mutation_pos_on_sequence+mutation_len, &mutated_allele);    
    }
    return mutated_sequence;
}




///Function that takes a genomic sequence and creates a given insertion of nucleotides at position X
/// 
/// # Arguments
/// 
/// * `sequence`         - Nucleotide sequence 
/// * `start`            - Start position of the sequence in the genome 
/// * `mutation_pos`     - Position of the mutation in the genome
/// * `mutated_allele`   - The mutated allele at position
/// * `reference_allele` - The reference allele at the position
/// 
/// # Example 
/// 
/// ```
/// let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
/// let start: usize     = 10;
/// let mutation_pos: usize = 14;
/// let mutated_allele: &str = "ATGCTCGA"; 
/// let reference_allele: &str ="A";
/// let mutated_sequence: String = insertion_mutator(&sequence, start, mutation_pos, mutated_allele, reference_allele);
/// let expected: String = "CCTAATGCTCGAATGCTCGAGGTTCCAAA".to_string();
/// assert_eq!(expected, mutated_sequence);
///    
/// let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
/// let start: usize     = 10;
/// let mutation_pos: usize = 10;
/// let mutated_allele: &str = "CTA"; 
/// let reference_allele: &str ="C";
/// let mutated_sequence: String = insertion_mutator(&sequence, start, mutation_pos, mutated_allele, reference_allele);
/// let expected: String = "CTACTAATGCTCGAGGTTCCAAA".to_string();
/// assert_eq!(expected, mutated_sequence);
/// 
/// let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
/// let start: usize     = 10;
/// let mutation_pos: usize = 30;
/// let mutated_allele: &str = "ACAAA"; 
/// let reference_allele: &str ="A";
/// let mutated_sequence: String = insertion_mutator(&sequence, start, mutation_pos, mutated_allele, reference_allele);
/// let expected: String = "CCTAATGCTCGAGGTTCCAAACAAA".to_string();
/// assert_eq!(expected, mutated_sequence);
/// ```
pub fn insertion_mutator(sequence: &String, mutation_pos_on_sequence: usize, mutated_allele: &str, reference_allele:&str) -> String{
    
    
    let mut mutated_sequence: String = String::new();
    if mutation_pos_on_sequence == 0
    {
        let tmp: &str = &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+1];
        if !compare_seq_slices(tmp, reference_allele)
        {
            panic!("First nuc of insertion is: {} and reference is: {}. There is something wrong.", tmp, reference_allele);
        }
        let slice1: &str = &sequence[mutation_pos_on_sequence+1..];
        mutated_sequence = mutated_allele.to_string() + slice1;
    }
    else if mutation_pos_on_sequence == sequence.len()-1
    {
        let tmp: &str = &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+1];
        if !compare_seq_slices(tmp, reference_allele)
        {
            panic!("First nuc of insertion is: {} and reference is: {}. There is something wrong.", tmp, reference_allele);
        }
        let slice1: &str = &sequence[..mutation_pos_on_sequence];
        mutated_sequence = slice1.to_string() + mutated_allele;
    }
    else 
    {
        let tmp: &str = &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence+1];
        if !compare_seq_slices(tmp, reference_allele)
        {
            panic!("First nuc of insertion is: {} and reference is: {}. There is something wrong.", tmp, reference_allele);
        }

        let slice1: &str = &sequence[..mutation_pos_on_sequence];
        let after: &str = &sequence[mutation_pos_on_sequence+1..];
        mutated_sequence = slice1.to_string() + mutated_allele + after;
    }
    return mutated_sequence;
}



/// Function that takes a genomic sequence and performs a deletion of X nucleotides 
/// 
/// # Arguments
/// 
/// * `sequence`         - Nucleotide sequence 
/// * `start`            - Start position of the sequence in the genome 
/// * `mutation_pos`     - Position of the mutation in the genome
/// * `mutated_allele`   - The mutated allele at position
/// * `reference_allele` - The reference allele at the position
/// 
/// # Example
/// 
/// ```
/// let sequence: String = "CCTAATGCTCGAGGTTCCAAA".to_string();
/// let start: usize     = 10;
/// let mutation_pos: usize = 14;
/// let mutated_allele: &str = "-"; 
/// let reference_allele: &str ="ATGCTCGA";
/// let mutated_sequence: String = deletion_mutator(&sequence, start, mutation_pos, reference_allele);
/// let expected_sequence: String = "CCTAGGTTCCAAA".to_string();
/// assert_eq!(mutated_sequence, expected_sequence);
/// ```

pub fn deletion_mutator(sequence: &String, mutation_pos_on_sequence: usize, reference_allele:&str) -> String 
{
    let mutation_len: usize = reference_allele.len();
    
    let mutated_slice: &str = &sequence[mutation_pos_on_sequence..mutation_pos_on_sequence + mutation_len];
    if !compare_seq_slices(reference_allele, mutated_slice)
    {
        panic!("The deletion is: {} but the extracted slice to remove is: {}. There is something wrong!", reference_allele, mutated_slice);
    }
    
    
    let mut mutated_sequence: String = String::new();
    
    if mutation_pos_on_sequence == 0
    {
        let slice1: &str = &sequence[0..1];
        let slice2: &str = &sequence[mutation_len..];// -1 is required here  because the first base on the ref allele for deletions is the ref pos base and the deletion is anything coming after it.
        mutated_sequence = slice1.to_string() + slice2;
    }
    else
    {
        let slice1: &str = &sequence[..mutation_pos_on_sequence+1];
        let slice2: &str = &sequence[mutation_pos_on_sequence+1 + mutation_len-1..];
        mutated_sequence = slice1.to_string() + slice2;

    }
    return mutated_sequence;
}



pub fn strand_aware_seq_orientation(sequence: &String) -> String
{
    let mut new_sequence: String = String::new();
    let tmp_seq: String = sequence.chars().rev().collect::<String>();
    for i in tmp_seq.chars()
    {
        if i == 'A'
        {
            new_sequence += "T";
        }
        else if i == 'T'
        {
            new_sequence += "A";
        }
        else if i == 'C'
        {
            new_sequence += "G";
        }
        else if i == 'G'
        {
            new_sequence += "C";
        }
        else
        {
            panic!("Found {} which is not a nucleotide!!!", i);    
        }
    }
    
    return new_sequence;
}
