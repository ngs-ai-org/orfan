/*use std::str;
use std::fs::File;
use std::io::{BufRead,BufReader};
use std::path::Path;

use rustc_hash::FxHashMap;
use core::panic;
use faimm::IndexedFasta;


/// Function that reads bed 3 file
/// generated from fetchICGCData.py
/// 
/// # Arguments: 
/// * `filename` - file path to be read. 
/// # Returns: 
///  FxhashMap < String, Vec < MutBedStruct >>

pub fn read_mut_bed_file (filename: &str) -> FxHashMap < String, Vec < MutBedStruct > >
{
    let my_file = Path::new(&filename);
    let input   = File::open(&my_file).expect("Couldn't open file");
    let buffer  = BufReader::new(input);

    let mut tmp_map: FxHashMap<String, Vec < MutBedStruct >> = FxHashMap::default();    
    
    for line in buffer.lines()
    {
        let l = line.expect("Error: Could not read the line!");
        let e: Vec < &str > = l.split("\t").collect();

        if e[0] == "#chr"
        {
            continue;
        }
        let entry_struct: MutBedStruct = MutBedStruct {
            chrom   : (*e[0]).to_string(),
            start   : (*e[1]).parse::<i64>().unwrap(),
            end     : (*e[2]).parse::<i64>().unwrap(),
            mutation_id : (*e[3]).to_string(),
            mutation : (*e[8]).to_string(),
            ref_allele : (*e[9]).to_string(),
            variant    : (*e[10]).to_string(),
            mutation_type : (*e[4]).to_string(),
            trans_mutation_type: (*e[11]).to_string()
        };
        let transcript_id = (*e[7]).to_string();
        if !tmp_map.contains_key(&transcript_id)
        {
            let mut tmp_vec: Vec < MutBedStruct > = Vec::new();
            tmp_vec.push(entry_struct.clone());
            tmp_map.insert(transcript_id, tmp_vec);
        }
        else
        {   
            let mut tmp_vec: Vec < MutBedStruct > = tmp_map[&transcript_id].clone();
            tmp_vec.push(entry_struct.clone());
            *tmp_map.get_mut(&transcript_id).unwrap() = tmp_vec; 
        }
    }

    return tmp_map;
}

*/