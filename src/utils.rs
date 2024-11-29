
use std::str;

use rustc_hash::FxHashMap;



pub fn codon_hash() -> FxHashMap < &'static str, &'static str >
{
    let mut codon: FxHashMap < &str, &str > = FxHashMap::default();
    
    codon.insert("GTT", "V");
    codon.insert("GTC", "V");
    codon.insert("GTA", "V");
    codon.insert("GTG", "V");

    codon.insert("GCT", "A");
    codon.insert("GCC", "A");
    codon.insert("GCA", "A");
    codon.insert("GCG", "A");

    codon.insert("GAT", "D");
    codon.insert("GAC", "D");

    codon.insert("GAA", "E");
    codon.insert("GAG", "E");

    codon.insert("GGT", "G");
    codon.insert("GGC", "G");
    codon.insert("GGG", "G");
    codon.insert("GGA", "G");

    codon.insert("TTT", "F");
    codon.insert("TTC", "F");

    codon.insert("TTA", "L");
    codon.insert("TTG", "L");

    codon.insert("TCT", "S");
    codon.insert("TCC", "S");
    codon.insert("TCA", "S");
    codon.insert("TCG", "S");

    codon.insert("TAT", "Y");
    codon.insert("TAC", "Y");
    codon.insert("TAA", "$"); // stop codon
    codon.insert("TAG", "$"); // stop codon

    codon.insert("TGT","C");
    codon.insert("TGC","C");
    codon.insert("TGA","$"); // Stop codon
    codon.insert("TGG","W");

    codon.insert("CTT", "L");
    codon.insert("CTC", "L");
    codon.insert("CTA", "L");
    codon.insert("CTG", "L");

    codon.insert("CCT", "P");
    codon.insert("CCC", "P");
    codon.insert("CCA", "P");
    codon.insert("CCG", "P");

    codon.insert("CAT", "H");
    codon.insert("CAC", "H");

    codon.insert("CAA", "Q");
    codon.insert("CAG", "Q");

    codon.insert("CGT", "R");
    codon.insert("CGC", "R");
    codon.insert("CGA", "R");
    codon.insert("CGG", "R");

    codon.insert("ATT", "I");
    codon.insert("ATC", "I");
    codon.insert("ATA", "I");

    codon.insert("ATG", "M");

    codon.insert("ACT", "T");
    codon.insert("ACC", "T");
    codon.insert("ACA", "T");
    codon.insert("ACG", "T");

    codon.insert("AAT", "N");
    codon.insert("AAC", "N");

    codon.insert("AAA", "K");
    codon.insert("AAG", "K");

    codon.insert("AGT", "S");
    codon.insert("AGC", "S");
    codon.insert("AGA", "R");
    codon.insert("AGG", "R");

    return codon 
}