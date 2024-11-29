# Open Reading Frame analysis (ORFan) toolkit

:heavy_exclamation_mark: ORFan is deprecated! Please use the version present here.

![](data/proteinator.jpeg)

ORFan is a suite of tools used for generating transcripts baring mutations or performing in-silico translation of cDNA sequences. It is written in Rust. 


ORFan has 2 main functions: 
1. Generation of in-silico mutated transcripts
2. Generation of in-silico mutated peptides

# transcript-mutator 

```bash 
Usage: orfan transcript-mutator [OPTIONS] --vcf-file <vcf_file> --gtf-file <gtf_file> --reference-genome <reference_genome> --output-file <output_file>

Options:
  -f, --vcf-file <vcf_file>
          input vcf file containing mutations
  -g, --gtf-file <gtf_file>
          path to gtf annotation file
  -r, --reference-genome <reference_genome>
          reference genome in fasta format
  -o, --output-file <output_file>
          Path to output file
  -s, --strand-aware-off
          Whether to turn OFF strand aware sequence report. NOT RECOMMENDED TO BE TURNED OFF if the transcripts are going to be used to generate ORFs
  -h, --help
          Print help
  -V, --version
          Print version
```

`transcript-mutator` requires you provide a vcf file containing coding mutations (mutations falling onto exonic regions), a gtf annotation file (can be found [here](https://www.gencodegenes.org/)), a reference genome in fasta format and an output file where to write the mutated transcripts. 


# orf-seeker 


```bash 
Usage: orfan orf-seeker [OPTIONS] --input-file <input_fasta_file> --normal-fasta-file <normal_fasta_file> --output-file <output_fasta_file>

Options:
  -i, --input-file <input_fasta_file>
          Input fasta file
  -n, --normal-fasta-file <normal_fasta_file>
          Fasta file containing transcript sequences as obtained by seqtk
  -o, --output-file <output_fasta_file>
          Output fasta file
  -m, --mutated-only
          Whether to only write ORFs containing mutation
  -l, --min-len <min_orf_len>
          Minimum length of the ORF [default: 6]
  -a, --all-codons
          Whether to print ORFs starting with any codon
  -w, --write-missing
          Whether to write header line of missing ORFs in output
  -h, --help
          Print help
  -V, --version
          Print version
```

`orf-seeker` takes a set of transcripts and performs in-silico translation in all 3 frames. The output are peptide sequences. `orf-seeker` is based on `tama-orf-seeker` (PUT GITHUB PAGE).