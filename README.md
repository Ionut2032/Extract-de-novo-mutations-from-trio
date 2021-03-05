# Extract-de-novo-mutations-from-trio
Project completed during my MSc in Bioinformatics at the University of Birmingham, United Kingdom, 2019-2020.

## Data pre-processing
1. **QC & FILTERING** 
   * pass threshold lod≥3 
   * exclude missing variants
   * exclude multi-allelic variants
   * exclude low no of reads - only keep ≥30
   * exclude low PHRED scores - only keep≥20

## Data analysis
* uses ***vcfR*** R package
1. **EXTRACT GENOTYPE** (reference allele vs. alternative allele e.g. 0/1)
2. **EXTRACT DNA BASES** (reference allele vs. alternative allele)
3. **FIND THE CHILD** (count the no of incompatible variants)
4. **EXTRACT DE NOVO VARIANTS FROM THE CHILD** 
5. **FURTHER RESEARCH**
   * validation of de novo variants 
    * **computational** e.g. calculate the probability of the mutation to be present in a family trio 
    * **experimental** e.g. SANGER sequencing
   * clinical impact of de novo varints (UCSC Genome Browser)
