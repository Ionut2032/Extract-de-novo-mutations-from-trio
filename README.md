# Extract-de-novo-mutations-from-trio
GENOMICS: finding the child from a family trio and extracting de novo mutations using vcfR R package.

Project completed during my MSc in Bioinformatics at the University of Birmingham, United Kingdom.

## Data pre-processing
1. **QC & FILTERING** 
   * pass threshold lod≥3 
   * exclude missing variants
   * exclude multi-allelic variants
   * exclude low no of reads - only keep ≥30
   * exclude low PHRED scores - only keep≥20

## Data analysis
1. **EXTRACT GENOTYPE** (reference allele vs. alternative allele e.g. 0/1)
2. **EXTRACT DNA BASES** (reference allele vs. alternative allele)
3. **FIND THE CHILD** (count the no of incompatible variants)
4. **EXTRACT DE NOVO VARIANTS FROM THE CHILD** 
5. **FURTHER RESEARCH**
   * validation of de novo variants (**computational** e.g. calculate the probability for the mutation to be present in a family trio vs. **experimental** e.g. SANGER sequencing)
   * clinical impact of de novo varints (UCSC Genome Browser)
