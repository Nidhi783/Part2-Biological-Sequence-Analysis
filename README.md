# Assessment 4 – Part 2 (R)
**Student ID:** 225044829  

## Purpose
This repository contains my solution for **Part 2: Examining Biological Sequence Diversity**.  
It compares genomic and proteomic sequence features between *Escherichia coli* and *Streptacidiphilus jiangxiensis* using R.  
The analysis includes coding sequence composition, nucleotide and amino-acid frequencies, codon usage bias, and k-mer profiling.

---

## Script
- `analysis_part2.R` — main script that runs all 10 steps of the assignment.  
- The script automatically downloads, extracts, and analyses CDS data for both organisms.  
- All figures and results are produced directly in the R console or Plots pane.

---

## Inputs
The script downloads the required CDS FASTA files from Ensembl Bacteria at runtime:
- `ecoli_cds.fa.gz` — *E. coli* (strain K-12 MG1655) coding sequences.  
- `sj_cds.fa.gz` — *S. jiangxiensis* coding sequences.  

After extraction:
- `ecoli_cds.fa`  
- `sj_cds.fa`

---

## Outputs
- **Console:**  
  - Number of coding sequences for both organisms.  
  - Total coding DNA length, mean, and median gene lengths.  
  - Top over- and under-represented 3-mer sequences and overlap comparison.  

- **Plots generated:**  
  - Boxplot of CDS length distributions.  
  - Barplots for nucleotide and amino-acid frequencies.  
  - Codon usage bias summary.  
  - K-mer frequency comparison plots for both organisms.

---

## Summary
*Streptacidiphilus jiangxiensis* shows higher GC content, longer coding sequences, and stronger hydrophobic amino-acid bias than *E. coli*, reflecting adaptation to acidic, nutrient-poor environments.  
*E. coli* displays balanced nucleotide composition and moderate codon bias, optimised for rapid growth and translational efficiency.  
Both species favour alanine- and leucine-rich motifs and avoid cysteine-, methionine-, and tryptophan-rich triplets.

---

## How to Run
From the project folder in RStudio or Terminal:

```bash
R -f analysis_part2.R
