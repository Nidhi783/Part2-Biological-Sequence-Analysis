####  Part 2: Examining Biological Sequence Diversity
####  Organisms: E. coli vs Streptacidiphilus jiangxiensis
####  Student ID: 225044829

#### Step 1. Load libraries ----
if (!require("seqinr")) install.packages("seqinr", dependencies=TRUE)
if (!require("R.utils")) install.packages("R.utils", dependencies=TRUE)

library(seqinr)
library(R.utils)

#### Step 2. Download and extract coding DNA sequences

### E. coli
url_ecoli <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(url_ecoli, destfile = "ecoli_cds.fa.gz", mode = "wb")
if (!file.exists("ecoli_cds.fa")) {
  gunzip("ecoli_cds.fa.gz", remove = FALSE)
}

### Streptacidiphilus jiangxiensis
url_sj <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_57_collection/streptacidiphilus_jiangxiensis_gca_900109465/cds/Streptacidiphilus_jiangxiensis_gca_900109465.IMG-taxon_2675903135_annotated_assembly.cds.all.fa.gz"
download.file(url_sj, destfile = "sj_cds.fa.gz", mode = "wb")
if (!file.exists("sj_cds.fa")) {
  gunzip("sj_cds.fa.gz", remove = FALSE)
}

list.files()

#### Step 3. Read FASTA files
ecoli_cds <- read.fasta("ecoli_cds.fa")
sj_cds    <- read.fasta("sj_cds.fa")

####Number of coding sequences
length(ecoli_cds)   # number of CDS in E. coli
length(sj_cds)      # number of CDS in S. jiangxiensis

######Step 5. Total coding DNA length
len_ecoli <- as.numeric(summary(ecoli_cds)[,1])
len_sj    <- as.numeric(summary(sj_cds)[,1])

sum(len_ecoli)   # total bp for E. coli
sum(len_sj)      # total bp for S. jiangxiensis

#####Step 6. Mean, median & boxplot of gene lengths
mean(len_ecoli); median(len_ecoli)
mean(len_sj);    median(len_sj)

boxplot(list(Ecoli=len_ecoli, Streptacidiphilus=len_sj),
        ylab="CDS length (bp)",
        main="CDS length comparison")

#####Step 7. Nucleotide frequencies
dna_ecoli <- unlist(ecoli_cds)
dna_sj    <- unlist(sj_cds)

freq_ecoli <- count(dna_ecoli, 1)
freq_sj    <- count(dna_sj, 1)

par(mfrow=c(1,2))
barplot(freq_ecoli, main="E. coli Nucleotide Frequency", xlab="Base", ylab="Count", col="skyblue")
barplot(freq_sj,    main="S. jiangxiensis Nucleotide Frequency", xlab="Base", ylab="Count", col="salmon")


####Step 8. Protein translation & amino acid frequency
prot_ecoli <- lapply(ecoli_cds, translate)
prot_sj    <- lapply(sj_cds, translate)

aa_ecoli <- unlist(prot_ecoli)
aa_sj    <- unlist(prot_sj)

aa_freq_ecoli <- count(aa_ecoli, 1)
aa_freq_sj    <- count(aa_sj, 1)

par(mfrow=c(1,2))
barplot(aa_freq_ecoli, main="E. coli AA Frequency", las=2)
barplot(aa_freq_sj,    main="S. jiangxiensis AA Frequency", las=2)

####Step 9. Codon usage bias
codon_ecoli <- uco(unlist(ecoli_cds), index="rscu", as.data.frame=TRUE)
codon_sj    <- uco(unlist(sj_cds),    index="rscu", as.data.frame=TRUE)

head(codon_ecoli)
head(codon_sj)


####Step 10. K-mer profiling

aa <- unique(aa_ecoli)
aa <- aa[aa != "*"]

####Protein sequences for SJ
prots_sj <- unlist(prot_sj)

kmer_sj <- count(prots_sj, wordsize=3, alphabet=aa, freq=TRUE)

####Top 10 over- and under-represented
sort(kmer_sj, decreasing=TRUE)[1:10]
sort(kmer_sj, decreasing=FALSE)[1:10]

####Protein sequences for E. coli
prots_ecoli <- unlist(prot_ecoli)

####Count 3-mers
kmer_ecoli <- count(prots_ecoli, wordsize=3, alphabet=aa, freq=TRUE)

####Top 10 over- and under-represented in E. coli
sort(kmer_ecoli, decreasing = TRUE)[1:10]
sort(kmer_ecoli, decreasing = FALSE)[1:10]

# Define amino acid alphabet (remove stop codon)
aa <- unique(aa_ecoli)
aa <- aa[aa != "*"]

# Flatten protein sequences for each organism
prots_ecoli <- unlist(prot_ecoli)
prots_sj    <- unlist(prot_sj)

# ----- Count 3-mer frequencies -----
kmer_ecoli <- count(prots_ecoli, wordsize = 3, alphabet = aa, freq = TRUE)
kmer_sj    <- count(prots_sj,    wordsize = 3, alphabet = aa, freq = TRUE)

# ----- Identify over- and under-represented sequences -----
top_over_ecoli  <- sort(kmer_ecoli, decreasing = TRUE)[1:10]
top_under_ecoli <- sort(kmer_ecoli, decreasing = FALSE)[1:10]
top_over_sj     <- sort(kmer_sj,    decreasing = TRUE)[1:10]
top_under_sj    <- sort(kmer_sj,    decreasing = FALSE)[1:10]

# ----- Barplots for visual comparison -----
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))  # Reduced margins to prevent "figure margins too large"

barplot(top_over_sj, las = 2, col = "tomato",
        main = "Top 10 Over-represented 3-mers (S. jiangxiensis)")
barplot(top_over_ecoli, las = 2, col = "skyblue",
        main = "Top 10 Over-represented 3-mers (E. coli)")
barplot(top_under_sj, las = 2, col = "tomato",
        main = "Top 10 Under-represented 3-mers (S. jiangxiensis)")
barplot(top_under_ecoli, las = 2, col = "skyblue",
        main = "Top 10 Under-represented 3-mers (E. coli)")

# ----- Compare overlap between organisms -----
intersect(names(top_over_sj), names(top_over_ecoli))
intersect(names(top_under_sj), names(top_under_ecoli))