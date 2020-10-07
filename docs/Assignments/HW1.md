# Homework 1

1. Write a script called `download_count.sh` which does the following.
   * Download the data file [https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core) from NCBI
   * Print out the count of the number of FASTA format sequences in this file - see [Wikipedia FASTA format](https://en.wikipedia.org/wiki/FASTA_format) - each record starts with a `>`

2. Write a script called `summary_exons.sh` which summarizes the total length of exons in the file [data/rice_random_exons.bed](https://raw.githubusercontent.com/biodataprog/GEN220_data/main/data/rice_random_exons.bed)
   * read in the file
   * use a loop structure to read each line
   * add up the values into a variable

3. Write a script called `strand_gene_count.sh` to calculate the number of genes that are on the + and - strand in the file.
  * [https://fungidb.org/common/downloads/release-48/ScerevisiaeS288c/gff/data/FungiDB-48_ScerevisiaeS288c.gff](https://fungidb.org/common/downloads/release-48/ScerevisiaeS288c/gff/data/FungiDB-48_ScerevisiaeS288c.gff)
