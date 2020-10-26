# Parse and Converting workshop

## BLAST report
- [Ecoli proteins](data/E_coli_K12.pep.gz)
- [S_enterica proteins](data/S_enterica.pep.gz)
- [Ecoli-vs-Senterica.BLASTP.tab.gz](data/Ecoli-vs-Senterica.BLASTP.tab.gz)

Write script to read in the BLAST report.
- Calculate for each alignment what the % aligned of the Ecoli query protein is?
- Calculate the % of proteome that was aligned (out of the total number of Ecoli proteins)

Parsing BLAST
[https://github.com/biodataprog/GEN220_2019_examples/blob/master/Bioinformatics_2/parse_blast.py](https://github.com/biodataprog/GEN220_2019_examples/blob/master/Bioinformatics_2/parse_blast.py)

## Orthofinder parsing

Here is a data file - [OrthoFinder result](data/Orthogroups.csv)

- Let's write a script which will summarize the data by Species.

The data look like
 ```
 ORTHOGROUP	GENENAME_SP1, GENENAME2_SP1	GENENAME_SP2
 ```

