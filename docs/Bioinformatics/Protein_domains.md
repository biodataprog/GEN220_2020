# Annotating Proteins

# Finding homologs

For Protein to Protein searches
BLASTP, phmmer ([HMMER](https://hmmer-web-docs.readthedocs.io/en/latest/algorithms.html#hmmer-algorithms)), FASTA

```bash
module load fasta
fasta36 query database > results.FASTA
fasta36 -m 8c -E 1e-3 query database > results.FASTA.tab
```

# To Find Domains

See Overview lecture [Domains lecture](Domains_lecture.pdf)

Searching with [HMMer](http://hmmer.org/) against [Pfam](http://pfam.xfam.org)

See the HMMER [tutorial](http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf)

Searching with [Interpro](https://www.ebi.ac.uk/interpro/search/)

## Searchin Intepro on HPCC

Note this can be slow.

```bash
#SBATCH -p batch -N 1 -n 8
module load iprscan
CPU=4
interproscan.sh  --goterms --pathways -f tsv -i PROTEINFILE.fa --cpu $CPU > SEARCH.log
```

The results will contain information like

Gene Ontology [http://geneontology.org/](http://geneontology.org/)

# Running Analyses on Biocluster

```bash
module load hmmer
module load db-pfam
hmmscan --domtbl domtbl_results.out $PFAM_DB/Pfam-A.hmm proteins.fa > proteins.hmmscan
hmmsearch --domtbl domtbl_results.out $HMM protein-db.fa > protein.hmmsearch
```


Pfam2GO - [http://current.geneontology.org/ontology/external2go/pfam2go](http://current.geneontology.org/ontology/external2go/pfam2go)

# Workshop

1. Searching for Pfam domains in sets of proteins - [https://github.com/biodataprog/GEN220_2019_examples/tree/master/Bioinformatics_4](Bioinformatics_4)
See [search_SOD1.sh](https://github.com/biodataprog/GEN220_2019_examples/blob/master/Bioinformatics_4/search_SOD1.sh)
2. Parsing report files
