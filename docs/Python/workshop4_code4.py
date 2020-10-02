import itertools, sys, re

# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'


def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

#filename="/bigdata/gen220/shared/data_files/S_cerevisiae_ORFs.fasta"
filename="S_cerevisiae_ORFs.fasta"
with open(filename,"r") as f:
   seqs = dict(aspairs(f))
   first_codon = {}
   last_codon  = {}
   sequence_count = 0
   strand_count = {}

   for seqname in seqs:
       sequence_count += 1
       
       firstcodon = seqs[seqname][0:3]
       lastcodon = seqs[seqname][-3:]
       if firstcodon in first_codon:
           first_codon[firstcodon] +=1
       else:
           first_codon[firstcodon] = 1

       if lastcodon in last_codon:
           last_codon[lastcodon] +=1
       else:
           last_codon[lastcodon] = 1

#       print(firstcodon,lastcodon)
       last_char = seqname[-1]
#       print(last_char, " in ", seqname)
       if last_char in strand_count:
           strand_count[last_char] += 1
       else:
           strand_count[last_char] = 1
print("1.")
print("There are %d sequences in the file"%(sequence_count))

print("2.")
print("The distribution of first codons is:")

for codon in first_codon:
    print("%s => %d (%.1f%%)" % (codon, first_codon[codon],
                                 100.0 * first_codon[codon] / sequence_count))

    
print("The distribution of last codons is:")
for codon in last_codon:
    print("%s => %d (%.1f%%)" % (codon, last_codon[codon],
                                 100.0 * last_codon[codon] / sequence_count))

print("There are %s genes on the Watson (+) strand"%(strand_count['W']))
print("          %s genes on the Crick  (-) strand"%(strand_count['C']))

print("3.")
print("There are %s genes on the Watson (+) strand"%(strand_count['W']))
print("          %s genes on the Crick  (-) strand"%(strand_count['C']))
