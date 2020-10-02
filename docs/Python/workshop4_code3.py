import csv

datafile = 'species.csv'

genera = {}
with open(datafile,"r") as fh:
    reader = csv.reader(fh,delimiter=",")
    for row in reader:
        if row[1] in genera:
            genera[row[1]] += 1
        else:
            genera[row[1]] = 1
for genus in genera:
    print(genus,genera[genus])
