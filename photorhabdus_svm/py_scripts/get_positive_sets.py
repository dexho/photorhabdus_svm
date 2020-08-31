import os
import csv
import random
import pandas as pd 
from Bio import SeqIO

'''MAKE A 2D LIST CONTAINING SEPARATED LOCI FOR GENES IN THE ORTHOGROUP. EXPORT GB FILE FROM ASSEMBLY.'''

orthofinder = '/Users/desho/Desktop/photorhabdus_svm/positives_intermediate/orthofinder_output.csv'
photorhabdus = '/Library/Photorhabdus/unzip'
os.chdir(photorhabdus)

df = pd.read_csv(orthofinder)

dfb = pd.DataFrame()
i = 0

effector = [[] for i in range(40)]

for column in df:
    df.rename(columns = {column:column[:-2]}, inplace = True)

def collect_genes(gb):
    global effector
    seq_record = SeqIO.parse(gb, "genbank")
    for record in seq_record:
        if record.name in df:
            print(record.name)

            # build 2D list of orthogroups and genes
            orthogroups = list(df[record.name])
            for x in range(40):
                if type(orthogroups[x]) == str:
                    effector[x] += orthogroups[x].split(',')

            # save contigs that contain the relevant genes
            with open(record.name + '.gb', 'w') as output_handle:
                SeqIO.write(record, output_handle, 'genbank')

for gb in os.listdir(photorhabdus):
    if 'gbff' in gb:
        collect_genes(gb)

df = pd.DataFrame(effector)
df.to_pickle('pickled_effector')