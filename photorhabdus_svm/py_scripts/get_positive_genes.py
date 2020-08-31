import os
import csv
import random
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''ADD ALL GENES FROM 108 GENBANK FILES INTO LIST. 
    PULL THE ONES CORRESPONDING TO EACH ORTHOGROUP AND EXPORT IN ONE FASTA. 
    ALL NEGATIVES GO TO ANOTHER FASTA.'''

photorhabdus = '/Library/Photorhabdus/genbanks'
os.chdir(photorhabdus)
df = pd.read_pickle('/Library/Photorhabdus/pickled_effector')

# print(df.iloc[0][1])
# print(df.iloc[0][2])
# print(list(df.iloc[0]))

all_genes = {}

def collect_genes(gb):
    global all_genes
    seq_record = SeqIO.parse(gb, "genbank")

    # add all sequences to all_genes
    for record in seq_record:
        for feature in record.features:
            if feature.type == 'CDS':
                try:
                    translation = feature.qualifiers['translation']
                except:
                    translation = None
                locus = feature.qualifiers['locus_tag']
                if translation and locus:
                    try:
                        desc = feature.qualifiers['product'][0]
                    except:
                        desc = ''
                    new_seqrecord = SeqRecord(
                        Seq(translation[0]),
                        id = locus[0],
                        name = feature.qualifiers['protein_id'][0],
                        description = desc,
                    )
                    all_genes[locus[0]] = new_seqrecord

def extract_orthogroups():
    global all_genes

    for i in range(40):
        orthogroup = []
        print(i, len(df.iloc[i]))
        for locus in df.iloc[i]:
            try:
                locus = locus.strip()
                orthogroup += [all_genes.pop(locus)]
            except:
                continue

        # with open('/library/Photorhabdus/orthogroup_' + str(i) + '.fasta', "a") as handle:
        #     SeqIO.write(orthogroup, handle, "fasta")

for gb in os.listdir(photorhabdus):
    if 'gb' in gb:
        print(gb)
        collect_genes(gb)

extract_orthogroups()
with open('/library/Photorhabdus/negative.fasta', 'a') as handle:
    negatives = random.sample(list(all_genes.values()), 605)
    SeqIO.write(negatives, handle, "fasta")
