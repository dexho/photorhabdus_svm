import pandas as pd
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''TAKES A GB GENOME FILE AND A CSV LIST OF LOCUS IDS TO AVOID. OUTPUT A FASTA FILE WITH ALL OTHER CDS'''

gb = 'FM162591'
filename = gb + "_negative.fasta"
genbank_location = "/Users/desho/Desktop/photorhabdus_svm/negatives/genbank/" + gb + ".gb"
avoid_positive = "/Users/desho/Desktop/photorhabdus_svm/negatives/avoid/" + gb + ".csv"

def get_negatives(gb, avoid_locus):
    all_genes = []
    seq_record = SeqIO.parse(genbank_location, "genbank")
    for record in seq_record:
        for feature in record.features:
            if feature.type == 'CDS':
                try:
                    translation = feature.qualifiers['translation']
                except:
                    translation = None
                locus = feature.qualifiers['locus_tag']
                if translation and locus:
                    if locus[0] in avoid_locus:
                        print(locus[0])
                    else:
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
                        all_genes.append(new_seqrecord)
                    

    
    with open(filename, "a") as handle:
        SeqIO.write(all_genes, handle, "fasta")
    #     #list of SeqRecord objects
    #     #add locus tag to SeqRecord description
    #     #SeqIO.write(record, handle, "fasta")

avoid_locus = []

with open(avoid_positive) as csvfile:
    avoid = csv.reader(csvfile)
    for row in avoid:
        avoid_locus += row

get_negatives(genbank_location, avoid_locus)