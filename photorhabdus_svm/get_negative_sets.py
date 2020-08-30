import os
import csv
import random
import pandas as pd 
from Bio import SeqIO

negatives_location = "/Users/desho/Desktop/photorhabdus_svm/negatives/fasta"
training_filename = "/Users/desho/Desktop/photorhabdus_svm/negatives/negative_training_set.fasta"
testing_filename = "/Users/desho/Desktop/photorhabdus_svm/negatives/negative_testing_set.fasta"

os.chdir(negatives_location)

all_negatives = []

# collect negative genes from all three genomes into list ALL_NEGATIVES # 
def collect_genes(fasta):
    global all_negatives

    seq_record = SeqIO.parse(fasta, "fasta")
    for record in seq_record:
        all_negatives += [record]

# split all_negatives into training and testing sets, each with 20*20=400 genes
def split(all_negatives):
    split = random.sample(all_negatives, 800)
    training = split[:400]
    testing = split[400:]

    with open(training_filename, "a") as handle:
        SeqIO.write(training, handle, "fasta")

    with open(testing_filename, "a") as handle:
        SeqIO.write(testing, handle, "fasta")

    #list of SeqRecord objects
    #add locus tag to SeqRecord description
    #SeqIO.write(record, handle, "fasta")


for fasta in os.listdir(negatives_location):
    collect_genes(fasta)
    print(fasta)

split(all_negatives)

