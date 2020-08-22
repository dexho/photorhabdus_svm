import sklearn
import pandas as pd
from sklearn import datasets
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

#reduce to N-terminal 1-30 region
# f = open("nterminal_30.fasta", "a")
# for record in SeqIO.parse("/Users/desho/Desktop/photorhabdus_svm/training_positives/x.fasta", "fasta"):
#     #print(str(record.seq))
#     record.seq = record.seq[0:30]
#     with open("nterminal_30.fasta", "a") as handle:
#         SeqIO.write(record, handle, "fasta")
#     print(record)
#     #analysed_seq = ProteinAnalysis(str(record.seq))
#     #print(analysed_seq.isoelectric_point)





#cancer = datasets.load_breast_cancer()






















