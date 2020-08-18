import sklearn
from sklearn import datasets
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

for record in SeqIO.parse("/Users/desho/Desktop/photorhabdus_training_positives/x.fasta", "fasta"):
    print(str(record.seq))
    analysed_seq = ProteinAnalysis(str(record.seq))
    print(analysed_seq.isoelectric_point)




#cancer = datasets.load_breast_cancer()






















