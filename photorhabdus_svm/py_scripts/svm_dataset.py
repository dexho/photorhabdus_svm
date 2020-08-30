import os
import random
import pandas as pd 
from sklearn.datasets import load_iris
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

'''READS A FASTA FILE AND GENERATES A PANDAS DATAFRAME CONTAINING FEATURES
    POSITIVE/NEGATIVE TRAINTING/TESTING SEQUENCES ARE ALL CONTAINED IN RESPECTIVE FASTA FILES'''


gene_location = "/Users/desho/Desktop/photorhabdus_svm/genes"
os.chdir(gene_location)

iris = load_iris()
i = 1

def make_dataset(fasta):
    global i
    ls_pI = []
    ls_aac = []
    # ls_locus = []
    # ls_desc = []
    # assign whether it's effector or non-effector
    if 'positive' in fasta:
        target = 0
    elif 'negative' in fasta:
        target = 1

    for record in SeqIO.parse(fasta, "fasta"):
        if i:
            i = 1
            # test with isoelectric point
            analysed_seq = ProteinAnalysis(str(record.seq))
            # aac = analysed_seq.get_amino_acids_percent()
            # ls_aac += [aac]
            ls_pI += [analysed_seq.isoelectric_point()]
            # ls_locus += [record.name]
            # ls_desc += [record.description]
    

    # df = pd.DataFrame(ls_aac)
    df = pd.DataFrame()
    df['pI'] = ls_pI
    df['target'] = target
    # df['locus_id'] = ls_locus
    # df['desc'] = ls_desc
    df.to_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/' + fasta[:-6] + '.pkl')



for fasta in os.listdir(gene_location):
    if fasta[-5:] == 'fasta':
        make_dataset(fasta)

df_X_0 = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/negative_training_genes.pkl')
df_X_1 = pd.read_pickle('/Users/desho/Desktop/photorhabdus_svm/testpI/positive_training_genes.pkl')


frames = [df_X_0, df_X_1]
result = pd.concat(frames)
print(df_X_0)