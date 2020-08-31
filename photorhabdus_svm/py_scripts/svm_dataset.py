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

def make_dataset(gb):
    ls_aac = []

    # assign whether it's effector or non-effector
    if 'positive' in gb:
        target = 0
    elif 'negative' in gb:
        target = 1

    helix = {'H_A': 0, 'H_C': 0, 'H_D': 0, 'H_E': 0, 'H_F': 0, 'H_G': 0, 'H_H': 0, 'H_I': 0, 'H_K': 0, 'H_L': 0, 'H_M': 0, 'H_N': 0, 'H_P': 0, 'H_Q': 0, 'H_R': 0, 'H_S': 0, 'H_T': 0, 'H_V': 0, 'H_W': 0, 'H_Y': 0}
    strand = {'E_A': 0, 'E_C': 0, 'E_D': 0, 'E_E': 0, 'E_F': 0, 'E_G': 0, 'E_H': 0, 'E_I': 0, 'E_K': 0, 'E_L': 0, 'E_M': 0, 'E_N': 0, 'E_P': 0, 'E_Q': 0, 'E_R': 0, 'E_S': 0, 'E_T': 0, 'E_V': 0, 'E_W': 0, 'E_Y': 0}

    for record in SeqIO.parse(gb, "genbank"):

        SS = []
        pointer = 0
        sequence = str(record.seq)
        length = min(100, len(sequence))

        if True:
            for feature in record.features:
                for x in range(feature.location.start - pointer):
                    if len(SS) >= length:
                        break
                    SS += ['n']
                for x in range(feature.location.end - feature.location.start):
                    if len(SS) >= length:
                        break
                    SS += [feature.type]
                pointer = feature.location.end
            
            # in case a strand/helix ends before the end of the actual sequence
            for x in range(len(SS), length):
                SS += ['n']
            
            # build a list SS that is a list of what kind of structure this residue is in

            helix_c, strand_c = 0, 0
            for r in range(length):
                if SS[r] == 'helix':
                    helix_c += 1
                    helix['H_' + sequence[r]] += 1
                elif SS[r] == 'strand':
                    strand_c += 1
                    strand['E_' + sequence[r]] += 1
            
            # increment helix and strand dictionaries for the amino acid in strand and helix
            # find frequency of each residue in this structure

            for k, v in helix.items():
                try:
                    helix[k] = v/helix_c
                except:
                    helix[k] = 0
            for k, v in strand.items():
                try:
                    strand[k] = v/strand_c
                except:
                    strand[k] = 0

        # add AAC composition of first 100
        analysed_seq = ProteinAnalysis(str(record.seq)[:length])
        aac = analysed_seq.get_amino_acids_percent()

        # merge all dictionaries into aac
        aac.update(helix)
        aac.update(strand)
        ls_aac += [aac]

    df = pd.DataFrame(ls_aac)
    df['target'] = target
    df.to_pickle('/Users/desho/Desktop/photorhabdus_svm/testSS/' + gb[:-3] + '.pkl')



for gb in os.listdir(gene_location):
    if 'gb' in gb:
        print(gb)
        make_dataset(gb)

