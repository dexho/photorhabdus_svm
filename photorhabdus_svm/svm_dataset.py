import csv
import random
import pandas as pd 
from Bio import SeqIO

'''READS A FASTA FILE AND GENERATES A PANDAS DATAFRAME CONTAINING FEATURES
    POSITIVE/NEGATIVE TRAINTING/TESTING SEQUENCES ARE ALL CONTAINED IN RESPECTIVE FASTA FILES'''


x = list(range(1, 31))
random.shuffle(x)
print(x)