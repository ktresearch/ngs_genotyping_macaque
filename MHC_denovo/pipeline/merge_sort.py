# -*- coding: utf-8 -*-
"""
Created on Tue July 21 2020

@author: TANIMOTO Kousuke
"""

import re
#import numpy as np
#from sklearn.cluster import KMeans
import sys
import glob
from tqdm import tqdm
import subprocess
from matplotlib import pyplot as plt


samplename = sys.argv[1]

alleles = {}
with open('count_uniquehit_150.txt') as f:
    header, *lines = f.readlines()
    for line in lines:
        cols = line.strip().split('\t')
        alleles.setdefault(cols[0], line.strip())

with open('count_uniquehit_150_denovo2.txt') as f:
    for line in f.readlines()[1:]:
        cols = line.strip().split('\t')
        if re.match('DeNovo_', cols[0]):
            alleles.setdefault(cols[0], line.strip())

with open(f'{samplename}_result.txt', 'w') as out:
    print(header.strip(), file=out)
    for allele, info in sorted(alleles.items(), key=lambda x:float(x[1].split('\t')[1]), reverse=True):
        print(info, file=out)
