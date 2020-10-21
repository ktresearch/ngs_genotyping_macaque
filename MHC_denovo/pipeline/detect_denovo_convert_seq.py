# -*- coding: utf-8 -*-
"""
Created on Wed July 22 2020

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


inputfile = 'detect_denovo_out3.txt'
outputfile1 = 'detect_denovo_out4.txt'
outputfile2 = 'denovo_candidate1.fasta'

referencefile = sys.argv[1]
blastfile = sys.argv[2]
fastafile = sys.argv[3]


ref = {}
id = {}
with open(referencefile) as f:
    for line in f:
        if re.match('>', line):
            id_temp = re.sub('^>', '', line.strip())
        else:
            ref.setdefault(id_temp, line.strip())


fasta = {}
with open(outputfile1, 'w') as out:
    with open(inputfile) as f:
        for line in f:
            cols = line.strip().split('\t')
            cutoff = int(cols[1]) * 0.1
            seq = list(ref[cols[0]]) #convert list
            mismatch = 0
            for var in cols[2:]:
                info = var.split(',')
                if int(info[1]) >= cutoff:
                    m = re.search('(\d+)([ATGC])', info[0])
                    seq[int(m.group(1))] = m.group(2)
                    mismatch += 1
            seq = ''.join(seq) #convert str
            if mismatch == 0:
                continue
            else:
                fasta.setdefault(seq, [])
                fasta[seq].append(cols[0])
    for seq_temp in fasta.keys():
        print('<>'.join(fasta[seq_temp]), seq_temp, sep='\t', file=out)


fasta = {}
with open(outputfile1) as f:
    for line in f:
        cols = line.strip().split('\t')
        fasta.setdefault(cols[0], cols[1])


with open(outputfile1) as f:
    denovo_temp = []
    for line in f:
        cols = line.strip().split('\t')
        denovo_temp.append(cols[1])

    denovo_temp = sorted(list(set(denovo_temp)), key=len)
    denovo = []
    for i in range(len(denovo_temp)):
        discard = 0
        for j in range(i + 1, len(denovo_temp)):
            if len(denovo_temp[i]) == len(denovo_temp[j]):
                continue
            elif denovo_temp[i] in denovo_temp[j]:
                discard = 1
            else:
                pass
        if discard == 0:
            denovo.append(denovo_temp[i])

with open(outputfile2, 'w') as out:
    denovo_id = 0
    for seq in denovo:
        denovo_id += 1
        print('>DeNovo_candidate_' + str(denovo_id) + '_Length' + str(len(seq)), seq, sep='\n', file=out)
