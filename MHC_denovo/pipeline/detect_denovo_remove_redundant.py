# -*- coding: utf-8 -*-
"""
Created on Fri July 31 2020

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


inputfile = 'count_multihit_150_denovo1.txt'
fastafile = 'denovo_candidate2.fasta'
outputfile = 'detect_denovo_out5.fasta'


fasta = {}
with open(fastafile) as f:
    header = ''
    for line in f:
        if re.match('>', line):
            header = re.sub('^>', '', line.strip())
            header = re.sub('_Length\d+', '', header)
        else:
            fasta.setdefault(header, line.strip())

seq_list = []
with open(outputfile, 'w') as out:
    with open(inputfile) as f:
        lines = f.readlines()
        for line in lines[1:]:
            cols = line.strip().split('\t')
            if float(cols[2]) < 0.5 or float(cols[1]) < 500:
                continue
            name = re.sub('_Length\d+', '', cols[0])
            with open(f'coverage_multihit_150_denovo1/{name}.txt') as fh_cov:
                #calculate ave coverage
                sum = 0
                lines_cov = fh_cov.readlines()
                for line_cov in lines_cov:
                    cols_cov = line_cov.strip().split('\t')
                    sum += int(cols_cov[1])
                ave = sum / len(lines_cov)

            with open(f'coverage_multihit_150_denovo1/{name}.txt') as fh_cov:
                seq = ''
                for line_cov in fh_cov:
                    cols_cov = line_cov.strip().split('\t')
                    if int(cols_cov[1]) > ave * 0.1:
                        seq += fasta[name][int(cols_cov[0]) - 1]
                    else:
                        if seq:
                            seq_list.append(seq)
                            seq = ''
                        else:
                            continue
                if seq:
                    seq_list.append(seq)

    candidateid = 0
    seq_list = sorted(list(set(seq_list)), key=len)
    for i in range(len(seq_list)):
        discard = 0
        for j in range(i + 1, len(seq_list)):
            if seq_list[i] == seq_list[j]:
                discard = 1
            elif seq_list[i] in seq_list[j]:
                discard = 1
            else:
                continue
        if discard == 0:
            candidateid += 1
            print('>DeNovo_candidate2_' + str(candidateid) + '_Length' + str(len(seq_list[i])), seq_list[i], sep='\n', file=out)
