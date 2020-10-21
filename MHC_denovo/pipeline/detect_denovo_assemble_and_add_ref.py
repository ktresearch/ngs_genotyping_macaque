# -*- coding: utf-8 -*-
"""
Created on Mon Aug 3 2020

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

reffile = sys.argv[1]
fastafile = sys.argv[2]
blastresultfile = sys.argv[3]
outputfile = sys.argv[4]
switch = sys.argv[5]


fasta = {}
header = ''
with open(fastafile) as f:
    for line in f:
        if re.match('>', line):
            header = re.sub('^>', '', line.strip())
        else:
            fasta.setdefault(header, line.strip())

remove = []
seq_list = []
combination = []
with open(blastresultfile) as f:
    for line in f:
        cols = line.strip().split('\t')
        if cols[0] == cols[1]:
            continue
        if float(cols[2]) != 100:
            continue

        combination_check = 0
        for combi in combination:
            if cols[0] in combi and cols[1] in combi:
                combination_check = 1
        if combination_check == 1:
            continue

        match = re.search('_Length(\d+)', cols[0])
        qulength = int(match.group(1))
        match = re.search('_Length(\d+)', cols[1])
        dblength = int(match.group(1))
        if int(cols[6]) == 1 and int(cols[9]) == dblength:
            seq = ''
            for i in range(int(cols[8]) - 1):
                seq += fasta[cols[1]][i]
            seq += fasta[cols[0]]
            seq_list.append(seq)
            remove.append(cols[0])
            remove.append(cols[1])
            combination.append([cols[0], cols[1]])
        elif int(cols[8]) == 1 and int(cols[7]) == qulength:
            seq = ''
            for i in range(int(cols[6]) - 1):
                seq += fasta[cols[0]][i]
            seq += fasta[cols[1]]
            seq_list.append(seq)
            remove.append(cols[0])
            remove.append(cols[1])
            combination.append([cols[0], cols[1]])
        else:
            continue
remove = list(set(remove))

seq_list = sorted(list(set(seq_list)), key=len)
remove_iter = []
for i in range(len(seq_list)):
    check = 0
    for j in range(i + 1, len(seq_list)):
        if seq_list[i] == seq_list[j]:
            check = 1
        elif seq_list[i] in seq_list[j]:
            check = 1
        else:
            continue
    if check == 1:
        remove_iter.append(i)

seq_list_res = []
for i in range(len(seq_list)):
    if i in remove_iter:
        continue
    else:
        seq_list_res.append(seq_list[i])


with open(outputfile, 'w') as out:
    if switch == 'on':
        with open(reffile) as f:
            for line in f:
                print(line.strip(), file=out)

    id = 0
    with open(fastafile) as f:
        header_out = ''
        skip = 0
        for line in f:
            if re.match('>', line):
                header = re.sub('^>', '', line.strip())
                if header in remove:
                    skip = 1
                else:
                    id += 1
                    header_out = f'>DeNovo_candidate2_{str(id)}_Length'
            else:
                if skip == 1:
                    skip = 0
                else:
                    header_out += str(len(line.strip()))
                    print(header_out, line.strip(), sep='\n', file=out)
                    header_out = ''

    for seq_out in seq_list_res:
        id += 1
        header_out = f'>DeNovo_candidate2_{str(id)}_Length{str(len(seq_out))}'
        print(header_out, seq_out, sep='\n', file=out)
