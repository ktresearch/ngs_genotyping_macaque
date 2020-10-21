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


min_length = int(sys.argv[1])
inputfile = sys.argv[2]
fastafile = sys.argv[3]
referencefile = sys.argv[4]
outputfile1 = 'detect_denovo_out1.txt'
outputfile2 = 'detect_denovo_out2.txt'
outputfile3 = 'detect_denovo_out3.txt'

seq = {}
with open(fastafile) as f:
    header = ''
    for line in f:
        if re.match('>', line):
            header = re.sub('^>', '', line.strip())
        else:
            seq.setdefault(header, line.strip())


blastmatch = {}
with open(outputfile1, 'w') as out:
    with open(inputfile) as f:
        for line in f:
            cols = line.strip().split('\t')
            length = int(cols[0].split('_')[-1])
            if length < min_length:
                continue
            if not cols[0] in blastmatch:
                blastmatch.setdefault(cols[0], float(cols[2]))
            else:
                if float(cols[2]) >= blastmatch[cols[0]]:
                    blastmatch[cols[0]] = float(cols[2])
            if blastmatch[cols[0]] == 100:
                continue
            else:
                print(line.strip(), seq[cols[0]], sep='\t', file=out)


ref = {}
with open(referencefile) as f:
    header = ''
    for line in f:
        if re.match('>', line):
            header = re.sub('^>', '', line.strip())
        else:
            ref.setdefault(header, line.strip())


with open(outputfile2, 'w') as out:
    with open(outputfile1) as f:
        for line in f:
            cols = line.strip().split()
            length = int(cols[0].split('_')[-1])
            if int(cols[8]) > int(cols[9]):
                query = cols[12].translate(str.maketrans('ATGC', 'TACG'))[::-1]
                refst = int(cols[9]) - 1
            else:
                query = cols[12]
                refst = int(cols[8]) - 1
            if not (length == int(cols[7]) and int(cols[6]) == 1):
                continue
            if int(cols[5]) != 0:
                continue
            db = ref[cols[1]][refst:refst + length]
            mismatchpos = []
            for i in range(length):
                if query[i] == db[i]:
                    continue
                else:
                    mismatchpos.append(str(refst + i) + query[i])
            print(line.strip(), ','.join(mismatchpos), sep='\t', file=out)


count = {}
pos = {}
with open(outputfile2) as f:
    for line in f:
        cols = line.strip().split('\t')
        count.setdefault(cols[1], 0)
        count[cols[1]] += 1
        pos.setdefault(cols[1], '')
        pos[cols[1]] += cols[13]


with open(outputfile3, 'w') as out:
    for allele in count.keys():
        pos_list = pos[allele].split(',')
        pos_count = {}
        for pos_temp in pos_list:
            pos_count.setdefault(pos_temp, 0)
            pos_count[pos_temp] += 1
        pos_out = []
        for pos_temp, pos_count_temp in sorted(pos_count.items(), key=lambda x:x[1], reverse=True):
            pos_out.append(pos_temp + ',' + str(pos_count_temp))
        print(allele, count[allele], '\t'.join(pos_out), sep='\t', file=out)
