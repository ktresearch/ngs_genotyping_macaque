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


candidatefile = 'denovo_candidate1.fasta'
blastfile = 'denovo_candidate1_blast_results.txt'
outputfile = 'denovo_candidate2.fasta'

remove = []
with open(blastfile) as f:
    for line in f:
        cols = line.strip().split('\t')
        match = re.search('Length(\d+)', cols[0])
        readlength = int(match.group(1))
        if float(cols[2]) == 100 and readlength == int(cols[3]):
            remove.append(cols[0])

with open(outputfile, 'w') as out:
    with open(candidatefile) as f:
        check = 0
        for line in f:
            if re.match('>', line):
                if re.sub('^>', '', line.strip()) in remove:
                    check = 1
                else:
                    print(line.strip(), file=out)
            else:
                if check == 1:
                    check = 0
                    continue
                else:
                    print(line.strip(), file=out)
