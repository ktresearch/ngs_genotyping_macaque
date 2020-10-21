# -*- coding: utf-8 -*-
"""
Created on Sun May 24 2020

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

input = sys.argv
min_length = sys.argv[1]
inputfile = sys.argv[2]
inputfastafile = sys.argv[3]
multihitfile = sys.argv[4]
outputname = sys.argv[5]
species = sys.argv[6]
reverse_undetermined_file = f'reference/{species}_reverse_undetermined_alleles.txt'

outputfile = f'count_uniquehit_{str(min_length)}_{outputname}.txt'
subprocess.run([f'mkdir coverage_uniquehit_{min_length}_{outputname}'], shell=True)
subprocess.run([f'mkdir coverage_uniquehit_{min_length}_{outputname}_plot'], shell=True)

#count a number of reads in a fasta file
fasta_count = 0
with open(inputfastafile) as f:
    for _ in f:
        fasta_count += 1
totalreads = fasta_count / 2

rev_undet = {}
with open(reverse_undetermined_file) as f:
    for line in f.readlines():
        rev_undet.setdefault(line.strip(), 1)

#for tqdm
line_count = 0
with open(inputfile) as f:
    for line in f:
        line_count += 1
pbar = tqdm(total=line_count)


#obtain a result of multihit
uniformity = []
with open(multihitfile) as f:
    for line in f.readlines()[1:]:
        cols = line.strip().split('\t')
        if re.search('DeNovo', cols[0]): #uniformity >= 0.9 in DeNovo detection
            if float(cols[2]) >= 0.9:
                uniformity.append(cols[0])
        else:
            if float(cols[2]) >= 1:
                uniformity.append(cols[0])

count = {}
count_tmp = {}
depth = {}
with open(inputfile) as f:
    for line in f:
        pbar.update(1) #tqdm
        cols = line.strip().split('\t')
        header_info = cols[0].split('_')
        read_length = int(header_info[1])
        if read_length < int(min_length):
            continue

        #skip if an uniformity of multihit is not 1
        if not cols[1] in uniformity:
            continue

        undet_check = 0
        if not re.search('DeNovo', cols[1]):
            allele_tmp = re.sub('//.+', '', cols[1])
            match = re.search('(.+)__Length(\d+)', allele_tmp)
            allele_length = match.group(2)
            for allele in match.group(1).split('<>'):
                if allele in rev_undet:
                    undet_check = 1

        if float(cols[2]) >= 100: #identity100
            if int(cols[3]) == read_length: #match entire read
                count_tmp.setdefault(cols[0], [])
                count_tmp[cols[0]].append(line.strip())
            else: #if not match entire read, consider reverse primer undetermined
                if undet_check == 1: #reverse undetermined
                    if int(cols[8]) > int(cols[9]):
                        dbst = int(cols[9])
                        dben = int(cols[8])
                    else:
                        dbst = int(cols[8])
                        dben = int(cols[9])
                    if dben == int(allele_length): #if match end of db
                        if int(cols[6]) == 1 or int(cols[7]) == read_length:
                            count_tmp.setdefault(cols[0], [])
                            count_tmp[cols[0]].append(line.strip())

#count uniquehit
for header in count_tmp.keys():
    if len(count_tmp[header]) == 1:
        cols = count_tmp[header][0].split('\t')
        count.setdefault(cols[1], 0)
        count[cols[1]] += 1
        if int(cols[8]) < int(cols[9]):
            dbst = int(cols[8])
            dben = int(cols[9])
        else:
            dbst = int(cols[9])
            dben = int(cols[8])
        for i in range(dbst,dben + 1):
            depth.setdefault(cols[1], {})
            depth[cols[1]].setdefault(i, 0)
            depth[cols[1]][i] += 1

#calculate uniformity
with open(outputfile,'w') as out:
    print('MHC_type', 'normalized_count', 'uniformity', 'TotalReads_' + str(totalreads), sep='\t', file=out)
    for allele in count.keys():
        if re.search('//', allele):
            allele_tmp = re.sub('//.+', '', allele)
        else:
            allele_tmp = allele
        if re.search('DeNovo', allele_tmp):
            match = re.search('_Length(\d+)', allele_tmp)
        else:
            match = re.search('__Length(\d+)', allele_tmp)
        allele_length = match.group(1)
        depth_sum = 0
        for i in range(1, int(allele_length) + 1):
            if i in depth[allele]:
                depth_sum += depth[allele][i]
        depth_ave = depth_sum / int(allele_length)
        over_ave_nt = 0
        for i in range(1, int(allele_length) + 1):
            if i in depth[allele]:
                if depth[allele][i] >= depth_ave / 10:
                    over_ave_nt += 1
        uniformity = over_ave_nt / int(allele_length)

        if re.search('<>', allele_tmp):
            filename_tmp = re.sub('<>.+', '', allele_tmp)
            if species == 'Mamu':
                match_name = re.search('Mamu.+', filename_tmp)
            elif species == 'Mafa':
                match_name = re.search('Mafa.+', filename_tmp)
            filename = match_name.group()
        else:
            if re.search('DeNovo', allele_tmp):
                match_name = re.search('(.+)_Length\d+', allele_tmp)
                filename = match_name.group(1)
            else:
                if species == 'Mamu':
                    match_name = re.search('(Mamu.+)__Length\d+', allele_tmp)
                elif species == 'Mafa':
                    match_name = re.search('(Mafa.+)__Length\d+', allele_tmp)
                filename = match_name.group(1)
        filename = re.sub('[*:]', '_', filename)

        plot_x = []
        plot_y = []
        outputfile_depth = f'coverage_uniquehit_{min_length}_{outputname}/{filename}.txt'
        with open(outputfile_depth, 'w') as depth_out:
            for i in range(1, int(allele_length) + 1):
                plot_x.append(i)
                if i in depth[allele]:
                    print(i, depth[allele][i], sep='\t', file=depth_out)
                    plot_y.append(depth[allele][i])
                else:
                    print(i, 0, sep='\t', file=depth_out)
                    plot_y.append(0)

        plt.figure()
        plt.stackplot(plot_x, plot_y, colors='gray')
        plt.title(filename, fontsize=10)
        plt.xlabel('position')
        plt.ylabel('count')
        plt.xlim(0, max(plot_x))
        plt.ylim(0,)
        plt.savefig(f'coverage_uniquehit_{min_length}_{outputname}_plot/{filename}.pdf')
        plt.close()

        normalized_count = (count[allele] / totalreads) * 100000
        print(allele, normalized_count, uniformity, sep='\t', file=out)
