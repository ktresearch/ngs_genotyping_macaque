#!/bin/sh

read -p "enter Fastqfile:" fastqfile
read -p "enter Mamu or Mafa:" species
read -p "enter runid (any characters):" runid
read -p "enter barcodeid (any characters):" barcodeid
read -p "enter samplename (any characters):" sampleid

length=150
qv=20
percent=90

if [ $species = "Mamu" ]; then
  blast_reference_fasta=./reference/190312_Mamu_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
else
  blast_reference_fasta=./reference/190304_Mafa_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
fi

cutadapt -b GTGGGCTACGTGGACGAC $fastqfile > $runid"_"$barcodeid"_removeFp.fastq"
cutadapt -b GTCGTCCACGTAGCCCAC $runid"_"$barcodeid"_removeFp.fastq" > $runid"_"$barcodeid"_removeFpm.fastq"

cutadapt -b CTTCTACCCTGCGGAGATCA $runid"_"$barcodeid"_removeFpm.fastq" > $runid"_"$barcodeid"_removeFRp.fastq"
cutadapt -b TGATCTCCGCAGGGTAGAAG $runid"_"$barcodeid"_removeFRp.fastq" > $runid"_"$barcodeid"_removeFR.fastq"

perl pipeline/fastq_remove_readlength.pl $length $runid"_"$barcodeid"_removeFR.fastq" $runid"_"$barcodeid"_"$length"bp.fastq"

perl pipeline/fastq_remove_quality.pl 20 90 $runid"_"$barcodeid"_"$length"bp.fastq" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fastq"

perl pipeline/fastq2fasta.pl $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fastq" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta"

blastn -db $blast_reference_fasta -query $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" -evalue 1e-10 -perc_identity 98 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98.txt"

python3 pipeline/count_multihit_v5.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" $species

python3 pipeline/count_uniquehit_v5.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" "count_multihit_"$length".txt" $species

python3 pipeline/sort.py $sampleid
