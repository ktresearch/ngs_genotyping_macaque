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
  blast_reference_fasta=reference/190312_Mamu_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
  blast_reference_fasta_filename=190312_Mamu_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
else
  blast_reference_fasta=reference/190304_Mafa_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
  blast_reference_fasta_filename=190304_Mafa_classI_Amplifiable_ver2_Realigned_Include_Rev1.fasta
fi

cutadapt -b GTGGGCTACGTGGACGAC $fastqfile > $runid"_"$barcodeid"_removeFp.fastq"
cutadapt -b GTCGTCCACGTAGCCCAC $runid"_"$barcodeid"_removeFp.fastq" > $runid"_"$barcodeid"_removeFpm.fastq"

cutadapt -b CTTCTACCCTGCGGAGATCA $runid"_"$barcodeid"_removeFpm.fastq" > $runid"_"$barcodeid"_removeFRp.fastq"
cutadapt -b TGATCTCCGCAGGGTAGAAG $runid"_"$barcodeid"_removeFRp.fastq" > $runid"_"$barcodeid"_removeFR.fastq"

perl pipeline/fastq_remove_readlength.pl $length $runid"_"$barcodeid"_removeFR.fastq" $runid"_"$barcodeid"_"$length"bp.fastq"

perl pipeline/fastq_remove_quality.pl 20 90 $runid"_"$barcodeid"_"$length"bp.fastq" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fastq"

perl pipeline/fastq2fasta.pl $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fastq" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta"

blastn -db $blast_reference_fasta -query $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" -evalue 1e-10 -perc_identity 95 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95.txt"

python3 pipeline/count_multihit_v5.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" $species

python3 pipeline/count_uniquehit_v5.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" "count_multihit_"$length".txt" $species


###denovo detection
mkdir ref_for_denovo
cp $blast_reference_fasta ref_for_denovo/
makeblastdb -in ref_for_denovo/$blast_reference_fasta_filename -dbtype nucl

python3 pipeline/detect_denovo.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" ref_for_denovo/$blast_reference_fasta_filename

python3 pipeline/detect_denovo_convert_seq.py ref_for_denovo/$blast_reference_fasta_filename $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta"

blastn -db $blast_reference_fasta -query denovo_candidate1.fasta -evalue 1e-10 -perc_identity 99 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out denovo_candidate1_blast_results.txt

python3 pipeline/detect_denovo_remove_known_alleles.py

###blast to denovo candidate
makeblastdb -in denovo_candidate2.fasta -dbtype nucl

blastn -db denovo_candidate2.fasta -query $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" -evalue 1e-10 -perc_identity 95 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95_DeNovo.txt"

python3 pipeline/count_multihit_v5_denovo.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc95_DeNovo.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" denovo1 $species

python3 pipeline/detect_denovo_remove_redundant.py

#assemble candidate
makeblastdb -in detect_denovo_out5.fasta -dbtype nucl

blastn -db detect_denovo_out5.fasta -query detect_denovo_out5.fasta -evalue 1e-10 -perc_identity 99 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out detect_denovo_out5_blast_results.txt

python3 pipeline/detect_denovo_assemble_and_add_ref.py $blast_reference_fasta detect_denovo_out5.fasta detect_denovo_out5_blast_results.txt reference_known_and_denovo.fasta on

###blast to known and denovo denovo_candidate1
makeblastdb -in reference_known_and_denovo.fasta -dbtype nucl

blastn -db reference_known_and_denovo.fasta -query $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" -evalue 1e-10 -perc_identity 98 -max_target_seqs 10000 -outfmt 6 -num_threads 4 -out $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98_DeNovo2.txt"

python3 pipeline/count_multihit_v5_denovo.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98_DeNovo2.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" denovo2 $species

python3 pipeline/count_uniquehit_v5_denovo.py $length $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent"_blast_eval1e10_perc98_DeNovo2.txt" $runid"_"$barcodeid"_"$length"bp_q"$qv"p"$percent".fasta" "count_multihit_"$length"_denovo2.txt" denovo2 $species

python3 pipeline/merge_sort.py $sampleid
