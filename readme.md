# README

---

This repository contains all code to reproduce the analysis in the article "Development and evaluation of a rapid and cost-efficient NGS-based MHC class I genotyping method fro macaques by using a prevalent short-read sequencer".

### requirements:

- blastn  
- cutadapt  
  https://cutadapt.readthedocs.io/en/stable/  
- python (version >= 3.7)  

### How to run codes:

(1) Please copy MHC_genotyping directory to your computer.  
(2) Please copy fastq file to MHC_genotyping directory.  
(3) Please move MHC_genotyping directory and run MHC_pipeline_v5.sh like below.  

> $ sh pipeline/MHC_pipeline_v5.sh  

(4) Please enter information following the instructions on the screen.  

`input Fastqfile:` Please enter Fastq file name.  
`input Mamu or Mafa:` Please enter species Mamu (Rhesus Macaque) or Mafa (Cynomolgus Macaque).  
`input runid (any characters):` Please enter any characters.  
`input barcodeid (any characters):` Please enter any characters, with a barcode ID number recommended for Ion sequencer.  
`input samplename (any characters):` Please enter any characters.  
