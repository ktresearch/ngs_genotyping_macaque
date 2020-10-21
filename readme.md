# README

---

This repository contains all code to reproduce the analysis in the article "Development and evaluation of a rapid and cost-efficient NGS-based MHC class I genotyping method for macaques by using a prevalent short-read sequencer".

### requirements:

- blastn  
- cutadapt  
  https://cutadapt.readthedocs.io/en/stable/  
- python (version >= 3.7)  
- perl (version5)  

### How to run codes:

(1) Please copy MHC_genotyping or MHC_denovo(for novel allele detection) directory to your computer.  
(2) Please copy fastq file to MHC_genotyping directory.  
(3) Please move MHC_genotyping or MHC_denovo directory and run .sh file like below.  

> $ sh pipeline/MHC_pipeline_v5.sh  
> $ sh pipeline/MHC_pipeline_denovo_v3.sh (for novel allele detection)    

(4) Please enter information following the instructions on the screen.  

`enter Fastqfile:` Please enter Fastq file name.  
`enter Mamu or Mafa:` Please enter species Mamu (Rhesus Macaque) or Mafa (Cynomolgus Macaque).  
`enter runid (any characters):` Please enter any characters.  
`enter barcodeid (any characters):` Please enter any characters, with a barcode ID number recommended for Ion sequencer.  
`enter samplename (any characters):` Please enter any characters.  

(5) result file name is "(samplename entered from terminal)_result.txt"  
