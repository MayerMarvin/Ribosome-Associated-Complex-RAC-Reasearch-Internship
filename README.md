# BPC Research internship 
Internship from 13.06.2022 - 29.07.2022 at Bukau lab at the University in Heidelberg
Supervisor Koji Ishikawa

Topic of the research intership is the analysis of the ribosome associated complex (RAC) which consits of Zuotin and Ssz1 binding to the ribosome. It is responsible for binding the emerging nascant chain and for the handover of the nascant chain towards Ssb chaperones.

### Cytosin Metagene Analysis
One hypothesis put forward by Koji Ishikawa is that RAC binds to hydrophobic amino acids around the first cysteine in the amino acid sequence. For this purpose, the script [aa_distr_first_cys.py](https://github.com/MayerMarvin/Ribosome-Associated-Complex-RAC-Reasearch-Internship/blob/4937c5910a7763a01b1b363998778b2297e4d139/Cytosine%20Metagene%20Analysis/aa_distr_first_cys.py) was created to determine the amino acids located around the first cysteine. Then, only the amino acids that are significantly increased are displayed. This is determined by comparing the abundance of the amino acid at this relative position to the first cysteine, with the normal distribution of the amino acid in the yeast proteome. Significantly elevated amino acids are then output in color coding, with positively charged amino acids in orange, negatively charged amino acids in yellow, polar and uncharged amino acids in red, special case amino acids in green, and hydrophobic amino acids in blue. 

![50_window_first_cys](https://github.com/MayerMarvin/Ribosome-Associated-Complex-RAC-Reasearch-Internship/blob/e38bca87ff728864db8b2690e785ccd5b61b6ade/Cytosine%20Metagene%20Analysis/50_window_first_cys.png)

### Bioinformatic Analysis of Selective Ribosome Profiling (SeRP)

For the analysis of the SeRP data, a customized pipeline by Koji Ishikawa was used. This pipeline runs the data analysis and quality assessment of Ribosome Profiling Data. It analyzes short unpaired raw reads requiring fastq files as input and runs following analysis:

  1. Trimming of 3' adapter with Cutadapt (v.1.13)

  2. Umis isolation with custom Julia script (by Ilia)

  3. rRNA sequences removal with Bowtie2 (v.2.3.5.1)

  4. Alignment of reads using STAR (on the Nextseq computer, v2.7.1a)

  5. Quality check with Fastqc (v.0.11.5)

  6. Reads assignemnt (5') with custom Julia script (by Ilia), including umi-collapsed assignment and the read length correction of soft-clipped reads with "Non-templated Nt additions" (updated version from 10th of Oct 2019)
  
  The produced results can be viewed in the XYZ folder
