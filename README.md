# NanoTel
### This program searches DNA sequences for telomeric patterns and includes 5 outputs:
 1. A folder named **reads** with all telomere-containing DNA reads in separate FASTA files. 
 2. A detailed CSV-format table named **summary** with the telomere length, location and ratio of telomeric patterns each telomere contains for each DNA sequence.
 3. A folder named **single_read_plots** containing plots of the telomere-pattern density throughout the DNA sequence for all sequences, normalized to 100kb, aiding telomere length comparison.
 4. A similar folder named **single_read_plots_adj** with the same plots, all spread out so that the x axis contains only the DNA read length (usually under 100kb).
 5. A folder named **log** with summary statistics of the analysis run.
 
 Example of hypothetical plot in the *single_read_plots* folder:
![plot_example](https://github.com/Tzfatilab/Telomere-Analyzer/blob/main/Example/graph_example.jpeg)



## Contents

- [Instructions](#instructions)
- [Run example](#example)
- [Changing default parameters](#changing-default-parameters) 
- [Required preinstallations](#preinstallations)
- [Citing NanoTel](#citing-nanotel)

### Instructions

If workig on a Linux OS, use the *NanoTel* file to run the program, otherwise use *nanotel_old.R* from the v1.0.0-beta. The difference between the two is the use of parallel computing in the former which speeds up the computation process and is not currently supported by Windows.  
All other files for now are just a draft with remarks for future planning.  
Tests were done on Ubuntu operating system version 22.04.1.

This script is suitable for running in Linux and demands 4 arguments: (1) The code file used for analysis, (2) DNA sequences to be tested in FASTA or FASTQ format, (3) An output directory, (4) Specification of sequences file type. 
  
To run it, open shell and use the command scaffold:  `Rscript --vanilla NanoTel.R input_dir output_dir file_format`  
Replace parameters as following:
- NanoTel.R: path for the code file (including file name).
- input_dir: path for the fastq/a file or directory containing fasta/fastq files.
- output_dir: path for the output directory.
- file_format: type *fasta* or *fastq* (default is fastq).  

**Make sure the output_dir is not a subdirectory of input_dir or vice versa**

After executing the code, a question will be asked `Use reverse complement?`. This refers to the default telomeric pattern which is searched. If the sequenced DNA is from the C-rich strand, type *yes*, otherwise type *no*.  
A second question will be asked `Use the filtration?`. This question refers to filteration of reads that are under 1000 nt and that have a low telomeric density at the edge of the read. Type *yes* if you want to filter, and type *yes* again after the next question if the the sequenced DNA is from the C-rich strand.

### Example 
`git clone https://github.com/Tzfatilab/Telomere-Analyzer.git`   
`cd Telomere-Analyzer`  
`Rscript --vanilla NanoTel.R Example/10_reads.fasta Example/Example_output fasta`  
Two questions will be asked. Type `no` in both.  
Expected run time is 11 seconds.
  
### Changing default parameters  
Before running the code, consider changing certain parameters in the code itself that by default have set values:
- The telomere pattern density is searched in segments of 100 consecutive bases. This could be changed in the `global_subseq_length` parameter.
- Each segment is classified as potentially telomeric or not depending if it passes a minimum density threshold. Changing the existing threshold (0.3) could be done in the `global_min_density` parameter. This could affect where the program sets the telomere beginning.  
- The telomeric sequence searched is TTAGGG. This could be changed under the `curr_patterns` variable.
- If you used guppy basecaller with the r9 kit it is best recommended to re-basecall with dorado due to problems in purine distinction by the sequencing technology or changing the `curr_patterns` variable to a list of "YYAGGG/CCCTRR" patterns.

  
### Preinstallations  
The code was written in R 4.2.2 version. The following R packages should be preinstalled (including versions used at the time of publishing): conflicted (1.2.0), tidyverse (2.0.0), logr (1.3.3), future (1.32.0), ggprism (1.0.4), IRanges (2.32.0), and Biostrings (2.66.0), while the last two are part of Bioconductor and should be installed through the BiocManager package. 

### Citing NanoTel 
If you use NanoTel in your work, please cite:  
Smoom R., May C.L., Ortiz V., Tigue, M., Kolev H., Reizel Y., Morgan A., Egyes, N., Lichtental D., Skordalakes E., Kaestner K.H., & Tzfati Y. A single amino acid in RTEL1 determines telomere length in mice. BioRxiv:  https://doi.org/10.1101/2021.06.06.447246
