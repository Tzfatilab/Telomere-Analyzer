# Telomere Analyzer
### This program searches DNA sequences for telomeric patterns and includes 5 outputs:
 1. A folder named **reads** containing all DNA reads possesing a telomere in separate FASTA files. 
 2. A detailed CSV-format table named **summary** with the telomere length, location and ratio of telomeric patterns it contains for each DNA sequence.
 3. A folder named **single_read_plots** containing plots of the telomere-pattern density throughout the DNA sequence for all sequences.
 4. A similar folder named **single_read_plots_adj** with the same plots, all adjusted to 100 kb, aiding telomere length comparison.
 5. A folder named **log** with summary statistics on the analysis run.
 
If workig on a Linux OS, use the *nanotel-multicore-10workers.R* to run the program, otherwise use *nanotel.R*. The difference between the two is the use of parallel computing in the former which speeds up the computation process.
All other files for now are just a draft with remarks for future planning.

This script is suitable for running in Linux and demands 4 arguments: (1) The code file used for analysis, (2) DNA sequences to be tested in FASTA or FASTQ format, (3) An output directory, (4) Specification of sequences file type. 
  
To run it, open shell and use the command scaffold:  *Rscript --vanilla nanotel-multicore-10workers.R input_dir output_dir file_format*  
Replace parameters as following:
- nanotel-multicore-10workers.R: path for the code file (including file name).
- input_dir: path for the fastq/a file or directory containing fasta/fastq files.
- output_dir: path for the output directory.
- file_format: type *fasta* or *fastq* (default is fastq).  

**Make sure the output_dir is not a subdirectory of input_dir or vice versa**




Before running the nanotel: check the default args of the functions: you may want to change some arguments such as the pattern, min_density etc...


The script is supported on Linux OS.
