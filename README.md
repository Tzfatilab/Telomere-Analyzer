# Telomere-Analyzer
 This program searches DNA sequences for telomeric patterns and includes 4 outputs:
 1. A FASTA file containing only DNA reads possesing a telomere
 2. A detailed CSV-format table with the telomere length, location and ratio of telomeric patterns it contains for each DNA sequence.
 3. A folder named **single_read_plots** containing plots of the telomere-pattern density throughout the DNA sequence for all sequences.
 4. A similar folder named **single_read_plots_adj** with the same plots, all adjusted to 100 kb, aiding telomere length comparison.
 
Use the nanotel-multicore-10workers.R to run the program.
All other files for now is just a draft with a lot of remarks for future planning.

The nanotel is a script suitable for running in Linux:
Accept 2 or 3 arguments:
arg[1] : input of fastq/fasta file/dir
arg[2]: The output directory
arg[3]: The format of the files: fastq (default) or fasta

The nanotel-multicore-10workers is the same script as nanotel except that it run in parallel using future::plan(multicore, workers = 10).
The nanotel-multicore-10workers is not supported on Windows OS.

Before running the nanotel: check the default args of the functions: you may want tochange some arguments such as the pattern, min_density ect...

To run it: use the commane:
" Rscript --vanilla nanotel-multicore-10workers.R input_dir output_dir file_format "
- input_dir: path for the fastq/a file or directory containing fasta/fastq files
- output_dir: path for the output directory
- file_format: fasta or fastq (default is fastq)
**** Make sure the output_dir is not a subdirectory of input_dir or vice versa.
 

The script is supported on Linux OS.
