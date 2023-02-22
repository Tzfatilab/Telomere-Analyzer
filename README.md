# Telomere-Analyzer
 Search over DNA sequences for telomeric patterns.
 
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
