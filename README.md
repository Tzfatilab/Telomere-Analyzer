# Telomere-Analyzer
 Search over DNA sequences for telomeric patterns.
 
 
TelomereAnalyzer.R for now is just a draft with a lot of remarks for future planning.
Use the TelomereAnalyzer_with_running_samples.R to run the program.

The nanotel is a script suitable for running in Linux:
Accept 2 or 3 arguments:
arg[1] : input of fastq/fasta file/dir
arg[2]: The output directory
arg[3]: The format of the files: fastq (default) or fasta

Before running the nanotel: check the default args of the functions: you may want tochange some arguments such as the pattern, min_density ect...

To run it: use the commane:
" Rscript --vanilla nanotel.R input_dir output_dir file_format "
 

The script is supported on Linux OS.
