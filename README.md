# NanoTel
### This program searches DNA sequences for telomeric patterns and includes 5 outputs:
 1. A folder named **reads** with all telomere-containing DNA reads in separate FASTA files. 
 2. A detailed CSV-format table named **summary** with the telomere length, location and ratio of telomeric patterns each telomere contains for each DNA sequence.
 3. A folder named **single_read_plots** containing plots of the telomere-pattern density throughout the DNA sequence for all sequences, normalized to 100kb, aiding telomere length comparison.
 4. A similar folder named **single_read_plots_adj** with the same plots + patterns density with 1 mismatch allowed , all spread out so that the x axis contains only the DNA read length (usually under 100kb).
 6. A folder named **log** with summary statistics of the analysis run.
 7. A file named **reads_ids.txt** with the ids of the telomeric reads.
 
 Example of hypothetical plot in the *single_read_plots* folder:
![plot_example](https://github.com/Tzfatilab/Telomere-Analyzer/blob/main/Example/graph_example.jpeg)



## Contents

- [Instructions](#instructions)
- [Run example](#example)
- [Changing default parameters](#changing-default-parameters) 
- [Required preinstallations](#preinstallations)
- [Citing NanoTel](#citing-nanotel)

### Instructions

The scripts s suitable for running in Linux OS.
Tests were done on Ubuntu operating system version 22.04.1.

This script i and demands 3 arguments: (1) The code file or dirctory containig files used for analysis, (2) An output directory, (3) Specification of the pattern/s to search. 
  
To run it, open shell and use the command scaffold:  `Rscript --vanilla NanoTel.R -i input_dir --save_path output_dir --patterns pattern_list`  
Replace parameters as following:
- NanoTel.R: path for the code file (including file name).
- input_dir: path for the fastq/a file or directory containing fasta/fastq files.
- output_dir: path for the output directory.
- pattern_list: A single pattern (for example: TTAGGG) or a list of patterns(Must be in double quotes) for example: "TTAGGG TTGGG "CCAGGG".
  It is possible to give a multi pattern as input for eample: "TTAGGG TCAGGG CTAGGG CCAGGG", but when it is possible we recommand to give it as a single pattern with the ambiguity letters:
  For eample: "TTAGGG TCAGGG CTAGGG CCAGGG" should be replaced by the single pattern YYAGGG.
  
### additional arguments: 
--format=FILE FORMAT
	input files format (Either "fastq" (the default) or "fasta", gzip is supported)

-n NREC, --nrec=NREC
	Single integer. The maximum of number of records to read in to memory for each iteration. Negative values are ignored.

-r, --rc
	Should we do reverse complement on the given reads.
	
--min_density=MINIMAL DENSITY.
	Minimal density to consider a subsequence as a pattern region.
	
--subseq_length=SUB-SEQUENCE LENGTH.
	The length of the sub-sequence.
	
--use_filter
	Filter reads accoding to the edge.
	
--check_right_edge
	When using the filter function, check the start of the sequence (left) or the end of the sequence (right), using this flag will tell the filter to check the right flag, the default will check the left edge!

--tvr_patterns= TELOMERE VARIANT REPEATS PATTERNS
	Space separated list of additional pattern(s) for Telomere variant repeats. Must be in double quotes.

	

	
 Sometimes we need to add  TVR pattersns for more accurate Telomere analysis. (see the different between 2 plots)
 ![plot_example](https://github.com/Tzfatilab/Telomere-Analyzer/blob/main/Example/noTVR.jpeg)
 ![plot_example](https://github.com/Tzfatilab/Telomere-Analyzer/blob/main/Example/withTVR.jpeg)
  

**Make sure the output_dir is not a subdirectory of input_dir or vice versa**
for additional available input flags: see Rscript --vanilla NanoTel.R --help


### Example 
`git clone https://github.com/Tzfatilab/Telomere-Analyzer.git`   
`cd Telomere-Analyzer`  
`Rscript --vanilla NanoTel.R -i Example/sample.fasta --save_path Example/Example_output patterns TTAGGG --min_density 0.6`
  
  
### Preinstallations  
The code was written in R 4.2.2 version. The following R packages should be preinstalled (including versions used at the time of publishing): conflicted (1.2.0), tidyverse (2.0.0), logr (1.3.3), future (1.32.0), ggprism (1.0.4), IRanges (2.32.0), optparse and Biostrings (2.66.0), while the last two are part of Bioconductor and should be installed through the BiocManager package. 

### Citing NanoTel 
If you use NanoTel in your work, please cite:  
Smoom R., May C.L., Ortiz V., Tigue, M., Kolev H., Reizel Y., Morgan A., Egyes, N., Lichtental D., Skordalakes E., Kaestner K.H., & Tzfati Y. A single amino acid in RTEL1 determines telomere length in mice. BioRxiv:  https://doi.org/10.1101/2021.06.06.447246


# chrMap
### This script help to filter NanoTel results with minimap2 alignment using The dorado aligner summary file. (works for dorado version 1.3.1 or newer +  and NanoTel version v1.1.4-beta or newer.)
See details using the --help flag.
Still under developments...






