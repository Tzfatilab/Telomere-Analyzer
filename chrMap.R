#!/usr/bin/env Rscript


# Given nanotel output + minimap2 dorado summaryfile -> extract filter with a
# dir for the mapped

# according to dorado aligner ( dorado version 1.3.1) +  and NanoTel version v1.1.4-beta

# last updated 2026-02-02


# todo: check if it has tvr_pattern columns....


#min_mapq <- 60 # add several parameters
# min_num_align <- 6000
# filter for each param : if not pass add why

#nanotel_path <- "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN"
library(optparse)


#' todo:
#' 1. Add log file
#' 2. Add vizualizations
#' 






# optional filterations:
#' 1.  alignment_genome != '*' This is must filter 
#' 2.  alignment_direction     v
#' 3.  alignment_genome_start  v
#' 4.  alignment_genome_end    v
#' 5.  alignment_strand_start  x
#' 6.  alignment_strand_end    x
#' 7.  alignment_identity      [0,1] - do not use 
#' 8.  alignment_accuracy      v [0,1]
#' 9.  alignment_coverage :    v ( claculate sub-telo coverage and see if it is ~ alignment_coverage(alignment_coverage = (aligned length of query) / (total query length) × 100)
#' 10. alignment_mapping_quality v [0,60] int 




# creat for each filter parameter a logical column if passed filter
#' for exmple: Head ~ 1, tail ~ -1 -> pass_direction  
#' 
#' 
#' 
#' todo: add function which do reverese Complement to reads which are mapped to 3'
  
###################  ###########################################################  


option_list = list(
  
  # telo summary input_path
  make_option("--telo_summary_path", action = "store", type = "character" , 
              default = NULL, help = "The path for the NanoTel summary.csv file. ", metavar = "NanoTel summary path"),
  
  # NanoTel input path 
  make_option("--nanotel_path", action = "store", type = "character" , 
              default = NULL, help = "The path for the NanoTel output dir. ", metavar = "NanoTel output path"),
  
  # dorado aligner --emit-summary data frame path
  make_option("--aligner_summary_path", action = "store", type = "character" , 
              default = NULL, help = "The path for alignments data frame. ", metavar = "alignments sequencing summary path"),
  
  # save_path
  make_option("--save_path", action = "store", 
              type = "character" , default = NULL, 
              help = "A path to a directory for storing the output files.", 
              metavar = "output dir path"),
  
  make_option("--filter_direction", action = "store_true", default = FALSE,
              help = "Filter alignments according to alignment_direction column: We asssume that reads mapped to Head are + and mapped to Tail are -.", 
              metavar = "USE filter_direction filter" ), 
  
  # if <0 ignore
  make_option("--filter_genome_position", action = "store", default = NULL,
              help = "Filter alignments according to alignment_genome_start/genome_end, with threshold of given integer >= 0:.", 
              type = 'integer',metavar = "USE filter_direction filter with a given threshold" ), 
  
  
  make_option("--min_alignment_accuracy", action = "store", default = NULL, 
              type = "double", 
              help = "Minimal alignment accuracy (alignment_accuracy = (number of matches) / (matches + mismatches + insertions + deletions) × 100).",
              metavar = "Minimal alignment accuracy score."), 
  
  
  # todo: need to take additional arg if to check with telo_mismach, telo or telo_tvr indices
  make_option("--min_alignment_coverage_thr", action = "store", default = NULL, 
              type = "double", 
              help = "Minimal threshold for The subtelomere coverage: abs(alignment_coverage - subtelomere_length/read_length) <= threshold.",
              metavar = "Alignmnet coverage threshold"), 
  
  # todo : add a logical which tells if it is right or left
  make_option("--telo_index", action = "store", default = "telomere", type = "character", 
              help = "The indices for computing the subtelomere length: choose from Telomere_start/end, Telomere_start_mismatch/end_mismatch or Telomere_start_mismatch_tvr/end_mismatch_tvr: the argument must be one of the 3: telomere, mismatch or tvr!",
              metavar = "Column for computing subtelomere length."),
  
  make_option("--telo_right", action = "store_true", default = FALSE,
              help = "Use this flag if the Telomere should be at the right edge of the read otherwise it will assume the telomere is at the left edge, helps for calculating the sub-telo length accordingly.",
              metavar = "Telomere position at the right edge."), 
  
  make_option("--min_alignment_mapping_quality", action = "store", default = NULL, 
              type = "integer", 
              help = "Minimal mapping quality (should be an integer in [0,60].",
              metavar = "Minimap alignmnet mapping quality."), 
  
  make_option("--genome_edges_length", action = "store", default = NULL, 
              type = "integer", 
              help = "The length of the reference genome edges(Heads/Tails) of the chromosomes.",
              metavar = "Genome edge length."), 
  
  make_option("--version", 
              action = "store_true", 
              default = FALSE,
              help = "Print version information and exit"), 
  
  # opt$file_extension for reads (now support fasta or fasta.gzip)
  make_option("--file_extension", action = "store", default = ".fasta", type = "character", 
              help = "File extension for HTS reads (for now only fasta(default) or fasta.gzip", 
              metavar = "fasta or fasta.gzip"), 
  # subtelo_length_thr
  make_option("--subtelo_length_thr", action = "store", default = 4000L, 
              type = "integer", 
              help = "The minimal subtelomeric length to pass mapping filteration with defaul length of 4000, insert 0 or lower for not using it.")
)   

opt = parse_args(OptionParser(option_list=option_list))  
  


# Handle --version flag
if (opt$version) {
  cat("Telomere Analyzer  version v1.1.7-beta 2026-02-18\n")
  quit(save = "no", status = 0)
}  



library(readr)
suppressPackageStartupMessages(require(dplyr))
library(stringr)
library(conflicted)
library(logr)

#' @description
#'  Create the input dir + subdirectories from the mapping
#' 
#' @param path the path for the dir we want to create.
#' @param alignments A vector of unique names for the mappings.
create_dirs <- function(path, alignments)
{
  if(!dir.exists(path)){
    dir.create(path, showWarnings = TRUE)
  }
  for(chr in alignments)
  {
    dir.create(path = str_c(path, chr, sep = "/"))
  }
  # create the unclassified
  dir.create(path = str_c(path, "unclassified", sep = "/"))
}

#' @description
#' copy reads from input path to output path according to certain alignment 
#' genome: This function will create a sub-dir for the current alignment and 
#' will add all the data of the reads which were mapped to it.
#'
#' @param df_merged: A data frame of minimap2 results (from dorado aligner 
#' combined with the summary file from NanoTel after passing mapping filteration). 
#' @param chrs: String rep for the alignments (vector of unique names for all the alignments).
#' @param nanotel_path: path for the NanoTel dir.
#' @param chrs_path: path for the inputs.
#' @param file_extension: The extension of the reads.  
#' @param unclassified_serial The reads which were not mapped or did not pass the filter.
#' (vector od serial numbers)
#'
#' @return 0 if all files were copied 
copy_reads <- function(df_merged, chrs, nanotel_path, chrs_path, file_extension = c(".fasta", ".fasta.gz") , unclassified_serial )
{
  file_extension <- match.arg(file_extension)
  
  res <- logical()
  for(chr in chrs) { ############################## todo: add pass_all column
    serial <- df_merged %>% dplyr::filter(alignment_genome == chr & pass_all) %>% 
      pull(Serial)
    
    # make reads dir
    base::dir.create(path = str_c(chrs_path,chr ,'reads', sep = '/') , showWarnings = TRUE)
    
    from_path <- str_c(nanotel_path, '/reads/',   as.character(serial), file_extension)
    to_path <- str_replace_all(string = from_path, pattern = nanotel_path, replacement = str_c(chrs_path,chr, sep = '/'))
    res <- c(res, file.copy(from = from_path  ,to = to_path, overwrite = TRUE) )
  }
  
  # unclassified
  base::dir.create(path = str_c(chrs_path,'unclassified' ,'reads', sep = '/') , showWarnings = TRUE)
  from_path <- str_c(nanotel_path, '/reads/',   as.character(unclassified_serial), file_extension)
  to_path <- str_replace_all(string = from_path, pattern = nanotel_path, replacement = str_c(chrs_path,'unclassified', sep = '/'))
  
  res <- c(res, file.copy(from = from_path  ,to = to_path, overwrite = TRUE) )
  return(sum(!res))
}



#' @description
#' copy plots from input path to output path according to certain alignment 
#' genome: This function will create a sub-dir for the current alignment and 
#' will add all the data of the reads which were mapped to it.
#'
#' @param df_merged: A data frame of minimap2 results (from dorado aligner 
#' combined with the summary file from NanoTel after passing mapping filteration). 
#' @param chrs: String rep for the alignments (vector of unique names for all the alignments).
#' @param nanotel_path: path for the NanoTel dir.
#' @param chrs_path: path for the inputs.
#' @param file_extension: The extension of the reads.  
#' @param unclassified_serial The reads which were not mapped or did not pass the filter.
#' (vector od serial numbers)
#'
#' @return 0 if all files were copied 
copy_plots <- function(df_merged, chrs, nanotel_path, chrs_path, 
                       file_extension =c('.jpeg', '.eps'), 
                       dir_name = c('single_read_plots', 'single_read_plots_adj'), unclassified_serial )
{
  file_extension<- match.arg(file_extension)
  dir_name <- match.arg(dir_name)
  res <- logical()
  
  
  for(chr in chrs) { ############################## todo: add pass_all column
    serial <- df_merged %>% dplyr::filter(alignment_genome == chr & pass_all) %>% 
      pull(Serial)
    
    # make reads dir
    if(!dir.exists(paths = str_c(chrs_path,chr, dir_name, sep = '/') )) {
      base::dir.create(path = str_c(chrs_path,chr, dir_name, sep = '/') , showWarnings = TRUE)  
    }
    
    
    from_path <- str_c(nanotel_path, '/',dir_name, '/read', as.character(serial), file_extension)
    to_path <- str_replace_all(string = from_path, pattern = nanotel_path, replacement = str_c(chrs_path,chr, sep = '/'))
    
    res <- c(res, file.copy(from = from_path  ,to = to_path, overwrite = TRUE) )
  }
  
  # unclassified
  if(!dir.exists(paths = str_c(chrs_path,'unclassified', dir_name, sep = '/'))) {
    base::dir.create(path = str_c(chrs_path,'unclassified', dir_name, sep = '/') , showWarnings = TRUE) 
  }
  
  
  from_path <- str_c(nanotel_path, '/',dir_name, '/read', as.character(unclassified_serial), file_extension)
  to_path <- str_replace_all(string = from_path, pattern = nanotel_path, replacement = str_c(chrs_path,'unclassified', sep = '/'))
  
  res <- c(res, file.copy(from = from_path  ,to = to_path, overwrite = TRUE) )
  return(sum(!res))
}



#' @description
#' Load the data frames from NanoTel and dorado aligner summary files and 
#' join them accordingly.
#' @param nanotel_summary_path: 
#' @param minimap_path: 
#' @return joined df of all the telomeric reads
join_df <- function(nanotel_summary_path, minimap_path){
  
  # for minimap
  col_types <- "ciccfiiiiiliiiidffffffiiiiiiiiddidiiiii"
  
  col_names <- c('read_id','alignment_genome', 
                 'alignment_direction', 'alignment_genome_start', 'alignment_genome_end', 
                 'alignment_strand_start', 'alignment_strand_end', 'alignment_num_insertions',
                 'alignment_num_deletions', 'alignment_num_aligned','alignment_num_correct', 
                 'alignment_identity', 'alignment_accuracy', 'alignment_score', 
                 'alignment_coverage',  'alignment_mapping_quality', 
                 'alignment_num_alignments', 'alignment_num_secondary_alignments'
  )           
  
  # first rename to read_id : I should change it in the new NanoTel
  df_telo <- read_csv(file = nanotel_summary_path, col_types = "icidiiidiiidiii")
  df_telo <- df_telo %>% rename(read_id = sequence_ID)
  df_telo$read_id <- str_sub(df_telo$read_id, start = 1 , end = 36)
  
  
  # Process chunks with a function
  process_chunk <- function(chunk, pos) {
    filtered <- chunk %>% dplyr::filter(read_id %in% df_telo$read_id)
    return(filtered)
  }
  
  df_minimap <- read_tsv_chunked(minimap_path, 
                                 DataFrameCallback$new(process_chunk),
                                 chunk_size = 10000, col_types = col_types) %>% 
    select(all_of( col_names))
  
  # adjust index fro R ( start <- start + 1)
  df_minimap$alignment_genome_start <- case_when(
    df_minimap$alignment_genome_start == -1 ~ -1,
    .default = df_minimap$alignment_genome_start + 1 )
  df_minimap$alignment_strand_start <- case_when(
    df_minimap$alignment_strand_start == -1 ~ -1, 
    .default = df_minimap$alignment_strand_start + 1 )
  
  df_join <- full_join(x = df_telo, y = df_minimap, by = "read_id")
  
  return(df_join)
}





# todo: load all the other packages after this line inorder to save time!
# todo: For -1 , find alternative: claculate instead with the valid indices 
# todo: check if it works!!!!

#' @description
#'  Create a sub-telomere column which gives the length of the subtelomeric part of the read.(length == -1 if the indices for the calculation are NA)
#' @param df The data frame (summary file of nanotel or more).
#' @param telo_index The columns to choose for the subtelomeric position:Telomere_start(end)/Telomere_startTelomere_start_mismatch_mismatch/Telomere_startTelomere_start_mismatch_mismatch_tvr. 
#' @param telo_position Is the telomere should be at the start or end(right) of the read.
#' @return The data frame with the subtelo_length column.
calculate_subtelo <- function(df, telo_index=c("telomere", "mismatch", "tvr"), telo_position=c("left", "right")) {
  telo_position <- match.arg(telo_position)
  telo_index <- match.arg(telo_index)
  
  if(telo_index == "telomere" && telo_position == "left") {
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_end), -1, sequence_length -Telomere_end ) )
    
  } else if(telo_index == "telomere"  && telo_position == "right") {
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_start), -1, sequence_length - Telomere_start + 1) )
    
  } else if(telo_index =="mismatch" && telo_position == "left") {
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_end_mismatch), -1, sequence_length -Telomere_end_mismatch) )
    
  } else if(telo_index ==  "mismatch"&& telo_position == "right") {
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_start_mismatch), -1, sequence_length - Telomere_start_mismatch + 1) )
    
  } else if(telo_index == "tvr" && telo_position =="left") {
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_end_mismatch_tvr), -1, sequence_length -Telomere_end_mismatch_tvr) )
    
  } else { # tvr + right
    df <- df %>% 
      mutate(subtelo_length = ifelse(is.na(Telomere_start_mismatch_tvr), -1, sequence_length - Telomere_start_mismatch_tvr + 1) )
  }
  return(df)
}  







# check flags: must flags
if (is.null(opt$save_path)) {
  stop("Missing required parameter:  --save_path", call.=FALSE)
}

if(!dir.exists(opt$save_path)) { # update  did it
  dir.create(opt$save_path)
}

tmp <- file.path(opt$save_path, "run.log")

# Open log
lf <- log_open(tmp)
log_print('Telomere Analyzer  version v1.1.7-beta 2026-02-18', hide_notes = TRUE, console = FALSE) 

# optional filterations:
#' 1.  alignment_genome != '*' This is must filter 
#' 2.  alignment_direction     v
#' 3.  alignment_genome_start  v
#' 4.  alignment_genome_end    v
#' 5.  alignment_accuracy      v [0,1]
#' 6.  alignment_coverage :    v ( claculate sub-telo coverage and see if it is ~ alignment_coverage(alignment_coverage = (aligned length of query) / (total query length) × 100)
#' 7.  alignment_mapping_quality v [0,60] int   
mapping_filter <- function(df_join, filter=NULL,filter_column='alignment_genome' , thr=NULL, genome_length=NULL, telo_index="telomere", telo_right=FALSE)
{
  if(is.null(filter)){ return(df_join)}
  if(filter==FALSE){ return(df_join)}
  
  if(filter_column == 'alignment_genome') {
    df_join <- df_join %>% 
      mutate(pass_alignment_genome = (alignment_genome != '*') )
    log_print(paste(sum(df_join$pass_alignment_genome), "reads pass the alignment filteration!"), console = FALSE, hide_notes = TRUE)
    return(df_join)
  }
  
  # We assume heads are + and atails are -
  if(filter_column == 'alignment_direction') {
    df_join <- df_join %>% 
      mutate(pass_alignment_direction = 
               ( (str_detect(string = alignment_genome, pattern = "Head") & alignment_direction == "+" ) | 
                   (str_detect(string = alignment_genome, pattern = "Tail") & alignment_direction == "-" ) ) )
    log_print(paste(sum(df_join$pass_alignment_direction), "reads pass the alignment direction filteration!"), console = FALSE, hide_notes = TRUE)
    return(df_join)
  }
  
  
  # This is to be detected!
  if(filter_column == 'alignment_genome_start') {
    if(is.null(genome_length)) {
      warning("No genome length givven so cab filter genome indices! Make sure to use the --genome_edges_length flag!")
      return(df_join)
    }
    df_join <- df_join %>% 
      mutate(pass_alignment_genome_start_end = (alignment_genome_start <= thr & str_detect(string = alignment_genome, pattern = "Head")) | 
               (abs(alignment_genome_end - genome_length) <= thr & str_detect(string = alignment_genome, pattern = "Tail")) ) 
    log_print(paste(sum(df_join$pass_alignment_genome_start_end), "reads pass the genome position filteration!"), console = FALSE, hide_notes = TRUE)
    return(df_join)
  }
  
  # same as alignment_genome_star
  if(filter_column == 'alignment_genome_end') {
    df_join <- df_join %>% 
      mutate(pass_alignment_genome_start_end = (alignment_genome_start <= thr & str_detect(string = alignment_genome, pattern = "Head")) | 
               (abs(alignment_genome_end - genome_length) <= thr & str_detect(string = alignment_genome, pattern = "Tail")) ) 
    
    return(df_join)
  }
  
  if(filter_column == 'alignment_accuracy') {
    df_join <- df_join %>% 
      mutate(pass_alignment_accuracy = (alignment_accuracy >= thr))
    log_print(paste(sum(df_join$pass_alignment_accuracy), "reads pass the alignment_accuracy filteration of", thr, "!"), console = FALSE, hide_notes = TRUE)
    return(df_join)
  }
  
  # todo: Need toa calculate the sub-telo length
  if(filter_column == 'alignment_coverage') {
    #df_join <- calculate_subtelo(df_join)
    df_join <- df_join %>% 
      mutate(pass_alignment_coverage = (abs(subtelo_length/sequence_length - alignment_coverage) < thr ))
    log_print(paste(sum(df_join$pass_alignment_coverage), "reads pass the alignment coverage filteration of", thr, " which is the diffrence between alignment coverage and sub-telomere coverage!"), console = FALSE, hide_notes = TRUE)
    # sum( abs( (subtelo_length2 / df_pass_genome$sequence_length) - df_pass_genome$alignment_coverage ) <0.2 )
    return(df_join)
  }
  
  if(filter_column == 'alignment_mapping_quality') {
    df_join <- df_join %>% 
      mutate(pass_alignment_mapping_quality = (alignment_mapping_quality  >= thr))
    log_print(paste(sum(df_join$pass_alignment_mapping_quality), "reads pass the alignment mapping quality filteration of", thr, "!"), console = FALSE, hide_notes = TRUE)
    return(df_join)
  }
  
} # End of mapping_filter









if (is.null(opt$telo_summary_path)) {
  log_error("Missing required parameter:  --telo_summary_path")
  log_close()
  writeLines(readLines(lf))
  stop("Missing required parameter:  --telo_summary_path", call.=FALSE)
} else {
  log_print(paste('NanoTel summary path:', opt$telo_summary_path), hide_notes = TRUE, console = FALSE) 
}



if (is.null(opt$nanotel_path)) {
  log_error("Missing required parameter:  --nanotel_path")
  log_close()
  writeLines(readLines(lf))
  stop("Missing required parameter:  --nanotel_path", call.=FALSE)
} else {
  log_print(paste('NanoTel output path:', opt$nanotel_path), hide_notes = TRUE, console = FALSE) 
}

if (is.null(opt$aligner_summary_path)) {
  log_error("Missing required parameter:  --aligner_summary_path")
  log_close()
  writeLines(readLines(lf))
  stop("Missing required parameter:  --aligner_summary_path", call.=FALSE)
} else {
  log_print(paste('Alignment summary path:', opt$aligner_summary_path), hide_notes = TRUE, console = FALSE)
}

# optional flags       min_alignment_mapping_quality min_alignment_accuracy min_alignment_coverage_th
if( !is.null(opt$min_alignment_mapping_quality)) {
  if( (opt$min_alignment_mapping_quality < 0) ||
      (opt$min_alignment_mapping_quality > 60)
  ){log_error("The alignment mapping quality threshold should be an integer in [0,60]!")
    log_close()
    writeLines(readLines(lf))
    stop("The alignment mapping quality threshold should be an integer in [0,60]!")
  } else 
  {
    log_print(paste('Alignment mapping quality threshold:', opt$min_alignment_mapping_quality), hide_notes = TRUE, console = FALSE)
  }
}

if( !is.null(opt$min_alignment_accuracy)) {
  if( (opt$min_alignment_accuracy < 0) ||
      (opt$min_alignment_accuracy > 1)
  ){
    log_error("The alignment accuracy threshold should be a float in [0,1]!")
    log_close()
    writeLines(readLines(lf))
    stop("The alignment accuracy threshold should be a float in [0,1]!")
  } else {
    log_print(paste('Alignment accuracy threshold:', opt$min_alignment_accuracy), hide_notes = TRUE, console = FALSE)
  }
} 

if( !is.null(opt$min_alignment_coverage_thr)) {
  if( (opt$min_alignment_coverage_thr < 0) ||
      (opt$min_alignment_coverage_thr > 1)
  ){
    log_error("The alignment coverage threshold should be a float in [0,1]!")
    log_close()
    writeLines(readLines(lf))
    stop("The alignment coverage threshold should be a float in [0,1]!")
  } else {
    log_print(paste('Alignment coverage threshold:', opt$min_alignment_coverage_thr), hide_notes = TRUE, console = FALSE)
  }
}


if( !is.null(opt$genome_edges_length)) {
  if( opt$genome_edges_length < 10000){
    log_error("The refrennce edges should be at least 10K length")
    log_close()
    writeLines(readLines(lf))
    stop("The refrennce edges should be at least 10K length!")
  } else {
    log_print(paste('refrennce edges length:', opt$genome_edges_length), hide_notes = TRUE, console = FALSE)
  }
}

if(! (opt$telo_index %in% c("telomere", "mismatch", "tvr"))) {
  log_error("The telomere index parameter should be telomere, mismatch or tvr!")
  log_close()
  writeLines(readLines(lf))
  stop("The telomere index parameter should be telomere, mismatch or tvr!")
} else {
  log_print(paste("Calculating the subtelomeric length using", opt$telo_index), , hide_notes = TRUE, console = FALSE)
}

 
# Null defaullts :   
df_join <- join_df(nanotel_summary_path = opt$telo_summary_path, minimap_path = opt$aligner_summary_path) 
log_print(paste("There are", nrow(df_join), "telomeric reads."), hide_notes = TRUE, console = FALSE) 


# calculate sub-telo accordingly
df_join <- calculate_subtelo(df_join, telo_index = opt$telo_index, telo_position = ifelse(opt$telo_right, yes = "right", no = "left"))




log_print("Arguments structure:", console = FALSE, hide_notes = TRUE)
log_print(capture.output(str(opt)), console = FALSE, hide_notes = TRUE)

if(opt$subtelo_length_thr > 0) { 
  df_join  <- df_join  %>% 
    mutate(pass_subtelo_length = subtelo_length >= opt$subtelo_length_thr)
  log_print(paste(sum(df_join$pass_subtelo_length), "reads pass the alignment subtelomeric length filteration of threshold", opt$subtelo_length_th, "!"), console = FALSE, hide_notes = TRUE)
}


# filterations  
df_join <- mapping_filter(df_join, filter = TRUE, filter_column = 'alignment_genome')  

df_join <- mapping_filter(df_join, filter = opt$min_alignment_mapping_quality,
           filter_column = 'alignment_mapping_quality', 
           thr = opt$min_alignment_mapping_quality)  

df_join <- mapping_filter(df_join, filter = opt$filter_genome_position , filter_column = 'alignment_genome_start', thr = opt$filter_genome_position, genome_length = opt$genome_edges_length) 

df_join <- mapping_filter(df_join, filter = opt$min_alignment_accuracy, filter_column = 'alignment_accuracy', thr = opt$min_alignment_accuracy)  

df_join <- mapping_filter(df_join, filter = opt$min_alignment_coverage_thr, filter_column = 'alignment_coverage', thr = opt$min_alignment_coverage_thr, telo_index = opt$telo_index,telo_right =telo_rightt$telo_right )    

df_join <- mapping_filter(df_join, filter = opt$filter_direction, filter_column = 'alignment_direction')   

write_csv(x = df_join, file = paste(opt$save_path, "summary_merged.csv", sep = '/'))

df_join <- df_join %>% 
  mutate(pass_all = (if_all(starts_with("pass"))== TRUE))
df_pass <- df_join %>% dplyr::filter(pass_all == TRUE)
log_print(paste(nrow(df_pass), "reads passed all alignment filterations!"), hide_notes = TRUE, console = FALSE)

create_dirs(path = opt$save_path, alignments = unique(df_pass$alignment_genome))

unclassified_serial <- df_join %>%  dplyr::filter(! (Serial %in% df_pass$Serial)) %>%  pull(Serial)

copy_reads(df_merged = df_join, chrs = unique(df_pass$alignment_genome), 
  nanotel_path = opt$nanotel_path, chrs_path = opt$save_path,file_extension = opt$file_extension, 
  unclassified_serial = unclassified_serial)

copy_plots(df_merged = df_join, chrs = unique(df_pass$alignment_genome), nanotel_path = opt$nanotel_path, chrs_path = opt$save_path, file_extension = '.jpeg', dir_name =  'single_read_plots_adj', unclassified_serial = unclassified_serial)

copy_plots(df_merged = df_join, chrs = unique(df_pass$alignment_genome), nanotel_path = opt$nanotel_path, chrs_path = opt$save_path, file_extension = '.eps', dir_name =  'single_read_plots_adj',  unclassified_serial = unclassified_serial)
copy_plots(df_merged = df_join, chrs = unique(df_pass$alignment_genome), nanotel_path = opt$nanotel_path, chrs_path = opt$save_path, file_extension = '.jpeg', dir_name =  'single_read_plots',  unclassified_serial = unclassified_serial)




log_close(footer = FALSE)
writeLines(readLines(lf))















# log ptint work inside










