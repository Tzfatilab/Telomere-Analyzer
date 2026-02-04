#!/usr/bin/env Rscript


# Given nanotel output + minimap2 dorado summaryfile -> extract filter with a
# dir for the mapped

# according to dorado aligner ( dorado version 1.3.1) +  and NanoTel version v1.1.4-beta

# last updated 2026-02-02


# todo: check if it has tvr_pattern columns....


min_mapq <- 60 # add several parameters
min_num_align <- 6000
# filter for each param : if not pass add why

nanotel_path <- "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN"

library(readr)
library(dplyr)
library(stringr)
library(conflicted)

library(optparse)
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
  res <- logical()
  for(chr in chrs) {
    serial <- df_merged %>% dplyr::filter(alignment_genome == chr) %>% 
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

copy_plots <- function(df_merged, chrs, nanotel_path, chrs_path, 
  file_extension =c('.jpeg', '.eps'), 
  dir_name = c('single_read_plots', 'single_read_plots_adj'), unclassified_serial )
{
  res <- logical()
  for(chr in chrs) {
    serial <- df_merged %>% dplyr::filter(alignment_genome == chr) %>% 
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

create_dirs(path = "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test3", alignments =  unique(df_merged$alignment_genome))

res2 <- copy_plots(df_merged = df_merged, chrs = unique(df_merged$alignment_genome) , nanotel_path = nanotel_path , 
                   chrs_path = "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test3",
                   file_extension = '.jpeg', dir_name =  'single_read_plots', unclassified_serial = unclassified$Serial)

res3 <- copy_plots(df_merged = df_merged, chrs = unique(df_merged$alignment_genome) , nanotel_path = nanotel_path , 
                   chrs_path = "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test3",
                   file_extension = '.jpeg', dir_name =  'single_read_plots_adj', unclassified_serial = unclassified$Serial)

res4 <- copy_plots(df_merged = df_merged, chrs = unique(df_merged$alignment_genome) , nanotel_path = nanotel_path , 
                   chrs_path = "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test3",
                   file_extension = '.eps', dir_name =  'single_read_plots_adj', unclassified_serial = unclassified$Serial)


res <- copy_reads(df_merged = df_merged, chrs = unique(df_merged$alignment_genome), 
           nanotel_path =nanotel_path , 
           chrs_path ="/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test3" 
           , file_extension = ".fasta", unclassified_serial = unclassified$Serial)

# get the df and return the merged one
get_df <- function(nanotel_summary_path, minimap_path)
{
  
}

minimap_path <- "/media/tzfati/Storage/Nanopore_runs/T83/T83_TGBrepeatedT81_8m3wTGB/20260101_1946_P2S-03064-A_PBK31479_c54f8d5a/fastq_pass/barcode09-minimap2/sequencing_summary.tsv"

setwd(dir_path)

nanotel_summary_path <- "~/Desktop/minknow_runs/T83/nanotel_output/barcode9/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/summary.csv"

col_names <- c('read_id', 'sequence_length_template', 'alignment_genome', 
              'alignment_direction', 'alignment_genome_start', 'alignment_genome_end', 'alignment_strand_start', 'alignment_strand_end', 'alignment_num_insertions', 'alignment_num_deletions',
              'alignment_num_aligned','alignment_num_correct', 'alignment_identity', 'alignment_accuracy', 'alignment_score', 
              'alignment_coverage',  'alignment_mapping_quality', 'alignment_num_alignments', 'alignment_num_secondary_alignments')                      

col_types <- "ciccfiiiiiliiiidffffffiiiiiiiiddidiiiii"
                                                                                               
df_minimap <- read_tsv(file = minimap_path, n_max = 100000, 
             col_select = col_names, col_types = col_types)
             
            
df_minimap_all <- read_tsv(file = minimap_path, n_max = 100000)
                                 




df_telo <- read_csv(file = nanotel_summary_path, col_types = "icidiiidiiidiii")

# change to new col : read_id
df_telo$read_id <- str_sub(string = df_telo$sequence_ID, start = 1, end = 36)

minimap_path <- "~/Desktop/minknow_runs/T83/sequencing_summary.tsv"

# Process chunks with a function
process_chunk <- function(chunk, pos) {
  filtered <- chunk %>% dplyr::filter(read_id %in% df_telo$read_id)
  return(filtered)
}

df_minimap <- read_tsv_chunked(minimap_path, 
                         DataFrameCallback$new(process_chunk),
                         chunk_size = 10000, col_types = col_types)


df_minimap <- df_minimap %>% 
  dplyr::select(c(4, 21:39))

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


df_join <- join_df(nanotel_summary_path = nanotel_summary_path, minimap_path = minimap_path)
# add 

df_join2 <- join_df(nanotel_summary_path = nanotel_summary_path, minimap_path = minimap_path)

# add summary to log file ....


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

filter_alignment <- function(df_join, )


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
  make_option("--filter_genome_position", action = "store", default = -1,
              help = "Filter alignments according to alignment_genome_start/genome_end, with threshold of given integer >= 0:.", 
              type = 'integer',metavar = "USE filter_direction filter with a given threshold" ), 
  
  
  make_option("--min_alignment_accuracy", action = "store", default = NULL, 
              type = "double", 
              help = "Minimal alignment accuracy (alignment_accuracy = (number of matches) / (matches + mismatches + insertions + deletions) × 100).",
              metavar = "Minimal alignment accuracy score."), 
  
  make_option("--min_alignment_coverage_thr", action = "store", default = NULL, 
              type = "double", 
              help = "Minimal threshold for The subtelomere coverage: abs(alignment_coverage - subtelomere_length/read_length) <= threshold.",
              metavar = "Alignmnet coverage threshold"), 
  
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
              help = "Print version information and exit")
  
)   

opt = parse_args(OptionParser(option_list=option_list))  
  


# Handle --version flag
if (opt$version) {
  cat("Telomere Analyzer  version v1.1.4-beta\n")
  quit(save = "no", status = 0)
}  

if( !is.null(opt$min_alignment_mapping_quality)) {
  if( (opt$min_alignment_mapping_quality < 0) ||
      (opt$min_alignment_mapping_quality > 60)
  ){
    stop("The alignment mapping quality threshold should be an integer in [0,60]!")
  }
}


if( !is.null(opt$min_alignment_accuracy)) {
  if( (opt$min_alignment_accuracy < 0) ||
      (opt$min_alignment_accuracy > 1)
  ){
    stop("The alignment accuracy threshold should be a float in [0,1]!")
  }
}

 
df_join <- join_df(nanotel_summary_path = opt$telo_summary_path, minimap_path = opt$aligner_summary_path)  
  
# filterations  
  
  
  
  
# optional filterations:
#' 1.  alignment_genome != '*' This is must filter 
#' 2.  alignment_direction     v
#' 3.  alignment_genome_start  v
#' 4.  alignment_genome_end    v
#' 5.  alignment_accuracy      v [0,1]
#' 6.  alignment_coverage :    v ( claculate sub-telo coverage and see if it is ~ alignment_coverage(alignment_coverage = (aligned length of query) / (total query length) × 100)
#' 7.  alignment_mapping_quality v [0,60] int   
mapping_filter <- function(df_join, filter_column='alignment_genome' , thr=NULL, genome_length=NULL)
{
  if(filter_column == 'alignment_genome') {
    df_join <- df_join %>% 
      mutate(pass_alignment_genome = (alignment_genome != '*') )
    
    return(df_join)
  }
  
  # We assume heads are + and atails are -
  if(filter_column == 'alignment_direction') {
    df_join <- df_join %>% 
      mutate(pass_alignment_direction = 
      ( (str_detect(string = alignment_genome, pattern = "Head") & alignment_direction == "+" ) | 
      (str_detect(string = alignment_genome, pattern = "Tail") & alignment_direction == "-" ) ) )
    
    return(df_join)
  }
  
  
  # This is to be detected!
  if(filter_column == 'alignment_genome_start') {
    df_join <- df_join %>% 
      mutate(pass_alignment_genome_start_end = (alignment_genome_start <= thr & str_detect(string = alignment_genome, pattern = "Head")) | 
               (abs(alignment_genome_end - genome_length) <= thr & str_detect(string = alignment_genome, pattern = "Tail")) ) 
    
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
    
    return(df_join)
  }
  
  # todo: Need toa calculate the sub-telo length
  if(filter_column == 'alignment_coverage') {
    df_join <- df_join %>% 
      mutate( )
    
    return(df_join)
  }
  
  if(filter_column == 'alignment_mapping_quality') {
    df_join <- df_join %>% 
      mutate(pass_alignment_mapping_quality = (alignment_mapping_quality  >= thr))
    
    return(df_join)
  }
}
  
  
# test alignment_genome filter
df_pass_genome <- df_join2 %>% 
  mutate(pass_alignment_mapping_quality = (alignment_mapping_quality  >= 50))

  
# test map quality 

# test direction
df_pass_genome <- df_join2 %>% 
  mutate(pass_alignment_direction =  ( (str_detect(string = alignment_genome, pattern = "Head") & alignment_direction == "+" ) | 
                 (str_detect(string = alignment_genome, pattern = "Tail") & alignment_direction == "-" ) ) )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# factorize
df_minimap$alignment_genome <- base::as.factor(df_minimap$alignment_genome)

# filter 
df_minimap_f <- df_minimap %>% 
  dplyr::filter(alignment_mapping_quality >= min_mapq) %>% 
  dplyr::filter( (str_detect(string = alignment_genome, pattern = "Head") & alignment_direction == "+" ) | 
            (str_detect(string = alignment_genome, pattern = "Tail") & alignment_direction == "-" ) ) %>% 
  dplyr::filter(alignment_num_aligned >= min_num_align)

unclassified <- df_telo %>% 
  dplyr::filter(!(sequence_ID %in% df_minimap_f$read_id))

# now merge with the telo df and make dirs and copy files accordingly

# merge 
df_merged <- left_join(x = df_minimap_f , y = df_telo, by = join_by(read_id==sequence_ID))
chr <- 'chromosome_25_Tail' 





df_chr <- df_merged %>% dplyr::filter(alignment_genome == chr)



dir_path <- nanotel_path

# make dir
base::dir.create(path = str_c(dir_path, "/", chr) )

# make reads dir
base::dir.create(path = str_c(dir_path, "/", chr, "/reads") )

from_path <- get_reads_path(serial = df_chr$Serial, reads_path = dir_path)
to_path <- str_replace_all(string = from_path, pattern = dir_path, replacement = str_c(dir_path, "/", chr))

file.copy(from = from_path  ,to = to_path)


# test copy_read function 
test_result <- copy_reads(df_merged = df_merged, chr = chr, nanotel_path = nanotel_path, chrs_path = "~/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN/test1", file_extension = ".fasta")

