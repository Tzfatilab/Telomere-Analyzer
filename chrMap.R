# Given nanotel output + minimap2 dorado summaryfile -> extract filter with a
# dir for the mapped

# according to dorado aligner ( dorado version 1.3.1) +  and NanoTel version v1.1.4-beta

# last updated 2026-01-29

min_mapq <- 60 # add several parameters
min_num_align <- 6000
# filter for each param : if not pass add why

nanotel_path <- "/home/tzfati/Downloads/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN-20260128T161200Z-3-001/barcode09-T83-Output-with-6TVR-CCCCTAGNNNNNNNN"

library(readr)
library(dplyr)
library(stringr)
library(conflicted)

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



dir_path <- "/media/tzfati/Storage/Nanopore_runs/T83/T83_TGBrepeatedT81_8m3wTGB/20260101_1946_P2S-03064-A_PBK31479_c54f8d5a/fastq_pass/barcode09-minimap2"

setwd(dir_path)

df_telo <- read_csv(file = "summary-bc09-Trial83.csv")


df_telo$sequence_ID <- str_sub(string = df_telo$sequence_ID, start = 1, end = 36)



# Process chunks with a function
process_chunk <- function(chunk, pos) {
  filtered <- chunk %>% dplyr::filter(read_id %in% df_telo$sequence_ID)
  return(filtered)
}

df_minimap <- read_tsv_chunked("sequencing_summary.tsv", 
                         DataFrameCallback$new(process_chunk),
                         chunk_size = 10000)


df_minimap <- df_minimap %>% 
  dplyr::select(c(4, 21:39))



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

