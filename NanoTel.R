#!/usr/bin/env Rscript

################################################################
# Author: Dan Lichtental
# Copyright (c) Dan Lichtental, 2022
# Email:  dan.lichtental@mail.huji.ac.il
# Last updated:  2024-Jun-19
# Script Name: Telomere pattern finder
# Script Description: In this script we search over fastq file sequences for
# telomeric patterns.
#
# Notes: Running the script on linux:
# "Rscript --vanilla NanoTel -i input_path --save_path output_path --patterns "some patterns"
# See Rscript --vanilla NanoTel.R  --help for more additional flags
# format can be either fasta or fastq, input_path can be either a full path for a
# single fastq/a file or a directory containing fastq/a files.
# Please try to give the output_path outside the input dir path and vice versa

suppressPackageStartupMessages(require(optparse))

########## start of the flags ###############
# global vars
# My last change: 22/10/2022 - change the global_min_density from 0.3 to 0.6
global_min_density <- 0.6

global_subseq_length <- 100



option_list = list(
  
  # input_path
  make_option(c("-i","--input_path"), action = "store", type = "character" , 
              default = NULL, help = "Path to input files.( dir or single file)", metavar = "input path"),
  
  # save_path
  make_option("--save_path", action = "store", 
              type = "character" , default = NULL, 
              help = "A path to a directory for storing the output files.", 
              metavar = "output dir path"),
  
  # format( fastq/fasta)
  make_option("--format", action = "store", 
              type = "character", default = "fastq" , 
              help = "input files format (Either \"fastq\" (the default) or \"fasta\", gzip is supported)", 
              metavar = "file format"),
  
  # nrec
  make_option(c("-n","--nrec"), action = "store", type = "integer", 
              default = 10000 , 
              help = "Single integer. The maximum of number of records to read in to memory for each iteration. Negative values are ignored."),
  
  # do rc ?
  make_option(c("-r","--rc"), action = "store_true", default = FALSE, 
              help = "Should we do reverse complement on the given reads.", 
              metavar = "Reverse complement on reads."),
  
  # pattern/s
  make_option("--patterns", action = "store", default = NULL, type = "character", 
              help = "Space separated list of pattern(s). Must be in double quotes.", 
              metavar = "pattern"), 
  
  make_option("--min_density", action = "store", default = 0.6, 
              type = "double", 
              help = "Minimal density to consider a subsequence as a pattern region.",
              metavar = "Minimal density."), 
  
  make_option("--subseq_length", action = "store", default = 100, 
              type = "integer", help = "The length of the sub-sequence.",
              metavar = "Sub-sequence length." ), 
  
  make_option("--use_filter", action = "store_true", default = FALSE,
              help = "Filter reads accoding to the edge.", 
              metavar = "USE FILTER" ), 
  
  make_option("--check_right_edge", action = "store_true", default = FALSE, 
              help = "This flag tells us the expected telomere position: helps with the filter and telo_position accuracy", 
              metavar = "Check right or left edge for filter and position"), 
  make_option("--tvr_patterns", action = "store", default = NULL, type = "character", 
              help = "Space separated list of additional pattern(s) for Telomere variant repeats. Must be in double quotes.", 
              metavar = " Telomere variant repeats patterns"), 
  
  make_option("--version", 
              action = "store_true", 
              default = FALSE,
              help = "Print version information and exit")
  
)   
opt = parse_args(OptionParser(option_list=option_list))

# Handle --version flag
if (opt$version) {
  cat("Telomere Analyzer  version v1.1.7-beta 2026-02-18 \n")
  quit(save = "no", status = 0)
}  



# The 'confilcted' package tries to make your function choice explicit.
#   it produce an error if a function name is found on 2 or more packages!!!
suppressPackageStartupMessages(require(S4Vectors))
suppressPackageStartupMessages(require(BiocGenerics))
suppressPackageStartupMessages(require(conflicted))
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(logr))
suppressPackageStartupMessages(require(future))
suppressPackageStartupMessages(require(ggprism))

suppressPackageStartupMessages(require(testit))

#utils::globalVariables(c("start_index"))




#' my changes: 5.11.2023
#' 1. Change thr for re-indexing  telomere_density < 0.85 
#'  (telomere_density2 < 0.85) this should exclude small islands from the telo.
#' 2. Added: print the input arguments to the log file.




# test matchPattern with latters !in c("A", "C", "G", "T")








#' 3. Add option to mismatch, max_mismatch (default == 1), or fixed list of mismatches 





#' TODO: replace all for loops to map with a private functio
#' instead of for(...{if...else...} : map(x, p_function)
#' map(numlist, ~.x %>% sqrt() %>% sin())

# Last change: 16:09 4/10/2023 max-matches <- 2 (insteadof 1)

#










global_min_density <- 0.5
global_subseq_length <- 100


#' bug fix of faild calculating the telomere indcies after re-telo_position:
#' Check: My last change : 22/10/2023



# try using this code: 
#' # Install and load the GenomicRanges package

#library(GenomicRanges)
# Set the total length and width
#total_length <- 1000
#width <- 100
# Create a GRanges object
#ranges <- GRanges(seqnames = 'chr1', ranges = IRanges(start = seq(1, total_length, width), end = seq(width, total_length, width)))
# Print the GRanges object
#print(ranges)

split_telo <- function(dna_seq, sub_length) {
  #' @title Splits a DNA sequence into subsequences.
  #' @description This function calculate the sequence len and creates IRanges
  #' objects of subseuences of a given length.
  #' if The dna_seq%%subseq != 0   and last' width < sub_length/2 Then we will
  #' remove this last index making the last subtelomere longer then by.
  #' (because we want to prevent a case where we have for
  #' example subtelomere of length < sub_length/2 which is too short to consider
  #'  -> then the last subtelomere is a bit longer)
  #' If the length of dna_seq is less then the sub_length it will return an
  #' empty Iranges Object.
  #' @usage
  #' @param dna_seq: DNAString object
  #' @param sub_length: The length of each subsequence
  #' @value An IRanges object of the subsequences.
  #' @return iranges_idx: IRanges obsject of the indices of each subsequence
  #' @examples
  idx_start <- seq(1, length(dna_seq), by = sub_length)
  idx_end <- idx_start + sub_length - 1
  idx_end[length(idx_end)] <- length(dna_seq)
  # if last subsequence is less then 50%
  if (length(dna_seq) - dplyr::last(idx_start) < (sub_length / 2)) {
    idx_start <- idx_start[1:length(idx_start) - 1]
    idx_end <- idx_end[1:length(idx_end) - 1]
    idx_end[length(idx_end)] <- length(dna_seq)
  }
  idx_ranges <- IRanges(start = idx_start, end = idx_end)
  return(idx_ranges)
} # IRanges object


# create a patterns which are create from 1 pattern using 1 
# mismatch + the original pattern
get_one_mismatch <- function(seq, abc = c("A", "C", "G", "T")) {
  seq_list <- list()
  seq_length <- str_length(seq)
  
  if(seq_length == 1) {
    return(as.list(abc))
  }
  
  if(seq_length == 0) {
    return(list())
  }
  
  for(i in seq_len(str_length(seq))) {
    for(l in abc ) {
      if(i == 1)
        curr_seq <- str_c(l, str_sub(seq, start = 2, end = seq_length))
      else if (i == seq_length)
        curr_seq <- str_c(str_sub(seq, start = 1, end = seq_length - 1), l)
      else
        curr_seq <- str_c(str_sub(seq, start = 1 , end = i -1), l , str_sub(seq, start = i + 1, end = seq_length))
      seq_list <- base::append(seq_list, curr_seq)
    }
  }
  
  return(seq_list)
}








# My last change: 15/10/2023 - use trim for out-of-bound matches"
#' This is a feature. The matchPattern/vmatchPattern family of functions treat 
#' out-of-bound nucleotide positions as mismatches (think of these positions as 
#' having letters that don't match any letter in DNA_ALPHABET).
#' 
#' Note that the out-of-bound matches are easy to remove:
#' 
#' library(Biostrings)
#' matches <- matchPattern("ATGG", subject, max.mismatch=1)
#' matches 

# Views on a 16-letter DNAString subject
# subject: AATGCGCGTGGATATG
# views:
#       start end width
#   [1]     2   5     4 [ATGC]
#   [2]     8  11     4 [GTGG]
#   [3]    14  17     4 [ATG ]
#' 
## Remove the out-of-bound matches:
#' keep_me <- (0 <= start(matches)) & (end(matches) <= length(subject(matches)))
#' matches[keep_me]
# Views on a 16-letter DNAString subject
# subject: AATGCGCGTGGATATG
# views:
#       start end width
#   [1]     2   5     4 [ATGC]
#   [2]     8  11     4 [GTGG]
#' Alternatively, you can trim the out-of-bound matches:
#' trim(matches)
# Views on a 16-letter DNAString subject
# subject: AATGCGCGTGGATATG
# views:
#       start end width
#   [1]     2   5     4 [ATGC]
#   [2]     8  11     4 [GTGG]
#   [3]    14  16     3 [ATG]

# my improvment: fiding the IRanges and making union for overlaps and then
#calculate according to sum(width of the IRanges)
# The calculation is on the full sequence and it is not fit for subsequences
#' use purrr::map to replace the for loops code
get_density_iranges <- function(sequence, patterns, with_mismatch = FALSE, tvr_patterns = NULL) {
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a
  #' list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
  #' @param tvr_patterns: Additional TVR's.(Telomeric varient repeats)
  #' @value A numeric for the total density of the pattern(patterns) in the
  #' sequences, a IRanges object of the indcies of the patterns found.
  #' @return a tuple of (density, IRanges) Total density of all the patterns in
  #' the list( % of the patterns in this sequence) and the IRanges of them
  #' @examples
  total_density <- 0
  max_mismatch <- 0
  if(with_mismatch) {
    max_mismatch <- 1
  }
  
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list
  if (is.list(patterns)) {
    patterns <- unique(patterns)  # make sure there are no dups
    
    #  use map instead and do 
    for (pat in patterns){
      #' letters that represent ambiguity which are used when more than one kind
      #'  of nucleotide could occur at that position - than it is not fixed pattern
      fixed <- !str_detect(string = pat, pattern = "[WSMKRYBDHVN]")
      curr_mp <- matchPattern(pattern = pat, subject =
              unlist(sequence),  fixed = fixed,  max.mismatch = max_mismatch)
      if( (fixed == FALSE) || (max_mismatch > 0) ) {
        curr_mp <- trim(curr_mp)
      }
      
      mp_all <- IRanges::union(mp_all, curr_mp)
      
      
    }
    mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
    
  } else {
    fixed <- !str_detect(string = patterns, pattern = "[WSMKRYBDHVN]")
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence),
                           max.mismatch = max_mismatch, fixed = fixed)
    if( (fixed == FALSE) || (max_mismatch > 0) ) {
      mp_all <- trim(mp_all)
      mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
    }
    
  }
  
  
  # additional TVR's
  if(!is.null(tvr_patterns)) {
    if (is.list(tvr_patterns)) {
      tvr_patterns <- unique(tvr_patterns)  # make sure there are no dups
      
      #  use map instead and do 
      for (pat in tvr_patterns){
        #' letters that represent ambiguity which are used when more than one kind
        #'  of nucleotide could occur at that position - than it is not fixed pattern
        fixed <- !str_detect(string = pat, pattern = "[WSMKRYBDHVN]")
        curr_mp <- matchPattern(pattern = pat, subject =
                                  unlist(sequence),  fixed = fixed )
        if( (fixed == FALSE) || (max_mismatch > 0) ) {
          curr_mp <- trim(curr_mp)
        }
        
        mp_all <- IRanges::union(mp_all, curr_mp)
      }
      
      
      
    mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
    
      
    }else {
      fixed <- !str_detect(string = tvr_patterns, pattern = "[WSMKRYBDHVN]")
      mp_curr <- matchPattern(pattern = tvr_patterns, subject = unlist(sequence),
                              fixed = fixed)
      if( (fixed == FALSE) || (max_mismatch > 0) ) {
        mp_curr <- trim(mp_curr)
        mp_all <- IRanges::union(mp_all, mp_curr) # incase there are overlaps
      }
      mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
    }
  }
  
  total_density <- sum(width(mp_all)) / nchar(sequence)
  return(list(total_density, mp_all))
}


# with a data frame we can further explore the different patterns
# my improvment: fiding the IRanges and making union for overlaps and then
# calculate according to sum(width of the IRanges)
# The calculation is on the full sequence and it is not fit for subsequences
#' use purrr::map to replace the for loops code
get_density_iranges_with_csv <- function(sequence, patterns, output_path = NA) {
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a
  #' list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
  #' @param output_path: the filename for the csv patterns summary.
  #' @value A numeric for the total density of the pattern(patterns) in the
  #' sequences, a IRanges object of the indcies of the patterns found.
  #' @return a list of (density, IRanges, data frame) Total density of all the
  #' patterns in the list( % of the patterns in this sequence) and the IRanges
  #' of them
  #' @examples
  total_density <- 0
  patterns_df <- data_frame(patterns = character(), start_idx = integer(),
                            end_idx = integer())
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list
  if (is.list(patterns)) {
    patterns <- unique(patterns)  # make sure there are no dups
    for (pat in patterns) {
      curr_match <- matchPattern(pattern = pat, subject = unlist(sequence),
                                 max.mismatch = 0, fixed = FALSE)
      patterns_df <- add_row(patterns_df, patterns = as.data.frame(curr_match)$x
                     , start_idx = start(curr_match), end_idx = end(curr_match))
      mp_all <- IRanges::union(mp_all, curr_match)
    }
  }else {
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence),
                           max.mismatch = 0, fixed = FALSE)
    patterns_df <- add_row(patterns_df, patterns = as.data.frame(mp_all)$x,
                           start_idx = start(mp_all), end_idx = end(mp_all))
    mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
  }
  total_density <- sum(width(mp_all)) / nchar(sequence)
  pat_list <- list(total_density, mp_all, patterns_df)
  names(pat_list) <- c("total density", "patterns IRanges",
                       "Patterns data frame")
  if (!isTRUE(is.na(output_path))) {
    write_csv(x = patterns_df, file = output_path)
  }
  return(pat_list)
}


get_sub_density <- function(sub_irange, ranges) {
  #' @title Calculate density of a subsequnce.
  #' @details with a given IRanges of a subsequence and IRanges of patterns,
  #' compute the density of the IRanges within the
  #'        subseuence range.
  #' @param sub_irange: the IRange of the subsequence
  #' @param ranges: The IRanges of the patterns found in the full sequence
  #' @value a numeric which is the density in range [0,1].
  #' @return The density of the patterns in the subseuence according to the
  #' IRanges.
  #' @description  sub_irange = (10, 30), ranges = {(2,8), (16,21), (29,56)} ->
  #' the intersect is {(16,21), (29, 30)} ->
  #                width = 6+2 = 8 , sub_irange width = 21 -> density = 8/21
  #' @examples sub_irange <- IRanges(start = 10, end = 30)
  #'           ranges <- IRanges(start = c(2,16,29), end = c(8,21,56))
  #'           get_sub_density(sub_irange =  sub_irange, ranges = ranges) # 0.38
  #' this compute the desity of iranges of patterns within a given irange of a
  #' subsequence
  return(sum(width(IRanges::intersect(sub_irange, ranges))) / width(sub_irange))
}















# I Can merge the 2 functions to 1 with additional paramater which tells me to return start or end index?
# ' Helper for search_left_patterns
#' I need to think about a pattern which is a list of patterns
#' -> return an start index ( maximal end index)
#' We assume patterns is a list of strings.
#' @param read - A dna seq DNAString
#' @param patterns - A list of Strings , each represen a pattern
#' @param subseq_start - Index for the start of the subsequence we want to check.
#' @param subseq_end - "               end ".
#' @param with_mismatches - If True, allow 1 mismatch.
#' @tvr_patterns - Additional patterns to check, with no mismatch for them.
#' @returns  - New index, or Inf if there is no updated start index.
multi_pattern_step_left <- function(read,patterns, subseq_start, subseq_end, with_mismatches = FALSE, tvr_patterns = NULL) {
  new_start <- Inf
  
  if( (is.null(tvr_patterns) || with_mismatches) == FALSE ) { # FF
    all_patterns <- unique(unlist(list(patterns, tvr_patterns)))
    for(pat in all_patterns) {
      curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end ) )
      if(length(curr_mp) > 0) {
        new_start <- min(new_start, min(start(curr_mp)))
      }
    }
    
    return(ifelse(new_start == Inf, new_start , as.integer(new_start + subseq_start - 1)) )
  }
  
  for(pat in patterns) {
    curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end), max.mismatch = as.numeric(with_mismatches) ) 
    if(length(curr_mp) > 0) {
      new_start <- min(new_start, min(start(curr_mp)))
    }
  }
  if(!is.null(tvr_patterns)) {
    for(pat in tvr_patterns) {
      curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end) )
      if(length(curr_mp) > 0) {
        new_start <- min(new_start, min(start(curr_mp)))
      }
    }
    
  }  
  
  return(ifelse(new_start == Inf, new_start , as.integer(new_start + subseq_start - 1)) )
}



######## Helper functions for more accurate telomere length ##############


#' I need to think about a pattern which is a list of patterns
#' -> return an end index ( maximal end index)
#' We assume patterns is a list of strings.
#' @param read - A dna seq DNAString
#' @param patterns - A list of Strings , each represen a pattern
#' @param subseq_start - Index for the start of the subsequence we want to check.
#' @param subseq_end - "               end ".
#' @param with_mismatches - If True, allow 1 mismatch.
#' @tvr_patterns - Additional patterns to check, with no mismatch for them.
multi_pattern_step_right <- function(read,patterns, subseq_start, subseq_end, with_mismatches = FALSE, tvr_patterns = NULL) {
  new_end <- -1
  if( (is.null(tvr_patterns) || with_mismatches) == FALSE ) { # FF
    all_patterns <- unique(unlist(list(patterns, tvr_patterns)))
    for(pat in all_patterns) {
      curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end ) )
      if(length(curr_mp) > 0) {
        new_end <- max(new_end,max(end(curr_mp)))
      }
    }
    
    return(ifelse(new_end == -1, new_end, new_end + subseq_start - 1) )
  }
  
  for(pat in patterns) {
    curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end), max.mismatch = as.numeric(with_mismatches) ) 
    if(length(curr_mp) > 0) {
      new_end <- max(new_end,max(end(curr_mp)))
    }
  }
  if(!is.null(tvr_patterns)) {
    for(pat in tvr_patterns) {
      curr_mp <- matchPattern(pattern = pat, subject = subseq(read, start = subseq_start, end = subseq_end) )
      if(length(curr_mp) > 0) {
        new_end <- max(new_end,max(end(curr_mp)))
      }
    }
    
  }  
  
  return(ifelse(new_end == -1, new_end, new_end + subseq_start - 1) )
}
search_left_patterns <- function(read, start_index, subseq_width = 18, step_size = 10, max_steps = 4, pattern , with_mismatch = FALSE, tvr_patterns = NULL) {
  #' This function advance to the left to find more patterns for more accurate 
  #' start index.
  #' @param read - A dna seq
  #' @param start_index - the index to start from
  #' @param subseq_width - The length of the subseq each time we search.
  #' @param step_size - The step for the subseq each time we move more to the left.
  #' @param max_steps - How many iterations.
  #' @param  pattern - The dna pattern to serach (a single string)
  #' @param tvr_patterns - As with right...
  #' @param with_mismatch - If TRUE allow 1 mismatch.
  
  
  assert("Subsequence lengt is not smaller than the Pattern length!", str_length(pattern) <= subseq_width)
  #  I need to set boundary checks: length(pattern) < subseq_width ....
  subseq_start <- max(start_index - subseq_width, 1)
  new_start <- start_index
  for( i in 1:max_steps) {
    curr_end <- min(subseq_start + subseq_width - 1, length(read)) # make sure no out of boundry
    
    if(is.list(pattern) ) {
      # max for incase returns Inf
      curr_start <- multi_pattern_step_left(read = read, patterns = pattern, 
                                            subseq_start = subseq_start, subseq_end = curr_end, 
                                            with_mismatches = as.numeric(with_mismatch), 
                                            tvr_patterns = tvr_patterns )
      if(curr_start == Inf) { break }
      new_start <- curr_start
      
    } else {
      if(!is.null(tvr_patterns)) {
        curr_start <- multi_pattern_step_left(read = read, patterns = as.list(pattern), 
                                              subseq_start = subseq_start, subseq_end = curr_end, 
                                              with_mismatches = with_mismatch, tvr_patterns = tvr_patterns)
        if(curr_start == Inf) { break }
        new_start <- curr_start
        
      } else { #No list, no TVR's
        curr_mp <- matchPattern(pattern = pattern, 
                                subject = subseq(read, start = subseq_start, end = curr_end), 
                                max.mismatch = as.numeric(with_mismatch) )
        if(length(curr_mp) == 0) { 
          break
        }
        new_start <- min(start(curr_mp)) + subseq_start - 1
        
      } # End of No list, no TVR's
    }
    
    
    # iteration part
    subseq_start_new <- max(subseq_start - step_size +1, 1) # check boundry
    if(subseq_start_new == subseq_start) {break} # we already checked it!
    subseq_start <- subseq_start_new
  }
  return(new_start)
  
}

search_right_patterns <- function(read, end_index, subseq_width = 18, step_size = 10, max_steps = 4, pattern , with_mismatch = FALSE, tvr_patterns = NULL) {
  #' This function advance to the left to find more patterns for more accurate 
  #' end index.
  #' @param read - A dna seq DNAString
  #' @param end_index - the index to the end to start from.
  #' @param subseq_width - The length of the subseq each time we search.
  #' @param step_size - The step for the subseq each time we move more to the left.
  #' @param max_steps - How many iterations.
  #' @param  pattern - The dna pattern to serach (a single string)
  #' @param with_mismatch - If TRUE allow 1 mismatch.
  
  
  assert("Subsequence length is not smaller than the Pattern length!", str_length(pattern) <= subseq_width)
  # avvoid of limit (end_index > length(read))
  subseq_end <- min(end_index + subseq_width, length(read))
  new_end <- end_index
  for( i in 1:max_steps) {
    # avoid of limit start index (<= 1)
    curr_start <- max(subseq_end - subseq_width + 1, 1)
    
    if(is.list(pattern) ) {
      # max for incase returns -1
      curr_end <- multi_pattern_step_right(read = read, patterns = pattern, 
                                           subseq_start = curr_start, subseq_end = subseq_end, 
                                           with_mismatches = with_mismatch, tvr_patterns = tvr_patterns)
      if(curr_end == -1) { break}
      new_end <- curr_end
      
      
      
      
    } else{
      if(!is.null(tvr_patterns)) {
        curr_end <- multi_pattern_step_right(read = read, patterns = as.list(pattern),
                                             subseq_start = curr_start, subseq_end = subseq_end, 
                                             with_mismatches = with_mismatch, tvr_patterns = tvr_patterns) 
        if(curr_end == -1) { break}
        new_end <- curr_end
        
        
      } else {
        curr_mp <- matchPattern(pattern = pattern, 
                                subject = subseq(read, start = curr_start, end = subseq_end ), 
                                max.mismatch = as.numeric(with_mismatch) )
        
        if(length(curr_mp) == 0) { 
          break
        }
        new_end <- max(end(curr_mp)) + curr_start - 1 #  adjust index to the full length read
        
        
        
      } 
    }
    # iteration part
    # avvoid of limit (end_index > length(read))
    subseq_end_new <- min(subseq_end + step_size +1, length(read))
    if(subseq_end_new == subseq_end) { break}
    subseq_end <- subseq_end_new
  }
  return(new_end)
  
}



######## End of Helper functions for more accurate telomere length ##############






#' future changes:
#'   1. remove the ID column
#'   2. class: give the class a name according to type "CCCTAA/TTAGGG"
#'   3. Add the middle_point as in the plotTelomere.R
#'   4. use purrr::map / future_map to replace the for loops code
#'   5. Give names to the elements in the list : for example:
#'      m_list[["subtelos]] <- subtelos
#'      ( see the course Foundations of Functional Programming with purrr)
###########3 My cahnge from prev - return a list 0f (df, total_density) ######
analyze_subtelos <- function(dna_seq, patterns, sub_length,
                       min_density, with_mismatch = FALSE, tvr_patterns = NULL) {
  # return list(subtelos, list_density_mp[1])
  #' @title Analyze the patterns for each subsequence.
  #' @description s split a dna sequence to subsequences and calculate the
  #' density of each subsequence
  #' @param dna_seq: a dna sequence (DNAString object)
  #' @param patterns: a list of patterns or a string of 1 pattern
  #' @param sub_length: The length of the subsequences for split_telo fuction.s
  #' @value a list(data frame of all subtelomeres and their properties ,numeric
  #' for total density)
  #' @return  a list of (a data frame, list(numeric: total density, IRanges for
  #' patterns)
  #' @examples
  # aother option is to create 5 vector s and then make a data.table from them

  # density and iranges of matchPattern
  list_density_mp <- get_density_iranges(dna_seq, patterns = patterns, with_mismatch = with_mismatch, tvr_patterns)
  mp_iranges <- list_density_mp[[2]]
  # get start indexes of "subtelo"
  idx_iranges <- split_telo(dna_seq, sub_length = sub_length)
  # create empty dataframe which will contain all subtelomeres and their
  # properties.
  subtelos <- data.frame(ID = as.integer(), start_index = as.integer(),
                         end_index = as.integer(), density = as.numeric(),
                         class = as.numeric())
  cur_id <- 1             # intitialize ID counter
  # loop through start indexs
  for (i in seq_along(idx_iranges)) {
    subtelo_density <- get_sub_density(sub_irange = idx_iranges[i],
                                      ranges = mp_iranges)
    #TODO : rename -5,1,0 to a facators
    classes <- list("CCCTAA" = -5, "NONE" = 1, "SKIP" = 0)
    ###############3 NEED TO CHANGE FOR AN ARGUMENT OF classes #################
    subtelo_class <- classes$CCCTAA
    if (subtelo_density < min_density) {
      if (subtelo_density < 0.1) {
        subtelo_class <- classes$SKIP
      }else {
        subtelo_class <- classes$NONE
      }
    }
    subtelos <- subtelos %>%
      add_row(ID = cur_id, start_index = start(idx_iranges[i]), end_index =
                end(idx_iranges[i]), density = subtelo_density, class =
                subtelo_class)
    cur_id <- cur_id + 1
  }
  return(list(subtelos, list_density_mp))
}
#  a a list of (a data frame, numeric: total density)

# find a hole inside a telomeric sequence
find_inner_hole <- function(subtelos, min_in_a_row = 3, max_density_score = 0.75) {
    classes <- list("CCCTAA" = -5, "NONE" = 1, "SKIP" = 0)
  # set score, start, in.a.row to 0,-1,0

  score <- 0.0
  start <- -1
  end <- -1
  in_a_row <- 0
  start_end_diff <- subtelos[1, "end_index"] - subtelos[1, "start_index"]

  # loop through subsequences
  end_position <- 0 # for end loop
  for (i in seq_len(nrow(subtelos))){
    subt <- subtelos[i, ]
    # if the subsequence's class is telomeric pattern reset values
    if (subt$class == classes$CCCTAA ) {
      score <- 0
      start <- -1
      in_a_row <- 0
    }else {
      # otherwise add one to in.a.row, update score and set start index
      in_a_row <- in_a_row + 1
      score <- score + subt$density
      if (start == -1) {
        start <- subt$start_index
      }
    }
    if (in_a_row >= min_in_a_row && score <= max_density_score) {
      end_position <- i + 1
      break
    }
  }
  if (end_position == 0) {
    return(IRanges(1, 1)) # no hole was found
  }


  #' search for end from the last subsequence (backward to finding stat incase
  #' there is island of non-telomeric subsequence)
  end <- -1
  score <- 0.0
  in_a_row <- 0
  for (i in nrow(subtelos):end_position){
    subt <- subtelos[i, ]
    if (subt$class ==  subt$CCCTAA) {
      score <- 0.0
      end <- -1
      in_a_row <- 0
    }else {
      # otherwise add one to in.a.row, update score and set start index
      in_a_row <- in_a_row + 1
      score <- score + subt$density
      if (end == -1) {
        end <- subt$end_index
      }
    }
    #' if more than MIN.IN.A.ROW subtelomeres were found and the overall score
    #' is high enough, return start index
    if (in_a_row >= min_in_a_row && score >= max_density_score) {
      break
    }
  }

  if (start > end) {
    end <- start + start_end_diff
  }

  return(IRanges(start = start, end = end))
}



# max_diff what is the max distance between the end of the read to the end of the telomere
find_right_telo <- function(seq_length, subtelos, max_diff = 200) {
  
  classes <- list("CCCTAA" = -5, "NONE" = 1, "SKIP" = 0)
  # set score, start, in.a.row to 0,-1,0
  
  
  # check also if there is the last were it is bigger.smaller
  start_end_diff <- subtelos[1, "end_index"] - subtelos[1, "start_index"]
  start <- 1
  end <- 1
  # loop through subsequences
  # I need to change it to fit the case when last subseq is in range of [101:150] 
  # and  the density is 36/122 < 0.3
  
  last_i <- 1
  
  for (i in nrow(subtelos):1){
    subt <- subtelos[i, ]
    if(subt$end_index < seq_length -max_diff) {
      return(IRanges(-1, -1)) # no telomere was found, change to NULL
    }
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == classes$SKIP || subt$class == classes$NONE ||
        is.na(subt$class)) {
      next
      
    } else {
      # otherwise add one to in.a.row, update score and set start index
      end <- subt$end_index
      last_i <- i
      break
    }
  }
  
  for(i in last_i:1) {
    subt <- subtelos[i, ]
    if (subt$class == classes$SKIP || subt$class == classes$NONE ||
        is.na(subt$class)) {
      break
      
    } else {
      # otherwise add one to in.a.row, update score and set start index
      start <- subt$start_index
      last_i <- i
    }
  }
  
  
  # check also if there is the last were it is bigger.smaller
  start_end_diff <- subtelos[last_i, "end_index"] - subtelos[last_i, "start_index"]
  if (start > end) {
    end <- start + start_end_diff
  }
  
  return(IRanges(start = start, end = end))
  
}  # end of function 'find_right_telo' 





# max_diff what is the max distance between the startof the read to the start of the telomere
find_left_telo <- function(seq_length, subtelos, max_diff = 200) {
  
  classes <- list("CCCTAA" = -5, "NONE" = 1, "SKIP" = 0)
  # set score, start, in.a.row to 0,-1,0
  
  
  
  start <- 1
  end <- 1
  
  last_i <- 1
  
  # start searching the telo patterns from the start (left)
  for (i in seq_len(nrow(subtelos))){
    subt <- subtelos[i, ]
    if(subt$start > max_diff) {
      return(IRanges(-1, -1)) # no telomere was found, change to NULL
    }
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == classes$SKIP || subt$class == classes$NONE ||
        is.na(subt$class)) {
      next
    }else {
      # otherwise add one to in.a.row, update score and set start index
      start <- subt$start
      last_i <- i
      break
    }
    
  }
  
  last_i_start <- last_i
  
  for (i in last_i:nrow(subtelos)){
    subt <- subtelos[i, ]
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if(subt$class == classes$SKIP || subt$class == classes$NONE ||
       is.na(subt$class)) {
      break
    }else {
      # otherwise add one to in.a.row, update score and set start index
      end<- subt$end
    }
  }
  # check also if there is the last were it is bigger.smaller
  start_end_diff <- subtelos[last_i_start, "end_index"] - subtelos[last_i_start, "start_index"]
  if (start > end) {
    end <- start + start_end_diff
  }
  
  return(IRanges(start = start, end = end))
  
  
} # end of function 'find_left_telo'


# I need to put The classes as an arument ?
#' Have to create help-functions to this function - it is too long"
#' help_function1 : find_the_start_of_telomere
#' help_function2: find_the_end_of_telomere
#' end/start should both use the find first/last pattern at the edges
#' check_edges: see if we can go forethere with the indice: using the
#' matchpattern with mismathes = 2 + indels.
#' If density < density_thr: 
#' First check if the low density is to the right or left.
#' Then rerun with 15,10 only to one of start/end accordingly.
#' Should solve problems like 
find_telo_position <- function(seq_length, subtelos, min_in_a_row = 3,
                        min_density_score = 2) { # 15,10 for sub_length == 20
  #' @title: Find the position of the Telomere(subsequence) within the seuence.
  #' @description:  Find the start and end indices of the subsequence within the
  #'  sequence according to a data frame.
  #' @usage
  #' @param seq_length: The length of the read.
  #' @param subtelos: data frame of a subseuences indices, density and class
  #' (from the analyze_subtelos)
  #' @value An IRanges object of length 1.
  #' @return (start, end) irange  of the Telomere.
  #' @examples

  ###############3 NEED TO CHANGE FOR AN ARGUMENT OF classes ###################
  #' we have a problem in this function: Error in if (subt$class == classes$SKIP
  #'  | subt$class == classes$NONE |  : argument is of length zero
  classes <- list("CCCTAA" = -5, "NONE" = 1, "SKIP" = 0)
  # set score, start, in.a.row to 0,-1,0

  score <- 0.0
  start <- -1
  end <- -1
  in_a_row <- 0
  # check also if there is the last were it is bigger.smaller
  start_end_diff <- subtelos[1, "end_index"] - subtelos[1, "start_index"]

  # loop through subsequences
  # I need to change it to fit the case when last subseq is in range of [101:150] 
  # and  the density is 36/122 < 0.3
  end_position <- 0 # for end loop
  for (i in seq_len(nrow(subtelos))){
    subt <- subtelos[i, ]
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == classes$SKIP || subt$class == classes$NONE ||
        is.na(subt$class)) {
      score <- 0
      start <- -1
      in_a_row <- 0
    }else {
      # otherwise add one to in.a.row, update score and set start index
      in_a_row <- in_a_row + 1
      score <- score + subt$density
      if (start == -1) {
        start <- subt$start_index
      }
    }
    #' if more than MIN.IN.A.ROW subtelomeres were found and the overall score
    #' is high enough, return start index
    if (in_a_row >= min_in_a_row && score >= min_density_score) {
      end_position <- i + 1
      break
    }
  }
  if (end_position == 0) {
    return(IRanges(-1, -1)) # no telomere was found, change to NULL
  }


  #' search for end from the last subsequence (backward to finding stat incase
  #' there is island of non-telomeric subsequence)
  end <- -1
  score <- 0.0
  in_a_row <- 0
  # what if end_position >= nrow(subtelos) - min_in_a_row +1
  if(end_position >= nrow(subtelos) - min_in_a_row +1 ) {
    i <- nrow(subtelos)
    subt <- subtelos[i, ] 
    while(subt$class != classes$CCCTAA && i > end_position) {
      i <- i-1
      subt <- subtelos[i, ] 
    }
    end <- subt$end_index
  } else {
    for (i in nrow(subtelos):end_position){
      subt <- subtelos[i, ]
      # if the subsequence's class is SKIP, NONE or NA, reset values
      if (subt$class == classes$SKIP || subt$class == classes$NONE ||
          is.na(subt$class)) {
        score <- 0.0
        end <- -1
        in_a_row <- 0
      } else {
        # otherwise add one to in.a.row, update score and set start index
        in_a_row <- in_a_row + 1
        score <- score + subt$density
        if (end == -1) {
          end <- subt$end_index
        }
      }
      
      #' if more than MIN.IN.A.ROW subtelomeres were found and the overall score
      #' is high enough, return start index
      if (in_a_row >= min_in_a_row && score >= min_density_score) {
        break
      }
    }
  
  }

  if (start > end) {
    end <- start + start_end_diff
  }

  return(IRanges(start = start, end = end))
}

# if right_edge is false then we check left edge instead
find_telo_position_wraper <- function(read, patterns, with_mismatch, tvr_patterns, seq_length, subtelos, analyze_list,right_edge = FALSE) {
  
  
  # I need to change the name to subseq_irange ( 1:100, 101:200 ,....)
  telo_position <- find_telo_position(seq_length = seq_length,
                                      subtelos = analyze_list[[1]],
                                      min_in_a_row = 3, min_density_score = 2)
  # TODO: duplicated code: use a helper functio....
  
  # change the name to iranges_patterns
  irange_telo <- analyze_list[[2]][[2]]
  
  
  
  # Go over the telo_index and add gray area to the plot ....
  
  
  #' Onc ew have the Telomere indices calculate the density of the patterns
  #' within it's range.
  telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  
  #' My last change: 22/10/2023
  #' Check first the length of the telomere, before setting min_in_a_row
  num_rows <- width(telo_position) %/% global_subseq_length
  if (telo_density < 0.85 && num_rows > 5) {
    min_rows <- ifelse(num_rows <= 7, yes = num_rows - 2, no = 7)
    min_density <- 0.6*min_rows
    telo_position <- find_telo_position(seq_length = length(current_seq_unlist),
                                        subtelos = analyze_list[[1]], min_in_a_row = min_rows,
                                        min_density_score = min_density)
    
    #' Onc ew have the Telomere indices calculate the density of the patterns
    #' within it's range.
    #telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  }
  
  
  
  # I should validate start <= end
  start_acc <- get_accurate_start(telo_start = start(telo_position) , irange_telo = irange_telo)
  end_acc <- get_accurate_end(telo_end = end(telo_position), irange_telo = irange_telo)
  
  if (start_acc > end_acc) {
    end_acc <- start_acc 
  }
  
  telo_position <- IRanges(start = start_acc, end = end_acc)
  
  
  if(width(telo_position) < 100) {
    # find at the telomere at the edge 
    if(right_edge) { # right edge TTAGGG
      telo_position <- find_right_telo(seq_length = seq_length, subtelos = subtelos)
    } else { # left edge CCCTAA
      telo_position <- find_left_telo(seq_length = seq_length, subtelos = subtelos)
    }
  }
  
  # new args: read, patterns, with_mismatch, tvr_patterns
  # use match pattern steps for more accurate results 
  if(end(telo_position) < length(read)) {
    end_acc <- search_right_patterns(read = read , end_index = end(telo_position) + 1, pattern = patterns , with_mismatch = with_mismatch, tvr_patterns = tvr_patterns)
  } else {
    end_acc <- end(telo_position) 
  }
  if(start(telo_position) > 1) {
    start_acc  <- search_left_patterns(read = read, start_index = start(telo_position) -1 , pattern =  patterns, with_mismatch = with_mismatch, tvr_patterns = tvr_patterns)
  } else {
    start_acc <- start(telo_position)
  }
  
  
  telo_position <- IRanges(start = start_acc, end = end_acc)
  
  return (telo_position)
}


# check the plot
# I need to adjust the plot to my code
# dna.seq is a DNAStringSet of length 1
# This is version updated at 6.01.2022
# add the posibility for the type of picture file (JPEGS, pdf, eps ....)
plot_single_telo <- function(x_length, seq_length, subs, serial_num, seq_start,
                    seq_end, save_it = TRUE, main_title = "", w = 750, h = 300,
                    output_jpegs, eps = FALSE) { # add output_jpegs as arg
  #' @title plot the density over a sequence
  #' @param x_length: The length of the x axis.
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos
  #' function
  #' @param serial_num: The serial number of the current sequence, used as the
  #' name of the file
  #' @param seq_start: The start of the Telomere.
  #' @param seq_end: The end " ".
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param output_jpegs: the output directory for saving the file
  subs <- na.omit(subs)
  # save file if specified
    if (save_it) {
      if (isTRUE(eps)) {
        eps_path <- file.path(output_jpegs, paste0("read", serial_num, ".eps"))
        setEPS()
        #naming the eps file
        postscript(eps_path)
      } else {
        jpeg_path <- file.path(output_jpegs, paste0("read", serial_num, ".jpeg"))
        jpeg(filename = jpeg_path, width = w, height = h)
      }
  }

  ## for eps
  #if(isTRUE(save_it)){
  #  eps_path <-paste(output_jpegs, paste0('read', serial_num, '.eps'),sep='/')

    #setEPS()
    # naming the eps file
    #postscript(eps_path)
  #}


  # 26-07: my addition: save the csv file subs
  #' write_csv(x = subs, file = paste(output_jpegs, paste('read', serial_num,
  #' '.csv',sep=''), sep='/') )

  # give extra x for the legend at the topRigth
  plot(subs$density ~ subs$start_index, type = "n", yaxt = "n", xaxt = "n",
       ylab = "", xlab = "", ylim = c(0, 1), xlim = c(1, x_length +
                                                       round(x_length / 4.15)))
  # create axes
  xpos <- seq(1, x_length, by = 1000) # I have cahnged from 0 to 1
  axis(1, at = xpos, labels = sprintf("%.1fkb", xpos / 1000))
  title(xlab = "Position", adj = 0)
  axis(2, at = seq(-0.1, 1, by = 0.1), las = 2)
  
  
  #' Last change: try more accurate polygon 13:06 27/9/2023 - change
  #' c(0, subs$density, 0) to c(0, subs$density, dplyr::last(subs$density), 0)
  #' x = c(1, subs$start_index, seq_length) to x = c(1, subs$start_index, seq_length, seq_length)
  
  
  # add polygon to plot for each variant repeat.
  # mychange: only comp_ttaggg
  suppressWarnings(polygon(y = c(0, subs$density, dplyr::last(subs$density), 0),
      x = c(1, subs$start_index, seq_length, seq_length), col = rgb(1, 0, 0, 0.5),
      lwd = 0.5))
  
  # no telo indices
  if(seq_start == -1) {
    rect(xleft = 1, ybottom = -0.1, xright = seq_length, ytop = 0,
         col = "blue")
    abline(h = 1, col = "black", lty = 2)
  abline(h = 0, col = "black", lty = 2)
  legend(x = x_length, y = 1, legend = c("telomere", "sub-telomere"),
         col = c("red", "blue"), lty = 1, lwd = 2, cex = 1.2)
  sub_title <- paste("read length:", seq_length, ", No telomere length")
  title(main = main_title, sub = sub_title, ylab = "Density")
  dev.off()
  return()
  }
  
  
  rect(xleft = seq_start, ybottom = -0.1, xright = seq_end, ytop = 0,
       col = "red")
  rect(xleft = seq_end + 1, ybottom = -0.1, xright = seq_length, ytop = 0,
       col = "blue")
  if (seq_start > 1) {
    rect(xleft = 1, ybottom = -0.1, xright = seq_start, ytop = 0, col = "blue")
  }

  abline(h = 1, col = "black", lty = 2)
  abline(h = 0, col = "black", lty = 2)
  legend(x = x_length, y = 1, legend = c("telomere", "sub-telomere"),
         col = c("red", "blue"), lty = 1, lwd = 2, cex = 1.2)
  sub_title <- paste("read length:", seq_length, ", telomere length:",
                     abs(seq_start - seq_end) + 1)
  title(main = main_title, sub = sub_title, ylab = "Density")
  dev.off()

}


# check the plot
# I need to adjust the plot to my code
# dna.seq is a DNAStringSet of length 1
# This is version updated at 6.01.2022
# add the posibility for the type of picture file (JPEGS, pdf, eps ....)
# Change the function to return the plot instead
plot_single_telo_with_gray_area <- function(x_length, seq_length, subs,subs_mismatch,  serial_num, seq_start,
  seq_end,gray_start, gray_end, save_it = TRUE, main_title = "", w = 750, h = 300,
  output_jpegs, eps = FALSE) { # add output_jpegs as arg
  #' @title plot the density over a sequence
  #' @param x_length: The length of the x axis.
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos
  #' function
  #' @param subs_mismatch: the data frame with the mismatch
  #' @param serial_num: The serial number of the current sequence, used as the
  #' name of the file
  #' @param seq_start: The start of the Telomere .
  #' @param seq_end: The end " ".
  #' @param gray_start: The start of the Telomere with mismatches.
  #' @param gray_ends: The end " ".
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param output_jpegs: the output directory for saving the file
  subs <- na.omit(subs)
  # save file if specified
  if (save_it) {
    if (isTRUE(eps)) {
      eps_path <- file.path(output_jpegs, paste0("read", serial_num, ".eps"))
      setEPS()
      #naming the eps file
      postscript(eps_path)
    } else {
      jpeg_path <- file.path(output_jpegs, paste0("read", serial_num, ".jpeg"))
      jpeg(filename = jpeg_path, width = w, height = h)
    }
  }
  
  ## for eps
  #if(isTRUE(save_it)){
  #  eps_path <-paste(output_jpegs, paste0('read', serial_num, '.eps'),sep='/')
  
  #setEPS()
  # naming the eps file
  #postscript(eps_path)
  #}
  
  
  # 26-07: my addition: save the csv file subs
  #' write_csv(x = subs, file = paste(output_jpegs, paste('read', serial_num,
  #' '.csv',sep=''), sep='/') )
  
  # give extra x for the legend at the topRigth
  plot(subs$density ~ subs$start_index, type = "n", yaxt = "n", xaxt = "n",
       ylab = "", xlab = "", ylim = c(0, 1), xlim = c(1, x_length +
                                                        round(x_length / 4.15)))
  # create axes
  xpos <- seq(1, x_length, by = 1000) # I have cahnged from 0 to 1
  axis(1, at = xpos, labels = sprintf("%.1fkb", xpos / 1000))
  title(xlab = "Position", adj = 0)
  axis(2, at = seq(-0.1, 1, by = 0.1), las = 2)
  # add polygon to plot for each variant repeat.
  # mychange: only comp_ttaggg
  
  # add the plot of the density with Mismatches
  suppressWarnings(polygon(y = c(0, subs_mismatch$density, dplyr::last(subs_mismatch$density), 0),
                           x = c(1, subs_mismatch$start_index, seq_length, seq_length), col = "orange",
                           lwd = 0.5))
  
  
  suppressWarnings(polygon(y = c(0, subs$density, dplyr::last(subs$density), 0),
                           x = c(1, subs$start_index, seq_length, seq_length), col = "salmon",
                           lwd = 0.5))
  
  
  
  if(seq_start > -1) {
    rect(xleft = seq_start, ybottom = -0.1, xright = seq_end, ytop = 0,
         col = "red")
    rect(xleft = seq_end + 1, ybottom = -0.1, xright = seq_length, ytop = 0,
         col = "blue")
    if (seq_start > 1) {
      rect(xleft = 1, ybottom = -0.1, xright = seq_start, ytop = 0, col = "blue")
    }
    
    
    #  My last change: 1/10/2023
    # I need to check if gray_start !=  -1
    if(gray_start == -1) {
      abline(h = 1, col = "black", lty = 2)
      abline(h = 0, col = "black", lty = 2)
      legend(x = x_length, y = 1, 
             legend = c("telomere", "gray area" ,"sub-telomere", "Density", 
                        "Density MM" ),
             col = c("red","yellow","blue", "salmon", "orange"), lty = 1, lwd = 2, 
             cex = 1.2)
      sub_title <- paste("Read length:", seq_length, ", Telomere length:",
                         abs(seq_start - seq_end) + 1, ", Faild to calculate Telomere length with mismatches")
      title(main = main_title, sub = sub_title, ylab = "Density")
      dev.off()
      return()
    }
    
    
    
    if(gray_start < seq_start) {
      rect(xleft = gray_start, ybottom = -0.1, xright = seq_start, ytop = 0,
           col = "yellow")
    }
    
    if(gray_end > seq_end) {
      rect(xleft = seq_end, ybottom = -0.1, xright = gray_end, ytop = 0,
           col = "yellow")
    }
  } else { # no telomere without mm
  
    rect(xleft = gray_start, ybottom = -0.1, xright = gray_end, ytop = 0,
         col = "yellow")
    rect(xleft = gray_end + 1, ybottom = -0.1, xright = seq_length, ytop = 0,
         col = "blue")
    if (gray_start > 1) {
      rect(xleft = 1, ybottom = -0.1, xright = gray_start, ytop = 0, col = "blue")
    }
  }
  
  
  abline(h = 1, col = "black", lty = 2)
  abline(h = 0, col = "black", lty = 2)
  legend(x = x_length, y = 1, 
    legend = c("telomere", "gray area" ,"sub-telomere", "Density", 
    "Density MM" ),
    col = c("red","yellow","blue", "salmon", "orange"), lty = 1, lwd = 2, 
    cex = 1.2)
  
  # My last change : 2/10/2023 11:40
  # incase no telo found
  string_telo <- ifelse(seq_start == -1, yes = ", No telomere length", no = paste(", Telomere length:",
                        abs(seq_start - seq_end) + 1))
  
  sub_title <- paste("Read length:", seq_length, string_telo , ", Telomere length with mismatches:", abs(gray_start - gray_end) + 1)
  title(main = main_title, sub = sub_title, ylab = "Density")
  dev.off()
  
} # plot_single_telo_with_gray_area



plot_single_telo_with_tvr <- function(x_length, seq_length, subs,subs_mismatch, subs_tvr,  serial_num, seq_start,
         seq_end,gray_start, gray_end,tvr_start, tvr_end, save_it = TRUE, main_title = "", w = 750, h = 300,
         output_jpegs, eps = FALSE) { # add output_jpegs as arg
  #' @title plot the density over a sequence
  #' @param x_length: The length of the x axis.
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos
  #' function
  #' @param subs_mismatch: the data frame with the mismatch
  #' @param  subs_tvr: The data frame with the mismatch + tvrs
  #' @param serial_num: The serial number of the current sequence, used as the
  #' name of the file
  #' @param seq_start: The start of the Telomere .
  #' @param seq_end: The end " ".
  #' @param gray_start: The start of the Telomere with mismatches.
  #' @param gray_end: The end " ".
  #' @param tvr_start: The start of the Telomere with mismatches + tvrs.
  #' @param tvr_end: The end " "
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param output_jpegs: the output directory for saving the file
  subs <- na.omit(subs)
  # save file if specified
  if (save_it) {
    if (isTRUE(eps)) {
      eps_path <- file.path(output_jpegs, paste0("read", serial_num, ".eps"))
      setEPS()
      #naming the eps file
      postscript(eps_path)
    } else {
      jpeg_path <- file.path(output_jpegs, paste0("read", serial_num, ".jpeg"))
      jpeg(filename = jpeg_path, width = w, height = h)
    }
  }
  
  ## for eps
  #if(isTRUE(save_it)){
  #  eps_path <-paste(output_jpegs, paste0('read', serial_num, '.eps'),sep='/')
  
  #setEPS()
  # naming the eps file
  #postscript(eps_path)
  #}
  
  
  # 26-07: my addition: save the csv file subs
  #' write_csv(x = subs, file = paste(output_jpegs, paste('read', serial_num,
  #' '.csv',sep=''), sep='/') )
  
  # give extra x for the legend at the topRigth
  plot(subs$density ~ subs$start_index, type = "n", yaxt = "n", xaxt = "n",
       ylab = "", xlab = "", ylim = c(0, 1), xlim = c(1, x_length +
                                                        round(x_length / 4.15)))
  # create axes
  xpos <- seq(1, x_length, by = 1000) # I have cahnged from 0 to 1
  axis(1, at = xpos, labels = sprintf("%.1fkb", xpos / 1000))
  title(xlab = "Position", adj = 0)
  axis(2, at = seq(-0.1, 1, by = 0.1), las = 2)
  # add polygon to plot for each variant repeat.
  # mychange: only comp_ttaggg
  
  # tvr 
  suppressWarnings(polygon(y = c(0, subs_tvr$density, dplyr::last(subs_tvr$density), 0),
                           x = c(1, subs_tvr$start_index, seq_length, seq_length), col = "orange3",
                           lwd = 0.5))
  
  
  # add the plot of the density with Mismatches
  suppressWarnings(polygon(y = c(0, subs_mismatch$density, dplyr::last(subs_mismatch$density), 0),
                           x = c(1, subs_mismatch$start_index, seq_length, seq_length), col = "orange",
                           lwd = 0.5))
  
  
  suppressWarnings(polygon(y = c(0, subs$density, dplyr::last(subs$density), 0),
                           x = c(1, subs$start_index, seq_length, seq_length), col = "salmon",
                           lwd = 0.5))
  
 
  
  if(seq_start > -1) {
    rect(xleft = seq_start, ybottom = -0.1, xright = seq_end, ytop = 0,
         col = "red")
    rect(xleft = seq_end + 1, ybottom = -0.1, xright = seq_length, ytop = 0,
         col = "blue")
    if (seq_start > 1) {
      rect(xleft = 1, ybottom = -0.1, xright = seq_start, ytop = 0, col = "blue")
    }
    
    
    #  My last change: 1/10/2023
    # I need to check if gray_start !=  -1
    if(gray_start == -1) {
      abline(h = 1, col = "black", lty = 2)
      abline(h = 0, col = "black", lty = 2)
      
      
      if(tvr_start == -1) {
        legend(x = x_length, y = 1, 
               legend = c("telomere", "gray area", "tvr" ,"sub-telomere", "Density", 
                          "Density MM", "DensMM+TVRs" ),
               col = c("red","yellow","yellow3","blue", "salmon", "orange", "orange3"), lty = 1, lwd = 2, 
               cex = 1.2)
        sub_title <- paste("Read length:", seq_length, ", Telomere length:",
                           abs(seq_start - seq_end) + 1, ", Faild to calculate Telomere length with mismatches/tvr")
        title(main = main_title, sub = sub_title, ylab = "Density")
        dev.off()
        return()
      } # add tvr 
      
      if(tvr_start < seq_start) {
        rect(xleft = tvr_start, ybottom = -0.1, xright = seq_start, ytop = 0,
             col = "yellow3")
      }
      
      if(tvr_end > seq_end) {
        rect(xleft = seq_end, ybottom = -0.1, xright = tvr_end, ytop = 0,
             col = "yellow3")
      }
      
      legend(x = x_length, y = 1, 
             legend = c("telomere", "gray area", "tvr" ,"sub-telomere", "Density", 
                        "Density MM", "DensMM+TVRs" ),
             col = c("red","yellow","yellow3","blue", "salmon", "orange", "orange3"), lty = 1, lwd = 2, 
             cex = 1.2)
      sub_title <- paste("Read length:", seq_length, ", Telomere length:",
                         abs(seq_start - seq_end) + 1, ", Faild to calculate Telomere length with mismatches", "Telomere length with tvr:", abs(tvr_start - tvr_end) + 1)
      title(main = main_title, sub = sub_title, ylab = "Density")
      dev.off()
      return()
    }
    
    
    
    if(gray_start < seq_start) {
      rect(xleft = gray_start, ybottom = -0.1, xright = seq_start, ytop = 0,
           col = "yellow")
    }
    
    if(gray_end > seq_end) {
      rect(xleft = seq_end, ybottom = -0.1, xright = gray_end, ytop = 0,
           col = "yellow")
    }
    
    # check tvr
    if(tvr_start != -1) {
      if(gray_start > tvr_start) {
        rect(xleft = tvr_start, ybottom = -0.1, xright = gray_start, ytop = 0,
             col = "yellow3")
      }
      
      if(gray_end < tvr_end) {
        rect(xleft = gray_end, ybottom = -0.1, xright = tvr_end, ytop = 0,
             col = "yellow3")
      }
    }
    
    
  } else { # no telomere without mm
    
    rect(xleft = gray_start, ybottom = -0.1, xright = gray_end, ytop = 0,
         col = "yellow")
    rect(xleft = gray_end + 1, ybottom = -0.1, xright = seq_length, ytop = 0,
         col = "blue")
    if (gray_start > 1) {
      rect(xleft = 1, ybottom = -0.1, xright = gray_start, ytop = 0, col = "blue")
    }
    
    # check tvr
    if(tvr_start != -1) {
      if(gray_start > tvr_start) {
        rect(xleft = tvr_start, ybottom = -0.1, xright = gray_start, ytop = 0,
             col = "yellow3")
      }
      
      if(gray_end < tvr_end) {
        rect(xleft = gray_end, ybottom = -0.1, xright = tvr_end, ytop = 0,
             col = "yellow3")
      }
    }
    
  }
  
  
  abline(h = 1, col = "black", lty = 2)
  abline(h = 0, col = "black", lty = 2)
  legend(x = x_length, y = 1, 
         legend = c("telomere", "gray area", "tvr", "sub-telomere", "Density", 
                    "Density MM", "DensMM+TVRs"),
         col = c("red","yellow", "yellow3","blue", "salmon", "orange", "orange3"), lty = 1, lwd = 2, 
         cex = 1.2)
  
  # My last change : 2/10/2023 11:40
  # incase no telo found
  string_telo <- ifelse(seq_start == -1, yes = ", No telomere length", no = paste(", Telomere length:",
                                                                                  abs(seq_start - seq_end) + 1))
  
  sub_title <- paste("Read length:", seq_length, string_telo , ", Telomere length with mismatches:", abs(gray_start - gray_end) + 1)
  
  if(tvr_start != -1) {
    sub_title <- paste(sub_title, ", with mismatch+tvr:", abs(tvr_start - tvr_end) + 1 )
  } else {
    sub_title <- paste(sub_title, ", failed to calculate Telomere length with mismatch+tvr")
  }
  
  
  title(main = main_title, sub = sub_title, ylab = "Density")
  dev.off()
  
} # plot_single_telo_with_gray_area




# add options for the image saved:
#'ggsave currently recognises the extensions eps/ps, tex (pictex), pdf, jpeg,
#'tiff, png, bmp, svg and wmf (windows only).
plot_single_telo_ggplot2 <- function(seq_length, subs, serial_num, seq_start,
                         seq_end, save_it = TRUE,
                         main_title = "", output_jpegs) {
  #' @title plot the density over a sequence
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos
  #' function
  #' @param serial_num: The serial number of the current sequence, used as the
  #' name of the file
  #' @param seq_start: The start of the Telomere.
  #' @param seq_end: The end " ".
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param output_jpegs: the output directory for saving the file
  subs <- na.omit(subs)
  # save the densities in a csv file for later use !
  write_csv(x = subs, file = file.path(output_jpegs, paste0(serial_num, ".csv")))
  sub_title <- paste("read length:", seq_length, ", telomere length:",
                     abs(seq_start - seq_end) + 1)

  my_ggplot <- ggplot2::ggplot(subs, aes(x = start_index, y = density)) +
    geom_area(color = NA, fill = "red") +
    # the 3 last lines are for the telomere-sub-telomre region ....
    geom_rect(aes(xmin = seq_start, ymin = -0.01, xmax = seq_end, ymax = 0),
              color = "green", fill = "green")  +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(title = main_title, x = "Position", y = "Density",
         subtitle = sub_title)

  if (seq_end < seq_length) {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = seq_end, ymin = -0.01, xmax = seq_length, ymax = 0),
                color = "blue", fill = "blue")
  }


  if (seq_start > 1) {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = 1, ymin = -0.01, xmax = seq_start, ymax = 0),
                color = "blue", fill = "blue")
  }


  # save file if specified
  if (isTRUE(save_it)) {
    eps_path <- file.path(output_jpegs, paste0("gg_read", serial_num, ".eps"))
    #'ggsave currently recognises the extensions eps/ps, tex (pictex), pdf,
    #'jpeg, tiff, png, bmp, svg and wmf (windows only).
    suppressMessages(ggsave(plot = my_ggplot, device = "eps",
                            filename = eps_path))
  }
}


# get the telomere end index, with more accuracy
# find_telo_position gives only the indices according to the subseqs
# When using the package : give option for a function as argument for the end/start
get_accurate_end <- function(telo_end, irange_telo) {
  #' try more accurate: take the max of which and also check new_end >=
  #' new_start before updating the IRange object
  if(telo_end == -1) { 
    return(-1)
  }
  e_index <- telo_end
  iranges_end <- which(end(irange_telo) %in% (e_index - 99):e_index)
  if (length(iranges_end) > 0) {
    #end(telo_position) <- end(irange_telo[iranges_end[1]])}
    # take the last pattern in range ( max(end))
    e_index <- max(end(irange_telo[iranges_end]))
    # make sure end >= start
    
  }
  # take to consideration if the end is above 100length for exmple 1:216 is (1:100, 102:216)
  # check next 100 block if partial ...
    # change 20 to subseq_length/3
  iranges_end <- which(end(irange_telo) %in% (telo_end + 1 ):(telo_end+50)  )
  if (length(iranges_end) > 0) {
      #end(telo_position) <- end(irange_telo[iranges_end[1]])}
      # take the last pattern in range ( max(end))
    e_index <- max(end(irange_telo[iranges_end]))
      
  }
  
  
  return(e_index)
  
}

#' This function adjust to a pattern of length ~6.
#' @param telo_start - the start position according to subseqs
#' @param irange_telo - The paterns iranges.
get_accurate_start <- function(telo_start, irange_telo) {
  if(telo_start == -1) {
    return(telo_start)
  }
  
  s_index <- telo_start
  first_50 <- get_sub_density(IRanges(start = telo_start, width = 50), ranges = irange_telo)
  # This means probably the subseq is half telomeric
  if(first_50 < 0.3) {
    iranges_start <- which(start(irange_telo) %in% (s_index+48):(s_index + 99))
    #message("start index if the first subseq is partial: ", iranges_start2) # to remove 2
    if (length(iranges_start) > 0) {
      telo_start <- min(start(irange_telo[iranges_start]))
      #message("The start index after iranges which check for -20: ", start(telo_position2)) # to remove 5
    }
    # search for pattern before, check range of 16 bases ( for the case when 2 patterns are at the end...)
    iranges_start <- which(start(irange_telo) %in% (s_index+33):(s_index + 48))
    if (length(iranges_start) > 0) {
      telo_start <- min(start(irange_telo[iranges_start]))
      #message("The start index after iranges which check for -20: ", start(telo_position2)) # to remove 5
    }
  } else {
    # check + option to look in the prev subseq
    iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 99))
    if (length(iranges_start) > 0) {
      telo_start <- min(start(irange_telo[iranges_start]))
    }
    if(first_50 >= 0.72){
      # check before (if prev 100 have at the end TTAGGG pattern) search last 20
      iranges_start <- which(start(irange_telo) %in% (s_index-36):(s_index-1) )
      if (length(iranges_start) > 0) {
        telo_start <- min(start(irange_telo[iranges_start]))
      }
    }
    
  }
  
  return(telo_start)
}


#' get a single read find and analyze the telomeric patterns in it.
#' Plot graphs and return a summary row.
#' This function is too long - make it more readable and shorter
#' get rid of the id, instead use a shortcut od the read_id as name: once done
#' with all the reads -> reaname all the files accoding to the final df.
#' Also make dir for each different type of plot
#' right_edge - if the telomere is expected to be to the right: NNNN...NNN(TTAGGG)n
analyze_read <- function(current_seq, current_serial, pattern_list, min_density,
                         output_dir, output_jpegs, output_jpegs_1, max_length, title, tvr_patterns, right_edge = FALSE ) {
  current_fastq_name <- names(current_seq)
  current_seq_unlist <- unlist(current_seq)
  #analyze_read(current_seq, pattern_list, min_density,title, output_jpegs, output_jpegs_1, max_length )

  # returns a a list of (a data frame, list(numeric: total density,iranges))
  analyze_list <- analyze_subtelos(dna_seq = current_seq_unlist, patterns =
                                     pattern_list, min_density = min_density, sub_length = global_subseq_length, with_mismatch = FALSE)
  
  
  # new args: read, patterns, with_mismatch, tvr_patterns
    # I need to change the name to subseq_irange ( 1:100, 101:200 ,....)
  telo_position <- find_telo_position_wraper(read = current_seq_unlist, patterns = pattern_list,with_mismatch = FALSE, tvr_patterns = NULL, seq_length = length(current_seq_unlist), subtelos = analyze_list[[1]], analyze_list = analyze_list , right_edge = right_edge)  
    
  
  
  
  analyze_list2 <- analyze_subtelos(dna_seq = current_seq_unlist, patterns =
                                     pattern_list, min_density = min_density, sub_length = global_subseq_length, with_mismatch = TRUE)
  
  telo_position2 <- find_telo_position_wraper(read = current_seq_unlist, patterns = pattern_list,with_mismatch = TRUE, tvr_patterns = NULL
    , seq_length = length(current_seq_unlist), subtelos = analyze_list2[[1]], analyze_list = analyze_list2, right_edge = right_edge)
  
  
  
 
  
  telo_position3 <- NULL
  irange_telo3 <- NULL
  analyze_list3 <- NULL
  telo_density3 <- NULL
  
  # Now check with TVR's 
  if(!is.null(tvr_patterns)) {
    analyze_list3 <- analyze_subtelos(dna_seq = current_seq_unlist, patterns =
                pattern_list, min_density = min_density, sub_length = global_subseq_length, with_mismatch = TRUE, tvr_patterns = tvr_patterns)
    telo_position3 <- find_telo_position_wraper(read = current_seq_unlist, patterns = pattern_list, with_mismatch = TRUE, tvr_patterns = tvr_patterns,
                                                seq_length = length(current_seq_unlist), subtelos = analyze_list3[[1]], analyze_list = analyze_list3, right_edge = right_edge )
    
   
  }
  
  
  

  df <- data.frame(Serial = integer(), sequence_ID = character(),
                   sequence_length = integer(), telo_density = double()
                   , Telomere_start = integer(), Telomere_end = integer(),
                   Telomere_length = integer(), telo_density_mismatch = double()
                   , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                   Telomere_length_mismatch = integer() )

  if(!is.null(tvr_patterns)) {
    df <- data.frame(Serial = integer(), sequence_ID = character(),
                     sequence_length = integer(), telo_density = double()
                     , Telomere_start = integer(), Telomere_end = integer(),
                     Telomere_length = integer(), telo_density_mismatch = double()
                     , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                     Telomere_length_mismatch = integer(), 
                     telo_density_mismatch_tvr = double()
                     , Telomere_start_mismatch_tvr = integer(), Telomere_end_mismatch_tvr = integer(),
                     Telomere_length_mismatch_tvr = integer())
  }
  
  # last update
  telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  telo_density2 <- get_sub_density(telo_position2, analyze_list2[[2]][[2]])
  if(!is.null(tvr_patterns)) {
    telo_density3 <- get_sub_density(telo_position3, analyze_list3[[2]][[2]])
  }
   
  
  if (max(width(telo_position), width(telo_position2) )  < 30 && is.null(tvr_patterns)) {
    df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                  , telo_density = NA, Telomere_start = NA, Telomere_end = NA,
                  Telomere_length = NA, telo_density_mismatch = NA, 
                  Telomere_start_mismatch = NA, Telomere_end_mismatch = NA, 
                  Telomere_length_mismatch = NA)
    
    return(df)
  } # not considered a Telomere
  if(!is.null(tvr_patterns)) {
    if (max(width(telo_position), width(telo_position2), width(telo_position3) )  < 30 ) {
      df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                    , telo_density = NA, Telomere_start = NA, Telomere_end = NA,
                    Telomere_length = NA, telo_density_mismatch = NA, 
                    Telomere_start_mismatch = NA, Telomere_end_mismatch = NA, 
                    Telomere_length_mismatch = NA, telo_density_mismatch_tvr = NA, 
                    Telomere_start_mismatch_tvr = NA, Telomere_end_mismatch_tvr = NA, 
                    Telomere_length_mismatch_tvr = NA )
      
      return(df)
    } # not consid
  }
  
  # changed to compressed file: 2026-01-29
  output_telo_fasta <- paste(output_dir, paste(toString(current_serial),
                       "fasta.gz", sep = "."),  sep = "/")
  writeXStringSet(current_seq, output_telo_fasta, compress = TRUE)

  
 if(is.null(tvr_patterns)) {
   plot_single_telo_with_gray_area(x_length = max_length, seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]],
                                   serial_num = current_serial, seq_start = start(telo_position), 
                                   seq_end = end(telo_position), gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2), save_it = TRUE, main_title = title, w = 750, 
                                   h = 300, output_jpegs = output_jpegs)
   
   plot_single_telo_with_gray_area(x_length = length(current_seq_unlist), seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]],
                                   serial_num = current_serial, seq_start = start(telo_position), 
                                   seq_end = end(telo_position), gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2), save_it = TRUE, main_title = title, w = 750, 
                                   h = 300, output_jpegs = output_jpegs_1)
   
   # eps
   plot_single_telo_with_gray_area(x_length = length(current_seq_unlist), seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]],  serial_num = current_serial,
                                   seq_start = start(telo_position), seq_end = end(telo_position),gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2) ,save_it = TRUE, main_title = title,  w = 750, 
                                   h = 300,output_jpegs = output_jpegs_1, eps = TRUE)
 } else {
   plot_single_telo_with_tvr(x_length = max_length, seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]], subs_tvr = analyze_list3[[1]], 
                                   serial_num = current_serial, seq_start = start(telo_position), 
                                   seq_end = end(telo_position), gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2),tvr_start = start(telo_position3) , tvr_end = end(telo_position3), save_it = TRUE, main_title = title, w = 750, 
                                   h = 300, output_jpegs = output_jpegs)
   
   plot_single_telo_with_tvr(x_length = length(current_seq_unlist), seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]], subs_tvr = analyze_list3[[1]],
                                   serial_num = current_serial, seq_start = start(telo_position), 
                                   seq_end = end(telo_position), gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2), tvr_start = start(telo_position3) , tvr_end = end(telo_position3),save_it = TRUE, main_title = title, w = 750, 
                                   h = 300, output_jpegs = output_jpegs_1)
   
   # eps
   plot_single_telo_with_tvr(x_length = length(current_seq_unlist), seq_length =
                                     length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]], subs_tvr = analyze_list3[[1]], serial_num = current_serial,
                                   seq_start = start(telo_position), seq_end = end(telo_position),gray_start = start(telo_position2), 
                                   gray_end = end(telo_position2),tvr_start = start(telo_position3) , tvr_end = end(telo_position3) ,save_it = TRUE, main_title = title,  w = 750, 
                                   h = 300,output_jpegs = output_jpegs_1, eps = TRUE)
 }
  
  

  # check for valid telomere iranges
  telo_length <- width(telo_position)
  telo_start <- start(telo_position)
  telo_end <- end(telo_position)
  if(start(telo_position) == -1) {
    telo_length <- NA
    telo_start <- NA
    telo_end <- NA
    telo_density <- NA
  }
  telo_length2 <- width(telo_position2)
  telo_start2 <- start(telo_position2)
  telo_end2 <- end(telo_position2)
  if(start(telo_position2) == -1) {
    telo_length2 <- NA
    telo_start2 <- NA
    telo_end2 <- NA
    telo_density2 <- NA
  }
  
  
  if(is.null(tvr_patterns)) {
  df <- df %>%
    add_row(Serial = current_serial, sequence_ID = current_fastq_name,
    sequence_length = length(current_seq_unlist), telo_density = telo_density
    , Telomere_start = telo_start, Telomere_end = telo_end, Telomere_length = 
      telo_length, telo_density_mismatch = telo_density2, 
    Telomere_start_mismatch = telo_start2, Telomere_end_mismatch = telo_end2, 
    Telomere_length_mismatch = telo_length2)
  } else {
    
    telo_length3 <- width(telo_position3)
    telo_start3 <- start(telo_position3)
    telo_end3 <- end(telo_position3)
    if(start(telo_position3) == -1) {
      telo_length3 <- NA
      telo_start3 <- NA
      telo_end3 <- NA
      telo_density3 <- NA
    }
    
    df <- df %>%
      add_row(Serial = current_serial, sequence_ID = current_fastq_name,
              sequence_length = length(current_seq_unlist), telo_density = telo_density
              , Telomere_start = telo_start, Telomere_end = telo_end, Telomere_length = 
                telo_length, telo_density_mismatch = telo_density2, 
              Telomere_start_mismatch = telo_start2, Telomere_end_mismatch = telo_end2, 
              Telomere_length_mismatch = telo_length2, 
              telo_density_mismatch_tvr = telo_density3
              , Telomere_start_mismatch_tvr = telo_start3, Telomere_end_mismatch_tvr = telo_end3,
              Telomere_length_mismatch_tvr = telo_length3)
  }
  return(df)
  
} # End of analyze_read


create_dirs <- function(output_dir) {
  if (!dir.exists(output_dir)) { # update  did it
    dir.create(output_dir)
  }
  output_jpegs <- paste(output_dir, "single_read_plots", sep = "/")
  if (!dir.exists(output_jpegs)) {
    dir.create(output_jpegs)
  }
  output_jpegs_1 <- paste(output_dir, "single_read_plots_adj", sep = "/")
  if (!dir.exists(output_jpegs_1)) {
    dir.create(output_jpegs_1)
  }

  reads_dir <- paste(output_dir, "reads", sep = "/")
  if (!dir.exists(reads_dir)) {
    dir.create(reads_dir)
  }
}


############## Running functions ###############################################

search_patterns <- function(sample_telomeres, pattern_list, max_length = 1e5,
     output_dir, serial_start = 1, min_density,  title = "Telomeric repeat density", tvr_patterns, right_edge = FALSE) {
  #'@title Search given Patterns over a DNA sequences.
  #'@param sample_telomeres: the DNAString Set of the reads.
  #'@param pattern_list: a list of patterns or a string of 1 pattern.
  #'@param max_length: The x-axis length for the plot.
  #'@param output_dir:
  #'@param serial_start: The first id serial number.
  #'@param min_density: The minimal density of the patterns in a sequence to
  #'be consider relevant.
  #'@param title: The title for the density plots.
  #'@param tvr_patterns : Additional patterns for TVR's.
  #'@param right_edge : if the expected telomere is to be to the rigth : NNNN....NNN(TTAGGG)n
  #'@param return: Returns a data frame of the reads.
  

  # output_telo_csv <- paste(output_dir, paste(csv_name, "csv", sep = "."), sep = "/")
  # output_telo_fasta <- paste(output_dir, paste("reads", "fasta", sep = "."),  sep = "/")
  output_jpegs <- paste(output_dir, "single_read_plots", sep = "/")
  output_jpegs_1 <- paste(output_dir, "single_read_plots_adj", sep = "/")
  output_reads <- paste(output_dir, "reads", sep = "/")

  #max_length <- max(width(sample_telomeres))
  #  I HAVE ADDED TLOMERE LENGTH, START @ END
  
  
  
  df <- data.frame(Serial = integer(), sequence_ID = character(),
                   sequence_length = integer(), telo_density = double()
                   , Telomere_start = integer(), Telomere_end = integer(),
                   Telomere_length = integer(), telo_density_mismatch = double()
                   , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                   Telomere_length_mismatch = integer() )
  if(!is.null(tvr_patterns)) {
    df <- data.frame(Serial = integer(), sequence_ID = character(),
                     sequence_length = integer(), telo_density = double()
                     , Telomere_start = integer(), Telomere_end = integer(),
                     Telomere_length = integer(), telo_density_mismatch = double()
                     , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                     Telomere_length_mismatch = integer(), 
                     telo_density_mismatch_tvr = double()
                     , Telomere_start_mismatch_tvr = integer(), Telomere_end_mismatch_tvr = integer(),
                     Telomere_length_mismatch_tvr = integer())
  }
  
  
  # add telo density : get_sub_density <- function(sub_irange, ranges){
  # For the fasta output of the reads which pass the filter
  #telomeres_dna_string_set <- DNAStringSet()
  current_serial <- serial_start
  # try map instead
  for (i in seq_along(sample_telomeres)) {
    #' now make corrections with mismax+indels patterns both for start and end
    #' of the telomere
    #telo_position <- gray_area_adges(telo_position, current_seq)
    curr_df <- analyze_read(current_seq = sample_telomeres[i], current_serial =
               current_serial, pattern_list = pattern_list, min_density =
               min_density, output_dir = output_reads, output_jpegs = output_jpegs,
               output_jpegs_1 = output_jpegs_1, max_length = max_length, title = title, tvr_patterns = tvr_patterns, right_edge = right_edge)
    # not a telomeric read
    if (is.na(curr_df[1,1])) {
      next

    } else {
      df <- dplyr::rows_append(df, curr_df)
    }

    #telomeres_dna_string_set <- append(telomeres_dna_string_set, values = sample_telomeres[i])
    current_serial <- current_serial + 1
  } # end of for loop
  # need to save the df in a file
  #write.csv(x = df, file = output_telo_csv, row.names = FALSE)
  #writeXStringSet(telomeres_dna_string_set, output_telo_fasta)
  #message("Done!")
  #'  now what's left is to extract the sequences from fasta to fasta using
  #'  the read_names list file ( 3 files )
  return(df)
} # end of the function search_patterns



# need to create parApply for filter
filter_density <- function(sequence, patterns, min_density = 0.18) {
  #'
  total_density <- 0
  # union of all the IRanges of all the patterns in the list
  mp_all <- IRanges::IRanges()
  if (is.list(patterns)) {
    patterns <- unique(patterns)  # make sure there are no dups
    for (pat in patterns){
      mp_all <- IRanges::union(mp_all, Biostrings::matchPattern(pattern = pat,
                subject = unlist(sequence), max.mismatch = 0, fixed = FALSE))
    }
  } else {
    mp_all <- Biostrings::matchPattern(pattern = patterns,
              subject = unlist(sequence), max.mismatch = 0, fixed = FALSE)
    mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
  }

  total_density <- sum(IRanges::width(mp_all)) / nchar(sequence)
  return(total_density >= min_density)

}


create_sample <- function(input_path, format = "fastq") {
  #' @title: Extract DA sampe from fasta/q files.
  #' @description By given input files(fastq otr fasta) creates a DNAStringSet
  #' object.
  #' @param input_path: path to the file or directory containing files.
  #' @param format: The file/files format should be either fastq format or
  #'        fasta format gz extension is supported.
  if (dir.exists(input_path)) {
    sample <- Biostrings::readDNAStringSet(filepath = dir(full.names = TRUE,
                                            path = input_path, recursive = TRUE, include.dirs = FALSE), format = format)
  }else { # it is a single file path
    sample <- Biostrings::readDNAStringSet(filepath = input_path,
                                           format = format)
  }
  return(sample)
}

filter_reads <- function(samples,  patterns, do_rc = TRUE, num_of_cores = 10, subread_width = 200, right_edge = TRUE, trimm_length = 70) {
  samps_1000 <- samples[width(samples) >= 1e3]
  if (isTRUE(do_rc)) {
    samps_1000 <- Biostrings::reverseComplement(samps_1000)
  }
  # FROK don't work on windows os 
  cl <- makeCluster(num_of_cores, type = "FORK")
  # change to -(61+ just incase there are indels ( barcode_telorette == 61))
  if(right_edge) {
    trimm_length <- trimm_length +1
    trimm_length <- -1*trimm_length
    samp_100 <- parLapply(cl, samps_1000, subseq, end = trimm_length, width = subread_width)
  } else { # check the start of the read (left)
    samp_100 <- parLapply(cl, samps_1000, subseq, start = trimm_length + 1, width = subread_width)
  }

  stopCluster(cl)


  logical_100 <- mclapply(X = samp_100, FUN = filter_density,
      patterns = patterns, min_density = global_min_density*0.8, mc.cores = num_of_cores)


  #test_filter2 <- samps_1000[unlist(logical_100)]

  #logical_100 <-sapply(X = samp_100, FUN = filter_density, patterns = curr_patterns, min_density = 0.175)


  rm(samp_100)
  names(logical_100) <- NULL

  samps_1000 <- samps_1000[unlist(logical_100)]

  if (length(samps_1000) < 1) {
    message("No read have passed the filteration at run_with_rc_and_filter!")
    return(NA)
  }else {

    return(samps_1000)
  }
}




################## final function ######################################################

# run The future par search_pattern on small chunks inorder to save memory.
run_future_worker_chuncks <- function(input_path, output_path, format = c("fasta","fastq"), nrec = 10000, 
                                      patterns, do_rc, use_filter = FALSE, right_edge = TRUE, tvr_patterns) {
  
  if (dir.exists(input_path)) {
    filepath <- dir(full.names = TRUE, path = input_path, recursive = TRUE, include.dirs = FALSE)
  } else { 
    filepath <- input_path 
  }
  
  files <- open_input_files(filepath)
  
  i <- 0
  dna_length <- integer(0)
  
  if(is.null(tvr_patterns)) {
    df_summary <- data.frame(Serial = integer(), sequence_ID = character(),
                             sequence_length = integer(), telo_density = double()
                             , Telomere_start = integer(), Telomere_end = integer(),
                             Telomere_length = integer(), telo_density_mismatch = double()
                             , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                             Telomere_length_mismatch = integer() )
  } else {
    df_summary <- data.frame(Serial = integer(), sequence_ID = character(),
                             sequence_length = integer(), telo_density = double()
                             , Telomere_start = integer(), Telomere_end = integer(),
                             Telomere_length = integer(), telo_density_mismatch = double()
                             , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                             Telomere_length_mismatch = integer(), 
                             
                             telo_density_mismatch_tvr = double()
                             , Telomere_start_mismatch_tvr = integer(), Telomere_end_mismatch_tvr = integer(),
                             Telomere_length_mismatch_tvr = integer() )
  }
  
  
  options(future.globals.maxSize = 1048576000*1.5) # 1.5 gigabyte max
  plan(multicore, workers = 8)
  serial_start <- 1
  while (TRUE) {
    i <- i + 1
    ## Load 40 records at a time. Each new call to readDNAStringSet()
    ## picks up where the previous call left.
    dna_reads <- readDNAStringSet(files, nrec=nrec, format = format)
    
    
    if (length(dna_reads) == 0L)
      break
    
    if(do_rc) {
      dna_reads <- reverseComplement(dna_reads)
    }
    print(Sys.time())
    cat("processing chunk", i, "...\n")
    ## do something with 'dna' ...
    dna_length <- c(dna_length, width(dna_reads))
    
    if(use_filter) {
      plan(sequential)
      dna_reads <- filter_reads(samples = dna_reads , patterns = patterns , do_rc = FALSE, right_edge = right_edge, num_of_cores = 8)
      plan(multicore, workers = 8)
      if(is.na(dna_reads[1])) {next}
    }
    
    groups_length <- 8
    seq_over_length <- seq.int(from = 1, by = 1, length.out = length(dna_reads))
    if( length(seq_over_length) < groups_length) {
      plan(sequential)
      curr_df <- search_patterns(sample_telomeres = dna_reads, pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = serial_start , tvr_patterns = tvr_patterns, right_edge = right_edge)
      df_summary <- union_all(df_summary, curr_df)
    } else {
      
      groups <- 1:groups_length
      split_seq <- suppressWarnings(split(x = seq_over_length, f = groups))
      
      df1 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`1`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge )
      df2 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`2`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(split_seq$`1`) + serial_start , tvr_patterns = tvr_patterns, right_edge = right_edge)
      df3 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`3`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:2])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      df4 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`4`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:3])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      df5 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`5`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:4])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      df6 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`6`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:5])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      df7 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`7`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:6])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      df8 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`8`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:7])) + serial_start, tvr_patterns = tvr_patterns, right_edge = right_edge)
      
      df_summary <- Reduce(union_all , list(df_summary, df1,df2,df3, df4, df5, df6, df7, df8))
      # end of parallel try
    }
    
    serial_start <- max(df_summary$Serial) + 1
    
  }
  
  plan(sequential)
  ans_list <- list(df_summary, dna_length)  
  names(ans_list) <- c("df_summary", "all_reads_length_vec")
  
  return(ans_list)
  
}


# TODO
# todo:
#' 2. Make a package
#' 3. upload to git
#' 4. Add a subtelo- only reads (Telomers are trimmed off)
#' 5. Give trimm barcode option (length of nt to trimm..1300/30)
#' 7. Create  a log file : https://cran.r-project.org/web/packages/logr/vignettes/logr.html
#' 8. Use the suppressMessages or suppressPackageStartupMessages(library_wich give as saving 7*7....) : https://campus.datacamp.com/courses/defensive-r-programming/early-warning-systems?ex=5
#' 9. Use message() to show to process of work % and how much time it takes ....
#' 10. Check the warrnings : see https://campus.datacamp.com/courses/defensive-r-programming/early-warning-systems?ex=10
#' 11. Add error messages using stop for bad arguments... see https://campus.datacamp.com/courses/defensive-r-programming/early-warning-systems?ex=13
#' 12. Use try for fixing bad arguments such as: 'fasa -> use fasta instead ....
#' 13. Make parallel option for searh_patterns - done using future ( dn't work on windows os)
#' 14. some reads (about 1 %) do,'t pass the filteration but are found to be telomeric with search_patterns: as the telomere start after more than 160 first bases
#' 15. Add interactive plotly plot with telomeric density and sequence letters on buttom + ability to zoom in/out ....
#' 16. Add the matchPatterns with indels and max.mis == 2 at the edges -> plot at the graph as yellow for pattern with error....
#' 17. Swith instersect.Vector/union.Vector to IRange intersect  - done.
#' 18. use assertive  https://www.rdocumentation.org/packages/assertive/versions/0.3-6 for checking assert_is_numeric(..) ....
#' 19. Use cuncks inorder to save memory
#' 20. Try use parApply instead of future: see Parallel Programming in R start of chapter 2  
#' ( with FORK it should work beacause of shared memo: mclapply is faster (uses fork))
#' 21. Use the microbenchmarck and profvis toos to find were we can make it more efficient.
#' (see Parallel Programming in R ,chapter 1 )
#' 22. Use another filter: count "TTAGGG" pattern and if it is 5-10% (depends 
#' on the read length and sample) count it as a suspect for telomeric read and
#' if is is filtered out check why. (could be read with a bug on indices)
#' 23. Perhaps splitting telomere into 2 or more telomeric region will solve 
#' problems of incorrect indices or [-1,-1] indices.





if(is.null(opt$patterns)) {
  stop("Missing required parameter:  --patterns", call.=FALSE)
}

if (is.null(opt$save_path)) {
  stop("Missing required parameter:  --save_path", call.=FALSE)
}

if (is.null(opt$i)) {
  stop("Missing required parameter:  --input_path", call.=FALSE)
}

if( is.null(opt$format)) {
  stop("Format should be a string fastq or fasta")
}


# check pattern
extract_patterns <- compose(unlist, partial(str_split, pattern = "\\s+") )
cur_patterns <- extract_patterns(opt$patterns)
if(length(cur_patterns) > 1) {
  cur_patterns <- as.list(cur_patterns)
}

cur_tvr_patterns <- NULL
if(!is.null(opt$tvr_patterns)) {
  cur_tvr_patterns <- extract_patterns(opt$tvr_patterns)
  if(length(cur_tvr_patterns) > 1) {
    cur_tvr_patterns <- as.list(cur_tvr_patterns)
  }
}


global_min_density <- opt$min_density
global_subseq_length <- opt$subseq_length

lockBinding("global_subseq_length", globalenv())
lockBinding("global_min_density", globalenv())


# log function blueprint : summary(sample), %telomeric_reads , summary(telo_read) ...
# I need to create a log funcion : with ifelse ( if telomeric patterns were found or not -> no one passed the filteration or df isempty ....)
# test log file
tmp <- file.path(opt$save_path, "run.log")

# Open log
lf <- log_open(tmp) 'Telomere Analyzer  version v1.1.7-beta 2026-02-18'
log_print('Telomere Analyzer  version v1.1.7-beta 2026-02-18', hide_notes = TRUE, console = FALSE) # Send message to log
t1 <- Sys.time()
log_print(base::paste("Work started at:", toString(t1)), hide_notes = TRUE, console = FALSE) # Send message to log

log_print("############### The input argumetns for this run: ################", hide_notes = TRUE, console = FALSE)
if(opt$r) {
  log_print("Reverse complement was applied on the input reads.", hide_notes = TRUE, console = FALSE)
}
log_print(base::paste("The patterns to search:", opt$patterns), hide_notes = TRUE, console = FALSE)
log_print(base::paste("The sub-sequence length  is:", 
                      toString(opt$subseq_length)), hide_notes = TRUE, console = FALSE)
log_print(base::paste("The minimal density for a telomeric subseq:", 
                      toString(opt$min_density)), hide_notes = TRUE, console = FALSE)

if(!is.null(cur_tvr_patterns) ) {
  log_print(base::paste("Additional Telomere variant repeats patterns were added:", 
                        opt$tvr_patterns), hide_notes = TRUE, console = FALSE)
}


log_print("##################################################################",hide_notes = TRUE, , console = FALSE)

log_print("The input files:", hide_notes = TRUE, console = FALSE)
# add the names of the files which we analyze.


if(dir.exists(opt$i) ){
  filepath = dir(full.names = TRUE,
                 path = opt$i, recursive = TRUE, include.dirs = FALSE)
  
  for(i in seq_along(filepath)) {
    log_print(filepath[i], hide_notes = TRUE, console = FALSE)
  }
} else {
  log_print(opt$i, hide_notes = TRUE, console = FALSE)
}

create_dirs(output_dir = opt$save_path)


  
ans_list <- run_future_worker_chuncks(input_path = opt$i, output_path = opt$save_path, 
    format = opt$format, nrec = opt$nrec, patterns = cur_patterns, do_rc = opt$r, 
    use_filter = opt$use_filter, right_edge = opt$check_right_edge, tvr_patterns = cur_tvr_patterns)

global_total_length <- length(ans_list$all_reads_length_vec)

log_print(base::paste("Total reads in sample:", global_total_length), hide_notes = TRUE, console = FALSE)
log_print("Summary statistics of the sample reads length:", hide_notes = TRUE, console = FALSE)
log_print(summary(ans_list$all_reads_length_vec), hide_notes = TRUE, console = FALSE)

# try FastqStreamerList for streaming, n records each yield ....

#now try parallel using future package
# prev cod
# df_summary <- search_patterns(sample_telomeres = dna_reads, pattern_list = curr_patterns, output_dir = args[2], min_density = 0.3, serial_start = )



log_print(base::paste0("Number of reads which identified as Telomeric: ", nrow(ans_list$df_summary)), hide_notes = TRUE, console = FALSE)
log_print(base::paste0("% of total reads: ", toString(round( (100*nrow(ans_list$df_summary)) / global_total_length, 2 ) ), "%" ), hide_notes = TRUE, console = FALSE)

# summary statistics of the Telomeric reads
log_print("Summary statistics for the Telomeric reads:", hide_notes = TRUE, console = FALSE)
log_print("reads length:", hide_notes = TRUE, console = FALSE)
log_print(summary(ans_list$df_summary$sequence_length), hide_notes = TRUE, console = FALSE)
log_print("Telomere length:", hide_notes = TRUE, console = FALSE)
log_print(summary(ans_list$df_summary$Telomere_length), hide_notes = TRUE, console = FALSE)

log_print("Telomere length with 1 mismatch allowed:", hide_notes = TRUE, console = FALSE)
log_print(summary(ans_list$df_summary$Telomere_length_mismatch), hide_notes = TRUE, console = FALSE)


if(!is.null(cur_tvr_patterns)) {
  log_print('Telomere length with 1 mismatch allowed + tvr patterns.:', hide_notes = TRUE, console = FALSE)
  log_print(summary(ans_list$df_summary$Telomere_length_mismatch_tvr), hide_notes = TRUE, console = FALSE)
}


write_csv(x = ans_list$df_summary, file = file.path(opt$save_path,"summary.csv"))
write_lines(x = ans_list$df_summary$sequence_ID  , file = file.path(opt$save_path, "reads_ids.txt"))
t2 <- Sys.time()


log_print(base::paste("Work ended at:", toString(t2)), hide_notes = TRUE, console = FALSE)
# Close log
log_close(footer = FALSE)
writeLines(readLines(lf))







