#!/usr/bin/env Rscript

################################################################
# Author: Dan Lichtental
# Copyright (c) Dan Lichtental, 2022
# Email:  dan.lichtental@mail.huji.ac.il
# Last updated:  2023-Oct-26
# Script Name: Telomere pattern finder
# Script Description: In this script we search over fastq file sequences for
# telomeric patterns.
#
# Notes: Running the script on linux:
# "Rscript --vanilla NanoTel -i input_path --save_path output_path --patterns "some patterns"
# See Rscript --vanilla nanotel-multicore-10workers.R --help for more additional flags
# format can be either fasta or fastq, input_path can be either a full path for a
# single fastq/a file or a directory containing fastq/a files.
# Please try to give the output_path outside the input dir path and vice versa

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
suppressPackageStartupMessages(require(optparse))
#utils::globalVariables(c("start_index"))

#' TODO: replace all for loops to map with a private functio
#' instead of for(...{if...else...} : map(x, p_function)
#' map(numlist, ~.x %>% sqrt() %>% sin())

# Last change: 16:09 4/10/2023 max-matches <- 2 (insteadof 1)


global_min_density <- 0.6
global_subseq_length <- 100


#' bug fix of faild calculating the telomere indcies after re-telo_position:
#' Check: My last change : 22/10/2023


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
get_density_iranges <- function(sequence, patterns, fixed = TRUE, with_mismatch = FALSE) {
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a
  #' list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
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
      curr_mp <- matchPattern(pattern = pat, subject =
              unlist(sequence), max.mismatch = max_mismatch, fixed = fixed )
      if(max_mismatch) { curr_mp <- trim(curr_mp)}
      mp_all <- IRanges::union(mp_all, curr_mp)
    }
  }else {
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence),
                           max.mismatch = max_mismatch, fixed = fixed)
    if( (fixed == FALSE) || (max_mismatch > 0) ) {
      mp_all <- trim(mp_all)
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
                       min_density, fixed = TRUE, with_mismatch = FALSE) {
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
  list_density_mp <- get_density_iranges(dna_seq, patterns = patterns, fixed = fixed , with_mismatch = with_mismatch)
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
  #' @param seq_start: The start of the Telomere with mismatches.
  #' @param gray_end: The end " ".
  #' @param gray_start: The start of the Telomere with mismatches.
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
  
}





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
analyze_read <- function(current_seq, current_serial, pattern_list, min_density,
                         output_dir, output_jpegs, output_jpegs_1, max_length, title) {
  current_fastq_name <- names(current_seq)
  current_seq_unlist <- unlist(current_seq)
  #analyze_read(current_seq, pattern_list, min_density,title, output_jpegs, output_jpegs_1, max_length )

  # returns a a list of (a data frame, list(numeric: total density,iranges))
  analyze_list <- analyze_subtelos(dna_seq = current_seq_unlist, patterns =
                                     pattern_list, min_density = min_density, sub_length = global_subseq_length, fixed = TRUE, with_mismatch = FALSE)
  
  # I need to change the name to subseq_irange ( 1:100, 101:200 ,....)
  telo_position <- find_telo_position(seq_length = length(current_seq_unlist),
                                    subtelos = analyze_list[[1]],
                                    min_in_a_row = 3, min_density_score = 2)
  # TODO: duplicated code: use a helper functio....
  
  # change the name to iranges_patterns
  irange_telo <- analyze_list[[2]][[2]]

  
  analyze_list2 <- analyze_subtelos(dna_seq = current_seq_unlist, patterns =
                                     pattern_list, min_density = min_density, sub_length = global_subseq_length, fixed = TRUE, with_mismatch = TRUE)
  telo_position2 <- find_telo_position(seq_length = length(current_seq_unlist),
                                    subtelos = analyze_list2[[1]],
                                    min_in_a_row = 3, min_density_score = 2)
  
  irange_telo2 <- analyze_list2[[2]][[2]]
  
  

  df <- data.frame(Serial = integer(), sequence_ID = character(),
                   sequence_length = integer(), telo_density = double()
                   , Telomere_start = integer(), Telomere_end = integer(),
                   Telomere_length = integer(), telo_density_mismatch = double()
                   , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                   Telomere_length_mismatch = integer() )

  
  
  
  if (max(width(telo_position), width(telo_position2) )  < 100) {
    df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
          , telo_density = NA, Telomere_start = NA, Telomere_end = NA,
          Telomere_length = NA, telo_density_mismatch = NA, 
          Telomere_start_mismatch = NA, Telomere_end_mismatch = NA, 
          Telomere_length_mismatch = NA)
    
    return(df)
  } # not considered a Telomere
  
  
  
  # Go over the telo_index and add gray area to the plot ....
  
  
  #' Onc ew have the Telomere indices calculate the density of the patterns
  #' within it's range.
  telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  
  #' My last change: 22/10/2023
  #' Check first the length of the telomere, before setting min_in_a_row
  num_rows <- width(telo_position) %/% global_subseq_length
  if (telo_density < 0.7 && num_rows > 5) {
    min_rows <- ifelse(num_rows <= 10, yes = num_rows - 2, no = 10)
    min_density <- 0.6*min_rows
    telo_position <- find_telo_position(seq_length = length(current_seq_unlist),
                                      subtelos = analyze_list[[1]], min_in_a_row = min_rows,
                                      min_density_score = min_density)
    
    #' Onc ew have the Telomere indices calculate the density of the patterns
    #' within it's range.
    #telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  }
  
  
  telo_density2 <- get_sub_density(telo_position2, analyze_list2[[2]][[2]])
  
  #' My last change: 22/10/2023
  #' Check first the length of the telomere, before setting min_in_a_row
  num_rows <- width(telo_position2) %/% global_subseq_length
  if (telo_density2 < 0.7 && num_rows > 5) {
    min_rows <- ifelse(num_rows <= 10, yes = num_rows - 2, no = 10)
    min_density <- 0.6*min_rows
    telo_position2 <- find_telo_position(seq_length = length(current_seq_unlist),
                                       subtelos = analyze_list2[[1]], min_in_a_row = min_rows,
                                       min_density_score = min_density)
    
    #' Onc ew have the Telomere indices calculate the density of the patterns
    #' within it's range.
    #telo_density2 <- get_sub_density(telo_position2, analyze_list2[[2]][[2]])
  }
  

  
  # I should validate start <= end
  start_acc <- get_accurate_start(telo_start = start(telo_position) , irange_telo = irange_telo)
  start2_acc <- get_accurate_start(telo_start = start(telo_position2) , irange_telo = irange_telo2 )
  end_acc <- get_accurate_end(telo_end = end(telo_position), irange_telo = irange_telo)
  end2_acc <- get_accurate_end(telo_end = end(telo_position2), irange_telo = irange_telo2 )
  
  
  if (start_acc > end_acc) {
    end_acc <- start_acc 
  }
  if (start2_acc > end2_acc) {
    end2_acc <- start2_acc 
  }
  telo_position <- IRanges(start = start_acc, end = end_acc)
  telo_position2 <- IRanges(start = start2_acc, end = end2_acc)
  
  
  # last update
  telo_density <- get_sub_density(telo_position, analyze_list[[2]][[2]])
  telo_density2 <- get_sub_density(telo_position2, analyze_list2[[2]][[2]])
  
  if (max(width(telo_position),width(telo_position2)) < 100) {
    df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                  , telo_density = NA, Telomere_start = NA, Telomere_end = NA,
                  Telomere_length = NA, telo_density_mismatch = NA, 
                  Telomere_start_mismatch = NA, Telomere_end_mismatch = NA, 
                  Telomere_length_mismatch = NA)
      return(df)
  }  # not considered a Telomere

  
  
  output_telo_fasta <- paste(output_dir, paste(toString(current_serial),
                       "fasta", sep = "."),  sep = "/")
  writeXStringSet(current_seq, output_telo_fasta)

  plot_single_telo(x_length = max(max_length, length(current_seq_unlist)),
        seq_length = length(current_seq_unlist), subs =  analyze_list[[1]],
        serial_num = current_serial, seq_start = start(telo_position),
        seq_end = end(telo_position), save_it = TRUE, main_title = title,
        w = 750, h = 300, output_jpegs = output_jpegs)

  plot_single_telo_with_gray_area(x_length = length(current_seq_unlist), seq_length =
    length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]],
    serial_num = current_serial, seq_start = start(telo_position), 
    seq_end = end(telo_position), gray_start = start(telo_position2), 
    gray_end = end(telo_position2), save_it = TRUE, main_title = title, w = 750, 
    h = 300, output_jpegs = output_jpegs_1)

  # eps
  plot_single_telo_with_gray_area(x_length = length(current_seq_unlist), seq_length =
    length(current_seq_unlist), subs =  analyze_list[[1]],subs_mismatch = analyze_list2[[1]],serial_num = current_serial,
    seq_start = start(telo_position), seq_end = end(telo_position),gray_start = start(telo_position2), 
    gray_end = end(telo_position2) ,save_it = TRUE, main_title = title,  w = 750, 
    h = 300,output_jpegs = output_jpegs_1, eps = TRUE)

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
  
  
  
  df <- df %>%
    add_row(Serial = current_serial, sequence_ID = current_fastq_name,
    sequence_length = length(current_seq_unlist), telo_density = telo_density
    , Telomere_start = telo_start, Telomere_end = telo_end, Telomere_length = 
      telo_length, telo_density_mismatch = telo_density2, 
    Telomere_start_mismatch = telo_start2, Telomere_end_mismatch = telo_end2, 
    Telomere_length_mismatch = telo_length2)

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
                            output_dir, serial_start = 1,
                           min_density,  title = "Telomeric repeat density") {
  #'@title Search given Patterns over a DNA sequences.
  #'@param sample_telomeres: the DNAString Set of the reads.
  #'@param pattern_list: a list of patterns or a string of 1 pattern.
  #'@param max_length: The x-axis length for the plot.
  #'@param output_dir:
  #'@param serial_start: The first id serial number.
  #'@param min_density: The minimal density of the patterns in a sequence to
  #'be consider relevant.
  #'@param title: The title for the density plots.
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
               output_jpegs_1 = output_jpegs_1, max_length = max_length, title = title)
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

run_with_filter <- function(samples,  patterns, output_dir, do_rc = TRUE,
                                   serial_start = 1, num_of_cores = 10) {
  #' @title: Run a search for Telomeric sequences on the reads
  #' @description use rc to adjust for the patterns and barcode/telorette l
  #' ocatio ( should be at the last ~ 60-70 bases), filter out reads with no
  #'              no telomeric pattern at the edge and run search_patterns.
  #'              We assume that the Telomeres are in the right edge, or at the
  #'              left edge and than we make the reverse comlement.
  #' @usage

  if (!dir.exists(output_dir)) { # update  did it
    dir.create(output_dir)
  }

  samples_filtered <- filter_reads(samples = samples, patterns = patterns,
                      do_rc = do_rc, num_of_cores = num_of_cores)

  if (is.na(samples_filtered[1])) {
    return(NA)
  }else {
    search_patterns(samples_filtered, pattern_list = patterns, max_length =
                  max(max(width(samples_filtered)), 150000), output_dir =
                  output_dir, min_density = global_min_density, serial_start = serial_start)
  }
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
                                      patterns, do_rc, use_filter = FALSE) {
  
  if (dir.exists(input_path)) {
    filepath <- dir(full.names = TRUE, path = input_path, recursive = TRUE, include.dirs = FALSE)
  } else { 
    filepath <- input_path 
  }
  
  files <- open_input_files(filepath)
  
  i <- 0
  dna_length <- integer(0)
  df_summary <- data.frame(Serial = integer(), sequence_ID = character(),
                           sequence_length = integer(), telo_density = double()
                           , Telomere_start = integer(), Telomere_end = integer(),
                           Telomere_length = integer(), telo_density_mismatch = double()
                           , Telomere_start_mismatch = integer(), Telomere_end_mismatch = integer(),
                           Telomere_length_mismatch = integer() )
  
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
      dna_reads <- filter_reads(samples = dna_reads , patterns = patterns , do_rc = FALSE, right_edge = TRUE, num_of_cores = 8)
      plan(multicore, workers = 8)
      if(is.na(dna_reads[1])) {next}
    }
    
    groups_length <- 8
    seq_over_length <- seq.int(from = 1, by = 1, length.out = length(dna_reads))
    if( length(seq_over_length) < groups_length) {
      plan(sequential)
      curr_df <- search_patterns(sample_telomeres = dna_reads, pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = 1 )
      df_summary <- union_all(df_summary, curr_df)
    } else {
      
      groups <- 1:groups_length
      split_seq <- suppressWarnings(split(x = seq_over_length, f = groups))
      
      df1 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`1`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = serial_start )
      df2 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`2`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(split_seq$`1`) + serial_start )
      df3 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`3`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:2])) + serial_start)
      df4 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`4`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:3])) + serial_start)
      df5 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`5`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:4])) + serial_start)
      df6 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`6`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:5])) + serial_start)
      df7 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`7`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:6])) + serial_start)
      df8 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`8`], pattern_list = patterns, output_dir = output_path, min_density = global_min_density, serial_start = length(unlist(split_seq[1:7])) + serial_start)
      
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
  make_option(c("-f","--format"), action = "store", 
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
              metavar = "Sub--sequence length." ), 
  
  make_option("--use_filter", action = "store_true", default = FALSE,
              help = "Filter reads accoding to the edge.", 
              metavar = "USE FILTER" )
  
)   
opt = parse_args(OptionParser(option_list=option_list))

if(is.null(opt$patterns)) {
  stop("Missing required parameter:  --patterns", call.=FALSE)
}

if (is.null(opt$save_path)) {
  stop("Missing required parameter:  --save_path", call.=FALSE)
}

if (is.null(opt$i)) {
  stop("Missing required parameter:  --input_path", call.=FALSE)
}


# check pattern
extract_patterns <- compose(unlist, partial(str_split, pattern = "\\s+") )
cur_patterns <- extract_patterns(opt$patterns)
if(length(cur_patterns) > 1) {
  cur_patterns <- as.list(cur_patterns)
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
lf <- log_open(tmp)

t1 <- Sys.time()
log_print(base::paste("Work started at:", toString(t1)), hide_notes = TRUE) # Send message to log

log_print(base::paste("The patterns to search:", opt$patterns), hide_notes = TRUE)

log_print("The input files:", hide_notes = TRUE)
# add the names of the files which we analyze.


if(dir.exists(opt$i) ){
  filepath = dir(full.names = TRUE,
                 path = opt$i, recursive = TRUE, include.dirs = FALSE)
  
  for(i in seq_along(filepath)) {
    log_print(filepath[i], hide_notes = TRUE)
  }
} else {
  log_print(opt$i, hide_notes = TRUE)
}

create_dirs(output_dir = opt$save_path)

  
ans_list <- run_future_worker_chuncks(input_path = opt$i, output_path = opt$save_path, 
    format = opt$f, nrec = opt$nrec, patterns = cur_patterns, do_rc = opt$r, 
    use_filter = opt$use_filter)

global_total_length <- length(ans_list$all_reads_length_vec)

log_print(base::paste("Total reads in sample:", global_total_length), hide_notes = TRUE)
log_print("Summary statistics of the sample reads length:", hide_notes = TRUE)
log_print(summary(ans_list$all_reads_length_vec), hide_notes = TRUE)

# try FastqStreamerList for streaming, n records each yield ....

#now try parallel using future package
# prev cod
#df_summary <- search_patterns(sample_telomeres = dna_reads, pattern_list = curr_patterns, output_dir = args[2], min_density = 0.3, serial_start = )



log_print(base::paste0("Number of reads which identified as Telomeric: ", nrow(ans_list$df_summary)), hide_notes = TRUE)
log_print(base::paste0("% of total reads: ", toString(round( (100*nrow(ans_list$df_summary)) / global_total_length, 2 ) ), "%" ), hide_notes = TRUE)

# summary statistics of the Telomeric reads
log_print("Summary statistics for the Telomeric reads:", hide_notes = TRUE)
log_print("reads length:", hide_notes = TRUE)
log_print(summary(ans_list$df_summary$sequence_length), hide_notes = TRUE)
log_print("Telomere length:", hide_notes = TRUE)
log_print(summary(ans_list$df_summary$Telomere_length), hide_notes = TRUE)

log_print("Telomere length with 1 mismatch allowed:", hide_notes = TRUE)
log_print(summary(ans_list$df_summary$Telomere_length_mismatch), hide_notes = TRUE)

write_csv(x = ans_list$df_summary, file = file.path(opt$save_path,"summary.csv"))
t2 <- Sys.time()


log_print(base::paste("Work ended at:", toString(t2)), hide_notes = TRUE)
# Close log
log_close()
writeLines(readLines(lf))







