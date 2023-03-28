#!/usr/bin/env Rscript
# args for input_path, output_path  and format(fastq/a)
args <- commandArgs(trailingOnly = TRUE) 
################################################################
# Author: Dan Lichtental
# Copyright (c) Dan Lichtental, 2022
# Email:  dan.lichtental@mail.huji.ac.il
# Last updated:  2023-03-09
# Script Name: Telomere pattern finder
# Script Description: In this script we search for 
# telomeric patterns in fastq/a file sequences.
#
# Notes: Running the script on linux: 
# "Rscript --vanilla NanoTel.R input_dir output_dir format"
# Format can be either fasta or fastq, input_dir can be either a full path for a
# single fastq/a file or a directory containing fastq/a files.
# Please try to give the output dir outside the input dir path and vice versa

# The 'confilcted' package tries to make your function choice explicit.
#  it produces an error if a function name is found in 2 or more packages.
library(conflicted)
library(Biostrings)
library(tidyverse)  
library(IRanges)
library(parallel)
library(logr)
library(future)
library(ggprism)
utils::globalVariables(c("start_index"))


# global vars
global_min_density <- 0.3
lockBinding("global_min_density", globalenv())
global_subseq_length <- 100
lockBinding("global_subseq_length", globalenv())



split_telo <- function(dna_seq, sub_length) {
  #' @title Splits a DNA sequence into subsequences. 
  #' @description This function calculates the sequence length and creates IRanges 
  #' objects of subseuences of a given length.
  #' if The dna_seq%%subseq != 0   and last' width < sub_length/2 Then we will 
  #' remove this last index making the last subtelomere longer then by. 
  #' (because we want to prevent a case where we have for 
  #' example subtelomere of length < sub_length/2 which is too short to consider
  #'  -> then the last subtelomere is a bit longer)
  #' If the length of dna_seq is less then the sub_length it will return an 
  #' empty Iranges object.            
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



get_density_iranges <- function(sequence, patterns) {
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a 
  #' list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
  #' @value A numeric for the total density of the pattern(patterns) in the 
  #' sequences, an IRanges object of the indices of the patterns found.
  #' @return a tupple of (density, IRanges) Total density of all the patterns in 
  #' the list( % of the patterns in this sequence) and the IRanges of them
  #' @examples 
  total_density <- 0 
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list 
  if (is.list(patterns)) {
    patterns <- unique(patterns)  # make sure there are no dups
    for (pat in patterns){
      mp_all <- IRanges::union(mp_all, matchPattern(pattern = pat, subject = 
                       unlist(sequence), max.mismatch = 0, fixed = FALSE))
    }
  }else {
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence), 
                           max.mismatch = 0, fixed = FALSE)
    mp_all <- IRanges::union(mp_all, mp_all) # incase there are overlaps
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
#'   4. use purrr::map to replace the for loops code
#'   5. Give names to the elements in the list : for example:
#'      m_list[["subtelos]] <- subtelos  
#'      ( see the course Foundations of Functional Programming with purrr)
###########3 My cahnge from prev - return a list 0f (df, total_density) ######
analyze_subtelos <- function(dna_seq, patterns, sub_length, 
                       min_density) { 
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
  list_density_mp <- get_density_iranges(dna_seq, patterns = patterns)
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


#' Have to create help-functions to this function - it is too long"
#' help_function1 : find_the_start_of_telomere
#' help_function2: find_the_end_of_telomere
#' end/start should both use the find first/last pattern at the edges
#' check_edges: see if we can go forethere with the indice: using the 
#' matchpattern with mismathes = 2 + indels.
#' 
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
  start_end_diff <- subtelos[1, "end_index"] - subtelos[1, "start_index"] 
  
  # loop through subsequences
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
    return(IRanges(1, 1)) # no telomere was found 
  }
  
  
  #' search for end from the last subsequence (backward to finding stat incase 
  #' there is island of non-telomeric subsequence)
  end <- -1 
  score <- 0.0
  in_a_row <- 0
  for (i in nrow(subtelos):end_position){ 
    subt <- subtelos[i, ]
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == classes$SKIP || subt$class == classes$NONE || 
        is.na(subt$class)) {
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
    if (in_a_row >= min_in_a_row && score >= min_density_score) {
      break
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
  # add polygon to plot for each variant repeat. 
  # mychange: only comp_ttaggg 
  suppressWarnings(polygon(y = c(0, subs$density, 0), 
      x = c(1, subs$start_index, seq_length), col = rgb(1, 0, 0, 0.5), 
      lwd = 0.5)) 
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

# add options for the image saved:
#'ggsave currently recognizes the extensions eps/ps, tex (pictex), pdf, jpeg, 
#'tiff, png, bmp, svg and wmf (windows only).
plot_single_telo_ggplot2 <- function(seq_length, subs, serial_num, seq_start, 
                                     seq_end, save_it = TRUE, 
                                     main_title = "", output_jpegs) { 
  #' @title plot the density over a sequence
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subsequences from the analyze_subtelos 
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
  #write_csv(x = subs, file = paste(output_jpegs, paste0(serial_num, ".csv"), sep = "/"))
  
  my_ggplot <- ggplot2::ggplot(subs, aes(x = start_index/1000, y = density)) +
    geom_area(fill = "lightcoral") +
    geom_line(color = "black", linewidth = 0.01) +
    scale_size_manual(values = 0.1) + 
    geom_rect(aes(xmin = seq_start/1000, ymin = -0.02, xmax = seq_end/1000, ymax = 0, fill = "red"), 
              color = "black", linewidth = 0.3, key_glyph = draw_key_dotplot)  +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
    labs(title = "Telomeric repeat density", x = "Position (kb)", y = "Density", 
         caption = paste("read length:", seq_length, ", telomere length:", abs(seq_start - seq_end) + 1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0.5, vjust = 0, size = 10),
          axis.text=element_text(size=11),
          axis.ticks.length.x = unit(5, "pt"),
          axis.ticks.length.y = unit(4.5, "pt"),
          legend.title = element_blank(),
          legend.background = element_rect(color = "black"),
          legend.position = c(1.1, 0.8),
          legend.spacing = unit(20, "mm"),
          panel.background = element_rect(color = "black", fill = NA),
          plot.margin = margin(r = 4, t = 0.5, b = 0.5, l = 0.5, unit = "cm")) +
    coord_cartesian(clip = "off") +
    annotation_ticks(sides = "bl", minor.length = unit(2.8, "pt"), outside = TRUE) +
    scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), guide = "prism_minor", minor_breaks = seq(0.0, 1.0, 0.1)) +
    scale_x_continuous(breaks = seq(0, 100, 5), limits = c(0,100), guide = "prism_minor", minor_breaks = 1:101) +
    scale_fill_identity(guide = "legend", labels = c("Sub-telomere", "Telomere")) 
  
  
  if (1821 < 1821) {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = seq_end/1000, ymin = -0.02, xmax = seq_length/1000, ymax = 0, fill = "blue"), 
                color = "black", linewidth = 0.3, key_glyph = draw_key_dotplot) }
  
  
  if (101 > 1) {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = 0, ymin = -0.02, xmax = seq_start/1000, ymax = 0, fill = "blue"), 
                color = "black", linewidth = 0.3, key_glyph = draw_key_dotplot) }
  
  
  # save file if specified
  if (isTRUE(save_it)) {
    eps_path <- paste(output_jpegs, paste0("gg_read", serial_num, ".eps"),
                      sep = "/")  
    #'ggsave currently recognizes the extensions eps/ps, tex (pictex), pdf, 
    #'jpeg, tiff, png, bmp, svg and wmf (windows only).
    ggsave(filename = eps_path, plot = my_ggplot, device = "eps", width = 2500, height = 900, units = "px", dpi = 250)
  }
  
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
                                     pattern_list, min_density = min_density, sub_length = global_subseq_length)
  telo_irange <- find_telo_position(seq_length = length(current_seq_unlist), 
                                    subtelos = analyze_list[[1]], 
                                    min_in_a_row = 3, min_density_score = 2)
  # TODO: duplicated code: use a helper functio....
  irange_telo <- analyze_list[[2]][[2]]
  
  
  df <- data.frame(Serial = integer(), sequence_ID = character(), 
                   sequence_length = integer(), telo_density = double()
                   , Telomere_start = integer(), Telomere_end = integer(), 
                   Telomere_length = integer())
  
  if (width(telo_irange) < 100) {
    df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                  , telo_density = NA, Telomere_start = NA, Telomere_end = NA, 
                  Telomere_length = NA)
    return(df)
  } # not considered a Telomere
  s_index <- start(telo_irange) 
  # make the strat/end more accurate (usethe IRanges for the patterns)
  iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) 
  if (length(iranges_start) > 0) { 
    start(telo_irange) <- start(irange_telo[iranges_start[1]])
  } 
  #' try more accurate: take the max of which and also check new_end >= 
  #' new_start before updating the IRange object
  e_index <- end(telo_irange) 
  iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
  if (length(iranges_end) > 0) {
    #end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
    # take the last pattern in range 
    new_end <- end(irange_telo[iranges_end[length(iranges_end)]])
    # make sure end >= start 
    if (new_end >= start(telo_irange)) {
      end(telo_irange) <- new_end
    }     
  }  
  #' Onc ew have the Telomere indices calculate the density of the patterns 
  #' within it's range.
  telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
  
  if (telo_density < 0.75) {
    telo_irange <- find_telo_position(seq_length = length(current_seq_unlist), 
                   subtelos = analyze_list[[1]], min_in_a_row = 10, 
                   min_density_score = 6)
    irange_telo <- analyze_list[[2]][[2]]
    if (width(telo_irange) < 100) {
      df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                  , telo_density = NA, Telomere_start = NA, Telomere_end = NA, 
                  Telomere_length = NA)
      return(df)
    } # not considered a Telomere
    s_index <- start(telo_irange) 
    # make the strat/end more accurate (usethe IRanges for the patterns)
    iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) 
    if (length(iranges_start) > 0) {
      start(telo_irange) <- start(irange_telo[iranges_start[1]])
    } 
    #' try more accurate: take the max of which and also check new_end >= 
    #' new_start before updating the IRange object
    e_index <- end(telo_irange) 
    iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
    if (length(iranges_end) > 0) {
      #end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
      # take the last pattern in range
      new_end <- end(irange_telo[iranges_end[length(iranges_end)]]) 
      # make sure end >= start 
      if (new_end >= start(telo_irange)) {
        end(telo_irange) <- new_end
      }     
    }  
    #' Onc ew have the Telomere indices calculate the density of the patterns 
    #' within it's range.
    telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
  } 
  
  if (width(telo_irange) < 100) {
      df <- add_row(df, Serial = NA, sequence_ID = NA, sequence_length = NA
                  , telo_density = NA, Telomere_start = NA, Telomere_end = NA, 
                  Telomere_length = NA)
      return(df)
  }  # not considered a Telomere
  
  output_telo_fasta <- paste(output_dir, paste(toString(current_serial), 
                       "fasta", sep = "."),  sep = "/")
  writeXStringSet(current_seq, output_telo_fasta)
  
  plot_single_telo(x_length = max(max_length, length(current_seq_unlist)),
        seq_length = length(current_seq_unlist), subs =  analyze_list[[1]], 
        serial_num = current_serial, seq_start = start(telo_irange), 
        seq_end = end(telo_irange), save_it = TRUE, main_title = title,   
        w = 750, h = 300, output_jpegs = output_jpegs)
  
  plot_single_telo(x_length = length(current_seq_unlist), seq_length = 
                     length(current_seq_unlist), subs =  analyze_list[[1]], 
                   serial_num = current_serial,
                   seq_start = start(telo_irange), seq_end = end(telo_irange), 
                   save_it = TRUE, main_title = title,  w = 750, h = 300, 
                   output_jpegs = output_jpegs_1)
  
  # eps
  plot_single_telo(x_length = length(current_seq_unlist), seq_length = 
                     length(current_seq_unlist), subs =  analyze_list[[1]], 
                   serial_num = current_serial,
                   seq_start = start(telo_irange), seq_end = end(telo_irange), 
                   save_it = TRUE, main_title = title,  w = 750, h = 300, 
                   output_jpegs = output_jpegs_1, eps = TRUE)
  
  plot_single_telo_ggplot2(seq_length = length(current_seq_unlist),
                           subs = analyze_list[[1]], serial_num = current_serial,
                           seq_start = start(telo_irange), seq_end = end(telo_irange),
                           save_it = TRUE, main_title = title, output_jpegs)
  
  
  
  df <- df %>%
    add_row(Serial = current_serial, sequence_ID = current_fastq_name, 
    sequence_length = length(current_seq_unlist), telo_density = telo_density
    , Telomere_start = start(telo_irange), Telomere_end = 
    end(telo_irange), Telomere_length = width(telo_irange))
  
  return(df)
}


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
                 Telomere_length = integer())
  # add telo density : get_sub_density <- function(sub_irange, ranges){
  # For the fasta output of the reads which pass the filter
  #telomeres_dna_string_set <- DNAStringSet() 
  current_serial <- serial_start
  
  for (i in seq_along(sample_telomeres)) {
    #' now make corrections with mismax+indels patterns both for start and end 
    #' of the telomere
    #telo_irange <- gray_area_adges(telo_irange, current_seq)
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
                                            path = input_path), format = format)
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
  
  cl <- makeCluster(num_of_cores)
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
      patterns = dna_rc_patterns, min_density = global_min_density*0.8, mc.cores = num_of_cores) 
  
  
  #test_filter2 <- samps_1000[unlist(logical_100)]
  
  #logical_100 <-sapply(X = samp_100, FUN = filter_density, patterns = dna_rc_patterns, min_density = 0.175)
  
  
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

bool_question <- function(question) {
  if(assertthat::is.string(question) == FALSE) {
    stop("The input for bool_question should be a string!", call. = TRUE)
  }
  ans <- ""
  repeat{ 
    writeLines(text = question, con = "stdin")
    ans <- readLines(con = "stdin", n = 1 )
    if(ans %in% list("Yes", "yes", "No", "no")) {
      break
    } 
  }
  return (ans %in% list("yes", "Yes"))
}


################## Arguments ######################################################

PATTERNS_LIST <- list("CCCTRR", "CCTRRC", "CTRRCC", "TRRCCC", "RRCCCT", "RCCCTR")  

patterns_dna <- lapply(PATTERNS_LIST, DNAString)
dna_rc_patterns <- lapply(patterns_dna, Biostrings::reverseComplement)
dna_rc_patterns <- lapply(dna_rc_patterns, toString)


####################################################
# TODO
#' 2. Make a package
#' 3. upload to git
#' 4. Add a subtelo- only reads (Telomers are trimmed off)
#' 5. Give trim barcode option (length of nt to trimm..1300/30)
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
#' 19. fixed file printing: exclude the input_dir


# code for running script on linux shell - use log file for statistics
if (length(args) < 2) {
  stop("At least 2 argument must be supplied (input file/dir and output dir).n", call.=FALSE)
} else if (length(args) == 2) {
  # default output file
  args[3] = "fastq"
}


# log function blueprint : summary(sample), %telomeric_reads , summary(telo_read) ...
# I need to create a log funcion : with ifelse ( if telomeric patterns were found or not -> no one passed the filteration or df isempty ....)
# test log file
tmp <- file.path(args[2], "run.log")

# Open log
lf <- log_open(tmp)

t1 <- Sys.time()
log_print(base::paste("Work started at:", toString(t1)), hide_notes = TRUE) # Send message to log

log_print("The input files:", hide_notes = TRUE)
# add the names of the files which we analyze.

if(dir.exists(args[1])) {
  filepath = dir(full.names = TRUE, 
                 path = args[1])
  
  for(i in seq_along(filepath)) {
    log_print(filepath[i], hide_notes = TRUE)
  }
} else {
  log_print(args[1], hide_notes = TRUE)
}


dna_reads <- create_sample(input_path = args[1], format = args[3])
global_total_length <- length(dna_reads)
# Print data to log: length(sample), nrow(df), summary sts of read_leangth, Telo_length
log_print(base::paste("Total reads in sample:", toString(length(dna_reads)), "\n"), hide_notes = TRUE)
log_print("Summary statistics of the sample reads length:", hide_notes = TRUE)
log_print(summary(width(dna_reads)), hide_notes = TRUE)
log_print("\n",  hide_notes = TRUE)


rc <- bool_question("Use reverse complement? (print yes/Yes or no/No)")
filter <- bool_question("Use the filteration? (print yes/Yes or no/No)\nIt filters the reads according to length (>= 1000 nt) and the density at the edge of the read).")
if(filter) {
  right_edge <- bool_question("Check the right edge? (print yes/Yes or no/No)\n Be aware if you chose to use the reverse complement!")
  dna_reads <- filter_reads(samples = dna_reads, patterns = dna_rc_patterns, 
                            do_rc = rc, num_of_cores = 10, subread_width = 200, right_edge = right_edge)
} else {
  if(rc) {
    dna_reads <- reverseComplement(dna_reads)
  }
}



create_dirs(output_dir = args[2])

if(Sys.info()[1] == "Linux"){
  #now try parallel using future package
  # prev cod
  #df_summary <- search_patterns(sample_telomeres = dna_reads, pattern_list = dna_rc_patterns, output_dir = args[2], min_density = 0.3, serial_start = )
  
  
  # now try to use future::plan(multicore, workers = 10) , I need to set the save_summary to false, and after reduce save it.
  # also save the reads as 1 read (serian.fasta : 1.fasta ...) per file instead of reads.fasta.
  options(future.globals.maxSize = 1048576000*1.5) # 1.5 gigabyte max
  plan(multicore, workers = 10)
  
  # create indices for each worker using split 
  
  groups_length <- 10
  seq_over_length <- seq.int(from = 1, by = 1, length.out = length(dna_reads))
  groups <- 1:groups_length
  split_seq <- suppressWarnings(split(x = seq_over_length, f = groups))
  
  df1 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`1`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = 1 )
  df2 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`2`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(split_seq$`1`) + 1 )
  df3 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`3`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:2])) + 1)
  df4 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`4`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:3])) + 1)
  df5 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`5`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:4])) + 1)
  df6 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`6`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:5])) + 1)
  df7 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`7`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:6])) + 1)
  df8 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`8`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:7])) + 1)
  df9 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`9`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:8])) + 1)
  df10 %<-% search_patterns(sample_telomeres = dna_reads[split_seq$`10`], pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density, serial_start = length(unlist(split_seq[1:9])) + 1)
  
  df_summary <- Reduce(union_all , list(df1,df2,df3, df4, df5, df6, df7, df8, df9, df10))
  # end of parallel try
} else {
  df_summary <- search_patterns(sample_telomeres = dna_reads, pattern_list = dna_rc_patterns, output_dir = args[2], min_density = global_min_density)
}




# add row number
df_summary <- df_summary %>% 
  mutate(row_number = seq_len(nrow(df_summary)), .before = "Serial")



log_print(base::paste0("Numer of reads which identified as Telomeric: ", nrow(df_summary)), hide_notes = TRUE)
log_print(base::paste0("% of total reads: ", toString(round( (100*nrow(df_summary)) / global_total_length, 2 ) ), "%\n" ), hide_notes = TRUE)

# summary statistics of the Telomeric reads
log_print("Summary statistics for the Telomeric reads:", hide_notes = TRUE)
log_print("reads length:", hide_notes = TRUE)
log_print(summary(df_summary$sequence_length), hide_notes = TRUE)
log_print("Telomere length:", hide_notes = TRUE)
log_print(summary(df_summary$Telomere_length), hide_notes = TRUE)


write_csv(x = df_summary, file = file.path(args[2],"summary.csv"))
t2 <- Sys.time()


log_print(base::paste("Work ended at:", toString(t2)), hide_notes = TRUE) 
# Close log
log_close()
writeLines(readLines(lf))

