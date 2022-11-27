################################################################
#
# Author: Dan Lichtental
# Copyright (c) Dan Lichtental, 2022
# Email:  dan.lichtental@mail.huji.ac.il
# 
# Date: 2022-07-05
#
# Script Name: Telomere pattern finder
#
# Script Description: In this script we search over fastq file sequences for telomeric patterns.
#
#
# Notes:




library(BiocManager)
library('BSgenome')
library('stringr')
library(tidyverse)  
library(IRanges)
library(purrr)
# library(data.table)
library("parallel")
library(logr)




split_telo <- function(dna_seq, sub_length = 100){
  #' @title Splits a DNA sequence nto subsequences. 
  #' @description This function calculate the sequence ength and creates IRanges objects of subseuences of a given length
  #' if The dna_seq%%subseq != 0   and last' width < sub_length/2 Then we will remove this last index making the last subtelomere longer then by. (because we want to prevent a case where we have for 
  #' example subtelomere of length < sub_length/2 which is too short to consider -> then the last subtelomere is a bit longer)
  #' If the length of dna_seq is less then the sub_length it will return an empty Iranges Object.            
  #' @usage 
  #' @param dna_seq: DNAString object 
  #' @param sub_length: The length of each subsequence
  #' @value An IRanges object of the subsequences.
  #' @return iranges_idx: IRanges obsject of the indices of each subsequence
  #' @examples 
  idx_start <- seq(1,length(dna_seq), by=sub_length)
  idx_end <- idx_start + sub_length -1
  idx_end[length(idx_end)] <- length(dna_seq)
  if(length(dna_seq) - last(idx_start) < (sub_length/2)){ # if last subsequence is less then 50% 
    idx_start <- idx_start[1:length(idx_start)-1]
    idx_end <- idx_end[1:length(idx_end)-1]
    idx_end[length(idx_end)] <- length(dna_seq)
  }
  idx_ranges <- IRanges(start = idx_start, end = idx_end)
  return(idx_ranges)
} # IRanges object



# my improvment: fiding the IRanges and making union for overlaps and then calculate according to sum(width of the IRanges)
# The calculation is on the full sequence and it is not fit for subsequences
get_densityIRanges <- function(sequence, patterns){
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
  #' @value A numeric for the total density of the pattern(patterns) in the sequences, a IRanges object of the indcies of the patterns found.
  #' @return a tuple of (density, IRanges) Total density of all the patterns in the list( % of the patterns in this sequence) and the IRanges of them
  #' @examples 
  total_density <- 0 
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list 
  if(is.list(patterns)){
    patterns <- unique(patterns)  # make sure there are no dups
    for( pat in patterns){
      mp_all <- union.Vector(mp_all, matchPattern(pattern = pat, subject = unlist(sequence), max.mismatch = 0) )
    }
  }
  else{
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence), max.mismatch = 0)
    mp_all <- union.Vector(mp_all, mp_all) # incase there are overlaps
  }
  total_density <-sum(width(mp_all)) / nchar(sequence)
  return(list(total_density, mp_all))
}  


# with a data frame we can further explore the different patterns
# my improvment: fiding the IRanges and making union for overlaps and then calculate according to sum(width of the IRanges)
# The calculation is on the full sequence and it is not fit for subsequences
  get_densityIRanges_with_csv <- function(sequence, patterns, output_path = NA){
  #' @title Pattern searching function.
  #' @description: get the density of a given pattern or a total density of a list of patterns, and IRanges of the patterns.
  #' @param pattern: a list of patterns or a string of 1 pattern.
  #' @param sequence: DNAString object
  #' @param output_path: the filename for the csv patterns summary.
  #' @value A numeric for the total density of the pattern(patterns) in the sequences, a IRanges object of the indcies of the patterns found.
  #' @return a list of (density, IRanges, data frame) Total density of all the patterns in the list( % of the patterns in this sequence) and the IRanges of them
  #' @examples 
  total_density <- 0 
  patterns_df <- data_frame(patterns = character(), start_idx = integer(), end_idx = integer())
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list 
  if(is.list(patterns))
  {
    patterns <- unique(patterns)  # make sure there are no dups
    for( pat in patterns)
    {
      curr_match <- matchPattern(pattern = pat, subject = unlist(sequence), max.mismatch = 0)
      patterns_df <- add_row(patterns_df, patterns = as.data.frame(curr_match)$x, start_idx = start(curr_match), end_idx = end(curr_match))
      mp_all <- union.Vector(mp_all, curr_match )
    }
  }
  else
  {
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence), max.mismatch = 0)
    patterns_df <- add_row(patterns_df, patterns = as.data.frame(mp_all)$x, start_idx = start(mp_all), end_idx = end(mp_all))
    mp_all <- union.Vector(mp_all, mp_all) # incase there are overlaps
  }
  total_density <-sum(width(mp_all)) / nchar(sequence)
  pat_list <- list(total_density, mp_all, patterns_df)
  names(pat_list) <- c("total density", "patterns IRanges", "Patterns data frame")
  if(!is.na(output_path)) 
  {
    write_csv(x = patterns_df, file = output_path)
  }
  return(pat_list)
}  


get_sub_density <- function(sub_irange, ranges){
  #' @title Calculate density of a subsequnce.
  #' @details with a given IRanges of a subsequence and IRanges of patterns, compute the density of the IRanges within the 
  #'        subseuence range.
  #' @param sub_irange: the IRange of the subsequence
  #' @param ranges: The IRanges of the patterns found in the full sequence
  #' @value a numeric which is the density in range [0,1].  
  #' @return The density of the patterns in the subseuence according to the IRanges.
  #' @description  sub_irange = (10, 30), ranges = {(2,8), (16,21), (29,56)} -> the intersect is {(16,21), (29, 30)} -> 
  #                width = 6+2 = 8 , sub_irange width = 21 -> density = 8/21
  #' @examples sub_irange <- IRanges(start = 10, end = 30)
  #'           ranges <- IRanges(start = c(2,16,29), end = c(8,21,56))
  #'           get_sub_density(sub_irange =  sub_irange, ranges = ranges) # 0.3809524
  # this compute the desity of iranges of patterns within a given irange of a subsequence
  return( sum(width( intersect.Vector(sub_irange, ranges))) / width(sub_irange) ) 
}


###########3 My cahnge from prev - return a list 0f (df, total_density) ###############
analyze_subtelos <- function(dna_seq, patterns , sub_length = 100, MIN_DENSITY = 0.3){ # return list(subtelos, list_density_mp[1])
  #' @title Analyze the patterns for each subsequence.
  #' @description s split a dna sequence to subsequences and calculate the density of each subsequence
  #' @param dna_seq: a dna sequence (DNAString object)
  #' @param patterns: a list of patterns or a string of 1 pattern
  #' @param sub_length: The length of the subsequences for split_telo fuction.s
  #' @value a list(data frame of all subtelomeres and their properties ,numeric for total density)
  #' @return  a list of (a data frame, list(numeric: total density, IRanges for patterns) 
  #' @examples 
  
  
  # aother option is to create 5 vector s and then make a data.table from them
  
  # density and iranges of matchPattern
  list_density_mp <-get_densityIRanges(dna_seq, patterns = patterns)
  #density <- list_density_mp[[1]]
  mp_iranges <- list_density_mp[[2]]
  # get start indexes of "subtelo"
  idx_iranges <- split_telo(dna_seq, sub_length = sub_length) 
  # create empty dataframe which will contain all subtelomeres and their properties
  subtelos <- data.frame(ID = as.integer(), start_index = as.integer(), end_index = as.integer(), density = as.numeric(), class = as.numeric() )
  cur.ID <- 1             # intitialize ID counter
  # loop through start indexs
  for(i in 1:length(idx_iranges)){
    #first.idx <- start(idxs_iranges[i])
    #last.idx <- end(idxs_iranges)
    subtelo.density <-get_sub_density(sub_irange= idx_iranges[i] , ranges = mp_iranges)
    #cur.seq <- subseq(dna_seq, start = first.idx, end = last.idx )  #instead_of  substr(dna.seq, first.idx, last.idx)         # sequence of the current subtelomere
    # subtelo.density <- get_density(cur.seq, patterns = patterns)        # calculate all density and classify each subtelo for the patern "CCCTRR"
    CLASSES <- list('CCCTAA'=-5, 'NONE'=1, 'SKIP'=0)
    
    
    ###############3 NEED TO CHANGE FOR AN ARGUMENT OF CLASSES #####################
    subtelo.class <- CLASSES$CCCTAA
    if(subtelo.density < MIN_DENSITY){
      if(subtelo.density < 0.1){
        subtelo.class <- CLASSES$SKIP
      }
      else{
        subtelo.class <- CLASSES$NONE
      }
    }
    #subtelos <- rbindlist(list(subtelos, list(cur.ID, start(idx_iranges[i]) , end(idx_iranges[i]), subtelo.density, subtelo.class)))
    subtelos <- subtelos %>%
      add_row(ID = cur.ID, start_index =start(idx_iranges[i]) , end_index = end(idx_iranges[i]), density = subtelo.density, class = subtelo.class)
    cur.id <- cur.ID + 1
  }
  return( list(subtelos, list_density_mp) ) # return the subtelos df and the list(total density, mp_iranges)
}
#  a a list of (a data frame, numeric: total density) 


# I need to put The CLASSES as an arument ?
find_telo_position <- function(seq_length, subtelos, min_in_a_row = 3, min_density_score = 2){ # 15,10 for sub_length == 20 
  #' @title: Find the position of the Telomere(subsequence) within the seuence.
  #' @description:  Find the start and end indices of the subsequence within the sequence according to a data frame.
  #' @usage 
  #' @param seq_length: The length of the read.
  #' @param subtelos: data frame of a subseuences indices, density and class (from the analyze_subtelos)
  #' @value An IRanges object of length 1.
  #' @return (start, end) irange  of the Telomere.
  #' @examples
  
  ###############3 NEED TO CHANGE FOR AN ARGUMENT OF CLASSES #####################
  CLASSES <- list('CCCTAA'=-5, 'NONE'=1, 'SKIP'=0)  # # we have a problem in this function: Error in if (subt$class == CLASSES$SKIP | subt$class == CLASSES$NONE |  : argument is of length zero
  # set score, start, in.a.row to 0,-1,0
  
  score <- 0.0
  start <- -1
  end <- -1
  in_a_row <- 0
  start_end_diff <- subtelos[1, "end_index"] - subtelos[1, "start_index"] 
  
  # loop through subsequences
  end_position <- 0 # for end loop
  for (i in 1:nrow(subtelos)){ # CLASSES <- list('CCCTAA'=-5, 'NONE'=1, 'SKIP'=0)
    subt <- subtelos[i,]
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == CLASSES$SKIP | subt$class == CLASSES$NONE | is.na(subt$class)){
      score <- 0
      start <- -1
      in_a_row <- 0
    }
    else{
      # otherwise add one to in.a.row, update score and set start index
      in_a_row <- in_a_row + 1
      score <- score + subt$density
      if (start == -1){
        start <- subt$start_index
      }
    }
    # if more than MIN.IN.A.ROW subtelomeres were found and the overall score is high enough, return start index
    if(in_a_row >= min_in_a_row && score >= min_density_score){
      j <- i+1
      end_position <- i+1
      break
    }
  }
  if(end_position == 0){
    return(IRanges(1,1)) # no telomere was found 
  }
  
  
  
  # search for end from the last subsequence (backward to finding stat incase there is island of non-telomeric subsequence)
  end <- -1 
  score <- 0.0
  in_a_row <- 0
  for (i in nrow(subtelos):end_position){ # CLASSES <- list('CCCTAA'=-5, 'NONE'=1, 'SKIP'=0)
    subt <- subtelos[i,]
    # if the subsequence's class is SKIP, NONE or NA, reset values
    if (subt$class == CLASSES$SKIP || subt$class == CLASSES$NONE || is.na(subt$class) ){
      score <- 0.0
      end <- -1
      in_a_row <- 0
    }
    else{
      # otherwise add one to in.a.row, update score and set start index
      in_a_row <- in_a_row + 1
      score <- score + subt$density
      if (end == -1){
        end <- subt$end_index
      }
    }
    # if more than MIN.IN.A.ROW subtelomeres were found and the overall score is high enough, return start index
    if(in_a_row >= min_in_a_row && score >= min_density_score){
      break
    }
  }
  
  
  
  if( start > end){
    end <- start + start_end_diff
    #browser()  ###### https://support.rstudio.com/hc/en-us/articles/205612627-Debugging-with-the-RStudio-IDE
  }
  
  return(IRanges(start = start, end = end))
}

# check the plot
# I need to adjust the plot to my code
# dna.seq is a DNAStringSet of length 1
# This is version updated at 6.01.2022
plot_single_telo <- function(x_length, seq_length, subs, serial_num, seq_start, seq_end,save.it=T, main_title = "", w=750, h=300, OUTPUT_JPEGS){ # add OUTPUT_JPEGS as arg
  #' @title plot the density over a sequence
  #' @param x_length: The length of the x axis.
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos function
  #' @param serial_num: The serial number of the current sequence, used as the name of the file
  #' @param seq_start: The start of the Telomere.
  #' @param seq_end: The end " ".
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param OUTPUT_JPEGS: the output directory for saving the file
  subs <- na.omit(subs)
  # save file if specified
  if(save.it){
    eps_path <- paste(OUTPUT_JPEGS, paste('read', serial_num, '.eps',sep=''), sep='/')  
    setEPS()
    # naming the eps file
    postscript(eps_path)
                                                  
  }
  
  
  # 26-07: my addition: save the csv file subs
  # write_csv(x = subs, file = paste(OUTPUT_JPEGS, paste('read', serial_num, '.csv',sep=''), sep='/') )
  
  
  # give extra x for the legend at the topRigth
  plot(subs$density ~ subs$start_index, type='n', yaxt='n', xaxt='n',ylab='', xlab='', ylim=c(0,1), xlim=c(1,x_length + round(x_length/4.15)) ) 
  # create axes
  xpos <- seq(1, x_length, by=1000) # I have cahnged from 0 to 1
  axis(1, at=xpos, labels=sprintf("%.1fkb", xpos/1000)); title(xlab = "Position", adj = 0)
  axis(2, at=seq(-0.1,1,by=0.1), las=2)
  # add polygon to plot for each variant repeat. 
  # mychange: only comp_ttaggg 
  polygon(y=c(0,subs$density,0), x=c(1,subs$start_index,seq_length), col=rgb(1,0,0,0.5), lwd=0.5) # change c(1,) instead of c(0, ) for x
  
  rect(xleft = seq_start, ybottom = -0.1, xright = seq_end, ytop = 0, col = "red") # 
  rect(xleft = seq_end+1, ybottom = -0.1, xright = seq_length, ytop = 0, col = "blue")
  if(seq_start > 1){
    rect(xleft = 1, ybottom = -0.1, xright = seq_start, ytop = 0, col = "blue")
  }
  
  abline(h=1, col="black", lty = 2)
  abline(h=0, col="black", lty = 2)
  legend(x = x_length, y = 1, legend=c("telomere", "sub-telomere"),col=c("red", "blue"), lty=1, lwd= 2,cex=1.2)
  sub_title <- paste("read length:", seq_length, ", telomere length:", abs(seq_start-seq_end)+1)
  title( main = main_title, sub = sub_title, ylab='Density')
  dev.off()
  
}


plot_single_telo_ggplot2 <- function( seq_length, subs, serial_num, seq_start, seq_end,save.it=T, main_title = "", OUTPUT_JPEGS){ # add OUTPUT_JPEGS as arg
  #' @title plot the density over a sequence
  #' @param seq_length: The length of the sequence
  #' @param subs: the data frame of subseuences from the analyze_subtelos function
  #' @param serial_num: The serial number of the current sequence, used as the name of the file
  #' @param seq_start: The start of the Telomere.
  #' @param seq_end: The end " ".
  #' @param sava.it: save the file if TRE
  #' @param main_title: ad a title.
  #' @param w: width of the jpeg
  #' @param h: height of the jpeg
  #' @param OUTPUT_JPEGS: the output directory for saving the file
  
  subs <- na.omit(subs)
  # save the densities in a csv file for later use !
  write_csv(x = subs, file = paste(OUTPUT_JPEGS, paste( serial_num, '.csv',sep=''), sep='/') )
  sub_title <- paste("read length:", seq_length, ", telomere length:", abs(seq_start-seq_end)+1)
  
  my_ggplot <- ggplot2::ggplot(subs, aes(x = start_index, y = density)) +
    geom_area(color = NA, fill = "red") +
    # the 3 last lines are for the telomere-sub-telomre region ....
    geom_rect(aes(xmin = seq_start, ymin= -0.01, xmax= seq_end, ymax = 0), color = "green", fill = "green")  +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1,linetype = "dashed") +
    labs(title = main_title, x = "Position", y = "Density", subtitle = sub_title)
  
  
  if(seq_end < seq_length)
  {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = seq_end , ymin= -0.01, xmax= seq_length, ymax = 0), color = "blue", fill = "blue")
  }
  
  if(seq_start > 1)
  {
    my_ggplot <- my_ggplot +
      geom_rect(aes(xmin = 1, ymin= -0.01, xmax= seq_start, ymax = 0), color = "blue", fill = "blue")
    
  }
  
  
  # save file if specified
  if(save.it){
    eps_path <- paste(OUTPUT_JPEGS, paste('gg_read', serial_num, '.eps',sep=''), sep='/')  
    ggsave(plot = my_ggplot, device = "eps", filename = eps_path)  
  }
  
}



############## Running functions #############################################################  


# removed the telorrete
# removed the if condition(filter), all sequences input are already came passing the filter 
# added a serial_start _ to work with:
searchPatterns <- function(sample_telomeres , pattern_list, max_length = 1e5, csv_name = "summary",output_dir,serial_start = 1, min_density,  title = "Telomeric repeat density"){
  #'@title Search given Patterns over a DNA sequences.
  #'@param sample_telomeres: the DNAString Set of the reads.
  #'@param pattern_list: a list of patterns or a string of 1 pattern.
  #'@param max_length: The x-axis length for the plot.
  #'@param csv_name: The name of the csv file
  #'@param output_dir: 
  #'@param serial_start: The first id serial number.
  #'@param min_density: The minimal density of the patterns in a sequence to be consider relevant.
  #'@param title: The title for the density plots.
  #'@param return: Returns a data frame of the reads.
  
  if(!dir.exists(output_dir)){ # update  did it 
    dir.create(output_dir)
  }
  
  OUTPUT_TELO_CSV <- paste(output_dir, paste(csv_name, 'csv', sep='.'), sep='/')
  OUTPUT_TELO_FASTA <- paste(output_dir, paste("reads", 'fasta', sep='.'), sep='/')
  OUTPUT_JPEGS <- paste(output_dir, 'single_read_plots', sep='/')
  dir.create(OUTPUT_JPEGS)
  OUTPUT_JPEGS.1 <- paste(output_dir, 'single_read_plots_adj', sep='/')
  dir.create(OUTPUT_JPEGS.1)
  
  #max_length <- max(width(sample_telomeres))
  #  I HAVE ADDED TLOMERE LENGTH, START @ END
  df<-data.frame(Serial = integer(), sequence_ID = character(), sequence_length = integer(), telo_density = double()
                 , Telomere_start = integer(), Telomere_end = integer(), Telomere_length = integer())
  # add telo density : get_sub_density <- function(sub_irange, ranges){
  LargeDNAStringSet <- DNAStringSet() # For the fasta output of the reads which pass the filter
  current_serial <- serial_start
  
  for( i in 1:length(sample_telomeres) ){
    current_fastq_name <- names(sample_telomeres[i])
    current_seq <- unlist(sample_telomeres[i])
    # we skip the adaptor and telorete so we start with base 57
    
    
    #################### NEED TO MAKE IT MORE RUBUST , MAYBE ADD THE FILTER FUNCTION AS AN INPUT ########################
    #current_filt_100 <-  get_densityIRanges( subseq(current_seq, start = length(current_seq)-99, end = length(current_seq) ), patterns = pattern_list) 
    
    # # returns a a list of (a data frame, list(numeric: total density,iranges)) 
    analyze_list <- analyze_subtelos(dna_seq = current_seq , patterns =  pattern_list, MIN_DENSITY = min_density)
    telo_irange <- find_telo_position(seq_length = length(current_seq), subtelos = analyze_list[[1]], min_in_a_row = 3, min_density_score = 2 )
    
    
    irange_telo <- analyze_list[[2]][[2]]
    if(width(telo_irange) < 100 ) {next} # not considered a Telomere
    s_index <- start(telo_irange) 
    # make the strat/end more accurate (usethe IRanges for the patterns)
    iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) 
    if(length(iranges_start) > 0){ start(telo_irange) <- start(irange_telo[iranges_start[1]])} 
    # try more accurate: take the max of which and also check new_end >= new_start before updating the IRange object
    e_index <- end(telo_irange) 
    iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
    if(length(iranges_end) >0 ) {
      #end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
      new_end <- end(irange_telo[iranges_end[length(iranges_end)]])# take the last pattern in range 
      if(new_end >= start(telo_irange)){  end(telo_irange) <- new_end} # make sure end >= start     
    }  
    
    # Onc ew have the Telomere indices calculate the density of the patterns within it's range.
    telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
    
    if(telo_density < 0.5)
    {
      telo_irange <- find_telo_position(seq_length = length(current_seq), subtelos = analyze_list[[1]], min_in_a_row = 10, min_density_score = 6 )
      irange_telo <- analyze_list[[2]][[2]]
      if(width(telo_irange) < 100 ) {next} # not considered a Telomere
      s_index <- start(telo_irange) 
      # make the strat/end more accurate (usethe IRanges for the patterns)
      iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) 
      if(length(iranges_start) > 0){ start(telo_irange) <- start(irange_telo[iranges_start[1]])} 
      # try more accurate: take the max of which and also check new_end >= new_start before updating the IRange object
      e_index <- end(telo_irange) 
      iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
      if(length(iranges_end) >0 ) {
        #end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
        new_end <- end(irange_telo[iranges_end[length(iranges_end)]])# take the last pattern in range 
        if(new_end >= start(telo_irange)){  end(telo_irange) <- new_end} # make sure end >= start     
      }  
      
      # Onc ew have the Telomere indices calculate the density of the patterns within it's range.
      telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
    }
    
    # now make corrections with mismax+indels patterns both for start and end of the telomere
    #telo_irange <- gray_area_adges(telo_irange, current_seq)
    
    df <- df %>%
      add_row(Serial = current_serial, sequence_ID = current_fastq_name, sequence_length = length(current_seq), 
              telo_density = telo_density, Telomere_start = start(telo_irange), Telomere_end = end(telo_irange), Telomere_length = width(telo_irange))
    
    if(max_length < length(current_seq)){
      max_length <- current_seq
    }
    
    plot_single_telo(x_length =max(max_length, length(current_seq)), seq_length = length(current_seq), subs =  analyze_list[[1]], serial_num = current_serial ,
                     seq_start = start(telo_irange),seq_end = end(telo_irange), save.it=T, main_title = title,  w=750, h=300, OUTPUT_JPEGS= OUTPUT_JPEGS)
    plot_single_telo(x_length = length(current_seq), seq_length = length(current_seq), subs =  analyze_list[[1]], serial_num = current_serial ,
                     seq_start = start(telo_irange),seq_end = end(telo_irange), save.it=T, main_title = title,  w=750, h=300, OUTPUT_JPEGS= OUTPUT_JPEGS.1)
    plot_single_telo_ggplot2(seq_length = length(current_seq),subs =  analyze_list[[1]], serial_num = current_serial ,
                             seq_start = start(telo_irange),seq_end = end(telo_irange), save.it=T, main_title = title,  OUTPUT_JPEGS= OUTPUT_JPEGS )
    
    LargeDNAStringSet <- append(LargeDNAStringSet, values = sample_telomeres[i])
    current_serial <- current_serial + 1
    
    
  } # end of for loop
  
  # need to save the df in a file
  write.csv(x=df, file=OUTPUT_TELO_CSV, col.names = F)
  writeXStringSet(LargeDNAStringSet, OUTPUT_TELO_FASTA)
  message("Done!") #  now what's left is to extract the sequences from fasta to fasta using the read_names list file ( 3 files )
  
  return(df)
  
} # end of the function searchPatterns




# need to create parApply for filter 
library("parallel")
filter_density <- function(sequence, patterns, min_density = 0.18){
  #'
  current_seq <- unlist(sequence)
  total_density <- 0 
  mp_all <- IRanges()# union of all the IRanges of all the patterns in the list 
  if(is.list(patterns)){
    patterns <- unique(patterns)  # make sure there are no dups
    for( pat in patterns){
      mp_all <- union.Vector(mp_all, matchPattern(pattern = pat, subject = unlist(sequence), max.mismatch = 0) )
    }
  }
  else{
    mp_all <- matchPattern(pattern = patterns, subject = unlist(sequence), max.mismatch = 0)
    mp_all <- union.Vector(mp_all, mp_all) # incase there are overlaps
  }
  total_density <-sum(width(mp_all)) / nchar(sequence)
  
  return( total_density >= min_density)
  
}


# need to correct the spelling for telorette
run_with_rc_and_filter <- function(samples,  patterns, output_dir,  do_rc = TRUE, serial_start = 1){
  #' @title: Run a search for Telomeric sequences on the reads   
  #' @description use rc to adjust for the patterns and barcode/telorette locatio ( should be at the last ~ 60-70 bases), filter out reads with no 
  #'              no telomeric pattern at the edge and run searchPatterns_withTelorette  
  #' @usage 
  #' @param samples: A DNAStringSet of reads
  #' @param patterns: The patterns for the telomere
  #' @param output_dir: The path for the output directory
  #' @param telorrete_pattern: The telorette pattern
  #' @return
  #' @examples
  
  
  
  if(!dir.exists(output_dir)){ # update  did it 
    dir.create(output_dir)
  }
  
  # create a log file
  tmp <- file.path(output_dir, 'run_summary.log')
  
  samps_1000 <- samples[width(samples) >= 1e3]
  if(do_rc){samps_1000 <- Biostrings::reverseComplement(samps_1000)}
  
  copies_of_r <- 10
  
  cl <- makeCluster(copies_of_r)
  samp_100 <- parLapply(cl, samps_1000, subseq, end = -67, width = 100)  # change to -(61+ just incase there are indels ( barcode_telorette == 61))
  stopCluster(cl)
  
  cl <- makeCluster(copies_of_r)
  logical_100 <- parSapply(cl, samp_100,  filter_density,patterns = patterns , min_density = 0.175)
  stopCluster(cl)
  names(logical_100) <- NULL
  
  samps_filtered <- samps_1000[logical_100]
  
  if(length(samps_filtered) < 1)
  {
    message('No read have passed the filteration at run_with_rc_and_filter function!')
  }else
  {
    searchPatterns(samps_filtered, pattern_list = patterns, max_length = 
                     max(max(width(samps_filtered)), 150000), output_dir = 
                     output_dir , min_density = 0.3, serial_start = serial_start )

  }
  
  
  
}


create_sample <- function(input_path, format = 'fastq')
{
  #' @title: Extract DA sampe from fasta/q files.
  #' @description By given input files(fastq otr fasta) creates a DNAStringSet object.
  #' @param input_path: path to the file or directory containing files.
  #' @param format: The file/files format should be either fastq format or fasta format
  #'                gz extension is supported.
  if(dir.exists(input_path))
  {
    sample <- Biostrings::readDNAStringSet(filepath = dir(full.names = T, path = input_path) , format = format)
  }else # it is a single file path
  {
    sample <- Biostrings::readDNAStringSet(filepath = input_path, format = format)
  }
  return(sample)
}




################## Arguments ######################################################

PATTERNS_LIST <- list("CCCTAA", "CCCTAG", "CCCTGA", "CCCTGG")  # CCCTRR
PATTERNS_LIST <- append(PATTERNS_LIST, list("CCTAAC", "CCTAGC", "CCTGAC", "CCTGGC")) # add CCTRRC
PATTERNS_LIST <- append(PATTERNS_LIST, list("CTAACC", "CTAGCC", "CTGACC", "CTGGCC")) # add CTRRCC
PATTERNS_LIST <- append(PATTERNS_LIST, list("TAACCC", "TAGCCC", "TGACCC", "TGGCCC")) # add TRRCCC
PATTERNS_LIST <- append(PATTERNS_LIST, list("AACCCT", "AGCCCT", "GACCCT", "GGCCCT")) # add RRCCCT
PATTERNS_LIST <- append(PATTERNS_LIST,list("ACCCTA", "GCCCTA", "ACCCTG", "GCCCTG"))  # add RCCCTR

patterns_dna <- lapply(PATTERNS_LIST, DNAString)
dna_rc_patterns <- lapply(patterns_dna, Biostrings::reverseComplement)
dna_rc_patterns <- lapply(dna_rc_patterns, toString)


####################################################

# I need to create a log funcion : with ifelse ( if telomeric patterns were found or not -> no one passed the filteration or df isempty ....)
sample_name <- "test1"
# test log file
tmp <- file.path('/media/lab/E/', "test.log")

# Open log
lf <- log_open(tmp)

# Send message to log
log_print(paste("Summary statistics of sample", sample_name))

# Perform operations
df1 <- read_csv('/home/lab/Downloads/Telomers/Trial7B/barcode03/summary.csv', show_col_types = F)
df1 <- select(df1, -c(1))
write.csv(x=df1, file= file.path('/media/lab/E/', "df.csv") , row.names = F)
# Print data to log: length(sample), nrow(df), summary sts of read_leangth, Telo_length
log_print()
log_print()
log_print()
# Close log
log_close()

# View results
writeLines(readLines(lf))










# todo: 
#' 1. Create a function for running from a filepath - Done: create_sample
#' 2. Make a package
#' 3. upload to git
#' 4. Add a subtelo- only reads (Telomers are trimmed off)
#' 5. Give trimm barcode option (length of nt to trimm..1300/30)
#' 6. Divide run_with_rc_and_filter to 2 function: one which do the initial filteration
#'    and a second one which do all the running.
#' 7. Create a csv files with a description...  a log file : https://cran.r-project.org/web/packages/logr/vignettes/logr.html  

# 24.11 - run on Trial15 unclassified
maindir <- '/home/lab/Downloads/Telomers/Trial15_bc1_2_3_6Telorette3Musmusculus/fastq_pass/rebarcoding_unclassified'
path_list <-dir(path = maindir)
path_list <- path_list[1:13]
path_list <- paste( maindir, path_list,sep = '/')
# I have a bug in  path_list[[4]]
for(i in 1:length(path_list))
{
  run_with_rc_and_filter(samples = create_sample(input_path = path_list[[i]]), 
                         patterns = dna_rc_patterns, output_dir = paste(path_list[[i]],'output', sep = '/')
                         , do_rc = T, serial_start = 1)
}

sample_bc04 <- create_sample(path_list[[4]])
sample_bc04 <- create_sample


# 24.11.2022 - rerun Trial 17 unclassified: 
maindir <- '/media/lab/E/Telomeres/Trial 17/Unclassified_rebarcoding-20221123T131555Z-001/Unclassified_rebarcoding/rebarcoding_output_Threshold_45'
path_list <- dir()
path_list <-path_list[c(1:12,18)]

path_list <- paste( maindir, path_list,sep = '/')


# 24.11 - run Trial 15 unclassified....
dir(full.names = T, path = path_list[[1]])
for(i in 7:13)
{
  run_with_rc_and_filter(samples = create_sample(input_path = path_list[[i]]), 
                         patterns = dna_rc_patterns, output_dir = paste(path_list[[i]],'output', sep = '/')
                         , do_rc = T, serial_start = 1)
}

# now run threshold 50
maindir <- '/media/lab/E/Telomeres/Trial 17/Unclassified_rebarcoding-20221123T131555Z-001/Unclassified_rebarcoding/rebarcoding_output_Threshold_50'
path_list <- dir(path = maindir)
path_list <-path_list[c(6:11,15)]

path_list <- paste( maindir, path_list,sep = '/')


# 24.11 - run Trial 15 unclassified....thteshold 50
for(i in 1:length(path_list))
{
  run_with_rc_and_filter(samples = create_sample(input_path = path_list[[i]]), 
                         patterns = dna_rc_patterns, output_dir = paste(path_list[[i]],'output', sep = '/')
                         , do_rc = T, serial_start = 1)
}






# 21.11.2022: Explore the T2T reference wrt Telomere start, end indexes
# upload the minimap2 results to compare

# minimap2 mapping results for trial 17





# Heads in T2T
heads_t2t <- readDNAStringSet(filepath = '/media/lab/E/Human  T2T-CHM13 Reference/HeadsAnalysis/reads.fasta')
heads_csv <- read_csv(file = '/media/lab/E/Human  T2T-CHM13 Reference/HeadsAnalysis/Heads.csv',show_col_types = F)

manual_corrected_heads_csv <- data.frame(sequence_ID = character(), Telomere_start = integer(), Telomere_end = integer())
#chr 1: 1,2705
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "1", Telomere_start = 1 , Telomere_end = 2705)

manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "2", Telomere_start = 1 , Telomere_end = 3617 )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "3", Telomere_start = 1 , Telomere_end = 2640)
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "4", Telomere_start = 1 , Telomere_end = 3271)
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "5", Telomere_start = 1 , Telomere_end = 2296)
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "6", Telomere_start = 1 , Telomere_end = 2899)
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "7", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "8", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "9", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "10", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "11", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "12", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "13", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "14", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "15", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "16", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "17", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "18", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "19", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "20", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "21", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "22", Telomere_start = 1 , Telomere_end = )
manual_corrected_heads_csv <- add_row(manual_corrected_heads_csv, sequence_ID = "23", Telomere_start = 1 , Telomere_end = )

#chr21 : 1,3015

# now make corrections with mismax+indels patterns both for start and end of the telomere
telo_irange <- gray_area_adges(telo_irange, current_seq, patterns)
{
  start <- start(telo_irange)
  end <- end(telo_irange)
}

# search at the right edge
end1 <- 2700
m1 <- Biostrings::matchPattern( pattern = "CCCTAA", subject = unlist(subseq(x = heads_t2t[1], start = end1 -5, end= end1 + 5)), max.mismatch = 2 , with.indels = T)
while(length(m1) > 0 && end1 < ( width(heads_t2t[1]) -50 ) )
{
  print(m1)
  end1 <- end1 + max(end(m1))
  print(end1)
  m1 <- Biostrings::matchPattern( pattern = "CCCTAA", subject = unlist(subseq(x = heads_t2t[1], start = end1 +1, end= end1 + 11)), max.mismatch = 2 , with.indels = T)
}

# search at the left edge
start1 <- 1101






# Running example
# I need to add jpeg file, improve creating csv , and another function which gets the dir of fastq files - and make a package
# 23.10.2022 - Trial 17
setwd("/media/lab/E/Telomeres/Trial 17/barcode02")
bc02_filepath <- dir()
t17_bc02 <- Biostrings::readDNAStringSet(filepath = bc02_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc02, patterns = dna_rc_patterns, output_dir = "bc02_output", do_rc = T)
rm(t17_bc02)


setwd("/media/lab/E/Telomeres/Trial 17/barcode01")
bc01_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc01_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)


setwd("/media/lab/E/Telomeres/Trial 17/barcode11")
bc011_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc011_filepath, format = "fastq")
t1 <- Sys.time()
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

setwd("/media/lab/E/Telomeres/Trial 17/barcode07")
bc07_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc07_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

setwd("/media/lab/E/Telomeres/Trial 17/barcode04")
bc04_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc04_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

setwd("/media/lab/E/Telomeres/Trial 17/barcode09")
bc09_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc09_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

setwd("/media/lab/E/Telomeres/Trial 17/barcode10")
bc010_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc010_filepath , format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

setwd("/media/lab/E/Telomeres/Trial 17/barcode08")
bc08_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = bc08_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)

Sys.time() - t1

# unclasiffied including miss classified - לא עבד, כנראה בגלל הזיכרון- כגודל הקבצים יותר מ- 4 גיגה..
t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied")
un_filepath <- dir(pattern =  "*.gz", recursive = T)
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T)
rm(t17_bc)
Sys.time() - t1



t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/1")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 1)
rm(t17_bc)
Sys.time() - t1

df1 <- read_csv("output/summary.csv")
# check serial end of prev before running
# check also the cores
ncores <- detectCores(logical = F)


t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/2")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 247)
rm(t17_bc)
Sys.time() - t1

df2 <- read_csv("output/summary.csv")


t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/3")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 510)
rm(t17_bc)
Sys.time() - t1

df3 <- read_csv("output/summary.csv")
max(df3$Serial)

t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/4")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 807)
rm(t17_bc)
Sys.time() - t1

df4 <- read_csv("output/summary.csv")
max(df4$Serial)

t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/5")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start= 1094)
rm(t17_bc)
Sys.time() - t1


df5 <- read_csv("output/summary.csv")
max(df5$Serial)

t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/6")
un_filepath <- dir()
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 1390)
rm(t17_bc)
Sys.time() - t1


df6 <- read_csv("output/summary.csv")
max(df6$Serial)

t1 <- Sys.time()
setwd("/media/lab/E/Telomeres/Trial 17/miss_clasiffied/unclassified/7")
un_filepath <- dir(pattern =  "*.gz", recursive = T)
t17_bc <- Biostrings::readDNAStringSet(filepath = un_filepath, format = "fastq")
run_with_rc_and_filter(samples = t17_bc, patterns = dna_rc_patterns, output_dir = "output", do_rc = T, serial_start = 1736)
rm(t17_bc)
Sys.time() - t1

df7 <- read_csv("output/summary.csv")
max(df7$Serial)
min(df7$Serial)



# megre all df's 1:7 ( 1902 reads)
df_list <- list(df1, df2, df3, df4, df5, df6,df7)
df_all <- bind_rows(df_list)
write_csv(x = df_all, 'summary_all.csv')

# merge all reads.fasta (1902 reads)
reads <- dir(pattern = "*.fasta")
all_reads <- Biostrings::readDNAStringSet(filepath = reads)
Biostrings::writeXStringSet(x = all_reads, filepath = 'all_reads.fasta')




# 31.08.2022
# comparisson of the reads to the reference ( buth tail and head)
# take the 30K sub-sequence
# for head( smooth sub-telomere ): chr 17  and Human read 63
# for tail (rough sub-telomere ): chr 10 and read 84

reads_t13 <-  readDNAStringSet(filepath = "/home/lab/Downloads/Telomers/Trial13hNL76telorettes/output_fastq_pass/reads.fasta")
# I will cut the 50 last bases of each read ( telorette + barcode)
# 63,728
read_84 <- reads_t13[84]
read_84_trimm_last50 <- Biostrings::subseq(x = read_84, start = 1, end = width(read_84) - 50)
read_84_30K <- Biostrings::subseq(x = read_84_trimm_last50, end = -1, width = 30000)
tails_human <- readDNAStringSet(filepath = "/media/lab/E/Human  T2T-CHM13 Reference/TailsAnalysis/reads.fasta")
tail10 <- tails_human[10] # The telomere length is 3140
# take 26K last bases
tail10_26K <- Biostrings::subseq(x = tail10, end = -1, width = 26000)

compare_tails <- append(x = read_84_30K, tail10_26K)

searchPatterns(sample_telomeres = compare_tails, pattern_list = dna_rc_patterns, max_length = 30000, csv_name = "readCompareTail10.csv",output_dir = "comparisson", min_density = 0.3) 




# 77,367, telo length 7181
# use rc , extract out the telorette , and take the first 30K
read_63 <- reads_t13[63]
read_63_trimm_last50 <- Biostrings::subseq(x = read_63, start = 1, end = width(read_63) - 50)
read_63_trimm_last50_rc <- Biostrings::reverseComplement(read_63_trimm_last50)
read_63_30K <- Biostrings::subseq(x = read_63_trimm_last50_rc, start = 1,  width = 30000 ) 
heads_human <- readDNAStringSet(filepath = "/media/lab/E/Human  T2T-CHM13 Reference/HeadsAnalysis/reads.fasta")
head17 <- heads_human[17] # 1690
head17_24500 <- Biostrings::subseq(x = head17, start = 1, width = 24500)

compare_heads <- append(x = read_63_30K, head17_24500)

searchPatterns(sample_telomeres = compare_heads, pattern_list = PATTERNS_LIST, max_length = 30000, csv_name = "readComparehead1",output_dir = "comparisson4", min_density = 0.2) 

# a 150K fisrt and 150K last nt of each chromosome
samples <- Biostrings::readDNAStringSet(filepath = "/media/lab/E/Human  T2T-CHM13 Reference/T2T_CHM13_150KchrHeadTail_Yexcluded.fasta")
# since we use the RRAGGG pattern we expect that only the tails will be idntfied as telomeric
run_with_rc_and_filter(samples =samples, patterns = dna_rc_patterns,output_dir = "Ref", do_rc = F)

# since we use the CCCTYY pattern we expect that only the heads will be idntfied as telomeric
run_with_rc_and_filter(samples =samples, patterns = PATTERNS_LIST,output_dir = "Ref", do_rc = F)


# let's test with ggplot2 gemo_area








# 12.09.2022 - test plot saved as eps, and another option using ggplot2
dna_samples <- Biostrings::readDNAStringSet("/home/lab/Downloads/Telomers/Trial13hNL76telorettes/output_fastq_pass/reads.fasta")

dna_samples <- dna_samples[c(12,13,37)]


searchPatterns(sample_telomeres = dna_samples, pattern_list = dna_rc_patterns, max_length = 140000, output_dir = "/home/lab/test_plots", min_density = 0.3) 


# check how many counts for each pattern
# 1. create a data frame according to pattern list: columns: pattern ,start(IRange), end(IRange)
sequence <- dna_samples[3]
pat_1 <- Biostrings::matchPattern(pattern = "TTAGGG", subject = unlist(sequence), max.mismatch = 0)
pat_2 <- Biostrings::matchPattern(pattern = "CCAGGG", subject = unlist(sequence), max.mismatch = 0)

patterns_df <- data_frame(patterns = character(), start_idx = integer(), end_idx = integer())
patterns_df <- add_row(patterns_df, patterns = as.data.frame(pat_1)$x, start_idx = start(pat_1), end_idx = end(pat_1))

Biostrings::matchPattern(pattern = "TTAGGG", subject = unlist(dna_samples[3]), max.mismatch = 0)


my_list <- get_densityIRanges_with_csv(sequence = dna_samples[3], patterns = dna_rc_patterns)

counts <- my_list$`Patterns data frame` %>% group_by(patterns) %>% count()
















































