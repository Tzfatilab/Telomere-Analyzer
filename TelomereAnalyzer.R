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
# Script Description: In this script we search over fastq file sequences for telomeric paterns.
#
#
# Notes:
#
#




################################################ TELOMERS patterns analyzer ################################################ 

########## changes for Jan.2022: try to find the teloret + filter from 55:155 ( avoid the adaptor and teloret)
#Telorette3

# P-5’-TGCTCCGTGCATCTCCAAGGTTCTAACC  # 28 nucs -> use countmatch Pattern with 11 mismatches (60%)




############# WORK FOR 12.01.2022 ##########################3
#' 1. make the programm more robust, not only for telomeric sequences
#' 
#' 
#' 2. make 2 plots: the second with total length of 100k (chec the max(seq_length on csv files)) : done
#' 3. for the plots add line for the end of the total read(vertical line).
#' 4. take the output fasta file and make it rc -> 5-3 order ( I can do it on linux-shell or on R): done using reverseComplement()
#' 5. adjust the patterns and telorrete to rc (5-3 order): done
#' 6. change the title to "Telomeric repeat density": done
#' 7. try to add "\n\n" for the sub-title (we don't want it to be so close to the x-title): does not work
#' 8. filter out the reads which thier subtelomeric regin is les then 50b: done with 100b
#' 9. the read after rc (5'-3'): [sub_telo ,   telo_start,    , telo_end,   telorrete,    , adaptor..s]
#' 9. adjust the patterns of the telomere and telorrete(rc - YYAGGG paterns , now since R->Y): done
#' 
#' #' Work for 20.11
#' 10. make it more efficient with data.table  library(data.table)  and library("parallel")
#' 11. Make the functions generic and place them in a project( the ones that don't have to be for telomeres)
#' 12. make the filter function an argument for searchPatten
#' 13. make parallel: get 1 fastq/a file as an input divide to 12 samples which all have seq which pass the filter and then
#'     run in parallel with foreach (i=1:12) %dopar%, make sure I give each samples serial accordingly:
#'     for exm: samp1 is [1:20] so serial[1:20], samp2[21:40] -> serial[21:40]
#' 14. need to find out how to filter with apply...
#' 15. Need to make the Telomere start,end , widht(length) more accurate (maybe use subseq of 20B and not include the telorette ...)
#' 15. chenge the tellorete: fnd only th unvarible bases ( no specific telorrete3)
#' 16. fix the bug: 
#' Error during wrapup: unimplemented type (29) in 'eval'

#Error: no more error handlers available (recursive errors?); invoking 'abort' restart
#Error during wrapup: INTEGER() can only be applied to a 'integer', not a 'unknown type #29'
#Error: no more error handlers available (recursive errors?); invoking 'abort' restart
#bc01 <- readDNAStringSet(filepath = "/home/lab/Downloads/Telomers/Trial 8/fastq_pass-20220508T085644Z-001/fastq_pass/barcode01/FAT23158_pass_barcode01_2c5a6cb7_0.fastq.gz", format = "fastq")
#run_with_rc_and_filter(samples = bc01, patterns = dna_rc_patterns,
#                       output_dir = "/home/lab/Downloads/Telomers/Trial 8/fastq_pass-20220508T085644Z-001/fastq_pass/barcode01/output",
#                       telorrete_pattern = the_telorete_pattern)


#'
#'
#' # plans for future:
#'  1. use plot_ly for interactive plots
#'  2. change the table creation by creating first the vectors and once all finish create from them the dt
#'  3. make adjustment and s..
#'  4. make searchPattern a parallel with foreach dopar
#'  5. create an efficient filter for rawData : I need to figure it out how to use vapply with my own function which returns logical  
#'     I need to extract a character vector of firts/last 100b for each seq in the DNAStringSet and then apply it to the filter which cheks for pattern density. 
#'  6.
#'  7. seperate to several files and create a project
#'  8. the speling of telorette is "telorette" 
#'  9. need to make the telomere length more accutrate: subtract from it the teloette indices...s

# 1. load fata/q file
# 2. filter reads according to MIN_LENGHT
# 3. find the telorette: GCTCCGTGCATCTCCAAGGTTCTAACC  # 28 nucs -> use countmatch Pattern with 11 mismatches
# 4. filter reads according to pattern density at the first subtelomere
# 5. summarise resullts and plot the telomeres
# 7. need to use the "cutter for barcode" and also cut the telorette from the sequence !!!!!( since it can be identified as subtelomeric region ): did it with telorette - need to add barcode option
# 8. add a sub-telomeric only fasta/fastq file output    
# 9. try maybe to use library(microseq) for fastq file output / Or ShortRead ( see Datacap intro to Bioconductor in R chapter 4)
# 10. try using library(microseq) readFastq and use the Quality for more accurate Or use ShortRead ( see intro to Bioconductor in R) .... or the XStringQuality-class {Biostrings} readDNAStringSet(with.qualities = TRUE, ....)
# 11. need to add option: search for all 6 types of telorette for a higher chance of finding one
# 12. add another fasta/fastq file -> telomere_trimed for the mapping to the genome 
# 13. Try QualityScaledDNAStringSet()
# 14.  # change to -(61+ ( barcode_telorette == 61)) at run_with_rc_and_filer - done !
# 15. There is a bug if the DNAstringSet is of length == 1
# 16. I need to adjust the sub-seq length ( 100 if length >= 40,000 , 20)

library(BiocManager)
library('BSgenome')
library('stringr')
library(tidyverse)  
library(IRanges)
library(purrr)
# library(data.table)
library("parallel")

#' @title  
#' @description 
#' @usage 
#' @param 
#' @param 
#' @return
#' @examples




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
analyze_subtelos <- function(dna_seq, patterns , sub_length = 100, MIN_DENSITY = 0.18){ # return list(subtelos, list_density_mp[1])
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
    jpeg_path <- paste(OUTPUT_JPEGS, paste('read', serial_num, '.jpeg',sep=''), sep='/')  
    jpeg(filename=jpeg_path, width=w, height=h)                                                            
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



# has bug, more accurate then the default of getting the density of subsequence , because it can calculate the edges which has the partial of the pattern
filter_first_100 <- function(sequence, patterns, min_density = 0.18, start = 1,end = 100){
  #' Title: filter only sequences which thier subseq has at least minimal density of the given patterns
  #' have density of at least min_density
  #' 
  #' Notice that it is not accurate since it is nly a subsequnce, does not count patterns which are at the edges ( starts before start or ends after end )
  #' 
  #' @param sequence: a dna sequence (DNAString object)
  #' @param pattern: a list of patterns or a string of 1 pattern
  #' @param min_density: threshold for the density of the pattern
  #' @param start: the start index of the subsequnce
  #' @param end: the end "      "     
  #' @return a list of (logical, numeric), the logical indicates if the density of the subsequnce is >= min_density, and the numeric
  end_2 <- end+100 # take a subseuence which is large enough for not missing patterns at the edges
  ranges <- IRanges(start = start, end = end)
  sub_ranges <- get_densityIRanges(sequence = unlist(subseq(sequence,start = 1, end = end2)), patterns = patterns)[2] # return(list(total_density, mp_all))
  densities <- get_sub_density(sub_irange = sub_ranges, ranges = ranges)
  if(densities >= min_density){
    return(list(TRUE, densities))
  }
  else{
    return(list(FALSE, densities))
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
  #'
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
    s_index <- start(telo_irange) 
    # make the strat/end more accurate (usethe IRanges for the patterns)
    iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) 
    if(length(iranges_start) > 0){ start(telo_irange) <- start(irange_telo[iranges_start[1]])} 
    
    e_index <- end(telo_irange) 
    iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
    if(length(iranges_end) >0 ) {end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
    
    
    if(width(telo_irange) < 100 ) {next} # not considered a Telomere
    telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
    
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
    
    LargeDNAStringSet <- append(LargeDNAStringSet, values = sample_telomeres[i])
    current_serial <- current_serial + 1
    
    
  } # end of for loop
  
  # need to save the df in a file
  write.csv(x=df, file=OUTPUT_TELO_CSV)
  writeXStringSet(LargeDNAStringSet, OUTPUT_TELO_FASTA)
  message("Done!") #  now what's left is to extract the sequences from fasta to fasta using the read_names list file ( 3 files )
} # end of the function searchPatterns






searchPatterns_withTelorette <- function(sample_telomeres , pattern_list, max_length = 1e5, csv_name = "summary",output_dir, serial_start = 1, min_density, telorete_pattern, title = "Telomeric repeat density"){
  #'@title Search given Patterns over a DNA sequences.
  #'
  #'
  #'@param sample_telomeres: the DNAString Set of the reads.
  #'@param pattern_list: a list of patterns or a string of 1 pattern.
  #'@param csv_name: The name of the csv file
  #'@param output_dir: 
  #'@param min_density: The minimal density of the patterns in a sequence to be consider relevant.
  #'@param telorete_pattern: The patern of the telorete
  #'@param title: The title for the density plots.
  #'
  OUTPUT_TELO_CSV <- paste(output_dir, paste(csv_name, 'csv', sep='.'), sep='/')
  OUTPUT_TELO_FASTA <- paste(output_dir, paste("reads", 'fasta', sep='.'), sep='/')
  OUTPUT_JPEGS <- paste(output_dir, 'single_read_plots', sep='/')
  dir.create(OUTPUT_JPEGS)
  OUTPUT_JPEGS.1 <- paste(output_dir, 'single_read_plots_adj', sep='/')
  dir.create(OUTPUT_JPEGS.1)
  
  #max_length <- max(width(sample_telomeres))
  
  
  #  I HAVE ADDED TLOMERE LENGTH, START @ END
  df<-data.frame(Serial = integer(), sequence_ID = character(), sequence_length = integer(),  telo_density = double(),
                 Telorette3 = logical(), Telorette3Start_index = integer(), Telorette_seq = character(), Telomere_start = integer(), Telomere_end = integer(), Telomere_length = integer())
  
  
  # add telo density : get_sub_density <- function(sub_irange, ranges){
  LargeDNAStringSet <- DNAStringSet() # For the fasta output of the reads which pass the filter
  current_serial <- serial_start
  
  for( i in 1:length(sample_telomeres) ){
    current_fastq_name <- names(sample_telomeres[i])
    current_seq <- unlist(sample_telomeres[i])
    # we skip the adaptor and telorete so we start with base 57
    
    
    # # returns a a list of (a data frame, list(numeric: total density,iranges)) 
    analyze_list <- analyze_subtelos(dna_seq = current_seq , patterns =  pattern_list, MIN_DENSITY = min_density)
    
    #list_CurrentDens_MP <- get_densityIRanges(current_seq, patterns = PATTERNS_LIST) # to remove: last_100_density = double(), total_density = double(),
    current_telorete <- matchPattern(pattern = telorete_pattern, 
                                     subject = subseq(current_seq, start =length(current_seq) -86, end = length(current_seq)), max.mismatch = 14, with.indels = TRUE)
    
    
    telo_irange <- find_telo_position(seq_length = length(current_seq), subtelos = analyze_list[[1]], min_in_a_row = 3, min_density_score = 2 )
    
    irange_telo <- analyze_list[[2]][[2]]
    s_index <- start(telo_irange) 
    # make the strat/end more accurate (usethe IRanges for the patterns)
    iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100))  # change to 20 if subseq == 20
    if(length(iranges_start) > 0){ start(telo_irange) <- start(irange_telo[iranges_start[1]])} 
    
    e_index <- end(telo_irange) 
    iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
    if(length(iranges_end) >0 ) {end(telo_irange) <- end(irange_telo[iranges_end[1]])} 
    
    
    
    if(width(telo_irange) < 100 ) {next} # not considered a Telomere
    
    if(length(current_telorete)){ # The telorrete was found
      telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
      " check if the telomre starts after the tellorete.. for meanwhile don't do it....
        # check if the telomre starts after the tellorete
        if( start(telo_irange) <= end(current_telorete[[1]]) ){ # check if the telorrete is before the start of the Telomere
          new_start <- end(current_telorete[[1]]) + 1
          if(end(telo_irange) <= new_start) { # The tellorete was found after the telomere
            
            ##################### CHECH IN THE FUTURE ###################
            if( width(telo_irange) < 200 ) {next}# won't be considered Telomere, COULD BE PROBLEMATIC IF WE HAVE TO SUBSEQUENCES WHICH CAN BE TELOMERES
            # else we leave it to be before the telorette
          }
          else{ # update the start of the Telo after the telorrete 
            telo_irange <- IRanges(start = new_start, end = end(telo_irange))
          }
        }
        "              # to remove: last_100_density = double(), total_density = double(),
      df <- df %>% add_row(Serial = current_serial, sequence_ID = current_fastq_name, sequence_length = length(current_seq), telo_density = telo_density,
                           Telorette3 = TRUE, Telorette3Start_index = start(current_telorete[1]) + length(current_seq) -87,
                           Telorette_seq = toString(unlist((current_telorete[1]))), Telomere_start = start(telo_irange), Telomere_end = end(telo_irange), Telomere_length = width(telo_irange))  
    }
    else{ # no telorrete
      df <- df %>% add_row(Serial = current_serial, sequence_ID = current_fastq_name, sequence_length = length(current_seq),telo_density = telo_density,Telorette3 = FALSE,
                           Telorette3Start_index = -1 , Telorette_seq = "",Telomere_start = start(telo_irange), Telomere_end = end(telo_irange), Telomere_length = width(telo_irange))  
    }
    # now save a plot 
    # plot_single_telo <- function(seq_length, subs, serial_num ,save.it=T, legend_string = c("TTAGGG density"),  w=500, h=200, OUTPUT_JPEGS)
    # need to add the start, end positions
    
    # plot_single_telo(seq_length, subs, serial_num, telo_irange,save.it=T, main_title = "", w=500, h=200, OUTPUT_JPEGS){ # add OUTPUT_JPEGS as arg
    # THE PLOT WITH AN ADJUST  X-AXIS LENGTH ACCORDING TO THE READ LENGTH
    
    plot_single_telo(x_length = length(current_seq),seq_length = length(current_seq), subs =  analyze_list[[1]], serial_num = current_serial,
                     seq_start = start(telo_irange),seq_end = end(telo_irange), save.it=T, main_title = title,  w=750, h=300, OUTPUT_JPEGS= OUTPUT_JPEGS.1)
    
    # THE PLOT WITH AN CONSTATNT X-AXIS LENGTH
    plot_single_telo(x_length =max_length,seq_length = length(current_seq), subs =  analyze_list[[1]], serial_num = current_serial ,
                     seq_start = start(telo_irange),seq_end = end(telo_irange), save.it=T, main_title = title,  w=750, h=300, OUTPUT_JPEGS= OUTPUT_JPEGS)
    
    LargeDNAStringSet <- append(LargeDNAStringSet, values = sample_telomeres[i])
    current_serial <- current_serial + 1
  }
  
  # need to save the df in a file
  write.csv(x=df, file=OUTPUT_TELO_CSV)
  writeXStringSet(LargeDNAStringSet, OUTPUT_TELO_FASTA)
  message("Done!") #  now what's left is to extract the sequences from fasta to fasta using the read_names list file ( 3 files )
} # end of the function searchPatterns  











############################ code for 30.01 #########################################################

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



run_with_rc_and_filter <- function(samples,  patterns, output_dir, telorrete_pattern){
  samps_1000 <- samples[width(samples) >= 1e3]
  samps_1000 <- reverseComplement(samps_1000)
  copies_of_r <- 10
  
  cl <- makeCluster(copies_of_r)
  samp_100 <- parLapply(cl, samps_1000, subseq, end = -67, width = 100)  # change to -(61+ just incase there are indels ( barcode_telorette == 61))
  stopCluster(cl)
  
  cl <- makeCluster(copies_of_r)
  logical_100 <- parSapply(cl, samp_100,  filter_density,patterns = patterns , min_density = 0.175)
  stopCluster(cl)
  names(logical_100) <- NULL
  
  samps_filtered <- samps_1000[logical_100]
  
  searchPatterns_withTelorette(samps_filtered, pattern_list = patterns, output_dir = output_dir , min_density = 0.18,
                               telorete_pattern = telorrete_pattern )
  
}




PATTERNS_LIST <- list("CCCTAA", "CCCTAG", "CCCTGA", "CCCTGG")  # CCCTRR
PATTERNS_LIST <- append(PATTERNS_LIST, list("CCTAAC", "CCTAGC", "CCTGAC", "CCTGGC")) # add CCTRRC
PATTERNS_LIST <- append(PATTERNS_LIST, list("CTAACC", "CTAGCC", "CTGACC", "CTGGCC")) # add CTRRCC
PATTERNS_LIST <- append(PATTERNS_LIST, list("TAACCC", "TAGCCC", "TGACCC", "TGGCCC")) # add TRRCCC
PATTERNS_LIST <- append(PATTERNS_LIST, list("AACCCT", "AGCCCT", "GACCCT", "GGCCCT")) # add RRCCCT
PATTERNS_LIST <- append(PATTERNS_LIST,list("ACCCTA", "GCCCTA", "ACCCTG", "GCCCTG"))  # add RCCCTR

patterns_dna <- lapply(PATTERNS_LIST, DNAString)
dna_rc_patterns <- lapply(patterns_dna,reverseComplement)
dna_rc_patterns <- lapply(dna_rc_patterns, toString)

the_telorete_pattern = "TGCTCCGTGCATCTCCAAGGTTCTAACC"
the_telorete_pattern <- toString(reverseComplement(DNAString(the_telorete_pattern)))



################## test new plots 26/07/2022 ################

## Work for today: 1. make it a plot fron subseq of 20, 2. take the Iranges and use it to find the first index of the patterns within the telomere start
# , 3.exclude the telorette+barcode from the plot/seq_length ,( 1+2+3 will make the telomere length more accurate: from 100B error to less then 5B error in length)
# 4.Make an option to find all the 6 different telorettes , 5. change the subs data frame: add indices of start/end of patterns within each sub-seq


filepath1 = "/home/lab/Downloads/Telomers/Trial10 2022_06_15_1335_MN34594_FAS36701_5116e444/output_bc03/reads__bc03_pass.fasta"
filepath2 = "/home/lab/Downloads/Telomers/Trial10 2022_06_15_1335_MN34594_FAS36701_5116e444/pass-20220707T125428Z-001/pass/barcode03/output/reads_bc03_un_mixed.fasta"
x_covered <- readDNAStringSet(filepath = c(filepath1, filepath2))
name_c <- "49f23ce2-8000-4c65-8790-91ecc042855f runid=7603663a2adec7cf9ecee66b7cd4a81aa9032d9a sampleid=no_sample read=69258 ch=104 start_time=2022-06-15T21:19:05Z model_version_id=2021-05-17_dna_r9.4.1_minion_768_2f1c8637 barcode=barcode03"
name_l = "c748e9a4-5166-4a97-9be7-3f9774df17c0 runid=7603663a2adec7cf9ecee66b7cd4a81aa9032d9a read=39099 ch=340 start_time=2022-06-16T00:42:45.655343+03:00 flow_cell_id=FAS36701 protocol_group_id=10KKbtMEFsMMpd110 sample_id=no_sample barcode=barcode03 barcode_alias=barcode03 parent_read_id=c748e9a4-5166-4a97-9be7-3f9774df17c0 basecall_model_version_id=2021-05-17_dna_r9.4.1_minion_384_d37a2ab9"
x_covered <- x_covered[which(names(x_covered) %in% c(name_c, name_l))] # this is the plot 55 in the mixed_bc03 from Trial10 which Dudy toled me 
# See this telomere. 130K! but the label is covering the graph !
x_covered <- reverseComplement(x_covered) # needed since it will be rc again !

setwd("/home/lab/Downloads/Telomers/plot_test")

dir.create("output_plotTest")
run_with_rc_and_filter(samples = x_covered, patterns = dna_rc_patterns,
                       output_dir = "output_plotTest",
                       telorrete_pattern = the_telorete_pattern )


# now try diff plots

plot_df <- read_csv(file = "/home/lab/Downloads/Telomers/plot_test/output_plotTest/output_sub20_AdjustTelomerePosition/single_read_plots_adj/read1.csv")
ggplot(plot_df , aes(x = start_index, y = density)) + geom_area(colour = "red", fill = "red")


plot_df100 <- read_csv(file = "/home/lab/Downloads/Telomers/plot_test/output_plotTest/old_output/single_read_plots/read1.csv")
ggplot(plot_df100 , aes(x = start_index, y = density )) + geom_area(colour = "red", fill = "red") +  
  scale_x_continuous(name="Position", limits=c(1, 123501), breaks = seq(10000, 123501, by = 10000)) +
  theme_light() + legend() # how to inser a legend ???????




# try to adjust the telomere start/end indices using the IRanges indices + telorete_startIndex
df<-data.frame(Serial = integer(), sequence_ID = character(), sequence_length = integer(),  telo_density = double(),
               Telorette3 = logical(), Telorette3Start_index = integer(), Telorette_seq = character(), Telomere_start = integer(), Telomere_end = integer(), Telomere_length = integer())


# add telo density : get_sub_density <- function(sub_irange, ranges){
LargeDNAStringSet <- DNAStringSet() # For the fasta output of the reads which pass the filter
current_serial <- serial_start


  current_fastq_name <- names(x_covered[1])
  current_seq <- unlist(x_covered[1])
  current_seq <- Biostrings::reverseComplement(current_seq)
  # we skip the adaptor and telorete so we start with base 57
  
  
  # # returns a a list of (a data frame, list(numeric: total density,iranges)) 
  analyze_list <- analyze_subtelos(dna_seq = current_seq , patterns =  dna_rc_patterns, MIN_DENSITY = 0.18)
  
  #list_CurrentDens_MP <- get_densityIRanges(current_seq, patterns = PATTERNS_LIST) # to remove: last_100_density = double(), total_density = double(),
  current_telorete <- matchPattern(pattern = the_telorete_pattern, 
                                   subject = subseq(current_seq, start =length(current_seq) -86, end = length(current_seq)), max.mismatch = 11, with.indels = T)
  
  
  # for the telorete : 3 stage: first whith indels, and max.mismatch == 5, ... last with indels , max.mismatch == 14 : use the width (max width ) to select the irange
  telo_irange <- find_telo_position(seq_length = length(current_seq), subtelos = analyze_list[[1]], min_in_a_row = 3, min_density_score = 2 )
  
  
  
  
  irange_telo <- analyze_list[[2]][[2]]
  s_index <- start(telo_irange) # 100401 
  iranges_start <- which(start(irange_telo) %in% s_index:(s_index + 100)) # 100 or 20 ( The length of the subseq)
  if(length(iranges_start) > 0){ s_index <- start(irange_telo[iranges_start[1]])} # fixed to 100439
  
  
  
  e_index <- end(telo_irange) # 124142
  iranges_end <- which(end(irange_telo) %in% (e_index - 100):e_index)
  if(length(iranges_end) >0 ) {e_index <- end(irange_telo[iranges_end[1]])} # fixed to 124060
  # There is improvment of 118B at Telomere length for this example.
  
  # chech the se itself : old_start_end = [ 100401, 124142], new_start_end = [100439, 124060 ]
  Biostrings::subseq(x = current_seq, start = 100401, width = 100)
  
  # IRanges: find IRange which is within other IRange
  
  
  
  
  
  which(start(irange_telo) <= start(telo_irange) + 20 && start(irange_telo) >= start(telo_irange) )
  
  
  if(width(telo_irange) < 100 ) {next} # not considered a Telomere
  
  if(length(current_telorete)){ # The telorrete was found
    telo_density <- get_sub_density(telo_irange, analyze_list[[2]][[2]])
    " check if the telomre starts after the tellorete.. for meanwhile don't do it....
        # check if the telomre starts after the tellorete
        if( start(telo_irange) <= end(current_telorete[[1]]) ){ # check if the telorrete is before the start of the Telomere
          new_start <- end(current_telorete[[1]]) + 1
          if(end(telo_irange) <= new_start) { # The tellorete was found after the telomere
            
            ##################### CHECH IN THE FUTURE ###################
            if( width(telo_irange) < 200 ) {next}# won't be considered Telomere, COULD BE PROBLEMATIC IF WE HAVE TO SUBSEQUENCES WHICH CAN BE TELOMERES
            # else we leave it to be before the telorette
          }
          else{ # update the start of the Telo after the telorrete 
            telo_irange <- IRanges(start = new_start, end = end(telo_irange))
          }
        }
        "              # to remove: last_100_density = double(), total_density = double(),
    df <- df %>% add_row(Serial = current_serial, sequence_ID = current_fastq_name, sequence_length = length(current_seq), telo_density = telo_density,
                         Telorette3 = TRUE, Telorette3Start_index = start(current_telorete[1]) + length(current_seq) -87,
                         Telorette_seq = toString(unlist((current_telorete[1]))), Telomere_start = start(telo_irange), Telomere_end = end(telo_irange), Telomere_length = width(telo_irange))  
  }



################## End of test new plots 26/07/2022 ################









#################### Trial 12  20.07.2022 #################
setwd("/home/lab/Downloads/Telomers/Trial 12/prelog/barcode01")
bc01 <-  readDNAStringSet(dir(), format = "fastq")
dir.create("output_bc01")
run_with_rc_and_filter(samples = bc01, patterns = dna_rc_patterns,
                       output_dir = "output_bc01",
                       telorrete_pattern = the_telorete_pattern)

setwd("/home/lab/Downloads/Telomers/Trial 12/prelog/barcode02")
bc02 <-  readDNAStringSet(dir(), format = "fastq")
dir.create("output_bc02")
run_with_rc_and_filter(samples = bc02, patterns = dna_rc_patterns,
                       output_dir = "output_bc02",
                       telorrete_pattern = the_telorete_pattern)

setwd("/home/lab/Downloads/Telomers/Trial 12/prelog/barcode03")
bc03 <-  readDNAStringSet(dir(), format = "fastq")
dir.create("output_bc03")
run_with_rc_and_filter(samples = bc03, patterns = dna_rc_patterns,
                       output_dir = "output_bc03",
                       telorrete_pattern = the_telorete_pattern)
###########################################################


#################### test 3 12.07.2022 #######################
test <- readDNAStringSet(filepath = "/home/lab/Downloads/Telomers/Trial7B/barcode01/test2/reads.fasta", format = 'fasta')
test <- reverseComplement(test)
setwd("/home/lab/Downloads/Telomers/Trial7B/barcode01/test2")
run_with_rc_and_filter(samples = test, patterns = dna_rc_patterns,
                       output_dir = "output",
                       telorrete_pattern = the_telorete_pattern)



######################## 16.06.2022: 1335_MN34594_FAS36701_5116e444 ###########################################

# bc01
bc1_0 <- readDNAStringSet(filepath = "barcode01/FAS36701_pass_barcode01_7603663a_0.fastq.gz", format = 'fastq') 
bc1_1 <- readDNAStringSet(filepath = "barcode01/FAS36701_pass_barcode01_7603663a_1.fastq.gz", format = 'fastq')   
bc1_2 <- readDNAStringSet(filepath =  "barcode01/FAS36701_pass_barcode01_7603663a_2.fastq.gz", format = 'fastq') 
bc1_3 <- readDNAStringSet(filepath = "barcode01/FAS36701_pass_barcode01_7603663a_3.fastq.gz", format = 'fastq') 
bc1_all <- append(x = append(bc1_0, bc1_1), append(bc1_2, bc1_3))

run_with_rc_and_filter(samples = bc1_all, patterns = dna_rc_patterns,
                       output_dir = "output_bc01",
                       telorrete_pattern = the_telorete_pattern)


#bc 02
bc2_0 <- readDNAStringSet(filepath = "barcode02/FAS36701_pass_barcode02_7603663a_0.fastq.gz", format = 'fastq') 
bc2_1 <- readDNAStringSet(filepath = "barcode02/FAS36701_pass_barcode02_7603663a_1.fastq.gz", format = 'fastq')
bc2_2 <- readDNAStringSet(filepath = "barcode02/FAS36701_pass_barcode02_7603663a_2.fastq.gz" , format = 'fastq')
bc2_3 <- readDNAStringSet(filepath ="barcode02/FAS36701_pass_barcode02_7603663a_3.fastq.gz" , format = 'fastq') 
bc2_4 <- readDNAStringSet(filepath ="barcode02/FAS36701_pass_barcode02_7603663a_4.fastq.gz" , format = 'fastq')
bc2_5 <- readDNAStringSet(filepath = "barcode02/FAS36701_pass_barcode02_7603663a_5.fastq.gz", format = 'fastq')
bc2_all <- append(x = append(bc2_0, bc2_1), append(append(bc2_2, bc2_3), append(bc2_4, bc2_5) ) )

run_with_rc_and_filter(samples = bc2_all, patterns = dna_rc_patterns,
                       output_dir = "output_bc02",
                       telorrete_pattern = the_telorete_pattern)




# bc 03
bc3_0 <- readDNAStringSet(filepath = "barcode03/FAS36701_pass_barcode03_7603663a_0.fastq.gz", format = 'fastq') 
bc3_1 <- readDNAStringSet(filepath = "barcode03/FAS36701_pass_barcode03_7603663a_1.fastq.gz", format = 'fastq')

bc3_all <- append(bc3_0, bc3_1)

run_with_rc_and_filter(samples = bc3_all, patterns = dna_rc_patterns,
                       output_dir = "output_bc03",
                       telorrete_pattern = the_telorete_pattern)

run_with_rc_and_filter(samples = bc1_all, patterns = dna_rc_patterns,
                       output_dir = "output_bc01",
                       telorrete_pattern = the_telorete_pattern)



# now run fail fastq
bc1_fail <-  readDNAStringSet(filepath = "barcode01_FAILED/FAS36701_fail_barcode01_7603663a_0.fastq.gz" , format = 'fastq') 
run_with_rc_and_filter(samples = bc1_fail, patterns = dna_rc_patterns,
                       output_dir = "failed_1" ,
                       telorrete_pattern = the_telorete_pattern)

bc2_fail <-  readDNAStringSet(filepath = "barcode02FAIL/FAS36701_fail_barcode02_7603663a_0.fastq.gz"   , format = 'fastq') 
run_with_rc_and_filter(samples = bc2_fail, patterns = dna_rc_patterns,
                       output_dir = "failed_2" ,
                       telorrete_pattern = the_telorete_pattern)

bc3_fail <-  readDNAStringSet(filepath = "barcode03FAIL/FAS36701_fail_barcode03_7603663a_0.fastq.gz" , format = 'fastq') 
run_with_rc_and_filter(samples = bc3_fail, patterns = dna_rc_patterns,
                       output_dir = "failed_3" ,
                       telorrete_pattern = the_telorete_pattern)



######################## END 16.06.2022: 1335_MN34594_FAS36701_5116e444 ###########################################




