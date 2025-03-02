# todo: create an Rmarkdown report on it...
library(Biostrings)
library(tidyverse)


index_adjust <- function(index, edge_name) {
  
}



# links for the internet:
# 1. https://github.com/yulab-ql/mhaESC_genome/releases  - the assembly
# 2. https://www.science.org/doi/10.1126/science.adq8191 -The paper



# load the full assembly
mm_file <-"~/Documents/Mouse-assembly/mouse.241018.v1.1.0.combined.fasta.gz"
ref_mouse <- readDNAStringSet(filepath = mm_file)


test1_head <- Biostrings::subseq(x =  ref_mouse[1], start = 1, end = 20000)


test1_tail <- Biostrings::subseq(x =  ref_mouse[1], end = width(ref_mouse[1]), width = 20000)


# create an ~50k assemly 
mm_50kref <- Biostrings::DNAStringSet()
for(i in 1:20) {
  curr_head <- subseq(ref_mouse[i], start = 1, end = 50000)
  curr_tail <- subseq(ref_mouse[i], end = width(ref_mouse[i]), width = 50000)
  # add names
  names(curr_head) <- str_c(names(curr_head), "_Head")
  names(curr_tail) <- str_c(names(curr_tail), "_Tail")
  
  mm_50kref <- append(mm_50kref, curr_head)
  mm_50kref <- append(mm_50kref, curr_tail)
  
}

Biostrings::writeXStringSet(x =  mm_50kref, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges.fasta")
# split to head and tail
mm_50k_ref_Head <- mm_50kref[str_ends(names(mm_50kref), pattern = "Head")]
mm_50k_ref_Tail <- mm_50kref[str_ends(names(mm_50kref), pattern = "Tail")]
Biostrings::writeXStringSet(x =  mm_50k_ref_Head, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Head.fasta")
Biostrings::writeXStringSet(x = mm_50k_ref_Tail, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Tail.fasta")


# run nanotel on it ....
Rscript --vanilla NanoTel.R -i /home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Head.fasta --save_path /home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Head --format fasta --patterns CCCTAA --min_density 0.5
;Rscript --vanilla NanoTel.R -i /home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Tail.fasta --save_path /home/tzfati/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Tail --format fasta --patterns CCCTAA --check_right_edge --min_density 0.5


# create short 40k edges

# Heads
head_1 <- mm_50k_ref_Head[1]  # chr1 left edge - as the summary file
df_head <- read_csv(file = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-50k_edges-Head/summary.csv")
df_head <- arrange(df_head, sequence_ID)

i <- 20 # I have a bug: the Telomere end is +1 than what it should be
df <- df_head[i, c(2, 6, 10)]
curr_head <- mm_50k_ref_Head[i]
df

# head1 without the telomeric part
head_1_trimmed <-  subseq(mm_50k_ref_Head[1],  start = 6138, width = 40000) # + 6137 for indices of mapping  (CCCTAA)n CCC
head_2_trimmed <- subseq(mm_50k_ref_Head[2], start = 7362 , width = 40000) # + 7361 for indices of mapping : (CCCTAA)nCCC -> should include CCC
head_3_trimmed <- subseq(mm_50k_ref_Head[3], start = 8055, width = 40000) # + 8054 for indices of mapping  : (CCCTAA)n
head_4_trimmed <- subseq(mm_50k_ref_Head[4], start = 5717 , width = 40000) # + 5716 for indices of mapping (CCCTAA)n CCC
head_5_trimmed <- subseq(mm_50k_ref_Head[5], start = 7362, width = 40000) # + 7361 for indices of mapping :(CCCTAA)n CCC
head_6_trimmed <- subseq(mm_50k_ref_Head[6], start = 7050 , width = 40000) # + 7049 for indices of mapping  (CCCTAA)n CCC
head_7_trimmed <- subseq(mm_50k_ref_Head[7], start = 6521 , width = 40000) # + 6520 for indices of mapping (CCCTAA)n CCC
head_8_trimmed <- subseq(mm_50k_ref_Head[8], start = 7889, width = 40000) # + 7888 for indices of mapping  (CCCTAA)n CCC
head_9_trimmed <- subseq(mm_50k_ref_Head[9], start = 8054, width = 40000) # + 8053 for indices of mapping ( CCCTAA)n CCCTA
head_10_trimmed <- subseq(mm_50k_ref_Head[10], start = 8053, width = 40000) # + 8052 for indices of mapping ( CCCTAA)n CCCT
head_11_trimmed <- subseq(mm_50k_ref_Head[11], start = 3464 , width = 40000) # + 3463 for indices of mapping (CCCTAA)n CCC
head_12_trimmed <- subseq(mm_50k_ref_Head[12], start = 8342, width = 40000) # + 8341 for indices of mapping (CCCTAA)n CCCT
head_13_trimmed <- subseq(mm_50k_ref_Head[13], start = 4684 , width = 40000) # + 4683 for indices of mapping (CCCTAA)n CCC
head_14_trimmed <- subseq(mm_50k_ref_Head[14], start = 7362, width = 40000) # + 7361 for indices of mapping :(CCCTAA)n CCC
head_15_trimmed <- subseq(mm_50k_ref_Head[15], start = 7362, width = 40000) # 7361 ? for indices of mapping :(CCCTAA)n CCC
head_16_trimmed <- subseq(mm_50k_ref_Head[16], start = 8057, width = 40000) # + 8056 for indices of mapping (CCCTAA)n CC
head_17_trimmed <- subseq(mm_50k_ref_Head[17], start = 7362, width = 40000) # + 7361 for indices of mapping (CCCTAA)n CCC
head_18_trimmed <- subseq(mm_50k_ref_Head[18], start = 7362, width = 40000) # + 7361 for indices of mapping (CCCTAA)n CCC
head_19_trimmed <- subseq(mm_50k_ref_Head[19], start = 7362, width = 40000) # + 7361 for indices of mapping (CCCTAA)n CCC
head_X_trimmed <- subseq(mm_50k_ref_Head[20], start = 7050, width = 40000) # + 7049 for indices of mapping (CCCTAA)n CCC


ref_40k_head_telo_trimmed <- DNAStringSet()

# rearrange accoridng to start index
head_11_trimmed <- subseq(mm_50k_ref_Head[11], start = 3464 , width = 40000) # + 3463 for indices of mapping (CCCTAA)n CCC
head_13_trimmed <- subseq(mm_50k_ref_Head[13], start = 4684 , width = 40000) # + 4683 for indices of mapping (CCCTAA)n CCC
head_4_trimmed <-  subseq(mm_50k_ref_Head[4],  start = 5717 , width = 40000) # + 5716 for indices of mapping (CCCTAA)n CCC
head_1_trimmed <-  subseq(mm_50k_ref_Head[1],  start = 6138, width = 40000) # + 6137 for indices of mapping  (CCCTAA)n CCC
head_7_trimmed <-  subseq(mm_50k_ref_Head[7],  start = 6521 , width = 40000) # + 6520 for indices of mapping (CCCTAA)n CCC

# check if the subtelomeres are equal : equal take one and call it 6_X
head_6_trimmed <- subseq(mm_50k_ref_Head[6],  start = 7050 , width = 40000) # + 7049 for indices of mapping (CCCTAA)n CCC
head_X_trimmed <- subseq(mm_50k_ref_Head[20], start = 7050, width = 40000) # + 7049 for indices of mapping  (CCCTAA)n CCC


# check if the subtelomeres are equal : all equal: take one and call it 2_5_14_15_17_18_19
head_2_trimmed <- subseq(mm_50k_ref_Head[2],   start = 7362 , width = 40000) # + 7361 for indices of mapping : (CCCTAA)nCCC -> should include CCC
head_5_trimmed <- subseq(mm_50k_ref_Head[5],   start = 7362, width = 40000) # + 7361 for indices of mapping  :(CCCTAA)n CCC
head_14_trimmed <- subseq(mm_50k_ref_Head[14], start = 7362, width = 40000) # + 7361 for indices of mapping  :(CCCTAA)n CCC
head_15_trimmed <- subseq(mm_50k_ref_Head[15], start = 7362, width = 40000) # 7361 ? for indices of mapping  :(CCCTAA)n CCC
head_17_trimmed <- subseq(mm_50k_ref_Head[17], start = 7362, width = 40000) # + 7361 for indices of mapping  :(CCCTAA)n CCC
head_18_trimmed <- subseq(mm_50k_ref_Head[18], start = 7362, width = 40000) # + 7361 for indices of mapping  :(CCCTAA)n CCC
head_19_trimmed <- subseq(mm_50k_ref_Head[19], start = 7362, width = 40000) # + 7361 for indices of mapping  :(CCCTAA)n CCC


head_8_trimmed <- subseq(mm_50k_ref_Head[8],   start = 7889, width = 40000) # + 7888 for indices of mapping   (CCCTAA)n CCC
head_10_trimmed <- subseq(mm_50k_ref_Head[10], start = 8053, width = 40000) # + 8052 for indices of mapping ( CCCTAA)n CCCT
head_9_trimmed <-  subseq(mm_50k_ref_Head[9],  start = 8054, width = 40000) # + 8053 for indices of mapping ( CCCTAA)n CCCTA
head_3_trimmed <-  subseq(mm_50k_ref_Head[3],  start = 8055, width = 40000) # + 8054 for indices of mapping ( CCCTAA)n
head_16_trimmed <- subseq(mm_50k_ref_Head[16], start = 8057, width = 40000) # + 8056 for indices of mapping (CCCTAA)n CC
head_12_trimmed <- subseq(mm_50k_ref_Head[12], start = 8342, width = 40000) # + 8341 for indices of mapping (CCCTAA)n CCCT


# todo: create a ref  seq for Head
Biostrings::writeXStringSet(x = ref_40k_head_telo_trimmed, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-40k-Head-teloTrimmed.fasta")



# Problem: alot of heads have the same row fro msummary file
df_head2 <- select(df_head, -c(1, 2))
unique(df_head2) # we have only 11 unique rows: check the seq : we have 13 uniqke 40k subtelomers 



# tail: I need to re-tkae with longer because some have tel length > 10,000
# create an ~50k assemly 
mm_60kref_tail <- Biostrings::DNAStringSet()
for(i in 1:20) {
  curr_tail <- subseq(ref_mouse[i], end = width(ref_mouse[i]), width = 60000)
  # add names
  names(curr_tail) <- str_c(names(curr_tail), "_Tail")
  mm_60kref_tail <- append(mm_60kref_tail, curr_tail)
  
}
Biostrings::writeXStringSet(x = mm_60kref_tail, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-60k_edges-Tail.fasta")


# trimm telomeres + tvr regions
df_tail <- read_csv("~/Documents//Mouse-assembly/mouse.241018.v1.1.0.-60k_edges-Tail/summary.csv")
df_tail <- arrange(df_tail, sequence_ID)


i <- 20 # I have a bug: the Telomere end is +1 than what it should be
df <- df_tail[i, c(1,2, 5, 9)]
curr_tail <- mm_60kref_tail[i]
df

# todo: fix partial TTAGGG at the start should be included ... also partial end....
tail_01_trimmed <- subseq(mm_60kref_tail[1],  end = 58948, width = 40000) # + ? for indices of mapping? need to check the length of crh  TAGGG (TTAGGG)n
tail_02_trimmed <- subseq(mm_60kref_tail[2],  end = 53582, width = 40000) # + ? for indices of mapping TTATGG (TTAGGG)n
tail_03_trimmed <- subseq(mm_60kref_tail[3],  end = 52750, width = 40000) # + ? for indices of mapping  GGG (TTAGGG)n
tail_04_trimmed <- subseq(mm_60kref_tail[4],  end = 58785, width = 40000) # + ? for indices of mapping GGG (TTAGGG)n
tail_05_trimmed <- subseq(mm_60kref_tail[5],  end = 58365 , width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_06_trimmed <- subseq(mm_60kref_tail[6],  end = 55593, width = 40000) # + ? for indices of mapping  TAAGGG (TTAGGG)n
tail_07_trimmed <- subseq(mm_60kref_tail[7],  end = 55685 , width = 40000) # + ? for indices of mapping AGGG (TTAGGG)n
tail_08_trimmed <- subseq(mm_60kref_tail[8],  end = 51949, width = 40000) # + ? for indices of mapping TTGGGG (TTAGGG)n
tail_09_trimmed <- subseq(mm_60kref_tail[9],  end = 55922, width = 40000) # + ? for indices of mapping   TTCAGGG (TTAGGG)n
tail_10_trimmed <- subseq(mm_60kref_tail[10], end = 55961 , width = 40000) # + ? for indices of mapping TAGGG (TTAGGG)n
tail_11_trimmed <- subseq(mm_60kref_tail[11], end = 51508 , width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_12_trimmed <- subseq(mm_60kref_tail[12], end = 50037, width = 40000) # + ? for indices of mapping TTCGGG (TTAGGG)n
tail_13_trimmed <- subseq(mm_60kref_tail[13], end = 46932, width = 40000) # + ? for indices of mapping GGG (TTAGGG)n
tail_14_trimmed <- subseq(mm_60kref_tail[14], end = 45508 , width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_15_trimmed <- subseq(mm_60kref_tail[15], end = 47734, width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_16_trimmed <- subseq(mm_60kref_tail[16], end = 52294, width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_17_trimmed <- subseq(mm_60kref_tail[17], end = 57764, width = 40000) # + ? for indices of mapping TCAGGG (TTAGGG)n
tail_18_trimmed <- subseq(mm_60kref_tail[18], end = 53568, width = 40000) # + ? for indices of mapping (TTAGGG)n
tail_19_trimmed <- subseq(mm_60kref_tail[19], end = 52148, width = 40000) # + ? for indices of mapping TTTGGG  (TTAGGG)n
tail_X_trimmed <- subseq(mm_60kref_tail[20],  end = 54590, width = 40000) # + ? for indices of mapping (TTAGGG)n

mm_40kref_tail_telo_trimmed <- Biostrings::DNAStringSet()
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_01_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_02_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_03_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_04_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_05_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_06_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_07_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_08_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_09_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_10_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_11_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_12_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_13_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_14_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_15_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_16_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_17_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_18_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_19_trimmed)
mm_40kref_tail_telo_trimmed <- append(mm_40kref_tail_telo_trimmed, tail_X_trimmed)

Biostrings::writeXStringSet(x = mm_40kref_tail_telo_trimmed, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-40k_edges-Tail-telo-trimmed.fasta")


# now merge all 40k : head + tail
mm_40k_telo_trimmed <- append(ref_40k_head_telo_trimmed, mm_40kref_tail_telo_trimmed)
Biostrings::writeXStringSet(x = mm_40k_telo_trimmed, filepath = "~/Documents/Mouse-assembly/mouse.241018.v1.1.0.-40k_edges-telo-trimmed.fasta")

# create function which inverte the indices of the edges to the real indices of the ref




# This ref do not have Y chr, check how similar Y to X using the other ref


# Now use the assembly to map mouse samples

# This is the same mouse diff age
# 1. Trial 18 bc01 (6 months) - maybe rebascalle? (r9)
# 2. Trial 36 bc01 (16 months) - maybe rebascalle? (r9)
# 3. Trial 47 bc01 (22 months) - maybe rebascalle? (r10)