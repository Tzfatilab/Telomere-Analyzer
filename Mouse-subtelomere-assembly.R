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

# run nanotel on it ....

# create short 40k edges



# trimm telomeres + tvr regions


# create function which inverte the indices of the edges to the real indices of the ref




# This ref do not have Y chr, check how similar Y to X using the other ref


# Now use the assembly to map mouse samples

# This is the same mouse diff age
# 1. Trial 18 bc01 (6 months) - maybe rebascalle? (r9)
# 2. Trial 36 bc01 (16 months) - maybe rebascalle? (r9)
# 3. Trial 47 bc01 (22 months) - maybe rebascalle? (r10)