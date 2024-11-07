library(Biostrings) # you need to google how to install this package if you don't have it installed already

# Read the fasta file containing the selected sequences
seqs <- readDNAStringSet("sequence_data/selected_seqs.fasta") # replace the path with the path to the link to the fasta file if you are using the online runtime

############# DATA CLEANING #############
############# 1. Load the provided SARS-CoV-2 FASTA file (20 sequences, each 29,903 bp) into the Colab notebook.
############# 2. Keep sequences with coverage over 90% (counting only A, C, T, and G bases). 


# Manually calculate the proportion of A, C, T, G in the sequences
## Try to extract one sequence and print it
seq1 <- seqs[[1]]
print(seq1)
### Try to change to a character vector
seq1_char <- as.character(seq1)
print(seq1_char)
### Split the sequence into a list of characters
seq1_char_list <- strsplit(seq1_char, "")
print(seq1_char_list)
### Count the number of each character
seq1_char_table <- table(seq1_char_list)
print(seq1_char_table)

## Repeat the above steps for all the sequences
seqs_char <- lapply(seqs, as.character)
seqs_char_list <- lapply(seqs_char, strsplit, "")
seqs_char_table <- lapply(seqs_char_list, table)
print(seqs_char_table)

## calculate the proportion of A, C, T, G in each sequence
seqs_char_table_prop <- lapply(seqs_char_table, function(x) x/sum(x))
print(seqs_char_table_prop)

## calculate the summed proportion of A, C, T, G in all the sequences
seqs_char_table_prop_actg <- lapply(seqs_char_table_prop, function(x) sum(x[names(x) %in% c("A", "C", "T", "G")]))
print(seqs_char_table_prop_actg)
print(seqs_char_table_prop_actg[seqs_char_table_prop_actg<0.9])
print(seqs[seqs_char_table_prop_actg<0.9])


# Using a built-in function to calculate the proportion of A, C, T, G in the sequences
seqs_nt_count <- alphabetFrequency(seqs)
print(seqs_nt_count)
seqs_nt_prop <- seqs_nt_count/rowSums(seqs_nt_count)
print(seqs_nt_prop)

## Way 1
seqs_nt_prop_actg_1 <- apply(seqs_nt_prop, 1, function(x) sum(x[names(x) %in% c("A", "C", "T", "G")]))
print(seqs_nt_prop_actg_1)
## Way 2
seqs_nt_prop_df <- as.data.frame(seqs_nt_prop)
seqs_nt_prop_actg_2 <- (seqs_nt_prop_df$A + seqs_nt_prop_df$C + seqs_nt_prop_df$T + seqs_nt_prop_df$G)
print(seqs_nt_prop_actg_2)

## Compare the two ways
print(seqs_nt_prop_actg_1 - seqs_nt_prop_actg_2)
(seqs_nt_prop_actg_1 - seqs_nt_prop_actg_2)<0.00000000001

# keeping sequences with at least 90% of A, C, T, G
seqs_filtered <- seqs[seqs_nt_prop_actg_1 >= 0.9]
print(seqs_filtered) # how many sequences are left after filtering?


############# SEQUENCE ANALYSIS #############
############# 1. Calculate GC content for each of two randomly selected sequences.
############# 2. Extract the spike gene region (positions 21,563 to 25,384) for both sequences.
############# 3. Calculate the codon usage for one of the extracted sequences.

# Calculate GC content for each of two randomly selected sequences
set.seed(2024) # set seed for reproducibility, in this way you will get the same random numbers as me every time
seqs_sample <- sample(seqs_filtered, 2)
seqs_sample_nt_count <- alphabetFrequency(seqs_sample)
print(seqs_sample_nt_count)
seqs_sample_nt_count_df <- as.data.frame(seqs_sample_nt_count)
seqs_sample_nt_count_df$length <- rowSums(seqs_sample_nt_count)
seqs_sample_nt_count_df$GC <- seqs_sample_nt_count_df$G + seqs_sample_nt_count_df$C
seqs_sample_nt_count_df$GC_content <- seqs_sample_nt_count_df$GC/seqs_sample_nt_count_df$length

print(names(seqs_sample))
print(seqs_sample_nt_count_df$GC_content)
print(paste0("GC content of sequence ", names(seqs_sample)[1], ": ", seqs_sample_nt_count_df$GC_content[1]))
print(paste0("GC content of sequence ", names(seqs_sample)[2], ": ", seqs_sample_nt_count_df$GC_content[2]))

# Extract the spike gene region (positions 21,563 to 25,384) for both sequences
## visit https://codon2nucleotide.theo.io for positions and annotations of the SARS-CoV-2 genome
spike_gene <- subseq(seqs_sample, start=21563, end=25384)
print(spike_gene)

# Calculate the codon usage for one of the extracted sequences
set.seed(2024)
spike_gene_selected <- spike_gene[sample(1:2, 1)]
length(spike_gene_selected)
width(spike_gene_selected)

## Manually calculate the codon usage
spike_gene_selected_char <- as.character(spike_gene_selected)
print(spike_gene_selected_char)
spike_gene_selected_char_list <- strsplit(spike_gene_selected_char, "")[[1]]
print(spike_gene_selected_char_list)
spike_gene_selected_codons <- lapply(seq(1, length(spike_gene_selected_char_list), by=3), function(x) paste(spike_gene_selected_char_list[x:(x+2)], collapse=""))
print(spike_gene_selected_codons)
codon_usage_manual <- table(unlist(spike_gene_selected_codons))
print(codon_usage_manual)

## Using a built-in function to calculate the codon usage
codon_usage_auto <- trinucleotideFrequency(spike_gene_selected, step = 3)
print(codon_usage_auto)

# check the reverse complement of the spike gene
spike_gene_selected_revcomp <- reverseComplement(spike_gene_selected)
print(spike_gene_selected_revcomp)
