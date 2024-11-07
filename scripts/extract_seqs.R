library(Biostrings)

# Read the fasta file containing all the sequences
fasta_all <- readDNAStringSet("/Volumes/Ext_28T/backup_2024-02-25/work/2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta")

# Select 20 sequences randomly
set.seed(20241108)
id_selected <- sample(1:length(fasta_all), 20)
fasta_selected <- fasta_all[id_selected]

# write the selected sequences to a new fasta file
writeXStringSet(fasta_selected, file = "sequence_data/selected_seqs.fasta")
