library(ape)
library(pegas)
library(RColorBrewer)
library(msa)
library(rstudioapi)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################ KELCH13 ##########
seqs <- readDNAStringSet("input/Pfk13_consensus_aligned.fa")
metadata<- read.csv("input/District_and_Sample_ID_Metadata_oct11.csv")



# Filter
dna_matrix <- as.matrix(seqs)
sample_threshold <- 1  # Allow up to 10% missing data in a sample
position_threshold <- 1  # Allow up to 10% missing data at a position

missing_per_sample <- rowMeans(dna_matrix == "N" | dna_matrix == "-")
filtered_matrix <- dna_matrix[missing_per_sample <= sample_threshold, ]
missing_per_position <- colMeans(filtered_matrix == "N" | filtered_matrix == "-")
filtered_matrix <- filtered_matrix[, missing_per_position <= position_threshold]

filtered_seqs <- DNAStringSet(apply(filtered_matrix, 1, paste0, collapse = ""))
writeXStringSet(filtered_seqs, "K13_filtered_aligned_sequences.fasta")

# Add district to sample name
current_names <- names(filtered_seqs)
metadata_ordered <- metadata[match(current_names, metadata$Sample_ID), ]
new_names <- paste0(metadata_ordered$Sample_ID, "_", metadata_ordered$District)
names(filtered_seqs) <- new_names

filtered_seqs_bin <- as.DNAbin(filtered_seqs)
haps <- haplotype(filtered_seqs_bin)

dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))
# Set a threshold for rare haplotypes
hap_threshold <- 2 #Filter out haplotypes that occur fewer than 2 times
# Identify haplotypes to keep (those occurring >= threshold)
haps_to_keep <- names(sz)[sz >= hap_threshold]
haps_to_keep_indices <- unlist(attr(haps, "index")[sz >= hap_threshold])

# Subset the sequences using indices for haplotypes to keep
filtered_seqs_bin <- filtered_seqs_bin[haps_to_keep_indices]

# Recompute the haplotypes after filtering
haps_filtered <- haplotype(filtered_seqs_bin)

# Recompute the distance matrix and network
dist <- dist.dna(haps_filtered, "N")
net <- rmst(dist, B = 100)

nt.labs <- attr(net, "labels")

# Update size summary for filtered haplotypes

sz_filtered <- log2(summary(haps_filtered)) * 100
sz_filtered <- sz_filtered[nt.labs]

# Create a matrix of filtered sequences with updated row names
matrix_filtered <- as.matrix(filtered_seqs_bin)

# Ensure row names match sequence names
rownames(matrix_filtered) <- names(filtered_seqs_bin)

# Recalculate region matrix
regions_filtered <- haploFreq(matrix_filtered, split = "_", what = 3)

regions_filtered <- regions_filtered[nt.labs , ]

# Plot with the filtered data
colors <- colorRampPalette(brewer.pal(8, "Accent"))(15)

png("K13_haplotype_network.png" , units = "in" ,width = 12,height = 7  , res = 600)
plot(net, size = sz_filtered, scale.ratio = 200, cex = 0.1 , pie = regions_filtered, threshold = 0, col.link = "black", bg = colors , show.mutation = 0
)
legend("topright", colnames(regions_filtered), col=colors, pch=20, cex=1.6, pt.cex=3)

dev.off()

##########CRT#############

seqs <- readDNAStringSet("input/Pfcrt_consensus_aligned.fa")
metadata<- read.csv("input/District_and_Sample_ID_Metadata_oct11.csv")

# Filter
dna_matrix <- as.matrix(seqs)
sample_threshold <- 0.1  # Allow up to 10% missing data in a sample
position_threshold <- 0.1  # Allow up to 10% missing data at a position

missing_per_sample <- rowMeans(dna_matrix == "N" | dna_matrix == "-")
filtered_matrix <- dna_matrix[missing_per_sample <= sample_threshold, ]
missing_per_position <- colMeans(filtered_matrix == "N" | filtered_matrix == "-")
filtered_matrix <- filtered_matrix[, missing_per_position <= position_threshold]

filtered_seqs <- DNAStringSet(apply(filtered_matrix, 1, paste0, collapse = ""))
writeXStringSet(filtered_seqs, "Crt_filtered_aligned_coding_sequences.fasta")


#Add district to sample name
current_names <- names(filtered_seqs)
metadata_ordered <- metadata[match(current_names, metadata$Sample_ID), ]
new_names <- paste0(metadata_ordered$Sample_ID, "_", metadata_ordered$District)
names(filtered_seqs) <- new_names

filtered_seqs_bin <- as.DNAbin(filtered_seqs)
haps <- haplotype(filtered_seqs_bin)

dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))

# Set a threshold for rare haplotypes
hap_threshold <- 2 #Filter out haplotypes that occur fewer than 2 times
# Identify haplotypes to keep (those occurring >= threshold)
haps_to_keep <- names(sz)[sz >= hap_threshold]
haps_to_keep_indices <- unlist(attr(haps, "index")[sz >= hap_threshold])

# Subset the sequences using indices for haplotypes to keep
filtered_seqs_bin <- filtered_seqs_bin[haps_to_keep_indices]

# Recompute the haplotypes after filtering
haps_filtered <- haplotype(filtered_seqs_bin)

# Recompute the distance matrix and network
dist <- dist.dna(haps_filtered, "N")
net <- rmst(dist, B = 100)

nt.labs <- attr(net, "labels")

# Update size summary for filtered haplotypes


sz_filtered <- log2(summary(haps_filtered)) * 70

sz_filtered <- sz_filtered[nt.labs]



# Create a matrix of filtered sequences with updated row names
matrix_filtered <- as.matrix(filtered_seqs_bin)

# Ensure row names match sequence names
rownames(matrix_filtered) <- names(filtered_seqs_bin)

# Recalculate region matrix
regions_filtered <- haploFreq(matrix_filtered, split = "_", what = 3)

regions_filtered <- regions_filtered[nt.labs , ]


# Plot with the filtered data
colors <- colorRampPalette(brewer.pal(8, "Accent"))(15)

plot(net, size = sz_filtered, scale.ratio = 50, cex = 0.1 , pie = regions_filtered, threshold = 0, col.link = "black", bg = colors)
legend("topright", colnames(regions_filtered), col=colors, pch=20, cex=1.5, pt.cex=3)
#legend("topleft" , legend = legendlable  , cex =1 , pt.cex =3  )

dev.off()

################ MDR1 ##########
seqs <- readDNAStringSet("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/Pfmdr1_consensus_aligned.fa")
metadata<- read.csv("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/District_and_Sample_ID_Metadata_oct11.csv")

# Filter
dna_matrix <- as.matrix(seqs)
sample_threshold <- 0.1  # Allow up to 10% missing data in a sample
position_threshold <- 0.1  # Allow up to 10% missing data at a position

missing_per_sample <- rowMeans(dna_matrix == "N" | dna_matrix == "-")
filtered_matrix <- dna_matrix[missing_per_sample <= sample_threshold, ]
missing_per_position <- colMeans(filtered_matrix == "N" | filtered_matrix == "-")
filtered_matrix <- filtered_matrix[, missing_per_position <= position_threshold]

filtered_seqs <- DNAStringSet(apply(filtered_matrix, 1, paste0, collapse = ""))
writeXStringSet(filtered_seqs, "Mdr1_filtered_aligned_sequences.fasta")


#Add district to sample name
current_names <- names(filtered_seqs)
metadata_ordered <- metadata[match(current_names, metadata$Sample_ID), ]
new_names <- paste0(metadata_ordered$Sample_ID, "_", metadata_ordered$District)
names(filtered_seqs) <- new_names

filtered_seqs_bin <- as.DNAbin(filtered_seqs)
haps <- haplotype(filtered_seqs_bin)

dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))

# Set a threshold for rare haplotypes
hap_threshold <- 2 #Filter out haplotypes that occur fewer than 2 times
# Identify haplotypes to keep (those occurring >= threshold)
haps_to_keep <- names(sz)[sz >= hap_threshold]
haps_to_keep_indices <- unlist(attr(haps, "index")[sz >= hap_threshold])

# Subset the sequences using indices for haplotypes to keep
filtered_seqs_bin <- filtered_seqs_bin[haps_to_keep_indices]

# Recompute the haplotypes after filtering
haps_filtered <- haplotype(filtered_seqs_bin)

# Recompute the distance matrix and network
dist <- dist.dna(haps_filtered, "N")
net <- rmst(dist, B = 100)

nt.labs <- attr(net, "labels")

# Update size summary for filtered haplotypes

sz_filtered <- log2(summary(haps_filtered)) * 50
sz_filtered <- sz_filtered[nt.labs]

# Create a matrix of filtered sequences with updated row names
matrix_filtered <- as.matrix(filtered_seqs_bin)

# Ensure row names match sequence names
rownames(matrix_filtered) <- names(filtered_seqs_bin)

# Recalculate region matrix
regions_filtered <- haploFreq(matrix_filtered, split = "_", what = 3)

regions_filtered <- regions_filtered[nt.labs , ]

# Plot with the filtered data
colors <- colorRampPalette(brewer.pal(8, "Accent"))(15)

#png("mdr1.png" , units = "in" ,width = 10,height = 7  , res = 600)
plot(net, size = sz_filtered, scale.ratio = 50, cex = 0.01 , pie = regions_filtered, threshold = 0, col.link = "black", bg = colors)
legend("topright", colnames(regions_filtered), col=colors, pch=20, cex=1.5, pt.cex=3)
#legend("topleft" , legend = legendlable  , cex =1 , pt.cex =3  )

dev.off()

######DHFR##########
seqs <- readDNAStringSet("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/Pfdhfr_consensus_aligned.fa")
metadata<- read.csv("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/District_and_Sample_ID_Metadata_oct11.csv")

# Filter
dna_matrix <- as.matrix(seqs)
sample_threshold <- 0.1  # Allow up to 10% missing data in a sample
position_threshold <- 0.1  # Allow up to 10% missing data at a position

missing_per_sample <- rowMeans(dna_matrix == "N" | dna_matrix == "-")
filtered_matrix <- dna_matrix[missing_per_sample <= sample_threshold, ]
missing_per_position <- colMeans(filtered_matrix == "N" | filtered_matrix == "-")
filtered_matrix <- filtered_matrix[, missing_per_position <= position_threshold]

filtered_seqs <- DNAStringSet(apply(filtered_matrix, 1, paste0, collapse = ""))
writeXStringSet(filtered_seqs, "Dhfr_filtered_aligned_sequences.fasta")


#Add district to sample name
current_names <- names(filtered_seqs)
metadata_ordered <- metadata[match(current_names, metadata$Sample_ID), ]
new_names <- paste0(metadata_ordered$Sample_ID, "_", metadata_ordered$District)
names(filtered_seqs) <- new_names

filtered_seqs_bin <- as.DNAbin(filtered_seqs)
haps <- haplotype(filtered_seqs_bin)

dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))

# Set a threshold for rare haplotypes
hap_threshold <- 2 #Filter out haplotypes that occur fewer than 2 times
# Identify haplotypes to keep (those occurring >= threshold)
haps_to_keep <- names(sz)[sz >= hap_threshold]
haps_to_keep_indices <- unlist(attr(haps, "index")[sz >= hap_threshold])

# Subset the sequences using indices for haplotypes to keep
filtered_seqs_bin <- filtered_seqs_bin[haps_to_keep_indices]

# Recompute the haplotypes after filtering
haps_filtered <- haplotype(filtered_seqs_bin)

# Recompute the distance matrix and network
dist <- dist.dna(haps_filtered, "N")
net <- rmst(dist, B = 100)

nt.labs <- attr(net, "labels")

# Update size summary for filtered haplotypes

sz_filtered <- log2(summary(haps_filtered)) * 20
sz_filtered <- sz_filtered[nt.labs]

# Create a matrix of filtered sequences with updated row names
matrix_filtered <- as.matrix(filtered_seqs_bin)

# Ensure row names match sequence names
rownames(matrix_filtered) <- names(filtered_seqs_bin)

# Recalculate region matrix
regions_filtered <- haploFreq(matrix_filtered, split = "_", what = 3)

regions_filtered <- regions_filtered[nt.labs , ]

# Plot with the filtered data
colors <- colorRampPalette(brewer.pal(8, "Accent"))(15)

#png("dhfr.png" , units = "in" ,width = 14,height = 7  , res = 600)
plot(net, size = sz_filtered, scale.ratio = 50, cex = 0.1 , pie = regions_filtered, threshold = 0, col.link = "black", bg = colors)
legend("topright", colnames(regions_filtered), col=colors, pch=20, cex=1.5, pt.cex=3)
#legend("topleft" , legend = legendlable  , cex =1 , pt.cex =3  )

dev.off()

################ DHPS ##########


seqs <- readDNAStringSet("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/Pfdhps_consensus_aligned.fa")

metadata<- read.csv("C:/Users/Fitsum/Desktop/Ale.Haplotype/New folder/Files for running R code/Files for running R code/District_and_Sample_ID_Metadata_oct11.csv")

# Filter
dna_matrix <- as.matrix(seqs)
sample_threshold <- 0.1  # Allow up to 10% missing data in a sample
position_threshold <- 0.1  # Allow up to 10% missing data at a position

missing_per_sample <- rowMeans(dna_matrix == "N" | dna_matrix == "-")
filtered_matrix <- dna_matrix[missing_per_sample <= sample_threshold, ]
missing_per_position <- colMeans(filtered_matrix == "N" | filtered_matrix == "-")
filtered_matrix <- filtered_matrix[, missing_per_position <= position_threshold]

filtered_seqs <- DNAStringSet(apply(filtered_matrix, 1, paste0, collapse = ""))
writeXStringSet(filtered_seqs, "Dhps_filtered_aligned_sequences.fasta")


#Add district to sample name
current_names <- names(filtered_seqs)
metadata_ordered <- metadata[match(current_names, metadata$Sample_ID), ]
new_names <- paste0(metadata_ordered$Sample_ID, "_", metadata_ordered$District)
names(filtered_seqs) <- new_names

filtered_seqs_bin <- as.DNAbin(filtered_seqs)
haps <- haplotype(filtered_seqs_bin)

dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))

# Set a threshold for rare haplotypes
hap_threshold <- 2 #Filter out haplotypes that occur fewer than 2 times
# Identify haplotypes to keep (those occurring >= threshold)
haps_to_keep <- names(sz)[sz >= hap_threshold]
haps_to_keep_indices <- unlist(attr(haps, "index")[sz >= hap_threshold])

# Subset the sequences using indices for haplotypes to keep
filtered_seqs_bin <- filtered_seqs_bin[haps_to_keep_indices]

# Recompute the haplotypes after filtering
haps_filtered <- haplotype(filtered_seqs_bin)

# Recompute the distance matrix and network
dist <- dist.dna(haps_filtered, "N")
net <- rmst(dist, B = 100)

nt.labs <- attr(net, "labels")

# Update size summary for filtered haplotypes

sz_filtered <- log2(summary(haps_filtered)) * 10
sz_filtered <- sz_filtered[nt.labs]

# Create a matrix of filtered sequences with updated row names
matrix_filtered <- as.matrix(filtered_seqs_bin)

# Ensure row names match sequence names
rownames(matrix_filtered) <- names(filtered_seqs_bin)

# Recalculate region matrix
regions_filtered <- haploFreq(matrix_filtered, split = "_", what = 3)



regions_filtered <- regions_filtered[nt.labs , ]



# Plot with the filtered data
colors <- colorRampPalette(brewer.pal(8, "Accent"))(15)

png("dhps.png" , units = "in" ,width = 15,height = 7  , res = 600)
plot(net, size = sz_filtered, scale.ratio = 50, cex = 0.1 , pie = regions_filtered, threshold = 0, col.link = "black", bg = colors)
legend("topright", colnames(regions_filtered), col=colors, pch=20, cex=1.5, pt.cex=3)

dev.off()