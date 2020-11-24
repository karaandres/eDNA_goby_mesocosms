# Created: Oct. 28, 2019 by KJA
# Last modified: Sept 21, 2020 by KJA (kja68@cornell.edu)

# This script accomplishes the following: 
# Part 1: Calculate tissue allele frequencies
# Part 2: Calculate eDNA allele frequencies
# Part 3: Check for HWE and identify problem loci
# Part 4: PCA of tissue allele frequencies and eDNA read frequencies 
# Part 5: Heatmap of allele/read frequencies
# Part 6: Correlation between population allele/read frequencies 
# Part 7: Contributor estimation
# Part 8: Supplementary figures 
# Part 9: Check blanks for contamination

library(adegenet) # 
library(pegas) #
library(dplyr) # 
library(tidyr) #
library(ggplot2) #
library(gridExtra)
library(splitstackshape) #
library(data.table)
library(RColorBrewer) #
library(gplots) #
library(viridis) # 

#############################################################################################
############## Part 1: Calculate tissue allele frequencies ##################################
#############################################################################################

# Recode genotypes for missing data (< 10 reads/locus)
# Import read count and genotype matrix: col = samples, row = loci, cells = allele1/allele2:reads1,reads2
hap_genotype <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/hap_genotype.csv", header=TRUE)
hap_genotype <- data.frame(lapply(hap_genotype[,-2], gsub, pattern=".*:", replacement=""))   # Remove haplotype column, remove everything before ":" in each cell
hap_genotype <- hap_genotype[,-grep(pattern="Locus|e|B|CAY", colnames(hap_genotype))] # subset to only individual genotypes (i.e. no locus names, eDNA, or blanks)
hap_genotype_mat <- t(as.matrix(hap_genotype[,]))   # Transpose, turn into matrix, remove col and row headers

# Create a new matrix w/ the total number of reads per locus per sample
counts_total <- matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat)) # nrow = samples, ncol = loci
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  counts_total[i] <- sum(as.numeric(strsplit(hap_genotype_mat[i], ",")[[1]]))
}

# Create a new matrix that treats total # reads/locus <=10 as missing data
counts_less10 <- matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat))
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  if (counts_total[i]<=10) {
    counts_less10[i]=NA} else {
      counts_less10[i]=counts_total[i]
    }
}

# avg and sd read depth 
mean(apply(counts_less10, 2, function(x) {mean(x, na.rm = TRUE)}))
sd(apply(counts_less10, 2, function(x) {mean(x, na.rm = TRUE)}))

# Read in genotypes
genotypes <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/hap_genotype_matrix.csv", header=TRUE)
genotypes <- genotypes[-grep(pattern="e|B|CAY", genotypes$X),] # only retain for individuals (e.g. not eDNA or blanks)
genotypes <- as.matrix(genotypes[,-1]) # remove column of individual IDs

# If total count <10, locus gets "0|0"
genotypes_edited <-  matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat))
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  if (is.na(counts_less10[i]==TRUE)) {
    genotypes_edited[i] <- paste0(0, "|", 0)} else {
      genotypes_edited[i] = genotypes[i]
    }
}
genotypes_edited_df <- as.data.frame(genotypes_edited, row.names = colnames(hap_genotype)) # rows=individual IDs, columns=loci
colnames(genotypes_edited_df) <- colnames(genotypes)

# Add mesocosm numbers
genotypes_edited_df$mesocosm <- c(rep(c(1,2), each = 10), rep(3, 11), rep(c(4,5,6), each = 5), rep(c(7,8,9), each = 3), 10,11,12)
# write.csv(genotypes_edited_df, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/tissue_genotypes_less10.csv")

# split into alleles
genotypes_edited_df_split <- as.data.frame(cSplit(genotypes_edited_df[,-36], 1:ncol(genotypes_edited_df[,-36]), '|'))  # split each locus by "|"
genotypes_edited_df_split$mesocosm <- genotypes_edited_df$mesocosm
# write.csv(genotypes_edited_df_split, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/tissue_genotypes_less10_split.csv")

# Total genotype frequencies
locus_positions <- seq(1, ncol(genotypes_edited_df_split)-1, 2)  # starting column number for each locus
lnames <- colnames(genotypes_edited_df_split)
total_tissue_allele_freqs <- NULL
for (j in locus_positions) { # for each locus (35)
  alleles <- c(genotypes_edited_df_split[,j], genotypes_edited_df_split[,j+1]) # Combine 2 columns per locus
  alleles2 <- as.data.frame(table(alleles)) # count each allele at locus x
  alleles3 <- alleles2[alleles2$Freq!=0 & alleles2$alleles!=0,] # remove missing data (otherwise 0 would be counted in total number of alleles)
  alleles4 <- cbind(alleles3,alleles3$Freq/sum(alleles3$Freq)) # calculate frequencies
  if (length(alleles4) > 0) {
    output <- cbind(rep(lnames[j], nrow(alleles4)),alleles4) #combine j, locus name, and frequencies
    total_tissue_allele_freqs <<- rbind(total_tissue_allele_freqs,output)          
  }
}

colnames(total_tissue_allele_freqs) <- c("locus","haplotypes","count","tissue_freq")
total_tissue_allele_freqs$locus <- gsub(total_tissue_allele_freqs$locus, pattern = "_1", replacement = "") # remove subscript from locus names 
# write.csv(total_tissue_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/total_tissue_allele_freqs.csv")

# Genotype frequencies per mesocosm
# Subset genotypes for individuals by mesocosm (mesocosm 1-12)
locus_positions <- seq(1, ncol(genotypes_edited_df_split)-1, 2)  # starting column number for each locus
lnames <- colnames(genotypes_edited_df_split)
mesocosm_tissue_allele_freqs <- NULL
for (i in unique(genotypes_edited_df_split$mesocosm)){ # for each mesocosm (12)
  mesocosm_subset <- genotypes_edited_df_split[genotypes_edited_df_split$mesocosm==i,]
  for (j in locus_positions) { # for each locus (35)
    alleles <- c(mesocosm_subset[,j], mesocosm_subset[,j+1]) # combine 2 alleles per locus
    alleles2 <- as.data.frame(table(alleles)) # count each allele at locus x
    alleles3 <- alleles2[alleles2$Freq!=0 & alleles2$alleles!=0,] # remove missing data (otherwise 0 would be counted in total number of alleles)
    alleles4 <- cbind(alleles3,alleles3[,2]/sum(alleles3[,2])) # calculate frequencies
    mesocosm <- rep(paste(i), nrow(alleles4)) # mesocosm name 
    if (length(mesocosm) > 0) {
      output <- cbind(lnames[j],alleles4, mesocosm) #combine locus name, frequencies
      mesocosm_tissue_allele_freqs <<- rbind(mesocosm_tissue_allele_freqs,output) # add the new rows to the bottom of the data frame        
    }
  }
}

colnames(mesocosm_tissue_allele_freqs) <- c("locus","haplotypes","count","tissue_freq", "mesocosm")
mesocosm_tissue_allele_freqs$locus <- gsub(mesocosm_tissue_allele_freqs$locus, pattern = "_1", replacement = "") # remove subscript from locus names 
# write.csv(mesocosm_tissue_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_tissue_allele_freqs.csv")

#############################################################################################
############## Part 2: Calculate eDNA read frequencies ######################################
#############################################################################################

# Total eDNA read frequencies for all frequencies > 0.01
# Import read counts for all alleles at all loci
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus

# Begin loop
j = 1 # individual locus indexed from "filenames"
# create empty dataframe to populate with allele frequencies for each mesocosm
total_edna_allele_freqs <- data.frame(NULL)
# Create matrix with total number of reads/haplotype <=10 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
    file_mat <- as.matrix(file[,grep("e", colnames(file))])
    file_mat <- file_mat[,-grep("Haplotype", colnames(file_mat))]
    file_mat[file_mat < 10] <- 0 # replace counts <10 with 0
    file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
    file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
    file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
    edna_allele_freqs_temp <- data.frame(locus=rep(filenames[[j]], nrow(file_mat_norm)), allele=file[,1], count=rowSums(file_mat_norm, na.rm = TRUE))
    edna_allele_freqs_temp$edna_freq <- edna_allele_freqs_temp$count/sum(edna_allele_freqs_temp$count, na.rm = TRUE)
    total_edna_allele_freqs <- rbind(total_edna_allele_freqs, edna_allele_freqs_temp)
    j <- j+1
  }
total_edna_allele_freqs <- total_edna_allele_freqs[total_edna_allele_freqs$edna_freq>0,]
# write.csv(total_edna_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/total_edna_allele_freqs.csv")

# Check blanks
# create empty dataframe to populate with max reads for each controls
edna_control_reads <- data.frame(locus=names(dataFiles), BL1=NA, BL2=NA, EB=NA)
# Create matrix that removes haplotypes w/ <10 reads
x <- 1
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,grep("BL|EB", colnames(file))])
  file_mat[file_mat < 10] <- 0 # replace counts <10 with 0
  counts.less10_df <- as.data.frame(file_mat)
  colnames(counts.less10_df) <- colnames(file_mat)
  edna_control_reads[x,2:4] <- apply(counts.less10_df,2,sum,na.rm=TRUE)
  x <- x+1
}

# Read frequencies per mesocosm, theshold 0.001-0.1 and by per-locus allelic richness 
mesocosm_tissue_allele_freqs <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_tissue_allele_freqs.csv")
total_tissue_allele_freqs <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/total_tissue_allele_freqs.csv")
allelic_richness <- data.frame(table(total_tissue_allele_freqs$locus))
bins <- seq(from=0.005,to=0.05,length.out=20)
calculate_allele_freqs <- function(threshold){
  j=1
  edna_allele_freqs <- data.frame(NULL)
  for (file in dataFiles) { # for the reads at each locus for all samples 
    file_mat <- as.matrix(file[,grep("e", colnames(file))])
    file_mat <- file_mat[,-grep("Haplotype", colnames(file_mat))]
    file_mat[file_mat < 10] <- 0 # replace counts <10 with 0
    #file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
    NumAl <- nrow(total_tissue_allele_freqs[total_tissue_allele_freqs$locus==paste(filenames[[j]]),])
    #if (threshold=="variable_threshold"){
    #  file_mat_norm[file_mat_norm<bins[NumAl]] <- 0
    #} else file_mat_norm[file_mat_norm<threshold] <- 0
    #file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
    file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE)*100)) 
    #file_mat <- apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE)*100) # read count normalized per sample to 100 counts
    edna_allele_freqs_temp <- data.frame(locus_allele = NA) # empty dataframe
    k=1 # number (label) for each mesocosm
    for (i in seq(1,24,2)){ # for each mesocosm
      mesocosm <- data.frame(haplotypes=as.factor(file$Haplotype),
                             edna_sums=rowSums(file_mat[,c(i, i+1)], na.rm = T)) # combine read counts for 2 eDNA samples
      mesocosm$edna_freq <- mesocosm$edna_sums/sum(mesocosm$edna_sums, na.rm = TRUE) # calculate read frequency for each allele in combined eDNA samples 
      #file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
      #NumAl <- nrow(total_tissue_allele_freqs[total_tissue_allele_freqs$locus==paste(filenames[[j]]),])
      if (threshold=="variable_threshold"){
        mesocosm$edna_freq[mesocosm$edna_freq<bins[NumAl]] <- 0
      } else mesocosm$edna_freq[mesocosm$edna_freq<threshold] <- 0 # remove alleles below threshold
      # file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
      mesocosm$edna_freq <- mesocosm$edna_freq/sum(mesocosm$edna_freq, na.rm = TRUE) # re-calculate read frequency
      mesocosm <- merge(mesocosm, mesocosm_tissue_allele_freqs[mesocosm_tissue_allele_freqs$locus==paste(filenames[[j]]) & 
                                                                 mesocosm_tissue_allele_freqs$mesocosm==k,], by = "haplotypes", all = TRUE)
      mesocosm$tissue_freq[is.na(mesocosm$tissue_freq)] <- 0
      mesocosm$edna_freq[is.na(mesocosm$edna_freq)] <- 0
      mesocosm <-  mesocosm[mesocosm$edna_freq!=0 | mesocosm$tissue_freq!=0,]
      if (nrow(mesocosm) > 0){
        mesocosm <- data.frame(locus_allele = paste(filenames[[j]],"_", mesocosm$haplotypes,sep = ""), mesocosm$tissue_freq, mesocosm$edna_freq) # combine locus/allele names, eDNA read frequencies
      } else {mesocosm <- data.frame(locus_allele = NA, 0,0)}
      colnames(mesocosm) <- c("locus_allele", paste("tissue_",k,sep=""), paste("edna_",k,sep=""))
      edna_allele_freqs_temp <- merge(edna_allele_freqs_temp, mesocosm, by = "locus_allele", all = TRUE)
      k=k+1
    }
    edna_allele_freqs <-  rbind(edna_allele_freqs, edna_allele_freqs_temp)
    j=j+1
  }
  edna_allele_freqs[is.na(edna_allele_freqs)] <- 0 # replace NAs with 0
  edna_allele_freqs <- edna_allele_freqs[which(rowSums(edna_allele_freqs[,2:25]) > 0),] # remove any rows where all frequencies == 0
  return(edna_allele_freqs)
}

mesocosm_edna_allele_freqs_0.1 <- calculate_allele_freqs(0.1)
# write.csv(mesocosm_edna_allele_freqs_0.1, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_edna_allele_freqs_0.1.csv")
mesocosm_edna_allele_freqs_0.01 <- calculate_allele_freqs(0.01)
# write.csv(mesocosm_edna_allele_freqs_0.01, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_edna_allele_freqs_0.01.csv")
mesocosm_edna_allele_freqs_0.001 <- calculate_allele_freqs(0.001)
# write.csv(mesocosm_edna_allele_freqs_0.001, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_edna_allele_freqs_0.001.csv")
mesocosm_edna_allele_freqs_NumAl <- calculate_allele_freqs("variable_threshold")
# write.csv(mesocosm_edna_allele_freqs_NumAl, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_edna_allele_freqs_NumAl.csv")

#############################################################################################
############## Part 3: Check for HWE and identify problem loci  #############################
#############################################################################################
# Exclude alleles not in HWE and with heterozygote excess 
geno <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/tissue_genotypes_less10.csv")
rownames(geno) <- geno$X # rownames are sample names
geno <- geno[,!names(geno) %in% "X"] # remove column of sample names
geno[geno == "0|0"] <- "NA|NA" # replace missing data "0|0" with NAs
geno.obj <- df2genind(geno[,1:35], sep="\\|", ploidy=2, type="codom") # "|" must be preceded by double backslashes
summary(geno.obj)
hw.output <- data.frame(hw.test(geno.obj, 10000))
hw.output$observed_heterozygosity <- summary(geno.obj)[[6]]
hw.output$expected_heterozygosity <- summary(geno.obj)[[7]]
hw.output[hw.output$Pr.exact< 0.05/35 & hw.output$observed>hw.output$expected,] # bonferroni corrected pvalue
# remove loci 155, 248, 361, 726 -- deviation from HWE and heterzygote excess
# also l1103 and 1531 because of 1.00 heterozygosity

apply(genotypes_edited_df_split, 2, function(x) {sum(x==0)}) # number of alleles with missing data per locus (column)
# also remove 185, maybe 344 and 422 due to poor amplification

patterns <- c("Nmel185","Nmel155","Nmel248","Nmel361","Nmel726", "Nmel1103","Nmel1531") # loci to remove

# Calculate number of alleles with missing data per individual and locus
genotypes_edited_df_split_HWE <- genotypes_edited_df_split[, !grepl(paste(patterns, collapse="|"), colnames(genotypes_edited_df_split))]
apply(genotypes_edited_df_split_HWE, 1, function(x) {sum(x==0)}) # number of alleles with missing data per individual (row)
# maximum of 2 alleles missing 

allelic_richness <- data.frame(allele = NA, allelic_rich = NA)
for (i in seq(1,ncol(genotypes_edited_df_split_HWE)-1,2)){
  uniq_alleles <- c(genotypes_edited_df_split_HWE[,i], genotypes_edited_df_split_HWE[,i+1])
  allele <- colnames(genotypes_edited_df_split_HWE)[i]
  allelic_rich <- length(unique(uniq_alleles[uniq_alleles>0]))
  allelic_richness_temp <- cbind(allele, allelic_rich)
  allelic_richness <- rbind(allelic_richness, allelic_richness_temp) 
}
allelic_richness$allelic_rich <- as.numeric(allelic_richness$allelic_rich)
mean(allelic_richness$allelic_rich, na.rm = TRUE) # mean alleic rich = 8.61

# Remove problem loci from datasets: 
remove_loci <- function(x) filter(x, !grepl(paste(patterns, collapse="|"), locus))
total_tissue_allele_freqs_HWE <- remove_loci(total_tissue_allele_freqs)
total_edna_allele_freqs_HWE <- remove_loci(total_edna_allele_freqs)
remove_loci <- function(x) filter(x, !grepl(paste(patterns, collapse="|"), locus_allele))
mesocosm_edna_allele_freqs_0.1_HWE <- remove_loci(mesocosm_edna_allele_freqs_0.1)
mesocosm_edna_allele_freqs_0.01_HWE <- remove_loci(mesocosm_edna_allele_freqs_0.01)
mesocosm_edna_allele_freqs_0.001_HWE <- remove_loci(mesocosm_edna_allele_freqs_0.001)
mesocosm_edna_allele_freqs_NumAl_HWE <- remove_loci(mesocosm_edna_allele_freqs_NumAl)

#############################################################################################
############## Part 4: PCA of tissue allele frequencies and eDNA read frequencies ###########
#############################################################################################
# PCA of genotype frequencies and eDNA read counts with threshold based on number of alleles 
PCA_data <- mesocosm_edna_allele_freqs_NumAl
n <- PCA_data$locus_allele # save vector of allele names
PCA_data <- as.data.frame(t(subset(PCA_data, select = -locus_allele))) # transpose data frame so each mesocosm is in a row, remove allele column
colnames(PCA_data) <- n # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], type = rownames(my.pca$x), 
                          mesocosm = rep(1:12, each = 2), sample = rep(c("Tissue", "eDNA"), 6), 
                          trtmt = rep(1:4, each=6), triplicate = rep(c("a","b","c"), each=2,2))
my_palette <- c("#762A83","#009E73","#F46D43","#0072B2")

my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=factor(trtmt), 
                                   shape = factor(triplicate),
                                   fill = factor(ifelse(sample == "Tissue", trtmt, "eDNA")))) +
  scale_color_manual(name = "trtmt", values = my_palette) +
  scale_shape_manual(name = "triplicate", values=c(24,22,21)) +
  scale_fill_manual(name = "sample", values = c(my_palette,"white")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold")) +
  geom_point(size = 5, stroke = 2)
# setEPS()
# postscript("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/Fig1_pca_num_alleles.eps",
#             width = 10, height = 8)
my_plot
# dev.off()


# PCA of genotype frequencies and eDNA read counts at threshold 0.01 
PCA_data <- mesocosm_edna_allele_freqs_0.01_HWE
n <- PCA_data$locus_allele # save vector of allele names
PCA_data <- as.data.frame(t(subset(PCA_data, select = -locus_allele))) # transpose data frame so each mesocosm is in a row, remove allele column
colnames(PCA_data) <- n # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], type = rownames(my.pca$x), 
                          mesocosm = rep(1:12, each = 2), sample = rep(c("Tissue", "eDNA"), 6), 
                          trtmt = rep(1:4, each=6), triplicate = rep(c("a","b","c"), each=2,2))
my_palette <- viridis(4)
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=factor(trtmt), 
                                   shape = factor(triplicate),
                                   fill = factor(ifelse(sample == "Tissue", trtmt, "eDNA")))) +
  scale_color_manual(name = "trtmt", values = my_palette) +
  scale_shape_manual(name = "triplicate", values=c(24,22,21)) +
  scale_fill_manual(name = "sample", values = c(my_palette,"white")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(size = 5, stroke = 2)
my_plot
#ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/Fig1_pca.eps"), plot=my_plot, dpi = 300,  width = 10, height = 8, units = "in")   


#############################################################################################
############# Part 5: Heatmap of allele/read frequencies ####################################
#############################################################################################
PCA_data <- mesocosm_edna_allele_freqs_0.01_HWE
n <- PCA_data$locus_allele # save vector of allele names
PCA_data <- as.data.frame(t(subset(PCA_data, select = -locus_allele))) # transpose data frame so each mesocosm is in a row, remove allele column
colnames(PCA_data) <- n # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data)
pca_summary <- summary(my.pca)
pca_heat <- my.pca$x[c(24:1), ]
rownames(pca_heat) <- c("Tissue 1a", "eDNA 1a", "Tissue 1b", "eDNA 1b", "Tissue 1c", "eDNA 1c",
                        "Tissue 3a", "eDNA 3a", "Tissue 3b", "eDNA 3b", "Tissue 3c", "eDNA 3c",
                        "Tissue 5a", "eDNA 5a", "Tissue 5b", "eDNA 5b", "Tissue 5c", "eDNA 5c",
                        "Tissue 10a", "eDNA 10a", "Tissue 10b", "eDNA 10b", "Tissue 11c", "eDNA 11c")
cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))(50)
cols2 <- rep(viridis(4), each = 6)
heat <- as.matrix(dist(pca_heat))
heat[upper.tri(heat)] <- NA

pdf('/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/Fig1_pca_heatmap.pdf',width = 6, height = 5.5)
heatmap.2(heat, Rowv = FALSE, Colv = FALSE, symm = TRUE,
          dendrogram = "none", trace = "none", scale = "none",
          col = cols, colRow = rev(cols2), colCol = rev(cols2), srtCol=45,
          density.info = "none", key.title = "", key.xlab = "Euclidian distance")
dev.off()

#############################################################################################
############## Part 6: Correlation between population allele/read frequencies ###############
#############################################################################################
genotype_freqs <- unite(total_tissue_allele_freqs_HWE, "allele", c(locus, haplotypes), sep = "_",remove = FALSE)
edna_freqs <- unite(total_edna_allele_freqs_HWE, "allele", c(locus, allele), sep = "_", remove = FALSE)
total_allele_freqs <- merge(genotype_freqs, edna_freqs, by = "allele", all = TRUE)
total_allele_freqs$tissue_freq[is.na(total_allele_freqs$tissue_freq)] <- 0
total_allele_freqs$edna_freq[is.na(total_allele_freqs$edna_freq)] <- 0
total_allele_freqs$locus <- ifelse(is.na(total_allele_freqs$locus.y), total_allele_freqs$locus.x, total_allele_freqs$locus.y)

cor.test(total_allele_freqs$tissue_freq, total_allele_freqs$edna_freq)
cors <- ddply(total_allele_freqs, "locus", summarise, cor=round(cor(tissue_freq,edna_freq,use="complete.obs"), 2))
min(cors$cor)
max(cors$cor)

# png('/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/mesocosm_correlation.png',width = 6, height = 6, units = "in", res = 400)
plot(total_allele_freqs$tissue_freq, total_allele_freqs$edna_freq, xlim = c(0,1), ylim = c(0,1),
     xlab = "Individual allele frequencies", ylab = "eDNA allele frequencies", cex.lab = 1.5, cex = 1.5, lwd = 1.5,
     col = "darkblue")
# dev.off()

p <- ggplot(total_allele_freqs, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylim(0,1) + xlim(0,1) + ylab("Mesocosm eDNA allele frequency") + xlab("Genotyped tissue allele frequency") +
  scale_color_manual(values=viridis(28)) + theme_bw() + theme(legend.position="none")
# + geom_text(aes(label=locus),hjust=0, vjust=0)
#ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/Fig1_mesocosm_allele_freq_corr.eps"), plot=p, dpi = 300,  width = 5, height = 5, units = "in")   

#############################################################################################
############## Part 7: Contributor estimation ###############################################
#############################################################################################

# following Weir, Haned: sum.p = 0 if # alleles exceeds 2x
# Suresh's code: speed up computation time when # of alleles is large
MixtureLikelihood_Parsed <- function(x,p.v){
  require(RcppAlgos)
  if(length(p.v)>(2*x)) {sum.p <- 0} 
  else {
    counter <- 0
    if (length(p.v)>0){
      sum.p <- sum(p.v)^(2*x)
      if(length(p.v)>1){ # only conduct if > 1 allele is observed
        for(i in 1:(length(p.v)-1)){
          counter <- counter+1
          temp.combo <- comboGeneral(v=length(p.v), m = length(p.v)-i, repetition = FALSE, Parallel=TRUE,nThreads=10) #; i; dim(temp.combo);i<-i+1
          if(nrow(temp.combo)>15e6){
            ix <- 999999		
            row.sums.pow2x <- 1:(floor(nrow(temp.combo)/ix)+1)
            for(m in 1:length(row.sums.pow2x)){
              temp.m <- temp.combo[ ((m-1)*ix+1) : min(m*ix,nrow(temp.combo)), ]
              temp.combo.v <- c(temp.m)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.m)[1],nc=dim(temp.m)[2],byrow=F)
              row.sums.pow2x[m] <- sum(rowSums(p.v.m)^(2*x)) # reference as folows
            } # end m loop
            sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) } else 
            {
              # non-parsed
              temp.combo.v <- c(temp.combo)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
              row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as folows
              sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) 
            } # end parsing option
        }# end i loop
      } # end if
    } else sum.p <- NA # if no alleles observed, assign NA
  }
  return(sum.p)
} # end function

# Example for a single locus
p.v. <- c(0.3, 0.1, 0.1, 0.05) # 4 distinct alleles observed at the locus from a DNA mixture with associated population allele frequencies
x. <- 3 # Putative number of contributors
MixtureLikelihood_Parsed(x=x.,p.v=p.v.)

# Multi-locus function across 1:x number of putative contributors
MixtureLikelihood2 <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Parsed(x=j,p.v=unlist(loci[i])) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

# Example for 1:100 putative contributors across 4 loci, 1 with no alleles 
y. <- 100
loci. <- list(locus1 =  c(0.3, 0.1, 0.1, 0.05, 0.05, 0.03, 0.01, 0.01),
              locus2 = c(0.3, 0.2, 0.1, 0.1, 0.05),
              locus3= c(0.4, 0.1, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01),
              locus4 = NA)
MixtureLikelihood2(y., loci.)

# 7b: Bias in estimator with varying thresholds, population = tissue allele freqs
# Total allele frequencies from all genotyped gobies 
pop_allele_freqs_tissue <- unite(total_tissue_allele_freqs_HWE, "locus_allele", c(locus, haplotypes), sep = "_")
pop_allele_freqs_tissue <- pop_allele_freqs_tissue[,c("locus_allele","tissue_freq")]
names(pop_allele_freqs_tissue)[2] <- "pop_freq"

contributor_estimation <- function(sample_alleles, population_alleles){
  sample_long <- gather(sample_alleles, "sample", "freq", tissue_1:edna_12)
  sample_long <- sample_long[sample_long$freq>0,]
  combined_allele_freqs <- merge(sample_long, population_alleles, by = "locus_allele")
  combined_allele_freqs_split <- cSplit(combined_allele_freqs, "locus_allele", '_', drop=FALSE)
  colnames(combined_allele_freqs_split)[c(5,6)] <- c("locus", "allele")
  z <- 1
  contrib_estim_mat <- matrix(nrow = 24, ncol = 2)
  for (i in unique(combined_allele_freqs_split$sample)){
    mesocosm_i <- combined_allele_freqs_split[combined_allele_freqs_split$sample == i]
    loci.list <- split(mesocosm_i, mesocosm_i$locus)
    loci.list <- lapply(loci.list, function(x) {x$pop_freq})
    contrib_estim_mat[z,1] <-  mesocosm_i$sample[1]
    contrib_estim_mat[z,2] <- MixtureLikelihood2(y = 50, loci = loci.list)$Contributors
    z <- z+1
  }
  contrib_estim <- as.data.frame(contrib_estim_mat)
  colnames(contrib_estim) = c("sample", "estimation")
  contrib_estim <- contrib_estim[order(contrib_estim$sample),]
  contrib_estim <- contrib_estim[c(13,17:24,14,15,16, 1,5:12,2,3,4),]
  return(contrib_estim)
}

mesocosm_edna_allele_freqs_0.001_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.001_HWE, pop_allele_freqs_tissue) 
mesocosm_edna_allele_freqs_0.01_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.01_HWE, pop_allele_freqs_tissue)
mesocosm_edna_allele_freqs_0.1_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.1_HWE, pop_allele_freqs_tissue)
mesocosm_edna_allele_freqs_NumAl_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_NumAl_HWE, pop_allele_freqs_tissue)

contrib_estim <- rbind(mesocosm_edna_allele_freqs_0.001_ContEst, mesocosm_edna_allele_freqs_0.01_ContEst[13:24,],
                       mesocosm_edna_allele_freqs_0.1_ContEst[13:24,], mesocosm_edna_allele_freqs_NumAl_ContEst[13:24,])
contrib_estim$true_cont <- as.numeric(rep(c(10,10,11,5,5,5,3,3,3,1,1,1), times = 5)) # true numbers of individuals in mesocosms
contrib_estim$bias <- as.numeric(as.character(contrib_estim$estimation)) - contrib_estim$true_cont # calculate bias in estimation
contrib_estim$num_cont <- as.numeric(rep(c(10,10,10,5,5,5,3,3,3,1,1,1), times = 5)) # mesocosm treatment
contrib_estim$threshold <- as.factor(rep(c("Genotypes", "0.001", "0.01", "0.1", "NumAl"), each = 12)) 
contrib_estim$threshold <- factor(contrib_estim$threshold, levels = c("Genotypes", "0.001", "0.01", "0.1", "NumAl"))
contrib_estim$origin = rep(c("Tissue","eDNA"), c(12, 48))
contrib_estim$triplicate = rep(c("a","b","c"), 20)

# make plot
my_palette <- rev(viridis(4))
p <- ggplot(contrib_estim, aes(factor(num_cont), bias, color = factor(num_cont), shape = factor(triplicate), 
                               fill = factor(ifelse(origin == "Tissue", factor(num_cont), "eDNA")))) +
  geom_hline(yintercept=0) + ylim(-11,11) + 
  geom_hline(yintercept=c(-5,5), linetype = "dashed", color = "grey87", lwd = 0.5) +
  scale_color_manual(name = "true_cont", values = my_palette) +
  scale_shape_manual(name = "triplicate", values=c(24,22,21)) +
  scale_fill_manual(name = "origin", values = c(my_palette,"white")) +
  geom_jitter(width = 0.1, height = 0, size = 3) + 
  xlab("Number of contributors") + ylab("Bias (Estimate - True)") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +
  facet_grid(threshold ~ .)
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/contrib_estimation_tissue.eps"), plot=p, dpi = 300,  width = 3.5, height = 8, units = "in")   

# 7c: Bias in estimator with varying thresholds, population = eDNA read freqs (relative abundance)
# Total allele frequencies from all genotyped gobies 
pop_allele_freqs_edna <- unite(total_edna_allele_freqs_HWE, "locus_allele", c(locus, allele), sep = "_")
pop_allele_freqs_edna <- pop_allele_freqs_edna[pop_allele_freqs_edna$edna_freq>0.01,]
pop_allele_freqs_edna <- pop_allele_freqs_edna[-1,c("locus_allele","edna_freq")]
names(pop_allele_freqs_edna)[2] <- "pop_freq"

mesocosm_edna_allele_freqs_0.001_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.001_HWE, pop_allele_freqs_edna) 
mesocosm_edna_allele_freqs_0.01_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.01_HWE, pop_allele_freqs_edna)
mesocosm_edna_allele_freqs_0.1_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_0.1_HWE, pop_allele_freqs_edna)
mesocosm_edna_allele_freqs_NumAl_ContEst <- contributor_estimation(mesocosm_edna_allele_freqs_NumAl_HWE, pop_allele_freqs_edna)

contrib_estim <- rbind(mesocosm_edna_allele_freqs_0.001_ContEst, mesocosm_edna_allele_freqs_0.01_ContEst[13:24,],
                       mesocosm_edna_allele_freqs_0.1_ContEst[13:24,], mesocosm_edna_allele_freqs_NumAl_ContEst[13:24,])
contrib_estim$true_cont <- as.numeric(rep(c(10,10,11,5,5,5,3,3,3,1,1,1), times = 5)) # true numbers of individuals in mesocosms
contrib_estim$bias <- as.numeric(as.character(contrib_estim$estimation)) - contrib_estim$true_cont # calculate bias in estimation
contrib_estim$num_cont <- as.numeric(rep(c(10,10,10,5,5,5,3,3,3,1,1,1), times = 5)) # mesocosm treatment
contrib_estim$threshold <- as.factor(rep(c("Genotypes", "0.001", "0.01", "0.1", "NumAl"), each = 12)) 
contrib_estim$threshold <- factor(contrib_estim$threshold, levels = c("Genotypes", "0.001", "0.01", "0.1", "NumAl"))
contrib_estim$origin = rep(c("Tissue","eDNA"), c(12, 48))
contrib_estim$triplicate = rep(c("a","b","c"), 20)

# make plot
my_palette <-  rev(viridis(4))
p <- ggplot(contrib_estim, aes(factor(num_cont), bias, color = factor(num_cont), shape = factor(triplicate), 
                               fill = factor(ifelse(origin == "Tissue", factor(num_cont), "eDNA")))) +
  geom_hline(yintercept=0) + ylim(-11,11) + 
  geom_hline(yintercept=c(-5,5), linetype = "dashed", color = "grey87", lwd = 0.5) +
  scale_color_manual(name = "true_cont", values = my_palette) +
  scale_shape_manual(name = "triplicate", values=c(24,22,21)) +
  scale_fill_manual(name = "origin", values = c(my_palette,"white")) +
  geom_jitter(width = 0.1, height = 0, size = 3) + 
  xlab("Number of contributors") + ylab("Bias (Estimate - True)") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15,face="bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +
  facet_grid(threshold ~ .)
p
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/contrib_estimation_edna.eps"), plot=p, dpi = 300,  width = 3.5, height = 8, units = "in")   

# 7d: Contributor estimator with simulated mesocosm mixtures
contributor_estimation_sim <- function(sample_alleles, population_alleles,threshold){
  contrib_estim_mat <- matrix(nrow = 1100, ncol = 2)
  z=1
  for (i in 2:12){
    for (j in 1:100){
      dat <- sample_alleles[,grep("locus_allele|edna",colnames(sample_alleles))]
      subsample <- data.frame(locus_allele = dat$locus_allele, dat[, sample(2:ncol(dat), i, replace = FALSE)])
      true_cont <- as.vector(colnames(subsample)[-1])
      true_cont <- gsub("edna_10|edna_11|edna_12",1,true_cont)
      true_cont <- gsub("edna_7|edna_8|edna_9",3,true_cont)
      true_cont <- gsub("edna_4|edna_5|edna_6",5,true_cont)
      true_cont <- gsub("edna_3",11,true_cont)
      true_cont <- gsub("edna_1|edna_2",10,true_cont)
      n_true_cont <- sum(as.numeric(true_cont))
      subsample$edna_freq <- rowSums(subsample[,-1])/sum(rowSums(subsample[,-1]))
      subsample <- subsample[subsample$edna_freq>0,]
      combined_allele_freqs <- merge(subsample[,grep("locus_allele|edna_freq",colnames(subsample))], population_alleles, by = "locus_allele")
      combined_allele_freqs_split <- cSplit(combined_allele_freqs, "locus_allele", '_')
      colnames(combined_allele_freqs_split)[3:4] <- c("locus", "allele")
      loci.list <- split(combined_allele_freqs_split, combined_allele_freqs_split$locus)
      loci.list <- lapply(loci.list, function(x) {x$pop_freq})
      contrib_estim_mat[z,1] <- n_true_cont
      contrib_estim_mat[z,2] <- (MixtureLikelihood2(y = 100, loci = loci.list))$Contributors
      print(contrib_estim_mat)
      z=z+1
    } 
  }
  return(contrib_estim_mat)
}

cont_est_sim_0.001 <- contributor_estimation_sim(mesocosm_edna_allele_freqs_0.001_HWE, pop_allele_freqs_tissue)
cont_est_sim_0.01 <- contributor_estimation_sim(mesocosm_edna_allele_freqs_0.01_HWE, pop_allele_freqs_tissue)
cont_est_sim_0.1 <- contributor_estimation_sim(mesocosm_edna_allele_freqs_0.1_HWE, pop_allele_freqs_tissue)
cont_est_sim_NumAl <- contributor_estimation_sim(mesocosm_edna_allele_freqs_NumAl_HWE, pop_allele_freqs_tissue)

# png(filename = "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/cont_estim_sim_0.001.png", width = 5, height = 5, res = 300, units = "in")
cont_est_sim_0.001 <- as.data.frame(cont_est_sim_0.001)
colnames(cont_est_sim_0.001) <- c("True", "Estimated")
par(mar=c(5,5,5,5))
plot(cont_est_sim_0.001$True, cont_est_sim_0.001$Estimated, col = "darkgray", 
     pch = 1, cex = 2, ylab = "Estimated", xlab = "True", ylim = c(0,75), 
     xlim = c(0,75), cex.lab = 2)
abline(0, 1)
# dev.off()

# png(filename = "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/cont_estim_sim_0.01.png", width = 5, height = 5, res = 300, units = "in")
cont_est_sim_0.01 <- as.data.frame(cont_est_sim_0.01)
colnames(cont_est_sim_0.01) <- c("True", "Estimated")
par(mar=c(5,5,5,5))
plot(cont_est_sim_0.01$True, cont_est_sim_0.01$Estimated, col = "darkgray", 
     pch = 1, cex = 2, ylab = "Estimated", xlab = "True", 
     ylim = c(0,75), xlim = c(0,75), cex.lab = 2)
abline(0, 1)
# dev.off()

# png(filename = "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/cont_estim_sim_0.1.png", width = 5, height = 5, res = 300, units = "in")
cont_est_sim_0.1 <- as.data.frame(cont_est_sim_0.1)
colnames(cont_est_sim_0.1) <- c("True", "Estimated")
par(mar=c(5,5,5,5))
plot(cont_est_sim_0.1$True, cont_est_sim_0.1$Estimated, col = "darkgray", 
     pch = 1, cex = 2, ylab = "Estimated", xlab = "True", 
     ylim = c(0,75), xlim = c(0,75), cex.lab = 2)
abline(0, 1)
# dev.off()

# png(filename = "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/cont_estim_sim_NumAl.png", width = 5, height = 5, res = 300, units = "in")
cont_est_sim_NumAl <- as.data.frame(cont_est_sim_NumAl)
colnames(cont_est_sim_NumAl) <- c("True", "Estimated")
par(mar=c(5,5,5,5))
plot(cont_est_sim_NumAl$True, cont_est_sim_NumAl$Estimated, col = "darkgray", 
     pch = 1, cex = 2, ylab = "Estimated", xlab = "True", 
     ylim = c(0,75), xlim = c(0,75), cex.lab = 2)
abline(0, 1)
# dev.off()

#############################################################################################
############# Part 8: Supplementary figures #################################################
#############################################################################################
# Import read counts for all alleles at all loci
# duplicate sample comparison: eDNA read frequencies for all frequencies > 0.01
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus

# Read frequencies per mesocosm replicate, theshold 0.01
j=1
edna_replicate_counts <- data.frame(locus = NA, allele = NA, eDNA10_1A = NA, eDNA10_1B = NA, eDNA10_2A = NA, eDNA10_2B = NA, eDNA10_3A = NA, eDNA10_3B = NA,
                              eDNA5_1A = NA, eDNA5_1B = NA, eDNA5_2A = NA, eDNA5_2B = NA, eDNA5_3A = NA, eDNA5_3B = NA, 
                              eDNA3_1A = NA, eDNA3_1B = NA, eDNA3_2A = NA, eDNA3_2B = NA, eDNA3_3A = NA, eDNA3_3B = NA, 
                              eDNA1_1A = NA, eDNA1_1B = NA, eDNA1_2A = NA, eDNA1_2B = NA, eDNA1_3A = NA, eDNA1_3B = NA)
edna_replicate_freqs <- data.frame(locus = NA, allele = NA, eDNA10_1A = NA, eDNA10_1B = NA, eDNA10_2A = NA, eDNA10_2B = NA, eDNA10_3A = NA, eDNA10_3B = NA,
                              eDNA5_1A = NA, eDNA5_1B = NA, eDNA5_2A = NA, eDNA5_2B = NA, eDNA5_3A = NA, eDNA5_3B = NA, 
                              eDNA3_1A = NA, eDNA3_1B = NA, eDNA3_2A = NA, eDNA3_2B = NA, eDNA3_3A = NA, eDNA3_3B = NA, 
                              eDNA1_1A = NA, eDNA1_1B = NA, eDNA1_2A = NA, eDNA1_2B = NA, eDNA1_3A = NA, eDNA1_3B = NA)
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,grep("e", colnames(file))])
  file_mat <- file_mat[,-grep("Haplotype", colnames(file_mat))]
  file_mat[file_mat <= 1] <- 0 # replace counts <10 with 0
  file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) if (sum(x)>0){
                                                      x/sum(x, na.rm = TRUE)} else {x=0}))
  file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
  file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
  edna_replicate_counts_temp <- data.frame(rep(filenames[[j]], nrow(file_mat)), file[,1], file_mat)
  edna_replicate_counts_temp <- edna_replicate_counts_temp[which(rowSums(edna_replicate_counts_temp[,-c(1:2)]) > 0),] # remove any rows where all frequencies == 0
  colnames(edna_replicate_counts_temp) <- colnames(edna_replicate_counts)
  edna_replicate_freqs_temp <- data.frame(rep(filenames[[j]], nrow(file_mat_norm)), file[,1], file_mat_norm)
  edna_replicate_freqs_temp <- edna_replicate_freqs_temp[which(rowSums(edna_replicate_freqs_temp[,-c(1:2)]) > 0),] # remove any rows where all frequencies == 0
  colnames(edna_replicate_freqs_temp) <- colnames(edna_replicate_freqs)
  edna_replicate_counts <- rbind(edna_replicate_counts, edna_replicate_counts_temp)
  edna_replicate_freqs <- rbind(edna_replicate_freqs, edna_replicate_freqs_temp)
  j<-j+1
}

# Figure S1: correlation between read counts
edna_replicate_counts_1 <- edna_replicate_counts[,c(1,2,3,5,7,9,11,13,15,17,19,21,23,25)]
edna_replicate_counts_2 <- edna_replicate_counts[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26)]
edna_replicate_counts_1_long <- gather(edna_replicate_counts_1, sample, reads_1, eDNA10_1A:eDNA1_3A, factor_key=TRUE)
edna_replicate_counts_2_long <- gather(edna_replicate_counts_2, sample, reads_2, eDNA10_1B:eDNA1_3B, factor_key=TRUE)
edna_replicate_counts_1_long <- as.data.frame(cSplit(edna_replicate_counts_1_long, "sample", '_', drop=FALSE))
edna_replicate_counts_2_long <- as.data.frame(cSplit(edna_replicate_counts_2_long, "sample", '_', drop=FALSE))
names(edna_replicate_counts_1_long) <- c("locus_1","allele_1","sample_1","reads_1","mesocosm_1","replicate_1")
names(edna_replicate_counts_2_long) <- c("locus_2","allele_2","sample_2","reads_2","mesocosm_2","replicate_2")
edna_replicate_counts_long <- cbind(edna_replicate_counts_1_long, edna_replicate_counts_2_long)
edna_replicate_counts_long$locus_1 <- as.factor(edna_replicate_counts_long$locus_1)
edna_replicate_counts_long$mesocosm_1 <- factor(edna_replicate_counts_long$mesocosm_1, levels=c("eDNA1","eDNA3","eDNA5","eDNA10"))
remove_loci <- function(x) filter(x, !grepl(paste(patterns, collapse="|"), locus_1))
edna_replicate_counts_long <- remove_loci(edna_replicate_counts_long)
edna_replicate_counts_long <- edna_replicate_counts_long[!is.na(edna_replicate_counts_long$locus_1),]
cors <- ddply(edna_replicate_counts_long, c("mesocosm_1", "replicate_1"), summarise, cor = round(cor(reads_1, reads_2, use="complete.obs"), 2))
cors$locus_1 <- "Nmel1132"
p <- ggplot(edna_replicate_counts_long, aes(reads_1, reads_2, col=locus_1)) +
        geom_point() + ylim(0,5000) + xlim(0,5000) + xlab("eDNA sample 1") + ylab("eDNA sample 2") +
        geom_abline(slope=1, intercept=0, colour="gray",linetype=2) +
        scale_color_manual(values=viridis(28)) + theme_bw()
FigS1 <- p + facet_grid(vars(mesocosm_1), vars(replicate_1)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=1000, y=4500)
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS1_duplicate_read_counts_corr.eps"), plot=FigS1, dpi = 600)   

# Figure S2: correlation between read frequencies 
edna_replicate_freqs_1 <- edna_replicate_freqs[,c(1,2,3,5,7,9,11,13,15,17,19,21,23,25)]
edna_replicate_freqs_2 <- edna_replicate_freqs[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26)]
edna_replicate_freqs_1_long <- gather(edna_replicate_freqs_1, sample, reads_1, eDNA10_1A:eDNA1_3A, factor_key=TRUE)
edna_replicate_freqs_2_long <- gather(edna_replicate_freqs_2, sample, reads_2, eDNA10_1B:eDNA1_3B, factor_key=TRUE)
edna_replicate_freqs_1_long <- as.data.frame(cSplit(edna_replicate_freqs_1_long, "sample", '_', drop=FALSE))
edna_replicate_freqs_2_long <- as.data.frame(cSplit(edna_replicate_freqs_2_long, "sample", '_', drop=FALSE))
names(edna_replicate_freqs_1_long) <- c("locus_1","allele_1","sample_1","reads_1","mesocosm_1","replicate_1")
names(edna_replicate_freqs_2_long) <- c("locus_2","allele_2","sample_2","reads_2","mesocosm_2","replicate_2")
edna_replicate_freqs_long <- cbind(edna_replicate_freqs_1_long, edna_replicate_freqs_2_long)
edna_replicate_freqs_long$locus_1 <- as.factor(edna_replicate_freqs_long$locus_1)
edna_replicate_freqs_long$mesocosm_1 <- factor(edna_replicate_freqs_long$mesocosm_1, levels=c("eDNA1","eDNA3","eDNA5","eDNA10"))
remove_loci <- function(x) filter(x, !grepl(paste(patterns, collapse="|"), locus_1))
edna_replicate_freqs_long <- remove_loci(edna_replicate_freqs_long)
edna_replicate_freqs_long <- edna_replicate_freqs_long[!is.na(edna_replicate_freqs_long$locus_1),]
options(scipen = 999)
cors <- ddply(edna_replicate_freqs_long, c("mesocosm_1", "replicate_1"), summarise, cor = round(cor(reads_1, reads_2, use="complete.obs"), 2))
cors$locus_1 <- "Nmel1132"
p <- ggplot(edna_replicate_freqs_long, aes(reads_1, reads_2, col=locus_1)) +
  geom_point() + ylim(0,1) + xlim(0,1) + xlab("eDNA sample 1") + ylab("eDNA sample 2") +
  geom_abline(slope=1, intercept=0, colour="gray",linetype=2) +
  scale_color_manual(values=viridis(28)) + theme_bw()
FigS2 <- p + facet_grid(vars(mesocosm_1), vars(replicate_1)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.2, y=0.9)
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS2_duplicate_read_freqs_corr.eps"), plot=FigS2, dpi = 300)   

# Figures S4-S5: variation in read depth among mesocosms
per_locus_sums <- rbind(aggregate(edna_replicate_counts_long$reads_1, by=list(edna_replicate_counts_long$locus_1,edna_replicate_counts_long$sample_1), FUN=sum),
                        aggregate(edna_replicate_counts_long$reads_2, by=list(edna_replicate_counts_long$locus_2,edna_replicate_counts_long$sample_2), FUN=sum))
per_locus_sums <- as.data.frame(cSplit(per_locus_sums, "Group.2", '_', drop=FALSE))
colnames(per_locus_sums) <- c("locus","sample","counts","treatment","replicate")
per_locus_sums$treatment <- factor(per_locus_sums$treatment, levels=c("eDNA1","eDNA3","eDNA5","eDNA10"))


total_count_sums <- rbind(aggregate(edna_replicate_counts_long$reads_1, by=list(edna_replicate_counts_long$sample_1), FUN=sum),
                        aggregate(edna_replicate_counts_long$reads_2, by=list(edna_replicate_counts_long$sample_2), FUN=sum))
total_count_sums <- as.data.frame(cSplit(total_count_sums, "Group.1", '_', drop=FALSE))
colnames(total_count_sums) <- c("sample","counts","treatment","replicate")
total_count_sums$treatment <- factor(total_count_sums$treatment, levels=c("eDNA1","eDNA3","eDNA5","eDNA10"))


# Figure S4: average read depth per eDNA sample
FigS4A <- ggplot(data=per_locus_sums, aes(x=treatment, y=counts, fill=treatment)) +
  geom_boxplot(show.legend=F) + geom_jitter(shape=16, position=position_jitter(0.2),colour="darkgray",show.legend=F) +
  xlab("Mesocosm density treatment") + ylab("Total reads per locus") +
  scale_fill_manual(values=rev(viridis(4))) + theme_bw() + 
  scale_x_discrete(labels=c("eDNA1"="1 fish","eDNA3"="3 fish","eDNA5"="5 fish","eDNA10"="10 fish")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS4A_per_locus_reads.eps"), plot=FigS4A, dpi = 300)   

FigS4B <- ggplot(data=total_count_sums, aes(x=treatment, y=counts, fill=treatment)) +
  geom_boxplot(show.legend=F) + geom_jitter(shape=16, position=position_jitter(0.2),colour="darkgray",show.legend=F) +
  xlab("Mesocosm density treatment") + ylab("Total reads per sample") +
  scale_fill_manual(values=rev(viridis(4))) + theme_bw() + 
  scale_x_discrete(labels=c("eDNA1"="1 fish","eDNA3"="3 fish","eDNA5"="5 fish","eDNA10"="10 fish")) +
  ylim(0,80000) + theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS4B_total_reads.eps"), plot=FigS4B, dpi = 300)   

