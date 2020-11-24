### Cayuga eDNA/tissue comparison
### Analysis added for Molecular Ecology revisions, round 2
### Last updated 9.15.2020 by Kara Andres

# This script accomplishes the following: 
# Calculate tissue allele frequencies
# Calculate eDNA allele frequencies
# Exclude loci with poor coverage or not in HWE
# Double check that the Cayuga gobies are panmictic 
# Correlation between population allele/read frequencies 

library(adegenet)
library(pegas)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(splitstackshape)
library(vegan)
library(data.table)
library(devtools)
library(ecodist)
library(RColorBrewer)
library(stringr)
library(viridis)

#############################################################################################
####################### Calculate tissue allele frequencies ########################
#############################################################################################

# Recode genotypes for missing data (< 10 reads/locus)
# Import read count and genotype matrices for both field experiment and mesocosm experiment
# read in MultAmp output: col = samples, row = loci, cells = allele1/allele2:reads1,reads2
hap_genotype_field <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/hap_genotype.csv", header=TRUE)
hap_genotype_field <- data.frame(lapply(hap_genotype_field[,-2], gsub, pattern=".*:", replacement=""))   # Remove haplotype column, remove everything before ":" in each cell
hap_genotype_field <- hap_genotype_field[,-grep("e", colnames(hap_genotype_field))] # remove eDNA samples
hap_genotype_mat <- t(as.matrix(hap_genotype_field[,-1]))   # Transpose, turn into matrix

# Re-code genotypes for missing data 
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
genotypes_field <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/hap_genotype_matrix.csv", header=TRUE)
genotypes_field <- genotypes_field[-grep("e_", genotypes_field$X),] # remove eDNA samples
genotypes <- as.matrix(genotypes_field[,-1]) # remove row and col names

# Loop through to re-code alleles appropriately
#If total count <10, locus gets "0|0"
genotypes_edited <-  matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat))
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  if (is.na(counts_less10[i]==TRUE)) {
    genotypes_edited[i] <- paste0(0, "|", 0)} else {
      genotypes_edited[i] = genotypes[i]
    }
}

genotypes_edited_df <- as.data.frame(genotypes_edited) # Turn genotypes_edited into df, name columns (loci), and split by "|"
colnames(genotypes_edited_df) <- hap_genotype_field$Locus # name columns (loci)
genotypes_edited_df_split <- as.data.frame(cSplit(genotypes_edited_df, 1:ncol(genotypes_edited_df), '|'))  # split each locus by "|"
# write.csv(genotypes_edited_df_split, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/total_population_genotypes_less10_split.csv")

# Calculate number of individuals and alleles with missing data
apply(genotypes_edited_df_split, 1, function(x) {sum(x==0)}) # number of alleles with missing data per individual (row)
apply(genotypes_edited_df_split, 2, function(x) {sum(x==0)}) # number of alleles with missing data per locus (column)

allelic_richness <- data.frame(allele = NA, allelic_rich = NA)
for (i in seq(1,ncol(genotypes_edited_df_split)-1,2)){
  uniq_alleles <- c(genotypes_edited_df_split[,i], genotypes_edited_df_split[,i+1])
  allele <- colnames(genotypes_edited_df_split)[i]
  allelic_rich <- length(unique(uniq_alleles[uniq_alleles>0]))
  allelic_richness_temp <- cbind(allele, allelic_rich)
  allelic_richness <- rbind(allelic_richness, allelic_richness_temp) 
}
allelic_richness <- allelic_richness[-1,]
allelic_richness$allelic_rich <- as.numeric(allelic_richness$allelic_rich)
mean(allelic_richness$allelic_rich, na.rm = TRUE) # mean alleic rich = 9.4

# Global genotype frequencies
genotypes_edited_df_total <- genotypes_edited_df_split
locus_positions <- seq(1, ncol(genotypes_edited_df_split)-1, 2)  # starting column number for each locus
lnames <- colnames(genotypes_edited_df_split)
OUT <- NULL
for (j in locus_positions) { # for each locus (31)
  alleles <- c(genotypes_edited_df_total[,j], genotypes_edited_df_total[,j+1]) # Combine 2 columns per locus
  alleles2 <- as.data.frame(table(alleles)) # count each allele at locus x
  alleles3 <- alleles2[alleles2$Freq!=0 & alleles2$alleles!=0,] # remove missing data (otherwise 0 would be counted in total number of alleles)
  alleles4 <- cbind(alleles3,alleles3$Freq/sum(alleles3$Freq)) # calculate frequencies
  if (length(alleles4) > 0) {
    output <- cbind(rep(lnames[j], nrow(alleles4)),alleles4) #combine j, locus name, and frequencies
    OUT <<- rbind(OUT,output)          
  }
}

colnames(OUT) <- c("locus","allele","count","inds_freq") #add column headers
total_tissue_allele_freqs <- OUT
total_tissue_allele_freqs$locus <- gsub(pattern = "_1", replacement = "", x = total_tissue_allele_freqs$locus) # remove subscript from locus names 
# write.csv(total_tissue_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/total_tissue_allele_freqs.csv")

#############################################################################################
########################## Calculate eDNA read frequencies ##########################
#############################################################################################

# Total normalized read frequencies for all counts > 10 and frequencies > 0.01
# Import read counts for all alleles at all loci
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus

# Begin loop
# create empty dataframe to populate with allele frequencies for each mesocosm
total_edna_allele_freqs <- data.frame(NULL)
j = 1 # individual locus indexed from "filenames"
# Create matrix with total number of reads/haplotype <=1 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,-1]) # remove "haplotype" column, turn into matrix
  file_mat <- file_mat[,grep("e_", colnames(file_mat))] # subset to eDNA samples
  file_mat <- file_mat[,-grep("_BL", colnames(file_mat))] # remove blank sample
  file_mat[file_mat <= 1] <- 0 # replace counts <1 with 0
  file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
  file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
  file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
  edna_allele_freqs_temp <- data.frame(locus=rep(filenames[[j]]), allele=file[,1], count=rowSums(file_mat_norm, na.rm = TRUE)) # read count normalized per sample
  edna_allele_freqs_temp$CAY_eDNA <- edna_allele_freqs_temp$count/sum(edna_allele_freqs_temp$count, na.rm = TRUE)
  
  edna_allele_freqs_temp[is.na(edna_allele_freqs_temp)] <- 0
  total_edna_allele_freqs <- rbind(total_edna_allele_freqs, edna_allele_freqs_temp)
  j <- j+1
}  

# write.csv(total_edna_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/total_edna_allele_freqs.csv")

# Check blanks
# create empty dataframe to populate with max reads for each controls
edna_control_reads <- data.frame(locus=names(dataFiles), BL=NA)
# Create matrix that removes haplotypes w/ <10 reads
x <- 1
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,grep("BL", colnames(file))])
  file_mat[file_mat <= 1] <- 0 # replace counts <=1 with 0
  edna_control_reads[x,2] <- sum(file_mat)
  x <- x+1
}

# Combine eDNA and tisuse allele freqs
total_allele_freqs <- merge(total_edna_allele_freqs, total_tissue_allele_freqs, by = c("locus","allele"), all=TRUE)
total_allele_freqs[is.na(total_allele_freqs)] <- 0
total_allele_freqs <-  total_allele_freqs[total_allele_freqs$CAY_eDNA>0|total_allele_freqs$inds_freq>0,]
# write.csv(total_allele_freqs, "/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/total_allele_freqs.csv")

# Read frequencies per eDNA sample
sample_edna_allele_freqs <- data.frame(NULL)
j = 1 # individual locus indexed from "filenames"
# Create matrix with total number of reads/haplotype <=1 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,-1]) # remove "haplotype" column, turn into matrix
  file_mat <- file_mat[,grep("e_", colnames(file_mat))] # subset to eDNA samples
  file_mat <- file_mat[,-grep("_BL", colnames(file_mat))] # remove blank sample
  file_mat[file_mat <= 1] <- 0 # replace counts <10 with 0
  file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) if (sum(x)>0){x/sum(x, na.rm = TRUE)} else {x})) # read count normalized per sample
  file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
  file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) if (sum(x)>0){x/sum(x, na.rm = TRUE)} else {x})) # overwrite with alleles < threshold removed
  edna_allele_freqs_temp <- data.frame(locus=rep(filenames[[j]]), allele=file[,1], file_mat_norm) # read count normalized per sample
  edna_allele_freqs_temp <- edna_allele_freqs_temp[rowSums(edna_allele_freqs_temp[3:5])>0,]
  sample_edna_allele_freqs <- rbind(sample_edna_allele_freqs, edna_allele_freqs_temp)
  j <- j+1
}  

#############################################################################################
#################### Check for HWE and identify problem loci  ##########################
#############################################################################################
geno <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/hap_genotype_matrix.csv")
geno <- geno[-grep("e_", geno$X),] # remove eDNA samples 
rownames(geno) <- geno$X # rownames are sample names
geno <- geno[,!names(geno) %in% "X"] # remove column of sample names
geno[geno == ".|."] <- "NA|NA" # replace missing data ".|." with NAs
geno.obj <- df2genind(geno, sep="\\|", ploidy=2, type="codom") # "|" must be preceded by double backslashes
summary(geno.obj)
hw.output <- data.frame(hw.test(geno.obj, 10000))
hw.output$observed_heterozygosity <- summary(geno.obj)[[6]]
hw.output$expected_heterozygosity <- summary(geno.obj)[[7]]
hw.output[hw.output$Pr.exact< 0.05/35 & hw.output$observed>hw.output$expected,] # bonferroni corrected pvalue
# remove loci 155, 248, 361, and 726 -- deviation from HWE and heterzygote excess
# also remove l1103 and 1531 because of 1.00 heterozygosity and 185 due to poor coverage

# Remove loci from datasets: 
patterns <- c("Nmel185","Nmel155","Nmel248","Nmel361","Nmel726", "Nmel1103","Nmel1531") # loci to remove
total_tissue_allele_freqs_HWE <- filter(total_tissue_allele_freqs, !grepl(paste(patterns, collapse="|"), locus))
total_edna_allele_freqs_HWE <- filter(total_edna_allele_freqs, !grepl(paste(patterns, collapse="|"), locus))
total_allele_freqs_HWE <- filter(total_allele_freqs, !grepl(paste(patterns, collapse="|"), locus))
sample_edna_allele_freqs_HWE <- filter(sample_edna_allele_freqs, !grepl(paste(patterns, collapse="|"), locus))

#############################################################################################
################## Double check that the Cayuga gobies are panmictic ####################
#############################################################################################
geno <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/hap_genotype_matrix.csv")
geno <- geno[, -which(names(geno) %in% patterns)]
geno <- geno[-grep("e_", geno$X),] # remove eDNA samples 
rownames(geno) <- geno$X
geno <- geno[,-1]
geno[geno == ".|."] <- "NA|NA"
geno[geno == "|"] <- "/"
geno.obj <- df2genind(geno, sep="/", ploidy=2, type="codom")
geno.cent <- scaleGen(geno.obj, center= T, scale= F, NA.method="mean")
geno.clusters <- find.clusters(geno.obj)
# retain 50 PC axes (~90% of variance) for a max of 7 clusters
# one cluster has the lowest BIC 

#############################################################################################
############## Correlation between population allele/read frequencies ###############
#############################################################################################
cor.test(total_allele_freqs_HWE$inds_freq, total_allele_freqs_HWE$CAY_eDNA)
cors <- ddply(total_allele_freqs_HWE, "locus", summarise, cor=round(cor(inds_freq,CAY_eDNA,use="complete.obs"), 2))
min(cors$cor)
max(cors$cor)

# png('/Users/kbja10/Downloads/Cayuga_eDNA_ind_corr.png',width = 6, height = 6, units = "in", res = 400)
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(total_allele_freqs_HWE$inds_freq, total_allele_freqs_HWE$CAY_eDNA, xlim = c(0,1), ylim = c(0,1),
     xlab = "Tissue-based allele frequencies", ylab = "eDNA allele frequencies", cex.lab = 1.5, cex = 1.5, lwd = 1.5,
     col = "darkblue")
# dev.off()

allele_recovery <- total_allele_freqs_HWE[total_allele_freqs_HWE$inds_freq>0,]
nrow(allele_recovery) # 253 alleles in genotyped individuals
nrow(allele_recovery[allele_recovery$CAY_eDNA>0,]) # 121/253 alleles recovered
edna_missed_alleles <- allele_recovery[allele_recovery$CAY_eDNA==0,]
max(edna_missed_alleles$inds_freq) # all alleles above 0.24 were detected in eDNA
mean(edna_missed_alleles$inds_freq) # mean undetected allele freq is 0.03
sd(edna_missed_alleles$inds_freq)
new_alleles <- total_allele_freqs_HWE[total_allele_freqs_HWE$CAY_eDNA>0&total_allele_freqs_HWE$inds_freq==0,]
nrow(new_alleles)
mean(new_alleles$CAY_eDNA)
sd(new_alleles$CAY_eDNA)
max(new_alleles$CAY_eDNA)
min(new_alleles$CAY_eDNA)

ggplot(total_allele_freqs_HWE, aes(x=inds_freq, y=CAY_eDNA, color=locus)) +
  xlim(0,1.1)+ylim(0,1.1) + geom_point() + geom_text(aes(label=locus),hjust=0, vjust=0)
p <- ggplot(total_allele_freqs_HWE, aes(x=inds_freq, y=CAY_eDNA, color=locus)) +
  geom_point() + ylab("Field eDNA allele frequency") + xlab("Genotyped tissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(28)) + theme_bw()
# + geom_text(aes(label=locus),hjust=0, vjust=0)
ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/Fig1_field_allele_freq_corr2.eps"), plot=p, dpi = 300,  width = 7, height = 5, units = "in")   

#############################################################################################
############################# Contributor estimator  #######################################
#############################################################################################
# contributor estimator 
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

pop_allele_freqs_edna <- unite(total_edna_allele_freqs_HWE, "locus_allele", c(locus, allele), sep = "_")
totalfreqs <- pop_allele_freqs_edna[,c("locus_allele","CAY_eDNA")]

pop_allele_freqs <- unite(total_tissue_allele_freqs_HWE, "locus_allele", c(locus, allele), sep = "_")
totalfreqs <- pop_allele_freqs[,c("locus_allele","inds_freq")]

sample_freqs <- unite(sample_edna_allele_freqs_HWE, "locus_allele", c(locus, allele), sep = "_")
sample_freqs <- sample_freqs[,c("locus_allele","e_CAY_03")]
sample_freqs <- sample_freqs[sample_freqs$e_CAY_03>0,]
  
# Combine eDNA and population allele frequencies
combined_allele_freqs <- merge(sample_freqs, totalfreqs, by = "locus_allele")
# combined_allele_freqs$inds_freq[is.na(combined_allele_freqs$inds_freq)] <- 0.01
# combined_allele_freqs <- combined_allele_freqs[combined_allele_freqs$inds_freq>0,]
combined_allele_freqs_split <- cSplit(combined_allele_freqs, "locus_allele", '_', drop = FALSE)
colnames(combined_allele_freqs_split)[c(4,5)] <- c("locus", "allele")

# Combine eDNA and population allele frequencies
combined_allele_freqs <- unite(total_allele_freqs_HWE, "locus_allele", c(locus, allele), sep = "_", remove=FALSE)
colnames(combined_allele_freqs)[c(2,3)] <- c("locus", "allele")
combined_allele_freqs_split <- combined_allele_freqs[combined_allele_freqs$CAY_eDNA>0&combined_allele_freqs$inds_freq>0,]

combined_allele_freqs <- totalfreqs[totalfreqs$CAY_eDNA>0.01,]
combined_allele_freqs_split <- cSplit(combined_allele_freqs, "locus_allele", '_', drop = FALSE)
colnames(combined_allele_freqs_split)[c(3,4)] <- c("locus", "allele")

# estimator 
loci.list <- split(combined_allele_freqs_split, combined_allele_freqs_split$locus)
loci.list <- lapply(loci.list, function(x) {x$inds_freq})
MixtureLikelihood2(y = 100, loci = loci.list)$Contributors

loci.list <- split(combined_allele_freqs_split, combined_allele_freqs_split$locus)
loci.list <- lapply(loci.list, function(x) {x$CAY_eDNA})
MixtureLikelihood2(y = 100, loci = loci.list)$Contributors

##############################################################################################
############################# Read depth per sample, Fig S3  #######################################
#############################################################################################

dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/cayuga_edna_tissue_comparison/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus
# Begin loop
j = 1 # individual locus indexed from "filenames"
total_read_depth <- data.frame(NULL)
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,-grep("Haplotype", colnames(file))])
  file_mat[file_mat < 10] <- 0 # replace counts <10 with 0
  file_mat_sum <- data.frame(apply(file_mat, 2, function(x) sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
  total_read_depth_temp <- data.frame(sample=rownames(file_mat_sum), locus=rep(filenames[[j]], nrow(file_mat_sum)), 
                                      total_reads=apply(file_mat, 2, function(x) sum(x, na.rm = TRUE)),
                                      mean_reads=apply(file_mat, 2, function(x) mean(x[x>0], na.rm = TRUE)))
  total_read_depth <- rbind(total_read_depth, total_read_depth_temp)
  j <- j+1
}

# add mesocosm eDNA samples 
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus
# Begin loop
j = 1 # individual locus indexed from "filenames"
total_read_depth2 <- data.frame(NULL)
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,grep("e", colnames(file))])
  file_mat <- file_mat[,-grep("Haplotype", colnames(file_mat))]
  file_mat[file_mat < 10] <- 0 # replace counts <10 with 0
  file_mat_sum <- data.frame(apply(file_mat, 2, function(x) sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
  total_read_depth2_temp <- data.frame(sample=rownames(file_mat_sum), locus=rep(filenames[[j]], nrow(file_mat_sum)), 
                                      total_reads=apply(file_mat, 2, function(x) sum(x, na.rm = TRUE)),
                                      mean_reads=apply(file_mat, 2, function(x) mean(x[x>0], na.rm = TRUE)))
  total_read_depth2 <- rbind(total_read_depth2, total_read_depth2_temp)
  j <- j+1
}

total_read_depth <- rbind(total_read_depth, total_read_depth2)
total_read_depth$sample_type <- ifelse(grepl("t", total_read_depth$sample, ignore.case = T), "tissue",
                                ifelse(grepl("e_", total_read_depth$sample, ignore.case = T), "field_eDNA","mesocosm_eDNA"))
patterns <- c("Nmel185","Nmel155","Nmel248","Nmel361","Nmel726", "Nmel1103","Nmel1531") # loci to remove
total_read_depth <- total_read_depth[!grepl(paste(patterns, collapse="|"), total_read_depth$locus),]
total_read_depth <- total_read_depth[-grep("BL", total_read_depth$sample),]
aggregate(total_read_depth$total, by=list(Category=total_read_depth$sample_type), FUN=function(x) c(mean = mean(x), sd = sd(x)))

sum_read_depth <- aggregate(total_read_depth$total_reads, by=list(Category=total_read_depth$sample), FUN=sum)
sum_read_depth$sample <- ifelse(grepl("t", sum_read_depth$Category, ignore.case = T), "tissue",
                        ifelse(grepl("e_", sum_read_depth$Category, ignore.case = T), "field_eDNA","mesocosm_eDNA"))
aggregate(sum_read_depth$x, by=list(Category=sum_read_depth$sample), FUN=function(x) c(mean = mean(x), sd = sd(x)))

mean_read_depth <- aggregate(total_read_depth$mean_reads, by=list(Category=total_read_depth$sample), FUN=mean, na.rm=T)
mean_read_depth$sample <- ifelse(grepl("t", mean_read_depth$Category, ignore.case = T), "tissue",
                                ifelse(grepl("e_", mean_read_depth$Category, ignore.case = T), "field_eDNA","mesocosm_eDNA"))

p <- ggplot(sum_read_depth, aes(x=sample, y=x)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2), colour="darkgray") + 
  scale_y_log10(limits = c(1000,100000)) + theme_classic() + 
  ylab("Total read depth") + xlab("Sample type")
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS3_totalreadspersample.eps"), plot=p, dpi = 300,  width = 5, height = 5, units = "in")   

ggplot(mean_read_depth, aes(x=sample, y=x)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic()

p <- ggplot(total_read_depth[total_read_depth$total_reads>0,], aes(x=sample_type, y=total_reads)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2), colour="darkgray") + 
  scale_y_continuous(trans='log10') + theme_classic() + 
    ylab("Total reads per locus") + xlab("Sample type")
#ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Mesocosm_experiment_2018/mesocosm_analysis/FigS3_totalreadsperlocus.eps"), plot=p, dpi = 300,  width = 5, height = 5, units = "in")   
