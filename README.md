# eDNA_goby_mesocosms

## This repository contains the files related to the bioinformatics, data processing, and data analysis for the eDNA round goby mesocosm experiment. 

### Overview

In this project, we evaluated intraspecific genetic diversity in the round goby (Neogobius melanostomus) using eDNA samples from experimental mesocosms. We tested the similarity between eDNA and individual tissue-based estimates of allele frequencies and used a likelihood-based DNA mixture framework to estimate the number of unique genetic contributors in mesocosm eDNA samples and in simulated mixtures of alleles. We also tested the similarity between eDNA and tissue-based estimates of allele frequencies in a field trial.

**Sequencing:**
  - Illumina MiSeq 2x250 bp kit
  - Raw read output: 

**Analysis:**
  - Calculate tissue allele frequencies: [tissue_allele_freqs.R](tissue_allele_freqs.R)
  - Calculate eDNA allele frequencies: [edna_allele_freqs.R](edna_allele_freqs.R)
  - Identify problem loci and remove from datasets: [inspect_loci.R](inspect_loci.R)
  - PCA, heatmap, and correlation of tissue allele frequencies and eDNA read frequencies: [PCA_heatmap_corr.R](PCA_heatmap_corr.R)
  - Contributor estimation: [contributor_estimation.R](contributor_estimation.R)
  - Supplementary figures: [supp_figs.R](supp_figs.R)
  - Field trial analysis: [Cayuga_field_trial.R](Cayuga_field_trial.R)
