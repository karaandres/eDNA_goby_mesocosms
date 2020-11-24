# eDNA_goby_mesocosms

## This repository contains the files related to the bioinformatics, data processing, and data analysis for the eDNA round goby mesocosm experiment. 

### Overview

In this project, we evaluated intraspecific genetic diversity in the round goby (Neogobius melanostomus) using eDNA samples from experimental mesocosms and a field trial. We tested the similarity between eDNA and individual tissue-based estimates of allele frequencies and used a likelihood-based DNA mixture framework to estimate the number of unique genetic contributors in mesocosm eDNA samples and in simulated mixtures of alleles. We also tested the similarity between eDNA and tissue-based estimates of allele frequencies in a field trial.

**Sequencing:**
  - Illumina MiSeq 2x250 bp kit
  - Raw read output: 47920390 reads
  
**Processing**
  - Mesocosm experiment data processinng: [mult_amp_mesocosms.txt](mult_amp_mesocosms.txt)
  - Field trial data processing: [mult_amp_field_trial.txt](mult_amp_field_trial.txt)

**Analysis:**
  - Mesocosm experiment analysis: [mesocosm_analysis.R](mesocosm_analysis.R)
  - Field trial analysis: [Cayuga_field_trial.R](Cayuga_field_trial.R)
