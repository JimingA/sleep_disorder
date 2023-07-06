# --------------------------------------------------------
# | Project: Sleep Disturbance Scale for Children GWAS   |
# | sarang seo                                           |
# | purpose: make phenotype
# | 
# --------------------------------------------------------

### set env
#setRepositories(ind=1:8)
setwd('/disk6/bissr/GWAS')

### load library
library(data.table)
library(tidyverse)

### load data
# fam
fam <- fread('0.data/genomics_sample03/ABCD_release_3.0_QCed.fam')
colnames(fam) <-c('FID','IID','PID','MID','Sex','P')

# sdsc
sdsc <- fread('0.data/abcd_sds01.txt')
sdsc.info <- sdsc[1, ]
sdsc <- sdsc[-1,]
colnames(sdsc)


### 1. sdsc cleansing -> make total
sdsc <- sdsc %>% 
  mutate_at(vars(c(8,9)),as.factor) %>%   # sex, eventname, 26 questionnaire answers
  mutate_at(vars(c(11:36)),as.numeric)

# 1) baseline
baseline_total <- sdsc %>% 
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  rowwise() %>% 
  mutate(total = sum(c_across(starts_with("sleepdisturb")), na.rm = T)) %>% 
  select(src_subject_id, eventname, total) %>% 
  rename(IID = src_subject_id)

sum(baseline_total$total<39) 

# 2) 1 year
oneyear_total <- sdsc %>%
  filter(eventname == '1_year_follow_up_y_arm_1') %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("sleepdisturb")), na.rm = T)) %>%
  select(src_subject_id, eventname, total) %>%
  rename(IID = src_subject_id)

sum(oneyear_total$total<39)


# 3) 2 year
twoyear_total <- sdsc %>%
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("sleepdisturb")), na.rm = T)) %>%
  select(src_subject_id, eventname, total) %>%
  rename(IID = src_subject_id)

sum(twoyear_total$total<39)


### 2. Phenotype file 
famID <- fam %>% filter(IID %in% baseline_total$IID) %>% select(FID, IID)


# 1) baseline
phe_sdsc <- merge(famID, baseline_total, all=T) %>% 
  select(-eventname) %>% na.omit() %>% 
  rename(TOTAL_SCORE=total)
phe_sdsc <- phe_sdsc[,c(2,1,3)]        

#write.table(phe_sdsc,file = '3.phenotype/phe_sdsc.txt',sep = '\t',col.names = T, row.names = F, quote = F)


# 2) 1 year
phe_sdsc_oneyear <- merge(famID, oneyear_total, all=T) %>% 
  select(-eventname) %>% na.omit() %>% 
  rename(TOTAL_SCORE=total)
phe_sdsc_oneyear <- phe_sdsc_oneyear[,c(2,1,3)]        

#write.table(phe_sdsc,file = '3.phenotype/phe_sdsc_oneyear.txt',sep = '\t',col.names = T, row.names = F, quote = F)


# 3) 2 year
phe_sdsc_twoyear <- merge(famID, oneyear_total, all=T) %>% 
  select(-eventname) %>% na.omit() %>% 
  rename(TOTAL_SCORE=total)
phe_sdsc_twoyear <- phe_sdsc_twoyear[,c(2,1,3)]        

#write.table(phe_sdsc,file = '3.phenotype/phe_sdsc_twoyear.txt',sep = '\t',col.names = T, row.names = F, quote = F)


### 3.
sdsc_total <- sdsc %>% 
  mutate(total_score = select(.,c(11:36)) %>% rowSums()) %>% 
  select(src_subject_id, eventname, total_score)
sdsc_disease <- sdsc_total %>% 
  mutate(disease = case_when(
    total_score >= 39 ~ "Sleep disturbed",
    TRUE ~ "Healthy"))

#write.table(sdsc_disease, file = '3.phenotype/0.phe_eventname_disease.txt',sep = '\t',col.names = T, row.names = F, quote = F)
