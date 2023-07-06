# --------------------------------------------------------
# | Project: Sleep Disturbance Scale for Children GWAS   |
# | sarang seo                                           |
# | purpose: 
# | 
# --------------------------------------------------------

### set env
setRepositories(ind=1:7)
setwd('/disk6/bissr/GWAS')

### load library
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

### load data
phe <- fread('3.phenotype/phe_sdsc.txt')
#phe_hy <- fread('/disk6/bissr/SDSC_gwas/3.Phenotype/phe_sdsc.txt')

white <- fread('2.filterlist/out_white_samplelist.txt', header=F)
trans <- fread('2.filterlist/out_trans_samplelist.txt')
trans_clear <- fread('2.filterlist/out_trans_samplelist_clear.txt')

pc_white <- fread('4.covariate/white/pca_white.eigenvec') %>% select(V2,V3,V4) %>% rename(IID=V2, PC1=V3, PC2=V4)
pc_trans <- fread('4.covariate/trans/pca_trans.eigenvec') %>% select(V2,V3,V4) %>% rename(IID=V2, PC1=V3, PC2=V4)
pc_trans_clear <- fread('4.covariate/trans_clear/pca_trans.eigenvec') %>% select(V2,V3,V4) %>% rename(IID=V2, PC1=V3, PC2=V4)


### 1. make data
# 1) white
phe_white <- phe %>% 
  filter(IID %in% white$V2) %>% 
  mutate(DISEASE = case_when(
    TOTAL_SCORE >= 39 ~ "Disease",
    TRUE ~ "Healthy",
  )) 

# 2 ) trans
phe_trans <- phe %>% 
  filter(IID %in% trans$IID) %>% 
  mutate(DISEASE = case_when(
    TOTAL_SCORE >= 39 ~ "Disease",
    TRUE ~ "Healthy",
  )) 


# 3) trnas - clear
phe_trans_clear <- phe %>% 
  filter(IID %in% trans_clear$IID) %>% 
  mutate(DISEASE = case_when(
    TOTAL_SCORE >= 39 ~ "Disease",
    TRUE ~ "Healthy",
  )) 

table(phe_white$DISEASE)
table(phe_trans$DISEASE)
table(phe_trans_clear$DISEASE)

#write.table(phe_white,file = '3.phenotype/phe_disease_white.txt',sep = '\t',col.names = T, row.names = F, quote = F)
#write.table(phe_trans,file = '3.phenotype/phe_disease_trans.txt',sep = '\t',col.names = T, row.names = F, quote = F)
#write.table(phe_trans_clear,file = '3.phenotype/phe_disease_trans_clear.txt',sep = '\t',col.names = T, row.names = F, quote = F)


### 2. plot -pca
plot_white <- merge(phe_white, pc_white, by='IID') %>% 
  mutate(DISEASE = factor(DISEASE))

ggplot(plot_white, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = DISEASE),alpha=0.4)+
  scale_color_brewer(palette = 'YlGnBu')+
  theme_classic()

plot_trans <- merge(phe_trans, pc_trans, by='IID') %>% 
  mutate(DISEASE = factor(DISEASE))

ggplot(plot_trans, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = DISEASE),alpha=0.4)+
  scale_color_brewer(palette = 'YlGnBu')+
  theme_classic()


