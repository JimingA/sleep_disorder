# --------------------------------------------------------
# | Project: Sleep Disturbance Scale for Children GWAS   |
# | sarang seo                                           |
# | purpose: make covariate file
# | save 4.covariate/cov_~~
# |                  dum_~
# --- -----------------------------------------------------

### set env
#setRepositories(ind=1:8)
setwd('/disk6/bissr/GWAS')


### load library
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(genio)
library(stringr)
library(janitor)
library(reshape2)
library(fastDummies)

### load data
# lnstitute
lt <- fread('0.data/abcd_lt01.txt')
lt.info <- lt[1, ]
lt <- lt[-1, ]

# cov demo
cov_demo <- fread('4.covariate/cov_demo.txt')
cov_demo.c <- fread('4.covariate/cov_demo_clear.txt')

# samplelist
list.w <- fread('2.filterlist/out_white_samplelist.txt', header=F) %>% rename(FID=V1, IID=V2)
list.t <- fread('2.filterlist/out_trans_samplelist.txt')

# eigenval eigenvalue
val.w <- fread('4.covariate/white/pca_white.eigenval')
vec.w <- fread('4.covariate/white/pca_white.eigenvec')

val.t <- fread('4.covariate/trans/pca_trans.eigenval')
vec.t <- fread('4.covariate/trans/pca_trans.eigenvec')

val.c <- fread('4.covariate/trans_clear/pca_trans.eigenval')
vec.c <- fread('4.covariate/trans_clear/pca_trans.eigenvec')


### 1. scree plot
# # white 
# var_explained <- val.w/sum(val.w)   #번째 pc가 전체 데이터의 분산 중 몇 퍼센트 설명한다.
# sum(var_explained[1:3])
# 
# qplot(c(1:10),var_explained$V1[1:10])+
#   geom_line() + 
#   xlab("Principal Component") + 
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot -white") +
#   ylim(0, 0.0015)+
#   theme_classic()
# 
# # trans 
# var_explained <- val.t/sum(val.t)   
# sum(var_explained[1:3])
# 
# qplot(c(1:10),var_explained$V1[1:10])+
#   geom_line() + 
#   xlab("Principal Component") + 
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot -trans") +
#   ylim(0, 0.03)+
#   theme_classic()
# 
# 
# # trans -clear 
# var_explained <- val.c/sum(val.c)   
# sum(var_explained[1:3])
# 
# qplot(c(1:10),var_explained$V1[1:10])+
#   geom_line() + 
#   xlab("Principal Component") + 
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot -trans_clear") +
#   ylim(0, 0.03)+
#   theme_classic()


### 2. data processing

# 1) cov : IID AGE SEX SITE

cov <- lt %>% 
  filter(eventname=='baseline_year_1_arm_1') %>% 
  select(src_subject_id, interview_age,sex, site_id_l) %>% 
  rename(IID = src_subject_id, SEX = sex, AGE = interview_age, SITE = site_id_l) %>% 
  mutate(IID = as.factor(IID), 
         SEX = sex_to_int(SEX), 
         AGE=as.numeric(AGE),
         SITE = as.numeric(str_remove(SITE, 'site')))


#list.w %>% filter(!(IID %in% cov$IID))  
#list.t %>% filter(!(IID %in% cov$IID)) 
#table(cov$SITE)


# 2) white : FID IID AGE SEX SITE -> dummy -> save
covariate.w <- merge(list.w, cov, all = T) %>% na.omit() %>% relocate(FID) 
  #dummy_cols(select_columns = 'SITE') %>% 
  #select(-SITE)

#write.table(covariate.w, file = '4.covariate/tmp_white.txt',sep = '\t',col.names = T, row.names = F, quote = F)

# 3)
dum.w <-  fread('4.covariate/z/cov_white.cov')
#write.table(merge(covariate.w %>% select(-SITE), dum.w), file = '4.covariate/cov_white.txt',sep = '\t',col.names = T, row.names = F, quote = F)


### 3. trans ) PC or RACE / clear ver

# 1) trans : FID IID PC1-3 AGE SEX SITE
cov_pc.t<- vec.t %>% rename(FID=V1, IID=V2, PC1=V3, PC2=V4, PC3=V5) %>% 
  select(FID, IID, PC1,PC2,PC3) 

covariate_pc <- merge(cov_pc.t, cov, all = T) %>% na.omit() %>% relocate(FID) 
#write.table(covariate_pc, file = '4.covariate/z/tmp_trans_pc.txt',sep = '\t',col.names = T, row.names = F, quote = F)

dum.t.pc <-  fread('4.covariate/z/cov_trans_pc.cov')
#write.table(merge(covariate_pc %>% select(-SITE), dum.t.pc), file = '4.covariate/cov_trans_pc.txt',sep = '\t',col.names = T, row.names = F, quote = F)


# 2) trans : FID IID AGE SEX SITE RACE 
covariate_demo <- merge(cov, cov_demo, all = T) %>% na.omit() %>% relocate(FID) 
#write.table(covariate_demo, file = '4.covariate/z/tmp_trans_race.txt',sep = '\t',col.names = T, row.names = F, quote = F)

dum.r.pc <-  fread('4.covariate/z/cov_trans_race.cov')
#write.table(merge(covariate_demo %>% select(-SITE,-RACE), dum.r.pc), file = '4.covariate/cov_trans_race.txt',sep = '\t',col.names = T, row.names = F, quote = F)


# 3) trans : FID IID PC1-3 AGE SEX SITE -clear
cov_pc.c<- vec.c %>% rename(FID=V1, IID=V2, PC1=V3, PC2=V4, PC3=V5) %>% 
  select(FID, IID, PC1,PC2,PC3) 

covariate_pc.c <- merge(cov_pc.c, cov, all = T) %>% na.omit() %>% relocate(FID) 
#write.table(covariate_pc.c, file = '4.covariate/z/tmp_trans_pc_clear.txt',sep = '\t',col.names = T, row.names = F, quote = F)

dum.tc.pc <-  fread('4.covariate/z/cov_trans_pc_clear.cov')
#write.table(merge(covariate_pc.c %>% select(-SITE), dum.tc.pc), file = '4.covariate/cov_trans_pc_clear.txt',sep = '\t',col.names = T, row.names = F, quote = F)


# 4) trans : FID IID AGE SEX SITE RACE -clear 
covariate_demo.c <- merge(cov, cov_demo.c, all = T) %>% na.omit() %>% relocate(FID)  
#write.table(covariate_demo.c, file = '4.covariate/z/tmp_trans_race_clear.txt',sep = '\t',col.names = T, row.names = F, quote = F)

dum.rc.pc <-  fread('4.covariate/z/cov_trans_race_clear.cov')
#write.table(merge(covariate_demo.c %>% select(-SITE,-RACE), dum.rc.pc), file = '4.covariate/cov_trans_race_clear.txt',sep = '\t',col.names = T, row.names = F, quote = F)



# # ================================================================================
# # 2) cov_pc : FID IID pc1-3
# cov_pc.w <- eigenvec.w %>% rename(FID=V1, IID=V2, PC1=V3, PC2=V4, PC3=V5) %>% 
#   select(FID, IID, PC1,PC2,PC3)
# cov_pc.t<- eigenvec.t %>% rename(FID=V1, IID=V2, PC1=V3, PC2=V4, PC3=V5) %>% 
#   select(FID, IID, PC1,PC2,PC3)
# 
# # 3) covarite  : FID IID PC1-3 AGE SEX SITE
# covariate.w <- merge(cov_pc.w, cov, all = T) %>% na.omit() %>% relocate(FID) 
# covariate.t <- merge(cov_pc.t, cov, all = T) %>% na.omit() %>% relocate(FID)
# 
# # 4) covarite_demo  : FID IID PC1-3 AGE SEX SITE RACE
# covariate_demo <- merge(covariate.t, cov_demo, all = T) %>% na.omit()
# 
# # +
# site_cnt <- covariate.w %>% group_by(SITE) %>% summarise(n=n())
