# --------------------------------------------------------
# | Project: Sleep Disturbance Scale for Children GWAS   |
# | sarang seo                                           |
# | purpose: 
# | 
# --------------------------------------------------------

### set env
setRepositories(ind=1:8)
setwd('/disk6/bissr/GWAS')

### load library
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

### load data

# race name
raceName <- as.matrix(fread('0.data/race.txt',header = F))

# fam
fam <- fread('0.data/genomics_sample03/ABCD_release_3.0_QCed.fam')
colnames(fam) <-c('FID','IID','PID','MID','Sex','P')

# sdsc
sdsc <- fread('0.data/abcd_sds01.txt')
sdsc.info <- sdsc[1, ]
sdsc <- sdsc[-1,]

# demo
demo <- fread('0.data/pdem02.txt')
demo.info <- demo[1, ]
demo <- demo[-1,]

# 2) demo -> long -> remove duplication
race_col <- which(str_detect(colnames(pdem),pattern='demo_race_a_p__'))
demo <- pdem %>% select(src_subject_id, race_col) 
colnames(demo) <- raceName

demo_long <- demo %>% 
  mutate_at(vars(-IID), as.numeric) %>% 
  gather(c(2:ncol(demo)), key="race", value = "race_val") %>% 
  filter(race_val == 1) %>%
  mutate(race = as.factor(race))

dup <- demo_long %>% 
  mutate(IID = as.factor(IID)) %>% 
  group_by(IID) %>% 
  summarise(sum = sum(race_val))%>% 
  filter(sum >= 2)  #1479

demo_long2 <- demo_long %>% filter (!(IID %in% dup$IID)) # 10375


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






### 2. modeling

# 1) data select & scale
data_ml <- data1 %>% 
  mutate_at(vars(-disease), scale) %>% 
  mutate_at(vars(-disease), as.numeric)

# 2) Data shuffling
randomIdx <- sample(1:nrow(data_ml))
data_ml <- data_ml[randomIdx,] 

# 3) Make train & test data
div <- createFolds(data_ml$disease, 10)
#div <- createDataPartition(y = data_ml$disease, p = 0.7, list = F)

train <- data_ml[-div[["Fold01"]], ]
test <-  data_ml[div[["Fold01"]], ]

# 4) modeling
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

ctrl <- trainControl(method = "cv", number = 10, savePredictions = "all", summaryFunction=twoClassSummary, classProbs=T )

model.svm <- train(disease~., data = train, method = "svmLinear", trainControl=ctrl)
model.rf <- train(disease~., data = train, method = "rf", trainControl=ctrl)
model.lb <- train(disease~., data = train, method = "LogitBoost", trainControl=ctrl)
pam
pda

stopCluster(cl)

# +)
# ctrl <- trainControl(method = "repeatedcv",number=10, repeats = 3)
# model.knn <- train(host_disease~., data=train_phylum,
#                                  method="knn", trControl=ctrl)





### 3. res compare
#
model.svm$resample
model.rf$resample

model.svm$results


## Prediction
pred.svm <- predict(model.svm, test)
pred.rf  <- predict(model.rf, test)

confu.svm <- confusionMatrix(data = pred.svm, reference = test$disease, positive = "Sleep disturbed")
confu.rf  <- confusionMatrix(data = pred.rf, reference = test$disease, positive = "Sleep disturbed")
confusionMatrix(model.svm, "none")


confu.svm$table
