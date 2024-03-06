# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in synchronous population dynamics 
# BBS climate synchrony analysis

rm(list=ls()) # clear R

## load packages
library(lme4)
library(lmerTest)
library(broom)
library(MuMIn)
library(ggplot2)
library(doSNOW)

# read in data
pair_attr <- readRDS(file = "pop_climate_synchrony_BBS_final.rds") # read in pair_attr file 

## make sure correct variables are factors ## 
str(pair_attr)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$pair.id <- as.character(pair_attr$pair.id)

length(unique(pair_attr$spp)) ## 24 species

# Run climate model with all 8 climate variables and covariates
climate_model1 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                         scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                         scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
summary(climate_model1) ## winter rain, spring rain and autumn temperature are significant
# save model output 
results_table_climate_bbs <- data.frame(summary(climate_model1)$coefficients[,1:5]) ## 31 species
write.csv(results_table_climate_bbs, file = "../Results/BBS/climate_synchrony.csv", row.names=TRUE)

## run model without winter rain, spring rain and autumn temperature to calculate R2 of significant variables
climate_model2 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + 
                         scale(autumn_rain) + scale(summer_rain) + scale(winter_temp) + 
                         scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
r.squaredGLMM(climate_model1)
r.squaredGLMM(climate_model2)
x = 0.0001192723 - 0.0000957853 
x ## 0.0000235
x2 = x*100 ## 0.00235%


###############################################################

# Next, run permutation tests to determine true significance of significant climate variables

# 1. Spring rainfall
true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with spring rainfall
true_result_table <- true_result_table[grep("spring_rain", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)

cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=8000, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

n_sims <- 999
para_start_time = Sys.time()
bbs_perm_spring_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4', .options.snow=opts) %dopar% {
  print(i)
  pair_attr$spring_rain_shuffle <- sample(pair_attr$spring_rain) ## randomly shuffle spring rainfall varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain_shuffle) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                  scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_perm_spring_rain$parameter <- paste(row.names(bbs_perm_spring_rain)) ## move row.names to parameter column
rownames(bbs_perm_spring_rain) <- 1:nrow(bbs_perm_spring_rain) ## change row names to numbers
bbs_perm_spring_rain <- bbs_perm_spring_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
bbs_perm_spring_rain <- bbs_perm_spring_rain[grep("spring_rain_shuffle", bbs_perm_spring_rain$parameter),]
final_results_table <- rbind(true_result_table, bbs_perm_spring_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/BBS/perm_spring_rain_bbs.csv", row.names=TRUE)
## read in file
perm_spring_rain_bbs <- read.csv("../Results/BBS/perm_spring_rain_bbs.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_spring_rain_bbs <- perm_spring_rain_bbs[perm_spring_rain_bbs$i==0,]
perm_spring_rain_bbs <- perm_spring_rain_bbs[!perm_spring_rain_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_spring_rain_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_spring_rain_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.003 significant


# 2. Winter rainfall

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with winter rainfall
true_result_table <- true_result_table[grep("winter_rain", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_perm_winter_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$winter_rain_shuffle <- sample(pair_attr$winter_rain) ## randomly shuffle winter rainfall varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain_shuffle) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                  scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_perm_winter_rain$parameter <- paste(row.names(bbs_perm_winter_rain)) ## move row.names to parameter column
rownames(bbs_perm_winter_rain) <- 1:nrow(bbs_perm_winter_rain) ## change row names to numbers
bbs_perm_winter_rain <- bbs_perm_winter_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
bbs_perm_winter_rain <- bbs_perm_winter_rain[grep("winter_rain_shuffle", bbs_perm_winter_rain$parameter),]
final_results_table <- rbind(true_result_table, bbs_perm_winter_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()
## save file
write.csv(final_results_table, file = "../Results/BBS/perm_winter_rain_bbs.csv", row.names=TRUE)
## read in file
perm_winter_rain_bbs <- read.csv("../Results/BBS/perm_winter_rain_bbs.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_winter_rain_bbs <- perm_winter_rain_bbs[perm_winter_rain_bbs$i==0,]
perm_winter_rain_bbs <- perm_winter_rain_bbs[!perm_winter_rain_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_winter_rain_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_winter_rain_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.001 significant


# 3. Autumn temperature

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with autumn temperature
true_result_table <- true_result_table[grep("autumn_temp", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_perm_autumn_temp <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$autumn_temp_shuffle <- sample(pair_attr$autumn_temp) ## randomly shuffle autumn temperature varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp_shuffle) + 
                  scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_perm_autumn_temp$parameter <- paste(row.names(bbs_perm_autumn_temp)) ## move row.names to parameter column
rownames(bbs_perm_autumn_temp) <- 1:nrow(bbs_perm_autumn_temp) ## change row names to numbers
bbs_perm_autumn_temp <- bbs_perm_autumn_temp[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
bbs_perm_autumn_temp <- bbs_perm_autumn_temp[grep("autumn_temp_shuffle", bbs_perm_autumn_temp$parameter),]
final_results_table <- rbind(true_result_table, bbs_perm_autumn_temp) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/BBS/perm_autumn_temp_bbs.csv", row.names=TRUE)
## read in file
perm_autumn_temp_bbs <- read.csv("../Results/BBS/perm_autumn_temp_bbs.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_autumn_temp_bbs <- perm_autumn_temp_bbs[perm_autumn_temp_bbs$i==0,]
perm_autumn_temp_bbs <- perm_autumn_temp_bbs[!perm_autumn_temp_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_autumn_temp_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_autumn_temp_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.002 significant


## Results summary: all climate variables are still significant 
# Spring rainfall: 0.003
# Winter rainfall: 0.001
# Autumn temperature: 0.002
