# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in synchronous population dynamics 
# UKBMS climate synchrony analysis

rm(list=ls()) # clear R

## load packages
library(lme4)
library(lmerTest)
library(broom)
library(MuMIn)
library(ggplot2)
library(parallel)
library(doParallel)

# read in data
pair_attr <- readRDS("pop_climate_synchrony_UKBMS_final.rds")  # pop synchrony for 32 species and climate synchrony for 8 variables

## make sure correct variables are factors ## 
str(pair_attr)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$pair.id <- as.character(pair_attr$pair.id)

length(unique(pair_attr$spp)) ## 32 species

# Run climate model with all 8 climate variables and covariates
climate_model1 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                         scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                         scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
summary(climate_model1) ## summer temperature non-significant, remaining climate variables are significant
# save model output 
results_table_climate_ukbms <- data.frame(summary(climate_model1)$coefficients[,1:5]) ## 31 species
write.csv(results_table_climate_ukbms, file = "../Results/UKBMS/climate_synchrony.csv", row.names=TRUE)

# Run second model with covariates and summer temperature only to calculate R2 of significant climate variables
climate_model2 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(summer_temp) +
                         (1|pair.id) + (1|spp), data = pair_attr)
summary(climate_model2)
## calc R2 and take difference (marginal == fixed effects)
r.squaredGLMM(climate_model1)
r.squaredGLMM(climate_model2)
x = 0.01260022 - 0.004692347   
x2 = x*100 ## 0.791% or 0.00791 proportion





###############################################################

# Next, run permutation tests to determine true significance of significant climate variables

# 1. Spring temperature
true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values from main climate model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with spring temperature
true_result_table <- true_result_table[grep("spring_temp", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_spring_temp <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$spring_temp_shuffle <- sample(pair_attr$spring_temp) ## randomly shuffle spring temperature varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                  scale(spring_temp_shuffle) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
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
ukbms_perm_spring_temp$parameter <- paste(row.names(ukbms_perm_spring_temp)) ## move row.names to parameter column
rownames(ukbms_perm_spring_temp) <- 1:nrow(ukbms_perm_spring_temp) ## change row names to numbers
ukbms_perm_spring_temp <- ukbms_perm_spring_temp[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
ukbms_perm_spring_temp <- ukbms_perm_spring_temp[grep("spring_temp_shuffle", ukbms_perm_spring_temp$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_spring_temp) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_spring_temp_ukbms.csv", row.names=TRUE)
## read in file
perm_spring_temp_ukbms <- read.csv("../Results/UKBMS/perm_spring_temp_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_spring_temp_ukbms <- perm_spring_temp_ukbms[perm_spring_temp_ukbms$i==0,]
perm_spring_temp_ukbms <- perm_spring_temp_ukbms[!perm_spring_temp_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_spring_temp_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_spring_temp_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant




# 2. Autumn temperature

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with autumn temperature
true_result_table <- true_result_table[grep("autumn_temp", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_autumn_temp <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
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
ukbms_perm_autumn_temp$parameter <- paste(row.names(ukbms_perm_autumn_temp)) ## move row.names to parameter column
rownames(ukbms_perm_autumn_temp) <- 1:nrow(ukbms_perm_autumn_temp) ## change row names to numbers
ukbms_perm_autumn_temp <- ukbms_perm_autumn_temp[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
ukbms_perm_autumn_temp <- ukbms_perm_autumn_temp[grep("autumn_temp_shuffle", ukbms_perm_autumn_temp$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_autumn_temp) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_autumn_temp_ukbms.csv", row.names=TRUE)
## read in file
perm_autumn_temp_ukbms <- read.csv("../Results/UKBMS/perm_autumn_temp_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_autumn_temp_ukbms <- perm_autumn_temp_ukbms[perm_autumn_temp_ukbms$i==0,]
perm_autumn_temp_ukbms <- perm_autumn_temp_ukbms[!perm_autumn_temp_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_autumn_temp_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_autumn_temp_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant


# 3. Winter temperature

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with winter temperature
true_result_table <- true_result_table[grep("winter_temp", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_winter_temp <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$winter_temp_shuffle <- sample(pair_attr$winter_temp) ## randomly shuffle winter temperature varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp_shuffle) + scale(autumn_temp) + 
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
ukbms_perm_winter_temp$parameter <- paste(row.names(ukbms_perm_winter_temp)) ## move row.names to parameter column
rownames(ukbms_perm_winter_temp) <- 1:nrow(ukbms_perm_winter_temp) ## change row names to numbers
ukbms_perm_winter_temp <- ukbms_perm_winter_temp[,-c(1:3)] ## remove unnecessary columns
## only interested in shuffle main effect
ukbms_perm_winter_temp <- ukbms_perm_winter_temp[grep("winter_temp_shuffle", ukbms_perm_winter_temp$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_winter_temp) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_winter_temp_ukbms.csv", row.names=TRUE)
## read in file
perm_winter_temp_ukbms <- read.csv("../Results/UKBMS/perm_winter_temp_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_winter_temp_ukbms <- perm_winter_temp_ukbms[perm_winter_temp_ukbms$i==0,]
perm_winter_temp_ukbms <- perm_winter_temp_ukbms[!perm_winter_temp_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_winter_temp_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_winter_temp_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant


# 4. Spring rainfall

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with spring rainfall
true_result_table <- true_result_table[grep("spring_rain", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_spring_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
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
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_perm_spring_rain$parameter <- paste(row.names(ukbms_perm_spring_rain)) ## move row.names to parameter column
rownames(ukbms_perm_spring_rain) <- 1:nrow(ukbms_perm_spring_rain) ## change row names to numbers
ukbms_perm_spring_rain <- ukbms_perm_spring_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in shuffle main effect
ukbms_perm_spring_rain <- ukbms_perm_spring_rain[grep("spring_rain_shuffle", ukbms_perm_spring_rain$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_spring_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_spring_rain_ukbms.csv", row.names=TRUE)
## read in file
perm_spring_rain_ukbms <- read.csv("../Results/UKBMS/perm_spring_rain_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_spring_rain_ukbms <- perm_spring_rain_ukbms[perm_spring_rain_ukbms$i==0,]
perm_spring_rain_ukbms <- perm_spring_rain_ukbms[!perm_spring_rain_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_spring_rain_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_spring_rain_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant



# 5. Summer rainfall

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with summer rainfall
true_result_table <- true_result_table[grep("summer_rain", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_summer_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$summer_rain_shuffle <- sample(pair_attr$summer_rain) ## randomly shuffle summer rainfall varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain_shuffle) + scale(winter_temp) + scale(autumn_temp) + 
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
ukbms_perm_summer_rain$parameter <- paste(row.names(ukbms_perm_summer_rain)) ## move row.names to parameter column
rownames(ukbms_perm_summer_rain) <- 1:nrow(ukbms_perm_summer_rain) ## change row names to numbers
ukbms_perm_summer_rain <- ukbms_perm_summer_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in shuffle main effect
ukbms_perm_summer_rain <- ukbms_perm_summer_rain[grep("summer_rain_shuffle", ukbms_perm_summer_rain$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_summer_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_summer_rain_ukbms.csv", row.names=TRUE)
## read in file
perm_summer_rain_ukbms <- read.csv("../Results/UKBMS/perm_summer_rain_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_summer_rain_ukbms <- perm_summer_rain_ukbms[perm_summer_rain_ukbms$i==0,]
perm_summer_rain_ukbms <- perm_summer_rain_ukbms[!perm_summer_rain_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_summer_rain_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_summer_rain_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant



# 6. Autumn rainfall

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with autumn rainfall
true_result_table <- true_result_table[grep("autumn_rain", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_autumn_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$autumn_rain_shuffle <- sample(pair_attr$autumn_rain) ## randomly shuffle autumn rainfall varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain_shuffle) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
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
ukbms_perm_autumn_rain$parameter <- paste(row.names(ukbms_perm_autumn_rain)) ## move row.names to parameter column
rownames(ukbms_perm_autumn_rain) <- 1:nrow(ukbms_perm_autumn_rain) ## change row names to numbers
ukbms_perm_autumn_rain <- ukbms_perm_autumn_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in shuffle main effect
ukbms_perm_autumn_rain <- ukbms_perm_autumn_rain[grep("autumn_rain_shuffle", ukbms_perm_autumn_rain$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_autumn_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_autumn_rain_ukbms.csv", row.names=TRUE)
## read in file
perm_autumn_rain_ukbms <- read.csv("../Results/UKBMS/perm_autumn_rain_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_autumn_rain_ukbms <- perm_autumn_rain_ukbms[perm_autumn_rain_ukbms$i==0,]
perm_autumn_rain_ukbms <- perm_autumn_rain_ukbms[!perm_autumn_rain_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_autumn_rain_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_autumn_rain_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant



# 7. Winter rainfall

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
ukbms_perm_winter_rain <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
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
ukbms_perm_winter_rain$parameter <- paste(row.names(ukbms_perm_winter_rain)) ## move row.names to parameter column
rownames(ukbms_perm_winter_rain) <- 1:nrow(ukbms_perm_winter_rain) ## change row names to numbers
ukbms_perm_winter_rain <- ukbms_perm_winter_rain[,-c(1:3)] ## remove unnecessary columns
## only interested in shuffle main effect
ukbms_perm_winter_rain <- ukbms_perm_winter_rain[grep("winter_rain_shuffle", ukbms_perm_winter_rain$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_winter_rain) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_winter_rain_ukbms.csv", row.names=TRUE)
## read in file
perm_winter_rain_ukbms <- read.csv("../Results/UKBMS/perm_winter_rain_ukbms.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_winter_rain_ukbms <- perm_winter_rain_ukbms[perm_winter_rain_ukbms$i==0,]
perm_winter_rain_ukbms <- perm_winter_rain_ukbms[!perm_winter_rain_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_winter_rain_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_winter_rain_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant


## Results summary ## 
# Spring temperature: p = 0
# Autumn temperature: p = 0
# Winter temperature: p = 0
# Spring rainfall: p = 0
# Summer rainfall: p = 0
# Autumn rainfall: p = 0.006
# Winter rainfall: p = 0

