# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in synchronous population dynamics 
# CBC climate synchrony analysis

rm(list=ls()) # clear R

## load packages
library(lme4)
library(lmerTest)
library(broom)
library(MuMIn)
library(ggplot2)
library(doSNOW)

# read in data
pair_attr <- readRDS("pop_climate_synchrony_CBC_final.rds") 

## make sure correct variables are factors ## 
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
str(pair_attr)

length(unique(pair_attr$spp)) ## 26 species

# Run climate model with all 8 climate variables and covariates
climate_model1 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(winter_rain) + 
                         scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                         scale(spring_temp) + scale(summer_temp) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
summary(climate_model1) ## summer temperature is significant
# save model output 
results_table_climate_cbc <- data.frame(summary(climate_model1)$coefficients[,1:5]) ## 31 species
write.csv(results_table_climate_cbc, file = "../Results/CBC/climate_synchrony.csv", row.names=TRUE)

## run model without summer temperature to calculate R2 
climate_model2 <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(winter_rain) + 
                         scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                         scale(spring_temp) + (1|pair.id) + (1|spp), data = pair_attr)
## calc R2 and take difference (marginal == fixed effects)
r.squaredGLMM(climate_model1)
r.squaredGLMM(climate_model2)
x = 0.00093626 - 0.000799783
x2 = x*100 ## 0.0136% or 0.000136


###############################################################

# Next, run permutation tests to determine true significance of significant climate variables

# 1. Summer temperature 

true_result_table <- data.frame(anova(climate_model1)[5]) ## save F values
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with spring temperature
true_result_table <- true_result_table[grep("summer_temp", true_result_table$parameter),]

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
cbc_perm_summer_temp <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4', .options.snow=opts) %dopar% {
  print(i)
  pair_attr$summer_temp_shuffle <- sample(pair_attr$summer_temp) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ scale(mid.year) + scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + 
                  scale(spring_temp) + scale(summer_temp_shuffle) + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 13.24 minutes for 999 runs

### save results
cbc_perm_summer_temp$parameter <- paste(row.names(cbc_perm_summer_temp)) ## move row.names to parameter column
rownames(cbc_perm_summer_temp) <- 1:nrow(cbc_perm_summer_temp) ## change row names to numbers
cbc_perm_summer_temp <- cbc_perm_summer_temp[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
cbc_perm_summer_temp <- cbc_perm_summer_temp[grep("summer_temp_shuffle", cbc_perm_summer_temp$parameter),]
final_results_table <- rbind(true_result_table, cbc_perm_summer_temp) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
## plot distribution of F values with vertical line (true F value)
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/CBC/perm_summer_temp_cbc.csv", row.names=TRUE)
## read in file
perm_summer_temp_cbc <- read.csv("../Results/CBC/perm_summer_temp_cbc.csv", header=TRUE)

## Calculate p value
number_of_permutations <- 1000
true_summer_temp_cbc <- perm_summer_temp_cbc[perm_summer_temp_cbc$i==0,]
perm_summer_temp_cbc <- perm_summer_temp_cbc[!perm_summer_temp_cbc$i==0,] ## remove true value to calc. p value
diff.observed <- true_summer_temp_cbc$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_summer_temp_cbc$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.028 significant

# Summary: summer temperature is still significant (p=0.028)



