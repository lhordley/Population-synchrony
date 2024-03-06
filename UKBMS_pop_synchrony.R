# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in synchronous population dynamics 
# UKBMS dataset analysis

rm(list=ls()) # clear R

library(plyr) # load packages
library(dplyr) 
library(lme4)
library(ggplot2)
library(lmerTest)
library(plotrix)
library(ggeffects)
library(DHARMa)
library(foreach)
library(doParallel)
library(stringr)
options(scipen=999)

## load data
pair_attr <- readRDS("pop_climate_synchrony_UKBMS_final.rds") # butterfly pair attribute data
str(pair_attr)
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$abund_change_00_12 <- as.factor(pair_attr$abund_change_00_12)
pair_attr$abund_change_85_00 <- as.factor(pair_attr$abund_change_85_00)
pair_attr$sig_abund_change_00_12 <- as.factor(pair_attr$sig_abund_change_00_12)
pair_attr$sig_abund_change_85_00 <- as.factor(pair_attr$sig_abund_change_85_00)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$specialism <- as.factor(pair_attr$specialism)

summary(pair_attr)

##########################################
### 1a. Average synchrony + specialism ### 
##########################################

length(unique(pair_attr$spp)) # 32 spp

spec_model_ukbms1 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(mid.year) + specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr)
summary(spec_model_ukbms1)
# non-significant

results_table_average_spec <- data.frame(summary(spec_model_ukbms1)$coefficients[,1:5]) ## 32 species
confidence_intervals <- data.frame(confint(spec_model_ukbms1, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_average_spec <- cbind(results_table_average_spec,confidence_intervals)

write.csv(results_table_average_spec, file = "../Results/UKBMS/average_spec_ukbms.csv", row.names=TRUE) # 32 species

############################################
### 1b. Change in synchrony + specialism ### 
############################################

pair_attr_early <- droplevels(pair_attr[pair_attr$mid.year==1984.5 | pair_attr$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr[pair_attr$mid.year==1999.5 | pair_attr$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

#### Early (1985-2000) ####
spec_model_ukbms2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(mid.year)*specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr_early)
summary(spec_model_ukbms2)
# significant
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms, file = "../Results/UKBMS/change_spec_ukbms_85_00.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(spec_model_ukbms2, terms = c("mid.year", "specialism"))
plot(dat) # generalists decline in synchrony more steeply compared to specialists

# Create figure 2a
# run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table_spec<-NULL
for (i in unique(pair_attr_early$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                              scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                              scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) +(1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                         scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                         scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                       scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                       scale(winter_temp) + scale(autumn_temp) + scale(spring_temp), data=pair_attr_early[pair_attr_early$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table_spec<-rbind(results_table_spec,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
specialism <- unique(pair_attr_early[,c("spp", "specialism")])
results_table_spec <- merge(results_table_spec, specialism, by.x="i", by.y="spp", all.x=TRUE)
## remove parameter column
results_table_spec <- results_table_spec[,-c(4)]
## re-name some columns
names(results_table_spec) <- c("species", "slope", "SE", "Specialism")
results_table_spec$Specialism <- revalue(results_table_spec$Specialism, c("Habitat specialist"="Specialist"))
results_table_spec$Specialism <- revalue(results_table_spec$Specialism, c("Wider countryside sp"="Generalist"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_spec_ukbms <- results_table_spec_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_spec_ukbms$Specialism <- paste(row.names(results_table_spec_ukbms))
rownames(results_table_spec_ukbms) <- 1:nrow(results_table_spec_ukbms)
results_table_spec_ukbms <- results_table_spec_ukbms[c(12,14),] # keep mid year and interaction rows
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("scale(mid.year)"="Specialist"))
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("scale(mid.year):specialismWider countryside sp"="Generalist"))

## create new dataframe
specialist_slope <- results_table_spec_ukbms[1,1]
generalist_slope <- sum(results_table_spec_ukbms$Estimate)
specialist_SE <- results_table_spec_ukbms[1,2]
interaction_SE <- results_table_spec_ukbms[2,2]
generalist_SE <- sqrt((specialist_SE)^2) + ((interaction_SE)^2)
model_summary_spec <- data.frame(Specialism=c("Specialist", "Generalist"), slope=c(specialist_slope, generalist_slope), SE=c(specialist_SE, generalist_SE))

model_summary_spec$Specialism <- factor(model_summary_spec$Specialism, levels=c("Generalist", "Specialist"))
results_table_spec$Specialism <- factor(results_table_spec$Specialism, levels=c("Generalist", "Specialist"))

png("../Graphs/New/Specialism_change_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
spec <- ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table_spec, size=2, color="grey71", position = myjit) +
  geom_errorbar(data=results_table_spec, color="grey71", position = myjit, width=0.1) +
  geom_point(data=model_summary_spec, size=4) +
  geom_errorbar(data=model_summary_spec, width=0.1) +
  labs(x="Specialism", y="Change in population \n synchrony 1985-2000") +
  #scale_y_continuous(breaks = seq(-1.8,0.8,0.4)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
spec
## generalists have a more negative change in syncrhony than specialists 
dev.off()

#### run permutations ####
## save true interaction results (to merge in with permutations later)
main_result_table <- data.frame(anova(spec_model_ukbms2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):specialism"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_early$spec_shuffle <- sample(pair_attr_early$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                  scale(autumn_temp) + scale(spring_temp) + scale(mid.year)*spec_shuffle + 
                  (1|pair.id) + (1|spp), data = pair_attr_early)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 31 minutes for 999 runs

### save results
ukbms_spec_para$parameter <- paste(row.names(ukbms_spec_para)) ## move row.names to parameter column
rownames(ukbms_spec_para) <- 1:nrow(ukbms_spec_para) ## change row names to numbers
ukbms_spec_para <- ukbms_spec_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spec_para <- ukbms_spec_para[grep("mid.year):spec_shuffle", ukbms_spec_para$parameter),]

final_results_table <- rbind(main_result_table, ukbms_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_spec_ukbms_85_00.csv", row.names=FALSE)
## read in file
perm_change_spec_ukbms <- read.csv("../Results/UKBMS/perm_change_spec_ukbms_85_00.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_spec_ukbms <- perm_change_spec_ukbms[perm_change_spec_ukbms$i==0,]
perm_change_spec_ukbms <- perm_change_spec_ukbms[!perm_change_spec_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant



#### Late (2000-2012) ####
spec_model_ukbms3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(mid.year)*specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr_late)
summary(spec_model_ukbms3)
# significant 
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms3)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms, file = "../Results/UKBMS/change_spec_ukbms_00_12.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(spec_model_ukbms3, terms = c("mid.year", "specialism"))
plot(dat) # generalists increase in synchrony more steeply compared to specialists

# run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table_spec<-NULL
for (i in unique(pair_attr_late$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                              scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                              scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) +(1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                         scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                         scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                       scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                       scale(winter_temp) + scale(autumn_temp) + scale(spring_temp), data=pair_attr_late[pair_attr_late$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table_spec<-rbind(results_table_spec,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
specialism <- unique(pair_attr_late[,c("spp", "specialism")])
results_table_spec <- merge(results_table_spec, specialism, by.x="i", by.y="spp", all.x=TRUE)
## remove parameter column
results_table_spec <- results_table_spec[,-c(4)]
## re-name some columns
names(results_table_spec) <- c("species", "slope", "SE", "Specialism")
results_table_spec$Specialism <- revalue(results_table_spec$Specialism, c("Habitat specialist"="Specialist"))
results_table_spec$Specialism <- revalue(results_table_spec$Specialism, c("Wider countryside sp"="Generalist"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_spec_ukbms <- results_table_spec_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_spec_ukbms$Specialism <- paste(row.names(results_table_spec_ukbms))
rownames(results_table_spec_ukbms) <- 1:nrow(results_table_spec_ukbms)
results_table_spec_ukbms <- results_table_spec_ukbms[c(12,14),] # keep mid year and interaction rows
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("scale(mid.year)"="Specialist"))
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("scale(mid.year):specialismWider countryside sp"="Generalist"))

## create new dataframe
specialist_slope <- results_table_spec_ukbms[1,1]
generalist_slope <- sum(results_table_spec_ukbms$Estimate)
specialist_SE <- results_table_spec_ukbms[1,2]
interaction_SE <- results_table_spec_ukbms[2,2]
generalist_SE <- sqrt((specialist_SE)^2) + ((interaction_SE)^2)
model_summary_spec <- data.frame(Specialism=c("Specialist", "Generalist"), slope=c(specialist_slope, generalist_slope), SE=c(specialist_SE, generalist_SE))

model_summary_spec$Specialism <- factor(model_summary_spec$Specialism, levels=c("Generalist", "Specialist"))
results_table_spec$Specialism <- factor(results_table_spec$Specialism, levels=c("Generalist", "Specialist"))

png("../Graphs/New/Specialism_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
spec2 <- ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table_spec, size=2, color="grey71", position = myjit) +
  geom_errorbar(data=results_table_spec, color="grey71", position = myjit, width=0.1) +
  geom_point(data=model_summary_spec, size=4) +
  geom_errorbar(data=model_summary_spec, width=0.1) +
  labs(x="Specialism", y="Change in population \n synchrony 2000-2012") +
  #scale_y_continuous(breaks = seq(-1.8,0.8,0.4)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
spec2
dev.off()

#### run permutations ####
## save true interaction results (to merge in with permutations later)
main_result_table <- data.frame(anova(spec_model_ukbms3)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):specialism"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$spec_shuffle <- sample(pair_attr_late$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                  scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                  scale(autumn_temp) + scale(spring_temp) + scale(mid.year)*spec_shuffle + 
                  (1|pair.id) + (1|spp), data = pair_attr_late)
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
ukbms_spec_para$parameter <- paste(row.names(ukbms_spec_para)) ## move row.names to parameter column
rownames(ukbms_spec_para) <- 1:nrow(ukbms_spec_para) ## change row names to numbers
ukbms_spec_para <- ukbms_spec_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spec_para <- ukbms_spec_para[grep("mid.year):spec_shuffle", ukbms_spec_para$parameter),]
final_results_table <- rbind(main_result_table, ukbms_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_spec_ukbms_00_12.csv", row.names=FALSE)
## read in file
perm_change_spec_ukbms <- read.csv("../Results/UKBMS/perm_change_spec_ukbms_00_12.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_spec_ukbms <- perm_change_spec_ukbms[perm_change_spec_ukbms$i==0,]
perm_change_spec_ukbms <- perm_change_spec_ukbms[!perm_change_spec_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant




########################################
### 2a. Average synchrony + mobility ### 
########################################

str(pair_attr)
pair_attr_mob <- pair_attr[!is.na(pair_attr$mobility_wil),]
length(unique(pair_attr_mob$spp)) # 31 spp

mobility_model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + scale(autumn_rain) + 
                         scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + 
                         scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_mob)
summary(mobility_model)
## significant (positive)
results_table_mob_ukbms <- data.frame(summary(mobility_model)$coefficients[,1:5]) ## 31 species
confidence_intervals <- data.frame(confint(mobility_model, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_mob_ukbms <- cbind(results_table_mob_ukbms,confidence_intervals)

write.csv(results_table_mob_ukbms, file = "../Results/UKBMS/average_mob_ukbms.csv", row.names=TRUE)

#### run permutations ####
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility main effect
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mobility_wil)"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para1 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_mob$mob_shuffle <- sample(pair_attr_mob$mobility_wil) ## randomly shuffle mobility varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + scale(autumn_rain) + 
                  scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + 
                  scale(mob_shuffle) + (1|spp) + (1|pair.id), data=pair_attr_mob)
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
ukbms_mob_para1$parameter <- paste(row.names(ukbms_mob_para1)) ## move row.names to parameter column
rownames(ukbms_mob_para1) <- 1:nrow(ukbms_mob_para1) ## change row names to numbers
ukbms_mob_para1 <- ukbms_mob_para1[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
ukbms_mob_para1 <- ukbms_mob_para1[grep("mob_shuffle", ukbms_mob_para1$parameter),]

final_results_table <- rbind(main_result_table, ukbms_mob_para1) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_average_mob_ukbms.csv", row.names=TRUE)
## read in file
perm_average_mob_ukbms <- read.csv("../Results/UKBMS/perm_average_mob_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_average_mob_ukbms <- perm_average_mob_ukbms[perm_average_mob_ukbms$i==0,]
perm_average_mob_ukbms <- perm_average_mob_ukbms[!perm_average_mob_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_average_mob_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_average_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.004 significant


##########################################
### 2b. Change in synchrony + mobility ### 
##########################################

pair_attr_early <- droplevels(pair_attr_mob[pair_attr_mob$mid.year==1984.5 | pair_attr_mob$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr_mob[pair_attr_mob$mid.year==1999.5 | pair_attr_mob$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

#### Early 1985-2000 ####
mobility_model2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(mid.year)*scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(mobility_model2)
## interaction non-significant
results_table_mobility_ukbms <- data.frame(summary(mobility_model2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/UKBMS/change_mobility_ukbms_85_00.csv", row.names=TRUE)

#### Late 2000-2012 ####
mobility_model3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(mid.year)*scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(mobility_model3)
## interaction significant
results_table_mobility_ukbms <- data.frame(summary(mobility_model3)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/UKBMS/change_mobility_ukbms_00_12.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(mobility_model3, terms = c("mid.year", "mobility_wil"))
plot(dat) # more mobile species increase in synchrony more rapidly compared to less mobile species

## run model without mobility for each species and extract the slope (mid.year) coefficient for each
results_table_mob<-NULL
for (i in unique(pair_attr_late$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({mob_ukbms <- lmer(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                           + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                           + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    mob_ukbms <- lmer(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                      + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                      + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    mob_ukbms <- lm(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                    + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                    + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(mob_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table_mob<-rbind(results_table_mob,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
mobility <- pair_attr[,c("spp", "mobility_wil")]
mobility <- unique(mobility)
results_table_mob <- merge(results_table_mob, mobility, by.x="i", by.y="spp")
## remove parameter column
results_table_mob <- results_table_mob[,-c(4)]
## re-name some columns
names(results_table_mob) <- c("species", "slope", "SE", "Mobility")


########## calculate slope and SE from main model (with mobility as x axis)
### plot points and SE ontop of raw data
results_table_mobility_ukbms <- results_table_mobility_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_mobility_ukbms$Mobility <- paste(row.names(results_table_mobility_ukbms))
rownames(results_table_mobility_ukbms) <- 1:nrow(results_table_mobility_ukbms)
results_table_mobility_ukbms <- results_table_mobility_ukbms[c(12,14),]

## extract slope and SE for interaction and mid.year
mid.year_slope <- results_table_mobility_ukbms[1,1]
interaction_slope <- results_table_mobility_ukbms[2,1]
se_mid.year <- results_table_mobility_ukbms[1,2]
se_interaction <- results_table_mobility_ukbms[2,2]

## get slopes for each value of mobility
mobility<-unique(results_table_mob$Mobility)
slopes <- mid.year_slope + interaction_slope * scale(mobility)
slopes 

## get standard errors for each value of mobility
SEs <- rep(NA, length(mobility))
for (i in 1:length(mobility)){
  j <- mobility[i]  
  SEs[i] <- sqrt((se_interaction * j)^2 + se_mid.year^2) 
}
model_summary_mob <- data.frame(Mobility=mobility, slope=slopes, SE=SEs)

png("../Graphs/New/FigureX_butterflies_mob.png", height = 150, width = 180, units = "mm", res = 300)
mob <- ggplot(mapping=aes(x=Mobility, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table_mob, size=2, color="grey71") +
  geom_errorbar(data=results_table_mob, color="grey71", width=2) +
  geom_ribbon(data=model_summary_mob, fill = "grey70", alpha=0.4) +
  geom_line(data=model_summary_mob, size=1) +
  labs(x="Mobility score", y="Change in population \n synchrony 2000-2012") +
  #scale_x_continuous(breaks = seq(0,4.5,0.5)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
mob
dev.off()

# plot both graphs together for figure 2
library(ggpubr)
png("../Graphs/New/Figure2_sync.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(spec, mob, 
          labels = c("(a)", "(b)"), font.label = list(size = 16, color ="black"),
          ncol = 2, nrow = 1)
dev.off()

#### run permutations ####

## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model3)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):scale(mobility_wil)"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para2 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$mob_shuffle <- sample(pair_attr_late$mobility_wil) ## randomly shuffle mobility varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                + scale(mid.year)*scale(mob_shuffle) + (1|spp) + (1|pair.id), data=pair_attr_late)
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
ukbms_mob_para2$parameter <- paste(row.names(ukbms_mob_para2)) ## move row.names to parameter column
rownames(ukbms_mob_para2) <- 1:nrow(ukbms_mob_para2) ## change row names to numbers
ukbms_mob_para2 <- ukbms_mob_para2[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility*midyear interaction
ukbms_mob_para2 <- ukbms_mob_para2 %>% filter_all(any_vars(grepl(":", .)))

final_results_table <- rbind(main_result_table, ukbms_mob_para2) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_mob_ukbms_00_12.csv", row.names=TRUE)
## read in file
perm_change_mob_ukbms <- read.csv("../Results/UKBMS/perm_change_mob_ukbms_00_12.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_ukbms <- perm_change_mob_ukbms[perm_change_mob_ukbms$i==0,]
perm_change_mob_ukbms <- perm_change_mob_ukbms[!perm_change_mob_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant




#################################################
### 3a. Average synchrony + average abundance ### 
#################################################

pair_attr_abund <- pair_attr[!is.na(pair_attr$average_abundance),]
length(unique(pair_attr_abund$spp)) # 31 spp
str(pair_attr_abund)

avg_abund_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(mid.year) + scale(average_abundance) + (1|pair.id) + (1|spp), data = pair_attr_abund)


summary(avg_abund_ukbms)
## non-significant
results_table_abund_ukbms <- data.frame(summary(avg_abund_ukbms)$coefficients[,1:5]) ## 31 species
confidence_intervals <- data.frame(confint(avg_abund_ukbms, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_abund_ukbms <- cbind(results_table_abund_ukbms,confidence_intervals)

write.csv(results_table_abund_ukbms, file = "../Results/UKBMS/average_abund_ukbms.csv", row.names=TRUE)

#############################################################
### 3b. Change in synchrony + non-sig change in abundance ### 
#############################################################

pair_attr_early <- droplevels(pair_attr[pair_attr$mid.year==1984.5 | pair_attr$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr[pair_attr$mid.year==1999.5 | pair_attr$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

pair_attr_early <- pair_attr_early[!is.na(pair_attr_early$abund_change_85_00),]
pair_attr_late <- pair_attr_late[!is.na(pair_attr_late$abund_change_00_12),]

length(unique(pair_attr_early$spp))
length(unique(pair_attr_late$spp))

#### Early 1985-2000 ####
abund_mod1 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(mid.year)*abund_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_early)

summary(abund_mod1)
## not significant (p=0.22)
results_table_abund_ukbms <- data.frame(summary(abund_mod1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/UKBMS/change_abund_85_00_ukbms.csv", row.names=TRUE)

#### Late 2000-2012 ####
abund_mod2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(mid.year)*abund_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_late)

summary(abund_mod2)
## significant
results_table_abund_ukbms <- data.frame(summary(abund_mod2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/UKBMS/change_abund_00_12_ukbms.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(abund_mod2, terms = c("mid.year", "abund_change_00_12"))
plot(dat) # species increasing in abundance increase in synchrony more rapidly

## run model without abundance for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_late$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                            + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                            + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with abundance info 
abundance <- pair_attr_late[,c("spp", "abund_change_00_12")]
abundance <- unique(abundance)
results_table <- merge(results_table, abundance, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Abundance")
results_table$Abundance <- revalue(results_table$Abundance, c("decrease"="Decrease"))
results_table$Abundance <- revalue(results_table$Abundance, c("increase"="Increase"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_abund_ukbms <- results_table_abund_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_abund_ukbms$Abundance <- paste(row.names(results_table_abund_ukbms))
rownames(results_table_abund_ukbms) <- 1:nrow(results_table_abund_ukbms)
results_table_abund_ukbms <- results_table_abund_ukbms[c(12,14),]
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year)"="Decrease"))
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year):abund_change_00_12increase"="Increase"))

## create new dataframe
decrease_slope <- results_table_abund_ukbms[1,1]
increase_slope <- sum(results_table_abund_ukbms$Estimate)
decrease_SE <- results_table_abund_ukbms[1,2]
interaction_SE <- results_table_abund_ukbms[2,2]
increase_SE <- sqrt((decrease_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Abundance=c("Decrease", "Increase"), slope=c(decrease_slope, increase_slope), SE=c(decrease_SE, increase_SE))

model_summary$Abundance <- factor(model_summary$Abundance, levels=c("Increase", "Decrease"))
results_table$Abundance <- factor(results_table$Abundance, levels=c("Increase", "Decrease"))

png("../Graphs/New/Abundance_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
abund_00_12 <- ggplot(mapping=aes(x=Abundance, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71", position=myjit) +
  geom_errorbar(data=results_table, color="grey71", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Change in abundance", y="Change in population \n synchrony 2000-2012") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
abund_00_12
dev.off()

#### run permutations ####
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(abund_mod2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with abundance*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):abund_change_00_12"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_abund_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$abund_shuffle <- sample(pair_attr_late$abund_change_00_12) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                + scale(mid.year)*abund_shuffle + (1|spp) + (1|pair.id), data=pair_attr_late)
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
ukbms_abund_para$parameter <- paste(row.names(ukbms_abund_para)) ## move row.names to parameter column
rownames(ukbms_abund_para) <- 1:nrow(ukbms_abund_para) ## change row names to numbers
ukbms_abund_para <- ukbms_abund_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
ukbms_abund_para <- ukbms_abund_para %>% filter_all(any_vars(grepl(":", .)))

final_results_table <- rbind(main_result_table, ukbms_abund_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_abund_ukbms.csv", row.names=TRUE)
## read in file
perm_change_abund_ukbms <- read.csv("../Results/UKBMS/perm_change_abund_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_abund_ukbms <- perm_change_abund_ukbms[perm_change_abund_ukbms$i==0,]
perm_change_abund_ukbms <- perm_change_abund_ukbms[!perm_change_abund_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_abund_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_abund_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant


#########################################################
### 3c. Change in synchrony + sig change in abundance ### 
#########################################################

pair_attr_early <- pair_attr_early[!pair_attr_early$sig_abund_change_85_00=="No change",]
pair_attr_late <- pair_attr_late[!pair_attr_late$sig_abund_change_00_12=="No change",] 

length(unique(pair_attr_early$spp)) # 9 species
length(unique(pair_attr_late$spp)) # 10 species

#### Early 1985-2000 ####
abund_mod3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(mid.year)*sig_abund_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_early)

summary(abund_mod3)
## significant
results_table_abund_ukbms <- data.frame(summary(abund_mod3)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/UKBMS/change_sig_abund_85_00_ukbms.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(abund_mod3, terms = c("mid.year", "sig_abund_change_85_00"))
plot(dat) # species decreasing in abundance decline in synchrony more rapidly

## run model without abundance for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_early$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                            + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                            + scale(mid.year) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(mid.year) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(mid.year), data=pair_attr_early[pair_attr_early$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with abundance info 
abundance <- pair_attr_early[,c("spp", "sig_abund_change_85_00")]
abundance <- unique(abundance)
results_table <- merge(results_table, abundance, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Abundance")

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_abund_ukbms <- results_table_abund_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_abund_ukbms$Abundance <- paste(row.names(results_table_abund_ukbms))
rownames(results_table_abund_ukbms) <- 1:nrow(results_table_abund_ukbms)
results_table_abund_ukbms <- results_table_abund_ukbms[c(12,14),]
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year)"="Decrease"))
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year):sig_abund_change_85_00Increase"="Increase"))

## create new dataframe
decrease_slope <- results_table_abund_ukbms[1,1]
increase_slope <- sum(results_table_abund_ukbms$Estimate)
decrease_SE <- results_table_abund_ukbms[1,2]
interaction_SE <- results_table_abund_ukbms[2,2]
increase_SE <- sqrt((decrease_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Abundance=c("Decrease", "Increase"), slope=c(decrease_slope, increase_slope), SE=c(decrease_SE, increase_SE))

model_summary$Abundance <- factor(model_summary$Abundance, levels=c("Increase", "Decrease"))
results_table$Abundance <- factor(results_table$Abundance, levels=c("Increase", "Decrease"))

png("../Graphs/New/Abundance_sig_change_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
abund_85_00 <- ggplot(mapping=aes(x=Abundance, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71", position=myjit) +
  geom_errorbar(data=results_table, color="grey71", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Change in abundance", y="Change in population \n synchrony 1985-2000") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
abund_85_00
dev.off()

## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(abund_mod3)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with abundance*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):sig_abund_change_85_00"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_abund_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_early$abund_shuffle <- sample(pair_attr_early$sig_abund_change_85_00) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                + scale(mid.year)*abund_shuffle + (1|spp) + (1|pair.id), data=pair_attr_early)
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
ukbms_abund_para$parameter <- paste(row.names(ukbms_abund_para)) ## move row.names to parameter column
rownames(ukbms_abund_para) <- 1:nrow(ukbms_abund_para) ## change row names to numbers
ukbms_abund_para <- ukbms_abund_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
ukbms_abund_para <- ukbms_abund_para %>% filter_all(any_vars(grepl(":", .)))

final_results_table <- rbind(main_result_table, ukbms_abund_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_abund_ukbms_85_00_3cat.csv", row.names=TRUE)
## read in file
perm_change_abund_ukbms <- read.csv("../Results/UKBMS/perm_change_abund_ukbms_85_00_3cat.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_abund_ukbms <- perm_change_abund_ukbms[perm_change_abund_ukbms$i==0,]
perm_change_abund_ukbms <- perm_change_abund_ukbms[!perm_change_abund_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_abund_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_abund_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.016 significant



#### Late 2000-2012 ####
abund_mod4 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(mid.year)*sig_abund_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_late)

summary(abund_mod4)
## significant
results_table_abund_ukbms <- data.frame(summary(abund_mod4)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/UKBMS/change_sig_abund_00_12_ukbms.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(abund_mod4, terms = c("mid.year", "sig_abund_change_00_12"))
plot(dat) # species declining in abundance increase in synchrony more rapidly (opposite of hypothesis)
# even non-sig increasing in abundance drive emigration of individuals


## run model without abundance for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_late$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                            + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                            + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with abundance info 
abundance <- pair_attr_late[,c("spp", "sig_abund_change_00_12")]
abundance <- unique(abundance)
results_table <- merge(results_table, abundance, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Abundance")

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_abund_ukbms <- results_table_abund_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_abund_ukbms$Abundance <- paste(row.names(results_table_abund_ukbms))
rownames(results_table_abund_ukbms) <- 1:nrow(results_table_abund_ukbms)
results_table_abund_ukbms <- results_table_abund_ukbms[c(12,14),]
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year)"="Decrease"))
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("scale(mid.year):sig_abund_change_00_12Increase"="Increase"))

## create new dataframe
decrease_slope <- results_table_abund_ukbms[1,1]
increase_slope <- sum(results_table_abund_ukbms$Estimate)
decrease_SE <- results_table_abund_ukbms[1,2]
interaction_SE <- results_table_abund_ukbms[2,2]
increase_SE <- sqrt((decrease_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Abundance=c("Decrease", "Increase"), slope=c(decrease_slope, increase_slope), SE=c(decrease_SE, increase_SE))

model_summary$Abundance <- factor(model_summary$Abundance, levels=c("Increase", "Decrease"))
results_table$Abundance <- factor(results_table$Abundance, levels=c("Increase", "Decrease"))

png("../Graphs/New/Abundance_sig_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
abund_00_12 <- ggplot(mapping=aes(x=Abundance, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71", position=myjit) +
  geom_errorbar(data=results_table, color="grey71", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Change in abundance", y="Change in population \n synchrony 2000-2012") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
abund_00_12
dev.off()

#### run permutations ####

## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(abund_mod4)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with abundance*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):sig_abund_change_00_12"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_abund_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$abund_shuffle <- sample(pair_attr_late$sig_abund_change_00_12) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                + scale(mid.year)*abund_shuffle + (1|spp) + (1|pair.id), data=pair_attr_late)
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
ukbms_abund_para$parameter <- paste(row.names(ukbms_abund_para)) ## move row.names to parameter column
rownames(ukbms_abund_para) <- 1:nrow(ukbms_abund_para) ## change row names to numbers
ukbms_abund_para <- ukbms_abund_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
ukbms_abund_para <- ukbms_abund_para %>% filter_all(any_vars(grepl(":", .)))

final_results_table <- rbind(main_result_table, ukbms_abund_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/UKBMS/perm_change_abund_ukbms_00_12_3cat.csv", row.names=TRUE)
## read in file
perm_change_abund_ukbms <- read.csv("../Results/UKBMS/perm_change_abund_ukbms_00_12_3cat.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_abund_ukbms <- perm_change_abund_ukbms[perm_change_abund_ukbms$i==0,]
perm_change_abund_ukbms <- perm_change_abund_ukbms[!perm_change_abund_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_abund_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_abund_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant







####### Figure 1 plot

# ukbms specialism
ukbms_spec <- read.csv("../Results/UKBMS/average_spec_ukbms.csv", header=TRUE)
ukbms_spec <- ukbms_spec[13,]
ukbms_spec <- ukbms_spec[,c(1:2,7:8)]
ukbms_spec$parameter <- "Specialism"
ukbms_spec$Scheme <- "UKBMS"
# ukbms mobility
ukbms_mob <- read.csv("../Results/UKBMS/average_mob_ukbms.csv", header=TRUE)
ukbms_mob <- ukbms_mob[13,]
ukbms_mob <- ukbms_mob[,c(1:2,7:8)]
ukbms_mob$parameter <- "Mobility"
ukbms_mob$Scheme <- "UKBMS"
# ukbms abundance
ukbms_abund <- read.csv("../Results/UKBMS/average_abund_ukbms.csv", header=TRUE)
ukbms_abund <- ukbms_abund[13,]
ukbms_abund <- ukbms_abund[,c(1:2,7:8)]
ukbms_abund$parameter <- "Abundance"
ukbms_abund$Scheme <- "UKBMS"

# cbc specialism
cbc_spec <- read.csv("../Results/CBC/average_spec_cbc.csv", header=TRUE)
cbc_spec <- cbc_spec[7,]
cbc_spec <- cbc_spec[,c(1:2,7:8)]
cbc_spec$parameter <- "Specialism"
cbc_spec$Scheme <- "CBC"
# cbc mobility
cbc_mob <- read.csv("../Results/CBC/average_mob_cbc.csv", header=TRUE)
cbc_mob <- cbc_mob[7,]
cbc_mob <- cbc_mob[,c(1:2,7:8)]
cbc_mob$parameter <- "Mobility"
cbc_mob$Scheme <- "CBC"
# cbc abundance
cbc_abund <- read.csv("../Results/CBC/average_abund_cbc.csv", header=TRUE)
cbc_abund <- cbc_abund[7,]
cbc_abund <- cbc_abund[,c(1:2,7:8)]
cbc_abund$parameter <- "Abundance"
cbc_abund$Scheme <- "CBC"

# bbs specialist
bbs_spec <- read.csv("../Results/BBS/average_spec_bbs.csv", header=TRUE)
bbs_spec <- bbs_spec[9,]
bbs_spec <- bbs_spec[,c(1:2,7:8)]
bbs_spec$parameter <- "Specialism"
bbs_spec$Scheme <- "BBS"
# bbs mobility
bbs_mob <- read.csv("../Results/BBS/average_mob_bbs.csv", header=TRUE)
bbs_mob <- bbs_mob[9,]
bbs_mob <- bbs_mob[,c(1:2,7:8)]
bbs_mob$parameter <- "Mobility"
bbs_mob$Scheme <- "BBS"
# bbs abundance
bbs_abund <- read.csv("../Results/BBS/average_abund_bbs.csv", header=TRUE)
bbs_abund <- bbs_abund[9,]
bbs_abund <- bbs_abund[,c(1:2,7:8)]
bbs_abund$parameter <- "Abundance"
bbs_abund$Scheme <- "BBS"


final <- rbind(ukbms_spec, ukbms_mob, ukbms_abund, cbc_spec, cbc_mob, cbc_abund, bbs_spec, bbs_mob, bbs_abund)
final$Scheme <- factor(final$Scheme, levels=c("UKBMS", "CBC", "BBS") )
final$parameter <- factor(final$parameter, levels=c("Specialism", "Mobility", "Abundance") )

plot <- ggplot(final, aes(x=parameter, y=Estimate))+
  geom_point(aes(shape=Scheme),  position=position_dodge(width=0.7), size=5)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, group=Scheme),  position=position_dodge(width=0.7), width=0.2)+
  geom_hline(yintercept = 0, linetype="dashed")+
  labs(x="Species attributes", y="Standardised coefficient")+
  theme_classic()+
  theme(text=element_text(size=25))
plot
ggsave(plot, file="../Graphs/New/Figure1.png", height=7, width=12)



### get trends in synchrony for each species
pair_attr <- readRDS("pop_climate_synchrony_UKBMS_final.rds") # butterfly pair attribute data
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)

spp <- unique(pair_attr$spp)
results_table_final <- NULL
for(i in spp){
  print(i)
mod <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
       scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
       scale(autumn_temp) + scale(spring_temp) + mid.year +
       (1|pair.id)-1, data = pair_attr[pair_attr$spp==i,])
results_table_temp <- data.frame(summary(mod)$coefficients[,1:5],i) ## save model results
results_table_final <-rbind(results_table_final,results_table_temp)

}

## save results
names(results_table_final) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table_final$parameter <- paste(row.names(results_table_final))
rownames(results_table_final) <- 1:nrow(results_table_final)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_final <- results_table_final[grep("mid.year", results_table_final$parameter),]
results_table_final <- results_table_final[,c(1,6:7)]
results_table_final$year <- rep(1985:2012)

for(i in spp){
  
  plot <- ggplot(results_table_final[results_table_final$species==i,], aes(y=Estimate, x=year))+
    geom_smooth(colour="black")+
    labs(x="Year", y="Residual population synchrony")+
    theme_classic()
}

############ jitter code #################
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.5,
                 dodge.width = 0.15,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )
