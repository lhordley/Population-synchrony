# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in synchronous population dynamics 
# BBS dataset analysis

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
pair_attr_BBS <- readRDS("pop_climate_synchrony_BBS_final.rds") # BBS bird pair attribute data
str(pair_attr_BBS)
pair_attr_BBS$pop_estimate <- as.numeric(pair_attr_BBS$pop_estimate)
pair_attr_BBS$ab_change_85_96 <- as.factor(pair_attr_BBS$ab_change_99_12)
pair_attr_BBS$sig_ab_change_85_96 <- as.factor(pair_attr_BBS$sig_ab_change_99_12)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)
pair_attr_BBS$specialism <- as.factor(pair_attr_BBS$specialism)

summary(pair_attr_BBS)

##########################################
### 1a. Average synchrony + specialism ### 
##########################################

length(unique(pair_attr_BBS$spp)) # 24 spp

strategy_model_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                             scale(spring_rain) + scale(autumn_temp) + specialism + (1|pair.id) + (1|spp), data = pair_attr_BBS)

summary(strategy_model_bbs)
anova(strategy_model_bbs)
## non-significant (p=0.83)
results_table_average_spec <- data.frame(summary(strategy_model_bbs)$coefficients[,1:5]) 
confidence_intervals <- data.frame(confint(strategy_model_bbs, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_average_spec <- cbind(results_table_average_spec,confidence_intervals)

write.csv(results_table_average_spec, file = "../Results/BBS/average_spec_bbs.csv", row.names=TRUE) # 24 species

############################################
### 1b. Change in synchrony + specialism ### 
############################################
pair_attr_bbs <- droplevels(pair_attr_BBS[pair_attr_BBS$mid.year==1998.5 | pair_attr_BBS$mid.year==2011.5,])

pair_attr_bbs$mid.year <- as.numeric(pair_attr_bbs$mid.year)
length(unique(pair_attr_bbs$spp)) # 24 spp

spec_model_bbs2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(spring_rain) + 
                          scale(autumn_temp) + scale(mid.year)*specialism + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(spec_model_bbs2)
## just significant! (p=0.045)
results_table_strategy_bbs <- data.frame(summary(spec_model_bbs2)$coefficients[,1:5]) ## 24 species
write.csv(results_table_strategy_bbs, file = "../Results/BBS/change_spec_bbs.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(spec_model_bbs2, terms = c("mid.year", "specialism"))
plot(dat) # generalists decline in synchrony, specialists show no change

## run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_bbs$spp)){
  print(i)
  # if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model
  loadError=F
  a=try({spec_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                              scale(autumn_rain) + scale(spring_rain) + (1|pair.id), data=pair_attr_bbs[pair_attr_bbs$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                         scale(autumn_rain) + scale(spring_rain) + (1|pair.id), data=pair_attr_bbs[pair_attr_bbs$spp==i,])
  }else{
    spec_bbs <- lm(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                       scale(autumn_rain) + scale(spring_rain), data=pair_attr_bbs[pair_attr_bbs$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_bbs)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info 
specialism <- pair_attr_bbs[,c("spp", "specialism")]
specialism <- unique(specialism)
results_table <- merge(results_table, specialism, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Specialism")
results_table$Specialism <- revalue(results_table$Specialism, c("generalist"="Generalist"))
results_table$Specialism <- revalue(results_table$Specialism, c("specialist"="Specialist"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_strategy_bbs <- results_table_strategy_bbs[,-c(3:5)]
## leave rows with year and interaction
results_table_strategy_bbs$Specialism <- paste(row.names(results_table_strategy_bbs))
rownames(results_table_strategy_bbs) <- 1:nrow(results_table_strategy_bbs)
results_table_strategy_bbs <- results_table_strategy_bbs[c(8,10),]
results_table_strategy_bbs$Specialism <- revalue(results_table_strategy_bbs$Specialism, c("scale(mid.year)"="Generalist"))
results_table_strategy_bbs$Specialism <- revalue(results_table_strategy_bbs$Specialism, c("scale(mid.year):specialismspecialist"="Specialist"))

## create new dataframe
generalist_slope <- results_table_strategy_bbs[1,1]
specialist_slope <- sum(results_table_strategy_bbs$Estimate)
generalist_SE <- results_table_strategy_bbs[1,2]
interaction_SE <- results_table_strategy_bbs[2,2]
specialist_SE <- sqrt((generalist_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Specialism=c("Generalist", "Specialist"), slope=c(generalist_slope, specialist_slope), SE=c(generalist_SE, specialist_SE))

model_summary$Specialism <- factor(model_summary$Specialism, levels=c("Generalist", "Specialist"))
results_table$Specialism <- factor(results_table$Specialism, levels=c("Generalist", "Specialist"))
results_table <- na.omit(results_table)

png("../Graphs/Specialism/Specialism_change_bbs.png", height = 150, width = 180, units = "mm", res = 300)
spec <- ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="lightgrey", position=myjit) +
  geom_errorbar(data=results_table, color="lightgrey", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Specialism", y="Change in population \n synchrony 1999-2012") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "lightgrey", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
spec
dev.off()

#### run permutations ####

## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_bbs2)[5]) ## save anova table from main model
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
bbs_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs$spec_shuffle <- sample(pair_attr_bbs$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(spring_rain) + 
                  scale(autumn_temp) + scale(mid.year)*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
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
bbs_spec_para$parameter <- paste(row.names(bbs_spec_para)) ## move row.names to parameter column
rownames(bbs_spec_para) <- 1:nrow(bbs_spec_para) ## change row names to numbers
bbs_spec_para <- bbs_spec_para[,-c(1:3)] ## remove unnecessary columns
## only interested in specialism interaction
bbs_spec_para <- bbs_spec_para[grep("mid.year):spec_shuffle", bbs_spec_para$parameter),]

final_results_table <- rbind(main_result_table, bbs_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/BBS/perm_change_spec_bbs.csv", row.names=TRUE)
## read in file
perm_change_spec_bbs <- read.csv("../Results/BBS/perm_change_spec_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_spec_bbs <- perm_change_spec_bbs[perm_change_spec_bbs$i==0,]
perm_change_spec_bbs <- perm_change_spec_bbs[!perm_change_spec_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.049 still *just* significant



########################################
### 2a. Average synchrony + mobility ### 
########################################

pair_attr_BBS_mob <- pair_attr_BBS[!is.na(pair_attr_BBS$Breeding_AM),]
str(pair_attr_BBS_mob)
length(unique(pair_attr_BBS_mob$spp)) # 17 spp

dispersal_model_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                               scale(spring_rain) + scale(autumn_temp) + scale(Breeding_AM) + (1|spp) + (1|pair.id), data=pair_attr_BBS_mob)
summary(dispersal_model_bbs)
## not significant 
results_table_mob_bbs <- data.frame(summary(dispersal_model_bbs)$coefficients[,1:5]) ## 17 species
confidence_intervals <- data.frame(confint(dispersal_model_bbs, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_mob_bbs <- cbind(results_table_mob_bbs,confidence_intervals)

write.csv(results_table_mob_bbs, file = "../Results/BBS/average_mob_bbs.csv", row.names=TRUE)

##########################################
### 2b. Change in synchrony + mobility ### 
##########################################

pair_attr_bbs_mob <- droplevels(pair_attr_BBS_mob[pair_attr_BBS_mob$mid.year==1998.5 | pair_attr_BBS_mob$mid.year==2011.5,])

pair_attr_bbs_mob$mid.year <- as.numeric(pair_attr_bbs_mob$mid.year)
length(unique(pair_attr_bbs_mob$spp)) # 17 spp

dispersal_model_bbs2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(spring_rain) + 
       scale(autumn_temp) + scale(mid.year)*scale(Breeding_AM) + (1|pair.id) + (1|spp), data = pair_attr_bbs_mob)
summary(dispersal_model_bbs2)
## mid.year*dispersal interaction is significant overall
results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs2)$coefficients[,1:5])
write.csv(results_table_dispersal_bbs, file = "../Results/BBS/change_mob_bbs.csv", row.names=TRUE)

qqnorm(resid(dispersal_model_bbs2))
qqline(resid(dispersal_model_bbs2))
plot(dispersal_model_bbs2, which = 1)

# plot result to show direction
dat <- ggpredict(dispersal_model_bbs2, terms = c("mid.year", "Breeding_AM"))
plot(dat) # less mobile species decline in synchrony over time, more mobile remain stable

simulationOutput <- simulateResiduals(fittedModel = dispersal_model_bbs2)
plot(simulationOutput)

## run model without mobility for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_bbs_mob$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({mob_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                           scale(spring_rain) + scale(autumn_rain) + scale(mid.year) +
                           (1|pair.id), data=pair_attr_bbs_mob[pair_attr_bbs_mob$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    mob_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                      scale(spring_rain) + scale(autumn_rain) + scale(mid.year) + (1|pair.id), data=pair_attr_bbs_mob[pair_attr_bbs_mob$spp==i,])
  }else{
    mob_bbs <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                    scale(spring_rain) + scale(autumn_rain) + scale(mid.year),
                  data=pair_attr_bbs_mob[pair_attr_bbs_mob$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(mob_bbs)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
mobility <- pair_attr_bbs_mob[,c("spp", "Breeding_AM")]
mobility <- unique(mobility)
results_table <- merge(results_table, mobility, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Mobility")


########## calculate slope and SE from main model (with mobility as x axis)
### plot points and SE ontop of raw data
results_table_dispersal_bbs <- results_table_dispersal_bbs[,-c(3:5)]
## leave rows with year and interaction
results_table_dispersal_bbs$Mobility <- paste(row.names(results_table_dispersal_bbs))
rownames(results_table_dispersal_bbs) <- 1:nrow(results_table_dispersal_bbs)
results_table_dispersal_bbs <- results_table_dispersal_bbs[c(8,10),]

## extract slope and SE for interaction and mid.year
mid.year_slope <- results_table_dispersal_bbs[1,1]
interaction_slope <- results_table_dispersal_bbs[2,1]
se_mid.year <- results_table_dispersal_bbs[1,2]
se_interaction <- results_table_dispersal_bbs[2,2]

## get slopes for each value of mobility
mobility<-unique(results_table$Mobility)
slopes <- mid.year_slope + interaction_slope * scale(mobility)
slopes

## get standard errors for each value of mobility
SEs <- rep(NA, length(mobility))
for (i in 1:length(mobility)){
  j <- mobility[i]  
  SEs[i] <- sqrt((se_interaction * j)^2 + se_mid.year^2) 
}
model_summary <- data.frame(Mobility=mobility, slope=slopes, SE=SEs)

png("../Graphs/New/FigureA5.png", height = 150, width = 180, units = "mm", res = 300)
ggplot(mapping=aes(x=Mobility, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71") +
  geom_errorbar(data=results_table, color="grey71", width=2) +
  geom_ribbon(data=model_summary, fill = "grey70", alpha=0.4) +
  geom_line(data=model_summary, size=1) +
  labs(x="Dispersal distance", y="Change in population \n synchrony 1999-2012") +
  #scale_x_continuous(breaks = seq(0,4.5,0.5)) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
dev.off()
## 4 species are driving this trend with high dispersal distances, try removing them and see if result is robust

## remove outlying species with high dispersal distances
pair_attr_bbs_mob2 <- pair_attr_bbs_mob[!(pair_attr_bbs_mob$species=="Blackcap" | pair_attr_bbs_mob$species=="Willow Warbler" |
                                            pair_attr_bbs_mob$species=="Robin" | pair_attr_bbs_mob$species=="Wren"),]

dispersal_model_bbs3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(spring_rain) + 
                               scale(autumn_temp) + scale(mid.year)*scale(Breeding_AM) + (1|pair.id) + (1|spp), data = pair_attr_bbs_mob2)
summary(dispersal_model_bbs3)
## interaction is now non-significant - outliers were driving the initial result
results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs3)$coefficients[,1:5])
write.csv(results_table_dispersal_bbs, file = "../Results/BBS/change_mob_bbs2.csv", row.names=TRUE)

#### run permutations on full model ####

## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(dispersal_model_bbs2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility*midyear interaction
main_result_table <- main_result_table[which(main_result_table$parameter=="scale(mid.year):scale(Breeding_AM)"),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_mob_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs_mob$mob_shuffle <- sample(pair_attr_bbs_mob$Breeding_AM) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(spring_rain) + 
                  scale(autumn_temp) + scale(mid.year)*scale(mob_shuffle) + (1|pair.id) + (1|spp), data = pair_attr_bbs_mob)
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
bbs_mob_para$parameter <- paste(row.names(bbs_mob_para)) ## move row.names to parameter column
rownames(bbs_mob_para) <- 1:nrow(bbs_mob_para) ## change row names to numbers
bbs_mob_para <- bbs_mob_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
bbs_mob_para <- bbs_mob_para %>% filter_all(any_vars(grepl(":", .)))

final_results_table <- rbind(main_result_table, bbs_mob_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
ggplot(final_results_table, aes(x=F.value)) +
  geom_histogram(fill= 'grey') +
  geom_vline(aes(xintercept= final_results_table$F.value[final_results_table$i==0]), colour="red")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()

## save file
write.csv(final_results_table, file = "../Results/BBS/perm_change_mob_bbs.csv", row.names=TRUE)
## read in file
perm_change_mob_bbs <- read.csv("../Results/BBS/perm_change_mob_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_bbs <- perm_change_mob_bbs[perm_change_mob_bbs$i==0,]
perm_change_mob_bbs <- perm_change_mob_bbs[!perm_change_mob_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant




#################################################
### 3a. Average synchrony + average abundance ### 
#################################################

common_model_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + scale(winter_rain) + 
                           scale(spring_rain) + scale(autumn_temp) + scale(pop_estimate) + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(common_model_bbs)
## non-significant
## save result
results_table_abund_bbs <- data.frame(summary(common_model_bbs)$coefficients[,1:5]) ## 24 species
confidence_intervals <- data.frame(confint(common_model_bbs, method="Wald"))
confidence_intervals <- na.omit(confidence_intervals)
colnames(confidence_intervals) <- c("lowerCI", "upperCI")
results_table_abund_bbs <- cbind(results_table_abund_bbs,confidence_intervals)

write.csv(results_table_abund_bbs, file = "../Results/BBS/average_abund_bbs.csv", row.names=TRUE)


#############################################################
### 3b. Change in synchrony + non-sig change in abundance ### 
#############################################################

pair_attr_bbs <- droplevels(pair_attr_BBS[pair_attr_BBS$mid.year==1998.5 | pair_attr_BBS$mid.year==2011.5,])
pair_attr_bbs <- pair_attr_bbs[!is.na(pair_attr_bbs$ab_change_99_12),]
pair_attr_bbs$mid.year <- as.numeric(pair_attr_bbs$mid.year)
length(unique(pair_attr_bbs$spp)) # 18 spp

abund_model_bbs <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                          scale(spring_rain) + scale(autumn_temp) + scale(mid.year)*ab_change_99_12 + (1|pair.id) + 
                          (1|spp), data = pair_attr_bbs)
summary(abund_model_bbs)
## non-significant with 3 categories
results_table_abund_bbs <- data.frame(summary(abund_model_bbs)$coefficients[,1:5]) ## 18 species
write.csv(results_table_abund_bbs, file = "../Results/BBS/change_abund_bbs.csv", row.names=TRUE)

#########################################################
### 3c. Change in synchrony + sig change in abundance ### 
#########################################################

pair_attr_bbs <- droplevels(pair_attr_BBS[pair_attr_BBS$mid.year==1998.5 | pair_attr_BBS$mid.year==2011.5,])
pair_attr_bbs <- pair_attr_bbs[!is.na(pair_attr_bbs$sig_ab_change_99_12),]
pair_attr_bbs$mid.year <- as.numeric(pair_attr_bbs$mid.year)

pair_attr_bbs <- pair_attr_bbs[!pair_attr_bbs$sig_ab_change_99_12=="No change",] 
length(unique(pair_attr_bbs$spp)) # 12 spp

abund_model_bbs2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                          scale(spring_rain) + scale(autumn_temp) + scale(mid.year)*ab_change_99_12 + (1|pair.id) + 
                          (1|spp), data = pair_attr_bbs)
summary(abund_model_bbs2)
## non-significant with 3 categories
results_table_abund_bbs <- data.frame(summary(abund_model_bbs2)$coefficients[,1:5]) ## 12 species
write.csv(results_table_abund_bbs, file = "../Results/BBS/change_sig_abund_bbs.csv", row.names=TRUE)


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
