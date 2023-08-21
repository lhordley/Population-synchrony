# --------------- Code --------------
# Disentangling how climate and dispersal drive temporal trends in spatial synchrony among populations
# CBC dataset analysis

rm(list=ls()) # clear R

library(plyr) # load packages
library(dplyr) 
library(lme4)
library(ggplot2)
library(lmerTest)
library(plotrix)
library(ggeffects)
library(DHARMa)
options(scipen=999)

# ## load data
# pair_attr_CBC <- readRDS("../Data/Bird_sync_data/pop_climate_synchrony_CBC.rds") # CBC pair attribute data
# bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)
# bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
# abundance_results_cbc1 <- read.csv("../Results/Bird_results/abundance_results_85_96.csv", header=TRUE)
# abundance_results_cbc2 <- read.csv("../Results/Bird_results/abundance_results_85_96_3cat.csv", header=TRUE)
# 
# bird_dispersal <- bird_dispersal[,c("Species_code", "Breeding_AM")]
# bird_common <- bird_common[,c("species_code", "pop_estimate")]
# abundance_results_cbc1 <- abundance_results_cbc1[,c("species_code", "ab_change_85_96")]
# abundance_results_cbc2 <- abundance_results_cbc2[,c("species_code", "abund_change")]
# colnames(abundance_results_cbc2)[2] <- "sig_ab_change_85_96"
# 
# attributes <- merge(bird_dispersal, bird_common, by.x="Species_code", by.y="species_code", all=TRUE)
# attributes <- merge(attributes, abundance_results_cbc1, by.x="Species_code", by.y="species_code", all=TRUE)
# attributes <- merge(attributes, abundance_results_cbc2, by.x="Species_code", by.y="species_code", all=TRUE)
# 
# pair_attr_CBC <- merge(pair_attr_CBC, attributes, by.x="spp", by.y="Species_code", all.x=TRUE)
# head(pair_attr_CBC)
# summary(pair_attr_CBC)
# 
# saveRDS(pair_attr_CBC, "../Data/Bird_sync_data/pop_climate_synchrony_CBC_final.rds")

## load data
pair_attr_CBC <- readRDS("../Data/Bird_sync_data/pop_climate_synchrony_CBC_final.rds") # CBC bird pair attribute data
str(pair_attr_CBC)
pair_attr_CBC$pop_estimate <- as.numeric(pair_attr_CBC$pop_estimate)
pair_attr_CBC$ab_change_85_96 <- as.factor(pair_attr_CBC$ab_change_85_96)
pair_attr_CBC$sig_ab_change_85_96 <- as.factor(pair_attr_CBC$sig_ab_change_85_96)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)
pair_attr_CBC$specialism <- as.factor(pair_attr_CBC$specialism)

summary(pair_attr_CBC)

### 1a. Average synchrony + specialism 
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)
length(unique(pair_attr_CBC$spp)) # 26 spp

strategy_model_cbc <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + mid.year + scale(summer_temp) + specialism + 
                             (1|family) + (1|pair.id) + (1|spp), data = pair_attr_CBC)
summary(strategy_model_cbc)
# non-significant
results_table_average_spec <- data.frame(summary(strategy_model_cbc)$coefficients[,1:5])
write.csv(results_table_average_spec, file = "../Results/Model_outputs/CBC/average_spec_cbc.csv", row.names=TRUE) # 26 species

### 1b. Change in synchrony + specialism
pair_attr_cbc <- droplevels(pair_attr_CBC[pair_attr_CBC$mid.year==1984.5 | pair_attr_CBC$mid.year==1995.5,])

pair_attr_cbc$mid.year <- as.numeric(pair_attr_cbc$mid.year)
length(unique(pair_attr_cbc$spp)) # 26 spp

spec_model_cbc2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(summer_temp) + 
                          scale(mid.year)*specialism + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(spec_model_cbc2)
# significant
results_table_strategy_cbc <- data.frame(summary(spec_model_cbc2)$coefficients[,1:5]) ## 26 species
write.csv(results_table_strategy_cbc, file = "../Results/Model_outputs/CBC/change_spec_cbc.csv", row.names=TRUE)

# plot result to find direction
dat <- ggpredict(spec_model_cbc2, terms = c("mid.year", "specialism"))
plot(dat) # specialists increase in synchrony, generalist small decline

## run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_cbc$spp)){
  print(i)
  # if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp) + (1|pair.id), data=pair_attr_cbc[pair_attr_cbc$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp) + (1|pair.id), data=pair_attr_cbc[pair_attr_cbc$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~  scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp), data=pair_attr_cbc[pair_attr_cbc$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info 
specialism <- pair_attr_cbc[,c("spp", "specialism")]
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
results_table_strategy_cbc <- results_table_strategy_cbc[,-c(3:5)]
## leave rows with year and interaction
results_table_strategy_cbc$Specialism <- paste(row.names(results_table_strategy_cbc))
rownames(results_table_strategy_cbc) <- 1:nrow(results_table_strategy_cbc)
results_table_strategy_cbc <- results_table_strategy_cbc[c(6,8),]
results_table_strategy_cbc$Specialism <- revalue(results_table_strategy_cbc$Specialism, c("scale(mid.year)"="Generalist"))
results_table_strategy_cbc$Specialism <- revalue(results_table_strategy_cbc$Specialism, c("scale(mid.year):specialismspecialist"="Specialist"))

## create new dataframe
generalist_slope <- results_table_strategy_cbc[1,1]
specialist_slope <- sum(results_table_strategy_cbc$Estimate)
generalist_SE <- results_table_strategy_cbc[1,2]
interaction_SE <- results_table_strategy_cbc[2,2]
specialist_SE <- sqrt((generalist_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Specialism=c("Generalist", "Specialist"), slope=c(generalist_slope, specialist_slope), SE=c(generalist_SE, specialist_SE))

model_summary$Specialism <- factor(model_summary$Specialism, levels=c("Generalist", "Specialist"))
results_table$Specialism <- factor(results_table$Specialism, levels=c("Generalist", "Specialist"))
results_table <- na.omit(results_table)

png("../Graphs/Specialism/Specialism_change_cbc.png", height = 150, width = 180, units = "mm", res = 300)
spec <- ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71", position=myjit) +
  geom_errorbar(data=results_table, color="grey71", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Specialism", y="Change in population \n synchrony 1985-1996") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
spec
dev.off()


#### 2a. Average synchrony + mobility
pair_attr_CBC_mob <- pair_attr_CBC[!is.na(pair_attr_CBC$Breeding_AM),]
str(pair_attr_CBC_mob)
length(unique(pair_attr_CBC_mob$spp)) # 20 spp

dispersal_model_cbc <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + mid.year + scale(summer_temp) + 
                              scale(Breeding_AM) + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_CBC_mob)
summary(dispersal_model_cbc)
# non-significant
results_table_mob_cbc <- data.frame(summary(dispersal_model_cbc)$coefficients[,1:5]) ## 20 species
write.csv(results_table_mob_cbc, file = "../Results/Model_outputs/CBC/average_mob_cbc.csv", row.names=TRUE)

qqnorm(resid(dispersal_model_cbc))
qqline(resid(dispersal_model_cbc))
plot(dispersal_model_cbc, which = 1)

### 2b. Change in synchrony + mobility
pair_attr_cbc_mob <- droplevels(pair_attr_CBC_mob[pair_attr_CBC_mob$mid.year==1984.5 | pair_attr_CBC_mob$mid.year==1995.5,])

pair_attr_cbc_mob$mid.year <- as.numeric(pair_attr_cbc_mob$mid.year)
length(unique(pair_attr_cbc_mob$spp)) # 20 spp

dispersal_model_cbc2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(summer_temp) + 
                               scale(mid.year)*scale(Breeding_AM) + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc_mob)

summary(dispersal_model_cbc2)
# non-significant
results_table_dispersal_cbc <- data.frame(summary(dispersal_model_cbc2)$coefficients[,1:5]) ## 20 species
write.csv(results_table_dispersal_cbc, file = "../Results/Model_outputs/CBC/change_mob_cbc.csv", row.names=TRUE)

#### 3a. Average synchrony + average abundance
common_model_cbc <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + mid.year + scale(summer_temp) + 
                           scale(pop_estimate) + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_CBC)

summary(common_model_cbc) ## pop estimate is significant (p=0.0128)
results_table_abund_cbc <- data.frame(summary(common_model_cbc)$coefficients[,1:5]) ## 26 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/average_abund_cbc.csv", row.names=TRUE)

qqnorm(resid(common_model_cbc))
qqline(resid(common_model_cbc))
plot(common_model_cbc, which = 1)

# plot result to show direction
dat <- ggpredict(common_model_cbc, terms = "pop_estimate")
plot(dat) # more common species have higher average levels of pop. synchrony

#### 3b. Change in synchrony + change in abundance
pair_attr_cbc <- droplevels(pair_attr_CBC[pair_attr_CBC$mid.year==1984.5 | pair_attr_CBC$mid.year==1995.5,])
pair_attr_cbc <- pair_attr_cbc[!is.na(pair_attr_cbc$ab_change_85_96),]
pair_attr_cbc$mid.year <- as.numeric(pair_attr_cbc$mid.year)
length(unique(pair_attr_cbc$spp)) # 22 spp

abund_model_cbc <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(summer_temp) + 
                          scale(mid.year)*ab_change_85_96 + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(abund_model_cbc)
# non-significant
results_table_abund_cbc <- data.frame(summary(abund_model_cbc)$coefficients[,1:5]) ## 22 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/change_abund_cbc.csv", row.names=TRUE)

qqnorm(resid(abund_model_cbc))
qqline(resid(abund_model_cbc))
plot(abund_model_cbc, which = 1)

#### 3b. Change in synchrony + change in abundance
pair_attr_cbc <- droplevels(pair_attr_CBC[pair_attr_CBC$mid.year==1984.5 | pair_attr_CBC$mid.year==1995.5,])
pair_attr_cbc <- pair_attr_cbc[!is.na(pair_attr_cbc$sig_ab_change_85_96),]
pair_attr_cbc$mid.year <- as.numeric(pair_attr_cbc$mid.year)

pair_attr_cbc <- pair_attr_cbc[!pair_attr_cbc$sig_ab_change_85_96=="No change",] 
length(unique(pair_attr_cbc$spp)) # 19 spp
# remove species 377 and 431 - v. high eror around the change in sync estimate
pair_attr_cbc2 <- pair_attr_cbc[!pair_attr_cbc$spp=="377",]
pair_attr_cbc2 <- pair_attr_cbc2[!pair_attr_cbc2$spp=="431",]
pair_attr_cbc2 <- droplevels(pair_attr_cbc2)
length(unique(pair_attr_cbc2$spp)) # 17 spp

abund_model_cbc2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(summer_temp) + 
                          scale(mid.year)*sig_ab_change_85_96 + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc2)
summary(abund_model_cbc2)
# significant (just)
results_table_abund_cbc <- data.frame(summary(abund_model_cbc2)$coefficients[,1:5]) ## 22 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/change_abund_cbc_3cat.csv", row.names=TRUE)

dat <- ggpredict(abund_model_cbc2, terms = c("mid.year", "sig_ab_change_85_96"))
plot(dat) # species significantly increasing in abundance increase in synchrony, whereas declining species decrease in synchrony

## run model without abundance for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_cbc2$spp)){
  print(i)
  # if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp) + (1|pair.id), data=pair_attr_cbc2[pair_attr_cbc2$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp) + (1|pair.id), data=pair_attr_cbc2[pair_attr_cbc2$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~  scale(mean_northing) + scale(distance) + scale(hab_sim) + scale(mid.year) + scale(summer_temp), data=pair_attr_cbc2[pair_attr_cbc2$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with abundance info 
abundance <- pair_attr_cbc2[,c("spp", "sig_ab_change_85_96")]
abundance <- unique(abundance)
results_table <- merge(results_table, abundance, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Abundance")

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_abund_cbc <- results_table_abund_cbc[,-c(3:5)]
## leave rows with year and interaction
results_table_abund_cbc$Abundance <- paste(row.names(results_table_abund_cbc))
rownames(results_table_abund_cbc) <- 1:nrow(results_table_abund_cbc)
results_table_abund_cbc <- results_table_abund_cbc[c(6,8),]
results_table_abund_cbc$Abundance <- revalue(results_table_abund_cbc$Abundance, c("scale(mid.year)"="Decrease"))
results_table_abund_cbc$Abundance <- revalue(results_table_abund_cbc$Abundance, c("scale(mid.year):sig_ab_change_85_96Increase"="Increase"))

## create new dataframe
decrease_slope <- results_table_abund_cbc[1,1]
increase_slope <- sum(results_table_abund_cbc$Estimate)
decrease_SE <- results_table_abund_cbc[1,2]
interaction_SE <- results_table_abund_cbc[2,2]
increase_SE <- sqrt((decrease_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Abundance=c("Decrease", "Increase"), slope=c(decrease_slope, increase_slope), SE=c(decrease_SE, increase_SE))

model_summary$Abundance <- factor(model_summary$Abundance, levels=c("Increase", "Decrease"))
results_table$Abundance <- factor(results_table$Abundance, levels=c("Increase", "Decrease"))
results_table <- na.omit(results_table)

png("../Graphs/Abundance/Abundance_change_cbc_85_96_3cat.png", height = 150, width = 180, units = "mm", res = 300)
abund2 <- ggplot(mapping=aes(x=Abundance, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="grey71", position=myjit) +
  geom_errorbar(data=results_table, color="grey71", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Change in abundance", y="Change in population \n synchrony 1985-1996") +
  #scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "grey71", size=1) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
abund2
dev.off()

library(ggpubr)
png("../Graphs/FINAL/Figure3_sync.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(abund, abund2, 
          labels = c("(a)", "(b)"), font.label = list(size = 16, color ="black"),
          ncol = 2, nrow = 1)
dev.off()
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
