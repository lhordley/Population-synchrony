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
options(scipen=999)

## load data
pair_attr <- readRDS("../Data/Butterfly_sync_data/pop_climate_synchrony_UKBMS_final.rds") # butterfly pair attribute data
str(pair_attr)
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$abund_change_00_12 <- as.factor(pair_attr$abund_change_00_12)
pair_attr$abund_change_85_00 <- as.factor(pair_attr$abund_change_85_00)
pair_attr$sig_abund_change_00_12 <- as.factor(pair_attr$sig_abund_change_00_12)
pair_attr$sig_abund_change_85_00 <- as.factor(pair_attr$sig_abund_change_85_00)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$specialism <- as.factor(pair_attr$specialism)

summary(pair_attr)

### 1a. Average synchrony + specialism 
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
length(unique(pair_attr$spp)) # 32 spp

spec_model_ukbms1 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) + mid.year + specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr)
summary(spec_model_ukbms1)
# non-significant

results_table_average_spec <- data.frame(summary(spec_model_ukbms1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_average_spec, file = "../Results/Model_outputs/UKBMS/average_spec_ukbms.csv", row.names=TRUE) # 32 species

### 1b. Change in synchrony + specialism
pair_attr_early <- droplevels(pair_attr[pair_attr$mid.year==1984.5 | pair_attr$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr[pair_attr$mid.year==1999.5 | pair_attr$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

# Early (1985-2000)
spec_model_ukbms2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) + scale(mid.year)*specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr_early)
summary(spec_model_ukbms2)
# significant
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms, file = "../Results/Model_outputs/UKBMS/change_spec_ukbms_85_00.csv", row.names=TRUE)

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
                              scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) +(1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                         scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                         scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                       scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                       scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp), data=pair_attr_early[pair_attr_early$spp==i,])
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
results_table_spec_ukbms <- results_table_spec_ukbms[c(13,15),] # keep mid year and interaction rows
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

png("../Graphs/Specialism/Specialism_change_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
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

# Late (2000-2012)
spec_model_ukbms3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + 
                            scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) +
                            scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) + scale(mid.year)*specialism + 
                            (1|pair.id) + (1|spp), data = pair_attr_late)
summary(spec_model_ukbms3)
# significant 
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms3)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms, file = "../Results/Model_outputs/UKBMS/change_spec_ukbms_00_12.csv", row.names=TRUE)

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
                              scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) +(1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                         scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                         scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(mid.year) + 
                       scale(winter_rain) + scale(autumn_rain) + scale(spring_rain) + scale(summer_rain) + 
                       scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + scale(summer_temp), data=pair_attr_late[pair_attr_late$spp==i,])
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
results_table_spec_ukbms <- results_table_spec_ukbms[c(13,15),] # keep mid year and interaction rows
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

png("../Graphs/Specialism/Specialism_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
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

#### 2a. Average synchrony + mobility
str(pair_attr)
pair_attr_mob <- pair_attr[!is.na(pair_attr$mobility_wil),]
pair_attr_mob$mid.year <- as.factor(pair_attr_mob$mid.year)
length(unique(pair_attr_mob$spp)) # 31 spp

mobility_model <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + mid.year + scale(winter_rain) + scale(autumn_rain) + 
                         scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp) + 
                         scale(summer_temp) + scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_mob)
summary(mobility_model)
## significant (positive)
results_table_mob_ukbms <- data.frame(summary(mobility_model)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mob_ukbms, file = "../Results/Model_outputs/UKBMS/average_mob_ukbms.csv", row.names=TRUE)

#### 2b. Change in synchrony + mobility

# Early
pair_attr_early <- droplevels(pair_attr_mob[pair_attr_mob$mid.year==1984.5 | pair_attr_mob$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr_mob[pair_attr_mob$mid.year==1999.5 | pair_attr_mob$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

mobility_model2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(summer_temp) + scale(mid.year)*scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(mobility_model2)
## interaction non-significant
results_table_mobility_ukbms <- data.frame(summary(mobility_model2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_85_00.csv", row.names=TRUE)

# Late
mobility_model3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(summer_temp) + scale(mid.year)*scale(mobility_wil) + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(mobility_model3)
## interaction significant
results_table_mobility_ukbms <- data.frame(summary(mobility_model3)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_00_12.csv", row.names=TRUE)

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
                           + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    mob_ukbms <- lmer(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                      + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                      + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    mob_ukbms <- lm(lag0 ~  scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                    + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                    + scale(summer_temp) + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
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
results_table_mobility_ukbms <- results_table_mobility_ukbms[c(13,15),]

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

png("../Graphs/FINAL/FigureX_butterflies_mob.png", height = 150, width = 180, units = "mm", res = 300)
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
png("../Graphs/FINAL/Figure2_sync.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(spec, mob, 
          labels = c("(a)", "(b)"), font.label = list(size = 16, color ="black"),
          ncol = 2, nrow = 1)
dev.off()


#### 3a. Average synchrony + average abundance
pair_attr_abund <- pair_attr[!is.na(pair_attr$average_abundance),]
pair_attr_abund$mid.year <- as.factor(pair_attr_abund$mid.year)
length(unique(pair_attr_abund$spp)) # 31 spp
str(pair_attr_abund)

avg_abund_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                        + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                        + scale(summer_temp) + mid.year + scale(average_abundance) + (1|pair.id) + (1|spp), data = pair_attr_abund)


summary(avg_abund_ukbms)
## non-significant
results_table_abund_ukbms <- data.frame(summary(avg_abund_ukbms)$coefficients[,1:5]) ## 31 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/average_abund_ukbms.csv", row.names=TRUE)

qqnorm(resid(avg_abund_ukbms))
qqline(resid(avg_abund_ukbms))
plot(avg_abund_ukbms, which = 1)

#### 3b. Change in synchrony + non-significant abundance change

pair_attr_early <- droplevels(pair_attr[pair_attr$mid.year==1984.5 | pair_attr$mid.year==1999.5,])
pair_attr_late <- droplevels(pair_attr[pair_attr$mid.year==1999.5 | pair_attr$mid.year==2011.5,])

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)

pair_attr_early <- pair_attr_early[!is.na(pair_attr_early$abund_change_85_00),]
pair_attr_late <- pair_attr_late[!is.na(pair_attr_late$abund_change_00_12),]

length(unique(pair_attr_early$spp))
length(unique(pair_attr_late$spp))

# Early
abund_mod1 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(summer_temp) + scale(mid.year)*abund_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_early)

summary(abund_mod1)
## not significant (p=0.22)
results_table_abund_ukbms <- data.frame(summary(abund_mod1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_85_00_ukbms.csv", row.names=TRUE)

qqnorm(resid(abund_mod1))
qqline(resid(abund_mod1))
plot(abund_mod1, which = 1)

# Late
abund_mod2 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(summer_temp) + scale(mid.year)*abund_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_late)

summary(abund_mod2)
## significant
results_table_abund_ukbms <- data.frame(summary(abund_mod2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_00_12_ukbms.csv", row.names=TRUE)

qqnorm(resid(abund_mod2))
qqline(resid(abund_mod2))
plot(abund_mod2, which = 1)

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
                            + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(summer_temp) + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
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
results_table_abund_ukbms <- results_table_abund_ukbms[c(13,15),]
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

png("../Graphs/Abundance/Abundance_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
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


#### 3c. Change in synchrony + significant abundance change
pair_attr_early <- pair_attr_early[!pair_attr_early$sig_abund_change_85_00=="No change",]
pair_attr_late <- pair_attr_late[!pair_attr_late$sig_abund_change_00_12=="No change",] 

length(unique(pair_attr_early$spp)) # 9 species
length(unique(pair_attr_late$spp)) # 10 species

# Early
abund_mod3 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(summer_temp) + scale(mid.year)*sig_abund_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_early)

summary(abund_mod3)
## significant
results_table_abund_ukbms <- data.frame(summary(abund_mod3)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/sig_change_abund_85_00_ukbms.csv", row.names=TRUE)

qqnorm(resid(abund_mod3))
qqline(resid(abund_mod3))
plot(abund_mod3, which = 1)

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
                            + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(summer_temp) + scale(mid.year), data=pair_attr_early[pair_attr_early$spp==i,])
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
results_table_abund_ukbms <- results_table_abund_ukbms[c(13,15),]
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

png("../Graphs/Abundance/Abundance_sig_change_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
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


# Late
abund_mod4 <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                   + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                   + scale(summer_temp) + scale(mid.year)*sig_abund_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_late)

summary(abund_mod4)
## not significant (p=0.22)
results_table_abund_ukbms <- data.frame(summary(abund_mod4)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/sig_change_abund_00_12_ukbms.csv", row.names=TRUE)

qqnorm(resid(abund_mod4))
qqline(resid(abund_mod4))
plot(abund_mod4, which = 1)

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
                            + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                       + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                       + scale(summer_temp) + scale(mid.year) + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ scale(mean_northing) + scale(distance) + scale(renk_hab_sim) + scale(winter_rain) + scale(autumn_rain)
                     + scale(spring_rain) + scale(summer_rain) + scale(winter_temp) + scale(autumn_temp) + scale(spring_temp)
                     + scale(summer_temp) + scale(mid.year), data=pair_attr_late[pair_attr_late$spp==i,])
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
results_table_abund_ukbms <- results_table_abund_ukbms[c(13,15),]
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

png("../Graphs/Abundance/Abundance_sig_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
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
