#Author: Marion Chevrinais
#Date: August 2024

# required packages
library(car)
library(dplyr)
library(bbmle)
library(MuMIn)
library(lmerTest)
library(lsmeans)
library(rstatix)
library(ARTool)
library(ggplot2)
library(effects)

# Updload raw data ####
raw.data <- read.csv2(file.path(here::here(), "Raw_data_010425.csv"))

# Experiment 1: effect of the sample preservation on the DNA copies ####
# Data subsetting
Sub1.s<-raw.data[raw.data$Date_stdcurve%in% c('25/11/2020'),]
Sub1.1FO<-Sub1.s[Sub1.s$Trt %in% c('FO', 'FE'),] 

# Data arrangment and summary statistics
Sub1.samples.cat <- Sub1.1FO %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                  Trt = factor(Trt),
                                  conservation_t = factor(conservation_t),
                                  DNA_log10 = as.numeric(as.character(DNA_log10))) #categoric variables

Sub1.samples.cont <- Sub1.1FO %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                  Trt = factor(Trt),
                                  conservation_t = as.numeric(as.character(conservation_t)),
                                  DNA_log10 = as.numeric(as.character(DNA_log10)))

Sub1.stats<-Sub1.samples.cat%>% group_by(Trt, conservation_t) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count))

#Detections in negative controls
Sub1.neg <-raw.data[raw.data$Trt %in% c('SNC0', 'ENC0', 'SNC2','FNC0','FNC2','ENC0','ENC2','QNC2'),]

Sub1.1neg<-Sub1.neg %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                               Trt = factor(Trt), conservation_t = as.numeric(as.character(conservation_t)))

Sub1.neg.stats<-Sub1.neg%>% group_by(Trt, conservation_t, rep_bio, ID) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm=TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count)) #results in Table S2

# Data visualization with an histogram
g1 <- ggplot(Sub1.samples.cat, aes(x=DNA_copy)) + geom_histogram()
g1

g2 <- ggplot(Sub1.samples.cat, aes(x=DNA_log10)) + geom_histogram()
g2

# Generalized linear mixed model
mBO.cat<- lmer(DNA_log10 ~ conservation_t + (1|rep_bio:ID), data= Sub1.samples.cat)
summary(mBO.cat) #results in Table S3

# plot residuals 
plot(mBO.cat)

# residual distribution  
hist(resid(mBO.cat))

# Pairwise comparisons
T.c <- pairs(lsmeans(mBO.cat, ~ conservation_t)) 
T.c #Results in Table S4

# Figure: 
pred.mBO <- effect("conservation_t", mBO.cat, confidence.level=0.95) #predictions of the model
pred.mBO
resp.mBO <- summary(pred.mBO, type="response") #response of the model 
resp.mBO

graph.mBO <- data.frame(pred.mBO$x, pred.mBO$fit, se=pred.mBO$se, lowerCI=pred.mBO$lower, upperCI=pred.mBO$upper)
graph.mBO <- graph.mBO %>%  mutate(conservation_t = as.numeric(as.character(conservation_t)))

# Experiment 2: compare the properties of different filter types to restitute eDNA ####

#Data subsetting
Sub2 <- raw.data[raw.data$Trt %in% c("FNC1", "GF", "NY", "ST", 
                                     "PES", "SNC1", "ENC1", "QNC1"), ]
Sub2.samples <- raw.data[raw.data$Trt %in% c("GF", "NY", "ST", 
                                             "PES"),]

#Data arrangment and summary stats
Sub2.samples <- Sub2.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        DNA_log10 = as.numeric(as.character(DNA_log10)),
                                        Trt = factor(Trt))

Sub2.stats<-Sub2.samples%>% group_by(Trt) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count))

# Data visualization with an histogram
g3 <- ggplot(Sub2.samples, aes(x=DNA_copy)) + geom_histogram()
g3

g4 <- ggplot(Sub2.samples, aes(x=DNA_log10)) + geom_histogram()
g4

# Generalized linear mixed model
m2<- lmer(DNA_copy ~ Trt + (1|rep_bio:ID), data= Sub2.samples)
summary(m2) #Table S5

# plot residuals 
plot(m2)

# residual distribution  
hist(resid(m2))

# Pairwise comparisons
T.c <- pairs(lsmeans(m2, ~ Trt))
T.c #Table S6

# Experiment 3: compare the effect of filter preservation by several methods during mid-time exposure ####

#Data subsetting
Sub3<-raw.data[raw.data$Date_stdcurve=='25/11/2020',]

Sub3.samples<-Sub3[Sub3$Trt %in% c('SI','ET', 'SP','FI', 'FE'),] 
Sub3.samples <- Sub3.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        DNA_log10 = as.numeric(as.character(DNA_log10)),
                                        Trt = factor(Trt))

# Data visualization with an histogram
g5 <- ggplot(Sub3.samples, aes(x=DNA_copy)) + geom_histogram()
g5

g6 <- ggplot(Sub3.samples, aes(x=DNA_log10)) + geom_histogram()
g6

#Data arrangement and summary statistics
Sub3.samples.cat <- Sub3.samples %>% mutate(DNA_log10 = as.numeric(as.character(DNA_log10)),
                                              Trt = factor(Trt),
                                              conservation_t = factor(conservation_t))

Sub3.samples.cont <- Sub3.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                         Trt = factor(Trt),
                                         conservation_t = as.numeric(as.character(conservation_t)))

Sub3.stats<-Sub3.samples%>% group_by(Trt, conservation_t) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count))

Sub3.neg <-raw.data[raw.data$Trt %in% c('SNC0', 'ENC0', 'SNC2','FNC0','FNC2','ENC0','ENC2','QNC2'),]

# Generalized linear mixed models
m0.cat <- lmer(DNA_log10 ~ conservation_t+Trt + conservation_t:Trt+ (1|rep_bio:ID), data= Sub3.samples.cat) 
m1.cat <- lmer(DNA_log10 ~ conservation_t + (1+conservation_t|Trt) + (1|rep_bio:ID), data= Sub3.samples.cat)

# Comparisons of GLM
MuMIn::AICc(m0.cat)
MuMIn::AICc(m1.cat)

# GLM output
summary(m0.cat) #Table S7

# plot residuals 
plot(m0.cat)

# residual distribution  
hist(resid(m0.cat))

# Pairwise comparisons
T.c <- pairs(lsmeans(m0.cat, ~ Trt|conservation_t))
c.T <- pairs(lsmeans(m0.cat, ~ conservation_t|Trt))
rbind(T.c,c.T) #Table S8

# Experiment 4: compare the effect of commercial DNA extraction kits on the detection of eDNA ####

#Data subsetting
Sub4.samples<-raw.data[raw.data$Trt %in% c('BT', "PW", "BTZ", 'PWZ','BTT'),]
Sub4.neg<-raw.data[raw.data$Trt %in% c('SNC3', "FNC3", "ENC3", 'QNC3'),]

# Data arrangment and summary statistics
Sub4.samples <- Sub4.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        DNA_log10 = as.numeric(as.character(DNA_log10)),
                                        Trt = factor(Trt))

Sub4.samples.stats<-Sub4.samples%>% group_by(Trt, conservation_t) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count))

Sub4.neg <- Sub4.neg %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        Trt = factor(Trt))

Sub4.neg.stats<-Sub4.neg%>% group_by(Trt, conservation_t, rep_bio, ID) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),    se   = sd / sqrt(count))

# Data visualization with an histogram
g5 <- ggplot(Sub4.samples, aes(x=DNA_copy)) + geom_histogram()
g5

g6 <- ggplot(Sub4.samples, aes(x=DNA_log10)) + geom_histogram()
g6

# Generalized linear mixed model
m4<- lmer(DNA_copy ~ Trt + (1|rep_bio:ID), data= Sub4.samples)
summary(m4) #Table S9

# plot residuals 
plot(m4)

# residual distribution  
hist(resid(m4))

# Pairwise comparisons
T.c <- pairs(lsmeans(m4, ~ Trt))
T.c #Table S10

# Correction of PW and PWZ values for 100 uL elution volume
# Data subsetting
Sub4.samples<-raw.data[raw.data$Trt %in% c('BT', 'BTT', 'BTZ',"PW", 'PWZ', 'PW_c', 'PWZ_c'),]

# Data arrangment and summary statistics
Sub4.samples <- Sub4.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        DNA_log10 = as.numeric(as.character(DNA_log10)),
                                        Trt = factor(Trt))

# Data visualization with an histogram
g7 <- ggplot(Sub4.samples, aes(x=DNA_copy)) + geom_histogram()
g7

g8<- ggplot(Sub4.samples, aes(x=DNA_log10)) + geom_histogram()
g8

# Generalized linear mixed model
m4<- lmer(DNA_copy ~ Trt + (1|rep_bio:ID), data= Sub4.samples)
summary(m4)

# plot residuals 
plot(m4)

# residual distribution  
hist(resid(m4))

# Pairwise comparisons
T.c <- pairs(lsmeans(m4, ~ Trt))
T.c
