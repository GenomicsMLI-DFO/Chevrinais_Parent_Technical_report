# code used for statistics in Chevrinais and Parent, Methods in Ecology and Evolution

# packages used
library(car)
library(dplyr)
library(bbmle)
library(MuMIn)
library(lmerTest)
library(lsmeans)
library(rstatix)
library(ARTool)
library(ggplot2)

# Updload raw data
raw.data <- read.csv(file.choose(), sep=';') #File Chevrinais_2022_raw_data.csv

# Step 1: sample preservation 

# Data subsetting
Sub1<-raw.data[raw.data$qPCR_mix=='GMM',]
Sub1<-Sub1[Sub1$Date_stdcurve=='25/11/2020',]
Sub1.1FO<-Sub1[Sub1$Trt %in% c('FO', 'FE'),] 

# Data arrangment and summary statistics
Sub1.samples.cat <- Sub1.1FO %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                  Trt = factor(Trt),
                                  conservation_t = factor(conservation_t))

Sub1.samples.cont <- Sub1.1FO %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                       Trt = factor(Trt),
                                       conservation_t = as.numeric(as.character(conservation_t)))

# Data visualization with an histogram
qplot(Sub1.1FO$DNA_copy, geom="histogram", bins=30)
qplot(Sub1.1FO$DNA_log10, geom="histogram",bins=30)

# Generalized linear mixed model
mBO.cat<- lmer(DNA_log10 ~ conservation_t + (1|rep_bio:ID), data= Sub1.samples.cat)
summary(mBO.cat)

mBO.cont <- lmer(DNA_log10 ~ conservation_t + (1|rep_bio:ID), data= Sub1.samples.cont)
summary(mBO.cont)

# plot residuals 
plot(mBO.cat)
plot(mBO.cont)

# residual distribution  
hist(resid(mBO.cat))
hist(resid(mBO.cont))

# Pairwise comparisons
T.c <- pairs(lsmeans(mBO.cat, ~ conservation_t))
T.c

# Step 2: compare the properties of different filter types to restitute eDNA

#Data subsetting
Sub2 <- raw.data[raw.data$Trt %in% c("FNC1", "GF", "NY", "ST", 
                                     "PES", "SNC1", "ENC1", "QNC1"), ]
Sub2.samples <- raw.data[raw.data$Trt %in% c("GF", "NY", "ST", 
                                             "PES"),]

#Data arrangment and summary stats
Sub2.samples <- Sub2.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        Trt = factor(Trt))
# Data visualization with an histogram
qplot(Sub2.samples$DNA_copy, geom="histogram", bins=30)
qplot(Sub2.samples$DNA_log10, geom="histogram",bins=30)

# Generalized linear mixed model
m2<- lmer(DNA_copy ~ Trt + (1|rep_bio:ID), data= Sub2.samples)
summary(m2)

# plot residuals 
plot(m2)

# residual distribution  
hist(resid(m2))

# Pairwise comparisons
T.c <- pairs(lsmeans(m2, ~ Trt))
T.c

# Step 3: compare the effect of filter preservation by several methods during mid-time exposure

#Data subsetting
Sub3<-raw.data[raw.data$qPCR_mix=='GMM',]
Sub3<-Sub3[Sub3$Date_stdcurve=='25/11/2020',]

Sub3.samples<-Sub3[Sub3$Trt %in% c('SI','ET', 'SP','FI', 'FE'),] 

# Data visualization with an histogram
qplot(Sub3.samples$DNA_copy, geom="histogram", bins=30)
qplot(Sub3.samples$DNA_log10, geom="histogram",bins=30)

#Data arrangement and summary statistics
Sub3.samples.cat <- Sub3.samples %>% mutate(DNA_log10 = as.numeric(as.character(DNA_log10)),
                                              Trt = factor(Trt),
                                              conservation_t = factor(conservation_t))

Sub3.samples.cont <- Sub3.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                         Trt = factor(Trt),
                                         conservation_t = as.numeric(as.character(conservation_t)))

# Generalized linear mixed models
m0.cat <- lmer(DNA_log10 ~ conservation_t+Trt + conservation_t:Trt+ (1|rep_bio:ID), data= Sub3.samples.cat) 
m1.cat <- lmer(DNA_log10 ~ conservation_t + (1+conservation_t|Trt) + (1|rep_bio:ID), data= Sub3.samples.cat)

m0.cont <- lmer(DNA_log10 ~ conservation_t+Trt + conservation_t:Trt+ (1|rep_bio:ID), data= Sub3.samples.cont) 
m1.cont <- lmer(DNA_log10 ~ conservation_t + (1+conservation_t|Trt) + (1|rep_bio:ID), data= Sub3.samples.cont)

# Comparisons of GLM
MuMIn::AICc(m0.cat)
MuMIn::AICc(m1.cat)
MuMIn::AICc(m0.cont)
MuMIn::AICc(m1.cont)

# GLM output
summary(m0.cat)
summary(m0.cont)

# plot residuals 
plot(m0.cont)
plot(m0.cat)

# residual distribution  
hist(resid(m0.cont))
hist(resid(m0.cat))

# Pairwise comparisons
T.c <- pairs(lsmeans(m0.cat, ~ Trt|conservation_t))
c.T <- pairs(lsmeans(m0.cat, ~ conservation_t|Trt))
rbind(T.c,c.T)

# Step 4: compare the effect of commercial DNA extraction kits on the detection of eDNA

#Data subsetting
Sub4.samples<-raw.data[raw.data$Trt %in% c('BT', "PW", "BTZ", 'PWZ','BTT'),]
Sub4.samples<-Sub4.samples[Sub4.samples$qPCR_mix %in% c('GMM'),]

# Data arrangment and summary statistics
Sub4.samples <- Sub4.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        Trt = factor(Trt))

# Data visualization with an histogram
qplot(Sub4.samples$DNA_copy, geom="histogram", bins=30)
qplot(Sub4.samples$DNA_log10, geom="histogram",bins=30)

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

# Correction of PW and PWZ values for 100 uL elution volume
# Data subsetting
Sub4.samples<-raw.data[raw.data$Trt %in% c('BT', 'BTT', 'BTZ',"PW", 'PWZ', 'PW_c', 'PWZ_c'),]
Sub4.samples<-Sub4.samples[Sub4.samples$qPCR_mix %in% c('GMM'),]

# Data arrangment and summary statistics
Sub4.samples <- Sub4.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        Trt = factor(Trt))

# Data visualization with an histogram
qplot(Sub4.samples$DNA_copy, geom="histogram", bins=30)
qplot(Sub4.samples$DNA_log10, geom="histogram",bins=30)

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

# Step 5: compare the effect of the qPCR master mix on the detection of eDNA

# Data subsetting
Sub5.samples<-raw.data[raw.data$Trt %in% c('FO', "ET", "SI", 'FI','SP'),]
Sub5.samples<-Sub5.samples[Sub5.samples$conservation_t %in% c('30'),]
Sub5.samples<-Sub5.samples[Sub5.samples$Date_stdcurve %in% c('19/10/2020','05/10/2020','09/11/2020'),]

# Data arrangment and summary statistics
Sub5.samples.model <- Sub5.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                              qPCR_mix = factor(qPCR_mix),
                                             Trt = factor(Trt))

Sub5.samples.graph <-Sub5.samples.model%>% group_by(Trt, qPCR_mix) %>% dplyr::summarise(count = n(), mean = mean(DNA_copy, na.rm = TRUE), sd=sd(DNA_copy, na.rm=TRUE),
                                                                                                               se   = sd / sqrt(count))

# Data visualization with an histogram
qplot(Sub5.samples$DNA_copy, geom="histogram", bins=30)
qplot(Sub5.samples$DNA_log10, geom="histogram",bins=30)

# ET treatment
# Generalized linear mixed model
m5.ET<- lmer(DNA_copy ~ qPCR_mix + (1|rep_bio:ID), data= Sub5.samples[Sub5.samples$Trt %in% c('ET'),])
summary(m5)

# plot residuals 
plot(m5.ET)

# residual distribution  
hist(resid(m5.ET))

# Pairwise comparisons
T.c <- pairs(lsmeans(m5.ET, ~ qPCR_mix))
T.c

# all other treatments
Sub5.samples<-Sub5.samples[Sub5.samples$Trt %in% c('FO', "SI", 'FI','SP'),]

# Generalized linear mixed model
m5<- lmer(DNA_copy ~ qPCR_mix+Trt + Trt:qPCR_mix + (1|rep_bio:ID), data= Sub5.samples)
summary(m5)

# plot residuals 
plot(m5)

# residual distribution  
hist(resid(m5))

# Pairwise comparisons
T.c <- pairs(lsmeans(m5, ~ Trt|qPCR_mix))
c.T <- pairs(lsmeans(m5, ~ qPCR_mix|Trt))
rbind(T.c,c.T)

# Extraction kits comparisons in environmental samples

#Data subsetting
Sub6<-raw.data[raw.data$Trt %in% c('BTTZ', 'BTT','PW'),]
Sub6<-Sub6[Sub6$Site %in% c('C', 'D','B', 'A'),]
SiteC<-Sub6[Sub6$Site %in% c('C'),]
SiteD<-Sub6[Sub6$Site %in% c('D'),]
SiteB<-Sub6[Sub6$Site %in% c('B'),]
SiteA<-Sub6[Sub6$Site %in% c('A'),]

#Data arrangment and summary statistics
Sub6.samples <- Sub6 %>% mutate(Cq = as.numeric(as.character(Cq)),
                                Site = factor(Site),
                                Trt = factor(Trt))


Sub6.stats<-Sub6.samples %>% group_by(Site, Trt) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(delta_cq, na.rm = TRUE),
    sd = sd(delta_cq, na.rm = TRUE),se   = sd / sqrt(count)
  )

# Two way ANOVA
res.aov<-aov(delta_cq~Trt*Site, data = Sub6) 
summary(res.aov)

#ANOVA for site A
res.aov<-aov(delta_cq~Trt, data = SiteA)
summary(res.aov)
#Tukey pairwise comparisons - site A
TukeyHSD(res.aov)
# Check homogeneity of residual variances - site A
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)
# check for normality of residuals - site A
plot(res.aov,1) 
leveneTest(delta_cq~Trt*Site, data = SiteA)

#ANOVA for site B
res.aov<-aov(delta_cq~Trt, data = SiteB)
summary(res.aov)
#Tukey pairwise comparisons - site B
TukeyHSD(res.aov)
# Check homogeneity of residual variances - site B
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)
# check for normality of residuals - site B
plot(res.aov,1) 
leveneTest(delta_cq~Trt*Site, data = SiteB)

#ANOVA for site C
res.aov<-aov(delta_cq~Trt, data = SiteC)
summary(res.aov)
#Tukey pairwise comparisons - site C
TukeyHSD(res.aov)
# Check homogeneity of residual variances - site C
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)
# check for normality of residuals - site C
plot(res.aov,1) 
leveneTest(delta_cq~Trt*Site, data = SiteC)

#ANOVA for site D
res.aov<-aov(delta_cq~Trt, data = SiteD)
summary(res.aov)
#Tukey pairwise comparisons - site D
TukeyHSD(res.aov)
# Check homogeneity of residual variances - site D
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)
# check for normality of residuals - site D
plot(res.aov,1) 
leveneTest(delta_cq~Trt*Site, data = SiteD)

# qPCR mix comparisons in environmental samples

#Data subsetting
Sub7<-raw.data[raw.data$Site %in% c('A', 'B','C', 'D'),]
Sub7<-Sub7[Sub7$qPCR_mix %in% c('GMM', 'EMM','TPM'),]

#Data arrangement and summary statistics
Sub7.stats<-group_by(Sub7,qPCR_mix) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(delta_cq, na.rm = TRUE),
    sd = sd(delta_cq, na.rm = TRUE), se = sd/sqrt(count)
  )

# two way ANOVA
res.aov2<-aov(delta_cq~qPCR_mix*Site, data = Sub7)
summary(res.aov2)

# Tukey multiple pairwise comparisons
TukeyHSD(res.aov2)

# Check homogeneity of residual variances
plot(res.aov2, 2)
aov_residuals <- residuals(object = res.aov2)
shapiro.test(x = aov_residuals)
# check for normality of residuals
plot(res.aov2,1) 
leveneTest(delta_cq ~ qPCR_mix*Site, data = Sub7)
