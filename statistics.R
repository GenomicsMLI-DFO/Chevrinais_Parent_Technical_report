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

# Updload raw data
raw.data <- read.csv(file.choose(), sep=';') #File Chevrinais_2022_raw_data.csv

# Step 1: sample preservation 

# Data subsetting
Sub1<-raw.data[raw.data$qPCR_mix=='GMM',]
Sub1<-Sub1[Sub1$Date_stdcurve=='25/11/2020',]
Sub1.1BO<-Sub1[Sub1$Trt %in% c('BO'),]
Sub1.1neg<-Sub1[Sub1$Trt %in% c('ENC2', 'ENC3', 'CONTROL'),]

# Data arrangment and summary statistics
Sub1.1neg<-Sub1.1neg %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                Trt = factor(Trt), conservation_t = as.numeric(as.character(conservation_t)))

data.model <- Sub1.1BO %>% mutate(DN_copy = as.numeric(as.character(DNA_copy)),
                                  Trt = factor(Trt),
                                  conservation_t = as.numeric(as.character(conservation_t)))

# Generalized linear model
mBO<- lmer(DNA_log10 ~ conservation_t + (1|rep_bio:ID), data= data.model)
summary(mBO)

# plot residuals 
plot(mBO)

# residual distribution  
hist(resid(mBO))

# Step 2: compare the properties of different filter types to restitute eDNA

#Data subsetting
Sub2 <- raw.data[raw.data$Trt %in% c("FNC1", "GF", "NY", "ST", 
                                   "PES", "SNC1", "ENC1", "QNC1"), ]
Sub2.samples <- raw.data[raw.data$Trt %in% c("GF", "NY", "ST", 
                                          "PES"),]
Sub2.neg<-Sub2[Sub2$Trt%in% c("SNC1"),]

#Data arrangment and summary stats
Sub2.neg<-Sub2.neg %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                              Trt = factor(Trt))

Sub2.samples <- Sub2samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                       Trt = factor(Trt))

# ANOVA to test differences between treatments
res.aov<-aov(DNA_copy~Trt, data = Sub2.samples)
summary(res.aov)

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Check homogeneity of residual variances
plot(res.aov,1)
leveneTest(DNA_copy ~ Trt, data = Sub2.samples)

# check for normality of residuals
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

# Step 3: compare the effect of filter preservation by several methods during mid-time exposure

#Data subsetting
Sub3<-raw.data[raw.data$qPCR_mix=='GMM',]
Sub3<-Sub3[Sub3$Date_stdcurve=='25/11/2020',]

Sub3.samples<-Sub3[Sub3$Trt %in% c('SI','ET', 'SP','FI'),]
Sub3.neg<-Sub3[Sub3$Trt %in% c('SNC2', 'SNC3', 'CONTROL'),]

#Data arrangement and summary statistics
Sub3.neglog<-Sub3.neg %>% mutate(DNA_log10 = as.numeric(as.character(DNA_log10)),
                                 Trt = factor(Trt),
                                 conservation_t = as.numeric(as.character(conservation_t)))

Sub3.samples.model <- Sub3.samples %>% mutate(DNA_log10 = as.numeric(as.character(DNA_log10)),
                                              Trt = factor(Trt),
                                              conservation_t = factor(conservation_t)
                                              
)
            
# Generalized linear models
m0 <- lmer(DNA_log10 ~ conservation_t+Trt + conservation_t:Trt+ (1|rep_bio:ID), data= Sub3.samples.model) 
m1 <- lmer(DNA_log10 ~ conservation_t + (1+conservation_t|Trt) + (1|rep_bio:ID), data= Sub3.samples.model)

# Comparisons of GLM
MuMIn::AICc(m0)
MuMIn::AICc(m1)

# GLM output
summary(m0)

# plot residuals 
plot(mBO)

# residual distribution  
hist(resid(mBO))

# Pairwise comparisons
T.c <- pairs(lsmeans(m0, ~ Trt|conservation_t))
c.T <- pairs(lsmeans(m0, ~ conservation_t|Trt))
rbind(T.c,c.T)

# Step 4: compare the effect of commercial DNA extraction kits on the detection of eDNA

#Data subsetting
Sub4.samples<-raw.data[raw.data$Trt %in% c('BT', "PW", "BTZ", 'PWZ','BTT'),]
Sub4.samples<-Sub4.samples[Sub4.samples$qPCR_mix %in% c('GMM'),]
Sub4.neg<-raw.data[raw.data$Trt %in% c('SNC3'),]


# Data arrangment and summary statistics
Sub4.neg<-Sub4.neg %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                              Trt = factor(Trt))

Sub4.samples <- Sub4.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                        Trt = factor(Trt))

# ANOVA to test differences between treatments
res.aov<-aov(DNA_copy~Trt, data = Sub4.samples)
summary(res.aov)

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Check homogeneity of residual variances
plot(res.aov,1) 
leveneTest(DNA_copy ~ Trt, data = Sub4.samples)

# check for normality of residuals
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

# Step 5: compare the effect of the qPCR master mix on the detection of eDNA

# Data subsetting
Sub5.samples<-raw.data[raw.data$Trt %in% c('BO', "ET", "SI", 'FI','SP'),]
Sub5.samples<-Sub5.samples[Sub5.samples$conservation_t %in% c('30'),]
Sub5.samples<-Sub5.samples[Sub5.samples$Date_stdcurve %in% c('19/10/2020','05/10/2020','09/11/2020'),]

# Data arrangment and summary statistics
Sub5.samples.model <- Sub5.samples %>% mutate(DNA_copy = as.numeric(as.character(DNA_copy)),
                                              qPCR_mix = factor(qPCR_mix),
                                              Trt = factor(Trt))

# ET treatment
# ANOVA for comparisons of qPCR master mixes
res.aov<-aov(DNA_log10~qPCR_mix, data = Sub5.samples[Sub5.samples$Trt %in% c('ET'),]) 
summary(res.aov)

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)

# Check homogeneity of residual variances
plot(res.aov,1) 
leveneTest(DNA_log10 ~ qPCR_mix, data = Sub5.samples[Sub5.samples$Trt %in% c('ET'),])

# check for normality of residuals
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

# Non parametric alternative to ANOVA
kruskal.test(DNA_copy ~ qPCR_mix, data = Sub5.samples[Sub5.samples$Trt %in% c('ET'),])
 
# Pairwise comparisons
pairwise.wilcox.test(Sub5.samples[Sub5.samples$Trt %in% c("ET"),]$DNA_copy, Sub5.samples[Sub5.samples$Trt %in% c("ET"),]$qPCR_mix,
                     p.adjust.method = "BH")

# all other treatments
Sub5.samples<-Sub5.samples[Sub5.samples$Trt %in% c('BO', "SI", 'FI','SP'),]

# Two-way ANOVA 
res.aov2<-aov(DNA_log10~qPCR_mix*Trt, data = Sub5.samples[Sub5.samples$Trt %in% c("SI","SP","BO","FI"),])
summary(res.aov2)

# Check homogeneity of residual variances
plot(res.aov2,1)
leveneTest(DNA_log10 ~ qPCR_mix*Trt, data = Sub5.samples[Sub5.samples$Trt %in% c("SI","SP","BO","FI"),]) 

# check for normality of residuals
plot(res.aov2, 2)
aov_residuals <- residuals(object = res.aov2 )
shapiro.test(x = aov_residuals)

# Non-parametric alternative to two-way ANOVA
Sub5.samples <- Sub5.samples[Sub5.samples$Trt %in% c("SI","SP","BO","FI"),]
Sub5.samples$DNA_log10= as.numeric(as.factor(Sub5.samples$DNA_log10)) 
Sub5.samples$Trt = factor(Sub5.samples$Trt) 
Sub5.samples$qPCR_mix = factor(Sub5.samples$qPCR_mix) 
m = art(DNA_log10 ~ Trt*qPCR_mix, data=Sub5.samples)
anova(m)

# Pairwise comparisons
art.con(m, ~ Trt*qPCR_mix, adjust="holm") %>% 
  summary() %>% 
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints= c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))

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

