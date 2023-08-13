rm(list=ls())
library(tidyverse); library(RCurl); library(ez); library(lme4); library(car); 
library(lmerTest);library(ggeffects); library(ggplot2);library(stringr);library(devtools);
library(ggiraphExtra);library(RColorBrewer);library(MuMIn)

DATA <- read.csv(str_c('fnfa246_6-13.csv'),
                 stringsAsFactors = TRUE)


data_TIME <- DATA
data_TIME$sex <- factor(data_TIME$sex)
data_TIME$centre <- factor(data_TIME$centre)

# Centering time on the first point and converting to years:
data_TIME$age.0 <- data_TIME$age
data_TIME$age.0_sq<-data_TIME$age.0^2


# Random slope random intercepts model ---- 
raneff_01<-lmer(get(varname)~
                # Fixed-effects
                1+sex+centre+tbv+age.0+
                (-1+age.0|id)+(1|id),
                data=data_TIME,REML=FALSE,
                control=lmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=5e6)))
Results1<-summary(raneff_01)[["coefficients"]]
summary(raneff_01)
confint(raneff_01)   # IC



# # Random quadratic slopes and intercepts model ---- 
raneff_02<-lmer(get(varname)~
                  # Fixed-effects
                  1+sex+centre+tbv+age.0+age.0_sq+
                  (-1+age.c_sq|id)+(-1+age.c|id)+(1|id), data=data_TIME, REML=FALSE,
                control=lmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=5e5)))
summary(raneff_02)

