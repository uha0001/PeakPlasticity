genotyping=read.csv("SNP_genotyping.csv",header = TRUE)
genotyping$Treatment=as.character(genotyping$Treatment)
genotyping$F1.Vial=as.numeric(genotyping$F1.Vial)
library(dplyr)
library(reshape2)
library(ggplot2)
library(emmeans)
library(lme4)
library(lmerTest)
library(doBy)

#get mean recombination
tapply(genotyping$rec_rate_total,genotyping$Treatment,mean)
tapply(genotyping$rec_rate_total,genotyping$Day,mean)
mean(genotyping_21$rec_rate_i3)
mean(genotyping_21$rec_rate_i2)
mean(genotyping_21$rec_rate_i1)
genotyping_21=subset(genotyping, genotyping$Treatment=="23")

fit_recrate=glm(rec_rate_total~Treatment*Day,data=genotyping,family=quasipoisson)
summary(fit_recrate)
anova(fit_recrate,test="Chisq")

fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="Day", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec

#i1 REGION
#logistic regression, similar to a t-test for count data
fit1=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit1)
coefs=coef(fit1)
coefs

anova(fit1,test="Chisq")
fit_contrast_i1 <- emmeans::emmeans(fit1, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i1 <- contrast(fit_contrast_i1, method="trt.vs.ctrl")

pheno_contr_i1 <- as.data.frame(summary(fit_contr_i1))
pheno_contr_i1
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day a
fit2a=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="a"))
summary(fit2a)
exp(coef(fit2a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit2b=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="b"))
summary(fit2b)
exp(coef(fit2b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b


#repeat, but only for day c
fit2c=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="c"))

summary(fit2c)
exp(coef(fit2c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day d
fit2d=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="d"))
summary(fit2d)
exp(coef(fit2d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit2e=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="e"))
summary(fit2e)
exp(coef(fit2e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e

#i2 region
#logistic regression, similar to a t-test for count data
fit3=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit3)
coefs=coef(fit3)
coefs

anova(fit3,test="Chisq")
fit_contrast_i2 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i2 <- contrast(fit_contrast_i2, method="trt.vs.ctrl")

pheno_contr_i2 <- as.data.frame(summary(fit_contr_i2))
pheno_contr_i2
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day a
fit4a=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="a"))
summary(fit4a)
exp(coef(fit4a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit4b=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="b"))
summary(fit4b)
exp(coef(fit4b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b


#repeat, but only for day c
fit4c=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="c"))
summary(fit4c)
exp(coef(fit4c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day d
fit4d=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="d"))
summary(fit4d)
exp(coef(fit4d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit4e=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="e"))
summary(fit4e)
exp(coef(fit4e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e

#i3 region
#logistic regression, similar to a t-test for count data
fit5=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit5)
coefs=coef(fit5)
coefs

anova(fit5,test="Chisq")

fit_contrast_i3 <- emmeans::emmeans(fit5, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i3 <- contrast(fit_contrast_i3, method="trt.vs.ctrl")

pheno_contr_i3 <- as.data.frame(summary(fit_contr_i3))
pheno_contr_i3

#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day a
fit6a=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="a"))
summary(fit6a)
exp(coef(fit6a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit6b=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="b"))
summary(fit6b)
exp(coef(fit6b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b


#repeat, but only for day c
fit6c=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="c"))
summary(fit6c)
exp(coef(fit6c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day d
fit6d=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="d"))
summary(fit6d)
exp(coef(fit6d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit6e=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="e"))
summary(fit6e)
exp(coef(fit6e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e

#i4 region
#logistic regression, similar to a t-test for count data
fit7=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit7)
coefs=coef(fit7)
coefs

anova(fit7,test="Chisq")
fit_contrast_i4 <- emmeans::emmeans(fit7, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i4 <- contrast(fit_contrast_i4, method="trt.vs.ctrl")

pheno_contr_i4 <- as.data.frame(summary(fit_contr_i4))
pheno_contr_i4
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day a
fit8a=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="a"))
summary(fit8a)
exp(coef(fit8a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit8b=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="b"))
summary(fit8b)
exp(coef(fit8b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b


#repeat, but only for day c
fit8c=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="c"))
summary(fit8c)
exp(coef(fit8c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day d
fit8d=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="d"))
summary(fit8d)
exp(coef(fit8d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit8e=glm(cbind(num_CO_4,num_NCO_4)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="e"))
summary(fit8e)
exp(coef(fit8e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e

#i5 region
#logistic regression, similar to a t-test for count data
fit9=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment*Day,data=genotyping,
         family=binomial(link="logit"))
summary(fit9)
coefs=coef(fit9)
coefs

anova(fit9,test="Chisq")
fit_contrast_i5 <- emmeans::emmeans(fit9, "Treatment", by="Day", mode="kenward-roger")
fit_contr_i5 <- contrast(fit_contrast_i5, method="trt.vs.ctrl")

pheno_contr_i5 <- as.data.frame(summary(fit_contr_i5))
pheno_contr_i5
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day a
fit10a=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="a"))
summary(fit10a)
exp(coef(fit10a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit10b=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="b"))
summary(fit10b)
exp(coef(fit10b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b


#repeat, but only for day c
fit10c=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="c"))
summary(fit10c)
exp(coef(fit10c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day d
fit10d=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="d"))
summary(fit10d)
exp(coef(fit10d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit10e=glm(cbind(num_CO_5,num_NCO_5)~(1|F1.Vial)+Treatment,data=genotyping,
          family=binomial(link="logit"),subset=c(Day=="e"))
summary(fit10e)
exp(coef(fit10e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e


#Odds ratios are entered excel manually

genotyping=read.csv(file="genotyping.csv", header=T)
genotyping$time=as.character(genotyping$time)
pdf("SNP.pdf")
genotyping_odds=ggplot(data=genotyping, aes(x=time, y=i1, group=1)) +
  geom_line(color= "black", size=2)+
  geom_errorbar(aes(ymin=i1-sd1, ymax=i1+sd1), width=.2, position=position_dodge(0.05))+
  geom_line(data=genotyping, aes(x=time, y=i2, group=1), color= "black",linetype="dashed",size=2)+
  geom_errorbar(aes(ymin=i2-sd2, ymax=i2+sd2), width=.2, position=position_dodge(0.05))+
  geom_line(data=genotyping, aes(x=time, y=i3, group=1), color= "black",linetype="dotted", size=2)+
  geom_errorbar(aes(ymin=i3-sd3, ymax=i3+sd3), width=.2, position=position_dodge(0.05))+
  geom_line(data=genotyping, aes(x=time, y=i4, group=1), color= "red", size=2)+
  geom_errorbar(aes(ymin=i4-sd4, ymax=i4+sd4), width=.2, position=position_dodge(0.05))+
  geom_line(data=genotyping, aes(x=time, y=i5, group=1), color= "red",linetype="dashed", size=2)+
  geom_errorbar(aes(ymin=i5-sd5, ymax=i5+sd5), width=.2, position=position_dodge(0.05))+
  #theme_bw()+
  theme_update(text = element_text(size=40))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(-1,8.0)+
  ylab("odd ratios")+
  scale_x_discrete(name="Days", labels=c("1-2","3-4","5-6","7-8","9-10"))+
  ggtitle("Odd Ratios for genotyping results")
genotyping_odds
dev.off()