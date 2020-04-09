#load the datasets
yv=read.csv("Exp2_updated.csv", header= TRUE)
backcross=read.csv("exp2_bc_setup.csv", header=T)
yv=na.omit(yv)
colnames(yv)

#The 24 hour transfers are aggregated into the same with experiment 3 due to insufficient sample sizes.
yv$new_day=ifelse(yv$Day=="A","A",
                  ifelse(yv$Day..letter.of.vial.=="B","A",
                         ifelse(yv$Day..letter.of.vial.=="C","B",
                                ifelse(yv$Day..letter.of.vial.=="D","B",
                                       ifelse(yv$Day..letter.of.vial.=="E","C",
                                              ifelse(yv$Day..letter.of.vial.=="F","C",
                                                     ifelse(yv$Day..letter.of.vial.=="G","D",
                                                            ifelse(yv$Day..letter.of.vial.=="H","D",
                                                                   ifelse(yv$Day..letter.of.vial.=="I","D",
                                                                          ifelse(yv$Day..letter.of.vial.=="J","E",
                                                                                 ifelse(yv$Day..letter.of.vial.=="K","E",
                                                                                        ifelse(yv$Day..letter.of.vial.=="L","E",
                                                                                               ifelse(yv$Day..letter.of.vial.=="M","F",
                                                                                                      ifelse(yv$Day..letter.of.vial.=="N","F",
                                                                                                             ifelse(yv$Day..letter.of.vial.=="o","F",NA)))))))))))))))
#backcross$Treatment=as.numeric(as.character(backcross$Treatment)) #to make sure R reads it as a character rather than a number.

yv$Wildtype=as.numeric(as.character(yv$Wildtype))
yv$Vellow.vermillion=as.numeric(as.character(yv$Vellow.vermillion))
yv$Yellow.only=as.numeric(as.character(yv$Yellow.only))
yv$Vermillion.only=as.numeric(as.character(yv$Vermillion.only))
yv$Day..letter.of.vial.=as.character(yv$Day..letter.of.vial.)
yv=na.omit(yv)

#necessary packages

library(ggplot2)
library(ggthemes)
library(devtools)
library(gridExtra)
library(gtable)
library(grid)
library(lsmeans)
library(lme4)
library(lmerTest)
library(doBy)
library(pbkrtest)
library(tidyr)

#CO groups were defined

yv$SCO1=yv$Yellow.only
yv$SCO2=yv$Vermillion.only
yv$NCO1=yv$Wildtype
yv$NCO2=yv$Vellow.vermillion

sco_count=sum(yv$SCO1, na.rm = TRUE)+sum(yv$SCO2, na.rm = TRUE)
nco_count=sum(yv$NCO1+yv$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#rate between yellow and vermillion

(sum(yv$SCO1, na.rm = TRUE)+sum(yv$SCO2, na.rm = TRUE))/num_samples


#merge with treatment data

yv_merged <- merge(yv, backcross, by.x = "Vial..", by.y = "ï..Vial.Number", all=T)
yv_merged = na.omit(yv_merged)
yv_merged=na.omit(yv_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+new_day+Treatment,data=yv_merged, FUN=sum,na.rm=T)

#add in a column for total offspring

dataset$total_offspring=dataset$SCO1.sum + dataset$SCO2.sum + dataset$NCO1.sum + dataset$NCO2.sum

num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
  #f1_vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
  #print(f1_vial)
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$Vial..)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms

dataset$fecundity=dataset$total_offspring/dataset$Num_moms

write.csv(dataset,file="Experiment2_data.csv")

pdf("Experiment2_fecundity.pdf")

#Fecundity figure
Fecund_figure=ggplot(aes(y=fecundity,x=new_day, col=as.factor(Treatment),label=Num_moms),data=dataset)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_update(text = element_text(size=20))
Fecund_figure= Fecund_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

dev.off()

#Statistics of the vy interval
#Let's start with cleanedup data for Exp 3:
dataset2=read.csv("Experiment2_data.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)

#get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$new_day,mean)

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*new_day,data=dataset2,family=quasipoisson)
summary(fit)
#Anova
anova(fit,test="Chisq")
#post-hoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="new_day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr

#Now do the same, but with recombination rate

#sum of crossovers in intervals 1
dataset2$num_CO=dataset2$SCO1.sum+dataset2$SCO2.sum
#sum of non-crossovers in intervals 1
dataset2$num_NCO=dataset2$NCO1.sum+dataset2$NCO2.sum
#total recombination rate
dataset2$rec_rate_total=dataset2$num_CO/dataset2$total_offspring

#get mean recombination
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_total,dataset2$new_day,mean)
fit_recrate=glm(rec_rate_total~Treatment*new_day,data=dataset2,family=quasipoisson)
summary(fit_recrate)
#Anova
anova(fit_recrate,test="Chisq")
#posthoc
fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="new_day", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec
#logistic regression, similar to a t-test for count data
dataset2$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset2$F1.Vial))
dataset2$F1.Vial=as.numeric(dataset2$F1.Vial)
fit2=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment*new_day,data=dataset2,
         family=binomial(link="logit"))
summary(fit2)
coefs=coef(fit2)
coefs

#repeat, but only for day a
fit2a=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="A"))
summary(fit2a)
exp(coef(fit2a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit2b=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="B"))
summary(fit2b)
exp(coef(fit2b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b

#repeat, but only for day c
fit2c=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="C"))
summary(fit2c)
exp(coef(fit2c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c

#repeat, but only for day d
fit2d=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="D"))
summary(fit2d)
exp(coef(fit2d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d

#repeat, but only for day e
fit2e=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="E"))
summary(fit2e)
exp(coef(fit2e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day e

#repeat, but only for day f
fit2f=glm(cbind(num_CO,num_NCO)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="F"))
summary(fit2f)
exp(coef(fit2f)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day f

#Summarization of the data

length(unique(dataset2$F1.Vial[dataset2$Treatment=="18"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="24"]))

median(dataset2$Num_moms[dataset2$Treatment=="18"])
median(dataset2$Num_moms[dataset2$Treatment=="24"])


sum(dataset2$num_CO[dataset2$Treatment=="18"])+sum(dataset2$num_NCO[dataset2$Treatment=="18"])

sum(dataset2$num_NCO[dataset2$Treatment=="24"])+sum(dataset2$num_CO[dataset2$Treatment=="24"])


#Odd ratios are extracted and entered to excel and plotted in r
exp1=read.csv(file="exp2.csv", header=T)
exp1$time=as.character(exp2$time)
pdf("odds_exp2.pdf")
# Basic line plot with points
exp2_odds=ggplot(data=exp2, aes(x=time, y=odds, group=1))+
  geom_line(size=2)+
  theme_bw()+
  theme_update(text = element_text(size=25))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point()+
  geom_errorbar(aes(ymin=odds-sd, ymax=odds+sd), width=.2, position=position_dodge(0.05))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(0.5,2.5)+
  ylab("Odds")+
  scale_x_discrete(name="Days", labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+
  ggtitle("Odd Ratios for experiment 2")
exp2_odds
dev.off()