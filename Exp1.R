#load the datasets
yv=read.csv("Exp1_rawdata.csv", header= TRUE)
backcross=read.csv("Exp1_bc_setup.csv", header=T)
yv=na.omit(yv)
backcross$Treatment=as.numeric(as.character(backcross$Treatment))

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

yv$wildtype=as.numeric(as.character(yv$wildtype))
yv$yellow.vermillion=as.numeric(as.character(yv$yellow.vermillion))
yv$yellow.only=as.numeric(as.character(yv$yellow.only))
yv$vermillion.only=as.numeric(as.character(yv$vermillion.only))
yv=na.omit(yv)

yv$SCO1=yv$yellow.only
yv$SCO2=yv$vermillion.only
yv$NCO1=yv$wildtype
yv$NCO2=yv$yellow.vermillion


sco_count=sum(yv$SCO1, na.rm = TRUE)
nco_count=sum(yv$NCO1, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#rate between yellow and vermillion

(sum(yv$SCO1, na.rm = TRUE))/num_samples

#merge with treatment data
yv_merged <- merge(yv, backcross, by.x = "Vial.number", by.y = "Vial.number", all=T)
yv_merged = na.omit(yv_merged)
colnames(yv_merged)

dataset=summaryBy(SCO1+NCO1~F1.Vial+Day..letter.of.vial.+Treatment,data=yv_merged, FUN=sum,na.rm=T)

#add in a column for total offspring

dataset$total_offspring=dataset$SCO1.sum + dataset$NCO1.sum

num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
    F1.vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
    #print(f1_vial)
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$Vial.number)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms
dataset$fecundity=dataset$total_offspring/dataset$Num_moms

write.csv(dataset,file="Experiment3_data.csv")
pdf("Experiment1_vialsnotremoved.pdf")

#Fecundity figure
Fecund_figure=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=as.factor(Treatment),label=Num_moms),data=dataset)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_update(text = element_text(size=20))
Fecund_figure= Fecund_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-5","6-10","11-15","16-20"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

dev.off()

#Let's start with cleanedup data for Exp 1:
dataset2=read.csv("Experiment3_data.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)

#get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$Day..letter.of.vial.,mean)

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
summary(fit)
#Anova
anova(fit,test="Chisq")
#posthoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr

#Now do the same, but with recombination rate

#sum of crossovers in intervals 1 
dataset2$num_CO=dataset2$SCO1.sum
#sum of non-crossovers in intervals 1
dataset2$num_NCO_1=dataset2$NCO1.sum
#total recombination rate
dataset2$rec_rate_total=dataset2$SCO1.sum/dataset2$total_offspring

#get mean recombination
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_total,dataset2$Day..letter.of.vial.,mean)
fit_recrate=glm(rec_rate_total~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
summary(fit_recrate)
#Anova
anova(fit_recrate,test="Chisq")
#Posthoc
fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec

#logistic regression, similar to a t-test for count data
dataset2$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset2$F1.Vial))
dataset2$F1.Vial=as.numeric(dataset2$F1.Vial)
fit2=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment*Day..letter.of.vial.,data=dataset2,
         family=binomial(link="logit"))
summary(fit2)
coefs=coef(fit2)
coefs

#repeat, but only for day a
fit2a=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
            family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="A"))
summary(fit2a)
exp(coef(fit2a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

#repeat, but only for day b
fit2b=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="B"))
summary(fit2b)
exp(coef(fit2b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day b

#repeat, but only for day c
fit2c=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="C"))
summary(fit2c)
exp(coef(fit2c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c

#repeat, but only for day d
fit2d=glm(cbind(num_CO,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="D"))
summary(fit2d)
exp(coef(fit2d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day d


#Summarizing the data

length(unique(dataset2$F1.Vial[dataset2$Treatment=="20"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="25"]))

median(dataset2$Num_moms[dataset2$Treatment=="20"])
median(dataset2$Num_moms[dataset2$Treatment=="25"])


sum(dataset2$total_offspring[dataset2$Treatment=="20"])
sum(dataset2$total_offspring[dataset2$Treatment=="25"])

#Odd ratio plots
exp1=read.csv(file="exp1.csv", header=T)
exp1$time=as.character(exp1$time)
pdf("odd1.pdf")
# Basic line plot with points
exp1_odds=ggplot(data=exp1, aes(x=time, y=odds, group=1))+
  geom_line(size=2)+
  theme_bw()+
  theme_update(text = element_text(size=40))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point()+
  geom_errorbar(aes(ymin=odds-sd, ymax=odds+sd), width=.2, position=position_dodge(0.05))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(0.5,2.5)+
  ylab("Odds")+
  scale_x_discrete(name="Days", labels=c("1-5","6-10","11-15","16-20"))+
  ggtitle("Odd Ratios for experiment 1")

exp1_odds
dev.off()
