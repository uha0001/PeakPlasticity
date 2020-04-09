#load the datasets
yv=read.csv("Exp3_rawdata.csv", header= TRUE, na.strings="-")
backcross=read.csv("Exp3bc_setup.csv", header=T, na.strings="-")



#backcross$Treatment=as.numeric(as.character(backcross$Treatment)) #to make sure R reads it as a character rather than a number.
yv$wildtype=as.numeric(as.character(yv$wildtype))
yv$yellow.vermillion=as.numeric(as.character(yv$yellow.vermillion))
yv$yellow.only=as.numeric(as.character(yv$yellow.only))
yv$vermillion.only=as.numeric(as.character(yv$vermillion.only))
yv$wildtype.1=as.numeric(as.character(yv$wildtype.1))
yv$yellow.vermillion.1=as.numeric(as.character(yv$yellow.vermillion.1))
yv$yellow.only.1=as.numeric(as.character(yv$yellow.only.1))
yv$vermillion.only.1=as.numeric(as.character(yv$vermillion.only.1))
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

yv$SCO1=yv$yellow.only+yv$yellow.only.1
yv$SCO2=yv$vermillion.only+yv$vermillion.only.1
yv$NCO1=yv$wildtype+yv$wildtype.1
yv$NCO2=yv$yellow.vermillion+yv$vermillion.only.1

sco_count=sum(yv$SCO1, na.rm = TRUE)+sum(yv$SCO2, na.rm = TRUE)
nco_count=sum(yv$NCO1+yv$NCO2, na.rm = TRUE)
num_samples=sum(nco_count+sco_count, na.rm = TRUE)

#rate between yellow and vermillion

(sum(yv$SCO1, na.rm = TRUE)+sum(yv$SCO2, na.rm = TRUE))/num_samples

#merge with treatment data
yv_merged <- merge(yv, backcross, by.x = "Vial..", by.y = "Vial.number", all=T)
yv_merged = na.omit(yv_merged)
colnames(yv_merged)

dataset=summaryBy(SCO1+SCO2+NCO1+NCO2~F1.Vial+Day..letter.of.vial.+Treatment,data=yv_merged, FUN=sum,na.rm=T)

#add in a column for total offspring

dataset$total_offspring=dataset$SCO1.sum + dataset$SCO2.sum + dataset$NCO1.sum + dataset$NCO2.sum

num_moms=vector(mode="numeric",length=length(dataset$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset$F1.Vial)) { 
  f1_vial=dataset$F1.Vial[h]
#  f1_vial=as.numeric(gsub("[^0-9\\.]","",f1_vial))
#  print(f1_vial)
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$Vial..)))
  #store result in the vector
  num_moms[h]=mom_ct
}
dataset$Num_moms=num_moms

dataset$fecundity=dataset$total_offspring/dataset$Num_moms

write.csv(dataset,file="Experiment4_data.csv")

#Fecundity_figure
pdf("Experiment4_vialsnotremoved.pdf")

Fecund_figure=ggplot(aes(y=fecundity,x=Day..letter.of.vial., col=Treatment,label=Num_moms),data=dataset)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_update(text = element_text(size=24))
#Now add points and change labels for x-axis
Fecund_figure= Fecund_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

dev.off()
#remove vials with few progeny

dataset=dataset[dataset$total_offspring>=10,]


dataset$rec_rate_total=(dataset$SCO1.sum+dataset$SCO2.sum)/(dataset$NCO1.sum+dataset$NCO2.sum+dataset$SCO1.sum+dataset$SCO2.sum)

pdf("Experiment3_fecundity.pdf")

#Recomb_figure_total

Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day..letter.of.vial., col=Treatment,label=total_offspring),data=dataset)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
#Recomb_figure

#Now add points and change labels for x-axis

Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
#Recomb_figure

#Add a line through the median of the points

Recomb_figure=Recomb_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)

#Add a label for sample size

Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0,1.2)
Recomb_figure

#Recomb_figure_y_and_vermillion

Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day..letter.of.vial., col=Treatment,label=total_offspring),data=dataset)+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))
Recomb_figure=Recomb_figure+ggtitle("y-v Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure


dev.off()
#summarizing the data
length(unique(dataset2$F1.Vial[dataset2$Treatment=="20°"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="25°"]))

median(dataset2$Num_moms[dataset2$Treatment=="20°"])
median(dataset2$Num_moms[dataset2$Treatment=="25°"])


sum(dataset2$total_offspring[dataset2$Treatment=="20°"])
sum(dataset2$total_offspring[dataset2$Treatment=="25°"])

#Odds ratios and statistics

#Let's start with Exp 3:
dataset2=read.csv("exp3_cleanedup.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)
recrate=subset(dataset2, dataset2$Treatment=="20°")
mean(recrate$rec_rate_yv)
#get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$Day,mean)
(sum(dataset2$SCO1.sum)+sum(dataset2$SCO2.sum))/sum(dataset2$SCO1.sum+dataset2$SCO2.sum+dataset2$NCO1.sum+dataset2$NCO2.sum)
#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
summary(fit)
#anova
anova(fit, test="Chisq")
#posthoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr


#Now do the same, but with recombination rate

#remove vials with few progeny
dataset2=dataset2[dataset2$total_offspring>=10,]

#sum of crossovers in yv
dataset2$num_CO_1=dataset2$SCO1.sum+dataset2$SCO2.sum


#sum of non-crossovers in yv
dataset2$num_NCO_1=dataset2$total_offspring-(dataset2$SCO1.sum+dataset2$SCO2.sum)

#total recombination rate
#dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_yv=(dataset2$SCO1.sum+dataset2$SCO2.sum)/dataset2$total_offspring

#get mean recombination
tapply(dataset2$rec_rate_yv,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_yv,dataset2$Day,mean)

fit_recrate=glm(rec_rate_yv~Treatment*Day..letter.of.vial.,data=dataset2,family=quasipoisson)
#Anova
summary(fit_recrate)
anova(fit_recrate, test="Chisq")
#posthoc
fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="Day..letter.of.vial.", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec
onlyf=subset(dataset2, dataset2$Day..letter.of.vial.=="F" & dataset2$Treatment=="20°")
#Y-St REGION

#dataset2$num_CO_1=as.character(as.numeric(dataset2$num_CO_1))
#dataset2$num_NCO_1=as.character(as.numeric(dataset2$num_NCO_1))
dataset2$F1.Vial=as.numeric(gsub("[^0-9\\.]","",dataset2$F1.Vial))
print(dataset2$F1.Vial)
#logistic regression, similar to a t-test for count data
fit2=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
         family=binomial(link="logit"))
summary(fit2)
anova(fit2)
coefs=coef(fit2)
coefs
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#repeat, but only for day A
fita=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
         family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="A"))
summary(fita)
exp(coef(fita)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day a

exp(-0.07490)
#repeat, but only for day B
fit2b=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="B"))
summary(fit2b)
exp(coef(fit2b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day w


#repeat, but only for day C
fit2c=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="C"))
summary(fit2c)

exp(coef(fit2c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day x


#repeat, but only for day D
fit2d=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="D"))
summary(fit2d)
exp(coef(fit2d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y


#repeat, but only for day E
fit2e=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="E"))
summary(fit2e)
exp(coef(fit2e)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y

#repeat, but only for day F
fit2f=glm(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day..letter.of.vial.=="F"))
summary(fit2f)
exp(coef(fit2f)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y

#odds ratios are entered excel manually
exp3=read.csv(file="exp3.csv", header=T)
exp3$time=as.character(exp3$time)

pdf("Odds3.pdf")
exp3_odds=ggplot(data=exp3, aes(x=time, y=odds, group=1)) +
  geom_line(size=2)+
  theme_bw()+
  theme_update(text = element_text(size=40))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point()+
  geom_errorbar(aes(ymin=odds-sd, ymax=odds+sd), width=.2, position=position_dodge(0.05))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(0.5,2.6)+
  ylab("odds")+
  scale_x_discrete(name="Days", labels=c("1-2","3-4","5-6","7-9","10-12","13-15"))+
  ggtitle("Odd Ratios for experiment 3")
exp3_odds
dev.off()

