#we have changed the file format to csv and make R read the dataset
yv=read.csv("Exp4_rawdata.csv",header=T)
bc_worksheet=read.csv("Exp_backcrosses.csv",header=T,stringsAsFactors = F) 

#recode phenotypes as characters we told R they were 1 s and 0 s so they were defined
#as characters R can understand.
#that our treatment is not a number and to treat it as a character.

bc_worksheet$Treatment=as.character(bc_worksheet$Treatment)

female_counts=read.csv(file="Exp4_female.csv",header=T,stringsAsFactors = F)
yv$sd=as.character(yv$sd)
yv$y=as.character(yv$y)
yv$se=as.character(yv$se)

#check for errors
yv[yv$co_class=="error_ctsd",]

#add the mandatory packages for analysis

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

#CO groups were defined, There should be 14 groups
yv$co_class=ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y==yv$se,"non_CO", 
                          ifelse(yv$sd==yv$ct & yv$sd!=yv$y & yv$y==yv$se,"single_CO_2",
                                 ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y!=yv$se,"single_CO_3",
                                        ifelse(yv$ct==yv$sd & yv$sd!=yv$y & yv$y!=yv$se,"double_CO_3",
                                            "error"))))
#check for errors

yv[yv$co_class=="error",]
#add columns for counting
yv$NCO=ifelse(yv$co_class=="non_CO",1,0)
yv$SCO_2=ifelse(yv$co_class=="single_CO_2",1,0)
yv$SCO_3=ifelse(yv$co_class=="single_CO_3",1,0)
yv$DCO_3=ifelse(yv$co_class=="double_CO_3",1,0)
yv$NCO=ifelse(yv$co_class=="non_CO",1,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.

sco_count=sum(yv$numbMales*yv$SCO_2, na.rm = TRUE)+sum(yv$numbMales*yv$SCO_3, na.rm = TRUE)
dco_count=sum(yv$numbMales*yv$DCO_1, na.rm = TRUE)+sum(yv$numbMales*yv$DCO_2, na.rm = TRUE)+sum(yv$numbMales*yv$DCO_3, na.rm = TRUE)
nco_count=sum(yv$numbMales*yv$NCO, na.rm = TRUE)
num_samples=sum(yv$numbMales , na.rm = TRUE)+181

#total rate

(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow

(sum(yv$numbMales*yv$SCO_2,na.rm = TRUE)+sum(yv$numbMales*yv$DCO_2,na.rm = TRUE)+sum(yv$numbMales*yv$DCO_1,na.rm = TRUE))/num_samples

#rate between sepia and yellow

(sum(yv$numbMales*yv$SCO_3, na.rm = TRUE)+sum(yv$numbMales*yv$DCO_3,na.rm = TRUE)+sum(yv$numbMales*yv$DCO_2,na.rm = TRUE)+sum(yv$numbMales*yv$TCO_1, na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s
#single crossovers 1s
#double crossovers 2s

yv$num_co=ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y==yv$se,0,
                 ifelse(yv$sd==yv$ct & yv$sd!=yv$y & yv$y==yv$se,yv$numbMales*1,
                               ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y!=yv$se,yv$numbMales*1,
                                      ifelse(yv$ct==yv$sd & yv$sd!=yv$y & yv$y!=yv$se,yv$numbMales*2,
                                                        "NA"))))
#This is for summarizing our data
yv$male=c(yv$numbMales)


#merge with treatment data
yv_merged <- merge(yv, bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#merge female data
female_merged=merge(female_counts, bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)
dataset=summaryBy(male+num_co+SCO_2+SCO_3+ DCO_3+NCO~F1.Vial+Day+Treatment,data=yv)
#summarize long form data

dataset=summaryBy(male+num_co+SCO_2+SCO_3+ DCO_3+NCO~F1.Vial+Day+Treatment,data=yv_merged, FUN=sum,na.rm=T)

#add in female data
female_short=summaryBy(Numbfemales~F1.Vial+Day+Treatment,data=female_merged, FUN=sum,na.rm=T)
female_short=na.omit(female_short)

#one more merge
dataset2=merge(female_short,dataset, by=c("F1.Vial","Day","Treatment"))

#add in a column for total offspring

dataset2$total_offspring=dataset2$Numbfemales.sum + dataset2$male.sum

#add column for male-to-female ratio

dataset2$mf_ratio=dataset2$male.sum/dataset2$Numbfemales.sum

#make a vector to store the data

num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial))

#loop through dataset to get the count
for (h in 1:length(dataset2$F1.Vial)) { 
  f1_vial=dataset2$F1.Vial[h]
  mom_ct=length(unique(sort(subset(yv_merged,yv_merged$F1.Vial==f1_vial)$ViaNumber)))
  #store result in the vector
  num_moms[h]=mom_ct
}
#add vector as a column in dataset2

dataset2$Num_moms=num_moms

#use data to get fecundity calculation

dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms

#We will write our data into a file so we can further explore it in excel. 

write.csv(dataset2,file="Experiment8_data.csv")


#Fecundity_figure
pdf("Experiment4_fecundity.pdf")

Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_update(text = element_text(size=30))
#Now add points and change labels for x-axis
Fecund_figure= Fecund_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

dev.off()

#remove vials with few progeny

dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1, 2, and 3 

#dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO_1.sum+dataset2$DCO_2.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO_3.sum
dataset2$num_CO_3=dataset2$SCO_3.sum+dataset2$DCO_3.sum

#sum of non-crossovers in intervals 1, 2, and 3

#dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO_1.sum+dataset2$DCO_2.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO_3.sum)
dataset2$num_NCO_3=dataset2$male.sum-(dataset2$SCO_3.sum+dataset2$DCO_3.sum)

#total recombination rate

dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_2.sum+dataset2$DCO_3.sum)/dataset2$male.sum
#dataset2$rec_rate_sdct=(dataset2$SCO_1.sum+dataset2$DCO_1.sum+dataset2$DCO_2.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_3.sum+dataset2$DCO_3.sum)/dataset2$male.sum

pdf("Experiment8_recrate.pdf")

#Recomb_figure_total

Recomb_figure=ggplot(aes(y=rec_rate_total,x=Day, col=Treatment,label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
Recomb_figure

#Now add points and change labels for x-axis

Recomb_figure=Recomb_figure+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure

#Add a line through the median of the points

Recomb_figure=Recomb_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)

#Add a label for sample size

Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0.5,1.2)
Recomb_figure
#Fecundity_figure

Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_base()

Fecund_figure=Fecund_figure+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

#Male female ratio figure (add horizontal line at 1)

mfratio_figure=ggplot(aes(y=mf_ratio,x=Day,col=Treatment,label=Num_moms),data=dataset2)+geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+ylab("Male-female ratio")+theme_base()+ggtitle("Male-female ratio vs. Days post-mating")

mfratio_figure=mfratio_figure+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)

mfratio_figure=mfratio_figure+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
mfratio_figure  

#Recomb_figure_yellow_scalloped

Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_scalloped_cut

Recomb_figure=ggplot(aes(y=rec_rate_sdct,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("sd-ct Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_yellow_sepia

Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("green", "purple"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("1-3","4-6","7-9","10-12"))
Recomb_figure=Recomb_figure+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

dev.off()

#rough mean rec rate per treatment calculations

mean(dataset2$rec_rate_total[dataset2$Treatment== "26"],na.rm=T)
mean(dataset2$rec_rate_total[dataset2$Treatment== "21"],na.rm=T)

mean(dataset2$total_offspring[dataset2$Treatment=="21"],na.rm=T)
mean(dataset2$total_offspring[dataset2$Treatment=="26"],na.rm=T)

mean(dataset2$rec_rate_total[dataset2$Treatment== "26" & dataset2$Day== "C"],na.rm=T)
mean(dataset2$rec_rate_total[dataset2$Treatment== "21" & dataset2$Day== "C"],na.rm=T)

mean(dataset2$rec_rate_ysd[dataset2$Treatment== "26"],na.rm=T)
mean(dataset2$rec_rate_ysd[dataset2$Treatment== "21"],na.rm=T)
sum(dataset2$total_offspring[dataset2$Treatment=="21"],na.rm=T)

mean(dataset2$rec_rate_ysd[dataset2$Treatment== "26" & dataset2$Day== "C"],na.rm=T)
mean(dataset2$rec_rate_ysd[dataset2$Treatment== "21" & dataset2$Day== "C"],na.rm=T)

mean(dataset2$rec_rate_sdct[dataset$Treatment== "26"],na.rm=T)
mean(dataset2$rec_rate_sdct[dataset$Treatment== "21"],na.rm=T)

mean(dataset2$rec_rate_sdct[dataset2$Treatment== "26" & dataset2$Day== "C"],na.rm=T)
mean(dataset2$rec_rate_sdct[dataset2$Treatment== "21" & dataset2$Day== "C"],na.rm=T)

mean(dataset2$rec_rate_yse[dataset$Treatment== "26"],na.rm=T)
mean(dataset2$rec_rate_yse[dataset$Treatment== "21"],na.rm=T)

mean(dataset2$rec_rate_yse[dataset2$Treatment== "26" & dataset2$Day== "C"],na.rm=T)
mean(dataset2$rec_rate_yse[dataset2$Treatment== "21" & dataset2$Day== "C"],na.rm=T)


length(unique(dataset2$F1.Vial[dataset2$Treatment=="21"]))
length(unique(dataset2$F1.Vial[dataset2$Treatment=="26"]))

median(dataset2$Num_moms[dataset2$Treatment=="21"])
median(dataset2$Num_moms[dataset2$Treatment=="26"])


sum(dataset2$male.sum[dataset2$Treatment=="21"])
sum(dataset2$male.sum[dataset2$Treatment=="26"])

###Statistics for experiment 4

#Let's start with Exp 3:
dataset2=read.csv("Exp4_cleanedup.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)
dataset2=na.omit(dataset2)

#get mean fecundity
tapply(dataset2$fecundity,dataset2$Treatment,mean)
tapply(dataset2$fecundity,dataset2$Day,mean)

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit)
anova(fit,test="Chisq")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr


#Now do the same, but with recombination rate

#remove vials with few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1, 2, and 3 
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO_3.sum
dataset2$num_CO_3=dataset2$SCO_3.sum+dataset2$DCO_3.sum

#sum of non-crossovers in intervals 1, 2, and 3
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO_3.sum)
dataset2$num_NCO_3=dataset2$male.sum-(dataset2$SCO_3.sum+dataset2$DCO_3.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_2.sum+dataset2$DCO_3.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_3.sum+dataset2$DCO_3.sum)/dataset2$male.sum


#get mean recombination
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_total,dataset2$Day,mean)
tapply(dataset2$rec_rate_yse,dataset2$Treatment,mean)
fit_recrate=glm(rec_rate_total~Treatment*Day,data=dataset2,family=quasipoisson)
summary(fit_recrate)
anova(fit_recrate,test="Chisq")

fit_contrast_rec <- emmeans::emmeans(fit_recrate, "Treatment", by="Day", mode="kenward-roger")
fit_contr_rec <- contrast(fit_contrast_rec, method="trt.vs.ctrl")

pheno_contr_rec <- as.data.frame(summary(fit_contr_rec))
pheno_contr_rec

#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*Day,data=dataset2,
         family=binomial(link="logit"))
coefs=coef(fit3)
coefs
anova(fit3, test="Chisq")
fit_contrast_3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr_3 <- contrast(fit_contrast_3, method="trt.vs.ctrl")

pheno_contr_3 <- as.data.frame(summary(fit_contr_3))
pheno_contr_3
#repeat, but only for day A
fit3a=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="A"))
summary(fit3a)
exp(coef(fit3a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day v

#repeat, but only for day B
fit3b=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="B"))
summary(fit3b)
exp(coef(fit3b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day w


#repeat, but only for day C
fit3c=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(new_day=="C"))
summary(fit3c)
exp(coef(fit3c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day c


#repeat, but only for day D
fit3d=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="D"))
summary(fit3d)
exp(coef(fit3d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y

#means contrast
fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr3))
pheno_contr


#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment*Day,data=dataset2,
         family=binomial(link="logit"))
summary(fit4)
coefs=coef(fit4)
coefs
anova(fit4, test="Chisq")

#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!


#means contrast
fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr4))
pheno_contr

#repeat, but only for day A
fit4a=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="A"))
summary(fit4a)
exp(coef(fit4a)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day v

#repeat, but only for day B
fit4b=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="B"))
summary(fit4b)
exp(coef(fit4b)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day w


#repeat, but only for day c
fit4c=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="C"))
summary(fit4c)
exp(coef(fit4c)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day x


#repeat, but only for day D
fit4d=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="D"))
summary(fit4d)
exp(coef(fit4d)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y


#Odds ratios entered manually into excel.
exp4=read.csv(file="exp4.csv", header=T)
exp4$time=as.character(exp4$time)

pdf("Odds4.pdf")
exp4_odds=ggplot(data=exp4, aes(x=time, y=sdy, group=1)) +
  geom_line(color= "red",size=2)+
  geom_errorbar(aes(ymin=sdy-sd1, ymax=sdy+sd1), width=.2, position=position_dodge(0.05))+
  geom_line(data=exp4, aes(x=time, y=yse, group=1), color= "blue",size=2)+
  geom_errorbar(aes(ymin=yse-sd2, ymax=yse+sd2), width=.2, position=position_dodge(0.05))+
  theme_bw()+
  theme_update(text = element_text(size=40))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(0.6,1.6)+
  ylab("odds")+
  scale_x_discrete(name="Days", labels=c("1-3","4-6","7-9","10-12"))+
  ggtitle("Odd Ratios for experiment 4")
exp4_odds
dev.off()