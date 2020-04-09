#we have changed the file format to csv and make R read the dataset
yv=read.csv("Exp5_rawdata.csv",header=T)
bc_worksheet=read.csv("Exp5_backcrosses.csv",header=T,stringsAsFactors = F) 
#recode phenotypes as characters we told R they were 1 s and 0 s so they were defined
#as characters R can understand.
#that our treatment is not a number and to treat it as a character.
female_counts=read.csv(file="Exp5_females.csv",header=T,stringsAsFactors = F)
bc_worksheet$Treatment=as.character(bc_worksheet$Treatment)

yv$sd=as.character(yv$sd)
yv$y=as.character(yv$y)
yv$se=as.character(yv$se)

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
                                        ifelse(yv$ct==yv$sd & yv$sd!=yv$y & yv$y!=yv$se,"double_CO_2",
                                                        "error"))))
#check for errors

yv[yv$co_class=="error",]

#add columns for counting
yv$SCO_2=ifelse(yv$co_class=="single_CO_2",1,0)
yv$SCO_3=ifelse(yv$co_class=="single_CO_3",1,0)
yv$DCO_2=ifelse(yv$co_class=="double_CO_2",1,0)
yv$NCO=ifelse(yv$co_class=="non_CO",1,0)
sum(yv$SCO_2*yv$numbMales)
sum(yv$SCO_3*yv$numbMales)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.

sco_count=sum(yv$numbMales*yv$SCO_2, na.rm = TRUE)+sum(yv$numbMales*yv$SCO_3, na.rm = TRUE)
dco_count=sum(yv$numbMales*yv$DCO_2, na.rm = TRUE)
num_samples=sum(yv$numbMales , na.rm = TRUE)
nco_count=sum(yv$numbMales*yv$NCO, na.rm = TRUE)
#total rate

(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow

(sum(yv$numbMales*yv$SCO_2,na.rm = TRUE)+sum(yv$numbMales*yv$DCO_2,na.rm = TRUE))/num_samples

#rate between sepia and yellow

(sum(yv$numbMales*yv$SCO_3, na.rm = TRUE)+sum(yv$numbMales*yv$DCO_2,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s
#single crossovers 1s
#double crossovers 2s

yv$num_co=ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y==yv$se,0,
                ifelse(yv$sd==yv$ct & yv$sd!=yv$y & yv$y==yv$se,yv$numbMales*1,
                               ifelse(yv$ct==yv$sd & yv$sd==yv$y & yv$y!=yv$se,yv$numbMales*1,
                                      ifelse(yv$ct==yv$sd & yv$sd!=yv$y & yv$y!=yv$se,yv$numbMales*2,
                                                           NA))))
#This is for summarizing our data
yv$male=c(yv$numbMales)

#merge with treatment data
yv_merged <- merge(yv, bc_worksheet, by.x = "ViaNumber", by.y = "Vial.Number", all=T)

#merge female data
female_merged=merge(female_counts, bc_worksheet, by.x = "Vial", by.y = "Vial.Number", all=T)

#summarize long form data

dataset=summaryBy(male+num_co+SCO_2+SCO_3+DCO_2~F1.Vial+Day+Treatment,data=yv_merged, FUN=sum,na.rm=T)


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

write.csv(dataset2,file="Exp5_cleanedup.csv")

#Fecundity_figure
pdf("Exp5_fecundity.pdf")

Fecund_figure=ggplot(aes(y=fecundity,x=Day, col=Treatment,label=Num_moms),data=dataset2)+ylab("# Progeny per mom")+ggtitle("Total Fecundity vs. Days post-mating")+theme_update(text = element_text(size=30))
#Now add points and change labels for x-axis
Fecund_figure= Fecund_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))


Fecund_figure=Fecund_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Fecund_figure=Fecund_figure+stat_summary(fun.y = median, geom="line",aes(group=Treatment),size=2)+geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Fecund_figure

dev.off()

#remove vials with few progeny

dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1, 2, and 3 

#dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO_1.sum+dataset2$DCO_3.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO_2.sum
dataset2$num_CO_3=dataset2$SCO_3.sum+dataset2$DCO_2.sum

#sum of non-crossovers in intervals 1, 2, and 3

#dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO_1.sum+dataset2$DCO_3.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO_2.sum)
dataset2$num_NCO_3=dataset2$male.sum-(dataset2$SCO_3.sum+dataset2$DCO_2.sum)

#total recombination rate

dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_2.sum+dataset2$DCO_2.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_3.sum+dataset2$DCO_2.sum)/dataset2$male.sum


pdf("Exp9_sdyboxplotonly.pdf")


Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_update(text = element_text(size=30))
#Now add points and change labels for x-axis
Recomb_figure= Recomb_figure+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("#f1a340", "#998ec3"))
Recomb_figure <- Recomb_figure + geom_boxplot() # Adds color
Recomb_figure <- Recomb_figure + scale_x_discrete(name="Days", labels=c("6","7","8","9","10")) # Adds kaveks
Recomb_figure # displays the boxplots

dev.off()

#Let's start with Exp 5 stats:
dataset2=read.csv("Exp5_cleanedup.csv",header=T,stringsAsFactors = F)
dataset2$Treatment=as.character(dataset2$Treatment)

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
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO_2.sum
dataset2$num_CO_3=dataset2$SCO_3.sum+dataset2$DCO_2.sum

#sum of non-crossovers in intervals 1, 2, and 3
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO_2.sum)
dataset2$num_NCO_3=dataset2$male.sum-(dataset2$SCO_3.sum+dataset2$DCO_2.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_2.sum+dataset2$DCO_2.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_3.sum+dataset2$DCO_2.sum)/dataset2$male.sum


#get mean recombination
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean)
tapply(dataset2$rec_rate_total,dataset2$Day,mean)
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
anova(fit3,test="Chisq")
fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
#repeat, but only for day V
fit3v=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="V"))
summary(fit3v)
exp(coef(fit3v)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day v

#repeat, but only for day W
fit3w=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="W"))
summary(fit3w)
exp(coef(fit3w)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day w


#repeat, but only for day X
fit3x=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="X"))
summary(fit3x)
exp(coef(fit3x)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day x


#repeat, but only for day Y
fit3y=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="Y"))
summary(fit3y)
exp(coef(fit3y)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y

#repeat, but only for day Z
fit3z=glm(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="Z"))
summary(fit3z)
exp(coef(fit3z)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day z


#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!


#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment*Day,data=dataset2,
           family=binomial(link="logit"))
summary(fit4)
coefs=coef(fit4)
coefs
anova(fit4,test="Chisq")
#negative 26 coefficent means odds of being recombinant is higher in other treatment (b/c it's negative)!

#means contrast
fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr4))
pheno_contr

#repeat, but only for day V
fit4v=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="V"))
summary(fit4v)
exp(coef(fit4v)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day v

#repeat, but only for day W
fit4w=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="W"))
summary(fit4w)
exp(coef(fit4w)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day w


#repeat, but only for day X
fit4x=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="X"))
summary(fit4x)
exp(coef(fit4x)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day x


#repeat, but only for day Y
fit4y=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="Y"))
summary(fit4y)
exp(coef(fit4y)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day y

#repeat, but only for day Z
fit4z=glm(cbind(num_CO_3,num_NCO_3)~(1|F1.Vial)+Treatment,data=dataset2,
          family=binomial(link="logit"),subset=c(Day=="Z"))
summary(fit4z)
exp(coef(fit4z)[3]) #At high temperature, odds of crossover are this times as likely as odds of crossover at low temperature for day z

#Odds ratios were entered manually

exp5=read.csv(file="exp5.csv", header=T)
exp5$time=as.character(exp5$time)
pdf("odds5.pdf")
exp5_odds=ggplot(data=exp5, aes(x=time, y=SDY, group=1)) +
  geom_line (color= "red", size=2)+
  geom_errorbar(aes(ymin=SDY-sd1, ymax=SDY+sd1), width=.2, position=position_dodge(0.05))+
  geom_line(data=exp5, aes(x=time, y=YSE, group=1), color= "blue",size=2)+
  geom_errorbar(aes(ymin=YSE-sd2, ymax=YSE+sd2), width=.2, position=position_dodge(0.05))+
  theme_bw()+
  theme_update(text = element_text(size=40))+
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=1, linetype="dashed", color = "grey",size=1)+
  ylim(0.6,1.6)+
  ylab("odds")+
  scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+
  ggtitle("Odd Ratios for experiment 5")
exp5_odds
dev.off()