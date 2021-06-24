#necessary packages are listed below, Everytime before running the scripts, the prerequisites
#should be uploaded.

library(ggplot2)
library(ggthemes)
library(emmeans)
library(lme4)
library(lmerTest)
library(doBy)
library(reshape2)
library(car)
#####
#Load the necessary data files.These should be in csv format.

ICE=read.csv("PhenotypingICExPlasticity.csv",header=T)
ICE_bc=read.csv("ICExPlasticity_backcrosses.csv",header=T,stringsAsFactors = F) 

#make separate datafiles for males and females. Remember only males should be mutant phenotype is present.

ICE$number.of.flies=ifelse(ICE$Sex.of.the.flies=="Male", ICE$Number.of.males, ICE$Number.of.females)
ICE_female=subset(ICE, ICE$Sex.of.the.flies=="Female", na.rm=TRUE)
ICE_male=subset(ICE, ICE$Sex.of.the.flies=="Male", na.rm=TRUE)
#Remove the repeated columns in the dataset

ICE_female <- ICE_female[,-15]
ICE_male <- ICE_male[,-15]
ICE_male <- ICE_male[,-16]
ICE_female <- ICE_female[,-7]

#Some of the vial names in F1 flies have V infront of them, the code velow removes those.
ICE_female$Vial.number=gsub("[^0-9\\.]","",ICE_female$Vial.number)
ICE_male$Vial.number=gsub("[^0-9\\.]","",ICE_male$Vial.number)

#Define Crossovers
ICE_male$co_class=ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia,"non_CO", 
                    ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia,"single_CO_1",
                           ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia,"single_CO_2",
                                  ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia,"double_CO",
                                         "error"))))

#check for errors, which are the removed cut phenotypes.
ICE[ICE$co_class=="error",]

#I am making sure the numbers are actually numbers. R doesn't like it
ICE_male$Number.of.males=as.numeric(ICE_male$Number.of.males)

#Haplotype analysis part 1 defining haplotype groups
#Define happlotype groups
ICE_male$gamete_class=ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Yellow=="Present","gt1a",
                             ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Yellow=="Absent","gt1b",
                                    ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Yellow=="Present","gt2a",
                                           ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Yellow=="Absent","gt2b",
                                                  ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia& ICE_male$Yellow=="Present","gt3a",
                                                         ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia& ICE_male$Yellow=="Absent","gt3b",
                                                                ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia& ICE_male$Yellow=="Present","gt4a",
                                                                       ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia& ICE_male$Yellow=="Absent","gt4b",
                                                                              "error"))))))))

#Sanity check, there should not be any error!
ICE_male[ICE_male$gamete_class=="error",]


#Observer bias
total=(tapply(ICE_male$Number.of.males,ICE_male$Name,sum, na.rm=T))
ICE_male$scallopedcnt=ifelse(ICE_male$Scalloped=="Present",ICE_male$Number.of.males,0)
ICE_male$yellowcnt=ifelse(ICE_male$Yellow=="Present",ICE_male$Number.of.males,0)
ICE_male$sepiacnt=ifelse(ICE_male$Sepia=="Present",ICE_male$Number.of.males,0)


ICE_male=subset(ICE_male, ICE_male$Letter.of.the.day==c("A","B","C","D","E","F"))
total=(tapply(ICE_male$Number.of.males,ICE_male$Vial.number,sd, na.rm=T))
ICE_male_badobs=(tapply(ICE_male$Number.of.males, ICE_male$Vial.number, sd,na.rm=T))
scalloped_sd=(tapply(ICE_male$scallopedcnt,ICE_male$Name,sd, na.rm=T)/total)*100
yellow_sd=(tapply(ICE_male$yellowcnt,ICE_male$Name,sd, na.rm=T)/total)*100
sepia_sd=(tapply(ICE_male$sepiacnt,ICE_male$Name,sd, na.rm=T)/total)*100
scalloped_sd=data.frame(scalloped_sd)
yellow_sd=data.frame(yellow_sd)
sepia_sd=data.frame(sepia_sd)
k=cbind(scalloped_sd,yellow_sd,sepia_sd)

df2 <- data.frame(phenotype=rep(c("scalloped","yellow","sepia"), each=10),
                  observer=rep(c("a","b","c","d","e","f","g","h","i","j"),3),
                  sd=c(0.24726979,1.23620614,0.52285211,0.21573412,0.38685277,0.17895677,0.89723581,0.51622318,0.12920541,0.05452998,0.23364624,1.09549989,0.52405423,0.21205176,0.40070764,0.18103302,0.92717161,0.54058514,0.12999266,0.05457178,0.22905714,1.21953985,0.47510231,0.21935828,0.38608879,0.18618752,0.92819169,0.63308597,0.11432229,0.05595982))

s2=ggplot(data=df2, aes(x=observer, y=sd, fill=phenotype)) +
  geom_bar(stat="identity")

pdf("observerstdsev.pdf")
s=ggplot(data=df2, aes(x=observer, y=sd, fill=phenotype)) +
  geom_bar(stat="identity",position=position_dodge())
s
dev.off()
scalloped=(tapply(ICE_male$scallopedcnt,ICE_male$Name,sum, na.rm=T)/total)*100
yellow=(tapply(ICE_male$yellowcnt,ICE_male$Name,sum, na.rm=T)/total)*100
sepia=(tapply(ICE_male$sepiacnt,ICE_male$Name,sum, na.rm=T)/total)*100

df1 <- data.frame(phenotype=rep(c("scalloped","yellow","sepia"), each=10),
                  observer=rep(c("a","b","c","d","e","f","g","h","i","j"),3),
                  mean=c(44.49339,44.55446,44.49153,50.24233,46.22951,45.60440,45.28302,42.67782,41.55154,45.64926,44.93392,47.52475,43.22034,50.72698,52.13115,49.03846,53.45912,44.35146,47.18385,47.61267,46.91630,52.47525,32.20339,50.08078,48.19672,49.45055,52.20126,47.69874,34.21892,47.74654),
                  sd=c(0.24726979,1.23620614,0.52285211,0.21573412,0.38685277,0.17895677,0.89723581,0.51622318,0.12920541,0.05452998,0.23364624,1.09549989,0.52405423,0.21205176,0.40070764,0.18103302,0.92717161,0.54058514,0.12999266,0.05457178,0.22905714,1.21953985,0.47510231,0.21935828,0.38608879,0.18618752,0.92819169,0.63308597,0.11432229,0.05595982))

pdf("observermean.pdf")

m<- ggplot(df1, aes(x=observer, y=mean, fill=phenotype)) + 
geom_bar(stat="identity", position=position_dodge()) +
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
position=position_dodge(.9))     
m
dev.off()

shapiro.test(as.numeric(df1$mean))
t.test()
#write.csv(a,"observerbias.csv")

#total=(tapply(ICE_male$Number.of.males,ICE_male$Name,sum, na.rm=T))
#sepia=tapply(ICE_male$sepiacnt,list(ICE_male$Name,ICE_male$Letter.of.the.day),sum, na.rm=T)
#yellow=tapply(ICE_male$yellowcnt,list(ICE_male$Name,ICE_male$Letter.of.the.day),mean, na.rm=T)
#scalloped=tapply(ICE_male$scallopedcnt,list(ICE_male$Name,ICE_male$Letter.of.the.day),mean, na.rm=T)
#b=rbind(sepia,yellow,scalloped)
#ICE_male=subset(ICE_male,ICE_male$Name==c("Ulku Huma Altindag", "Abigail Burnett", "Taylor Papstein-Novak","Hannah Taylor","Maxwell McLaughlin"))

#add columns for counting the haplotype groups
ICE_male$gt1a_111=ifelse(ICE_male$gamete_class=="gt1a",ICE_male$Number.of.males,0)
ICE_male$gt1b_000=ifelse(ICE_male$gamete_class=="gt1b",ICE_male$Number.of.males,0)
ICE_male$gt2a_011=ifelse(ICE_male$gamete_class=="gt2a",ICE_male$Number.of.males,0)
ICE_male$gt2b_100=ifelse(ICE_male$gamete_class=="gt2b",ICE_male$Number.of.males,0)
ICE_male$gt3a_110=ifelse(ICE_male$gamete_class=="gt3a",ICE_male$Number.of.males,0)
ICE_male$gt3b_001=ifelse(ICE_male$gamete_class=="gt3b",ICE_male$Number.of.males,0)
ICE_male$gt4a_010=ifelse(ICE_male$gamete_class=="gt4a",ICE_male$Number.of.males,0)
ICE_male$gt4b_101=ifelse(ICE_male$gamete_class=="gt4b",ICE_male$Number.of.males,0)

#define them as vectors so we can do a binomial test.
gt1a_111=sum(ICE_male$Number.of.males*ICE_male$gt1a_111)
gt1b_000=sum(ICE_male$Number.of.males*ICE_male$gt1b_000)
gt2a_011=sum(ICE_male$Number.of.males*ICE_male$gt2a_011)
gt2b_100=sum(ICE_male$Number.of.males*ICE_male$gt2b_100)
gt3a_110=sum(ICE_male$Number.of.males*ICE_male$gt3a_110)
gt3b_001=sum(ICE_male$Number.of.males*ICE_male$gt3b_001)
gt4a_010=sum(ICE_male$Number.of.males*ICE_male$gt4a_010)
gt4b_101=sum(ICE_male$Number.of.males*ICE_male$gt4b_101)

ICE_haplotype=cbind(gt1a_111,gt1b_000,gt2a_011,gt2b_100,gt3a_110,gt3b_001,gt4a_010,gt4b_101)

#haplotype analysis - statistical test to see if they are different.
binom.test(687,687+1097,0.5)
binom.test(198,198+210,0.5)
binom.test(479,479+510,0.5)
binom.test(124,124+80,0.5)
binom.test(gt1a_111,,p=0.5)
binom.test(gt2a_011,gt2b_100,p=0.5)
binom.test(gt3a_110,gt3b_001,p=0.5)



#BACK to THE PHENOTYPING DATA and CALCULATING CROSSOVER FREQUENCIES
#add columns for counting

ICE_male$NCO=ifelse(ICE_male$co_class=="non_CO",ICE_male$Number.of.males,0)
ICE_male$SCO_1=ifelse(ICE_male$co_class=="single_CO_1",ICE_male$Number.of.males,0)
ICE_male$SCO_2=ifelse(ICE_male$co_class=="single_CO_2",ICE_male$Number.of.males,0)
ICE_male$DCO=ifelse(ICE_male$co_class=="double_CO",ICE_male$Number.of.males,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(as.numeric(ICE_male$NCO),na.rm=T)
sco_count=sum(as.numeric(ICE_male$SCO_1), na.rm = TRUE)+sum(as.numeric(ICE_male$SCO_2), na.rm = TRUE)
dco_count=sum(as.numeric(ICE_male$DCO), na.rm = TRUE)
num_samples=sum(as.numeric(ICE_male$Number.of.males), na.rm = TRUE) 

#total rate
(sco_count+(2*dco_count))/num_samples

#rate between scalloped and yellow
(sum(as.numeric(ICE_male$SCO_1),na.rm = TRUE)+(sum(as.numeric(ICE_male$DCO),na.rm = TRUE)))/num_samples

#rate between sepia and yellow
(sum(as.numeric(ICE_male$SCO_2), na.rm = TRUE)+(sum(as.numeric(ICE_male$DCO),na.rm = TRUE)))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
ICE_male$num_co=ifelse(ICE_male$Yellow==ICE_male$Scalloped & ICE_male$Yellow==ICE_male$Sepia,0, 
                  ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia,1*ICE_male$Number.of.males, 
                         ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia,1*ICE_male$Number.of.males,  
                                ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia,2*ICE_male$Number.of.males, 
                                       NA))))
#need to manually add in treatment. This happened because I forgot to start a new backcross sheet for ICE treatment and continued from 75 NoICE backcross. 
#ICE: VWXYZ; NO-ICE: BCDEF
ICE_male$Treatment1=ifelse(grepl("[BCDEF]",ICE_male$Letter.of.the.day),"no-ICE","ICE")

#add a new vial number that renumbers the ICE to 76-150
ICE_male$newVial=ifelse(ICE_male$Treatment1=="ICE",as.numeric(ICE_male$Vial.number)+75,as.numeric(ICE_male$Vial.number))


#repeat for females
#ICE_female: VWXYZ; NO-ICE_female: BCDEF
ICE_female$Treatment1=ifelse(grepl("[BCDEF]",ICE_female$Letter.of.the.day),"no-ICE","ICE")

#add a new vial number that renumbers the ICE_female to 76-150
ICE_female$newVial=ifelse(ICE_female$Treatment1=="ICE",as.numeric(ICE_female$Vial.number)+75,as.numeric(ICE_female$Vial.number))

#Changes to make it similar with the peak plasticity code, also remove na.s and columns that are unnnecessary.
ICE_bc= ICE_bc[,-15:-21]
ICE_bc=na.omit(ICE_bc)
colnames(ICE_bc)[1]="Vial.number"

###Subset for the ones with more standard deviation

ICE_male1=subset(ICE_male,ICE_male$Name==c("Andrew Mayer", "Hannah Forget", "Natalia Rivera Rincon","Nicholas Morbidelli","Joshua Thomas"))
x=unique(ICE_male1$Vial.number)
x=as.numeric(x)
sort(x)
ICE_male<- ICE_male[!(ICE_male$Vial.number %in% x),]
vials_sd=(tapply(ICE_male$Number.of.males,ICE_male$Vial.number,sd, na.rm=T))

if(vials_sd > 0.05){
  remove("Non-negative number")
} else {
  print("Negative number")
}

#finally we can continue analyzing the recombination data.
#merge with treatment data
ICE_merged <- merge(ICE_male, ICE_bc, by.x = "newVial", by.y = "Vial.number", all=T)

keep=c("Ulku Huma Altindag", "Abigail Burnett", "Taylor Papstein-Novak","Hannah Taylor","Maxwell McLaughlin")

ICE_merged2=ICE_merged[ICE_merged$Name.x==keep,]


#ICE_male=subset(ICE_male,ICE_male$Name==c("Ulku Huma Altindag", "Abigail Burnett", "Taylor Papstein-Novak","Hannah Taylor","Maxwell McLaughlin"))

#merge female data
ICE_female_merged=merge(ICE_female, ICE_bc, by.x = "newVial", by.y = "Vial.number", all=T)
ICE_merged$F1.Vial.Number=as.numeric(gsub("[^0-9\\.]","",ICE_merged$F1.Vial.Number))

#summarize long form data
dataset=summaryBy(Number.of.males.x+num_co+gt1a_111+gt1b_000+gt2a_011+gt2b_100+gt3a_110+gt3b_001+gt4a_010+gt4b_101+SCO_1+SCO_2+DCO~+Letter.of.the.day+Treatment+F1.Vial.Number,data=ICE_merged, FUN=sum,na.rm=T)
dataset=na.omit(dataset)
ICE_female_merged$F1.Vial.Number=as.numeric(gsub("[^0-9\\.]","",ICE_female_merged$F1.Vial.Number))

#add in female data
ICE_female_merged$Number.of.females.x=as.numeric(ICE_female_merged$Number.of.females.x)
ICE_female_short=summaryBy(Number.of.females.x~+F1.Vial.Number+Letter.of.the.day+Treatment,data=ICE_female_merged, FUN=sum,na.rm=T)
ICE_female_short=na.omit(ICE_female_short)

#one more merge
dataset2=merge(ICE_female_short,dataset, by=c("F1.Vial.Number","Letter.of.the.day","Treatment"))
dataset2=na.omit(dataset2)

#add in a column for total offspring
dataset2$total_offspring=dataset2$Number.of.females.x.sum + dataset2$Number.of.males.x.sum

#make a vector to store the data
num_moms=vector(mode="numeric",length=length(dataset2$F1.Vial.Number))

#loop through dataset to get the count
for (h in 1:length(dataset2$F1.Vial.Number)) { 
  f1_vial=dataset2$F1.Vial.Number[h]
  mom_ct=length(unique(sort(subset(ICE_merged,ICE_merged$F1.Vial.Number==f1_vial)$Vial.number)))
  #store result in the vector
  num_moms[h]=mom_ct
}

#add vector as a column in dataset2
dataset2$Num_moms=num_moms

#use data to get fecundity calculation
dataset2$fecundity=dataset2$total_offspring/dataset2$Num_moms

#We will write our data into a file so we can read it in for later analysis. 
write.csv(dataset2,file="ICE_cleanedup.csv")

#additional corrections of the dataset
# this is because R doesn't recognize that day b and day v are the same time point.
dataset2$Day=dataset2$Letter.of.the.day
levels(as.factor(dataset2$Letter.of.the.day))
dataset2$Letter.of.the.day=as.factor(dataset2$Letter.of.the.day)
levels(dataset2$Letter.of.the.day)=c("B", "C","D", "E", "F","B", "C","D","E","F")

###Haplotype analysis
#gamete analysis by treatment

gt111=tapply(dataset2$gt1a_111.sum,dataset2$Treatment,sum,na.rm=T)

gt000=tapply(dataset2$gt1b_000.sum,dataset2$Treatment,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Treatment,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Treatment,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Treatment,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Treatment,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Treatment,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Treatment,sum,na.rm=T)
i=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(i,"gamete_class_by_treatment.csv")
#gamete analysis by day
gt111=tapply(dataset2$gt1a_111.sum,dataset2$Day,sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,dataset2$Day,sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,dataset2$Day,sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,dataset2$Day,sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,dataset2$Day,sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,dataset2$Day,sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,dataset2$Day,sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,dataset2$Day,sum,na.rm=T)
s=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(s,"gamete_class_byday.csv")
gt111=tapply(dataset2$gt1a_111.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt000=tapply(dataset2$gt1b_000.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt011=tapply(dataset2$gt2a_011.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt100=tapply(dataset2$gt2b_100.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt110=tapply(dataset2$gt3a_110.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt001=tapply(dataset2$gt3b_001.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt010=tapply(dataset2$gt4a_010.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
gt101=tapply(dataset2$gt4b_101.sum,list(dataset2$Treatment,dataset2$Day),sum,na.rm=T)
t=rbind(gt111, gt000, gt011, gt100, gt110, gt001, gt010, gt101)
write.csv(t,"gamete_class_byday_and_treatment.csv")


#Fecundity stats

#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Letter.of.the.day,data=dataset2,family=quasipoisson)
#summary(fit)
anova_fec=anova(fit, test="Chisq")
write.csv(anova_fec,"ICE_fecundity_model_table.csv")

#posthoc
fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Letter.of.the.day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"ICE_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))

pdf("ICE_Fecundity.pdf")

Fecund_figure_ICE=ggplot(dataset2, aes(x=Letter.of.the.day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment ICE")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot(lwd=2)+scale_x_discrete(name="Days", labels=c("6","7","8","9","10","6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("green","orange","red","blue"))
Fecund_figure_ICE

dev.off()

tapply(dataset2$fecundity,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$fecundity,dataset2$Letter.of.the.day,mean,na.rm=T)

#remove vials with few progeny
dataset2=dataset2[dataset2$Number.of.males.x.sum>=10,]
dataset2=dataset2[dataset2$Number.of.females.x.sum>=10,]

#sum of crossovers in intervals 1 & 2 
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$Number.of.males.x.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$Number.of.males.x.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$Number.of.males.x.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$Number.of.males.x.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$Number.of.males.x.sum
mean(dataset2$rec_rate_ysd)
mean(dataset2$rec_rate_yse)

#Kosambi correction
dataset2$kosambi_rec_rate_ysd=((2.71^(4*dataset2$rec_rate_ysd)-1)/2*(2.71^(-4*dataset2$rec_rate_ysd)+1))/10
dataset2$kosambi_rec_rate_yse=((2.71^(4*dataset2$rec_rate_yse)-1)/2*(2.71^(-4*dataset2$rec_rate_yse)+1))/10
levels(as.factor(dataset2$Treatment))

sum(dataset2$total_offspring[dataset2$Treatment=="Control (21)"])
sum(dataset2$total_offspring[dataset2$Treatment=="temp (26)"])
sum(dataset2$total_offspring[dataset2$Treatment=="ICE (ICE only at 21)"])
sum(dataset2$total_offspring[dataset2$Treatment=="ICExTemp (ICE in 26)"])

#Recomb_figure_total
Recomb_figure=ggplot(aes(y=rec_rate_total,x=Letter.of.the.day, col=as.factor(Treatment),label=Number.of.males.x.sum),data=dataset2)+ylab("% recombination")+ggtitle("Total Recombination rate vs. Days post-mating")+theme_base()
Recomb_figure=Recomb_figure+scale_colour_manual(values=c("blue", "red","yellow","green"))+geom_point(alpha=0.6,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)
Recomb_figure=Recomb_figure+ggtitle("Total Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3) +ylim(0.2,1)
Recomb_figure


pdf("ICEysd and yse combined.pdf")
#Recomb_figure_yellow_scalloped
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_ysd,x=Letter.of.the.day, col=as.factor(Treatment),label=Number.of.males.x.sum),data=dataset2)+scale_colour_manual(values=c("blue", "green", "orange","red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day")
Recomb_figure=Recomb_figure+ylim(0,0.5)+ggtitle("y-sd Recombination rate") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10","6","7","8","9","10"))
Recomb_figure

#Recomb_figure_yellow_sepia
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_yse,x=Letter.of.the.day, col=as.factor(Treatment),label=Number.of.males.x.sum),data=dataset2)+scale_colour_manual(values=c("blue", "green", "orange","red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day")
Recomb_figure=Recomb_figure+ylim(0,1)+ggtitle("y-se Recombination rate") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10","6","7","8","9","10"))
Recomb_figure
dev.off()
#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_ysd,list(dataset2$Treatment, dataset2$Letter.of.the.day=="F"),mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Letter.of.the.day,mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_yse,dataset2$Letter.of.the.day,mean,na.rm=T)

#under construction I need to make a pllot without the day effect just treatment
####
pdf("ICEmain.pdf")
ICE_figure=ggplot(aes(y=rec_rate_total,x=as.factor(Treatment), col=as.factor(Treatment),label=Number.of.males.x.sum),data=dataset2)+scale_colour_manual(values=c("blue", "orange", "red","green"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Treatment")
ICE_figure=ICE_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)+scale_x_discrete(name="Treatment",labels=c("control","ICEonly","ICExHeat","Heat"))
ICE_figure
dev.off()
#These can ve used to summarize the results in the manuscript.
tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Letter.of.the.day,mean,na.rm=T)

tapply(dataset2$kosambi_rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$kosambi_rec_rate_yse,dataset2$Letter.of.the.day,mean,na.rm=T)

###Odds ratios
#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial.Number)+Treatment*Letter.of.the.day,data=dataset2,
           family=binomial(link="logit"))
#coefs=coef(fit3)
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"ICE_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Letter.of.the.day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"ICE_recrate_sd-y_posthoc_table.csv")

#we can extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
y_sd_or=exp(pheno_contr3$estimate)
y_sd_error=(pheno_contr3$SE)

#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial.Number)+Treatment*Letter.of.the.day,data=dataset2,
           family=binomial(link="logit"))
#summary(fit4)
#coefs=coef(fit4)
#coefs
anova_yse=Anova(fit4,test="Chisq")
write.csv(anova_yse,"ICE_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Letter.of.the.day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"ICE_recrate_y-se_posthoc_table.csv")

#Odds ratio doesn't work for the four catagories, we need to use the actual recomvination data,
#can ignore till line 362
y_se_or=exp(pheno_contr4$estimate)
y_se_error=(pheno_contr4$SE)


#setup data frame
odds_ratios=cbind(y_sd_or,y_se_or)
colnames(odds_ratios)=c("sd-y","y-se")
odds_ratios=na.omit(odds_ratios)
rownames(odds_ratios)=c("6","7","8","9","10")

#melt dataframe into long form
odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(y_sd_error,y_se_error)
x=cbind(odds_ratios,SE)

#convert p-values to stars for plot
ysd_sig=ifelse(pheno_contr3$p.value<0.001,"***",ifelse(pheno_contr3$p.value<0.01,"**",ifelse(pheno_contr3$p.value<0.05,"*","")))
yse_sig=ifelse(pheno_contr4$p.value<0.001,"***",ifelse(pheno_contr4$p.value<0.01,"**",ifelse(pheno_contr4$p.value<0.05,"*","")))

#add significance to table
sig=c(ysd_sig,yse_sig)
y=cbind(x,sig)
odds_ratios=y

pdf("ICE_odds_ratio.pdf")
odds_figure_ICE=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.75,1.76)+
  geom_line(size=2)+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ggtitle("Experiment ICE")+theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=6, y=1.65, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=7, y=1.65, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=8, y=1.65, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=ysd_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=yse_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=6, y=1.6, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=7, y=1.6, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=8, y=1.6, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_ICE
dev.off()
####

#boxplot for total recombination rate
pdf("ICE_sdyboxplotonly.pdf")

Recomb_figure_ysd=ggplot(aes(y=kosambi_rec_rate_ysd,x=Letter.of.the.day, col=Treatment,label=Number.of.males.x.sum),data=dataset2)+scale_colour_manual(values=c("green", "orange", "red","blue"))+geom_boxplot(size=1.5)+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure_ysd=Recomb_figure_ysd+ggtitle("sd-y Interval")+ylim(0,0.6)
Recomb_figure_ysd=Recomb_figure_ysd+annotate(geom="text", x=1, y=0.58, label=ysd_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=ysd_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=ysd_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=ysd_sig[4],size=10)
Recomb_figure_ysd

dev.off()


pdf("ICE_yseboxplotonly.pdf")


Recomb_figure_yse=ggplot(aes(y=kosambi_rec_rate_yse,x=Letter.of.the.day, col=Treatment,label=),data=dataset2)+scale_colour_manual(values=c("green", "orange", "red","blue"))+geom_boxplot(size=1.5)+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure_yse=Recomb_figure_yse+ggtitle("y-se Interval")+ylim(0.1,0.6)
Recomb_figure_yse=Recomb_figure_yse+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
Recomb_figure_yse

dev.off()

#Crossover interference
dataset2$Exp_DCO=(((dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$Number.of.males.x.sum)*((dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$Number.of.males.x.sum))
dataset2$Obs_DCO=dataset2$DCO.sum/dataset2$Number.of.males.x.sum

dataset2$COC=dataset2$Obs_DCO/dataset2$Exp_DCO
dataset2$Interference=1-dataset2$COC

fit5=glm(Interference~Treatment*Letter.of.the.day,data=dataset2)
summary(fit5)
anova_coi=anova(fit5, test="Chisq")
anova_coi
write.csv(anova_coi,"ICE_interference_model_table.csv")


fit_contrast5 <- emmeans::emmeans(fit5, "Treatment", by="Letter.of.the.day", mode="kenward-roger")
fit_contr5 <- contrast(fit_contrast5, method="trt.vs.ctrl")

pheno_contr5 <- as.data.frame(summary(fit_contr5))
pheno_contr5
write.csv(pheno_contr5,"ICE_interference_posthoc_table.csv")

or=exp(pheno_contr5$estimate)
or

#convert p-values to stars for plot
sig2=ifelse(pheno_contr5$p.value<0.001,"***",ifelse(pheno_contr5$p.value<0.01,"**",ifelse(pheno_contr5$p.value<0.05,"*","")))

#Interference figure
pdf("ICECOI.pdf")
COI_figure=ggplot(aes(y=Interference,x=Letter.of.the.day, col=Treatment,label=Num_moms),data=dataset2)+ylab("COI")+ggtitle("Interference")+theme_base()+ylim(-1,1)
COI_figure=COI_figure+scale_colour_manual(values=c("blue","green","orange","red"))+geom_point(size=1.5)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
COI_figure=COI_figure+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+
  annotate(geom="text", x=1, y=0.25, label=sig2[1],size=5)+annotate(geom="text", x=2, y=0.25, label=sig2[2],size=5)+annotate(geom="text", x=3, y=0.25, label=sig2[3],size=5)+annotate(geom="text", x=4, y=0.25, label=sig2[4],size=5)
COI_figure
dev.off()
#box plot

COI_box=ggplot(aes(y=Interference,x=Letter.of.the.day, col=Treatment,label=),data=dataset2)+scale_colour_manual(values=c("blue", "orange", "red","green"))+geom_boxplot()+ylab("COI")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
COI_box=COI_box+ggtitle("interference")+ylim(-1,1)
COI_box=COI_box+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
COI_box

#####
#recode phenotypes as characters
#ICE_male$cut=as.character(ICE_male$cut)
ICE_male$Scalloped=as.character(ICE_male$Scalloped)
ICE_male$Yellow=as.character(ICE_male$Yellow)
ICE_male$Sepia=as.character(ICE_male$Sepia)
ICE_male$co_class=ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow=="Absent" & ICE_male$Scalloped=="Absent" & ICE_male$Yellow!=ICE_male$Sepia,"A", 
                   ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow=="Present" & ICE_male$Scalloped=="Present" & ICE_male$Yellow!=ICE_male$Sepia,"B",
                          ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow=="Absent" & ICE_male$Scalloped=="Absent" & ICE_male$Yellow==ICE_male$Sepia,"A",
                                 ifelse(ICE_male$Scalloped==ICE_male$Yellow & ICE_male$Yellow=="Present" & ICE_male$Scalloped=="Present" & ICE_male$Yellow==ICE_male$Sepia,"B",
                                        ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Scalloped=="Absent","A",
                                               ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow==ICE_male$Sepia & ICE_male$Scalloped=="Present","B",
                                                      ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia & ICE_male$Yellow=="Absent","B",
                                                             ifelse(ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Yellow!=ICE_male$Sepia & ICE_male$Yellow=="Present","A",
                                                                    "error_ctScalloped"))))))))
#check for errors
ICE_male[ICE_male$co_class=="error_ctsd3",]
A=subset(ICE_male,ICE_male$co_class=="A")
B=subset(ICE_male,ICE_male$co_class=="B")
#add columns for counting
ICE_male$NCO_1=ifelse(ICE_male$co_class=="non_CO_1_ctsd1",1,0)
ICE_male$NCO_2=ifelse(ICE_male$co_class=="non_CO_1_ctsd2",1,0)
ICE_male$NCO_3=ifelse(ICE_male$co_class=="non_CO_2_ctsd3",1,0)
ICE_male$NCO_4=ifelse(ICE_male$co_class=="non_CO_2_ctsd4",1,0)
ICE_male$SCO_1=ifelse(ICE_male$co_class=="single_CO_1_ctsd1",1,0)
ICE_male$SCO_2=ifelse(ICE_male$co_class=="single_CO_1_ctsd2",1,0)
ICE_male$SCO_3=ifelse(ICE_male$co_class=="single_CO_2_ctsd3",1,0)
ICE_male$SCO_4=ifelse(ICE_male$co_class=="single_CO_2_ctsd4",1,0)
ICE_male$error=ifelse(ICE_male$co_class=="error_ctsd3",1,0)



sco_count_1_ctsd1=sum(ICE_male$SCO_1)
sco_count_1_ctsd2=sum(ICE_male$SCO_2)
sco_count_2_ctsd3=sum(ICE_male$SCO_3)
sco_count_2_ctsd4=sum(ICE_male$SCO_4)
nco_count_1_ctsd1=sum(ICE_male$NCO_1)
nco_count_1_ctsd2=sum(ICE_male$NCO_2)
nco_count_2_ctsd3=sum(ICE_male$NCO_3)
nco_count_2_ctsd4=sum(ICE_male$NCO_4)
error_count_ctsd3=sum(ICE_male$error)


######Seems to be working####

error=subset(ICE_male,ICE_male$Scalloped!=ICE_male$Yellow)
error1=subset(ICE_male,ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Scalloped=="Present")
error2=subset(ICE_male,ICE_male$Scalloped!=ICE_male$Yellow &ICE_male$Scalloped=="Absent")
View(error)
tapply(error$Number.of.males,error$Name,FUN= function(x) length(unique(x)))
error$uniqueID=paste(error$Start.time,error$Vial.number,error$Letter.of.the.day)
length(unique(error$uniqueID))
error$uniqueID=rownames(error)
length(unique(error$uniqueID))
tapply(error$uniqueID,error$Name,FUN= function(x) length(unique(x)))
ICE_male$uniqueID=rownames(ICE_male)
tapply(ICE_male$uniqueID,ICE_male$Name,FUN= function(x) length(unique(x)))
x=tapply(ICE_male$uniqueID,ICE_male$Name,FUN= function(x) length(unique(x)))
y=tapply(error$uniqueID,error$Name,FUN= function(x) length(unique(x)))
y/x
error1=subset(ICE_male,ICE_male$Scalloped!=ICE_male$Yellow & ICE_male$Scalloped=="Present")
error2=subset(ICE_male,ICE_male$Scalloped!=ICE_male$Yellow &ICE_male$Scalloped=="Absent")
tapply(ICE_male$uniqueID,ICE_male$Name,FUN= function(x) length(unique(x)))
tapply(error1$uniqueID,error1$Name,FUN= function(x) length(unique(x)))
tapply(error2$uniqueID,error2$Name,FUN= function(x) length(unique(x)))
tapply(error1$uniqueID,list(error1$Name,error1$Letter.of.the.day),FUN= function(x) length(unique(x)))
tapply(ICE_male$uniqueID,list(ICE_male$Name,ICE_male$Letter.of.the.day),FUN= function(x) length(unique(x)))
q=tapply(ICE_male$uniqueID,list(ICE_male$Name,ICE_male$Letter.of.the.day),FUN= function(x) length(unique(x)))
w=tapply(error1$uniqueID,list(error1$Name,error1$Letter.of.the.day),FUN= function(x) length(unique(x)))
w/q
cutzero=subset(ICE_male,ICE_male$ct==0)
tapply(cutzero$uniqueID,list(cutzero$Initials,cutzero$Day),FUN= function(x) length(unique(x)))
p=tapply(cutzero$uniqueID,list(cutzero$Initials,cutzero$Day),FUN= function(x) length(unique(x)))
w/p
p/q
tapply(cutzero$uniqueID,cutzero$Day,FUN= function(x) length(unique(x)))
a=tapply(cutzero$uniqueID,cutzero$Day,FUN= function(x) length(unique(x)))
c=tapply(ICE_male$uniqueID,ICE_male$Day,FUN= function(x) length(unique(x)))
a/c
sdzero=subset(ICE_male,ICE_male$sd==0)
d=tapply(sdzero$uniqueID,sdzero$Day,FUN= function(x) length(unique(x)))
d/c
yzero=subset(ICE_male,ICE_male$y==0)
e=tapply(yzero$uniqueID,yzero$Day,FUN= function(x) length(unique(x)))
e/c
sezero=subset(ICE_male,ICE_male$se==0)
f=tapply(sezero$uniqueID,sezero$Day,FUN= function(x) length(unique(x)))
f/c




# less biased ------NOT FULLY DONE-------------------------------------------------------


#I subsetted the less biased haplotype groups, so we could reanalyze our previous data. 
yv1=subset(ICE_male,ICE_male$gamete_class=="gt2a")
yv2=subset(ICE_male,ICE_male$gamete_class=="gt4a")
yv3=subset(ICE_male,ICE_male$gamete_class=="gt1b")
yv4=subset(ICE_male,ICE_male$gamete_class=="gt3b")
yv=rbind(yv1,yv2,yv3,yv4)

#Define Crossovers
yv$co_class=ifelse(yv$Scalloped==yv$Yellow & yv$Yellow==yv$Sepia,"non_CO", 
                   ifelse(yv$Scalloped!=yv$Yellow & yv$Yellow==yv$Sepia,"single_CO_1",
                          ifelse(yv$Scalloped==yv$Yellow & yv$Yellow!=yv$Sepia,"single_CO_2",
                                 ifelse(yv$Scalloped!=yv$Yellow & yv$Yellow!=yv$Sepia,"double_CO",
                                        "error"))))

#add columns for counting
yv$NCO=ifelse(yv$co_class=="non_CO",yv$Number.of.males,0)
yv$SCO_1=ifelse(yv$co_class=="single_CO_1",yv$Number.of.males,0)
yv$SCO_2=ifelse(yv$co_class=="single_CO_2",yv$Number.of.males,0)
yv$DCO=ifelse(yv$co_class=="double_CO",yv$Number.of.males,0)

#get rough crossover rate as a sanity check the numbers of crossovers at the intervals 
#should equal to the total crossover rate.
nco_count=sum(yv$Number.of.males*yv$NCO)
sco_count=sum(yv$Number.of.males*yv$SCO_1, na.rm = TRUE)+sum(yv$Number.of.males*yv$SCO_2, na.rm = TRUE)
dco_count=sum(yv$Number.of.males*yv$DCO, na.rm = TRUE)
num_samples=sum(nco_count+sco_count+dco_count, na.rm = TRUE) 


#rate between scalloped and yellow
(sum(yv$Number.of.males*yv$SCO_1,na.rm = TRUE)+sum(yv$Number.of.males*yv$DCO,na.rm = TRUE))/num_samples

#rate between sepia and yellow
(sum(yv$Number.of.males*yv$SCO_2, na.rm = TRUE)+sum(yv$Number.of.males*yv$DCO,na.rm = TRUE))/num_samples

#we defined the crossover groups as if non crossovers will be identified as 0s; this is same as co_class column, but numerical
#single crossovers 1s
#double crossovers 2s
yv$num_co=ifelse(yv$Yellow==yv$Scalloped & yv$Yellow==yv$Sepia,0, 
                 ifelse(yv$sd==yv$Yellow & yv$Yellow!=yv$Sepia,1*yv$Number.of.males, 
                        ifelse(yv$Scalloped!=yv$Yellow & yv$Yellow==yv$Sepia,1*yv$Number.of.males,  
                               ifelse(yv$Scalloped!=yv$Yellow & yv$yellow!=yv$Sepia,2*yv$Number.of.males, 
                                      NA))))

#This is for summarizing our data
yv$male=c(yv$Number.of.males)

yv_merged <- merge(yv, ICE_bc, by.x = "newVial", by.y = "Vial.number", all=T)


#summarize long form data
dataset=summaryBy(male+num_co+SCO_1+SCO_2+DCO~F1.Vial+Day+Treatment,data=yv_merged, FUN=mean,na.rm=T)

#one more merge
dataset2=merge(yv_female_short,dataset, by=c("F1.Vial","Day","Treatment"))


#poisson regression, similar to a t-test for count data
fit=glm(fecundity~Treatment*Day,data=dataset2,family=quasipoisson)
#summary(fit)
anova_fec=anova(fit, test="Chisq")
write.csv(anova_fec,"yv_fecundity_model_table.csv")

fit_contrast <- emmeans::emmeans(fit, "Treatment", by="Day", mode="kenward-roger")
fit_contr <- contrast(fit_contrast, method="trt.vs.ctrl")

pheno_contr <- as.data.frame(summary(fit_contr))
pheno_contr
write.csv(pheno_contr,"yv_fecundity_posthoc_table.csv")

#convert p-values to stars for plot
sig=ifelse(pheno_contr$p.value<0.001,"***",ifelse(pheno_contr$p.value<0.01,"**",ifelse(pheno_contr$p.value<0.05,"*","")))


#boxplots for fecundity analysis
Fecund_figure_yv=ggplot(dataset2, aes(x=Day, y=fecundity, col=Treatment)) + theme_base()+ylab("# Progeny per mom")+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  geom_boxplot()+scale_x_discrete(name="Days", labels=c("6","7","8","9","10"))+ylim(0,50)+
  annotate(geom="text", x=1, y=38, label=sig[1],size=10)+annotate(geom="text", x=2, y=38, label=sig[2],size=10)+annotate(geom="text", x=3, y=38, label=sig[3],size=10)+annotate(geom="text", x=4, y=38, label=sig[4],size=10)+
  scale_color_manual(values = c("blue","red"))
Fecund_figure_yv


#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$fecundity,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$fecundity,dataset2$Day,mean,na.rm=T)

#remove vials with few progeny
dataset2=dataset2[dataset2$male.sum>=10,]
dataset2=dataset2[dataset2$Numbfemales.sum>=10,]

#sum of crossovers in intervals 1 & 2 
dataset2$num_CO_1=dataset2$SCO_1.sum+dataset2$DCO.sum
dataset2$num_CO_2=dataset2$SCO_2.sum+dataset2$DCO.sum

#sum of non-crossovers in intervals 1 & 2
dataset2$num_NCO_1=dataset2$male.sum-(dataset2$SCO_1.sum+dataset2$DCO.sum)
dataset2$num_NCO_2=dataset2$male.sum-(dataset2$SCO_2.sum+dataset2$DCO.sum)

#total recombination rate
dataset2$rec_rate_total=dataset2$num_co.sum/dataset2$male.sum
dataset2$rec_rate_ysd=(dataset2$SCO_1.sum+dataset2$DCO.sum)/dataset2$male.sum
dataset2$rec_rate_yse=(dataset2$SCO_2.sum+dataset2$DCO.sum)/dataset2$male.sum
mean(dataset2$rec_rate_total)
mean(dataset2$rec_rate_ysd)
mean(dataset2$rec_rate_yse)

#Kosambi correction
dataset2$kosambi_rec_rate_ysd=((2.71^(4*dataset2$rec_rate_ysd)-1)/2*(2.71^(-4*dataset2$rec_rate_ysd)+1))/10
dataset2$kosambi_rec_rate_yse=((2.71^(4*dataset2$rec_rate_yse)-1)/2*(2.71^(-4*dataset2$rec_rate_yse)+1))/10



###Recombination figures

#Recomb_figure_yellow_scalloped
Recomb_figure=ggplot(aes(y=rec_rate_ysd,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue", "red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-sd Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure

#Recomb_figure_yellow_sepia
Recomb_figure=ggplot(aes(y=rec_rate_yse,x=Day, col=as.factor(Treatment),label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_point(alpha=0.6,size=3)+ylab("% recombination")+theme_base()+stat_summary(fun = median, geom="line",aes(group=Treatment),size=2)+scale_x_discrete(name="Day",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-se Recombination rate vs. Days post-mating") +geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=3)
Recomb_figure


#get averages for text; These are the numbers we put in the paper!
tapply(dataset2$rec_rate_total,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_total,dataset2$Day,mean,na.rm=T)

d=subset(dataset2, dataset2$Day=="Y")
tapply(d$rec_rate_total,d$Treatment,mean,na.rm=T)
tapply(d$rec_rate_ysd,d$Treatment,mean,na.rm=T)
tapply(d$rec_rate_yse,d$Treatment,mean,na.rm=T)
e=subset(dataset2,dataset2$Day=="V")
tapply(e$rec_rate_yse,e$Treatment,mean,na.rm=T)
tapply(e$rec_rate_yse,e$Day,mean,na.rm=T)

tapply(dataset2$kosambi_rec_rate_ysd,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_ysd,dataset2$Day,mean,na.rm=T)

tapply(dataset2$kosambi_rec_rate_yse,dataset2$Treatment,mean,na.rm=T)
tapply(dataset2$rec_rate_yse,dataset2$Day,mean,na.rm=T)

###Odds ratios
#SD-Y REGION

#logistic regression, similar to a t-test for count data
fit3=glmer(cbind(num_CO_1,num_NCO_1)~(1|F1.Vial)+Treatment*Day,data=dataset2,
           family=binomial(link="logit"))
#coefs=coef(fit3)
anova_sdy=Anova(fit3,test="Chisq")
write.csv(anova_sdy,"yv_recrate_sd-y_model_table.csv")

fit_contrast3 <- emmeans::emmeans(fit3, "Treatment", by="Day", mode="kenward-roger")
fit_contr3 <- contrast(fit_contrast3, method="trt.vs.ctrl")

pheno_contr3 <- as.data.frame(summary(fit_contr3))
pheno_contr3
write.csv(pheno_contr3,"yv_recrate_sd-y_posthoc_table.csv")

#we can extract the odd ratios and SE from the posthoc table, which is a LOT cleaner and does not require additional model fits for each time point!
y_sd_or=exp(pheno_contr3$estimate)
y_sd_error=(pheno_contr3$SE)

#Y-SE REGION

#logistic regression, similar to a t-test for count data
fit4=glmer(cbind(num_CO_2,num_NCO_2)~(1|F1.Vial)+Treatment*Day,data=dataset2,
           family=binomial(link="logit"))
#summary(fit4)
#coefs=coef(fit4)
#coefs
anova_yse=Anova(fit4,test="Chisq")
write.csv(anova_yse,"yv_recrate_y-se_model_table.csv")

fit_contrast4 <- emmeans::emmeans(fit4, "Treatment", by="Day", mode="kenward-roger")
fit_contr4 <- contrast(fit_contrast4, method="trt.vs.ctrl")

pheno_contr4 <- as.data.frame(summary(fit_contr4))
pheno_contr4
write.csv(pheno_contr4,"yv_recrate_y-se_posthoc_table.csv")

y_se_or=exp(pheno_contr4$estimate)
y_se_error=(pheno_contr4$SE)

#setup data frame
odds_ratios=cbind(y_sd_or,y_se_or)
colnames(odds_ratios)=c("sd-y","y-se")
rownames(odds_ratios)=c("6","7","8","9","10")

#melt dataframe into long form
odds_ratios=melt(odds_ratios)
colnames(odds_ratios)=c("Day","Interval","OR")
odds_ratios

#add in standard error
SE=c(y_sd_error,y_se_error)
x=cbind(odds_ratios,SE)

#convert p-values to stars for plot
ysd_sig=ifelse(pheno_contr3$p.value<0.001,"***",ifelse(pheno_contr3$p.value<0.01,"**",ifelse(pheno_contr3$p.value<0.05,"*","")))
yse_sig=ifelse(pheno_contr4$p.value<0.001,"***",ifelse(pheno_contr4$p.value<0.01,"**",ifelse(pheno_contr4$p.value<0.05,"*","")))

#add significance to table
sig=c(ysd_sig,yse_sig)
y=cbind(x,sig)
odds_ratios=y

#Odds ratio plot
pdf("yv_odds_ratio.pdf")
odds_figure_yv=ggplot(aes(y=OR,x=Day, col=Interval,group=Interval),data=odds_ratios)+scale_colour_manual(values=c("#f1a340", "#998ec3"))+
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.50,1.76)+
  geom_line(size=2)+geom_errorbar(aes(ymin=OR-SE,ymax=OR+SE))+ggtitle("Experiment 5")+theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=6, y=1.65, label=ysd_sig[1],color="#f1a340",size=10)+annotate(geom="text", x=7, y=1.65, label=ysd_sig[2],color="#f1a340",size=10)+annotate(geom="text", x=8, y=1.65, label=ysd_sig[3],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=ysd_sig[4],color="#f1a340",size=10)+annotate(geom="text", x=9, y=1.65, label=yse_sig[4],color="#f1a340",size=10)+
  annotate(geom="text", x=6, y=1.6, label=yse_sig[1],color="#998ec3",size=10)+annotate(geom="text", x=7, y=1.6, label=yse_sig[2],color="#998ec3",size=10)+annotate(geom="text", x=8, y=1.6, label=yse_sig[3],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)+annotate(geom="text", x=9, y=1.6, label=yse_sig[4],color="#998ec3",size=10)
odds_figure_yv
dev.off()


#boxplot for total recombination rate
pdf("yv_sdyboxplotonly.pdf")
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_ysd,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("sd-y Interval")+ylim(0,0.6)
Recomb_figure=Recomb_figure+annotate(geom="text", x=1, y=0.58, label=ysd_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=ysd_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=ysd_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=ysd_sig[4],size=10)
Recomb_figure

dev.off()


pdf("yv_yseboxplotonly.pdf")
Recomb_figure=ggplot(aes(y=kosambi_rec_rate_yse,x=Day, col=Treatment,label=male.sum),data=dataset2)+scale_colour_manual(values=c("blue","red"))+geom_boxplot()+ylab("% recombination")+theme_base()+scale_x_discrete(name="Days post-mating",labels=c("6","7","8","9","10"))
Recomb_figure=Recomb_figure+ggtitle("y-se Interval")+ylim(0.1,0.7)
Recomb_figure=Recomb_figure+annotate(geom="text", x=1, y=0.58, label=yse_sig[1],size=10)+annotate(geom="text", x=2, y=0.58, label=yse_sig[2],size=10)+annotate(geom="text", x=3, y=0.58, label=yse_sig[3],size=10)+annotate(geom="text", x=4, y=0.58, label=yse_sig[4],size=10)
Recomb_figure

dev.off()

