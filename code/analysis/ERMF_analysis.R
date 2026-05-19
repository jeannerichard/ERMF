# Written by Ben Meuleman
# Edited by Jeanne Richard
# Last update: 2026-05-19

library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(foreach)
library(visreg)
library(lcmm)
library(plyr)
library(dplyr)
library(DHARMa)
library(clubSandwich)
library(MuMIn)
library(permuco)
library(effectsize)
library(tidyverse)
library(xtable)
library(forcats)
library(cowplot)
library(ggpubr)

raw_dir <- here::here("rawdata")
fig_dir <- here::here("figures")
data_file <- here::here(raw_dir, "ERMF_dataframe.csv")

#############################################
## FINAL ANALYSIS
#############################################

## DATA SUBSETTING, AGGREGATION AND WIDE FORMAT
subset <- data[,c("site","task_group","id","trial","MFcoeff","MFindex","video_cond","age","gender","BMI_V1","rscale","rating")]
head(subset) ; dim(subset)

meanlong <- aggregate(cbind(rating,MFcoeff)~age+gender+site+task_group+BMI_V1+video_cond+rscale+id,data=subset,FUN=mean,na.rm=TRUE)
head(meanlong) ; dim(meanlong)
str(meanlong)
xtabs(~id,data=meanlong) ; range(xtabs(~id,data=meanlong))

meanshort <- foreach(i=levels(meanlong$id),.combine="rbind",.errorhandling="remove") %do% { 
  subs <- meanlong[which(meanlong$id==i),]
  subs <- rbind(data.frame(subs[subs$video_cond=="positif" & subs$rscale=="amu_resp",]),
                data.frame(subs[subs$video_cond=="mixte" & subs$rscale=="rep_resp",]),
                data.frame(subs[subs$video_cond=="negatif" & subs$rscale=="rep_resp",]))
  subs
}
head(meanshort)
meanshort$rating <- ifelse(meanshort$video_cond=="mixte",meanshort$MFcoeff,meanshort$rating)
meanshort$rating[meanshort$video_cond=="positif"] <- scale(meanshort$rating[meanshort$video_cond=="positif"])
meanshort$rating[meanshort$video_cond=="mixte"] <- scale(meanshort$rating[meanshort$video_cond=="mixte"])
meanshort$rating[meanshort$video_cond=="negatif"] <- scale(meanshort$rating[meanshort$video_cond=="negatif"])
head(meanshort) ; dim(meanshort)

meanwide <- foreach(i=unique(meanlong$id),.combine="rbind") %do% {
  sub <- rbind(meanlong[meanlong$id==i,])
  labels <- c(paste(as.character(sub$video_cond),as.character(sub$rscale),"rating",sep="_"),paste(as.character(sub$video_cond),as.character(sub$rscale),"MFcoeff",sep="_"))
  widetmp <- data.frame(sub[1,c("id","age","gender","site","task_group","BMI_V1")],rbind(sub$rating),rbind(sub$MFcoeff))
  colnames(widetmp)[7:18] <- labels
  widetmp
}
head(meanwide) ; dim(meanwide)
str(meanwide)

#############################################
## RAW DESCRIPTIVES
#############################################

desc <- rbind(aggregate(cbind(rating,MFcoeff)~task_group*video_cond*rscale,data=meanlong,FUN=mean),
aggregate(cbind(rating,MFcoeff)~task_group*video_cond*rscale,data=meanlong,FUN=sd))

#write.csv(desc,"desc_tmp.csv")



#############################################
## AMUSEMENT AND REPULSION ANALYSIS
#############################################

## BASIC DESIGN EFFECTS
#--------------------------------------------

### MODELLING
AIC(lmer(rating~task_group*video_cond*rscale+(1|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(1+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(0+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(1+rscale|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(0+rscale|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(1+rscale+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(0+rscale+video_cond|id),data=meanlong)) # WINNER
AIC(lmer(rating~task_group*video_cond*rscale+(1+rscale+video_cond||id),data=meanlong))
AIC(lmer(rating~task_group*video_cond*rscale+(0+rscale+video_cond||id),data=meanlong))

model <- lmer(rating~task_group*video_cond*rscale+(0+rscale+video_cond|id),data=meanlong)
model_scaled <- lmer(scale(rating)~task_group*video_cond*rscale+(0+rscale+video_cond|id),data=meanlong)
summary(model)
anova(model,type=2)
eta_squared(anova(model,type=2),alternative="two.sided")

### FOLLOW-UP TESTS
joint_tests(model,by="video_cond",lmer.df="satterthwaite")
joint_tests(model,by="task_group",lmer.df="satterthwaite")

rg <- ref_grid(model,specs="rscale",by=c("video_cond","task_group"),lmer.df="satterthwaite")
emmeans(rg,specs="rscale",by=c("video_cond","task_group"),infer=TRUE,adjust="none")
tabS2 <- data.frame(pairs(emmeans(rg,specs="rscale",by=c("video_cond","task_group"),adjust="none"))) # Table S2

rg_scaled <- ref_grid(model_scaled,specs="rscale",by=c("video_cond","task_group"),lmer.df="satterthwaite")
tabS2_scaled <- data.frame(pairs(emmeans(rg_scaled,specs="rscale",by=c("video_cond","task_group"),adjust="none"))) # Table S2_scaled
tabS2_dz <- tabS2_scaled$estimate
S2 <- tabS2 |> 
  mutate(dz = tabS2_dz) |> 
  relocate(task_group, .before = video_cond) |> 
  relocate(contrast, .after = video_cond) |> 
  mutate(task_group = fct_recode(task_group, "Look" = "evalemo", "Regulate" = "regulemo")) |> 
  mutate(video_cond = fct_recode(video_cond, "Positive" = "positif", "Negative" = "negatif", "Mixed" = "mixte")) |> 
  mutate(video_cond = fct_relevel(video_cond, "Positive")) |> 
  mutate(contrast = fct_recode(contrast, "Amusement - Repulsion" = "amu_resp - rep_resp")) |> 
  dplyr::rename(Group = task_group, Video = video_cond, "Scale contrast" = contrast, D = estimate, DF = df, t = t.ratio, p = p.value) |> 
  arrange(Group, Video)

print(xtable(S2), include.rownames=FALSE)


rg <- ref_grid(model,specs=c("rscale","video_cond"),by=c("task_group"),lmer.df="satterthwaite")
emmeans(rg,specs=c("rscale","video_cond"),by=c("task_group"),infer=TRUE,adjust="none")
pairs(emmeans(rg,specs=c("rscale","video_cond"),by=c("task_group"),adjust="none"))

rg <- ref_grid(model,specs=c("rscale"),by=c("video_cond"),lmer.df="satterthwaite")
emmeans(rg,specs=c("rscale"),by=c("video_cond"),infer=TRUE,adjust="none")
pairs(emmeans(rg,specs=c("rscale"),by=c("video_cond"),adjust="none"))


### EFFECT PLOTS
clrs <- c("steelblue","goldenrod")
par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="evalemo"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=clrs), points.par=list(pch=4,col=adjustcolor(clrs,alpha=0.4),cex=0.8),fill.par=list(col=adjustcolor(clrs,alpha=0.3)))
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="regulemo"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=clrs), points.par=list(pch=4,col=adjustcolor(clrs,alpha=0.4),cex=0.8),fill.par=list(col=adjustcolor(clrs,alpha=0.3)))

# FIGURE 2
### RELABEL
meanlong2 <- meanlong
meanlong2$video_cond <- factor(meanlong2$video_cond)
meanlong2$video_cond <- relevel(meanlong2$video_cond, ref = "positif")
meanlong2$rscale <- factor(meanlong2$rscale)
meanlong2$task_group <- factor(meanlong2$task_group)

meanlong2$video_cond 
meanlong2$rscale
meanlong2$task_group

meanlong2$rscale2 <- factor(meanlong2$rscale, labels=c("amu_resp" = "Amusement", "rep_resp" = "Repulsion"))
meanlong2$video_cond2 <- factor(meanlong2$video_cond, labels=c("positif" = "Positive", "mixte" = "Mixed", "negatif" = "Negative"))
meanlong2$task_group2 <- factor(meanlong2$task_group, labels=c("evalemo" = "No regulation", "regulemo" = "Regulation"))

model2 <- lmer(rating~task_group2*video_cond2*rscale2+(0+rscale2+video_cond2|id),data=meanlong2)

par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model2,xvar="video_cond2",by="rscale2",cond=list(task_group2="No regulation"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
       line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model2,xvar="video_cond2",by="rscale2",cond=list(task_group2="Regulation"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
       line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))

vg1 <- visreg(model2,
              xvar="video_cond2",
              by="rscale2",
              cond=list(task_group2="No regulation"),
              overlay=TRUE,
              ylab="Emotion rating",
              xlab="Video",
              gg = TRUE,
              line.par=list(), 
              fill.par=list(alpha=0.2),
              points=list(size=2, pch=4)) + 
  ylim(-1, 5) +
  scale_fill_manual(values = c("steelblue","goldenrod")) +
  scale_color_manual(values = c("steelblue","goldenrod")) + 
  labs(color = "", fill = "") + 
  guides(fill = "none") +
  theme_cowplot(font_size = 24)  +
  theme(legend.position="top")
vg2 <- visreg(model2,
              xvar="video_cond2",
              by="rscale2",
              cond=list(task_group2="Regulation"),
              overlay=TRUE,
              ylab="Emotion rating",
              xlab="Video",
              gg = TRUE,
              line.par=list(), 
              fill.par=list(alpha=0.2),
              points=list(size=2, pch=4)) + 
  ylim(-1, 5) +
  scale_fill_manual(values = c("steelblue","goldenrod")) +
  scale_color_manual(values = c("steelblue","goldenrod")) + 
  labs(color = "", fill = "") + 
  guides(fill = "none") +
  theme_cowplot(font_size = 24)  +
  theme(legend.position="top")

figure_2 <- ggarrange(vg1, vg2, labels = c("A", "B"))

png(file= here::here(fig_dir,"figure-2-2.png"),width=40,height=30,units="cm",res=300)
figure_2
dev.off()


### DIAGNOSTICS
r.squaredGLMM(model)

vif(model)

P <- length(fixef(model))
N <- length(residuals(model))
threshold <- qf(0.5,P,N-P)
plot(cooks.distance(model),type="h",col="slategray3",ylab="Cook's distance")
abline(h=threshold,col="firebrick3")

influential <- as.numeric(names(which(cooks.distance(model)>threshold)))
influential

Boxplot(ranef(model)$id)

set.seed(1999)
simulateResiduals(model,n=1000,plot=TRUE)


## BMI AND COVARIATES
#--------------------------------------------

### BMI cutoffs
cuts <- c(18,22,26,35)

### CENTERED VARIABLES (INSERT THESE IN MODELS THAT NEED TO COLLAPSE OVER THE CENTERED VARIABLE!)
meanlong$BMI_c <- scale(meanlong$BMI_V1)
meanlong$task_group_c <- meanlong$task_group
meanlong$rscale_c <- meanlong$rscale
contrasts(meanlong$task_group_c) <- contr.sum(levels(meanlong$task_group_c))
contrasts(meanlong$rscale_c) <- contr.sum(levels(meanlong$rscale_c))

### MODELLING
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(1|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(1+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(0+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(1+rscale|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(0+rscale|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(1+rscale+video_cond|id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(0+rscale+video_cond|id),data=meanlong)) # WINNER
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(1+rscale+video_cond||id),data=meanlong))
AIC(lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(0+rscale+video_cond||id),data=meanlong))

model <- lmer(rating~task_group*BMI_V1*video_cond*rscale+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
summary(model)
anova(model,type=2)
eta_squared(anova(model,type=2),alternative="two.sided")

model <- lmer(rating~(task_group+BMI_V1+video_cond+rscale)^3+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
summary(model)
anova(model,type=2)
eta_squared(anova(model,type=2),alternative="two.sided")


### FOLLOW-UP TESTS
rg <- ref_grid(model,at=list(BMI_V1=cuts),lmer.df="satterthwaite",nuisance=c("gender","site"))
rg
joint_tests(rg,by=c("BMI_V1"))
joint_tests(rg,by=c("video_cond"))
desc <- data.frame(emmeans(rg,specs="rscale",by=c("video_cond","BMI_V1"),infer=TRUE,adjust="none"))
#write.csv2(data.frame(emmeans(rg,specs="rscale",by=c("video_cond","BMI_V1"),infer=TRUE,adjust="none")),"temp_desc.csv") # EXPORT MODEL-BASED DESCRIPTIVES

tabS3_a <- data.frame(desc[,1:4],SD=desc$SE*sqrt(desc$df))
  


rg <- ref_grid(model,specs="rscale",by=c("video_cond","task_group"),lmer.df="satterthwaite")
emmeans(rg,specs="rscale",by=c("video_cond","task_group"),infer=TRUE,adjust="none")
data.frame(pairs(emmeans(rg,specs="rscale",by=c("video_cond","task_group"),adjust="none")))

emmeans(rg,specs=c("video_cond","rscale"),by="BMI_V1",infer=TRUE,adjust="none")
data.frame(pairs(emmeans(rg,specs=c("rscale"),by=c("video_cond","BMI_V1")),adjust="none"))

emtrends(model,var="BMI_V1",specs=c("video_cond","rscale"),infer=TRUE,lmer.df = "satterthwaite")
pairs(emtrends(model,var="BMI_V1",specs=c("video_cond","rscale")),lmer.df = "satterthwaite",adjust="none")

emmeans(model,specs=c("site"),infer=TRUE,adjust="none")
data.frame(pairs(emmeans(model,specs=c("site"),infer=TRUE,adjust="none"),adjust="none"))


### MODEL PLOTS
par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="evalemo"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="regulemo"),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))

par(mfrow=c(2,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_V1=18),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_V1=22),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_V1=46),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_V1=35),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))

par(mfrow=c(1,2))
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="evalemo"),overlay=TRUE,ylim=c(-1,5))
visreg(model,xvar="video_cond",by="rscale",cond=list(task_group="regulemo"),overlay=TRUE,ylim=c(-1,5))

par(mfrow=c(1,2))
visreg(model,xvar="video_cond",by="BMI_V1",cond=list(rscale="amu_resp"),breaks=c(17,44),overlay=TRUE,ylim=c(-3,5),ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("bisque3","deeppink3")), points.par=list(pch=4,col=c("bisque3","deeppink3"),cex=0.8),fill.par=list(col=adjustcolor(c("bisque3","deeppink3"),alpha=0.2)))
visreg(model,xvar="video_cond",by="BMI_V1",cond=list(rscale="rep_resp"),breaks=c(17,44),overlay=TRUE,ylim=c(-3,5),ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("bisque3","deeppink3")), points.par=list(pch=4,col=c("bisque3","deeppink3"),cex=0.8),fill.par=list(col=adjustcolor(c("bisque3","deeppink3"),alpha=0.2)))

par(mfrow=c(1,2))
visreg(model,xvar="BMI_V1",by="video_cond",cond=list(rscale="amu_resp"),overlay=TRUE,ylim=c(-3,5))
visreg(model,xvar="BMI_V1",by="video_cond",cond=list(rscale="rep_resp"),overlay=TRUE,ylim=c(-3,5))


### DIAGNOSTICS
r.squaredGLMM(model)

vif(model)

P <- length(fixef(model))
N <- length(residuals(model))
threshold <- qf(0.5,P,N-P)
plot(cooks.distance(model),type="h",col="slategray3",ylab="Cook's distance")
abline(h=threshold,col="firebrick3")

influential <- as.numeric(names(which(cooks.distance(model)>threshold)))
influential

Boxplot(ranef(model)$id)

set.seed(1999)
simulateResiduals(model,n=1000,plot=TRUE)


## SENSITIVITY CHECK
#---------------------------------------------

### REPEATED MEASURES MANOVA
design <- data.frame(rscale=factor(c("amu","amu","amu","rep","rep","rep")),video_cond=factor(c("pos","mix","neg","pos","mix","neg"),levels=c("pos","mix","neg")))
model <- lm(cbind(positif_amu_resp_rating,mixte_amu_resp_rating,negatif_amu_resp_rating,positif_rep_resp_rating,mixte_rep_resp_rating,negatif_rep_resp_rating)~BMI_V1*task_group+age+gender+site,data=meanwide)
Anova(model,type=2,test="Pillai",idesign=~rscale*video_cond,idata=design)

### HETEROSCEDASTIC MLM
model <- lmer(rating~(task_group+BMI_V1+video_cond+rscale)^2+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
model <- lmer(rating~(task_group+BMI_V1+video_cond+rscale)^3+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
model <- lmer(rating~(task_group+BMI_V1+video_cond+rscale)^4+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
summary(model)
coef_test(model, vcov = "CR2")

Wald_test(model, constraints = constrain_zero(c("video_condmixte:rscalerep_resp","video_condnegatif:rscalerep_resp")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("task_groupregulemo:video_condmixte:rscalerep_resp","task_groupregulemo:video_condnegatif:rscalerep_resp")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("BMI_V1:video_condmixte:rscalerep_resp","BMI_V1:video_condnegatif:rscalerep_resp")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("task_groupregulemo:BMI_V1:video_condmixte:rscalerep_resp","task_groupregulemo:BMI_V1:video_condnegatif:rscalerep_resp")), vcov = "CR2", test = "HTZ")

### PERMUTATION REPEATED MEASURES ANOVA
set.seed(1985)
model <- aovperm(rating~(task_group+BMI_V1+video_cond+rscale)^2+age+gender+site+ Error(id/(video_cond*rscale)),data=meanlong, method="Rde_kheradPajouh_renaud") ; model
model <- aovperm(rating~(task_group+BMI_V1+video_cond+rscale)^3+age+gender+site+ Error(id/(video_cond*rscale)),data=meanlong, method="Rde_kheradPajouh_renaud") ; model
model <- aovperm(rating~(task_group+BMI_V1+video_cond+rscale)^4+age+gender+site+ Error(id/(video_cond*rscale)),data=meanlong, method="Rde_kheradPajouh_renaud") ; model

### POLYNOMIAL BMI
model <- lmer(rating~(task_group_c+poly(BMI_c,2)+video_cond+rscale)^3+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
model <- lmer(rating~(task_group_c+poly(BMI_c,2)+video_cond+rscale)^4+age+gender+site+(0+rscale+video_cond|id),data=meanlong)
anova(model,type=2)

par(mfrow=c(1,2))
visreg(model,xvar="BMI_c",by="video_cond",cond=list(rscale="amu_resp"),overlay=TRUE,ylim=c(-3,5))
visreg(model,xvar="BMI_c",by="video_cond",cond=list(rscale="rep_resp"),overlay=TRUE,ylim=c(-3,5))

par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_c=-1.4),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_c=2.2),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))


#############################################
## MIXED EMOTIONS ANALYSIS
#############################################

## BASIC DESIGN EFFECTS
#--------------------------------------------

### MODELLING
model <- lmer(MFcoeff~task_group*video_cond+(1|id),data=meanlong,subset=rscale=="amu_resp")
summary(model)
anova(model,type=2)
eta_squared(anova(model,type=2),alternative="two.sided")

### FOLLOW-UP TESTS
emmeans(model,specs="video_cond",infer=TRUE,adjust="none",lmer.df="satterthwaite")
tabS4 <- data.frame(pairs(emmeans(model,specs="video_cond",lmer.df="satterthwaite"),adjust="none")) # Table S4

# Scale MFcoeff to get standardized mean differences
model2 <- lmer(scale(MFcoeff)~task_group*video_cond+(1|id),data=meanlong,subset=rscale=="amu_resp")
tabS4_scaled <- pairs(emmeans(model2,specs="video_cond",lmer.df="satterthwaite"),adjust="none") # Table S4 (estimate = dz)

tabS4_dz <- data.frame(tabS4_scaled)$estimate

S4 <- tabS4 |> 
  mutate(dz = tabS4_dz) |> 
  mutate(contrast = fct_recode(contrast, 
                               "Mixed - Positive" = "mixte - positif",
                               "Mixed - Negative" = "mixte - negatif",
                               "Positive - Negative" = "negatif - positif" ), # reversed
         contrast = fct_relevel(contrast, "Mixed - Positive")) |> 
  arrange(contrast) |> 
  rename("Video contrast" = contrast, D = estimate, DF = df, t = t.ratio, p = p.value)

print(xtable(S4), include.rownames=FALSE)


### EFFECT PLOTS
clrs <- c("burlywood3","darkcyan")
par(cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="task_group",overlay=TRUE,ylim=c(-1,2),ylab="Mixed emotion score",xlab="Video",
  line.par=list(col=clrs), points.par=list(pch=4,col=adjustcolor(clrs,alpha=0.4),cex=0.8),fill.par=list(col=adjustcolor(clrs,alpha=0.3)))


# FIGURE 3
### RELABEL
meanlong2 <- meanlong
meanlong2$video_cond <- factor(meanlong2$video_cond)
meanlong2$video_cond 
meanlong2$video_cond <- relevel(meanlong2$video_cond, ref = "positif")
meanlong2$video_cond 

meanlong2$task_group_c <- factor(meanlong2$task_group_c)
meanlong2$task_group_c

meanlong2$video_cond2 <- factor(meanlong2$video_cond, labels=c("positif" = "Positive", "mixte" = "Mixed", "negatif" = "Negative"))
meanlong2$task_group_c2 <- factor(meanlong2$task_group_c, labels=c("evalemo" = "No regulation", "regulemo" = "Regulation"))

meanlong2_amu <- meanlong2 |> 
  filter(rscale == "amu_resp") # because the same MFscore is calculated on both scales 

model2 <- lmer(MFcoeff~task_group_c2*video_cond2+(1|id),data=meanlong2_amu)


par(mfrow=c(1,1),cex.axis=1.5,cex.lab=1.5)
figure_3 <- visreg(model2,
       xvar="video_cond2",
       by="task_group_c2",
       overlay=TRUE,
       ylab="Mixed feelings score",
       xlab="Video",
       gg = TRUE,
       line.par=list(), 
       fill.par=list(alpha=0.2),
       points=list(size=2, pch=4)) + 
  ylim(-1, 2) +
  scale_fill_manual(values = c("burlywood3","darkcyan")) +
  scale_color_manual(values = c("burlywood3","darkcyan")) + 
  labs(color = "", fill = "") + 
  guides(fill = "none") +
  theme_cowplot(font_size = 24)  +
  theme(legend.position="top")

png(file= here::here(fig_dir,"figure-3.png"),width=30,height=30,units="cm",res=300)
figure_3
dev.off()

### DIAGNOSTICS
r.squaredGLMM(model)

vif(model)

P <- length(coef(model))
N <- length(residuals(model))
threshold <- qf(0.5,P,N-P)
plot(cooks.distance(model),type="h",col="slategray3",ylab="Cook's distance")
abline(h=threshold,col="firebrick3")

influential <- as.numeric(names(which(cooks.distance(model)>threshold)))
influential

Boxplot(ranef(model)$id)

set.seed(1999)
simulateResiduals(model,n=1000,plot=TRUE)


## BMI AND COVARIATES
#--------------------------------------------

### CENTERED VARIABLES (INSERT THESE IN MODELS THAT NEED TO COLLAPSE OVER THE CENTERED VARIABLE!)
meanlong$BMI_c <- scale(meanlong$BMI_V1)
meanlong$task_group_c <- meanlong$task_group
contrasts(meanlong$task_group_c) <- contr.sum(levels(meanlong$task_group_c))

### MODELLING
model <- lmer(MFcoeff~(task_group+video_cond+BMI_V1)^3+age+gender+site+(1|id),data=meanlong,subset=rscale=="amu_resp")
summary(model)
anova(model,type=2)
eta_squared(anova(model,type=2),alternative="two.sided")

### FOLLOW-UP TESTS
rg <- ref_grid(model,at=list(BMI_V1=cuts),lmer.df="satterthwaite")
rg
joint_tests(rg,by=c("BMI_V1"))
joint_tests(rg,by=c("video_cond"),adjust="none")
desc <- data.frame(emmeans(rg,specs="video_cond",by=c("BMI_V1"),infer=TRUE,adjust="none"))
#write.csv2(data.frame(emmeans(rg,specs="video_cond",by=c("BMI_V1"),infer=TRUE,adjust="none")),"temp_desc.csv") # EXPORT MODEL-BASED DESCRIPTIVES
tabS3_b <- data.frame(desc[,1:2],rscale=rep("MFcoeff",12),emmean=desc[,3],SD=desc$SE*sqrt(desc$df))
S3 <- rbind(tabS3_a,tabS3_b) 

emmeans(rg,specs=c("video_cond"),by="BMI_V1",infer=TRUE,adjust="none")
pairs(emmeans(rg,specs=c("video_cond"),by=c("BMI_V1")),adjust="none")

emtrends(model,var="BMI_c",specs=c("video_cond"),infer=TRUE,lmer.df = "satterthwaite")
pairs(emtrends(model,var="BMI_c",specs=c("video_cond")),lmer.df = "satterthwaite",adjust="none")

emmeans(model,specs="video_cond",infer=TRUE,adjust="none",lmer.df="satterthwaite")
pairs(emmeans(model,specs="video_cond",lmer.df="satterthwaite"),adjust="none")


### MODEL PLOTS
par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="task_group",cond=list(BMI_V1=17),overlay=TRUE,ylim=c(-1,2),ylab="Mixed emotion score",xlab="Video",
  line.par=list(col=c("burlywood3","darkcyan")), points.par=list(pch=4,col=c("burlywood3","darkcyan"),cex=0.8),fill.par=list(col=adjustcolor(c("burlywood3","darkcyan"),alpha=0.2)))
visreg(model,xvar="video_cond",by="task_group",,cond=list(BMI_V1=44),overlay=TRUE,ylim=c(-1,2),ylab="Mixed emotion score",xlab="Video",
  line.par=list(col=c("burlywood3","darkcyan")), points.par=list(pch=4,col=c("burlywood3","darkcyan"),cex=0.8),fill.par=list(col=adjustcolor(c("burlywood3","darkcyan"),alpha=0.2)))

par(mfrow=c(1,2))
visreg(model,xvar="BMI_V1",by="video_cond",cond=list(task_group_c="evalemo"),overlay=TRUE,ylim=c(-1,2))
visreg(model,xvar="BMI_V1",by="video_cond",cond=list(task_group_c="regulemo"),overlay=TRUE,ylim=c(-1,2))

### DIAGNOSTICS
r.squaredGLMM(model)

vif(model)

P <- length(coef(model))
N <- length(residuals(model))
threshold <- qf(0.5,P,N-P)
plot(cooks.distance(model),type="h",col="slategray3",ylab="Cook's distance")
abline(h=threshold,col="firebrick3")

influential <- as.numeric(names(which(cooks.distance(model)>threshold)))
influential

Boxplot(ranef(model)$id)

set.seed(1999)
simulateResiduals(model,n=1000,plot=TRUE)


## SENSITIVITY CHECK
#---------------------------------------------

### REPEATED MEASURES MANOVA
design <- data.frame(video_cond=factor(c("pos","mix","neg"),levels=c("pos","mix","neg")))
model <- lm(cbind(positif_amu_resp_MFcoeff,mixte_amu_resp_MFcoeff,negatif_amu_resp_MFcoeff)~BMI_V1*task_group+age+gender+site,data=meanwide)
Anova(model,type=2,test="Pillai",idesign=~video_cond,idata=design)

### HETEROSCEDASTIC MLM
model <- lmer(MFcoeff~(task_group_c+video_cond+BMI_V1)^3+age+gender+site+(1|id),data=meanlong,subset=rscale=="amu_resp")
model <- lmer(MFcoeff~(task_group_c+video_cond+BMI_V1)^2+age+gender+site+(1|id),data=meanlong,subset=rscale=="amu_resp")
summary(model)
coef_test(model, vcov = "CR2")

Wald_test(model, constraints = constrain_zero(c("task_group_c1:video_condmixte","task_group_c1:video_condnegatif")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("task_group_c1:BMI_V1")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("video_condmixte:BMI_V1","video_condnegatif:BMI_V1")), vcov = "CR2", test = "HTZ")
Wald_test(model, constraints = constrain_zero(c("task_group_c1:video_condmixte:BMI_V1","task_group_c1:video_condnegatif:BMI_V1")), vcov = "CR2", test = "HTZ")


### PERMUTATION REPEATED MEASURES ANOVA
set.seed(1985)
model <- aovperm(MFcoeff~(task_group_c+video_cond+BMI_V1)^3+age+gender+site+Error(id/video_cond),data=meanlong[meanlong$rscale=="amu_resp",], method="Rde_kheradPajouh_renaud") ; model
model <- aovperm(MFcoeff~(task_group_c+video_cond+BMI_V1)^2+age+gender+site+Error(id/video_cond),data=meanlong[meanlong$rscale=="amu_resp",], method="Rde_kheradPajouh_renaud") ; model

### ROBUST MLM
model <- lmer(MFcoeff~task_group_c*video_cond+(1|id),data=meanlong,subset=rscale=="amu_resp")
influential <- as.numeric(rownames(model.matrix(model)[c(387,398,446),]))
meanlong[influential,]

model <- lmer(MFcoeff~(task_group_c+video_cond+BMI_V1)^3+age+gender+site+(1|id),data=meanlong[-influential,],subset=rscale=="amu_resp")
summary(model)
anova(model,type=2)

### POLYNOMIAL BMI
model <- lmer(MFcoeff~(task_group_c+video_cond+poly(BMI_V1,2))^3+age+gender+site+(1|id),data=meanlong,subset=rscale=="amu_resp")
anova(model,type=2)

par(mfrow=c(1,2))
visreg(model,xvar="BMI_c",by="video_cond",cond=list(rscale="amu_resp"),overlay=TRUE,ylim=c(-3,5))
visreg(model,xvar="BMI_c",by="video_cond",cond=list(rscale="rep_resp"),overlay=TRUE,ylim=c(-3,5))

par(mfrow=c(1,2),cex.axis=1.5,cex.lab=1.5)
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_c=-1.4),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))
visreg(model,xvar="video_cond",by="rscale",cond=list(BMI_c=2.2),overlay=TRUE,ylim=c(-1,5),ylab="Emotion rating",xlab="Video",
  line.par=list(col=c("steelblue","goldenrod")), points.par=list(pch=4,col=c("steelblue","goldenrod"),cex=0.8),fill.par=list(col=adjustcolor(c("steelblue","goldenrod"),alpha=0.2)))