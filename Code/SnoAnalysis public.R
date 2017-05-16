#Code for Stevens 2017, Ecological Applications
#Jens Stevens, stevensjt@gmail.com

library(lme4) #for lmer(). version 1.1-12
library(bbmle) #for BICtab(). version 1.0.18
library(Kmisc) #for write.cb(). version 0.5.0
library(ggplot2) #for ggplot(). version 2.2.1
library(plyr) #for ddply(). version 1.8.4
library(reshape2) #for melt(). version 1.4.2
library(gridExtra) #for grid.arrange(). version 2.2.1
library(grid) #for textGrob(). version 3.3.2

GetME_PVals=function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs #Everything is significantly different from UB
  return(coefs)
}


####1. Setup data####
d <- read.csv("Data/Stevens_2017_EAP_data.csv")
d$pt.ID <- as.character(d$pt.ID) #Format plot ID to character
d$Week <- factor(d$Week) #Format week to factor
d$BurnSev <- factor(d$BurnSev) #Format burn severity to factor
d$C.Pos <- factor(d$C.Pos,levels=c("","U","E","O")) #Reorder factor levels in canopy position

####2. Mixed model####

#2a. Compare models, pick best, generate Table 2, Table S2
m0 <- lmer(Depth~ (1|Week/Site),data=d) #Null model
m1a <- #Aspect effect: SW slopes have less
  lmer(Depth~ AspectCat + (1|Week/Site),data=d) 
m1b <- #Burn effect: Significant negative effect of severity on depth.
  lmer(Depth~ BurnSev + (1|Week/Site),data=d) 
m1c <- #Overhead canopy effect: Significant negative effect of canopy on depth
  lmer(Depth~ C.Pos + (1|Week/Site),data=d) 
m2a <- #Burn effect controlling for C.Pos
  lmer(Depth~ BurnSev + C.Pos + (1|Week/Site),data=d)
m2b <- #Burn effect controlling for C.Pos as a random effect
  lmer(Depth~ BurnSev + (1|Week/Site) + (1|C.Pos) ,data=d) 
m2c <- #Overhead canopy effect controlling for burn severity
  lmer(Depth~ C.Pos + (1|Week/Site) + (1|BurnSev),data=d) 
m2d <- #Burn effect controlling for Aspect
  lmer(Depth~ BurnSev + Aspect + (1|Week/Site),data=d) 
m3b <- #Burn effect controlling for C.Pos and Aspect as random effects
  lmer(Depth~ BurnSev + (1|Week/Site) + (1|C.Pos) +(1|AspectCat),data=d) 
m3c <- #Canopy effect controlling for burn severity and Aspect as random effects
  lmer(Depth~ C.Pos + (1|Week/Site) + (1|BurnSev) +(1|AspectCat),data=d) 
m4b <- #Interaction between burn severity and aspect
  lmer(Depth~ BurnSev*AspectCat + (1|Week/Site),data=d)
m4c <- #Interaction between canopy position and aspect
  lmer(Depth~ C.Pos*AspectCat + (1|Week/Site),data=d) 
m.full <- #Full model with all individual terms and no interactions
  lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Site/Week),data=d)
m.inter <- #Full model with interaction between canopy and burn severity. Not fully crossed, can't include aspect.
  lmer(Depth~ C.Pos * BurnSev + AspectCat + (1|Week/Site),data=d)

BICt <- BICtab(m0,m1a,m1b,m1c,m2a,m2b,m2c,m2d,m3b,m3c,m4b,m4c,m.full,m.inter)# m.full best; 

#Calculate R2 to append to Table S2
m.list <- list(m.full,m3b,m3c,m.inter,m4b,m2d,m2a,m2b,m2c,m4c,m1b,m1a,m0,m1c)#Arrange in order of BITc
#Copy and run code from https://github.com/jslefche/rsquared.glmm/blob/master/rsquaredglmm.R
rsquares <- rsquared.glmm(m.list)
#write.cb(cbind(attr(BICt,"row.names"),list2df(BICt),round(rsquares$Marginal,2),rsquares$Conditional)) #Copy for TableS2

#Get the coefficients from m.full
Table2 <- GetME_PVals(m.full) #All the trends hold up in the best model, everything is still pretty significant.
Table2 <- round(Table2,3)
#write.cb(Table2) #Copy for Table 2

#2b: Make parameter distribution figure (Figure S3)
d.full <- d[-which(is.na(d$Depth)),] #Remove NA's. ONLY DO ONCE
for(s in 1:1000){ #Subset data by half, run model 1000 times
  d.sub <- #Subset to half the data
    d.full[-sample(1:nrow(d.full), nrow(d.full)/2, replace=FALSE),]
  m.full <- #Full model with all individual terms and no interactions
    lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Site/Week),data=d.sub)
  Table2 <- GetME_PVals(m.full) 
  Table2 <- round(Table2,3)
  if(s==1){ #Optional if subsetting
    param.dist=data.frame(Sim=1,Var=rownames(Table2),Estimate=Table2$Estimate)} else 
  { #Fill out parameter table
    param.dist[((s-1)*nrow(Table2)+1):(s*nrow(Table2)+nrow(Table2)),"Sim"] <- s
    param.dist[((s-1)*nrow(Table2)+1):(s*nrow(Table2)+nrow(Table2)),"Var"] <- rownames(Table2)
    param.dist[((s-1)*nrow(Table2)+1):(s*nrow(Table2)+nrow(Table2)),"Estimate"] <- Table2$Estimate
  }
  print(s) 
} #End subsetting loop

levels(param.dist$Var) <- #Rename factor levels for plot
  c("Intercept","NE Aspect", "SW Aspect", "Severity Class 1"," Severity Class 2", "Severity Class 3", "Severity Class 4", "Canopy Condition: Edge", "Canopy Condition: Open")
ggplot(param.dist) +
  geom_density(aes(x=Estimate)) +
  geom_vline(xintercept=0,linetype="longdash") +
  facet_wrap(~Var,ncol=3) +
  theme_bw()

#dev.copy2pdf(file="Figures/FigS3.pdf",width=7,height=7)

#2c:Get the coefficients from m.full for a single fire (Table S3)
d.site <- d[d$Site=="Angora",] #Set to a single fire
m.full <- lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Week),data=d.site) #Run model
TableS3 <- GetME_PVals(m.full) #All the trends hold up for each fire, everything is still pretty significant.
TableS3 <- round(TableS3,3)
#write.cb(TableS3) #Copy for Table S3 (do for each fire)

#2d: Summarize actual mean snow depths across scenarios (Table S4)
d.S4 <- d
d.S4$Depth.Angora[d.S4$Site=="Angora"] <- d.S4[d.S4$Site=="Angora","Depth"]
d.S4$Depth.Reading[d.S4$Site=="Reading"] <- d.S4[d.S4$Site=="Reading","Depth"]
d.S4$Depth.Showers[d.S4$Site=="Showers"] <- d.S4[d.S4$Site=="Showers","Depth"]
TableS4 <- ddply(d.S4,.(Week,BurnSev,C.Pos),summarize,
                 MeanSnowDepth.Angora=round(mean(Depth.Angora,na.rm=T),2),
                 sdSnowDepth.Angora=round(sd(Depth.Angora,na.rm=T),2),
                 MeanSnowDepth.Reading=round(mean(Depth.Reading,na.rm=T),2),
                 sdSnowDepth.Reading=round(sd(Depth.Reading,na.rm=T),2),
                 MeanSnowDepth.Showers=round(mean(Depth.Showers,na.rm=T),2),
                 sdSnowDepth.Showers=round(sd(Depth.Showers,na.rm=T),2))
#write.cb(TableS4)

#Multi-panel figure
#For clearer visualization, take the 3 points at Showers that were the only representative from their canopy class-burn severity combination in weeks 13+14 and assign to next nearest class.
d.S4[which(d.S4$Site=="Showers"&as.character(d.S4$BurnSev)=="4" & 
             as.character(d.S4$C.Pos)=="U"),"C.Pos"] = factor("E")
d.S4[which(d.S4$Site=="Showers"&as.character(d.S4$BurnSev)=="1" & 
             as.character(d.S4$C.Pos)=="U"),"BurnSev"] = factor("2")
d.S4[which(d.S4$Site=="Showers"&as.character(d.S4$BurnSev)=="1" & 
             as.character(d.S4$C.Pos)=="O"),"BurnSev"] = factor("2")
d.FigS2 <- ddply(d.S4,.(Week,BurnSev,C.Pos),summarize,
                 MeanSnowDepth.Angora=round(mean(Depth.Angora,na.rm=T),2),
                 sdSnowDepth.Angora=round(sd(Depth.Angora,na.rm=T),2),
                 MeanSnowDepth.Reading=round(mean(Depth.Reading,na.rm=T),2),
                 sdSnowDepth.Reading=round(sd(Depth.Reading,na.rm=T),2),
                 MeanSnowDepth.Showers=round(mean(Depth.Showers,na.rm=T),2),
                 sdSnowDepth.Showers=round(sd(Depth.Showers,na.rm=T),2))


FS2.melt <- melt(d.FigS2[-which(is.na(d.FigS2$Week)),c(1,2,3,4,6,8)],
                 id=c("Week", "BurnSev", "C.Pos"))
FS2.melt <- FS2.melt[-which(is.na(FS2.melt$value)),] #Remove NA's
levels(FS2.melt$variable) = c("Angora","Reading","Showers") #Shorten labels
ggplot(FS2.melt)+
  geom_line(aes(x=as.integer(as.character(Week)),y=value,lty=C.Pos,col=BurnSev))+
  facet_wrap(~variable,ncol=1)+
  scale_color_manual(values=c("#2b83ba","#66c2a5","#ffc425","#e0301e","#740001"))+
  scale_linetype_manual(values=c("solid","longdash","dotted"))+
  labs(y= "mean snow depth (cm)", x= "Week",col="Burn \nSeverity \nClass",lty="C.Pos")+
  theme_bw()
#dev.copy2pdf(file="Figures/FigS2.pdf",width=7,height=7)

#2e: Investigating SWE (Table S1) for Reading fire in week 15
TableS1 <- d.S4[d.S4$Site=="Reading" & !is.na(match(d.S4$Week,"15")),]
TableS1 <- TableS1[which(!is.na(TableS1$SWE)),]
TableS1 <- 
  TableS1[,colnames(TableS1) %in% 
            c("pt.ID", "Depth", "SWE", "SWE2", "Week", "C.Pos", "BurnSev")]
TableS1$SWE2 <- as.character(TableS1$SWE2)
#write.cb(TableS1) #Copy for Table S1
sd(TableS1$SWE) #34.9, standard deviation in mean density
SWE_sd <- function(SWE2){
  return(sd(scan(text=SWE2, sep="/", what = numeric(), quiet=TRUE)))
}
#mean of standard deviation in within-point density at different depths:
mean(sapply(TableS1$SWE2,SWE_sd),na.rm = TRUE) #18.03
#write.cb(TableS1) #Copy for Table S1

####3. Plots of Models####
#3a: Burn Severity
m.pred.2a <- #Model of burn severity holding other variables constant
  lmer(Depth~ 0 + BurnSev + (1|Week/Site) + (1|C.Pos) +(1|AspectCat),data=d) 
coefs.2a <- data.frame(coef(summary(m.pred.2a)))
coefs.2a$BurnSev <- c("Class 0","Class 1","Class 2", "Class 3", "Class 4")
coefs.2a$Std..Error <- with(coefs.2a,ifelse(Estimate-Std..Error<0,Estimate,Std..Error)) #Truncate SE at 0.
p.2a <- ggplot(coefs.2a,aes(x=BurnSev,y=Estimate))+
  geom_bar(stat="identity",fill="gray") +
  geom_errorbar(aes(ymax=Estimate+Std..Error,ymin=Estimate-Std..Error,width=0.2)) +
  theme_bw() +
  labs(x="Burn severity class",y="") +
  theme(plot.margin=unit(c(0.2,0.2,0.3,-0.5),"cm"))
#p.2a

#3b: Canopy position
m.pred.2b <- #Model of canopy position holding other variables constant
  lmer(Depth~ 0 + C.Pos + (1|Week/Site) + (1|BurnSev) +(1|AspectCat),data=d) 
coefs.2b <- data.frame(coef(summary(m.pred.2b))); 
coefs.2b$C.Pos=factor(c("Under","Edge","Open"),levels=c("Under","Edge","Open"))
p.2b=ggplot(coefs.2b,aes(x=C.Pos,y=Estimate))+
  geom_bar(stat="identity",fill="gray")+
  geom_errorbar(aes(ymax=Estimate+Std..Error,ymin=Estimate-Std..Error),width=0.2)+
  theme_bw()+
  labs(x="Overhead canopy condition",y="") +
  theme(plot.margin=unit(c(0.3,0.2,0.2,-0.5),"cm"))
#p.2b

ylabel=textGrob("Model estimate of snow depth (cm)",rot=90,vjust=0.5)
grid.arrange(ylabel,arrangeGrob(p.2a,p.2b,
                                top = textGrob("A", x=0, hjust = -0.1),
                                left = textGrob("B", x=0, y=0.48, hjust = -1.7),
                                padding = unit(0.1, "line"),
                                nrow=2),nrow=1,widths=c(0.05,0.95) )

#dev.copy2pdf(file="./Figures/Fig2.pdf",width=3.75,height=3.75)

#3c Time since snowfall effects 
m.full <- #Full model with all individual terms and no interactions
  lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Site/Week),data=d)
df.3d <- data.frame(RE= #Random effects for Site/Week, i.e. mean snow depth for Site/Week
                      ranef(m.full)$`Week:Site`[,1],  
  Site=c(rep("Angora",3),"Showers","Reading",rep("Angora",2),"Showers","Angora","Showers","Reading"),
  TSS=c(15,4,9,9,15,1,1,2,6,8,15),
  StormTot=c(40.5,36.6,7.4,18.8,20.4,15.7,36.9,46.3,36.9,46.3,42.6))
ggplot(df.3d,aes(x=TSS,y=RE,col=Site))+
  geom_point()+
  geom_smooth(aes(col=Site))
summary(lm(RE~TSS,data=df.3d)) #Change predictor variable, either TSS (time since snow) or StormTot.
#In conclusion, neither the time since fire or the storm total amount were significant predictors of mean snow depth at a given site. Week random effects account for overall site effects (they are nested in 'site', thus are positive or negative adjustments to the overall site effect)

####4. Reading Patch Area analysis####
#Depth~Patch Size, by severity
pd <- d[d$Site=="Reading" & d$BurnSev!=0 & !is.na(d$Depth),] #Edit min() or max()
#Descriptive stats on pd:
#ddply(pd[pd$Week==max(pd$Week,na.rm=T),],.(BurnSev),summarize,count=length(unique(PatchAreaHA_Regular)), mean=mean(PatchAreaHA_Regular)) #Only run on second week because all points were sampled then.
m.pa0 <- lmer(Depth~1+BurnSev+(1|Week),data=pd)
m.pa1 <- lmer(Depth~1+PatchAreaHA_Regular+BurnSev+(1|Week),data=pd)
m.pa2 <- lmer(Depth~PatchAreaHA_Regular*BurnSev+(1|Week),data=pd)
GetME_PVals(m.pa2)
BICtab(m.pa0,m.pa1,m.pa2)

w1 <- ggplot(pd[pd$Week==11,], 
             aes(y=Depth,x=PatchAreaHA_Regular,col=BurnSev, shape=BurnSev))+
  geom_point(size=1.8)+stat_smooth(method="lm",alpha=0.2)+
  scale_color_manual(values=c("#66c2a5","#ffc425","#e0301e","#740001"),guide=F)+
  scale_shape_manual(values=c(16,17,18,3),guide=F)+
  #scale_linetype_manual(values=c(1,2,3,4),guide=F)+
  theme_bw()+
  labs(y= "snow depth (cm)",title="Week 11")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.1,0.2,0.2,-0.5),"cm"))

w2 <- ggplot(pd[pd$Week==15,], 
          aes(y=Depth,x=PatchAreaHA_Regular,col=BurnSev,shape=BurnSev))+
  geom_point()+stat_smooth(method="lm",alpha=0.2)+
  scale_color_manual(values=c("#66c2a5","#ffc425","#e0301e","#740001"))+
  scale_shape_manual(values=c(16,17,18,3))+
  #scale_linetype_manual(values=c(1,2,3,4))+
  theme_bw()+
  labs(title="Week 15",col="Burn \nSeverity \nClass",
       shape="Burn \nSeverity \nClass",lty="Burn \nSeverity \nClass") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(0.1,0.2,0.2,0.2),"cm"))

xlabel=textGrob("patch area (ha)")
grid.arrange(arrangeGrob(w1,w2,
                         top = textGrob("B", x=0.45, vjust = 2.5),
                         left = textGrob("A", y=0.9, hjust = -1),
                         padding = unit(0.5, "line"),nrow=1, widths = c(0.435,0.565)),
             xlabel,
             nrow=2,heights=c(0.95,0.05) )

#dev.copy2pdf(file="./Figures/Fig3.pdf",width=7.5,height=4)

####5. Variance in snow depth at different severity levels####
#Variance~Severity
d.var=ddply(d[!is.na(d$Week) & d$Site=="Reading",],
  .(BurnSev,Week),summarize,Var=var(Depth,na.rm=T))
m.var=lm(Var~BurnSev, data=d.var)
ggplot(d.var,aes(x=factor(BurnSev),y=Var,pch=factor(Week)))+geom_point()+
  stat_smooth(method="lm",col="black",fill="gray",aes(group=1))+
  theme_bw()+
  labs(y= "variance in snow depth \nwithin severity class", 
       x= "burn severity class", pch="Week")
#dev.copy2pdf(file="./Figures/Fig4.pdf",width=3.75,height=3.5)
#For both weeks, the variance is highest at low severity

####6. Estimate SWE and water volume for Reading Week 15####
swe.d <- d[d$Week%in%c(15),c("Site","pt.ID","Depth","C.Pos","BurnSev")]
swe.d$SWE_m <- ((swe.d$Depth/100) * 439)/1000 #Snow water equivalent at the sample point, in m. Convert depth in cm to depth in m, multiply by density of snow (439 kg/m3), and divide by density of water (1000 kg/m3)
swe.d$Volume_m3 <- swe.d$SWE_m * 3600 #Volume of water contained in snowpack represented by sample point, in cubic m. Each sample point represents 0.36 ha, or 3600 m2.

agg.d <- ddply(swe.d,.(BurnSev,C.Pos), summarize,
               TotalVolume=sum(Volume_m3),
               TotalSWE=round(sum(SWE_m),3),
               AreaCoveredM2=round(length(Depth)*3600,3),
               AreaCoveredHA=round(length(Depth)*0.36,3),
               FractionCovered=round(length(Depth)/520,3)
               )
agg.d$ComboClass <- paste(agg.d$BurnSev,agg.d$C.Pos,sep="_")
agg.d <- mutate(agg.d,AreaHACum = cumsum(AreaCoveredHA))
agg.d$FractionVolume <- round(agg.d$TotalVolume/sum(agg.d$TotalVolume),2)
agg.d$TotalVolumeCalc <- agg.d$TotalSWE*agg.d$AreaCoveredM2

#Set up positioning for bars in plot:
agg.d$xpos[1] <- agg.d$AreaHACum[1]/2
for(i in 2:nrow(agg.d)){
  agg.d$xpos[i] <- agg.d$AreaHACum[i-1]+(agg.d$AreaCoveredHA[i]/2)
}
agg.d$C.Pos.size <- c(rep(3,9),1.5,1.5,3,3)
agg.d$FractionVolume.size <- c(rep(3,6),2,2,3,0,0,3,3)
agg.d$FractionVolumeAngle <- c(0,89,0,89,89,0,89,89,0,90,90,0,0)
agg.d$FractionVolume[agg.d$FractionVolume==0]=""

#Draw plot
ggplot(agg.d)+
  geom_bar(aes(y=TotalVolume,x=xpos,fill=BurnSev),
           width=agg.d$AreaCoveredHA,
           stat="identity",position="identity",col="black")+
  scale_fill_manual(values=c("#2b83ba","#66c2a5","#ffc425","#e0301e","#740001"))+
  geom_point(aes(x=xpos,y=TotalVolume+800,pch=C.Pos),size=agg.d$C.Pos.size)+
  scale_shape_manual(values=c("U","E","O"),guide=F)+
  geom_text(aes(x=xpos,y=TotalVolume/2,label=FractionVolume),
            size=agg.d$FractionVolume.size,fontface="bold",angle=agg.d$FractionVolumeAngle)+
  labs(y= expression(atop("total water", paste("volume in snowpack (",m^3,")"))), 
       x= "cumulative landscape area (ha)", fill= "Burn \nSeverity \nClass")+
  theme_bw()+
  theme(legend.position=c(0.8,0.6))

#dev.copy2pdf(file="./Figures/Fig5.pdf",width=7.5,height=4)


####7. Calculate Canopy Cover####
tmp=ddply(d[d$Site=="Reading" & !is.na(d$Depth) & d$Week==15,],
          .(Week, BurnSev, C.Pos),summarize, L=length(na.exclude(Depth)))
ddply(tmp,.(BurnSev),summarize,L2=1-L[3]/sum(L)) 
#43% of points in the unburned were in the open.
#Canopy cover estimates are 57, 42, 25 and 9, respectively.