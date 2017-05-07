#Snowpack Analysis: Mixed model
#Jens Stevens 6/25/15

library(ggplot2);library(lme4);library(plyr);library(gridExtra);library(bbmle);
#library(ez);library(lattice);library(cowplot)
library(Kmisc);library(dplyr)
library(grid)
#setwd("/Users/Jens/Documents/Davis/Research/6.Snowpack Fuels Project/FTE Snowpack/Analysis/CSVs for R")
GetME_PVals=function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs #Everything is significantly different from UB
  return(coefs)
}#Get approximate p-vals (http://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html)


####1. Setup data####
Sites=c("Reading","Angora","Showers")
Angora=read.csv("CSVs for R/Angora.csv");Reading=read.csv("CSVs for R/Reading.csv");Showers=read.csv("CSVs for R/Showers.csv")
Angora$pt.ID=as.character((Angora$pt.ID));Reading$pt.ID=as.character((Reading$pt.ID));Showers$pt.ID=as.character((Showers$pt.ID))
#Temp: C.Pos2 might be the right one to use for Reading, C.Pos1 seems to have possibly counted "under dead tree" as "under" when it should have been "open". I think this is right, but play around with it some more.
Reading$C.Pos=Reading$C.Pos2
d=rbind(Angora,Reading,Showers) 
d$Week=factor(d$Week);d$BurnSev=factor(d$BurnSev)
d$C.Pos=factor(d$C.Pos,levels=c("","U","E","O"))

####2. Mixed model####
#2a: Gaussian model (although snow depth is not normally distributed).

#d.full=d[-which(is.na(d$Depth)),] #Optional if subsetting; ONLY DO ONCE
#for(s in 1:1000){ #Optional if subsetting
#d=d.full[-sample(1:nrow(d.full), nrow(d.full)/2, replace=FALSE),] #Subset to half the data; Optional if subsetting
m0=lmer(Depth~ (1|Week/Site),data=d)
m1a=lmer(Depth~ AspectCat + (1|Week/Site),data=d) #Aspect effect: SW slopes have less
m1b=lmer(Depth~ BurnSev + (1|Week/Site),data=d) #Burn effect: Significant negative effect of severity on depth.
m1c=lmer(Depth~ C.Pos + (1|Week/Site),data=d) #Overhead canopy effect: Significant negative effect of canopy on depth
m2a=lmer(Depth~ BurnSev + C.Pos + (1|Week/Site),data=d) #Burn effect controlling for C.Pos: Stronger than before!
m2b=lmer(Depth~ BurnSev + (1|Week/Site) + (1|C.Pos) ,data=d) #Burn effect controlling for C.Pos as a random effect: Stronger than before!
m2c=lmer(Depth~ C.Pos + (1|Week/Site) + (1|BurnSev),data=d) #Overhead canopy effect controlling for burn: Stronger!
m2d=lmer(Depth~ BurnSev + Aspect + (1|Week/Site),data=d) #Burn effect controlling for Aspect
m3b=lmer(Depth~ BurnSev + (1|Week/Site) + (1|C.Pos) +(1|AspectCat),data=d) 
m3c=lmer(Depth~ C.Pos + (1|Week/Site) + (1|BurnSev) +(1|AspectCat),data=d) 
m4b=lmer(Depth~ BurnSev*AspectCat + (1|Week/Site),data=d)
m4c=lmer(Depth~ C.Pos*AspectCat + (1|Week/Site),data=d) #Interaction here; Benefit of open canopy is much greater on NE aspects.
m.full=lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Site/Week),data=d)
m.inter=lmer(Depth~ C.Pos * BurnSev + AspectCat + (1|Week/Site),data=d)#Not fully crossed so can't include Aspect interaction.
BICt=BICtab(m0,m1a,m1b,m1c,m2a,m2b,m2c,m2d,m3b,m3c,m.full,m.inter)# m.full best; 

#Calculate R2 to append to table
#Based on 
#1) https://www.r-bloggers.com/mixed-model-r2-updated/
#2) https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/
#3) Nakagawa and Schielzeth 2013
m.list=list(m.full,m3b,m3c,m.inter,m2d,m2a,m2b,m2c,m1b,m1a,m0,m1c)#Arrange in order of BITc
#Read code from https://github.com/jslefche/rsquared.glmm/blob/master/rsquaredglmm.R
rsquares=rsquared.glmm(m.list)
#write.cb(cbind(attr(BICt,"row.names"),list2df(BICt),round(rsquares$Marginal),rsquares$Conditional))

#Optional: Get the coefficients from m.full for a single fire (Table S4)
d.site=d[d$Site=="Showers",]
m.full=lmer(Depth~ C.Pos + BurnSev + AspectCat + (1|Week),data=d.site)

Table1=GetME_PVals(m.full) #All the trends hold up in the best model, everything is still pretty significant.
Table1=round(Table1,3);#write.cb(Table1) #Error here in write.cb is apparently ok.
#if(s==1){ #Optional if subsetting
#  param.dist=data.frame(Sim=1,Var=rownames(Table1),Estimate=Table1$Estimate)} else {
#    param.dist[((s-1)*nrow(Table1)+1):(s*nrow(Table1)+nrow(Table1)),"Sim"]=s
#    param.dist[((s-1)*nrow(Table1)+1):(s*nrow(Table1)+nrow(Table1)),"Var"]=rownames(Table1)
#    param.dist[((s-1)*nrow(Table1)+1):(s*nrow(Table1)+nrow(Table1)),"Estimate"]=Table1$Estimate
#  }
#print(s) #Optional if subsetting
#Take-away: Negative effect of severity and positive effect of overhead canopy openness on snow depth, controlling for aspect and variation among sites and weeks.
#}#End subsetting loop (optional)

#2b: Make parameter distribution figure (optional, if subsetting)
#levels(param.dist$Var)=c("Intercept","NE Aspect", "SW Aspect", "Severity Class 1"," Severity Class 2", "Severity Class 3", "Severity Class 4", "Canopy Condition: Edge", "Canopy Condition: Open")
#AppendixParams=ggplot(param.dist)+
#  geom_density(aes(x=Estimate))+
#  geom_vline(xintercept=0,linetype="longdash")+
#  facet_wrap(~Var,ncol=3)+
#  theme_bw()

#dev.copy2pdf(file="Figures/FigA2.pdf",width=7,height=7)
  
#2c: Summarize actual mean snow depths across scenarios (Optional, for a)
d.tmp=d
d.tmp$Depth.Angora[d.tmp$Site=="Angora"]=d.tmp[d.tmp$Site=="Angora","Depth"]
d.tmp$Depth.Reading[d.tmp$Site=="Reading"]=d.tmp[d.tmp$Site=="Reading","Depth"]
d.tmp$Depth.Showers[d.tmp$Site=="Showers"]=d.tmp[d.tmp$Site=="Showers","Depth"]
#For clearer visualization, take the 3 points at Showers that were the only representative from their canopy class-burn severity combination in weeks 13+14 and assign to next nearest class.
d.tmp[which(d.tmp$Site=="Showers"&as.character(d.tmp$BurnSev)=="4"&as.character(d.tmp$C.Pos)=="U"),"C.Pos"]=factor("E")
d.tmp[which(d.tmp$Site=="Showers"&as.character(d.tmp$BurnSev)=="1"&as.character(d.tmp$C.Pos)=="U"),"BurnSev"]=factor("2")
d.tmp[which(d.tmp$Site=="Showers"&as.character(d.tmp$BurnSev)=="1"&as.character(d.tmp$C.Pos)=="O"),"BurnSev"]=factor("2")


AppendixA2=ddply(d.tmp,.(Week,BurnSev,C.Pos),summarize,
                 MeanSnowDepth.Angora=round(mean(Depth.Angora,na.rm=T),2),
                 sdSnowDepth.Angora=round(sd(Depth.Angora,na.rm=T),2),
                 MeanSnowDepth.Reading=round(mean(Depth.Reading,na.rm=T),2),
                 sdSnowDepth.Reading=round(sd(Depth.Reading,na.rm=T),2),
                 MeanSnowDepth.Showers=round(mean(Depth.Showers,na.rm=T),2),
                 sdSnowDepth.Showers=round(sd(Depth.Showers,na.rm=T),2))
#write.cb(AppendixA2)

#Multi-panel figure
library(reshape2)
A2.melt=melt(AppendixA2[-which(is.na(AppendixA2$Week)),c(1,2,3,4,6,8)],id=c("Week", "BurnSev", "C.Pos"))
A2.melt=A2.melt[-which(is.na(A2.melt$value)),] #Remove NA's
levels(A2.melt$variable)=c("Angora","Reading","Showers") #Shorten labels
ggplot(A2.melt)+
  geom_line(aes(x=as.integer(as.character(Week)),y=value,lty=C.Pos,col=BurnSev))+
  facet_wrap(~variable,ncol=1)+
  scale_color_manual(values=c("#2b83ba","#66c2a5","#ffc425","#e0301e","#740001"))+
  scale_linetype_manual(values=c("solid","longdash","dotted"))+
  labs(y= "mean snow depth (cm)", x= "Week",col="Burn \nSeverity \nClass",lty="C.Pos")+
  theme_bw()
dev.copy2pdf(file="Figures/FigA3.pdf",width=7,height=7)

#2d: Investigating SWE
d.Reading.15=d.tmp[d.tmp$Site=="Reading" & !is.na(match(d.tmp$Week,"15")),]
d.Reading.15.SWE=d.Reading.15[which(!is.na(d.Reading.15$SWE)),]
d.Reading.15.SWE=d.Reading.15.SWE[,colnames(d.Reading.15.SWE)%in%c("pt.ID", "Depth", "SWE", "SWE2", "Week", "C.Pos", "BurnSev")]

AppendixA3=d.Reading.15.SWE
write.cb(AppendixA3)
sd(AppendixA3$Depth)

####3. Plots of Models####
#3a: Burn Severity (thresholds are generally 5%, 25%, 75%, from MTBS report)
m.pred.3b=lmer(Depth~ 0+ BurnSev + (1|Week/Site) + (1|C.Pos) +(1|AspectCat),data=d) 
coefs.3b <- data.frame(coef(summary(m.pred.3b)))
coefs.3b$BurnSev=c("Class 0","Class 1","Class 2", "Class 3", "Class 4")
coefs.3b$Std..Error=with(coefs.3b,ifelse(Estimate-Std..Error<0,Estimate,Std..Error))
p.3b=ggplot(coefs.3b,aes(x=BurnSev,y=Estimate))+geom_bar(stat="identity",fill="gray")+
  geom_errorbar(aes(ymax=Estimate+Std..Error,ymin=Estimate-Std..Error,width=0.2))+
  theme_bw()+labs(x="Burn severity class",y="")
#p.3b

#3b: Canopy position
m.pred.3c=lmer(Depth~ 0+ C.Pos + (1|Week/Site) + (1|BurnSev) +(1|AspectCat),data=d) 
coefs.3c <- data.frame(coef(summary(m.pred.3c))); 
coefs.3c$C.Pos=factor(c("Under","Edge","Open"),levels=c("Under","Edge","Open"))
p.3c=ggplot(coefs.3c,aes(x=C.Pos,y=Estimate))+geom_bar(stat="identity",fill="gray")+
  geom_errorbar(aes(ymax=Estimate+Std..Error,ymin=Estimate-Std..Error),width=0.2)+
  theme_bw()+labs(x="Overhead canopy condition",y="")
#p.3c

ylabel=textGrob("Model estimate of snow depth (cm)",rot=90,vjust=0.5)
grid.arrange(ylabel,arrangeGrob(p.3b,p.3c,nrow=2),nrow=1,widths=c(0.05,0.95))
#plot_grid(wv1, wv2, nrow=2,labels = c("A", "B")) Cowplot might not work here because of the wierd way I did this.

#grid.arrange(p.3b,p.3c,ncol=1,widths=c(3,3),left="Model estimate of snow depth (cm)",clip=F) #Deprecated?
#dev.copy2pdf(file="/Users/Jens/Documents/Davis/Research/6.Snowpack Fuels Project/FTE Snowpack/Manuscript/Figures/Fig2.pdf",width=3.75,height=3.75)

#3d Time Since Snow effects
df.3d=data.frame(RE=ranef(m.full)$`Week:Site`[,1],
  Site=c(rep("Angora",3),"Showers","Reading",rep("Angora",2),"Showers","Angora","Showers","Reading"),
  TSS=c(15,4,9,9,15,1,1,2,6,8,15),
  StormTot=c(40.5,36.6,7.4,18.8,20.4,15.7,36.9,46.3,36.9,46.3,42.6))
ggplot(df.3d,aes(x=TSS,y=RE,col=Site))+
  geom_point()+
  geom_smooth(aes(col=Site))
summary(lm(RE~TSS,data=df.3d)) #Change predictor variable
#In conclusion, neither the time since fire or the storm total amount were significant predictors of mean snow depth at a given site. Week random effects account for overall site effects (they are nested in 'site', thus are positive or negative adjustments to the overall site effect)

####4. Reading Patch Area analysis####
#Depth~Patch Size, by severity
pd=Reading[Reading$BurnSev!=0,]#Edit min() or max()
#Descriptive stats on pd:
#ddply(pd[pd$Week==max(pd$Week,na.rm=T),],.(BurnSev),summarize,count=length(unique(PatchAreaHA_Regular)), mean=mean(PatchAreaHA_Regular)) #Only run on second week because all points were sampled then.
m.pa0=lmer(Depth~1+as.factor(BurnSev)+(1|Week),data=pd)
m.pa1=lmer(Depth~1+PatchAreaHA_Regular+as.factor(BurnSev)+(1|Week),data=pd)
m.pa2=lmer(Depth~PatchAreaHA_Regular*as.factor(BurnSev)+(1|Week),data=pd)
GetME_PVals(m.pa2)
BICtab(m.pa0,m.pa1,m.pa2)
w1=ggplot(pd[pd$Week==11,], aes(y=Depth,x=PatchAreaHA_Regular,col=factor(BurnSev),shape=factor(BurnSev)))+
  geom_point(size=1.8)+stat_smooth(method="lm",alpha=0.2)+
  scale_color_manual(values=c("#66c2a5","#ffc425","#e0301e","#740001"),guide=F)+
  scale_shape_manual(values=c(16,17,18,3),guide=F)+
  #scale_linetype_manual(values=c(1,2,3,4),guide=F)+
  theme_bw()+labs(y= "snow depth (cm)", x= "patch area (ha)",title="Week 11")
w2=ggplot(pd[pd$Week==15,], aes(y=Depth,x=PatchAreaHA_Regular,col=factor(BurnSev),shape=factor(BurnSev)))+
  geom_point()+stat_smooth(method="lm",alpha=0.2)+
  scale_color_manual(values=c("#66c2a5","#ffc425","#e0301e","#740001"))+
  scale_shape_manual(values=c(16,17,18,3))+
  #scale_linetype_manual(values=c(1,2,3,4))+
  #scale_color_manual(values=c("#2b83ba","#66c2a5","#fdae61","#d7191c"))+ #original colors
  theme_bw()+labs(y= "snow depth (cm)", x= "patch area (ha)",title="Week 15",col="Burn \nSeverity \nClass",shape="Burn \nSeverity \nClass",lty="Burn \nSeverity \nClass")
f3=grid.arrange(w1,w2,nrow=1,widths=c(0.5,0.6))
#dev.copy2pdf(file=paste0("Figures/Fig3",Sys.Date(),".pdf"),width=7.5,height=4)

####5. Variance in snow depth at different severity levels####
#Variance~Severity
d.var=ddply(Reading[-which(is.na(Reading$Week)),],
  .(BurnSev,Week),summarize,Var=var(Depth,na.rm=T))
m.var=lm(Var~BurnSev, data=d.var)
f4=ggplot(d.var,aes(x=factor(BurnSev),y=Var,pch=factor(Week)))+geom_point()+
  stat_smooth(method="lm",col="black",fill="gray",aes(group=1))+theme_bw()+
  labs(y= "variance in snow depth \nwithin severity class", x= "burn severity class", pch="Week")
f4
#dev.copy2pdf(file="/Users/Jens/Documents/Davis/Research/6.Snowpack Fuels Project/FTE Snowpack/Manuscript/Figures/Fig4.pdf",width=3.75,height=3.5)

#For both weeks, the variance is highest at low severity

####6. Plots of each variable combination of snow depth for the appendix####
#Deprecated now?
df.base=expand.grid(BurnSev=as.factor(c(0:4)),C.Pos=as.factor(c("O","U","E")))
l=list();p=list()
l[["ac"]]=factor(rep(c("FLAT","SW","NE"),3))
l[["s"]]=factor(rep(c("Angora","Reading","Showers"),each=3))
l[["w"]]=factor(rep(c(13,15,13),each=3))
l[["label"]]=paste(l[["ac"]],l[["s"]],l[["w"]])
for(i in 1:length(l[["ac"]])){
  df.base$AspectCat=l[["ac"]][i];df.base$Site=l[["s"]][i];df.base$Week=l[["w"]][i]
  pred.d=ezPredict(m.full,df.base)
  pd=cbind(df.base,Depth=pred.d$cells$value)
  p[[i]]=wireframe(Depth ~ BurnSev*C.Pos, data = pd,
          xlab = "Burn Severity Class", ylab = "Overhead Canopy Class",
          main = l[["label"]][i],
          scales=list(arrows=FALSE),
          #drape = TRUE,
          #colorkey = TRUE,
          screen = list(z = -60, x = -60)
  )
}
do.call("grid.arrange",c(p,ncol=3))

####7a. Estimate SWE and water volume for Reading Week 15####
swe.d=d[d$Week%in%c(15),c("Site","pt.ID","Depth","C.Pos","BurnSev")]
swe.d$SWE_m=((swe.d$Depth/100) * 439)/1000 #Snow water equivalent at the sample point, in m. Convert depth in cm to depth in m, multiply by density of snow (439 kg/m3), and divide by density of water (1000 kg/m3)
swe.d$Volume_m3=swe.d$SWE_m * 3600 #Volume of water contained in snowpack represented by sample point, in cubic m. Each sample point represents 0.36 ha, or 3600 m2.

agg.d=ddply(swe.d,.(BurnSev,C.Pos), summarize,
            TotalVolume=sum(Volume_m3),
            TotalSWE=round(sum(SWE_m),3),
            AreaCoveredM2=round(length(Depth)*3600,3),
            AreaCoveredHA=round(length(Depth)*0.36,3),
            FractionCovered=round(length(Depth)/520,3)
            )
agg.d$ComboClass=paste(agg.d$BurnSev,agg.d$C.Pos,sep="_")
agg.d=agg.d %>% mutate(AreaHACum = cumsum(AreaCoveredHA))
agg.d$FractionVolume=round(agg.d$TotalVolume/sum(agg.d$TotalVolume),2)
agg.d$TotalVolumeCalc=agg.d$TotalSWE*agg.d$AreaCoveredM2

#Set up positioning for bars in plot:
agg.d$xpos[1]=agg.d$AreaHACum[1]/2
for(i in 2:nrow(agg.d)){
  agg.d$xpos[i]=agg.d$AreaHACum[i-1]+(agg.d$AreaCoveredHA[i]/2)
}
agg.d$C.Pos.size=c(rep(3,9),1.5,1.5,3,3)
agg.d$FractionVolumeAngle=c(0,89,0,89,89,0,89,89,0,90,90,0,0)
#Draw plot
p5=ggplot(agg.d)+
  geom_bar(aes(y=TotalVolume,x=xpos,width=AreaCoveredHA,fill=BurnSev),stat="identity",position="identity",col="black")+
  scale_fill_manual(values=c("#2b83ba","#66c2a5","#ffc425","#e0301e","#740001"))+
  geom_point(aes(x=xpos,y=TotalVolume+300,pch=C.Pos),size=agg.d$C.Pos.size)+
  scale_shape_manual(values=c("U","E","O"),guide=F)+
  geom_text(aes(x=xpos,y=TotalVolume/2,label=FractionVolume),size=3,fontface="bold",angle=agg.d$FractionVolumeAngle)+
  labs(y= expression(atop("total water", paste("volume in snowpack (",m^3,")"))), x= "cumulative landscape area (ha)", fill= "Burn \nSeverity \nClass")+
  theme_bw()+
  theme(legend.position=c(0.8,0.6))
#p5
dev.copy2pdf(file=paste0("Figures/Fig5",Sys.Date(),".pdf"),width=7.5,height=5)


####8. Calculate Canopy Cover####
tmp=ddply(Reading,.(Week, BurnSev, C.Pos),summarize, L=length(na.exclude(Depth)))
ddply(tmp,.(Week,BurnSev),summarize,L2=L[2]/sum(L))
