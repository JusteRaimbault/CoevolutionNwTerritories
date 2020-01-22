setwd(paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Models/MacroCoevol/MacroCoevol'))

library(dplyr)
library(ggplot2)
library(reshape2)

source(paste0(Sys.getenv('CN_HOME'),'/Models/Utils/R/plots.R'))

resprefix = '20191228_175942_HIERARCHIES_SYNTHETICVIRTUAL_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))

resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir)

params = c("synthRankSize","nwExponent","nwThreshold","gravityWeight","gravityGamma","gravityDecay")
indics = c("hierarchiesPop","hierarchiesCloseness","hierarchiesAccessibility",
           "segHierarchiesPop","segHierarchiesCloseness","segHierarchiesAccessibility",
           "rankCorrPop","rankCorrCloseness","rankCorrAccessibility",
           "rankCorrsPopCloseness","rankCorrsPopAccessibility","rankCorrsClosenessAccessibility"
           )

finalTime = 30

# reshape the data to long format
# specific indic format for multidim indics (alpha + rsquared for example)
timeindicsstruct = list("hierarchiesPop"=c("Alpha","RSquared"),"hierarchiesCloseness"=c("Alpha","RSquared"),"hierarchiesAccessibility"=c("Alpha","RSquared"),
                    "segHierarchiesPop"=c("Alpha1","Alpha2","Psi","RSQuared"),"segHierarchiesCloseness"=c("Alpha1","Alpha2","Psi","RSQuared"),"segHierarchiesAccessibility"=c("Alpha1","Alpha2","Psi","RSQuared"),
                    "rankCorrsPopCloseness"=c(""),"rankCorrsPopAccessibility"=c(""),"rankCorrsClosenessAccessibility"=c("")
                    )

# the arrayOnRows option may not be the right one - kind of a hell to switch to a proper long format
# local function - depends on parameter names but not correct level of code in group_by
sumdata <- function(p){
timedata = data.frame()
sdata = data.frame()
for(timeindic in names(timeindicsstruct)){
  k = length(timeindicsstruct[[timeindic]])
  for(subindicind in 0:(k-1)){
    varyinginds = (0:finalTime)*k+subindicind
    finalindicname=paste0(timeindic,timeindicsstruct[[timeindic]][subindicind+1])
    # reshape is not handy to use
    #reshaped = reshape(res[,c(params,paste0(timeindic,varyinginds))],varying = (length(params)+1):(length(params)+length(varyinginds)),direction="long",v.names = c(paste0(timeindic,timeindicsstruct[[timeindic]][subindicind])))
    toreshape=res[,c(params,paste0(timeindic,varyinginds))]
    names(toreshape)[(length(params)+1):(length(params)+length(varyinginds))]<-0:finalTime
    reshaped = melt(toreshape,id.vars=params,variable.name="time",value.name = "value")
    if (p=='virtual'){
     sdelta=left_join(reshaped[reshaped$time==finalTime,],reshaped[reshaped$time==0,],by=c('synthRankSize'='synthRankSize','nwExponent'='nwExponent','nwThreshold'='nwThreshold','gravityWeight'='gravityWeight','gravityGamma'='gravityGamma','gravityDecay'='gravityDecay')) %>%
       group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>% summarize(valueSd=sd(value.x-value.y),value=mean(value.x-value.y))
     #sdelta = as.tbl(reshaped[reshaped$time==finalTime,]-reshaped[reshaped$time==0,]) %>% group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>% summarize(valueSd=sd(value),value=mean(value))
     #sf = as.tbl(reshaped[reshaped$time==finalTime,]) %>% group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>% summarize(value=mean(value))
    }else{
      sdelta=left_join(reshaped[reshaped$time==finalTime,],reshaped[reshaped$time==0,],by=c('synthRankSize'='synthRankSize','nwExponent'='nwExponent','nwPhysQuantile'='nwPhysQuantile','gravityWeight'='gravityWeight','gravityGamma'='gravityGamma','gravityDecay'='gravityDecay')) %>%
        group_by(synthRankSize,nwExponent,nwPhysQuantile,gravityWeight,gravityGamma,gravityDecay) %>% summarize(valueSd=sd(value.x-value.y),value=mean(value.x-value.y))
      #s0 = as.tbl(reshaped[reshaped$time==0,]) %>% group_by(synthRankSize,nwExponent,nwPhysQuantile,gravityWeight,gravityGamma,gravityDecay) %>% summarize(value=mean(value))
      #sf = as.tbl(reshaped[reshaped$time==finalTime,]) %>% group_by(synthRankSize,nwExponent,nwPhysQuantile,gravityWeight,gravityGamma,gravityDecay) %>% summarize(value=mean(value))
    }
    if(ncol(timedata)==0){
      timedata = reshaped
      colnames(timedata)[ncol(timedata)]<-finalindicname
      sdata = sdelta;colnames(sdata)[(ncol(sdata)-1):ncol(sdata)]<-c(paste0(finalindicname,'Sd'),finalindicname)
    }else{
      timedata=cbind(timedata,reshaped[,"value"])
      colnames(timedata)[ncol(timedata)]<-finalindicname
      sdata = cbind(sdata,valueSd=sdelta$valueSd,value=sdelta$value)
      colnames(sdata)[(ncol(sdata)-1):ncol(sdata)]<-c(paste0(finalindicname,'Sd'),finalindicname)
    }
  }
}
return(list(timedata=timedata,sdata=sdata))
}
d = sumdata(p="virtual")
timedata = d$timedata;sdata = d$sdata

allindics = names(timedata)[-(1:(length(params)+1))]

###
# indicator plots - in time

filter<-function(d,values){
  rows = rep(T,nrow(d))
  for(indic in names(values)){rows=rows&d[,indic]==values[indic]}
  return(d[rows,])
}

fixedparams = list(
#nwExponent=1,
  synthRankSize=1,
  nwThreshold=2.5,
  #gravityGamma=1,
  gravityWeight=1e-03 #5e-04
)
indic="rankCorrsPopCloseness"

g = ggplot(filter(timedata,fixedparams),aes_string(x="time",y=indic,color="gravityDecay",group="gravityDecay"))
g+geom_point(pch='.')+geom_smooth()+facet_grid(nwExponent~gravityGamma,scales='free')

####
# nothing really interesting
# -> summarize with Delta indic (tf - t0)
# sdata already averaged on realizations

dir.create(paste0(resdir,'summarizedindics'))

for(nwExp in unique(sdata$nwExponent)){
  for(wg in unique(sdata$gravityWeight)){
for(indic in allindics){
    
fixedparams = list(
  nwExponent=nwExp, #0.5 1.0 1.5,
  #synthRankSize=1,
  #nwThreshold=2.5,
  #gravityGamma=1,
  gravityWeight=wg #1e-03 #5e-04
)
#indic="rankCorrsPopCloseness"

g = ggplot(filter(sdata,fixedparams),aes_string(x="gravityDecay",y=indic,color="gravityGamma",group="gravityGamma"))
g+geom_point()+geom_line()+facet_grid(synthRankSize~nwThreshold,scales='free')+stdtheme
ggsave(file=paste0(resdir,"summarizedindics/",indic,'_nwExp',nwExp,'_wG',wg,'_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png'),width=20,height=18,units='cm')

  }
  }
}


#####
# targeted plots

# rankCorrsPopCloseness_nwExp1_wG0.001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png
# - rank corr pop access: how the territory and network self-reinforce
# - qual inversion as a function of initial hierarchy
# in the same setting hierarchiesClosenessAlpha_nwExp1_wG0.001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png
# hierarchy of closenesses is interesting also
# and psi population and/or closenesses

# => rerun new experiments with more repets for standard devs and more gravity_gamma
# 

resprefix = '20200120_175200_HIERARCHIES_SYNTHETICVIRTUAL_TARGETED_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))
resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir);dir.create(paste0(resdir,'targeted'))
params = c("synthRankSize","nwExponent","nwThreshold","gravityWeight","gravityGamma","gravityDecay")
d = sumdata(p='virtual')
timedata = d$timedata;sdata = d$sdata

# rank correlation pop closeness
g = ggplot(sdata,aes(x=gravityDecay,y=rankCorrsPopCloseness,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+
  # error bars are ~ significant but unease reading
  #geom_errorbar(aes(ymin=rankCorrsPopCloseness-rankCorrsPopClosenessSd,ymax=rankCorrsPopCloseness+rankCorrsPopClosenessSd))+
  facet_grid(synthRankSize~nwThreshold,scales='free')+
  xlab(expression(d[G]))+ylab(expression(rho[r*","*Delta]*"[P,C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png"),width=20,height=18,units='cm')

# hierarchy closenesses
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesClosenessAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThreshold,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png"),width=20,height=18,units='cm')

# hierarchy pops
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesPopAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThreshold,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png"),width=20,height=18,units='cm')

# Psi pops
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesPopPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThreshold,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesPopPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png"),width=20,height=18,units='cm')

# Psi closeness
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesClosenessPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThreshold,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png"),width=20,height=18,units='cm')






######
# similar plots with physical network

resprefix = '20200115_184543_HIERARCHIES_SYNTHETICPHYSICAL_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))
resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir)
# one param is different (but ~ equivalent) for the physical network
params = c("synthRankSize","nwExponent","nwPhysQuantile","gravityWeight","gravityGamma","gravityDecay")
# indics are the same

d = sumdata(p='physical')
timedata = d$timedata;sdata = d$sdata

dir.create(paste0(resdir,'summarizedindics'))
for(nwExp in unique(sdata$nwExponent)){
  for(wg in unique(sdata$gravityWeight)){
    for(indic in allindics){

      fixedparams = list(
        nwExponent=nwExp, #0.5 1.0 1.5,
        #synthRankSize=1,
        #nwThreshold=2.5,
        #gravityGamma=1,
        gravityWeight=wg #1e-03 #5e-04
      )
      #indic="rankCorrsPopCloseness"
      
      g = ggplot(filter(sdata,fixedparams),aes_string(x="gravityDecay",y=indic,color="gravityGamma",group="gravityGamma"))
      g+geom_point()+geom_line()+facet_grid(synthRankSize~nwPhysQuantile,scales='free')+stdtheme
      ggsave(file=paste0(resdir,"summarizedindics/",indic,'_nwExp',nwExp,'_wG',wg,'_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThreshold.png'),width=20,height=18,units='cm')
    }
  }
}


####
# similar targeted plots ?


#resprefix = '20200120_194245_HIERARCHIES_SYNTHETICPHYSICAL_TARGETED_GRID' # one param value missing in this grid
resprefix = '20200121_114422_HIERARCHIES_SYNTHETICPHYSICAL_TARGETED_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))
resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir);dir.create(paste0(resdir,'targeted'))
params = c("synthRankSize","nwExponent","nwPhysQuantile","gravityWeight","gravityGamma","gravityDecay")
d = sumdata(p='physical')
timedata = d$timedata;sdata = d$sdata

# rank correlation pop closeness
g = ggplot(sdata,aes(x=gravityDecay,y=rankCorrsPopCloseness,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+
  facet_grid(synthRankSize~nwPhysQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(rho[r*","*Delta]*"[P,C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwPhysQuantile.png"),width=20,height=18,units='cm')

# hierarchy closenesses
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesClosenessAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwPhysQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwPhysQuantile.png"),width=20,height=18,units='cm')

# hierarchy pops
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesPopAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwPhysQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwPhysQuantile.png"),width=20,height=18,units='cm')

# Psi pops
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesPopPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwPhysQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesPopPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwPhysQuantile.png"),width=20,height=18,units='cm')

# Psi closeness
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesClosenessPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwPhysQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwPhysQuantile.png"),width=20,height=18,units='cm')










