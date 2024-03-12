setwd(paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/CoevolutionNwTerritories/Models/MacroCoevol/MacroCoevol'))

library(dplyr, warn.conflicts = F)
library(ggplot2)
library(reshape2)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

resprefix = '20191228_175942_HIERARCHIES_SYNTHETICVIRTUAL_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))

resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir)

params = c("synthRankSize","nwExponent","nwThresholdQuantile","gravityWeight","gravityGamma","gravityDecay")
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
sumdata <- function(){
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
    
     sdelta=left_join(reshaped[reshaped$time==finalTime,],reshaped[reshaped$time==0,],by=c('synthRankSize'='synthRankSize','nwExponent'='nwExponent','nwThresholdQuantile'='nwThresholdQuantile','gravityWeight'='gravityWeight','gravityGamma'='gravityGamma','gravityDecay'='gravityDecay')) %>%
       group_by(synthRankSize,nwExponent,nwThresholdQuantile,gravityWeight,gravityGamma,gravityDecay) %>% summarize(valueSd=sd(value.x-value.y),value=mean(value.x-value.y))
    
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
d = sumdata()
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

#resprefix = '20200120_175200_HIERARCHIES_SYNTHETICVIRTUAL_TARGETED_GRID' # rerun with quantile param
resprefix = '20200123_183857_HIERARCHIES_SYNTHETICVIRTUAL_TARGETED_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))
resdir=paste0(Sys.getenv('CS_HOME'),'/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir);dir.create(paste0(resdir,'targeted'))
params = c("synthRankSize","nwExponent","nwThresholdQuantile","gravityWeight","gravityGamma","gravityDecay")
d = sumdata()
timedata = d$timedata;sdata = d$sdata

# rank correlation pop closeness
g = ggplot(sdata,aes(x=gravityDecay,y=rankCorrsPopCloseness,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+
  # error bars are ~ significant but unease reading
  #geom_errorbar(aes(ymin=rankCorrsPopCloseness-rankCorrsPopClosenessSd,ymax=rankCorrsPopCloseness+rankCorrsPopClosenessSd))+
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(rho[r*","*Delta]*"[P,C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')

# hierarchy closenesses
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesClosenessAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')

# hierarchy pops
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesPopAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')

# Psi pops
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesPopPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesPopPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')

# Psi closeness
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesClosenessPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')


##
## distribution of initial hierarchies of closeness
# (for pse indicators: in order to know if alpha_C of [alpha_C - alpha_C0] is better)

#timedata$synthRankSizeF = cut(timedata$synthRankSize,5) # note: no param influences distrib of closeness
# except for the physical if would have changed nw generation params - but initial pop rank size no as only distances, not access!
g=ggplot(timedata[as.numeric(as.character(timedata$time))%%5==0,],aes(x=hierarchiesClosenessAlpha,color=time,group=time))
g+geom_density()+scale_color_discrete(name='Time')+xlab(expression(alpha[C]))+ylab('Density')+stdtheme
ggsave(file=paste0(resdir,"targeted/alphaClosenessDistributions.png"),width=20,height=18,units='cm')




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
#resprefix = '20200121_114422_HIERARCHIES_SYNTHETICPHYSICAL_TARGETED_GRID' # not correct implementation
#resprefix = '20200122_202902_HIERARCHIES_SYNTHETICPHYSICAL_TARGETED_GRID' # do similar DOE than virtual
resprefix = '20200123_194611_HIERARCHIES_SYNTHETICPHYSICAL_TARGETED_GRID'
res <- as.tbl(read.csv(paste0('exploration/',resprefix,'.csv'),stringsAsFactors = FALSE))
resdir=paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/CoevolutionNwTerritories/Results/MacroCoevol/',resprefix,'/');dir.create(resdir);dir.create(paste0(resdir,'targeted'))
params = c("synthRankSize","nwExponent","nwThresholdQuantile","gravityWeight","gravityGamma","gravityDecay")
d = sumdata()
timedata = d$timedata;sdata = d$sdata

# rank correlation pop closeness
g = ggplot(sdata,aes(x=gravityDecay,y=rankCorrsPopCloseness,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(rho[r*","*Delta]*"[P,C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/physical_rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')
# pdf for paper
ggsave(file=paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/physical_rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf"),width=20,height=18,units='cm')


# hierarchy closenesses
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesClosenessAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/physical_hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')
# pdf for paper
ggsave(file=paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/physical_hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf"),width=20,height=18,units='cm')


# hierarchy pops
g = ggplot(sdata,aes(x=gravityDecay,y=hierarchiesPopAlpha,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(alpha[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/physical_hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')
# pdf for paper
ggsave(file=paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/physical_hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf"),width=20,height=18,units='cm')


# Psi pops
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesPopPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[P]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/physical_segHierarchiesPopPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')
# pdf for paper - not needed
#ggsave(file=paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/physical_segHierarchiesPopPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf"),width=20,height=18,units='cm')


# Psi closeness
g = ggplot(sdata,aes(x=gravityDecay,y=segHierarchiesClosenessPsi,color=gravityGamma,group=gravityGamma))
g+geom_point()+geom_line()+ 
  facet_grid(synthRankSize~nwThresholdQuantile,scales='free')+
  xlab(expression(d[G]))+ylab(expression(Psi[Delta]*"[C]"))+scale_color_continuous(name=expression(gamma[G]))+
  stdtheme
ggsave(file=paste0(resdir,"targeted/physical_segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png"),width=20,height=18,units='cm')
# pdf for paper
ggsave(file=paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/physical_segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf"),width=20,height=18,units='cm')


# distrib closeness
g=ggplot(timedata[as.numeric(as.character(timedata$time))%%5==0,],aes(x=hierarchiesClosenessAlpha,color=time,group=time))
g+geom_density()+scale_color_discrete(name='Time')+xlab(expression(alpha[C]))+ylab('Density')+stdtheme
ggsave(file=paste0(resdir,"targeted/alphaClosenessDistributions.png"),width=20,height=18,units='cm')








