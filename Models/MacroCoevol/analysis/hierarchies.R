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
    s0 = as.tbl(reshaped[reshaped$time==0,]) %>% group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>% summarize(value=mean(value))
    sf = as.tbl(reshaped[reshaped$time==finalTime,]) %>% group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>% summarize(value=mean(value))
    if(ncol(timedata)==0){
      timedata = reshaped
      colnames(timedata)[ncol(timedata)]<-finalindicname
      sdata = s0; sdata$value = sf$value - s0$value; colnames(sdata)[ncol(sdata)]<-finalindicname
    }else{
      timedata=cbind(timedata,reshaped[,"value"])
      colnames(timedata)[ncol(timedata)]<-finalindicname
      sdata = cbind(sdata,value=sf$value - s0$value); colnames(sdata)[ncol(sdata)]<-finalindicname
    }
  }
}

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

#sres = data.frame()
#for(indic in allindics){
#  currentdata = timedata[,c(params,indic)];names(currentdata)[ncol(currentdata)]<-"indic"
#  scurrentdata = as.tbl(currentdata) %>% group_by(synthRankSize,nwExponent,nwThreshold,gravityWeight,gravityGamma,gravityDecay) %>%
#    summarize(sdindic=sd(indic[time==0]),indic=mean(indic[time==0]))
#  if(ncol(sres)==0){sres=scurrentdata}else{sres=cbind(sres,scurrentdata[,c("indic","sdindic")])}
#  colnames(sres)[(ncol(sres)-2):ncol(sres)]=c(indic,paste0(sd,indic))
#}










