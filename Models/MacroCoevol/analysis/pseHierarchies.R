setwd(paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/CoevolutionNwTerritories/Models/MacroCoevol/MacroCoevol'))

library(dplyr)
library(ggplot2)
library(GGally)

source(paste0(Sys.getenv('CS_HOME'),'/Organisation/Models/Utils/R/plots.R'))

resprefix = '20200124_PSE_SYNTHETICPHYSICAL_GRID'

resdir = paste0(Sys.getenv('CS_HOME'),'/NetworksTerritories/CoevolutionNwTerritories/Results/MacroCoevol/PSE/',resprefix,'/');dir.create(resdir,showWarnings = F)

###
# check convergence (number of patterns)

counts = Reduce(rbind,
lapply(list.files(paste0('pse/',resprefix)),function(s){
  gen = as.numeric(strsplit(substring(s,11),'.',fixed=T)[[1]][1])
  patterns = read.csv(paste0('pse/',resprefix,'/population',gen,'.csv'))
  return(data.frame(generation=rep(gen,3),patterns=c(length(which(patterns$evolution.samples>=1)),length(which(patterns$evolution.samples>=5)),length(which(patterns$evolution.samples>=10))),samples=c(1,5,10)))
})
)

g=ggplot(counts,aes(x=generation,y=patterns,group=samples,color=as.factor(samples)))
g+geom_point()+geom_line()+scale_color_discrete(name=('Samples'))+xlab('Generation')+ylab('Patterns')+stdtheme
ggsave(filename = paste0(resdir,'patterncount.png'),width=20,height=18,units='cm')


###
# plot results
latestgen = max(as.integer(sapply(strsplit(sapply(strsplit(list.files(paste0('pse/',resprefix)),"population"),function(s){s[2]}),".csv"),function(s){s[1]})))
res <- as.tbl(read.csv(paste0('pse/',resprefix,'/population',latestgen,'.csv')))

res$DeltaHierarchyPop=res$hierarchyPopFinal+res$synthRankSize # ! synthranksize is in absolute value
res$DeltaHierarchyCloseness=res$hierarchyClosenessFinal+0.2 # remove average value - if rerun, effectively take the delta
res$gravityDecayF = cut(res$gravityDecay,6)
res$synthRankSizeF = cut(res$synthRankSize,6)
res$gravityGammaF = cut(res$gravityGamma,6)

# Q: why values outside bounds for pop hierarchy? outlier - fail algo ? -> filter
res = res[res$DeltaHierarchyPop>=-0.2,]


for(samples in c(1,5,10)){
  ggsave(
plot=ggpairs(data=res[res$evolution.samples>=samples,],columns = c("hierarchyPopFinal","hierarchyClosenessFinal","rankCorrsPopClosenessFinal"),
        aes(colour=gravityDecayF,alpha=0.4)
)+stdtheme,
filename = paste0(resdir,'scatterobjs_colorgravityDecay_samples',samples,'.png'),width=30,height=30,units='cm')
}

for(samples in c(1,5,10)){
  ggsave(
    plot=ggpairs(data=res[res$evolution.samples>=samples,],columns = c("DeltaHierarchyPop","DeltaHierarchyCloseness","rankCorrsPopClosenessFinal"),
                 aes(colour=gravityDecayF,alpha=0.4)
    )+stdtheme,
    filename = paste0(resdir,'scatterdeltaobjs_colorgravityDecay_samples',samples,'.png'),width=30,height=30,units='cm')
}


for(samples in c(1,5,10)){
  ggsave(
plot=ggpairs(data=res[res$evolution.samples>=samples,],columns = c("hierarchyPopFinal","hierarchyClosenessFinal","rankCorrsPopClosenessFinal"),
        aes(colour=synthRankSizeF,alpha=0.4)
)+stdtheme,
filename = paste0(resdir,'scatterobjs_colorsynthRankSize_samples',samples,'.png'),width=30,height=30,units='cm')
}



for(samples in c(1,5,10)){
  ggsave(
    plot=ggpairs(data=res[res$evolution.samples>=samples,],columns = c("DeltaHierarchyPop","DeltaHierarchyCloseness","rankCorrsPopClosenessFinal"),
                 aes(colour=gravityGammaF,alpha=0.4)
    )+stdtheme,
    filename = paste0(resdir,'scatterdeltaobjs_colorgravityGamma_samples',samples,'.png'),width=30,height=30,units='cm')
}


# pdf figure for paper
samples=10
ggsave(
  plot=ggpairs(data=res[res$evolution.samples>=samples,],columns = c("DeltaHierarchyPop","DeltaHierarchyCloseness","rankCorrsPopClosenessFinal"),
               aes(colour=gravityDecayF,alpha=0.4)
  )+stdtheme,
  filename =paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Docs/Papers/HierarchyCoevolution/chapter/figuresraw/scatterdeltaobjs_colorgravityDecay_samples",samples,".pdf"),width=30,height=30,units='cm')



# values

summary(res[res$DeltaHierarchyPop>=-0.2&res$evolution.samples>=10,])

summary(res[res$DeltaHierarchyPop>=-0.2&res$evolution.samples>=10&res$rankCorrsPopClosenessFinal<0,])


####
# test linear model

# note: this would be useful for the sake of the statistcial analysis to have all individuals and not medians ffrom replications
# -> ! must weight them in the regression ?

### pop hierarchy

summary(lm(data=res,formula=DeltaHierarchyPop~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay))

# weighted
summary(lm(data=res,
           formula=DeltaHierarchyPop~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay,
           weights = evolution.samples
           )
        )
# basically the same - but more accurate


### closeness hierarchy

summary(lm(data=res,formula=DeltaHierarchyCloseness~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay))
# weighted
summary(lm(data=res,formula=DeltaHierarchyCloseness~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay,weights = evolution.samples))

### rank correlation

summary(lm(data=res,formula=rankCorrsPopClosenessFinal~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay))
# weighted
summary(lm(data=res,formula=rankCorrsPopClosenessFinal~synthRankSize+nwThresholdQuantile+nwGmax+nwExponent+gravityWeight+gravityGamma+gravityDecay,weights = evolution.samples))









