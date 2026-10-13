
setwd(paste0(Sys.getenv('CS_HOME'),"/NetworksTerritories/CoevolutionNwTerritories/Models/MesoCoevolCA/"))

library(readr)
library(dplyr)

resprefix = '20261013_100155_EXPLORATION'
res <- read_csv(paste0('exploration/',resprefix,'.csv'))
resdir = paste0('results/',resprefix,'/');dir.create(resdir,recursive = T)

TEvariables = c("TECentrRoad","TERoadCentr","TECentreDens","TEDensCentre","TEDensRoad","TERoadDens")

# scatterplot
library(GGally)

ggsave(filename = paste0(resdir,'TE-scatterplot.png'),
  plot = ggpairs(res, columns = TEvariables),
  height = 25,width=25,units = 'cm' 
)


# compute number of reciprocal pairs (rule above median)
res$coevolCentreRoad = as.numeric((res$TECentrRoad > median(res$TECentrRoad))&(res$TERoadCentr>median(res$TERoadCentr)))
res$coevolDensRoad = as.numeric((res$TEDensRoad > median(res$TEDensRoad))&(res$TERoadDens>median(res$TERoadDens)))
res$coevolDensCentre =  as.numeric((res$TEDensCentre > median(res$TEDensCentre))&(res$TECentreDens>median(res$TECentreDens)))

res$coevol = res$coevolCentreRoad+res$coevolDensRoad+res$coevolDensCentre

table(res$coevol)*100/nrow(res)


# lagged correlations

res <- read_csv('exploration/20261013_180315_EXPLORATION_LAGGEDCORRS.csv')

# proportions with median on absolute correlations
res$coevolCentreRoad = as.numeric((abs(res$rhoCentrRoad) > median(abs(res$rhoCentrRoad)))&(abs(res$rhoRoadCentr)>median(abs(res$rhoRoadCentr))))
res$coevolDensRoad = as.numeric((abs(res$rhoDensRoad) > median(abs(res$rhoDensRoad)))&(abs(res$rhoRoadDens)>median(abs(res$rhoRoadDens))))
res$coevolDensCentre =  as.numeric((abs(res$rhoDensCentre) > median(abs(res$rhoDensCentre)))&(abs(res$rhoCentreDens)>median(abs(res$rhoCentreDens))))

res$coevol = res$coevolCentreRoad+res$coevolDensRoad+res$coevolDensCentre

table(res$coevol)*100/nrow(res)





