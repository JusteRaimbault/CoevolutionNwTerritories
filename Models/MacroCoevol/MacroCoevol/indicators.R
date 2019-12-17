
library(segmented)

# test for double array json import
# library("jsonlite")
# fromJSON('[[1,2],[3,4]]',simplifyMatrix = FALSE)
# => Array[Array[T]] ~ list[c[T]]

#'
#' basic rank size hierarchy
hierarchy <- function(values){
  x = values[values>0]
  reg = lm(data = data.frame(x = log(1:length(x)),y = log(x)),y~x)
  return(list(slope = reg$coefficients[2],rsquared = summary(reg)$adj.r.squared))
}

#'
#' Hierarchy with one breakpoint
segHierarchy <- function(values){
  tryCatch({
    x = values[values>0]
    reg = lm(data = data.frame(x = log(1:length(x)),y = log(x)),y~x)
    seg<-segmented(reg,seg.Z=~x,npsi=1)
    return(list(slope1 = seg$coefficients[2], slope2 = seg$coefficients[3], breakRank = exp(seg$psi[2]),rsquared = summary(seg)$adj.r.squared))
  },error = function(e){return(list(slope1=NA,slope2=NA,breakRank=NA,rsquared=NA))})
}

hierarchies <- function(trajectories){
  return(unlist(sapply(trajectories,hierarchy)))
}

segHierarchies <- function(trajectories){
  return(unlist(sapply(trajectories,segHierarchy)))
}


#'
#' spearman rank corrs between initial/final
rankCorrelation <- function(trajectories){
  cortest = cor.test(x = trajectories[[1]],y=trajectories[[length(trajectories)]],method="spearman")
  # more accurate test not needed as we do not return the pvalue
  #cortest = spearman.test(x = trajectories[[1]],y=trajectories[[length(trajectories)]])
  return(cortest$estimate)
}


#'
#' time serie of rank correlations between variables; assume traj have the same length 
rankCorrelationsTS <- function(traj1,traj2){
  res = c()
  for(i in 1:length(traj1)){
    cortest = cor.test(x = traj1[[i]],y=traj2[[i]],method="spearman")
    res = append(res,cortest$estimate)
  }
  return(res)
}






# trajectories = sapply(1:30,function(i){runif(30)},simplify=F)
# traj2 = sapply(1:30,function(i){runif(30)},simplify=F)
# hierarchies(trajectories)
# segHierarchies(trajectories)
# rankCorrelation(trajectories)
# rankCorrelationsTS(trajectories,traj2)



