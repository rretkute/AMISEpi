#' Calculates Effective Sample Sizes and weights
#'
#' AMIS.
#'
#' @param prev Matrix with observed prevalences
#' @param sim Vector of simulated prevalences
#' @param w1 Weight matrix 1
#' @param delta Threshold value to compare similarity between prevalences
#' @return Vector with ESS and weight mnatrix
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

calculate_ESS<-function(prev, sim, w1, delta){

  n.pixels<-nrow(prev)
  n.param<-length(sim)

  WW.cnt<-matrix(0, nrow=n.pixels, ncol=n.param)
  g<-rep(0, n.param)
  for(j in 1:n.param){
    x1<-H(delta/2-abs(prev-sim[j]))
    x2<-matrix(x1, ncol=ncol(prev), nrow=nrow(prev))
    x3<-rowSums(x2)
    WW.cnt[,j]<-x3
    g[j]<-sum(w1[which(abs(sim-sim[j])<=delta/2)])/sum(w1)
  }

  ess<-c()
  WW<-matrix(0, nrow=n.pixels, ncol=n.param)
  for(i in 1:n.pixels){
    f<-WW.cnt[i,]
    ww<-w1*(f/g)
    if(sum(ww) >0)
      ww<-ww/sum(ww)
    WW[i,]<-ww
    if(sum(ww)>0) {
      www<-(sum((ww)^2))^(-1)
    } else {
      www<-0
    }
    ess[i]<- www
  }

  return(list(ess=ess, WW=WW))
}
