#' Mixture distribution
#'
#' AMIS.
#'
#' @param G Sequence of clusters
#' @return Parameters of mixture
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
mclustMix<-function(G=1:10){
  if (any(as.numeric(G)<=0)) stop("G must be positive")
  function(xx){
    clustering <- fitMclust(xx,modelName="VVV",G= G)
    G <- clustering$G
    if(G==1) clustering$parameters$pro <- 1
    return(list(alpha=clustering$parameters$pro, muHat=t(clustering$parameters$mean), SigmaHat=clustering$parameters$variance$sigma,G=G,cluster=clustering$classification))
  }
}


