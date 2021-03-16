#' Fit clusters of mixtures
#'
#' AMIS.
#'
#' @param xx Sampled parameters
#' @param G Sequence of clusters
#' @return Parameters of mixture
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

fitMclust<-function(xx, modelName="VVV", G= G){
  options(warn=-1)
  control <- emControl(eps=sqrt(.Machine$double.eps))
  n <- nrow(xx)
  p <- ncol(xx)
  clustering <-Gout <- BIC <- NA
  if (G[1] == 1) {
    clustering <- mvn(modelName = modelName, data = xx)
    BIC <- bic(modelName=modelName,loglik=clustering$loglik,n=n,d=p,G=1)
    Gout <- 1
    G <- G[-1]
  }
  if (p != 1) {
    if (n > p) {
      hcPairs <- hc(modelName="VVV",data=xx)
    }else {
      hcPairs <- hc(modelName="EII",data=xx)
    }
  }else hcPairs <- NULL
  if (p > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)

  for (g in G) {

    if (p > 1 || !is.null(hcPairs)) {
      cl <- clss[, as.character(g)]
    }else {
      cl <- .qclass(data[subset], as.numeric(g))
    }
    z <- unmap(cl, groups = 1:max(cl))
    new <- me(modelName=modelName,data=xx,z=z,control=control)
    if(!is.na(new$loglik)){
      BICnew <- bic(modelName=modelName,loglik=new$loglik,n=n,d=p,G=g,equalPro=control$equalPro)
      if(is.na(BIC)){
        clustering <- new
        BIC <- BICnew
        Gout <- g
      }else{
        if(BICnew>BIC){
          clustering <- new
          BIC <- BICnew
          Gout <- g
        }
      }
    }
  }
  options(warn=0)
  return(c(clustering,G=Gout))
}




