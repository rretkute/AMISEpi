#' Ascaris model
#'
#' Mechanistic model of Ascaris prevalence.
#'
#' @param M Number of worms
#' @param k Work clumping coeficient
#' @return Prevalence
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
ascaris_model<-function(M,k){
  sM<-sum(exp(M))
  p<-(1-(1+sM/k)^(-k))
  return(p)
}


