#' Proposal for mixtur distribution
#'
#' Mechanistic model of Ascaris prevalence.
#'
#' @param df Number of degrees of freedom, df=3 for t Student distribution.
#' @return Distribution
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
mvtComp<-function(df=3){
  list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
    dmt(xx,mean=mu,S=Sig,df=df,log=log)
  },"r"=function(n=1,mu=0,Sig=1){
    rmt(n,mean=mu,S=Sig,df=df)
  })
}


