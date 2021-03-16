#' Runs EPIONCHO model under two alternative control scenarious
#'
#'
#' @param abr Annual biting rate
#' @param k Host variability contant
#' @return Prevalence
#' @import Rcpp
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
Oncho_model<-function(abr, k){
  sourceCpp(file="EPIONCHOv2.cpp")
  load(file="theta.Rdata")
  thetaA<-theta; thetaB<-theta
  thetaB$ntr5<-1; thetaB$ftrt1<-0.25; thetaB$ftrt2<-0.25; thetaB$ftrt3<-0.25; thetaB$ftrt4<-0.25; thetaB$ftrt5<-0.25
  thetaA<-theta;
  thetaA$ABR<-exp(abr)
  thetaA$kW0<-k
  out1 <- runEPIONCHO(theta = as.double(thetaA), itervtn = 1)
  sourceCpp(file="EPIONCHOv2.cpp")
  load(file="theta.Rdata")
  thetaB<-theta
  thetaB$ABR<-exp(new.param[1])
  thetaB$kW0<-new.param[2]
  thetaB$cov1<-0.65; thetaB$cov2<-0.65; thetaB$cov3<-0.65; thetaB$cov4<-0.65; thetaB$cov5<-0.65
  out2 <- runEPIONCHO(theta = as.double(thetaB), itervtn = 1)
  return(list(t1=out1$time, prev1=out1$Mp5, time1=out2$time, prev2=out2$Mp5))
}
