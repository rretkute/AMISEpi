#' Plot Ascaris AMIS results
#'
#'
#' @param N0 Initial population size
#' @param f0 Fraction of initially infeceted
#' @param t0 Year HIV outreak started
#' @return Prevalence
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
HIV_model<-function(N0, f0, t0, phi, mu, lambda0, dt, TT){
  x<-data.frame(X=N0* (1-f0), Y=N0*f0, Z=0, t=t0)
  t<-dt+t0
  while(t<=TT){
    new_x<-modelHIV(t, x)
    x<-rbind(x, data.frame(X=new_x[1], Y=new_x[2], Z=new_x[3], t=t))
    t<-t+dt
  }
  return(x)
}

# The population at timet is divided into three groups,
# a not-at-risk group X(t),
# an at-risk group Z(t) and
# an in-fected group Y(t).

modelHIV <- function(t, x, phi, mu, lambda0){
  n<-nrow(x)
  X <- x[n,1]
  Z <- x[n,2]
  Y <- x[n,3]
  tt<-x[,4]
  N <- X+Z+Y
  wh<-which(tt<=t-15)
  if(length(wh)>0){
    E<-x[max(wh),1]+ x[max(wh),2] + x[max(wh),3]
  } else {
    E<-x[1,1]+x[1,2]+x[1,3]
  }
  dX <- (1-F1(X,N,f0,phi))*b*E-mu*X
  if(n==1){
    dZ <- F1(X,N,f0,phi)*b*E-(mu+r*Y/N +lambda0)*Z
    dY <- (r*Y/N +lambda0)*Z
  } else {
    dZ <- F1(X,N,f0,phi)*b*E-(mu+r*Y/N)*Z
    dY <- (r*Y/N)*Z -dt*sum((r*x[2:n,3]/(x[2:n,1]+x[2:n,2]+x[2:n,3]))*x[2:n,2]*g(t-tt[2:n])) - dt*(r*x[1,3]/(x[1,1]+x[1,2]+x[1,3])+lambda0)*x[1,2]*g(t-tt[1])
  }
  res <- pmax(0,c(dt*dX+X, dt*dZ+Z, dt*dY+Y))
  return(res)
}

F1 <- function(X, N, f0, phi) {
  x<-X/N-(1-f0)
  return(exp(phi*x)/(exp(phi*x)-1+1/f0))
}

# Death rate without treatments
g<-function(t){
  dweibull(t, shape=2.4, scale = 10.5, log = FALSE)
}

