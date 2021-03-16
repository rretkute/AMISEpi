#' Prior density for Ascaris
#'
#' AMIS.
#'
#' @param a Proposed value of M
#' @param b Proposed value of k
#' @param fit.v Fitted relationship between M and k
#' @param fit.hess Fitted Hesian matrix
#' @return Probability
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

dprop0_ascaris<-function(a, b, fit.v, fit.hess, fit.inflation.factor){ # (M,k)
  A<-prod(sapply(1:length(a), function(i) dunif(a[i], min=log(0.01), max=log(60))))
  sa<-log(sum(exp(a)))
  B<-prod(sapply(1:length(b), function(i) dnorm(b[i], mean=fit.v[1]+fit.v[2]*exp(sa), sd = fit.inflation.factor*sqrt(c(1, exp(sa))%*%solve(fit.hess,c(1,exp(sa)))))))
  return(A*B)
}
