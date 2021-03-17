#' Run AMIS for Ascaris
#'
#' Parameter estimation for Adcaris prevalence distributions.
#'
#' @param prev Matrix with observed prevalences: row for a separate pixel
#' @param n.param Number of parameters
#' @param NN Vector containing number of samples for each iteration
#' @param delta Threshold for prevalence
#' @param ESS.R Minimum Effective Sample Size
#' @return Matrix with sampled vectpors and weigts for each pixel
#' @import mnormt
#' @import mclust
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#'
#' @examples
#' run_AMIS_ascaris(prev=prev, n.param=2, NN=rep(1000, 100), delta=5, ESS.R=2000)
#'
#' @export
#'

run_AMIS_ascaris<-function(prev, n.param, NN, delta, ESS.R){
  start_time <- Sys.time()
  TT<-length(NN)# Max number of iterations
  n.pixels<-nrow(prev) # Number of pixels

  #  Set up for mixture function
  proposal=mvtComp(df=3); mixture=mclustMix();
  dprop <- proposal$d
  rprop <- proposal$r

  # Fitted dependancy between W and k
  fit.v<-c(0.33371009, 0.01719462)
  fit.inflation.factor<-5
  fit.hess<-matrix(0, nrow = 2, ncol = 2)
  fit.hess[[1,1]]<-5138.97
  fit.hess[[1,2]]<-49499.4
  fit.hess[[2,1]]<-49499.40
  fit.hess[[2,2]]<-677831.0

  param<-matrix(NA, ncol=n.param+1, nrow=0)
  Sigma <- list(NA)
  Mean<-list(NA)
  PP<-list(NA)
  GG<-list(NA)

  # Iteration 1
  it<-1
  cat(c("Started iteration ", it,"\n"))

  for(i in 1:NN[it]){
    M<-c(); k<-c()
    while(length(k)==0){
      a<-runif(n.param-1, min=log(0.01), max=log(60))
      sa<-sum(a)
      b<-rnorm(1, mean=fit.v[1]+fit.v[2]*exp(sa), sd = fit.inflation.factor*sqrt(c(1, exp(sa))%*%solve(fit.hess,c(1,exp(sa)))))
      if(b>0 & b<3){
        M<-c(M, a)
        k<-c(k, b)
      }
    }
    pr<- ascaris_model(M,k)
    param<-rbind(param, c(M,k,100*pr))
  }
  sim<-param[,n.param+1]

  #Calculate ESS and weighst
  w1<-rep(1, length(sim))
  tmp<-calculate_ESS(prev, sim, w1, delta)
  ess<-tmp$ess
  WW<-tmp$WW
  rm(tmp)

  # Iterations 2+
  stop<-0
  while(stop==0){
    it<-it+1
    cat(c("Started iteration ", it,"\n"))
    wh<-which(ess>=ESS.R)
    W1<-WW; W1[wh,]<-0
    wh<-which(ess<ESS.R)
    w1<- c(colSums(W1))
    J<-sample(1:sum(NN[1:(it-1)]), NN[it], prob= w1, replace=TRUE)
    xx<-param[J,1:n.param]
    clustMix <- mixture(xx)
    G <- clustMix$G
    cluster <- clustMix$cluster
    ### Components of the mixture
    ppt <- clustMix$alpha
    muHatt <- clustMix$muHat
    varHatt <- clustMix$SigmaHat
    GG[[it-1]]<-G
    G1<-0; G2<-G
    if(it>2) {
      G1<-sum(sapply(1:(it-2), function(a) GG[[a]]))
      G2<-sum(sapply(1:(it-1), function(a) GG[[a]]))
    }
    for(i in 1:G){
      Sigma[[i+G1]] <- varHatt[,,i]
      Mean[[i+G1]] <- muHatt[i,]
      PP[[i+G1]]<-ppt[i]
    }

    # Draw new parameters and calculate prevalences
    for(i in 1:NN[it]){
      while(nrow(param)<sum(NN[1:it])){
        compo <- sample(1:G,1,prob=ppt)
        x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
        new.param<-as.numeric(x1)
        M<-new.param[1:(n.param-1)];
        k<-new.param[n.param]
        if(dprop0_ascaris(M,k,fit.v, fit.hess, fit.inflation.factor)>0 & k>0 & k<3){
          pr<- ascaris_model(M,k)
          if(!(is.na(pr))){
            param<-rbind(param, c(M,k,100*pr))
          }
        }
      }
    }
    sim<-param[1:sum(NN[1:it]), n.param+1]
    w1 <- sapply(1:sum(NN[1:it]), function(b)  dprop0_ascaris(param[b,1:(n.param-1)], param[b,n.param], fit.v, fit.hess, fit.inflation.factor))/
      (sapply(1:sum(NN[1:it]), function(b)  dprop0_ascaris(param[b,1:(n.param-1)], param[b, n.param], fit.v, fit.hess, fit.inflation.factor) +
                sum(sapply(1:G2, function(a) PP[[a]] * dprop(param[b,1:n.param],mu= Mean[[a]], Sig=Sigma[[a]])))))

    #Calculate ESS and weighst
    tmp<-calculate_ESS(prev, sim, w1, delta)
    ess<-tmp$ess
    WW<-tmp$WW
    rm(tmp)

    if(min(ess)>=ESS.R) stop<-1
    if(it>= TT) stop<-1
  }
  end_time <- Sys.time()
  total_time<-difftime(end_time, start_time, units ="hours")
  return(list(param=param, WW=WW, ESS=ess, total_time=total_time))
}
