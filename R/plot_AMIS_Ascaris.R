#' Plot Ascaris AMIS results
#'
#'
#' @param coord Coordinates od pixels in Ethiopia
#' @param p Fraction of initially infeceted
#' @param t0 Year HIV outreak started
#' @return Prevalence
#' @import ggplot2
#' @import RColorBrewer
#' @import gridExtra
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
plot_AMIS_Ascaris<-function(coord, ans){
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  data1<-coord
  data1$Prevalence<-data1$prevalence/100
  hm1 <- ggplot(data1,
                aes(x = latitude, y = longitude, fill = Prevalence)) + geom_tile()
  hm1 <- hm1 + scale_fill_gradientn(limits = c(0, 1), name = "Mean prevalence", colours = myPalette(100))
  hm1 <- hm1 + scale_x_continuous(expand = c(0,0), breaks = seq(0.05,0.25,0.05))
  hm1 <- hm1 + scale_y_continuous(expand = c(0, 0))
  hm1 <- hm1 + theme_bw() + theme(legend.position="bottom")  + coord_fixed(ratio = 1)
  f1<-hm1 + ggtitle("(a)")
  par<-ans$param
  par<-par[c(1,2,ncol(par)),]
  n.iter<-nrow(par)/1000-1
  colnames(par)<-c("logW","k","Prevalence")
  par<-as.data.frame(par)
  par$Prevalence<-par$Prevalence/100
  par$iteration<-1
  for(i in 1:n.iter) par$iteration[(i*1000+1):((i+1)*1000)]<-i+1
  par$iteration<-as.factor(par$iteration)
  f2<- ggplot(par, aes(x=Prevalence, col=I('white'), fill = forcats::fct_rev(iteration))) +
    geom_histogram(bins=50) +xlab("Prevalence") + ylab("") + labs(fill = "") +
    ggtitle("(b)") + theme_bw() + theme(legend.position="right")
  f3<-ggplot(par, aes(logW, k)) + geom_point(aes(colour = Prevalence)) +
    scale_colour_gradientn(limits = c(0,1), name = "Prevalence",
                           colours =  myPalette(100)) +xlab("log(W)") + ylab("k") + ggtitle("(c)") +
    theme_bw() + theme(legend.position="right")
  f4 <- ggplot(par, aes(logW, k)) + geom_bin2d(bins=50) +
    xlab("log(W)") + ylab("k") + ggtitle("(d)") +
    scale_fill_gradient(low="gray", high = "black")  +
    theme_bw() + theme(legend.position="right")
  return(grid.arrange(f1, f2, f3, f4, ncol = 2))
}
