# AMISEpi
 
This package contains the code and data to reproduce analysis in Retkute et al. "Integrating geostatistical maps and infectious disease transmission models using adaptive multiple importance sampling".

Run following:

```r
library(AMISEpi)

data(prev)

data(coord)

ans<-run_AMIS_ascaris(prev=prev, n.param=2, NN=rep(1000, 100), delta=5, ESS.R=2000)

plot_AMIS_Ascaris(coord, ans)
``` 


Figure. Results  for  Ascaris  lumbricoides  in  Ethiopia.  

![](pkg_img.png)
