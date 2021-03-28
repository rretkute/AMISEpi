# AMISEpi
 
This package contains the code and data to reproduce analysis in Retkute et al. "Integrating geostatistical maps and infectious disease transmission models using adaptive multiple importance sampling" https://www.medrxiv.org/content/10.1101/2020.08.03.20146241v1.

#### Installation

```r
library(devtools)

install_github("rretkute/AMISEpi")
``` 

#### Example use: Ascariasis in Ethiopia

```r
library(AMISEpi)

data(prev)

data(coord)

ans<-run_AMIS_ascaris(prev=prev, n.param=2, NN=rep(1000, 100), delta=5, ESS.R=2000)

plot_AMIS_Ascaris(coord, ans)
``` 

Results  for  Ascaris  lumbricoides  in  Ethiopia.  

![](pkg_img.png)
