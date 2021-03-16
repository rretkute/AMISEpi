#' Runs OpenMalaria model
#'
#'
#' @param ibra Infectious bites rate per year
#' @param sigma2i Host inter-variability
#' @return Prevalence
#' @import xml2
#' @import methods
#' @import plyr
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
Malaria_model<-function(ibra, sigma2i){
  d.obs.1.s<-round(39*365/10)
  d.obs.1.e<-round(40*365/10)
  d.obs.2.s<-round(44*365/10)
  d.obs.2.e<-round(45*365/10)
  d.obs.3.s<-round(49*365/10)
  d.obs.3.e<-round(50*365/10)
  x <- read_xml("scenario_7B.xml")
  y <- xml_find_all(x, ".//seasonality")
  xml_set_attr(y, "annualEIR", paste(exp(ibra)))
  y<-xml_find_all(x, ".//parameter")
  xml_set_attr(y[6], "value", sigma2i)
  write_html(x, 'scenario_8B.xml', options = "format")
  system('./openMalaria --scenario scenario_8B.xml')
  output<-read.table("output.txt", header=F)
  time<-unique(output[,1])
  NN<-sapply(1:length(time), function(a) output[output[,1]==time[a] & output[,3]==0 & output[,2]==2,4])
  II<-sapply(1:length(time), function(a) output[output[,1]==time[a] & output[,3]==1 & output[,2]==2,4])
  pr1<-mean(II[d.obs.1.s:d.obs.1.e]/NN[d.obs.1.s:d.obs.1.e])
  pr2<-mean(II[d.obs.2.s:d.obs.2.e]/NN[d.obs.2.s:d.obs.2.e])
  pr3<-mean(II[d.obs.3.s:d.obs.3.e]/NN[d.obs.3.s:d.obs.3.e])
  return(c(pr1, pr2,pr3))
}
