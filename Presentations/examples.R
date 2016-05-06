
library(devtools)
dir.create("libraries")
.libPaths("libraries")
devtools::install_github("markmfredrickson/RItools@randomization-distribution")

library(RItools)

Z <- c(0,1,0,1)
Y <- c(16,22,7,14)
Om <- matrix(0,ncol=choose(4,2),nrow=length(Z))
whotrted <- combn(1:4,2)
for(i in 1:choose(4,2)){ Om[cbind(whotrted[,i],i)]<-1 }
meandifftz <- function(y,z){ mean(y[z==1]) - mean(y[z==0]) }
thedist<-apply(Om,2, function(z){ meandifftz(Y,z) }) 
theobs <- meandifftz(Y,Z)
mean(thedist >= theobs)




