
set.seed(20150313)
## Bias refers to a relationship between the repeated operation of a procedure and a truth. So we have to invent a truth.
newsdf$latenty0<-rnorm(nrow(newsdf))
newsdf$fakey0<-as.numeric(newsdf$latenty0 >= qnorm(.55,lower.tail=FALSE) )
prop.table(table(newsdf$fakey0))
## newsdf$fakey0<-sample(rep(c(0,1),c(nrow(newsdf)-numhlpers,numhlpers)))
trueATE<-.25 ## posit a true average treatment effect
## In the context of a binary outcome such a treatment effect is a difference of proportions
## that is, we should change 25\% of the 0s in fakey0 to 1.
newsdf$latenty1<-newsdf$latenty0+trueATE
newsdf$fakey1<-as.numeric(newsdf$latenty1 > qnorm(.8,mean=mean(newsdf$latenty1),lower.tail=FALSE))

newsdf$obsy<-with(newsdf, postbomb*fakey1+(1-postbomb)*fakey0 ) ## what we observe

## calculate the true ATE and the $\hat{\bar{\tau}}$
trueATEfake<-with(newsdf,mean(fakey1)-mean(fakey0))
trueTotal<-with(newsdf,sum(fakey1))
trueDelta<-with(newsdf, log( mean(fakey1)/(1-mean(fakey1))) - log( mean(fakey0)/(1-mean(fakey0))))
## true Logit?
## estimate the true ATE using the data that we would observe in this fake experiment
estATEfake<-coef(lm(obsy~postbomb,newsdf))["postbomb"] ## same as a mean difference on obsy
estTotal<-with(newsdf,mean(obsy[postbomb==1])*length(obsy))
estDelta1<-coef(glm(obsy~postbomb,newsdf,family=binomial(link="logit")))[["postbomb"]]
estDelta2<-with(newsdf, log( mean(obsy[postbomb==1])/(1-mean(obsy[postbomb==1]))) -
		 log( mean(obsy[postbomb==0])/(1-mean(obsy[postbomb==0])))
	      )
## Notice that estDelta1 and estDelta2 are the same.

# define a function which reveals a difference in observed outcome and calculates
## estimates of the ATE given a different treatment vector
makeNewObsyAndEst<-function(thez){
    newobsy<-with(newsdf, thez*fakey1+(1-thez)*fakey0 )
    lmATE<-coef(lm(newobsy~thez))[["thez"]]
    totalEffect<-mean(newobsy[thez==1])*length(newobsy)
    logitglm<-glm(newobsy~thez,family=binomial(link="logit"))
    haty0<-predict(logitglm,newdata=data.frame(thez=0),type="response")
    haty1<-predict(logitglm,newdata=data.frame(thez=1),type="response")
    logitDelta<-log( mean(haty1)/(1-mean(haty1))) - log( mean(haty0)/(1-mean(haty0)))
    logitglmATE<-haty1-haty0
    logitcoef<-coef(logitglm)[["thez"]]
    return(c(lmATE=lmATE,totalTE=totalEffect,logitcoef=logitcoef,logitglmATE=logitglmATE,logitDelta=logitDelta))
}

## Does the pair of functions do what we want them to do?
makeNewObsyAndEst(sample(newsdf$postbomb))

nsims<-10000
## For many of the possible ways to run the experiment, calculate this mean difference
### The slow way:
## dist.sample.est<-replicate(nsims,makeNewObsyAndEst(sample(newsdf$postbomb)))

### The fast way uses all of the cores on your unix-based machine (mac or linux):
require(parallel)
ncores<-detectCores()
system.time(
dist.sample.est<-simplify2array(
                                mclapply(1:nsims,function(i){
                                         makeNewObsyAndEst(sample(newsdf$postbomb))
                                 },mc.cores=ncores)
                                )
)

str(dist.sample.est)
apply(dist.sample.est,1,summary)

## Compare to
trueATEfake
trueTotal
trueDelta

## And recall that we have simulation error on the order of 1/sqrt(nsims)
SEsims<-apply(dist.sample.est,1,function(x){ sqrt(var(x)/nsims) })



