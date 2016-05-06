rm(list=ls())

setwd(COPY AND PASTE YOUR OWN PATHWAY TO YOUR EXPERIMENT>ASSIGNMENTS>HYPOTHESIS FOLDER HERE)
data  <- read.csv("hypothesis.csv", header = T)

#Question 1(c)
#Whole Sample
mean(data$testscore)
sd(data$testscore)

#Treatment
tmu <- with(data, mean(testscore[treat==1]))
tmu
with(data, sd(testscore[treat==1]))


#Control
cmu <- with(data, mean(testscore[treat==0]))
cmu
with(data, sd(testscore[treat==0]))


#Question 1(d)
with(data,t.test(testscore[treat==1], testscore[treat==0]))


