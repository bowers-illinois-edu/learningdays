rm(list = ls())

##IF YOU HAVE NOT INSTALLED BLOCKTOOLS:
install.packages("blockTools")
library("blockTools")

##########################################
##randomization assignment R code
##########################################


##Administration
setwd(COPY AND PASTE YOUR OWN PATHWAY TO YOUR EXPERIMENT>ASSIGNMENTS>RANDOMIZATION FOLDER HERE)
randomize<-read.csv("randomize.csv", header = TRUE)
TestScore<-runif(nrow(randomize), min = 35, max = 100)
randomize$TestScore<-TestScore
randomize$town<-rep(1:25, each = 5 , len = 123)

#Question 3(a)
#number of schools by gender (note; there are three gender classifications!)
table(randomize$Gender)
#number of schools by language type
table(randomize$Language)


by(randomize$TestScore, list(randomize$Gender), FUN = mean)
by(randomize$TestScore, list(randomize$Language), FUN = mean)

#Question 3(b)
set.seed (21042015)
#simple randomization over all schools
randomize$treat3b<-rbinom(n=length(randomize$SchoolID), size=1, p=0.5)
#number of schools in treatment and control group by gender
table(randomize$treat3b, randomize$Gender)
#number of schools in treatment and control group by language
table(randomize$treat3b, randomize$Language)

#Question 3(c)
set.seed (04212015)
#simple randomization over all schools
randomize$treat3c<-rbinom(n=length(randomize$SchoolID), size=1, p=0.5)
#number of schools in treatment and control group by gender
table(randomize$treat3c, randomize$Gender)
#number of schools in treatment and control group by language
table(randomize$treat3c, randomize$Language)

#Question 3(e)
set.seed(21042015)
#create a numeric variable that corresponds to Gender (since gender is a character variable)
randomize$gender<-as.integer(randomize$Gender)
#block randomize by gender with two treatment types (treatment & control)
treat3e <-assignment(block(data=randomize, block.vars=c("gender"),n.tr=2, id.vars=c("SchoolID")))
#create a variable in which treatment == 1 and control == 0 based on above block randomization
randomize$treat3e <- ifelse(
  # Is the person assigned to treatment by blocktools?
  test = randomize$SchoolID %in% 
    as.numeric(as.character(as.data.frame(
      treat3e[[1]])[,"X1.Treatment.2"]
    )),
  # If yes, return 1
  yes = 1,
  # If no, return 0
  no = 0
)

#total observations in treatment and control groups
table(randomize$treat3e)
#observations in treatment and control groups by gender
table(randomize$treat3e, randomize$Gender)

#Question 3(f)
set.seed(4212015)
#block randomize by gender, with two treatment types (treatment & control)
treat3f<-block(data=randomize, block.vars=c("gender"),n.tr=2, id.vars=c("SchoolID"))
#assign treatment using the block randomization from the previous command
blocktools.treat3f<-assignment(treat3f)

#generate a variable that treatment = 1 and control = 0 based on the block assignment in previous command
randomize$treat3f <- ifelse(
  # Is the person assigned to treatment by blocktools?
  test = randomize$SchoolID %in% 
    as.numeric(as.character(as.data.frame(
      blocktools.treat3f[[1]])[,"X1.Treatment.2"]
    )),
  # If yes, return 1
  yes = 1,
  # If no, return 0
  no = 0
)

#summarize treatment numbers by gender
table(randomize$treat3f, randomize$Gender)

#Question 3(g)
#create a variable that calculates many times was a school assigned to treatment over the 4 randomizations
randomize$treat.assign<-rowSums(randomize[,c(6,7,9,10)])
#create a variable that calculates the proportion of times each school was assigned to treatment
randomize$treat.assign.prop<-randomize$treat.assign/4
#Table that shows how many schools were assigned with each probability
table(randomize$treat.assign.prop)
#average probability of assignment to treatment
mean(randomize$treat.assign.prop)

#Question 3(i)
set.seed (21042015)
#how many towns are there
max(randomize$town)

#Use simple randomization (think coin flipping) for each town to assign to treatment or control
town.treat  <- rbinom(n=length(unique(randomize$town)), size=1, p=0.5)
we multiply town.treat (a vector of 0 and 1's) with a vector of the town ID values (1 through 25) so that what is
  #returned is a vector where if treated, the town ID is returned but if in control a zero is returned.
town.treat  <- 1:25*town.treat

#The vector created above is only a vector of 25 (the number of towns) but we have 124 schools,
#the next step is to determine if a school's town is in the treatment group or not, and assign treatment at the school
#level in the main dataframe.

randomize$cluster.treat<- ifelse(
  # This allows us
  #to test at the school observation level whether a school's town is treated or not.
  test = randomize$town %in% town.treat ,
  # If yes, return 1
  yes = 1,
  # If no, return 0
  no = 0
)

#total number of schools treated
sum(randomize$cluster.treat)
#total number of treated TOWNS
length(unique(randomize$town[randomize$cluster.treat==1]))
#Proportion of treated TOWNS (treated towns/total # of towns)
length(unique(randomize$town[randomize$cluster.treat==1]))/length(unique(randomize$town))

#Question 3(j)
u<-sample(x=1:length(unique(randomize$town)), size=13, replace = FALSE)
randomize$cluster.control<-ifelse(
  #is the town one of the randomly selected to be treated
  test = randomize$town %in% u,
  #if yes, return 1
  yes = 1, 
  #if no, return 0
  no = 0)
#proportion of treated TOWNS (treated towns/total # of towns)
length(unique(randomize$town[randomize$cluster.control == 1]))/length(unique(randomize$town))

#Question 3(k)
set.seed (21042015)

t<-sample(x=randomize$SchoolID, size = 40, replace = F)
randomize$limited.treat<-ifelse(
  #is the school ID one of the randomly selected to be treated
  test = randomize$SchoolID %in% t, 
  #if yes, return 1
  yes = 1 , 
  #if no, return 0
  no = 0)
table(randomize$limited.treat)
table(randomize$limited.treat, randomize$Gender)