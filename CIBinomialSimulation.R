Ns = c(10:210)
p1 <- 0.1; p2 <- 0.2; p3 <- 0.3; p4 <- 0.4; p5 <- 0.5
nrep <- 10^5

################
##### WALD #####
################

# Calculate Binomial CI: Wald formula

liw <- function(x, conf=0.95) { 
  n <- length(x)
  pe <- mean(x)
  z = abs(qnorm((1-conf)/2))
  return((pe)  - z*sqrt(pe*(1-pe)/n))
}

lsw <- function(x, conf=0.95) { 
  n <- length(x)
  pe <- mean(x)
  z = abs(qnorm((1-conf)/2))
  return((pe)  + z*sqrt(pe*(1-pe)/n))
}

# Simulation results: Wald
liwres1 = list(); liwres2 = list(); liwres3 = list(); liwres4 = list(); liwres5 = list()
lswres1 = list(); lswres2 = list(); lswres3 = list(); lswres4 = list(); lswres5 = list()
exw1=list(); exw2=list(); exw3=list(); exw4=list(); exw5=list()
totw1=NULL; totw2=NULL; totw3=NULL; totw4=NULL; totw5=NULL

set.seed(1234)

for(i in 1:length(Ns)){
  n <- Ns[i]
  dat1 <- rbinom(nrep*ns,1,p1)
  dat2 <- rbinom(nrep*ns,1,p2)
  dat3 <- rbinom(nrep*ns,1,p3)
  dat4 <- rbinom(nrep*ns,1,p4)
  dat5 <- rbinom(nrep*ns,1,p5)
  x1 <- matrix(dat1, ncol=nrep)
  x2 <- matrix(dat2, ncol=nrep)
  x3 <- matrix(dat3, ncol=nrep)
  x4 <- matrix(dat4, ncol=nrep)
  x5 <- matrix(dat5, ncol=nrep)
  liwres1[[i]] <- apply(x1, 2, liw)
  lswres1[[i]] <- apply(x1, 2, lsw)
  liwres2[[i]] <- apply(x2, 2, liw)
  lswres2[[i]] <- apply(x2, 2, lsw)
  liwres3[[i]] <- apply(x3, 2, liw)
  lswres3[[i]] <- apply(x3, 2, lsw)
  liwres4[[i]] <- apply(x4, 2, liw)
  lswres4[[i]] <- apply(x4, 2, lsw)
  liwres5[[i]] <- apply(x5, 2, liw)
  lswres5[[i]] <- apply(x5, 2, lsw)
  exw1[[i]] <- liwres1[[i]] <= p1 & p1 <= lswres1[[i]]
  exw2[[i]] <- liwres2[[i]] <= p2 & p2 <= lswres2[[i]]
  exw3[[i]] <- liwres3[[i]] <= p3 & p3 <= lswres3[[i]]
  exw4[[i]] <- liwres4[[i]] <= p4 & p4 <= lswres4[[i]]
  exw5[[i]] <- liwres5[[i]] <= p5 & p5 <= lswres5[[i]]
  totw1[i] = sum(exw1[[i]]==T)/nrep
  totw2[i] = sum(exw2[[i]]==T)/nrep
  totw3[i] = sum(exw3[[i]]==T)/nrep
  totw4[i] = sum(exw4[[i]]==T)/nrep
  totw5[i] = sum(exw5[[i]]==T)/nrep
}

resuls.wald <- cbind(totw1, totw2, totw3, totw4, totw5); results.wald

rm(list = ls(pattern = "^liwres"))
rm(list = ls(pattern = "^lswres"))
rm(list = ls(pattern = "^exw"))
rm(dat1, dat2, dat3, dat4, dat5)
rm(x1, x2, x3, x4, x5)

################
#### SCORE #####
################

# Calculate Binomial CI: Score formula

lis <- function(x, conf=0.95){
  n <- length(x)
  pe <- mean(x)
  z = abs(qnorm((1-conf)/2))
  reslis <- (2*n*pe+(z^2)-z*(sqrt((4*n*pe*(1-pe))+z^2)))/ (2*(n+z^2))
  return(reslis)
}

lss <- function(x, conf=0.95){
  n <- length(x)
  pe <- mean(x)
  z = abs(qnorm((1-conf)/2))
  reslss <- (2*n*pe+(z^2)+z*(sqrt((4*n*pe*(1-pe))+z^2)))/ (2*(n+z^2))
  return(reslss)
}

# Simulation results: Score
lisres1 = list(); lisres2 = list(); lisres3 = list(); lisres4 = list(); lisres5 = list()
lssres1 = list(); lssres2 = list(); lssres3 = list(); lssres4 = list(); lssres5 = list()
exs1=list(); exs2=list(); exs3=list(); exs4=list(); exs5=list()
tots1=NULL; tots2=NULL; tots3=NULL; tots4=NULL; tots5=NULL

set.seed(1234)

for(i in 1:length(Ns)){
  ns <- Ns[i]
  dat1 <- rbinom(nrep*ns,1,p1)
  dat2 <- rbinom(nrep*ns,1,p2)
  dat3 <- rbinom(nrep*ns,1,p3)
  dat4 <- rbinom(nrep*ns,1,p4)
  dat5 <- rbinom(nrep*ns,1,p5)
  x1 <- matrix(dat1, ncol=nrep)
  x2 <- matrix(dat2, ncol=nrep)
  x3 <- matrix(dat3, ncol=nrep)
  x4 <- matrix(dat4, ncol=nrep)
  x5 <- matrix(dat5, ncol=nrep)
  lisres1[[i]] <- apply(x1, 2, lis)
  lssres1[[i]] <- apply(x1, 2, lss)
  lisres2[[i]] <- apply(x2, 2, lis)
  lssres2[[i]] <- apply(x2, 2, lss)
  lisres3[[i]] <- apply(x3, 2, lis)
  lssres3[[i]] <- apply(x3, 2, lss)
  lisres4[[i]] <- apply(x4, 2, lis)
  lssres4[[i]] <- apply(x4, 2, lss)
  lisres5[[i]] <- apply(x5, 2, lis)
  lssres5[[i]] <- apply(x5, 2, lss)
  exs1[[i]] <- lisres1[[i]] <= p1 & p1 <= lssres1[[i]]
  exs2[[i]] <- lisres2[[i]] <= p2 & p2 <= lssres2[[i]]
  exs3[[i]] <- lisres3[[i]] <= p3 & p3 <= lssres3[[i]]
  exs4[[i]] <- lisres4[[i]] <= p4 & p4 <= lssres4[[i]]
  exs5[[i]] <- lisres5[[i]] <= p5 & p5 <= lssres5[[i]]
  tots1[i] = sum(exs1[[i]]==T)/nrep
  tots2[i] = sum(exs2[[i]]==T)/nrep
  tots3[i] = sum(exs3[[i]]==T)/nrep
  tots4[i] = sum(exs4[[i]]==T)/nrep
  tots5[i] = sum(exs5[[i]]==T)/nrep
}

results.score <- cbind(tots1, tots2, tots3, tots4, tots5); results.score

rm(list = ls(pattern = "^lisres"))
rm(list = ls(pattern = "^lssres"))
rm(list = ls(pattern = "^exs"))
rm(dat1, dat2, dat3, dat4, dat5)
rm(x1, x2, x3, x4, x5)

#################################
#### EXACT/ Clopper Pearson ####
#################################

lie <- function(x, conf=0.95) { 
  n <- length(x)
  alfa <- (1-conf)
  exitos <- sum(x)
  if (exitos == 0) {
    reslie <- 0
  }
  else if (exitos == n) {
    reslie <- (alfa/2)^(1/n)
  }
  else {
    reslie <- 1/(1+(n-exitos+1)/(exitos*qf(alfa/2, 2*exitos, 2*(n-exitos+1))))}
  return(reslie)
}

lse <- function(x, conf=0.95) { 
  n <- length(x)
  alfa <- (1-conf)
  exitos <- sum(x)
  if (exitos == 0) {
    reslse <- 1 - (alfa/2)^(1/n)
  }
  else if (exitos == n) {
    reslse <- 1
  }
  else {
    reslse <- 1/(1+(n-exitos)/((exitos+1)*qf(1-alfa/2, 2*(exitos+1), 2*(n-exitos))))}
  return(reslse)
}

# Simulation results: Exact
lieres1 = list(); lieres2 = list(); lieres3 = list(); lieres4 = list(); lieres5 = list()
lseres1 = list(); lseres2 = list(); lseres3 = list(); lseres4 = list(); lseres5 = list()
exe1=list(); exe2=list(); exe3=list(); exe4=list(); exe5=list()
tote1=NULL; tote2=NULL; tote3=NULL; tote4=NULL; tote5=NULL

set.seed(1234)

for(i in 1:length(Ns)){
  ns <- Ns[i]
  dat1 <- rbinom(nrep*ns,1,p1)
  dat2 <- rbinom(nrep*ns,1,p2)
  dat3 <- rbinom(nrep*ns,1,p3)
  dat4 <- rbinom(nrep*ns,1,p4)
  dat5 <- rbinom(nrep*ns,1,p5)
  x1 <- matrix(dat1, ncol=nrep)
  x2 <- matrix(dat2, ncol=nrep)
  x3 <- matrix(dat3, ncol=nrep)
  x4 <- matrix(dat4, ncol=nrep)
  x5 <- matrix(dat5, ncol=nrep)
  lieres1[[i]] <- apply(x1, 2, lie)
  lseres1[[i]] <- apply(x1, 2, lse)
  lieres2[[i]] <- apply(x2, 2, lie)
  lseres2[[i]] <- apply(x2, 2, lse)
  lieres3[[i]] <- apply(x3, 2, lie)
  lseres3[[i]] <- apply(x3, 2, lse)
  lieres4[[i]] <- apply(x4, 2, lie)
  lseres4[[i]] <- apply(x4, 2, lse)
  lieres5[[i]] <- apply(x5, 2, lie)
  lseres5[[i]] <- apply(x5, 2, lse)
  exe1[[i]] <- lieres1[[i]] <= p1 & p1 <= lseres1[[i]]
  exe2[[i]] <- lieres2[[i]] <= p2 & p2 <= lseres2[[i]]
  exe3[[i]] <- lieres3[[i]] <= p3 & p3 <= lseres3[[i]]
  exe4[[i]] <- lieres4[[i]] <= p4 & p4 <= lseres4[[i]]
  exe5[[i]] <- lieres5[[i]] <= p5 & p5 <= lseres5[[i]]
  tote1[i] = sum(exe1[[i]]==T)/nrep
  tote2[i] = sum(exe2[[i]]==T)/nrep
  tote3[i] = sum(exe3[[i]]==T)/nrep
  tote4[i] = sum(exe4[[i]]==T)/nrep
  tote5[i] = sum(exe5[[i]]==T)/nrep
}

resultados.exact <- cbind(tote1, tote2, tote3, tote4, tote5); resultados.exact

rm(list = ls(pattern = "^lieres"))
rm(list = ls(pattern = "^lseres"))
rm(list = ls(pattern = "^exe"))
rm(dat1, dat2, dat3, dat4, dat5)
rm(x1, x2, x3, x4, x5)

################
#### T WALD ####
################

liwt <- function(x, conf=0.95) { 
  n <- length(x)
  pe <- mean(x)
  t <- abs(qt((1-conf)/2, df=n-1))
  return((pe)  - t*sqrt(pe*(1-pe)/n))
}

lswt <- function(x, conf=0.95) { 
  n <- length(x)
  pe <- mean(x)
  t <- abs(qt((1-conf)/2, df=n-1))
  return((pe)  + t*sqrt(pe*(1-pe)/n))
}

# Simulation results: T-Wald

liwtres1 = list(); liwtres2 = list(); liwtres3 = list(); liwtres4 = list(); liwtres5 = list()
lswtres1 = list(); lswtres2 = list(); lswtres3 = list(); lswtres4 = list(); lswtres5 = list()
exwt1=list(); exwt2=list(); exwt3=list(); exwt4=list(); exwt5=list()
totwt1=NULL; totwt2=NULL; totwt3=NULL; totwt4=NULL; totwt5=NULL

set.seed(1234)

for(i in 1:length(Ns)){
  ns <- Ns[i]
  dat1 <- rbinom(nrep*ns,1,p1)
  dat2 <- rbinom(nrep*ns,1,p2)
  dat3 <- rbinom(nrep*ns,1,p3)
  dat4 <- rbinom(nrep*ns,1,p4)
  dat5 <- rbinom(nrep*ns,1,p5)
  x1 <- matrix(dat1, ncol=nrep)
  x2 <- matrix(dat2, ncol=nrep)
  x3 <- matrix(dat3, ncol=nrep)
  x4 <- matrix(dat4, ncol=nrep)
  x5 <- matrix(dat5, ncol=nrep)
  liwtres1[[i]] <- apply(x1, 2, liwt)
  lswtres1[[i]] <- apply(x1, 2, lswt)
  liwtres2[[i]] <- apply(x2, 2, liwt)
  lswtres2[[i]] <- apply(x2, 2, lswt)
  liwtres3[[i]] <- apply(x3, 2, liwt)
  lswtres3[[i]] <- apply(x3, 2, lswt)
  liwtres4[[i]] <- apply(x4, 2, liwt)
  lswtres4[[i]] <- apply(x4, 2, lswt)
  liwtres5[[i]] <- apply(x5, 2, liwt)
  lswtres5[[i]] <- apply(x5, 2, lswt)
  exwt1[[i]] <- liwtres1[[i]] <= p1 & p1 <= lswtres1[[i]]
  exwt2[[i]] <- liwtres2[[i]] <= p2 & p2 <= lswtres2[[i]]
  exwt3[[i]] <- liwtres3[[i]] <= p3 & p3 <= lswtres3[[i]]
  exwt4[[i]] <- liwtres4[[i]] <= p4 & p4 <= lswtres4[[i]]
  exwt5[[i]] <- liwtres5[[i]] <= p5 & p5 <= lswtres5[[i]]
  totwt1[i] = sum(exwt1[[i]]==T)/nrep
  totwt2[i] = sum(exwt2[[i]]==T)/nrep
  totwt3[i] = sum(exwt3[[i]]==T)/nrep
  totwt4[i] = sum(exwt4[[i]]==T)/nrep
  totwt5[i] = sum(exwt5[[i]]==T)/nrep
}

results.waldt <- cbind(totwt1, totwt2, totwt3, totwt4, totwt5); results.waldt

rm(list = ls(pattern = "^liwtres"))
rm(list = ls(pattern = "^lswtres"))
rm(list = ls(pattern = "^exwt"))
rm(dat1, dat2, dat3, dat4, dat5)
rm(x1, x2, x3, x4, x5)

################
#### A WALD ####
################

#  Wald (Agresti-Coull)

liwa <- function(x, conf=0.95) { 
  n <- length(x)
  ex <- sum(x) 
  pt <- (ex+2)/(n+4)
  z = abs(qnorm((1-conf)/2))
  se <- sqrt(pt*(1-pt)/(n+4))
  resliwa = pt - z * se
  if(resliwa < 0) resliwa = 0
  return(resliwa)
}


lswa <- function(x, conf=0.95) { 
  n <- length(x)
  ex <- sum(x)
  pt <- (ex+2)/(n+4)
  z = abs(qnorm((1-conf)/2))
  se <- sqrt(pt*(1-pt)/(n+4))
  reslswa = pt + z * se
  if(reslswa > 1) reslswa = 1
  return(reslswa)
}

#Simulation results: Agresti-Coull 

liwares1 = list(); liwares2 = list(); liwares3 = list(); liwares4 = list(); liwares5 = list()
lswares1 = list(); lswares2 = list(); lswares3 = list(); lswares4 = list(); lswares5 = list()
exwa1=list(); exwa2=list(); exwa3=list(); exwa4=list(); exwa5=list()
totwa1=NULL; totwa2=NULL; totwa3=NULL; totwa4=NULL; totwa5=NULL

set.seed(1234)

for(i in 1:length(Ns)){
  ns <- Ns[i]
  dat1 <- rbinom(nrep*ns,1,p1)
  dat2 <- rbinom(nrep*ns,1,p2)
  dat3 <- rbinom(nrep*ns,1,p3)
  dat4 <- rbinom(nrep*ns,1,p4)
  dat5 <- rbinom(nrep*ns,1,p5)
  x1 <- matrix(dat1, ncol=nrep)
  x2 <- matrix(dat2, ncol=nrep)
  x3 <- matrix(dat3, ncol=nrep)
  x4 <- matrix(dat4, ncol=nrep)
  x5 <- matrix(dat5, ncol=nrep)
  liwares1[[i]] <- apply(x1, 2, liwa)
  lswares1[[i]] <- apply(x1, 2, lswa)
  liwares2[[i]] <- apply(x2, 2, liwa)
  lswares2[[i]] <- apply(x2, 2, lswa)
  liwares3[[i]] <- apply(x3, 2, liwa)
  lswares3[[i]] <- apply(x3, 2, lswa)
  liwares4[[i]] <- apply(x4, 2, liwa)
  lswares4[[i]] <- apply(x4, 2, lswa)
  liwares5[[i]] <- apply(x5, 2, liwa)
  lswares5[[i]] <- apply(x5, 2, lswa)
  exwa1[[i]] <- liwares1[[i]] <= p1 & p1 <= lswares1[[i]]
  exwa2[[i]] <- liwares2[[i]] <= p2 & p2 <= lswares2[[i]]
  exwa3[[i]] <- liwares3[[i]] <= p3 & p3 <= lswares3[[i]]
  exwa4[[i]] <- liwares4[[i]] <= p4 & p4 <= lswares4[[i]]
  exwa5[[i]] <- liwares5[[i]] <= p5 & p5 <= lswares5[[i]]
  totwa1[i] = sum(exwa1[[i]]==T)/nrep
  totwa2[i] = sum(exwa2[[i]]==T)/nrep
  totwa3[i] = sum(exwa3[[i]]==T)/nrep
  totwa4[i] = sum(exwa4[[i]]==T)/nrep
  totwa5[i] = sum(exwa5[[i]]==T)/nrep
}

results.wald_a <- cbind(totwa1, totwa2, totwa3, totwa4, totwa5); results.wald_a
rm(list = ls(pattern = "^liwares"))
rm(list = ls(pattern = "^lswares"))
rm(list = ls(pattern = "^exwa"))
rm(dat1, dat2, dat3, dat4, dat5)
rm(x1, x2, x3, x4, x5)


### Plots ###


# Wald
plot(Ns, totw1, type = "l", xlab = "Sample size", ylab = "Coverage", main="Wald, p=0.1", ylim=c(0.65, 1))
abline(h = 0.95, col = "firebrick")

plot(Ns, totw2, type = "l", xlab = "Sample sizel", ylab = "Coverage", main="Wald, p=0.2", ylim=c(0.65, 1))
abline(h = 0.95, col = "firebrick")

plot(Ns, totw3, type = "l", xlab = "Sample sizel", ylab = "Coverage", main="Wald, p=0.3", ylim=c(0.65, 1))
abline(h = 0.95, col = "firebrick")

plot(Ns, totw4, type = "l", xlab = "Sample sizel", ylab = "Coverage", main="Wald, p=0.4", ylim=c(0.65, 1))
abline(h = 0.95, col = "firebrick")

plot(Ns, totw5, type = "l", xlab = "Sample size", ylab = "Coverage", main="Wald, p=0.5", ylim=c(0.65, 1))
abline(h = 0.95, col = "firebrick")
