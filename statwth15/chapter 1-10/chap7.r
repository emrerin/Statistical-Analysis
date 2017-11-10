#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#               Kapitel 7: Schliessende Statistik               #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------
# r01.txt
# sodium.txt
# cholesterol.txt
# euroweight.txt
# lifetimes.txt

# Packages:
#-----------
# base, MASS


# Bsp 7.3
#=========
dist.exp <- function(DATA, tau=NA, plotit=TRUE) {
 n <- length(DATA)
 mu <- ifelse(is.na(tau),mean(DATA),tau)
 d.s <- sort(DATA)
 e.1 <- (1:n)/n
 e.2 <- (0:(n-1))/n
 Fx <- pexp(d.s,1/mu)
 ABS <- c(abs(Fx-e.1),abs(Fx-e.2))
 MAX <- max(ABS)
 ind <- which.max(ABS)
 if ( plotit ) {
   plot(ecdf(DATA),verticals=TRUE,do.p=FALSE,
      main="Emp. VF - Exponential",xlab="x",
      ylab=expression(hat(F)[n](x)))
   t <- seq(0,d.s[n]+max(diff(d.s))/5,length=100)
   lines(t,pexp(t,1/mu),lty=1,col="red",lwd=2)
   IND <- ifelse(ind <= n,ind,ind-n)
   x.max <- sort(DATA)[IND]
   abline(v=x.max,lty=2)
   points(c(x.max,x.max),c(pexp(x.max,1/mu),pexp(x.max,1/mu) +
        MAX*ifelse(ind<=n,1,-1)),pch=21,cex=1.1)
   mtext(paste("n =",n,"     D =",round(MAX,4),"     x =",
        round(d.s[IND],4)),side=3,line=0,cex=0.8)
   }
   return(list(data.frame(x=d.s,dist.1=abs(Fx-e.1),
     dist.2=abs(Fx-e.2)),maxDist=MAX))
}

dist.exp(x <- rexp(10, rate=1/2), tau=2)
dist.exp(x <- rexp(100, rate=1/2), tau=2)


# Bsp 7.8
#=========
require(MASS)
x <- c(11.96,5.03,67.40,16.07,31.50,7.73,11.10,22.38)
# avoid spurious accuracy
op <- options(digits=3)
fitdistr(x, "gamma")
# now do this with more control
fitdistr(x, dgamma, start=list(shape=1, rate=0.1),
  lower=0.001)
options(op)


# Bsp 7.12
#==========
n <- 100
N <- 10000
eVF <- numeric(N)
quant <- 0.5
qx <- qexp(quant, rate=1)
for (k in 1:N) {
  x <- rexp(n, rate=1)
  eVF[k] <- sum( x <= qx )/n
}

par(mar=c(5, 5, 4, 2) + 0.1)
hist(eVF, breaks=30, prob=T, col="grey70", axes=F,,
  xlab=expression(hat(F)[n](tilde(x))), xlim=c(0.2,0.8),
  main="")
axis(1, pos=0)
axis(2, pos=0.2)
y <- seq(0, 1, length=300)
mx <- pexp(qx); sx <- sqrt(pexp(qx)*(1-pexp(qx))/n)
lines(y, dnorm(y, mean=mx, sd=sx), type="l", lty=1,
  lwd=3, col=2)
par(mar=c(5, 4, 4, 2) + 0.1)


# Bsp 7.13
#==========
# Abb 7.5
#---------
mu <- 0; sx <- 1
alph <- 0.05
n <- 10
N <- 100
ki <- c()
for (i in 1:N) {
  x <- rnorm(n, mean=mu, sd=sx)
  ki <- rbind(ki, mean(x)+c(-1,1)*qnorm(1-alph/2)/sqrt(n))
}
f <- function(x) {
  if ((mu<x[1])|(mu>x[2])) y <- 2
  else y <- 1
  y
}
ind <- apply(ki, 1, f)
plot(c(ki[1,1],ki[1,2]), c(1,1), type="l", lwd=2, col=ind[1],
  xlim=c(min(ki),max(ki)), ylim=c(-3,N), xlab="", ylab="",
  axes=FALSE)
for (i in 2:N) {
  lines(c(ki[i,1],ki[i,2]), c(i,i), lty=1, lwd=2, col=ind[i])
}
lines(c(mu,mu), c(-0.5,N+1), lty=1, lwd=3)
text(locator(1), expression(mu), cex=1.5)


# Bsp 7.14
#==========
n <- 10^6
x <- rexp(n, rate=1)
y <- sqrt(x)
alph <- 0.05
op <- options(digits=5)
(Ihat <- mean(y))
(Ihat + c(-1,1)*qnorm(1-alph/2)*sd(y)/sqrt(n))
options(op)


# Bsp 7.15
#==========
x <- c(81,87,86,82,90,86,96,73,74,75,72,80,66,72,56,82)
y <- c(78,91,78,78,84,67,92,70,58,62,70,58,66,60,65,73)
d <- x-y

boxplot(d, vertical=TRUE, main="Boxplot",
  ylab="Difference", col="grey50")
lines(rep(1, length(d)), d, type="p", pch=19)

t.test(x, y, paired=TRUE)
t.test(x, y, var.equal=TRUE)


# Bsp 7.16
#==========
# Abb 7.7
#---------
# Standard (Wald)
cp1 <- function(p, n, alph) {
  y <- 0
  for (x in 0:n) {
    T1 <- x/n - qnorm(1-alph/2)*sqrt((x/n)*(1-x/n)/n)
    T2 <- x/n + qnorm(1-alph/2)*sqrt((x/n)*(1-x/n)/n)
    if ((T1 < p)&(p < T2)) {
      y <- y + dbinom(x, n, p)
    }
  }
  return(y)
}

# Score
cp2 <- function(p, n, alph) {
  y <- 0
  a <- qnorm(1-alph/2)
  for (x in 0:n) {
    T1 <- (a^2+2*n*(x/n)-a*sqrt(a^2-4*n*(x/n-1)*x/n))/(2*(a^2+n))
    T2 <- (a^2+2*n*(x/n)+a*sqrt(a^2-4*n*(x/n-1)*x/n))/(2*(a^2+n))
    if ((T1 < p)&(p < T2)) {
      y <- y + dbinom(x, n, p)
    }
  }
  return(y)
}

alph <- 0.05
n <- 35
p <- seq(0.001, 0.999, by=0.001)
uew1 <- sapply(p, cp1, n=n, alph=alph)
uew2 <- sapply(p, cp2, n=n, alph=alph)
uew <- cbind(uew1,uew2)
matplot(p, uew, type="l", lty=1, lwd=2, col=c("grey50",1),
  ylim=c(0.8,1), ylab="Coverage Probability")
abline(h=1-alph, lty=2, lwd=2)
legend("topright", c("standard","score"), lty=1, lwd=2,
  col=c("grey50",1))


# Bsp 7.17
#==========
# Abb 7.8
#---------
# Standard (Wald)
cp1 <- function(lamb, n, alph) {
  y <- 0
  a <- qnorm(1-alph/2)
  for (x in 0:100) {
    T1 <- x/n - a*sqrt((x/n)/n)
    T2 <- x/n + a*sqrt((x/n)/n)
    if ((T1 < lamb)&(lamb < T2)) {
      y <- y + dpois(x, n*lamb)
    }
  }
  return(y)
}

# Score
cp2 <- function(lamb, n, alph) {
  y <- 0
  a <- qnorm(1-alph/2)
  for (x in 0:100) {
    T1 <- (a^2+2*n*(x/n)-a*sqrt(a^2+4*n*x/n))/(2*n)
    T2 <- (a^2+2*n*(x/n)+a*sqrt(a^2+4*n*x/n))/(2*n)
    if ((T1 < lamb)&(lamb < T2)) {
      y <- y + dpois(x, n*lamb)
    }
  }
  return(y)
}

alph <- 0.05
n <- 35
lamb <- seq(0.001, 2, by=0.001)
uew1 <- sapply(lamb, cp1, n=n, alph=alph)
uew2 <- sapply(lamb, cp2, n=n, alph=alph)
uew <- cbind(uew1,uew2)
matplot(lamb, uew, type="l", lty=1, lwd=2, col=c("grey50",1),
  ylim=c(0.8,1), xlab=expression(lambda),
  ylab="Coverage Probability")
abline(h=1-alph, lty=2, lwd=2)
legend("topright", c("standard","score"), lty=1, lwd=2,
  col=c("grey50",1))
  
  
# Bsp 7.18
#==========
percentciboot <- function(x, B, alpha) {
  # x: vector containing the original sample
  # B: desired number of resamples
  # alpha: (1 - alpha) is the confidence coefficient
  #
  # theta: point estimate.
  # lower: lower end of the percentile confidence interval
  # upper: upper end of the percentile confidence interval
  # thetastar: vector of bootstrapped theta^*s.
  #
  theta <- mean(x)
  thetastar <- rep(0, B)
  n <- length(x)
  for (i in 1:B) {
    xstar <- sample(x, n, replace=TRUE)
    thetastar[i] <- mean(xstar)
  }
  thetastar <- sort(thetastar)
  m <- floor((alpha/2)*B)
  lower <- thetastar[m]
  upper <- thetastar[B-m+1]
  list(theta=theta, lower=lower, upper=upper, thetastar=thetastar)
}

# Exponential
tau <- 3
n <- 25
x <- rexp(n, rate=1/tau)
tauhat <- mean(x)

op <- options(digits=4)

# KIe
alpha <- 0.05
exact <- 2*n*tauhat/qchisq(c(1-alpha/2, alpha/2), 2*n)
wald <- tauhat + c(-1,1)*qnorm(1-alpha/2)*tauhat/sqrt(n)
score <- tauhat/(1 + c(1,-1)*qnorm(1-alpha/2)/sqrt(n))
erg <- percentciboot(x, 3000, alpha)
ci.boot <- c(erg$lower, erg$upper)

# Histogramm
y <- erg$thetastar
hist(erg$thetastar, breaks=20, prob=TRUE,
 xlab=expression(bar(x)^" *"), main="", col="gray50")
abline(v=c(erg$lower,erg$theta,erg$upper), lty=2, lwd=2)

# Tabelle
ci.all <- rbind(exact, wald, score, ci.boot)
len <- ci.all[,2] - ci.all[,1]
ci <- data.frame(ci.all[,1], ci.all[,2], len)
names(ci) <- c(paste(c("Lower","Upper"), 100*c(alpha/2,1-alpha/2),
  "%", sep=""),"Length")
row.names(ci) <- c("Exact", "Wald", "Score", "Boot")
ci

options(op)


# Bsp 7.23
#==========
x <- c(52.1,49.0,51.4,50.0,50.3,49.6,50.6,50.8,51.0,51.7)
t.test(x, mu=50)


# Bsp 7.25
#==========
p0 <- 0.02
n <- 300
alpha <- 0.05
g <- function(k, n, p) { sum(dbinom(k:n, n, p)) }
k <- 0
while (g(k, n, p0) > alpha) {
  k <- k+1
}
(kstar <- k)
# p-value
(pv <- g(10, n, p0))


# Bsp 7.26/7.27
#===============
p0 <- 0.04
n <- 500
p1 <- sum(dbinom(0:16, n, p0))
p2 <- sum(dbinom(16:500, n, p0))
(pv <- 2*min(p1, p2))

# binom.test()
binom.test(16, 500, p=0.04)

k <- 0:500
ind <- which(dbinom(k, n, p0) <= dbinom(16, n, p0))
sum(dbinom(ind-1, n, p0))

# Normalapproximation
z0 <- (16-n*p0)/(sqrt(n*p0*(1-p0)))
2*(1-pnorm(abs(z0)))


# Bsp 7.28
#==========
x <- c(3,7,25,10,15,6,12,25,15,7)
y <- c(48,44,40,38,33,21,20,12,1,18)
mean(x); mean(y)
sd(x); sd(y)
t.test(x, y, var.equal=TRUE)
t.test(x, y, var.equal=FALSE)


# Bsp 7.29
#==========
x <- c(3,7,25,10,15,6,12,25,15,7)
y <- c(48,44,40,38,33,21,20,12,1,18)
var.test(x, y)
var.test(x, y, alternative="less")


# Bsp 7.31
#==========
x <- c(176,183,185,190,191,192,201,205,214,220)
n <- length(x)
x <- sort(x)
p <- (1:n-0.5)/n
z <- qnorm((1:n-0.5)/n)
nnetz <- data.frame(x=x, p=p, z=z)
round(nnetz, 3)

net.normal1 <- function(dat, ml.line=TRUE, h.line=FALSE, txt=TRUE) {
  require(MASS)
  MLE <- fitdistr(dat, "normal")$estimate
  mu <- MLE[1]
  sig <- MLE[2]

  ld <- length(dat)

  # pp = (i-0.5)/n
  i <- 1:ld
  pp <- (i-0.5)/ld
  plot(sort(dat), qnorm(pp), pch=19, cex=1.2,
    xlab=expression(x[(i)]), ylab=expression(z[i]))

  if (ml.line) {
    abline(-mu/sig, 1/sig, lwd=2, col=2)
  }

  if (h.line) {
    abline(h=c(-1,0,1), lty=2, lwd=1)
    abline(v=sig*c(-1,0,1)+mu, lty=2, lwd=1)
  }

  grid()

  if (txt) {
    mtext("Normal QQ-Plot", line=2, cex=1.3, font=1)
    mtext(paste("MLE:   mu =",round(mu,2),"  sig =",
      round(sig,2)), line=0.5, cex=1, font=1)
  }

  list(mu, sig)
}

net.normal1(x, h.line=TRUE)


# Bsp 7.32
#==========
H <- c(13,19,11,8,5,4)
sum(H)
w <- rep(1/6,6)
(H-60*w)^2/(60*w)
sum((H-60*w)^2/(60*w))

chisq.test(H, p=w)


# Bsp 7.33
#==========
chi2.normal <- function(DAT, nk) {
  n <- length(DAT)
  if ( nk < 4 ) return( cat("Zu wenig Klassen!", "\n") )
  if ( nk > n/5 ) return( cat("Zu viele Klassen!", "\n") )
  m <- mean(DAT)
  s <- sd(DAT)*sqrt((n-1)/n)
  cl <- cut(DAT, breaks=c(-Inf,m+s*qnorm((1:(nk-1))/nk),Inf))
  result <- chisq.test(table(cl), p=rep(1/nk,nk))
  obs <- as.data.frame(result$observed)
  rownames(obs) <- 1:nk
  names(obs)[1] <- "class"
  tab <- data.frame(
      obs,
      expec = result$expected,
      xsquared = (as.vector(result$residuals))^2)
  return( list(
      tab,
      result$statistic,
      df = nk-3,
      pvalue = as.vector(pchisq(result$statistic,
        df=nk-3, lower.tail=FALSE)) ))
}

dat <- scan("rn01.txt")
m <- mean(dat)
s2 <- var(dat)*(length(dat)-1)/length(dat)
class.2 <- cut(dat, breaks=c(-Inf,m+sqrt(s2)*qnorm((1:7)/8),Inf))
result <- chisq.test(table(class.2), p=rep(1/8,8))
result

# Zus. Chi^2-Test
chi2.normal(dat, 8)

# QQ-Plot
net.normal1(dat)


# UE - Aufgaben
#===============

# ue7.1
#=======
dist.norm <- function(DATA, mu=NA, sig=NA, plotit=TRUE) {
 n <- length(DATA)
 mu <- ifelse(is.na(mu),mean(DATA),mu)
 sig <- ifelse(is.na(sig),sd(DATA),sig)
 d.s <- sort(DATA)
 e.1 <- (1:n)/n
 e.2 <- (0:(n-1))/n
 Fx <- pnorm(d.s,mu,sig)
 ABS <- c(abs(Fx-e.1),abs(Fx-e.2))
 MAX <- max(ABS)
 ind <- which.max(ABS)
 if ( plotit ) {
   plot(ecdf(DATA),verticals=TRUE,do.p=FALSE,
      main="Emp. VF - Normal",xlab="DATA",
      ylab=expression(Fx[n](x)))
   t <- seq(d.s[1]-max(diff(d.s))/10,d.s[n]+max(diff(d.s))/10,
     length=100)
   lines(t,pnorm(t,mu,sig),lty=1,col="red",lwd=2)
   IND <- ifelse(ind <= n,ind,ind-n)
   x.max <- sort(DATA)[IND]
   abline(v=x.max,lty=2)
   points(c(x.max,x.max),c(pnorm(x.max,mu,sig),pnorm(x.max,mu,sig) +
        MAX*ifelse(ind<=n,1,-1)),pch=21,cex=1.1)
   mtext(paste("N =",n,"     D =",round(MAX,4),"     x =",
        round(d.s[IND],4)),side=3,line=0,cex=0.8) }
   return(list(data.frame(x=d.s,dist.1=abs(Fx-e.1),
     dist.2=abs(Fx-e.2)),maxDist=MAX))
}

dist.norm(rnorm(10), mu=0, sig=1)
dist.norm(rnorm(100), mu=0, sig=1)


# ue7.8
#=======
n <- 10^6
x <- runif(n)
y <- 4*sqrt(1-x^2)
alph <- 0.05
opt <- options(digits=5)
(Ihat <- mean(y))
(Ihat + c(-1,1)*qnorm(1-alph/2)*sd(y)/sqrt(n))
options(opt)


# ue7.15
#========
n <- 300
x <- 13
alph <- 0.05
a <- qnorm(1-alph/2)
wald <- x/n + c(-1,1)*a*sqrt((x/n)*(1-x/n)/n)
score <- (a^2+2*n*(x/n) + c(-1,1)*a*sqrt(a^2-4*n*(x/n-1)*x/n))/(2*(a^2+n))

op <- options(digits=3)
(phat <- x/n)
ci.all <- rbind(wald, score)
len <- ci.all[,2] - ci.all[,1]
ci <- data.frame(ci.all[,1], ci.all[,2], len)
names(ci) <- c(paste(c("Lower","Upper"),
  100*c(alph/2,1-alph/2), "%", sep=""),"Length")
row.names(ci) <- c("Standard", "Score")
ci
options(op)


# ue7.17
#========
# s. Bsp 7.18
x <- c(79,88,39,17,40,27,45,100,50,71)
erg <- percentciboot(x, 3000, alpha=0.05)
(ci.boot <- c(erg$lower, erg$upper))
# Histogram
y <- erg$thetastar
hist(erg$thetastar, breaks=20, prob=TRUE,
  xlab=expression(bar(x)^" *"), main="", col="gray50")
abline(v=c(erg$lower,erg$theta,erg$upper), lty=2, lwd=2)


# ue7.22
#========
x <- scan("sodium.txt")
t.test(x, mu=130)


# ue7.25
#========
dat <- read.table("cholesterol.txt", header=TRUE)
t.test(dat$vor, dat$nach, paired=TRUE, alternative="greater")


# ue7.26
#========
x <- c(10.33,9.53,9.82,10.11,8.99,10.37,9.99,9.01)
y <- c(9.75,8.44,10.03,10.67,9.30,10.68,11.14,7.87)
cor.test(x, y, method="pearson", alternative="greater")


# ue7.28
#========
x <- scan("sodium.txt")
net.normal1(x)


# ue7.30
#========
# s. Bsp 7.31
dat <- read.table("euroweight.txt", header=TRUE)
attach(dat)
old.par <- par(mfrow=c(4,2), pin=c(2.1,1.3))
for (i in 1:8) {
 subdat <- weight[batch==i]
 net.normal1(subdat, txt=FALSE)
 mtext(paste("Batch ",i), line=1, cex=0.7, font=1)
 }
par(old.par)
detach(dat)


# ue7.33
#========
dat <- scan("rn01.txt")
class.1 <- cut(dat, breaks=c(-Inf,qnorm((1:7)/8),Inf))
result <- chisq.test(table(class.1), p=rep(1/8,8))
result


# ue 7.35
#=========
chi2.exp <- function(DAT, nk) {
  n <- length(DAT)
  if ( nk < 3 ) return( cat("Zu wenig Klassen!", "\n") )
  if ( nk > n/5 ) return( cat("Zu viele Klassen!", "\n") )
  m <- mean(DAT)
  cl <- cut(DAT, breaks=c(0,qexp((1:(nk-1))/nk, 1/m),Inf))
  result <- chisq.test(table(cl), p=rep(1/nk,nk))
  obs <- as.data.frame(result$observed)
  rownames(obs) <- 1:nk
  names(obs)[1] <- "class"
  tab <- data.frame(
      obs,
      expec = result$expected,
      xsquared = (as.vector(result$residuals))^2)
  return( list(
      tab,
      result$statistic,
      df = nk-2,
      pvalue = as.vector(pchisq(result$statistic,
        df=nk-2, lower.tail=FALSE)) ))
}

x <- scan("lifetimes.txt")
chi2.exp(x, 4)
