#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#          Kapitel 6: Folgen von stochastischen Größen          #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# Bsp 6.7
#=========
# Abb 6.2
#---------
x <- 0:4
px <- c(0.1,0.2,0.3,0.35,0.05)
(EX <- sum(x*px))
(EX2 <- sum(x^2*px))
(VX <- EX2 - EX^2)

Nsim <- 10000
xx <- sample(x, Nsim, replace=TRUE, prob=px)
N <- 1:Nsim
xbar <- cumsum(xx)/N
par(mar=c(5, 5, 4, 2) + 0.1)
plot(N, xbar, type="l", lty=1, lwd=2, col="grey50",
  xlab="n", ylab=expression(bar(X)[n]), log="x")
abline(h=EX, lty=2)
text(locator(1), "E(X)", cex=1.5)
par(mar=c(5, 4, 4, 2) + 0.1)


# Bsp 6.9
#=========
# Abb 6.3
#---------
conv <- function(px, n) {
  a <- b <- px
  if (n > 1) {
    for (k in 2:n) {
      b <- convolve(b, a, type="open")
    }
  }
  return(b)
}

x <- 0:4
px <- c(0.25,0.15,0.1,0.2,0.3)
n <- c(1,5,10,25)

old.par <- par(mfrow=c(2,2))
for (i in 1:length(n)) {
  xn <- 0:(n[i]*4)
  pn <- conv(px, n[i])
  plot(xn, pn, type="l", lty=1, ylim=c(0,max(pn)),
    xlab="x", ylab="p(x)", , cex.axis=0.8)
  lines(xn, pn, type="p", pch=21, col=2)
  mtext(line=0.5, paste("n = ", n[i]))
}
par(old.par)


# Bsp 6.10
#==========
# Abb 6.4
#---------
old.par <- par(mfrow=c(2,2))
for (n in c(3,5,10,50)) {
  results <- c()
  mu <- sigma <- 1
  for (i in 1:10000) {
    X <-  rexp(n, 1/mu)
    results[i] <- (mean(X)-mu)/(sigma/sqrt(n))
  }
  his <- hist(results, breaks=seq(min(results)-1, max(results)+1,
           by=0.25), plot=FALSE)
  ylim <- range(his$density, dnorm(0))
  hist(results, breaks=seq(min(results)-1, max(results)+1, by=0.25),
    prob=TRUE, xlab="y", xlim=c(-4,4), ylim=ylim,
    main="", col="grey80", cex.axis=0.8)
  x <- seq(-4, 4, by=0.1)
  lines(x, dnorm(x), lty=1, lwd=2, col="red")
  mtext(line=0.5, paste("n = ", n))
}
par(old.par)


# UE - Aufgaben
#===============

# ue6.3
#=======
x <- 1:4
px <- rep(1,4)/4
conv <- function(px, n) {
  a <- b <- px
  if (n > 1) {
    for (k in 2:n) {
      b <- convolve(b, a, type="open")
    }
  }
  return(b)
}

conv(px, 2) *16
conv(px, 3) *64


# ue6.4
#=======
fxpy <- function(a) {
  ifelse(a<0, 0,
  ifelse(a<1, a/2,
  ifelse(a<2, 1/2,
  ifelse(a<3, (3-a)/2, 0))))
}

N <- 10000
x <- runif(N, min=0, max=1)
y <- runif(N, min=0, max=2)
xpy <- x+y
a <- seq(-0.1, 3.1, length=300)
hist(xpy, breaks=30, freq=F, col="grey70", axes=F,
  xlab="a", ylab=expression(f[X+Y](a)), main="")
axis(1, pos=0)
axis(2, pos=-0.1)
lines(a, fxpy(a), type="l", lty=1, lwd=3, col=2)


# ue6.5
#=======
fX <- function(x) {
  ifelse(x<0, 0, exp(-x/4)-exp(-x/3))
}
N <- 10000
xa <- rexp(N, rate=1/4)
xr <- rexp(N, rate=1/3)
x <- xa + xr
mean(x)
sd(x)

y <- seq(0, 30, length=500)
hist(x, breaks=30, freq=F, col="grey70", axes=F, xlab="x",
  ylab=expression(f[X](x)), xlim=c(0,30), main="")
axis(1, pos=0)
axis(2, pos=0)
lines(y, fX(y), type="l", lty=1, lwd=3, col=2)


# ue6.10
#========
n <- 10^6
x <- runif(n, -1, 1)
y <- runif(n, -1, 1)
estim.pi <- 4*sum(x^2+y^2<1)/n
matrix(c(estim.pi, sqrt(pi*(4-pi)/n)), ncol=2,
  dimnames=list(c(""), c("estim","sd")))
  
# Zusatz: Grafische Darstellung
#-------------------------------
# Note: Adapted from J. Verzani
simpi3 <- function(n=1000) {
  plot(0, 0, pch="", xlim=c(-1,1), ylim=c(-1,1), xlab="x",
    ylab="y", main=substitute(n == list(x), list(x=n)),
    cex=2, axes=FALSE)
  axis(1, pos=-1)
  axis(2, pos=-1)
  axis(3, pos=1)
  axis(4, pos=1)
  theta <- seq(0, 2*pi, length=100)
  polygon(cos(theta), sin(theta))
  x <- runif(n, min=-1, max=1)
  y <- runif(n, min=-1, max=1)
  inorout <- x^2+y^2 < 1
  points(x, y, pch=21, col=as.numeric(inorout)+1, cex=0.5)
  text(0, 0, substitute(hat(pi)==list(x),
    list(x=4*sum(inorout)/n)), cex=2)
  return(as.numeric(inorout))
}

par(pty="s")
n <- 10000
simpi <- simpi3(n)
simpi <- 4*cumsum(simpi)/(1:n)
plot(1:n, simpi, type="l", lty=1)
lines(1:n, pi+2*sqrt(4*pi*(1-pi/4)/(1:n)), lty=2)
lines(1:n, pi-2*sqrt(4*pi*(1-pi/4)/(1:n)), lty=2)
abline(h=pi)


# ue6.14
#========
n <- 1
f <- function(n) pnorm((2000-100*n)/(30*sqrt(n)))
while ( f(n) > 0.05 ) n <- n+1
n