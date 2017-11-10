#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#        Kapitel 3: Stochastische Größen u. Verteilungen        #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# Bsp 3.2
#=========

# Abb 3.2
#---------
Fx <- function(x) {
   ifelse(x < 0 , 0, ifelse(x < 3/4, 2*x/3,
     ifelse(x < 3/2, 1/2, 1)))
}

x <- seq(-0.2, 1.7, by=0.005)
FF <- Fx(x)
FF[(1.49<x)&(x<1.51)] <- NA
old.par <- par(mar=c(5, 5, 4, 2) + 0.1)
plot(x, FF, type="l", lty=1, lwd=3, ylim=c(0,1),
  ylab=expression(F[X](x)), cex.lab=2)
lines(c(3/2,3/2), c(1/2,1), type="p", pch=c(21,19), cex=1.3)
abline(h=c(0,1), lty=2)
par(old.par)


# Bsp 3.3
#=========

# Abb 3.3
#---------
x <- seq(0, 1.5, by=0.001)
p <- seq(0, 1, by=0.001)
xp <- function(p) { min(x[Fx(x) >= p]) }
Q <- sapply(p, xp)
Q[(0.49<p)&(p<0.51)] <- NA
old.par <- par(mar=c(5, 5, 4, 2) + 0.1)
plot(p, Q, type="l", lty=1, lwd=3, xlab="p",
  ylab=expression(x[p]), cex.lab=2)
lines(c(0.5,0.5), c(3/4,1.5), type="p",
  pch=c(19,21), cex=1.3)
par(old.par)


# Bsp 3.4
#=========
x <- 2:12
px <- (6-abs(7-x))/36
Fx <- cumsum(px)

# Funktion zur Darstellung einer diskreten Verteilung
px.dis <- function(Merk, Px, lw=2, point=T, ...) {
  m <- length(Merk)
  Px <- Px /sum(Px)
  d.u <- diff(Merk)[1]
  d.o <- diff(Merk)[m-1]
  Merk.aug <- c(Merk[1]-d.u/3,Merk,Merk[m]+d.o/3)
  old.part <- par(mfrow=c(2,1), mar=c(4, 5, 2, 2) + 0.1)
  plot(Merk, Px, type="h", lty=1, lwd=lw, xlab="x",
    ylab=expression(p[X]), xlim=c(min(Merk.aug),max(Merk.aug)),
    cex.main=0.8, ...)
  if ( point ) {
    lines(Merk, Px, type="p", pch=19,
    cex=0.8, ...)
  }
  Fx.aug <- c(0,cumsum(Px),1)
  plot(Merk.aug, Fx.aug, type="s", lty=1, lwd=lw, xlab="x",
    ylab=expression(F[X]), xlim=c(min(Merk.aug),max(Merk.aug)),
    ylim=c(0,1), cex.main=0.8, ...)
  if ( point ) {
    lines(Merk, cumsum(Px), type="p", pch=19,
    cex=0.8, ...)
  }
  abline(h=c(0,1), lty=2)
  par(old.par)
}

# Abb 3.4
#---------
px.dis(2:12, px, lw=3, col="grey50")


# Bsp 3.6
#=========

# Abb 3.7
#---------
Fs <- function(x) {
  ifelse(x <= 0, 0, ifelse(x <= 3/4, 4*x/3, 1))
}
Fd <- function(x) {
  ifelse(x < 3/2, 0, 1)
}
Fx <- function(x) {
   ifelse(x < 0 , 0, ifelse(x < 3/4, 2*x/3,
     ifelse(x < 3/2, 1/2, 1)))
}

x <- seq(0, 2, by=0.001)
old.part <- par(mfrow=c(1,3), mar=c(4, 5, 2, 2) + 0.1, xpd=T)
plot(x, Fd(x), type="l", lty=1, lwd=3, ylab=expression(F[d]),
  cex.lab=1.8)
plot(x, Fs(x), type="l", lty=1, lwd=3, ylab=expression(F[s]),
  cex.lab=1.8)
plot(x, Fx(x), type="l", lty=1, lwd=3, ylab=expression(F),
  cex.lab=1.8)
par(old.par)


# Bsp 3.17
#==========
m <- 23^4
a <- 7200
c <- 1
g <- function(n, m, a, c, x0=1) {
   x <- numeric(n)
   x[1] <- (a*x0+1)%%m
   for (i in 2:n) {
     x[i] <- (a*x[i-1]+c)%%m
   }
  x
}

Nsim <- 10^4
x <- g(Nsim, m, a, c)
u <- x/m
hist(u, breaks=(0:10)/10, prob=T, main="",
  col="grey60")
lines(c(0,1), c(1,1), lty=2, lwd=2)
u1 <- x[-Nsim]
u2 <- u[-1]
plot(u1, u2)
acf(u)


# Bsp 3.18
#==========
lam <- 1
Nsim <- 10^4
U <- runif(Nsim)
X <- -log(1-U)/lam
hist(X, breaks=20, freq=F, main="", col="grey60",
  xlim=c(0,6))
x <- seq(0, 6, by=0.01)
lines(x, lam*exp(-lam*x), lty=1, lwd=2)


# Bsp 3.19
#==========
x <- 1:6
px <- x/sum(x)
Fx <- cumsum(px)
Nsim <- 10^4
u <- runif(Nsim)
xx <- sapply(u, function(u) x[min(which(Fx >= u))])
round(table(xx)/Nsim, 3)
round(px, 3)

px.xx <- rbind(px, table(xx)/Nsim)
rownames(px.xx) <- c("P(X=x)","rel. Häufig.")
barplot(px.xx, beside=T, axis.lty=1, xlab="x",
  col=c("grey50","grey80"), legend=rownames(px.xx),
  args.legend=list(x=5, y=0.25))


# UE - Aufgaben
#===============

# ue3.1
#=======
x <- 1:6
p.max <- (2*x-1)/36
Fmax <- cumsum(p.max)
p.min <- (13-2*x)/36
Fmin <- cumsum(p.min)
plot(c(0.5,x,6.5), c(0,Fmax,1), type="s", lty=1, lwd=3,
  xlab="x", ylab="F(x)", col="grey50",
  main=expression(paste("Verteilungsfunktion für ", X[min],
       " und ", X[max])))
lines(c(0.5,x,6.5), c(0,Fmin,1), type="s", lty=1, lwd=3,
  col="grey50")
abline(h=c(0,1), lty=2)
abline(v=c(1,6), lty=2)
text(locator(1), "Max", cex=1.5)
text(locator(1), "Min", cex=1.5)


# ue3.2
#=======
N <- 20; M <- 10
p <- M/(N+M)
x <- 0:15
Fx <- 1-(1-p)^x
plot(x, Fx, type="s", lty=1, lwd=3, xlab="x", ylab="F(x)",
  col="grey50", main=paste("Verteilungsfunktion (",
  "N = ", N, ", M = ", M, ")"))
abline(h=c(0,1), lty=2)


# ue3.4
#=======
x <- seq(-0.2, 3.2, by=0.001)
Fx <- function(x) { ifelse(x<0, 0, ifelse(x<=1, x^2/5,
  ifelse(x<=3, (-x^2+6*x-4)/5, 1))) }
fx <- function(x) { ifelse(x<0, 0, ifelse(x<=1, 2*x/5,
  ifelse(x<=3, (-2*x+6)/5, 0))) }
par(mfrow=c(2,1))
plot(x, Fx(x), type="l", lty=1, lwd=3, col="grey50", ylab="F(x)",
  main="Verteilungsfunktion")
xd <- c(0,1,3)
lines(xd, Fx(xd), type="p", pch=21, col="grey50", cex=1.5)
abline(h=c(0,1), lty=2)
plot(x, fx(x), type="l", lty=1, lwd=3, col="grey50", ylab="f(x)",
  main="Dichtefunktion")
xd <- c(0,1,3)
lines(xd, c(0,0.6,0), type="p", pch=21, col="grey50", cex=1.5)
abline(h=0, lty=2)
par(mfrow=c(1,1))


# ue3.5
#=======
Fx <- function(x, g, r) {
  ifelse(x<0, 0, ifelse(x>=r, 1, (x+g)/(g+r)))
}
x <- seq(-3, 70, by=0.05)
plot(x, Fx(x, g=25, r=65), type="l", lty=1, lwd=3, col="grey50",
  ylab="F(x)", main="Verteilungsfunktion")
abline(h=c(0,1), lty=2)
abline(v=c(0,65), lty=2)


# ue3.6
#=======
x <- seq(-5, 5, by=0.01)
Fx <- function(x) exp(x)/(1+exp(x))
fx <- function(x) exp(x)/(1+exp(x))^2
old.par <- par(mfrow=c(2,1))
plot(x, Fx(x), type="l", lty=1, lwd=3, col="grey50",
  ylab="F(x)", main="Verteilungsfunktion")
abline(h=c(0,1), lty=2)
plot(x, fx(x), type="l", lty=1, lwd=3, col="grey50",
  ylab="f(x)", ylim=c(0,0.25), main="Dichtefunktion")
abline(h=0, lty=2)
par(old.par)


# ue3.8
#=======
x <- seq(-5, 5, by=0.01)
fx <- function(x, d) 1/(pi*(1+(x/d)^2))
plot(x, fx(x, d=1), type="l", lty=1, lwd=3, col="grey50", xlab="y",
  ylab="f(y)", ylim=c(0,0.35), main="Cauchyverteilung")
abline(h=0, lty=2)


# ue3.9
#=======
x <- seq(0, 3, by=0.01)
fx <- function(x) ifelse(x<0, 0, 2*x*exp(-x^2))
plot(x, fx(x), type="l", lty=1, lwd=3, col="grey50",
  xlab="y", ylab=expression(f[Y](y)), main="Rayleighverteilung",
  axes=F)
axis(1, pos=0)
axis(2, pos=0)


# ue3.22
#========
x <- seq(0, 1, by=0.005)
fx <- function(x) { ifelse(x<0, 0, ifelse(x>1, 0, 3*x^2)) }
Fx <- function(x) { ifelse(x<0, 0, ifelse(x>1, 1, x^3)) }

par(mfrow=c(2,2))

plot(x, fx(x), type="l", lty=1, lwd=2, main="Dichte", ylab="f(x)")
abline(h=0, lty=2)
plot(x, Fx(x), type="l", lty=1, lwd=2, main="Verteilungsfunktion",
  ylab="F(x)")
abline(h=c(0,1), lty=2)
u <- runif(1)
xx <- u^(1/3)
arrows(0, u, xx, u, length=0.15, angle=15,
  lwd=2, col="grey60")
arrows(xx, u, xx, 0, length=0.15, angle=15,
  lwd=2, col="grey60")

Nsim <- 1000
brk <- seq(0, 1, by=0.1)
u <- runif(Nsim)
xx <- u^(1/3)
hist(xx, breaks=brk, prob=T, col="grey60", xlab="x",
  main=paste("Simulation","(N = ",Nsim,")"))
lines(x, fx(x), type="l", lty=2, lwd=2)

# GGZ
N <- 1:Nsim
xbar <- cumsum(xx)/N
plot(N, xbar, type="l", lty=1, lwd=2, col="grey50",
  xlab="n", ylab=expression(bar(x)[n]), log="x",
  main="Gesetz der großen Zahlen")
abline(h=3/4, lty=2)

par(mfrow=c(1,1))


# ue3.23
#========
x <- seq(-5, 5, by=0.01)
fx <- function(x) exp(x)/(1+exp(x))^2
Fx <- function(x) exp(x)/(1+exp(x))

par(mfrow=c(2,2))

plot(x, fx(x), type="l", lty=1, lwd=2, main="Dichte", ylab="f(x)",
  ylim=c(0,1/4))
abline(h=0, lty=2)
plot(x, Fx(x), type="l", lty=1, lwd=2, main="Verteilungsfunktion",
  ylab="F(x)")
abline(h=c(0,1), lty=2)
abline(v=0, lty=2)
u <- runif(1)
xx <- log(u/(1-u))
arrows(0, u, xx, u, length=0.15, angle=15,
  lwd=2, col="grey60")
arrows(xx, u, xx, 0, length=0.15, angle=15,
  lwd=2, col="grey60")

Nsim <- 10000
u <- runif(Nsim)
xx <- log(u/(1-u))
hist(xx, breaks=25, prob=T, col="grey60", xlab="x", xlim=c(-5,5),
  ylim=c(0,0.25), main=paste("Simulation","(N = ",Nsim,")"))
lines(x, fx(x), type="l", lty=2, lwd=2)

# GGZ
N <- 1:Nsim
xbar <- cumsum(xx)/N
plot(N, xbar, type="l", lty=1, lwd=2, col="grey50",
  xlab="n", ylab=expression(bar(x)[n]), log="x",
  main="Gesetz der großen Zahlen")
abline(h=0, lty=2)

par(mfrow=c(1,1))


# ue3.24
#========
Nsim <- 100
u <- runif(Nsim)
xx <- ifelse(u<=25/90, 0, 90*u-25)
round(xx)

# GGZ
N <- 1:Nsim
xbar <- cumsum(xx)/N
plot(N, xbar, type="l", lty=1, lwd=2, col="grey50",
  xlab="n", ylab=expression(bar(x)[n]), log="x",
  main="Gesetz der großen Zahlen")
abline(h=65^2/180, lty=2)


# ue3.25
#========
x <- seq(0, 1, by=0.005)
fx <- function(x) {30*(x^2-2*x^3+x^4)}
Fx <- function(x, a=0) {30*(x^3/3-x^4/2+x^5/5) - a}

par(mfrow=c(2,2))

plot(x, fx(x), type="l", lty=1, lwd=2, main="Dichte")
abline(h=0, lty=2)
plot(x, Fx(x), type="l", lty=1, lwd=2, main="Verteilungsfunktion")
abline(h=c(0,1), lty=2)
u <- runif(1)
xx <- uniroot(Fx, interval=c(0,1), a=u)$root
arrows(0, u, xx, u, length=0.15, angle=15,
  lwd=2, col="grey60")
arrows(xx, u, xx, 0, length=0.15, angle=15,
  lwd=2, col="grey60")

Nsim <- 10000
brk <- seq(0, 1, by=0.1)
u <- runif(Nsim)
xx <- sapply(u, function(u) {uniroot(Fx, interval=c(0,1), a=u)$root})
hist(xx, breaks=brk, prob=T, col="grey60", xlab="x",
  main=paste("Simulation","(N = ",Nsim,")"))
lines(x, fx(x), type="l", lty=2, lwd=2)

# GGZ
N <- 1:Nsim
xbar <- cumsum(xx)/N
plot(N, xbar, type="l", lty=1, lwd=2, col="grey50",
  xlab="n", ylab=expression(bar(x)[n]), log="x",
  main="Gesetz der großen Zahlen")
abline(h=0.5, lty=2)

par(mfrow=c(1,1))
