#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#               Kapitel 4: Spezielle Verteilungen               #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# Funktion zur grafischen Darstellung 
#     von diskreten Verteilungen
#=====================================
px.dis <- function(Merk, Px, lw=2, point=T, ...) {
  m <- length(Merk)
  d.u <- diff(Merk)[1]
  d.o <- diff(Merk)[m-1]
  Merk.aug <- c(Merk[1]-d.u/2,Merk,Merk[m]+d.o/2)
  old.par <- par(mfrow=c(2,1), mar=c(4, 5, 2, 2) + 0.1)
  plot(Merk, Px, type="h", lty=1, lwd=lw, xlab="x", 
    ylab=expression(p[X]), xlim=c(min(Merk.aug),max(Merk.aug)), 
    ylim=c(0,max(Px)*(1+0.2)), cex.main=0.8, ...)
  abline(h=0, lty=2)
  if ( point ) {
    lines(Merk, Px, type="p", pch=19, 
    cex=0.8, ...) 
  }
  Fx.aug <- c(0,cumsum(Px),cumsum(Px)[m])
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


# Bsp 4.1: Diskrete uniforme Verteilung
#=======================================
x <- 0:10
px <- rep(1, length(x))/length(x)
Fx <- cumsum(px)
px.dis(x, px, lw=3, col="grey50")

# Zufallszahlen
M <- 0:10
sample(M, size=25, replace=TRUE)


# Bsp 4.2: Bernoulli-Verteilung
#===============================
x <- 0:1
px <- c(0.3,0.7)
Fx <- cumsum(px)
px.dis(x, px, lw=3, col="grey50")


# Bsp 4.3: Binomialverteilung
#=============================
x <- 0:10
px <- dbinom(x, 10, 0.7)
Fx <- cumsum(px)
px.dis(x, px, lw=3, col="grey50")


# Bsp 4.4: Geometrische Verteilung
#==================================
x <- 0:20
px <- dgeom(x, prob=1/10)
Fx <- cumsum(px)
px.dis(x, px, lw=3, col="grey50")


# Bsp 4.5: Hypergeometrische Verteilung
#=======================================
N <- 100
n <- 22
v <- 2
p <- (0:N)/N
hyp <- sapply(p, function(x) phyper(v, N*x, N*(1-x), n))
bin <- sapply(p, function(x) pbinom(v, n, x))
hyp.bin <- cbind(hyp, bin)
matplot(p, hyp.bin, type="o", lty=c(1,2), pch=19, lwd=2, cex=0.7, 
  col="grey50", ylab="P(Los akzeptiert)", xlim=c(0,0.3))
legend("topright", legend=c("Hypergeometrisch","Binomial"), 
  lty=c(1,2), pch=19, lwd=2, col="grey50")


# Bsp 4.6: Poisson-Verteilung
#=============================
binom.pois <- function(n, p) {
 int <- 0:n
 binom <- dbinom(int, n, p)
 pois <- dpois(int, n*p)
 int1 <- ( (binom >0.001)|(pois >0.001) )
 binom.pois <- data.frame(Binomial=binom[int1], Poisson=pois[int1])
 rownames(binom.pois) <- int[int1]
 binom.pois.m <- t(as.matrix(binom.pois))
 binom.pois.t <- as.table(binom.pois.m)
 barplot(binom.pois.t, beside=TRUE, axis.lty=1, col=c("grey50","grey80"), 
   xlab="x", ylab="P(X=x)", legend=colnames(binom.pois), 
   ylim=c(0,max(binom.pois))) 
 invisible(binom.pois)
}

binom.pois(50, 1/10)


# Exponentialverteilung (Abb 4.8)
#=================================
old.par <- par(mfrow=c(2,1), las=1, mar=c(4, 5, 2, 3)+0.1)
lamb <- 1/c(0.5,1.0,2)
x <- seq(0, 3.0, length=300)
matplot(x, outer(x, lamb, pexp), lty=1:3, type="l", lwd=3,
  xlim=c(0,3), ylim = c(0,1), ylab = "F(x)", xlab="x",
  main="(a) Verteilungsfunktion", axes=FALSE, cex.main=0.75,
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3))
axis(2, pos=0, at=c(0,0.5,1))
legend(x=2.2, y=0.5, expression(tau == 0.5,
                                tau == 1.0,
                                tau == 2.0),
                                lty=1:3, col="grey50", lwd=3)

matplot(x, outer(x, lamb, dexp), lty=1:3, type="l", lwd=3,
  xlim=c(0,3), ylim=c(0,2), ylab="f(x)", xlab="x",
  main="(b) Dichte", axes=FALSE, cex.main=0.75, 
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3))
axis(2, pos=0, at=c(0,0.5,1,1.5,2))
legend(x=2.2, y=1.2, expression(tau == 0.5,
                              tau == 1.0,
                              tau == 2.0),
                              lty=1:3, col="grey50", lwd=3)
par(old.par)


# Normalverteilung (Abb 4.10)
#=============================
old.par <- par(mfrow=c(2,1), las=1, mar=c(4, 5, 2, 3)+0.1)
sig <- 1
mu <- -2:2
x <- seq(-6, 6, length=500)
matplot(x, outer(x, mu, pnorm, sd=sig), lty=1:5, type="l", lwd=3,
  xlim=c(-6,6), ylim = c(0,1), ylab = "", xlab="x",
  main="(a) Verteilungsfunktion", axes=FALSE, cex.main=0.75,
  col="grey50")
axis(1, pos=0, at=-6:6)
axis(2, pos=0, at=c(0,0.5,1), labels=c("",0.5,1))
legend(x=3.5, y=0.8, expression(mu == -2,
                                mu == -1,
                                mu == -0,
                                mu ==  1,
                                mu ==  2),
                                lty=1:5, col="grey50", lwd=3)

matplot(x, outer(x, mu, dnorm, sd=sig), lty=1:5, type="l", lwd=3,
  xlim=c(-6,6), ylim=c(0,0.4), ylab="", xlab="x",
  main="(b) Dichte", axes=FALSE, cex.main=0.75, 
  col="grey50")
axis(1, pos=0, at=-6:6)
axis(2, pos=0, at=c(0,0.2,0.45), labels=c("",0.2,""))
legend(x=3.5, y=0.4, expression(mu == -2,
                                mu == -1,
                                mu == -0,
                                mu ==  1,
                                mu ==  2),
                                lty=1:5, col="grey50", lwd=3)
par(old.par)


# UE - Aufgaben
#===============

# ue4.2
#=======
# (a)
#-----
n <- 15
p <- 0.2
prob1 <- dbinom(0:n, n, p)
plot(0:n, prob1, type="h", lty=1, lwd=5, col="grey50", 
  pch=19, axes=FALSE, xlab="x", ylab="Binomialwahrscheinlichkeiten")
axis(1, pos=0)
axis(2, pos=0)

# (b)
#-----
n <- 15
p <- (1:9)/10
prob2 <- matrix(nrow=n+1, ncol=9)
for (i in 1:9) {
  prob2[,i] <- dbinom(0:n, n, p[i])
}
matplot(0:n, prob2, type="o", lty=1, lwd=2, col="grey50", 
  pch=19, axes=FALSE, xlab="x", ylab="Binomialwahrscheinlichkeiten")
axis(1, pos=0)
axis(2, pos=0)

# (c)
#-----
p <- 0.05
n <- c(10,20,50,200)
prob3 <- matrix(nrow=201, ncol=4)
for (i in 1:4) {
  prob3[,i] <- dbinom(0:200, n[i], p)
}
matplot(0:20, prob3[1:21,], type="o", lty=1, lwd=2, col="grey50", 
  pch=19, axes=FALSE, xlab="x", ylab="Binomialwahrscheinlichkeiten")
axis(1, pos=0)
axis(2, pos=0)


# ue4.4
#=======
polyroot(c(-3,12,-15,6))

# Zusatz: Abbildung
#-------------------
k <- 1:6
K <- 2*k-1
p <- seq(0, 1, by=0.01)
f <- function(x, y) {
  sum(dbinom(y:(2*y-1), 2*y-1, x))
}
vf <- Vectorize(f)
intakt <- outer(p, k, vf)
matplot(p, intakt, type="l", lty=1:length(k), 
  lwd=2, col="grey50", ylab="P(System intakt)",
  axes=FALSE)
axis(1, pos=0)
axis(2, pos=0)
axis(3, pos=1)
axis(4, pos=1)
abline(h=0.5, lty=2)
abline(v=0.5, lty=2)
legend(locator(1), paste(K, "K", sep=""),
  lty=1:length(k), lwd=2, col="grey50")


# ue4.5
#=======
k <- 0:10
(pa <- sum(dbinom(k, 10, 0.8)*pbinom(k-1, 10, 0.85)))
(pb <- sum(dbinom(k, 10, 0.85)*pbinom(k-1, 10, 0.8)))
(pc <- sum(dbinom(k, 10, 0.8)*dbinom(k, 10, 0.85)))


# ue4.7
#=======
n <- 5; N <- 15; A <- 6
round(dhyper(0:n, A, N-A, n), 4)
p <- A/N
(E <- n*p)
(V <- (N-n)/(N-1)*n*p*(1-p))
1-phyper(3, A, N-A, n)

# Zusatz1
Personen <- c(paste("M",1:6, sep=""),paste("F",1:9, sep=""))
sort(sample(Personen, 5))

# Zusatz2
1-phyper(3, A, N-A, n)


# ue4.13
#========
N <- 500
n <- 50
d <- 1
# Hypergeometrisch
#A <- N*c(0.008,0.09)
A <- 0:50 # für Abbildung
ocH <- phyper(d, A, N-A, n)
abH <- 1-ocH
# Binomial
p <- A/N
ocB <- pbinom(d, n, p)
abB <- 1-ocB
# Poisson
mu <- n*A/N
ocP <- ppois(d, mu)
abP <- 1-ocP

ab <- cbind(abH, abB, abP)
oc <- cbind(ocH, ocB, ocP)

matplot(p, oc, type="l", lty=1:3, lwd=2, col=1:3,
  xlab="p", ylab="P(Losannahme)")
legend("topright", c("Hypergeometrisch","Binomial","Poisson"),
  lty=1:3, lwd=2, col=1:3)
  
  
# ue4.20
#========
old.par <- par(mfrow=c(2,1), las=1, mar=c(4, 5, 2, 3)+0.1)
y <- seq(0, 3.0, length=300)
plot(y, 2*pnorm(y)-1, lty=1, type="l", lwd=3,
  xlim=c(0,3), ylim=c(0,1), ylab=expression(F[Y](y)), xlab="y",
  main="(a) Verteilungsfunktion", axes=FALSE, cex.main=0.75,
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3))
axis(2, pos=0, at=c(0,0.5,1))

plot(y, 2*dnorm(y), lty=1, type="l", lwd=3,
  xlim=c(0,3), ylab=expression(f[Y](y)), xlab="y",
  main="(b) Dichte", axes=FALSE, cex.main=0.75, 
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3))
axis(2, pos=0)
par(old.par)


# ue4.22
#========
old.par <- par(mfrow=c(2,1), las=1, mar=c(4, 5, 2, 3)+0.1)
y <- seq(0, 4.0, length=300)
plot(y, plnorm(y), lty=1, type="l", lwd=3,
  xlim=c(0,4), ylim=c(0,1), ylab=expression(F[Y](y)), xlab="y",
  main="(a) Verteilungsfunktion", axes=FALSE, cex.main=0.75,
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3,4))
axis(2, pos=0, at=c(0,0.5,1))

plot(y, dlnorm(y), lty=1, type="l", lwd=3,
  xlim=c(0,4), ylab=expression(f[Y](y)), xlab="y",
  main="(b) Dichte", axes=FALSE, cex.main=0.75, 
  col="grey50")
axis(1, pos=0, at=c(0,1,2,3,4))
axis(2, pos=0)
par(old.par)


# ue4.23
#========
n <- 5
(c <- qchisq(0.025, n))
(d <- qchisq(0.975, n))
# check
pchisq(d, n) - pchisq(c, n)


# ue4.24
#========
m <- 5; n <- 10
(a <- qf(0.05, m, n))
(b <- qf(0.95, m, n))
# check
pf(b, m, n) - pf(a, m, n)


# ue4.25
#========
n <- 14
(b <- qt(0.95, n))
# check
pt(b, n) - pt(-b, n)