#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#              Kapitel 5: Multivariate Verteilungen             #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base, MASS


# Bsp 5.3
#=========
# Abb 5.1 (gemeinsame Dichte)
#-----------------------------
x1 <- x2 <- seq(0, 1, length=30)
f <- function(x1, x2) 6*x1^2*x2
z <- outer(x1, x2, f)
persp(x1, x2, z, theta=30, phi=30, expand=0.5, col="lightblue",
  ticktype="detailed", zlab="f(x1,x2)")


# Bsp 5.4
#=========
# Abb 5.2 (Randdichten)
#-----------------------
old.par <- par(mar=c(5, 5, 4, 2) + 0.1)
x1 <- x2 <- seq(0, 1, length=100)
plot(x1, 3*x1^2, type="l", lty=1, lwd=5, col="grey50",
  xlab=expression(x[1]), ylab=expression(f[1](x[1])),
  pty="s", cex.lab=1.5)
abline(h=0, lty=2)
mtext(text=expression(paste("Randdichte von  ", X[1])),
  line=0.5, cex=1.5)
plot(x2, 2*x2, type="l", lty=1, lwd=5, col="grey50",
  xlab=expression(x[2]), ylab=expression(f[2](x[2])),
  pty="s", cex.lab=1.5)
abline(h=0, lty=2)
mtext(text=expression(paste("Randdichte von  ", X[2])),
  line=0.5, cex=1.5)
par(old.par)


# Bsp 5.11
#==========
# Abb 5.4 (gemeinsame Dichte)
#-----------------------------
x1 <- x2 <- seq(0, 1, length=30)
f <- function(x1, x2) x1 + x2
z <- outer(x1, x2, f)
persp(x1, x2, z, theta=30, phi=30, expand=0.5, col="lightblue",
  ticktype="detailed", zlab="f(x1,x2)")
  
  
# Bsp 5.18
#==========
# Abb 5.6 (k-aus-n System)
#--------------------------
f <- function(x, lamb, k, n) {
  p <- 1-pexp(x, rate=lamb)
  pbinom(k-1, size=n, prob=p)
}

n <- 10
lamb <- 1/100
x <- seq(0, 500, length=500)
Fkn <- c()
for (k in 1:n) {
  z <- sapply(x, f, lamb=lamb, k=k, n=n)
  Fkn <- cbind(Fkn, z)
}

matplot(x, Fkn, type="l", lty=1:n, lwd=3, col="grey50")
legend("bottomright", paste(1:n, "-aus-", n, " ", sep=""),
  lty=1:n, lwd=3, col="grey50")
  
  
# Bsp 5.21
#==========
# Abb 5.7 (bivariate Standardnormaldichte)
#------------------------------------------
f <- function(x, y) {
  z <- (1/(2*pi)) * exp(-.5 * (x^2 + y^2))
}

y <- x <- seq(-3, 3, length=50)
z <- outer(x, y, f)
persp(x, y, z, theta=45, phi=30, expand=0.6, ltheta=120,
  shade=0.75, ticktype="detailed", col="lightblue",
  xlab="x", ylab="y", zlab="f(x, y)")
  
# Abb 5.8 (Contourplots)
#------------------------
norm.2d <- function(x, y, m.x=0, m.y=0, s.x=1, s.y=1, r.xy=0) {
     C.1 <- 1/(2*pi*s.x*s.y*sqrt(1-r.xy^2))
     C.2 <- 1/(2*(1-r.xy^2))
     e <- (((x-m.x)/s.x)^2 -
         2*r.xy*((x-m.x)/s.x)*((y-m.y)/s.y) +
         ((y-m.y)/s.y)^2)
     return(C.1*exp(-C.2*e))
     }

x <- y <- seq(-3, 3, by=0.05)

old.par <- par(mfrow=c(2,2), pin=c(2.1,2.1))

s.x <- 1
s.y <- 2
r.xy <- 0
z <- outer(x, y, norm.2d, s.x=s.x, s.y=s.y, r.xy=r.xy)
contour(x, y, z, nlevels=8, xlab="x", ylab="y", col="darkblue", lwd=2)
mtext(line=0.5, substitute(group("(",list(sigma[1],sigma[2],rho),")") ==
  group("(",list(a,b,c),")"), list(a=s.x, b=s.y, c=r.xy)))

s.x <- 1
s.y <- 2
r.xy <- 0.5
z <- outer(x, y, norm.2d, s.x=s.x, s.y=s.y, r.xy=r.xy)
contour(x, y, z, nlevels=8, xlab="x", ylab="y", col="darkblue", lwd=2)
mtext(line=0.5, substitute(group("(",list(sigma[1],sigma[2],rho),")") ==
  group("(",list(a,b,c),")"), list(a=s.x, b=s.y, c=r.xy)))

s.x <- 1
s.y <- 2
r.xy <- 0.8
z <- outer(x, y, norm.2d, s.x=s.x, s.y=s.y, r.xy=r.xy)
contour(x, y, z, nlevels=8, xlab="x", ylab="y", col="darkblue", lwd=2)
mtext(line=0.5, substitute(group("(",list(sigma[1],sigma[2],rho),")") ==
  group("(",list(a,b,c),")"), list(a=s.x, b=s.y, c=r.xy)))

s.x <- 1
s.y <- 1
r.xy <- -0.5
z <- outer(x, y, norm.2d, s.x=s.x, s.y=s.y, r.xy=r.xy)
contour(x, y, z, nlevels=8, xlab="x", ylab="y", col="darkblue", lwd=2)
mtext(line=0.5, substitute(group("(",list(sigma[1],sigma[2],rho),")") ==
  group("(",list(a,b,c),")"), list(a=s.x, b=s.y, c=r.xy)))

par(old.par)


# UE - Aufgaben
#===============

# ue5.5
#=======
# (a)
#-----
fxy <- function(x, y) 6/7*(x^2+x*y/2)
x <- seq(0, 1, by=0.05)
y <- seq(0, 2, by=0.1)
zf <- outer(x, y, fxy)
persp(x, y, zf, theta=30, phi=30, col="lightblue", expand=0.8,
  ltheta=120, shade=0.3, ticktype="detailed", xlab="x", ylab="y",
  zlab="f(x,y)")
  
  
# ue5.12
#========
# (a)
#-----
fxy <- function(x, y) 2*exp(-(x+2*y))
x <- seq(0, 2, by=0.05)
y <- seq(0, 1, by=0.05)
zf <- outer(x, y, fxy)
persp(x, y, zf, theta=30, phi=30, col="lightblue", expand=0.8,
  ltheta=120, shade=0.3, ticktype="detailed", xlab="x", ylab="y",
  zlab="f(x,y)")
  
  
# ue5.15
#========
# (b) (Simulation)
#------------------
n <- 100000
S <- numeric(n)
for (i in 1:n) {
  K1 <- max(rexp(3, rate=1/1000))
  K2 <- max(rexp(2, rate=1/3000))
  K3 <- rexp(1, rate=1/5000)
  S[i] <- min(K1,K2,K3)
}

hist(S, breaks=50, prob=TRUE, main="Systemlebensdauer",
  xlab="Stunden", col="lightgrey", xlim=c(0,6000))
abline(v=c(mean(S),median(S)), lty=c(1,2), lwd=2)
legend("topright", legend=c("Mittelwert","Median"),
  lty=c(1,2), lwd=2)


# ue5.16
#========
a <- c(1,-2,1)
Sig <- matrix(c(3,2,1,2,2,1,1,1,3), ncol=3, byrow=TRUE)
t(a) %*% Sig %*% a


# ue5.17
#========
# inch
mu.in <- c(70,64)
sigx.in <- 2
sigy.in <- 1.5
rho <- 0.7
Sig.in <- matrix(c(sigx.in^2,rep(sigx.in*sigy.in*rho,2),sigy.in^2), ncol=2)
# cm
A <- matrix(c(2.54,rep(0,2),2.54), ncol=2)
(mu.cm <- mu.in*2.54)
(Sig.cm <- A %*% Sig.in %*% t(A))

# Prognose
(male.cm <- 6*12*2.54)
(progn <- mu.cm[2]+0.7*sqrt(Sig.cm[2,2]/Sig.cm[1,1])*(male.cm-mu.cm[1]))

# Regressionsschere
x <- seq(165, 190, by=0.1)
y <- seq(155, 170, by=0.1)
m.x <- mu.cm[1]
m.y <- mu.cm[2]
s.x <- sqrt(Sig.cm[1,1])
s.y <- sqrt(Sig.cm[2,2])
r.xy <- rho
a1 <- mu.cm[2]-r.xy*(s.y/s.x)*mu.cm[1]
a2 <- mu.cm[2]-(r.xy)^(-1)*(s.y/s.x)*mu.cm[1]
z <- outer(x, y, norm.2d, m.x=m.x, m.y=m.y, s.x=s.x, s.y=s.y, r.xy=r.xy)
contour(x, y, z, nlevels=10, xlab="x", ylab="y", col="grey80", lwd=2)
lines(m.x, m.y, type="p", pch=21, cex=1.3, col="darkblue")
abline(a=a1, b=r.xy*s.y/s.x, lty=1, lwd=3, col="darkblue")
text(locator(1), "E(Y|x)", cex=1.8, col="darkblue", pos=4)
abline(a=a2, b=(r.xy)^(-1)*s.y/s.x, lty=1, lwd=3, col="darkblue")
text(locator(1), "E(X|y)", cex=1.8, col="darkblue", pos=2)
abline(v=male.cm, lty=2)
abline(h=progn, lty=2)


# ue5.19
#========
require(MASS)
biv.rnorm <- function(n=1, mux=0, muy=0, sigx=1, sigy=1, rho=0) {
  mu <- c(mux,muy)
  Sigma <- matrix(c(sigx^2,rep(rho*sigx*sigy,2),sigy^2), ncol=2)
  rn <- mvrnorm(n, mu, Sigma)
  return(rn)
}

# Note: The following functions were adapted
# from J. Verzani (Package: UsingR)
scatter.with.hist <- function(x, y, hist.col="grey50", trend.line="lm", ...) {
  on.par <- par(no.readonly=TRUE)
  on.exit(par(on.par))
  nf <- layout(matrix(c(1,0,
                        3,2),
                        2,2, byrow=TRUE),
                      widths=c(3,1),
                      heights=c(1,3),
                      respect=TRUE)
  layout.show(nf)
  n <- length(x)
  no.breaks <- max(nclass.scott(x),nclass.scott(y))
  xhist <- hist(x, breaks=no.breaks, plot=FALSE)
  yhist <- hist(y, breaks=no.breaks, plot=FALSE)
  top <- max(c(xhist$counts, yhist$counts))
  par(mar=c(0,4,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0,top),
    space=0, col=hist.col)
  #box()
  par(mar=c(4,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0,top),
    space=0, col=hist.col, horiz=TRUE)
  #box()
  par(mar=c(4,4,1,1))
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  plot(x, y, pch=21, xlab=x.name, ylab=y.name, ...)
  if(!is.null(trend.line) && !is.na(trend.line)) {
    switch(trend.line,
      "lm"=abline(lm(y ~ x)),
      "supsmu"=lines(supsmu(x, y)),
      "lowess"=lines(lowess(x, y)),
      NULL)
    }
  invisible()
}

scatter.with.box <- function(x, y, box.col="grey50", ...) {
  on.par <- par(no.readonly=TRUE)
  on.exit(par(on.par))
  nf <- layout(matrix(c(1,0,
                        3,2),
                        2,2, byrow=TRUE),
                      widths=c(3,1),
                      heights=c(1,3),
                      respect=TRUE)
  layout.show(nf)
  n <- length(x)
  par(mar=c(0,4,1,1))
  boxplot(x, col=box.col, axes=FALSE, horizontal=TRUE)
  par(mar=c(4,0,1,1))
  boxplot(y, col=box.col, axes=FALSE)
  par(mar=c(4,4,1,1))
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  plot(x, y, pch=21, xlab=x.name, ylab=y.name, ...)
  invisible()
}

# (a)
#-----
BNa <- biv.rnorm(500)
x <- BNa[,1]
y <- BNa[,2]
scatter.with.hist(x, y, trend.line=FALSE)
scatter.with.box(x, y)

# (b)
#-----
BNb <- biv.rnorm(500, 100, 200, 5, 6, 0.8)
x <- BNb[,1]
y <- BNb[,2]
scatter.with.hist(x, y, trend.line=FALSE)
scatter.with.box(x, y)

# (c)
#-----
BNc <- biv.rnorm(500, 100, 200, 5, 6, -0.6)
x <- BNc[,1]
y <- BNc[,2]
scatter.with.hist(x, y, trend.line=FALSE)
scatter.with.box(x, y)


# ue5.20
#========
require(MASS)
Sig <- matrix(c(3,2,1,2,2,1,1,1,3), ncol=3, byrow=TRUE)
rtrinorm <- mvrnorm(500, mu=c(0,0,0), Sigma=Sig)
colnames(rtrinorm) <- c("X1","X2","X3")
panel.reg <- function(x, y) {
  points(x, y, pch=21, cex=1.5, col="lightblue")
  abline(lm(y ~ x), col="darkblue", lwd=2)
}
pairs(rtrinorm, panel=panel.reg, gap=0)


# Anhang: Methode der Annahme und Verwerfung
#============================================

# Bsp1
#------
f <- function(x, a, b) { dbeta(x, a, b) }
g <- function(x) { dunif(x) }

a <- 2; b <- 4
x <- seq(0, 1, by=0.01)
plot(x, f(x, a, b), type="l", lty=1, lwd=3,
  xlab="x", ylab="", ylim=c(0,2.5), axes=FALSE)
axis(1, pos=0)
axis(2, pos=0, at=seq(0, 2.5, by=0.5))

cc <- 135/64

N <- 1000
rndf <- numeric(N)
k <- 0
while( k < N ) {
  Y <- runif(1)
  U <- runif(1, min=0, max=cc*g(Y))
  if ( U <= f(Y, a, b) ) {
    rndf[k+1] <- Y
    k <- k + 1
  }
}

hist(rndf, breaks=seq(0, 1, by=0.1), prob=TRUE,
  col="grey70", main="", xlab="x", ylim=c(0, 2.5), 
  axes=FALSE)
axis(1, pos=0)
axis(2, pos=0, at=seq(0, 2.5, by=0.5))
lines(x, f(x, a, b), lty=1, lwd=3)


# Bsp2
#------
f <- function(x, a, b) { 2*dnorm(x) }
g <- function(x) { dexp(x) }

x <- seq(0, 4, by=0.01)
plot(x, f(x), type="l", lty=1, lwd=3,
  xlab="x", ylab="", ylim=c(0,1), axes=FALSE)
axis(1, pos=0)
axis(2, pos=0, at=seq(0, 1, by=0.2))

cc <- sqrt(2*exp(1)/pi)

N <- 1000
rndf <- numeric(N)
k <- 0
while( k < N ) {
  Y <- rexp(1)
  U <- runif(1, min=0, max=cc*g(Y))
  if ( U <= f(Y) ) {
    rndf[k+1] <- Y
    k <- k + 1
  }
}

hist(rndf, breaks=seq(0, 4, by=0.2), prob=TRUE,
  col="grey70", main="", xlab="x", ylim=c(0, 1), 
  axes=FALSE)
axis(1, pos=0)
axis(2, pos=0, at=seq(0, 1, by=0.2))
lines(x, f(x), lty=1, lwd=3)
