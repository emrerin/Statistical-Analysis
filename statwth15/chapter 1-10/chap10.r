#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#                    Kapitel 10: Ergänzungen                    #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# 10.1 Überraschung und Entropie
#================================


# 10.3 Bayes'sche Netze
#=======================

# Bsp 10.5 (Klassifikation)
#===========================
x <- c(1.67, 2.00, 4.23)
x <- matrix(rep(x, 4), nrow=4, byrow=T)
mu <- matrix(
  c(1.6, 2.4, 4.3,
    1.5, 2.9, 6.1,
    1.8, 2.5, 4.2,
    1.1, 3.1, 5.6),
  ncol=3, byrow=T)

sig <- matrix(
  c(0.1, 0.5, 0.2,
    0.2, 0.6, 0.9,
    0.3, 0.3, 0.3,
    0.2, 0.7, 0.3),
  ncol=3, byrow=T)
  
NN <- exp((-0.5)*apply((x-mu)^2/sig^2, 1, sum))
DD <- apply(sig, 1, prod)
NN/DD


# 10.3 Markow-Ketten
#====================

# Bsp 10.6 (Wettermodell)
#=========================
P <- matrix(c(0.6,0.4,0.3,0.7), nrow=2, byrow=TRUE)
(P2 <- P %*% P)
(P3 <- P2 %*% P)

p0 <- c(0.8,0.2)
t(p0) %*% P
t(p0) %*% P2
t(p0) %*% P3


# Bsp 10.7 (Ehrenfestmodell)
#============================
efmatrix <- function(M) {
  P <- matrix(rep(0, (M+1)^2), nrow=M+1, ncol=M+1)
  for (i in 1:M) { 
    P[i,i+1] <- (M-(i-1))/M 
  }
  for (i in 2:(M+1)) { 
    P[i,i-1] <- (i-1)/M
  }
  return(P)
}

(P <- efmatrix(4))


# Bsp 10.8 (Random Walk, 1-dim)
#===============================
n <- 100
N <- 5
RW <- matrix(nrow=n, ncol=N)
for (i in 1:N) {
  rb <- 2*rbinom(n, size=1, prob=1/2) - 1
  rw <- cumsum(rb)
  RW[,i] <- rw  
}
RW0 <- rbind(rep(0, N), RW)
matplot(c(0:n), RW0, type="l", lty=1,
  xlab="Zeitpunkt", ylab="Zustand",
  main="Random Walk (symmetrisch)")
abline(h=0, lty=2, lwd=2)
lines(0, 0, type="p", pch=19, cex=1)


# Abb 10.8 (Random Walk, 2-dim)
#-------------------------------
N <- 10000
xy <- matrix(nrow=N+1, ncol=2)
Hx <- 0; Hy <- 0
H <- c(Hx,Hy)
xy[1,] <- H
for (i in 2:(N+1)) {
  st <- sample(1:4, 1)
  Hx <- Hx + ifelse(st==1, -1, ifelse(st==2, 1, 0))
  Hy <- Hy + ifelse(st==3, 1, ifelse(st==4, -1, 0))
  xy[i,] <- c(Hx,Hy)
}

plot(xy[,1], xy[,2], type="l", lty=1, lwd=0.5,
  xlab="x", ylab="y")
lines(xy[1,1], xy[1,2], type="p", pch=19, cex=0.7)
lines(xy[N+1,1], xy[N+1,2], type="p", pch=19, cex=0.7)


# Bsp 10.11
#===========
MC.sim <- function(n, P, x1) {
  sim <- as.numeric(n)
  m <- ncol(P)
  if (missing(x1)) {
    sim[1] <- sample(1:m, 1) # random start
  } else { sim[1] <- x1 }
  for (i in 2:n) {
    newstate <- sample(1:m, 1, prob=P[sim[i-1],])
    sim[i] <- newstate
  }
 sim
}

P <- efmatrix(4)

N <- 10000
SIM <- MC.sim(N, P)-1
head(SIM, n=100)
table(SIM)/N
dbinom(0:4, 4, 1/2)


# 10.5 Zählprozesse
#===================

# Bsp 10.12
#===========
tt <- c(0:7,7)
x <- c(0,177,242,293,336,368,395,410,420)
plot(x, tt, type="s", lty=1, lwd=2, xlab="Kalenderzeit t", 
  ylab="Anzahl der Ausfälle X(t)")
  
  
# Bsp 10.13 (Bernoulli-Prozess)
#===============================
Delta <- 0.05
p <- 0.1
tt <- seq(Delta, 10, by=Delta)
n <- length(tt)
N <- 3
X <- matrix(nrow=n, ncol=N)
for (i in 1:N) {
  X[,i] <- cumsum(rbinom(n, 1, p))
}
X <- rbind(rep(0,N), X)

matplot(c(0,tt), X, type="s", lty=1, lwd=2,
  xlab="Zeit (Minuten)", ylab="X(n)")
x <- seq(0, 10, by=0.05)
lines(x, x/Delta*p, lty=2, lwd=2)
lines(x, x/Delta*p + 2*sqrt(x/Delta*p*(1-p)), lty=1, lwd=2)
lines(x, x/Delta*p - 2*sqrt(x/Delta*p*(1-p)), lty=1, lwd=2)


# Bsp 10.14 (Poisson-Prozess)
#=============================
pprocess <- function(lamb, T0) {
  x <- 0
  while (sum(x) < T0) {
    ex <- rexp(1, rate=lamb)
    x <- c(x, ex)
  }
  x
}

lamb <- 2 #/min
T0 <- 30
x <- pprocess(lamb, T0)
xx <- cumsum(x)
yy <- 0:(length(x)-1)
plot(xx, yy, type="s", lty=1, lwd=2,
  xlab="Zeit (Minuten)", ylab="X(t)",
  xlim=c(0,T0), ylim=c(0,60))

for (k in 2:5) {
  x <- pprocess(lamb, T0)
  xx <- cumsum(x)
  yy <- 0:(length(x)-1)
  lines(xx, yy, type="s", lty=1, lwd=2, col=k)
}
  
z <- seq(0, T0, by=0.01)
lines(z, z*lamb, lty=2, lwd=2)
lines(z, z*lamb + 2*sqrt(z*lamb), lty=1, lwd=2)
lines(z, z*lamb - 2*sqrt(z*lamb), lty=1, lwd=2)


# Bsp 10.16
#===========
lamb <- 2 # /min
T0 <- 30
x <- 0
while (sum(x) < T0) {
 ex <- rexp(1, rate=lamb)
 x <- c(x, ex)
}

round(x, 2)
xx <- cumsum(x)
round(xx, 2)


# 10.5 Warteschlangen
#=====================

# Bsp 10.18 (M/M/1)
#===================
lambdaA <- 5         # arrival rate
lambdaS <- 6         # service rate
duration <- 10       # total duration of the simulation
t <- 0               # current time in the simulation
queue <- 0           # start with empty queue
s <- 0               # running sum for computing average
                     # number of jobs in the system
N <- Inf             # capacity of the system (set N=Inf 
                     # if no limitation)

rho <- lambdaA/lambdaS
a <- 1:N
(L <- (1-rho)/(1-rho^(N+1))*sum(a * rho^a))

# first arrival to start process
T1 <- rexp(1, rate=lambdaA)
currentqueue <- 1
eventsTime <- T1
t <- T1
nEvents <- 1  # total number of events that have occurred

while (t < duration) {
  nEvents <- nEvents + 1
  if (currentqueue > 0) { # nonempty queue
     T1 <- rexp(1, rate=lambdaA + lambdaS) # time until next event
  # is event an arrival or departure?
  p <- runif(1)
  queue[nEvents] <- currentqueue # how many are in the system 
                                 # before this new event?
  currentqueue <- ifelse(p < lambdaA/(lambdaA + lambdaS),
                         min(currentqueue + 1, N), # arrival
                         currentqueue - 1)         # departure
  } 
  else { # empty system
    T1 <- rexp(1, rate=lambdaA)
    queue[nEvents] <- currentqueue
    currentqueue <- 1
  }
  t <- t + T1                # time of next arrival
  eventsTime[nEvents] <- T1  # inter-event time
  s <- s + T1*queue[nEvents] # spent T1 at nth number 
                             # of jobs in the system
}

# Schätzwert für E(X)
s/t

lq <- length(queue)
tt <- c(0,cumsum(eventsTime)[-lq])
bServ <- ifelse(queue > 0, 1, 0)
nWait <- queue - bServ
mq <- max(queue)

pp <- NA
cl <- "grey60"
old.par <- par(mfrow=c(3,1), mar=c(5, 5, 4, 3))
plot(tt, queue, type="s", main="Jobs in the system",
  xlab="Time (min)", ylab="# (jobs in the system)", 
  xlim=c(0,duration), ylim=c(0,mq), col=cl, axes=F)
lines(c(tt[lq],duration), rep(queue[lq], 2), col=cl)
axis(1, pos=pp); axis(2, pos=pp)
plot(tt, nWait, type="s", main="Jobs waiting", 
  xlab="Time (min)", ylab="# (jobs waiting)", 
  xlim=c(0,duration), ylim=c(0,mq), col=cl, axes=F)
lines(c(tt[lq],duration), rep(nWait[lq], 2), col=cl)
axis(1, pos=pp); axis(2, pos=pp)
plot(tt, bServ, type="s", main="Busy periods",
  xlab="Time (min)", ylab="# (busy servers)", 
  xlim=c(0,duration), col=cl, axes=F)
lines(c(tt[lq],duration), rep(bServ[lq], 2), col=cl)
axis(1, pos=pp); axis(2, at=c(0,1), pos=pp)
par(old.par)


# Bsp 10.19 (M/M/1/N)
#=====================
# Setze N <- 3 im Code zu Bsp 10.18


# Bsp 10.20 (M/M/k)
#===================
lambdaA <- 10        # arrival rate
lambdaS <- 6         # service rate
duration <- 200      # total duration of the simulation
t <- 0               # current time in the simulation
queue <- 0           # start with empty queue
s <- 0               # running sum for computing average
                     # number of jobs in the system
k <- 3               # number of servers

# first arrival to start process
T1 <- rexp(1, rate=lambdaA)
currentqueue <- 1
eventsTime <- T1
t <- T1
nEvents <- 1  # total number of events that have occurred
bServ <- 1    # number of busy servers (<= k)

while (t < duration) {
  nEvents <- nEvents + 1
  if (currentqueue > 0) { # nonempty queue
     T1 <- rexp(1, rate=lambdaA + bServ*lambdaS) # time until next event
  # is event an arrival or departure?
  p <- runif(1)
  queue[nEvents] <- currentqueue # how many are in the system 
                                 # before this new event?
  
  currentqueue <- ifelse(p < lambdaA/(lambdaA + bServ*lambdaS),
                         currentqueue + 1, # arrival
                         currentqueue - 1) # departure
  bServ <- ifelse(p < lambdaA/(lambdaA + bServ*lambdaS),
                  min(bServ + 1, k), # arrival
                  ifelse(queue[nEvents] <= k, bServ-1, bServ)) # departure
  } 
  else { # empty system
    T1 <- rexp(1, rate=lambdaA)
    queue[nEvents] <- currentqueue
    currentqueue <- 1
    bServ <- 1
  }
  
  t <- t + T1                # time of next arrival
  eventsTime[nEvents] <- T1  # inter-event time
  s <- s + T1*queue[nEvents] # spent T1 at nth number 
                             # of jobs in the system
}

# Schätzwert für E(X)
s/t

lq <- length(queue)
tt <- c(0,cumsum(eventsTime)[-lq])
bServ <- ifelse(queue > 0, 1, 0)
nWait <- queue - bServ
mq <- max(queue)

pp <- NA
cl <- "grey30"
old.par <- par(mar=c(5, 5, 4, 3))
plot(tt, queue, type="s", xlab="Time (min)", 
  ylab="# (jobs in the system)", xlim=c(0,duration), 
  ylim=c(0,mq), col=cl, axes=F)
lines(c(tt[lq],duration), rep(queue[lq], 2), lty=1)
axis(1, pos=pp); axis(2, pos=pp)
par(old.par)


# UE - Aufgaben
#===============

# ue10.1
#========
x <- 2:12
px <- c(1:5,6,5:1)/36
(HX <- -sum(px*log2(px)))


# ue10.2
#========
(HY <- log2(6))
(H.YX <- (5/6)*(-(5/6)*log2(5/6)-(1/6)*log2(1/6)))
(HXY <- -5*(5/36)*log2(5/36)-(6/36)*log2(6/36)-5*(1/36)*log2(1/36))
# check
HY + H.YX; HXY


# ue10.3
#========
p <- seq(0, 1, by=0.001)
n <- 6
HX <- -n*(p*log2(p)+(1-p)*log2(1-p))
plot(p, HX, type="l", lty=1, lwd=3, 
  ylim=c(0,6), xlab="p", ylab="H(X)")
  

# ue10.7
#========
x <- seq(0, 5, by=0.01)
plot(x, exp(-x), type="l", lty=1, lwd=3,
  xlab="x", ylab="")
lines(x, 1/(1+x)^2, lty=2, lwd=3)
legend("topright", legend=c(expression(e^-x), 
  expression(frac(1,(1+x)^2))), bty="n", 
  lty=c(1,2), lwd=3)
  

# ue10.10
#=========
X <- matrix(c(9,-1, 0, 0,
              0, 0,-1, 9,
              9,-5, 4, 0,
              0, 4,-5, 9,
              1, 1, 1, 1),
              ncol=4, byrow=TRUE)
Y <- c(0,0,0,0,1)

qr.solve(X, Y)

  
# ue10.11
#=========
P <- matrix(c(0,1/2,1/2,
              1/3,1/3,1/3,
              1/3,1/3,1/3),
              ncol=3, byrow=TRUE)
N <- 30
plot(1:N, MC.sim(N, P, 1+1), type="o", lty=1, 
  pch=19, cex=1, xlab="Tag", ylab="", axes=FALSE) 
abline(h=1:3, lty=2)
axis(1, pos=1, at=c(1,5,10,15,20,25,30))
axis(2, pos=0, at=c(0,1,2,3), 
  labels=c("","regnerisch","sonnig","bewölkt"))


# ue10.12
#=========
P <- matrix(c(0,1,0,0,
              1/9,4/9,4/9,0,
              0,4/9,4/9,1/9,
              0,0,1,0),
              ncol=4, byrow=TRUE)
              
N <- 500
SIM <- MC.sim(N, P)-1
head(SIM, n=100)
table(SIM)/N


# ue10.13
#=========
# Übergangsmatrix
#-----------------
P <- matrix(c(0,1/2,1/2,
              1/2,0,1/2,
              1,0,0),
              3, 3, byrow=TRUE)

# regulär/ergodisch
#-------------------
P %*% P
P %*% P %*% P
P %*% P %*% P %*% P

# Simulation
#------------
N <- 10000
SIM <- MC.sim(N, P)
head(SIM, n=100)
table(SIM)/N

# Lösung mittels qr.solve()
#---------------------------
X <- matrix(c(1,-1/2,-1,
              1/2,-1,0,
              1/2,1/2,-1,
              1,1,1),
              ncol=3, byrow=TRUE)
Y <- c(0,0,0,1)
 
qr.solve(X, Y)


# ue10.15
#=========
Delta <- 0.05  # min
p <- 0.1
tt <- seq(Delta, 60, by=Delta)
n <- length(tt)
N <- 5
X <- matrix(nrow=n, ncol=N)
for (i in 1:N) {
  X[,i] <- cumsum(rbinom(n, 1, p))
}
X <- rbind(rep(0,N), X)

matplot(c(0,tt), X, type="s", lty=1, lwd=2,
  xlab="Zeit (Minuten)", ylab="X(n)")
x <- seq(0, 60, by=0.05)
lines(x, x/Delta*p, lty=2, lwd=2)
lines(x, x/Delta*p + 2*sqrt(x/Delta*p*(1-p)), lty=1, lwd=2)
lines(x, x/Delta*p - 2*sqrt(x/Delta*p*(1-p)), lty=1, lwd=2)


# ue10.22
#=========
mm1 <- function(lambdaA, lambdaB) {
  rho <- lambdaA/lambdaS
  EX <- rho/(1-rho)
  EXw <- rho^2/(1-rho)
  EXs <- rho
  EW <- rho/(lambdaS*(1-rho))
  ER <- 1/(lambdaS*(1-rho))
  round(data.frame(EX=EX, EXw=EXw, EXs=EXs, EW=EW, ER=ER), 2)
}

lambdaA <- 5 # pro min
lambdaS <- 6 # pro min
mm1(lambdaA, lambdaB)

lambdaA <- 5.5 # pro min
lambdaS <- 6 # pro min
mm1(lambdaA, lambdaB)

# Vgl. für die Simulation der Warteschlange
# den Code zu Bsp 10.18 (M/M/1).
