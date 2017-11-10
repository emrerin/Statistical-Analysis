#===============================================================# 
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#                Kapitel 2: Wahrscheinlichkeit                  #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# Bsp 2.1 (eGGZ)
#================
N <- 1000
K <- 10
M <- matrix(ncol=K, nrow=N)
for (i in 1:N) {
  for (j in 1:K) {
    Augens <- sum(sample(1:6, size=2, replace=T))
    M[i,j] <- ifelse(Augens == 7, 1, 0)
  }
}
Mcum <- apply(M, 2, cumsum)
matplot(1:N, Mcum/(1:N), type="l", lty=1, lwd=2, log="x",
  xlab="Anzahl der Würfe n", 
  ylab=expression(paste(h[n],"(Augensumme = 7)")))
abline(h=1/6, lty=2)


# Bsp 2.5 (Rendevousproblem)
#============================

# Abb 2.3
#---------
# A: 10 min, B: 20 min
x <- c(0,1); y <- c(0,1)
plot(x, y, type="l", lwd=2, xlim=c(0,1), ylim=c(0,1), 
  pty="s", asp=1, xlab="Eintreffzeitpunkt A", 
  ylab="Eintreffzeitpunkt B", axes=F)
axis(1, pos=0)
axis(2, pos=0)
axis(3, pos=1)
axis(4, pos=1)
lines(c(0,0,1,1,0), c(0,1,1,0,0), type="l", lty=1, lwd=2)
polygon(c(0,0,5/6,1,1,1/3,0), c(0,1/6,1,1,2/3,0,0), lwd=2,
  col="lightgrey")
lines(c(0,1), c(0,1), lty=2, lwd=2)
text(locator(1), expression(y == x - 1/3), srt=45, 
  pos=4, cex=1.5)
text(locator(1), expression(y == x + 1/6), srt=45, 
  pos=4, cex=1.5)    
  
  
# Bsp 2.13
#==========

# Abb 2.6
#---------
x <- 0:4
P <- matrix(ncol=4, nrow=5)
for (i in 1:4) {
  P[,i] <- dbinom(x, size=i, prob=1/2)
}
matplot(x, P, type="o", lty=1, lwd=3, pch=21, 
  cex=2, col=1:4, xlab="Zahl der Köpfe", 
  ylab="Bedingte Wahrscheinlichkeiten")
legend("topright", paste("Augenzahl = ", 1:4),
  lty=1, col=1:4, lwd=3, pch=21, cex=1.1)  


# 2.16 Beispiele
#================

# Geburtstagsproblem
#====================
n <- 23 # Personen
nn <- 0:(n-1)
x <- (365-nn)/365
1-prod(x)
(odds <- (1-prod(x))/prod(x))


# Matchingproblem
#=================

# Exakte Verteilung
#-------------------
f <- function(N, k) {
  NN <- 0:(N-k)
  x <- (-1)^NN/factorial(NN)
  sum(x)/factorial(k)
}

# N = 5
round(sapply(0:5, f, N=5), 4)

# Simulation
#------------
N <- 25
B <- 100000
ueber <- numeric(B)
for (i in 1:B) {
  x <- sample.int(N, N)
  ueber[i] <- sum((1:N) == x)
}
(tab <- table(ueber)/B)
bpl <- barplot(tab, axis.lty=1, space=0.5, col="darkgrey",
  xlab="k", ylab="P( {exakt  k  Übereinstimmungen} )")


# Diagnostische Tests
#=====================
PPV <- function(prv, sens, spez) {
  (prv*sens)/(prv*sens+(1-prv)*(1-spez))  # P(D+|T+)
}
NPV <- function(prv, sens, spez) {
  ((1-prv)*spez)/((1-prv)*spez+prv*(1-sens))  # P(D-|T-)
}

# Abb 2.10
#----------
prv <- seq(0, 50, by=0.5)
plot(prv, PPV(prv/100, 0.95, 0.99)*100, type="l", lty=1, 
  lwd=2, col=1, xlab="Prävalenz (%)", ylab="PPV (%)")
  

# UE - Aufgaben
#===============

# ue2.1
#=======
demere <- function(B, plotit=TRUE) {
  versuche <- B
  n <- versuche*4
  x <- matrix(sample(1:6, n, replace=TRUE), nrow=versuche, ncol=4)
  sechs.in.4 <- apply(x==6, 1, any)
  freq.6.in.4 <- cumsum(sechs.in.4)/(1:versuche)
  n <- versuche*48
  x <- matrix(sample(1:6, n, TRUE), nrow=versuche, ncol=48)
  doppel.6.in.24 <- apply(x==6, 1, function(x) any(x[1:24] & x[25:48]))
  freq.doppel.6.in.24 <- cumsum(doppel.6.in.24)/(1:versuche)
  if (plotit) {
    plot(freq.6.in.4, ylim=0:1, log="x", type="l", bty="n", lwd=2,
      ylab="Rel. Häufigkeit", xlab="Versuche", col="black")
    lines(1:versuche, freq.doppel.6.in.24, lty=1, lwd=2, col="red")
    lines(c(1,versuche), c(0.5,0.5), lty=1)
    legend("topright", legend=c("6 bei 4 Würfen (1 Würfel)", 
      "Doppel-6 bei 24 Würfen (2 Würfel)"), lty=1, lwd=2, 
      col=c("black","red"))
  }
  invisible(list(freq.6.in.4=freq.6.in.4,
    freq.doppel.6.in.24=freq.doppel.6.in.24))
}

demere(10000)


# ue2.8
#=======
pp <- numeric(13)
for (i in 1:13) {
  n <- i
  nn <- 0:(n-1)
  x <- (12-nn)/12
  pp[i] <- prod(x)
}

round(data.frame(pn=pp, pnc=1-pp), 4)
  
  
# ue2.9
#=======
f <- function(N, k) {
  NN <- 0:(N-k)
  x <- (-1)^NN/factorial(NN)
  sum(x)/factorial(k)
}

N <- 8
round(sapply(0:8, f, N=8), 4)
factorial(8)*sapply(0:8, f, N=8)


# ue2.10
#========
# Rekursive Funktion
T.rek <- function(n) {
  if (n == 1) return(1)
  else { M <- 0
    for (k in 1:(n-1)) {
      M <- M + choose(n-1, k)*T.rek(k) }
    M <- M + 1 
    return(M) }
}

nm <- matrix(1:12, ncol=1)
data.frame(T = apply(nm, 1, T.rek))


# ue2.16
#========
# Simulation
#------------
# Es genügt zu überprüfen, ob die Summe der beiden 
# kleinsten Stücke größer als das größte Stück ist.
N <- 10000
DK <- numeric(N)
for (i in 1:N) {
  x <- sort(runif(2))
  y <- sort(c(x[1],x[2]-x[1],1-x[2]))
  DK[i] <- ifelse(y[1]+y[2] > y[3], 1, 0)
}
DKum <- cumsum(DK)
plot(1:N, DKum/(1:N), type="l", lty=1, lwd=2, 
  log="x", xlab="n", col="blue3",
  ylab=expression(paste(h[n],"(Dreieck konstruierbar)")))
abline(h=1/4, lty=2)


# ue2.23
#========
p <- seq(0, 1, by=0.001)
p.intakt <- 2*p^2+2*p^3-5*p^4+2*p^5
plot(p, p.intakt, type="l", lty=1, lwd=2, asp=1,
  xlab="p", ylab="P( System intakt )", axes=F)
lines(c(0,1), c(0,1), lty=2)
axis(1, pos=0)
axis(2, pos=0)
axis(3, pos=1)
axis(4, pos=1)


# ue 2.25
#=========
P <- numeric(7)
j <- 1:6
for (k in 1:7) {
  P[k] <- sum(dbinom(k-1, size=j, prob=1/2)) /6 
}

round(P, 4)
P *6*2^6  # P *384
