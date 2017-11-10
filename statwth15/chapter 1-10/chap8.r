#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#                 Kapitel 8: Bayes - Statistik                  #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------

# Packages:
#-----------
# base


# Bsp 8.6
#=========
# Equal-Tails
#-------------
qgamma(c(0.025,0.975), shape=6, rate=3)

# HPD
#-----
gam.hpd <- function(cpr, shape, rate, ...) {
  
  f <- function(x) {
    qgamma(x + cpr, shape, rate) - qgamma(x, shape, rate)
  }
  
  op <- optimize(f, c(0,1 - cpr), ...)
  list(cpr = cpr, param = c(shape, rate),
    lower = qgamma(op$minimum, shape, rate),
    upper = qgamma(op$minimum + cpr, shape, rate),
    length.min = op$objective,
    alpha = op$minimum)
}

gam.hpd(0.95, 6, 3)


# Bsp 8.7
#=========
x <- 0:20
p0 <- 0.5; p1 <- 0.5
b0 <- dbinom(x, 20, 0.05)
b1 <- dbinom(x, 20, 0.10)
post0 <- b0*p0/(b0*p0+b1*p1)
post1 <- b1*p1/(b0*p0+b1*p1)
round(data.frame(post0, post1), 4)
# Fehler 1. Art
sum(dbinom(3:20, 20, 0.05))


# UE - Aufgaben
#===============

# ue8.4
#=======
beta.hpd <- function(cpr, shape1, shape2, ...) {
  
  f <- function(x) {
    qbeta(x + cpr, shape1, shape2) - qbeta(x, shape1, shape2)
  }
  
  op <- optimize(f, c(0,1 - cpr), ...)
  list(cpr = cpr, param = c(shape1, shape2),
    lower = qbeta(op$minimum, shape1, shape2),
    upper = qbeta(op$minimum + cpr, shape1, shape2),
    length.min = op$objective,
    alpha = op$minimum)
}

n <- 8; x <- 3; a <- 1; b <- 1
# Equal-Tails
qbeta(c(0.025,0.975), shape1=a+x, shape2=b+n-x)
# HPD
beta.hpd(0.95, a+x, b+n-x)

n <- 8; x <- 3; a <- 1; b <- 2
# Equal-Tails
qbeta(c(0.025,0.975), shape1=a+x, shape2=b+n-x)
# HPD
beta.hpd(0.95, a+x, b+n-x)


# ue8.6
#=======
x <- c(11, 7,11, 6, 5, 9,14,10, 9, 5,
        8,10, 8,10,12, 9, 3,12,14, 4)
alphstar <- 10 + sum(x)
lambstar <- 5/6 + 20

# Equal-Tails
qgamma(c(0.025,0.975), shape=alphstar, rate=lambstar)  
 # HPD
gam.hpd(0.95, alphstar, lambstar)
# Test
(alpha0 <- pgamma(10, shape=alphstar, rate=lambstar))
alpha0/(1-alpha0)