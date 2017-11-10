#===============================================================#
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#                 Kapitel 9: Regressionsanalyse                 #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------
# carsnew.txt
# anscombe.txt
# wirebond.txt
# normtemp.txt
# wine.txt

# Packages:
#-----------
# base


# Bsp 9.2
#=========
carsn <- read.table("carsnew.txt", header=TRUE)

attach(carsn)

plot(100/MPGHwy ~ CurbWeight, type="p", pch=19,
  xlab="weight", ylab="gpm (Highway)", col="grey50")
mod <- lm(100/MPGHwy ~ CurbWeight)
abline(mod, lwd=2)
coef(mod)
fit <- data.frame(y=100/MPGHwy, yhat=fitted(mod), e=resid(mod))
round(fit, 4)

# Check
sum(fit$e)
sum(fit$e*CurbWeight)

detach(carsn)


# Bsp 9.3
#=========
SST <- var(100/carsn$MPGHwy)*(dim(carsn)[1]-1)
anova(mod)


# Bsp 9.4
#=========
anova(mod)


# Bsp 9.5
#=========
summary(mod)


# Bsp 9.6
#=========
attach(carsn)

plot(100/MPGHwy ~ CurbWeight, type="p", pch=19,
  xlab="weight", ylab="gpm (Highway)", col="grey50")
mod <- lm(100/MPGHwy ~ CurbWeight)
x <- CurbWeight
CuWe.new <- data.frame(CurbWeight=pretty(range(x), 100))
pred.c <- predict(mod, CuWe.new, interval="confidence")
pred.p <- predict(mod, CuWe.new, interval="prediction")
matplot(CuWe.new, cbind(pred.c, pred.p[,-1]), type="l",
  lty=c(1,2,2,3,3), lwd=3, col=c(1,2,2,3,3), add=TRUE)
legend("bottomright", c("LS - line", "95% confidence band",
  "95% prediction band"), lty=1:3, lwd=3, col=1:3)

detach(carsn)


# Bsp 9.8
#=========
attach(carsn)

mod <- lm(100/MPGHwy ~ CurbWeight)

x <- CurbWeight
y <- 100/MPGHwy

# Ohne Beob. # 27
mod1 <- lm(y[-27] ~ x[-27])

# Abb 9.6
#---------
old.par <- par(mfrow=c(2,2))
plot(y ~ x, type="p", pch=19, col="grey50",
  xlab="weight", ylab="gpm (Highway)", main="(a)")
abline(mod, lwd=2)
abline(mod1, lty=2, lwd=2, col="red")
identify(x, y)
plot(rstandard(mod) ~ fitted(mod), type="p",
  pch=19, col="grey50", ylim=c(-4,4), main="(b)")
abline(h=0, lty=2, lwd=2)
identify(fitted(mod), rstandard(mod))
plot(rstandard(mod) ~ x, type="p", pch=19, col="grey50",
  ylim=c(-4,4), xlab="weight", main="(c)")
abline(h=0, lty=2, lwd=2)
identify(x, rstandard(mod))
par(old.par)

detach(carsn)


# Abb 9.7
#=========
dat <- read.table("anscombe.txt", header=TRUE, sep=",")
cx <- 0.8
old.par <- par(mfrow=c(2,2))
for (i in c(1,3,5,7)) {
  plot(dat[,i+1] ~ dat[,i], type="p", pch=19, xlab="x",
    ylab="y", xlim=c(0,20), ylim=c(2,13), cex=1.2,
    col="grey50")
  xax <- par("xaxp")
  fit.ls <- lm(dat[,i+1] ~ dat[,i])
  abline(fit.ls, lwd=2)
  mtext(bquote(hat(beta)[0] == .(round(coef(fit.ls)[1], 2))),
    line=0.5, at=xax[1]+(xax[2]-xax[1])*(0.5/5), cex=cx)
  mtext(bquote(hat(beta)[1] == .(round(coef(fit.ls)[2], 2))),
    line=0.5, at=xax[1]+(xax[2]-xax[1])*(2/5), cex=cx)
  mtext(bquote(R^2 == .(paste(round(summary(fit.ls)$r.sq*100, 1),
    "%"))), line=0.5, at=xax[1]+(xax[2]-xax[1])*(4/5), cex=cx)
 }
par(old.par)


# 9.2.4 Beispiele
#=================

# 1. Quadratisches Modell
#=========================
x <- 0:8
y <- c(1.2,16.5,28.9,23.1,81.7,120.3,132.5,197.6,283.8)
plot(y ~ x, type="p", pch=19, col="grey50", cex=1.3)

mod1 <- lm(y ~ x)
abline(mod1, lty=2, lwd=3)

mod2 <- lm(y ~ x + I(x^2))
x.new <- data.frame(x=pretty(range(x), 100))
pred2 <- predict(mod2, x.new, interval="none")
lines(x.new$x, pred2, type="l", lty=1, lwd=3)

mod3 <- lm(y ~ x + I(x^2)+ I(x^3))
pred3 <- predict(mod3, x.new, interval="none")
lines(x.new$x, pred3, type="l", lty=3, lwd=3)
legend("bottomright", c("linear","quadratic","cubic"),
  lty=c(2,1,3), lwd=3)

summary(mod1)
summary(mod2)
summary(mod3)


# 2. Multiples Modell
#=====================
bond <- read.table("wirebond.txt", header=TRUE)[,-1]

pairs(bond, pch=19, cex=1.8, col="grey50", gap=0)

mod <- lm(Strength ~ Length + Height, data=bond)
summary(mod)

old.par <- par(mfrow=c(2,2), oma=c(0,1,2,0))
plot(mod)
par(old.par)


# 3. Dummy - Variable
#=====================
# Celsius/Fahrenheit: T_C = (T_F-32)*5/9

(98.6-32)*5/9  # = 37°C

normtemp <- read.table("normtemp.txt", header=TRUE)

old.par <- par(oma=c(0,1,0,1))

plot(temp ~ hr, type="p", pch=(3-gender)+20, cex=1.3,
  col=ifelse(gender==1, "blue", "red"), data=normtemp,
  xlab="heart rate", ylab="temperature (°F)")
CC <- seq(35.6, 39, by=0.2)
FF <- 9*CC/5+32
axis(4, at=FF, labels=CC, cex.axis=0.7)

mod0 <- lm(temp ~ hr, data=normtemp)
mod1 <- lm(temp ~ hr + factor(gender), data=normtemp)
summary(mod0)
summary(mod1)

cf1 <- coef(mod1)
abline(a=cf1[1], b=cf1[2], lty=1, lwd=3, col="blue")
abline(a=cf1[1]+cf1[3], b=cf1[2], lty=2, lwd=3, col="red")

legend(locator(1), c("male","female"), lty=c(1,2),
  col=c("blue","red"), lwd=3)
legend(locator(1), c("male","female"), pch=c(22,21),
  col=c("blue","red"))

par(old.par)

anova(mod0, mod1)


# UE - Aufgaben
#===============

# ue9.1
#=======
n <- 14
sumx <- 43; sumx2 <- 157.42
sumy <- 572; sumy2 <- 23530
sumxy <- 1697.80
xbar <- sumx/n
ybar <- sumy/n

opt <- options(digits=4)

(bet1 <- (sumxy - n*xbar*ybar)/(sumx2 - n*xbar^2))
(bet0 <- ybar - bet1*xbar)
(yprog <- bet0 + (3.7)*bet1)

(SST <- sumy2 - n*ybar^2)
(SSR <- bet1^2*(sumx2-n*xbar^2))
(SSE <- SST - SSR)
(MSE <- SSE/(n-2))


ANOVAtab <- data.frame(
   df=c(1,n-2,n-1),
   SS=c(SSR,SSE,SST),
   MS=c(SSR,SSE,NA)/c(1,n-2,n-1),
   Fval=c(SSR/MSE,NA,NA),
   pval=c(pf(SSR/MSE, 1, n-2, lower.tail=F),NA,NA),
   row.names=c("Regression","Error","Total"))

print.anova(ANOVAtab)

(R2 <- (SSR/SST)*100)

options(opt)


# ue9.4
#=======
y <- c(1,0,1,2,5,1,4,6,2,3,
       5,4,6,8,4,5,7,9,7,6)
x <- c(60,63,65,70,70,70,80,90,80,80,
       85,89,90,90,90,90,94,100,100,100)

# Scatterplot
plot(y ~ x, type="p", pch=19, cex=1.3, col="grey50",
  xlab="x [db]", ylab="y [mmHg]")

# Einfaches lineares Modell
mod <- lm(y ~ x)

# KQ--Gerade zeichnen
abline(mod, lty=1, lwd=2)

# F-Test, R^2, ...
summary(mod)

# Schätzwert für sig^2
summary(mod)$sigma^2

# Residualanalyse
old.par <- par(mfrow=c(2,2), oma=c(0,1,2,0))
plot(mod)
par(old.par)


# ue9.5
#=======
plot(y ~ x, type="p", pch=19, cex=1.3, col="grey50",
  xlab="x [db]", ylab="y [mmHg]" )
x.new <- data.frame(x=pretty(range(x), 100))
pred.c <- predict(mod, x.new, interval="confidence")
pred.p <- predict(mod, x.new, interval="prediction")
matplot(x.new, cbind(pred.c, pred.p[,-1]), type="l",
  lty=c(1,2,2,3,3), lwd=3, col=c(1,2,2,3,3), add=TRUE)
legend("bottomright", c("LS - line", "95% confidence band",
  "95% prediction band"), lty=1:3, lwd=3, col=1:3)


# ue9.7
#=======
y <- c(1.81,1.70,1.65,1.55,1.48,1.40,
       1.30,1.26,1.24,1.21,1.20,1.18)
x <- c(20,25,30,35,40,50,60,65,70,75,80,90)

# Scatterplot
plot(y ~ x, type="p", pch=19, cex=1.3, col="grey50",
  xlab="x", ylab="y", ylim=c(1.1,1.9))

# Quadratisches Modell
mod <- lm(y ~ x + I(x^2))
summary(mod)

# 95% Konfidenz- und Prognosebänder
x.new <- data.frame(x=pretty(range(x), 100))
pred.c <- predict(mod, x.new, interval="confidence")
pred.p <- predict(mod, x.new, interval="prediction")
matplot(x.new, cbind(pred.c, pred.p[,-1]), type="l",
  lty=c(1,2,2,3,3), lwd=3, col=c(1,2,2,3,3), add=TRUE)
legend("topright", c("LS - line", "95% confidence band",
  "95% prediction band"), lty=1:3, lwd=3, col=1:3)


# ue9.8
#=======
bond <- read.table("wirebond.txt", header=TRUE)[,-1]
mod1 <- lm(Strength ~ Length + Height, data=bond)

# Residualplots
plot(rstandard(mod1) ~ Length, type="p", pch=19,
  col="grey50", ylim=c(-3,3), data=bond)
x1 <- bond$Length; y <- rstandard(mod1)
lines(lowess(x1, y), lty=1, lwd=2, col="red")
abline(h=0, lty=2, lwd=2)
plot(rstandard(mod1) ~ Height, type="p", pch=19,
  col="grey50", ylim=c(-3,3), data=bond)
x2 <- bond$Height
lines(lowess(x2, y), lty=1, lwd=2, col="red")
abline(h=0, lty=2, lwd=2)

# Quadratische Abhängigkeit bez. Length
mod2 <- lm(Strength ~ Length + I(Length^2) + Height, data=bond)
summary(mod2)

# Residualanalyse
old.par <- par(mfrow=c(2,2), oma=c(0,1,2,0))
plot(mod2)
par(old.par)


# ue9.9
#=======
wine <- read.table("wine.txt", header=TRUE)

# Scatterplots
pairs(wine[1:29,-1], gap=0, panel=panel.smooth,
  pch=19, cex=1.3, lwd=2, col="darkblue")

# Regressionsmodell
mod <- lm(log(Price) ~ ., data=wine[1:29,-1])
summary(mod)

# Residualanalyse
old.par <- par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(mod)
par(old.par)

# Prognose
x.new <- wine[30:34, -c(1,6)]
PricePred <- predict(mod, newdata=x.new, conf="none")
y.new <- wine[30:34,6]
y.pred <- exp(PricePred)
data.frame(y.new, y.pred)

# Vergleich Price/PriceHat
ind <- complete.cases(wine)
comp <- data.frame(Price=wine$Price[ind],
  PriceHat=c(exp(fitted(mod)),y.pred))
plot(wine$Year[ind], comp$Price, type="o", pch=19, lty=1,
  lwd=2, xlab="", ylab="Price", cex.axis=0.8, col="blue1")
lines(wine$Year[ind], comp$PriceHat, type="o", pch=21,
  lty=2, lwd=2, col="blue4")
abline(v=1986, lty=2)
legend(locator(1), c("Price","PriceHat"), lty=1:2,
  lwd=2, pch=c(19,21), col=c("blue1","blue4"))
  
  
# ue9.10
#========
# Numerisches Beispiel
#----------------------
g <- function(x, xstar) ifelse(x > xstar, 1, 0)

# Regressionsfunktion
xstar <- 6
x <- seq(0, 10, by=0.01)
y <- 1 + 1.8*x - 1*(x - xstar)*g(x, xstar)

# Simulierte Beobachtungen
x1 <- seq(0, 10, by=0.1)
y1 <- 1 + 1.8*x1 - 1*(x1-xstar)*g(x1, xstar) + rnorm(length(x1), sd=1.3)
plot(y1 ~ x1, type="p", pch=19, cex=1.2, xlab="x", ylab="y")
lines(x, y, type="l", lty=2, lwd=2, col="red")
abline(v=xstar, lty=2)

# Modellanpassung
mod <- lm(y1 ~ x1 + I((x1-xstar)*g(x1, xstar)))
summary(mod)

# Prognose
x.new <- data.frame(x1=x)
y.pred <- predict(mod, newdata=x.new, conf="none")
lines(x, y.pred, lty=1, lwd=2, col="red")
legend("bottomright", c("true","LS"),
  lty=c(2,1), lwd=2, col="red")
