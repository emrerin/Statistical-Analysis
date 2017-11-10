#===============================================================# 
#                                                               #
#                STAT u. WTH  f. INF    WS 2015                 #
#                                                               #
#===============================================================#
#                                                               #
#        Kapitel 1: Deskriptive u. Explorative Statistik        #
#                                                               #
#===============================================================#

# c/ W. Gurker


# Datenfiles:
#-------------
# TempWien1951-2012.txt
# body.txt
# beginner.txt
# IllmitzMai-Sept2010.txt
# IllmitzMai-Sept2011.txt
# IllmitzMai-Sept2012.txt
# UsingR::galton
# pkw-neuzul11.txt
# euroweight.txt

# Packages:
#-----------
# base, zoo, qcc, vioplot, UsingR, aplpack


# Bsp 1.4
#=========
library(zoo)
  
dat <- read.table("TempWien1951-2012.txt", header=T)
temp <- zoo(dat[,2:4], dat[,1])
ma <- rollmean(temp, k=10, align="right")

# Abb 1.1
#---------
plot(cbind(temp[,2],ma[,2]), screens=1, type="l", lwd=c(2,1), 
  lty=c(1,2), col="red2", xlab="", ylab="Temperatur °C",
  main="Jahreshöchsttemperaturen  Wien / Hohe Warte")


# Bsp 1.5
#=========
dat <- read.table("body.txt", header=T)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (0 - female, 1 - male)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
head(dat, n=10)
tail(dat, n=10)


# Bsp 1.6
#=========
dat <- read.table("beginner.txt", header=T, sep=";")
a <- c(3,5,6,8,10,12,15,21)
dat1 <- dat[-a,]
anf.w13 <- sort(dat1$W2013, decreasing=T, index.return=T)
anf.w13.sum <- sum(anf.w13$x)
stud <- dat1$Studienrichtung

# Häufigkeitstabelle
#--------------------
w13 <- data.frame(Studienrichtung=stud[anf.w13$ix],
  Absolut=anf.w13$x,
  Relativ=round(anf.w13$x/anf.w13.sum*100, 2),
  Kumuliert=round(cumsum(anf.w13$x/anf.w13.sum*100), 2))
rownames(w13) <- rownames(dat1)[anf.w13$ix]
w13

# Abb 1.2 (Kreisdiagramm)
#-------------------------
pie(anf.w13$x, labels=stud[anf.w13$ix], 
  main="Neuinskriptionen  TU-Wien  WS 2013",
  clockwise=FALSE, col=gray((10:26)/26),
  radius=0.8, cex=0.5)
  
# Abb 1.3 (Barplot)
#-------------------
dat2 <- dat[-a,c(8,6,4,2)]
dat3 <- as.matrix(dat2[anf.w13$ix,])
dat4 <- apply(dat3, 2, function(x) 100*x/sum(x)) 
barplot(dat3[,4], beside=T, col=gray((10:26)/26),
  legend.text=paste(stud[anf.w13$ix]," (",rownames(dat3),") ", sep=""),
  args.legend=list(cex=0.8), xlab="Studienrichtung",
  main="Neuinskriptionen  TU-Wien  WS 2013")

# Abb 1.4 (Barplot)
#-------------------
old.par <- par(no.readonly=T)
nf <- layout(matrix(c(1,2,1,0), 2, 2, byrow = TRUE),
  widths=c(3,2))
layout.show(nf)

par(mar=c(5, 4, 4, 2) + 0.1)
barplot(dat3, beside=F, col=gray((10:26)/26), xpd=T)
par(mar=c(1,1,2,1) + 0.1)
plot(x=0, y=0, type="n", axes=F, xlab="", ylab="",
  xlim=c(-1,1), ylim=c(-1,1), xpd=T)
legend("center", legend=stud[anf.w13$ix], cex=0.8,
 fill=gray((10:26)/26))
par(old.par)

# Abb 1.5 (Mosaikplot)
#----------------------
mosaicplot(dat3, col=gray((4:1)/5),
  main="Neuinskriptionen  TU-Wien  WS 2013")


# Bsp 1.7
#=========
install.packages("qcc")
require(qcc)
defects <- c(54,39,20,9,7,5,2)
names(defects) <- c("Insulating varnish","Loose leads",
  "Solder joint A", "Solder joint B","Resistor 1","Resistor 2",
  "Capacitor")

# Abb 1.6
#---------
pareto.chart(defects, ylab="Error frequency", xlab="", 
  main="Paretodiagramm", col="darkgray", 
  cumperc = seq(0, 100, by=10))


# Bsp 1.8
#=========
x <- round(rnorm(10), 2)  # Zufallszahlen (normal)
sort(x)
rank(x)


# Bsp 1.9
#=========
dat <- read.table("body.txt", header=T)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

ecdf.mf <- ecdf(dat[,1])
ecdf.m <- ecdf(datm[,1])
ecdf.f <- ecdf(datf[,1])

# Abb 1.7
#---------
old.par <- par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(ecdf.mf, do.p=F, verticals=T, xlab="Biacromial diameter [cm]",
  ylab=expression(hat(F)[n](x)), lwd=2, 
  main="Empirische Verteilungsfunktion: Biacromial")
lines(ecdf.m, do.p=F, verticals=T, col="darkblue", lwd=2)
lines(ecdf.f, do.p=F, verticals=T, col="red2", lwd=2)
legend(45, 0.2, c("female","male","both"),
  lty=1, col=c("red2","darkblue","black"), lwd=2)
par(old.par)


# Bsp 1.10
#==========
# Abb 1.8
#---------
stem(datm[,1], scale=2)
stem(datm[,1])


# Bsp 1.11
#==========
his <- hist(datm[,1], freq=F, plot=F)
table(cut(datm[,1], breaks=his$breaks))


# Bsp 1.12
#==========
# Abb 1.9 (Gender = 0 und Gender = 1)
#-------------------------------------
brk <- seq(29.5, 48.5, by=1)
old.par <- par(mfrow=c(1,2))
hist(datm[,1], breaks=brk, prob=T, col="darkgrey", right=F,
  xlab="Biacromial diameter [cm]", main="Male", xlim=c(30,50))
hist(datf[,1], breaks=brk, prob=T, col="darkgrey", right=F,
  xlab="Biacromial diameter [cm]", main="Female", xlim=c(30,50))
par(old.par)


# Abb 1.10
#==========
K1 <- function(x) dunif(x, min=-1, max=1)
K2 <- function(x) (1-abs(x))*dunif(x, min=-1, max=1)*2
K3 <- function(x) dnorm(x)
K4 <- function(x) (3/4)*(1-x^2)*dunif(x, min=-1, max=1)*2

x <- seq(-3, 3, by=0.01)
kern <- cbind(K1(x), K2(x), K3(x), K4(x))
matplot(x, kern, type="l", lty=1, xlab="z", ylab="K(z)",
  col=1:4, lwd=2)
legend("topright", c("Rechteck","Dreieck", "Normal","Epanechnikov"),
  lty=1, lwd=2, col=1:4)


# Bsp 1.13
#==========
dat <- rnorm(6)  # Zufallszahlen (normal)
n <- length(dat)
bb <- 1  
x <- seq(-5, 5, by=0.1)
kern <- c()
for (i in 1:n) {
  kern <- cbind(kern, dnorm(x, dat[i], bb))
}

# Abb 1.11
#----------
plot(density(dat, bw=bb), type="l", lty=1, lwd=2,
  main="Prinzip der Kerndichteschätzung")
rug(dat, ticksize=0.35, lty=2)
matplot(x, kern/(n*bb), type="l", lty=1, 
  col="darkgrey", add=T)


# Bsp 1.14
#==========
dat <- read.table("body.txt", header=T)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

# Abb 1.12 (Gender = 0 und Gender = 1)
#--------------------------------------
brk <- seq(29.5, 48.5, by=1)
old.par <- par(mfrow=c(1,2))
hist(datm[,1], breaks=brk, prob=T, col="lightgrey", right=F,
  xlab="Biacromial diameter [cm]", main="Male", xlim=c(30,50))
lines(density(datm[,1]), lty=1, lwd=2, col="darkblue")
hist(datf[,1], breaks=brk, prob=T, col="lightgrey", right=F,
  xlab="Biacromial diameter [cm]", main="Female", xlim=c(30,50))
lines(density(datf[,1]), lty=1, lwd=2, col="red3")
par(old.par)


# Bsp 1.15
#==========
x <- c(3,1,7,2,4,5,4,10,6,9)
quantile(x, c(5,25,50,75)/100, type=1)
quantile(x, c(5,25,50,75)/100, type=2)
quantile(x, c(5,25,50,75)/100, type=4)
quantile(x, c(5,25,50,75)/100, type=7)


# Bsp 1.16
#==========
fivenum(x)[c(2,4)]  # Hinges
quantile(x, c(25,50,75)/100)  # Quartile (Typ = 7)


# Bsp 1.18
#==========
x2010 <- read.table("IllmitzMai-Sept2010.txt")
x2011 <- read.table("IllmitzMai-Sept2011.txt")
x2012 <- read.table("IllmitzMai-Sept2012.txt")
x101112  <- cbind(x2010[,2], x2011[,2], x2012[,2])

# Abb 1.14
#----------
qqplot(x2010[,2], x2011[,2], type="o", pch=19, cex=0.8,
  xlim=c(30,190), ylim=c(30,190), pty="s",
  xlab="Quantile (2010)", ylab="Quantile (2011)",
  main="QQ - Plot")
abline(a=0, b=1, lty=2)


# Bsp 1.19
#==========
# Abb 1.15 (Boxplots)
#---------------------
boxplot(x101112, col="darkgrey", boxwex=0.5, 
  names=paste("20",c("10","11","12"), sep=""), 
  notch=T, main="Boxplots")
mtext(expression(paste(mu,g/m^3)), side=2, line=2.5)

# Abb 1.16 (Violinplots)
#------------------------
install.packages("vioplot")
require(vioplot)
vioplot(x2010[,2][!is.na(x2010[,2])], x2011[,2], x2012[,2], 
  names=paste("20",c("10","11","12"), sep=""), col="lightgrey")
mtext(expression(paste(mu,g/m^3)), side=2, line=2.5)
mtext("Violinplots", side=3, line=1)


# Zusätzliches Bsp (Violinplot)
#===============================
# Vergleich zweier Sortierungsalgorithmen
#-----------------------------------------
N <- 1000
Sim <- 50
rep <- 5000   # <<-- Eventuell verkleinern!
c1 <- c2 <- numeric(Sim)
for(is in 1:Sim){
  x <- rnorm(N)
  c1[is] <- system.time(for(i in 1:rep) sort(x, method="shell"))[1]
  c2[is] <- system.time(for(i in 1:rep) sort(x, method="quick"))[1]
  stopifnot(sort(x, method = "s") == sort(x, method = "q"))
}
erg <- cbind(ShellSort = c1, QuickSort = c2)
require(vioplot)
vioplot(erg[,1], erg[,2], names=c("Shell","Quick"), 
  col="lightgrey")
mtext("Seconds", side=2, line=2.5)


# Bsp 1.21
#==========
x <- c( 77, 87, 87,114,151,210,219,246,253,262,
       296,299,306,376,428,515,666,1310,2611)
boxplot(x, col="darkgrey", boxwex=0.8)
mean(x, trim=0.2)


# Bsp 1.22
#==========
(x <- 1:8)
(m <- median(x))
sort(abs(x-m))
mad(x, constant=1)
mad(x)
mad(x, constant=1, low=TRUE)
mad(x, constant=1, high=TRUE)


# Bsp 1.23
#==========
x2011 <- read.table("IllmitzMai-Sept2011.txt")
x2012 <- read.table("IllmitzMai-Sept2012.txt")
x101112  <- data.frame(x2010[,2], x2011[,2], x2012[,2])
colnames(x101112) <- c("2010","2011","2012")
x <- stack(x101112)
colnames(x) <- c("Ozone","Year")
attach(x)
by(Ozone, Year, summary)
detach(x)


# Bsp 1.24
#==========
dat <- read.table("body.txt", header=T)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)
nf <- dim(datf)[1]; nm <- dim(datm)[1]
brk <- seq(32, 50, by=1.2)
hist(dat[,1], breaks=brk, prob=T, col="lightgrey", right=F,
  xlab="Biacromial diameter [cm]", xlim=c(30,50),
  main="Mischverteilung")
bw <- 0.75
lines(density(dat[,1], bw=bw), lty=1, lwd=2, col=1)
lines(density(datm[,1], bw=bw)$x, density(datm[,1], bw=bw)$y*(nm/(nm+nf)), 
  lty=1, lwd=1, col="darkblue")
lines(density(datf[,1], bw=bw)$x, density(datf[,1], bw=bw)$y*(nf/(nm+nf)), 
  lty=1, lwd=1, col="red2")


# Bsp 1.25
#==========
# Abb 1.23
#----------
plot(Weight ~ Height, type="p", pch=24, data=dat,
  bg=ifelse(dat$Gender==0, "lightpink", "lightblue"),
  xlab="Height [cm]", ylab="Weight [kg]", main="Scatterplot")
legend("topright", c("Female","Male"), pch=24, 
  pt.bg=c("lightpink", "lightblue"))
  
# Abb 1.24
#----------
pairs(dat[,-6], gap=0, col=ifelse(dat$Gender==0, "lightpink3", 
  "lightblue3"))
pairs(dat[,-6], gap=0, col=ifelse(dat$Gender==0, "lightpink3", 
  "white"))
pairs(dat[,-6], gap=0, col=ifelse(dat$Gender==0, "white", 
  "lightblue3"))

# Histogramme in der Diagonalen
#-------------------------------
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="lightgrey", 
    border=1, ...)
}
pairs(datf[,-6], gap=0, diag.panel=panel.hist,
  pch=24, bg="lightpink3")
pairs(datm[,-6], gap=0, diag.panel=panel.hist,
  pch=24, bg="lightblue3")
pairs(dat[,-6], gap=0, diag.panel=panel.hist,
  pch=24, bg=ifelse(dat$Gender==0, "lightpink3",
  "lightblue"))


# Zusätzliches Bsp (Scatterplot)
#================================
install.packages("UsingR")
require(UsingR)
data(galton)
sunflowerplot(galton, main="Sunflowerplot")


# Bsp 1.26
#==========
dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

wx <- range(dat[,"Height"])
wxm <- range(datm[,"Height"])
wxf <- range(datf[,"Height"])
wy <- range(dat[,"Weight"])
wym <- range(datm[,"Weight"])
wyf <- range(datf[,"Weight"])

# Abb 1.25
#----------
par(mfrow=c(2,2))
plot(Weight ~ Height, type="p", pch=21, data=dat,
  bg=ifelse(dat$Gender==0, "lightpink3", "lightblue3"),
  xlim=wx, ylim=wy, xlab="Height [cm]", ylab="Weight [kg]",
  main="Female & Male", cex=0.8)

dens.estim <- kde2d(dat[,"Height"], dat[,"Weight"], n=50,
  lims=c(wx,wy))
image(dens.estim, xlab="Height [cm]", ylab="Weight [kg]",
  col=grey(5:10 /10), main="Female & Male")
contour(dens.estim, add=TRUE, drawlabels=FALSE)
lines(Weight ~ Height, type="p", pch=21, data=dat, cex=0.5,
  bg=ifelse(dat$Gender==0, "lightpink3", "lightblue3"))
  
dens.estim.f <- kde2d(datf[,"Height"], datf[,"Weight"], n=50,
  lims=c(wxf,wyf))
image(dens.estim.f, xlab="Height [cm]", ylab="Weight [kg]",
  col=grey(5:10 /10), main="Female")
contour(dens.estim.f, add=TRUE, drawlabels=FALSE)
lines(Weight ~ Height, type="p", pch=21, data=datf, cex=0.5,
  bg="lightpink3")
  
dens.estim.m <- kde2d(datm[,"Height"], datm[,"Weight"], n=50,
  lims=c(wxm,wym))
image(dens.estim.m, xlab="Height [cm]", ylab="Weight [kg]",
  col=grey(5:10 /10), main="Male")
contour(dens.estim.m, add=TRUE, drawlabels=FALSE)
lines(Weight ~ Height, type="p", pch=21, data=datm, cex=0.5,
  bg="lightblue3")
par(mfrow=c(1,1))


# Abb 1.26
#==========
# Gender = 0 und Gender = 1
par(mfrow=c(1,2))

plot(Weight ~ Height, type="p", pch=21, data=datf,
  bg=ifelse(datf$Gender==0, "lightpink3", "lightblue3"),
  xlab="Height [cm]", ylab="Weight [kg]", main="Female", cex=0.8)
abline(h=mean(datf[,"Weight"]), lty=1, col="darkgrey")
abline(v=mean(datf[,"Height"]), lty=1, col="darkgrey")
modf <- lm(Weight ~ Height, dat=datf)
abline(modf, lty=1, lwd=2, col="lightpink3")
identify(datf[,5], datf[,4])
  
plot(Weight ~ Height, type="p", pch=21, data=datm,
  bg=ifelse(datm$Gender==0, "lightpink3", "lightblue3"),
  xlab="Height [cm]", ylab="Weight [kg]", main="Male", cex=0.8)
abline(h=mean(datm[,"Weight"]), lty=1, col="darkgrey")
abline(v=mean(datm[,"Height"]), lty=1, col="darkgrey")
modm <- lm(Weight ~ Height, dat=datm)
abline(modm, lty=1, lwd=2, col="lightblue3")
identify(datm[,5], datm[,4])

par(mfrow=c(1,1))


# Bsp 1.27
#==========
round(cor(datm[,1:5]), 4)


# Bsp 1.28
#==========
x <- c(41, 37, 38, 39, 49, 47, 42, 34, 36, 48, 29)
y <- c(36, 20, 31, 24, 37, 35, 42, 26, 27, 29, 23)
cor(x, y)  # Pearson
cor(x, y, method="spearman")
(rx <- rank(x))
(ry <- rank(y))
cor(rx, ry)


# Ausprobieren!!
#================
g1 <- function(closed=T, xL, yL) {
  plot(c(0), c(0), type="n", xlim=xL, ylim=yL,
    xlab="x", ylab="y")
  abline(h=seq(xL[1], xL[2], by=0.1), lty=2, lwd=1, 
    col="lightgrey")
  abline(v=seq(yL[1], yL[2], by=0.1), lty=2, lwd=1, 
    col="lightgrey")
  z <- locator()
  if (closed) {
    plot(z$x, z$y, type="p", pch=19, cex=0.6, 
      xlim=xL, ylim=yL, xlab="x", ylab="y")
    polygon(z$x, z$y, lty=1, lwd=2)
    } else {
    plot(z$x, z$y, type="o", lty=1, pch=19, cex=0.6, 
      lwd=2, xlim=xL, ylim=yL, xlab="x", ylab="y")
  }
  z
}
g2 <- function(x, y, mt=T) {
  rp <- cor(x, y)
  rs <- cor(x, y, method="spearman")
  if (mt) {
    mtext(paste("r.pearson = ", round(rp, 3),
        "       r.spearman = ", round(rs, 3)),
      cex=0.8, line=1)
  }
}
g <- function(mt=T, closed=F, xL, yL) {
  z1 <- g1(closed, xL, yL)
  z2 <- locator()
  z <- rbind(data.frame(z1), data.frame(z2))
  lines(z$x, z$y, type="p", pch=19, col=2)
  abline(lm(z$y ~ z$x), lty=2, lwd=2, col=2)
  g2(z$x, z$y, mt)
}

# Wählen Sie einige Punkte! (Jeweils mit 
# rechter Maustaste beenden.)
# Polygonzug (offen)
g(mt=T, closed=F, xL=c(0,1), yL=c(0,1))
# Polygonzug (geschlossen)
g(mt=T, closed=T, xL=c(0,1), yL=c(0,1))


# KQ-Gerade Abb 1.26
#====================
modm <- lm(Weight ~ Height, dat=datm)
coef(modm)



# UE - Aufgaben
#===============

# ue1.1
#=======
require(zoo)
dat <- read.table("TempWien1951-2012.txt", header=T)
temp <- zoo(dat[,2:4], dat[,1])
ma <- rollmean(temp, k=10, align="right")

old.par <- par(mfrow=c(2,2), mar=c(3, 4, 3, 2) + 0.1)
plot(temp, screens=1, type="l", lwd=2, lty=1, col=c(1,2,4), 
  xlab="", ylab="Temperatur °C")
plot(cbind(temp[,1],ma[,1]), screens=1, type="l", lwd=c(2,1), 
  lty=c(1,2), col=1, xlab="", ylab="Temperatur °C")
legend("bottomright", "Durchschnitt")
plot(cbind(temp[,2],ma[,2]), screens=1, type="l", lwd=c(2,1), 
  lty=c(1,2), col=2, xlab="", ylab="Temperatur °C")
legend("bottomright", "Maximum")
plot(cbind(temp[,3],ma[,3]), screens=1, type="l", lwd=c(2,1), 
  lty=c(1,2), col=4, xlab="", ylab="Temperatur °C")
legend("bottomright", "Minimum")
par(old.par)
title("Lufttemperatur Wien/Hohe Warte 1951-2012", 
  outer=T, line=-0.8)


# ue1.2
#=======
pkw <- read.table("pkw-neuzul11.txt", header=TRUE, sep=";")
TO <- sum(pkw$TOTAL)
TO.other <- pkw$TOTAL[pkw$GROUP == "OTHER"]
pkw2 <- subset(pkw, (GROUP != "OTHER")&(100*pkw$TOTAL/TO >= 3))
TO2 <- sum(pkw2$TOTAL)
pkw2.other <- data.frame(GROUP="OTHER", TOTAL=TO-TO2)
(pkw3 <- rbind(pkw2, pkw2.other))
ra <- sort(pkw3$TOTAL, index.return=TRUE)
pie(pkw3$TOTAL[ra$ix], labels=pkw3$GROUP[ra$ix],
  col=gray(seq(0.4, 1.0, length=12)),
  main="PKW Neuzulassungen 2011 (Western Europe)")


# ue1.3
#=======
require(qcc)

# Neuinskriptionen W2013
#------------------------
dat <- read.table("beginner.txt", header=TRUE, sep=";")
neu13 <- dat[,2]
names(neu13) <- dat$Studienrichtung
prta <- pareto.chart(neu13, ylab="Neusinskiptionen W2013", 
  cumperc = seq(0, 100, by=10), cex.names=0.7)
round(prta, 1)

# Pareto für Neuzulassungen
#---------------------------
pkw <- read.table("pkw-neuzul11.txt", header=TRUE, sep=";")
neuzul11 <- pkw$TOTAL
names(neuzul11) <- pkw$GROUP
prtb <- pareto.chart(neuzul11, cex.names=0.7, 
  cumperc=seq(0, 100, by=10))
round(prtb, 1)


# ue1.4
#=======
x <- c(0,2,0,0,1,3,0,3,1,1,0,0,1,2,0,
       0,0,1,1,3,0,1,0,0,0,5,1,0,2,0)
n <- length(x)
tab <- table(x)
mr <- as.numeric(names(tab))
n1 <- min(mr)
n2 <- max(mr)
nk <- numeric(n2-n1+1)
nk[mr-n1+1] <- tab
names(nk) <- n1:n2

old.par <- par(mar=c(5, 4.5, 4, 2) + 0.1)
barplot(nk/n, axis.lty=1, space=2, xlab="Zahl der Fehler",
  ylab="Relative  Häufigkeit", main="Balkendiagramm", 
  col="darkgray")
plot(ecdf(x), verticals=TRUE, do.points=FALSE, lwd=2,
  xlab="Zahl der Fehler", ylab=expression(hat(F)(x)),
  main="Summentreppe")
par(old.par)


# ue1.5
#=======
dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

ecdf.mf <- ecdf(dat[,2])
ecdf.m <- ecdf(datm[,2])
ecdf.f <- ecdf(datf[,2])

old.par <- par(mar=c(5, 4.5, 4, 2) + 0.1)
plot(ecdf.mf, do.p=F, verticals=T, xlab="Waist girth [cm]",
  ylab=expression(hat(F)[n](x)), lwd=1.5, 
  main="Empirische Verteilungsfunktion")
lines(ecdf.m, do.p=F, verticals=T, col="darkblue", lwd=1.5)
lines(ecdf.f, do.p=F, verticals=T, col="red2", lwd=1.5)
legend("bottomright", c("female","male","both"),
  lty=1, col=c("red2","darkblue","black"), lwd=1.5)
par(old.par)

# Stem-and-leaf-Plot
#--------------------
stem(datm[,1], scale=2)
stem(datf[,1], scale=2)
stem(datf[,1])

# Back-to-Back Stem-and-Leaf Plot
#---------------------------------
install.packages("aplpack")
require(aplpack)
stem.leaf.backback(datf[,1], datm[,1], m=2, depth=F)


# ue1.7
#=======
dat <- read.table("euroweight.txt", header=TRUE, skip=1)
attach(dat)
brk <- seq(7.200, 7.760, by=0.010)
par(mfrow=c(4,2))
for (i in 1:8) {
 subdat <- weight[batch==i] 
 hist(subdat, breaks=brk, freq=FALSE, main=paste("Batch",i), 
    xlim=c(7.35,7.65), xlab="weight [g]", col="lightgrey",
    cex.main=1)
 lines(density(subdat), lwd=2)
 }
par(mfrow=c(1,1))
detach(dat)


# ue1.8
#=======
dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

# Boxplot
boxplot(Biacromial ~ Gender, data=dat, col="darkgrey",
  names=c("female","male"), pars=list(boxwex=0.6))
mtext("Biacromial [cm]", side=2, line=2.5)
mtext("Boxplots", side=3, line=1)

# Violinplot
require(vioplot)
vioplot(datf$Biacromial, datm$Biacromial, col="lightgrey",
  names=c("female","male"))
mtext("Biacromial [cm]", side=2, line=2.5)
mtext("Violinplots", side=3, line=1)

# Kennzahlen
kennz <- function(x, ro=4) {
  param <- c(mean(x), median(x), max(x)-min(x), var(x), 
    sd(x), sd(x)/mean(x), IQR(x), mean(abs(x-mean(x))), 
    mean(abs(x-median(x)))) 
  param.m <- matrix(param, ncol=1)
  dimnames(param.m) <- list(c("Mittel","Median","Spannweite",
    "Varianz","Streuung","VarKoef","IQR","MAD.Mittel",
    "MAD.Median"), "")
  round(param.m, ro)
}

attach(dat)
by(Biacromial, Gender, summary)
by(Biacromial, Gender, fivenum)
by(Biacromial, Gender, kennz, ro=3)
detach(dat)


# ue1.9
#=======
fivenum(1:100)
fivenum(c(1:75, (76:100)*10))


# ue1.10
#========
yr <- seq(1951, 2011, by=10)
x <- c(7568,9315,10062,10102,10349,11334,12995)
bpl <- barplot(x, axis.lty=1, space=2, xlab="", 
  ylab="Bevölkerung", names.arg=yr, main="Eisenstadt 1951 - 2011", 
  col="darkgray", ylim=c(0,13500))

# 10-jährliche Zunahme
#----------------------
(zu10 <- ((x[7]/x[1])^(1/6) - 1)*100)
x[1]*(1 + zu10/100)^6
lines(bpl, x[1]*(1 + zu10/100)^(0:6), type="o",
  lty=1, pch=19, lwd=2)

# 1-jährliche Zunahme
#---------------------
(zu1 <- ((x[7]/x[1])^(1/60) - 1)*100)
x[1]*(1 + zu1/100)^60

# Prognose 2030
#---------------
(zu1a <- ((x[7]/x[6])^(1/10) - 1)*100)
x[7]*(1 + zu1a/100)^19
x[6]*(1 + zu1a/100)^29

# ue1.11

avg <- (10+20+50+20) / (10/35000 + 20/42000 + 50/52500 + 20/28000)

# ue1.14
#========
attach(dat)
by(Height, Gender, kennz, ro=3)
detach(dat)


# ue1.15
#========
dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
dat1 <- within(dat, { BMI <- Weight/(Height/100)^2 })
dat1f <- subset(dat1, Gender==0)
dat1m <- subset(dat1, Gender==1)

# Boxplots
#----------
boxplot(BMI ~ Gender, data=dat1, col="darkgrey",
  names=c("female","male"), pars=list(boxwex=0.6))
mtext("BMI [ kg/m² ]", side=2, line=2.5)
mtext("Boxplots", side=3, line=1)

# Histogramme
#-------------
brk <- seq(16, 40, by=2)
old.par <- par(mfrow=c(1,2))
hist(dat1m[,7], breaks=brk, prob=T, col="darkgrey", right=F,
  xlab="kg/m²", main="BMI: Men", xlim=c(15,40))
lines(density(dat1m[,7]), lwd=2)
abline(v=25, lty=2, lwd=2)
x25m <- sum(dat1m[,7] > 25)/length(dat1m[,7])*100
text(25, 0.15, labels=paste("> 25: ", round(x25m,1), "%"), 
  pos=4, cex=1.1) 
hist(dat1f[,7], breaks=brk, prob=T, col="darkgrey", right=F,
  xlab="kg/m²", main="BMI: Women", xlim=c(15,40))
lines(density(dat1f[,7]), lwd=2)
abline(v=25, lty=2, lwd=2)
x25f <- sum(dat1f[,7] > 25)/length(dat1f[,7])*100
text(25, 0.15, labels=paste("> 25: ", round(x25f, 1), "%"), 
  pos=4, cex=1.1) 
par(old.par)

# Kennzahlen
#------------
bowley <- function(x) {
  Q <- quantile(x, (1:3)/4)
  (Q[3]-2*Q[2]+Q[1])/(Q[3]-Q[1])
}
moors <- function(x) {
  A <- quantile(x, (1:7)/8)
  ((A[7]-A[5])+(A[3]-A[1]))/(A[6]-A[2])
} 

kennz2 <- function(x, ro=4) {
  param <- c(mean(x), median(x), max(x)-min(x), var(x), 
    sd(x), sd(x)/mean(x), IQR(x), mean(abs(x-median(x))),
    e1071::skewness(x), e1071::kurtosis(x), bowley(x),
    moors(x)) 
  param.m <- matrix(param, ncol=1)
  dimnames(param.m) <- list(c("Mittel","Median","Spannweite",
    "Varianz","Streuung","VarKoef","IQR","MAD","Schiefe",
    "Exzess","Bowley","Moors"), "")
  round(param.m, ro)
}

attach(dat1)
by(BMI, Gender, kennz2, ro=3)
detach(dat1)


# ue1.16
#========
require(UsingR)
data(brightness)
hist(brightness, breaks=25, prob=T, xlab="mag",
  main="", col="darkgrey")
lines(density(brightness), lwd=2)
  
kennz2(brightness)


# ue1.17
#========
dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)

attach(dat)
print(by(dat[,1:5], Gender, cor), digits=4)
detach(dat)


# ue1.18
#========

dat <- read.table("body.txt", header=TRUE)[,c(1,12,22,23,24,25)]
#  V1 = Biacromial diameter (Schulter) (cm)
# V12 = Waist girth (Taille) (cm)
# V22 = Age (years)
# V23 = Weight (kg)
# V24 = Height (cm)
# V25 = Gender (1 - male, 0 - female)
names(dat) <- c("Biacromial","Waist","Age","Weight","Height","Gender")
datf <- subset(dat, Gender==0)
datm <- subset(dat, Gender==1)


xR <- range(dat$Waist)
yR <- range(dat$Weight)

par(mfrow=c(1,2))

plot(Weight ~ Waist, type="p", pch=21, data=datf,
  bg=ifelse(datf$Gender==0, "lightpink3", "lightblue3"),
  xlab="Waist [cm]", ylab="Weight [kg]",
  main="Female", cex=0.8, xlim=xR, ylim=yR)
modf <- lm(Weight ~ Waist, dat=datf)
abline(modf, lty=1, lwd=2, col="lightpink3")
  
plot(Weight ~ Waist, type="p", pch=21, data=datm,
  bg=ifelse(datm$Gender==0, "lightpink3", "lightblue3"),
  xlab="Waist [cm]", ylab="Weight [kg]",
  main="Male", cex=0.8, xlim=xR, ylim=yR)
modm <- lm(Weight ~ Waist, dat=datm)
abline(modm, lty=1, lwd=2, col="lightblue3")

par(mfrow=c(1,1))

coef(modf)
coef(modm)


# ue1.19
#========
x <- c(3,1,5,6,3,4)
y <- c(4,2,4,8,6,5)
plot(y ~ x, type="p", pch=21, bg="lightblue3",
  xlim=c(0,6), ylim=c(0,8), cex=2,
  main="KQ - Gerade durch den Nullpunkt")
mod1 <- lm(y ~ -1 + x)
coef(mod1)
abline(mod1, lwd=2, col="lightblue3")
mod2 <- lm(y ~ x)
coef(mod2)
abline(mod2, lty=2, lwd=2, col="lightblue3")


# ue1.20
#========
x <- c(-2,3,-1,0,-3,1,5,-3)
y <- c(7,15,3,1,11,6,20,16)
plot(y ~ x, type="p", pch=21, bg="lightblue3",
  xlim=c(-5,5), ylim=c(0,25), cex=2,
  main="KQ - Parabel")
mod <- lm(y ~ I(x^2))
coef(mod)
xnew <- data.frame(x=seq(-5, 5, by=0.1))
ypred <- predict(mod, newdata=xnew, interval="n")
lines(xnew$x, ypred, lwd=2, col="lightblue3")






























































