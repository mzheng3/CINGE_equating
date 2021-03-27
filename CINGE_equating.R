#-----------------------------------------------------------------
# PSQF 7358 Equating
# Final Project_Mingying Zheng




library(irtoys)  # for IRT calibration
library(equate)	# for traditional equating
library(plink)	# for IRT linking and equating
detach("package:equate", unload = TRUE)
#-----------------------------------------------------------------
# IRT simulation

# Example research questions:
# How well does the 1PL perform with data generated under
# a 3PL with non-negligible c parameters?
# Does 1PL performance differ by ability?

# Generate fixed item parameters for 3PL
# Modify seed, ni, nvi, and parameter distributions
# but keep the 3PL, that is, don't set c = 0 or a = 1
set.seed(10000)
ni <- 80
nvi1 <- 16
nvi2<-20
nvi3<-24

xabc <- data.frame(a = 1, b = rnorm(ni), c = 0)
yabc <- data.frame(a = 1, b = rnorm(ni), c = 0)

#xpar<-cbind(a=rnorm(), b=rnorm(), c=rnorm())
#Ypar<-cbind(a=rnorm(), b=rnorm(), c=rnorm())
#A<-12
#B<-7
#ypar[1-10, "6"]<-xpars[1-10, "6"]
# Set linking constants
# Linearly transform scales
# Optional - if you want to look at recovery of A and B
# A <- 1.2
# B <- .8
# yabc[1:nvi, 1] <- xabc[1:nvi, 1]/A
# yabc[1:nvi, 2] <- xabc[1:nvi, 2]*A + B

# Set n = 1000, ability fixed
nj <- 1000
xm <- .5
#xm<- rnorm(0,1)
ym <- .5
#xm<-rnorm(0,1)
xs <- ys <- 1

# True equating functions_nvi1=16
irtlink1 <- plink(as.irt.pars(
  x = list(yabc, xabc),
  common = cbind(1:nvi1, 1:nvi1),
  poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
  cat = list(rep(2, ni), rep(2, ni))),
  rescale = "SL")
irtscale1 <- equate(irtlink1, method = c("TSE", "OSE"),
                    base.grp = 2)
tscale1<- irtscale1$tse[, 3]
oscale1 <- irtscale1$ose$scores[, 2]

# True equating functions_nvi2=20
irtlink2 <- plink(as.irt.pars(
  x = list(yabc, xabc),
  common = cbind(1:nvi2, 1:nvi2),
  poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
  cat = list(rep(2, ni), rep(2, ni))),
  rescale = "SL")
irtscale2 <- equate(irtlink2, method = c("TSE", "OSE"),
                    base.grp = 2)
tscale2<- irtscale2$tse[, 3]
oscale2 <- irtscale2$ose$scores[, 2]

# True equating functions_nvi3=24
irtlink3 <- plink(as.irt.pars(
  x = list(yabc, xabc),
  common = cbind(1:nvi3, 1:nvi3),
  poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
  cat = list(rep(2, ni), rep(2, ni))),
  rescale = "SL")
irtscale3 <- equate(irtlink3, method = c("TSE", "OSE"),
                    base.grp = 2)
tscale3<- irtscale3$tse[, 3]
oscale3 <- irtscale3$ose$scores[, 2]

# A small simulation study
# Comparing separate and concurrent linking
# Modify as needed
reps <-500

stscale1 <- soscale1 <- ctscale1 <- coscale1 <-
  stscale2 <- soscale2 <- ctscale2 <- coscale2 <-
  stscale3 <- soscale3 <- ctscale3 <- coscale3 <-
  matrix(nrow = ni + 1, ncol = reps)
for(i in 1:reps) {
  
  # Simulate item responses
  # Note that ability is a random effect_nvi1=20
  xresp1 <- sim(xabc, rnorm(nj, xm, xs))
  yresp1 <- sim(yabc, rnorm(nj, ym, ys))
  colnames(xresp1) <-
    paste("x", 1:ni, sep = "")
  colnames(yresp1) <-
    paste("y", 1:ni, sep = "")
  colnames(xresp1)[1:nvi1] <- 
    colnames(yresp1)[1:nvi1] 
  paste("c", 1:nvi1, sep = "")
  
  # Note that ability is a random effect_nvi2=16
  xresp2 <- sim(xabc, rnorm(nj, xm, xs))
  yresp2 <- sim(yabc, rnorm(nj, ym, ys))
  colnames(xresp2) <-
    paste("x", 1:ni, sep = "")
  colnames(yresp2) <-
    paste("y", 1:ni, sep = "")
  colnames(xresp2)[1:nvi2] <- 
    colnames(yresp2)[1:nvi2] 
  paste("c", 1:nvi2, sep = "")
  
  
  xresp3 <- sim(xabc, rnorm(nj, xm, xs))
  yresp3 <- sim(yabc, rnorm(nj, ym, ys))
  colnames(xresp3) <-
    paste("x", 1:ni, sep = "")
  colnames(yresp3) <-
    paste("y", 1:ni, sep = "")
  colnames(xresp3)[1:nvi3] <- 
    colnames(yresp3)[1:nvi3] 
  paste("c", 1:nvi3, sep = "")
  # Sparse matrices for concurrent calibration
  sresp1 <- testEquatingData(list(xresp1, yresp1))
  sresp2 <- testEquatingData(list(xresp2, yresp2))
  sresp3 <- testEquatingData(list(xresp3, yresp3))
  # Models
  xmodel1 <- rasch(xresp1, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  ymodel1 <- rasch(yresp1, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  yxmodel1 <- rasch(sresp1, constraint = cbind(ncol(sresp1) + 1, 1),
                    IRT.param = T)
  
  xmodel2 <- rasch(xresp2, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  ymodel2 <- rasch(yresp2, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  yxmodel2 <- rasch(sresp2, constraint = cbind(ncol(sresp2) + 1, 1),
                    IRT.param = T)
  
  xmodel3 <- rasch(xresp3, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  ymodel3 <- rasch(yresp3, constraint = cbind(ni + 1, 1),
                   IRT.param = T)
  yxmodel3 <- rasch(sresp3, constraint = cbind(ncol(sresp3) + 1, 1),
                    IRT.param = T)
  # Parameters
  xparams1 <- data.frame(coef(xmodel1)[, 2:1], 0)
  yparams1 <- data.frame(coef(ymodel1)[, 2:1], 0)
  xyparams1 <- data.frame(coef(yxmodel1)[, 2:1], 0)
  
  xparams2 <- data.frame(coef(xmodel2)[, 2:1], 0)
  yparams2 <- data.frame(coef(ymodel2)[, 2:1], 0)
  xyparams2 <- data.frame(coef(yxmodel2)[, 2:1], 0)
  
  xparams3 <- data.frame(coef(xmodel3)[, 2:1], 0)
  yparams3<- data.frame(coef(ymodel3)[, 2:1], 0)
  xyparams3 <- data.frame(coef(yxmodel3)[, 2:1], 0)
  
  # Separate linking with SL_nvi1=20
  seplink1 <- plink(as.irt.pars(
    x = list(yparams1, xparams1),
    common = cbind(1:nvi1, 1:nvi1),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni))),
    rescale = "SL")
  # Separate linking with SL_nvi2=16
  seplink2 <- plink(as.irt.pars(
    x = list(yparams1, xparams1),
    common = cbind(1:nvi2, 1:nvi2),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni))),
    rescale = "SL")
  #nvi=24
  seplink3 <- plink(as.irt.pars(
    x = list(yparams3, xparams3),
    common = cbind(1:nvi3, 1:nvi3),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni))),
    rescale = "SL")
  # Concurrent linking
  conlink1 <- as.irt.pars(
    x = list(xyparams1[colnames(yresp1), ], xyparams1[colnames(xresp1), ]),
    common = cbind(1:ni, 1:ni),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni)))
  
  conlink2 <- as.irt.pars(
    x = list(xyparams1[colnames(yresp1), ], xyparams1[colnames(xresp1), ]),
    common = cbind(1:ni, 1:ni),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni)))
  
  conlink3 <- as.irt.pars(
    x = list(xyparams3[colnames(yresp3), ], xyparams1[colnames(xresp3), ]),
    common = cbind(1:ni, 1:ni),
    poly.mod = list(as.poly.mod(ni), as.poly.mod(ni)),
    cat = list(rep(2, ni), rep(2, ni)))
  
  # Equating, mean equating should work best.
  # smoothing should be made, such as loglinear
  sscale1 <- equate(seplink1, method = c("TSE", "OSE"), base.grp = 2)
  cscale1 <- equate(conlink1, method = c("TSE", "OSE"), base.grp = 2)
  
  
  stscale1[, i] <- sscale1$tse[, 3]
  soscale1[, i] <- sscale1$ose$scores[, 2]
  ctscale1[, i] <- cscale1$tse[, 3]
  coscale1[, i] <- cscale1$ose$scores[, 2]
  
  
  sscale2 <- equate(seplink2, method = c("TSE", "OSE"), base.grp = 2)
  cscale2 <- equate(conlink2, method = c("TSE", "OSE"), base.grp = 2)
  
  
  stscale2[, i] <- sscale2$tse[, 3]
  soscale2[, i] <- sscale2$ose$scores[, 2]
  ctscale2[, i] <- cscale2$tse[, 3]
  coscale2[, i] <- cscale2$ose$scores[, 2]
  
  sscale3 <- equate(seplink3, method = c("TSE", "OSE"), base.grp = 2)
  cscale3 <- equate(conlink3, method = c("TSE", "OSE"), base.grp = 2)
  
  
  stscale3[, i] <- sscale3$tse[, 3]
  soscale3[, i] <- sscale3$ose$scores[, 2]
  ctscale3[, i] <- cscale3$tse[, 3]
  coscale3[, i] <- cscale3$ose$scores[, 2]
  
  # Capture other information before moving on?
  # Item parameter recovery?
  # Correlations with generating parameters?
}


# Look at results
# Here's a start to some plots

# Ploting SEE_nvi1-nvi2)TSE vs. OSE; Separate vs. Concurrent 
#plotting SEE_nvi1=16
stsee1 <- apply(stscale1, 1, sd)
sosee1 <- apply(soscale1, 1, sd)
ctsee1 <- apply(ctscale1, 1, sd)
cosee1 <- apply(coscale1, 1, sd)
plot(main="SEE for Four Rasch Equating Methods (v1=16)", c(0, 81), c(0, .5), type = "n", xlab="item")
lines(0:80, stsee1, col = 2)
lines(0:80, sosee1, col = 3)
lines(0:80, ctsee1, col = 4)
lines(0:80, cosee1, col = 5)
legend("topright",legend=c("stsee1","sosee1","ctsee1", "cosee1"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)

meansee1 <- lapply(list(stsee1, sosee1, ctsee1,cosee1), function(x)
  mean(x, na.rm = T))
#plotting SEE_nvi2=20
stsee2 <- apply(stscale2, 1, sd)
sosee2 <- apply(soscale2, 1, sd)
ctsee2 <- apply(ctscale2, 1, sd)
cosee2 <- apply(coscale2, 1, sd)

plot(main="SEE for for Four Rasch Equating Methods (v2=20)", c(0, 81), c(0, .5), type = "n", xlab="item")
lines(0:80, stsee2, col = 2)
lines(0:80, sosee2, col = 3)
lines(0:80, ctsee2, col = 4)
lines(0:80, cosee2, col = 5)
legend("topright",legend=c("stsee2", "sosee2","ctsee2","cosee2"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)

meansee2 <- lapply(list(stsee2, sosee2, ctsee2,cosee2), function(x)
  mean(x, na.rm = T))
#nvi3=24
stsee3 <- apply(stscale3, 1, sd)
sosee3 <- apply(soscale3, 1, sd)
ctsee3 <- apply(ctscale3, 1, sd)
cosee3 <- apply(coscale3, 1, sd)
plot(main="SEE for Four Rasch Equating Methods (v3=24)", c(0, 81), c(0, 1.0), type = "n", xlab="item")
lines(0:80, stsee3, col = 2)
lines(0:80, sosee3, col = 3)
lines(0:80, ctsee3, col = 4)
lines(0:80, cosee3, col = 5)
legend("topright",legend=c("stsee3","sosee3","ctsee3", "cosee3"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)
meansee3 <- lapply(list(stsee3, sosee3, ctsee3,cosee3), function(x)
  mean(x, na.rm = T))
# Ploting Bias_nvi1-nvi2)TSE vs. OSE; Separate vs. Concurrent
#nvi1=16

stbias1 <- apply(stscale1, 1, mean) - tscale1
sobias1 <- apply(soscale1, 1, mean) - oscale1
ctbias1 <- apply(ctscale1, 1, mean) - tscale1
cobias1 <- apply(coscale1, 1, mean) - oscale1

summary (stbias1)
plot(main="Bias for for Four Rasch Equating Methods (v1=16)", c(0, 81), c(0, .5), type = "n", xlab="item")
lines(0:80, stbias1, col = 2)
lines(0:80, sobias1, col = 3)
lines(0:80, ctbias1, col = 4)
lines(0:80, cobias1, col = 5)
legend("topright",legend=c("stbias1", "sobias1","ctbias1","cobias1"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)
meanbias1 <- lapply(list(stbias1,  sobias1, ctbias1, cobias1), function(x)
  mean(x, na.rm = T))
# Ploting Bias_nvi1-nvi2)TSE vs. OSE; Separate vs. Concurrent
#nvi2=16
stbias2 <- apply(stscale2, 1, mean) - tscale2
sobias2 <- apply(soscale2, 1, mean) - oscale2
ctbias2 <- apply(ctscale2, 1, mean) - tscale2
cobias2 <- apply(coscale2, 1, mean)- oscale2

plot(main="Bias for for Four Rasch Equating Methods (v2=20)", c(0, 81), c(0, 0.5), type = "n", xlab="item")

lines(0:80, stbias2, col = 2)
lines(0:80, sobias2, col = 3)
lines(0:80, ctbias2, col = 4)
lines(0:80, cobias2, col = 5)

legend("topright",legend=c("stbias2", "sobias2","ctbias2","cobias2"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)

meanbias2 <- lapply(list(stbias2,  sobias2, ctbias2, cobias2), function(x)
  mean(x, na.rm = T))
#nvi3=24
stbias3 <- apply(stscale3, 1, mean) - tscale3
sobias3 <- apply(soscale3, 1, mean) - oscale3
ctbias3 <- apply(ctscale3, 1, mean) - tscale3
cobias3 <- apply(coscale3, 1, mean) - oscale3


plot(main="Bias for for Four Rasch Equating Methods (v3=24)", c(0, 81), c(0, .5), type = "n", xlab="item")
lines(0:80, stbias3, col = 2)
lines(0:80, sobias3, col = 3)
lines(0:80, ctbias3, col = 4)
lines(0:80, cobias3, col = 5)
legend("topright",legend=c("stbias3", "sobias3","ctbias3","cobias3"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)
meanbias3 <- lapply(list(stbias3,  sobias3, ctbias3, cobias3), function(x)
  mean(x, na.rm = T))

# Ploting RMSE_nvi1-nvi2)TSE vs. OSE; Separate vs. Concurrent
#nvi1=16

strmse1 <- sqrt(stbias1^2 + stsee1^2)
sormse1 <- sqrt(sobias1^2 + sosee1^2)
ctrmse1 <- sqrt(ctbias1^2 + ctsee1^2)
cormse1 <- sqrt(cobias1^2 + cosee1^2)

plot(main="RMSE for for Four Rasch Equating Methods (v1=16)", c(0, 81), c(0, 3), type = "n", xlab="item")
lines(0:80, strmse1, col = 2)
lines(0:80, sormse1, col = 3)
lines(0:80, ctrmse1, col = 4)
lines(0:80, cormse1, col = 5)

legend("topright",legend=c("strmse1","sormse1","ctrmse1","cormse1"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)

meanrmse1 <- lapply(list(strmse1,sormse1, ctrmse1, sormse1), function(x)
  mean(x, na.rm = T))
#plotting RMSE_nvi1-nvi5_OSE
#nvi2=20
strmse2 <- sqrt(stbias2^2 + stsee2^2)
sormse2 <- sqrt(sobias2^2 + sosee2^2)
ctrmse2 <- sqrt(ctbias2^2 + ctsee2^2)
cormse2 <- sqrt(cobias2^2 + cosee2^2)


plot(main="RMSE for for Four Rasch Equating Methods (v2=20)", c(0, 81), c(0, 3), type = "n", xlab="item")
lines(0:80, strmse2, col = 2)
lines(0:80, sormse2, col = 3)
lines(0:80, ctrmse2, col = 4)
lines(0:80, cormse2, col = 5)

legend("topright",legend=c("strmse2","sormse2","ctrmse2", "cormse2"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)
meanrmse2 <- lapply(list(strmse2,sormse2, ctrmse2, sormse2), function(x)
  mean(x, na.rm = T))
#nvi3=24
strmse3 <- sqrt(stbias3^2 + stsee3^2)
sormse3 <- sqrt(sobias3^2 + sosee3^2)
ctrmse3 <- sqrt(ctbias3^2 + ctsee3^2)
cormse3 <- sqrt(cobias3^2 + cosee3^2)

plot(main="RMSE for for Four Rasch Equating Methods (v3=24)", c(0, 81), c(0, 3), type = "n", xlab="item")
lines(0:80, strmse3, col = 2)
lines(0:80, sormse3, col = 3)
lines(0:80, ctrmse3, col = 4)
lines(0:80, cormse3, col = 5)

legend("topright",legend=c("strmse3","sormse3","ctrmse3","cormse3"),col=2:5, bty = "n",horiz=F,cex=.75,lty=c(1,1),pt.cex=1)
meanrmse3 <- lapply(list(strmse3,sormse3, ctrmse3, sormse3), function(x)
  mean(x, na.rm = T))