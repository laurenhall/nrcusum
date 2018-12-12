#####################################
# NRCUSUM Correction Code           #
#####################################
source("nrcusum_functions.R")

# We'll start by selecting values for Lambda[0], Lambda[a], and the Lambda to detect (Lambda[1]).
# Lambda[a] is represented by "lambdat", the (t)rue Lambda value.
# We'll simulate 100 single CUSUM streams with these values.
lambda0 = 5
lambda1 = 7.5
lambdat = 10
outbreak = c(rep(0,50),rep(1,25),rep(0,50)) # outbreak for 25 periods in the middle

# Simulate 100 data sets of 125 time periods.
nsets = 100
set.seed(6)
simpois = lapply(seq_len(nsets), function(x) rpois(length(outbreak), ifelse(outbreak, lambdat, lambda0)))

##########################################################
# NRCUSUM Correction Tests                                #
# These apply the given correction to the simulated data. #
###########################################################

# No correction - Lambda[1] to detect = 7.5
set.seed(12)
outnc = pblapply(simpois, FUN = function(x) nrc.csave(x, lambda0, lambda1, correction="n"))
decnc = sapply(outnc, getElement, "decision")
pvnc = sapply(outnc, getElement, "pvalue")
# Correction with Lambda[a] specified as the known true value
set.seed(6)
outs = pblapply(simpois, FUN = function(x) nrc.csave(x, lambda0, lambda1, lambdat=lambdat, correction="s"))
decs = sapply(outs, getElement, "decision")
pvs = sapply(outs, getElement, "pvalue")
# Correction with Lambda[a] specified to be Lambda[1]
set.seed(4)
outs2 = pblapply(simpois, FUN = function(x) nrc.csave(x, lambda0, lambda1, correction="s"))
decs2 = sapply(outs2, getElement, "decision")
pvs2 = sapply(outs2, getElement, "pvalue")
# Correction where Lambda is estimated from the mean of outbreak time periods
set.seed(1)
oute = pblapply(simpois, FUN = function(x) nrc.csave(x, lambda0, lambda1, correction="e"))
dece = sapply(oute, getElement, "decision")
pve = sapply(oute, getElement, "pvalue")
# Correction with bootstrapped samples
set.seed(3)
outb = pblapply(simpois, FUN = function(x) nrc.csave(x, lambda0, lambda1, correction="b"))
decb = sapply(outb, getElement, "decision")
pvb = sapply(outb, getElement, "pvalue")

#######################################################################
# CUSUM Correction Results                                            #
# This generates a text file with the mean alarm rate before, during, #
# and after the simulated outbreak for each method.                   #
#######################################################################
sink(file="nrcusum_correction_test.txt")
cat("Results of simulated tests - The numbers below represent the following: \n 1: Pre-outbreak false alarm rate \n 2: Power during outbreak \n 3: Post-outbreak false alarm rate \n \n")
cat("No Correction \n")
mean(decnc[1:50,]) # Pre-outbreak false positives
mean(decnc[51:75,]) # power during outbreak
mean(decnc[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Select Lambda - Lambda[a] \n")
mean(decs[1:50,]) # Pre-outbreak false positives
mean(decs[51:75,]) # power during outbreak
mean(decs[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Select Lambda - Lambda[1] \n")
mean(decs2[1:50,]) # Pre-outbreak false positives
mean(decs2[51:75,]) # power during outbreak
mean(decs2[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Estimate Lambda \n")
mean(dece[1:50,]) # Pre-outbreak false positives
mean(dece[51:75,]) # power during outbreak
mean(dece[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Bootstrap \n")
mean(decb[1:50,]) # Pre-outbreak false positives
mean(decb[51:75,]) # power during outbreak
mean(decb[76:125,]) # post outbreak false alarm rate
cat("\n")
sink()

###############################################
# AUC - Calculate the AUC-ROC for the         #
# trade-off between TPR and post-outbreak FPR #
###############################################
# We start with the AUC
library(pROC)
# Remove the pre-outbreak time periods, which are well controlled
pvnc_auc = pvnc[-(1:50),]
pvs_auc = pvs[-(1:50),]
pvs2_auc = pvs2[-(1:50),]
pve_auc = pve[-(1:50),]
pvb_auc = pvb[-(1:50),]

truth_auc = outbreak[-(1:50)]

auc_nc = apply(pvnc_auc, 2, function(x) auc(roc(truth_auc, x)))
auc_s = apply(pvs_auc, 2, function(x) auc(roc(truth_auc, x)))
auc_s2 = apply(pvs2_auc, 2, function(x) auc(roc(truth_auc, x)))
auc_e = apply(pve_auc, 2, function(x) auc(roc(truth_auc, x)))
auc_b = apply(pvb_auc, 2, function(x) auc(roc(truth_auc, x)))

mean(auc_nc);mean(auc_s);mean(auc_s2);mean(auc_e);mean(auc_b)

##########################################
# Figures - Code to generate the figures #
# used in the body of the manuscript.    #
##########################################
# setwd("C:/Users/Lauren/Dropbox/NRCUSUM/paper_v1/Figures/")
# setwd("~/Dropbox/NRCUSUM/paper_v1/Figures/")
# NRCUSUM Problem Example
x = 1:100
set.seed(12)
y = c(rpois(30, 4), rpois(40, 6), rpois(30, 4)) # Ran 4 times after setting seed to get easy-to-read shape
k = (6 - 4)/(log(6) - log(4))
stat = csm.s(y, k)
pdf(file="nrcusumexample.pdf", height=5, width=10, bg="white")
#setEPS()
#postscript(file="nrcusumexample.eps", height=5, width=7, bg="white")
plot(x, stat, lwd=2, type="n", yaxs="i", ylim=c(-1, 35), xlab="Time", ylab="CUSUM Value")
rect(31, -1, 70, 130, col="#EEEEEE", xpd=0.25)
lines(x, stat, lwd=2)
abline(h = 5, lty=2, lwd=2, col="darkgrey")
box(which="plot")
axis(4, at=5, labels=expression(bolditalic("h")), cex.axis = 1.5, las=2)
text(x = c(14, 50, 86), y = 33, labels=c("In Control", "Out of Control", "In Control"), cex=1.3)
dev.off()

# Select a random simulated data set to generate figures for
set.seed(1)
pick = sample(1:100,1) #27

# Cases over time

setEPS()
postscript(file="simcases.eps", bg="white", width=6, height=5)
out = simpois[[pick]]
plot(1:125, out, type="b", pch=19, xlab="Time", ylab="Cases")
dev.off()

############
# Figure 3 #
############

# Figure 3-A, uncorrected
#pdf(file="lambdaa.pdf", width=10, height=5)
setEPS()
postscript(file="Fig3a.eps", width=5, height=5)
out = outnc[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(a) No Correction:   ", lambda[1], " to Detect = 7.5,   ", italic(k), " = 6.17")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#EEEEEE", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[which(out$decision==1)],
       out$stat[which(out$decision==1)], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4),
       lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-B, known lambda
setEPS()
postscript(file="Fig3b.eps", width=5, height=5)
out = outs[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(b) Simulation Correction with ", lambda[a], " = 10")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#EEEEEE", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[which(out$decision==1)],
       out$stat[which(out$decision==1)], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-C, lambda 1
#pdf(file="lambda1.pdf", width=10, height=5)
setEPS()
postscript(file="Fig3c.eps", width=5, height=5)
out = outs2[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(c) Simulation Correction with ", lambda[a], " = ", lambda[1])),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#EEEEEE", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[which(out$decision==1)],
       out$stat[which(out$decision==1)], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

#Figure 3-D, estimated lambda
#pdf(file="lambdam.pdf", width=10, height=5)
setEPS()
postscript(file="Fig3d.eps", width=5, height=5)
out = oute[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(d) Simulation Correction with ", lambda[a], " = ", hat(lambda)[m])),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#EEEEEE", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[which(out$decision==1)],
       out$stat[which(out$decision==1)], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-E, bootstrap
#pdf(file="lambdab.pdf", width=10, height=5)
setEPS()
postscript(file="Fig3e.eps", width=5, height=5)
out = outb[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(e) Simulation Correction with Bootstrap")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#EEEEEE", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[which(out$decision==1)],
       out$stat[which(out$decision==1)], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

####################################
# Figure 4 - AUC curves for Set 27 #
####################################

setEPS()
postscript(file="Fig4.eps", width=5, height=5)
plot(roc(truth_auc, pvs[51:125,75]), xaxs="i",legacy.axes=T, ylab="True Positive Rate",
     xlab = "False Positive Rate", col="#f4a582", lty=4, identity.lty=2)
plot(roc(truth_auc, pvb[51:125,75]), col="#92c5de", lty=5, add=TRUE)
plot(roc(truth_auc, pve[51:125,75]), col="#0571b0", lty=3, type="o", pch=20, cex=0.5, add=TRUE)
plot(roc(truth_auc, pvs2[51:125,75]), col="#ca0020", lty=3, add=TRUE)
plot(roc(truth_auc, pvnc[51:125,75]), col="#777777", add=TRUE)
legend("bottomright", col=c("#777777","#f4a582","#ca0020","#0571b0","#92c5de"), lty=c(1,4,3,3,5), lwd=2,
       pch=c(NA,NA,NA,20,NA), pt.cex=0.8, legend=c("Uncorrected", expression(lambda[a]),
                expression(lambda[1]),
                expression(hat(lambda)[m]),
                "Bootstrap"))
dev.off()

###############################################
# Figure 5 - average AUC for various lambda 1 #
###############################################
l1grid = seq(5, 15, 0.1)
outlg = pbsapply(l1grid, function(l){
  auc = sapply(simpois, function(x){
    tmp = nrc.csave(x, lambda0, lambda1, lambdat=l, correction="s")$pvalue
    auc(roc(truth_auc, tmp[51:125]))})
    return(mean(auc))})
auc = apply(outlg, 2, mean)
setEPS()
postscript(file="Fig5.eps", width=5, height=5)
plot(l1grid, outlg, pch=16, xlab=expression(lambda[1]), ylab="AUC")
dev.off()
