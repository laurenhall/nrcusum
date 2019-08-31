###########################################################
# NRCUSUM Correction Code - Expanded Version              #
# Analyses are separated for easier exploration of 		  #
# individual sets of tests, if desired. 				  #
########################################################### 
source("nrcusum_functions.R")

# Set the values for Lambda[0], Lambda[1], and the true Lambda[a].
lambda0 = 5
lambda1 = 7.5
lambdaa = 10

# Simulate 100 data sets using the above values. 
# Each data set consists of 125 time periods, with a 25 period outbreak in the middle. 
outbreak = rep(c(0,1,0), times=c(50,25,50)) 

# Simulate 100 data sets of 125 time periods.
nsets = 100
set.seed(6)
simpois = lapply(seq_len(nsets), function(x) rpois(length(outbreak), ifelse(outbreak, lambdaa, lambda0)))

###########################################################
# NRCUSUM Correction Tests                                #
# These apply the given correction to the simulated data. #
###########################################################

# (1) No correction - Lambda[1] to detect = 7.5
set.seed(1)
out_n = pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, correction="n"))
dec_n = sapply(out_n, getElement, "decision")
pv_n = sapply(out_n, getElement, "pvalue")

# (2) Correction simulating at the known Lambda[a] value
set.seed(6)
out_a = pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, lambdasim=lambdaa, correction="s"))
dec_a = sapply(out_a, getElement, "decision")
pv_a = sapply(out_a, getElement, "pvalue")

# (3) Correction simulating at Lambda[1]
set.seed(4)
out_s = pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, correction="s"))
dec_s = sapply(out_s, getElement, "decision")
pv_s = sapply(out_s, getElement, "pvalue")

# (4) Correction where Lambda is estimated from the mean of outbreak time periods
set.seed(1)
out_e = pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, correction="e"))
dec_e = sapply(out_e, getElement, "decision")
pv_e = sapply(out_e, getElement, "pvalue")

# (5) Correction with bootstrapped samples
set.seed(3)
out_b = pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, correction="b"))
dec_b = sapply(out_b, getElement, "decision")
pv_b = sapply(out_b, getElement, "pvalue")

#######################################################################
# CUSUM Correction Results                                            #
# This generates a text file with the mean alarm rate before, during, #
# and after the simulated outbreak for each method.                   #
#######################################################################

sink(file="nrcusum_correction_test.txt")
cat("Results of simulated tests - The numbers below represent the following: \n 1: Pre-outbreak false alarm rate \n 2: Power during outbreak \n 3: Post-outbreak false alarm rate \n \n")
cat("No Correction \n")
mean(dec_n[1:50,]) # Pre-outbreak false positives
mean(dec_n[51:75,]) # power during outbreak
mean(dec_n[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Select Lambda - Lambda[a] \n")
mean(dec_a[1:50,]) # Pre-outbreak false positives
mean(dec_a[51:75,]) # power during outbreak
mean(dec_a[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Select Lambda - Lambda[1] \n")
mean(dec_s[1:50,]) # Pre-outbreak false positives
mean(dec_s[51:75,]) # power during outbreak
mean(dec_s[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Estimate Lambda \n")
mean(dec_e[1:50,]) # Pre-outbreak false positives
mean(dec_e[51:75,]) # power during outbreak
mean(dec_e[76:125,]) # post outbreak false alarm rate
cat("\n")
cat("Bootstrap \n")
mean(dec_b[1:50,]) # Pre-outbreak false positives
mean(dec_b[51:75,]) # power during outbreak
mean(dec_b[76:125,]) # post outbreak false alarm rate
cat("\n")
sink()

###############################################
# AUC - Calculate the AUC-ROC for the         #
# trade-off between TPR and post-outbreak FPR #
###############################################
index = 51:125
truth_auc = outbreak[index]
# Remove the pre-outbreak time periods from the p-value data 
pvn_auc = pv_n[index,]
pva_auc = pv_a[index,]
pvs_auc = pv_s[index,]
pve_auc = pv_e[index,]
pvb_auc = pv_b[index,]

# Compute the AUC for each simulated data set 
auc_n = apply(pvn_auc, 2, function(x) auc(roc(truth_auc, x, quiet=TRUE)))
auc_a = apply(pva_auc, 2, function(x) auc(roc(truth_auc, x, quiet=TRUE)))
auc_s = apply(pvs_auc, 2, function(x) auc(roc(truth_auc, x, quiet=TRUE)))
auc_e = apply(pve_auc, 2, function(x) auc(roc(truth_auc, x, quiet=TRUE)))
auc_b = apply(pvb_auc, 2, function(x) auc(roc(truth_auc, x, quiet=TRUE)))

# Compute the mean AUC 
mean(auc_n); mean(auc_a); mean(auc_s); mean(auc_e); mean(auc_b)

##########################################
# Figures - Code to generate the figures #
# used in the body of the manuscript.    #
##########################################

##########################################
# Figure 1 - NRCUSUM Problem Example     #
# This sample data was generated with    #
# Lambda[0] = 4 and Lambda[1] = 6.       #
##########################################
load("cusum_example.rda")

runs = cut(x, c(0, which(stat > h), 100))
spx = split(x, runs)
sps = split(stat, runs)

# Generate the figure 
setEPS()
postscript(file="Fig1.eps", height=5, width=7, bg="white")
plot(x, statnr, lwd=2, type="n", yaxs="i", ylim=c(-1, 35), xlab="Time", ylab="CUSUM Value")
rect(31, -1, 70, 130, col="#eeeeee", xpd=0.25)
for(i in seq_along(spx)){
	lines(spx[[i]], sps[[i]], lwd=2)
}
lines(x, statnr, lwd=2, lty=2)
abline(h = 5, lty=2, lwd=2, col="darkgrey")
box(which="plot")
axis(4, at=5, labels=expression(bolditalic("h")), cex.axis = 1.5, las=2)
text(x = c(14, 50, 86), y = 33, labels=c("In Control", "Out of Control", "In Control"), cex=1.3)
dev.off()

################################
# Figure 2-3 : Sample Data Set #
################################
# Select a random simulated data set to generate figures for
set.seed(1)
pick = sample(1:100,1) #27

# Figure 2 - Cases over time 
setEPS()
postscript(file="Fig2.eps", bg="white", width=6, height=5)
out = simpois[[pick]]
plot(1:125, out, type="b", pch=19, xlab="Time", ylab="Cases")
dev.off()

# Figure 3-A, No Correction 
setEPS()
postscript(file="Fig3a.eps", width=5, height=5)
out = out_nc[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(a) No Correction:   ", lambda[1], " to Detect = 7.5,   ", italic(k), " = 6.17")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#eeeeee", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[out$decision==1],
       out$stat[out$decision==1], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4),
       lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-B, Lambda[a]
setEPS()
postscript(file="Fig3b.eps", width=5, height=5)
out = out_a[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(b) Simulation Correction with ", lambda[a], " = 10")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#eeeeee", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[out$decision==1],
       out$stat[out$decision==1], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-C, Lambda[1]
setEPS()
postscript(file="Fig3c.eps", width=5, height=5)
out = out_s[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(c) Simulation Correction with ", lambda[a], " = ", lambda[1])),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#eeeeee", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[out$decision==1],
       out$stat[out$decision==1], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

#Figure 3-D, Lambda estimated from the data 
setEPS()
postscript(file="Fig3d.eps", width=5, height=5)
out = out_e[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(d) Simulation Correction with ", lambda[a], " = ", hat(lambda)[m])),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#eeeeee", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[out$decision==1],
       out$stat[out$decision==1], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

# Figure 3-E, Bootstrap
setEPS()
postscript(file="Fig3e.eps", width=5, height=5)
out = out_b[[pick]]
plot(1:125, out$stat, xlab="Time", ylab="CUSUM",
     main=expression(paste("(e) Simulation Correction with Bootstrap")),
     type="n", yaxs="i", ylim=c(-1, 130))
rect(51, -1, 75, 130, col="#eeeeee", xpd=0.25, border=NA)
lines(1:125, out$stat, lwd=1)
abline(v = c(51,75), lwd=1, col="darkgrey", lty=2)
points((1:125)[out$decision==1],
       out$stat[out$decision==1], col="red", pch=4)
box(which="plot")
legend("topleft", col=c("grey","red"), pch=c(15,4), lwd=c(1,1), legend=c("Outbreak", "Alarm"), lty=NA, bg="white")
dev.off()

####################################
# Figure 4 - ROC curves for Set 27 #
####################################
# NOTE: The version of Figure 4 which appears in the manuscript incorrectly shows the ROC curves for data set number 75 rather than data set 27. Code to generate both the original figure and the correct figure appear below. 

# Original Figure - Data set 75
setEPS()
postscript(file="Fig4.eps", width=5, height=5)
plot(roc(truth_auc, pv_a[51:125, 75], quiet=TRUE), legacy.axes=T, ylab="True Positive Rate", xlab = "False Positive Rate", col="#f4a582", lty=4, identity.lty=2)
plot(roc(truth_auc, pv_b[51:125, 75], quiet=TRUE), col="#92c5de", lty=5, add=TRUE)
plot(roc(truth_auc, pv_e[51:125, 75], quiet=TRUE), col="#0571b0", lty=3, type="o", pch=20, cex=0.5, add=TRUE)
plot(roc(truth_auc, pv_s[51:125, 75], quiet=TRUE), col="#ca0020", lty=3, add=TRUE)
plot(roc(truth_auc, pv_n[51:125, 75], quiet=TRUE), col="#777777", add=TRUE)
legend("bottomright", col=c("#777777","#f4a582","#ca0020","#0571b0","#92c5de"), lty=c(1,4,3,3,5), lwd=2, pch=c(NA,NA,NA,20,NA), pt.cex=0.8, legend=c("Uncorrected", expression(lambda[a]), expression(lambda[1]), expression(hat(lambda)[m]), "Bootstrap"))
dev.off()

# Correct Figure - Data set 27
setEPS()
postscript(file="Fig4_Corrected.eps", width=5, height=5)
plot(roc(truth_auc, pv_a[51:125, pick], quiet=TRUE), legacy.axes=T, ylab="True Positive Rate", xlab = "False Positive Rate", col="#f4a582", lty=4, identity.lty=2)
plot(roc(truth_auc, pv_b[51:125, pick], quiet=TRUE), col="#92c5de", lty=5, add=TRUE)
plot(roc(truth_auc, pv_e[51:125, pick], quiet=TRUE), col="#0571b0", lty=3, type="o", pch=20, cex=0.5, add=TRUE)
plot(roc(truth_auc, pv_s[51:125, pick], quiet=TRUE), col="#ca0020", lty=3, add=TRUE)
plot(roc(truth_auc, pv_n[51:125, pick], quiet=TRUE), col="#777777", add=TRUE)
legend("bottomright", col=c("#777777","#f4a582","#ca0020","#0571b0","#92c5de"), lty=c(1,4,3,3,5), lwd=2, pch=c(NA,NA,NA,20,NA), pt.cex=0.8, legend=c("Uncorrected", expression(lambda[a]), expression(lambda[1]), expression(hat(lambda)[m]), "Bootstrap"))
dev.off()

##################################################
# Figure 5 - average AUC for varying Lambda[sim] #
##################################################
# Compute the average AUC with different selected levels of Lambda[sim]
# Please be aware that this may take considerable time to run. 
l1grid = seq(5, 15, 0.1)
outlg = pbsapply(l1grid, function(L){
  auc = sapply(simpois, function(x){
    tmp = csm_corrected(x, lambda0, lambda1, lambdasim=L, correction="s")$pvalue
    auc(roc(truth_auc, tmp[51:125],quiet=TRUE))})
    return(mean(auc))})

setEPS()
postscript(file="Fig5.eps", width=5, height=5)
plot(l1grid, outlg, pch=16, xlab=expression(lambda[1]), ylab="AUC")
dev.off()
