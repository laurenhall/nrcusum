#####################################
# NRCUSUM Correction Code           #
#####################################
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

# A total of five tests are performed: 
# (1) No correction - Lambda[1] to detect = 7.5
# (2) Correction simulating at the known Lambda[a] value
# (3) Correction simulating at Lambda[1]
# (4) Correction where Lambda is estimated from the mean of outbreak time periods
# (5) Correction with bootstrapped samples

corrections = c("n","s","s","e","b")
lambdas = rep(lambda1, 5); lambdas[2] = lambdaa
seeds = c(12, 6, 4, 1, 3)
labs = corrections; labs[2] = "a"
index = rep(c(0,1,2), times=c(50,25,50))

# Apply the uncorrected and corrected tests to the simulated data  
test_out = lapply(seq_along(corrections), function(i){
	set.seed(seeds[i])
	pblapply(simpois, FUN = function(x) csm_corrected(x, lambda0, lambda1, lambdasim = lambdas[i], correction=corrections[i]))
	})
names(test_out) = labs 

# Calculate the AUC-ROC for the trade-off between true positive rate and post-outbreak false positive rate 

auc_out = sapply(test_out, function(x){
	tmp = sapply(x, getElement, "pvalue")[index > 0,]
	apply(tmp, 2, function(p) auc(roc(outbreak[index > 0], p, quiet=T)))
})

# Generate the table of test results 
test_results = matrix(nrow=5, ncol=4, dimnames= list(c("No Correction","Lambda[a]","Lambda[1]", "Estimation","Bootstrap"), c("Pre-Outbreak","Outbreak","Post-Outbreak","AUC")))
test_results[,4] = colMeans(auc_out)

for(i in seq_along(test_out)){
	tmp = sapply(test_out[[i]], getElement, "decision")
	test_results[i,-4] = tapply(tmp, matrix(index, 125, 100), mean)
}

print(test_results, digits=3)


##########################################
# Figures - Code to generate the figures #
# used in the body of the manuscript.    #
##########################################

##########################################
# Figure 1 - NRCUSUM Problem Example     #
# This sample data was generated with    #
# Lambda[0] = 4 and Lambda[1] = 6.       #
##########################################
load(example.rda)

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
out = test_out$n[[pick]]
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
out = test_out$a[[pick]]
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
out = test_out$s[[pick]]
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
out = test_out$e[[pick]]
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
out = test_out$b[[pick]]
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

# Order to draw the lines to maximize visibility 
ord = c(2,5,4,3,1)
cols = c("#777777","#f4a582","#ca0020","#0571b0","#92c5de")
ltys = c(1, 4, 3, 3, 5)
types = c("l","l","l","o","l")

# Original Figure - Data set 75
setEPS()
postscript(file="Fig4.eps", width=5, height=5)
plot(roc(outbreak[index > 0], test_out[[2]][[75]]$pvalue[index > 0], quiet=TRUE), ylab="True Positive Rate", xlab="False Positive Rate", col=cols[2], lty=ltys[2], identity.lty = 2, legacy.axes=TRUE)
for(i in ord[-1]){
	plot(roc(outbreak[index > 0], test_out[[i]][[75]]$pvalue[index > 0], quiet=TRUE), col=cols[i], lty=ltys[i], type=types[i], pch=20, add=TRUE)
}
legend("bottomright", col=cols, lty=ltys, lwd=2, pch=c(NA,NA,NA,20,NA), pt.cex=0.8, legend=c("Uncorrected", expression(lambda[a]), expression(lambda[1]), expression(hat(lambda)[m]), "Bootstrap"))
dev.off()

# Correct Figure - Data set 27
setEPS()
postscript(file="Fig4_Corrected.eps", width=5, height=5)
plot(roc(outbreak[index > 0], test_out[[2]][[27]]$pvalue[index > 0], quiet=TRUE), ylab="True Positive Rate", xlab="False Positive Rate", col=cols[2], lty=ltys[2], identity.lty = 2, legacy.axes=TRUE)
for(i in ord[-1]){
	plot(roc(outbreak[index > 0], test_out[[i]][[27]]$pvalue[index > 0], quiet=TRUE), col=cols[i], lty=ltys[i], type=types[i], pch=20, add=TRUE)
}
legend("bottomright", col=cols, lty=ltys, lwd=2, pch=c(NA,NA,NA,20,NA), pt.cex=0.8, legend=c("Uncorrected", expression(lambda[a]), expression(lambda[1]), expression(hat(lambda)[m]), "Bootstrap"))
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
