##########################################
# NRCUSUM Correction Case Study          #
# Cases of Salmonella Newport in Germany #
##########################################
source("nrcusum_functions.R")
# Load in the data, freely available through the Surveillance package
library(surveillance)
data(salmNewport)
# Isolate the observed counts
dat = observed(salmNewport); dim(dat)

# Figure 6 - Total S. Newport cases over time. 
setEPS()
postscript(file="Fig6.eps", width=7, height=5, bg="white")
plot(apply(dat, 1, sum), type="l", xlab="Time (in Weeks)",
     ylab="Total Cases", main="")
dev.off()

# The outbreak we want to detect is the large spike from about t=408 to t=410.
# Region 12 was the only region not to report any cases of S. Newport between weeks 400 and 415, and so is excluded from the analysis.
dat = dat[,-12]

# Estimate the in-control mean for each region using three years of data.
lambdas = sapply(1:15, function(x) mean(dat[1:156,x]))

# Run the CUSUM on 15 states that recorded cases during the outbreak.
test_n = pblapply(1:15, function(x) {set.seed(x); csm_corrected(dat[,x], lambdas[x], lambdas[x]*1.5, correction="n", sig=0.05, nsim = 1999)})
test_s = pblapply(1:15, function(x) {set.seed(x); csm_corrected(dat[,x], lambdas[x], lambdas[x]*1.5, correction="s", sig=0.05, nsim = 1999)})
test_e = pblapply(1:15, function(x) {set.seed(x); csm_corrected(dat[,x], lambdas[x], lambdas[x]*1.5, correction="e", sig=0.05, nsim = 1999)})
test_b = pblapply(1:15, function(x) {set.seed(x); csm_corrected(dat[,x], lambdas[x], lambdas[x]*1.5, correction="b", sig=0.05, nsim = 1999)})


firstalarm = alarms = matrix(NA, nrow=15, ncol=4, dimnames = list(colnames(dat), c("No Correction","Lambda[1]", "Estimation", "Bootstrap")))
# Time of first alarm after outbreak
for(i in 1:15){
  firstalarm[i,1] = (408:528)[test_n[[i]]$decision[408:528] == 1][1]
  firstalarm[i,2] = (408:528)[test_s[[i]]$decision[408:528] == 1][1]
  firstalarm[i,3] = (408:528)[test_e[[i]]$decision[408:528] == 1][1]
  firstalarm[i,4] = (408:528)[test_b[[i]]$decision[408:528] == 1][1]
}

print(firstalarm)

# Calculate the number of alarms that occured between week 412 and the end of the study period.  
for(i in 1:15){
  alarms[i,1] = sum(which(test_n[[i]]$decision == 1) > 412)
  alarms[i,2] = sum(which(test_s[[i]]$decision == 1) > 412)
  alarms[i,3] = sum(which(test_e[[i]]$decision == 1) > 412)
  alarms[i,4] = sum(which(test_b[[i]]$decision == 1) > 412)
}

print(alarms)
