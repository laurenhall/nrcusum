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

# Plot the total cases over time to isolate the outbreak we want to detect
#pdf(file="snewportcases.pdf", width=7, height=5)
setEPS()
postscript(file="snewportcases.eps", width=7, height=5, bg="white")
plot(apply(dat, 1, sum), type="l", xlab="Time (in Weeks)", 
     ylab="Total Cases", main="")
dev.off()
# The outbreak we want to detect is the large spike from about t=408 to t=410.
# Region 12 was the only region not to report any cases of S. Newport between weeks 400 and 415, and so can be excluded.
dat = dat[,-12] 
# Estimate the in-control mean for each region using three years of data. 
lambda = sapply(1:15, function(x) mean(dat[1:156,x])) 
# Run the CUSUM on 15 states that recorded cases during the outbreak. 
testn = pblapply(1:15, function(x) {set.seed(x); nrc.csave(dat[,x], lambda[x], lambda[x]*1.5, correction="n", sig=0.05, nsim = 1999)})
tests = pblapply(1:15, function(x) {set.seed(x); nrc.csave(dat[,x], lambda[x], lambda[x]*1.5, correction="s", sig=0.05, nsim = 1999)})
teste = pblapply(1:15, function(x) {set.seed(x); nrc.csave(dat[,x], lambda[x], lambda[x]*1.5, correction="e", sig=0.05, nsim = 1999)})
testb = pblapply(1:15, function(x) {set.seed(x); nrc.csave(dat[,x], lambda[x], lambda[x]*1.5, correction="b", sig=0.05, nsim = 1999)})

# Time of first alarm after outbreak
firstalarm = matrix(NA, nrow=15, ncol=4); rownames(firstalarm) = colnames(dat); colnames(firstalarm) = c("No Correction","Lambda to Detect", "Estimation", "Bootstrap")
for(i in 1:15){
  firstalarm[i,1] = (408:528)[which(testn[[i]]$decision[408:528] == 1)][1]
  firstalarm[i,2] = (408:528)[which(tests[[i]]$decision[408:528] == 1)][1]
  firstalarm[i,3] = (408:528)[which(teste[[i]]$decision[408:528] == 1)][1]
  firstalarm[i,4] = (408:528)[which(testb[[i]]$decision[408:528] == 1)][1]
}

firstalarm

# calculate number of alarms after detection: total alarms from 412
alarms = matrix(NA, nrow=15, ncol=4); rownames(alarms) = colnames(dat); colnames(alarms) = c("No Correction","Lambda to Detect", "Estimation", "Bootstrap") 
for(i in 1:15){
  alarms[i,1] = max(0,sum(which(testn[[i]]$decision == 1) > 412))
  alarms[i,2] = max(0,sum(which(tests[[i]]$decision == 1) > 412))
  alarms[i,3] = max(0,sum(which(teste[[i]]$decision == 1) > 412))
  alarms[i,4] = max(0,sum(which(testb[[i]]$decision == 1) > 412))
}

alarms
