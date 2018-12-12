# Libraries
library(pbapply)

# A function to simulate under the alternative hypothesis.
pickc = function(cases, dec, lambda1, lambdat, nsim, correction){
if(correction == "s"){out = rpois(nsim, ifelse(is.null(lambdat), lambda1, lambdat))}
if(correction == "e"){out = rpois(nsim, mean(cases[which(dec==1)]))}
if(correction == "b"){
  if(sum(dec) == 1){
    out = rep(cases[which(dec==1)], times=nsim)
  } else {
    out = sample(cases[which(dec==1)], nsim, replace=TRUE)}}
return(out)
  }

# A function to implement the corrected cusum.
# Correction key:
# "n" = no correction, "s" = simulate at lambda1, or a specified lambdat (if different from l1),
# "e" = estimate lambda from data, "b" = bootstrap.
nrc.csave = function(cases, lambda0, lambda1, lambdat=NULL, nsim=499, sig=0.05, start=0, correction=c("n","s","e","b")){
  if(lambda0 == lambda1){k=0
  } else {
  k = (lambda1 - lambda0)/(log(lambda1) - log(lambda0))}
stat = csm.s(cases, k, start)
dec = pval = numeric(length(cases))
# start value - a constant (first time period) or a vector of length nsim (last time period)
c0 = start
for(i in 1:length(cases)){
  # Sim under null
  simi = matrix(rpois(nsim, lambda0))
  # Calculate the cusum for this step
  ci = matrix(apply(simi + c0 - k, 1, function(x) max(0,x)))
  # calculate p-value and make decision
  pval[i] = (sum(as.numeric(ci>= stat[i]))+1)/(nsim+1)
  dec[i] = as.numeric(pval[i] < sig)
  if(dec[i] == 1 & correction!="n"){
    # re-simulate the CUSUMs under H1 if appropriate
    sim1 = matrix(pickc(cases, dec, lambda1, lambdat, nsim, correction))
    c0 = matrix(apply(sim1 + c0 - k, 1, function(x) max(0,x)))
  } else {
    # If no outbreak or no correction, just save the CUSUMs simulated under the null
    c0 = ci
   }
  }
return(list("stat" = stat, "decision" = dec, "pvalue" = pval))}
