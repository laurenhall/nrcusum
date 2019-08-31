# Libraries
library(pbapply)
library(pROC)

# A function to compute the CUSUM for a single stream of data and a given value of k. 
csm = function(cases, k, start=0){
	ct = numeric(length(cases))
	ct0 = start
	for(i in seq_along(cases)){
		ct0 = max(0, ct0 + cases[i] - k)
		ct[i] = ct0
	}
	return(ct)
}

# A function to simulate under the alternative hypothesis.
# Note: the sample() function behaves differently if a single integer is given as the data, so the bootstrap correction requires an extra line of code to function correctly.  
sim_alt = function(cases, dec, lambdasim, correction, nsim){
	if(correction == "s"){out = rpois(nsim, lambdasim)}
	if(correction == "e"){out = rpois(nsim, mean(cases[dec == 1]))}
	if(correction == "b"){
		if(sum(dec) == 1){
			out = rep(cases[dec == 1], times=nsim)
		} else {
			out = sample(cases[dec == 1], nsim, replace=TRUE)}}
	return(matrix(out))
  }

# The primary function to implement the corrected CUSUM. 
# Correction must be one of the following:
# "n" - No correction.
# "s" - Simulate Poisson draws at the specified Lambda[sim]. 
# "e" - Simulate Poisson draws by estimating Lambda from the observed data. 
# "b" - Bootstrap from the observed data. 
csm_corrected = function(cases, lambda0, lambda1, lambdasim=lambda1, start=0, correction="n", sig=0.05, nsim=499){
	dec = pval = numeric(length(cases))
	ct0 = start
	# Compute k using the standard formula for a Poisson process 
	if(lambda0 == lambda1){
		k = 0
	} else {
		k = (lambda1 - lambda0)/(log(lambda1) - log(lambda0))
	}	
	# Compute the CUSUM statistic 
	stat = csm(cases, k, start)
	# Perform the hypothesis test 
	for(i in seq_along(cases)){
		simi = matrix(rpois(nsim, lambda0))
		cti = matrix(apply(simi + ct0 - k, 1, function(x) max(0,x)))
		pval[i] = (sum(cti >= stat[i]) + 1)/(nsim + 1)
		dec[i] = as.numeric(pval[i] < sig)
		# re-simulate the CUSUMs under H1 if appropriate, otherwise save the CUSUMs simulated under the null hypothesis 
		if(dec[i] == 1 & correction %in% c("s", "b", "e")){
			sim1 = sim_alt(cases, dec, lambdasim, correction, nsim)
			ct0 = matrix(apply(sim1 + ct0 - k, 1, function(x) max(0, x)))
		} else {
			ct0 = cti
		}
	}
	return(list("stat" = stat, "decision" = dec, "pvalue" = pval))
}
