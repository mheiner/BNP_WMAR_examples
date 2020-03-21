# rm(list=ls())
library("coda")

id = commandArgs(trailingOnly = TRUE)



### user input

simregexp = paste0(".*", id, ".*\\.rda")
dctory = "./postsim"
allfiles = list.files(dctory)
(selfiles = allfiles[grep(simregexp, allfiles)])
(nfiles = length(selfiles))

### end user input

for ( i in 1:nfiles ) {
  load( paste0(dctory, "/", selfiles[i]) )

  nuse = 2000
  whichiter = floor(seq(1, nsim, length=nuse))

  pdf(paste0("plots/trace/", mesg, ".pdf"), height=4, width=9)

  traceplot(as.mcmc(sims_llik[whichiter]), main="log-likelihood")

  traceplot(as.mcmc(sims_noccup[whichiter]), main="occupied clusters")

  cM_Scounts = colMeans(sims_Scounts[whichiter,])
  (maxclust_single = which.max(cM_Scounts))
  traceplot(as.mcmc(sims_Scounts[whichiter,maxclust_single]), main="size of most occupied single index")
  maxclust_indx = apply(sims_Scounts[whichiter,], 1, which.max)
  traceplot(as.mcmc(sims_Scounts[cbind(whichiter,maxclust_indx)]), main="size of most occupied cluster")


  if (gam_type == "local") {

    for (k in 1:K) {
        traceplot(as.mcmc(sims_gam_occ[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main=paste0("Inclusion of lag ", k, "\n p_global_hypothesis: ", round(mean(sims_gam_glob[whichiter,k]),3)))
        axis(1)
        axis(2, at=c(0,1))
        traceplot(as.mcmc(sims_gam_w[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main="", add=TRUE, col="red")
        traceplot(as.mcmc(sims_gam_glob[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main="", add=TRUE, col="black", lty=2)
        traceplot(as.mcmc(sims_gam_alt_occup[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main="", add=TRUE, col="blue", lty=3)

        legend("topright", bty="n", lty=c(1,1,2,3), col=c("black", "red", "black", "blue"), legend=c("occupied", "weights", "global hypoth", "slope"))
    }

    for (k in 1:K) {
        traceplot(as.mcmc(sims_pigam[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main=paste0("pi_gamma, lag ", k))
        axis(1)
        axis(2, at=c(0,1))
    }

  } else if (gam_type == "global") {

    if (typeof(sims_fc_on) == "character") {
      for (k in 1:K) {
        traceplot(as.mcmc(sims_gam[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main=paste0("Inclusion of lag ", k, "\n p: ", round(mean(sims_gam[whichiter,k]),3)))
        axis(1)
        axis(2, at=c(0,1))
        traceplot(as.mcmc(sims_gam_alt_occup[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main="", add=TRUE, col="blue", lty=3)
        legend("topright", bty="n", lty=c(1,3), col=c("black", "blue"), legend=c("standard", "slope"))
      }
    } else {
      for (k in 1:K) {
        traceplot(as.mcmc(sims_gam[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main=paste0("Inclusion of lag ", k, "\n p: ", round(mean(sims_gam[whichiter,k]),3), ", p_fc: ", round(mean(sims_fc_on[whichiter,k]),3)))
        axis(1)
        axis(2, at=c(0,1))
        traceplot(as.mcmc(sims_gam_alt_occup[whichiter,k]), ylim=c(0,1), axes=FALSE,
                main="", add=TRUE, col="blue", lty=3)
        legend("topright", bty="n", lty=c(1,3), col=c("black", "blue"), legend=c("standard", "slope"))
      }
    }  
  }


  traceplot(as.mcmc(sims_w[whichiter,1]), main="w[1]")
  traceplot(as.mcmc(sims_w[whichiter,H]), main="w[H]")
  traceplot(as.mcmc(sims_w[whichiter,maxclust_single]), main=bquote("w"*.(maxclust_single)~": most occupied single index"))
  traceplot(as.mcmc(sims_w[cbind(whichiter,maxclust_indx)]), main="w : most occupied cluster")

  traceplot(as.mcmc(sims_alpha[whichiter]), main=expression(alpha))

  traceplot(as.mcmc(sims_mu[whichiter,maxclust_single]), main=bquote(mu[y]~": most occupied single index"))
  traceplot(as.mcmc(sims_mu[cbind(whichiter,maxclust_indx)]), main=bquote(mu[y]~": most occupied cluster"))

  traceplot(as.mcmc(sims_int[whichiter,maxclust_single]), main=bquote("intercept: most occupied single index"))
  traceplot(as.mcmc(sims_int[cbind(whichiter,maxclust_indx)]), main=bquote("intercept: most occupied cluster"))

  for (k in 1:K) {
    traceplot(as.mcmc(sims_beta[whichiter,k,maxclust_single]), main=bquote(beta[y]~"lag"~.(k)~": most occupied single index"))
    abline(h=0, lty=2)
    traceplot(as.mcmc(sims_beta[cbind(whichiter,rep(k,nuse),maxclust_indx)]), main=bquote(beta[y]~"lag"~.(k)~": most occupied cluster"))
    abline(h=0, lty=2)
  }

  traceplot(as.mcmc(sims_sig2[whichiter,maxclust_single]), main=bquote(sigma^2~": most occupied single index"))
  traceplot(as.mcmc(sims_sig2[cbind(whichiter,maxclust_indx)]), main=bquote(sigma^2~": most occupied cluster"))


  dev.off()

  cat(i, " of ", nfiles, "\r")
}
