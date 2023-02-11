# Function that returns all different qs, not only the selected one.
# Used to correlate WS.dist with difference in CV - true error.

knndmW <- function(tpoints, modeldomain = NULL, ppoints = NULL,
                  k = 10, maxp = NULL, linkf = "ward.D2",
                  samplesize = 1000, sampling = "regular"){

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(tpoints)
  distmat <- sf::st_distance(tpoints)
  units(distmat) <- NULL
  Gj <- distclust(distmat, clust)

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(ppoints, tpoints)
  Gij <- apply(Gij, 1, min)
  units(Gij) <- NULL

  # Check if there's is clustering in the data in the first place
  testks <- stats::ks.test(Gj, Gij, alternative = "great")
  if(testks$p.value >= 0.05){

    clust <- sample(rep(1:k, ceiling(nrow(tpoints)/k)), size = nrow(tpoints), replace=F)
    Gjstar <- distclust(distmat, clust)
    stat_final <- twosamples::wass_stat(Gjstar, Gij)
    warning("No evidence of clustering has been found, a random CV assignment is returned")

  }else{

    # Hierarchical clustering
    hc <- stats::hclust(d = stats::as.dist(distmat), method=linkf)

    # Build grid of number of clusters to try - we sample low numbers more intensively
    clustgrid <- data.frame(nk = as.integer(round(exp(seq(log(k), log(nrow(tpoints)),
                                                          length.out = 100)))))
    clustgrid$stat <- NA
    clustgrid <- clustgrid[!duplicated(clustgrid$nk),]
    clustgroups <- list()

    # We test each number of clusters
    for(nk in clustgrid$nk){

      # Cut nk clusters
      clust_nk <- stats::cutree(hc, k=nk)
      tabclust <- as.data.frame(table(clust_nk))
      tabclust <- tabclust[order(tabclust$Freq, decreasing=T),]
      tabclust$clust_k <- NA

      # We don't merge big clusters
      clust_i <- 1
      for(i in 1:nrow(tabclust)){
        if(tabclust$Freq[i] >= nrow(tpoints)/k){
          tabclust$clust_k[i] <- clust_i
          clust_i <- clust_i + 1
        }
      }
      rm("clust_i")
      clust_i <- setdiff(1:k, unique(tabclust$clust_k))
      tabclust$clust_k[is.na(tabclust$clust_k)] <- rep(clust_i, ceiling(nk/length(clust_i)))[1:sum(is.na(tabclust$clust_k))]
      tabclust2 <- data.frame(ID = 1:length(clust_nk), clust_nk = clust_nk)
      tabclust2 <- merge(tabclust2, tabclust, by = "clust_nk")
      tabclust2 <- tabclust2[order(tabclust2$ID),]
      clust_k <- tabclust2$clust_k

      # Compute statistic if not exceeding limit
      if(!any(table(clust_k)/length(clust_k)>maxp)){
        Gjstar_i <- distclust(distmat, clust_k)
        clustgrid$stat[clustgrid$nk==nk] <- twosamples::wass_stat(Gjstar_i, Gij)
        clustgroups[[paste0("nk", nk)]] <- clust_k
      }
    }

    stat_final <- clustgrid[!is.na(clustgrid$stat),"stat"]
    clust <- clustgroups[!is.na(clustgrid$stat)]

  }

  # Output
    res <- list(clusters = clust, stat=stat_final)
    class(res) <- c("knndm", "list")
    res

}

# Helper function: Compute out-of-fold NN distance
distclust <- function(distm, folds){

  alldist <- c()
  for(f in unique(folds)){
    alldist <- c(alldist, apply(distm[f == folds, f != folds, drop=FALSE], 1, min))
  }
  alldist
}
