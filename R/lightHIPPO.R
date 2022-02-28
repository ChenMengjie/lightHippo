cut_off_zscore <- function(num.gene, significance_level = 0.1, method = "fdr"){

  find_cutoff <- function(x){ abs(p.adjust(x, method = method, n = num.gene) - significance_level)}
  pvalue_before_correction <- optim(significance_level/num.gene, find_cutoff, method= "Brent",  lower = 0.1/num.gene, upper = 0.1/num.gene*20)$par
  zscore <- qnorm(pvalue_before_correction, lower.tail=FALSE)
  return(zscore)

}

compute_test_statistic_for_one_gene <- function(gene_mean, zero_proportion, cell.number) {

  expected_pi <- pmin(exp(-gene_mean), 1-1e-10)
  expected_pi <- pmax(1e-10, expected_pi)
  est_se <- sqrt(expected_pi*(1-expected_pi)/cell.number)
  pvalue <- pnorm(zero_proportion, expected_pi, est_se, lower.tail = FALSE)
  return(pvalue)

}

compute_zscore_for_one_gene <- function(gene_mean, zero_proportion, cell.number) {
  expected_pi <- pmin(exp(-gene_mean), 1-1e-10)
  expected_pi <- pmax(1e-10, expected_pi)
  est_se <- sqrt(expected_pi*(1-expected_pi)/cell.number)
  propdiff <- zero_proportion - expected_pi
  zvalue <- propdiff/est_se
  return(zvalue)
}

select_features_light <- function(X, stop_at = 500, Zscore.cutoff = 2) {

  gene.number <- nrow(X)
  cell.number <- ncol(X)

  selected <- NULL
  selected.Zscore <- NULL
  seleced.num <- 0
  permutated <-  sample(1:gene.number)

  for(i in permutated){
    gene_mean <- sum(X[i, ])/cell.number
    zero_proportion <- sum(X[i, ]==0)/cell.number
    if(gene_mean > 0){
      zvalue_i <- compute_zscore_for_one_gene(gene_mean, zero_proportion, cell.number)
      if(zvalue_i >= Zscore.cutoff){
        selected <- c(selected, i)
        selected.Zscore <- c(selected.Zscore, zvalue_i)
        seleced.num <- seleced.num + 1
      }
      if(seleced.num >= stop_at) break
    }
  }

  res <- list(selected = selected, Zscore = selected.Zscore)
  return(res)
}


select_features_full <- function(X, Zscore.cutoff = 2) {

  gene.number <- nrow(X)
  cell.number <- ncol(X)
  zvalues <- rep(NA, gene.number)

  for(i in 1:gene.number){
    gene_mean <- sum(X[i, ])/cell.number
    zero_proportion <- sum(X[i, ]==0)/cell.number
    if(gene_mean > 0){
      zvalues[i] <- compute_zscore_for_one_gene(gene_mean, zero_proportion, cell.number)
    }
  }

  selected <- which(zvalues >= Zscore.cutoff)
  selected.Zscore <- zvalues[selected]

  res <- list(selected = selected, Zscore = selected.Zscore)
  return(res)
}


check_zero_inflation_numbers <- function(X, Zscore.cutoff = 2) {

  gene.number <- nrow(X)
  cell.number <- ncol(X)
  zvalues <- rep(NA, gene.number)
  for(i in 1:gene.number){
    gene_mean <- sum(X[i, ])/cell.number
    zero_proportion <- sum(X[i, ]==0)/cell.number
    if(gene_mean > 0){
      zvalues[i] <- compute_zscore_for_one_gene(gene_mean, zero_proportion, cell.number)
    }
  }

  selected.len <- length(zvalues[zvalues >= Zscore.cutoff & !is.na(zvalues)])
  return(selected.len)
}


run_kmeans_clustering <- function(subset.dat, npc = 10, km_nstart = 20, km_iter.max = 25){

   pcs <- tryCatch(expr = {
    irlba::irlba(log1p(subset.dat), npc)$v
  }, error = function(e) NA, warning = function(w) NA)
   if (is.na(pcs[1])) {
    return(null_result)
  }

  km_cluster <- kmeans(pcs, 2, nstart = km_nstart, iter.max = km_iter.max)$cluster
  return(km_cluster)

}


initialize_HIPPO <- function(dat, initial.round = 5, stop_at = 500, Zscore.cutoff = 2){

  res <- NULL
  res$sequence <- NULL

  for(i.round in 1:initial.round) {

    if(i.round == 1){
      selected.ID <- select_features_light(dat, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)$selected
      subset.dat <- dat[selected.ID, ]
      clusterID <- run_kmeans_clustering(subset.dat)
      next_round_IDs <- clusterID
    } else {
      tab.IDs <- table(next_round_IDs)
      go_with_larger_number <- as.numeric(names(tab.IDs)[which.max(tab.IDs)])
      selected.dat <- dat[, next_round_IDs%in%go_with_larger_number]
      selected.ID <- select_features_light(selected.dat, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)$selected
      new.subset.dat <- selected.dat[selected.ID, ]
      clusterID <- run_kmeans_clustering(new.subset.dat)
      next_round_IDs[next_round_IDs%in%go_with_larger_number][clusterID == 2] <- i.round + 1
      res$sequence <- c(res$sequence, go_with_larger_number)
    }
    res$next_round_IDs <- next_round_IDs
  }
  return(res)
}

selectCluster_to_proceed <- function(var.list, IDs, cluster.size.cutoff = 100){
  cluster.sizes <- table(IDs)
  passed.clusters <- which(cluster.sizes >= cluster.size.cutoff)
  go_with_larger_variance <- which.max(var.list[passed.clusters])
  selected.cluster <- passed.clusters[go_with_larger_variance]
  return(selected.cluster)
}

selectCluster_to_proceed_inflation <- function(inflation.list, IDs, cluster.size.cutoff = 100){
     cluster.sizes <- table(IDs)
     passed.clusters <- which(cluster.sizes >= cluster.size.cutoff)
     go_with_higher_inflation <- which.max(inflation.list[passed.clusters])
     selected.cluster <- passed.clusters[go_with_higher_inflation]
     return(selected.cluster)
}



#' lightHIPPO
#'
#' @param dat input data matrix (dense or sparse)
#' @param K.round maximum number of rounds
#' @param initial.labels The initial group labels to start with. The default is NULL. This is useful for the situation that you wish to perform sub-clustering on existing group labels. When initial labels are provided, initial.round will be forced to 0.
#' @param initial.round HIPPO rounds using a subset of features, the default is 0. The initial rounds of lighter feature selection will not be supported when given initial group labels.
#' @param stop_at when initial.round > 0, each round will select up to this number of features for clustering and it will not go over all the features.  The initial rounds of lighter feature selection will not be supported when given initial group labels.
#' @param correctByK whether taking into account the number of clusters when calculating the z-score cut-off, default is FALSE
#' @param override.Zscore.cutoff a pre-specified cut-off on Z-scores, the default is NULL. If the default z-score cut-off returns too few features, consider to relax the cut-off here.
#' @param smallest.cluster.num smallest number of cells required in the cluster. Clusters with cells smaller than this number will not be selected for further clustering.
#' @param random.num if smaller than the number of features, only this random set of features will be used to track numbers of inflated genes in each cluster, which will be used to determine the next cluster for further breakdown.
#' @param move.by.inflation the metric used to determine the order in the hierarchical clustering, when this option is TRUE, the procedure moves with the cluster with the largest number of inflated genes, when this option is FALSE, it moves with the cluster with large (variance \times number of cells). The latter one favors larger clusters while considering the within-cluster variance. The default if TRUE.
#' @return A list with clustering results, selected features and z-scores.
#' \itemize{
#'   \item next_round_IDs - the current cluster labels
#'   \item sequence - the sequence of hierarchical clusters
#'   \item selected.gene.list - a list of IDs of selected genes  at each round
#'   \item selected.gene.Zscore - a list of z-scores of selected genes at each round
#'   \item type - "Rooted" if initial.labels = NULL or "Truncated" if generated with initial.labels
#' }
#' @export

lightHIPPO <- function(dat, K.round = 10, initial.labels = NULL, initial.round = 0, stop_at = 500, correctByK = FALSE, override.Zscore.cutoff = NULL, smallest.cluster.num = 200, random.num = 2500, move.by.inflation = TRUE){

  require(irlba)
  total.num.cell <- ncol(dat)
  total.num.gene <- nrow(dat)

  if(!is.null(initial.labels)){
    if(length(initial.labels) != total.num.cell){
      stop("Length of initial group labels doesn't match the number of cell.")
    }
    initial.round <- 0
  }

  if(!is.null(override.Zscore.cutoff)) {
    Zscore.cutoff <- override.Zscore.cutoff
  } else if(correctByK == FALSE){
    Zscore.cutoff <- cut_off_zscore(total.num.gene)
  } else {
    Zscore.cutoff <- cut_off_zscore(total.num.gene*K)
  }

  if(move.by.inflation == TRUE){

    ### calculate the inflation number for each cluster based on a random set of genes ###
    set.seed(1234567)
    randomIDs <- sample(1:total.num.gene, random.num)

    if(is.null(initial.labels) & initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      inflation.tracking <- NULL
      for(i in 1:c(initial.round+1)){
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
      }
      names(inflation.tracking) <- 1:c(initial.round + 1)

      for(i.round in (initial.round + 1):K.round){

        go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in%go_with_higher_inflationID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)
        next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_higher_inflationID)

        inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
        names(inflation.tracking)[i.round+1] <- i.round+1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore
      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else if(is.null(initial.labels) & initial.round == 0) {

      res <- NULL
      res$sequence <- NULL
      inflation.tracking <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      for(i.round in (initial.round + 1):K.round){

        if(i.round == 1){

          selected.res <- select_features_full(dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs <- clusterID
          inflation.tracking[1] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==1], Zscore.cutoff = Zscore.cutoff)
          inflation.tracking[2] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==2], Zscore.cutoff = Zscore.cutoff)
          names(inflation.tracking) <- c(1:2)

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore

        } else {

          go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
          selected.dat <- dat[, next_round_IDs%in% go_with_higher_inflationID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore

          new.subset.dat <- selected.dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
          res$sequence <- c(res$sequence, go_with_higher_inflationID)

          inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
          inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
          names(inflation.tracking)[i.round+1] <- i.round+1

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore
        }

      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else if(!is.null(initial.labels)) {

      res <- NULL
      res$sequence <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      next_round_IDs <- initial.labels

      inflation.tracking <- NULL
      for(i in 1:max(initial.labels)){
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
      }
      names(inflation.tracking) <- 1:max(initial.labels)

      for(i.round in (max(initial.labels)):(K.round+max(initial.labels)-1)){

        go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in% go_with_higher_inflationID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)
        next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_higher_inflationID)

        inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
        names(inflation.tracking)[i.round+1] <- i.round+1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Truncated"
    }

  }

  if(move.by.inflation == FALSE) {

    ### calculating first 10 PCs
    first_pcs <- tryCatch(expr = {
      irlba::irlba(log1p(dat), 10)$v
    }, error = function(e) NA, warning = function(w) NA)

    if(is.null(initial.labels) & initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      cluster.heterogeneity.tracking <- NULL

      kkk <- apply(first_pcs, 2, function(x){
        tapply(x, next_round_IDs, var)
      })
      cluster.heterogeneity.tracking <- apply(kkk, 1, sum)*table(next_round_IDs)

      for(i.round in (initial.round + 1):K.round){

        go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in%go_with_larger_varianceID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)

        new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
          tapply(x, clusterID, var)
        })

        next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_larger_varianceID)

        new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
        cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
        cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
        names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else if(is.null(initial.labels) & initial.round == 0) {

      res <- NULL
      res$sequence <- NULL
      cluster.heterogeneity.tracking <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      for(i.round in (initial.round + 1):K.round){

        if(i.round == 1){

          selected.res <- select_features_full(dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs <- clusterID
          new.clusters.var <- apply(first_pcs, 2, function(x){
            tapply(x, clusterID, var)
          })
          cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
          names(cluster.heterogeneity.tracking) <- c(1:2)

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore

        } else {

          go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
          selected.dat <- dat[, next_round_IDs%in%go_with_larger_varianceID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          selected.Zscore <- selected.res$Zscore

          new.subset.dat <- selected.dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)

          new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
            tapply(x, clusterID, var)
          })

          next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
          res$sequence <- c(res$sequence, go_with_larger_varianceID)

          new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
          cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
          cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
          names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore
        }

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else if(!is.null(initial.labels)) {

      res <- NULL
      res$sequence <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL
      cluster.heterogeneity.tracking <- NULL

      next_round_IDs <- initial.labels

      kkk <- apply(first_pcs, 2, function(x){
        tapply(x, next_round_IDs, var)
      })
      cluster.heterogeneity.tracking <- apply(kkk, 1, sum)*table(next_round_IDs)

      for(i.round in (initial.round + 1):K.round){

        go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
        selected.dat <- dat[, next_round_IDs%in%go_with_larger_varianceID]
        selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
        selected.ID <- selected.res$selected
        selected.Zscore <- selected.res$Zscore

        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)

        new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
          tapply(x, clusterID, var)
        })

        next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_larger_varianceID)

        new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
        cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
        cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
        names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

        selected.gene.list[[i.round]] <- selected.ID
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Truncated"
    }


  }

  return(res)

}


#' lightHIPPO nested
#'
#' The feature selection in this function is nested, i.e., only the features selected from the previous round will be passed to next round. Will not add new features for testing.
#' Nested procedure will significantly reduce the computing time. However if the purpose is to select cluster-specific features, we recommend to use the complete procedure.
#' @param dat input data matrix (dense or sparse)
#' @param K.round maximum number of rounds
#' @param initial.round HIPPO rounds using a subset of features, default is 0
#' @param stop_at when initial.round > 0, each round will select up to this number of features for clustering and it will not go over all the features
#' @param correctByK whether taking into account the number of clusters when calculating the z-score cut-off, default is FALSE
#' @param override.Zscore.cutoff a pre-specified cut-off on Z-scores, the default is NULL.
#' @param smallest.cluster.num smallest number of cells required in the cluster. Clusters with cells smaller than this number will not be selected for further clustering.
#' @param random.num if smaller than the number of features, only this random set of features will be used to track numbers of inflated genes in each cluster, which will be used to determine the next cluster for further breakdown.
#' @param move.by.inflation the metric used to determine the order in the hierarchical clustering, when this option is TRUE, the procedure moves with the cluster with the largest number of inflated genes, when this option is FALSE, it moves with the cluster with large (variance \times number of cells). The latter one favors larger clusters while considering the within-cluster variance. The default if TRUE.
#' @return A list with clustering results, selected features and z-scores.
#' \itemize{
#'   \item next_round_IDs - the current cluster labels
#'   \item sequence - the sequence of hierarchical clusters
#'   \item selected.gene.list - a list of IDs of selected genes  at each round
#'   \item selected.gene.Zscore - a list of z-scores of selected genes at each round
#'   \item type - "Rooted", nested procedure only implements for runs without initial labels

#' }
#' @export

lightHIPPO_nested <- function(dat, K.round = 10, initial.round = 0, stop_at = 500, correctByK = FALSE, override.Zscore.cutoff = NULL, smallest.cluster.num = 200, random.num = 2500, move.by.inflation = TRUE){

  require(irlba)
  total.num.cell <- ncol(dat)
  total.num.gene <- nrow(dat)

  if(!is.null(override.Zscore.cutoff)) {
    Zscore.cutoff <- override.Zscore.cutoff
  } else if(correctByK == FALSE){
    Zscore.cutoff <- cut_off_zscore(total.num.gene)
  } else {
    Zscore.cutoff <- cut_off_zscore(total.num.gene*K)
  }

  if(move.by.inflation){

    ### calculate the inflation number for each cluster based on a random set of genes ###
    set.seed(1235467)
    randomIDs <- sample(1:total.num.gene, random.num)

    if(initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL


      inflation.tracking <- NULL
      for(i in 1:c(initial.round+1)){
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
      }
      names(inflation.tracking) <- 1:c(initial.round + 1)

      nested.feature.list <- NULL
      for(i.round in (initial.round + 1):K.round){

        go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)

        if(i.round == initial.round + 1){

          selected.dat <- dat[, next_round_IDs%in%go_with_higher_inflationID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- selected.ID

        } else {

          selected.dat <- dat[nested.feature.list, next_round_IDs%in%go_with_higher_inflationID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- nested.feature.list[selected.ID]
        }

        selected.Zscore <- selected.res$Zscore
        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)
        next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_higher_inflationID)

        inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
        inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
        names(inflation.tracking)[i.round+1] <- i.round+1

        selected.gene.list[[i.round]] <- nested.feature.list
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else {

      res <- NULL
      res$sequence <- NULL
      inflation.tracking <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      for(i.round in (initial.round + 1):K.round){

        if(i.round == 1){

          selected.res <- select_features_full(dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- selected.ID
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs <- clusterID
          inflation.tracking[1] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==1], Zscore.cutoff = Zscore.cutoff)
          inflation.tracking[2] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs==2], Zscore.cutoff = Zscore.cutoff)
          names(inflation.tracking) <- c(1:2)

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore

        } else {

          go_with_higher_inflationID <- selectCluster_to_proceed_inflation(inflation.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)
          selected.dat <- dat[nested.feature.list, next_round_IDs%in%go_with_higher_inflationID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- nested.feature.list[selected.ID]
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- selected.dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs[next_round_IDs%in%go_with_higher_inflationID][clusterID == 2] <- i.round + 1
          res$sequence <- c(res$sequence, go_with_higher_inflationID)

          inflation.tracking[go_with_higher_inflationID] <- check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%go_with_higher_inflationID], Zscore.cutoff = Zscore.cutoff)
          inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%(i.round+1)], Zscore.cutoff = Zscore.cutoff))
          names(inflation.tracking)[i.round+1] <- i.round+1

          selected.gene.list[[i.round]] <- nested.feature.list
          selected.gene.Zscore[[i.round]] <- selected.Zscore
        }

      }
      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    }


  } else {

    ### calculating first 10 PCs
    first_pcs <- tryCatch(expr = {
      irlba::irlba(log(dat + 1), 10)$v
    }, error = function(e) NA, warning = function(w) NA)

    if(initial.round > 0) {

      initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
      next_round_IDs <- initial_clusters$next_round_IDs
      res <- initial_clusters

      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      kkk <- apply(first_pcs, 2, function(x){
        tapply(x, next_round_IDs, var)
      })
      cluster.heterogeneity.tracking <- apply(kkk, 1, sum)*table(next_round_IDs)

      for(i.round in (initial.round + 1):K.round){

        go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)

        if(i.round == initial.round + 1){

          selected.dat <- dat[, next_round_IDs%in%go_with_larger_varianceID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- selected.ID

        } else {

          selected.dat <- dat[nested.feature.list, next_round_IDs%in%go_with_larger_varianceID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- nested.feature.list[selected.ID]
        }

        selected.Zscore <- selected.res$Zscore
        new.subset.dat <- selected.dat[selected.ID, ]
        clusterID <- run_kmeans_clustering(new.subset.dat)

        new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
          tapply(x, clusterID, var)
        })

        next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
        res$sequence <- c(res$sequence, go_with_larger_varianceID)

        new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
        cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
        cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
        names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

        selected.gene.list[[i.round]] <- nested.feature.list
        selected.gene.Zscore[[i.round]] <- selected.Zscore

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    } else {

      res <- NULL
      res$sequence <- NULL
      cluster.heterogeneity.tracking <- NULL
      selected.gene.list <- NULL
      selected.gene.Zscore <- NULL

      for(i.round in (initial.round + 1):K.round){

        if(i.round == 1){

          selected.res <- select_features_full(dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- selected.ID
          selected.Zscore <- selected.res$Zscore
          new.subset.dat <- dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)
          next_round_IDs <- clusterID
          new.clusters.var <- apply(first_pcs, 2, function(x){
            tapply(x, clusterID, var)
          })
          cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
          names(cluster.heterogeneity.tracking) <- c(1:2)

          selected.gene.list[[i.round]] <- selected.ID
          selected.gene.Zscore[[i.round]] <- selected.Zscore

        } else {

          go_with_larger_varianceID <- selectCluster_to_proceed(cluster.heterogeneity.tracking, next_round_IDs, cluster.size.cutoff = smallest.cluster.num)

          selected.dat <- dat[nested.feature.list, next_round_IDs%in%go_with_larger_varianceID]
          selected.res <- select_features_full(selected.dat, Zscore.cutoff = Zscore.cutoff)
          selected.ID <- selected.res$selected
          nested.feature.list <- nested.feature.list[selected.ID]
          selected.Zscore <- selected.res$Zscore

          new.subset.dat <- selected.dat[selected.ID, ]
          clusterID <- run_kmeans_clustering(new.subset.dat)

          new.clusters.var <- apply(first_pcs[next_round_IDs%in%go_with_larger_varianceID, ], 2, function(x){
            tapply(x, clusterID, var)
          })

          next_round_IDs[next_round_IDs%in%go_with_larger_varianceID][clusterID == 2] <- i.round + 1
          res$sequence <- c(res$sequence, go_with_larger_varianceID)

          new.cluster.heterogeneity.tracking <- apply(new.clusters.var, 1, sum)*table(clusterID)
          cluster.heterogeneity.tracking[go_with_larger_varianceID] <- new.cluster.heterogeneity.tracking[1]
          cluster.heterogeneity.tracking <- c(cluster.heterogeneity.tracking, new.cluster.heterogeneity.tracking[2])
          names(cluster.heterogeneity.tracking)[i.round + 1] <-  i.round + 1

          selected.gene.list[[i.round]] <- nested.feature.list
          selected.gene.Zscore[[i.round]] <- selected.Zscore
        }

      }

      res$next_round_IDs <- next_round_IDs
      names(res$sequence) <- 1:length(res$sequence)+2
      res$selected.gene.list <- selected.gene.list
      res$selected.gene.Zscore <- selected.gene.Zscore
      res$type <- "Rooted"

    }

  }

  return(res)

}

identify_private_features <- function(hippo.res){

  private_features_ID <- NULL
  private_features_Zscore <- NULL
  num_list <- length(hippo.res$selected.gene.list)
  all_num <- 1:num_list
  for(i in all_num){
    geneIDs <- hippo.res$selected.gene.list[[i]]
    zscores <- hippo.res$selected.gene.Zscore[[i]]
    if(length(geneIDs) > 0){
      after.i <- all_num[all_num>i]
      all.other.IDs <- NULL
      for(jj in after.i){
        all.other.IDs <- c(all.other.IDs, unlist(hippo.res$selected.gene.list[[jj]]))
      }
      others <- unique(all.other.IDs)
      private_features_ID[[i]] <- geneIDs[!geneIDs%in%others]
      private_features_Zscore[[i]] <- zscores[!geneIDs%in%others]
    }
  }

  private_features <- list(ID = private_features_ID, Zscore = private_features_Zscore)
  return(private_features)
}



identify_common_features <- function(hippo.res){

  common_features_ID <- NULL
  num_list <- length(hippo.res$selected.gene.list)
  all_num <- 1:num_list
  for(i in all_num){
    geneIDs <- hippo.res$selected.gene.list[[i]]
    if(length(geneIDs) > 0){
      if(is.null(common_features_ID)) {
        common_features_ID <- geneIDs
      } else {
        common_features_ID <- intersect(common_features_ID, geneIDs)
      }
    }
  }

  common_features <- list(ID = common_features_ID)
  return(common_features)
}

#' Post-process HIPPO Features
#'
#'This function will take the HIPPO result as input and return selected features for each round.
#'It will first identify genes that are selected at each round and define those as common features. These genes remain heterogeneous at each round, usually are not interesting. We remove them from the final curated gene list.
#'The function will then identify genes that are private to each round, defined as genes that are only inflated at round $k$, but no longer inflated in later rounds. These genes carry information for separation at round $k$, but their heterogeneity got reconciled by the newly introduced cluster.
#' @param hippo.res output from lightHIPPO
#' @return A list of selected features and z-scores that are ranked by their significance.
#' \itemize{
#'   \item ID - a list of IDs of selected genes at each round
#'   \item Zscore - a list of z-scores of selected genes at each round
#' }
#' @export


organizing_hippo_features <- function(hippo.res){

  private_features <- identify_private_features(hippo.res)
  common_ids <- identify_common_features(hippo.res)$ID
  num.rounds <- length(private_features$ID)
  num.features <- unlist(lapply(private_features$ID, length))

  re_organized <- NULL
  for(i in 1:(num.rounds-1)) {
    if(num.features[i] > 0) {
      de_order <- order(private_features$Zscore[[i]], decreasing = TRUE)
      re_organized$ID[[i]] <- private_features$ID[[i]][de_order]
      re_organized$Zscore[[i]] <- private_features$Zscore[[i]][de_order]
    }
  }

  if(num.features[num.rounds] > 0){
    aa <- private_features$ID[[num.rounds]]
    re_organized_lastround_ID <- aa[which(!aa%in%common_ids)]
    re_organized_lastround_Zscore <- private_features$Zscore[[num.rounds]][which(!aa%in%common_ids)]
    de_order <- order(re_organized_lastround_Zscore, decreasing = TRUE)
    re_organized$ID[[num.rounds]] <- re_organized_lastround_ID[de_order]
    re_organized$Zscore[[num.rounds]] <- re_organized_lastround_Zscore[de_order]
  }
  return(re_organized)
}


#' Clustering pruning
#'
#'This function will take the HIPPO result as input and return clustering labels with the desired number of clusters.
#' @param cluster_res output from lightHIPPO
#' @param K maximum number of clusters
#' @return A vector of cluster labels.
#' @export

cut_hierarchy <- function(cluster_res,  K){
  cluster_sequence <- cluster_res$sequence
  tracebackIDs <- cluster_res$next_round_IDs
  final_cluster_num <- length(cluster_sequence) + 2
  for(rev_ii in c(final_cluster_num:(K+1))){
    tracebackIDs[tracebackIDs%in%rev_ii] <- cluster_sequence[names(cluster_sequence)%in%rev_ii]
  }
  return(tracebackIDs)
}

summarize_for_feature_dot <- function(dataset, clusters){
  ttt1 <- apply(dataset, 1, function(x){
    x <- expm1(x)
    tapply(x, clusters, mean)
  })
  ttt2 <- apply(dataset, 1, function(x){
    tapply(x==0, clusters, mean)
  })
  gene.res <- list(Cluster_mean = ttt1, Cluster_zero_proportion = ttt2)
  return(gene.res)
}


summarize_current_zero_proportions <- function(dataset, clusters){
  ttt1 <- apply(dataset, 1, function(x){
    tapply(x, clusters, mean)
  })
  ttt2 <- apply(dataset, 1, function(x){
    tapply(x==0, clusters, mean)
  })

  melt.ttt1 <- reshape2::melt(ttt1)
  melt.ttt2 <- reshape2::melt(ttt2)
  df <- data.frame(Cluster = melt.ttt1[, 1], gene_mean = melt.ttt1[, 3], zero_proportion = melt.ttt2[, 3])
  return(df)

}

visualize_current_zero_proportions <- function(df){

   require(ggplot2)
   x_intervals <- seq(0, 10, by = 0.1)
    expected_line <- data.frame(aa = x_intervals, bb = exp(-x_intervals))
    g = ggplot2::ggplot(df, ggplot2::aes(x = gene_mean,
                                      y = zero_proportion)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::facet_wrap(~.data$Cluster, ncol = 4) +
    ggplot2::geom_line(data = expected_line, aes(aa, bb), col = "black") +
    ggplot2::xlim(c(0, 10)) +
    ggplot2::theme(legend.position = "none") + ggplot2::theme_bw() +
    ggplot2::ylab("Zero Proportion") +
    ggplot2::xlab("Gene Mean") +
    ggplot2::guides(colour =
                      ggplot2::guide_legend(override.aes =
                                              list(size = 5,alpha = 1),
                                            shape = 19)) +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = 45,hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none", strip.placement = "inside") +
    ggplot2::scale_color_manual(values = c("black", "red"))
    return(g)

}



create_edge_for_sequence <- function(sequence_of_number){
  aa <- sort(sequence_of_number)
  aa.len <- length(aa)
  new.edges <- NULL
  for(i in 1:(aa.len-1)){
    new.edges <- rbind(new.edges, c(aa[i+1], aa[i]))
  }
  return(new.edges)
}

trace_back_hierarchy <- function(cluster_res){

  cluster_sequence <- cluster_res$sequence
  parent.nodes <- c(1, cluster_sequence)
  leaf.nodes <- c(2, as.numeric(names(cluster_sequence)))
  max.child <- max(leaf.nodes)
  n.round <- length(leaf.nodes)
  renumbered.parent.nodes <- seq(n.round+max.child, max.child+1, by = -1)
  parent.layers <- unique(parent.nodes)

  renumbered.child.nodes <- NULL
  for(ii in 1:n.round){
    the.child <- leaf.nodes[ii]
    if(the.child%in%parent.layers){
      renumbered.child.nodes <- c(renumbered.child.nodes, max(renumbered.parent.nodes[parent.nodes%in%the.child]))
    } else {
      renumbered.child.nodes <- c(renumbered.child.nodes, the.child)
    }
  }

  pairs <- cbind(renumbered.parent.nodes, parent.nodes, renumbered.child.nodes)

  added <- NULL
  for(jj in parent.layers){
    the.parent <- jj
    multi.layer.numbers <- pairs[parent.nodes%in%the.parent, 1]
    if(length(multi.layer.numbers) == 1) {
      added <- rbind(added, c(pairs[parent.nodes%in%the.parent, 1], pairs[parent.nodes%in%the.parent, 2]), c(pairs[parent.nodes%in%the.parent, 1], pairs[parent.nodes%in%the.parent, 3]))
    } else{
      bottom.layer <- min(multi.layer.numbers)
      added <- rbind(added, c(pairs[renumbered.parent.nodes%in%bottom.layer, 1], pairs[renumbered.parent.nodes%in%bottom.layer, 2]), c(pairs[renumbered.parent.nodes%in%bottom.layer, 1], pairs[renumbered.parent.nodes%in%bottom.layer, 3]))
      rest.layers <- multi.layer.numbers[multi.layer.numbers!=bottom.layer]
      for(kk in rest.layers){
        added <- rbind(added, c(pairs[renumbered.parent.nodes%in%kk, 1], pairs[renumbered.parent.nodes%in%kk, 3]))
      }
      added <- rbind(added, create_edge_for_sequence(multi.layer.numbers))
    }
  }

  connectivity <- as.data.frame(added)
  connectivity <- connectivity[order(connectivity[, 1], decreasing = T), ]
  colnames(connectivity) <- c("parent", "child")
  siblings_sequence <- connectivity[connectivity$child%in%c(1, leaf.nodes), "child"]

  res <- list(connectivity = connectivity, siblings_sequence = siblings_sequence)
  return(res)
}



#' Visualize hippo hierarchy
#'
#'This function will take the HIPPO result as input and visualize the hierarchy.
#' @param cluster_res output from lightHIPPO
#' @param vertical layout of the plot, vertial if set to be TURE, horizontal otherwise
#' @param colorlist a palette to color different clusters, when set to NULL, a pre-selected list of color will be used.
#' @export

visualize_hippo_hierarchy <- function(cluster_res, vertical = TRUE, colorlist = NULL){


  if(cluster_res$type == "Truncated"){
    stop("Now only support the visualization of a rooted tree.")
  }

  res <- trace_back_hierarchy(cluster_res)
  connectivity <- res$connectivity
  siblings_sequence <- res$siblings_sequence

  des.list <- .dendro.node.descendant(connectivity)
  ans.list <- .dendro.node.ancestor(connectivity)

  one.tree <- des.list[[which.max(unlist(lapply(des.list, length))[-max(connectivity$parent)])]]
  siblings_sequence <- c(siblings_sequence[siblings_sequence%in%one.tree], siblings_sequence[!siblings_sequence%in%one.tree])

  if(vertical){
    xx <- .dendro.node.dis(connectivity)
    yy <- .dendro.node.height(connectivity, siblings_sequence)
  } else {
    yy <- .dendro.node.dis(connectivity)
    xx <- .dendro.node.height(connectivity, siblings_sequence)
  }

  plot.new();plot.window(ylim=c(min(xx) - 0.1, max(xx) + 0.1), xlim=c(min(yy) - 0.1, max(yy) + 0.1))
  par(mar=c(0, 0, 0, 0))
  ##put in background
  for(k in 1:nrow(connectivity)){
    lines(c(yy[connectivity$child[k]],
            yy[connectivity$parent[k]]),
          c(xx[connectivity$child[k]],
            xx[connectivity$parent[k]]),
          col='gray', lwd=60)
  }

  leaf.node <- connectivity$child[!connectivity$child%in%connectivity$parent]

  if(is.null(colorlist)){
    colorlist <- c("turquoise4", "lavender", "slateblue1", "violet", "skyblue1",  "gold", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "cyan", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
  }
  if(length(leaf.node) > length(colorlist)){
    color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    colorlist <- sample(color, length(leaf.node))
  }

  width.set <- 1
  root <- max(connectivity$parent)
  plotted <- NULL
  for(jj in leaf.node){
    ans.jj <- c(jj, ans.list[[jj]])
    if(width.set > 1){
      with.common.ans <- intersect(plotted, ans.jj)
      with.common.ans <- with.common.ans[!with.common.ans%in%root]
      if(length(with.common.ans)==0){
        width.set <- 1
      }
    }
    lines.jj <- create_edge_for_sequence(ans.jj)
    for(kk in 1:nrow(lines.jj)){
      lines(c(yy[lines.jj[kk, 1]],
              yy[lines.jj[kk, 2]]),
            c(xx[lines.jj[kk, 1]],
              xx[lines.jj[kk, 2]]),
            col=colorlist[jj], lwd=60-width.set*5)

    }
    plotted <- c(plotted, ans.jj)
    width.set <- width.set + 1
  }

  for(jj in 1:c(max(connectivity$parent)-1)){
    if(is.null(des.list[[jj]])){
      text(yy[jj], xx[jj], jj, col = "black")
    } else {
      aa <- unlist(des.list[[jj]])
      aa <- aa[aa%in%leaf.node]
      str.jj <- paste(aa,  collapse= ", ")
      text(yy[jj], xx[jj], str.jj, col = "black")
    }
  }



}





.dendro.node.height <- function(connectivity, siblings_sequence){
  cood <- rep(0, max(connectivity$parent))
  all.len <- length(siblings_sequence)
  for(i in 1:all.len){
    cood[siblings_sequence[i]] <- i
  }
  parents <- sort(unique(connectivity$parent))
  for(j in parents){
    child <- connectivity$child[connectivity$parent==j]
    cood[j] <- mean(cood[child])
  }
  return(cood)
}

.dendro.node.dis <- function(connectivity){
  ancestor <- max(connectivity$parent)
  edge.length <- function(z, node){
    count <- 1
    parent <- node
    while(parent <= ancestor){
      parent <- z$parent[z$child==parent]
      if(any(z$child==parent)){
        count <- 1 + count
      }  else {
        break
      }
    }
    return(count)
  }
  edge.degree <- NULL
  for(i in 1:(max(connectivity$child))){
    edge.degree <- c(edge.degree, edge.length(connectivity, i))
  }
  edge.degree <- c(edge.degree, 0)
  return(edge.degree)
}

.dendro.node.descendant <- function(connectivity){

  find.descendant <- function(z, node){
    parent <- node
    descendant.list <- NULL
    while(1){
      parent <- z$child[z$parent%in%parent]
      if(length(parent)==0){
        break
      } else {
        descendant.list <- c(descendant.list, parent)
      }
    }
    return(descendant.list)
  }

  num.of.node <- max(connectivity$parent)
  descendant <- NULL
  for(i in 1:num.of.node){
    descendant[[i]] <- sort(find.descendant(connectivity, i))
  }

  return(descendant)
}


.dendro.node.ancestor <- function(connectivity){

  find.ancestor <- function(z, node){
    child <- node
    ancestor.list <- NULL
    while(1){
      child <- z$parent[z$child%in%child]
      if(length(child)==0){
        break
      } else {
        ancestor.list <- c(ancestor.list, child)
      }
    }
    return(ancestor.list)
  }

  num.of.node <- max(connectivity$parent)
  ancestor <- NULL
  for(i in 1:(num.of.node-1)){
    ancestor[[i]] <- sort(find.ancestor(connectivity, i))
  }

  return(ancestor)
}









