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
#' @param K maximum number of clusters
#' @param initial.round HIPPO rounds using a subset of features, default is 0
#' @param stop_at when initial.round > 0, each round will select up to this number of features for clustering and it will not go over all the features
#' @param correctByK whether taking into account the number of clusters when calculating the z-score cut-off, default is FALSE
#' @param override.Zscore.cutoff a pre-specified cut-off on Z-scores, the default is NULL.
#' @param smallest.cluster.num smallest number of cells required in the cluster. Clusters with cells smaller than this number will not be selected for further clustering.
#' @param random.num if smaller than the number of features, only this random set of features will be used to track numbers of inflated genes in each cluster, which will be used to determine the next cluster for further breakdown.
#' @return A list with clustering results, selected features and z-scores.
#' \itemize{
#'   \item next_round_IDs - the current cluster labels
#'   \item sequence - the sequence of hierarchical clusters
#'   \item selected.gene.list - a list of IDs of selected genes  at each round
#'   \item selected.gene.Zscore - a list of z-scores of selected genes at each round
#' }
#' @export

lightHIPPO <- function(dat, K = 10, initial.round = 0, stop_at = 500, correctByK = FALSE, override.Zscore.cutoff = NULL, smallest.cluster.num = 200, random.num = 2500){

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

  if(initial.round > 0) {

    initial_clusters <- initialize_HIPPO(dat, initial.round = initial.round, stop_at = stop_at, Zscore.cutoff = Zscore.cutoff)
    next_round_IDs <- initial_clusters$next_round_IDs

    res <- initial_clusters

    selected.gene.list <- NULL
    selected.gene.Zscore <- NULL
    ### calculate the variance for each cluster based on a random set of genes ###
    randomIDs <- sample(1:total.num.gene, random.num)

    inflation.tracking <- NULL
    for(i in 1:c(initial.round+1)){
      inflation.tracking <- c(inflation.tracking, check_zero_inflation_numbers(dat[randomIDs, next_round_IDs%in%i], Zscore.cutoff = Zscore.cutoff))
    }
    names(inflation.tracking) <- 1:c(initial.round + 1)

    for(i.round in (initial.round + 1):K){

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

  } else {


    res <- NULL
    res$sequence <- NULL
    inflation.tracking <- NULL
    selected.gene.list <- NULL
    selected.gene.Zscore <- NULL
    randomIDs <- sample(1:total.num.gene, random.num)

    for(i.round in (initial.round + 1):K){

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
  parent.nodes <- cluster_sequence
  leaf.nodes <- as.numeric(names(cluster_sequence))
  max.child <- max(leaf.nodes)
  n.round <- length(leaf.nodes)
  renumbered.parent.nodes <- seq(n.round+max.child, max.child+1, by = -1)
  pairs <- cbind(renumbered.parent.nodes, leaf.nodes)
  parent.layers <- table(parent.nodes)
  more.than.one.layers <- which(parent.layers>1)
  added <- NULL
  for(jj in more.than.one.layers){
    the.parent <- as.numeric(names(parent.layers)[jj])
    multi.layer.numbers <- pairs[parent.nodes%in%the.parent, 1]
    added <- rbind(added, create_edge_for_sequence(multi.layer.numbers))
  }

  nodes <- rbind(pairs, added)
  nodes <- nodes[order(nodes[, 1], decreasing = T), ]
  return(nodes)
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








