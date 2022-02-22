## lightHippo


**lightHippo** is a light implementation for **HIPPO**. 

**HIPPO** is a tool that performs single cell UMI data pre-processing using zero-inflation tests. The main differences of **HIPPO** over existing methods include 1) recommending no normalization or imputation and 2) an iterative hierarchical clustering and feature selection procedure. In order to apply **HIPPO** to large data sets, we make following changes to the original release. 

- **lightHippo** only computes z-scores for inflation tests. It will not compute p-values or deviance.
- **lightHippo** only uses kmeans to perform clustering. It can not call alternative clustering algorithms. 
- **lightHippo** performs SVD using the default of `irlba` in package **irlba**. It will not compute other version of PCs.
- **lightHippo** tracks the number of inflated gene for each cluster based on a random set (much smaller), not on all of the genes. 
- **lightHippo** allows early termination of the feature selection at first several rounds. Users can specify the number of rounds with lighter feature selection and the number of features when to stop inflation testing.

In addition, we add the following new characteristics to **lightHippo** procedure: 
- introducing a function that computes z-score cut-off based on the number of input genes. The default is computed using a significance level of 0.1 after correcting for FDR. 
- at each round, all eligible genes will be tested for zero-inflation, compared to the original **HIPPO** release that only genes selected in the previous rounds will be subsetted in the following rounds. 
- implementing a post-processing procedure for selected features. This new function will first identify genes that are selected at each round and define those as common features. These genes remain heterogeneous at each round and at the end. They will be removed from the final curated feature list. The function will then identify genes that are private to each round. These are genes only inflated at round $k$, but no longer inflated in later rounds. These genes carry information for separation at round $k$, but their heterogeneity got reconciled by the newly introduced cluster.
- introducing a function that can prune clustering results to any given $K$. 


## Additional notes 
The burden of **lightHippo** procedure mainly comes from the number of input genes. It is recommended to filter out genes that are only expressed in very few cells. 
We will add a procedure that operates on features selected by the previous round like **HIPPO**, which will significantly reduce computing time. Empirical analysis show very rare new features will be added by the complete procedure. 

## Installation

**lightHippo** can be installed from github directly as follows:

```r
devtools::install_github("ChenMengjie/lightHippo")
```

## **lightHippo** Procedure Step 1: iterative clustering and feature selection 

The main function is `lightHIPPO`, with following arguments:
- `dat` input data matrix (dense or sparse)
- `K` maximum number of clusters
- `initial.round` HIPPO rounds using a subset of features, default is 0
- `stop_at` when initial.round > 0, each round will select up to this number of features for clustering and it will not go over all the features
- `correctByK` whether taking into account the number of clusters when calculating the z-score cut-off, default is FALSE
- `override.Zscore.cutoff` a pre-specified cut-off on Z-scores, the default is NULL.
- `smallest.cluster.num` smallest number of cells required in the cluster. Clusters with cells smaller than this number will not be selected for further clustering.
- `random.num` if smaller than the number of features, only this random set of features will be used to track numbers of inflated genes in each cluster, which will be used to determine the next cluster for further breakdown.

`lightHIPPO` will return a list with clustering lab, selected features and z-scores.

- next_round_IDs : the current cluster labels
- sequence : the sequence of hierarchical clusters
- selected.gene.list : a list of IDs of selected genes  at each round
- selected.gene.Zscore : a list of z-scores of selected genes at each round

```r       
check_ttt <- lightHIPPO(dat, K = 10, initial.round = 0)   
```
You can use the following command to run a lighter **HIPPO** procedure on first 5 rounds, where the zero-inflation tests stop when over 500 features are selected based on the z-score cut-offs. The clustering for first 5 round will run on 500 features. 

```r       
check_ttt_2 <- lightHIPPO(dat, K = 10, initial.round = 5, stop_by = 500)   
```

## **lightHippo** Procedure Step 2: selected feature post-processing

The function `organizing_hippo_features` will take the HIPPO result as input and return selected features for each round. In brief, this new function will remove common features that appear at each round and identify features that are private to each round. These are genes only inflated at round $k$, but no longer inflated in later rounds. These genes carry information for separation at round $k$, but their heterogeneity got reconciled by the newly introduced cluster.

```r       
final_feature_list <- organizing_hippo_features(check_ttt)
```

## **lightHippo** Procedure Step 3:  cluster pruning 

You can use function `cut_hierarchy` to prune the clustering results. This function will take the `lightHIPPO` result as input and return clustering labels with the desired number of clusters.
 
```r  
new_group_labels <- cut_hierarchy(check_ttt, K=4)
```


## **lightHippo** Visualizations 

### Check zero-inflation for each cluster
To speed up, we recommend to check the inflation on a random subset. 
```r 
total.num.gene <- nrow(dat)
randomIDs <- sample(1:total.num.gene, 5000)
summarizing_dat <- summarize_current_zero_proportions(dat[randomIDs, ], check_ttt$next_round_IDs)
plot_dat_per_cluster_inflation <- visualize_current_zero_proportions(summarizing_dat)     
png("lightHIPPO_counts_inflation_check.png")
print(plot_dat_per_cluster_inflation)
dev.off()
```

### Check mean and non-zero percentage for each cluster on top features or a list of markers

```r 
ttID <- cut_hierarchy(check_ttt, K = 4)
check_these <- final_feature_list$ID[[3]]
tt.selected <- summarize_for_feature_dot(dat[check_these , ], ttID)
p <- makeDotplot(input = tt.selected, topN = 50)
ggsave(filename = "test_feature_round3.pdf", plot = p, width = 20, height = 6)
```


### Author

**Mengjie Chen** (U Chicago)

Bug report, comments or questions please send to mengjiechen@uchicago.edu.
