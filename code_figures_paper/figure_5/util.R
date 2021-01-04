merge_signatures = function(mixtures) {
  found_signatures = as.character(mixtures[,1])
  
  summarise_signatures = function(found_signatures, mixtures, to_match_signatures, combined_name) {
    sign_index = which(found_signatures %in% to_match_signatures)
    combined = as.data.frame(t(data.matrix(colSums(mixtures[sign_index,2:ncol(mixtures)]))))
    combined = cbind(data.frame(X=combined_name), combined)
    mixtures = mixtures[-sign_index,]
    mixtures = rbind(mixtures, combined)
    return(mixtures)
  }
  
  # merge 10a and 10b into 10
  if (any(found_signatures %in% c("10a", "10b"))) {
    mixtures = summarise_signatures(found_signatures, mixtures, c("10a", "10b"), "10")
    found_signatures = as.character(mixtures[,1])
  }
  
  
  # 2 and 13 into 2+13
  if (any(found_signatures %in% c("2", "13"))) {
    mixtures = summarise_signatures(found_signatures, mixtures, c("2", "13"), "2+13")
    found_signatures = as.character(mixtures[,1])
  }
  
  # 17a and 17b into 17
  if (any(found_signatures %in% c("17a", "17b"))) {
    mixtures = summarise_signatures(found_signatures, mixtures, c("17a", "17b"), "17")
    found_signatures = as.character(mixtures[,1])
  }
  
  # 7a-d into 7
  if (any(found_signatures %in% c("7a", "7b", "7c", "7d"))) {
    mixtures = summarise_signatures(found_signatures, mixtures, c("7a", "7b", "7c", "7d"), "7")
    found_signatures = as.character(mixtures[,1])
  }
  return(mixtures)
}

########################################################################################################
# functions to find transition points in signature trace
########################################################################################################
get_transition_points_from_consensus_format <- function(clusters_sample, mixtures, sliding_window=T) {
  n_timesteps = ncol(mixtures)-1
  if (!"perc_snvs" %in% colnames(clusters_sample)) {
    clusters_sample$perc_snvs = (clusters_sample$n_snvs / sum(clusters_sample$n_snvs))*100
  }
  n_clusters <- 0
  transition_points <- assignments <- NULL
  
  n_clusters <- nrow(clusters_sample)
  if (sliding_window) {
    number_of_mutation_100 = n_timesteps + 3
  } else {
    number_of_mutation_100 = n_timesteps
  }
  
  transition_points = round(number_of_mutation_100 * cumsum(clusters_sample$perc_snvs / 100))[-length(clusters_sample$perc_snvs)]
  if (sliding_window) {
    transition_points = transition_points - 3
    transition_points[transition_points <= 0] <- 1
  }
  transition_points <- transition_points[transition_points != 0]
  
  assignments = c()
  previous_tp = 1
  if (length(transition_points) > 0) {
    for (i in 1:length(transition_points)) {
      assignments <- c(assignments, rep(toString(i), transition_points[i] - previous_tp))
      previous_tp = transition_points[i]
    }
  }
  assignments <- c(assignments, rep(toString(length(transition_points) + 1), n_timesteps - previous_tp + 1))
  
  # cluster_range = data.frame()
  # for (cluster in unique(assignments)) {
  #   cluster_min = min(bins[assignments==cluster])
  #   cluster_max = max(bins[assignments==cluster])
  #   mean_ccf = mean(ccf$ccf[ccf$mcn <= cluster_max & ccf$mcn >= cluster_min], na.rm=T)
  #   median_ccf = median(ccf$ccf[ccf$mcn <= cluster_max & ccf$mcn >= cluster_min], na.rm=T)
  #   cluster_range = rbind(cluster_range, data.frame(cluster=cluster,
  #                                                   min_ccf=cluster_min,
  #                                                   max_ccf=cluster_max,
  #                                                   mean_ccf=mean_ccf,
  #                                                   median_ccf=median_ccf))
  # }
  
  found_clusters = unique(assignments)
  if (length(found_clusters) > 1) {
    boundary = c()
    for (cluster in found_clusters[-length(found_clusters)]) {
      boundary = c(boundary, max(which(assignments==cluster)))
    }
  } else {
    boundary = 1
  }
  
  
  # stopifnot(n_timesteps == length(assignments))
  
  # return(list(n_clusters, transition_points, assignments))
  return(boundary)
}

get_clone_boundary_clonal_clustsize = function(clusters_sample, bins, method="prop_ccf") {
  
  clone_id = which.min(1-clusters_sample$fraction_cancer_cells)
  highest_ccf_subclone_id = clone_id + 1 # assuming clusters are ordered by CCF
  
  clone_size = clusters_sample$n_snvs[clone_id]
  subclone_size = clusters_sample$n_snvs[highest_ccf_subclone_id]
  
  proportion = clone_size / (clone_size+subclone_size)
  
  clone_bin = which.min(abs(bins - clusters_sample$fraction_cancer_cells[clone_id]))
  subclone_bin = which.min(abs(bins - clusters_sample$fraction_cancer_cells[highest_ccf_subclone_id]))
  
  # Find boundary by obtaining the CCF in between the clone&subclone and then find the bin by matching it's location
  if (method=="prop_ccf") {
    boundary = clusters_sample$fraction_cancer_cells[highest_ccf_subclone_id] + (clusters_sample$fraction_cancer_cells[clone_id]-clusters_sample$fraction_cancer_cells[highest_ccf_subclone_id])*proportion
    boundary = which.min(abs(bins-boundary))
  }
  # Find boundary by assigning bins in between the clone&subclone proportionately by the cluster size
  else if (method=="prop_bins") {
    boundary = clone_bin + ceiling((subclone_bin-clone_bin) * proportion)
  }
  
  return(boundary)
}
