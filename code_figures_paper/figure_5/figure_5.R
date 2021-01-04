#
# Creates the signatures overview plot consisting of a number of bar charts, 
# radial plots, pre-created trajectories, further summary figures and legends
#
# and creates supplementary summary trajectories
#

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(dplyr)
library(reshape2)

source("encoding.R")
source("summarise.R")
source("util.R")
source("make_summary_plots.R")
source("barplot_supported_boudaries.r")
source("pcawg.colour.palette.R")

signature_annotations_file = "signature_annotations.txt"

# what to do with clonal mutations
only_subclonal_changes = F
only_bins_below_ccf_1 = T # < most figures so far
only_clone_subclone_boundary = F
only_subclone_subclone_boundary = F

# parameters to set limits of the axis of subplots
if (only_subclone_subclone_boundary) {
  radial_y_max = 40 #35
  
  most_change_min_samples_threshold = 30 #14
  both_change_min_samples_threshold = 30 #14
  
  both_change_y_max = 1
  most_change_y_max = 1
  
} else if (only_clone_subclone_boundary) {
  radial_y_max = 40
  
  most_change_min_samples_threshold = 30
  both_change_min_samples_threshold = 30
  
  both_change_y_max = 1
  most_change_y_max = 1
} else {
  radial_y_max = 30
  
  most_change_min_samples_threshold = 30 
  both_change_min_samples_threshold = 30
  
  both_change_y_max = 0.6
  most_change_y_max = 0.64
}

# select figures how the space is used
use_half = F
use_third = T

if (use_half) {
  postfix = "_half"
} else if (use_third) {
  postfix = "_third"
}

# clone boundary finding settings
boundary_snv_sum = T # Yulias approach
boundary_clonal_clustsize = F

# encoding settings
encode_use_trace = F
encode_use_change_start_end = F
encode_use_max_change_before_after = T # use max change before and after a boundary

change_detection_threshold = 0.06

# Load data
summary_table_file = "consensus_subclonal_reconstruction_v1.1_20181121_summary_table.txt"
drivers = readr::read_tsv("TableS3_panorama_driver_mutations_pcawg_annotated_v1.1.tsv.gz")
summary_table = readr::read_tsv(summary_table_file)
annotations = readr::read_tsv("pcawg.wg11.final_sample_list.txt")

# Yulias figures
# Load these before colours as the traj_means contain that variable as well
max_size = 5
load("traj_means_all_cancer_types_list.RData")
load("individual_tumours_CLL_with_bootstraps.RData")
load("traj_CLL.RData")

load("sig_colors.sigProfiler.RData")
names(sig_colors) = gsub("SBS", "", names(sig_colors))

# Signature colours - these are subset once the full set of known signatures is obtained because sig_colors does not have a colour for every signature
known_signatures_order = c("1", "2+13", 3:12, 14:60)
plot_title_theme = plot_list[[1]]$theme$plot.title

adjust_labels = function(z, clonal_label_y, title_label_y, title_label) {
    # move plot title
    z = z + annotate(geom="text", x=1, y=title_label_y, size=8, label=title_label, hjust=0)
    z = z + theme(plot.title=element_blank(), legend.position="none")
    return(z)
}

# now custom per plot as the y-axis is different
plot_list[["Lymph-CLL"]] = adjust_labels(plot_list[["Lymph-CLL"]], -15, 13, "Lymph-CLL")
plot_list[["Lung-AdenoCA"]] = adjust_labels(plot_list[["Lung-AdenoCA"]], -19, 22, "Lung-AdenoCA")
plot_list[["Eso-AdenoCA"]] = adjust_labels(plot_list[["Eso-AdenoCA"]], -14, 4, "Eso-AdenoCA")
plot_list[["Thy-AdenoCA"]] = adjust_labels(plot_list[["Thy-AdenoCA"]], -15, 15, "Thy-AdenoCA")

# older figures still have the clonal/subclonal labels, these should be removed manually
for (i in 1:length(cll_sig_plot_list)) {
  z = cll_sig_plot_list[[i]]
  z = z + theme(plot.title=element_blank(), legend.position="none")
  cll_sig_plot_list[[i]] = z
}

cll_sig_plot_list[[1]] = cll_sig_plot_list[[1]] + annotate(geom="text", x=1, y=48.0, size=8, label="Lymph-CLL - Sign. 1", hjust=0)
cll_sig_plot_list[[2]] = cll_sig_plot_list[[2]] + annotate(geom="text", x=1, y=48.0, size=8, label="Lymph-CLL - Sign. 5", hjust=0)
cll_sig_plot_list[[3]] = cll_sig_plot_list[[3]] + annotate(geom="text", x=1, y=48.0, size=8, label="Lymph-CLL - Sign. 8", hjust=0)
cll_sig_plot_list[[4]] = cll_sig_plot_list[[4]] + annotate(geom="text", x=1, y=48.0, size=8, label="Lymph-CLL - Sign. 9", hjust=0)


######### uses clusters BEFORE mtimer #########
clusters = readr::read_tsv("consensus_subclonal_reconstruction_v1.1_20181121_structures.txt.gz")
samplenames_tokeep = annotations$samplename[annotations$multi_rep] # & annotations$trinuc_signif_diff

clusters = clusters[clusters$samplename %in% samplenames_tokeep,]

summary_table$histology_abbreviation = annotations$histology_detailed[match(summary_table$samplename, annotations$samplename)]

## FILTER THRESHOLD
# remove cancer types with low counts
# establish counts for each cancer type, later used for filtering
summary_table = summary_table[summary_table$samplename %in% unique(clusters$samplename),]
cancer_type_counts = table(summary_table$histology_abbreviation)

########################################################################################################
# Signature colours
########################################################################################################

get_sig_colors <- function(sig_names) {
  library(metafolio)
  library(RColorBrewer)
  n <- length(sig_names)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual_col_pals <- qual_col_pals[c(1,2,3,6,7),]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-c(4)]
  indices <- c(seq(1,length(col_vector), 3), seq(2,length(col_vector), 3), seq(3,length(col_vector), 3))
  col_vector <- col_vector[indices][-c(19, 24,41)]
  col_vector <- c(gg_color_hue(10), col_vector)
  
  sig_colors <- list()
  for (i in 1:length(sig_names)) {
    sig <- sig_names[i]
    sig_colors[[sig]] <- col_vector[i]
  }
  return(sig_colors)
}

########################################################################################################
# get legend
########################################################################################################
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}

########################################################################################################
# Start script
########################################################################################################
change_inventory = data.frame()
all_encodings = list()
known_signatures = c() # capture all known signatures
num_changepoints_inventory = data.frame() # count changepoints and save bin in which lowest changepoint occurs for power figure
for (samplename in unique(clusters$samplename)) {
  print(samplename)
  
  histology = summary_table$histology_abbreviation[summary_table$samplename==samplename]

    # reset kidney-rcc or melanoma as the directory is named as such on disk
  histology = ifelse(grepl("Kidney-RCC", histology), "Kidney-RCC", histology)
  histology = ifelse(grepl("Skin-Melanoma", histology), "Skin-Melanoma", histology)
  path_prefix = file.path("signature_trajectories", 
                          histology,
                          samplename)
  
  mixtures_file = file.path(path_prefix, "mixtures.csv")
  changepoints_file = file.path(path_prefix, "changepoints.txt")
  
  if (file.exists(mixtures_file) & file.info(changepoints_file)$size > 1) {
    mixtures = read.csv(mixtures_file, header=T, stringsAsFactors=F)
    # Remove SBS denomination
    mixtures[,1] = gsub("SBS", "", mixtures[,1])
    changes = read.table(changepoints_file, header=F, stringsAsFactors=F)
    
    bins = as.numeric(unlist(lapply(stringr::str_split(colnames(mixtures), "X"), function(x) x[2])))
    bins = bins[!is.na(bins)]
    
    # adjust for purity
    bins = bins / summary_table$purity[summary_table$samplename==samplename]
    
    # merge a couple of signatures into a single one - disabled for now
    mixtures = merge_signatures(mixtures)
    
    # find cluster boundaries
    clusters_sample = clusters[clusters$samplename==samplename,]
    # get the boundary between clone and subclone
    if (boundary_clonal_clustsize) {
      boundary = get_clone_boundary_clonal_clustsize(clusters_sample, bins, method="prop_ccf")
    } else if (boundary_snv_sum) {
      boundary = get_transition_points_from_consensus_format(clusters_sample, mixtures, sliding_window=T)
    }
    
    # remove clonal bins    
    if (only_subclonal_changes) {
      # Skip clonal tumours when selecting only subclonal mutations
      if (nrow(clusters_sample)==1) { next() }
      # mixtures data frame contains signature as first column, so ofset boundary by 1
      mixtures = mixtures[,c(1, (boundary[1]+1):ncol(mixtures))]
      bins = bins[boundary[1]:length(bins)]
      # Remove change points no longer in range and adjust the remaining to match the mixtures
      changes = changes[,changes >= boundary[1], drop=F] - boundary[1] + 1
      
    } else if (only_bins_below_ccf_1) {
      selection = which(bins <= 1)
      
      # if there are no subclonal bins, move to the next sample
      if (length(selection)==0) { next() }
      
      # take bins from the first until the first selected bin (i.e. the clone). These will all need to be merged.
      # the bins are ofset by 1 against the mixtures as the mixtures contain the signature as the first column
      # if the clonal bin (min(selection)+1) is bin 2, then there is nothing to merge and this step is skipped
      if (min(selection)+1 > 2) {
        clonal_to_merge = mixtures[,2:(min(selection)+1)]
        clonal_merged = apply(clonal_to_merge, 1, mean)
        mixtures[,min(selection)+1] = clonal_merged
      }
      
      # mixtures data frame contains signature as first column, so ofset selection by 1
      mixtures = mixtures[,c(1, (selection[1]+1):ncol(mixtures))]
      bins = bins[selection[1]:length(bins)]
      # Remove change points no longer in range and adjust the remaining to match the mixtures
      # changes = changes[,changes >= boundary[1], drop=F] - boundary + 1
      changes = changes[,changes %in% selection, drop=F] - min(selection) + 1
      
      
    } else if (only_clone_subclone_boundary) {
      
      # ______X_______1_______2______3
      # X is the clonal peak center
      # want what is before boundary 1 collapsed
      # want what is between boundaries 1 and 2 collapsed
      
      changepoint_bins = bins[as.numeric(changes)]
      subclonal_changepoints = which(changepoint_bins < 1)
      if (length(subclonal_changepoints)==0) { next() }
      
      changepoint_bins = changepoint_bins[subclonal_changepoints]
      # get the subclonal changepoint with CCF closest to 1 as the clone/subclone boundary
      clone_subclone_changepoint = as.numeric(changes[subclonal_changepoints[which.min(abs(changepoint_bins-1))]])
      
      # collapse everything before and after this changepoint
      if (clone_subclone_changepoint > 2) {
        clonal_to_merge = mixtures[,2:clone_subclone_changepoint]
        clonal_merged = apply(clonal_to_merge, 1, mean)
      } else {
        clonal_merged = mixtures[,2]
      }
      
      if ((clone_subclone_changepoint+1) < ncol(mixtures)) {
        subclonal_to_merge = mixtures[,(clone_subclone_changepoint+1):ncol(mixtures)]
        subclonal_merged = apply(subclonal_to_merge, 1, mean)
      } else {
        subclonal_merged = mixtures[,ncol(mixtures)]
      }
      
      # merge the mixtures before and after
      mixtures = mixtures[,1:3]
      mixtures[,2] = clonal_merged
      mixtures[,3] = subclonal_merged
      
      
      # setup dummy bins and changes vectors
      bins = bins[(clone_subclone_changepoint-1):clone_subclone_changepoint]
      bins[1] = 1
      names(bins) = paste0("X",as.numeric(bins))
      colnames(mixtures)[2:3] = names(bins)
      changes = changes[1]
      changes[1] = 2
      
    } else if (only_subclone_subclone_boundary) {
      
      # ______X_______1_______2______3
      # X is the clonal peak center
      # want what is between boundaries 1 and 2 collapsed
      # want what is between boundaries 2 and 3 collapsed
      
      
      changepoint_bins = bins[as.numeric(changes)]
      subclonal_changepoints = which(changepoint_bins < 1)
      # skip if there are no subclonal changepoints
      if (length(subclonal_changepoints)==0) { next() }
      
      changepoint_bins = changepoint_bins[subclonal_changepoints]
      # get the subclonal changepoint with CCF closest to 1 as the clone/subclone boundary
      clone_subclone_changepoint = subclonal_changepoints[which.min(abs(changepoint_bins-1))]
      
      if (length(changes) > clone_subclone_changepoint) {
        selection = as.numeric(changes[(clone_subclone_changepoint):length(changes)])
        mixtures = mixtures[,c(1, (selection[1]+1):ncol(mixtures))]
        bins = bins[(selection[1]):length(bins)]
        # Remove change points no longer in range and adjust the remaining to match the mixtures
        changes = changes[,changes %in% selection, drop=F] - min(selection)
        changes = changes[changes!=0]
      } else {
        # no changepoints beyond the clone/subclone boundary
        next()
      }
    }
    
    
    # store changes for a separate file
    for (index in unlist(changes)) {
      if (index==1) {
        signatures_change_a = rep(0, nrow(mixtures))
      } else {
        signatures_change_a = mixtures[,index+1]-mixtures[,index]
      }
      signatures_change_b = mixtures[,index+2]-mixtures[,index+1]
      signif_change = which(abs(signatures_change_a) > 0.06 | signatures_change_b > 0.06)
      if (length(signif_change) > 0) {
        signatures_that_change = mixtures$X[signif_change]
        change_inventory = rbind(change_inventory, data.frame(samplename=samplename, ccf=bins[index], signature=signatures_that_change, abs_change_before=signatures_change_a[signif_change], abs_change_after=signatures_change_b[signif_change]))
      }
    }
    
    colnames(mixtures)[1] = "Signature"    
    mixtures_m = melt(mixtures, id.vars="Signature")
    names(bins) = colnames(mixtures)[-1]
    mixtures_m$variable = bins[mixtures_m$variable]
    known_signatures = c(known_signatures, unique(mixtures_m$Signature))
    mixtures_m$Signature = factor(mixtures_m$Signature, levels=sort(unique(mixtures_m$Signature)))
    
    # Check if there is data left, cannot detect a changepoint when there is a single bin
    if (ncol(mixtures) > 2) {
      
      if (encode_use_trace) {
        res = encode_trace(mixtures, samplename)
      }
      else if (encode_use_change_start_end) {
        # bins don't have the signature included, so add 1 to the startbin
        res = encode_change_start_end(mixtures, samplename, startbin=which.min(abs(1-bins))+1, endbin=ncol(mixtures))
      }
      else if (encode_use_max_change_before_after) {
        res = encode_max_change_before_after(mixtures, samplename, startbin=which.min(abs(1-bins))+1, endbin=ncol(mixtures), changes=as.numeric(changes))
      }
      
      all_encodings[[samplename]] = res$encoding
      num_changepoints_inventory = rbind(num_changepoints_inventory, res$inventory)
    } else {
      
      # put no detected changepoints in the inventory
      num_changepoints_inventory = rbind(num_changepoints_inventory, data.frame(samplename=samplename,
                                                                                num_changepoints=0,
                                                                                lowest_changepoint_bin=NA,
                                                                                mcn_lowest_bin=NA))
    }
  }
}
known_signatures = sort(unique(known_signatures))
known_signatures_order = known_signatures_order[known_signatures_order %in% known_signatures]
sig_colors = sig_colors[known_signatures_order]

# store just the changed bin and write out the data
change_inventory$abs_change_before[change_inventory$abs_change_before <= 0.06] = 0
change_inventory$abs_change_after[change_inventory$abs_change_after <= 0.06] = 0
write.table(change_inventory, file="signature_change_inventory.txt", quote=F, sep="\t", row.names=F)

####################################################################################################################################
# obtain inventory of signatures going up and own
####################################################################################################################################
num_types = 4
histologies = sort(unique(summary_table$histology_abbreviation))
signatures = unlist(lapply(known_signatures, rep, num_types))
encodings_histology = matrix(0, nrow=length(histologies), ncol=length(known_signatures)*num_types)
row.names(encodings_histology) = histologies
signature_active_count = matrix(NA, nrow=length(histologies), ncol=length(known_signatures))
signature_up_count = matrix(0, nrow=length(histologies), ncol=length(known_signatures))
signature_down_count =matrix(0, nrow=length(histologies), ncol=length(known_signatures))

row.names(encodings_histology) = histologies
num_samples = rep(0, length(known_signatures)) # count how often signatures are found
average_change_size = matrix(NA, nrow=length(histologies), ncol=length(known_signatures))
average_change_size_up = matrix(NA, nrow=length(histologies), ncol=length(known_signatures))
average_change_size_down = matrix(NA, nrow=length(histologies), ncol=length(known_signatures))
for (i in 1:length(histologies)) {
  print(histologies[i])
  samplenames = summary_table$samplename[summary_table$histology_abbreviation==histologies[i]]
  samplenames = samplenames[samplenames %in% names(all_encodings)]
  num_samples[i] = length(samplenames)
  if (length(samplenames) > 0) {
    if (encode_use_trace) {
      res = summarise_trace(all_encodings, samplenames, known_signatures, encodings_histology, i, signatures)
    } else if (encode_use_change_start_end) {
      res = summarise_change_start_end(all_encodings, encodings_histology, known_signatures, samplenames, i, signatures)
    } else if (encode_use_max_change_before_after) {
      res = summarise_max_change_before_after(all_encodings, encodings_histology, known_signatures, samplenames, i, signatures)
    }
    
    samples_with_signature_count = res$samples_with_signature_count
    samples_with_signature_count_change_up = res$samples_with_signature_count_change_up
    samples_with_signature_count_change_down = res$samples_with_signature_count_change_down
    encodings_histology = res$encodings_histology
    average_change_size[i,] = res$average_change_size
    average_change_size_up[i,] = res$change_size_up
    average_change_size_down[i,] = res$change_size_down
    
    for (j in 1:length(samples_with_signature_count)) {
      if (samples_with_signature_count[j] > 0) {
        # If fewer than 5 samples have the signature active, dont keep the data (i.e. doenst show up in the plot)
        if (samples_with_signature_count[j] < 5) {
          encodings_histology[i, signatures==known_signatures[j]] = 0
          average_change_size[i, j] = NA
          average_change_size_up[i, j] = NA
          average_change_size_down[i, j] = NA
        } else {
          encodings_histology[i, signatures==known_signatures[j]] = encodings_histology[i, signatures==known_signatures[j]] / samples_with_signature_count[j]
          signature_active_count[i, j] = samples_with_signature_count[j]
          signature_up_count[i, j] = samples_with_signature_count_change_up[j]
          signature_down_count[i, j] = samples_with_signature_count_change_down[j]
        }
      }
    }
  }
}

all_dat = list()
for (i in 1:length(histologies)) {
  dat = data.frame(signature=factor(unlist(lapply(known_signatures, rep, num_types)), levels=rev(known_signatures)),
                   type=rep(c("same", "up", "down", "both"), length(known_signatures)),
                   value=encodings_histology[i,])
  all_dat[[histologies[i]]] = dat
}

# plot data without 'same' category
data_combined_up = matrix(NA, nrow=length(known_signatures), ncol=length(histologies))
data_combined_down = matrix(NA, nrow=length(known_signatures), ncol=length(histologies))
data_category_up = matrix("up", nrow=length(known_signatures), ncol=length(histologies))
data_category_down = matrix("down", nrow=length(known_signatures), ncol=length(histologies))
for (i in 1:length(all_dat)) {
  dat = all_dat[[i]]
  dat = dat[as.character(dat$type)!="same",]
  r = reshape2::dcast(dat, formula=signature~type, value.var="value", fun.aggregate=sum)
  r = r[match(r$signature, known_signatures),]
  max_category = apply(r[,-1], 1, which.max) + 1
  
  data_combined_up[,i] = r$up
  data_combined_down[,i] = r$down
}

prepare_data = function(data_category_selected, data_combined, change_size, known_signatures, histologies, signature_active_count, signature_change_count) {
  data_category_selected = as.data.frame(data_category_selected)
  colnames(data_category_selected) = histologies
  data_category_selected$signature = known_signatures
  data_category_selected = data_category_selected[, c(ncol(data_category_selected), 1:(ncol(data_category_selected)-1))]
  
  # set x-axis ordering
  known_signatures_order = known_signatures_order[known_signatures_order %in% known_signatures]
  data_category_selected = data_category_selected[match(known_signatures_order, data_category_selected$signature),]
  data_category_selected$signature = factor(data_category_selected$signature, levels=known_signatures_order)
  
  data_category_selected_m = melt(data_category_selected, id.vars="signature")
  data_category_selected_m$value = factor(data_category_selected_m$value, levels=c("up", "down"))
  
  # transform the sizes of the squares
  average_change_size_selected = cbind(data.frame(signature=known_signatures),t(change_size))
  colnames(average_change_size_selected) = c("signature", histologies)
  average_change_size_selected = average_change_size_selected[match(known_signatures_order, average_change_size_selected$signature),]
  average_change_size_selected$signature = factor(average_change_size_selected$signature, levels=known_signatures_order)
  average_change_size_selected_m = melt(average_change_size_selected, id.vars="signature")
  
  data_category_selected_m$changesize_x = average_change_size_selected_m$value #average_change_size_selected_m_normalised
  data_category_selected_m$changesize_y = average_change_size_selected_m$value #average_change_size_selected_m_normalised
  
  # add the proportion of samples that changed
  data_combined = cbind(data.frame(signature=known_signatures), data_combined)
  colnames(data_combined) = c("signature", histologies)
  data_combined = data_combined[match(known_signatures_order, data_combined$signature),]
  data_combined$signature = factor(data_combined$signature, levels=known_signatures_order)
  
  data_combined_m = melt(data_combined, id.vars="signature")
  data_category_selected_m$prop_samples = data_combined_m$value
  
  data_category_selected_m$num_samples = 0
  data_category_selected_m$num_samples_change = 0
  for (i in 1:nrow(data_category_selected_m)) {
    data_category_selected_m$num_samples[i] = signature_active_count[histologies==as.character(data_category_selected_m$variable[i]), 
                                                                     known_signatures==as.character(data_category_selected_m$signature[i])]
    data_category_selected_m$num_samples_change[i] = signature_change_count[histologies==as.character(data_category_selected_m$variable[i]), 
                                                                            known_signatures==as.character(data_category_selected_m$signature[i])]
  }
  
  return(data_category_selected_m)
}

data_category_selected_m_up = prepare_data(data_category_up, data_combined_up, average_change_size_up, known_signatures, histologies, signature_active_count, signature_up_count)
data_category_selected_m_down = prepare_data(data_category_down, data_combined_down, average_change_size_down, known_signatures, histologies, signature_active_count, signature_down_count)
data_category_selected_m = rbind(data_category_selected_m_up, data_category_selected_m_down)

# assign a bin to determine alpha on
data_category_selected_m$prop_samples_category = NA
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples < 0.1] = "x < 0.1"
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples >= 0.1 & data_category_selected_m$prop_samples < 0.2] = "0.1 <= x < 0.2"
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples >= 0.2 & data_category_selected_m$prop_samples < 0.3] = "0.2 <= x < 0.3"
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples >= 0.3 & data_category_selected_m$prop_samples < 0.4] = "0.3 <= x < 0.4"
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples >= 0.4 & data_category_selected_m$prop_samples < 0.5] = "0.4 <= x < 0.5"
data_category_selected_m$prop_samples_category[data_category_selected_m$prop_samples >= 0.5] = "0.5 <= x"
data_category_selected_m$prop_samples_category = factor(data_category_selected_m$prop_samples_category, levels=c("x < 0.1", "0.1 <= x < 0.2", "0.2 <= x < 0.3", "0.3 <= x < 0.4", "0.4 <= x < 0.5", "0.5 <= x"))

data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples < 0.1] = 0.1
data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples >= 0.1 & data_category_selected_m$prop_samples < 0.2] = 0.2
data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples >= 0.2 & data_category_selected_m$prop_samples < 0.3] = 0.3
data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples >= 0.3 & data_category_selected_m$prop_samples < 0.4] = 0.4
data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples >= 0.4 & data_category_selected_m$prop_samples < 0.5] = 0.5
data_category_selected_m$prop_samples_value[data_category_selected_m$prop_samples >= 0.5] = 0.6

####################################################################################################################################
# make the pancan figure
####################################################################################################################################


################################ Signature change summary ################################
most_changing = data_category_selected_m %>% group_by(signature) %>% summarise(num_samples_change=sum(num_samples_change, na.rm=T), num_samples=(sum(num_samples, na.rm=T)/2))
most_changing$fraction = most_changing$num_samples_change / most_changing$num_samples
most_changing = most_changing[most_changing$num_samples > most_change_min_samples_threshold & is.finite(most_changing$fraction),]
most_changing_ordered = most_changing$signature[order(most_changing$fraction, decreasing=T)]
most_changing_ordered = most_changing_ordered[1:10]

plot_data = data_category_selected_m %>% group_by(signature, value) %>% summarise(num_samples_change=sum(num_samples_change, na.rm=T), num_samples=sum(num_samples, na.rm=T))
plot_data$fraction = plot_data$num_samples_change / plot_data$num_samples
plot_data$fraction[as.character(plot_data$value)=="down"] = plot_data$fraction[as.character(plot_data$value)=="down"] * -1
plot_data = plot_data[plot_data$signature %in% most_changing_ordered,]
plot_data$signature = factor(as.character(plot_data$signature), levels=rev(as.character(most_changing_ordered)))

y_axis_ticks = round(unique(c(rev(seq(0, most_change_y_max, 0.3))*-1, seq(0, most_change_y_max, 0.3))), 1)
y_axis_labels = abs(y_axis_ticks)

signature_change_figure = ggplot(plot_data[as.character(plot_data$value)=="up",]) +
  aes(x=signature, y=fraction) +
  geom_bar(position="stack", stat="identity", fill="red") +
  geom_bar(data=plot_data[as.character(plot_data$value)=="down",],position="stack", stat="identity", fill="blue") +
  geom_hline(yintercept=0, size=2) +
  ggtitle("Frac. tumours change") +
  scale_y_continuous(breaks=y_axis_ticks, labels=y_axis_labels, limits=c(most_change_y_max*-1, most_change_y_max)) +
  theme_bw() +
  theme(axis.text = element_text(colour="black",size=16,face="plain"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(colour="black",size=18,face="plain"),
        strip.text.y = element_text(colour="black",size=18,face="plain"),
        legend.position = "none",
        plot.title = plot_title_theme,
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(0.05,0.5,0.5,0.5,"cm")) +
  coord_flip()

################################ Inventory of joint changing signatures ################################
# matrices to keep track of the data
both_active = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))
both_change = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))
both_change_up = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))
both_change_down = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))
both_change_up_down = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))
both_change_down_up = matrix(0, ncol=length(known_signatures_order), nrow=length(known_signatures_order))

for (encoding in all_encodings) {
  # account for pairs of active signatures that change
  active = encoding$signatures[!is.na(encoding$encoding$direction)]
  
  if (length(active) > 1) {
    all_pairs = t(combn(active, 2))
    
    for (i in 1:nrow(all_pairs)) {
      first = as.character(all_pairs[i,1])
      second = as.character(all_pairs[i,2])
      
      first_index = which(first==known_signatures_order)
      second_index = which(second==known_signatures_order)
      
      # make sure the first active signature is always mentioned first in the known_signatures list - ensures all counts end up in the upper triangle
      if (first_index > second_index) {
        first = as.character(all_pairs[i,2])
        second = as.character(all_pairs[i,1])
        
        first_index = which(first==known_signatures_order)
        second_index = which(second==known_signatures_order)
      }
      
      # get direction of change for both signatures
      first_direction = encoding$encoding$direction[encoding$signatures==first]
      second_direction = encoding$encoding$direction[encoding$signatures==second]
      
      if (first_direction=="up" & second_direction=="up") { both_change_up[first_index, second_index] = both_change_up[first_index, second_index] + 1
      } else if (first_direction=="down" & second_direction=="down") { both_change_down[first_index, second_index] = both_change_down[first_index, second_index] + 1
      } else if (first_direction=="up" & second_direction=="down") { both_change_up_down[first_index, second_index] = both_change_up_down[first_index, second_index] + 1
      } else if (first_direction=="down" & second_direction=="up") { both_change_down_up[first_index, second_index] = both_change_down_up[first_index, second_index] + 1 }
      
      if (first_direction!="both" & second_direction!="both") {
        both_change[first_index, second_index] = both_change[first_index, second_index] + 1
      }
    }
  }
  
  # account for pairs of signatures active
  active = encoding$signatures
  if (length(active) > 1) {
    all_pairs = t(combn(active, 2))
    
    for (i in 1:nrow(all_pairs)) {
      first = as.character(all_pairs[i,1])
      second = as.character(all_pairs[i,2])
      
      first_index = which(first==known_signatures_order)
      second_index = which(second==known_signatures_order)
      
      # make sure the first active signature is always mentioned first in the known_signatures list - ensures all counts end up in the upper triangle
      if (first_index > second_index) {
        first = as.character(all_pairs[i,2])
        second = as.character(all_pairs[i,1])
        
        first_index = which(first==known_signatures_order)
        second_index = which(second==known_signatures_order)
      }
      
      # get direction of change for both signatures
      first_direction = encoding$encoding$direction[encoding$signatures==first]
      second_direction = encoding$encoding$direction[encoding$signatures==second] 
      
      both_active[first_index, second_index] = both_active[first_index, second_index] + 1
    }
  }
}

transform_matrix = function(both_active) {
  # explicitly assume data went into the upper triangle
  both_active[lower.tri(both_active)] = NA
  both_active = as.data.frame(both_active)
  both_active$signature = known_signatures_order
  colnames(both_active) = c(known_signatures_order, "signature")
  both_active = melt(both_active, id.vars="signature", na.rm=T)
  return(both_active)
}
both_change_up = transform_matrix(both_change_up)
both_change_down = transform_matrix(both_change_down)
both_change_down_up = transform_matrix(both_change_down_up)
both_change_up_down = transform_matrix(both_change_up_down)

both_change = transform_matrix(both_change)
both_active = transform_matrix(both_active)

colnames(both_change) = c("signature_a", "signature_b", "num_change")
both_change$num_active = both_active$value
both_change$fraction = both_change$num_change / both_change$num_active
both_change$both_up = both_change_up$value
both_change$both_down = both_change_down$value
both_change$up_down = both_change_up_down$value
both_change$down_up = both_change_down_up$value

both_change = both_change[both_change$num_active > both_change_min_samples_threshold,]
both_change = head(both_change[order(both_change$fraction, decreasing=T),], 10)

both_change$frac_both_up = both_change$both_up / both_change$num_active
both_change$frac_both_down = both_change$both_down / both_change$num_active
both_change$frac_up_down = both_change$up_down / both_change$num_active
both_change$frac_down_up = both_change$down_up / both_change$num_active

both_change$signatures = paste0(both_change$signature_a, " / ", both_change$signature_b)
both_change$signatures = factor(both_change$signatures, levels=rev(both_change$signatures))
both_change = both_change[, c("signatures", "frac_both_up", "frac_both_down", "frac_up_down", "frac_down_up")]
both_change_m = melt(both_change)

mycolours = c("red", "blue", "#00B0F6", "purple") #"#D95F02"
both_change_figure = ggplot(both_change_m) + aes(x=signatures, y=value, fill=variable) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=mycolours, labels=c("Up", "Down", "Up a / Down b", "Down a / Up b"), drop=FALSE) +
  ggtitle("Frac. pair change") +
  ylim(0, both_change_y_max) +
  theme_bw() +
  theme(axis.text = element_text(colour="black",size=16,face="plain"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(colour="black",size=18,face="plain"),
        strip.text.y = element_text(colour="black",size=18,face="plain"),
        legend.position = "right",
        legend.text=plot_title_theme,
        legend.title=plot_title_theme,
        plot.title = plot_title_theme,
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(0.05,0.5,0.5,0.5,"cm")) + 
  guides(fill=guide_legend(title="Type of change", 
                           override.aes=list(size=6), 
                           ncol=1, 
                           keywidth=0.3,
                           keyheight=0.3,
                           default.unit="inch",
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5)) +
  coord_flip()

bar_legend = g_legend(both_change_figure)
both_change_figure = both_change_figure + theme(legend.position="none")

####################################################################################################################################
# Radial plots
####################################################################################################################################

################################ Apply filters to show most confident signal ################################
## FILTER THRESHOLD
data_category_selected_m_up = data_category_selected_m[as.character(data_category_selected_m$value)=="up",]
data_category_selected_m_down = data_category_selected_m[as.character(data_category_selected_m$value)=="down",]
# remove data that is not finite and where fewer than 3 samples are changing in a particular direction
to_remove = !is.finite(data_category_selected_m_up$changesize_x) & !is.finite(data_category_selected_m_down$changesize_x) | (data_category_selected_m_up$num_samples_change < 3 & data_category_selected_m_down$num_samples_change < 3)
data_category_selected_m_up = data_category_selected_m_up[!to_remove,]
data_category_selected_m_down = data_category_selected_m_down[!to_remove,]
data_category_selected_m = rbind(data_category_selected_m_up, data_category_selected_m_down)

# mask all remaining signatures where less than 3 are changing (combined with the above filter, this removes signal from cases where )
data_category_selected_m$changesize_x[data_category_selected_m$num_samples_change < 3] = NaN

# Remove cancer types with few samples, threshold set to fill the space reserved for radials
too_few_samples = names(cancer_type_counts)[cancer_type_counts < 35]
data_category_selected_m = data_category_selected_m[!data_category_selected_m$variable %in% too_few_samples,]

################################ establish radials order ################################
# establish an order for the cancer types
signal_sum = array(0, length(unique(data_category_selected_m$variable)))
names(signal_sum) = unique(data_category_selected_m$variable)

signal_diff = array(0, length(unique(data_category_selected_m$variable)))
names(signal_diff) = unique(data_category_selected_m$variable)
for (histology in unique(data_category_selected_m$variable)) {
  test = data_category_selected_m[!is.na(data_category_selected_m$value) & data_category_selected_m$variable==histology,]
  signal_sum[histology] = sum(abs(test$num_samples_change) * test$changesize_x, na.rm=T) / (nrow(test))^2 #sum(test$num_samples) #
  
  a = test[test$value=="up",]
  b = test[test$value=="down",]
  # set the reset signatures due to too few samples active proportions to 0 temporarily
  # this makes sure only the visible signatures in the plot are used for ranking
  a$prop_samples[!is.finite(a$changesize_x)] = 0
  b$prop_samples[!is.finite(b$changesize_x)] = 0
  
  signal_diff[histology] = mean(abs(a$prop_samples - b$prop_samples), na.rm=T)
}
histology_order = names(signal_diff)[order(signal_diff, decreasing=T)]

################################ make plots ################################
# for radial text size and legend text
x_axis_text_theme = plot_list[[1]]$theme$axis.text

all_plots = list()
all_colours = list()
for (i in 1:length(histology_order)) {
  histology = histology_order[i]
  test = data_category_selected_m[!is.na(data_category_selected_m$value) & data_category_selected_m$variable==histology,]
  print(histology)
  print(sum(abs(test$num_samples_change)))
  
  ## FILTER THRESHOLD
  # Require at lmore than one signature to change
  if (length(unique(test$signature)) > 1) {
    
    signs = as.character(test$signature)[as.character(test$value)=="up"]
    # set the order in which the signatures are plot along the radial
    test_up = test[as.character(test$value)=="up",]
    test_down = test[as.character(test$value)=="down",]
    signs_order = signs[order(unlist(lapply(1:nrow(test_up), function(i) { max(test_up$changesize_x[i], test_down$changesize_x[i], na.rm=T) })), decreasing=T)]
    # signs_order = signs[order(test_up$num_samples_change+(test_down$num_samples_change*-1), decreasing=T)]
    test$x_ticks = factor(paste0(test$value, "_", test$signature), levels=c(paste0("up_", signs_order), paste0("down_", rev(signs_order))))
    test$changesize_x = test$changesize_x*100
    test$prop_samples = test$prop_samples*100

    # get colour for strip / border
    if (grepl("Kidney-RCC", histology)) {
      strip_background = element_rect(fill=pcawg.colour.palette(x="Kidney-RCC", scheme='tumour.subtype'), colour="white")
    } else {
      strip_background = element_rect(fill=pcawg.colour.palette(x=histology, scheme='tumour.subtype'), colour="white")
    }
    
    if (histology %in% c("Lung-AdenoCA")) {
      strip_background = element_rect(fill=pcawg.colour.palette(x=histology, scheme='tumour.subtype'), colour="black")
    }
    
    # set colour of the strip text
    strip_text = element_text(colour="black", size=plot_title_theme$size, face="plain", margin=margin(0.1,3,0.1,3,"cm")) #t, r, b, l
    
    signature_labels = c(signs_order, rev(signs_order))
    
    signature_colours = unlist(sig_colors[as.character(test$signature)])
    names(signature_colours) = test$x_ticks
    all_colours[[histology]] = signature_colours
    
    ymax = radial_y_max
    background_coord_red = data.frame(xmin=0.5, xmax=0.5+(length(signature_labels)/2), ymin=0, ymax=ymax)
    background_coord_blue = data.frame(xmin=0.5+(length(signature_labels)/2), xmax=0.5+(length(signature_labels)), ymin=0, ymax=ymax)
    
    all_plots[[histology]] = ggplot_build(ggplot() +  #+  #, width=prop_samples_value
                                            geom_rect(data=background_coord_red, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.08) +
                                            geom_rect(data=background_coord_blue, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.08) +
                                            geom_hline(yintercept=seq(10, ymax, 10), linetype=2, colour="black", alpha=0.4) +
                                            geom_vline(xintercept=1:(length(unique(test$signature))*2), colour="grey", alpha=0.5) +
                                            geom_vline(xintercept=c(0.5, 0.5+length(unique(test$signature))), colour="black") +
                                            geom_bar(data=test, mapping=aes(x=x_ticks, y=changesize_x, fill=x_ticks, alpha=prop_samples_category), stat="identity", colour="black") + #, position="dodge"  
                                            geom_hline(yintercept=6, linetype=2, colour="black", alpha=0.7) +
                                            facet_wrap(~variable, ncol=10) + #ylim(0,60) + 
                                            scale_fill_manual(values=signature_colours, drop=F) +
                                            # scale_colour_manual(values=mycolours) +
                                            scale_x_discrete(labels=signature_labels) +
                                            scale_y_log10() + #limits=c(1.1,60)
                                            scale_alpha_discrete(drop=F) +
                                            coord_cartesian(ylim=c(1,ymax)) +
                                            coord_polar(theta="x") +
                                            theme_bw() +
                                            theme(axis.text.x = element_text(colour="black",size=15,face="plain"),
                                                  axis.title.x = element_blank(),
                                                  axis.text.y = element_blank(),
                                                  axis.ticks = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  strip.text.x = element_text(colour="black",size=20,face="plain"),
                                                  legend.position = "none",
                                                  panel.border = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid = element_blank(),
                                                  plot.margin = margin(0.05,0.5,0.5,0.5,"cm")))$plot #t, r, b, l
  } else {
    print(paste0("Found not enough signature changes for ", histology))
  }
}

####################################################################################################################################
# single cancer type figure for an explanation
####################################################################################################################################

################################ Example plot top left ################################
test_up = data_category_selected_m[as.character(data_category_selected_m$value)=="up" & data_category_selected_m_down$variable=="Lung-AdenoCA",]
signs = as.character(test$signature)[as.character(test$value)=="up"]
test_down = data_category_selected_m[as.character(data_category_selected_m$value)=="down" & data_category_selected_m_down$variable=="Lung-AdenoCA",]
test = rbind(test_up, test_down)

sign_empty = test_up$signature[!(is.finite(test_up$changesize_x) | is.finite(test_down$changesize_x))]
test = test[!test$signature %in% sign_empty,]

signs = as.character(test$signature)[as.character(test$value)=="up"]
# set the order in which the signatures are plot along the radial
test_up = test[as.character(test$value)=="up",]
test_down = test[as.character(test$value)=="down",]
signs_order = signs[order(unlist(lapply(1:nrow(test_up), function(i) { max(test_up$changesize_x[i], test_down$changesize_x[i], na.rm=T) })), decreasing=T)]
# signs_order = signs[order(test_up$num_samples_change+(test_down$num_samples_change*-1), decreasing=T)]
test$x_ticks = factor(paste0(test$value, "_", test$signature), levels=c(paste0("up_", signs_order), paste0("down_", rev(signs_order))))
test$changesize_x = test$changesize_x*100
test$prop_samples = test$prop_samples*100

# get colour for strip / border
if (grepl("Kidney-RCC", histology)) {
  strip_background = element_rect(fill=pcawg.colour.palette(x="Kidney-RCC", scheme='tumour.subtype'), colour="white")
} else {
  strip_background = element_rect(fill=pcawg.colour.palette(x=histology, scheme='tumour.subtype'), colour="white")
}

if (histology %in% c("Lung-AdenoCA")) {
  strip_background = element_rect(fill=pcawg.colour.palette(x=histology, scheme='tumour.subtype'), colour="black")
}

# set colour of the strip text
strip_text = element_text(colour="black",size=12,face="plain", margin=margin(0.1,3,0.1,3,"cm")) #t, r, b, l
if (histology %in% c("ColoRect-AdenoCA", "Skin-Melanoma", "Head-SCC", "Panc-AdenoCA", "Kidney-ChRCC", "CNS-GBM", "CNS-Oligo", "Liver-HCC", "Ovary-AdenoCA", "Kidney-RCC.clearcell", "Kidney-RCC.papillary", "Breast-AdenoCA", "Lymph-BNHL", "Thy-AdenoCA")) {
  strip_text = element_text(colour="white",size=12,face="plain", margin=margin(0.1,3,0.1,3,"cm"))
}

signature_labels = c(signs_order, rev(signs_order))
# signature_labels[signature_labels=="2+13"] = "2/13"

signature_colours = unlist(sig_colors[as.character(test$signature)])
names(signature_colours) = test$x_ticks

ymax = radial_y_max
background_coord_red = data.frame(xmin=0.5, xmax=0.5+(length(signature_labels)/2), ymin=0, ymax=ymax)
background_coord_blue = data.frame(xmin=0.5+(length(signature_labels)/2), xmax=0.5+(length(signature_labels)), ymin=0, ymax=ymax)


p = ggplot() +  #+  #, width=prop_samples_value
  geom_rect(data=background_coord_red, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.08) +
  geom_rect(data=background_coord_blue, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.08) +
  geom_hline(yintercept=seq(10, ymax, 10), linetype=2, colour="black", size=1.2, alpha=0.4) +
  geom_vline(xintercept=1:(length(unique(test$signature))*2), colour="grey", alpha=0.5) +
  geom_vline(xintercept=c(0.5, 0.5+length(unique(test$signature))), colour="black") +
  geom_bar(data=test, mapping=aes(x=x_ticks, y=changesize_x, fill=x_ticks, alpha=prop_samples_category), stat="identity", colour="black") + #, position="dodge"  
  geom_hline(yintercept=6, linetype=2, colour="black", size=1.2) +
  geom_vline(xintercept=1, colour="black", size=1, arrow(length = unit(0.5, "cm"))) +
  facet_wrap(~variable, ncol=10) + #ylim(0,60) + 
  scale_fill_manual(values=signature_colours, drop=F, guide = F) +
  # scale_colour_manual(values=mycolours) +
  scale_x_discrete(labels=signature_labels) +
  scale_y_log10() + #limits=c(1.1,60)
  scale_alpha_discrete(drop=F) +
  coord_cartesian(ylim=c(1,ymax)) +
  coord_polar(theta="x") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=plot_title_theme$size, face="plain"),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(colour="black", size=plot_title_theme$size, face="plain"),
        strip.text.y = element_text(colour="black", size=plot_title_theme$size, face="plain"),
        legend.position = "right",
        legend.text=plot_title_theme,
        legend.title=plot_title_theme,
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0.15,2,2,2,"cm")) +
  guides(alpha=guide_legend(title="Frac. samples", 
                           ncol=1, 
                           keywidth=0.3,
                           keyheight=0.3,
                           default.unit="inch",
                           label.hjust=0,
                           label.vjust=0.5,
                           title.hjust=0.5))

alpha_legend = g_legend(p)
p = p + theme(legend.position="none")

################################ Legend ################################
# dummy figure for the legend
signatures_plot = known_signatures_order[known_signatures_order %in% data_category_selected_m$signature]
colours_plot = unlist(sig_colors[known_signatures_order])
colours_plot = colours_plot[signatures_plot]

temp = data.frame(x=rnorm(length(signatures_plot)), y=rnorm(length(signatures_plot)), col=factor(signatures_plot, levels=signatures_plot))

signature_names_extended = readr::read_tsv(signature_annotations_file)
signature_names_extended$signature = gsub("SBS", "", signature_names_extended$signature)

# add a space in front of low numbers to line the dashes up nicely
signature_names_extended$combined = NA
for (i in 1:nrow(signature_names_extended)) {
  if (signature_names_extended$signature[i] %in% as.character(1:9)) {
    signature_names_extended$combined[i] = paste0("  ", signature_names_extended$signature[i], " - ", signature_names_extended$annotation[i])
  } else {
    signature_names_extended$combined[i] = paste0(signature_names_extended$signature[i], " - ", signature_names_extended$annotation[i])
  }
}

temp$col = signature_names_extended$combined[match(temp$col, signature_names_extended$signature)]

colours_plot_annotated = colours_plot
names(colours_plot_annotated) = signature_names_extended$combined[match(names(colours_plot), signature_names_extended$signature)]

temp$col = factor(temp$col, levels=names(colours_plot_annotated))

smaller_fontsize = plot_title_theme
smaller_fontsize$size = 16

l = ggplot(temp) + aes(x=x, y=y, colour=col) + geom_point(shape=15) + 
  scale_colour_manual(values=colours_plot_annotated,
                      labels=names(colours_plot_annotated), drop=F) +
  theme_bw() + 
  theme(legend.text=smaller_fontsize,
        legend.title=plot_title_theme) +
  guides(colour=guide_legend(title="Signatures", 
                             override.aes=list(size=5), 
                             ncol=2, 
                             keywidth=0.3,
                             keyheight=0,
                             default.unit="inch",
                             label.hjust=0.0,
                             label.vjust=0.5,
                             title.hjust=0.5))
legend <- g_legend(l) 


################################ Build plot ################################
bar_summary = arrangeGrob(signature_change_figure, both_change_figure, ncol=2, widths=c(42/90, 48/90))
row1 = arrangeGrob(bar_summary, bar_legend, legend, alpha_legend, ncol=4, widths=c(2/6, 1/6, 2/6, 1/6))
row2 = arrangeGrob(grobs=all_plots[1:6], ncol=6)

if (only_subclone_subclone_boundary) {
  row3 = arrangeGrob(grobs=all_plots[7:14], ncol=8)
  row4 = arrangeGrob(grobs=all_plots[15:20], ncol=8)
  
} else {
  row3 = arrangeGrob(grobs=all_plots[7:14], ncol=8)
  row4 = arrangeGrob(grobs=all_plots[15:22], ncol=8)
}

row1_2 = arrangeGrob(p, arrangeGrob(row1, row2), ncol=2, widths=c(1,3))
top_plot = arrangeGrob(row1_2, row3, row4, ncol=1, heights=c(2,1,1))

####################################################################################################################################
# Yulias figure
####################################################################################################################################
temp = lapply(plot_list[c("Lymph-CLL", "Lung-AdenoCA", "Eso-AdenoCA", "Thy-AdenoCA")], function(p) p + theme(legend.position="none", plot.margin=margin(0.2,2,0.2,0.2,"cm")))
l = arrangeGrob(grobs=temp, ncol=1)

# summary figures
change_summary = make_summary_plots(summary_table, annotations, summary_table_file)

# CLL summary
temp = lapply(cll_sig_plot_list[1:4], function(p) p + theme(legend.position="none", plot.margin=margin(0.2,0.2,0.2,2,"cm")))
r = arrangeGrob(grobs=temp, ncol=1)

# CLL individual tumours
temp = lapply(plots_cll[c("0b6cd7df-6970-4d60-b7b5-85002a7d8781", "5b4b2312-acb5-4329-8d46-7f93213e3daf", "132f7f2a-b902-4343-aa08-cf6a7af10b9a", "ebc1a26b-9582-4756-acd5-b02d1152319d")], function(p) p$plot + theme(legend.position="none") + scale_x_reverse(expand = c(0,0)))
temp[[4]] = temp[[4]] + scale_x_continuous(expand=c(0,0), breaks=c(1.0, 0.8, 0.6, 0.4, 0.2), labels=c("1.0", "0.8", "0.6", "0.4", "0.2"), trans="reverse")
# temp[[4]] = change_summary
r_samples = arrangeGrob(grobs=temp, ncol=1)

# empty plot as spacer in between rows
blank = ggplot(data.frame(wt=rnorm(10), mpg=rnorm(10)), aes(x = wt, y = mpg)) + geom_blank() + theme_bw() + theme(line = element_blank(),
                                                                                                    text = element_blank(),
                                                                                                    title = element_blank(),
                                                                                                    panel.border = element_blank())

ggsave(arrangeGrob(top_plot, blank, arrangeGrob(l, r, ncol=2), blank, change_summary, ncol=1, heights=c(3.63/9, 0.18/9, 3.63/9, 0.18/9, 1.38/9)), file=paste0("figure_main", postfix, ".pdf"), height=25, width=20*1.125, dpi=300)
ggsave(arrangeGrob(top_plot, blank, arrangeGrob(l, r_samples, ncol=2), blank, change_summary, ncol=1, heights=c(3.63/9, 0.18/9, 3.63/9, 0.18/9, 1.38/9)), file=paste0("figure_main_samples", postfix, ".pdf"), height=25, width=20*1.125, dpi=300)

########################################################################################################################################
# Supplementary figure with all trajectories
########################################################################################################################################

load("traj_means_all_cancer_types_list.RData")
plot_order = names(plot_list)[order(names(plot_list))]
plot_order = c(plot_order[plot_order!="Other" & plot_order!="Kidney-RCC"], "Other")
plot_list = plot_list[plot_order]

for (i in 1:length(plot_order)) {
  temp_plot = plot_list[[i]]
  min_y = min(temp_plot$coordinates$limits$y)
  print(min_y)
  temp_plot = temp_plot + annotate(geom="text", label="Clonal", x=17, y=min_y+0.4, size=4) + annotate(geom="text", label="Subclonal", x=67, y=min_y+0.4, size=4)
  plot_list[[i]] = temp_plot
}

single_row_height = 23.38 / 9
ggsave(arrangeGrob(grobs=plot_list[1:18], ncol=2), file=paste0("figure_supplement_1", postfix, ".pdf"), height=single_row_height*9, width=16.54, dpi=300)
ggsave(arrangeGrob(grobs=plot_list[19:32], ncol=2), file=paste0("figure_supplement_2", postfix, ".pdf"), height=single_row_height*8, width=16.54, dpi=300)

