########################################################################################################
# Summarise encodings
########################################################################################################
summarise_trace = function(all_encodings, encodings_histology, samplenames, known_signatures, i, signatures) {
  # count number of samples per signature active
  samples_with_signature_count = rep(0, length(known_signatures))
  all_signs = c()
  for (samplename in samplenames) {
    encoding = all_encodings[[samplename]]
    
    # add the signatures to the list of signatures found in the dataset
    all_signs = c(all_signs, encoding$signatures)
    
    for (j in 1:length(encoding$signatures)) {
      # count another sample with this signature active
      samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] = 
        samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] + 1
      
      # count proportion of bins
      counts = table(encoding$encoding[j,]) / ncol(encoding$encoding)
      
      # Original -> proportion of bins with this same/up/down      
      # encodings_histology[i, which(signatures==encoding$signatures[j])[1]] = encodings_histology[i, which(signatures==encoding$signatures[j])[1]] + ifelse("same" %in% names(counts), counts["same"], 0)
      # encodings_histology[i, which(signatures==encoding$signatures[j])[2]] = encodings_histology[i, which(signatures==encoding$signatures[j])[2]] + ifelse("up" %in% names(counts), counts["up"], 0)
      # encodings_histology[i, which(signatures==encoding$signatures[j])[3]] = encodings_histology[i, which(signatures==encoding$signatures[j])[3]] + ifelse("down" %in% names(counts), counts["down"], 0)
      
      # Select single up/down/both/same from the counts inventory
      if ("up" %in% names(counts) & !"down" %in% names(counts)) {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] + 1
      } else if (!"up" %in% names(counts) & "down" %in% names(counts)) {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] + 1
      } else if ("up" %in% names(counts) & "down" %in% names(counts)) {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[4]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[4]] + 1        
      } else {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] + 1
      }
      
      
    }
    # normalise for number of samples per histology - deactivated as counting number of tumours with signature active
    # encodings_histology[i, ] = encodings_histology[i, ] / num_samples[i]
  }
  return(list(encodings_histology=encodings_histology,
              all_signs=all_signs,
              samples_with_signature_count=samples_with_signature_count,
              encodings_histology=encodings_histology,
              average_change_size=NA))
}

summarise_change_start_end = function(all_encodings, encodings_histology, known_signatures, samplenames, i, signatures) {
  # count number of samples per signature active
  samples_with_signature_count = rep(0, length(known_signatures))
  all_change_sizes = matrix(NA, nrow=length(known_signatures), ncol=length(samplenames))
  all_signs = c()
  for (k in 1:length(samplenames)) {
    samplename = samplenames[k]
    encoding = all_encodings[[samplename]]
    
    # add the signatures to the list of signatures found in the dataset
    all_signs = c(all_signs, encoding$signatures)
    
    for (j in 1:length(encoding$signatures)) {
      # count another sample with this signature active
      samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] = 
        samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] + 1
      
      # Select single up/down/both/same from the counts inventory
      if (is.na(encoding$encoding$direction[j])) {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] + 1
      } else if (encoding$encoding$direction[j]=="up") {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] + 1
      } else if (encoding$encoding$direction[j]=="down") {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] + 1        
      } else {
        # this cannot occur
      }
      
      all_change_sizes[known_signatures==as.character(encoding$signatures[j]), k] = encoding$encoding$changesize[j]
    }
  }
  
  average_change_size = apply(all_change_sizes, 1, mean, na.rm=T)
  
  return(list(encodings_histology=encodings_histology,
              all_signs=all_signs,
              samples_with_signature_count=samples_with_signature_count,
              encodings_histology=encodings_histology,
              average_change_size=average_change_size))
}

# Summarise by taking the largest difference before and after a boundary
summarise_max_change_before_after = function(all_encodings, encodings_histology, known_signatures, samplenames, i, signatures) {
  # count number of samples per signature active
  samples_with_signature_count = rep(0, length(known_signatures))
  samples_with_signature_count_change_up = rep(0, length(known_signatures))
  samples_with_signature_count_change_down = rep(0, length(known_signatures))
  all_change_sizes = matrix(NA, nrow=length(known_signatures), ncol=length(samplenames))
  all_change_sizes_up = matrix(NA, nrow=length(known_signatures), ncol=length(samplenames))
  all_change_sizes_down = matrix(NA, nrow=length(known_signatures), ncol=length(samplenames))
  all_signs = c()
  for (k in 1:length(samplenames)) {
    print(samplenames[k])
    samplename = samplenames[k]
    encoding = all_encodings[[samplename]]
    
    # add the signatures to the list of signatures found in the dataset
    all_signs = c(all_signs, encoding$signatures)
    
    for (j in 1:length(encoding$signatures)) {
      # if (encoding$signatures[j]=="21") {
      # print(encoding$signatures[j])
      # print(paste0("Count before: ", samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])]))
      # print(encoding$encoding$direction[j])
      # }
      # count another sample with this signature active
      samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] = 
        samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])] + 1
      # if (encoding$signatures[j]=="21") {
      # print(paste0("Count after: ", samples_with_signature_count[known_signatures==as.character(encoding$signatures[j])]))
      # }
      
      # Select single up/down/both/same from the counts inventory
      if (is.na(encoding$encoding$direction[j])) {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[1]] + 1
      } else if (encoding$encoding$direction[j]=="up") {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[2]] + 1
      } else if (encoding$encoding$direction[j]=="down") {
        encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] = encodings_histology[i, which(signatures==as.character(encoding$signatures[j]))[3]] + 1        
      } else {
        # this cannot occur
      }
      
      all_change_sizes[known_signatures==as.character(encoding$signatures[j]), k] = encoding$encoding$changesize[j]
      if (is.na(encoding$encoding$direction[j])) {
        # do nothing
      } else if (encoding$encoding$direction[j]=="up") {
        all_change_sizes_up[known_signatures==as.character(encoding$signatures[j]), k] = encoding$encoding$changesize[j]
        samples_with_signature_count_change_up[known_signatures==as.character(encoding$signatures[j])] = samples_with_signature_count_change_up[known_signatures==as.character(encoding$signatures[j])] + 1
      } else if (encoding$encoding$direction[j]=="down") {
        all_change_sizes_down[known_signatures==as.character(encoding$signatures[j]), k] = encoding$encoding$changesize[j]
        samples_with_signature_count_change_down[known_signatures==as.character(encoding$signatures[j])] = samples_with_signature_count_change_down[known_signatures==as.character(encoding$signatures[j])] + 1
      }
    }
  }
  
  average_change_size = apply(all_change_sizes, 1, mean, na.rm=T)
  change_size_up = apply(all_change_sizes_up, 1, mean, na.rm=T)
  change_size_down = apply(all_change_sizes_down, 1, mean, na.rm=T)
  
  # average_change_size = apply(all_change_sizes, 1, median, na.rm=T)
  # change_size_up = apply(all_change_sizes_up, 1, median, na.rm=T)
  # change_size_down = apply(all_change_sizes_down, 1, median, na.rm=T)
  
  return(list(encodings_histology=encodings_histology,
              all_signs=all_signs,
              samples_with_signature_count=samples_with_signature_count,
              encodings_histology=encodings_histology,
              average_change_size=average_change_size,
              change_size_up=change_size_up,
              change_size_down=change_size_down,
              samples_with_signature_count_change_up=samples_with_signature_count_change_up,
              samples_with_signature_count_change_down=samples_with_signature_count_change_down))
}

