########################################################################################################
# functions to summarise a trace
########################################################################################################
encode_trace = function(mixtures, samplename) {
  encoding = matrix(NA, nrow=nrow(mixtures), ncol=ncol(mixtures))
  for (j in 1:nrow(mixtures)) {
    signature = mixtures[j,]
    for (i in 2:(length(signature)-1)) {
      if (abs(signature[i]-signature[i+1]) <= change_detection_threshold) {
        encoding[j,i] = "same"
      } else if ((signature[i]-signature[i+1]) < 0) {
        encoding[j,i] = "up"
      } else {
        encoding[j,i] = "down"
      }
    }
  }
  
  # count the number of changepoints in this sample
  num_changepoints = sum(apply(encoding, 1, function(x) sum(x!="same" & !is.na(x))))
  # get the bin with the lowest ccf/mcn value to save
  lowest_changepoint_bin = max(unlist(apply(encoding, 1, function(x) which(x!="same" & !is.na(x)))))
  lowest_changepoint_bin = ifelse(!is.finite(lowest_changepoint_bin), NA, lowest_changepoint_bin)
  mcn_lowest_bin = ifelse(lowest_changepoint_bin > 0, bins[lowest_changepoint_bin], NA)
  
  return(list(encoding=list(signatures=mixtures[,1], encoding=encoding[,2:(ncol(encoding)-1), drop=F]),
              inventory=data.frame(samplename=samplename,
                                   num_changepoints=num_changepoints,
                                   lowest_changepoint_bin=lowest_changepoint_bin,
                                   mcn_lowest_bin=mcn_lowest_bin)))
}


encode_change_start_end = function(mixtures, samplename, startbin, endbin) {
  endbin = ncol(mixtures)
  encoding = array(F, nrow(mixtures))
  direction = array(NA, nrow(mixtures))
  changesize = array(NA, nrow(mixtures))
  for (j in 1:nrow(mixtures)) {
    if (abs(mixtures[j,endbin]-mixtures[j,startbin]) > change_detection_threshold) {
      encoding[j] = T
      direction[j] = ifelse((mixtures[j,endbin]-mixtures[j,startbin]) > 0, "up", "down")
      changesize[j] = abs(mixtures[j,endbin]-mixtures[j,startbin])
    }
  }
  return(list(encoding=list(signatures=mixtures[,1], encoding=data.frame(encoding=encoding, direction=direction, changesize=changesize)),
              inventory=data.frame()))
}

encode_max_change_before_after = function(mixtures, samplename, startbin, endbin, changes) {
  endbin = ncol(mixtures)
  encoding = array(F, nrow(mixtures))
  direction = array(NA, nrow(mixtures))
  changesize = array(NA, nrow(mixtures))
  if (length(changes) > 0 && changes >= startbin) {
    for (j in 1:nrow(mixtures)) {
      # for (changepoint in changes) {
      # initiate with first and last bins as boundaries
      startpoint = startbin
      endpoint = endbin
      
      for (k in 1:length(changes)) {
        changepoint = changes[k]
        
        # startpoint is previous boundary+1, unless its the first changepoint
        if (k > 1) {
          startpoint = changes[k-1] + 1
        }
        
        # the endpoint is either the next changepoint or the endbin
        if (k < length(changes)) {
          endpoint = changes[k+1]
        } else {
          endpoint = endbin
        }
        
        # print(changepoint)
        # pull all bins before and after the boundary
        before_boundary = mixtures[j, startpoint:changepoint]
        after_boundary = mixtures[j,(changepoint+1):endpoint]
        
        # print(before_boundary)
        # print(after_boundary)
        
        min_max = min(before_boundary)-max(after_boundary)
        max_min = max(before_boundary)-min(after_boundary)
        
        # if (min_max > max_min & min_max > change_detection_threshold) {
        if (min_max < 0 & abs(min_max) > change_detection_threshold) {
          # print("1")
          # new change in this signature
          if (!encoding[j]) {
            # print("1.1")
            encoding[j] = T
            direction[j] = "up"
            changesize[j] = abs(min_max)
          # not new change
          } else if (encoding[j] & direction[j] == "up") {
            # print("1.2")
            # adapt the changesize, if this is larger
            changesize[j] = ifelse(abs(min_max) > changesize[j], abs(min_max), changesize[j])
          } else if (encoding[j] & direction[j] == "down") {
            # print("1.3")
            # change in other direction..
            print(paste0(samplename, " multiple direction changes in signature ", mixtures[j,1]))
            direction[j] = "both"
            changesize[j] = ifelse(abs(min_max) > changesize[j], abs(min_max), changesize[j])
          }
          
        } else if (max_min > 0 & abs(max_min) > change_detection_threshold) {
          # print("2")
          # new change in this signature
          if (!encoding[j]) {
            # print("2.1")
            encoding[j] = T
            direction[j] = "down"
            changesize[j] = abs(max_min)
          } else if (encoding[j] & direction[j] == "down") {
            # print("2.2")
            # adapt the changesize, if this is larger
            changesize[j] = ifelse(max_min > changesize[j], max_min, changesize[j])
          } else if (encoding[j] & direction[j] == "down") {
            # print("2.3")
            # change in other direction..
            print(paste0(samplename, " multiple direction changes in signature ", mixtures[j,1]))
            direction[j] = "both"
            changesize[j] = ifelse(abs(max_min) > changesize[j], abs(max_min), changesize[j])
          }
        }
      }
    }
  }
  return(list(encoding=list(signatures=mixtures[,1], encoding=data.frame(encoding=encoding, direction=direction, changesize=changesize)),
              inventory=data.frame()))
}