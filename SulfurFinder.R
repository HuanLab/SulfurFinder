#SulfurFinder ver. 1.0
#This R program will 1) clean LC-HRMS data, 2) recognize S-containing features, and
#3) predict the number of S
#Brian Low, Feb 12, 2025
#Copyright @ The University of British Columbia

################################################################################
################################################################################
################################################################################

#Load libraries

library("xcms")
library("foreach")
library("doParallel")
library("ranger")
library("xlsx")

#Set directory

setwd("C:/Users/User/Desktop")

#Read in feature table

ft = read.xlsx("demo_feature_table.xlsx", sheetIndex = 1)

################################################################################

#Set parameters

polarity = "negative" #"positive" or "negative" for adduct annotation
mz_tol = 0.01 #m/z tolerance (Da) for MS1 assignment 
rt_tol = 6 #retention time (RT) tolerance (in seconds) for MS1 assignment 
peak_cor = 0.7 #minimum threshold for peak-peak correlation
min_int = 1000 #minimum intensity for isotopic peak extraction
blank_threshold = 3 #minimum intensity for blank filtering, NULL to skip
rt_range = c(900,2340) #RT range (in seconds) to use for filtering, NULL to skip
ms2_tol = 0.05 #m/z tolerance (Da) for MS/MS assignment

annotate_isotopes = T #TRUE to annotate isotopes, FALSE to skip
annotate_adducts = T #TRUE to annotate adducts, FALSE to skip
annotate_isf = T #TRUE to annotate in-source fragments, FALSE to skip. Set FALSE if no MS/MS data was collected.
save = T #TRUE to save results, FALSE to skip

################################################################################
################################################################################
################################################################################

#Clean feature table with RT filter

if(!is.null(rt_range)){
  ft = subset(ft, rt >= rt_range[1] & rt <= rt_range[2])
}

################################################################################

#Clean feature table with blank filter

if(!is.null(blank_threshold)){
  mb_idx = grep("MB", colnames(ft))
  sx_idx = 4:ncol(ft)
  sx_idx = setdiff(sx_idx, mb_idx)
  
  mb_int = apply(ft, 1, function(w) {
    mean(w[mb_idx])
  })
  sx_int = apply(ft, 1, function(w) {
    mean(w[sx_idx])
  })
  
  ft = ft[which(sx_int > blank_threshold*mb_int),]
}

smooth = T
smoothing_level = 2
rt_tol = rt_tol/60

################################################################################

#Sanity check to make sure all feature ids are unique

if(length(unique(ft$featureID)) != nrow(ft)){
  stop("Feature IDs in feature table must all be unique!")
}

#Clean up sample names in ft

clean_names = colnames(ft)[4:ncol(ft)]
clean_names = sub("^X", "", clean_names)
clean_names = gsub("\\.", " ", clean_names)
clean_names = gsub(" ", "-", clean_names)
colnames(ft)[4:ncol(ft)] = clean_names

#Set up backend and read all raw data

numCores = detectCores() - 5
cl = makeCluster(numCores)
registerDoParallel(cl)
invisible(clusterEvalQ(cl, library("xcms")))

ms_files = list.files(pattern = ".mzML")
raw_data = foreach(i = 1:length(ms_files)) %dopar% {
  return(xcmsRaw(ms_files[i], profstep = 0))
}
names(raw_data) = gsub(".mzML", "", ms_files)
names(raw_data) = gsub(" ", "-", names(raw_data))

#Sanity check 

if(!identical(clean_names, names(raw_data))){
  stop("Ordering of raw file names and feature table names do not match...")
}

################################################################################

#General functions

#Extract EIC
get_EIC = function(data, mz, rt){
  
  eic = tryCatch({rawEIC(data, mzrange = c(mz - mz_tol, mz + mz_tol),
                         rtrange = c(rt - rt_tol*60, rt + rt_tol*60))}, 
                 error = function(e) {
                   return(NA)
                 })
  if(anyNA(eic)){return(NULL)}
  
  scan_rt = data@scantime[eic$scan]/60
  
  #Get mz; sometimes the mz does not perfectly match the raw scan
  
  scan_mz = c()
  
  for(i in 1:length(eic$scan)){
    
    #If no intensity, no mz detected
    if(eic$intensity[i] == 0){
      scan_mz[i] = NA
      next
    }
    
    temp_scan = data.frame(getScan(data, eic$scan[i]))
    int_candidates = temp_scan[which(abs(temp_scan$mz - mz) <= mz_tol),]
    int_index = which.min(abs(int_candidates$intensity - eic$intensity[i]))
    
    scan_mz[i] = int_candidates$mz[int_index]
    
  }
  
  eic$rt = scan_rt
  eic$mz = scan_mz
  
  return(eic)
  
}
#Smooth peaks
peak_smooth <- function(x,level = smoothing_level){
  n <- level
  if(length(x) <= 2*n){
    return(x)
  } else if(length(unique(x))==1){
    return(x)
  } else{
    y <- vector(length=length(x))
    for(i in 1:n){
      y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
    }
    for(i in (n+1):(length(y)-n)){
      y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
    }
    for(i in (length(y)-n+1):length(y)){
      y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
    }
    return(y)
  }
} 

######################################################################################

#Annotate isotopes

if(annotate_isotopes){
  
  search_isotopes = foreach(i = 1:nrow(ft), .combine = cbind) %dopar% {
    
    temp_results = rep(";", nrow(ft))
    
    query_mz = ft$mz[i]
    query_rt = ft$rt[i]
    #Use the sample with the highest intensity when constructing EICs
    raw_idx = which.max(ft[i, 4:(3+length(ms_files))])
    
    query_eic = get_EIC(raw_data[[raw_idx]], query_mz, query_rt)
    
    if(smooth){
      query_eic$intensity = peak_smooth(query_eic$intensity)
    }
    
    #Consider only M+1 and M+2
    
    m1 = query_mz + 1.003355
    m2 = query_mz + 1.9958
    
    #Look for M+1
    
    mz_idx = which(abs(ft$mz - m1) <= mz_tol)
    rt_idx = which(abs(ft$rt - query_rt) <= rt_tol*60)
    m1_idx = intersect(mz_idx, rt_idx)
    
    #For each M+1 candidate, check whether it is an isotope
    
    if(length(m1_idx) != 0){
      
      for(j in 1:length(m1_idx)){
        #For each candidate for M+1, construct the EIC and get a correlation
        temp_eic = get_EIC(raw_data[[raw_idx]], ft$mz[m1_idx[j]], query_rt)
        
        if(smooth){
          temp_eic$intensity = peak_smooth(temp_eic$intensity)
        }
        
        temp_cor = suppressWarnings(cor(query_eic$intensity, temp_eic$intensity))
        
        if(is.na(temp_cor)){temp_cor = 0}
        
        #Must be well correlated and have an intensity less than M
        if(temp_cor >= peak_cor & ft[i, raw_idx+3] > ft[m1_idx[j], raw_idx+3]){
          #Write that it is an isotope
          temp_results[m1_idx[j]] = paste0(c(temp_results[m1_idx[j]], paste0(ft$featureID[i],"[M+1]")), 
                                           collapse = ";")
          
        }
        
      }
      
    }
    
    #Consider M+2
    
    mz_idx = which(abs(ft$mz - m2) <= mz_tol)
    rt_idx = which(abs(ft$rt - query_rt) <= rt_tol*60)
    m2_idx = intersect(mz_idx, rt_idx)
    
    if(length(m2_idx) != 0){
      
      for(j in 1:length(m2_idx)){
        #For each candidate for M+1, construct the EIC and get a correlation
        temp_eic = get_EIC(raw_data[[raw_idx]], ft$mz[m2_idx[j]], query_rt)
        
        if(smooth){
          temp_eic$intensity = peak_smooth(temp_eic$intensity)
        }
        
        temp_cor = suppressWarnings(cor(query_eic$intensity, temp_eic$intensity))
        
        if(is.na(temp_cor)){temp_cor = 0}
        
        #Must be well correlated and have an intensity less than M
        if(temp_cor >= peak_cor & ft[i, raw_idx+3] > ft[m2_idx[j], raw_idx+3]){
          #Write that it is an isotope
          temp_results[m2_idx[j]] = paste0(c(temp_results[m2_idx[j]], paste0(ft$featureID[i],"[M+2]")), 
                                           collapse = ";")
          
        }
        
      }
    }
    
    return(temp_results)
    
  }
  
  search_isotopes = data.frame(search_isotopes)
  
  #Dereplicate isotope results
  
  dereplicate_isotopes = foreach(i = 1:nrow(search_isotopes), .combine = c) %dopar% {
    
    temp_result = as.character(search_isotopes[i,])
    
    #If nothing detected for that feature, end
    #Else, combine
    
    if(!any(temp_result != ";")){
      return("")
    } else {
      temp_result = temp_result[temp_result != ";"]
      temp_result = gsub(";", "", temp_result)
      temp_result = paste0(temp_result, collapse = ";")
      return(temp_result)
    }
    
  }
  
  ft$isotopes = dereplicate_isotopes
  
}

################################################################################

#Annotate adducts

if(polarity == "positive"){
  adduct_list = data.frame(c("[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+H-H2O]+", "[M+2Na-H]+"),
                           c(17.0265, 21.9819, 37.9559, -18.0106, 43.9639))
  colnames(adduct_list) = c("adduct", "mz_diff")
} else {
  adduct_list = data.frame(c("[M-H2O-H]-", "[M+Na-2H]-", "[M+Cl]-", "[M+K-2H]-"),
                           c(-18.0106, 21.9819, 35.9767, 37.9559))
  colnames(adduct_list) = c("adduct", "mz_diff")
}

if(annotate_adducts){
  
  search_adducts = foreach(i = 1:nrow(ft), .combine = cbind) %dopar% {
    
    temp_results = rep(";", nrow(ft))
    
    query_mz = ft$mz[i]
    query_rt = ft$rt[i]
    
    #Use sample with the highest intensity
    raw_idx = which.max(ft[i, 4:(3+length(ms_files))])
    
    query_eic = get_EIC(raw_data[[raw_idx]], query_mz, query_rt)
    
    if(smooth){
      query_eic$intensity = peak_smooth(query_eic$intensity)
    }
    
    for(j in 1:nrow(adduct_list)){
      
      adduct_mz = query_mz + adduct_list$mz_diff[j]
      mz_idx = which(abs(ft$mz - adduct_mz) <= mz_tol)
      rt_idx = which(abs(ft$rt - query_rt) <= rt_tol*60)
      adduct_idx = intersect(mz_idx, rt_idx)
      
      if(length(adduct_idx) != 0){
        
        for(k in 1:length(adduct_idx)){
          
          temp_eic = get_EIC(raw_data[[raw_idx]], adduct_mz, query_rt)
          
          if(smooth){
            temp_eic$intensity = peak_smooth(temp_eic$intensity)
          }
          
          temp_cor = suppressWarnings(cor(query_eic$intensity, temp_eic$intensity))
          
          if(is.na(temp_cor)){temp_cor = 0}
          
          if(temp_cor >= peak_cor){
            temp_results[adduct_idx[k]] = paste0(c(temp_results[adduct_idx[k]], paste0(ft$featureID[i], adduct_list$adduct[j])), 
                                                 collapse = ";")
          }
          
        }
        
      } 
      
    }
    
    return(temp_results)
    
  }
  
  search_adducts = data.frame(search_adducts)
  
  #Dereplicate adduct results
  
  dereplicate_adducts = foreach(i = 1:nrow(search_adducts), .combine = c) %dopar% {
    
    temp_result = as.character(search_adducts[i,])
    
    #If nothing detected for that feature, end
    #Else, combine
    
    if(!any(temp_result != ";")){
      return("")
    } else {
      temp_result = temp_result[temp_result != ";"]
      temp_result = gsub(";", "", temp_result)
      temp_result = paste0(temp_result, collapse = ";")
      return(temp_result)
    }
    
  }
  
  ft$adducts = dereplicate_adducts
  
}

################################################################################

#Extract MS/MS and annotate ISFs

stopCluster(cl)
cl = makeCluster(numCores)
registerDoParallel(cl)
invisible(clusterEvalQ(cl, library("xcms")))

extract_msms = function(msms){
  
  msms = strsplit(msms, split = " ")[[1]]
  mz = c()
  int = c()
  
  for(i in 1:length(msms)){
    mz[i] = strsplit(msms[i], split = ":")[[1]][1]
    int[i] = strsplit(msms[i], split = ":")[[1]][2]
  }
  
  msms_df = data.frame(as.numeric(mz), as.numeric(int))
  
  #Remove entries with 0
  
  msms_df = subset(msms_df, int > 0)
  colnames(msms_df) = c("mz", "int")
  
  return(msms_df)
  
}
buddy_denoise = function(msms, max_rsd = 0.25, max_ratio = 0.8, max_int = 5000){
  
  if(nrow(msms) > 10){
    
    msms = msms[order(msms$int),]
    
    #Start with first 3 lowest intense ions
    
    temp_int = msms$int[1:3]
    
    temp_rsd = sd(temp_int)/mean(temp_int)
    temp_ratio = length(temp_int)/nrow(msms)
    temp_threshold = mean(temp_int) + 3*sd(temp_int)
    
    #Iterate until one of the following conditions are met
    
    counter = 1
    
    while(temp_rsd <= max_rsd & temp_ratio <= max_ratio & temp_threshold <= max_int){
      
      ion_grouping = round(0.05*counter*nrow(msms), 0)
      temp_int = msms$int[1:ion_grouping]
      
      if(length(temp_int) < 3){
        temp_int = msms$int[1:3]
      }
      
      temp_rsd = sd(temp_int)/mean(temp_int)
      temp_ratio = length(temp_int)/nrow(msms)
      temp_threshold = mean(temp_int) + 3*sd(temp_int)
      
      counter = counter + 1
      
    }
    
    #Filter
    
    msms = subset(msms, int > temp_threshold)
    msms = msms[order(msms$mz),]
    
  }
  
  return(msms)
  
}
r_dp = function(isf_msms, precursor_msms){
  
  #Normalize intensities
  
  isf_msms$int = isf_msms$int/max(isf_msms$int)
  precursor_msms$int = precursor_msms$int/max(precursor_msms$int)
  
  aligned_table = isf_msms
  aligned_table$precursor_mz = NA
  aligned_table$precursor_int = NA
  
  #Each fragment in the ISF candidate must match to only one precursor_ref
  for(i in 1:nrow(isf_msms)){
    frag_mz = isf_msms$mz[i]
    precursor_idx = any(abs(frag_mz - precursor_msms$mz) <= ms2_tol)
    
    if(precursor_idx){
      precursor_idx = which.min(abs(frag_mz - precursor_msms$mz))
      aligned_table$precursor_mz[i] = precursor_msms$mz[precursor_idx]
      aligned_table$precursor_int[i] = precursor_msms$int[precursor_idx]
      #Remove the matched fragment in the precursor so that it will not
      #be considered again
      precursor_msms = precursor_msms[-precursor_idx,]
    } else {
      aligned_table$precursor_mz[i] = 0
      aligned_table$precursor_int[i] = 0
    }
    
  }
  
  #Match ratio 
  match_ratio = sum(aligned_table$precursor_int != 0)/nrow(isf_msms)
  
  #Cosine
  num = as.numeric(aligned_table$int %*% aligned_table$precursor_int)
  denom = sqrt(sum(aligned_table$int^2)) * sqrt(sum(aligned_table$precursor_int^2))
  score = num/denom
  
  if(is.na(score)){score = 0}
  
  results = data.frame(score, match_ratio)
  
  return(results)
  
  
}

if(annotate_isf){
  
  #Add MS/MS
  
  raw_ms2 = foreach(i = 1:length(ms_files)) %dopar% {
    
    temp_raw = readMSData(ms_files[i], msLevel. = 2)
    precursor_mz = c()
    precursor_rt = c()
    ms2 = c()
    
    for(j in 1:length(temp_raw)){
      
      precursor_mz[j] = temp_raw[[j]]@precursorMz
      precursor_rt[j] = temp_raw[[j]]@rt
      ms2_mz = temp_raw[[j]]@mz
      ms2_int = temp_raw[[j]]@intensity
      temp_ms2 = c()
      
      for(k in 1:length(ms2_mz)){
        temp_ms2[k] = paste0(c(round(ms2_mz[k],4), ms2_int[k]), collapse = ":")
      }
      
      ms2[j] = paste0(temp_ms2, collapse = " ")
      
    }
    
    ms2_spectra = data.frame(precursor_mz, precursor_rt, ms2)
    
    return(ms2_spectra)
    
  }
  
  names(raw_ms2) = gsub(".mzML", "", ms_files)
  names(raw_ms2) = gsub(" ", "-", names(raw_ms2))
  
  msms = foreach(i = 1:nrow(ft), .combine = c) %dopar% {
    
    query_mz = ft$mz[i]
    query_rt = ft$rt[i]
    
    #Retrieve MS/MS from the samples
    
    int_order = order(as.numeric(ft[i,4:(length(ms_files)+3)]), decreasing = T)
    msms_found = F
    counter = 0
    
    while(msms_found == F){
      msms_idx = which(abs(raw_ms2[[int_order[counter+1]]]$precursor_mz - query_mz) <= mz_tol)
      rt_idx = which(abs(raw_ms2[[int_order[counter+1]]]$precursor_rt - query_rt) <= rt_tol*60)
      msms_idx = intersect(msms_idx, rt_idx)
      
      if(length(msms_idx) != 0){
        msms_found = T
      } else {
        counter = counter + 1
      }
      
      #Stop when the last sample is reached
      
      if(counter > length(ms_files)){break}
      
    }
    
    if(msms_found){
      #Consider the case where there are multiple candidates; use the one
      #with the closest RT
      msms_idx = msms_idx[which.min(abs(raw_ms2[[int_order[counter+1]]]$precursor_rt[msms_idx] - 
                                          query_rt))]
      return(raw_ms2[[int_order[counter+1]]]$ms2[msms_idx])
    } else {
      return("")
    }
    
  }
  
  ft$msms = msms
  
  #Reformat some MS/MS output
  
  ft$msms[grep("NA", ft$msms)] = ""
  
  #Remove raw_ms2 data 
  
  rm(raw_ms2)
  
  #Begin ISF annotation
  #For each feature mz, look for its potential ISFs
  
  isf = foreach(i = 1:nrow(ft), .combine = cbind) %dopar% {
    
    temp_ft = ft
    temp_ft$isf = ""
    
    query_mz = temp_ft$mz[i]
    query_rt = temp_ft$rt[i]
    raw_idx = which.max(temp_ft[i, 4:(3+length(ms_files))])
    query_eic = get_EIC(raw_data[[raw_idx]], query_mz, query_rt)
    
    if(smooth){
      query_eic$intensity = peak_smooth(query_eic$intensity)
    }
    
    query_msms = temp_ft$msms[i]
    
    if(!identical(query_msms, "")){
      query_msms = extract_msms(query_msms)
      #Denoise
      query_msms = buddy_denoise(query_msms)
      #Normalize
      query_msms$int = query_msms$int/max(query_msms$int)*100
    }
    
    #Level 3 ISF -> only look at peak correlation
    lvl3_idx = which(abs(temp_ft$rt - query_rt) <= rt_tol * 60)
    lvl3_idx = lvl3_idx[-which(lvl3_idx == i)]
    
    #If nothing, next
    
    if(length(lvl3_idx) == 0){return(temp_ft$isf)}
    
    #Index the candidates that could be level 2 ISF
    lvl2_idx = c()
    
    for(j in 1:length(lvl3_idx)){
      
      temp_mz = temp_ft$mz[lvl3_idx[j]]
      
      #Cannot be ISF if higher mass than the query
      if(temp_mz >= query_mz){next}
      
      temp_eic = get_EIC(raw_data[[raw_idx]], temp_mz, query_rt)
      
      if(smooth){
        temp_eic$intensity = peak_smooth(temp_eic$intensity)
      }
      
      temp_cor = suppressWarnings(cor(query_eic$intensity, temp_eic$intensity))
      
      if(is.na(temp_cor)){temp_cor = 0}
      
      if(temp_cor >= peak_cor){
        temp_ft$isf[lvl3_idx[j]] = paste0(c(temp_ft$isf[lvl3_idx[j]], paste0(temp_ft$featureID[i], "[L3]")), 
                                          collapse = ";")
        lvl2_idx = c(lvl2_idx, lvl3_idx[j])
      }
      
    }
    
    #Return what we have so far if cannot proceed
    if(is.null(lvl2_idx) | identical(query_msms, "")){
      return(temp_ft$isf)
    }
    
    #Level 2 ISF --> the precursor exists in the query's MS/MS
    lvl1_idx = c()
    
    for(j in 1:length(lvl2_idx)){
      
      temp_mz = temp_ft$mz[lvl2_idx[j]]
      if(any(abs(query_msms$mz - temp_mz) <= ms2_tol)){
        temp_ft$isf[lvl2_idx[j]] = paste0(c(temp_ft$isf[lvl2_idx[j]], paste0(temp_ft$featureID[i], "[L2]")), 
                                          collapse = ";")
        lvl1_idx = c(lvl1_idx, lvl2_idx[j])
      }
      
    }
    
    #Stop if we cannot proceed further
    if(is.null(lvl1_idx)){
      return(temp_ft$isf)
    }
    
    #Level 1 ISF --> cosine score and match ratio
    
    for(j in 1:length(lvl1_idx)){
      
      temp_msms = temp_ft$msms[lvl1_idx[j]]
      
      if(temp_msms == ""){
        next
      } else {
        temp_msms = extract_msms(temp_msms)
        temp_msms = buddy_denoise(temp_msms)
        temp_msms$int = temp_msms$int/max(temp_msms$int)*100
      }
      
      #Not enough fragments in the experimental MS2 spectra
      if(nrow(temp_msms) < 5){next}
      
      msms_dp = r_dp(temp_msms, query_msms)
      
      if(msms_dp$score >= 0.7 | msms_dp$match_ratio >= 0.7){
        temp_ft$isf[lvl1_idx[j]] = paste0(c(temp_ft$isf[lvl1_idx[j]], paste0(temp_ft$featureID[i], "[L1]")), 
                                          collapse = ";")
      }
      
    }
    
    return(temp_ft$isf)
    
  }
  
  isf = data.frame(isf)
  
  #Clean up and dereplicate ISF results; only consider level 1 ISF
  isf_results = foreach(i = 1:nrow(isf), .combine = c) %dopar% {
    
    query = as.character(isf[i,])
    query = query[query != ""]
    
    if(length(query) == 0){
      return("")
    }
    
    query = sub("^;", "", query)
    query = strsplit(query, split = ";")
    query = unlist(query)
    
    id = c()
    level = c()
    
    for(j in 1:length(query)){
      id[j] = gsub("\\[.*\\]", "", query[j])
      level[j] = gsub(".*\\[L([0-9]+)\\].*", "\\1", query[j])
    }
    
    temp_table = data.frame(id, level, query)
    temp_table$id = as.numeric(temp_table$id)
    temp_table$level = as.numeric(temp_table$level)
    temp_table = temp_table[order(temp_table$level),]
    
    #Consider level 1 ISF and above only
    
    temp_table = subset(temp_table, level == 1)
    
    if(nrow(temp_table) == 0){
      return("")
    }
    
    #Only keep the highest level for each unique entry
    
    unique_id = unique(temp_table$id)
    temp_result = c()
    
    #In cases where we want more than level 1 ISF
    
    for(j in 1:length(unique_id)){
      temp_id = unique_id[j]
      id_idx = which(temp_table$id == temp_id)
      grouped_id = temp_table[id_idx,]
      temp_result[j] = grouped_id$query[which.min(grouped_id$level)]
    }
    
    temp_result = paste0(temp_result, collapse = ";")
    
    return(temp_result)
    
  }
  
  ft$isf = isf_results
  
}

################################################################################

#For each feature, extract information to use for sulfur model

clean_int = function(x){
  x = x[is.finite(x)]
  if(any(x == 0)){
    x = x[-which(x == 0)]
  }
  return(x)
}

out = foreach(i = 1:nrow(ft), .combine = rbind) %dopar% {
  
  query_mz = ft$mz[i]
  query_rt = ft$rt[i]
  
  #Extract chromatograms for M and isotopes
  
  m1 = query_mz + 1.003355
  m2 = query_mz + 1.9958
  m3 = query_mz + 1.003355 + 1.9958
  m4 = query_mz + 1.9958 + 1.9958
  
  #Only consider the isotope pattern of the sample with the highest intensity
  
  raw_idx = which.max(as.numeric(ft[i,4:(length(ms_files)+3)]))
  
  m0_eic = get_EIC(raw_data[[raw_idx]], mz = query_mz, rt = query_rt)
  m1_eic = get_EIC(raw_data[[raw_idx]], mz = m1, rt = query_rt)
  m2_eic = get_EIC(raw_data[[raw_idx]], mz = m2, rt = query_rt)
  m3_eic = get_EIC(raw_data[[raw_idx]], mz = m3, rt = query_rt)
  m4_eic = get_EIC(raw_data[[raw_idx]], mz = m4, rt = query_rt)
  
  #Need to consider the possibility that there is no scans within the time range
  #for small scan windows
  
  if(is.null(m0_eic)){
    return(rep(0, 58))
  }
  
  #Find which m0 scan to use considering relative and absolute intensities
  
  scan_index = which(m0_eic$intensity > min_int & m0_eic$intensity > 
                       (max(m0_eic$intensity)/10))
  
  if(length(scan_index) == 0){
    return(rep(0, 58))
  }
  
  #Smooth
  
  if(smooth){
    m0_eic$intensity = peak_smooth(m0_eic$intensity)
    m1_eic$intensity = peak_smooth(m1_eic$intensity)
    m2_eic$intensity = peak_smooth(m2_eic$intensity)
    m3_eic$intensity = peak_smooth(m3_eic$intensity)
    m4_eic$intensity = peak_smooth(m4_eic$intensity)
  }
  
  #Take correlation
  
  m1_cor = round(suppressWarnings(cor(m0_eic$intensity, m1_eic$intensity)),2)
  m2_cor = round(suppressWarnings(cor(m0_eic$intensity, m2_eic$intensity)),2)
  m3_cor = round(suppressWarnings(cor(m0_eic$intensity, m3_eic$intensity)),2)
  m4_cor = round(suppressWarnings(cor(m0_eic$intensity, m4_eic$intensity)),2)
  
  if(is.na(m1_cor)){m1_cor = 0}
  if(is.na(m2_cor)){m2_cor = 0}
  if(is.na(m3_cor)){m3_cor = 0}
  if(is.na(m4_cor)){m4_cor = 0}
  
  mz0 = mean(na.omit(m0_eic$mz[scan_index]))
  int0 = m0_eic$intensity[scan_index]
  
  if(sum(int0 == 0) == length(scan_index)){
    return(rep(0, 58))
  }
  
  #Isotope information
  
  if(m1_cor >= peak_cor & sum(!is.na(m1_eic$mz[scan_index])) >= 1){
    mz1 = mean(na.omit(m1_eic$mz[scan_index]))
    int1 = m1_eic$intensity[scan_index]
  } else {
    mz1 = 0
    int1 = 0
  }
  
  if(m2_cor >= peak_cor & sum(!is.na(m2_eic$mz[scan_index])) >= 1){
    mz2 = mean(na.omit(m2_eic$mz[scan_index]))
    int2 = m2_eic$intensity[scan_index]
  } else {
    mz2 = 0
    int2 = 0
  }
  
  if(m3_cor >= peak_cor & sum(!is.na(m3_eic$mz[scan_index])) >= 1){
    mz3 = mean(na.omit(m3_eic$mz[scan_index]))
    int3 = m3_eic$intensity[scan_index]
  } else {
    mz3 = 0
    int3 = 0
  }
  
  if(m4_cor >= peak_cor & sum(!is.na(m4_eic$mz[scan_index])) >= 1){
    mz4 = mean(na.omit(m4_eic$mz[scan_index]))
    int4 = m4_eic$intensity[scan_index]
  } else {
    mz4 = 0
    int4 = 0
  }
  
  int1 = mean(clean_int(int1/int0)) * 100
  int2 = mean(clean_int(int2/int0)) * 100
  int3 = mean(clean_int(int3/int0)) * 100
  int4 = mean(clean_int(int4/int0)) * 100
  int0 = 100
  
  if(is.na(int1)){int1 = 0}
  if(is.na(int2)){int2 = 0}
  if(is.na(int3)){int3 = 0}
  if(is.na(int4)){int4 = 0}
  
  #Clean up intensity variables
  
  if(int1 >= 100){
    mz1 = 0
    int1 = 0
  }
  
  if(int2 >= 100){
    mz2 = 0
    int2 = 0
  }
  
  if(int3 >= 100){
    mz3 = 0
    int3 = 0
  }
  
  if(int4 >= 100){
    mz4 = 0
    int4 = 0
  }
  
  #Do not allow int3 or in4 to be larger than int2
  
  if(int3 > int2){
    mz3 = 0
    int3 = 0
  }
  
  if(int4 > int2){
    mz4 = 0
    int4 = 0
  }
  
  if(int2 >= 25){
    notes = "Cl flag"
  } else {
    notes = ""
  }
  
  #Compile isotope pattern
  
  iso_mz = round(c(mz0, mz1, mz2, mz3, mz4), 4)
  iso_mz = paste0(iso_mz, collapse = ",")
  iso_int = round(c(int0, int1, int2, int3, int4),2)
  iso_int = paste0(iso_int, collapse = ",")  
  iso_pattern = paste0(c(iso_mz, iso_int), collapse = ";")
  
  #Select a model
  
  if(int4 > 0){
    model = 3 
  } else if(int3 > 0){
    model = 2
  } else {
    model = 1
  }
  
  #Calculate mz variables
  
  mz1_mz0 = mz1 - mz0
  mz2_mz0 = mz2 - mz0
  mz3_mz0 = mz3 - mz0
  mz4_mz0 = mz4 - mz0
  mz2_mz1 = mz2 - mz1
  mz3_mz1 = mz3 - mz1
  mz4_mz1 = mz4 - mz1
  mz3_mz2 = mz3 - mz2
  mz4_mz2 = mz4 - mz2
  mz4_mz3 = mz4 - mz3
  
  #Calculate intensity variables
  
  int1_int0 = int1 - int0
  int2_int0 = int2 - int0
  int3_int0 = int3 - int0
  int4_int0 = int4 - int0
  int2_int1 = int2 - int1
  int3_int1 = int3 - int1
  int4_int1 = int4 - int1
  int3_int2 = int3 - int2
  int4_int2 = int4 - int2
  int4_int3 = int4 - int3
  
  #Do not predict if int2 if too low as most likely no sulfur
  
  if(int2 < 2.6 | int1 == 0){
    to_predict = F
  } else {
    to_predict = T
  }
  
  #Save first layer variables
  
  first_layer = c(int1, int2, int3, int4, mz0, mz1_mz0 , mz2_mz0 , mz3_mz0,
                  mz4_mz0, mz2_mz1, mz3_mz1, mz4_mz1, mz3_mz2, mz4_mz2, mz4_mz3, 
                  int1_int0, int2_int0, int3_int0, int4_int0, int2_int1, int3_int1,
                  int4_int1, int3_int2, int4_int2, int4_int3, iso_pattern, notes,
                  to_predict, model, m1_cor, m2_cor, m3_cor, m4_cor)
  
  #Additional filters for coeluting ions
  
  if(int2 != 0 & int3 != 0 & int3/int2 > 0.27){
    mz3 = 0
    int3 = 0
  }
  
  if(int2 != 0 & int4 != 0 & int4/int2 > 0.23){
    mz4 = 0
    int4 = 0
  }
  
  #Calculate mz variables
  
  mz1_mz0 = mz1 - mz0
  mz2_mz0 = mz2 - mz0
  mz3_mz0 = mz3 - mz0
  mz4_mz0 = mz4 - mz0
  mz2_mz1 = mz2 - mz1
  mz3_mz1 = mz3 - mz1
  mz4_mz1 = mz4 - mz1
  mz3_mz2 = mz3 - mz2
  mz4_mz2 = mz4 - mz2
  mz4_mz3 = mz4 - mz3
  
  #Calculate intensity variables
  
  int1_int0 = int1 - int0
  int2_int0 = int2 - int0
  int3_int0 = int3 - int0
  int4_int0 = int4 - int0
  int2_int1 = int2 - int1
  int3_int1 = int3 - int1
  int4_int1 = int4 - int1
  int3_int2 = int3 - int2
  int4_int2 = int4 - int2
  int4_int3 = int4 - int3
  
  #Compile second layer variables
  
  second_layer = c(int1, int2, int3, int4, mz0, mz1_mz0 , mz2_mz0 , mz3_mz0,
                   mz4_mz0, mz2_mz1, mz3_mz1, mz4_mz1, mz3_mz2, mz4_mz2, mz4_mz3, 
                   int1_int0, int2_int0, int3_int0, int4_int0, int2_int1, int3_int1,
                   int4_int1, int3_int2, int4_int2, int4_int3)
  
  return(c(first_layer, second_layer))
  
}

out = data.frame(out)
names = c("int1", "int2", "int3", "int4", "mz0", "mz1_mz0", "mz2_mz0", "mz3_mz0",
          "mz4_mz0", "mz2_mz1", "mz3_mz1", "mz4_mz1", "mz3_mz2", "mz4_mz2", "mz4_mz3",
          "int1_int0", "int2_int0", "int3_int0", "int4_int0", "int2_int1", "int3_int1",
          "int4_int1", "int3_int2", "int4_int2", "int4_int3", "iso_pattern", "notes", "to_predict", "model",
          "cor_1", "cor_2", "cor_3", "cor_4")

colnames(out) = c(names, names[1:25])

#Add iso_pattern to main feature table

ft$iso_pattern = out$iso_pattern

#If did not meet intensity cutoff, no need to consider S

if(any(out$notes == "0")){
  out$to_predict[which(out$notes == "0")] = "FALSE"
}

#################################################################################

#Load models

model1 = readRDS("M+2_S_recog.RDS")
model2 = readRDS("M+3_S_recog.RDS")
model3 = readRDS("M+4_S_recog.RDS")

model1_no = readRDS("M+2_S_number.RDS")
model2_no = readRDS("M+3_S_number.RDS")
model3_no = readRDS("M+4_S_number.RDS")

################################################################################

#S recognition

model1_pred = predict(model1, data.frame(sapply(out[,1:25], as.numeric)))
model1_class = ifelse(model1_pred$predictions[,2] >= 0.5, T, F)
model2_pred = predict(model2, data.frame(sapply(out[,1:25], as.numeric)))
model2_class = ifelse(model2_pred$predictions[,2] >= 0.5, T, F)
model3_pred = predict(model3, data.frame(sapply(out[,1:25], as.numeric)))
model3_class = ifelse(model3_pred$predictions[,2] >= 0.5, T, F)

s_recog = c()
s_recog_prob = c()

for(i in 1:nrow(out)){
  
  if(out$model[i] == "1"){
    s_recog[i] = model1_class[i]
    s_recog_prob[i] = model1_pred$predictions[,2][i]
  } else if(out$model[i] == "2"){
    s_recog[i] = model2_class[i]
    s_recog_prob[i] = model2_pred$predictions[,2][i]
  } else {
    s_recog[i] = model3_class[i]
    s_recog_prob[i] = model3_pred$predictions[,2][i]
  }
  
}

s_recog_prob = round(s_recog_prob,2)

#Clean up S recognition predictions

s_recog[which(out$to_predict == "FALSE")] = "FALSE"
s_recog_prob[which(out$to_predict == "FALSE")] = 0

################################################################################

#S number prediction

model1_pred = predict(model1_no, data.frame(sapply(out[,34:58], as.numeric)))
model1_class = ifelse(model1_pred$predictions[,2] >= 0.5, "2", "1")
model2_pred = predict(model2_no, data.frame(sapply(out[,34:58], as.numeric)))
model2_class = ifelse(model2_pred$predictions[,2] >= 0.5, "2", "1")
model3_pred = predict(model3_no, data.frame(sapply(out[,34:58], as.numeric)))
model3_class = ifelse(model3_pred$predictions[,2] >= 0.5, "2", "1")

s_no = c()
s_no_prob = c()

for(i in 1:nrow(out)){
  
  if(out$model[i] == "1"){
    s_no[i] = model1_class[i]
    s_no_prob[i] = model1_pred$predictions[,2][i]
  } else if(out$model[i] == "2"){
    s_no[i] = model2_class[i]
    s_no_prob[i] = model2_pred$predictions[,2][i]
  } else {
    s_no[i] = model3_class[i]
    s_no_prob[i] = model3_pred$predictions[,2][i]
  }
  
}

#If S was not recognized previously, put 0 in the S number predictions

s_no[which(s_recog == "FALSE")] = 0
s_no_prob = round(as.numeric(s_no_prob),2)
s_no_prob[which(s_recog == "FALSE")] = ""

ft$sulfur = s_recog
ft$sulfur_prob = s_recog_prob
ft$sulfur_no = s_no
ft$multi_sulfur_prob = s_no_prob
ft$cl_flag = ""
ft$cl_flag[which(out$notes == "Cl flag")] = "Cl flag"

#Clean up isotope pattern results

ft$iso_pattern[ft$iso_pattern == 0] = ""

#Clean up backend

stopCluster(cl)

#Save results

if(save){
  write.xlsx(ft, "SulfurFinder_results.xlsx", row.names = F)
}

