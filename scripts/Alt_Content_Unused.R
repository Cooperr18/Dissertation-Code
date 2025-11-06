# ALTERNATIVE CONTENT BIAS VERSIONS (UNUSED)


# Pipeline hard-coding focal
content_bias_snapshot <- function(N, mu, s_true, burnin, timesteps, p_value_lvl, n_runs) {
  
  # Assign selected variant
  sel_variant_snap <- 400  # we choose the focal variant
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = rep(sel_variant_snap, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop, replace = T) # we keep it neutral, as we don't count it
      # Add innovations
      innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
    }
    
    # Observation period after equilibrium
    for (i in 1:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == sel_variant_snap, 1 + s_true, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      traitmatrix[i,] <- pop #record the variants of each individual
    }
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to long format
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- unique_variants # give names
    
    # Prepare FIT input
    freq_long <- as.data.frame(freq_mat) %>%
      mutate(time = 1:timesteps) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(freq > 0) %>% # remove zeros
      mutate(variant = as.integer(variant))
    
    # Filter data to 3 time points
    freq_long_filtered <- freq_long %>%
      group_by(variant) %>%
      filter(n_distinct(time) >= 3) %>%
      ungroup()
    
    # FIT application and storing
    fit_results <- freq_long_filtered %>%
      group_split(variant) %>%
      map_dfr(~ {
        df <- as.data.frame(.x)  # .x is a dataframe
        
        # check if we have enough data points
        if (nrow(df) < 3) {
          return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
        }
        
        # safely apply FIT
        res <- tryCatch(
          fit(time = df$time, v = df$freq),
          error = function(e) list(fit_p = NA)  # handle errors
        )
        
        data.frame(
          variant = df$variant[1],
          time_points = nrow(df),
          fit_p = res$fit_p,
          stringsAsFactors = FALSE
        )
      }) %>%
      mutate(
        sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"), # neutral if >0.05
        sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
      )
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    this_fit <- fit_results %>% 
      filter(variant == sel_variant_ta) # call the focal variant
    
    # Avoid length-zero 
    p_value <- this_fit$fit_p[1] 
    if (is.null(p_value) || length(p_value) == 0) {
      p_value <- NA_real_
    }
    
    # assign inference
    inf <- case_when(
      is.na(p_value) ~ NA_character_,
      p_value > 0.05 ~ "neutral",
      TRUE ~ "selection" # if non of the previous, selection
    ) 
    
    # Store both across runs
    results_snap$fit_p[run]  <- p_value
    results_snap$inference[run] <- inf
  }
  
  # Compute output
  sel <- sum(results_snap$inference == "selection", na.rm = T)
  neut <- sum(results_snap$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results_snap$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)
  
  # Export output
  return(list(
    sel_variant_snap = sel_variant_snap, # variant
    
    # parameters
    N = N,
    mu = mu,
    s_true = s_true,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    
    # outputs
    all_pvals = all_pvals,
    SSR = SSR,
    sumNA = sumNA,
    FNR = FNR,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    fit_results = fit_results, # global results
    results_snap = results_snap # focal variant results
    )
  )
}


# Pipeline modal sampling
content_bias_snapshot <- function(N, mu, s_true, burnin, timesteps, p_value_lvl, n_runs) {
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = rep(NA_integer_, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop, replace = T) # we keep it neutral, as we don't record it
      # Add innovations
      innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
    }
    
    # Record the focal variant as the modal at equilibrium
    pop_counts <- table(pop)
    count_max <- max(pop_counts)
    if (count_max < 3) next
    foc_variant_snap <- as.integer(names(pop_counts[which.max(pop_counts)]))
    
    # Store the focal variant
    results_snap$variant[run] <- foc_variant_snap
    
    # Observation period after equilibrium
    for (i in 1:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == foc_variant_snap, 1 + s_true, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      traitmatrix[i,] <- pop #record the variants of each individual
    }
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to long format
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- unique_variants # give names
    
    # Prepare FIT input
    freq_long <- as.data.frame(freq_mat) %>%
      mutate(time = 1:timesteps) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(freq > 0) %>% # remove zeros
      mutate(variant = as.integer(variant))
    
    # Filter data to 3 time points
    freq_long_filtered <- freq_long %>%
      group_by(variant) %>%
      filter(n_distinct(time) >= 3) %>%
      ungroup()
    
    # FIT application and storing
    fit_results <- freq_long_filtered %>%
      group_split(variant) %>%
      map_dfr(~ {
        df <- as.data.frame(.x)  # .x is a dataframe
        
        # check if we have enough data points
        if (nrow(df) < 3) {
          return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
        }
        
        # safely apply FIT
        res <- tryCatch(
          fit(time = df$time, v = df$freq),
          error = function(e) list(fit_p = NA)  # handle errors
        )
        
        data.frame(
          variant = df$variant[1],
          time_points = nrow(df),
          fit_p = res$fit_p,
          stringsAsFactors = FALSE
        )
      }) %>%
      mutate(
        sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"), # neutral if >0.05
        sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
      )
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    this_fit <- fit_results %>% 
      filter(variant == foc_variant_snap) # call the focal variant
    
    # Avoid length-zero 
    p_value <- this_fit$fit_p[1] 
    if (is.null(p_value) || length(p_value) == 0) {
      p_value <- NA_real_
    }
    
    # assign inference
    inf <- case_when(
      is.na(p_value) ~ NA_character_,
      p_value > 0.05 ~ "neutral",
      TRUE ~ "selection" # if non of the previous, selection
    ) 
    
    # Store both across runs
    results_snap$fit_p[run]  <- p_value
    results_snap$inference[run] <- inf
  }
  
  # Compute output
  sel <- sum(results_snap$inference == "selection", na.rm = T)
  neut <- sum(results_snap$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results_snap$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)
  
  # Export output
  return(list(
    foc_variant_snap = foc_variant_snap, # variant
    
    # parameters
    N = N,
    mu = mu,
    s_true = s_true,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    
    # outputs
    all_pvals = all_pvals,
    SSR = SSR,
    sumNA = sumNA,
    FNR = FNR,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    fit_results = fit_results, # global results
    results_snap = results_snap # focal variant results
    )
  )
}


# Pipeline unbiased sampling
content_bias_snapshot <- function(N, mu, s_true, burnin, timesteps, p_value_lvl, n_runs) {
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = rep(NA_integer_, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  # Initialize population
  for (run in 1:n_runs) {
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop, replace = T) # we keep it neutral, as we don't record it
      # Add innovations
      innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
    }
    
    # Record the focal variant sampling randomly >= 3 time points
    pop_counts <- table(pop)
    viable <- as.integer(names(pop_counts[pop_counts>=3]))
    if (length(viable)==0) next
    foc_variant_snap <- sample(viable, 1)
    
    # Store the focal variant
    results_snap$variant[run] <- foc_variant_snap
    
    # Observation period after equilibrium
    for (i in 1:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == foc_variant_snap, 1 + s_true, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      traitmatrix[i,] <- pop #record the variants of each individual
    }
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to long format
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- unique_variants # give names
    
    # Prepare FIT input
    freq_long <- as.data.frame(freq_mat) %>%
      mutate(time = 1:timesteps) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(freq > 0) %>% # remove zeros
      mutate(variant = as.integer(variant))
    
    # Filter data to 3 time points
    freq_long_filtered <- freq_long %>%
      group_by(variant) %>%
      filter(n_distinct(time) >= 3) %>%
      ungroup()
    
    # FIT application and storing
    fit_results <- freq_long_filtered %>%
      group_split(variant) %>%
      map_dfr(~ {
        df <- as.data.frame(.x)  # .x is a dataframe
        
        # check if we have enough data points
        if (nrow(df) < 3) {
          return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
        }
        
        # safely apply FIT
        res <- tryCatch(
          fit(time = df$time, v = df$freq),
          error = function(e) list(fit_p = NA)  # handle errors
        )
        
        data.frame(
          variant = df$variant[1],
          time_points = nrow(df),
          fit_p = res$fit_p,
          stringsAsFactors = FALSE
        )
      }) %>%
      mutate(
        sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"), # neutral if >0.05
        sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
      )
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    this_fit <- fit_results %>% 
      filter(variant == foc_variant_snap) # call the focal variant
    
    # Avoid length-zero 
    p_value <- this_fit$fit_p[1] 
    if (is.null(p_value) || length(p_value) == 0) {
      p_value <- NA_real_
    }
    
    # assign inference
    inf <- case_when(
      is.na(p_value) ~ NA_character_,
      p_value > 0.05 ~ "neutral",
      TRUE ~ "selection" # if non of the previous, selection
    ) 
    
    # Store both across runs
    results_snap$fit_p[run]  <- p_value
    results_snap$inference[run] <- inf
  }
  
  # Compute output
  sel <- sum(results_snap$inference == "selection", na.rm = T)
  neut <- sum(results_snap$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results_snap$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)
  
  # Export output
  return(list(
    foc_variant_snap = foc_variant_snap, # focal variant
    
    # parameters
    N = N,
    mu = mu,
    s_true = s_true,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    
    # outputs
    all_pvals = all_pvals,
    SSR = SSR,
    sumNA = sumNA,
    FNR = FNR,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    fit_results = fit_results, # global results
    results_snap = results_snap # focal variant results
    )
  )
}
