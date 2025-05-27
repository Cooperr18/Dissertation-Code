################### CONTENT BIAS SNAPSHOT ###################

pkgs <- c(
  "signatselect","dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl"
)
lapply(pkgs, library, character.only = TRUE)

# Pipeline with tsinfer ------------------------------------------
content_bias_snapshot <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs, sample_size) {
  
  p_sel_vec <- numeric(n_runs) # empty vector to store the p-value of foc_variant_snap
  detected_sel <- logical(n_runs) # true if we called selection
  fit_p_count <- vector("list", n_runs) # store p-values
  s_est_list <- vector("list", n_runs) # store estimated coefficient of selection of the focal variant
  N_est_list <- vector("list", n_runs) # store estimated population size
  
  # Assign selected variant and S
  foc_variant_snap <- 15 # we choose the focal variant
  s_true <- 0.1 # and true coefficient of selection
  
  # how many points max we hand to tsinfer
  max_ts <- 100
  
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
    # *** REPLACE map_dfr with loop to avoid deep stack
    # Prepare FIT input
    variants     <- unique(freq_long_filtered$variant)
    results_list <- vector("list", length(variants))
    
    for (j in seq_along(variants)) {
      v    <- variants[j]
      df   <- freq_long_filtered %>% filter(variant == v)
      
      if (nrow(df) < 3 || any(df$freq <= 0) || any(df$freq >= 1)) {
        warning("Skipping variant ", v, " due to extreme frequencies")
        results_list[[j]] <- tibble(
          variant     = v,
          time_points = nrow(df),
          fit_p       = NA_real_,
          s_est       = NA_real_,
          N_est       = NA_real_
        )
        next
      }
      
      if (nrow(df) > max_ts) {
        idx    <- round(seq(1, nrow(df), length.out = max_ts))
        df_sub <- df[idx, ]
      } else {
        df_sub <- df
      }
      
      if (nrow(df_sub) == 0) {
        warning("Variant ", v, ": no data points after subsampling; skipping")
        results_list[[j]] <- tibble(
          variant     = v,
          time_points = nrow(df),
          fit_p       = NA_real_,
          s_est       = NA_real_,
          N_est       = NA_real_
        )
        next
      }
      
      df_sub$count <- as.numeric(rbinom(nrow(df_sub), sample_size, df_sub$freq))
      nvec         <- rep(sample_size, nrow(df_sub))
      
      # ALWAYS run fit() so res_fit is defined
      res_fit <- tryCatch(
        fit(time = df_sub$time, v = df_sub$freq),
        error = function(e) {
          warning("fit() error for variant ", v, ": ", conditionMessage(e))
          list(fit_p = NA_real_)
        }
      )
      
      # Guard for tsinfer
      if (
        nrow(df_sub) < 3 ||
        any(df_sub$freq <= 0)   || any(df_sub$freq >= 1) ||
        any(df_sub$count <= 0)  || any(df_sub$count >= nvec) ||
        any(!is.finite(df_sub$freq)) ||
        any(!is.finite(df_sub$count)) ||
        any(!is.finite(df_sub$time))
      ) {
        warning("Skipping tsinfer for variant ", v, ": bad inputs")
        s_est <- NA_real_
        N_est <- NA_real_
      } else {
        res_ts <- tryCatch(
          tsinfer(tvec = df_sub$time, bvec = df_sub$count, nvec = nvec),
          error = function(e) {
            warning("tsinfer() error for variant ", v, ": ", conditionMessage(e))
            list(s = NA_real_, N = NA_real_)
          }
        )
        s_est <- res_ts$s
        N_est <- res_ts$N
      } # *** ADDED/UPDATED BRACE: end of tsinfer guard
      
      results_list[[j]] <- tibble(
        variant     = v,
        time_points = nrow(df),
        fit_p       = res_fit$fit_p,
        s_est       = s_est,
        N_est       = N_est
      )
    } # *** ADDED BRACE: end of for (j ...)
    
    fit_results <- bind_rows(results_list) %>%
      mutate(sig = case_when(
        is.na(fit_p)         ~ "NA",
        fit_p <= p_value_lvl ~ "selection",
        TRUE                 ~ "neutral"
      ))
    
    fit_p_count[[run]] <- fit_results$fit_p
    s_est_list[[run]] <- fit_results$s_est
    N_est_list[[run]] <- fit_results$N_est
    
    if (foc_variant_snap %in% fit_results$variant) {
      this_p           <- fit_results$fit_p[fit_results$variant == foc_variant_snap]
      p_sel_vec[run]   <- this_p
      detected_sel[run] <- !is.na(this_p) && (this_p <= p_value_lvl)
    } else {
      p_sel_vec[run]    <- NA
      detected_sel[run] <- FALSE
    }
  } # end of for (run ...)
  
  SSR          <- mean(detected_sel, na.rm = TRUE)
  FNR          <- mean(!detected_sel, na.rm = TRUE)
  all_pvals    <- unlist(fit_p_count)
  sumNA        <- sum(is.na(all_pvals))
  proportionNA <- sumNA / length(all_pvals) * 100
  
  return(list(
    SSR            = SSR,
    FNR            = FNR,
    sumNA          = sumNA,
    proportionNA   = proportionNA,
    p_sel_vec      = p_sel_vec,
    detected_sel   = detected_sel,
    fit_p_count    = fit_p_count,
    s_est_list     = s_est_list,
    N_est_list     = N_est_list,
    N              = N,
    mu             = mu,
    burnin         = burnin,
    timesteps      = timesteps,
    p_value_lvl    = p_value_lvl,
    n_runs         = n_runs,
    all_pvals      = all_pvals
  ))
}


cb_snap_sim <- content_bias_snapshot(N = 100, mu = 0.02, burnin = 1000,
                                     timesteps = 1000, p_value_lvl = 0.05, n_runs = 5,
                                     sample_size = 60)

# Store output across runs in a table
results_table_cb_snapshot <- tibble(
  N  = cb_snap_sim$N,
  mu = cb_snap_sim$mu,
  burnin = cb_snap_sim$burnin,
  timesteps = cb_snap_sim$timesteps,
  p_value_lvl = cb_snap_sim$p_value_lvl,
  n_runs = cb_snap_sim$n_runs,
  mean_accuracy = cb_snap_sim$mean_accuracy,
  proportionNA = cb_snap_sim$proportionNA,
  mean_p_value = mean(cb_snap_sim$all_pvals, na.rm = TRUE),
  SSR = cb_snap_sim$SSR,
  FNR = cb_snap_sim$FNR
)

results_table_neutral_snapshot


# Pipeline without tsinfer ----------------------------
set.seed(1234)

# Pipeline hard coding focal
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

# Check results
cb_snap_sim <- content_bias_snapshot(N = 100, mu = 0.02, s_true = 0.3,
                                     burnin = 100, timesteps = 100, 
                                     p_value_lvl = 0.05, n_runs = 100)

results_snap <- cb_snap_sim$results_snap
results_snap
fit_results <- cb_snap_sim$fit_results

# Check table output
results_table_cb_snapshot <- tibble(
  "Focal Variant" = cb_snap_sim$foc_variant_snap,
  N  = cb_snap_sim$N,
  mu = cb_snap_sim$mu,
  "s" = cb_snap_sim$s_true,
  burnin = cb_snap_sim$burnin,
  timesteps = cb_snap_sim$timesteps,
  "α" = cb_snap_sim$p_value_lvl,
  SSR = cb_snap_sim$SSR,
  FNR = cb_snap_sim$FNR,
  "%NA" = cb_snap_sim$proportionNA,
  n_runs = cb_snap_sim$n_runs
)
results_table_cb_snapshot

# Run multiple parameter settings
params_cb_snapshot <- list(
  list(N=100, mu=0.02, s_true = 0.1, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.02, s_true = 0.25, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.02, s_true = 0.5, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.02, s_true = 0.75, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.02, s_true = 0.9, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100)
)

all_results_cb_snapshot <- map_dfr(params_cb_snapshot, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble("Focal Variant" = sim$sel_variant_snap,
         N  = .x$N,
         mu = .x$mu,
         "S" = .x$s_true,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         α = .x$p_value_lvl,
         SSR = sim$SSR,
         FNR = sim$FNR,
         "%NA" = round(sim$proportionNA, 2)*100,
         n_runs = .x$n_runs,
  )
})

all_results_cb_snapshot 

# Export excel
write_xlsx(all_results_cb_snapshot, "multiple_S_cb_snapshot.xlsx")


sum(fit_results$sig == "neutral")
sum(fit_results$sig == "selection")
sum(fit_results$sig == "NA")
nrow(fit_results)
