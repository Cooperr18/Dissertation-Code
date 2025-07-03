#########################################################
################ CONTENT BIAS SNAPSHOT ##################
#########################################################

# Reading packages
pkgs <- c(
  "pak","dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl"
)
lapply(pkgs, library, character.only = TRUE)

# Install signatselect
pak::pkg_install("benmarwick/signatselect")
library(signatselect)

# Modeled as in Boyd & Richerson (1985, p. 140)

# PIPELINE -------------------------------------------------------------------
content_bias_snapshot <- function(N, mu, b, burnin, timesteps,
                                  p_value_lvl, n_runs,
                                  p_extinct = 0.05, inject_override=NULL) {
  
  # compute minimal inject count to ensure focal survives >= 3 time points
  if (is.null(inject_override)) {
    inject_fvar <- ceiling(N * (1 - p_extinct^(1 / (2 * N))))
  } else {
    inject_fvar <- inject_override
  }
  message(sprintf("[Diagnostics] inject_fvar = %d (override=%s)",
                  inject_fvar, ifelse(is.null(inject_override), "F", "T")))
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = NA_integer_, # number of variant
    fit_p = NA_real_, # default NA
    inference = NA_character_ # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    
    # Initialize population and matrix
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Only sampled times, not whole sequence
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage (neutral transmission and innovation)
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop, replace = T) # keep it neutral, as we don't count it
      # Add innovations
      innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
    }
    
    # Define focal variant label after burn-in to avoid innovation collision
    sel_variant_snap <- maxtrait + 1
    results_snap$variant[run] <- sel_variant_snap
    
    # Inject focal variant in t = 1
    pop[1:inject_fvar] <- sel_variant_snap        
    maxtrait  <- sel_variant_snap
    traitmatrix[1, ] <- pop 
    
    # Observation period after equilibrium
    for (gen in 2:timesteps) {
      
      # apply content biased transmission
      a <- ifelse(pop == sel_variant_snap, 1 + b, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = a) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      traitmatrix[gen, ] <- pop
    }
    
    # Build count matrix
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab) / N # counts to frequencies
    }))
    colnames(freq_mat) <- as.character(unique_variants) # give names
    
    # DIAGNOSTIC: check focal frequencies
    
    # Prepare FIT input
    freq_long <- as_tibble(freq_mat) %>%
      mutate(time = seq_len(nrow(freq_mat))) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(variant == sel_variant_snap | freq > 0) %>%
      mutate(variant = as.integer(variant))
    
    # Filter data to 3 time points
    freq_long_filtered <- freq_long %>%
      group_by(variant) %>%
      filter(n_distinct(time) >= 3) %>%
      ungroup()
    
    # If it's empty (0 time points) return NA
    if (nrow(freq_long_filtered) == 0) {
      fit_results <- tibble(
        variant = sel_variant_snap,
        time_points = 0L,
        fit_p = NA_real_,
        sig = "NA"
      )
    } else {
      
      # Apply FIT and store
      fit_results <- freq_long_filtered %>%
        group_split(variant) %>%
        map_dfr(~ {
          df <- as.data.frame(.x)  # .x is a dataframe
          
          # check if we have enough data points
          if (nrow(df) < 3) {
            return(data.frame(variant = df$variant[1], 
                              time_points = nrow(df), 
                              fit_p = NA_real_, 
                              stringsAsFactors = FALSE))
          }
          
          # safely apply FIT
          res <- tryCatch(
            fit(time = df$time, v = df$freq),
            error = function(e) list(fit_p = NA_real_)
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
    }
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    this_fit <- fit_results %>% 
      filter(variant == sel_variant_snap) # call the focal variant
    
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
  
  # Summarise output
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
    b = b,
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

# BASELINE SIMULATION
cb_snap_sim <- content_bias_snapshot(N = 100, mu = 0.02, burnin = 500, b = 0.25,
                                     timesteps = 100, p_value_lvl = 0.05, n_runs = 100)

# Store output across runs in a table
results_table_cb_snapshot <- tibble(
  N  = cb_snap_sim$N,
  "µ" = cb_snap_sim$mu,
  b = cb_snap_sim$b,
  "Burn-in" = cb_snap_sim$burnin,
  "Time steps" = cb_snap_sim$timesteps,
  "α" = cb_snap_sim$p_value_lvl,
  SSR = cb_snap_sim$SSR,
  FNR = cb_snap_sim$FNR,
  proportionNA = cb_snap_sim$proportionNA,
  "Runs" = cb_snap_sim$n_runs
)

print(results_table_cb_snapshot)


# PARAMETER SWEEP ------------------------------------------------------------------

# Innovation rate ------------------------------------------------------------------
cb_snap_mu_params <- list(
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.025, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.05, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.075, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.1, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.125, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.15, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.175, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.2, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100)
)

# Population size ------------------------------------------------------------------
cb_snap_N_params <- list(
  list(N=10, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=50, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=150, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=200, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=250, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=300, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=350, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=400, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100)
)

# Content bias
cb_snap_b_params <- list(
  list(N=100, mu=0.01, b = 0, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.05, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.1, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.2, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.25, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.3, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.35, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.4, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100)
)

# Time series ----------------------------------------------------------------------
cb_snap_time_params <- list(
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=10, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=50, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=150, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=200, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=250, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=300, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=350, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=400, p_value_lvl=0.05, n_runs=100)
)

# Runs ----------------------------------------------------------------------
cb_snap_time_params <- list(
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=50),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=500),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=5000),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=10000),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=15000),
  list(N=100, mu=0.01, b = 0.15, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=100000)
)


# Run and store --------------------------------------------------------------------

# INNOVATION RATE (MU)
cb_snap_mu_results <- map_dfr(cb_snap_mu_params, ~ {
  sim <- do.call(neutral_snaphot, args = .x)
  tibble(N  = .x$N,
         "µ" = .x$mu,
         b = .x$b,
         "Burn-in" = .x$burnin,
         "Time steps" = .xm$timesteps,
         "α" = .x$p_value_lvl,
         SSR = sim$SSR,
         FNR = sim$FNR,
         proportionNA = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})

# POPULATION SIZE (N)
cb_snap_N_results <- map_dfr(cb_snap_N_params, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble(N  = .x$N,
         "µ" = .x$mu,
         b = .x$b,
         "Burn-in" = .x$burnin,
         "Time steps" = .xm$timesteps,
         "α" = .x$p_value_lvl,
         SSR = sim$SSR,
         FNR = sim$FNR,
         proportionNA = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})

# CONTENT BIAS (b)
cb_snap_N_results <- map_dfr(cb_snap_b_params, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble(N  = .x$N,
         "µ" = .x$mu,
         b = .x$b,
         "Burn-in" = .x$burnin,
         "Time steps" = .xm$timesteps,
         "α" = .x$p_value_lvl,
         SSR = sim$SSR,
         FNR = sim$FNR,
         proportionNA = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})

# POPULATION SIZE (N)
cb_snap_N_results <- map_dfr(cb_snap_b_params, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble(N  = .x$N,
         "µ" = .x$mu,
         b = .x$b,
         "Burn-in" = .x$burnin,
         "Time steps" = .xm$timesteps,
         "α" = .x$p_value_lvl,
         SSR = sim$SSR,
         FNR = sim$FNR,
         proportionNA = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})

# TIME SERIES
cb_snap_time_results <- map_dfr(cb_snap_time_params, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         b = .x$b,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = sim$mean_accuracy,
         "FPR" = 1 - sim$mean_accuracy,
         "%NA" = sim$proportionNA,
         "Runs" = .x$n_runs,
  )
})

# PRINT RESULTS
print(cb_snap_mu_results)
print(cb_snap_N_results)
print(cb_snap_b_results)
print(cb_snap_time_results)

# EXPORT TO SPREADSHEET
write_xlsx(cb_snap_mu_results, "tables/cb_snap_output/cb_snap_mu_params.xlsx") # mu

write_xlsx(cb_snap_N_results, "tables/cb_snap_output/cb_snap_N_params.xlsx") # N

write_xlsx(cb_snap_b_results, "tables/cb_snap_output/cb_snap_b_params.xlsx") # b

write_xlsx(cb_snap_time_params, "tables/cb_snap_output/cb_snap_time_params.xlsx") # time series
