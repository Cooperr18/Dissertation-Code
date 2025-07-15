#########################################################
################ CONTENT BIAS INTERVALS #################
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


# PARAMETERS -----------------------------------------------------------------
# N = Number of individuals
# mu = innovation rate (numeric, between 0 and 1, inclusive)
# b = direct/content bias
# w = normalized weight (1 + d)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs
# sample_int = sample interval

set.seed(1234)

# PIPELINE -------------------------------------------------------------------
content_bias_snapshot <- function(N, mu, b, burnin, timesteps, p_value_lvl, 
                                  n_runs, inject_fvar, sample_int) {
  
  # compute generations recorded
  sample_times <- seq(1, timesteps, by = sample_int)
  n_samples <- length(sample_times)
  
  # Assign selected variant
  sel_variant_snap <- N + 1 # reserve this slot
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = rep(sel_variant_snap, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    
    # Initialize population and matrix
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=n_samples,ncol=N) # Only sampled times, not whole sequence
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
    
    # Inject focal variant in t = 1
    pop[1:inject_fvar] <- sel_variant_snap        
    maxtrait  <- sel_variant_snap
    traitmatrix[1, ] <- pop 
    
    # Observation period after equilibrium
    sample_row <- 1 # start of sampling intervals
    
    for (i in 2:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == sel_variant_snap, 1 + b, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      if (i %in% sample_times) {
        sample_row  <- sample_row + 1
        traitmatrix[sample_row, ] <- pop
      }
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
      mutate(time = sample_times) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(freq > 0) %>% # remove zeros
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
            return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
          }
          
          # safely apply FIT
          res <- tryCatch(
            fit(time = df$time, v = df$freq),
            error = function(e) list(fit_p = NA)
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
    b = b,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    inject_fvar = inject_fvar,
    sample_int = sample_int,
    
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

# Optimal sampling intervals
interval_chunks <- list(
  `1-200`   = 1:200,
  `201-400` = 201:400,
  `401-600` = 401:600,
  `601-800` = 601:800,
  `801-1000`= 801:1000
)

chunked_osi_results <- map(interval_chunks, function(chunk_vec) {
  map_dfr(chunk_vec, function(si) {
    sim <- content_bias_sampl_int(N=100, mu=0.01, b=0.5,
                                  burnin=1000, timesteps=1000,
                                  p_value_lvl=0.05, n_runs=50,
                                  inject_fvar=20, sample_int=si)
    tibble(N  = sim$N, mu = sim$mu, b = sim$b,
           burnin = sim$burnin, timesteps = sim$timesteps,
           sample_int = si, SSR = sim$SSR, FNR = sim$FNR,
           "%NA" = sim$proportionNA, "Runs" = sim$n_runs)
  })
})

print(chunked_osi_results)

osi_list <- list(
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, sample_int = 1:200),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, sample_int = 201:400),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, sample_int = 401:600),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, sample_int = 601:800),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, sample_int = 801:1000)
)

# Run simulation for each sampling interval
osi_result <- map_dfr(osi_list, ~ {
  args <- .x
  args$sample_int <- as.integer(args$sample_int)
  sim <- do.call(content_bias_snapshot, args = args)
  tibble(N  = .x$N,
         mu = .x$mu,
         b = .x$b,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "SSR" = round(sim$SSR, 3),
         "FNR" = round(sim$FNR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Sampling interval" = .x$sample_int,
         "Runs" = .x$n_runs
  )
})

print(osi_result)
