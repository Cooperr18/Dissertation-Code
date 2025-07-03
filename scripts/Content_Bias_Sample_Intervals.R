#########################################################
############ CONTENT BIAS SAMPLE INTERVALS ##############
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
# mu = innovation rate (between 0 and 1, inclusive)
# b = direct/content bias
# a = normalized weight (1 + b)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs
# sample_int = sample interval (L)

set.seed(1234)

# Pipeline ------------------------------------------------------------------
content_bias_sampl_int <- function(N, mu, b, burnin, timesteps,
                                  p_value_lvl, n_runs, sample_int,
                                  p_extinct = 0.05) {
  
  # compute generations recorded
  sample_times <- seq(1, timesteps, by = sample_int)
  n_samples <- length(sample_times)
  
  # compute minimal inject count to ensure focal survives >= 3 time points
  inject_fvar <- ceiling(N * (1 - p_extinct^(1 / (2 * N))))
  
  # Table to store results for the focal variant 
  results_si <- tibble(
    run = seq_len(n_runs),
    variant = NA_integer_, # number of variant
    fit_p = NA_real_, # default NA
    inference = NA_character_ # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    
    # Initialize population and matrix
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=n_samples,ncol=N) # Only sampled times, not whole sequence
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage (neutral transmission + innovation)
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
    sel_variant_si <- maxtrait + 1
    results_si$variant[run] <- sel_variant_si
    
    # Inject focal variant in t = 1
    pop[1:inject_fvar] <- sel_variant_si        
    maxtrait  <- sel_variant_si
    traitmatrix[1, ] <- pop 
    
    # Observation period after equilibrium
    sample_row <- 1 # start of sampling intervals
    
    for (i in 2:timesteps) {
      
      # apply content biased transmission
      a <- ifelse(pop == sel_variant_si, 1 + b, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = a) # add the weights as probabilities
      
      # Random innovation
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
        variant = sel_variant_si,
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
      filter(variant == sel_variant_si) # call the focal variant
    
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
    results_si$fit_p[run]  <- p_value
    results_si$inference[run] <- inf
  }
  
  # Summarise output
  sel <- sum(results_si$inference == "selection", na.rm = T)
  neut <- sum(results_si$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results_si$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)
  
  # Export output
  return(list(
    sel_variant_si = sel_variant_si, # variant
    
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
    results_si = results_si # focal variant results
    )
  )
}

# DEFINE SAMPLING INTERVALS
interval_chunks <- list(
  "1-200"   = 1:200,
  "201-400" = 201:400,
  "401-600" = 401:600,
  "601-800" = 601:800,
  "801-1000"= 801:1000
)

# RUN 
chunked_osi_results <- map(interval_chunks, function(chunk_vec) {
  map_dfr(chunk_vec, function(si) {
    sim <- content_bias_sampl_int(N=100, mu=0.01, b=0.5,
                                  burnin=1000, timesteps=1000,
                                  p_value_lvl=0.05, n_runs=5,
                                  sample_int=si)
    tibble(N  = sim$N, mu = sim$mu, b = sim$b,
           burnin = sim$burnin, timesteps = sim$timesteps,
           sample_int = si, SSR = sim$SSR, FNR = sim$FNR,
           "%NA" = sim$proportionNA, "Runs" = sim$n_runs)
  })
})

print(chunked_osi_results)