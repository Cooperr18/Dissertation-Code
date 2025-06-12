#########################################################
############ CONTENT BIAS TIME AVERAGING ################
#########################################################

# Install signatselect
install.packages("pak")
pak::pkg_install("benmarwick/signatselect")

# Reading packages
pkgs <- c(
  "signatselect","dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl"
)
lapply(pkgs, library, character.only = TRUE)

# Modeled as in Boyd & Richerson (1985, p. 140)

# PARAMETERS -----------------------------------------------------------------
# N = Number of individuals
# mu = innovation rate (between 0 and 1, inclusive)
# b = direct/content bias
# w = normalized weight (1 + b)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs
# time_window = time window size

set.seed(1234)

# PIPELINE -------------------------------------------------------------------

# Without sampling interval
content_bias_ta <- function(N, mu, b, burnin, timesteps, p_value_lvl, 
                                  n_runs, time_window, inject_fvar) {
  
  # Assign selected variant
  sel_variant_ta <- N + 1 # reserve this slot
  
  # Table to store results for the focal variant 
  results_ta <- tibble(
    run = seq_len(n_runs),
    variant = rep(sel_variant_ta, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    
    # Initialize population and matrix
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # record every time step
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
    pop[1:inject_fvar] <- sel_variant_ta        
    maxtrait  <- sel_variant_ta
    traitmatrix[1, ] <- pop 
    
    # Observation period after equilibrium
    for (i in 1:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == sel_variant_ta, 1 + b, 1) # selected variant weighs 1 + s, neutral = 1
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }

        traitmatrix[i, ] <- pop
    }
    
    # Time averaging step
    averaged_rows <- floor(timesteps / time_window) # round down to nearest whole number
    
    # list of time-averaged samples (each is a vector of N * time_window variants)
    averaged_samples <- vector("list", averaged_rows)
    for (j in 1:averaged_rows) {
      start <- (j - 1) * time_window + 1
      end <- j * time_window
      # Flatten all traits in the time window into one vector
      averaged_samples[[j]] <- as.vector(traitmatrix[start:end, ])
    }
    
    unique_variants <- sort(unique(unlist(averaged_samples))) # store unique variants across all bins
    
    # Trait matrix to frequency matrix (row = bins, col = variants)
    freq_mat <- t(sapply(averaged_samples, function(variants) {
      tab <- table(factor(variants, levels = unique_variants))
      as.numeric(tab) / length(variants)  # Proportions relative to N * time_window
    }))
    colnames(freq_mat) <- unique_variants
    
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
    results_ta$fit_p[run]  <- p_value
    results_ta$inference[run] <- inf
  }
  
  # Compute output
  sel <- sum(results_ta$inference == "selection", na.rm = T)
  neut <- sum(results_ta$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results_ta$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)
  
  # Export output
  return(list(
    sel_variant_ta = sel_variant_ta, # variant
    
    # parameters
    N = N,
    mu = mu,
    b = b,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    time_window = time_window,
    inject_fvar = inject_fvar,
    
    # outputs
    all_pvals = all_pvals,
    SSR = SSR,
    sumNA = sumNA,
    FNR = FNR,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    fit_results = fit_results, # global results
    results_ta = results_ta # focal variant results
  )
  )
}

# Initial frequency of the variant (inject_fvar) ------------------------------
fvar_list_ta <- list(
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 5),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 10),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 25),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 50),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 60),
  list(N=100, mu=0.01, b = 0.5, burnin=1000, timesteps=1000, 
       p_value_lvl=0.05, n_runs=50, inject_fvar = 75),
)

# Run and store in a table
fvar_result_ta <- map_dfr(fvar_list_ta, ~ {
  sim <- do.call(content_bias_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         b = .x$b,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "SSR" = round(sim$SSR, 3),
         "FNR" = round(sim$FNR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Initial freq" = .x$inject_fvar,
         "Runs" = .x$n_runs
  )
})

print(fvar_result_ta) 



