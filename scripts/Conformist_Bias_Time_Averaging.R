#########################################################
########## CONFORMIST BIAS TIME AVERAGING ###############
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

# Modeled as in Boyd & Richerson (1985)

# PARAMETERS -----------------------------------------------------------------
# N = Number of individuals
# mu = innovation rate (between 0 and 1, inclusive)
# c = conformist bias
# a = normalized weight (1 + c)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs
# time_window = time window size

set.seed(1234)

# PIPELINE ---------------------------------------------------------------------
conformist_bias_ta <- function(N, mu, c, burnin, timesteps, 
                               p_value_lvl, n_runs, time_window) {
  
  # Table to store results for the focal variant 
  results_ta <- tibble(
    run = seq_len(n_runs),
    variant = NA_integer_, # number of variant
    fit_p = NA_real_, # default NA
    inference = NA_character_ # default NA character
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
    
    # Define focal variant label after burn-in to avoid innovation collision
    tab <- table(pop)
    sel_variant_ta <- as.integer(names(tab)[ which.max(tab) ])
    results_ta$variant[run] <- sel_variant_ta
    
    # record t = 1
    traitmatrix[1, ] <- pop  
    
    # Observation period after equilibrium
    for (gen in 2:timesteps) {
      
      # conformist biased transmission
      counts <- tabulate(pop, nbins = maxtrait) # get counts for all variants
      w_var <- counts^(1 + c) # assign normalized weight to all variants:
      a <- w_var[pop] # give each individual its variant's weight
      # sample next generation
      pop <- sample(pop,N,replace=T,prob = a) # add the weights as probabilities
      
      # Neutral innovation
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      traitmatrix[gen, ] <- pop
    }
    
    # Time averaging step
    averaged_rows <- floor(timesteps / time_window) # round down to nearest whole number
    
    # list of time-averaged samples (each is a vector of N * time_window variants)
    averaged_samples <- vector("list", averaged_rows)
    if (averaged_rows >= 1) {
      for (j in seq_len(averaged_rows)) {
        start <- (j - 1) * time_window + 1
        end <- j * time_window
        # Flatten all traits in the time window into one vector
        averaged_samples[[j]] <- as.vector(traitmatrix[start:end, ])
      }
      
      # Build frequency matrix
      unique_variants <- sort(unique(unlist(averaged_samples))) # store unique variants across all bins
      freq_mat <- t(sapply(averaged_samples, function(variants) {
        tab <- table(factor(variants, levels = unique_variants))
        as.numeric(tab) / length(variants)  # Proportions relative to pool of variants
      }))
      colnames(freq_mat) <- unique_variants
      
      # Prepare FIT input
      freq_long <- as.data.frame(freq_mat) %>%
        mutate(time = 1:nrow(.)) %>%
        pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
        filter(freq > 0) %>% # remove zeros
        mutate(variant = as.integer(variant))
      
      # Filter data to 3 time points
      freq_long_filtered <- freq_long %>%
        group_by(variant) %>%
        filter(n_distinct(time) >= 3) %>%
        ungroup()
      
      if (nrow(freq_long_filtered) == 0) {
        # No data at all -> one row of defaults
        fit_results <- tibble(
          variant = sel_variant_ta,
          time_points = 0L,
          fit_p = NA_real_,
          sig = "NA"
        )
      } else {
        
        # Apply FIT and store
        fit_results <- freq_long_filtered %>%
          group_split(variant) %>%
          map_dfr(function(df) {
            pval <- tryCatch(
              fit(time = df$time, v = df$freq)$fit_p,
              error = function(e) NA_real_
            )
            # classify immediately
            sig_flag <- if (is.na(pval)) {
              "NA"
            } else if (pval > p_value_lvl) {
              "neutral"
            } else {
              "selection"
            }
            tibble(
              variant = df$variant[1],
              time_points = nrow(df),
              fit_p = pval,
              sig = sig_flag)
          })
      }  
    } else {
      # no full windows → skip all of the above
      fit_results <- tibble(
        variant = sel_variant_ta,
        time_points = 0,
        fit_p = NA_real_,
        sig = "NA")
    }
    
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
    sel_variant_ta = sel_variant_ta, # focal variant
    
    # parameters
    N = N,
    mu = mu,
    c = c,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    time_window = time_window,
    
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


# BASELINE SIMULATION
conf_ta_sim <- conformist_bias_ta(N = 750, mu = 0.02, burnin = 1000, c = 0.1,
                             timesteps = 100, p_value_lvl = 0.05,n_runs = 100,
                             time_window = 10)

# Store output across runs in a table
results_table_conf_ta <- tibble(
  N  = conf_ta_sim$N,
  "µ" = conf_ta_sim$mu,
  c = conf_ta_sim$c,
  "Burn-in" = conf_ta_sim$burnin,
  "Time steps" = conf_ta_sim$timesteps,
  "α" = conf_ta_sim$p_value_lvl,
  "w" = conf_ta_sim$time_window,
  SSR = conf_ta_sim$SSR,
  FNR = conf_ta_sim$FNR,
  proportionNA = conf_ta_sim$proportionNA,
  "Runs" = conf_ta_sim$n_runs
)
print(results_table_conf_ta)
write_xlsx(results_table_conf_ta, "conf_ta_output/conf_ta_baseline.xlsx")



# ADD THIS TO EXAMINE FREQUENCIES FOR EACH RUN (best when low n_runs)
# Add it just before "# Prepare FIT input"
fv <- freq_mat[, as.character(sel_variant_ta)]
message(sprintf("Run %d: focal counts first 20 gens = %s", run,
                paste(head(fv,20), collapse=",")))