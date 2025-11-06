#########################################################
############## CONFORMIST BIAS SNAPSHOT #################
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
# a = normalized weight (1 + b)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs


set.seed(1234)


# PIPELINE -------------------------------------------------------------------
conformist_bias_snapshot <- function(N, mu, c, burnin, timesteps,
                                     p_value_lvl, n_runs) {
  
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
    
    # Define focal variant as the most common
    tab <- table(pop)
    sel_variant_snap <- as.integer(names(tab)[ which.max(tab) ])
    results_snap$variant[run] <- sel_variant_snap
    
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
    
    # Build count matrix
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab) / N # counts to frequencies
    }))
    colnames(freq_mat) <- as.character(unique_variants) # give names
    
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
    c = c,
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
conf_snap_sim <- conformist_bias_snapshot(N = 2000, mu = 0.02, burnin = 1000, c = 0.1,
                                          timesteps = 10, p_value_lvl = 0.05, n_runs = 100)

# Store output across runs in a table
results_table_conf_snapshot <- tibble(
  N  = conf_snap_sim$N,
  "µ" = conf_snap_sim$mu,
  c = conf_snap_sim$c,
  "Burn-in" = conf_snap_sim$burnin,
  "Time steps" = conf_snap_sim$timesteps,
  "α" = conf_snap_sim$p_value_lvl,
  SSR = conf_snap_sim$SSR,
  FNR = conf_snap_sim$FNR,
  proportionNA = conf_snap_sim$proportionNA,
  "Runs" = conf_snap_sim$n_runs
)

print(results_table_conf_snapshot)
write_xlsx(results_table_conf_snapshot, "tables/conf_snap_output/conf_snap_baseline.xlsx")


# ADD THIS TO EXAMINE FREQUENCIES FOR EACH RUN (best when low n_runs)
# Add it just before "# Prepare FIT input"
fv <- freq_mat[, as.character(sel_variant_snap)]
message(sprintf("Run %d: focal counts first 20 gens = %s", run,
                paste(head(fv,20), collapse=",")))