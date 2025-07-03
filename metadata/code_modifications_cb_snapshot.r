sel_variant <- 15 # we choose the focal variant
s_true <- 0.1 # and true coefficient of selection

# apply content biased transmission
w <- ifelse(pop == sel_variant, 1 + s_true, 1) # selected variant weighs 1 + s, neutral = 1
pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities 

# Table to store results for the focal variant 
  results <- tibble(
    run = seq_len(n_runs),
    variant = rep(sel_variant, n_runs),
    fit_p = rep(NA_real_, n_runs),
    inference = rep("absent", n_runs) # default = not detected
  )

this_fit <- fit_results %>% filter(variant == sel_variant) # call the focal variant
    
    if (nrow(this_fit) == 1) {
      p <- this_fit$fit_p # p-value of focal
      sig <- this_fit$sig # either "neutral" or "selection"
      
      # Store across runs
      results$fit_p[run]  <- p
      results$inference[run] <- sig
    }
    # else leave: fit_p = NA, inference = "absent"

# Compute output
  SSR <- mean(results$inference == "selection", na.rm = TRUE)
  FNR <- 1 - SSR

# 25/05/25
# Changes SSR and FNR formulas, excluding NAs i.e. non-tested variants
this_fit <- fit_results %>% 
      filter(variant == sel_variant) # call the focal variant
    
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
      results$fit_p[run]  <- p_value
      results$inference[run] <- inf
  }
  
  # Compute output
  sel <- sum(results$inference == "selection", na.rm = T)
  neut <- sum(results$inference == "neutral", na.rm = T)
  sumNA <- sum(is.na(results$inference))
  
  SSR <- sel / (n_runs - sumNA)
  FNR <- neut / (n_runs - sumNA)
  proportionNA <- sumNA / n_runs
  all_pvals <- unlist(fit_p_count)

# Focal variant problem 26-27/05/25
# SSR unbiased
pop_counts <- table(pop)
    viable <- as.integer(names(pop_counts[pop_counts>=3]))
    if (length(viable)==0) next
    foc_variant_snap <- sample(viable, 1)
	
# SSR modal
pop_counts <- table(pop)
    count_max <- max(pop_counts)
    if (count_max < 3) next
    foc_variant_snap <- as.integer(names(pop_counts[which.max(pop_counts)]))
	
# Hard coding focal variant:
  sel_variant_snap <- 400  # we choose the focal variant
  
  # Table to store results for the focal variant 
  results_snap <- tibble(
    run = seq_len(n_runs),
    variant = rep(sel_variant_snap, n_runs), # number of variant
    fit_p = rep(NA_real_, n_runs), # default NA
    inference = rep(NA_character_, n_runs) # default NA character
  )
  
sel_variant_snap <- N + 1 # reserve this slot

pop[1:inject_cnt] <- sel_variant_snap        
    maxtrait         <- sel_variant_snap
    traitmatrix[1, ] <- pop 
    
    # Observation period after equilibrium
    for (i in 2:timesteps) {
      
      # apply content biased transmission
      w <- ifelse(pop == sel_variant_snap, 1 + b, 1) # normalized weights
      pop <- sample(pop,replace=T,prob = w) # add the weights as probabilities
	  
	  
sample_times <- seq(1, timesteps, by = sample_int)
n_samples    <- length(sample_times)

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
	
# Check results with one parameter setting:
cb_snap_sim <- content_bias_snapshot(N=100, mu=0.01, b=0.5,
                                     burnin=1000, timesteps=1000,
                                     p_value_lvl=0.05, n_runs=100,
                                     sample_int=50)
									 
# Multiple parameter settings:
# Optimal sampling intervals
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
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "Runs" = .x$n_runs,
         "SSR" = round(sim$SSR, 3),
         "FNR" = round(sim$FNR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Sampling interval" = .x$sample_int
  )
})

print(osi_result)

# Sample REAL INTERVALS:
interval_chunks <- list(
  `1-200`   = 1:200,
  `201-400` = 201:400,
  `401-600` = 401:600,
  `601-800` = 601:800,
  `801-1000`= 801:1000
)

chunked_osi_results <- map(interval_chunks, function(chunk_vec) {
  map_dfr(chunk_vec, function(si) {
    sim <- content_bias_snapshot(N=100, mu=0.01, b=0.5,
                                 burnin=1000, timesteps=1000,
                                 p_value_lvl=0.05, n_runs=50,
                                 inject_fvar=20, sample_int=si)
    tibble(N  = sim$N, mu = sim$mu, b = sim$b,
           burnin = sim$burnin, timesteps = sim$timesteps,
           sample_int = si, SSR = sim$SSR, FNR = sim$FNR,
           "%NA" = sim$proportionNA, "Runs" = sim$n_runs)
  })
})  


# PLOTS ----------------------------------------------

# %NA vs SAMPLE_INT
# bind into one big tibble with a `chunk` column
chunked_results_df <- imap_dfr(chunked_results, ~ .x %>% mutate(chunk = .y))

# plot %NA vs. sample_int
ggplot(chunked_results_df, aes(x = sample_int, y = pct_NA)) +
  geom_line(alpha = 0.2) +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, colour = "steelblue") +
  labs(
    x     = "Sampling Interval (time steps)",
    y     = "% of runs with NA",
    title = "%NA vs. Sampling Interval (smoothed)"
  ) +
  theme_minimal(base_size = 16) +   # increase the base text size
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18)
  )
  
# compute minimal inject count to ensure focal survives >= 3 time points
  inject_fvar <- ceiling(N * (1 - p_extinct^(1 / (2 * N))))
  
# Inject variant using the formula, otherwise inject it calling an integer in the function argument
if (is.null(inject_override)) {
    inject_fvar <- ceiling(N * (1 - p_extinct^(1 / (2 * N))))
  } else {
    inject_fvar <- inject_override
  }
  message(sprintf("[Diagnostics] inject_fvar = %d (override = %s)",
                  inject_fvar,
                  ifelse(is.null(inject_override), "false", "true")))
				  
# From frequencies to raw counts per generation
as.numeric(tab)/N # Convert to frequencies
# to
as.integer(tab)  # raw counts per generation

# Went back to frequencies
# Added a message to examine the frequencies of the focal variants across runs
# DIAGNOSTIC: check focal frequencies
    fv <- freq_mat[, as.character(sel_variant_snap)]
    message(sprintf("Run %d: focal counts first 10 gens = %s", run,
                    paste(head(fv,25), collapse=",")))