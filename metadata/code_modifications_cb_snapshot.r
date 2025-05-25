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