#########################################################
############ MODEL SIMULATIONS & PLOTS ##################
#########################################################


# INITIALIZE PACKAGES

pkgs <- c(
  "dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl", "tidyverse",
  "viridis", "stringr", "viridisLite",
  "readxl", "patchwork", "scales"
)
lapply(pkgs, library, character.only = TRUE)
library(pak)

pak::pkg_install("benmarwick/signatselect")
pak::pkg_install("benmarwick/evoarchdata")
library(signatselect)

set.seed(1234)

# FIGURE 1 ---------------------------------------------------------------------

# BORROWED FROM https://github.com/benmarwick/signatselect

# Data setup
library(evoarchdata)
library(signatselect)
data("ceramics_lbk_merzbach")
ceramics_lbk_merzbach


# Decoration types are ordered numerically and plotted each one of them. 
decoration_types <- 
  names(ceramics_lbk_merzbach)[-1] %>%
  enframe() %>% 
  separate(value, into = c('a', 'b'), 2) %>% 
  mutate(b = parse_number(b)) %>% 
  arrange(b) %>% 
  unite(decorations, c(a,b), sep = "") %>% 
  pull(decorations)

# see how the freqs of each change over time
# gather reorganizes the data to a long format with its columns, 
# and mutate and fct_relevel indicates the order in which you want them to be showed when you plot them. 
ceramics_lbk_merzbach_long <-
  ceramics_lbk_merzbach %>%
  gather(variable, value, -Phase) %>% 
  mutate(Phase = fct_relevel(Phase, ceramics_lbk_merzbach$Phase)) %>% 
  mutate(variable = fct_relevel(variable, decoration_types))

max_n <- 50
ceramics_lbk_merzbach_long_subset <-
  ceramics_lbk_merzbach_long %>% 
  group_by(variable) %>% 
  filter(max(value) > max_n)

# keep these decorations
decorations_to_keep <- unique(as.character(ceramics_lbk_merzbach_long_subset$variable))


# reshape data for each decoration type to go into test:
ceramics_lbk_merzbach_prop <- 
  ceramics_lbk_merzbach %>% 
  select(Phase, decorations_to_keep)

df <- ceramics_lbk_merzbach_prop[ , 2:ncol(ceramics_lbk_merzbach_prop)]
time <- utils:::.roman2numeric(ceramics_lbk_merzbach_prop$Phase)

list_of_dfs <- vector("list", ncol(df))
names(list_of_dfs) <- names(df)

for(i in 1:ncol(df)){
  tmp <-
    data.frame(time = time,
               count_this_one = df[[i]],
               count_others = rowSums(df[, (seq_len(ncol(df)))[-i]   ]))
  
  tmp$frequency = with(tmp, count_this_one / count_others)
  
  # collect results and exclude rows with zero counts for this type i
  list_of_dfs[[i]] <- tmp[which(tmp$count_this_one != 0 ), ]
}

# we need a min of three time points to compute the FIT, so drop decoration types with less than 3
list_of_dfs_three_or_more <- 
  keep(list_of_dfs, ~nrow(.x) >= 3)

fit_safely <- 
  safely(fit, 
         otherwise = data.frame(fit_stat = NA,
                                fit_p = NA))

# apply test to each decoration type
df_fit_test_results <-
  list_of_dfs_three_or_more %>%
  bind_rows(.id = "type") %>%
  nest(data = -type) %>%
  mutate(fit_test = map(data,
                        ~fit_safely(time = .x$time,
                                    v =    .x$frequency))) %>%
  mutate(fit_p = map(fit_test, ~.x$result %>% bind_rows)) %>%
  unnest(fit_p) %>%
  mutate(sig = ifelse(fit_p <= 0.05, "selection", "neutral"))

ceramics_lbk_merzbach_long_sig <-
  ceramics_lbk_merzbach_long_subset %>%
  ungroup %>% 
  left_join(df_fit_test_results %>% 
              select(type, sig), by = c("variable" = "type")) %>%
  mutate(Phase_num = utils:::.roman2numeric(as.character(Phase))) %>% 
  mutate(variable = fct_relevel(factor(variable, levels = decoration_types))) %>% 
  arrange(variable, Phase_num)

# Run tsinfer() and prepare the summary table
decoration_types_tsinfer_output <- 
  list_of_dfs_three_or_more %>% 
  purrr::map_dfr(
    ~ tsinfer(
      tvec = .x$time,
      bvec = .x$count_this_one,
      nvec = .x$count_others,
      verbose = FALSE),
    .id = "type") %>% 
  left_join(df_fit_test_results %>%  # join p-value and "sig" tag
              select(type, fit_p, sig), by = "type") %>% 
  rename(p_value = fit_p, inference_tag = sig) # rename for better clarity

# Check results
print(decoration_types_tsinfer_output)

# Add p-value and inference and round to three decimal places
decoration_types_tsinfer_rounded <- 
  decoration_types_tsinfer_output %>%
  transmute(Type = type, `Coefficient of selection (s)` = s,
    `Initial frequency (f0)` = f0, LL = LL, `p-value` = p_value,
    `Population size (α)` = alpha, Inference = inference_tag) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Print table and export excel
print(decoration_types_tsinfer_rounded)
write_xlsx(decoration_types_tsinfer_rounded, "tables/Merzbach_types_ts_output.xlsx")

# RELATIVE FREQUENCIES
df_rel <- ceramics_lbk_merzbach_long_sig_to_plot_with_others %>%
  group_by(time) %>% 
  mutate(
    total_counts = sum(count_this_one, na.rm = TRUE),
    frequency = count_this_one / total_counts
  ) %>% 
  ungroup()

# PLOT INFERENCE BY COLOUR AND VARIANT BY SHAPE
ggplot(df_rel, aes(x = time, y = frequency, group = type,
                   shape = type, # variant different shapes
                   colour = sig )) + # inference by colour
  geom_line(aes(colour = sig), linewidth = 1) + 
  geom_point(size = 4) + # points by variant
  scale_shape_manual(name = "Decoration\nVariant",
                     values = seq_along(unique(df_rel$type))) + # different shapes
  scale_colour_viridis_d(name  = "Inference", begin = 0.25, end = 0.75) +
  guides(shape  = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_minimal() +
  labs(x = "Time series", y = "Relative frequency") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# FIGURE 2 ---------------------------------------------------------------------

# define parameters
neutral_df <- data.frame(freq = seq(1, 100)/100) %>% 
  mutate(probability = freq)

# neutral transmission adoption probability based on Wright-Fisher neutral model
ggplot(neutral_df, aes(x = freq, y = probability)) +
  geom_line(colour = "grey80", size = 1) +
  geom_point(aes(colour = probability), size = 4) +
  scale_colour_viridis_c(guide = "none") +
  coord_equal() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid = element_blank()
  ) +
  labs(x = "Frequency of Variant ", y = "Adoption Probability")


# FIGURES 3a & 3b --------------------------------------------------------------

# FIGURE 3a
neutral_snapshot <- function(N, mu, timesteps, p_value_lvl, n_runs) {
  runs <- vector("list", length = n_runs)
  
  for (run in 1:n_runs) {
    # Initialize population: 33% start with variant 1, others have unique variants
    initial_variant_count <- round(N * 0.33)
    pop <- c(rep(1, initial_variant_count), 2:(N - initial_variant_count + 1))
    
    traitmatrix <- matrix(NA, nrow = timesteps, ncol = N)
    maxtrait <- max(pop)
    
    for (i in 1:timesteps) {
      pop <- sample(pop, replace = TRUE)  # Neutral drift
      innovate <- which(runif(N) < mu)
      if (length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      traitmatrix[i, ] <- pop
    }
    
    # Track ONLY variant 1 (initial frequency = 0.33)
    freq_variant1 <- rowMeans(traitmatrix == 1, na.rm = TRUE)
    
    runs[[run]] <- data.frame(
      time = 1:timesteps,
      freq = freq_variant1,
      run = run,
      N = N,
      mu = mu,
      p_value_lvl = p_value_lvl,
      variant = "Variant 1 (Initial freq = 0.33)"
    )
  }
  bind_rows(runs)
}

n_snap_N_params <- list(
  list(N=10, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=50, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=100, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=250, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=500, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=1000, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4)
)

# Run all simulations and bind into one big tibble
all_results <- map_dfr(n_snap_N_params, ~
                         neutral_snapshot(.x$N, .x$mu, .x$timesteps, 
                                          .x$p_value_lvl, .x$n_runs))

ggplot(all_results, aes(x = time, y = freq, color = factor(run))) +
  geom_line(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0.33, linetype = "dashed", color = "gray30", size = 0.8) +
  facet_wrap(~ N, scales = "fixed", ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Run") +
  labs(
    x = "Time Step",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 20) +  # Base font size increased to 20
  theme(
    axis.title.x = element_text(size = 22, margin = margin(t = 10)),  # Larger X-axis label
    axis.title.y = element_text(size = 22, margin = margin(r = 10)),  # Larger Y-axis label
    axis.text = element_text(size = 18),  # Larger tick labels
    strip.text = element_text(size = 20, face = "bold"),  # Larger facet labels
    legend.title = element_text(size = 20),  # Larger legend title
    legend.text = element_text(size = 18),  # Larger legend labels
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")  # Extra space between facets
  )


# FIGURE 3b
n_snap_init_freq <- function(N, mu, timesteps,n_variants, n_runs) {
  
  all_runs <- vector("list", n_runs)
  ini <- seq_len(n_variants) # Initial cultural variants
  
  # Initial frequencies
  p0 <- 0.6  
  init_prob <- rep((1 - p0) / (n_variants - 1), n_variants)
  init_prob[5] <- p0
  
  for (run in 1:n_runs) {
    traitmatrix <- matrix(NA, nrow = timesteps, ncol = N)
    
    # t = 1 initial frequencies
    pop <- sample(1:n_variants, N, replace = TRUE, prob = init_prob)
    traitmatrix[1, ] <- pop
    
    for(t in 2:timesteps) {
      pop <- sample(pop, N, replace = TRUE) # neutrality
      innovate <- which(runif(N) < mu) # random innovation
      if(length(innovate)) {
        pop[innovate] <- sample(ini, length(innovate), replace = TRUE) 
      }
      traitmatrix[t, ] <- pop
    }
    
    # Trait matrix to frequency matrix
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = ini)) 
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- as.character(ini)
    df_long <- as.data.frame(freq_mat) %>%
      mutate(time = seq_len(nrow(freq_mat)),
             run  = run) %>%
      pivot_longer(cols = -c(time, run),
                   names_to = "variant",
                   values_to = "frequency")
    all_runs[[run]] <- df_long
  }
  bind_rows(all_runs)
}

neutral_sim_init_freq <- n_snap_init_freq(N = 100, mu = 0.01, timesteps = 100, 
                                               n_runs = 1, n_variants = 6)

ggplot(neutral_sim_init_freq, aes(x = time, y = frequency,
                            group  = interaction(variant, run),
                            colour = variant)) +
  geom_line(linewidth = 0.75) +
  ylim(c(0, 1)) +
  scale_colour_viridis_d(name = "Variant") +
  theme_minimal() +
  theme(axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  labs(x = "Time step", y = "Relative frequency") +
  facet_wrap(~ variant, ncol = 2)


# CONTENT BIASED TRANSMISSION
# FIGURE 4 ------------------------------------------------------------
# Build a data.frame of x, y and all b‐values
df <- tibble(x = 1:100) %>%
  mutate(x_scaled = x/100, y = 100 - x) %>%
  crossing(b = seq(0.2, 1.4, by = 0.2)) %>%
  mutate(p = x*(1 + b) / (x*(1 + b) + y))

df_neutral <- tibble(x_scaled = seq(0, 1, length.out = 100), 
                     p = seq(0, 1, length.out = 100))

ggplot(df, aes(x = x_scaled, y = p, color = as.factor(b))) +
  geom_line(linewidth = 1) +
  geom_line(data = df_neutral,
            aes(x = x_scaled, y = p),
            linetype = "dashed",
            color = "black",
            linewidth = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_viridis_d(name = "b") +
  labs(x = "Frequency of Variant",
       y = "Adoption probability (p)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid.minor = element_blank())

# FIGURE 5 ---------------------------------------------------------------------
df_negative <- tibble(x = 1:100) %>%
  mutate(x_scaled = x/100, y = 100 - x) %>%
  crossing(b = seq(-1, -0.2, by = 0.2)) %>%
  mutate(p = x*(1 + b) / ( x*(1 + b) + y ))

# Plot
ggplot(df_negative, aes(x = x_scaled, y = p, color = as.factor(b))) +
  # biased curves
  geom_line(size = 1) +
  # neutral curve
  geom_line(data = df_neutral, aes(x = x_scaled, y = p),
            linetype = "dashed", color = "black", size = 1) +
  # scales & labels
  scale_color_viridis_d(name = "Content bias (b)",
                        option = "D") +
  labs(
    x = "Frequency of Variant",
    y = "Probability of adoption") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position  = "top"
  )

# FIGURE 6 ---------------------------------------------------------------------
dff <- tibble(x = 1:100) %>%
  mutate(x_scaled = x/100, y = 100 - x) %>%
  crossing(b = seq(-1, 1.4, by = 0.2)) %>%
  mutate(p = x^(1 + b) / ( x^(1 + b) + y ))

# Separate neutral for dashed line
df_neutral <- tibble(x_scaled = seq(0, 1, length.out = 100), 
                     p = seq(0, 1, length.out = 100))
# Plot
ggplot(dff, aes(x = x_scaled, y = p, color = as.factor(b))) +
  # biased curves
  geom_line(size = 1) +
  # neutral curve
  geom_line(data = df_neutral, aes(x = x_scaled, y = p),
            linetype = "dashed", color = "black", size = 1) +
  # scales & labels
  scale_color_viridis_d(name = "Conformist bias (c)",
                        option = "D") +
  labs(
    x = "Frequency of variant",
    y = "Probability of adoption") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position  = "top"
  )

# FIGURE 7 ---------------------------------------------------------------------
neutral_snapshot_fig <- function(N, mu, burnin, timesteps, n_variants = 4, n_runs) {
  
  all_runs <- vector("list", n_runs)
  
  for (run in 1:n_runs) {
    ini <- 1:n_variants # initial variants
    traitmatrix <- matrix(NA, nrow = timesteps, ncol = N)
    
    # Create initial population with EXACT equal frequencies
    pop <- rep(ini, each = N/n_variants)[1:N] # Perfect 0.25 each
    
    # Burn-in stage - no recording
    for(i in 1:burnin) {
      # Neutral transmission with frequency correction
      current_freqs <- tabulate(pop, nbins = n_variants)/N
      weights <- (1/current_freqs)[pop] # Inverse frequency weighting
      pop <- sample(pop, size = N, replace = TRUE, prob = weights)
      
      # Innovation with equal probability for all variants
      innovate <- which(runif(N) < mu)
      if(length(innovate) > 0) {
        pop[innovate] <- sample(ini, length(innovate), 
                                replace = TRUE, 
                                prob = rep(0.25, n_variants))
      }
    }
    
    # Observation period after equilibrium
    for (t in 1:timesteps) {
      # Frequency-corrected sampling
      current_freqs <- tabulate(pop, nbins = n_variants)/N
      weights <- (1/current_freqs)[pop]
      pop <- sample(pop, size = N, replace = TRUE, prob = weights)
      
      # Balanced innovation
      innovate <- which(runif(N) < mu)
      if(length(innovate) > 0) {
        pop[innovate] <- sample(ini, length(innovate), 
                                replace = TRUE,
                                prob = rep(0.25, n_variants))
      }
      
      traitmatrix[t,] <- pop # record the variants
    }
    
    # Calculate frequencies
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = ini))
      as.numeric(tab)/N
    }))
    colnames(freq_mat) <- as.character(ini)
    
    # Convert to long format
    df_long <- as.data.frame(freq_mat) %>%
      mutate(time = seq_len(nrow(freq_mat)),
             run = run) %>%
      pivot_longer(cols = -c(time, run),
                   names_to = "variant",
                   values_to = "frequency")
    
    all_runs[[run]] <- df_long
  }
  bind_rows(all_runs)
}


neutral_snap_traj <- neutral_snapshot_fig(N = 100, mu = 0.02, burnin = 200,
                                          timesteps = 200, n_runs = 1)

ggplot(neutral_snap_traj, aes(x = time, y = frequency,
                            group  = interaction(variant, run))) +
  geom_line(linewidth = 0.8) +
  ylim(c(0, 0.4)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 18)) +
  labs(x = "Time step", y = "Frequency") +
  facet_wrap(~ variant, ncol = 2)

# FIGURE 8 ---------------------------------------------------------------------
neutral_ta_fig <- function(N, mu, burnin, timesteps, n_variants = 4, n_runs,
                           time_window = 20) {
  
  all_runs <- vector("list", n_runs)
  
  for (run in 1:n_runs) {
    ini <- 1:n_variants # initial variants
    traitmatrix <- matrix(NA, nrow = timesteps, ncol = N)
    
    # Create initial population with EXACT equal frequencies
    pop <- rep(ini, each = N/n_variants)[1:N] # Perfect 0.25 each
    
    # Burn-in stage - no recording
    for(i in 1:burnin) {
      # Neutral transmission with frequency correction
      current_freqs <- tabulate(pop, nbins = n_variants)/N
      weights <- (1/current_freqs)[pop] # Inverse frequency weighting
      pop <- sample(pop, size = N, replace = TRUE, prob = weights)
      
      # Innovation with equal probability for all variants
      innovate <- which(runif(N) < mu)
      if(length(innovate) > 0) {
        pop[innovate] <- sample(ini, length(innovate), 
                                replace = TRUE, 
                                prob = rep(0.25, n_variants))
      }
    }
    
    # Observation period after equilibrium
    for (t in 2:timesteps) {
      # Frequency-corrected sampling
      current_freqs <- tabulate(pop, nbins = n_variants)/N
      weights <- (1/current_freqs)[pop]
      pop <- sample(pop, size = N, replace = TRUE, prob = weights)
      
      # Balanced innovation
      innovate <- which(runif(N) < mu)
      if(length(innovate) > 0) {
        pop[innovate] <- sample(ini, length(innovate), 
                                replace = TRUE,
                                prob = rep(0.25, n_variants))
      }
      
      traitmatrix[t,] <- pop # record the variants
    }
    
    # Time averaging step
    n_bins <- floor(timesteps / time_window)
    bin_labels <- 1:n_bins  # Integer bin numbers
    
    averaged_samples <- lapply(seq_len(n_bins), function(j) {
      rows <- ((j - 1) * time_window + 1):(j * time_window)
      as.vector(traitmatrix[rows, ])
    })
    
    unique_variants <- sort(unique(unlist(averaged_samples)))
    freq_mat <- t(sapply(averaged_samples, function(x) {
      tab <- table(factor(x, levels = ini))  # Track only original variants
      as.numeric(tab) / (N * time_window)
    }))
    colnames(freq_mat) <- as.character(ini)
    
    # Convert to long format with integer bins
    df_long <- as.data.frame(freq_mat) %>%
      mutate(bin = bin_labels,  # Use integer sequence
             run = run) %>%
      pivot_longer(cols = -c(bin, run),
                   names_to = "variant",
                   values_to = "frequency")
    
    all_runs[[run]] <- df_long
  }
  bind_rows(all_runs)
}

neutral_ta_traj <- neutral_ta_fig(N = 100, mu = 0.02, burnin = 200,
                                  timesteps = 200, 
                                  n_runs = 1, time_window = 20)

ggplot(neutral_ta_traj, aes(x = bin, y = frequency,
                            group = interaction(variant, run))) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = scales::pretty_breaks(),  # Integer breaks
                     labels = scales::number_format(accuracy = 1)) +  # No decimals
  ylim(c(0, 0.4)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  labs(x = "Time bin", y = "Frequency") + 
  facet_wrap(~ variant, ncol = 2)


# FIGURE 9 - BORROWED FROM FEDER ET AL. (2014) ---------------------------------
# FIGURE 10 --------------------------------------------------------------------

# Import results from baseline simulation
baseline_df <- read_excel("tables/n_snap_output/n_snap_baseline.xlsx")
glimpse(baseline_df)

binwidth <- 0.05

# PLOT NDR DISTRIBUTION WHEN NON-NAs ARE COUNTED
ggplot(baseline_df, aes(x = NDR)) +
  geom_histogram(binwidth = 0.0025, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 1 - baseline_df$α[1], linetype   = "dashed", 
             color = "red", linewidth  = 1) +
  annotate("text",x  = 0.95,y = Inf,label = "Expected NDR",
           hjust = -0.1,vjust = 2, size = 6) +
  labs(x = "Neutral Detection Rate",y = "Frequency",
    caption  = paste0(
      "Average = ", round(baseline_df$NDR[1], 3), 
      " | Runs ≥ 95% = ", baseline_df$`>95% runs`[1], "%", 
      " | Number of runs = ", baseline_df$Runs[1], 
      " | % NA = ", round(baseline_df$`%NA`[1], 1), "%")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.caption = element_text(size = 14))

# FIGURE 11 ---------------------------------------------------------------------
# RUN NEUTRAL SNAPSHOT COUNTING NAs WHEN COMPUTING NDR
neutral_snapshot_NA <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs) {
  
  accuracy_snapshot <- numeric(n_runs) # empty vector for accuracy tracking each run 
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop,replace=T)
      
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
      pop <- sample(pop,replace=T)
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      traitmatrix[i,] <- pop #record the variants of each individual
    }
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to frequency matrix
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
    
    # Store metrics for this run
    total_variants <- nrow(fit_results)
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_snapshot[run] <- NDR  # track for each run
  }
  
  # Store results
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_snapshot >= 0.95) / n_runs * 100
  
  all_pvals <- unlist(fit_p_count)
  
  # Proportion of NA from the simulation
  sumNA <- sum(is.na(all_pvals))
  proportionNA <- sumNA / length(all_pvals) * 100
  if(sumNA == 0) {
    proportionNA <- 0
  }
  
  # Store output
  return(list(
    
    # PARAMETERS
    N = N,
    mu = mu,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    all_pvals = all_pvals,
    
    # OUTPUTS
    accuracy_snapshot = accuracy_snapshot,
    mean_accuracy = mean_accuracy,
    high_accuracy_runs = high_accuracy_runs,
    sumNA = sumNA,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count
    )
  )
}

# BASELINE SIMULATION
n_snap_NA <- neutral_snapshot_NA(N = 100, mu = 0.01, burnin = 1000,
                               timesteps = 1000, p_value_lvl = 0.05, n_runs = 100)

results_n_snap_NA <- tibble(
  N  = n_snap_NA$N,
  "µ" = n_snap_NA$mu,
  burnin = n_snap_NA$burnin,
  timesteps = n_snap_NA$timesteps,
  "α" = n_snap_NA$p_value_lvl,
  NDR = n_snap_NA$mean_accuracy,
  FPR = n_snap_NA$mean_FPR,
  ">95% runs" = n_snap_NA$high_accuracy_runs,
  "%NA" = n_snap_NA$proportionNA,
  mean_p_value = mean(n_snap_NA$all_pvals, na.rm = TRUE),
  "Runs" = n_snap_NA$n_runs
)
print(results_n_snap_NA)

plot_n_snap_NA <- function(results_n_snap_NA, binwidth = 0.005) {
  ggplot(data = data.frame(NDR = results_n_snap_NA$accuracy_snapshot), aes(x = NDR)) +
    geom_histogram(aes(y = ..count.. / sum(..count..)),
                   binwidth = binwidth, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
    labs(x = "Neutral Detection Rate", y = "Frequency",
         caption = paste("Average =", round(mean(results_n_snap_NA$accuracy_snapshot), 3), "|", 
                         "Runs ≥ 95% =", round(results_n_snap_NA$high_accuracy_runs, 1), "%", "|",
                         "Number of runs =", results_n_snap_NA$n_runs, "|",
                         "% NA =", round(results_n_snap_NA$proportionNA, 2), "%")) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.caption = element_text(size = 14)
    )
}
plot_n_snap_NA(n_snap_NA)


# FIGURE 12 -----------------------------------------------------------------------
# BUILD DF OF %NA PER RUN
na_by_run <- tibble(
  run = seq_along(n_snap_NA$fit_p_count),
  propNA = map_dbl(n_snap_NA$fit_p_count, ~ mean(is.na(.x)) * 100))

# Define maximum %NA
max_idx <- na_by_run$run[which.max(na_by_run$propNA)]
max_prop <- max(na_by_run$propNA)

# PLOT %NA BY RUN, HIGHLIGHTING HIGHER %NAs
ggplot(na_by_run, aes(x = run, y = propNA)) +
  geom_area(fill = "skyblue", alpha = 0.3) +
  geom_line(color = "#2C3E50", size = 1) +
  scale_color_viridis(option = "C",name = "% missing",
                      labels = function(x) paste0(x, "%")) +
  annotate("text", x = max_idx, y = max_prop,
    label = paste0("Peak: ", round(max_prop, 3), "%"), vjust = -1.2,
    size = 6, fontface = "italic") +
  # breaks & limits
  scale_x_continuous(
    breaks = pretty(na_by_run$run, n = 10)) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Simulation run",
       y = "% Missing values (NA)",
       caption  = paste0("Overall avg: ", round(mean(na_by_run$propNA), 1), "%")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.caption = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank())


# FIGURE 13 --------------------------------------------------------------------
# Import results from baseline simulation
baseline_ta <- read_excel("tables/n_ta_output/n_ta_baseline.xlsx")
glimpse(baseline_ta)

# PLOT NDR DISTRIBUTION
ggplot(baseline_ta, aes(x = NDR)) +
  geom_histogram(binwidth = 0.0025, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 1 - baseline_ta$α[1], linetype   = "dashed", 
             color = "red", linewidth  = 1) +
  annotate("text",x  = 0.95,y = Inf,label = "Expected NDR",
           hjust = -0.1,vjust = 2, size = 6) +
  labs(x = "Neutral Detection Rate",y = "Frequency",
       caption  = paste0(
         "Average = ", round(baseline_ta$NDR[1], 3), 
         " | Runs ≥ 95% = ", baseline_ta$`>95% runs`[1], "%", 
         " | Number of runs = ", baseline_ta$Runs[1], 
         " | % NA = ", round(baseline_ta$`%NA`[1], 1), "%")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.caption  = element_text(size = 14))


# FIGURE 14 --------------------------------------------------------------------
# import sweep
n_snap_mu_df <- read_excel("tables/n_snap_output/n_snap_mu_params2.xlsx")
df <- n_snap_mu_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = `mu`)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.0125) +
  geom_line(aes(y = NDR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = expression(mu), breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )


# FIGURE 15 --------------------------------------------------------------------
n_snap_N_df <- read_excel("tables/n_snap_output/n_snap_N_params2.xlsx")
df <- n_snap_N_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )

# FIGURE 16 --------------------------------------------------------------------
n_snap_time_df <- read_excel("tables/n_snap_output/n_snap_time_params2.xlsx")
df <- n_snap_time_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = timesteps)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 150) +
  geom_line(aes(y = NDR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = "Time steps", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )


# FIGURE 17 --------------------------------------------------------------------
n_ta_mu_df <- read_excel("tables/n_ta_output/n_ta_mu_params2.xlsx")
df <- n_ta_mu_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = mu)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[1], size = 3) +
  scale_x_continuous(name = expression(mu), breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Purple = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )


# FIGURE 18 --------------------------------------------------------------------
n_ta_N_df <- read_excel("tables/n_ta_output/n_ta_N_params2.xlsx")
df <- n_ta_N_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[1], size = 3) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Purple = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )

# FIGURE 19 --------------------------------------------------------------------
n_ta_time_df <- read_excel("tables/n_ta_output/n_ta_time_params2.xlsx")
df <- n_ta_time_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = timesteps)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[1], size = 3) +
  scale_x_continuous(name = "Time steps", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Purple = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )


# FIGURE 20 --------------------------------------------------------------------
n_ta_tw_df <- read_excel("tables/n_ta_output/n_ta_tw_params2.xlsx")
df <- n_ta_tw_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = `Time window size`)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[1], size = 3) +
  scale_x_continuous(name = "w", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Purple = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )

# FIGURE 21 --------------------------------------------------------------------
mu_df <- read_excel("tables/cb_snap_output/cb_snap_overall.xlsx", sheet = "mu")

ggplot(mu_df, aes(x = `µ`)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.0125) +
  geom_line(aes(y = SSR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = expression(mu), breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "SSR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = SSR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )

# FIGURE 22 --------------------------------------------------------------------
N_df <- read_excel("tables/cb_snap_output/cb_snap_overall.xlsx", sheet = "N")

ggplot(N_df, aes(x = N)) +
  geom_col(aes(y = `%NA`),
           fill  = "grey80") +
  geom_line(aes(y = SSR), color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = SSR), color = plasma(5)[4], size = 3) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "SSR")
  ) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Orange = SSR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x       = element_text(size = 20),
    axis.text          = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)
  )

# FIGURE 23 --------------------------------------------------------------------
cb_df <- read_excel("tables/cb_snap_output/cb_snap_overall.xlsx", sheet = "b")

ggplot(cb_df, aes(x = b)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.025) +
  geom_line(aes(y = SSR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[4], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "b", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Orange = SSR"))


# FIGURE 24 --------------------------------------------------------------------
t_df <- read_excel("tables/cb_snap_output/cb_snap_overall.xlsx", sheet = "t")

ggplot(t_df, aes(x = timesteps)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 30) +
  geom_line(aes(y = SSR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[4], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Orange = SSR")),
       x = "Time Series Length (t)")

# FIGURE 25 --------------------------------------------------------------------
mu_df <- read_excel("tables/cb_ta_output/cb_ta_overall.xlsx", sheet = "mu")

ggplot(mu_df, aes(x = `µ`)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.012) +
  geom_line(aes(y = SSR),   color = plasma(5)[2], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[2], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = expression(mu), breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Blue = SSR")))

# FIGURE 26 --------------------------------------------------------------------
N_df <- read_excel("tables/cb_ta_output/cb_ta_overall.xlsx", sheet = "N")

ggplot(N_df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 35) +
  geom_line(aes(y = SSR),   color = plasma(5)[2], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[2], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Blue = SSR")))


# FIGURE 27 --------------------------------------------------------------------
cb_df <- read_excel("tables/cb_ta_output/cb_ta_overall.xlsx", sheet = "b")

ggplot(cb_df, aes(x = b)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.025) +
  geom_line(aes(y = SSR),   color = plasma(5)[2], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[2], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "b", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)) +
  labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Blue = SSR"))

# FIGURE 28 --------------------------------------------------------------------
t_df <- read_excel("tables/cb_ta_output/cb_ta_overall.xlsx", sheet = "t")

ggplot(t_df, aes(x = timesteps)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 30) +
  geom_line(aes(y = SSR),   color = plasma(5)[2], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[2], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Blue = SSR")),
    x = "Time Series Length (t)")


# FIGURE 29 --------------------------------------------------------------------
tw_df <- read_excel("tables/cb_ta_output/cb_ta_overall.xlsx", sheet = "tw")

ggplot(tw_df, aes(x = w)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 7.5) +
  geom_line(aes(y = SSR),   color = plasma(5)[2], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[2], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Blue = SSR")),
    x = "Time Window Width (w)")


# FIGURE 30 --------------------------------------------------------------------
c_df <- read_excel("tables/conf_snap_output/conf_snap_overall.xlsx", sheet = "c")

ggplot(c_df, aes(x = c)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.025) +
  geom_line(aes(y = SSR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[1], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "c", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Purple = SSR")),
    x = "c")

# FIGURE 31 --------------------------------------------------------------------
c_df <- read_excel("tables/conf_snap_output/conf_snap_overall.xlsx", sheet = "N")

ggplot(c_df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 25) +
  geom_line(aes(y = SSR),   color = plasma(5)[1], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[1], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Purple = SSR")),
    x = "N")

# FIGURE 32 --------------------------------------------------------------------
c_df <- read_excel("tables/conf_ta_output/conf_ta_overall.xlsx", sheet = "c")

ggplot(c_df, aes(x = c)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.012) +
  geom_line(aes(y = SSR),   color = plasma(5)[3], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[3], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "c", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Magenta = SSR")),
    x = "c")

# FIGURE 33 --------------------------------------------------------------------
c_df <- read_excel("tables/conf_ta_output/conf_ta_overall.xlsx", sheet = "N")

ggplot(c_df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 30) +
  geom_line(aes(y = SSR),   color = plasma(5)[3], linewidth = 1.2) +
  geom_point(aes(y = SSR),  color = plasma(5)[3], size = 3) +
  scale_y_continuous(
    name = "NA",
    sec.axis = sec_axis(
      ~ ., name   = "SSR")) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)) +
  labs(caption = labs(caption = paste0(
    "Runs = ", 1000, "    |    ",
    "Grey = NA, Magenta = SSR")),
    x = "N")

