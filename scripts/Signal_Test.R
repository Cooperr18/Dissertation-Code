# TSINFER ---------------------------

library(here)
library(signatselect)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)


# Data setup

library(evoarchdata)
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
# gather reorganizes the data to a long format with its columns, and mutate and fct_relevel indicates the order in which you want them to be showed when you plot them. 
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

n <- 5
df_with_rolling_idx <- function(df, window = n) {
  nr <- nrow(df)
  w <- window       # window size
  i <- 1:nr         # indices of the rows
  iw <-
    embed(i, w)[, w:1]   # matrix of rolling-window indices of length w
  wnum <- rep(1:nrow(iw), each = w)   # window number
  inds <-
    i[c(t(iw))]         # the indices flattened, to use below
  zw <- sapply(df, '[', inds)
  zw <- transform(data.frame(zw), w = wnum)
  return(zw)
}

merzbach_long_sig_mid_time_point <- 
  list_of_dfs %>% 
  bind_rows(.id = "type") %>% # to get rolling window of n
  nest(-type) %>%
  mutate(rolled = map(data, df_with_rolling_idx)) %>% 
  unnest(rolled) %>% 
  mutate(unid = str_glue('{type}_{w}')) %>% 
  nest(-unid) %>% 
  mutate(fit_test = map(data,
                        ~fit_safely(time = .x$time,
                                    v =    .x$frequency))) %>%
  mutate(fit_p = map(fit_test, ~.x$result %>% bind_rows)) %>%
  unnest(fit_p) %>%
  mutate(sig = ifelse(fit_p <= 0.05, "selection", "neutral")) %>% 
  unnest(data) 

# make type a factor so we can order the plots nicely 
merzbach_long_sig_mid_time_point$type <- 
  fct_relevel(merzbach_long_sig_mid_time_point$type,
              decoration_types[decoration_types %in% merzbach_long_sig_mid_time_point$type])

# plot with overall time-series results also
# harmonize some variable names first

ceramics_lbk_merzbach_long_sig_to_plot_with_others <- 
  ceramics_lbk_merzbach_long_sig %>% 
  rename( time = Phase_num,
          count_this_one = value,
          type = variable) %>% 
  filter(count_this_one != 0) %>% 
  arrange(type, time) %>% 
  mutate(type = fct_relevel(type, decoration_types))

sig_decorations <- 
  ceramics_lbk_merzbach_long_sig_to_plot_with_others %>% 
  filter(sig == "selection") %>% 
  pull(type) %>% 
  as.character() %>% 
  unique()

ggplot()  +
  geom_line(data = merzbach_long_sig_mid_time_point %>% 
              filter(sig == "selection"),
            aes(time,
                count_this_one,
                group = type),
            size = 5,
            colour = "grey80",
            lineend = "round") +
  geom_point(data = merzbach_long_sig_mid_time_point %>% 
               filter(sig == "selection"),
             aes(time,
                 count_this_one,
                 group = type),
             size = 5,
             colour = "grey80") +
  geom_line(data = ceramics_lbk_merzbach_long_sig_to_plot_with_others,
            aes(time,
                count_this_one, 
                group = type,
                colour = sig)) +
  geom_point(data = ceramics_lbk_merzbach_long_sig_to_plot_with_others,
             aes(time,
                 count_this_one, 
                 group = type,
                 colour = sig,
                 shape = sig))  +
  scale_color_viridis_d(name = "", 
                        begin = 0.25, 
                        end = 0.75) + 
  guides(shape = FALSE) +
  facet_wrap( ~ type, scales = "free_y") +
  theme_minimal(base_size = 8) +
  ggtitle(str_glue('Application of the FIT to decoration frequency data from Merzbach.\nShading highlights the data points where FIT identifies selection'))+
  labs(y="Number of counts",x="Time series")

decoration_types_tsinfer_output <- 
  list_of_dfs_three_or_more %>% 
  map(~tsinfer(
    tvec = .x$time,
    bvec = .x$count_this_one,
    nvec = .x$count_others,
    verbose = FALSE
  )) %>% 
  bind_rows(.id = "type")

# inspect the output
decoration_types_tsinfer_output

decoration_types_tsinfer_output <- 
  list_of_dfs_three_or_more %>% 
  purrr::map_dfr(
    ~ tsinfer(
      tvec    = .x$time,
      bvec    = .x$count_this_one,
      nvec    = .x$count_others,
      verbose = FALSE
    ),
    .id = "type"
  ) %>% 
  left_join(
    df_fit_test_results %>% 
      select(type, fit_p, sig),
    by = "type"
  ) %>% 
  # Rename
  rename(
    p_value       = fit_p,
    inference_tag = sig
  )

# Print result
print(decoration_types_tsinfer_output)

# Add p-value and inference
decoration_types_tsinfer_selected <- 
  decoration_types_tsinfer_output %>%
  select(
    type, 
    s, 
    f0, 
    LL, 
    p_value,        
    inference_tag
  )

# Print table
print(decoration_types_tsinfer_selected)

# Install packages

install.packages(c("officer", "flextable"))

library(officer)
library(flextable)

# flextable from data
ft <- flextable(decoration_types_tsinfer_selected)
ft <- theme_vanilla(ft) %>%
  autofit()

# Table to Word
doc <- read_docx() %>%
  body_add_par("Results of FIT Analysis for Decoration Types", style = "heading 1") %>%
  body_add_flextable(ft) 

# save in directoru
print(doc, target = "decoration_types_tsinfer_results.docx")


# Visualize by Relative frequencies
library(viridis)
library(stringr)

# Convert frequencies
df_rel <- ceramics_lbk_merzbach_long_sig_to_plot_with_others %>%
  group_by(time) %>% 
  mutate(
    total_counts = sum(count_this_one, na.rm = TRUE),
    frequency    = count_this_one / total_counts
  ) %>% 
  ungroup()

ggplot(df_rel, aes(
  x      = time,
  y      = frequency,
  group  = type,
  shape  = type,    # variant different shapes
  colour = sig      # inference by colour
)) +
  geom_line(aes(colour = sig), size = 1) +  # variants by colour
  geom_point(size = 3) +                    # points by variant
  scale_shape_manual(
    name   = "Decoration\nVariant",
    values = seq_along(unique(df_rel$type))  # different shapes
  ) +
  scale_colour_viridis_d(
    name  = "Inference",
    begin = 0.25,
    end   = 0.75
  ) +
  guides(
    shape  = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 8) +
  labs(
    title = str_glue(
      "Relative frequencies by variant over time\n",
      "Application of the FIT to Merzbach decoration data"
    ),
    x = "Time series",
    y = "Relative frequency"
  )

