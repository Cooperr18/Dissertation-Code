#########################################################
################ CONTENT BIAS SNAPSHOT ##################
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

