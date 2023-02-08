
#==============================================================================#
# Libraries and functions ----
#==============================================================================#

library(data.table)
library(ggplot2)
library(zoo)

source("functions.R")

#==============================================================================#
# Hardcodes ----
#==============================================================================#

# where we're pulling data and storing plots 
data_path <- "../output/"
ests_path <- paste0(data_path, "impute_ests/")
plot_path <- paste0(data_path, "plots/")

# input files
data_file <- paste0(data_path, "/bcsignals/bcsignals_none.csv")

# dates to get beta distributions for
dates <- as.yearmon(c("Jun 1990", "Jun 2000", "Jun 2010"))

#==============================================================================#
# Read in data ----
#==============================================================================#

# Read in the covariance matrices
covs <- lapply(dates, function(dd) {

    # Match file format
    dd_file <- gsub(" ", "", dd)

    # Return the covariance matrix to list 
    as.matrix(read.csv(paste0(ests_path, "estR_", dd_file, ".csv"), 
        row.names = 1))
})
names(covs) <- dates

# Read in the observation dataset to find missing obs
data <- fread(data_file)[as.yearmon(yyyymm) %in% dates][order(permno)]
data[, yyyymm := as.yearmon(yyyymm)]

#==============================================================================#
# Extract betas ----
#==============================================================================#

# for each month, get each observation's coefficients
beta_lists <- list()
for (dd in 1:length(dates)) {
    this_cov <- covs[[as.character(dates[dd])]]
    vars <- colnames(this_cov) # Ensure variables line up correctly
    beta_lists[[dd]] <- getBetas(this_cov[vars, vars], 
        is.na(data[yyyymm == dates[dd], ..vars])) # is.na(data) returns mat of bools
}
names(beta_lists) <- dates

# basically just reducing the betas down to one huge vector per month 
# to plot a histogram
beta_lists <- lapply(beta_lists, function(x) lapply(x, as.vector))
all_the_betas <- lapply(beta_lists, function(x) Reduce(c, x))
all_the_betas <- lapply(all_the_betas, as.data.table)
all_the_betas <- lapply(names(all_the_betas), function(x) all_the_betas[[x]][,month := x])

# finally, combine all the betas into a plotable data.frame
beta_dt <- Reduce(rbind, all_the_betas)
beta_dt[, beta_winsor := winsorize(V1)]

rm(beta_lists, all_the_betas)

#==============================================================================#
# Plot ----
#==============================================================================#

# make the density plot of the beta values
beta_dist <- ggplot(data = beta_dt, aes(x = V1, group = month)) + 
  geom_density(aes(colour = month, linetype = month), size = 1) +
  labs(
    x = "Imputation Slope",
    y = "Density",
    color = "Month",
    linetype = "Month"
  ) + 
  scale_x_continuous(
    limits = c(-1,1),
    breaks = seq(-1,1,0.25)
  ) + 
  scale_y_continuous(
    expand = expand_scale(mult = c(0,0.05), add = 0),
    breaks = seq(0,10,2)
  ) +
  theme_bw() + 
  theme(
    legend.position = c(0.8, 0.8),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    panel.grid = element_line(colour = "grey", size = 0.2),
    axis.text = element_text(size = 13),
    text = element_text(size = 15)
  ) 

# beta_dist


# Output the plots
ggsave(plot = beta_dist
       , filename = paste0(plot_path, "beta_dist.pdf")
       , height = 4, width = 8, scale = 0.9
     )



# Numbers for the text ---------------------------------------------------

mean(abs(beta_dt$V1))

quantile(beta_dt$V1, c(2.5,5,50,95,97.5)/100)

hist(beta_dt$beta_winsor)
