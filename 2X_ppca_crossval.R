# 2023 09 andrew's check on ppca cross val
# doc
# https://www.bioconductor.org/packages/release/bioc/vignettes/pcaMethods/inst/doc/missingValues.pdf

#===============================================================================
# Packages 
#===============================================================================

library(data.table)
library(optparse)
library(pcaMethods)
library(tidyverse)
library(zoo)

source("functions.R")

#===============================================================================
# Set up
#===============================================================================

option_list <- list(
    optparse::make_option(c("--yearm_select"),
        type = "character",
        default = "Jun 1985",
        help = "yearmon to select in format `Mmm YYYY`"),
    optparse::make_option(c("--max_numPCs"),
        type = "numeric",
        default = 80,
        help = "maximum number of PCs to try in imputation crossval"),
    optparse::make_option(c("--num_numPCs"),
        type = "numeric",
        default = 20,
        help = "number of different PCs to try, ranging from 1 to `--max_numPCs`"),
    optparse::make_option(c("--num_folds"),
        type = "numeric",
        default = 5,
        help = "number of folds in crossval"),
    optparse::make_option(c("--maxiter"),
        type = "numeric",
        default = 500,
        help = "max number of iterations in pPCA")
)

# Unpack the options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# File paths for the project
getFilePaths()

#===============================================================================
# Hardcoded settings
#===============================================================================

seed <- 864

#===============================================================================
# Read in data
#===============================================================================

# Read in each dataset we need
# immediately filter here for mem
signals <- fread(paste0(FILEPATHS$data_path, "bcsignals/bcsignals_none.csv"))[
    as.yearmon(yyyymm) == as.yearmon(opt$yearm_select)
]
crsp <- fread(paste0(FILEPATHS$data_path, "raw/crsp_data.csv"))[
    as.yearmon(yyyymm) == as.yearmon(opt$yearm_select)
]
me_breaks <- fread(paste0(FILEPATHS$data_path, "raw/ME_Breakpoints_20_50_clean.csv"))[
    as.yearmon(yyyymm) == as.yearmon(opt$yearm_select)
]

# Pretty dates
signals[, yyyymm := as.yearmon(yyyymm)] 
crsp[, yyyymm := as.yearmon(yyyymm)]
me_breaks[, yyyymm := as.yearmon(yyyymm)]

# Get the good signals
signals_good <- names(signals %>% select(-c("permno", "yyyymm")))

# add me to signals. select necessary columns
signals <- merge(
    signals, 
    crsp[, .(permno, yyyymm, me)],
    by = c("permno", "yyyymm"),
    all.x = TRUE
)[ # Duplicative for now, but in case we put better check on signals_good
    , .SD,
    .SDcols = c(signals_good, "permno", "yyyymm", "me")
]

#===============================================================================
# prep data 
#===============================================================================

# Group by me breakpoints in corresponding month.
# Divide me by 1000 because it's in thousands, but breakpoint data in millions
# Using Fama-French (2008) group definitions
signals[
    me_breaks,
    group := dplyr::case_when(
        me/1000 > p50 ~ "Big",
        p20 < me/1000 & me/1000 < p50 ~ "Small",
        TRUE ~ "Micro"
    ),
    on = .(yyyymm)
]

# make signals matrix
xsignal <- as.matrix(signals[, .SD, .SDcols = c("permno", signals_good)],
    rownames = "permno") # permno as rowname not col. see `as.matrix.data.table`

# shuffle observed stock-signal indexes
iobs <- which(!is.na(xsignal))
set.seed(seed)
iobs <- sample(iobs, length(iobs))

# split into folds
ifold <- split(iobs, ceiling(seq_along(iobs) / length(iobs) * opt$num_folds))

#===============================================================================
# declare function
#===============================================================================

# function for finding prediction error for a fold
fold_err <- function(itest, numPCs){
    
    itrain <- setdiff(iobs, itest)
    
    # force test ids to missing
    train_data <- xsignal
    train_data[itest] <- NA
    test_data <- xsignal
    test_data[itrain] <- NA
    
    # remove stocks with no obs
    badstock <- apply(is.na(train_data), 1, all)
    train_data <- train_data[!badstock, ]
    test_data <- test_data[!badstock, ]
    
    # remove signals with no obs
    badsignal <- apply(is.na(train_data), 2, all)
    train_data <- train_data[, !badsignal]
    test_data <- test_data[, !badsignal]
    
    # estimate
    pc <- pca(train_data,
        method = "ppca",
        nPcs = numPCs,
        seed = 0, # Consistency
        maxIterations = opt$maxiter,
        scale = "none",
        center = FALSE
    )
    
    # measure rmse (only measured on obs data, hence why force NAs in 2nd line)
    test_pred <- scores(pc) %*% t(loadings(pc))
    test_pred[which(is.na(test_data))] <- NA

    errsum <- melt(
        as.data.table(test_data - test_pred, keep.rownames = "permno")[
            signals[, .(permno = as.character(permno), group)],
            on = .(permno)
        ],
        id.vars = c("permno", "group")
    )[
        !is.na(value),
        .(
            rmse = sqrt(mean(value^2, na.rm = TRUE)),
            nobs = sum(!is.na(value)),
            numPCs = numPCs,
            badstock_all = sum(badstock)
        ),
        by = .(group)
    ]
    
    return(errsum)
    
} # end function fold_err

#===============================================================================
# rmse vs num_PCs
#===============================================================================

# loop over numPCs
numPCs_grid <- unique(round(seq(1, opt$max_numPCs, length.out = opt$num_numPCs)))
errsum <- tibble()
for (i in 1:opt$num_numPCs) {
    cat("numPCs =", numPCs_grid[i], "\n")
    numPCs <- numPCs_grid[i]

    # evaluate prediction error for each fold
    for (j in 1:opt$num_folds) {
        cat("\tfold =", j, "\n")
        temperr <- fold_err(ifold[[j]], numPCs)[, fold := j]
        errsum <- errsum %>% rbind(temperr)
    } # for j
    
} # for i in 1:num_numPCs

stop()

# Aggregate over folds and signals
rmsedat <- errsum[,
    .(
        mean = mean(rmse, na.rm = TRUE),
        sd = sd(rmse, na.rm = TRUE),
        SE = sd(rmse, na.rm = TRUE) / sqrt(opt$num_folds),
        nobs = sum(nobs, na.rm = TRUE)
    ),
    by = .(numPCs, group)
][, # Better for plot
    group := factor(group, levels = c("Micro", "Small", "Big")) 
]

#===============================================================================
# Plots
#===============================================================================

# plot error by group
cval_plot <- ggplot(rmsedat, aes(x = numPCs, y = mean, group = group)) +
    geom_line(aes(color = group, linetype = group), size = 0.5) +
    geom_point(aes(color = group, shape = group)) +
    geom_errorbar(aes(ymin = mean - 2*SE, ymax = mean + 2*SE, color = group)) +
    theme_bw() +
    theme(
        legend.position = c(77, 75) / 100,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing.y = unit(0.005, 'cm'),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.box = 'horizontal',
        legend.key.width = unit(0.55, 'cm'),
        legend.key.height = unit(0.38, 'cm'),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)
    ) + 
    labs(
        colour = "Size Group",
        shape = "Size Group",
        linetype = "Size Group",
        #title = paste0(
        #    "ppca oos error ",
        #    as.yearmon(opt$yearm_select),
        #    " num_folds = ",
        #    opt$num_folds
        #),
        x = "Number of PCs",
        y = "RMSE",
        caption = paste0('error bars are 2 SE')
    )

dir.create(paste0(FILEPATHS$out_path, "plots/"))

# Output
ggsave(
    plot = cval_plot,
    filename = paste0(
        FILEPATHS$out_path,
        "plots/",
        "ppca_cval_",
        gsub(" ", "", as.character(opt$yearm_select)), 
        "_", opt$num_folds, "folds",
        ".pdf"
    ),
    width = 8, 
    height = 7,
    unit = "in",
    scale = 0.55
)

# size counts for reference
signals[,
    .(nstock = .N, pct = .N/nrow(signals)),
    by = .(group)
]

