
library(raster)
library(scales)

# function to create values with known correlation
complement <- function(y, rho, x) {
  if (missing(x)) x <- (rnorm(length(y), mean=2, sd=10))^2 # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

# function to split values into bins that have the same amount of elements 
# (grid cells)
EqualFreq <-function(x,n,include.lowest=TRUE,...){
  nx <- length(x)    
  id <- round(c(1,(1:(n-1))*(nx/n),nx))
  
  breaks <- sort(x)[id]
  if (sum(duplicated(breaks)) > 0 ) {
    stop("n is too large.")
  }  
  cut(x,breaks,include.lowest=include.lowest,...)
}




#### create a made-up tenure layer using the continuous cost layers as templates ####
cost_mask <- raster("Baseline_costs.tif")
pri <- raster("Baseline_sp_rich.tif")
pri_v <- na.omit(getValues(pri))
cost_files <- c("Baseline_costs.tif", "LargeCV_Neg.tif", "SmallCV_Neg.tif",  
                "LargeCV_Pos.tif", "SmallCV_Pos.tif")
cost_names <- paste0("tenure_cost_", c("baseline", "neg_high", "neg_low", 
                                       "pos_high", "pos_low"))

tenure_summary <- data.frame(name = cost_names, rho = NA, CV = NA, min = NA, 
                             max = NA)
cost_bins <- list()
tenure_cost_values <- matrix(NA, 7,6, 
                             dimnames = list(NULL, 
                                             c("tenure", "baseline", 
                                               "neghigh", "neglow", 
                                               "poshigh", "poslow")))
tenure_cost_values[,1] <- 1:7


for (i in seq_along(cost_files)){
  r <- raster(cost_files[i])
  
  # unify grid cells
  r <- mask(r, cost_mask)
  r <- mask(r, pri)
  
  # find the min and max values of those bins that split the data into
  # equal area tenure categories
  r_values <- na.omit(getValues(r))
  y <- EqualFreq(r_values, 7)

  bins <- unlist( strsplit( as.character( levels(y) ), ',') )
  bins <- gsub('\\[', '', gsub('\\]', '', gsub('\\(', '', bins)))
  
  # create a reclassification table
  rcl_table <- matrix(as.numeric(bins), 7, 2, byrow = TRUE)
  rcl_table[1,1] <- 0
  rcl_table[7,2] <- max(r_values)
  rcl_table <- cbind(rcl_table, 1:7)
  
  cost_bins <- append(cost_bins, list(rcl_table))
  
  # create the dummy tenure layer
  tenure <- reclassify(r, rcl_table)
  
  # gather mean values and create tenure averaged cost layer
  cell.data <- data.frame(tenure = na.omit(getValues(tenure)),
                          cost = r_values)
  cell.data <- aggregate(.~tenure, data = cell.data, mean)
  tenure_cost_values[,i+1] <- cell.data$cost
  tenure_cost <- reclassify(tenure, rcl = as.matrix(cell.data[c('tenure','cost')]))
  
  # collect data
  new_cost_v <- na.omit(getValues(tenure_cost))
  
  tenure_summary$rho[i] <- cor(pri_v, new_cost_v, method = "spearman")
  tenure_summary$CV[i] <- sd(new_cost_v) / mean(new_cost_v)
  tenure_summary$min[i] <- min(new_cost_v)
  tenure_summary$max[i] <- max(new_cost_v)
  
  writeRaster(tenure, 
              paste0(gsub("cost_", "", cost_names[i]), "_binned.tif"), 
              format = 'GTiff', overwrite = TRUE) 
}

names(cost_bins) <- cost_names
tenure_cost_values <- as.data.frame(tenure_cost_values)

# print summaries
tenure_summary
tenure_cost_values





#### change CV ####
mod_tenure_cost_values <- tenure_cost_values
mod_tenure_summary <- data.frame(name = cost_names, rho = NA, CV = NA, min = NA, 
                                 max = NA)
mod_tenure_summary[1,] <- tenure_summary[1,]
tenure_files <- grep("_binned.tif", list.files(), value = TRUE)[-1]

# set target CV
t_cv_low <- 0.3

#-------------------------------------------------------------------------------
# negative correlation

tenure_neghigh <- raster(tenure_files[1])

## high CV
mod_tenure_cost_values$neghigh <- 10 ^ (rescale(tenure_cost_values$neghigh, 
                                                to = c(0, 6)))
tenure_cost_neghigh <- reclassify(tenure_neghigh, 
                                  rcl = as.matrix(mod_tenure_cost_values[
                                    c("tenure", "neghigh")]))

# checks
mod_tenure_summary$rho[2] <-
  cor(na.omit(getValues(pri)), na.omit(getValues(tenure_cost_neghigh)), method = "spearman")

new_cost_v <- na.omit(getValues(tenure_cost_neghigh))
mod_tenure_summary$CV[2] <- sd(new_cost_v) / mean(new_cost_v)
mod_tenure_summary$min[2] <- min(new_cost_v)
mod_tenure_summary$max[2] <- max(new_cost_v)

writeRaster(tenure_cost_neghigh, 
            paste0(cost_names[2], "_v3.tif"), 
            format = 'GTiff', overwrite = TRUE)


## low CV
tenure_neglow <- raster(tenure_files[2])
t_mean <- mean(tenure_cost_values$neglow)
t_sd <- t_mean * t_cv_low

tmp <- t_sd * ((tenure_cost_values$neglow - mean(tenure_cost_values$neglow)) / sd(tenure_cost_values$neglow)) + t_mean
mod_tenure_cost_values$neglow <- tmp

tenure_cost_neglow <- reclassify(tenure_neglow, 
                                 rcl = as.matrix(mod_tenure_cost_values[
                                   c("tenure", "neglow")]))

# checks
mod_tenure_summary$rho[3] <- 
  cor(na.omit(getValues(pri)), na.omit(getValues(tenure_cost_neglow)), method = "spearman")

new_cost_v <- na.omit(getValues(tenure_cost_neglow))
mod_tenure_summary$CV[3] <- sd(new_cost_v) / mean(new_cost_v)
mod_tenure_summary$min[3] <- min(new_cost_v)
mod_tenure_summary$max[3] <- max(new_cost_v)

writeRaster(tenure_cost_neglow, 
            paste0(cost_names[3], "_v3.tif"), 
            format = 'GTiff', overwrite = TRUE)


#-------------------------------------------------------------------------------
# Positive correlation

tenure_poshigh <- raster(tenure_files[3])

## high CV
mod_tenure_cost_values$poshigh <- 10 ^ (rescale(tenure_cost_values$poshigh, 
                                                to = c(0, 6)))
tenure_cost_poshigh <- reclassify(tenure_poshigh, 
                                  rcl = as.matrix(mod_tenure_cost_values[
                                    c("tenure", "poshigh")]))

# checks
mod_tenure_summary$rho[4] <-
  cor(na.omit(getValues(pri)), na.omit(getValues(tenure_cost_poshigh)), method = "spearman")

new_cost_v <- na.omit(getValues(tenure_cost_poshigh))
mod_tenure_summary$CV[4] <- sd(new_cost_v) / mean(new_cost_v)
mod_tenure_summary$min[4] <- min(new_cost_v)
mod_tenure_summary$max[4] <- max(new_cost_v)

writeRaster(tenure_cost_poshigh, 
            paste0(cost_names[4], "_v3.tif"), 
            format = 'GTiff', overwrite = TRUE)

## low CV
tenure_poslow <- raster(tenure_files[4])
t_mean <- mean(tenure_cost_values$poslow)
t_sd <- t_mean * t_cv_low

tmp <- t_sd * ((tenure_cost_values$poslow - mean(tenure_cost_values$poslow)) / sd(tenure_cost_values$poslow)) + t_mean
mod_tenure_cost_values$poslow <- tmp

tenure_cost_poslow <- reclassify(tenure_poslow, 
                                 rcl = as.matrix(mod_tenure_cost_values[
                                   c("tenure", "poslow")]))

# checks
mod_tenure_summary$rho[5] <-
  cor(na.omit(getValues(pri)), na.omit(getValues(tenure_cost_poslow)), method = "spearman")

new_cost_v <- na.omit(getValues(tenure_cost_poslow))
mod_tenure_summary$CV[5] <- sd(new_cost_v) / mean(new_cost_v)
mod_tenure_summary$min[5] <- min(new_cost_v)
mod_tenure_summary$max[5] <- max(new_cost_v)

writeRaster(tenure_cost_poslow, 
            paste0(cost_names[5], "_v3.tif"),
            format = 'GTiff', overwrite = TRUE)


# print summaries
mod_tenure_summary
mod_tenure_cost_values