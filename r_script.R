# Retrieve data from text file
df <- read.table("household_power_consumption.txt", header = TRUE, sep = ";")
# Checked to see if dataframe has any missing values using print(any(is.na(df))) -> returns true

###########################
# --- Linear Interpolation Function ---
# Returns the resulting linear interpolation on a given feature
linearInterpolation <- function(feature) {
  x <- seq_along(df[[feature]])
  y <- df[[feature]]
  result <- NULL
  
  if(any(is.na(y))) { # Check if there is any NA values in the column
    result <- approx(x, y, method = "linear", n = length(x))
    return(result$y) # Return the resulting linear interpolation of the feature
  }
  return(y) # return the original value if no missing values
}
###########################
# --- Z Score ---
# Computes the z score
computeZScore <- function(feature) {
  if(is.numeric(df[[feature]])) {
    x <- df[[feature]]
    u <- mean(x, trim = 0, na.rm = FALSE)
    sd <- sd(x, na.rm = FALSE)
    z_score <- (x - u) / sd
    return(z_score)
  }
  return(NULL)
}
###########################
# ---- Script Start ----- #
exclude_features <- c("Date", "Time")
for(feature in names(df)) {
  if(is.character(df[[feature]]) && !(feature %in% exclude_features)) {
    df[[feature]] <- as.numeric(df[[feature]])
  }
  df[[feature]] <- linearInterpolation(feature)
}
z_scores <- df
for(feature in names(df)) {
  z_scores[[feature]] <- computeZScore(feature)
}

# ----- PCA ----- #
library(ggbiplot)

pca <- prcomp(z_scores, scale. = FALSE)
summary(pca)
#ggbiplot(pca, obs.scale = 1, var.scale = 1)
#loadings <- pca$rotation
#print(loadings)

# ----- HMM Training and Testing ----- #

# Formatting into POSIXct
library(lubridate)
library(depmixS4) # For HMM modeling

# Assuming df is already loaded
df$DateTime <- as.POSIXct(paste(df$Date, df$Time), format = "%d/%m/%Y %H:%M:%S")
df$Date <- as.Date(df$DateTime) # Ensure 'Date' is in Date format

# Filter data
df_filtered <- subset(df, year(Date) >= 2006 & year(Date) <= 2009 & wday(Date) == 1)
df_filtered <- subset(df_filtered, hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) >= 16 & hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) < 19)

# Preprocess: Assuming 'Global_active_power' is numeric. If it's factor or character due to missing values or placeholders like "?", convert it and handle NAs.
df_filtered$Global_active_power <- as.numeric(as.character(df_filtered$Global_active_power))
df_filtered <- na.omit(df_filtered) # Removing rows with NAs after conversion

# Setup for HMM
set.seed(10000)
states <- 1:16
logLikelihoods <- numeric(length(states))
BICs <- numeric(length(states))


# Assuming 'ntimes' needs to reflect the actual sequence of observations per year or another logical segmentation
# For simplicity, let's assume an equal distribution across the segments. Adjust this according to your specific needs.
# The number of states is varied from 4 to 20 as per your instructions.

# Train data HMM
for (state in 4:13) {
  model <- depmix(list(Global_active_power ~ 1, Global_intensity ~ 1), data = df_filtered, 
                  nstates = state, ntimes = nrow(df_filtered), family = list(gaussian(), gaussian())) # Adjusted to Gaussian assuming continuous data
  fitModel <- fit(model) # Fit the model
  print("------- Fit Model -------")
  print(summary(fitModel))
  # Store log-likelihood and BIC for each model
  logLikelihoods[state] <- logLik(fitModel)
  BICs[state] <- BIC(fitModel)
} 

# Test Data HMM
test_data <- subset(df, year(Date) == 2010 & wday(Date) == 1)
best_model <- depmix(list(Global_active_power ~ 1, Global_intensity ~ 1), data = df_filtered, 
                              nstates = 13, ntimes = nrow(df_filtered), family = list(gaussian(), gaussian()))
fit_best_model <- fit(best_model)
test_logLikelihood <- forwardbackward(fit_best_model, newdata = test_data)
print(test_logLikelihood)
print(logLikelihoods)
print(BICs)
print(states)
# Plotting results
plot(states, logLikelihoods, type = "l", main = "Log-Likelihoods and BICs vs. Number of States",
     xlab = "Number of States", ylab = "Log-Likelihood / BIC", col = "blue", ylim = range(c(logLikelihoods, BICs)))

# Add BICs to the plot
lines(states, BICs, type = "l", col = "red")

# Add legend
legend("topright", legend = c("Log-Likelihood", "BIC"),
       col = c("blue", "red"), lty = 1, cex = 0.8)

# Normalize log-likelihoods
train_normalized_likelihoods <- (logLikelihoods - min(logLikelihoods)) / (max(logLikelihoods) - min(logLikelihoods))
test_normalized_likelihoods <- (test_logLikelihood$logLike - min(logLikelihoods)) / (max(logLikelihoods) - min(logLikelihoods))

# Filter data
df_filtered <- subset(df, year(Date) == 2010 & wday(Date) == 1)
df_filtered <- subset(df_filtered, hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) >= 12 & hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) < 15)
df_filtered$week <- week(df_filtered$Date)

# Divide test data into 10 equal parts
partitionlist <- list()
for (partition in 1:10) {
  start_week <- (partition - 1) * 5 + 1
  end_week <- partition * 5
  
  if (partition == 10) {
    end_week <- start_week + 1 
  }
  currentpartition <- df_filtered[df_filtered$week >= start_week & df_filtered$week <= end_week, ]
  partitionlist[[partition]] <- currentpartition
}

# Compute the log likelihood using the best fit model of each subset
partition_likelihoods <- vector("list", length = length(partitionlist))
for (i in 1:length(partitionlist)) {
  fit_best_model <- fit(best_model)
  partition_likelihood <- forwardbackward(fit_best_model, newdata = partitionlist[[i]])
  partition_likelihoods[[i]] <- partition_likelihood
}
print(partition_likelihoods)  # Print the entire list after computation


