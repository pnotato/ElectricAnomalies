df <- read.table("household_power_consumption.txt", header = TRUE, sep = ";")

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

# ----- HMM Training and Testing ----- #

library(lubridate)
library(depmixS4) # For HMM modeling

# Assuming df is already loaded
df$DateTime <- as.POSIXct(paste(df$Date, df$Time), format = "%d/%m/%Y %H:%M:%S")
df$Date <- as.Date(df$DateTime) # Ensure 'Date' is in Date format

# Filter data
df_filtered <- subset(df, year(Date) >= 2006 & year(Date) <= 2009 & wday(Date) == 1)
df_filtered <- subset(df_filtered, hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) >= 12 & hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) < 15)

# Preprocess: Assuming 'Global_active_power' is numeric. If it's factor or character due to missing values or placeholders like "?", convert it and handle NAs.
df_filtered$Global_active_power <- as.numeric(as.character(df_filtered$Global_active_power))
df_filtered$Sub_metering_1 <- as.numeric(as.character(df_filtered$Sub_metering_1))
df_filtered <- na.omit(df_filtered) # Removing rows with NAs after conversion


# Setup for HMM
logLikelihoods <- numeric(17)
BICs <- numeric(17)
states <- 4:20


for (nstates in states) {
  model <- depmix(list(Global_active_power ~ 1, Sub_metering_1 ~ 1), data = df_filtered, 
                  nstates = nstates, family = list(gaussian(), multinomial())) 
  fitModel <- fit(model) 
  print(summary(fitModel))

  logLikelihoods[nstates - 3] <- logLik(fitModel)
  BICs[nstates - 3] <- BIC(fitModel)
}

# Plotting results
plot(states, logLikelihoods, type = "l", main = "Log-Likelihoods vs. Number of States",
     xlab = "Number of States", ylab = "Log-Likelihood", col = "blue")
plot(states, BICs, type = "l", main = "BIC vs. Number of States",
     xlab = "Number of States", ylab = "BIC", col = "red")