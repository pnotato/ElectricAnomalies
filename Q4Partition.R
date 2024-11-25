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


# ----- Partition ----- #

# Formatting into POSIXct
library(lubridate)
library(depmixS4) # For HMM modeling

# Assuming df is already loaded
df$DateTime <- as.POSIXct(paste(df$Date, df$Time), format = "%d/%m/%Y %H:%M:%S")
df$Date <- as.Date(df$DateTime) # Ensure 'Date' is in Date format

# Filter data
df_filtered <- subset(df, year(Date) == 2010 & wday(Date) == 1)
df_filtered <- subset(df_filtered, hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) >= 12 & hour(as.POSIXct(df_filtered$Time, format = "%H:%M:%S", tz = "UTC")) < 15)

df_filtered$week <- week(df_filtered$Date)


# There are 329 / 7 = 47 weeks in our 2010 data set. To divide into 10 partitions, we will make each partition 5 weeks long, 
# with the last partition being 2 weeks long.

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


