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
ggbiplot(pca, obs.scale = 1, var.scale = 1)
loadings <- pca$rotation
print(loadings)