#install.packages("corrr")

library(raster)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(GGally)

# Define the directory containing the TIFF files
directory_path <- "/path/to/wc2.1_30s_bio"

# List all TIFF files in the directory
tiff_files <- list.files(directory_path, pattern = "\\.tif$", full.names = TRUE)

# Load the raster files into a stack
bioclim_raster <- stack(tiff_files)

##get LAT LON from diarrhea preprocessed data
# Create cluster_id and extract distinct LAT LON for each cluster
distinct_under5child_diarrhea <- under5children_filter3 %>%
  distinct(V000, V001, LATNUM, LONGNUM) %>%  # Extract distinct latitude and longitude pairs
  filter(LATNUM != 0.00 & LONGNUM != 0.00) %>%  # Remove rows where LATNUM or LONGNUM is 0.00
  mutate(cluster_id = row_number())  # Add unique cluster_id



# Function to generate 8 additional points around a central point
generate_additional_points <- function(lon, lat, distance_km, cluster_id) {
  directions <- data.frame(
    dlon = c(0, 0, distance_km, -distance_km, distance_km, -distance_km, distance_km, -distance_km),
    dlat = c(distance_km, -distance_km, 0, 0, distance_km, distance_km, -distance_km, -distance_km)
  )
  directions <- directions * 0.009  # Convert km to degrees approximately
  additional_points <- data.frame(
    LONGNUM = lon + directions$dlon,
    LATNUM = lat + directions$dlat,
    cluster_id = cluster_id,
    is_original = FALSE
  )
  return(additional_points)
}





# Generate additional points for all coordinates
additional_points_list <- lapply(1:nrow(distinct_under5child_diarrhea), function(i) {
  cluster_id <- i
  original_point <- data.frame(
    LONGNUM = distinct_under5child_diarrhea$LONGNUM[i], 
    LATNUM = distinct_under5child_diarrhea$LATNUM[i], 
    cluster_id = cluster_id,
    is_original = TRUE
  )
  neighbors <- generate_additional_points(distinct_under5child_diarrhea$LONGNUM[i], distinct_under5child_diarrhea$LATNUM[i], 1.5, cluster_id)
  rbind(original_point, neighbors)
})

# Combine original and additional points
all_points <- do.call(rbind, additional_points_list)

# Convert to spatial points
coordinates <- SpatialPoints(all_points[, c("LONGNUM", "LATNUM")], proj4string = CRS(projection(bioclim_raster)))

# Extract raster values for all points
raster_values <- raster::extract(bioclim_raster, coordinates)
raster_values_df <- data.frame(all_points, raster_values)

# Remove rows with NAs in any of the raster values
raster_values_df <- raster_values_df %>% filter(if_all(starts_with("wc2.1_30s_bio"), ~ !is.na(.)))

# Calculate mean, variance, and standard deviation for each cluster
summary_stats <- raster_values_df %>%
  group_by(cluster_id) %>%
  summarise(across(starts_with("wc2.1_30s_bio"), list(
    mean = ~mean(., na.rm = TRUE),
    variance = ~var(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE)
  )), .groups = 'drop')

# Remove clusters with any NA values in summary statistics
summary_stats <- summary_stats %>% filter(if_all(everything(), ~ !is.na(.)))

# Rename columns in summary_stats to avoid duplicate column names
summary_stats <- summary_stats %>%
  rename_with(~ paste0(.x, "_summary"), starts_with("wc2.1_30s_bio"))

summary(summary_stats)

#Most appropriate variables will be selected based on vairance values,
#collinearity and VIF, and biological raitoal, since we are looking at the 
#long-term impacts of temperature and precipitation on diarrhea we will take variables assessing long-term impacts only



#####Checking variance
# Summary statistics of variance columns
variance_cols <- grep("_variance", names(summary_stats), value = TRUE)
summary(summary_stats[, variance_cols])

# Visualizing distribution of variances using histograms
par(mfrow=c(2, 4))  # Adjust to plot 8 histograms
for (col in variance_cols) {
  hist(summary_stats[[col]], main = paste("Variance:", col), xlab = "Variance", col = "lightblue", breaks = 20)
}
par(mfrow=c(1, 1))  # Reset the plot layout

# Boxplot to visualize potential outliers
boxplot(summary_stats[, variance_cols], main="Variance Boxplots", col="lightblue", outline=TRUE)


# Calculate the interquartile range (IQR) for each variance column
outliers <- summary_stats %>%
  select(cluster_id, all_of(variance_cols)) %>%
  pivot_longer(cols = all_of(variance_cols), names_to = "bioclim_var", values_to = "variance") %>%
  group_by(bioclim_var) %>%
  mutate(
    Q1 = quantile(variance, 0.25),
    Q3 = quantile(variance, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = (variance < lower_bound) | (variance > upper_bound)
  )

# View the outliers
outliers_filtered <- outliers %>% filter(is_outlier == TRUE)
print(outliers_filtered)

# Count of outliers by variable
outliers_filtered %>%
  group_by(bioclim_var) %>%
  summarise(outlier_count = n())


# Remove clusters with outliers
non_outlier_clusters <- outliers %>% filter(is_outlier == FALSE) %>% pull(cluster_id) %>% unique()

# Filter the summary_stats dataset to keep only non-outlier clusters
summary_stats_filtered <- summary_stats %>% filter(cluster_id %in% non_outlier_clusters)

# Visualize the filtered dataset summary
summary(summary_stats_filtered)



# Summary statistics of variance columns
variance_cols <- grep("_variance", names(summary_stats), value = TRUE)
summary(summary_stats[, variance_cols])

# Visualizing distribution of variances using histograms
par(mfrow=c(2, 4))  # Adjust to plot 8 histograms
for (col in variance_cols) {
  hist(summary_stats[[col]], main = paste("Variance:", col), xlab = "Variance", col = "lightblue", breaks = 20)
}
par(mfrow=c(1, 1))  # Reset the plot layout

# Boxplot to visualize potential outliers
boxplot(summary_stats[, variance_cols], main="Variance Boxplots", col="lightblue", outline=TRUE)




# Load necessary libraries
library(car)  # For VIF calculation

# Step 1: Subset the selected Bioclim variables (bio2, bio10, bio11, bio15)
bioclim_selected_vars <- summary_stats %>%
  dplyr::select(wc2.1_30s_bio_1_mean_summary, 
                wc2.1_30s_bio_7_mean_summary,
                wc2.1_30s_bio_12_mean_summary, 
                wc2.1_30s_bio_15_mean_summary,
                wc2.1_30s_bio_17_mean_summary)

# Step 2: Remove rows with NA or infinite values
bioclim_selected_vars_clean <- bioclim_selected_vars %>%
  filter_all(all_vars(!is.na(.))) %>%
  filter_all(all_vars(is.finite(.)))

# Step 3: Fit a linear model with dummy response for VIF calculation (since VIF requires a model context)
dummy_response <- rnorm(nrow(bioclim_selected_vars_clean))  # Create a random dummy response

# Step 4: Fit the model using the selected Bioclim variables as predictors
lm_model <- lm(dummy_response ~ ., data = bioclim_selected_vars_clean)

# Step 5: Calculate VIF for the selected Bioclim variables
vif_values <- vif(lm_model)

# Step 6: Print VIF values
print(vif_values)


##based on above mentioned test, we will now filter the data to selected set fo bioclim variables

# Step 1: Filter summary_stats to keep only relevant Bioclim variables (1, 7, 12, 15, and 17)
summary_stats_filtered <- summary_stats %>%
  select(cluster_id,
         starts_with("wc2.1_30s_bio_1_mean_summary"), starts_with("wc2.1_30s_bio_1_variance_summary"), starts_with("wc2.1_30s_bio_1_sd_summary"),
         starts_with("wc2.1_30s_bio_7_mean_summary"), starts_with("wc2.1_30s_bio_7_variance_summary"), starts_with("wc2.1_30s_bio_7_sd_summary"),
         starts_with("wc2.1_30s_bio_12_mean_summary"), starts_with("wc2.1_30s_bio_12_variance_summary"), starts_with("wc2.1_30s_bio_12_sd_summary"),
         starts_with("wc2.1_30s_bio_15_mean_summary"), starts_with("wc2.1_30s_bio_15_variance_summary"), starts_with("wc2.1_30s_bio_15_sd_summary"),
         starts_with("wc2.1_30s_bio_17_mean_summary"), starts_with("wc2.1_30s_bio_17_variance_summary"), starts_with("wc2.1_30s_bio_17_sd_summary"))

# Step 2: Merge with LATNUM and LONGNUM based on cluster_id
# Assuming that `all_points` contains the original LATNUM, LONGNUM, and cluster_id data

final_data <- all_points %>%
  filter(is_original == TRUE) %>%  # Keep only original points (no additional neighbors)
  select(LATNUM, LONGNUM, cluster_id) %>%  # Select relevant columns
  left_join(summary_stats_filtered, by = "cluster_id")  # Merge summary statistics using cluster_id

# Step 3: View the final data
print(head(final_data))

