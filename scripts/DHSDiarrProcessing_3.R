# Load necessary libraries
library(dplyr)
library(mice)
library(ggplot2)
library(gridExtra)

# Filter data for children who are alive and under 5 years old
under5children_filter3 <- allcountries_BR_data %>%
  filter(age_categories != "Older kids", B5 == 1)

# Create new variable for the probability of diarrhea by cluster (V000 and V001)
under5children_filter3 <- under5children_filter3 %>%
  group_by(V000, V001) %>%
  mutate(prob_diarrhea = mean(diarrhea, na.rm = TRUE)) %>%
  ungroup()

# Convert placeholders to NA for consistency
under5children_filter3 <- under5children_filter3 %>%
  mutate(
    women_weight = ifelse(women_weight > 999, NA, women_weight),
    women_bmi = ifelse(women_bmi > 99, NA, women_bmi),
    child_weight = ifelse(child_weight > 99, NA, child_weight),
    child_height = ifelse(child_height > 999, NA, child_height),
    vitamin_a = ifelse(vitamin_a == 8, NA, vitamin_a),
    iron_supplement = ifelse(iron_supplement == 8, NA, iron_supplement),
    vaccinations = ifelse(vaccinations == 8, NA, vaccinations),
    postnatal_checks = ifelse(postnatal_checks == 8, NA, postnatal_checks),
    stunting_cont = ifelse(stunting_cont > 9000, NA, stunting_cont),
    underweight_cont = ifelse(underweight_cont > 9000, NA, underweight_cont),
    wasting_cont = ifelse(wasting_cont > 9000, NA, wasting_cont)
  )

#Keep orignal var before imputation
under5children_original <- under5children_filter3



#FOR ALL VARIABLES

#imputing missing values
if (!is.data.frame(under5children_filter3)) {
  under5children_filter3 <- as.data.frame(under5children_filter3)
}


# Load necessary library
library(mice)

# Create a copy of the original data
imputation_data <- standard_data[, sapply(standard_data, is.numeric)]

# Identify variables with missing values
vars_with_missing <- c("drinking_water", "sanitation_facility", "displace_status", "access_healthcare", "women_bmi", "delivery_by_csection", "postnatal_checks", 
                       "place_of_delivery",  "stunting_cont", "underweight_cont", "wasting_cont", "drugs_intestinal_parasites", "birth_size")

# Define your categorical and continuous variables
categorical_vars <- c("drinking_water", "sanitation_facility", "displace_status", "access_healthcare",  "delivery_by_csection", "postnatal_checks", 
                      "place_of_delivery",   "drugs_intestinal_parasites", "birth_size")
continuous_vars <- c("women_bmi", "stunting_cont", "underweight_cont", "wasting_cont")

# Ensure all categorical variables are factors
for (var in categorical_vars) {
  if (var %in% names(imputation_data) && length(imputation_data[[var]]) > 0) {
    imputation_data[[var]] <- as.factor(imputation_data[[var]])
  } else {
    warning(paste("Variable", var, "is missing or empty."))
  }
}

# Ensure all continuous variables are numeric
for (var in continuous_vars) {
  if (var %in% names(imputation_data) && length(imputation_data[[var]]) > 0) {
    imputation_data[[var]] <- as.numeric(imputation_data[[var]])
  } else {
    warning(paste("Variable", var, "is missing or empty."))
  }
}


# Subset the data to include only variables of interest
subset_data <- imputation_data[c(continuous_vars, categorical_vars)]

# Specify the methods for imputation for each variable
methods <- make.method(subset_data)
methods[categorical_vars] <- "polyreg"  # Polytomous regression for categorical variables
methods[continuous_vars] <- "pmm"       # Predictive mean matching for continuous variables

# Check predictor matrix
predictorMatrix <- make.predictorMatrix(subset_data)
print(predictorMatrix)

# Perform imputation with detailed warning handling
mice_object <- tryCatch(
  {
    mice(subset_data, method = methods, predictorMatrix = predictorMatrix, m = 1, printFlag = FALSE)
  },
  warning = function(w) {
    print("Warning occurred during imputation:")
    print(w)
    return(mice(subset_data, method = methods, predictorMatrix = predictorMatrix, m = 1, printFlag = FALSE, maxit = 0))  # Return mice object with diagnostics
  },
  error = function(e) {
    print("Error occurred during imputation:")
    print(e)
    return(NULL)
  }
)

if (!is.null(mice_object)) {
  # Print diagnostics if warnings were captured
  if (inherits(mice_object, "mids")) {
    warning_messages <- mice_object$loggedEvents
    print("Detailed warning messages during imputation:")
    print(warning_messages)
  }
  
  # Extract imputed data
  imputed_data <- complete(mice_object)
  
  # Merge imputed values back into the original dataset
  for (var in vars_with_missing) {
    if (var %in% names(under5children_filter3) && var %in% names(imputed_data)) {
      under5children_filter3[[var]] <- imputed_data[[var]]
    } else {
      warning(paste("Variable", var, "is missing in the original or imputed dataset."))
    }
  }
  
  # View the updated dataset
  head(under5children_filter3)
}




#FOR ANEMIA_LEVEL


# Ensure the dataset is a data frame
if (!is.data.frame(under5children_filter3)) {
  under5children_filter3 <- as.data.frame(under5children_filter3)
}

# Step 1: Filter the data to exclude rows where anemia_level is equal to 99, keeping NA values
filtered_data <- under5children_filter3 %>% filter(anemia_level != 99 | is.na(anemia_level))

# Select the anemia_level and predictor variables
imputation_data <- filtered_data[, c("anemia_level", "child_weight", "child_height", "hemoglobin_level", "underweight", "Wasting", "Stunting")]

# Ensure anemia_level is a factor for polyreg imputation
imputation_data$anemia_level <- as.factor(imputation_data$anemia_level)

# Step 2: Specify the methods for imputation
methods <- make.method(imputation_data)
methods["anemia_level"] <- "polyreg"  # Polytomous regression for anemia_level
methods["child_weight"] <- ""  # No imputation for predictor variables
methods["child_height"] <- ""
methods["hemoglobin_level"] <- ""
methods["underweight"] <- ""
methods["Wasting"] <- ""
methods["Stunting"] <- ""

# Create the predictor matrix, ensuring only anemia_level is imputed
predictorMatrix <- make.predictorMatrix(imputation_data)
predictorMatrix[, "anemia_level"] <- 1  # Use all variables to predict anemia_level
predictorMatrix["anemia_level", ] <- 0  # anemia_level itself is not used to predict other variables

# Perform imputation
mice_object <- mice(imputation_data, method = methods, predictorMatrix = predictorMatrix, m = 1, printFlag = FALSE)

# Extract imputed data
imputed_data <- complete(mice_object)

# Replace imputed values back into the filtered data
filtered_data$anemia_level <- imputed_data$anemia_level

# Step 3: Create a logical index of the rows to update in the original data
rows_to_update <- which(under5children_filter3$anemia_level != 99 | is.na(under5children_filter3$anemia_level))

# Ensure there are no NAs in the index
rows_to_update <- rows_to_update[!is.na(rows_to_update)]

# Replace the values in the original dataset with the imputed values
under5children_filter3$anemia_level[rows_to_update] <- filtered_data$anemia_level

# Verify the unique values in anemia_level
unique_anemia_level <- table(under5children_filter3$anemia_level, useNA = "always")
print(unique_anemia_level)

# Check logged events in the mice object
if (!is.null(mice_object)) {
  logged_events <- mice_object$loggedEvents
  print("Detailed warning messages during imputation:")
  print(logged_events)
}




#FOR ROTAVIRUS

# Load necessary libraries
library(mice)
library(dplyr)

# Ensure the dataset is a data frame
if (!is.data.frame(under5children_filter3)) {
  under5children_filter3 <- as.data.frame(under5children_filter3)
}

# Step 1: Filter the data to exclude rows where received_rotavirus is equal to 3, keeping NA values
filtered_data <- under5children_filter3 %>% filter(received_rotavirus != 3 | is.na(received_rotavirus))

# Select the received_rotavirus and predictor variables
imputation_data <- filtered_data[, c("received_rotavirus", "child_age", "child_height", "Stunting", "Wasting", 
                                     "received_rotavirus1", "received_rotavirus2", "received_polio1", "received_polio2", "vaccinations")]

# Ensure received_rotavirus is a factor for polyreg imputation
imputation_data$received_rotavirus <- as.factor(imputation_data$received_rotavirus)

# Specify the methods for imputation
methods <- make.method(imputation_data)
methods["received_rotavirus"] <- "polyreg"  # Polytomous regression for received_rotavirus
methods["child_age"] <- ""  # No imputation for predictor variables
methods["child_height"] <- ""
methods["Stunting"] <- ""
methods["Wasting"] <- ""
methods["received_rotavirus1"] <- ""
methods["received_rotavirus2"] <- ""
methods["received_polio1"] <- ""
methods["received_polio2"] <- ""
methods["vaccinations"] <- ""

# Create the predictor matrix, ensuring only received_rotavirus is imputed
predictorMatrix <- make.predictorMatrix(imputation_data)
predictorMatrix[, "received_rotavirus"] <- 1  # Use all variables to predict received_rotavirus
predictorMatrix["received_rotavirus", ] <- 0  # received_rotavirus itself is not used to predict other variables

# Perform imputation
mice_object <- mice(imputation_data, method = methods, predictorMatrix = predictorMatrix, m = 1, printFlag = FALSE)

# Extract imputed data
imputed_data <- complete(mice_object)

# Replace imputed values back into the filtered data
filtered_data$received_rotavirus <- imputed_data$received_rotavirus

# Step 3: Create a logical index of the rows to update in the original data
rows_to_update <- which(under5children_filter3$received_rotavirus != 3 | is.na(under5children_filter3$received_rotavirus))

# Ensure there are no NAs in the index
rows_to_update <- rows_to_update[!is.na(rows_to_update)]

# Replace the values in the original dataset with the imputed values
under5children_filter3$received_rotavirus[rows_to_update] <- filtered_data$received_rotavirus

# Verify the unique values in received_rotavirus
unique_received_rotavirus <- table(under5children_filter3$received_rotavirus, useNA = "always")
print(unique_received_rotavirus)

# Check logged events in the mice object
if (!is.null(mice_object)) {
  logged_events <- mice_object$loggedEvents
  print("Detailed warning messages during imputation:")
  print(logged_events)
}




#FOR REGION

# Load necessary libraries
library(mice)
library(dplyr)

# Ensure the dataset is a data frame
if (!is.data.frame(under5children_filter3)) {
  under5children_filter3 <- as.data.frame(under5children_filter3)
}

# Step 1: Filter the data to exclude rows where region is equal to 7, keeping NA values
filtered_data_region <- under5children_filter3 %>% filter(region != 7 | is.na(region))

# Select the region and predictor variables
imputation_data_region <- filtered_data_region[, c("region", "women_education", "working_mother", "wealth_index", "HH_members")]

# Ensure region is a factor for polyreg imputation
imputation_data_region$region <- as.factor(imputation_data_region$region)

# Specify the methods for imputation
methods_region <- make.method(imputation_data_region)
methods_region["region"] <- "polyreg"  # Polytomous regression for region
methods_region["women_education"] <- ""  # No imputation for predictor variables
methods_region["working_mother"] <- ""
methods_region["wealth_index"] <- ""
methods_region["HH_members"] <- ""

# Create the predictor matrix, ensuring only region is imputed
predictorMatrix_region <- make.predictorMatrix(imputation_data_region)
predictorMatrix_region[, "region"] <- 1  # Use all variables to predict region
predictorMatrix_region["region", ] <- 0  # region itself is not used to predict other variables

# Perform imputation
mice_object_region <- mice(imputation_data_region, method = methods_region, predictorMatrix = predictorMatrix_region, m = 1, printFlag = FALSE)

# Extract imputed data
imputed_data_region <- complete(mice_object_region)

# Replace imputed values back into the filtered data
filtered_data_region$region <- imputed_data_region$region

# Step 3: Create a logical index of the rows to update in the original data
rows_to_update_region <- which(under5children_filter3$region != 7 | is.na(under5children_filter3$region))

# Ensure there are no NAs in the index
rows_to_update_region <- rows_to_update_region[!is.na(rows_to_update_region)]

# Replace the values in the original dataset with the imputed values
under5children_filter3$region[rows_to_update_region] <- filtered_data_region$region

# Verify the unique values in region
unique_region <- table(under5children_filter3$region, useNA = "always")
print(unique_region)

# Check logged events in the mice object
if (!is.null(mice_object_region)) {
  logged_events <- mice_object_region$loggedEvents
  print("Detailed warning messages during imputation:")
  print(logged_events)
}


#CREATE CSV for under5children_filter 3 file which will be further used for other analysis as well

write.csv(under5children_filter3, "under5children_imputed.csv", row.names = FALSE)

#data visualization
#create histograms for different variables to see their performance after imputation


# Generate histograms for continuous and categorical variables after imputation
generate_histogram_continuous <- function(data, variable) {
  ggplot(data, aes(x = .data[[variable]])) + 
    geom_histogram(bins = 30, fill = 'blue', alpha = 0.5) + 
    ggtitle(paste(variable, "Distribution")) + 
    theme_minimal()
}

generate_histogram_categorical <- function(data, variable) {
  ggplot(data, aes(x = .data[[variable]])) + 
    geom_bar(fill = 'blue', alpha = 0.5) + 
    ggtitle(paste(variable, "Distribution")) + 
    theme_minimal()
}

# Example for generating histograms for BMI and drinking water
p1 <- generate_histogram_continuous(under5children_filter3, "women_bmi")
p2 <- generate_histogram_categorical(under5children_filter3, "drinking_water")

# Display histograms
grid.arrange(p1, p2, ncol = 2)

# Clean up unnecessary datasets and save workspace
rm(list = setdiff(ls(), c("under5children_filter3", "allcountries_BR_data", "under5children_original")))
save.image()

