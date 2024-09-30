# Load required libraries
library(dplyr)

# Read data
# Replace the file path with the actual path where your CSV is stored
allcountries_BR_data <- read.csv("path_to_file/allcountries_BR_data_unprocessed.csv", row.names = NULL)

# Estimate child's age in months
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(estimated_child_age = V008 - B3)

# Categorize children based on estimated age
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(age_categories = cut(estimated_child_age,
                              breaks = c(-1, 5, 11, 17, 23, 35, 47, 60, Inf),
                              labels = c("Infants (up to 5 months)", "Babies (6 to 11 months)", 
                                         "Toddlers (12 to 17 months)", "Young children (18 to 23 months)", 
                                         "Children (24 to 35 months)", "Preschoolers (36 to 47 months)", 
                                         "Early school age (48 to 60 months)", "Older kids")))

# Identify children alive under 5 years old
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(alive_children_under_5 = ifelse(HW1 < 60 & B5 == 1, 1, 
                                         ifelse(HW1 >= 60, 0, NA)))

# Assign urban/rural regions based on V140
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(region = case_when(
    V140 == 1 ~ "Urban", 
    V140 == 2 ~ "Rural", 
    V140 == 7 ~ "Not a resident"
  ))



# Define safe or unsafe drinking water based on V113
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(drinking_water = case_when(
    V113 %in% c(11,12,13,14,21,31,41,51,61,62,71) ~ 1,  # Improved source
    V113 %in% c(32,42,43,44,63,72,92,96,97) ~ 2,  # Poor source
    TRUE ~ NA_integer_))

# Define sanitation facilities and shared/safe status
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(sanitation_facility = case_when(
    V116 %in% c(11,12,13,21,22,41) ~ 1,  # Improved source
    V116 %in% c(14,15,23,31,42,43,44,45,71,96,97) ~ 2,  # Poor source
    TRUE ~ NA_integer_),
    shared = ifelse(V160 %in% c(1, 7), "Shared", "Unshared"),
    safe = ifelse(shared == "Unshared" & sanitation_facility == 1, 1, 2))

# Additional variable transformations for maternal and child health metrics
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(women_weight = V437 / 10,
         women_bmi = V445 / 100,
         cooking_fuel = case_when(V161 %in% c(1:4, 12) ~ 1, V161 %in% c(6:11, 96) ~ 2, V161 == 95 ~ 3, TRUE ~ NA_integer_),
         women_occupation = case_when(V717 == 0 ~ 0, V717 %in% c(1, 2) ~ 1, V717 %in% c(3, 7) ~ 3, V717 == 8 ~ 4, TRUE ~ NA_integer_),
         marital_status = case_when(V501 %in% c(1, 2) ~ 1, TRUE ~ 0),
         access_healthcare = ifelse(V467B == 1 | V467C == 1 | V467D == 1 | V467F == 1, 0, 1),
         total_pregnancies_terminated = ifelse(V228 == 1, 1, 0))

# Miscarriages and stillbirths based on V233 and V228
allcountries_BR_data <- allcountries_BR_data %>%
  mutate(miscarriages = ifelse(V228 == 1 & V233 < 5, 1, 0),
         stillbirths = ifelse(V228 == 1 & V233 >= 5, 1, 0),
         miscarriage_date = ifelse(miscarriages == 1, paste0(V229, "-", V230), NA),
         stillbirth_date = ifelse(stillbirths == 1, paste0(V229, "-", V230), NA))

# Clean and rename variables for vaccinations, stunting, underweight, etc.
allcountries_BR_data <- allcountries_BR_data %>%
  rename(delivery_by_csection = M17, working_mother = V714, ANC_visits = M14, 
         LBW = M18, postnatal_checks = M70, vaccinations = H10, 
         stunting_cont = HW70, underweight_cont = HW71, wasting_cont = HW72)

# Calculate summary statistics for LBW, diarrhea, ARI, child mortality, etc.
total_lbw_cases <- sum(allcountries_BR_data$LBW == 1, na.rm = TRUE)
diarrhea_stats <- allcountries_BR_data %>%
  group_by(country_code) %>%
  summarise(total_diarrhea_cases = sum(diarrhea == 1, na.rm = TRUE))


# View outputs
View(diarrhea_stats)

