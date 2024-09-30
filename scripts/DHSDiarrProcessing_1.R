# Reading and processing DHS data for child health-climate relationships
# Hira Fatima & Corey Bradshaw
# September 2024

# 1. we will define the list of all plausible vairables for the analysis based on previous extensive literature review and desinged hypothesis
# 2. we will read the files for survey data files and GPS data
# 3. we will merge each survey file with its GPS file using cluster number variable
# 4. we will create dummy variables for DHS 6 or 7 or 8, where vairable was not present in any of these, assigning NA values
# 5. cleanup the environment to keep necessary datasets

# Load required packages
library(foreign)
library(haven)
library(sf)
library(sp)
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)
library(leaflet)
library(mice)


# Define common variables for all datasets
common_vars <- c(
  "CASEID", "BORD", "V000", "V001", "V002", "V003", "V005", "V008", "V012", "V024", 
  "V025", "V127", "V155", "V140", "V104", "V106", "V107", "V151", "V701", "V501", 
  "V714", "V113", "V116", "V160", "V161", "V136", "V137", "V190", "V437", "V445", 
  "V463Z", "V464", "V465", "V717", "V705", "V467B", "V467C", "V467D", "V467F", 
  "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "V201", "V206", "V207", 
  "V208", "V209", "V218", "V219", "V228", "V229", "V230", "V231", "V233", "V234", 
  "V238", "M14", "M15", "M17", "M70", "M18", "M19", "H0", "H2", "H3", "H4", "H5", 
  "H6", "H7", "H8", "H9", "H10", "H51", "H52", "H53", "H54", "H55", "H56", "H57", 
  "H58", "H59", "H33", "H43", "H11", "H22", "H31", "H31B", "H31C", "H42", "HW1", 
  "HW2", "HW3", "HW56", "HW57", "HW70", "HW71", "HW72"
)

# Load datasets with error handling
BR_PK17 <- try(read_sas("PKBR71FL.SAS7BDAT"), silent = TRUE)
BR_BD17 <- try(read_sas("BDBR7RFL.SAS7BDAT"), silent = TRUE)
BR_BD14 <- try(read_sas("BDBR72FL.SAS7BDAT"), silent = TRUE)
BR_BD11 <- try(read_sas("BDBR61FL.SAS7BDAT"), silent = TRUE)
BR_KH21 <- try(read_sas("KHBR82FL.SAS7BDAT"), silent = TRUE)
BR_KH14 <- try(read_sas("KHBR73FL.SAS7BDAT"), silent = TRUE)
BR_KH10 <- try(read_sas("KHBR61FL.SAS7BDAT"), silent = TRUE)
BR_IA19 <- try(read.csv("IABR7EFL.SAS7BDAT"), silent = TRUE)
BR_IA15 <- try(read.csv("IABR74FL.SAS7BDAT"), silent = TRUE)
BR_MM15 <- try(read_sas("MMBR71FL.SAS7BDAT"), silent = TRUE)
BR_NP22 <- try(read_sas("NPBR82FL.SAS7BDAT"), silent = TRUE)
BR_NP16 <- try(read_sas("NPBR7HFL.SAS7BDAT"), silent = TRUE)
BR_NP11 <- try(read_sas("NPBR61FL.SAS7BDAT"), silent = TRUE)
BR_PH22 <- try(read_sas("PHBR82FL.SAS7BDAT"), silent = TRUE)
BR_PH17 <- try(read_sas("PHBR71FL.SAS7BDAT"), silent = TRUE)
BR_TL16 <- try(read_sas("TLBR71FL.SAS7BDAT"), silent = TRUE)
BR_TL10 <- try(read_sas("TLBR61FL.SAS7BDAT"), silent = TRUE)

# set base directory path
base_dir <- "/path/to/your/DHS_data/DHS_allcountries/"

# Define file paths fr GPS data using file.path()
file_paths <- c(
  file.path(base_dir, "BD_2011_DHS_04042024_20_194714/BDGE61FL/BDGE61FL.shp"),
  file.path(base_dir, "BD_2014_DHS_04042024_158_194714/BDGE71FL/BDGE71FL.shp"),
  file.path(base_dir, "BD_2017-18_DHS_04042024_156_194714/BDGE7SFL/BDGE7SFL.shp")
  # Add other paths as needed
)


gps_names <- c(
  "GPS_BD11", "GPS_BD14", "GPS_BD17", # Add other names as needed
)

# Read shapefiles and assign names
for (i in seq_along(file_paths)) {
  assign(gps_names[i], st_read(file_paths[i]))
}

# Define a list of GPS datasets
gps_datasets <- list(GPS_BD11, GPS_BD14, GPS_BD17) # Add other GPS datasets as needed

# Merge datasets with GPS data
merged_datasets <- lapply(gps_datasets, function(gps_data) {
  left_join(BR_PK17, gps_data, by = c("V001" = "DHSCLUST"))
})

# Mutate necessary variables to create dummy columns for each dataset
BR_BD14_merged <- BR_BD14 %>%
  mutate(H51 = NA, H52 = NA, H53 = NA, H54 = NA, H55 = NA, H56 = NA, H57 = NA, H58 = NA, H59 = NA)

BR_BD11_merged <- BR_BD11 %>%
  mutate(H51 = NA, H52 = NA, H53 = NA, H54 = NA, H55 = NA, H56 = NA, H57 = NA, H58 = NA, H59 = NA)

# Subset datasets to include only common variables
subsetted_datasets <- lapply(merged_datasets, function(dataset) {
  select(dataset, all_of(common_vars))
})

# Combine all subsetted datasets into one
allcountries_BR_data <- bind_rows(subsetted_datasets)

# Save the final dataset
write.csv(allcountries_BR_data, "allcountries_BR_data_unprocessed.csv", row.names = FALSE)

# Clean up environment, keep only necessary datasets
rm(list = setdiff(ls(), "allcountries_BR_data"))

# View the first few rows of the combined dataset
head(allcountries_BR_data)

# Save the R workspace
save.image()
