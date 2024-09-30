# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(dplyr)

# Ensure all variables are numeric where necessary

# Create the new subset with the selected variables, removing rows where LATNUM or LONGNUM are 0.00
under5children_filter5 <- under5children_filter3 %>%
  select(V000, V001, LATNUM, LONGNUM, 
         HH_members, women_bmi, women_age, women_edu_cont, 
         stunting_cont, underweight_cont, wasting_cont, 
         delivery_by_csection, postnatal_checks, household_head_sex, 
         displace_status, access_healthcare, region, child_gender, 
         drugs_intestinal_parasites, drinking_water, sanitation_facility, 
         anemia_level, diarrhea, 
         wealth_index, birth_size) %>%
  filter(LATNUM != 0.00 & LONGNUM != 0.00)  # Remove rows where LATNUM or LONGNUM are 0.00



#### data processing at cluster level
u5sub <- under5children_filter5
head(u5sub)
table(u5sub$V000)
u5sub$clust <- paste(u5sub$V000,".", u5sub$V001, sep="")
table(u5sub$clust)

# Get the unique LATNUM and LONGNUM for each cluster
latlong_summ <- u5sub %>%
  group_by(clust) %>%
  summarise(LATNUM = first(LATNUM), 
            LONGNUM = first(LONGNUM))

clust.vec <- attr(table(u5sub$clust), "names")
lclust <- length(clust.vec)
clustn <- as.data.frame(table(u5sub$clust))[,2]

# delivery by c-section (Bernoulli)
table(u5sub$delivery_by_csection)
u5sub$delivery_by_csection <- as.numeric(as.character(u5sub$delivery_by_csection))
csectp <- as.vector(xtabs(u5sub$delivery_by_csection ~ u5sub$clust) / clustn)
csectv <- csectp * (1 - csectp) # Bernoulli variance

# drinking water (Bernoulli; 1 = improved; 0 = unimproved)
table(u5sub$drinking_water)
u5sub$dwimp <- ifelse(u5sub$drinking_water == 1, 1, 0)
dwimprp <- as.vector(xtabs(u5sub$dwimp ~ u5sub$clust) / clustn)
dwimprv <- dwimprp * (1 - dwimprp) # Bernoulli variance

# child gender (Bernoulli; 0 = male; 1 = female)
table(u5sub$child_gender)
u5sub$cgend <- ifelse(u5sub$child_gender == 2, 1, 0)
femp <- as.vector(xtabs(u5sub$cgend ~ u5sub$clust) / clustn)
femv <- femp * (1 - femp)


# postnatal_checks (Bernoulli; 0 = no check; 1 = yes check)
table(u5sub$postnatal_checks)
u5sub$pnc <- ifelse(u5sub$postnatal_checks == 1, 1, 0)
pncp <- as.vector(xtabs(u5sub$pnc ~ u5sub$clust) / clustn)
pncv <- pncp * (1 - pncp)

# household_head_sex (Bernoulli; 0 = male; 1 = female)
table(u5sub$household_head_sex)
u5sub$household_head_sex[u5sub$household_head_sex == 3] <- NA
u5sub$hhsex <- ifelse(u5sub$household_head_sex == 2, 1, 0)
hhsexp <- as.vector(xtabs(u5sub$hhsex ~ u5sub$clust) / clustn)
hhsexv <- hhsexp * (1 - hhsexp)

# displacement_status (Bernoulli; 0 = no dplace; 1 = yes dplace)
table(u5sub$displace_status)
u5sub$dplace <- ifelse(u5sub$displace_status == 1, 1, 0)
dplacep <- as.vector(xtabs(u5sub$dplace ~ u5sub$clust) / clustn)
dplacev <- dplacep * (1 - dplacep)

# access_healthcare (Bernoulli; 0 = no access; 1 = yes access)
table(u5sub$access_healthcare)
u5sub$acchc <- ifelse(u5sub$access_healthcare == 1, 1, 0)
acchcp <- as.vector(xtabs(u5sub$acchc ~ u5sub$clust) / clustn)
acchcv <- acchcp * (1 - acchcp)

# region (Bernoulli; 0 = urban; 1 = rural)

u5sub$region[u5sub$region == 7] <- NA
u5sub$regionc <- ifelse(u5sub$region == 1, 1, 0)
clustn_region <- xtabs(~ u5sub$clust) - xtabs(is.na(u5sub$regionc) ~ u5sub$clust)
regionp <- as.vector(xtabs(u5sub$regionc ~ u5sub$clust, na.rm = TRUE) / clustn_region)
regionv <- regionp * (1 - regionp)

# Ensure that regionp is correctly stored in the final dataset
u5sub.clustsumm$regionp <- regionp
u5sub.clustsumm$regionv <- regionv

# Check the distribution of regionp
table(regionp)




# drugs intestinal parasites (Bernoulli; 0 = no drugs; 1 =  yes drugs)
table(u5sub$drugs_intestinal_parasites)
u5sub$drugip <- ifelse(u5sub$drugs_intestinal_parasites == 1, 1, 0)
drugipp <- as.vector(xtabs(u5sub$drugip ~ u5sub$clust) / clustn)
drugipv <- drugipp * (1 - drugipp)

# sanitation_facility (Bernoulli; 0 = unimporved; 1 = improved)
table(u5sub$sanitation_facility)
u5sub$sfimp <- ifelse(u5sub$sanitation_facility == 1, 1, 0)
sfimpp <- as.vector(xtabs(u5sub$sfimp ~ u5sub$clust) / clustn)
sfimpv <- sfimpp * (1 - sfimpp) # Bernoulli variance



# anaemia (recode to binary; (1 = none; 2 = mild; 3 = severe; 99 = NA) => 0 = no anemia, 1 = yes anemia)
table(u5sub$anemia_level)
u5sub$anemia_level[u5sub$anemia_level == 99] <- NA
u5sub$anaemrc <- ifelse(u5sub$anemia_level == 1, 0, 1)
anaemp <- as.vector(xtabs(u5sub$anaemrc ~ u5sub$clust, na.rm=T) / (clustn -
                                                                     as.vector(xtabs(is.na(u5sub$anaemrc) ~ u5sub$clust))))
anaemv <- anaemp * (1 - anaemp)

# diarrhoea (Bernouilli)
table(u5sub$diarrhea)
diarp <- as.vector(xtabs(u5sub$diarrhea ~ u5sub$clust, na.rm=T) / (clustn -
                                                                     as.vector(xtabs(is.na(u5sub$diarrhea) ~ u5sub$clust))))
diarv <- diarp * (1 - diarp)

#ordinal variables
# wealth index (1 = lowest; 2 = 2nd-lowest; 3 = middle; 4 = 2nd-highest; 5 = highest)
#u5sub$wealth_index <- as.ordered(u5sub$wealth_index)
#str(u5sub)
wlthsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(median=median(wealth_index), lo75=quantile(wealth_index, probs=0.125),
            up75=quantile(wealth_index, probs=0.875))
head(wlthsumm)
wlthmd <- as.integer(round(wlthsumm$median, 0))

# for implementation in stochastic resampling:
# wlth.rnd <- rpoistrunc(length(wlthmd), wlthmd, minimum = 1, method="harding")
# wlth.rnd <- ifelse(wlth.rnd > 5, 5, wlth.rnd)

# birth size (1 = small; 2 = average; 3 = large)
table(u5sub$birth_size)
u5sub$birth_size <- as.numeric(as.character(u5sub$birth_size))
birthsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(median=median(birth_size), lo75=quantile(birth_size, probs=0.125),
            up75=quantile(birth_size, probs=0.875))
head(birthsumm)
birthsmd <- as.integer(round(birthsumm$median, 0))
# for implementation in stochastic resampling:
# births.rnd <- rpoistrunc(length(birthsmd), birthsmd, minimum = 1, method="harding")
# births.rnd <- ifelse(births.rnd > 3, 3, wlth.rnd)

#continous variables
# woman bmi (just take average and SD)
hist(u5sub$women_bmi)
wbmisumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(women_bmi), sd=sd(women_bmi))
head(wbmisumm)
wbmimn <- wbmisumm$mean
wbmisd <- wbmisumm$sd

# HH_members
hist(u5sub$HH_members)
hhmemsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(HH_members), sd=sd(HH_members))
head(hhmemsumm)
hhmemmn <- hhmemsumm$mean
hhmemsd <- hhmemsumm$sd


# women_age
hist(u5sub$women_age)
wagesumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(women_age), sd=sd(women_age))
head(wagesumm)
wagemn <- wagesumm$mean
wagesd <- wagesumm$sd



# women_edu_cont
hist(u5sub$women_edu_cont)
wedusumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(women_edu_cont), sd=sd(women_edu_cont))
head(wedusumm)
wedumn <- wedusumm$mean
wedusd <- wedusumm$sd


# stunting_cont (just take average and SD)
hist(u5sub$stunting_cont)
stuntsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(stunting_cont), sd=sd(stunting_cont))
head(stuntsumm)
stuntmn <- stuntsumm$mean
stuntsd <- stuntsumm$sd

# underweight_cont (just take average and SD)
hist(u5sub$underweight_cont)
uweighsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(underweight_cont), sd=sd(underweight_cont))
head(uweighsumm)
uweighmn <- uweighsumm$mean
uweighsd <- uweighsumm$sd


# wasting_cont (just take average and SD)
hist(u5sub$wasting_cont)
wastsumm <- u5sub %>%
  group_by(clust) %>%
  summarise(mean=mean(wasting_cont), sd=sd(wasting_cont))
head(wastsumm)
wastmn <- wastsumm$mean
wastsd <- wastsumm$sd


# combine swomen_bmi# combine summary parameters by cluster
u5sub.clustsumm <- data.frame("clust"=clust.vec, LATNUM = latlong_summ$LATNUM, LONGNUM = latlong_summ$LONGNUM,  csectp, csectv, dwimprp, dwimprv, femp, femv,pncp,pncv, hhsexp,hhsexv, dplacep,dplacev,acchcp,acchcv,regionp,regionv, drugipp,drugipv, sfimpp,sfimpv, anaemp, anaemv,  diarp, diarv, wlthmd, birthsmd, hhmemmn,hhmemsd, wbmimn, wbmisd, wagemn,wagesd, wedumn,wedusd, stuntmn, stuntsd, uweighmn, uweighsd, wastmn, wastsd)
head(u5sub.clustsumm)



###plot variance and probabilities


# Create individual plots
plot1 <- ggplot(u5sub.clustsumm, aes(x = csectp, y = csectv)) +
  geom_point(color = "blue", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "C-Section", x = "Probability", y = "Variance") +
  theme_minimal()

plot2 <- ggplot(u5sub.clustsumm, aes(x = dwimprp, y = dwimprv)) +
  geom_point(color = "green", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Improved Drinking Water", x = "Probability", y = "Variance") +
  theme_minimal()

plot3 <- ggplot(u5sub.clustsumm, aes(x = femp, y = femv)) +
  geom_point(color = "purple", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Female Child", x = "Probability", y = "Variance") +
  theme_minimal()

plot4 <- ggplot(u5sub.clustsumm, aes(x = pncp, y = pncv)) +
  geom_point(color = "orange", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Postnatal Check", x = "Probability", y = "Variance") +
  theme_minimal()

plot5 <- ggplot(u5sub.clustsumm, aes(x = hhsexp, y = hhsexv)) +
  geom_point(color = "red", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Household Head (Female)", x = "Probability", y = "Variance") +
  theme_minimal()

plot6 <- ggplot(u5sub.clustsumm, aes(x = dplacep, y = dplacev)) +
  geom_point(color = "yellow", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Displaced Status", x = "Probability", y = "Variance") +
  theme_minimal()

plot7 <- ggplot(u5sub.clustsumm, aes(x = acchcp, y = acchcv)) +
  geom_point(color = "cyan", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Access to Healthcare", x = "Probability", y = "Variance") +
  theme_minimal()

plot8 <- ggplot(u5sub.clustsumm, aes(x = regionp, y = regionv)) +
  geom_point(color = "pink", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Region (Urban/Rural)", x = "Probability", y = "Variance") +
  theme_minimal()

plot9 <- ggplot(u5sub.clustsumm, aes(x = drugipp, y = drugipv)) +
  geom_point(color = "gray", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Intestinal Parasite Drugs", x = "Probability", y = "Variance") +
  theme_minimal()

plot10 <- ggplot(u5sub.clustsumm, aes(x = sfimpp, y = sfimpv)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Improved Sanitation", x = "Probability", y = "Variance") +
  theme_minimal()

plot11 <- ggplot(u5sub.clustsumm, aes(x = anaemp, y = anaemv)) +
  geom_point(color = "darkblue", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Anemia", x = "Probability", y = "Variance") +
  theme_minimal()

plot12 <- ggplot(u5sub.clustsumm, aes(x = diarp, y = diarv)) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Diarrhea", x = "Probability", y = "Variance") +
  theme_minimal()

# Arrange all plots into a single figure (grid layout)
combined_plot <- grid.arrange(
  plot1, plot2, plot3, plot4,
  plot5, plot6, plot7, plot8,
  plot9, plot10, plot11, plot12,
  ncol = 3, nrow = 4
)

# Print the combined plot
combined_plot

u5sub.clustsumm$regionp <- as.numeric(u5sub.clustsumm$regionp)
u5sub.clustsumm$regionv <- as.numeric(u5sub.clustsumm$regionv)


# Check the structure of both datasets
str(final_data)
str(u5sub.clustsumm)

# Convert LATNUM and LONGNUM to numeric if they are not, to make sure merge is done correctly
final_data$LATNUM <- as.numeric(final_data$LATNUM)
final_data$LONGNUM <- as.numeric(final_data$LONGNUM)

u5sub.clustsumm$LATNUM <- as.numeric(u5sub.clustsumm$LATNUM)
u5sub.clustsumm$LONGNUM <- as.numeric(u5sub.clustsumm$LONGNUM)


# Merge the processed survey data with processed raster data from DHSDiarrProcessing_4 file
# Perform the left join, ensuring exact matching on LATNUM and LONGNUM
Final_Cluster_level_data <- final_data %>%
  left_join(u5sub.clustsumm, by = c("LATNUM", "LONGNUM"))

# View the first few rows of the merged data
head(Final_Cluster_level_data)

# Save the merged data if needed
write.csv(Final_Cluster_level_data, "DHSclusterLevelDiarrData.csv", row.names = FALSE)

#some additional checks
# Check the structure and summary after the join
str(Final_Cluster_level_data)
summary(Final_Cluster_level_data)



# Step 1: Identify clusters with any NAs in the summary statistics
na_clusters <- u5sub.clustsumm %>%
  filter(if_any(everything(), is.na))  # Filters rows where any column has an NA

# Step 2: View the clusters that have NAs
print(na_clusters$clust)

# Optional: View the entire data for clusters with NAs
print(na_clusters)


##check mean and variance for bioclim vars

# List of Bioclim variables (mean and variance pairs)
bioclim_vars <- c("wc2.1_30s_bio_1", "wc2.1_30s_bio_7", "wc2.1_30s_bio_12", 
                  "wc2.1_30s_bio_15", "wc2.1_30s_bio_17")

# Step 1: Select the mean and variance columns for the chosen bioclim variables
bioclim_mean_variance <- Final_Cluster_level_data %>%
  select(LATNUM, LONGNUM, 
         ends_with("mean_summary"), 
         ends_with("variance_summary"))

# Step 2: Reshape the data to long format for easier plotting
bioclim_long <- bioclim_mean_variance %>%
  pivot_longer(cols = contains("wc2.1_30s_bio"), 
               names_to = c("bioclim", "type"), 
               names_pattern = "(wc2.1_30s_bio_\\d+)_(mean_summary|variance_summary)")

# Step 3: Pivot the table wider to separate mean and variance columns
bioclim_long_wide <- bioclim_long %>%
  pivot_wider(names_from = type, values_from = value)

# Step 4: Plot variance vs mean for each Bioclim variable
ggplot(bioclim_long_wide, aes(x = mean_summary, y = variance_summary, color = bioclim)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear fit lines
  labs(title = "Variance vs Mean for Bioclim Variables",
       x = "Mean",
       y = "Variance") +
  theme_minimal() +
  facet_wrap(~ bioclim, scales = "free")  # Create separate plots for each Bioclim variable



save.image()