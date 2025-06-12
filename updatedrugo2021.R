#SCRIPT 2

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(readxl)       # Reading Excel files
library(mgcv)         # Generalized Additive Models (GAMs)
library(performance)  # Model diagnostics
library(ggeffects)    # Predictions and effect plots
library(ggpubr)       # Arranging ggplots
library(broom)        # Tidy model output
library(MASS)         # For Negative Binomial GLM
library(visreg)
library(DHARMa)
library(lsmeans)
library(emmeans)
library(multcomp)
library(multcompView)
library(MuMIn)
library(glmmTMB)

#Christine Model 2 = quantify relationship with rugosity and location

##########################Data processing#######################################

# # Part 1: Combining CSV files from all sites
# folder_path <- "data/rugo csv"
# 
# # Get a list of all CSV files in the folder
# file_list <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
# 
# # Initialize an empty list to store data frames
# data_list <- list()
# 
# # Loop through each file, read it, and ensure consistent column names
# for (file in file_list) {
#   data <- read.csv(file)  # Read the file
#   colnames(data) <- tolower(colnames(data))  # Standardize column names
#   data_list[[length(data_list) + 1]] <- data  # Append data to list
# }
# 
# # Combine all data frames into one
# rug_combined_data <- bind_rows(data_list)

# Save the combined dataset to a CSV (optional)
#write.csv(rug_combined_data, file = "~/Desktop/Rproject/rugosity.csv", row.names = FALSE)

#################################################################################

rug_combined_data = read_csv("data/old/combined_sites.csv")

# Part 2: Calculate rugosity index for each transect
rugosity <- rug_combined_data %>%
  filter(z != 0) %>%  # Exclude points where z = 0
  group_by(site, transect) %>%
  summarize(
    diff_x = diff(x),
    diff_y = diff(y),
    diff_z = diff(z),
    standard_length = sum(sqrt(diff_x^2 + diff_y^2)),  # 2D Euclidean distance
    true_length = sum(sqrt(diff_x^2 + diff_y^2 + diff_z^2)),  # 3D distance
    rugosity_index = ifelse(standard_length > 0, true_length / standard_length, NA)  # Rugosity calculation
  ) %>%
  ungroup()

# Histogram to visualize rugosity distribution
hist(rugosity$rugosity_index)

# Part 3: Calculate average rugosity index for each site, handling NAs
Avg_rugosity_per_site <- rugosity %>%
  group_by(site) %>%
  summarize(
    Avg_rugosity = mean(rugosity_index, na.rm = TRUE),    # Average rugosity
    rugosity_variance = var(rugosity_index, na.rm = TRUE),  # Variance per site
  )

# View the calculated rugosity indices per transect
head(Avg_rugosity_per_site)

# Separate 'site' column into 'site', 'plot', and 'timepoint'
Avg_rugosity_per_site <- Avg_rugosity_per_site %>%
  separate(site, into = c("site", "plot", "timepoint"), sep = "-", remove = TRUE, fill = "right")

# Save processed data
#write.csv(rugosity, file = "data/rugosity_transects.csv", row.names = FALSE)
#write.csv(Avg_rugosity_per_site, file = "data/Avg_rugosity_sites.csv", row.names = FALSE)

# Ensure column names match for joins
Avg_rugosity_per_site <- Avg_rugosity_per_site %>% rename_with(~ stringr::str_to_title(.))

# Print the summary of average rugosity per site
print(Avg_rugosity_per_site)

#join rugosity-site data to snapper and grouper MaxN and Occurrence respectively **new_gs2 file
new_gs2 <- new_gs2 %>% mutate(Plot = as.character(Plot))
new_gs2$Site <- trimws(new_gs2$Site) # Trim whitespace from both Site and Plot columns
new_gs2$Plot <- trimws(new_gs2$Plot)

# Create SitePlot column by combining Site and Plot
Avg_rugosity_per_site <- Avg_rugosity_per_site %>%
  mutate(SitePlot = paste(Site, Plot, sep = ""))

# SitePlotNew inserted into table: The site plots did not match
newkeys<- read_excel("data/old/Trip_Records.xlsx", sheet = "All sites")
newkeys <- subset(newkeys, select = c(SitePlotOld, SitePlotNew))

write_csv(newkeys, "data/newkeys.csv")

# Merge new site plots into data set
new_gs2_joined <- new_gs2 %>%
  left_join(newkeys, by = c("SitePlot" = "SitePlotOld"))

new_gs2_joined <- new_gs2_joined %>%
  left_join(
    subset(Avg_rugosity_per_site, select = c(SitePlot))
    , by = c("SitePlotNew" = "SitePlot"))

#write.csv(new_gs2_joined, file = "data/new_gs2_joined.csv", row.names = FALSE)
#yes


# Start Here --------------------------------------------------------------


new_gs2_joined = read_csv("data/old/new_gs2_joined.csv")
snapper_oc = read_csv("data/old/snapper_oc.csv")
grouper_oc = read_csv("data/old/grouper_oc.csv")

unique(grouper_oc$Common)

# Join rugosity data with occurence datasets
snapper_oc_joined <- snapper_oc %>%  left_join(newkeys, by = c("SitePlot" = "SitePlotOld"))
grouper_oc_joined <- grouper_oc %>%  left_join(newkeys, by = c("SitePlot" = "SitePlotOld"))

#remove MaxN since it adds 0s to the data.
snapper_oc_rugo_joined <- snapper_oc_joined %>%
  left_join(
    subset(Avg_rugosity_per_site, select = c(SitePlot,Avg_rugosity,Rugosity_variance))
    , by = c("SitePlotNew" = "SitePlot"))

grouper_oc_rugo_joined <- grouper_oc_joined %>%
  left_join(
    subset(Avg_rugosity_per_site, select = c(SitePlot,Avg_rugosity,Rugosity_variance))
    , by = c("SitePlotNew" = "SitePlot"))

# Save joined datasets
#write.csv(snapper_oc_rugo_joined, file = "data/snapper_oc_rugo_joined.csv", row.names = FALSE)
#write.csv(grouper_oc_rugo_joined, file = "data/grouper_oc_rugo_joined.csv", row.names = FALSE)

# Remove rows with missing values in relevant columns
snapper_oc_rugo_joined <- snapper_oc_rugo_joined %>%
  filter(
    !is.na(Avg_rugosity),
    !is.na(SitePlotNew),
    !is.na(Rugosity_variance)
  )

grouper_oc_rugo_joined <- grouper_oc_rugo_joined %>%
  filter(
    !is.na(Avg_rugosity),
    !is.na(SitePlotNew),
    !is.na(Rugosity_variance)
  )

# Prerpare datasets for modeling.
snapper_abundance_r <- snapper_oc_rugo_joined %>% #Drop PA column
  dplyr::select(SitePlot, Site, Plot, Location, Common, MaxN, SitePlotNew, Avg_rugosity, Rugosity_variance) %>%
  filter(MaxN != 0) |> 
  filter(Common %in% c("Yellowtail","Schoolmaster"))

write_csv(snapper_abundance_r, "data/snapper_abundance_r.csv")

snapper_occurrence_r <- snapper_oc_rugo_joined %>% #Drop MaxN column
  dplyr::select(SitePlot, Site, Plot, Location, Common, PA, SitePlotNew, Avg_rugosity, Rugosity_variance)|> 
  filter(Common %in% c("Yellowtail","Schoolmaster"))

write_csv(snapper_occurrence_r, "data/snapper_occurrence_r.csv")

grouper_abundance_r <- grouper_oc_rugo_joined %>% #Drop PA column
  dplyr::select(SitePlot, Site, Plot, Location, Common, MaxN, SitePlotNew, Avg_rugosity, Rugosity_variance) %>%
  filter(MaxN != 0)|> 
  filter(Common %in% c("Graysby","Red Hind"))

write_csv(grouper_abundance_r, "data/grouper_abundance_r.csv")

grouper_occurrence_r <- grouper_oc_rugo_joined %>% #Drop MaxN column
  dplyr::select(SitePlot, Site, Plot, Location, Common, PA, SitePlotNew, Avg_rugosity, Rugosity_variance)|> 
  filter(Common %in% c("Graysby","Red Hind"))
write_csv(grouper_occurrence_r, "data/grouper_occurrence_r.csv")

###########################################################################
#### Model Check:
# TO DO: Run GLM/Neg Binomial, compare NB vs Poisson, check model structure
# Plan: Pairwise test to identify mean groups, gg predict, ggplot, ggarrange

# Model 1: MaxN ~ Rugosity * Location (Stratified by "Common" Species)
#Snapper abundance: MaxN
snapp_maxn_pois <- glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = poisson, data = snapper_abundance_r)
snapp_maxn_nb <- glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = nbinom2,data = snapper_abundance_r)

check_model(snapp_maxn_pois)
check_model(snapp_maxn_nb)
compare_performance(snapp_maxn_pois, snapp_maxn_nb) 

# Compare models using AIC
AIC(snapp_maxn_pois, snapp_maxn_nb) #NB lower AIC
check_overdispersion(snapp_maxn_nb)  # no overdispersion

summary(snapp_maxn_nb)     # Model coefficients
car::Anova(snapp_maxn_nb)  # ANOVA for significance
performance::check_collinearity(snapp_maxn_nb)  # high corelation.

options(na.action = "na.fail")
dredge(snapp_maxn_nb) # Model selection, looks like location is unnecessaty

snapp_maxn_nb <- glmmTMB(MaxN ~ Avg_rugosity  + Common  + (1 | Site), family = nbinom2,data = snapper_abundance_r)
check_model(snapp_maxn_nb)
car::Anova(snapp_maxn_nb) 
performance::check_collinearity(snapp_maxn_nb)


#I did also run an interaction between Avg_rugosity and common but it was not signifcant and did not inform on anything more
# em_snapp_maxn <- emmeans(snapp_maxn_nb, ~ Common | Avg_rugosity, 
#                          at = list(Avg_rugosity = c(1.0, 1.5, 2.0, 2.5))) 
# pairs(em_snapp_maxn, simple = "Common")


# Snapper model prediction
snapp_maxn_commonpreds <- ggpredict(snapp_maxn_nb, terms = c("Common")) %>% 
  rename(Common = x, MaxN = predicted) 


em_snapp_maxn <- emmeans(snapp_maxn_nb, ~ Common)
pairs(em_snapp_maxn, simple = "Common")

cld_snapp_maxn <- cld(em_snapp_maxn, 
                      Letters = letters,
                      adjust = "bonf", 
                      sort = FALSE, 
                      alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

gg.snapp_maxn = left_join(snapp_maxn_commonpreds, cld_snapp_maxn, by = c("Common"))

ggplot(gg.snapp_maxn, aes(x = Common, y = MaxN)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = code), vjust = -0.5, size = 3) +
  labs(title = "Snapper Abundance (MaxN) by Species",
       x = "Species",
       y = "MaxN") +
  theme_classic()

snapp_maxn_rugopreds = ggpredict(snapp_maxn_nb, terms = c("Avg_rugosity [1:2.5 by = 0.05]")) %>% 
  rename(Avg_rugosity = x, MaxN = predicted)

plot(snapp_maxn_rugopreds) +
  theme_classic()


##Grouper abundance: MaxN
grouper_maxn_pois <- glm(MaxN ~ Avg_rugosity * Location * Common, family = poisson, data = grouper_abundance_r)
grouper_maxn_nb   <- glm.nb(MaxN ~ Avg_rugosity * Location * Common, data = grouper_abundance_r)

grouper_maxn_pois1 = glmmTMB(MaxN ~ Avg_rugosity * Location * Common  * (1 | Site), family = poisson, data = grouper_abundance_r)
grouper_maxn_nb1 = glmmTMB(MaxN ~ Avg_rugosity * Location * Common  * (1 | Site), family = nbinom2,data = grouper_abundance_r)
grouper_maxn_pois2 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = poisson, data = grouper_abundance_r)
grouper_maxn_nb2 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = nbinom2,data = grouper_abundance_r)

# Model diagnostics.
plot(grouper_maxn_pois)   # Check Poisson model structure
plot(grouper_maxn_nb)     # Check NB model structure

# Compare models using AIC
AIC(grouper_maxn_pois, maxn_nb) #NB lower AIC
check_overdispersion(grouper_maxn_nb)  # no overdispersion

summary(grouper_maxn_nb)     # Model coefficients
anova(grouper_maxn_nb, test = "Chisq")  # ANOVA for significance
performance::check_collinearity(grouper_maxn_nb)  # Check multicollinearity. high correlation.

# grouper model prediction
grouper_maxn_preds <- ggpredict(grouper_maxn_nb, terms = c("Common")) %>% 
  rename(Common = x, MaxN = predicted)

#help!
# Run pairwise comparisons using emmeans

# Get compact letter display (CLD) for grouper

# Merge CLD with predictions


# Create the grouper abundance (MaxN) bar graph

###########################################################################
# Model 2: Occurrence ~ Rugosity * Location (Stratified by "Common" Species)

#Snapper occurrence: PA
snapp_occ_pois <- glm(PA ~ Avg_rugosity * Location * Common, family = poisson, data = snapper_occurrence_r)
snapp_occ_nb   <- glm.nb(PA ~ Avg_rugosity * Location * Common, data = snapper_occurrence_r)

compare_performance(snapp_occ_pois, snapp_occ_nb) #POIS lower

# Model diagnostics.
plot(snapp_occ_pois)   # Check Poisson model structure
check_overdispersion(snapp_occ_pois)  # no overdispersion
summary(snapp_occ_pois)     # Model coefficients
anova(snapp_occ_pois, test = "Chisq")  # ANOVA for significance
performance::check_collinearity(snapp_occ_pois)  # high corelation.

# Snapper model prediction
snapp_occ_preds <- ggpredict(snapp_occ_pois, terms = c("Avg_rugosity [mean]", "Common")) %>% 
  rename(Common = group, Occurrence = predicted) %>% 
  mutate(Common = as.factor(Common))


# Run pairwise comparisons using emmeans


# Get compact letter display (CLD) for grouper

#HELP! :( where do I go from here??

#Grouper occurrence: PA
grouper_occ_pois <- glm(PA ~ Avg_rugosity * Location * Common, family = poisson, data = grouper_occurrence_r)
grouper_occ_nb   <- glm.nb(PA ~ Avg_rugosity * Location * Common, data = grouper_occurrence_r)

compare_performance(grouper_occ_pois, grouper_occ_nb)

# Model diagnostics.
plot(grouper_occ_pois)   # Check Poisson model structure
check_overdispersion(grouper_occ_pois)  # no overdispersion
summary(grouper_occ_pois)     # Model coefficients
anova(grouper_occ_pois, test = "Chisq")  # ANOVA for significance
performance::check_collinearity(grouper_occ_pois)  # Check multicollinearity. high correlation.

# grouper model prediction
grouper_occ_preds <- ggpredict(grouper_occ_pois, terms = c("Avg_rugosity", "Location", "Common"))

# Run pairwise comparisons using emmeans


# Get compact letter display (CLD) for grouper

#help! how do I graph predicted values showing the species vs occurrence vs rugosity vs location?

### COMBINE THE PLOTS ###


