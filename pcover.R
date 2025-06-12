#SCRIPT 3:
#Pcover data

##########################Data processing#######################################

# Read the pcover data from CSV
pcover <- read_csv("~/Desktop/Rproject/pcover.csv")

# Inspect and clean the pcover data if needed
pcover <- pcover %>%
  mutate(
    Site = trimws(Site),
    Plot = trimws(Plot),
    Location = trimws(Location)
  )


pcover <- pcover %>%
  mutate(Site = str_pad(Site, width = 2, pad = "0"),
         Plot = str_pad(Plot, width = 2, pad = "0"),
         SitePlotNew = paste(Location, Site, Plot, sep = ""))

#bruvs data is just baseline from 2021 so filtering pcover
pcover <- pcover %>% 
  filter(Year == 2021)

# Now join pcover to each of your datasets by SitePlotNew and Location

# For Snapper Abundance
snapper_abundance_pcover <- snapper_abundance_r %>%
  left_join(pcover, by = c("SitePlotNew", "Location")) %>%
  mutate(Category = as.factor(Category),
         Location = as.factor(Location))

# For Snapper Occurrence
snapper_occurrence_pcover <- snapper_occurrence_r %>%
  left_join(pcover, by = c("SitePlotNew", "Location")) %>%
  mutate(Category = as.factor(Category),
         Location = as.factor(Location))

# For Grouper Abundance
grouper_abundance_pcover <- grouper_abundance_r %>%
  left_join(pcover, by = c("SitePlotNew", "Location")) %>%
  mutate(Category = as.factor(Category),
         Location = as.factor(Location))

# For Grouper Occurrence
grouper_occurrence_pcover <- grouper_occurrence_r %>%
  left_join(pcover, by = c("SitePlotNew", "Location")) %>%
  mutate(Category = as.factor(Category),
         Location = as.factor(Location))

#to do: run glm?neg binomial, compare model structure NB vs Poisson, check_model, pairwise test to identify unique mean groups (gg predict), ggplot, then ggarrange.
#~ rugosity * location + coral + algae + sponge

#data |> filter(species == 'species name')

##(MAXN ~ LOCATION * AVG_RUGO * CATEGORY * PCOVER)
#SNAPPER - ABUNDANCE
snapper_model2 <- glm(MaxN ~ Location * Avg_rugosity * Category * Pcover, 
                      family = poisson, 
                      data = snapper_abundance_pcover)

# Check the model summary and ANOVA table
summary(snapper_model2)
anova(snapper_model2, test = "Chisq")
plot(snapper_model2)

#Location and Avg_rugosity are the strongest predictors of MaxN.
#Category and Pcover are not significant, but their interactions (especially with each other and with Location and Avg_rugosity) are significant.
#The significant higher-order interactions indicate that the relationship between rugosity and snapper abundance depends on both the type of habitat (Category) and its percent cover (Pcover), and that these relationships vary by location.

# Get predictions from the model.
# Here we let Avg_rugosity vary across its full range, and set Pcover to representative values (0, 1)
snapper_preds <- ggpredict(snapper_model2, 
                           terms = c("Avg_rugosity [all]", "Pcover [0, 1]", "Location", "Category [all]"))
head(snapper_preds)
colnames(snapper_preds)

unique(snapper_preds$facet)
#HELP: how do I model predict and plot?

#GROUPER - ABUNDANCE
grouper_model2 <- glm(MaxN ~ Location * Avg_rugosity * Category * Pcover, 
                      family = poisson, 
                      data = grouper_abundance_pcover)
summary(grouper_model2)
anova(grouper_model2, test = "Chisq")
plot(grouper_model2)
#Location:avg_rugosity significant. Avg_rugosity significant. Loocation significant. Pcover and category NOT.

grouper_preds <- ggpredict(grouper_model2, 
                           terms = c("Avg_rugosity [all]", "Pcover [0, 1]", "Location", "Category [all]"))
head(grouper_preds)

#help!
# Run pairwise comparisons using emmeans

# Get compact letter display (CLD) for grouper

# Merge CLD with predictions


# Create graph

##(OCC ~ LOCATION * AVG_RUGO * CATEGORY * PCOVER)
#SNAPPER - OCCURRENCE
snapper_model3 <- glm(PA ~ Location * Avg_rugosity * Category * Pcover, 
                      family = poisson, 
                      data = snapper_occurrence_pcover)

# Check the model summary and ANOVA table
summary(snapper_model3)
anova(snapper_model3, test = "Chisq") #avg_rugo and location:avg_rugo significant
plot(snapper_model3)

# Get predictions from the model.
# Here we let Avg_rugosity vary across its full range, and set Pcover to representative values (0, 1)
snapper_preds3 <- ggpredict(snapper_model3, 
                           terms = c("Avg_rugosity [all]", "Pcover [0, 1]", "Location", "Category [all]"))
head(snapper_preds3)

#how do I plot?

#GROUPER - OCCURRENCE
grouper_model3 <- glm(PA ~ Location * Avg_rugosity * Category * Pcover, 
                      family = poisson, 
                      data = grouper_occurrence_pcover)
summary(grouper_model3)
anova(grouper_model3, test = "Chisq") #location and avg_rugosity significant
plot(grouper_model3)

grouper_preds3 <- ggpredict(grouper_model2, 
                           terms = c("Avg_rugosity [all]", "Pcover [0, 1]", "Location", "Category [all]"))
head(grouper_preds3)

# Run pairwise comparisons using emmeans

# Get compact letter display (CLD) for grouper

# Merge CLD with predictions


# Create graph
