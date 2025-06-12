#' """ Data analyses for BRUVs dataset: Grouper vs Snapper paper
#'     @authors: Christine Nation, Nicolas Rivas, W. Ryan James, Rolando O. Santos
#'     date: 03/27/24"""

library(tidyverse)
library(viridis)
library(ggpubr)

#Models libraries
library(mgcv)
library(glmmTMB)

# library(gam)
# library(gamm4)
# library(nlme)

#Model check libraries
library(visreg)
library(DHARMa)
library(performance)
library(lsmeans)
library(emmeans)
library(multcomp)
library(multcompView)
library(MuMIn)
library(ggeffects)

# Load required packages
library(MASS)

###############################################################################


##########################Data processing#######################################

#new_gs_data is the bruvs masterlist already filtered with just groupers and snappers (gs)
new_gs_data<- read_csv("data/old/NEW_GS_DATA.csv") %>%   unite(SitePlot, Site, Plot, sep = "", remove = FALSE) |> 
  mutate(Common = case_when(Common == "Red_hind" ~ "Red Hind", TRUE ~ Common))

list(unique(new_gs_data$Common))
unique(new_gs_data$Site)
length(unique(new_gs_data)) 
length(unique(new_gs_data$SitePlot))
       
#Change MaxN to numeric 
new_gs_data$MaxN<-as.numeric(new_gs_data$MaxN)

# Replace "MAG" with "M" and "TAM" with "T" in the Site column (to match the rugosity viscore data)
new_gs_data$Site <- gsub("^MAG", "M", new_gs_data$Site)
new_gs_data$Site <- gsub("^TAM", "T", new_gs_data$Site)

# Add leading zeros to Plot in new_gs2
new_gs_data$Plot <- sprintf("%02d", as.numeric(new_gs_data$Plot))

#Cut videos to 60 minutes and flipped Time before flipped TBF
new_gs2 <- new_gs_data %>%
  group_by(Site, Plot, Location, Common, Taxa) %>%
  filter(DifT <= 60, TBF >= 60) %>%
  summarize(MaxN = max(MaxN)) %>%
  unite(SitePlot, Site, Plot, sep = "", remove = FALSE) %>%
  distinct() %>%
  ungroup()

#Filter by just snapper data
snapper = new_gs2 |> 
  filter(Taxa %in% c("Snapper")) 
#Checking frequency of observation to eliminate species with <5% observations
snapper |> count(Common) #no elimination based on the 5% threshold

#Filter by just grouper data
grouper = new_gs2 |> 
  filter(Taxa %in% c("Grouper")) 
#Checking frequency of observation to eliminate species with <5% observations
grouper |> count(Common) #Nassau eliminated for analysis purposes 

grouper = grouper |> 
  filter(Common != c("Nassau"))

#write_csv(snapper, "data/snapper.csv")
#write_csv(grouper, "data/grouper.csv")

#################################################################################



# Question 1 Analysis - Species differences -------------------------------
#Model 1 = Identify difference between species to help id species for submodels

###
#Abundance - MaxN
##
library(MASS)

#snappers
# (Rolo/Ryan did this part) task: run glm(MaxN ~ Species, family = negative binomial) 

#Comparing model structure = NB vs Poisson
glm_snap_MaxN.nb <- glm.nb(MaxN ~ Common, data = snapper)
glm_snap_MaxN.poi <- glm(MaxN ~ Common, family = poisson(link = "log"), data = snapper)
compare_performance(glm_snap_MaxN.nb, glm_snap_MaxN.poi)
#nb model better than poisson

#Checking for the model assumptions and results
check_model(glm_snap_MaxN.nb) #looks ok
summary(glm_snap_MaxN.nb)
anova(glm_snap_MaxN.nb, test = "Chisq")
hist(snapper$MaxN)
ggplot(snapper,aes(Common,MaxN))+geom_boxplot()+scale_y_log10()
#plot(glm_snap_MaxN.nb)

#Pairwise tests to identify unique mean groups
gg.snapper.MaxN = ggpredict(glm_snap_MaxN.nb, terms = c("Common")) |> 
  rename(Common = x, MaxN = predicted)


em_snapper_MaxN = emmeans(glm_snap_MaxN.nb, "Common")
pairs(em_snapper_MaxN, simple = "Common") |> capture.output(file = "./tables/q1.emmean.pairs.SnapperSummary.txt")

#Add letter code to unique mean groups
meanmodel_means_cld = cld(object = em_snapper_MaxN,
                          adjust = "bonf",
                          Letters = letters,
                          sort = F,
                          alpha = 0.05) |> 
  mutate(code = str_replace_all(.group, ' ', ''))

#Merging fitted models with unique groups code
custom_colors1 <- c("#FEE08B","#E69F00","#009E73", "#56B4E9","#0072B2")  # color blind friendly

gg.snapper.MaxN2 = left_join(gg.snapper.MaxN, meanmodel_means_cld, by = "Common")

snapper_MaxN = ggplot(gg.snapper.MaxN2, aes(Common, MaxN, fill = Common), colour = "black")+
  geom_col()+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black")+
  geom_text(aes(label = code, y = conf.high + 0.5), size = 5)+
  labs(x = "Snapper Species", y = "MaxN")+
  scale_fill_manual(values = custom_colors1) +  # Use the custom color palette
  # scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

##grouper
grouper$Common = as.factor(grouper$Common)
glm_grp_MaxN.nb <- glm.nb(MaxN ~ Common, data = grouper)
glm_grp_MaxN.poi <- glm(MaxN ~ Common, family = poisson(link = "log"), data = grouper)
compare_performance(glm_grp_MaxN.nb, glm_snap_MaxN.poi)
#nb model better than poisson

check_model(glm_grp_MaxN.nb) #looks ok
summary(glm_grp_MaxN.nb)
anova(glm_grp_MaxN.nb, test = "Chisq")
#plot(glm_grp_MaxN.nb)

gg.grouper.MaxN = ggpredict(glm_grp_MaxN.nb, terms = c("Common"))  %>% 
  rename(Common = x, MaxN = predicted)


em_grouper_MaxN = emmeans(glm_grp_MaxN.nb, "Common")
pairs(em_grouper_MaxN, simple = "Common") |> capture.output(file = "./tables/q1.emmean.pairs.GrouperSummary.txt")

meanmodel_means_cld2 = cld(object = em_grouper_MaxN,
                           adjust = "bonf",
                           Letters = letters,
                           sort = F,
                           alpha = 0.05) |> 
  mutate(code = str_replace_all(.group, ' ', ''))


gg.grouper.MaxN2 = left_join(gg.grouper.MaxN, meanmodel_means_cld2, by = "Common")

g_colors1 <- c("#5AB4AC","#F46D43","#B5A1D3")

grouper_MaxN = ggplot(gg.grouper.MaxN2, aes(Common, MaxN, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.2), size = 5) +
  labs(x = "Grouper Species", y = "MaxN") +
  scale_fill_manual(values = g_colors1) +  # Use the custom color palette
 # scale_x_discrete(labels = c("Graysby", "Red Hind", "Coney")) +
  theme(
    axis.text = element_text(size = 14, face = "bold", colour = "black"),
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )


ggarrange(snapper_MaxN, grouper_MaxN,
          labels = c('a)','b)'),
          ncol = 2, vjust = 1, align = "h")


#ggsave("./figs/MaxN.Fit.GrouperSnapper.png", units = "in", width = 14, height = 8, dpi =  600)


####
####
###
#Occurrence (Christine did this- needs checking) task: run glm(occurrence ~ Species, family = binomial) 
###---------------------------------------------
# Part 2: Occurrence Analysis (Matching Part 1 Workflow)
#---------------------------------------------

#---------------------------------------------
# Data Preparation
#---------------------------------------------
# Pivot and calculate presence/absence (PA)
gs <- new_gs2 |>
  dplyr::select(-Taxa) |>
  pivot_wider(names_from = Common, values_from = MaxN, values_fill = 0) |>
  pivot_longer(cols = Graysby:Nassau, names_to = "Common", values_to = "MaxN") |>
  mutate(PA = if_else(MaxN > 0, 1, 0))

#write_csv(gs, "data/groupersnapper_pivot.csv")

#---------------------------------------------
# Section 1: Occurrence GLM Analysis
#---------------------------------------------

### Snapper Occurrence Analysis ###
snapper_oc <- gs |>
  filter(Common %in% c("Yellowtail", "Lane", "Mahogany", "Schoolmaster", "Mutton"))


# Fit GLM (binomial) for snapper occurrence
glm_snap_occur <- glm(PA ~ Common, family = binomial(link = "logit"), data = snapper_oc)
check_model(glm_snap_occur)   # Model diagnostics
summary(glm_snap_occur)
anova(glm_snap_occur, test = "Chisq")

# Obtain predicted occurrence probabilities for each species
gg.snapper.occ <- ggpredict(glm_snap_occur, terms = "Common") |>
  rename(Common = x, Occurrence = predicted)

# Pairwise comparisons using emmeans
em_snapper_occ <- emmeans(glm_snap_occur, "Common")
pairs(em_snapper_occ, simple = "Common") |>
  capture.output(file = "./tables/q2.emmean.pairs.SnapperOccurrenceSummary.txt")

# Add compact letter display (CLD) for significance groups
meanmodel_means_cld_occ <- cld(em_snapper_occ,
                               adjust = "bonf",
                               Letters = letters,
                               sort = FALSE,
                               alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

# Merge predicted values with CLD letters
gg.snapper.occ2 <- left_join(gg.snapper.occ, meanmodel_means_cld_occ, by = "Common")

# Define custom colors (as used in Part 1)
custom_colors1 <- c("#FEE08B","#E69F00","#009E73", "#56B4E9","#0072B2")

# Create the Snapper Occurrence Plot
snapper_occ_plot <- ggplot(gg.snapper.occ2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Snapper Species", y = "Predicted Probability \n of Occurrence") +
  scale_fill_manual(values = custom_colors1) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

### Grouper Occurrence Analysis ###
grouper_oc <- gs |>
  filter(Common %in% c("Graysby", "Red Hind", "Coney"))

# Fit GLM (binomial) for grouper occurrence
glm_grouper_occur <- glm(PA ~ Common, family = binomial(link = "logit"), data = grouper_oc)
check_model(glm_grouper_occur)
summary(glm_grouper_occur)
anova(glm_grouper_occur, test = "Chisq")

# Obtain predicted occurrence probabilities for grouper
gg.grouper.occ <- ggpredict(glm_grouper_occur, terms = "Common") |>
  rename(Common = x, Occurrence = predicted)

# Pairwise comparisons for grouper occurrence
em_grouper_occ <- emmeans(glm_grouper_occur, "Common")
pairs(em_grouper_occ, simple = "Common") |>
  capture.output(file = "./tables/q2.emmean.pairs.GrouperOccurrenceSummary.txt")

# Add CLD for grouper occurrence
meanmodel_means_cld_occ_grp <- cld(em_grouper_occ,
                                   adjust = "bonf",
                                   Letters = letters,
                                   sort = FALSE,
                                   alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

# Merge predicted values with CLD letters
gg.grouper.occ2 <- left_join(gg.grouper.occ, meanmodel_means_cld_occ_grp, by = "Common")

# Define grouper color palette (as in Part 1)
g_colors1 <- c("#5AB4AC","#F46D43","#B5A1D3")

# Create the Grouper Occurrence Plot
grouper_occ_plot <- ggplot(gg.grouper.occ2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Grouper Species", y = "Predicted Probability \n of Occurrence") +
  scale_fill_manual(values = g_colors1) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

snapgroupoccur = rbind(gg.snapper.occ2, gg.grouper.occ2) |> 
  mutate(code = case_when(Common == "Graysby" ~ "d", Common == "Coney" ~ "e", TRUE ~ code)) |> 
  drop_na()

ggplot(snapgroupoccur, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "red", linewidth = 1)+
  labs(x = "Species", y = "Predicted Occurrence Probability") +
  scale_fill_manual(values = c(custom_colors1, g_colors1)) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
#ggsave("./figs/Occurrence.Fit.GrouperSnapperNICO.png", units = "in", width = 12, height = 8, dpi = 600)

# Combine Snapper and Grouper Occurrence Plots
ggarrange(snapper_occ_plot, grouper_occ_plot,
          labels = c('a)', 'b)'),
          ncol = 1, vjust = 1, align = "v")



#ggsave("./figs/Occurrence.Fit.GrouperSnapper.png",units = "in", width = 8, height = 12, dpi = 600)

#---------------------------------------------
# Section 2 (Revised): Occurrence Analysis by Location – Model Based
#---------------------------------------------

# For Snapper:
# Fit a GLM including Location as an interaction
glm_snap_occur_loc <- glm(PA ~ Common * Location, 
                          family = binomial(link = "logit"), 
                          data = snapper_oc)
#plot(glm_snap_occur_loc)
summary(glm_snap_occur_loc)
#yellowtail high positive coefficient -> higher occurrence. problematic due to high std error
#The interaction for Schoolmaster is significant -> the effect of location on occurrence is different for Schoolmaster compared to the baseline species.
#Other species effects and the main effect of Location (Tampico vs. Maguey) are not statistically significant.

anova(glm_snap_occur_loc, test = "Chisq")
#

# Generate model predictions for each combination of Common and Location
gg.snapper.occ_loc <- ggpredict(glm_snap_occur_loc, terms = c("Common", "Location")) |>
  rename(Common = x, Occurrence = predicted, Location = group)

#gg.snapper.occ_loc <- ggpredict(glm_snap_occur_loc, terms = c("Common", "Location"))
#head(gg.snapper.occ_loc)

gg.snapper.occ <- ggpredict(glm_snap_occur, terms = "Common") |>
  rename(Common = x, Occurrence = predicted)

# Perform pairwise comparisons within each Location using emmeans
em_snap_occur_loc <- emmeans(glm_snap_occur_loc, ~ Common | Location)
pairs(em_snap_occur_loc, simple = "Common") |>
  capture.output(file = "./tables/q2.emmean.pairs.SnapperOccurrence_byLocation.txt")

# Get compact letter displays for significance groups
cld_snap_occur_loc <- cld(em_snap_occur_loc,
                          adjust = "bonf",
                          Letters = letters,
                          sort = FALSE,
                          alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

# Merge the CLD letters with the predicted values
gg.snapper.occ_loc2 <- left_join(gg.snapper.occ_loc, cld_snap_occur_loc, 
                                 by = c("Common", "Location"))

# Define the custom color palette (reusing Part 1 colors)
custom_colors1 <- c("#FEE08B","#E69F00","#009E73", "#56B4E9","#0072B2")

gg.snapper.occ_loc2 = gg.snapper.occ_loc2 |> 
  mutate(Location = case_when(Location == "M" ~ "Maguey",Location == "T" ~ "Tampico", TRUE ~ Location),
         code = case_when(Location == "Tampico" & Common == "Yellowtail" ~ "*", TRUE ~ code)) |> 
  drop_na()

# Create the Snapper Occurrence Plot by Location
snapper_occ_loc_plot <- ggplot(gg.snapper.occ_loc2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(data = subset(gg.snapper.occ_loc2, !(Common == "Yellowtail" & Location == "Tampico")),
                aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  scale_fill_manual(values = custom_colors1) +
  facet_wrap(~ Location) +
  labs(x = "Snapper Species", y = "Predicted Probability \n of Occurrence") +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold", colour = "black"))

snapper_occ_loc_plot

# For Grouper:
# Fit a GLM including Location as an interaction
glm_grouper_occur_loc <- glm(PA ~ Common * Location, 
                             family = binomial(link = "logit"), 
                             data = grouper_oc)
#plot(glm_grouper_occur_loc)
summary(glm_grouper_occur_loc)
anova(glm_grouper_occur_loc, test = "Chisq")

# Generate predictions
gg.grouper.occ_loc <- ggpredict(glm_grouper_occur_loc, terms = c("Common", "Location"))|>
  rename(Common = x, Occurrence = predicted, Location = group)

# Pairwise comparisons within each Location
em_grouper_occur_loc <- emmeans(glm_grouper_occur_loc, ~ Common | Location)
pairs(em_grouper_occur_loc, simple = "Common") |>
  capture.output(file = "./tables/q2.emmean.pairs.GrouperOccurrence_byLocation.txt")

# Obtain CLD letters for grouper
cld_grouper_occur_loc <- cld(em_grouper_occur_loc,
                             adjust = "bonf",
                             Letters = letters,
                             sort = FALSE,
                             alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

# Merge CLD letters with predictions
gg.grouper.occ_loc2 <- left_join(gg.grouper.occ_loc, cld_grouper_occur_loc, 
                                 by = c("Common", "Location"))

# Define grouper custom colors (as used in Part 1)
g_colors1 <- c("#5AB4AC","#F46D43","#B5A1D3")

gg.grouper.occ_loc2 = gg.grouper.occ_loc2 |> 
  mutate(Location = case_when(Location == "M" ~ "Maguey",Location == "T" ~ "Tampico", TRUE ~ Location)) |> 
  drop_na()

# Create the Grouper Occurrence Plot by Location
grouper_occ_loc_plot <- ggplot(gg.grouper.occ_loc2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Grouper Species", y = "Predicted Probability \n of Occurrence") +
  scale_fill_manual(values = g_colors1) +
  facet_wrap(~ Location) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold", colour = "black"))
grouper_occ_loc_plot

# Optionally, combine and save the plots
ggarrange(snapper_occ_loc_plot, grouper_occ_loc_plot,
          labels = c('a)', 'b)'),
          ncol = 1, vjust = 1, align = "v")

#ggsave("./figs/Occurrence_byLocation.Fit.GrouperSnapper.png", 
#        units = "in", width = 14, height = 12, dpi = 600)

#---------------------------------------------
# Section 3: Additional - Interaction Tests on Abundance
# (These replicate extra tests from your provided code)
#---------------------------------------------

# Snapper abundance with interaction (MaxN ~ Common * Site)
glm_snap_MaxN.nb_interaction <- glm.nb(MaxN ~ Common * Site, data = snapper)
summary(glm_snap_MaxN.nb_interaction)

# Grouper abundance with interaction (MaxN ~ Common * Site)
glm_grp_MaxN.nb_interaction <- glm.nb(MaxN ~ Common * Site, data = grouper)
summary(glm_grp_MaxN.nb_interaction)

