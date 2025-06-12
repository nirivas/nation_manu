#' """ Data analyses for BRUVs dataset: Grouper vs Snapper paper
#'     @authors: Christine Nation, Nicolas Rivas, W. Ryan James, Rolando O. Santos
#'     date: 03/27/24"""


# Load Libraries ----------------------------------------------------------

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


# Load Data ----------------------------------------------------------------

snapper = read_csv("data/snapper.csv") |> 
  mutate(Common = as.factor(Common))
grouper = read_csv("data/grouper.csv") # not including Nassau grouper
gs = read_csv("data/groupersnapper_pivot.csv") #For occurrence analysis

### For Rugosity ###

#Only top grouper and snapper

snapper_abundance_r = read_csv("data/snapper_abundance_r.csv")
grouper_abundance_r = read_csv("data/grouper_abundance_r.csv")
snapper_occurrence_r = read_csv("data/snapper_occurrence_r.csv")
grouper_occurrence_r = read_csv("data/grouper_occurrence_r.csv")

### For P cover ###

pcover = read_csv("data/old/pcover.csv")


# Q1: Species Differences -------------------------------------------------

## Abundance - MaxN ####

## snappers ####

#Comparing model structure = NB vs Poisson
glm_snap_MaxN.nb <- glmmTMB(MaxN ~ Common, data = snapper, family = nbinom2)

glmm_nb <- glmmTMB(MaxN ~ Common + (1 | Site), data = snapper, family = nbinom2)

glm_snap_MaxN.poi <- glm(MaxN ~ Common, family = poisson(link = "log"), data = snapper)

compare_performance(glm_snap_MaxN.nb,glm_snap_MaxN.poi, glmm_nb)

check_model(glmm_nb) # Model diagnostics
check_model(glm_snap_MaxN.nb) # Model diagnostics
#nb model better than poisson
rm(glm_snap_MaxN.poi,glm_snap_MaxN.nb)


#Checking for the model assumptions and results
x11()
check_model(glmm_nb) #looks ok

plot(glmm_nb, which = 3)
abline(h = 0, lty = 2)

summary(glmm_nb)
car::Anova(glmm_nb, test = "Chisq")
hist(snapper$MaxN)
ggplot(snapper,aes(Common,MaxN))+geom_boxplot()
#plot(glm_snap_MaxN.nb)

#Pairwise tests to identify unique mean groups
gg.snapper.MaxN = ggpredict(glmm_nb, terms = c("Common")) |> 
  rename(Common = x, MaxN = predicted)


em_snapper_MaxN = emmeans(glmm_nb, "Common")
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
snapper_MaxN

##grouper####

grouper$Common = as.factor(grouper$Common)

glm_grp_MaxN.nb <- glmmTMB(MaxN ~ Common, 
                           family = nbinom2(link = "log"), 
                           data = grouper)

glm_grp_MaxN.poi <- glmmTMB(MaxN ~ Common , 
                        family = poisson(link = "log"), 
                        data = grouper)

glmm_grp_MaxN_nb <- glmmTMB(MaxN ~ Common + (1 | Site), 
                             family = nbinom2(), 
                             data = grouper)

glmm_grp_MaxN_poi <- glmmTMB(MaxN ~ Common + (1 | Site), 
                              family = poisson(link = "log"), 
                              data = grouper)

compare_performance(glm_grp_MaxN.poi,glmm_snap_MaxN_poi, 
                    glm_grp_MaxN.nb, glmm_grp_MaxN_nb)

check_model(glmm_snap_MaxN_nb) # Model diagnostics


#Poisson is better than nb for grouper
rm(glm_grp_MaxN.nb,glmm_snap_MaxN_poi,glmm_snap_MaxN_nb)
x11()
check_model(glm_grp_MaxN.poi) #looks ok
summary(glm_grp_MaxN.poi)
anova(glm_grp_MaxN.poi, test = "Chisq")
#plot(glm_grp_MaxN.poi)

gg.grouper.MaxN = ggpredict(glm_grp_MaxN.poi, terms = c("Common"))  %>% 
  rename(Common = x, MaxN = predicted)


em_grouper_MaxN = emmeans(glm_grp_MaxN.poi, "Common")
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

grouper_MaxN

ggarrange(snapper_MaxN, grouper_MaxN,
          labels = c('a)','b)'),
          ncol = 2, vjust = 1, align = "h")

#ggsave("./figs/MaxN.Fit.GrouperSnapper.png", units = "in", width = 14, height = 8, dpi =  600)


# Q2: Occurrence  ---------------------------------------------------------

## snapper ####

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

# Create the Snapper Occurrence Plot
snapper_occ_plot <- ggplot(gg.snapper.occ2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Snapper Species", y = "Predicted Probability /n of Occurrence") +
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

snapper_occ_plot

## grouper ####

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

# Create the Grouper Occurrence Plot
grouper_occ_plot <- ggplot(gg.grouper.occ2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Grouper Species", y = "Predicted Probability /n of Occurrence") +
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

grouper_occ_plot

# Combine Snapper and Grouper Occurrence Plots
ggarrange(snapper_occ_plot, grouper_occ_plot,
          labels = c('a)', 'b)'),
          ncol = 1, vjust = 1, align = "v")



#ggsave("./figs/Occurrence.Fit.GrouperSnapper.png",units = "in", width = 8, height = 12, dpi = 600)

##combo plot ####

snapgroupoccur = rbind(gg.snapper.occ2, gg.grouper.occ2) |> 
  mutate(code = case_when(Common == "Graysby" ~ "d", Common == "Red Hind" ~ "e", Common == "Coney" ~ "f", TRUE ~ code)) |> 
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



##Occurrence by Loc ####

### For Snapper ####

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

check_model(glmmtmb_snap_occur_loc) # Model diagnostics
#

# Generate model predictions for each combination of Common and Location
gg.snapper.occ_loc <- ggpredict(glm_snap_occur_loc, terms = c("Common", "Location")) |>
  rename(Common = x, Occurrence = predicted, Location = group)


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
  labs(x = "Snapper Species", y = "Predicted Probability /n of Occurrence") +
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

### For Grouper ####

# Fit a GLM including Location as an interaction
glm_grouper_occur_loc <- glm(PA ~ Common * Location, 
                             family = binomial(link = "logit"), 
                             data = grouper_oc)

check_model(glm_grouper_occur_loc) # Model diagnostics
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


gg.grouper.occ_loc2 = gg.grouper.occ_loc2 |> 
  mutate(Location = case_when(Location == "M" ~ "Maguey",Location == "T" ~ "Tampico", TRUE ~ Location)) |> 
  drop_na()

# Create the Grouper Occurrence Plot by Location
grouper_occ_loc_plot <- ggplot(gg.grouper.occ_loc2, aes(Common, Occurrence, fill = Common), colour = "black") +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, color = "black") +
  geom_text(aes(label = code, y = conf.high + 0.05), size = 5) +
  labs(x = "Grouper Species", y = "Predicted Probability /n of Occurrence") +
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

#ggsave("./figs/Occurrence_byLocation.Fit.GrouperSnapper.png", units = "in", width = 14, height = 12, dpi = 600)


# Interaction Tests -------------------------------------------------------

# Snapper abundance with interaction (MaxN ~ Common * Site)
glm_snap_MaxN.nb_interaction <- glm.nb(MaxN ~ Common * Site, data = snapper)
summary(glm_snap_MaxN.nb_interaction)

# Grouper abundance with interaction (MaxN ~ Common * Site)
glm_grp_MaxN.nb_interaction <- glm.nb(MaxN ~ Common * Site, data = grouper)
summary(glm_grp_MaxN.nb_interaction)


# Rugosity Analysis -------------------------------------------------------

## Snapper Abundance ####

snapp_maxn_pois <- glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = poisson, data = snapper_abundance_r)
snapp_maxn_nb <- glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = nbinom2,data = snapper_abundance_r)

check_model(snapp_maxn_pois)
check_model(snapp_maxn_nb)
compare_performance(snapp_maxn_pois, snapp_maxn_nb) 
AIC(snapp_maxn_pois, snapp_maxn_nb) #NB lower AIC

#nb is better model
rm(snapp_maxn_pois)

check_overdispersion(snapp_maxn_nb)  # no overdispersion

summary(snapp_maxn_nb)     # Model coefficients
car::Anova(snapp_maxn_nb)  # ANOVA for significance
performance::check_collinearity(snapp_maxn_nb) #low correlations

options(na.action = "na.fail")
dredge(snapp_maxn_nb) # Model selection, looks like location is unnecessaty


snapper_abundance_r$Common = as.factor(snapper_abundance_r$Common)

snapp_maxn_nb <- glmmTMB(MaxN ~ Avg_rugosity  + Common  + (1 | Site), family = nbinom2,data = snapper_abundance_r)

check_model(snapp_maxn_nb)
car::Anova(snapp_maxn_nb) 
performance::check_collinearity(snapp_maxn_nb)


#I did also run an interaction between Avg_rugosity and common but it was not signifcant and did not inform on anything more


### Snapper  prediction ####


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

snappcommon_pred = ggplot(gg.snapp_maxn, aes(x = Common, y = MaxN, fill = Common)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = code), vjust = -0.5, size = 3) +
  labs(title = "Snapper Abundance (MaxN) by Species",
       x = "Species",
       y = "MaxN") +
  scale_fill_manual(values = custom_colors1[c(1, 4)])+
  theme_classic() +
  theme(legend.position = "none")

snappcommon_pred

snapp_maxn_rugopreds = ggpredict(snapp_maxn_nb, terms = c("Avg_rugosity [1:2.5 by = 0.05]"))

snaprugo_pred = ggplot(snapp_maxn_rugopreds, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D62", alpha = 0.3) +
  geom_line( linewidth = 1) +  
  theme_classic() +
  labs(title = "Predicted Snapper Abundance (MaxN) by Rugosity",
       x = "Rugosity",
       y = "MaxN")

snaprugo_pred

ggarrange(snappcommon_pred, snaprugo_pred,
          labels = c('a)', 'b)'),
          ncol = 2, vjust = 1, align = "v")

#ggsave("./figs/Snapper_Abundance.Rugosity.png", units = "in", width = 12, height = 8, dpi = 600)

## Grouper Abundance ####

grouper_maxn_pois1 = glmmTMB(MaxN ~ Avg_rugosity * Location * Common  * (1 | Site), family = poisson, data = grouper_abundance_r)
grouper_maxn_nb1 = glmmTMB(MaxN ~ Avg_rugosity * Location * Common  * (1 | Site), family = nbinom2,data = grouper_abundance_r)
grouper_maxn_pois2 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = poisson, data = grouper_abundance_r)
grouper_maxn_nb2 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = nbinom2,data = grouper_abundance_r)

compare_performance(grouper_maxn_pois1, grouper_maxn_nb1,grouper_maxn_pois2,grouper_maxn_nb2) #grouper_maxn_pois2 lower AIC, better fit
rm(grouper_maxn_pois1,grouper_maxn_nb1,grouper_maxn_pois2,grouper_maxn_nb2)

grouper_maxn_pois1 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = poisson, data = grouper_abundance_r)

grouper_maxn_pois2 = glmmTMB(MaxN ~ Avg_rugosity + Location + Common  + (1 | Site), family = compois(), data = grouper_abundance_r)

compare_performance(grouper_maxn_pois1,grouper_maxn_pois2)
check_model(grouper_maxn_pois2) #looks ok, misspecified doesnt look great
check_overdispersion(grouper_maxn_pois)  # no overdispersion

summary(grouper_maxn_pois)     # Model coefficients
car::Anova(grouper_maxn_pois, test = "Chisq")  # ANOVA for significance
performance::check_collinearity(grouper_maxn_pois)  # Check multicollinearity. high correlation.
rm(grouper_maxn_pois1,grouper_maxn_pois2)

grouper_abundance_r$Common = as.factor(grouper_abundance_r$Common)

grouper_maxn_pois_r = glmmTMB(MaxN ~ Avg_rugosity  + Common  + (1 | Site), family = compois(), data = grouper_abundance_r)

dredge(grouper_maxn_pois) 

summary(grouper_maxn_pois_r)     # Model coefficients

### Grouper prediction####
grouper_maxn_preds_r <- ggpredict(grouper_maxn_pois_r, terms = c("Common")) %>% 
  rename(Common = x, MaxN = predicted)

em_grouper_maxn_r <- emmeans(grouper_maxn_pois_r, ~ Common)
pairs(em_grouper_maxn_r, simple = "Common")

cld_grouper_maxn_r <- cld(em_grouper_maxn_r, 
                        Letters = letters,
                        adjust = "bonf", 
                        sort = FALSE, 
                        alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))
gg.grouper_maxn_r = left_join(grouper_maxn_preds_r, cld_grouper_maxn_r, by = c("Common"))

groupercommon_pred = ggplot(gg.grouper_maxn_r, aes(x = Common, y = MaxN, fill = Common)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = code), vjust = -0.5, size = 3) +
  labs(title = "Grouper Abundance (MaxN) by Species",
       x = "Species",
       y = "MaxN") +
  scale_fill_manual(values = g_colors1[c(1, 2, 3)])+
  theme_classic() +
  theme(legend.position = "none")

groupercommon_pred

grouper_maxn_rugopreds = ggpredict(grouper_maxn_pois_r, terms = c("Avg_rugosity [1:2.5 by = 0.05]"))

grouper_rugo_pred = ggplot(grouper_maxn_rugopreds, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D62", alpha = 0.3) +
  geom_line( linewidth = 1) +  
  theme_classic() +
  labs(title = "Predicted Grouper Abundance (MaxN) by Rugosity",
       x = "Rugosity",
       y = "MaxN")
grouper_rugo_pred

ggarrange(groupercommon_pred, grouper_rugo_pred,
          labels = c('a)', 'b)'),
          ncol = 2, vjust = 1, align = "v")
#ggsave("./figs/Grouper_Abundance.Rugosity.png", units = "in", width = 12, height = 8, dpi = 600)

## Snapper Occuurence ####

snapp_occur_binom = glmmTMB(PA ~ Avg_rugosity  + Location + Common  + (1 | Site), family = binomial, data = snapper_occurrence_r)
check_model(snapp_occur_binom1) 
dredge(snapp_occur_binom1) 

snapp_occur_binom = glmmTMB(PA ~ Avg_rugosity + Common  + (1 | Site), family = binomial, data = snapper_occurrence_r)

check_model(snapp_occur_binom) 

snapper_occurrence_r$Common = as.factor(snapper_occurrence_r$Common)
snapp_occur_binom = glmmTMB(PA ~ Avg_rugosity  + Common  + (1 | Site), family = binomial, data = snapper_occurrence_r)

summary(snapp_occur_binom)

car::Anova(snapp_occur_binom)

snapp_paoccur_preds = ggpredict(snapp_occur_binom, terms = c("Common")) %>% 
  rename(Common = x, Occurrence = predicted)

em_snapp_paoccur = emmeans(snapp_occur_binom, ~ Common)
pairs(em_snapp_paoccur, simple = "Common")

cld_snapp_paoccur = cld(em_snapp_paoccur, 
                      Letters = letters,
                      adjust = "bonf", 
                      sort = FALSE, 
                      alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

gg.snapp_paoccur = left_join(snapp_paoccur_preds, cld_snapp_paoccur, by = c("Common"))

snappcommon_paoccur = ggplot(gg.snapp_paoccur, aes(x = Common, y = Occurrence, fill = Common)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = code), vjust = -0.5, size = 3) +
  labs(title = "Snapper Occurrence by Species",
       x = "Species",
       y = "Occurrence") +
  scale_fill_manual(values = custom_colors1[c(1, 4)])+
  theme_classic() +
  theme(legend.position = "none")

snappcommon_paoccur

snapp_parugopreds = ggpredict(snapp_occur_binom, terms = c("Avg_rugosity [1:2.5 by = 0.05]"))

snap_pa_rugo_pred = ggplot(snapp_parugopreds, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D62", alpha = 0.3) +
  geom_line( linewidth = 1) +  
  theme_classic() +
  labs(title = "Predicted Snapper Occurrence by Rugosity",
       x = "Rugosity",
       y = "Occurrence")

snap_pa_rugo_pred


ggarrange(snappcommon_paoccur, snap_pa_rugo_pred,
          labels = c('a)', 'b)'),
          ncol = 2, vjust = 1, align = "v")

#ggsave("./figs/Snapper_Occurrence.Rugosity.png", units = "in", width = 12, height = 8, dpi = 600)

## Grouper Occurrence ####

grouper_occur_binom1 = glmmTMB(PA ~ Avg_rugosity  + Location + Common  + (1 | Site), family = binomial, data = grouper_occurrence_r)
check_model(grouper_occur_binom1)
dredge(grouper_occur_binom1)

grouper_occur_binom2 = glmmTMB(PA ~ Avg_rugosity + Common  + (1 | Site), family = binomial, data = grouper_occurrence_r)
check_model(grouper_occur_binom2)
compare_performance(grouper_occur_binom1, grouper_occur_binom2) # 2 has a higher weight, will continue with this model

rm(grouper_occur_binom1, grouper_occur_binom2)

grouper_occurrence_r$Common = as.factor(grouper_occurrence_r$Common)
grouper_occur_binom = glmmTMB(PA ~ Avg_rugosity + Common  + (1 | Site), family = binomial, data = grouper_occurrence_r)
summary(grouper_occur_binom)
car::Anova(grouper_occur_binom)

grouper_paoccur_preds = ggpredict(grouper_occur_binom, terms = c("Common")) %>% 
  rename(Common = x, Occurrence = predicted)

em_grouper_paoccur = emmeans(grouper_occur_binom, ~ Common)
pairs(em_grouper_paoccur, simple = "Common")

cld_grouper_paoccur = cld(em_grouper_paoccur, 
                      Letters = letters,
                      adjust = "bonf", 
                      sort = FALSE, 
                      alpha = 0.05) |>
  mutate(code = str_replace_all(.group, ' ', ''))

gg.grouper_paoccur = left_join(grouper_paoccur_preds, cld_grouper_paoccur, by = c("Common"))

groupercommon_paoccur = ggplot(gg.grouper_paoccur, aes(x = Common, y = Occurrence, fill = Common)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = code), vjust = -0.5, size = 3) +
  labs(title = "Grouper Occurrence by Species",
       x = "Species",
       y = "Occurrence") +
  scale_fill_manual(values = g_colors1[c(1, 2, 3)])+
  theme_classic() +
  theme(legend.position = "none")

groupercommon_paoccur

grouper_parugopreds = ggpredict(grouper_occur_binom, terms = c("Avg_rugosity [1:2.5 by = 0.05]"))

grouper_pa_rugo_pred = ggplot(grouper_parugopreds, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FC8D62", alpha = 0.3) +
  geom_line( linewidth = 1) +  
  theme_classic() +
  labs(title = "Predicted Grouper Occurrence by Rugosity",
       x = "Rugosity",
       y = "Occurrence")

grouper_pa_rugo_pred

grouper_parugopreds = ggpredict(grouper_occur_binom, terms = c("Avg_rugosity [1:2.5 by = 0.05]")) |> 
  rename(Rugosity = x, Occurrence = predicted)


ggarrange(groupercommon_paoccur, grouper_pa_rugo_pred,
          labels = c('a)', 'b)'),
          ncol = 2, vjust = 1, align = "v")


test = as.data.frame(gg.grouper_paoccur)
test
test2 = as.data.frame(grouper_parugopreds)
test2

test3 = ggpredict(grouper_occur_binom, terms = c("Avg_rugosity [1:2.5 by = 0.05]", "Common")) |>
  rename(Rugosity = x, Occurrence = predicted, Common = group)

ggplot(test3, aes(x = Rugosity, y = Occurrence, color = Common, fill = Common)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +  
  theme_classic() +
  labs(title = "Predicted Grouper Occurrence by Rugosity",
       x = "Rugosity",
       y = "Probability of Occurrence")



# Pcover ------------------------------------------------------------------

grouper_occurrence_r = read_csv("data/grouper_occurrence_r.csv") |> 
  dplyr::select(-SitePlot) |> 
  rename(SitePlot = SitePlotNew) 

snapper_occurrence_r = read_csv("data/snapper_occurrence_r.csv") |> 
  dplyr::select(-SitePlot) |> 
  rename(SitePlot = SitePlotNew)

pcover = read_csv("data/old/pcover.csv") |> 
  mutate(Pcover = Pcover*100,
         Plot = str_pad(Plot, width = 2, side = "left", pad = "0"),
         Site = str_pad(Site, width = 2, side = "left", pad = "0")) |>
  unite("SitePlot", Location, Site, Plot, sep = "", remove = FALSE) |> 
  filter(Category %in% c("coral", "gorgonian", "algae", "sponge","sand-rubble"),
         Period == "baseline") |> 
  dplyr::select(-Period, -Year,-Location,-Site,-Plot) 


lsmetrics = read_csv("data/lsmetrics.csv")


dft = left_join(grouper_occurrence_r, pcover, by = c("SitePlot")) 

dft_wide = dft |> 
  pivot_wider(
    names_from = Category,
    values_from = Pcover
  ) |> 
  dplyr::select(-14) |> 
  drop_na() |> 
  rename(
    sand = "sand-rubble") |> 
  left_join(lsmetrics, by = c("SitePlot"))

df_gray_O = dft_wide |> 
  filter(Common == "Graysby") 


grayO1 = glmmTMB(PA ~ Avg_rugosity + Location +  coral + gorgonian + sand + sponge + algae  + (1 | Site), family = binomial, data = df_gray_O)


check_model(grayO1)
options(na.action = "na.fail")
dredge(grayO1) # Model selecection. Going with most pasimonius that only inclues Avg_rugosity 

grayo1 = glmmTMB(PA ~ Avg_rugosity + (1 | Site), family = binomial, data = df_gray_O)


check_model(grayo1) # Model diagnostics

grayo2_1 = glmmTMB(PA ~ pd_AggregateReef + pd_AggregatedPatchReefs + pd_SandwithScatteredCoralandRock + pd_Seagrass + 
                   pland_AggregateReef + pland_AggregatedPatchReefs + pland_SandwithScatteredCoralandRock + pland_Seagrass
                   + (1|Site), data = df_gray_O, family = binomial)


grayo2_2 = glmmTMB(PA ~ pd_AggregateReef + pd_AggregatedPatchReefs + pd_SandwithScatteredCoralandRock + pd_Seagrass + 
                   pland_AggregateReef + pland_AggregatedPatchReefs + pland_SandwithScatteredCoralandRock + pland_Seagrass, 
                 data = df_gray_O, family = binomial)

compare_performance(grayo2_1, grayo2_2) # Model comparison, grayo2_2 is better
rm(grayo2_1)

check_model(grayo2_2) # Model diagnostics, looks okayish

dr_invT = dredge(grayo2_2) |> 
  filter(delta < 4)

grayo2_top = which.min(dr_invT$df)

grayo2_tops = get.models(dr_invT, subset = top)[[1]]

summary(grayo2_tops)

grayo2 = glmmTMB(PA ~ pland_AggregateReef + pland_SandwithScatteredCoralandRock + pland_Seagrass, data = df_gray_O, family = binomial)
check_model(grayo2) # Model diagnostics
summary(grayo2)

grayo3_1 = glmmTMB(PA ~ mpa_distance + mang + land + inlet + (1 | Site), data = df_gray_O, family = binomial)
grayo3_2 = glmmTMB(PA ~ mpa_distance + mang + land + inlet, data = df_gray_O, family = binomial)
compare_performance(grayo3_1, grayo3_2) # Model comparison, grayo3_1 is better
rm(grayo3_1)
check_model(grayo3_2) # Model diagnostics

dr_invT = dredge(grayo3_2) |> 
  filter(delta < 4)

grayo3_top = which.min(dr_invT$df)

grayo3_tops = get.models(dr_invT, subset = top)[[1]]

summary(grayo3_tops)

grayo3 = glmmTMB(PA ~  mang +  inlet, data = df_gray_O, family = binomial)

summary(grayo3) # Model coefficients

grayofinal_1 = glmmTMB(PA ~ pland_AggregateReef + pland_SandwithScatteredCoralandRock + pland_Seagrass + mang + inlet
                     + (1|Site), data = df_gray_O, family = binomial)
grayofinal_2 = glmmTMB(PA ~ pland_AggregateReef + pland_SandwithScatteredCoralandRock + pland_Seagrass + mang + inlet,
                      data = df_gray_O, family = binomial)
compare_performance(grayofinal_1, grayofinal_2) # Model comparison, grayofinal_2 is better

check_model(grayofinal_2) # Model diagnostics, looks okayish
summary(grayofinal_2) # Model coefficients


# Abundance ---------------------------------------------------------------


grouper_abundance_r = read_csv("data/grouper_abundance_r.csv") |> 
  dplyr::select(-SitePlot) |> 
  rename(SitePlot = SitePlotNew) 


pcover = read_csv("data/old/pcover.csv") |> 
  mutate(Pcover = Pcover*100,
         Plot = str_pad(Plot, width = 2, side = "left", pad = "0"),
         Site = str_pad(Site, width = 2, side = "left", pad = "0")) |>
  unite("SitePlot", Location, Site, Plot, sep = "", remove = FALSE) |> 
  filter(Category %in% c("coral", "gorgonian", "algae", "sponge","sand-rubble"),
         Period == "baseline") |> 
  dplyr::select(-Period, -Year,-Location,-Site,-Plot) 

lsmetrics = read_csv("data/lsmetrics.csv")


dft = left_join(grouper_abundance_r, pcover, by = c("SitePlot")) 

dft_wide = dft |> 
  pivot_wider(
    names_from = Category,
    values_from = Pcover
  ) |> 
  dplyr::select(-14) |> 
  drop_na() |> 
  rename(
    sand = "sand-rubble") |> 
  left_join(lsmetrics, by = c("SitePlot"))

df_gray_O = dft_wide |> 
  filter(Common == "Graysby") 

grayA1 = glmmTMB(MaxN ~ Avg_rugosity + Location +  coral + gorgonian + sand + sponge + algae  + (1 | Site), family = poisson(), data = df_gray_O)
check_model(grayA1)
dredge(grayA1)

dr_invT = dredge(grayA1) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)

grayA1 = glmmTMB(MaxN ~ sand + (1 | Site), family = poisson(), data = df_gray_O) #only good one

grayA2 = glmmTMB(MaxN ~ pland_AggregateReef + pland_AggregatedPatchReefs + pland_SandwithScatteredCoralandRock + pland_Seagrass
                   + (1|Site), data = df_gray_O, family = poisson())

check_model(grayA2) # Model diagnostics, looks okayish
dredge(grayA2) 

dr_invT = dredge(grayA2) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)

grayA3 = glmmTMB(MaxN ~ mpa_distance + mang + land + inlet, data = df_gray_O, family = nbinom2())
check_model(grayA3) # Model diagnostics, looks okayish
dredge(grayA3)


dr_invT = dredge(grayA3) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)

grayAF1 = glmmTMB(MaxN ~ sand + mang , data = df_gray_O, family = nbinom2(link = "sqrt"))
grayAF2 = glmmTMB(MaxN ~ sand + mang  + (1|Site), data = df_gray_O, family = nbinom2(link = "sqrt"))
compare_performance(grayAF1, grayAF2) # Model comparison, grayAF2 is better
AIC(grayAF1, grayAF2) # AIC comparison, grayAF2 is better

check_model(grayAF2) # Model diagnostics, looks okayish
summary(grayAF2) # Model coefficients
