#' ---
#' title: "R code for the analysis conducted in 'BrainAGE and regional volumetric analysis of a Buddhist monk: a longitudinal MRI case study'"
#' author: "Nagesh Adluru and Derek L. Norton"
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#' ---
#'
#' # Initialization
# Loading the libraries ====
library(ggsci)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(latex2exp)
library(forcats)
library(broom)
library(purrr)
library(scales)
library(ggpol)
library(ggrepel)
library(lemon)
library(ggsignif)
library(gmodels)
library(margins)
library(lmPerm)
library(permuco)
library(permute)
library(modelr)
library(cowplot)
library(tidyverse)

# Initializing variables ====
rm(list = ls(all = TRUE))
csvroot = 'CSVs/'
figroot = 'Figures/'

# ggplot theme ====
dodge = position_dodge(width = 0.9)
txtSize = 12
gtheme = theme(legend.key = element_rect(colour = "black"),
               legend.title = element_text(size = txtSize),
               legend.text = element_text(size = txtSize),
               legend.background = element_blank(),
               legend.position = "top",
               strip.text.x = element_text(size = txtSize),
               strip.text.y = element_text(size = txtSize),
               axis.text = element_text(colour = "black", size = txtSize),
               plot.title = element_text(size = txtSize),
               axis.title = element_text(size = txtSize),
               axis.line = element_line(),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray"),
               panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "gray"))

#' # Demographics (Figure 1)
# Ages ====
dfAges = read.csv(paste0(csvroot, 'Ages.csv'))
dfAgesControls = dfAges %>% filter(Group %>% str_detect('Control'))
dfAgesYMR = dfAges %>% filter(Group %>% str_detect('YMR'))
agedensity = dfAgesControls %>%  pull(Age) %>% density
#+ fig.width=5.5, fig.height=4, warning=F
p = ggplot() + geom_density(data = dfAgesControls, aes(x = Age, fill = 'Control'),
    color = 'blue', alpha = 0.2) +
  geom_histogram(data = dfAgesControls, aes(x = Age, y = ..density.., fill = 'Control'),
    color = 'blue', alpha = 0.2, binwidth = 3) +
  geom_segment(data = dfAgesYMR, aes(x = Age, y = 0, xend = Age, yend = 0.03, color = 'YMR'), size = 1.5) +
  scale_fill_manual(values = c('blue')) + scale_color_manual(values = c('red')) +
  xlim(agedensity$x %>% range) +
  labs(x = 'Age [y]', y = TeX('Density of samples \\[y$^{-1}$\\]')) +
  gtheme + theme(legend.title = element_blank(),
                 legend.position = c(0.2, 0.9),
                 legend.box = 'horizontal')
p
pdf(paste0(figroot, 'Ages', '.pdf'), width = 5.5,  height = 4.0)
print(p)
dev.off()

# Delta ages =====
dfDeltaAges = read.csv(paste0(csvroot, 'DeltaAges.csv'))
dfDeltaAgesControls = dfDeltaAges %>% filter(Group %>% str_detect('Control'))
dfDeltaAgesYMR = dfDeltaAges %>% filter(Group %>% str_detect('YMR'))
#+ fig.width=5.5, fig.height=4, warning=FALSE
p = ggplot() + geom_density(data = dfDeltaAgesControls, aes(x = dAge, fill = 'Control'),
    color = 'blue', alpha = 0.2) +
  geom_histogram(data = dfDeltaAgesControls, aes(x = dAge, y = ..density.., fill = 'Control'),
    color = 'blue', alpha = 0.2, binwidth = 0.05) +
  geom_segment(data = dfDeltaAgesYMR, aes(x = dAge, xend = dAge, y = 0, yend = 2.5, color = 'YMR'), size = 1.5) +
  scale_fill_manual(values = c('blue')) + scale_color_manual(values = c('red')) +
  scale_x_log10(minor_breaks = pretty_breaks(n = 5), breaks = c(pretty_breaks()(0:1), pretty_breaks()(1:10)), limits = c(0.1, 10)) +
  annotation_logticks(sides = 'b') +
  labs(x = 'Age difference between consecutive visits [y]',
       y = TeX('Density of samples \\[y$^{-1}$\\]')) +
  gtheme + theme(legend.title = element_blank(),
                 legend.position = c(0.8, 0.9),
                 legend.box = 'horizontal')
p
pdf(paste0(figroot, 'DeltaAges', '.pdf'), width = 5.5,  height = 4.0)
print(p)
dev.off()

#' # Noise (Figure 2)
varcsv = read.csv(paste0(csvroot, 'NoiseVariability.csv')) %>% mutate(Year = Year %>% as.factor) %>%
  filter(!(ROIName %in% c('GM', 'WM'))) %>% group_by(ROIName, Year) %>%
  summarise(MeanVolume = ROIVolume %>% mean, SDVolume = ROIVolume %>% sd)

#+ fig.width=4.7, fig.height=4.5, warning=FALSE
p = varcsv %>% mutate(Year = Year %>% as.factor) %>%
  ggplot(aes(x = ROIName, y = MeanVolume, fill = Year)) +
  geom_bar(alpha = 0.5, position = dodge, stat = 'identity') +
  geom_errorbar(aes(x = ROIName, ymin = MeanVolume - SDVolume, ymax = MeanVolume + SDVolume), width = 0.6, position = dodge) +
  labs(x = '', y = TeX('Mean volume \\[mm$^3$\\]')) +
  gtheme + theme(axis.text.x = element_text(angle = 90))
p
pdf(paste0(figroot, 'NoiseVariability', '.pdf'), width = 4.7,  height = 4.5)
print(p)
dev.off()

#' # BrainAGE (Figures 6(b,c), 8, 9)
#' ## Calculations
#' ### Permutation testing function
# Functions ====
GetPVal = function(dset, cg, aform, onesided){
  nperm = 10000 # set number of permutations
  set.seed(987) # set the seed
  perms = shuffleSet(nrow(dset), nperm) # get the permutations for group reshuffling
  null_est_vector =
    apply(perms, 1, function(xx, fset = dset, fform = aform) {
      fset %<>% mutate(Group = Group[xx] %>% factor(levels = c(cg, 'YMR')))
      mod = lm(fform, data = fset)
      null_est = mod$coefficients[4]}) # generate and collect the null distributions of how the YMR slope differs from controls
  
  actmod = lm(aform, data = dset)
  act_stat = actmod$coefficients[[4]]
  return(case_when(onesided ~ sum(null_est_vector <= act_stat)/nperm,
                   T ~ sum(abs(null_est_vector) >= abs(act_stat)) / nperm))
}

#' ### RVM data
# BrainAGE data =============
# generated from Matlab code
df = read.table(paste0(csvroot, 'BrainAGE_Results.csv'), header = FALSE) %>%
  rename(CalendarAge = V3, EstimatedAge = V2, BrainAGE = V1) %>%
  mutate(Group = c('Controls' %>% rep(n() - 4), 'YMR' %>% rep(4)) %>% as.factor)

#' ### Bias correction
# Estimating the bias correction factor ====
set.seed(1)
folds = crossv_loo(df)
folds %<>% mutate(model = map(train, ~ lm(BrainAGE ~ CalendarAge, data = .)))
predicted = folds %>% mutate(predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>% unnest(predicted)
df$BiasCorrectedEstimatedAge = df$EstimatedAge - predicted$.fitted
df$BiasCorrectedBrainAGE = df$BiasCorrectedEstimatedAge - df$CalendarAge

# Performance comparison MAE, RMSE R^2 =====
# MAE
paste('MAE (uncorrected):', abs(df$EstimatedAge - df$CalendarAge) %>% mean)
paste('MAE (bias corrected):', abs(df$BiasCorrectedEstimatedAge - df$CalendarAge) %>% mean)

# RMSE
paste('RMSE (uncorrected):', (df$EstimatedAge - df$CalendarAge) ^ 2 %>% mean %>% sqrt)
paste('RMSE (bias corrected):', (df$BiasCorrectedEstimatedAge - df$CalendarAge) ^ 2 %>% mean %>% sqrt)

# R^2
paste('R^2 (uncorrected):', cor(df$EstimatedAge, df$CalendarAge) ^ 2)
paste('R^2 (bias corrected):', cor(df$BiasCorrectedEstimatedAge, df$CalendarAge) ^ 2)

#' ### CA vs. BrainAGE ($\Delta$)
# CA vs BAGE (uncorrected and bias corrected) models ======
caseq = seq(min(df$CalendarAge), max(df$CalendarAge), length.out = 100)
predvis = data.frame(CalendarAge = caseq)
# uncorrected
lmfit = df %>% lm(BrainAGE ~ CalendarAge, data = .)
pred = predict(lmfit, newdata = predvis, interval = 'prediction')
predvis %<>% mutate(fit = pred[, 1], lwr = pred[, 2], upr = pred[, 3])
# bias corrected
lmfit = df %>% lm(BiasCorrectedBrainAGE ~ CalendarAge, data = .)
pred = predict(lmfit, newdata = predvis, interval = 'prediction')
predvis %<>% mutate(fit_bc = pred[, 1], lwr_bc = pred[, 2], upr_bc = pred[, 3])

# Scatter CA vs. BAGE (uncorrected) prediction intervals ======================
#+ fig.width=4, fig.height=4
p = df %>% mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = BrainAGE, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', se = F, aes(x = CalendarAge, y = BrainAGE),
              size = 1.2, color = 'gray61', fullrange = T, inherit.aes = F) +
  geom_line(data = predvis, aes(x = CalendarAge, y = lwr),
            size = 1.2, color = 'gray61', linetype = 'dashed', inherit.aes = F) +
  geom_line(data = predvis, aes(x = CalendarAge, y = upr),
            size = 1.2, color = 'gray61', linetype = 'dashed', inherit.aes = F) +
  labs(x = "Calendar age [y]", y = TeX("Brain AGE ($\\Delta$) (uncorrected) \\[y\\]")) +
  scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10)) +
  gtheme + theme(legend.position = c(0.85, 0.9))
p
pdf(paste0(figroot, 'CAvsBAGE_uncorrected', '.pdf'), width = 4.0,  height = 4.0)
print(p)
dev.off()

# Scatter CA vs. BAGE (bias corrected) prediction intervals ======================
#+ fig.width=4, fig.height=4
p = df %>% mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = BiasCorrectedBrainAGE, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', se = F, aes(x = CalendarAge, y = BiasCorrectedBrainAGE),
              size = 1.2, color = 'gray61', fullrange = T, inherit.aes = F) +
  geom_line(data = predvis, aes(x = CalendarAge, y = lwr_bc),
            size = 1.2, color = 'gray61', linetype = 'dashed', inherit.aes = F) +
  geom_line(data = predvis, aes(x = CalendarAge, y = upr_bc),
            size = 1.2, color = 'gray61', linetype = 'dashed', inherit.aes = F) +
  labs(x = "Calendar age [y]", y = TeX("Brain AGE ($\\Delta$) (bias corrected) \\[y\\]")) +
  scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10)) +
  gtheme + theme(legend.position = c(0.85, 0.9))
p
pdf(paste0(figroot, 'CAvsBAGE_biascorrected', '.pdf'), width = 4.0,  height = 4.0)
print(p)
dev.off()

# Control subgroup data =====
controls = read.csv(paste0(csvroot, 'SubjIDVisit.csv')) %>% filter(!(SubjectID %>% str_detect('MYR')))
controlsgrp = read.csv(paste0(csvroot, 'ControlSubjectGroups.csv')) %>% mutate(SubjectID = paste0('WB', studyID)) %>% rename(ControlGroup = Group)
dfControls = df %>% filter(Group %>% str_detect('Controls')) %>% cbind(inner_join(controls, controlsgrp)) %>%
  mutate(ControlGroup = ifelse(Visit == 'Visit1', 'WL', ControlGroup %>% as.character) %>% as.factor)
dfYMR = df %>% filter(Group == 'YMR') %>% mutate(SubjectID = 'YMRID', Visit = c('Visit1', 'Visit2', 'Visit3', 'Visit4'), studyID = 'YMRSID')
dfWithControlGroups = dfControls %>% group_by(ControlGroup) %>% group_modify(~ rbind(.x, dfYMR)) %>%
  group_by(ControlGroup) %>%
  mutate(CGSamp = paste0(ControlGroup, ' n=(', n() - 4, ')'))
dfWithControlGroups$CGSamp %<>% as.factor() %>% fct_shift(-1)

#' ### YMR vs. controls (uncorrected)
# CA vs BA (uncorrected) models ======
mdls = df %>% group_by(Group) %>% do(fit = lm(EstimatedAge ~ CalendarAge, data = .) %>% tidy) %>%
  rowwise() %>% do(estimates = .$fit %>% select(term, estimate) %>% spread(term, estimate) %>%
                     rename(`Starting brain age [y]` = `(Intercept)`, `Rate of brain age change` = CalendarAge),
                   SEs = .$fit %>% select(term, std.error) %>% spread(term, std.error) %>%
                     rename(BACAInterceptSE = `(Intercept)`, BACASlopeSE = CalendarAge),
                   Group = .$Group)

dfGroupBAMdls = cbind(do.call(rbind, mdls$estimates), do.call(rbind, mdls$SEs), Group = mdls$Group %>% unlist) %>%
  mutate(x1 = c(0, 0), y1 = `Starting brain age [y]`,
         x2 = c(27, 27), y2 = `Rate of brain age change` * 27 + `Starting brain age [y]`,
         x3 = c(75, 75), y3 = `Rate of brain age change` * 75 + `Starting brain age [y]`)

#' ### YMR vs. controls (bias corrected)
# CA vs BA (bias corrected) models ======
mdls = df %>% group_by(Group) %>% do(fit = lm(BiasCorrectedEstimatedAge ~ CalendarAge, data = .) %>% tidy) %>%
  rowwise() %>% do(estimates = .$fit %>% select(term, estimate) %>% spread(term, estimate) %>% 
                     rename(`Starting brain age [y]` = `(Intercept)`, `Rate of brain age change` = CalendarAge),
                   SEs = .$fit %>% select(term, std.error) %>% spread(term, std.error) %>%
                     rename(BACAInterceptSE = `(Intercept)`, BACASlopeSE = CalendarAge),
                   Group = .$Group)

dfGroupBABCMdls = cbind(do.call(rbind, mdls$estimates), do.call(rbind, mdls$SEs), Group = mdls$Group %>% unlist) %>%
  mutate(x1 = c(0, 0), y1 = `Starting brain age [y]`,
         x2 = c(27, 27), y2 = `Rate of brain age change` * 27 + `Starting brain age [y]`,
         x3 = c(75, 75), y3 = `Rate of brain age change` * 75 + `Starting brain age [y]`)

#' ### YMR vs. subgroups (uncorrected)
# CA vs BA (uncorrected) models (by control subgroups) ======
ctrlmdls = dfControls %>% group_by(ControlGroup) %>%
  group_modify( ~ lm(EstimatedAge ~ CalendarAge, data = .x) %>% tidy) %>%
  select(ControlGroup, term, estimate) %>% spread(term, estimate) %>%
  mutate(x1 = c(0), y1 = `(Intercept)`, 
         x2 = c(27), y2 = CalendarAge * 27 + `(Intercept)`,
         x3 = c(75), y3 = CalendarAge * 75 + `(Intercept)`)

ymrmdl = dfYMR %>% lm(EstimatedAge ~ CalendarAge, data = .) %>% tidy %>% select(term, estimate) %>% spread(term, estimate) %>% 
  mutate(x1 = c(0), y1 = `(Intercept)`,
         x2 = c(27), y2 = CalendarAge * 27 + `(Intercept)`,
         x3 = c(75), y3 = CalendarAge * 75 + `(Intercept)`,
         Group = 'YMR')

dfmdlsUC = ctrlmdls %>% mutate(Group = "Controls") %>% group_by(ControlGroup) %>% group_modify(~ rbind(.x, ymrmdl)) %>%
  inner_join(dfWithControlGroups %>% distinct(ControlGroup, Group, CGSamp))

#' ### YMR vs. subgroups (bias corrected)
# CA vs BA (bias corrected) models (by control subgroups) ======
ctrlmdls = dfControls %>% group_by(ControlGroup) %>%
  group_modify( ~ lm(BiasCorrectedEstimatedAge ~ CalendarAge, data = .x) %>% tidy) %>%
  select(ControlGroup, term, estimate) %>% spread(term, estimate) %>%
  mutate(x1 = c(0), y1 = `(Intercept)`,
         x2 = c(27), y2 = CalendarAge * 27 + `(Intercept)`,
         x3 = c(75), y3 = CalendarAge * 75 + `(Intercept)`)

ymrmdl = dfYMR %>% lm(BiasCorrectedEstimatedAge ~ CalendarAge, data = .) %>% tidy %>% select(term, estimate) %>% spread(term, estimate) %>% 
  mutate(x1 = c(0), y1 = `(Intercept)`,
         x2 = c(27), y2 = CalendarAge * 27 + `(Intercept)`,
         x3 = c(75), y3 = CalendarAge * 75 + `(Intercept)`,
         Group = 'YMR')

dfmdls = ctrlmdls %>% mutate(Group = "Controls") %>% group_by(ControlGroup) %>% group_modify(~ rbind(.x, ymrmdl)) %>%
  inner_join(dfWithControlGroups %>% distinct(ControlGroup, Group, CGSamp))

#' ### Permutation testing (YMR vs. controls)
# CA vs. BA (uncorrected) permutation testing ===================
slopepval_onesidedUC = GetPVal(df, 'Controls', as.formula("EstimatedAge ~ CalendarAge * Group"), T)

# CA vs. BA (bias corrected) permutation testing ===================
slopepval_onesided = GetPVal(df, 'Controls', as.formula("BiasCorrectedEstimatedAge ~ CalendarAge * Group"), T)

#' ### Permutation testing (YMR vs. subgroups)
# CA vs. BA (uncorrected) permutation testing (by control subgroups) ===================
dfPValUC = dfWithControlGroups %>% group_by(ControlGroup) %>% group_nest() %>%
  mutate(pval = map2_dbl(data, ControlGroup, ~GetPVal(.x %>% mutate(Group = case_when(Group == 'Controls' ~ as.character(.y), T ~ 'YMR') %>% as.factor),
                                                      as.character(.y), as.formula("EstimatedAge ~ CalendarAge * Group"), T)))

# CA vs. BA (bias corrected) permutation testing (by control subgroups) ===================
dfPVal = dfWithControlGroups %>% group_by(ControlGroup) %>% group_nest() %>%
  mutate(pval = map2_dbl(data, ControlGroup, ~GetPVal(.x %>% mutate(Group = case_when(Group == 'Controls' ~ as.character(.y), T ~ 'YMR') %>% as.factor),
                                                      as.character(.y), as.formula("BiasCorrectedEstimatedAge ~ CalendarAge * Group"), T)))

#' ## Visualizations
#' ### YMR vs. controls (uncorrected)
# Scatter CA vs. BA (uncorrected) ==================
#+ fig.width=4, fig.height=4
p = df %>% mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = EstimatedAge, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', se = F, size = 0.8, fullrange = F) +
  geom_segment(aes(x = 25, y = 25, xend = 68, yend = 68),
               color = 'gray61', size = 0.75, alpha = 1.0) +
  geom_segment(aes(x = 0, y = 0, xend = 25, yend = 25),
               color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 68, y = 68, xend = 75, yend = 75), arrow = arrow(length = unit(0.03, "npc"), type = "closed"),
               color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(data = dfGroupBAMdls, aes(x = x1, y = y1, xend = x2, yend = y2, color = Group),
               linetype = 2, size = 0.75, alpha = 0.75) +
  geom_segment(data = dfGroupBAMdls, aes(x = x2, y = y2, xend = x3, yend = y3, color = Group),
               linetype = 2, size = 0.75, alpha = 0.75) +
  geom_text(x = 10, y = 60, label = paste('p =', slopepval_onesidedUC %>% round(digits = 4)),
            size = txtSize / 2, color = 'black', hjust = -0.1, vjust = -1, inherit.aes = F) +
  scale_x_continuous(breaks = pretty_breaks(10), limits = c(0, 75)) + scale_y_continuous(breaks = pretty_breaks(10), limits = c(0, 75)) +
  labs(x = "Calendar age [y]", y = "Estimated brain age (uncorrected) [y]") +
  gtheme + theme(legend.position = c(0.8, 0.2),
                 legend.box = 'horizontal')
p
pdf(paste0(figroot, 'CAvsBA_uncorrected', '.pdf'), width = 4.0,  height = 4.0)
print(p)
dev.off()

#' ### YMR vs. controls (bias corrected)
# Scatter CA vs. BA (bias corrected) ========
#+ fig.width=4, fig.height=4
p = df %>% mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = BiasCorrectedEstimatedAge, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', size = 0.8, se = F, fullrange = F) +
  geom_segment(aes(x = 25, y = 25, xend = 68, yend = 68),
               color = 'gray61', size = 0.75, alpha = 1.0) +
  geom_segment(aes(x = 0, y = 0, xend = 25, yend = 25),
               color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 68, y = 68, xend = 75, yend = 75),
               color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2,
               arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +
  geom_segment(data = dfGroupBABCMdls, aes(x = x1, y = y1, xend = x2, yend = y2, color = Group),
               linetype = 2, size = 0.75, alpha = 0.75) +
  geom_segment(data = dfGroupBABCMdls, aes(x = x2, y = y2, xend = x3, yend = y3, color = Group),
               size = 0.75, alpha = 0.75, linetype = 2) +
  geom_text(x = 10, y = 60, label = paste('p =', slopepval_onesided %>% round(digits = 4)),
            inherit.aes = F, size = txtSize / 2, color = 'black', hjust = -0.1, vjust = -1) +
  scale_x_continuous(breaks = pretty_breaks(10), limits = c(0, 75)) + scale_y_continuous(breaks = pretty_breaks(10), limits = c(0, 75)) +
  labs(x = "Calendar age [y]", y = "Estimated brain age (bias corrected) [y]") +
  gtheme + theme(legend.position = c(0.8, 0.2),
                 legend.box = 'horizontal')
p
pdf(paste0(figroot, 'CAvsBA_biascorrected', '.pdf'), width = 4.0,  height = 4.0)
print(p)
dev.off()

#' ### YMR vs. subgroups (uncorrected)
# Scatter CA vs. BA (uncorrected) by control subgroups ============
#+ fig.width=8, fig.height=3.15
p = dfWithControlGroups %>% ungroup() %>% mutate(ControlGroup = fct_relevel(ControlGroup, "WL")) %>%
  mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = EstimatedAge, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', size = 0.8, se = F, fullrange = F) +
  geom_text(data = dfPValUC %>% mutate(CGSamp = map2_chr(data, ControlGroup, ~ 
                                                           paste0(.y, ' n=(', nrow(.x) - 4, ')')) %>% as.factor %>% fct_shift(-1)),
    aes(x = rep(30, 3), y = rep(5, 3), label = paste('p =', pval %>% round(digits = 4))),
    inherit.aes = F, size = txtSize / 2, color = 'black', hjust = -0.1, vjust = -1) +
  geom_segment(data = dfmdlsUC, aes(x = x1, y = y1, xend = x2, yend = y2, color = Group),
    linetype = 2, size = 0.75, alpha = 0.75) +
  geom_segment(data = dfmdlsUC, aes(x = x2, y = y2, xend = x3, yend = y3, color = Group),
    size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 25, y = 25, xend = 68, yend = 68), color = 'gray61', size = 0.75, alpha = 1.0) +
  geom_segment(aes(x = 0, y = 0, xend = 25, yend = 25), color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 68, y = 68, xend = 75, yend = 75), arrow = arrow(length = unit(0.03, "npc"), type = "closed"), 
               color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  facet_rep_wrap(. ~ CGSamp) +
  scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10)) +
  labs(x = "Calendar age [y]", y = "Estimated brain age [y]") +
  gtheme + theme(legend.position = c(0.1, 0.9))
p
pdf(paste0(figroot, 'CAvsBA_uncorrected_subgroups', '.pdf'), width = 8.0,  height = 3.15)
print(p)
dev.off()

#' ### YMR vs. subgroups (bias corrected)
# Scatter CA vs. BA (bias corrected) by control subgroups =======
#+ fig.width=8, fig.height=3.15
p = dfWithControlGroups %>% ungroup() %>%
  mutate(ControlGroup = ControlGroup %>% fct_relevel("WL"),
         Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = CalendarAge, y = BiasCorrectedEstimatedAge, color = Group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm', size = 0.8, se = F, fullrange = F) +
  geom_text(data = dfPVal %>% mutate(CGSamp = map2_chr(data, ControlGroup, ~ 
                                                         paste0(.y, ' n=(', nrow(.x) - 4, ')')) %>% as.factor %>% fct_shift(-1)),
    aes(x = rep(30, 3), y = rep(5, 3), label = paste('p =', pval %>% round(digits = 4))),
    inherit.aes = F, size = txtSize / 2, color = 'black', hjust = -0.1, vjust = -1) +
  geom_segment(data = dfmdls, aes(x = x1, y = y1, xend = x2, yend = y2, color = Group),
    linetype = 2, size = 0.75, alpha = 0.75) +
  geom_segment(data = dfmdls, aes(x = x2, y = y2, xend = x3, yend = y3, color = Group),
    size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 25, y = 25, xend = 68, yend = 68),
    color = 'gray61', size = 0.75, alpha = 1.0) +
  geom_segment(aes(x = 0, y = 0, xend = 25, yend = 25),
    color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  geom_segment(aes(x = 68, y = 68, xend = 75, yend = 75), arrow = arrow(length = unit(0.03, "npc"), type = "closed"),
    color = 'gray61', size = 0.75, alpha = 0.75, linetype = 2) +
  facet_rep_wrap(. ~ CGSamp) +
  scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10)) +
  labs(x = "Calendar age [y]", y = "Estimated brain age (bias corrected) [y]") +
  gtheme + theme(legend.position = c(0.1, 0.9))
p
pdf(paste0(figroot, 'CAvsBA_biascorrected_subgroups', '.pdf'), width = 8.0,  height = 3.15)
print(p)
dev.off()

#' # Brain resemblance (Figure 10, Table 2)
kdf = data.frame(SubjectID = 'YMRID', k = c(1, 3, 5, 7, 11, 13, 17))
knnage = dfYMR %>% nest_join(kdf, name = 'kdf') %>% unnest(kdf) %>%
  mutate(ControlAge = map2_dbl(k, BiasCorrectedEstimatedAge, ~
                                 dfControls %>% mutate(AbsEstDiff = abs(BiasCorrectedEstimatedAge - .y)) %>% 
                                 arrange(AbsEstDiff) %>% head(.x) %>% pull(CalendarAge) %>% mean))
knnagesummary = knnage %>% group_by(CalendarAge) %>% summarise(Mean = mean(ControlAge) %>% round(2) %>% format(nsmall = 2), 
                                                      SD = sd(ControlAge) %>% round(2) %>% format(nsmall = 2))
knnagesummary
# kNN calendar age histogram (Figure 9) ======
#+ fig.width=5, fig.height=4.5
knnagevis = knnage %>% filter(CalendarAge == 41) %>% select(k, ControlAge)
p = ggplot() + geom_histogram(data = df %>% filter(Group %>% str_detect('Controls')), aes(x = CalendarAge), binwidth = 3,
                              color = 'black', fill = 'white', alpha = 0.5) +
  geom_segment(data = knnagevis, aes(x = ControlAge, xend = ControlAge, y = 0, yend = 20, color = k)) +
  geom_point(data = knnagevis, aes(x = ControlAge, y = 20), color = "red", size = 0.85) +
  geom_text_repel(data = knnagevis, aes(x = ControlAge, y = 20, label = ControlAge %>% round(2) %>% format(nsmall = 2)),
                  nudge_y = 1.0, direction = "x", angle = 90, vjust = 0, segment.size = 0.2) +
  guides(color = guide_legend(nrow = 1, label.position = 'bottom')) +
  annotate("segment", x = 24, xend = 41, y = 25, yend = 25, colour = "black", size = 1.25) +
  annotate("text", x = 33, y = 26.25, label = paste(knnagesummary$Mean[[4]], '+/-', knnagesummary$SD[[4]], '[y]'), parse = F) +
  scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10)) +
  labs(x = 'Calendar age [y]', y = 'Number of control samples') +
  gtheme + theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_line(size = 0.0))
p
pdf(paste0(figroot, 'knnHist', '.pdf'), width = 5.0,  height = 4.5)
print(p)
dev.off()

# kNN projection ====
k = 11
knnprojvis = dfControls %>% mutate(AbsEstDiff = abs(BiasCorrectedEstimatedAge - dfYMR$BiasCorrectedEstimatedAge[[4]])) %>% arrange(AbsEstDiff) %>% head(k)
#+ fig.width=4, fig.height=4
p = dfControls %>% ggplot() + geom_point(aes(x = CalendarAge, y = BiasCorrectedEstimatedAge),
                                         color = 'blue', size = 2, alpha = 0.5) +
  geom_point(aes(x = df$CalendarAge[[243]], y = df$BiasCorrectedEstimatedAge[[243]]),
             size = 5, color = 'orange', alpha = 0.25) +
  geom_vline(xintercept = df$CalendarAge[[243]]) +
  geom_segment(data = knnprojvis, aes(x = CalendarAge, y = BiasCorrectedEstimatedAge, xend = dfYMR$CalendarAge[[4]], yend = BiasCorrectedEstimatedAge),
               linetype = 1, size = 0.2) +
  geom_segment(data = knnprojvis, aes(x = CalendarAge, y = BiasCorrectedEstimatedAge, xend = CalendarAge, yend = 23),
               linetype = 1, size = 0.2, arrow = arrow(length = unit(0.03, "npc"))) +
  geom_label(aes(x = 53, y = 35, label = "k = 11\n nearest neighbors"),
             fill = 'gray',alpha = 0.25) +
  annotate('rect', xmin = dfYMR$CalendarAge[[4]] - 1, xmax = dfYMR$CalendarAge[[4]] + 1, ymin = knnprojvis$BiasCorrectedEstimatedAge %>% min, ymax = knnprojvis$BiasCorrectedEstimatedAge %>% max, alpha = 0.5) +
  geom_curve(aes(x = 52, y = 37.3, xend = dfYMR$CalendarAge[[4]] + 1, yend = knnprojvis$BiasCorrectedEstimatedAge %>% max),
             curvature = 0.3, arrow = arrow(length = unit(0.03, "npc"))) +
  xlim(23, 68) + ylim(23, 68) + scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks(25)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Calendar age [y]', y = 'Estimated brain age [y]') +
  gtheme + theme(axis.text.x = element_text(angle = 90),
                 panel.grid.minor.x = element_blank())
p
pdf(paste0(figroot, 'knnProj', '.pdf'), width = 4.0,  height = 4.0)
print(p)
dev.off()

#' # ROI volumetry (Figure 11)
# ROI data =================
dfROI = read.csv(paste0(csvroot, 'NormalizedVolumes.csv'), na.strings = c('', ' ', 'NA')) %>%
  spread(MeasureName, MeasureValue) %>% gather(ROIName, NormROIVolume, Cingulum_Ant_L:Temporal_Inf_R) %>% 
  mutate(ROIName = ROIName %>% as.factor) %>% filter(Age >= 27 & Age <= 41)

# Vol vs. age permutation testing =======
dfPValROI = dfROI %>% group_by(ROIName) %>% group_nest() %>% mutate(pval = map_dbl(data, ~GetPVal(.x, 'Controls', as.formula("NormROIVolume ~ Age * Group"), F)))

# ROIVol ~ Age ====
#+ fig.width=9.5, fig.height=9.5
p = dfROI %>% mutate(Group = Group %>% fct_relevel('Controls', after = 2)) %>%
  ggplot(aes(x = Age, y = NormROIVolume, color = Group)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = 'lm', se = F, fullrange = F) +
  geom_text(data = dfPValROI, aes(x = -Inf, y = -Inf, label = paste('p = ', pval)),
            color = 'black', size = 6, hjust = -0.1, vjust = -1) +
  facet_rep_wrap(~ROIName, scales = 'free', ncol = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(x = 'Calendar age [y]', y = 'Regional volume fraction of the gray matter') +
  gtheme + theme(legend.title = element_blank(),
                 axis.title = element_text(size = txtSize + 4),
                 legend.position = c(0.7, 0.1),
                 legend.background = element_blank())
p
pdf(paste0(figroot, 'ROIScatter', '.pdf'), width = 9.5,  height = 9.5)
print(p)
dev.off()
