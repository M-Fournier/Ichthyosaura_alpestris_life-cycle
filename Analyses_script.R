## ------------------------------------------------------------------------- ##
##        Testing the impact of life cycle, size and their interaction       ##
##                     on body morphology and bones shape                    ##
##                        using multivariate analyses                        ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(mvMORPH)
library(geomorph)
library(dplyr)

# Data import
#_______________________________________________________________________________
df_linear <-read.csv("Data/Linear_measurements/df_linear_filtered.csv", 
                     sep=",",
                     row.names = 1) # ordi pro

#

load(file = "Data/Landmarks/pts/Y.gpa_cr.Rdata")
Y_cr <- Y.gpa_cr$rotated
library(geomorph)
Y_cr_2D <- two.d.array(Y_cr)

df_Csize_cr <- read.csv("Data/Landmarks/Ident/ident_&_Csize_cr.csv",
                        sep=",",header = T, row.names = 1)

df_Csize_cr <- df_Csize_cr[match(rownames(Y_cr_2D), rownames(df_Csize_cr)), ]

all(rownames(df_Csize_cr) == rownames(Y_cr_2D))  

#

load(file = "Data/Landmarks/pts/Y.gpa_md.Rdata")
Y_md <- Y.gpa_md$rotated
library(geomorph)
Y_md_2D <- two.d.array(Y_md)

df_Csize_md <- read.csv("Data/Landmarks/Ident/ident_&_Csize_md.csv",
                        sep=",",header = T, row.names = 1)

all(rownames(Y_md_2D) == rownames(df_Csize_md))

# Fitting of models 
#_______________________________________________________________________________

Y_linear <- as.matrix(df_linear[, c(7:11, 21:30)])
Y_cr <- as.matrix(Y_cr_2D)
Y_md <- as.matrix(Y_md_2D)

library(mvMORPH)

fit_linear <- mvols(
  Y_linear ~ life_cycle * svl + country,
  data = df_linear,  
  method = "PL-LOOCV",
  iter = 999
)

fit_cr <- mvols(
  Y_cr  ~ life_cycle * Csize + country,
  data = df_Csize_cr,
  method = "PL-LOOCV", 
  iter = 999 
)

fit_md <- mvols(
  Y_md  ~ life_cycle * Csize + country,
  data = df_Csize_md,  
  method = "PL-LOOCV", 
  iter = 999 
)

{
save(fit_linear, file = "//fit_linear.Rdata")
save(fit_cr, file = "/fit_cr.Rdata")
save(fit_md, file = "/fit_md.Rdata")
}

# Manovas
#_______________________________________________________________________________
results_linear <- manova.gls(fit_linear, 
                             test = "Pillai", 
                             type="II", 
                             nperm = 999,
                             seed = 123)
results_linear
# Type II MANOVA Tests with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat)   
# life_cycle        0.9356     0.001 **
#   svl               0.9152     0.001 **
#   country           0.5822     0.001 **
#   life_cycle:svl    0.3560     0.001 **
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

effect_size_linear <- effectsize(results_linear)
print(effect_size_linear)
# ## Multivariate measure(s) of association ## 
# life_cycle    svl country life_cycle:svl
# ξ^2 [Pillai]     0.4227 0.9082   0.231          0.108

pairwise_linear <- pairwise.glh(fit_linear, term="life_cycle",
             test="Pillai",
             type="II",
             adjust="holm",
             nperm=999, 
             verbose=TRUE,
             seed = 123)

pairwise_linear
# General Linear Hypothesis Test with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat) adjusted   
# metamorphic - paedomorphic                0.29518     0.001    0.003 **
#   metamorphic - undergoing metamorphosis    0.04594     0.909    0.909   
# paedomorphic - undergoing metamorphosis   0.17511     0.006    0.012  *
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

effect_size_pairwise_linear <- effectsize(pairwise_linear)
print(effect_size_pairwise_linear)
# ## Multivariate measure(s) of association ## 
# metamorphic - paedomorphic | metamorphic - undergoing metamorphosis | paedomorphic - undergoing metamorphosis |
#   ξ^2 [Pillai]                       0.2353                                 -0.03752                                    0.1023
# ## Values < 0 represent no association    ##

results_cr <- manova.gls(fit_cr, 
                         test = "Pillai", 
                         type="II", 
                         nperm = 999,
                         seed = 123,
                         verbose=TRUE)
results_cr
# Type II MANOVA Tests with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat)   
# life_cycle          1.6971     0.001 **
#   Csize               0.9307     0.001 **
#   country             1.6869     0.001 **
#   life_cycle:Csize    1.3098     0.001 **
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

effect_size_cr <- effectsize(results_cr)
print(effect_size_cr)
# ## Multivariate measure(s) of association ## 
# life_cycle  Csize country life_cycle:Csize
# ξ^2 [Pillai]     0.6856 0.8528  0.6685           0.2622

pairwise_cr <- pairwise.glh(fit_cr, term="life_cycle",
                                test="Pillai",
                                type="II",
                                adjust="holm",
                                nperm=999, 
                                verbose=TRUE,
                                seed = 123)
pairwise_cr
# General Linear Hypothesis Test with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat) adjusted   
# metamorphic - paedomorphic                 0.7183     0.001    0.003 **
#   metamorphic - undergoing metamorphosis     0.5605     0.237    0.237   
# paedomorphic - undergoing metamorphosis    0.6642     0.001    0.003 **
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

effect_size_pairwise_cr <- effectsize(pairwise_cr)
print(effect_size_pairwise_cr)
# ## Multivariate measure(s) of association ## 
# metamorphic - paedomorphic | metamorphic - undergoing metamorphosis | paedomorphic - undergoing metamorphosis |
#   ξ^2 [Pillai]                       0.3962                                  0.05644                                    0.2778

#

results_md <- manova.gls(fit_md, 
                         test = "Pillai", 
                         type="II", 
                         nperm = 999,
                         seed = 123)
results_md
# Type II MANOVA Tests with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat)   
# life_cycle          1.5150     0.001 **
#   Csize               0.7668     0.001 **
#   country             1.3358     0.001 **
#   life_cycle:Csize    1.2223     0.001 **
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

effect_size_md <- effectsize(results_md)
print(effect_size_md)
# ## Multivariate measure(s) of association ## 
# life_cycle  Csize country life_cycle:Csize
# ξ^2 [Pillai]     0.5613 0.5776  0.3806           0.2772

pairwise_md <- pairwise.glh(fit_md, term="life_cycle",
                            test="Pillai",
                            type="II",
                            adjust="holm",
                            nperm=999, 
                            verbose=TRUE,
                            seed = 123)
pairwise_md
# General Linear Hypothesis Test with 999 permutations: Pillai test statistic 
# Test stat Pr(>Stat) adjusted   
# metamorphic - paedomorphic                 0.6771     0.001    0.003 **
#   metamorphic - undergoing metamorphosis     0.5200     0.166    0.166   
# paedomorphic - undergoing metamorphosis    0.5853     0.021    0.042  *
#   --- 
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1  

effect_size_pairwise_md <- effectsize(pairwise_md)
print(effect_size_pairwise_md)
# ## Multivariate measure(s) of association ## 
# metamorphic - paedomorphic | metamorphic - undergoing metamorphosis | paedomorphic - undergoing metamorphosis |
#   ξ^2 [Pillai]                       0.3966                                  0.09394                                    0.2217

{
save(results_linear, file = "/results_linear.Rdata")
save(pairwise_linear, file = "/pairwise_linear.Rdata")
save(results_cr, file = "/results_cr.Rdata")
save(pairwise_cr, file = "/pairwise_cr.Rdata")
save(results_md, file = "/results_md.Rdata")
save(pairwise_md, file = "/pairwise_md.Rdata")
}

# Export results
#_______________________________________________________________________________
{
res_linear <- data.frame(
  term = as.character(results_linear$terms),
  Pillai = as.numeric(results_linear$stat),
  p_value = as.numeric(results_linear$pvalue),
  effect_size = as.numeric(effect_size_linear$effect)
)
res_cr <- data.frame(
  term = as.character(results_cr$terms),
  Pillai = as.numeric(results_cr$stat),
  p_value = as.numeric(results_cr$pvalue),
  effect_size = as.numeric(effect_size_cr$effect)
)
res_md <- data.frame(
  term = as.character(results_md$terms),
  Pillai = as.numeric(results_md$stat),
  p_value = as.numeric(results_md$pvalue),
  effect_size = as.numeric(effect_size_md$effect)
)

res_linear_pairwise <- data.frame(
  comparison = colnames(effect_size_pairwise_linear$effect),
  Pillai = as.numeric(pairwise_linear$stat),
  p_value_adjusted = as.numeric(pairwise_linear$adjust),
  effect_size = as.numeric(effect_size_pairwise_linear$effect)
)
res_cr_pairwise <- data.frame(
  comparison = colnames(effect_size_pairwise_cr$effect),
  Pillai = as.numeric(pairwise_cr$stat),
  p_value_adjusted = as.numeric(pairwise_cr$adjust),
  effect_size = as.numeric(effect_size_pairwise_cr$effect)
)
res_md_pairwise <- data.frame(
  comparison = colnames(effect_size_pairwise_md$effect),
  Pillai = as.numeric(pairwise_md$stat),
  p_value_adjusted = as.numeric(pairwise_md$adjust),
  effect_size = as.numeric(effect_size_pairwise_md$effect)
)
}

{
  write.csv(res_linear, file = "/res_linear.csv",
          row.names = F)
  write.csv(res_cr, file = "/res_cr.csv",
            row.names = F)
  write.csv(res_md, file = "/res_md.csv",
            row.names = F)
  write.csv(res_linear_pairwise, file = "/res_linear_pairwise.csv",
            row.names = F)
  write.csv(res_cr_pairwise, file = "/res_cr_pairwise.csv",
            row.names = F)
  write.csv(res_md_pairwise, file = "/res_md_pairwise.csv",
            row.names = F)
}

###############################################################################
## ------------------------------------------------------------------------- ##
##                  Testing for allometric relationship                      ##
##                   between body/bones shape and size                       ##
##                       depending on life cycle                             ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(mvMORPH)
library(ggplot2)
library(ggrepel)

# Data import
#_______________________________________________________________________________
df_linear <-read.csv("Data/Linear_measurements/df_linear_filtered.csv", 
                     sep=",",
                     row.names = 1) 

#

load(file = "Data/Landmarks/pts/Y.gpa_cr.Rdata")
Y_cr <- Y.gpa_cr$rotated
library(geomorph)
Y_cr_2D <- two.d.array(Y_cr)

df_Csize_cr <- read.csv("Data/Landmarks/Ident/ident_&_Csize_cr.csv",
                        sep=",",header = T, row.names = 1)


df_Csize_cr <- df_Csize_cr[match(rownames(Y_cr_2D), rownames(df_Csize_cr)), ]

all(rownames(df_Csize_cr) == rownames(Y_cr_2D))  


#

load(file = "Data/Landmarks/pts/Y.gpa_md.Rdata")
Y_md <- Y.gpa_md$rotated
library(geomorph)
Y_md_2D <- two.d.array(Y_md)

df_Csize_md <- read.csv("Data/Landmarks/Ident/ident_&_Csize_md.csv",
                        sep=",",header = T, row.names = 1)


all(rownames(Y_md_2D) == rownames(df_Csize_md))


# Fitting of models 
#_______________________________________________________________________________
Y_linear <- as.matrix(df_linear[, c(7:11, 21:30)])
Y_cr <- as.matrix(Y_cr_2D)
Y_md <- as.matrix(Y_md_2D)

fit_linear_allo <- mvols(
  Y_linear ~ life_cycle * svl,
  data = df_linear,  
  method = "PL-LOOCV",
  iter = 999
)

fit_cr_allo <- mvols(
  Y_cr  ~ life_cycle * Csize,
  data = df_Csize_cr,
  method = "PL-LOOCV", 
  iter = 999 
)

fit_md_allo <- mvols(
  Y_md  ~ life_cycle * Csize,
  data = df_Csize_md,  
  method = "PL-LOOCV", 
  iter = 999 
)

{
save(fit_linear_allo, file = "/fit_linear_allo.Rdata")
save(fit_cr_allo, file = "/fit_cr_allo.Rdata")
save(fit_md_allo, file = "/fit_md_allo.Rdata")
}

# Simple plots and extraction of allometric data
#_______________________________________________________________________________
plot_linear <- plot(fit_linear_allo,
                    term = "svl",
                    residuals = FALSE,
                    fitted = TRUE,
                    pch = 19)
plot_cr <- plot(fit_cr_allo,
                term = "Csize",
                residuals = FALSE,
                fitted = TRUE,
                pch = 19)
plot_md <- plot(fit_md_allo,
                term = "Csize",
                residuals = FALSE,
                fitted = TRUE,
                pch = 19)

{
  mv_scores_linear <- (plot_linear$scores)
  fitted <- fitted(fit_linear_allo)%*%plot_linear$standBeta
  mv_scores_linear <- as.data.frame(mv_scores_linear)
  fitted <- as.data.frame(fitted)
  mv_scores_linear$fitted_values <- fitted[, 1]
  library(dplyr)
  mv_scores_linear <- merge(mv_scores_linear,
                            df_linear,
                            by = 'row.names',
                            all = T)
  row.names(mv_scores_linear) <- mv_scores_linear$Row.names
  mv_scores_linear$Row.names <- NULL
}
{
  mv_scores_cr <- (plot_cr$scores)
  fitted <- fitted(fit_cr_allo)%*%plot_cr$standBeta
  mv_scores_cr <- as.data.frame(mv_scores_cr)
  fitted <- as.data.frame(fitted)
  mv_scores_cr$fitted_values <- fitted[, 1]
  library(dplyr)
  mv_scores_cr <- merge(mv_scores_cr,
                        df_Csize_cr,
                        by = 'row.names',
                        all = T)
  row.names(mv_scores_cr) <- mv_scores_cr$Row.names
  mv_scores_cr$Row.names <- NULL
}
{
  mv_scores_md <- (plot_md$scores)
  fitted <- fitted(fit_md_allo)%*%plot_md$standBeta
  mv_scores_md <- as.data.frame(mv_scores_md)
  fitted <- as.data.frame(fitted)
  mv_scores_md$fitted_values <- fitted[, 1]
  library(dplyr)
  mv_scores_md <- merge(mv_scores_md,
                        df_Csize_md,
                        by = 'row.names',
                        all = T)
  row.names(mv_scores_md) <- mv_scores_md$Row.names
  mv_scores_md$Row.names <- NULL
}

# Plots
#_______________________________________________________________________________
plot_allo_linear <- ggplot(mv_scores_linear, aes(x = svl, y = V1, color = life_cycle)) +
  geom_point() +
  labs(x = "Log (svl)", y = "mv scores") +
  geom_line(aes(y = fitted_values, color = life_cycle),
            linewidth = 1) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "metamorphic" = "#009900", "undergoing metamorphosis" = "#00CCFF")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid = element_line(color = "gray90"),
        axis.line = element_line(color = "black", size = 0.5),
        aspect.ratio = 1)
plot_allo_linear

plot_allo_cr <- ggplot(mv_scores_cr, aes(x = Csize, y = V1, color = life_cycle)) +
  geom_point() +
  labs(x = "Csize", y = "mv scores") +
  geom_line(aes(y = fitted_values, color = life_cycle),
            linewidth = 1) +
  # geom_text_repel(aes(label = rownames(mv_scores_cr)),
  #                 size = 3,
  #                 max.overlaps = Inf) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "metamorphic" = "#009900", "undergoing metamorphosis" = "#00CCFF")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid = element_line(color = "gray90"),
        axis.line = element_line(color = "black", size = 0.5),
        aspect.ratio = 1)
plot_allo_cr

plot_allo_md <- ggplot(mv_scores_md, aes(x = Csize, y = V1, color = life_cycle)) +
  geom_point() +
  labs(x = "Csize", y = "mv scores") +
  #theme_minimal() +
  geom_line(aes(y = fitted_values, color = life_cycle),
            linewidth = 1) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "metamorphic" = "#009900", "undergoing metamorphosis" = "#00CCFF")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid = element_line(color = "gray90"),
        axis.line = element_line(color = "black", size = 0.5),
        aspect.ratio = 1)
plot_allo_md

{
ggsave("/Plot_allometry_linear.png",
       plot = plot_allo_linear, width = 8, height = 6, dpi = 300)
ggsave("/Plot_allometry_cr.png",
       plot = plot_allo_cr, width = 8, height = 6, dpi = 300) 
ggsave("/Plot_allometry_md.png",
       plot = plot_allo_md, width = 8, height = 6, dpi = 300) 
}


# Contrast coding
#_______________________________________________________________________________

# IMPORTANT: check the order of the coefficient to be sure to look at the interaction between the group and size
df_linear$life_cycle <- as.factor(df_linear$life_cycle)
df_Csize_cr$life_cycle <- as.factor(df_Csize_cr$life_cycle)
df_Csize_md$life_cycle <- as.factor(df_Csize_md$life_cycle)
contrasts(df_linear$life_cycle)
contrasts(df_Csize_cr$life_cycle)
contrasts(df_Csize_md$life_cycle)

pairwise.contrasts(fit_linear_allo, term = "life_cycle")
pairwise.contrasts(fit_cr_allo, term = "life_cycle")
pairwise.contrasts(fit_md_allo, term = "life_cycle")

# allows to test hypothesis of allometry in the different groups
allom_m <- matrix(c(0, 0, 0, 1, 0, 0),nrow = 1)
allom_p <- matrix(c(0, 0, 0, 1, 1, 0),nrow = 1)
allom_u <- matrix(c(0, 0, 0, 1, 0, 1),nrow = 1)

contrast_list <- list(
  allom_m = allom_m,
  allom_p = allom_p,
  allom_u = allom_u
)
names_contrast <- names(contrast_list)

### function to test allometry
test_allometry <- function(model, contrast) {
  manova.gls(model,
             test = "Pillai",
             nperm = 999,
             L = contrast,
             seed = 123,
             verbose = TRUE,
             type = "II")
}

### run analysis
allometry_results_linear <- list()
allometry_results_cr <- list()
allometry_results_md <- list()

for (j in 1:length(contrast_list)) {
  result <- test_allometry(fit_linear_allo, contrast_list[[j]])
  
  # store the results
  allometry_results_linear[[paste(names_contrast[j], sep = "_")]] <- result
}
for (j in 1:length(contrast_list)) {
  result <- test_allometry(fit_cr_allo, contrast_list[[j]])

  # store the results
  allometry_results_cr[[paste(names_contrast[j], sep = "_")]] <- result
}
for (j in 1:length(contrast_list)) {
  result <- test_allometry(fit_md_allo, contrast_list[[j]])

  # store the results
  allometry_results_md[[paste(names_contrast[j], sep = "_")]] <- result
}


### test for effect size

# create a list to store the results of the analyses
names_allometry_results_linear <- names(allometry_results_linear)
names_allometry_results_cr <- names(allometry_results_cr)
names_allometry_results_md <- names(allometry_results_md)
effectsize_allo_linear <- list()
effectsize_allo_cr <- list()
effectsize_allo_md <- list()

# loop to compute the analyses
for (i in 1:length(allometry_results_linear)) {
  model <- allometry_results_linear[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_allo_linear[[names_allometry_results_linear[i]]] <- result
}
for (i in 1:length(allometry_results_cr)) {
  model <- allometry_results_cr[[i]]
  result <- effectsize(model)

  # store the results
  effectsize_allo_cr[[names_allometry_results_cr[i]]] <- result
}
for (i in 1:length(allometry_results_md)) {
  model <- allometry_results_md[[i]]
  result <- effectsize(model)

  # store the results
  effectsize_allo_md[[names_allometry_results_md[i]]] <- result
}

### extract the results

allometry_results_linear <- data.frame(
  Model = names(allometry_results_linear),
  Test_stat = sapply(allometry_results_linear, function(x) x$stat),
  Effect_size = sapply(effectsize_allo_linear, function(x) x$effect),
  Pvalue = sapply(allometry_results_linear, function(x) x$pvalue),
  row.names = NULL
)
allometry_results_cr <- data.frame(
  Model = names(allometry_results_cr),
  Test_stat = sapply(allometry_results_cr, function(x) x$stat),
  Effect_size = sapply(effectsize_allo_cr, function(x) x$effect),
  Pvalue = sapply(allometry_results_cr, function(x) x$pvalue),
  row.names = NULL
)
allometry_results_md <- data.frame(
  Model = names(allometry_results_md),
  Test_stat = sapply(allometry_results_md, function(x) x$stat),
  Effect_size = sapply(effectsize_allo_md, function(x) x$effect),
  Pvalue = sapply(allometry_results_md, function(x) x$pvalue),
  row.names = NULL
)

{
write.csv(allometry_results_linear, file = "/results_allometry_linear.csv",
          row.names = FALSE)
write.csv(allometry_results_cr, file = "/results_allometry_cr.csv",
          row.names = FALSE)
write.csv(allometry_results_md, file = "/results_allometry_md.csv",
          row.names = FALSE)
}


# Allometric slope sifferences scross the sife cycle categories

diff_allom_p_m <- matrix(
  c(0, 0, 0, 0, 1, 0),
  nrow = 1
)

diff_allom_u_m <- matrix(
  c(0, 0, 0, 0, 0, 1),
  nrow = 1
)

diff_allom_u_p <- matrix(
  c(0, 0, 0, 0, -1, 1),
  nrow = 1
)

contrast_list <- list(
  diff_allom_p_m = diff_allom_p_m,
  diff_allom_u_m = diff_allom_u_m,
  diff_allom_u_p = diff_allom_u_p
)

test_diff_allometry <- function(model, contrast) {
  manova.gls(model,
             test = "Pillai",
             nperm = 999,
             L = contrast,
             seed = 123,
             verbose = TRUE)
}

diff_allometry_results_linear <- list()
diff_allometry_results_cr <- list()
diff_allometry_results_md <- list()

for (j in 1:length(contrast_list)) {
  result <- test_diff_allometry(fit_linear_allo, contrast_list[[j]])
  
  diff_allometry_results_linear[[names(contrast_list)[j]]] <- result
}
for (j in 1:length(contrast_list)) {
  result <- test_diff_allometry(fit_cr_allo, contrast_list[[j]])
  
  diff_allometry_results_cr[[names(contrast_list)[j]]] <- result
}
for (j in 1:length(contrast_list)) {
  result <- test_diff_allometry(fit_md_allo, contrast_list[[j]])
  
  diff_allometry_results_md[[names(contrast_list)[j]]] <- result
}

names_diff_allometry_results_linear <- names(diff_allometry_results_linear)
names_diff_allometry_results_cr <- names(diff_allometry_results_cr)
names_diff_allometry_results_md <- names(diff_allometry_results_md)
effectsize_diff_allo_linear <- list()
effectsize_diff_allo_cr <- list()
effectsize_diff_allo_md <- list()

# loop to compute the analyses
for (i in 1:length(diff_allometry_results_linear)) {
  model <- diff_allometry_results_linear[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_diff_allo_linear[[names_diff_allometry_results_linear[i]]] <- result
}
for (i in 1:length(diff_allometry_results_cr)) {
  model <- diff_allometry_results_cr[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_diff_allo_cr[[names_diff_allometry_results_cr[i]]] <- result
}
for (i in 1:length(diff_allometry_results_md)) {
  model <- diff_allometry_results_md[[i]]
  result <- effectsize(model)
  
  # store the results
  effectsize_diff_allo_md[[names_diff_allometry_results_md[i]]] <- result
}

diff_allometry_results_linear <- data.frame(
  Model = names(diff_allometry_results_linear),
  Test_stat = sapply(diff_allometry_results_linear, function(x) x$stat),
  Effect_size = sapply(effectsize_diff_allo_linear, function(x) x$effect),
  Pvalue = sapply(diff_allometry_results_linear, function(x) x$pvalue),
  row.names = NULL
)
diff_allometry_results_cr <- data.frame(
  Model = names(diff_allometry_results_cr),
  Test_stat = sapply(diff_allometry_results_cr, function(x) x$stat),
  Effect_size = sapply(effectsize_diff_allo_cr, function(x) x$effect),
  Pvalue = sapply(diff_allometry_results_cr, function(x) x$pvalue),
  row.names = NULL
)
diff_allometry_results_md <- data.frame(
  Model = names(diff_allometry_results_md),
  Test_stat = sapply(diff_allometry_results_md, function(x) x$stat),
  Effect_size = sapply(effectsize_diff_allo_md, function(x) x$effect),
  Pvalue = sapply(diff_allometry_results_md, function(x) x$pvalue),
  row.names = NULL
)

{
write.csv(diff_allometry_results_linear, file = "/diff_allometry_results_linear.csv",
          row.names = FALSE) 
write.csv(diff_allometry_results_cr, file = "/diff_allometry_results_cr.csv",
          row.names = FALSE) 
write.csv(diff_allometry_results_md, file = "/diff_allometry_results_md.csv",
          row.names = FALSE)
}

###############################################################################
## ------------------------------------------------------------------------- ##
##          Testing for difference in integration and modularity             ##
##                 between paedomorphic and metamorphic                      ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(geomorph)

# Body
#_______________________________________________________________________________
df_linear <-read.csv("Data/Linear_measurements/df_linear_mp.csv", 
                     sep=",",
                     row.names = 1)


# partition for functional hypothesis
partition_linear <- c("loco", "loco",
                      "hd", "hd", "hd",
                      "hd",
                      "loco", "loco", "loco", "loco",
                      "loco", "loco", "loco", "loco",
                      "loco"
)

# loco = locomotion
# hd = head --> prey capture

Y_linear <- as.matrix(df_linear[, c(7:11, 21:30)]) # You need at least 2 measurements per module, so remove the tail.
group <- df_linear$life_cycle
names(group) <- row.names(df_linear)

sp_m <- which(group == "metamorphic")
sp_p <- which(group == "paedomorphic")

min_size <- min(length(sp_m), length(sp_p))

# Perform bootstrap resampling
n_bootstrap <- 100  # Define the number of bootstrap samples
results_m <- vector("list", n_bootstrap)
results_p <- vector("list", n_bootstrap)
set.seed(123)  # For reproducibility

for (i in 1:n_bootstrap) {
  # Sample with replacement to balance the groups
  boot_m <- sample(sp_m, min_size, replace = TRUE)
  boot_p <- sample(sp_p, min_size, replace = TRUE)
  
  # Run the integration test on the resampled data
  int_m <- integration.test(A = Y_linear[boot_m, ],
                            partition.gp = NULL,
                            iter = 999,
                            seed = "random",
                            print.progress = T)
  int_p <- integration.test(A = Y_linear[boot_p, ],
                            partition.gp = NULL,
                            iter = 999,
                            seed = "random",
                            print.progress = T)
  
  # Store the results
  results_m[[i]] <- int_m
  results_p[[i]] <- int_p
}


# Example to access the first result of the bootstrap
summary(results_m[[1]])
summary(results_p[[1]])

### plot distributions of points per group

# create list of results

int_results <- list(int_m = int_m,
                    int_p = int_p)

# save results
int_rslts <- list(int_m = summary(int_m),
                  int_p = summary(int_p))

table_integration_linear <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  r.pls = c(int_m$r.pls,
            int_p$r.pls),
  Effect_size = c(int_m$Z,
                  int_p$Z),
  P_value = c(int_m$P.value,
              int_p$P.value)
)
table_integration_linear
# life_cycle     r.pls Effect_size P_value
# metamorphic   metamorphic 0.8430562    4.648706   0.001
# paedomorphic paedomorphic 0.7826485    3.799484   0.001

write.csv(table_integration_linear,
          "/results_integration_linear.csv",
          row.names = FALSE)


compare_integration <- compare.pls(int_m, 
                                   int_p)
compare_integration
# Effect sizes
# 
# int_m    int_p 
# 4.648706 3.799484 
# 
# Effect sizes for pairwise differences in PLS effect size
# 
# int_m     int_p
# int_m 0.0000000 0.4565851
# int_p 0.4565851 0.0000000
# 
# P-values
# 
# int_m     int_p
# int_m 1.0000000 0.6479693
# int_p 0.6479693 1.0000000


##### Modularity  
{
  modul_m <- modularity.test(A = Y_linear[boot_m, ],
                             partition.gp = partition_linear,
                             CI = TRUE,
                             iter = 999,
                             seed = "random",
                             print.progress = TRUE)
  
  modul_p <- modularity.test(A = Y_linear[boot_p, ],
                             partition.gp = partition_linear,
                             CI = TRUE,
                             iter = 999,
                             seed = "random",
                             print.progress = TRUE)
  
}


# create list of results
modul_results <- list(modul_m = modul_m,
                      modul_p = modul_p)

# save results
modul_rslts <- list(modul_m = summary(modul_m),
                    modul_p = summary(modul_p))

table_modularity_linear <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  cr = c(modul_m$CR,
         modul_p$CR),
  effect_size = c(modul_m$Z,
                  modul_p$Z),
  p_value = c(modul_m$P.value,
              modul_p$P.value)
)
table_modularity_linear
# life_cycle        cr effect_size p_value
# metamorphic   metamorphic 0.9896830   -1.992008  0.0230
# paedomorphic paedomorphic 0.9303156   -2.269787  0.0125

write.csv(table_modularity_linear,
          "/results_modularity_linear.csv",
          row.names = FALSE)


compare_modularity <- compare.CR(modul_m, 
                                 modul_p)
compare_modularity
# NOTE: more negative effects represent stronger modular signal! 
#   
#   
# Effect sizes
# 
# No_Modules    modul_m    modul_p 
# 0.000000  -1.992008  -2.269787 
# 
# Effect sizes for pairwise differences in CR effect size
# 
# No_Modules   modul_m   modul_p
# No_Modules   0.000000 1.9920077 2.2697871
# modul_m      1.992008 0.0000000 0.5314757
# modul_p      2.269787 0.5314757 0.0000000
# 
# P-values
# 
# No_Modules    modul_m   modul_p
# No_Modules 1.00000000 0.04637021 0.0232205
# modul_m    0.04637021 1.00000000 0.5950892
# modul_p    0.02322050 0.59508921 1.0000000

# Mandible
#_______________________________________________________________________________
load(file = "Data/Landmarks/pts/Y.gpa_md_mp.Rdata")
Y_md_mp <- Y.gpa_md_mp$rotated

ident_md <- read.csv("Data/Landmarks/ident/ident_md_mp.csv",
                     sep=",",header = T, row.names = 1)
table(ident_md$life_cycle)

specimens <- dimnames(Y_md_mp)[[3]]
ident_md_ordered <- ident_md[specimens, ]
sum(is.na(ident_md_ordered$life_cycle)) # must be 0
is_meta  <- ident_md_ordered$life_cycle == "metamorphic"
is_paedo <- ident_md_ordered$life_cycle == "paedomorphic"
Y_md_m  <- Y_md_mp[, , is_meta]
Y_md_p <- Y_md_mp[, , is_paedo]
all(specimens %in% rownames(ident_md))
all(rownames(ident_md) %in% specimens)

# Create partition for each bone 
partition_md <- c("d", "d", "d",
                  "pa", "pa",
                  "d", "d", "d", "d", "d", "d", "d", "d", "d", "d","d", "d", "d", "d", "d", "d", "d", "d", "d", "d","d", "d", "d", "d", "d", "d", "d", "d", "d", "d","d", "d", "d", "d", "d", "d", "d", "d", "d", "d",
                  "pa", "pa", "pa","pa", "pa", "pa","pa", "pa", "pa", "pa","pa", "pa", "pa","pa", "pa", "pa","pa", "pa", "pa", "pa"
)

# d = dentary
# pa = prearticular

{
  spm <- dimnames(Y_md_m)[[3]]
  spp <- dimnames(Y_md_p)[[3]]
}

# Determine the minimum group size
min_size <- min(length(spm), length(spp))
# Perform bootstrap resampling
n_bootstrap <- 100  # Define the number of bootstrap samples
results_m <- vector("list", n_bootstrap)
results_p <- vector("list", n_bootstrap)
set.seed(123)  # For reproducibility

# ##### Integration measure  

for (i in 1:n_bootstrap) {
  
  # 1️⃣ Bootstrap individus
  boot_m <- sample(spm, min_size, replace = TRUE)
  boot_p <- sample(spp, min_size, replace = TRUE)
  
  # 2️⃣ Sous-ensembles morpho
  Y_m_boot <- Y_md_m[, , boot_m]
  Y_p_boot <- Y_md_p[, , boot_p]
  
  # 6️⃣ Integration
  int_m <- integration.test(
    A = Y_m_boot,
    partition.gp = partition_md,
    iter = 999,
    seed = "random",
    print.progress = T
  )
  
  int_p <- integration.test(
    A = Y_p_boot,
    partition.gp = partition_md,
    iter = 999,
    seed = "random",
    print.progress = T
  )
  
  # 7️⃣ Stockage
  results_m[[i]] <- int_m
  results_p[[i]] <- int_p
}


# Example to access the first result of the bootstrap
summary(results_m[[1]])
summary(results_p[[1]])

### plot distributions of points per group

# create list of results

int_results <- list(int_m = int_m,
                    int_p = int_p)

# save results
int_rslts <- list(int_m = summary(int_m),
                  int_p = summary(int_p))

table_integration_md <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  r.pls = c(int_m$r.pls,
            int_p$r.pls),
  Effect_size = c(int_m$Z,
                  int_p$Z),
  P_value = c(int_m$P.value,
              int_p$P.value)
)
table_integration_md
# life_cycle     r.pls Effect_size P_value
# metamorphic   metamorphic 0.9063347    5.393192   0.001
# paedomorphic paedomorphic 0.9477788    6.275334   0.001

write.csv(table_integration_md,
          "/results_integration_md.csv",
          row.names = FALSE)


compare_integration <- compare.pls(int_m, 
                                   int_p)
compare_integration
# Effect sizes
# 
# int_m    int_p 
# 5.393192 6.275334 
# 
# Effect sizes for pairwise differences in PLS effect size
# 
# int_m    int_p
# int_m 0.000000 1.153188
# int_p 1.153188 0.000000
# 
# P-values
# 
# int_m     int_p
# int_m 1.0000000 0.2488332
# int_p 0.2488332 1.0000000


##### Modularity  
{
  modul_m <- modularity.test(A = Y_m_boot,
                             partition.gp = partition_md,
                             CI = TRUE,
                             iter = 999,
                             seed = "random",
                             print.progress = TRUE)
  
  modul_p <- modularity.test(A = Y_p_boot,
                             partition.gp = partition_md,
                             CI = TRUE,
                             iter = 999,
                             seed = "random",
                             print.progress = TRUE)
  
}


# create list of results
modul_results <- list(modul_m = modul_m,
                      modul_p = modul_p)

# save results
modul_rslts <- list(modul_m = summary(modul_m),
                    modul_p = summary(modul_p))

table_modularity_md <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  cr = c(modul_m$CR,
         modul_p$CR),
  effect_size = c(modul_m$Z,
                  modul_p$Z),
  p_value = c(modul_m$P.value,
              modul_p$P.value)
)
table_modularity_md
# life_cycle        cr effect_size p_value
# metamorphic   metamorphic 0.8304171   -4.104008   0.001
# paedomorphic paedomorphic 0.8478613   -4.576659   0.001

write.csv(table_modularity_md,
          "/results_modularity_md.csv",
          row.names = FALSE)

compare_modularity <- compare.CR(modul_m, 
                                 modul_p)
compare_modularity
# NOTE: more negative effects represent stronger modular signal! 
#   
#   
#   Effect sizes
# 
# No_Modules    modul_m    modul_p 
# 0.000000  -4.104008  -4.576659 
# 
# Effect sizes for pairwise differences in CR effect size
# 
# No_Modules   modul_m   modul_p
# No_Modules   0.000000 4.1040078 4.5766587
# modul_m      4.104008 0.0000000 0.6887441
# modul_p      4.576659 0.6887441 0.0000000
# 
# P-values
# 
# No_Modules      modul_m      modul_p
# No_Modules 1.000000e+00 4.060538e-05 4.724614e-06
# modul_m    4.060538e-05 1.000000e+00 4.909843e-01
# modul_p    4.724614e-06 4.909843e-01 1.000000e+00


# Cranium
#_______________________________________________________________________________
load(file = "Data/Landmarks/pts/Y.gpa_cr_mp.Rdata")
Y_cr_mp <- Y.gpa_cr_mp$rotated

ident_cr <- read.csv("Data/Landmarks/ident/ident_cr_mp.csv",
                     sep=",",header = T, row.names = 1)
table(ident_cr$life_cycle)

specimens <- dimnames(Y_cr_mp)[[3]]
ident_cr_ordered <- ident_cr[specimens, ]
sum(is.na(ident_cr_ordered$life_cycle)) 
is_meta  <- ident_cr_ordered$life_cycle == "metamorphic"
is_paedo <- ident_cr_ordered$life_cycle == "paedomorphic"
Y_cr_m  <- Y_cr_mp[, , is_meta]
Y_cr_p <- Y_cr_mp[, , is_paedo]
all(specimens %in% rownames(ident_cr))
all(rownames(ident_cr) %in% specimens)

# partition for developmental hypothesis
partition_cr <- c("early", "mid", "mid",
                  "early", "early", "early",
                  "late", "late", "late", "late", "late",
                  "late", "late", "late", "late",
                  "late", "late", "late", "late",
                  "mid", "mid", "mid", "mid", "mid", "mid", "mid", 
                  "mid", "mid", "mid", "mid",
                  "met", "met", "met", "met", "met",
                  "mid", "mid", "mid", "mid", "mid", "mid", "mid", "mid", 
                  "mid", "mid", "mid",
                  "mid", "mid", "mid", "mid",
                  "mid", "mid", "mid", "mid", "mid", "mid", "mid", "mid", "mid", "mid", "mid",
                  "early", "early", "early", "early", "early", 
                  "late", "late", "late", "late", "late", "late", 
                  "met", "met", "met") 

{
  spm <- dimnames(Y_cr_m)[[3]]
  spp <- dimnames(Y_cr_p)[[3]]
}

# Determine the minimum group size
min_size <- min(length(spm), length(spp))
# Perform bootstrap resampling
n_bootstrap <- 100  # Define the number of bootstrap samples
results_m <- vector("list", n_bootstrap)
results_p <- vector("list", n_bootstrap)
set.seed(123)  # For reproducibility

# ##### Integration measure  
library(geomorph)

for (i in 1:n_bootstrap) {
  
  # 1️⃣ Bootstrap individus
  boot_m <- sample(spm, min_size, replace = TRUE)
  boot_p <- sample(spp, min_size, replace = TRUE)
  
  # 2️⃣ Sous-ensembles morpho
  Y_m_boot <- Y_cr_m[, , boot_m]
  Y_p_boot <- Y_cr_p[, , boot_p]
  
  # 6️⃣ Integration
  int_m <- integration.test(
    A = Y_m_boot,
    partition.gp = partition_cr,
    iter = 999,
    seed = "random",
    print.progress = T
  )
  
  int_p <- integration.test(
    A = Y_p_boot,
    partition.gp = partition_cr,
    iter = 999,
    seed = "random",
    print.progress = T
  )
  
  # 7️⃣ Stockage
  results_m[[i]] <- int_m
  results_p[[i]] <- int_p
}


# Example to access the first result of the bootstrap
summary(results_m[[1]])
summary(results_p[[1]])

### plot distributions of points per group

# create list of results

int_results <- list(int_m = int_m,
                    int_p = int_p)

# save results
int_rslts <- list(int_m = summary(int_m),
                  int_p = summary(int_p))

table_integration_cr <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  r.pls = c(int_m$r.pls,
            int_p$r.pls),
  Effect_size = c(int_m$Z,
                  int_p$Z),
  P_value = c(int_m$P.value,
              int_p$P.value)
)
table_integration_cr
# life_cycle     r.pls Effect_size P_value
# metamorphic   metamorphic 0.8035619    8.767140   0.001
# paedomorphic paedomorphic 0.8744739    6.826211   0.001

write.csv(table_integration_cr,
          "/results_integration_cr.csv",
          row.names = FALSE)

compare_integration <- compare.pls(int_m, 
                                   int_p)
compare_integration
# Effect sizes
# 
# int_m    int_p 
# 5.084311 5.268103 
# 
# Effect sizes for pairwise differences in PLS effect size
# 
# int_m     int_p
# int_m 0.0000000 0.9561079
# int_p 0.9561079 0.0000000
# 
# P-values
# 
# int_m     int_p
# int_m 1.0000000 0.3390177
# int_p 0.3390177 1.0000000


{
  modul_m <- modularity.test(A = Y_m_boot,
                              partition.gp = partition_cr,
                              CI = TRUE,
                              iter = 999,
                              seed = "random",
                              print.progress = TRUE)
  
  modul_p <- modularity.test(A = Y_p_boot,
                              partition.gp = partition_cr,
                              CI = TRUE,
                              iter = 999,
                              seed = "random",
                              print.progress = TRUE)
  
}


# create list of results
modul_results <- list(modul_m = modul_m,
                      modul_p = modul_p)

# save results
modul_rslts <- list(modul_m = summary(modul_m),
                    modul_p = summary(modul_p))

table_modularity_cr <- data.frame(
  row.names =  c("metamorphic",
                 "paedomorphic"),
  life_cycle = c("metamorphic",
                 "paedomorphic"),
  cr = c(modul_m$CR,
         modul_p$CR),
  effect_size = c(modul_m$Z,
                  modul_p$Z),
  p_value = c(modul_m$P.value,
              modul_p$P.value)
)
table_modularity_cr
# life_cycle        cr effect_size p_value
# metamorphic   metamorphic 0.9166347   -1.526055   0.072
# paedomorphic paedomorphic 0.9169588   -2.502965   0.008

write.csv(table_modularity_cr,
          "/results_modularity_cr.csv",
          row.names = FALSE)

compare_modularity <- compare.CR(modul_m, 
                                 modul_p)
compare_modularity
# NOTE: more negative effects represent stronger modular signal! 
#   
#   
#   Effect sizes
# 
# No_Modules    modul_m    modul_p 
# 0.000000  -1.526055  -2.502965 
# 
# Effect sizes for pairwise differences in CR effect size
# 
# No_Modules   modul_m   modul_p
# No_Modules   0.000000 1.5260552 2.5029646
# modul_m      1.526055 0.0000000 0.8479685
# modul_p      2.502965 0.8479685 0.0000000
# 
# P-values
# 
# No_Modules   modul_m    modul_p
# No_Modules 1.00000000 0.1269961 0.01231579
# modul_m    0.12699612 1.0000000 0.39645551
# modul_p    0.01231579 0.3964555 1.00000000


###############################################################################
## ------------------------------------------------------------------------- ##
##         Testing for difference in morphology and shape diversity          ##
##                    between paedomorphic and metamorphic                   ##
##                       using disparity analyses                            ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(dispRity)

# Data import
#_______________________________________________________________________________
{
  df_linear <-read.csv("Data/Linear_measurements/df_linear_mp.csv", 
                     sep=",",
                     row.names = 1)

load(file = "Data/Landmarks/pts/Y.gpa_cr_mp.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_md_mp.Rdata")

load(file = "Data/Landmarks/pts/Y.gpa_cr_pmx.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_max.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_psph.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_nas.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_prf.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_fr.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_par.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_vom.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_osph.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_oticocc.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_sq.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_qd.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_cr_pt.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_md_d.Rdata")
load(file = "Data/Landmarks/pts/Y.gpa_md_pa.Rdata")
}

# List all Y.gpa objects
gpa_objs <- ls(pattern = "^Y\\.gpa_")
# Loop through all objects
for (obj in gpa_objs) {
  gpa <- get(obj)
  # Extract the aligned coordinates
  Y_ <- gpa$rotated
  # Conversion to a 2D matrix
  Y_2d <- two.d.array(Y_)
  # Output name: Y_...
  out_name <- sub("^Y\\.gpa_", "Y_", obj)
  assign(out_name, Y_2d, envir = .GlobalEnv)
}


ident_cr <- read.csv("Data/Landmarks/Ident/ident_cr_mp.csv",
                     sep=",",header = T, row.names = 1)
ident_md <- read.csv("Data/Landmarks/Ident/ident_md_mp.csv",
                     sep=",",header = T, row.names = 1)

Y_cr <- Y_cr_mp[match(rownames(ident_cr),
                             rownames(Y_cr_mp)), ]
all(rownames(Y_cr) == rownames(ident_cr))

Y_md <- Y_md_mp[match(rownames(ident_md),
                             rownames(Y_md_mp)), ]
all(rownames(Y_md) == rownames(ident_md))

{
all(rownames(Y_cr_fr) == rownames(ident_cr))
all(rownames(Y_cr_max) == rownames(ident_cr))
all(rownames(Y_cr_nas) == rownames(ident_cr))
all(rownames(Y_cr_osph) == rownames(ident_cr))
all(rownames(Y_cr_oticocc) == rownames(ident_cr))
all(rownames(Y_cr_par) == rownames(ident_cr))
all(rownames(Y_cr_pmx) == rownames(ident_cr))
all(rownames(Y_cr_prf) == rownames(ident_cr))
all(rownames(Y_cr_psph) == rownames(ident_cr))
all(rownames(Y_cr_pt) == rownames(ident_cr))
all(rownames(Y_cr_qd) == rownames(ident_cr))
all(rownames(Y_cr_sq) == rownames(ident_cr))
all(rownames(Y_cr_vom) == rownames(ident_cr))
}

Y_md_d <- Y_md_d[match(rownames(ident_md),
                      rownames(Y_md_d)), ]
all(rownames(Y_md_d) == rownames(ident_md))

Y_md_pa <- Y_md_pa[match(rownames(ident_md),
                       rownames(Y_md_pa)), ]
all(rownames(Y_md_pa) == rownames(ident_md))

group_factor_lin <- as.factor(df_linear$life_cycle)
group_factor_cr <- as.factor(ident_cr$life_cycle)
group_factor_md <- as.factor(ident_md$life_cycle)

# Disparity tests
#_______________________________________________________________________________

# For each structure

disparity_result_linear <- dispRity.per.group(
  data  =df_linear[, c(7:11, 21:30)],
  group = group_factor_lin,
  metric = c(sum, variances)
)

disparity_result_cr <- dispRity.per.group(
  data  = Y_cr,
  group = group_factor_cr,
  metric = c(sum, variances)
)

disparity_result_md <- dispRity.per.group(
  data  = Y_md,
  group = group_factor_md,
  metric = c(sum, variances)
)

summary(disparity_result_linear)
summary(disparity_result_cr)
summary(disparity_result_md)

{
save(disparity_result_linear, file = "/disparity_result_linear.Rdata")
save(disparity_result_cr, file = "/disparity_result_cr.Rdata")
save(disparity_result_md, file = "/disparity_result_md.Rdata")
}

# For each part of the body

disparity_result_head <- dispRity.per.group(
  data  = df_linear[, c(9:11, 21)],
  group = group_factor_lin,
  metric = c(sum, variances)
)
save(disparity_result_head, file = "/disparity_result_head.Rdata")
summary(disparity_result_head)

disparity_result_trunk <- dispRity.per.group(
  data  =df_linear[, 7:8],
  group = group_factor_lin,
  metric = c(sum, variances)
)
save(disparity_result_trunk, file = "/disparity_result_trunk.Rdata")
summary(disparity_result_trunk)

disparity_result_legs <- dispRity.per.group(
  data  = df_linear[, c(22:29)],
  group = group_factor_lin,
  metric = c(sum, variances)
)
save(disparity_result_legs, file = "/disparity_result_legs.Rdata")
summary(disparity_result_legs)

# For each bone

disparity_fr <- dispRity.per.group(
  data  = Y_cr_fr,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_fr, file = "/disparity_fr.Rdata")
summary(disparity_fr)

disparity_max <- dispRity.per.group(
  data  = Y_cr_max,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_max, file = "/disparity_max.Rdata")
summary(disparity_max)

disparity_nas <- dispRity.per.group(
  data  = Y_cr_nas,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_nas, file = "/disparity_nas.Rdata")
summary(disparity_nas)

disparity_osph <- dispRity.per.group(
  data  = Y_cr_osph,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_osph, file = "/disparity_osph.Rdata")
summary(disparity_osph)

disparity_oticocc <- dispRity.per.group(
  data  = Y_cr_oticocc,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_oticocc, file = "/disparity_oticocc.Rdata")
summary(disparity_oticocc)

disparity_par <- dispRity.per.group(
  data  = Y_cr_par,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_par, file = "/disparity_par.Rdata")
summary(disparity_par)

disparity_pmx <- dispRity.per.group(
  data  = Y_cr_pmx,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_pmx, file = "/disparity_pmx.Rdata")
summary(disparity_pmx)

disparity_prf <- dispRity.per.group(
  data  = Y_cr_prf,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_prf, file = "/disparity_prf.Rdata")
summary(disparity_prf)

disparity_psph <- dispRity.per.group(
  data  = Y_cr_psph,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_psph, file = "/disparity_psph.Rdata")
summary(disparity_psph)

disparity_pt <- dispRity.per.group(
  data  = Y_cr_pt,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_pt, file = "/disparity_pt.Rdata")
summary(disparity_pt)

disparity_qd <- dispRity.per.group(
  data  = Y_cr_qd,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_qd, file = "/disparity_qd.Rdata")
summary(disparity_qd)

disparity_sq <- dispRity.per.group(
  data  = Y_cr_sq,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_sq, file = "/disparity_sq.Rdata")
summary(disparity_sq)

disparity_vom <- dispRity.per.group(
  data  = Y_cr_vom,
  group = group_factor_cr,
  metric = c(sum, variances)
)
save(disparity_vom, file = "/disparity_vom.Rdata")
summary(disparity_vom)

disparity_d <- dispRity.per.group(
  data  = Y_md_d,
  group = group_factor_md,
  metric = c(sum, variances)
)
save(disparity_d, file = "/disparity_d.Rdata")
summary(disparity_d)

disparity_pa <- dispRity.per.group(
  data  = Y_md_pa,
  group = group_factor_md,
  metric = c(sum, variances)
)
save(disparity_pa, file = "/disparity_pa.Rdata")
summary(disparity_pa)


## To test if the disparity is different between the two life cycle categories

test_disp_linear <- test.dispRity(disparity_result_linear,
                                  test = wilcox.test,
                                  comparisons = "pairwise")
# [[1]]
# statistic: W
# metamorphic : paedomorphic         1645
# 
# [[2]]
# p.value
# metamorphic : paedomorphic 2.478002e-16
save(test_disp_linear, file = "/test_disp_linear.Rdata")

test_disp_cr <- test.dispRity(disparity_result_cr,
                                  test = wilcox.test,
                                  comparisons = "pairwise")
# [[1]]
# statistic: W
# metamorphic : paedomorphic           12
# 
# [[2]]
# p.value
# metamorphic : paedomorphic 3.672819e-34
save(test_disp_cr, file = "/test_disp_cr.Rdata")

test_disp_md <- test.dispRity(disparity_result_md,
                              test = wilcox.test,
                              comparisons = "pairwise")
# [[1]]
# statistic: W
# metamorphic : paedomorphic           14
# 
# [[2]]
# p.value
# metamorphic : paedomorphic 3.899684e-34
save(test_disp_md, file = "/test_disp_md.Rdata")

test_disp_head <- test.dispRity(disparity_result_head,
                              test = wilcox.test,
                              comparisons = "pairwise")
save(test_disp_head, file = "/test_disp_head.Rdata")

test_disp_trunk <- test.dispRity(disparity_result_trunk,
                                test = wilcox.test,
                                comparisons = "pairwise")
save(test_disp_trunk, file = "/test_disp_trunk.Rdata")

test_disp_legs <- test.dispRity(disparity_result_legs,
                                 test = wilcox.test,
                                 comparisons = "pairwise")
save(test_disp_legs, file = "/test_disp_legs.Rdata")

test_disp_fr <- test.dispRity(disparity_fr,
                                  test = wilcox.test,
                                  comparisons = "pairwise")
save(test_disp_fr, file = "/test_disp_fr.Rdata")

test_disp_max <- test.dispRity(disparity_max,
                              test = wilcox.test,
                              comparisons = "pairwise")
save(test_disp_max, file = "/test_disp_max.Rdata")

test_disp_nas <- test.dispRity(disparity_nas,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_nas, file = "/test_disp_nas.Rdata")

test_disp_osph <- test.dispRity(disparity_osph,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_osph, file = "/test_disp_osph.Rdata")

test_disp_oticocc <- test.dispRity(disparity_oticocc,
                                test = wilcox.test,
                                comparisons = "pairwise")
save(test_disp_oticocc, file = "/test_disp_oticocc.Rdata")

test_disp_par <- test.dispRity(disparity_par,
                                   test = wilcox.test,
                                   comparisons = "pairwise")
save(test_disp_par, file = "/test_disp_par.Rdata")

test_disp_pmx <- test.dispRity(disparity_pmx,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_pmx, file = "/test_disp_pmx.Rdata")

test_disp_prf <- test.dispRity(disparity_prf,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_prf, file = "/test_disp_prf.Rdata")

test_disp_psph <- test.dispRity(disparity_psph,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_psph, file = "/test_disp_psph.Rdata")

test_disp_pt <- test.dispRity(disparity_pt,
                                test = wilcox.test,
                                comparisons = "pairwise")
save(test_disp_pt, file = "/test_disp_pt.Rdata")

test_disp_qd <- test.dispRity(disparity_qd,
                              test = wilcox.test,
                              comparisons = "pairwise")
save(test_disp_qd, file = "/test_disp_qd.Rdata")

test_disp_sq <- test.dispRity(disparity_sq,
                              test = wilcox.test,
                              comparisons = "pairwise")
save(test_disp_sq, file = "/test_disp_sq.Rdata")

test_disp_vom <- test.dispRity(disparity_vom,
                              test = wilcox.test,
                              comparisons = "pairwise")
save(test_disp_vom, file = "/test_disp_vom.Rdata")

test_disp_d <- test.dispRity(disparity_d,
                               test = wilcox.test,
                               comparisons = "pairwise")
save(test_disp_d, file = "/test_disp_d.Rdata")

test_disp_pa <- test.dispRity(disparity_pa,
                             test = wilcox.test,
                             comparisons = "pairwise")
save(test_disp_pa, file = "/test_disp_pa.Rdata")



# Simple plot
plot(disparity_result_linear,
     type = "box",
     observed = list(pch = 21, 
                     col = "black",
                     bg = "#e74c3c",
                     cex = 1.4))
plot(disparity_result_cr,
     type = "box",
     observed = list(pch = 21, 
                     col = "black",
                     bg = "#e74c3c",
                     cex = 1.4))
plot(disparity_result_md,
     type = "box",
     observed = list(pch = 21, 
                     col = "black",
                     bg = "#e74c3c",
                     cex = 1.4))


extract_disp <- function(disp_obj, structure_name) {
  
  disp <- get.disparity(disp_obj, observed = TRUE)
  disp <- as.data.frame(disp)
  disp <- as.data.frame(t(disp))
  
  disp$life_cycle <- c("metamorphic", "paedomorphic")
  disp$structure  <- structure_name
  
  colnames(disp)[1] <- "disparity_scores"
  
  disp
}

disp_objs <- ls(pattern = "^disparity_")

disp_list <- lapply(disp_objs, function(obj) {
  extract_disp(get(obj), structure_name = obj)
})

names(disp_list) <- disp_objs
disp_all <- do.call(rbind, disp_list)
rownames(disp_all) <- NULL


# List all objects named disparity_
disp_objs <- ls(pattern = "^disparity_")

# Apply the t-test to each one
test_list <- lapply(disp_objs, function(obj) {
  test.dispRity(get(obj), test = t.test)
})

names(test_list) <- disp_objs


# Plots
#_______________________________________________________________________________
disp_plot <- disp_all %>%
  # dplyr::filter(structure != "disparity_result_linear") %>%
  dplyr::mutate(
    region = dplyr::case_when(
      structure %in% c("disparity_result_linear", "disparity_result_head", "disparity_result_trunk", "disparity_result_legs") ~ "1Linear",
      structure %in% c("disparity_result_md", "disparity_d", "disparity_pa") ~ "3Mandible",
      TRUE                                              ~ "2Cranium"
    )
  )

library(ggplot2)
plot_disp_regions <- ggplot(
  disp_plot,
  aes(
    x     = region,
    y     = disparity_scores,
    color = life_cycle,
    shape = structure
  )
) +
  geom_point(
    size = 3,
   # position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(
    values = c(
      "metamorphic"  = "#009900",
      "paedomorphic" = "#6633FF"
    )
  ) +
  scale_shape_manual(
    values = c(
      "disparity_result_cr"     = 16,
      "disparity_result_linear" = 16,
      "disparity_result_md"     = 16
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_line(color = "gray90"),
    axis.line        = element_line(color = "black", size = 0.5),
    aspect.ratio     = 1.1
  ) +
  xlab("Anatomical region") +
  ylab("Disparity score")
plot_disp_regions

ggsave("/plot_disp_regions.png",
       plot = plot_disp_regions, width = 8, height = 6, dpi = 300) 

plot_disp_regions_zoom <- ggplot(
  disp_plot,
  aes(
    x     = region,
    y     = disparity_scores,
    color = life_cycle,
    shape = structure
  )
) +
  geom_point(
    size = 3,
    # position = position_dodge(width = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_color_manual(
    values = c(
      "metamorphic"  = "#009900",
      "paedomorphic" = "#6633FF"
    )
  ) +
  scale_shape_manual(
    values = c(
      "disparity_result_cr"     = 16,
      "disparity_result_linear" = 16,
      "disparity_result_md"     = 16
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_line(color = "gray90"),
    axis.line        = element_line(color = "black", size = 0.5),
    aspect.ratio     = 1.1
  ) +
  xlab("Anatomical region") +
  ylab("Disparity score")
plot_disp_regions_zoom
ggsave("/plot_disp_regions_zoom.png",
       plot = plot_disp_regions_zoom, width = 8, height = 6, dpi = 300) 

plot_disp_body <- ggplot(
  subset(disp_plot, region == "1Linear"& 
           structure != "disparity_result_linear"),
  aes(
    x     = structure,
    y     = disparity_scores,
    color = life_cycle
  )
) +
  geom_point(size = 3) +
  coord_cartesian(ylim = c(0, 0.04)) +
  scale_color_manual(
    values = c(
      "metamorphic"  = "#009900",
      "paedomorphic" = "#6633FF"
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_line(color = "gray90"),
    axis.line        = element_line(color = "black", size = 0.5),
    aspect.ratio     = 4,
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  xlab("Body structure") +
  ylab("Disparity score")
plot_disp_body

ggsave("/plot_disp_body.png",
       plot = plot_disp_body, width = 8, height = 6, dpi = 300)


plot_disp_skull <- ggplot(
  subset(disp_plot,
         (region == "3Mandible" & structure != "disparity_result_md") |
           (region == "2Cranium"  & structure != "disparity_result_cr")),
  aes(
    x     = structure,
    y     = disparity_scores,
    color = life_cycle
  )
) +
  geom_point(size = 3) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_color_manual(
    values = c(
      "metamorphic"  = "#009900",
      "paedomorphic" = "#6633FF"
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_line(color = "gray90"),
    axis.line        = element_line(color = "black", size = 0.5),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  xlab("Skull structure") +
  ylab("Disparity score")

plot_disp_skull
ggsave("/plot_disp_cr_and_md.png",
       plot = plot_disp_skull, width = 8, height = 6, dpi = 300) 


###############################################################################
## ------------------------------------------------------------------------- ##
##                   Testing for differences in body size                    ##
##                         depending on life cycles                          ##
#### ----------------------------------------------------------------------####
library(nlme)
library(car)
library(emmeans)
library(dplyr)

df_linear <-read.csv("Data/Linear_measurements/df_svl.csv", 
                     sep=",",
                     row.names = 1) # ordi pro

df_linear2 <-read.csv("Data/Linear_measurements/Linear_measurements_Ichthyosaura_alpestris_lc.csv", 
                     sep=";",
                     row.names = 1) # ordi pro

df_linear2$svl <- as.numeric(as.character(df_linear2$svl))
df_linear2 %>%
  group_by(life_cycle) %>%
  summarise(mean_svl = mean(svl, na.rm = TRUE),
            sd_svl   = sd(svl, na.rm = TRUE),
            n        = n())
# # A tibble: 3 × 4
# life_cycle               mean_svl sd_svl     n
# <chr>                       <dbl>  <dbl> <int>
#   1 metamorphic                  43.7   5.63    65
# 2 paedomorphic                 40.5   5.21   105
# 3 undergoing metamorphosis     40.3   4.07    27

# Comparison between life cycle
#_______________________________________________________________________________
fit_svl <- lm(
  svl ~ life_cycle + country,
  data = df_linear
)

summary(fit_svl)

res_anova_svl <- Anova(fit_svl, type=c("II"))
res_anova_svl
# Anova Table (Type II tests)
# 
# Response: svl
# Sum Sq  Df F value    Pr(>F)    
# life_cycle 0.04233   2  7.4935 0.0007351 ***
#   country    0.03661   2  6.4801 0.0018907 ** 
#   Residuals  0.54234 192                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


emm <- emmeans(fit_svl, ~ life_cycle)
res_svl_pairwise <- pairs(emm, adjust = "holm")
res_svl_pairwise
# contrast                                estimate     SE  df t.ratio p.value
# metamorphic - paedomorphic                0.0319 0.0084 192   3.797  0.0006
# metamorphic - undergoing metamorphosis    0.0283 0.0123 192   2.301  0.0449
# paedomorphic - undergoing metamorphosis  -0.0036 0.0116 192  -0.311  0.7564
# 
# Results are averaged over the levels of: ss_country 
# P value adjustment: holm method for 3 tests 

save(res_anova_svl, file = "/res_anova_svl.Rdata")
save(res_svl_pairwise, file = "/res_svl_pairwise.Rdata")

# Extractdata as tables
#_______________________________________________________________________________
table_anova_svl <- data.frame(
  term = rownames(res_anova_svl)[1:2],
  Sum_Sq = res_anova_svl$`Sum Sq`[1:2],
  Df = res_anova_svl$Df[1:2],
  F = res_anova_svl$`F value`[1:2],
  p_value = res_anova_svl$`Pr(>F)`[1:2]
)
table_anova_svl

table_pairwise_svl <- as.data.frame(summary(res_svl_pairwise))
table_pairwise_svl

write.csv(table_anova_svl,
          "/table_anova_svl.csv",
          row.names = FALSE)
write.csv(table_pairwise_svl,
          "/table_pairwise_svl.csv",
          row.names = FALSE)


###############################################################################
## ------------------------------------------------------------------------- ##
##                                  PCA                                      ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(ggplot2)
library(factoextra)
library(dplyr)
library(tidyr)
library(ggrepel)
library(geomorph)
library(Rvcg)
library(Morpho)


# Body
#_______________________________________________________________________________
df_linear <-read.csv("Data/Linear_measurements/df_linear_filtered.csv",
                     sep=",",
                     row.names = 1)

# PCA
pca_res_linear <- prcomp(df_linear[, c(7:11, 21:30)], center = TRUE, scale = TRUE)

df_pca_res_linear <- data.frame(pca_res_linear$x[,1:4])

df_pca_res_linear <- cbind(df_pca_res_linear, 
                           df_linear$sub.species, 
                           df_linear$country, 
                           df_linear$life_cycle)
colnames(df_pca_res_linear) <- c("PC1", "PC2", "PC3", "PC4", "subspecies", "country", "life_cycle")


write.csv(df_pca_res_linear,
          "/df_pca_res_linear.csv",
          row.names = TRUE) 



# Visualization of % variance and variable contributions

fviz_eig(pca_res_linear, addlabels = TRUE)
# Dim1 = 53.6 %
# Dim2 = 13.3 %
# Dim3 = 6.4 %
# Dim4 = 4.4 %

fviz_pca_var(pca_res_linear, axes = c(1,2), col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal(),
             repel = TRUE)

fviz_contrib(pca_res_linear, choice = "var", axes = c(1))
var_contrib <- fviz_contrib(pca_res_linear, choice = "var", axes = c(2))
fviz_contrib(pca_res_linear, choice = "var", axes = c(3))
fviz_contrib(pca_res_linear, choice = "var", axes = c(4))

ggsave("/Analyses/Results/PCA/var_contrib.png",
       plot = var_contrib, width = 8, height = 6, dpi = 300) 

loadings <- pca_res_linear$rotation
loadings_axis2 <- sort(loadings[,2], decreasing = TRUE)
loadings_axis2
# l_metatarsus  l_metacarpus    l_toe_hind    l_toe_fore       l_tibia l_lower_jaw_l     l_humerus      l_radius       l_femur l_interlimb_l          tail        head_w        head_l        body_h        head_h 
# 0.34699509    0.33172996    0.26188835    0.16540375    0.11664558    0.06981181    0.06812485    0.06616035    0.02471869   -0.01145585   -0.13153688   -0.15160564   -0.35738893   -0.46737584   -0.50883412 


vars <- c("head_h", "body_h", "head_l",
          "l_metacarpus", "r_metatarsus")

means_multi <- df_linear %>%
  select(life_cycle, all_of(vars)) %>%
  pivot_longer(cols = -life_cycle, names_to = "variable", values_to = "value") %>%
  group_by(life_cycle, variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

# Plots

plot_mean_measurements <- ggplot(means_multi, aes(x = variable, y = mean_value, fill = life_cycle)) +
  geom_col(position = position_dodge(width = 0.7)) +
  labs(y = "Mean", x = "Variables") +
  scale_fill_manual(
    name = "Life cycle",
    values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave("/plot_mean_measurements.png",
       plot = plot_mean_measurements, width = 8, height = 6, dpi = 300) 

# PCA plot
hulls_lin <- df_pca_res_linear %>%
  group_by(life_cycle) %>%
  slice(chull(PC1, PC2))
plot_PCA_ext <- ggplot(df_pca_res_linear) +
  aes(PC1, PC2, color = life_cycle, shape = country) +
  # geom_text_repel(aes(label = rownames(df_pca_res_linear)),
  #                 size = 3,
  #                 max.overlaps = Inf) +
  geom_polygon(
    data = hulls_lin,
    aes(fill = life_cycle, group = interaction(life_cycle)),
    alpha = 0.2,
    color = "gray80",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  ) +
  scale_fill_manual(
    name = "Life cycle",
    values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  ) +
  scale_shape_manual(values = c("Italy" = 17, "Bosnia" = 1, "Switzerland" = 8)) +
  geom_point(size = 2) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "gray90"),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)) +
  xlab("PC1: 53.6%") +
  ylab("PC2: 13.3%") +
  coord_fixed(ratio = 1)
plot_PCA_ext


ggsave("/plot_PCA_linear.png",
       plot = plot_PCA_ext, width = 8, height = 6, dpi = 300)


# Cranium
#_______________________________________________________________________________
load(file = "Data/Landmarks/pts/Y.gpa_cr.Rdata")
ident_cr <- read.csv("Data/Landmarks/Ident/ident_&_Csize_cr.csv",
                     sep=";",header = T, row.names = 1)
Y_cr <- Y.gpa_cr$rotated

# PCA
PCA_cr <- gm.prcomp(Y_cr)

# % Variance Explained
summary(PCA_cr)
# Comp1 = 71.7 %
# Comp2 = 3.3 %


Pcscores_cr <- data.frame(PCA_cr$x[,1:4])

Pcscores_cr <- merge(Pcscores_cr,
                       ident_cr,
                       by = 'row.names',
                       all = T)

write.csv(Pcscores_cr,
          "/Pcscores_cr.csv",
          row.names = TRUE)

# PCA plot
hulls_lin <- Pcscores_cr %>%
  group_by(life_cycle) %>%
  slice(chull(Comp1, Comp2))
plot_PCA_cr <- ggplot(Pcscores_cr) +
  aes(Comp1, Comp2, color = life_cycle, shape = country) +
  # geom_text_repel(aes(label = rownames(Pcscores_cr)),
  #                 size = 3,
  #                 max.overlaps = Inf) +
  geom_polygon(
    data = hulls_lin,
    aes(fill = life_cycle, group = interaction(life_cycle)),
    alpha = 0.2,
    color = "gray80",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  ) +
  scale_fill_manual(
    name = "Life cycle",
    values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  ) +
  scale_shape_manual(values = c("Italy" = 17, "Bosnia" = 1, "Switzerland" = 8)) +
  geom_point(size = 2) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "gray90"),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)) +
  xlab("PC1: 71.7 %") +
  ylab("PC2: 3.3 %") +
  coord_fixed(ratio = 1)
plot_PCA_cr


ggsave("/plot_PCA_cr.png",
       plot = plot_PCA_cr, width = 8, height = 6, dpi = 300) 


# Calculate the distance from the origin for each point
Pcscores_cr$distance_to_origin <- sqrt(Pcscores_cr$Comp1^2 + Pcscores_cr$Comp2^2)
# Find the line with the minimum distance
consensus_specimen <- Pcscores_cr[which.min(Pcscores_cr$distance_to_origin), ]
# Display the reference specimen label
consensus_specimen$Row.names #[1] "A_Ichthyosaura_alpestris_MNHN_1992-2685_cranium"


Consensus_mesh_cr <- vcgPlyRead("Data/Landmarks/Consensus_specimens/A_Ichthyosaura_alpestris_MNHN_1992-2685_cranium.ply",
                                 clean = TRUE)
load(file = "Data/Landmarks/pts/landmarks_cranium.Rdata")

dimnames(landmarks_cranium)
Consensus_cr <- landmarks_cranium[,, 14]   

## Display the minimum and maximum deformations along the axes

#PC min and max
PC1max<-PCA_cr$shapes$shapes.comp1$max
PC1min<-PCA_cr$shapes$shapes.comp1$min
# PC2max<-PCA_cr$shapes$shapes.comp2$max
# PC2min<-PCA_cr$shapes$shapes.comp2$min

# max PC1
warpcr1 <- tps3d(Consensus_mesh_cr, Consensus_cr, PC1max)
shade3d(warpcr1, col="bisque", specular=1)
spheres3d(PC1max,col="red", radius = 0.003)
# min PC1
warpcr2 <- tps3d(Consensus_mesh_cr, Consensus_cr, PC1min)
shade3d(warpcr2, col="bisque", specular=1)
spheres3d(PC1min,col="skyblue1", radius = 0.003)
# # max PC2
# warpcr3 <- tps3d(Consensus_mesh_cr, Consensus_cr, PC2max)
# shade3d(warpcr3, col="bisque", specular=1)
# spheres3d(PC2max,col="red",  radius = 0.003)
# # min PC2
# warpcr4 <- tps3d(Consensus_mesh_cr, Consensus_cr, PC2min)
# shade3d(warpcr4, col="bisque", specular=1)
# spheres3d(PC2min,col="skyblue1",  radius = 0.003)


# Mandible
#_______________________________________________________________________________
load(file = "Data/Landmarks/pts/Y.gpa_md.Rdata")
ident_md <- read.csv("Data/Landmarks/Ident/ident_&_Csize_md.csv",
                     sep=";",header = T, row.names = 1)
Y_md<- Y.gpa_md$rotated

# PCA
PCA_md <- gm.prcomp(Y_md)

# % Variance Explained
summary(PCA_md)
# Comp1 = 66 %
# Comp2 = 10 %
# Comp3 = 5.4 %
# Comp4 = 3.5 %

Pcscores_md <- data.frame(PCA_md$x[,1:4])

Pcscores_md <- merge(Pcscores_md,
                     ident_md,
                     by = 'row.names',
                     all = T)

write.csv(Pcscores_md,
          "/Pcscores_md.csv",
          row.names = TRUE)

# PCA plot
hulls_md <- Pcscores_md %>%
  group_by(life_cycle) %>%
  slice(chull(Comp1, Comp2))
plot_PCA_md <- ggplot(Pcscores_md) +
  aes(Comp1, Comp2, color = life_cycle, shape = country) +
  # geom_text_repel(aes(label = rownames(Pcscores_md)),
  #                 size = 3,
  #                 max.overlaps = Inf) +
  geom_polygon(
    data = hulls_md,
    aes(fill = life_cycle, group = interaction(life_cycle)),
    alpha = 0.2,
    color = "gray80",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF", "metamorphic" = "#009900")
  ) +
  scale_fill_manual(
    name = "Life cycle",
    values = c("paedomorphic" = "#6633FF", "undergoing metamorphosis" = "#00CCFF",  "metamorphic" = "#009900")
  ) +
  scale_shape_manual(values = c("Italy" = 17, "Bosnia" = 1, "Switzerland" = 8)) +
  geom_point(size = 2) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_line(color = "gray90"),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)) +
  xlab("PC1: 66 %") +
  ylab("PC2: 10 %") +
  coord_fixed(ratio = 1)
plot_PCA_md


ggsave("/plot_PCA_md.png",
       plot = plot_PCA_md, width = 8, height = 6, dpi = 300) 


# Calculate the distance from the origin for each point
Pcscores_md$distance_to_origin <- sqrt(Pcscores_md$Comp1^2 + Pcscores_md$Comp2^2)
# Find the line with the minimum distance
consensus_specimen <- Pcscores_md[which.min(Pcscores_md$distance_to_origin), ]
# Display the reference specimen label
consensus_specimen$Row.names #[1] "A_Ichthyosaura_alpestris_MNHN_1992-2823_mandible"


Consensus_mesh_md <- vcgPlyRead("Data/Landmarks/Consensus_specimens/A_Ichthyosaura_alpestris_MNHN_1992-2823_mandible.ply",
                                clean = TRUE)
coronoid <- vcgPlyRead("Data/Landmarks/Consensus_specimens/A_Ichthyosaura_alpestris_MNHN_1992-2823_coronoid.ply",
                                clean = TRUE)
load(file = "Data/Landmarks/pts/landmarks_mandible_slided.Rdata")

dimnames(landmarks_mandible_slided$dataslide)
Consensus_md <- landmarks_mandible_slided$dataslide[,, 27]   

## Display the minimum and maximum deformations along the axes

#PC min and max
PC1max<-PCA_md$shapes$shapes.comp1$max
PC1min<-PCA_md$shapes$shapes.comp1$min
# PC2max<-PCA_md$shapes$shapes.comp2$max
# PC2min<-PCA_md$shapes$shapes.comp2$min

# max PC1
warpmd1 <- tps3d(Consensus_mesh_md, Consensus_md, PC1max)
warpcoro <- tps3d(coronoid, Consensus_md, PC1max)
shade3d(warpmd1, col="bisque", specular=1)
shade3d(warpcoro, col="#CCCCCC", specular=1, alpha = 0.5)
spheres3d(PC1max,col="red", radius = 0.003)
# min PC1
warpmd2 <- tps3d(Consensus_mesh_md, Consensus_md, PC1min)
shade3d(warpmd2, col="bisque", specular=1)
spheres3d(PC1min,col="skyblue1", radius = 0.003)
# # max PC2
# warpmd3 <- tps3d(Consensus_mesh_md, Consensus_md, PC2max)
# warpcoro <- tps3d(coronoid, Consensus_md, PC2max)
# shade3d(warpmd3, col="bisque", specular=1)
# shade3d(warpcoro, col="#CCCCCC", specular=1, alpha = 0.5)
# spheres3d(PC2max,col="red",  radius = 0.003)
# # min PC2
# warpmd4 <- tps3d(Consensus_mesh_md, Consensus_md, PC2min)
# warpcoro <- tps3d(coronoid, Consensus_md, PC2min)
# shade3d(warpmd4, col="bisque", specular=1)
# shade3d(warpcoro, col="#CCCCCC", specular=1, alpha = 0.5)
# spheres3d(PC2min,col="skyblue1",  radius = 0.003)

###############################################################################
## ------------------------------------------------------------------------- ##
##                   Percentage of hyoid bones presence                      ##
#### ----------------------------------------------------------------------####
# Packages
#_______________________________________________________________________________
library(dplyr)
library(tidyr)

# Data import
#_______________________________________________________________________________
df_hyoid <-read.csv("Data/list_hyoid.csv", 
                     sep=";")

unique(df_hyoid$life.cycle)
unique(df_hyoid$Bb)
unique(df_hyoid$mandible.coronoid.presence)

df_by_cycle <- split(df_hyoid, df_hyoid$life.cycle)

names(df_by_cycle)

bones <- c("Uh", "Bb", "Ch", "Hb.1", "Hb.2", "Cb.1", "Cb.2", "mandible.coronoid.presence")

# Percentage
#_______________________________________________________________________________
calc_presence <- function(df) {
  sapply(bones, function(b) {
    mean(df[[b]] == "yes", na.rm = TRUE) * 100
  })
}

presence_by_cycle <- lapply(df_by_cycle, calc_presence)

presence_by_cycle

presence_df <- do.call(rbind, presence_by_cycle)
presence_df <- as.data.frame(presence_df)

presence_df
# Uh        Bb       Ch       Hb.1     Hb.2     Cb.1     Cb.2 mandible.coronoid.presence
# paedomorphic             100.000000  1.904762 19.04762   1.904762 0.000000 17.14286 3.809524                  97.142857
# triphasic                  0.000000 96.969697 98.48485  96.969697 3.030303 98.48485 0.000000                   1.515152
# undergoing metamorphosis   3.703704 77.777778 96.29630 100.000000 0.000000 88.88889 3.703704                   3.703704

write.csv(presence_df,
          "/Analyses/Results/presence_hyoid.csv",
          row.names = FALSE)

###############################################################################


