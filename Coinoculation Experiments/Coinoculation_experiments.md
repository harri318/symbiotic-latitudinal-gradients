Co-inoculation experiments
================
Tia Harrison
09/04/2021

# Set up

## Load the packages for analysis

``` r
# Packages for analysis
library(car) 
library(lsmeans)
library(multcomp)
library(lme4)
library(performance)
library(see)
library(gridExtra)
library(wesanderson)
library(RColorBrewer)
library(R2admb)
library(glmmADMB)
library(dplyr)
library(tidyverse)
library(nationalparkcolors)
library(ggpubr)
```

## Set contrasts

When performing an Anova, should set contrasts to use effects coding
instead of the default (treatment contrasts).

``` r
options(contrasts = c("contr.sum","contr.poly")) 
```

## Overdispersion function

Ben Bolker’s function for testing overdispersion in poisson models
<https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion>

``` r
overdisp_fun <- function(model) {
rdf <- df.residual(model)
rp <- residuals(model,type="pearson")
Pearson.chisq <- sum(rp^2)
prat <- Pearson.chisq/rdf
pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
```

# Data

## Load and explore the two datasets for analysis

The Bacillus data set was collected from the growth chamber experiment
where plants that are tropical in origin (Costa Rica) were used in the
experiment. The Pseudomonas data set was collected from the greenhouse
experiment and the seeds and strains used in that experiment were from a
temperate population in Virginia.

``` r
bacillus <- read_csv("Bacillus_data.csv") # This data already has the initial dead plants removed from it and has 95 plants in it
pseudomonas <- read_csv("Pseudomonas_data.csv") # This data set is still the full dataset and has 270 plants in it 

# Convert biomass to mg in bacillus dataset 
# Pseudomonas data set is already in mg 
# Convert to factor 
bacillus <- bacillus%>% 
  mutate(Aboveground_mg= Aboveground_g*1000) %>%
  mutate(Belowground_mg=Belowground_g*1000) %>%
  mutate(Rhizobia=as.factor(Rhizobia), 
         Bacillus=as.factor(Bacillus), 
         Block=as.factor(Block), 
         Plant_ID=as.factor(Plant_ID))

# Remove plants that died before treatment took effect 
# Convert relevant data to factors
pseudomonas2 <- pseudomonas %>% 
  filter(Mortality1 == 1) %>%
  mutate(Bradyrhizobia=as.factor(Bradyrhizobia), 
         Pseudomonas=as.factor(Pseudomonas), 
         Variovorax=as.factor(Variovorax),
         Block=as.factor(Block), 
         Plant_ID=as.factor(Plant_ID)) # Now the experiment has 216 plants in it 

# Break down of how many plants alive at the end of the experiment 
bacillus3 <- bacillus %>%
  filter(Dead ==0) # 64 plants but if they are in the treatments with rhizobia might be enough to compare 

pseudomonas3 <- pseudomonas2 %>%
  filter(Dead ==0) # 164 plants alive at the end of the experiment 
```

# Analysis

## Tropical non-rhizobia strains

Test if the Bacillus strain have an impact on plant growth, nodule
number, and nodule weight belowground, and whether or not there is an
interaction with rhizobia inoculations.

``` r
# histogram of aboveground data 
hist(bacillus3$Aboveground_mg) 
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-1.png)<!-- -->

``` r
# aboveground model 
mod_full_s <- lmer(log10(Aboveground_mg) ~ Bacillus*Rhizobia + (1|Family), data=bacillus3)

# Check assumptions 
plot(mod_full_s)
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-2.png)<!-- -->

``` r
qqnorm(resid(mod_full_s))
qqline(resid(mod_full_s))
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-3.png)<!-- -->

``` r
plot(fitted(mod_full_s), sqrt(abs(resid(mod_full_s))), main="Scale-location")
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-4.png)<!-- -->

``` r
# Tried a bunch of things but the results came out all the same and this is the best fitting model 

# Test significance and look at summary 
summary(mod_full_s)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log10(Aboveground_mg) ~ Bacillus * Rhizobia + (1 | Family)
    ##    Data: bacillus3
    ## 
    ## REML criterion at convergence: 111.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5309 -0.3107 -0.0253  0.7237  1.6745 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Family   (Intercept) 0.0000   0.000   
    ##  Residual             0.2851   0.534   
    ## Number of obs: 64, groups:  Family, 5
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)          1.60441    0.06761  23.729
    ## Bacillus1            0.07622    0.06761   1.127
    ## Rhizobia1           -0.92241    0.06761 -13.642
    ## Bacillus1:Rhizobia1 -0.06830    0.06761  -1.010
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Bclls1 Rhizb1
    ## Bacillus1   0.108               
    ## Rhizobia1   0.108  0.083        
    ## Bclls1:Rhz1 0.083  0.108  0.108 
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
Anova(mod_full_s, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log10(Aboveground_mg)
    ##                      Chisq Df Pr(>Chisq)    
    ## (Intercept)       563.0776  1     <2e-16 ***
    ## Bacillus            1.2707  1     0.2596    
    ## Rhizobia          186.1172  1     <2e-16 ***
    ## Bacillus:Rhizobia   1.0205  1     0.3124    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Nodule analysis 
hist(bacillus3$Nodule_no) 
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-5.png)<!-- -->

``` r
# Include random term for overdispersion 
mod_nod_s <- glmer(Nodule_no ~ Bacillus*Rhizobia + (1|Family) + (1|Plant_ID), family=poisson(link="sqrt"),  data=bacillus3)

# Test fit 
plot(mod_nod_s)
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-6.png)<!-- -->

``` r
qqnorm(resid(mod_nod_s))
qqline(resid(mod_nod_s))
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-7.png)<!-- -->

``` r
plot(fitted(mod_nod_s), sqrt(abs(resid(mod_nod_s))), main="Scale-location")
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-8.png)<!-- -->

``` r
# Check overdispersion 
overdisp_fun(mod_nod_s) # if non-significant then no overdispersion 
```

    ##       chisq       ratio         rdf           p 
    ##  4.34953082  0.07499191 58.00000000  1.00000000

``` r
# Test significance and look at summary 
summary(mod_nod_s)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( sqrt )
    ## Formula: Nodule_no ~ Bacillus * Rhizobia + (1 | Family) + (1 | Plant_ID)
    ##    Data: bacillus3
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##      470      483     -229      458       58 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -0.99780 -0.10156 -0.05417  0.19828  0.59896 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant_ID (Intercept) 3.429    1.852   
    ##  Family   (Intercept) 0.000    0.000   
    ## Number of obs: 64, groups:  Plant_ID, 64; Family, 5
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)           3.9489     0.2460  16.052   <2e-16 ***
    ## Bacillus1             0.1837     0.2460   0.747    0.455    
    ## Rhizobia1            -3.3371     0.2460 -13.567   <2e-16 ***
    ## Bacillus1:Rhizobia1  -0.3699     0.2460  -1.504    0.133    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Bclls1 Rhizb1
    ## Bacillus1   0.113               
    ## Rhizobia1   0.129  0.088        
    ## Bclls1:Rhz1 0.087  0.128  0.113 
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
Anova(mod_nod_s, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: Nodule_no
    ##                      Chisq Df Pr(>Chisq)    
    ## (Intercept)       257.6721  1     <2e-16 ***
    ## Bacillus            0.5579  1     0.4551    
    ## Rhizobia          184.0669  1     <2e-16 ***
    ## Bacillus:Rhizobia   2.2612  1     0.1327    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Nodule weight 
# Only performed on plants that formed nodules and so we removed Control and treatments without rhizobia
bacillus_r <- bacillus3 %>%
  filter(Rhizobia=="YES")

# Look at the data 
hist(bacillus_r$Nodule_weight_ug) 
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-9.png)<!-- -->

``` r
# Run model 
mod_weight <- lmer(Nodule_weight_ug ~ Bacillus + (1|Family), data=bacillus_r)

# Test significance and look at summary 
Anova(mod_weight, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: Nodule_weight_ug
    ##               Chisq Df Pr(>Chisq)    
    ## (Intercept) 64.8345  1  8.146e-16 ***
    ## Bacillus     0.0477  1     0.8271    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(mod_weight)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: Nodule_weight_ug ~ Bacillus + (1 | Family)
    ##    Data: bacillus_r
    ## 
    ## REML criterion at convergence: 575.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.93527 -0.70300  0.08072  0.75096  2.01390 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Family   (Intercept) 3468343  1862    
    ##  Residual             8371619  2893    
    ## Number of obs: 32, groups:  Family, 5
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   7915.1      983.0   8.052
    ## Bacillus1      112.8      516.5   0.218
    ## 
    ## Correlation of Fixed Effects:
    ##           (Intr)
    ## Bacillus1 0.028

``` r
# Check assumptions 
plot(mod_weight)
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-10.png)<!-- -->

``` r
qqnorm(resid(mod_weight))
qqline(resid(mod_weight))
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-11.png)<!-- -->

``` r
plot(fitted(mod_weight), sqrt(abs(resid(mod_weight))), main="Scale-location") # everything looks good 
```

![](Coinoculation_experiments_files/figure-gfm/tropical%20nr_model-12.png)<!-- -->

## Temperate non-rhizobia strains

Test whether the non-rhizobia strains Pseudomonas and Variovorax have
impacts on plant growth and nodule formation when coinoculation with
rhizobia. Look for interactions between all three strains inoculated on
the plants.

``` r
# histogram of aboveground data 
hist(pseudomonas3$Aboveground_mg) 
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-1.png)<!-- -->

``` r
# Model on full data set with all treatments
mod_full_ps <- lmer(log10(Aboveground_mg) ~ Bradyrhizobia*Pseudomonas*Variovorax + (1|Family) + (1|Block), data=pseudomonas3)

# Check assumptions 
plot(mod_full_ps)
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-2.png)<!-- -->

``` r
qqnorm(resid(mod_full_ps))
qqline(resid(mod_full_ps))
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-3.png)<!-- -->

``` r
plot(fitted(mod_full_ps), sqrt(abs(resid(mod_full_ps))), main="Scale-location") # Pretty good fit! 
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-4.png)<!-- -->

``` r
# Test significance and look at summary 
summary(mod_full_ps)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log10(Aboveground_mg) ~ Bradyrhizobia * Pseudomonas * Variovorax +  
    ##     (1 | Family) + (1 | Block)
    ##    Data: pseudomonas3
    ## 
    ## REML criterion at convergence: 230.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6066 -0.6221 -0.1351  0.4066  2.8944 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Family   (Intercept) 0.01083  0.1041  
    ##  Block    (Intercept) 0.02785  0.1669  
    ##  Residual             0.17818  0.4221  
    ## Number of obs: 163, groups:  Family, 11; Block, 10
    ## 
    ## Fixed effects:
    ##                                          Estimate Std. Error t value
    ## (Intercept)                              1.044769   0.070162  14.891
    ## Bradyrhizobia1                          -0.105277   0.034273  -3.072
    ## Pseudomonas1                            -0.011420   0.034140  -0.334
    ## Variovorax1                              0.020999   0.034019   0.617
    ## Bradyrhizobia1:Pseudomonas1              0.012986   0.033645   0.386
    ## Bradyrhizobia1:Variovorax1              -0.007474   0.034114  -0.219
    ## Pseudomonas1:Variovorax1                -0.002204   0.034178  -0.064
    ## Bradyrhizobia1:Pseudomonas1:Variovorax1  0.011069   0.033906   0.326
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Brdyr1 Psdmn1 Vrvrx1 Br1:P1 Br1:V1 Ps1:V1
    ## Bradyrhizb1  0.038                                          
    ## Pseudomons1 -0.013 -0.086                                   
    ## Variovorax1  0.006  0.067  0.006                            
    ## Brdyrhz1:P1 -0.043 -0.021  0.058  0.071                     
    ## Brdyrhz1:V1  0.033  0.031  0.062  0.077  0.019              
    ## Psdmns1:Vr1  0.012  0.083  0.008 -0.020  0.062 -0.090       
    ## Brdy1:P1:V1  0.030  0.024  0.077 -0.094  0.016 -0.022  0.071

``` r
Anova(mod_full_ps, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log10(Aboveground_mg)
    ##                                         Chisq Df Pr(>Chisq)    
    ## (Intercept)                          221.7347  1  < 2.2e-16 ***
    ## Bradyrhizobia                          9.4354  1   0.002128 ** 
    ## Pseudomonas                            0.1119  1   0.738008    
    ## Variovorax                             0.3810  1   0.537067    
    ## Bradyrhizobia:Pseudomonas              0.1490  1   0.699522    
    ## Bradyrhizobia:Variovorax               0.0480  1   0.826592    
    ## Pseudomonas:Variovorax                 0.0042  1   0.948576    
    ## Bradyrhizobia:Pseudomonas:Variovorax   0.1066  1   0.744084    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Nodules 
# Data is overdispersed and zero -inflated but when we tried zero inflated models the model was overparametized
# Ttesting the rhizobia treatments only worked best 
pseudomonas_r <- pseudomonas3 %>%
  filter(Bradyrhizobia=="Y")

# look at the data 
hist(pseudomonas3$Nodule_no) 
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-5.png)<!-- -->

``` r
# Fit the model 
mod_nod_glm_full <- glmer(Nodule_no ~ Pseudomonas*Variovorax + (1|Family) + (1|Block) + (1|Plant_ID), family=poisson(link="sqrt"), data=pseudomonas_r)

# Check assumptions 
plot(mod_nod_glm_full)
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-6.png)<!-- -->

``` r
qqnorm(resid(mod_nod_glm_full))
qqline(resid(mod_nod_glm_full))
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-7.png)<!-- -->

``` r
plot(fitted(mod_nod_glm_full), sqrt(abs(resid(mod_nod_glm_full))), main="Scale-location")
```

![](Coinoculation_experiments_files/figure-gfm/temperate%20nr_model-8.png)<!-- -->

``` r
overdisp_fun(mod_nod_glm_full)
```

    ##       chisq       ratio         rdf           p 
    ##  6.36833412  0.07960418 80.00000000  1.00000000

``` r
# Check significance  
summary(mod_nod_glm_full)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( sqrt )
    ## Formula: Nodule_no ~ Pseudomonas * Variovorax + (1 | Family) + (1 | Block) +  
    ##     (1 | Plant_ID)
    ##    Data: pseudomonas_r
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    492.1    509.4   -239.1    478.1       80 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -0.3828 -0.2193 -0.1266  0.2395  0.6125 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant_ID (Intercept) 2.8197   1.6792  
    ##  Family   (Intercept) 0.1431   0.3783  
    ##  Block    (Intercept) 0.6198   0.7873  
    ## Number of obs: 87, groups:  Plant_ID, 87; Family, 11; Block, 10
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               1.57453    0.33781   4.661 3.15e-06 ***
    ## Pseudomonas1             -0.13483    0.20096  -0.671    0.502    
    ## Variovorax1              -0.07619    0.19632  -0.388    0.698    
    ## Pseudomonas1:Variovorax1 -0.12215    0.19672  -0.621    0.535    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Psdmn1 Vrvrx1
    ## Pseudomons1  0.046              
    ## Variovorax1 -0.027 -0.039       
    ## Psdmns1:Vr1 -0.026 -0.083  0.084

``` r
Anova(mod_nod_glm_full, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: Nodule_no
    ##                          Chisq Df Pr(>Chisq)    
    ## (Intercept)            21.7247  1  3.147e-06 ***
    ## Pseudomonas             0.4502  1     0.5023    
    ## Variovorax              0.1506  1     0.6980    
    ## Pseudomonas:Variovorax  0.3855  1     0.5347    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Tests of variance

Testing whether variance in nodule number in the rhiozbia only
treatments is different from in coinoculation treatments in both the
tropical and temperate experiment data.

``` r
# tropical experiment 
trop_variance <- bartlett.test(Nodule_no ~ Bacillus, bacillus_r)
trop_variance
```

    ## 
    ##  Bartlett test of homogeneity of variances
    ## 
    ## data:  Nodule_no by Bacillus
    ## Bartlett's K-squared = 0.11394, df = 1, p-value = 0.7357

``` r
# temperate experiment 
temp_variance <- leveneTest(Nodule_no ~ Pseudomonas*Variovorax, pseudomonas_r)
temp_variance
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  3  0.5045 0.6802
    ##       83

# Figures

Plotting the data for the tropical and temperate datasets

## Plotting aboveground biomass

``` r
# Specify the colour palette 
treatments3 <- park_palette("RockyMountains")
grey_col <- brewer.pal(9, "Greys") 
 
# Prep the data 
bacillus34 <- bacillus3 %>%
  mutate(rhizobia_new=ifelse(Rhizobia=="YES", "Rhizobia", "Control")) %>%
  mutate(bacillus_new=ifelse(Bacillus=="YES", "          Bacillus", "Control")) %>%
  mutate(bacillus_new=fct_relevel(bacillus_new, "Control", "          Bacillus"))

# Plot tropical experiment full dataset 
bacillus_mass <- ggplot(bacillus34, aes(x=bacillus_new, y=log10(Aboveground_g), fill=Rhizobia)) + 
  geom_boxplot(fatten=1) +
  (ylab("Log biomass (mg)"))+
  (xlab(""))+
  labs(tag = "(a)") +
  facet_grid(~rhizobia_new, switch="x") +
  scale_fill_manual(values = c("white", grey_col[4])) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_text(size=9, hjust = 0.98, angle=75),
        strip.background =element_rect(fill=grey_col[2]))

# Prep temperate data 
pseudomonas34 <- pseudomonas3 %>%
  mutate(rhizobia_new=ifelse(Bradyrhizobia=="Y", "Rhizobia", "Control")) %>%
  mutate(treatment_new=ifelse(Pseudomonas =="N" & Variovorax =="N", "Control", Treatment)) %>%
  mutate(treatment_new=ifelse(Pseudomonas=="Y" & Variovorax=="N", "Pseudomonas", treatment_new)) %>%
  mutate(treatment_new=ifelse(Pseudomonas=="N" & Variovorax=="Y", "Variovorax", treatment_new)) %>%
  mutate(treatment_new=ifelse(Pseudomonas=="Y" & Variovorax=="Y", "All strains", treatment_new)) %>%
  mutate(treatment_new=fct_relevel(treatment_new, "Control", "Pseudomonas", "Variovorax", "All strains"))

# Plot the temperate data 
pseud_mass <- ggplot(pseudomonas34, aes(x=treatment_new, y=log10(Aboveground_mg), fill=rhizobia_new)) + 
  geom_boxplot(fatten=1) +
  (ylab("Log biomass (mg)"))+
  (xlab(""))+
  labs(tag = "(b)") +
  facet_wrap(~rhizobia_new, switch="x") +
  scale_fill_manual(values = c("white", grey_col[4])) +
  theme_classic() +
  theme(legend.position="none", 
        axis.text.x=element_text(size=9, hjust = 0.98, angle=75),
        strip.background =element_rect(fill=grey_col[2]))

# Look at the plot 
grid.arrange(bacillus_mass, pseud_mass,  
             ncol=3, nrow=1, 
             layout_matrix=rbind(c(1,2,2)))
```

![](Coinoculation_experiments_files/figure-gfm/biomass%20plot-1.png)<!-- -->

``` r
# Save the plot 
# pdf(file="Biomass_plot.pdf", width=7, height=5)

test3<-grid.arrange(bacillus_mass, pseud_mass,  
             ncol=3, nrow=1, 
             layout_matrix=rbind(c(1,2,2)))
```

![](Coinoculation_experiments_files/figure-gfm/biomass%20plot-2.png)<!-- -->

``` r
# dev.off()
ggsave(filename="Figures/Biomass_plot.pdf", plot=test3)
```

## Plotting nodule number

``` r
# Plot tropical experiment full dataset 
bacillus_nod <- ggplot(bacillus34, aes(x=bacillus_new, y=Nodule_no, fill=Rhizobia)) + 
  geom_boxplot(fatten=1) +
  (ylab("Nodules (no)"))+
  (xlab(""))+
  labs(tag = "(a)") +
  facet_grid(~rhizobia_new, switch="x") +
  scale_fill_manual(values = c("white", grey_col[4])) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x=element_text(size=9, hjust = 0.98, angle=75),
        strip.background =element_rect(fill=grey_col[2]))

# Plot the temperate data 
pseud_nod <- ggplot(pseudomonas34, aes(x=treatment_new, y=Nodule_no, fill=rhizobia_new)) + 
  geom_boxplot(fatten=1) +
  (ylab("Nodules (no)"))+
  (xlab(""))+
  labs(tag = "(b)") +
  facet_wrap(~rhizobia_new, switch="x") +
  scale_fill_manual(values = c("white", grey_col[4])) +
  theme_classic() +
  theme(legend.position="none", 
        axis.text.x=element_text(size=9, hjust = 0.98, angle=75),
        strip.background =element_rect(fill=grey_col[2]))

# Look at the plot 
grid.arrange(bacillus_nod, pseud_nod,  
             ncol=3, nrow=1, 
             layout_matrix=rbind(c(1,2,2)))
```

![](Coinoculation_experiments_files/figure-gfm/nodule%20plot-1.png)<!-- -->

``` r
# Save the plot 
# pdf(file="Nodule_plot.pdf", width=7, height=5)

test4<-grid.arrange(bacillus_nod, pseud_nod,  
             ncol=3, nrow=1, 
             layout_matrix=rbind(c(1,2,2)))
```

![](Coinoculation_experiments_files/figure-gfm/nodule%20plot-2.png)<!-- -->

``` r
#dev.off() 
ggsave(filename="Figures/Nodule_plot.pdf", plot=test4)
```
