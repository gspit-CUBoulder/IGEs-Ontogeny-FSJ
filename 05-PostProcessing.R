############# Table and Figure Generation for Indirect Genetic Effects Across Ontogeny #####################

#### Script 5/5 in Florida Scrub Jay IGE Analysis
#### Original creator: Gladiana Spitz

#### Goal: This script processes modeling, testing, and dilution parameter data 
#### into tables and figures for the main text and supplemental

#### Output: Main text tables (Table2_MainText.csv) and
#### figures (Figure1_MainText.png, Figure2_MainText.png)
#### Supplement tables (Suppl_Table1.csv, Suppl_Table2.csv, Suppl_Table3.csv, Suppl_Table4.csv)
#### and figures (Figure1_Supplement.png)


# Package Set Up
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(viridis)
library(tidyverse)
library(kableExtra)
library(asremlPlus)
library(emdbook)

# Load output data 

## Modeling data
load("ModelResults.rdata")

## Testing data 
load("TestingResults.Rdata")
load("TestingModels.Rdata")

#### Data Post Processing ####

# clean total effects table - contains social effect information
totaleffects_all <- totaleffects %>% 
  filter(!(Trait %in% c("d11.c","juv.c"))) %>%
  # add column of complete trait names
  mutate(Trait_full = ifelse(Trait == "d11.w", "Nestling Weight",Trait),
         Trait_full = ifelse(Trait_full == "juv.w", "Juvenile Weight",Trait_full),
         Trait_full = ifelse(Trait_full == "d11.t", "Nestling Tarsus",Trait_full),
         Trait_full = ifelse(Trait_full == "juv.t", "Juvenile Tarsus",Trait_full),
         Trait_full = factor(Trait_full, levels = c("Nestling Weight", "Juvenile Weight",
                                                    "Nestling Tarsus", "Juvenile Tarsus",
                                                    "Nestling Condition", "Juvenile Condition")),
         ## Rename models
         Model = ifelse(Model == "Maternal", "Mom", Model),
         Model = ifelse(Model == "Paternal", "Dad", Model),
         Model = ifelse(Model == "Helper", "Helper", Model),
         Model = ifelse(Model == "MatPat", "Mom + Dad", Model),
         Model = ifelse(Model == "MatHelper", "Mom + Helper", Model),
         Model = ifelse(Model == "PatHelper", "Dad + Helper", Model),
         Model = ifelse(Model == "MatPatHelper", "Mom + Dad + Helper", Model),
         # factor model column for better plotting
         Model = factor(Model, levels=c("DGE-Only","Mom","Dad","Helper","Mom Cov","Dad Cov","Helper Cov","Mom + Dad","Mom + Helper","Dad + Helper","Mom + Dad + Helper")),
         # make indirect effect name column1
         `Indirect Effect` = ifelse(type == "GenetICS","Genetic",type),
         `Indirect Effect` = ifelse(`Indirect Effect` == "environmental","Environmental",`Indirect Effect`),
         fullsig = Estimate - 2*SE,
         sigmark = ifelse(fullsig > 0, " *", ""),
         full = paste(round(Estimate,2)," (",
                      round(SE,2),")",sigmark,sep=""))

# Split total effects into different traits
totaleffects_weight <- subset(totaleffects_all, Trait_full %in% c("Nestling Weight","Juvenile Weight"))

totaleffects_tarsus <- subset(totaleffects_all, Trait_full %in% c("Nestling Tarsus","Juvenile Tarsus"))


# Clean heritability table - contains heritabilty parameter for each model and trait
herit_all <- herit %>%
  filter(!(Trait %in% c("d11.c","juv.c"))) %>%
  # add column of complete trait names
  mutate(Trait_full = ifelse(Trait == "d11.w", "Nestling Weight",Trait),
         Trait_full = ifelse(Trait_full == "juv.w", "Juvenile Weight",Trait_full),
         Trait_full = ifelse(Trait_full == "d11.t", "Nestling Tarsus",Trait_full),
         Trait_full = ifelse(Trait_full == "juv.t", "Juvenile Tarsus",Trait_full),
         type = "h2",
         Trait_full = factor(Trait_full, levels = c("Nestling Weight", "Juvenile Weight",
                                                    "Nestling Tarsus", "Juvenile Tarsus",
                                                    "Nestling Condition", "Juvenile Condition")))
# add column for the type of parameter

# Clean social heritability parameter table
params_all <- params %>% 
  filter(!(Trait %in% c("d11.c","juv.c"))) %>%
  mutate(Trait_full = ifelse(Trait == "d11.w", "Nestling Weight",Trait),
         Trait_full = ifelse(Trait_full == "juv.w", "Juvenile Weight",Trait_full),
         Trait_full = ifelse(Trait_full == "d11.t", "Nestling Tarsus",Trait_full),
         Trait_full = ifelse(Trait_full == "juv.t", "Juvenile Tarsus",Trait_full),
         Trait_full = factor(Trait_full, levels = c("Nestling Weight", "Juvenile Weight",
                                                    "Nestling Tarsus", "Juvenile Tarsus",
                                                    "Nestling Condition", "Juvenile Condition")))
# Make heritability parameters table and remove TBV (same as tau in no covariance model)
Params_all <- rbind(params_all,herit_all) %>% 
  mutate(Model = ifelse(Model == "Maternal", "Mom", Model),
         Model = ifelse(Model == "Paternal", "Dad", Model),
         Model = ifelse(Model == "Helper", "Helper", Model),
         Model = ifelse(Model == "MatPat", "Mom + Dad", Model),
         Model = ifelse(Model == "MatHelper", "Mom + Helper", Model),
         Model = ifelse(Model == "PatHelper", "Dad + Helper", Model),
         Model = ifelse(Model == "MatPatHelper", "Mom + Dad + Helper", Model),
         Model = factor(Model, levels=c("DGE-Only","Mom","Dad","Helper","Mom Cov","Dad Cov","Helper Cov","Mom + Dad","Mom + Helper","Dad + Helper","Mom + Dad + Helper")),
         fullsig = Estimate - 2*SE,
         sigmark = ifelse(fullsig > 0, " *", ""),
         full = paste(round(Estimate,2)," (",
                      round(SE,2),")",sigmark,sep=""))
# split parameter table into different traits
params_weight <- subset(Params_all, Trait_full %in% c("Nestling Weight","Juvenile Weight"))
params_tarsus <- subset(Params_all, Trait_full %in% c("Nestling Tarsus","Juvenile Tarsus"))

# Clean varcomps table - table containing all variance components
Varcomps_all <- Varcomps %>%
  mutate(varcomp = rownames(Varcomps),
         varcomp = ifelse(grepl("!diag\\(2\\)_1",varcomp), "vm(USFWS, ainv)",varcomp),
         varcomp = ifelse(grepl("!diag\\(2\\)_2",varcomp), "vm(Help1, ainv)",varcomp),
         varcomp = ifelse(grepl("!us\\(2\\)_1:1",varcomp), "vm(USFWS, ainv)",varcomp),
         varcomp = ifelse(grepl("vm\\(mom, ainv\\)!us\\(2\\)_2:1",varcomp), "Maternal Cov",varcomp),
         varcomp = ifelse(grepl("vm\\(dad, ainv\\)!us\\(2\\)_2:1",varcomp), "Paternal Cov",varcomp),
         varcomp = ifelse(grepl("vm\\(Help5, ainv\\)\\)!us\\(2\\)_2:1",varcomp), "Helper Cov",varcomp),
         varcomp = ifelse(grepl("vm\\(mom, ainv\\)!us\\(2\\)_2:2",varcomp), "Maternal Genetic",varcomp),
         varcomp = ifelse(grepl("vm\\(dad, ainv\\)!us\\(2\\)_2:2",varcomp), "Paternal Genetic",varcomp),
         varcomp = ifelse(grepl("vm\\(Help5, ainv\\)\\)!us\\(2\\)_2:2",varcomp), "Helper Genetic",varcomp),
         varcomp = gsub("[0-9]*$", "", varcomp),
         ## Rename variance components
         VC_full = ifelse(varcomp == "NatalYear", "Natal Year", varcomp),
         VC_full = ifelse(VC_full == "NatalNest", "Natal Nest", VC_full),
         VC_full = ifelse(VC_full %in% c("vm(USFWS, ainv)"), "Direct Genetic", VC_full),
         VC_full = ifelse(VC_full == "vm(mom, ainv)", "Maternal Genetic", VC_full),
         VC_full = ifelse(VC_full == "ide(mom)", "Maternal Environment", VC_full),
         VC_full = ifelse(VC_full == "units!units", "Residual", VC_full),
         VC_full = ifelse(VC_full == "vm(dad, ainv)", "Paternal Genetic", VC_full),
         VC_full = ifelse(VC_full == "ide(dad)", "Paternal Environment", VC_full),
         VC_full = ifelse(VC_full == "vm(Help1, ainv)", "Helper Genetic", VC_full),
         VC_full = ifelse(VC_full == "d:ide(Help1)", "Helper Environment", VC_full),
         VC_full = factor(VC_full, levels = c("Residual","Natal Year", "Natal Nest", "Direct Genetic",
                                              "Maternal Genetic","Maternal Environment", 
                                              "Paternal Genetic", "Paternal Environment", 
                                              "Helper Genetic", "Helper Environment")),
         ## Rename traits
         Trait_full = ifelse(Trait == "d11.w", "Nestling Weight",Trait),
         Trait_full = ifelse(Trait_full == "juv.w", "Juvenile Weight",Trait_full),
         Trait_full = ifelse(Trait_full == "d11.t", "Nestling Tarsus",Trait_full),
         Trait_full = ifelse(Trait_full == "juv.t", "Juvenile Tarsus",Trait_full),
         Trait_full = factor(Trait_full, levels = c("Nestling Weight", "Juvenile Weight",
                                                    "Nestling Tarsus", "Juvenile Tarsus")),
         ## Rename models
         Model = ifelse(Model == "Maternal", "Mom", Model),
         Model = ifelse(Model == "Paternal", "Dad", Model),
         Model = ifelse(Model == "Helper", "Helper", Model),
         Model = ifelse(Model == "MatPat", "Mom + Dad", Model),
         Model = ifelse(Model == "MatHelper", "Mom + Helper", Model),
         Model = ifelse(Model == "PatHelper", "Dad + Helper", Model),
         Model = ifelse(Model == "MatPatHelper", "Mom + Dad + Helper", Model),
         Model = factor(Model, levels=c("DGE-Only","Mom","Dad","Helper","Mom Cov","Dad Cov","Helper Cov","Mom + Dad","Mom + Helper","Dad + Helper","Mom + Dad + Helper")),
         ## Match Effect to sig_table Effect
         Effect = ifelse(Model == "Mom + Dad" & 
                           VC_full %in% c("Maternal Genetic","Maternal Environment",
                                          "Paternal Genetic","Paternal Environment"), 
                         as.character(Model), as.character(VC_full)),
         Effect = ifelse(Model == "Mom + Helper" & 
                           VC_full %in% c("Maternal Genetic","Maternal Environment",
                                          "Helper Genetic","Helper Environment"), 
                         as.character(Model), Effect),
         Effect = ifelse(Model == "Dad + Helper" & 
                           VC_full %in% c("Paternal Genetic","Paternal Environment",
                                          "Helper Genetic","Helper Environment"), 
                         as.character(Model), Effect),
         Effect = ifelse(Model == "Mom + Dad + Helper" & 
                           VC_full %in% c("Maternal Genetic","Maternal Environment",
                                          "Paternal Genetic","Paternal Environment",
                                          "Helper Genetic","Helper Environment"),
                         as.character(Model), Effect),
         # add pretty standard error
         std.error.chr = ifelse(is.na(std.error),"",paste(" (", round(std.error,2),")", sep = "")),
         fullsig = component - 2*std.error,
         fullsig = ifelse(is.na(fullsig),"",fullsig),
         sigmark = ifelse(fullsig > 0, " *", ""),
         full = paste(round(component,2),std.error.chr,sigmark, sep="")) %>% 
  filter(component != 1.00,
         !is.na(Effect))

# Split into different traits
var_weight <- subset(Varcomps_all, Trait_full %in% c("Nestling Weight","Juvenile Weight"))

var_tarsus <- subset(Varcomps_all, Trait_full %in% c("Nestling Tarsus",
                                                     "Juvenile Tarsus"))

# Total variance as calculated by summing all variance components
Vps <- Varcomps_all %>% group_by(Trait_full,Model) %>% dplyr::summarise(VP=round(sum(component),2))


#### Clean up modelFit to label model IDs with tested variance component

# named vectors of variance components
vcs <- c("1","2","3","4","5","6","7", "8a","8b","9")
VarianceComp_mat <- c("No RE","Natal Year","Natal Nest","Direct Genetic",
                      "Maternal Environment", "Maternal Genetic", "Cov",
                      "Mom + Dad", "Mom + Helper", "Mom + Dad + Helper")
names(VarianceComp_mat) <- vcs

VarianceComp_pat <- c("No RE","Natal Year","Natal Nest","Direct Genetic",
                      "Paternal Environment", "Paternal Genetic", "Cov",
                      "Mom + Dad", "Dad + Helper", "Mom + Dad + Helper")
names(VarianceComp_pat) <- vcs

VarianceComp_help <-  c("No RE","Natal Year","Natal Nest","Direct Genetic",
                        "Helper Environment", "Helper Genetic", "Cov",
                        "Mom + Helper", "Dad + Helper", "Mom + Dad + Helper")
names(VarianceComp_help) <- vcs

# traits
traits <- c("Nestling Weight", "Juvenile Weight", "Nestling Tarsus","Juvenile Tarsus")
names(traits) <-  c("d11.w","juv.w","d11.t","juv.t")

# social partners
social <- c("m","m","m","p","h")
names(social) <- c("none","all","maternal", "paternal", "helper")


df <- modelFit %>%
  mutate(Trait = traits[trait], # Cleaned up trait name
         # Assign variance component names based on social partner and model number in LRT
         # Maternal effects
         VComp = ifelse(socialPartner %in% c("maternal","none","all"), # none and all are same across models
                        VarianceComp_mat[model],
                        NA),
         # Paternal effects
         VComp = ifelse(socialPartner %in% c("paternal"),
                        VarianceComp_pat[model],
                        VComp),
         # Helper Effects
         VComp = ifelse(socialPartner %in% c("helper"),
                        VarianceComp_help[model],
                        VComp)) %>%
  rename(Model = model) %>%
  dplyr::select(-trait)

##### LRTs Extraction
sig_table <- data.frame()
# social partners
social <- c("maternal", "paternal", "helper","none","all")
names(social) <- c("m","p","h","n","a")


for(i in 1:length(LRTs)){
  # pull single LRT from list
  x <- LRTs[[i]]
  # get model name
  y <- names(LRTs)[i]
  x <- x %>%
    # Pull model comparison IDs from rownames
    mutate(mod_compared = as.character(row.names(x))) %>%
    # Make model tested vs compared to
    separate_wider_delim(mod_compared, delim = "/", names = c("Model","ComparedTo")) %>%
    # add row for m1 models
    add_row(Model = "m1") %>%
    # split model name into social partner and trait, assign each
    mutate(Model = str_remove_all(Model, "m"),
           ComparedTo = str_remove_all(ComparedTo, "m"),
           soc = str_split(y, "_")[[1]][1],
           socialPartner = social[soc],
           trt = str_split(y, "_")[[1]][2],
           Trait = traits[trt],
           # assign tested variance components based on social partner
           VComp = ifelse(soc == "m", VarianceComp_mat[Model], NA),
           VComp = ifelse(soc == "p", VarianceComp_pat[Model], VComp),
           VComp = ifelse(soc == "h", VarianceComp_help[Model],VComp),
           # reassign none/all
           socialPartner = ifelse(Model %in% c("1","2","3","4"), "none",socialPartner),
           socialPartner = ifelse(Model =="9", "all",socialPartner),
           # assign compared variance components based on social partner
           VC_Compare = ifelse(soc == "m", VarianceComp_mat[ComparedTo], NA),
           VC_Compare = ifelse(soc == "p", VarianceComp_pat[ComparedTo], VC_Compare),
           VC_Compare = ifelse(soc == "h", VarianceComp_help[ComparedTo],VC_Compare)) %>%
    dplyr::select(-soc, -trt)
  
  sig_table  <- bind_rows(sig_table,x)
}

sig_table <- sig_table %>% distinct()

#### Clean, round, and adjust table with signficance values

sig_table <- sig_table %>% 
  left_join(df) %>% 
  mutate(roundicant = ifelse(`Pr(Chisq)` <= 0.05, "*",""),
         across(names(.), as.character),
         socialPartner = factor(socialPartner, levels = social),
         Trait = factor(Trait, levels = traits),
         LogLik = round(as.numeric(LogLik), 2),
         AIC = round(as.numeric(AIC), 2),
         BIC = round(as.numeric(BIC), 2),
         `LR-statistic` = round(as.numeric(`LR-statistic`), 4),
         `Pr(Chisq)` = round(as.numeric(`Pr(Chisq)`), 2)) %>%
  arrange(Trait, Model, socialPartner) %>%
  dplyr::select(Trait, Model, socialPartner, Effect = VComp, LogLik, AIC, BIC, Reference_Model = VC_Compare,
                `LR-statistic`, `Pr(Chisq)`,roundicant)

All_RE_sig <- sig_table %>%
  group_by(Trait) %>%
  mutate(minAIC = min(AIC),
         deltaAIC = AIC - minAIC) %>%
  ungroup() %>% dplyr::select(-minAIC)

sig_table[is.na(sig_table)] <- ""

#### Main Text Tables and Figures ####

#### Tables

## Table 2 - Variance components for all models

# Function to make variance components names into more readable names that match manuscript
translator <- function(x){
  # Residual rename
  row.names(x)[which(row.names(x) == "units!units")] <- "VR"
  
  # Year Rename
  row.names(x)[which(row.names(x) == "NatalYear")] <- "VK"
  
  # Nest Rename
  row.names(x)[which(row.names(x) == "NatalNest")] <- "VN"
  
  # Permanent environment rename
  row.names(x)[which(row.names(x) == "ide(dad)")] <- "VPE"
  row.names(x)[which(row.names(x) == "ide(dad, ainv)")] <- "VPE"
  row.names(x)[which(row.names(x) == "ide(mom)")] <- "VME"
  row.names(x)[which(row.names(x) == "ide(mom, ainv)")] <- "VME"
  row.names(x)[which(row.names(x) == "d:ide(Help1)")] <- "VHE"
  row.names(x)[which(row.names(x) == "d:ide(Help1, ainv)")] <- "VHE"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:ide(Help1, ainv)+and(d:ide(Help2, ainv))+and(d:ide(Help3, ainv))+and(d:ide(Help4, ainv))+and(d:ide(Help5, ainv))!diag(2)_2")] <- "VHE"
  
  # Genetics rename
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)")] <- "VDG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:vm(Help1, ainv)+and(d:vm(Help2, ainv))+and(d:vm(Help3, ainv))+and(d:vm(Help4, ainv))+and(d:vm(Help5, ainv))!us(2)_1:1")] <- "VDG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:vm(Help1, ainv)+and(d:vm(Help2, ainv))+and(d:vm(Help3, ainv))+and(d:vm(Help4, ainv))+and(d:vm(Help5, ainv))!diag(2)_1")] <- "VDG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:ide(Help1, ainv)+and(d:ide(Help2, ainv))+and(d:ide(Help3, ainv))+and(d:ide(Help4, ainv))+and(d:ide(Help5, ainv))!diag(2)_1")] <- "VDG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(mom, ainv)!us(2)_1:1")] <- "VDG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(dad, ainv)!us(2)_1:1")] <- "VDG"
  
  # IGE rename
  row.names(x)[which(row.names(x) == "vm(mom, ainv)")] <- "VMG"
  row.names(x)[which(row.names(x) == "vm(dad, ainv)")] <- "VPG"
  row.names(x)[which(row.names(x) == "vm(Help1, ainv)")] <- "VHG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:vm(Help1, ainv)+and(d:vm(Help2, ainv))+and(d:vm(Help3, ainv))+and(d:vm(Help4, ainv))+and(d:vm(Help5, ainv))!us(2)_2:2")] <- "VHG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:vm(Help1, ainv)+and(d:vm(Help2, ainv))+and(d:vm(Help3, ainv))+and(d:vm(Help4, ainv))+and(d:vm(Help5, ainv))!diag(2)_2")] <- "VHG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(mom, ainv)!us(2)_2:2")] <- "VMG"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(dad, ainv)!us(2)_2:2")] <- "VPG"
  
  # Covariance Rename
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+d:vm(Help1, ainv)+and(d:vm(Help2, ainv))+and(d:vm(Help3, ainv))+and(d:vm(Help4, ainv))+and(d:vm(Help5, ainv))!us(2)_2:1")] <- "Cov"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(mom, ainv)!us(2)_2:1")] <- "Cov"
  row.names(x)[which(row.names(x) == "vm(USFWS, ainv)+vm(dad, ainv)!us(2)_2:1")] <- "Cov"
  return(x)
}

# Assign trait and model names to extract from testing lists
allmodels_traits <- c("nestweight", "juvweight","nestTarsus", "juvTarsus")
trait_names <- c("Nestling Weight", "Juvenile Weight", "Nestling Tarsus","Juvenile Tarsus")
names(trait_names) <- allmodels_traits

allmodels_Models <- c("mat","pat","help")
model_names <- c("maternal","paternal","helper")
names(model_names) <- allmodels_Models

# Extract all model variance components and match to model testing assessment

results <- data.frame()
for (j in 1:4) {
  trait <- allmodels_traits[j]
  for(k in 1:3){
    model <- allmodels_Models[k]
    working <- eval(parse(text=paste(trait, model,sep="_")))
    for(i in 1:(length(working)/2)){
      if(model == "mat"){
        df <- working[[10+i]]
      }else{
        df <- working[[5+i]]
      }
      df <- translator(df)
      df <- df %>% 
        mutate(component = as.numeric(df$component),
               std.error = ifelse(is.na(as.numeric(std.error)), "",as.numeric(df$std.error)),
               full = ifelse(is.na(as.numeric(std.error)),
                             round(df$component,3),
                             paste(round(df$component,3)," (", round(df$std.error,3),")",sep = ""))) %>%
        dplyr::select(full)
      df1 <- as.data.frame(t(df)) %>%
        mutate(socialPartner = ifelse(model == "mat" & i %in% 1:4, "none", model_names[model]),
               socialPartner = ifelse(model == "mat" & i == 10, "all", socialPartner),
               Trait = trait_names[trait],
               Model = as.character(ifelse(model == "mat", i, i+4)),
               Model = ifelse(Model == "8", "8a", Model),
               Model = ifelse(Model == "9", "8b", Model),
               Model = ifelse(Model == "10", "9", Model))
      results<- bind_rows(results,df1)
      rownames(results) <- NULL
    }
  }
}

# Final cleaning of table 2 to match naming conventions in main text, remove extra models and testing results

Table2 <- left_join(results, All_RE_sig %>%
                      mutate(deltaAIC = round(deltaAIC, 2)) %>%
                      dplyr::select(Trait, socialPartner, `Effect Tested` = Effect, Model,`Reference Model` =Reference_Model, LogLik, `LR-statistic`,`Pr(Chisq)`,AIC, deltaAIC)) %>%
  mutate(Trait = factor(Trait, levels = c("Nestling Weight", "Juvenile Weight",
                                          "Nestling Tarsus", "Juvenile Tarsus")),
         Model = ifelse(socialPartner == "maternal" & Model %in% c("5","6","7"), 
                        paste(Model," - M", sep = ""), Model),
         Model = ifelse(socialPartner == "paternal" & Model %in% c("5","6","7"), 
                        paste(Model," - D", sep = ""), Model),
         Model = ifelse(socialPartner == "helper" & Model %in% c("5","6","7"), 
                        paste(Model," - H", sep = ""), Model),
         Model = ifelse(socialPartner %in% c("maternal","paternal") & Model == "8a",
                        "M + D", Model),
         Model = ifelse(socialPartner %in% c("maternal") & Model == "8b",
                        "M + H", Model),
         Model = ifelse(socialPartner %in% c("paternal","helper") & Model == "8b",
                        "D + H", Model),
         Model = ifelse(socialPartner %in% c("helper") & Model == "8a",
                        "M + H", Model),
         Model = ifelse(Model == "9", "M + D + H", Model),
         Model = factor(Model, 
                        levels = c("1","2","3","4","5 - M","5 - D","5 - H",
                                   "6 - M","6 - D","6 - H",
                                   "7 - M","7 - D","7 - H",
                                   "M + D","M + H","D + H",
                                   "M + D + H")),
         across(c("LR-statistic", "Pr(Chisq)"),
                ~ ifelse(Model %in% c("M + D","M + H","D + H","M + D + H"),NA_real_, .))) %>%
  arrange(Trait,Model) %>%
  distinct(Trait, Model, .keep_all = T) %>%
  dplyr::select(Trait, Model,
                VDG, Cov, VMG, VPG, VHG, 
                VME, VPE, VHE,
                VK, VN, VR,
                LogLik,`LR-statistic`, `Pr(Chisq)`, 
                AIC, deltaAIC) %>%
  replace(is.na(.),"")


write_csv(Table2, "Table2_MainText.csv")

#### Figures

## Figure 1 - Variance components

varcomps_weight_plot <- ggplot(var_weight %>% filter(!(Model %in% c("Mom Cov","Dad Cov","Helper Cov"))),aes(y = Model, x = component, fill=VC_full)) +
  geom_bar(position = "fill",stat = "identity", size=1) +
  
  facet_wrap(~Trait_full, nrow=1)+
  coord_flip()+
  scale_fill_manual(values = c(viridis(5)[1:4],plasma(6)),
                    name = "Variance Component",
                    labels = rev(c(bquote(V[HE]),bquote(V[HG]),
                                   bquote(V[PE]),bquote(V[PG]),
                                   bquote(V[ME]),bquote(V[MG]),
                                   bquote(V[DG]),bquote(V[N]),bquote(V[K]),bquote(V[R])))) +
  labs(x="Proportion of Variance", y="") +
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1))

varcomps_tarsus_plot <- ggplot(var_tarsus%>% filter(!(Model %in% c("Mom Cov","Dad Cov","Helper Cov"))),aes(y = Model, x = component, fill=VC_full)) +
  geom_bar(position = "fill",stat = "identity", size=1) +
  
  facet_wrap(~Trait_full, nrow=1)+
  coord_flip()+
  scale_fill_manual(values = c(viridis(5)[1:4],plasma(6)),
                    name = "Variance Component",
                    labels = rev(c(bquote(V[HE]),bquote(V[HG]),
                                   bquote(V[PE]),bquote(V[PG]),
                                   bquote(V[ME]),bquote(V[MG]),
                                   bquote(V[DG]),bquote(V[N]),bquote(V[K]),bquote(V[R])))) +
  labs(x="Proportion of Variance", y="") +
  theme_classic() + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))

varcomps_all <- ggarrange(varcomps_weight_plot, varcomps_tarsus_plot, nrow = 1, align = "hv", common.legend = T, legend = "bottom")

ggsave("Figure1_MainText.png", plot = varcomps_all, width = 10, height =5)


## Figure 2 - Heritabilities

# Only including weight and tarsus in main text
heritability_params <- ggplot(Params_all%>% filter(!(Trait_full %in% c("Nestling Condition","Juvenile Condition")),
                                                   type != "TBV"),
                              # Plot estimates across each model colored by type of parameter
                              aes(x=Model,y=Estimate, col=type))+
  geom_point(position = position_dodge(0.9), size=2)+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE),
                width=0.25, size=0.5,
                position = position_dodge(0.9))+
  coord_cartesian(ylim=c(0,0.6))+
  labs(x="Trait",y="Heritabiltiy Estimate")+
  facet_wrap(~Trait_full, nrow = 2)+
  theme_bw() +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))+
  scale_color_manual(values = c("#35B779FF","#FCA636FF"), name="Parameter",labels = c(bquote(h^2),bquote(tau^2)))

ggsave("Figure2_MainText.png", plot = heritability_params, width = 6, height = 6)


#### Supplemental Tables and Figures

#### Tables

## Table s1 - Parameters

# Rename variables from all both parameter estimate tables

Tables1a <- Params_all %>% 
  dplyr::select(-Estimate, -SE, -fullsig,-sigmark,-Trait) %>%
  pivot_wider(names_from = type, values_from = full)

Tables1b <- totaleffects_all %>% 
  dplyr::select(-Estimate, -SE, -fullsig,-sigmark,-Trait,-type) %>%
  pivot_wider(names_from = `Indirect Effect`, values_from = full)

Tables1 <- left_join(Tables1a, Tables1b) %>%
  arrange(Trait_full, Model) %>%
  dplyr::select(Trait = Trait_full, Model, h2, tau, 
                Genetic, Environmental, Total)

write_csv(Tables1, "Suppl_Table1.csv")

## Table s2 - Dilution parameter without relatedness

# Helpers only

load("DilutionParameter-Helper-Tables.rdata")

# Adjust table for nestling weight
df1 <- d11.w_df %>%
  # calculate difference in AIC from best model
  mutate(deltaAIC = round(AIC - min(AIC),2),
         # Clean estimates and names
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Weight",
         `Social Partner` = "Helper") %>%
  arrange(AIC)%>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

# juvenile weight dilution
df2 <- juv.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

# nestling tarsus dilution
df3 <- d11.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Tarsus",
         `Social Partner` = "Helper") %>%
  arrange(AIC)%>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

# juvenile tarsus dilution
df4 <- juv.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

# Combine tables for overall helper dilution parameters
help_dilution <- bind_cols(bind_rows(df1,df3),bind_rows(df2,df4))

# Helpers and Mom

load("DilutionParameter-MomHelper-Tables.rdata")

df1 <- d11.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Weight",
         `Social Partner` = "Mom + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df2 <- juv.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

df3 <- d11.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Tarsus",
         `Social Partner` = "Mom + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df4 <- juv.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

mathelp_dilution <- bind_cols(bind_rows(df1,df3),bind_rows(df2,df4))

# Dad and helper

load("DilutionParameter-DadHelper-Tables.rdata")

df1 <- d11.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Weight",
         `Social Partner` = "Dad + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df2 <- juv.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

df3 <- d11.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Tarsus",
         `Social Partner` = "Dad + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df4 <- juv.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

pathelp_dilution <- bind_cols(bind_rows(df1,df3),bind_rows(df2,df4))

### All Three social partners

load("DilutionParameter-MomDadHelper-Tables.rdata")

df1 <- d11.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Weight",
         `Social Partner` = "Mom + Dad + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df2 <- juv.w_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

df3 <- d11.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL),
         Trait = "Tarsus",
         `Social Partner` = "Mom + Dad + Helper") %>%
  arrange(AIC) %>%
  dplyr::select(Trait, `Social Partner`, Model, LogL, AIC, deltaAIC)

df4 <- juv.t_df %>%
  mutate(deltaAIC = round(AIC - min(AIC),2),
         AIC = round(AIC, 2),
         LogL = round(LogL)) %>%
  arrange(AIC) %>%
  dplyr::select(LogL, AIC, deltaAIC)

matpathelp_dilution <- bind_cols(bind_rows(df1,df3),bind_rows(df2,df4)) 

# Combine different social partner tables and clean names

all_dilution <- bind_rows(help_dilution, pathelp_dilution, mathelp_dilution, matpathelp_dilution) %>%
  mutate(Relatedness = ifelse(grepl("_rel",Model),"y","n")) %>%
  select(Trait, `Social Partner`, Model, 
         Nest.LogLik = LogL...4, Nest.AIC = AIC...5, Nest.deltaAIC = deltaAIC...6,  
         Juv.LogLik = LogL...7, Juv.AIC = AIC...8, Juv.deltaAIC = deltaAIC...9, Relatedness)

# split overall table into no relatedness
dilution_norel <- all_dilution %>%
  filter(Relatedness == 'n') %>%
  select(-Relatedness)

write_csv(dilution_norel, "SupplTable2.csv")

## Table s3 - Dilution Parameter with relatedness

# split overall table into relatedness dilution
dilution_rel <- all_dilution %>%
  filter(Relatedness == 'y') %>%
  select(-Relatedness)

write_csv(dilution_rel, "SupplTable3.csv")

## Table s4 - Fixed Effect Estimates

# pull walds information
walds_out <- data.frame()
trts <- c("Nestling Weight", "Juvenile Weight",
          "Nestling Tarsus", "Juvenile Tarsus")
for(i in 1:length(walds)){
  df <- as.data.frame(walds[[i]]) %>%
    dplyr::select(P_val = Pr) %>%
    mutate(Trait = trts[i],
           FE = rownames(.),
           FE = ifelse(FE == "Sex", "Sex_M",FE))
  rownames(df) <- NULL
  walds_out <- bind_rows(walds_out, df)
}

# rename fixed effect estimates in list
fixed_effects <- data.frame()
names(fixed_list) <- c("Simple_Nestling Weight","Simple_Juvenile Weight",
                       "Simple_Nestling Tarsus","Simple_Juvenile Tarsus",
                       "Maternal_Nestling Weight","Maternal_Juvenile Weight",
                       "Maternal_Nestling Tarsus","Maternal_Juvenile Tarsus",
                       "Paternal_Nestling Weight","Paternal_Juvenile Weight",
                       "Paternal_Nestling Tarsus","Paternal_Juvenile Tarsus",
                       "Helper_Nestling Weight","Helper_Juvenile Weight",
                       "Helper_Nestling Tarsus","Helper_Juvenile Tarsus",
                       "MatPat_Nestling Weight","MatPat_Juvenile Weight",
                       "MatPat_Nestling Tarsus","MatPat_Juvenile Tarsus",
                       "MatHelper_Nestling Weight","MatHelper_Juvenile Weight",
                       "MatHelper_Nestling Tarsus","MatHelper_Juvenile Tarsus",
                       "PatHelper_Nestling Weight","PatHelper_Juvenile Weight",
                       "PatHelper_Nestling Tarsus","PatHelper_Juvenile Tarsus",
                       "MatPatHelper_Nestling Weight","MatPatHelper_Juvenile Weight",
                       "MatPatHelper_Nestling Tarsus","MatPatHelper_Juvenile Tarsus")

# pull fixed effect estimates from list and clean for making into a table
for(i in 1:length(fixed_list)){
  df <- as.data.frame(fixed_list[[i]])
  name <- str_split(names(fixed_list)[i], pattern="_")
  df <- df %>% 
    mutate(Trait = as.character(name[[1]][2]),
           FE = row.names(df),
           solution = as.numeric(df$solution),
           `std error` = as.numeric(df$`std error`)) %>%
    left_join(walds_out) %>%
    mutate(sigmark = ifelse(P_val < 0.05, " *", ""),
           full = paste(round(df$solution,3)," (", round(df$`std error`,3),")", sigmark,sep = ""))
  rownames(df) <- df$FE
  df <- df %>% dplyr::select(full)
  df1 <- as.data.frame(t(apply(df,2,rev)))
  df1$Model <- as.character(name[[1]][1])
  df1$Trait <- as.character(name[[1]][2])
  fixed_effects <- bind_rows(fixed_effects,df1)
}

# Create table and switch names to match text
tables4 <- fixed_effects %>%
  mutate(
    ## Rename models
    Model = ifelse(Model == "Maternal", "Mom", Model),
    Model = ifelse(Model == "Paternal", "Dad", Model),
    Model = ifelse(Model == "Helper", "Helper", Model),
    Model = ifelse(Model == "MatPat", "Mom + Dad", Model),
    Model = ifelse(Model == "MatHelper", "Mom + Helper", Model),
    Model = ifelse(Model == "PatHelper", "Dad + Helper", Model),
    Model = ifelse(Model == "MatPatHelper", "Mom + Dad + Helper", Model),
    Model = factor(Model, levels=c("DGE-Only","Mom","Dad","Helper","Mom Cov","Dad Cov","Helper Cov","Mom + Dad","Mom + Helper","Dad + Helper","Mom + Dad + Helper")),
    Trait = factor(Trait, levels = c("Nestling Weight", "Juvenile Weight",
                                     "Nestling Tarsus", "Juvenile Tarsus"))) %>%
  arrange(Trait, Model) %>%
  dplyr::select(Trait, Model,
                Intercept = `(Intercept)`, 
                `Sex (M)` = Sex_M, 
                `Hatching Date` = HatchDOYS,
                `Brood Size` = HatchNumS, 
                `Number of Helpers` = NumHelpS, 
                `Inbreeding Coefficient` = ICS, 
                `Age Measured` = AgeMeasS,
                `Territory Area` = TerrSizeS) %>%
  replace(is.na(.),"")

rownames(tables4) <- NULL

write_csv(tables4, "SupplTable4.csv")


#### Figures

## Figure s1 - Variance components with covariance models
varcomps_weight_plot <- ggplot(var_weight,aes(y = Model, x = component, fill=VC_full)) +
  geom_bar(position = "fill",stat = "identity", size=1) +
  
  facet_wrap(~Trait_full, nrow=1)+
  coord_flip()+
  scale_fill_manual(values = c(viridis(5)[1:4],plasma(6)),
                    name = "Variance Component",
                    labels = rev(c(bquote(V[HE]),bquote(V[HG]),
                                   bquote(V[PE]),bquote(V[PG]),
                                   bquote(V[ME]),bquote(V[MG]),
                                   bquote(V[DG]),bquote(V[N]),bquote(V[K]),bquote(V[R])))) +
  labs(x="Proportion of Variance", y="") +
  theme_classic() + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1))

varcomps_tarsus_plot <- ggplot(var_tarsus,aes(y = Model, x = component, fill=VC_full)) +
  geom_bar(position = "fill",stat = "identity", size=1) +
  
  facet_wrap(~Trait_full, nrow=1)+
  coord_flip()+
  scale_fill_manual(values = c(viridis(5)[1:4],plasma(6)),
                    name = "Variance Component",
                    labels = rev(c(bquote(V[HE]),bquote(V[HG]),
                                   bquote(V[PE]),bquote(V[PG]),
                                   bquote(V[ME]),bquote(V[MG]),
                                   bquote(V[DG]),bquote(V[N]),bquote(V[K]),bquote(V[R])))) +
  labs(x="Proportion of Variance", y="") +
  theme_classic() + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1))

varcomps_all <- ggarrange(varcomps_weight_plot, varcomps_tarsus_plot, nrow = 1, align = "hv", common.legend = T, legend = "bottom")

ggsave("Figure1_Supplement.png", plot = varcomps_all, width = 10, height =5)
