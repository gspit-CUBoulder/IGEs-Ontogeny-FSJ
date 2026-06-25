############# Significance Testing for Indirect Genetic Effects Across Ontogeny #####################

#### Script 3/5 in Florida Scrub Jay IGE Analysis
#### Original creator: Gladiana Spitz

#### Goal: This script tests statistical significance of all estimates. LRTs assess variance
#### component significance and Wald test fixed effect significance

#### Output: Wald results for each model. Hierarchical likelihood ratio tests across all models - (TestingResults.rdata)
#### Models and model variance components (TestingModels.rdata)


# Load Packages
library(asreml)
# Larger workspace makes this much easier
asreml.options(workspace="16gb") #if R crashes, increase this
#Filtered Data from 01-Cleaning.r
load("cleandata_FSJ_IGE_Analysis.RData") 

#Invert pedigree as usual
ainv <- ainverse(pedPrunedASREML)

# Inititalize lists and dataframes to catch results
LRTs <- list()
walds <- list()

# NC save log likelihoods and AIC/BIC values
# orig models
modelFit <- data.frame(matrix(ncol=6,nrow=80))
colnames(modelFit) <- c("trait", "socialPartner", "model", "LogLik", "AIC", "BIC")
modelFit$trait <- c(rep("d11.w",20),rep("juv.w",20),rep("d11.t",20),rep("juv.t",20))
modelFit$socialPartner <- c(rep("none",4),rep("maternal",5),"all",rep("paternal",5),rep("helper",5))
modelFit$model <- c(c("1","2","3","4","5","6","7","8a","8b","9"),rep(c("5","6","7","8a","8b"),2))

# General order of testing for each single IGE model
# m1 - no random effects
# m2 - only year
# m3 - year and nest
# m4 - year, nest, and additive genetic
# m5 - year, nest, additive genetic, and social environment for each social partner a -> maternal, b -> paternal, c -> helper
# m6 - year, nest, additive genetic, social environment, social genetic 
# m7 - year, nest, social environment, covariance structure between social genetic and additive genetic

# testing for multiple IGE models
# m8a - year, nest, additive genetic, social environment, social genetic, social 2 environment, social 2 genetic  
# m8b - year, nest, additive genetic, social environment, social genetic, social 3 environment, social 3 genetic  
# m9 - year, nest, additive genetic, social environment, social genetic, social 2 environment, social 2 genetic, social 3 environment, social 3 genetic  
# for mat: m8a = mat + pat, m8b = mat + helper; can compare m8a & m6
# for pat: m8a = mat + pat, m8b = pat + helper; can compare m8a & m6
# for helper: m8a = mat + helper, m8b = pat + helper; can compare m8a & m6, m8b & m6, m8ab & m9
# lrt.asreml(m6, m7, m8a, m9) and lrt.asreml(m6, m7, m8b, m9)

#### Nestling Weight Testing Models ####

# Random effect significance tests of maternal effect models on nestling weight

# Nestling weight dilution parameter
d11.w$GroupSize <- d11.w$NumHelp + 3
d <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^1

# models for testing
m1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit = 1000, data=d11.w)
m2 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m3 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m4 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+ide(mom,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+vm(mom,ainv)+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)

m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~str(~vm(USFWS,ainv)+vm(mom,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=d11.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)
m9 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
              vm(mom,ainv)+vm(dad,ainv)+
              ide(mom)+ide(dad)+
              d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
              NatalYear+NatalNest,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.w)

m1_VC <- summary(m1)$varcomp
m2_VC <- summary(m2)$varcomp
m3_VC <- summary(m3)$varcomp
m4_VC <- summary(m4)$varcomp
m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp
m9_VC <- summary(m9)$varcomp

# Save walds results
walds[["d11.w_m1"]] <- wald.asreml(m1, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m2"]] <- wald.asreml(m2, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m3"]] <- wald.asreml(m3, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m4"]] <- wald.asreml(m4, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m5_m"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m6_m"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m7_m"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m8a_mp"]] <- wald.asreml(m8a, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m8b_mh"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m9"]] <- wald.asreml(m9, ssType = 'conditional', denDF = 'numeric')$Wald

# Save results
LRTs[["m_d11.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["m_d11.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[1, 4:6] <- c(m1$loglik, summary(m1)$aic, summary(m1)$bic)
modelFit[2, 4:6] <- c(m2$loglik, summary(m2)$aic, summary(m2)$bic)
modelFit[3, 4:6] <- c(m3$loglik, summary(m3)$aic, summary(m3)$bic)
modelFit[4, 4:6] <- c(m4$loglik, summary(m4)$aic, summary(m4)$bic)
modelFit[5, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[6, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[7, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[8, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[9, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)
modelFit[10, 4:6] <- c(m9$loglik, summary(m9)$aic, summary(m9)$bic)

# NC save model output 
nestweight_mat <- list(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
                       
                       m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


# Random effect significance tests of paternal effect models on nestling weight

m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+ide(dad,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+vm(dad,ainv)+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)
m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~str(~vm(USFWS,ainv)+vm(dad,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=d11.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random =  ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

# Save walds results
walds[["d11.w_m5_p"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m6_p"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m7_p"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m8b_ph"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald

# Save results
LRTs[["p_d11.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["p_d11.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[11, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[12, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[13, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[14, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[15, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
nestweight_pat <- list(m5, m6, m7, m8a, m8b,
                       m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

# Random effect significance tests of helper effect models on nestling weight
m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:ide(Help1,ainv)+and(d:ide(Help2,ainv))+and(d:ide(Help3,ainv))+and(d:ide(Help4,ainv))+and(d:ide(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)
m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random =  ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.w)

# Save walds results
walds[["d11.w_m5_h"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m6_h"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.w_m7_h"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald

# save variance components
m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

LRTs[["h_d11.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["h_d11.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)
# NC save logL/AIC/BIC
modelFit[16, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[17, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[18, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[19, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[20, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
nestweight_help <- list(m5, m6, m7, m8a, m8b,
                        m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

rm(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
   m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


#### Juvenile Weight Testing Models ####

# Random effect significance tests of maternal effect models on juvenile weight

# Juvenile weight dilution parameter
juv.w$GroupSize <- juv.w$NumHelp + 3
d <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^1

m1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit = 1000, data=juv.w)
m2 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m3 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m4 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+ide(mom,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+vm(mom,ainv)+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~str(~vm(USFWS,ainv)+vm(mom,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=juv.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)
m9 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+vm(dad,ainv)+
               ide(mom)+ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)

m1_VC <- summary(m1)$varcomp
m2_VC <- summary(m2)$varcomp
m3_VC <- summary(m3)$varcomp
m4_VC <- summary(m4)$varcomp
m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp
m9_VC <- summary(m9)$varcomp

# Save walds results
walds[["juv.w_m1"]] <- wald.asreml(m1, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m2"]] <- wald.asreml(m2, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m3"]] <- wald.asreml(m3, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m4"]] <- wald.asreml(m4, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m5_m"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m6_m"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m7_m"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m8a_mp"]] <- wald.asreml(m8a, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m8b_mh"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m9"]] <- wald.asreml(m9, ssType = 'conditional', denDF = 'numeric')$Wald

LRTs[["m_juv.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["m_juv.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC 
modelFit[21, 4:6] <- c(m1$loglik, summary(m1)$aic, summary(m1)$bic)
modelFit[22, 4:6] <- c(m2$loglik, summary(m2)$aic, summary(m2)$bic)
modelFit[23, 4:6] <- c(m3$loglik, summary(m3)$aic, summary(m3)$bic)
modelFit[24, 4:6] <- c(m4$loglik, summary(m4)$aic, summary(m4)$bic)
modelFit[25, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[26, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[27, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[28, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[29, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)
modelFit[30, 4:6] <- c(m9$loglik, summary(m9)$aic, summary(m9)$bic)

# NC save model output 
juvweight_mat <- list(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
                      m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


# Random effect significance tests of paternal effect models on juvenile weight
m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+ide(dad,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+vm(dad,ainv)+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~str(~vm(USFWS,ainv)+vm(dad,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=juv.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)


# Save walds results
walds[["juv.w_m5_p"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m6_p"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m7_p"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m8b_ph"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

# Save results
LRTs[["p_juv.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["p_juv.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[31, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[32, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[33, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[34, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[35, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
juvweight_pat <- list(m5, m6, m7, m8a, m8b,
                      m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

# Random effect significance tests of helper effect models on juvenile weight
m5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:ide(Help1,ainv)+and(d:ide(Help2,ainv))+and(d:ide(Help3,ainv))+and(d:ide(Help4,ainv))+and(d:ide(Help5,ainv)),
                           ~diag(2):vm(USFWS, ainv))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)
m6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                           ~diag(2):vm(USFWS, ainv))+NatalYear+NatalNest+d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5)),
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.w)
m7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                           ~us(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)

m8a <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)
m8b <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.w)

# Save walds results
walds[["juv.w_m5_h"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m6_h"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.w_m7_h"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

LRTs[["h_juv.w_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["h_juv.w_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[36, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[37, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[38, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[39, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[40, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
juvweight_help <- list(m5, m6, m7, m8a, m8b,
                      m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

rm(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
   m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


#### Nestling Tarsus Model Testing ####

# nestling tarsus dilution parameter
d11.t$GroupSize <- d11.t$NumHelp + 3
d <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0

# Random effect significance tests of maternal effect models on nestling Tarsus

m1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit = 1000, data=d11.t)
m2 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m3 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m4 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+ide(mom,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+vm(mom,ainv)+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~str(~vm(USFWS,ainv)+vm(mom,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=d11.t)
m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)
m9 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
              vm(mom,ainv)+vm(dad,ainv)+
              ide(mom)+ide(dad)+
              d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
              NatalYear+NatalNest,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.t)

# Save walds results
walds[["d11.t_m1"]] <- wald.asreml(m1, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m2"]] <- wald.asreml(m2, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m3"]] <- wald.asreml(m3, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m4"]] <- wald.asreml(m4, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m5_m"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m6_m"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m7_m"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m8a_mp"]] <- wald.asreml(m8a, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m8b_mh"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m9"]] <- wald.asreml(m9, ssType = 'conditional', denDF = 'numeric')$Wald

m1_VC <- summary(m1)$varcomp
m2_VC <- summary(m2)$varcomp
m3_VC <- summary(m3)$varcomp
m4_VC <- summary(m4)$varcomp
m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp
m9_VC <- summary(m9)$varcomp

LRTs[["m_d11.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["m_d11.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC - need to edit
modelFit[41, 4:6] <- c(m1$loglik, summary(m1)$aic, summary(m1)$bic)
modelFit[42, 4:6] <- c(m2$loglik, summary(m2)$aic, summary(m2)$bic)
modelFit[43, 4:6] <- c(m3$loglik, summary(m3)$aic, summary(m3)$bic)
modelFit[44, 4:6] <- c(m4$loglik, summary(m4)$aic, summary(m4)$bic)
modelFit[45, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[46, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[47, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[48, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[49, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)
modelFit[50, 4:6] <- c(m9$loglik, summary(m9)$aic, summary(m9)$bic)

# NC save model output 
nestTarsus_mat <- list(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
                       m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)

# Random effect significance tests of paternal effect models on nestling Tarsus
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+ide(dad,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~ vm(USFWS,ainv)+vm(dad,ainv)+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,
             random = ~str(~vm(USFWS,ainv)+vm(dad,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=d11.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=d11.t)
m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)

# Save walds results
walds[["d11.t_m5_p"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m6_p"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m7_p"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m8b_ph"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

# Save results

LRTs[["p_d11.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["p_d11.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[51, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[52, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[53, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[54, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[55, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
nestTarsus_pat <- list(m5, m6, m7, m8a, m8b,
                       m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

# Random effect significance tests of helper effect models on nestling Tarsus
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:ide(Help1,ainv)+and(d:ide(Help2,ainv))+and(d:ide(Help3,ainv))+and(d:ide(Help4,ainv))+and(d:ide(Help5,ainv)),
                           ~diag(2):vm(USFWS, ainv))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                           ~diag(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                           ~us(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)
m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
             random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=d11.t)

# Save walds results
walds[["d11.t_m5_h"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m6_h"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["d11.t_m7_h"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

LRTs[["h_d11.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["h_d11.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[56, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[57, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[58, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[59, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[60, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
nestTarsus_help <- list(m5, m6, m7, m8a, m8b,
                        m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)


rm(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
   m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


#### Juvenile Tarsus Testing Models ####

# juvenile tarsus dilution parameter

juv.t$GroupSize <- juv.t$NumHelp + 3
d <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^1

# Random effect significance tests of maternal effect models on juvenile Tarsus

m1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit = 1000, data=juv.t)
m2 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m3 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m4 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+ide(mom,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+vm(mom,ainv)+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~str(~vm(USFWS,ainv)+vm(mom,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(mom),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=juv.t)
m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)
m9 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
            random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
              vm(mom,ainv)+vm(dad,ainv)+
              ide(mom)+ide(dad)+
              d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
              NatalYear+NatalNest,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=juv.t)

m1_VC <- summary(m1)$varcomp
m2_VC <- summary(m2)$varcomp
m3_VC <- summary(m3)$varcomp
m4_VC <- summary(m4)$varcomp
m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp
m9_VC <- summary(m9)$varcomp

# Save walds results
walds[["juv.t_m1"]] <- wald.asreml(m1, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m2"]] <- wald.asreml(m2, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m3"]] <- wald.asreml(m3, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m4"]] <- wald.asreml(m4, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m5_m"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m6_m"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m7_m"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m8a_mp"]] <- wald.asreml(m8a, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m8b_mh"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m9"]] <- wald.asreml(m9, ssType = 'conditional', denDF = 'numeric')$Wald

LRTs[["m_juv.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["m_juv.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC - need to edit
modelFit[61, 4:6] <- c(m1$loglik, summary(m1)$aic, summary(m1)$bic)
modelFit[62, 4:6] <- c(m2$loglik, summary(m2)$aic, summary(m2)$bic)
modelFit[63, 4:6] <- c(m3$loglik, summary(m3)$aic, summary(m3)$bic)
modelFit[64, 4:6] <- c(m4$loglik, summary(m4)$aic, summary(m4)$bic)
modelFit[65, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[66, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[67, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[68, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[69, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)
modelFit[70, 4:6] <- c(m9$loglik, summary(m9)$aic, summary(m9)$bic)

# NC save model output 
juvTarsus_mat <- list(m1, m2, m3, m4, m5, m6, m7, m8a, m8b, m9,
                      m1_VC, m2_VC, m3_VC, m4_VC, m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC, m9_VC)


# Random effect significance tests of paternal effect models on juvenile Tarsus
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+ide(dad,ainv)+NatalYear+NatalNest,
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~ vm(USFWS,ainv)+vm(dad,ainv)+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,
             random = ~str(~vm(USFWS,ainv)+vm(dad,ainv),~us(2):vm(USFWS,ainv))+NatalYear+NatalNest+ide(dad),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit = 1000, data=juv.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
             residual=~idv(units), na.action = na.method(x="include", y="include"),		
             maxit=1000, data=juv.t)
m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

# Save walds results
walds[["juv.t_m5_p"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m6_p"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m7_p"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m8b_ph"]] <- wald.asreml(m8b, ssType = 'conditional', denDF = 'numeric')$Wald

LRTs[["p_juv.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["p_juv.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[71, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[72, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[73, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[74, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[75, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
juvTarsus_pat <- list(m5, m6, m7, m8a, m8b,
                      m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

# Random effect significance tests of helper effect models on juvenile Tarsus
m5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:ide(Help1,ainv)+and(d:ide(Help2,ainv))+and(d:ide(Help3,ainv))+and(d:ide(Help4,ainv))+and(d:ide(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)
m6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)
m7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)

m8a <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(mom,ainv)+
               ide(mom)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)

m8b <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
             random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
               vm(dad,ainv)+
               ide(dad)+
               d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
               NatalYear+NatalNest,
             equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
             residual=~idv(units), na.action = na.method(x="include", y="include"),
             maxit=1000,
             data=juv.t)

m5_VC <- summary(m5)$varcomp
m6_VC <- summary(m6)$varcomp
m7_VC <- summary(m7)$varcomp
m8a_VC <- summary(m8a)$varcomp
m8b_VC <- summary(m8b)$varcomp

# Save walds results
walds[["juv.t_m5_h"]] <- wald.asreml(m5, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m6_h"]] <- wald.asreml(m6, ssType = 'conditional', denDF = 'numeric')$Wald
walds[["juv.t_m7_h"]] <- wald.asreml(m7, ssType = 'conditional', denDF = 'numeric')$Wald

LRTs[["h_juv.t_1"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8a, m9)
LRTs[["h_juv.t_2"]] <- lrt.asreml(m1, m2, m3, m4, m5, m6, m7, m8b, m9)

# NC save logL/AIC/BIC
modelFit[76, 4:6] <- c(m5$loglik, summary(m5)$aic, summary(m5)$bic)
modelFit[77, 4:6] <- c(m6$loglik, summary(m6)$aic, summary(m6)$bic)
modelFit[78, 4:6] <- c(m7$loglik, summary(m7)$aic, summary(m7)$bic)
modelFit[79, 4:6] <- c(m8a$loglik, summary(m8a)$aic, summary(m8a)$bic)
modelFit[80, 4:6] <- c(m8b$loglik, summary(m8b)$aic, summary(m8b)$bic)

# NC save model output 
juvTarsus_help <- list(m5, m6, m7, m8a, m8b,
                       m5_VC, m6_VC, m7_VC, m8a_VC, m8b_VC)

# save model objects

save(nestweight_mat, juvweight_mat, nestweight_pat, juvweight_pat, nestweight_help, juvweight_help,
     nestTarsus_mat, juvTarsus_mat, nestTarsus_pat, juvTarsus_pat, nestTarsus_help, juvTarsus_help, file = "TestingModels.Rdata")

# save testing objects
save(LRTs, walds, modelFit,file="TestingResults.Rdata")