##### Estimating Best Dilution Parameter for Indirect Genetic Effects Across Ontogeny #####

#### Script 2c/5 in Florida Scrub Jay IGE Analysis
#### Original creator: Gladiana Spitz

#### Goal: This script tests the best group size dilution parameter for models that include 
#### up to 5 helper social partners and dads. Diluting the effect of each helper by the size of the group
#### gives a more accurate estimate of variation in juvenile metrics that occurs due to differences in helpers

#### Output: Model results for each dilution parameter (DilutionParameter-DadHelper-Models.rdata) and 
#### data tables of model fit results (DilutionParameter-DadHelper-Tables.rdata) for comparison and selection of the appropriate dilution parameter


# Load ASReml-R
library(asreml)

# Larger workspace makes this much easier
asreml.options(workspace="8gb") #if R crashes, increase this
#Filtered Data from 01-Cleaning.r
load("cleandata_FSJ_IGE_Analysis.RData") 
#Invert pedigree as usual
ainv <- ainverse(pedPrunedASREML)


#### 0 Nestling Weight ####
d11.w$GroupSize <- d11.w$NumHelp + 3

# group size vector
w0 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0

# each helper relatedness vector
w0_h1 <- w0 * d11.w$Help1_rel
w0_h2 <- w0 * d11.w$Help2_rel
w0_h3 <- w0 * d11.w$Help3_rel
w0_h4 <- w0 * d11.w$Help4_rel
w0_h5 <- w0 * d11.w$Help5_rel


# group size model
d0 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS, ainv)+w0:vm(Help1,ainv)+and(w0:vm(Help2,ainv))+and(w0:vm(Help3,ainv))+and(w0:vm(Help4,ainv))+and(w0:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w0:ide(Help1)+and(w0:ide(Help2))+and(w0:ide(Help3))+and(w0:ide(Help4))+and(w0:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.w)

# relatedness and group size model
d0_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS, ainv)+w0_h1:vm(Help1,ainv)+and(w0_h2:vm(Help2,ainv))+and(w0_h3:vm(Help3,ainv))+and(w0_h4:vm(Help4,ainv))+and(w0_h5:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w0_h1:ide(Help1)+and(w0_h2:ide(Help2))+and(w0_h3:ide(Help3))+and(w0_h4:ide(Help4))+and(w0_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.w)

#### 0.1 Nestling Weight ####
w0.1 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.1

w0.1_h1 <- w0.1 * d11.w$Help1_rel
w0.1_h2 <- w0.1 * d11.w$Help2_rel
w0.1_h3 <- w0.1 * d11.w$Help3_rel
w0.1_h4 <- w0.1 * d11.w$Help4_rel
w0.1_h5 <- w0.1 * d11.w$Help5_rel


d0.1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.1:vm(Help1,ainv)+and(w0.1:vm(Help2,ainv))+and(w0.1:vm(Help3,ainv))+and(w0.1:vm(Help4,ainv))+and(w0.1:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.1:ide(Help1)+and(w0.1:ide(Help2))+and(w0.1:ide(Help3))+and(w0.1:ide(Help4))+and(w0.1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.1_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                random = ~str(~vm(USFWS, ainv)+w0.1_h1:vm(Help1,ainv)+and(w0.1_h2:vm(Help2,ainv))+and(w0.1_h3:vm(Help3,ainv))+and(w0.1_h4:vm(Help4,ainv))+and(w0.1_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w0.1_h1:ide(Help1)+and(w0.1_h2:ide(Help2))+and(w0.1_h3:ide(Help3))+and(w0.1_h4:ide(Help4))+and(w0.1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=d11.w)

#### 0.2 Nestling Weight ####
w0.2 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.2

w0.2_h1 <- w0.2 * d11.w$Help1_rel
w0.2_h2 <- w0.2 * d11.w$Help2_rel
w0.2_h3 <- w0.2 * d11.w$Help3_rel
w0.2_h4 <- w0.2 * d11.w$Help4_rel
w0.2_h5 <- w0.2 * d11.w$Help5_rel


d0.2 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.2:vm(Help1,ainv)+and(w0.2:vm(Help2,ainv))+and(w0.2:vm(Help3,ainv))+and(w0.2:vm(Help4,ainv))+and(w0.2:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.2:ide(Help1)+and(w0.2:ide(Help2))+and(w0.2:ide(Help3))+and(w0.2:ide(Help4))+and(w0.2:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.2_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.2_h1:vm(Help1,ainv)+and(w0.2_h2:vm(Help2,ainv))+and(w0.2_h3:vm(Help3,ainv))+and(w0.2_h4:vm(Help4,ainv))+and(w0.2_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.2_h1:ide(Help1)+and(w0.2_h2:ide(Help2))+and(w0.2_h3:ide(Help3))+and(w0.2_h4:ide(Help4))+and(w0.2_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.3 Nestling Weight ####
w0.3 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.3

w0.3_h1 <- w0.3 * d11.w$Help1_rel
w0.3_h2 <- w0.3 * d11.w$Help2_rel
w0.3_h3 <- w0.3 * d11.w$Help3_rel
w0.3_h4 <- w0.3 * d11.w$Help4_rel
w0.3_h5 <- w0.3 * d11.w$Help5_rel


d0.3 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.3:vm(Help1,ainv)+and(w0.3:vm(Help2,ainv))+and(w0.3:vm(Help3,ainv))+and(w0.3:vm(Help4,ainv))+and(w0.3:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.3:ide(Help1)+and(w0.3:ide(Help2))+and(w0.3:ide(Help3))+and(w0.3:ide(Help4))+and(w0.3:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.3_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.3_h1:vm(Help1,ainv)+and(w0.3_h2:vm(Help2,ainv))+and(w0.3_h3:vm(Help3,ainv))+and(w0.3_h4:vm(Help4,ainv))+and(w0.3_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.3_h1:ide(Help1)+and(w0.3_h2:ide(Help2))+and(w0.3_h3:ide(Help3))+and(w0.3_h4:ide(Help4))+and(w0.3_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.4 Nestling Weight ####
w0.4 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.4

w0.4_h1 <- w0.4 * d11.w$Help1_rel
w0.4_h2 <- w0.4 * d11.w$Help2_rel
w0.4_h3 <- w0.4 * d11.w$Help3_rel
w0.4_h4 <- w0.4 * d11.w$Help4_rel
w0.4_h5 <- w0.4 * d11.w$Help5_rel


d0.4 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.4:vm(Help1,ainv)+and(w0.4:vm(Help2,ainv))+and(w0.4:vm(Help3,ainv))+and(w0.4:vm(Help4,ainv))+and(w0.4:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.4:ide(Help1)+and(w0.4:ide(Help2))+and(w0.4:ide(Help3))+and(w0.4:ide(Help4))+and(w0.4:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.4_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.4_h1:vm(Help1,ainv)+and(w0.4_h2:vm(Help2,ainv))+and(w0.4_h3:vm(Help3,ainv))+and(w0.4_h4:vm(Help4,ainv))+and(w0.4_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.4_h1:ide(Help1)+and(w0.4_h2:ide(Help2))+and(w0.4_h3:ide(Help3))+and(w0.4_h4:ide(Help4))+and(w0.4_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.5 Nestling Weight ####
w0.5 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.5

w0.5_h1 <- w0.5 * d11.w$Help1_rel
w0.5_h2 <- w0.5 * d11.w$Help2_rel
w0.5_h3 <- w0.5 * d11.w$Help3_rel
w0.5_h4 <- w0.5 * d11.w$Help4_rel
w0.5_h5 <- w0.5 * d11.w$Help5_rel


d0.5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.5:vm(Help1,ainv)+and(w0.5:vm(Help2,ainv))+and(w0.5:vm(Help3,ainv))+and(w0.5:vm(Help4,ainv))+and(w0.5:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.5:ide(Help1)+and(w0.5:ide(Help2))+and(w0.5:ide(Help3))+and(w0.5:ide(Help4))+and(w0.5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.5_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.5_h1:vm(Help1,ainv)+and(w0.5_h2:vm(Help2,ainv))+and(w0.5_h3:vm(Help3,ainv))+and(w0.5_h4:vm(Help4,ainv))+and(w0.5_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.5_h1:ide(Help1)+and(w0.5_h2:ide(Help2))+and(w0.5_h3:ide(Help3))+and(w0.5_h4:ide(Help4))+and(w0.5_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.6 Nestling Weight ####
w0.6 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.6

w0.6_h1 <- w0.6 * d11.w$Help1_rel
w0.6_h2 <- w0.6 * d11.w$Help2_rel
w0.6_h3 <- w0.6 * d11.w$Help3_rel
w0.6_h4 <- w0.6 * d11.w$Help4_rel
w0.6_h5 <- w0.6 * d11.w$Help5_rel


d0.6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
               random = ~str(~vm(USFWS, ainv)+w0.6:vm(Help1,ainv)+and(w0.6:vm(Help2,ainv))+and(w0.6:vm(Help3,ainv))+and(w0.6:vm(Help4,ainv))+and(w0.6:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.6:ide(Help1)+and(w0.6:ide(Help2))+and(w0.6:ide(Help3))+and(w0.6:ide(Help4))+and(w0.6:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=d11.w)

d0.6_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.6_h1:vm(Help1,ainv)+and(w0.6_h2:vm(Help2,ainv))+and(w0.6_h3:vm(Help3,ainv))+and(w0.6_h4:vm(Help4,ainv))+and(w0.6_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.6_h1:ide(Help1)+and(w0.6_h2:ide(Help2))+and(w0.6_h3:ide(Help3))+and(w0.6_h4:ide(Help4))+and(w0.6_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.7 Nestling Weight ####
w0.7 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.7

w0.7_h1 <- w0.7 * d11.w$Help1_rel
w0.7_h2 <- w0.7 * d11.w$Help2_rel
w0.7_h3 <- w0.7 * d11.w$Help3_rel
w0.7_h4 <- w0.7 * d11.w$Help4_rel
w0.7_h5 <- w0.7 * d11.w$Help5_rel


d0.7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
               random = ~str(~vm(USFWS, ainv)+w0.7:vm(Help1,ainv)+and(w0.7:vm(Help2,ainv))+and(w0.7:vm(Help3,ainv))+and(w0.7:vm(Help4,ainv))+and(w0.7:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.7:ide(Help1)+and(w0.7:ide(Help2))+and(w0.7:ide(Help3))+and(w0.7:ide(Help4))+and(w0.7:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=d11.w)

d0.7_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.7_h1:vm(Help1,ainv)+and(w0.7_h2:vm(Help2,ainv))+and(w0.7_h3:vm(Help3,ainv))+and(w0.7_h4:vm(Help4,ainv))+and(w0.7_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.7_h1:ide(Help1)+and(w0.7_h2:ide(Help2))+and(w0.7_h3:ide(Help3))+and(w0.7_h4:ide(Help4))+and(w0.7_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.8 Nestling Weight ####
w0.8 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.8

w0.8_h1 <- w0.8 * d11.w$Help1_rel
w0.8_h2 <- w0.8 * d11.w$Help2_rel
w0.8_h3 <- w0.8 * d11.w$Help3_rel
w0.8_h4 <- w0.8 * d11.w$Help4_rel
w0.8_h5 <- w0.8 * d11.w$Help5_rel


d0.8 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.8:vm(Help1,ainv)+and(w0.8:vm(Help2,ainv))+and(w0.8:vm(Help3,ainv))+and(w0.8:vm(Help4,ainv))+and(w0.8:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.8:ide(Help1)+and(w0.8:ide(Help2))+and(w0.8:ide(Help3))+and(w0.8:ide(Help4))+and(w0.8:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)


d0.8_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.8_h1:vm(Help1,ainv)+and(w0.8_h2:vm(Help2,ainv))+and(w0.8_h3:vm(Help3,ainv))+and(w0.8_h4:vm(Help4,ainv))+and(w0.8_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.8_h1:ide(Help1)+and(w0.8_h2:ide(Help2))+and(w0.8_h3:ide(Help3))+and(w0.8_h4:ide(Help4))+and(w0.8_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 0.9 Nestling Weight ####
w0.9 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^0.9

w0.9_h1 <- w0.9 * d11.w$Help1_rel
w0.9_h2 <- w0.9 * d11.w$Help2_rel
w0.9_h3 <- w0.9 * d11.w$Help3_rel
w0.9_h4 <- w0.9 * d11.w$Help4_rel
w0.9_h5 <- w0.9 * d11.w$Help5_rel


d0.9 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.9:vm(Help1,ainv)+and(w0.9:vm(Help2,ainv))+and(w0.9:vm(Help3,ainv))+and(w0.9:vm(Help4,ainv))+and(w0.9:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.9:ide(Help1)+and(w0.9:ide(Help2))+and(w0.9:ide(Help3))+and(w0.9:ide(Help4))+and(w0.9:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.w)

d0.9_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.9_h1:vm(Help1,ainv)+and(w0.9_h2:vm(Help2,ainv))+and(w0.9_h3:vm(Help3,ainv))+and(w0.9_h4:vm(Help4,ainv))+and(w0.9_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.9_h1:ide(Help1)+and(w0.9_h2:ide(Help2))+and(w0.9_h3:ide(Help3))+and(w0.9_h4:ide(Help4))+and(w0.9_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### 1 Nestling Weight ####
w1 <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^1

w1_h1 <- w1 * d11.w$Help1_rel
w1_h2 <- w1 * d11.w$Help2_rel
w1_h3 <- w1 * d11.w$Help3_rel
w1_h4 <- w1 * d11.w$Help4_rel
w1_h5 <- w1 * d11.w$Help5_rel



d1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS, ainv)+w1:vm(Help1,ainv)+and(w1:vm(Help2,ainv))+and(w1:vm(Help3,ainv))+and(w1:vm(Help4,ainv))+and(w1:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w1:ide(Help1)+and(w1:ide(Help2))+and(w1:ide(Help3))+and(w1:ide(Help4))+and(w1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.w)

d1_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w1_h1:vm(Help1,ainv)+and(w1_h2:vm(Help2,ainv))+and(w1_h3:vm(Help3,ainv))+and(w1_h4:vm(Help4,ainv))+and(w1_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w1_h1:ide(Help1)+and(w1_h2:ide(Help2))+and(w1_h3:ide(Help3))+and(w1_h4:ide(Help4))+and(w1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.w)

#### Save Nestling Weight Output ####

d11.w_df <- data.frame(Model = c("d0","d0.1","d0.2","d0.3","d0.4","d0.5","d0.6","d0.7","d0.8","d0.9","d1",
                                 "d0_rel","d0.1_rel","d0.2_rel","d0.3_rel","d0.4_rel","d0.5_rel","d0.6_rel","d0.7_rel","d0.8_rel","d0.9_rel","d1_rel"), 
LogL = c(d0$loglik,d0.1$loglik,d0.2$loglik,d0.3$loglik,d0.4$loglik,d0.5$loglik,d0.6$loglik,d0.7$loglik,d0.8$loglik,d0.9$loglik,d1$loglik,
         d0_rel$loglik,d0.1_rel$loglik,d0.2_rel$loglik,d0.3_rel$loglik,d0.4_rel$loglik,d0.5_rel$loglik,d0.6_rel$loglik,d0.7_rel$loglik,d0.8_rel$loglik,d0.9_rel$loglik,d1_rel$loglik),
AIC = c(summary(d0)$aic,summary(d0.1)$aic,summary(d0.2)$aic,summary(d0.3)$aic,summary(d0.4)$aic,summary(d0.5)$aic,summary(d0.6)$aic,summary(d0.7)$aic,summary(d0.8)$aic,summary(d0.9)$aic,summary(d1)$aic,
        summary(d0_rel)$aic,summary(d0.1_rel)$aic,summary(d0.2_rel)$aic,summary(d0.3_rel)$aic,summary(d0.4_rel)$aic,summary(d0.5_rel)$aic,summary(d0.6_rel)$aic,summary(d0.7_rel)$aic,summary(d0.8_rel)$aic,summary(d0.9_rel)$aic,summary(d1_rel)$aic),
BIC = c(summary(d0)$bic,summary(d0.1)$bic,summary(d0.2)$bic,summary(d0.3)$bic,summary(d0.4)$bic,summary(d0.5)$bic,summary(d0.6)$bic,summary(d0.7)$bic,summary(d0.8)$bic,summary(d0.9)$bic,summary(d1)$bic,
        summary(d0_rel)$bic,summary(d0.1_rel)$bic,summary(d0.2_rel)$bic,summary(d0.3_rel)$bic,summary(d0.4_rel)$bic,summary(d0.5_rel)$bic,summary(d0.6_rel)$bic,summary(d0.7_rel)$bic,summary(d0.8_rel)$bic,summary(d0.9_rel)$bic,summary(d1_rel)$bic))

d11.w_models <- list(d0,d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9,d1,
                     d0_rel,d0.1_rel,d0.2_rel,d0.3_rel,d0.4_rel,d0.5_rel,d0.6_rel,d0.7_rel,d0.8_rel,d0.9_rel,d1_rel)
#### 0 Juvenile Weight ####


juv.w$GroupSize <- juv.w$NumHelp + 3

# group size vector
w0 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0

# each helper relatedness vector
w0_h1 <- w0 * juv.w$Help1_rel
w0_h2 <- w0 * juv.w$Help2_rel
w0_h3 <- w0 * juv.w$Help3_rel
w0_h4 <- w0 * juv.w$Help4_rel
w0_h5 <- w0 * juv.w$Help5_rel


# group size model
d0 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
            random = ~str(~vm(USFWS, ainv)+w0:vm(Help1,ainv)+and(w0:vm(Help2,ainv))+and(w0:vm(Help3,ainv))+and(w0:vm(Help4,ainv))+and(w0:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w0:ide(Help1)+and(w0:ide(Help2))+and(w0:ide(Help3))+and(w0:ide(Help4))+and(w0:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=juv.w)

# relatedness and group size model
d0_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                random = ~str(~vm(USFWS, ainv)+w0_h1:vm(Help1,ainv)+and(w0_h2:vm(Help2,ainv))+and(w0_h3:vm(Help3,ainv))+and(w0_h4:vm(Help4,ainv))+and(w0_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w0_h1:ide(Help1)+and(w0_h2:ide(Help2))+and(w0_h3:ide(Help3))+and(w0_h4:ide(Help4))+and(w0_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=juv.w)

#### 0.1 Juvenile Weight ####
w0.1 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.1

w0.1_h1 <- w0.1 * juv.w$Help1_rel
w0.1_h2 <- w0.1 * juv.w$Help2_rel
w0.1_h3 <- w0.1 * juv.w$Help3_rel
w0.1_h4 <- w0.1 * juv.w$Help4_rel
w0.1_h5 <- w0.1 * juv.w$Help5_rel


d0.1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.1:vm(Help1,ainv)+and(w0.1:vm(Help2,ainv))+and(w0.1:vm(Help3,ainv))+and(w0.1:vm(Help4,ainv))+and(w0.1:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.1:ide(Help1)+and(w0.1:ide(Help2))+and(w0.1:ide(Help3))+and(w0.1:ide(Help4))+and(w0.1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)

d0.1_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.1_h1:vm(Help1,ainv)+and(w0.1_h2:vm(Help2,ainv))+and(w0.1_h3:vm(Help3,ainv))+and(w0.1_h4:vm(Help4,ainv))+and(w0.1_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.1_h1:ide(Help1)+and(w0.1_h2:ide(Help2))+and(w0.1_h3:ide(Help3))+and(w0.1_h4:ide(Help4))+and(w0.1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.2 Juvenile Weight ####
w0.2 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.2

w0.2_h1 <- w0.2 * juv.w$Help1_rel
w0.2_h2 <- w0.2 * juv.w$Help2_rel
w0.2_h3 <- w0.2 * juv.w$Help3_rel
w0.2_h4 <- w0.2 * juv.w$Help4_rel
w0.2_h5 <- w0.2 * juv.w$Help5_rel


d0.2 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.2:vm(Help1,ainv)+and(w0.2:vm(Help2,ainv))+and(w0.2:vm(Help3,ainv))+and(w0.2:vm(Help4,ainv))+and(w0.2:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.2:ide(Help1)+and(w0.2:ide(Help2))+and(w0.2:ide(Help3))+and(w0.2:ide(Help4))+and(w0.2:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)

d0.2_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.2_h1:vm(Help1,ainv)+and(w0.2_h2:vm(Help2,ainv))+and(w0.2_h3:vm(Help3,ainv))+and(w0.2_h4:vm(Help4,ainv))+and(w0.2_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.2_h1:ide(Help1)+and(w0.2_h2:ide(Help2))+and(w0.2_h3:ide(Help3))+and(w0.2_h4:ide(Help4))+and(w0.2_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.3 Juvenile Weight ####
w0.3 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.3

w0.3_h1 <- w0.3 * juv.w$Help1_rel
w0.3_h2 <- w0.3 * juv.w$Help2_rel
w0.3_h3 <- w0.3 * juv.w$Help3_rel
w0.3_h4 <- w0.3 * juv.w$Help4_rel
w0.3_h5 <- w0.3 * juv.w$Help5_rel


d0.3 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.3:vm(Help1,ainv)+and(w0.3:vm(Help2,ainv))+and(w0.3:vm(Help3,ainv))+and(w0.3:vm(Help4,ainv))+and(w0.3:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.3:ide(Help1)+and(w0.3:ide(Help2))+and(w0.3:ide(Help3))+and(w0.3:ide(Help4))+and(w0.3:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)

d0.3_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.3_h1:vm(Help1,ainv)+and(w0.3_h2:vm(Help2,ainv))+and(w0.3_h3:vm(Help3,ainv))+and(w0.3_h4:vm(Help4,ainv))+and(w0.3_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.3_h1:ide(Help1)+and(w0.3_h2:ide(Help2))+and(w0.3_h3:ide(Help3))+and(w0.3_h4:ide(Help4))+and(w0.3_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.4 Juvenile Weight ####
w0.4 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.4

w0.4_h1 <- w0.4 * juv.w$Help1_rel
w0.4_h2 <- w0.4 * juv.w$Help2_rel
w0.4_h3 <- w0.4 * juv.w$Help3_rel
w0.4_h4 <- w0.4 * juv.w$Help4_rel
w0.4_h5 <- w0.4 * juv.w$Help5_rel


d0.4 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.4:vm(Help1,ainv)+and(w0.4:vm(Help2,ainv))+and(w0.4:vm(Help3,ainv))+and(w0.4:vm(Help4,ainv))+and(w0.4:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.4:ide(Help1)+and(w0.4:ide(Help2))+and(w0.4:ide(Help3))+and(w0.4:ide(Help4))+and(w0.4:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)

d0.4_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.4_h1:vm(Help1,ainv)+and(w0.4_h2:vm(Help2,ainv))+and(w0.4_h3:vm(Help3,ainv))+and(w0.4_h4:vm(Help4,ainv))+and(w0.4_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.4_h1:ide(Help1)+and(w0.4_h2:ide(Help2))+and(w0.4_h3:ide(Help3))+and(w0.4_h4:ide(Help4))+and(w0.4_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.5 Juvenile Weight ####
w0.5 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.5

w0.5_h1 <- w0.5 * juv.w$Help1_rel
w0.5_h2 <- w0.5 * juv.w$Help2_rel
w0.5_h3 <- w0.5 * juv.w$Help3_rel
w0.5_h4 <- w0.5 * juv.w$Help4_rel
w0.5_h5 <- w0.5 * juv.w$Help5_rel


d0.5 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.5:vm(Help1,ainv)+and(w0.5:vm(Help2,ainv))+and(w0.5:vm(Help3,ainv))+and(w0.5:vm(Help4,ainv))+and(w0.5:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.5:ide(Help1)+and(w0.5:ide(Help2))+and(w0.5:ide(Help3))+and(w0.5:ide(Help4))+and(w0.5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.w)

d0.5_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.5_h1:vm(Help1,ainv)+and(w0.5_h2:vm(Help2,ainv))+and(w0.5_h3:vm(Help3,ainv))+and(w0.5_h4:vm(Help4,ainv))+and(w0.5_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.5_h1:ide(Help1)+and(w0.5_h2:ide(Help2))+and(w0.5_h3:ide(Help3))+and(w0.5_h4:ide(Help4))+and(w0.5_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.6 Juvenile Weight ####
w0.6 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.6

w0.6_h1 <- w0.6 * juv.w$Help1_rel
w0.6_h2 <- w0.6 * juv.w$Help2_rel
w0.6_h3 <- w0.6 * juv.w$Help3_rel
w0.6_h4 <- w0.6 * juv.w$Help4_rel
w0.6_h5 <- w0.6 * juv.w$Help5_rel


d0.6 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.6:vm(Help1,ainv)+and(w0.6:vm(Help2,ainv))+and(w0.6:vm(Help3,ainv))+and(w0.6:vm(Help4,ainv))+and(w0.6:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.6:ide(Help1)+and(w0.6:ide(Help2))+and(w0.6:ide(Help3))+and(w0.6:ide(Help4))+and(w0.6:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.w)

d0.6_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.6_h1:vm(Help1,ainv)+and(w0.6_h2:vm(Help2,ainv))+and(w0.6_h3:vm(Help3,ainv))+and(w0.6_h4:vm(Help4,ainv))+and(w0.6_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.6_h1:ide(Help1)+and(w0.6_h2:ide(Help2))+and(w0.6_h3:ide(Help3))+and(w0.6_h4:ide(Help4))+and(w0.6_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.7 Juvenile Weight ####
w0.7 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.7

w0.7_h1 <- w0.7 * juv.w$Help1_rel
w0.7_h2 <- w0.7 * juv.w$Help2_rel
w0.7_h3 <- w0.7 * juv.w$Help3_rel
w0.7_h4 <- w0.7 * juv.w$Help4_rel
w0.7_h5 <- w0.7 * juv.w$Help5_rel


d0.7 <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.7:vm(Help1,ainv)+and(w0.7:vm(Help2,ainv))+and(w0.7:vm(Help3,ainv))+and(w0.7:vm(Help4,ainv))+and(w0.7:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.7:ide(Help1)+and(w0.7:ide(Help2))+and(w0.7:ide(Help3))+and(w0.7:ide(Help4))+and(w0.7:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.w)

d0.7_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.7_h1:vm(Help1,ainv)+and(w0.7_h2:vm(Help2,ainv))+and(w0.7_h3:vm(Help3,ainv))+and(w0.7_h4:vm(Help4,ainv))+and(w0.7_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.7_h1:ide(Help1)+and(w0.7_h2:ide(Help2))+and(w0.7_h3:ide(Help3))+and(w0.7_h4:ide(Help4))+and(w0.7_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.8 Juvenile Weight ####
w0.8 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.8

w0.8_h1 <- w0.8 * juv.w$Help1_rel
w0.8_h2 <- w0.8 * juv.w$Help2_rel
w0.8_h3 <- w0.8 * juv.w$Help3_rel
w0.8_h4 <- w0.8 * juv.w$Help4_rel
w0.8_h5 <- w0.8 * juv.w$Help5_rel


d0.8 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.8:vm(Help1,ainv)+and(w0.8:vm(Help2,ainv))+and(w0.8:vm(Help3,ainv))+and(w0.8:vm(Help4,ainv))+and(w0.8:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.8:ide(Help1)+and(w0.8:ide(Help2))+and(w0.8:ide(Help3))+and(w0.8:ide(Help4))+and(w0.8:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)


d0.8_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.8_h1:vm(Help1,ainv)+and(w0.8_h2:vm(Help2,ainv))+and(w0.8_h3:vm(Help3,ainv))+and(w0.8_h4:vm(Help4,ainv))+and(w0.8_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.8_h1:ide(Help1)+and(w0.8_h2:ide(Help2))+and(w0.8_h3:ide(Help3))+and(w0.8_h4:ide(Help4))+and(w0.8_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 0.9 Juvenile Weight ####
w0.9 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^0.9

w0.9_h1 <- w0.9 * juv.w$Help1_rel
w0.9_h2 <- w0.9 * juv.w$Help2_rel
w0.9_h3 <- w0.9 * juv.w$Help3_rel
w0.9_h4 <- w0.9 * juv.w$Help4_rel
w0.9_h5 <- w0.9 * juv.w$Help5_rel


d0.9 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.9:vm(Help1,ainv)+and(w0.9:vm(Help2,ainv))+and(w0.9:vm(Help3,ainv))+and(w0.9:vm(Help4,ainv))+and(w0.9:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.9:ide(Help1)+and(w0.9:ide(Help2))+and(w0.9:ide(Help3))+and(w0.9:ide(Help4))+and(w0.9:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.w)

d0.9_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.9_h1:vm(Help1,ainv)+and(w0.9_h2:vm(Help2,ainv))+and(w0.9_h3:vm(Help3,ainv))+and(w0.9_h4:vm(Help4,ainv))+and(w0.9_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.9_h1:ide(Help1)+and(w0.9_h2:ide(Help2))+and(w0.9_h3:ide(Help3))+and(w0.9_h4:ide(Help4))+and(w0.9_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.w)

#### 1 Juvenile Weight ####
w1 <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^1

w1_h1 <- w1 * juv.w$Help1_rel
w1_h2 <- w1 * juv.w$Help2_rel
w1_h3 <- w1 * juv.w$Help3_rel
w1_h4 <- w1 * juv.w$Help4_rel
w1_h5 <- w1 * juv.w$Help5_rel



d1 <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
            random = ~str(~vm(USFWS, ainv)+w1:vm(Help1,ainv)+and(w1:vm(Help2,ainv))+and(w1:vm(Help3,ainv))+and(w1:vm(Help4,ainv))+and(w1:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w1:ide(Help1)+and(w1:ide(Help2))+and(w1:ide(Help3))+and(w1:ide(Help4))+and(w1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=juv.w)

d1_rel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                random = ~str(~vm(USFWS, ainv)+w1_h1:vm(Help1,ainv)+and(w1_h2:vm(Help2,ainv))+and(w1_h3:vm(Help3,ainv))+and(w1_h4:vm(Help4,ainv))+and(w1_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w1_h1:ide(Help1)+and(w1_h2:ide(Help2))+and(w1_h3:ide(Help3))+and(w1_h4:ide(Help4))+and(w1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=juv.w)

#### Save Juvenile Weight Output ####

juv.w_df <- data.frame(Model = c("d0","d0.1","d0.2","d0.3","d0.4","d0.5","d0.6","d0.7","d0.8","d0.9","d1",
                                 "d0_rel","d0.1_rel","d0.2_rel","d0.3_rel","d0.4_rel","d0.5_rel","d0.6_rel","d0.7_rel","d0.8_rel","d0.9_rel","d1_rel"), 
                       LogL = c(d0$loglik,d0.1$loglik,d0.2$loglik,d0.3$loglik,d0.4$loglik,d0.5$loglik,d0.6$loglik,d0.7$loglik,d0.8$loglik,d0.9$loglik,d1$loglik,
                                d0_rel$loglik,d0.1_rel$loglik,d0.2_rel$loglik,d0.3_rel$loglik,d0.4_rel$loglik,d0.5_rel$loglik,d0.6_rel$loglik,d0.7_rel$loglik,d0.8_rel$loglik,d0.9_rel$loglik,d1_rel$loglik),
                       AIC = c(summary(d0)$aic,summary(d0.1)$aic,summary(d0.2)$aic,summary(d0.3)$aic,summary(d0.4)$aic,summary(d0.5)$aic,summary(d0.6)$aic,summary(d0.7)$aic,summary(d0.8)$aic,summary(d0.9)$aic,summary(d1)$aic,
                               summary(d0_rel)$aic,summary(d0.1_rel)$aic,summary(d0.2_rel)$aic,summary(d0.3_rel)$aic,summary(d0.4_rel)$aic,summary(d0.5_rel)$aic,summary(d0.6_rel)$aic,summary(d0.7_rel)$aic,summary(d0.8_rel)$aic,summary(d0.9_rel)$aic,summary(d1_rel)$aic),
                       BIC = c(summary(d0)$bic,summary(d0.1)$bic,summary(d0.2)$bic,summary(d0.3)$bic,summary(d0.4)$bic,summary(d0.5)$bic,summary(d0.6)$bic,summary(d0.7)$bic,summary(d0.8)$bic,summary(d0.9)$bic,summary(d1)$bic,
                               summary(d0_rel)$bic,summary(d0.1_rel)$bic,summary(d0.2_rel)$bic,summary(d0.3_rel)$bic,summary(d0.4_rel)$bic,summary(d0.5_rel)$bic,summary(d0.6_rel)$bic,summary(d0.7_rel)$bic,summary(d0.8_rel)$bic,summary(d0.9_rel)$bic,summary(d1_rel)$bic))

juv.w_models <- list(d0,d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9,d1,
                     d0_rel,d0.1_rel,d0.2_rel,d0.3_rel,d0.4_rel,d0.5_rel,d0.6_rel,d0.7_rel,d0.8_rel,d0.9_rel,d1_rel)



#### 0 Nestling Tarsus ####
d11.t$GroupSize <- d11.t$NumHelp + 3

# group size vector
w0 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0

# each helper relatedness vector
w0_h1 <- w0 * d11.t$Help1_rel
w0_h2 <- w0 * d11.t$Help2_rel
w0_h3 <- w0 * d11.t$Help3_rel
w0_h4 <- w0 * d11.t$Help4_rel
w0_h5 <- w0 * d11.t$Help5_rel


# group size model
d0 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS, ainv)+w0:vm(Help1,ainv)+and(w0:vm(Help2,ainv))+and(w0:vm(Help3,ainv))+and(w0:vm(Help4,ainv))+and(w0:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w0:ide(Help1)+and(w0:ide(Help2))+and(w0:ide(Help3))+and(w0:ide(Help4))+and(w0:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.t)

# relatedness and group size model
d0_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                random = ~str(~vm(USFWS, ainv)+w0_h1:vm(Help1,ainv)+and(w0_h2:vm(Help2,ainv))+and(w0_h3:vm(Help3,ainv))+and(w0_h4:vm(Help4,ainv))+and(w0_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w0_h1:ide(Help1)+and(w0_h2:ide(Help2))+and(w0_h3:ide(Help3))+and(w0_h4:ide(Help4))+and(w0_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=d11.t)

#### 0.1 Nestling Tarsus ####
w0.1 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.1

w0.1_h1 <- w0.1 * d11.t$Help1_rel
w0.1_h2 <- w0.1 * d11.t$Help2_rel
w0.1_h3 <- w0.1 * d11.t$Help3_rel
w0.1_h4 <- w0.1 * d11.t$Help4_rel
w0.1_h5 <- w0.1 * d11.t$Help5_rel


d0.1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.1:vm(Help1,ainv)+and(w0.1:vm(Help2,ainv))+and(w0.1:vm(Help3,ainv))+and(w0.1:vm(Help4,ainv))+and(w0.1:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.1:ide(Help1)+and(w0.1:ide(Help2))+and(w0.1:ide(Help3))+and(w0.1:ide(Help4))+and(w0.1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)

d0.1_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.1_h1:vm(Help1,ainv)+and(w0.1_h2:vm(Help2,ainv))+and(w0.1_h3:vm(Help3,ainv))+and(w0.1_h4:vm(Help4,ainv))+and(w0.1_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.1_h1:ide(Help1)+and(w0.1_h2:ide(Help2))+and(w0.1_h3:ide(Help3))+and(w0.1_h4:ide(Help4))+and(w0.1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.2 Nestling Tarsus ####
w0.2 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.2

w0.2_h1 <- w0.2 * d11.t$Help1_rel
w0.2_h2 <- w0.2 * d11.t$Help2_rel
w0.2_h3 <- w0.2 * d11.t$Help3_rel
w0.2_h4 <- w0.2 * d11.t$Help4_rel
w0.2_h5 <- w0.2 * d11.t$Help5_rel


d0.2 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.2:vm(Help1,ainv)+and(w0.2:vm(Help2,ainv))+and(w0.2:vm(Help3,ainv))+and(w0.2:vm(Help4,ainv))+and(w0.2:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.2:ide(Help1)+and(w0.2:ide(Help2))+and(w0.2:ide(Help3))+and(w0.2:ide(Help4))+and(w0.2:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)

d0.2_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.2_h1:vm(Help1,ainv)+and(w0.2_h2:vm(Help2,ainv))+and(w0.2_h3:vm(Help3,ainv))+and(w0.2_h4:vm(Help4,ainv))+and(w0.2_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.2_h1:ide(Help1)+and(w0.2_h2:ide(Help2))+and(w0.2_h3:ide(Help3))+and(w0.2_h4:ide(Help4))+and(w0.2_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.3 Nestling Tarsus ####
w0.3 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.3

w0.3_h1 <- w0.3 * d11.t$Help1_rel
w0.3_h2 <- w0.3 * d11.t$Help2_rel
w0.3_h3 <- w0.3 * d11.t$Help3_rel
w0.3_h4 <- w0.3 * d11.t$Help4_rel
w0.3_h5 <- w0.3 * d11.t$Help5_rel


d0.3 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.3:vm(Help1,ainv)+and(w0.3:vm(Help2,ainv))+and(w0.3:vm(Help3,ainv))+and(w0.3:vm(Help4,ainv))+and(w0.3:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.3:ide(Help1)+and(w0.3:ide(Help2))+and(w0.3:ide(Help3))+and(w0.3:ide(Help4))+and(w0.3:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)

d0.3_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.3_h1:vm(Help1,ainv)+and(w0.3_h2:vm(Help2,ainv))+and(w0.3_h3:vm(Help3,ainv))+and(w0.3_h4:vm(Help4,ainv))+and(w0.3_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.3_h1:ide(Help1)+and(w0.3_h2:ide(Help2))+and(w0.3_h3:ide(Help3))+and(w0.3_h4:ide(Help4))+and(w0.3_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.4 Nestling Tarsus ####
w0.4 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.4

w0.4_h1 <- w0.4 * d11.t$Help1_rel
w0.4_h2 <- w0.4 * d11.t$Help2_rel
w0.4_h3 <- w0.4 * d11.t$Help3_rel
w0.4_h4 <- w0.4 * d11.t$Help4_rel
w0.4_h5 <- w0.4 * d11.t$Help5_rel


d0.4 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.4:vm(Help1,ainv)+and(w0.4:vm(Help2,ainv))+and(w0.4:vm(Help3,ainv))+and(w0.4:vm(Help4,ainv))+and(w0.4:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.4:ide(Help1)+and(w0.4:ide(Help2))+and(w0.4:ide(Help3))+and(w0.4:ide(Help4))+and(w0.4:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)

d0.4_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.4_h1:vm(Help1,ainv)+and(w0.4_h2:vm(Help2,ainv))+and(w0.4_h3:vm(Help3,ainv))+and(w0.4_h4:vm(Help4,ainv))+and(w0.4_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.4_h1:ide(Help1)+and(w0.4_h2:ide(Help2))+and(w0.4_h3:ide(Help3))+and(w0.4_h4:ide(Help4))+and(w0.4_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.5 Nestling Tarsus ####
w0.5 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.5

w0.5_h1 <- w0.5 * d11.t$Help1_rel
w0.5_h2 <- w0.5 * d11.t$Help2_rel
w0.5_h3 <- w0.5 * d11.t$Help3_rel
w0.5_h4 <- w0.5 * d11.t$Help4_rel
w0.5_h5 <- w0.5 * d11.t$Help5_rel


d0.5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
               random = ~str(~vm(USFWS, ainv)+w0.5:vm(Help1,ainv)+and(w0.5:vm(Help2,ainv))+and(w0.5:vm(Help3,ainv))+and(w0.5:vm(Help4,ainv))+and(w0.5:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.5:ide(Help1)+and(w0.5:ide(Help2))+and(w0.5:ide(Help3))+and(w0.5:ide(Help4))+and(w0.5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=d11.t)

d0.5_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.5_h1:vm(Help1,ainv)+and(w0.5_h2:vm(Help2,ainv))+and(w0.5_h3:vm(Help3,ainv))+and(w0.5_h4:vm(Help4,ainv))+and(w0.5_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.5_h1:ide(Help1)+and(w0.5_h2:ide(Help2))+and(w0.5_h3:ide(Help3))+and(w0.5_h4:ide(Help4))+and(w0.5_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.6 Nestling Tarsus ####
w0.6 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.6

w0.6_h1 <- w0.6 * d11.t$Help1_rel
w0.6_h2 <- w0.6 * d11.t$Help2_rel
w0.6_h3 <- w0.6 * d11.t$Help3_rel
w0.6_h4 <- w0.6 * d11.t$Help4_rel
w0.6_h5 <- w0.6 * d11.t$Help5_rel


d0.6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
               random = ~str(~vm(USFWS, ainv)+w0.6:vm(Help1,ainv)+and(w0.6:vm(Help2,ainv))+and(w0.6:vm(Help3,ainv))+and(w0.6:vm(Help4,ainv))+and(w0.6:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.6:ide(Help1)+and(w0.6:ide(Help2))+and(w0.6:ide(Help3))+and(w0.6:ide(Help4))+and(w0.6:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=d11.t)

d0.6_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.6_h1:vm(Help1,ainv)+and(w0.6_h2:vm(Help2,ainv))+and(w0.6_h3:vm(Help3,ainv))+and(w0.6_h4:vm(Help4,ainv))+and(w0.6_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.6_h1:ide(Help1)+and(w0.6_h2:ide(Help2))+and(w0.6_h3:ide(Help3))+and(w0.6_h4:ide(Help4))+and(w0.6_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.7 Nestling Tarsus ####
w0.7 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.7

w0.7_h1 <- w0.7 * d11.t$Help1_rel
w0.7_h2 <- w0.7 * d11.t$Help2_rel
w0.7_h3 <- w0.7 * d11.t$Help3_rel
w0.7_h4 <- w0.7 * d11.t$Help4_rel
w0.7_h5 <- w0.7 * d11.t$Help5_rel


d0.7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
               random = ~str(~vm(USFWS, ainv)+w0.7:vm(Help1,ainv)+and(w0.7:vm(Help2,ainv))+and(w0.7:vm(Help3,ainv))+and(w0.7:vm(Help4,ainv))+and(w0.7:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.7:ide(Help1)+and(w0.7:ide(Help2))+and(w0.7:ide(Help3))+and(w0.7:ide(Help4))+and(w0.7:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=d11.t)

d0.7_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.7_h1:vm(Help1,ainv)+and(w0.7_h2:vm(Help2,ainv))+and(w0.7_h3:vm(Help3,ainv))+and(w0.7_h4:vm(Help4,ainv))+and(w0.7_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.7_h1:ide(Help1)+and(w0.7_h2:ide(Help2))+and(w0.7_h3:ide(Help3))+and(w0.7_h4:ide(Help4))+and(w0.7_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.8 Nestling Tarsus ####
w0.8 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.8

w0.8_h1 <- w0.8 * d11.t$Help1_rel
w0.8_h2 <- w0.8 * d11.t$Help2_rel
w0.8_h3 <- w0.8 * d11.t$Help3_rel
w0.8_h4 <- w0.8 * d11.t$Help4_rel
w0.8_h5 <- w0.8 * d11.t$Help5_rel


d0.8 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.8:vm(Help1,ainv)+and(w0.8:vm(Help2,ainv))+and(w0.8:vm(Help3,ainv))+and(w0.8:vm(Help4,ainv))+and(w0.8:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.8:ide(Help1)+and(w0.8:ide(Help2))+and(w0.8:ide(Help3))+and(w0.8:ide(Help4))+and(w0.8:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)


d0.8_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.8_h1:vm(Help1,ainv)+and(w0.8_h2:vm(Help2,ainv))+and(w0.8_h3:vm(Help3,ainv))+and(w0.8_h4:vm(Help4,ainv))+and(w0.8_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.8_h1:ide(Help1)+and(w0.8_h2:ide(Help2))+and(w0.8_h3:ide(Help3))+and(w0.8_h4:ide(Help4))+and(w0.8_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 0.9 Nestling Tarsus ####
w0.9 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0.9

w0.9_h1 <- w0.9 * d11.t$Help1_rel
w0.9_h2 <- w0.9 * d11.t$Help2_rel
w0.9_h3 <- w0.9 * d11.t$Help3_rel
w0.9_h4 <- w0.9 * d11.t$Help4_rel
w0.9_h5 <- w0.9 * d11.t$Help5_rel


d0.9 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
              random = ~str(~vm(USFWS, ainv)+w0.9:vm(Help1,ainv)+and(w0.9:vm(Help2,ainv))+and(w0.9:vm(Help3,ainv))+and(w0.9:vm(Help4,ainv))+and(w0.9:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.9:ide(Help1)+and(w0.9:ide(Help2))+and(w0.9:ide(Help3))+and(w0.9:ide(Help4))+and(w0.9:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=d11.t)

d0.9_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                  random = ~str(~vm(USFWS, ainv)+w0.9_h1:vm(Help1,ainv)+and(w0.9_h2:vm(Help2,ainv))+and(w0.9_h3:vm(Help3,ainv))+and(w0.9_h4:vm(Help4,ainv))+and(w0.9_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.9_h1:ide(Help1)+and(w0.9_h2:ide(Help2))+and(w0.9_h3:ide(Help3))+and(w0.9_h4:ide(Help4))+and(w0.9_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=d11.t)

#### 1 Nestling Tarsus ####
w1 <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^1

w1_h1 <- w1 * d11.t$Help1_rel
w1_h2 <- w1 * d11.t$Help2_rel
w1_h3 <- w1 * d11.t$Help3_rel
w1_h4 <- w1 * d11.t$Help4_rel
w1_h5 <- w1 * d11.t$Help5_rel



d1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
            random = ~str(~vm(USFWS, ainv)+w1:vm(Help1,ainv)+and(w1:vm(Help2,ainv))+and(w1:vm(Help3,ainv))+and(w1:vm(Help4,ainv))+and(w1:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w1:ide(Help1)+and(w1:ide(Help2))+and(w1:ide(Help3))+and(w1:ide(Help4))+and(w1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=d11.t)

d1_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                random = ~str(~vm(USFWS, ainv)+w1_h1:vm(Help1,ainv)+and(w1_h2:vm(Help2,ainv))+and(w1_h3:vm(Help3,ainv))+and(w1_h4:vm(Help4,ainv))+and(w1_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w1_h1:ide(Help1)+and(w1_h2:ide(Help2))+and(w1_h3:ide(Help3))+and(w1_h4:ide(Help4))+and(w1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=d11.t)

#### Save Nestling Tarsus Output ####

d11.t_df <- data.frame(Model = c("d0","d0.1","d0.2","d0.3","d0.4","d0.5","d0.6","d0.7","d0.8","d0.9","d1",
                                 "d0_rel","d0.1_rel","d0.2_rel","d0.3_rel","d0.4_rel","d0.5_rel","d0.6_rel","d0.7_rel","d0.8_rel","d0.9_rel","d1_rel"), 
                       LogL = c(d0$loglik,d0.1$loglik,d0.2$loglik,d0.3$loglik,d0.4$loglik,d0.5$loglik,d0.6$loglik,d0.7$loglik,d0.8$loglik,d0.9$loglik,d1$loglik,
                                d0_rel$loglik,d0.1_rel$loglik,d0.2_rel$loglik,d0.3_rel$loglik,d0.4_rel$loglik,d0.5_rel$loglik,d0.6_rel$loglik,d0.7_rel$loglik,d0.8_rel$loglik,d0.9_rel$loglik,d1_rel$loglik),
                       AIC = c(summary(d0)$aic,summary(d0.1)$aic,summary(d0.2)$aic,summary(d0.3)$aic,summary(d0.4)$aic,summary(d0.5)$aic,summary(d0.6)$aic,summary(d0.7)$aic,summary(d0.8)$aic,summary(d0.9)$aic,summary(d1)$aic,
                               summary(d0_rel)$aic,summary(d0.1_rel)$aic,summary(d0.2_rel)$aic,summary(d0.3_rel)$aic,summary(d0.4_rel)$aic,summary(d0.5_rel)$aic,summary(d0.6_rel)$aic,summary(d0.7_rel)$aic,summary(d0.8_rel)$aic,summary(d0.9_rel)$aic,summary(d1_rel)$aic),
                       BIC = c(summary(d0)$bic,summary(d0.1)$bic,summary(d0.2)$bic,summary(d0.3)$bic,summary(d0.4)$bic,summary(d0.5)$bic,summary(d0.6)$bic,summary(d0.7)$bic,summary(d0.8)$bic,summary(d0.9)$bic,summary(d1)$bic,
                               summary(d0_rel)$bic,summary(d0.1_rel)$bic,summary(d0.2_rel)$bic,summary(d0.3_rel)$bic,summary(d0.4_rel)$bic,summary(d0.5_rel)$bic,summary(d0.6_rel)$bic,summary(d0.7_rel)$bic,summary(d0.8_rel)$bic,summary(d0.9_rel)$bic,summary(d1_rel)$bic))

d11.t_models <- list(d0,d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9,d1,
                     d0_rel,d0.1_rel,d0.2_rel,d0.3_rel,d0.4_rel,d0.5_rel,d0.6_rel,d0.7_rel,d0.8_rel,d0.9_rel,d1_rel)
#### 0 Juvenile Tarsus ####


juv.t$GroupSize <- juv.t$NumHelp + 3

# group size vector
w0 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0

# each helper relatedness vector
w0_h1 <- w0 * juv.t$Help1_rel
w0_h2 <- w0 * juv.t$Help2_rel
w0_h3 <- w0 * juv.t$Help3_rel
w0_h4 <- w0 * juv.t$Help4_rel
w0_h5 <- w0 * juv.t$Help5_rel


# group size model
d0 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
            random = ~str(~vm(USFWS, ainv)+w0:vm(Help1,ainv)+and(w0:vm(Help2,ainv))+and(w0:vm(Help3,ainv))+and(w0:vm(Help4,ainv))+and(w0:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w0:ide(Help1)+and(w0:ide(Help2))+and(w0:ide(Help3))+and(w0:ide(Help4))+and(w0:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=juv.t)

# relatedness and group size model
d0_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                random = ~str(~vm(USFWS, ainv)+w0_h1:vm(Help1,ainv)+and(w0_h2:vm(Help2,ainv))+and(w0_h3:vm(Help3,ainv))+and(w0_h4:vm(Help4,ainv))+and(w0_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w0_h1:ide(Help1)+and(w0_h2:ide(Help2))+and(w0_h3:ide(Help3))+and(w0_h4:ide(Help4))+and(w0_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=juv.t)

#### 0.1 Juvenile Tarsus ####
w0.1 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.1

w0.1_h1 <- w0.1 * juv.t$Help1_rel
w0.1_h2 <- w0.1 * juv.t$Help2_rel
w0.1_h3 <- w0.1 * juv.t$Help3_rel
w0.1_h4 <- w0.1 * juv.t$Help4_rel
w0.1_h5 <- w0.1 * juv.t$Help5_rel


d0.1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.1:vm(Help1,ainv)+and(w0.1:vm(Help2,ainv))+and(w0.1:vm(Help3,ainv))+and(w0.1:vm(Help4,ainv))+and(w0.1:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.1:ide(Help1)+and(w0.1:ide(Help2))+and(w0.1:ide(Help3))+and(w0.1:ide(Help4))+and(w0.1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)

d0.1_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.1_h1:vm(Help1,ainv)+and(w0.1_h2:vm(Help2,ainv))+and(w0.1_h3:vm(Help3,ainv))+and(w0.1_h4:vm(Help4,ainv))+and(w0.1_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.1_h1:ide(Help1)+and(w0.1_h2:ide(Help2))+and(w0.1_h3:ide(Help3))+and(w0.1_h4:ide(Help4))+and(w0.1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.2 Juvenile Tarsus ####
w0.2 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.2

w0.2_h1 <- w0.2 * juv.t$Help1_rel
w0.2_h2 <- w0.2 * juv.t$Help2_rel
w0.2_h3 <- w0.2 * juv.t$Help3_rel
w0.2_h4 <- w0.2 * juv.t$Help4_rel
w0.2_h5 <- w0.2 * juv.t$Help5_rel


d0.2 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.2:vm(Help1,ainv)+and(w0.2:vm(Help2,ainv))+and(w0.2:vm(Help3,ainv))+and(w0.2:vm(Help4,ainv))+and(w0.2:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.2:ide(Help1)+and(w0.2:ide(Help2))+and(w0.2:ide(Help3))+and(w0.2:ide(Help4))+and(w0.2:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)

d0.2_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.2_h1:vm(Help1,ainv)+and(w0.2_h2:vm(Help2,ainv))+and(w0.2_h3:vm(Help3,ainv))+and(w0.2_h4:vm(Help4,ainv))+and(w0.2_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.2_h1:ide(Help1)+and(w0.2_h2:ide(Help2))+and(w0.2_h3:ide(Help3))+and(w0.2_h4:ide(Help4))+and(w0.2_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.3 Juvenile Tarsus ####
w0.3 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.3

w0.3_h1 <- w0.3 * juv.t$Help1_rel
w0.3_h2 <- w0.3 * juv.t$Help2_rel
w0.3_h3 <- w0.3 * juv.t$Help3_rel
w0.3_h4 <- w0.3 * juv.t$Help4_rel
w0.3_h5 <- w0.3 * juv.t$Help5_rel


d0.3 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.3:vm(Help1,ainv)+and(w0.3:vm(Help2,ainv))+and(w0.3:vm(Help3,ainv))+and(w0.3:vm(Help4,ainv))+and(w0.3:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.3:ide(Help1)+and(w0.3:ide(Help2))+and(w0.3:ide(Help3))+and(w0.3:ide(Help4))+and(w0.3:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)

d0.3_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.3_h1:vm(Help1,ainv)+and(w0.3_h2:vm(Help2,ainv))+and(w0.3_h3:vm(Help3,ainv))+and(w0.3_h4:vm(Help4,ainv))+and(w0.3_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.3_h1:ide(Help1)+and(w0.3_h2:ide(Help2))+and(w0.3_h3:ide(Help3))+and(w0.3_h4:ide(Help4))+and(w0.3_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.4 Juvenile Tarsus ####
w0.4 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.4

w0.4_h1 <- w0.4 * juv.t$Help1_rel
w0.4_h2 <- w0.4 * juv.t$Help2_rel
w0.4_h3 <- w0.4 * juv.t$Help3_rel
w0.4_h4 <- w0.4 * juv.t$Help4_rel
w0.4_h5 <- w0.4 * juv.t$Help5_rel


d0.4 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.4:vm(Help1,ainv)+and(w0.4:vm(Help2,ainv))+and(w0.4:vm(Help3,ainv))+and(w0.4:vm(Help4,ainv))+and(w0.4:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.4:ide(Help1)+and(w0.4:ide(Help2))+and(w0.4:ide(Help3))+and(w0.4:ide(Help4))+and(w0.4:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)

d0.4_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.4_h1:vm(Help1,ainv)+and(w0.4_h2:vm(Help2,ainv))+and(w0.4_h3:vm(Help3,ainv))+and(w0.4_h4:vm(Help4,ainv))+and(w0.4_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.4_h1:ide(Help1)+and(w0.4_h2:ide(Help2))+and(w0.4_h3:ide(Help3))+and(w0.4_h4:ide(Help4))+and(w0.4_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.5 Juvenile Tarsus ####
w0.5 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.5

w0.5_h1 <- w0.5 * juv.t$Help1_rel
w0.5_h2 <- w0.5 * juv.t$Help2_rel
w0.5_h3 <- w0.5 * juv.t$Help3_rel
w0.5_h4 <- w0.5 * juv.t$Help4_rel
w0.5_h5 <- w0.5 * juv.t$Help5_rel


d0.5 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.5:vm(Help1,ainv)+and(w0.5:vm(Help2,ainv))+and(w0.5:vm(Help3,ainv))+and(w0.5:vm(Help4,ainv))+and(w0.5:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.5:ide(Help1)+and(w0.5:ide(Help2))+and(w0.5:ide(Help3))+and(w0.5:ide(Help4))+and(w0.5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.t)

d0.5_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.5_h1:vm(Help1,ainv)+and(w0.5_h2:vm(Help2,ainv))+and(w0.5_h3:vm(Help3,ainv))+and(w0.5_h4:vm(Help4,ainv))+and(w0.5_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.5_h1:ide(Help1)+and(w0.5_h2:ide(Help2))+and(w0.5_h3:ide(Help3))+and(w0.5_h4:ide(Help4))+and(w0.5_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.6 Juvenile Tarsus ####
w0.6 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.6

w0.6_h1 <- w0.6 * juv.t$Help1_rel
w0.6_h2 <- w0.6 * juv.t$Help2_rel
w0.6_h3 <- w0.6 * juv.t$Help3_rel
w0.6_h4 <- w0.6 * juv.t$Help4_rel
w0.6_h5 <- w0.6 * juv.t$Help5_rel


d0.6 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.6:vm(Help1,ainv)+and(w0.6:vm(Help2,ainv))+and(w0.6:vm(Help3,ainv))+and(w0.6:vm(Help4,ainv))+and(w0.6:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.6:ide(Help1)+and(w0.6:ide(Help2))+and(w0.6:ide(Help3))+and(w0.6:ide(Help4))+and(w0.6:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.t)

d0.6_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.6_h1:vm(Help1,ainv)+and(w0.6_h2:vm(Help2,ainv))+and(w0.6_h3:vm(Help3,ainv))+and(w0.6_h4:vm(Help4,ainv))+and(w0.6_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.6_h1:ide(Help1)+and(w0.6_h2:ide(Help2))+and(w0.6_h3:ide(Help3))+and(w0.6_h4:ide(Help4))+and(w0.6_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.7 Juvenile Tarsus ####
w0.7 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.7

w0.7_h1 <- w0.7 * juv.t$Help1_rel
w0.7_h2 <- w0.7 * juv.t$Help2_rel
w0.7_h3 <- w0.7 * juv.t$Help3_rel
w0.7_h4 <- w0.7 * juv.t$Help4_rel
w0.7_h5 <- w0.7 * juv.t$Help5_rel


d0.7 <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
               random = ~str(~vm(USFWS, ainv)+w0.7:vm(Help1,ainv)+and(w0.7:vm(Help2,ainv))+and(w0.7:vm(Help3,ainv))+and(w0.7:vm(Help4,ainv))+and(w0.7:vm(Help5,ainv)),
                             ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                 w0.7:ide(Help1)+and(w0.7:ide(Help2))+and(w0.7:ide(Help3))+and(w0.7:ide(Help4))+and(w0.7:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
               equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
               residual=~idv(units), na.action = na.method(x="include", y="include"),
               maxit=1000,
               data=juv.t)

d0.7_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.7_h1:vm(Help1,ainv)+and(w0.7_h2:vm(Help2,ainv))+and(w0.7_h3:vm(Help3,ainv))+and(w0.7_h4:vm(Help4,ainv))+and(w0.7_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.7_h1:ide(Help1)+and(w0.7_h2:ide(Help2))+and(w0.7_h3:ide(Help3))+and(w0.7_h4:ide(Help4))+and(w0.7_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.8 Juvenile Tarsus ####
w0.8 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.8

w0.8_h1 <- w0.8 * juv.t$Help1_rel
w0.8_h2 <- w0.8 * juv.t$Help2_rel
w0.8_h3 <- w0.8 * juv.t$Help3_rel
w0.8_h4 <- w0.8 * juv.t$Help4_rel
w0.8_h5 <- w0.8 * juv.t$Help5_rel


d0.8 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.8:vm(Help1,ainv)+and(w0.8:vm(Help2,ainv))+and(w0.8:vm(Help3,ainv))+and(w0.8:vm(Help4,ainv))+and(w0.8:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.8:ide(Help1)+and(w0.8:ide(Help2))+and(w0.8:ide(Help3))+and(w0.8:ide(Help4))+and(w0.8:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)


d0.8_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.8_h1:vm(Help1,ainv)+and(w0.8_h2:vm(Help2,ainv))+and(w0.8_h3:vm(Help3,ainv))+and(w0.8_h4:vm(Help4,ainv))+and(w0.8_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.8_h1:ide(Help1)+and(w0.8_h2:ide(Help2))+and(w0.8_h3:ide(Help3))+and(w0.8_h4:ide(Help4))+and(w0.8_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 0.9 Juvenile Tarsus ####
w0.9 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^0.9

w0.9_h1 <- w0.9 * juv.t$Help1_rel
w0.9_h2 <- w0.9 * juv.t$Help2_rel
w0.9_h3 <- w0.9 * juv.t$Help3_rel
w0.9_h4 <- w0.9 * juv.t$Help4_rel
w0.9_h5 <- w0.9 * juv.t$Help5_rel


d0.9 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
              random = ~str(~vm(USFWS, ainv)+w0.9:vm(Help1,ainv)+and(w0.9:vm(Help2,ainv))+and(w0.9:vm(Help3,ainv))+and(w0.9:vm(Help4,ainv))+and(w0.9:vm(Help5,ainv)),
                            ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                w0.9:ide(Help1)+and(w0.9:ide(Help2))+and(w0.9:ide(Help3))+and(w0.9:ide(Help4))+and(w0.9:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
              equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
              residual=~idv(units), na.action = na.method(x="include", y="include"),
              maxit=1000,
              data=juv.t)

d0.9_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                  random = ~str(~vm(USFWS, ainv)+w0.9_h1:vm(Help1,ainv)+and(w0.9_h2:vm(Help2,ainv))+and(w0.9_h3:vm(Help3,ainv))+and(w0.9_h4:vm(Help4,ainv))+and(w0.9_h5:vm(Help5,ainv)),
                                ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                    w0.9_h1:ide(Help1)+and(w0.9_h2:ide(Help2))+and(w0.9_h3:ide(Help3))+and(w0.9_h4:ide(Help4))+and(w0.9_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                  equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                  residual=~idv(units), na.action = na.method(x="include", y="include"),
                  maxit=1000,
                  data=juv.t)

#### 1 Juvenile Tarsus ####
w1 <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^1

w1_h1 <- w1 * juv.t$Help1_rel
w1_h2 <- w1 * juv.t$Help2_rel
w1_h3 <- w1 * juv.t$Help3_rel
w1_h4 <- w1 * juv.t$Help4_rel
w1_h5 <- w1 * juv.t$Help5_rel



d1 <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
            random = ~str(~vm(USFWS, ainv)+w1:vm(Help1,ainv)+and(w1:vm(Help2,ainv))+and(w1:vm(Help3,ainv))+and(w1:vm(Help4,ainv))+and(w1:vm(Help5,ainv)),
                          ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
              w1:ide(Help1)+and(w1:ide(Help2))+and(w1:ide(Help3))+and(w1:ide(Help4))+and(w1:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
            equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
            residual=~idv(units), na.action = na.method(x="include", y="include"),
            maxit=1000,
            data=juv.t)

d1_rel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                random = ~str(~vm(USFWS, ainv)+w1_h1:vm(Help1,ainv)+and(w1_h2:vm(Help2,ainv))+and(w1_h3:vm(Help3,ainv))+and(w1_h4:vm(Help4,ainv))+and(w1_h5:vm(Help5,ainv)),
                              ~diag(2):vm(USFWS, ainv))+ vm(dad,ainv) +
                  w1_h1:ide(Help1)+and(w1_h2:ide(Help2))+and(w1_h3:ide(Help3))+and(w1_h4:ide(Help4))+and(w1_h5:ide(Help5)) + ide(dad) +NatalNest+NatalYear,
                equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                residual=~idv(units), na.action = na.method(x="include", y="include"),
                maxit=1000,
                data=juv.t)

#### Save Juvenile Tarsus Output ####

juv.t_df <- data.frame(Model = c("d0","d0.1","d0.2","d0.3","d0.4","d0.5","d0.6","d0.7","d0.8","d0.9","d1",
                                 "d0_rel","d0.1_rel","d0.2_rel","d0.3_rel","d0.4_rel","d0.5_rel","d0.6_rel","d0.7_rel","d0.8_rel","d0.9_rel","d1_rel"), 
                       LogL = c(d0$loglik,d0.1$loglik,d0.2$loglik,d0.3$loglik,d0.4$loglik,d0.5$loglik,d0.6$loglik,d0.7$loglik,d0.8$loglik,d0.9$loglik,d1$loglik,
                                d0_rel$loglik,d0.1_rel$loglik,d0.2_rel$loglik,d0.3_rel$loglik,d0.4_rel$loglik,d0.5_rel$loglik,d0.6_rel$loglik,d0.7_rel$loglik,d0.8_rel$loglik,d0.9_rel$loglik,d1_rel$loglik),
                       AIC = c(summary(d0)$aic,summary(d0.1)$aic,summary(d0.2)$aic,summary(d0.3)$aic,summary(d0.4)$aic,summary(d0.5)$aic,summary(d0.6)$aic,summary(d0.7)$aic,summary(d0.8)$aic,summary(d0.9)$aic,summary(d1)$aic,
                               summary(d0_rel)$aic,summary(d0.1_rel)$aic,summary(d0.2_rel)$aic,summary(d0.3_rel)$aic,summary(d0.4_rel)$aic,summary(d0.5_rel)$aic,summary(d0.6_rel)$aic,summary(d0.7_rel)$aic,summary(d0.8_rel)$aic,summary(d0.9_rel)$aic,summary(d1_rel)$aic),
                       BIC = c(summary(d0)$bic,summary(d0.1)$bic,summary(d0.2)$bic,summary(d0.3)$bic,summary(d0.4)$bic,summary(d0.5)$bic,summary(d0.6)$bic,summary(d0.7)$bic,summary(d0.8)$bic,summary(d0.9)$bic,summary(d1)$bic,
                               summary(d0_rel)$bic,summary(d0.1_rel)$bic,summary(d0.2_rel)$bic,summary(d0.3_rel)$bic,summary(d0.4_rel)$bic,summary(d0.5_rel)$bic,summary(d0.6_rel)$bic,summary(d0.7_rel)$bic,summary(d0.8_rel)$bic,summary(d0.9_rel)$bic,summary(d1_rel)$bic))

juv.t_models <- list(d0,d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9,d1,
                     d0_rel,d0.1_rel,d0.2_rel,d0.3_rel,d0.4_rel,d0.5_rel,d0.6_rel,d0.7_rel,d0.8_rel,d0.9_rel,d1_rel)


# Output of all models saved to rdata
save(d11.w_df, juv.w_df, 
     d11.t_df, juv.t_df,
     file = "DilutionParameter-DadHelper-Tables.rdata")

save(d11.w_models, juv.w_models,
     d11.t_models, d11.t_models,
     file = "DilutionParameter-DadHelper-Models.rdata")
