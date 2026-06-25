############# Main Model Script for Indirect Genetic Effects Across Ontogeny #####################

#### Script 4/5 in Florida Scrub Jay IGE Analysis
#### Original creator: Gladiana Spitz

#### Goal: This script runs IGE models with single and multiple social partners (moms, dads, and helpers) 
#### models on FSJ nestling/juvenile weight and tarsus in ASReml-R.

#### Output: data file containing model variance component table (Varcomps), 
#### heritability estimates for DGEs (herit) and total heritability (params)
#### social effect estimates (totaleffects)
#### fixed effect estimates (fixed_list)
#### All tables are in file ModelResults.rdata

# Load ASReml-R
library(asreml)

# Larger workspace makes this much easier
asreml.options(workspace="8gb") #if R crashes, increase this
#Filtered Data from 01-Cleaning.r
load("cleandata_FSJ_IGE_Analysis.RData")
#Invert pedigree as usual
ainv <- ainverse(pedPrunedASREML)

# General form of the model
# input object <- response variable~fixed effects
#                 random effects
#                 residual handling, missing value handling
#                 dataframe used


#### Nestling Weight ####

# Group size dilution parameter for use in helper models
HelpNum <- mean(d11.w$NumHelp)
d11.w$GroupSize <- d11.w$NumHelp + 3
d <- ((mean(d11.w$GroupSize) - 1)/(d11.w$GroupSize-1))^1

#### No IGE model

s_d11.wmodel <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~vm(USFWS,ainv)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.w)

s_d11.wvarcomp <- summary(s_d11.wmodel)$varcomp
s_d11.wvarcomp$trait <- "Nestling Weight"
s_d11.w_h2 <- vpredict(s_d11.wmodel, h2~V3/(V1+V2+V3+V4))

#### Single partner models

# Maternal Effect, Nestling Weight
m_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~vm(USFWS,ainv)+vm(mom,ainv)+ide(mom)+NatalNest+NatalYear,
                      residual=~idv(units), 
                      na.action = na.method(x="include", y="include"),		
                      maxit=1000,			
                      data=d11.w)
m_d11.wvarcomp <- summary(m_d11.wmodel)$varcomp 
# Pulling parameter
# narrow-sense heritability: Additive genetic variance/Total Variance
m_d11.w_h2 <- vpredict(m_d11.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
# total social effect: Indirect genetic variance + Indirect Environmental variance / Total Variance
m_d11.w_Tot <- vpredict(m_d11.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
# genetic social effect: Indirect genetic variance / Total Variance
m_d11.w_g <- vpredict(m_d11.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
# environmental social effect: Indirect Environmental variance / Total Variance
m_d11.w_e <- vpredict(m_d11.wmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
# Total heritable variation: Additive genetic variance + 2(num in group-1)*cov(V_A, V_IGE) + indirect genetic variance
m_d11.w_TBV <- vpredict(m_d11.wmodel, TBV~V3+V4)
# relative total heritable variation: TBV / total variance
m_d11.w_tau <- vpredict(m_d11.wmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))


# Maternal Effect with covariance structure
mc_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(vm(USFWS,ainv)+vm(mom,ainv), ~us(2):vm(USFWS, ainv))+ide(mom)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=d11.w)

mc_d11.wvarcomp <- summary(mc_d11.wmodel)$varcomp 
mc_d11.w_h2 <- vpredict(mc_d11.wmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
mc_d11.w_Tot <- vpredict(mc_d11.wmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
mc_d11.w_g <- vpredict(mc_d11.wmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
mc_d11.w_e <- vpredict(mc_d11.wmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
mc_d11.w_TBV <- vpredict(mc_d11.wmodel, TBV~V4+2*V5+V6)
mc_d11.w_tau <- vpredict(mc_d11.wmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

# Paternal Effect
p_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~vm(USFWS,ainv)+vm(dad,ainv)+ide(dad)+NatalNest+NatalYear,
                      residual=~idv(units), 
                      na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=d11.w)
p_d11.wvarcomp <- summary(p_d11.wmodel)$varcomp 
p_d11.w_h2 <- vpredict(p_d11.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
p_d11.w_Tot <- vpredict(p_d11.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
p_d11.w_g <- vpredict(p_d11.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
p_d11.w_e <- vpredict(p_d11.wmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
p_d11.w_TBV <- vpredict(p_d11.wmodel, TBV~V3+V4)
p_d11.w_tau <- vpredict(p_d11.wmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

# Paternal Effect with covariance structure
pc_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~vm(USFWS,ainv)+vm(dad,ainv)+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.w)
pc_d11.wvarcomp <- summary(pc_d11.wmodel)$varcomp 
pc_d11.w_h2 <- vpredict(pc_d11.wmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
pc_d11.w_Tot <- vpredict(pc_d11.wmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
pc_d11.w_g <- vpredict(pc_d11.wmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
pc_d11.w_e <- vpredict(pc_d11.wmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
pc_d11.w_TBV <- vpredict(pc_d11.wmodel, TBV~V4+2*V5+V6)
pc_d11.w_tau <- vpredict(pc_d11.wmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

# Helper Effect

h_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                        d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                      equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=d11.w)
h_d11.wvarcomp <- summary(h_d11.wmodel)$varcomp
h_d11.w_h2 <- vpredict(h_d11.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
h_d11.w_Tot <- vpredict(h_d11.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
h_d11.w_g <- vpredict(h_d11.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
h_d11.w_e <- vpredict(h_d11.wmodel, Se~V5/(V3+V5+V1+V2+V4+V6))
h_d11.w_TBV <- vpredict(h_d11.wmodel, TBV~V3+(HelpNum^2)*V4)
h_d11.w_tau <- vpredict(h_d11.wmodel, tau~(V3+(HelpNum^2)*V4)/(V1+V2+V3+V4+V5+V6))

# Helper Effect with covariance structure

hc_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.w)
hc_d11.wvarcomp <- summary(hc_d11.wmodel)$varcomp
hc_d11.w_h2 <- vpredict(hc_d11.wmodel, h2~V3/(V1+V2+V3+V5+V6+V7))
hc_d11.w_Tot <- vpredict(hc_d11.wmodel, TotS~(V4+V5)/(V1+V2+V3+V5+V6+V7))
hc_d11.w_g <- vpredict(hc_d11.wmodel, Sg~V4/(V1+V2+V3+V5+V6+V7))
hc_d11.w_e <- vpredict(hc_d11.wmodel, Se~V5/(V1+V2+V3+V5+V6+V7))
hc_d11.w_TBV <- vpredict(hc_d11.wmodel, TBV~V3+2*HelpNum*V4+(HelpNum^2)*V5)
hc_d11.w_tau <- vpredict(hc_d11.wmodel, tau~(V3+2*HelpNum*V4+(HelpNum^2)*V5)/(V1+V2+V3+V5+V6+V7))


#### Multiple partner models

# Mat + Pat
mp_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=d11.w)

mp_d11.wvarcomp <- summary(mp_d11.wmodel)$varcomp 
mp_d11.w_h2 <- vpredict(mp_d11.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.w_Tot <- vpredict(mp_d11.wmodel, TotS~(V4+V5+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.w_g <- vpredict(mp_d11.wmodel, Sg~(V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.w_e <- vpredict(mp_d11.wmodel, Se~(V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.w_TBV <- vpredict(mp_d11.wmodel, TBV~V3+V4+V5)
mp_d11.w_tau <- vpredict(mp_d11.wmodel, tau~(V3+V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + helpers

mh_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(mom,ainv)+ide(mom)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.w)

mh_d11.wvarcomp <- summary(mh_d11.wmodel)$varcomp
mh_d11.w_h2 <- vpredict(mh_d11.wmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.w_Tot <- vpredict(mh_d11.wmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.w_g <- vpredict(mh_d11.wmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.w_e <- vpredict(mh_d11.wmodel, Se~(V4+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.w_TBV <- vpredict(mh_d11.wmodel, TBV~V5+V3+(HelpNum^2)*V6)
mh_d11.w_tau <- vpredict(mh_d11.wmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Pat + helpers

ph_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random =  ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(dad,ainv)+
                         ide(dad)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.w)

ph_d11.wvarcomp <- summary(ph_d11.wmodel)$varcomp
ph_d11.w_h2 <- vpredict(ph_d11.wmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.w_Tot <- vpredict(ph_d11.wmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.w_g <- vpredict(ph_d11.wmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.w_e <- vpredict(ph_d11.wmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.w_TBV <- vpredict(ph_d11.wmodel, TBV~V5+V3+(HelpNum^2)*V6)
ph_d11.w_tau <- vpredict(ph_d11.wmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + Pat + helpers

mph_d11.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                        random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                          vm(mom,ainv)+vm(dad,ainv)+
                          ide(mom)+ide(dad)+
                          d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                          NatalYear+NatalNest,
                        equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                        residual=~idv(units), na.action = na.method(x="include", y="include"),
                        maxit=1000,
                        data=d11.w)

mph_d11.wvarcomp <- summary(mph_d11.wmodel)$varcomp
mph_d11.w_h2 <- vpredict(mph_d11.wmodel, h2~V7/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.w_Tot <- vpredict(mph_d11.wmodel, TotS~(V3+V4+V5+V6+V8+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.w_g <- vpredict(mph_d11.wmodel, Sg~(V3+V4+V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.w_e <- vpredict(mph_d11.wmodel, Se~(V5+V6+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.w_TBV <- vpredict(mph_d11.wmodel, TBV~V7+V3+V4+(HelpNum^2)*V8)
mph_d11.w_tau <- vpredict(mph_d11.wmodel, tau~(V7+V3+V4+(HelpNum^2)*V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))

#### Juvenile Weight ####

# Group size dilution parameter for use in helper models
HelpNum <- mean(juv.w$NumHelp)
juv.w$GroupSize <- juv.w$NumHelp + 3
d <- ((mean(juv.w$GroupSize) - 1)/(juv.w$GroupSize-1))^1

#### No IGE Model

s_juv.wmodel <- asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~vm(USFWS,ainv)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,		
                       data=juv.w)
s_juv.wvarcomp <- summary(s_juv.wmodel)$varcomp
s_juv.wvarcomp$trait <- "Juvenile Weight"
s_juv.w_h2 <- vpredict(s_juv.wmodel, h2~V3/(V1+V2+V3+V4))

#### Single social partners

## Maternal
m_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~vm(USFWS,ainv)+vm(mom,ainv)+ide(mom)+NatalNest+NatalYear,
                      residual=~idv(units), 
                      na.action = na.method(x="include", y="include"),
                      maxit=1000,		
                      data=juv.w)
m_juv.wvarcomp <- summary(m_juv.wmodel)$varcomp 
m_juv.w_h2 <- vpredict(m_juv.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
m_juv.w_Tot <- vpredict(m_juv.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
m_juv.w_g <- vpredict(m_juv.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
m_juv.w_e <- vpredict(m_juv.wmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
m_juv.w_TBV <- vpredict(m_juv.wmodel, TBV~V3+V4)
m_juv.w_tau <- vpredict(m_juv.wmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

## Maternal effect with covariance structure
mc_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(vm(USFWS,ainv)+vm(mom,ainv), ~us(2):vm(USFWS, ainv))+ide(mom)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.w)

mc_juv.wvarcomp <- summary(mc_juv.wmodel)$varcomp 
mc_juv.w_h2 <- vpredict(mc_juv.wmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
mc_juv.w_Tot <- vpredict(mc_juv.wmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
mc_juv.w_g <- vpredict(mc_juv.wmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
mc_juv.w_e <- vpredict(mc_juv.wmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
mc_juv.w_TBV <- vpredict(mc_juv.wmodel, TBV~V4+2*V5+V6)
mc_juv.w_tau <- vpredict(mc_juv.wmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

## Paternal
p_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~vm(USFWS,ainv)+vm(dad,ainv)+ide(dad)+NatalNest+NatalYear,
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=juv.w)
p_juv.wvarcomp <- summary(p_juv.wmodel)$varcomp 
p_juv.w_h2 <- vpredict(p_juv.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
p_juv.w_Tot <- vpredict(p_juv.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
p_juv.w_g <- vpredict(p_juv.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
p_juv.w_e <- vpredict(p_juv.wmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
p_juv.w_TBV <- vpredict(p_juv.wmodel, TBV~V3+V4)
p_juv.w_tau <- vpredict(p_juv.wmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

## Paternal effect with covariance structure
pc_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(vm(USFWS,ainv)+vm(dad,ainv), ~us(2):vm(USFWS, ainv))+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.w)

pc_juv.wvarcomp <- summary(pc_juv.wmodel)$varcomp 
pc_juv.w_h2 <- vpredict(pc_juv.wmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
pc_juv.w_Tot <- vpredict(pc_juv.wmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
pc_juv.w_g <- vpredict(pc_juv.wmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
pc_juv.w_e <- vpredict(pc_juv.wmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
pc_juv.w_TBV <- vpredict(pc_juv.wmodel, TBV~V4+2*V5+V6)
pc_juv.w_tau <- vpredict(pc_juv.wmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

## Helper
h_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                                    ~diag(2):vm(USFWS, ainv))+
                        d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                      equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=juv.w)
h_juv.wvarcomp <- summary(h_juv.wmodel)$varcomp
h_juv.w_h2 <- vpredict(h_juv.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
h_juv.w_Tot <- vpredict(h_juv.wmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
h_juv.w_g <- vpredict(h_juv.wmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
h_juv.w_e <- vpredict(h_juv.wmodel, Se~V5/(V3+V5+V1+V2+V4+V6))
h_juv.w_TBV <- vpredict(h_juv.wmodel, TBV~V3+(HelpNum^2)*V4)
h_juv.w_tau <- vpredict(h_juv.wmodel, tau~(V3+(HelpNum^2)*V4)/(V1+V2+V3+V4+V5+V6))

# Helper Effect with covariance structure
hc_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.w)
hc_juv.wvarcomp <- summary(hc_juv.wmodel)$varcomp
hc_juv.w_h2 <- vpredict(hc_juv.wmodel, h2~V3/(V1+V2+V3+V5+V6+V7))
hc_juv.w_Tot <- vpredict(hc_juv.wmodel, TotS~(V4+V5)/(V1+V2+V3+V5+V6+V7))
hc_juv.w_g <- vpredict(hc_juv.wmodel, Sg~V4/(V1+V2+V3+V5+V6+V7))
hc_juv.w_e <- vpredict(hc_juv.wmodel, Se~V5/(V1+V2+V3+V5+V6+V7))
hc_juv.w_TBV <- vpredict(hc_juv.wmodel, TBV~V3+2*HelpNum*V4+(HelpNum^2)*V5)
hc_d11.w_tau <- vpredict(hc_juv.wmodel, tau~(V3+2*HelpNum*V4+(HelpNum^2)*V5)/(V1+V2+V3+V5+V6+V7))

#### Multiple Social Partners
# Mat + Pat
mp_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.w)
mp_juv.wvarcomp <- summary(mp_juv.wmodel)$varcomp 
mp_juv.w_h2 <- vpredict(mp_juv.wmodel, h2~V3/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.w_Tot <- vpredict(mp_juv.wmodel, TotS~(V4+V5+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.w_g <- vpredict(mp_juv.wmodel, Sg~(V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.w_e <- vpredict(mp_juv.wmodel, Se~(V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.w_TBV <- vpredict(mp_juv.wmodel, TBV~V3+V4+V5)
mp_juv.w_tau <- vpredict(mp_juv.wmodel, tau~(V3+V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + helpers

mh_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(mom,ainv)+
                         ide(mom)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.w)

mh_juv.wvarcomp <- summary(mh_juv.wmodel)$varcomp
mh_juv.w_h2 <- vpredict(mh_juv.wmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.w_Tot <- vpredict(mh_juv.wmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.w_g <- vpredict(mh_juv.wmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.w_e <- vpredict(mh_juv.wmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.w_TBV <- vpredict(mh_juv.wmodel, TBV~V5+V3+(HelpNum^2)*V6)
mh_juv.w_tau <- vpredict(mh_juv.wmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Pat + helpers
ph_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(dad,ainv)+
                         ide(dad)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.w)

ph_juv.wvarcomp <- summary(ph_juv.wmodel)$varcomp
ph_juv.w_h2 <- vpredict(ph_juv.wmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.w_Tot <- vpredict(ph_juv.wmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.w_g <- vpredict(ph_juv.wmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.w_e <- vpredict(ph_juv.wmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.w_TBV <- vpredict(ph_juv.wmodel, TBV~V5+V3+(HelpNum^2)*V6)
ph_juv.w_tau <- vpredict(ph_juv.wmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + Pat + helpers

mph_juv.wmodel <-asreml(fixed=Weight~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                        random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                          vm(mom,ainv)+vm(dad,ainv)+
                          ide(mom)+ide(dad)+
                          d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                          NatalYear+NatalNest,
                        equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                        residual=~idv(units), na.action = na.method(x="include", y="include"),
                        maxit=1000,
                        data=juv.w)

mph_juv.wvarcomp <- summary(mph_juv.wmodel)$varcomp
mph_juv.w_h2 <- vpredict(mph_juv.wmodel, h2~V7/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.w_Tot <- vpredict(mph_juv.wmodel, TotS~(V3+V4+V5+V6+V8+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.w_g <- vpredict(mph_juv.wmodel, Sg~(V3+V4+V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.w_e <- vpredict(mph_juv.wmodel, Se~(V5+V6+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.w_TBV <- vpredict(mph_juv.wmodel, TBV~V7+V3+V4+(HelpNum^2)*V8)
mph_juv.w_tau <- vpredict(mph_juv.wmodel, tau~(V7+V3+V4+(HelpNum^2)*V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))

#### Nestling Tarsus ####
HelpNum <- mean(d11.t$NumHelp)
d11.t$GroupSize <- d11.t$NumHelp + 3
d <- ((mean(d11.t$GroupSize) - 1)/(d11.t$GroupSize-1))^0

#### No IGE Model
s_d11.tmodel <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~vm(USFWS,ainv)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.t)

s_d11.tvarcomp <- summary(s_d11.tmodel)$varcomp
s_d11.tvarcomp$trait <- "Nestling Tarsus"
s_d11.t_h2 <- vpredict(s_d11.tmodel, h2~V3/(V1+V2+V3+V4))

#### Single partner models

## Maternal model
m_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~vm(USFWS,ainv)+vm(mom,ainv)+ide(mom)+NatalNest+NatalYear,
                      residual=~idv(units), na.action = na.method(x="include", y="include"),		
                      maxit=1000,			
                      data=d11.t)
m_d11.tvarcomp <- summary(m_d11.tmodel)$varcomp 
m_d11.t_h2 <- vpredict(m_d11.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
m_d11.t_Tot <- vpredict(m_d11.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
m_d11.t_g <- vpredict(m_d11.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
m_d11.t_e <- vpredict(m_d11.tmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
m_d11.t_TBV <- vpredict(m_d11.tmodel, TBV~V3+V4)
m_d11.t_tau <- vpredict(m_d11.tmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

# Maternal Effect with covariance structure
mc_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(vm(USFWS,ainv)+vm(mom,ainv), ~us(2):vm(USFWS, ainv))+ide(mom)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=d11.t)

mc_d11.tvarcomp <- summary(mc_d11.tmodel)$varcomp 
mc_d11.t_h2 <- vpredict(mc_d11.tmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
mc_d11.t_Tot <- vpredict(mc_d11.tmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
mc_d11.t_g <- vpredict(mc_d11.tmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
mc_d11.t_e <- vpredict(mc_d11.tmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
mc_d11.t_TBV <- vpredict(mc_d11.tmodel, TBV~V4+2*V5+V6)
mc_d11.t_tau <- vpredict(mc_d11.tmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))


## Paternal model
p_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~vm(USFWS,ainv)+vm(dad,ainv)+ide(dad)+NatalNest+NatalYear,
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=d11.t)
p_d11.tvarcomp <- summary(p_d11.tmodel)$varcomp 
p_d11.t_h2 <- vpredict(p_d11.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
p_d11.t_Tot <- vpredict(p_d11.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
p_d11.t_g <- vpredict(p_d11.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
p_d11.t_e <- vpredict(p_d11.tmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
p_d11.t_TBV <- vpredict(p_d11.tmodel, TBV~V3+V4)
p_d11.t_tau <- vpredict(p_d11.tmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

# Maternal Effect with covariance structure
pc_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(vm(USFWS,ainv)+vm(dad,ainv), ~us(2):vm(USFWS, ainv))+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=d11.t)

pc_d11.tvarcomp <- summary(pc_d11.tmodel)$varcomp 
pc_d11.t_h2 <- vpredict(pc_d11.tmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
pc_d11.t_Tot <- vpredict(pc_d11.tmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
pc_d11.t_g <- vpredict(pc_d11.tmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
pc_d11.t_e <- vpredict(pc_d11.tmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
pc_d11.t_TBV <- vpredict(pc_d11.tmodel, TBV~V4+2*V5+V6)
pc_d11.t_tau <- vpredict(pc_d11.tmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

## Helper models
h_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                      random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)),
                                    ~diag(2):vm(USFWS, ainv))+
                        d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                      equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=d11.t)

h_d11.tvarcomp <- summary(h_d11.tmodel)$varcomp
h_d11.t_h2 <- vpredict(h_d11.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
h_d11.t_Tot <- vpredict(h_d11.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
h_d11.t_g <- vpredict(h_d11.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
h_d11.t_e <- vpredict(h_d11.tmodel, Se~V5/(V3+V5+V1+V2+V4+V6))
h_d11.t_TBV <- vpredict(h_d11.tmodel, TBV~V3+(HelpNum^2)*V4)
h_d11.t_tau <- vpredict(h_d11.tmodel, tau~(V3+(HelpNum^2)*V4)/(V1+V2+V3+V4+V5+V6))

# Helper Effect with covariance structure

hc_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.t)
hc_d11.vtarcomp <- summary(hc_d11.tmodel)$varcomp
hc_d11.t_h2 <- vpredict(hc_d11.tmodel, h2~V3/(V1+V2+V3+V5+V6+V7))
hc_d11.t_Tot <- vpredict(hc_d11.tmodel, TotS~(V4+V5)/(V1+V2+V3+V5+V6+V7))
hc_d11.t_g <- vpredict(hc_d11.tmodel, Sg~V4/(V1+V2+V3+V5+V6+V7))
hc_d11.t_e <- vpredict(hc_d11.tmodel, Se~V5/(V1+V2+V3+V5+V6+V7))
hc_d11.t_TBV <- vpredict(hc_d11.tmodel, TBV~V3+2*HelpNum*V4+(HelpNum^2)*V5)
hc_d11.t_tau <- vpredict(hc_d11.tmodel, tau~(V3+2*HelpNum*V4+(HelpNum^2)*V5)/(V1+V2+V3+V5+V6+V7))

#### Multiple partner models

# Mat + Pat
mp_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=d11.t)

mp_d11.tvarcomp <- summary(mp_d11.tmodel)$varcomp 
mp_d11.t_h2 <- vpredict(mp_d11.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.t_Tot <- vpredict(mp_d11.tmodel, TotS~(V4+V5+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.t_g <- vpredict(mp_d11.tmodel, Sg~(V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.t_e <- vpredict(mp_d11.tmodel, Se~(V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_d11.t_TBV <- vpredict(mp_d11.tmodel, TBV~V3+V4+V5)
mp_d11.t_tau <- vpredict(mp_d11.tmodel, tau~(V3+V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + helpers
mh_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(mom,ainv)+
                         ide(mom)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.t)
mh_d11.tvarcomp <- summary(mh_d11.tmodel)$varcomp
mh_d11.t_h2 <- vpredict(mh_d11.tmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.t_Tot <- vpredict(mh_d11.tmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.t_g <- vpredict(mh_d11.tmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.t_e <- vpredict(mh_d11.tmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_d11.t_TBV <- vpredict(mh_d11.tmodel, TBV~V5+V3+(HelpNum^2)*V6)
mh_d11.t_tau <- vpredict(mh_d11.tmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Pat + helpers

ph_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                       random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(dad,ainv)+
                         ide(dad)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=d11.t)
ph_d11.tvarcomp <- summary(ph_d11.tmodel)$varcomp
ph_d11.t_h2 <- vpredict(ph_d11.tmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.t_Tot <- vpredict(ph_d11.tmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.t_g <- vpredict(ph_d11.tmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.t_e <- vpredict(ph_d11.tmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_d11.t_TBV <- vpredict(ph_d11.tmodel, TBV~V5+V3+(HelpNum^2)*V6)
ph_d11.t_tau <- vpredict(ph_d11.tmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + Pat + helpers

mph_d11.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS,	
                        random = ~str(~vm(USFWS,ainv) + d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                          vm(mom,ainv)+vm(dad,ainv)+
                          ide(mom)+ide(dad)+
                          d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                          NatalYear+NatalNest,
                        equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                        residual=~idv(units), na.action = na.method(x="include", y="include"),
                        maxit=1000,
                        data=d11.t)

mph_d11.tvarcomp <- summary(mph_d11.tmodel)$varcomp
mph_d11.t_h2 <- vpredict(mph_d11.tmodel, h2~V7/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.t_Tot <- vpredict(mph_d11.tmodel, TotS~(V3+V4+V5+V6+V8+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.t_g <- vpredict(mph_d11.tmodel, Sg~(V3+V4+V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.t_e <- vpredict(mph_d11.tmodel, Se~(V5+V6+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_d11.t_TBV <- vpredict(mph_d11.tmodel, TBV~V7+V3+V4+(HelpNum^2)*V8)
mph_d11.t_tau <- vpredict(mph_d11.tmodel, tau~(V7+V3+V4+(HelpNum^2)*V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))

#### Juvenile Tarsus ####
HelpNum <- mean(juv.t$NumHelp)
juv.t$GroupSize <- juv.t$NumHelp + 3
d <- ((mean(juv.t$GroupSize) - 1)/(juv.t$GroupSize-1))^1

#### No IGE Model
s_juv.tmodel <- asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~vm(USFWS,ainv)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,		
                       data=juv.t)

s_juv.tvarcomp <- summary(s_juv.tmodel)$varcomp
s_juv.tvarcomp$trait <- "Juvenile Tarsus"
s_juv.t_h2 <- vpredict(s_juv.tmodel, h2~V3/(V1+V2+V3+V4))

#### Single partner models

## Maternal model
m_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~vm(USFWS,ainv)+vm(mom,ainv)+ide(mom)+NatalNest+NatalYear,
                      residual=~idv(units), 
                      na.action = na.method(x="include", y="include"),
                      maxit=1000,		
                      data=juv.t)

m_juv.tvarcomp <- summary(m_juv.tmodel)$varcomp 
m_juv.t_h2 <- vpredict(m_juv.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
m_juv.t_Tot <- vpredict(m_juv.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
m_juv.t_g <- vpredict(m_juv.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
m_juv.t_e <- vpredict(m_juv.tmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
m_juv.t_TBV <- vpredict(m_juv.tmodel, TBV~V3+V4)
m_juv.t_tau <- vpredict(m_juv.tmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

## Maternal effect with covariance structure
mc_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(vm(USFWS,ainv)+vm(mom,ainv), ~us(2):vm(USFWS, ainv))+ide(mom)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.t)

mc_juv.tvarcomp <- summary(mc_juv.tmodel)$varcomp 
mc_juv.t_h2 <- vpredict(mc_juv.tmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
mc_juv.t_Tot <- vpredict(mc_juv.tmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
mc_juv.t_g <- vpredict(mc_juv.tmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
mc_juv.t_e <- vpredict(mc_juv.tmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
mc_juv.t_TBV <- vpredict(mc_juv.tmodel, TBV~V4+2*V5+V6)
mc_juv.t_tau <- vpredict(mc_juv.tmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

## Paternal model
p_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~vm(USFWS,ainv)+vm(dad,ainv)+ide(dad)+NatalNest+NatalYear,
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=juv.t)

p_juv.tvarcomp <- summary(p_juv.tmodel)$varcomp 
p_juv.t_h2 <- vpredict(p_juv.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
p_juv.t_Tot <- vpredict(p_juv.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
p_juv.t_g <- vpredict(p_juv.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
p_juv.t_e <- vpredict(p_juv.tmodel, Se~V5/(V1+V2+V3+V4+V5+V6))
p_juv.t_TBV <- vpredict(p_juv.tmodel, TBV~V3+V4)
p_juv.t_tau <- vpredict(p_juv.tmodel, tau~(V3+V4)/(V1+V2+V3+V4+V5+V6))

## Paternal effect with covariance structure
pc_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(vm(USFWS,ainv)+vm(dad,ainv), ~us(2):vm(USFWS, ainv))+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), 
                       na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.t)

pc_juv.tvarcomp <- summary(pc_juv.tmodel)$varcomp 
pc_juv.t_h2 <- vpredict(pc_juv.tmodel, h2~V4/(V1+V2+V3+V4+V6+V7))
pc_juv.t_Tot <- vpredict(pc_juv.tmodel, TotS~(V3+V6)/(V1+V2+V3+V4+V6+V7))
pc_juv.t_g <- vpredict(pc_juv.tmodel, Sg~V6/(V1+V2+V3+V4+V6+V7))
pc_juv.t_e <- vpredict(pc_juv.tmodel, Se~V3/(V1+V2+V3+V4+V6+V7))
pc_juv.t_TBV <- vpredict(pc_juv.tmodel, TBV~V4+2*V5+V6)
pc_juv.t_tau <- vpredict(pc_juv.tmodel, tau~(V4+2*V5+V6)/(V1+V2+V3+V4+V6+V7))

## Helper model
h_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                      random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                        d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                      equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                      residual=~idv(units), na.action = na.method(x="include", y="include"),
                      maxit=1000,
                      data=juv.t)
h_juv.tvarcomp <- summary(h_juv.tmodel)$varcomp
h_juv.t_h2 <- vpredict(h_juv.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6))
h_juv.t_Tot <- vpredict(h_juv.tmodel, TotS~(V4+V5)/(V1+V2+V3+V4+V5+V6))
h_juv.t_g <- vpredict(h_juv.tmodel, Sg~V4/(V1+V2+V3+V4+V5+V6))
h_juv.t_e <- vpredict(h_juv.tmodel, Se~V5/(V3+V5+V1+V2+V4+V6))
h_juv.t_TBV <- vpredict(h_juv.tmodel, TBV~V3+(HelpNum^2)*V4)
h_juv.t_tau <- vpredict(h_juv.tmodel, tau~(V3+(HelpNum^2)*V4)/(V1+V2+V3+V4+V5+V6))

# Helper Effect with covariance structure
hc_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~us(2):vm(USFWS, ainv))+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+NatalNest+NatalYear,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.t)
hc_juv.tvarcomp <- summary(hc_juv.tmodel)$varcomp
hc_juv.t_h2 <- vpredict(hc_juv.tmodel, h2~V3/(V1+V2+V3+V5+V6+V7))
hc_juv.t_Tot <- vpredict(hc_juv.tmodel, TotS~(V4+V5)/(V1+V2+V3+V5+V6+V7))
hc_juv.t_g <- vpredict(hc_juv.tmodel, Sg~V4/(V1+V2+V3+V5+V6+V7))
hc_juv.t_e <- vpredict(hc_juv.tmodel, Se~V5/(V1+V2+V3+V5+V6+V7))
hc_juv.t_TBV <- vpredict(hc_juv.tmodel, TBV~V3+2*HelpNum*V4+(HelpNum^2)*V5)
hc_d11.t_tau <- vpredict(hc_juv.tmodel, tau~(V3+2*HelpNum*V4+(HelpNum^2)*V5)/(V1+V2+V3+V5+V6+V7))


#### Multiple partner models

# Mat + Pat
mp_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~vm(USFWS,ainv)+vm(mom,ainv)+vm(dad,ainv)+ide(mom)+ide(dad)+NatalNest+NatalYear,
                       residual=~idv(units), na.action = na.method(x="include", y="include"),		
                       maxit=1000,			
                       data=juv.t)

mp_juv.tvarcomp <- summary(mp_juv.tmodel)$varcomp 
mp_juv.t_h2 <- vpredict(mp_juv.tmodel, h2~V3/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.t_Tot <- vpredict(mp_juv.tmodel, TotS~(V4+V5+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.t_g <- vpredict(mp_juv.tmodel, Sg~(V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.t_e <- vpredict(mp_juv.tmodel, Se~(V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mp_juv.t_TBV <- vpredict(mp_juv.tmodel, TBV~V3+V4+V5)
mp_juv.t_tau <- vpredict(mp_juv.tmodel, tau~(V3+V4+V5)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + helpers

mh_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(mom,ainv)+
                         ide(mom)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.t)

mh_juv.tvarcomp <- summary(mh_juv.tmodel)$varcomp
mh_juv.t_h2 <- vpredict(mh_juv.tmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.t_Tot <- vpredict(mh_juv.tmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.t_g <- vpredict(mh_juv.tmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.t_e <- vpredict(mh_juv.tmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
mh_juv.t_TBV <- vpredict(mh_juv.tmodel, TBV~V5+V3+(HelpNum^2)*V6)
mh_juv.t_tau <- vpredict(mh_juv.tmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Pat + helpers

ph_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                       random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                         vm(dad,ainv)+
                         ide(dad)+
                         d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                         NatalYear+NatalNest,
                       equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                       residual=~idv(units), na.action = na.method(x="include", y="include"),
                       maxit=1000,
                       data=juv.t)
ph_juv.tvarcomp <- summary(ph_juv.tmodel)$varcomp
ph_juv.t_h2 <- vpredict(ph_juv.tmodel, h2~V5/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.t_Tot <- vpredict(ph_juv.tmodel, TotS~(V3+V4+V6+V7)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.t_g <- vpredict(ph_juv.tmodel, Sg~(V3+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.t_e <- vpredict(ph_juv.tmodel, Se~(V4+V6)/(V1+V2+V3+V4+V5+V6+V7+V8))
ph_juv.t_TBV <- vpredict(ph_juv.tmodel, TBV~V5+V3+(HelpNum^2)*V6)
ph_juv.t_tau <- vpredict(ph_juv.tmodel, tau~(V5+V3+(HelpNum^2)*V6)/(V1+V2+V3+V4+V5+V6+V7+V8))

# Mat + Pat + helpers

mph_juv.tmodel <-asreml(fixed=Tarsus~1+Sex+HatchDOYS+HatchNumS+NumHelpS+TerrSizeS+ICS+AgeMeasS,	
                        random = ~str(~vm(USFWS, ainv)+d:vm(Help1,ainv)+and(d:vm(Help2,ainv))+and(d:vm(Help3,ainv))+and(d:vm(Help4,ainv))+and(d:vm(Help5,ainv)), ~diag(2):vm(USFWS, ainv))+
                          vm(mom,ainv)+vm(dad,ainv)+
                          ide(mom)+ide(dad)+
                          d:ide(Help1)+and(d:ide(Help2))+and(d:ide(Help3))+and(d:ide(Help4))+and(d:ide(Help5))+
                          NatalYear+NatalNest,
                        equate.levels = c("Help1","Help2","Help3","Help4","Help5"),
                        residual=~idv(units), na.action = na.method(x="include", y="include"),
                        maxit=1000,
                        data=juv.t)
mph_juv.tvarcomp <- summary(mph_juv.tmodel)$varcomp
mph_juv.t_h2 <- vpredict(mph_juv.tmodel, h2~V7/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.t_Tot <- vpredict(mph_juv.tmodel, TotS~(V3+V4+V5+V6+V8+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.t_g <- vpredict(mph_juv.tmodel, Sg~(V3+V4+V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.t_e <- vpredict(mph_juv.tmodel, Se~(V5+V6+V9)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))
mph_juv.t_TBV <- vpredict(mph_juv.tmodel, TBV~V7+V3+V4+(HelpNum^2)*V8)
mph_juv.t_tau <- vpredict(mph_juv.tmodel, tau~(V7+V3+V4+(HelpNum^2)*V8)/(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10))

## You can export all data at this point. However, I prefer to pull the pieces I need and organize them

library(tidyverse)

# Vector of trait names
traits <- c("d11.w","juv.w","d11.t","juv.t")
#vector of model names 
socialpartner <- c("DGE-Only",
                   "Maternal","Maternal Cov",
                   "Paternal","Paternal Cov",
                   "Helper","Helper Cov",
                   "MatPat", "MatHelper", "PatHelper", 
                   "MatPatHelper")
names(socialpartner) <- c("s","m","mc","p","pc","h","hc","mp","mh","ph","mph")
#initialize data frames used in 04-Plotting.R
herit <- data.frame() #heritabilities
Varcomps <- data.frame() #Variance components
totaleffects <- data.frame() #Social effects
params <- data.frame() #TBV and tau measures
fixed_list <- list() # fixed effects table

# This loop parses information from each parameter to add to a single data frame for plotting
for(j in 1:length(socialpartner)){
  for(i in 1:length(traits)){
    # get model name
    mod = paste(names(socialpartner)[j], traits[i], sep = "_")
    
    # get variance component table
    varcomp = eval(parse(text=paste(mod,"varcomp",sep="")))
    # add to variance component table
    varcomp$Trait = traits[i]
    varcomp$Model = socialpartner[j]
    Varcomps = rbind(Varcomps,varcomp)	
    
    # get heritability
    h2 = eval(parse(text=paste(mod,"h2",sep="_")))
    # add to heritability table
    h2$Trait = traits[i]
    h2$Model = socialpartner[j]
    herit = rbind(herit,h2)
    if(names(socialpartner)[j] != "s"){
      # get total social effect
      Tot = eval(parse(text=paste(mod,"Tot",sep="_")))
      Tot$Trait = traits[i]
      Tot$Model = socialpartner[j]
      Tot$type = "Total"
      
      # get social genetic effect
      g = eval(parse(text=paste(mod,"g",sep="_")))
      g$Trait = traits[i]
      g$Model = socialpartner[j]
      g$type = "Genetic"
      
      # get social environmental effect
      e = eval(parse(text=paste(mod,"e",sep="_")))
      e$Trait = traits[i]
      e$Model = socialpartner[j]
      e$type = "environmental"
      
      # make social effects table
      totaleffects = rbind(totaleffects,Tot,g,e)
      
      # get total heritable variance
      TBV = eval(parse(text=paste(mod,"TBV",sep="_")))
      TBV$Trait = traits[i]
      TBV$Model = socialpartner[j]
      TBV$type = "TBV"
      
      # get total heritability
      tau = eval(parse(text=paste(mod,"tau",sep="_")))
      tau$Trait = traits[i]
      tau$Model =socialpartner[j]
      tau$type = "tau"
      
      # make social parameter table
      params = rbind(params,TBV,tau)
    }
    fixed_list[[paste(names(socialpartner)[j], traits[i], sep = "_")]] <- summary(mod, coef=T)$coef.fixed
  }
}

# Output of all models saved to rdata
save(herit, Varcomps, totaleffects, params, fixed_list, file = "ModelResults.rdata")
