# Analysis for Indirect genetic effect across ontengy in Florida scrub-jays

This repository stores raw data and scripts used in the submitted manuscript for Spitz et al. 2026 studying the effect of social partners on juvenile morphometrics in Florida scrub-jays.

The input data is found in cleandata_FSJ_IGE_Analysis.RData which contains the following dataframes:

d11.w - Juvenile individual ID (USFWS), nest of origin of the focal juvenile (NatalNest), Sex, the weight of the juvenile at day 11 (Weight), standardized hatch date (HatchDOYS), standardized clutch size (HatchNumS), number of helpers at the nest (NumHelp) that is standardized (NumHelpS), standardized territory size (TerrSizeS), standardized inbreeding coefficient (ICS), identity of the father (dad), identity of the mother (mom), identity of up to 5 helpers (Help1-Help5), relatedness of up to 5 helpers (Help1_rel-Help5_rel)

juv.w - Juvenile individual ID (USFWS), nest of origin of the focal juvenile (NatalNest), Sex, the weight of the juvenile between day 60 and day 90 (Weight), standardized hatch date (HatchDOYS), standardized clutch size (HatchNumS), number of helpers at the nest (NumHelp) that is standardized (NumHelpS), standardized territory size (TerrSizeS), standardized inbreeding coefficient (ICS), standardized age of weight measurement (AgeMeasS), identity of the father (dad), identity of the mother (mom), identity of up to 5 helpers (Help1-Help5), relatedness of up to 5 helpers (Help1_rel-Help5_rel)

d11.t - Juvenile individual ID (USFWS), nest of origin of the focal juvenile (NatalNest), Sex, tarsus length of the juvenile at day 11 (Tarsus), standardized hatch date (HatchDOYS), standardized clutch size (HatchNumS), number of helpers at the nest (NumHelp) that is standardized (NumHelpS), standardized territory size (TerrSizeS), standardized inbreeding coefficient (ICS), identity of the father (dad), identity of the mother (mom), identity of up to 5 helpers (Help1-Help5), relatedness of up to 5 helpers (Help1_rel-Help5_rel)

juv.t - Juvenile individual ID (USFWS), nest of origin of the focal juvenile (NatalNest), Sex, tarsus length of the juvenile between day 60 and day 90 (Tarsus), standardized hatch date (HatchDOYS), standardized clutch size (HatchNumS), number of helpers at the nest (NumHelp) that is standardized (NumHelpS), standardized territory size (TerrSizeS), standardized inbreeding coefficient (ICS), standardized age of tarsus measurement (AgeMeasS), identity of the father (dad), identity of the mother (mom), identity of up to 5 helpers (Help1-Help5), relatedness of up to 5 helpers (Help1_rel-Help5_rel)

pedPrunedASREML - Juvenile individual ID (USFWS), paternal ID (Dad), maternal ID (Mom)

02a-DilutionParameter-Helpers.R - Dilution parameter assessment of models with only helpers

02b-DilutionParameter-MomHelpers.R - Dilution parameter assessment of models with helpers and moms

02c-DilutionParameter-DadHelpers.R - Dilution parameter assessment of models with helpers and dads

02d-DilutionParameter-MomDadHelpers.R - Dilution parameter assessment of models with helpers, dads, and moms

03-SignificanceTesting.R - variance component testing by LRTs and fixed effect testing by walds

04-Modeling.R - primary modeling and calculating variance parameters
