#==========================#
# 0 - Cleaning             #
# Marquisseau Anaïs        #
#==========================#

# 2022
# Internship INRAE M2 BI
# Relationship between feed efficiency and metabolomics in sheep

#### Goals ####

# Cleaning of sheep zootechnical and metabolomic datasets from 2018 to 2020 (Bourges's experimental unit)

#### 0. Working directory and Session information ####

#Clear global environment
rm(list = ls(all.names = TRUE))

#free up memory and report the memory usage
gc()

#### 1. Load libraries and data ####

##### a) Libraries #####

library(tidyverse)
library(zCompositions) #Compositional data 
library(mixOmics)

#theme for ggplot
theme_set(theme_minimal())

##### b) Data #####

#File : variable's description in variable_description.txt file
DAC = readxl::read_xlsx("Data/zootechnie\ DAC.xlsx")

# Rearrange variable's types
DAC = DAC %>%
  
  #factor
  mutate(across(annee:lignee_gen, as_factor)) %>%
  mutate(across(COMNAI:mnal, as_factor)) %>%
  mutate(across(num_passage:NULOBA_annee, as_factor)) %>%
  
  #numeric
  mutate(across(age_sevrage:age.prelevements, as.numeric)) %>%
  mutate(across(MUMC:PAT, as.numeric)) %>%
  mutate_at(c("PODEC1","PODEC2", "POMIC1", "POFIC1", "POFIC2", "index", "cd", "RFI","DAPDC1","DAPDC2", "DAFIC1", "DAFIC2","date_prelevement","date_sevr", "DANAOV"), as.numeric) 

# Reorder factor's levels

DAC = DAC %>%
  mutate(annee = fct_relevel(annee, "2018", "2019", "2020")) 

DAC = DAC %>%
  mutate(generation = fct_relevel(generation, "G2", "G3")) 

DAC = DAC %>%
  mutate(lignee_gen = fct_relevel(lignee_gen, "rfi-_G2", "rfi+_G2", "rfi-_G3", "rfi+_G3")) 

DAC = DAC %>%
  mutate(NULOBA = fct_relevel(NULOBA, "1", "2", "3", "4", "5", "6")) 

DAC = DAC %>%
  mutate(NULOBA = fct_relevel(NULOBA, "1", "2", "3", "4", "5", "6")) 

DAC = DAC %>%
  mutate(NULOBA_annee = fct_relevel(NULOBA_annee, "1_2018", "2_2018","3_2018","4_2018","5_2018","6_2018","1_2019","2_2019","3_2019","4_2019","5_2019","1_2020","2_2020","3_2020","4_2020")) 

DAC = DAC %>%
  mutate(COMNAI = fct_relevel(COMNAI, "1", "2", "3", "4")) 

DAC = DAC %>%
  mutate(COMOEL = fct_relevel(COMOEL, "0","1", "2")) 

DAC = DAC %>%
  mutate(MEALLA = fct_relevel(MEALLA, "0","1")) 

DAC = DAC %>%
  mutate(mnal = fct_relevel(mnal, "1", "2", "3", "4", "6")) 

#### 2. Data checking ####

##### a. Missing values ##### 
na = apply(is.na(DAC), 2, which) 

# 1 not usable for the following steps : filter
`%notin%`=negate(`%in%`)
DAC = DAC %>% filter(row_number() %notin% na$RFI)

#### 3. New variables and outliers ####

###### a. FCR ######

DAC = DAC %>%
  add_column(FCR = DAC$CONMOY/DAC$GMQ)

##### b. Outliers #####

#Considered as outliers : 3sd to the mean per default

zscore = function(x, delimiter = 3, na.rm = TRUE) {
  centered = x - mean(x, na.rm = na.rm)
  limit = delimiter * sd(x, na.rm = na.rm)
  io = ifelse(abs(centered) > limit, 0, 1)
  return(io)
}

quanti = DAC %>% 
  dplyr::select_if(is.numeric)%>%
  dplyr::select(-date_sevr, -date_prelevement, -DAPDC1, -DAPDC2, -DAFIC1, -DAFIC2, -DANAOV, -cd)

#Filter the outliers
DAC = DAC %>%
  filter(animal %notin% c("188530", "188714", "104549", "188741", "104364", "104376", "194452", "188737"))

##### c. new RFI #####

# first, add a metabolic weight column :
DAC = DAC %>%
  add_column(POMET = DAC$POFIC**0.75)

#new rfi in corrected for year and pen 
lmCONMOY_2 = lm(CONMOY~GMQ+GRFC+MUFC+POMET+annee+NULOBA_annee %in% annee, data=DAC)

DAC = DAC %>%
  add_column(newRFI = residuals(lmCONMOY_2))

##### d. factRFI #####

#Create a factor version of the RFI
DAC = DAC %>%
  mutate(factRFI = as_factor(if_else(newRFI>0, "rfi+", "rfi-")))

#### 4. Corresponding buckets ####

corresp_RMN_nuofov = readxl::read_xlsx("Data/corresp_RMN_nuofov.xlsx")
corresp_RMN_nuofov$nuofov=gsub("^20000","",as.character(corresp_RMN_nuofov$nuofov))

#left join to add nmr id's only for animals with zootechnical values
DAC_nmr = left_join(
  x= DAC, 
  y= corresp_RMN_nuofov %>% 
    filter(regime == "DAC") %>%
    dplyr::select(-regime, -fluide,-annee),
  by = c("animal" = "nuofov")
)

#buckets
buckets_DAC_v2 = read_table("Data/RMN/buckets_DAC_v3_row")

#filter animals that don't have metabolomic data
DAC_nmr = DAC_nmr %>%
  filter(nom_anim_RMN %in% names(buckets_DAC_v2))

#filter animals that don't have zootechnical values
buckets_DAC_filtered = buckets_DAC_v2 %>%
  dplyr::select("ppm",DAC_nmr[["nom_anim_RMN"]])

#### 5. Corresponding metabolites ####

#metabolites quantifications
metabolites = readxl::read_xlsx("Data/quantifications/PlasmaDACv4.xlsx")
metabolites = metabolites %>%
  dplyr::select(Id, Citrate,`L-Leucine`, `beta-Hydroxyisovalerate`)
metabolites$Id = gsub("^20000","",as.character(metabolites$Id))

summary(metabolites)

DAC_nmr = left_join(DAC_nmr,
                    metabolites,
                    by = c("animal"="Id")
                    ) 

#### 6. Zero handling and transformations (buckets) ####

##### a. Water zone filtering #####
buckets_DAC_filtered = buckets_DAC_filtered %>%
  filter (ppm < 4.500 | ppm > 5.100)

##### b. 100% zeros buckets filtering #####

#filter the 100% zeros buckets
buckets_DAC_filtered %>%
  filter(rowSums(buckets_DAC_filtered[,2:ncol(buckets_DAC_filtered)])<= 0)%>%
  count()

buckets_DAC_filtered = buckets_DAC_filtered %>%
  filter(rowSums(buckets_DAC_filtered[,2:ncol(buckets_DAC_filtered)])> 0)

##### c. Zero imputation #####
#hp : zeros are values below the detection threshold (min non nul : 1.49e-09, 1.795ppm : nmr = low Sensitivity)

#Transpose the buckets matrix 
tbuckets_DAC = t(buckets_DAC_filtered) %>%
  as.data.frame()%>%
  set_names(as.character(slice(., 1))) %>%
  rownames_to_column("animal") %>%
  filter(animal != "ppm")

#Sort by animal nmr name
DAC_nmr = DAC_nmr %>%
  arrange(nom_anim_RMN)

tbuckets_DAC = tbuckets_DAC %>%
  arrange(animal)

#2.29% zeros in df
zPatterns(tbuckets_DAC[,2:ncol(tbuckets_DAC)], label =0)

#detection limit for each buckets : min for each bucket
dl = map_dfc(tbuckets_DAC[,2:878], function (x) {
  if_else(min(x)==0, min(replace(x, x==0,max(x))), min(x))
})

#imputed_buckets = multRepl(tbuckets_DAC[,2:ncol(tbuckets_DAC)], label = 0, frac = 0.65, dl = rep(1.4e-09,ncol(tbuckets_DAC)-1)) #dl = minimum of the dataframe
imputed_buckets = multRepl(tbuckets_DAC[,2:ncol(tbuckets_DAC)], label = 0, frac = 0.65, dl = dl)

#still zeros ? No Label 0 was not found in the data set
#zPatterns(imputed_buckets, label =0)

##### d. Transformation CLR #####

#geometric mean : applied on a sample x 
fct_gm = function(x) {exp(mean(log(x)))} 

#clr transformation : applied on the dataframe
fct_clr = function (x) {log(x / apply(x,1,fct_gm))}

clr_buckets = fct_clr(imputed_buckets)

#### 7. Correct phenotypes and buckets for fixed effects ####

quanti = DAC_nmr %>% 
  dplyr::select_if(is.numeric)%>%
  dplyr::select(-date_sevr, -date_prelevement, -DAPDC1, -DAPDC2, -DAFIC1, -DAFIC2, -DANAOV) %>%
  dplyr::select(-cd, -RFI, -PAT, -age_sevrage, -p_sevr, -PODEC1, -PODEC2, -POFIC1, -POFIC2, -POMIC1, -MUMC, -GRMC)

##### b. Correct #####
#linear model
fct_correct_fixedEffects_lm = function(x) {
  lm = lm(x~DAC_nmr$annee + DAC_nmr$NULOBA_annee%in%DAC_nmr$annee)
  return(residuals(lm))
}

corr_DAC_lm = map_dfc(quanti %>%
                        dplyr::select(-newRFI),
                      fct_correct_fixedEffects_lm)

corr_DAC_lm = corr_DAC_lm %>%
  add_column("newRFI" = quanti$newRFI)

corr_Buckets = map_dfc(clr_buckets, fct_correct_fixedEffects_lm)

#### 8. Standardisation ####

scaled_DAC = scale(corr_DAC_lm)
scaled_DAC = as_tibble(scaled_DAC)

scaled_not_cor = scale(quanti)
scaled_not_cor = as_tibble(scaled_not_cor)

#scaled_buckets = scale(corr_Buckets)
#scaled_buckets = as_tibble(scaled_buckets)

#### 9. Bucket selection ####

#Elimination 25% or 50% buckets according to variation coefficient
coeff_var = function (x) {(sd(x)/mean(x))*100}

select_buckets = reshape2::melt(purrr::map_dfc(imputed_buckets, coeff_var))%>%
  arrange(value)

select_buckets_25 = select_buckets[round(length(select_buckets$variable)*0.25):length(select_buckets$variable),]
select_buckets_50 = select_buckets[round(length(select_buckets$variable)*0.5):length(select_buckets$variable),]

select_buckets_25 = as.character(select_buckets_25$variable)
select_buckets_50= as.character(select_buckets_50$variable)

#### 10. Divergent RFI animal selection ####

select_rfi = DAC_nmr %>%
  arrange(newRFI)

select_rfi = select_rfi[c(1:round(length(select_rfi$newRFI)*0.25),round(length(select_rfi$newRFI)*0.75):length(select_rfi$newRFI)),]

select_rfi %>%
ggplot(aes(factRFI, newRFI))+
  geom_boxplot(aes(col = factRFI))

select_rfi = select_rfi$nom_anim_RMN

#### 11. Prepare environment for analysis ####

#add column animal
corr_Buckets = corr_Buckets %>%
  add_column("animal" = tbuckets_DAC$animal) %>%
  relocate(animal)

clr_buckets = clr_buckets %>%
  add_column("animal" = tbuckets_DAC$animal) %>%
  relocate(animal)

##### a) df for prediction SC #####
df_pred_SC_V3 = DAC_nmr %>% 
  unite("annee_lignee", c(annee,lignee), sep= "_", 
        remove = FALSE) %>%
  dplyr::select(nom_anim_RMN, annee_lignee, annee, NULOBA_annee, 
                lignee,factRFI, CONMOY, newRFI, FCR, 
                Citrate, `L-Leucine`,`beta-Hydroxyisovalerate`,
                MUFC, GRFC, GMQ, PODEC, POFIC,
                ageDC)

annee_dummy =unmap(df_pred_SC_V3$annee)
colnames(annee_dummy) = levels(df_pred_SC_V3$annee)

df_pred_SC_V3 = cbind(df_pred_SC_V3, annee_dummy)

lot_dummy =unmap(df_pred_SC_V3$NULOBA_annee)
colnames(lot_dummy) = levels(df_pred_SC_V3$NULOBA_annee)

df_pred_SC_V3 = cbind(df_pred_SC_V3, lot_dummy)

df_pred_SC_V3 = merge(
  df_pred_SC_V3%>%rename(animal=nom_anim_RMN),
  clr_buckets,
  by = "animal"
)

##### b) df for prediction corr #####
df_pred_corr_V3 = DAC_nmr %>% 
  unite("annee_lignee", c(annee,lignee), sep= "_", 
        remove = FALSE) %>%
  dplyr::select(nom_anim_RMN, annee_lignee, annee, NULOBA_annee, 
                lignee, factRFI, newRFI, 
                ageDC)

df_pred_corr_V3 =  cbind(df_pred_corr_V3, corr_DAC_lm %>%
                           dplyr::select(CONMOY, FCR, Citrate,`L-Leucine`,`beta-Hydroxyisovalerate`, MUFC, GRFC, GMQ, PODEC, POFIC)) %>% 
  relocate(.after = factRFI, CONMOY)%>%
  relocate(.after = newRFI, FCR)%>%
  relocate(.after =POFIC, ageDC)

df_pred_corr_V3 = cbind(df_pred_corr_V3, annee_dummy)
df_pred_corr_V3 = cbind(df_pred_corr_V3, lot_dummy)

df_pred_corr_V3 = merge(
  df_pred_corr_V3%>%rename(animal=nom_anim_RMN),
  corr_Buckets,
  by = "animal"
)

##### c) df corr and SC eliminate buckets 50 and 25% ##### 

#25%
df_pred_corr_25_V3 = cbind(
  df_pred_corr_V3[1:36],
  df_pred_corr_V3 %>% dplyr::select_if(names(df_pred_corr_V3) %in% select_buckets_25)
)

df_pred_SC_25_V3 = cbind(
  df_pred_SC_V3[1:36],
  df_pred_SC_V3 %>% dplyr::select_if(names(df_pred_SC_V3) %in% select_buckets_25)
)

#50%
df_pred_corr_50_V3 = cbind(
  df_pred_corr_V3[1:36],
  df_pred_corr_V3 %>% dplyr::select_if(names(df_pred_corr_V3) %in% select_buckets_50)
)

df_pred_SC_50_V3 = cbind(
  df_pred_SC_V3[1:36],
  df_pred_SC_V3 %>% dplyr::select_if(names(df_pred_SC_V3) %in% select_buckets_50)
)


##### d) df corr and SC for discriminant RFI #####
df_pred_corr_divRFI_V3 = df_pred_corr_V3 %>%
  filter(animal %in% select_rfi)

df_pred_SC_divRFI_V3 = df_pred_SC_V3 %>%
  filter(animal %in% select_rfi)

#### 12. write tables ####

ifelse(!dir.exists("for_cluster/V3/"), dir.create("for_cluster/V3/", recursive=TRUE), "folder exists already")

#All
write_tsv(df_pred_corr_V3, file = "for_cluster/V3/df_corr_V3")
write_tsv(df_pred_SC_V3, file = "for_cluster/V3/df_SC_V3")

#25
write_tsv(df_pred_corr_25_V3, file = "for_cluster/V3/df_corr_25_V3")
write_tsv(df_pred_SC_25_V3, file = "for_cluster/V3/df_SC_25_V3")

#50
write_tsv(df_pred_corr_50_V3, file = "for_cluster/V3/df_corr_50_V3")
write_tsv(df_pred_SC_50_V3, file = "for_cluster/V3/df_SC_50_V3")

#divRFI
write_tsv(df_pred_corr_divRFI_V3, file = "for_cluster/V3/df_corr_rfi_V3")
write_tsv(df_pred_SC_divRFI_V3, file = "for_cluster/V3/df_SC_rfi_V3")