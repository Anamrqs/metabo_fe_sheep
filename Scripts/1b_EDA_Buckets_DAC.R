#==========================#
# EDA Buckets              #
# Marquisseau Anaïs        #
#==========================#

# 2022
# Internship INRAE M2 BI
# Relationship between feed efficiency and metabolomics in sheep

#### Goals ####

# Description of metabolomic dataset


#### Libraries ####

library(tidyverse)
library(ggpubr)
library(plotly)

#### Files ####

source("0_Cleaning.R")


#Plot raw
#V1
ggplot(buckets_DAC_v1_long, aes(ppm,value, color=variable)) + 
  geom_line(show.legend = FALSE) +
  theme_minimal()

#V2
ggplot(buckets_DAC_long, aes(ppm,value, color=variable)) + 
  geom_line(show.legend = FALSE) +
  theme_minimal()

#All in one
ggarrange(ggplot(buckets_DAC_v1_long, aes(ppm,value, color=variable)) + 
            geom_line(show.legend = FALSE) +
            theme_minimal(),
          
          ggplot(buckets_DAC_long, aes(ppm,value, color=variable)) + 
            geom_line(show.legend = FALSE) +
            theme_minimal(),
          
          labels = c("DAC before baseline",  "DAC after baseline"),
          ncol = 2, nrow = 1)

##### where are zeros #####
ggarrange(
  ggplot(buckets_DAC_v1_long, aes(variable, ppm)) +
    geom_tile(aes(fill = as.factor(value<=0)), show.legend = FALSE)+
    scale_fill_manual(values=c("#999999", "#E83A14"))
  ,
  ggplot(buckets_DAC_long, aes(variable, ppm)) +
    geom_tile(aes(fill = as.factor(value<=0)), show.legend = FALSE)+
    scale_fill_manual(values=c("#999999", "#E83A14"))
  , labels = c("DAC v1",  "DAC v2"),
  ncol = 2, nrow = 1)


#### Caracterisation des zéros dans buckets ####
View(buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value) %>%
  arrange(desc(n)))

#DAC
#100%
a= buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value)%>%
  filter(n==ncol(buckets_DAC_filtered)-1)%>%
  ungroup()%>%
  summarise(n = n())

#>85%
b= buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value)%>%
  filter(n>0.85*(ncol(buckets_DAC_filtered)-1))%>%
  ungroup()%>%
  summarise(n = n())

#>50%
c= buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value)%>%
  filter(n>0.5*(ncol(buckets_DAC_filtered)-1))%>%
  ungroup()%>%
  summarise(n = n())

#<15%
e= buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value)%>%
  filter(n<0.15*(ncol(buckets_DAC_filtered)-1))%>%
  ungroup()%>%
  summarise(n = n())

#<50
d= buckets_DAC_long %>%
  filter(value == 0) %>%
  group_by(ppm) %>%
  count(value)%>%
  filter(n<0.5*(ncol(buckets_DAC_filtered)-1))%>%
  ungroup()%>%
  summarise(n = n())

prop_zeros_DAC = data.frame(
  proportion_zeros = c("100%","100 > n > 85%","85 > n > 50%","50 > n > 15%", "< 15%"), 
  
  nombre_buckets = c( a[[1]], b[[1]]-a[[1]],c[[1]]-b[[1]],d[[1]]-e[[1]],e[[1]])
)

prop_zeros_DAC

rm(a,b,c,d,e)



#### Prop zeros

#rfi+
prop_zeros_per_bucket_rfiP = buckets_DAC_long %>%
  filter(lignee == "rfi+") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(desc(prop_zeros))

prop_zeros_per_bucket_rfiP$ppm = factor(prop_zeros_per_bucket_rfiP$ppm, levels = prop_zeros_per_bucket_rfiP$ppm)

#rfi-
prop_zeros_per_bucket_rfiN = buckets_DAC_long %>%
  filter(lignee == "rfi-")%>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(desc(prop_zeros))

prop_zeros_per_bucket_rfiN$ppm = factor(prop_zeros_per_bucket_rfiN$ppm, levels = prop_zeros_per_bucket_rfiP$ppm)

#filter when there are no 0
prop_zeros_per_bucket = full_join(
  prop_zeros_per_bucket_rfiP,
  prop_zeros_per_bucket_rfiN,
  by ="ppm"
) %>%
  dplyr::select(ppm, prop_zeros.x, prop_zeros.y) %>%
  filter(prop_zeros.x != 0 & prop_zeros.y != 0)

#plot prop zeros all ppm 
ggplot() +
  geom_point(data =prop_zeros_per_bucket_rfiP, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "blue")+
  
  geom_point(data =prop_zeros_per_bucket_rfiN, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "red")+
  theme_classic() 

#plot prop zeros filtered ppm
ggplot(data=prop_zeros_per_bucket) +
  geom_point( aes(
    x = ppm,
    y= prop_zeros.x, 
  ), col = "blue")+
  
  geom_point(aes(
    x = ppm,
    y= prop_zeros.y, 
  ), col = "red")+
  theme_minimal() 

#plot rearranged by bucket

#rfi+
prop_zeros_per_bucket_rfiP = buckets_DAC_long %>%
  filter(lignee == "rfi+") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(ppm)

prop_zeros_per_bucket_rfiP$ppm = factor(prop_zeros_per_bucket_rfiP$ppm, levels = prop_zeros_per_bucket_rfiP$ppm)

#rfi-
prop_zeros_per_bucket_rfiN = buckets_DAC_long %>%
  filter(lignee == "rfi-")%>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100)

prop_zeros_per_bucket_rfiN$ppm = factor(prop_zeros_per_bucket_rfiN$ppm, levels = prop_zeros_per_bucket_rfiP$ppm)

#filter when there are no 0
prop_zeros_per_bucket = full_join(
  prop_zeros_per_bucket_rfiP,
  prop_zeros_per_bucket_rfiN,
  by ="ppm"
) %>%
  dplyr::select(ppm, prop_zeros.x, prop_zeros.y) %>%
  filter(prop_zeros.x != 0 & prop_zeros.y != 0)

ggplotly(
ggplot(data=prop_zeros_per_bucket) +
  geom_point( aes(
    x = ppm,
    y= prop_zeros.x, 
  ), col = "blue")+
  
  geom_point(aes(
    x = ppm,
    y= prop_zeros.y, 
  ), col = "red")+
  theme_minimal() 
)

ggplot(data=prop_zeros_per_bucket, aes(x=ppm, group=1)) +
  geom_line( aes(
    x = ppm,
    y= prop_zeros.x, 
  ), col = "blue")+
  
  geom_line(aes(
    x = ppm,
    y= prop_zeros.y, 
  ), col = "red")+
  
  geom_point(aes(y= prop_zeros.y), col = "red")+
  geom_point(aes(y= prop_zeros.x), col = "blue")+
  theme_minimal() 

##prop zero par bucket par annee
#2018
prop_zeros_per_bucket_2018 = buckets_DAC_long %>%
  filter(annee == "2018") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(desc(prop_zeros))

prop_zeros_per_bucket_2018$ppm = factor(prop_zeros_per_bucket_2018$ppm, levels = prop_zeros_per_bucket_2018$ppm)

#2019
prop_zeros_per_bucket_2019 = buckets_DAC_long %>%
  filter(annee == "2019") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(desc(prop_zeros))

prop_zeros_per_bucket_2019$ppm = factor(prop_zeros_per_bucket_2019$ppm, levels = prop_zeros_per_bucket_2018$ppm)

#2020
prop_zeros_per_bucket_2020 = buckets_DAC_long %>%
  filter(annee == "2020") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(desc(prop_zeros))

prop_zeros_per_bucket_2020$ppm = factor(prop_zeros_per_bucket_2020$ppm, levels = prop_zeros_per_bucket_2018$ppm)

prop_zeros_per_bucket_annee = full_join(full_join(
  prop_zeros_per_bucket_2018,
  prop_zeros_per_bucket_2019,
  by ="ppm"),
  prop_zeros_per_bucket_2020,
  by ="ppm"
) %>%
  dplyr::select(ppm, prop_zeros.x, prop_zeros.y, prop_zeros) %>%
  filter(prop_zeros.x != 0 & prop_zeros.y != 0 & prop_zeros != 0 )

#Plot all
ggplot() +
  geom_point(data =prop_zeros_per_bucket_2018, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "blue")+
  
  geom_point(data =prop_zeros_per_bucket_2019, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "red")+
  geom_point(data =prop_zeros_per_bucket_2020, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "dark green")+
  theme_classic() 

#plot without same prop 0% zeros
ggplotly(
  ggplot(data=prop_zeros_per_bucket_annee) +
    geom_point( aes(
      x = ppm,
      y= prop_zeros.x, 
    ), col = "blue")+
    
    geom_point(aes(
      x = ppm,
      y= prop_zeros.y, 
    ), col = "red")+
    
    geom_point(aes(
      x = ppm,
      y= prop_zeros, 
    ), col = "dark green")+
    theme_minimal() 
)

#plot rearranged by bucket
#2018
prop_zeros_per_bucket_2018 = buckets_DAC_long %>%
  filter(annee == "2018") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) %>%
  arrange(ppm)

prop_zeros_per_bucket_2018$ppm = factor(prop_zeros_per_bucket_2018$ppm, levels = prop_zeros_per_bucket_2018$ppm)

#2019
prop_zeros_per_bucket_2019 = buckets_DAC_long %>%
  filter(annee == "2019") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) 

prop_zeros_per_bucket_2019$ppm = factor(prop_zeros_per_bucket_2019$ppm, levels = prop_zeros_per_bucket_2018$ppm)

#2020
prop_zeros_per_bucket_2020 = buckets_DAC_long %>%
  filter(annee == "2020") %>%
  group_by(ppm) %>%
  summarise( n=n(), n_zeros = sum(value == 0)) %>% 
  mutate(prop_zeros = (n_zeros/n)*100) 

prop_zeros_per_bucket_2020$ppm = factor(prop_zeros_per_bucket_2020$ppm, levels = prop_zeros_per_bucket_2018$ppm)

prop_zeros_per_bucket_annee = full_join(full_join(
  prop_zeros_per_bucket_2018,
  prop_zeros_per_bucket_2019,
  by ="ppm"),
  prop_zeros_per_bucket_2020,
  by ="ppm"
) %>%
  dplyr::select(ppm, prop_zeros.x, prop_zeros.y, prop_zeros) %>%
  filter(prop_zeros.x != 0 & prop_zeros.y != 0 & prop_zeros != 0 )


#plot without same prop 0% zeros
ggplotly(
  ggplot(data=prop_zeros_per_bucket_annee) +
    geom_point( aes(
      x = ppm,
      y= prop_zeros.x, 
    ), col = "blue")+
    
    geom_point(aes(
      x = ppm,
      y= prop_zeros.y, 
    ), col = "red")+
    
    geom_point(aes(
      x = ppm,
      y= prop_zeros, 
    ), col = "dark green")+
    theme_minimal() 
)

ggplotly(
ggplot(data=prop_zeros_per_bucket_annee, aes(x=ppm, group=1)) +
  geom_line( aes(
    y= prop_zeros.x
  ), color = "blue")+
  
  geom_line(aes(
    y= prop_zeros.y
  ), color = "red")+
  
  geom_line(aes(
    y= prop_zeros
  ), color = "dark green")+
  
  geom_point( aes(
    y= prop_zeros.x
  ), color = "blue")+
  
  geom_point(aes(
    y= prop_zeros.y
  ), color = "red")+
  
  geom_point(aes(
    y= prop_zeros
  ), color = "dark green")+
  
   theme_minimal()
)

ggplot() +
  geom_line(data =prop_zeros_per_bucket_2018, aes(
    x = ppm,
    y= prop_zeros, 
    group=1
  ), color = "blue")+
  
  geom_line(data =prop_zeros_per_bucket_2019, aes(
    x = ppm,
    y= prop_zeros,
    group=1
  ), color = "red")+
  geom_line(data =prop_zeros_per_bucket_2020, aes(
    x = ppm,
    y= prop_zeros, 
    group=1 
  ), color = "dark green")+
  geom_point(data =prop_zeros_per_bucket_2018, aes(
    x = ppm,
    y= prop_zeros, 
  ), col= "blue")+
  
  geom_point(data =prop_zeros_per_bucket_2019, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "red")+
  geom_point(data =prop_zeros_per_bucket_2020, aes(
    x = ppm,
    y= prop_zeros, 
  ), col = "dark green")+
  theme_classic()

#### Correction effets : lesquels prendre en compte ? ####
#Which effects affect metabolomic variables ?

#Input : a dataframe of predictors, an effect to test (lm)
#Output : a df with a pvalue, r2 and 1-log(pvalue) for each predictor
recup_pval_r2 = function (X,y) {
  
  r2 = apply(X,2, function (x) {summary(lm(x~y))$adj.r.squared})
  pval = apply(X,2, function (x) {pf(summary(lm(x~y))$fstatistic[1],summary(lm(x~y))$fstatistic[2],summary(lm(x~y))$fstatistic[3],lower.tail=F)[[1]]})
  pval_adj = apply(X,2, function (x) {p.adjust(pf(summary(lm(x~y))$fstatistic[1],summary(lm(x~y))$fstatistic[2],summary(lm(x~y))$fstatistic[3],lower.tail=F)[[1]], method = "BH", n = 877)})
  df = data.frame (noms=names(X), r2_adj = as.vector(r2), pvalue = as.vector(pval), pval_adj = as.vector(pval_adj), log_pval = as.vector(1-log(pval_adj)) )
  
  return(df)
}

X = clr_buckets %>%
  dplyr::select(-animal)

DAC_nmr = DAC_nmr %>%
  add_column(regime = ifelse(DAC_nmr$annee != "2019", "J", "NJ" ))

effet_annee = recup_pval_r2(X, DAC_nmr$annee) %>% add_column(effet = rep("annee", ncol(X)))
effet_NULOBA = recup_pval_r2(X, DAC_nmr$NULOBA_annee) %>% add_column(effet = rep("lot", ncol(X)))
effet_age = recup_pval_r2(X, DAC_nmr$ageDC) %>% add_column(effet = rep("age", ncol(X)))
effet_lignee = recup_pval_r2(X, DAC_nmr$lignee) %>% add_column(effet = rep("lignee", ncol(X)))
effet_MEALLA = recup_pval_r2(X, DAC_nmr$MEALLA) %>% add_column(effet = rep("MEALLA", ncol(X)))
effet_COMNAI = recup_pval_r2(X, DAC_nmr$COMNAI) %>% add_column(effet = rep("COMNAI", ncol(X)))
effet_mnal = recup_pval_r2(X, DAC_nmr$mnal) %>% add_column(effet = rep("mnal", ncol(X)))
effet_COMOEL = recup_pval_r2(X, DAC_nmr$COMOEL) %>% add_column(effet = rep("COMOEL", ncol(X)))
effet_regime = recup_pval_r2(X, DAC_nmr$regime) %>% add_column(effet = rep("regime", ncol(X)))

effets = rbind (
  effet_annee, effet_NULOBA, effet_age, effet_lignee,effet_MEALLA,effet_COMNAI,effet_mnal,effet_COMOEL, effet_regime 
)

effets %>%
  filter(effet %in% c("annee", "lot", "age"))%>%
  ggplot(aes(noms, log_pval, color = effet, group =1)) +
  geom_point()+ 
  geom_path()+
  geom_hline(yintercept = 1-log(0.05), linetype = "dashed")+
  labs(x = "buckets (ppm)", y = "1-log(p-value)", color = "Effet")+
  ggtitle("Effets sur les variables zootechniques") +
  theme_void()

ggplotly(
effets %>%
  filter(effet %in% c("annee", "lot", "age"))%>%
  ggplot(aes(noms, r2_adj, color = effet, group =1)) +
  geom_point()+ 
  geom_path()+
  labs(x = "buckets (ppm)", y = "r2 ajusté", color = "Effet")+
  ggtitle("Effets annee, lot et age sur le metabolome") +
  theme_void()
)


r2 = apply(X,2, function (x) {summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$adj.r.squared})
pval = apply(X,2, function (x) {pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[3],lower.tail=F)[[1]]})
pval_adj = apply(X,2, function (x) {p.adjust(pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee))$fstatistic[3],lower.tail=F)[[1]], method = "BY", n = 877)})

effet_annee_NULOBA = data.frame (noms=names(X), r2_adj = as.vector(r2), pvalue = as.vector(pval), pval_adj = as.vector(pval_adj), log_pval = as.vector(1-log(pval_adj))) %>% add_column(effet = rep("annee_lot", ncol(X)))

r2 = apply(X,2, function (x) {summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$adj.r.squared})
pval = apply(X,2, function (x) {pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[3],lower.tail=F)[[1]]})
pval_adj = apply(X,2, function (x) {p.adjust(pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC))$fstatistic[3],lower.tail=F)[[1]], method = "BY", n = 877)})


effet_annee_NULOBA_age = data.frame (noms=names(X), r2_adj = as.vector(r2), pvalue = as.vector(pval), pval_adj = as.vector(pval_adj), log_pval = as.vector(1-log(pval_adj)) ) %>% add_column(effet = rep("annee_lot_age", ncol(X)))

r2 = apply(X,2, function (x) {summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$adj.r.squared})
pval = apply(X,2, function (x) {pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[3],lower.tail=F)[[1]]})
pval_adj = apply(X,2, function (x) {p.adjust(pf(summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[1],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[2],summary(lm(x~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime))$fstatistic[3],lower.tail=F)[[1]], method = "BY", n = 877)})

View(lm(X$`0.505`~DAC_nmr$annee+DAC_nmr$NULOBA_annee%in%DAC_nmr$annee+DAC_nmr$ageDC+DAC_nmr$regime)[3]$effects)
summary(lm(X$`0.615`~DAC_nmr$regime+DAC_nmr$annee))

effet_annee_NULOBA_age_regime = data.frame (noms=names(X), r2_adj = as.vector(r2), pvalue = as.vector(pval), pval_adj = as.vector(pval_adj), log_pval = as.vector(1-log(pval_adj)) ) %>% add_column(effet = rep("annee_lot_age_regime", ncol(X)))


effets_buckets = list("annee" = effet_annee,
                      "lot" =effet_NULOBA,
                      "age" = effet_age,
                      "lignee" =effet_lignee,
                      "mode allaitement" =effet_MEALLA ,
                      "mode naissance" =effet_COMNAI , 
                      "synthese mnal" =effet_mnal ,
                      "mode elevage" =effet_COMOEL,
                      "regime" = effet_regime,
                      "annee + lot" =effet_annee_NULOBA, 
                      "annee, lot et age" =effet_annee_NULOBA_age, 
                      "annee_NULOBA_age_regime" = effet_annee_NULOBA_age_regime)
effets = rbind (
  effet_annee, effet_NULOBA, effet_age, effet_lignee,effet_MEALLA,effet_COMNAI,effet_mnal,effet_COMOEL, effet_regime,
  effet_annee_NULOBA, effet_annee_NULOBA_age, effet_annee_NULOBA_age_regime
)

effets %>%
  filter(effet %in% c("annee", "annee_lot", "annee_lot_age"))%>%
  ggplot(aes(noms, log_pval, color = effet, group =1)) +
  geom_point()+ 
  geom_path()+
  geom_hline(yintercept = 1-log(0.05), linetype = "dashed")+
  labs(x = "buckets (ppm)", y = "1-log(p-value)", color = "Effet")+
  ggtitle("Effets sur les variables zootechniques") +
  theme_void()

effets %>%
  filter(effet %in% c("annee", "annee_lot", "annee_lot_age"))%>%
  dplyr::select(r2_adj)
  

Y = clr_buckets$`1.925`
lm0 =lm(Y ~ DAC_nmr$annee )
lm1 = lm(Y ~ DAC_nmr$annee +DAC_nmr$NULOBA_annee %in% DAC_nmr$annee)
lm2 = lm(Y ~ DAC_nmr$annee + DAC_nmr$NULOBA_annee %in% DAC_nmr$annee + DAC_nmr$regime + DAC_nmr$ageDC)
summary(lm1)
plot(lm1 )

anova(lm0,lm1,lm2)

y = clr_buckets$`1.925`
a =anova(lm2)
a$`Pr(>F)`


recup_pval_anova = function(y) {
  lm = lm(y ~ DAC_nmr$annee + DAC_nmr$NULOBA_annee %in% DAC_nmr$annee + DAC_nmr$ageDC + DAC_nmr$regime)
  a = anova(lm)
  recup = p.adjust(a$`Pr(>F)`, method = "BH", n = 877)
  return(recup)
}

res_anova = purrr::map_dfc(X, recup_pval_anova)
view(res_anova)
plot(res_anova[1,])

res_anova = cbind(
gather(res_anova[1,], "bucket", "annee"),
gather(res_anova[2,], "bucket", "lot")[,2],
gather(res_anova[3,], "bucket", "ageDC")[,2],
gather(res_anova[4,], "bucket", "regime")[,2]
)

#Regime est NA
res_anova %>%
  filter(!is.na(regime))

#plot
res_anova %>%
  filter(annee < 0.05) %>%
  count()

res_anova %>%
  filter(lot < 0.05) %>%
  count()

res_anova %>%
  filter(ageDC < 0.05) %>%
  count()



an =anova(lm2)
an2=anova(lm0, lm1, lm2)
an$`Pr(>F)`
an 
an2

summary(lm1)$adj.r.squared
summary(lm2)$adj.r.squared

  
View(effets_buckets$annee)
ifelse(!dir.exists("Output/DAC/EDA_Buckets/Effets_buckets"), dir.create("Output/DAC/EDA_Buckets/Effets_buckets", recursive=FALSE), "folder exists already")

#Plot pvalue and R2
pdf("Output/DAC/EDA_Buckets/Effets_buckets/Effets_buckets_pvaladjR2.pdf", paper = "USr", width = 11, height = 8)
lapply(1:length(effets_buckets), function (x) { 
  df = as.data.frame(effets_buckets[x])
  names(df) = c("noms", "r2_adj", "pvalue", "pvalue_adj", "log_pval", "effet")
  signif = df %>% filter(pval_adj < 0.05) %>% count()
  df %>%
    ggplot(aes(as.numeric(noms),pval_adj))+
    geom_point(aes(color = ifelse(pval_adj < 0.05, paste("pvalue < 0.05 (", round(100*(signif$n[1])/877), "%)", sep=""), "pvalue > 0.05"),
                   size = r2_adj))+
    labs(x = "buckets (ppm)", y = "p-value", color = "P-value", size = "R2 ajusté")+
    ggtitle(paste("Effet ",names(effets_buckets)[x]," sur les buckets", sep=""))
})
dev.off()

#Plot 1-log(pvalue) and R2
pdf("Output/DAC/EDA_Buckets/Effets_buckets/Effets_buckets_pvaladjlogR2.pdf", paper = "USr", width = 11, height = 8)
lapply(1:length(effets_buckets), function (x) { 
  df = as.data.frame(effets_buckets[x])
  names(df) = c("noms", "r2_adj", "pvalue", "pval_adj", "log_pval", "effet")
  signif = df %>% filter(pval_adj < 0.05) %>% count()
  df %>%
    ggplot(aes(as.numeric(noms),log_pval))+
    geom_point(aes(color = ifelse(log_pval > 1-log(0.05), paste("pvalue < 0.05 (", round(100*(signif$n[1])/877), "%)", sep=""), "pvalue > 0.05"),
                   size = r2_adj))+
    labs(x = "buckets (ppm)", y = "1-log(p-value)", color = "P-value", size = "R2 ajusté")+
    ggtitle(paste("Effet ",names(effets_buckets)[x]," sur les buckets", sep=""))
})
dev.off()

#Plot 1-log(pvalue) and line
pdf("Output/DAC/EDA_Buckets/Effets_buckets/Effets_buckets_logLine.pdf", paper = "USr", width = 11, height = 8)
lapply(1:length(effets_buckets), function (x) { 
  df = as.data.frame(effets_buckets[x])
  names(df) = c("noms", "r2_adj", "pvalue", "pval_adj", "log_pval", "effet")
  signif = df %>% filter(pval_adj < 0.05) %>% count()
  df %>%
    ggplot(aes(as.numeric(noms),log_pval))+
    geom_point(aes(color = ifelse(log_pval > 1-log(0.05), paste("pvalue < 0.05 (", round(100*(signif$n[1])/877), "%)", sep=""), "pvalue > 0.05")))+
    geom_line()+
    labs(x = "buckets (ppm)", y = "1-log(p-value)", color = "P-value")+
    ggtitle(paste("Effet ",names(effets_buckets)[x]," sur les buckets", sep=""))
})
dev.off()

#plotly for interactive plots

recup_pval_anova_lm = function(y) {
  lm0 = lm(y ~ DAC_nmr$annee)
  lm1 = lm(y ~ DAC_nmr$annee + DAC_nmr$NULOBA_annee %in% DAC_nmr$annee)
  lm2 = lm(y ~ DAC_nmr$annee + DAC_nmr$NULOBA_annee %in% DAC_nmr$annee + DAC_nmr$ageDC)
  lm3 = lm(y ~ DAC_nmr$annee + DAC_nmr$NULOBA_annee %in% DAC_nmr$annee + DAC_nmr$ageDC + DAC_nmr$regime)
  a = anova(lm0,lm1,lm2,lm3)
  recup = a$`Pr(>F)`
  return(recup)
}

res_anova = purrr::map_dfc(X, recup_pval_anova_lm)
res_anova[1:4,1:4]
res_anova = cbind(
  gather(res_anova[1,], "bucket", "model_0"),
  gather(res_anova[2,], "bucket", "model_1")[,2],
  gather(res_anova[3,], "bucket", "model_2")[,2],
  gather(res_anova[4,], "bucket", "model_3")[,2]
)

view(res_anova)

res_anova %>%
  filter(p.adjust(model_1, method = "BH", n=length(model_1)) < 0.05)%>%
  summarise(n=n())

res_anova %>%
  filter(p.adjust(model_2, method = "BH", n=length(model_1)) < 0.05)%>%
  summarise(n=n())

res_anova %>%
  filter(p.adjust(model_3, method = "BH", n=length(model_1)) < 0.05)%>%
  summarise(n=n())


#### ACP ####
quali = quali %>%
  dplyr::select (-num_passage, -NULOBA)

tbuckets_DAC[,2:878]


library(mixOmics)

pca = pca(
  tbuckets_DAC[,2:878], 
  scale = FALSE
)

plotIndiv(pca)

apply(quali, 2, function (x) {
plotIndiv(
  pca,
  legend = TRUE,
  group = x
)
})


#### 1. PCA ####

ifelse(!dir.exists("Output/DAC/EDA_Buckets/PCA"), dir.create("Output/DAC/EDA_Buckets/PCA", recursive=TRUE), "folder exists already")

###### i. raw DAC (all) ######
pca_raw_DAC = pca (tbuckets_DAC[,2:878], center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_raw_DAC_2.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_raw_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA raw DAC Buckets", 
                                            ellipse = F))
dev.off()

plotVar(pca_raw_DAC)
plotVar(pca_scaled_DAC)

###### ii. imputed DAC (all) ######
pca_imp_DAC = pca (imputed_buckets, center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_imp_DLcol_DAC_2.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_imp_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA corr DAC Buckets", 
                                            ellipse =F))
dev.off()


###### iii. clr DAC (all) ######
pca_clr_DAC = pca (clr_buckets, center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_clr_DAC_2.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_clr_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA corr DAC Buckets", 
                                            ellipse =F))
dev.off()

###### iv. corr DAC (all) ######

#lm
pca_corr_DAC_lm = pca (corr_Buckets_lm, center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_corr_DAC_lm.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_corr_DAC_lm, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA corr DAC Buckets lm", 
                                            ellipse =F))
dev.off()


#lmrob
pca_corr_DAC = pca (corr_Buckets, center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_corr_DAC_2.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_corr_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA corr DAC Buckets", 
                                            ellipse =F))
dev.off()


###### v. scaled DAC (all) ###### 
pca_scaled_DAC = pca (scaled_buckets, center = FALSE, scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PCA/pca_scaled_DAC_2.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_scaled_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA scaled DAC Buckets", 
                                            scale = FALSE,
                                            ellipse = F))
dev.off()

#### 1. PLSDA ####

ifelse(!dir.exists("Output/DAC/EDA_Buckets/PLSDA"), dir.create("Output/DAC/EDA_Buckets/PLSDA", recursive=TRUE), "folder exists already")

###### i. raw DAC (all) ######
plsda_raw_DAC = plsda (tbuckets_DAC[,2:878], 
                       group = quali$lignee, 
                       scale = FALSE)

pdf("Output/DAC/EDA_Buckets/PLSDA/plsda_raw_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda (X = tbuckets_DAC[,2:878], 
                                                   Y = quali[[x]], 
                                                   scale = FALSE), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PLSDA raw DAC Buckets", 
                                            ellipse = TRUE))
dev.off()
tbuckets_DAC[46,1]
DAC_nmr[46,]


plotVar(pca_raw_DAC)
plotVar(pca_scaled_DAC)

###### ii. corr DAC (all) ###### 
pdf("Output/DAC/EDA_Buckets/PLSDA/plsda_corr_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda (X = corr_Buckets, 
                                                   Y = quali[[x]], 
                                                   scale = FALSE),
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PLSDA corr DAC Buckets", 
                                            ellipse = TRUE))
dev.off()

###### iii. scaled DAC (all) ###### 
pdf("Output/DAC/EDA_Buckets/PLSDA/plsda_scaled_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda (X = scaled_buckets, 
                                                   Y = quali[[x]], 
                                                   scale = FALSE),
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PLSDA scaled DAC Buckets", 
                                            ellipse = TRUE))
dev.off()



###### V3 ######
#Reshape df long version
buckets_DAC_v3_long = reshape2::melt(buckets_DAC_filtered, id.vars="ppm")
buckets_DAC_long = reshape2::melt(buckets_DAC_filtered, id.vars="ppm")

#### test tSNE ####
library(M3C)

pp = tsne(df_pred_SC_V3[37:ncol(df_pred_SC_V3)])
plot(pp$data, col=df_pred_SC_V3$annee, pch=20)
plot(pp$data, col=df_pred_SC_V3$lignee, pch=20)

