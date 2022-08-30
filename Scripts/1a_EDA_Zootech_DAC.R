#=====================#
# Zootechny DAC - EDA #
# Marquisseau Anaïs   #
#=====================#

# 2022
# Internship INRAE M2 BI
# Relationship between feed efficiency and metabolomics in sheep

#### Goals ####

# Description of zootechnical dataset
# Output : pdf files, Ctrl+F to search variables of interest

#### 1. Load libraries and data ####

##### a) Libraries #####

library(corrplot) #correlation plot
library(GGally)
library(rcompanion) #links between variables
library(corrr) #links between variables
#library(plotly) #interactive plots
library(mixOmics) # multivariate analysis # /!\ need to specify dplyr::select or mixomics::select

##### b) Data #####

#source("0_Cleaning.R")

#### 2. Data description ####

##### a) Summary #####

summary(DAC_nmr)
names(DAC_nmr)

# /!\ cv for variables with negative values are false (rfi and index)
DAC_nmr %>%
  dplyr::select(animal,names(quanti))%>%
  gather(variable, value, ageDC:`beta-Hydroxyisovalerate`, factor_key=TRUE)%>%
  group_by(variable)%>%
  summarise(min = min(value), mean = mean(value), max = max(value), sd=sd(value), CV = coeff_var(value))

#per line
DAC_nmr %>%
  dplyr::select(animal,lignee, names(quanti))%>%
  gather(variable, value, ageDC:`beta-Hydroxyisovalerate`, factor_key=TRUE)%>%
  group_by(variable, lignee)%>%
  summarise(min = min(value), mean = mean(value), max = max(value), sd=sd(value), CV = coeff_var(value))

#To check a specific variable (ex : "lignee") : summary(DAC$lignee)

##### b) Factors #####

ifelse(!dir.exists("Output/DAC/EDA_Zootech"), dir.create("Output/DAC/EDA_Zootech", recursive=TRUE), "folder exists already")

#-- barplots  
quali = DAC_nmr %>%
  select_if(is.factor)%>%
  dplyr::select(-num_passage, -NULOBA)

pdf("Output/DAC/EDA_Zootech/barplots_zootech_DAC.pdf")
lapply(1:ncol(quali), 
       function(x) quali %>%
         ggplot(aes(quali[[x]]))+
            geom_bar()+
            geom_text(
              aes(label=after_stat(..count..)),
              stat="count",
              nudge_y=2) + #position count
            xlab(names(quali)[x])+
            ggtitle(names(quali)[x])+
            theme_minimal()
)
dev.off()

##### c) Numeric variables #####

###### i. Histograms ######

pdf("Output/DAC/EDA_Zootech/histograms_DAC.pdf")
lapply(1:ncol(quanti), function(x) hist(quanti[[x]], main=names(quanti)[x], xlab=names(quanti)[x]))
dev.off()

pdf("Output/DAC/EDA_Zootech/hist_density_DAC.pdf")
lapply(1:ncol(quanti), function(x) ggplot(quanti, aes(x= quanti[[x]])) + 
         geom_histogram(aes(y=..density..), fill="white", color="grey")+
         geom_density(alpha=.2, fill="#FF6666")+
         labs(title=names(quanti)[x], x=names(quanti)[x], y="density") +
         theme_minimal())
dev.off()

###### ii. Normality ######

shapiro = quanti %>% apply(2, shapiro.test)

pdf("Output/DAC/EDA_Zootech/qqplots_DAC.pdf")
lapply(1:ncol(quanti), function(x) ggplot(data = quanti)+
         geom_qq(aes(sample = quanti[[x]]))+
         geom_qq_line(aes(sample = quanti[[x]]), color ="red")+
         ggtitle(names(quanti)[x], paste("shapiro p-value =",shapiro[[x]]$p.value))+
         ylab(names(quanti)[x])+
         xlab("Theorique")+
         theme_minimal())
dev.off()

###### iii. Boxplot by factors ######

pdf("Output/DAC/EDA_Zootech/boxplots_quali_all.pdf")
lapply(1:ncol(quali), function(x) lapply(1:ncol(quanti), function(y) 
  ggplot(data = DAC_nmr, aes(quali[[x]], quanti[[y]])) +
    geom_boxplot()+
    geom_jitter(aes(color=quali[[x]]), size=1.5, alpha=0.5) +
    theme_minimal()+
    ggtitle(names(quali[x]), names(quanti[y]))+
    xlab(names(quali[x])) + ylab(names(quanti[y]))+
    labs(color=names(quali[x]))
))
dev.off()

#tests : ttest or wilcox for lines
lapply(1:length(shapiro), 
       function(x) ifelse(shapiro[[x]]$p.value < 0.05, 
                          paste(names(shapiro)[x],"; t.test =", t.test(DAC_nmr %>% 
                                                                                  dplyr::select(names(shapiro)[x])%>% 
                                                                                  as_vector() ~ DAC_nmr$lignee)$p.value), 
                          paste(names(shapiro)[x],"; wilcox.test =",wilcox.test(DAC_nmr %>% 
                                                                                           dplyr::select(names(shapiro)[x])%>% 
                                                                                           as_vector() ~ DAC_nmr$lignee)$p.value)
       ))

lapply(1:length(shapiro), 
       function(x) ifelse(shapiro[[x]]$p.value < 0.05, 
                          paste(names(shapiro)[x],"; t.test =", p.adjust(t.test(DAC_nmr %>% 
                                                                                  dplyr::select(names(shapiro)[x])%>% 
                                                                                  as_vector() ~ DAC_nmr$lignee)$p.value, method="BH")), 
                          paste(names(shapiro)[x],"; wilcox.test =",p.adjust(wilcox.test(DAC_nmr %>% 
                                                                                  dplyr::select(names(shapiro)[x])%>% 
                                                                                  as_vector() ~ DAC_nmr$lignee)$p.value, method="BH"))
))

#test for year
names(shapiro)

lapply(1:length(shapiro), 
       function(x)  {
         kruskal.test(DAC_nmr %>% dplyr::select(names(shapiro)[x])%>% as_vector() ~ DAC_nmr$annee)
       })

lapply(1:length(shapiro), 
       function(x)  {
         pairwise.wilcox.test(DAC_nmr %>% dplyr::select(names(shapiro)[x])%>% as_vector(),DAC_nmr$annee, p.adjust.method = "BH")
       })

#test for pen
names(shapiro)

lapply(1:length(shapiro), 
       function(x)  {
         kruskal.test(DAC_nmr %>% dplyr::select(names(shapiro)[x])%>% as_vector() ~ DAC_nmr$NULOBA_annee)
       })

lapply(1:length(shapiro), 
       function(x)  {
         pairwise.wilcox.test(DAC_nmr %>% dplyr::select(names(shapiro)[x])%>% as_vector(),DAC_nmr$NULOBA_annee, p.adjust.method = "BH")
       })

###### iv. Violin plots per factor ######

pdf("Output/DAC/EDA_Zootech/violin_quali_all.pdf")
lapply(1:ncol(quali), function(x) lapply(1:ncol(quanti), function(y)
  ggplot(data = DAC_nmr, aes(quali[[x]], quanti[[y]])) +
         geom_violin() +
         geom_boxplot(aes(fill=quali[[x]]), width=0.1) +
    theme_minimal()+
    ggtitle(names(quali[x]), names(quanti[y]))+
    xlab(names(quali[x])) + ylab(names(quanti[y]))+
    labs(fill=names(quali[x]))
))
dev.off()

###### v. density per factor ######
pdf("Output/DAC/EDA_Zootech/density_quali_all.pdf")
lapply(1:ncol(quali), function(x) lapply(1:ncol(quanti), function(y)
  ggplot(data = DAC_nmr, aes(quanti[[y]])) + 
         geom_histogram(aes(y=..density.., fill=quali[[x]]), color="grey", position="identity", alpha=0.5)+
         geom_density(alpha=.2, fill="#FF6666")+
    theme_minimal()+
    ggtitle(names(quali[x]), names(quanti[y]))+
    xlab(names(quali[x])) + ylab(names(quanti[y]))+
    labs(fill=names(quali[x]))
))
dev.off()

###### vi. Correlations ######

#-- all numeric

pdf("Output/DAC/EDA_Zootech/corrplot_DAC.pdf")
quanti %>%
  dplyr::select(-`beta-Hydroxyisovalerate`, -Citrate, -`L-Leucine`, -POMET, -age.prelevements, -index)%>%
  cor(use = "complete.obs") %>%
  corrplot(type="upper", order="hclust", tl.col="black", tl.srt=45)
dev.off()

quanti %>%
  dplyr::select(-`beta-Hydroxyisovalerate`, -Citrate, -`L-Leucine`, -POMET, -age.prelevements, -index)%>%
  cor(use = "complete.obs") %>%
  corrplot(type="upper", order="hclust", tl.col="black", tl.srt=45, method = "circle")

#-- selected variables
names(quanti)

pdf("Output/DAC/EDA_Zootech/corrplot_10selectedVars_DAC.pdf", width=10)
quanti %>%
  dplyr::select(-"Citrate",-"L-Leucine",-"beta-Hydroxyisovalerate",-"age.prelevements", -"POMET")%>%
  ggpairs(ggplot2::aes(colour=DAC_nmr$lignee))
dev.off()

pdf("Output/DAC/EDA_Zootech/corrplot_selectedVars_predictors_DAC.pdf", width=10)
DAC_nmr %>%
  dplyr::select(MUFC, GRFC, POFIC, GMQ, FCR, newRFI, CONMOY )%>%
  ggpairs(aes(colour=DAC_nmr$lignee))+
  ggtitle("Correlations between MUFC,GRFC, GMQ, PODEC and POFIC")
dev.off()


####### vii. tables #######

lapply(1:ncol(quali), function(x) lapply(1:ncol(quali), function(y) 
  table(quali[[x]], quali[[y]])
))

pdf("Output/DAC/EDA_Zootech/contingenceTiles_all.pdf")
lapply(1:ncol(quali), function(x) lapply(1:ncol(quali), function(y) quali %>%
         mutate(tmpX = quali[[x]], tmpY = quali[[y]]) %>%
         count(tmpX, tmpY ) %>%  
  ggplot(aes(x = tmpX, y = tmpY)) +
  geom_tile(mapping = aes(fill = n))+
    #scale_fill_gradientn(limits=c(0,268), colours=c("navyblue","darkmagenta", "darkorange1"),breaks=c(0,45,90,134,178,223,268), labels=format(c(0,45,90,134,178,223,268)))+
  geom_text(aes(label= n), col="white")+
    ggtitle(names(quali[x]), names(quali[y]))+
    xlab(names(quali[x]))+
    ylab(names(quali[y]))
))
dev.off()

###### viii. scatterplots ######

pdf("Output/DAC/EDA_Zootech/scatterplot_lines_all.pdf")
lapply(1:ncol(quanti), function(x) lapply(1:ncol(quanti), function(y) 
  ggplot(data = DAC_nmr, aes(quanti[[x]], quanti[[y]])) +
    geom_point(aes(col=DAC_nmr$lignee))+
    ggtitle(names(quanti[x]), names(quanti[y]))+
    xlab(names(quanti[x])) + ylab(names(quanti[y]))+
    labs(color=names(quanti[x]))+
    theme_minimal()
))
dev.off()

pdf("Output/DAC/EDA_Zootech/scatterplot_year_all.pdf")
lapply(1:ncol(quanti), function(x) lapply(1:ncol(quanti), function(y) 
  ggplot(data = DAC_nmr, aes(quanti[[x]], quanti[[y]])) +
    geom_point(aes(col=DAC_nmr$annee))+
    ggtitle(names(quanti[x]), names(quanti[y]))+
    xlab(names(quanti[x])) + ylab(names(quanti[y]))+
    labs(color=names(quanti[x]))+
    theme_minimal()
))
dev.off()

##### d) Association between variables #####

# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
# Adopted from https://stackoverflow.com/a/52557631/590437
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
  # https://github.com/r-lib/rlang/issues/781
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

liaisons_DAC = DAC_nmr %>%
  dplyr::select(-animal, -nuofmg, -NUOFPE, -DAENAT, -DAPMC1,-date_sevr, -date_prelevement, -DAPDC1, -DAPDC2, -DAFIC1, -DAFIC2, -DANAOV, -cd) %>%
  mixed_assoc()

pdf("Output/DAC/EDA_Zootech/liaison_0.3_all.pdf")
DAC_nmr %>%
  dplyr::select(names(quanti), names(quali)) %>%
  dplyr::select( -POMET, -factRFI) %>%
  mixed_assoc() %>%
  dplyr::select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf %>%
  network_plot(min_cor = 0.3, legend = TRUE)
dev.off()

pdf("Output/DAC/EDA_Zootech/liaison_0.5_all.pdf")
DAC_nmr %>%
  dplyr::select(names(quanti), names(quali)) %>%
  dplyr::select(-POMET, -factRFI) %>%
  mixed_assoc() %>%
  dplyr::select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf %>%
  network_plot(min_cor = 0.5, legend = TRUE)
dev.off()

gc()

#### Corrections ####
#Which effects affect zootechnical variables ?

#Input : a dataframe of predictors, an effect to test (lm)
#Output : a df with a pvalue, r2 and 1-log(pvalue) for each predictor
recup_pval_r2 = function (X,y) {
  
  r2 = apply(X,2, function (x) {summary(lm(x~y))$adj.r.squared})
  pval = apply(X,2, function (x) {pf(summary(lm(x~y))$fstatistic[1],summary(lm(x~y))$fstatistic[2],summary(lm(x~y))$fstatistic[3],lower.tail=F)[[1]]})
  df = data.frame (noms=names(X), r2_adj = as.vector(r2), pvalue = as.vector(pval), log_pval = as.vector(1-log(pval)) )
  return(df)
}

X = DAC_nmr %>% 
  dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")

effet_annee = recup_pval_r2(X, DAC_nmr$annee) %>% add_column(effet = rep("annee", ncol(X)))
effet_NULOBA = recup_pval_r2(X, DAC_nmr$NULOBA_annee) %>% add_column(effet = rep("lot", ncol(X)))
effet_age = recup_pval_r2(X, DAC_nmr$ageDC) %>% add_column(effet = rep("age", ncol(X)))
effet_lignee = recup_pval_r2(X, DAC_nmr$lignee) %>% add_column(effet = rep("lignee", ncol(X)))
effet_MEALLA = recup_pval_r2(X, DAC_nmr$MEALLA) %>% add_column(effet = rep("MEALLA", ncol(X)))
effet_COMNAI = recup_pval_r2(X, DAC_nmr$COMNAI) %>% add_column(effet = rep("COMNAI", ncol(X)))
effet_mnal = recup_pval_r2(X, DAC_nmr$mnal) %>% add_column(effet = rep("mnal", ncol(X)))
effet_COMOEL = recup_pval_r2(X, DAC_nmr$COMOEL) %>% add_column(effet = rep("COMOEL", ncol(X)))

effets = rbind (
  effet_annee, effet_NULOBA, effet_age, effet_lignee,effet_MEALLA,effet_COMNAI,effet_mnal,effet_COMOEL
)

effets %>%
  filter(effet %notin% c("lignee", "MEALLA", "COMNAI"))%>%
  ggplot(aes(noms, log_pval, color = effet)) +
  geom_point(aes(size = r2_adj))+  
  geom_hline(yintercept = 1-log(0.05), linetype = "dashed")+
  labs(x = "buckets (ppm)", y = "1-log(p-value)", color = "Effet", size = "R2 ajusté")+
  ggtitle("Effets sur les variables zootechniques") 

effets %>%
  filter(effet %notin% c("lignee", "MEALLA", "COMNAI"))%>%
  ggplot(aes(noms, pvalue, color = effet)) +
    geom_point(aes(size = r2_adj))+  
    geom_hline(yintercept = 0.05, linetype = "dashed")+
    labs(x = "buckets (ppm)", y = "p-value", color = "Effet", size = "R2 ajusté")+
    ggtitle("Effets sur les variables zootechniques") 

#### 3. Multivariate ####
#mixomics

quali = quali %>%
  dplyr::select(-num_passage, -NULOBA)

sel_Vars = scaled_DAC %>%
  dplyr::select(MUFC, GRFC, POFIC, PODEC, GMQ)

# When variable selection : MUFC,GRFC,POFIC,PODEC,GMQ

##### a) PCA #####

ifelse(!dir.exists("Output/DAC/EDA_Zootech/PCA"), dir.create("Output/DAC/EDA_Zootech/PCA", recursive=FALSE), "folder exists already")

###### i. raw DAC (all) ######
pca_raw_DAC = pca (
  quanti %>%
    dplyr::select(names(scaled_DAC)), center = TRUE, scale = TRUE)

pdf("Output/DAC/EDA_Zootech/PCA/pca_raw_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca_raw_DAC, group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "PCA raw DAC", 
                                            ellipse = TRUE))
dev.off()

plotVar(pca_raw_DAC)
plotVar(pca_scaled_DAC)

# Define factors for colours matching plotIndiv above
lignee <- factor(quali$lignee, 
                 levels = c("rfi-", "rfi+"))
annee <- factor(quali$annee, 
                levels = c("2018", "2019", "2020"))

# Set up colours and symbols
col.lignee <- color.mixo(lignee)
pch.annee <- as.numeric(annee)

plot(pca_raw_DAC$variates$X, col = col.lignee, pch = pch.annee, 
     xlab = round((pca_raw_DAC$prop_expl_var$X[1])*100), 
     ylab = round((pca_raw_DAC$prop_expl_var$X[2])*100))
legend("topleft", col = color.mixo(1:2), legend = levels(lignee),
       lty = 1, title = "Lignee")
legend("bottomright", legend = levels(annee), pch = 1:3,
       title = "Annee")
title("PCA zootechnie brute (MUFC,GRFC,POFIC,PODEC,GMQ)")

plot(pca_raw_DAC)
plot(pca_corr)

png("Output/DAC/EDA_Zootech/PCA/pca_raw_DAC_ligneeAnnee.png",  width = 960, height = 960, units = "px")
plotIndiv(pca_raw_DAC, group = DAC_nmr$lignee, 
          pch = DAC_nmr$annee,
          legend = TRUE, legend.title = 'Lignée',
          legend.title.pch = 'Année',
          title = 'ACP données brutes zootechnie')
dev.off()

png("Output/DAC/EDA_Zootech/PCA/pca_corrScaled_DAC_ligneeAnnee.png", width = 960, height = 960, units = "px")
plotIndiv(pca_corr, group = DAC_nmr$lignee, 
          pch = DAC_nmr$annee,
          legend = TRUE, legend.title = 'Lignée',
          legend.title.pch = 'Année',
          title = 'ACP données corrigées zootechnie')
dev.off()

###### ii. scaled DAC (all) ###### 

pdf("Output/DAC/EDA_Zootech/PCA/pca_scaled_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca(scaled_DAC, center = FALSE, scale = FALSE), group=quali[[x]], 
          legend = TRUE, legend.title = names(quali[x]), 
          title = "PCA scaled DAC", 
          scale = FALSE,
          ellipse = TRUE))
dev.off()

###### i. raw DAC (selected) ######

pdf("Output/DAC/EDA_Zootech/PCA/pca_rawSelec_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca(X = quanti %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Pca raw DAC select", 
                                            ellipse = TRUE))
dev.off()

###### ii. corrected DAC (selected) ###### 

#lm
pdf("Output/DAC/EDA_Zootech/PCA/pca_corrSelec_DAC_lm.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca(X = corr_DAC_lm %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Pca corrected DAC select lm", 
                                            ellipse = F))
dev.off()

#lm_rob
pdf("Output/DAC/EDA_Zootech/PCA/pca_corrSelec_DAC_lmrob.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca(X = corr_DAC %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Pca corrected DAC select lmrob", 
                                            ellipse = F))
dev.off()

#lm is sufficient to correct year and pen 
pca_corr = pca(X = corr_DAC_lm %>% 
                 dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), 
               center = TRUE, scale = TRUE)

plot(pca_corr$variates$X, col = col.lignee, pch = pch.annee, 
     xlab = round((pca_corr$prop_expl_var$X[1])*100), 
     ylab = round((pca_corr$prop_expl_var$X[2])*100))
legend("topleft", col = color.mixo(1:2), legend = levels(lignee),
       lty = 1, title = "Lignee")
legend("bottomright", legend = levels(annee), pch = 1:3,
       title = "Annee")
title("PCA zootechnie corrigée pour année et n° lot (MUFC,GRFC,POFIC,PODEC,GMQ)")

plotIndiv(pca_raw_DAC, group = DAC_nmr$annee)
plotIndiv(pca_corr, group = DAC_nmr$annee)
plotIndiv(pca(X = scaled_DAC %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")), 
          group=DAC_nmr$annee, scale = FALSE)

cov(scaled_DAC_lm %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"))
cov(DAC_nmr %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"))

corrplot(cor(DAC_nmr %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), method = "pearson"))
corrplot(cor(scaled_DAC_lm %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), method = "pearson")
)

DAC_nmr %>% 
  dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ") %>%
  ggpairs(aes(colour=DAC_nmr$lignee))

scaled_DAC_lm %>% 
  dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ") %>%
  ggpairs(aes(colour=DAC_nmr$lignee))

DAC_nmr %>% 
  dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ") %>%
  ggpairs(aes(colour=DAC_nmr$annee))

scaled_DAC_lm %>% 
  dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ") %>%
  ggpairs(aes(colour=DAC_nmr$annee))

###### iii. scaled DAC (selected) ###### 

pdf("Output/DAC/EDA_Zootech/PCA/pca_scaledSelect_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(pca(X = scaled_DAC %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ")), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Pca scaled DAC select", 
                                            scale = FALSE,
                                            ellipse = TRUE))
dev.off()

##### b) PLS-DA #####

ifelse(!dir.exists("Output/DAC/EDA_Zootech/PLSDA"), dir.create("Output/DAC/EDA_Zootech/PLSDA", recursive=FALSE), "folder exists already")


###### i. raw DAC (all) ######

pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_raw_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = quanti %>% dplyr::select(names(scaled_DAC)), 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda raw DAC", 
                                            ellipse = TRUE))
dev.off()

###### ii. scaled DAC (all) ###### 

pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_scaled_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = scaled_DAC, 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda scaled DAC",
                                            scale = FALSE,
                                            ellipse = TRUE))
dev.off()

###### iii. raw DAC (selected) ######

pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_rawSlect_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = quanti %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda raw select DAC", 
                                            ellipse = TRUE))
dev.off()

###### iv. corrected DAC ######
#lm
pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_corrSelect_DAC_lm.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = corr_DAC_lm %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda corr lm select DAC", 
                                            scale = FALSE,
                                            ellipse = TRUE))
dev.off()

pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_corrSelect_DAC_lmrob.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = corr_DAC %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda corr lmrob select DAC", 
                                            scale = FALSE,
                                            ellipse = TRUE))
dev.off()

###### v. scaled DAC (selected) ###### 

pdf("Output/DAC/EDA_Zootech/PLSDA/plsda_scaledSelect_DAC.pdf")
lapply(1:ncol(quali), function(x) plotIndiv(plsda(X = scaled_DAC %>% dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"), 
                                                  Y = quali[[x]]), 
                                            group=quali[[x]], 
                                            legend = TRUE, legend.title = names(quali[x]), 
                                            title = "Plsda scaled select DAC", 
                                            scale = FALSE,
                                            ellipse = TRUE))
dev.off()

###### vi. lignee #######
plsda_lignee_raw = plsda(X = DAC_nmr %>% 
                           dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"),
                         Y = DAC_nmr$lignee,
                         scale = TRUE)

plsda_lignee_corr = plsda(X = corr_DAC %>% 
                            dplyr::select("MUFC","GRFC","POFIC","PODEC", "GMQ"),
                          Y = DAC_nmr$lignee,
                          scale = TRUE)

plotIndiv(plsda_lignee_raw, group = DAC_nmr$lignee, 
          pch = DAC_nmr$annee,
          legend = TRUE, legend.title = 'Lignée',
          legend.title.pch = 'Année',
          title = 'PLSDA données brutes zootechnie')

plotIndiv(plsda_lignee_corr, group = DAC_nmr$lignee, 
          pch = DAC_nmr$annee,
          legend = TRUE, legend.title = 'Lignée',
          legend.title.pch = 'Année',
          title = 'PLSDA données corrigées zootechnie')

