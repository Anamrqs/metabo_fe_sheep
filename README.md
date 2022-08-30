# Relationship between feed efficiency and metabolome in sheep
Master's degree internship in INRAE Toulouse-Occitanie, [GenPhySE](https://genphyse.toulouse.inra.fr/), GESPR team.  
Supervisors : Flavie Tortereau, Christel Marie-Etancelin, Quentin Le Graverand 

Part of the [SMARTER](https://www.smarterproject.eu/) project.

## Experiment 
The experiment was performed from 2018 to 2020.  
During their growth phase, lambs from two divergent lines for RFI (Residual Feed Intake) were fed with a 100% concentrate diet.
Their characteristics were recorded among which :
- Body weight
- Average daily feed intake
- Back fat thickness
- Muscle depth

Their average daily gain, feed conversion ratio and phenotypic residual feed intake were calculated. 
Their blood was collected for untargeted metabolomics (NMR). 

## Rproject

The Rproject was run with R v4.1.2.
It contains :
- 3 folders (Data, for_cluster and Output)
- 4 Scripts
- 2 html Reports (and their Rmarkdown version)  
  
html reports available [here](https://filesender.renater.fr/?s=download&token=cee0fbd8-e0ce-4e98-b762-6a643f89b5cd)

### Data (Private folder) 

- The "corresp_RMN_nuofov.xlsx" file contains correspondences between animals's ID and their metabolomic samples. 

#### Zootechnical

- "variables_dac_2018_2019_2020.txt" is an explicative file about all the available variables written by Quentin Le Graverand. 
- The "zootechnie_DAC.csv" file contains all the informations and performances of the sheeps (variables in columns, animals in rows). 

#### Metabolome 

- The "RMN" folder contains NMR buckets, acquired by Flavie Tortereau with the ASICS R package.

columns = Animals, rows = buckets (chemical shifts), buckets's ID = first column, animals's ID = first row
files format : txt, separator : space

- The "quantifications" folder contains metabolite quantifications of the experiment, acquired by Florian Touitou with the ASICS R package. 

columns = Metabolites's quantifications (and some informations about the animal's line, year...), rows = animals 
files format : xlsx

- The "Metabolites_pics_florian.tsv" file corresponds to peaks manually observed and annotated by Florian Touitou, more precisions in the "metabolites_pics.pptx" file.

### Scripts

- 0_Cleaning.R
- 1a_EDA_Zootech_DAC.R
- 1b_EDA_Buckets_DAC.R
- t_tests.R (private file)

The cleaning script allows to clean zootechnical and metabolomic data. Among the cleaning steps : searching for outliers, zero imputation of the buckets and CLR transformation, and write tables for the prediction scripts.
The EDA scripts permit the exploration of both zootechnical and metabolomics datasets, and their outputs are recorded in the Output folder. 
The t_tests.R script was developped by Quentin Le Graverand to compare prediction results. Comparisons are made using Student's t-tests corrected with Nadeau and Bengio's method, which takes into account the dependence between samples during repeated cross-validations.

### for_cluster 

Folder with scripts and outputs obtained using the high-performance computing resources of the [GenoToul bioinformatics facility](http://bioinfo.genotoul.fr/). 

It contains 3 folders and 1 R script.
- "V3" folder : contains tables from "0_Cleaning.R" script. 
- "Scripts_V3"
- "output_V3"
- "template_V3.R"

### Output

Code outputs


