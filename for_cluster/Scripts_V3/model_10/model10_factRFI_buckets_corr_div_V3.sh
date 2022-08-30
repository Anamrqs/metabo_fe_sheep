#!/bin/bash
#SBATCH -J BfactRC_predict
#SBATCH --nodes=1
#SBATCH -c 40
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.stdout
#SBATCH --error=slurm-%j.stderr
#SBATCH --mail-user=anais.marquisseau@inrae.fr
#SBATCH --mail-type=ALL
#SBATCH --mem=32G #by node
#SBATCH --partition=workq

# set up a scratch directory
SCRATCHDIR=/work/$USER/scratchdir/job-$SLURM_JOB_ID
mkdir -m 700 $SCRATCHDIR
cd $SCRATCHDIR

# set up job environnement
module purge
module load system/R-4.1.2_gcc-9.3.0/

# copy template
cp /save/amarquissea/analysis/Scripts_V3/template_V3.R .

# Modif template script
sed -i 's/nbCores/39/g' template_V3.R
sed -i 's/repetsInner/50/g' template_V3.R
sed -i 's/foldsInner/5/g' template_V3.R
sed -i 's/repetsOutter/2/g' template_V3.R
sed -i 's/foldsOutter/5/g' template_V3.R
sed -i 's/PREDICTIONS/factrfi/g' template_V3.R
sed -i 's/PREDICTEURS/buckets/g' template_V3.R
sed -i 's/DATASET/df_corr_rfi_V3/g' template_V3.R
sed -i 's/NCOMP/10/g' template_V3.R
sed -i 's/LIST_KEEPX/c(seq(1,5,1),10,25,50,100)/g' template_V3.R
sed -i 's/MODEL/splsda/g' template_V3.R
sed -i 's/RANGE_X/c(37:length(df_prediction))/g' template_V3.R
sed -i 's/RANGE_Y/6/g' template_V3.R
sed -i 's/CORRECTION/corr_div/g' template_V3.R

# Modif script name
mv template_V3.R model10_factRFI_Buckets_corr_V3.R

# Run Rscript, output a logfile -- prediction ingestion --
Rscript --vanilla --verbose model10_factRFI_Buckets_corr_V3.R > slurm-${SLURM_JOB_ID}.Rout 2>&1 

# append logfile to this script logfile
cat slurm-${SLURM_JOB_ID}.Rout >> slurm-${SLURM_JOB_ID}.out

# remove Rout log
rm slurm-${SLURM_JOBID}.Rout

# Move results in the initial directory
mv model10_factRFI_Buckets_corr_V3.R /work/$USER/Scripts_V3/
mv $SCRATCHDIR ~/work/
