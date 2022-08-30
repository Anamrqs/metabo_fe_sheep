#!/bin/bash
#SBATCH -J BLSC_predict
#SBATCH --nodes=1
#SBATCH -c 40
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%j.stdout
#SBATCH --error=slurm-%j.stderr
#SBATCH --mail-user=anais.marquisseau@inrae.fr
#SBATCH --mail-type=ALL
#SBATCH --mem=64G #by node
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
sed -i 's/repetsOutter/50/g' template_V3.R
sed -i 's/foldsOutter/5/g' template_V3.R
sed -i 's/PREDICTIONS/lignee/g' template_V3.R
sed -i 's/PREDICTEURS/zootech/g' template_V3.R
sed -i 's/DATASET/df_SC_V3/g' template_V3.R
sed -i 's/NCOMP/5/g' template_V3.R
sed -i 's/LIST_KEEPX/c(seq(1,5,1),10,22)/g' template_V3.R
sed -i 's/MODEL/splsda/g' template_V3.R
sed -i 's/RANGE_X/c(13:36)/g' template_V3.R
sed -i 's/RANGE_Y/5/g' template_V3.R
sed -i 's/CORRECTION/SC/g' template_V3.R

# Modif script name
mv template_V3.R model250_Lignee_Zootech_SC_V3.R

# Run Rscript, output a logfile -- prediction ingestion --
Rscript --vanilla --verbose model250_Lignee_Zootech_SC_V3.R > slurm-${SLURM_JOB_ID}.Rout 2>&1 

# append logfile to this script logfile
cat slurm-${SLURM_JOB_ID}.Rout >> slurm-${SLURM_JOB_ID}.out

# remove Rout log
rm slurm-${SLURM_JOBID}.Rout

# Move results in the initial directory
mv model250_Lignee_Zootech_SC_V3.R/work/$USER/Scripts_V3/
mv $SCRATCHDIR ~/work/
