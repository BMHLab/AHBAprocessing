#!/bin/bash

#SBATCH --account=monash076
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=aurina.arnatkeviciute@monash.edu
#SBATCH --mail-type=ALL

# massive usage (more than one subject)
#  1. [code]$ SUBJECTIDS="1"; SUBJECT_CODE=H0351_2001
#  2. [....]$ for SUBJECTID in $SUBJECTIDS; do sbatch --job-name="${SUBJECTID}" --output="/gpfs/M2Scratch/Monash076/aurina/AllenInstitute/slurm-${SUBJECTID}_surfaces.out" --error="/gpfs/M2Scratch/Monash076/aurina/AllenInstitute/slurm-${SUBJECTID}_surfaces.err" runSurf2surf.sh $SUBJECTID $SUBJECT_CODE; done

SUBJECTID=$1
SUBJECT_CODE=$2

module load freesurfer
cd code
source SetupEnv.sh
cd ..
WHEREISMYSCRIPT="/home/aurina/Monash076_scratch/aurina/AllenInstitute"
WHEREISMYDATA="/home/aurina/Monash076_scratch/aurina/HumanExpression/data/genes/forFreesurfer/S0${SUBJECTID}"
SAMPLELIST=${WHEREISMYDATA}/S${SUBJECTID}sampleList.txt; SAMPLEIDS=`cat ${SAMPLELIST}`;
for SAMPLEID in $SAMPLEIDS; do 
mri_surf2surf --srcsubject ${SUBJECT_CODE} --trgsubject fsaverage --hemi lh --sval ${WHEREISMYDATA}/S${SUBJECTID}sample${SAMPLEID}_singleVert.mgz --nsmooth-in 1 --tval ${WHEREISMYDATA}/S${SUBJECTID}sample${SAMPLEID}ONfsaveragePial.mgz
done
