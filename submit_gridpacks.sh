#!/bin/bash
#
#SBATCH --account=t3
#SBATCH --time 01:00:00
#SBATCH --chdir /work/areimers/workdir_slurm
#SBATCH -J gridpack_array
#SBATCH -e %x-%A-%a.err
#SBATCH -o %x-%A-%a.out
#SBATCH --mem=4000
#SBATCH --mail-type FAIL
#SBATCH --mail-user arne.reimers@physik.uzh.ch
#####SBATCH --export ALL

ulimit -a

echo HOME: $HOME
echo USER: $USER
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME
pwd

# each worker node has local /scratch space to be used during job run
export TMPDIR=/scratch/$USER/test_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
echo TMPDIR: $TMPDIR
mkdir -p $TMPDIR

# actual job
# source /t3home/areimers/genprod.sh
# source /t3home/areimers/gridpacks.sh
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
# cd /work/areimers/CMSSW_10_2_10/src/genproductions/bin/MadGraph5_aMCatNLO


cd /work/areimers/GENSIM
export TASKID=$SLURM_ARRAY_TASK_ID
echo $TASKID
JOBLIST=$1
TASKCMD=$(cat $JOBLIST | sed "${TASKID}q;d")
cd /work/areimers/CMSSW_10_2_10/src/genproductions/bin/MadGraph5_aMCatNLO

echo $TASKCMD
eval $TASKCMD


# cleaning of temporal working dir when job was completed:
rm -rf $TMPDIR
echo Removed TMPDIR: $TMPDIR
echo Done.
sstat -j $SLURM_JOB_ID.batch --format=jobid,MaxRSS,MaxVMSize
