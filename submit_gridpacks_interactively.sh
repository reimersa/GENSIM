#!/bin/bash
#

ulimit -a

echo HOME: $HOME
echo USER: $USER
echo HOSTNAME: $HOSTNAME
pwd


source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
# cd /work/areimers/CMSSW_10_2_10/src/genproductions/bin/MadGraph5_aMCatNLO

TASKCMD="./gridpack_generation.sh ScalarLQ_Single_M1500_L2p5_LOPDF31_FIRSTGEN ../../../../../GENSIM/cards/ScalarLQ_Single local"
cd /work/areimers/CMSSW_10_2_10/src/genproductions/bin/MadGraph5_aMCatNLO

echo $TASKCMD
eval $TASKCMD


# cleaning of temporal working dir when job was completed:
echo Done.
