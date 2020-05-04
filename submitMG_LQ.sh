#! /bin/bash
#$ -l h=!(t3wn46|t3wn5*|t3wn6*)


# Settings for SLURM
#
#SBATCH --account=t3
#SBATCH --time 01:00:00
#SBATCH --chdir /work/areimers/workdir_slurm
#SBATCH -e %x-%A.err
#SBATCH -o %x-%A.out
#SBATCH --mail-type FAIL
#SBATCH --mail-user arne.reimers@physik.uzh.ch
##########SBATCH --mem-per-cpu=4gb
##########SBATCH --export ALL

printf "###############################################\n"
printf "##   Run MadGraph for sample %-16s##\n" "$1"
printf "###############################################\n"



SAMPLE="$1"
BASESAMPLE="$SAMPLE" # SLQ, VLQ

DBG=2
LOGDIR="/work/areimers/workdir_slurm"
MGDIR="/work/areimers/MG5_aMC_v2_7_2"
BASEDIR="/work/areimers/GENSIM"
CMSSWDIR="/work/areimers/CMSSW_10_2_10"
XROOTD="root://t3dcachedb.psi.ch:1094"
GFAL="gsiftp://t3se01.psi.ch"
SE_HOME="/pnfs/psi.ch/cms/trivcat/store/user/$USER"
JOBDIR="LQCrossSections"

# optional parameters (always after non-optional ones!)
OPTIND=2 # look for optional arguments after first two non-optional ones
XSEC=0
MASS=1000
LAMBDA=1.0
BWCUTOFF=15
KAPPA=-1.0
NEVENT=10000
EBEAMS="6500.0"
TAG=""
# OPTS="-f "
OPTS=""
echo ">>> arguments are $@"
while getopts E:B:M:L:K:N:T:j:c:l:Pt:x option; do
case "${option}" in
  E) EBEAMS=${OPTARG};;
  B) BWCUTOFF=${OPTARG};;
  M) MASS=${OPTARG};;
  L) LAMBDA=${OPTARG};;
  K) KAPPA=${OPTARG};;
  N) NEVENTS=${OPTARG};;
  T) TAG=${OPTARG};;
  j) OPTS+="--multicore --nb_core=${OPTARG} ";;
  c) OPTS+="--multicore --nb_core=${OPTARG} ";;
  P) OPTS+="-p ";; # Stop the run after the parton level file generation
  x) XSEC=1;;
esac
done
LAMBDASTR=${LAMBDA/./p}
KAPPASTR=${KAPPA/./p}

# ISNLO=false
# if [[ "$SAMPLE" = *"Scalar"*"Pair"* ]] || [[ "$SAMPLE" = *"Scalar"*"Single" ]]
# then
#   ISNLO=true
# fi
# echo $ISNLO

[[ ($SAMPLE = *Scalar*Pair*) || ($SAMPLE = *Scalar*Single*) ]] && ISNLO=true || ISNLO=false
# [[ $SAMPLE = *Scalar*Pair* ]] && ISNLO=true || ISNLO=false
echo $ISNLO
[[ $SAMPLE = *Vector* ]] && ISVECTOR=true || ISVECTOR=false
echo $ISVECTOR

SCALE=$MASS
[[ "$ISVECTOR" = true ]] && SAMPLE="${SAMPLE}_M${MASS}_L${LAMBDASTR}_K${KAPPASTR}" || SAMPLE="${SAMPLE}_M${MASS}_L${LAMBDASTR}"
[[ "$ISNLO" = true ]] && RUNNAME="aMC@NLO -n $SAMPLE" || RUNNAME="$SAMPLE"
[[ $OPTS != *"multicore"* ]] && OPTS+="--nb_core=1 "

cat <<EOF

###########################################
##            JOB PARAMETERS:            ##
###########################################
  \$BASESAMPLE=$BASESAMPLE
  \$SAMPLE=$SAMPLE
  \$SAMPLE=$SAMPLE
  \$BWCUTOFF=$BWCUTOFF
  \$MASS=$MASS
  \$LAMBDA=$LAMBDA
  \$KAPPA=$KAPPA
  \$EBEAMS=$EBEAMS
  \$SCALE=$SCALE
  \$NEVENTS=$NEVENTS
  \$RUNNAME=$RUNNAME
  \$OPTS=$OPTS
EOF

WORKDIR="/scratch/$USER/$JOBDIR/$SAMPLE"  # local workdir for slurm
WORKMGDIR="$WORKDIR/$BASESAMPLE"     # MG directory generated with proc_card
OUTDIR="/work/$USER/LQCrossSections" # where summary.txt files are stored


CARDS="run param"
CARDDIR="$BASEDIR/cards/CrossSections"  # GENSIM dir / cards (contains templates of cards for x-sec computation)
WORKCARDDIR="$WORKMGDIR/Cards"



##### MONITORING/DEBUG INFORMATION #######################################################

DATE_START=`date +%s`
echo "Job started at " `date`
function peval { echo ">>> $@"; eval "$@"; }
cat <<EOF

###########################################
##       QUEUEING SYSTEM SETTINGS:       ##
###########################################
  \$HOME=$HOME
  \$USER=$USER
  \$JOB_ID=$JOB_ID
  \$JOB_NAME=$JOB_NAME
  \$HOSTNAME=$HOSTNAME
  \$TASK_ID=$TASK_ID
  \$QUEUE=$QUEUE
EOF



##### SET ENVIRONMENT ####################################################################

if test -e "$WORKDIR"; then
   echo "ERROR: WORKDIR ($WORKDIR) already exists!" >&2
   peval "ls $WORKDIR" >&2
   #exit 1
fi
peval "mkdir -p $WORKDIR"
if [ ! -d $WORKDIR ]; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi


cat <<EOF

###########################################
##             JOB SETTINGS:             ##
###########################################
  \$MGDIR=$MGDIR
  \$WORKDIR=$WORKDIR
  \$WORKMGDIR=$WORKMGDIR
  \$WORKCARDDIR=$WORKCARDDIR
  \$VO_CMS_SW_DIR=$VO_CMS_SW_DIR
  \$CMSSWDIR=$CMSSWDIR
EOF



##### MAIN FUNCTIONALITY CODE ############################################################

echo " "
echo "###########################################"
echo "##           SETUP ENVIRONMENT           ##"
echo "###########################################"

peval 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MGDIR/HEPTools/hepmc/lib'
peval 'export CXXFLAGS="$CXXFLAGS -I$MGDIR/HEPTools/boost/include"' # silly lhapdf doesn't find boost includes otherwise... seems to be a MG/LHAPDF error
peval 'export PYTHIA8DATA=/work/areimers/MG5_aMC_v2_7_2/HEPTools/pythia8/share/Pythia8/xmldoc'

echo ">>> setting up CMSSW"
peval "source $VO_CMS_SW_DIR/cmsset_default.sh" || exit 1
peval "export SCRAM_ARCH=slc6_amd64_gcc630" || exit 1
peval "cd $CMSSWDIR/src"  || exit 1
eval `scramv1 runtime -sh`
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi



echo " "
echo "###########################################"
echo "##       SETUP FUNCTIONALITY CODE        ##"
echo "###########################################"

# CREATE process dir
peval "cd $WORKDIR"
THISPOSTFIX=""
if [[ $SAMPLE == *Scalar*Single* ]]; then
  ### TODO: Automatize this, give all template-cards lo or nlo postfix, as well as outputs. Then this could be done automatically instead of by hand...
  # THISPOSTFIX="_lo"
  # THISPOSTFIX="_lobw"
  THISPOSTFIX=""
fi
  THISPROCCARD="$CARDDIR/${BASESAMPLE}_proc_card$THISPOSTFIX.dat"

peval "$MGDIR/bin/mg5_aMC $THISPROCCARD"

# COPY cards
for thiscard in $CARDS; do
  peval "cp $CARDDIR/${BASESAMPLE}_${thiscard}_card$THISPOSTFIX.dat $WORKCARDDIR/${thiscard}_card.dat"
done
if [ -e $CARDDIR/${BASESAMPLE}_madspin_card$THISPOSTFIX.dat ]; then
  peval "cp $CARDDIR/${BASESAMPLE}_madspin_card$THISPOSTFIX.dat $WORKCARDDIR/madspin_card.dat"
fi

# REPLACE parameters

PARAM_CARD="$WORKCARDDIR/param_card.dat"
RUN_CARD="$WORKCARDDIR/run_card.dat"
[ $NEVENTS -le 0 ] && echo "ERROR! NEVENTS=$NEVENTS<=0!" >&2 && exit 1
[ $MASS    -le 0 ] && echo "ERROR! MASS=$MASS<=0!" >&2 && exit 1

# NEVENTS
grep -q '$NEVENTS' $RUN_CARD || echo "WARNING! '\$NEVENTS' not found in $RUN_CARD"
echo ">>> replacing \$NEVENTS to $NEVENTS in $RUN_CARD"
sed -i "s/\$NEVENTS/$NEVENTS/g" $RUN_CARD

# EBEAMS
grep -q '$EBEAMS' $RUN_CARD || echo "WARNING! '\$EBEAMS' not found in $RUN_CARD"
echo ">>> replacing \$EBEAMS to $EBEAMS in $RUN_CARD"
sed -i "s/\$EBEAMS/$EBEAMS/g" $RUN_CARD

# BWCUTOFF
grep -q '$BWCUTOFF' $RUN_CARD || echo "WARNING! '\$BWCUTOFF' not found in $RUN_CARD"
echo ">>> replacing \$BWCUTOFF to $BWCUTOFF in $RUN_CARD"
sed -i "s/\$BWCUTOFF/$BWCUTOFF/g" $RUN_CARD

# MASS
grep -q '$MASS' $PARAM_CARD || echo "WARNING! '\$MASS' not found in $PARAM_CARD"
echo ">>> replacing \$MASS to $MASS in $PARAM_CARD"
sed -i "s/\$MASS/$MASS/g" $PARAM_CARD

# LAMBDA
grep -q '$LAMBDA' $PARAM_CARD || echo "WARNING! '\$LAMBDA' not found in $PARAM_CARD"
echo ">>> replacing \$LAMBDA to $LAMBDA in $PARAM_CARD"
sed -i "s/\$LAMBDA/$LAMBDA/g" $PARAM_CARD

# KAPPA
grep -q '$KAPPA' $PARAM_CARD || echo "WARNING! '\$KAPPA' not found in $PARAM_CARD"
echo ">>> replacing \$KAPPA to $KAPPA in $PARAM_CARD"
sed -i "s/\$KAPPA/$KAPPA/g" $PARAM_CARD

# SCALE
grep -q '$SCALE' $RUN_CARD || echo "WARNING! '\$SCALE' not found in $RUN_CARD"
echo ">>> replacing \$SCALE to $SCALE in $RUN_CARD"
sed -i "s/\$SCALE/$SCALE/g" $RUN_CARD
peval "ls"


echo " "
echo "###########################################"
echo "##         MY FUNCTIONALITY CODE         ##"
echo "###########################################"


if [ $XSEC -gt 0 ]; then

  # CALCULATE cross section
  peval "$WORKMGDIR/bin/calculate_xsect $RUNNAME $OPTS"

  # RENAME output
  peval "cd $OUTDIR"
  peval "ls"
  peval "cat *_banner.txt summary.txt > ${SAMPLE}_summary.txt"
  peval "ls"

else

  # GENERATE events
  peval "pwd"
  peval "ls"
  peval "echo $RUNNAME"
  peval "echo $OPTS"
  # peval "$WORKMGDIR/bin/generate_events $RUNNAME $OPTS"
  if [ "$ISNLO" = true ]; then
    echo "order=NLO" >> commands.dat
    echo "fixed_order=OFF" >> commands.dat
    echo "madanalysis=OFF" >> commands.dat
    # if [[ $SAMPLE = *Scalar*Single* ]]; then
    #   echo "madspin=ON" >> commands.dat
    # else
    #  echo "madspin=OFF" >> commands.dat
    echo "madspin=OFF" >> commands.dat
    #fi
  else
    echo "analysis=OFF" >> commands.dat
  fi
  # echo "shower=PYTHIA8" >> commands.dat
  echo "shower=OFF" >> commands.dat
  echo "done" >> commands.dat
  echo "done" >> commands.dat
  # peval "cat commands.dat"
  peval "cat commands.dat | $WORKMGDIR/bin/generate_events $RUNNAME $OPTS"


fi



##### RETRIEVAL OF OUTPUT FILES AND CLEANING UP ##########################################
peval "ls"
# peval "cat ${BASESAMPLE}/SubProcesses/results.dat"
# peval "ls ${BASESAMPLE}/Events/${SAMPLE}"
# peval "cat ${BASESAMPLE}/Events/${SAMPLE}/parton_systematics.log"
peval "mkdir -p $OUTDIR"

if [ "$ISNLO" = true ]
then
  peval "mv ${BASESAMPLE}/Events/${SAMPLE}/summary.txt $OUTDIR/${SAMPLE}${TAG}_output_nlo.txt" # for NLO samples
else
  # peval "mv ${BASESAMPLE}/Events/${SAMPLE}/parton_systematics.log $OUTDIR/${SAMPLE}_output$THISPOSTFIX.txt" # for LO samples
  peval "cp $LOGDIR/$SLURM_JOB_NAME-$SLURM_JOB_ID.out $OUTDIR/${SAMPLE}${TAG}_output$THISPOSTFIX.txt"
fi

printf "\nRemoving $WORKDIR (on ${HOSTNAME})\n"
peval "rm -rf $WORKDIR"



##########################################################################################

DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
printf "\n#####################################################"
printf "\n    Job finished at %s" "$(date)"
printf "\n    Wallclock running time: %02d:%02d:%02d" "$(( $RUNTIME / 3600 ))" "$(( $RUNTIME % 3600 /60 ))" "$(( $RUNTIME % 60 ))"
printf "\n#####################################################\n\n"

exit 0
