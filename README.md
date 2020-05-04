# CRAB

Generate events or process `miniAOD` files with `CRAB3`.

#### Table of Contents  
* [Installation](#Installation)<br>
  * [DeepTau2017v2p1](#DeepTau2017v2p1)
  * [Environment setup](#environment-setup)
* [Event processing](#event-processing)<br>
  * [Local run](#local-run)<br>
  * [CRAB submission](#CRAB-submission)<br>
* [Event generation](#event-generation)<br>
  * [MadGraph gridpack generation](#madgraph-gridpack-generation)
  * [Local event generation](#local-event-generation)
* [Notes](#Notes)<br>
  * [CRAB3](#CRAB3)
  * [NanoAOD](#NanoAOD)
  * [Samples](#Samples)
  * [PSet fragments](#PSet-fragments)


## Installation

First, get this repository:
```
git clone https://github.com/IzaakWN/CRAB CRAB
cd CRAB
```
then install the relevant `CMSSW` in this directory, e.g.:
```
CMSSW=CMSSW_10_2_16_patch1
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel $CMSSW
cd $CMSSW/src
cmsenv
scram b -j4
cd ../..
```
When you want to use CRAB, you need to source the relevant setup script:
```
source $VO_CMS_SW_DIR/crab3/crab_slc6.sh  # for SLC6
source $VO_CMS_SW_DIR/crab3/crab.sh       # for SLC7
```

Make sure you have a valid VOMS proxy,
```
voms-proxy-info --timeleft              # check how many seconds you have left
voms-proxy-init -voms cms -valid 400:0  # renew proxy if it is too short
source setupVOMS.sh                     # OR, source this script instead of the above two lines
```

Before submitting jobs to CRAB, make sure that in [`submit_crab.py`](submit_crab.py) you have specified the correct tier storage element that you have writing permissions to. For example with `'T2_CH_CSCS'` for PSI's T2:
```
sed "s/\(site\s*=\s*\)'\w+'/\1'T2_CH_CSCS'/" submit_crab.py
``` 


### DeepTau2017v2p1

In case you want to use `DeepTau2017v2p1` in `102X` samples, see this [this TWiki page](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Running_of_the_DeepTauIDs_ver_20), and setup `CMSSW` as
```
CMSSW=CMSSW_10_2_16_patch1
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel $CMSSW
cd $CMSSW/src/
cmsenv
git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2
git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_deepTauVetoPCA
sed 's/idDeepTau2017v2/idDeepTau2017v2p1/g' PhysicsTools/NanoAOD/python/taus_cff.py -i
sed 's/rawDeepTau2017v2/rawDeepTau2017v2p1/g' PhysicsTools/NanoAOD/python/taus_cff.py -i
scram b -j 4
cd ../..
```


### Environment setup

With each new shell session, do something like
```
cd CRAB
source setupVOMS.sh
source $VO_CMS_SW_DIR/cmsset_default.sh
source $VO_CMS_SW_DIR/crab3/crab_slc6.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cd CMSSW_10_2_16_patch1/src
cmsenv
cd ../..
```


## Event processing

### Local run

Get some files with e.g.
```
dasgoclient --limit=0 --query="dataset=/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM file" | head -n1
xrdcp -f root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/042C8EE9-9431-5443-88C8-77F1D910B3A5.root ./test.root
```

Reprocess `miniAOD` to add `DeepTau2017v2p1` with
```
cmsRun pset_miniAOD_rerun.py
```

Process `miniAOD` to `nanoAODv5` with
```
cmsRun pset_nanoAODv5.py
```


### Check output

Check the content of a `miniAOD` file with
```
edmProvDump input/VLQ-p_M1100_2017_rerun.root > dump.log
```
or use as an example (see [this TWiki page](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Example_code_accessing_all_high))
```
./checkMiniAOD.py
```



### CRAB submission

Collect your favorite samples in the respective python file, sorted in a dictionary per year:
```
samples_miniAOD.py
samples_nanoAOD.py
```

Then you can submit a specific parameter-set configuration file as
```
./submit_crab.py pset_nanoAODv5.py -y 2017
```

To filter a specific sample, do e.g.
```
./submit_crab.py -y 2017 -s DY*Jets
```

To run test job(s), do e.g.
```
./submit_crab.py -y 2017 -t    # one test job
./submit_crab.py -y 2017 -t 2  # two test jobs
```

Check the task status, and resubmit if needed as:
```
crab status -d crab_projects/crab_<request>
crab resubmit -d crab_projects/crab_<request>
```
replacing `<request>`. Or, use the local [`crab.sh`](crab.sh) script to check many directories at once
```
crab.sh status crab_projects/crab_*
crab.sh resubmit crab_projects/crab_*
```

If you publish your output files (default, if it is not a test job), you can retrieve them in DAS with (replacing `<user>`)
```
dasgoclient --limit=0 --query="dataset=/*/<user>*/USER instance=prod/phys03"
```



## Event generation

### MadGraph gridpack generation
Before you can generate MadGraph events with CMSSW, it's useful to create a gridpack, as described [here](https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuickGuideMadGraph5aMCatNLO), or install
```
git clone git@github.com:cms-sw/genproductions.git genproductions
```
and in a clean shell session, run
```
source $VO_CMS_SW_DIR/cmsset_default.sh
cd genproductions/bin/MadGraph5_aMCatNLO/
./gridpack_generation.sh <process name without _proc_card.dat> <card dir>
```
To produce a large set of MadGraph cards with varying parameters (mass, coupling strenghts, ...), you can prepare some template cards like the examples in [`cards/ScalarLQ_Single/ScalarLQ_Single_template_*.dat`](cards/ScalarLQ_Single), and use `create_cards.py` as e.g.
```
./create_cards.py cards/ScalarLQ_Single/ScalarLQ_Single_template_*.dat -m 1000 -p LAMBDA=0.1,1.0,1.5
```
Produce the cards and gridpacks in series with [`generate_gridpacks.py`](generate_gridpacks.py), e.g.
```
./generate_gridpacks.py cards/ScalarLQ_Single ScalarLQ_Single -m 1000 -p LAMBDA=0.1,1.0,1.5
```

### Local event generation
Produce events with a gridpack with [`pset_GENSIM.py`](pset_GENSIM.py) as e.g.
```
cmsRun pset_GENSIM.py nevents=100 gridpack=ScalarLQ_Single_M500_slc6_amd64_gcc630_CMSSW_9_3_16_tarball.tar.xz
```



## Notes

### CRAB3

* Tutorial: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial
* Configuration: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
* Commands: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands
* FAQ: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ
* Dashboard: https://monit-grafana.cern.ch (Jobs > CMS Tasks Monitoring GlobalView > Select User > `<username>`)


### NanoAOD

* **working book**: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
* **2016 `9_4_X`**: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X2016_doc.html
* **2017 `9_4_X`**: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html
* **2018 `10_2_X`**: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html

More [notes](https://www.evernote.com/l/Ac8PKYGpaJxJArj4eng5ed95_wvpzwSNTgc).


### Samples

[PPD Run II summary table](https://docs.google.com/presentation/d/1YTANRT_ZeL5VubnFq7lNGHKsiD7D3sDiOPNgXUYVI0I/edit#slide=id.g4dfd66f53d_1_7)
* **2016**: [list](samples_2016.cfg), [DAS](https://cmsweb.cern.ch/das/request?view=plain&limit=50&instance=prod%2Fglobal&input=dataset%3D%2F*%2FRunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic*%2FNANOAODSIM), [notes](https://www.evernote.com/l/Ac9nVeF2tcdJI7R-is1KPT2Ukv7A260zNX0)
* **2017**: [list](samples_2017.cfg), [DAS](https://cmsweb.cern.ch/das/request?view=plain&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*94X*_realistic_v14*%2FNANOAOD*), [notes](https://www.evernote.com/l/Ac8WfL3Mzx1MrKdm1LfIOl-F-j7NeScPKxs)
* **2018**: [list](samples_2018.cfg), [DAS](https://cmsweb.cern.ch/das/request?view=plain&limit=50&instance=prod%2Fglobal&input=%2F*%2FRunIIAutumn18NanoAODv4-Nano14Dec2018*%2FNANOAODSIM), [notes](https://www.evernote.com/l/Ac9yyi7wtg9LaYgxOIz11jFyzLV0ztkemtE)


### PSet fragments

To get a fragment of a `PSet` configuration file with `cmsDriver.py` is explained [in the NanoAOD Workbook](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Running_on_various_datasets_from).
```
cmsDriver.py myNanoProdMc -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions MyGlobalTag --era MyEraModifiers
cmsDriver.py myNanoProdData -s NANO --data --eventcontent NANOAOD --datatier NANOAOD --no_exec --conditions MyGlobalTag --era MyEraModifiers
```

Get a fragment of a `PSet` configuration file for a particular sample:
1. Go to [MCM](https://cms-pdmv.cern.ch/mcm/).
2. Go to 'Request' > 'Output Dataset', and type in the DAS path of your favorite sample, e.g.
```
/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
```
3. In the 'Actions' column, click on the 'Get setup command' symbol (circle with down arrow).<br>
   https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_setup/MUO-RunIIAutumn18MiniAOD-00016<br>
   https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/MUO-RunIIAutumn18MiniAOD-00016
4. Copy-paste the `cmsDriver.py` the command line (in a `CMSSW` environment).
5. You will find a configuration file called `<prep-id>.py`.
