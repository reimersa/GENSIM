#! /usr/bin/env python
# Author: Izaak Neutelings (October, 2019)
from CRABClient.UserUtilities import getUsernameFromSiteDB
from CRABClient.UserUtilities import config as crabconfig
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from utils import filterSamplesWithPattern, getSampleSites, formatTag, bold
from eras import globaltags
import re


def main(args):
  
  years    = args.years
  pset     = args.pset
  samples  = args.samples
  vetoes   = args.vetoes
  priority = args.priority
  force    = args.force
  test     = args.test
  tag      = "DeepTau2017v2p1"
  
  
  # SAMPLES
  if 'nanoaod' in pset.lower():
    import samples_nanoAOD
    samplesets = samples_nanoAOD.samples
  else:
    import samples_miniAOD
    samplesets = samples_miniAOD.samples
  
  # SUBMIT
  for year in years:
    datasets = samplesets.get(year,[])
    
    if samples:
      datasets = filterSamplesWithPattern(datasets,samples)
    if vetoes:
      datasets = filterSamplesWithPattern(datasets,vetoes,veto=True)
    submitSampleToCRAB(pset,year,datasets,tag=tag,priority=priority,
                       test=test,force=force)
  


def submitSampleToCRAB(pset,year,samples,**kwargs):
  """Create a CRAB configuration and submit a given list of samples."""
  
  assert isinstance(samples,list), "Samples list should be a list or tuple! Given %s"%samples
  
  # USER OPTIONS
  year          = year
  test          = kwargs.get('test',       0)
  force         = kwargs.get('force',      False)
  datatier      = 'nanoAOD' if 'nanoaod' in pset.lower() else 'miniAOD'
  version       = re.findall("(?<=AOD)v\d+",pset)
  version       = version[0] if version else ""
  pluginName    = 'Analysis' #'PrivateMC'
  splitting     = kwargs.get('split',      'FileBased' if year==2018 or datatier=='nanoAOD' else 'Automatic')
  tag           = kwargs.get('tag',        "")
  instance      = kwargs.get('instance',   'global')
  nevents       = -1
  unitsPerJob   = 1 # files per job for 'FileBased'
  njobs         = -1
  ncores        = kwargs.get('ncores',     1) # make sure nCores > nThreads in pset.py
  maxRunTime    = kwargs.get('maxRunTime', 6*60) #1250 # minutes
  priority      = kwargs.get('priority',   10)
  workArea      = "crab_tasks" #"crab_projects"
  outdir        = '/store/user/%s/%s_%s%s'%(getUsernameFromSiteDB(),datatier,year,formatTag(tag))
  publish       = True #and False
  site          = 'T2_CH_CSCS'
  
  # OVERRIDE
  if test>0:
    splitting    = 'FileBased'
    unitsPerJob  = 1 # files per job
    njobs        = int(test)
    outdir      += '_test'
    publish      = False
  if splitting=='Automatic':
    unitsPerJob  = -1
    njobs        = -1
    maxRunTime   = -1
  
  # PRINT
  print ">>> "+'='*70
  print ">>> year        = %s"%year
  print ">>> pset        = '%s'"%bold(pset)
  print ">>> pluginName  = '%s'"%pluginName
  print ">>> splitting   = '%s'"%splitting
  print ">>> tag         = '%s'"%bold(tag)
  print ">>> nevents     = %s"%nevents
  print ">>> unitsPerJob = %s"%unitsPerJob
  print ">>> njobs       = %s"%njobs
  print ">>> nCores      = %s"%ncores
  print ">>> maxRunTime  = %s"%maxRunTime
  print ">>> priority    = %s"%priority
  print ">>> workArea    = '%s'"%workArea
  print ">>> site        = '%s'"%site
  print ">>> outdir      = '%s'"%outdir
  print ">>> publish     = %r"%publish
  print ">>> "+'='*70
  
  if len(samples)==0:
    print ">>> No samples given..."
    print ">>> "
    return
  
  # CRAB CONFIGURATION
  config = crabconfig()
  config.General.workArea           = workArea
  config.General.transferOutputs    = True
  config.General.transferLogs       = False
  
  config.JobType.pluginName         = pluginName
  config.JobType.psetName           = pset
  config.JobType.pyCfgParams        = [ "year=%s"%year, "nThreads=%s"%ncores ]
  config.JobType.numCores           = ncores
  if maxRunTime>0:
    config.JobType.maxJobRuntimeMin = maxRunTime # minutes
  config.JobType.priority           = priority
  
  config.Data.splitting             = splitting
  if unitsPerJob>0:
    config.Data.unitsPerJob         = unitsPerJob
    if njobs>0:
      config.Data.totalUnits        = unitsPerJob * njobs
  config.Site.storageSite           = site
  config.Data.outLFNDirBase         = outdir
  config.Data.publication           = publish
  
  for dataset in samples:
    
    # INDIVIDUAL CONFIG
    request       = (datatier.lower().replace('aod','')+'_'+shortenDASPath(dataset))[:100]
    private       = dataset.endswith('/USER')
    sites         = getSampleSites(dataset,instance=None)
    if private:
      ignoreLocal = True
      inputDBS    = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/"
      whitelist   = getOptimalWhitelist(sites,instance=instance)
      #whitelist   = ['T2_CH_*','T2_DE_*','T2_IT_*']
    else:
      ignoreLocal = False
      inputDBS    = "https://cmsweb.cern.ch/dbs/prod/%s/DBSReader/"%instance
      whitelist   = [ ]
    outtag        = createDatasetOutTag(dataset,tag=tag,datatier=datatier,version=version,year=year)
    
    # PRINT
    print ">>> "+'-'*5+" Submitting... "+'-'*50
    print ">>> request     = '%s'"%bold(request)
    print ">>> dataset     = '%s'"%bold(dataset)
    print ">>> inputDBS    = '%s'"%inputDBS
    print ">>> sites       = %s"%sites
    print ">>> whitelist   = %s"%whitelist
    print ">>> ignoreLocal = %s"%ignoreLocal
    print ">>> outtag      = '%s'"%outtag
    print ">>> "+'-'*70
    
    # INDIVIDUAL CONFIG
    config.General.requestName    = request # max. 100 characters
    config.Data.inputDataset      = dataset
    config.Data.inputDBS          = inputDBS
    #config.Data.outputPrimaryDataset = 'LQ_test' # only for 'PrivateMC'
    config.Data.outputDatasetTag  = outtag
    config.Data.ignoreLocality    = ignoreLocal # do not run on same site the dataset is stored on
    if whitelist:
      config.Site.whitelist       = whitelist
    print str(config).rstrip('\n')
    print ">>> "+'-'*70
    
    # SUBMIT
    if force:
      print ">>> Do you want to submit this job to CRAB? [y/n]? force"
      print ">>> Submitting..."
      submitCRABConfig(config)
    else:
      while True:
        submit = raw_input(">>> Do you want to submit this job to CRAB? [y/n]? ").strip().lower()
        if any(s in submit for s in ['quit','exit']):
          print ">>> Exiting..."
          exit(0)
        elif 'force' in submit:
          submit = 'y'
          force = True
        if 'y' in submit:
          print ">>> Submitting..."
          submitCRABConfig(config)
          break
        elif 'n' in submit:
          print ">>> Not submitting."
          break
        else:
          print ">>> '%s' is not a valid answer, please choose 'y' or 'n'."%submit
    
    print ">>> "
  


def executeCRABCommand(command,**kwargs):
  """Submit a command and catch exceptions."""
  try:
    crabCommand(command,**kwargs)
  except HTTPException as hte:
    print "Failed submitting task: %s"%(hte.headers)
  except ClientException as cle:
    print "Failed submitting task: %s"%(cle)
  
def submitCRABConfig(config):
  """Submit a single config."""
  executeCRABCommand('submit',config=config)
  
countrypattern = re.compile(r"T2_([A-Z]{2})_.+")
def getOptimalWhitelistForDataset(dataset,instance='global'):
  """Get optimal whitelist for a DAS dataset."""
  # https://dashb-ssb.cern.ch/dashboard/request.py/siteviewhome
  sites = getSampleSites(dataset,instance=None)
  return getOptimalWhitelist(sites)
  
def getOptimalWhitelist(sites,instance='global'):
  """Get optimal whitelist for a list of sites."""
  countries = [ ]
  for site in sites:
    countries += countrypattern.findall(site)
  whitelist = [ ]
  if any(c in ['US','CA','MX'] for c in countries): # America
    whitelist += ['T2_US_*','T2_CA_*','T2_MX_*']
  if any(c in ['CH','DE','IT','FR','BE','UK','FI','ES','PT'] for c in countries): # West-Europe
    whitelist += ['T2_CH_*','T2_DE_*','T2_IT_*','T2_FR_*','T2_BE_*','T2_UK_*','T2_FI_*','T2_ES_*','T2_PT_*']
  if any(c in ['EE','HR','GR','PL','RU','UA','HU','BY','BG'] for c in countries): # East-Europe
    whitelist += ['T2_EE_*','T2_HR_*','T2_GR_*','T2_PL_*','T2_RU_*','T2_UA_*','T2_HU_*','T2_BY_*','T2_BG_*']
  if any(c in ['TR','IR','IN'] for c in countries): # Middle-East
    whitelist += ['T2_TR_*','T2_IR_*','T2_IN_*']
  if any(c in ['TW','CN','KR','TH'] for c in countries): # Asia
    whitelist += ['T2_TW_*','T2_CN_*','T2_KR_*','T2_IN_*','T2_TH_*']
  if not whitelist:
    whitelist = sites or ['T2_*']
  return whitelist
  
  
hashpattern = re.compile(r"-(?!v)[0-9a-z]{5,}$")
minipattern = re.compile(r"[Mm]iniAOD(?:v\d+)?")
def shortenDASPath(daspath):
  """Shorten a long DAS path."""
  replace = [
    ('Aut18',          "RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15"),
    ('Aut18FlatPU',    "RunIIAutumn18MiniAOD-FlatPU0to70_102X_upgrade2018_realistic_v15"),
    ('Aut18FlatPURAW', "RunIIAutumn18MiniAOD-FlatPU0to70RAW_102X_upgrade2018_realistic_v15"),
    ('Fall17',         "RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ('Fall17newPMX',   "RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14"),
    ('Fall17RECOSIM',  "RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14"),
    ('Fall17v2',       "RunIIFall17MiniAODv2-PU2017_12Apr2018_v2_94X_mc2017_realistic_v14"),
    ('Sum16',          "RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3"),
    ("MG",             "madgraph"),
    ("PY",             "pythia"),
    ("",               "_13TeV"),
    ("",               "Tune"),
    ("Incl",           "Inclusive"),
  ]
  daspath = '_'.join(daspath.split('/')[1:-1])
  for short, long in replace:
    daspath = daspath.replace(long,short)
  daspath = hashpattern.sub("",daspath)
  daspath = daspath.replace(getUsernameFromSiteDB()+'-',"")
  return daspath
  
def getCampaign(daspath):
  """Return a simple campaign label."""
  for campaign in [ 'Summer16', 'Fall17', 'Autumn18' ]:
    if campaign in daspath:
      return campaign
  return daspath
  
def createDatasetOutTag(dataset,tag='',datatier=None,year=None,version=""):
  """Create dataset output tag from DAS path."""
  outtags = dataset.strip('/').split('/')
  assert len(outtags)>=2, "Invalid DAS path '%s'!"%(dataset)
  outtag   = outtags[1]
  outtag   = outtag.replace(getUsernameFromSiteDB()+'-',"")
  outtag   = hashpattern.sub("",outtag)
  dtype    = 'data' if '/Run201' in dataset else 'mc'
  if datatier=='nanoAOD':
    if 'miniaod' in outtag.lower():
      newtier = 'NanoAOD'+version
      outtag  = minipattern.sub(newtier,outtag)
    globaltag = globaltags[dtype]['nanoAOD'].get(year,False)
    if globaltag:
      outtag = outtag[:outtag.index(newtier)+len(newtier)]+'_'+globaltag
  if tag and not outtag.endswith(tag):
    outtag += formatTag(tag)
  return outtag
  
def createRequest(string,year=None,tag=""):
  """Create request."""
  request = string.split('/')[-1].replace('pset_','').replace('.py','')
  if year and str(year) not in request:
    request += "_%s"%(year)
  return request
  


if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('pset',             type=str, action='store',
                      metavar='FILE',     help="parameter-set configuration file to submit" )
  parser.add_argument('-f', '--force',    dest='force', action='store_true', default=False,
                                          help="submit jobs without asking confirmation" )
  parser.add_argument('-y', '--year',     dest='years', choices=[2016,2017,2018], type=int, nargs='+', default=[2017], action='store',
                      metavar='YEAR',     help="select year" )
  parser.add_argument('-s', '--sample',   dest='samples', type=str, nargs='+', default=[ ], action='store',
                      metavar='DATASET',  help="samples to submit" )
  parser.add_argument('-x', '--veto',     dest='vetoes', nargs='+', default=[ ], action='store',
                      metavar="DATASET",  help="exclude/veto this sample" )
  parser.add_argument('-p', '--priority', dest='priority', type=int, default=10, action='store',
                      metavar='PRIORITY', help="submit with priority (default=10)" )
  parser.add_argument('-t', '--test',     dest='test', type=int, nargs='?', default=-1, const=1,
                      metavar='NJOBS',    help="submit test job(s)" )
  args = parser.parse_args()
  print ">>> "
  main(args)
  #print ">>> "
  

