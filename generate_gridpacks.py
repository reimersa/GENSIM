#! /usr/bin/env python
# Author: Izaak Neutelings (November, 2019)
# Source:
#   git clone git@github.com:cms-sw/genproductions.git genproductions
#   cd genproductions/bin/MadGraph5_aMCatNLO/
#   ./gridpack_generation.sh <process name without _proc_card.dat> <card dir>
# Examples
#   ./generate_gridpacks.py cards/ScalarLQ_Pair ScalarLQ_Pair -m 800 1100 1500 -p 'LAMBDA=1.5,2.0,2.5'
#   ./generate_gridpacks.py cards/VectorLQ_Single VectorLQ_Single -m 1000 1500 2000 -p 'LAMBDA=1.5,2.0,2.5:KAPPA=0.0,1.0'
import os, re, glob, shutil
import itertools
import subprocess
from utils import formatTag, bold, error, warning, ensureDirectory
from create_cards import makeCardLabel, makeCardName, makeCard, getCards



def main(args):
  
  carddir    = args.carddir
  sample     = args.sample
  years      = args.years
  cardlabel  = args.cardlabel
  masses     = args.masses
  params     = args.params
  tag        = args.tag
  copy       = args.copy
  remove     = not args.keep
  test       = args.test
  cmsswdir   = "/work/areimers/CMSSW_10_2_10"
  genproddir = "genproductions/bin/MadGraph5_aMCatNLO"
  workdir    = "%s/src/%s"%(cmsswdir,genproddir)
  assert os.path.isdir(cmsswdir), error("CMSSW directory '%s' does not exists!"%(cmsswdir))
  assert os.path.isdir(workdir), error("Working directory '%s' does not exists!"%(workdir))
  oldcarddir = carddir[:]
  
  # CHECK environment
  assert not os.getenv('CMSSW_BASE'), error("A CMSSW environment is set. Please retry in a clean shell session.")
  assert os.getenv('CMS_PATH'), error("No CMS default environment set! Please do 'source $VO_CMS_SW_DIR/cmsset_default.sh' first.")
  
  # CREATE POINTS
  if masses:
    keys   = ['MASS']
    params = [('MASS',masses)] #{ 'MASS': masses }
  else:
    keys   = [ ]
    params = [ ]
  if args.params:
    for param in args.params.split(':'):
      assert '=' in param, "Invalid format '%s'; no '=' for '%s'"%(args.params,param)
      param, values = param[:param.index('=')], param[param.index('=')+1:].split(',')
      assert param not in keys, error("Key '%s' defined multiple times!"%param)
      keys.append(param)
      params.append((param,values))
      #params[param] = values
  if not cardlabel:
    cardlabel = '_'.join(k[0]+'$'+k for k in keys)
  if 'OUTPUT' not in keys:
    keys.append('OUTPUT')
    params.append(('OUTPUT',["$SAMPLE_%s"%cardlabel]))
  if params:
    points    = list(itertools.product(*[v for k, v in params]))
    pattern   = os.path.join(carddir,"%s_template*.dat"%(sample))
    templates = glob.glob(pattern)
    assert templates, error("Did not find any template cards '%s' in %s!"%(os.path.basename(pattern),carddir))
  else:
    points    = [ ]
    templates = [ ]
  
  # PRINT
  print ">>> "+'='*90
  print ">>> cmsswdir    = '%s'"%cmsswdir
  print ">>> genproddir  = '%s'"%genproddir
  print ">>> workdir     = '%s'"%workdir
  print ">>> sample      = '%s'"%sample
  print ">>> carddir     = '%s'"%bold(carddir)
  print ">>> params      = %s"%', '.join("%s: %s"%(bold(k),l) for k, l in params)
  print ">>> templates   = '%s'"%"', '".join(bold(os.path.basename(t)) for t in templates)
  print ">>> "+'='*90
  
  # GENERATE
  if points:
    samplenames = [ ]
    for values in points:
      kwargs = { }
      for key, value in zip(keys,values):
        kwargs[key] = value
      for template in templates:
        cardname = makeCardName(template,cardlabel,**kwargs)
        makeCard(template,cardname,**kwargs)
      samplenames.append("%s_%s"%(sample,makeCardLabel(cardlabel,**kwargs)))
    carddir    = os.path.relpath(carddir,workdir)
    os.chdir(workdir)
    for samplename in samplenames:
      generateGridpack(carddir,samplename,remove=remove,copy=copy)
  else:
    carddir    = os.path.relpath(carddir,workdir)
    os.chdir(workdir)
    generateGridpack(carddir,sample,remove=remove,copy=copy)
  

def generateGridpack(carddir, sample, **kwargs):
  """Create a CRAB configuration and submit a given list of samples."""
  
  cards  = getCards(carddir,sample) #[os.path.basename(f) for f in glob.glob("%s/%s_*.dat"%(carddir,sample))]
  copy   = kwargs.get('copy',   False )
  remove = kwargs.get('remove', True  )
  
  # COPY
  if copy:
    newdir = ensureDirectory("%s_InputCards"%sample)
    print ">>> copying cards to '%s'..."%newdir
    for card in cards:
      shutil.copy(os.path.join(carddir,card),newdir)
    carddir = newdir
  
  # PRINT
  print ">>> "
  print ">>> "+'-'*100
  #print ">>> year        = %s"%year
  print ">>> sample      = '%s'"%bold(sample)
  print ">>> cards       = '%s'"%"', '".join(c.replace(sample,'*') for c in cards)
  print ">>> carddir     = '%s'"%bold(carddir)
  
  # CLEAN
  if os.path.join(sample):
    print ">>> "+warning("Directory '%s' already exists! Removing..."%sample)
    rmcommand  = "rm -rf %s"%(sample)
    print ">>> "+bold(rmcommand)
    os.system(rmcommand)
  
  # GENERATE
  extraopts  = ""#%()
  gencommand = "./gridpack_generation.sh %s %s"%(sample,carddir)
  gencommand = gencommand.rstrip()
  print ">>> "+bold(gencommand)
  os.system(gencommand)
  
  # REMOVE
  if remove:
    rmdirs = [sample,carddir] if copy else [sample]
    for dir in rmdirs:
      rmcommand  = "rm -rf %s"%(dir)
      print ">>> "+bold(rmcommand)
      os.system(rmcommand)
  
  print ">>> "+'-'*100
  


if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('carddir',           type=str, action='store',
                       metavar='CARDDIR',  help="directoy with cards" )
  parser.add_argument('sample',            type=str, action='store',
                       metavar='SAMPLE',   help="name of sample" )
  parser.add_argument('-k', '--keep',      action='store_true', default=False,
                                           help="do not remove the gridpack directory" )
  #parser.add_argument('-f', '--force',     action='store_true', default=False,
  #                                         help="submit jobs without asking confirmation" )
  parser.add_argument('-c', '--copy',      action='store_true', default=False,
                                           help="copy cards into new directory to avoid stowaway cards" )
  parser.add_argument('-y', '--year',      dest='years', choices=[2016,2017,2018], type=int, nargs='+', default=[2017], action='store',
                                           help="select year" )
  parser.add_argument('-m', '--mass',      dest='masses', type=int, nargs='+', default=[ ], action='store',
                                           help="generate this mass points" )
  parser.add_argument('-t', '--tag',       type=int, nargs='?', default=-1, const=1,
                                           help="tag" )
  parser.add_argument('-T', '--test',      type=int, nargs='?', default=-1, const=1,
                      metavar='NJOBS',     help="submit test job(s)" )
  parser.add_argument('-n', '--cardlabel', type=str, action='store', default=None,
                                           help="card label replacing" )
  parser.add_argument('-p', '--param',     dest='params', default="",
                      metavar='PARAMS',    help="single string of parameters separated by colons,"+\
                                                "each with a list of values separated commas after an equal sign" )
  args = parser.parse_args()
  print ">>> "
  main(args)
  print ">>> "
  

