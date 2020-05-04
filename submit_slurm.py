#! /usr/bin/env python
# Author: Arne Reimers (March, 2020)



import os, sys, math
from os.path import isfile, join
import parse
import subprocess
import time
from argparse import ArgumentParser
from utils import ensureDirectory

m_min = 200
m_max = 3000
stepsize = 25


## What to run on
samples       = ['ScalarLQ_Single', 'ScalarLQ_Pair', 'ScalarLQ_NonRes', 'VectorLQ_Single', 'VectorLQ_Pair', 'VectorLQ_NonRes']
# samples       = ['ScalarLQ_Single']
# lambdas       = [0.5, 1.0, 1.5, 2.0, 2.5]
lambdas       = [2.0, 2.5]
masses_scalar = [800,  1100, 1500]
# masses_scalar = [1500]
masses_vector = [1000, 1500, 2000]
# masses_vector = [1500]
masses_xsec   = [25*i+200 for i in range((m_max-m_min)/stepsize + 1)]
# masses_xsec   = [575]
kappas        = [1.0, 0.0]
# kappas        = [1.0]
# tag           = '_FIRSTGEN'
tag           = '_BW15'



## How to run on it
maxindex        = 50    # Number of samples per configuration
nevents         = 3000  # Events per sample
nevents_xsec    = 10000 # Events per sample (for cross section calculation)
ncores          = 4     # Number of cores per job
submit_to_queue = None  # Override the queue to submit jobs to

gensim          = False # Submit GENSIM jobs
resubmit        = False # Submit missing GENSIM jobs
tuplize         = False # Tuplize GENSIM samples
add             = False # Hadd tuples
crosssections   = True  # calculate cross sections
resubmit_xsec   = False  # calculate cross sections

submit          = True # Actually execute the job


## Usually default options
pset  = 'pset_GENSIM.py'
queue = submit_to_queue or ('quick' if nevents<=3000 else 'wn')
if tuplize or add or crosssections:
    queue = submit_to_queue or ('quick' if nevents <= 15000 else 'wn')
max_time_option = '' if queue == 'wn' else '-t 01:00:00'




director     = 'root://t3dcachedb03.psi.ch:1094/'
WORKPATH     = "/work/areimers/GENSIM"
MGPATH       = "/work/areimers/MG5_aMC_v2_7_2"
DATAPATH     = director+"/pnfs/psi.ch/cms/trivcat/store/user/areimers/GENSIM/LQ/"
first        = 1                     # first index of a sample
indices      = range(first,maxindex+1)
nJobs        = 0
mg_bwcutoff  = 15.
lambda_ref   = 1.5



def main():
    global pset, samples, ncores, nevents, masses_scalar, masses_vector, lambdas, indices, nJobs

    INPUTDIR                = "%s/gridpacks" % (WORKPATH)
    REPORTDIR               = "%s/../workdir_slurm"%(WORKPATH)
    gpversion               = "LOPDF31%s_slc7_amd64_gcc700_CMSSW_9_3_16" % (tag)

    # ensure directory
    ensureDirectory(REPORTDIR)
    if crosssections: ensureDirectory(REPORTDIR)


    # Create commands
    for sample in samples:
      masses = masses_scalar if 'Scalar' in sample else masses_vector
      if 'Vector' in sample and -1 in kappas: raise ValueError('Sample is vector but no list of kappas is provided')

      if 'Scalar' in sample: mykappas = [kappas[0]]
      else: mykappas = kappas

      for mass in masses:
        for lamb in lambdas:
          for kappa in mykappas:
            lambstr = str(lamb).replace('.', 'p')
            kappastr = str(kappa).replace('.', 'p')
            gridpack = "%s/%s_M%s_L%s_%s_tarball.tar.xz"%(INPUTDIR,sample,mass,lambstr,gpversion)
            jobname  = "%s_M%s_L%s"%(sample,mass,lambstr)
            if 'Vector' in sample:
              gridpack = "%s/%s_M%s_L%s_K%s_%s_tarball.tar.xz"%(INPUTDIR,sample,mass,lambstr,kappastr,gpversion)
              jobname  = "%s_M%s_L%s_K%s"%(sample,mass,lambstr,kappastr)
            jobname += tag
            nJobs = 0
            f = open('commands/commands_submit_slurm_%s.txt'%(jobname), 'w')
            for index in indices:
              command = getCommandLine(pset,gridpack,jobname,index,N=nevents)
              f.write(command + '\n')
              nJobs += 1
            f.close
            ft = open('commands/commands_submit_slurm_%s_tuplize.txt'%(jobname), 'w')
            nJobs_tup = 0
            for index in indices:
              infilename = '/pnfs/psi.ch/cms/trivcat/store/user/areimers/GENSIM/LQ/' + jobname + '/GENSIM_%d.root'%(index)
              outfilename = 'simple_files/' + jobname + '_GENSIM_simple_%d.root'%(index)
              command = './convertGENSIM.py %s -o %s' % (infilename, outfilename)
              # print command
              ft.write(command + '\n')
              nJobs_tup += 1
    ft.close

    for sample in samples:

      if crosssections or resubmit_xsec:
        masses = masses_xsec
      else:
        masses = masses_scalar if 'Scalar' in sample else masses_vector

      if 'Scalar' in sample: mykappas = [kappas[0]]
      else: mykappas = kappas

      for mass in masses:
        samplemass = "%s_M%s" % (sample, mass)
        for lamb in lambdas:

          this_cutoff    = getNewCutoff(sample=sample, mg_bwcutoff=mg_bwcutoff, lamb_ref=lambda_ref, lamb=lamb, mass=mass)
          if 'BW15' in tag:
            this_cutoff    = mg_bwcutoff

          # if sample has smaller lambda than reference-sample, only calculate the cross section in it's default window, as the lambda-reweighting factors only cover the region where the narrower sample is generated (all samples are generated at LO inside their default windows)
          if lamb < lambda_ref: this_cutoff = mg_bwcutoff

          for kappa in mykappas:

            lambstr = str(lamb).replace('.', 'p')
            kappastr = str(kappa).replace('.', 'p')
            jobname  = "%s_M%s_L%s"%(sample,mass,lambstr)
            kappaopt = ''
            if 'Vector' in sample:
              jobname  = "%s_M%s_L%s_K%s"%(sample,mass,lambstr,kappastr)
              kappaopt = "-K %s" % (kappa)
            jobname += tag



            # MAIN FUNCTIONS #
            # ============== #

            # command_xsec = 'sbatch -J xsecs_%s -p %s --cpus-per-task 1 submitMG_LQ.sh %s 1 -M %s -N %s' % (samplemass, queue, sample, mass, nevents_xsec)
            # command_xsec = 'sbatch -J xsecs_%s -p %s --cpus-per-task 2 submitMG_LQ.sh %s -c 2 -M %s -L %s %s -N %s' % (jobname, queue, sample, mass, lamb, kappaopt, nevents_xsec)
            # command_xsec = 'sbatch -J xsecs_%s -p %s --cpus-per-task 2 submitMG_LQ.sh %s -c 2 -M %s -L %s %s -N %s' % (jobname, queue, sample, mass, lamb, kappaopt, nevents_xsec)
            command_xsec = 'sbatch -J xsecs_%s -p %s --cpus-per-task 2 submitMG_LQ.sh %s -M %s -L %s -B %.2f %s -N %s -T %s' % (jobname, queue, sample, mass, lamb, this_cutoff, kappaopt, nevents_xsec, tag)
            if crosssections:
              print command_xsec
              if submit: os.system(command_xsec)
            elif resubmit_xsec:
              do_resubmit = isMissingCrossSections(jobname)
              if do_resubmit:
                print 'need to resubmit jobname %s with command %s.' % (jobname, command_xsec)
                if submit: os.system(command_xsec)




            elif resubmit:
              # Check GENSIM rootfiles that are stored, compare with expectation, resubmit missing ones
              missing_indices = findMissingFilesGENSIM(jobname, maxindex)
              nJobs_re = 0
              fr = open('commands/commands_submit_slurm_%s_resubmit.txt'%(jobname), 'w')
              for index in missing_indices:
                command = getCommandLine(pset,gridpack,jobname,index,N=nevents)
                fr.write(command + '\n')
                nJobs_re += 1
              fr.close
              # submit, this time make sure it's to the 'wn' partition
              command = 'sbatch -a 1-%s -J resubmit_%s -p wn --cpus-per-task %i submit_slurm_interface.sh commands/commands_submit_slurm_%s_resubmit.txt' % (str(nJobs_re), jobname, ncores, jobname)
              print command
              if submit: os.system(command)
              print ">>> re-submitted an array of %s jobs"%(str(nJobs_re))
            elif tuplize:
              command = 'sbatch -a 1-%s -J tuplize_%s -p %s --cpus-per-task 1 tuplize_gensim.sh commands/commands_submit_slurm_%s_tuplize.txt' % (str(nJobs_tup), jobname, queue, jobname)
              print command
              if submit: os.system(command)
              print ">>> submitted an array of %s jobs for tuplizing gensim. "%(str(nJobs_tup))
            elif add:
              command = 'sbatch -J add_%s -p quick --cpus-per-task 1 add_tuples.sh %s' % (jobname, jobname)
              print command
              if submit: os.system(command)
              print ">>> submitted job for adding gensim samples: %s. "%(jobname)
            elif gensim:
              command = 'sbatch -a 1-%s -J gensim_%s -p %s %s --cpus-per-task %i submit_slurm_interface.sh commands/commands_submit_slurm_%s.txt' % (str(nJobs), jobname, queue, max_time_option, ncores, jobname)
              print command
              if submit: os.system(command)
              print ">>> submitted an array of %s jobs"%(str(nJobs))
            elif not crosssections:
              raise ValueError("No run-option specified - what should be submitted?")




def getSampleNameFromGridpack(gridpack):
  gridpack = os.path.basename(gridpack)
  for pattern in ['_slc6_','_slc7_','_CMSSW_',]:
    if pattern in gridpack:
      gridpack = gridpack[:gridpack.index(pattern)]
      break
  return gridpack


def findMissingFilesGENSIM(jobname, maxindex):
  missing_indices = []
  filename_base = '/pnfs/psi.ch/cms/trivcat/store/user/areimers/GENSIM/LQ/%s/GENSIM_' % (jobname)
  for idx in range(1,maxindex+1):
    filename = filename_base + str(idx) + '.root'
    if not os.path.isfile(filename):
      missing_indices.append(idx)
  return missing_indices


def isMissingCrossSections(jobname):
    # if shortfile exists, job was successful, otherwise it wasn't
    crosssection_filepath = '/work/areimers/LQCrossSections/'

    # look for everything that doesn't have '_short' in it
    filenames = [f for f in os.listdir(crosssection_filepath) if isfile(join(crosssection_filepath, f)) and not '_short' in f]
    filenames_short = [f for f in os.listdir(crosssection_filepath) if isfile(join(crosssection_filepath, f)) and '_short' in f]
    # print filenames_short

    resub = False if len(filenames_short) > 0 else True
    if not resub:
        found_in_filenames = False
        for f in filenames_short:
            if jobname in f:
                found_in_filenames = True
        resub = not found_in_filenames


    # print jobname, resub
    return resub


def getCommandLine(pset,gridpack,sample,index=-1,N=-1,ncores=ncores):
    """Submit PSet config file and gridpack to SLURM batch system."""
    global nJobs
    jobname  = sample
    gridpack = os.path.abspath(gridpack)
    options  = ""
    if index>0:
      jobname = "%s_%d"%(jobname,index)
      options = "%s -i %d"%(options,index)
    if N>0:
      options = "%s -N %s"%(options,N)
    if ncores>1:
      options  = "%s -c %s"%(options,ncores)
    options = options.lstrip(' ')
    command = "./submit_slurm.sh %s %s %s %s"%(pset,gridpack,sample,options.strip())
    return command

def getBWCutoffs(mg_bwcutoff=15., lamb=1.5, mass=1000, spin="Scalar"):
    m_upper = (1+mg_bwcutoff*lamb*lamb/(16*math.pi)) * mass;
    if spin == "Vector": m_upper = (1+mg_bwcutoff*lamb*lamb/(24*math.pi)) * mass;
    m_lower = (1-mg_bwcutoff*lamb*lamb/(16*math.pi)) * mass;
    if spin == "Vector": m_lower = (1-mg_bwcutoff*lamb*lamb/(24*math.pi)) * mass;

    return(m_lower, m_upper)

def getNewCutoff(sample, mg_bwcutoff, lamb_ref, lamb, mass):
    spin = "Scalar" if "Scalar" in sample else "Vector"
    (default_cutoff_lo, default_cutoff_hi) = getBWCutoffs(mg_bwcutoff=mg_bwcutoff, lamb=lamb_ref, mass=mass, spin=spin)
    bw_window_width_default                = default_cutoff_hi - default_cutoff_lo
    (this_cutoff_lo, this_cutoff_hi)       = getBWCutoffs(mg_bwcutoff=mg_bwcutoff, lamb=lamb, mass=mass, spin=spin)
    bw_window_width_default_this           = this_cutoff_hi - this_cutoff_lo
    bw_window_factor                       = float(bw_window_width_default) / float(bw_window_width_default_this)
    this_cutoff                            = mg_bwcutoff * bw_window_factor
    return this_cutoff

def is_file_empty(file_path):
    """ Check if file is empty by confirming if its size is 0 bytes"""
    # Check if file exist and it is empty
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0

if __name__ == '__main__':
  print "\n>>> "
  main()
  print ">>> done\n"
