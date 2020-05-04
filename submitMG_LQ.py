#! /usr/bin/env python
#
# Creating dir:
#   uberftp t3se01.psi.ch 'mkdir /pnfs/psi.ch/cms/trivcat/store/user/ineuteli/samples/LowMassDiTau_madgraph'
#
# Multicore jobs:
#   to submit multicore job:    qsub -pe smp 8 ...
#   in mg_configuration.txt:    run_mode=2 # multicore
#                               nb_core=8
#   note it might wait longer in queue
#

import os, sys
import subprocess
import time
from argparse import ArgumentParser

argv = sys.argv
usage = """This script will submit jobs to generate LQ samples."""
parser = ArgumentParser(prog="submitMG_LQ",description=usage,epilog="Succes!")
parser.add_argument( '-m', '--mock-sub', dest="submit", default=True, action='store_false',
                                         help="do not submit job to batch (mock submit)")
parser.add_argument( '-s', "--signal",   dest="signals", choices=['SLQ-s','SLQ-s_NLO','SLQ-p','SLQ-p_NLO','SLQ-t','VLQ-s','VLQ-p','VLQ-t'], nargs='+', default=[ 'SLQ-s' ], action='store',
                     metavar="SIGNAL",   help="list of signals to get output for" )
parser.add_argument( '-N', '--n-events', dest="n_events", default=10000, action='store',
                                         help="number of event to be generated in each job")
parser.add_argument( '-c', '--n-cores',  dest="n_cores", type=int, default=4, action='store',
                                         help="number of core in each job")
parser.add_argument( '-q', "--queue",    dest="queue", choices=['all.q','short.q','long.q'], default='short.q', action='store',
                     metavar="QUEUE",    help="job queue" )
parser.add_argument( '-E', '--ebeam',    dest="ebeams", type=float, default=-1, action='store',
                                         help="energy per beam")
parser.add_argument( '-M', '--mass',     dest="masses", type=int, nargs='+', default=[ 1000 ], action='store',
                     metavar="MASS",     help="LQ mass point(s) to run over" )
parser.add_argument( '-p', '--parton',   dest="parton", default=False, action='store_true',
                                         help="stop the run after the parton level file generation")
parser.add_argument( '-x', '--xsec',     dest="xsec", default=False, action='store_true',
                                         help="only calculate cross section, without event generation")
parser.add_argument( '-i', '--index',    dest="indices", type=int, nargs='+', default=[ ], action='store',
                     metavar="INDEX",    help="indices to run over" )
parser.add_argument( '-t', '--tag',      dest="tag", type=str, default="", action='store',
                                         help="extra tag for output")
args = parser.parse_args()
# if len(argv) == 1:
#     parser.print_help()
#     sys.exit()

WORKPATH    = "/work/areimers/MG5_aMC_v2_7_2"
# DATAPATH    = "root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/areimers/GENSIM/"
queue       = args.queue            # all.q (10h), short.q (90 min.)
n_cores     = args.n_cores          # "nb_core" in input/mg5_configuration.txt should be the same number!
n_events    = args.n_events         # number of events to be generated
xsec        = args.xsec
signals     = args.signals
ebeams      = args.ebeams
masses      = args.masses
parton      = args.parton
tag         = args.tag
first       = 1                     # first index of a sample
last        = 1                     # last index of a sample
indices     = range(first,last+1)
nJobs       = 0
if args.indices: indices = args.indices



def main():
    global signals, n_cores, n_events, ebeams, masses, indices, parton, nJobs, xsec, tag

    # ensure directory
    REPORTDIR = "%s/LQCrossSections"%(WORKPATH)
    if not os.path.exists(REPORTDIR):
      os.makedirs(REPORTDIR)
      print ">>> made directory %s\n>>>"%(REPORTDIR)

    for sample in signals:
      for mass in masses:
        for index in indices:
          submitSample(sample,index,E=ebeams,M=mass,N=n_events,parton=parton,xsec=xsec,tag=tag)

    print ">>> submitted %s jobs"%(nJobs)



def submitSample(sample,index,M=-1,N=-1,E=-1,parton=parton,xsec=xsec,tag=tag):
    global nJobs
    jobname = "%s"%(sample)
    options = ""
    coptions = ""
    if M>0:
      jobname = "%s_M-%s"%(jobname,M)
      options = "%s -M %s"%(options,M)
    if E>0:
      options = "%s -E %s"%(options,E)
    if N>0:
      options = "%s -N %s"%(options,N)
    if tag:
      options = "%s -t %s"%(options,tag)
    if n_cores>1:
      coptions = "-pe smp %d"%(n_cores)
      options  = "%s -c %s"%(options,int(round(n_cores/2.)))
    if parton:
      options = "%s -P"%(options)
    if xsec:
      options = "%s -x"%(options)
    jobname = "%s_%d"%(jobname,index)
    options = options.lstrip(' ')
    command = "qsub -q %s %s -N %s submitMG_LQ.sh %s %s %s" % (queue,coptions,jobname,sample,index,options)
    print ">>> %s"%(command.replace(jobname,"\033[;1m%s\033[0;0m"%jobname,1))
    nJobs  += 1
    if not args.submit:
      return
    sys.stdout.write(">>> ")
    sys.stdout.flush()
    os.system(command)
    print ">>> "



def printColums(list):
    N = len(list)
    if N%4: list.extend( [" "]*(4-N%4) ); N = len(list)
    for row in zip(list[:N/4],list[N/4:N/2],list[N/2:N*3/4],list[N*3/4:]):
      print ">>> %18s %18s %18s %18s" % row



def proceed(prompt=">>> proceed?",proceed_message=">>> proceeding...",stop_message=">>> stop"):
    proceed_ = False
    while True:
      answer = raw_input(prompt+" (y or n) ").lower()
      if answer.lower() in [ 'y', 'ye', 'yes', 'yeah', 'yep', 'jep' ]:
        print proceed_message
        proceed_ = True
        break
      elif answer.lower() in [ 'n', 'no', 'na', 'nah', 'nee', 'neen', 'nop' ]:
        print stop_message
        proceed_ = False
        break
      else:
        print ">>> incorrect input"
        continue
    return proceed_



if __name__ == '__main__':
    print "\n>>> "
    main()
    print ">>> done\n"
