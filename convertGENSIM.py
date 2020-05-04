#! /usr/bin/env python
# Author: Yuta Takahashi
# Description: converts GENSIM sample to a simplified format for generator-level analysis
# http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
# ./convertGENSIM.py /pnfs/psi.ch/cms/trivcat/store/user/ineuteli/samples/SLQ-s_M-1000_MGP_GENSIM/*GENSIM*.root -o checkGENSIM/SLQ-s_M-1000_MGP_LOPDF_CP2_GENSIM_simple2.root -N 100
# ./convertGENSIM.py /pnfs/psi.ch/cms/trivcat/store/user/ineuteli/samples/SLQ-s_M-1000_MGonly_GENSIM/*GENSIM*.root -o checkGENSIM/SLQ-s_M-1000_MGonly_LOPDF_CP2_GENSIM_simple2.root -N 100
import time
start = time.time()
import os, sys, glob, copy, math
from math import log, ceil, exp
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TFile, TH1F, TTree, gStyle, TH2F
from DataFormats.FWLite import Events, Handle
import numpy as num
from argparse import ArgumentParser
usage = """Converts GENSIM sample to a simplified format for generator-level analysis."""
parser = ArgumentParser(prog="fragment_GS", description=usage, epilog="Success!")
parser.add_argument('infiles',         nargs='*',
                    metavar='INFILE',  help="signal to generate events for" )
parser.add_argument('-o', '--out',     dest='outfile', type=str, default=None, action='store',
                    metavar='OUTFILE', help="output file" )
parser.add_argument('-d', '--outdir',  type=str, default=None, action='store',
                                       help="output file" )
parser.add_argument('-n', '--nmax',    dest='Nmax', type=int, default=-1, action='store',
                                       help="number of event to be generated in each job" )
parser.add_argument('-p', '--prod',    dest='prods', choices=[p+'-'+t for p in ['SLQ','VLQ'] for t in ['s','p','t']], nargs='*', default=[ ], action='store',
                                       help="production modes" )
parser.add_argument('-m', '--mass',    dest='masses', type=int, nargs='*', default=[ ], action='store',
                                       help="mass values" )
parser.add_argument('-l', '--lambda',  dest='lambdas', type=float, nargs='*', default=[1.0], action='store',
                                       help="lambda values" )
parser.add_argument('-k', '--kappa',   dest='kappas', type=float, nargs='*', default=[1.0], action='store',
                                       help="kappa values" )
parser.add_argument('-t', '--tag',     type=str, default="", action='store',
                                       help="extra tag for output" )
args      = parser.parse_args()

gStyle.SetOptStat(11111)
decaydict = { 'ele': 0, 'muon': 1, 'tau': 2 }



def main():

  # SETTING
  director  = 'root://t3dcachedb.psi.ch:1094/'
  indir     = "/pnfs/psi.ch/cms/trivcat/store/user/areimers/GENSIM/LQ"
  outdir    = ensureDir(args.outdir or "GENSIM")
  Nmax      = args.Nmax
  tag       = args.tag
  infileset = [ ]

  # USER FILES
  if args.infiles:
    infiles = [director+f if f.startswith('/pnfs/psi.ch/cms') else f for f in args.infiles]
    if args.outfile:
      outfile = args.outfile
    else:
      outfile = infile.split('/')[-1]
      if 'GENSIM' in outfile:
        outfile = outfile[:outfile.rfind('GENSIM')]+"GENSIM_simple%s.root"%tag
      else:
        outfile = outfile.replace(".root","_simple%s.root"%tag)
    infileset.append(("User input",infiles,outfile))

  # USER PARAMETERS
  elif args.masses and args.prods:
    for prod in args.prods:
      kappas = [''] if 'SLQ' in prod else args.kappas
      model  = 'Scalar' if 'SLQ' in prod else 'Vector'
      prod   = 'NonRes' if '-t' in prod else 'Pair' if '-p' in prod else 'Single'
      sample = "%sLQ_%s"%(model,prod)
      for mass in args.masses:
        if model=='Scalar' and mass in [1000,2000]: continue
        if model=='Vector' and mass in [ 800,1100]: continue
        for kappa in kappas:
          for lambd in args.lambdas:
            title    = "%s: M=%s, lambda=%s"%(sample,mass,lambd)
            if kappa!='':
              title += ", kappa=%s"%(kappa)
            lambd    = str(lambd).replace('.','p')
            kappa    = str(kappa).replace('.','p')
            if kappa:
              infile = "%s/%s_M%s_L%s_K%s_GENSIM/GENSIM_*.root"%(indir,sample,mass,lambd,kappa)
              outfile  = "%s/%s_M%s_L%s_K%s_GENSIM_simple%s.root"%(outdir,sample,mass,lambd,kappa,tag)
            else:
              infile = "%s/%s_M%s_L%s_GENSIM/GENSIM_*.root"%(indir,sample,mass,lambd)
              outfile  = "%s/%s_M%s_L%s_GENSIM_simple%s.root"%(outdir,sample,mass,lambd,tag)
            infiles  = [director+f for f in glob.glob(infile) if 'LHE' not in f]
            if not infiles:
              print ">>> Warning! No files found for '%s' !!!"%(infile)
            infileset.append((title,infiles,outfile))
  else:
    raise IOError("No input given to run on.")

  # CONVERT
  for title, infiles, outfile in infileset:
    print ">>> %s"%bold(title)
    print ">>>   input files:"
    for file in infiles:
      print ">>>     %s"%(file)
    print ">>>   output file: %s"%(outfile)
    isPythia  = 'Pythia' in infiles[0]
    convertGENSIM(infiles,outfile,Nmax=Nmax,isPythia=isPythia)
    print ">>> "

  print ">>>   done in in %s"%(formatTime(time.time()-start))


def convertGENSIM(infiles,outfilename,Nmax=-1,isPythia=False):
  """Loop over GENSIM events and save custom trees."""
  start1 = time.time()

  lqids = [46] if isPythia else [9000002,9000006]

  print ">>>   loading files..."
  events  = Events(infiles)
  outfile = TFile(outfilename, 'RECREATE')

  print ">>>   creating trees and branches..."
  tree_event  = TTree('event', 'event')
  tree_jet    = TTree('jet',   'jet')
  tree_mother = TTree('mother','mother')
  tree_decay  = TTree('decay', 'decay')
  tree_assoc  = TTree('assoc', 'assoc')

  # EVENT
  tree_event.addBranch('nbgen',       'i')
  tree_event.addBranch('nbcut',       'i')
  tree_event.addBranch('nbcut50',     'i')
  tree_event.addBranch('ntgen',       'i')
  tree_event.addBranch('njet',        'i')
  tree_event.addBranch('nlepton',     'i')
  tree_event.addBranch('ntau',        'i')
  tree_event.addBranch('ntaucut',     'i')
  tree_event.addBranch('ntaucut50',   'i')
  tree_event.addBranch('ntaucut_vis', 'i')
  tree_event.addBranch('ntaucut50_vis','i')
  tree_event.addBranch('nnu',         'i')
  tree_event.addBranch('nlq',         'i')
  tree_event.addBranch('ntau_assoc',  'i')
  tree_event.addBranch('ntau_decay',  'i')
  tree_event.addBranch('nbgen_decay', 'i')
  tree_event.addBranch('met',         'f')
  tree_event.addBranch('jpt1',        'f')
  tree_event.addBranch('jpt2',        'f')
  tree_event.addBranch('sumjet',      'f')
  tree_event.addBranch('dphi_jj',     'f')
  tree_event.addBranch('deta_jj',     'f')
  tree_event.addBranch('dr_jj',       'f')
  tree_event.addBranch('ncentral',    'i')
  tree_event.addBranch('mjj',         'f')
  tree_event.addBranch('lq1_mass',    'f')
  tree_event.addBranch('lq2_mass',    'f')
  tree_event.addBranch('lq1_pt',      'f')
  tree_event.addBranch('lq2_pt',      'f')
  tree_event.addBranch('m_lqlq',      'f')
  tree_event.addBranch('tau1_pt',     'f')
  tree_event.addBranch('tau1_eta',    'f')
  tree_event.addBranch('tau1_y',      'f')
  tree_event.addBranch('tau1_ptvis',  'f')
  tree_event.addBranch('tau1_etavis', 'f')
  tree_event.addBranch('tau2_pt',     'f')
  tree_event.addBranch('tau2_eta',    'f')
  tree_event.addBranch('tau2_y',      'f')
  tree_event.addBranch('tau2_ptvis',  'f')
  tree_event.addBranch('tau2_etavis', 'f')
  tree_event.addBranch('ditau_dy',    'f')
  tree_event.addBranch('ditau_chi',   'f')
  tree_event.addBranch('st',          'f') # scalar sum pT
  tree_event.addBranch('st_met',      'f') # scalar sum pT with MET
  tree_event.addBranch('weight',      'f')

  # LQ DECAY
  tree_mother.addBranch('pid',        'i')
  tree_mother.addBranch('moth',       'i')
  tree_mother.addBranch('status',     'i')
  tree_mother.addBranch('pt',         'f')
  tree_mother.addBranch('eta',        'f')
  tree_mother.addBranch('phi',        'f')
  tree_mother.addBranch('mass',       'f')
  tree_mother.addBranch('inv',        'f')
  tree_mother.addBranch('ndau',       'i')
  tree_mother.addBranch('dau',        'i')
  tree_mother.addBranch('dphi_ll',    'f')
  tree_mother.addBranch('deta_ll',    'f')
  tree_mother.addBranch('dr_ll',      'f')
  tree_mother.addBranch('st',         'f') # scalar sum pT
  tree_mother.addBranch('st_met',     'f') # scalar sum pT with MET
  tree_mother.addBranch('weight',     'f')

  # FROM LQ DECAY
  tree_decay.addBranch('pid',         'i')
  tree_decay.addBranch('pt',          'f')
  tree_decay.addBranch('eta',         'f')
  tree_decay.addBranch('phi',         'f')
  tree_decay.addBranch('lq_mass',     'f')
  tree_decay.addBranch('ptvis',       'f')
  tree_decay.addBranch('type',        'i')
  tree_decay.addBranch('isBrem',      'i')
  tree_decay.addBranch('weight',      'f')

  # NOT FROM LQ DECAY (ASSOCIATED)
  tree_assoc.addBranch('pid',         'i')
  tree_assoc.addBranch('moth',        'i')
  tree_assoc.addBranch('pt',          'f')
  tree_assoc.addBranch('ptvis',       'f')
  tree_assoc.addBranch('eta',         'f')
  tree_assoc.addBranch('phi',         'f')
  tree_assoc.addBranch('weight',      'f')

  # JETS
  tree_jet.addBranch('pt',            'f')
  tree_jet.addBranch('eta',           'f')
  tree_jet.addBranch('phi',           'f')
  tree_jet.addBranch('weight',        'f')

  hist_LQ_decay = TH1F('LQ_decay',"LQ decay",60,-30,30)
  handle_gps,    label_gps    = Handle('std::vector<reco::GenParticle>'), 'genParticles'
  handle_jets,   label_jets   = Handle('std::vector<reco::GenJet>'), 'ak4GenJets'
  handle_met,    label_met    = Handle('vector<reco::GenMET>'), 'genMetTrue'
  handle_weight, label_weight = Handle('GenEventInfoProduct'), 'generator'

  evtid = 0
  sec_per_evt = 0.023 # seconds per event
  Ntot = Nmax if Nmax>0 else events.size()
  print ">>>   start processing %d events, ETA %s..."%(Ntot,formatTimeShort(sec_per_evt*Ntot))
  step = stepsize(Ntot)
  start_proc = time.time()

  # LOOP OVER EVENTS
  for event in events:
      # print ' --- NEW EVENT'
      # print '='*30
      # print evtid
      if Nmax>0 and evtid>=Nmax: break
      if evtid>0 and evtid%step==0: print ">>>     processed %4s/%d events, ETA %s"%(evtid,Ntot,ETA(start_proc,evtid+1,Ntot))
      evtid += 1

      event.getByLabel(label_gps,handle_gps)
      gps = handle_gps.product()

      event.getByLabel(label_jets,handle_jets)
      jets = handle_jets.product()

      event.getByLabel(label_met,handle_met)
      met = handle_met.product()

      event.getByLabel(label_weight,handle_weight)
      gweight = handle_weight.product()
      weight = gweight.weight()

      # GEN PARTICLES
      gps_mother = [p for p in gps if isFinal(p) and abs(p.pdgId()) in [42]]
      gps_final  = [p for p in gps if isFinal(p) and abs(p.pdgId()) in [5,6,15,16]+lqids]
      gps_mother = [p for p in gps_final if abs(p.pdgId()) in lqids and p.status()>60] #not(moth.numberOfDaughters()==2 and abs(moth.daughter(0).pdgId()) in lqids)
      gps_bgen   = [p for p in gps_final if abs(p.pdgId())==5 and p.status()==71]
      gps_bcut   = [p for p in gps_bgen  if p.pt()>20 and abs(p.eta())<2.5]
      gps_bcut50 = [p for p in gps_bgen  if p.pt()>50 and abs(p.eta())<2.5]
      gps_tgen   = [p for p in gps_final if abs(p.pdgId())==6] #[-1:]
      gps_nugen  = [p for p in gps_final if abs(p.pdgId())==16]
      gps_tau    = [p for p in gps_final if abs(p.pdgId())==15 and p.status()==2]
      gps_tau.sort(key=lambda p: p.pt(), reverse=True)
      gps_taucut = [p for p in gps_tau   if p.pt()>20 and abs(p.eta())<2.5]
      gps_taucut50 = [p for p in gps_tau   if p.pt()>50 and abs(p.eta())<2.5]

      gps_tau_vis = []
      gps_taucut_vis = []
      gps_taucut50_vis = []
      # find taus that survive ptvis > 20 and >50 cuts
      for p in gps_tau:
        while p.status()!=2 :
          p = p.daughter(0)
        findau = finalDaughters(p, [])
        thisptvis = p4sumvis(findau).pt()
        thisetavis = p4sumvis(findau).eta()
        gps_tau_vis.append(p)
        if thisptvis > 20 and abs(thisetavis) < 2.5: gps_taucut_vis.append(p)
        if thisptvis > 50 and abs(thisetavis) < 2.5: gps_taucut50_vis.append(p)


      #print '-'*10
      #for p in gps_tgen:
      #  printParticle(p)
      #if gps_tgen:
      #  print "has top"
      #for p in gps_nugen:
      #  printParticle(p)

      # REMOVE TOP QUARK if its final daughter is also in the list
      for top in gps_tgen[:]:
        dau = top
        while abs(dau.daughter(0).pdgId())==6:
          dau = dau.daughter(0)
        if dau!=top and dau in gps_tgen:
          gps_tgen.remove(top)

      # REMOVE JET-LEPTON OVERLAP
      jets, dummy = cleanObjectCollection(jets,gps_tau,dRmin=0.5)
      njets  = 0
      sumjet = 0
      jets30  = [ ]
      for jet in jets:
        if jet.pt()>30 and abs(jet.eta())<5:
          sumjet            += jet.pt()
          njets             += 1
          tree_jet.pt[0]     = jet.pt()
          tree_jet.eta[0]    = jet.eta()
          tree_jet.phi[0]    = jet.phi()
          tree_jet.weight[0] = weight
          tree_jet.Fill()
          jets30.append(jet)

      # MULTIPLICITIES
      tree_event.nlq[0]           = len(gps_mother)
      tree_event.nbcut[0]         = len(gps_bcut)
      tree_event.nbcut50[0]       = len(gps_bcut50)
      tree_event.nbgen[0]         = len(gps_bgen)
      tree_event.ntgen[0]         = len(gps_tgen)
      tree_event.njet[0]          = njets
      tree_event.nlepton[0]       = len(gps_tau)
      tree_event.ntau[0]          = len(gps_tau)
      tree_event.ntaucut[0]       = len(gps_taucut)
      tree_event.ntaucut50[0]     = len(gps_taucut50)
      tree_event.ntaucut_vis[0]   = len(gps_taucut_vis)
      tree_event.ntaucut50_vis[0] = len(gps_taucut50_vis)
      tree_event.nnu[0]           = len(gps_nugen)

      # JETS
      tree_event.met[0]     = met[0].pt()
      tree_event.sumjet[0]  = sumjet
      if len(jets30)>=2:
        centrajpt1s = findCentrajpt1s(jets30[:2],jets30[2:])
        tree_event.ncentral[0] = len(centrajpt1s)
      else:
        tree_event.ncentral[0] = -9
      if(len(jets30)>=2):
          tree_event.jpt1[0]    = jets30[0].pt()
          tree_event.jpt2[0]    = jets30[1].pt()
          tree_event.dphi_jj[0] = deltaPhi(jets30[0].phi(), jets30[1].phi())
          tree_event.deta_jj[0] = jets30[0].eta() - jets30[1].eta()
          tree_event.dr_jj[0]   = deltaR(jets30[0].eta(),jets30[0].phi(),jets30[1].eta(),jets30[1].phi())
          dijetp4               = jets30[0].p4() + jets30[1].p4()
          tree_event.mjj[0]     = dijetp4.M()
      elif (len(jets30)==1):
          tree_event.jpt1[0]    = jets30[0].pt()
          tree_event.jpt2[0]    = -1
          tree_event.dphi_jj[0] = -9
          tree_event.deta_jj[0] = -9
          tree_event.dr_jj[0]   = -1
          tree_event.mjj[0]     = -1
      else:
          tree_event.jpt1[0]    = -1
          tree_event.jpt2[0]    = -1
          tree_event.dphi_jj[0] = -9
          tree_event.deta_jj[0] = -9
          tree_event.dr_jj[0]   = -1
          tree_event.mjj[0]     = -1

      # SCALAR SUM PT
      # if len(gps_taucut)>=2 and len(gps_bcut)>=1:
      if len(gps_tau_vis)>=2 and len(gps_bcut)>=1:
        st = 0
        #gps_taucut.sort(key=lambda p: p.pt(), reverse=True)
        gps_bcut.sort(key=lambda p: p.pt(), reverse=True)
        #taus_assoc.sort(key=lambda p: p.pt(), reverse=True)
        #taus_decay.sort(key=lambda p: p.pt(), reverse=True)
        #bgen_decay.sort(key=lambda p: p.pt(), reverse=True)
        for tau in gps_tau_vis[:2]:
          st += p4sumvis(finalDaughters(tau, [])).pt()
        st += gps_bcut[0].pt()
        # for part in gps_taucut[:2]+gps_bcut[:1]:
        #   st += part.pt()
        stmet = st + met[0].pt()
      else:
        st    = -1
        stmet = -1

      if len(gps_tau) > 0:
        tree_event.tau1_pt[0]      = gps_tau[0].pt()
        tree_event.tau1_eta[0]     = gps_tau[0].eta()
        tree_event.tau1_y[0]       = gps_tau[0].p4().Rapidity()
        tree_event.tau1_ptvis[0]   = p4sumvis(finalDaughters(gps_tau_vis[0], [])).pt()
        tree_event.tau1_etavis[0]  = p4sumvis(finalDaughters(gps_tau_vis[0], [])).eta()
      if len(gps_tau) > 1:
        tree_event.tau2_pt[0]      = gps_tau[1].pt()
        tree_event.tau2_eta[0]     = gps_tau[1].eta()
        tree_event.tau2_y[0]       = gps_tau[1].p4().Rapidity()
        tree_event.tau2_ptvis[0]   = p4sumvis(finalDaughters(gps_tau_vis[1], [])).pt()
        tree_event.tau2_etavis[0]  = p4sumvis(finalDaughters(gps_tau_vis[1], [])).eta()
        dy = abs(gps_tau[0].p4().Rapidity() - gps_tau[1].p4().Rapidity())
        tree_event.ditau_dy[0]     = dy
        tree_event.ditau_chi[0]    = exp(dy)
        
      tree_event.st[0]           = st
      tree_event.st_met[0]       = stmet
      tree_mother.st[0]          = st
      tree_mother.st_met[0]      = stmet

      tree_event.weight[0] = weight

      #print 'len, gps_mother = ', len(gps_mother)
      #if len(gps_mother)==1:
      #    print gps_mother[0].pdgId(), gps_mother[0].status(), gps_mother[0].pt(), gps_mother[0].eta(), gps_mother[0].phi()
      #    print '1 (ndaughter, daughter pdgid) =', gps_mother[0].numberOfDaughters(), gps_mother[0].daughter(0).pdgId(), '(pdgId, status, pt, eta, phi) = ', gps_mother[0].pdgId(), gps_mother[0].status(), gps_mother[0].pt(), gps_mother[0].eta(), gps_mother[0].phi()
      #if len(gps_mother)>=2:
      #    print '2 (ndaughter, daughter 1/2 pdgid) =', gps_mother[0].numberOfDaughters(), gps_mother[0].daughter(0).pdgId(), gps_mother[0].daughter(1).pdgId(), '(pdgId, status, pt, eta, phi) = ', gps_mother[0].pdgId(), gps_mother[0].status(), gps_mother[0].pt(), gps_mother[0].eta(), gps_mother[0].phi()
      #    print '2 (ndaughter, daughter 1/2 pdgid) =', gps_mother[1].numberOfDaughters(), gps_mother[1].daughter(0).pdgId(), gps_mother[1].daughter(1).pdgId(), '(pdgId, status, pt, eta, phi) = ', gps_mother[1].pdgId(), gps_mother[1].status(), gps_mother[1].pt(), gps_mother[1].eta(), gps_mother[1].phi()

      # TAU
      taus_assoc = [ ]
      for gentau in gps_tau:

          while gentau.status()!=2 :
            gentau = gentau.daughter(0)
          genfinDaughters = finalDaughters(gentau, [])
          genptvis = p4sumvis(genfinDaughters).pt()

          # CHECK MOTHER
          taumoth  = gentau.mother(0)
          mothpid  = abs(taumoth.pdgId())
          from_LQ  = False
          #from_had = False # from hadron decay
          #print '-'*30
          while mothpid!=2212:
            #print taumoth.pdgId()
            if mothpid in lqids:
              from_LQ = True
              break
            elif 100<mothpid<10000: #and mothpid!=2212:
              #from_had = True
              break
            taumoth = taumoth.mother(0)
            mothpid = abs(taumoth.pdgId())

          # ASSOC
          if not from_LQ:
            tree_assoc.pt[0]     = gentau.pt()
            tree_assoc.ptvis[0]  = genptvis
            tree_assoc.eta[0]    = gentau.eta()
            tree_assoc.phi[0]    = gentau.phi()
            tree_assoc.pid[0]    = gentau.pdgId()
            tree_assoc.moth[0]   = taumoth.pdgId()
            tree_assoc.weight[0] = weight
            tree_assoc.Fill()
            #if not from_had:
            taus_assoc.append(gentau)

      # B QUARK
      for genb in gps_bgen:
          bmoth   = genb.mother(0)
          mothpid = abs(bmoth.pdgId())
          from_LQ = False
          while mothpid!=2212:
            if mothpid in lqids:
              from_LQ = True
              break
            bmoth   = bmoth.mother(0)
            mothpid = abs(bmoth.pdgId())
          if not from_LQ:
            tree_assoc.pt[0]     = genb.pt()
            tree_assoc.ptvis[0]  = -1
            tree_assoc.eta[0]    = genb.eta()
            tree_assoc.phi[0]    = genb.phi()
            tree_assoc.pid[0]    = genb.pdgId()
            tree_assoc.moth[0]   = bmoth.pdgId()
            tree_assoc.weight[0] = weight
            tree_assoc.Fill()

      # MOTHER LQ
      #print '-'*80
      taus_decay = [ ]
      bgen_decay = [ ]
      gps_mother.sort(key=lambda p: p.pt(), reverse=True)
      for moth in gps_mother:

          dau_pid = 0
          pair    = [ ]

          if moth.numberOfDaughters()==2:
            if moth.daughter(0).pdgId() in [21,22] or moth.daughter(1).pdgId() in [21,22]:
              continue
            if abs(moth.daughter(0).pdgId()) in lqids: # single production with t-channel LQ
              continue

          lq_moth = moth.mother(0)
          while abs(lq_moth.pdgId()) in lqids:
            lq_moth = lq_moth.mother(0)

          for i in range(moth.numberOfDaughters()):
              #print '\t', dau.pdgId()
              dau = moth.daughter(i)

              # TAU
              isBrem = False
              if abs(dau.pdgId())==15:
                while dau.status()!=2:
                  dau = dau.daughter(0)
                if dau.numberOfDaughters()==2 and abs(dau.daughter(0).pdgId())==15 and dau.daughter(1).pdgId()==22:
                    #print "This is brems !?!"
                    isBrem = True
                else:
                  taus_decay.append(dau)

              # BOTTOM QUARK
              elif abs(dau.pdgId())==5:
                dau_pid = dau.pdgId()
                bgen_decay.append(dau)

              # TOP QUARK
              elif abs(dau.pdgId())==6:
                dau_pid = dau.pdgId()
                newdau = dau
                while abs(newdau.daughter(0).pdgId())==6:
                  newdau = newdau.daughter(0)
                if isFinal(newdau):
                  dau = newdau

              pair.append(dau.p4())
              tree_decay.lq_mass[0] = moth.mass()
              tree_decay.pid[0]     = dau.pdgId()
              tree_decay.pt[0]      = dau.pt()
              tree_decay.eta[0]     = dau.eta()
              tree_decay.phi[0]     = dau.phi()
              tree_decay.isBrem[0]  = isBrem

              if abs(dau.pdgId())==15:
                finDaughters        = finalDaughters(dau, [ ])
                ptvis               = p4sumvis(finDaughters).pt()
                tree_decay.ptvis[0] = ptvis
                decaymode           = tauDecayMode(dau)
                tree_decay.type[0]  = decaydict[decaymode]
                #print decaymode, 'vis pt = ', ptvis , 'tau pt = ', dau.pt()
                if ptvis > dau.pt():
                  print "%s, vis pt = %s, tau pt = %s "%(decaymode,ptvis,dau.pt())+'!'*30
              else:
                tree_decay.ptvis[0] = dau.pt()
                tree_decay.type[0]  = -1
              tree_decay.weight[0]  = weight
              tree_decay.Fill()

              if abs(moth.pdgId()) in lqids:
                hist_LQ_decay.Fill(dau.pdgId())

          if len(pair)==2:
            tree_mother.inv[0]     = (pair[0] + pair[1]).mass()
            tree_mother.dphi_ll[0] = deltaPhi(pair[0].phi(), pair[1].phi())
            tree_mother.deta_ll[0] = pair[0].eta() - pair[1].eta()
            tree_mother.dr_ll[0]   = deltaR(pair[0].eta(),pair[0].phi(),pair[1].eta(),pair[1].phi())
          else:
            tree_mother.inv[0]     = -1
            tree_mother.dphi_ll[0] = -99
            tree_mother.deta_ll[0] = -99
            tree_mother.dr_ll[0]   = -99

          tree_mother.pid[0]       = moth.pdgId()
          tree_mother.moth[0]      = lq_moth.pdgId()
          tree_mother.status[0]    = moth.status()
          tree_mother.mass[0]      = moth.mass()
          tree_mother.pt[0]        = moth.pt()
          tree_mother.eta[0]       = moth.eta()
          tree_mother.phi[0]       = moth.phi()
          tree_mother.ndau[0]      = len(pair)
          tree_mother.dau[0]       = dau_pid # save PDG ID for quark daughter
          tree_mother.weight[0]    = weight
          tree_mother.Fill()

      if len(gps_mother)==1:
        tree_event.lq1_mass[0]  = gps_mother[0].mass()
        tree_event.lq1_pt[0]    = gps_mother[0].pt()
        tree_event.lq2_mass[0]  = -1
        tree_event.lq2_pt[0]    = -1
        tree_event.m_lqlq[0]    = -1
      elif len(gps_mother)>=2:
        tree_event.lq1_mass[0]  = gps_mother[0].mass()
        tree_event.lq1_pt[0]    = gps_mother[0].pt()
        tree_event.lq2_mass[0]  = gps_mother[1].mass()
        tree_event.lq2_pt[0]    = gps_mother[1].pt()
        dilqp4                  = gps_mother[0].p4() + gps_mother[1].p4()
        tree_event.m_lqlq[0]    = dilqp4.M()
      else:
        tree_event.lq1_mass[0]  = -1
        tree_event.lq1_pt[0]    = -1
        tree_event.lq2_mass[0]  = -1
        tree_event.lq2_pt[0]    = -1
        tree_event.m_lqlq[0]    = -1
      tree_event.ntau_assoc[0]  = len(taus_assoc)
      tree_event.ntau_decay[0]  = len(taus_decay)
      tree_event.nbgen_decay[0] = len(bgen_decay)
      tree_event.Fill()

  print ">>>   processed %4s events in %s"%(evtid,formatTime(time.time()-start_proc))
  print ">>>   writing to output file %s..."%(outfilename)
  outfile.Write()
  outfile.Close()
  print ">>>   done in in %s"%(formatTime(time.time()-start1))


root_dtype = {
  float: 'D',  int: 'I',  bool: 'O',
  'f':   'D',  'i': 'I',  '?':  'O',  'b': 'b', }

def addBranch(self, name, dtype='f', default=None):
   """Add branch with a given name, and create an array of the same name as address."""
   if hasattr(self,name):
     print "ERROR! TreeProducerCommon.addBranch: Branch of name '%s' already exists!"%(name)
     exit(1)
   if isinstance(dtype,str):
     if dtype.lower()=='f': # 'f' is only a 'float32', and 'F' is a 'complex64', which do not work for filling float branches
       dtype = float        # float is a 'float64' ('f8')
     elif dtype.lower()=='i': # 'i' is only a 'int32'
       dtype = int            # int is a 'int64' ('i8')
   setattr(self,name,num.zeros(1,dtype=dtype))
   self.Branch(name, getattr(self,name), '%s/%s'%(name,root_dtype[dtype]))
   if default!=None:
     getattr(self,name)[0] = default
TTree.addBranch = addBranch

def findCentrajpt1s(leadJets, otherJets ):
  """Finds all jets between the 2 leading jets, for central jet veto."""
  if not len(otherJets):
      return []
  etamin = leadJets[0].eta()
  etamax = leadJets[1].eta()
  def isCentral(jet):
    if jet.pt()<30.:
      return False
    eta = jet.eta()
    return etamin<eta and eta<etamax
  if etamin > etamax:
      etamin, etamax = etamax, etamin
  centrajpt1s = filter(isCentral,otherJets)
  return centrajpt1s

def p4sumvis(particles):
  #import pdb; pdb.set_trace()
  visparticles = copy.deepcopy([p for p in particles if abs(p.pdgId()) not in [12, 14, 16]])
  p4 = visparticles[-1].p4() if particles else 0.
  visparticles.pop()
  for p in visparticles:
      p4 += p.p4()
  return p4

def finalDaughters(particle, daughters):
  """Fills daughters with all the daughters of particle recursively."""
  if particle.numberOfDaughters()==0:
    daughters.append(particle)
  else:
    foundDaughter = False
    for i in range( particle.numberOfDaughters() ):
      dau = particle.daughter(i)
      if dau.status()>=1:
        daughters = finalDaughters( dau, daughters )
        foundDaughter = True
    if not foundDaughter:
      daughters.append(particle)
  return daughters

def tauDecayMode(tau):
  unstable = True
  dm = 'tau'
  final_daughter = tau
  while unstable:
    nod = tau.numberOfDaughters()
    for i in range(nod):
      dau = tau.daughter(i)
      if abs(dau.pdgId())==11 and dau.status()==1:
        dm = 'ele'
        final_daughter = dau
        unstable = False
        break
      elif abs(dau.pdgId())==13 and dau.status()==1:
        dm = 'muon'
        final_daughter = dau
        unstable = False
        break
      elif abs(dau.pdgId())==15: #taus may do bremsstrahlung
        dm = 'tau'
        final_daughter = dau
        tau = dau # check its daughters
        break
      elif abs(dau.pdgId()) not in (12, 14, 16):
        unstable = False
        break
  return dm

def diTauDecayMode(decay1,decay2):
  if decay1=='ele' and decay2=='ele':
    return 5, 'ee'
  elif decay1=='muon' and decay2=='muon':
    return 3, 'mm'
  elif (decay1=='ele' and decay2=='muon') or (decay1=='muon' and decay2=='ele'):
    return 4, 'em'
  elif (decay1=='tau' and decay2=='ele') or (decay1=='ele' and decay2=='tau'):
    return 2, 'et'
  elif (decay1=='tau' and decay2=='muon') or (decay1=='muon' and decay2=='tau'):
    return 1, 'mt'
  elif decay1=='tau' and decay2=='tau':
    return 0, 'tt'
  else:
    return -1, 'anything'

def deltaR(eta1, phi1, eta2, phi2):
  deta = eta1 - eta2
  dphi = deltaPhi(phi1, phi2)
  return math.sqrt(deta*deta+dphi*dphi)

def deltaPhi(phi1, phi2):
  """Computes delta phi, handling periodic limit conditions."""
  res = phi1 - phi2
  while res > math.pi:
    res -= 2*math.pi
  while res<-math.pi:
    res += 2*math.pi
  return res

def cleanObjectCollection(objects, masks, dRmin):
  """Masks objects using a deltaR cut."""
  if len(objects)==0 or len(masks)==0:
    return objects, []
  cleanObjects = [ ]
  dirtyObjects = [ ]
  for object in objects:
    overlap = True
    for mask in masks:
        dR = deltaR(object.eta(),object.phi(),mask.eta(),mask.phi())
        if dR<dRmin:
          overlap = False
          break
    if overlap:
      dirtyObjects.append(object)
    else:
      cleanObjects.append(object)
  return cleanObjects, dirtyObjects

def isFinal(p):
  # TODO: check if one daughter is final and has same PID
  return not (p.numberOfDaughters()==1 and p.daughter(0).pdgId()==p.pdgId())

def isFinalM(p):
  #print p.numberOfDaughters(), p.daughter(0).pdgId(), p.pdgId(), p.status()
  #return p.numberOfDaughters()==3
  return not (p.numberOfDaughters()==3 and p.daughter(0).pdgId()==p.pdgId())

def printParticle(p):
  string = "%9d: status=%2d, pt=%7.2f, eta=%5.2f, phi=%5.2f, final=%5s"%(p.pdgId(),p.status(),p.pt(),p.eta(),p.phi(),isFinal(p))
  if p.numberOfMothers()>=2:
    string += ", mothers %s, %s"%(p.mother(0).pdgId(),p.mother(1).pdgId())
  elif p.numberOfMothers()==1:
    string += ", mother %s"%(p.mother(0).pdgId())
  if p.numberOfDaughters()>=2:
    string += ", daughters %s, %s"%(p.daughter(0).pdgId(),p.daughter(1).pdgId())
  elif p.numberOfDaughters()==1:
    string += ", daughter %s"%(p.daughter(0).pdgId())
  print string

def formatTime(seconds):
  minutes, seconds = divmod(seconds, 60)
  hours,   minutes = divmod(minutes, 60)
  if   hours:   return "%d hours, %d minutes and %.3g seconds"%(hours,minutes,seconds)
  elif minutes: return "%d minutes and %.3g seconds"%(minutes,seconds)
  return "%.3g seconds"%(seconds)

def formatTimeShort(seconds):
  minutes, seconds = divmod(seconds, 60)
  hours,   minutes = divmod(minutes, 60)
  return "%02d:%02d:%02d"%(hours,minutes,seconds)

def bold(string):
  return "\033[1m%s\033[0m"%(string)

def ETA(start,iJob,nJobs):
  return formatTimeShort((time.time()-start)*(nJobs-iJob)/iJob)

def stepsize(total):
  return min(max(10**(ceil(log(total,10))-1),20),100)

def ensureDir(dir):
  if not os.path.exists(dir):
    os.makedirs(dir)
    print ">>> made directory %s\n>>>"%(dir)
  return dir


if __name__=='__main__':
  print ">>> "
  main()
  print ">>> done\n"
