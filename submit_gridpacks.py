#! /usr/bin/env python

import os, sys
import subprocess
import time

def is_file_empty(file_path):
    """ Check if file is empty by confirming if its size is 0 bytes"""
    # Check if file exist and it is empty
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0


make_cards = False
submit_gridpacks = True
copy_gridpacks = False
clean_mg_area = False


# masses_scalar = [800,  1100, 1500]
masses_scalar = [1500]
masses_vector = [1000, 1500, 2000]
# masses_vector = [1000]
# lambdas = [0.5, 1.0, 1.5, 2.0, 2.5]
lambdas = [2.5]
kappas = [1.0, 0.0]
# kappas = [1.0]
# productions = ['Single', 'Pair', 'NonRes']
productions = ['Single']
# spins = ['Scalar', 'Vector']
spins = ['Scalar']
# masses = [1100]
# lambdas = [1.0]
# tag = 'LOPDF31'
tag = 'LOPDF31_FIRSTGEN'


# Make cards first, if needed
for prod in productions:
    for spin in spins:
        masses = masses_scalar if spin=='Scalar' else masses_vector
        massstring = ' '.join([str(m) for m in masses])
        lambdastring = ','.join([str(l) for l in lambdas])
        kappastring = ','.join([str(k) for k in kappas])
        command = "./create_cards.py cards/%sLQ_%s/%sLQ_%s_template_*.dat -n 'M$MASS_L$LAMBDA_%s' -m " % (spin, prod, spin, prod, tag)
        if spin == 'Vector':
            command = command.replace("-n 'M$MASS_L$LAMBDA_%s'" % (tag), "-n 'M$MASS_L$LAMBDA_K$KAPPA_%s'" % (tag))
        command += massstring + " -p LAMBDA=" + lambdastring
        if spin == 'Vector':
            command+= ':KAPPA=' + kappastring
        print command
        if make_cards: os.system(command)



# Submit gridpacks based on cards created above
for prod in productions:
    for spin in spins:
        masses = masses_scalar if spin=='Scalar' else masses_vector
        base_command = './gridpack_generation.sh %sLQ_%s_MMYMASS_LMYLAMBDA_%s ../../../../../GENSIM/cards/%sLQ_%s local' % (spin, prod, tag, spin, prod)
        if spin == 'Vector':
            base_command = base_command.replace("_MMYMASS_LMYLAMBDA_%s" % (tag), "_MMYMASS_LMYLAMBDA_KMYKAPPA_%s" % (tag))

        commands = []
        idx = 1
        for mass in masses:
            for lamb in lambdas:
                command = base_command
                command = command.replace('MYMASS', str(mass))
                command = command.replace('MYLAMBDA', str(lamb).replace('.', 'p'))
                if spin == 'Vector':
                    for kap in kappas:
                        thiscommand = command
                        thiscommand = thiscommand.replace('MYKAPPA', str(kap).replace('.', 'p'))
                        commands.append(thiscommand)
                        print thiscommand
                        idx += 1
                else:
                    commands.append(command)
                    print command
                    idx += 1

        f = open('commands/commands_submit_gridpacks_%s.txt'%(prod), 'w')
        for comm in commands:
            f.write(comm + '\n')
        f.close()

        njobs = len(commands)
        print njobs
        jobname = spin + "LQ_" + prod
        command = 'sbatch -a 1-%s --mem 4000 -J gridpacks_%s submit_gridpacks.sh commands/commands_submit_gridpacks_%s.txt' % (str(njobs), jobname, prod)
        # command = 'sbatch -a 1-%s -J gridpacks_%s submit_gridpacks.sh commands/commands_submit_gridpacks_%s.txt' % (str(njobs), jobname, prod)
        print command
        if submit_gridpacks: os.system(command)

# Check job status
if (copy_gridpacks or clean_mg_area) and submit_gridpacks:
    while True:
        fout = open('.logfile_slurmsubmission.out', 'w')
        ferr = open('.logfile_slurmsubmission.err', 'w')
        p = subprocess.Popen("squeue -u areimers | grep gridpack", stdout=fout, stderr=ferr, shell=True)
        fout.close()
        ferr.close()

        time.sleep(30)
        if is_file_empty('.logfile_slurmsubmission.out'):
            print 'No gridpack jobs running anymore.'
            break



# Copy gridpacks to new dir
gpnames = []
for spin in spins:
    masses = masses_scalar if spin=='Scalar' else masses_vector
    for prod in productions:
        time.sleep(5)
        mgfolder = '/work/areimers/CMSSW_10_2_10/src/genproductions/bin/MadGraph5_aMCatNLO'
        gpname_template = '%sLQ_%s_MMYMASS_LMYLAMBDA_%s_slc7_amd64_gcc700_CMSSW_9_3_16_tarball.tar.xz' % (spin, prod, tag)
        if spin == 'Vector':
            gpname_template = gpname_template.replace('_MMYMASS_LMYLAMBDA', '_MMYMASS_LMYLAMBDA_KMYKAPPA')
        for mass in masses:
            for lamb in lambdas:
                gpname = gpname_template
                gpname = gpname.replace('MYMASS', str(mass))
                gpname = gpname.replace('MYLAMBDA', str(lamb).replace('.', 'p'))

                if spin == 'Vector':
                    for kap in kappas:
                        thisgpname = gpname
                        thisgpname = thisgpname.replace('MYKAPPA', str(kap).replace('.', 'p'))
                        gpnames.append(thisgpname)
                        print thisgpname
                else:
                    gpnames.append(gpname)
                    print gpname

if copy_gridpacks:
    procs = []
    for gp in gpnames:
        command = 'cp ' + mgfolder + '/' + gp + ' gridpacks'
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        procs.append(p)

    for p in procs:
        p.wait()

    print 'Done moving gridpacks.'

# Clean up mg folder from remnants of gp generation
if clean_mg_area:
    n_running = 0
    n_completed = 0
    n_jobs = len(gpnames)
    processes = []
    for gp in gpnames:
        b_wait = (n_running >= 15)
        while b_wait:
            n_running = 0
            n_completed = 0
            for proc in processes:
                if proc.poll() == None: n_running += 1
                else:
                    n_completed += 1
            percent = round(float(n_completed)/float(n_jobs)*100, 1)
            sys.stdout.write( '{0:d} of {1:d} ({2:4.2f}%) jobs done.\r'.format(n_completed, n_jobs, percent))
            sys.stdout.flush()
            time.sleep(10)
            b_wait = (n_running >= 15)

        n_running += 1
        gp = gp.replace('_slc7_amd64_gcc700_CMSSW_9_3_16_tarball.tar.xz', '')
        command = 'rm -rf ' + mgfolder + '/' + gp + '*'
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        processes.append(p)

    b_wait = (n_completed < n_jobs)
    while b_wait:
        n_running = 0
        n_completed = 0
        for proc in processes:
            if proc.poll() == None: n_running += 1
            else:
                n_completed += 1
        percent = float(n_completed)/float(n_jobs)*100
        sys.stdout.write( '{0:d} of {1:d} ({2:4.2f} %) jobs done.\r'.format(n_completed, n_jobs, percent))
        sys.stdout.flush()
        time.sleep(10)
        b_wait = (n_completed < n_jobs)
    print "Done cleaning up MadGraph folder."
