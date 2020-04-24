#!/usr/bin/python2.7

#   =================== aMD/MD/MSM Analysis STEP I v1.0 ===================	#
#    Copyright (C) 2016-2018  Jordi Juarez-Jimenez, PhD				#
#    Michel lab, The University of Edinburgh.					#
#										#
#    This program is free software: you can redistribute it and/or modify	#
#    it under the terms of the GNU General Public License as published by	#
#    the Free Software Foundation, either version 3 of the License, or		#
#    (at your option) any later version.				        #
#										#
#    This program is distributed in the hope that it will be useful, 		#
#    but WITHOUT ANY WARRANTY; without even the implied warranty of		#
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		#
#    GNU General Public License for more details.				#
#										#
#    You should have received a copy of the GNU General Public License		#
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.	#
#################################################################################

# This script was used to obtain some of the results reported in:
# Jordi Juárez-Jiménez*, Arun A. Gupta, Gogulan Karunanithy, Antonia S. J. S. Mey, 
# Charis Georgiou, Harris Ioannidis, Alessio De Simone, Paul N. Barlow, Alison N. Hulme,
# Malcolm D. Walkinshaw, Andrew J. Baldwin and  Julien Michel*
# "Dynamic design: manipulation of millisecond timescale motions on the energy landscape 
# of cyclophilin A" , Chem. Sci., 2020, 11, 2670
# Please cite that work if you find this script useful.

import os, sys, pickle
import MDAnalysis as mda
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import mdtraj as md
import numpy as np


def gather_trajfiles(myrootpath,longtrajs, seeded_folders, folder_names,trajname, sim_range):
    '''(str, list, list,list,list, list) -> list

    Gathers all the paths to trajectory files in longtrajs and seeded_folders and returns
    them in a list

    '''
    trajfiles = []
    print 'Adding long trajectories'
    for traj in longtrajs:
        trajfiles.append(traj)

    first = sim_range[0]
    last = sim_range[1]
    print 'Adding seeded trajectories...'
    for seed in seeded_folders:
        for i in range(first,last):
            for sysname in folder_names:
                for o in range(len(trajname)):
                    filename = myrootpath+seed+'/'+sysname+'_st%d/'%i+sysname+trajname[o]
                    if os.path.isfile(filename): 
                        trajfiles.append(filename)

    return trajfiles

def calculate_loops_RMSD(topfile, refile, trajfiles):
    '''

    Calculate the RMSD of the 65-75 region and 100-110 region after removing translation/rotation movements
    by fitting the more rigid protein regions.

    '''

    ref=mda.Universe(refile)
    Y = []
    for inputraj in trajfiles:
        MD=mda.Universe(topfile,inputraj)
        R = mda.analysis.rms.RMSD(MD,ref, select='resid 5-60 and name CA or resid 80-90 and name CA or resid 120-160 and name CA', 
                            groupselections=[" resid 5-60 and name CA or resid 80-90 and name CA or resid 120-160 and name CA",   # CORE
                                   "name CA and resid 65-75",                                   # 70s
                                   "name CA and resid 100-110"], # 100s
                                    filename="__temp_rmsd_.dat")
        R.run()
        Y.append(R.rmsd[:,[3,4,5]])

    return Y
    
    

##

__author__='Jordi Juarez-Jimenez, PhD'

if __name__ == '__main__':
    
    # This section should be completed with paths towards different trajectories.
    # The objective is to pupulate the trajfiles list with all the trajectories

    #=== CYPA ===# 
    print 'Loading trajectories for CypA...'
    myrootpath = '/media/data4/jjuarez/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/DRY_CAT_TRAJ/'
    refile='/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/SYS_PREP/1ak4/1ak4.pdb'
    topfile = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/cypa.dry.gro'
    trajname = ['.cattraj.xtc']
    seeded_folders = ['100s_seeding','100s_seeding_2w', '100s_seeding_3w', '70s_seeding', '70s_seeding_2w', '70s_seeding_3w','100s-2_seeding','100s-2_seeding_2w', '100s-2_seeding_3w']
    folder_names = ['cypa_100s', 'cypa_70s', 'cypa_100s-2']
    # Long trajectories
    longtraj1 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/cypa.cattraj.xtc'
    longtraj2 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/cypa_100s_st86.cattraj.xtc'
    longtraj3 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/cypa_70s.cattraj.xtc'

    longtrajs = [longtraj1, longtraj2, longtraj3]
    
    # This object will hold the paths to all the trajectories keeping the order for subsequent analysis.
    # i.e: trajfile[0] will  point to '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/cypa.cattraj.xtc'
    # in this and subsequent scripts

    trajfiles = gather_trajfiles(myrootpath, longtrajs, seeded_folders, folder_names,trajname, [1,200])
    
    print 'Done!'
    print 'Total number of trajectories for CypA: %d' % len(trajfiles)

    # The cypa_traj_cattraj.pickle will be used in subsequent scripts to track the trajectories

    f=open('Intermediate_pickle_files/cypa_traj_cattraj.pickle', 'w')
    pickle.dump(trajfiles, f)
    f.close()

    print 'Calculating RMSD...'
    sys.stdout.flush()
    Y = calculate_loops_RMSD(topfile, refile, trajfiles)

    # Save the pickle containing the rmsd values
    picfile = open ( 'Intermediate_pickle_files/cypa_cattraj_rmsd_array.pickle', 'w')
    pickle.dump(Y,picfile)
    picfile.close()

    print "Done!"
    #====================================================#

    #=== D66A ===#
    print 'Loading trajectories for D66A...'
    myrootpath = '/media/data4/jjuarez/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/DRY_CAT_TRAJ/'
    refile='/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/SYS_PREP/1ak4/1ak4.pdb'
    topfile = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/d66a.dry.gro'
    trajname = ['.cattraj.xtc']
    seeded_folders = ['d66a_100s_seeding','d66a_100s_seeding_2w', 'd66a_100s_seeding_3w', 'd66a_70s_seeding', 'd66a_70s_seeding_2w', 'd66a_70s_seeding_3w','d66a_100s-2_seeding','d66a_100s-2_seeding_2w', 'd66a_100s-2_seeding_3w']
    folder_names = ['d66a_100s', 'd66a_70s', 'd66a_100s-2']
    #===#
    longtraj1 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/d66a.cattraj.xtc'
    longtraj2 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/d66a_100s.cattraj.xtc'
    longtraj3 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/d66a_70s.cattraj.xtc'

    longtrajs = [longtraj1, longtraj2, longtraj3]    
    trajfiles = gather_trajfiles(myrootpath, longtrajs, seeded_folders, folder_names,trajname, [1,200])

    print 'Done!'
    print 'Total number of trajectories for D66A: %d' % len(trajfiles)
    
    f=open('Intermediate_pickle_files/d66a_traj_cattraj.pickle', 'w')
    pickle.dump(trajfiles, f)
    f.close()

    print 'Calculating RMSD...'
    sys.stdout.flush()
    Y = calculate_loops_RMSD(topfile, refile, trajfiles)

    # Save the pickle containing the rmsd values
    picfile = open ( 'Intermediate_pickle_files/d66a_cattraj_rmsd_array.pickle', 'w')
    pickle.dump(Y,picfile)
    picfile.close()

    print "Done!"
    #====================================================#

    #=== H70A ===#
    myrootpath = '/media/data4/jjuarez/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/DRY_CAT_TRAJ/'
    refile='/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/SYS_PREP/1ak4/1ak4.pdb'
    topfile = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/h70a.dry.gro'
    trajname = ['.cattraj.xtc']
    seeded_folders = ['h70a_100s_seeding','h70a_100s_seeding_2w', 'h70a_100s_seeding_3w', 'h70a_70s_seeding', 'h70a_70s_seeding_2w', 'h70a_70s_seeding_3w','h70a_100s-2_seeding','h70a_100s-2_seeding_2w', 'h70a_100s-2_seeding_3w']
    folder_names = ['h70a_100s', 'h70a_70s', 'h70a_100s-2']
    trajfiles = []
    #===#
    longtraj1 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/h70a.cattraj.xtc'
    longtraj2 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/h70a_100s.cattraj.xtc'
    longtraj3 = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/long_traj/h70a_70s.cattraj.xtc'
    
    longtrajs = [longtraj1, longtraj2, longtraj3]
    trajfiles = gather_trajfiles(myrootpath, longtrajs, seeded_folders, folder_names,trajname, [1,200])

    print 'Done!'
    print 'Total number of trajectories H70A: %d' % len(trajfiles)
    f=open('Intermediate_pickle_files/h70a_traj_cattraj.pickle', 'w')
    pickle.dump(trajfiles, f)
    f.close()

    print 'Calculating RMSD...'
    sys.stdout.flush()
    Y = calculate_loops_RMSD(topfile, refile, trajfiles)

    # Save the pickle containing the rmsd values
    picfile = open ( 'Intermediate_pickle_files/h70a_cattraj_rmsd_array.pickle', 'w')
    pickle.dump(Y,picfile)
    picfile.close()

    print 'Done!'
    print "Completed the first aMD/MD/MSM analysis step"

