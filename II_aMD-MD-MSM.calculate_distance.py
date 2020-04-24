#!/usr/bin/python2.7

#   =================== aMD/MD/MSM Analysis STEP II v1.0 ===================    #
#    Copyright (C) 2016-2018  Jordi Juarez-Jimenez, PhD                         #
#    Michel lab, The University of Edinburgh                                    #
#                                                                               #
#    This program is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by       #
#    the Free Software Foundation, either version 3 of the License, or          #
#    (at your option) any later version.                                        #
#                                                                               #
#    This program is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#    GNU General Public License for more details.                               #
#                                                                               #
#    You should have received a copy of the GNU General Public License          #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
#################################################################################

# This script was used to obtain some of the results reported in:
# Jordi Juárez-Jiménez*, Arun A. Gupta, Gogulan Karunanithy, Antonia S. J. S. Mey, 
# Charis Georgiou, Harris Ioannidis, Alessio De Simone, Paul N. Barlow, Alison N. Hulme,
# Malcolm D. Walkinshaw, Andrew J. Baldwin and  Julien Michel*
# "Dynamic design: manipulation of millisecond timescale motions on the energy landscape 
# of cyclophilin A" , Chem. Sci., 2020, 11, 2670
# Please cite that work if you find this script useful.

import pickle,os, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis import *
from MDAnalysis.analysis.distances import *
import time
import pickle

__author__='Jordi Juarez-Jimenez, PhD'

def module_vector(p1,p2):
    """
    by Marina Corbella
    University of Barcelona
    """

    v1 = [0,0,0]

    v1[0] = p2[0] - p1[0]

    v1[1] = p2[1] - p1[1]

    v1[2] = p2[2] - p1[2]

    module = np.sqrt((v1[0]**2) + (v1[1]**2) + (v1[2]**2))

    return module

def calculate_100s2core_dist(topfile, trajfiles):
    '''

    '''

    mda_dist = []
    
    j = 0.
    control = np.arange(0,len(trajfiles),round(len(trajfiles)*0.1))
    print "Starting mdanalysis block..."
    for inputtraj in trajfiles:
        if j in control:
            print "%.2f" % (j/len(trajfiles)*100) + "% " +"completed"
            sys.stdout.flush()
        MD = mda.Universe(topfile,inputtraj)
        A= MD.select_atoms('resid 102 108 and name CA ')
        B= MD.select_atoms('resid 24 and name CA or resid 30 and name CA or resid 129 and name CA')
        D = []
        for ts in MD.trajectory:
            com1=A.center_of_mass()
            com2=B.center_of_mass()
            d = module_vector(com1,com2)
            D.append(d)
        D = np.array(D)
        mda_dist.append(D)
        j = j + 1.

    return mda_dist

    


if __name__ == '__main__':

    
    #=== CYPA ===#
    # Load the topology file as a gromacs gro file
    myrootpath = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/'
    topfile = myrootpath + 'cypa.dry.gro'
    # Load the pickle containing the list of and path to trajectories
    print 'Loading the trajectory paths for CypA...'
    
    f=open('Intermediate_pickle_files/cypa_traj_cattraj.pickle', 'r')
    trajfiles = pickle.load(f)
    f.close()

    print 'Done!'
    print 'Total number of trajectories: %d' % len(trajfiles)
    
    print ' Calculating 100s loop distance to the protein core...'

    mda_dist = calculate_100s2core_dist(topfile, trajfiles)

    f=open('Intermediate_pickle_files/cypa_cattraj_dist4id.pickle', 'w')
    pickle.dump(mda_dist,f)
    f.close()
    print "Done!"
    #====================================================#

    #=== D66A ===#
    # Load the topology file as a gromacs gro file
    myrootpath = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/'
    topfile = myrootpath + 'd66a.dry.gro'
    
    # Load the pickle containing the list of and path to trajectories
    print 'Loading the trajectory paths for D66A...'
    
    f=open('Intermediate_pickle_files/d66a_traj_cattraj.pickle', 'r')
    trajfiles = pickle.load(f)
    f.close()

    print 'Done!'
    print 'Total number of trajectories: %d' % len(trajfiles)
    
    print ' Calculating 100s loop distance to the protein core...'

    mda_dist = calculate_100s2core_dist(topfile, trajfiles)

    f=open('Intermediate_pickle_files/d66a_cattraj_dist4id.pickle', 'w')
    pickle.dump(mda_dist,f)
    f.close()
    print "Done!"
    #====================================================#

    #=== H70A ===#
    # Load the topology file as a gromacs gro file
    myrootpath = '/home/jjuarez/Work/UEDIN/ISOTRAPSS/CONF_ANALYSIS/CYPA/AMD_PROD_RUNS/dry_grofiles/'
    topfile = myrootpath + 'h70a.dry.gro'
    
    # Load the pickle containing the list of and path to trajectories
    print 'Loading the trajectory paths for H70A...'
    
    f=open('Intermediate_pickle_files/h70a_traj_cattraj.pickle', 'r')
    trajfiles = pickle.load(f)
    f.close()

    print 'Done!'
    print 'Total number of trajectories: %d' % len(trajfiles)
    
    print ' Calculating 100s loop distance to the protein core...'

    mda_dist = calculate_100s2core_dist(topfile, trajfiles)

    f=open('Intermediate_pickle_files/h70a_cattraj_dist4id.pickle', 'w')
    pickle.dump(mda_dist,f)
    f.close()
    print "Done!"
    #====================================================#
    
    print "Completed the second aMD/MD/MSM analysis step"
