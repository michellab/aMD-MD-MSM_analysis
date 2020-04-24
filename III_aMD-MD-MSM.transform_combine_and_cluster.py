#!/usr/bin/python2.7

#   =================== aMD/MD/MSM Analysis STEP III v1.0 ===================   #
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


import pyemma
import pyemma.coordinates as coor
import os, sys, pickle
import MDAnalysis as mda
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import mdtraj as md
import numpy as np

__author__='Jordi Juarez-Jimenez, PhD'

if __name__ == '__main__':
        

    #=== CYPA ===#
    print "Processing CypA..."

    picfile = open ('Intermediate_pickle_files/cypa_cattraj_rmsd_array.pickle', 'r')
    wt = pickle.load(picfile)
    picfile.close()
    picfile = open ('Intermediate_pickle_files/cypa_cattraj_dist4id.pickle', 'r')
    wt_D = pickle.load(picfile)
    picfile.close()

    trans_wt = []
    tres100s = 14.0
    for traj in wt_D:
        temp = []
        for snapshot in traj:
            if snapshot >= tres100s:
                temp.append(1)
            else:
                temp.append(-1)
        trans_wt.append(temp)
    np.shape(trans_wt)

    wt_dir_rmsd = []
    for i in range(len(wt)):
        traj = []
        for j in range(len(wt[i])):
            new_rmsd_70s = wt[i][j][1]*1 # No modification of the 70s rmsd
            new_rmsd_100s = wt[i][j][2]*trans_wt[i][j] # therefore each snapshot in trans contains a single value
            traj.append([new_rmsd_70s, new_rmsd_100s])
        wt_dir_rmsd.append(traj)

    print "Done!"
    
    #=== D66A ===#
    print "Processing D66A..."

    picfile = open ('Intermediate_pickle_files/d66a_cattraj_rmsd_array.pickle', 'r')
    d66a = pickle.load(picfile)
    picfile.close()
    picfile = open ('Intermediate_pickle_files/d66a_cattraj_dist4id.pickle', 'r')
    d66a_D = pickle.load(picfile)
    picfile.close()

    trans_d66a = []
    tres100s = 14.0
    for traj in d66a_D:
        temp = []
        for snapshot in traj:
            if snapshot >= tres100s:
                temp.append(1)
            else:
                temp.append(-1)
        trans_d66a.append(temp)
    np.shape(trans_d66a)

    d66a_dir_rmsd = []
    for i in range(len(d66a)):
        traj = []
        for j in range(len(d66a[i])):
            new_rmsd_70s = d66a[i][j][1]*1 # No modification of the 70s rmsd
            new_rmsd_100s = d66a[i][j][2]*trans_d66a[i][j] # therefore each snapshot in trans contains a single value
            traj.append([new_rmsd_70s, new_rmsd_100s])
        d66a_dir_rmsd.append(traj)

    print "Done!"

    #=== H70A ===#
    print "Processing H70A..."
    
    picfile = open ('Intermediate_pickle_files/h70a_cattraj_rmsd_array.pickle', 'r')
    h70a = pickle.load(picfile)
    picfile.close()
    picfile = open ('Intermediate_pickle_files/h70a_cattraj_dist4id.pickle', 'r')
    h70a_D = pickle.load(picfile)
    picfile.close()

    trans_h70a = []
    tres100s = 14.0
    for traj in h70a_D:
        temp = []
        for snapshot in traj:
            if snapshot >= tres100s:
                temp.append(1)
            else:
                temp.append(-1)
        trans_h70a.append(temp)
    np.shape(trans_h70a)

    h70a_dir_rmsd = []
    for i in range(len(h70a)):
        traj = []
        for j in range(len(h70a[i])):
            new_rmsd_70s = h70a[i][j][1]*1 # No modification of the 70s rmsd
            new_rmsd_100s = h70a[i][j][2]*trans_h70a[i][j] # therefore each snapshot in trans contains a single value
            traj.append([new_rmsd_70s, new_rmsd_100s])
        h70a_dir_rmsd.append(traj)

    print "Done!"


    print "Combining dir RMSD data..."

    combined = np.concatenate((wt_dir_rmsd,d66a_dir_rmsd, h70a_dir_rmsd), axis=0)


    # In[20]:

    histodataX=[]
    histodataY=[]
    for i in range(len(combined)):
        for p in range(len(combined[i])):
            histodataX.append(combined[i][p][0])
            histodataY.append(combined[i][p][1])
            
    z,x,y = np.histogram2d(histodataX,histodataY,bins=200)
    F = -np.log(z)

    zfile = open('Intermediate_pickle_files/wt-h70a-d66a_cattraj_contourmap.pickle', 'w')
    pickle.dump(F,zfile)
    zfile.close()

    test =[]
    for t in combined:
        test.append(t)
    np.shape(test)

    # ###  100 K-means clusters
    nclusters = 100

    kmean_cluster100=coor.cluster_kmeans(data=test, k=nclusters, max_iter=1000, tolerance=1e-6)

    print "Done!"
    print "Saving cluster centers..."
    ccenters100 = kmean_cluster100.clustercenters
    f=open('Intermediate_pickle_files/wt-h70a-d66a_cattraj_dirrmsd_ccenter-100.pickle', 'w')
    pickle.dump(ccenters100,f)
    f.close()
    
    wt_dtrajs=coor.assign_to_centers(data=wt_dir_rmsd, centers=ccenters100)
    f=open('Intermediate_pickle_files/cypa_wt-d66a_cattraj_dirrmsd_dtrajs.pickle', 'w')
    pickle.dump(wt_dtrajs,f)
    f.close()
    
    d66a_dtrajs=coor.assign_to_centers(data=d66a_dir_rmsd, centers=ccenters100)
    f=open('Intermediate_pickle_files/d66a_cattraj_dirrmsd_dtrajs.pickle', 'w')
    pickle.dump(wt_dtrajs,f)
    f.close()
    
    h70a_dtrajs=coor.assign_to_centers(data=h70a_dir_rmsd, centers=ccenters100)
    f=open('Intermediate_pickle_files/h70a_cattraj_dirrmsd_dtrajs.pickle', 'w')
    pickle.dump(wt_dtrajs,f)
    f.close()
    
    
    print "Done!"
    print "Completed the third aMD/MD/MSM analysis step"
    

    print "Done!"
    
