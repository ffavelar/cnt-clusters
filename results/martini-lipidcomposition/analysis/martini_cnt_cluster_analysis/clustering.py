#!/usr/bin/env python

#-------------------------#
# CNT CLUSTERING ANALYSIS #
#-------------------------#
#
# Reads the COM of all carbon nanotubes and writes to a file the number of clusters as well as the size of the biggest one
#
# Martin Voegele, MPI of Biophysics, 2016-09-18
#



#-----------------#
# IMPORT PACKAGES #
#-----------------#


import sys
import argparse
import numpy as np
import scipy as sp
import MDAnalysis as mda
import time




#----------------------#
# FUNCTION DEFINITIONS #
#----------------------#


def distance(x0, x1, dimensions):
    """Calculate the distance concerning periodic boundary conditions."""
    
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    
    return np.sqrt((delta ** 2).sum(axis=-1))


def connection_list(coord, dim, cutoff):
    """Create a connection matrix for the coordinates based on the cutoff radius"""
    
    n = len(coord)
    
    connections = []
    
    for i in xrange(n): 
        newconnect = []        
        for j in xrange(n):
            dist = distance(coord[i,:3],coord[j,:3],dim[:3])
            if dist <= cutoff:
                newconnect.append(j)
        connections.append(newconnect)
    
    return connections
    

def clustering(coord, dim, cutoff):
    """Calculate the number of clusters and the size of the biggest cluster."""
    
    # Create the connection list
    connections          = connection_list(coord,dim,cutoff)

    # Initialize the counting variables
    n                    = len(coord)
    number_of_clusters   = 0
    biggest_cluster_size = 0 
    cluster_index        = np.zeros(n, dtype=np.int)
    cluster_size         = []
    
    # Count the clusters and number of elements per cluster
    for i in xrange(n): 
        if cluster_index[i] == 0:
            number_of_clusters += 1
            cluster_index[i] = number_of_clusters
            cluster_size.append(1)
            queue = connections[i]
            while len(queue) > 0:
                newqueue = []
                for j in queue:
                    if cluster_index[j] == 0:
                        cluster_index[j] = number_of_clusters
                        cluster_size[cluster_index[i]-1] += 1
                        for k in connections[j]:
                            newqueue.append(k)
                queue = newqueue         
                
    biggest_cluster_size = np.max(cluster_size)
            
    return number_of_clusters, biggest_cluster_size




#----------------#
#  MAIN ROUTINE  #
#----------------#


start = time.time()


# PARSE THE ARGUMENTS

parser = argparse.ArgumentParser(description='Reads the COM of all Carbon nanotubes and writes to a file the number of clusters as well as the size of the biggest one.')
parser.add_argument( '-c', '--cutoff',  type=float, default=22.5,             help='cut-off radius for clustering analysis' )
parser.add_argument( '-b', '--begin',   type=int,   default=0,                help='frame to begin with'                    )
parser.add_argument( '-e', '--end',     type=int,   default=None,             help='last frame'                             )
parser.add_argument( '-s', '--skip',    type=int,   default=1,                help='use only every sth frame'               )
parser.add_argument( '-t', '--topol',   type=str,   default='martini_md.tpr', help='topology file'                          )
parser.add_argument( '-x', '--traj',    type=str,   default='martini_md.xtc', help='trajectory file'                        )
parser.add_argument( '-o', '--outfile', type=str,   default='clustering.dat', help='output file'                            ) 
args = parser.parse_args()


# LOAD NECESSARY DATA

universe     = mda.Universe(args.topol,args.traj)
cnt_atoms    = universe.select_atoms("resname CNT")
cnt_residues = cnt_atoms.residues


# GO THROUGH THE TRAJECTORY AND DO ANALYSIS

with open(args.outfile, 'w') as output:

    # If no end frame is given, the last frame will be the end
    length = len(universe.trajectory)
    if args.end==None:
        end = length
        
    # Go through the trajectory
    for ts in universe.trajectory[args.begin:end:args.skip]:
        
        # Timestep properties
        print "Frame %4d" % ts.frame
        t   = ts.time # in ps    
        dim = ts.dimensions[0:3]
        
        # Make a list of centers of mass of all CNTs
        com_list = []
        for tube in cnt_residues:
            com_list.append(tube.center_of_mass())
        com_list = np.array(com_list)
        
        # Calculate the number of clusters and the size of the biggest one
        num, size = clustering(com_list, dim, args.cutoff)
    
        # Write data to file
        output.write("%12.3f \t %8d \t %8d\n" % (t, num, size))


# Calculate the elapsed time
end  = time.time()
elap = end - start
print 'Run time:', elap, 's'
