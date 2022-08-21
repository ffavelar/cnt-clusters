#!/usr/bin/env python

# Reads the COM of all Carbon naotubes and writes to a file the number of clusters as well as the size of the biggest one
# Martin Voegele, MPI of Biophysics, 2015-08-18



# IMPORT PACKAGES

import sys
import MDAnalysis
import numpy
import scipy
import argparse



# PARSE THE ARGUMENTS

parser = argparse.ArgumentParser(description='Reads the average angle (in degrees) with respect to the z axis of all carbon nanotubes and writes it to a file.')
parser.add_argument('-s', '--skip', type=int, default=1, help='use only every sth frame')
parser.add_argument('-t', '--topol', type=str, default='martini_md.tpr', help='topology file')
parser.add_argument('-x', '--traj', type=str, default='martini_md.xtc', help='trajectory file')
parser.add_argument('-o', '--outfile', type=str, default='angles.dat', help='output file')

args = parser.parse_args()



# LOAD NECESSARY DATA

universe  = MDAnalysis.Universe(args.topol,args.traj)
cntatoms  = universe.selectAtoms("resname CNT")
cntresids = cntatoms.resids()



# FUNCTION DEFINITIONS


def distance(x0, x1, dimensions):
    """Calculate the distance concerning periodic boundary conditions."""
    delta = numpy.abs(x0 - x1)
    delta = numpy.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return numpy.sqrt((delta ** 2).sum(axis=-1))


   
# MAIN ROUTINE

with open(args.outfile, 'w') as output:

    for ts in universe.trajectory[::args.skip]:

        print "Frame %4d" % ts.frame
        time = ts.time # in ps
    
        allAngles = numpy.zeros(len(cntresids))
        allTiltcos = numpy.zeros(len(cntresids))        

        for i,resid in enumerate(cntresids):
            ring1 = universe.selectAtoms("resid "+str(resid)+" and name F00*")
            ring2 = universe.selectAtoms("resid "+str(resid)+" and (name F08* or name F09*)")
            vector = ring2.centerOfMass() - ring1.centerOfMass()
            tiltcos = vector[2]/numpy.linalg.norm(vector)
            angle = numpy.arccos(tiltcos)
            allAngles[i] = angle*180/numpy.pi
            allTiltcos[i] = tiltcos
            
        averAngle = numpy.average(allAngles)
        averTiltcos = numpy.average(allTiltcos)        

        output.write("%12.3f \t %8d \t %8d \n" % (time, averAngle, averTiltcos))

