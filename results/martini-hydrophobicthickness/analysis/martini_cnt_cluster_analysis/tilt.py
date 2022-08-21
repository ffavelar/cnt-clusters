#!/usr/bin/env python

# Calculates the cosine of the tilt angle of all Carbon naotubes and writes it to a file
# Martin Voegele, MPI of Biophysics, 2016-09-12



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
parser.add_argument('-rs', '--ringsize', type=int, default=8, help='number of CG beads in one ring')


args = parser.parse_args()



# LOAD NECESSARY DATA

universe  = MDAnalysis.Universe(args.topol,args.traj)
cntatoms  = universe.select_atoms("resname CNT")
cntresids = cntatoms.resids


   
# MAIN ROUTINE

with open(args.outfile, 'w') as output:

    # go through the trajectory
    for ts in universe.trajectory[::args.skip]:
        
        # Get number, time and box size of the current time step
        print "Frame %4d" % ts.frame
        time = ts.time # in ps
        dimensions = ts.dimensions[:3]
        
	# initialize the array for the tilt cosine
        allTiltcos = numpy.zeros(len(cntresids))

	# Go through all residues
        for i,resid in enumerate(cntresids):
            
            # The main axis of the CNT is calculated as the vector between the COM of the uppermost and the lowermost ring
            ring1 = cntatoms.select_atoms("resid "+str(resid)).atoms[:args.ringsize]
            ring2 = cntatoms.select_atoms("resid "+str(resid)).atoms[-args.ringsize:]
            vector = numpy.abs( ring2.center_of_mass() - ring1.center_of_mass() )
            numpy.where(vector > 0.5 * dimensions, dimensions - vector, vector)

            # Calculate the cosine of the tilting angle
            tiltcos = vector[2]/numpy.linalg.norm(vector)
            allTiltcos[i] = tiltcos

        # Calculate average and standard deviation
        averTiltcos = numpy.average(allTiltcos)        
        stdvTiltcos = numpy.std(allTiltcos)

        # Write the results to the output file
        output.write("%12.3f \t %8f \t %8f \n" % (time, averTiltcos, stdvTiltcos))

