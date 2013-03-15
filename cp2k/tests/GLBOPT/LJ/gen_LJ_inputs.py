#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import floor, ceil
from os import mkdir
from os import path

#===============================================================================
def read_references():
    known_minima = {}
    f = open("known_LJ_minima.txt")
    for line in f.readlines():
        parts = line.split()
        known_minima[int(parts[0])] = parts[1]
    return known_minima

#===============================================================================
def main():
    size_min, size_max = 2, 5
    run_min, run_max   = 1, 3


    dirname = "Farming_LJ%.3d-%.3d_RUN%.3d-%.3d"%(size_min, size_max, run_min, run_max)
    print "Creating dir: "+dirname
    mkdir(dirname)

    jobs = []
    known_minima = read_references()
    for r in range(run_min, run_max):
        for s in range(size_min, size_max):
            fn = "LJ%.3d_RUN%.3d.inp"%(s,r)
            jobs.append(fn)
            Emin = 0.001*float(known_minima[s]) + 1.0e-6
            write2file(dirname+"/"+fn, gen_glbopt_input(size=s, Emin=Emin, run=r))

    write2file(dirname+"/"+dirname+".inp", gen_framing_input(jobs))

#===============================================================================
def write2file(fn, content):
    print "Writting: "+fn
    f = open(fn, "w")
    f.write(content)
    f.close()

#===============================================================================
def gen_framing_input(jobs):
    output = ""
    output += "&GLOBAL\n"
    output += "   PROGRAM_NAME FARMING\n"
    output += "   RUN_TYPE NONE\n"
    output += "   PROJECT_NAME LJ_framing\n"
    output += "&END GLOBAL\n"

    output += "&FARMING\n"
    output += "   GROUP_SIZE 1\n"
    output += "   MASTER_SLAVE .TRUE.\n"

    for j in jobs:
        output += "   &JOB\n"
        output += "       INPUT_FILE_NAME  %s\n"%path.basename(j)
        output += "       OUTPUT_FILE_NAME %s\n"%path.basename(j).replace("inp","out")
        output += "       DIRECTORY .\n"
        output += "   &END JOB\n"

    output += "&END FARMING\n"

    return(output)

#===============================================================================
def gen_glbopt_input(size, Emin, run):
    output = ""
    output += "&GLOBAL\n"
    output += "   PROGRAM_NAME GLOBAL_OPT\n"
    output += "   RUN_TYPE NONE\n"
    output += "   PROJECT_NAME LJ%.3d_RUN%.3d\n"%(size, run)
    output += "   SEED %d\n"%(100*run)
    output += "&END GLOBAL\n"

    output += "&GLOBAL_OPT\n"
    output += "   NUMBER_OF_WALKERS  1\n"
    output += "   Emin %.10f\n"%Emin
    output += "&END GLOBAL_OPT\n"

    output += """
&MOTION
  &PRINT  ! IO is expensive, turning everything off
    &RESTART OFF
    &END RESTART
    &RESTART_HISTORY OFF
    &END RESTART_HISTORY
    &TRAJECTORY OFF
    &END TRAJECTORY
  &END PRINT

  &MD
    ENSEMBLE NVE
    STEPS 100
    TIMESTEP 1.0
    TEMPERATURE 10.0
    &PRINT   ! IO is expensive, turning everything off
      &ENERGY OFF
      &END ENERGY
    &END PRINT
  &END MD

  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 300
    &BFGS
     USE_RAT_FUN_OPT  ! otherwise LJ particle sth. get too close.
     &RESTART OFF     ! IO is expensive, turning everything off
     &END RESTART
    &END BFGS
  &END GEO_OPT

&END MOTION

&FORCE_EVAL
  METHOD FIST
  &MM
    &FORCEFIELD
     &SPLINE
        R0_NB 0.0     ! solely MAX_SPLINE shall control spline range
        EMAX_SPLINE   [hartree] 1000    ! yields r_min = 0.66 bohr
        EMAX_ACCURACY [hartree] 1000
        EPS_SPLINE    [hartree] 1.0E-10   ! yields 1698 spline points
     &END SPLINE
      &NONBONDED
        &LENNARD-JONES
          atoms X X
          EPSILON [hartree] 0.001
          SIGMA 1.0
          RCUT 25.0
        &END LENNARD-JONES
      &END NONBONDED
      &CHARGE
        ATOM X
        CHARGE 0.0
      &END CHARGE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE none
      &END EWALD
    &END POISSON
    &PRINT
      &FF_INFO
      &END FF_INFO
    &END PRINT
  &END MM
  &SUBSYS
    &CELL
    ABC [angstrom] 50.0 50.0 50.0
      !PERIODIC NONE
   &END CELL

   &COORD
      """

    output += gen_coords(size=size, lattice_const=1.5)

    output += """
    &END COORD

    &TOPOLOGY
      CONNECTIVITY OFF
    &END TOPOLOGY

    &COLVAR
      &U
      &END U
    &END COLVAR

    !&COLVAR
    !  &DISTANCE
    !    ATOMS 1 2
    !  &END DISTANCE
    !&END COLVAR


     &KIND X
        ELEMENT H
        MASS 1.0
     &END KIND
  &END SUBSYS
  STRESS_TENSOR ANALYTICAL
&END FORCE_EVAL
"""
    return output

#===============================================================================
def gen_coords(size, lattice_const):
    #output = "%d\n\n"%size
    output = ""
    nx = floor(size**(1.0/3))
    ny = floor((size/nx)**(1.0/2))
    #print "nx: ",nx, "ny: ",ny

    n = 0
    for i in range(int(ceil(size/nx/ny))+1):
        z = i*lattice_const
        for j in range(int(nx)):
            y = j*lattice_const
            for k in range(int(ny)):
                x = k*lattice_const
                n += 1
                if(n>size): return output
                #print n, i, j, k
                output += "X %f %f %f\n"%(x, y, z)

    assert(False) #should always exit through the return-statement

#===============================================================================
if(__name__ == "__main__"):
    main()
#EOF
