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
    run_min, run_max   = 1, 1024

    method = "MinhopIndep"
    dirname = "LJ38_%s_RUN%.3d-%.3d"%(method,run_min, run_max)
    print "Creating dir: "+dirname
    mkdir(dirname)

    jobs = []
    known_minima = read_references()
    for r in range(run_min, run_max+1):
       fn = "LJ38_RUN%.3d.inp"%(r,)
       jobs.append(fn)
       Emin = 0.001*float(known_minima[38]) + 1.0e-6
       inp = gen_glbopt_input(size=38, Emin=Emin, run=r, method=method)
       write2file(dirname+"/"+fn, inp)

    write2file(dirname+"/"+dirname+".inp", gen_framing_input(jobs))
    write2file(dirname+"/"+dirname+".job", gen_jobfile(dirname))


#===============================================================================
def main2():
    size_min, size_max = 38,38
    run_min, run_max   = 1, 50


    #dirname = "Farming_LJ%.3d-%.3d_RUN%.3d-%.3d"%(size_min, size_max, run_min, run_max)
    dirname = "Farming_LJ38_RUN%.3d-%.3d"%(run_min, run_max)
    print "Creating dir: "+dirname
    mkdir(dirname)

    jobs = []
    known_minima = read_references()
    for r in range(run_min, run_max+1):
        for s in range(size_min, size_max+1):
            fn = "LJ%.3d_RUN%.3d.inp"%(s,r)
            jobs.append(fn)
            Emin = 0.001*float(known_minima[s]) + 1.0e-6
            write2file(dirname+"/"+fn, gen_glbopt_input(size=s, Emin=Emin, run=r))

    write2file(dirname+"/"+dirname+".inp", gen_framing_input(jobs))
    write2file(dirname+"/"+dirname+".job", gen_jobfile(dirname))

#===============================================================================
def write2file(fn, content):
    print "Writting: "+fn
    f = open(fn, "w")
    f.write(content)
    f.close()


#===============================================================================
def gen_jobfile(name):
   output =  "#! /bin/bash\n"
   output += "#SBATCH --job-name=cp2k\n"
   output += "#SBATCH --ntasks=1024\n"
   output += "#SBATCH --partition=day\n"
#   output += "#SBATCH --partition=night\n"
   output += "#SBATCH --time=0:59:00\n"
   output += "#SBATCH --account=h05\n\n"

   output += "cd $SLURM_SUBMIT_DIR\n"
   output += "export OMP_NUM_THREADS=1\n\n"

   output += "aprun -n 1024  -d 1  -N 16 -cc cpu  ../../../../exe/todi/cp2k.popt -i %s.inp -o %s.out\n"%(name,name)
   return(output)  



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
    output += "   MASTER_SLAVE .FALSE.\n"

    for j in jobs:
        output += "   &JOB\n"
        output += "       INPUT_FILE_NAME  %s\n"%path.basename(j)
        output += "       OUTPUT_FILE_NAME %s\n"%path.basename(j).replace("inp","out")
        output += "       DIRECTORY .\n"
        output += "   &END JOB\n"

    output += "&END FARMING\n"

    return(output)

#===============================================================================
def gen_glbopt_input(size, Emin, run, method):
    assert(method == "MinhopIndep")

    output = ""
    output += "&GLOBAL\n"
    output += "   PROGRAM_NAME GLOBAL_OPT\n"
    output += "   RUN_TYPE NONE\n"
    output += "   PROJECT_NAME LJ%.3d_RUN%.3d\n"%(size, run)
    output += "   SEED %d\n"%(100*run)
    output += "   &TIMINGS\n"
    output += "      THRESHOLD 0.0\n"
    output += "   &END TIMINGS\n"
    output += "&END GLOBAL\n"

    output += "&GLOBAL_OPT\n"
    output += "   NUMBER_OF_WALKERS  1\n"
    #output += "   MAX_ITER 1000\n"
    output += "   E_MIN %.10f\n"%Emin
    output += "   &MINIMA_HOPPING\n"
    output += "      SHARE_HISTORY FALSE\n"
    output += "   &END MINIMA_HOPPING\n"
    output += "&END GLOBAL_OPT\n"

    output += """
&MOTION
  &PRINT  ! IO is expensive, turning everything off
    &RESTART OFF
    &END RESTART
    &RESTART_HISTORY OFF
    &END RESTART_HISTORY
    &TRAJECTORY
      ADD_LAST NUMERIC
      &EACH
        GEO_OPT -1
        MD -1
      &END EACH
    &END TRAJECTORY
  &END PRINT

  &MD
    ENSEMBLE NVE
    STEPS 100
    TIMESTEP 1.0
    TEMPERATURE 10.0

    &PRINT
      &ENERGY OFF
      &END ENERGY
      &CENTER_OF_MASS OFF
      &END CENTER_OF_MASS
      &COEFFICIENTS OFF
      &END COEFFICIENTS
      &PROGRAM_RUN_INFO OFF
      &END PROGRAM_RUN_INFO
      &ROTATIONAL_INFO OFF
      &END ROTATIONAL_INFO
      &SHELL_ENERGY OFF
      &END SHELL_ENERGY
      &TEMP_KIND OFF
      &END TEMP_KIND
      &TEMP_SHELL_KIND OFF
      &END TEMP_SHELL_KIND
      FORCE_LAST .TRUE.
    &END PRINT
  &END MD

  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 300
    &BFGS
     TRUST_RADIUS [angstrom] 0.1
     USE_RAT_FUN_OPT  ! otherwise LJ particle sth. get too close.
     &RESTART OFF
     &END RESTART
    &END BFGS
    &PRINT
      &PROGRAM_RUN_INFO OFF
      &END PROGRAM_RUN_INFO
    &END PRINT
  &END GEO_OPT

&END MOTION

&FORCE_EVAL
  &PRINT
    &DISTRIBUTION OFF
    &END DISTRIBUTION
    &DISTRIBUTION1D OFF
    &END DISTRIBUTION1D
    &DISTRIBUTION2D OFF
    &END DISTRIBUTION2D
    &FORCES OFF
    &END FORCES
    &GRID_INFORMATION OFF
    &END GRID_INFORMATION
    &PROGRAM_RUN_INFO OFF
    &END PROGRAM_RUN_INFO
    &STRESS_TENSOR OFF
    &END STRESS_TENSOR
    &TOTAL_NUMBERS OFF
    &END TOTAL_NUMBERS
  &END PRINT


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
      &DERIVATIVES OFF
      &END DERIVATIVES
      &DIPOLE OFF
      &END DIPOLE
      &EWALD_INFO OFF
      &END EWALD_INFO
      &FF_INFO OFF
      &END FF_INFO
      &FF_PARAMETER_FILE OFF
      &END FF_PARAMETER_FILE
      &ITER_INFO OFF
      &END ITER_INFO
      &NEIGHBOR_LISTS OFF
      &END NEIGHBOR_LISTS
      &PROGRAM_BANNER OFF
      &END PROGRAM_BANNER
      &PROGRAM_RUN_INFO OFF
      &END PROGRAM_RUN_INFO
      &SUBCELL OFF
      &END SUBCELL
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
