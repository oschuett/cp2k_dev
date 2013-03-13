#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import floor, ceil


#===============================================================================
def read_references():
    known_minima = {}
    f = open("LJ_known_minima.txt")
    for line in f.readlines():
        parts = line.split()
        known_minima[int(parts[0])] = parts[1]
    return known_minima

#===============================================================================
def main():
    known_minima = read_references()

    for s in range(2,40):
        fn = "LJ%.3d.inp"%s
        print "Writting: "+fn
        f = open(fn, "w")
        Emin = 0.00099999*float(known_minima[s])
        f.write(gen_input(size=s, Emin=Emin))
        f.close()

#===============================================================================
def gen_input(size, Emin):
    output = ""
    output += "&GLOBAL\n"
    output += "   PROGRAM_NAME GLOBAL_OPT\n"
    output += "   RUN_TYPE NONE\n"
    output += "   PROJECT_NAME LJ%.3d\n"%size
#    output += "   SEED 200\n"
    output += "&END GLOBAL\n"

    output += "&GLOBAL_OPT\n"
    output += "   NUMBER_OF_WALKERS  1\n"
    output += "   Emin %f\n"%Emin
    output += "&END GLOBAL_OPT\n"

    output += """
&MOTION
  &PRINT
    &RESTART ! writting restarts is expensive, turning them off
       &EACH
         MD -1
         GEO_OPT -1
       &END EACH
       ADD_LAST NO
    &END RESTART
  &END PRINT

  &MD
    ENSEMBLE NVE
    STEPS 1000
    TIMESTEP 1.0
    TEMPERATURE 300
    STEP_START_VAL 1 !otherwise md_energies::md_write_output flushes trajectory
  &END MD

  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 3000
    !MAX_DR 0.0001
    &BFGS
     USE_RAT_FUN_OPT  ! otherwise LJ particle sth. get too close.

      &RESTART  ! writting restarts is expensive, turning them off
        &EACH
          GEO_OPT -1
        &END EACH
        ADD_LAST NO
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
