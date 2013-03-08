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
        Estop = 0.99999*float(known_minima[s])
        f.write(gen_input(size=s, Estop=Estop))
        f.close()

#===============================================================================
def gen_input(size, Estop):
    output = ""
    output += "&GLOBAL\n"
    output += "   PROGRAM_NAME GLOBAL_OPT\n"
    output += "   RUN_TYPE NONE\n"
    output += "   PROJECT_NAME LJ%.3d\n"%size
    output += "&END GLOBAL\n"

    output += "&GLOBAL_OPT\n"
    output += "   NUMBER_OF_WALKERS  10\n"
    output += "   ENERGY_STOP %f\n"%Estop
    output += "&END GLOBAL_OPT\n"

    output += """
&MOTION

  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 3000
    !MAX_DR 0.0001
  &END GEO_OPT

  &MD
    ENSEMBLE NVE
    STEPS 3000
    TIMESTEP 0.5
    TEMPERATURE 1.0
  &END MD
&END MOTION

&FORCE_EVAL
  METHOD FIST
  &MM
    &FORCEFIELD
     &SPLINE
        EMAX_ACCURACY 1.0E12
        EMAX_SPLINE 1.0E12
        EPS_SPLINE 1.0E-12
      &END
      &NONBONDED
        &LENNARD-JONES
          atoms Ar Ar
          EPSILON [hartree] 1.0
          SIGMA 1.0
          RCUT 50.0
        &END LENNARD-JONES
      &END NONBONDED
      &CHARGE
        ATOM Ar
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
      &END
    &END
  &END MM
  &SUBSYS
    &CELL
      ABC 100.0 100.0 100.0
      PERIODIC NONE
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
                output += "Ar %f %f %f\n"%(x, y, z)

    assert(False) #should always exit through the return-statement

#===============================================================================
if(__name__ == "__main__"):
    main()
#EOF
