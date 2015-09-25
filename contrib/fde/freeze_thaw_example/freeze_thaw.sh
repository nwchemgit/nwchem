#!/bin/bash

if [ ! -f h2o1.movecs ] || [ ! -f h2o2.movecs ]; then 
  mpirun -np 4 ~/Programs/nwchem-dev-new-build/bin/LINUX64/./nwchem dimer.nw > dimer.out.rhf;
  mpirun -np 4 ~/Programs/nwchem-dev-new-build/bin/LINUX64/./nwchem first.nw > first.out.rhf;
  rm h2o*.out.rhf
fi

export how_many_iterations=0
while [ "$how_many_iterations" != "4" ]; 
do 
#  ~/Programs/nwchem-dev/bin/LINUX64/./nwchem h2o2.nw > h2o2.out.rhf; 
#  ~/Programs/nwchem-dev/bin/LINUX64/./nwchem h2o1.nw > h2o1.out.rhf; 
  mpirun -np 4 ~/Programs/nwchem-dev-new-build/bin/LINUX64/./nwchem h2o2.nw > h2o2.out.rhf;
  mpirun -np 4 ~/Programs/nwchem-dev-new-build/bin/LINUX64/./nwchem h2o1.nw > h2o1.out.rhf;
  export how_many_iterations=`grep " d=" h2o1.out.rhf h2o2.out.rhf | wc -l`
  echo "Last self-consistent loop has $how_many_iterations iterations in total. Trying to reach 4!";
done

