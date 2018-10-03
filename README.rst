==========================================================================
 CT-INT-TG : Continuous time interaction expansion Tremblay Group
==========================================================================

:Authors: Charles-David Hébert, Maxime Charlebois, Patrick Sémon 
:Date: $Date: 2017-09-21 01:10:53 +0000 (Fri, 21 Sept 2017) $
:Revision: $Revision: 1.01 $
:Description: Description

Graham
-------

g++
^^^^^^

* module reset 
* module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi
* mkdir build && cd build && cmake -Dhome=OFF -Dgraham=ON .. && make

icpc (mpic++)
^^^^^^^^^^^^^^
* module reset
* module load intel/2017.5 armadillo boost-mpi
* mkdir build && cd buil && cmake -Dhome=OFF -Dgraham=ON .. && make

Mp2
------

g++ and icpc (mpic++)
^^^^^^^^^^^^^^^^^^^^^^
* module reset
* module load cmake/3.6.1  gcc/6.1.0  intel64/17.4  boost64/1.65.1_intel17 openmpi/1.8.4_intel17  armadillo/8.300.0
* mkdir build && cd buil && cmake -Dhome=OFF -Dmp2=ON .. && make



    