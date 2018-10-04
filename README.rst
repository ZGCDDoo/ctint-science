==========================================================================
 CT-INT-TG : Continuous time interaction expansion Tremblay Group
==========================================================================

:Authors: Charles-David Hébert, Maxime Charlebois, Patrick Sémon 

Installation
================================


**Note :**
If build problems,
please remove the build directory if it exists, then retry :
    
    $ rm -rf build

Dependencies
--------------
1. Armadillo
2. boost (mpi, serialization, filesystem, system)


Pre-Steps
----------
1. Make sure you have a "bin" directory in your home folder
2. Append the bin folder to your path. Add the following line to your ~/.bashrc:  export PATH="$PATH:~/bin"
3. $ source ~/.bashrc

Linux (Ubuntu 16.04 and 18.04)
--------------------------------
This installation procedure should work for many recent Linux flavors. For the following
we present the instructions specific for Ubuntu or derivatives.

1. Install the Dependencies
    $ sudo apt-get install libarmadillo-dev libboost-all-dev cmake liblapack-dev
2. | $ mkdir build && cd build && cmake .. && make install



Example
================================
1. Go to the examples/U3_b10 directory:
   $ cd examples/U3_b10
   $ bash run_dmft.sh

2. The important output is the selfenergy, given at each dmft iteration, by self+{$ITER}.dat .

    