# CTIN-SCIENCE


[![Build Status](https://travis-ci.com/ZGCDDoo/ctint-science.svg?branch=master)](https://travis-ci.com/ZGCDDoo/ctint-science)

[![Snap Status](https://build.snapcraft.io/badge/ZGCDDoo/ctint-science.svg)](https://build.snapcraft.io/user/ZGCDDoo/ctint-science)

[![Get it from the Snap Store](https://snapcraft.io/static/images/badges/en/snap-store-black.svg)](https://snapcraft.io/ctint-science)


## Notes

This program reproduces the DMFT results that can be found in the following article: https://arxiv.org/abs/1802.09456 .
However, this program is not optimized for speed. The main reasons being ease of use and installation. Please contact the developper to discuss obtaining the optimized version.


## Installation


### Preferred way: Snap package

The preferred and easiest way to install is to use the snap package. Follow the instructions to install snap capability here: [Snap Installation](https://docs.snapcraft.io/installing-snapd/6735). The snap package can be obtained by clicking the snap store button at the top.


### Manual installation

**Note :**
If build problems,
please remove the build directory if it exists, then retry :
    
    $ rm -rf build

#### Dependencies
1. Armadillo
2. boost (serialization, filesystem, system)


#### Pre-Steps
1. Make sure you have a "bin" directory in your home folder
2. Append the bin folder to your path. Add the following line to your ~/.bashrc:  export PATH="$PATH:~/bin"
3. $ source ~/.bashrc


#### Mac (Tested on macOS 10.13.6)

1. Install the Dependencies (with Homebrew : https://brew.sh/)
      * $ brew install armadillo
      * $ brew install boost

2. $ mkdir build && cd build && cmake .. && make install


#### Linux (Ubuntu 16.04 and 18.04)

This installation procedure should work for many recent Linux flavors. For the following
we present the instructions specific for Ubuntu or derivatives.

1. Install the Dependencies
    $ sudo apt-get install libarmadillo-dev libboost-all-dev cmake liblapack-dev
2. | $ mkdir build && cd build && cmake .. && make install



## Example
1. Go to the examples/U3_b10 directory:
   ```bash
    $ cd examples/U3_b10
    $ bash run_dmft.sh
   ```

2. The important output is the selfenergy, given at each dmft iteration, by self+\{\$ITER\}.dat .   
   The columns are the matsubara frequencies, the real and imaginary parts of the self-energy.  
   The parameter files is given by params{$ITER}.json, which is in big part self-explainatory.


### Parameters


We use json as the parameter file. Please keep the same structure and the same names (with case) as in the examples.

* SEED: The seed for the random number generator
* beta: Inverse temperature
* mu : Chemical potentiel
* U : Hubbard Interaction
* NMAT : The number of matsubara frequencies
* NTAU : The discretization in imaginary time. Should be at least 1000
* UPDATESMEAS : The number of updates bewteen each measure
* THERMALIZATION : The number of UPDATESMEAS for the thermalization, i.e , do UPDATESMEAS*THERMALIZATION updates before measuring 
* TOTALNMEAS : The number of measures
* CLEANUPDATE : The number of UPDATESMEAS before a cleanupdate
* delta : Value for the auxiliary Ising Spins
* HybFile: the name of the hybridization file for the current iteration

## Used third-party tools

* Json Libraray: https://github.com/nlohmann/json
* armadillo : http://arma.sourceforge.net/
* boost: https://www.boost.org/
    
   
## Contact and help

To get help, please leave an issue or contact me by email at charles-david.hebert@usherbrooke.ca

## Authors:
Charles-David Hébert.

## Acknowledgments

Maxime charlebois and Patrick Sémon have inspired some portions of the code.
