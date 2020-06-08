# Quandary - Quantum control on HPC clusters
This project implements a parallel-in-time optimization solver for quantum control. The underlying quantum dynamics model open quantum systems, using the Lindblad master equation to evolve the density matrix in time. The control problem aims to find control pulses that realize a certain gate, i.e. drive the system to a desired target state

## Requirements:
To build this project, you need to have the following packages installed:
* Petsc [https://www.mcs.anl.gov/petsc/]
* Xbraid [https://github.com/XBraid/xbraid], on branch 'solveadjointwithxbraid'

## Installation
* Download XBraid, switch to the 'solveadjointwithxbraid' branch and build the shared library:
    - git clone https://github.com/XBraid/xbraid.git
    - cd xbraid
    - git checkout solveadjointwithxbraid
    - make braid
* Install Petsc (see Petsc manual for installation guide). Set the `PETSC_DIR` and `PETSC_ARCH` variables.
* In the main directory of this project, adapt the beginning of the Makefile to set the path to XBraid and Petsc. 
* Type `make cleanup` to clean the build directory.
* Type `make -j main` to build the code. 


## Running
The code builds into the executable `./main`. It takes one argument being the name of the config file. The config file `AxC.cfg`, lists all possible config options. It is currently set to run the Alice-Cavity testcase (3x20 levels).
* `./main AxC.cfg` for serial run
* `srun -n36 ./main AxC.cfg` for parallel run using 36 cores

Make sure to include the Petsc directory in the `LD_LIBRARY_PATH`:
* `export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH`

### Notes to install and run on Quartz (LLNL, LC):
On quartz, you'll have to load modules first. In a terminal, you can type `module avail` to see available modules, `module list` to see what is currently loaded, and `module load <nameofmodule>` to load a module. For this project, you'll need to load a newer gcc compiler, and a module for MPI.
> `module load gcc/8.1.0
> `module load mvapich2/2.3

Petsc is already installed on Quartz, under `/usr/tce/packages/petsc/petsc-3.12.4-mvapich2-2.3-gcc-4.8-redhat`. Set the `PETSC_DIR` variable with
> `export PETSC_DIR = /usr/tce/packages/petsc/petsc-3.12.4-mvapich2-2.3-gcc-4.8-redhat`
and the `LD_LIBRARY_PATH` 
> `export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:$PETSC_DIR`

You'll have to do this in each new terminal on quartz. Alternatively, you can put all this into the `~/.bashrc` file, which should be loaded whenever you start a new terminal. You can also load the file yourself with `source ~/.bashrc`.


