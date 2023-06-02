#!/bin/bash

# Defines how scratch directory is made, might change with different systems
export scratch=/scratch

# Full path to directory with Chemshell binaries
export chemshellpath=/users/home/ragnarbj/chemshell-360/chemshell-ifort-new-polfix-PAR-DLFINDtaskfarmfix-openmpi1-10-icc-xtbaddition2
# Full path to Psi4 directory
export PYTHONHOME=/users/work/ragnarbj/psi4conda_1.1
# Full path to main ORCA directory, only used in serial Chemshell/parallel ORCA calculations
export orcapath=/users/work/ragnarbj/orca_4_0_1_linux_x86-64_openmpi202

# Printing paths to output
echo "Chemshell path is: $chemshellpath"
echo "Psi4 path is $PYTHONHOME"
echo "ORCA path is $orcapath"
