#!/usr/bin/env csh
#
# Use this script to assure that modules are properly loaded before compiling
#
setenv MPILIB mpi-serial
setenv COMPILER intel
source env_mach_specific
gmake $*
