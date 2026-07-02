help([==[

Description
===========
Firedrake is a FEM solver for a range of applications


More information
================
 - Homepage: https://firedrakeproject.org/
]==])
whatis([==[Description: Firedrake - FEM sovler]==])
whatis([==[Homepage: https://firedrakeproject.org/]==])
load("GDAL/3.11.3-foss-2025b")
load("Python/3.13.5-GCCcore-14.3.0")
load("OpenMPI/5.0.8-GCC-14.3.0")
load("Bison/3.8.2-GCCcore-14.3.0")
load("M4/1.4.20-GCCcore-14.3.0")
load("CMake/3.31.8-GCCcore-14.3.0")
load("OpenBLAS/0.3.30-GCC-14.3.0")
load("FFTW/3.3.10-GCC-14.3.0")
setenv("OMP_NUM_THREADS","1")
setenv("OPENBLAS_NUM_THREADS","1")
setenv("OMPI_MCA_orte_precondition_transports","13241234acffedeb-abcdefabcdef1233")
setenv("PETSC_DIR","/mnt/scratch/projects/env-tsunami-2019/petsc-2026.4.1")
setenv("PETSC_ARCH","arch-firedrake-default")
