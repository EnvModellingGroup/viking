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
load("GCC/10.3.0")
load("Autoconf/2.71-GCCcore-10.3.0")
load("Automake/1.16.3-GCCcore-10.3.0")
load("Python/3.9.5-GCCcore-10.3.0")
load("CMake/3.20.1-GCCcore-10.3.0")
load("git/2.32.0-GCCcore-10.3.0-nodocs")
load("libtool/2.4.6-GCCcore-10.3.0")
load("Bison/3.7.6-GCCcore-10.3.0")
load("UCX/1.10.0-GCCcore-10.3.0")
load("zlib/1.2.11-GCCcore-10.3.0")
load("MPICH/4.1.2-GCC-10.3.0")
load("PROJ/8.0.1-GCCcore-10.3.0")
load("JasPer/2.0.28-GCCcore-10.3.0")
prepend_path("PATH","/mnt/scratch/projects/env-tsunami-2019/gdal/bin")
prepend_path("LD_LIBRARY_PATH","/mnt/scratch/projects/env-tsunami-2019/gdal/lib")
setenv("GDAL_DATA","-L/mnt/scratch/projects/env-tsunami-2019/gdal/share/gdal")
