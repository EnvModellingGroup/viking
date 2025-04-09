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
load("Python/3.11.3-GCCcore-12.3.0")
load("OpenMPI/4.1.5-GCC-12.3.0")
load("Bison/3.8.2-GCCcore-12.3.0")
load("M4/1.4.19-GCCcore-12.3.0")
load("CMake/3.26.3-GCCcore-12.3.0")
load("OpenBLAS/0.3.23-GCC-12.3.0")
load("FFTW/3.3.10-GCC-12.3.0")
load("hwloc/2.9.1-GCCcore-12.3.0")
load("libjpeg-turbo/2.1.5.1-GCCcore-12.3.0")
load("LibTIFF/4.5.0-GCCcore-12.3.0")
load("libdeflate/1.18-GCCcore-12.3.0")
prepend_path("PATH","/mnt/scratch/projects/env-tsunami-2019/gdal-3.10/bin")
prepend_path("LD_LIBRARY_PATH","/mnt/scratch/projects/env-tsunami-2019/gdal-3.10/lib64")
setenv("GDAL_DATA","-L/mnt/scratch/projects/env-tsunami-2019/gdal/share/gdal-3.10")
