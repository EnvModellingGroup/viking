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
load(" tools/Mercurial/3.8.3-foss-2018b-Python-2.7.15")
unload("toolchain/foss/2018b")
unload("tools/numactl/2.0.11-GCCcore-7.3.0")
unload("tools/XZ/5.2.4-GCCcore-7.3.0")
unload("lib/libxml2/2.9.8-GCCcore-7.3.0")
unload("system/libpciaccess/0.14-GCCcore-7.3.0")
unload("system/hwloc/1.11.10-GCCcore-7.3.0")
unload("mpi/OpenMPI/3.1.1-GCC-7.3.0-2.30")
unload("numlib/OpenBLAS/0.3.1-GCC-7.3.0-2.30")
unload("toolchain/gompi/2018b")
unload("numlib/FFTW/3.3.8-gompi-2018b")
unload("numlib/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1")
unload("tools/bzip2/1.0.6-foss-2018b")
unload("devel/ncurses/6.1-GCCcore-7.3.0")
unload("lib/libreadline/7.0-GCCcore-7.3.0")
unload("lang/Tcl/8.6.8-GCCcore-7.3.0")
unload("devel/SQLite/3.24.0-GCCcore-7.3.0")
unload("math/GMP/6.1.2-GCCcore-7.3.0 lib/libffi/3.2.1-foss-2018b")
load("compiler/GCC/10.3.0")
load("devel/Autoconf/2.71-GCCcore-10.3.0")
load("devel/Automake/1.16.3-GCCcore-10.3.0")
load("lang/Python/3.9.5-GCCcore-10.3.0")
load("devel/CMake/3.20.1-GCCcore-10.3.0")
load("tools/git/2.32.0-GCCcore-10.3.0-nodocs")
load("lib/libtool/2.4.6-GCCcore-10.3.0")
load("lang/Bison/3.7.6-GCCcore-10.3.0")
load("lib/UCX/1.10.0-GCCcore-10.3.0")
load("lib/zlib/1.2.11-GCCcore-10.3.0")
load("system/hwloc/2.4.1-GCCcore-10.3.0")
load("mpi/MPICH/3.4.2-GCC-10.3.0") 
prepend_path("PATH","/mnt/lustre/groups/env-tsunami-2019/gdal/bin")
prepend_path("LD_LIBRARY_PATH","/mnt/lustre/groups/env-tsunami-2019/gdal/lib")
prepend_path("LD_LIBRARY_PATH","/mnt/lustre/groups/env-tsunami-2019/lapack")
setenv("GDAL_DATA","-L/mnt/lustre/groups/env-tsunami-2019/gdal/share/gdal")
prepend_path("LIBRARY_PATH","/mnt/lustre/groups/env-tsunami-2019/lapack")
setenv("LAPACK","/mnt/lustre/groups/env-tsunami-2019/lapack")
