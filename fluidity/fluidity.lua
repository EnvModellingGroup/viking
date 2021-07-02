help([==[

Description
===========
Fluidity is a CFD solver for a range of applications


More information
================
 - Homepage: http://fluidityproject.github.io/
]==])
whatis([==[Description: Fluidity - CFD sovler]==])
whatis([==[Homepage: http://fluidityproject.github.io/]==])
load("compiler/GCCcore/7.3.0")
load("mpi/OpenMPI/3.1.1-GCC-7.3.0-2.30")
load("lib/libreadline/8.0-GCCcore-7.3.0")
load("numlib/OpenBLAS/0.3.7-GCC-7.3.0-2.30")
load("tools/binutils/2.30-GCCcore-7.3.0")
load("lang/Tcl/8.6.8-GCCcore-7.3.0")
load("tools/numactl/2.0.11-GCCcore-7.3.0")
load("numlib/FFTW/3.3.8-gompi-2018b")
load("numlib/SuiteSparse/5.1.2-foss-2018b-METIS-5.1.0")
load("devel/SQLite/3.24.0-GCCcore-7.3.0")
load("tools/XZ/5.2.4-GCCcore-7.3.0")
load("numlib/ScaLAPACK/2.0.2-gompi-2018b-OpenBLAS-0.3.1")
load("numlib/Hypre/2.14.0-foss-2018b")
load("tools/bzip2/1.0.6-GCCcore-7.3.0")
load("lib/libxml2/2.9.8-GCCcore-7.3.0")
load("lib/zlib/1.2.11-GCCcore-7.3.0")
load("devel/CMake/3.12.1-GCCcore-7.3.0")
load("system/libpciaccess/0.14-GCCcore-7.3.0")
load("devel/Boost/1.66.0-foss-2018b")
load("devel/ncurses/6.1-GCCcore-7.3.0")
load("system/hwloc/1.11.10-GCCcore-7.3.0")
setenv("PETSC_DIR","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/petsc-3.13.6/")
prepend_path("PYTHONPATH","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/lib/python3.6/site-packages")
prepend_path("PATH","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/bin/")
setenv("LDFLAGS","-L/mnt/lustre/groups/env-tsunami-2019/fluidity-build/lib/ -L/mnt/lustre/groups/env-tsunami-2019/fluidity-build/lib64/")
prepend_path("LD_LIBRARY_PATH","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/lib/")
prepend_path("LD_LIBRARY_PATH","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/lib64/")
setenv("PYTHONPATH","/mnt/lustre/groups/env-tsunami-2019/fluidity-build/python")



