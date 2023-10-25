#!/bin/sh
#
# Fluidity buildenv script written for Viking with local PETSc, Zoltan, NetCDF, and conda
#
# NOTE: There are no niceties here. If something fails, the script carries on!

# Prefix for build
PREFIX=$HOME/scratch/fluidity-build/

# We want everything python to use conda
unset PYTHONPATH


# get and install vtk
curl -s https://vtk.org/files/release/8.2/VTK-8.2.0.tar.gz | tar -zxf -


# Get Conda
CONDA_SERVER=https://repo.continuum.io/
CONDA_SERVERDIR=miniconda/
# Python 3.6
CONDA_INSTALLER=Miniconda3-4.3.31-Linux-x86_64.sh

#curl -O -L ${CONDA_SERVER}${CONDA_SERVERDIR}${CONDA_INSTALLER}

# Install conda without prompting
bash ./${CONDA_INSTALLER} -b -f -p ${PREFIX}

# Point to the conda binary
export PATH="${PREFIX}/bin:$PATH"


# Install other packages needed by Fluidity
conda install -y numpy scipy 
conda install -y sympy matplotlib
conda install -c scitools -y udunits2 


# hdf5
HDF5_VER="1.10.7"
curl -s "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${HDF5_VER}/src/CMake-hdf5-${HDF5_VER}.tar.gz" | tar -zxf -
pushd CMake-hdf5-${HDF5_VER}/hdf5-${HDF5_VER}
./configure --prefix=${PREFIX} --enable-parallel --enable-fortran
make -j 5 install

popd

# Get the PETSc tarball
PETSC_VERSION=3.13.6
PETSC_SERVER=" https://codeload.github.com/petsc/petsc/tar.gz/v${PETSC_VERSION}"

curl -s ${PETSC_SERVER} | tar -zxf - 

PETSC_SOURCEDIR=petsc-${PETSC_VERSION}

# Change into the PETSc sourcedir
pushd ${PETSC_SOURCEDIR}

export PETSC_DIR=${PWD}

# Set up the PETSc build environment
PETSC_ARCH=linux-gnu-opt

# Configure and build PETSc
./configure \
 --with-debugging=0 \
 --PETSC_ARCH=${PETSC_ARCH} \
 --download-hypre=1 \
 --download-metis=1 \
 --download-ml=1 \
 --download-parmetis=1 \
 --with-sowing=1 \
 --download-triangle=1 \
 --download-ptscotch=1 \
 --download-suitesparse=1 \
 --download-ctetgen=1 \
 --download-chaco=1 \
 --download-scalapack=1 \
 --download-blacs=1 \
 --download-fblaslapack=1 \
 --known-mpi-shared-libraries=1 \
 --prefix=${PREFIX} \
 --with-fortran-interfaces=1 

make PETSC_DIR=$PWD PETSC_ARCH=${PETSC_ARCH} all
make PETSC_DIR=$PWD PETSC_ARCH=${PETSC_ARCH} install

# Revert to our previous location
popd

# link blas and lapack libraries to expected names
rm ${PREFIX}/lib/libblas.a ${PREFIX}/lib/liblapack.a
ln -s ${PREFIX}/lib/libfblas.a ${PREFIX}/lib/libblas.a 
ln -s ${PREFIX}/lib/libflapack.a ${PREFIX}/lib/liblapack.a 

# NetCDF-C build
# Obtain source
NETCDF_VERSION="4.5.0"
NETCDF_TARBALL="netcdf-${NETCDF_VERSION}.tar.gz"
NETCDF_SERVER="ftp://ftp.unidata.ucar.edu/"
NETCDF_SERVERDIR="pub/netcdf/"

curl -s ${NETCDF_SERVER}${NETCDF_SERVERDIR}${NETCDF_TARBALL} | tar -zxf -

NETCDF_SOURCEDIR=netcdf-${NETCDF_VERSION}

# Change into the source directory
pushd ${NETCDF_SOURCEDIR}

# Configure and build
LDFLAGS="-L${PREFIX}/lib" CPPFLAGS="-I${PREFIX}/include" \
   LIBS="-L${PREFIX}/lib" \
CC=mpicc ./configure --prefix=${PREFIX} --includedir=${PREFIX}/include

make install

# Return to our previous directory
popd

# NetCDF-Fortran build

# Obtain source
NETCDFF_VERSION="4.5.0"
NETCDFF_TARBALL="netcdf-fortran-${NETCDFF_VERSION}.tar.gz"

curl -s -L ${NETCDF_SERVER}${NETCDF_SERVERDIR}${NETCDFF_TARBALL} | tar -zxf -

NETCDFF_SOURCEDIR=netcdf-fortran-${NETCDFF_VERSION}

# Change into the source directory
pushd ${NETCDFF_SOURCEDIR}

# Configure and build
LDFLAGS="-L${PREFIX}/lib" CPPFLAGS="-I${PREFIX}/include" \
   LIBS="-L${PREFIX}/lib" \
CC=mpicc FC=mpif90   ./configure --prefix=${PREFIX} --includedir=${PREFIX}/include

make install

# Return to our previous directory
popd

# Obtain source for Zoltan
ZOLTAN_VERSION=v3.83
ZOLTAN_TAR=https://raw.githubusercontent.com/sandialabs/Zoltan/main/zoltan_distrib_${ZOLTAN_VERSION}.tar.gz

curl -s ${ZOLTAN_TAR} | tar --warning=no-unknown-keyword -zxf -

# Create the Zoltan build directory and change into it
ZOLTAN_SOURCEDIR=${PWD}/Zoltan_${ZOLTAN_VERSION}
ZOLTAN_BUILDDIR=${PWD}/zoltan-build

mkdir ${ZOLTAN_BUILDDIR}

pushd ${ZOLTAN_BUILDDIR}

# Set up the Zoltan build environment
ZOLTAN_LIBDIR=${PREFIX}/lib
ZOLTAN_INCDIR=${PREFIX}/include

# With some Zoltan-build-local environment, configure Zoltan
CC=mpicc FC=mpif90 F90=mpif90 CXX=mpic++ \
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ZOLTAN_LIBDIR} \
LDFLAGS="-L${PREFIX}/lib" CPPFLAGS="-I${PREFIX}/include" \
  ${ZOLTAN_SOURCEDIR}/configure \
  --prefix=${PREFIX} \
  --libdir=${ZOLTAN_LIBDIR} \
  --enable-mpi \
  --with-mpi-compilers \
  --with-gnumake \
  --enable-f90interface \
  --enable-zoltan-cppdriver \
  --disable-examples \
  --with-parmetis \
  --with-parmetis-libdir=${ZOLTAN_LIBDIR} \
  --with-parmetis-incdir=${ZOLTAN_INCDIR} \
  --with-scotch \
  --with-scotch-libdir=${ZOLTAN_LIBDIR} \
  --with-scotch-incdir=${ZOLTAN_INCDIR}

# Build and install Zoltan
make install

# Return to our previous location
popd

# We also need VTK.
# You have to configure by hand, but this will download and extract
curl -s https://vtk.org/files/release/8.2/VTK-8.2.0.tar.gz | tar -zxf -
cd VTK-8.2.0/
mkdir build
cd build

# from here you need to then run:
# ccmake ..
# Edit some options. Turn python to python 3.6, turn on parallel IO, set your PREFIX to:
# /users/jh1889/scratch/fluidity-build
# (Swampping my username for yours)
# Then press c to configuree. Assuming it succeeds, press g to generate and that should exit the config
# Then run:
# make -j 8 install
#
#
#
#
# That should be the build environment done
# Check out fluidity:
# git clone https://github.com/FluidityProject/fluidity.git
#
# This bit is horrid. For some reason libtool doesn't make a file that some other bit needs
# So copy it:
# cd fluidity
# cp /usr/share/libtool/config/ltmain.sh .
# Now configure:
# VTK_INSTALL_PREFIX=/users/jh1889/scratch/fluidity-build/include/ ./configure --prefix /users/jh1889/scratch/fluidity-build/ --with-lapack=/users/jh1889/scratch/fluidity-build/lib/liblapack.a --enable-2d-adaptivity
# make -j 5
# make install
#
# Fluidity is then in /users/jh1889/scratch/fluidity-build/bin/, but as this is in your path from the module, simply typing fluidity should run it.
