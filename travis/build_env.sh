#!/bin/bash
os=`uname`
arch=`uname -m`
 if [[ "$os" == "Darwin" ]]; then 
#  HOMEBREW_NO_AUTO_UPDATE=1 brew cask uninstall oclint || true  
#  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" openblas python3 ||true
  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" python3 ||true
#  if [[ "$MPI_IMPL" == "openmpi" ]]; then
#      HOMEBREW_NO_INSTALL_CLEANUP=1 HOMEBREW_NO_AUTO_UPDATE=1 brew install scalapack
#  fi
fi
if [[ "$os" == "Linux" ]]; then 
    if [[ "$MPI_IMPL" == "openmpi" ]]; then
	mpi_bin="openmpi-bin" ; mpi_libdev="libopenmpi-dev" scalapack_libdev="libscalapack-openmpi-dev"
    fi
    if [[ "$MPI_IMPL" == "mpich" ]]; then
        mpi_bin="mpich" ; mpi_libdev="libmpich-dev" scalapack_libdev="libscalapack-mpich-dev"
    fi
    cat /etc/apt/sources.list
    sudo add-apt-repository universe && sudo apt update
#    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin" "$scalapack_libdev"  make perl  libopenblas-dev python3 rsync
    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin"  make perl  python3 rsync
    if [[ "$MPI_IMPL" == "mpich" ]]; then
#fix for github actions	
	sudo ln -sf /usr/bin/mpifort.mpich /etc/alternatives/mpif90
	sudo ln -sf /etc/alternatives/mpif90  /usr/bin/mpif90
    fi
fi
