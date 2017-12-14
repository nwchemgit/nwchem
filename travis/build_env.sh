#!/bin/bash
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then 
  brew cask uninstall oclint || true  
  brew install gcc "$MPI_IMPL" openblas ||true
  if [[ "$MPI_IMPL" == "openmpi" ]]; then
     brew install scalapack
  fi
fi
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then 
    if [[ "$MPI_IMPL" == "openmpi" ]]; then
	mpi_bin="openmpi-bin" ; mpi_libdev="libopenmpi-dev"
    fi
    if [[ "$MPI_IMPL" == "mpich" ]]; then
        mpi_bin="mpich" ; mpi_libdev="libmpich-dev"
    fi
    sudo apt-get -y install gfortran gcc python-dev  cmake "$mpi_libdev" "$mpi_bin" tcsh make perl subversion 
fi
