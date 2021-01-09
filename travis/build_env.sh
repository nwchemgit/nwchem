#!/bin/bash
os=`uname`
arch=`uname -m`
dist="ubuntu"
if test -f " /usr/lib/fedora-release"; then
    dist="fedora"
fi
echo dist is "$dist"
echo DISTR is "$DISTR"
 if [[ "$os" == "Darwin" ]]; then 
#  HOMEBREW_NO_AUTO_UPDATE=1 brew cask uninstall oclint || true  
#  HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" openblas python3 ||true
     HOMEBREW_NO_INSTALL_CLEANUP=1  HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc "$MPI_IMPL" python3 ||true
     #hack to fix Github actions mpif90
     ln -sf /usr/local/bin/$FC /usr/local/bin/gfortran
     $FC --version
     gfortran --version
#  if [[ "$MPI_IMPL" == "openmpi" ]]; then
#      HOMEBREW_NO_INSTALL_CLEANUP=1 HOMEBREW_NO_AUTO_UPDATE=1 brew install scalapack
#  fi
fi
 if [[ "$os" == "Linux" ]]; then
     if [[ "$DISTR" == "fedora" ]];then
	 sudo dnf udate;  sudo dnf -y install perl perl python3-devel time patch openblas-serial64 openmpi-devel cmake gcc-gfortran unzip
	 #	 module load mpi
	 export PATH=/usr/lib64/openmpi/bin:$PATH
	 export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
	 which mpif90
	 mpif90 -show
     else
    if [[ "$MPI_IMPL" == "openmpi" ]]; then
	mpi_bin="openmpi-bin" ; mpi_libdev="libopenmpi-dev" scalapack_libdev="libscalapack-openmpi-dev"
    fi
    if [[ "$MPI_IMPL" == "mpich" ]]; then
        mpi_bin="mpich" ; mpi_libdev="libmpich-dev" scalapack_libdev="libscalapack-mpich-dev"
    fi
    if [[ "$MPI_IMPL" == "intel" ]]; then
	wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB \
            && sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB  \
            && rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB  \
            && echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list \
            && sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"  \
	    && sudo apt update
	sudo apt-cache search intel-oneapi-mpi
        mpi_bin="  " ; mpi_libdev="intel-oneapi-mpi-devel" scalapack_libdev="intel-oneapi-mkl"
    fi
    sudo add-apt-repository universe && sudo apt update
#    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin" "$scalapack_libdev"  make perl  libopenblas-dev python3 rsync
    sudo apt-get -y install gfortran python3-dev python-dev cmake "$mpi_libdev" "$mpi_bin"  make perl  python3 rsync
    if [[ "$FC" == "ifort" ]]; then
	sudo apt-get -y install intel-oneapi-ifort intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic  intel-oneapi-mkl
	sudo apt-get -y install intel-oneapi-mpi-devel
    fi
    fi
fi
