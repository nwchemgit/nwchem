# NWChem SPEC file
# ================
#
# This spec file was created for a tutorial at Imperial College London. On
# the training machines at the site no compilers nor MPI commands were available
# requiring the generation of an RPM that installs mpirun (from OpenMPI) and
# nwchem (the patches NWChem 6.5 release). Advice was obtained from David Brown
# of EMSL operations team, in addition information was found online at:
#
#   http://www.rpm.org/max-rpm/ch-rpm-build.html
#   http://www.ibm.com/developerworks/library/l-rpm1/
#
# The process starts by creating the following directory structure:
#
#   ~/rpmbuild              -- the top-level directory for the rpmbuild tool
#   ~/rpmbuild/BUILD        -- the directory where stuff is compiled
#   ~/rpmbuild/BUILDROOT    -- the RPM contents stage area
#   ~/rpmbuild/RPMS         -- the directory for the binary .rpm files
#   ~/rpmbuild/SOURCES      -- the directory for source code .tar.gz files
#   ~/rpmbuild/SPECS        -- the RPM specifications directory 
#   ~/rpmbuild/SRPMS        -- the directory for the sources .rpm files
#
# The SPECS directory is where this file and other files like it live. 
# A .spec file contains a number of sections most of which are labeled 
# %<name> where <name> indicates a particular function/purpose. The first
# section is the preamble which does not have a particular name. The preamble
# sets a number of variables that say something about the RPM, things like
# the name of the code, the version number, the names of the source code tar
# files, etc.
#
# The labeled sections start of with %description which provides a few lines
# of information about the code.
#
# %prep lists the commands that are needed to put the source code in the right
# place. The %setup macro can automate most of this. The first %setup is applied
# to Source0 from the preamble, the next %setup to Source1, etc.
#
# %build contains a shell script (using Bash) that compiles the code. The
# script launches in the directory containing the source code from Source0. In
# this particular case we needed to compile OpenMPI (with a specific --prefix)
#first so we could build NWChem with MPI, then we build NWChem, and then we
# build OpenMPI again but this time to generate an mpirun executable that will
# be installed on the target machine (i.e. with the default --prefix). 
#
# %install is the tricky part as this script installs the RPM on the target
# machine, but is also used to collect all the files that will be rolled up
# into the RPM. When the RPM is generated this script is launched in the
# directory containing Source0, on the target machine it is launched in the
# unpacked RPM directory. The way I got this to work was to have the %build
# script create a structure of files within the nwchem-6.5 build directory
# that matches the installation directory structure on the target machine.
# Note that $RPM_BUILD_ROOT/usr/local are replaced by the rpmbuild and rpm
# commands with locations that are particular to the RPM generation and RPM
# installation process. Potentially these replacements are different at the
# installation time depending on user specified flags on the rpm command line.
# Unfortunately the publically available documentation is extremely sketchy
# on what the expectations for file names are during RPM generation and RPM
# installation (most documentation simply states "make install" as the
# install script which is as good as useless). Hence most of what I did here is
# based on guess work.
#
# %files must contain a list of all files (or directories) that will be included
# in the RPM. If a file is not listed it will not be included.
#
# To generate the RPM run:
#
#   cd ~/rpmbuild
#   rpmbuild -ba SPECS/nwchem-6-5.spec
#
# If something did not work then running this command again will start over
# from scratch. Using the --short-circuit flag allows you to skip prior steps
# and skip directly to the -b<step> build step. E.g. rpmbuild --short-circuit
# -bi skips directly to the %install step. The rpmbuild command creates the
# .rpm file in the RPMS directory.
#
# On the target machine use:
#
#   rpm -i nwchem-6.5-1.el6.x86_64.rpm
#
# to install the code. 
#
# Hopefully these brief notes provide a leg up the next time we need to roll
# an RPM.
#
# Huub van Dam, April 01, 2015
# 
#
# === Start ===
#
# This is a spec file for compiling an RPM for NWChem-6.5 including dependencies
# like OpenMPI and Global Arrays. Of course the various components have their
# own legacy build systems. So we just want to run a bunch of commands and
# tar-up some resulting executables. The fewer assumptions are made what is 
# where the better.
#
Name:		nwchem
Version:	6.5
Release:	1%{?dist}
Summary:	Parallel Ab-Initio Quantum Chemistry Code

Group:		Applications/Scientific
License:	Various: NWChem: ECL2.0; Global Arrays: OpenSource; OpenMPI: New BSD
URL:		http://www.nwchem-sw.org
Source0:	nwchem-6.5.tar.gz
Source1:	openmpi-1.8.4.tar.gz
BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)

BuildRequires:	gcc == 4.4.7

%description
NWChem parallel ab-initio quantum chemistry code targetting computer systems. \
For more details see: http://www.nwchem-sw.org


%prep
%setup 
%setup -T -D -b 1


%build
cd ..
export BUILD_DIR=`pwd`
export MRCC_METHODS=y
export EACCSD=y
export IPCCSD=y
export BLASOPT=" "
export NWCHEM_MODULES="all python"
export ARMCI_NETWORK=MPI-TS
# Someone added a -pthread flag to gcc in mpicc which loads libpthread through
# a non-standard back door. So we now need to explicitly say what we need
# including libpthread as we cannot infer the right information from 
# "mpicc -show" as we cannot know about these non-standard compiler flags.
export LIBMPI="-L/people/d3y133/rpmbuild/BUILD/lib -lmpi_usempi -lmpi_mpifh -lmpi -lopen-rte -lopen-pal -lm -lnuma -ldl -lrt -lutil -lpthread"
export SLURMOPT=" "
env | sort
echo "Building in " `pwd`
echo "Directory structure"
ls -l
# Build OpenMPI so we can use it to compile NWChem
cd openmpi-1.8.4
./configure --without-verbs --disable-mpi-profile --disable-shared --enable-static --disable-java --disable-mca-dso --disable-vt --disable-libompitrace --without-slurm --without-loadleveler --prefix=$BUILD_DIR
make all
make install 
cd ..
export PATH=$BUILD_DIR/bin:$PATH
# Build NWChem using OpenMPI mpif90, mpicc, etc.
cd nwchem-6.5
./contrib/distro-tools/build_nwchem 
echo "=== done building nwchem ==="
cp -a bin/LINUX64/nwchem bin
mkdir -p lib/basis
cp -a src/basis/libraries lib/basis
cp -a src/nwpw/libraryps lib/basis
cp -a src/data           lib/ff
cd ..
# Rebuild OpenMPI for the target machine as we need mpirun to run NWChem
echo "current directory" `pwd`
ls -l
cd openmpi-1.8.4
make clean
./configure --disable-mpi-profile --disable-shared --enable-static --disable-java --disable-mca-dso --disable-vt --disable-libompitrace --without-slurm --without-loadleveler --without-verbs
make all
cd ..
cp -a openmpi-1.8.4/orte/tools/orterun/orterun nwchem-6.5/bin/mpirun
mkdir -p nwchem-6.5/share/openmpi
cp -a openmpi-1.8.4/orte/tools/orterun/help-orterun.txt nwchem-6.5/share/openmpi
cp -a openmpi-1.8.4/opal/runtime/help-opal-runtime.txt nwchem-6.5/share/openmpi

%install

echo "Build root     = " %{buildroot}
echo "RPM_BUILD_ROOT = " $RPM_BUILD_ROOT
pwd; env
mkdir -p $RPM_BUILD_ROOT/usr/local/bin
cp bin/nwchem $RPM_BUILD_ROOT/usr/local/bin
cp bin/mpirun $RPM_BUILD_ROOT/usr/local/bin
mkdir -p $RPM_BUILD_ROOT/usr/local/lib
cp -a lib/basis $RPM_BUILD_ROOT/usr/local/lib
cp -a lib/ff $RPM_BUILD_ROOT/usr/local/lib
mkdir -p $RPM_BUILD_ROOT/usr/local/share/openmpi
cp -a share/openmpi/help-orterun.txt $RPM_BUILD_ROOT/usr/local/share/openmpi
cp -a share/openmpi/help-opal-runtime.txt $RPM_BUILD_ROOT/usr/local/share/openmpi
mkdir -p $RPM_BUILD_ROOT/etc
cp -a etc/nwchemrc $RPM_BUILD_ROOT/etc/nwchemrc

%clean
rm -rf %{buildroot}


%files
%defattr(-,root,root,-)
/usr/local/bin/mpirun
/usr/local/bin/nwchem
/usr/local/lib/basis
/usr/local/lib/ff
/usr/local/share/openmpi
%config /etc/nwchemrc

%changelog

