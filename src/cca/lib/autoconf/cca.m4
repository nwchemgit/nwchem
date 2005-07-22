dnl AC_CHECK_CCA()
dnl detects the cca tools environment
dnl author: Joseph Kenny, jpkenny@sandia.gov
dnl
dnl predefined variables:
dnl   CCAFE_CONFIG        full path to ccafe-config (optional)
dnl   ENABLE_PYTHON       enable python bindings (yes/no)
dnl
dnl output variables:
dnl   CCAFE_CONFIG CCAFE_INCLUDE CCAFE_LIB CCAFE_SHARE CCAFE_BIN
dnl   CCAFE_MPI_ENABLE CCAFE_MPI_INCLUDE CCAFE_MPI_LIB CCAFE_MPI_BIN
dnl   CCA_SPEC_BABEL_CONFIG CCA_SPEC_BABEL_INCLUDE
dnl   CCA_SPEC_BABEL_LIB CCA_SPEC_BABEL_SHARE
dnl   CCA_SPEC_CLASSIC_CONFIG CCA_SPEC_CLASSIC_INCLUDE
dnl   CCA_SPEC_CLASSIC_LIB CCA_SPEC_CLASSIC_SHARE
dnl   BABEL_CONFIG BABEL_INCLUDE BABEL_LIB BABEL_SHARE BABEL_BIN
dnl   BABEL_CC BABEL_CFLAGS BABEL_CXX BABEL_CXXFLAGS
dnl   BABEL_LIBTOOL
dnl   BABEL_PYTHON_ENABLE CCAFE_PYTHON_ENABLE
dnl   BABEL_PYTHON BABEL_PYTHON_VERSION
dnl   BABEL_PYTHON_LIB BABEL_PYTHON_INCLUDE

AC_DEFUN([AC_CHECK_CCA],[

  # ccaffeine gives us everything else
  AC_ARG_WITH(ccafe-config,
  [  --with-ccafe-config     path to the ccafe-config script.],
  [ CCAFE_CONFIG=$withval ],
  [ if test -z $CCAFE_CONFIG || test ! -x $CCAFE_CONFIG; then
      AC_PATH_PROG(CCAFE_CONFIG,ccafe-config,"not-found")
    fi
  ]
  )
  if ! test -x $CCAFE_CONFIG; then
    AC_MSG_ERROR([ccaffeine not found, use --with-ccafe-config])
  fi
  CCAFE_INCLUDE=`$CCAFE_CONFIG --var CCAFE_pkgincludedir`
  CCAFE_LIB=`$CCAFE_CONFIG --var CCAFE_pkglibdir`
  CCAFE_SHARE=`$CCAFE_CONFIG --var CCAFE_pkgdatadir`
  CCAFE_BIN=`$CCAFE_CONFIG --var CCAFE_bindir`
  AC_SUBST(CCAFE_CONFIG)
  AC_SUBST(CCAFE_INCLUDE)
  AC_SUBST(CCAFE_LIB)
  AC_SUBST(CCAFE_SHARE)
  AC_SUBST(CCAFE_BIN)

  # check for cca-spec-babel
  CCA_SPEC_BABEL_CONFIG=`$CCAFE_CONFIG --var CCAFE_CCA_SPEC_BABEL_CONFIG`
  if test -z $CCA_SPEC_BABEL_CONFIG || ! test -x $CCA_SPEC_BABEL_CONFIG; then
    AC_MSG_ERROR([can't find cca-spec-babel-config])
  fi
  CCA_SPEC_BABEL_INCLUDE=`$CCA_SPEC_BABEL_CONFIG --var CCASPEC_pkgincludedir`
  CCA_SPEC_BABEL_LIB=`$CCA_SPEC_BABEL_CONFIG --var CCASPEC_pkglibdir`
  CCA_SPEC_BABEL_SHARE=`$CCA_SPEC_BABEL_CONFIG --var CCASPEC_pkgdatadir`
  AC_SUBST(CCA_SPEC_BABEL_CONFIG)
  AC_SUBST(CCA_SPEC_BABEL_INCLUDE)
  AC_SUBST(CCA_SPEC_BABEL_LIB)
  AC_SUBST(CCA_SPEC_BABEL_SHARE)

  # check for cca-spec-classic
  CCA_SPEC_CLASSIC_ROOT=`$CCAFE_CONFIG --var CCAFE_CLASSIC_CCA_ROOT`
  CCA_SPEC_CLASSIC_CONFIG="$CCA_SPEC_CLASSIC_ROOT/bin/cca-spec-classic-config"
  if test -z $CCA_SPEC_CLASSIC_CONFIG || test ! -e $CCA_SPEC_CLASSIC_CONFIG; then
    AC_MSG_ERROR([can't find cca-spec-classic-config])
  fi
  CCA_SPEC_CLASSIC_INCLUDE=`$CCA_SPEC_CLASSIC_CONFIG --var CLASSIC_pkgincludedir`
  CCA_SPEC_CLASSIC_LIB=`$CCA_SPEC_CLASSIC_CONFIG --var CLASSIC_pkglibdir`
  CCA_SPEC_CLASSIC_SHARE=`$CCA_SPEC_CLASSIC_CONFIG --var CLASSIC_pkgdatadir`
  AC_SUBST(CCA_SPEC_CLASSIC_CONFIG)
  AC_SUBST(CCA_SPEC_CLASSIC_INCLUDE)
  AC_SUBST(CCA_SPEC_CLASSIC_LIB)
  AC_SUBST(CCA_SPEC_CLASSIC_SHARE)

  # check for babel
  BABEL_CONFIG=`$CCA_SPEC_BABEL_CONFIG --var CCASPEC_BABEL_BABEL_CONFIG`
  if test ! -x $BABEL_CONFIG; then
    AC_MSG_ERROR([can't find babel-config])
  fi
  BABEL_INCLUDE=`$BABEL_CONFIG --includedir`
  BABEL_LIB=`$BABEL_CONFIG --libdir`
  BABEL_SHARE=`$BABEL_CONFIG --datadir`
  BABEL_BIN=`$BABEL_CONFIG --bindir`
  AC_SUBST(BABEL_CONFIG)
  AC_SUBST(BABEL_INCLUDE)
  AC_SUBST(BABEL_LIB)
  AC_SUBST(BABEL_SHARE)
  AC_SUBST(BABEL_BIN)

  # check for babel compilers
  BABEL_CC=`$BABEL_CONFIG --query-var=CC`
  BABEL_CFLAGS=`$BABEL_CONFIG --query-var=CFLAGS`
  BABEL_CXX=`$BABEL_CONFIG --query-var=CXX`
  BABEL_CXXFLAGS=`$BABEL_CONFIG --query-var=CXXFLAGS`
  AC_SUBST(BABEL_CC)
  AC_SUBST(BABEL_CFLAGS)
  AC_SUBST(BABEL_CXX)
  AC_SUBST(BABEL_CXXFLAGS)

  # might as well use babel's libtool
  BABEL_BIN=`$BABEL_CONFIG --bindir`
  BABEL_LIBTOOL=$BABEL_BIN/babel-libtool
  if test -z $BABEL_LIBTOOL || ! test -x $BABEL_LIBTOOL; then
    AC_MSG_ERROR([can't find babel-libtool])
  fi
  AC_SUBST(BABEL_LIBTOOL)

  # check mpi configuration
  CCAFE_MPI_INCLUDE=`$CCAFE_CONFIG --var CCAFE_MPI_INC`
  CCAFE_MPI_LIB=`$CCAFE_CONFIG --var CCAFE_MPI_LIBDIR`
  CCAFE_MPI_BIN=`$CCAFE_CONFIG --var CCAFE_MPI_BIN`
  if test -z "$CCAFE_MPI_INCLUDE"; then
    CCAFE_MPI_ENABLE="no"
    AC_MSG_WARN([Ccaffeine not configured for MPI])
  else
   CCAFE_MPI_ENABLE="yes"
   CCAFE_MPI_INCLUDE=`echo $CCAFE_MPI_INCLUDE | sed 's/^\-I//'`
  fi
  AC_SUBST(CCAFE_MPI_ENABLE)
  AC_SUBST(CCAFE_MPI_INCLUDE)
  AC_SUBST(CCAFE_MPI_LIB)
  AC_SUBST(CCAFE_MPI_BIN)

  if test $ENABLE_PYTHON == "yes"; then
    # check for babel python
    BABEL_PYTHON_ENABLE=`$BABEL_CONFIG --query-var=SUPPORT_PYTHON`
    if test $BABEL_PYTHON_ENABLE == "false"; then
      AC_MSG_ERROR([Babel not properly configured for python])
    fi
    # check that ccafe is configured for python
    if ! test -d $CCAFE_ROOT/lib/python$PYTHON_VERSION/site-packages/ccaffeine; then
      AC_MSG_ERROR([Ccaffeine not properly configured for Python])
    else
      CCAFE_PYTHON_ENABLE="yes"
    fi
    BABEL_PYTHON=`$BABEL_CONFIG --query-var=WHICH_PYTHON`
    BABEL_PYTHON_VERSION=`$BABEL_CONFIG --query-var=PYTHON_VERSION`
    BABEL_PYTHON_LIB=`$BABEL_CONFIG --query-var=PYTHONLIB`/site-packages
    BABEL_PYTHON_INCLUDE=`$BABEL_CONFIG --query-var=PYTHONINC`
    AC_SUBST(BABEL_PYTHON)
    AC_SUBST(BABEL_PYTHON_VERSION)
    AC_SUBST(BABEL_PYTHON_LIB)
    AC_SUBST(BABEL_PYTHON_INCLUDE)
  else
    BABEL_PYTHON_ENABLE="no"
    CCAFE_PYTHON_ENABLE="no"
  fi

  echo -e "\nCCA Tools Configuration:"
  echo -e "---------------------------------------------------------------"
  echo -e "ccafe config:\n  $CCAFE_CONFIG"
  echo -e "ccafe include:\n  $CCAFE_INCLUDE"
  echo -e "ccafe lib:\n  $CCAFE_LIB"
  echo -e "ccafe share:\n  $CCAFE_SHARE"
  echo -e "ccafe bin:\n  $CCAFE_BIN"
  echo -e "ccafe python enabled:\n  $CCAFE_PYTHON_ENABLE"
  echo -e "ccafe mpi enabled\n  $CCAFE_MPI_ENABLE"
  if test $CCAFE_MPI_ENABLE == "yes"; then
    echo -e "ccafe mpi include:\n  $CCAFE_MPI_INCLUDE"
    echo -e "ccafe mpi lib:\n  $CCAFE_MPI_LIB"
    echo -e "ccafe mpi bin:\n  $CCAFE_MPI_BIN"
  fi
  echo -e "cca-spec-babel-config:\n $CCA_SPEC_BABEL_CONFIG"
  echo -e "cca-spec-babel include:\n  $CCA_SPEC_BABEL_INCLUDE"
  echo -e "cca-spec-babel lib:\n  $CCA_SPEC_BABEL_LIB"
  echo -e "cca-spec-babel share:\n  $CCA_SPEC_BABEL_SHARE"
  echo -e "cca-spec-classic-config:\n  $CCA_SPEC_CLASSIC_CONFIG"
  echo -e "cca-spec-classic include:\n  $CCA_SPEC_CLASSIC_INCLUDE"
  echo -e "cca-spec-classic lib:\n  $CCA_SPEC_CLASSIC_LIB"
  echo -e "cca-spec-classic share:\n  $CCA_SPEC_CLASSIC_SHARE"
  echo -e "babel-config:\n  $BABEL_CONFIG"
  echo -e "babel include:\n  $BABEL_INCLUDE"
  echo -e "babel lib:\n  $BABEL_LIB"
  echo -e "babel share:\n  $BABEL_SHARE"
  echo -e "babel bin:\n  $BABEL_BIN"
  echo -e "babel C compiler:\n  $BABEL_CC"
  echo -e "babel C++ compiler:\n  $BABEL_CXX"
  echo -e "babel CFLAGS:\n  $BABEL_CFLAGS"
  echo -e "babel CXXFLAGS:\n  $BABEL_CXXFLAGS"
  echo -e "babel libtool:\n  $BABEL_LIBTOOL"
  echo -e "babel python enabled:\n  $BABEL_PYTHON_ENABLE\n"

])

