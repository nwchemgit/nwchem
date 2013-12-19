C> \brief The Fortran interface routine to get the NWCHEM source
C> directory
C>
C> In order to compile the source code location into the code we have
C> to choose between two evils: 
C> - Either we struggle with the random Fortran line length limits
C>   (because they were such a disaster in Fortran77 compiler developers
C>   have kindly carried these over into Fortran90 and further :-(.
C> - Or we have to generate the code is C and use to ISO_C_BINDING
C>   module to integrate this into Fortran.
C> Either approach is a bit of a mess, so we encapsulate all the messing
C> around in this routine. This routine then provides a straightforward
C> interface that can simply be called from the Fortran code without
C> having to pull any stunts.
C>
      subroutine util_nwchem_srcdir(pathname)
#ifdef NWCHEM_LONG_PATHS
      USE ISO_C_BINDING
#endif
      implicit none
#ifdef NWCHEM_LONG_PATHS
#include "utilc_nwchem_srcdir.fh"
      character (KIND=C_CHAR,LEN=256) :: compiled_name
      integer (C_INT) :: length
#else
      character*256 compiled_name
      integer length
#endif
      character*(*) pathname !< [Output] The compiled in pathname
#ifdef NWCHEM_LONG_PATHS
      length = len(compiled_name)
      call utilc_nwchem_srcdir(compiled_name,length)
#else
      compiled_name =
     &NWCHEM_SRCDIR
#endif
      pathname = compiled_name
      end
c $Id$
