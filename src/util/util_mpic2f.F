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
      integer function util_mpicommrank()
      USE ISO_C_BINDING
      implicit none
      interface
        subroutine utilc_mpicommrank(rank) bind(C)
          use, intrinsic                  :: ISO_C_BINDING
          integer (C_INT)                 :: rank
        end subroutine
      end interface
      integer (C_INT) :: rank
      call utilc_mpicommrank(rank)
      util_mpicommrank=rank
      end
      integer function util_mpicommsize()
      USE ISO_C_BINDING
      implicit none
      interface
        subroutine utilc_mpicommsize(size) bind(C)
          use, intrinsic                  :: ISO_C_BINDING
          integer (C_INT)                 :: size
        end subroutine
      end interface
      integer (C_INT) :: size
      call utilc_mpicommsize(size)
      util_mpicommsize=size
      end
c $Id$
