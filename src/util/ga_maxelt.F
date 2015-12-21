      subroutine ga_maxelt(g_a, value)
C     $Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "errquit.fh"
      integer g_a
      double precision value
c     
c     Determine the absolute value of the element with largest 
c     absolute magnitude
c     
      integer ilo, ihi, jlo, jhi
      integer itype, ndim, dims(ga_max_dim)
      integer ga_maxelt_dbl_doit, ga_maxelt_dcpl_doit
      external ga_maxelt_dbl_doit, ga_maxelt_dcpl_doit
c     
c     Each process first determines max of data local to it
c     
      call ga_sync()
      call nga_inquire(g_a,itype,ndim,dims)
      if (ndim.ne.2) call errquit("ga_maxelt: invalid no. dims",
     +                            ndim,UERR)
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      value = 0.0d0    ! for nodes with no data
      if (ihi.gt.0 .and. ihi.ge.ilo) then
        if (itype.eq.MT_DBL) then
          call ga_access_callback_release(
     $         g_a, ilo, ihi, jlo, jhi,ga_maxelt_dbl_doit, value, 
     $         0, 0, 0, 0, 0,0)
        else if (itype.eq.MT_DCPL) then
          call ga_access_callback_release(
     $         g_a, ilo, ihi, jlo, jhi,ga_maxelt_dcpl_doit, value, 
     $         0, 0, 0, 0, 0,0)
        else
          call errquit("ga_maxelt: invalid data type",itype,UERR)
        endif
      endif
c     
c     Finally a global max via message passing
c     
      call ga_sync()
#ifdef NWCHEM_USE_GOP_ABSMAX
      call ga_dgop(Msg_MxElt, value, 1, 'absmax')
#else
      value = abs(value)
      call ga_dgop(Msg_MxElt, value, 1, 'max')
#endif
      call ga_sync()
c     
      end
      integer function ga_maxelt_dbl_doit(
     $     g_a, ilo, ihi, jlo, jhi, a, ld, value, 
     $     idum2, idum3, idum4, idum5, idum6,idum7)
c
c     Double precision version
c
      implicit none
      integer g_a, ilo, ihi, jlo, jhi, ld, 
     $     idum2, idum3, idum4, idum5, idum6,idum7
      double precision a(ilo:ilo+ld-1,jlo:jhi), value
      integer i, j
c
      value = 0.0d0
      do j = jlo, jhi
         do i = ilo, ihi
            value = max(abs(a(i,j)),value)
         end do
      end do
c
      ga_maxelt_dbl_doit = 0    ! no update necessary
c
      end
      integer function ga_maxelt_dcpl_doit(
     $     g_a, ilo, ihi, jlo, jhi, a, ld, value, 
     $     idum2, idum3, idum4, idum5, idum6,idum7)
c
c     Double complex version
c
      implicit none
      integer g_a, ilo, ihi, jlo, jhi, ld, 
     $     idum2, idum3, idum4, idum5, idum6,idum7
      double complex a(ilo:ilo+ld-1,jlo:jhi)
      double precision value
      integer i, j
c
      value = 0.0d0
      do j = jlo, jhi
         do i = ilo, ihi
            value = max(abs(a(i,j)),value)
         end do
      end do
c
      ga_maxelt_dcpl_doit = 0    ! no update necessary
c
      end
