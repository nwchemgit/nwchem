*
* $Id$
*
*======================================================================
*
* DISCLAIMER
*
* This material was prepared as an account of work sponsored by an
* agency of the United States Government.  Neither the United States
* Government nor the United States Department of Energy, nor Battelle,
* nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
* COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
* SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
* INFRINGE PRIVATELY OWNED RIGHTS.
*
* ACKNOWLEDGMENT
*
* This software and its documentation were produced with Government
* support under Contract Number DE-AC06-76RLO-1830 awarded by the United
* States Department of Energy.  The Government retains a paid-up
* non-exclusive, irrevocable worldwide license to reproduce, prepare
* derivative works, perform publicly and display publicly by or for the
* Government, including the right to distribute to other Government
* contractors.
*
*======================================================================
*
*  -- PFFT routine (version 1.0) --
*     Pacific Northwest Laboratory
*     April 5, 1995
*
*======================================================================
      SUBROUTINE WRKSPC2(JOB, N1, N2, n3,
     $     MAP1, MAP2, LIWORK, LWORK,
     $     INFO )
	implicit none
      INTEGER           JOB, N1, N2, n3, LIWORK, LWORK, INFO
      INTEGER           MAP1(n3), MAP2(n2)
*-----------------------------------------------------------------------
*     Compute workspace requirements of various fft routines:
c     
*     JOB = 1, Workspace for BLKTRN
c     complex to complex , itype is used in the call to 3d fft
c     
c     m = number of panels
c     n = number of 3rd dim
c     
c     returns
c     liwork = length of integer work array
c     lwork = length of double precision work array
c     info
c
c
*-----------------------------------------------------------------------
      integer nelebk
      INTEGER           NAPROC, ME, N1MAX, N2MAX, I
      INTEGER           LSTCNT
      integer           max
      INTEGER           MXNPRC, MXMYND
      integer m2, mmax,nmax, iii, nele
      intrinsic         max
c     
      EXTERNAL          LSTCNT
      EXTERNAL          MXNPRC, MXMYND
*-----------------------------------------------------------------------
c     
      NAPROC = MXNPRC()
      ME     = MXMYND()
      INFO = 0
c     
c     computes the large number of panels that is stored on map1(1:n3) and
c     map2(1:n2)
c     

      N1MAX = 0
      DO 10 I = 0, NAPROC-1
         iii = LSTCNT( I, MAP1, n3 )
         N1MAX = MAX( N1MAX, iii)
 10   CONTINUE
      
      N2MAX = 0
      DO 20 I = 0, NAPROC-1
         iii = lstcnt(i, map2, n2)
         N2MAX = MAX( N2MAX, iii)
 20   CONTINUE
c
c     
      nelebk = 2*n1          ! satisfy compiler warnings
      if ( job .eq. 1 ) then
         nelebk = 2*n1
         m2 = n2
      endif
      if ( job .eq. 2 ) then
         nelebk = 2*(n1/2 + 1)
         if ( n1 .eq. 1 ) then
            m2 = n2/2 + 1
         else
            m2 = n2
         endif
      endif
c     
c     
*     BLKTRN
*     ======
c     
      LIWORK = MAX( N1MAX, N2MAX )
      LWORK = NELEBK * N1MAX * LSTCNT( ME, MAP2, N2 ) + 10 +
     $     NELEBK * N2MAX * LSTCNT( ME, MAP1, n3 ) + 10
      
c     IF( JOB .EQ. 1 )
c     $     RETURN
c     
c     
      mmax = 0
      nmax = 0
      do iii = 0, naproc-1
         mmax = max(mmax, lstcnt(iii, map1, n3))
      enddo
      do iii = 0, naproc-1
         nmax = max(nmax, lstcnt(iii, map2, n2))
      enddo
c     
      nele = nelebk * max(n2*mmax, n3*nmax)
c     
      lwork = 5*max(nele, lwork)/2
      liwork = 3*max(2*liwork, 3*naproc)
c     
      RETURN
      END
