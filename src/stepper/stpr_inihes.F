C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_inihes ( HESIAN,rtdb )
c $Id$
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
#include "errquit.fh"
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
#include "stdio.fh"
#include "chow.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stpr_clean.fh"
#include "util.fh"
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION HESIAN(NAT3TR)
      logical oprint_inithess
#if defined(WEIGHT_HESS)
      integer indx, indx_x, indx_y, indx_z
      integer nactive,nat,ma_type
      integer h_actlist, k_actlist
      double precision valinf
#endif
c
* zeroed matrix
      call dfill (NAT3TR,0.0d00, HESIAN,1)
* print initial hessian?
      oprint_inithess = util_print('inithess',print_debug)
* check for active atom list
#if defined(WEIGHT_HESS)
      if (rtdb_ma_get(rtdb,'geometry:actlist', ma_type, 
     $     nactive, h_actlist)) then 
       if (MWC) then
          write(luout,*)' Unknown how to fill hessian for NWC ',
     &        'stepper walking'
          call errquit('stpr_inihes: fatal logical error',911,
     &       UNKNOWN_ERR)
        endif
        nat = nat3/3
        if (nactive.le.0 .or. nactive.gt.nat) call errquit
     &      ('stpr_inihes: invalid number of active atoms',nactive,
     &       INPUT_ERR)
        if (.not.MA_get_index(h_actlist,k_actlist)) call errquit
     &      ('stpr_inihes: ma_get_index failed for actlist',911, MA_ERR)
* have active list of atoms now fill matrix with infinite values        
        valinf = 300.0d00/thrcln
        do i=1,nat3
          indx = i*(i-1)/2 + i
          hesian(indx) = valinf
        enddo
* modify hesian such that values for active atoms are 1.0
        do j = 0,(nactive-1)
          indx = int_mb(k_actlist+j)   ! = atom number
          indx = 3*(indx-1) + 1        ! = first index in nat3
          indx_x = indx*(indx-1)/2 + indx
          indx = indx + 1
          indx_y = indx*(indx-1)/2 + indx
          indx = indx + 1
          indx_z = indx*(indx-1)/2 + indx
          hesian(indx_x) = ONE
          hesian(indx_y) = ONE
          hesian(indx_z) = ONE
        enddo
        if (.not.MA_free_heap(h_actlist)) call errquit
     &      ('stpr_inihes: ma_free_heap failed for actlist',911, MA_ERR)
        WRITE(luout,*)'  Nuclear Hessian initialized with UNIT',
     &      ' matrix for active atoms only.'
        write(luout,*)'  Inactive atoms given pseudo-infinite value.'
      else
#endif
*---
*...... i   : 1  2  3  4  5  6
*...... indx: 1  3  6 10 15 21
*---
        if (MWC) then
          do i = 1,nat3
            indx = i*(i-1)/2 + i
            hesian(indx) = 4.0d-04
          enddo
        else
          do i = 1,nat3
            indx = i*(i-1)/2 + i
            hesian(indx) = one
          enddo
        endif
        WRITE(luout,*)'  Nuclear Hessian initialized with UNIT',
     &      ' matrix. '
#if defined(WEIGHT_HESS)
      endif
#endif
      if (oprint_inithess)
     &    call stpr_prntpd(hesian,nat3tr,nat3,luout)
#if defined(OLD_STPR_CODE)
      IJ = 0
      DO 2 I = 1,NAT3
        DO 1 J = 1,I
          IJ = IJ + 1
          IF (I.EQ.J) THEN
            IF (MWC) THEN
              HESIAN(IJ) = 4.0D-04
            ELSE
              HESIAN(IJ) = ONE
            ENDIF
          ELSE
            HESIAN(IJ) = ZERO
          ENDIF
    1   CONTINUE
    2 CONTINUE
#endif
*debug:      write(6,*)'hesian stpr_inihes', nat3tr
*debug:      do i=1,nat3tr
*debug:        if (abs(hesian(i)).gt.1.0d-5)
*debug:     &      write(6,*)' hesian (',i,') =',hesian(i)
*debug:      enddo
      RETURN
      END
C> @}
