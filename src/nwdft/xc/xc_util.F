c
c     xc utility functions
c
      logical function xc_chkgrad()
*
* $Id$
*
      implicit none
c
#include "cdft.fh"
#include "steric.fh"
c
c     see if any Gradient Corrected Functional is not null
c
      double precision tot,eps
      parameter(eps=1.d-10)
c
      tot = xfac(3) + xfac(4) + xfac(5) + xfac(6) + xfac(7) +
     + xfac(8)  + xfac(9)  + xfac(10) + xfac(11) +xfac(12)+ xfac(13)+
     + xfac(14) + xfac(15) + xfac(16)+xfac(17)+xfac(18) +xfac(19) +
     + xfac(22) +
     + xfac(30) + xfac(31) + xfac(33) + xfac(39) + cfac(29) +
     + cfac(2)  + cfac(4)  + cfac(5)  + cfac(12) +
     + cfac(13) + cfac(14) + cfac(15) + cfac(16) + cfac(17)+cfac(18)+
     + cfac(19) + cfac(20) + cfac(21) + cfac(22) + cfac(23) + cfac(24) +
     + xfac(25) + cfac(25) + cfac(26) + xfac(20) + xfac(21) + cfac(27) +
     + xfac(26) + xfac(27) + xfac(28) + xfac(29) + cfac(31) + xfac(32) +
     + cfac(32) + cfac(33) + cfac(34) + cfac(35) + cfac(36) +
     + xfac(34) + xfac(35) + xfac(36) + xfac(37) +
     + cfac(30) + cfac(38) + cfac(39) + cfac(40) + cfac(41) + cfac(42) +
     + xfac(46) + cfac(46) + xfac(47) + cfac(47) + xfac(40) + xfac(41) +
     + xfac(48) + xfac(49) + xfac(50) + xfac(51) +
     + cfac(48) + cfac(49) + cfac(50) + cfac(51) +
     + xfac(52) + xfac(53) + xfac(54) +
     + cfac(52) + cfac(53) + cfac(54)
c
      if (abs(tot).gt.eps.or.lb94.or.cs00)then
         xc_chkgrad = .true.
      else if (lsteric) then ! turn on gradients for steric calculations
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      end if
c
      return
      end
c
      logical function xc_gotxc()
      implicit none
#include "cdft.fh"
      double precision tot
      integer i
c
c     return true if there is any XC functional to compute 
c     (other than HF Exch)
c
      tot = -abs(xfac(1))
      do i = 1, numfunc
         tot = tot + abs(xfac(i)) + abs(cfac(i))
      enddo
      xc_gotxc = tot.gt.1d-8.or.lb94.or.cs00
      return
      end
c
      integer function xc_getnqmx()
      implicit none
#include "cdft.fh"
      xc_getnqmx=nqmax
      return
      end
c
      integer function xc_getipol()
      implicit none
#include "cdft.fh"
      xc_getipol=ipol
      return
      end
c
c     return true for functional without 2nds
c
      logical function xc_got2nd()
      implicit none
c
#include "cdft.fh"
c
      xc_got2nd=.not.(
c
     .     xfac(4).ne.0d0.or.
     .     xfac(5).ne.0d0.or.
     .     xfac(6).ne.0d0.or.
     .     xfac(10).ne.0d0.or.
     .     xfac(11).ne.0d0.or.
     .     xfac(12).ne.0d0.or.
     .     xfac(13).ne.0d0.or.
     .     xfac(14).ne.0d0.or.
     .     xfac(15).ne.0d0.or.
     .     xfac(16).ne.0d0.or.  ! optx
     .     xfac(17).ne.0d0.or.
     .     xfac(18).ne.0d0.or.
     .     xfac(21).ne.0d0.or.
     .     xfac(22).ne.0d0.or.
     .     xfac(28).ne.0d0.or.
     .     xfac(29).ne.0d0.or.
     .     xfac(32).ne.0d0.or.
     .     xfac(33).ne.0d0.or.
     .     xfac(34).ne.0d0.or.
     .     xfac(35).ne.0d0.or.
     .     xfac(36).ne.0d0.or.
     .     xfac(37).ne.0d0.or.
     .     xfac(48).ne.0d0.or.
     .     xfac(49).ne.0d0.or.
     .     xfac(50).ne.0d0.or.
     .     xfac(51).ne.0d0.or.
     .     xfac(52).ne.0d0.or.
     .     xfac(53).ne.0d0.or.
     .     xfac(54).ne.0d0.or.
c
     .     cfac(13).ne.0d0.or.
     .     cfac(14).ne.0d0.or.
     .     cfac(15).ne.0d0.or.
     .     cfac(16).ne.0d0.or.
     .     cfac(17).ne.0d0.or.
     .     cfac(18).ne.0d0.or.
     .     cfac(19).ne.0d0.or.
     .     cfac(20).ne.0d0.or.
     .     cfac(21).ne.0d0.or.
     .     cfac(22).ne.0d0.or.
     .     cfac(25).ne.0d0.or.
     .     cfac(27).ne.0d0.or.
     .     cfac(29).ne.0d0.or.
     .     cfac(30).ne.0d0.or.
     .     cfac(31).ne.0d0.or.
     .     cfac(32).ne.0d0.or.
     .     cfac(33).ne.0d0.or.
     .     cfac(34).ne.0d0.or.
     .     cfac(35).ne.0d0.or.
     .     cfac(38).ne.0d0.or.
     .     cfac(39).ne.0d0.or.
     .     cfac(40).ne.0d0.or.
     .     cfac(42).ne.0d0.or.
     .     cfac(48).ne.0d0.or.
     .     cfac(49).ne.0d0.or.
     .     cfac(50).ne.0d0.or.
     .     cfac(51).ne.0d0.or.
     .     cfac(52).ne.0d0.or.
     .     cfac(53).ne.0d0.or.
     .     cfac(54).ne.0d0.or.
     .     cfac(36).ne.0d0)
c
      return
      end
c
      logical function xc_chktau()
C
C $Id$ 
C
      implicit none
 
c
c     see if any metaGGA Functional using TAU is not null
c
#include "cdft.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
c 
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25) + xfac(21) + cfac(27) 
     +       + xfac(28) + xfac(29) + xfac(33) + xfac(34) + xfac(35)
     +       + xfac(36) + cfac(29) + cfac(31) + cfac(32) + cfac(33) 
     +       + cfac(34) + cfac(35) + cfac(30) + cfac(38)
     +       + cfac(39) + cfac(40) + xfac(32) + cfac(42)
     +       + cfac(48) + cfac(49) + cfac(50) + cfac(51)
c
      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      return
      end
c
      logical function xc_getxcfac(rtdb)
      implicit none
      integer rtdb
c     load cfac and xfac into cdft.fh
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      xc_getxcfac=rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,
     &xfac)
      xc_getxcfac=xc_getxcfac.and.
     A     rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc,
     &     cfac)
      return
      end
c
c     return true for functionals with cam
c
      logical function xc_hascam()
      implicit none
c
#include "cdft.fh"
#include "case.fh"
c
      xc_hascam=(
     .     xfac(38).ne.0d0.or.
     .     xfac(40).ne.0d0.or.
     .     xfac(41).ne.0d0.or.
     .     xfac(42).ne.0d0.or.
     .     xfac(43).ne.0d0.or.
     .     xfac(50).ne.0d0.or.
     .     xfac(51).ne.0d0)
c
      xc_hascam = xc_hascam.or.cam_exch
c
      return
      end
c
c     return true for functionals with short-range HF
c
      logical function xc_hassrhf()
      implicit none
c
#include "cdft.fh"
#include "case.fh"
c
      xc_hassrhf = cam_srhf
c
      return
      end
c
c     set cam parameters
c
      subroutine xc_setcamparam(rtdb,dosrhf,omega,alpha,beta)
      implicit none
c
      integer rtdb
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      logical cam_exch,dosrhf
      double precision omega,alpha,beta
c
c     set cam parameters
      cam_exch = .true.    
c
c     set cam flag 
      if (.not.rtdb_put(rtdb,'dft:cam_exch', mt_log, 1, cam_exch))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set short-range HF (srhf) flag 
      if (.not.rtdb_put(rtdb,'dft:cam_srhf', mt_log, 1, dosrhf))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_omega
      if (.not.rtdb_put(rtdb,'dft:cam_omega',mt_dbl,1,omega))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_alpha
      if (.not.rtdb_put(rtdb,'dft:cam_alpha',mt_dbl,1,alpha))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_beta
      if (.not.rtdb_put(rtdb,'dft:cam_beta',mt_dbl,1,beta))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     Force texas integrals to false for CAM calculations
      if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &   call errquit('cannot set int:cando_txs',0, RTDB_ERR)
c
      return
      end
