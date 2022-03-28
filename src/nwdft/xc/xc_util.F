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
#include "util.fh"
c
c     see if any Gradient Corrected Functional is not null
c
      logical out1,out2
      logical,external :: nwchem_libxc_family
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
     + cfac(45) + xfac(45) + 
     + xfac(43) + xfac(64) + cfac(64) +
     + xfac(46) + cfac(46) + xfac(47) + cfac(47) + xfac(40) + xfac(41) +
     + xfac(48) + xfac(49) + xfac(50) + xfac(51) +
     + cfac(48) + cfac(49) + cfac(50) + cfac(51) +
     + xfac(52) + xfac(53) + xfac(54) + xfac(55) + xfac(56) +
     + cfac(52) + cfac(53) + cfac(54) +
     + xfac(60) + xfac(61) + xfac(62) + xfac(63) + xfac(68) + xfac(69)+
     +     xfac(66) + cfac(66) + xfac(67) + cfac(67) +
     +     cfac(68)+  cfac(69) + xfac(70) + xfac(71) + cfac(71) +
     +     xfac(72) + xfac(73) + xfac(74) +
     +     xfac(75) + cfac(75) +
     w     xfac(76) + xfac(77) + cfac(76) + cfac(77) +
     w     xfac(78) + cfac(78) + xfac(79) + cfac(79) +
     x     xfac(80) + cfac(80) + xfac(81) + cfac(81) +
     n     cfac(82) + xfac(82) + cfac(83) + xfac(83) +
     c     cfac(84) + xfac(84) + cfac(85) + xfac(85) +
     L     lxdm +
cc AJL/Begin/FDE
     + xfac_fde(3) + xfac_fde(4) + xfac_fde(5) + xfac_fde(6) + 
     + xfac_fde(7) +
     + xfac_fde(8)  + xfac_fde(9)  + xfac_fde(10) + xfac_fde(11) +
     + xfac_fde(12)+ xfac_fde(13)+
     + xfac_fde(14) + xfac_fde(15) + xfac_fde(16)+xfac_fde(17)+
     + xfac_fde(18) +xfac_fde(19) +
     + xfac_fde(22) +
     + xfac_fde(30) + xfac_fde(31) + xfac_fde(33) + xfac_fde(39) + 
     + cfac_fde(29) +
     + cfac_fde(2)  + cfac_fde(4)  + cfac_fde(5)  + cfac_fde(12) +
     + cfac_fde(13) + cfac_fde(14) + cfac_fde(15) + cfac_fde(16) + 
     + cfac_fde(17)+cfac_fde(18)+
     + cfac_fde(19) + cfac_fde(20) + cfac_fde(21) + cfac_fde(22) + 
     + cfac_fde(23) + cfac_fde(24) +
     + xfac_fde(25) + cfac_fde(25) + cfac_fde(26) + xfac_fde(20) + 
     + xfac_fde(21) + cfac_fde(27) +
     + xfac_fde(26) + xfac_fde(27) + xfac_fde(28) + xfac_fde(29) + 
     + cfac_fde(31) + xfac_fde(32) +
     + cfac_fde(32) + cfac_fde(33) + cfac_fde(34) + cfac_fde(35) + 
     + cfac_fde(36) +
     + xfac_fde(34) + xfac_fde(35) + xfac_fde(36) + xfac_fde(37) +
     + cfac_fde(30) + cfac_fde(38) + cfac_fde(39) + cfac_fde(40) + 
     + cfac_fde(41) + cfac_fde(42) +
     + xfac_fde(43) +
     + xfac_fde(46) + cfac_fde(46) + xfac_fde(47) + cfac_fde(47) + 
     + xfac_fde(40) + xfac_fde(41) +
     + xfac_fde(48) + xfac_fde(49) + xfac_fde(50) + xfac_fde(51) +
     + cfac_fde(48) + cfac_fde(49) + cfac_fde(50) + cfac_fde(51) +
     + xfac_fde(52) + xfac_fde(53) + xfac_fde(54) + xfac_fde(55) + 
     + xfac_fde(56) +
     + cfac_fde(52) + cfac_fde(53) + cfac_fde(54) +
     + xfac_fde(60) + xfac_fde(61) + xfac_fde(62) + xfac_fde(63) +
     + tsfac(2)
cc AJL/End
c
      if (abs(tot).gt.eps.or.lb94.or.cs00) then
         xc_chkgrad = .true.
      else if (lsteric) then ! turn on gradients for steric calculations
         xc_chkgrad = .true.
      else
         xc_chkgrad = .false.
      end if
      if (util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_is_gga",out1)
         call nwxc_getvals("nwxc_is_mgga",out2)
         xc_chkgrad = xc_chkgrad .or. out1.or.out2
      endif
      out1 = nwchem_libxc_family("gga")
      out2 = nwchem_libxc_family("mgga")
      xc_chkgrad = xc_chkgrad .or. out1 .or. out2
c
      return
      end
c
      logical function xc_gotxc()
      implicit none
#include "cdft.fh"
#include "util.fh"
      double precision tot
      integer i
      logical out1
c
c     return true if there is any XC functional to compute 
c     (other than HF Exch)
c
cc AJL/Begin/FDE
      tot = -abs(xfac(1)) - abs(xfac_fde(1))
      do i = 1, numfunc
         tot = tot + abs(xfac(i)) + abs(cfac(i))
         tot = tot + abs(xfac_fde(i)) + abs(cfac_fde(i))
cc AJL/End
      enddo
      xc_gotxc = tot.gt.1d-8.or.lb94.or.cs00
      if (util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_has_xc",out1)
        xc_gotxc = xc_gotxc.or.out1
      endif
      if (libxcon) xc_gotxc = .true.

      return
      end
      logical function xc_gothfx()
      implicit none
#include "cdft.fh"
#include "util.fh"
c
c     return true if there is any HF exch component
      xc_gothfx = abs(xfac(1)).gt.1d-8
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
#include "util.fh"
      logical out1
      logical, external :: nwchem_libxc_family
c
      xc_got2nd=.not.(
c
cXhcth     .     xfac(4).ne.0d0.or.
cXbecke97     .     xfac(5).ne.0d0.or.
cXb971     .     xfac(6).ne.0d0.or.
cXh120     .     xfac(10).ne.0d0.or.
cXh147     .     xfac(11).ne.0d0.or.
cXB98     .     xfac(12).ne.0d0.or.
cXB97gga1     .     xfac(13).ne.0d0.or.
ch407     .     xfac(14).ne.0d0.or.
chp14     .     xfac(15).ne.0d0.or.
     .     xfac(16).ne.0d0.or.  ! optx
cXB972     .     xfac(17).ne.0d0.or.
     .     xfac(18).ne.0d0.or.
     .     xfac(21).ne.0d0.or.
cb973     .     xfac(22).ne.0d0.or.
     .     xfac(28).ne.0d0.or.
     .     xfac(29).ne.0d0.or.
     .     xfac(32).ne.0d0.or.
     .     xfac(33).ne.0d0.or.
     .     xfac(34).ne.0d0.or.
     .     xfac(35).ne.0d0.or.
     .     xfac(36).ne.0d0.or.
     .     xfac(37).ne.0d0.or.
     .     xfac(45).ne.0d0.or.
     .     xfac(48).ne.0d0.or.
     .     xfac(49).ne.0d0.or.
     .     xfac(50).ne.0d0.or.
     .     xfac(51).ne.0d0.or.
     .     xfac(52).ne.0d0.or.
     .     xfac(53).ne.0d0.or.
     .     xfac(54).ne.0d0.or.
     .     xfac(55).ne.0d0.or.
     .     xfac(56).ne.0d0.or.
     .     xfac(64).ne.0d0.or.
     .     xfac(66).ne.0d0.or.
     .     xfac(68).ne.0d0.or.
     .     xfac(69).ne.0d0.or.
     .     xfac(71).ne.0d0.or.
     .     xfac(73).ne.0d0.or.
     .     xfac(74).ne.0d0.or.
cwb97 should have 2nds
     .     xfac(76).ne.0d0.or.
     .     xfac(77).ne.0d0.or.
     .     xfac(78).ne.0d0.or.
     .     xfac(79).ne.0d0.or.
     .     xfac(80).ne.0d0.or.
     .     xfac(81).ne.0d0.or.
     .     xfac(82).ne.0d0.or.
     .     xfac(83).ne.0d0.or.
     .     xfac(84).ne.0d0.or.
     .     xfac(85).ne.0d0.or.
c
chcth     .     cfac(13).ne.0d0.or.
cCbecke97     .     cfac(14).ne.0d0.or.
cCB97-1     .     cfac(15).ne.0d0.or.
cCh120     .     cfac(16).ne.0d0.or.
cCh147     .     cfac(17).ne.0d0.or.
CCB98     .     cfac(18).ne.0d0.or.
cCB97gga1     .     cfac(19).ne.0d0.or.
cCh407     .     cfac(20).ne.0d0.or.
chp414     .     cfac(21).ne.0d0.or.
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
     .     cfac(68).ne.0d0.or.
     .     cfac(69).ne.0d0.or.
     .     cfac(71).ne.0d0.or.
     .     cfac(73).ne.0d0.or.
     .     cfac(74).ne.0d0.or.
cwb97 should have 2nds
     .     cfac(76).ne.0d0.or.
     .     cfac(77).ne.0d0.or.
     .     cfac(78).ne.0d0.or.
     .     cfac(79).ne.0d0.or.
     .     cfac(80).ne.0d0.or.
     .     cfac(81).ne.0d0.or.
     .     cfac(83).ne.0d0.or.
     .     cfac(84).ne.0d0.or.
     .     cfac(85).ne.0d0.or.
     .     cfac(36).ne.0d0)
      if (util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_has_2nd",out1)
        xc_got2nd = xc_got2nd .or. out1
      endif

      if (libxcon .and. xc_got2nd) then
        xc_got2nd = xc_got2nd .and. nwchem_libxc_family("fxc")
      endif
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
#include "util.fh"
C  note in cdft.fh,they define cfac(numfunc) and xfac(numfunc)
c 
      logical out1
      logical,external :: nwchem_libxc_family
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(18) + cfac(25) + xfac(21) + cfac(27) +
     +      xfac(28) + xfac(29) + xfac(33) + xfac(34) + xfac(35) +
     +      xfac(36) + xfac(37) + xfac(48) + xfac(49) + xfac(50) +
     +      xfac(51) + lxdm + xfac(64) +
     +      cfac(29) + cfac(31) + cfac(32) + cfac(33) +
     +      cfac(34) + cfac(35) + cfac(30) + cfac(38) +
     +      cfac(39) + cfac(40) + xfac(32) + cfac(42) +
     +      cfac(48) + cfac(49) + cfac(50) + cfac(51) +
     +      xfac(66) + cfac(66) + xfac(67) + cfac(67) +
     +      xfac(68) + cfac(68) + xfac(69) + cfac(69) +
     +      xfac(71) + cfac(71) + xfac(73) + cfac(73) +
     +      xfac(79) + cfac(79) + xfac(80) + cfac(80) +
     +      xfac(81) + cfac(81) + xfac(83) + cfac(83) +
     +      xfac(84) + cfac(84) + xfac(85) + cfac(85) +
cc AJL/Begin/FDE
     +      xfac_fde(18) + cfac_fde(25) + xfac_fde(21) + cfac_fde(27) +
     +      xfac_fde(28) + xfac_fde(29) + xfac_fde(33) + xfac_fde(34) + 
     +      xfac_fde(35) +
     +      xfac_fde(36) + xfac_fde(37) + xfac_fde(48) + xfac_fde(49) + 
     +      xfac_fde(50) +
     +      xfac_fde(51) + lxdm +
     +      cfac_fde(29) + cfac_fde(31) + cfac_fde(32) + cfac_fde(33) +
     +      cfac_fde(34) + cfac_fde(35) + cfac_fde(30) + cfac_fde(38) +
     +      cfac_fde(39) + cfac_fde(40) + xfac_fde(32) + cfac_fde(42) +
     +      cfac_fde(48) + cfac_fde(49) + cfac_fde(50) + cfac_fde(51)
cc AJL/End
c
      if (abs(tot).gt.eps)then
         xc_chktau = .true.
      else
         xc_chktau = .false.
      endif
      if (util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_is_mgga",out1)
         xc_chktau = xc_chktau .or. out1
      endif
      out1 = nwchem_libxc_family("mgga")
      xc_chktau = xc_chktau .or. out1
      return
      end
c
      logical function xc_getxcfac(rtdb)
      implicit none
      integer rtdb
      logical dolibxc
c     load cfac and xfac into cdft.fh
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
      logical use_nwxc, lhashfxc
      use_nwxc=util_module_avail("nwxc")
      if (use_nwxc) then
        call nwxc_rtdb_load(rtdb,"dft",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_has_hfx",lhashfxc)
        if (lhashfxc) then
          call nwxc_getwght("nwxc_wght_hfx",xfac(1))
        endif
      endif
      xc_getxcfac=use_nwxc
      if (use_nwxc) return
c
      xc_getxcfac=rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,
     &xfac)
      xc_getxcfac=xc_getxcfac.and.
     A     rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc,
     &     cfac)
      if(.not.rtdb_get(rtdb, 'dft:libxcon', mt_log, 1, dolibxc))
     D     dolibxc=.false.
      xc_getxcfac=xc_getxcfac.or.dolibxc
      if (dolibxc) then
        call nwchem_libxc_rdinput(rtdb,'dft')
      endif

      return
      end
c
c     return true for functionals with cam
c
      logical function xc_hascam(rtdb)
      implicit none
c
#include "cdft.fh"
#include "case.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb              ! [in]
c      
      logical out1
c
      xc_hascam=(
     .     xfac(38).ne.0d0.or.
     .     xfac(40).ne.0d0.or.
     .     xfac(41).ne.0d0.or.
     .     xfac(42).ne.0d0.or.
     .     xfac(43).ne.0d0.or.
     .     xfac(50).ne.0d0.or.
     .     xfac(76).ne.0d0.or.
     .     xfac(77).ne.0d0.or.
     .     xfac(78).ne.0d0.or.
     .     xfac(79).ne.0d0.or.
     .     xfac(80).ne.0d0.or.
     .     xfac(81).ne.0d0.or.
     .     xfac(85).ne.0d0.or.
     .     xfac(51).ne.0d0)
      if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1, cam_exch)) 
     &     cam_exch=.false.
c
      xc_hascam = xc_hascam.or.cam_exch
      if (util_module_avail("nwxc")) then
         call nwxc_getvals("nwxc_has_cam",out1)
        xc_hascam = xc_hascam.or.out1
      endif
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
c     set range-separation parameters
c
      subroutine xc_setcamparam(rtdb,cam_exch,cam_srhf,
     &   cam_omega,cam_alpha,cam_beta)
c
      implicit none
c
      integer rtdb
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
c
      logical cam_exch,cam_srhf
      double precision cam_omega,cam_alpha,cam_beta
c
      logical direct
      integer memsize, filesize
c
c     set cam flag 
      if (.not.rtdb_put(rtdb,'dft:cam_exch', mt_log, 1, cam_exch))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set short-range HF (srhf) flag 
      if (.not.rtdb_put(rtdb,'dft:cam_srhf', mt_log, 1, cam_srhf))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_omega
      if (.not.rtdb_put(rtdb,'dft:cam_omega',mt_dbl,1,cam_omega))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_alpha
      if (.not.rtdb_put(rtdb,'dft:cam_alpha',mt_dbl,1,cam_alpha))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
c     set cam_beta
      if (.not.rtdb_put(rtdb,'dft:cam_beta',mt_dbl,1,cam_beta))
     & call errquit('xc_setcamparam: rtdb_put failed', 2700, RTDB_ERR)
c
      if (cam_exch) then
c
c       Force texas integrals to false for CAM calculations
        if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &     call errquit('cannot set int:cando_txs',0, RTDB_ERR)
c     unset cando_nw if defined to avoid Simint
        if (rtdb_get(rtdb,'int:cando_nw',mt_log,1,direct)) then
           if (.not.rtdb_delete(rtdb,'int:cando_nw'))
     &          call errquit('cannot del int:cando_nw',0, RTDB_ERR)
           if(ga_nodeid().eq.0) write(luout,*) ' cando_nw unset'
        endif
c
c       Force use of direct
c
        direct = .true.
        filesize = -1
        memsize  = -1
        if (.not. rtdb_put(rtdb, 'dft:direct',mt_log,1,direct))
     &     call errquit('xc_setcamparam: rtdb_put failed',400,RTDB_ERR)
        if (.not. rtdb_put(rtdb,'int2e:filesize',mt_int,1,filesize))
     $     call errquit('xc_setcamparam: failed to store filesize',401,
     &                  RTDB_ERR)
        if (.not. rtdb_put(rtdb,'int2e:memsize',mt_int,1,memsize))
     $     call errquit('xc_setcamparam: failed to store filesize',402,
     &                  RTDB_ERR)
      endif
c
      return
      end
c
c     return true for functionals without 3rds
c
      logical function xc_ngot3rd()
      implicit none
      logical,external :: nwchem_libxc_family
c
#include "cdft.fh"
c
      xc_ngot3rd = (
c
     .     xfac(4).ne.0d0.or.   ! HCTH
     .     xfac(5).ne.0d0.or.   ! HCTH (b970)
     .     xfac(6).ne.0d0.or.   ! HCTH (b971)
     .     xfac(8).ne.0d0.or.   ! Gill 96
     .     xfac(9).ne.0d0.or.   ! PW91
     .     xfac(10).ne.0d0.or.  ! HCTH (h120)
     .     xfac(11).ne.0d0.or.  ! HCTH (h147)
     .     xfac(12).ne.0d0.or.  ! HCTH (b980)
     .     xfac(13).ne.0d0.or.  ! HCTH (b97g)
     .     xfac(14).ne.0d0.or.  ! HCTH (h407)
     .     xfac(15).ne.0d0.or.  ! HCTH (hp14)
     .     xfac(16).ne.0d0.or.  ! OPTX
     .     xfac(17).ne.0d0.or.  ! HCTH (b972)
     .     xfac(18).ne.0d0.or.  ! PKZB99-EXC
     .     xfac(19).ne.0d0.or.  ! mPW91
     .     xfac(20).ne.0d0.or.  ! HCTH (407p) 
     .     xfac(21).ne.0d0.or.  ! TPSS03-EXC
     .     xfac(22).ne.0d0.or.  ! HCTH (b973)
     .     xfac(25).ne.0d0.or.  ! ft97 
     .     xfac(26).ne.0d0.or.  ! PW6B95
     .     xfac(27).ne.0d0.or.  ! PW6BK
     .     xfac(28).ne.0d0.or.  ! M05 
     .     xfac(29).ne.0d0.or.  ! M05-2X
     .     xfac(32).ne.0d0.or.  ! dlDF
     .     xfac(33).ne.0d0.or.  ! VSXC
     .     xfac(34).ne.0d0.or.  ! M06-L
     .     xfac(35).ne.0d0.or.  ! M06-HF
     .     xfac(36).ne.0d0.or.  ! M06
     .     xfac(37).ne.0d0.or.  ! M06-2X
     .     xfac(39).ne.0d0.or.  ! HCTH (b97d)
     .     xfac(45).ne.0d0.or.  ! n12
     .     xfac(46).ne.0d0.or.  ! SSB-D
     .     xfac(47).ne.0d0.or.  ! kt1
     .     xfac(48).ne.0d0.or.  ! M08-HX 
     .     xfac(49).ne.0d0.or.  ! M08-SO
     .     xfac(50).ne.0d0.or.  ! M11
     .     xfac(51).ne.0d0.or.  ! M11-L
     .     xfac(52).ne.0d0.or.  ! SOGGA
     .     xfac(53).ne.0d0.or.  ! SOGGA11
     .     xfac(54).ne.0d0.or.  ! SOGGA11-X
     .     xfac(64).ne.0d0.or.  ! MVS15-EXC
     .     xfac(66).ne.0d0.or.  ! SCAN
     .     xfac(67).ne.0d0.or.  ! SCAN-L
     .     xfac(68).ne.0d0.or.  ! revM06
     .     xfac(69).ne.0d0.or.  ! revM06-L
     .     xfac(71).ne.0d0.or.  ! rSCAN
     .     xfac(73).ne.0d0.or.  ! r^2SCAN
     .     xfac(74).ne.0d0.or.  ! r^2SCAN-L
     .     xfac(75).ne.0d0.or.  ! hcth147@tz2p
     .     xfac(76).ne.0d0.or.  ! wb97
     .     xfac(77).ne.0d0.or.  ! wb97x
     .     xfac(78).ne.0d0.or.  ! wb97x-d3
     .     xfac(79).ne.0d0.or.  ! revm11
     .     xfac(80).ne.0d0.or.  ! mn12-l
     .     xfac(81).ne.0d0.or.  ! mn12-sx
     .     xfac(82).ne.0d0.or.  ! n12-sx
     .     xfac(83).ne.0d0.or.  ! mn15-l
     .     xfac(84).ne.0d0.or.  ! mn15
     .     xfac(85).ne.0d0.or.  ! m06-sx
c
     .     cfac(5).ne.0d0.or.   ! PW91
     .     cfac(13).ne.0d0.or.  ! HCTH
     .     cfac(14).ne.0d0.or.  ! HCTH (b970)
     .     cfac(15).ne.0d0.or.  ! HCTH (b971)
     .     cfac(16).ne.0d0.or.  ! HCTH (h120)
     .     cfac(17).ne.0d0.or.  ! HCTH (h147)
     .     cfac(18).ne.0d0.or.  ! HCTH (b980)
     .     cfac(19).ne.0d0.or.  ! HCTH (b97g)
     .     cfac(20).ne.0d0.or.  ! HCTH (h407)
     .     cfac(21).ne.0d0.or.  ! HCTH (hp14)
     .     cfac(22).ne.0d0.or.  ! OPTC
     .     cfac(23).ne.0d0.or.  ! HCTH (b972)
     .     cfac(24).ne.0d0.or.  ! ft97 
     .     cfac(25).ne.0d0.or.  ! PKZB99-COR
     .     cfac(26).ne.0d0.or.  ! HCTH (407p)
     .     cfac(27).ne.0d0.or.  ! TPSS03-COR
     .     cfac(28).ne.0d0.or.  ! HCTH (b973)
     .     cfac(29).ne.0d0.or.  ! cvs98
     .     cfac(30).ne.0d0.or.  ! M06-L
     .     cfac(31).ne.0d0.or.  ! Bc95-COR
     .     cfac(32).ne.0d0.or.  ! PW6B95
     .     cfac(33).ne.0d0.or.  ! PWB6K
     .     cfac(34).ne.0d0.or.  ! M05
     .     cfac(35).ne.0d0.or.  ! M05-2X
     .     cfac(36).ne.0d0.or.  ! Be88 
     .     cfac(37).ne.0d0.or.  ! OP 
     .     cfac(38).ne.0d0.or.  ! M06-HF
     .     cfac(39).ne.0d0.or.  ! M06
     .     cfac(40).ne.0d0.or.  ! M06-2X
     .     cfac(41).ne.0d0.or.  ! HCTH (b97d)
     .     cfac(42).ne.0d0.or.  ! dlDF
     .     cfac(45).ne.0d0.or.  ! n12
     .     cfac(46).ne.0d0.or.  ! SSB-D
     .     cfac(48).ne.0d0.or.  ! M08-HX
     .     cfac(49).ne.0d0.or.  ! M08-SO
     .     cfac(50).ne.0d0.or.  ! M11
     .     cfac(51).ne.0d0.or.  ! M11-L
     .     cfac(52).ne.0d0.or.  ! SOGGA
     .     cfac(53).ne.0d0.or.  ! SOGGA11
     .     cfac(54).ne.0d0.or.  ! SOGGA11-X
     .     cfac(66).ne.0d0.or.  ! SCAN
     .     cfac(67).ne.0d0.or.  ! SCAN-L
     .     cfac(68).ne.0d0.or.  ! revM06
     .     cfac(69).ne.0d0.or.  ! revM06-L
     .     cfac(71).ne.0d0.or.  ! rSCAN
     .     cfac(73).ne.0d0.or.  ! r^2SCAN
     .     cfac(74).ne.0d0.or.  ! r^2SCAN-L
     .     cfac(75).ne.0d0.or.  ! hcth147@tz2p
     .     cfac(76).ne.0d0.or.  ! wb97
     .     cfac(77).ne.0d0.or.  ! wb97x
     .     cfac(78).ne.0d0.or.  ! wb97x-d3
     .     cfac(79).ne.0d0.or.  ! revm11
     .     cfac(80).ne.0d0.or.  ! mn12-l
     .     cfac(81).ne.0d0.or.  ! mn12-sx
     .     cfac(82).ne.0d0.or.  ! n12-sx
     .     cfac(83).ne.0d0.or.  ! mn15-l
     .     cfac(84).ne.0d0.or.  ! mn15
     .     cfac(85).ne.0d0      ! m06-sx
     .                )
c
      if (libxcon .and. .not.xc_ngot3rd) then
        xc_ngot3rd = xc_ngot3rd .and. (.not.nwchem_libxc_family("kxc"))
      endif

      return
      end
C
      logical function xc_chklap()
      implicit none
#include "cdft.fh"
#include "util.fh"
      logical out1,nwchem_libxc_family
      double precision tot,eps
      parameter(eps=1.d-10)
      external nwchem_libxc_family
      tot = xfac(67) + cfac(67) + xfac(74) + cfac(74)

      if (abs(tot).gt.eps) then
        xc_chklap = .true.
      else
        xc_chklap = .false.
      end if
      out1 = nwchem_libxc_family("lapl")
      xc_chklap = xc_chklap .or. out1
      return
      end
      logical function xc_chkdelrq()
      implicit none
#include "cdft.fh"
#include "util.fh"
      double precision tot,eps
      parameter(eps=1.d-10)
      tot = xfac(70) 

      if (abs(tot).gt.eps) then
        xc_chkdelrq = .true.
      else
        xc_chkdelrq = .false.
      end if
      return
      end

cc DMR/End

      subroutine xc_os2cs(funcxc,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     amat, cmat, nq, ipol, ex, 
     &                     qwght, ldew, func)
      implicit none
c     convert open-shell fxc to closed-shell
#include "errquit.fh"
#include "dft2drv.fh"
c      
      external funcx
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision delrho(nq,3,ipol)
      double precision qwght(nq)
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision tol_rho, pi,fact_cs
      integer r_a,r_b,dr_a,dr_b
      integer a_a,a_b,c_a,c_b
      if(ipol.eq.1) then
         r_a=1
         r_b=1
         dr_a=1
         dr_b=1
         a_a=D1_RA
         a_b=D1_RA
         c_a=D1_GAA
         c_b=D1_GAA
         fact_cs=2d0
         call dscal(nq,0.5d0,rho,1)
         call dscal(nq*3,0.5d0,delrho,1)
      else
         r_a=2
         r_b=3
         dr_a=1
         dr_b=2
         a_a=D1_RA
         a_b=D1_RB
         c_a=D1_GAA
         c_b=D1_GBB
         fact_cs=1d0
      endif
c     alpha
      call funcxc(tol_rho, fac,lfac,nlfac,
     C     rho(1,r_a), delrho(1,1,dr_a), 
     &     amat(1,a_a), cmat(1,c_a), nq,  ex, 
     &     qwght, ldew, func, fact_cs)
c     beta
      if (ipol.eq.2) then
      call funcxc(tol_rho, fac,lfac,nlfac,
     C     rho(1,r_b), delrho(1,1,dr_b), 
     &     amat(1,a_b), cmat(1,c_b), nq,  ex, 
     &        qwght, ldew, func, fact_cs)
      endif
      if(ipol.eq.1) then
         call dscal(nq,2d0,rho,1)
         call dscal(nq*3,2d0,delrho,1)
      endif
      return
      end
      subroutine xc_os2cs_mgga(funcxc,
     &     tol_rho, fac,lfac,nlfac, rho, delrho, tau,
     &     amat, cmat, mmat,
     N     nq, ipol, ex, 
     &     qwght, ldew, func)
      implicit none
c     convert open-shell fxc to closed-shell
#include "errquit.fh"
#include "dft2drv.fh"
c      
      external funcx
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision delrho(nq,3,ipol)
      double precision tau(nq,ipol)
      double precision qwght(nq)
      double precision Amat(nq,ipol), Cmat(nq,*),mmat(nq,*)
      double precision tol_rho, pi,fact_cs
      integer r_a,r_b,dr_a,dr_b
      integer a_a,a_b,c_a,c_b
      if(ipol.eq.1) then
         r_a=1
         r_b=1
         dr_a=1
         dr_b=1
         a_a=D1_RA
         a_b=D1_RA
         c_a=D1_GAA
         c_b=D1_GAA
         fact_cs=2d0
         call dscal(nq,0.5d0,rho,1)
         call dscal(nq*3,0.5d0,delrho,1)
      else
         r_a=2
         r_b=3
         dr_a=1
         dr_b=2
         a_a=D1_RA
         a_b=D1_RB
         c_a=D1_GAA
         c_b=D1_GBB
         fact_cs=1d0
         call dscal(nq*ipol,2d0,tau,1)
      endif
c     alpha
      call funcxc(tol_rho, fac,lfac,nlfac,
     C     rho(1,r_a), delrho(1,1,dr_a),
     t     tau(1,dr_a),
     &     amat(1,a_a), cmat(1,c_a), mmat(1,dr_a),
     N     nq,  ex, 
     &     qwght, ldew, func, fact_cs)
c     beta
      if (ipol.eq.2) then
      call funcxc(tol_rho, fac,lfac,nlfac,
     C        rho(1,r_b), delrho(1,1,dr_b),
     t        tau(1,dr_b),
     &        amat(1,a_b), cmat(1,c_b), mmat(1,dr_b),
     n        nq,  ex, 
     &        qwght, ldew, func, fact_cs)
      endif
      if(ipol.eq.1) then
         call dscal(nq,2d0,rho,1)
         call dscal(nq*3,2d0,delrho,1)
      else
         call dscal(nq*ipol,0.5d0,tau,1)
      endif
      return
      end
