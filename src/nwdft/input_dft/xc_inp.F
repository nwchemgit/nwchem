      subroutine xc_input(rtdb, module)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "xc.fh"
c
      integer rtdb
      character*(*) module
c     
c     Parse the xc directive which specifies how to construct the
c     exchange-correlation potential.  Insert options into the 
c     database prefixed with module:...
c     
c     Possible variables are:
c
c     acm
c     b3lyp
c     becke88
c     beckehandh
c     HFexch
c     lyp
c     perdew81
c     perdew86
c     perdew91
c     pw91lda
c     slater
c     vwn_1
c     vwn_2
c     vwn_3
c     vwn_4
c     vwn_5
c     vwn_1_rpa
c     hcth
c     becke97
c     becke97-1
c     cpbe96 => 2100
c     cvpbe09 => 2109
c     xpbe96
c     xmpw91
c     ft97
c     xpkzb99 => 3600
c     cpkzb99 => 3650
c     xtpss03 => 4100
c     xmvs15  => 4115
c     ctpss03 => 4150
c     xctpssh => 4175
c     becke97-3 => 4500
c     ssb-d => 4600 (numfunc=46)
c     kt1 => 4700 (numfunc=47)
c     becke97-d => 4800
c     bc95 => 5000
c     mpw1b95=> 5100
c     mpwb1k => 5200
c     b1b95  => 5300
c     bb1k   => 5400
c     pw6b95 => 5500
c     pwb6k  => 5600
c     m05    => 5700
c     m05-2x => 5800
c     xpw6b95 => 5900
c     cpw6b95 => 6000
c     xpwb6k  => 6100
c     cpwb6k  => 6200
c     xm05    => 6300
c     cm05    => 6400
c     xm05-2x => 6500
c     cm05-2x => 6600
c     vs98    => 6700
c     xvs98   => 6800
c     cvs98   => 6900
c     m06-L    => 7000
c     xm06-L   => 7100
c     cm06-L   => 7200
c     m06-hf   => 7300
c     xm06-hf  => 7400
c     cm06-hf  => 7500
c     m06      => 7600
c     xm06     => 7700
c     cm06     => 7800
c     revm06   => 7610
c     revm06-L => 7010
c     m06-2x   => 7900
c     xm06-2x  => 8000
c     cm06-2x  => 8100
c     xbnl07   => 8200
c     xcamb88  => 8300      
c     xcampbe96 => 8400
c     mp2       => 8500
c     b2plyp    => 8510
c     xcamlsd   => 8520
c     xwpbe     => 8530
c     dldf      => 9300
c
c     m08-HX   => 8600
c     xm08-HX  => 8610
c     cm08-HX  => 8620
c
c     m08-SO   => 8700
c     xm08-SO  => 8710
c     cm08-SO  => 8720
c
c     m11      => 8800
c     xm11     => 8810
c     cm11     => 8820
c
c     revm11   => 8850
c
c     m11-L    => 8900
c     xm11-L   => 8910
c     cm11-L   => 8920
c      
c     mn12-L   => 8950
c     mn12-sx   => 8960
c     m06-sx   => 8962
c      
c     mn15-L   => 8970
c     mn15     => 8980
c     
c
c     sogga    => 9000
c     xsogga   => 9010
c     csogga   => 9020
c
c     sogga11  => 9100
c     xsogga11 => 9110
c     csogga11 => 9120
c
c     sogga11-x  => 9200
c     xsogga11-x => 9210
c     csogga11-x => 9220
c
c
c     n12 numfunc=45 > 9245
c     
c     n12-sx -> 9260
c      
c     srexch ==> 9310
c     hse03 ==> 9320
c
c     wb97 ==> 9410
c     wb97x ==> 9420
c     wb97x-d3 ==> 9430
c
c     becke86b ==> 9510
c     xperde86 ==> 9520
c
c     hcth147@tz2p ==> 1947
c
c     r2scan0 ==> 19
c
      integer num_dirs, ind, mlen, iline, n, itemp
      parameter (num_dirs=168)
c
      character*15 dirs(num_dirs)
      character*255 test
c
c     Define logical flags for local (l) and nonlocal (nl)
c     parts of exchange-correlation functionals.
c
      logical  lcfac(numfunc),  lxfac(numfunc), ck_local
      logical nlcfac(numfunc), nlxfac(numfunc)
      logical xccomb(numfunc)
      logical success
      logical,external :: nwchem_libxc_init
      logical,external :: is_libxcon
c
c     XC functionals coefficients
c     
      double precision mp2fac
      double precision cfac(numfunc),xfac(numfunc)
      data dirs /'acm', 'b3lyp', 'beckehandh', 'HFexch', 'becke88',  
     &           'lyp', 'perdew81', 'perdew86', 'perdew91', 'pw91lda', 
     &           'slater', 'vwn_1', 'vwn_2', 'vwn_3', 'vwn_4', 
     &           'vwn_5','vwn_1_rpa','hcth','becke97','becke97-1',
     &           'cpbe96','xpbe96','gill96','xperdew91',
     ,     'xhcth','xbecke97','xbecke97-1',
     ,     'chcth','cbecke97','cbecke97-1','hcth120','hcth147',
     ,     'becke98','pbe0','becke97gga1','hcth407','hcthp14',
     ,     'optx','optc','becke97-2','mpw91',
     ,     'ft97','xft97','cft97',
     &     'xpkzb99','mpw1k','cpkzb99','hcth407p','xtpss03',
     &     'ctpss03', 'xctpssh' , 'becke97-3',
     E     'op', 'bop', 'pbeop','rpbe','revpbe',
     &     'bc95','mpw1b95',
     &     'mpwb1k','b1b95','bb1k','pw6b95','pwb6k','m05','m05-2x',
     &     'xpw6b95','cpw6b95','xpwb6k','cpwb6k','xm05','cm05',
     &     'xm05-2x','cm05-2x','vs98','xvs98','cvs98','m06-L',
     &     'xm06-L','cm06-L','m06-hf','xm06-hf','cm06-hf','m06',
     &     'revm06','revm06-L',
     &     'xm06','cm06','m06-2x','xm06-2x','cm06-2x',
     &     'm08-hx','xm08-hx','cm08-hx',
     &     'm08-so','xm08-so','cm08-so',
     &     'm11','xm11','cm11','revm11',
     &     'm11-l','xm11-l','cm11-l','mn12-l','mn12-sx','m06-sx',
     &     'mn15-l','mn15',
     &     'sogga','xsogga','csogga',
     &     'sogga11','xsogga11','csogga11',
     &     'sogga11-x','xsogga11-x','csogga11-x',
     %     'n12','xn12','cn12','n12-sx',
     &     'xbnl07',
     &     'becke97-d','xcamb88','xcampbe96','mp2','b2plyp',
     &     'xcamlsd','xwpbe','ssb-d','kt1','kt2','bhlyp','b3p86',
     &     'b3pw91','pbe96','dldf','srhfexch','hse03',
     &     's12g','s12h','cam-s12g','cam-s12h',
     &     'becke86b','xperdew86','xmvs15','cvpbe09','mvs',
     &     'hle16', 'scan', 'scanl','scan0','lb94', 'rscan',
     &     'xncap', 'ncap', 'r2scan', 'r2scanl',
     &     'xscan', 'cscan',
     &     'hcth147@tz2p',
     W     'wb97','wb97x','wb97x-d3', 
     &     'xr2scan', 'cr2scan','r2scan0' /
c     becke97-d, ssb-d
      logical disp, dumdirect
      integer ivdw
      integer ip
c
      logical cam_exch,cam_srhf
      double precision cam_omega,cam_alpha,cam_beta
c
      logical use_hcth147
c
c     Grab the current input line and dump the functional specification
c     onto the RTDB. This way external programs can pick the actual
c     functional definition up as it was entered.
c
      success = nwchem_libxc_init()
c
c     reset cam parameters
c
      call xc_setcamparam(rtdb,.false.,.false.,0d0,0d0,0d0)

      if (.not.inp_line(test)) then
         write(*,*)'xc_input: Could not read functional specification'
      endif
      if (.not.inp_contains(.false.,"xc",test,ip)) then
         call errquit("xc_input: Not a functional specification",0,
     +                INPUT_ERR)
      else
         ip = ip + 2
         do while (test(ip:ip).eq." ")
            ip = ip + 1
         enddo
         if (.not.rtdb_cput(rtdb,"dft:xc_spec",1,
     +            test(ip:inp_strlen(test)))) then
            call errquit("xc_input: rtdb_cput failed",0,RTDB_ERR)
         endif
      endif
c
c     Only these few arrays are initialized (not defaults).
c
      do iline = 1, numfunc
          xfac(iline)=0.d0
          cfac(iline)=0.d0
         lcfac(iline)=.false.
        nlcfac(iline)=.false.
         lxfac(iline)=.false.
        nlxfac(iline)=.false.
        xccomb(iline)=.false.
      enddo
c     
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 1999
c     
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
        call nwchem_libxc_read(rtdb,test,success,xfac,lxfac,
     $                         nlxfac)
        if (success) goto 10
c     
c        Does not match a keyword ... 
c     
         goto 10000
      endif
c
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 
     &      1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     ,      2100,2200,2300, 2400, 
     ,     2500, 2600, 2700, 2510, 2610, 2710,1812,1847,2800,2900,
     , 3000, 3100,3200,3300,3400,3500,3700,3800,3810,3820,
     &  3600,4000, 3650, 3110,4100,4150,4175,4500,
     E     1997, 2046, 2047, 2201, 2202,
     &  5000,5100,5200,5300,5400,5500,5600,5700,5800,
     &  5900,6000,6100,6200,6300,6400,6500,6600,6700,6800,6900,
     &     7000,7100,7200,7300,7400,7500,7600,7610,7010,
     &     7700,7800,7900,8000,8100,
     &  8600,8610,8620,8700,8710,8720,8800,8810,8820,8850,
     M  8900,8910,8920,8950,8960,8962,8970,8980,
     &  9000,9010,9020,9100,9110,9120,9200,9210,9220,
     N  9245,9246,9247,9260,
     &  8200,4800,8300,8400,8500,8510,8520,8530,4600,4700,4701,8531,
     &  8532,8533,8534,9300,9310,9320,4610,4620,4630,4640,
     &  9510,9520,4115,2109,4124,3101,
     &  9610,9620,9630,9494,9640,
     &     9710,9720,9650,9660,9615,9616,
     &  1947,
     &  9410,9420,9430,
     &  9661,9662,9670, 
     &     1999 ) ind


      call errquit('xc_inp: unimplemented directive', ind, INPUT_ERR)
c     
c     acm; adiabatic connection method
c          AD Becke JCP 98, 5648 (1993)
c          Exc(S+VWN) + A*(Ex(HF)-Ex(S)) + B*Ex(B88(nl)) + 
c          C*Ec(P91(nl)) 
c          A=0.20, B=0.72, C=0.81               
c             
  100 xccomb(1) = .true.
        cfac(1) = 1.00d0
       lcfac(1) = .true.
      nlcfac(1) = .false.
        cfac(5) = 0.81d0
       lcfac(5) = .false.
      nlcfac(5) = .true.
        xfac(1) = 0.20d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(2) = 0.80d0
       lxfac(2) = .true.
      nlxfac(2) = .false.
        xfac(3) = 0.72d0
       lxfac(3) = .false.
      nlxfac(3) = .true.
c
      goto 10
c     
c     b3lyp; Becke3LYP
c            AD Becke JCP 98, 5648 (1993)
c            A*Ex(HF) + (1-A)*Ex(S) + B*Ex(B88(nl)) + 
c            C*Ec(LYP(total)) + (1-C)*Ec(VWN_1_RPA)
c            
c            A=0.20, B=0.72, C=0.81               
c             
  200 xccomb(2) = .true.
        cfac(2) = 0.81d0
       lcfac(2) = .true.
      nlcfac(2) = .true.
        cfac(7) = 0.19d0
       lcfac(7) = .true.
      nlcfac(7) = .false.
        xfac(1) = 0.20d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(2) = 0.80d0
       lxfac(2) = .true.
      nlxfac(2) = .false.
        xfac(3) = 0.72d0
       lxfac(3) = .false.
      nlxfac(3) = .true.
c
      goto 10
c     
c     beckehandh; Becke Half and Half
c                 AD Becke JCP 98, 1372
c                 0.5*Exc(LSD) + 0.5*Ex(HF) + 0.5*Ec(PW91LDA)
c             
  300 xccomb(3) = .true.
        xfac(1) = 0.5d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(2) = 0.5d0
       lxfac(2) = .true.
      nlxfac(2) = .false.
        cfac(6) = 0.5d0
       lcfac(6) = .true.
      nlcfac(6) = .false.
c
      goto 10
c     
c     bhlyp; 0.5*Ex(HF)+0.5*Ex(LSDA)+0.5*Ex(B88(nl)+Ec(LYP(total)) 
c             
 8531 xccomb(55) = .true.
        xfac(1)  = 0.5d0
       lxfac(1)  = .true.
      nlxfac(1)  = .true.
        xfac(2)  = 0.5d0
       lxfac(2)  = .true.
      nlxfac(2)  = .false.
        xfac(3)  = 0.5d0
       lxfac(3)  = .false.
      nlxfac(3)  = .true.
        cfac(2)  = 1.0d0
       lcfac(2)  = .true.
      nlcfac(2)  = .true.
c
      goto 10
c
c     b3p86; A*Ex(HF) + (1-A)*Ex(S) + B*Ex(B88(nl)) + 
c            C*Ec(Perdew86(total)) + Ec(VWN_1_RPA)
c            
c            A=0.20, B=0.72, C=0.81               
c             
 8532 xccomb(56) = .true.
        cfac(4)  = 0.81d0
       lcfac(4)  = .false.
      nlcfac(4)  = .true.
        cfac(7)  = 1d0
       lcfac(7)  = .true.
      nlcfac(7)  = .false.
        xfac(1)  = 0.20d0
       lxfac(1)  = .true.
      nlxfac(1)  = .true.
        xfac(2)  = 0.80d0
       lxfac(2)  = .true.
      nlxfac(2)  = .false.
        xfac(3)  = 0.72d0
       lxfac(3)  = .false.
      nlxfac(3)  = .true.
c
      goto 10
c
c     b3pw91; A*Ex(HF) + (1-A)*Ex(S) + B*Ex(B88(nl)) + 
c            C*Ec(PW91(nl)) + 1*Ec(pw91lda)
c            
c            A=0.20, B=0.72, C=0.81               
c             
 8533 xccomb(57) = .true.
        cfac(5)  = 0.81d0
       lcfac(5)  = .false.
      nlcfac(5)  = .true.
        cfac(6)  = 1d0
       lcfac(6)  = .true.
      nlcfac(6)  = .false.
        xfac(1)  = 0.20d0
       lxfac(1)  = .true.
      nlxfac(1)  = .true.
        xfac(2)  = 0.80d0
       lxfac(2)  = .true.
      nlxfac(2)  = .false.
        xfac(3)  = 0.72d0
       lxfac(3)  = .false.
      nlxfac(3)  = .true.
c
      goto 10
c
c     PBE96; xpbe96 + cpbe96
c             
 8534 xccomb(58) = .true.
        xfac(7)  =  1.0d0
       lxfac(7)  = .true.
      nlxfac(7)  = .true.
        cfac(12) =  1.0d0
       lcfac(12) = .false.
      nlcfac(12) = .true.
       lcfac(6) = .true.
       cfac(6) = 1.0d0
c
      goto 10
c     
c     HFexch; HF exact exchange
c             
  400 if (inp_a(test) ) then
         if (inp_compare(.false.,'local', test) ) then
            write(LuOut,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total HF exchange is assumed desired. '
            lxfac(1) = .true.
            nlxfac(1) = .true.
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total HF exchange is assumed desired. '
            lxfac(1) = .true.
            nlxfac(1) = .true.
         else
c
c           assume total wanted
c
            lxfac(1) = .true.
            nlxfac(1) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(1) = .true.
         nlxfac(1) = .true.
      endif
      if (.not. inp_f(xfac(1)) ) then
         write(LuOut,*)' xc_inp: hfexch multiplicative factor not ',
     &                 'found.'
         xfac(1) = 1.0d0
      endif            
c
      goto 10
c     
c     becke88; Becke '88 functional for exchange
c     
  500 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' becke88 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
            lxfac(3) = .true.
            nlxfac(3) = .false.
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(3) = .true.
         else
c
c           assume total wanted
c
            lxfac(3) = .true.
            nlxfac(3) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(3) = .true.
         nlxfac(3) = .true.
      endif
      if (.not. inp_f(xfac(3))) then
         xfac(3) = 1.0d0
      endif
       goto 10
c     
c     Handys OPTX 2001
c     
 3300  continue
       if (.not. inp_f(xfac(16))) then
          xfac(16) = 1.43169d0
       endif
       nlxfac(16)=.true.
       lxfac(2)=.true.
       if(abs(xfac(2)).lt.1d-6) xfac(2)=1.05151d0
       goto 10
c     
c     Handys OPTC 2001
c     
 3400  continue
       if (.not. inp_f(cfac(22))) then
          cfac(22) = 1.d0
       endif
         lcfac(22) = .true.
         nlcfac(22) = .true.
       goto 10
c     
c     HCTH for exchange & correlation
c     
 1800  continue
c
c        assume total wanted
c     
       xccomb(6) = .true.
         lxfac(4) = .true.
         nlxfac(4) = .true.
         xfac(4) = 1.0d0
         lcfac(13) = .true.
         nlcfac(13) = .true.
         cfac(13) = 1.0d0
       goto 10
c     
c     HCTH120 for exchange & correlation
c     
 1812    continue
c
c        assume total wanted
c     
       xccomb(7) = .true.
         lxfac(10) = .true.
         nlxfac(10) = .true.
         xfac(10) = 1.0d0
         lcfac(16) = .true.
         nlcfac(16) = .true.
         cfac(16) = 1.0d0
       goto 10
c     
c     HCTH147 for exchange & correlation
c     
 1847    continue
c
c        assume total wanted
c
         if (.not. rtdb_get(rtdb, 'dft:use_hcth147', mt_log, 1,
     c        use_hcth147)) use_hcth147=.false.
         if(.not.use_hcth147) then
            write(luout,*) ' ********** WARNING *********'
            write(luout,*) ' Use of the hcth147 is deprecated,'
            write(luout,*) ' the hcth147@tz2p should be used instead.'
            write(luout,*) ' Details at '
        write(luout,*) ' https://github.com/nwchemgit/nwchem/pull/301'
            write(luout,*) ' ****************************'
            write(luout,*) ' If you still wish to use hcth147, '
            write(luout,*) ' add the folllowing input line '
            write(luout,*) ' before the dft input section'
            write(luout,*) '    set dft:use_hcth147 t'
            write(luout,*) ' ****************************'
            call errquit(' use of hcth147 deprecated ',0, INPUT_ERR)
         endif
       xccomb(8) = .true.
         lxfac(11) = .true.
         nlxfac(11) = .true.
         xfac(11) = 1.0d0
         lcfac(17) = .true.
         nlcfac(17) = .true.
         cfac(17) = 1.0d0
       goto 10
c
c
c     HCTH147@TZ2P for exchange & correlation
c
 1947    continue
c
c        assume total wanted
c
       xccomb(75) = .true.
         lxfac(75) = .true.
         nlxfac(75) = .true.
         xfac(75) = 1.0d0
         lcfac(75) = .true.
         nlcfac(75) = .true.
         cfac(75) = 1.0d0
       goto 10
c
c     becke 97 gga-1 by  Cohen and Handy
c     
 3000  continue
c
c        assume total wanted
c     
       xccomb(11) = .true.
         lxfac(13) = .true.
         nlxfac(13) = .true.
         xfac(13) = 1.0d0
         lcfac(19) = .true.
         nlcfac(19) = .true.
         cfac(19) = 1.0d0
       goto 10
c     
c     HCTH407 by Boese
c     
 3100  continue
c
c        assume total wanted
c     
       xccomb(12) = .true.
         lxfac(14) = .true.
         nlxfac(14) = .true.
         xfac(14) = 1.0d0
         lcfac(20) = .true.
         nlcfac(20) = .true.
         cfac(20) = 1.0d0
       goto 10
c     
c     HLE16 aka modified HCTH407 by Truhlar and Verma
c     
 3101  continue
c
c        assume total wanted
c     
       xccomb(65) = .true.
         lxfac(14) = .true.
         nlxfac(14) = .true.
         xfac(14) = 1.25d0
         lcfac(20) = .true.
         nlcfac(20) = .true.
         cfac(20) = 0.5d0
       goto 10
c     
c     HCTH407p by Boese
c     
 3110  continue
c
c        assume total wanted
c     
       xccomb(17) = .true.
         lxfac(20) = .true.
         nlxfac(20) = .true.
         xfac(20)= 1.0d0
         lcfac(26) = .true.
         nlcfac(26) = .true.
         cfac(26) = 1.0d0
       goto 10
c     
c     HCTHP14 by Tozer
c     
 3200  continue
c
c        assume total wanted
c     
       xccomb(13) = .true.
         lxfac(15) = .true.
         nlxfac(15) = .true.
         xfac(15) = 1.0d0
         lcfac(21) = .true.
         nlcfac(21) = .true.
         cfac(21) = 1.0d0
       goto 10
c     
c     becke1998 for exchange & correlation
c     
 2800  continue
c
c        assume total wanted
c     
       xccomb(9) = .true.
       xfac(1)=0.219847d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
         lxfac(12) = .true.
         nlxfac(12) = .true.
         xfac(12) = 1.0d0
         lcfac(18) = .true.
         nlcfac(18) = .true.
         cfac(18) = 1.0d0
       goto 10
c     
c     PBE0 for exchange & correlation
C     Adamo & Barone JCP 110, 6158 (1998)
c     
 2900  continue
c
c        assume total wanted
c     
       xccomb(10) = .true.
       xfac(1)=0.25d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
         lxfac(7) = .true.
         nlxfac(7) = .true.
         xfac(7) = .75d0
         lcfac(12) = .false.
         nlcfac(12) = .true.
         cfac(12) = 1.0d0
         lcfac(6) = .true.
         cfac(6) = 1.0d0
       goto 10
c     
c     mpw1k for exchange & correlation
C     Lynch, Fast, Harris, and Truhlar, J. Phys. Chem. A 104, 4811 (2000)
c     
 4000  continue
c
c        assume total wanted
c     
       xccomb(16) = .true.
       xfac(1)=0.428d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
         lxfac(19) = .true.
         nlxfac(19) = .true.
         xfac(19) = .572d0
         lcfac(5) = .false.
         nlcfac(5) = .true.
         cfac(5) = 1.0d0
         lcfac(6) = .true.
         cfac(6) = 1.0d0
       goto 10
c     
c     HCTH for exchange 
c     
 2500    continue
c

c     
         lxfac(4) = .true.
         nlxfac(4) = .true.
         xfac(4) = 1.0d0

c
      goto 10
c     
c     HCTH for correlation
c     
 2510 continue
c

c     
         lcfac(13) = .true.
         nlcfac(13) = .true.
         cfac(13) = 1.0d0

c
      goto 10
 2600 continue
c
c     Becke97 X
c
        xfac(1) = 0.1943d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(5) = 1.d0
       lxfac(5) = .true.
      nlxfac(5) = .true.
      goto 10
 2610 continue
c
c     Becke97 C
c
        cfac(14) = 1.d0
       lcfac(14) = .true.
      nlcfac(14) = .true.
      goto 10
 2700 continue
c
c     Becke97-1 X
c
        xfac(1) = 0.2100d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(6) = 1.d0
       lxfac(6) = .true.
      nlxfac(6) = .true.
      goto 10
 2710 continue
c
c     Becke97-1
c
        cfac(15) = 1.d0
       lcfac(15) = .true.
      nlcfac(15) = .true.
      goto 10
c     
c     ft97 Filatov Thiel
c     
 3800  continue
c
c        assume total wanted
c     
       xccomb(15) = .true.

         nlxfac(25) = .true.
         xfac(25) = 1.0d0
c     Slater XC
         lxfac(2) = .true.
         nlxfac(2) = .false.
         xfac(2)=1d0
         lcfac(24)=.true.
         nlcfac(24)=.true.
         cfac(24)=1d0
       goto 10
c     
c    X ft97 Filatov Thiel
c     
 3810  continue
c
c        assume total wanted
c     
         nlxfac(25) = .true.
         xfac(25) = 1.0d0
c     Slater XC
         xfac(2) = 1d0
         lxfac(2) = .true.
         nlxfac(2) = .false.
       goto 10
c     
c   C  ft97 Filatov Thiel
c     
 3820  continue
c
c        assume total wanted
c     
         lcfac(24)=.true.
         nlcfac(24)=.true.
         cfac(24)=1d0
       goto 10
c     
c     lyp; LYP functional for correlation
c             
  600 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' The lyp functional cannot be conveniently',
     &                ' split into local and nonlocal pieces; '
            write(LuOut,*)' total lyp is assumed desired. ' 
            lcfac(2) = .true.
            nlcfac(2) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' The lyp functional cannot be conveniently',
     &                ' split into local and nonlocal pieces; '
            write(LuOut,*)' total lyp is assumed desired. ' 
            lcfac(2) = .true.
            nlcfac(2) = .true.
         else
c
c           assume total wanted
c
            lcfac(2) = .true.
            nlcfac(2) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(2) = .true.
         nlcfac(2) = .true.
      endif
      if (.not. inp_f(cfac(2))) then
         cfac(2) = 1.0d0
      endif            
c
      goto 10
c     
c     perdew81; P81 functional for correlation
c             
  700 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(3) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' perdew81 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(3) = .false.
         else
c
c           assume total wanted
c
            lcfac(3) = .true.
            nlcfac(3) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(3) = .true.
         nlcfac(3) = .false.
      endif
      if (.not. inp_f(cfac(3))) then
         cfac(3) = 1.0d0
      endif            
c
      goto 10
c     
c     perdew86; P86 functional for correlation
c     default local functional (perdew81 LDA) is set later in this file
c             
  800 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' perdew86 is a nonlocal functional;',
     &                ' "local" descriptor is ignored. '                   
            lcfac(4) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(4) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(4) = .false.
            nlcfac(4) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(4) = .false.
         nlcfac(4) = .true.
      endif
      if (.not. inp_f(cfac(4))) then
         cfac(4) = 1.0d0
      endif            
c
      goto 10
c     
c     perdew91; P91 functional for correlation
c     default local functional (PW91 LDA) is set later in this file
c             
  900 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' perdew91 is a nonlocal functional;',
     &                ' "local" descriptor is ignored. '                   
c
c           This is also done by default inside xc_perdew91.F
c           when lcfac(5) is .true.
c
            lcfac(5) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(5) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(5) = .false.
            nlcfac(5) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(5) = .false.
         nlcfac(5) = .true.
      endif
      if (.not. inp_f(cfac(5))) then
         cfac(5) = 1.0d0
      endif            
c
      goto 10
c     
c     pw91lda; P91LDA functional for correlation
c     
 1000 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(6) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' pw91lda is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(6) = .false.
         else
c
c           assume total wanted
c
            lcfac(6) = .true.
            nlcfac(6) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(6) = .true.
         nlcfac(6) = .false.
      endif
      if (.not. inp_f(cfac(6))) then
         cfac(6) = 1.0d0
      endif            
c
      goto 10
c     
c     slater; Slater functional for exchange
c     
 1100 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lxfac(2) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlxfac(2) = .false.
            write(LuOut,*)' slater is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
         else
c
c           assume total wanted
c
            lxfac(2) = .true.
            nlxfac(2) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(2) = .true.
         nlxfac(2) = .false.
      endif
      if (.not. inp_f(xfac(2))) then
         xfac(2) = 1.0d0
      endif            
c
      goto 10
c     
c     vwn_1; VWN_1 functional for correlation
c     
 1200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(8) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_1 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(8) = .false.
         else
c
c           assume total wanted
c
            lcfac(8) = .true.
            nlcfac(8) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(8) = .true.
         nlcfac(8) = .false.
      endif
      if (.not. inp_f(cfac(8))) then
         cfac(8) = 1.0d0
      endif            
c
      goto 10
c     
c     vwn_2; VWN_2 functional for correlation
c     
 1300 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(9) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_2 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(9) = .false.
         else
c
c           assume total wanted
c
            lcfac(9) = .true.
            nlcfac(9) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(9) = .true.
         nlcfac(9) = .false.
      endif
      if (.not. inp_f(cfac(9))) then
         cfac(9) = 1.0d0
      endif            
c
      goto 10
c     
c     vwn_3; VWN_3 functional for correlation
c     
 1400 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(10) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_3 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(10) = .false.
         else
c
c           assume total wanted
c
            lcfac(10) = .true.
            nlcfac(10) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(10) = .true.
         nlcfac(10) = .false.
      endif
      if (.not. inp_f(cfac(10))) then
         cfac(10) = 1.0d0
      endif    
c
      goto 10
c     
c     vwn_4; VWN_4 functional for correlation
c     
 1500 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(11) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_4 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(11) = .false.
         else
c
c           assume total wanted
c
            lcfac(11) = .true.
            nlcfac(11) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(11) = .true.
         nlcfac(11) = .false.
      endif
      if (.not. inp_f(cfac(11))) then
         cfac(11) = 1.0d0
      endif    
c
      goto 10
c     
c     vwn_5; VWN_5 functional for correlation
c     
 1600 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(1) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_5 is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '                   
            nlcfac(1) = .false.
         else
c
c           assume total wanted
c
            lcfac(1) = .true.
            nlcfac(1) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(1) = .true.
         nlcfac(1) = .false.
      endif
      if (.not. inp_f(cfac(1))) then
         cfac(1) = 1.0d0
      endif
c
      goto 10
c     
c     vwn_1_rpa; VWN_1_RPA functional for correlation
c     
 1700 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lcfac(7) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' vwn_1_rpa is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(7) = .false.
         else
c
c           assume total wanted
c
            lcfac(7) = .true.
            nlcfac(7) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(7) = .true.
         nlcfac(7) = .false.
      endif
      if (.not. inp_f(cfac(7))) then
         cfac(7) = 1.0d0
      endif            
c
      goto 10
c     
c     
c     becke97; Becke 1997 Hybrid functional
c            AD Becke JCP 107, 8554 (1997)
c            c*Ex(HF) + Exc
c            
c            c=0.1943
c             
 1900 xccomb(4) = .true.
        xfac(1) = 0.1943d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(5) = 1.d0
       lxfac(5) = .true.
      nlxfac(5) = .true.
        cfac(14) = 1.d0
       lcfac(14) = .true.
      nlcfac(14) = .true.
      goto 10
c     
c     
c     becke97; Becke 1997 Hybrid functional
c            AD Becke JCP 107, 8554 (1997)
c     parameters re-fitted by
c     Hamprecht et al, JCP 109, 6264 (1998)
c            c*Ex(HF) + Exc
c            
c            c=0.2100
c             
 2000 xccomb(5) = .true.
        xfac(1) = 0.2100d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(6) = 1.d0
       lxfac(6) = .true.
      nlxfac(6) = .true.
        cfac(15) = 1.d0
       lcfac(15) = .true.
      nlcfac(15) = .true.
c
      goto 10
c     
c     becke97-2; Becke 1997 Hybrid functional
c            AD Becke JCP 107, 8554 (1997)
c     parameters re-fitted by
c     Wilson Bradley & Tozer JCP 115, 9233 (2001)
c            c*Ex(HF) + Exc
c            
c            c=0.2100
c             
 3500 xccomb(14) = .true.
        xfac(1) = 0.2100d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(17) = 1.d0
       lxfac(17) = .true.
      nlxfac(17) = .true.
        cfac(23) = 1.d0
       lcfac(23) = .true.
      nlcfac(23) = .true.
c
      goto 10
c     
c     becke97-3; Becke 1997 Hybrid functional
c            AD Becke JCP 107, 8554 (1997)
c     parameters re-fitted by
c     Keal & Tozer JCP 123, 121103 (2005)
c            c*Ex(HF) + Exc
c            
c            c=2.692880D-01
c             
 4500 xccomb(19) = .true.
        xfac(1) = 2.692880D-01
       lxfac(1) = .true.
      nlxfac(1) = .true.
        xfac(22) = 1.d0
       lxfac(22) = .true.
      nlxfac(22) = .true.
        cfac(28) = 1.d0
       lcfac(28) = .true.
      nlcfac(28) = .true.
c
      goto 10
c
c     becke97-d; Becke 1997 Hybrid functional
c            AD Becke JCP 107, 8554 (1997)
c     parameters re-fitted by
c     Grimme JCC 27, 1787 (2006)
c     to be used with empirical vdw contrib
c            c*Ex(HF) + Exc
c
c            c=0.0000
c
 4800 xccomb(37) = .true.
c       xfac(1) = 0.0d0
c      lxfac(1) = .false.
c     nlxfac(1) = .false.
        xfac(39) = 1.d0
       lxfac(39) = .true.
      nlxfac(39) = .true.
        cfac(41) = 1.d0
       lcfac(41) = .true.
      nlcfac(41) = .true.
c
      goto 10
c     
c     
c     X mpw91
c     by Adamo & Barone
c             
 3700  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' mpw91 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
            lxfac(19) = .true.
            nlxfac(19) = .false.
          
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(19) = .true.
         else
c
c           assume total wanted
c
            lxfac(19) = .true.
            nlxfac(19) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(19) = .true.
         nlxfac(19) = .true.
      endif
      if (.not. inp_f(xfac(19))) then
         xfac(19) = 1.0d0
      endif
      goto 10
c
c     
c     cpbe96 Perdew-Burke-Ernzerhof functional for correlation
c     default local functional (PW91 LDA) is set later in this file
c             
 2100 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' cpbe96 is a nonlocal functional;',
     &                ' "local" descriptor is ignored. '
c
c           This is also done by default inside xc_pbe96.F
c           when lcfac(12) is .true.
c
            lcfac(12) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(12) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(12) = .false.
            nlcfac(12) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(12) = .false.
         nlcfac(12) = .true.
      endif
      if (.not. inp_f(cfac(12))) then
         cfac(12) = 1.0d0
      endif            
c
      goto 10

c
c     cvpbe09 Perdew-Burke-Ernzerhof functional for correlation
c     but then the revTPSS variant version of it
c     default local functional (PW91 LDA) is set later in this file
c
 2109 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' cvpbe09 is a nonlocal functional;',
     &                ' "local" descriptor is ignored. '
c
c           This is also done by default inside xc_pbe96.F
c           when lcfac(64) is .true.
c
            lcfac(64) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(64) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(64) = .false.
            nlcfac(64) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(64) = .false.
         nlcfac(64) = .true.
      endif
      if (.not. inp_f(cfac(64))) then
         cfac(64) = 1.0d0
      endif
c
      goto 10

c     
c     xpbe96; Perdew-Burke-Ernzerhof functional for exchange
c     
 2200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xpbe96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(7) = .true.
         else
c
c           assume total wanted
c
            lxfac(7) = .true.
            nlxfac(7) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(7) = .true.
         nlxfac(7) = .true.
      endif
      if (.not. inp_f(xfac(7))) then
         xfac(7) = 1.0d0
      endif
       goto 10
c     
c     rpbe; Perdew-Burke-Ernzerhof + Norskov exchange 
c     
 2201  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' rpbe "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(30) = .true.
         else
c
c           assume total wanted
c
            lxfac(30) = .true.
            nlxfac(30) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(30) = .true.
         nlxfac(30) = .true.
      endif
      if (.not. inp_f(xfac(30))) then
         xfac(30) = 1.0d0
      endif
       goto 10
c     
c     revpbe; Perdew-Burke-Ernzerhof + W Yang exchange 
c     
 2202  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' rpbe "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(31) = .true.
         else
c
c           assume total wanted
c
            lxfac(31) = .true.
            nlxfac(31) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(31) = .true.
         nlxfac(31) = .true.
      endif
      if (.not. inp_f(xfac(31))) then
         xfac(31) = 1.0d0
      endif
       goto 10

c            META GGA

c     
c     xpkzb99; Perdew Kurth Blaha Zupan metaGGA functional for exchange
c     
 3600 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xpkzb "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(18) = .true.
         else
c
c           assume total wanted
c
            lxfac(18) = .true.
            nlxfac(18) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(18) = .true.
         nlxfac(18) = .true.
      endif
      if (.not. inp_f(xfac(18))) then
         xfac(18) = 1.0d0
      endif
       goto 10


c     cpkzb99 Perdew Kurth Zupan Blaha metaGGA functional for correlation
c     default local functional (PW91 LDA) is set in cpkzb99.F code
c             
 3650 continue
c
c           assume total wanted (but total here = nonlocal)
c
      lcfac(25) = .true.
      nlcfac(25) = .true.
      call xcinp_fact_noloc(test,'cpkzb99',cfac(25))
      if (.not. inp_f(cfac(25))) then
         cfac(25) = 1.0d0
      endif            
c
      goto 10


c     
c     xtpss03  metaGGA functional for EXCHANGE
c     
 4100 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xtpss "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using only slater exchange. '
            call errquit(' xc_inp - functional description ',0,0)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' xtpss "nonlocal" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without nonlocal descriptor. '
            call errquit(' xc_inp - functional description ',0,0)



c            nlxfac(21) = .true.


         else
c
c  assume total wanted but total is all that is allowed for metaGGAs
c
            lxfac(21) = .true.
            nlxfac(21) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(21) = .true.
         nlxfac(21) = .true.
      endif
      if (.not. inp_f(xfac(21))) then
         xfac(21) = 1.0d0
      endif
       goto 10


c
c     xmvs15  metaGGA functional for EXCHANGE
c
 4115 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xmvs15 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using only slater exchange. '
            call errquit(' xc_inp - functional description ',0,0)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' xmvs15 "nonlocal" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without nonlocal descriptor. '
            call errquit(' xc_inp - functional description ',0,0)



c            nlxfac(64) = .true.


         else
c
c  assume total wanted but total is all that is allowed for metaGGAs
c
            lxfac(64) = .true.
            nlxfac(64) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(64) = .true.
         nlxfac(64) = .true.
      endif
      if (.not. inp_f(xfac(64))) then
         xfac(64) = 1.0d0
      endif
       goto 10

c     MVS (Made Very Simple) metagga functional
c     Sun, Perdew, Ruzsinszky, PNAS 2015, 112, 685-689

 4124 xccomb(64) = .true.
       lcfac(64) = .false.
      nlcfac(64) = .true.
        cfac(64) = 1.0d0
       lcfac(6)  = .true.
      nlcfac(6)  = .false.
        cfac(6)  = 1.0d0
       lxfac(64) = .true.
      nlxfac(64) = .true.
        xfac(64) = 1.0d0

      goto 10

c     TPSS  metaGGA functional for CORRELATION
c     default local functional (PW91 LDA) 
c     and PBE96  is set in tpss03.F code
c     LOCAL part is always false
 4150 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' ctpss03 is a nonlocal functional;',
     &                ' "local" descriptor is ignored. '                   
           call inp_prev_field()
         elseif (inp_compare(.false.,'nonlocal', test)) then
           write(LuOut,*)' ctpss03 cannot be conveniently',
     &           ' split into local and nonlocal pieces; '
            write(LuOut,*)' total ctpss03 is used. '
            call inp_prev_field()
         endif
      endif
      lcfac(27) = .true.
      nlcfac(27) = .true.
      if (.not. inp_f(cfac(27))) then
         cfac(27) = 1.0d0
      endif            
c
      goto 10




c     
c TPSSH      hybrid meta GGA
c            A*EXX + (1-a)TPSSX + TPSSC
c            A=0.10               
c             
4175  xccomb(18) = .true.
         xfac(21)= 0.90d0
        lxfac(21)=.true.
       nlxfac(21)=.true.
        cfac(27) = 1.0d0
       lcfac(27) = .true.
      nlcfac(27) = .true.
        xfac(1) = 0.10d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c
c SCAN meta GGA
c
9610  xccomb(66) = .true.
       lcfac(66) = .true.
      nlcfac(66) = .true.
        cfac(66) = 1.0d0
       lxfac(66) = .true.
      nlxfac(66) = .true.
        xfac(66) = 1.0d0
      goto 10
c
c X-only SCAN meta GGA
c
 9615 lxfac(66) = .true.
      nlxfac(66) = .true.
      call xcinp_fact_noloc(test,'xscan',xfac(66))
      goto 10
c
c C-only SCAN meta GGA
c
 9616 lcfac(66) = .true.
      nlcfac(66) = .true.
      call xcinp_fact_noloc(test,'cscan',cfac(66))
      goto 10
c
c SCAN0 hybrid meta GGA
c
9630  xccomb(70) = .true.
       lcfac(66) = .true.
      nlcfac(66) = .true.
        cfac(66) = 1.0d0
       lxfac(66) = .true.
      nlxfac(66) = .true.
        xfac(66) = .75d0
       xfac(1)=0.25d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
      goto 10
c
c regularized SCAN meta GGA
c
9640  xccomb(71) = .true.
       lcfac(71) = .true.
      nlcfac(71) = .true.
        cfac(71) = 1.0d0
       lxfac(71) = .true.
      nlxfac(71) = .true.
        xfac(71) = 1.0d0
      goto 10
c
c SCAN-L laplacian-dependent mGGA
c
9620  xccomb(67) = .true.
       lcfac(67) = .true.
      nlcfac(67) = .true.
        cfac(67) = 1.0d0
       lxfac(67) = .true.
      nlxfac(67) = .true.
        xfac(67) = 1.0d0
      goto 10
c
c NCAP functional
c
9710   lxfac(72) = .true.
      nlxfac(72) = .true.
        xfac(72) = 1.0d0
      goto 10

9720  xccomb(72) = .true.
      cfac(3) = 1d0
      lcfac(3) = .true.
      cfac(4) = 1d0
      nlcfac(4) = .true.
      lxfac(72) = .true.
      nlxfac(72) = .true.
      xfac(72) = 1.0d0
      goto 10
c
c r2SCAN metaGGA
c
9650  xccomb(73) = .true.
       lcfac(73) = .true.
      nlcfac(73) = .true.
        cfac(73) = 1.0d0
       lxfac(73) = .true.
      nlxfac(73) = .true.
        xfac(73) = 1.0d0
      goto 10
c
c r2SCAN-L deorbitalized metaGGA
c
9660  xccomb(74) = .true.
       lcfac(74) = .true.
      nlcfac(74) = .true.
        cfac(74) = 1.0d0
       lxfac(74) = .true.
      nlxfac(74) = .true.
        xfac(74) = 1.0d0
      goto 10
c
c X-only r2SCAN meta GGA
c
 9661 lxfac(73) = .true.
      nlxfac(73) = .true.
      call xcinp_fact_noloc(test,'xr2scan',xfac(73))
      goto 10
c
c C-only r2SCAN meta GGA
c
 9662 lcfac(73) = .true.
      nlcfac(73) = .true.
      call xcinp_fact_noloc(test,'cr2scan',cfac(73))
      goto 10
c
c r2SCAN0 hybrid meta GGA
c
9670  xccomb(86) = .true.
       lcfac(73) = .true.
      nlcfac(73) = .true.
        cfac(73) = 1.0d0
       lxfac(73) = .true.
      nlxfac(73) = .true.
        xfac(73) = .75d0
       xfac(1)=0.25d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
      goto 10
c
c lb94 Vxc function
c fill only x here, no c
c      
 9494 continue
       lxfac(70) = .true.
      nlxfac(70) = .false.
        xfac(70) = 1.0d0
      goto 10

c     bc95


 5000 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' bc95 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' bc95 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(31) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'bc95 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(31) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(31) = .false.
            nlcfac(31) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(31) = .false.
         nlcfac(31) = .true.
      endif
            write(LuOut,*)'Bc95 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(31))) then
         cfac(31) = 1.0d0
      endif
c
      goto 10

c     mpw1b95
 5100  xccomb(21) = .true.
         xfac(19)= 0.69d0
        lxfac(19)=.true.
       nlxfac(19)=.true.
        cfac(31) = 1.0d0
       lcfac(31) = .true.
      nlcfac(31) = .true.
        xfac(1) = 0.31d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c     mpwb1k
 5200  xccomb(22) = .true.
         xfac(19)= 0.56d0
        lxfac(19)=.true.
       nlxfac(19)=.true.
        cfac(31) = 1.0d0
       lcfac(31) = .true.
      nlcfac(31) = .true.
        xfac(1) = 0.44d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c     
      goto 10
c     b1b95
 5300  xccomb(23) = .true.
         xfac(3)= 0.72d0
        lxfac(3)=.true.
       nlxfac(3)=.true.
        cfac(31) = 1.0d0
       lcfac(31) = .true.
      nlcfac(31) = .true.
        xfac(1) = 0.28d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c     bb1k
 5400  xccomb(24) = .true.
         xfac(3)= 0.58d0
        lxfac(3)=.true.
       nlxfac(3)=.true.
        cfac(31) = 1.0d0
       lcfac(31) = .true.
      nlcfac(31) = .true.
        xfac(1) = 0.42d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c     pw6b95
 5500  xccomb(25) = .true.
         xfac(26)= 0.72d0
        lxfac(26)=.true.
       nlxfac(26)=.true.
        cfac(32) = 1.0d0
       lcfac(32) = .true.
      nlcfac(32) = .true.
        xfac(1) = 0.28d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c     pwb6k
 5600  xccomb(26) = .true.
         xfac(27)= 0.54d0
        lxfac(27)=.true.
       nlxfac(27)=.true.
        cfac(33) = 1.0d0
       lcfac(33) = .true.
      nlcfac(33) = .true.
        xfac(1) = 0.46d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

 5700  xccomb(29) = .true.
         xfac(28)= 0.72d0
        lxfac(28)=.true.
       nlxfac(28)=.true.
        cfac(34) = 1.0d0
       lcfac(34) = .true.
      nlcfac(34) = .true.
        xfac(1) = 0.28d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

 5800  xccomb(30) = .true.
         xfac(29)= 0.44d0
        lxfac(29)=.true.
       nlxfac(29)=.true.
        cfac(35) = 1.0d0
       lcfac(35) = .true.
      nlcfac(35) = .true.
        xfac(1) = 0.56d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

 5900  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xpw6b95 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(26) = .true.
         else
c
c           assume total wanted
c
            lxfac(26) = .true.
            nlxfac(26) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(26) = .true.
         nlxfac(26) = .true.
      endif
      if (.not. inp_f(xfac(26))) then
         xfac(26) = 1.0d0
      endif
       goto 10
  
 6000 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cpw6b95 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cpw6b95 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)
            

c
c
c
            lcfac(32) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cpw6b95 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(32) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(32) = .false.
            nlcfac(32) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(32) = .false.
         nlcfac(32) = .true.
      endif
            write(LuOut,*)'cpw6b95 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(32))) then
         cfac(32) = 1.0d0
      endif
c
      goto 10


 6100  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xpwb6k "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(27) = .true.
         else
c
c           assume total wanted
c
            lxfac(27) = .true.
            nlxfac(27) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(27) = .true.
         nlxfac(27) = .true.
      endif
      if (.not. inp_f(xfac(27))) then
         xfac(27) = 1.0d0
      endif
       goto 10

 6200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cpwb6k is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cpwb6k "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(33) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cpwb6k is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(33) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(33) = .false.
            nlcfac(33) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(33) = .false.
         nlcfac(33) = .true.
      endif
            write(LuOut,*)'cpwb6k uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(33))) then
         cfac(33) = 1.0d0
      endif
c
      goto 10




 6300  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm05 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(28) = .true.
         else
c
c           assume total wanted
c
            lxfac(28) = .true.
            nlxfac(28) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(28) = .true.
         nlxfac(28) = .true.
      endif
      if (.not. inp_f(xfac(28))) then
         xfac(28) = 1.0d0
      endif
       goto 10

 6400 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm05 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm05 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(34) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm05 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(34) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(34) = .false.
            nlcfac(34) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(34) = .false.
         nlcfac(34) = .true.
      endif
            write(LuOut,*)'cm05 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(34))) then
         cfac(34) = 1.0d0
      endif
c
      goto 10


 6500  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm05-2x "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(29) = .true.
         else
c
c           assume total wanted
c
            lxfac(29) = .true.
            nlxfac(29) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(29) = .true.
         nlxfac(29) = .true.
      endif
      if (.not. inp_f(xfac(29))) then
         xfac(29) = 1.0d0
      endif
       goto 10

 6600 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm05-2x is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm05-2x "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(35) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm05-2x is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(35) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(35) = .false.
            nlcfac(35) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(35) = .false.
         nlcfac(35) = .true.
      endif
            write(LuOut,*)'cm05-2x uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(35))) then
         cfac(35) = 1.0d0
      endif
c
      goto 10


c
c VS98, also called VSXC
c
 6700  xccomb(32) = .true.
         xfac(33)= 1.0d0
        lxfac(33)=.true.
       nlxfac(33)=.true.
        cfac(29) = 1.0d0
       lcfac(29) = .true.
      nlcfac(29) = .true.
c
      goto 10
c
c VS98 exchange part
c
 6800  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xVS98 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(33) = .true.
         else
c
c           assume total wanted
c
            lxfac(33) = .true.
            nlxfac(33) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(33) = .true.
         nlxfac(33) = .true.
      endif
      if (.not. inp_f(xfac(33))) then
         xfac(33) = 1.0d0
      endif
       goto 10

 6900 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cvs98 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cvs98 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(29) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cvs98 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(29) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(29) = .false.
            nlcfac(29) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(29) = .false.
         nlcfac(29) = .true.
      endif
            write(LuOut,*)'cvs98 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(29))) then
         cfac(29) = 1.0d0
      endif
c
      goto 10

c
cM06 suite
c

c
c M06-L
c
 7000  xccomb(33) = .true.
         xfac(34)= 1.0d0
        lxfac(34)=.true.
       nlxfac(34)=.true.
        cfac(30) = 1.0d0
       lcfac(30) = .true.
      nlcfac(30) = .true.
c
      goto 10
c
c revM06-L
c
 7010  xccomb(69) = .true.
         xfac(69)= 1.0d0
        lxfac(69)=.true.
       nlxfac(69)=.true.
        cfac(69) = 1.0d0
       lcfac(69) = .true.
      nlcfac(69) = .true.
c
      goto 10
c
c M06-L exchange part
c
 7100  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm06-L "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(34) = .true.
         else
c
c           assume total wanted
c
            lxfac(34) = .true.
            nlxfac(34) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(34) = .true.
         nlxfac(34) = .true.
      endif
      if (.not. inp_f(xfac(34))) then
         xfac(34) = 1.0d0
      endif
       goto 10
c
c M06-L correlation part
c
 7200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm06-L is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm06-L "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(30) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm06-L is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(30) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(30) = .false.
            nlcfac(30) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(30) = .false.
         nlcfac(30) = .true.
      endif
            write(LuOut,*)'cm06-L uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(30))) then
         cfac(30) = 1.0d0
      endif
c
      goto 10
c
c  M06-HF
c
 7300  xccomb(34) = .true.
         xfac(35)= 1.0d0
        lxfac(35)=.true.
       nlxfac(35)=.true.
        cfac(38) = 1.0d0
       lcfac(38) = .true.
      nlcfac(38) = .true.
        xfac(1) = 1.0d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

c
c M06-HF exchange part
c
 7400  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm06-hf "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(35) = .true.
         else
c
c           assume total wanted
c
            lxfac(35) = .true.
            nlxfac(35) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(35) = .true.
         nlxfac(35) = .true.
      endif
      if (.not. inp_f(xfac(35))) then
         xfac(35) = 1.0d0
      endif
       goto 10
c
c M06-HF correlation part
c
 7500 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm06-hf is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm06-hf "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)

c
c
c
            lcfac(38) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm06-hf is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(38) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(38) = .false.
            nlcfac(38) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(38) = .false.
         nlcfac(38) = .true.
      endif
            write(LuOut,*)'cm06-hf uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(38))) then
         cfac(38) = 1.0d0
      endif
c
      goto 10

c
c M06
c
 7600  xccomb(35) = .true.
         xfac(36)= 1.0d0
        lxfac(36)=.true.
       nlxfac(36)=.true.
        cfac(39) = 1.0d0
       lcfac(39) = .true.
      nlcfac(39) = .true.
        xfac(1) = 0.27d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

c
c revM06
c
 7610  xccomb(68) = .true.
         xfac(68)= 1.0d0
        lxfac(68)=.true.
       nlxfac(68)=.true.
        cfac(68) = 1.0d0
       lcfac(68) = .true.
      nlcfac(68) = .true.
        xfac(1) = 0.4041d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10
c
c M06 exchange part
c
 7700  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm06 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(36) = .true.
         else
c
c           assume total wanted
c
            lxfac(36) = .true.
            nlxfac(36) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(36) = .true.
         nlxfac(36) = .true.
      endif
      if (.not. inp_f(xfac(36))) then
         xfac(36) = 1.0d0
      endif
       goto 10
c
c M06 correlation part
c
 7800 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm06 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm06 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)

c
c
c
            lcfac(39) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm06 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(39) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(39) = .false.
            nlcfac(39) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(39) = .false.
         nlcfac(39) = .true.
      endif
            write(LuOut,*)'cm06 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(39))) then
         cfac(39) = 1.0d0
      endif
c
      goto 10

c
c  M06-2X
c
 7900  xccomb(36) = .true.
         xfac(37)= 1.0d0
        lxfac(37)=.true.
       nlxfac(37)=.true.
        cfac(40) = 1.0d0
       lcfac(40) = .true.
      nlcfac(40) = .true.
        xfac(1) = 0.54d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

c
c M06-2X exchange part
c
 8000  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm06-2x "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(37) = .true.
         else
c
c           assume total wanted
c
            lxfac(37) = .true.
            nlxfac(37) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(37) = .true.
         nlxfac(37) = .true.
      endif
      if (.not. inp_f(xfac(37))) then
         xfac(37) = 1.0d0
      endif
cadded HF X bit
         xfac(1) = 0.54d0*xfac(37)
         lxfac(1) = .true.
         nlxfac(1) = .true.
       goto 10
c
c M06-2X correlation part
c
 8100 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm06-2x is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm06-2x "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)

c
c
c
            lcfac(40) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm06-2x is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(40) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(40) = .false.
            nlcfac(40) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(40) = .false.
         nlcfac(40) = .true.
      endif
            write(LuOut,*)'cm06-2x uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(40))) then
         cfac(40) = 1.0d0
      endif
c
      goto 10
c
cM08 suite
c
c M08-HX
c
 8600  xccomb(48) = .true.
         xfac(48)= 1.0d0
        lxfac(48)=.true.
       nlxfac(48)=.true.
        cfac(48) = 1.0d0
       lcfac(48) = .true.
      nlcfac(48) = .true.
        xfac(1) = 0.5223d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c    
      goto 10
c
c M08-HX exchange part
c
 8610  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm08-HX "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(48) = .true.
         else
c
c           assume total wanted
c
            lxfac(48) = .true.
            nlxfac(48) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(48) = .true.
         nlxfac(48) = .true.
      endif
      if (.not. inp_f(xfac(48))) then
         xfac(48) = 1.0d0
      endif
cadded HF X bit
        xfac(1) = 0.5223d0*xfac(48)
       lxfac(1) = .true.
      nlxfac(1) = .true.
       goto 10
c
c M08-HX correlation part
c
 8620 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm08-hx is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm08-hx "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(48) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm08-hx is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(48) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(48) = .false.
            nlcfac(48) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(48) = .false.
         nlcfac(48) = .true.
      endif
            write(LuOut,*)'cm08-hx uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(48))) then
         cfac(48) = 1.0d0
      endif
c
      goto 10
c
c  M08-SO
c
 8700  xccomb(49) = .true.
         xfac(49)= 1.0d0
        lxfac(49)=.true.
       nlxfac(49)=.true.
        cfac(49) = 1.0d0
       lcfac(49) = .true.
      nlcfac(49) = .true.
        xfac(1) = 0.5679d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
      goto 10

c
c M08-SO exchange part
c
 8710  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm08-so "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(49) = .true.
         else
c
c           assume total wanted
c
            lxfac(49) = .true.
            nlxfac(49) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(49) = .true.
         nlxfac(49) = .true.
      endif
      if (.not. inp_f(xfac(49))) then
         xfac(49) = 1.0d0
      endif
cadded HF X bit
        xfac(1) = 0.5679d0*xfac(49)
       lxfac(1) = .true.
      nlxfac(1) = .true.
       goto 10
c
c M08-SO correlation part
c
 8720 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm08-so is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm08-so "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)

c
c
c
            lcfac(49) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm08-so is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(49) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(49) = .false.
            nlcfac(49) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(49) = .false.
         nlcfac(49) = .true.
      endif
            write(LuOut,*)'cm08-so uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(49))) then
         cfac(49) = 1.0d0
      endif
c
      goto 10
      
c
cM11 suite
c
c M11
c
 8800  xccomb(50) = .true.
         xfac(50)= 1.0d0
        lxfac(50)=.true.
       nlxfac(50)=.true.
        cfac(50) = 1.0d0
       lcfac(50) = .true.
      nlcfac(50) = .true.
        xfac(1) = 1d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.25
      cam_alpha = 0.428
      cam_beta = 0.572
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &    cam_beta)
c
      goto 10
c
c M11 exchange part
c
 8810  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm11 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(50) = .true.
         else
c
c           assume total wanted
c
            lxfac(50) = .true.
            nlxfac(50) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(50) = .true.
         nlxfac(50) = .true.
      endif
      if (.not. inp_f(xfac(50))) then
         xfac(50) = 1.0d0
      endif
cadded HF X bit
        xfac(1) = 1d0*xfac(50)
       lxfac(1) = .true.
      nlxfac(1) = .true.
c    
      goto 10
c
c M11 correlation part
c
 8820 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm11 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm11 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(50) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm11 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(50) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(50) = .false.
            nlcfac(50) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(50) = .false.
         nlcfac(50) = .true.
      endif
            write(LuOut,*)'cm11 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(50))) then
         cfac(50) = 1.0d0
      endif
c
      goto 10
c
c revM11
c
 8850  xccomb(79) = .true.
         xfac(79)= 1.0d0
        lxfac(79)=.true.
       nlxfac(79)=.true.
        cfac(79) = 1.0d0
       lcfac(79) = .true.
      nlcfac(79) = .true.
        xfac(1) = 1d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.40
      cam_alpha = 0.225
      cam_beta = 0.775
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &    cam_beta)
c
      goto 10
c
c  M11-L
c
 8900  xccomb(51) = .true.
         xfac(51)= 1.0d0
        lxfac(51)=.true.
       nlxfac(51)=.true.
        cfac(51) = 1.0d0
       lcfac(51) = .true.
      nlcfac(51) = .true.
c
      goto 10

c
c M11-L exchange part
c
 8910  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xm11-l "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(51) = .true.
         else
c
c           assume total wanted
c
            lxfac(51) = .true.
            nlxfac(51) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(51) = .true.
         nlxfac(51) = .true.
      endif
      if (.not. inp_f(xfac(51))) then
         xfac(51) = 1.0d0
      endif
       goto 10
c
c M11-L correlation part
c
 8920 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm11-l is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' cm11-l "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)

c
c
c
            lcfac(51) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm11-l is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(51) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(51) = .false.
            nlcfac(51) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(51) = .false.
         nlcfac(51) = .true.
      endif
            write(LuOut,*)'cm11-l uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(51))) then
         cfac(51) = 1.0d0
      endif
c
      goto 10
c
c  MN12-L
c
 8950 xccomb(80) = .true.
         xfac(80)= 1.0d0
        lxfac(80)=.true.
       nlxfac(80)=.true.
        cfac(80) = 1.0d0
       lcfac(80) = .true.
      nlcfac(80) = .true.
c
      goto 10
c
c  MN12-SX
c
 8960 xccomb(81) = .true. 
c
c       Ex(HF-SR)
c
        xfac(1) = 0.25d0
       lxfac(1) = .true.
       nlxfac(1) = .true.
c       
        xfac(81)= 1.0d0
        lxfac(81)=.true.
       nlxfac(81)=.true.
        cfac(81) = 1.0d0
       lcfac(81) = .true.
       nlcfac(81) = .true.
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.11d0
      cam_alpha = 0.d0
      cam_beta =  1.0
      cam_exch = .true.
      cam_srhf = .true. ! for short-range HF
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &     cam_beta)
c
      goto 10
c
c  M06-SX
c
 8962 xccomb(81) = .true. 
c
c       Ex(HF-SR)
c
        xfac(1) = 0.335d0
       lxfac(1) = .true.
       nlxfac(1) = .true.
c       
        xfac(85)= 1.0d0
        lxfac(85)=.true.
       nlxfac(85)=.true.
        cfac(85) = 1.0d0
       lcfac(85) = .true.
       nlcfac(85) = .true.
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.10d0
      cam_alpha = 0.d0
      cam_beta =  1.0
      cam_exch = .true.
      cam_srhf = .true. ! for short-range HF
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &     cam_beta)
c
      goto 10
c
c  MN15-L
c
 8970 xccomb(83) = .true.
         xfac(83)= 1.0d0
        lxfac(83)=.true.
       nlxfac(83)=.true.
        cfac(83) = 1.0d0
       lcfac(83) = .true.
      nlcfac(83) = .true.
c
      goto 10
c
c  MN15
c
 8980 xccomb(84) = .true.
         xfac(84)= 1.0d0
        lxfac(84)=.true.
       nlxfac(84)=.true.
        xfac(1) = 0.44d0
       lxfac(1) = .true.
       nlxfac(1) = .true.
        cfac(84) = 1.0d0
       lcfac(84) = .true.
      nlcfac(84) = .true.
c
      goto 10
c
cSOGGA suite
c

c
c SOGGA
c
 9000  xccomb(52) = .true.
         xfac(52)= 1.0d0
        lxfac(52)=.true.
       nlxfac(52)=.true.
         lcfac(12) = .false.
         nlcfac(12) = .true.
         cfac(12) = 1.0d0
         lcfac(6) = .true.
         cfac(6) = 1.0d0
c    
      goto 10
c
c SOGGA exchange part
c
 9010  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' sogga "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(52) = .true.
         else
c
c           assume total wanted
c
            lxfac(52) = .true.
            nlxfac(52) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(52) = .true.
         nlxfac(52) = .true.
      endif
      if (.not. inp_f(xfac(52))) then
         xfac(52) = 1.0d0
      endif
       goto 10
c
c SOGGA correlation part
c
 9020 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm11 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' sogga "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(52) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'cm11 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(52) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(52) = .false.
            nlcfac(52) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(52) = .false.
         nlcfac(52) = .true.
      endif
            write(LuOut,*)'cm11 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(52))) then
         cfac(52) = 1.0d0
      endif
c
      goto 10      
c
c SOGGA11
c
 9100  xccomb(53) = .true.
         xfac(53)= 1.0d0
        lxfac(53)=.true.
       nlxfac(53)=.true.
        cfac(53) = 1.0d0
       lcfac(53) = .true.
      nlcfac(53) = .true.
c    
      goto 10
c
c SOGGA11 exchange part
c
 9110  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' sogga11 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(53) = .true.
         else
c
c           assume total wanted
c
            lxfac(53) = .true.
            nlxfac(53) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(53) = .true.
         nlxfac(53) = .true.
      endif
      if (.not. inp_f(xfac(53))) then
         xfac(53) = 1.0d0
      endif
       goto 10
c
c SOGGA11 correlation part
c
 9120 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm11 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' sogga11 "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(53) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'csogga11 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(53) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(53) = .false.
            nlcfac(53) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(53) = .false.
         nlcfac(53) = .true.
      endif
            write(LuOut,*)'csogga11 uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(53))) then
         cfac(53) = 1.0d0
      endif
c
      goto 10  
c
c SOGGA11-X
c

 9200  xccomb(54) = .true.
         xfac(54)= 1.0d0
        lxfac(54)=.true.
       nlxfac(54)=.true.
        cfac(54) = 1.0d0
       lcfac(54) = .true.
      nlcfac(54) = .true.
        xfac(1) = 0.4015d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
      goto 10
c
c SOGGA11-X exchange part
c
 9210  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' sogga11-x "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(54) = .true.
         else
c
c           assume total wanted
c
            lxfac(54) = .true.
            nlxfac(54) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(54) = .true.
         nlxfac(54) = .true.
      endif
      if (.not. inp_f(xfac(54))) then
         xfac(54) = 1.0d0
      endif
cadded HF X bit
        xfac(1) = 0.4015d0*xfac(54)
       lxfac(1) = .true.
      nlxfac(1) = .true.
       goto 10
c
c SOGGA11-X correlation part
c
 9220 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
c            write(LuOut,*)' cm11 is a semi-local functional;',
c     &                ' "local" descriptor is ignored. '



            write(LuOut,*)' sogga11-x "local" is not defined',
     &                ' for the metaGGA '
            write(LuOut,*)' Resubmit without local descriptor. '
            call errquit(' xc_inp - functional description ',0,0)


c
c
c
            lcfac(54) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)'csogga11 is already a nonlocal functional;',
     &                ' "nonlocal" descriptor is ignored. '
            nlcfac(54) = .true.
         else
c
c           assume total wanted (but total here = nonlocal)
c
            lcfac(54) = .false.
            nlcfac(54) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(54) = .false.
         nlcfac(54) = .true.
      endif
            write(LuOut,*)'csogga11-x uses PW91LDA as',
     &                ' defaults.These defaults cannot be changed '

      if (.not. inp_f(cfac(54))) then
         cfac(54) = 1.0d0
      endif
      goto 10  
c
c  N12 xc
c
 9245 xccomb(45) = .true.
         xfac(45)= 1.0d0
        lxfac(45)=.true.
       nlxfac(45)=.true.
        cfac(45) = 1.0d0
       lcfac(45) = .true.
      nlcfac(45) = .true.
c
      goto 10

c
c  N12 x
c
 9246  continue
        xfac(45) = 1.0d0
       lxfac(45) = .true.
      nlxfac(45) = .true.
c
      goto 10

c
c  N12 c
c
 9247 continue
         cfac(45)= 1.0d0
        lcfac(45)=.true.
       nlcfac(45)=.true.
c
      goto 10

c
c  N12-SX
c
 9260 xccomb(82) = .true. 
c
c       Ex(HF-SR)
c
        xfac(1) = 0.25d0
       lxfac(1) = .true.
       nlxfac(1) = .true.
c       
        xfac(82)= 1.0d0
        lxfac(82)=.true.
        nlxfac(82)=.true.
        cfac(82) = 1.0d0
       lcfac(82) = .true.
       nlcfac(82) = .true.
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.11d0
      cam_alpha = 0.d0
      cam_beta =  1.0
      cam_exch = .true.
      cam_srhf = .true. ! for short-range HF
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &     cam_beta)
c
      goto 10
       

c     
c     becke86b; Perdew-Burke-Ernzerhof exchange functional 
c     
 9510  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xbecke86b "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(55) = .true.
         else
c
c           assume total wanted
c
            lxfac(55) = .true.
            nlxfac(55) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(55) = .true.
         nlxfac(55) = .true.
      endif
      if (.not. inp_f(xfac(55))) then
         xfac(55) = 1.0d0
      endif

      goto 10
c     
c     xperdew86: Perdew-Wang 86 exchange
c     
 9520  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xperdew86 exchange "local" not possible.'
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' xperdew86 exchange "nonlocal" not possible.'
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         else
c           total wanted
            lxfac(56) = .true.
            nlxfac(56) = .true.
            call inp_prev_field()
         endif
      else
c        total wanted
         lxfac(56) = .true.
         nlxfac(56) = .true.
      endif
      if (.not. inp_f(xfac(56))) then
         xfac(56) = 1.0d0
      endif
c
      goto 10  
c
c     xbnl07: Long-range corrected Baer-Neuhauser-Livshits (BNL) 2007 exchange functional
c     
 8200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            lxfac(38) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlxfac(38) = .false.
            write(LuOut,*)' BNL is a local functional;',
     &                ' "nonlocal" descriptor is ignored. '
         else
c
c           assume total wanted
c
            lxfac(38) = .true.
            nlxfac(38) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(38) = .true.
         nlxfac(38) = .false.
      endif
      if (.not. inp_f(xfac(38))) then
         xfac(38) = 1.0d0
      endif
c
      goto 10
c
c     xcamb88: Coulomb-attenuating Becke88 exchange functional
c
 8300 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' becke88 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
            lxfac(40) = .true.
            nlxfac(40) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlxfac(40) = .true.
         else
c
c           assume total wanted
c
            lxfac(40) = .true.
            nlxfac(40) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(40) = .true.
         nlxfac(40) = .true.
      endif
      if (.not. inp_f(xfac(40))) then
         xfac(40) = 1.0d0
      endif
c
      goto 10
c     
c     xcampbe96; CAM-Perdew-Burke-Ernzerhof 96 functional for exchange
c     
 8400 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' cam-xpbe96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(41) = .true.
         else
c
c           assume total wanted
c
            lxfac(41) = .true.
            nlxfac(41) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(41) = .true.
         nlxfac(41) = .true.
      endif
      if (.not. inp_f(xfac(41))) then
         xfac(41) = 1.0d0
      endif
c
      goto 10
c
c     ssb-d; spin-state and sn2 barrier based functional
c            by Swart, Sola, Bickelhaupt, JCP 131, 094103 (2009)
c            includes Grimme's dispersion correction
c        
 4600 xccomb(46) = .true.
        xfac(46) = 1.d0
      nlxfac(46) = .true.
        xfac(47) = 0.358830d0
      nlxfac(47) = .true.
        cfac(46) = 1.d0
      nlcfac(46) = .true.
       lxfac(2)  = .true.
        xfac(2)  = 1.079966d0
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
c
c     this functional must be used with vdw contrib
c     (Grimme) with factor 0.847455
c
      goto 10
c
c     s12g functional
c
 4610 xccomb(60) = .true.
        xfac(60) = 1.d0
      nlxfac(60) = .true.
        xfac(2)  = 1.03842032d0
       lxfac(2)  = .true.
c
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
      nlcfac(12) = .true.
        cfac(12) = 1.0d0
c
c       this functional must be used with vdw contrib
c       (Grimme3) with factors 1.17755954 & 0.84432515
c
      ivdw = 3
      if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, ivdw))
     &    call errquit('xc_inp: rtdb_put failed vdw', 0, RTDB_ERR)
c
      goto 10
c
c     s12h functional
c
 4620 xccomb(61) = .true.
         xfac(1) = 0.25d0
        lxfac(1) = .true.
       nlxfac(1) = .true.
        xfac(61) = 0.75d0
      nlxfac(61) = .true.
        xfac(2)  = 1.02543951d0*xfac(61)
       lxfac(2)  = .true.
c
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
      nlcfac(12) = .true.
        cfac(12) = 1.0d0
c
c       this functional must be used with vdw contrib
c       (Grimme3) with factors 1.07735222 & 0.37705816
c
      ivdw = 3
      if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, ivdw))
     &    call errquit('xc_inp: rtdb_put failed vdw', 0, RTDB_ERR)
c
      goto 10
c
 4630 xccomb(62) = .true.
c
c     cam-s12g
c
      cam_omega = 1.52420731d0
      cam_alpha = 0.0d0
      cam_beta  = 0.34485046d0
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &    cam_beta)
c
           xfac(1) = 1.0d0
          lxfac(1) = .true.
         nlxfac(1) = .true.
          xfac(62) = 1.0d0
         lxfac(62) = .true.
        nlxfac(62) = .true.
c
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
      nlcfac(12) = .true.
        cfac(12) = 1.0d0
c
c       this functional must be used with vdw contrib
c       (Grimme3) with factors 1.20250451 & 0.86124355
c
      ivdw = 3
      if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, ivdw))
     &    call errquit('xc_inp: rtdb_put failed vdw', 0, RTDB_ERR)
c
      goto 10
c
 4640 xccomb(63) = .true.
c
c     cam-s12h
c
      cam_omega = 0.48516891d0
      cam_alpha = 0.25d0
      cam_beta  = 0.10897845d0
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &    cam_beta)
c
           xfac(1) = 1.0d0
          lxfac(1) = .true.
         nlxfac(1) = .true.
          xfac(63) = 1.0d0
         lxfac(63) = .true.
        nlxfac(63) = .true.
c
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
      nlcfac(12) = .true.
        cfac(12) = 1.0d0
c
c       this functional must be used with vdw contrib
c       (Grimme3) with factors 1.08034183 & 0.37999939
c
      ivdw = 3
      if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, ivdw))
     &    call errquit('xc_inp: rtdb_put failed vdw', 0, RTDB_ERR)
c
      goto 10
c
c     kt1: Keal-Tozer functional nr. 1
c
 4700   xfac(47) = 1.d0
      nlxfac(47) = .true.
       lxfac(2)  = .true.
        xfac(2)  = 1.0d0
       lcfac(6)  = .true.
        cfac(6)  = 1.0d0
c
      goto 10
c
c     kt1: Keal-Tozer functional nr. 2
c
 4701 xccomb(47) = .true.
        xfac(47) = 1.d0
      nlxfac(47) = .true.
       lxfac(2)  = .true.
        xfac(2)  = 1.07173d0
       lcfac(6)  = .true.
        cfac(6)  = 0.576727d0
c
      goto 10
c     
c     xcamlsd; CAM-LSD for exchange
c     
 8520 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' cam-xlsd "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(42) = .false.
         else
c
c           assume total wanted
c
            lxfac(42) = .true.
            nlxfac(42) = .false.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(42) = .true.
         nlxfac(42) = .false.
      endif
      if (.not. inp_f(xfac(42))) then
         xfac(42) = 1.0d0
      endif
c
      goto 10
c
c     xwpbe; HSE screened exchange
c
 8530 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xwpbe "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(43) = .true.
         else
c
c           assume total wanted
c
            lxfac(43) = .true.
            nlxfac(43) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(43) = .true.
         nlxfac(43) = .true.
      endif
      if (.not. inp_f(xfac(43))) then
         xfac(43) = 1.0d0
      endif
c
      goto 10
c
c     mp2; MP2 ``functional'' for correlation
c
 8500 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test) ) then
            write(LuOut,*)' MP2 correlation is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total MP2 correlation is assumed desired. '
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' MP2 correaltion is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total MP2 correlation is assumed desired. '
         else
            call inp_prev_field()
         endif
      endif
      if (.not. inp_f(mp2fac) ) then
         write(LuOut,*)' xc_inp: mp2 multiplicative factor not found.'
         mp2fac = 1.0d0
      endif
      if (.not. rtdb_put(rtdb,'dft:mp2fac', mt_dbl, 1, mp2fac))
     1     call errquit('xc_inp: failed to store dft:mp2fac in rtdb',
     2                  0,RTDB_ERR)
c
      goto 10
c
c    b2plyp; B2PLYP fifth-rung double hybrid of Grimme 
c
c            S. Grimme JCP 124, 034108 (2006)
c
c            A*Ex(HF) + (1-A)*Ex(S) + B*Ex(B88(nl)) +
c            C*Ec(LYP(total)) + (1-C)*Ec(VWN_1_RPA)
c
c            A=0.20, B=0.72, C=0.81
c
c     b3lyp  = HFexch 0.20 becke88 nonlocal 0.72 lyp 0.81          vwn_1_rpa 0.19 slater 0.80
c     b2plyp = HFexch 0.53 becke88 0.47          lyp 0.73 mp2 0.27
c
 8510 xccomb(2) = .true.
c HF exchange
        xfac(1) = 0.53d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c Becke88
        xfac(3) = 0.47d0
       lxfac(3) = .true.
      nlxfac(3) = .true.
c LYP
        cfac(2) = 0.73d0
       lcfac(2) = .true.
      nlcfac(2) = .true.
c MP2
      mp2fac = 0.27d0
      if (.not. rtdb_put(rtdb,'dft:mp2fac', mt_dbl, 1, mp2fac))
     1     call errquit('xc_inp: failed to store dft:mp2fac in rtdb',
     2                  0,RTDB_ERR)
c
      goto 10

c     dlDF      
 9300 xccomb(31) = .true.
        xfac(32) = 0.3855871d0
       lxfac(32) = .true.
      nlxfac(32) = .true.
        cfac(42) = 1.0d0
       lcfac(42) = .true.
      nlcfac(42) = .true.
         xfac(1) = 0.6144129d0
        lxfac(1) = .true.
       nlxfac(1) = .true.
c
      goto 10
c
c     SRHFexch; Short-Range HF exact exchange
c
 9310 if (inp_a(test) ) then
         if (inp_compare(.false.,'local', test) ) then
            write(LuOut,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total HF exchange is assumed desired. '
            lxfac(1) = .true.
            nlxfac(1) = .true.
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(LuOut,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(LuOut,*)' into local and nonlocal pieces;',
     &                ' total HF exchange is assumed desired. '
            lxfac(1) = .true.
            nlxfac(1) = .true.
         else
c
c           assume total wanted
c
            lxfac(1) = .true.
            nlxfac(1) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(1) = .true.
         nlxfac(1) = .true.
      endif
      if (.not. inp_f(xfac(1)) ) then
         write(LuOut,*)' xc_inp: hfexch multiplicative factor not ',
     &                 'found.'
         xfac(1) = 1.0d0
      endif
c
      cam_exch = .true.
      cam_srhf = .true.
      if (.not.rtdb_put(rtdb,'dft:cam_exch', mt_log, 1, cam_exch))
     & call errquit('xc_inp: rtdb_put failed', 2700, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'dft:cam_srhf', mt_log, 1, cam_srhf))
     & call errquit('xc_inp: rtdb_put failed', 2700, RTDB_ERR)
c
      goto 10
c
c     hse03; HSE03
c     0.25*Ex(HF-SR) -0.25*Ex(PBE-SR) + Ex(PBE) + Ec(PBE) 
c
 9320 xccomb(59) = .true.
c
c       Ex(HF-SR)
c
        xfac(1) = 0.25d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
c       Ex(PBE-SR)
c
        xfac(41) = -0.25
       lxfac(41) = .true.
      nlxfac(41) = .true.
c
c       Ex(PBE) + Ec(PBE)
c
        xfac(7)  =  1.0d0
       lxfac(7)  = .true.
      nlxfac(7)  = .true.
        cfac(12) =  1.0d0
       lcfac(12) = .false.
      nlcfac(12) = .true.
       lcfac(6) = .true.
       cfac(6) = 1.0d0
c    
c     set cam parameters since this is a range-separated functional 
c
      cam_omega = 0.33
      cam_alpha = 0.d0
      cam_beta =  1.0
      cam_exch = .true.
      cam_srhf = .true. ! for short-range HF
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &   cam_beta)
c
      goto 10
c
c
c wb97
c
 9410  xccomb(76) = .true.
         xfac(76)= 1.0d0
        lxfac(76)=.true.
       nlxfac(76)=.true.
        cfac(76) = 1.d0
       lcfac(76) = .true.
       nlcfac(76) = .true.
c
c       Ex(HF-LR)=Ex(HF)-Ex(HF-SR)
c
        xfac(1) = 1d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
c     set cam parameters since this is a range-separated functional
c
      cam_omega = 0.4d0
      cam_alpha = 0d0
      cam_beta = 1d0
      cam_exch = .true.
      cam_srhf = .false.

      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &     cam_beta)

c
      goto 10
c
c wb97x
c
 9420 xccomb(77) = .true.
         xfac(77)= 1.0d0
        lxfac(77)=.true.
       nlxfac(77)=.true.
        cfac(77) = 1.d0
       lcfac(77) = .true.
      nlcfac(77) = .true.
        xfac(1) = 1d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
c     set cam parameters since this is a range-separated functional
c
      cam_omega = 0.3d0
      cam_alpha = 0.157706d0
      cam_beta =  0.842294d0
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &    cam_beta)
c
      goto 10
c
c wb97x-d3
c
 9430 xccomb(78) = .true.
         xfac(78)= 1.0d0
        lxfac(78)=.true.
       nlxfac(78)=.true.
        cfac(78) = 1.d0
       lcfac(78) = .true.
      nlcfac(78) = .true.
        xfac(1) = 1d0
       lxfac(1) = .true.
      nlxfac(1) = .true.
c
c     set cam parameters since this is a range-separated functional
c
      cam_omega = 0.25d0
      cam_alpha = 0.195728d0
      cam_beta =  1d0-cam_alpha
      cam_exch = .true.
      cam_srhf = .false.
      call xc_setcamparam(rtdb,cam_exch,cam_srhf,cam_omega,cam_alpha,
     &     cam_beta)
c D3
      ivdw = 3
      if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, ivdw))
     &    call errquit('xc_inp: rtdb_put failed vdw', 0, RTDB_ERR)
c
      goto 10
c
c     gill96; Gill '96 functional for exchange
c     
 2300  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' gill96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(8) = .true.
         else
c
c           assume total wanted
c
            lxfac(8) = .true.
            nlxfac(8) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(8) = .true.
         nlxfac(8) = .true.
      endif
      if (.not. inp_f(xfac(8))) then
         xfac(8) = 1.0d0
      endif
       goto 10
c     
c     xperdew91; Perdew '91 functional for exchange
c     
 2400  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xperdew91 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0,
     &       INPUT_ERR)
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            nlxfac(9) = .true.
         else
c
c           assume total wanted
c
            lxfac(9) = .true.
            nlxfac(9) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lxfac(9) = .true.
         nlxfac(9) = .true.
      endif
      if (.not. inp_f(xfac(9))) then
         xfac(9) = 1.0d0
      endif
      goto 10
c     
c     OP; OP Tsuneda Suzumura Hirao functional for correlation
c             
 1997 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' OP is a nonlocal functional;',
     &           ' "local" descriptor is ignored. '                   
            call inp_prev_field()
c     
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' The OP functional cannot be conveniently',
     &           ' split into local and nonlocal pieces; '
            write(LuOut,*)' total OP is assumed desired. '
            call inp_prev_field()
         endif
      endif
      lcfac(36) = .true.
      nlcfac(36) = .true.
      if (.not. inp_f(cfac(36))) then
         cfac(36) = 1.0d0
      endif            
c     
      goto 10
c     
c     BOP Tsuneda Suzumura Hirao functional for XC
c             
 2046 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' OP is a nonlocal functional;',
     &           ' "local" descriptor is ignored. '                   
            call inp_prev_field()
c     
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' The OP functional cannot be conveniently',
     &           ' split into local and nonlocal pieces; '
            write(LuOut,*)' total OP is assumed desired. '
            call inp_prev_field()
         endif
      endif
      lcfac(36) = .true.
      nlcfac(36) = .true.
      lxfac(3) = .true.
      nlxfac(3) = .true.
      if (.not. inp_f(cfac(36))) then
         cfac(36) = 1.0d0
      endif            
      xfac(3)=cfac(36)
      xccomb(27) = .true.
c     
      goto 10
c     
c     PBEOP Tsuneda Suzumura Hirao functional for XC
c             
 2047 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' OP is a nonlocal functional;',
     &           ' "local" descriptor is ignored. '                   
            call inp_prev_field()
c     
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(LuOut,*)' The OP functional cannot be conveniently',
     &           ' split into local and nonlocal pieces; '
            write(LuOut,*)' total OP is assumed desired. '
            call inp_prev_field()
         endif
      endif
      lcfac(37) = .true.
      nlcfac(37) = .true.
      lxfac(7) = .true.
      nlxfac(7) = .true.
      if (.not. inp_f(cfac(37))) then
         cfac(37) = 1.0d0
      endif            
      xfac(7)=cfac(37)
      xccomb(28) = .true.
c     
      goto 10
 1999 continue
c
c     Check nonlocal COR functionals:
c     perdew86, perdew91, pbe96, cvpbe96
c     to see if user has specified an accompanying local
c     functional.  If not, set appropriately.
c
      ck_local = .false.
      if(nlcfac(4).or.nlcfac(5).or.nlcfac(12).or.nlcfac(64))then
        do n = 1, numfunc
           if(lcfac(n))ck_local = .true.
        enddo
c
c       if ck_local is true than local functional has
c       been specified.  If false than a default local
c       functional is needed.
c
        if(.not.ck_local)then
          if(nlcfac(4))then
            write(LuOut,*)' perdew86 is a nonlocal functional;',
     &                ' adding perdew81 local functional. '
            lcfac(3) = .true.
            cfac(3) = cfac(4)
          elseif(nlcfac(5))then
            write(LuOut,*)' perdew91 is a nonlocal functional;',
     &                ' adding pw91lda local functional. '
c
c           This is also done by default inside xc_perdew91.F
c           when lcfac(5) is .true.
c
            lcfac(6) = .true.
            cfac(6) = 1.0d0
          elseif(nlcfac(12))then
            write(LuOut,*)' pbe96 is a nonlocal functional;',
     &                ' adding pw91lda local functional. '
c
c           This is also done by default inside xc_perdew91.F
c           when lcfac(5) is .true.
c
            lcfac(6) = .true.
            cfac(6) = 1.0d0

          elseif(nlcfac(64))then
            write(LuOut,*)' cvpbe09 is a nonlocal functional;',
     &                ' adding pw91lda local functional. '
c
            lcfac(6) = .true.
            cfac(6) = 1.0d0

          endif
        endif
      endif
      
c     
c    
c      write(*,*) ' lcfac '
c      write(*,*) (lcfac(iline),iline=1,numfunc) 
c      write(*,*) ' nlcfac '
c      write(*,*) (nlcfac(iline),iline=1,numfunc) 
c      write(*,*) ' lxfac '
c      write(*,*) (lxfac(iline),iline=1,numfunc) 
c      write(*,*) ' nlxfac '
c      write(*,*) (nlxfac(iline),iline=1,numfunc) 
c      write(*,*) ' xccomb '
c      write(*,*) (xccomb(iline),iline=1,numfunc) 
c      write(*,*) ' cfac '
c      write(*,*) (cfac(iline),iline=1,numfunc) 
c      write(*,*) ' xfac '
c      write(*,*) (xfac(iline),iline=1,numfunc) 
c
c     No alternative but to dump all array's associated
c     with type of XC into rtdb (whether specified or not).
c     
c>>>  place logical XC parameters in rtdb
c
cc AJL/Begin/FDE
c
c    This has been made slightly more dynamic so we can get information
c    in to the RTDB for e.g. frozen density embedding
c
c      if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, numfunc, lcfac))
c     &   call errquit('xc_inp: rtdb_put failed', 1, RTDB_ERR)
c      if (.not. rtdb_put(rtdb, 'dft:nlcfac', mt_log, numfunc, nlcfac))
c     &   call errquit('xc_inp: rtdb_put failed', 2, RTDB_ERR)
c      if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, numfunc, lxfac))
c     &   call errquit('xc_inp: rtdb_put failed', 3, RTDB_ERR)
c      if (.not. rtdb_put(rtdb, 'dft:nlxfac', mt_log, numfunc, nlxfac))
c     &   call errquit('xc_inp: rtdb_put failed', 4, RTDB_ERR)
c      if (.not. rtdb_put(rtdb, 'dft:xccomb', mt_log, numfunc, xccomb))
c     &   call errquit('xc_inp: rtdb_put failed', 5, RTDB_ERR)
c     
c>>>  place real XC parameters in rtdb
c     
c      if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, numfunc, cfac))
c     &   call errquit('xc_inp: rtdb_put failed', 6, RTDB_ERR)
c      if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
c     &   call errquit('xc_inp: rtdb_put failed', 7, RTDB_ERR)

      if (.not. rtdb_put(rtdb, module(1:mlen)//':lcfac', 
     &                   mt_log, numfunc, lcfac))
     &   call errquit('xc_inp: rtdb_put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, module(1:mlen)//':nlcfac', 
     &                   mt_log, numfunc, nlcfac))
     &   call errquit('xc_inp: rtdb_put failed', 2, RTDB_ERR)
      if (.not. rtdb_put(rtdb, module(1:mlen)//':lxfac', 
     &                   mt_log, numfunc, lxfac))
     &   call errquit('xc_inp: rtdb_put failed', 3, RTDB_ERR)
      if (.not. rtdb_put(rtdb, module(1:mlen)//':nlxfac', 
     &                   mt_log, numfunc, nlxfac))
     &   call errquit('xc_inp: rtdb_put failed', 4, RTDB_ERR)
      if (.not. rtdb_put(rtdb, module(1:mlen)//':xccomb', 
     &                   mt_log, numfunc, xccomb))
     &   call errquit('xc_inp: rtdb_put failed', 5, RTDB_ERR)
c     
c>>>  place real XC parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, module(1:mlen)//':cfac', 
     &                   mt_dbl, numfunc, cfac))
     &   call errquit('xc_inp: rtdb_put failed', 6, RTDB_ERR)
      if (.not. rtdb_put(rtdb, module(1:mlen)//':xfac',
     &                   mt_dbl, numfunc, xfac))
     &   call errquit('xc_inp: rtdb_put failed', 7, RTDB_ERR)

cc AJL/End

      if (is_libxcon()) then
        if (.not.(rtdb_put(rtdb, module(1:mlen)//':libxcon',
     $                     mt_log,1,.true.)))
     $   call errquit('xc_inp: rtdb_put failed', 8, RTDB_ERR)
        call nwchem_libxc_rtdbput(rtdb,module(1:mlen))
        else
        if (.not.(rtdb_put(rtdb, module(1:mlen)//':libxcon',
     $                     mt_log,1,.false.)))
     $   call errquit('xc_inp: rtdb_put failed', 9, RTDB_ERR)
      endif   
c
      return
c     
10000 write(luout,777)
 777  format(/,' xc [acm],')
      do n=2,num_dirs
         write(luout,778) dirs(n)
 778     format(' [',A10,'],')
      enddo

      call util_flush(LuOut)
      call errquit('xc_input: invalid format', 0, INPUT_ERR)
c     
      end
      subroutine sic_input(rtdb, module)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the sic directive which specifies how to construct the
c     sic/oep approach is used.
c     
c     Possible variables are:
c
c     perturbative
c     oep
c     oep-loc
c     
      integer num_dirs, ind, mlen, test_sic
      parameter (num_dirs = 3)
      character*12 dirs(num_dirs)
      character*255 test
c
c
      data dirs /'perturbative', 'oep', 'oep-loc'/
c     
      mlen = inp_strlen(module)
      test_sic = 1
c     
 10   if (.not.inp_a(test)) goto 1999
c     
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c     
c        Does not match a keyword ... 
c     
         goto 10000
      endif
c
      goto (100, 200, 300, 1999) ind
      call errquit('sic_inp: unimplemented directive', ind, INPUT_ERR)
c     
  100 if (.not.inp_a(test)) then
        test_sic = 1
      else
        goto 10000
      end if            
c
      goto 10
  200 if (.not.inp_a(test)) then
        test_sic = 2
      else
        goto 10000
      end if
c
      goto 10
  300 if (.not.inp_a(test)) then
        test_sic = 4
      else
        goto 10000
      end if
c
      goto 10
c             
 1999 continue
c
      if (.not. rtdb_put(rtdb, 'dft:test_sic', mt_int, 1,test_sic))
     &     call errquit('sic_inp: rtdb_put failed', 1, RTDB_ERR)
c     
      return
c     
10000 write(LuOut,10001)

10001 format(/,' sic [[perturbative], [oep], [oep-loc]]')
      call util_flush(LuOut)
      call errquit('sic_input: invalid format', 0, INPUT_ERR)
c     
      end
      subroutine forcedirect(rtdb)
      implicit none

      integer rtdb

#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer memsize, filesize
c
c     Force use of direct
c
      direct = .true.
      filesize = -1
      memsize  = -1
      if (.not. rtdb_put(rtdb, 'dft:direct',mt_log,1,direct))
     &   call errquit('cam-s12h: rtdb_put failed',400,RTDB_ERR)
      if (.not. rtdb_put(rtdb,'int2e:filesize',mt_int,1,filesize))
     $     call errquit('cam-s12h: failed to store filesize',401,
     &                  RTDB_ERR)
      if (.not. rtdb_put(rtdb,'int2e:memsize',mt_int,1,memsize))
     $     call errquit('cam-s12h: failed to store filesize',402,
     &                  RTDB_ERR)
      end
c
      subroutine ts_input(rtdb, module)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the kinetic energy directive which specifies how to deal 
c     with kinetic energy for frozen-density embedding calculations
c     
c     Possible variables are currently:
c
c     lda:
c     Thomas-Fermi (tf) (-> 100)
c     General Expansion Approx., 2nd Order (gea0) (-> 100)
c     pbe:
c     von-Weizsacker (vw) (-> 200)
c     General Expansion Approx., 2nd Order (gea2) (-> 300)
c     Thomas-Fermi plus von Weisacker (tfw) (-> 400)
c     
      integer num_dirs, ind, mlen
      parameter (num_dirs = 2)
      double precision tsfac(num_dirs)
      character*12 dirs(5)
      character*255 test
c
c
      data dirs /'tf', 'vw', 'gea0', 'gea2', 'tfw'/
c
c Benchmarking of combinatorial functionals: Ts atomic energies
c Previous results in brackets, our results following.
c See the TF subroutine for more detailed discussion of test settings.
c
c - spin-paired: Jones and Gunnarsson, 10.1103/RevModPhys.61.689
c --> He : GEA2 (2.8785), 2.8778 
c --> Ne : GEA2 (127.79), 127.7898
c
c - spin-unpaired: Iyrngar et al., 10.1103/PhysRevA.63.052508,
c --> H: GEA2 (0.511) , 0.512
c        TFW  (0.956) , 0.972
c --> C: GEA2 (37.56) , 37.639
c        TFW  (66.179), 66.276
c     
      mlen = inp_strlen(module)
      do ind=1,num_dirs
        tsfac(ind) = 0.0d0
      enddo
c     
 10   if (.not.inp_a(test)) goto 1999
c     
      if (.not. inp_match(5, .false., test, dirs, ind)) then
c     
c        Does not match a keyword ... 
c     
         goto 10000
      endif
c
      goto (100, 200, 100, 300, 400, 1999) ind
      call errquit('ts_input: unimplemented directive', ind, INPUT_ERR)
c     
  100 if (.not. inp_f(tsfac(1))) then
        tsfac(1) = 1.0d0
      endif
c
      goto 10
  200 if (.not. inp_f(tsfac(2))) then
        tsfac(2) = 1.d0
      end if
c
      goto 10
c
  300 continue
      tsfac(1) = 1.d0
      tsfac(2) = 1.d0/9.d0
c
      goto 10
c
  400 continue
      tsfac(1) = 1.d0
      tsfac(2) = 1.d0
c
      goto 10
c             
 1999 continue
c
      if (.not. rtdb_put(rtdb, module(1:mlen)//':ts', mt_dbl, num_dirs,
     &          tsfac))
     &     call errquit('ts_input: rtdb_put failed', 1, RTDB_ERR)
c     
      return
c     
10000 write(LuOut,10001)

10001 format(/,' frozemb_ts [[tf], [vw], [gea0], [gea2], [tfw]]')
      call util_flush(LuOut)
      call errquit('ts_input: invalid format', 0, INPUT_ERR)
c     
      end
      subroutine xcinp_fact_noloc(test,funcname,fact_out)
      implicit none
      character*(*) test        ! [in]
      character*(*) funcname    ! [in]
      double precision fact_out ! [out]
#include "inp.fh"
#include "stdio.fh"
 1    if (.not. inp_f(fact_out)) then
         if (inp_a(test)) then
            if (inp_compare(.false.,'local', test)) then
               write(LuOut,*) funcname(1:inp_strlen(funcname)),
     &              ' is a nonlocal functional;',
     &              ' "local" descriptor is ignored. '
               goto 1
            elseif (inp_compare(.false.,'nonlocal', test)) then
               write(LuOut,*) funcname(1:inp_strlen(funcname)),
     &              ' cannot be conveniently',
     &              ' split into local and nonlocal pieces; '
               write(LuOut,*)' total ',funcname(1:inp_strlen(funcname)),
     &              ' is used. '
               goto 1
            endif
         else
            fact_out=1d0
         endif
      endif
      return
      end
