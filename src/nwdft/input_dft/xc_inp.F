      subroutine xc_input(rtdb, module)
C     $Id: 
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
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
c     cpbe96
c     xpbe96
c     
      integer num_dirs, ind, mlen, iline, n,numfunc
      parameter (num_dirs = 36,numfunc=20)
      character*15 dirs(num_dirs)
      character*255 test
c
c     Define logical flags for local (l) and nonlocal (nl)
c     parts of exchange-correlation functionals.
c
      logical  lcfac(numfunc),  lxfac(numfunc), ck_local
      logical nlcfac(numfunc), nlxfac(numfunc)
      logical xccomb(numfunc)
c
c     XC functionals coefficients
c     
      double precision cfac(numfunc),xfac(numfunc)
      data dirs /'acm', 'b3lyp', 'beckehandh', 'HFexch', 'becke88',  
     &           'lyp', 'perdew81', 'perdew86', 'perdew91', 'pw91lda', 
     &           'slater', 'vwn_1', 'vwn_2', 'vwn_3', 'vwn_4', 
     &           'vwn_5','vwn_1_rpa','hcth','becke97','becke97-1',
     &           'cpbe96','xpbe96','gill96','xperdew91',
     ,     'xhcth','xbecke97','xbecke97-1',
     ,     'chcth','cbecke97','cbecke97-1','hcth120','hcth147',
     ,     'becke98','pbe0','becke97gga1','hcth402'/
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
     , 3000, 3100,
     ,     1999) ind
      call errquit('xc_inp: unimplemented directive', ind)
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
c                 0.5*Exc(LSD) + 0.5*Ex(HF) 
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
            call errquit(' xc_inp - functional description ',0)
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
       xccomb(8) = .true.
         lxfac(11) = .true.
         nlxfac(11) = .true.
         xfac(11) = 1.0d0
         lcfac(17) = .true.
         nlcfac(17) = .true.
         cfac(17) = 1.0d0
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
c     HCTH402 by Boese
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
c     cpbe96 Perdew-Burke-Ernzerhof functional for correlation
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
c     xpbe96; Perdew-Burke-Ernzerhof functional for exchange
c     
 2200 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' xpbe96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0)
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
c     gill96; Gill '96 functional for exchange
c     
 2300  if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(LuOut,*)' gill96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(LuOut,*)' Resubmit using slater exchange. '
            call errquit(' xc_inp - functional description ',0)
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
            call errquit(' xc_inp - functional description ',0)
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
 1999 continue
c
c     Check nonlocal functionals: perdew86 and perdew91 to
c     see if user has specified and accompanying local
c     functional.  If not, set appropriately.
c
      ck_local = .false.
      if(nlcfac(4).or.nlcfac(5).or.nlcfac(12))then
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
      if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, numfunc, lcfac))
     &   call errquit('xc_inp: rtdb_put failed', 1)
      if (.not. rtdb_put(rtdb, 'dft:nlcfac', mt_log, numfunc, nlcfac))
     &   call errquit('xc_inp: rtdb_put failed', 2)
      if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, numfunc, lxfac))
     &   call errquit('xc_inp: rtdb_put failed', 3)
      if (.not. rtdb_put(rtdb, 'dft:nlxfac', mt_log, numfunc, nlxfac))
     &   call errquit('xc_inp: rtdb_put failed', 4)
      if (.not. rtdb_put(rtdb, 'dft:xccomb', mt_log, numfunc, xccomb))
     &   call errquit('xc_inp: rtdb_put failed', 5)
c     
c>>>  place real XC parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, numfunc, cfac))
     &   call errquit('xc_inp: rtdb_put failed', 6)
      if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     &   call errquit('xc_inp: rtdb_put failed', 7)
      return
c     
10000 write(LuOut,10001)

10001 format(/,' xc [[acm], [b3lyp], [beckehandh], [HFexch rval], ',
     &       /,' [becke88 local/nonlocal rval], ',
     &       /,' [lyp local/nonlocal rval], ',
     &       /,' [perdew81 local rval], ',
     &       /,' [perdew86 nonlocal rval], ',
     &       /,' [perdew91 nonlocal rval], ',
     &       /,' [pw91lda local rval], ',
     &       /,' [slater local rval], ',
     &       /,' [vwn_1 local rval], ',
     &       /,' [vwn_2 local rval], ',
     &       /,' [vwn_3 local rval], ',
     &       /,' [vwn_4 local rval], ',
     &       /,' [vwn_5 local rval], ',
     &       /,' [vwn_1_rpa local rval]] ',
     &       /,' [becke97] ',
     &       /,' [becke97-1] ',
     &       /,' [cpbe96 nonlocal rval]',
     &       /,' [xpbe96] ',
     &       /,' [gill96] ',
     &       /,' [xperdew91] ',
     &       /,' [cbecke97] ',
     &       /,' [xbecke97] ',
     &       /,' [cbecke97-1',
     &       /,' [xbecke97-1',
     &       /,' [xhcth] ',
     &       /,' [chcth] ',
     &       /,' [hcth120] ',
     &       /,' [hcth147] ',
     &       /,' [hcth402] ',
     &       /,' [becke97gga1]')
      call util_flush(LuOut)
      call errquit('xc_input: invalid format', 0)
c     
      end
      subroutine sic_input(rtdb, module)
      implicit none
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
      call errquit('sic_inp: unimplemented directive', ind)
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
     &     call errquit('sic_inp: rtdb_put failed', 1)
c     
      return
c     
10000 write(LuOut,10001)

10001 format(/,' sic [[perturbative], [oep], [oep-loc]]')
      call util_flush(LuOut)
      call errquit('sic_input: invalid format', 0)
c     
      end



