      subroutine xc_input(rtdb, module)
C     $Id: 
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
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
c     
      integer num_dirs, ind, mlen, iline
      parameter (num_dirs = 17)
      character*10 dirs(num_dirs)
      character*255 test
c
c     Define logical flags for local (l) and nonlocal (nl)
c     parts of exchange-correlation functionals.
c
      logical  lcfac(15),  lxfac(15)
      logical nlcfac(15), nlxfac(15)
      logical xccomb(15)
c
c     XC functionals coefficients
c     
      double precision cfac(15),xfac(15)
      data dirs /'acm', 'b3lyp', 'beckehandh', 'HFexch', 'becke88',  
     &           'lyp', 'perdew81', 'perdew86', 'perdew91', 'pw91lda', 
     &           'slater', 'vwn_1', 'vwn_2', 'vwn_3', 'vwn_4', 'vwn_5',
     &           'vwn_1_rpa'/
c
c     Only these few arrays are initialized (not defaults).
c
      do iline = 1, 15
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
         goto 10
      endif
c
      write(6,*)' test = ',test     
      write(6,*)' ind = ',ind     
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 
     &      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1999) ind
      call errquit('xc_inp: unimplemented directive', ind)
c     
c     acm; adiabatic connection method
c          AD Becke JCP 98, 5648 (1993)
c          Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + C*Ec(P91) 
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
c            Exc(S+VWN) + A*(Ex(HF)-Ex(S) + B*Ex(B88) + (1-C)*Ec(VWN)
c            + C*Ec(VWN) 
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
            write(6,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(6,*)' into local and nonlocal pieces;',
     &                ' total HF exchange is assumed desired. '                   
            lxfac(1) = .true.
            nlxfac(1) = .true.
         elseif (inp_compare(.false.,'nonlocal', test) ) then
            write(6,*)' HF exchange is by nature nonlocal',
     &                ' and cannot be conveniently split'
            write(6,*)' into local and nonlocal pieces;',
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
         write(6,*)' xc_inp: hfexch multiplicative factor not found.'
         xfac(1) = 1.0d0
      endif            
c
      goto 10
c     
c     becke88; Becke '88 functional for exchange
c     
  500 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(6,*)' becke88 "local" is equivalent to the',
     &                ' slater exchange local functional. '
            write(6,*)' Resubmit using slater exchange. '
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
c
      goto 10
c     
c     lyp; LYP functional for correlation
c             
  600 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(6,*)' The lyp functional cannot be conveniently',
     &                ' split into local and nonlocal pieces; '
            write(6,*)' total lyp is assumed desired. '                   
            lcfac(2) = .true.
            nlcfac(2) = .true.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            write(6,*)' The lyp functional cannot be conveniently',
     &                ' split into local and nonlocal pieces; '
            write(6,*)' total lyp is assumed desired. '                   
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
            write(6,*)' perdew81 is a local functional;',
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
            write(6,*)' perdew86 is a nonlocal functional;',
     &                ' adding perdew81 local functional. '                   
            lcfac(3) = .true.
            lcfac(4) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(3) = .false.
            nlcfac(4) = .true.
         else
c
c           assume total wanted
c
            write(6,*)' perdew86 is a nonlocal functional;',
     &                ' adding perdew81 local functional. '                   
            lcfac(3) = .true.
            nlcfac(3) = .false.
            lcfac(4) = .false.
            nlcfac(4) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         write(6,*)' perdew86 is a nonlocal functional;',
     &             ' adding perdew81 local functional. '                   
         lcfac(3) = .true.
         nlcfac(3) = .false.
         lcfac(4) = .false.
         nlcfac(4) = .true.
      endif
      if (.not. inp_f(cfac(4))) then
         cfac(4) = 1.0d0
         if(lcfac(3))cfac(3) = cfac(4)
      endif            
c
      goto 10
c     
c     perdew91; P91 functional for correlation
c             
  900 if (inp_a(test)) then
         if (inp_compare(.false.,'local', test)) then
            write(6,*)' perdew91 is a nonlocal functional;',
     &                ' adding pw91lda local functional. '
c
c           This is also done by default inside xc_perdew91.F
c           when lcfac(5) is .true.
c
            lcfac(6) = .true.
            lcfac(5) = .false.
         elseif (inp_compare(.false.,'nonlocal', test)) then
            nlcfac(5) = .true.
         else
c
c           assume total wanted
c
            lcfac(5) = .true.
            nlcfac(5) = .true.
            call inp_prev_field()
         endif
      else
c
c        assume total wanted
c
         lcfac(5) = .true.
         nlcfac(5) = .true.
      endif
      if (.not. inp_f(cfac(5))) then
         cfac(5) = 1.0d0
         if(lcfac(6))cfac(6) = cfac(5)
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
            write(6,*)' pw91lda is a local functional;',
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
            write(6,*)' slater is a local functional;',
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
            write(6,*)' vwn_1 is a local functional;',
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
            write(6,*)' vwn_2 is a local functional;',
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
            write(6,*)' vwn_3 is a local functional;',
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
            write(6,*)' vwn_4 is a local functional;',
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
            write(6,*)' vwn_5 is a local functional;',
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
            write(6,*)' vwn_1_rpa is a local functional;',
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
 1999 continue
c     
c
c     No alternative but to dump all array's associated
c     with type of XC into rtdb (whether specified or not).
c     
c>>>  place logical XC parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, 15, lcfac))
     &   call errquit('xc_inp: rtdb_put failed', 1)
      if (.not. rtdb_put(rtdb, 'dft:nlcfac', mt_log, 15, nlcfac))
     &   call errquit('xc_inp: rtdb_put failed', 2)
      if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, 15, lxfac))
     &   call errquit('xc_inp: rtdb_put failed', 3)
      if (.not. rtdb_put(rtdb, 'dft:nlxfac', mt_log, 15, nlxfac))
     &   call errquit('xc_inp: rtdb_put failed', 4)
      if (.not. rtdb_put(rtdb, 'dft:xccomb', mt_log, 15, xccomb))
     &   call errquit('xc_inp: rtdb_put failed', 5)
c     
c>>>  place real XC parameters in rtdb
c     
      if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, 15, cfac))
     &   call errquit('xc_inp: rtdb_put failed', 6)
      if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, 15, xfac))
     &   call errquit('xc_inp: rtdb_put failed', 7)
      return
c     
10000 write(6,10001)

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
     &       /,' [vwn_1_rpa local rval]] ')
      call util_flush(6)
      call errquit('grid_input: invalid format', 0)
c     
      end
