      subroutine dft_input(rtdb)
c
c     $Id: dft_input.F,v 1.55 2009-02-06 00:25:25 niri Exp $
c
c     The main purpose of this subroutine is to parse the input
c     and place any "found" parameters into the rtdb.
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "case.fh"
c
c     read input for DFT calculation.  All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
c
      logical status
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 35)
      character*12 dirs(num_dirs)
      character*255 test
      logical disp
      integer vdw
      integer alpha
      integer s6
      integer nopen
      logical lcdft
      logical pseudo
      logical dftmp2
      logical molden
      logical noscf
      integer which_mp2
      character*20 mp2char
      double precision smear_sigma
      double precision shift
      integer nel_fon(2),nmo_fon(2),ispin
      logical steric
      logical efield
      double precision efield_x
      double precision efield_y
      double precision efield_z
c
      data dirs /'convergence','decomp','dft','direct','semidirect',
     $           'ecp','grid','vectors','incore','iterations', 
     $           'max_ovl','mulliken','mult','noio','noprint','odft', 
     $           'print','tolerances', 'xc','sym','adapt','sic', 
     $           'smear','fon','lb94','cs00','cdft','cam','disp',
     $           'dftmp2','molden','noscf','steric','efield','end'/
c
c Qin / cdft
      integer nconstr
      integer dtype(max_constr), iatms1(max_constr), fatms1(max_constr)
      integer iatms2(max_constr), fatms2(max_constr)
      double precision constr(max_constr)
      integer memsize, filesize
      nconstr = 0
c Qin / end
c
c     Any DFT input causes calculation to be marked as not converged
c
      status = rtdb_delete(rtdb,'dft:converged')
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &   ipol=1
      nel_fon(1)=0
      nmo_fon(1)=0
      nel_fon(2)=0
      nmo_fon(2)=0
c
c     Lebedev on by default
c
c
 10   if (.not. inp_read()) call errquit
     $     ('dft_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('dft_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('dft_input: unknown directive', 0,
     &       INPUT_ERR)
c     
      goto (100,  200,  300,  400,  450,  500,  600,  700,  800,
     $      900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
     $     1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 
     $     2700, 2800, 2900, 3000, 3100, 3200, 3300, 9999) ind

      call errquit('dft_input: unimplemented directive', ind,
     &       INPUT_ERR)
c
c     convergence;
c             
  100 call conv_input(rtdb, 'dft')
c
      goto 10
c             
c     decomp;
c             
  200 idecomp = 1
      if (.not. rtdb_put(rtdb, 'dft:idecomp', mt_int, 1, idecomp))
     &   call errquit('dft_input: rtdb_put failed', 200, RTDB_ERR)
c
      goto 10
c
c     dft;
c             
  300 ipol=1
c
      goto 10
c
c     direct;
c             
  400 direct = .true.
      filesize = -1
      memsize  = -1
      if (.not. rtdb_put(rtdb, 'dft:direct',mt_log,1,direct))
     &   call errquit('dft_input: rtdb_put failed',400,RTDB_ERR)
         if (.not. rtdb_put(rtdb,'int2e:filesize',mt_int,1,filesize))
     $        call errquit('dft_input: failed to store filesize',401,
     &                     RTDB_ERR)
         if (.not. rtdb_put(rtdb,'int2e:memsize',mt_int,1,memsize))
     $        call errquit('dft_input: failed to store filesize',402,
     &                     RTDB_ERR)
c
      goto 10
c
c     semidirect;
c
  450 direct = .false.
c
c     this part is copied directly from ddscf/scf_input.F
c
         filesize = -33
         memsize  = -33
 3601    if (inp_a(test)) then
            if (inp_compare(.false., test, 'filesize')) then
               if (.not. inp_i(filesize)) goto 3602
            else if (inp_compare(.false., test, 'memsize')) then
               if (.not. inp_i(memsize)) goto 3602
            else if (inp_compare(.false., test, 'filename')) then
               if (.not. inp_a(test)) goto 3602
               if (.not. rtdb_cput(rtdb, 'int2e:filename', 1, test))
     $              call errquit('dft_input: bad int2e:filename',450,
     $                           RTDB_ERR)
            else
               goto 3602
            endif
            goto 3601
         endif
         call inp_clear_err()
c
c     User should see a value of zero meaning zero not use default
c
         if (filesize.eq.0) filesize = -1
         if (memsize.eq.0) memsize = -1
         if (filesize.eq.-33) filesize=0
         if (memsize.eq.-33) memsize=0
         if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, filesize))
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, memsize))
     $        call errquit('scf_input: failed to store filesize',0,
     &       RTDB_ERR)
c
c     ecp;
c             
  500 pseudo = .true.
c
      goto 10
c
c     grid;
c             
  600 call grid_input(rtdb, 'dft')
c
      goto 10
c             
c     vectors; get input vectors, output vectors, initial guess, etc.
c             
  700 call vectors_input(rtdb,'dft')
c
      goto 10
c
c     incore; try to do an incore coulomb construction
c             
  800 incore = .true.
      if (.not.rtdb_put(rtdb,'dft:incore', mt_log, 1, incore)) 
     &   call errquit('dft_input: rtdb_put failed', 800, RTDB_ERR)
c
      goto 10
c
c     iterations;
c             
  900 if (.not. inp_i(iterations)) then
         write(LuOut,*)' dft_input: iterations value not found;',
     &             ' defaulting to 30'
         iterations = 30
      endif 
      if (.not. rtdb_put(rtdb, 'dft:iterations', mt_int, 1, iterations))
     &   call errquit('dft_input: rtdb_put failed', 900, RTDB_ERR)
c
      goto 10
c
c     max_ovl; Move orbitals around in movecs to maintain max 
c              overlap between iterations
c             
 1000 lmaxov = .true.
      if (.not. rtdb_put(rtdb, 'dft:max_ovl', mt_log, 1, lmaxov))
     &   call errquit('dft_input: rtdb_put failed', 1000, RTDB_ERR)
c
      goto 10
c
c     mulliken; do Mulliken pop analysis
c             
 1100 imull=1
      if (.not. rtdb_put(rtdb, 'dft:imull', mt_int, 1, imull))
     &   call errquit('dft_input: rtdb_put failed', 1100, RTDB_ERR)
c
      goto 10
c             
c     mult;
c             
 1200 if (.not. inp_i(mult)) then
         write(LuOut,*)' dft_input: spin multiplicity value not found;',
     &             ' ignoring option using default of singlet'
         mult = 1
      endif 
      if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
     &   call errquit('dft_input: rtdb_put failed', 1200, RTDB_ERR)
      nopen = mult - 1
      if (.not.rtdb_put(rtdb, 'dft:nopen',mt_int,1,nopen))
     &   call errquit('dft_input: rtdb_put failed', 1201, RTDB_ERR)
c
      goto 10
c             
c     noio;
c             
 1300 noio = 1
      if (.not. rtdb_put(rtdb, 'dft:noio', mt_int, 1, noio))
     &   call errquit('dft_input: rtdb_put failed', 1300, RTDB_ERR)
c
      goto 10
c
c     noprint;
c
 1400 call util_print_input(rtdb, 'dft')
c
      goto 10
c
c     odft;
c             
 1500 ipol=2
c
      goto 10
c
c     print;
c
 1600 call util_print_input(rtdb, 'dft')
c
      goto 10
c
c     tolerances;
c             
 1700 call tol_input(rtdb, 'dft')
c
      goto 10
c
c     xc;
c             
 1800 call xc_input(rtdb, 'dft')
c
      goto 10
c
c     sym [off|on] (defaults to on)
c
 1900 if (.not. inp_a(test)) test = 'on'
      status = inp_compare(.false.,test,'on')
      if (.not. rtdb_put(rtdb,'dft:skeleton', MT_LOG, 1, status))
     &   call errquit('dft_input: rtdb_put failed', 1900, RTDB_ERR)
      goto 10
c
c     adapt [off|on] (defaults to on)
c
 2000 if (.not. inp_a(test)) test = 'on'
      status = inp_compare(.false.,test,'on')
      if (.not. rtdb_put(rtdb,'dft:adapt', MT_LOG, 1, status))
     &   call errquit('dft_input: rtdb_put failed', 2000, RTDB_ERR)
      goto 10
c
 2100 call sic_input(rtdb, 'dft')
c
      goto 10
c
c     smearing
c
 2200 if (.not. inp_f(smear_sigma)) then
         write(LuOut,*)'  smearing value not found;',
     &             ' defaulting to 1d-3 Hartree.'
         smear_sigma = 1.d-3
      endif
      if (.not.rtdb_put(rtdb, 'dft:smear_sigma',mt_dbl,1,smear_sigma))
     &   call errquit('dft_inp: rtdb_put failed', 2200, RTDB_ERR)
      goto 10
c
c     fon
c
 2300 continue
      ispin=1
 2301 if (inp_a(test)) then
         if (inp_compare(.false.,'alpha', test)) then
            ispin=1
            goto 2301
         else if (inp_compare(.false.,'beta', test)) then
            ispin=2
            goto 2301
         endif
         if (inp_compare(.false.,'electrons', test)) then
            if(.not.inp_i(nel_fon(ispin))) call errquit(
     .           ' fon: required integer after electrons keyword',0,
     &       INPUT_ERR)
            goto 2301
         else if (inp_compare(.false.,'orbitals', test)) then
            if(.not.inp_i(nmo_fon(ispin))) call errquit(
     .           ' fon: required integer after orbitals keyword',0,
     &       INPUT_ERR)
            goto 2301
         endif
      endif
      if (.not.rtdb_put(rtdb, 'dft:nmo_fon', mt_int, 2,
     & nmo_fon)) call errquit('dft_inp: rtdb_put failed', 2300,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'dft:nel_fon', mt_int, 2,
     & nel_fon)) call errquit('dft_inp: rtdb_put failed', 2300,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'dft:fon', mt_log, 1,
     & .true.)) call errquit('dft_inp: rtdb_put failed', 2300,
     &       RTDB_ERR)
      goto 10
c
c     lb94; adds van-Leeuwen-Baerends correction to potential
c             
 2400 lb94 = .true.
      if (.not.rtdb_put(rtdb,'dft:lb94', mt_log, 1, lb94)) 
     &   call errquit('dft_input: rtdb_put failed', 2400, RTDB_ERR)
c
      goto 10
c
c     cs00; adds Casida-Salahub correction to potential
c             
 2500 cs00 = .true.
      if (.not.rtdb_put(rtdb,'dft:cs00', mt_log, 1, cs00)) 
     &   call errquit('dft_input: rtdb_put failed', 2500, RTDB_ERR)
      if (.not.inp_f(shift)) then
         write(LuOut,*)' shift for CS00 correction not found;',
     &   ' Zhan-Nichols-Dixon shift (B3LYP) will be used.'
         shift = 1.0d99
      endif
      if (.not.rtdb_put(rtdb,'dft:delta_ac',mt_dbl,1,shift))
     &   call errquit('dft_inp: rtdb_put failed', 2500, RTDB_ERR)
c
      goto 10
c
c     constrained dft
c
 2600 lcdft = .true.
      if (.not. rtdb_put(rtdb,'dft:lcdft', mt_log, 1, lcdft))
     &   call errquit('dft_input: rtdb_put failed', 2600, RTDB_ERR)
      nconstr = nconstr + 1
      call cdft_input(rtdb, nconstr, dtype, iatms1, iatms2, fatms1,
     $                fatms2, constr)
c
c Qin
      if (.not. rtdb_put(rtdb,'dft:nconstr', mt_int, 1, nconstr))
     &   call errquit('dft_input: rtdb_put failed', 2600, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'dft:dtype', mt_int, nconstr, dtype(1)))
     &   call errquit('dft_input: rtdb_put failed', 2600, RTDB_ERR)
      if (.not. (
     $          rtdb_put(rtdb,'dft:iatms1', mt_int, nconstr, iatms1(1))
     $    .and. rtdb_put(rtdb,'dft:iatms2', mt_int, nconstr, iatms2(1))
     $    .and. rtdb_put(rtdb,'dft:fatms1', mt_int, nconstr, fatms1(1))
     $    .and. rtdb_put(rtdb,'dft:fatms2', mt_int, nconstr, fatms2(1))
     $    ) )
     &   call errquit('dft_input: rtdb_put atms failed', 2600, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'dft:constr', mt_dbl, nconstr, constr(1)))
     &   call errquit('dft_input: rtdb_put failed', 2600, RTDB_ERR)
c Qin/end
c
      goto 10
c
c     == Coulomb attenuation method for Exchange ==
c
 2700 cam_exch = .true.
c      
c     persist the information in the runtime database
      if (.not.rtdb_put(rtdb,'dft:cam_exch', mt_log, 1, cam_exch)) 
     &   call errquit('dft_input: rtdb_put failed', 2700, RTDB_ERR)
      if (.not.inp_f(cam_omega)) cam_omega = 0.5d0
c           write(luout,*) "cam_omega:",cam_omega
      if (.not.rtdb_put(rtdb,'dft:cam_omega',mt_dbl,1,cam_omega))
     &   call errquit('dft_inp: rtdb_put failed cam_omega', 2700, 
     &     RTDB_ERR)
c
 2710 if (.not. inp_a(test)) goto 10
c  get cam_alpha
      if (inp_compare(.false.,'cam_alpha', test)) then
       if (.not. inp_f(cam_alpha)) cam_alpha = 1.d0
       if (.not.rtdb_put(rtdb,'dft:cam_alpha', mt_dbl, 1, cam_alpha))
     &  call errquit('dft_input: rtdb_put failed cam_alpha', 2710, 
     &     RTDB_ERR)
c           write(luout,*) "cam_alpha:",cam_alpha
      end if
c
 2720 if (.not. inp_a(test)) goto 10
c  get cam_beta
      if (inp_compare(.false.,'cam_beta', test)) then
       if (.not. inp_f(cam_beta)) cam_beta = 1.d0
       if (.not.rtdb_put(rtdb,'dft:cam_beta', mt_dbl, 1, cam_beta))
     &  call errquit('dft_input: rtdb_put failed cam_beta', 2720, 
     &     RTDB_ERR)
c           write(luout,*) "cam_beta:",cam_beta
      end if
c
c     Force texas integrals to false for CAM calculations
      if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &   call errquit('cannot set int:cando_txs',0, RTDB_ERR)
c
      goto 10
c
 2800 disp = .true.
c
      if (.not.rtdb_put(rtdb,'dft:disp', mt_log, 1, disp))
     &   call errquit('dft_input: rtdb_put failed', 2800, RTDB_ERR)

 2810 if (.not. inp_a(test)) goto 10
c version vdw
      if (inp_compare(.false.,'vdw', test)) then
          if (.not. inp_i(vdw))
     &     call errquit('dft_input: vdw case failed', 0, INPUT_ERR)
          if (.not.rtdb_put(rtdb,'dft:ivdw', mt_int, 1, vdw))
     &    call errquit('dft_input: rtdb_put failed vdw', 2810, RTDB_ERR)
c alpha damping parameter
      elseif (inp_compare(.false.,'alpha', test)) then
        if (.not. inp_f(alpha))
     &   call errquit('dft_input: alpha failed', 0, INPUT_ERR)
        if (.not.rtdb_put(rtdb,'dft:alpha',mt_dbl,1,alpha))
     &  call errquit('dft_input: rtdb_put failed alpha', 2810, RTDB_ERR)
c s6 contribution
      elseif (inp_compare(.false.,'s6', test)) then
        if (.not. inp_f(s6))
     &   call errquit('dft_input: s6 failed', 0, INPUT_ERR)
        if (.not.rtdb_put(rtdb,'dft:vdw',mt_dbl,1,s6))
     &   call errquit('dft_input: rtdb_put failed', 2810, RTDB_ERR)
c turn off vdw contrib
      elseif (inp_compare(.false.,'off', test)) then
        if (.not.rtdb_put(rtdb,'dft:disp',mt_dbl,1,.false.))
     &   call errquit('dft_input: rtdb_put failed', 2810, RTDB_ERR)
      else
c
        write(LuOut,2819)
c
 2819   format
     $  (/' valid disp structure input : '/,
     $    '  DISP  '/,
     $    '   VDW <real vdw integer default 2> '/,
     $    '   S6 <real s6 default depends on XC functional>  '/,
     $    '   ALPHA <real s6 default 20.0d0> '/ ,
     $    '   OFF '/)
c
        call errquit('dft_input: unknown directive,disp', 911,INPUT_ERR)
      endif
c
      goto 2810
c
c     dftmp2: enables fifth-rung functionals by modifying parameters
c               for the MP2 code so that it uses DFT reference orbitals, etc.
c               which_mp2 is necessary for task selection
c
 2900 dftmp2 = .true.
c
      if (.not.rtdb_put(rtdb,'dft:dftmp2', mt_log, 1, dftmp2))
     &   call errquit('dft_input: rtdb_put failed', 2901, RTDB_ERR)
c
        if (.not.inp_a(mp2char)) then
           write(LuOut,*)' dft_input: which_mp2 not specified,',
     1                   ' using semidirect MP2 by default.'
           which_mp2 = 1
        else
          if (mp2char.eq.'semidirect') then
            which_mp2 = 1
          else if (mp2char.eq.'semi') then
            which_mp2 = 1
          else if (mp2char.eq.'disk') then
            which_mp2 = 1
          else if (mp2char.eq.'direct') then
            which_mp2 = 2
          else if (mp2char.eq.'ri') then
            which_mp2 = 3
          else if (mp2char.eq.'rimp2') then
            which_mp2 = 3
          else
            write(LuOut,*)' dft_input: invalid which_mp2 input,',
     1                    ' using semidirect MP2 by default.'
            which_mp2 = 1
          endif
        endif
c
      if (.not. rtdb_put(rtdb, 'dft:which_mp2', mt_int, 1, which_mp2))
     &   call errquit('dft_input: rtdb_put failed', 2902, RTDB_ERR)
c
      goto 10
c
c     molden; 
c
 3000 molden = .true.
      if (.not. rtdb_put(rtdb, 'dft:molden', mt_log, 1, molden))
     &   call errquit('dft_input: rtdb_put failed', 3000, RTDB_ERR)
c
      goto 10
c
c     noscf; 
c
 3100 noscf = .true.
      if (.not. rtdb_put(rtdb, 'dft:noscf', mt_log, 1, noscf))
     &   call errquit('dft_input: rtdb_put failed', 3100, RTDB_ERR)
c
      goto 10
c
c     steric
c
 3200 steric = .true.
      if (.not.rtdb_put(rtdb,'dft:steric', mt_log, 1, steric))
     &   call errquit('dft_input: rtdb_put failed', 3200, RTDB_ERR)
c
      goto 10
c
c     efield
c
 3300 efield = .true.
      if (.not.rtdb_put(rtdb,'dft:efield', mt_log, 1, efield))
     &   call errquit('dft_input: rtdb_put failed', 3300, RTDB_ERR)
      if (.not. inp_f(efield_x)) 
     &   call errquit('dft_input: efield_x not found', 3300, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:efield_x', mt_dbl, 1, efield_x))
     &   call errquit('dft_input: rtdb_put failed', 3300, RTDB_ERR)
      if (.not. inp_f(efield_y)) 
     &   call errquit('dft_input: efield_y not found', 3300, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:efield_y', mt_dbl, 1, efield_y))
     &   call errquit('dft_input: rtdb_put failed', 3300, RTDB_ERR)
      if (.not. inp_f(efield_z)) 
     &   call errquit('dft_input: efield_z not found', 3300, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:efield_z', mt_dbl, 1, efield_z))
     &   call errquit('dft_input: rtdb_put failed', 3300, RTDB_ERR)
c
      goto 10
c
 9999 continue
c
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &   call errquit('dft_input: rtdb_put failed', 300, RTDB_ERR)
      return 
      write(6,10001)
10001 format(/' fon [alpha|beta] orbitals <norb> electrons <nelec>'/)
      call util_flush(6)
      call errquit('vectors_input: invalid format', 0, INPUT_ERR)
c
c     for semidirect
c
 3602 call errquit('semi-direct [filesize <integer filesize'//
     $     ' = disksize>] [memsize <integer memsize = available>]'//
     $     ' [filename <string filename>] ',450,INPUT_ERR)
c       
      end
