c
c $Id$
c

*     ***************************************
*     *					    *
*     *          psp_generator_input        *
*     *	  				    *
*     ***************************************
      subroutine psp_generator_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "btdb.fh"
c

      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 35)
      integer NMAX              ! maximum wavfunctions
      parameter (NMAX = 100)
 
      character*30 dirs(num_dirs)
      character*255 test,dir_name,sdir_name

      character*255 pspgenerator
      character*4  element
      double precision charge
      double precision mass_number
      character*25  pseudopotential_type
      character*25  solver_type
      character*25  correlation_type
      character*25  exchange_type,semicore_type
      character*255 pseudopotential_filename
      integer number_core,number_valence
      integer n_number(NMAX)
      character l_char
      integer l_number(NMAX)
      real*8 filling(NMAX)
      real*8 rcut(10)
      integer max_rcut,l_rcut(10)
      real*8 rlocal,clocal
      integer vstates
      integer v_n(NMAX)
      integer v_l(NMAX)
      real*8  v_eig(NMAX)
      real*8  v_rcut(NMAX)
      real*8  r_semicore
      integer lmax00,locp00,efg00
      real*8  rlocal00,screening_cut

      character lmax,locp
      integer h,i,j,k,debug,print
      data dirs / 'element:',
     >            'charge:',
     >            'mass_number:',
     >            'pseudopotential_type:',
     >            'atomic_filling:',
     >            'cutoff:',
     >            'maximum_angular_momentum:',
     >            'pseudopotential_filename:',
     >            'solver_type:',
     >            'correlation_type:',
     >            'exchange_type:',
     >            'vanderbilt_local:',
     >            'vanderbilt_states:',
     >            'semicore_radius:',
     >            'element',
     >            'charge',
     >            'mass_number',
     >            'pseudopotential_type',
     >            'atomic_filling',
     >            'cutoff',
     >            'maximum_angular_momentum',
     >            'pseudopotential_filename',
     >            'solver_type',
     >            'correlation_type',
     >            'exchange_type',
     >            'vanderbilt_local',
     >            'vanderbilt_states',
     >            'semicore_radius',
     >            'lmax',
     >            'locp',
     >            'aperiodic_cutoff',
     >            'screening_cutoff',
     >            'semicore_type',
     >            'efg_type',
     >            'end'/



*     ***** initializations ****
      correlation_type     = 'vosko'
      exchange_type        = 'dirac'
      solver_type          = 'pauli'
      pseudopotential_type = 'hamann'
      element = 'H'
      charge = 1.0d0
      mass_number = 1.0d0
      number_core =0
      number_valence =1
      n_number(1) = 1
      l_number(1) = 0
      filling(1) = 1.0d0
      rlocal     = -1.0d0
      clocal     = 0.0d0
      vstates    = -1
      max_rcut = 0
      lmax = '?'
      locp = '?'
      r_semicore = -1.0d0
      semicore_type = 'fuchs'
      rlocal00  = 1.0d0
      screening_cut = 0.0d0
      efg00 = 0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'psp_generator_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'psp_generator_input: unknown directive', 0, INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 600, 700,800,900,
     >       1000,1100,1200,1300,1400,
     >       100,  200, 300, 400, 500, 600, 700,800,900,
     >       1000,1100,1200,1300,1400,
     >       700,710,1500,1600,1450,1700,
     >      9999) ind
      call errquit(
     >      'psp_generator_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  element
c     
 100  if (.not. inp_a(element))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c charge
c
 200  if (.not. inp_f(charge))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c  mass_number
c 
 300  if (.not. inp_f(mass_number))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)

      goto 10

c
c  pseudopotential_type
c 
 400  if (.not. inp_a(pseudopotential_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10

c
c  atomic_filling
c 
 500  if (.not. inp_i(number_core))
     >     call errquit(
     >          'psp_generator_input: failed to read ncore', 0,
     &       INPUT_ERR)

      if (.not. inp_i(number_valence))
     >     call errquit(
     >          'psp_generator_input: failed to read nvalence', 0,
     &       INPUT_ERR)

      do i=1,(number_core+number_valence)
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(n_number(i)))
     >     call errquit(
     >          'psp_generator_input: inp_i(n_number) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0,
     &       INPUT_ERR)

           if (l_char.eq.'s') l_number(i) = 0
           if (l_char.eq.'p') l_number(i) = 1
           if (l_char.eq.'d') l_number(i) = 2
           if (l_char.eq.'f') l_number(i) = 3
           if (l_char.eq.'g') l_number(i) = 4
           if (l_char.eq.'S') l_number(i) = 0
           if (l_char.eq.'P') l_number(i) = 1
           if (l_char.eq.'D') l_number(i) = 2
           if (l_char.eq.'F') l_number(i) = 3
           if (l_char.eq.'G') l_number(i) = 4
           if (l_char.eq.'0') l_number(i) = 0
           if (l_char.eq.'1') l_number(i) = 1
           if (l_char.eq.'2') l_number(i) = 2
           if (l_char.eq.'3') l_number(i) = 3
           if (l_char.eq.'4') l_number(i) = 4

           if (.not. inp_f(filling(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(filling) failed', 0,
     &       INPUT_ERR)

      end do
      goto 10
 
c
c  cutoff 
c 
 600  if (.not. inp_i(max_rcut))
     >     call errquit(
     >          'psp_generator_input: failed to read max_rcut', 0,
     &       INPUT_ERR)

      do i=1,(max_rcut)
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0, INPUT_ERR)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0,
     &       INPUT_ERR)

           if (l_char.eq.'s') l_rcut(i) = 0
           if (l_char.eq.'p') l_rcut(i) = 1
           if (l_char.eq.'d') l_rcut(i) = 2
           if (l_char.eq.'f') l_rcut(i) = 3
           if (l_char.eq.'g') l_rcut(i) = 4
           if (l_char.eq.'S') l_rcut(i) = 0
           if (l_char.eq.'P') l_rcut(i) = 1
           if (l_char.eq.'D') l_rcut(i) = 2
           if (l_char.eq.'F') l_rcut(i) = 3
           if (l_char.eq.'G') l_rcut(i) = 4
           if (l_char.eq.'0') l_rcut(i) = 0
           if (l_char.eq.'1') l_rcut(i) = 1
           if (l_char.eq.'2') l_rcut(i) = 2
           if (l_char.eq.'3') l_rcut(i) = 3
           if (l_char.eq.'4') l_rcut(i) = 4
           if (.not. inp_f(rcut(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(rcut) failed', 0,
     &       INPUT_ERR)
      end do

      goto 10
 
 
c
c  maximum_angular_momentum - lmax
c 
 700  if (.not. inp_a(lmax))
     >     call errquit(
     >          'psp_generator_input: failed to read lmax', 0,
     &       INPUT_ERR)
      goto 10

c
c  local potential - locp
c 
 710  if (.not. inp_a(locp))
     >     call errquit(
     >          'psp_generator_input: failed to read locp',0,0)
      goto 10

 
c
c  pseudopotential_filename
c 
 800  if (.not. inp_a(pseudopotential_filename))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10

c
c  solver_type
c 
 900  if (.not. inp_a(solver_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10

c
c  correlation_type
c 
 1000 if (.not. inp_a(correlation_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10

c
c  exchange_type
c 
 1100 if (.not. inp_a(exchange_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c  vanderbilt_local
c 
 1200 if (.not. inp_f(rlocal))
     >     call errquit(
     >          'psp_generator_input: failed to read rlocal', 0,
     &       INPUT_ERR)
      if (.not. inp_f(clocal))
     >     call errquit(
     >          'psp_generator_input: failed to read clocal', 0,
     &       INPUT_ERR)
      goto 10


c
c  vanderbilt_states
c
 1300 if (.not. inp_i(vstates))
     >  call errquit(
     >   'psp_generator_input: failed to read vstates', 0,
     &       INPUT_ERR)
  
      do i=1,vstates
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0,
     &       INPUT_ERR)

           if (.not. inp_i(v_n(i)))
     >     call errquit(
     >          'psp_generator_input: inp_i(v_n) failed', 0,
     &       INPUT_ERR)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0,
     &       INPUT_ERR)

           if (l_char.eq.'s') v_l(i) = 0
           if (l_char.eq.'p') v_l(i) = 1
           if (l_char.eq.'d') v_l(i) = 2
           if (l_char.eq.'f') v_l(i) = 3
           if (l_char.eq.'g') v_l(i) = 4
           if (l_char.eq.'S') v_l(i) = 0
           if (l_char.eq.'P') v_l(i) = 1
           if (l_char.eq.'D') v_l(i) = 2
           if (l_char.eq.'F') v_l(i) = 3
           if (l_char.eq.'G') v_l(i) = 4
           if (l_char.eq.'0') v_l(i) = 0
           if (l_char.eq.'1') v_l(i) = 1
           if (l_char.eq.'2') v_l(i) = 2
           if (l_char.eq.'3') v_l(i) = 3
           if (l_char.eq.'4') v_l(i) = 4

           if (.not. inp_f(v_eig(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(v_eig) failed', 0,
     &       INPUT_ERR)
         if (.not. inp_f(v_rcut(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(v_rcut) failed', 0,
     &       INPUT_ERR)
         
      end do

      goto 10
c
c  semicore_radius
c 
 1400 if (.not. inp_f(r_semicore))
     >     call errquit(
     >     'psp_generator_input: failed to read semicore_radius',0,
     &       INPUT_ERR)
      goto 10

c
c  semicore_type
c 
 1450 if (.not. inp_a(semicore_type))
     >     call errquit(
     >     'psp_generator_input: failed to read semicore_type',0,
     &       INPUT_ERR)
      goto 10
 
 
c
c  aperiodic cutoff radius
c 
 1500 if (.not. inp_f(rlocal00)) rlocal00 = 1.0d0
      goto 10
 
c
c  screening cutoff radius
c
 1600 if (.not. inp_f(screening_cut)) screening_cut = 0.0d0
      goto 10

c
c  efg_type
c
 1700 if (.not. inp_i(efg00)) efg00 = 1
      goto 10


*     ***** generate temporary data file  and run pspsolve ****
 9999 continue
      k = index(pseudopotential_filename,' ') -1
      pspgenerator = pseudopotential_filename(1:k)//'.generator'
      call util_file_name_noprefix(pspgenerator, 
     >                             .true.,
     >                             .false.,
     >                             test)
      h = index(test,' ')-1
      open(unit=17,file=test(1:h),form='formatted')
      write(17,*) "<atom>"
      write(17,*) element
      write(17,111) charge,mass_number,number_core,number_valence
 111  format(2E14.6,2I4)
      do i=1,(number_core+number_valence)
         write(17,112) n_number(i),l_number(i),filling(i)
      end do
 112  format(2I4,E14.6)
      write(17,*) "<end>"
      write(17,*) "<solver> ",solver_type," <end>"
      write(17,*) "<correlation>  ",correlation_type,"  <end>"
      write(17,*) "<exchange>     ",exchange_type,   "  <end>"
      write(17,*) "<pseudopotential> ",pseudopotential_type," <end>"
      if (r_semicore.gt. 0.0d0) then
        write(17,*) "<semicore> ",r_semicore," <end>" 
        write(17,*) "<semicore_type> ",semicore_type," <end>" 
      end if
      if (lmax.ne.'?') write(17,*) "<lmax> ",lmax," <end>"
      if (max_rcut.gt.0) then
         write(17,*) "<rcut>"
         do i=1,max_rcut
            write(17,113) l_rcut(i),rcut(i)
         end do
         write(17,*) "<end>"
      end if
 113  format(I4,E14.6)

      if (rlocal.gt.0.0) then
         write(17,*) "<vanderbilt-local>", rlocal,clocal," <end>"
      end if
      if (vstates.gt.0) then
         write(17,*) "<npsp-states>",vstates," <end>"
         write(17,*) "<vanderbilt-states>"
         do i=1,vstates
            write(17,114) (i-1),v_n(i),v_l(i),v_eig(i),v_rcut(i)
         end do
         write(17,*) "<end>"
      end if
 114  format(3I4,2E14.6)

      if (screening_cut.gt.0.0)
     >  write(17,*) "<screening_cut> ",screening_cut," <end>" 

      close(17)

      write(*,*) 
     >   ">>>> PSPW Serial Module - 1d pseudopotential solver <<<<"

      call util_directory_name(sdir_name,.true.,0)
      call util_directory_name(dir_name,.false.,0)
      h = index(sdir_name,' ') -1
      i = index(dir_name,' ') -1
      j = index(pspgenerator,' ')-1
      k = index(pseudopotential_filename,' ') -1

      debug = 1
      print = 1
      lmax00   = -1
      if (lmax.eq.'s') lmax00 = 0
      if (lmax.eq.'p') lmax00 = 1
      if (lmax.eq.'d') lmax00 = 2
      if (lmax.eq.'f') lmax00 = 3
      if (lmax.eq.'g') lmax00 = 4
      if (lmax.eq.'S') lmax00 = 0
      if (lmax.eq.'P') lmax00 = 1
      if (lmax.eq.'D') lmax00 = 2
      if (lmax.eq.'F') lmax00 = 3
      if (lmax.eq.'G') lmax00 = 4
      if (lmax.eq.'0') lmax00 = 0
      if (lmax.eq.'1') lmax00 = 1
      if (lmax.eq.'2') lmax00 = 2
      if (lmax.eq.'3') lmax00 = 3
      if (lmax.eq.'4') lmax00 = 4

      locp00   = -1
      if (locp.eq.'s') locp00 = 0
      if (locp.eq.'p') locp00 = 1
      if (locp.eq.'d') locp00 = 2
      if (locp.eq.'f') locp00 = 3
      if (locp.eq.'g') locp00 = 4
      if (locp.eq.'S') locp00 = 0
      if (locp.eq.'P') locp00 = 1
      if (locp.eq.'D') locp00 = 2
      if (locp.eq.'F') locp00 = 3
      if (locp.eq.'G') locp00 = 4
      if (locp.eq.'0') locp00 = 0
      if (locp.eq.'1') locp00 = 1
      if (locp.eq.'2') locp00 = 2
      if (locp.eq.'3') locp00 = 3
      if (locp.eq.'4') locp00 = 4

      !rlocal00 = 1.0d0 ! 
      call pspsolve(print,debug,
     >              lmax00,locp00,rlocal00,efg00,
     >              sdir_name(1:h),   h,
     >              dir_name(1:i),    i,
     >              pspgenerator(1:j),j,
     >              pseudopotential_filename(1:k),k)
      call nwpw_message(4) 

c      write(6,*)"done with psp gen pspsolve!"
c      call flush(6)

      return
      end

*     ***************************************
*     *					    *
*     *          psp_generator_auto         *
*     *					    *
*     ***************************************
      subroutine psp_generator_auto(iopt,atom,lmax,locp,rlocal)
      implicit none

#include "inp.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"
#include "stdio.fh"

      integer iopt
      character*4 atom
      double precision rlocal
      integer	    locp,lmax
c
      integer MASTER,taskid,msglen
      parameter (MASTER=0)

      logical value,hprint
      character*255 dir_name,sdir_name
      character*16 psp_in,psp_out
      integer h,i,j,k,ka,debug,psp_type,print,efg_type
      character ch_tmp

*     **** external functions ****
      logical  control_print
      integer  nwpw_libglmax,nwpw_libglocp,control_gga
      integer  nwpw_libgpsp_type
      integer  nwpw_libgefg_type
      double precision nwpw_libgrlocal
      external control_print
      external nwpw_libglmax,nwpw_libglocp,control_gga
      external nwpw_libgpsp_type
      external nwpw_libgefg_type
      external nwpw_libgrlocal
c
      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)

      value = btdb_parallel(.false.)
      if (taskid.eq.MASTER) then

      call util_directory_name(sdir_name,.true.,0)
      call util_directory_name(dir_name,.false.,0)
      h = index(sdir_name,' ') - 1
      i = index(dir_name,' ') - 1
c      
c     clean crap
c
      open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='unknown')
      close(unit=99,status='delete')
      psp_in = "junk.inp "
      ka = index(atom,' ') - 1
      if (ka.lt.0) ka = 4
      if (iopt.eq.1) psp_out = atom(1:ka)//".psp "
      if (iopt.eq.2) psp_out = atom(1:ka)//"_basis "
      j = index(psp_in,' ')  - 1
      k = index(psp_out,' ') - 1

      call nwpw_libgeninp(iopt,atom,sdir_name(1:h)//'/junk.inp')
      lmax=nwpw_libglmax(atom)
      locp=nwpw_libglocp(atom)
      rlocal=nwpw_libgrlocal(atom)
      psp_type=nwpw_libgpsp_type(atom)
      efg_type=nwpw_libgefg_type(atom)
       
c      if (hprint) then
      write(luout,*) "Generating 1d pseudopotential for ",atom

c      end if

      if ((psp_type.eq.0).or.(psp_type.eq.7)) then


c        *** add GGA options to junk.inp ***
*  
c        *** Vosko ***
         if (control_gga().eq.0) then
           if (hprint) write(*,*) "+ generated using Vosko et al."

c        *** PBE96 ***
         else if ((control_gga().eq.10).or.(control_gga().eq.110)) then
           if (hprint) write(*,*) "+ generated using PBE96"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=30,END=30) ch_tmp
           end do
   30      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    pbe96 <end>'
           write(99,*) '<correlation> pbe96 <end>'
           close(unit=99)

c        *** BLYP - use b3lyp gga until fixed! ***
         else if ((control_gga().eq.11).or.(control_gga().eq.115)) then
           !if (hprint) write(*,*) "+ generated using PBE96 (BLYP psp)"
           if (hprint) write(*,*) "+ generated using BLYP psp"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=31,END=31) ch_tmp
           end do
   31      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    becke <end>'
           write(99,*) '<correlation> lyp   <end>'
           !write(99,*) '<exchange>    pbe96 <end>'
           !write(99,*) '<correlation> pbe96 <end>'
           close(unit=99)

c        *** revPBE ***
         else if ((control_gga().eq.12).or.(control_gga().eq.112)) then
           if (hprint) write(*,*) "+ generated using revPBE"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=32,END=32) ch_tmp
           end do
   32      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    revpbe <end>'
           write(99,*) '<correlation> revpbe <end>'
           close(unit=99)
         end if

         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call pspsolve(print,debug,
     >        lmax,locp,rlocal,efg_type,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k)

      !*** HGH pseudopotential format ***
      else if (psp_type.eq.1) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call hgh_parse(debug,    
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k)

      !*** CPI pseudopotential format ***
      else if (psp_type.eq.2) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call cpi_parse(debug,    
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** TETER pseudopotential format ***
      else if (psp_type.eq.3) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call teter_parse(debug,    
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** PAW basis format ***
      else if (psp_type.eq.4) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call paw_atom_driver(debug,
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** QMMM format ***
      else if (psp_type.eq.5) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call qmmm_parse(debug,    
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** CARTER pseudopotential format ***
      else if (psp_type.eq.6) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call carter_parse(debug,    
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)
      end if

      open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
      close(unit=99,status='delete')


      end if !** taskid.eq.MASTER **
      value = btdb_parallel(.true.)

      call ga_sync()

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp)
      call Parallel_Brdcst_values(MASTER,msglen,rlocal)

      return
      end



*     ***************************************
*     *					    *
*     *          psp_generator_auto1        *
*     *					    *
*     ***************************************
      subroutine psp_generator_auto1(atom,lmax,locp,rlocal)
      implicit none

#include "inp.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"
#include "stdio.fh"

      character*4 atom
      double precision rlocal
      integer	    locp,lmax
c
      integer MASTER,taskid,msglen
      parameter (MASTER=0)

      logical value,hprint
      character*255 dir_name,sdir_name
      character*16 psp_in,psp_out
      integer h,i,j,k,ka,debug,psp_type,print
      character ch_tmp

*     **** external functions ****
      logical  control_print
      integer  nwpw_libglmax,nwpw_libglocp,control_gga
      integer  nwpw_libgpsp_type
      double precision nwpw_libgrlocal
      external control_print
      external nwpw_libglmax,nwpw_libglocp,control_gga
      external nwpw_libgpsp_type
      external nwpw_libgrlocal
c
      call Parallel_taskid(taskid)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)

      value = btdb_parallel(.false.)
      if (taskid.eq.MASTER) then

      call util_directory_name(sdir_name,.true.,0)
      call util_directory_name(dir_name,.false.,0)
      h = index(sdir_name,' ') - 1
      i = index(dir_name,' ') - 1
c      
c     clean crap
c
      open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='unknown')
      close(unit=99,status='delete')
      psp_in = "junk.inp "
      ka = index(atom,' ') - 1
      if (ka.lt.0) ka = 4
      psp_out = atom(1:ka)//".psp1"
      j = index(psp_in,' ')  - 1
      k = index(psp_out,' ') - 1

      call nwpw_libgeninp(1,atom,sdir_name(1:h)//'/junk.inp')
      lmax=nwpw_libglmax(atom)
      locp=nwpw_libglocp(atom)
      rlocal=nwpw_libgrlocal(atom)
      psp_type=nwpw_libgpsp_type(atom)
       
      if (hprint) then
      write(luout,*) "Generating 1d atomic orbitals for ",atom
      end if

      !*** CPI pseudopotential format ***
      if (psp_type.eq.2) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call cpi_parse(debug,
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** TETER pseudopotential format ***
      else if (psp_type.eq.3) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call teter_parse(debug,
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      !*** PAW pseudopotential format ***
      else if (psp_type.eq.4) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call paw_atom_driver(debug,
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      else if (psp_type.eq.6) then
         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call carter_parse(debug,
     >        lmax,locp,rlocal,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k,
     >        atom(1:ka),ka)

      else

c        *** add GGA options to junk.inp ***
*  
c        *** Vosko ***
         if (control_gga().eq.0) then
           if (hprint) write(*,*) "+ generated using Vosko et al."

c        *** PBE96 ***
         else if ((control_gga().eq.10).or.(control_gga().eq.110)) then
           if (hprint) write(*,*) "+ generated using PBE96"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=30,END=30) ch_tmp
           end do
   30      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    pbe96 <end>'
           write(99,*) '<correlation> pbe96 <end>'
           close(unit=99)

c        *** BLYP - use until blyp vxc fixed ***
         else if (control_gga().eq.11) then
           !if (hprint) write(*,*) "+ generated using PBE96 (BLYP psp)"
           if (hprint) write(*,*) "+ generated using BLYP"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=31,END=31) ch_tmp
           end do
   31      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    becke <end>'
           write(99,*) '<correlation> lyp   <end>'
           !write(99,*) '<exchange>    pbe96 <end>'
           !write(99,*) '<correlation> pbe96 <end>'
           close(unit=99)

c        *** revPBE ***
         else if ((control_gga().eq.12).or.(control_gga().eq.112)) then
           if (hprint) write(*,*) "+ generated using revPBE"
           open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
           do while(.true.)
             read(99,*,ERR=32,END=32) ch_tmp
           end do
   32      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(MACX) || defined(GCC46)
           backspace 99
#endif
           write(99,*) '<exchange>    revpbe <end>'
           write(99,*) '<correlation> revpbe <end>'
c           write(99,*) '<correlation> off   <end>'
           close(unit=99)


         end if


         debug = 0
         print = 0
         if (control_print(print_debug)) debug = 1
         if (control_print(print_high)) print  = 1
         call pspsolve(print,debug,
     >        lmax,locp,rlocal,0,
     >        sdir_name(1:h),h,
     >        dir_name(1:i),i,
     >        psp_in(1:j),j,
     >        psp_out(1:k),k)


      end if
      open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
      close(unit=99,status='delete')


      end if !** taskid.eq.MASTER **
      value = btdb_parallel(.true.)

      call ga_sync()

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,lmax)
      call Parallel_Brdcst_ivalues(MASTER,msglen,locp)
      call Parallel_Brdcst_values(MASTER,msglen,rlocal)

      return
      end
