*     ***************************************
*     *										*
*     *          psp_generator_input    	*
*     *										*
*     ***************************************
      subroutine psp_generator_input(rtdb)
c
c $Id: psp_generator_input.F,v 1.12 2002-02-13 20:02:06 edo Exp $
c
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c

      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 15)
      integer NMAX              ! maximum wavfunctions
      parameter (NMAX = 100)
 
      character*30 dirs(num_dirs)
      character*255 test,dir_name,sdir_name

      character*255 pspgenerator
      character*2  element
      double precision charge
      double precision mass_number
      character*25  pseudopotential_type
      character*25  solver_type
      character*25  correlation_type
      character*25  exchange_type
      character*25  pseudopotential_filename
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
      integer lmax00,locp00
      real*8  rlocal00

      character lmax
      integer h,i,j,k,debug
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
      r_semicore = -1.0d0

 10   if (.not. inp_read()) 
     >     call errquit(
     >           'psp_generator_input: inp_read failed', 0)
      if (.not. inp_a(test))
     >     call errquit(
     >           'psp_generator_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'psp_generator_input: unknown directive', 0)


      goto ( 100,  200, 300, 400, 500, 600, 700,800,900,
     >       1000,1100,1200,1300,1400,
     >      9999) ind
      call errquit(
     >      'psp_generator_input: unimplemented directive', ind)


c     
c  element
c     
 100  if (.not. inp_a(element))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10


c
c charge
c
 200  if (.not. inp_f(charge))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10


c
c  mass_number
c 
 300  if (.not. inp_f(mass_number))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)

      goto 10

c
c  pseudopotential_type
c 
 400  if (.not. inp_a(pseudopotential_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10

c
c  atomic_filling
c 
 500  if (.not. inp_i(number_core))
     >     call errquit(
     >          'psp_generator_input: failed to read ncore', 0)

      if (.not. inp_i(number_valence))
     >     call errquit(
     >          'psp_generator_input: failed to read nvalence', 0)

      do i=1,(number_core+number_valence)
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0)

           if (.not. inp_i(n_number(i)))
     >     call errquit(
     >          'psp_generator_input: inp_i(n_number) failed', 0)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0)

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

           if (.not. inp_f(filling(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(filling) failed', 0)

      end do
      goto 10
 
c
c  cutoff 
c 
 600  if (.not. inp_i(max_rcut))
     >     call errquit(
     >          'psp_generator_input: failed to read max_rcut', 0)

      do i=1,(max_rcut)
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0)

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
           if (.not. inp_f(rcut(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(rcut) failed', 0)
      end do

      goto 10
 
 
c
c  maximum_angular_momentum
c 
 700  if (.not. inp_a(lmax))
     >     call errquit(
     >          'psp_generator_input: failed to read lmax', 0)
      goto 10
 
c
c  pseudopotential_filename
c 
 800  if (.not. inp_a(pseudopotential_filename))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10

c
c  solver_type
c 
 900  if (.not. inp_a(solver_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10

c
c  correlation_type
c 
 1000 if (.not. inp_a(correlation_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10

c
c  exchange_type
c 
 1100 if (.not. inp_a(exchange_type))
     >     call errquit(
     >          'psp_generator_input: failed to read keyword', 0)
      goto 10


c
c  vanderbilt_local
c 
 1200 if (.not. inp_f(rlocal))
     >     call errquit(
     >          'psp_generator_input: failed to read rlocal', 0)
      if (.not. inp_f(clocal))
     >     call errquit(
     >          'psp_generator_input: failed to read clocal', 0)       
      goto 10


c
c  vanderbilt_states
c
 1300 if (.not. inp_i(vstates))
     >  call errquit(
     >   'psp_generator_input: failed to read vstates', 0)
  
      do i=1,vstates
         if (.not. inp_read()) 
     >        call errquit(
     >           'psp_generator_input: inp_read failed', 0)

           if (.not. inp_i(v_n(i)))
     >     call errquit(
     >          'psp_generator_input: inp_i(v_n) failed', 0)

           if (.not. inp_a(l_char))
     >     call errquit(
     >          'psp_generator_input: inp_i(l_char) failed', 0)

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

           if (.not. inp_f(v_eig(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(v_eig) failed', 0)
         if (.not. inp_f(v_rcut(i)))
     >     call errquit(
     >          'psp_generator_input: inp_f(v_rcut) failed', 0)
         
      end do

      goto 10
c
c  semicore_radius
c 
 1400 if (.not. inp_f(r_semicore))
     >     call errquit(
     >     'psp_generator_input: failed to read semicore_radius',0)
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
      if (r_semicore.gt. 0.0d0) 
     >  write(17,*) "<semicore> ",r_semicore," <end>" 
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
      locp00   = -1
      rlocal00 = 1.0d0
      call pspsolve(debug,
     >              lmax00,locp00,rlocal00,
     >              sdir_name(1:h),   h,
     >              dir_name(1:i),    i,
     >              pspgenerator(1:j),j,
     >              pseudopotential_filename(1:k),k)
      call nwpw_message(4) 

      return
      end

*     ***************************************
*     *										*
*     *          psp_generator_auto     	*
*     *										*
*     ***************************************
      subroutine psp_generator_auto(atom,lmax,locp,rlocal)
      implicit none

#include "inp.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

      character*2 atom
      double precision rlocal
      integer	    locp,lmax
c
      integer MASTER,taskid,msglen
      parameter (MASTER=0)

      character*255 dir_name,sdir_name
      character*16 psp_in,psp_out
      integer h,i,j,k,debug
      integer nwpw_libglmax,nwpw_libglocp
      double precision nwpw_libgrlocal
      external nwpw_libglmax,nwpw_libglocp,nwpw_libgrlocal
c
      call Parallel_taskid(taskid)

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
      k = index(atom,' ') - 1
      if (k.lt.0) k = 2
      psp_out = atom(1:k)//".psp "
      j = index(psp_in,' ')  - 1
      k = index(psp_out,' ') - 1

      call nwpw_libgeninp(atom,sdir_name(1:h)//'/junk.inp')
      lmax=nwpw_libglmax(atom)
      locp=nwpw_libglocp(atom)
      rlocal=nwpw_libgrlocal(atom)
      write(*,*) "Generating 1d pseudopotential for ",atom
      debug = 0
      call pspsolve(debug,
     >     lmax,locp,rlocal,
     >     sdir_name(1:h),h,
     >     dir_name(1:i),i,
     >     psp_in(1:j),j,
     >     psp_out(1:k),k)
!      stop
      open(unit=99,file=sdir_name(1:h)//'/junk.inp',status='old')
      close(unit=99,status='delete')

      end if
      call ga_sync()

      msglen = 1
      call BRDCST(9+MSGINT,lmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,locp,mitob(msglen),MASTER)
      call BRDCST(9+MSGDBL,rlocal,mitob(msglen),MASTER)

      return
      end
