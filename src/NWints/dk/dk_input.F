      subroutine dk_input(rtdb)
*
* $Id$
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rel_consts.fh"
#include "util.fh"
#include "errquit.fh"
c     
c     Read input for calculations with Douglas-Kroll
c     
c     All input is of the form:
c     <keyword>                    or   
c     <keyword> <on||off>          or
c     <keyword> <fpp||dkh||dkfull> or a combination of the last two
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
      logical dk_run
      integer dk_type
c
c     Read input data beyond douglas-kroll; store in rtdb.
c     
      dk_run  = .true.
      dk_type = 2
c
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            dk_run = .false.
         elseif (inp_compare(.false., 'on', test)) then
            dk_run  = .true.
         elseif (inp_compare(.false.,'fpp',test)) then
            dk_type = 1
         elseif (inp_compare(.false.,'dkh',test)) then
            dk_type = 2
         elseif (inp_compare(.false.,'dkfull',test)) then
            dk_type = 3
         elseif (inp_compare(.false.,'dk1',test)) then
            dk_type = 1
         elseif (inp_compare(.false.,'dk1full',test)) then
            dk_type = 1
         elseif (inp_compare(.false.,'dk2',test)) then
            dk_type = 2
         elseif (inp_compare(.false.,'dk2full',test)) then
            dk_type = 3
         elseif (inp_compare(.false.,'dk3',test)) then
            dk_type = 4
         elseif (inp_compare(.false.,'dk3full',test)) then
            dk_type = 5
         elseif (inp_compare(.false.,'resc',test)) then
            dk_type = 0
         else
            call errquit('dk_input: unknown directive',0, UNKNOWN_ERR)
         endif
         goto 10
      endif
c
c     If we are going to do Douglas-Kroll, write out a line informing user.
c
      if ((ga_nodeid() .eq. 0) .and .dk_run .and. 
     $     util_print('douglas-kroll', print_low)) then
         write(luout,20)
  20     format(/'Douglas-Kroll integrals are switched on !')
         if (dk_type .eq. 1) then
            write(luout,21)
  21        format(/'Integrals for the no-pair equation with ',
     $             'free-particle projection operators'/
     $             'will be computed.'/)
         elseif (dk_type .eq. 2) then
            write(luout,22)
  22        format(/'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/
     $             'will be computed (Douglas-Kroll-Hess).'/)
         elseif (dk_type .eq. 3) then
            write(luout,23)
  23        format(/'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/,
     $             'including cross-product terms will be computed.'/)
         elseif (dk_type .eq. 4) then
            write(luout,24)
  24        format(/'Integrals for the no-pair equation with '/
     $             'Douglas-Kroll Hamiltonian ',
     $             'correct to third order in the external potential'/
     $             'will be computed (3rd-order Douglas-Kroll).'/)
         elseif (dk_type .eq. 5) then
            write(luout,25)
  25        format(/'Integrals for the no-pair equation with '/
     $             'Douglas-Kroll Hamiltonian ',
     $             'correct to third order in the external potential'/
     $             'including cross-product terms will be computed.'/)
         elseif (dk_type .eq. 0) then
            write(luout,200)
 200        format(/'Integrals for the relativistic scheme by ',
     $             'eliminating small-components of the Dirac equation'/
     $             'will be computed (RESC).'/)
         endif
      endif
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $   call errquit('dk_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'doug_kroll:type', mt_int, 1, dk_type))
     $   call errquit('dk_input: rtdb_put failed', 0, RTDB_ERR)
c
      return
c       
      end
