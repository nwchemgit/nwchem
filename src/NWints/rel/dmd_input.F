      subroutine dmd_input(rtdb)
*
* $Id: dmd_input.F,v 1.5 2000-01-25 19:37:20 mg201 Exp $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "rel_consts.fh"
c     
c     Read input for calculations with Dyall-Modified-Dirac or Dyall-Mod-Dirac
c     
c     All input is of the form:
c     <keyword>                          or   
c     <keyword> <on||off>                or
c     <keyword> <sf>  !future option     or
c     <keyword> <afw> !future option     or a combination of the last three
c     
      integer rtdb ! [input] handle to database
c    
      character*255 test
      logical dmd_run, atomic_nesc, one_e_approx, one_cent
c
c     Read input data beyond Dyall-Modified-Dirac or Dyall-Mod-Dirac; store in rtdb.
c     
      dmd_run  = .true.
      atomic_nesc = .true.
      one_e_approx = .true.
c
  10  if (inp_a(test)) then
        if (inp_compare(.false., 'off', test)) then
          dmd_run = .false.
          if (.not. rtdb_put(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run))
     &        call errquit('dmd_input: rtdb_put failed', 0)
          return
        elseif (inp_compare(.false., 'on', test)) then
          dmd_run  = .true.
        else if (inp_compare(.false., 'nesc1e', test)) then
          one_e_approx = .true.
        else if (inp_compare(.false., 'nesc2e', test)) then
          one_e_approx = .false.
          one_cent = .false.
          if (inp_a(test)) then
            if (inp_compare(.false., 'ss1cent', test)) then
              one_cent = .true.
            else
              call errquit('dmd_input: unknown option for nesc2e',0)
            end if
          end if
        else
          call errquit('dmd_input: unknown directive',0)
        end if
        goto 10
      endif
c
c     If we are going to do Dyall-Modified-Dirac then print out information 
c
      if ((ga_nodeid() .eq. 0) .and . dmd_run .and. 
     $    util_print('dyall-mod-dirac', print_low)) then
        write(luout,20)
   20   format(//' Relativistic calculation using Dyall''s modification'
     &      ,' to the Dirac equation,',
     &      /5x,'K. G. Dyall, J. Chem. Phys. 100, 2118 (1994) ',/)
        if (atomic_nesc) write (luout,30)
   30   format(' Calculation performed with normalized elimination of ',
     &      'small components (NESC) and atomic 4-spinor expansion,',/,
     &      5x,'K. G. Dyall, J. Chem. Phys. 106, 9618 (1997); K. G. Dy',
     &      'all and T. Enevoldsen, J. Chem. Phys. 111, 10000 (1999).')
        if (one_e_approx) then
          write (luout,40)
   40     format (' Electron-electron interaction approximated by ',
     &        'neglecting the commutator of 1/r12 with the atomic FW ',
     &        'transformation,',/5x,
     &        'K. G. Dyall, in preparation.')
        else if (one_cent) then
          write (luout,50)
   50     format (' Only one-center small component densities used',
     &        'in two-electron integrals.')
        end if
        write (luout,*)
      endif
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb,'dyall_mod_dir',mt_log,1,dmd_run))
     &    call errquit('dmd_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
     &    call errquit('dmd_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb,'nesc_1e_approx',mt_log,1,one_e_approx))
     &    call errquit('dmd_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb,'ss_one_cent',mt_log,1,one_cent)) 
     &    call errquit('dmd_input: rtdb_put failed', 0)
c
      return
c       
      end
