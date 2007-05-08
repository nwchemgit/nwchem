*
* $Id: nwpw_timing.F,v 1.6 2007-05-08 17:55:17 bylaska Exp $
*

* nwpw_timing.F
* Author - Eric Bylaska
*
*   These routines are to be used to time the nwpw module
*
*
*   1 - total FFT 
*   2 - total dot products
*   3 - lagrange multipliers
*   4 - exchange correlation
*   5 - local pseudopotentials
*   6 - non-local pseudopotentials
*   7 - hartree potentials
*   8 - structure factors
*   9 - masking and packing
*   10 - geodesic time
*   11 - gen psi_r and dn
*   12 - allocating memory from stack
*   13 - miscellaneous steepest descent update
*   15 - ffm_dgemm
*   16 - fmf_dgemm
*   17 - m_diagonalize

*
*   20 - phase factors
*   21 - ewald /ion-ion

*   22 - tredq
*   23 - getdiags
*   24 - tqliq
*   25 - eigsrt

*   30 - queue fft
*   31 - queue fft serial
*   32 - queue fft parallel
*   33 - HFX
*   40 - qmmm LJ
*   41 - qmmm residual Q


      subroutine nwpw_timing_init()
      implicit none

#include "nwpw_timing_common.fh"
      call dcopy(50,0.0d0,0,times,1)
      return 
      end

      subroutine nwpw_timing_start(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"

	  call current_second(nwpw_tim1(counter))
      return
      end

      subroutine nwpw_timing_end(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"

      
	  call current_second(nwpw_tim2(counter))
     
      times(counter) = times(counter) 
     >               + (nwpw_tim2(counter)-nwpw_tim1(counter))
      return
      end

      real*8 function nwpw_timing(counter)
      implicit none
      integer counter
#include "nwpw_timing_common.fh"

      nwpw_timing = times(counter)
      return
      end


      subroutine nwpw_timing_print_final(oprint,counter)
      implicit none
      logical  oprint
      integer counter

#include "stdio.fh"

c     **** external functions ****
      real*8   nwpw_timing
      external nwpw_timing

      if (oprint) then
      write(luout,1809) 'Time spent doing               ',
     >                  'total','step'
         write(luout,1808) '  FFTs                       : ',
     >                 nwpw_timing(1),
     >                 nwpw_timing(1)/dble(counter)
         write(luout,1808) '  dot products               : ',
     >                 nwpw_timing(2),
     >                 nwpw_timing(2)/dble(counter)
         write(luout,1808) '  geodesic                   : ',
     >                 nwpw_timing(10),
     >                 nwpw_timing(10)/dble(counter)
         write(luout,1808) '  ffm_dgemm                  : ',
     >                 nwpw_timing(15),
     >                 nwpw_timing(15)/dble(counter)
         write(luout,1808) '  fmf_dgemm                  : ',
     >                 nwpw_timing(16),
     >                 nwpw_timing(16)/dble(counter)
         write(luout,1808) '  m_diagonalize              : ',
     >                 nwpw_timing(17),
     >                 nwpw_timing(17)/dble(counter)

         write(luout,1808) '    - m_tredq                : ',
     >                 nwpw_timing(22),
     >                 nwpw_timing(22)/dble(counter)
         write(luout,1808) '    - m_getdiags             : ',
     >                 nwpw_timing(23),
     >                 nwpw_timing(23)/dble(counter)
         write(luout,1808) '    - m_tqliq                : ',
     >                 nwpw_timing(24),
     >                 nwpw_timing(24)/dble(counter)
         write(luout,1808) '    - m_eigsrt               : ',
     >                 nwpw_timing(25),
     >                 nwpw_timing(25)/dble(counter)

         write(luout,1808) '  exchange correlation       : ',
     >                 nwpw_timing(4),
     >                 nwpw_timing(4)/dble(counter)
         write(luout,1808) '  local pseudopotentials     : ',
     >                 nwpw_timing(5),
     >                 nwpw_timing(5)/dble(counter)
         write(luout,1808) '  non-local pseudopotentials : ',
     >                 nwpw_timing(6),
     >                 nwpw_timing(6)/dble(counter)
         write(luout,1808) '  hartree potentials         : ',
     >                 nwpw_timing(7),
     >                 nwpw_timing(7)/dble(counter)
         write(luout,1808) '  ion-ion interaction        : ',
     >                 nwpw_timing(21),
     >                 nwpw_timing(21)/dble(counter)
         write(luout,1808) '  structure factors          : ',
     >                 nwpw_timing(8),
     >                 nwpw_timing(8)/dble(counter)
         write(luout,1808) '  phase factors              : ',
     >                 nwpw_timing(20),
     >                 nwpw_timing(20)/dble(counter)
         write(luout,1808) '  masking and packing        : ',
     >                 nwpw_timing(9),
     >                 nwpw_timing(9)/dble(counter)
         write(luout,1808) '  queue fft                  : ',
     >                 nwpw_timing(30),
     >                 nwpw_timing(30)/dble(counter)
         write(luout,1808) '  queue fft (serial)         : ',
     >                 nwpw_timing(31),
     >                 nwpw_timing(31)/dble(counter)
         write(luout,1808) '  queue fft (message passing): ',
     >                 nwpw_timing(32),
     >                 nwpw_timing(32)/dble(counter)
         write(luout,1808) '  HFX potential              : ',
     >                 nwpw_timing(33),
     >                 nwpw_timing(33)/dble(counter)
         write(luout,1808) '  qmmm LJ                    : ',
     >                 nwpw_timing(40),
     >                 nwpw_timing(40)/dble(counter)
         write(luout,1808) '  qmmm residual Q            : ',
     >                 nwpw_timing(41),
     >                 nwpw_timing(41)/dble(counter)
      end if

      return 
 1808 FORMAT(A,E14.6,E14.6)
 1809 FORMAT(//A,2A14)
      end
