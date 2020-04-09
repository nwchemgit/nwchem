*
* $Id$
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
*   18 - mmm_dgemm

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

*   34 - paw gaussian integrals
*   35 - paw atomic coulomb
*   36 - paw atomic xc
*   37 - paw gen dEmult/dQlm
*   38 - paw gen dElocal/dQlm
*   39 - paw cmp operations

*   40 - qmmm LJ
*   41 - qmmm residual Q

*   42 - MATHIAS InnerLoop
*   43 - MATHIAS Phaze
*   44 - MATHIAS Pipelined FFTs
*   45 - MATHIAS Lagrange
*   46 - MATHIAS Exch Corr
*   47 - MATHIAS Hpsi


*   50 - io time

*   52 - HFX localization
*   53 - HFX DM columns
*   54 - HFX DM Cholesky
*   55 - re-gridding


      subroutine nwpw_timing_init()
      implicit none

#include "nwpw_timing_common.fh"
      call current_second(t0)
      call dcopy(nwpw_tim_max,0.0d0,0,times,1)
      call dcopy(272*nwpw_tim_max,0.0d0,0,thr_times,1)
      return 
      end

      subroutine nwpw_timing_start(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"
!$OMP MASTER
      call current_second(nwpw_tim1(counter))
!$OMP END MASTER
      return
      end

      subroutine nwpw_timing_end(counter)
      implicit none
      integer counter

#include "nwpw_timing_common.fh"

      
!$OMP MASTER
      call current_second(nwpw_tim2(counter))
     
      times(counter) = times(counter) 
     >               + (nwpw_tim2(counter)-nwpw_tim1(counter))
!$OMP END MASTER
      return
      end

      subroutine nwpw_timing_start_thr(counter)
      USE omp_lib
      implicit none
      integer counter
      integer tid 
#include "nwpw_timing_common.fh"
#ifdef USE_OPENMP
      tid = omp_get_thread_num()
#else
      tid = 0
#endif
      call current_second(thr_nwpw_tim1(counter,tid+1))
      return
      end

      subroutine nwpw_timing_end_thr(counter)
      USE omp_lib
      implicit none
      integer counter
#include "nwpw_timing_common.fh"
      integer tid

#ifdef USE_OPENMP
      tid = omp_get_thread_num()
#else
      tid = 0
#endif
      call current_second(thr_nwpw_tim2(counter,tid+1))
     
      thr_times(counter,tid+1) = thr_times(counter,tid+1) 
     >     + (thr_nwpw_tim2(counter,tid+1)-thr_nwpw_tim1(counter,tid+1))
      return
      end


      real*8 function nwpw_timing(counter)
      implicit none
      integer counter
#include "nwpw_timing_common.fh"

      nwpw_timing = times(counter)
      return
      end

      subroutine nwpw_timing_print(msg,time,counter,ttime)
      implicit none
      character*(*) msg
      real*8 time,ttime
      integer counter
#include "stdio.fh"
      if (time>1.0d-9) then
         write(luout,1708) msg,time,time/dble(counter),100*time/ttime
      end if
 1708 FORMAT(A,E14.6,E14.6,F12.1,' %')
      return
      end

      subroutine nwpw_timing_print_thr(msg,id,counter,ttime)
      USE omp_lib
      implicit none
#include "nwpw_timing_common.fh"
      character*(*) msg
      real*8 time,ttime
      integer counter,tid,nthr,used_threads,id

      used_threads = 0
#ifdef USE_OPENMP
      nthr = omp_get_max_threads()
#else
      nthr = 1
#endif
      time = 0
      do tid=1,nthr
        if (thr_times(id,tid)>1.0d-9) then
          time = time + thr_times(id,tid)
          used_threads=used_threads+1
        end if
      end do
      
      if (used_threads>0) then
        time = time / used_threads
        if (time>1.0d-9) then
          write(*,1708) msg,time,time/dble(counter),100*time/ttime
        end if
      end if
 1708 FORMAT(A,E14.6,E14.6,F12.1,' %')
      return
      end

      subroutine nwpw_timing_print_thr_max(msg,id,counter,ttime)
      USE omp_lib
      implicit none
#include "nwpw_timing_common.fh"
      character*(*) msg
      real*8 time,ttime
      integer counter,tid,nthr,id

#ifdef USE_OPENMP
      nthr = omp_get_max_threads()
#else
      nthr = 1
#endif
      time = 0
      do tid=1,nthr
        if (thr_times(id,tid)>1.0d-9) then
          time = max(time,thr_times(id,tid))
        end if
      end do
      
      if (time>1.0d-9) then
        write(*,1708) msg,time,time/dble(counter),100*time/ttime
      end if
 1708 FORMAT(A,E14.6,E14.6,F12.1,' %')
      return
      end




      subroutine nwpw_timing_print_final(oprint,counter)
      implicit none
      logical  oprint
      integer counter

#include "stdio.fh"

      real*8 ttime
c     **** external functions ****
      real*8   nwpw_timing
      external nwpw_timing

#include "nwpw_timing_common.fh"

      call current_second(tf)
      ttime = tf-t0
      if (oprint) then
      write(luout,1809) 'Time spent doing               ',
     >                  'total','step', 'percent'
         call nwpw_timing_print(
     >                 '  total time                 : ',
     >                 ttime,counter,ttime)
         call nwpw_timing_print(
     >                 '  i/o time                   : ',
     >                 nwpw_timing(50),counter,ttime)
         call nwpw_timing_print(
     >                 '  FFTs                       : ',
     >                 nwpw_timing(1),counter,ttime)
         call nwpw_timing_print(
     >                 '  dot products               : ',
     >                 nwpw_timing(2),counter,ttime)
         call nwpw_timing_print(
     >                 '  geodesic                   : ',
     >                 nwpw_timing(10),counter,ttime)
         call nwpw_timing_print(
     >                 '  two-electron Gaussian      : ',
     >                 nwpw_timing(13),counter,ttime)
         call nwpw_timing_print(
     >                 '  ffm_dgemm                  : ',
     >                 nwpw_timing(15),counter,ttime)
         call nwpw_timing_print(
     >                 '  fmf_dgemm                  : ',
     >                 nwpw_timing(16),counter,ttime)
         call nwpw_timing_print(
     >                 '  mmm_dgemm                  : ',
     >                 nwpw_timing(18),counter,ttime)
         call nwpw_timing_print(
     >                 '  m_diagonalize              : ',
     >                 nwpw_timing(17),counter,ttime)

         call nwpw_timing_print(
     >                 '    - m_tredq                : ',
     >                 nwpw_timing(22),counter,ttime)
         call nwpw_timing_print(
     >                 '       - m_tredq_houseq      : ',
     >                 nwpw_timing(26),counter,ttime)
         call nwpw_timing_print(
     >                 '       - m_tredq_houseq_dgemm: ',
     >                 nwpw_timing(28),counter,ttime)
         call nwpw_timing_print(
     >                 '       - m_tredq_dgemm1      : ',
     >                 nwpw_timing(27),counter,ttime)

         call nwpw_timing_print(
     >                 '    - m_getdiags             : ',
     >                 nwpw_timing(23),counter,ttime)
         call nwpw_timing_print(
     >                 '    - m_tqliq                : ',
     >                 nwpw_timing(24),counter,ttime)
         call nwpw_timing_print(
     >                 '    - m_eigsrt               : ',
     >                 nwpw_timing(25),counter,ttime)

         call nwpw_timing_print(
     >                 '  exchange correlation       : ',
     >                 nwpw_timing(4),counter,ttime)
         call nwpw_timing_print(
     >                 '  local pseudopotentials     : ',
     >                 nwpw_timing(5),counter,ttime)
         call nwpw_timing_print(
     >                 '  non-local pseudopotentials : ',
     >                 nwpw_timing(6),counter,ttime)
         call nwpw_timing_print(
     >                 '  hartree potentials         : ',
     >                 nwpw_timing(7),counter,ttime)
         call nwpw_timing_print(
     >                 '  ion-ion interaction        : ',
     >                 nwpw_timing(21),counter,ttime)
         call nwpw_timing_print(
     >                 '  structure factors          : ',
     >                 nwpw_timing(8),counter,ttime)
         call nwpw_timing_print(
     >                 '  phase factors              : ',
     >                 nwpw_timing(20),counter,ttime)
         call nwpw_timing_print(
     >                 '  masking and packing        : ',
     >                 nwpw_timing(9),counter,ttime)
         call nwpw_timing_print(
     >                 '  queue fft                  : ',
     >                 nwpw_timing(30),counter,ttime)
         call nwpw_timing_print(
     >                 '  queue fft (serial)         : ',
     >                 nwpw_timing(31),counter,ttime)
         call nwpw_timing_print(
     >                 '  queue fft (message passing): ',
     >                 nwpw_timing(32),counter,ttime)
         call nwpw_timing_print(
     >                 '  HFX potential              : ',
     >                 nwpw_timing(33),counter,ttime)

         call nwpw_timing_print(
     >                 '  paw gaussian integrals     : ',
     >                 nwpw_timing(34),counter,ttime)
         call nwpw_timing_print(
     >                 '  paw cgaussian integrals    : ',
     >                 nwpw_timing(64),counter,ttime)

         call nwpw_timing_print(
     >                 '  paw atomic coulomb         : ',
     >                 nwpw_timing(35),counter,ttime)
         call nwpw_timing_print(
     >                 '  paw atomic xc              : ',
     >                 nwpw_timing(36),counter,ttime)
         call nwpw_timing_print(
     >                 '  paw gen dEmult/dQlm        : ',
     >                 nwpw_timing(37),counter,ttime)
         call nwpw_timing_print(
     >                 '  paw gen dElocal/dQlm       : ',
     >                 nwpw_timing(38),counter,ttime)
         call nwpw_timing_print(
     >                 '  paw cmp operations         : ',
     >                 nwpw_timing(38),counter,ttime)
         call nwpw_timing_print(
     >                 '  qmmm LJ                    : ',
     >                 nwpw_timing(40),counter,ttime)
         call nwpw_timing_print(
     >                 '  qmmm residual Q            : ',
     >                 nwpw_timing(41),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS InnerLoop          : ',
     >                 nwpw_timing(42),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS Phaze              : ',
     >                 nwpw_timing(43),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS Pipelined FFTs     : ',
     >                 nwpw_timing(44),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS Lagrange           : ',
     >                 nwpw_timing(45),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS Exch Corr          : ',
     >                 nwpw_timing(46),counter,ttime)
         call nwpw_timing_print(
     >                 '  MATHIAS Hpsi               : ',
     >                 nwpw_timing(47),counter,ttime)
         call nwpw_timing_print(
     >                 '  nwpw_ugauss                : ',
     >                 nwpw_timing(48),counter,ttime)
         call nwpw_timing_print(
     >                 '  nwpw_wgauss                : ',
     >                 nwpw_timing(49),counter,ttime)
         call nwpw_timing_print(
     >                 '  nwpw_dwgauss               : ',
     >                 nwpw_timing(55),counter,ttime)
         call nwpw_timing_print(
     >                 '  nwpw_gaunt                 : ',
     >                 nwpw_timing(51),counter,ttime)
         call nwpw_timing_print(
     >                 '  HFX localization           : ',
     >                 nwpw_timing(52),counter,ttime)
         call nwpw_timing_print(
     >                 '  HFX Finding DM columns     : ',
     >                 nwpw_timing(53),counter,ttime)
         call nwpw_timing_print(
     >                 '  HFX DM Cholesky            : ',
     >                 nwpw_timing(54),counter,ttime)
         call nwpw_timing_print(
     >                 '  HFX localized re-gridding  : ',
     >                 nwpw_timing(55),counter,ttime)
         call nwpw_timing_print(
     >                 '  non-local psp FFM          : ',
     >                 nwpw_timing(56),counter,ttime)
         call nwpw_timing_print(
     >                 '  non-local psp FMF          : ',
     >                 nwpw_timing(57),counter,ttime)
         call nwpw_timing_print(
     >                 '  non-local psp FFM A        : ',
     >                 nwpw_timing(58),counter,ttime)
         call nwpw_timing_print(
     >                 '  non-local psp FFM B        : ',
     >                 nwpw_timing(59),counter,ttime)

      end if


      return 
 1808 FORMAT(A,E14.6,E14.6)
 1809 FORMAT(//A,3A14)
      end

