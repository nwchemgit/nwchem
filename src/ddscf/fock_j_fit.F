      subroutine fock_j_fit(nmat, g_dens, g_j)
c
c     $Id$
c
c     Wrapper routine for AO-basis J matrices with fitting
c
c     BGJ - 2/99
c
      implicit none
c
      integer nmat              ! [input] number of J matrices
      integer g_dens(nmat)      ! [input] array of handles to densities
      integer g_j(nmat)         ! [input] array of handles to J matrices
c
c     Local declarations
c
      call fock_j_fit_gen(nmat, g_dens, g_j, 0d0, .false., .true.)
c
      return
      end
c
c
      subroutine fock_j_fit_gen(nmat, g_dens, g_j, Tvec, have_Tvec,
     &     do_init)
c
c     $Id$
c
c     Wrapper routine for AO-basis J matrices with fitting
c
c     BGJ - 2/99
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "../nwdft/include/cdft.fh"
c!!! BGJ
#include "bgj.fh"
#include "rtdb.fh"
c!!! BGJ
c
      integer nmat              ! [input] number of J matrices
      integer g_dens(nmat)      ! [input] array of handles to densities
      integer g_j(nmat)         ! [input] array of handles to J matrices
      logical derfit            ! [input] whether derivative density fit
c                               !  *** from rtdb
c                               !  ******* !!! document here !!!
      double precision Tvec(*)  ! [input] right side vectors for fits
      logical have_Tvec         ! [input] whether right side vectors are given
      logical do_init           ! [input] whether to zero J matrices
c
c     Local declarations
c
      integer imat
      integer rtdb
      integer ntotel, n_batch, n3c_int, n3c_dbl, iwhat_max, 
     &     n_semi_bufs, fd
      integer l_cd_coef, k_cd_coef
      integer idum
      double precision Ecoul
      double precision ddum
      logical IOLGC, old_incore, old_direct
c
      if (bgj_print() .gt. 0)
     &     write(*,*)'--------Entered fock_j_fit-------------'
      rtdb = bgj_get_rtdb_handle()
c      write(*,*)'nmat rtdb geom',nmat,rtdb,geom
c
c     Get fock_j_fit variables
c
      if (.not. rtdb_get(rtdb, 'fock_j:derfit', mt_log, 1,
     &     derfit)) then
         if (bgj_print() .gt. 0)
     &      write(*,*)'*** fock_j_fit: derfit not set: setting to false'
         derfit = .false.       !!!
      endif
      if (bgj_print() .gt. 0)
     &     write(*,*)'fock_j_fit: derfit =',derfit
c
c     Initialize J matrices
c
c      write(*,*)'*** do_init',do_init
      if (do_init) then
         do imat = 1, nmat
c     write(*,*)'Density matrix',imat
c     call ga_print(g_dens(imat))
c            write(*,*)'*** zeroing g_j'
            call ga_zero(g_j(imat))
         enddo
      endif
c
c     Set up scratch space for dft_fitcd and dft_fitvc
c
      if (.not.MA_Push_Get(MT_Dbl,nbf_cd*nmat,'cd_coef',l_cd_coef,
     &     k_cd_coef))
     &     call errquit('fock_j_fit: cannot allocate cd_coef',0, MA_ERR)
c
c     Prepare to call fitting routines
c
      idum = -1
      ddum = -1d0
c!!! Set up for cdinv to be allocated in dft_cdfit
      IOLGC = .true.
      if (noio.eq.1)IOLGC = .FALSE.
      if (ipol.eq.1)then
         nTotEl = 2*noc(1)
      else
         nTotEl = noc(1) + noc(2)
      endif
c!!! Set things up for direct to begin with - can refine to semi-direct
c!!! after things are working
      old_incore = incore
      old_direct = direct
      incore = .false.          !!!
      direct = .true.           !!!
      n_batch = 0               !!!
      n3c_int = 0               !!!
      n3c_dbl = 0               !!!
      fd = 0                    !!!
      n_semi_bufs = 0           !!!
      iwhat_max = -1            !!!
c
c     Determine the fitting coefficients
c
      call dft_fitcd(nmat,dbl_mb(k_cd_coef), idum, Ecoul, 
     &     g_dens, nTotEl, n_batch, n3c_int,
     &     idum, n3c_dbl, iwhat_max, n_semi_bufs, fd, IOLGC, 
     &     idum, derfit,  Tvec,
     &     have_Tvec)
c
c     Compute the J matrices using the fitting coefficients
c
      call dft_fitvc(dbl_mb(k_cd_coef), idum, Ecoul,  g_j,
     &     n_batch, n3c_int, idum, n3c_dbl,
     &     iwhat_max, n_semi_bufs, fd, IOLGC, 
     &     derfit, nmat, do_init)
c
c     Clean up
c
      if (.not.ma_pop_stack(l_cd_coef))
     &     call errquit('fock_j_fit: cannot pop stack',0, MA_ERR)
      direct = old_direct
      incore = old_incore
c
      if (bgj_print() .gt. 0) then
         do imat = 1, nmat
            write(*,*)'Fitted J matrix',imat
            call ga_print(g_j(imat))
         enddo
         write(*,*)'--------Leaving fock_j_fit-------------'
      endif
c
      return
      end
