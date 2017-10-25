      subroutine dimqm_polarizability(rtdb, om)
c  ==================================================================
c  Calculate the standalone polarizability of the DIM system.  Here,
c  we apply a unit field in each cartesian direction and determine
c  the DIM system's response.
c
c  Called from: property/aoresponse_driver.F
c  ==================================================================
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input variables
      integer rtdb        ! RTDB handle
      double precision om ! Frequency of incident light (zero if static calculation)
c
c     MA handles
      integer k_fld, l_fld
      integer l_dimxyz, k_dimxyz ! DIM coordinates handles
      integer l_muind,  k_muind  ! DIM dipoles handles
c
c     Local variables
      integer idir, m, i3
      integer icmplx
      double precision pol(3,3,2)
      double precision time
      logical stat
      integer id
      double precision dsum
      external dsum
      logical seed_save
c
c   Common variables used
c      logical lcpim, lcmplx, ldebug
c      integer nOrder, nDIM
c
c   ==========
c   Initialize 
c   ==========
c
      i3 = 3 * nDIM
      id = ga_nodeid()
      icmplx = 1
      if(om > ZERO) icmplx = 2
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "Begin Standalone Polarizability Routine"
c
c   ========================
c   Allocate arrays and fill
c   ========================
c
c     Coords
      if(.not.ma_push_get(mt_dbl,i3,'dim coords',l_dimxyz,k_dimxyz))
     $   call errquit('dimqm_pol malloc k_xyz failed',1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimpar:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $   call errquit('dimpar:coords rtdb get failed',1,RTDB_ERR)
c
c     Induced Dipoles
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dim muind',
     $                                            l_muind,k_muind))
     $   call errquit('dimqm_pol malloc k_muind failed',1,MA_ERR)
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dim:muold',
     $                                            l_muold,k_muold))
     $   call errquit('dimqm_pol malloc k_muind failed',1,MA_ERR)
c
c     Applied Field
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dim fld',l_fld,k_fld))
     $   call errquit('dimqm_pol malloc k_fld failed',1,MA_ERR)
c
c   Zero polarizability
      pol = ZERO
c
c   ===================================================
c   Loop over each unit field and sum up the dipoles to
c   get the polarizability
c   ===================================================
c
      time = util_timer()
      lfirst = .true.
      seed_save = dimqm_noseed
      dimqm_noseed = .true.
      do idir = 1, 3
c
c           Zero field
          call dfill(i3*icmplx, ZERO, dbl_mb(k_fld), 1)
          call dfill(i3*icmplx, ZERO, dbl_mb(k_muold), 1)
c
c           Unit field in specified direction
          call dfill(nDIM, ONE, dbl_mb(k_fld+idir-1), 3)
          call dimqm_f2d(rtdb, dbl_mb(k_fld),
     $                   dbl_mb(k_muind), dbl_mb(k_muold),
     $                   dbl_mb(k_dimxyz), icmplx, 's', ' ',.false.)
c           Sum up dipoles to get polarizability
          pol(1,idir,1) = dsum(nDIM, dbl_mb(k_muind), 3)
          pol(2,idir,1) = dsum(nDIM, dbl_mb(k_muind+1), 3)
          pol(3,idir,1) = dsum(nDIM, dbl_mb(k_muind+2), 3)
          if(icmplx > 1) then
            pol(1,idir,2) = dsum(nDIM, dbl_mb(k_muind+nOrder), 3)
            pol(2,idir,2) = dsum(nDIM, dbl_mb(k_muind+1+nOrder), 3)
            pol(3,idir,2) = dsum(nDIM, dbl_mb(k_muind+2+nOrder), 3)
          end if
      end do
c
c   Print polarizability
      if(id.eq.0) then
        if(om > ZERO) then
          write(luout,5000) "Frequency-Dependent Real",
     $                      pol(1,1,1), pol(1,2,1), pol(1,3,1),
     $                      pol(2,1,1), pol(2,2,1), pol(2,3,1),
     $                      pol(3,1,1), pol(3,2,1), pol(3,3,1) 
          write(luout,*)
          write(luout,5000) "Frequency-Dependent Imag",
     $                      pol(1,1,2), pol(1,2,2), pol(1,3,2),
     $                      pol(2,1,2), pol(2,2,2), pol(2,3,2),
     $                      pol(3,1,2), pol(3,2,2), pol(3,3,2)
          write(luout,*)
        else
          write(luout,5000) "Static",
     $                      pol(1,1,1), pol(1,2,1), pol(1,3,1),
     $                      pol(2,1,1), pol(2,2,1), pol(2,3,1),
     $                      pol(3,1,1), pol(3,2,1), pol(3,3,1)
          write(luout,*)
        end if
        call util_flush(LuOut)
      end if
      dimqm_noseed = seed_save
c
c   ==============================
c   Deallocate arrays and clean up
c   ==============================
c
      stat = ma_chop_stack(l_dimxyz)
c  jbecca START: dont think we need to create and destroy int mat here.
c        this is called from ao_driver_new, where it is already made
c      call dimqm_tidy()
c  jbecca END
      if(id .eq. 0) write(LuOut,*)
     $  "Total time:", util_timer() - time
 5000 format("Polarizability tensor for DIM system: ", A ,/
     $       1x,'             X              Y              Z',/
     $       1x,'-----------------------------------------------',/
     $       1x,'X ',3f15.7,/
     $       1x,'Y ',3f15.7,/
     $       1x,'Z ',3f15.7,/
     $       1x,'-----------------------------------------------')


      end subroutine dimqm_polarizability

