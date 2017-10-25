      subroutine dimqm_main(rtdb, basis, geom, g_dens, lastcy)
c
c     Author: Justin Moore
c             Jeff Becca
c
c     Called from: dft_scf.F, rohf_fock.F 
c
c     Subroutines called: dimqm_prep, dimqm_EqmE, dimqm_EqmN, dimqm_f2d
c                         dimqm_seed_init_tolerance, dimqm_check_dipoles,
c                         dimqm_energy, dimqm_printAtoms, dimqm_printDipole,
c                         dimqm_printAtomicDipoles, dimqm_tidy 
c
c     Main driver for the DIM/QM module.  Called every SCF cycle.
c    
c

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c   
c     Input Variables
      integer rtdb      ! RTDB handle
      integer basis     ! Basis set handle
      integer geom      ! Geometry handle
      integer g_dens    ! QM density GA handle (we only work with the total density
      logical lastcy    ! Last SCF iteration flag
c
c     Local variables
      integer i3, m, i                  ! Looping and convenience variables
      double precision EqmE(3, nDIM)    ! Electric field from QM electrons felt at each DIM atom
      double precision Etot(3, nDIM)    ! Total QM electric field (electrons + nuclei) felt at each DIM atom
      integer l_fin,  k_fin             ! Finite field
      double precision dimxyz(3, nDIM)  ! DIM coordinates
      double precision muind(3, nDIM)   ! DIM induced dipoles
      double precision max_err, calcErr ! Error in dipoles
      logical stat                      ! Dummy for function returns
      integer id                        ! Processor ID
      double precision dime             ! DIM/QM interaction energy function
      integer a
      double precision dx, dy, dz       ! Total induced dipoles moment for the DIM system
      double precision efld_int
c
c     Common variables used from dimqm.fh 
c
c      integer l_eqmn, k_eqmn
c      integer l_muold, k_muold
c      logical lfirst, ldebug
c      integer nDIM, dimqm_rtdb
c      logical lefield
c      double precision finiteFieldX, finiteFieldY, finiteFieldZ
c
c   ===========
c   Initialize
c   =========== 
      i3 = 3 * nDIM
      id = ga_nodeid()
      dimqm_rtdb = rtdb
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Entering DIM/QM Main"
      end if
c
c   ================
c   Allocate arrays
c   ================
c
      EqmE  = ZERO
      muind = ZERO
      Etot  = ZERO
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Allocated local arrays"
      end if
c
c      Read in coords from RTDB
      if(.not. rtdb_get(rtdb,'dimpar:coords',mt_dbl,i3,dimxyz))
     $   call errquit('dimqm:xyz rtdb get failed',1,RTDB_ERR)
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Pulled coords"
c        write(*,*) "GEOM:", geom
      end if
c     Read in whether or not this is first cycle
c     It seems the house cleaning made all of this rtdb dependent now
      if(.not. rtdb_get(rtdb,'dimqm:lfirst',mt_log, 1, lfirst))
     $   call errquit('dimqm_main:lfirst get failed', 1, RTDB_ERR)
c
c   ===========
c   First Cycle
c   ===========
c
      if(lfirst) then
c
c        Allocate memory for QM nuclear contribution on the heap (PERSISTANT)
        if(.not.ma_alloc_get(mt_dbl,i3,'dimqm eqmn',l_eqmn,k_eqmn))
     $    call errquit('dimqm_main malloc k_eqmn failed', 1, MA_ERR)
c
c        Allocate memory for last cycle's dipoles on the heap (PERSISTANT)
        if(.not.ma_alloc_get(mt_dbl,i3,'dimqm muold',l_muold,k_muold))
     $    call errquit('dimqm_main malloc k_muold failed',1,MA_ERR)
c
        if(id.eq.0.and.ldebug) then
          write(luout,*) "Allocated heap arrays"
        end if
c        Zero both arrays
        call dfill(i3, ZERO, dbl_mb(k_muold), 1)
        call dfill(i3, ZERO, dbl_mb(k_eqmn), 1)
c        Build interaction matrix
        call dimqm_prep(rtdb, -ONE)
c        Generate electric field due to QM nuclei.  This doesn't change so we only do it once.
        call dimqm_EqmN(rtdb, geom, dimxyz,
     $                  dbl_mb(k_eqmn))
        if(id.eq.0.and.ldebug) then
          write(luout,*) "DIM/QM First Cycle Prep Done"
        end if
        call dimqm_seed_init_tolerance(1.0d-4)
      end if
c
c      Generate E-field from QM electrons
      call dimqm_EqmE(rtdb, g_dens, geom, basis, EqmE,
     $                dimxyz)
c
c      Add finite field if present
      if(lefield .and. lclfld) then
        !if(id.eq.0 .and. .not. lastcy) then 
        !  write(luout,*) "Adding finite field to Eqm"
        !  write(luout,*) finiteFieldX, finiteFieldY, finiteFieldZ
        !end if
        if(.not.ma_push_get(mt_dbl,i3,'dimqm finite',l_fin,k_fin))
     $   call errquit('dimqm_main malloc k_finite failed',1,MA_ERR)
        call dfill(nDIM, finiteFieldX, dbl_mb(k_fin)  , 3)
        call dfill(nDIM, finiteFieldY, dbl_mb(k_fin+1), 3)
        call dfill(nDIM, finiteFieldZ, dbl_mb(k_fin+2), 3)
c
c        Add finite field for finite-difference local field calculations
c  jbecca TODO: make a separate variable for the finitefield so 
c              the energy can be separated out for static local field
        call daxpy(i3, ONE, dbl_mb(k_fin), 1, EqmE, 1)
        if(.not. ma_pop_stack(l_fin)) call
     $       errquit('ma_pop failed l_fin',1,MA_ERR)
      end if
c
c      Add electronic and nuclear efields into etot
      Etot = EqmE
      call daxpy(i3, ONE, dbl_mb(k_eqmn), 1, Etot, 1)
c      
c   =======================================================
c   Generate induced dipoles (and charges) on the DIM atoms
c   =======================================================
c
c  jbecca: Adding in option to add external field induced dipoles
      if (lefield .and. lclfld) then
         call dimqm_f2d(rtdb, Etot, muind, dbl_mb(k_muold), dimxyz, 1,
     $               'g', ' ',.true.)
      else
         call dimqm_f2d(rtdb, Etot, muind, dbl_mb(k_muold), dimxyz, 1,
     $               'g', ' ',.false.)
      endif
c
c      Write induced dipoles and charges to disk
      if(.not.rtdb_put(rtdb,'dimqm:muind',mt_dbl,i3, muind))
     $  call errquit('dimqm_main muind put failed', 1, RTDB_ERR)
c
c      Calculate DIM/QM interaction energy
c      This is needed in the RTDBdd
      edimqm = dime(EqmE, dbl_mb(k_eqmn), muind, nDIM)
c  jbecca START: need to add in the interaction energy with the efield
c              if the efield is present
c  NOTE: This is not needed for now, since the finitefield is added to
c        EqmE, the interaction energy will also be calculated as part
c        of -HALF * SUM(EqmE * muind)
c
c      if (lefield .and. lastcy) then
c         write(luout,*)'adding in DIM/QM interaction energy with efield'
c         efld_int    =  ZERO
c         efld_int    =  SUM(muind(1,:) * finiteFieldX)
c     $               +  SUM(muind(2,:) * finiteFieldY)
c     $               +  SUM(muind(3,:) * finiteFieldZ)
c         efld_int    =  -HALF * efld_int
c         write(luout,*)'edimqm before efld_int',edimqm
c         write(luout,*)'efld_int',efld_int
c         edimqm      =  edimqm + efld_int
c         write(luout,*)'edimqm after efld_int',edimqm
c      endif       !lefield
c  jbecca END

      if (.not. rtdb_put(rtdb, 'dimqm:edimqm', MT_DBL, 1, edimqm))
     $  call errquit('dimqm_main edimqm put failed', 1, RTDB_ERR)
c
c
      dx = SUM(muind(1,:))
      dy = SUM(muind(2,:))
      dz = SUM(muind(3,:))
      if(id.eq.0 .and. ldebug .and. .not. lastcy) then
        write(luout,*) "Induced dipoles:"
        write(luout,*) dx, dy, dz
        write(luout,*) ''
      end if
c
c      Determine error in dipoles
      max_err = calcErr(i3, dbl_mb(k_muold), muind)
      if(id.eq.0 .and. ldebug .and. .not. lastcy) then
        write(LuOut,*) "Max error in dipoles:", max_err
        write(LuOut,*) ''
      end if
c     call dimqm_check_dipoles(2.0d-1, max_err)
      call dimqm_check_dipoles(dimtol, max_err)
c
c      Copy new dipoles into muold
      call dcopy(i3, muind, 1, dbl_mb(k_muold), 1)
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $  "Saved dipoles as old dipoles"
c
c      No longer the first cycle!
      lfirst = .false.
      if (.not. rtdb_put(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
     $   call errquit('dimqm_main: lfirst put failed', 1, RTDB_ERR)
c
c      Last cycle only
      if(lastcy) then
        if(id.eq.0) then
          write(luout,*)
     $    "                                    DIM/QM Results"
          write(luout,*)
     $    "                                    --------------"
          write(luout,*) ""
          write(luout,*) ""
        end if

c        Calculate energy terms
        call dimqm_energy(nDIM, muind, EqmE, dbl_mb(k_eqmn))
c        Print coordinates
        call dimqm_printAtoms(rtdb, dimxyz)
c        Print total dipole moment
        call dimqm_printDipole(rtdb, muind, 1, .true.)
c        Print atomic dipoles (and charges)
        call dimqm_printAtomicDipoles(rtdb, .false., .false., '')
c        Destroy persistant arrays
        if(.not.ma_free_heap(l_muold))
     $    call errquit("dimqm_main: can't free l_muold", MA_ERR)
        if(.not.ma_free_heap(l_eqmn))
     $    call errquit("dimqm_main: can't free l_eqmn", MA_ERR)
        call dimqm_tidy()
      end if
c
      if(id.eq.0.and.ldebug) write(LuOut,*) "End dimqm_main"
      call util_flush(LuOut)
      call ga_sync()
      end subroutine dimqm_main

      function calcErr(n, old, new)
      implicit none
#include "stdio.fh"
      integer n
      double precision old(n)
      double precision new(n)
      double precision calcErr
      double precision temp(n)
c
      temp(:) = new(:) - old(:)
      calcErr = MAXVAL(temp)
      if (ABS(MINVAL(temp)) > calcErr) then
         calcErr = ABS(MINVAL(temp))
      end if
      return
c
      end function calcErr

      function dime(eqme, eqmn, muind, nAtoms)
      implicit none
#include "dimqm_constants.fh"
      integer nAtoms
      double precision eqme(3,nAtoms)
      double precision eqmn(3,nAtoms)
      double precision muind(3,nAtoms)
      double precision dime
      
      dime = -HALF * SUM(muind * eqme)
      dime = dime +(-HALF)* SUM(muind * eqmn)
      return
      end function dime

      subroutine print_dim(array, nel, label)
      implicit none
#include "stdio.fh"
      integer nel
      double precision array(3,nel)
      character*4 label
      write(luout,*) label
      write(luout,*) array
      call util_flush(LuOut)
      end subroutine print_dim
