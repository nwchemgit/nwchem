      subroutine dimqm_iterative(rtdb, fld, muold, xyz, icmplx, iscmplx,
     $                            dir, pm)
c      use constants
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input Variables
      integer rtdb
      integer dtype 
      integer icmplx
      double precision xyz(3,nDIM)
      double precision fld(nOrder, icmplx)
      double precision muold(nOrder, icmplx)
      logical iscmplx
      character dir, pm
c
c     GA handles
      integer g_work(8)
c
c     Common variables
c      double precision dimtol
c      logical ldebug
c      integer nDIM, nOrder
c
c     Local variables
      integer ierr
      external matvecReal
      external matvecComplex
      double precision time
      logical stat
      integer i, iOrder
      integer id
      double precision err
      integer k_cmplxFld, l_cmplxFld
      integer k_cmplxMuold, l_cmplxMuold
      double precision seed(nOrder, icmplx)
c
c     Get node ID
      id = ga_nodeid()
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "Begin dimqm_iterative"
c
c     Get data type of the external field
      dtype = mt_dbl
      if(iscmplx) dtype = mt_dcpl
c
c     Allocate the work arrays
      if(.not. ga_create(dtype, nOrder, 1, 'g_work1', -1, 1,
     $  g_work(1)))
     $  call errquit('g_work1: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work2', -1, 1,
     $  g_work(2)))
     $  call errquit('g_work2: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work3', -1, 1,
     $  g_work(3)))
     $  call errquit('g_work3: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work4', -1, 1,
     $  g_work(4)))
     $  call errquit('g_work4: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work5', -1, 1,
     $  g_work(5)))
     $  call errquit('g_work5: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work6', -1, 1,
     $  g_work(6)))
     $  call errquit('g_work6: ga create failed',0,GA_ERR)
      if(iscmplx) then
        if(.not. ga_create(dtype, nOrder, 1, 'g_work7', -1, 1,
     $    g_work(7)))
     $    call errquit('g_work7: ga create failed',0,GA_ERR)
      end if
c     
c     Some prep before calling the solver
      nIter = MAX(nOrder, 1000)
      err = dimtol

c
c     Call the solver with the proper matrix vector multiplication
c     routine
      time = util_timer()
c      write(luout,*) "noseed?", dimqm_noseed
      if(iscmplx) then
c       
c           Transform fld vector from real to complex
        if(.not.ma_push_get(mt_dcpl, nOrder,'dimqm:cmplxFld',
     $                      l_cmplxFld, k_cmplxFld))
     $   call errquit('dimqm_f2d malloc k_cmplxArr failed', 1, MA_ERR)
        if(.not.ma_push_get(mt_dcpl, nOrder,'dimqm:cmplxMuold',
     $                      l_cmplxMuold, k_cmplxMuold))
     $   call errquit('dimqm_f2d malloc k_cmplxMuold failed', 1, MA_ERR)
        do iOrder = 1, nOrder
          dcpl_mb(k_cmplxFld+iOrder-1) = DCMPLX(fld(iOrder,1),
     $                                          fld(iOrder,2))
          dcpl_mb(k_cmplxMuold+iOrder-1) = DCMPLX(muold(iOrder,1),
     $                                            muold(iOrder,2))
        end do
        if(dimqm_noseed) then
          if(id.eq.0 .and. ldebug) write(luout,*) "DIM/QM QMR"
          call qmr_complex(rtdb, nOrder, g_work, nIter, err, ierr, 
     $                     dcpl_mb(k_cmplxFld), matvecComplex, xyz,
     $                     ldebug)
        else if(dimqm_seeded) then
          if(id.eq.0.and.ldebug) write(luout,*) "Using previous dipoles"
          call seeded_qmr_complex(rtdb, nOrder, g_work, nIter, err,
     $                            ierr, dcpl_mb(k_cmplxFld),
     $                            matvecComplex, xyz, ldebug, dir, pm)
        else
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM Seeding QMR"
          call qmr_seed_complex(rtdb, nOrder, g_work, nIter, err, ierr, 
     $                          dcpl_mb(k_cmplxFld), matvecComplex, xyz,
     $                          ldebug, dir, pm)
        end if
c
c       Transform fld vector from complex to real
        do iOrder = 1, nOrder
          fld(iOrder, 1) =  REAL(dcpl_mb(k_cmplxFld+iOrder-1))
          fld(iOrder, 2) = AIMAG(dcpl_mb(k_cmplxFld+iOrder-1))
        end do
      else
        if(dimqm_noseed) then
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM QMR"
          call qmr_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $                  fld(:,1), matvecReal, xyz, ldebug)
        else if(dimqm_seeded) then 
          if(id.eq.0.and.ldebug) write(luout,*) "Using previous dipoles"
          call seeded_qmr_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $                         fld(:,1), matvecReal, xyz, ldebug, dir)
        else
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM Seeding QMR"
          call qmr_seed_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $                        fld(:,1), matvecReal, xyz, ldebug, dir)
        end if
      end if

      time = util_timer() - time
      select case (ierr)
      case (0)
c       Success
        if(id .eq. 0 .and. ldebug) then
          write(LuOut,*)"DIM/QM dipoles converged in",nIter,"iterations"
          write(LuOut,*) "Took", time, "seconds"
          write(LuOut,*) ''
        end if
      case (2)
        call errquit('Error in arguments for QMR routine',ierr,GA_ERR)
      case (3)
        call errquit('Iterative solver did not converge',ierr,GA_ERR)
      case (4)
        call errquit('Nonsingular system in iterative solver',0,GA_ERR)
      case (5)
        call errquit('A-invarient subspace found in iterative solver',0,
     $                  GA_ERR)
      end select
c
c     Clean up
      do i=1,6 
        stat = ga_destroy(g_work(i))
      end do
      if(iscmplx) then
        stat = ga_destroy(g_work(7))
      end if
      if(iscmplx) then
        if(.not.ma_chop_stack(l_cmplxFld))
     $      call errquit('pop cmplxFld in iterative', 1, MA_ERR)
      end if
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "End dimqm_iterative"
      end subroutine dimqm_iterative
