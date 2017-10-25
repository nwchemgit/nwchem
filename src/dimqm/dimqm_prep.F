      subroutine dimqm_prep(rtdb, om)
c
c     Author: Justin Moore
c
c     Called from: dimqm_main.F cphf_solve2.F cphf_solve3.F 
c                  aoresponse_driver.F aoresponse_driver_new.F
c                  dimqm_polarizability.F
c                  
c
c     Subroutines called: None
c
c     Builds the interaction matrix for the DIM system.   Called
c     at the start of every DIM/QM calculation.
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
c     Input variables
c
      integer rtdb        ! RTDB handle
      double precision om ! Frequency, zero or negative means static
c
c     Local variables
      integer stat
      integer nLength
      integer g_matrix
      integer g_indx
      integer m, i, itype
      integer mLo, mUp, mOnNode
      double precision polr
      double complex   polc
      double complex   die
      integer tx(nDIM) ! Type index
      double precision sPol(nDIMTypes)
      double precision bound(nDIMTypes)
      double precision drudeterms(2, nDIMTypes)
      double complex Drude
      double complex ExperimentalDielectric
      integer id
      double complex expdies(nDIMTypes)
      logical expdiesRead(nDIMTypes)
      logical hasExpDie(nDimTypes)
c
c     Common variables used from dimqm.fh
c
c      integer alg, g_DIM_diag
c      logical ldebug, lpim
c      integer nOrder, nDIM, nDIMTypes
c
      expdies = ZERO_C
      expdiesRead = .false.
      id = ga_nodeid()
      if(id.eq.0 .and. ldebug) write(luout,*) "Start dimqm_prep"
      call pphilo(id, nDIM, mLo, mUp, mOnNode)
c
c   Pull needed info from RTDB
      if(.not. rtdb_get(rtdb,'dimpar:sPol', mt_dbl, nDIMTypes, sPol))
     $  call errquit('dimpar:sPol rtdb get failed', 0, RTDB_ERR)
c
      if(.not. rtdb_get(rtdb,'dimpar:typeindex', mt_int, nDIM, tx))
     $  call errquit('dimpar:type index rtdb get failed', 0, RTDB_ERR)
c
      if(.not. rtdb_get(rtdb,'dimpar:bound', mt_dbl, nDIMTypes, bound))
     $  call errquit('dimpar:bound rtdb get failed', 0, RTDB_ERR)
c
      if(.not. rtdb_get(rtdb,'dimpar:drude', mt_dbl, nDIMTypes*2,
     $                                    drudeterms))
     $  call errquit('dimpar:drude rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:hasexp', mt_log, nDIMTypes,
     $                                    hasExpDie))
     $  call errquit('dimpar:hasexp rtdb get failed', 0, RTDB_ERR)
     
c
c   =========================================
c    Allocate GAs for the interaction matrix
c   =========================================
c
c     Exact Solver
      if(alg ==  0) then   
        nLength = nOrder*(nOrder+1)/2
c
c        Complex matrix allocation
        if(om > ZERO) then   
          if(.not.ga_create(mt_dcpl, nLength, 1, 'dimMatrix_FD',
     $                      -1, 1, g_matrix))
     $      call errquit('dimMatrix_FD: ga_create failed',9,GA_ERR)
c
c        DP matrix allocation
        else
          if(.not.ga_create(mt_dbl, nLength, 1, 'dimMatrix_S', -1, 1, 
     $                      g_matrix))
     $      call errquit('dimMatrix_S: ga_create failed',9,GA_ERR)
        end if
c
c        Allocate index array
        if(.not. ga_create(mt_int, nOrder, 1, 'dimMatrix_indx', -1, 1,
     $                     g_indx))
     $    call errquit('dimMatrix_indx: ga_create failed',0,GA_ERR)
c
c        Zero both arrays
        call ga_zero(g_matrix)
        call ga_zero(g_indx)
c
c        Place handles into RTDB
        if(.not. rtdb_put(rtdb,'dimqm:g_matrix',mt_int,1,g_matrix))
     $    call errquit('g_matrix: rtdb put failed',0,RTDB_ERR)
        if(.not. rtdb_put(rtdb,'dimqm:g_indx',mt_int,1,g_indx))
     $    call errquit('g_indx: rtdb put failed',0,RTDB_ERR)
c
c      Iterative Solver
      else
c        Complex diagonal allocation
        if(om > ZERO) then
          if(.not. ga_create(mt_dcpl, nOrder, 1, 'dimqm:diagFD',
     $                       -1, 1, g_DIM_diag))
     $      call errquit('dimqm:alphainvFD: create failed',9,GA_ERR)
c        DP diagonal allocation
        else
          if(.not. ga_create(mt_dbl, nOrder, 1, 'dimqm:diagS',
     $                       -1, 1, g_DIM_diag))
     $      call errquit('dimqm:alphainvFD: create failed',9,GA_ERR)
        end if
c
c        Zero the diagonal
        call ga_zero(g_DIM_diag)
      end if
c
c   =========================================================================
c   Build the interaction matrix (only the diagonal for the iterative solver) 
c   =========================================================================
c
c      Iterative solver
      if(alg > 0) then
        call ga_init_fence()
c
c        Frequency dependent calculation
        if(om > ZERO) then
c          Loop over DIM atoms
          do m = mLo, mUp 
            itype = tx(m)
            die = ExperimentalDielectric(rtdb, itype, om, hasExpDie,
     $                                   expDiesRead, expDies)
            if (die /= ZERO_C) then
              ! Size correction
            else
              die = CMPLX(bound(itype), ZERO)
              if(drudeterms(1,itype) > ZERO) then ! Subtract Drude term if given
              die = die -
     $              Drude(om, drudeterms(1,itype), drudeterms(2,itype))
              endif
            end if
            polc = ONE_C/(((die - ONE_C) / (die + TWO_C)) * sPol(itype))
            i = 3*(m-1)
            call ga_put(g_DIM_diag, i+1, i+1, 1, 1, polc, 1)
            call ga_put(g_DIM_diag, i+2, i+2, 1, 1, polc, 1)
            call ga_put(g_DIM_diag, i+3, i+3, 1, 1, polc, 1)
          end do
c
c        Static calculation
        else
c          Loop over DIM atoms
          do m = mLo, mUp
            itype = tx(m)
            i = 3*(m-1)
            polr = ONE / sPol(itype)
            call ga_put(g_DIM_diag, i+1, i+1, 1, 1, polr, 1)
            call ga_put(g_DIM_diag, i+2, i+2, 1, 1, polr, 1)
            call ga_put(g_DIM_diag, i+3, i+3, 1, 1, polr, 1)
          end do
        end if
        call ga_fence()
c
c      Exact Solver
      else
        call errquit("dimqm direct solver not implemented", 1,
     $               INPUT_ERR)
      end if
      dimqm_on = .true.
      lfirst   = .true.
      end subroutine dimqm_prep

      function Drude (omega, plasma, damp) result (drd)
      implicit none
#include "dimqm_constants.fh"
c
c   =====================================
c   Calculates the complex Drude function
c   =====================================
c
      double precision omega
      double precision plasma
      double precision damp
      double complex drd

      drd = plasma**2 / (omega * (omega + I_C * damp))
      end function Drude

      function ExperimentalDielectric(rtdb, atomtype, om, hasExpDie,
     $                                expDieRead, expDies) result (die)
      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "dimqm.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"
#include "errquit.fh"

      integer rtdb
      integer atomtype
      double precision om
      double complex die
      logical hasExpDie(nDIMTypes)
      double complex expDies(nDIMTypes)
      logical expDieRead(nDIMTypes)
      
      integer iexp
      double precision rawexpdie(3, 10000)
      double precision splined(2, 10000)
      double complex interpolate_complex1

      if (.not. hasExpDie(atomtype)) then
        die = ZERO_C
      else if (expDieRead(atomtype)) then
        die = expDies(atomtype)
      else
        if(.not.rtdb_get(rtdb,'dimpar:iexp'//char(atomtype),
     $                   mt_int, 1, iexp))
     $      call errquit('iexp: rtdb get failed', atomtype, RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimpar:expdata'//char(atomtype), mt_dbl,
     $                   3*iexp, rawexpdie(1:3,1:iexp)))
     $      call errquit('expdata: rtdb get failed', atomtype, RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimpar:splinedata'//char(atomtype),
     $                   mt_dbl, 2*iexp, splined(1:2,1:iexp)))
     $     call errquit('splinedata: rtdb get failed',atomtype,RTDB_ERR)
        expDies(atomtype) = interpolate_complex1(iexp,
     $                                           rawexpdie(1,1:iexp),
     $                                           rawexpdie(2,1:iexp),
     $                                           rawexpdie(3,1:iexp),
     $                                           splined(1,1:iexp),
     $                                           splined(2,1:iexp),
     $                                           HART2NM / om)
        expDieRead(atomtype) = .true.
        die = expDies(atomtype)
      end if
      end function ExperimentalDielectric
