#ifdef GWCMPLX
#define ISOGW complex(kind=c_double_complex)
#define GWTYP double complex
#else
#define ISOGW real(kind=c_double)
#define GWTYP double precision
#endif 

#define USE_ALLOC 1
      subroutine gw_cdgw_buildiw(pars,eia,x,wgt,iw,iwidx,nmo,nri,ngl,
     $                           mynmo,maxpoles,maxqp,ipol)

C**********************************************************************
C*     This subroutine precomputes the screened Coulomb matrix        *
C*     elements over the imaginary axis. This corresponds to          *
C*     Equation (41) of:                                              *
C*                                                                    *
C*     J. Chem. Theory Comput. 17, 7504 (2021).                       *
C*                                                                    *
C*     The array iW stores all matrix elements corresponding to       *
C*     the MO pairs given in iWidx.                                   *
C**********************************************************************

#ifdef USE_ALLOC
      USE ISO_C_BINDING
#endif
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nmo,nri,ngl,mynmo,maxpoles,maxqp,ipol

      double precision,intent(in)    :: eia(maxpoles,*),x(ngl),wgt(ngl)
      integer,intent(in)             :: iwidx(0:mynmo,maxqp,*)
      GWTYP,intent(out)   :: iw(mynmo,ngl,maxqp,*)

      character(len=12),parameter :: pname = 'gw_buildiw: '

      logical ok,main,diag
      integer molo,mohi,ldpi,sizepi,l_sol,k_sol
      integer ilo,ihi,me,ld,mynpoles,igl,ipole,info
      integer pilo,pihi,k_pi,l_pi,idx,idy,l_tmp,tmp,l_vx,k_vx
      integer imo,jmo,kmo,lmo,nprocs,kOV1,kOV2,isp,addr,idz,kOO(2)
      integer kVV(2),kOV(2),mynvir,mynocc,nlocal,nmod
      double precision w2,ecrpa,ecgm
      double precision, external :: gw_ddot
      double complex, external :: gw_ddot_cmplx
#ifdef USE_ALLOC
      ISOGW, allocatable    :: factor(:,:)
      integer(kind=c_int)   :: stat_alloc
#else
      GWTYP :: factor(maxpoles,ipol)
#endif
      double complex :: imw
      double complex :: ieta
      double precision t0, t1


#ifdef USE_OPENMP
      integer iMaxThreads,ithread
      integer,external :: omp_get_max_threads,omp_get_thread_num
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif
#ifdef USE_ALLOC
      allocate(factor(maxpoles,ipol),STAT = stat_alloc)
      if(stat_alloc.ne.0)
     &     call errquit(' gw_buildiw: allocate failed for size',
     M     maxpoles*ipol, MA_ERR)
#endif
      call ga_sync()

      ieta = dcmplx(0d0,pars%eta)      
      me = ga_nodeid()
      nprocs = ga_nnodes()
      main = me.eq.0
      pars%ecrpa = 0d0
      pars%ecgm = 0d0
C
C     Decide the leading dimension of the \Pi array according to
C     the rectangular full-packed (RFP) format
C     
      if (mod(nri,2).eq.0) then
        ldpi = nri + 1
      else
        ldpi = nri
      endif
      sizepi = (nri*(nri+1))/2
C
C     If diag is .true., the inverse of the dielectric matrix
C     will be obatined using eigendecomposition in GA.

      diag = pars%diag

C
C     Allocate local fields
C  
#ifdef GWCMPLX
      if (.not.ma_alloc_get(mt_dcpl,nri**2,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dcpl,nri*mynmo,'sol',l_sol,k_sol))
     $  call errquit(pname//'can''t allocate sol',0,MA_ERR) 
#else
      if (.not.ma_alloc_get(mt_dbl,sizepi,'Pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate Pi',0,MA_ERR) 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'sol',l_sol,k_sol))
     $  call errquit(pname//'can''t allocate sol',0,MA_ERR) 
      if (pars%rpa) then
        if (.not.ma_alloc_get(mt_dbl,sizepi,'vx',l_vx,k_vx))
     $  call errquit(pname//'can''t allocate vx',0,MA_ERR)
      endif
#endif 
      if (.not.ma_alloc_get(mt_dbl,nri*mynmo,'temp',l_tmp,tmp))
     $  call errquit(pname//'can''t allocate tmp',0,MA_ERR) 

      kOO = 1
      kOV = 1
      kVV = 1

C
C     Pointers to alpha spin MO ERIs stored locally for a given MPI rank.
C
      if (pars%ovlo(1).le.pars%ovhi(1) .and.pars%ovlo(1).gt.0) then
        call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $               kOV(1),ld)
      endif
      if (pars%oolo(1).le.pars%oohi(1) .and.pars%oolo(1).gt.0) then
        call ga_access(pars%g_erioo(1),1,nri,pars%oolo(1),pars%oohi(1),
     $                 kOO(1),ld) 
      endif
      if (pars%nvqp(1).gt.0. and. pars%vvlo(1).le.pars%vvhi(1) .and.
     $    pars%vvlo(1).gt.0) then
        call ga_access(pars%g_erivv(1),1,nri,pars%vvlo(1),pars%vvhi(1),
     $                 kVV(1),ld) 
      end if

C
C     Pointers to beta spin MO ERIs stored locally for a given MPI rank.
C
      if (ipol.gt.1 .and. pars%mynpoles(2).gt.0 ) then
        if (pars%ovlo(2).le.pars%ovhi(2) .and. pars%ovlo(2).gt.0) then
        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV(2),ld)
        endif
        if (pars%oolo(2).le.pars%oohi(2) .and. pars%oolo(2).gt.0) then
        call ga_access(pars%g_erioo(2),1,nri,pars%oolo(2),pars%oohi(2),
     $                 kOO(2),ld) 
        endif
        if (pars%nvqp(2).gt.0 .and. pars%vvlo(2).le.pars%vvhi(2) .and.
     $      pars%vvlo(2).gt.0) then
          call ga_access(pars%g_erivv(2),1,nri,pars%vvlo(2),
     $                   pars%vvhi(2),kVV(2),ld) 
        endif
      endif

      !
      ! Start loop over all imaginary grid points `ngl`
      !
      do igl=ngl,1,-1
c
#ifdef GWDEBUG
        if (main) write(luout,*) ' ## igl ',igl, ' diag ',diag
#endif

        !
        ! Load imaginary omega
        !
        w2 = x(igl)
        imw = dcmplx(0d0,dsqrt(w2))

        !
        ! Precompute square brackets of Equation (39).
        !
        do isp=1,ipol
!$omp      parallel do simd schedule(static)
           do ipole=1,pars%mynpoles(isp)
#ifdef GWCMPLX
             factor(ipole,isp) = sqrt(
     $                           0.5d0/(eia(ipole,isp) - imw + ieta) +
     $                           0.5d0/(eia(ipole,isp) + imw + ieta) )
#else
             factor(ipole,isp) = 
     $            dsqrt(eia(ipole,isp)/(eia(ipole,isp)**2 + w2))
#endif
          enddo
!$omp     end parallel do simd
        enddo

        !
        ! Build the Polarizability matrix \Pi (Equation (39)) and
        ! transform it into the dieletric matrix
        ! \epsilon = 1 - \Pi (in-line text just after Equation (39)).
        !
        ! In real arithmetic, the matrix is stored in RFP format.
        ! In complex arithmetic, the matrix is stored in full format.
        !
#ifdef GWDEBUG
        t0 = util_wallsec()
#endif

#ifdef GWCMPLX
        call gw_cdgw_buildpi_cmplx('w',dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $                         dcpl_mb(k_pi),factor,nri,nri,
     $                         pars%mynpoles,main,ipol,maxpoles)
#else
        call gw_cdgw_buildpi('w',dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $           dbl_mb(k_pi),factor,ldpi,nri,pars%mynpoles,main,ipol,
     $           maxpoles) 
#endif

#ifdef GWDEBUG
        if (main) write(luout,*) '     ## buildpi ',util_wallsec()-t0
#endif

        !
        ! Sum contributions from all MPI ranks
        !
#ifdef GWDEBUG
        t0 = util_wallsec()
#endif
        call ga_mask_sync(.true.,.false.)

#ifdef GWCMPLX
        call ga_zgop(1038,dcpl_mb(k_pi),nri**2,'+')
#else
        call ga_dgop(1038,dbl_mb(k_pi),sizepi,'+')
#endif

#ifdef GWDEBUG
        if (main) write(luout,*) '     ## gop ',util_wallsec()-t0
#endif

        !
        ! Copy dielectric matrix in case we need to compute the
        ! RPA correlation energy.
        !
        if (pars%rpa) call dcopy(sizepi,dbl_mb(k_pi),1,dbl_mb(k_vx),1)

        !
        ! The dielectric matrices along all the imaginary axes are 
        ! positive definite. The default behavior is to use a 
        ! Cholesky decomposition followed by a triangular solve
        ! in order to get \epsilon_{PQ}^{-1}(Q|mn) in Equation (41). 
        ! The user can request to perform a diagonalization 
        ! instead (using dft_scaleig) which could be faster than 
        ! the serial Cholesky decomposition.
        !
        ! In real arithmetic, the decomposition is performed directly in
        ! RFP format.
        !
#ifdef GWDEBUG
        t0 = util_wallsec()
#endif

#ifdef GWCMPLX
        call zpotrf('l',nri,dcpl_mb(k_pi),nri,info)
#else
        if (diag) then
          call gw_invdiag(dbl_mb(k_vx),ldpi,nri,ecrpa,.true.)
        else
          call dpftrf('n','l',nri,dbl_mb(k_pi),info)
        endif
#endif
#ifdef GWDEBUG
        if (main) write(luout,*) '     ## inverse ',util_wallsec()-t0
#endif

        !
        ! Compute the RPA correlation energy, if requested
        !
#ifndef GWCMPLX
        if (me.eq.0 .and. pars%rpa) then
          if (.not.diag) then
            call gw_cdgw_ecrpa(dbl_mb(k_pi),ldpi,nri,ecrpa)
          endif
          pars%ecrpa = pars%ecrpa + wgt(igl)*ecrpa
        endif
#endif


        ! So far we have either \epsilon^{-1} or the lower Cholesky
        ! factors of \epsilon. Here, we actually compute Equation (41)
        ! with the data at hand, first by gathering all ERIs into a
        ! contiguous array, and then by performing the appropriate
        ! contractions or solvers.
        !
        ! All the work perfomed here is done locally in each MPI rank.
        !
#ifdef GWDEBUG
        t0 = util_wallsec()
#endif

        do isp=1,ipol
          if (pars%nstates(isp).eq.0) cycle

          do imo=pars%llmo(isp),pars%llmo(isp)+pars%nqp(isp)-1
            jmo = imo - pars%llmo(isp) + 1
            
            !
            ! Deactive threaded BLAS since the work below involve
            ! only copies of small vectors.
            !
            call util_blas_set_num_threads(1)

            if (imo.le.pars%nocc(isp)) then
              
              !
              ! Gather all needed ERIs for a QP below the Fermi level.
              !
!$omp         parallel do private(kmo,lmo,idx,idy) schedule(static)
              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = (imo-1)*pars%nocc(isp)+lmo-1
                  idy = idy/nprocs*nri + kOO(isp)
                else
                  idy = ((imo-1)*pars%nvir(isp)+lmo-pars%nocc(isp)-1)
                  idy = idy/nprocs*nri + kOV(isp)
                endif
                call dcopy(nri,dbl_mb(idy),1,dbl_mb(idx),1)
              enddo
!$omp         end parallel do

            else

              !
              ! Gather all needed ERIs for a QP above the Fermi level.
              !
!$omp         parallel do private(kmo,lmo,idx,idy) schedule(static)
              do kmo=1,iwidx(0,jmo,isp)
                lmo = iwidx(kmo,jmo,isp)
                idx = tmp + (kmo-1)*nri
                if (lmo.le.pars%nocc(isp)) then
                  idy = ((lmo-1)*pars%nvir(isp)+imo-pars%nocc(isp)-1)
                  idy = idy/nprocs*nri + kOV(isp)
                else
                  idy = ((imo-pars%nocc(isp)-1)*pars%nvir(isp)+lmo-
     $                    pars%nocc(isp)-1)/nprocs*nri + kVV(isp)
                endif
                call dcopy(nri,dbl_mb(idy),1,dbl_mb(idx),1)
              enddo
!$omp         end parallel do
            endif

            !
            ! Reactivate threaded BLAS for the matrix-matrix operations
            ! below
            !
#ifdef USE_OPENMP
            call util_blas_set_num_threads(iMaxThreads)
#endif            

            !
            ! If diagonalization was requested, we directly transform
            ! the ERIs with the inverse dielectric matrix:
            !
            ! \epsilon_{PQ}^{-1} (Q|mn)
            !
            ! Note that the matrix is stored in RFP format and some more
            ! logic is needed.
            !

            if (diag) then
              call gw_rfpmm(dbl_mb(k_vx),dbl_mb(tmp),dbl_mb(k_sol),
     $                      ldpi,nri,iwidx(0,jmo,isp),.true.)
            else
            !
            ! If Cholesky decomposition was used, we first need to copy
            ! the ERIs into the solution arrays and then call the
            ! triangular solvers from LAPACK
            !
            ! Note that in real arithmetic, we compute 
            ! L_{PQ}^{-1} (Q|mn), where L is the lower Cholesky factor.
            ! In complex arithmetic, we compute
            ! \epsilon^{-1} (Q|mn) directly.
            !
#ifdef GWCMPLX
              call zlacp2('A',nri,iwidx(0,jmo,isp),dbl_mb(tmp),nri,
     $                   dcpl_mb(k_sol),nri)
              call zpotrs('l',nri,iwidx(0,jmo,isp),dcpl_mb(k_pi),nri,
     $                   dcpl_mb(k_sol),nri,info)
#else
              ! Copy the ERIs into the solution array
              call gw_copy(nri*iwidx(0,jmo,isp),dbl_mb(tmp),
     $                     dbl_mb(k_sol))

              ! Obtain \epsilon^{-1/2} (Q|mn)
              call dtfsm('n','l','l','n','n',nri,iwidx(0,jmo,isp),1d0,
     $                   dbl_mb(k_pi),dbl_mb(k_sol),nri)       
#endif
            endif

            !
            ! This is the last part to assembly the screened Coulomb
            ! matrix elements given in Equation (41):

            ! iW(mn,mn) = (mn|P)\epsilon^{-1}(Q|mn) - (mn|P)(P|mn)
            !
            ! for GWCMPLX and diag=.true.  `idy` will point to
            ! \epsilon^{-1} (Q|mn)
            !
            ! For Cholesky decomposition over reals, `idy` will point
            ! to \epsilon^{-1/2} (Q|mn)
            !
            ! `idx` points to (P|mn) in all cases
            !
!$omp       parallel do private(kmo,idx,idy) schedule(static)
            do kmo=1,iwidx(0,jmo,isp)
              idx = tmp + (kmo-1)*nri
              idy = k_sol + (kmo-1)*nri
              if (diag) then
                iw(kmo,igl,jmo,isp) = 
     $           gw_ddot(nri,dbl_mb(idx),dbl_mb(idy)) - 
     $           gw_ddot(nri,dbl_mb(idx),dbl_mb(idx))
              else
#ifdef GWCMPLX
                iw(kmo,igl,jmo,isp) = 
     $             gw_ddot_cmplx(nri,dbl_mb(idx),dcpl_mb(idy)) -
     $             gw_ddot(nri,dbl_mb(idx),dbl_mb(idx))
#else
                iw(kmo,igl,jmo,isp) = 
     $             gw_ddot(nri,dbl_mb(idy),dbl_mb(idy)) - 
     $             gw_ddot(nri,dbl_mb(idx),dbl_mb(idx))
#endif
              endif
            enddo
!$omp       end parallel do
          enddo
        enddo

#ifdef GWDEBUG
        if (main) write(luout,*) '     ## transf. ',util_wallsec()-t0
#endif

#ifndef GWCMPLX
        ! Compute the RPA correlation energy according to the
        ! Galitskii-Migdal formula.
        if (me.eq.0 .and. pars%rpa) then
          call gw_cdgw_ecgm(dbl_mb(k_vx),dbl_mb(k_pi),nri,ldpi,ecgm,
     $                      ecrpa,diag)
          pars%ecgm = pars%ecgm + wgt(igl)*ecgm
          pars%ecrpa = pars%ecrpa + wgt(igl)*ecrpa
        endif
#endif

      enddo

      !
      ! Free all local arrays
      !
      deallocate(factor,STAT = stat_alloc)
      if(stat_alloc.ne.0)
     &     call errquit(' gw_buildiw: deallocate failed for size',
     M     maxpoles*ipol, MA_ERR)
      if (.not.(ma_free_heap(l_sol)))
     $  call errquit(pname//'can''t free SOL heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_tmp)))
     $  call errquit(pname//'can''t free TMP heap',me,MA_ERR) 
      if (.not.(ma_free_heap(l_pi)))
     $  call errquit(pname//'can''t free PI heap',me,MA_ERR) 
#ifndef GWCMPLX
      if (pars%rpa) then
        if (.not.(ma_free_heap(l_vx)))
     $  call errquit(pname//'can''t free vX heap',me,MA_ERR) 
      endif
#endif

      pars%ecrpa = pars%ecrpa*0.5d0
      pars%ecgm  = pars%ecgm*0.5d0


      !
      ! Release modified GA arrays
      !
      if (pars%ovlo(1).le.pars%ovhi(1) .and.pars%ovlo(1).gt.0) then
        call ga_release_update(pars%g_eriov(1),1,nri,pars%ovlo(1),
     $                         pars%ovhi(1))
      endif
      if (pars%oolo(1).le.pars%oohi(1) .and.pars%oolo(1).gt.0) then
        call ga_release(pars%g_erioo(1),1,nri,pars%oolo(1),pars%oohi(1))
      endif
      if (pars%nvqp(1).gt.0. and. pars%vvlo(1).le.pars%vvhi(1) .and.
     $    pars%vvlo(1).gt.0) then
        call ga_release(pars%g_erivv(1),1,nri,pars%vvlo(1),pars%vvhi(1))
      end if

      if (ipol.gt.1 .and. pars%mynpoles(2).gt.0 ) then
        if (pars%ovlo(2).le.pars%ovhi(2) .and. pars%ovlo(2).gt.0) then
        call ga_release_update(pars%g_eriov(2),1,nri,pars%ovlo(2),
     $                         pars%ovhi(2))
        endif
        if (pars%oolo(2).le.pars%oohi(2) .and. pars%oolo(2).gt.0) then
        call ga_release(pars%g_erioo(2),1,nri,pars%oolo(2),pars%oohi(2))
        endif
        if (pars%nvqp(2).gt.0 .and. pars%vvlo(2).le.pars%vvhi(2) .and.
     $      pars%vvlo(2).gt.0) then
          call ga_release(pars%g_erivv(2),1,nri,pars%vvlo(2),
     $                   pars%vvhi(2))
        endif
      endif

      ! Go back to 1 thread
#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      
c
      end subroutine



      subroutine gw_copy(n,from,to)
      integer n
      double precision,dimension(n) :: from, to
      to(:) = from(:)
      end subroutine

      double precision function gw_ddot(n,a,b)
      integer n
      double precision,dimension(n) :: a,b
      gw_ddot = dot_product(a,b)
      end function

      double complex function gw_ddot_cmplx(n,a,b)
      integer n
      double precision,dimension(n) :: a
      double complex,dimension(n) :: b
      gw_ddot_cmplx = sum(a(:)*b(:))
      end function
