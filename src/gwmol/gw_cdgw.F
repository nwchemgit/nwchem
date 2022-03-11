#ifdef GWCMPLX
#define GWMT mt_dcpl
#define GWMB dcpl_mb
#else
#define GWMT mt_dbl
#define GWMB dbl_mb
#endif

      subroutine gw_cdgw(pars)
c
c     Contour Deformation GW
c
c     D. Mejia-Rodriguez, A. Kunitsa, E. Apra, N. Govind
c     J. Chem. Theory Comput. 17, 7504 (2021).
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"      
#include "gw.fh"

      type(gw_params_t) :: pars
      character(len=9), parameter :: pname = 'gw_cdgw: '

      integer,dimension(2) :: nocc,nvir,npoles,nstates,mynpoles

      integer nri, nmo
      integer nri_me, istart, iend, ld
      logical main, oactive

      integer k_eriij,k_eriia,kSigma,kSigmax,kEref,kEnew,kEmod,kValues
      integer l_eriij,l_eriia,lSigma,lSigmax,lEref,lEnew,lEmod,lValues
      integer kErrors,lErrors,ndiis, offset, isp
      integer k_omega,kEia
      integer l_omega,lEia
      integer g_epsilon(pars%ngl)
      integer igl,ngl,nqp

      integer nprocs, maxqp, lsize, maxpoles, temp
      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4, mynmo
      integer oolo,oohi,jlo,jhi,vvlo,vvhi,kiW,liW,kiWidx,liWidx

      double precision efermi(ipol),efermi_new(ipol)
      double precision timer

      logical converged,ok
      double precision, external :: gw_efermi
c
c     -----------------------------------------------------------------
c
      ngl = pars%ngl
      nri = nbf_cd
      nmo = pars%nmo
     

      do i=1,ipol
        nocc(i) = pars%nocc(i)
        nvir(i) = pars%nvir(i)
        npoles(i) = nocc(i)*nvir(i)
        nstates(i) = pars%nstates(i)
        mynpoles(i) = pars%mynpoles(i)
      enddo

      main = pars%me .eq. 0
      nprocs = ga_nnodes()
      maxqp = maxval(pars%nqp(1:ipol))
      maxpoles = maxval(pars%mynpoles(1:ipol))
c
c     Obtain Vxc
c
      if (main) write(luout,9000)
 9000 format(10x,'Computing Vxc')

      timer = util_wallsec()

      lSize = ipol*nmo
      if (.not.ma_push_get(mt_dbl,lSize,'Sigma',lSigma,kSigma))
     $    call errquit(pname//'can''t create Sigma',0,MA_ERR)
      call dfill(lSize,0d0,dbl_mb(kSigma),1)

      call gw_vxc(dbl_mb(kSigma),nmo,pars)

      pars%timer%vxc = pars%timer%vxc + util_wallsec() - timer
c
c     Obtain Sigma_x
c
      if (main) write(luout,9010)
 9010 format(10x,'Computing Sigma_x')

      timer = util_wallsec()

      lSize = ipol*nmo
      if (.not.ma_push_get(mt_dbl,lSize,'Sigma_x',lSigmax,kSigmax))
     $  call errquit(pname//'cannot allocate Sigma_x',0,MA_ERR)

      call gw_sigmax(pars,dbl_mb(kSigmax),nri,pars%nmo)

      pars%timer%vmo = pars%timer%vmo + util_wallsec() - timer
c
c     Allocations
c
      mynmo = pars%nmo/nprocs + maxval(pars%nocc(1:pars%ipol)) + 1

      lSize = ipol*mynmo*maxqp*ngl
      if (.not.ma_alloc_get( GWMT ,lSize,'W(igl)',liW,kiW))
     $  call errquit(pname//'can''t create iw array',0,MA_ERR)

      lSize = ipol*(mynmo+1)*maxqp
      if (.not.ma_alloc_get(mt_int,lSize,'iW idx',liWidx,kiWidx))
     $  call errquit(pname//'can''t create iw index',0,MA_ERR)

      lSize = ipol*pars%nmo
      if (.not.ma_push_get(mt_dbl,lSize,'evals ref',lEref,kEref))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)

      if (.not.ma_push_get(mt_dbl,lSize,'evals mod',lEmod,kEmod))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)

      if (.not.ma_push_get(mt_dbl,lSize,'evals new',lEnew,kEnew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)
      call dfill(lSize,0d0,dbl_mb(kEnew),1)

      if (pars%evgw .or. pars%evgw0) then
        lSize = lSize*max(pars%eviter,1)
        if (.not.ma_push_get(mt_dbl,lSize,'values',lValues,kValues))
     $    call errquit(pname//'failed to allocate values',0,MA_ERR)
        lSize = lSize*max(pars%eviter,1)
        if (.not.ma_push_get(mt_dbl,lSize,'errors',lErrors,kErrors))
     $    call errquit(pname//'failed to allocate errors',0,MA_ERR)
      endif

      lSize = ipol*maxpoles
      if (.not.ma_push_get(mt_dbl,lSize,'eval_ia',lEia,kEia))
     $  call errquit(pname//'cannot allocate eval_ia',0,MA_ERR)
c
c     Get fermi energy
c
      do i=1,ipol
        lSize = pars%k_evals + (i-1)*nbf_ao
        efermi(i) = gw_efermi(dbl_mb(lSize),pars%nocc(i))
      enddo
c
c     Shift eigenvalues
c
      call gw_shift(dbl_mb(pars%k_evals),dbl_mb(kEref),efermi,
     $              nbf_ao,pars%nmo,pars%ipol)

      call dcopy(nmo*ipol,dbl_mb(kEref),1,dbl_mb(kEmod),1)
c
      if(.not.(ga_create( GWMT ,nri,pars%nmo,'sols',nri,0,
     $    pars%g_sols)))
     $  call errquit('gw_cdgw: cannot create GA',0,GA_ERR) 
c
      eviter = 0
      ndiis = 0
c
  100 continue ! evGW calculations
c
      if (pars%evgw.and.main) then
        if (eviter.lt.10) then
          write(luout,9020) eviter,eviter
        else
          write(luout,9021) eviter,eviter
        endif
      endif
 9020 format(/,2x,"*** G",I1,"W",I1)     
 9021 format(/,2x,"*** G",I2,"W",I2)     

      if (eviter.eq.0 .and. pars%ncap .and. pars%evgw) then
        do i=1,ipol
          do j=0,pars%nmo-1
            if (j.lt.pars%nocc(i)) then
              dbl_mb(kEmod+(i-1)*nmo+j) = 
     $        dbl_mb(kEmod+(i-1)*nmo+j) + pars%vxddp(i)
            else
              dbl_mb(kEmod+(i-1)*nmo+j) = 
     $        dbl_mb(kEmod+(i-1)*nmo+j) + pars%vxddm(i)
            endif
          enddo
        enddo
      endif
c
c     Obtain eigenvalue difference
c
      do i=1,ipol
        lSize = (i-1)*maxpoles
        call gw_get_eia(dbl_mb(kEmod+(i-1)*pars%nmo),dbl_mb(kEia+lSize),
     $                  pars%nocc(i),pars%nvir(i))
      enddo

c
c     Experimental bit, might be useful for BSE calculations
c
      if (eviter.eq.0 .and. pars%ncap .and. (.not.pars%evgw)) then
        do i=1,ipol
          do j=0,pars%nmo-1
            if (j.lt.pars%nocc(i)) then
              dbl_mb(kEmod+(i-1)*nmo+j) = 
     $        dbl_mb(kEmod+(i-1)*nmo+j) + pars%vxddp(i)
            else
              dbl_mb(kEmod+(i-1)*nmo+j) = 
     $        dbl_mb(kEmod+(i-1)*nmo+j) + pars%vxddm(i)
            endif
          enddo
        enddo
      endif
c
c      if (eviter.eq.0 .and. pars%ncap) then
c        call dcopy(nmo*ipol,dbl_mb(kEmod),1,dbl_mb(kEref),1)
c      endif
c
c     Obtain Screened Coulomb matrices on imaginary axis
c
      if (main) write(luout,9030)
 9030 format(10x,'Computing W(iw) on the imaginary grid')
c
      timer = util_wallsec()

      do i=1,ipol
        lSize = (i-1)*(mynmo+1)*maxqp
        call gw_get_iwidx(int_mb(kiWidx+lSize),pars%nmo,pars%nocc(i),
     $                    pars%nvir(i),mynmo,pars%nqp(i),pars%llmo(i))
      enddo

      call gw_cdgw_buildiw(pars,dbl_mb(kEia),dbl_mb(pars%k_glx),
     $          dbl_mb(pars%k_glw), GWMB(kiW), int_mb(kiWidx),
     $          pars%nmo,nri,pars%ngl,mynmo,maxpoles,maxqp,pars%ipol)


      pars%timer%iw = pars%timer%iw + util_wallsec() - timer
c
  200 continue ! evGW0
c
      if (pars%evgw0.and.main) then
        if (eviter.lt.10) then
          write(luout,9040) eviter
        else
          write(luout,9041) eviter
        endif
      endif
 9040 format(/,2x,"*** G",I1,"W0")     
 9041 format(/,2x,"*** G",I2,"W0")     


      if (main) write(luout,9050)
 9050 format(10x,'Solving quasiparticle equations')

      if (pars%graph) then
        call gw_cdgw_graph(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEref),dbl_mb(kEmod),dbl_mb(kEnew),
     $                      dbl_mb(kEia),
     $                      GWMB(kiW), int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles,eviter)
      elseif (pars%newton) then
        call gw_cdgw_newton(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEref),dbl_mb(kEmod),dbl_mb(kEnew),
     $                      dbl_mb(kEia),
     $                      GWMB(kiW), int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles,eviter)
      endif

      if (pars%evgw0.or.pars%evgw) then
c
c       Apply shift to the eigenvalues that were not modified by GW
c
        call gw_scissor(dbl_mb(kEnew),dbl_mb(kEmod),pars%nmo,pars%nocc,
     $                  pars%nvir,pars%noqp,pars%nvqp,pars%ipol,pars%me)

        call ga_brdcst(1038,dbl_mb(kEnew),
     $       ma_sizeof(mt_dbl,nmo*pars%ipol,mt_byte),0)
        call ga_brdcst(1038,dbl_mb(kSigma),
     $       ma_sizeof(mt_dbl,nmo*pars%ipol,mt_byte),0)
        !
        ! Update Fermi energy
        !
        do i=1,ipol
          lSize = kEnew + (i-1)*nbf_ao
          efermi_new(i) = gw_efermi(dbl_mb(lSize),pars%nocc(i))
        enddo
        call gw_shift(dbl_mb(kEnew),dbl_mb(kemod),efermi_new,
     $                nbf_ao,pars%nmo,pars%ipol)
        efermi(:) = efermi(:) + efermi_new(:)
        call dcopy(pars%ipol*nmo,dbl_mb(kemod),1,dbl_mb(kEnew),1)

        eviter = eviter + 1
        
        if (eviter.lt.pars%eviter) then
          if (pars%evgw0) goto 200
          if (pars%evgw)  goto 100
        endif  
      endif

      if (.not.ga_destroy(pars%g_sols))
     $  call errquit('gw_cdgw_newton: cannot destroy GA',0,GA_ERR) 

      if (.not.(ma_free_heap(liW).and.ma_free_heap(liWidx)))
     $  call errquit(pname//'cannot destroy iW arrays',0,MA_ERR)

c
      end subroutine

