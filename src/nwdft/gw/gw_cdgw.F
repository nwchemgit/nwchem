      subroutine gw_cdgw(pars)
c
c     Contour Deformation GW
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

      integer k_eriij,k_eriia,kSigma,kSigmax,kEold,kEnew
      integer l_eriij,l_eriia,lSigma,lSigmax,lEold,lEnew
      integer k_omega,kEia
      integer l_omega,lEia
      integer g_epsilon(pars%ngl)
      integer igl,ngl,nqp

      integer nprocs, maxqp, lsize, maxpoles
      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4, mynmo
      integer oolo,oohi,jlo,jhi,vvlo,vvhi,kiW,liW,kiWidx,liWidx

      double precision efermi(ipol)
      double precision timer


      logical converged,ok

      double precision, parameter :: fdstep = 0.002d0
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
      if (.not.ma_alloc_get(mt_dbl,lSize,'W(igl)',liW,kiW))
     $  call errquit(pname//'can''t create iw array',0,MA_ERR)

      lSize = ipol*(mynmo+1)*maxqp
      if (.not.ma_alloc_get(mt_int,lSize,'iW idx',liWidx,kiWidx))
     $  call errquit(pname//'can''t create iw index',0,MA_ERR)

      lSize = ipol*pars%nmo
      if (.not.ma_push_get(mt_dbl,lSize,'evals old',lEold,kEold))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)

      if (.not.ma_push_get(mt_dbl,lSize,'evals new',lEnew,kEnew))
     $  call errquit(pname//'failed to allocate old evals',0,MA_ERR)

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
      call gw_shift(dbl_mb(pars%k_evals),dbl_mb(kEold),efermi,
     $              nbf_ao,pars%nmo,pars%ipol)
c
      eviter = 0
c
  100 continue ! evGW calculations
c
      if (pars%evgw.and.main) then
        write(luout,9020) eviter,eviter
      endif
 9020 format(/,2x,"*** G",I1,"W",I1)     
c
c     Obtain eigenvalue difference
c
      do i=1,ipol
        lSize = (i-1)*maxpoles
        call gw_get_eia(dbl_mb(kEold+(i-1)*pars%nmo),dbl_mb(kEia+lSize),
     $                  pars%nocc(i),pars%nvir(i))
      enddo
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

      if(.not.(ga_create(mt_dbl,nri,pars%nmo,'sols',nri,0,pars%g_sols)))
     $  call errquit('gw_cdgw_newton: cannot create GA',0,GA_ERR) 

      call gw_buildiw(pars,dbl_mb(kEia),dbl_mb(pars%k_glx),dbl_mb(kiW),
     $    int_mb(kiWidx),pars%nmo,nri,pars%ngl,mynmo,maxpoles,maxqp,
     $    pars%ipol)

c      call gw_buildiw_minres(pars,dbl_mb(kEia),dbl_mb(pars%k_glx),
c     $     dbl_mb(kiW),
c     $    int_mb(kiWidx),pars%nmo,nri,pars%ngl,mynmo,maxpoles,maxqp,
c     $    pars%ipol)

      pars%timer%iw = pars%timer%iw + util_wallsec() - timer
c
  200 continue ! evGW0
c
      if (pars%evgw0.and.main) then
        write(luout,9040) eviter
      endif
 9040 format(/,2x,"*** G",I1,"W0")     

      if (main) write(luout,9050)
 9050 format(10x,'Solving quasiparticle equations')

      if (pars%graph) then
        call gw_cdgw_graph(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                          dbl_mb(kEold),dbl_mb(kEnew),
     $                          dbl_mb(kEia),dbl_mb(kiW),int_mb(kiWidx),
     $                          efermi,nri,mynmo,maxqp,maxpoles)
      elseif (pars%newton) then
        call gw_cdgw_newton(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEold),dbl_mb(kEnew),dbl_mb(kEia),
     $                      dbl_mb(kiW),int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles,eviter)
      else
        call gw_cdgw_linear(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEold),dbl_mb(kEnew),dbl_mb(kEia),
     $                      dbl_mb(kiW),int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles)
      endif

      if (pars%evgw0.or.pars%evgw) then
        call ga_brdcst(mt_dbl,dbl_mb(kEnew),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call ga_brdcst(mt_dbl,dbl_mb(kSigma),
     $       ma_sizeof(mt_dbl,nmo,mt_byte),0)
        call dcopy(nmo,dbl_mb(kEnew),1,dbl_mb(kEold),1)
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

c     *****************************************************************

      logical function is_near(test,targ,tol)
      implicit none
      double precision, intent(in) :: test, targ, tol

      if (abs(test-targ).le.tol) then
        is_near = .true.
      else
        is_near = .false.
      endif

      return
      end function

c     *****************************************************************

      subroutine gw_sigmax(pars,sigmax,nri,nmo)
      implicit none
#include "mafdecls.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: nri, nmo
      double precision,intent(out),dimension(nmo,*) :: sigmax

      double precision :: temp
      double precision,external :: ddot

      integer idx, idum, isp
      integer imo, jmo, kmo, ipol
      integer addrOO, addrOV, addrVV
      integer,external :: omp_get_num_threads,omp_get_thread_num

      do ipol=1,pars%ipol
        if (pars%nqp(ipol).eq.0) cycle

        call ga_access(pars%g_erioo(ipol),1,nri,pars%oolo(ipol),
     $                 pars%oohi(ipol),addrOO,idum)
        call ga_access(pars%g_eriov(ipol),1,nri,pars%ovlo(ipol),
     $                 pars%ovhi(ipol),addrOV,idum)


        do imo=pars%llmo(ipol),pars%llmo(ipol)+pars%nqp(ipol)-1
          jmo = imo - pars%llmo(ipol) + 1
          temp = 0d0

          if (imo.le.pars%nocc(ipol)) then
c
!$omp       parallel do private(idx) reduction(+:temp)
            do kmo=1,pars%nocc(ipol)
              idx = (imo-1)*pars%nocc(ipol)+kmo-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOO
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
!$omp       end parallel do
c
          else
c
!$omp       parallel do private(idx) reduction(+:temp)
            do kmo=1,pars%nocc(ipol)
              idx = (kmo-1)*pars%nvir(ipol)+imo-pars%nocc(ipol)-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOV
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
!$omp       end parallel do
c
          endif
          sigmax(imo,ipol) = temp

        enddo

      enddo

      call dgop(mt_dbl,sigmax,nmo*pars%ipol,'+')

      end subroutine

c     *****************************************************************

      subroutine gw_buildpi(which,ovup,ovdw,pi,factor,ldpi,nri,npoles,
     $                      main,ipol,maxpoles)
      implicit none
#include "mafdecls.fh"      
      logical,intent(in) :: main
      integer,intent(in) :: ldpi, nri, maxpoles, ipol
      character(len=1),intent(in) :: which
      integer,intent(in),dimension(ipol) :: npoles
      double precision,intent(in),dimension(maxpoles,ipol) :: factor
      double precision,intent(inout),dimension(nri,*) :: ovup,ovdw
      double precision,intent(out),dimension(ldpi,*) :: pi

      double precision, dimension(nri) :: vector, scaling, scal
      integer ipole,iri,ma_dbl,l_scr,k_scr,mynpoles,llpole,mymax
      integer ulpole,nbatch,ibatch
      logical ok

      scal = 2d0
      if (ipol.eq.1) scal = 2d0*scal

      if (which.eq.'w') then

!$omp   parallel do
        do ipole=1,npoles(1)
          ovup(1:nri,ipole) = ovup(1:nri,ipole)*factor(ipole,1)
        enddo
!$omp   end parallel do

        call dsfrk('n','l','n',nri,npoles(1),scal,ovup,nri,0d0,pi)

!$omp   parallel do
        do ipole=1,npoles(1)
          ovup(1:nri,ipole) = ovup(1:nri,ipole)/factor(ipole,1)
        enddo
!$omp   end parallel do

        if (ipol.gt.1) then

!$omp     parallel do
          do ipole=1,npoles(2)
            ovdw(1:nri,ipole) = ovdw(1:nri,ipole)*factor(ipole,2)
          enddo
!$omp     end parallel do

          call dsfrk('n','l','n',nri,npoles(2),scal,ovdw,nri,1d0,pi)

!$omp     parallel do
          do ipole=1,npoles(2)
            ovdw(1:nri,ipole) = ovdw(1:nri,ipole)/factor(ipole,2)
          enddo
!$omp     end parallel do

        endif

        if (main) then
          if ( mod(nri,2) .eq. 0 ) then
!$omp       parallel do
            do iri=1,nri/2
              pi(iri,iri) = pi(iri,iri) + 1d0
              pi(iri+1,iri) = pi(iri+1,iri) + 1d0
            enddo
!$omp       end parallel do
          else
            pi(1,1) = pi(1,1) + 1d0
!$omp       parallel do
            do iri=2,nri/2+1
              pi(iri-1,iri) = pi(iri-1,iri) + 1d0
              pi(iri,iri) = pi(iri,iri) + 1d0
            enddo
!$omp       end parallel do
          endif
        endif

      else if (which.eq.'r') then

        scal = scal/2d0
        ma_dbl = ma_inquire_avail(mt_dbl)/nri
        nbatch = npoles(1)/ma_dbl + 1
        mymax = min(npoles(1)/nbatch + 1,npoles(1))
        ok = ma_alloc_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

        call dcopy(nri*mymax,ovup,1,dbl_mb(k_scr),1)
        do ipole=1,mymax
          call dscal(nri,factor(ipole,1),dbl_mb(k_scr+(ipole-1)*nri),1)
        enddo
        call dsyr2k('l','n',nri,mymax,scal,ovup,nri,dbl_mb(k_scr),nri,
     $               0d0,pi,nri)

        do ibatch=2,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(1))
          mynpoles = ulpole - llpole + 1
          call dcopy(nri*mynpoles,ovup(1,llpole),1,dbl_mb(k_scr),1)
          do ipole=1,mynpoles
            call dscal(nri,factor(llpole+ipole-1,1),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call dsyr2k('l','n',nri,mynpoles,scal,ovup(1,llpole),nri,
     $                 dbl_mb(k_scr),nri,1d0,pi,nri)     
        enddo

        ok = ma_free_heap(l_scr)

        if (ipol.gt.1) then
          nbatch = npoles(2)/ma_dbl + 1
          mymax = min(npoles(2)/nbatch + 1,npoles(2))
          ok = ma_alloc_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

          call dcopy(nri*mymax,ovdw,1,dbl_mb(k_scr),1)
          do ipole=1,mymax
            call dscal(nri,factor(ipole,2),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call dsyr2k('l','n',nri,mymax,scal,ovdw,nri,dbl_mb(k_scr),nri,
     $                 1d0,pi,nri)

          do ibatch=2,nbatch
            llpole = (ibatch-1)*mymax + 1
            ulpole = min(ibatch*mymax,npoles(2))
            mynpoles = ulpole - llpole + 1
            call dcopy(nri*mynpoles,ovdw(1,llpole),1,dbl_mb(k_scr),1)
            do ipole=1,mynpoles
              call dscal(nri,factor(llpole+ipole-1,2),
     $                   dbl_mb(k_scr+(ipole-1)*nri),1)
            enddo
            call dsyr2k('l','n',nri,mynpoles,scal,ovdw(1,llpole),nri,
     $                   dbl_mb(k_scr),nri,1d0,pi,nri)     
          enddo

          ok = ma_free_heap(l_scr)

        endif

        if (main) then
!$omp     parallel do
          do iri=1,nri
            pi(iri,iri) = pi(iri,iri) + 1d0
          enddo
!$omp     end parallel do
        endif

      endif

      end subroutine

c     *****************************************************************

      double precision function gw_efermi(eig,nocc)
      implicit none
      integer,intent(in) :: nocc
      double precision,intent(in),dimension(*) :: eig
      gw_efermi = (eig(nocc) + eig(nocc+1))/2d0
      return
      end function

      subroutine gw_shift(eig,shifted,efermi,nbf,nmo,npol)
      implicit none
      integer,intent(in) :: nbf, nmo, npol
      double precision,intent(in) :: efermi(npol),eig(nbf,npol)
      double precision,intent(out) :: shifted(nmo,npol)
      integer imo,ipol
!$omp parallel do collapse(2)      
      do ipol=1,npol
        do imo=1,nmo
          shifted(imo,ipol) = eig(imo,ipol) - efermi(ipol)
        enddo
      enddo
!$omp end parallel do      
      return
      end subroutine

      subroutine gw_get_eia(eig,eia,nocc,nvir)
      implicit none
#include "global.fh"      
      integer,intent(in) :: nocc,nvir
      double precision,intent(in),dimension(nocc+nvir) :: eig
      double precision,intent(out),dimension(*) :: eia
      integer ipole,imo,jmo,me,nprocs

      me = ga_nodeid()
      nprocs = ga_nnodes()
c
      ipole = 0
      do imo=1,nocc
        do jmo=nocc+1,nocc+nvir
          if (me.ne.mod((imo-1)*nvir+jmo-nocc-1,nprocs)) cycle
          ipole = ipole + 1
          eia(ipole) = eig(jmo) - eig(imo)
        enddo
      enddo
c      
      return
      end subroutine

      subroutine gw_get_iwidx(idx,nmo,nocc,nvir,mynmo,nqp,llmo)
      implicit none
      integer,intent(in) :: nmo, nocc, nvir, mynmo, nqp, llmo
      integer,intent(out),dimension(0:mynmo,*) :: idx

      integer me, imo, jmo, kmo, lmo, nprocs
      integer,external :: ga_nodeid, ga_nnodes

      me = ga_nodeid()
      nprocs = ga_nnodes()

      do imo=llmo,llmo+nqp-1
        jmo = imo - llmo + 1
        lmo = 0

        if (imo.le.nocc) then
          do kmo=1,nocc
            if (me.ne.mod((imo-1)*nocc+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        else
          do kmo=1,nocc
            if (me.ne.mod((kmo-1)*nvir+imo-nocc-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo
          enddo
          do kmo=1,nvir
            if (me.ne.mod((imo-nocc-1)*nvir+kmo-1,nprocs)) cycle
            lmo = lmo + 1
            idx(lmo,jmo) = kmo+nocc
          enddo
        endif
        idx(0,jmo) = lmo
      enddo
      return
      end subroutine

      subroutine gw_symmat(matrix,n,lo,hi)
      implicit none
      integer,intent(in) :: n, lo, hi
      double precision,intent(inout),dimension(n,n) :: matrix
      integer i,ilo

      do i=1,n-1
        matrix(i,i+1:) = matrix(i+1:,i)
      enddo

      return
      end subroutine
