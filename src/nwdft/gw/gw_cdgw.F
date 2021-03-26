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

      integer k_eriij,k_eriia,kSigma,kSigmax,kEref,kEnew,kEmod
      integer l_eriij,l_eriia,lSigma,lSigmax,lEref,lEnew,lEmod
      integer k_omega,kEia
      integer l_omega,lEia
      integer g_epsilon(pars%ngl)
      integer igl,ngl,nqp

      integer nprocs, maxqp, lsize, maxpoles, temp
      integer eviter, i, j, ibf, ilo, ihi, lo2, ipole, iqp, lo4, mynmo
      integer oolo,oohi,jlo,jhi,vvlo,vvhi,kiW,liW,kiWidx,liWidx

      double precision efermi(ipol)
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
      if (.not.ma_alloc_get(mt_dbl,lSize,'W(igl)',liW,kiW))
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
#ifdef GWCMPLX
      if(.not.(ga_create(mt_dcpl,nri,pars%nmo,'sols',nri,0,
     $    pars%g_sols)))
#else
      if(.not.(ga_create(mt_dbl,nri,pars%nmo,'sols',nri,0,
     $    pars%g_sols)))
#endif
     $  call errquit('gw_cdgw: cannot create GA',0,GA_ERR) 
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
        call gw_get_eia(dbl_mb(kEmod+(i-1)*pars%nmo),dbl_mb(kEia+lSize),
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

      call gw_buildiw(pars,dbl_mb(kEia),dbl_mb(pars%k_glx),
     $          dbl_mb(pars%k_glw),dbl_mb(kiW),int_mb(kiWidx),
     $          pars%nmo,nri,pars%ngl,mynmo,maxpoles,maxqp,pars%ipol)

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
     $                      dbl_mb(kEref),dbl_mb(kEmod),dbl_mb(kEnew),
     $                      dbl_mb(kEia),
     $                      dbl_mb(kiW),int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles,eviter)
      elseif (pars%newton) then
        call gw_cdgw_newton(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEref),dbl_mb(kEmod),dbl_mb(kEnew),
     $                      dbl_mb(kEia),
     $                      dbl_mb(kiW),int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles,eviter)
      else
        call gw_cdgw_linear(pars,dbl_mb(kSigmax),dbl_mb(kSigma),
     $                      dbl_mb(kEref),dbl_mb(kEmod),dbl_mb(kEnew),
     $                      dbl_mb(kEia),
     $                      dbl_mb(kiW),int_mb(kiWidx),efermi,nri,
     $                      mynmo,maxqp,maxpoles)
      endif

      if (pars%evgw0.or.pars%evgw) then
        call gw_scissor(dbl_mb(kEnew),dbl_mb(kEmod),pars%nmo,pars%nocc,
     $                  pars%nvir,pars%noqp,pars%nvqp,pars%ipol,pars%me)

        call ga_brdcst(mt_dbl,dbl_mb(kEnew),
     $       ma_sizeof(mt_dbl,nmo*pars%ipol,mt_byte),0)
        call ga_brdcst(mt_dbl,dbl_mb(kSigma),
     $       ma_sizeof(mt_dbl,nmo*pars%ipol,mt_byte),0)
        call dcopy(nmo*pars%ipol,dbl_mb(kEnew),1,dbl_mb(kEmod),1)
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

      pars%ex = 0d0

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
            do kmo=1,pars%nocc(ipol)
              idx = (imo-1)*pars%nocc(ipol)+kmo-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOO
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
c
          else
c
            do kmo=1,pars%nocc(ipol)
              idx = (kmo-1)*pars%nvir(ipol)+imo-pars%nocc(ipol)-1
              if (pars%me.ne.mod(idx,pars%nprocs)) cycle
              idx = (idx/pars%nprocs)*nri + addrOV
              temp = temp - ddot(nri,dbl_mb(idx),1,dbl_mb(idx),1)
            enddo
c
          endif
          sigmax(imo,ipol) = temp

        enddo

      enddo

      call ga_dgop(mt_dbl,sigmax,nmo*pars%ipol,'+')

      do ipol=1,pars%ipol
        pars%ex = pars%ex+sum(sigmax(1:pars%nocc(ipol),ipol))/pars%ipol
      enddo

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
      integer ipole,imo,jmo,me,nprocs,jpole

      me = ga_nodeid()
      nprocs = ga_nnodes()
c
      ipole = 0
      jpole = 0
      do imo=1,nocc
        do jmo=nocc+1,nocc+nvir
          jpole = jpole + 1
          if (me.ne.mod(jpole-1,nprocs)) cycle
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

      subroutine gw_scissor(enew,eold,nmo,nocc,nvir,noqp,nvqp,ipol,me)
      implicit none
      integer,intent(in) :: ipol, nmo, me
      integer,intent(in),dimension(ipol) :: nocc,nvir,noqp,nvqp
      double precision,intent(in),dimension(nmo,ipol) :: eold
      double precision,intent(out),dimension(nmo,ipol) :: enew

      integer isp,imo
      double precision shift

      do isp=1,ipol

        ! shift occupied states as HOMO
        if (noqp(isp).lt.nocc(isp) .and.noqp(isp).gt.0) then
          shift = 0d0
          do imo=nocc(isp),nocc(isp)-noqp(isp)+1,-1
            shift = shift + enew(nocc(isp),isp) - eold(nocc(isp),isp)
          enddo
          shift = shift/noqp(isp)
          if (me.eq.0) write(*,9000) shift*27.2114d0
          do imo=nocc(isp)-noqp(isp),1,-1
            enew(imo,isp) = eold(imo,isp) + shift
          enddo
        endif

        ! shift unoccupied states as LUMO
        if (nvqp(isp).lt.nvir(isp) .and. nvqp(isp).gt.0) then
          shift = 0d0
          do imo=nocc(isp)+1,nocc(isp)+nvqp(isp)
            shift = shift + enew(nocc(isp)+1,isp)-eold(nocc(isp)+1,isp)
          enddo
          shift = shift/nvqp(isp)
          if (me.eq.0) write(*,9010) shift*27.2114d0
          do imo=nocc(isp)+nvqp(isp)+1,nmo
            enew(imo,isp) = eold(imo,isp) + shift
          enddo
        endif
      enddo

 9000 format(/,3X,'Applying ',F8.4,' eV scissor-shift to rest of ',
     $         'particle states',/)
 9010 format(/,3X,'Applying ',F8.4,' eV scissor-shift to rest of ',
     $         'hole states',/)
      end subroutine
