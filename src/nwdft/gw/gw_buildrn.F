#ifdef GWCMPLX
#define GWTYP double complex
#else
#define GWTYP double precision
#endif

      subroutine gw_buildrn(pars,Eig,Eia,omega,R,dR,imo,maxpoles,
     $                      nri,me,nmo,ipol,isp)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: imo, maxpoles, nri, me, nmo, ipol, isp
      double precision, intent(in) :: omega
      double precision, intent(in), dimension(nmo) :: Eig
      double precision, intent(in), dimension(maxpoles,ipol) :: Eia
      double precision, intent(out) :: R,dR

      character(len=12), parameter :: pname = 'gw_buildrn: '

      integer mynpoles, root2
      integer root, oolo, oohi, vvlo, vvhi, nrisize, ovlo, ovhi, nspace
      integer g_a, g_b, g_c, g_d, g_e, g_epsilon, k_eps, k_sol, nprocs
      integer jmo, ipole, ilo, ihi, jlo, jhi, info, ld, mynri, idx, jsp
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, k_ov, lscr
      integer g_pi, g_sol, addrPi, addrSol, llpi, off
      integer pihi, pilo, kOV(2), ga, first, second, step
      integer spacelo, addr, gb, owner, ints, nris
      double precision fac, arg, w, arg2, dfac, sgn
      double precision factor(maxpoles,ipol), vector(nri), vector2(nri)
      double precision vector3(nri)
      double precision dfactor(maxpoles,ipol)

      GWTYP,dimension(nri) ::   solvector

      logical main, ok, parallel_diag
      double precision, parameter :: shift = 0.00001d0
      double precision, parameter :: eta = 0.001d0
      double complex, parameter :: ieta = dcmplx(0d0,eta)
      double complex :: temp
      logical, external :: is_near
#ifdef USE_OPENMP      
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif


      main = me.eq.0
      nprocs = ga_nnodes()
      ints = ma_sizeof(mt_int,1,mt_byte)
      nris = ma_sizeof(mt_dbl,nri,mt_byte)


      if (.not.ma_alloc_get(mt_dbl,nri**2,'pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate pi matri',0,MA_ERR) 
      if (.not.ma_push_get(mt_int,nri,'ipiv',l_ipiv,k_ipiv))
     $  call errquit(pname//'can''t allocate ipiv',0,MA_ERR) 
      call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $            factor,-1,info)
      lscr = int(factor(1,1))
      if (.not.ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $  call errquit(pname//'can''t allocate scratch',0,MA_ERR) 

      do jsp=1,ipol
        call ga_access(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                 pars%ovhi(jsp),kOV(jsp),ld)   
      enddo
      if (ipol.eq.1) kOV(2) = 1

#if PARALLEL_DIAG
      parallel_diag = nri .gt. 2000
#else
      parallel_diag = .false.
#endif


      if (parallel_diag) then
        pilo = 1
        pihi = nri
#ifndef SCALAPACK
        ok = ga_create(mt_dbl,nri,nri,'Pi',nri,0,g_pi)
        if (.not.ok) call errquit(pname//'can''t create Pi',0,GA_ERR)
        call ga_distribution(g_pi,me,ilo,ihi,pilo,pihi)
        call ga_access(g_pi,1,nri,pilo,pihi,addrPi,ld)
        mynri = pihi - pilo + 1
        llpi = k_pi + (pilo-1)*nri
#endif
      endif

      R = 0d0
      dR = 0d0

      dfac = 4d0
      if (ipol.eq.1) dfac = 2d0*dfac

      gb = pars%g_sols

      if (imo.le.pars%nocc(isp)) then
        sgn = -1d0
        ga = pars%g_erioo(isp)
        spacelo = pars%oolo(isp)
        first = pars%nocc(isp)
        second = 1
        step = -1
        off = 0
        nspace = pars%nocc(isp)
      else
        sgn = 1d0
        ga = pars%g_erivv(isp)
        spacelo = pars%vvlo(isp)
        first = pars%nocc(isp) + 1
        second = nmo
        step = 1
        off = pars%nocc(isp)
        nspace = pars%nvir(isp)
      endif

      if (omega.lt.0d0) then
        sgn = -1d0
      else
        sgn = 1d0
      endif

       
      do jmo=first,second,step

        if (omega.lt.0d0) then
          if (eig(jmo) .lt. omega-shift) cycle
          if (eig(jmo) .gt. 0d0) cycle
        else
          if (eig(jmo) .gt. omega+shift) cycle
          if (eig(jmo) .lt. 0d0) cycle
        endif

        if (is_near(eig(jmo),omega,shift)) then
          fac = sgn*0.5d0
          arg = 0d0
          arg2 = 0d0
        else
          fac = sgn
          arg = eig(jmo) - omega
          arg2 = arg**2
        endif

!$omp   parallel do collapse(2) private(temp)
        do jsp=1,ipol
          do ipole=1,pars%mynpoles(jsp)
            temp = (Eia(ipole,jsp) - ieta)/
     $      ((Eia(ipole,jsp)**2-arg2)-2d0*ieta*(Eia(ipole,jsp)-sgn*arg))
            factor(ipole,jsp) = dble(temp)

            temp = temp*(arg+ieta*sgn)/
     $      ((Eia(ipole,jsp)**2-arg2)-2d0*ieta*(Eia(ipole,jsp)-sgn*arg))
            dfactor(ipole,jsp) = -dfac*dble(temp)
          enddo
        enddo
!$omp   end parallel do        

        call gw_buildpi('r',dbl_mb(kOV(1)),dbl_mb(kOV(2)),dbl_mb(k_pi),
     $                   factor,nri,nri,pars%mynpoles,main,ipol,
     $                   maxpoles)
        call ga_dgop(mt_dbl,dbl_mb(k_pi),nri**2,'+')

        root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
        idx = ((imo-off-1)*nspace+jmo-off-1)/nprocs + spacelo

        if (parallel_diag) then
          call gw_symmat(dbl_mb(k_pi),nri,pilo,pihi)
          call ga_brdcst(mt_int,idx,ints,root)
          call ga_get(ga,1,nri,idx,idx,vector,nri)
#ifdef SCALAPACK
          call dcopy(nri,vector,1,vector2,1)
          call gw_lu_solve(dbl_mb(k_pi),vector2,nri)
#else
          vector2(:) = 0d0
          call dcopy(mynri*nri,dbl_mb(llpi),1,dbl_mb(addrPi),1)
          call gw_invdiag(g_pi,g_pi,nri)
          call dgemv('t',nri,mynri,1d0,dbl_mb(addrPi),nri,vector,1,
     $                0d0,vector2(pilo),1)
          call ga_dgop(mt_dbl,vector2,nri,'+')
#endif
        else
          if (me.eq.root) then
            call ga_get(ga,1,nri,idx,idx,vector,nri)
            call dcopy(nri,vector,1,vector2,1)
            call dsytrf('l',nri,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $                   dbl_mb(k_scr),lscr,info)
            call dsytrs('l',nri,1,dbl_mb(k_pi),nri,int_mb(k_ipiv),
     $                   vector2,nri,info)
          end if
          call ga_brdcst(mt_dbl,vector2,nris,root)
        endif

        if (.not.ga_locate(gb,1,jmo,owner))
     $    call errquit('gw_buildrn: cannot locate owner',0,GA_ERR)
#ifdef GWCMPLX        
        solvector = dcmplx(vector2,0d0)
#else
        solvector = vector2
#endif
        if (owner.eq.me) call ga_put(gb,1,nri,jmo,jmo,solvector,nri)

        if (arg.ne.0d0) then
          do jsp=1,ipol
            addr = kOV(jsp)
            call dgemv('t',nri,pars%mynpoles(jsp),1d0,dbl_mb(addr),nri,
     $                  vector2,1,0d0,factor(1,jsp),1)

!$omp       parallel do reduction(+:dR)
            do ipole=1,pars%mynpoles(jsp)
              dR = dR + dfactor(ipole,jsp)*factor(ipole,jsp)**2
            enddo
!$omp       end parallel do

          enddo
        endif

        if (root.ne.me) cycle

        w = dot_product(vector,vector2) - dot_product(vector,vector)
        R = R + fac*w

        if (abs(fac).lt.1d0) dR = dR + w

      enddo

      call ga_dgop(mt_dbl,R,1,'+')
      call ga_dgop(mt_dbl,dR,1,'+')


      if (.not.ma_chop_stack(l_ipiv))
     $   call errquit(pname//'can''t chop stack',0,MA_ERR)
      if (.not.ma_free_heap(l_pi))
     $   call errquit(pname//'can''t free heap',0,MA_ERR)
c     

      if (parallel_diag) then
#ifndef SCALAPACK
        ok = ga_destroy(g_pi)
        if (.not.ok) call errquit(pname//'can''t destroy GAs',0,GA_ERR)
#endif
      endif

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      return
      end subroutine
