#ifdef GWCMPLX
#define ISOGW complex(kind=C_DOUBLE_COMPLEX)
#define GWTYP double complex
#define SOLVER cs_minres_qlp

#else

#define ISOGW real(kind=C_DOUBLE)
#define GWTYP double precision
#ifdef GWEN
#define SOLVER gw_en
#else
#define SOLVER gw_minres
#endif
#endif

#define USE_ALLOC 1

      subroutine gw_buildrn_minres(pars,Eig,Eia,omega,R,dR,imo,
     $                      maxpoles,nri,me,nmo,ipol,isp)
      USE ISO_C_BINDING
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: imo, maxpoles, nri, me, nmo, ipol, isp
      double precision, intent(in) :: omega
      double precision, intent(in), dimension(nmo) :: Eig
      double precision, intent(in), dimension(maxpoles,ipol) :: Eia
      GWTYP, intent(out) :: R,dR

      character(len=12), parameter :: pname = 'gw_buildrn: '
      integer, parameter :: maxthread = 1024

      integer root, oolo, vvlo, ovlo, ovhi, nspace
      integer nprocs
      integer jmo, ipole, info, ld, mynri, idx, jsp
      integer k_ipiv, l_scr, k_scr, lscr
      integer g_pi, g_sol, off
      integer kOV(2), ga, first, second, step
      integer spacelo, addr, gb, owner, lo, hi
      double precision fac, arg, arg2, dfac, sgn, t1
      double precision vector(nri)
      double complex :: temp
      integer(kind=c_int)         :: stat_alloc
      ISOGW,allocatable,dimension(:,:) :: factor, dfactor
      ISOGW,dimension(nri) :: xvec

      logical main, ok, parallel_diag
      logical, external :: is_near

#ifdef GWCMPLX
      double precision, parameter :: eta = 0.001d0
#else
      double precision, parameter :: eta = 0.00001d0
#endif

      double precision, parameter :: shift = 0.00001d0
      double complex :: w, R1, dR1
      double complex, parameter :: ieta = dcmplx(0d0,eta)

#ifdef USE_OPENMP      
      integer iMaxThreads
      integer,external :: omp_get_max_threads, omp_get_num_threads
      integer,external :: omp_get_thread_num
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif

      call ga_sync()

      main = me.eq.0
      nprocs = ga_nnodes()

      allocate(factor(maxpoles,ipol),dfactor(maxpoles,ipol),
     $         stat=stat_alloc)
      if (stat_alloc.ne.0) then
        call errquit('gw_buildrn: allocation failed',0,0)
      endif

      do jsp=1,ipol
        if (pars%mynpoles(jsp).eq.0) cycle
        call ga_access(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                 pars%ovhi(jsp),kOV(jsp),ld)   
      enddo
      if (ipol.eq.1) kOV(2) = 1

      R1 = (0d0,0d0)
      dR1 = (0d0,0d0)

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
        else
          fac = sgn
          arg = eig(jmo) - omega
        endif


        do jsp=1,ipol
!$omp   parallel do simd private(temp) schedule(static)
          do ipole=1,pars%mynpoles(jsp)
            temp = Eia(ipole,jsp) + ieta

#ifdef GWCMPLX
            factor(ipole,jsp) = 0.5d0/(arg + temp) + 0.5d0/(temp - arg)
            dfactor(ipole,jsp) = dfac*(-(0.5d0/(arg + temp))**2 + 
     $                                  (0.5d0/(temp - arg))**2)
#else
            factor(ipole,jsp) = dble( 0.5d0/(arg + temp) + 
     $                                0.5d0/(temp - arg) )
            dfactor(ipole,jsp) = dfac*dble( (-(0.5d0/(arg + temp))**2 +
     $                                      (0.5d0/(temp - arg))**2) )
#endif

#ifdef GWDEBUG
            if ((factor(ipole,jsp).ne.factor(ipole,jsp)).or.
     $           dfactor(ipole,jsp).ne.dfactor(ipole,jsp)) then
              write(*,*) " ## Singular factor in buildrn. Shift denom"
              temp = temp + dcmplx(0.0d0,0.001d0)
#ifdef GWCMPLX
            factor(ipole,jsp) = 0.5d0/(arg + temp) + 0.5d0/(temp - arg)
            dfactor(ipole,jsp) = dfac*(-(0.5d0/(arg + temp))**2 + 
     $                                  (0.5d0/(temp - arg))**2)
#else
            factor(ipole,jsp) = dble( 0.5d0/(arg + temp) + 
     $                                0.5d0/(temp - arg) )
            dfactor(ipole,jsp) = dfac*dble( (-(0.5d0/(arg + temp))**2 +
     $                                      (0.5d0/(temp - arg))**2) )
#endif
            endif
#endif
          enddo
!$omp   end parallel do simd
        enddo

        root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
        idx = ((imo-off-1)*nspace+jmo-off-1)/nprocs + spacelo

        if (.not.ga_locate(gb,1,jmo,owner))
     $   call errquit('gw_buildrn_minres: cannot locate owner',0,GA_ERR)

        if (me.eq.owner) call ga_get(gb,1,nri,jmo,jmo,xvec,nri)
        if (me.eq.root) call ga_get(ga,1,nri,idx,idx,vector,nri)
        call ga_brdcst(mt_dbl,vector,ma_sizeof(mt_dbl,nri,mt_byte),root)
        call ga_brdcst(mt_dbl,xvec,ma_sizeof(mt_dbl,nri,mt_byte),owner)


        call SOLVER(dbl_mb(kOV(1)),dbl_mb(kOV(2)),factor,xvec,
     $              vector,pars%mynpoles,maxpoles,nri,ipol)

        if (owner.eq.me) call ga_put(gb,1,nri,jmo,jmo,xvec,nri)

        if (arg.ne.0d0) then
          do jsp=1,ipol
            if (pars%mynpoles(jsp).eq.0) cycle
            addr = kOV(jsp)
#ifdef GWCMPLX            
            call dgemm('n','n',2,pars%mynpoles(jsp),nri,1d0,xvec,2,
     $                  dbl_mb(addr),nri,0d0,factor,2)
#else
            call dgemv('t',nri,pars%mynpoles(jsp),1d0,dbl_mb(addr),nri,
     $                  xvec,1,0d0,factor,1)
#endif            

!$omp  parallel do simd reduction(+:dR1) schedule(static)
            do ipole=1,pars%mynpoles(jsp)
#ifdef GWCMPLX
               dR1 = dR1 + 
     $               fac*dfactor(ipole,jsp)*factor(ipole,1)**2
#else
               dR1 = dR1 + 
     $               fac*dble(dfactor(ipole,jsp)*factor(ipole,1)**2)
#endif
            enddo
!$omp  end parallel do simd

          enddo
        endif

        if (root.ne.me) cycle

c        w = ddot(nri,vector,1,xvec,1) - ddot(nri,vector,1,vector,1)
        w = sum(vector(:)*xvec(:)) - dot_product(vector,vector)

        R1 = R1 + fac*w
        if (abs(fac).lt.0.6d0) dR1 = dR1 + sign(1d0,fac)*w

      enddo

      call ga_zgop(mt_dbl,R1,1,'+')
      call ga_zgop(mt_dbl,dR1,1,'+')

      deallocate(factor,dfactor)

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

#ifdef GWCMPLX
      R = R1
      dR = dR1
#else
      R = dble(R1)
      dR = dble(dR1)
#endif

      return
      end subroutine
