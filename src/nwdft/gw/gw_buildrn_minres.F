#ifdef GWCMPLX
#define GWTYP double complex
#define SOLVER cs_minres_qlp
#else
#define GWTYP double precision
#ifdef GWEN
#define SOLVER gw_en
#else
#define SOLVER gw_minres
#endif
#endif

      subroutine gw_buildrn_minres(pars,Eig,Eia,omega,R,dR,imo,
     $                      maxpoles,nri,me,nmo,ipol,isp)
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
      double precision, intent(out) :: R,dR

      character(len=12), parameter :: pname = 'gw_buildrn: '

      integer root, oolo, vvlo, ovlo, ovhi, nspace
      integer nprocs
      integer jmo, ipole, info, ld, mynri, idx, jsp
      integer k_ipiv, l_scr, k_scr, lscr
      integer g_pi, g_sol, off
      integer kOV(2), ga, first, second, step
      integer spacelo, addr, gb, owner
      double precision fac, arg, w, arg2, dfac, sgn
      double precision vector(nri)
      double complex :: temp
      GWTYP,dimension(maxpoles,ipol) :: factor, dfactor
      GWTYP,dimension(nri) :: xvec

      logical main, ok, parallel_diag
      logical, external :: is_near

      double precision, parameter :: eta = 0.001d0
      double precision, parameter :: shift = 0.00001d0
      double complex, parameter :: ieta = dcmplx(0d0,eta)

#ifdef USE_OPENMP      
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif


      main = me.eq.0
      nprocs = ga_nnodes()

      do jsp=1,ipol
        if (pars%mynpoles(jsp).eq.0) cycle
        call ga_access(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                 pars%ovhi(jsp),kOV(jsp),ld)   
      enddo
      if (ipol.eq.1) kOV(2) = 1

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
        else
          fac = sgn
          arg = eig(jmo) - omega
        endif

!$omp   parallel do collapse(2) private(temp)
        do jsp=1,ipol
          do ipole=1,pars%mynpoles(jsp)
            temp = Eia(ipole,jsp) - ieta
            temp = 0.5d0/(arg + temp) + 0.5d0/(temp - arg)
#ifdef GWCMPLX
            factor(ipole,jsp) = temp
#else            
            factor(ipole,jsp) = dble(temp)
#endif            

            temp = Eia(ipole,jsp) - ieta
            temp = -(0.5d0/(arg + temp))**2 + (0.5d0/(temp - arg))**2
#ifdef GWCMPLX
            dfactor(ipole,jsp) = dfac*temp
#else            
            dfactor(ipole,jsp) = dfac*dble(temp)
#endif            
          enddo
        enddo
!$omp   end parallel do        

        root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
        idx = ((imo-off-1)*nspace+jmo-off-1)/nprocs + spacelo


        call ga_brdcst(mt_int,idx,ma_sizeof(mt_int,1,mt_byte),root)

        if (.not.ga_locate(gb,1,jmo,owner))
     $  call errquit('gw_buildrn_minres: cannot locate owner',0,GA_ERR)


        call ga_get(gb,1,nri,jmo,jmo,xvec,nri)
        call ga_get(ga,1,nri,idx,idx,vector,nri)

        call SOLVER(dbl_mb(kOV(1)),dbl_mb(kOV(2)),factor,xvec,vector,
     $              pars%mynpoles,maxpoles,nri,ipol)

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

!$omp       parallel do reduction(+:dR)
            do ipole=1,pars%mynpoles(jsp)
               dR = dR + fac*dble(dfactor(ipole,jsp)*factor(ipole,1)**2)
            enddo
!$omp       end parallel do

          enddo
        endif

        if (root.ne.me) cycle

        w = dot_product(vector,dble(xvec)) - dot_product(vector,vector)
        R = R + fac*w

        if (abs(fac).lt.0.6d0) dR = dR + sign(1d0,fac)*w

      enddo

      call ga_dgop(mt_dbl,R,1,'+')
      call ga_dgop(mt_dbl,dR,1,'+')


#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      return
      end subroutine
