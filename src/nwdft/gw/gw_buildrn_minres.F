      subroutine gw_buildrn_minres(pars,Eig,Eia,omega,R,dR,d2R,imo,
     $                      maxpoles,nri,me,nmo,ipol,isp)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: imo, maxpoles, nri, me, nmo, ipol, isp
      double precision, intent(in) :: omega
      double precision, intent(in), dimension(nmo) :: Eig
      double precision, intent(in), dimension(maxpoles,ipol) :: Eia
      double precision, intent(out) :: R,dR,d2R

      character(len=12), parameter :: pname = 'gw_buildrn: '

      integer mynpoles, root2
      integer root, oolo, oohi, vvlo, vvhi, nrisize, ovlo, ovhi, nspace
      integer g_a, g_b, g_c, g_d, g_e, g_epsilon, k_eps, k_sol, nprocs
      integer jmo, ipole, ilo, ihi, jlo, jhi, info, ld, mynri, idx, jsp
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, k_ov, lscr
      integer g_pi, g_sol, addrPi, addrSol, llpi, off
      integer pihi, pilo, kOV(2), ga, first, second, step
      integer spacelo, addr, gb, owner
      double precision fac, arg, w, arg2, dfac, sgn
      double precision vector(nri), vector3(nri)
#ifdef GWCMPLX
      double complex factor(maxpoles,ipol),dfactor(maxpoles,ipol)
      double complex xvec(nri)
      double precision, parameter :: eta = 0.003d0
      double complex, parameter :: ieta = dcmplx(0d0,eta)
#else
      double precision factor(maxpoles,ipol),dfactor(maxpoles,ipol)
      double precision xvec(nri)
      double precision,parameter :: ieta = 0d0
#endif      
      double precision, parameter :: shift = 0.00001d0

      logical main, ok, parallel_diag
      double precision, external :: ddot
      logical, external :: is_near, nga_create, ga_duplicate, ga_destroy
      logical, external :: ga_create, ga_copy, ga_locate
      integer,external :: ga_nodeid,ga_nnodes
#ifdef USE_OPENMP      
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif


      main = me.eq.0
      nprocs = ga_nnodes()

      do jsp=1,ipol
        call ga_access(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                 pars%ovhi(jsp),kOV(jsp),ld)   
      enddo
      if (ipol.eq.1) kOV(2) = 1

      R = 0d0
      dR = 0d0
      d2R = 0d0

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

       
      do jmo=first,second,step
        if ( abs(eig(jmo)) .gt. abs(omega+sgn*shift) ) exit
        if (is_near(eig(jmo),omega,shift)) then
          fac = sgn*0.5d0
          arg = 0d0
          arg2 = 0d0
        else
          fac = sgn
          arg = sgn*(omega - eig(jmo))
          arg2 = arg**2
        endif

!$omp   parallel do collapse(2)
        do jsp=1,ipol
          do ipole=1,pars%mynpoles(jsp)
            factor(ipole,jsp) = (Eia(ipole,jsp)+ieta)/
     $                          ((Eia(ipole,jsp)+ieta)**2-arg2)
            dfactor(ipole,jsp) = -dfac*arg*factor(ipole,jsp)/
     $                          ((Eia(ipole,jsp)+ieta)**2-arg2)
          enddo
        enddo
!$omp   end parallel do        

        root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
        idx = ((imo-off-1)*nspace+jmo-off-1)/nprocs + spacelo

        call ga_brdcst(mt_int,idx,mitob(1),root)

        if (.not.ga_locate(gb,1,jmo,owner))
     $  call errquit('gw_buildrn_minres: cannot locate owner',0,GA_ERR)

        call ga_get(gb,1,nri,jmo,jmo,xvec,nri)
        call ga_get(ga,1,nri,idx,idx,vector,nri)

#ifdef GWCMPLX
        call cs_minres_qlp(dbl_mb(kOV(1)),dbl_mb(kOV(2)),factor,xvec,
     $                     vector,pars%mynpoles,maxpoles,nri,ipol)
#elif GWMINRES
        call gw_minres(dbl_mb(kOV(1)),dbl_mb(kOV(2)),factor,xvec,
     $            vector,pars%mynpoles,maxpoles,nri,ipol,20,.true.)
#else
        call gw_en(dbl_mb(kOV(1)),dbl_mb(kOV(2)),factor,xvec,
     $                 vector,pars%mynpoles,maxpoles,nri,ipol,
     $                 20,.true.)
#endif        
        if (owner.eq.me) call ga_put(gb,1,nri,jmo,jmo,xvec,nri)

        if (arg.ne.0d0) then
          do jsp=1,ipol
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
               dR = dR + dble(dfactor(ipole,jsp)*factor(ipole,1)**2)
            enddo
!$omp       end parallel do

          enddo
        endif

        if (root.ne.me) cycle

        w = dot_product(vector,dble(xvec)) - dot_product(vector,vector)
        R = R + fac*w

        if (abs(fac).lt.1d0) dR = dR + w

      enddo

      call ga_dgop(mt_dbl,R,1,'+')
      call ga_dgop(mt_dbl,dR,1,'+')


#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      return
      end subroutine
