#ifdef GWCMPLX
#define GWTYPE double complex
#define GWMT mt_dcpl
#define GWMB dcpl_mb      
#define GWSYTRF zsytrf
#define GWSYSV zsysv
#define GWGOP ga_zgop
#define GWGEMV zgemv      
#else
#define GWTYPE double precision
#define GWMT mt_dbl
#define GWMB dbl_mb
#define GWSYTRF dsytrf
#define GWSYSV dsysv      
#define GWGOP ga_dgop
#define GWGEMV dgemv      
#endif

      subroutine gw_buildrn_graph(pars,Eig,Eia,omega,output,llmo,
     $   maxpoles,nri,me,nmo,ipol,isp,nqp)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "gw.fh"
      type(gw_params_t) :: pars
      integer,intent(in) :: llmo,maxpoles,nri,me,nmo,ipol,isp,nqp
      double precision, intent(in) :: omega
      double precision, intent(in), dimension(nmo) :: Eig
      double precision, intent(in), dimension(maxpoles,ipol) :: Eia
      double precision, intent(out),dimension(2,nqp) :: output

      character(len=12), parameter :: pname = 'gw_buildrn_graph: '

      logical, external :: is_near

      logical parallel_diag
      logical main, ok, active
      integer imo, jmo, kmo, lmo, ints
      integer mynpoles, ipole, kov(2), ga
      integer root, idx, nspace, nprocs, first, second, step
      integer info, ld, jsp, off, spacelo, owner,mynmo
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, lscr
      integer mynri, llpi, pilo, myimo
      double precision fac, arg, w, arg2, dfac, sgn, dummy

      GWTYPE :: factor(maxpoles,ipol),dfactor(maxpoles,ipol)
      GWTYPE :: vectors2(nri,nqp)
      double complex :: temp
      double precision :: vectors(nri,nqp)

      double precision, parameter :: shift = 0.0001d0
      double complex :: ieta

#ifdef USE_OPENMP      
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif

      ieta = dcmplx(0d0,pars%eta)

      main = me.eq.0
      nprocs = ga_nnodes()
      ints = ma_sizeof(mt_int,1,mt_byte)

      if (.not.ma_alloc_get( GWMT ,nri**2,'pi',l_pi,k_pi))
     $  call errquit(pname//'can''t allocate pi matri',0,MA_ERR) 
      if (.not.ma_push_get(mt_int,nri,'ipiv',l_ipiv,k_ipiv))
     $  call errquit(pname//'can''t allocate ipiv',0,MA_ERR) 
      call GWSYTRF ('l',nri,GWMB(k_pi),nri,int_mb(k_ipiv),
     $               factor,-1,info)
      lscr = int(dble(factor(1,1)))

      if (.not.ma_push_get( GWMT ,lscr,'scr',l_scr,k_scr))
     $  call errquit(pname//'can''t allocate scratch',0,MA_ERR) 

      do jsp=1,ipol
        call ga_access(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                 pars%ovhi(jsp),kOV(jsp),ld)   
      enddo
      if (ipol.eq.1) kOV(2) = 1

#if PARALLEL_DIAG
      parallel_diag = nri .gt. 4000
#else
      parallel_diag = .false.
#endif

      if (parallel_diag) then
        mynri = nri/nprocs
        pilo = me*mynri + min(me,mod(nri,nprocs)) + 1
        llpi = k_pi + (pilo-1)*nri
        if (me.lt.mod(nri,nprocs)) mynri = mynri + 1
      endif

      output(:,:) = 0d0

      dfac = 4d0
      if (ipol.eq.1) dfac = 2d0*dfac

      if (llmo.le.pars%nocc(isp)) then
        ga = pars%g_erioo(isp)
        spacelo = pars%oolo(isp)
        first = pars%nocc(isp)
        second = 1
        step = -1
        off = 0
        nspace = pars%nocc(isp)
      else
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

        ! Loop over omega < eig(jmo) < 0d0
        !   or over 0d0 < eig(jmo) < omega
        if (omega.lt.0d0) then
          if (eig(jmo) .lt. omega-shift) cycle
          if (eig(jmo) .gt. 0d0) cycle
        else
          if (eig(jmo) .gt. omega+shift) cycle
          if (eig(jmo) .lt. 0d0) cycle
        endif

        ! Check for degeneracies
        if (is_near(eig(jmo),omega,shift)) then
          fac = sgn*0.5d0
          arg = 0d0
          arg2 = 0d0
        else
          fac = sgn
          arg = eig(jmo) - omega
          arg2 = arg**2
        endif

        !obtain the factor to build Pi matrix and its derivative
        do jsp=1,ipol
!$omp     parallel do simd private(temp) schedule(static)
          do ipole=1,pars%mynpoles(jsp)
            temp = (Eia(ipole,jsp) - ieta)/
     $      ((Eia(ipole,jsp)**2-arg2)-2d0*ieta*(Eia(ipole,jsp)-sgn*arg))
#ifdef GWCMPLX            
            factor(ipole,jsp) = temp
#else            
            factor(ipole,jsp) = dble(temp)
#endif            
          enddo
!$omp   end parallel do simd
        enddo

        !get all active vectors
        mynmo = 0
        active = .false.
        do imo=llmo,llmo+nqp-1
          root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
          if (root.ne.me) cycle
          active = .true.
          mynmo = mynmo+1
          idx = ((imo-off-1)*nspace+jmo-off-1)/nprocs + spacelo
          call ga_get(ga,1,nri,idx,idx,vectors(1,mynmo),nri)
        enddo

        !build Pi matrix
#ifdef GWCMPLX
        call gw_buildpi_cmplx('r',dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $              GWMB(k_pi),factor,nri,nri,pars%mynpoles,main,ipol,
     $              maxpoles)
#else        
        call gw_buildpi('r',dbl_mb(kOV(1)),dbl_mb(kOV(2)),GWMB(k_pi),
     $                   factor,nri,nri,pars%mynpoles,main,ipol,
     $                   maxpoles)
#endif        
        call GWGOP(GWMT,GWMB(k_pi),nri**2,'+')


        if (parallel_diag) then
          call gw_invdiag(dbl_mb(k_pi),nri,nri,dummy,.false.)
          call GWGOP(GWMT,GWMB(k_pi),nri**2,'+')
          if (.not.active) cycle
          call dgemm('t',nri,mynmo,nri,1d0,GWMB(k_pi),nri,vectors,
     $                nri,0d0,vectors2,nri)
        else
          if (.not.active) cycle
#ifdef GWCMPLX        
          vectors2(:,1:mynmo) = dcmplx(vectors(:,1:mynmo),0d0)
#else
          vectors2(:,1:mynmo) = vectors(:,1:mynmo)
#endif
          !solve the linear equation systems
          if (mynmo.gt.5) then
            call GWSYSV('l',nri,mynmo,GWMB(k_pi),nri,int_mb(k_ipiv),
     $                  vectors2,nri,GWMB(k_scr),lscr,info)
          else
            do myimo=1,mynmo
              vectors2(:,myimo) = 0d0
              call minres(GWMB(k_pi),vectors2(1,myimo),
     $                    vectors(1,myimo),nri)
            enddo
          endif
        endif

        lmo = 0
        do imo=llmo,llmo+nqp-1
          root = mod((imo-off-1)*nspace+jmo-off-1,nprocs)
          kmo = imo - llmo + 1
          if (root.eq.me) then            
            lmo = lmo + 1
            w = dot_product(vectors(:,lmo),dble(vectors2(:,lmo))) -
     $          dot_product(vectors(:,lmo),vectors(:,lmo))
            output(1,kmo) = output(1,kmo) + fac*w
          endif
        enddo

      enddo
       
      call ga_dgop(mt_dbl,output,2*nqp,'+')


      if (.not.ma_chop_stack(l_ipiv))
     $   call errquit(pname//'can''t chop stack',0,MA_ERR)
      if (.not.ma_free_heap(l_pi))
     $   call errquit(pname//'can''t free heap',0,MA_ERR)
c     

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      return
      end subroutine
