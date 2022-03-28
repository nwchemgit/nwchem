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

      subroutine gw_cdgw_buildrn_graph(pars,Eig,Eia,omega,output,llmo,
     $   maxpoles,nri,me,nmo,ipol,isp,nqp)
C
C     Computes the contribution to the self-energy matrix elements, and
C     their derivatives with respect to omega, from the residues inside 
C     the contours. 
C
C     Either MINRES or the Bunch-Kaufmann triangular factorization can
C     be used, depending on the number of eigenvalues within the current
C     cluster.
C
C     See Equations (42)-(47) of
C     JCTC 17, 7504 (2021)
C
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
      logical main, ok, active, skip
      integer imo, jmo, kmo, lmo, ints
      integer mynpoles, ipole, kov(2), ga
      integer root, idx, nspace, nprocs, first, second, step
      integer info, ld, jsp, off, spacelo, owner,mynmo
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, lscr
      integer mynri, llpi, pilo, myimo
      double precision fac, arg, w, arg2, dfac, sgn, dummy

      GWTYPE :: factor(maxpoles,ipol)
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

      !
      ! Initialization
      !
      ieta = dcmplx(0d0,pars%eta)
      main = me.eq.0
      nprocs = ga_nnodes()
      ints = ma_sizeof(mt_int,1,mt_byte)
      output(:,:) = 0d0

      !
      ! Allocate local fields
      !
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

      !
      ! If the size of the dielectric matrix gets too large, opt for a
      ! parallel diagonalization instead
      !
      parallel_diag = nri .gt. 20000
#ifdef GWCMPLX
      paralle_diag = .false.
#endif

      if (parallel_diag) then
        mynri = nri/nprocs
        pilo = me*mynri + min(me,mod(nri,nprocs)) + 1
        llpi = k_pi + (pilo-1)*nri
        if (me.lt.mod(nri,nprocs)) mynri = mynri + 1
      endif

      !
      ! Set the arrays according to the occ/vir label of the current
      ! orbital. For occupied orbitals, the loop will go from the Fermi
      ! level down, while for unocuppied orbitals the loop will go from
      ! the Fermi level up.
      !
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

      !
      ! The sign of the contribution will depend on wether omega is
      ! below or above the Fermi level (first or second term of Equation
      ! (42)).
      !
      if (omega.lt.0d0) then
        sgn = -1d0
      else
        sgn = 1d0
      endif

      do jmo=first,second,step

        !
        ! Only orbitals closer to the Fermi level than omega contribute
        ! to the resiude term (the Heaviside step functions in Equation
        ! (42) take care of this). 
        ! 
        ! Skip the rest, but leave some space for possible degeneracies.
        ! 
        skip = .false.
        if (omega.lt.0d0) then
          if (eig(jmo) .lt. omega-shift) skip = .true.
          if (eig(jmo) .gt. 0d0) skip = .true.
        else
          if (eig(jmo) .gt. omega+shift) skip = .true.
          if (eig(jmo) .lt. 0d0) skip = .true.
        endif

        if (is_near(eig(jmo),omega,shift)) then
          skip = .true.
          fac = sgn*0.5d0
          arg = 0d0
        else
          fac = sgn
          arg = eig(jmo) - omega
        endif

        !
        ! Make sure all MPI ranks skip the same set of orbitals
        !
        call ga_brdcst(1038,skip,ma_sizeof(mt_log,1,mt_byte),0)
        if (skip) cycle

        !
        ! Precompute the factor in square brackets of Equation (39)
        !
        do jsp=1,ipol
!$omp     parallel do simd private(temp) schedule(static)
          do ipole=1,pars%mynpoles(jsp)
            temp = Eia(ipole,jsp) + ieta
#ifdef GWCMPLX
            factor(ipole,jsp) = 0.5d0/(arg + temp) + 0.5d0/(temp - arg)
#else
            factor(ipole,jsp) = dble( 0.5d0/(arg + temp) + 
     $                                0.5d0/(temp - arg) )
#endif
          enddo
!$omp     end parallel do simd
        enddo

        !
        ! Gather all ERI vectors which are local to the current MPI rank
        !
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

        !
        ! Build the polarizability matrix \Pi and transform it into the
        ! dielectric matrix \epsilon. Equation (39).
        !
        ! In real arithmetics, the dielectric matrix is obtained in RFP
        ! format, so we need to take care to bring it back to full
        ! format for subsequent steps.
        !
#ifdef GWCMPLX
        call gw_cdgw_buildpi_cmplx('r',dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $              GWMB(k_pi),factor,nri,nri,pars%mynpoles,main,ipol,
     $              maxpoles)
        call GWGOP(GWMT,GWMB(k_pi),nri**2,'+')
#else        
        call gw_cdgw_buildpi('r',dbl_mb(kOV(1)),dbl_mb(kOV(2)),
     $         GWMB(k_pi),factor,nri,nri,pars%mynpoles,main,ipol,
     $         maxpoles)

        call GWGOP(GWMT,GWMB(k_pi),(nri*(nri+1))/2,'+')
        call rfp2std(GWMB(k_pi),nri)
#endif        

        !
        ! If parallel_diag = .true., invert the dielectric matrix using
        ! a parallel eigensolver.
        !
        if (parallel_diag) then
          call gw_invdiag(dbl_mb(k_pi),nri,nri,dummy,.false.)
          call GWGOP(GWMT,GWMB(k_pi),nri**2,'+')
          if (.not.active) cycle
          call dgemm('t',nri,mynmo,nri,1d0,GWMB(k_pi),nri,vectors,
     $                nri,0d0,vectors2,nri)
        else
          if (.not.active) cycle

          !
          ! Get \epsilon_{PQ}^{-1} (Q|mn) for all active mn pairs. This
          ! step uses MINRES for a small number of pairs, and switches
          ! to Baunch-Kauffmann LDL^T factorization otherwise.
          !
          ! Obtains the `x` vectors of Equations (43), (45)
          !
#ifdef GWCMPLX        
          vectors2(:,1:mynmo) = dcmplx(vectors(:,1:mynmo),0d0)
#else
          vectors2(:,1:mynmo) = vectors(:,1:mynmo)
#endif
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

        !
        ! Complete the contribution in Equations (43), (45) by
        ! contracting with the untransformed ERIs (vectors `b`).
        !
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
      
      !
      ! Sum contributions from all MPI ranks.
      !
      call ga_dgop(1038,output,2*nqp,'+')


      if (.not.ma_chop_stack(l_ipiv))
     $   call errquit(pname//'can''t chop stack',0,MA_ERR)
      if (.not.ma_free_heap(l_pi))
     $   call errquit(pname//'can''t free heap',0,MA_ERR)
c     

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      do jsp=1,ipol
        call ga_release(pars%g_eriov(jsp),1,nri,pars%ovlo(jsp),
     $                  pars%ovhi(jsp))
      enddo

      return
      end subroutine
