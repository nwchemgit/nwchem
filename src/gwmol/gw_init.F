      subroutine gw_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "bas.fh"
#include "case.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"      
#include "global.fh"      
#include "gw.fh"
#ifdef SCALAPACK
#include "scaleig.fh"
#endif
c
      type(gw_params_t) :: pars
c
      integer i, ilo, ihi, j
      integer nbf_temp, ipol_temp, nocc, nvir, nmo_temp(2)
c
      character(len=9), parameter :: pname = 'gw_init: '
      character(len=255) :: basisname, scftype

      logical init
c
      logical int_normalize, int_norm_2c
      external int_normalize, int_norm_2c
c
      logical movecs_read, movecs_read_header, atom_tag_check
      external movecs_read, movecs_read_header, atom_tag_check
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      integer resultlen,tid,cpu,ierr,iproc,itid
      character(len=80) names
      double precision t0, ehomo, vxddp, vxddm, aq2, zeta

#ifdef USE_OPENMP
      integer,external :: omp_get_max_threads, omp_get_thread_num
      integer,external :: omp_get_num_threads
#endif      
c
c     -----------------------------------------------------------------
c
      zeta = 0.304121d0
      aq2 = 2d0/9d0 * (1d0-zeta)**2

      pars%timer%init  = util_wallsec()
      pars%timer%total = pars%timer%init
      pars%timer%vxc   = 0d0
      pars%timer%eris  = 0d0
      pars%timer%distr = 0d0
      pars%timer%sigmax= 0d0
      pars%timer%sigmac= 0d0
      pars%timer%rpa   = 0d0
      pars%timer%wm    = 0d0
      pars%timer%rn    = 0d0
      pars%timer%in    = 0d0
      pars%timer%iw    = 0d0
      pars%timer%vmo   = 0d0
      pars%eviter      = 0

      pars%ecrpa = 0d0
      pars%ecgm  = 0d0
      pars%ex    = 0d0

      pars%me = ga_nodeid()
      pars%nprocs = ga_nnodes()

#ifdef USE_OPENMP
      pars%iMaxThreads = omp_get_max_threads()
#else
      pars%iMaxThreads = 1
#endif

      if (.not.rtdb_get(pars%rtdb,'gw:init',mt_log,1,init))
     $  call gw_defaults(pars%rtdb) 

      if (.not.rtdb_get(pars%rtdb,'dft:energy',mt_dbl,1,pars%escf))
     $  call errquit(pname//'failed to read DFT energy',0,RTDB_ERR) 

      if (.not.rtdb_get(pars%rtdb, 'gw:eta',mt_dbl,1,pars%eta))
     $  call errquit(pname//'failed to read eta',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:evgw0',mt_log,1,pars%evgw0))
     $  call errquit(pname//'failed to read evgw0',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:evgw',mt_log,1,pars%evgw))
     $  call errquit(pname//'failed to read evgw',0,RTDB_ERR)
      if (pars%evgw0 .or. pars%evgw) then
        if (.not.rtdb_get(pars%rtdb,'gw:eviter',mt_int,1,pars%eviter))
     $    call errquit(pname//'failed to read eviter',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(pars%rtdb,'gw:cdgw',mt_log,1,pars%cdgw))
     $  call errquit(pname//'failed to read cdgw',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:acgw',mt_log,1,pars%acgw))
     $  call errquit(pname//'failed to read acgw',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:ncap',mt_log,1,pars%ncap))
     $  call errquit(pname//'failed to read ncap',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:rpa',mt_log,1,pars%rpa))
     $  call errquit(pname//'failed to read rpa',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:diag',mt_log,1,pars%diag))
     $  call errquit(pname//'failed to read diag',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:core',mt_log,1,pars%core))
     $  call errquit(pname//'failed to read diag',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:first',mt_int,1,pars%first))
     $  pars%first = 1
c
      if (.not.rtdb_get(pars%rtdb,'gw:graph',mt_log,1,pars%graph))
     $  call errquit(pname//'failed to read graph',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:newton',mt_log,1,pars%newton))
     $  call errquit(pname//'failed to read newton',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:linear',mt_log,1,pars%linear))
     $  call errquit(pname//'failed to read linear',0,RTDB_ERR)

      if (pars%linear .or. pars%graph ) then
        pars%newton = .false.
      endif

      if (pars%newton) then
        if (.not.rtdb_get(pars%rtdb,'gw:qpiter',mt_int,1,pars%qpiter))
     $    call errquit(pname//'failed to read qpiter',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(pars%rtdb,'gw:threshold',mt_dbl,1,pars%thresh))
     $  call errquit(pname//'failed to read threshold',0,RTDB_ERR)
c
      if (.not.rtdb_get(pars%rtdb,'gw:noqp_alpha',mt_int,1,
     $                  pars%noqp(1)))
     $  call errquit(pname//'failed to read noqp',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:nvqp_alpha',mt_int,1,
     $                  pars%nvqp(1)))
     $  call errquit(pname//'failed to read nvqp',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:noqp_beta',mt_int,1,pars%noqp(2)))
     $  call errquit(pname//'failed to read noqp',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:nvqp_beta',mt_int,1,pars%nvqp(2)))
     $  call errquit(pname//'failed to read nvqp',0,RTDB_ERR)
c
      if (.not.rtdb_get(pars%rtdb,'gw:analyticw',mt_log,1,pars%anaw))
     $  call errquit(pname//'failed to read analyticw',0,RTDB_ERR)
c
      call dft_rdinput(pars%rtdb)
c      
      if (cdfit) then
        if (.not.bas_destroy(cd_bas_han))      
     $    call errquit(pname//'failed to destroy cd_bas',0,BASIS_ERR)
      endif

      if (.not.rtdb_get(pars%rtdb, 'dft:itol2e',mt_int,1,itol2e))
     $  call errquit(pname//'failed to read itol2e',0,RTDB_ERR)
      pars%tol2e = 10.0d0**(-itol2e-2)
c
c     RI basis set
c
      if (.not.bas_create(cd_bas_han,'ri basis'))
     $  call errquit(pname//'bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(pars%rtdb,geom,cd_bas_han,'ri basis'))then
        if (.not.bas_rtdb_load(pars%rtdb,geom,cd_bas_han,'cd basis'))
     $  call errquit(pname//'a "cd basis" or "ri basis" is needed',0,
     $  BASIS_ERR)
      endif
      call int_init(pars%rtdb, 1, cd_bas_han)
      if (.not.int_norm_2c(pars%rtdb, cd_bas_han))
     $  call errquit(pname//'int_norm_2c failed',0,INT_ERR)
      call int_terminate()
      if (.not. bas_numbf(cd_bas_han, nbf_cd))
     $  call errquit(pname//'basis set error',0,BASIS_ERR)
      if (.not. bas_nprim_cn_max(cd_bas_han,nbf_cd_mxprim))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(cd_bas_han,nbf_cd_mxang))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(cd_bas_han,nbf_cd_mxcont))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(cd_bas_han,nbf_cd_mxnbf_cn))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(cd_bas_han,nbf_cd_mxnbf_ce))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(cd_bas_han,nshells_cd))
     $  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (bas_is_spherical(ao_bas_han).and.
     $   (.not.bas_is_spherical(cd_bas_han)))
     $        call int_app_set_no_texas(pars%rtdb)
c
c     Spin multiplicity, occupations
c
      if (.not.rtdb_get(pars%rtdb,'dft:ipol',mt_int,1,ipol))
     $  call errquit(pname//'failed to read ipol',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb, 'dft:noc',mt_int,2,noc))
     $  call errquit(pname//'failed to read noc',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:mult',mt_int,1,mult))
     $  call errquit(pname//'failed to read mult',0,RTDB_ERR)
      pars%ipol = ipol
c
c     MO vectors
c
      if (.not. rtdb_cget(pars%rtdb, 'dft:output vectors',1,movecs_out))
     $     call errquit(pname//'no final MO found',0,RTDB_ERR)
c
c     Allocate arrays
c
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,
     $    'eigenvalues', pars%l_evals, pars%k_evals))
     $  call errquit(pname//'failed to allocate evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,'occupancies', pars%l_occ,
     $  pars%k_occ))
     $  call errquit(pname//'failed to allocate occ',0,MA_ERR) 
      if (.not.movecs_read_header(movecs_out,title,basisname,
     $    scftype,nbf_temp,ipol_temp,nmo_temp,2))
     $    call errquit(pname//'failed to read MO header',0,DISK_ERR)
      if (nbf_ao.ne.nbf_temp)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      if (ipol.ne.ipol_temp)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)

      pars%nmo = nmo_temp(1)

      pars%nocc(:) = 0
      pars%mynpoles(:) = 0
      pars%nvir(:) = 0
      pars%nqp(:) = 0
c
      do i=1,ipol
        pars%g_movecs(i) = ga_create_atom_blocked(geom,
     $    ao_bas_han,'MO eigenvectors')
        if (.not.movecs_read(movecs_out,i,
     $      dbl_mb(pars%k_occ+(i-1)*nbf_ao),
     $      dbl_mb(pars%k_evals+(i-1)*nbf_ao),
     $      pars%g_movecs(i)))
     $    call errquit(pname//'failed to read MO vectors',0,DISK_ERR)

        nocc = 0
        nvir = 0
        do j=1,pars%nmo
          if (dbl_mb(pars%k_occ+(i-1)*nbf_ao+j-1).gt.0.d0) then
            nocc = nocc + 1
            cycle
          else
            nvir = nvir + 1
          endif
        enddo
        pars%nocc(i) = nocc
        pars%nvir(i) = nvir
        if (nocc+nvir.ne.pars%nmo)
     $    call errquit(pname//'something went wrong',0,0)
      
c
        if ((pars%noqp(i).gt.nocc).or.(pars%noqp(i).eq.-1)) 
     $     pars%noqp(i) = nocc

        if (pars%noqp(i).lt.-1)
     $    call errquit(pname//'invalid number noqp specified',0,0)
c
        if ((pars%nvqp(i).gt.nvir).or.(pars%nvqp(i).eq.-1)) 
     $     pars%nvqp(i) = nvir

        if (pars%nvqp(i).lt.-1)
     $    call errquit(pname//'invalid number nvqp specified',0,0)

        if (pars%core) then
          pars%llmo(i) = pars%first
        else
          pars%llmo(i) = pars%nocc(i) - pars%noqp(i) + 1
        endif
        pars%nqp(i) = pars%noqp(i) + pars%nvqp(i)
c      

      enddo
c
c      if (pars%evgw.or.pars%evgw0) then
c        pars%noqp(:) = pars%nocc(:)
c        pars%nvqp(:) = pars%nvir(:)
c      endif
c
      if (sum(pars%noqp)+sum(pars%nvqp).eq.0)
     $  call errquit(pname//'no quasiparticle energies to obtain',0,0)
c
      do i=1,ipol
        pars%nstates(i) = pars%nocc(i) + pars%nvqp(i)
        pars%npoles(i) = pars%nocc(i)*pars%nvir(i)
      enddo
c
 1001 format(2X,A12,F8.3,A3)
      if (pars%ncap) then
        do i=1,ipol
          ehomo = dbl_mb(pars%k_evals+(i-1)*nbf_ao+pars%nocc(i)-1)
          vxddp = -0.5d0*aq2*(1d0 + sqrt(1d0-4d0*ehomo/aq2))
          vxddm = -0.5d0*aq2*(1d0 - sqrt(1d0-4d0*ehomo/aq2))
          pars%vxddp(i) = vxddp
          pars%vxddm(i) = vxddm
          if (pars%me.eq.0) then
            write(luout,1001) ' v_x^{DD+}: ',vxddp*ha2ev, ' eV'
            write(luout,1001) ' v_x^(DD-): ',vxddm*ha2ev, ' eV'
          endif
        enddo
      endif
c
      if (pars%evgw .or. pars%evgw0) then
        if (pars%eviter.lt.1)
     $    call errquit(pname//'invalid eviter',0,INPUT_ERR)
      endif
c      
      if (pars%newton .and. pars%qpiter.lt.1)
     $  call errquit(pname//'invalid qpiter',0,INPUT_ERR)
c
c     Gauss-Legendre quadrature
c
      if (pars%cdgw) then
        if (.not.rtdb_get(pars%rtdb,'gw:ngl',mt_int,1,pars%ngl))
     $    call errquit(pname//'failed to read ngl',0,RTDB_ERR)
        if (.not.ma_push_get(mt_dbl,pars%ngl,'GL roots',
     $    pars%l_glx,pars%k_glx))
     $    call errquit(pname//'can''t allocate GL roots',0,MA_ERR) 
        if (.not.ma_push_get(mt_dbl,pars%ngl,'GL weights',
     $    pars%l_glw,pars%k_glw))
     $    call errquit(pname//'can''t allocate GL roots',0,MA_ERR) 
        call gw_gaussleg(pars%ngl,dbl_mb(pars%k_glx),dbl_mb(pars%k_glw))
c        call gw_dequad(pars%ngl,dbl_mb(pars%k_glx),dbl_mb(pars%k_glw))
      endif
c
c     3-center ERIs 
c
      call int_init(pars%rtdb, 2, (/ao_bas_han, cd_bas_han/))
      call print_integrals((/ao_bas_han,cd_bas_han/),.false.)
      call schwarz_init(geom, ao_bas_han)
      call scf_get_fock_param(pars%rtdb, pars%tol2e)

      pars%timer%init = util_wallsec() - pars%timer%init
c
      call grid_cleanup(.false.)
      call gw_memory(pars)
      call gw_ri_init(pars)
c
      pars%fdstep = 0.0005d0

#if defined(PARALLEL_DIAG) && defined(SCALAPACK)
      call slexit
      call slexit2
      call slexit3
      call slexit4
      if (.not. rtdb_cget(pars%rtdb, 'dft:scaladiag',1,yev))
     &     yev='d'
      if (.not. rtdb_get(pars%rtdb, 'dft:scalanb',mt_int,1,nb))
     &     nb=32

      if (pars%cdgw) then
        call slinit2(nbf_cd)
      else
        call slinit2(sum(pars%npoles(1:pars%ipol)))
      endif

#endif

      return
      end subroutine gw_init

      subroutine gw_gaussleg(ngl,x,w)
      implicit none
      integer,intent(in) :: ngl
      double precision,intent(inout),dimension(ngl) :: x, w

      integer igl,jgl,m
      double precision pi
      double precision p1,p2,p3,pp,z,z1

      pi = acos(-1d0)

      m = (ngl + 1)/2
      do igl=1,m
        z = dcos(pi*(igl-0.25d0)/(ngl+0.5d0))
 100    p1 = 1d0
        p2 = 0d0
        do jgl=1,ngl
          p3 = p2
          p2 = p1
          p1 = ((2d0*jgl-1d0)*z*p2 - (jgl-1d0)*p3)/jgl
        enddo
        pp = ngl*(z*p1-p2)/(z*z-1d0)
        z1 = z
        z = z1 - p1/pp
        if (abs(z-z1).gt.1d-15) goto 100

        x(igl) = -z
        x(ngl+1-igl) = z
        w(igl) = 2d0/((1d0-z*z)*pp*pp)
        w(ngl+1-igl) = w(igl)
      enddo

      do igl=1,ngl
        w(igl) = 0.5d0 * 2d0 * w(igl) / (1d0-x(igl))**2 / pi
        x(igl) = (0.5d0 * ((1d0 + x(igl))/(1d0 - x(igl))))**2
      enddo


      end subroutine

      subroutine gw_dequad(ngl,x,w)
      implicit none
      integer,intent(in) :: ngl
      double precision,intent(out),dimension(ngl) :: x, w

      integer imin,imax,igl,jgl
      double precision alpha, h, s, arcsinh, rmax, rmin, pi
      double precision,dimension(ngl) :: tmp
      arcsinh(s) = dlog(s + dsqrt(1d0 + s**2))

      alpha = 2.2d0
      rmin = 0.01d0
      rmax = 1d6
      pi = dacos(-1d0)

      imin = -ceiling(dble(ngl-1)/(1d0 - arcsinh(dlog(rmax)/alpha)/
     $               arcsinh(dlog(rmin)/alpha)))
      imax = ngl + imin - 1
      h = arcsinh(dlog(rmin)/alpha)/dble(imin)

      igl = 0
      do jgl=imin,imax
        igl = igl + 1
        x(igl) = dexp(alpha*dsinh(jgl*h))**2
        w(igl) = alpha*h*dexp(3d0*alpha*dsinh(jgl*h))*dcosh(jgl*h)/
     $           pi/x(igl)**2
        write(*,*) x(igl)
      enddo
      return
      end subroutine



