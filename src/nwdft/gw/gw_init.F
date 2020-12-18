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
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      integer i, ilo, ihi
      integer nbf_temp, ipol_temp, nocc, nvir
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
      logical ga_create_atom_blocked
      external ga_create_atom_blocked
c
c     -----------------------------------------------------------------
c
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

      if (.not.rtdb_get(pars%rtdb,'gw:init',mt_log,1,init))
     $  call gw_defaults(pars%rtdb) 

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
c
      if (.not.rtdb_get(pars%rtdb,'gw:graph',mt_log,1,pars%graph))
     $  call errquit(pname//'failed to read graph',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:newton',mt_log,1,pars%newton))
     $  call errquit(pname//'failed to read newton',0,RTDB_ERR)
      if (pars%newton) then
        if (.not.rtdb_get(pars%rtdb,'gw:qpiter',mt_int,1,pars%qpiter))
     $    call errquit(pname//'failed to read qpiter',0,RTDB_ERR)
      endif
c
      if (.not.rtdb_get(pars%rtdb,'gw:threshold',mt_dbl,1,pars%thresh))
     $  call errquit(pname//'failed to read threshold',0,RTDB_ERR)
c
      if (.not.rtdb_get(pars%rtdb,'gw:noqp',mt_int,1,pars%noqp))
     $  call errquit(pname//'failed to read noqp',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'gw:nvqp',mt_int,1,pars%nvqp))
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
      pars%tol2e = 10.0d0**(-itol2e)
c
c     RI basis set
c
      if (.not.bas_create(cd_bas_han,'ri basis'))
     $  call errquit(pname//'bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(pars%rtdb,geom,cd_bas_han,'ri basis'))
     $  call errquit(pname//'an "ri basis" is needed',0,BASIS_ERR)
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
     $    scftype,nbf_temp,ipol_temp,pars%nmo,2))
     $    call errquit(pname//'failed to read MO header',0,DISK_ERR)
      if (nbf_ao.ne.nbf_temp)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      if (ipol.ne.ipol_temp)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
c
      do i=1,ipol
        pars%g_movecs(i) = ga_create_atom_blocked(geom,
     $    ao_bas_han,'MO eigenvectors')
        if (.not.movecs_read(movecs_out,i,
     $      dbl_mb(pars%k_occ+(i-1)*nbf_ao),
     $      dbl_mb(pars%k_evals+(i-1)*nbf_ao),
     $      pars%g_movecs(i)))
     $    call errquit(pname//'failed to read MO vectors',0,DISK_ERR)
      enddo
c
      nocc = 0
      nvir = 0
      do i=1,pars%nmo(1)
        if (dbl_mb(pars%k_occ+i-1).gt.0.d0) then
          nocc = nocc + 1
          cycle
        else
          nvir = nvir + 1
        endif
      enddo
      pars%nocc(1) = nocc
      pars%nvir(1) = nvir
      pars%nov(1)  = nocc*nvir
      if (nocc+nvir.ne.pars%nmo(1))
     $  call errquit(pname//'something went wrong',0,0)
c
      if ((pars%noqp.gt.nocc).or.(pars%noqp.eq.-1)) pars%noqp = nocc
      if (pars%noqp.lt.-1)
     $  call errquit(pname//'invalid number noqp specified',0,0)
c
      if ((pars%nvqp.gt.nvir).or.(pars%nvqp.eq.-1)) pars%nvqp = nvir
      if (pars%nvqp.lt.-1)
     $  call errquit(pname//'invalid number nvqp specified',0,0)
c
      if (pars%evgw.or.pars%evgw0) then
        pars%noqp = nocc
        pars%nvqp = nvir
      endif
c
      if (pars%noqp+pars%nvqp.eq.0)
     $  call errquit(pname//'no quasiparticle energies to obtain',0,0)
c
      pars%nstates = nocc + pars%nvqp
c
      if (pars%evgw .or. pars%evgw0) then
        if (pars%eviter.lt.1)
     $    call errquit(pname//'invalid eviter',0,0)
      endif
c      
      if (pars%qpiter.lt.1)
     $  call errquit(pname//'invalid qpiter',0,0)
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
      call gw_memory(pars)
      call gw_ri_init(pars)
c      
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
        w(igl) = 2d0 * w(igl) / (1d0-x(igl))**2 / pi
        x(igl) = ((1d0 + x(igl))/(1d0 - x(igl)))**2
      enddo

      end subroutine

