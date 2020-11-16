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
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      integer i
      integer nbf_temp, ipol_temp, nocc, nvir
c
      character(len=9), parameter :: pname = 'gw_init: '
      character(len=255) :: basisname, scftype
c
      logical int_normalize, int_norm_2c
      external int_normalize, int_norm_2c
c
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      logical ga_create_atom_blocked
      external ga_create_atom_blocked
c
c     -----------------------------------------------------------------
c
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
      if (.not.rtdb_get(pars%rtdb,'gw:acgw',mt_log,1,pars%cdgw))
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
      if (.not.rtdb_get(pars%rtdb,'gw:noqp',mt_int,1,pars%nvqp))
     $  call errquit(pname//'failed to read nvqp',0,RTDB_ERR)
c
      if (.not.rtdb_get(pars%rtdb,'gw:analyticw',mt_log,1,pars%anaw))
     $  call errquit(pname//'failed to read analyticw',0,RTDB_ERR)
c
c     Geometry
c
      if (.not.geom_create(geom, 'geometry'))
     $  call errquit(pname//'geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(pars%rtdb, geom, 'geometry'))
     $  call errquit(pname//'geom_rtdb_load failed',0,RTDB_ERR)
c
c     Coulomb Attenuation Method
c
      if (.not.rtdb_get(pars%rtdb, 'dft:cam_exch',mt_log,1,cam_exch))
     $  cam_exch = .false.
      if (.not.rtdb_get(pars%rtdb, 'dft:cam_omega',mt_dbl,1,cam_omega))
     $  cam_omega = 0.0d0
      if (.not.rtdb_get(pars%rtdb, 'dft:cam_alpha',mt_dbl,1,cam_omega))
     $  cam_alpha = 0.0d0
      if (.not.rtdb_get(pars%rtdb, 'dft:cam_beta',mt_dbl,1,cam_omega))
     $  cam_beta = 0.0d0
c
c     XC functionals
c
      if (.not.rtdb_get(pars%rtdb,'dft:xfac',mt_dbl,numfunc,xfac))
     1  call errquit(pname//'failed to read xfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:cfac',mt_dbl,numfunc,cfac))
     1  call errquit(pname//'failed to read cfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:lxfac',mt_log,numfunc,lxfac))
     1  call errquit(pname//'failed to read lxfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:nlxfac',mt_log,numfunc,nlxfac))
     1  call errquit(pname//'failed to read nlxfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:lcfac',mt_log,numfunc,lcfac))
     1  call errquit(pname//'failed to read lcfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:nlcfac',mt_log,numfunc,nlcfac))
     1  call errquit(pname//'failed to read nlcfac',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:xccomb',mt_log,numfunc,xccomb))
     1  call errquit(pname//'failed to read xccomb',0, RTDB_ERR)
      if (.not.rtdb_cget(pars%rtdb,'dft:cname',numfunc,cname))
     1  call errquit(pname//'failed to read xname',0, RTDB_ERR)
      if (.not.rtdb_cget(pars%rtdb,'dft:xcname',numfunc,xcname))
     1  call errquit(pname//'failed to read cname',0, RTDB_ERR)
      if (.not.rtdb_cget(pars%rtdb,'dft:xname',numfunc,xname))
     1  call errquit(pname//'failed to read xcname',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:lb94',mt_log,1,lb94))
     1  call errquit(pname//'failed to read lb94',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:cs00',mt_log,1,cs00))
     1  call errquit(pname//'failed to read cs00',0, RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:ncap',mt_log,1,ncap))
     1  call errquit(pname//'failed to read ncap',0, RTDB_ERR)
c
c     2-e integral tolerance
c
      if (.not.rtdb_get(pars%rtdb, 'dft:itol2e',mt_int,1,itol2e))
     $  call errquit(pname//'failed to read itol2e',0,RTDB_ERR)
      pars%tol2e = 10.0d0**(-itol2e)
c
c     Basis set
c
      if (.not.bas_create(ao_bas_han,'ao basis'))
     $  call errquit(pname//'bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(pars%rtdb,geom,ao_bas_han,'ao basis'))
     $  call errquit(pname//'no ao basis found',0,BASIS_ERR)
      if (.not.int_normalize(pars%rtdb, ao_bas_han))
     $  call errquit(pname//'int_normalize failed',0,BASIS_ERR)
      if (.not.bas_numbf(ao_bas_han,nbf_ao))
     $  call errquit(pname//'failed to read nbf_ao',0,BASIS_ERR)
      if (.not. bas_nprim_cn_max(ao_bas_han,nbf_ao_mxprim))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(ao_bas_han,nbf_ao_mxang))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(ao_bas_han,nbf_ao_mxcont))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(ao_bas_han,nbf_ao_mxnbf_cn))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(ao_bas_han,nbf_ao_mxnbf_ce))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(ao_bas_han,nshells_ao))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
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
      call scf_get_fock_param(pars%rtdb, pars%tol2e)
c
      call int_init(pars%rtdb, 2, (/ao_bas_han, cd_bas_han/))
      call print_integrals((/ao_bas_han,cd_bas_han/),.false.)
      call schwarz_init(geom, ao_bas_han)
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
      if (pars%noqp+pars%nvqp.eq.0)
     $  call errquit(pname//'no quasiparticle energies to obtain',0,0)
c
      if (pars%eviter.lt.1)
     $  call errquit(pname//'invalid eviter',0,0)
      if (pars%qpiter.lt.1)
     $  call errquit(pname//'invalid qpiter',0,0)
c
c     3-center ERIs 
c
      call gw_memory(pars)
      call gw_ri_init(pars)
c
      return
      end subroutine gw_init
