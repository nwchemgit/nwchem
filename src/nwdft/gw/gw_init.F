      subroutine gw_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "case.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
      integer i
      integer nbf_ao, ipol, itol2e, nocc, nvir, g_2ceri, g_2cinv
c
      character(len=9), parameter :: pname = 'gw_init: '
      character(len=255) :: title, basisname, scftype
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
      if (.not.geom_create(pars%geom, 'geometry'))
     $  call errquit(pname//'geom_create failed',0, GEOM_ERR)
      if (.not.geom_rtdb_load(pars%rtdb, pars%geom, 'geometry'))
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
c     2-e integral tolerance
c
      if (.not.rtdb_get(pars%rtdb, 'dft:itol2e',mt_int,1,itol2e))
     $  call errquit(pname//'failed to read itol2e',0,RTDB_ERR)
      pars%tol2e = 10.0d0**(-itol2e)
c
c     Basis set
c
      if (.not.bas_create(pars%ao_bas_han,'ao basis'))
     $  call errquit(pname//'bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(pars%rtdb,pars%geom,pars%ao_bas_han,
     $  'ao basis'))
     $  call errquit(pname//'no ao basis found',0,BASIS_ERR)
      if (.not.int_normalize(pars%rtdb, pars%ao_bas_han))
     $  call errquit(pname//'int_normalize failed',0,BASIS_ERR)
      if (.not.bas_numbf(pars%ao_bas_han,pars%nbf_ao))
     $  call errquit(pname//'failed to read nbf_ao',0,BASIS_ERR)
      if (.not. bas_nprim_cn_max(pars%ao_bas_han,pars%nbf_ao_mxprim))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(pars%ao_bas_han,pars%nbf_ao_mxang))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(pars%ao_bas_han,pars%nbf_ao_mxcont))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(pars%ao_bas_han,pars%nbf_ao_mxnbf_cn))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(pars%ao_bas_han,pars%nbf_ao_mxnbf_ce))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(pars%ao_bas_han,pars%nshells_ao))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
c
c     RI basis set
c
      if (.not.bas_create(pars%ri_bas_han,'ri basis'))
     $  call errquit(pname//'bas_create failed',0,BASIS_ERR)
      if (.not.bas_rtdb_load(pars%rtdb,pars%geom,pars%ri_bas_han,
     $  'ri basis'))
     $  call errquit(pname//'an "ri basis" is needed',0,BASIS_ERR)
      call int_init(pars%rtdb, 1, pars%ri_bas_han)
      if (.not.int_norm_2c(pars%rtdb, pars%ri_bas_han))
     &  call errquit(pname//'int_norm_2c failed',0,INT_ERR)
      call int_terminate()
      if (.not. bas_numbf(pars%ri_bas_han, pars%nbf_ri))
     &  call errquit(pname//'basis set error',0,BASIS_ERR)
      if (.not. bas_nprim_cn_max(pars%ri_bas_han,pars%nbf_ri_mxprim))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(pars%ri_bas_han,pars%nbf_ri_mxang))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(pars%ri_bas_han,pars%nbf_ri_mxcont))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(pars%ri_bas_han,pars%nbf_ri_mxnbf_cn))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(pars%ri_bas_han,pars%nbf_ri_mxnbf_ce))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(pars%ri_bas_han,pars%nshells_ri))
     &  call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (bas_is_spherical(pars%ao_bas_han).and.
     &   (.not.bas_is_spherical(pars%ri_bas_han)))
     &        call int_app_set_no_texas(pars%rtdb)
c
      call scf_get_fock_param(pars%rtdb, pars%tol2e)
c
      call int_init(pars%rtdb, 2, (/pars%ao_bas_han, pars%ri_bas_han/))
      call print_integrals((/pars%ao_bas_han,pars%ri_bas_han/),.false.)
      call schwarz_init(pars%geom, pars%ao_bas_han)
c
c     Spin multiplicity, occupations
c
      if (.not.rtdb_get(pars%rtdb,'dft:ipol',mt_int,1,pars%ipol))
     $  call errquit(pname//'failed to read ipol',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb, 'dft:noc',mt_int,2,pars%nocc))
     $  call errquit(pname//'failed to read noc',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:mult',mt_int,1,pars%mult))
     $  call errquit(pname//'failed to read mult',0,RTDB_ERR)
c
c     MO vectors
c
      if (.not. rtdb_cget(pars%rtdb, 'dft:output vectors',1,
     $     pars%movecs_out))
     $     call errquit(pname//'no final MO found',0,RTDB_ERR)
c
c     Allocate arrays
c
      if (.not.ma_push_get(mt_dbl,pars%nbf_ao*pars%ipol,
     $    'eigenvalues', pars%l_evals, pars%k_evals))
     $  call errquit(pname//'failed to allocate evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,pars%nbf_ao*pars%ipol,
     $    'occupancies', pars%l_occ, pars%k_occ))
     $  call errquit(pname//'failed to allocate occ',0,MA_ERR) 
      if (.not.movecs_read_header(pars%movecs_out,title,basisname,
     $    scftype,nbf_ao,ipol,pars%nmo,2))
     $    call errquit(pname//'failed to read MO header',0,DISK_ERR)
      if (pars%nbf_ao.ne.nbf_ao)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      if (pars%ipol.ne.ipol)
     $  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
c
      do i=1,ipol
        pars%g_movecs(i) = ga_create_atom_blocked(pars%geom,
     $    pars%ao_bas_han,'MO eigenvectors')
        if (.not.movecs_read(pars%movecs_out,i,
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
      if (pars%noqp.gt.nocc)
     $  call errquit(pname//'invalid number noqp specified',0,0)
      if (pars%noqp .eq. -1) pars%noqp = nocc
c
      if (pars%nvqp.gt.nvir)
     $  call errquit(pname//'invalid number nvqp specified',0,0)
      if (pars%nvqp .eq. -1) pars%nvqp = nvir
c
      if (pars%eviter.lt.1)
     $  call errquit(pname//'invalid eviter',0,0)
      if (pars%qpiter.lt.1)
     $  call errquit(pname//'invalid qpiter',0,0)
c
c     3-center ERIs 
c
      call gw_memory(pars)
      call gw_coul_init(pars)
c
      return
      end subroutine gw_init
