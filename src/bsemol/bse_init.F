      subroutine bse_init(pars)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "bas.fh"
#include "case.fh"
#include "geom.fh"
#include "global.fh"
#include "bse.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"

#ifdef SCALAPACK
#include "scaleig.fh"
#endif

      type(bse_params_t) :: pars
      character(*), parameter :: pname = 'bse_init: '

      character(len=255) :: basisname, scftype, gwevals
      integer i, ilo, ihi, j, g_work
      integer nbf_temp, ipol_temp, nocc, nvir, nmo_temp(2)
      logical ldum, skipgw

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

      integer,parameter :: unitno = 66
      integer :: nmix(2)

      logical ncapr

#ifdef USE_OPENMP
      integer,external :: omp_get_max_threads, omp_get_thread_num
      integer,external :: omp_get_num_threads
#endif      


      ! Initialize GA variables
      pars%me = ga_nodeid()
      pars%nprocs = ga_nnodes()

      ! Initialize OMP variables
#ifdef USE_OPENMP
      pars%iMaxThreads = omp_get_max_threads()
#else
      pars%iMaxThreads = 1
#endif

      pars%npoles(:) = 0

      ! Get defaults in case input block is not present
      if (.not.rtdb_get(pars%rtdb,'bse:init',mt_log,1,ldum))
     &  call bse_defaults(pars%rtdb)

      ! Check method to solve BSE
      if (.not.rtdb_get(pars%rtdb,'bse:analytic',mt_log,1,
     &                  pars%analytic))
     &  call errquit(pname//'could not read from rtdb',1,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'bse:davidson',mt_log,1,
     &                  pars%davidson))
     &  call errquit(pname//'could not read from rtdb',2,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'bse:lanczos',mt_log,1,
     &                  pars%lanczos))
     &  call errquit(pname//'could not read from rtdb',3,RTDB_ERR)

      if ((pars%analytic.and.(pars%davidson.or.pars%lanczos)) .or.
     &    (pars%davidson.and.(pars%analytic.or.pars%lanczos)) .or.
     &    (pars%lanczos.and.(pars%analytic.or.pars%davidson)))
     &  call errquit(pname//'more than one method',4,INPUT_ERR)

      if (.not.rtdb_get(pars%rtdb, 'bse:skipgw', mt_log, 1, skipgw))
     &  call errquit(pname//'could not read skipgw',5, RTDB_ERR)

      ! Check wether we are computing singlet or triplet excitations
      if(.not.rtdb_get(pars%rtdb,'bse:singlet',mt_log,1,pars%singlet))
     &  call errquit(pname//'could not read from rtdb',1,RTDB_ERR)
      if(.not.rtdb_get(pars%rtdb,'bse:triplet',mt_log,1,pars%triplet))
     &  call errquit(pname//'could not read from rtdb',1,RTDB_ERR)
      if((pars%singlet .and. pars%triplet).and.(pars%analytic)) then
        call errquit(pname//'select either triplet or singlet',0,
     &               INPUT_ERR)
      endif   

      ! Check if we are doing Tamm-Dancoff approximation
      if(.not.rtdb_get(pars%rtdb,'bse:tda',mt_log,1,pars%tda))
     &  call errquit(pname//'could not read from rtdb',101,RTDB_ERR)

      ! Check if Output is restricted to an energy window
      if(.not.rtdb_get(pars%rtdb,'bse:dowindow',mt_log,1,pars%dowindow))
     &  call errquit(pname//'could not read from rtdb',105,RTDB_ERR)
      if(pars%dowindow) then
        if(.not.rtdb_get(pars%rtdb,'bse:elower',mt_dbl,1,pars%elower))
     &    call errquit(pname//'could not read from rtdb',108,RTDB_ERR)
        if(.not.rtdb_get(pars%rtdb,'bse:eupper',mt_dbl,1,pars%eupper))
     &    call errquit(pname//'could not read from rtdb',110,RTDB_ERR)
      endif

      ! Number of roots
      if(.not.rtdb_get(pars%rtdb,'bse:nroots',mt_int,1,pars%nroots))
     &  call errquit(pname//'could not read from rtdb',115,RTDB_ERR)

      ! Starting space size
      if(.not.rtdb_get(pars%rtdb, 'bse:nspace', mt_int, 1, pars%nspace))
     &  call errquit(pname//'could not read from rtdb', 119, RTDB_ERR) 

      ! Maximum Number of roots
      if(.not.rtdb_get(pars%rtdb,'bse:nmax',mt_int,1,pars%nmax))
     &  call errquit(pname//'could not read from rtdb',123,RTDB_ERR)
      if ((pars%nmax.le.5*pars%nroots).and.pars%davidson) then
        call errquit(pname//'maxroots should be at least 5*nroots',
     &               126,INPUT_ERR)
      endif
      if((pars%nroots*pars%nspace.gt.pars%nmax).and.pars%davidson)
     &  call errquit(pname//'maxroots should be at least nroots*nspace',
     &               130, INPUT_ERR) 

      ! Maximum number of iterations
      if(.not.rtdb_get(pars%rtdb,'bse:maxiter',mt_int,1,pars%maxiter))
     &  call errquit(pname//'could not read from rtdb',134,RTDB_ERR)
      ! Get some info from DFT module
      call dft_rdinput(pars%rtdb)
      if (cdfit) then
        if (.not.bas_destroy(cd_bas_han))      
     $    call errquit(pname//'failed to destroy cd_bas',0,BASIS_ERR)
      endif
      if (.not.rtdb_get(pars%rtdb, 'dft:itol2e',mt_int,1,itol2e))
     $  call errquit(pname//'failed to read itol2e',0,RTDB_ERR)
      pars%tol2e = 10.0d0**(-itol2e-2)

      ! Fitting basis
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
      pars%nri = nbf_cd

      !Spin multiplicity, occupations
      if (.not.rtdb_get(pars%rtdb,'dft:ipol',mt_int,1,ipol))
     $  call errquit(pname//'failed to read ipol',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb, 'dft:noc',mt_int,2,noc))
     $  call errquit(pname//'failed to read noc',0,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'dft:mult',mt_int,1,mult))
     $  call errquit(pname//'failed to read mult',0,RTDB_ERR)
      pars%ipol = ipol

      ! MO vectors
      if (.not. rtdb_cget(pars%rtdb, 'dft:output vectors',1,movecs_out))
     $     call errquit(pname//'no final MO found',0,RTDB_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,
     $    'mf evals', pars%l_mf_evals, pars%k_mf_evals))
     $  call errquit(pname//'failed to allocate evals',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao*ipol,
     $    'gw evals', pars%l_gw_evals, pars%k_gw_evals))
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
      do i=1,ipol
        pars%g_movecs(i) = ga_create_atom_blocked(geom,
     $    ao_bas_han,'MO eigenvectors')
        if (.not.movecs_read(movecs_out,i,
     $      dbl_mb(pars%k_occ+(i-1)*nbf_ao),
     $      dbl_mb(pars%k_mf_evals+(i-1)*nbf_ao),
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
        pars%npoles(i) = pars%nocc(i)*pars%nvir(i)
      enddo

      if ((pars%nroots.eq.-1).or.(pars%nroots.gt.sum(pars%npoles))) then
        pars%nroots = sum(pars%npoles)
      endif

      if (pars%nmax.gt.sum(pars%npoles)) then
        pars%nmax = sum(pars%npoles)
      endif


      ! Obtain NCAP DD shift
 1001 format(2X,A12,F8.3,A3)
      if (.not.rtdb_get(pars%rtdb,'bse:ncap',mt_log,1,pars%ncap))
     &  call errquit(pname//'rtdb_get failed',176,RTDB_ERR)
      if (.not.rtdb_get(pars%rtdb,'bse:ncapr',mt_log,1,ncapr))
     &  call errquit(pname//'rtdb_get failed',176,RTDB_ERR)

      pars%ncap = ncapr .or. pars%ncap

      if (pars%ncap) then
        zeta = 0.304121d0
        if (abs(xfac(87)).gt.1d-10) zeta=0.5d0
        aq2 = 2d0/9d0 * (1d0-zeta)**2
        do i=1,ipol
          ehomo = dbl_mb(pars%k_mf_evals+(i-1)*nbf_ao+pars%nocc(i)-1)
          vxddp = -0.5d0*aq2*(1d0 + sqrt(1d0-4d0*ehomo/aq2))!*0.85d0
          vxddm = -0.5d0*aq2*(1d0 - sqrt(1d0-4d0*ehomo/aq2))!*0.85d0
          pars%vxddp(i) = vxddp
          pars%vxddm(i) = vxddm
          if (pars%me.eq.0) then
            write(luout,*)    ' Shifting eigvals by 100% NCAP DD'
            write(luout,1001) ' v_x^{DD-}: ',vxddp*ha2ev, ' eV'
            write(luout,1001) ' v_x^(DD+): ',vxddm*ha2ev, ' eV'
          endif
          call daxpy(pars%nocc(i),1d0,vxddp,0,
     &       dbl_mb(pars%k_mf_evals+(i-1)*nbf_ao),1)
          call daxpy(pars%nvir(i),1d0,vxddm,0,
     &       dbl_mb(pars%k_mf_evals+(i-1)*nbf_ao+pars%nocc(i)),1)
        enddo
      endif

      ! Read GW QP energies
      !if (.not.pars%ncap) then
      if(skipgw) then
        call dcopy(nbf_ao*pars%ipol,dbl_mb(pars%k_mf_evals),1,
     &             dbl_mb(pars%k_gw_evals),1)
      else   
        call util_file_name('gwevals',.false.,.false.,gwevals)
        call util_file_name_resolve(gwevals,.false.)
        if (pars%me.eq.0) then
          open(unit=unitno,status='old',form='unformatted',file=gwevals)
          do i=1,pars%ipol
            call sread(unitno, dbl_mb(pars%k_gw_evals+(i-1)*nbf_ao),
     &                 pars%nmo)     
          enddo
          close(unit=unitno)
        endif
        call ga_brdcst(mt_dbl,dbl_mb(pars%k_gw_evals),
     &                 ma_sizeof(mt_dbl,nbf_ao*pars%ipol,mt_byte),0)
      endif  

      !3-center ERIs 
      call int_init(pars%rtdb, 2, (/ao_bas_han, cd_bas_han/))
      call print_integrals((/ao_bas_han,cd_bas_han/),.false.)
      call schwarz_init(geom, ao_bas_han)
      call scf_get_fock_param(pars%rtdb, pars%tol2e)

      ! Obtain irreps information
      do i=1,ipol
        if(.not.ma_push_get(mt_int,nbf_ao,'irreps',pars%l_irs(i),
     $          pars%k_irs(i))) call errquit
     $          ('bse_init: failed to allocate irs',0,MA_ERR)
        if(.not.ga_create(mt_dbl,nbf_ao,pars%nmo,'Work',nbf_ao,-1,
     $                    g_work))
     $    call errquit('bse_init: failed to allocate work',0,GA_ERR)
        call ga_copy_patch('n',pars%g_movecs(i),1,nbf_ao,1,pars%nmo,
     $                      g_work,1,nbf_ao,1,pars%nmo)
        call sym_movecs_adapt(ao_bas_han, 1.0d-8, g_work,
     $                        int_mb(pars%k_irs(i)),nmix(i))
        if(.not.ga_destroy(g_work)) call errquit
     $    ('bse_init: failed to destroy g_work',0, GA_ERR)
      enddo

      ! Get dipole in CIS basis
      call bse_dipole(pars)

      !Clean previous grid
      call grid_cleanup(.false.)

      !Print expected memory requirements

      !Compute three-center integrals
      call bse_ri_init(pars)
     
      end subroutine bse_init
