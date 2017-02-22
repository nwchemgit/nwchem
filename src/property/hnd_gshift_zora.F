      subroutine hnd_gshift_zora(rtdb,basis,geom)
c Modified from hnd_giaox by F. Aquino 11-17-10
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"
#include "zora.fh" 
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer nclosed(2), nocc(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, ic, iatom, iocc, ifld, ioff
      integer sh_atom1
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer l_dia, k_dia, l_para, k_para
      integer l_xyz, k_xyz, l_zan, k_zan
      integer l_tmp, k_tmp
      integer g_dens(3),g_h01,g_d1,
     &        g_rhs,g_rhs0,g_fock
      integer g_fock_Coul,g_fock_Exch,
     &        g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_eSji,g_rhs_1e
      integer g_d1_oo,g_d1_ov,type_NMR

      integer g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji,g_d1_ov_ExchSFB
      integer vectors(2), geomnew, i, j, ij
      double precision atn, tol2e, val, isotr, aniso
      double precision val_oo,val_ov

      double precision jfac(12),kfac(12),a(6),axs(3,3),eig(3),xfac
      character*255 zorafilename
      integer type_nmrdata,g_AtNr1,icalczora,icalczora_AB
      integer ga_dia_epr,ga_para1_epr,ga_Fji,ga_h01_epr
      integer ga_dia_epr_AB,ga_para1_epr_AB ! for A,B tensor contrib.
      logical dft_zoraNMR_read,dft_zoraNMR_read_AB
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      integer ld(2),cbuf,ind,cbuf1,cbuf2
      logical  cphf2, file_write_ga, file_read_ga, cphf,
     &         switch_gshift_analysis,switch_skip_cphf
      double precision coeff_exch,ac_occ(2)
c ----- Definitions for NLMO analysis ---- START
      integer i1,j1,acc_vec,l_tvec,k_tvec,gshiftfile          
      integer g_tvec ! for nbo analysis
      external create_munu4nbo_gshift
c ----- Definitions for NLMO analysis ---- END      
      logical status
      double precision ppt,coeffpol
      double precision par_arr(9,3,3) 
      double precision val_ov_Coul,val_ov_Exch,val_ov_noJK,
     &                 val_ov_1e,val_ov_eSji

      data ppt     /1.0d+03/  
      data tol2e   /1.0d-10/
      integer npol,ntot,ispin,nind_jk,cdens
      external cphf2, file_write_ga, file_read_ga, cphf
      external giao_aotomo,mat_transpose,shell_fock_build 
      external get_prelim_fock,get_P10,get_P10_1,
     &         add_H10,add_fock,
     &         skip_cphf,skip_cphf_JK,
     &         new_giao_2e,new_giao_2e_JK,
     &         dft_zoraNMR_read,dft_zoraNMR_read_AB,
     &         get_gshift_tensor_AorB,
     &         dft_zoraCPHF_write,dft_zoraCPHF_read
      integer debug_gshift
      integer shift,disp
      logical skip_cphf_ev_gshift
c
      debug_gshift=0 ! =1 for debugging print outs of matrices
c     switch_skip_cphf=.true. ! For skipping cphf or cpks
      switch_skip_cphf=.false. ! For NOT skipping cphf or cpks
      switch_gshift_analysis=.false. ! For default mode gshift tensor
c     switch_gshift_analysis=.true.  ! For analysis of gshift tensor
c
      npol=2 ! gshift works ONLY for npol=2
      if (ga_nodeid().eq.0) then
        write(LuOut,*)
        call util_print_centered(LuOut,
     $   'g-Shift Tensor (in ppm)', 23, .true.)
        write(LuOut,*)
      end if
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb   
c
c     If DFT get part of the exact exchange defined
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
      nind_jk=12
c ------- using different settings (j,k)  03-14-11--- START
      do ifld = 1,nind_jk
        jfac(ifld) =  0.0d0       ! used in shell_fock_build()
        kfac(ifld) = -1.0d0*xfac  ! used in shell_fock_build()
        if (ga_nodeid().eq.0.and.debug_gshift.eq.1) then
          write(*,145) ifld,jfac(ifld),kfac(ifld)
 145   format('(j,k)(',i3,')=(',f15.8,',',f15.8,')')
        endif
      enddo
c     Integral initialization
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
c     Find out from rtdb which atoms we need to calculate shielding for
c     Get number of atoms (all or number from rtdb)
c     Get which atoms (all or some read from rtdb)
c     Allocate arrays which will hold atomic information (k_zan and k_xyz)
      status = rtdb_parallel(.true.)
c      if (.not.rtdb_get(rtdb, 'giao:natoms', MT_INT, 1,sh_atom)) then  
c         if (.not.geom_ncent(geom,sh_atom)) call                   
c     &       errquit('hnd_giaox: geom_ncent failed ',0, GEOM_ERR)
c      endif
         if (.not.geom_ncent(geom,sh_atom1)) call                   
     &       errquit('hnd_giaox: geom_ncent failed ',0, GEOM_ERR)
      sh_atom=1 ! for EPR
c     if (ga_nodeid().eq.0) then
c       write(*,116) sh_atom1,sh_atom
c116    format('(sh_atom1,sh_atom)=(',i3,',',i3,')')
c     endif

      if (.not. ma_push_get(mt_int,sh_atom1,'nmr tmp',l_tmp,k_tmp)) 
     &    call errquit('hnd_gshift: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,3*sh_atom1,'nmr at',l_xyz,k_xyz)) 
     &    call errquit('hnd_gshift: ma_push_get failed k_xyz',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,sh_atom1,'nmr zan',l_zan,k_zan)) 
     &    call errquit('hnd_gshift: ma_push_get failed k_zan',0,MA_ERR)

      gshiftfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
      if (gshiftfile.eq.1) then ! ------- g-shift-if++++ START
c  Allocate memory for l_tvec,k_tvec --- start 
       if (.not.ma_alloc_get(mt_dbl,3*3,'tvec',l_tvec,k_tvec)) 
     &     call errquit('hnd_gshift: ma failed',0,MA_ERR) 
       call dcopy(3*3,0.0d0,0,dbl_mb(k_tvec),1) ! reset 
c  Allocate memory for l_tvec,k_tvec --- end
      endif ! --------------------------- g-shift-if++++ END

c     By default set array from 1 to sh_atom (which is all when giao:natoms is not there)
      do ixy = 1, sh_atom1
         int_mb(k_tmp+ixy-1) = ixy
      enddo
c     Try to read the atom list from rtdb. If it is not there, we still have the default list
c     status = rtdb_get(rtdb,'giao:atom list',mt_int,sh_atom,
c     &                  int_mb(k_tmp))
      do ixy = 0, sh_atom1-1
         if (.not. geom_cent_get(geom, int_mb(k_tmp+ixy), tag, 
     &                           dbl_mb(k_xyz+3*ixy),dbl_mb(k_zan+ixy)))
     &   call errquit('hnd_gshift: geom_cent_tag failed',0,GEOM_ERR)
c       if (ga_nodeid().eq.0) then
c        write(*,119) ixy,dbl_mb(k_xyz+3*ixy),dbl_mb(k_xyz+3*ixy+1),
c    &                dbl_mb(k_xyz+3*ixy+2),
c    &                dbl_mb(k_zan+ixy)
c119     format('atom(',i3,')=(',f15.8,',',f15.8,',',
c    &          f15.8,',',f15.8,')')
c       endif
      enddo 
c
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_gshift: could not get nbf',0, BASIS_ERR)
c ++++++ Reading gshift data from file ++++++ START
c       Note.- lbl_nmrgshift defined in zora.fh
        call util_file_name(lbl_nmrgshift,.false.,.false.,
     &                      zorafilename)
      icalczora = 0  ! initialize the flag
      type_nmrdata=3 ! =1,2,3=shieldings,hyperfine,gshift      
      if (.not.dft_zoraNMR_read(
     &     zorafilename,
     &     type_nmrdata,
     &     nbf,sh_atom,
     &     g_AtNr1,
     &     ga_dia_epr,
     &     ga_para1_epr,
     &     ga_h01_epr,
     &     ga_Fji)) icalczora=1 
c Note.- If I print the GAs here it gets freezed
c ++++++ Reading gshift data from file ++++++ END
c ++++++ Reading gshift-AB data from file ++++++ START
c       Note.- lbl_nmrgshift defined in zora.fh
        call util_file_name(lbl_nmrgshift_AB,.false.,.false.,
     &                      zorafilename)
      icalczora_AB = 0  ! initialize the flag
      type_nmrdata=3 ! =1,2,3=shieldings,hyperfine,gshift      
      if (.not.dft_zoraNMR_read_AB(
     &     zorafilename,
     &     type_nmrdata,
     &     nbf,
     &     sh_atom,
     &     g_AtNr1,
     &     ga_dia_epr_AB,
     &     ga_para1_epr_AB)) icalczora_AB=1 
c Note.- If I print the GAs here it gets freezed
c ++++++ Reading gshift-AB data from file ++++++ END

      if (.not. ma_push_get(mt_dbl,npol*nbf,'occ num',
     &                      l_occ,k_occ)) 
     &   call errquit('hnd_gshift: ma_push_get failed k_occ',
     &            0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,npol*nbf,'eigenval',
     &                      l_eval,k_eval)) 
     &   call errquit('hnd_gshift: ma_push_get failed k_eval',
     &             0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nocc,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
      if (ga_nodeid().eq.0) then
        write(*,10) nocc(1),nocc(2),
     &              nclosed(1),nclosed(2),
     &              nvirt(1),nvirt(2),scftyp
 10     format('nocc=(',i3,',',i3,') ',
     &         'nclos=(',i3,',',i3,') ',
     &         'nvirt=(',i3,',',i3,') scftyp=',a,')')
      endif

c ---- calculate total occupation alpha and beta --- START
       do i=1,npol
        val=0.0d0
        ind=nbf*(i-1)
        do j=1,nocc(i)         
         val=val+dbl_mb(k_occ+ind+j-1)
        enddo
        ac_occ(i)=val
       enddo
c      if (ga_nodeid().eq.0) then
c       write(*,2) ac_occ(1),ac_occ(2),nocc(1),nocc(2)
c2      format('In hnd_gshift_zora:ac_occ=(',
c    &         f15.8,',',f15.8,')  nocc=(',
c    &         i4,',',i4,')')
c      endif
c ---- calculate total occupation alpha and beta --- END
c      if (noc(1) .ne. noc(2)) then
      if (ac_occ(1) .ne. ac_occ(2)) then
c      coeffpol=4.0d0/(noc(1)-noc(2))
       coeffpol=4.0d0/(ac_occ(1)-ac_occ(2))
      else
       write(*,1) nocc(1),nocc(2)
 1     format('Error in hnd_gshift_zora(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       stop
      endif 
      if (ga_nodeid().eq.0) 
     &  write(*,*) 'coeffpol=',coeffpol

c ------ Store nopen in rtdb so that CPHF routine is happy ---- START
          if (.not. rtdb_put(rtdb, 'scf:nopen', 
     &         MT_INT, 1, nocc(1)-nocc(2)))
     *         call errquit('scfinit:rtdbput nopen failed',
     &         nocc(1)-nocc(2),
     &       RTDB_ERR)
c ------ Store nopen in rtdb so that CPHF routine is happy ---- END
c
c     Get Unperturbed Density Matrix
      call hnd_prp_get_dens(rtdb,geom,basis,
     &                      g_dens,ndens, ! out
     &                      scftyp,nclosed,nocc,nvirt)
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
      enddo
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs))
     &   call errquit('hnd_gshift: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)

c ----- TEST: for testing Coulomb and Exchange contrib --- START
      if (switch_gshift_analysis) then
       if(.not.ga_create(MT_DBL,ntot,3,'rhs_Coul',-1,-1,g_rhs_Coul))
     &   call errquit('hnd_gshift: ga_create failed g_rhsJ',0,GA_ERR)
       call ga_zero(g_rhs_Coul)
       if(.not.ga_create(MT_DBL,ntot,3,'rhs_Exch',-1,-1,g_rhs_Exch))
     &   call errquit('hnd_gshift: ga_create failed g_rhsK',0,GA_ERR)
       call ga_zero(g_rhs_Exch)
       if(.not.ga_create(MT_DBL,ntot,3,'rhs_noJK',-1,-1,g_rhs_noJK))
     &   call errquit('hnd_gshift: ga_create failed g_rhsnJK',0,GA_ERR)
       call ga_zero(g_rhs_noJK)
       if(.not.ga_create(MT_DBL,ntot,3,'rhs_eSji',-1,-1,g_rhs_eSji))
     &   call errquit('hnd_gshift: ga_create failed g_rhseSji',0,GA_ERR)
       call ga_zero(g_rhs_eSji)
       if(.not.ga_create(MT_DBL,ntot,3,'rhs_1e',-1,-1,g_rhs_1e))
     &   call errquit('hnd_gshift: ga_create failed g_rhs1e',0,GA_ERR)
       call ga_zero(g_rhs_1e)
      endif
c ----- TEST: for testing Coulomb and Exchange contrib --- END

      if (debug_gshift.eq.1.and.ga_nodeid().eq.0) then
      write(*,*) 'after-creating g_rhs ...'
      write(*,102) npol,nocc(1),nocc(2),
     &            nclosed(1),nclosed(2),
     &            nvirt(1),nvirt(2),scftyp
 102   format('BEF pre-fock::npol=',i3,' nocc=(',i3,',',i3,') ',
     &       'nclos=(',i3,',',i3,') ',
     &       'nvirt=(',i3,',',i3,') scftyp=',a,')')
      endif

      if (switch_gshift_analysis) then
       call get_prelim_fock_debug(
     &                      g_d1,  ! out: 
     &                      g_rhs, ! out: rhs expression
     &                      g_rhs0,     ! out: used in get_d1()
     &                      g_rhs_eSji, ! out
     &                      vectors,dbl_mb(k_eval),
     &                      dbl_mb(k_xyz),sh_atom1,basis, 
     &                      nbf,nmo,npol,nocc,nvirt) 
      else
       call get_prelim_fock(g_d1,  ! out: 
     &                      g_rhs, ! out: rhs expression
     &                      g_rhs0,! out: used in get_d1()
     &                      vectors,dbl_mb(k_eval),
     &                      dbl_mb(k_xyz),sh_atom1,basis, 
     &                      nbf,nmo,npol,nocc,nvirt) 
      endif

      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-aft-get_prelim_fock -------- START'
       call ga_print(g_rhs)
       call ga_print(g_rhs0)
       call ga_print(g_d1)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-aft-get_prelim_fock -------- START'
        write(*,103) npol,nocc(1),nocc(2),
     &               nclosed(1),nclosed(2),
     &               nvirt(1),nvirt(2),scftyp
 103    format('AFT pre-fock:: npol=',i3,' nocc=(',i3,',',i3,') ',
     &         'nclos=(',i3,',',i3,') ',
     &         'nvirt=(',i3,',',i3,') scftyp=',a,')')
      endif 
c
c     Build "first order fock matrix"
      if (use_theory.eq.'dft') then
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
      endif
      clo(1) = 3*npol*2
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',chi,g_fock)) call 
     &    errquit('hnd_gshift: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)

c     if (ga_nodeid().eq.0) then
c        write(*,101) npol,nocc(1),nocc(2),
c    &                nclosed(1),nclosed(2),
c    &                nvirt(1),nvirt(2),scftyp
c101     format('BEF. shell_fock_build: npol=',i3,
c    &          ' nocc=(',i3,',',i3,') ',
c    &          'nclos=(',i3,',',i3,') ',
c    &          'nvirt=(',i3,',',i3,') scftyp=',a,')')
c     endif

c
c     Note: Just the exchange: jfac = 0.d0 (see above)

      if (.not.cam_exch) then
         call shell_fock_build(geom, basis,0,3*npol*2,
     $                         jfac,kfac,tol2e,
     &                         g_d1,  ! input
     &                         g_fock,! output
     &                         .false.)
      else
         call shell_fock_build_cam(geom, basis,0,3*npol*2,
     $                             jfac,kfac,tol2e,
     &                             g_d1,  ! input
     &                             g_fock,! output
     &                             .false.)
      end if
       if (debug_gshift.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '----- AFT shell_fock_build():g_fock ----- START'
         call ga_print(g_fock)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '----- AFT shell_fock_build():g_fock ----- END'
       endif
c      if (ga_nodeid().eq.0) then
c        write(*,104) npol,nocc(1),nocc(2),
c    &                nclosed(1),nclosed(2),
c    &                nvirt(1),nvirt(2),scftyp
c104     format('AFT. shell_fock_build: npol=',i3,
c    &          ' nocc=(',i3,',',i3,') ',
c    &          'nclos=(',i3,',',i3,') ',
c    &          'nvirt=(',i3,',',i3,') scftyp=',a,')')
c      endif

      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0) then
        write(*,*) '---- g_rhs-bef-add_fock -------- START'
        write(*,*) 'Bef add_fock()'
       endif
        call ga_print(g_rhs)
        call ga_print(g_fock)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-bef-add_fock -------- END'
      endif
      call add_fock(g_rhs, ! out: g_rhs=g_rhs+g_fock
     &              g_fock,vectors,
     &              nbf,nmo,npol,nocc,nvirt) 

       if (debug_gshift.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,*) 'Aft add_fock()'
         write(*,*) '---- g_rhs-aft-add_fock -------- START'
        endif
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0)
     &    write(*,*) '---- g_rhs-aft-add_fock -------- END'
       endif

c     Cleanup of g_d1 and g_fock, not needed for now
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('hnd_gshift: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_gshift: ga_destroy failed g_fock',0,GA_ERR)

      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-bef-add_H10 -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-bef-add_H10 -------- END'
      endif
      if (switch_gshift_analysis) then
       call ga_zero(g_rhs_1e)
       call add_H10_debug(
     &             g_rhs, ! out: ga_rhs(a,i)=ga_rhs(a,i)+H10(a,i)
     &             g_rhs_1e, ! out
     &             ga_Fji,vectors, 
     &             dbl_mb(k_xyz),sh_atom1,basis, 
     &             nbf,nmo,npol,nocc,nvirt,do_zora,rtdb) 
      else
       call add_H10(g_rhs, ! out: ga_rhs(a,i)=ga_rhs(a,i)+H10(a,i)
     &             ga_Fji,vectors, 
     &             dbl_mb(k_xyz),sh_atom1,basis, 
     &             nbf,nmo,npol,nocc,nvirt,do_zora,rtdb) 
      endif 

       if (debug_gshift.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,*) 'Aft. add_H10()'
         write(*,*) '---- g_rhs-aft-add_H10 -------- START'
        endif
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-aft-add_H10 -------- END'
       endif

c     Remaining term is Perturbed (GIAO) two-electron term times Unperturbed density
c     Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in AO basis
      alo(1) = -1 
      alo(2) = -1
      alo(3) =  1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*npol
      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
     &    errquit('hnd_gshift: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
c -- TEST: create g_fock_Coul,g_fock_Exch ----------- START
      if (switch_gshift_analysis) then
       if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix g_fock_Coul',
     &    alo,g_fock_Coul)) call 
     &    errquit('hnd_gshift: nga_create failed g_fockJ',0,GA_ERR)
       call ga_zero(g_fock_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix g_fock_Exch',
     &    alo,g_fock_Exch)) call 
     &    errquit('hnd_gshift: nga_create failed g_fockX',0,GA_ERR)
       call ga_zero(g_fock_Exch)
      endif
c -- TEST: create g_fock_Coul,g_fock_Exch ----------- END

      if(use_theory.eq.'dft') then
       ifld = 4
         if (.not. rtdb_put(rtdb,'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('hnd_gshift: rtdb_put failed',0,RTDB_ERR)

      endif

       if (debug_gshift.eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_dens bef new_giao_2e  -------- START'
         do cdens=1,npol
          call ga_print(g_dens(cdens))
         enddo
        if (ga_nodeid().eq.0)
     &    write(*,*) '---- g_dens bef new_giao_2e  -------- END'    
       endif

      if (switch_gshift_analysis) then
c +++++ Test: getting Coulomb and Exchange contrib separate --- START
       call new_giao_2e_JK(geom,basis,nbf,tol2e,
     &                 g_dens, !  in: e-denstiy 
     &                 g_fock, ! out: fock matrix
     &                 g_fock_Coul,
     &                 g_fock_Exch,
     &                 xfac,
     &                 npol)
c +++++ Test: getting Coulomb and Exchange contrib separate --- END
      else
       call new_giao_2e(geom,basis,nbf,tol2e,
     &                 g_dens, !  in: e-denstiy 
     &                 g_fock, ! out: fock matrix
     &                 xfac,
     &                 npol)
      endif
       if (debug_gshift.eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-new_giao -------- START'
         call ga_print(g_fock)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-new_giao -------- END'
       endif
c ++++++++ calling new_giao_2e() ++++++
c +++++++++++++++++++++++++++++++++++++

      if(use_theory.eq.'dft') then
      ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('hnd_gshift: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('hnd_gshift: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      endif
c
c     Transform to MO basis and add to right-hand-side
      call giao_aotomo(g_fock,vectors,nocc,nvirt,npol,3,nbf)
c ----- TEST: AO2MO for g_fock_Coul,g_fock_Exch ------ START
      if (switch_gshift_analysis) then
       call giao_aotomo(g_fock_Coul,vectors,nocc,nvirt,npol,3,nbf)
       call giao_aotomo(g_fock_Exch,vectors,nocc,nvirt,npol,3,nbf)
      endif
c ----- TEST: AO2MO for g_fock_Coul,g_fock_Exch ------ END

      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0) then
        write(*,*) 'Aft add_fock1-giao-aotomo()'
        write(*,*) '---- g_fock-aft-giao_aotomo -------- START'
       endif
       call ga_print(g_fock)
       if (ga_nodeid().eq.0)        
     &   write(*,*) '---- g_fock-aft-giao_aotomo -------- END'
      endif
c ------- Store g_rhs without Coulomb and Exchange ---- START
      if (switch_gshift_analysis) then
       call ga_copy(g_rhs,g_rhs_noJK)
       do ispin=1,npol
c ------ definitions for g_rhs -------- START
        disp = nocc(1)*nvirt(1)*(ispin-1)
        shift=3*(ispin-1)
        blo(1) = disp+1
        bhi(1) = disp+nocc(ispin)*nvirt(ispin)
        blo(2) = 1
        bhi(2) = 3
c ------ definitions for g_rhs -------- END
        alo(1) = nocc(ispin)+1
        ahi(1) = nmo
        alo(2) = 1
        ahi(2) = nocc(ispin)
        alo(3) = shift+1
        ahi(3) = shift+3
        if      (npol.eq.1) then
         call nga_scale_patch(g_rhs_noJK,blo,bhi,-4.0d0)
         call nga_scale_patch(g_rhs_1e  ,blo,bhi,-4.0d0)
         call nga_scale_patch(g_rhs_eSji,blo,bhi,-4.0d0)
        else if (npol.eq.2) then
         call nga_scale_patch(g_rhs_noJK,blo,bhi,-2.0d0)
         call nga_scale_patch(g_rhs_1e  ,blo,bhi,-2.0d0)
         call nga_scale_patch(g_rhs_eSji,blo,bhi,-2.0d0)
        endif
       enddo ! end-loop-ispin
      endif
c ------- Store g_rhs without Coulomb and Exchange ---- END

      call add_fock1(g_rhs, ! out: accumulated rhs expression
     &               g_fock,! in 
     &               nmo,npol,nocc,nvirt) 
      if (switch_gshift_analysis) then
       call add_fock1(g_rhs_Coul, ! out: accumulated rhs expression
     &                g_fock_Coul,! in 
     &                nmo,npol,nocc,nvirt) 
        if (.not.ga_destroy(g_fock_Coul)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_fock_Coul',
     &           0,GA_ERR)
       call add_fock1(g_rhs_Exch, ! out: accumulated rhs expression
     &                g_fock_Exch,! in 
     &                nmo,npol,nocc,nvirt) 
        if (.not.ga_destroy(g_fock_Exch)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_fock_Exch',
     &           0,GA_ERR)
      endif
      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0) then
        write(*,*) 'Aft add_fock1()'
        write(*,*) '---- g_rhs-aft-add_fock1 -------- START'
       endif
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-aft-add_fock1 -------- END'
      endif

      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_gshift: ga_destroy failed g_fock',0,GA_ERR)
c ----> SKIPPING CPHF ------------------ START
      if (switch_skip_cphf) then
       if (switch_gshift_analysis) then
        if (ga_nodeid().eq.0)
     &     write(*,*) 'WARNING : Using skip_cphf_JK ...'
           call skip_cphf_JK(
     &                    g_rhs,          ! IN/OUT
     &                    g_rhs_Coul,
     &                    g_rhs_Exch,
     &                    g_rhs_noJK,
     &                    g_rhs_1e,
     &                    g_rhs_eSji,
     &                    dbl_mb(k_eval), ! IN: energy eigenvalues
     &                    nocc,
     &                    nvirt,
     &                    nbf,  ! FA-04-28-12: replacing nmo by nbf
     &                    npol)
       else
         if (ga_nodeid().eq.0)
     &     write(*,*) 'WARNING : Using skip_cphf ...'
           call skip_cphf(g_rhs,          ! IN/OUT
     &                    dbl_mb(k_eval), ! IN: energy eigenvalues
     &                    nocc,
     &                    nvirt,
     &                    nbf,  ! FA-04-28-12: replacing nmo by nbf
     &                    npol)
       endif
      endif ! END-skip-cphf
c ----> SKIPPING CPHF ------------------ END
c -------free allocated memory -------------- START
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_occ',0,MA_ERR)
c -------free allocated memory -------------- END
c     OUTPUT: g_rhs

      if (debug_gshift.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- END'
      endif
      if(.not.rtdb_get(rtdb,'zora:skip_cphf_ev_gshift',
     &                 mt_log,1,skip_cphf_ev_gshift))        
     &  skip_cphf_ev_gshift = .false.       
      if (.not.(switch_skip_cphf) .and.
     &    .not.(skip_cphf_ev_gshift)) then ! START-check-skip-cphf
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'COMPUTE cphf g-shift data ...'
c     Write ga_rhs to disk 
       call cphf_fname('cphf_rhs',cphf_rhs)
       call cphf_fname('cphf_sol',cphf_sol)
       if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_gshift: could not write cphf_rhs',0, DISK_ERR)
       call schwarz_tidy()
       call int_terminate()
c     Call the CPHF routine    
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
       if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('hnd_gshift: failed to write skew ', 0, RTDB_ERR)
       if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_gshift: failure in cphf ',0, RTDB_ERR)
       if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('hnd_gshift: rtdb_delete failed ', 0, RTDB_ERR)
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
       call ga_zero(g_rhs)
       if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_gshift: could not read cphf_rhs',0, DISK_ERR) 
c ----- write CPHF data to file ----------------- START
c       Note.- lbl_cphfhyp defined in zora.fh
       call util_file_name(lbl_cphfgshift,
     &                     .false.,.false.,zorafilename)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-gshift --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-gshift --------- END'
       call dft_zoraCPHF_write(
     &           zorafilename, ! in: filename
     &           npol,         ! in: nr polarization
     &           nocc,         ! in: nr occupied MOs
     &           nvirt,        ! in: nr virtual  MOs
     &           nbf,          ! in: nr basis functions
     &           vectors,      ! in: MOs
     &           g_rhs0,       ! out: (ntot,3)       GA matrix
     &           g_rhs)        ! in: (nocc*nvirt,3) GA matrix

c ----- write CPHF data to file ----------------- END
      else !--------------------------- else-check-skip-cphf
       if (ga_nodeid().eq.0)
     &  write(*,*) 'WARNING: SKIP cphf ...'
       if (skip_cphf_ev_gshift) then
        call ga_zero(g_rhs)
        call ga_zero(g_rhs0)
        do i=1,npol
         call ga_zero(vectors(i))
        enddo
        if (ga_nodeid().eq.0)
     &   write(*,*) 'READ cphf g-shift data from file ...'
        call util_file_name(lbl_cphfgshift,.false.,.false.,zorafilename)
        call dft_zoraCPHF_read(
     &           zorafilename, ! in: filename
     &           npol,         ! in: nr polarization
     &           nocc,         ! in: nr occupied MOs
     &           nvirt,        ! in: nr virtual  MOs
     &           nbf,          !  in: nr basis functions
     &           vectors,      ! out: MOs
     &           g_rhs0,       ! out: (ntot,3)      GA matrix
     &           g_rhs)        ! in: (nocc*nvirt,3) GA matrix
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-gshift-read --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-gshift-read --------- END'
       endif
      endif !--------------------------- END-check-skip-cphf
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh para',l_para,k_para)) 
     &    call errquit('hnd_gshift: ma_push_get failed k_para',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_gshift: ma_push_get failed k_dia',0,MA_ERR)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++ Obtain: Alpha or Beta contributions ---- START
      if (icalczora_AB .eq. 0 ) then
c      if (ga_nodeid().eq.0)
c    &   write(*,*) 'BEF. get_gshift_tensor_AorB ...'
      call get_gshift_tensor_AorB(
     &                  rtdb,geom,basis,
     &                  g_rhs,
     &                  g_rhs0,
     &                  g_rhs_Coul,
     &                  g_rhs_Exch,
     &                  g_rhs_noJK,
     &                  g_rhs_1e,
     &                  g_rhs_eSji,
     &                  ga_h01_epr,
     &                  ga_para1_epr_AB,
     &                  ga_dia_epr_AB,
     &                  vectors,
     &                  g_CiFull,
     &                  coeffpol,
     &                  nbf,nmo,
     &                  npol,nocc,nvirt,
     &                  sh_atom, ! in : nr. of atoms
     &                  do_zora,do_NonRel,
     &                  not_zora_scale,
     &                  switch_gshift_analysis)
       if (.not.ga_destroy(ga_dia_epr_AB)) call 
     &    errquit('hnd_gshift_zora: ga_destroy failed g_d1_epr_AB',
     &             0,GA_ERR)
       if (.not.ga_destroy(ga_para1_epr_AB)) call 
     &    errquit('hnd_gshift_zora: ga_destroy failed g_para1_epr_AB',
     &             0,GA_ERR)
      endif
c +++++++++ Obtain: Alpha or Beta contributions ---- END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
c      if (debug_gshift.eq.1) then
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- START'
c      call ga_print(g_rhs)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- END'  
c      endif
      type_NMR=3 ! =1,2,3=shieldings,hyperfine,gshift
       if (switch_gshift_analysis) then !-- START-if-switch_gshift_analysis
        call get_P10_JK(
     &          g_d1_oo,        ! out: Perturbed density matrix occ-occ  contrib
     &          g_d1_ov,        ! out: Perturbed density matrix occ-virt contrib
     &          type_NMR,       ! in : =1,2,3=shieldings,hyperfine,gshift
     &          g_d1_ov_Coul,   ! out
     &          g_d1_ov_Exch,   ! out
     &          g_d1_ov_noJK,   ! out
     &          g_d1_ov_1e,     ! out
     &          g_d1_ov_eSji,   ! out
     &          g_rhs,          ! in: accumulated rhs expression
     &          g_rhs_Coul,     ! in
     &          g_rhs_Exch,     ! in
     &          g_rhs_noJK,     ! in
     &          g_rhs_1e,       ! in
     &          g_rhs_eSji,     ! in
     &          g_rhs0,         ! in: from get_prelim_fock()
     &          vectors,g_CiFull, 
     &          nbf,nmo,npol,nocc,nvirt,
     &          do_zora,do_NonRel,not_zora_scale,
     &          lbl_nlmogshift, ! in: for g-shift NLMO analysis
     &          lbl_nlmoshield, ! in: for shield  NLMO analysis
     &          rtdb) 
        if (.not.ga_destroy(g_rhs_Coul)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_rhs_Coul',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_rhs_Exch)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_rhs_Exch',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_rhs_noJK)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_rhs_noJK',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_rhs_1e)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_rhs_1e',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_rhs_eSji)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_rhs_eSji',
     &           0,GA_ERR)
       else
c      call get_P10( g_d1, ! out: Perturbed density matrix
c     &             g_rhs, ! in: accumulated rhs expression
c     &            g_rhs0, ! in: from get_prelim_fock()
c     &            vectors,g_CiFull, 
c     &            nbf,nmo,npol,nocc,nvirt,
c     &            do_zora,do_NonRel,not_zora_scale,rtdb) 
c -------- checking input to get_P10_1 --------- START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------g_rhs0 ------- START'
c        call ga_print(g_rhs0)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------g_rhs0 ------- END'
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------vectors ------- START'
c        call ga_print(vectors(1))
c        call ga_print(vectors(2))
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------vectors ------- END'
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------g_CiFull ------- START'
c        call ga_print(g_CiFull)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '------g_CiFull ------- END'
c -------- checking input to get_P10_1 --------- END
        call get_P10_1( 
     &         g_d1_oo,        ! out: Perturbed density matrix occ-occ  contrib
     &         g_d1_ov,        ! out: Perturbed density matrix occ-virt contrib
     &         type_NMR,       ! in : =1,2,3=shieldings,hyperfine,gshift
     &         g_rhs,          ! in: accumulated rhs expression
     &         g_rhs0,         ! in: from get_prelim_fock()
     &         vectors,g_CiFull, 
     &         nbf,nmo,npol,nocc,nvirt,
     &         do_zora,do_NonRel,not_zora_scale,
     &         lbl_nlmogshift, ! in: for g-shift NLMO analysis
     &         lbl_nlmoshield, ! in: for shield  NLMO analysis
     &         rtdb) 
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '------- g_d1_ov --------- START'
c         call ga_print(g_d1_ov)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '------- g_d1_ov --------- END'

c      if (ga_nodeid().eq.0)
c     &  write(*,*) '--BEF: get_P10_JK ---g_rhs ------ START'
c       call ga_print(g_rhs)
c      if (ga_nodeid().eq.0)
c     &  write(*,*) '--BEF: get_P10_JK ---g_rhs ------ END'
c      if (ga_nodeid().eq.0)
c     &  write(*,*) '--BEF: get_P10_JK ---g_rhs_noJK ------ START'
c       call ga_print(g_rhs_noJK)
c      if (ga_nodeid().eq.0)
c     &  write(*,*) '--BEF: get_P10_JK ---g_rhs_noJK ------ END'
       endif !--------------- END-if-switch_gshift_analysis
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the shielding tensor
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
c     s(para)xy = Sum(n,l) D1x(n,l) * H01y(n,l)
      do ixy = 1, 9*sh_atom
       dbl_mb(k_para+ixy-1) = 0.0d0  ! initialize the paramagnetic part
      enddo
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++ Initialize dbl_mb(k_para) with ga_para1 ++++ START
c ---- STORE: g_para1_epr --> dbl_mb(k_para)
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=1
       ld(1)=3
       ld(2)=3
      call nga_get(ga_para1_epr,alo,ahi,dbl_mb(k_para),ld)
c +++++ Initialize dbl_mb(k_para) with g_para1 ++++ END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      ixy = 0
      blo(3) = 0
      bhi(3) = 0
      do iy = 1, 3
       blo(3) = blo(3) + 1
       bhi(3) = bhi(3) + 1
       do ix = 1, 3
        alo(3) = ix
        ahi(3) = ix
        ixy = ixy + 1
        if (switch_gshift_analysis) then !-- START-if-switch_gshift_analysis
         val_oo = nga_ddot_patch(g_d1_oo,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov = nga_ddot_patch(g_d1_ov,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_Coul = nga_ddot_patch(g_d1_ov_Coul,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_Exch = nga_ddot_patch(g_d1_ov_Exch,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_noJK = nga_ddot_patch(g_d1_ov_noJK,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)  
         val_ov_1e = nga_ddot_patch(g_d1_ov_1e,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)   
         val_ov_eSji = nga_ddot_patch(g_d1_ov_eSji,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)             
        else
         val_oo = nga_ddot_patch(g_d1_oo,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov = nga_ddot_patch(g_d1_ov,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
        endif !-- END-if-switch_gshift_analysis
        cbuf=3*(ix-1)+iy-1 ! transpose
c ----- store in par_arr ------- START
        par_arr(1,iy,ix)=dbl_mb(k_para+cbuf)
        par_arr(2,iy,ix)=val_oo*coeffpol*ppt*(-0.25d0)
        par_arr(3,iy,ix)=val_ov*coeffpol*ppt*(-0.25d0)
        par_arr(4,iy,ix)=dbl_mb(k_para+cbuf)+
     &                   val_oo*coeffpol*ppt*(-0.25d0)+
     &                   val_ov*coeffpol*ppt*(-0.25d0)
        if (switch_gshift_analysis) then !-- START-if-switch_gshift_analysis
         par_arr(5,iy,ix)=val_ov_Coul*coeffpol*ppt*(-0.25d0)
         par_arr(6,iy,ix)=val_ov_Exch*coeffpol*ppt*(-0.25d0)
         par_arr(7,iy,ix)=val_ov_noJK*coeffpol*ppt*(-0.25d0)
         par_arr(8,iy,ix)=val_ov_1e*coeffpol*ppt*(-0.25d0)
         par_arr(9,iy,ix)=val_ov_eSji*coeffpol*ppt*(-0.25d0)
        endif
c ----- store in par_arr ------- END      
        dbl_mb(k_para+cbuf)=dbl_mb(k_para+cbuf)+
     &                      val_oo*coeffpol*ppt*(-0.25d0)+
     &                      val_ov*coeffpol*ppt*(-0.25d0)
       enddo
      enddo
c ------------- symmetrize total para ------------ START
      do iy = 1, 3
       do ix = iy+1, 3
        cbuf1=k_para+3*(ix-1)+iy-1 ! transpose
        cbuf2=k_para+3*(iy-1)+ix-1 ! transpose
        val=(dbl_mb(cbuf1)+dbl_mb(cbuf2))/2.0d0
        dbl_mb(cbuf1)=val
        dbl_mb(cbuf2)=val        
       enddo
      enddo
c ------------- symmetrize total para ------------ END
c
c     s(dia)xy = Sum(n,l) D0(n,l) * H11(dia)xy(n,l) 
c     Ordering of H11 blocks is Bxx,Bxy,Bxz,Byx,etc 
      do ixy = 1,9*sh_atom
         dbl_mb(k_dia+ixy-1)  = 0.0d0  ! initialize the diamagnetic part
      enddo
c ---- STORE: ga_dia_epr --> dbl_mb(k_dia)
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=sh_atom
      ld(1)=3
      ld(2)=3 
      call nga_get(ga_dia_epr,alo,ahi,dbl_mb(k_dia),ld)
c -------- symmetrize total dia ------------ START
      do ix = 1,3
       do iy = ix+1,3
        cbuf1=k_dia+3*(ix-1)+iy-1 ! transpose
        cbuf2=k_dia+3*(iy-1)+ix-1 ! transpose
        val=(dbl_mb(cbuf1)+dbl_mb(cbuf2))/2.0d0
        dbl_mb(cbuf1)=val
        dbl_mb(cbuf2)=val          
       enddo
      enddo
c -------- symmetrize total dia ------------ END
c +++++++++ print-total-pardia-transferred +++ START
       if (ga_nodeid().eq.0) then
        ic=1
        do ix = 1, 3
         do iy = 1, 3
          if (switch_gshift_analysis) then !-- START-if-switch_gshift_analysis
           write(*,179) ix,iy,dbl_mb(k_dia+ic-1),
     &                 par_arr(1,ix,iy),par_arr(2,ix,iy),
     &                 par_arr(3,ix,iy),
     &                 par_arr(5,ix,iy),par_arr(6,ix,iy),
     &                 par_arr(7,ix,iy),
     &                 par_arr(8,ix,iy),par_arr(9,ix,iy),
     &                 par_arr(4,ix,iy)
 179        format('NW-T:(dia,gauge,OO,OV,',
     &            'OV_Coul,OV_Exch,OV_nJK,',
     &            'OV_1e,OV_eSji,'
     &            'Totpar)(',i1,',',i1,')=(',
     &             f12.6,' ',f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' )')
          else
           write(*,19) ix,iy,dbl_mb(k_dia+ic-1),
     &                 par_arr(1,ix,iy),par_arr(2,ix,iy),
     &                 par_arr(3,ix,iy),par_arr(4,ix,iy)
 19        format('NW-T:(dia,gauge,OO,OV,Totpar)(',i1,',',i1,')=(',
     &             f12.6,' ',f12.6,' ',f12.6,' ',f12.6,' ',
     &             f12.6,' )')
          endif
          ic=ic+1
         enddo ! end-loop-iy
        enddo ! end-loop-ix
       endif ! end-if-ga_nodeid-eq-0
c +++++++++ print-total-dia-transferred +++ END
       if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_gshift_zora: ga_destroy failed g_rhs',
     &            0,GA_ERR)
       if (.not.ga_destroy(g_rhs0)) call 
     &    errquit('hnd_gshift_zora: ga_destroy failed g_rhs0',
     &            0,GA_ERR)
       do ispin=1,ndens
        if (.not.ga_destroy(g_dens(ispin))) call 
     &    errquit('hnd_gshift_zora: ga_destroy failed g_dens',
     &    0,GA_ERR)
       enddo
c
c     Print out tensor information, and write to Ecce file if necessary
      status = rtdb_parallel(.false.)   
      if (ga_nodeid().gt.0) goto 300
      acc_vec=0 ! For NMLO analysis
      call ecce_print_module_entry('nmr')
      do iatom = 1, sh_atom
       ioff = (iatom-1)*9
       if (.not. geom_cent_get(geom, int_mb(k_tmp+iatom-1), tag, 
     &     dbl_mb(k_xyz), dbl_mb(k_zan))) call 
     &     errquit('hnd_gshift: geom_cent_tag failed',0,GEOM_ERR)
       if (.not. geom_tag_to_element(tag, symbol, element, atn)) call
     &     errquit('hnd_gshift: geom_tag_to_element failed',0,GEOM_ERR)
c
c      Print tensor pieces and sum for total shielding tensor

         if (ga_nodeid().eq.0) then
c            write(luout,9700) iatom,symbol
            write(luout,9800) (dbl_mb(k_dia+ioff+ix-1),ix=1,9)
            write(luout,9801) (dbl_mb(k_para+ioff+ix-1),ix=1,9)
         endif
         do ix = 0, 8 
            dbl_mb(k_para+ioff+ix) = dbl_mb(k_dia+ioff+ix) + 
     &                               dbl_mb(k_para+ioff+ix)
         enddo
c
c     Print total shielding tensor
         if (ga_nodeid().eq.0) then
            write(luout,9802) (dbl_mb(k_para+ioff+ix-1),ix=1,9)
c
c     Diagonalize total tensor
c     Order in a: xx xy yy xz yz zz 
            a(1) = dbl_mb(k_para+ioff)     
            a(2) = dbl_mb(k_para+ioff+1)
            a(3) = dbl_mb(k_para+ioff+4)
            a(4) = dbl_mb(k_para+ioff+2)
            a(5) = dbl_mb(k_para+ioff+5)
            a(6) = dbl_mb(k_para+ioff+8)
            ij = 0
            do 241 i = 1, 3
            do 241 j = 1, i
               ij = ij + 1
               axs(i,j) = a(ij)
               axs(j,i) = a(ij)
  241       continue
            call hnd_diag(axs,eig,3,.true.,.true.)
            isotr =(eig(1) + eig(2) + eig(3))/3.0d0
            aniso = eig(1) -(eig(2) + eig(3))/2.0d0
c ++++++++++ get eigenvectors for NMLO analysis +++++ START
            gshiftfile=0 ! not doing NLMO analysis by default
            status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
            if (gshiftfile.eq.1) then
             do i1=1,3
              do j1=1,3
               dbl_mb(k_tvec+acc_vec)=axs(i1,j1)
               acc_vec=acc_vec+1
              enddo 
             enddo 
            endif
c ++++++++++ get eigenvectors for NMLO analysis +++++ END
            write(luout,9987) isotr,aniso
            write(luout,9986) (ix,ix=1,3)
            write(luout,9985) (eig(ix),ix=1,3)
            do iy=1,3
              write(luout,9983) iy,(axs(iy,ix),ix=1,3)
            enddo
            write(luout,'(//)')
c
c     Print Ecce information

            call ecce_print1_char('atom name',symbol,1)
            call ecce_print2('g-shift tensor',MT_DBL,
     &                       dbl_mb(k_para+ioff),3,3,3)
            call ecce_print1('g-shifts isotropic',MT_DBL,isotr,1)
            call ecce_print1('g-shifts anisotropy',MT_DBL,aniso,1)
            call ecce_print1('g-shifts eigenvalues',MT_DBL,eig,3)
            call ecce_print2('g-shifts eigenvectors',MT_DBL,axs,
     &                       3,3,3)
         endif
      enddo
      call ecce_print_module_exit('nmr','ok')
300   call ga_sync()
      status = rtdb_parallel(.true.)   
      gshiftfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
      if (gshiftfile.eq.1) then ! ------- gshiftfile-if++++ START
         if (.not. ga_create(mt_dbl,1,9,
     &       'munu4nbo: g_tvec',0,0,g_tvec)) 
     $       call errquit('munu4nbo: g_tvec', 0,GA_ERR)
        call ga_dgop(msg_efgs_col,dbl_mb(k_tvec),9,'+')
        call ga_put(g_tvec,1,1,1,9,dbl_mb(k_tvec),1)      
        call create_munu4nbo_gshift(
     &                           rtdb,g_tvec,
     &                           basis,npol,nocc,nvirt,nmo) 
         if (.not. ga_destroy(g_tvec)) call errquit( ! destroy GA
     &    'hnd_hyperfine_zora: ga_destroy failed ',0, GA_ERR)   
       if (.not.ma_free_heap(l_tvec)) call
     &     errquit('hnd_hyperfine_zora: ma_free_heap l_tvec',0,MA_ERR)    
      endif ! ------------------------ gshiftfile-if++++ END
c ---- Destroy stored ga arrays ------ START
c      if (.not.ga_destroy(g_d1)) call 
c     &    errquit('hnd_gshift: ga_destroy failed g_d1',0,GA_ERR)
c ----- destroy created GA arrays in get_P10_JK ---- START
       if (.not.ga_destroy(g_d1_oo)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_oo',
     &           0,GA_ERR)
       if (.not.ga_destroy(g_d1_ov)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov',
     &           0,GA_ERR)
       if (switch_gshift_analysis) then
        if (.not.ga_destroy(g_d1_ov_Coul)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov_Coul',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_Exch)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov_Exch',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_noJK)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov_noJK',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_1e)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov_1e',
     &           0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_eSji)) call 
     &  errquit('hnd_gshift_zora: ga_destroy failed g_d1_ov_eSji',
     &           0,GA_ERR)
       endif
c ----- destroy created GA arrays in get_P10_JK ---- END
           if (.not. ga_destroy(ga_dia_epr)) call errquit(
     &    'hnd_gshift_zora: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_para1_epr)) call errquit(
     &    'hnd_gshift_zora: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_epr)) call errquit(
     &    'hnd_gshift_zora: ga_destroy failed ',0, GA_ERR)   
        if (.not. ga_destroy(ga_Fji)) call errquit(
     &    'hnd_gshift_zora: ga_destroy failed ',0, GA_ERR)    
c        if (.not. ga_destroy(g_AtNr1)) call errquit(
c     &    'hnd_gshift_zora: ga_destroy failed ',0, GA_ERR)          
c ---- Destroy stored ga arrays ------ END
c
c     Clean up all remaining memory
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_dia',0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_para',0,MA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('giao_aotomo: ga_destroy failed vectors',0,GA_ERR)
      enddo
      if (.not.ma_pop_stack(l_zan)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_zan',0,MA_ERR)
      if (.not.ma_pop_stack(l_xyz)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_xyz',0,MA_ERR)
      if (.not.ma_pop_stack(l_tmp)) call
     &    errquit('hnd_gshift: ma_pop_stack failed l_tmp',0,MA_ERR)
      call schwarz_tidy()
      call int_terminate()
c
      return
 7000 format(/,10x,'EPR g-shifts cannot be calculated for UHF',
     1      ' or ROHF wave functions: use RHF')
 9700 format(6x,'Atom: ',i4,2x,a2)
 9800 format(8x,'Diamagnetic',/,3(3F12.4,/))
 9801 format(8x,'Paramagnetic',/,3(3F12.4,/))
 9802 format(8x,'Total g shift Tensor',/,3(3F12.4,/))
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9986 format(10x,'Principal Components and Axis System',/,10x,
     1       3(7x,i1,4x))
 9987 format(10x,' isotropic = ',f12.4,/,
     1       10x,'anisotropy = ',f12.4,/)
 9999 format(
     1 /,10x,41(1h-),/,
     2 10x,'g shift Tensors (GIAO, in ppt)',/,
     3 10x,41(1h-),/)
      end
 
      subroutine skip_cphf(g_rhs, ! IN/OUT
     &                     eval,  ! IN: energy eigenvalues
     &                     nocc,  ! IN: nr. occupied MOs
     &                     nvirt, ! IN: nr. virtual  MOs
     &                     nbf,   ! IN: nr. basis functions
     &                     npol)  ! IN: nr. polarizations  
c     Purpose : Avoid usage of cphf() routine in
c               the case of BP functionals.
c               It produces the output ready to
c               be used to comput 1st derivative
c               of density matrix (P^{10})
c               done by routine get_P10_1()
c     Author : Fredy W. Aquino
c     Date   : 02-22-11
c     Update-04-28-12: Replacing nmo by nbf
c -------- Calculation g_rhs/(ej-ei) ----- START
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer g_rhs
      integer disp,disp1,iocc,ivirt,ind
      integer ispin,alo(3),ahi(3) 
      integer nbf,npol,nocc(2),nvirt(2)
      double precision eval(nbf*npol),tosclji,coeff

      if (npol.eq.1) then
        coeff=0.25d0
      else if (npol.eq.2) then
        coeff=0.5d0
      else
       write(*,*) 'Error in skip_cphf: npol not 1 or 2!'
       stop
      endif         
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- BEF scaling g_rhs -------- START'
c       call ga_print(g_rhs)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- BEF scaling g_rhs -------- END'

       if (ga_nodeid().eq.0) then
        write(*,97) npol,nocc(1),nocc(2),
     &              nvirt(1),nvirt(2),nbf    
 97     format('(npol,nocc,nvirt,nbf)=[',i5,',(',
     &          i8,',',i8,'),(',i8,',',i8,'),',
     &          i8,']')               
       endif
       do ispin=1,npol 
         disp=nbf*(ispin-1) 
         disp1=nocc(1)*nvirt(1)*(ispin-1)
         do ivirt = 1, nvirt(ispin)
          do iocc = 1, nocc(ispin)
           tosclji=coeff/(eval(disp+nocc(ispin)+ivirt)-
     &                    eval(disp+iocc)) 
c           if (ga_nodeid().eq.0) then
c            write(*,16) iocc,ivirt,tosclji
c 16          format('1/Diff-E(',i5,',',i5,')=',f15.8) 
c           endif
           ind=disp1+nvirt(ispin)*(iocc-1)+ivirt 
           alo(1) = ind
           ahi(1) = ind
           alo(2) = 1
           ahi(2) = 3
           alo(3) = 0
           ahi(3) = 0 
           call nga_scale_patch(g_rhs,alo,ahi,tosclji)      
          enddo ! end-loop-iocc    
         enddo ! end-loop-ivirt
       enddo  ! end-loop-ispin
c -------- Calculation g_rhs/(ej-ei) ----- END
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- AFTER scaling g_rhs -------- START'
c       call ga_print(g_rhs)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- AFTER scaling g_rhs -------- END'
      return
      end

c ------ For testing Coulomb,Exchange contrib ----- START
      subroutine skip_cphf_JK(
     &                     g_rhs, ! IN/OUT
     &                     g_rhs_Coul,
     &                     g_rhs_Exch,
     &                     g_rhs_noJK,
     &                     g_rhs_1e,
     &                     g_rhs_eSji,
     &                     eval,  ! IN: energy eigenvalues
     &                     nocc,  ! IN: nr. occupied MOs
     &                     nvirt, ! IN: nr. virtual  MOs
     &                     nbf,   ! IN: nr. basis functions
     &                     npol)  ! IN: nr. polarizations  
c     Purpose : Avoid usage of cphf() routine in
c               the case of BP functionals.
c               It produces the output ready to
c               be used to comput 1st derivative
c               of density matrix (P^{10})
c               done by routine get_P10_1()
c     Author : Fredy W. Aquino
c     Date   : 02-22-11
c -------- Calculation g_rhs/(ej-ei) ----- START
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer g_rhs,g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_1e,g_rhs_eSji
      integer disp,disp1,iocc,ivirt,ind
      integer ispin,alo(3),ahi(3) 
      integer nbf,npol,nocc(2),nvirt(2)
      double precision eval(nbf*npol),tosclji,coeff
      logical debug_skip

      debug_skip= .false. ! =1 for debugging 

      if (npol.eq.1) then
        coeff=0.25d0
      else if (npol.eq.2) then
        coeff=0.5d0
      else
       write(*,*) 'Error in skip_cphf: npol not 1 or 2!'
       stop
      endif         
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- BEF scaling g_rhs -------- START'
c       call ga_print(g_rhs)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- BEF scaling g_rhs -------- END'

       if (ga_nodeid().eq.0) then
        write(*,97) npol,nocc(1),nocc(2),
     &              nvirt(1),nvirt(2),nbf    
 97     format('(npol,nocc,nvirt,nbf)=[',i5,',(',
     &          i8,',',i8,'),(',i8,',',i8,'),',
     &          i8,']')               
       endif
       do ispin=1,npol 
         disp=nbf*(ispin-1)
         disp1=nocc(1)*nvirt(1)*(ispin-1)
         do ivirt = 1, nvirt(ispin)
          do iocc = 1, nocc(ispin)
           tosclji=coeff/(eval(disp+nocc(ispin)+ivirt)-
     &                    eval(disp+iocc)) 
           if (debug_skip) then
            if (ga_nodeid().eq.0) then
              write(*,16) iocc,ivirt,tosclji
 16           format('1/Diff-E(',i5,',',i5,')=',f15.8) 
            endif
           endif
           ind=disp1+nvirt(ispin)*(iocc-1)+ivirt 
           alo(1) = ind
           ahi(1) = ind
           alo(2) = 1
           ahi(2) = 3
           alo(3) = 0
           ahi(3) = 0 
           call nga_scale_patch(g_rhs,alo,ahi,tosclji)     
           call nga_scale_patch(g_rhs_Coul,alo,ahi,tosclji)       
           call nga_scale_patch(g_rhs_Exch,alo,ahi,tosclji) 
           call nga_scale_patch(g_rhs_noJK,alo,ahi,tosclji) 
           call nga_scale_patch(g_rhs_1e,alo,ahi,tosclji) 
           call nga_scale_patch(g_rhs_eSji,alo,ahi,tosclji)
          enddo ! end-loop-iocc    
         enddo ! end-loop-ivirt
       enddo  ! end-loop-ispin
c -------- Calculation g_rhs/(ej-ei) ----- END
       goto 111
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs -------- END'
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs_Coul ---- START'
       call ga_print(g_rhs_Coul)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs_Coul ---- END'
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs_Exch ---- START'
       call ga_print(g_rhs_Exch)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- AFTER scaling g_rhs_Exch ---- END'
 111   continue
      return
      end
c ------ For testing Coulomb,Exchange contrib ----- END

      subroutine get_prelim_fock( g_d2, ! out: 
     &                           g_rhs, ! out: rhs expression
     &                          g_rhs0, ! out: to be used in get_d1()
     &                         vectors, !  in: MO  coeffs
     &                            eval, !  in: energy vals
     &                             pos, !  in: Nuclear positions (x,y,z)
     &                          natoms, !  in: nr. selected nuclei (atoms)
     &                           basis, !  in: basis handle
     &                             nbf, !  in: nr. basis functions
     &                             nmo, !  in: nr. MOs (occ+virt)
     &                            npol, !  in: nr. of polarizations
     &                            nocc, !  in: nr. occ     MOs
     &                           nvirt) !  in: nr. virtual MOs
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol  ! nr. of polarizations
      integer g_rhs, ! OUT-1
     &        g_d1,  ! tmp
     &        g_d2,  ! OUT-2
     *        g_rhs0 ! OUT-3
      integer g_u,g_s10,g_s10_1
      integer ntot,nbf,nmo,ispin,iocc,
     &        nocc(npol),nvirt(npol)
      integer shift,disp,xyz
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        dlo(3), dhi(3)
      integer natoms
      double precision pos(3*natoms)
      double precision eval(nmo*npol),toscl
      integer basis ! basis handle
      integer vectors(npol)
      logical oskel
      external int_giao_1ega,giao_aotomo
      integer debug_prelim
      integer ndens

      debug_prelim=0 ! =1 for debugging

      oskel = .false.
c ------------ creating ga arrays ---------- START
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nocc(ispin)
      enddo
      if(.not.ga_create(MT_DBL,ntot,3,'rhs0',-1,-1,g_rhs0))
     &   call errquit('get_prelim_fock: ga_create failed g_rhs0',
     &                 0,GA_ERR)
      call ga_zero(g_rhs0)
      clo(1) = 3*npol
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'g_d1 matrix',
     &                    chi,g_d1)) 
     &  call errquit('gprelim_fock: nga_create failed g_d1',
     &                0,GA_ERR)
      call ga_zero(g_d1)

c     Get S10 in GA and transform to MO set (virt,occ)
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10_1)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10_1',
     &               0,GA_ERR)
       call ga_zero(g_s10_1)
       ahi(3) = 3*npol
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10',
     &               0,GA_ERR)
       call ga_zero(g_s10)

c ------------ creating ga arrays ---------- END
       call int_giao_1ega(basis,basis,g_s10_1,'s10',
     &                    pos(1),natoms,oskel)    

       blo(1) = 1
       bhi(1) = nbf ! nmo-fixing4lineardependency
       blo(2) = 1
       bhi(2) = nbf ! nmo-fixing4lineardependency
       blo(3) = 1
       bhi(3) = 3
      do ispin=1,npol  
       disp=3*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf ! nmo-fixing4lineardependency
       alo(2) = 1
       ahi(2) = nbf ! nmo-fixing4lineardependency
       alo(3) = disp+1
       ahi(3) = disp+3   
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin

      if (debug_prelim.eq.1) then
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10 ----START'
        call ga_print(g_s10)
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10 ----END'
      endif

      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,3,nbf)
      if (debug_prelim.eq.1) then
        do ispin=1,npol
          if (ga_nodeid().eq.0) 
     &     write(*,*) '-------- MO vect(',ispin,')----START'
         call ga_print(vectors(ispin))
          if (ga_nodeid().eq.0) 
     &     write(*,*) '-------- MO vect(',ispin,')----END'
        enddo
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ----gprelim: ---- MO:g_s10 ---- START'
       call ga_print(g_s10)
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ----gprelim: ---- MO:g_s10 ---- END'
      endif

      do ispin=1,npol ! ++++++  START-loop-ispin
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) =  3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &                     alo,g_u)) 
     &   call errquit('gprelim_fock: nga_create failed g_u',
     &                0,GA_ERR)
       call ga_zero(g_u)
       disp=3*(ispin-1)
c     COPY : g_s10 -> g_s10_1 ---------- START
       blo(1) = 1
       bhi(1) = nmo
       blo(2) = 1
       bhi(2) = nmo
       blo(3) = 1
       bhi(3) = 3
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nmo
       alo(3) = disp+1
       ahi(3) = disp+3   
       call ga_zero(g_s10_1) ! FA-04-28-12-added
       call nga_copy_patch('n',g_s10,alo,ahi,
     &                       g_s10_1,blo,bhi) 
c     COPY : g_s10 -> g_s10_1 ---------- END    
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi    
c ----definitions for g_rhs -------- START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c --- definitions for g_rhs -------- END   
c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block g_s10 with - (minus) eigenvalues 
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(3) = 1
       ahi(3) = 3
       disp=nbf*(ispin-1)  ! FA-04-28-12-energy-fix

       if (debug_prelim.eq.1) then
         write(*,*) 'get_prelim_fock: nmo=',nmo,
     &              ' nbf=',nbf
         if (ga_nodeid().eq.0) 
     &    write(*,*) ' ---gprelim: ---- BEF-scl:g_s10_1 ----START'
         call ga_print(g_s10_1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) ' ---gprelim: ---- BEF-scl:g_s10_1 ----END'    
       endif

        do iocc = 1, nocc(ispin)
         alo(2) = iocc
         ahi(2) = iocc
         toscl=-eval(disp+iocc) 

          if (debug_prelim.eq.1) then
           if (ga_nodeid().eq.0) then
           write(*,1) iocc,toscl
 1         format('E2scl(',i5,')=',f15.8) 
           endif
          endif

         call nga_scale_patch(g_s10_1,alo,ahi,toscl)          
        enddo ! end-loop-iocc
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- AFT-scl:g_s10_1 ----START'
        call ga_print(g_s10_1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- AFT-scl:g_s10_1 ----END'    
      endif

c     Copy to g_rhs 
c     alo(1) and ahi(1) the same as before
       alo(2) = 1
       ahi(2) = nocc(ispin)

      if (debug_prelim.eq.1.and.ga_nodeid().eq.0) then
         write(*,16) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3)
 16      format('alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
      endif

       call nga_copy_patch('n',g_s10_1,alo,ahi,
     &                           g_rhs,blo,bhi)
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- g_rhs(',ispin,') ----START'
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- g_rhs(',ispin,') ----END'    
      endif
c
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
       alo(1) = 1
       ahi(1) = nocc(ispin)
       call nga_scale_patch(g_s10_1,alo,ahi,-0.5d0)
       call nga_copy_patch('n',g_s10_1,alo,ahi,
     &                             g_u,alo,ahi)
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u(',ispin,') ----START'
        call ga_print(g_u)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u(',ispin,') ----END'    
      endif
c
c     We also need the occupied-occupied contribution of g_u contributing
c     to the first order density matrix. As this block does not change 
c     during the CPHF we can calculate it once and subtract it from the
c     RHS. We will reuse g_s10 as scratch space.
        call ga_zero(g_s10_1)
        alo(1) = 1
        alo(2) = 1
        blo(1) = 1
        blo(2) = 1
        bhi(1) = nbf
        clo(2) = 1
        clo(3) = 1
        chi(2) = nbf
        chi(3) = nbf
        dlo(1) = 1
        dlo(2) = 1
        dhi(1) = nbf
        dhi(2) = nocc(ispin)
c     Create "perturbed density matrix" for closed-closed g_u block
       do xyz = 1,3 ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
        bhi(2) = nmo 
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                    vectors(ispin),blo,bhi,  
     &                               g_u,alo,ahi,
     &                           g_s10_1,dlo,dhi)  
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u-1(',ispin,',',xyz,') ----START'
        call ga_print(g_s10_1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u-1(',ispin,',',xyz,') ----END'    
      endif
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
        bhi(2) = nocc(ispin)
c     Minus sign as we subtract it from the RHS as we do not include 
c     it in the LHS
        disp=3*(ispin-1)
        clo(1) = disp+xyz
        chi(1) = disp+xyz
        call nga_matmul_patch('n','t',-1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,
     &                               g_s10_1,alo,ahi,
     &                                  g_d1,clo,chi)  
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_d1(',ispin,',',xyz,') ----START'
        call ga_print(g_d1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_d1(',ispin,',',xyz,') ----END'    
      endif
       enddo ! end-loop-xyz
c ------------ back-up g_u --> g_rhs0 ---- START
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       call nga_copy_patch('n',g_u   ,alo,ahi,
     &                         g_rhs0,blo,bhi) ! copy to g_rhs0
c ------------ back-up g_u --> g_rhs0 ---- END
c ----- Remove scratch ga arrays ------
       if (.not.ga_destroy(g_u)) call 
     &     errquit('gprelim_fock: ga_destroy failed g_u',0,GA_ERR)
      enddo ! end-loop-ispin

c -------- Creating g_d2 ---------------START
      ndens=3*npol
      clo(1) = ndens*2
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'g_d2 matrix',
     &                    chi,g_d2)) 
     &  call errquit('gprelim_fock: nga_create failed g_d2',
     &                0,GA_ERR)
       call ga_zero(g_d2)
       blo(1) = 1
       bhi(1) = ndens   
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 1
       bhi(3) = nbf
      do ispin=1,npol 
       disp=ndens*(ispin-1) 
       alo(1) = disp+1
       ahi(1) = disp+ndens   
       alo(2) = 1
       ahi(2) = nbf
       alo(3) = 1
       ahi(3) = nbf
       call nga_copy_patch('n',g_d1,blo,bhi,
     &                         g_d2,alo,ahi) 
      enddo ! end-loop-ispin

 145  continue

      if (debug_prelim.eq.1) then
       if (ga_nodeid().eq.0) then
        write(*,*) '----------- g_d2------------ START'
       endif
       call ga_print(g_d1)
       call ga_print(g_d2)
       if (ga_nodeid().eq.0) then
        write(*,*) '----------- g_d2------------ END'
       endif
      endif
c -------- Creating g_dens1 ---------------END
       if (.not.ga_destroy(g_s10)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
       if (.not.ga_destroy(g_s10_1)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10_1',0,GA_ERR)
       if (.not.ga_destroy(g_d1)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_d1',0,GA_ERR)
      return
      end

      subroutine add_fock(g_rhs, ! out: accumulated rhs expression
     &                   g_fock, !  in: Fock-term
     &                  vectors, !  in: MO  coeffs
     &                      nbf, !  in: nr. basis functions
     &                      nmo, !  in: nr. MOs (occ+virt)
     &                     npol, !  in: nr. of polarizations
     &                     nocc, !  in: nr. occ     MOs
     &                    nvirt) !  in: nr. virtual MOs
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol  ! nr. of polarizations
      integer g_rhs
      integer g_fock
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,disp,xyz
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3)
      integer vectors(npol)
      integer g_s10_1,g_fck ! scratch ga arrays
c     Transform to the occ-virt MO basis and add to RHS
      clo(1) = 3
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',
     &                    chi,g_fck)) 
     &    call errquit('add_fock: nga_create failed g_fck',
     &                 0,GA_ERR)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10_1)) 
     &    call errquit('add_fock: nga_create failed g_s10_1',
     &                  0,GA_ERR)
      do ispin=1,npol
c     Transfer: g_fock --> g_fck
       shift=3*(ispin-1)
       alo(1) = shift+1
       ahi(1) = shift+3   
       alo(2) = 1
       ahi(2) = nbf ! nmo-fixing4lineardependency
       alo(3) = 1
       ahi(3) = nbf ! nmo-fixing4lineardependency
       blo(1) = 1
       bhi(1) = 3
       blo(2) = 1
       bhi(2) = nbf ! nmo-fixing4lineardependency
       blo(3) = 1
       bhi(3) = nbf ! nmo-fixing4lineardependency
        call nga_copy_patch('n',g_fock,alo,ahi,
     &                           g_fck,blo,bhi) 
c       if (ga_nodeid().eq.0)
c     &  write(*,17)  ispin
c 17     format('------g_fck(',i3,')--------- START')
c        call ga_print(g_fck)
c       if (ga_nodeid().eq.0)
c     &  write(*,18)  ispin
c 18     format('------g_fck(',i3,')--------- END')

       call ga_zero(g_s10_1)
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nocc(ispin)
       clo(2) = 1
       chi(2) = nbf
       clo(3) = 1
       chi(3) = nbf
       do xyz = 1,3 ! =x,y,z
         alo(3) = xyz
         ahi(3) = xyz
         clo(1) = xyz
         chi(1) = xyz
        call nga_matmul_patch('n','n',2.0d0,0.0d0,
     $                            g_fck,clo,chi,
     $                   vectors(ispin),alo,ahi,
     $                          g_s10_1,alo,ahi)
       enddo ! end-loop-xyz
       call ga_zero(g_fck)
       clo(2) = nocc(ispin)+1
       clo(3) = 1
       chi(2) = nmo
       chi(3) = nocc(ispin)
       do xyz = 1,3 ! = x,y,z
         blo(1) = nocc(ispin)+1
         bhi(1) = nmo
         blo(2) = 1
         bhi(2) = nbf
         alo(3) = xyz
         ahi(3) = xyz
         clo(1) = xyz
         chi(1) = xyz
         call nga_matmul_patch('t','n',1.0d0,0.0d0,
     $                     vectors(ispin),blo,bhi,
     $                            g_s10_1,alo,ahi,
     $                              g_fck,clo,chi)
c ------ definitions for g_rhs -------- START
         disp=nocc(1)*nvirt(1)*(ispin-1)
         blo(1) = disp+1
         bhi(1) = disp+nocc(ispin)*nvirt(ispin)
         blo(2) = xyz
         bhi(2) = xyz
         blo(3) = 0
         bhi(3) = 0
c ------ definitions for g_rhs -------- END
         call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                      1.0d0,g_fck,clo,chi,
     &                            g_rhs,blo,bhi)
       enddo ! end-loop-xyz
      enddo ! end-loop-ispin
      if (.not.ga_destroy(g_s10_1)) call 
     &    errquit('add_fock: ga_destroy failed g_s10_1',0,GA_ERR)
      if (.not.ga_destroy(g_fck)) call 
     &    errquit('add_fock: ga_destroy failed g_fck',0,GA_ERR)
      return
      end

      subroutine add_fock1(g_rhs, ! out: accumulated rhs expression
     &                    g_fock, !  in: Fock-term
     &                       nmo, !  in: nr. MOs (occ+virt)
     &                      npol, !  in: nr. of polarizations
     &                      nocc, !  in: nr. occ     MOs
     &                     nvirt) !  in: nr. virtual MOs
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol  ! nr. of polarizations
      integer g_rhs
      integer g_fock 
      integer nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,disp
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3)

      do ispin=1,npol
c ------ definitions for g_rhs -------- START
       disp = nocc(1)*nvirt(1)*(ispin-1)
       shift=3*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c ------ definitions for g_rhs -------- END
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+3
       call nga_add_patch(1.0d0, g_rhs,blo,bhi,
     &                    1.0d0,g_fock,alo,ahi,
     &                           g_rhs,blo,bhi)

       if      (npol.eq.1) then
        call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)
       else if (npol.eq.2) then
        call nga_scale_patch(g_rhs,blo,bhi,-2.0d0)
       endif
      enddo ! ++++++  END-loop-ispin

      return
      end

      subroutine add_H10( g_rhs, ! out: accumulated rhs expression
     &                   ga_Fji, !  in: Fock 1st-deriv without V (pot.) contrib.
     &                  vectors, !  in: MO  coeffs
     &                      pos, !  in: Nuclear positions (x,y,z)
     &                   natoms, !  in: nr. selected nuclei (atoms)
     &                    basis, !  in: basis handle
     &                      nbf, !  in: nr. basis functions
     &                      nmo, !  in: nr. MOs (occ+virt)
     &                     npol, !  in: nr. of polarizations
     &                     nocc, !  in: nr. occ     MOs
     &                    nvirt, !  in: nr. virtual MOs
     &                  do_zora, !  in: .false. if doing nonrel calc
     &                     rtdb) !  in: for COSMO
                                 !      therefore: HF or DFT (restricted/unrestricted)
                                 !      Then evaluate analytical GIAO AOs. instead of numerical
c =============  Computing  ==============
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c ========================================
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer rtdb
      integer npol    ! nr. of polarizations
      logical do_zora ! logical to check if NOT doing zora calc
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,disp
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3)
      integer vectors(npol)
      integer g_rhs,ga_Fji
      integer g_s10_1,g_s10 ! scratch ga arrays
      integer natoms
      double precision pos(3*natoms)
      integer basis ! basis handle
      logical oskel
      integer NoKinetic
      common /skipKinetic/NoKinetic ! goes to int_giao_1ega()
      external int_giao_1ega,giao_aotomo
      integer debug_addH10
c     bq charges
      integer nbq,nextbq,ncosbq,
     &        l_xyz,k_xyz ! dummy variables
      nbq = 0    ! COSMO-variables
      nextbq = 0 ! COSMO-variables
      ncosbq = 0 ! COSMO-variables
      oskel = .false.

      debug_addH10=0 ! 0= no debug 1= debug

c ----------- Create scratch ga-arrays ------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                   alo,g_s10_1)) 
     &    call errquit('add_H10: nga_create failed g_s10',
     &                  0,GA_ERR)
       call ga_zero(g_s10_1)
      ahi(3) = 3*npol
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                   alo,g_s10)) 
     &    call errquit('add_H10: nga_create failed g_s10',
     &                  0,GA_ERR)
       call ga_zero(g_s10)

      if (.not.(do_zora)) then ! do: HF or DFT (norel calc.)

c       if (ga_nodeid().eq.0)
c    &   write(*,*) 'enter-NonRel-addH10'

        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'l10' ,pos(1),natoms,oskel)
        NoKinetic=0 ! =0 DO-kinetic, =1 SKIP-kinetic
        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'tv10',pos(1),natoms,oskel)
      else                            ! do: zora (relativistic calc.)
c ----------- Create scratch ga-arrays ------- END
        NoKinetic=1 ! =0 DO-kinetic, =1 SKIP-kinetic
        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'tv10',pos(1),natoms,oskel)
        if (debug_addH10.eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '---------g_s10-AO --------- START'
          call ga_print(g_s10_1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---------g_s10-AO --------- END'
        endif

        call ga_add(1.0d0,g_s10_1,
     &              1.0d0, ga_Fji, ! update g_s10_1 with ga_Fji
     &              g_s10_1)       ! out
      endif
c
c     Get external and cosmo bq contribution
      nbq = 0
      nextbq = 0
      ncosbq = 0
      if(geom_extbq_on()) nextbq = geom_extbq_ncenter()
      nbq = nextbq ! external bq's
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &    nbq = ncosbq ! cosmo bq's
      if (nextbq.gt.0.and.ncosbq.gt.0)
     &    nbq = nextbq + ncosbq  ! tally up cosmo and external bqs
      if (nbq.gt.0) then
        call int_giao_1ega(basis,basis,g_s10_1,'bq10',dbl_mb(k_xyz),
     &                     natoms,oskel)
      end if
c
c     Transform H10 to MO and add to g_rhs
c     Copy: g_s10_1 --> g_s10
       blo(1) = 1
       bhi(1) = nbf ! nmo-fixing4lineardependency
       blo(2) = 1
       bhi(2) = nbf ! nmo-fixing4lineardependency
       blo(3) = 1
       bhi(3) = 3
      do ispin=1,npol   
       shift=3*(ispin-1)
       alo(1) = 1
       ahi(1) = nbf ! nmo-fixing4lineardependency
       alo(2) = 1
       ahi(2) = nbf ! nmo-fixing4lineardependency
       alo(3) = shift+1
       ahi(3) = shift+3   
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,3,nbf)
        if (debug_addH10.eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '---------g_s10-MO --------- START'
          call ga_print(vectors)
          call ga_print(g_s10)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---------g_s10-MO --------- END'
        endif
      do ispin=1,npol
       shift=3*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+3   
c ------ definitions for g_rhs -------- START
       disp=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c ------ definitions for g_rhs -------- END
       if (debug_addH10.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,69) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 69      format('addH10-indices::blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
        endif
       endif
       call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs,blo,bhi)
      enddo ! end-loop-ispin
      if (.not.ga_destroy(g_s10_1)) call 
     &    errquit('add_H10: ga_destroy failed g_s10',0,GA_ERR)
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('add_H10: ga_destroy failed g_s10',0,GA_ERR)
      return
      end

      subroutine get_P10(  g_p1, ! out: Perturbed (spin)density matrix
     &                 type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift              
     &                    g_rhs, ! in: accumulated rhs expression
     &                   g_rhs0, ! in: from get_prelim_fock()         
     &                  vectors, ! in: MO                coeffs
     &                 g_CiFull, ! in: MO zora weighting coeffs
     &                      nbf, ! in: nr. basis functions
     &                      nmo, ! in: nr. MOs (occ+virt)
     &                     npol, ! in: nr. of polarizations
     &                     nocc, ! in: nr. occ     MOs
     &                    nvirt, ! in: nr. virtual MOs
     &                  do_zora, ! in: .true.  if doing zora calc
     &                do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &           not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &                     rtdb,
     &              lbl_nlmohyp)
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol   ! nr. of polarizations
      integer g_p1   ! OUT: Perturbed density matrix
      integer g_rhs  ! IN : accumulated right-hand-side expression
      integer g_rhs0 ! IN : from get_prelim_fock()
      integer g_u,g_d1 ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,type_NMR,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      double precision coeff(2)
      logical status,do_zora,do_NonRel,
     &        not_zora_scale,skip_hypAOev
c ------ for Hyperfine NMLO analysis --------- START
      integer hypfile,ntot,g_c1, ! g_c1 , collects pertrbed MO coeffs C1
     &        plo(3),phi(3),qlo(3),qhi(3),
     &        ndata,ndir,ndir1,nlist,
     &        g_munuFCSD,g_munuPSOSO,g_sdens  ! dummy variables not used 
      logical dft_zoraHYP_NLMOAnalysis_write
      character lbl_nlmohyp*(*)
      character*255 zorafilename
      external dft_zoraHYP_NLMOAnalysis_write,
     &         util_file_name
c ------ for Hyperfine NMLO analysis --------- END          
      integer debug_p10

      debug_p10=0 ! =1 for debugging

       if      (type_NMR.eq.1) then ! Shieldings
         coeff(1)= 1.0d0
         coeff(2)= 1.0d0
       else if (type_NMR.eq.2 .or.  ! Hyperfine
     &          type_NMR.eq.3) then ! g-shifts
         coeff(1)=  1.0d0
         coeff(2)= -1.0d0
       else
        write(*,*) 'Error in get_P10_1:',
     &             ' Calc. should be giao, gshift or hyperfine.'
        stop
       endif
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'g_d1-a matrix',alo,g_d1)) call 
     &     errquit('g_d1: nga_create failed g_d1',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1 matrix',alo,g_p1)) call 
     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
       call ga_zero(g_p1)
c +++++++++ store g_u for NMLO analysis +++++++++ START
        if (type_NMR.eq.2) then ! store only if Hyperfine calc.
         if(.not.rtdb_get(rtdb,'zora:skip_hypAOev',
     &                    mt_log,1,skip_hypAOev))        
     &   skip_hypAOev = .false. 
         hypfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
         if (hypfile.eq.1 .and. .not.(skip_hypAOev)) then
          alo(1) = nbf
          alo(2) = -1
          alo(3) = -1
          ahi(1) = nbf
          ntot=nocc(1)+nocc(2)
          ahi(2) = ntot
          ahi(3) = 3
          if (.not.nga_create(MT_DBL,3,ahi,'g_c1 matrix',alo,g_c1)) call 
     &     errquit('g_c1: nga_create failed g_c1',0,GA_ERR)
          call ga_zero(g_c1)
         endif
        endif
c +++++++++ store g_u for NMLO analysis +++++++++ END
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('g_d1: nga_create failed g_u',0,GA_ERR)
       call ga_zero(g_u)
c ----------- copy occ-occ info  ------------ START
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       call nga_copy_patch('n',g_rhs0,blo,bhi,
     &                            g_u,alo,ahi) 
c ----------- copy occ-occ info  ------------ END
c ----------- copy occ-virt info ------------ START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3

       if (debug_p10 .eq.1) then
         write(*,16) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 16      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
       endif
       call nga_copy_patch('n',g_rhs,blo,bhi,
     &                           g_u,alo,ahi)  
c ----------- copy occ-virt info ------------ START
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '------- g_u------- START'
        call ga_print(g_u)
        if (ga_nodeid().eq.0)
     &   write(*,*) '------- g_u------- END'
       endif

c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dhi(1) = nbf
       dlo(2) = 1
       dhi(2) = nocc(ispin)

c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then

c       if (ga_nodeid().eq.0) write(*,*) 'FA-enter-scaling'
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---g_CiFull(',ispin,')------ START'
c       call ga_print(g_CiFull(ispin))
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---g_CiFull(',ispin,')------ END'

        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
        if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
          write(*,17) blo(1),bhi(1),blo(2),bhi(2),
     &                blo(3),bhi(3),
     &                alo(1),ahi(1),alo(2),ahi(2),
     &                alo(3),ahi(3),
     &                dlo(1),dhi(1),dlo(2),dhi(2),
     &                dlo(3),dhi(3)
 17      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
         endif
        endif
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u,alo,ahi,
     &                                  g_d1,dlo,dhi) 
        call nga_copy_patch('n',g_d1,dlo,dhi,
     &                           g_u,dlo,dhi)
c ++++++ store g_u for NMLO analysis in Hyperfine calc.++++ START
        if (type_NMR.eq.2) then ! store only if Hyperfine calc.
         if(.not.rtdb_get(rtdb,'zora:skip_hypAOev',
     &                    mt_log,1,skip_hypAOev))        
     &   skip_hypAOev = .false.    
         hypfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
         if (hypfile.eq.1 .and. .not.(skip_hypAOev)) then
c -------- g_u --> g_c1 ----------START
          shift=nocc(1)*(ispin-1)
          plo(1) = 1
          phi(1) = nbf
          plo(2) = 1
          phi(2) = nocc(ispin)
          plo(3) = xyz
          phi(3) = xyz
          qlo(1) = 1
          qhi(1) = nbf
          qlo(2) = shift+1
          qhi(2) = shift+nocc(ispin)
          qlo(3) = xyz
          qhi(3) = xyz
          call nga_copy_patch('n',g_u ,plo,phi,
     &                            g_c1,qlo,qhi)     
c -------- g_u --> g_c1 ----------END     
         endif
        endif
c ++++++ store g_u for NMLO analysis in Hyperfine calc.++++ END
        if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,30) ispin,xyz
  30     format('------ g_u(',i3,',',i3,')----------- START')
         call ga_print(g_u)
         if (ga_nodeid().eq.0)
     &    write(*,31) ispin,xyz
  31     format('------ g_u(',i3,',',i3,')----------- END')
        endif

        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1
       
       if (debug_p10 .eq.1) then
         write(*,18) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
       endif

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u,alo,ahi,
     &                                  g_d1,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,32) ispin,xyz
  32    format('------ g_d1-1(',i3,',',i3,')----------- START')
        call ga_print(g_d1)
        if (ga_nodeid().eq.0)
     &   write(*,33) ispin,xyz
  33    format('------ g_d1-1(',i3,',',i3,')----------- END')
       endif

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,34) ispin,xyz
  34    format('------ g_d1-2(',i3,',',i3,')----------- START')
        call ga_print(g_d1)
        if (ga_nodeid().eq.0)
     &   write(*,35) ispin,xyz
  35    format('------ g_d1-2(',i3,',',i3,')----------- END')
       endif

       enddo ! end-loop-xyz

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,1) ispin
 1       format('------g_d1(',i3,')---------- START')
         call ga_print(g_d1)
        if (ga_nodeid().eq.0) 
     &   write(*,2) ispin
 2       format('------g_d1(',i3,')---------- END')     
       endif

       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
       call nga_add_patch(1.0d0       ,g_p1,elo,ehi,
     &                    coeff(ispin),g_d1,elo,ehi,
     &                                 g_p1,elo,ehi)
       if (.not.ga_destroy(g_u)) call 
     &    errquit('get_d1: ga_destroy failed g_u',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
      enddo ! end-loop-ispin

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_p1-1---------- START'
         call ga_print(g_p1)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-1---------- END'
       endif

      if (npol.eq.1 .and. type_NMR.eq.1) then ! this happens ONLY for NMR-restricted calc.
        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1,elo,ehi,2.0d0)
      endif

        if (type_NMR.eq.2) then ! store only if Hyperfine calc.
         if(.not.rtdb_get(rtdb,'zora:skip_hypAOev',
     &                    mt_log,1,skip_hypAOev))        
     &   skip_hypAOev = .false.    
         hypfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:hypfile',mt_int,1,hypfile) ! for NLMO analysis
         if (hypfile.eq.1 .and. .not.(skip_hypAOev)) then
c -------- g_c1 --> file ----------START
         ndir1=3
         ndir =6
         ndata=2 !  =1 write FCSD,PSOSO,sdens =2 write g_c1
         call util_file_name(lbl_nlmohyp,.false.,.false.,zorafilename)
         if (.not.dft_zoraHYP_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munuPSOSO
     &              nlist, ! in: dummy
     &               nocc, ! in: used here
     &              ndata, ! in: =1 write FCSD,PSOSO,sdens =2 write vectors,g_c1
     &         g_munuFCSD, ! in: dummy
     &        g_munuPSOSO, ! in: dummy
     &               npol, ! in: nr of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: used here
     &            g_sdens))! in: dummy
     &   call errquit('get_Htensor_fast: dft_zoraHYPNLMO_write failed',
     &                0,DISK_ERR)
c -------- g_c1 --> file ----------END
         if (.not.ga_destroy(g_c1)) call 
     &       errquit('get_P10: ga_destroy failed g_c1',0,GA_ERR)    
         endif
        endif
      
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- START'
        call ga_print(g_p1)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- END'
       endif

       if (.not.ga_destroy(g_d1)) call 
     &    errquit('get_d1: ga_destroy failed g_d1',0,GA_ERR)      
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_P10_1(
     &                     g_p1_oo, ! out: Perturbed (spin)density matrix-occ-occ       contrib
     &                     g_p1_ov, ! out: Perturbed (spin)density matrix-occ-virt      contrib
     &                    type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                       g_rhs, ! in: accumulated rhs expression
     &                      g_rhs0, ! in: from get_prelim_fock()         
     &                     vectors, ! in: MO                coeffs
     &                    g_CiFull, ! in: MO zora weighting coeffs
     &                         nbf, ! in: nr. basis functions
     &                         nmo, ! in: nr. MOs (occ+virt)
     &                        npol, ! in: nr. of polarizations
     &                        nocc, ! in: nr. occ     MOs
     &                       nvirt, ! in: nr. virtual MOs
     &                     do_zora, ! in: .true.  if doing zora calc
     &                   do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &              not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &              lbl_nlmogshift, ! in: for g-shift nlmo analysis
     &              lbl_nlmoshield, ! in: for shield  nlmo analysis
     &                        rtdb)
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      
      integer g_p1 ! for debugging

      integer npol   ! nr. of polarizations
      integer g_p1_oo,   ! OUT: Perturbed density matrix occ-occ  contrib
     &        g_p1_ov    ! OUT: Perturbed density matrix occ-virt contrib 
      integer g_rhs  ! IN : accumulated right-hand-side expression
      integer g_rhs0 ! IN : from get_prelim_fock()
      logical do_zora,do_NonRel,not_zora_scale
      integer g_u_oo,g_u_ov,g_d1_oo,g_d1_ov ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz,type_NMR
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      double precision coeff(2)
      logical status,
     &        skip_gshiftAOev,
     &        skip_csAOev
c ------ for g-shift NLMO analysis --------- START
      integer gshiftfile,ntot,g_c1,  ! g_c1 , collects perturbed MO coeffs C1
     &        plo(3),phi(3),qlo(3),qhi(3),
     &        ndata,ndir,ndir1,nlist,
     &        g_munuEPRdia,g_munuEPRpar1, ! dummy variables not used here
     &        g_munuEPRHpar,              ! dummy variables not used here
     &        g_sdens                     ! dummy variables not used here
      logical dft_zoraGshift_NLMOAnalysis_write
      character lbl_nlmogshift*(*)
      character*255 zorafilename
      external dft_zoraGshift_NLMOAnalysis_write,
     &         util_file_name
c ------ for g-shift NLMO analysis --------- END    
c ------ for shield  NLMO analysis --------- START
      character lbl_nlmoshield*(*)
      integer shldfile
      integer g_munudia,  ! in: dummy variable not used
     &        g_munupar1, ! in: dummy variable not used
     &        g_munu_h01, ! in: dummy variable not used
     &        g_dens      ! in: dummy not used here
      logical dft_zoraShield_NLMOAnalysis_write
      external dft_zoraShield_NLMOAnalysis_write
c ------ for shield  NLMO analysis --------- END  
         
      integer debug_p10

       debug_p10=0 ! =1 for debugging

       if      (type_NMR.eq.1) then ! Shieldings
         coeff(1)=  1.0d0
         coeff(2)=  1.0d0
       else if (type_NMR.eq.2 .or.  ! Hyperfine
     &          type_NMR.eq.3) then ! g-shifts
         coeff(1)=  1.0d0
         coeff(2)= -1.0d0
       else
        write(*,*) 'Error in get_P10_1:',
     &             ' Calc. should be giao, gshift or hyperfine.'
        stop
       endif
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'g_d1_oo matrix',
     &      alo,g_d1_oo)) call 
     &     errquit('g_d1_oo: nga_create failed g_d1_oo',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1_oo matrix',
     &      alo,g_p1_oo)) call 
     &     errquit('g_p1_oo: nga_create failed g_p1_oo',0,GA_ERR)
       call ga_zero(g_p1_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'g_d1_ov matrix',
     &      alo,g_d1_ov)) call 
     &     errquit('g_d1_ov: nga_create failed g_d1_ov',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1_ov matrix',
     &      alo,g_p1_ov)) call 
     &     errquit('g_p1_ov: nga_create failed g_p1_ov',0,GA_ERR)
       call ga_zero(g_p1_ov)
c ------ For debugging ---- START
c       if (debug_p10 .eq.1) then
c        if (.not.nga_create(MT_DBL,3,ahi,'g_p1-b matrix',
c     &      alo,g_p1)) call 
c     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
c        call ga_zero(g_p1)
c       endif
c ------ For debugging ---- END
c +++++++++ store g_u for NMLO analysis +++++++++ START
        if (type_NMR.eq.1 .or. type_NMR.eq.3) then ! store only if Hyperfine calc.
         if(.not.rtdb_get(rtdb,'zora:skip_csAOev',
     &                    mt_log,1,skip_csAOev))        
     &   skip_csAOev = .false.    
         if(.not.rtdb_get(rtdb,'zora:skip_gshiftAOev',
     &                    mt_log,1,skip_gshiftAOev))        
     &   skip_gshiftAOev = .false.    
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if ((gshiftfile.eq.1 .and. .not.(skip_gshiftAOev)) .or. 
     &       (shldfile  .eq.1 .and. .not.(skip_csAOev))) then
          alo(1) = nbf
          alo(2) = -1
          alo(3) = -1
          ahi(1) = nbf
          ntot=nocc(1)+nocc(2)
          ahi(2) = ntot
          ahi(3) = 3
          if (.not.nga_create(MT_DBL,3,ahi,'g_c1 matrix',alo,g_c1)) call 
     &     errquit('g_c1: nga_create failed g_c1',0,GA_ERR)
          call ga_zero(g_c1)
         endif
        endif
c +++++++++ store g_u for NMLO analysis +++++++++ END
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_oo)) call 
     &    errquit('g_u_oo: nga_create failed g_u_oo',0,GA_ERR)
       call ga_zero(g_u_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_ov)) call 
     &    errquit('g_u_ov: nga_create failed g_u_ov',0,GA_ERR)
       call ga_zero(g_u_ov)
c ----------- copy occ-occ info  ------------ START
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       call nga_copy_patch('n',g_rhs0,blo,bhi,
     &                         g_u_oo,alo,ahi) 
c ----------- copy occ-occ info  ------------ END
c ----------- copy occ-virt info ------------ START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,16) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 16      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
        endif
       endif
       call nga_copy_patch('n',g_rhs,blo,bhi,
     &                           g_u_ov,alo,ahi)  
c ----------- copy occ-virt info ------------ START
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_u_oo------- START'
        call ga_print(g_u_oo)
         if (ga_nodeid().eq.0) then
          write(*,*) '------- g_u_oo------- END'
          write(*,*) '------- g_u_ov------- START'
         endif
        call ga_print(g_u_ov)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------- g_u_ov------- END'
       endif

c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dlo(2) = 1
       dhi(1) = nbf
       dhi(2) = nocc(ispin)
c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then
        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,*) 'bef-1-matmul-patch'
         write(*,17) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               dlo(1),dhi(1),dlo(2),dhi(2),
     &               dlo(3),dhi(3)
 17      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
         endif
        endif
c ------ get:  g_u_oo ------------------------ START
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,dlo,dhi) 
        call nga_copy_patch('n',g_d1_oo,dlo,dhi,
     &                           g_u_oo,dlo,dhi)
c ------ get:  g_u_oo ------------------------ END
c ------ get:  g_u_ov ------------------------ START
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov,dlo,dhi,
     &                           g_u_ov,dlo,dhi)
c ------ get:  g_u_ov ------------------------ END
c ++++++ store g_u for NMLO analysis in g-shift calc.++++ START
        if (type_NMR.eq.1 .or. type_NMR.eq.3) then ! store only if g-shift calc.
         if(.not.rtdb_get(rtdb,'zora:skip_csAOev',
     &                    mt_log,1,skip_csAOev))        
     &   skip_csAOev = .false.    
         if(.not.rtdb_get(rtdb,'zora:skip_gshiftAOev',
     &                    mt_log,1,skip_gshiftAOev))        
     &   skip_gshiftAOev = .false.    
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if ((gshiftfile.eq.1 .and. .not.(skip_gshiftAOev)) .or. 
     &       (shldfile  .eq.1 .and. .not.(skip_csAOev))) then
c -------- g_u --> g_c1 ----------START
          shift=nocc(1)*(ispin-1)
          plo(1) = 1
          phi(1) = nbf
          plo(2) = 1
          phi(2) = nocc(ispin)
          plo(3) = xyz
          phi(3) = xyz
          qlo(1) = 1
          qhi(1) = nbf
          qlo(2) = shift+1
          qhi(2) = shift+nocc(ispin)
          qlo(3) = xyz
          qhi(3) = xyz
          call nga_copy_patch('n',g_u_oo,plo,phi,
     &                            g_c1  ,qlo,qhi)     
          call nga_add_patch(1.0d0,g_c1  ,qlo,qhi,
     &                       1.0d0,g_u_ov,plo,phi,
     &                             g_c1  ,qlo,qhi)
c -------- g_u --> g_c1 ----------END     
         endif
        endif
c ++++++ store g_u for NMLO analysis in g-shift calc.++++ END

       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &     write(*,30) ispin,xyz
  30     format('------ g_u(',i3,',',i3,')----------- START')
         call ga_print(g_u_oo)
         call ga_print(g_u_ov)
         if (ga_nodeid().eq.0)
     &    write(*,31) ispin,xyz
  31     format('------ g_u(',i3,',',i3,')----------- END')
        endif

        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1
       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,18) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
        endif
       endif

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,32) ispin,xyz
  32    format('------ g_d1-1(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,33) ispin,xyz
  33    format('------ g_d1-1(',i3,',',i3,')----------- END')
       endif

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_oo,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,34) ispin,xyz
  34    format('------ g_d1-2(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,35) ispin,xyz
  35    format('------ g_d1-2(',i3,',',i3,')----------- END')
       endif

       enddo ! end-loop-xyz

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,1) ispin
 1       format('------g_d1(',i3,')---------- START')
         call ga_print(g_d1_oo)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- START'
         call ga_print(g_d1_ov)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- END'
        if (ga_nodeid().eq.0) 
     &   write(*,2) ispin
 2       format('------g_d1(',i3,')---------- END')     
       endif

       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
c       write(*,*) 'ispin=',ispin,'coeff=',coeff(ispin)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--BEF--ST'
c        call ga_print(g_p1_oo)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--BEF--EN'
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--BEF--ST'
c        call ga_print(g_p1_ov)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--BEF--EN'
       call nga_add_patch(1.0d0       ,g_p1_oo,elo,ehi,
     &                    coeff(ispin),g_d1_oo,elo,ehi,
     &                                 g_p1_oo,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov,elo,ehi,
     &                    coeff(ispin),g_d1_ov,elo,ehi,
     &                                 g_p1_ov,elo,ehi)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--AFT--ST'
c        call ga_print(g_p1_oo)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--AFT--EN'
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--AFT--ST'
c        call ga_print(g_p1_ov)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--AFT--EN'

       if (.not.ga_destroy(g_u_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oo',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
      enddo ! end-loop-ispin

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_p1-1---------- START'
         call ga_print(g_p1_oo)
         call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-1---------- END'
       endif

      if (npol.eq.1 .and. type_NMR.eq.1) then ! this happens ONLY for NMR-restricted calc.
        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1_oo,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov,elo,ehi,2.0d0)
      endif

        if (type_NMR.eq.1) then ! store only if shield calc.
         if(.not.rtdb_get(rtdb,'zora:skip_csAOev',
     &                    mt_log,1,skip_csAOev))        
     &   skip_csAOev = .false.    
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         if (shldfile.eq.1 .and. .not.(skip_csAOev)) then
c -------- g_c1 --> file ----------START
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '----------BEF-writing: g_c1----------- START'
c         call ga_print(g_c1)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '----------BEF-writing: g_c1----------- END'
         ndir1=3
         ndir =6
         ndata=2 !  =1 write g_munudia,g_munupar1,g_munu_h01 sdens =2 write g_c1
         call util_file_name(lbl_nlmoshield,.false.,.false.,
     &                       zorafilename)
         if (.not.dft_zoraShield_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munudia,g_munupar1
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munu_h01
     &              nlist, ! in: list of selected atoms
     &               nocc, ! in: used here
     &              ndata, ! in: =1 write dia,par1,sdens =2 write g_c1
     &          g_munudia, ! in: dummy variable not used
     &         g_munupar1, ! in: dummy variable not used
     &         g_munu_h01, ! in: dummy variable not used
     &               npol, ! in: nr. of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: used here
     &             g_dens))! in: dummy not used here
     &   call errquit('get_P10_1: dft_zoraShieldNLMO_write failed',
     &                0,DISK_ERR)
c -------- g_c1 --> file ----------END
         if (.not.ga_destroy(g_c1)) call 
     &       errquit('get_P10_1: ga_destroy failed g_c1',0,GA_ERR)    
         endif
        endif
        if (type_NMR.eq.3) then ! store only if g-shift calc.
         if(.not.rtdb_get(rtdb,'zora:skip_gshiftAOev',
     &                    mt_log,1,skip_gshiftAOev))        
     &   skip_gshiftAOev = .false.    
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if (gshiftfile.eq.1 .and. .not.(skip_gshiftAOev)) then
c -------- g_c1 --> file ----------START
          ndir1=3
          ndir =6
          ndata=2 !  =1 write g_munuEPRdia,g_munuEPRpar1,g_munuEPRHpar sdens =2 write g_c1
          call util_file_name(lbl_nlmogshift,.false.,.false.,
     &                        zorafilename)
          if (.not.dft_zoraGshift_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuEPRdia,g_munuEPRpar1
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munuEPRHpar
     &               nocc, ! in: used here
     &              ndata, ! in: =1 write dia,par1,sdens =2 write g_c1
     &       g_munuEPRdia, ! in: dummy not used here
     &      g_munuEPRpar1, ! in: dummy not used here
     &      g_munuEPRHpar, ! in: dummy not used here
     &               npol, ! in: nr. of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: used here
     &            g_sdens))! in: dummy not used here
     &    call errquit('get_P10_1: dft_zoraGshiftNLMO_write failed',
     &                 0,DISK_ERR)
c -------- g_c1 --> file ----------END
         if (.not.ga_destroy(g_c1)) call 
     &       errquit('get_P10_1: ga_destroy failed g_c1',0,GA_ERR)    
         endif
        endif
      
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- START'
        call ga_print(g_p1_oo)
        call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- END'
c       if (ga_nodeid().eq.0) 
c     &    write(*,*) '------g_p1-tot---------- START'
c        call ga_add(1.0d0,g_p1_oo,1.0d0,g_p1_ov,g_p1)
c        call ga_print(g_p1)
c       if (ga_nodeid().eq.0) 
c     &    write(*,*) '------g_p1-tot---------- END'
       endif

       if (.not.ga_destroy(g_d1_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oo',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov',0,GA_ERR)    
      return
      end
c ----------------------------------------
c ---- 04-26-11 ----- get A or B --- START
      subroutine get_P10_1_AorB(
     &                     g_p1_oo, ! out: Perturbed (spin)density matrix-occ-occ       contrib
     &                     g_p1_ov, ! out: Perturbed (spin)density matrix-occ-virt      contrib
     &                    type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                       g_rhs, ! in: accumulated rhs expression
     &                      g_rhs0, ! in: from get_prelim_fock()         
     &                     vectors, ! in: MO                coeffs
     &                    g_CiFull, ! in: MO zora weighting coeffs
     &                         nbf, ! in: nr. basis functions
     &                         nmo, ! in: nr. MOs (occ+virt)
     &                       ispin, ! in: select polarization = 1(A) or 2(B)
     &                        npol, ! in: nr. of polarizations
     &                        nocc, ! in: nr. occ     MOs
     &                       nvirt, ! in: nr. virtual MOs
     &                     do_zora, ! in: .true.  if doing zora calc
     &                   do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &              not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &                        rtdb)
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      
      integer g_p1 ! for debugging
      integer npol   ! nr. of polarizations
      integer g_p1_oo,   ! OUT: Perturbed density matrix occ-occ  contrib
     &        g_p1_ov    ! OUT: Perturbed density matrix occ-virt contrib 
      integer g_rhs  ! IN : accumulated right-hand-side expression
      integer g_rhs0 ! IN : from get_prelim_fock()
      logical do_zora,do_NonRel,not_zora_scale
      integer g_u_oo,g_u_ov,g_d1_oo,g_d1_ov ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz,type_NMR
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      logical status
         
      integer debug_p10

       debug_p10=0 ! =1 for debugging
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'g_d1_ooAB matrix',
     &       alo,g_d1_oo)) call 
     &     errquit('g_d1_oo: nga_create failed g_d1_oo',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1_ooAB matrix',
     &       alo,g_p1_oo)) call 
     &     errquit('g_p1_oo: nga_create failed g_p1_oo',0,GA_ERR)
       call ga_zero(g_p1_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'g_d1_ovAB matrix',
     &       alo,g_d1_ov)) call 
     &     errquit('g_d1_ov: nga_create failed g_d1_ov',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1_ovAB matrix',
     &       alo,g_p1_ov)) call 
     &     errquit('g_p1_ov: nga_create failed g_p1_ov',0,GA_ERR)
       call ga_zero(g_p1_ov)
c ------ For debugging ---- START
       if (.not.nga_create(MT_DBL,3,ahi,'g_p1AB matrix',
     &       alo,g_p1)) call 
     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
       call ga_zero(g_p1)
c ------ For debugging ---- END

c      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_oo)) call 
     &    errquit('g_u_oo: nga_create failed g_u_oo',0,GA_ERR)
       call ga_zero(g_u_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_ov)) call 
     &    errquit('g_u_ov: nga_create failed g_u_ov',0,GA_ERR)
       call ga_zero(g_u_ov)
c ----------- copy occ-occ info  ------------ START
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       call nga_copy_patch('n',g_rhs0,blo,bhi,
     &                            g_u_oo,alo,ahi) 
c ----------- copy occ-occ info  ------------ END
c ----------- copy occ-virt info ------------ START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,16) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 16      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
        endif
       endif
       call nga_copy_patch('n',g_rhs,blo,bhi,
     &                           g_u_ov,alo,ahi)  
c ----------- copy occ-virt info ------------ START
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_u_oo------- START'
        call ga_print(g_u_oo)
         if (ga_nodeid().eq.0) then
          write(*,*) '------- g_u_oo------- END'
          write(*,*) '------- g_u_ov------- START'
         endif
        call ga_print(g_u_ov)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------- g_u_ov------- END'
       endif

c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dlo(2) = 1
       dhi(1) = nbf
       dhi(2) = nocc(ispin)
c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then
        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,*) 'bef-1-matmul-patch'
         write(*,17) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               dlo(1),dhi(1),dlo(2),dhi(2),
     &               dlo(3),dhi(3)
 17      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
         endif
        endif
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,dlo,dhi) 
        call nga_copy_patch('n',g_d1_oo,dlo,dhi,
     &                           g_u_oo,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov,dlo,dhi,
     &                           g_u_ov,dlo,dhi)

       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &     write(*,30) ispin,xyz
  30     format('------ g_u(',i3,',',i3,')----------- START')
         call ga_print(g_u_oo)
         call ga_print(g_u_ov)
         if (ga_nodeid().eq.0)
     &    write(*,31) ispin,xyz
  31     format('------ g_u(',i3,',',i3,')----------- END')
        endif

        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1
       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,18) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
        endif
       endif

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,32) ispin,xyz
  32    format('------ g_d1-1(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,33) ispin,xyz
  33    format('------ g_d1-1(',i3,',',i3,')----------- END')
       endif

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_oo,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,34) ispin,xyz
  34    format('------ g_d1-2(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,35) ispin,xyz
  35    format('------ g_d1-2(',i3,',',i3,')----------- END')
       endif

       enddo ! end-loop-xyz

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,1) ispin
 1       format('------g_d1(',i3,')---------- START')
         call ga_print(g_d1_oo)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- START'
         call ga_print(g_d1_ov)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- END'
        if (ga_nodeid().eq.0) 
     &   write(*,2) ispin
 2       format('------g_d1(',i3,')---------- END')     
       endif

       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
c       write(*,*) 'ispin=',ispin,'coeff=',coeff(ispin)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--BEF--ST'
c        call ga_print(g_p1_oo)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--BEF--EN'
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--BEF--ST'
c        call ga_print(g_p1_ov)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--BEF--EN'
       call nga_add_patch(1.0d0,g_p1_oo,elo,ehi,
     &                    1.0d0,g_d1_oo,elo,ehi,
     &                          g_p1_oo,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov,elo,ehi,
     &                    1.0d0,g_d1_ov,elo,ehi,
     &                          g_p1_ov,elo,ehi)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--AFT--ST'
c        call ga_print(g_p1_oo)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_oo--AFT--EN'
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--AFT--ST'
c        call ga_print(g_p1_ov)
c       if (ga_nodeid().eq.0) write(*,*) '---g_p1_ov--AFT--EN'

       if (.not.ga_destroy(g_u_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oo',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
c      enddo ! end-loop-ispin

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_p1-1---------- START'
         call ga_print(g_p1_oo)
         call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-1---------- END'
       endif

      if (npol.eq.1 .and. type_NMR.eq.1) then ! this happens ONLY for NMR-restricted calc.
        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1_oo,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov,elo,ehi,2.0d0)
      endif
      
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- START'
        call ga_print(g_p1_oo)
        call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- END'
       if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-tot---------- START'
        call ga_add(1.0d0,g_p1_oo,1.0d0,g_p1_ov,g_p1)
        call ga_print(g_p1)
       if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-tot---------- END'
       endif

       if (.not.ga_destroy(g_d1_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oo',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov',0,GA_ERR)    
       if (.not.ga_destroy(g_p1)) call 
     &    errquit('get_d1: ga_destroy failed g_p1AB',0,GA_ERR)   
      return
      end
c ---- 04-26-11 ----- get A or B --- END
c ----------------------------------------
c ++++++++++++++++++++++++++++++++++++++++++++
c =======================================================
c +++++++++++++++ for Coul,Exch test ++++++++++++++ START
      subroutine get_P10_JK(
     &                     g_p1_oo, ! out: Perturbed (spin)density matrix-occ-occ       contrib
     &                     g_p1_ov, ! out: Perturbed (spin)density matrix-occ-virt      contrib
     &                    type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                g_p1_ov_Coul, ! out: Perturbed (spin)density matrix with Coulomb  contrib
     &                g_p1_ov_Exch, ! out: Perturbed (spin)density matrix with Exchange contrib (fock_xc)
     &                g_p1_ov_noJK, ! out: Perturbed (spin)density matrix remaining terms not J or K
     &                  g_p1_ov_1e, ! out: Perturbed (spin)density matrix 1e-operator contrib
     &                g_p1_ov_eSji, ! out: Perturbed (spin)density matrix perturbed overlap matrix
     &                       g_rhs, ! in: accumulated rhs expression
     &                  g_rhs_Coul,
     &                  g_rhs_Exch,
     &                  g_rhs_noJK,
     &                    g_rhs_1e,
     &                  g_rhs_eSji,
     &                      g_rhs0, ! in: from get_prelim_fock()         
     &                     vectors, ! in: MO                coeffs
     &                    g_CiFull, ! in: MO zora weighting coeffs
     &                         nbf, ! in: nr. basis functions
     &                         nmo, ! in: nr. MOs (occ+virt)
     &                        npol, ! in: nr. of polarizations
     &                        nocc, ! in: nr. occ     MOs
     &                       nvirt, ! in: nr. virtual MOs
     &                     do_zora, ! in: .true.  if doing zora calc
     &                   do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &              not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &              lbl_nlmogshift, ! in: for g-shift nlmo analysis
     &              lbl_nlmoshield, ! in: for shield  nlmo analysis
     &                        rtdb)
c  g_p1_ov      = g_p1_ov_Coul + g_p1_ov_Exch + g_p1_ov_noJK
c  g_p1_ov_noJK = g_p1_ov_1e   + g_p1_ov_eSji
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      
      integer g_p1 ! for debugging

      integer npol   ! nr. of polarizations
      integer g_p1_oo,   ! OUT: Perturbed density matrix occ-occ  contrib
     &        g_p1_ov    ! OUT: Perturbed density matrix occ-virt contrib 
      integer g_rhs  ! IN : accumulated right-hand-side expression
      integer g_rhs0 ! IN : from get_prelim_fock()
      integer g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_1e,g_rhs_eSji,
     &        g_p1_ov_Coul,g_p1_ov_Exch,g_p1_ov_noJK,
     &        g_p1_ov_1e,g_p1_ov_eSji,
     &        g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji,
     &        g_u_ov_Coul,g_u_ov_Exch,g_u_ov_noJK,
     &        g_u_ov_1e,g_u_ov_eSji
      logical do_zora,do_NonRel,not_zora_scale
      integer g_u_oo,g_u_ov,g_d1_oo,g_d1_ov ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,type_NMR,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz,slc_spinpolAO
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      double precision coeff(2)
      logical status
c ------ for g-shift NMLO analysis --------- START
      integer gshiftfile,ntot,g_c1,  ! g_c1 , collects perturbed MO coeffs C1
     &        plo(3),phi(3),qlo(3),qhi(3),
     &        ndata,ndir,ndir1,
     &        g_munuEPRdia,g_munuEPRpar1, ! dummy variables not used here
     &        g_munuEPRHpar,              ! dummy variables not used here
     &        g_sdens                     ! dummy variables not used here
      logical dft_zoraGshift_NLMOAnalysis_write
      character lbl_nlmogshift*(*)
      character*255 zorafilename
      external dft_zoraGshift_NLMOAnalysis_write,
     &         util_file_name
c ------ for g-shift NMLO analysis --------- END   
c ------ for shield  NLMO analysis --------- START
      character lbl_nlmoshield*(*)
      integer shldfile,nlist
      integer g_munudia,  ! in: dummy variable not used
     &        g_munupar1, ! in: dummy variable not used
     &        g_munu_h01, ! in: dummy variable not used
     &        g_dens      ! in: dummy not used here
      logical dft_zoraShield_NLMOAnalysis_write
      external dft_zoraShield_NLMOAnalysis_write
c ------ for shield  NLMO analysis --------- END  
      integer nogshift
      data nogshift/1/            
      integer nogiao             
      data nogiao/1/             
      integer debug_p10

      if(.not.rtdb_get(rtdb,'zora:slc_spinpolAO',        ! FA
     &                 mt_int,1,slc_spinpolAO))          ! FA
     &  slc_spinpolAO= 0 ! 0=A-B,1=A,2=-B selecting spin-density matrix      

      debug_p10=0 ! =1 for debugging

       if      (type_NMR.eq.1) then ! Shieldings
              if (slc_spinpolAO .eq. 0) then !----- START-slc-spindensm
          coeff(1)=  1.0d0
          coeff(2)=  1.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC A+B in P^{10}'
         else if (slc_spinpolAO .eq. 1) then
          coeff(1)=  1.0d0
          coeff(2)=  0.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC Alpha in P^{10}'
         else if (slc_spinpolAO .eq. 2) then
          coeff(1)=  0.0d0
          coeff(2)=  1.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC Beta in P^{10}'
         else
          write(*,22) slc_spinpolAO
 22       format('Error in get_P10_JK:: slc_spinpolAO=',i3,
     &           'It should be: 0,1 or 2')
          stop
        endif !------- ----------------------------- END-slc-spindensm
       else if (type_NMR.eq.2 .or.  ! Hyperfine
     &          type_NMR.eq.3) then ! g-shifts
              if (slc_spinpolAO .eq. 0) then !----- START-slc-spindensm
          coeff(1)=  1.0d0
          coeff(2)= -1.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC A-B in P^{10}'
         else if (slc_spinpolAO .eq. 1) then
          coeff(1)=  1.0d0
          coeff(2)=  0.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC Alpha in P^{10}'
         else if (slc_spinpolAO .eq. 2) then
          coeff(1)=  0.0d0
          coeff(2)= -1.0d0
          if (ga_nodeid().eq.0) 
     &     write(*,*) 'WARNING: SLC neg Beta in P^{10}'
         else
          write(*,12) slc_spinpolAO
 12       format('Error in get_P10_JK:: slc_spinpolAO=',i3,
     &           'It should be: 0,1 or 2')
          stop
        endif !------- ----------------------------- END-slc-spindensm
       else
        write(*,*) 'Error in get_P10_1:',
     &             ' Calc. should be giao, gshift or hyperfine.'
        stop
       endif
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'d1oo matrix',alo,g_d1_oo)) call 
     &     errquit('g_d1_oo: nga_create failed g_d1_oo',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'p1oo matrix',alo,g_p1_oo)) call 
     &     errquit('g_p1_oo: nga_create failed g_p1_oo',0,GA_ERR)
       call ga_zero(g_p1_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'d1ov matrix',alo,g_d1_ov)) call 
     &     errquit('g_d1_ov: nga_create failed g_d1_ov',0,GA_ERR)

       if (.not.nga_create(MT_DBL,3,ahi,'d1ov_Coul matrix',
     &     alo,g_d1_ov_Coul)) call 
     &     errquit('g_d1_ovJ: nga_create failed g_d1_ovJ',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'d1ov_Exch matrix',
     &     alo,g_d1_ov_Exch)) call 
     &     errquit('g_d1_ovK: nga_create failed g_d1_ovK',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'d1ov_noJK matrix',
     &     alo,g_d1_ov_noJK)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_ovnJK',
     &     0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'d1ov_1e matrix',
     &     alo,g_d1_ov_1e)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_ov1e',
     &     0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'d1ov_eSji matrix',
     &     alo,g_d1_ov_eSji)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_oveSji',
     &     0,GA_ERR)

       if (.not.nga_create(MT_DBL,3,ahi,'p1ov matrix',alo,g_p1_ov)) call 
     &     errquit('g_p1_ov: nga_create failed g_p1_ov',0,GA_ERR)
       call ga_zero(g_p1_ov)

       if (.not.nga_create(MT_DBL,3,ahi,'p1ov_Coul matrix',
     &     alo,g_p1_ov_Coul)) call 
     &     errquit('g_p1_ovJ: nga_create failed g_p1_ovJ',0,GA_ERR)
       call ga_zero(g_p1_ov_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'p1ov_Exch matrix',
     &     alo,g_p1_ov_Exch)) call 
     &     errquit('g_p1_ovK: nga_create failed g_p1_ovK',0,GA_ERR)
       call ga_zero(g_p1_ov_Exch)
       if (.not.nga_create(MT_DBL,3,ahi,'p1ov_noJK matrix',
     &     alo,g_p1_ov_noJK)) call 
     &     errquit('g_p1_ovnJK: nga_create failed g_p1_ovnJK',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_noJK)
       if (.not.nga_create(MT_DBL,3,ahi,'p1ov_1e matrix',
     &     alo,g_p1_ov_1e)) call 
     &     errquit('g_p1_ov1e: nga_create failed g_p1_ov1e',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_1e)
       if (.not.nga_create(MT_DBL,3,ahi,'p1ov_eSji matrix',
     &     alo,g_p1_ov_eSji)) call 
     &     errquit('g_p1_oveSji: nga_create failed g_p1_oveSji',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_eSji)
c ------ For debugging ---- START
       if (.not.nga_create(MT_DBL,3,ahi,'p1 matrix',alo,g_p1)) call 
     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
       call ga_zero(g_p1)
c ------ For debugging ---- END
c +++++++++ store g_u for NMLO analysis +++++++++ START
        if (type_NMR.eq.1 .or. type_NMR.eq.3) then ! store only if Hyperfine calc.
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if (gshiftfile.eq.1 .or. shldfile.eq.1) then
          alo(1) = nbf
          alo(2) = -1
          alo(3) = -1
          ahi(1) = nbf
          ntot=nocc(1)+nocc(2)
          ahi(2) = ntot
          ahi(3) = 3
          if (.not.nga_create(MT_DBL,3,ahi,'g_c1 matrix',alo,g_c1)) call 
     &     errquit('g_c1: nga_create failed g_c1',0,GA_ERR)
          call ga_zero(g_c1)
         endif
        endif
c +++++++++ store g_u for NMLO analysis +++++++++ END
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_oo)) call 
     &    errquit('g_u_oo: nga_create failed g_u_oo',0,GA_ERR)
       call ga_zero(g_u_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_ov)) call 
     &    errquit('g_u_ov: nga_create failed g_u_ov',0,GA_ERR)
       call ga_zero(g_u_ov)

       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_Coul)) call 
     &    errquit('g_u_ovJ: nga_create failed g_u_ovJ',0,GA_ERR)
       call ga_zero(g_u_ov_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_Exch)) call 
     &    errquit('g_u_ovK: nga_create failed g_u_ovK',0,GA_ERR)
       call ga_zero(g_u_ov_Exch)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_noJK)) call 
     &    errquit('g_u_ovnJK: nga_create failed g_u_ovnJK',0,GA_ERR)
       call ga_zero(g_u_ov_noJK)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_1e)) call 
     &    errquit('g_u_ov1e: nga_create failed g_u_ov1e',0,GA_ERR)
       call ga_zero(g_u_ov_1e)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_eSji)) call 
     &    errquit('g_u_oveSji: nga_create failed g_u_oveSji',
     &    0,GA_ERR)
       call ga_zero(g_u_ov_eSji)
c ----------- copy occ-occ info  ------------ START
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       call nga_copy_patch('n',g_rhs0,blo,bhi,
     &                         g_u_oo,alo,ahi) 
c ----------- copy occ-occ info  ------------ END
c ----------- copy occ-virt info ------------ START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,16) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 16      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
        endif
       endif
       call nga_copy_patch('n',g_rhs,blo,bhi,
     &                           g_u_ov,alo,ahi)  

       call nga_copy_patch('n',g_rhs_Coul,blo,bhi,
     &                           g_u_ov_Coul,alo,ahi)  
       call nga_copy_patch('n',g_rhs_Exch,blo,bhi,
     &                           g_u_ov_Exch,alo,ahi)  
       call nga_copy_patch('n',g_rhs_noJK,blo,bhi,
     &                           g_u_ov_noJK,alo,ahi)  
       call nga_copy_patch('n',g_rhs_1e,blo,bhi,
     &                           g_u_ov_1e,alo,ahi)  
       call nga_copy_patch('n',g_rhs_eSji,blo,bhi,
     &                           g_u_ov_eSji,alo,ahi)  
c ----------- copy occ-virt info ------------ START
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_u_oo------- START'
        call ga_print(g_u_oo)
         if (ga_nodeid().eq.0) then
          write(*,*) '------- g_u_oo------- END'
          write(*,*) '------- g_u_ov------- START'
         endif
        call ga_print(g_u_ov)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------- g_u_ov------- END'
       endif

c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dlo(2) = 1
       dhi(1) = nbf
       dhi(2) = nocc(ispin)
c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then

c       if (ga_nodeid().eq.0) write(*,*) 'FA-enter-scaling'
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---g_CiFull(',ispin,')------ START'
c       call ga_print(g_CiFull(ispin))
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---g_CiFull(',ispin,')------ END'

        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,*) 'bef-1-matmul-patch'
         write(*,17) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               dlo(1),dhi(1),dlo(2),dhi(2),
     &               dlo(3),dhi(3)
 17      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
         endif
        endif
c ---------- get g_u_oo --------------------START
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,dlo,dhi) 
        call nga_copy_patch('n',g_d1_oo,dlo,dhi,
     &                           g_u_oo,dlo,dhi)
c ---------- get g_u_oo --------------------END
c ---------- get g_u_ov --------------------START
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov,dlo,dhi,
     &                           g_u_ov,dlo,dhi)
c ---------- get g_u_ov --------------------END
c ++++++ store g_u for NMLO analysis in shield,g-shift calc.++++ START
        if (type_NMR.eq.1 .or. type_NMR.eq.3) then ! store only if g-shift calc.
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if (gshiftfile.eq.1 .or. shldfile.eq.1) then
c -------- g_u --> g_c1 ----------START
          shift=nocc(1)*(ispin-1)
          plo(1) = 1
          phi(1) = nbf
          plo(2) = 1
          phi(2) = nocc(ispin)
          plo(3) = xyz
          phi(3) = xyz
          qlo(1) = 1
          qhi(1) = nbf
          qlo(2) = shift+1
          qhi(2) = shift+nocc(ispin)
          qlo(3) = xyz
          qhi(3) = xyz
          call nga_copy_patch('n',g_u_oo,plo,phi,
     &                            g_c1  ,qlo,qhi)     
          call nga_add_patch(1.0d0,g_c1  ,qlo,qhi,
     &                       1.0d0,g_u_ov,plo,phi,
     &                             g_c1  ,qlo,qhi)
c -------- g_u --> g_c1 ----------END     
         endif
        endif
c ++++++ store g_u for NMLO analysis in shield,g-shift calc.++++ END

        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_Coul,alo,ahi,
     &                                  g_d1_ov_Coul,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_Coul,dlo,dhi,
     &                           g_u_ov_Coul,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_Exch,alo,ahi,
     &                                  g_d1_ov_Exch,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_Exch,dlo,dhi,
     &                           g_u_ov_Exch,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_noJK,alo,ahi,
     &                                  g_d1_ov_noJK,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_noJK,dlo,dhi,
     &                           g_u_ov_noJK,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_1e,alo,ahi,
     &                                  g_d1_ov_1e,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_1e,dlo,dhi,
     &                           g_u_ov_1e,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_eSji,alo,ahi,
     &                                  g_d1_ov_eSji,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_eSji,dlo,dhi,
     &                           g_u_ov_eSji,dlo,dhi)
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &     write(*,30) ispin,xyz
  30     format('------ g_u(',i3,',',i3,')----------- START')
         call ga_print(g_u_oo)
         call ga_print(g_u_ov)
         if (ga_nodeid().eq.0)
     &    write(*,31) ispin,xyz
  31     format('------ g_u(',i3,',',i3,')----------- END')
        endif

        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1
       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,18) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
        endif
       endif

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,clo,chi)

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_Coul,alo,ahi,
     &                                  g_d1_ov_Coul,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_Exch,alo,ahi,
     &                                  g_d1_ov_Exch,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_noJK,alo,ahi,
     &                                  g_d1_ov_noJK,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_1e,alo,ahi,
     &                                  g_d1_ov_1e,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_eSji,alo,ahi,
     &                                  g_d1_ov_eSji,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,32) ispin,xyz
  32    format('------ g_d1-1(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,33) ispin,xyz
  33    format('------ g_d1-1(',i3,',',i3,')----------- END')
       endif

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_oo,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov,clo,chi)

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_Coul,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_Coul,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_Exch,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_Exch,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_noJK,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_noJK,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_1e,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_1e,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_eSji,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_eSji,clo,chi)
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,34) ispin,xyz
  34    format('------ g_d1-2(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,35) ispin,xyz
  35    format('------ g_d1-2(',i3,',',i3,')----------- END')
       endif

       enddo ! end-loop-xyz

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,1) ispin
 1       format('------g_d1(',i3,')---------- START')
         call ga_print(g_d1_oo)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- START'
         call ga_print(g_d1_ov)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- END'
        if (ga_nodeid().eq.0) 
     &   write(*,2) ispin
 2       format('------g_d1(',i3,')---------- END')     
       endif

       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
       call nga_add_patch(1.0d0       ,g_p1_oo,elo,ehi,
     &                    coeff(ispin),g_d1_oo,elo,ehi,
     &                                 g_p1_oo,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov,elo,ehi,
     &                    coeff(ispin),g_d1_ov,elo,ehi,
     &                                 g_p1_ov,elo,ehi)

       call nga_add_patch(1.0d0       ,g_p1_ov_Coul,elo,ehi,
     &                    coeff(ispin),g_d1_ov_Coul,elo,ehi,
     &                                 g_p1_ov_Coul,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov_Exch,elo,ehi,
     &                    coeff(ispin),g_d1_ov_Exch,elo,ehi,
     &                                 g_p1_ov_Exch,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov_noJK,elo,ehi,
     &                    coeff(ispin),g_d1_ov_noJK,elo,ehi,
     &                                 g_p1_ov_noJK,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov_1e,elo,ehi,
     &                    coeff(ispin),g_d1_ov_1e,elo,ehi,
     &                                 g_p1_ov_1e,elo,ehi)
       call nga_add_patch(1.0d0       ,g_p1_ov_eSji,elo,ehi,
     &                    coeff(ispin),g_d1_ov_eSji,elo,ehi,
     &                                 g_p1_ov_eSji,elo,ehi)
       if (.not.ga_destroy(g_u_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oo',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_Coul)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovJ',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_Exch)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovK',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_noJK)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovnJK',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_1e)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov1e',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_eSji)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oveSji',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
      enddo ! end-loop-ispin

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_p1-1---------- START'
         call ga_print(g_p1_oo)
         call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-1---------- END'
        if (ga_nodeid().eq.0) 
     &   write(*,11) npol,nogiao
 11   format('(npol,nogiao)=(',i3,',',i3,')')
       endif

      if (npol.eq.1 .and. type_NMR.eq.1) then ! this happens ONLY for NMR-restricted calc.
        if (ga_nodeid().eq.0)
     &    write(*,*) 'enter-scaling by 2 g_p1'
        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1_oo,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_Coul,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_Exch,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_noJK,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_1e,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_eSji,elo,ehi,2.0d0)
      endif
        if (type_NMR.eq.1) then ! store only if shield calc.
         shldfile  =0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile'  ,mt_int,1,shldfile)   ! for NLMO analysis
         if (shldfile.eq.1) then
c -------- g_c1 --> file ----------START
         ndir1=3
         ndir =6
         ndata=2 !  =1 write FCSD,PSOSO,sdens =2 write g_c1
         call util_file_name(lbl_nlmoshield,.false.,.false.,
     &                       zorafilename)
         if (.not.dft_zoraShield_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munudia,g_munupar1
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munu_h01
     &              nlist, ! in: list of selected atoms
     &               nocc, ! in: used here
     &              ndata, ! in: =1 write dia,par1,sdens =2 write g_c1
     &          g_munudia, ! in: dummy variable not used
     &         g_munupar1, ! in: dummy variable not used
     &         g_munu_h01, ! in: dummy variable not used
     &               g_c1, ! in: used here
     &             g_dens))! in: dummy not used here
     &   call errquit('get_P10_JK: dft_zoraShieldNLMO_write failed',
     &                0,DISK_ERR)
c -------- g_c1 --> file ----------END
         if (.not.ga_destroy(g_c1)) call 
     &       errquit('get_P10_JK: ga_destroy failed g_c1',0,GA_ERR)    
         endif
        endif
        if (type_NMR.eq.3) then ! store only if g-shift calc.
         gshiftfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
         if (gshiftfile.eq.1) then
c -------- g_c1 --> file ----------START
         ndir1=3
         ndir =6
         ndata=2 !  =1 write FCSD,PSOSO,sdens =2 write g_c1
         call util_file_name(lbl_nlmogshift,.false.,.false.,
     &                       zorafilename)
         if (.not.dft_zoraGshift_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuEPRdia,g_munuEPRpar1
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munuEPRHpar
     &               nocc, ! in: used here
     &              ndata, ! in: =1 write dia,par1,sdens =2 write g_c1
     &       g_munuEPRdia, ! in: dummy not used here
     &      g_munuEPRpar1, ! in: dummy not used here
     &      g_munuEPRHpar, ! in: dummy not used here
     &               g_c1, ! in: used here
     &            g_sdens))! in: dummy not used here
     &   call errquit('get_P10_1: dft_zoraGshiftNLMO_write failed',
     &                0,DISK_ERR)
c -------- g_c1 --> file ----------END
         if (.not.ga_destroy(g_c1)) call 
     &       errquit('get_P10_1: ga_destroy failed g_c1',0,GA_ERR)    
         endif
        endif      
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- START'
        call ga_print(g_p1_oo)
        call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- END'
       if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-tot---------- START'
        call ga_add(1.0d0,g_p1_oo,1.0d0,g_p1_ov,g_p1)
        call ga_print(g_p1)
       if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-tot---------- END'
       endif
       if (.not.ga_destroy(g_p1)) call 
     &    errquit('get_p1: ga_destroy failed g_p1',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oo',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov',0,GA_ERR)    

       if (.not.ga_destroy(g_d1_ov_Coul)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovJ',0,GA_ERR)          
       if (.not.ga_destroy(g_d1_ov_Exch)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovK',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov_noJK)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovnJK',0,GA_ERR)    
       if (.not.ga_destroy(g_d1_ov_1e)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov1e',0,GA_ERR)     
       if (.not.ga_destroy(g_d1_ov_eSji)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oveSji',0,GA_ERR)      
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Debugging routines ----------------------------------------START
      subroutine get_prelim_fock_debug( 
     &                            g_d2, ! out: 
     &                           g_rhs, ! out: rhs expression
     &                          g_rhs0, ! out: to be used in get_d1()
     &                      g_rhs_eSji, ! out: -Sji^k \epsilon_i
     &                         vectors, !  in: MO  coeffs
     &                            eval, !  in: energy vals
     &                             pos, !  in: Nuclear positions (x,y,z)
     &                          natoms, !  in: nr. selected nuclei (atoms)
     &                           basis, !  in: basis handle
     &                             nbf, !  in: nr. basis functions
     &                             nmo, !  in: nr. MOs (occ+virt)
     &                            npol, !  in: nr. of polarizations
     &                            nocc, !  in: nr. occ     MOs
     &                           nvirt) !  in: nr. virtual MOs
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer npol  ! nr. of polarizations
      integer g_rhs, ! OUT-1
     &        g_d1,  ! tmp
     &        g_d2,  ! OUT-2
     *        g_rhs0 ! OUT-3
      integer g_rhs_eSji
      integer g_u,g_s10,g_s10_1
      integer ntot,nbf,nmo,ispin,iocc,
     &        nocc(npol),nvirt(npol)
      integer shift,disp,xyz
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3),
     &        clo(3), chi(3),
     &        dlo(3), dhi(3)
      integer natoms
      double precision pos(3*natoms)
      double precision eval(nbf*npol),toscl
      integer basis ! basis handle
      integer vectors(npol)
      logical oskel
      external int_giao_1ega,
     &         giao_aotomo
      integer debug_prelim
      integer ndens

      debug_prelim=0 ! =1 for debugging

      oskel = .false.
c ------------ creating ga arrays ---------- START
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nocc(ispin)
      enddo

      if(.not.ga_create(MT_DBL,ntot,3,'prelim_fock_debug g_rhs0',
     &                  -1,-1,g_rhs0))
     &   call errquit('get_prelim_fock: ga_create failed g_rhs0',
     &                 0,GA_ERR)
      call ga_zero(g_rhs0)
      clo(1) = 3*npol
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'prelim_fock_debug g_d1',
     &                    chi,g_d1)) 
     &  call errquit('gprelim_fock: nga_create failed g_d1',
     &                0,GA_ERR)
      call ga_zero(g_d1)

c     Get S10 in GA and transform to MO set (virt,occ)
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10_1)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10_1',
     &               0,GA_ERR)
       call ga_zero(g_s10_1)
       ahi(3) = 3*npol
       if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                    alo,g_s10)) 
     &  call errquit('gprelim_fock: nga_create failed g_s10',
     &               0,GA_ERR)
c ------------ creating ga arrays ---------- END
       call ga_zero(g_s10)
       call int_giao_1ega(basis,basis,g_s10_1,'s10',
     &                    pos,natoms,oskel)   
      if (debug_prelim.eq.1) then
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10_1 ----START'
        call ga_print(g_s10_1)
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10_1 ----END'
      endif

       blo(1) = 1
       bhi(1) = nbf ! nmo-fixing4lineardependency
       blo(2) = 1
       bhi(2) = nbf ! nmo-fixing4lineardependency
       blo(3) = 1
       bhi(3) = 3
      do ispin=1,npol  
       disp=3*(ispin-1) 
       alo(1) = 1
       ahi(1) = nbf ! nmo-fixing4lineardependency
       alo(2) = 1
       ahi(2) = nbf ! nmo-fixing4lineardependency
       alo(3) = disp+1
       ahi(3) = disp+3   
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin

      if (debug_prelim.eq.1) then
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10 ----START'
        call ga_print(g_s10)
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ---gprelim: ---- g_s10 ----END'
      endif
c Note.- Output g_s10 is (nmo,nmo) not (nbf,nbf) matrix
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,3,nbf)

      if (debug_prelim.eq.1) then
        do ispin=1,npol
          if (ga_nodeid().eq.0) 
     &     write(*,*) '-------- MO vect(',ispin,')----START'
         call ga_print(vectors(ispin))
          if (ga_nodeid().eq.0) 
     &     write(*,*) '-------- MO vect(',ispin,')----END'
        enddo
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ----gprelim: ---- MO:g_s10 ---- START'
       call ga_print(g_s10)
       if (ga_nodeid().eq.0) 
     &  write(*,*) ' ----gprelim: ---- MO:g_s10 ---- END'
      endif

      do ispin=1,npol ! ++++++  START-loop-ispin
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) =  3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &                     alo,g_u)) 
     &   call errquit('gprelim_fock: nga_create failed g_u',
     &                0,GA_ERR)
       call ga_zero(g_u)
       disp=3*(ispin-1)
c     COPY : g_s10 -> g_s10_1 ---------- START
       blo(1) = 1
       bhi(1) = nmo
       blo(2) = 1
       bhi(2) = nmo
       blo(3) = 1
       bhi(3) = 3
       alo(1) = 1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nmo
       alo(3) = disp+1
       ahi(3) = disp+3  
       call ga_zero(g_s10_1) ! FA-04-28-12-added
       call nga_copy_patch('n',g_s10,alo,ahi,
     &                       g_s10_1,blo,bhi) 
c     COPY : g_s10 -> g_s10_1 ---------- END    
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi    
c ----definitions for g_rhs -------- START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c --- definitions for g_rhs -------- END   
c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block g_s10 with - (minus) eigenvalues 
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(3) = 1
       ahi(3) = 3
       disp=nbf*(ispin-1)  ! FA-04-28-12-energy-fix

       if (debug_prelim.eq.1) then
         if (ga_nodeid().eq.0) 
     &    write(*,*) ' ---gprelim: ---- BEF-scl:g_s10_1 ----START'
         call ga_print(g_s10_1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) ' ---gprelim: ---- BEF-scl:g_s10_1 ----END'    
       endif

        do iocc = 1, nocc(ispin)
         alo(2) = iocc
         ahi(2) = iocc
         toscl=-eval(disp+iocc) 

          if (debug_prelim.eq.1) then
           if (ga_nodeid().eq.0) then
            write(*,1) ispin,iocc,toscl
 1          format('E2scl(',i4,',',i5,')=',f15.8) 
           endif
          endif

         call nga_scale_patch(g_s10_1,alo,ahi,toscl)          
        enddo ! end-loop-iocc
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: - AFT-scl:g_s10_1(',ispin,') --START'
        call ga_print(g_s10_1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: - AFT-scl:g_s10_1(',ispin,') --END'    
      endif

c     Copy to g_rhs 
c     alo(1) and ahi(1) the same as before
       alo(2) = 1
       ahi(2) = nocc(ispin)

      if (debug_prelim.eq.1.and.ga_nodeid().eq.0) then
         write(*,16) alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3)
 16      format('alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
      endif

       call nga_copy_patch('n',g_s10_1,alo,ahi,
     &                           g_rhs,blo,bhi)
c ----- output: g_rhs_eSji --- START
       call nga_copy_patch('n',g_s10_1,alo,ahi,
     &                      g_rhs_eSji,blo,bhi)
c ----- output: g_rhs_eSji --- END
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- g_rhs(',ispin,') ----START'
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: ---- g_rhs(',ispin,') ----END'    
      endif
c
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
       alo(1) = 1
       ahi(1) = nocc(ispin)
       call nga_scale_patch(g_s10_1,alo,ahi,-0.5d0)
       call nga_copy_patch('n',g_s10_1,alo,ahi,
     &                             g_u,alo,ahi)
c         write(*,13) ispin,xyz,nmo,nbf,
c     &              alo(1),ahi(1),alo(2),ahi(2),
c     &              alo(3),ahi(3)
c 13      format('prelim-x: (ispin,xyz,nmo,nbf)=(',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'alo-ahi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,')')
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u(',ispin,') ----START'
        call ga_print(g_u)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u(',ispin,') ----END'    
      endif
c
c     We also need the occupied-occupied contribution of g_u contributing
c     to the first order density matrix. As this block does not change 
c     during the CPHF we can calculate it once and subtract it from the
c     RHS. We will reuse g_s10 as scratch space.
        call ga_zero(g_s10_1)
        alo(1) = 1
        alo(2) = 1
        blo(1) = 1
        bhi(1) = nbf
        blo(2) = 1
        clo(2) = 1
        chi(2) = nbf
        clo(3) = 1
        chi(3) = nbf
        dlo(1) = 1
        dhi(1) = nbf
        dlo(2) = 1
        dhi(2) = nocc(ispin)
c     Create "perturbed density matrix" for closed-closed g_u block
       do xyz = 1,3 ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
        bhi(2) = nmo 

c         write(*,3) ispin,xyz,nmo,nbf,
c     &              blo(1),bhi(1),blo(2),bhi(2),
c     &              blo(3),bhi(3),
c     &              alo(1),ahi(1),alo(2),ahi(2),
c     &              alo(3),ahi(3),
c     &              dlo(1),dhi(1),dlo(2),dhi(2),
c     &              dlo(3),dhi(3)
c 3      format('prelim-1: (ispin,xyz,nmo,nbf)=(',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'blo-bhi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'alo-ahi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'dlo-dhi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,')')

        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                    vectors(ispin),blo,bhi,  
     &                               g_u,alo,ahi,
     &                           g_s10_1,dlo,dhi)  
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u-1(',ispin,',',xyz,') ----START'
        call ga_print(g_s10_1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_u-1(',ispin,',',xyz,') ----END'    
      endif
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
        bhi(2) = nocc(ispin)
c     Minus sign as we subtract it from the RHS as we do not include 
c     it in the LHS
        disp=3*(ispin-1)
        clo(1) = disp+xyz
        chi(1) = disp+xyz

c         write(*,4) ispin,xyz,nmo,nbf,
c     &              blo(1),bhi(1),blo(2),bhi(2),
c     &              blo(3),bhi(3),
c     &              alo(1),ahi(1),alo(2),ahi(2),
c     &              alo(3),ahi(3),
c     &              clo(1),chi(1),clo(2),chi(2),
c     &              clo(3),chi(3)
c 4      format('prelim-2: (ispin,xyz,nmo,nbf)=(',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'blo-bhi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'alo-ahi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,') ',
c     &          'clo-chi=(',i3,',',i3,',',
c     &          i3,',',i3,',',i3,',',i3,')')


        call nga_matmul_patch('n','t',-1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,
     &                               g_s10_1,alo,ahi,
     &                                  g_d1,clo,chi)  
      if (debug_prelim.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_d1(',ispin,',',xyz,') ----START'
        call ga_print(g_d1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) ' ---gprelim: -- g_d1(',ispin,',',xyz,') ----END'    
      endif
       enddo ! end-loop-xyz
c ------------ back-up g_u --> g_rhs0 ---- START
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       call nga_copy_patch('n',g_u   ,alo,ahi,
     &                         g_rhs0,blo,bhi) ! copy to g_rhs0
c ---- back-up g_u --> g_rhs0 ---- END
c ---- Remove scratch ga arrays ------
       if (.not.ga_destroy(g_u)) call 
     &     errquit('gprelim_fock: ga_destroy failed g_u',0,GA_ERR)
      enddo ! end-loop-ispin

c -------- Creating g_d2 ---------------START
      ndens=3*npol
      clo(1) = ndens*2
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'prelim_fock_debug g_d2',
     &                    chi,g_d2)) 
     &  call errquit('gprelim_fock: nga_create failed g_d2',
     &                0,GA_ERR)
       call ga_zero(g_d2)
       blo(1) = 1
       bhi(1) = ndens   
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 1
       bhi(3) = nbf
      do ispin=1,npol 
       disp=ndens*(ispin-1) 
       alo(1) = disp+1
       ahi(1) = disp+ndens   
       alo(2) = 1
       ahi(2) = nbf
       alo(3) = 1
       ahi(3) = nbf
       call nga_copy_patch('n',g_d1,blo,bhi,
     &                         g_d2,alo,ahi) 
      enddo ! end-loop-ispin

 145  continue

      if (debug_prelim.eq.1) then
       if (ga_nodeid().eq.0) then
        write(*,*) '----------- g_d2------------ START'
       endif
       call ga_print(g_d1)
       call ga_print(g_d2)
       if (ga_nodeid().eq.0) then
        write(*,*) '----------- g_d2------------ END'
       endif
      endif
c -------- Creating g_dens1 ---------------END
       if (.not.ga_destroy(g_s10)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
       if (.not.ga_destroy(g_s10_1)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10_1',0,GA_ERR)
       if (.not.ga_destroy(g_d1)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_d1',0,GA_ERR)
      return
      end

      subroutine add_H10_debug( 
     &                   g_rhs, ! out: accumulated rhs expression
     &                g_rhs_1e, ! out : Fji^{k,1e}
     &                   ga_Fji, !  in: Fock 1st-deriv without V (pot.) contrib.
     &                  vectors, !  in: MO  coeffs
     &                      pos, !  in: Nuclear positions (x,y,z)
     &                   natoms, !  in: nr. selected nuclei (atoms)
     &                    basis, !  in: basis handle
     &                      nbf, !  in: nr. basis functions
     &                      nmo, !  in: nr. MOs (occ+virt)
     &                     npol, !  in: nr. of polarizations
     &                     nocc, !  in: nr. occ     MOs
     &                    nvirt, !  in: nr. virtual MOs
     &                  do_zora, !  in: .false. if doing nonrel calc
     &                     rtdb) !  in: for COSMO
                                 !      therefore: HF or DFT (restricted/unrestricted)
                                 !      Then evaluate analytical GIAO AOs. instead of numerical
c =============  Computing  ==============
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c ========================================
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      integer rtdb
      integer npol    ! nr. of polarizations
      logical do_zora ! logical to check if NOT doing zora calc
      integer nbf,nmo,ispin,
     &        nocc(npol),nvirt(npol)
      integer shift,disp
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3)
      integer vectors(npol)
      integer g_rhs,ga_Fji
      integer g_rhs_1e
      integer g_s10_1,g_s10 ! scratch ga arrays
      integer natoms
      double precision pos(3*natoms)
      integer basis ! basis handle
      logical oskel
      integer NoKinetic
      common /skipKinetic/NoKinetic ! goes to int_giao_1ega()
      external int_giao_1ega,
     &         giao_aotomo
c     bq charges
      integer nbq,nextbq,ncosbq,
     &        l_xyz,k_xyz ! dummy variables
      nbq = 0    ! COSMO-variables
      nextbq = 0 ! COSMO-variables
      ncosbq = 0 ! COSMO-variables
      oskel = .false.
c ----------- Create scratch ga-arrays ------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                   alo,g_s10_1)) 
     &    call errquit('add_H10: nga_create failed g_s10',
     &                  0,GA_ERR)
       call ga_zero(g_s10_1)
      ahi(3) = 3*npol
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',
     &                   alo,g_s10)) 
     &    call errquit('add_H10: nga_create failed g_s10',
     &                  0,GA_ERR)
       call ga_zero(g_s10)

      if (.not.(do_zora)) then ! do: HF or DFT (norel calc.)
c        call int_giao_1ega(basis,basis,g_s10_1,
c     &                     'l10' ,pos(1),natoms,oskel)
        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'l10' ,pos,natoms,oskel)
        NoKinetic=0 ! =0 DO-kinetic, =1 SKIP-kinetic
c        call int_giao_1ega(basis,basis,g_s10_1,
c     &                     'tv10',pos(1),natoms,oskel)
        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'tv10',pos,natoms,oskel)
      else                            ! do: zora (relativistic calc.)
c ----------- Create scratch ga-arrays ------- END
        NoKinetic=1 ! =0 DO-kinetic, =1 SKIP-kinetic
c        call int_giao_1ega(basis,basis,g_s10_1,
c     &                     'tv10',pos(1),natoms,oskel)
        call int_giao_1ega(basis,basis,g_s10_1,
     &                     'tv10',pos,natoms,oskel)
        call ga_add(1.0d0,g_s10_1,
     &              1.0d0,ga_Fji, ! update g_s10_1 with ga_Fji
     &                    g_s10_1)! out
      endif
c
c     Get external and cosmo bq contribution
      nbq = 0
      nextbq = 0
      ncosbq = 0
      if(geom_extbq_on()) nextbq = geom_extbq_ncenter()
      nbq = nextbq ! external bq's
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &    nbq = ncosbq ! cosmo bq's
      if (nextbq.gt.0.and.ncosbq.gt.0)
     &    nbq = nextbq + ncosbq  ! tally up cosmo and external bqs
      if (nbq.gt.0) then
        call int_giao_1ega(basis,basis,g_s10_1,'bq10',dbl_mb(k_xyz),
     &                     natoms,oskel)
      end if
c
c     Transform H10 to MO and add to g_rhs
c     Copy: g_s10_1 --> g_s10
       blo(1) = 1
       bhi(1) = nbf ! nmo-fixing4lineardependency
       blo(2) = 1
       bhi(2) = nbf ! nmo-fixing4lineardependency
       blo(3) = 1
       bhi(3) = 3
      do ispin=1,npol   
       shift=3*(ispin-1)
       alo(1) = 1
       ahi(1) = nbf ! nmo-fixing4lineardependency
       alo(2) = 1
       ahi(2) = nbf ! nmo-fixing4lineardependency
       alo(3) = shift+1
       ahi(3) = shift+3   
       call nga_copy_patch('n',g_s10_1,blo,bhi,
     &                           g_s10,alo,ahi) 
      enddo ! end-loop-ispin
      call giao_aotomo(g_s10,vectors,nocc,nvirt,npol,3,nbf)

      do ispin=1,npol
       shift=3*(ispin-1)
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = shift+1
       ahi(3) = shift+3   
c ------ definitions for g_rhs -------- START
       disp=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = disp+1
       bhi(1) = disp+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
c ------ definitions for g_rhs -------- END
       call nga_add_patch(1.0d0,g_rhs,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs,blo,bhi)
c ----- out: g_rhs_1e ------- START
       call nga_add_patch(1.0d0,g_rhs_1e,blo,bhi,
     &                    1.0d0,g_s10,alo,ahi,
     &                          g_rhs_1e,blo,bhi)
c ----- out: g_rhs_1e ------- END
      enddo ! end-loop-ispin
      if (.not.ga_destroy(g_s10_1)) call 
     &    errquit('add_H10: ga_destroy failed g_s10',0,GA_ERR)
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('add_H10: ga_destroy failed g_s10',0,GA_ERR)
      return
      end

      subroutine get_gshift_tensor_AorB(
     &                  rtdb,geom,basis,
     &                  g_rhs,
     &                  g_rhs0,
     &                  g_rhs_Coul,
     &                  g_rhs_Exch,
     &                  g_rhs_noJK,
     &                  g_rhs_1e,
     &                  g_rhs_eSji,
     &                  ga_h01_epr,
     &                  ga_para1_epr_AB,
     &                  ga_dia_epr_AB,
     &                  vectors,
     &                  g_CiFull,
     &                  coeffpol,
     &                  nbf,nmo,
     &                  npol,nocc,nvirt,
     &                  sh_atom, ! in : nr. of atoms
     &                  do_zora,do_NonRel,
     &                  not_zora_scale,
     &                  switch_gshift_analysis)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer nclosed(2), nocc(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, ic, iatom, iocc, ifld, ioff
      integer sh_atom1
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer l_dia, k_dia, l_para, k_para
      integer l_xyz, k_xyz, l_zan, k_zan
      integer l_tmp, k_tmp
      integer g_dens(3),g_h01,g_d1,
     &        g_rhs,g_rhs0,g_fock
      integer g_fock_Coul,g_fock_Exch,
     &        g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_eSji,g_rhs_1e
      integer g_d1_oo,g_d1_ov,type_NMR
      integer g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji,g_d1_ov_ExchSFB
      logical do_zora,do_NonRel,not_zora_scale,
     &        switch_gshift_analysis

      integer vectors(2), geomnew, i, j, ij
      double precision atn, tol2e, val, isotr, aniso
      double precision val_oo,val_ov
      double precision a(6),axs(3,3),eig(3),xfac
      character*255 zorafilename
      integer type_nmrdata,g_AtNr1,icalczora
      integer ga_dia_epr_AB,ga_para1_epr_AB,ga_Fji,ga_h01_epr
      integer g_CiFull(2)
      integer dft_zoraNMR_read
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      integer ld(2),cbuf,ind
      logical  cphf2, file_write_ga, file_read_ga, cphf
      double precision coeff_exch,ac_occ(2)
      
      logical status
      double precision ppt,coeffpol
      double precision par_arr(9,3,3) 
      double precision val_ov_Coul,val_ov_Exch,val_ov_noJK,
     &                 val_ov_1e,val_ov_eSji

      data ppt     /1.0d+03/  
      data tol2e   /1.0d-10/
      integer npol,ntot,ispin,nind_jk,cdens
      external cphf2, file_write_ga, file_read_ga, cphf
      external giao_aotomo
      external get_P10,get_P10_1,get_P10_JK
     &         dft_zoraNMR_read,get_P10_1_AorB,
     &         get_P10_JK_AorB
      integer debug_gshift
      integer shift,disp

      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh para',l_para,k_para)) 
     &    call errquit('hnd_gshift: ma_push_get failed k_para',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_gshift: ma_push_get failed k_dia',0,MA_ERR)

      do ispin=1,npol ! loop in spin----- START
       type_NMR=3 ! =1,2,3=shieldings,hyperfine,gshift
       if (switch_gshift_analysis) then ! --- START-if-switch_gshift_analysis
        call get_P10_JK_AorB(
     &               g_d1_oo,      ! out: Perturbed density matrix occ-occ  contrib
     &               g_d1_ov,      ! out: Perturbed density matrix occ-virt contrib
     &               type_NMR,     ! in : =1,2,3=shieldings,hyperfine,gshift
     &               g_d1_ov_Coul, ! out
     &               g_d1_ov_Exch, ! out
     &               g_d1_ov_noJK, ! out
     &               g_d1_ov_1e,   ! out
     &               g_d1_ov_eSji, ! out
     &                 g_rhs,      ! in: accumulated rhs expression
     &                 g_rhs_Coul, ! in
     &                 g_rhs_Exch, ! in
     &                 g_rhs_noJK, ! in
     &                   g_rhs_1e, ! in
     &                 g_rhs_eSji, ! in
     &                g_rhs0,      ! in: from get_prelim_fock()
     &               vectors,g_CiFull, 
     &               nbf,nmo,
     &               ispin,
     &               npol,nocc,nvirt,
     &               do_zora,do_NonRel,not_zora_scale,rtdb) 
      else
       call get_P10_1_AorB( 
     &           g_d1_oo, ! out: Perturbed density matrix occ-occ  contrib
     &           g_d1_ov, ! out: Perturbed density matrix occ-virt contrib
     &           type_NMR,! in : =1,2,3=shieldings,hyperfine,gshift
     &             g_rhs, ! in: accumulated rhs expression
     &            g_rhs0, ! in: from get_prelim_fock()
     &            vectors,g_CiFull, 
     &            nbf,nmo,
     &            ispin,
     &            npol,nocc,nvirt,
     &            do_zora,do_NonRel,not_zora_scale,rtdb) 
      endif ! ------------- END-if-switch_gshift_analysis
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the shielding tensor
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
c     s(para)xy = Sum(n,l) D1x(n,l) * H01y(n,l)
      do ixy = 1, 9*sh_atom
       dbl_mb(k_para+ixy-1) = 0.0d0  ! initialize the paramagnetic part
      enddo
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++ Initialize dbl_mb(k_para) with ga_para1 ++++ START
c ---- STORE: g_para1_epr --> dbl_mb(k_para)
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=ispin
      ahi(3)=ispin
       ld(1)=3
       ld(2)=3
      call nga_get(ga_para1_epr_AB,alo,ahi,dbl_mb(k_para),ld)
c +++++ Initialize dbl_mb(k_para) with g_para1 ++++ END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      ixy = 0
      blo(3) = 0
      bhi(3) = 0
      do iy = 1, 3
       blo(3) = blo(3) + 1
       bhi(3) = bhi(3) + 1
       do ix = 1, 3
        alo(3) = ix
        ahi(3) = ix
        ixy = ixy + 1
        if (switch_gshift_analysis) then ! --- START-if-switch_gshift_analysis
         val_oo = nga_ddot_patch(g_d1_oo,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov = nga_ddot_patch(g_d1_ov,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_Coul = nga_ddot_patch(g_d1_ov_Coul,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_Exch = nga_ddot_patch(g_d1_ov_Exch,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov_noJK = nga_ddot_patch(g_d1_ov_noJK,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)  
         val_ov_1e = nga_ddot_patch(g_d1_ov_1e,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)   
         val_ov_eSji = nga_ddot_patch(g_d1_ov_eSji,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)             
        else
         val_oo = nga_ddot_patch(g_d1_oo,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
         val_ov = nga_ddot_patch(g_d1_ov,'n',alo,ahi,
     &                          ga_h01_epr,'n',blo,bhi)
        endif
        cbuf=3*(ix-1)+iy-1 ! transpose
c ----- store in par_arr ------- START
        par_arr(1,iy,ix)=dbl_mb(k_para+cbuf)
        par_arr(2,iy,ix)=val_oo*coeffpol*ppt*(-0.25d0)
        par_arr(3,iy,ix)=val_ov*coeffpol*ppt*(-0.25d0)
        par_arr(4,iy,ix)=dbl_mb(k_para+cbuf)+
     &                   val_oo*coeffpol*ppt*(-0.25d0)+
     &                   val_ov*coeffpol*ppt*(-0.25d0)
        if (switch_gshift_analysis) then ! --- START-if-switch_gshift_analysis
         par_arr(5,iy,ix)=val_ov_Coul*coeffpol*ppt*(-0.25d0)
         par_arr(6,iy,ix)=val_ov_Exch*coeffpol*ppt*(-0.25d0)
         par_arr(7,iy,ix)=val_ov_noJK*coeffpol*ppt*(-0.25d0)
         par_arr(8,iy,ix)=val_ov_1e*coeffpol*ppt*(-0.25d0)
         par_arr(9,iy,ix)=val_ov_eSji*coeffpol*ppt*(-0.25d0)
        endif
c ----- store in par_arr ------- END      
        dbl_mb(k_para+cbuf)=dbl_mb(k_para+cbuf)+
     &                      val_oo*coeffpol*ppt*(-0.25d0)+
     &                      val_ov*coeffpol*ppt*(-0.25d0)
       enddo
      enddo
c
c     s(dia)xy = Sum(n,l) D0(n,l) * H11(dia)xy(n,l) 
c     Ordering of H11 blocks is Bxx,Bxy,Bxz,Byx,etc 
      do ixy = 1,9*sh_atom
         dbl_mb(k_dia+ixy-1)  = 0.0d0  ! initialize the diamagnetic part
      enddo
c ---- STORE: ga_dia_epr --> dbl_mb(k_dia)
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=ispin ! assuming for g-shifts natoms=1
      ahi(3)=ispin ! assuming for g-shifts natoms=1
      ld(1)=3
      ld(2)=3 
      call nga_get(ga_dia_epr_AB,alo,ahi,dbl_mb(k_dia),ld)
c +++++++++ print-total-pardia-transferred +++ START
       if (ga_nodeid().eq.0) then
        ic=1
        do ix = 1, 3
         do iy = 1, 3
          if (switch_gshift_analysis) then ! --- START-if-switch_gshift_analysis
           write(*,179) ispin, 
     &                 ix,iy,dbl_mb(k_dia+ic-1),
     &                 par_arr(1,ix,iy),par_arr(2,ix,iy),
     &                 par_arr(3,ix,iy),
     &                 par_arr(5,ix,iy),par_arr(6,ix,iy),
     &                 par_arr(7,ix,iy),
     &                 par_arr(8,ix,iy),par_arr(9,ix,iy),
     &                 par_arr(4,ix,iy)
 179        format('NW-',i1,':',
     &             '(dia,gauge,OO,OV,',
     &             'OV_Coul,OV_Exch,OV_nJK,',
     &             'OV_1e,OV_eSji,'
     &             'Totpar)(',i1,',',i1,')=(',
     &             f12.6,' ',f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' ',f12.6,' ',
     &             f12.6,' )')
          else
           write(*,19) ispin,
     &                 ix,iy,dbl_mb(k_dia+ic-1),
     &                 par_arr(1,ix,iy),par_arr(2,ix,iy),
     &                 par_arr(3,ix,iy),par_arr(4,ix,iy)
 19        format('NW-',i1,':',
     &             '(dia,gauge,OO,OV,Totpar)(',i1,',',i1,')=(',
     &             f12.6,' ',f12.6,' ',f12.6,' ',f12.6,' ',
     &             f12.6,' )')
          endif
          ic=ic+1
         enddo ! end-loop-iy
        enddo ! end-loop-ix
       endif ! ---- END-if-ga_nodeid-eq-0
c ----- destroy created GA arrays in get_P10_JK ---- START
       if (.not.ga_destroy(g_d1_oo)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
       if (.not.ga_destroy(g_d1_ov)) call 
     &  errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
       if (switch_gshift_analysis) then ! --- START-if-switch_gshift_analysis
        if (.not.ga_destroy(g_d1_ov_Coul)) call 
     &   errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_Exch)) call 
     &   errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_noJK)) call 
     &   errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_1e)) call 
     &   errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
        if (.not.ga_destroy(g_d1_ov_eSji)) call 
     &   errquit('gprelim_fock: ga_destroy failed g_s10',0,GA_ERR)
       endif ! ---- END-if-switch_gshift_analysis
c ----- destroy created GA arrays in get_P10_JK ---- END
      enddo ! ! loop in spin----- END
c     Clean up all remaining memory
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_dia',0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_gshift: ma_pop_stack failed k_para',0,MA_ERR)
      return
      end
      subroutine get_P10_JK_AorB(
     &                     g_p1_oo, ! out: Perturbed (spin)density matrix-occ-occ       contrib
     &                     g_p1_ov, ! out: Perturbed (spin)density matrix-occ-virt      contrib
     &                    type_NMR, ! in: =1,2,3=shieldings,hyperfine,gshift
     &                g_p1_ov_Coul, ! out: Perturbed (spin)density matrix with Coulomb  contrib
     &                g_p1_ov_Exch, ! out: Perturbed (spin)density matrix with Exchange contrib (fock_xc)
     &                g_p1_ov_noJK, ! out: Perturbed (spin)density matrix remaining terms not J or K
     &                  g_p1_ov_1e, ! out: Perturbed (spin)density matrix 1e-operator contrib
     &                g_p1_ov_eSji, ! out: Perturbed (spin)density matrix perturbed overlap matrix
     &                       g_rhs, ! in: accumulated rhs expression
     &                  g_rhs_Coul,
     &                  g_rhs_Exch,
     &                  g_rhs_noJK,
     &                    g_rhs_1e,
     &                  g_rhs_eSji,
     &                      g_rhs0, ! in: from get_prelim_fock()         
     &                     vectors, ! in: MO                coeffs
     &                    g_CiFull, ! in: MO zora weighting coeffs
     &                         nbf, ! in: nr. basis functions
     &                         nmo, ! in: nr. MOs (occ+virt)
     &                       ispin, ! in: select polarization = 1(A) or 2(B)
     &                        npol, ! in: nr. of polarizations
     &                        nocc, ! in: nr. occ     MOs
     &                       nvirt, ! in: nr. virtual MOs
     &                     do_zora, ! in: .true.  if doing zora calc
     &                   do_NonRel, ! in: .true.  if doing nonrel within zora scheme
     &              not_zora_scale, ! in: .true.  not scaling perturbed density matrix
     &                        rtdb)
c  g_p1_ov      = g_p1_ov_Coul + g_p1_ov_Exch + g_p1_ov_noJK
c  g_p1_ov_noJK = g_p1_ov_1e   + g_p1_ov_eSji
c Note: If shieldings calc --> g_p1= Perturbed        density matrix
c                                  = P^(1,0)_A +  P^(1,0)_B
c       ==> This calc. could be npol=1 (  restricted shell calc.)
c                               npol=2 (unrestricted shell calc.)       
c       If gshift     calc --> g_p1= Perturbed (spin) density matrix
c                                  = P^(1,0)_A -  P^(1,0)_B
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
      
      integer g_p1 ! for debugging

      integer npol   ! nr. of polarizations
      integer g_p1_oo,   ! OUT: Perturbed density matrix occ-occ  contrib
     &        g_p1_ov    ! OUT: Perturbed density matrix occ-virt contrib 
      integer g_rhs  ! IN : accumulated right-hand-side expression
      integer g_rhs0 ! IN : from get_prelim_fock()
      integer g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_1e,g_rhs_eSji,
     &        g_p1_ov_Coul,g_p1_ov_Exch,g_p1_ov_noJK,
     &        g_p1_ov_1e,g_p1_ov_eSji,
     &        g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji,
     &        g_u_ov_Coul,g_u_ov_Exch,g_u_ov_noJK,
     &        g_u_ov_1e,g_u_ov_eSji
      logical do_zora,do_NonRel,not_zora_scale
      integer g_u_oo,g_u_ov,g_d1_oo,g_d1_ov ! scratch ga array
      integer vectors(npol),vectors_scl(npol),
     &        g_CiFull(npol)
      integer rtdb
      integer nbf,nmo,ispin,type_NMR,
     &        nocc(npol),nvirt(npol)
      integer shift,xyz,slc_spinpolAO
      integer alo(3), ahi(3), 
     &        blo(3), bhi(3), 
     &        clo(3), chi(3),
     &        dlo(3), dhi(3),
     &        elo(3), ehi(3)
      logical status
      integer nogshift
      data nogshift/1/            
      integer nogiao             
      data nogiao/1/             
      integer debug_p10

      if(.not.rtdb_get(rtdb,'zora:slc_spinpolAO',        ! FA
     &                 mt_int,1,slc_spinpolAO))          ! FA
     &  slc_spinpolAO= 0 ! 0=A-B,1=A,2=-B selecting spin-density matrix      

       debug_p10=0 ! =1 for debugging
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_d1_oo',
     &                     alo,g_d1_oo)) call 
     &     errquit('g_d1_oo: nga_create failed g_d1_oo',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_p1_oo',
     &                     alo,g_p1_oo)) call 
     &     errquit('g_p1_oo: nga_create failed g_p1_oo',0,GA_ERR)
       call ga_zero(g_p1_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_d1_ov',
     &                     alo,g_d1_ov)) call 
     &     errquit('g_d1_ov: nga_create failed g_d1_ov',0,GA_ERR)

       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_d1_ov_Coul',
     &     alo,g_d1_ov_Coul)) call 
     &     errquit('g_d1_ovJ: nga_create failed g_d1_ovJ',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_d1_ovJ',
     &     alo,g_d1_ov_Exch)) call 
     &     errquit('g_d1_ovK: nga_create failed g_d1_ovK',0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_d1_ovnoJK',
     &     alo,g_d1_ov_noJK)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_ovnJK',
     &     0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_d1_ov_1e)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_ov1e',
     &     0,GA_ERR)
       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_d1_ov_eSji)) call 
     &     errquit('g_d1_ovnJK: nga_create failed g_d1_oveSji',
     &     0,GA_ERR)

       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_p1_ov)) call 
     &     errquit('g_p1_ov: nga_create failed g_p1_ov',0,GA_ERR)
       call ga_zero(g_p1_ov)

       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_p1_ov_Coul)) call 
     &     errquit('g_p1_ovJ: nga_create failed g_p1_ovJ',0,GA_ERR)
       call ga_zero(g_p1_ov_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_p1_ov_Exch)) call 
     &     errquit('g_p1_ovK: nga_create failed g_p1_ovK',0,GA_ERR)
       call ga_zero(g_p1_ov_Exch)
       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_p1_ov_noJK)) call 
     &     errquit('g_p1_ovnJK: nga_create failed g_p1_ovnJK',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_noJK)
       if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',
     &     alo,g_p1_ov_1e)) call 
     &     errquit('g_p1_ov1e: nga_create failed g_p1_ov1e',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_1e)
       if (.not.nga_create(MT_DBL,3,ahi,'P10_JK_AorB g_p1_ov_eSji',
     &     alo,g_p1_ov_eSji)) call 
     &     errquit('g_p1_oveSji: nga_create failed g_p1_oveSji',
     &     0,GA_ERR)
       call ga_zero(g_p1_ov_eSji)
c ------ For debugging ---- START
c       if (.not.nga_create(MT_DBL,3,ahi,
c     &          'P10_JK_AorB g_p1',alo,g_p1)) call 
c     &     errquit('g_p1: nga_create failed g_p1',0,GA_ERR)
c       call ga_zero(g_p1)
c ------ For debugging ---- END

c      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_oo)) call 
     &    errquit('g_u_oo: nga_create failed g_u_oo',0,GA_ERR)
       call ga_zero(g_u_oo)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u_ov)) call 
     &    errquit('g_u_ov: nga_create failed g_u_ov',0,GA_ERR)
       call ga_zero(g_u_ov)

       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_Coul)) call 
     &    errquit('g_u_ovJ: nga_create failed g_u_ovJ',0,GA_ERR)
       call ga_zero(g_u_ov_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_Exch)) call 
     &    errquit('g_u_ovK: nga_create failed g_u_ovK',0,GA_ERR)
       call ga_zero(g_u_ov_Exch)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_noJK)) call 
     &    errquit('g_u_ovnJK: nga_create failed g_u_ovnJK',0,GA_ERR)
       call ga_zero(g_u_ov_noJK)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_1e)) call 
     &    errquit('g_u_ov1e: nga_create failed g_u_ov1e',0,GA_ERR)
       call ga_zero(g_u_ov_1e)
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',
     &    alo,g_u_ov_eSji)) call 
     &    errquit('g_u_oveSji: nga_create failed g_u_oveSji',
     &    0,GA_ERR)
       call ga_zero(g_u_ov_eSji)
c ----------- copy occ-occ info  ------------ START
       shift=nocc(1)*nocc(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nocc(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = 1
       ahi(1) = nocc(ispin)
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3
       call nga_copy_patch('n',g_rhs0,blo,bhi,
     &                            g_u_oo,alo,ahi) 
c ----------- copy occ-occ info  ------------ END
c ----------- copy occ-virt info ------------ START
       shift=nocc(1)*nvirt(1)*(ispin-1)
       blo(1) = shift+1
       bhi(1) = shift+nocc(ispin)*nvirt(ispin)
       blo(2) = 1
       bhi(2) = 3
       alo(1) = nocc(ispin)+1
       ahi(1) = nmo
       alo(2) = 1
       ahi(2) = nocc(ispin)
       alo(3) = 1
       ahi(3) = 3

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,16) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3)
 16      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ')
        endif
       endif
       call nga_copy_patch('n',g_rhs,blo,bhi,
     &                           g_u_ov,alo,ahi)  

       call nga_copy_patch('n',g_rhs_Coul,blo,bhi,
     &                           g_u_ov_Coul,alo,ahi)  
       call nga_copy_patch('n',g_rhs_Exch,blo,bhi,
     &                           g_u_ov_Exch,alo,ahi)  
       call nga_copy_patch('n',g_rhs_noJK,blo,bhi,
     &                           g_u_ov_noJK,alo,ahi)  
       call nga_copy_patch('n',g_rhs_1e,blo,bhi,
     &                           g_u_ov_1e,alo,ahi)  
       call nga_copy_patch('n',g_rhs_eSji,blo,bhi,
     &                           g_u_ov_eSji,alo,ahi)  
c ----------- copy occ-virt info ------------ START
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '------- g_u_oo------- START'
        call ga_print(g_u_oo)
         if (ga_nodeid().eq.0) then
          write(*,*) '------- g_u_oo------- END'
          write(*,*) '------- g_u_ov------- START'
         endif
        call ga_print(g_u_ov)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------- g_u_ov------- END'
       endif

c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
       alo(1) = 1
       alo(2) = 1
       blo(1) = 1
       blo(2) = 1
       clo(1) = 1
       chi(1) = nbf
       clo(2) = 1
       chi(2) = nbf
       dlo(1) = 1
       dlo(2) = 1
       dhi(1) = nbf
       dhi(2) = nocc(ispin)
c --------- zora scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(ispin),
     &                        vectors_scl(ispin),'vscl 1'))
     &  call errquit('g_d1: ga_duplicate failed',1,GA_ERR)
       call ga_copy(vectors(ispin),vectors_scl(ispin))
       if (do_zora .and. .not.(do_NonRel) .and.
     &     .not.(not_zora_scale)) then
        if (ga_nodeid().eq.0)
     &   write(*,*) 'get_P10_JK_AorB:: scaling perturbed density matrix'
        call ga_scale_cols(vectors_scl(ispin),g_CiFull(ispin))
       endif
c --------- zora scaling of MO vectors(1) ----- END
       do xyz = 1,3  ! = x,y,z
        alo(3) = xyz
        ahi(3) = xyz
        blo(3) = xyz
        bhi(3) = xyz
        clo(3) = xyz
        chi(3) = xyz
        dlo(3) = xyz
        dhi(3) = xyz
        bhi(1) = nbf
        bhi(2) = nmo 
        ahi(1) = nmo
        ahi(2) = nocc(ispin)
c     Make C1       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,*) 'bef-1-matmul-patch'
         write(*,17) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               dlo(1),dhi(1),dlo(2),dhi(2),
     &               dlo(3),dhi(3)
 17      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'dlo-dhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
         endif
        endif
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,dlo,dhi) 
        call nga_copy_patch('n',g_d1_oo,dlo,dhi,
     &                           g_u_oo,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov,dlo,dhi,
     &                           g_u_ov,dlo,dhi)

        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_Coul,alo,ahi,
     &                                  g_d1_ov_Coul,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_Coul,dlo,dhi,
     &                           g_u_ov_Coul,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_Exch,alo,ahi,
     &                                  g_d1_ov_Exch,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_Exch,dlo,dhi,
     &                           g_u_ov_Exch,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_noJK,alo,ahi,
     &                                  g_d1_ov_noJK,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_noJK,dlo,dhi,
     &                           g_u_ov_noJK,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_1e,alo,ahi,
     &                                  g_d1_ov_1e,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_1e,dlo,dhi,
     &                           g_u_ov_1e,dlo,dhi)
        call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                        vectors(ispin),blo,bhi,  
     &                                   g_u_ov_eSji,alo,ahi,
     &                                  g_d1_ov_eSji,dlo,dhi) 
        call nga_copy_patch('n',g_d1_ov_eSji,dlo,dhi,
     &                           g_u_ov_eSji,dlo,dhi)
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0)
     &     write(*,30) ispin,xyz
  30     format('------ g_u(',i3,',',i3,')----------- START')
         call ga_print(g_u_oo)
         call ga_print(g_u_ov)
         if (ga_nodeid().eq.0)
     &    write(*,31) ispin,xyz
  31     format('------ g_u(',i3,',',i3,')----------- END')
        endif

        bhi(1) = nbf
        bhi(2) = nocc(ispin)
        ahi(1) = nocc(ispin)
        ahi(2) = nbf
c     Make D1
       
       if (debug_p10 .eq.1) then
         if (ga_nodeid().eq.0) then
         write(*,18) blo(1),bhi(1),blo(2),bhi(2),
     &               blo(3),bhi(3),
     &               alo(1),ahi(1),alo(2),ahi(2),
     &               alo(3),ahi(3),
     &               clo(1),chi(1),clo(2),chi(2),
     &               clo(3),chi(3)
 18      format('blo-bhi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'alo-ahi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,') ',
     &          'clo-chi=(',i3,',',i3,',',
     &          i3,',',i3,',',i3,',',i3,')')
        endif
       endif

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_oo,alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov,alo,ahi,
     &                                  g_d1_ov,clo,chi)

        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_Coul,alo,ahi,
     &                                  g_d1_ov_Coul,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_Exch,alo,ahi,
     &                                  g_d1_ov_Exch,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_noJK,alo,ahi,
     &                                  g_d1_ov_noJK,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_1e,alo,ahi,
     &                                  g_d1_ov_1e,clo,chi)
        call nga_matmul_patch('n','t',1.0d0,0.0d0,
     &                    vectors_scl(ispin),blo,bhi,  
     &                                   g_u_ov_eSji,alo,ahi,
     &                                  g_d1_ov_eSji,clo,chi)

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,32) ispin,xyz
  32    format('------ g_d1-1(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,33) ispin,xyz
  33    format('------ g_d1-1(',i3,',',i3,')----------- END')
       endif

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_oo,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_oo,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov,clo,chi)

        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_Coul,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_Coul,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_Exch,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_Exch,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_noJK,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_noJK,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_1e,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_1e,clo,chi)
        call nga_matmul_patch('n','t',-1.0d0,1.0d0,
     &                                   g_u_ov_eSji,blo,bhi,  
     &                    vectors_scl(ispin),alo,ahi,
     &                                  g_d1_ov_eSji,clo,chi)
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,34) ispin,xyz
  34    format('------ g_d1-2(',i3,',',i3,')----------- START')
        call ga_print(g_d1_oo)
        call ga_print(g_d1_ov)
        if (ga_nodeid().eq.0)
     &   write(*,35) ispin,xyz
  35    format('------ g_d1-2(',i3,',',i3,')----------- END')
       endif

       enddo ! end-loop-xyz

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,1) ispin
 1       format('------g_d1(',i3,')---------- START')
         call ga_print(g_d1_oo)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- START'
         call ga_print(g_d1_ov)
         if (ga_nodeid().eq.0)
     &    write(*,*) '---g_d1_ov ---- END'
        if (ga_nodeid().eq.0) 
     &   write(*,2) ispin
 2       format('------g_d1(',i3,')---------- END')     
       endif

       elo(1) = 1
       ehi(1) = nbf
       elo(2) = 1
       ehi(2) = nbf
       elo(3) = 1
       ehi(3) = 3
       call nga_add_patch(1.0d0,g_p1_oo,elo,ehi,
     &                    1.0d0,g_d1_oo,elo,ehi,
     &                          g_p1_oo,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov,elo,ehi,
     &                    1.0d0,g_d1_ov,elo,ehi,
     &                          g_p1_ov,elo,ehi)

       call nga_add_patch(1.0d0,g_p1_ov_Coul,elo,ehi,
     &                    1.0d0,g_d1_ov_Coul,elo,ehi,
     &                          g_p1_ov_Coul,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov_Exch,elo,ehi,
     &                    1.0d0,g_d1_ov_Exch,elo,ehi,
     &                          g_p1_ov_Exch,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov_noJK,elo,ehi,
     &                    1.0d0,g_d1_ov_noJK,elo,ehi,
     &                          g_p1_ov_noJK,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov_1e,elo,ehi,
     &                    1.0d0,g_d1_ov_1e,elo,ehi,
     &                          g_p1_ov_1e,elo,ehi)
       call nga_add_patch(1.0d0,g_p1_ov_eSji,elo,ehi,
     &                    1.0d0,g_d1_ov_eSji,elo,ehi,
     &                          g_p1_ov_eSji,elo,ehi)
       if (.not.ga_destroy(g_u_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oo',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_Coul)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovJ',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_Exch)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovK',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_noJK)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ovnJK',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_1e)) call 
     &    errquit('get_d1: ga_destroy failed g_u_ov1e',0,GA_ERR)
       if (.not.ga_destroy(g_u_ov_eSji)) call 
     &    errquit('get_d1: ga_destroy failed g_u_oveSji',0,GA_ERR)
       if (.not.ga_destroy(vectors_scl(ispin))) call 
     &    errquit('get_d1: ga_destroy failed vscl',0,GA_ERR)
c      enddo ! end-loop-ispin

       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_p1-1---------- START'
         call ga_print(g_p1_oo)
         call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-1---------- END'
        if (ga_nodeid().eq.0) 
     &   write(*,11) npol,nogiao
 11   format('(npol,nogiao)=(',i3,',',i3,')')
       endif

      if (npol.eq.1 .and. nogiao.eq.0) then ! this happens ONLY for NMR-restricted calc.
        if (ga_nodeid().eq.0)
     &    write(*,*) 'enter-scaling by 2 g_p1'

        elo(1) = 1
        ehi(1) = nbf
        elo(2) = 1
        ehi(2) = nbf
        elo(3) = 1
        ehi(3) = 3
        call nga_scale_patch(g_p1_oo,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_Coul,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_Exch,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_noJK,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_1e,elo,ehi,2.0d0)
        call nga_scale_patch(g_p1_ov_eSji,elo,ehi,2.0d0)
      endif
      
       if (debug_p10 .eq.1) then
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- START'
        call ga_print(g_p1_oo)
        call ga_print(g_p1_ov)
        if (ga_nodeid().eq.0) 
     &    write(*,*) '------g_p1-2---------- END'
c       if (ga_nodeid().eq.0) 
c     &    write(*,*) '------g_p1-tot---------- START'
c        call ga_add(1.0d0,g_p1_oo,1.0d0,g_p1_ov,g_p1)
c        call ga_print(g_p1)
c       if (ga_nodeid().eq.0) 
c     &    write(*,*) '------g_p1-tot---------- END'
       endif
c       if (.not.ga_destroy(g_p1)) call 
c     &    errquit('get_p1: ga_destroy failed g_d1_oo',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_oo)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oo',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov',0,GA_ERR)    

       if (.not.ga_destroy(g_d1_ov_Coul)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovJ',0,GA_ERR)          
       if (.not.ga_destroy(g_d1_ov_Exch)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovK',0,GA_ERR)   
       if (.not.ga_destroy(g_d1_ov_noJK)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ovnJK',0,GA_ERR)    
       if (.not.ga_destroy(g_d1_ov_1e)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_ov1e',0,GA_ERR)     
       if (.not.ga_destroy(g_d1_ov_eSji)) call 
     &    errquit('get_d1: ga_destroy failed g_d1_oveSji',0,GA_ERR)      
      return
      end
c Debugging routines ----------------------------------------START
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c -------------- for g-shift NMLO analysis ----------------- START
      subroutine create_munu4nbo_gshift(
     &             rtdb,    ! in: rtdb handle
     &             g_tvec,  ! in: eigenvectors or T diagonalizing matrix
     &             basis,   ! in: basis handle
     &             npol,    ! in: nr. polarizaitons
     &             nocc,    ! in: nr. occ   nocc(i) i=1,npol
     &             nvirt,   ! in: nr. virt nvirt(i) i=1,npol
     &             nmo)     ! in: nr. MO

      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh" 
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "msgids.fh" 
#include "zora.fh"   
c FA: Revised on 06-22-11
c ------Main outputs -------- START
      integer g_munu_rot,         ! hyp-FCSD  dia  part
     &        g_munu_rot1,
     &        g_munu_rot2,
     &        g_acc2  ! hyp-PSOSO para part
c ------Main outputs -------- END
      integer rtdb,basis
      integer g_sdens,g_tvec
      integer g_munuEPRdia,g_munuEPRpar1,g_munuEPRHpar,
     &        g_tnp,g_acc,vectors(2),
     &        g_tnp1,g_acc1,g_acc3     
      integer ispin,iat1  
      double precision ac_val,val1,sign    
      integer i,j,k,m,n,ndir,ndir1
      integer jlo,jhi,s,nbf,nmo,nsize,nsize1,nsize2
      integer npol,ntot
      integer ind,nlst,count,nocc(npol),nvirt(npol)
      integer Ndir_munu
c     Ndir_munu, Nr. of directions stored
c                =3  xx yy zz
      double precision coeff,fact,para_rot(9)
      double precision tmn(2),chcdata(3)
      integer jlo1,jhi1,jlo2,jhi2
      integer g_sdens1,g_c1,
     &        g_p10,g_munuEPRHpar2d,g_munuEPRpar12d
      integer iind(2),jind(2),icalczora,type_NMR  
      integer alo(3),ahi(3),elo(3),ehi(3),flo(3),fhi(3)
      logical dft_zoraGshift_NLMOAnalysis_read ! for read-nlmo-mat
      character*255 zorafilename              ! for read-nlmo-mat
      integer arr_ind(6,2)
       data ((arr_ind(j,i),i=1,2),j=1,6)
     &  /1,1,2,2,3,3,1,2,1,3,2,3/
      external dft_zoraHYP_NLMOAnalysis_read,get_P10_rot,
     &         fill_munuPSOSO_1,get_par_gshift_rot,
     &         get2dmat
c     --> To store ONLY munu principal components xx,yy,zz 
c     g_munuEPRdia    is created in dft_zora_EPR.F
c     size(g_munuEPRdia)=nlst*ndir (linear vector)
c     g_sdens, spin density matrix
c     nbf x nbf elements (bidimensional matrix)
c     Legend:
c     ndir=6 = xx, yy, zz, xy, xz, yz
c     nbf, Nr of basis functions
c     nlst=nbf*(nbf+1)/2

      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('munu: bas_numbf failed',555, BASIS_ERR)
      Ndir_munu=3
      nlst=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
c ++++++ Read NLMO matrices +++++++++ START
      ndir=6
      ndir1=3
      call util_file_name(lbl_nlmogshift,.false.,.false.,
     &                    zorafilename)    
      icalczora = 0  ! initialize the flag
      if (.not.dft_zoraGshift_NLMOAnalysis_read(
     &       zorafilename, ! in : filename
     &                nbf, ! in : nr basis functions
     &               ndir, ! in : nr of directions: 6 = xx yy zz xy xz yz
     &              ndir1, ! in : nr of directions: 3 = x y z
     &               nocc, ! in : nocc(i) i=1,2 nr. occupations in alpha and beta   
     &               npol, ! in : nr polarizations
     &       g_munuEPRdia, ! out: munu matrix of dia
     &      g_munuEPRpar1, ! out: munu matrix of par 1st term
     &      g_munuEPRHpar, ! out: munu matrix of H10
     &            vectors, ! out: MOs
     &               g_c1, ! out: perturbed MO
     &            g_sdens)) icalczora=1 

      goto 10

      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuEPRdia -- START'
        call ga_print(g_munuEPRdia)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuEPRdia -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &   'AFT-reading-HYP-NLMO matrices: g_munuEPRHpar -- START'
        call ga_print(g_munuEPRHpar)
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &  'AFT-reading-HYP-NLMO matrices: g_munuEPRHpar -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_sdens ----- START'
        call ga_print(g_sdens)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_sdens ----- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- START'
        call ga_print(g_c1)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- END'
 10   continue
c ++++++ Read NLMO matrices +++++++++ END

      call get_unique_elmat(g_sdens,g_sdens1,nlst,nbf)   ! out: g_sdens1
      ndir=6 ! Nr. of directions: xx,yy,zz,xy,xz,yz
      nsize=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
      nsize1=nsize*ndir   ! size of whole munu per atom
      nsize2=nsize*ndir1  ! size of whole munu per atom
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_tnp',0,0,g_tnp)) 
     $    call errquit('munu4nbo: g_tnp', 0,GA_ERR)
        call ga_zero(g_tnp)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_tnp1',0,0,g_tnp1)) 
     $    call errquit('munu4nbo: g_tnp1', 0,GA_ERR)
        call ga_zero(g_tnp1)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_acc',0,0,g_acc)) 
     $    call errquit('munu4nbo: g_acc', 0,GA_ERR)
        call ga_zero(g_acc)

        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_acc3',0,0,g_acc3)) 
     $    call errquit('munu4nbo: g_acc3', 0,GA_ERR)
        call ga_zero(g_acc3)
        if (.not. ga_create(mt_dbl,1,nsize,
     &                      'munu4nbo: g_acc1',0,0,g_acc1)) 
     $    call errquit('munu4nbo: g_acc1', 0,GA_ERR)
        call ga_zero(g_acc1)       
         alo(1) = nbf
         alo(2) = -1
         alo(3) = -1
         ahi(1) = nbf
         ntot=nocc(1)+nocc(2)
         ahi(2) = ntot
         ahi(3) = 3
         if (.not.nga_create(MT_DBL,3,ahi,'g_acc2 matrix',
     &       alo,g_acc2)) call 
     &        errquit('g_acc2: nga_create failed g_acc2',0,GA_ERR)
         call ga_zero(g_acc2)
        if (.not. ga_create(mt_dbl,1,nlst*3,
     &                       'munu4nbo: g_munu_rot',0,0,g_munu_rot)) 
     $    call errquit('munu4nbo: g_munu_rot', 0,GA_ERR)
        if (.not. ga_create(mt_dbl,1,nlst*3,
     &                       'munu4nbo: g_munu_rot2',0,0,g_munu_rot2)) 
     $    call errquit('munu4nbo: g_munu_rot2', 0,GA_ERR)
        if (.not. ga_create(mt_dbl,1,nlst*3,
     &                       'munu4nbo: g_munu_rot1',0,0,g_munu_rot1)) 
     $    call errquit('munu4nbo: g_munu_rot1', 0,GA_ERR)
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
      if (.not.nga_create(mt_dbl,3,ahi,'g_munuEPRHpar2d matrix',
     &    alo,g_munuEPRHpar2d)) call 
     &    errquit('munu4nbo: nga_create failed g_munuEPRHpar2d',
     &            0,GA_ERR)
        call ga_zero(g_munuEPRHpar2d)
      if (.not.nga_create(mt_dbl,3,ahi,'g_munuEPRpar12d matrix',
     &    alo,g_munuEPRpar12d)) call 
     &    errquit('munu4nbo: nga_create failed g_munuEPRpar12d',
     &            0,GA_ERR)
        call ga_zero(g_munuEPRpar12d)
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &              ' NW-g-shift: Summary C+HC data [ppt] ',
     &              'oooooooooooooooooo START'
        do n=1,3  ! xx,yy,zz
         m=n ! For principal components ONLY
         call ga_zero(g_acc)
         call ga_zero(g_acc3)
c ----- Do: A'= T^t A T, calculate only [A']_pp --> (do n=1,3 m=n)
c       a_pp'=    \sum_i t_ip a_ii t_ip + 
c               2 \sum_{j=2}^n \sum_{i=1}^{j-1} t_jp a_ji t_ip
c       g_munu_rot = A'
c       WARNING: g_munu_rot, contains several rotated matrices
c                since the matrices are symmetric I store only
c                the main diagonal + lower (upper) triangular 
c                matrix in a format that looks like:
c                a_11 a_22 ... a_nn 
c                a_21
c                a_31 a_32
c                a_41 a_42 a_43
c                ...
c                a_n1 a_n2 ... a_{n(n-1)}
c      There are two additional transformations on g_munu_rot
c      before leaving this routine and entering wefgfile()
c      1. I make the diagonalized matrix traceless
c ===== Transform xx_munu to 2xx_munu-(yy_munu+zz_munu) = START
c       or                   3xx_munu-(xx_munu+yy_munu+zz_munu)
c      2. I need to do a reordering of elements so that it is
c         compatible with wefgfile()
c        call reorder_munu(g_munu_rot,nlist,nlst,nbf,Ndir_munu)
c --------------------------------------------------------------
         do s=1,6
c ------- get coeff() --- START
          iind(1)=1
          iind(2)=1
          jind(1)=3*(arr_ind(s,1)-1)+m
          jind(2)=3*(arr_ind(s,2)-1)+n   
          call ga_gather(g_tvec,tmn,iind,jind,2)
          fact=1.0d0
          if (s.gt.3) fact=2.0d0
          coeff=fact*tmn(1)*tmn(2)      
c ------- get coeff() --- END
c Note.- g_munuFCSD will be the (hyp-diag)_uv matrix
           jlo=nsize*(s-1)+1
           jhi=jlo+nsize-1
           call ga_copy_patch('n',g_munuEPRdia,1,1,jlo,jhi,
     &                            g_tnp       ,1,1,1  ,nsize)
           call ga_add(1.0d0,g_acc,coeff,g_tnp,g_acc)  

           call ga_copy_patch('n',g_munuEPRpar1,1,1,jlo,jhi,
     &                            g_tnp        ,1,1,1  ,nsize)
           call ga_add(1.0d0,g_acc3,coeff,g_tnp,g_acc3)    
         enddo ! end-loop-s
c Note.- g_acc = \widetilde{H}_{mu nu}^{(m,m)}
c        it is the rotated munu matrix using:  T^t H T
         call ga_zero(g_acc1)
         do s=1,3
c ------- get coeff() --- START
          iind(1)=1
          jind(1)=3*(s-1)+m
          call ga_gather(g_tvec,tmn,iind,jind,1)
c          if (ga_nodeid().eq.0) then
c           write(*,1) m,s,tmn(1)
c 1         format('(m,s,tvec)=(',i5,',',i5,',',f15.8,')')
c          endif
          coeff=tmn(1) 
c ------- get coeff() --- END
c Note.- g_munuEPRHpar will be the (g-shift-para)_uv matrix
          jlo=nsize*(s-1)+1
          jhi=jlo+nsize-1
          call ga_copy_patch('n',g_munuEPRHpar,1,1,jlo,jhi,
     &                           g_tnp1       ,1,1,1  ,nsize)
          call ga_add(1.0d0,g_acc1,coeff,g_tnp1,g_acc1)   
c ----- Calculate rotated perturbed MO: g_acc2 ----- START
c  \sum_{s=1,3} t_{sm} C_{ri}^{(s) sigma} --> g_acc2
          elo(1) = 1
          ehi(1) = nbf
          elo(2) = 1
          ehi(2) = ntot
          elo(3) = s
          ehi(3) = s
          flo(1) = 1
          fhi(1) = nbf
          flo(2) = 1
          fhi(2) = ntot
          flo(3) = m
          fhi(3) = m
         call nga_add_patch(1.0d0,g_acc2,flo,fhi,
     &                       coeff,g_c1  ,elo,ehi,
     &                             g_acc2,flo,fhi)
c          if (s.eq.m) then
c           call nga_add_patch(1.0d0,g_acc2,flo,fhi,
c     &                        1.0d0,g_c1  ,elo,ehi,
c     &                              g_acc2,flo,fhi)
c          endif
c ----- Calculate rotated perturbed MO: g_acc2 ----- END 
         enddo ! end-loop-s
c Note: g_acc1 = \widetilde{H}_{mu nu}^{(m)}  m=x,y,z
c       it is the rotated munu matrix using: T H
c ====== Store final munu matrices === START
         jlo2=nlst*(n-1)+1
         jhi2=jlo2+nlst-1
         call ga_copy_patch('n',g_acc     ,1,1,   1,nlst,
     &                          g_munu_rot,1,1,jlo2,jhi2)   
         call ga_copy_patch('n',g_acc3    ,1,1,   1,nlst,
     &                         g_munu_rot2,1,1,jlo2,jhi2)         
         call ga_copy_patch('n',g_acc1    ,1,1,   1,nlst,
     &                         g_munu_rot1,1,1,jlo2,jhi2)     
c ====== Store final munu matrices === END
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== START
c ==== sum (g_acc .* g_sdens1 + Nuclear CONTRIB) 
c      = TOTAL HYP diagonalized   
         jlo1=1+nbf
         jhi1=nsize
         call ga_scale_patch(g_acc,1,1,jlo1,jhi1,2.0d0)
         chcdata(m)=ga_ddot(g_acc,g_sdens1)
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
        enddo ! end-loop-n
        if (ga_nodeid().eq.0) then
          write(*,23) 1,
     &                chcdata(1),            ! dia-x
     &                chcdata(2),chcdata(3)  ! dia-y,z
 23       format(' CHC   dia(xx,yy,zz)(',i3,')=(', 
     &           f15.8,',',f15.8,',',f15.8,')') 
        endif

c ++++++++++ CHECK diagonalization in para hyperfine ++++++++ START
c Note.- Variables defined in zora.fh:
c        g_CiFull, zora scaling factors 
c                 filled out with values in dft_zora_scale
c        zora switches: do_zora,do_NonRel,not_zora_scale
        type_NMR=3 ! =1,2,3=shieldings,hyperfine,gshift
        iat1=1
        call get_P10_rot(
     &       g_p10,            ! out: Perturbed density matrix (munu nbf x nbf x 3 square matrix)
     &       type_NMR,         !  in: =1,2,3=shieldings,hyperfine,gshift
     &       g_acc2,           !  in: rotated perturbed MO vector
     &       vectors,g_CiFull, !  in: to build zora scaled MO vector 
     &       iat1,             !  in: index for selected atom nr =1,nlist
     &       nbf,nmo,npol,nocc,nvirt,
     &       do_zora,do_NonRel,not_zora_scale,rtdb) 
      call fill_munuPSOSO_1(   ! g_munuEPRHpar --> g_munuEPRHpar2d
     &        g_munu_rot1    , ! in: array with unique elements
     &        g_munuEPRHpar2d, !out: nbf x nbf x 3 munu matrix for ith atom
     &        1,               ! in: one single unit
     &        2,               ! in: type_symm = 1 symm  = 2 antisymm
     &        nbf) 
      call fill_munuPSOSO_1(   ! g_munuEPRpar1 --> g_munuEPRHpar2d
     &        g_munu_rot2    , ! in: array with unique elements
     &        g_munuEPRpar12d, !out: nbf x nbf x 3 munu matrix for ith atom
     &        1,               ! in: one singlel unit
     &        1,               ! in: type_symm = 1 symm  = 2 antisymm
     &        nbf) 

c +++++++++ NOW: do ddot product to get diagonalized tensor +++ START
      call get_par_gshift_rot(
     &                   g_sdens,         ! IN : spin-density
     &                   g_munuEPRpar12d, ! IN : par 1st term
     &                   g_munuEPRHpar2d, ! IN : h01 matrix
     &                   g_p10,           ! IN : Perturbed density matrix
     &                   basis,nbf,1,rtdb)
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &              ' NW-g-shift: Summary C+HC data [ppt] ',
     &              'oooooooooooooooooo END'
c +++++++++ NOW: do ddot product to get diagonalized tensor +++ END
c ++++++++++ CHECK diagonalization in para hyperfine ++++++++ END
c --> Main outputs: g_acc2     ,rotated perturbed MO nbf*ntot*ndir*nlist
c                               ndir=1,2,3=x,y,z
c                               nbf, nr of basis functions
c                               ntot=nocc(1)+nocc(2)
c                               nlist, nr of selected atoms
c                   g_munu_rot1,rotated perturbed AO matrix
c                               storing only diag + off-diag elements
c                               Reminder: this comes from an antisymmetrix matrix
c                                         in case we want to pull back the 2d munu-matrix
c                   g_munu_rot, rotated AO matrix for dia part
c                               storing only diag + off-diag elements
c                               Reminder: this comes from a  symmetric matrix
c                                         in case we want to pull back the 2d munu-matrix
      nlist=1 ! ONLY for g-shift
      call reorder_munu(g_munu_rot ,nlist,nlst,nbf,Ndir_munu) ! reoder-munu matrix
      call reorder_munu(g_munu_rot1,nlist,nlst,nbf,Ndir_munu) ! reoder-munu matrix
      call reorder_munu(g_munu_rot2,nlist,nlst,nbf,Ndir_munu) ! reoder-munu matrix
c ------ destroy unnecessary GAs 
      if (.not. ga_destroy(g_munuEPRdia)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuEPRpar1)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuEPRHpar)) call errquit(
     &  'create_munu4nbo-2: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuEPRHpar2d)) call errquit(
     &  'create_munu4nbo-7a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuEPRpar12d)) call errquit(
     &  'create_munu4nbo-7: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp)) call errquit(
     &  'create_munu4nbo-5: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc)) call errquit(
     &  'create_munu4nbo-6: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc3)) call errquit(
     &  'create_munu4nbo-6: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_tnp1)) call errquit(
     &  'create_munu4nbo-5a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc1)) call errquit(
     &  'create_munu4nbo-6a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_sdens)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_sdens1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_p10)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_c1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('create_munu4nbo: ga_destroy failed vectors',0,GA_ERR)
      enddo
      call wgshiftfile(rtdb,
     &                 g_munu_rot,  ! dia term 
     &                 g_munu_rot1, ! perturbed AO operator g_munuEPRHpar x,y,z
     &                 g_acc2,      ! perturbed MO vector x,y,z
     &                 g_munu_rot2, ! 1st term in para      g_munuEPRpar1
     &                 nlst,       
     &                 Ndir_munu)

      if (.not. ga_destroy(g_munu_rot)) call errquit( ! destroy GA - FA
     &  'wgshiftfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munu_rot1)) call errquit( ! destroy GA - FA
     &  'wgshiftfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munu_rot2)) call errquit( ! destroy GA - FA
     &  'wgshiftfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc2)) call errquit(      ! destroy GA - FA
     &  'wgshiftfile: ga_destroy failed ',0, GA_ERR)
      return
      end  

      subroutine get_par_gshift_rot(  ! for checking rotated Pert. MO and rotated munu AO
     &                   g_sdens,     ! IN  : spin density
     &                   g_par1,      ! IN  : 1st term in para 
     &                   ga_h01,      ! IN  :
     &                   g_d1,        ! IN  : Perturbed density matrix 
     &                   basis,       ! IN  : basis handle
     &                   nbf,         ! IN  : nr. basis functions
     &                   iat,         ! IN  : atom nr
     &                   rtdb)
c Purpose : Assemble NMR Hyperfine: paramagnetic tensor
c Author  : Fredy Aquino
c Date    : 07-08-11
c Note    : This routine is taken from get_par_hfine()
c           The only difference is that here natoms=1
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "bas.fh"
#include "stdio.fh"
      integer g_d1 ! input: Perturbed density matrices 
      integer ga_h01,g_par1,g_sdens
      integer rtdb
      integer basis
      integer alo(3),ahi(3),ld(2)
      integer iat,ix,iy,ind,nbf
      double precision val(3),val1(3),chcdata(3)
      logical oskel

       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
        do iy = 1, 3
         alo(3) = iy
         ahi(3) = iy
         val1(iy)= nga_ddot_patch(g_sdens,'n',alo,ahi,
     &                            g_par1 ,'n',alo,ahi) 
          val(iy)= nga_ddot_patch(g_d1   ,'n',alo,ahi,
     &                            ga_h01 ,'n',alo,ahi) 
         chcdata(iy)=val(iy)+val1(iy)
        enddo ! end-loop-iy
        if (ga_nodeid().eq.0) then
          write(*,1) iat,
     &               chcdata(1),            ! par-x
     &               chcdata(2),chcdata(3)  ! par-y,z
 1        format(' CHC   par(xx,yy,zz)(',i3,')=(', 
     &            f15.8,',',f15.8,',',f15.8,')') 
        endif
      return
      end
c -------------- for g-shift NMLO analysis ----------------- END
c $Id$
