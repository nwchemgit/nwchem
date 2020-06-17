      subroutine hnd_giaox_zora(rtdb,basis,geom)
c $Id$
c
      implicit none
c
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
#include "inp.fh"
#include "zora.fh"
c
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer l_dia, k_dia, l_para, k_para
      integer l_xyz, k_xyz, l_zan, k_zan
      integer l_AtNr, k_AtNr,type_NMR
      integer l_fukui,k_fukui ! to include Fukui term
                              ! allows dia,para tensors 
                              ! to be origin invariant separately
      integer g_dens(3), g_s10, g_h01, g_h11,
     &        g_d1, g_rhs, g_rhs0,g_fock, g_u
      integer g_fock_Coul,g_fock_Exch,
     &        g_rhs_Coul,g_rhs_Exch,g_rhs_noJK,
     &        g_rhs_eSji,g_rhs_1e      
      integer g_d1_oo,g_d1_ov ! split occ-occ,occ-virt contrib to 
                              ! perturbed density matrix, P10 01-26-11
      integer g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji,g_d1_ov_ExchSFB

      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer vectors_scl(2) ! to be scaled with g_CiFull 
      double precision atn, tol2e, val, isotr, aniso
      double precision val_oo,val_ov
      double precision a(6),axs(3,3),eig(3),xfac
      double precision jfac(12),kfac(12)
      character*255 zorafilename
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      integer ld(2),cbuf ! FA-08-19-10
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
      external giao_aotomo,hnd_giasym,mat_transpose, ! FA-09-16-10
     &         dft_zoraNMR_read,get_slctd_atoms
      logical dft_zoraNMR_read
      double precision val1
      logical  oskel, status,debug_cs
      double precision val_ov_Coul,val_ov_Exch,val_ov_noJK,
     &                 val_ov_1e,val_ov_eSji

      double precision ppm
      data ppm     /26.62566914d+00/ 
      data tol2e   /1.0d-10/
      integer ic,npol 
c ----- Definitions for NLMO analysis ---- START
      integer i1,j1,acc_vec,l_tvec,k_tvec,shldfile          
      integer g_tvec ! for nbo analysis
c ----- Definitions for NLMO analysis ---- END
      external get_prelim_fock,get_P10,add_H10,add_fock,
     &         get_P10_1,skip_cphf,get_par,get_dia,
     &         skip_cphf_JK,get_P10_JK,new_giao_2e_JK,
     &         get_par_JK,
     &         dft_zoraCPHF_write,dft_zoraCPHF_read
      integer ntot,ispin,indA,indB,nocc(2),nind_jk,
     &        typeprop,nat_slc,nat
      integer cdens,l_pararr,k_pararr,icalczora,type_nmrdata
      integer ga_h01,npar_analysis
      integer debug_giaox
      integer shift,disp,nbf_dat,nat_dat
      logical switch_nmrcs_analysis,switch_skip_cphf
      integer g_AtNr1
      integer ga_dia,    ! OUTPUT
     &        ga_para1,  ! OUTPUT
     &        ga_h01_num,! OUTPUT
     &        ga_Fji     ! OUTPUT
      logical skip_cphf_ev_shield
      oskel=.false.
c Note.- switch_nmrcs_analysis=.true. has to go together with
c        switch_skip_cphf=.true.
c -----------------------------------
      debug_giaox=0 ! =1 for debugging print outs of matrices
c     debug_giaox=1 ! =1 for debugging print outs of matrices
c     switch_skip_cphf=.true. ! For skipping cphf or cpks
      switch_skip_cphf=.false. ! For NOT skipping cphf or cpks
      switch_nmrcs_analysis=.false. ! For default mode NO nmrcs analysis
c      switch_nmrcs_analysis=.true.  ! For analysis of nmrcs tensor
        if (ga_nodeid().eq.0.and.debug_giaox.eq.1)
     &    write(*,*) 'switch_skip_cphf=',switch_skip_cphf
        if (ga_nodeid().eq.0.and.debug_giaox.eq.1)
     &    write(*,*) 'switch_nmrcs_analysis=',switch_nmrcs_analysis
c
c     Print NMR shielding header
c
      if (ga_nodeid().eq.0) write(luout,9999)
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
c     
c     If DFT get part of the exact exchange defined
c
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()    
      nind_jk=12
      do ifld = 1,nind_jk
        jfac(ifld) =  0.0d0       ! used in shell_fock_build()
        kfac(ifld) = -1.0d0*xfac  ! used in shell_fock_build()
        if (ga_nodeid().eq.0.and.debug_giaox.eq.1) then
         write(*,144) ifld,jfac(ifld),kfac(ifld)
  144    format('(j,k)(',i3,')=(',f15.8,',',f15.8,')')
        endif
      enddo
c
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
c ------- Read (nat,atmnr) --------- START
         status=geom_ncent(geom,nat)   
      if (.not.ma_alloc_get(
     &       mt_int,nat,'nmt tmp',l_AtNr,k_AtNr))
     &    call errquit('hnd_giaox_zora: ma failed',0,MA_ERR)
         typeprop=2 ! =1 EFG =2 Shieldings =3 Hyperfine  
         call get_slctd_atoms(nat_slc,       ! out: selected atoms
     &                        int_mb(k_AtNr), ! out: list of selected atom nr.     
     &                        nat,           ! in : total nr atoms in molecule            
     &                        rtdb,          ! in : rdt  handle
     &                        typeprop)      ! in : =1,2,3=EFG,Shieldings,Hyperfine
      if (ga_nodeid().eq.0.and.debug_giaox.eq.1) then
       write(*,*) 'nat_slc=',nat_slc
       do i=1,nat_slc
        write(*,7) i,int_mb(k_AtNr+i-1)
 7      format('atomnr(',i3,')=',i5)
       enddo
      endif
c ------- Read (nat,atmnr) --------- END
      if (.not. ma_push_get(mt_dbl,3*nat_slc,'nmr at',l_xyz,k_xyz)) 
     &    call errquit('hnd_giaox_zora: ma_push_get failed k_xyz',
     &                  0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat_slc,'nmr zan',l_zan,k_zan)) 
     &    call errquit('hnd_giaox_zora: ma_push_get failed k_zan',
     &                  0,MA_ERR)
c
c     Try to read the atom list from rtdb. If it is not there, we still have the default list
      do ixy = 0, nat_slc-1
         if (.not. geom_cent_get(geom, int_mb(k_AtNr+ixy), tag, 
     &                           dbl_mb(k_xyz+3*ixy),dbl_mb(k_zan+ixy)))
     &       call errquit('hnd_giaox_zora: geom_cent_tag failed',
     &                     0, GEOM_ERR)
      enddo 
c
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_giaox_zora: could not get nbf',0, BASIS_ERR)

c ++++++ Reading shieldings data from file ++++++ START
      if (do_zora) then ! ==========if-do_zora-START
c       Note.- lbl_nmrcs defined in zora.fh
        call util_file_name(lbl_nmrcs,.false.,.false.,zorafilename)
       icalczora = 0  ! initialize the flag
       type_nmrdata=1 ! =1,2,3=shieldings,hyperfine,gshift      
       if (.not.dft_zoraNMR_read(zorafilename,
     &     type_nmrdata,
     &     nbf,nat_slc,
     &     g_AtNr1,ga_dia,
     &     ga_para1,ga_h01_num,ga_Fji)) icalczora=1 
c Note.- If I print the GAs here it gets freezed
      endif ! ======================if-do_zora-END
c ++++++ Reading shieldings data from file ++++++ END
      if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_para1---------- START'
         call ga_print(ga_para1)
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_para1---------- END'
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_h01_num---------- START'
         call ga_print(ga_h01_num)
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_h01_num---------- END'
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_Fji---------- START'
         call ga_print(ga_Fji)
        if (ga_nodeid().eq.0)
     &   write(*,*) '----ga_Fji---------- END'
      endif
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('hnd_giaox_zora: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('hnd_giaox_zora: ma_push_get failed k_eval',0,MA_ERR)
c
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c
c
c ------ define npol ----- START
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, npol)) then
        if      (scftyp .eq. 'UHF') then
         npol=2
        else if (scftyp .eq. 'RHF') then
         npol=1
        endif
      endif
c ------ define npol ----- END
c     if (ga_nodeid().eq.0)
c    & write(*,*) 'npol=',npol

      if      (npol.eq.1) then
       nocc(1)=nclosed(1)
       nocc(2)=0
      else if (npol.eq.2) then
       nocc(1)=nopen(1)   
       nocc(2)=nopen(2)     
      endif
      if (debug_giaox.eq.1) then
       write(*,10) nocc(1),nocc(2),
     &             nopen(1),nopen(2),
     &             nclosed(1),nclosed(2),
     &             nvirt(1),nvirt(2),scftyp,nmo
 10    format('nocc =(',i3,',',i3,') ',
     &        'nopen=(',i3,',',i3,') ',
     &        'nclos=(',i3,',',i3,') ',
     &        'nvirt=(',i3,',',i3,') scftyp=',a,', nmo=',i3)
      endif
c ------ Store nopen in rtdb so that CPHF routine is happy ---- START
c WARNING: For restricted calc nocc(1)=9 nocc(2)=0 <=== PROBLEM
        if (npol.eq.2) then   
          if (.not. rtdb_put(rtdb, 'scf:nopen', 
     &         MT_INT, 1, nocc(1)-nocc(2)))
     *         call errquit('hnd_giaox_zora:rtdbput nopen failed',
     &         nocc(1)-nocc(2),
     &       RTDB_ERR)
        endif
c ------ Store nopen in rtdb so that CPHF routine is happy ---- END
c
c     Get Unperturbed Density Matrix
      call hnd_prp_get_dens(rtdb,geom,basis,
     &                       g_dens, ! out
     &                       ndens,scftyp,
     &                       nclosed,nopen,nvirt)
c
      ntot=0
      do ispin=1,npol
        ntot=ntot+nocc(ispin)*nvirt(ispin)
      enddo
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)

      if (switch_nmrcs_analysis) then
        if (ga_nodeid().eq.0)
     &    write(*,*) 'ENTER-define g_rhs(J,K,nJK,eSji'
c ----- TEST: for testing Coulomb and Exchange contrib --- START
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs_Coul))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhsJ',
     &                 0,GA_ERR)
      call ga_zero(g_rhs_Coul)
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs_Exch))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhsK',
     &                 0,GA_ERR)
      call ga_zero(g_rhs_Exch)
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs_noJK))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhsnJK',
     &                 0,GA_ERR)
      call ga_zero(g_rhs_noJK)
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs_eSji))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhseSji',
     &                 0,GA_ERR)
      call ga_zero(g_rhs_eSji)
      if(.not.ga_create(MT_DBL,ntot,3,'RHS',-1,-1,g_rhs_1e))
     &   call errquit('hnd_giaox_zora: ga_create failed g_rhs1e',
     &                 0,GA_ERR)
      call ga_zero(g_rhs_1e)
c ----- TEST: for testing Coulomb and Exchange contrib --- END
      endif

      if (debug_giaox.eq.1) then
      write(*,*) 'after-creating g_rhs ...'
      write(*,102) npol,nocc(1),nocc(2),
     &            nclosed(1),nclosed(2),
     &            nvirt(1),nvirt(2),scftyp,ntot
 102   format('BEF pre-fock::npol=',i3,' nocc=(',i3,',',i3,') ',
     &       'nclos=(',i3,',',i3,') ',
     &       'nvirt=(',i3,',',i3,') scftyp=',a,' ntot=',i3)
      endif

      if (switch_nmrcs_analysis) then
       call get_prelim_fock_debug( 
     &                            g_d1, ! out: 
     &                           g_rhs, ! out: rhs expression
     &                          g_rhs0, ! out: to be used in get_d1()
     &                      g_rhs_eSji, ! out: -Sji^k \epsilon_i
     &                         vectors, !  in: MO  coeffs
     &                  dbl_mb(k_eval), !  in: energy vals
     &                   dbl_mb(k_xyz), !  in: Nuclear positions (x,y,z)
     &                         nat_slc, !  in: nr. selected nuclei (atoms)
     &                           basis, !  in: basis handle
     &                             nbf, !  in: nr. basis functions
     &                             nmo, !  in: nr. MOs (occ+virt)
     &                            npol, !  in: nr. of polarizations
     &                            nocc, !  in: nr. occ     MOs
     &                           nvirt) !  in: nr. virtual MOs
      else ! default
       call get_prelim_fock(g_d1,  ! out: 
     &                      g_rhs, ! out: rhs expression
     &                      g_rhs0,! out: to be used in get_d1()
     &                      vectors,dbl_mb(k_eval),
     &                      dbl_mb(k_xyz),nat_slc,basis, 
     &                      nbf,nmo,npol,nocc,nvirt) 
      endif
c
      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-aft-get_prelim_fock -------- START'
       call ga_print(g_rhs)
       call ga_print(g_rhs0)
       call ga_print(g_d1)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-aft-get_prelim_fock -------- START'
      endif 
c
c     Build "first order fock matrix"
c
      if (use_theory.eq.'dft') then
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hnd_giaox_zora: rtdb_put of xc_active failed',
     &       0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hnd_giaox_zora: rtdb_put of calc_type failed',
     &       0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hnd_giaox_zora: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
      endif
      clo(1) = 3*npol*2  
      clo(2) = nbf
      clo(3) = nbf
      chi(1) =  1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',chi,g_fock)) call 
     &    errquit('hnd_giaox_zora: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
c
c     Note: Just the exchange: jfac = 0.d0 (see above)
c
      if (.not.cam_exch) then
         call shell_fock_build(geom, basis,0,3*npol*2,
     $                         jfac,kfac,tol2e,
     &                         g_d1,  ! input
     &                         g_fock,! output
     &                         .false.)
      else
          call shell_fock_build_cam(geom, basis,0,3*npol*2,
     $                              jfac,kfac,tol2e,
     &                              g_d1,  ! input
     &                              g_fock,! output
     &                              .false.)
      end if
c
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('hnd_giaox_zora: rtdb_put failed',0,RTDB_ERR)
      endif
      if      (npol.eq.1) then
       nocc(1)=nclosed(1)
       nocc(2)=0
      else if (npol.eq.2) then
       nocc(1)=nopen(1)   
       nocc(2)=nopen(2)     
      endif
      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-bef-add_fock -------- START'
        call ga_print(g_rhs)
        call ga_print(g_fock)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-bef-add_fock -------- END'
      endif
      call add_fock(g_rhs, ! out: g_rhs=g_rhs+g_fock
     &              g_fock,vectors,
     &              nbf,nmo,npol,nocc,nvirt) 

       if (debug_giaox.eq.1) then
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
     &    errquit('hnd_giaox_zora: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed g_fock',0,GA_ERR)

      if (debug_giaox.eq.1) then
       write(*,*) '---- g_rhs-bef-add_H10 -------- START'
       call ga_print(g_rhs)
       write(*,*) '---- g_rhs-bef-add_H10 -------- END'
      endif
c ------ WARNING: skip new_giao_2d() ---- START
c       if (ga_nodeid().eq.0)
c    &   write(*,*) 'WARNING: SKIP add_H10() ...'
c       goto 178
c ------ WARNING: skip new_giao_2d() ---- END

      if (switch_nmrcs_analysis) then
       call  add_H10_debug( 
     &                    g_rhs, ! out: accumulated rhs expression
     &                 g_rhs_1e, ! out : Fji^{k,1e}
     &                   ga_Fji, !  in: Fock 1st-deriv without V (pot.) contrib.
     &                  vectors, !  in: MO  coeffs
     &            dbl_mb(k_xyz), !  in: Nuclear positions (x,y,z)
     &                  nat_slc, !  in: nr. selected nuclei (atoms)
     &                    basis, !  in: basis handle
     &                      nbf, !  in: nr. basis functions
     &                      nmo, !  in: nr. MOs (occ+virt)
     &                     npol, !  in: nr. of polarizations
     &                     nocc, !  in: nr. occ     MOs
     &                    nvirt, !  in: nr. virtual MOs
     &                  do_zora, !  in: switch = .true. zora on
     &                     rtdb) !  in: rtdb handle
      else
       call add_H10(g_rhs, ! out: ga_rhs(a,i)=ga_rhs(a,i)+H10(a,i)
     &             ga_Fji,vectors, 
     &             dbl_mb(k_xyz),nat_slc,basis, 
     &             nbf,nmo,npol,nocc,nvirt,do_zora,rtdb) 
      endif
c178  continue

       if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,*) 'Aft. add_H10()'
         write(*,*) '---- g_rhs-aft-add_H10 -------- START'
        endif
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-aft-add_H10 -------- END'
       endif
c ------------ added-4 ------------ END
c     Remaining term is Perturbed (GIAO) two-electron term times Unperturbed density
c     Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in AO basis
      alo(1) = -1 
      alo(2) = -1
      alo(3) = 1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*npol
      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
     &    errquit('hnd_giaox_zora: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
      if (switch_nmrcs_analysis) then    
        if (ga_nodeid().eq.0)
     &    write(*,*) 'ENTER-define g_fock_Coul,g_fock_Exch'
c -- TEST: create g_fock_Coul,g_fock_Exch ----------- START
       if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',
     &    alo,g_fock_Coul)) call 
     &    errquit('hnd_giaox_zora: nga_create failed g_fockJ',0,GA_ERR)
       call ga_zero(g_fock_Coul)
       if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',
     &    alo,g_fock_Exch)) call 
     &    errquit('hnd_giaox_zora: nga_create failed g_fockX',0,GA_ERR)
       call ga_zero(g_fock_Exch)
c -- TEST: create g_fock_Coul,g_fock_Exch ----------- END
      endif
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('hnd_giaox_zora: rtdb_put failed',0,RTDB_ERR)
      endif
c +++++++++++++++++++++++++++++++++++++
c ++++++++ calling new_giao_2e() ++++++
       if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_dens bef new_giao_2e  -------- START'
         do cdens=1,npol
          call ga_print(g_dens(cdens))
         enddo
        if (ga_nodeid().eq.0)
     &    write(*,*) '---- g_dens bef new_giao_2e  -------- END'    
       endif
       call ga_sync()
c ------ WARNING: skip new_giao_2d() ---- START
c       if (ga_nodeid().eq.0)
c    &   write(*,*) 'WARNING: SKIP new_giao_2e() ...'
c      goto 176
c ------ WARNING: skip new_giao_2d() ---- END

      if (switch_nmrcs_analysis) then
c +++++ Test: getting Coulomb and Exchange contrib separate --- START
       call new_giao_2e_JK(geom,basis,nbf,tol2e,
     &                 g_dens, !  in: e-denstiy 
     &                 g_fock, ! out: fock matrix
     &                 g_fock_Coul,
     &                 g_fock_Exch,
     &                 xfac,
     &                 npol)
c +++++ Test: getting Coulomb and Exchange contrib separate --- END
      else  ! default
       call new_giao_2e(geom,basis,nbf,tol2e,
     &                 g_dens, !  in: e-denstiy 
     &                 g_fock, ! out: fock matrix
     &                 xfac,
     &                 npol)
      endif
  176   continue

c ++++++++ calling new_giao_2e() ++++++
c +++++++++++++++++++++++++++++++++++++

       if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-new_giao -------- START'
         call ga_print(g_fock)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-new_giao -------- END'
       endif

      if(use_theory.eq.'dft') then
         ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('hnd_giaox_zora: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $      call errquit('hnd_giaox_zora: rtdb_put of xc_active failed',
     &      0,RTDB_ERR)
      endif
      if (debug_giaox.eq.1) then
       write(*,114) npol,nocc(1),nocc(2),
     &              nclosed(1),nclosed(2),
     &              nvirt(1),nvirt(2),scftyp
 114   format('BEF. giao_aotomo: npol=',i3,
     &       ' nocc=(',i3,',',i3,') ',
     &       'nclos=(',i3,',',i3,') ',
     &       'nvirt=(',i3,',',i3,') scftyp=',a,')')
      endif
c     Transform to MO basis and add to right-hand-side
      call giao_aotomo(g_fock,vectors,nocc,nvirt,npol,3,nbf)
      if (switch_nmrcs_analysis) then
c ----- TEST: AO2MO for g_fock_Coul,g_fock_Exch ------ START
       call giao_aotomo(g_fock_Coul,vectors,nocc,nvirt,npol,3,nbf)
       call giao_aotomo(g_fock_Exch,vectors,nocc,nvirt,npol,3,nbf)
c ----- TEST: AO2MO for g_fock_Coul,g_fock_Exch ------ END
      endif
      if (debug_giaox.eq.1) then
       write(*,*) 'Aft add_fock1-giao-aotomo()'
       write(*,*) 'BEF. giao_aotomo: npol='
       write(*,*) '---- g_fock-aft-giao_aotomo -------- START'
       call ga_print(g_fock)
       write(*,*) '---- g_fock-aft-giao_aotomo -------- END'
      endif

      if (switch_nmrcs_analysis) then
c ------- Store g_rhs without Coulomb and Exchange ---- START
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
c ------- Store g_rhs without Coulomb and Exchange ---- END
      endif

      call add_fock1(g_rhs, ! out: accumulated rhs expression
     &               g_fock,! in 
     &               nmo,npol,nocc,nvirt) 

      if (switch_nmrcs_analysis) then
       call add_fock1(g_rhs_Coul, ! out: accumulated rhs expression
     &                g_fock_Coul,! in 
     &                nmo,npol,nocc,nvirt) 
       call add_fock1(g_rhs_Exch, ! out: accumulated rhs expression
     &                g_fock_Exch,! in 
     &                nmo,npol,nocc,nvirt) 
       if (.not. ga_destroy(g_fock_Coul)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed g_fock_Coul',0, GA_ERR) 
       if (.not. ga_destroy(g_fock_Exch)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed g_fock_Exch',0, GA_ERR) 
      endif
      if (debug_giaox.eq.1) then
       write(*,*) 'Aft add_fock1()'
       write(*,*) '---- g_rhs-aft-add_fock1 -------- START'
       call ga_print(g_rhs)
       write(*,*) '---- g_rhs-aft-add_fock1 -------- END'
      endif
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed g_fock',0,GA_ERR)

      if (switch_skip_cphf) then
       if (switch_nmrcs_analysis) then    
        if (ga_nodeid().eq.0)
     &    write(*,*) 'WARNING : Using skip_cphf_JK ...'
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
     &   write(*,*) 'WARNING : Using skip_cphf ...'
            call skip_cphf(g_rhs,          ! IN/OUT
     &                     dbl_mb(k_eval), ! IN: energy eigenvalues
     &                     nocc,
     &                     nvirt,
     &                     nbf,  ! FA-04-28-12: replacing nmo by nbf
     &                     npol)
       endif
      endif

c -------free allocated memory -------------- START
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_occ',0,MA_ERR)
c -------free allocated memory -------------- END
      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- END'
      endif
      if(.not.rtdb_get(rtdb,'zora:skip_cphf_ev_shield',
     &                 mt_log,1,skip_cphf_ev_shield))        
     &  skip_cphf_ev_shield = .false.       
      if (.not.(switch_skip_cphf) .and.
     &    .not.(skip_cphf_ev_shield)) then
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'COMPUTE cphf shield data ...'
c     Write ga_rhs to disk 
       call cphf_fname('cphf_rhs',cphf_rhs)
       call cphf_fname('cphf_sol',cphf_sol)
       if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_giaox_zora: could not write cphf_rhs',0, DISK_ERR)
       call schwarz_tidy()
       call int_terminate()
c
c     Call the CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.

       if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('hnd_giaox_zora: failed to write skew ', 0, RTDB_ERR)
       if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_giaox_zora: failure in cphf ',0, RTDB_ERR)
       if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('hnd_giaox_zora: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
        call ga_zero(g_rhs)
       if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_giaox_zora: could not read cphf_rhs',0, DISK_ERR)  
c ----- write CPHF data to file ----------------- START
c       Note.- lbl_cphfhyp defined in zora.fh
       call util_file_name(lbl_cphfshield,
     &                     .false.,.false.,zorafilename)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs0-cs --------- START'
c        call ga_print(g_rhs0)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs0-cs --------- END'
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-cs --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-cs --------- END'
       call dft_zoraCPHF_write(
     &           zorafilename, ! in: filename
     &           npol,         ! in: nr polarization
     &           nocc,         ! in: nr occupied MOs
     &           nvirt,        ! in: nr virtual  MOs
     &           nbf,          ! in: nr basis functions
     &           vectors,      ! in: MOs
     &           g_rhs0,       ! in: (ntot,3)       GA matrix
     &           g_rhs)        ! in: (nocc*nvirt,3) GA matrix
c ----- write CPHF data to file ----------------- END
      else
       if (ga_nodeid().eq.0)
     &  write(*,*) 'WARNING: SKIP cphf ...'
       if (skip_cphf_ev_shield) then
        call ga_zero(g_rhs)
        do i=1,npol
         call ga_zero(vectors(i))
        enddo
        if (ga_nodeid().eq.0)
     &   write(*,*) 'READ cphf shieldings data from file ...'
        call util_file_name(lbl_cphfshield,.false.,.false.,zorafilename)
        call dft_zoraCPHF_read(
     &           zorafilename, ! in: filename
     &           npol,         ! in: nr polarization
     &           nocc,         ! in: nr occupied MOs
     &           nvirt,        ! in: nr virtual  MOs
     &           nbf,          !  in: nr basis functions
     &           vectors,      ! out: MOs
     &           g_rhs0,       ! out: (ntot,3)       GA matrix
     &           g_rhs)        ! out: (nocc*nvirt,3) GA matrix
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-cs-read --------- START'
c        call ga_print(g_rhs)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '---------- g_rhs-cs-read --------- END'
       endif
      endif ! end-if-switch-skip_cphf

      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- START'
      call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- END'  
      endif
      type_NMR=1 ! =1,2,3=shieldings,hyperfine,gshift
       if (switch_nmrcs_analysis) then
        call get_P10_JK(
     &          g_d1_oo, ! out: Perturbed density matrix occ-occ  contrib
     &          g_d1_ov, ! out: Perturbed density matrix occ-virt contrib
     &          type_NMR,
     &          g_d1_ov_Coul,! out: Perturbed (spin)density matrix with Coulomb  contrib
     &          g_d1_ov_Exch,! out: Perturbed (spin)density matrix with Exchange contrib (fock_xc)
     &          g_d1_ov_noJK,! out: Perturbed (spin)density matrix remaining terms not J or K
     &          g_d1_ov_1e,  ! out: Perturbed (spin)density matrix 1e-operator contrib
     &          g_d1_ov_eSji,! out: Perturbed (spin)density matrix perturbed overlap matrix
     &          g_rhs,       ! in: accumulated rhs expression
     &          g_rhs_Coul,
     &          g_rhs_Exch,
     &          g_rhs_noJK,
     &          g_rhs_1e,
     &          g_rhs_eSji,
     &          g_rhs0,      ! in: from get_prelim_fock()
     &          vectors,g_CiFull, 
     &          nbf,nmo,npol,nocc,nvirt,
     &          do_zora,do_NonRel,not_zora_scale,
     &          lbl_nlmogshift, ! in: for g-shift NLMO analysis
     &          lbl_nlmoshield, ! in: for shield  NLMO analysis
     &          rtdb) 
      if (.not. ga_destroy(g_rhs_Coul)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_rhs_Exch',0, GA_ERR)
      if (.not. ga_destroy(g_rhs_Exch)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_rhs_Exch',0, GA_ERR)
      if (.not. ga_destroy(g_rhs_noJK)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_rhs_noJK',0, GA_ERR)
      if (.not. ga_destroy(g_rhs_1e)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_rhs_1e',0, GA_ERR)
      if (.not. ga_destroy(g_rhs_eSji)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_rhs_1e',0, GA_ERR)
      else  ! default
c      call get_P10(  g_d1, ! out: Perturbed density matrix
c     &              g_rhs, !  in: accumulated rhs expression
c     &             g_rhs0, !  in: from get_prelim_fock()
c     &            vectors,g_CiFull, 
c     &            nbf,nmo,npol,nocc,nvirt,rtdb) 
         call get_P10_1( 
     &         g_d1_oo,        ! out: Perturbed density matrix occ-occ  contrib
     &         g_d1_ov,        ! out: Perturbed density matrix occ-virt contrib
     &         type_NMR,       ! in : =1,2,3=shieldings,hyperfine,gshift
     &         g_rhs,          ! in: accumulated rhs expression
     &         g_rhs0,         ! in: from get_prelim_fock()
     &         vectors,g_CiFull, 
     &         nbf,nmo,npol,nocc,nvirt,
     &         do_zora,do_NonRel,not_zora_scale,
     &         lbl_nlmogshift, ! in : label for g-shift NLMO analysis
     &         lbl_nlmoshield, ! in : label for shield  NLMO analysis
     &         rtdb) 
      endif

       if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed g_rhs',0,GA_ERR)
       if (.not.ga_destroy(g_rhs0)) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed g_rhs0',0,GA_ERR)

c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the shielding tensor
      if (.not. ma_push_get(mt_dbl,9*nat_slc,'sh para',l_para,k_para)) 
     &    call errquit('hnd_giaox_zora: ma_push_get failed k_para',
     &    0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*nat_slc,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_giaox_zora: ma_push_get failed k_dia',0,MA_ERR)
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)

c ------- TEST: get_par() get_dia()------ START    
      if (switch_nmrcs_analysis) then 
       npar_analysis=9
      else 
       npar_analysis=4
      endif
      if (.not. ma_push_get(mt_dbl,
     &          npar_analysis*3*3*nat_slc,'sh pararr',
     &          l_pararr,k_pararr)) 
     &    call errquit('hnd_giaox_zora: ma_push_get failed k_pararr',
     &          0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,
     &          3*3*nat_slc,'sh fukui',
     &          l_fukui,k_fukui)) 
     &    call errquit('hnd_giaox_zora: ma_push_get failed l_fukui',
     &          0,MA_ERR)
      call ga_sync()
      if (switch_nmrcs_analysis) then 
       call get_par_JK(
     &   dbl_mb(k_para),  ! OUT : paramagnetic tensor
     &   ga_para1,        ! IN  : paramagnetic tensor (GA) 
     &   ga_h01_num,      ! IN  : H01 matrix
     &   dbl_mb(k_pararr),! OUT : par tensor gauge,OO,OV
     &   g_d1_oo,         ! IN  : Perturbed density matrix (OO part)
     &   g_d1_ov,         ! IN  : Perturbed density matrix (OV part) 
     &   g_d1_ov_Coul,    ! IN  : Perturbed density matrix (OV part-Coul only) 
     &   g_d1_ov_Exch,    ! IN  : Perturbed density matrix (OV part-Exch only) 
     &   g_d1_ov_noJK,    ! IN  : Perturbed density matrix (OV part-all -(Coul,Exch)) 
     &   g_d1_ov_1e,      ! IN  : Perturbed density matrix (OV part-1e   contrib)
     &   g_d1_ov_eSji,    ! IN  : Perturbed density matrix (OV part-eSji contrib)
     &   g_dens,
     &   npol,        ! IN  : nr. polarizations
     &   dbl_mb(k_fukui),
     &   basis,           ! IN  : basis handle
     &   nbf,             ! IN  : nr. basis functions
     &   dbl_mb(k_xyz),   ! IN  : Nuclear positions (x,y,z)
     &   nat_slc,         ! IN  : nr. atoms selected
     &   nat,             ! IN  : nr. atoms total
     &   geom,
     &   rtdb,
     &   oskel)
      if (.not. ga_destroy(g_d1_ov_Coul)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_d1_ov_Coul',0, GA_ERR)
      if (.not. ga_destroy(g_d1_ov_Exch)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_d1_ov_Exch',0, GA_ERR)
      if (.not. ga_destroy(g_d1_ov_noJK)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_d1_ov_noJK',0, GA_ERR)
      if (.not. ga_destroy(g_d1_ov_1e)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_d1_ov_noJK',0, GA_ERR)
      if (.not. ga_destroy(g_d1_ov_eSji)) call errquit(  
     &  'hnd_giaox_zora: ga_destroy failed g_d1_ov_noJK',0, GA_ERR)
      else
       call get_par(
     &   dbl_mb(k_para),  ! OUT : paramagnetic tensor
     &   ga_para1,        ! IN  : paramagnetic tensor (GA) 
     &   ga_h01_num,      ! IN  : H01 matrix
     &   dbl_mb(k_pararr),! OUT : par tensor gauge,OO,OV
     &   g_d1_oo,         ! IN  : Perturbed density matrix (OO part)
     &   g_d1_ov,         ! IN  : Perturbed density matrix (OV part) 
     &   g_dens,          ! IN  : e-density
     &   npol,            ! IN  : nr. polarizations
     &   dbl_mb(k_fukui), ! OUT : Fukui term
     &   basis,           ! IN  : basis handle
     &   nbf,             ! IN  : nr. basis functions
     &   dbl_mb(k_xyz),   ! IN  : Nuclear positions (x,y,z)
     &   nat_slc,         ! IN  : nr. atoms selected
     &   nat,             ! IN  : nr. atoms total
     &   geom,
     &   rtdb,
     &   oskel)           ! IN  : = .false.
      endif

      call get_dia(dbl_mb(k_dia),   ! OUT: dia tensor
     &             ga_dia,          ! IN : ga dia tensor   
     &             basis,           ! IN : basis handle
     &             g_dens,          ! IN : e-density
     &             nbf,             ! IN : nr. basis functions
     &             npol,            ! IN : nr. polarizations
     &             dbl_mb(k_fukui), ! IN : Fukui term
     &             dbl_mb(k_xyz),   ! IN : nuclear positions
     &             nat_slc,         ! IN : selected atoms
     &             rtdb,            ! IN : rtdb handle
     &             oskel)           ! IN : = .false
c ------- TEST: get_par() get_dia()------ END
c +++++++++ print-total-pardia-transferred +++ START
      debug_cs=.false.
      if (debug_cs) then
       if (ga_nodeid().eq.0) then
      ic=1
       do iatom = 1, nat_slc
        do ix = 1, 3
         do iy = 1, 3
           cbuf=k_pararr+
     &          3*3*npar_analysis*(iatom-1)+
     &          3*npar_analysis*(iy-1)+npar_analysis*(ix-1)-1
          if (switch_nmrcs_analysis) then !-- START-if-switch_gshift_analysis
           write(*,179) ix,iy,iatom,
     &                 dbl_mb(k_dia+ic-1),
     &                 dbl_mb(cbuf+1),dbl_mb(cbuf+2),
     &                 dbl_mb(cbuf+3),
     &                 dbl_mb(cbuf+5),dbl_mb(cbuf+6),
     &                 dbl_mb(cbuf+7),
     &                 dbl_mb(cbuf+8),dbl_mb(cbuf+9),
     &                 dbl_mb(cbuf+4)
 179        format('NW:(dia,gauge,OO,OV,',
     &            'OV_Coul,OV_Exch,OV_nJK,',
     &            'OV_1e,OV_eSji,'
     &            'Totpar)(',i1,',',i1,',',i1,')=(',
     &             f18.6,' ',f18.6,' ',f18.6,' ',
     &             f18.6,' ',f18.6,' ',
     &             f18.6,' ',f18.6,' ',
     &             f18.6,' ',f18.6,' ',
     &             f18.6,' )')
          else
           write(*,19) ix,iy,iatom,
     &                 dbl_mb(k_dia+ic-1),
     &                 dbl_mb(cbuf+1),dbl_mb(cbuf+2),
     &                 dbl_mb(cbuf+3),dbl_mb(cbuf+4),
     &                 dbl_mb(k_dia+ic-1)+dbl_mb(cbuf+4)
 19        format('NW:(dia,gauge,OO,OV,Totpar,dia+par)(',
     &             i1,',',i1,',',i1,')=(',
     &             f18.6,' ',f18.6,' ',f18.6,' ',f18.6,' ',
     &             f18.6,' ',f18.6,' )')
          endif
          ic=ic+1
         enddo ! end-loop-iy
        enddo ! end-loop-ix
       enddo ! end-loop-iatom
       endif ! end-if-ga_nodeid-eq-0
      endif ! end-if-debug_cs
c +++++++++ print-total-dia-transferred +++ END
      if (.not.ma_pop_stack(l_fukui)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_fukui',
     &    0,MA_ERR)
      if (.not.ma_pop_stack(l_pararr)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_para',0,MA_ERR)
       do ispin=1,ndens
        if (.not.ga_destroy(g_dens(ispin))) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed g_dens',
     &    0,GA_ERR)
       enddo

c --------- allocating array for NLMO analysis --------------- START
      shldfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
      if (shldfile.eq.1) then
       if (.not. ma_alloc_get(mt_dbl,nat_slc*9,'tvec',l_tvec,k_tvec)) 
     &    call
     &    errquit('hnd_giaox_zora: ma_push_get failed tvec',0,MA_ERR)   
      endif ! end-if-shldfile
c --------- allocating array for NLMO analysis --------------- END
c
c     Print out tensor information, and write to Ecce file if necessary
      status = rtdb_parallel(.false.)   
      if (ga_nodeid().gt.0) goto 300
      acc_vec=0 ! For NMLO analysis
      call ecce_print_module_entry('nmr')
      do iatom = 1, nat_slc
         ioff = (iatom-1)*9
         if (.not. geom_cent_get(geom, int_mb(k_AtNr+iatom-1), tag, 
     &       dbl_mb(k_xyz), dbl_mb(k_zan))) call 
     &       errquit('hnd_giaox_zora: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(tag, symbol, element, atn)) then
           if (.not. inp_compare(0,tag(1:2),'bq')) call ! check for bq as a fall back
     &       errquit('hnd_giaox_zora: geom_tag_to_element failed',
     &               0,GEOM_ERR)
         endif
c
c      Print tensor pieces and sum for total shielding tensor
         if (ga_nodeid().eq.0) then
            write(luout,9700) iatom,symbol ! Showing original atm nr.
            write(luout,9800) (dbl_mb(k_dia+ioff+ix-1),ix=1,9)
            write(luout,9801) (dbl_mb(k_para+ioff+ix-1),ix=1,9)
         endif
         do ix = 0, 8 
            dbl_mb(k_para+ioff+ix) = dbl_mb(k_dia +ioff+ix) + 
     &                               dbl_mb(k_para+ioff+ix)
         enddo
c
c     Print total shielding tensor
c
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
            shldfile=0 ! not doing NLMO analysis by default
            status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
            if (shldfile.eq.1) then
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
c
            call ecce_print1_char('atom name',symbol,1)
            call ecce_print2('shielding tensor',MT_DBL,
     &                       dbl_mb(k_para+ioff),3,3,3)
            call ecce_print1('shielding isotropic'   ,MT_DBL,isotr,1)
            call ecce_print1('shielding anisotropy'  ,MT_DBL,aniso,1)
            call ecce_print1('shielding eigenvalues' ,MT_DBL,eig,3)
            call ecce_print2('shielding eigenvectors',MT_DBL,axs,
     &                       3,3,3)
         endif
      enddo
      call ecce_print_module_exit('nmr','ok')
300   call ga_sync()

      status = rtdb_parallel(.true.)   
      shldfile=0 ! not doing NLMO analysis by default
      status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
      if (shldfile.eq.1) then ! ------- hypfile-if++++ START
         if (.not. ga_create(mt_dbl,1,nat_slc*9,
     &       'munu4nbo: g_tvec',0,0,g_tvec)) 
     $       call errquit('hnd_giaox_zora: g_tvec', 0,GA_ERR)
        call ga_dgop(msg_efgs_col,dbl_mb(k_tvec),nat_slc*9,'+')
        call ga_put(g_tvec,1,1,1,nat_slc*9,dbl_mb(k_tvec),1)      
        call create_munu4nbo_shield(
     &                           rtdb,
     &                           g_tvec,
     &                           nat_slc,int_mb(k_AtNr),
     &                           basis,npol,nocc,nvirt,nmo) 
        if (.not. ga_destroy(g_tvec)) call errquit( ! destroy GA
     &      'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)   
        if (.not.ma_free_heap(l_tvec)) call
     &      errquit('hnd_giaox_zora: ma_free_heap l_tvec',0,MA_ERR)    
      endif ! ------------------------ hypfile-if++++ END
c ---- Destroy stored ga arrays ------ START
           if (.not. ga_destroy(g_d1_oo)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)  
           if (.not. ga_destroy(g_d1_ov)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR) 
       if (do_zora) then 
           if (.not. ga_destroy(ga_dia)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_para1)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_num)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)   
        if (.not. ga_destroy(ga_Fji)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR)    
        if (.not. ga_destroy(g_AtNr1)) call errquit(
     &    'hnd_giaox_zora: ga_destroy failed ',0, GA_ERR) 
       endif           
c ---- Destroy stored ga arrays ------ END
c
c     Clean up all remaining memory
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_dia',0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_para',0,MA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('hnd_giaox_zora: ga_destroy failed vectors',
     &             0,GA_ERR)
      enddo
c      if (.not.ga_destroy(vectors_scl(1))) call 
c     &    errquit('giao_aotomo: ga_destroy failed vectors',0,GA_ERR)
      if (.not.ma_pop_stack(l_zan)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_zan',0,MA_ERR)
      if (.not.ma_pop_stack(l_xyz)) call
     &    errquit('hnd_giaox_zora: ma_pop_stack failed k_xyz',0,MA_ERR)
       if (.not.ma_free_heap(l_AtNr)) call
     &     errquit('hnd_giaox_zora: ma_free_heap l_AtNr',0,MA_ERR)    
      call schwarz_tidy()
      call int_terminate()
      return
 7000 format(/,10x,'NMR shielding cannot be calculated for UHF',
     1      ' or ROHF wave functions: use RHF')
 9700 format(6x,'Atom: ',i4,2x,a2)
 9800 format(8x,'Diamagnetic',/,3(3F12.4,/))
 9801 format(8x,'Paramagnetic',/,3(3F12.4,/))
 9802 format(8x,'Total Shielding Tensor',/,3(3F12.4,/))
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9986 format(10x,'Principal Components and Axis System',/,10x,
     1       3(7x,i1,4x))
 9987 format(10x,' isotropic = ',f12.4,/,
     1       10x,'anisotropy = ',f12.4,/)
 9999 format(
     1 /,10x,41(1h-),/,
     2 10x,'Chemical Shielding Tensors (GIAO, in ppm)',/,
     3 10x,41(1h-),/)
      end

      subroutine get_par(par,       ! OUT : paramagnetic tensor
     &                   ga_para1,  ! IN  : paramagnetic tensor (GA) 
     &                   ga_h01_num,! IN  : H01 matrix
     &                   par_arr,   ! OUT : par tensor gauge,OO,OV
     &                   g_d1_oo,   ! IN  : Perturbed density matrix (OO part)
     &                   g_d1_ov,   ! IN  : Perturbed density matrix (OV part) 
     &                   g_dens,    ! IN  : e-density
     &                   npol,      ! IN  : nr. polarizations
     &                   Fukui,     ! OUT : Fukui term
     &                   basis,     ! IN  : basis handle
     &                   nbf,       ! IN  : nr. basis functions
     &                   pos,       ! IN  : Nuclear positions (x,y,z)
     &                   natoms_slc,! IN  : nr. atoms selected
     &                   natoms_tot,! IN  : nr. atoms total
     &                   geom,
     &                   rtdb,
     &                   oskel)     ! IN  : =.false.
c Purpose : Assemble NMR Chemical Shielding: paramagnetic tensor
c Author  : Fredy Aquino
c Date    : 03-03-11, 11-23-12 (adding Fukui term)
c Note.- Adaptation of hnd_giaox.F to hold 
c        - Nonrel(HF or DFT)   in unrestricted/restricted calculations.
c        - Relativistic (zora) in unrestricted/restricted calculations.
c
c Date    : 02-05-14: Fixed bug when selecting atoms
c           BEFORE :  call get_chi_centers_ga(g_R,basis,nbf,geom,natoms_slc)
c           AFTER  :  call get_chi_centers_ga(g_R,basis,nbf,geom,natoms_tot)
      implicit none
c
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
#include "zora.fh"
c
c      Global arrays variables defined in zora.fh
c      --> ga_dia,ga_para1,ga_Fji,ga_h01_num
      integer rtdb,geom
      integer basis
      integer nbf,natoms_slc,natoms_tot
      integer g_dens(3)          ! IN : electronic density
      double precision Fukui(3,3,natoms_slc) ! used in analyt shield
      integer npol  ! OUT: Used only in analytical shieldings
      integer alo(3), ahi(3), blo(3), bhi(3),
     &        dlo(3), dhi(3)
      integer ld(2),cbuf,cbuf1,cbuf2
      integer i,dims(3),chunk(3),
     &        g_R(3),ispin,a,b
      double precision pos(3*natoms_slc),val
      integer iatom,ix,iy,ixy,ind
      double precision val_oo,val_ov
      integer g_h01, ! local
     &        ga_h01_num
      integer g_d1_oo,g_d1_ov ! input: Perturbed density matrics OO and OV
                              ! OO, due to occupied-occupied MOs
                              ! OV, due to occupied-viritual MOs
      integer ga_para1,g_t1,g_t2,g_t3
      integer debug_get_par
      double precision par(*)
      logical oskel
      double precision ppm,par_arr(*)
      integer ind_tmn(2,3)
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      external get_chi_centers_ga
      data ppm     /26.62566914d+00/ 

c Warning: If I activate debugging (=1)
c          I could get trouble in final outputs
c          the content is affected if I attempt
c          to print out variables. (FA-11-24-11)

      debug_get_par=0

      do ixy = 1, 9*natoms_slc
         par(ixy) = 0.0d0  ! initialize
      enddo
      if (do_zora) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'Calc. par tensor-> zora'
c +++++ Initialize dbl_mb(k_para) with ga_para1 ++++ START
c ---- STORE: g_para1 --> dbl_mb(k_para)
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=natoms_slc
       ld(1)=3
       ld(2)=3
       call nga_get(ga_para1,alo,ahi,par(1),ld)
c +++++ Initialize dbl_mb(k_para) with g_para1 ++++ END
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
       ind=0
       do iatom = 1, natoms_slc
        do iy = 1, 3
         blo(3) = blo(3) + 1
         bhi(3) = bhi(3) + 1
         do ix = 1, 3
          alo(3) = ix
          ahi(3) = ix
          ixy = ixy + 1
          val_oo = nga_ddot_patch(g_d1_oo   ,'n',alo,ahi,
     &                            ga_h01_num,'n',blo,bhi)
          val_ov = nga_ddot_patch(g_d1_ov   ,'n',alo,ahi,
     &                            ga_h01_num,'n',blo,bhi)
         
          cbuf=9*(iatom-1)+3*(ix-1)+iy !-1 transpose
c ----- store in par_arr ------- START
          par_arr(ind+1)=par(cbuf) 
          par_arr(ind+2)=val_oo*ppm
          par_arr(ind+3)=val_ov*ppm
          par_arr(ind+4)=par(cbuf)+ 
     &                   val_oo*ppm+
     &                   val_ov*ppm

         if (debug_get_par.eq.1) then
               if (ga_nodeid().eq.0) then
                write(*,14) iatom,iy,ix,
     &                      alo(1),alo(2),alo(3),
     &                      ahi(1),ahi(2),ahi(3),
     &                      blo(1),blo(2),blo(3),
     &                      bhi(1),bhi(2),bhi(3),
     &                      par_arr(ind+1),
     &                      par_arr(ind+2),par_arr(ind+3),
     &                      par_arr(ind+4)
 14             format('(iatom,iy,ix)=(',i3,',',i3,',',i3,') ',
     &                  'alo=(',i3,',',i3,',',i3,') ',
     &                  'ahi=(',i3,',',i3,',',i3,') ',
     &                  'blo=(',i3,',',i3,',',i3,') ',
     &                  'bhi=(',i3,',',i3,',',i3,') ',
     &                  'para=(',f15.8,',',f15.8,',',
     &                   f15.8,',',f15.8,')')
               endif
         endif
c ----- store in par_arr ------- END      
          par(cbuf)=par(cbuf)+ 
     &              val_oo*ppm+val_ov*ppm
          ind=ind+4
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
      else ! Nonrel calc
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_t1))
     $ call errquit('get_par: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t2',0,0,g_t2))
     $ call errquit('get_par: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t3',0,0,g_t3))
     $ call errquit('get_par: g_t3',0,GA_ERR)
       call ga_zero(g_t3)
c -------- get R_nu ------------- START
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array R",chunk,g_R(i)))
     $    call errquit('get_par: g_R', 0,GA_ERR)
       enddo

       call get_chi_centers_ga(g_R,basis,nbf,geom,natoms_tot)

!       if (ga_nodeid().eq.0)
!     &   write(*,*) '----------g_R--------- START'
!       call ga_print(g_R(1))
!       call ga_print(g_R(2))
!       call ga_print(g_R(3))
!       if (ga_nodeid().eq.0)
!     &   write(*,*) '----------g_R--------- END'
c -------- get R_nu ------------- END
       if (ga_nodeid().eq.0)
     &  write(*,*) 'Calc. par tensor-> nonrel'
c ----- Calculate g_h01 --------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*natoms_slc
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,g_h01)) call 
     &    errquit('get_par: nga_create failed g_h01',0,GA_ERR)
      call ga_zero(g_h01)
      call int_giao_1ega(basis,basis,g_h01,'h01',pos(1),
     &                   natoms_slc,oskel)

         if (debug_get_par.eq.1) then
          if (ga_nodeid().eq.0)
     &     write(*,*) '------ g_h01 --------- START'
           call ga_print(g_h01)
          if (ga_nodeid().eq.0)
     &     write(*,*) '------ g_h01 --------- END'
         endif
c ----- Calculate g_h01 --------- END
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 0
       bhi(3) = 0
       dlo(1) = 1
       dhi(1) = nbf
       dlo(2) = 1
       dhi(2) = nbf
       ind=0

c       if (ga_nodeid().eq.0) 
c     &   write(*,*) '-------g_dens(1)----START'
c       call ga_print(g_dens(1))
c       if (ga_nodeid().eq.0) 
c     &   write(*,*) '-------g_dens(1)----END'
c       if (ga_nodeid().eq.0) 
c     &   write(*,*) '-------g_dens(2)----START'
c       call ga_print(g_dens(2))
c       if (ga_nodeid().eq.0) 
c     &   write(*,*) '-------g_dens(2)----END'

       do iatom = 1, natoms_slc
        do iy = 1, 3
         blo(3) = blo(3) + 1
         bhi(3) = bhi(3) + 1
         do ix = 1, 3
          alo(3) = ix
          ahi(3) = ix
          val_oo = nga_ddot_patch(g_d1_oo ,'n',alo,ahi,
     &                               g_h01,'n',blo,bhi)
          val_ov = nga_ddot_patch(g_d1_ov ,'n',alo,ahi,
     &                               g_h01,'n',blo,bhi)        
          cbuf=9*(iatom-1)+3*(ix-1)+iy ! transpose
c ----- Calc. Fukui term ------- START
          a=ind_tmn(1,ix)
          b=ind_tmn(2,ix)
          call nga_copy_patch('n',g_h01,blo,bhi,
     &                            g_t1 ,dlo,dhi)
          call ga_scale_cols(g_t1,g_R(b))          ! R_{nu,b} g_N
          call ga_scale_rows(g_t1,g_R(a))          ! R_{mu,a} [R_{nu,b} g_N] -> g_t1
          call nga_copy_patch('n',g_h01,blo,bhi,
     &                            g_t2 ,dlo,dhi)
          call ga_scale_cols(g_t2,g_R(a))          ! R_{nu,a} g_N
          call ga_scale_rows(g_t2,g_R(b))          ! R_{mu,b} [R_{nu,a} g_N] -> g_t2
          call ga_add(1.0d0,g_t1,-1.0d0,g_t2,g_t3) ! g_t3=-4c^2 <chi_mu|i/(2c)(R_mu x R_nu)_k h_t^{01}|chi_nu>

c          if (ga_nodeid().eq.0) then
c           write(*,2) iatom,ix,iy
c 2         format('-------g_t3(',i4,',',i4,',',i4,')----START')
c          endif
c           call ga_print(g_t3)
c          if (ga_nodeid().eq.0) then
c           write(*,3) iatom,ix,iy
c 3         format('-------g_t3(',i4,',',i4,',',i4,')----END')
c          endif

          val=0.0d0
          do ispin=1,npol
           val=val+nga_ddot_patch(g_dens(ispin),'n',dlo,dhi,
     &                                     g_t3,'n',dlo,dhi)         
          enddo ! end-loop-ispin
          Fukui(ix,iy,iatom)=val*ppm

c WARNING: If I attempt to print lines below the values for dia
c          are affected somehow.  This is a weird problem
c          Otherwise, everything is fine.
c         if (ga_nodeid().eq.0) then
c          write(*,1) ix,iy,iatom,Fukui(ix,iy,iatom)
c 1        format('Fukui(',i4,',',i4,',',i4,')=',f15.8)
c         endif

c ----- Calc. Fukui term ------- END
c ----- store in par_arr ------- START
          par_arr(ind+1)=Fukui(ix,iy,iatom)
          par_arr(ind+2)=val_oo*ppm
          par_arr(ind+3)=val_ov*ppm
          par_arr(ind+4)=par_arr(ind+1)+ 
     &                   val_oo*ppm+
     &                   val_ov*ppm

          if (debug_get_par.eq.1) then
               if (ga_nodeid().eq.0) then
                write(*,12) iatom,iy,ix,
     &                      alo(1),alo(2),alo(3),
     &                      ahi(1),ahi(2),ahi(3),
     &                      blo(1),blo(2),blo(3),
     &                      bhi(1),bhi(2),bhi(3),
     &                      par_arr(ind+1),
     &                      par_arr(ind+2),par_arr(ind+3),
     &                      par_arr(ind+4)
 12             format('(iatom,iy,ix)=(',i3,',',i3,',',i3,') ',
     &                  'alo=(',i3,',',i3,',',i3,') ',
     &                  'ahi=(',i3,',',i3,',',i3,') ',
     &                  'blo=(',i3,',',i3,',',i3,') ',
     &                  'bhi=(',i3,',',i3,',',i3,') ',
     &                  'para=(',f15.8,',',f15.8,',',
     &                   f15.8,',',f15.8,')')
               endif
          endif

c ----- store in par_arr ------- END  
          par(cbuf)=par_arr(ind+1)+
     &              val_oo*ppm+val_ov*ppm
          ind=ind+4
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
        if (.not.ga_destroy(g_h01)) call 
     &    errquit('get_par: ga_destroy failed g_h01',0,GA_ERR)
        if (.not.ga_destroy(g_t1)) call 
     &    errquit('get_par: ga_destroy failed g_t1',0,GA_ERR)
        if (.not.ga_destroy(g_t2)) call 
     &    errquit('get_par: ga_destroy failed g_t2',0,GA_ERR)
        if (.not.ga_destroy(g_t3)) call 
     &    errquit('get_par: ga_destroy failed g_t3',0,GA_ERR)
        do i=1,3
          if (.not.ga_destroy(g_R(i))) call 
     &    errquit('get_par: ga_destroy failed g_R',0,GA_ERR)
        enddo
      endif
c -------- symmetrize total par ------------ START
       do iatom = 1, natoms_slc
        do iy = 1, 3
         do ix = iy+1, 3     
          cbuf1=9*(iatom-1)+3*(ix-1)+iy ! transpose 
          cbuf2=9*(iatom-1)+3*(iy-1)+ix 
          val=(par(cbuf1)+par(cbuf2))/2.0d0
          par(cbuf1)=val
          par(cbuf2)=val 
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
c -------- symmetrize total par ------------ END
      return
      end
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ============== get_par_JK() ===========================START
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_par_JK(
     &                   par,         ! OUT : paramagnetic tensor
     &                   ga_para1,    ! IN  : paramagnetic tensor (GA) 
     &                   ga_h01_num,  ! IN  : H01 matrix
     &                   par_arr,     ! OUT : par tensor gauge,OO,OV
     &                   g_d1_oo,     ! IN  : Perturbed density matrix (OO part)
     &                   g_d1_ov,     ! IN  : Perturbed density matrix (OV part) 
     &                   g_d1_ov_Coul,! IN  : Perturbed density matrix (OV part-Coul only) 
     &                   g_d1_ov_Exch,! IN  : Perturbed density matrix (OV part-Exch only) 
     &                   g_d1_ov_noJK,! IN  : Perturbed density matrix (OV part-all -(Coul,Exch)) 
     &                   g_d1_ov_1e,  ! IN  : Perturbed density matrix (OV part-1e   contrib)
     &                   g_d1_ov_eSji,! IN  : Perturbed density matrix (OV part-eSji contrib)
     &                   g_dens,      ! IN  : e-density
     &                   npol,        ! IN  : nr. polarizations
     &                   Fukui,       ! OUT : Fukui term
     &                   basis,       ! IN  : basis handle
     &                   nbf,         ! IN  : nr. basis functions
     &                   pos,         ! IN  : Nuclear positions (x,y,z)
     &                   natoms_slc,  ! IN  : nr. atoms selected
     &                   natoms_tot,  ! IN  : nr. atoms total
     &                   geom,
     &                   rtdb,
     &                   oskel)
c Purpose : Assemble NMR Chemical Shielding: paramagnetic tensor
c Author  : Fredy Aquino
c Date    : 03-03-11, 11-23-12 (adding Fukui term)
c Note.- Adaptation of hnd_giaox.F to hold 
c        - Nonrel(HF or DFT)   in unrestricted/restricted calculations.
c        - Relativistic (zora) in unrestricted/restricted calculations.
c
      implicit none
c
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
#include "zora.fh"
c
      integer rtdb,geom
      integer basis
      integer alo(3), ahi(3), blo(3), bhi(3),
     &        dlo(3),dhi(3)
      integer ld(2),cbuf,cbuf1,cbuf2
      integer nbf,natoms_slc,natoms_tot
      integer g_dens(3)          ! IN : electronic density
      double precision Fukui(3,3,natoms_slc) ! used in analyt shield
      integer npol  ! OUT: Used only in analytical shieldings
      integer i,dims(3),chunk(3),
     &        g_R(3),ispin,a,b
      double precision pos(3*natoms_slc),val
      integer iatom,ix,iy,ixy,ind
      double precision val_oo,val_ov
      double precision par_arr(*) 
      double precision val_ov_Coul,val_ov_Exch,val_ov_noJK,
     &                 val_ov_1e,val_ov_eSji
      integer g_h01, ! local
     &        ga_h01_num
      integer g_d1_oo,g_d1_ov ! input: Perturbed density matrics OO and OV
                              ! OO, due to occupied-occupied MOs
                              ! OV, due to occupied-viritual MOs
      integer g_d1_ov_Coul,g_d1_ov_Exch,g_d1_ov_noJK,
     &        g_d1_ov_1e,g_d1_ov_eSji
      integer ga_para1,g_t1,g_t2,g_t3
      integer debug_get_par
      double precision par(*)
      logical oskel
      integer ind_tmn(2,3)
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      external get_chi_centers_ga
      double precision ppm
      data ppm     /26.62566914d+00/ 

      debug_get_par=0

      do ixy = 1, 9*natoms_slc
         par(ixy) = 0.0d0  ! initialize
      enddo
      if (do_zora) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'Calc. par tensor-> zora'
c +++++ Initialize dbl_mb(k_para) with ga_para1 ++++ START
c ---- STORE: g_para1 --> dbl_mb(k_para)
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=natoms_slc
       ld(1)=3
       ld(2)=3
       call nga_get(ga_para1,alo,ahi,par(1),ld)
c +++++ Initialize dbl_mb(k_para) with g_para1 ++++ END
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
       ind=0
       do iatom = 1, natoms_slc
        do iy = 1, 3
         blo(3) = blo(3) + 1
         bhi(3) = bhi(3) + 1
         do ix = 1, 3
          alo(3) = ix
          ahi(3) = ix
          ixy = ixy + 1
c         val = nga_ddot_patch(g_d1 ,'n',alo,ahi,
c     &                        ga_h01_num,'n',blo,bhi)
          val_oo      = nga_ddot_patch(g_d1_oo   ,'n',alo,ahi,
     &                         ga_h01_num,'n',blo,bhi)
          val_ov      = nga_ddot_patch(g_d1_ov   ,'n',alo,ahi,
     &                         ga_h01_num,'n',blo,bhi)
          val_ov_Coul = nga_ddot_patch(g_d1_ov_Coul,'n',alo,ahi,
     &                          ga_h01_num,'n',blo,bhi)
          val_ov_Exch = nga_ddot_patch(g_d1_ov_Exch,'n',alo,ahi,
     &                          ga_h01_num,'n',blo,bhi)
          val_ov_noJK = nga_ddot_patch(g_d1_ov_noJK,'n',alo,ahi,
     &                          ga_h01_num,'n',blo,bhi)  
          val_ov_1e   = nga_ddot_patch(g_d1_ov_1e,'n',alo,ahi,
     &                          ga_h01_num,'n',blo,bhi)   
          val_ov_eSji = nga_ddot_patch(g_d1_ov_eSji,'n',alo,ahi,
     &                          ga_h01_num,'n',blo,bhi)      
         
          cbuf=9*(iatom-1)+3*(ix-1)+iy !-1 transpose
c ----- store in par_arr ------- START
          par_arr(ind+1)=par(cbuf)
          par_arr(ind+2)=val_oo*ppm
          par_arr(ind+3)=val_ov*ppm
          par_arr(ind+4)=par(cbuf)+
     &                   val_oo*ppm+
     &                   val_ov*ppm
          par_arr(ind+5)=val_ov_Coul*ppm
          par_arr(ind+6)=val_ov_Exch*ppm
          par_arr(ind+7)=val_ov_noJK*ppm
          par_arr(ind+8)=val_ov_1e*ppm
          par_arr(ind+9)=val_ov_eSji*ppm
c ----- store in par_arr ------- END      
          par(cbuf)=par(cbuf)+ 
     &              val_oo*ppm+val_ov*ppm
          ind=ind+9
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
      else ! Nonrel calc
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_t1))
     $ call errquit('get_par_JK: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t2',0,0,g_t2))
     $ call errquit('get_par_JK: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t3',0,0,g_t3))
     $ call errquit('get_par_JK: g_t3',0,GA_ERR)
       call ga_zero(g_t3)
c -------- get R_nu ------------- START
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array R",chunk,g_R(i)))
     $    call errquit('get_par_JK: g_R', 0,GA_ERR)
       enddo
       call get_chi_centers_ga(g_R,basis,nbf,geom,natoms_tot)
c -------- get R_nu ------------- END
       if (ga_nodeid().eq.0)
     &  write(*,*) 'Calc. par tensor-> nonrel'
c ----- Calculate g_h01 --------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*natoms_slc
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,g_h01)) call 
     &    errquit('get_par_JK: nga_create failed g_h01',0,GA_ERR)
      call ga_zero(g_h01)
      call int_giao_1ega(basis,basis,g_h01,'h01',pos(1),
     &                   natoms_slc,oskel)
         if (debug_get_par.eq.1) then
          if (ga_nodeid().eq.0)
     &     write(*,*) '------ g_h01 --------- START'
           call ga_print(g_h01)
          if (ga_nodeid().eq.0)
     &     write(*,*) '------ g_h01 --------- END'
         endif
c ----- Calculate g_h01 --------- END
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 0
       bhi(3) = 0
       ind=0
       do iatom = 1, natoms_slc
        do iy = 1, 3
         blo(3) = blo(3) + 1
         bhi(3) = bhi(3) + 1
         do ix = 1, 3
          alo(3) = ix
          ahi(3) = ix
          val_oo = nga_ddot_patch(g_d1_oo ,'n',alo,ahi,
     &                               g_h01,'n',blo,bhi)
          val_ov = nga_ddot_patch(g_d1_ov ,'n',alo,ahi,
     &                               g_h01,'n',blo,bhi)        
          cbuf=9*(iatom-1)+3*(ix-1)+iy ! transpose
c ----- Calc. Fukui term ------- START
         a=ind_tmn(1,ix)
         b=ind_tmn(2,ix)
         call nga_copy_patch('n',g_h01,blo,bhi,
     &                           g_t1 ,dlo,dhi)
         call ga_scale_cols(g_t1,g_R(b))          ! R_{nu,b} g_N
         call ga_scale_rows(g_t1,g_R(a))          ! R_{mu,a} [R_{nu,b} g_N] -> g_t1
         call nga_copy_patch('n',g_h01,blo,bhi,
     &                           g_t2 ,dlo,dhi)
         call ga_scale_cols(g_t2,g_R(a))          ! R_{nu,b} g_N
         call ga_scale_rows(g_t2,g_R(b))          ! R_{mu,a} [R_{nu,b} g_N] -> g_t2
         call ga_add(1.0d0,g_t1,-1.0d0,g_t2,g_t3) ! g_t3=-4c^2 <chi_mu|i/(2c)(R_mu x R_nu)_k h_t^{01}|chi_nu>
         val=0.0d0
         do ispin=1,npol
          val=val+nga_ddot_patch(g_dens(ispin),'n',dlo,dhi,
     &                                    g_t3,'n',dlo,dhi)         
         enddo ! end-loop-ispin

         Fukui(ix,iy,iatom)=val*ppm

c         write(*,1) ix,iy,iatom,Fukui(ix,iy,iatom)
c 1       format('Fukui(',i4,',',i4,',',i4,')=',f15.8)

c ----- Calc. Fukui term ------- END
c ----- store in par_arr ------- START
          par_arr(ind+1)=Fukui(ix,iy,iatom)
          par_arr(ind+2)=val_oo*ppm
          par_arr(ind+3)=val_ov*ppm
          par_arr(ind+4)=par_arr(ind+1)+
     &                   val_oo*ppm+
     &                   val_ov*ppm
          if (debug_get_par.eq.1) then
               if (ga_nodeid().eq.0) then
                write(*,12) iatom,iy,ix,
     &                      alo(1),alo(2),alo(3),
     &                      ahi(1),ahi(2),ahi(3),
     &                      blo(1),blo(2),blo(3),
     &                      bhi(1),bhi(2),bhi(3),
     &                      par_arr(ind+1),
     &                      par_arr(ind+2),par_arr(ind+3),
     &                      par_arr(ind+4)
 12             format('(iatom,iy,ix)=(',i3,',',i3,',',i3,') ',
     &                  'alo=(',i3,',',i3,',',i3,') ',
     &                  'ahi=(',i3,',',i3,',',i3,') ',
     &                  'blo=(',i3,',',i3,',',i3,') ',
     &                  'bhi=(',i3,',',i3,',',i3,') ',
     &                  'para=(',f15.8,',',f15.8,',',
     &                   f15.8,',',f15.8,')')
               endif
          endif
c ----- store in par_arr ------- END  
          par(cbuf)=par_arr(ind+1)+
     &              val_oo*ppm+val_ov*ppm
          ind=ind+4
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
        if (.not.ga_destroy(g_h01)) call 
     &    errquit('get_par_JK: ga_destroy failed g_h01',0,GA_ERR)
        if (.not.ga_destroy(g_t1)) call 
     &    errquit('get_par_JK: ga_destroy failed g_t1',0,GA_ERR)
        if (.not.ga_destroy(g_t2)) call 
     &    errquit('get_par_JK: ga_destroy failed g_t2',0,GA_ERR)
        if (.not.ga_destroy(g_t3)) call 
     &    errquit('get_par_JK: ga_destroy failed g_t3',0,GA_ERR)
        do i=1,3
          if (.not.ga_destroy(g_R(i))) call 
     &    errquit('get_par_JK: ga_destroy failed g_R',0,GA_ERR)
        enddo
      endif
c -------- symmetrize total par ------------ START
       do iatom = 1, natoms_slc
        do iy = 1, 3
         do ix = iy+1, 3     
          cbuf1=9*(iatom-1)+3*(ix-1)+iy ! transpose 
          cbuf2=9*(iatom-1)+3*(iy-1)+ix 
          val=(par(cbuf1)+par(cbuf2))/2.0d0
          par(cbuf1)=val
          par(cbuf2)=val 
         enddo ! end-loop-ix
        enddo ! end-loop-iy
       enddo ! end-loop-iatom
c -------- symmetrize total par ------------ END
      return
      end
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ============== get_par_JK() =============================END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_dia(dia,    ! OUT : paramagnetic tensor
     &                   ga_dia, ! IN  : dia tensor
     &                   basis,  ! IN  : basis handle
     &                   g_dens, ! IN  : e-density
     &                   nbf,    ! IN  : nr. basis functions
     &                   npol,   ! IN  : nr. of polarizations
     &                   Fukui,  ! IN  : Fukui term
     &                   pos,    ! IN  : Nuclear positions (x,y,z)
     &                   natoms, ! IN  : nr. atoms
     &                   rtdb,
     &                   oskel)

c Purpose : Assemble NMR Chemical Shielding: diamagnetic tensor
c Author  : Fredy Aquino
c Date    : 03-03-11, 11-23-12 (adding Fukui term)
c Note.- Adaptation of hnd_giaox.F to hold 
c        - Nonrel(HF or DFT)   in unrestricted/restricted calculations.
c        - Relativistic (zora) in unrestricted/restricted calculations.
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "stdio.fh"
#include "zora.fh"
c
c      Global arrays variables defined in zora.fh
c      --> ga_dia,ga_para1,ga_Fji,ga_h01_num
      integer rtdb
      integer basis     ! [input] basis handle
      integer g_dens(3) ! IN: electronic density
      integer ga_dia
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer ld(2),cbuf,
     &        ix1,iy1
      integer nbf,natoms,npol,ispin
      double precision ac,pos(3*natoms)
      double precision Fukui(3,3,natoms) ! used in analyt shield
      integer iatom,ix,iy,ixy,ind1,ind2
      integer g_h11 ! local ga
      double precision dia(*)
      logical oskel
      double precision ppm,val
      data ppm     /26.62566914d+00/ 
      do ixy = 1, 9*natoms
         dia(ixy) = 0.0d0  ! initialize
      enddo
      if (do_zora) then  ! get-dia---- START
c ---- STORE: ga_dia --> dbl_mb(k_dia)
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=natoms
       ld(1)=3
       ld(2)=3 
       call nga_get(ga_dia,alo,ahi,dia(1),ld)
      else
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 9*natoms
       if (.not.nga_create(MT_DBL,3,ahi,'H11 matrix',alo,g_h11)) call 
     &    errquit('get_dia: nga_create failed g_h11 all',0,GA_ERR)
       call ga_zero(g_h11)
       call int_giao_1ega(basis,basis,g_h11,'h11 para',
     &                    pos(1),natoms,oskel)
       blo(1) = 1
       bhi(1) = nbf
       blo(2) = 1
       bhi(2) = nbf
       blo(3) = 0
       bhi(3) = 0
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       ixy = 0
       do iatom = 1, natoms
         ix1=1
         iy1=1
         do ix = 1, 9
            if (ix1.gt.3) then
             ix1=1
             iy1=iy1+1
            endif
            ixy = ixy + 1
            alo(3) = ixy
            ahi(3) = ixy
            do ispin=1,npol
             val=nga_ddot_patch(g_dens(ispin),'n',blo,bhi,
     &                                  g_h11,'n',alo,ahi)
             dia(ixy)=dia(ixy) + val * ppm            
            enddo ! end-loop-ispin
            dia(ixy)=dia(ixy)-Fukui(ix1,iy1,iatom)
c Note.- Printing lines below does not affect output
c        It is ok to print it out but I will comment to
c        reduce printouts.

c             if (ga_nodeid().eq.0) then
c              write(*,1) ix1,iy1,iatom,Fukui(ix1,iy1,iatom)
c  1           format('Fukui(ix1,iy1,iatom)(',i4,',',i4,',',i4,')=',
c     &               f15.8)
c             endif

            ix1=ix1+1
         enddo ! end-loop-ix
       enddo ! end-loop-atoms
c       
c     s(dia)xy  = s(dia)xy + Sum(n,l) D0(n,l) * H11(dia)xy(n,l)

       call ga_zero(g_h11)
       call int_giao_1ega(basis,basis,g_h11,'h11 dia',
     &                    pos(1),natoms,oskel)
       alo(1) = 1
       ahi(1) = nbf
       alo(2) = 1
       ahi(2) = nbf
       ixy=0
       do iatom = 1, natoms
         ix1=1
         iy1=1
         do ix = 1, 9
            if (ix1.gt.3) then
             ix1=1
             iy1=iy1+1
            endif
            ixy = ixy + 1
            alo(3) = ixy
            ahi(3) = ixy
            ac=0.0d0
            do ispin=1,npol
             val=nga_ddot_patch(g_dens(ispin),'n',blo,bhi,
     &                                  g_h11,'n',alo,ahi)
             ac=ac + val * ppm 
            enddo ! end-loop-ispin
            dia(ixy)=dia(ixy)+ac 
c
c WARNING: If I try to print ac, uncommenting lines below
c          the dia(ixy) values are affected somehow. (FA-11-24-12)
c          This is a weird problem, otherwise everything comes
c          fine.
c             if (ga_nodeid().eq.0) then
c              write(*,2) ix1,iy1,iatom,ac
c 2            format('h11-dia(',i4,',',i4,',',i4,')=',
c     &               f15.8)
c             endif

            ix1=ix1+1
         enddo
       enddo
       if (.not.ga_destroy(g_h11)) call 
     &    errquit('get_dia: ga_destroy failed g_h11',0,GA_ERR)
      endif ! get-dia---- END
c ---------- symmetrize dia ----------- START
       do iatom=1,natoms
        do ix=1,3
         do iy=ix+1,3
           ind1=9*(iatom-1)+3*(iy-1)+ix
           ind2=9*(iatom-1)+3*(ix-1)+iy
           val=(dia(ind1)+dia(ind2))/2.0d0
           dia(ind1)=val
           dia(ind2)=val        
         enddo ! end-loop-iy 
        enddo ! end-loop ix
       enddo ! end-loop-iat
c ---------- symmetrize dia ----------- END
      return
      end
c -------------- for shieldings NMLO analysis ----------------- START
      subroutine create_munu4nbo_shield(
     &             rtdb,    ! in: rtdb handle
     &             g_tvec,  ! in: eigenvectors or T diagonalizing matrix
     &             nat,     ! in: nr. atoms
     &             atmnr,   ! in: selected atoms
     &             basis,   ! in: basis handle
     &             npol,    ! in: nr. polarizaitons
     &             nocc,    ! in: nr. occ   nocc(i) i=1,npol
     &             nvirt,   ! in: nr. virt nvirt(i) i=1,npol
     &             nmo)     ! in: nr. MO
c
      implicit none
c
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
c
c FA: Revised on 06-22-11
c ------Main outputs -------- START
      integer g_munu_rot,         ! hyp-FCSD  dia  part
     &        g_munu_rot1,
     &        g_munu_rot2,
     &        g_acc2  ! hyp-PSOSO para part
c ------Main outputs -------- END
      integer rtdb,basis,nat
      integer g_dens,g_tvec,atmnr(nat)
      integer g_munuCSdia,g_munuCSpar1,g_munuCSHpar,
     &        g_tnp,g_acc,vectors(2),
     &        g_tnp1,g_acc1,g_acc3   
      integer vectors_scl(2),ispin,iat1  
      double precision ac_val,val1,sign    
      integer i,j,k,m,n,ndir,ndir1
      integer jlo,jhi,s,nbf,nmo,nsize,nsize1,nsize2
      integer npol,npol_munu,ntot
      integer ind,nlst,count,nocc(npol),nvirt(npol)
      integer Natoms_munu,atmnr_munu(nat)
      integer Ndir_munu
c     Ndir_munu, Nr. of directions stored
c                =3  xx yy zz
      double precision coeff,fact,para_rot(9)
      double precision tmn(2),chcdata(3)
      integer jlo1,jhi1,jlo2,jhi2
      integer g_dens1,g_c1,
     &        g_p10,g_munuCSHpar2d,g_munuCSpar12d
      integer iind(2),jind(2),icalczora,type_NMR,iat
      integer alo(3),ahi(3),elo(3),ehi(3),flo(3),fhi(3)
      logical dft_zoraShield_NLMOAnalysis_read ! for read-nlmo-mat
      character*255 zorafilename              ! for read-nlmo-mat
      integer arr_ind(6,2)
       data ((arr_ind(j,i),i=1,2),j=1,6)
     &  /1,1,2,2,3,3,1,2,1,3,2,3/
      external dft_zoraShield_NLMOAnalysis_read,get_P10_rot,
     &         fill_munuPSOSO_1,get_par_gshift_rot,
     &         get2dmat,wshldfile
c     --> To store ONLY munu principal components xx,yy,zz 
c     g_munuCSdia    is created in dft_zora_EPR.F
c     size(g_munuCSdia)=nlst*ndir (linear vector)
c     g_dens, spin density matrix
c     nbf x nbf elements (bidimensional matrix)
c     Legend:
c     ndir=6 = xx, yy, zz, xy, xz, yz
c     nbf, Nr of basis functions
c     nlst=nbf*(nbf+1)/2

      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('munu: bas_numbf failed',555, BASIS_ERR)
      Natoms_munu=nat
      do i=1,Natoms_munu
       atmnr_munu(i)=atmnr(i)
      enddo
      npol_munu=npol
      Ndir_munu=3
      nlst=nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
c ++++++ Read NLMO matrices +++++++++ START
      ndir =6
      ndir1=3
      call util_file_name(lbl_nlmoshield,.false.,.false.,
     &                    zorafilename)    
      icalczora = 0  ! initialize the flag
      if (.not.dft_zoraShield_NLMOAnalysis_read(
     &       zorafilename, ! in : filename
     &                nbf, ! in : nr basis functions
     &               ndir, ! in : nr of directions: 6 = xx yy zz xy xz yz
     &              ndir1, ! in : nr of directions: 3 = x y z
     &                nat, ! in : list of selected atoms
     &               nocc, ! in : nocc(i) i=1,2 nr. occupations in alpha and beta   
     &               npol, ! in: nr polarizations
     &        g_munuCSdia, ! out: munu matrix of dia
     &       g_munuCSpar1, ! out: munu matrix of par 1st term
     &       g_munuCSHpar, ! out: munu matrix of H10
     &            vectors, ! out: MOs
     &               g_c1, ! out: perturbed MO
     &             g_dens)) icalczora=1 
 
       goto 10
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuCSdia -- START'
        call ga_print(g_munuCSdia)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_munuCSdia -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &   'AFT-reading-HYP-NLMO matrices: g_munuCSHpar -- START'
        call ga_print(g_munuCSHpar)
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &  'AFT-reading-HYP-NLMO matrices: g_munuCSHpar -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &   'AFT-reading-HYP-NLMO matrices: g_munuCSpar1 -- START'
        call ga_print(g_munuCSpar1)
      if (ga_nodeid().eq.0)
     & write(*,*) 
     &  'AFT-reading-HYP-NLMO matrices: g_munuCSpar1 -- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_dens ----- START'
        call ga_print(g_dens)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_dens ----- END'
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- START'
        call ga_print(g_c1)
      if (ga_nodeid().eq.0)
     & write(*,*) 'AFT-reading-HYP-NLMO matrices: g_c1 ----- END'
 10   continue
c ++++++ Read NLMO matrices +++++++++ END
      call get_unique_elmat(g_dens,g_dens1,nlst,nbf)   ! out: g_dens1
      ndir =6 ! Nr. of directions: xx,yy,zz,xy,xz,yz
      ndir1=3 ! Nr. of directions: x,y,z
      nsize =nbf*(nbf+1)/2 ! size of xx,yy,zz,xy,xz,yz chunk
      nsize1=nsize*ndir    ! size of whole munu per atom
      nsize2=nsize*ndir1   ! size of whole munu per atom
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
         ahi(3) = 3*nat

         if (.not.nga_create(MT_DBL,3,ahi,'g_acc2 matrix',
     &       alo,g_acc2)) call 
     &        errquit('g_acc2: nga_create failed g_acc2',0,GA_ERR)
         call ga_zero(g_acc2)
        if (.not. ga_create(mt_dbl,1,nlst*3*nat,
     &                       'munu4nbo: g_munu_rot',0,0,g_munu_rot)) 
     $    call errquit('munu4nbo: g_munu_rot', 0,GA_ERR)
        if (.not. ga_create(mt_dbl,1,nlst*3*nat,
     &                       'munu4nbo: g_munu_rot2',0,0,g_munu_rot2)) 
     $    call errquit('munu4nbo: g_munu_rot2', 0,GA_ERR)
        if (.not. ga_create(mt_dbl,1,nlst*3*nat,
     &                       'munu4nbo: g_munu_rot1',0,0,g_munu_rot1)) 
     $    call errquit('munu4nbo: g_munu_rot1', 0,GA_ERR)
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nbf
       ahi(3) = 3
      if (.not.nga_create(mt_dbl,3,ahi,'g_munuCSHpar2d matrix',
     &    alo,g_munuCSHpar2d)) call 
     &    errquit('munu4nbo: nga_create failed g_munuCSHpar2d',
     &            0,GA_ERR)
        call ga_zero(g_munuCSHpar2d)
      if (.not.nga_create(mt_dbl,3,ahi,'g_munuCSpar12d matrix',
     &    alo,g_munuCSpar12d)) call 
     &    errquit('munu4nbo: nga_create failed g_munuCSpar12d',
     &            0,GA_ERR)
        call ga_zero(g_munuCSpar12d)
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &              ' NW-Shieldings: Summary C+HC data [ppm] ',
     &              'ooooooooooooooo START'
      do iat1=1,nat
        call ga_zero(g_munuCSpar12d)
        iat=atmnr(iat1)
        do n=1,3  ! xx,yy,zz
         m=n ! For principal components ONLY
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
c        call reorder_munu(g_munu_rot,nat,nlst,nbf,Ndir_munu)
c --------------------------------------------------------------
         call ga_zero(g_acc)
         call ga_zero(g_acc3)
         do s=1,6
c ------- get coeff() --- START
          iind(1)=1
          iind(2)=1
          jind(1)=9*(iat1-1)+3*(arr_ind(s,1)-1)+m
          jind(2)=9*(iat1-1)+3*(arr_ind(s,2)-1)+n   
          call ga_gather(g_tvec,tmn,iind,jind,2)
          fact=1.0d0
          if (s.gt.3) fact=2.0d0
          coeff=fact*tmn(1)*tmn(2)      
c ------- get coeff() --- END
c Note.- g_munuFCSD will be the (hyp-diag)_uv matrix
           jlo=nsize1*(iat1-1)+nsize*(s-1)+1
           jhi=jlo+nsize-1
           call ga_copy_patch('n',g_munuCSdia,1,1,jlo,jhi,
     &                            g_tnp      ,1,1,1  ,nsize)
           call ga_add(1.0d0,g_acc,coeff,g_tnp,g_acc)  

           call ga_copy_patch('n',g_munuCSpar1,1,1,jlo,jhi,
     &                            g_tnp       ,1,1,1  ,nsize)
           call ga_add(1.0d0,g_acc3,coeff,g_tnp,g_acc3)    
         enddo ! end-loop-s

c Note.- g_acc = \widetilde{H}_{mu nu}^{(m,m)}
c        it is the rotated munu matrix using:  T^t H T
         call ga_zero(g_acc1)
         do s=1,3
c ------- get coeff() --- START
          iind(1)=1
          jind(1)=9*(iat1-1)+3*(s-1)+m
          call ga_gather(g_tvec,tmn,iind,jind,1)
c          if (ga_nodeid().eq.0) then
c           write(*,1) m,s,tmn(1)
c 1         format('(m,s,tvec)=(',i5,',',i5,',',f15.8,')')
c          endif
          coeff=tmn(1) 
c ------- get coeff() --- END
c Note.- g_munuCSHpar will be the (g-shift-para)_uv matrix
          jlo=nsize2*(iat1-1)+nsize*(s-1)+1
          jhi=jlo+nsize-1
          call ga_copy_patch('n',g_munuCSHpar,1,1,jlo,jhi,
     &                           g_tnp1      ,1,1,1  ,nsize)
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
          flo(3) = 3*(iat1-1)+m
          fhi(3) = 3*(iat1-1)+m
          call nga_add_patch(1.0d0,g_acc2,flo,fhi,
     &                       coeff,g_c1  ,elo,ehi,
     &                             g_acc2,flo,fhi)
c ----- Calculate rotated perturbed MO: g_acc2 ----- END 
         enddo ! end-loop-s
c Note: g_acc1 = \widetilde{H}_{mu nu}^{(m)}  m=x,y,z
c       it is the rotated munu matrix using: T H
c ====== Store final munu matrices === START
         jlo2=nlst*Ndir_munu*(iat1-1)+
     &        nlst*(n-1)+1
         jhi2=jlo2+nlst-1
         call ga_copy_patch('n',g_acc     ,1,1,   1,nlst,
     &                          g_munu_rot,1,1,jlo2,jhi2)   
         call ga_copy_patch('n',g_acc3    ,1,1,   1,nlst,
     &                         g_munu_rot2,1,1,jlo2,jhi2)         
         call ga_copy_patch('n',g_acc1    ,1,1,   1,nlst,
     &                         g_munu_rot1,1,1,jlo2,jhi2)     
c ====== Store final munu matrices === END
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== START
c ==== sum (g_acc .* g_dens1 + Nuclear CONTRIB) 
c      = TOTAL Shieldings diagonalized   
         jlo1=1+nbf
         jhi1=nsize
         call ga_scale_patch(g_acc,1,1,jlo1,jhi1,2.0d0)
         chcdata(m)=ga_ddot(g_acc,g_dens1)
c ++++++++++++++++++CHECK++++ DIAGONALIZATION ==== END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++
        enddo ! end-loop-n
        if (ga_nodeid().eq.0) then
          write(*,23) iat,
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
        type_NMR=1 ! =1,2,3=shieldings,hyperfine,gshift
        call get_P10_rot(
     &          g_p10,            ! out: Perturbed density matrix (munu nbf x nbf x 3 square matrix)
     &          type_NMR,         !  in: =1,2,3=shieldings,hyperfine,gshift
     &          g_acc2,           !  in: rotated perturbed MO vector
     &          vectors,g_CiFull, !  in: to build zora scaled MO vector 
     &          iat1,             !  in: index for selected atom nr =1,nat
     &          nbf,nmo,npol,nocc,nvirt,
     &          do_zora,do_NonRel,not_zora_scale,rtdb) 
      call fill_munuPSOSO_1(   ! g_munuCSHpar --> g_munuCSHpar2d
     &        g_munu_rot1,     ! in: array with unique elements
     &        g_munuCSHpar2d,  !out: nbf x nbf x 3 munu matrix for ith atom
     &        iat1,            ! in: = 1,2,...,nat
     &        2,               ! in: type_symm = 1 symm  = 2 antisymm
     &        nbf) 
      call fill_munuPSOSO_1(   ! g_munuCSpar1 --> g_munuCSHpar2d
     &        g_munu_rot2,     ! in: array with unique elements
     &        g_munuCSpar12d,  !out: nbf x nbf x 3 munu matrix for ith atom
     &        iat1,            ! in: = 1,2,...,nat
     &        1,               ! in: type_symm = 1 symm  = 2 antisymm
     &        nbf) 
c +++++++++ NOW: do ddot product to get diagonalized tensor +++ START
      call get_par_gshift_rot(
     &                   g_dens,         ! IN : spin-density
     &                   g_munuCSpar12d, ! IN : par 1st term
     &                   g_munuCSHpar2d, ! IN : h01 matrix
     &                   g_p10,          ! IN : Perturbed density matrix
     &                   basis,nbf,iat,rtdb)
c +++++++++ NOW: do ddot product to get diagonalized tensor +++ END
c ++++++++++ CHECK diagonalization in para hyperfine ++++++++ END
      if (.not. ga_destroy(g_p10)) call errquit(
     &  'create_munu4nbo_shield: ga_destroy failed g_p10',0, GA_ERR)
      enddo ! end-loop-iat
        if (ga_nodeid().eq.0)
     &   write(*,*) 'CHCooooooooooooo',
     &              ' NW-Shieldings: Summary C+HC data [ppm] ',
     &              'ooooooooooooooo END'
c --> Main outputs: g_acc2     ,rotated perturbed MO nbf*ntot*ndir*nat
c                               ndir=1,2,3=x,y,z
c                               nbf, nr of basis functions
c                               ntot=nocc(1)+nocc(2)
c                               nat, nr of selected atoms
c                   g_munu_rot1,rotated perturbed AO matrix
c                               storing only diag + off-diag elements
c                               Reminder: this comes from an antisymmetrix matrix
c                                         in case we want to pull back the 2d munu-matrix
c                   g_munu_rot, rotated AO matrix for dia part
c                               storing only diag + off-diag elements
c                               Reminder: this comes from a  symmetric matrix
c                                         in case we want to pull back the 2d munu-matrix
      call reorder_munu(g_munu_rot ,nat,nlst,nbf,Ndir_munu) ! reoder-munu matrix
      call reorder_munu(g_munu_rot1,nat,nlst,nbf,Ndir_munu) ! reoder-munu matrix
      call reorder_munu(g_munu_rot2,nat,nlst,nbf,Ndir_munu) ! reoder-munu matrix
c ------ destroy unnecessary GAs 
      if (.not. ga_destroy(g_munuCSdia)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuCSpar1)) call errquit(
     &  'create_munu4nbo-1: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuCSHpar)) call errquit(
     &  'create_munu4nbo-2: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuCSHpar2d)) call errquit(
     &  'create_munu4nbo-7a: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munuCSpar12d)) call errquit(
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
      if (.not. ga_destroy(g_dens)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_dens1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_c1)) call errquit(
     &  'create_munu4nbo: ga_destroy failed ',0, GA_ERR)
      do i=1,npol
 911  if (.not.ga_destroy(vectors(i))) call 
     &    errquit('create_munu4nbo: ga_destroy failed vectors',
     &             0,GA_ERR)
      enddo
      call wshldfile(rtdb,
     &               g_munu_rot,  ! dia term 
     &               g_munu_rot2, ! 1st term in para      g_munuEPRpar1
     &               g_acc2,      ! perturbed MO vector x,y,z
     &               g_munu_rot1, ! perturbed AO operator g_munuEPRHpar x,y,z
     &               nlst,npol_munu,       
     &               Ndir_munu,   ! OUTPUT: used in wgshiftfile(rtdb)
     &               Natoms_munu,
     &               atmnr_munu)
      if (.not. ga_destroy(g_munu_rot)) call errquit(  
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munu_rot1)) call errquit( 
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_munu_rot2)) call errquit( 
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_acc2)) call errquit(     
     &  'wshldfile: ga_destroy failed ',0, GA_ERR)
      return
      end  
c -------------- for shieldings NMLO analysis ----------------- END
