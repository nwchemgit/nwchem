      subroutine hnd_giaox_ZORA(rtdb,basis,geom)
c $Id: hnd_giaox.F,v 1.10 2007-10-15 23:31:27 bert Exp $
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
#include "zora.fh"  ! for num NMRCS-Z4
c      Global arrays variables defined in zora.fh
c      --> ga_dia,ga_para1,ga_Fji,ga_h01_num

      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
c
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer l_dia, k_dia, l_para, k_para
      integer l_xyz, k_xyz, l_zan, k_zan
      integer l_tmp, k_tmp
      integer g_dens(3), g_s10, g_h01, g_d1, g_rhs, g_fock, g_u
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer vectors_scl(2) ! to be scaled with g_CiFull 
      double precision atn, tol2e, val, isotr, aniso
      double precision jfac(3),kfac(3),a(6),axs(3,3),eig(3),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      integer ld(2),cbuf ! FA-08-19-10
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
      external giao_aotomo,hnd_giasym,mat_transpose ! FA-09-16-10

      double precision val1

      logical     oskel, status
      double precision ppm
      data ppm     /26.62566914d+00/ 
      data tol2e   /1.0d-10/
      integer NoKinetic
      common /skipKinetic/NoKinetic

      if (ga_nodeid().eq.0) write(luout,9999)
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
c     
      oskel = .false.
c
c     If DFT get part of the exact exchange defined
c
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
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
      if (.not.rtdb_get(rtdb, 'giao:natoms', MT_INT, 1,sh_atom)) then  
         if (.not.geom_ncent(geom,sh_atom)) call                   
     &       errquit('hnd_giaox: geom_ncent failed ',0, GEOM_ERR)
      endif
      if (.not. ma_push_get(mt_int,sh_atom,'nmr tmp',l_tmp,k_tmp)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,3*sh_atom,'nmr at',l_xyz,k_xyz)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_xyz',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,sh_atom,'nmr zan',l_zan,k_zan)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_zan',0,MA_ERR)
c     By default set array from 1 to sh_atom (which is all when giao:natoms is not there)

      do ixy = 1, sh_atom
         int_mb(k_tmp+ixy-1) = ixy
      enddo
c
c     Try to read the atom list from rtdb. If it is not there, we still have the default list
      status = rtdb_get(rtdb,'giao:atom list',mt_int,sh_atom,
     &                  int_mb(k_tmp))
      do ixy = 0, sh_atom-1
         if (.not. geom_cent_get(geom, int_mb(k_tmp+ixy), tag, 
     &                           dbl_mb(k_xyz+3*ixy),dbl_mb(k_zan+ixy)))
     &       call errquit('hnd_giaox: geom_cent_tag failed',0, GEOM_ERR)
      enddo 
c
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('hnd_giaox: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('hnd_giaox: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('hnd_giaox: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c
c     Get Unperturbed Density Matrix
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Return if scftyp equals UHF (= ROHF)
      if (scftyp.eq.'UHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          goto 911
      endif
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nclosed(1)
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('hnd_giaox: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u)
c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
      if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,'RHS',-1,-1,g_rhs))
     &   call errquit('hnd_giaox: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
c    
c     Get S10 in GA and transform to MO set (virt,occ)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10)) call 
     &    errquit('hnd_giaox: nga_create failed g_s01',0,GA_ERR)
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,g_s10,'s10',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      call giao_aotomo(g_s10,vectors,nclosed,nvirt,1,3,nbf)
c
c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block g_s10 with - (minus) eigenvalues 
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(3) = 1
      ahi(3) = 3
      do iocc = 1, nclosed(1)
         alo(2) = iocc
         ahi(2) = iocc
         call nga_scale_patch(g_s10,alo,ahi,-dbl_mb(k_eval+iocc-1)) 
      enddo
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_occ',0,MA_ERR)
c
c     Copy to ga_rhs 
c     alo(1) and ahi(1) the same as before
      alo(2) = 1
      ahi(2) = nclosed(1)
      call nga_copy_patch('n',g_s10,alo,ahi,g_rhs,blo,bhi)
c
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
      alo(1) = 1
      ahi(1) = nclosed(1)
      call nga_scale_patch(g_s10,alo,ahi,-0.5d0)
      call nga_copy_patch('n',g_s10,alo,ahi,g_u,alo,ahi)
c
c     We also need the occupied-occupied contribution of g_u contributing
c     to the first order density matrix. As this block does not change 
c     during the CPHF we can calculate it once and subtract it from the
c     RHS. We will reuse g_s10 as scratch space.
      call ga_zero(g_s10)
      clo(1) = 3
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',chi,g_fock)) call 
     &    errquit('hnd_giaox: nga_create failed g_fock',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,clo,'D10 matrix',chi,g_d1)) call 
     &    errquit('hnd_giaox: nga_create failed g_d1',0,GA_ERR)
      call ga_zero(g_fock)
      call ga_zero(g_d1)
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
      dhi(2) = nclosed(1)
c
c     Create "perturbed density matrix" for closed-closed g_u block
      do ifld = 1, 3
         jfac(ifld) = 0.0d0
         kfac(ifld) = -1.0d0*xfac
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         dlo(3) = ifld
         dhi(3) = ifld
         ahi(1) = nmo
         ahi(2) = nclosed(1)
         bhi(2) = nmo 
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(1),blo,bhi,  
     &                         g_u,alo,ahi,g_s10,dlo,dhi)  
         ahi(2) = nbf
         ahi(1) = nclosed(1)
         bhi(2) = nclosed(1)
c
c     Minus sign as we subtract it from the RHS as we do not include 
c     it in the LHS
         call nga_matmul_patch('n','t',-1.0d0,0.0d0,vectors(1),blo,bhi,
     &                         g_s10,alo,ahi,g_d1,clo,chi)  
      enddo
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
      call shell_fock_build(geom, basis, 0, 3,
     $     jfac, kfac,tol2e, g_d1, g_fock,.false.)
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
c
c     Transform to the occ-virt MO basis and add to RHS
      call ga_zero(g_s10)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nclosed(1)
      do ifld = 1, 3
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         call nga_matmul_patch('n','n',2.0d0,0.0d0,
     $        g_fock,clo,chi,
     $        vectors(1),alo,ahi,
     $        g_s10,alo,ahi)
      enddo
      call ga_zero(g_fock)
      clo(2) = nclosed(1)+1
      clo(3) = 1
      chi(2) = nmo
      chi(3) = nclosed(1)
      do ifld = 1, 3
         blo(1) = nclosed(1)+1
         blo(2) = 1
         bhi(1) = nmo
         bhi(2) = nbf
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         call nga_matmul_patch('t','n',1.0d0,0.0d0,
     $        vectors(1), blo,bhi,
     $        g_s10, alo,ahi,
     $        g_fock, clo,chi )
         blo(1) = 1
         bhi(1) = nclosed(1)*nvirt(1)
         blo(2) = ifld
         bhi(2) = ifld
         call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,clo,chi,
     &                      g_rhs,blo,bhi)
      enddo
c
c     Cleanup of g_d1 and g_fock, not needed for now
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_fock',0,GA_ERR)
c
c     Get H10 in GA, reusing g_s10 array
      call ga_zero(g_s10)
      NoKinetic=1 ! =0 DO-kinetic, =1 SKIP-kinetic
      call int_giao_1ega(basis,basis,g_s10,'tv10',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
c ++++++++++++++++++++++++++++++++++++++++
c +++++++++++ update g_s10 +++++++++ START
c g_s10=gs10+ga_Fji
      call ga_add(1.0d0,g_s10,1.0d0,ga_Fji,g_s10)
c +++++++++++ update g_s10 +++++++++ END
c ++++++++++++++++++++++++++++++++++++++
c
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs
      call giao_aotomo(g_s10,vectors,nclosed,nvirt,1,3,nbf)
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_s10,alo,ahi,
     &                   g_rhs,blo,bhi)
c
c     Cleanup g_s10 as we do not need it right now
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_s10',0,GA_ERR)
c
c     Remaining term is Perturbed (GIAO) two-electron term times Unperturbed density
c     Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in AO basis
      alo(1) = -1 
      alo(2) = -1
      alo(3) = 1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
     &    errquit('hnd_giaox: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
      call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac)
      if(use_theory.eq.'dft') then
         ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('hess_cphf: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      endif
c
c     Transform to MO basis and add to right-hand-side
      call giao_aotomo(g_fock,vectors,nclosed,nvirt,1,3,nbf)
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,alo,ahi,
     &                   g_rhs,blo,bhi)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_fock',0,GA_ERR)
      call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)
c
c     Write ga_rhs to disk 
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_giaox: could not write cphf_rhs',0, DISK_ERR)
      call schwarz_tidy()
      call int_terminate()
c
c     Call the CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
c
      if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('hnd_giaox: failed to write skew ', 0, RTDB_ERR)
c
      if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_giaox: failure in cphf ',0, RTDB_ERR)
c
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('hnd_giaox: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_giaox: could not read cphf_rhs',0, DISK_ERR)      
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u,alo,ahi)
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_rhs',0,GA_ERR)
c
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1)) call 
     &    errquit('hnd_giaox: nga_create failed g_d1',0,GA_ERR)
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
      dhi(2) = nclosed(1)
c --------- ZORA scaling of MO vectors(1) ----- START
c Note.- g_CiFull is defined in dft_zora_scale() (source dft_zora_utils.F)
         if(.not.ga_duplicate(vectors(1),
     &                        vectors_scl(1),'vscl 1'))
     &  call errquit('hnd_giaox_Z4: ga_duplicate failed',1,GA_ERR)
        call ga_copy(vectors(1),vectors_scl(1))
        call ga_scale_cols(vectors_scl(1),g_CiFull(1))
c --------- ZORA scaling of MO vectors(1) ----- END
      do ifld = 1, 3
         alo(3) = ifld
         ahi(3) = ifld
         blo(3) = ifld
         bhi(3) = ifld
         clo(3) = ifld
         chi(3) = ifld
         dlo(3) = ifld
         dhi(3) = ifld
         bhi(1) = nbf
         bhi(2) = nmo 
         ahi(1) = nmo
         ahi(2) = nclosed(1)
c     Make C1
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(1),blo,bhi,  
     &                         g_u,alo,ahi,g_d1,dlo,dhi)  
         call nga_copy_patch('n',g_d1,dlo,dhi,g_u,dlo,dhi)
         bhi(1) = nbf
         bhi(2) = nclosed(1)
         ahi(1) = nclosed(1)
         ahi(2) = nbf
c     Make D1
         call nga_matmul_patch('n','t',2.0d0,0.0d0,
     &                         vectors_scl(1),blo,bhi,  
     &                         g_u,alo,ahi,g_d1,clo,chi)
         call nga_matmul_patch('n','t',-2.0d0,1.0d0,g_u,blo,bhi,  
     &                         vectors_scl(1),alo,ahi,g_d1,clo,chi)
      enddo
      if (.not.ga_destroy(g_u)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_u',0,GA_ERR)
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the shielding tensor
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh para',l_para,k_para)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_para',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_giaox: ma_push_get failed k_dia',0,MA_ERR)
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
c     s(para)xy = Sum(n,l) D1x(n,l) * H01y(n,l)
      do ixy = 1, 9*sh_atom
         dbl_mb(k_para+ixy-1) = 0.0d0  ! initialize the paramagnetic part
      enddo
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++ Initialize dbl_mb(k_para) with ga_para1 ++++ START
c ---- STORE: g_para1 --> dbl_mb(k_para)
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=sh_atom
      ld(1)=3
      ld(2)=3
      call nga_get(ga_para1,alo,ahi,dbl_mb(k_para),ld)
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
      do iatom = 1, sh_atom
       do iy = 1, 3
        blo(3) = blo(3) + 1
        bhi(3) = bhi(3) + 1
        do ix = 1, 3
         alo(3) = ix
         ahi(3) = ix
         ixy = ixy + 1
         val = nga_ddot_patch(g_d1 ,'n',alo,ahi,
     &                        ga_h01_num,'n',blo,bhi)
         cbuf=9*(iatom-1)+3*(ix-1)+iy-1 ! transpose
         dbl_mb(k_para+cbuf)=dbl_mb(k_para+cbuf)+val*ppm
        enddo
       enddo
      enddo
cf
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_d1',0,GA_ERR)
c
c     s(dia)xy = Sum(n,l) D0(n,l) * H11(dia)xy(n,l) 
c     Ordering of H11 blocks is Bxx,Bxy,Bxz,Byx,etc 
      do ixy = 1, 9*sh_atom
         dbl_mb(k_dia+ixy-1)  = 0.0d0  ! initialize the diamagnetic part
      enddo
c ---- STORE: ga_dia --> dbl_mb(k_dia)
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=sh_atom
       ld(1)=3
       ld(2)=3 
       call nga_get(ga_dia,alo,ahi,dbl_mb(k_dia),ld)
c +++++++++ print-total-dia-transferred +++ START
c      do iatom = 1, sh_atom
c         do ix = 1, 9
c          if (ga_nodeid().eq.0) then
c           write(*,17) iatom,ixy,dbl_mb(k_dia+ix-1)
c 17        format('dia-NUM(',i3,',',i3,')=',f15.8)
c          endif
c         enddo
c      enddo
c +++++++++ print-total-dia-transferred +++ END
      if (.not.ga_destroy(g_dens(1))) call 
     &    errquit('hnd_giaox: ga_destroy failed g_dens',0,GA_ERR)
c
c     Print out tensor information, and write to Ecce file if necessary

      call ecce_print_module_entry('nmr')
      do iatom = 1, sh_atom
         ioff = (iatom-1)*9
         if (.not. geom_cent_get(geom, int_mb(k_tmp+iatom-1), tag, 
     &       dbl_mb(k_xyz), dbl_mb(k_zan))) call 
     &       errquit('hnd_giaox: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(tag, symbol, element, atn)) call
     &       errquit('hnd_giaox: geom_tag_to_element failed',0,GEOM_ERR)
c
c      Print tensor pieces and sum for total shielding tensor
c
         if (ga_nodeid().eq.0) then
            write(luout,9700) iatom,symbol
            write(luout,9800) (dbl_mb(k_dia+ioff+ix-1),ix=1,9)
            write(luout,9801) (dbl_mb(k_para+ioff+ix-1),ix=1,9)
         endif
         do ix = 0, 8 
            dbl_mb(k_para+ioff+ix) = dbl_mb(k_dia+ioff+ix) + 
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
c
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
            call ecce_print1('shielding isotropic',MT_DBL,isotr,1)
            call ecce_print1('shielding anisotropy',MT_DBL,aniso,1)
            call ecce_print1('shielding eigenvalues',MT_DBL,eig,3)
            call ecce_print2('shielding eigenvectors',MT_DBL,axs,
     &                       3,3,3)
         endif
      enddo
      call ecce_print_module_exit('nmr','ok')
c
      call ga_sync()
c
c     Clean up all remaining memory
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_dia',0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_para',0,MA_ERR)
 911  if (.not.ga_destroy(vectors(1))) call 
     &    errquit('giao_aotomo: ga_destroy failed vectors',0,GA_ERR)
      if (.not.ga_destroy(vectors_scl(1))) call 
     &    errquit('giao_aotomo: ga_destroy failed vectors',0,GA_ERR)
c
      if (.not.ma_pop_stack(l_zan)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_zan',0,MA_ERR)
      if (.not.ma_pop_stack(l_xyz)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_xyz',0,MA_ERR)
      if (.not.ma_pop_stack(l_tmp)) call
     &    errquit('hnd_giaox: ma_pop_stack failed l_tmp',0,MA_ERR)
      call schwarz_tidy()
      call int_terminate()
c
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
c
