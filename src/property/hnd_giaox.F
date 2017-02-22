      subroutine hnd_giaox(rtdb,basis,geom)
c $Id$
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
#include "util.fh"
#include "inp.fh"
c
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
      integer g_dens(3), g_s10, g_h01, g_h11, g_d1, g_rhs, g_fock, g_u
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val, isotr, aniso
      double precision jfac(3),kfac(3),a(6),axs(3,3),eig(3),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
c
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
c
      logical     oskel, status, debug
      double precision ppm
      data ppm     /26.62566914d+00/ 
      data tol2e   /1.0d-10/

      integer debug_giaox

      integer NoKinetic,npol
      common /skipKinetic/NoKinetic
c     bq charges
      integer nbq,nextbq,ncosbq
      nbq = 0
      nextbq = 0
      ncosbq = 0

      debug_giaox=0 ! =1 for debugging print outs of matrices
c
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, npol))
     $   npol = 1 ! Default nr of polarizations
      if (ga_nodeid().eq.0) write(luout,9999)

      debug = .false. .and. (ga_nodeid().eq.0) ! special debugging
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
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
c     Find out from rtdb which atoms we need to calculate shielding for
c     Get number of atoms (all or number from rtdb)
c     Get which atoms (all or some read from rtdb)
c     Allocate arrays which will hold atomic information (k_zan and k_xyz)
c
      status = rtdb_parallel(.true.)
      if (.not.rtdb_get(rtdb, 'giao:natoms', MT_INT, 1,sh_atom)) then  
         if (.not.geom_ncent(geom, sh_atom)) call                   
     &       errquit('hnd_giaox: geom_ncent failed ',0, GEOM_ERR)
      endif
      if (.not. ma_push_get(mt_int,sh_atom,'nmr tmp',l_tmp,k_tmp)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_tmp',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,3*sh_atom,'nmr at',l_xyz,k_xyz)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_xyz',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,sh_atom,'nmr zan',l_zan,k_zan)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_zan',0,MA_ERR)
c
c     By default set array from 1 to sh_atom (which is all when giao:natoms is not there)
c
      do ixy = 1, sh_atom
         int_mb(k_tmp+ixy-1) = ixy
      enddo
c
c     Try to read the atom list from rtdb. If it is not there, we still have the default list
c
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
c
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
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Return if scftyp equals UHF (= ROHF)
c
      if (scftyp.eq.'UHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          goto 911
      endif
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c
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
c
      if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,'RHS',-1,-1,g_rhs))
     &   call errquit('hnd_giaox: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
c
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
c    
c     Get S10 in GA and transform to MO set (virt,occ)
c
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
c
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
c
      alo(2) = 1
      ahi(2) = nclosed(1)
      call nga_copy_patch('n',g_s10,alo,ahi,g_rhs,blo,bhi)
c
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
c
      alo(1) = 1
      ahi(1) = nclosed(1)
      call nga_scale_patch(g_s10,alo,ahi,-0.5d0)
      call nga_copy_patch('n',g_s10,alo,ahi,g_u,alo,ahi)
c
c     We also need the occupied-occupied contribution of g_u contributing
c     to the first order density matrix. As this block does not change 
c     during the CPHF we can calculate it once and subtract it from the
c     RHS. We will reuse g_s10 as scratch space.
c
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
c
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
c
         call nga_matmul_patch('n','t',-1.0d0,0.0d0,vectors(1),blo,bhi,
     &                         g_s10,alo,ahi,g_d1,clo,chi)  
      enddo
c
c     Build "first order fock matrix"
c
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
c
c     Note: Just the exchange: jfac = 0.d0 (see above)
c
      if (.not.cam_exch) then
         call shell_fock_build(geom, basis, 0, 3,
     $     jfac, kfac, tol2e, g_d1, g_fock, .false.)
      else
         call shell_fock_build_cam(geom, basis, 0, 3,
     $     jfac, kfac, tol2e, g_d1, g_fock, .false.)
      end if
c
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
c
c     Transform to the occ-virt MO basis and add to RHS
c
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
     $             g_s10, alo,ahi,
     $            g_fock, clo,chi )
         blo(1) = 1
         bhi(1) = nclosed(1)*nvirt(1)
         blo(2) = ifld
         bhi(2) = ifld
         call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,clo,chi,
     &                      g_rhs,blo,bhi)

      if (debug) then
        write (luout,*) 'Fock matrix from S(1) in U-vector'
        call nga_print_patch(g_fock,clo,chi,1)
      end if

      enddo
c
c     Cleanup of g_d1 and g_fock, not needed for now
c
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_fock',0,GA_ERR)
c
c     Get H10 in GA, reusing g_s10 array
c
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,g_s10,'l10',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      NoKinetic=0 ! =0 DO-kinetic, =1 SKIP-kinetic
      call int_giao_1ega(basis,basis,g_s10,'tv10',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
c
c     Get external and cosmo bq contribution
c
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
        call int_giao_1ega(basis,basis,g_s10,'bq10',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      end if
c
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs
c
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
c
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_s10',0,GA_ERR)
c
c     Remaining term is Perturbed (GIAO) two-electron term times Unperturbed density
c     Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in AO basis
c
      alo(1) = -1 
      alo(2) = -1
      alo(3) = 1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
     &    errquit('hnd_giaox: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
c
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif
c

       if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0) then
         write(*,*) 'Aft. add_H10()'
         write(*,*) '---- g_rhs-BEF-new_giao_2e -------- START'
        endif
        call ga_print(g_rhs)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_rhs-BEF-new_giao_2e -------- END'
       endif

      call new_giao_2e(geom,basis,nbf,tol2e,g_dens,g_fock,xfac,npol)
c
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
c
      call giao_aotomo(g_fock,vectors,nclosed,nvirt,1,3,nbf)

      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-giao_aotomo -------- START'
       call ga_print(g_fock)
       if (ga_nodeid().eq.0)
     &   write(*,*) '---- g_fock-aft-giao_aotomo -------- END'
      endif

      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,alo,ahi,
     &                   g_rhs,blo,bhi)


      if (debug) then
        write (luout,*) 'Fock matrix from giao G(1)'
        call nga_print_patch(g_fock,alo,ahi,1)
      end if

      if (.not.ga_destroy(g_fock)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_fock',0,GA_ERR)
      call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)

      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- START'
       call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading INPUT-cphf: g_rhs -------- END'
      endif

c
c     Write ga_rhs to disk 
c
      call cphf_fname('cphf_rhs',cphf_rhs)
      call cphf_fname('cphf_sol',cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_giaox: could not write cphf_rhs',0, DISK_ERR)
c
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
c
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_giaox: could not read cphf_rhs',0, DISK_ERR) 

      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- START'
      call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- Reading sol-cphf: g_rhs -------- END'  
      endif

     
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u,alo,ahi)
c
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_rhs',0,GA_ERR)
c
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1)) call 
     &    errquit('hnd_giaox: nga_create failed g_d1',0,GA_ERR)
c
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
c
c     Make C1
c
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(1),blo,bhi,  
     &                         g_u,alo,ahi,g_d1,dlo,dhi)  
         call nga_copy_patch('n',g_d1,dlo,dhi,g_u,dlo,dhi)

c ... jochen: this patch of g_u now has the perturbed MO
c        coefficients. let's print that for debug purposes:
c         if (debug) then
c           write (luout,*) 'giaox: perturbed C, direction ',ifld
c           call nga_print_patch(g_u,dlo,dhi,1)
c         end if

         bhi(1) = nbf
         bhi(2) = nclosed(1)
         ahi(1) = nclosed(1)
         ahi(2) = nbf
c
c     Make D1
c
         call nga_matmul_patch('n','t',2.0d0,0.0d0,vectors(1),blo,bhi,  
     &                      g_u,alo,ahi,g_d1,clo,chi)
         call nga_matmul_patch('n','t',-2.0d0,1.0d0,g_u,blo,bhi,  
     &                      vectors(1),alo,ahi,g_d1,clo,chi)
      enddo

      if (debug_giaox.eq.1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_d1---------- START'
        call ga_print(g_d1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '------g_d1---------- END'
      endif

      if (.not.ga_destroy(g_u)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_u',0,GA_ERR)
c
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the shielding tensor
c
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh para',l_para,k_para)) 
     &    call errquit('hnd_giaox: ma_push_get failed k_para',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,9*sh_atom,'sh dia',l_dia,k_dia)) call
     &    errquit('hnd_giaox: ma_push_get failed k_dia',0,MA_ERR)
c
c     Before we start getting the integrals we need to reinitialize the
c     integrals. They were terminated by the cphf.
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
c
c     s(para)xy = Sum(n,l) D1x(n,l) * H01y(n,l)
c
      do ixy = 1, 9*sh_atom
         dbl_mb(k_para+ixy-1) = 0.0d0  ! initialize the paramagnetic part
      enddo
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*sh_atom
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,g_h01)) call 
     &    errquit('hnd_giaox: nga_create failed g_h01',0,GA_ERR)
      call ga_zero(g_h01)
      call int_giao_1ega(basis,basis,g_h01,'h01',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      if (debug_giaox.eq.1) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '------ g_h01 --------- START'
       call ga_print(g_h01)
       if (ga_nodeid().eq.0)
     &  write(*,*) '------ g_h01 --------- END'
      endif
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

               dbl_mb(k_para+ixy-1) = nga_ddot_patch(g_d1 ,'n',alo,ahi,
     &                                               g_h01,'n',blo,bhi)
               dbl_mb(k_para+ixy-1) = dbl_mb(k_para+ixy-1) * ppm

      if (debug_giaox.eq.1) then
               if (ga_nodeid().eq.0) then
                write(*,12) iatom,iy,ix,
     &                      alo(1),alo(2),alo(3),
     &                      ahi(1),ahi(2),ahi(3),
     &                      blo(1),blo(2),blo(3),
     &                      bhi(1),bhi(2),bhi(3),dbl_mb(k_para+ixy-1)
 12             format('(iatom,iy,ix)=(',i3,',',i3,',',i3,') ',
     &                  'alo=(',i3,',',i3,',',i3,') ',
     &                  'ahi=(',i3,',',i3,',',i3,') ',
     &                  'blo=(',i3,',',i3,',',i3,') ',
     &                  'bhi=(',i3,',',i3,',',i3,') ',
     &                  'para=',f15.8)
               endif
       endif

            enddo
         enddo
         call mat_transpose(dbl_mb(k_para+ixy-9)) ! transpose the paramagnetic part
      enddo
cf
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_h01)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_h01',0,GA_ERR)
c
c     s(dia)xy = Sum(n,l) D0(n,l) * H11(dia)xy(n,l) 
c     Ordering of H11 blocks is Bxx,Bxy,Bxz,Byx,etc 
c
      do ixy = 1, 9*sh_atom
         dbl_mb(k_dia+ixy-1)  = 0.0d0  ! initialize the diamagnetic part
      enddo
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 9*sh_atom
      if (.not.nga_create(MT_DBL,3,ahi,'H11 matrix',alo,g_h11)) call 
     &    errquit('hnd_giaox: nga_create failed g_h11 all',0,GA_ERR)
c
      call ga_zero(g_h11)
      call int_giao_1ega(basis,basis,g_h11,'h11 para',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      ixy = 0
      do iatom = 1, sh_atom
         do ix = 1, 9
            ixy = ixy + 1
            alo(3) = ixy
            ahi(3) = ixy
            val=nga_ddot_patch(g_dens(1),'n',blo,bhi,g_h11,'n',alo,ahi)
            dbl_mb(k_dia+ixy-1) = dbl_mb(k_dia+ixy-1) + val * ppm 
         enddo
      enddo
c       
c     s(dia)xy  = s(dia)xy + Sum(n,l) D0(n,l) * H11(dia)xy(n,l)
c
      call ga_zero(g_h11)
      call int_giao_1ega(basis,basis,g_h11,'h11 dia',dbl_mb(k_xyz),
     &                   sh_atom,oskel)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      ixy = 0
      do iatom = 1, sh_atom
         do ix = 1, 9
            ixy = ixy + 1
            alo(3) = ixy
            ahi(3) = ixy
            val=nga_ddot_patch(g_dens(1),'n',blo,bhi,g_h11,'n',alo,ahi)
            dbl_mb(k_dia+ixy-1) =  dbl_mb(k_dia+ixy-1) + val * ppm 
         enddo
      enddo
c
      if (.not.ga_destroy(g_dens(1))) call 
     &    errquit('hnd_giaox: ga_destroy failed g_dens',0,GA_ERR)
      if (.not.ga_destroy(g_h11)) call 
     &    errquit('hnd_giaox: ga_destroy failed g_h11',0,GA_ERR)
c
c     Print out tensor information, and write to Ecce file if necessary
c
      call ecce_print_module_entry('nmr')
c
      do iatom = 1, sh_atom
         ioff = (iatom-1)*9
         if (.not. geom_cent_get(geom, int_mb(k_tmp+iatom-1), tag, 
     &     dbl_mb(k_xyz), dbl_mb(k_zan))) call 
     &     errquit('hnd_giaox: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(tag, symbol, element, atn)) then
           if (.not. inp_compare(0,tag(1:2),'bq')) call ! check for bq as a fall back
     &     errquit('hnd_giaox: geom_tag_to_element failed',0,GEOM_ERR)
         end if  ! tag check
c
c     Symmetrize the diamagnetic and paramagnetic tensors
c
         call hnd_giasym(dbl_mb(k_dia+ioff),dbl_mb(k_para+ioff))
c
c     Print tensor pieces and sum for total shielding tensor
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
      enddo ! iatom
c
      call ecce_print_module_exit('nmr','ok')
c
      call ga_sync()
c
c     Clean up all remaining memory
c
      if (.not.ma_pop_stack(l_dia)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_dia',0,MA_ERR)
      if (.not.ma_pop_stack(l_para)) call
     &    errquit('hnd_giaox: ma_pop_stack failed k_para',0,MA_ERR)
c
 911  if (.not.ga_destroy(vectors(1))) call 
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
c     symmetrize matrix
      subroutine hnd_giasym(dia,para)
      implicit none
      double precision dia(3,3), para(3,3)
      double precision val
      integer i,j
c
      do i = 1, 2
         do j = i+1, 3
            val = (dia(i,j) + dia(j,i)) / 2.0d0
            dia(i,j) = val
            dia(j,i) = val
            val = (para(i,j) + para(j,i)) / 2.0d0
            para(i,j) = val
            para(j,i) = val
         enddo
      enddo
      return
      end
c
c     transpose matrix
      subroutine mat_transpose(mat)
      implicit none
      double precision mat(3,3)
      double precision val
      integer i,j
c
      do i = 1, 2
         do j = i+1, 3
            val = mat(i,j)
            mat(i,j) = mat(j,i)
            mat(j,i) = val
         enddo
      enddo
      return
      end
