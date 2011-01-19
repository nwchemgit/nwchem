      subroutine fiao_b1_movecs(rtdb,basis,geom, g_vecB1)

c $Id: $

c     This routine is a modification of fiao_b1_movecs.F. We calculate
c     the MO vectors perturbed by a STATIC external magnetic field but
c     not with GIAOs. This is a veavily stripped-down version of
c     hnd_giaox.

c     g_vecB1 should have at least dimensions (nbf, nocc, 3)

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
c
      integer rtdb    ! [in] rtdb handle
      integer basis   ! [in] basis handle
      integer geom    ! [in] geometry handle
      integer g_vecB1 ! [out] B-field perturbed MO coefficients GIAO

      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval

      integer g_dens(3), g_s10, g_d1, g_rhs, g_fock, g_u
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val
      double precision jfac(3),kfac(3),a(6),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol

      double precision origin(3)
      data origin/0d0,0d0,0d0/

      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf

      logical     oskel, status, debug
      double precision ppm
      data tol2e   /1.0d-10/

c     ==================================================================
c
      if (ga_nodeid().eq.0) write(luout,9999)

      call ga_zero(g_vecB1) ! initialize array for main result

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
      status = rtdb_parallel(.true.)

c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues
c
      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('fiao_b1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('fiao_b1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('fiao_b1: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('fiao_b1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('fiao_b1: ma_pop_stack failed k_occ',0,MA_ERR)
c
c     Get Unperturbed Density Matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Error exit if scftyp equals UHF (= ROHF)
c
      if (scftyp.eq.'UHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          call errquit('fiao_b1: incompatible SCF type for Response',
     &       0,INPUT_ERR)
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
     &    errquit('fiao_b1: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u)
c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
c
      if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,'RHS',-1,-1,g_rhs))
     &   call errquit('fiao_b1: ga_create failed g_rhs',0,GA_ERR)
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
     &    errquit('fiao_b1: nga_create failed g_s01',0,GA_ERR)


c     Get H10 in GA, using g_s10 array

      call ga_zero(g_s10)

      call int_giao_1ega(basis,basis,g_s10,'angmom',origin,1,
     &   oskel)
      call ga_scale (g_s10, 0.5d0)

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
c
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('fiao_b1: ga_destroy failed g_s10',0,GA_ERR)


c     there is a division by -4 somewhere in the
c     CPKS solver so we scale the rhs by -4 in order to
c     get the correct magnitude of the result back
      
      call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)

c     Write ga_rhs to disk 

      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('fiao_b1: could not write cphf_rhs',0, DISK_ERR)
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
     $   errquit('fiao_b1: failed to write skew ', 0, RTDB_ERR)
c
      if (.not.cphf2(rtdb)) call errquit
     $  ('fiao_b1: failure in cphf ',0, RTDB_ERR)
c
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('fiao_b1: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
c
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('fiao_b1: could not read cphf_rhs',0, DISK_ERR)      
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u,alo,ahi)
c
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('fiao_b1: ga_destroy failed g_rhs',0,GA_ERR)
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
     &    errquit('fiao_b1: nga_create failed g_d1',0,GA_ERR)
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
         call nga_copy_patch('n',g_d1,dlo,dhi,g_vecB1,dlo,dhi)

c        This patch of g_vecB1 now has the perturbed MO
c        coefficients. let's print them for debug purposes:

         if (debug) then
           write (luout,*) 'fiao_b1: perturbed C, direction ',ifld
           call nga_print_patch(g_vecB1,dlo,dhi,1)
         end if

       enddo ! ifld

      if (.not.ga_destroy(g_u)) call 
     &    errquit('fiao_b1: ga_destroy failed g_u',0,GA_ERR)
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('fiao_b1: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(vectors(1))) call 
     &   errquit('fiao_b1: ga_destroy failed vectors',0,GA_ERR)
      if (.not.ga_destroy(g_dens(1))) call 
     &    errquit('fiao_b1: ga_destroy failed g_dens',0,GA_ERR)

c     Alll done.
c     At this point, we don't need to terminate the integrals.
c     They were terminated by the cphf.

      call ga_sync()

      return

 7000 format(/,10x,'B-field perturbed MOs cannot be calculated for',
     1      ' UHF or ROHF wave functions at present')
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9999 format(
     1 /,10x,54(1h-),/,
     2 10x,'Calculating magnetic field perturbed MO vectors (FIAO)',/,
     3 10x,54(1h-),/)

      end

