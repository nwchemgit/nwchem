      subroutine hnd_psoten(rtdb,geom,basis,vectors,nclosed,nopen,nvirt,
     &                      nbf,nmo,pairlist,translate,ipairs,unique,
     &                      i_pert,i_resp,tensor,coords,npso)
c $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "apiP.fh"
c
      integer rtdb                             ! [input] rtdb handle
      integer basis                            ! [input] basis handle
      integer geom                             ! [input] geometry handle
      integer vectors(2)                       ! [input] vectors
      integer nclosed(2), nopen(2), nvirt(2)   ! [input] occupation info
      integer nbf, nmo                         ! [input] orbital info
      integer ipairs                           ! [input] number of spin-spin pairs
      integer pairlist(2*ipairs)               ! [input] list of the pairs
      integer translate(2*ipairs)              ! [input] translation to unique list
      integer i_pert, i_resp                   ! [input] # of unique responding and perturbing atoms
      integer unique(i_pert+i_resp)            ! [input] list of unique atoms
      double precision tensor(3,3,5,ipairs)    ! [output] spin-spin tensor, one for each term
      double precision coords(3,i_pert+i_resp) ! [input] coordinates of unique atoms
      double precision npso                    ! [input] pso prefactor
c
      integer ix, iy, ifld
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer g_fca, g_fcb, ii
      integer g_d1(3),g_rhs,g_u(2)
      integer i, iat, jat
      double precision tol2e
      character*256 cphf_rhs, cphf_sol
c
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
c
      logical     oskel
      double precision valuea, valueb
      data tol2e   /1.0d-16/
c
      oskel = .false.
c
c     Integral initialization
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nopen(1)
      ahi(3) = 3*i_pert
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u(1))) call 
     &    errquit('hnd_pso: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u(1))
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u(2))) call 
     &    errquit('hnd_pso: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u(2))
c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,atom)
c
      if(.not.ga_create(MT_DBL,2*nopen(1)*nvirt(1),3*i_pert,'RHS',
     &                  -1,-1,g_rhs))
     &   call errquit('hnd_pso: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)
c
c     Get PSO integrals
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*i_pert
c
c     Allocate memory
c
      if (.not.nga_create(MT_DBL,3,ahi,'fca matrix',alo,g_fca)) call
     &    errquit('hnd_giaox: nga_create failed g_fca',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'fcb matrix',alo,g_fcb)) call
     &    errquit('hnd_giaox: nga_create failed g_fcb',0,GA_ERR)
c
c     Calculate integrals
c
      call ga_zero(g_fca)
      call int_giao_1ega(basis,basis,g_fca,'pso',coords,i_pert,oskel)
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
c
      blo(1) = 1
      bhi(1) = nopen(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3*i_pert
c
c     ga_rhs(a,i) = ga_rhs(a,i) + FC(a,i)
c     Transform FC and SD to MO and add to g_rhs
c
      call ga_copy(g_fca,g_fcb)
      call giao_aotomo(g_fca,vectors(1),nopen(1),nvirt(1),1,3*i_pert,
     &                 nbf)
      call giao_aotomo(g_fcb,vectors(2),nopen(2),nvirt(2),1,3*i_pert,
     &                 nbf)
      alo(1) = nopen(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nopen(1)
      alo(3) = 1
      ahi(3) = 3*i_pert
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fca,alo,ahi,
     &                   g_rhs,blo,bhi)
      alo(1) = nopen(2)+1
      ahi(2) = nopen(2)
      blo(1) = blo(1)+nopen(1)*nvirt(1)
      bhi(1) = bhi(1)+nopen(2)*nvirt(2)
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fcb,alo,ahi,
     &                   g_rhs,blo,bhi)
c
      call ga_scale(g_rhs,-4.0d0)
c
c     Some cleanup
c
      if (.not.ga_destroy(g_fcb)) call 
     &    errquit('hnd_pso: ga_destroy failed g_fcb',0,GA_ERR)
c
c     Write ga_rhs to disk 
c
      call cphf_fname('cphf_rhs',cphf_rhs)
      call cphf_fname('cphf_sol',cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_pso: could not write cphf_rhs',0, DISK_ERR)
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
     $   errquit('hnd_pso: failed to write skew ', 0, RTDB_ERR)
c
      if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_pso: failure in cphf ',0, RTDB_ERR)
c
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('hnd_pso: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_pso: could not read cphf_rhs',0, DISK_ERR)      
      alo(1) = nopen(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nopen(1)
      alo(3) = 1
      ahi(3) = 3*i_pert
      blo(1) = 1
      bhi(1) = nopen(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3*i_pert
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u(1),alo,ahi)
      blo(1) = blo(1)+nopen(1)*nvirt(1)
      bhi(1) = bhi(1)+nopen(2)*nvirt(2)
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u(2),alo,ahi)
c
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_pso: ga_destroy failed g_rhs',0,GA_ERR)
c
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = [(C1*C0+) + (C0*C1+)]
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*i_pert
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1(1))) call 
     &    errquit('hnd_pso: nga_create failed g_u',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1(2))) call 
     &    errquit('hnd_pso: nga_create failed g_u',0,GA_ERR)
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
      dhi(1) = nbf
      dlo(2) = 1
      do ii = 1, 2
         dhi(2) = nopen(ii)
      do ifld = 1, 3*i_pert
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
         ahi(2) = nopen(ii)
c
c     Make C1
c
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(ii),blo,bhi,
     &                         g_u(ii),alo,ahi,g_d1(ii),dlo,dhi)  
         call nga_copy_patch('n',g_d1(ii),dlo,dhi,g_u(ii),dlo,dhi)
         bhi(1) = nbf
         bhi(2) = nopen(ii)
         ahi(1) = nopen(ii)
         ahi(2) = nbf
c
c     Make D1
c
         call nga_matmul_patch('n','t',1.0d0,0.0d0,vectors(ii),blo,bhi,
     &                      g_u(ii),alo,ahi,g_d1(ii),clo,chi)
         call nga_matmul_patch('n','t',-1.0d0,1.0d0,g_u(ii),blo,bhi,  
     &                      vectors(ii),alo,ahi,g_d1(ii),clo,chi)
      enddo
      enddo
c
      if (.not.ga_destroy(g_u(1))) call 
     &    errquit('hnd_pso: ga_destroy failed g_u(1)',0,GA_ERR)
      if (.not.ga_destroy(g_u(2))) call 
     &    errquit('hnd_pso: ga_destroy failed g_u(2)',0,GA_ERR)
c
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the FC term
c
      if (.not.ga_destroy(g_fca)) call 
     &    errquit('hnd_pso: ga_destroy failed g_fca2',0,GA_ERR)
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*ipairs
c
c     Allocate memory
c
      if (.not.nga_create(MT_DBL,3,ahi,'fca matrix',alo,g_fca)) call
     &    errquit('hnd_psoten: nga_create failed g_fca',0,GA_ERR)

      call ga_zero(g_fca)
      call int_giao_1ega(basis,basis,g_fca,'pso',coords(1,i_pert+1),
     &                   i_resp,oskel)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      do i = 1, ipairs
c        iat = unique(translate(i))
c        jat = unique(translate(i+ipairs))
         iat = translate(i)
         jat = translate(i+ipairs)
         do ix = 1, 3
            alo(3) = (iat-1)*3+ix
            ahi(3) = (iat-1)*3+ix
            do iy = 1, 3
               blo(3) = (jat-1)*3+iy 
               bhi(3) = (jat-1)*3+iy
          valuea=nga_ddot_patch(g_d1(1),'n',alo,ahi,g_fca,'n',blo,bhi)
          valueb=nga_ddot_patch(g_d1(2),'n',alo,ahi,g_fca,'n',blo,bhi)
            tensor(ix,iy,4,i)=-(valuea+valueb)*npso*0.5d0
c           tensor(ix,iy,4,i)=-(valuea)*npso
c           if(ga_nodeid().eq.0) write(6,'(A,9I4,5F12.6)') 'PSO ij',i,
c    & pairlist(i),pairlist(i+ipairs),iat,jat,alo(3),blo(3),
c    &                 ix,iy,tensor(ix,iy,4,i),-valuea*npso,-valueb*npso
            enddo
         enddo
      enddo
c
      call ga_sync()
c
c     Clean up all remaining memory
c
      if (.not.ga_destroy(g_d1(1))) call 
     &    errquit('hnd_pso: ga_destroy failed g_d1(1)',0,GA_ERR)
      if (.not.ga_destroy(g_d1(2))) call 
     &    errquit('hnd_pso: ga_destroy failed g_d1(2)',0,GA_ERR)
      if (.not.ga_destroy(g_fca)) call 
     &    errquit('hnd_pso: ga_destroy failed g_fca',0,GA_ERR)
c
      call schwarz_tidy()
      call int_terminate()
c
      return
      end
