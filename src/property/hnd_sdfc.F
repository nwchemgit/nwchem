      subroutine hnd_sdfc(rtdb,geom,basis,vectors,nclosed,nopen,nvirt,
     &                    nbf,nmo,pairlist,translate,ipairs,unique,
     &                    i_pert,i_resp,tensor,coords,nfc,nsd,nsdfc)
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
      double precision nfc, nsd, nsdfc         ! [input] prefactors for the three terms
c
      integer ixy, ix, iy, iz, ifld
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer g_fca, g_fcb, g_sda, g_sdb, ii
      integer g_d1(3),g_rhs,g_u(2)
      integer i, j, i_total
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
      double precision pifac, froth
c
      integer ilist(3,3)
      data ilist /1,4,5, 4,2,6, 5,6,3/
c
      parameter(froth=4.0d0/3.0d0) 
      parameter(pifac=froth*3.14159265358979323846264338327950288419d0)
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
c     Total number of responses to be calculated
c     i_pert FC + 6 * i_pert SD + 6 * i_resp SD-FC
c
      i_total = 7*i_pert + 6*i_resp
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nopen(1)
      ahi(3) = i_total
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u(1))) call 
     &    errquit('hnd_sdfc: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u(1))
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u(2))) call 
     &    errquit('hnd_sdfc: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u(2))
c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,atom)
c
      if(.not.ga_create(MT_DBL,2*nopen(1)*nvirt(1),i_total,'RHS',
     &                  -1,-1,g_rhs))
     &   call errquit('hnd_sdfc: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)
c
c     Get FC and SD+FC in GA 
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = i_pert+i_resp
c
c     Allocate memory
c
      if (.not.nga_create(MT_DBL,3,ahi,'fca matrix',alo,g_fca)) call
     &    errquit('hnd_giaox: nga_create failed g_fca',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'fcb matrix',alo,g_fcb)) call
     &    errquit('hnd_giaox: nga_create failed g_fca',0,GA_ERR)
      ahi(3) = 6*(i_pert+i_resp)
      if (.not.nga_create(MT_DBL,3,ahi,'sda matrix',alo,g_sda)) call
     &    errquit('hnd_giaox: nga_create failed g_sda',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'sda matrix',alo,g_sdb)) call
     &    errquit('hnd_giaox: nga_create failed g_sda',0,GA_ERR)
c
c     Calculate integrals
c
      call ga_zero(g_fca)
      call ga_zero(g_sda)
      call int_giao_1ega(basis,basis,g_fca,'fc',coords,
     &                   i_pert+i_resp,oskel)
      call int_giao_1ega(basis,basis,g_sda,'sd+fc',coords,
     &                   i_pert+i_resp,oskel)
c
c     Take out FC part from SD+FC (i.e. from xx, yy and zz components)
c
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      do ixy = 1, i_pert+i_resp
         alo(3) = 1+(ixy-1)*6
         ahi(3) = 1+(ixy-1)*6
         blo(3) = ixy
         bhi(3) = ixy
         do i = 1, 3
            call nga_add_patch(1.0d0,g_sda,alo,ahi,pifac,g_fca,blo,bhi,
     &                         g_sda,alo,ahi)
            alo(3) = alo(3)+1
            ahi(3) = ahi(3)+1
         enddo
      enddo
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
c
      blo(1) = 1
      bhi(1) = nopen(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = i_pert
c
c     Now that all the work is done on the integrals in Alpha, copy to Beta
c
      call ga_copy(g_fca,g_fcb)
      call ga_copy(g_sda,g_sdb)
c
c     ga_rhs(a,i) = ga_rhs(a,i) + FC(a,i)
c     Transform FC and SD to MO and add to g_rhs
c
      call giao_aotomo(g_fca,vectors(1),nopen(1),nvirt(1),1,i_pert,nbf)
      call giao_aotomo(g_fcb,vectors(2),nopen(2),nvirt(2),1,i_pert,nbf)
      call giao_aotomo(g_sda,vectors(1),nopen(1),nvirt(1),1,
     &                 6*(i_pert+i_resp),nbf)
      call giao_aotomo(g_sdb,vectors(2),nopen(2),nvirt(2),1,
     &                 6*(i_pert+i_resp),nbf)
c
c     Add to g_rhs
c
      alo(1) = nopen(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nopen(1)
      alo(3) = 1
      ahi(3) = i_pert 
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,0.5d0,g_fca,alo,ahi,
     &                   g_rhs,blo,bhi)
      blo(1) = blo(1)+nopen(1)*nvirt(1)
      bhi(1) = bhi(1)+nopen(2)*nvirt(2)
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,-0.5d0,g_fcb,alo,ahi,
     &                   g_rhs,blo,bhi)
c
      blo(1) = 1
      bhi(1) = nopen(1)*nvirt(1)
      blo(2) = i_pert+1
      bhi(2) = i_pert+6*(i_pert+i_resp)
      ahi(3) = 6*(i_pert+i_resp)
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,0.5d0,g_sda,alo,ahi,
     &                   g_rhs,blo,bhi)
      blo(1) = blo(1)+nopen(1)*nvirt(1)
      bhi(1) = bhi(1)+nopen(2)*nvirt(2)
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,-0.5d0,g_sdb,alo,ahi,
     &                   g_rhs,blo,bhi)
      call ga_scale(g_rhs,-4.0d0)
c
c     Some memory cleanup here
c
      if (.not.ga_destroy(g_fcb)) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_fcb',0,GA_ERR)
      if (.not.ga_destroy(g_sdb)) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_sdb',0,GA_ERR)
c
c     Write ga_rhs to disk 
c
      call cphf_fname('cphf_rhs',cphf_rhs)
      call cphf_fname('cphf_sol',cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('hnd_sdfc: could not write cphf_rhs',0, DISK_ERR)
c
      call schwarz_tidy()
      call int_terminate()
c
c     Call the CPHF routine
c     
      if (.not.cphf2(rtdb)) call errquit
     $  ('hnd_sdfc: failure in cphf ',0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
c
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('hnd_sdfc: could not read cphf_rhs',0, DISK_ERR)      
      alo(1) = nopen(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nopen(1)
      alo(3) = 1
      ahi(3) = i_total
      blo(1) = 1
      bhi(1) = nopen(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = i_total
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u(1),alo,ahi)
      blo(1) = blo(1)+nopen(1)*nvirt(1)
      bhi(1) = bhi(1)+nopen(2)*nvirt(2)
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u(2),alo,ahi)
c
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_rhs',0,GA_ERR)
c
c     From U matrices, generate the perturbed density matrices D1
c     C1 = C0 * U10
c     D1 = [(C1*C0+) + (C0*C1+)]
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = i_total
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1(1))) call 
     &    errquit('hnd_sdfc: nga_create failed g_d1(1)',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1(2))) call 
     &    errquit('hnd_sdfc: nga_create failed g_d1(2)',0,GA_ERR)
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
      do ifld = 1, i_total
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
         call nga_matmul_patch('n','t',1.0d0,0.0d0,g_u(ii),blo,bhi,  
     &                      vectors(ii),alo,ahi,g_d1(ii),clo,chi)
      enddo
      enddo
c
      call ga_sync()
c
      if (.not.ga_destroy(g_u(1))) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_u(1)',0,GA_ERR)
      if (.not.ga_destroy(g_u(2))) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_u(2)',0,GA_ERR)
c
c     Now we have in g_d1(nmo,nmo,3) the derivative densities and
c     hence we can calculate the contributions to the FC term
c
      call ga_zero(g_fca)
      call int_giao_1ega(basis,basis,g_fca,'fc',coords(1,i_pert+1),
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
c        alo(3) = unique(translate(i))
c        ahi(3) = unique(translate(i))
c        blo(3) = unique(translate(i+ipairs))
c        bhi(3) = unique(translate(i+ipairs))
         alo(3) = translate(i)
         ahi(3) = translate(i)
         blo(3) = translate(i+ipairs)
         bhi(3) = translate(i+ipairs)
         valuea=nga_ddot_patch(g_d1(1),'n',alo,ahi,g_fca,'n',blo,bhi)
         valueb=nga_ddot_patch(g_d1(2),'n',alo,ahi,g_fca,'n',blo,bhi)
c              if (ga_nodeid().eq.0) 
c    &            write(6,'(A,2I4,4F12.6)') 'FC tens i j com',
c    &            pairlist(ixy),pairlist(ixy+ipairs),
c    &            tensor(1,1,1,ixy),(valuea-valueb)*0.5d0*nfc,
c    &             valuea,valueb
         do ixy = 1, 3
            tensor(ixy,ixy,1,i)=(valuea-valueb)*nfc*0.5d0
         enddo
      enddo
c
c     Calculate the SD-SD tensor. First calculate the sd+fc integrals
c     and subtract the fc contribution
c
      call ga_zero(g_sda)
      call int_giao_1ega(basis,basis,g_sda,'sd+fc',coords(1,i_pert+1),
     &                   i_resp,oskel)
      alo(1) = 1
      alo(2) = 1
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf
      do i = 1, i_resp
         alo(3) = 1+(i-1)*6
         ahi(3) = 1+(i-1)*6
         blo(3) = i
         bhi(3) = i
         do ixy = 1, 3
            call nga_add_patch(1.0d0,g_sda,alo,ahi,pifac,g_fca,blo,bhi,
     &                         g_sda,alo,ahi)
            alo(3) = alo(3)+1
            ahi(3) = ahi(3)+1
         enddo
      enddo
c
c     Calculate the tensor contributions to the SD-SD term
c
      do ixy = 1, ipairs
c        i = unique(translate(ixy))
c        j = unique(translate(ixy+ipairs))
         i = translate(ixy)
         j = translate(ixy+ipairs)
c        RAW contributions for reference
c        do ix = 1, 3
c           do iy = 1, 3 
c              alo(3) = i_pert + (i-1)*6 + ilist(ix,iy)
c              ahi(3) = i_pert + (i-1)*6 + ilist(ix,iy)
c              blo(3) = (j-1)*6 + ilist(ix,iy)
c              bhi(3) = (j-1)*6 + ilist(ix,iy)
c              valuea = nga_ddot_patch(g_sda,'n',blo,bhi,
c    &                                 g_d1(1),'n',alo,ahi)
c              valueb = nga_ddot_patch(g_sda,'n',blo,bhi,
c    &                                 g_d1(2),'n',alo,ahi)
c              if (ga_nodeid().eq.0) 
c    &            write(6,'(A,4I4,5F12.6)') 'SD raw i j com',
c    &            pairlist(ixy),pairlist(ixy+ipairs),ix,iy,
c    &            (valuea-valueb)*nsd*0.5d0,
c    &             valuea,valueb,nsd
c           enddo
c        enddo
         do ix = 1, 3
            do iy = 1, 3 
               valuea = 0.0d0
               valueb = 0.0d0
               do iz = 1, 3
                  alo(3) = i_pert + (i-1)*6 + ilist(iz,iy)
                  ahi(3) = i_pert + (i-1)*6 + ilist(iz,iy)
                  blo(3) = (j-1)*6 + ilist(iz,ix)
                  bhi(3) = (j-1)*6 + ilist(iz,ix)
                  valuea = valuea + nga_ddot_patch(g_sda,'n',blo,bhi,
     &                                             g_d1(1),'n',alo,ahi)
                  valueb = valueb + nga_ddot_patch(g_sda,'n',blo,bhi,
     &                                             g_d1(2),'n',alo,ahi)
               enddo
               tensor(iy,ix,2,ixy)=(valuea-valueb)*nsd*0.5d0
            enddo
         enddo
      enddo
c
c     The final term to be calculated is the SD-FC tensorterm , which consists
c     of two contributions:
c     D1(SD)_N * FC_M + D1(SD)_M * FC_N  N=perturbing, M=responding
c     We have all the contributions for the first term in the matrices, 
c     hence lets do the ddots
c
      do ixy = 1, ipairs
c        i = unique(translate(ixy))
c        j = unique(translate(ixy+ipairs))
         i = translate(ixy)
         j = translate(ixy+ipairs)
         do ix = 1, 3
            do iy = 1, 3 
               alo(3) = i_pert + (i-1)*6 + ilist(ix,iy)
               ahi(3) = i_pert + (i-1)*6 + ilist(ix,iy)
               blo(3) = j
               bhi(3) = j
               valuea = nga_ddot_patch(g_fca,'n',blo,bhi,
     &                                 g_d1(1),'n',alo,ahi)
               valueb = nga_ddot_patch(g_fca,'n',blo,bhi,
     &                                 g_d1(2),'n',alo,ahi)
               tensor(ix,iy,3,ixy)=(valuea-valueb)*0.5d0
c              if (ga_nodeid().eq.0) 
c    &            write(6,'(A,4I4,3F12.6)') 'SDFC NM tens i j com',
c    &            pairlist(ixy),pairlist(ixy+ipairs),
c    &            ix,iy,tensor(ix,iy,3,ixy)*nsdfc,valuea,valueb
            enddo
         enddo
c        if (ga_nodeid().eq.0) print*,''
      enddo
c
c     For the D1(SD)_M * FC_N  N=perturbing, M=responding we need to calculate the
c     fc integrals for the perturbing atoms. We already calculated the D1(SD) for
c     the responding atoms as extra terms.
c
      call ga_zero(g_fca)
      call int_giao_1ega(basis,basis,g_fca,'fc',coords,
     &                   i_pert,oskel)
      do ixy = 1, ipairs
c        i = unique(translate(ixy))
c        j = unique(translate(ixy+ipairs))
         i = translate(ixy)
         j = translate(ixy+ipairs)
         do ix = 1, 3
            do iy = 1, 3 
               valuea = 0.0d0
               valueb = 0.0d0
               alo(3) = 7*i_pert + (j-1)*6 + ilist(ix,iy)
               ahi(3) = 7*i_pert + (j-1)*6 + ilist(ix,iy)
               blo(3) = i
               bhi(3) = i
               valuea = valuea + nga_ddot_patch(g_fca,'n',blo,bhi,
     &                                          g_d1(1),'n',alo,ahi)
               valueb = valueb + nga_ddot_patch(g_fca,'n',blo,bhi,
     &                                          g_d1(2),'n',alo,ahi)
               tensor(ix,iy,3,ixy)=tensor(ix,iy,3,ixy)+
     &                             (valuea-valueb)*0.5d0
               tensor(ix,iy,3,ixy)=tensor(ix,iy,3,ixy)*nsdfc
c              if (ga_nodeid().eq.0) 
c    &            write(6,'(A,4I4,4F12.6)') 'SDFC MN tens i j com',
c    &            pairlist(ixy),pairlist(ixy+ipairs),
c    &            ix,iy,tensor(ix,iy,3,ixy),(valuea-valueb)*0.5d0*nsdfc,
c    &             valuea,valueb
            enddo
         enddo
c        if (ga_nodeid().eq.0) print*,''
      enddo
c
      call ga_sync()
c
c     Clean up all remaining memory
c
      if (.not.ga_destroy(g_d1(1))) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_d1(1)',0,GA_ERR)
      if (.not.ga_destroy(g_d1(2))) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_d1(2)',0,GA_ERR)
      if (.not.ga_destroy(g_sda)) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_sda',0,GA_ERR)
      if (.not.ga_destroy(g_fca)) call 
     &    errquit('hnd_sdfc: ga_destroy failed g_fca',0,GA_ERR)
c
      call schwarz_tidy()
      call int_terminate()
c
      return
      end
