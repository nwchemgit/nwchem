       subroutine mcscf_hdiag( nbf, nclosed, nact, dm1, dm2,
     $                         g_coul, g_exch, g_ifock,
     $                         g_afock, g_gfock, g_hdiag )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
       integer nbf, nclosed, nact
       double precision dm1(nact,nact)
       double precision dm2(nact,nact,nact,nact)
       integer g_coul, g_exch
       integer g_ifock, g_afock, g_gfock
       integer g_hdiag
c
       integer nvir, i, j, ioff, cbase, abase, vbase
       integer jjj, ii, ij
       integer numnodes, myid
       double precision xx, yy, zz
       integer l_di, k_di
       integer l_da, k_da
       integer l_dg, k_dg
       integer l_tmp, k_tmp

       INTEGER L_TMP2, K_TMP2, VLEN
       double precision di, da, dg
       di(i) = dbl_mb(k_di + i - 1)
       da(i) = dbl_mb(k_da + i - 1)
       dg(i) = dbl_mb(k_dg + i - 1)
c     
       nvir = nbf - nclosed - nact
       cbase = 0 
       abase = nclosed
       vbase = nclosed + nact
       numnodes = ga_nnodes()
       myid = ga_nodeid()
c     
c     Get diagonals of each matrix
c
       if (.not.ma_push_get(MT_DBL,nbf,'mcscf_hdiag: tmp',l_di,k_di))
     $      call errquit('rohf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,nbf,'mcscf_hdiag: tmp',l_da,k_da))
     $      call errquit('rohf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,nbf,'mcscf_hdiag: tmp',l_dg,k_dg))
     $      call errquit('rohf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,nbf,'mcscf_hdiag: tmp',l_tmp,k_tmp))
     $      call errquit('rohf_hdiag_scale: cannot allocate',0)
c     
       call ga_get_diagonal(g_ifock, dbl_mb(k_di))
       call ga_get_diagonal(g_afock, dbl_mb(k_da))
       call ga_get_diagonal(g_gfock, dbl_mb(k_dg))
c
c  Inactive-virtual piece
c
       do i = myid+1, nclosed, numnodes
         ioff = (i-1)*nvir + 1
         do j=1,nvir
           dbl_mb(k_tmp+j-1) = 4.d0*(di(j+vbase)+da(j+vbase)
     $                              -di(i+cbase)-da(i+cbase))
         enddo
         call ga_put(g_hdiag,ioff,ioff+nvir-1,1,1,dbl_mb(k_tmp),1)
       enddo
c     
c  Active-virtual piece
c     
       do i = myid+1, nact, numnodes
         ioff = (nclosed+i-1)*nvir + 1
         do j=1,nvir
           dbl_mb(k_tmp+j-1) = 2.d0*(dm1(i,i)*di(j+vbase)-dg(i+abase))
         enddo
         call ga_put(g_hdiag, ioff, ioff+nvir-1, 1, 1, dbl_mb(k_tmp),1)
       enddo
c     
c  Inactive-active piece
c     
       do i = myid+1, nclosed, numnodes
         ioff = (nclosed+nact)*nvir + (i-1)*nact + 1
         do j=1,nact
           dbl_mb(k_tmp+j-1) = 2.d0*(dm1(j,j)*di(i+cbase) - dg(j+abase)
     $                       + 2.d0*(di(j+abase)+da(j+abase))
     $                       - 2.d0*(di(i+cbase)+da(i+cbase)))
         enddo
         call ga_put(g_hdiag, ioff, ioff+nact-1, 1, 1, dbl_mb(k_tmp),1)
       enddo
c
c
C       CALL GA_ROWPRINT('1e Hessian diagonal',G_HDIAG)
c
c  Part of AV piece from 2-pdm contraction against
c  MO integerals, is approximated by
c         
c         gamma  sum [ (vv|aa) - (va|va)]
c              tt   v
c
c  putting gamma_tt = 1 reduces to ROHF expression
c  Compute 
c         z(a) = sum [ (vv|aa) - (va|va) ]
c                   v
c  (this intermediate term can potentially be precomputed
c  and stored in the unused blocks of the gen. Fock
c  inside the Fock build)
c
       do j=1,nvir
         jjj = (vbase+j-1)*nbf + (vbase+j)
         zz = 0.d0
         do i=1,nact
           ii = ((i-1)*i)/2 + i
           call ga_get(g_coul,jjj,jjj,ii,ii,xx,1)
           call ga_get(g_exch,jjj,jjj,ii,ii,yy,1)
           zz = zz + xx - yy
         enddo
         dbl_mb(k_tmp+j-1) = zz
       enddo
c
c  Now contract against density
c
       ioff = nclosed*nvir
       do i=myid+1, nact, numnodes
         do j=1,nvir
           ij = ioff + (i-1)*nvir + j
           xx = dm1(i,i)*dbl_mb(k_tmp+j-1)*2.d0
           call ga_acc(g_hdiag,ij,ij,1,1,xx,1,1.d0)
         enddo
       enddo
c
c
c  Part of CA term from 2-pdm contraction
c  (equivalent to above), also approximate
c
c       gamma  sum [ (vv|ii) - (vi|vi) ]
c            tt   v
c
c  This is exact for ROHF
c
       do j=1,nclosed
         jjj = (cbase+j-1)*nbf + (cbase+j)
         zz = 0.d0
         do i=1,nact
           ii = ((i-1)*i)/2 + i
           call ga_get(g_coul,jjj,jjj,ii,ii,xx,1)
           call ga_get(g_exch,jjj,jjj,ii,ii,yy,1)
           zz = zz + xx - yy
         enddo
         dbl_mb(k_tmp+j-1) = zz
       enddo
c
c  Now contract against density
c
       ioff = nvir*(nclosed+nact)
       do i=myid+1, nact, numnodes
         do j=1,nclosed
           ij = ioff + (j-1)*nact + i
           xx = dm1(i,i)*dbl_mb(k_tmp+j-1)*2.d0
           call ga_acc(g_hdiag,ij,ij,1,1,xx,1,1.d0)
         enddo
       enddo
c
c  Clean up
c
       if (.not.ma_pop_stack(l_tmp))
     $      call errquit('mcscf_hdiag: cannot pop stack',0)
       if (.not.ma_pop_stack(l_dg))
     $      call errquit('mcscf_hdiag: cannot pop stack',0)
       if (.not.ma_pop_stack(l_da))
     $      call errquit('mcscf_hdiag: cannot pop stack',0)
       if (.not.ma_pop_stack(l_di))
     $      call errquit('mcscf_hdiag: cannot pop stack',0)
c
c
c
c$$$       CALL GA_ROWPRINT('1e Hessian diagonal',G_HDIAG)
c
c
c
       return
       end




       subroutine mcscf_hdiag_exact( geom, basis, nbf, nclosed, nact,
     $                               oskel, pflg, lshift, dm1, dm2,
     $                               g_coul, g_exch,
     $                               g_movecs, g_ifock, g_afock,
     $                               g_gfock, g_hdiag )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
       integer geom, basis
       integer nbf, nclosed, nact
       logical oskel
       integer pflg
       double precision lshift
       double precision dm1(nact,nact)
       double precision dm2(nact,nact,nact,nact)
       integer g_coul, g_exch, g_movecs
       integer g_ifock, g_afock, g_gfock
       integer g_hdiag
c
       integer nvir, i, vlen
       integer g_x, g_prod
       logical status
       double precision xx,yy
c     
       nvir = nbf - nclosed - nact
       vlen = nvir*(nclosed+nact) + nact*nclosed
       status = ga_duplicate(g_hdiag,g_x,'crap I')
       status = ga_duplicate(g_hdiag,g_prod,'crap II')
       xx = 1.d0
       do i=1,vlen
         call ga_zero(g_x)
         call ga_put(g_x,i,i,1,1,xx,1)
         call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                     oskel, pflg, lshift, dm1, dm2,
     $                     g_movecs, g_ifock, g_afock, g_gfock,
     $                     g_coul, g_exch, g_x, g_prod )
         call ga_get(g_prod,i,i,1,1,yy,1)
         call ga_put(g_hdiag,i,i,1,1,yy,1)
       enddo
       status = ga_destroy(g_prod)
       status = ga_destroy(g_x)
       return
       end
