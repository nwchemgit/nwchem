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








       subroutine mcscf_hdiag_scale( nbf, nclosed, nact, lshift, dm1,
     $                               g_ifock, g_afock, g_gfock,
     $                               g_hdiag, g_x )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
       integer nbf, nclosed, nact
       double precision lshift
       double precision dm1(nact,nact)
       integer g_ifock, g_afock, g_gfock
       integer g_hdiag
       integer g_x
c
c ********** NOT correct for MCSCF, based on ROHF equivalent **********
c Scale the vector x by the ROHF 1e-Hessian diagonal
c
c             x                    x
c              pq                   ia
c   x   =  --------           = ---------    [ closed-virtual ]
c    pq      (1e)                 cv    cv
c           H                   4F  - 4F
c            pq,pq                aa    ii
c
c                                  x
c                                   pa
c                             = ---------    [ partial-virtual ]
c                                 pv    pv
c                               2F  - 2F
c                                 aa   pp
c
c                                  x
c                                   ip
c                             = ---------    [ closed-partial ]
c                                 cp    cp
c                               2F  - 2F
c                                 pp    ii
c
c
c
c
       integer nvir, voff, ooff, xoff, vlen
       integer i, j, ii, jj, ij
       double precision xxi, xx, denom, yy
       integer l_di, k_di
       integer l_da, k_da
       integer l_dg, k_dg
       integer l_v, k_v
       integer numnodes, myid
c
       double precision di, da, dg
       di(i) = dbl_mb(k_di + i - 1)
       da(i) = dbl_mb(k_da + i - 1)
       dg(i) = dbl_mb(k_dg + i - 1)
c
c
       nvir = nbf - nclosed - nact
       voff = nclosed + nact + 1
       ooff = nclosed + 1
       vlen = (nclosed+nact)*nvir + nclosed*nact
       myid = ga_nodeid()
       numnodes = ga_nnodes()
c
c
c
       do i=myid+1,vlen,numnodes
         call ga_get(g_x,i,i,1,1,xx,1)
         call ga_get(g_hdiag,i,i,1,1,yy,1)
C         if (abs(yy).lt.1.d-8) yy=1.d-8
         xx = xx/(yy + lshift)
         call ga_put(g_x,i,i,1,1,xx,1)
       enddo
       goto 1000                !  Hessian diag is stored
c
c
c
       if (.not.ma_push_get(MT_DBL,nbf,'Fi diag',l_di,k_di))
     $      call errquit('mcscf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,nbf,'Fa diag',l_da,k_da))
     $      call errquit('mcscf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,nbf,'Fg diag',l_dg,k_dg))
     $      call errquit('mcscf_hdiag_scale: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL,vlen,'Fg diag',l_v,k_v))
     $      call errquit('mcscf_hdiag_scale: cannot allocate',0)
c
c
c
       do i=1,nbf
         call ga_get(g_ifock,i,i,i,i,dbl_mb(k_di+i-1),1)
         call ga_get(g_afock,i,i,i,i,dbl_mb(k_da+i-1),1)
         call ga_get(g_gfock,i,i,i,i,dbl_mb(k_dg+i-1),1)
       enddo
c
c
c
       do i=1,nclosed
         xxi = di(i) + da(i)
         do j=1,nvir
           jj = voff + j - 1
           ij = (i-1)*nvir + j
           call ga_get(g_x,ij,ij,1,1,xx,1)
           dbl_mb(k_v+ij-1) = 4.d0*(di(jj)+da(jj)-xxi)
           xx = xx/(4.d0*(di(jj)+da(jj)-xxi)+lshift)
           call ga_put(g_x,ij,ij,1,1,xx,1)
         enddo
       enddo

       write(6,901) (di(voff+j-1),j=1,nvir)
       xoff = nclosed*nvir
       do i=1,nact
         ii = ooff + i - 1
         do j=1,nvir
           jj = voff + j - 1
           ij = xoff + (i-1)*nvir + j
           call ga_get(g_x,ij,ij,1,1,xx,1)
           dbl_mb(k_v+ij-1) = 2.d0*(dm1(i,i)*di(jj)-dg(ii))
           xx = xx/(2.d0*(dm1(i,i)*di(jj)-dg(ii)) + lshift)
           call ga_put(g_x,ij,ij,1,1,xx,1)
         enddo
       enddo

       xoff = xoff + nact*nvir
       do i=1,nclosed
         xxi = di(i) + da(i)
         do j=1,nact
           jj = ooff + j - 1
           ij = xoff + (i-1)*nact + j
           call ga_get(g_x,ij,ij,1,1,xx,1)
           dbl_mb(k_v+ij-1) = 2.d0*(dm1(j,j)*di(i) - dg(jj) +
     $                        2.d0*(di(jj)+da(jj)) - 2.d0*xxi)
           xx = xx/(2.d0*(dm1(j,j)*di(i) - dg(jj) +
     $             2.d0*(di(jj)+da(jj)) - 2.d0*xxi) + lshift)
           call ga_put(g_x,ij,ij,1,1,xx,1)
         enddo
       enddo

       WRITE(6,900)
 900   FORMAT('1E HESSIAN DIAGONAL')
       WRITE(6,901) (DBL_MB(K_V+I-1),I=1,VLEN)
 901   FORMAT(10F12.6)
c
c
c
10000  continue
       if (.not.ma_pop_stack(l_v))
     $    call errquit('mcscf_hdiag_scale: cannot pop stack',0)
       if (.not.ma_pop_stack(l_dg))
     $    call errquit('mcscf_hdiag_scale: cannot pop stack',0)
       if (.not.ma_pop_stack(l_da))
     $    call errquit('mcscf_hdiag_scale: cannot pop stack',0)
       if (.not.ma_pop_stack(l_di))
     $    call errquit('mcscf_hdiag_scale: cannot pop stack',0)

 1000  continue
c

       return
       end







       subroutine mcscf_hdiag_scale00( nbf, nclosed, nact, lshift, dm1,
     $                               g_ifock, g_afock, g_gfock,
     $                               g_hdiag, g_x )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
c
       integer nbf, nclosed, nact
       double precision lshift
       double precision dm1(nact,nact)
       integer g_ifock, g_afock, g_gfock
       integer g_hdiag
       integer g_x
c
c ********** NOT correct for MCSCF, based on ROHF equivalent **********
c Scale the vector x by the ROHF 1e-Hessian diagonal
c
c             x                    x
c              pq                   ia
c   x   =  --------           = ---------    [ closed-virtual ]
c    pq      (1e)                 cv    cv
c           H                   4F  - 4F
c            pq,pq                aa    ii
c
c                                  x
c                                   pa
c                             = ---------    [ partial-virtual ]
c                                 pv    pv
c                               2F  - 2F
c                                 aa   pp
c
c                                  x
c                                   ip
c                             = ---------    [ closed-partial ]
c                                 cp    cp
c                               2F  - 2F
c                                 pp    ii
c
c
c
c
       integer nvir, voff, ooff, xoff, vlen, xend
       integer i, j, ii, jj, ij
       double precision xxi, xx, denom, yy
       integer l_di, k_di
       integer l_da, k_da
       integer l_dg, k_dg
       integer l_v, k_v
c
       double precision di, da, dg
       di(i) = dbl_mb(k_di + i - 1)
       da(i) = dbl_mb(k_da + i - 1)
       dg(i) = dbl_mb(k_dg + i - 1)
c
c
       nvir = nbf - nclosed - nact
       voff = nclosed + nact + 1
       ooff = nclosed + 1
       vlen = (nclosed+nact)*nvir + nclosed*nact
c
c  CV bit
c
       xend = nvir*nclosed
       do i=1,xend
         call ga_get(g_x,i,i,1,1,xx,1)
         call ga_get(g_hdiag,i,i,1,1,yy,1)
         xx = xx/(yy + lshift)
         call ga_put(g_x,i,i,1,1,xx,1)
       enddo
c
c  AV bit
c     
       xoff = nvir*nclosed+1
       xend = nvir*(nclosed+nact)
       do i=xoff,xend
         call ga_get(g_x,i,i,1,1,xx,1)
         call ga_get(g_hdiag,i,i,1,1,yy,1)
         xx = xx/(yy + lshift + 0.5d0)
         call ga_put(g_x,i,i,1,1,xx,1)
       enddo
c
c  CA bit
c
       xoff = nvir*(nclosed+nact) + 1
       do i=xoff,vlen
         call ga_get(g_x,i,i,1,1,xx,1)
         call ga_get(g_hdiag,i,i,1,1,yy,1)
         xx = xx/(yy + lshift)
         call ga_put(g_x,i,i,1,1,xx,1)
       enddo
       return
       end

