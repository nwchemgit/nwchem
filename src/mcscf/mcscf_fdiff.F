c
c  Form gradient via finite difference (for debugging only)
c
c
       subroutine mcscf_fdiff_grad( geom, basis, nbf, nclosed, nact,
     $                              oskel, tol2e, dm1, dm2, g_movecs,
     $                              g_coul, g_grad )
*
* $Id$
*
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
       integer geom, basis                            ! [input] Geometry and basis handles
       integer nbf                                    ! [input] Number of basis functions
       integer nclosed                                ! [input] Number of closed shells
       integer nact                                   ! [input] Number of open shells
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       double precision dm1(nact,nact)                ! [input] Active space 1PDM
       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
       integer g_movecs                               ! [input] MO coefficients
       integer g_coul                                 ! [input] Coulomb integrals
       integer g_grad                                 ! [output] Gradient
c
c
c     
       logical osym
       integer npt
       parameter(npt=2)
       integer g_kvec, g_tmpmo, g_exch
       integer vlen, nvir, voff, mo_lo, mo_hi, i, j, ij, ipt
       double precision eone, etwo, ecore, energy0
       double precision g1
       double precision del(4), ee(4)
       data del/0.001d0,-0.001d0,0.01d0,-0.01d0/
       data osym/.false./
c
c
c
       if (.not.ga_duplicate(g_grad,g_kvec,' crap '))
     $   call errquit('mcscf: cannot dup grad vector',0, GA_ERR)
       if (.not.ga_duplicate(g_movecs,g_tmpmo,' MO crap '))
     $   call errquit('mcscf: cannot dup grad vector',0, GA_ERR)
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       nvir = nbf - nact - nclosed
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact
       g_exch = 0
c
c
c
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    .true., oskel, tol2e, dm1, dm2,
     $                    g_movecs, g_coul, eone, etwo, ecore )
       energy0 = eone + etwo
c
c
c
       do i=1,(nclosed+nact)
         do j=1,nvir
           ij = (i-1)*nvir + j
           do ipt=1,npt
             call ga_zero(g_kvec)
             call ga_put(g_kvec,ij,ij,1,1,del(ipt),1)
             call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_kvec,
     $                       g_movecs, g_tmpmo )
             call moints_build_6x( basis, osym,
     $                             mo_lo, mo_lo, mo_hi, 1, nbf,
     $                             g_tmpmo, g_coul, .true.,
     $                             g_exch, .false., 16, .false.)
             call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          .true., oskel, tol2e, dm1, dm2,
     $                          g_tmpmo, g_coul, eone, etwo, ecore )
             ee(ipt) = eone + etwo
           enddo
           g1 = (ee(1) - ee(2))/(2.d0*del(1))
           call ga_put(g_grad,ij,ij,1,1,g1,1)
c$$$           write(6,901) i,(voff+j),ij,g1
 901       format(5x,2i4,2x,i5,2(5x,f14.10),5x,e14.6)
         enddo
       enddo
       do i=1,nclosed
         do j=1,nact
           ij = (nclosed+nact)*nvir+(i-1)*nact + j
           do ipt=1,npt
             call ga_zero(g_kvec)
             call ga_put(g_kvec,ij,ij,1,1,del(ipt),1)
             call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_kvec,
     $                       g_movecs, g_tmpmo )
             call moints_build_6x( basis, osym,
     $                             mo_lo, mo_lo, mo_hi, 1, nbf,
     $                             g_tmpmo, g_coul, .true.,
     $                             g_exch, .false., 16, .false. )
             call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          .true., oskel, tol2e, dm1, dm2,
     $                          g_tmpmo, g_coul, eone, etwo, ecore )
             ee(ipt) = eone + etwo
           enddo
           g1 = (ee(1) - ee(2))/(2.d0*del(1))
           call ga_put(g_grad,ij,ij,1,1,g1,1)
c$$$           write(6,901) i,(nclosed+j),ij,g1
         enddo
       enddo
c
c  Restore old MO integrals
c
       call moints_build_6x( basis, osym,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .false., 16, .false. )
c
c
c
       if (.not. ga_destroy(g_tmpmo))
     $   call errquit('mcscf: cannot destroy ga',0, GA_ERR)
       if (.not. ga_destroy(g_kvec))
     $   call errquit('mcscf: cannot destroy ga',0, GA_ERR)
       return
       end


c
c
c
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c  Form Hessian via finite difference (for debugging only)
c
c
       subroutine mcscf_fdiff_hess( geom, basis, nbf, nclosed, nact,
     $                              oskel, tol2e, dm1, dm2, g_movecs,
     $                              g_coul, g_k1 )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
       integer geom, basis                            ! [input] Geometry and basis handles
       integer nbf                                    ! [input] Number of basis functions
       integer nclosed                                ! [input] Number of closed shells
       integer nact                                   ! [input] Number of open shells
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       double precision dm1(nact,nact)                ! [input] Active space 1PDM
       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
       integer g_movecs                               ! [input] MO coefficients
       integer g_coul                                 ! [input] Coulomb integrals
       integer g_k1                                   ! [output] KVec
c
c
c
       logical osym
       integer npt
       parameter(npt=2)
       integer g_k2, g_tmp1, g_tmp2, g_exch
       integer vlen, nvir, voff, mo_lo, mo_hi
       integer i, j, ij, ji, ipt, jpt
       integer vlo, vhi
       integer l_v, k_v
       double precision eone, etwo, h, ecore
       double precision del(4), ee(4,4)
       data osym/.false./
       data del/0.001d0,-0.001d0,0.01d0,-0.01d0/
c
c
c
       if (.not.ga_duplicate(g_k1,g_k2,' crap '))
     $   call errquit('mcscf: cannot dup grad vector',0, GA_ERR)
       if (.not.ga_duplicate(g_movecs,g_tmp1,' MO crap '))
     $   call errquit('mcscf: cannot dup grad vector',0, GA_ERR)
       if (.not.ga_duplicate(g_movecs,g_tmp2,' MO crap 1'))
     $   call errquit('mcscf: cannot dup grad vector',0, GA_ERR)
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       nvir = nbf - nact - nclosed
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact
       g_exch = 0
       vlo = 1
       vhi = vlen
       if (.not.ma_push_get(MT_DBL, vlen, 'h', l_v, k_v))
     $      call errquit('mcscf_fdiff: cannot allocate',0, MA_ERR)
c
       open(unit=90,file='hessFF',form='unformatted',status='unknown')
       write(90) nclosed,nact,nvir
       write(90) (nclosed*nvir),((nclosed+nact)*nvir),vlen
c
       PRINT*,' Finite difference Hessian start '
       do i=1,vlen
         write(6,711) i
 711     format('Row:',i4)
         do j=1,i
           ij = (i-1)*vlen + j
           ji = (j-1)*vlen + i
           do ipt=1,npt
             call ga_zero(g_k1)
             call ga_put(g_k1,i,i,1,1,del(ipt),1)
             call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_k1,
     $                     g_movecs, g_tmp1 )
             do jpt=1,npt
               call ga_zero(g_k2)
               call ga_put(g_k2,j,j,1,1,del(jpt),1)
               call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_k2,
     $                         g_tmp1, g_tmp2 )
               call moints_build_6x( basis, osym,
     $                               mo_lo, mo_lo, mo_hi, 1, nbf,
     $                               g_tmp2, g_coul, .true.,
     $                               g_exch, .false., 16, .false. )
               call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                            .true., oskel, tol2e, dm1, dm2,
     $                            g_tmp2, g_coul, eone, etwo, ecore )
               ee(ipt,jpt) = eone + etwo
             enddo
           enddo
           h = (ee(1,1)-ee(1,2)-ee(2,1)+ee(2,2))/(4.d0*del(1)*del(1))
           dbl_mb(k_v+j-1) = h
         enddo
         call hv_writev(90,i,dbl_mb(k_v))
       enddo
       close(90)
c
c  Restore old MO integrals
c
       call moints_build_6x( basis, osym,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .false., 16, .false. )
c
c
c
       if (.not.ma_pop_stack(l_v))
     $      call errquit('mcscf_fdiff: cannot pop stack?',0, MA_ERR)
       if (.not. ga_destroy(g_tmp1))
     $   call errquit('mcscf: cannot destroy ga',0, GA_ERR)
       if (.not. ga_destroy(g_tmp2))
     $   call errquit('mcscf: cannot destroy ga',0, GA_ERR)
       if (.not. ga_destroy(g_k2))
     $   call errquit('mcscf: cannot destroy ga',0, GA_ERR)
       return
       end










