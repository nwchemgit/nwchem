c
c  Form gradient via finite difference (for debugging only)
c
c
       subroutine mcscf_fdiff_grad( geom, basis, nbf, nclosed, nact,
     $                              oskel, tol2e, dm1, dm2, g_movecs,
     $                              g_coul, g_grad )
       implicit none
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
       integer npt
       parameter(npt=2)
       integer g_kvec, g_tmpmo, g_exch
       integer vlen, nvir, voff, mo_lo, mo_hi, i, j, ij, ipt
       double precision eone, etwo, energy0
       double precision g1
       double precision del(4), ee(4)
       data del/0.001d0,-0.001d0,0.01d0,-0.01d0/
c
c
c
       if (.not.ga_duplicate(g_grad,g_kvec,' crap '))
     $   call errquit('mcscf: cannot dup grad vector')
       if (.not.ga_duplicate(g_movecs,g_tmpmo,' MO crap '))
     $   call errquit('mcscf: cannot dup grad vector')
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
     $                    oskel, tol2e, dm1, dm2,
     $                    g_movecs, g_coul, eone, etwo )
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
             call rohf_k2cf( basis, nbf, nclosed, nact, g_kvec,
     $                       g_movecs, g_tmpmo )
             call moints_build_6x( basis, mo_lo, mo_lo, mo_hi, 1, nbf,
     $                             g_tmpmo, g_coul, .true.,
     $                             g_exch, .false., 16 )
             call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          oskel, tol2e, dm1, dm2,
     $                          g_tmpmo, g_coul, eone, etwo )
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
             call rohf_k2cf( basis, nbf, nclosed, nact, g_kvec,
     $                       g_movecs, g_tmpmo )
             call moints_build_6x( basis, mo_lo, mo_lo, mo_hi, 1, nbf,
     $                             g_tmpmo, g_coul, .true.,
     $                             g_exch, .false., 16 )
             call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          oskel, tol2e, dm1, dm2,
     $                          g_tmpmo, g_coul, eone, etwo )
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
       call moints_build_6x( basis, mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .false., 16 )
c
c
c
       if (.not. ga_destroy(g_tmpmo))
     $   call errquit('mcscf: cannot destroy ga')
       if (.not. ga_destroy(g_kvec))
     $   call errquit('mcscf: cannot destroy ga')
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
     $                              g_coul, g_k1, hess )
       implicit none
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
       double precision hess(*)                        ! [output] Hessian
c
c
c
       integer npt
       parameter(npt=2)
       integer g_k2, g_tmp1, g_tmp2, g_exch
       integer vlen, nvir, voff, mo_lo, mo_hi
       integer i, j, ij, ji, ipt, jpt
       integer vlo, vhi
       double precision eone, e2ii, e2ia, e2aa, h
       double precision del(4), ee(4,4)
       double precision ee2ii(4,4), ee2ia(4,4), ee2aa(4,4)
       data del/0.001d0,-0.001d0,0.01d0,-0.01d0/
c
c
c
       if (.not.ga_duplicate(g_k1,g_k2,' crap '))
     $   call errquit('mcscf: cannot dup grad vector')
       if (.not.ga_duplicate(g_movecs,g_tmp1,' MO crap '))
     $   call errquit('mcscf: cannot dup grad vector')
       if (.not.ga_duplicate(g_movecs,g_tmp2,' MO crap 1'))
     $   call errquit('mcscf: cannot dup grad vector')
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       nvir = nbf - nact - nclosed
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact
       g_exch = 0
       vlo = 1
       vhi = vlen
c
c
c
       PRINT*,' Finite difference Hessian start '
       call dfill((vlen*vlen),0.d0,hess,1)
       do i=vlo,vhi
         PRINT*,'I ',I
         do j=1,i-1
           ij = (i-1)*vlen + j
           ji = (j-1)*vlen + i
           do ipt=1,npt
             call ga_zero(g_k1)
             call ga_put(g_k1,i,i,1,1,del(ipt),1)
             call rohf_k2cf( basis, nbf, nclosed, nact, g_k1,
     $                     g_movecs, g_tmp1 )
             do jpt=1,npt
               call ga_zero(g_k2)
               call ga_put(g_k2,j,j,1,1,del(jpt),1)
               call rohf_k2cf( basis, nbf, nclosed, nact, g_k2,
     $                         g_tmp1, g_tmp2 )
               call moints_build_2x( basis, .true.,
     $                               mo_lo, mo_lo, mo_hi, 1, nbf,
     $                               g_tmp2, g_coul, .true.,
     $                               g_exch, .false., 16 )
               call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                            oskel, tol2e, dm1, dm2,
     $                            g_tmp2, g_coul, eone, e2ii,
     $                            e2ia, e2aa )
               ee(ipt,jpt) = eone + e2ii + e2ia + e2aa
               ee2ii(ipt,jpt) = e2ii
               ee2ia(ipt,jpt) = e2ia
               ee2aa(ipt,jpt) = e2aa
             enddo
           enddo
           h = (ee(1,1)-ee(1,2)-ee(2,1)+ee(2,2))/(4.d0*del(1)*del(1))
           hess(ij) = h
           hess(ji) = h
         enddo
       enddo
c
c  Restore old MO integrals
c
       call moints_build_6x( basis, mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .false., 16 )
c
c  Write Hessian to file
c
       open(unit=90,file='hessFF',form='unformatted',status='unknown')
       do i=1,vlen
         ij = (i-1)*vlen
         write(90) (hess(ij+j),j=1,vlen)
       enddo
       close(90)
c
c
c
       if (.not. ga_destroy(g_tmp1))
     $   call errquit('mcscf: cannot destroy ga')
       if (.not. ga_destroy(g_tmp2))
     $   call errquit('mcscf: cannot destroy ga')
       if (.not. ga_destroy(g_k2))
     $   call errquit('mcscf: cannot destroy ga')
       return
       end











