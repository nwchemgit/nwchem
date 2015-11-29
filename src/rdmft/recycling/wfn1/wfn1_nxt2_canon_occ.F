C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a new state canonicalizing everything
C>
C> As with regular effective 1-electron models there are redundant
C> rotations in this model also. For the occupation functions those
C> redundant rotations are easily identified as mixings in the 
C> occupied-occupied block and the virtual-virtual block. For the
C> natural orbitals this distinction is not so clear cut as they have
C> fractional occupations. Hence we construct a regular Fock matrix
C> and diagonalize that to canonicalize the natural orbitals.
C>
C> A problem with canonicalizing the orbitals in this setting is 
C> that because of the degeneracies due to Gilbert's theorem even the
C> canonical orbitals may be arbitrarily mixed. A double
C> canonicalization may be needed to fix this. 
C>
      subroutine wfn1_nxt2_canon_all(nbf,nea,neb,statein,fock,stateout)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "wfn1_nxt.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision statein(8*nbf*nbf+2*nbf) !< [Input] The initial state
      double precision fock(8*nbf*nbf)    !< [Input] The fock matrices
c
      double precision stateout(8*nbf*nbf+2*nbf) !< [Input] The final state
c
      double precision tt
      double precision levelshift
      parameter(levelshift = -1.0d0)
c
      integer info  !< The status flag for DSYEV
      integer lwork !< The work space size for DSYEV
      integer iocn  !< Offset short hand for occupation numbers
      integer iocf  !< Offset short hand for occupation functions
      integer icnf  !< Offset short hand for electron configurations
c
      integer l_m1,   k_m1
      integer l_m2,   k_m2
      integer l_exp,  k_exp
      integer l_occa, k_occa
      integer l_occb, k_occb
      integer l_gg,   k_gg
c
      integer ir !< Counter
      integer iq !< Counter
      integer ii !< Counter
      integer jj !< Counter
c
c     Allocate memory
c
      lwork = max(nbf*nbf,3*nbf)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp1",l_m1,k_m1))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,lwork,"tmp2",l_m2,k_m2))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"exp",l_exp,k_exp))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf,"occa",l_occa,k_occa))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate occa",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occb",l_occb,k_occb))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate occb",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"gershgorin",l_gg,k_gg))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate gg",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
c
c     Compute the occupation numbers
c
      call dfill(nbf,0.0d0,dbl_mb(k_occa),1)
      call dfill(nbf,0.0d0,dbl_mb(k_occb),1)
c
      icnf = wfn1_nxt_aconf()
      iocf = wfn1_nxt_aocc()
      iocn = k_occa
      do ir = 0, nbf-1
        do ii = 0, nbf-1
          dbl_mb(iocn+ii) = dbl_mb(iocn+ii)
     &       + statein(iocf+ii+ir*nbf)*statein(iocf+ii+ir*nbf)
     &         *statein(icnf+ir)
        enddo 
      enddo
c
      icnf = wfn1_nxt_bconf()
      iocf = wfn1_nxt_bocc()
      iocn = k_occb
      do ir = 0, nbf-1
        do ii = 0, nbf-1
          dbl_mb(iocn+ii) = dbl_mb(iocn+ii)
     &       + statein(iocf+ii+ir*nbf)*statein(iocf+ii+ir*nbf)
     &         *statein(icnf+ir)
        enddo 
      enddo
c
c     Initialize state out
c
      call dcopy(wfn1_nxt_size1(),statein,1,stateout,1)
c
c     Alpha natural orbitals (occupied part)
c     --------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_anat()),1,dbl_mb(k_exp),1)
c
c     scale columns
c
      do ii = 0, nbf-1
c       call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_exp+ii*nbf),1)
c       just diagonal
        dbl_mb(k_exp+ii*nbf+ii) = dbl_mb(k_exp+ii*nbf+ii)
     +                          * dbl_mb(k_occa+ii)
      enddo
c
c     scale rows
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occa+ii)),dbl_mb(k_exp+ii),nbf)
c     enddo
c
c     Alpha natural orbitals (virtual part)
c     -------------------------------------
c
c     call dcopy(nbf*nbf,fock(wfn1_nxt_anat()),1,dbl_mb(k_m1),1)
c     do ii = 0, nbf-1
c       dbl_mb(k_occa+ii) = 1.0d0 - dbl_mb(k_occa+ii)
c     enddo
c
c     scale columns
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occa+ii)),dbl_mb(k_m1+ii*nbf),1)
c     enddo
c
c     scale rows
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occa+ii)),dbl_mb(k_m1+ii),nbf)
c     enddo
c
c     Add occupied-occupied and virtual-virtual blocks
c
c     call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m1),1,dbl_mb(k_exp),1)
c
c     Print diagonal elements of Fock matrix
c
      write(*,*)'wfn1_nxt2_canon_all'
      write(*,*)'Alpha natural orbital Fock elements'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+1,dbl_mb(k_occa+ii),
     +                             dbl_mb(k_exp+ii+nbf*ii)
      enddo
c     write(*,*)'Alpha natural orbital original Fock'
c     call hess_hssout(fock(wfn1_nxt_anat()),nbf,nbf,nbf)
c     write(*,*)'Alpha natural orbital occ-occ+virt-virt Fock'
c     call hess_hssout(dbl_mb(k_exp),nbf,nbf,nbf)
c
c     - symmetrize
c
      do jj = 0, nbf-1
        do ii = jj+1, nbf-1
          tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
          dbl_mb(k_exp+ii+jj*nbf) = tt
          dbl_mb(k_exp+jj+ii*nbf) = tt
        enddo
      enddo
c
c     - sort stuff
c
c     call wfn1_nxt2_sort_diag(nbf,dbl_mb(k_exp),
c    &     statein(wfn1_nxt_anat()),statein(wfn1_nxt_aocc()))
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
     &                      dbl_mb(k_occa))
c
c     - establish rotation and apply it to both natural orbitals
c       and occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
cDEBUG
c     write(*,*)'Alpha natural orbital Fock eigenvalues'
c     do ii = 0, nbf-1
c       write(*,*)ii+1,dbl_mb(k_m1+ii)
c     enddo
c     call hess_hssout(dbl_mb(k_exp),nbf,nbf,nbf)
cDEBUG
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed a: ",info,UERR)
      endif
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_anat()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_anat()),nbf)
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_aocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta natural orbitals (occupied part)
c     -------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_bnat()),1,dbl_mb(k_exp),1)
c
c     Scale columns
c
      do ii = 0, nbf-1
c       call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_exp+ii*nbf),1)
c       just diagonal
        dbl_mb(k_exp+ii*nbf+ii) = dbl_mb(k_exp+ii*nbf+ii)
     +                          * dbl_mb(k_occb+ii)
      enddo
c
c     Scale rows
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occb+ii)),dbl_mb(k_exp+ii),nbf)
c     enddo
c
c     Beta natural orbitals (virtual part)
c     ------------------------------------
c
c     call dcopy(nbf*nbf,fock(wfn1_nxt_bnat()),1,dbl_mb(k_m1),1)
c     do ii = 0, nbf-1
c       dbl_mb(k_occb+ii) = 1.0d0 - dbl_mb(k_occb+ii)
c     enddo
c
c     Scale columns
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occb+ii)),dbl_mb(k_m1+ii*nbf),1)
c     enddo
c
c     Scale rows
c
c     do ii = 0, nbf-1
c       call dscal(nbf,sqrt(dbl_mb(k_occb+ii)),dbl_mb(k_m1+ii),nbf)
c     enddo
c
c     Add occupied-occupied and virtual-virtual blocks
c
c     call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m1),1,dbl_mb(k_exp),1)
c
c     Print diagonal elements of Fock matrix
c
      write(*,*)'Beta natural orbital Fock elements'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+1,dbl_mb(k_occb+ii),
     +                             dbl_mb(k_exp+ii+nbf*ii)
      enddo
c
c     - symmetrize
c
      do jj = 0, nbf-1
        do ii = jj+1, nbf-1
          tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
          dbl_mb(k_exp+ii+jj*nbf) = tt
          dbl_mb(k_exp+jj+ii*nbf) = tt
        enddo
      enddo
c
c     - sort stuff
c
c     call wfn1_nxt2_sort_diag(nbf,dbl_mb(k_exp),
c    &     statein(wfn1_nxt_bnat()),statein(wfn1_nxt_bocc()))
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
     &                      dbl_mb(k_occb))
c
c     - establish rotation and apply it to both natural orbitals
c       and occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed c: ",info,UERR)
      endif
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bnat()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bnat()),nbf)
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_bocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Alpha occupation functions (occupied part)
c     ------------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_exp),1)
      call dcopy(nbf,statein(wfn1_nxt_aconf()),1,dbl_mb(k_occa),1)
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_exp+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_exp+ii),nbf)
      enddo
c
c     Alpha occupation functions (virtual part)
c     ------------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_m2),1)
      do ii = 0, nbf-1
        dbl_mb(k_occa+ii) = 1.0d0 - dbl_mb(k_occa+ii)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_m2+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_m2+ii),nbf)
      enddo
      do ii = 0, nbf-1
        dbl_mb(k_occa+ii) = 1.0d0 - dbl_mb(k_occa+ii)
      enddo
c
c     Sum occupied-occupied and virtual-virtual blocks
c
      call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m2),1,dbl_mb(k_exp),1)
c
c     Print diagonal elements of Fock matrix
c
      write(*,*)'Alpha occupation function Fock elements'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+1,dbl_mb(k_occa+ii),
     &                             dbl_mb(k_exp+ii+nbf*ii)
      enddo
c
c     - symmetrize
c
c     do jj = 0, nbf-1
c       do ii = jj+1, nbf-1
c         tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
c         dbl_mb(k_exp+ii+jj*nbf) = tt
c         dbl_mb(k_exp+jj+ii*nbf) = tt
c       enddo
c     enddo
c
c     - sort stuff
c
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
     &                      dbl_mb(k_occa))
c
c     - establish rotation and apply it to the occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed e: ",info,UERR)
      endif
      call dcopy(nbf*nbf,stateout(wfn1_nxt_aocc()),1,dbl_mb(k_m2),1)
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_m2),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta occupation functions (occupied part)
c     -----------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_exp),1)
      call dcopy(nbf,statein(wfn1_nxt_bconf()),1,dbl_mb(k_occb),1)
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_exp+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_exp+ii),nbf)
      enddo
c
c     Beta occupation functions (virtual part)
c     ----------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_m2),1)
      do ii = 0, nbf-1
        dbl_mb(k_occb+ii) = 1.0d0 - dbl_mb(k_occb+ii)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_m2+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_m2+ii),nbf)
      enddo
      do ii = 0, nbf-1
        dbl_mb(k_occb+ii) = 1.0d0 - dbl_mb(k_occb+ii)
      enddo
c
c     Sum occupied-occupied and virtual-virtual blocks
c
      call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m2),1,dbl_mb(k_exp),1)
c
c     Print diagonal elements of Fock matrix
c
      write(*,*)'Beta occupation function Fock elements'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+11,dbl_mb(k_occb+ii),
     &                              dbl_mb(k_exp+ii+nbf*ii)
      enddo
c
c     - symmetrize
c
c     do jj = 0, nbf-1
c       do ii = jj+1, nbf-1
c         tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
c         dbl_mb(k_exp+ii+jj*nbf) = tt
c         dbl_mb(k_exp+jj+ii*nbf) = tt
c       enddo
c     enddo
c
c     - sort stuff
c
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
      call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
     &                      dbl_mb(k_occb))
c
c     - establish rotation and apply it to both natural orbitals
c       and occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1+nbf),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed g: ",info,UERR)
      endif
      call dcopy(nbf*nbf,stateout(wfn1_nxt_bocc()),1,dbl_mb(k_m2),1)
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_m2),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Print the occupation function eigenvalues
c
c     write(*,*)"occupation function energies"
c     do ir = 1, nbf
c       write(*,'(i4,2f14.6)')ir,dbl_mb(k_m1+(ir-1)),
c    &                           dbl_mb(k_m1+nbf+(ir-1))
c     enddo
c     write(*,*)
c
c     Alpha natural orbitals
c
      call dcopy(nbf    ,statein(wfn1_nxt_aconf()),1,
     &                   stateout(wfn1_nxt_aconf()),1)
c
c     Beta natural orbitals
c
      call dcopy(nbf    ,statein(wfn1_nxt_bconf()),1,
     &                   stateout(wfn1_nxt_bconf()),1)
c
c     Cleanup memory
c
      if (.not.ma_pop_stack(l_gg))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate gg",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occb))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate occb",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occa))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate occa",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_pop_stack(l_exp))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m2))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      end
C>
C> \brief Create a new state canonicalizing the occupation functions
C>
C> As with regular effective 1-electron models there are redundant
C> rotations in this model also. For the occupation functions those
C> redundant rotations are easily identified as mixings in the 
C> occupied-occupied block and the virtual-virtual block. 
C>
C> A problem with canonicalizing the orbitals in this setting is 
C> that because of the degeneracies due to Gilbert's theorem even the
C> canonical orbitals may be arbitrarily mixed. A double
C> canonicalization may be needed to fix this. 
C>
      subroutine wfn1_nxt2_canon_occ(nbf,nea,neb,statein,fock,stateout)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "wfn1_nxt.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision statein(8*nbf*nbf+2*nbf) !< [Input] The initial state
      double precision fock(8*nbf*nbf)    !< [Input] The fock matrices
c
      double precision stateout(8*nbf*nbf+2*nbf) !< [Input] The final state
c
      double precision tt
      double precision levelshift
      parameter(levelshift = -1.0d0)
c
      integer info  !< The status flag for DSYEV
      integer lwork !< The work space size for DSYEV
      integer iocn  !< Offset short hand for occupation numbers
      integer iocf  !< Offset short hand for occupation functions
      integer icnf  !< Offset short hand for electron configurations
c
      integer l_m1,   k_m1
      integer l_m2,   k_m2
      integer l_exp,  k_exp
      integer l_occa, k_occa
      integer l_occb, k_occb
      integer l_gg,   k_gg
c
      integer ir !< Counter
      integer iq !< Counter
      integer ii !< Counter
      integer jj !< Counter
c
c     Allocate memory
c
      lwork = max(nbf*nbf,3*nbf)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp1",l_m1,k_m1))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,lwork,"tmp2",l_m2,k_m2))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"exp",l_exp,k_exp))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf,"occa",l_occa,k_occa))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate occa",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"occb",l_occb,k_occb))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate occb",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf,"gershgorin",l_gg,k_gg))
     &  call errquit("wfn1_nxt2_canon_occ: could not allocate gg",
     &  ma_sizeof(MT_DBL,nbf,MT_BYTE),MA_ERR)
c
c     Initialize state out
c
      call dcopy(wfn1_nxt_size1(),statein,1,stateout,1)
c
c     Alpha occupation functions (occupied part)
c     ------------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_exp),1)
      call dcopy(nbf,statein(wfn1_nxt_aconf()),1,dbl_mb(k_occa),1)
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_exp+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_exp+ii),nbf)
      enddo
c
c     Alpha occupation functions (virtual part)
c     ------------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_aocc()),1,dbl_mb(k_m2),1)
      do ii = 0, nbf-1
        dbl_mb(k_occa+ii) = 1.0d0 - dbl_mb(k_occa+ii)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_m2+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occa+ii),dbl_mb(k_m2+ii),nbf)
      enddo
c
c     Sum occupied-occupied and virtual-virtual blocks
c
      call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m2),1,dbl_mb(k_exp),1)
c
c     Print diagonal elements of Fock matrix
c
      write(*,*)'wfn1_nxt2_canon_occ'
      write(*,*)'Alpha occupation function Fock elements'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+1,1.0d0-dbl_mb(k_occa+ii),
     &                             dbl_mb(k_exp+ii+nbf*ii)
      enddo
c
c     - symmetrize
c
c     do jj = 0, nbf-1
c       do ii = jj+1, nbf-1
c         tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
c         dbl_mb(k_exp+ii+jj*nbf) = tt
c         dbl_mb(k_exp+jj+ii*nbf) = tt
c       enddo
c     enddo
c
c     - sort stuff
c
c     call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
c     call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
c    &                      dbl_mb(k_occa))
c
c     - establish rotation and apply it to the occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed e: ",info,UERR)
      endif
      call dcopy(nbf*nbf,stateout(wfn1_nxt_aocc()),1,dbl_mb(k_m2),1)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_m2),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Print occupation function energies
c
      write(*,*)'Alpha occupation function Hartree-Fock energies'
      do ii = 0, nbf-1
        write(*,'(i4,2f14.6)')ii+1,1.0d0-dbl_mb(k_occa+ii),
     &                             dbl_mb(k_m1+ii)
      enddo
c
c     Beta occupation functions (occupied part)
c     -----------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_exp),1)
      call dcopy(nbf,statein(wfn1_nxt_bconf()),1,dbl_mb(k_occb),1)
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_exp+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_exp+ii),nbf)
      enddo
c
c     Beta occupation functions (virtual part)
c     ----------------------------------------
c
      call dcopy(nbf*nbf,fock(wfn1_nxt_bocc()),1,dbl_mb(k_m2),1)
      do ii = 0, nbf-1
        dbl_mb(k_occb+ii) = 1.0d0 - dbl_mb(k_occb+ii)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_m2+ii*nbf),1)
      enddo
      do ii = 0, nbf-1
        call dscal(nbf,dbl_mb(k_occb+ii),dbl_mb(k_m2+ii),nbf)
      enddo
c
c     Sum occupied-occupied and virtual-virtual blocks
c
      call daxpy(nbf*nbf,1.0d0,dbl_mb(k_m2),1,dbl_mb(k_exp),1)
c
c     - symmetrize
c
c     do jj = 0, nbf-1
c       do ii = jj+1, nbf-1
c         tt = 0.5d0*(dbl_mb(k_exp+ii+jj*nbf)+dbl_mb(k_exp+jj+ii*nbf))
c         dbl_mb(k_exp+ii+jj*nbf) = tt
c         dbl_mb(k_exp+jj+ii*nbf) = tt
c       enddo
c     enddo
c
c     - sort stuff
c
c     call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_gg))
c     call wfn1_nxt2_lshift(nbf,dbl_mb(k_exp),dbl_mb(k_gg),levelshift,
c    &                      dbl_mb(k_occb))
c
c     - establish rotation and apply it to both natural orbitals
c       and occupation functions
c
      info = 0 
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1+nbf),
     &           dbl_mb(k_m2),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_nxt2_canon_occ: dsyev failed g: ",info,UERR)
      endif
      call dcopy(nbf*nbf,stateout(wfn1_nxt_bocc()),1,dbl_mb(k_m2),1)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_m2),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Print occupation function energies
c
      write(*,*)'Alpha-Beta occupation function Hartree-Fock energies'
      do ii = 0, nbf-1
        write(*,'(i4,4f14.6)')ii+1,
     &     1.0d0-dbl_mb(k_occa+ii),
     &     1.0d0-dbl_mb(k_occb+ii),
     &     dbl_mb(k_m1+ii),
     &     dbl_mb(k_m1+nbf+ii)
      enddo
cDEBUG
      write(*,*)'Alpha occupation functions after HF canonicalization'
      call hess_hssout(stateout(wfn1_nxt_aocc()),nbf,nbf,nbf)
      write(*,*)'Beta occupation functions after HF canonicalization'
      call hess_hssout(stateout(wfn1_nxt_bocc()),nbf,nbf,nbf)
cDEBUG
c
c     Print the occupation function eigenvalues
c
c     write(*,*)"occupation function energies"
c     do ir = 1, nbf
c       write(*,'(i4,2f14.6)')ir,dbl_mb(k_m1+(ir-1)),
c    &                           dbl_mb(k_m1+nbf+(ir-1))
c     enddo
c     write(*,*)
c
c     Alpha natural orbitals
c
      call dcopy(nbf    ,statein(wfn1_nxt_aconf()),1,
     &                   stateout(wfn1_nxt_aconf()),1)
c
c     Beta natural orbitals
c
      call dcopy(nbf    ,statein(wfn1_nxt_bconf()),1,
     &                   stateout(wfn1_nxt_bconf()),1)
c
c     Cleanup memory
c
      if (.not.ma_pop_stack(l_gg))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate gg",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occb))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate occb",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_occa))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate occa",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      if (.not.ma_pop_stack(l_exp))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m2))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     &  call errquit("wfn1_nxt2_canon_occ: could not deallocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      end
C> @}
