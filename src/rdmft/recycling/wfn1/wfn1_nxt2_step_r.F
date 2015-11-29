C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Create a new state by moving one step from the old one
C>
C> The new state is generated from the old one by applying appropriate
C> rotations.
C>
      subroutine wfn1_nxt2_step_r(nea,neb,nbf,it,statein,fock,res,
     &                            stateout)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "wfn1_nxt.fh"
c
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
      integer nbf !< [Input] The number of basis functions
      integer it  !< [Input] The iteration number
c
      double precision statein(8*nbf*nbf) !< [Input] The initial state
      double precision fock(8*nbf*nbf)    !< [Input] The fock matrices
      double precision res(8*nbf*nbf)     !< [Input] The residuals
c
      double precision stateout(8*nbf*nbf) !< [Input] The final state
c
      integer info  !< The status flag for DSYEV
      integer lwork !< The work space size for DSYEV
c
      integer l_m1,  k_m1
      integer l_m2,  k_m2
      integer l_exp, k_exp
      integer ir
c
c     Allocate memory
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp1",l_m1,k_m1))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"tmp2",l_m2,k_m2))
     &  call errquit("wfn1_nxt2_step: could not allocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"exp",l_exp,k_exp))
     &  call errquit("wfn1_nxt2_step: could not allocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
c     Alpha occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_aocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,res(wfn1_nxt_aocc()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
cDEBUG
      do ir = 1, nbf
        write(*,*)"alpha occ eig in:",ir,dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
      enddo
      write(*,*)
cDEBUG
      info = 0 
      lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
cDEBUG
      do ir = 1, nbf
        write(*,*)"alpha occ eig ou:",ir,dbl_mb(k_m1+(ir-1))
      enddo
      write(*,*)
      write(*,*)
cDEBUG
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step_r: dsyev_nsrt failed a: ",
     &               info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: alpha occupation functions"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_aocc()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta occupation functions
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_bocc()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,res(wfn1_nxt_bocc()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
cDEBUG
      do ir = 1, nbf
        write(*,*)"beta occ eig in:",ir,dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
      enddo
      write(*,*)
cDEBUG
      info = 0 
      lwork = nbf*nbf
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
cDEBUG
      do ir = 1, nbf
        write(*,*)"beta occ eig ou:",ir,dbl_mb(k_m1+(ir-1))
      enddo
      write(*,*)
      write(*,*)
cDEBUG
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step_r: dsyev_nsrt failed b: ",
     &               info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: beta occupation functions"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bocc()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bocc()),nbf,
c    &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Alpha natural orbitals
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_anat()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,res(wfn1_nxt_anat()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      info = 0 
      lwork = nbf*nbf
cDEBUG
      do ir = 1, nbf
        write(*,*)"alpha nat eig in:",ir,dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
      enddo
      write(*,*)
cDEBUG
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
cDEBUG
      do ir = 1, nbf
        write(*,*)"alpha nat eig ou:",ir,dbl_mb(k_m1+(ir-1))
      enddo
      write(*,*)
      write(*,*)
cDEBUG
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step_r: dsyev_nsrt failed a: ",
     &               info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: alpha natural orbitals"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_anat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_anat()),nbf)
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_aocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_aocc()),nbf)
c
c     Beta natural orbitals
c
c     call wfn1_exp(nbf,grad(wfn1_nxt_bnat()),dbl_mb(k_m1),dbl_mb(k_m2),
c    &              dbl_mb(k_exp))
      call dcopy(nbf*nbf,res(wfn1_nxt_bnat()),1,dbl_mb(k_exp),1)
      call wfn1_nxt2_gershgorin(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call wfn1_nxt2_shift(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      info = 0 
      lwork = nbf*nbf
cDEBUG
      do ir = 1, nbf
        write(*,*)"beta nat eig in:",ir,dbl_mb(k_exp+(ir-1)+(ir-1)*nbf)
      enddo
      write(*,*)
cDEBUG
      call dsyev('v','u',nbf,dbl_mb(k_exp),nbf,dbl_mb(k_m1),
     &           dbl_mb(k_m2),lwork,info)
cDEBUG
      do ir = 1, nbf
        write(*,*)"beta nat eig ou:",ir,dbl_mb(k_m1+(ir-1))
      enddo
      write(*,*)
      write(*,*)
cDEBUG
      if (info.ne.0) then
        call errquit("wfn1_nxt2_step_r: dsyev_nsrt failed b: ",
     &               info,UERR)
      endif
cDEBUG
c     write(*,*)"wfn1_nxt2_step: beta natural orbitals"
cDEBUG
c     call wfn1_nxt2_sort_columns(nbf,dbl_mb(k_exp),dbl_mb(k_m1))
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,statein(wfn1_nxt_bnat()),nbf,
     &           dbl_mb(k_exp),nbf,0.0d0,stateout(wfn1_nxt_bnat()),nbf)
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,dbl_mb(k_exp),nbf,
c    &           statein(wfn1_nxt_bocc()),nbf,
c    &           0.0d0,stateout(wfn1_nxt_bocc()),nbf)
c
c     Cleanup memory
c
      if (.not.ma_pop_stack(l_exp))
     &  call errquit("wfn1_nxt2_step: could not deallocate exp",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m2))
     &  call errquit("wfn1_nxt2_step: could not deallocate tmp2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     &  call errquit("wfn1_nxt2_step: could not deallocate tmp1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      end
C>
C> \brief Shift the matrix diagonal elements to ensure non-overlapping
C> Gershgorin discs
C>
C> LAPACK will sort the eigenvalues (and eigenvectors) of a matrix in
C> increasing order. In our application the order of the eigenvectors
C> is physically significant and hence we cannot tolerate this 
C> behavior. Therefore we need to modify the diagonal elements such 
C> the eigenvalues of the matrix get sorted into the same order as
C> diagonal elements of the input matrix. I.e. we need to change the
C> diagonal elements such that the Gershgorin discs of the resulting
C> matrix do not overlap.
C>
      subroutine wfn1_nxt2_shift(nbf,a,r)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision a(nbf,nbf) !< [In/Output] The matrix in question
      double precision r(nbf)     !< [Input] The Gershgorin radii
c
      integer i !< Counter
c
      do i = 2, nbf
        a(i,i) = max(a(i-1,i-1)+r(i-1)+r(i)+0.0d0,a(i,i))
      enddo
c
      end
C>
C> \brief Shift the matrix diagonal using the Gershgorin discs and
C> an additional level shifter
C>
      subroutine wfn1_nxt2_lshift(nbf,a,r,s,o)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision a(nbf,nbf) !< [In/Output] The matrix in question
      double precision r(nbf)     !< [Input] The Gershgorin radii
      double precision s          !< [Input] The level shifter
      double precision o(nbf)     !< [Input] The occupation numbers
c
      integer i !< Counter
c
      if (s.ge.0.0d0) then
        do i = 2, nbf
c         a(i,i) = max(a(i-1,i-1)+r(i-1)+r(i)+s*(1.0d0-o(i)),a(i,i))
          a(i,i) = max(a(i-1,i-1)+r(i-1)+r(i),a(i,i))
        enddo
        do i = 1, nbf
          a(i,i) = a(i,i)+s*(1.0d0-o(i))
        enddo
      else
        do i = 2, nbf
          a(i,i) = max(a(i-1,i-1)+r(i-1)+r(i),a(i,i))
        enddo
      endif
c
      end
C>
C> @}
