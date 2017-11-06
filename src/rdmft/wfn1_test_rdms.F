c-----------------------------------------------------------------------
c
!> \ingroup wfn1
!> @{
!>
!> \file wfn1_test_rdms.F
!> Test 1RDM and 2RDM against N-representability conditions
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Test the RDM traces
!>
!> The 1RDMs are supposed to satisfy 
!> \f{eqnarray*}{
!>   n_e^\alpha &=& \tr(D^\alpha) \\\\
!>              &=& \sum_{a,b=1}^{n_b}D^\alpha_{ab}S_{ab} \\\\
!>   n_e^\beta  &=& \tr(D^\beta) \\\\
!>              &=& \sum_{a,b=1}^{n_b}D^\beta_{ab}S_{ab}
!> \f}
!> The 2RDM consists of 4 different blocks: the \f$\alpha\alpha\f$
!> block, the \f$\alpha\beta\f$ block, the \f$\beta\alpha\f$ block,
!> and the \f$\beta\beta\f$ block.
!> For these blocks the following equalities must hold
!> \f{eqnarray*}{
!>   \left(\begin{array}{c}n_e^\alpha \\ 2\end{array}\right)
!>   &=& \tr(\Gamma^{\alpha\alpha}) \\\\
!>   &=& \sum_{a,b,c,d=1}^{n_b}
!>       \Gamma^{\alpha\alpha}_{abcd}S_{ac}S_{bd} \\\\
!>   \left(\begin{array}{c}n_e^\alpha \\ 1\end{array}\right)
!>   \left(\begin{array}{c}n_e^\alpha \\ 1\end{array}\right)
!>   &=& \tr(\Gamma^{\alpha\beta}+\Gamma_{\beta\alpha}) \\\\
!>   &=& \sum_{a,b,c,d=1}^{n_b}\left(
!>         \Gamma^{\alpha\beta}_{abcd}
!>       + \Gamma^{\beta\alpha}_{abcd}\right)
!>         S_{ac}S_{bd} \\\\
!>   \left(\begin{array}{c}n_e^\beta \\ 2\end{array}\right)
!>   &=& \tr(\Gamma^{\beta\beta}) \\\\
!>   &=& \sum_{a,b,c,d=1}^{n_b}
!>       \Gamma^{\beta\beta}_{abcd}S_{ac}S_{bd} 
!> \f}
!> This subroutine evaluates the relevant traces and checks whether
!> the results match the required data.
!>
      subroutine wfn1_test_rdms_trace(wfn1_param,s,da,db,gaa,gab,gbb,
     &                                tol)
      implicit none
#include "wfn1_param.fh"
      type(wfn1_prm),   intent(in) :: wfn1_param
c
c     Overlap matrix
c
      double precision, intent(in) :: s(wfn1_param%nbf,wfn1_param%nbf)
c
c     1RDMs
c
      double precision, intent(in) :: da(wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: db(wfn1_param%nbf,wfn1_param%nbf)
c
c     2RDMs
c
      double precision, intent(in) :: gaa(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: gab(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: gbb(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
c
      double precision, intent(in) :: tol ! tolerance
c
c     Local variables
c
      integer :: nbf ! the number of basis functions
      integer :: nea ! the number of alpha electrons
      integer :: neb ! the number of beta electrons
c
      double precision :: tr_da ! trace of alpha 1RDM
      double precision :: tr_db ! trace of beta 1RDM
c
      double precision :: tr_gaa ! trace of alpha-alpha 2RDM block
      double precision :: tr_gab ! trace of alpha-beta 2RDM block
      double precision :: tr_gbb ! trace of beta-beta 2RDM block
c
      integer :: ia, ib, ic, id  ! counters over atomic basis functions
c
      character*20 pname
      parameter(pname = "wfn1_test_rdms_trace")
c
c     Code
c
      nbf = wfn1_param%nbf
      nea = wfn1_param%nea
      neb = wfn1_param%neb
c
c     Calculate tr(da)
c
      tr_da = 0.0d0
      do ib = 1, nbf
        do ia = 1, nbf
          tr_da = tr_da + da(ia,ib)*s(ia,ib)
        enddo
      enddo
      if (abs(tr_da-1.0d0*nea).gt.tol) then
        write(*,*)pname//" tr(da) failed ",nea,tr_da
      endif
c
c     Calculate tr(db)
c
      tr_db = 0.0d0
      do ib = 1, nbf
        do ia = 1, nbf
          tr_db = tr_db + db(ia,ib)*s(ia,ib)
        enddo
      enddo
      if (abs(tr_db-1.0d0*neb).gt.tol) then
        write(*,*)pname//" tr(db) failed ",neb,tr_db
      endif
c
c     Calculate tr(gaa)
c
      tr_gaa = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ia = 1, nbf
              tr_gaa = tr_gaa + gaa(ia,ib,ic,id)*s(ia,ib)*s(ic,id)
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
      if (abs(tr_gaa-0.5d0*nea*(nea-1)).gt.tol) then
        write(*,*)pname//" tr(gaa) failed ",nea*(nea-1)/2,tr_gaa
      endif
c
c     Calculate tr(gab)
c
      tr_gab = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ia = 1, nbf
              tr_gab = tr_gab + gab(ia,ib,ic,id)*s(ia,ib)*s(ic,id)
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
      if (abs(tr_gab-0.5d0*nea*neb).gt.tol) then
        write(*,*)pname//" tr(gab) failed ",nea,neb,tr_gab
      endif
c
c     Calculate tr(gbb)
c
      tr_gbb = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ia = 1, nbf
              tr_gbb = tr_gbb + gbb(ia,ib,ic,id)*s(ia,ib)*s(ic,id)
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
      if (abs(tr_gbb-0.5d0*neb*(neb-1)).gt.tol) then
        write(*,*)pname//" tr(gbb) failed ",neb*(neb-1)/2,tr_gbb
      endif
c
      end
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Test the RDM permutation symmetries
!>
!> The 1RDM and the 2RDM have to satisfy certain permutation symmetries.
!> Here we consider real valued RDMs. In that case the 1RDM is 
!> symmetric as in
!> \f{eqnarray*}{
!>   D^\sigma_{ab} &=& D^\sigma_{ba}
!> \f}
!> The 2RDM is symmetric under interchanging the labels that belong to
!> electron 1 as well as under interchanging the labels that belong to
!> electron 2. However, interchanging one of the labels belonging to
!> electron 1 with one of the labels belonging to electron 2 changes
!> sign. The way the 2RDM is stored in the program the first two labels
!> associated with electron 1 and the second two labels with electron 2.
!> This routine we first bring the labels in a sensible order of
!> r_1,r_2,r'_1,r'_2. After that the permutation symmetry is checked.
!> This means that the 2RDM must satisfy all conditions:
!> \f{eqnarray*}{
!> \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{abdc} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{bacd} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{badc} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{cdab} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{cdba} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{dcab} \\\\
!> \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{dcba}
!> \f}
!>
      subroutine wfn1_test_rdms_permute(nbf,rdm1,rdm2,tol)
      implicit none
c
      integer, intent(in) :: nbf ! the number of basis functions
c
c     rdm1 is a sigma spin 1RDM
c     rdm2 is a sigma,sigma spin 2RDM (no sigma,sigma' spin)
c
      double precision, intent(in) :: rdm1(nbf,nbf)
      double precision, intent(in) :: rdm2(nbf,nbf,nbf,nbf)
c
      double precision, intent(in) :: tol ! the tolerance
c
c     Local
c
      integer :: ia ! counter over labels a
      integer :: ib ! counter over labels b
      integer :: ic ! counter over labels c
      integer :: id ! counter over labels d
c
      double precision, allocatable :: gaa(:,:,:,:)
c
c     Code
c
      do ia = 1, nbf
        do ib = 1, ia
          if (abs(rdm1(ia,ib)-rdm1(ib,ia)).gt.tol) then
            write(*,*)'failed 1rdm not symmetric'
            write(*,*)ia,ib
            write(*,*)rdm1(ia,ib),rdm1(ib,ia)
          endif
        enddo
      enddo
c
      allocate(gaa(nbf,nbf,nbf,nbf))
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ia = 1, nbf
              gaa(ia,ib,ic,id) = rdm2(ia,ic,ib,id)
            enddo
          enddo
        enddo
      enddo
      do ia = 1, nbf
        do ib = 1, nbf
          do ic = 1, ia
            do id = 1, ib
c
              if(abs(gaa(ia,ib,ic,id)+( 1)*gaa(ia,ib,id,ic)).gt.tol)
     &        then
                write(*,*)'failed -1 abdc'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(ia,ib,id,ic)
              endif
              if(abs(gaa(ia,ib,ic,id)+( 1)*gaa(ib,ia,ic,id)).gt.tol)
     &        then
                write(*,*)'failed -1 bacd'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(ib,ia,ic,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+(-1)*rdm2(ib,ia,id,ic)).gt.tol)
     &        then
                write(*,*)'failed  1 badc'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(ib,ia,id,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+(-1)*rdm2(ic,id,ia,ib)).gt.tol)
     &        then
                write(*,*)'failed  1 cdab'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(ic,id,ia,ib)
              endif
              if(abs(gaa(ia,ib,ic,id)+( 1)*gaa(ic,id,ib,ia)).gt.tol)
     &        then
                write(*,*)'failed -1 cdba'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(ic,id,ib,ia)
              endif
              if(abs(gaa(ia,ib,ic,id)+( 1)*gaa(id,ic,ia,ib)).gt.tol)
     &        then
                write(*,*)'failed  1 dcab'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(id,ic,ia,ib)
              endif
              if(abs(gaa(ia,ib,ic,id)+(-1)*gaa(id,ic,ib,ia)).gt.tol)
     &        then
                write(*,*)'failed  1 dcba'
                write(*,*)ia,ib,ic,id
                write(*,*)gaa(ia,ib,ic,id),gaa(id,ic,ib,ia)
              endif
            enddo ! id
          enddo ! ic
        enddo ! ib
      enddo ! ia
      deallocate(gaa)
      end
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Test the integration of RDMs
!>
!> Integrating the coordinates of one electron out of the 2RDM should
!> produce \f$(n^\sigma_e-1)/2\f$ times the 1RDM (assuming electron
!> 1 and 2 are of same spin). This subroutine tests this requirement.
!>
      subroutine wfn1_test_rdms_integrate(nbf,nea,s,da,gaa,tol)
      implicit none
c
      integer, intent(in) :: nbf ! the number of basis functions
      integer, intent(in) :: nea ! the number of electrons
c
c     Overlap matrix
c
      double precision, intent(in) :: s(nbf,nbf)
c
c     1RDM
c
      double precision, intent(in) :: da(nbf,nbf)
c
c     2RDM
c
      double precision, intent(in) :: gaa(nbf,nbf,nbf,nbf)
c
      double precision, intent(in) :: tol ! tolerance
c
c     Local
c
      double precision, allocatable :: dd(:,:)
c
      integer :: ia, ib, ic, id
c
c     Code
c
      allocate(dd(nbf,nbf))
c
      dd = 0.0d0
      do id = 1, nbf
        do ib = 1, nbf
          do ic = 1, nbf
            do ia = 1, nbf
              dd(ia,ic) = dd(ia,ic)
     +          + gaa(ia,ic,ib,id)*s(ib,id)
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
c
      if (nea.gt.1) then
        do ic = 1, nbf
          do ia = 1, nbf
            dd(ia,ic) = dd(ia,ic)*(1.0d0/(1.0d0*nea-1.0d0))
          enddo ! ia
        enddo ! ic
      endif
c
      do ic = 1, nbf
        do ia = 1, nbf
          if (abs(da(ia,ic)-dd(ia,ic)).gt.tol) then
            write(*,*)'integrate failure ',ia,ic,da(ia,ic),dd(ia,ic)
          endif
        enddo ! ia
      enddo ! ic
c
      deallocate(dd)
c
      end
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Test the 1RDM eigenvalues
!>
!> The eigenvalues of the 1RDM must all be non-negative. Hence
!> express the 1RDM in an orthogonal basis, diagonalize it, and
!> check the eigenvalues.
!>
      subroutine wfn1_test_rdm_eigval1(nbf,s12,da,tol)
      implicit none
c
      integer, intent(in) :: nbf ! the number of basis functions
c
      double precision, intent(in) :: s12(nbf,nbf) ! S^{1/2}
c
c     The 1RDM
c
      double precision, intent(in) :: da(nbf,nbf) 
c
      double precision, intent(in) :: tol ! tolerance
c
c     Local
c
      double precision, allocatable :: da1(:,:)
      double precision, allocatable :: da2(:,:)
      double precision, allocatable :: eig(:)
      double precision, allocatable :: work(:)
c
      double precision :: epsilon = 1.0d-15
c
c     LAPACK variables
c
      integer :: info
      integer :: lwork
c
      integer :: ia, ib, ic, id, ig ! counters
c
c     Code
c
      allocate(da2(nbf,nbf))
      allocate(da1(nbf,nbf))
c
c     Reorder elements
c
c     On input the 2RDM elements are ordered as:
c     - r_1 r'_1; r_2 r'_2
c     we need them in the order
c     - r_1 r_2; r'_1 r'_2
c
      do ib = 1, nbf
        do ia = 1, nbf
          da2(ia,ib) = da(ia,ib)
        enddo ! ia
      enddo ! ic
c
c     Transform the 1RDM into an orthogonal basis
c
c     - Transform index ia
c
      da1 = 0.0d0
      do ib = 1, nbf
        do ia = 1, nbf
          do ig = 1, nbf
            da1(ia,ib) = da1(ia,ib)
     &        + da2(ig,ib)*s12(ig,ia)
          enddo ! ig
        enddo ! ia
      enddo ! ib
c
c     - Transform index ib
c
      da2 = 0.0d0
      do ib = 1, nbf
        do ig = 1, nbf
          do ia = 1, nbf
            da2(ia,ib) = da2(ia,ib)
     &        + da1(ia,ig)*s12(ig,ib)
          enddo ! ig
        enddo ! ia
      enddo ! ib
c
      deallocate(da1)
c
c     Diagonalize the 2RDM
c
      allocate(eig(nbf))
      lwork = -1
      info  = 0
      call dsyev('N','L',nbf,da2,nbf,eig,eig,lwork,info)
      if (info.lt.0) then
        write(*,*)'eigval0: dsyev: invalid argument no. ',-info
      else if (info.gt.0) then
        write(*,*)'eigval0: dsyev: error ',info
      endif
      lwork = eig(1)+1
c
      allocate(work(lwork))
      info  = 0
      call dsyev('N','L',nbf,da2,nbf,eig,work,lwork,info)
      if (info.lt.0) then
        write(*,*)'eigval: dsyev: invalid argument no. ',-info
      else if (info.gt.0) then
        write(*,*)'eigval: dsyev: failed to converge ',info
      endif
c
      do ia = 1, nbf
        write(*,*)'eigval: ',ia,eig(ia)
        if (eig(ia).lt.-epsilon) then
c
c         If the eigenvalue is less than negative the machine
c         precision then it is definitely wrong
c
          write(*,*)'1RDM eigenvalue < 0 ',ia,eig(ia)
        endif
      enddo
c
      deallocate(work,eig,da2)
c
      end
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Test the 2RDM eigenvalues
!>
!> The eigenvalues of the 2RDM must all be non-negative. Hence
!> express the 2RDM in an orthogonal basis, diagonalize it, and
!> check the eigenvalues.
!>
      subroutine wfn1_test_rdm_eigval2(nbf,s12,gaa,tol)
      implicit none
c
      integer, intent(in) :: nbf ! the number of basis functions
c
      double precision, intent(in) :: s12(nbf,nbf) ! S^{1/2}
c
c     The same spin block of the 2RDM
c
      double precision, intent(in) :: gaa(nbf,nbf,nbf,nbf) 
c
      double precision, intent(in) :: tol ! tolerance
c
c     Local
c
      double precision, allocatable :: gaa1(:,:,:,:)
      double precision, allocatable :: gaa2(:,:,:,:)
      double precision, allocatable :: eig(:)
      double precision, allocatable :: work(:)
c
      double precision :: epsilon = 1.0d-15
c
c     LAPACK variables
c
      integer :: info
      integer :: lwork
c
      integer :: ia, ib, ic, id, ig ! counters
      integer :: nbf2 ! dimension of 2RDM as a matrix
c
c     Code
c
      allocate(gaa2(nbf,nbf,nbf,nbf))
      allocate(gaa1(nbf,nbf,nbf,nbf))
c
c     Reorder elements
c
c     On input the 2RDM elements are ordered as:
c     - r_1 r'_1; r_2 r'_2
c     we need them in the order
c     - r_1 r_2; r'_1 r'_2
c
      do id = 1, nbf
        do ib = 1, nbf
          do ic = 1, nbf
            do ia = 1, nbf
              gaa2(ia,ib,ic,id) = gaa(ia,ic,ib,id)
            enddo ! ia
          enddo ! ic
        enddo ! ib
      enddo ! id
c
c     Transform the 2RDM into an orthogonal basis
c
c     - Transform index ia
c
      gaa1 = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ia = 1, nbf
              do ig = 1, nbf
                gaa1(ia,ib,ic,id) = gaa1(ia,ib,ic,id)
     &            + gaa2(ig,ib,ic,id)*s12(ig,ia)
              enddo ! ig
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
c
c     - Transform index ib
c
      gaa2 = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ib = 1, nbf
            do ig = 1, nbf
              do ia = 1, nbf
                gaa2(ia,ib,ic,id) = gaa2(ia,ib,ic,id)
     &            + gaa1(ia,ig,ic,id)*s12(ig,ib)
              enddo ! ig
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
c
c     - Transform index ic
c
      gaa1 = 0.0d0
      do id = 1, nbf
        do ic = 1, nbf
          do ig = 1, nbf
            do ib = 1, nbf
              do ia = 1, nbf
                gaa1(ia,ib,ic,id) = gaa1(ia,ib,ic,id)
     &            + gaa2(ia,ib,ig,id)*s12(ig,ic)
              enddo ! ig
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
c
c     - Transform index id
c
      gaa2 = 0.0d0
      do id = 1, nbf
        do ig = 1, nbf
          do ic = 1, nbf
            do ib = 1, nbf
              do ia = 1, nbf
                gaa2(ia,ib,ic,id) = gaa2(ia,ib,ic,id)
     &            + gaa1(ia,ib,ic,ig)*s12(ig,id)
              enddo ! ig
            enddo ! ia
          enddo ! ib
        enddo ! ic
      enddo ! id
c
      deallocate(gaa1)
c
c     Diagonalize the 2RDM
c
      nbf2 = nbf*nbf
      allocate(eig(nbf2))
      lwork = -1
      info  = 0
      call dsyev('N','L',nbf2,gaa2,nbf2,eig,eig,lwork,info)
      if (info.lt.0) then
        write(*,*)'eigval0: dsyev: invalid argument no. ',-info
      else if (info.gt.0) then
        write(*,*)'eigval0: dsyev: error ',info
      endif
      lwork = eig(1)+1
c
      allocate(work(lwork))
      info  = 0
      call dsyev('N','L',nbf2,gaa2,nbf2,eig,work,lwork,info)
      if (info.lt.0) then
        write(*,*)'eigval: dsyev: invalid argument no. ',-info
      else if (info.gt.0) then
        write(*,*)'eigval: dsyev: failed to converge ',info
      endif
c
      do ia = 1, nbf2
        write(*,*)'eigval: ',ia,eig(ia)
        if (eig(ia).lt.-epsilon) then
c
c         If the eigenvalue is less than negative the machine
c         precision then it is definitely wrong
c
          write(*,*)'2RDM eigenvalue < 0 ',ia,eig(ia)
        endif
      enddo
c
      deallocate(work,eig,gaa2)
c
      end
!>
!> @}
c
c-----------------------------------------------------------------------
c
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Put a set of density matrices to the test
!>
!> RDMs have to satisfy a number of N-representability conditions. 
!> This routine takes 1RDM and 2RMDs and puts them to the various
!> tests. The tests include:
!> - checking the traces,
!> - checking the permutation symmetry
!> - checking that higher order matrices integrate to lower order
!>   ones
!> - checking that the eigenvalues are non-negative
!>
      subroutine wfn1_test_rdms(wfn1_param,wfn1_mat,da,db,gaa,gab,
     &                          gbb,tol)
#include "global.fh"
#include "errquit.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_param.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param
      type(wfn1_matcnst), intent(in) :: wfn1_mat
c
      double precision, intent(in) :: da(wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: db(wfn1_param%nbf,wfn1_param%nbf)
c
      double precision, intent(in) :: gaa(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: gab(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: gbb(wfn1_param%nbf,wfn1_param%nbf,
     &                                    wfn1_param%nbf,wfn1_param%nbf)
      double precision, intent(in) :: tol ! the tolerance
c
c     Local
c
      integer :: nbf
      integer :: nea
      integer :: neb
      integer :: iproc
c
      double precision, allocatable :: s(:,:)
      double precision, allocatable :: s12(:,:)
c
      iproc = ga_nodeid()
c
      nbf = wfn1_param%nbf
      nea = wfn1_param%nea
      neb = wfn1_param%neb
c
      call ga_sync()
      if (iproc.eq.0) then
        allocate(s(nbf,nbf),s12(nbf,nbf))
        call ga_get(wfn1_mat%g_s,1,nbf,1,nbf,s,nbf)
        call ga_get(wfn1_mat%g_s12,1,nbf,1,nbf,s12,nbf)
        write(*,*)'BEGIN RDM tests'
        write(*,*)'Print overlap matrix'
        call wfn1_print_mat(s,nbf,nbf)
        write(*,*)'Print alpha 1RDM'
        call wfn1_print_mat(da,nbf,nbf)
        write(*,*)'Print alpha-alpha 2RDM'
        call wfn1_print_mat2(gaa,nbf)
        write(*,*)'Start: test trace'
        call wfn1_test_rdms_trace(wfn1_param,s,da,db,gaa,gab,gbb,tol)
        write(*,*)'Done : test trace'
        write(*,*)'Start: test permutations alpha-alpha'
        call wfn1_test_rdms_permute(nbf,da,gaa,tol)
        write(*,*)'Do beta-beta'
        call wfn1_test_rdms_permute(nbf,db,gbb,tol)
        write(*,*)'Done : test permutations'
        write(*,*)'Start: test positivity alpha'
        call wfn1_test_rdm_eigval1(nbf,s12,da,tol)
        write(*,*)'Start: test positivity beta'
        call wfn1_test_rdm_eigval1(nbf,s12,db,tol)
        write(*,*)'Start: test positivity alpha-alpha'
        call wfn1_test_rdm_eigval2(nbf,s12,gaa,tol)
        write(*,*)'Do beta-beta'
        call wfn1_test_rdm_eigval2(nbf,s12,gbb,tol)
        write(*,*)'Done : test positivity'
        write(*,*)'Start: test integration alpha-alpha'
        call wfn1_test_rdms_integrate(nbf,nea,s,da,gaa,tol)
        write(*,*)'Do beta-beta'
        call wfn1_test_rdms_integrate(nbf,neb,s,db,gbb,tol)
        write(*,*)'Done : test integration'
        write(*,*)'LEAVE RDM tests'
        deallocate(s,s12)
      endif
      call ga_sync()

      end
!>
!> @}
c
c-----------------------------------------------------------------------
