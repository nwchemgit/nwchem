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
              tr_gaa = tr_gaa + gaa(ia,ib,ic,id)*s(ia,ic)*s(ib,id)
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
              tr_gab = tr_gab + gab(ia,ib,ic,id)*s(ia,ic)*s(ib,id)
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
              tr_gbb = tr_gbb + gbb(ia,ib,ic,id)*s(ia,ic)*s(ib,id)
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
!> sign. This means that the 2RDM must satisfy all conditions:
!> \f{eqnarray*}{
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{abdc} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{acbd} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{adbc} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{acdb} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{adcb} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{bacd} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{badc} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{cabd} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{dabc} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{cadb} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{dacb} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{bcad} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{bdac} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{cbad} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{dbac} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{cdab} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{dcab} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{bcda} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{bdca} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{cbda} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{dbca} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=& -\Gamma^{\sigma\sigma}_{cdba} \\\\
!>  \Gamma^{\sigma\sigma}_{abcd} &=&  \Gamma^{\sigma\sigma}_{dcba} 
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
      do ia = 1, nbf
        do ib = 1, ia
          do ic = 1, ia
            do id = 1, ib
c
c             The following tests were generated by wfn1_2rdm.F
c
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ia,ib,id,ic)).gt.tol)then
                write(*,*)'failed -1 abdc'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ia,ib,id,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ia,ic,ib,id)).gt.tol)then
                write(*,*)'failed -1 acbd'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ia,ic,ib,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ia,id,ib,ic)).gt.tol)then
                write(*,*)'failed  1 adbc'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ia,id,ib,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ia,ic,id,ib)).gt.tol)then
                write(*,*)'failed -1 acdb'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ia,ic,id,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ia,id,ic,ib)).gt.tol)then
                write(*,*)'failed  1 adcb'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ia,id,ic,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ib,ia,ic,id)).gt.tol)then
                write(*,*)'failed -1 bacd'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,ia,ic,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ib,ia,id,ic)).gt.tol)then
                write(*,*)'failed  1 badc'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,ia,id,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ic,ia,ib,id)).gt.tol)then
                write(*,*)'failed  1 cabd'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,ia,ib,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(id,ia,ib,ic)).gt.tol)then
                write(*,*)'failed -1 dabc'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ia,ib,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ic,ia,id,ib)).gt.tol)then
                write(*,*)'failed  1 cadb'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,ia,id,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(id,ia,ic,ib)).gt.tol)then
                write(*,*)'failed -1 dacb'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ia,ic,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ib,ic,ia,id)).gt.tol)then
                write(*,*)'failed -1 bcad'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,ic,ia,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ib,id,ia,ic)).gt.tol)then
                write(*,*)'failed  1 bdac'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,id,ia,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ic,ib,ia,id)).gt.tol)then
                write(*,*)'failed  1 cbad'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,ib,ia,id)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(id,ib,ia,ic)).gt.tol)then
                write(*,*)'failed -1 dbac'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ib,ia,ic)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ic,id,ia,ib)).gt.tol)then
                write(*,*)'failed  1 cdab'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,id,ia,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(id,ic,ia,ib)).gt.tol)then
                write(*,*)'failed -1 dcab'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ic,ia,ib)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(ib,ic,id,ia)).gt.tol)then
                write(*,*)'failed  1 bcda'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,ic,id,ia)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ib,id,ic,ia)).gt.tol)then
                write(*,*)'failed -1 bdca'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ib,id,ic,ia)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ic,ib,id,ia)).gt.tol)then
                write(*,*)'failed -1 cbda'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,ib,id,ia)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(id,ib,ic,ia)).gt.tol)then
                write(*,*)'failed  1 dbca'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ib,ic,ia)
              endif
              if(abs(rdm2(ia,ib,ic,id)+ 1*rdm2(ic,id,ib,ia)).gt.tol)then
                write(*,*)'failed -1 cdba'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(ic,id,ib,ia)
              endif
              if(abs(rdm2(ia,ib,ic,id)+-1*rdm2(id,ic,ib,ia)).gt.tol)then
                write(*,*)'failed  1 dcba'
                write(*,*)ia,ib,ic,id
                write(*,*)rdm2(ia,ib,ic,id),rdm2(id,ic,ib,ia)
              endif
            enddo ! id
          enddo ! ic
        enddo ! ib
      enddo ! ia
      end
!>
!> @}
c
c-----------------------------------------------------------------------
