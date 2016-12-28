c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_2el(wfn1_param,wfn1_inst,
     &                        wfn1_wave,wfn1_deriv,
     &                        v_dca,v_dcb,g_2aa,g_2bb,
     &                        g_waa,g_wab,g_wba,g_wbb,
     &                        g_bfmo,g_momo,E_2el)
      implicit none
C>
C> \brief Calculate the 2-electron terms
C>
C> This subroutine implements the Hartree-Fock 2-electron terms within
C> the wfn1 formulation. An electron correlation model will be added
C> later. As a guiding principle we will not store any quantities that
C> are larger than N^2. Larger quantities will be recomputed as needed.
C>
C> In the algorithm used here we assume that the 2-electron integrals
C> are the most expensive part of the calculation. The cost of the
C> calculation can be reduced by using the Gauchy-Schwarz inequality.
C> The good thing is that we do not need to calculate factors for 
C> integrals we do not calculate. At the moment the Gauchy-Schwarz
C> inequality is not implemented but is planned to be added later.
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
      type(wfn1_prm), intent(in)     :: wfn1_param ! the job parameters
      type(wfn1_prminst), intent(in) :: wfn1_inst  ! the job instances
      type(wfn1_wfn), intent(in)     :: wfn1_wave  ! the wavefunction
      type(wfn1_drv), intent(inout)  :: wfn1_deriv ! the derivative
c
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer, intent(inout) :: g_2aa  ! nmo x nmo matrix for aa-block
                                       ! of the 2-electron density
                                       ! matrix diagonal
      integer, intent(inout) :: g_2bb  ! nmo x nmo matrix for bb-block
                                       ! of the 2-electron density
                                       ! matrix diagonal
      integer, intent(inout) :: g_waa  ! nea x nea matrix for aa-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wab  ! nea x neb matrix for ab-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wba  ! neb x nea matrix for ba-block
                                       ! of the 2-electron interaction
      integer, intent(inout) :: g_wbb  ! neb x nea matrix for ba-block
                                       ! of the 2-electron interaction
c
      integer, intent(inout) :: g_bfmo ! nbf x nmo matrix
      integer, intent(inout) :: g_momo ! nmo x nmo matrix
c
      double precision, intent(out) :: E_2el ! the two-electron energy
c
c     Local
c
      double precision, external :: ddot
c
      double precision, allocatable :: occa(:)  ! a-occupation numbers
      double precision, allocatable :: occb(:)  ! b-occupation numbers
c
      integer :: nproc
      integer :: iproc
      integer :: jproc
c
      integer :: nmo            ! number of MOs
      integer :: nea            ! number of alpha electrons
      integer :: neb            ! number of beta electrons
      integer :: ir, is, it     ! correlation function labels
      integer :: ii, jj, kk, ll ! natural orbital labels
      integer :: ia, ib, ic, id ! atomic orbital labels
      integer :: la, lb, lc, ld ! atomic orbital shell labels
      integer :: na, nb, nc, nd ! number of shell functions
      integer :: nia, nja       ! number of natural orbitals (alpha)
      integer :: nib, njb       ! number of natural orbitals (beta)
c
      integer :: ialo, iahi     ! atomic orbital limits
      integer :: iblo, ibhi     ! atomic orbital limits
      integer :: iclo, ichi     ! atomic orbital limits
      integer :: idlo, idhi     ! atomic orbital limits
c
      integer :: iialo, iiahi   ! natural orbital limits (alpha)
      integer :: ijalo, ijahi   ! natural orbital limits (alpha)
      integer :: iiblo, iibhi   ! natural orbital limits (beta)
      integer :: ijblo, ijbhi   ! natural orbital limits (beta)
c
      integer :: max2e          ! quartet maximum no. 2-el integrals
      integer :: mem2           ! amount of scratch space needed
      integer :: maxsbf         ! maximum no. of shell basis functions
      integer :: nshell         ! the number of shells
c
      integer, parameter :: msg_null = 0
c
      double precision, pointer :: buf(:)
      double precision, pointer :: eri(:,:,:,:)
      double precision, pointer :: buf2(:)
      double precision, pointer :: den(:,:,:,:)
      double precision, pointer :: buf3a(:)
      double precision, pointer :: d3a(:,:,:)
      double precision, pointer :: buf3b(:)
      double precision, pointer :: d3b(:,:,:)
      double precision, pointer :: bufaa(:)
      double precision, pointer :: naia(:,:)
      double precision, pointer :: bufba(:)
      double precision, pointer :: nbia(:,:)
      double precision, pointer :: bufca(:)
      double precision, pointer :: ncia(:,:)
      double precision, pointer :: bufda(:)
      double precision, pointer :: ndia(:,:)
      double precision, pointer :: bufab(:)
      double precision, pointer :: naib(:,:)
      double precision, pointer :: bufbb(:)
      double precision, pointer :: nbib(:,:)
      double precision, pointer :: bufcb(:)
      double precision, pointer :: ncib(:,:)
      double precision, pointer :: bufdb(:)
      double precision, pointer :: ndib(:,:)
      double precision, pointer :: bufija(:)
      double precision, pointer :: occija(:,:)
      double precision, pointer :: bufijb(:)
      double precision, pointer :: occijb(:,:)
      double precision, pointer :: bufaca(:)
      double precision, pointer :: daca(:,:) ! a,c 1-RDM alpha
      double precision, pointer :: bufacb(:)
      double precision, pointer :: dacb(:,:) ! a,c 1-RDM beta
      double precision, pointer :: bufbda(:)
      double precision, pointer :: dbda(:,:) ! b,d 1-RDM alpha
      double precision, pointer :: bufbdb(:)
      double precision, pointer :: dbdb(:,:) ! b,d 1-RDM beta
      double precision, allocatable :: scr(:)
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      E_2el = 0.0d0
      nmo = wfn1_param%nmo
      nea = wfn1_param%nea
      neb = wfn1_param%neb
      if (.not. bas_numcont(wfn1_inst%basis, nshell)) call errquit
     &     ('wfn1_e_g_2el: bas_numcont failed for basis', 
     &      wfn1_inst%basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(wfn1_inst%basis, maxsbf)) call errquit
     &     ('wfn1_e_g_2el: bas_nbf_cn_max failed for basis', 
     &      wfn1_inst%basis, BASIS_ERR)
      call ga_distribution(g_2aa,0,iialo,iiahi,ijalo,ijahi)
      call ga_distribution(g_2bb,0,iiblo,iibhi,ijblo,ijbhi)
      nia = iiahi-iialo+1
      nja = ijahi-ijalo+1
      nib = iibhi-iiblo+1
      njb = ijbhi-ijblo+1
      call int_mem_2e4c(max2e,mem2)
      allocate(occa(1:nmo),occb(1:nmo))
      allocate(buf(max2e),buf2(max2e),scr(mem2))
      allocate(bufaa(maxsbf*nmo),bufba(maxsbf*nmo))
      allocate(bufca(maxsbf*nmo),bufda(maxsbf*nmo))
      allocate(bufab(maxsbf*nmo),bufbb(maxsbf*nmo))
      allocate(bufcb(maxsbf*nmo),bufdb(maxsbf*nmo))
      allocate(bufija(nia*nja),bufijb(nib*njb))
      allocate(bufaca(maxsbf*maxsbf),bufacb(maxsbf*maxsbf))
      allocate(bufbda(maxsbf*maxsbf),bufbdb(maxsbf*maxsbf))
      allocate(buf3a(maxsbf*maxsbf*max(nia,nja)))
      allocate(buf3b(maxsbf*maxsbf*max(nib,njb)))
      occa = 0.0d0
      occb = 0.0d0
c
c     Calculate occupation numbers
c     ----------------------------
c
c     Calculate alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,occa)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,occb)
c
c     Calculate alpha-alpha-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                    g_2aa)
c
c     Calculate beta-beta-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_cb,nmo,nea,wfn1_wave%icnfb,
     &                    g_2bb)
c
c     Shell loops
c     -----------
c
      do ld = 1, nshell
        if (.not. bas_cn2bfr(wfn1_inst%basis, ld, idlo, idhi))
     &      call errquit('wfn1_e_g_2el: bas_cn2bfr', 
     &                   wfn1_inst%basis, BASIS_ERR)
        nd = idhi-idlo+1
        ndia(idlo:idhi,1:nmo)=>bufda
        ndib(idlo:idhi,1:nmo)=>bufdb
        call ga_get(wfn1_wave%g_na,idlo,idhi,1,nmo,ndia,nd)
        call ga_get(wfn1_wave%g_nb,idlo,idhi,1,nmo,ndib,nd)
        do lc = 1, nshell
          if (.not. bas_cn2bfr(wfn1_inst%basis, lc, iclo, ichi))
     &        call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                     wfn1_inst%basis, BASIS_ERR)
          nc = ichi-iclo+1
          ncia(iclo:ichi,1:nmo)=>bufca
          ncib(iclo:ichi,1:nmo)=>bufcb
          call ga_get(wfn1_wave%g_na,iclo,ichi,1,nmo,ncia,nc)
          call ga_get(wfn1_wave%g_nb,iclo,ichi,1,nmo,ncib,nc)
          do lb = 1, nshell
            if (.not. bas_cn2bfr(wfn1_inst%basis, lb, iblo, ibhi))
     &          call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                       wfn1_inst%basis, BASIS_ERR)
            nb = ibhi-iblo+1
            nbia(iblo:ibhi,1:nmo)=>bufba
            nbib(iblo:ibhi,1:nmo)=>bufbb
            call ga_get(wfn1_wave%g_na,iblo,ibhi,1,nmo,nbia,nb)
            call ga_get(wfn1_wave%g_nb,iblo,ibhi,1,nmo,nbib,nb)
            do la = 1, nshell
              if (.not. bas_cn2bfr(wfn1_inst%basis, la, ialo, iahi))
     &            call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                         wfn1_inst%basis, BASIS_ERR)
              na = iahi-ialo+1
              naia(ialo:iahi,1:nmo)=>bufaa
              naib(ialo:iahi,1:nmo)=>bufab
              call ga_get(wfn1_wave%g_na,ialo,iahi,1,nmo,naia,na)
              call ga_get(wfn1_wave%g_nb,ialo,iahi,1,nmo,naib,na)
              eri(ialo:iahi,iclo:ichi,iblo:ibhi,idlo:idhi)=>buf
              den(ialo:iahi,iclo:ichi,iblo:ibhi,idlo:idhi)=>buf2
              call int_2e4c(wfn1_inst%basis, ld, lb,
     &                      wfn1_inst%basis, lc, la,
     &                      mem2, scr, max2e, eri)
              den = 0.0d0
c
c             calculate the 1RDM terms
c
              daca(ialo:iahi,iclo:ichi)=>bufaca
              dacb(ialo:iahi,iclo:ichi)=>bufacb
              dbda(iblo:ibhi,idlo:idhi)=>bufbda
              dbdb(iblo:ibhi,idlo:idhi)=>bufbdb
              dbda = 0.0d0
              dbdb = 0.0d0
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ii = 1, nmo
                    dbda(ib,id) = dbda(ib,id)
     &                          + nbia(ib,ii)*ndia(id,ii)*occa(ii)
                    dbdb(ib,id) = dbdb(ib,id)
     &                          + nbib(ib,ii)*ndib(id,ii)*occb(ii)
                  enddo
                enddo
              enddo
              daca = 0.0d0
              dacb = 0.0d0
              do ic = iclo, ichi
                do ia = ialo, iahi
                  do ii = 1, nmo
                    daca(ia,ic) = daca(ia,ic)
     &                          + naia(ia,ii)*ncia(ic,ii)*occa(ii)
                    dacb(ia,ic) = dacb(ia,ic)
     &                          + naib(ia,ii)*ncib(ic,ii)*occb(ii)
                  enddo
                enddo
              enddo
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ic = iclo, ichi
                    do ia = ialo, iahi
                      den(ia,ic,ib,id) = den(ia,ic,ib,id) 
     &                                 + daca(ia,ic)*dbdb(ib,id)
     &                                 + dacb(ia,ic)*dbda(ib,id)
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
c
c             calculate the 2RDM terms
c
              do jproc = 0, nproc-1
                call ga_distribution(g_2aa,jproc,iialo,iiahi,
     &                                           ijalo,ijahi)
                call ga_distribution(g_2bb,jproc,iiblo,iibhi,
     &                                           ijblo,ijbhi)
                nia = iiahi-iialo+1
                nib = iibhi-iiblo+1
                occija(iialo:iiahi,ijalo:ijahi)=>bufija
                occijb(iiblo:iibhi,ijblo:ijbhi)=>bufijb
                call ga_get(g_2aa,iialo,iiahi,ijalo,ijahi,occija,nia)
                call ga_get(g_2bb,iiblo,iibhi,ijblo,ijbhi,occijb,nib)
                d3a(ialo:iahi,iclo:ichi,ijalo:ijahi)=>buf3a
                d3b(ialo:iahi,iclo:ichi,ijblo:ijbhi)=>buf3b
                d3a = 0.0d0
                d3b = 0.0d0
                do ic = iclo, ichi
                  do ia = ialo, iahi
                    do jj = ijalo, ijahi
                      do ii = iialo, iiahi
                        d3a(ia,ic,jj) = d3a(ia,ic,jj) 
     &                                + naia(ia,ii)*ncia(ic,ii)*
     &                                  occija(ii,jj)
                      enddo
                    enddo
                    do jj = ijblo, ijbhi
                      do ii = iiblo, iibhi
                        d3b(ia,ic,jj) = d3b(ia,ic,jj) 
     &                                + naib(ia,ii)*ncib(ic,ii)*
     &                                  occijb(ii,jj)
                      enddo
                    enddo
                  enddo
                enddo
                do id = idlo, idhi
                  do ib = iblo, ibhi
                    do ic = iclo, ichi
                      do ia = ialo, iahi
                        do jj = ijalo, ijahi
                          den(ia,ic,ib,id) = den(ia,ic,ib,id) 
     &                                     + nbia(ib,jj)*ndia(id,jj)*
     &                                       d3a(ia,ic,jj)
                        enddo ! jj
                        do jj = ijblo, ijbhi
                          den(ia,ic,ib,id) = den(ia,ic,ib,id) 
     &                                     + nbib(ib,jj)*ndib(id,jj)*
     &                                       d3b(ia,ic,jj)
                        enddo ! jj
                      enddo ! ia
                    enddo ! ic
                  enddo ! ib
                enddo ! id
              enddo ! jproc
c
c             Add the 2-electron energy
c
              E_2el = E_2el + ddot(na*nb*nc*nd,eri,1,den,1)
c
            enddo ! la
          enddo ! lb
        enddo ! lc
      enddo ! ld
c
c     Clean up memory
c     ---------------
c
      deallocate(buf3a,  buf3b)
      deallocate(bufbda, bufbdb)
      deallocate(bufaca, bufacb)
      deallocate(bufija, bufijb)
      deallocate(bufcb,  bufdb)
      deallocate(bufab,  bufbb)
      deallocate(bufca,  bufda)
      deallocate(bufaa,  bufba)
      deallocate(buf,    buf2, scr)
      deallocate(occa,   occb)
c
c     Global sums
c     -----------
c
      call ga_dgop(msg_null,E_2el,1,'+')
c
      end subroutine wfn1_e_g_2el
c
c-----------------------------------------------------------------------
