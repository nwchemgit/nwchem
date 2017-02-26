c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_2el(wfn1_param,wfn1_inst,
     &                        wfn1_wave,wfn1_deriv,
     &                        v_dca,v_dcb,g_2aa,g_2bb,E_2el)
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
      double precision, allocatable :: ocxa(:)  ! a-average exchange 
                                                ! occupation numbers
                                                ! of 2RDM
      double precision, allocatable :: ocxb(:)  ! b-average exchange
                                                ! occupation numbers
                                                ! of 2RDM
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
      double precision, pointer :: denx(:,:,:,:)
c
      double precision, pointer :: buf3aad(:)
      double precision, pointer :: d3aad(:,:,:)
      double precision, pointer :: buf3abd(:)
      double precision, pointer :: d3abd(:,:,:)
      double precision, pointer :: buf3bad(:)
      double precision, pointer :: d3bad(:,:,:)
      double precision, pointer :: buf3bbd(:)
      double precision, pointer :: d3bbd(:,:,:)
c
      double precision, pointer :: buf3aac(:)
      double precision, pointer :: d3aac(:,:,:)
      double precision, pointer :: buf3abc(:)
      double precision, pointer :: d3abc(:,:,:)
      double precision, pointer :: buf3bac(:)
      double precision, pointer :: d3bac(:,:,:)
      double precision, pointer :: buf3bbc(:)
      double precision, pointer :: d3bbc(:,:,:)
c
      double precision, pointer :: bufaa(:)
      double precision, pointer :: naia(:,:)
      double precision, pointer :: bufba(:)
      double precision, pointer :: nbia(:,:)
      double precision, pointer :: bufca(:)
      double precision, pointer :: ncia(:,:)
      double precision, pointer :: bufda(:)
      double precision, pointer :: ndia(:,:)
c
      double precision, pointer :: bufab(:)
      double precision, pointer :: naib(:,:)
      double precision, pointer :: bufbb(:)
      double precision, pointer :: nbib(:,:)
      double precision, pointer :: bufcb(:)
      double precision, pointer :: ncib(:,:)
      double precision, pointer :: bufdb(:)
      double precision, pointer :: ndib(:,:)
c
      double precision, pointer :: bufija(:)
      double precision, pointer :: occija(:,:)
      double precision, pointer :: bufijb(:)
      double precision, pointer :: occijb(:,:)
c
      double precision, pointer :: bufaca(:)
      double precision, pointer :: daca(:,:) ! a,c 1-RDM alpha
      double precision, pointer :: bufacb(:)
      double precision, pointer :: dacb(:,:) ! a,c 1-RDM beta
      double precision, pointer :: bufbda(:)
      double precision, pointer :: dbda(:,:) ! b,d 1-RDM alpha
      double precision, pointer :: bufbdb(:)
      double precision, pointer :: dbdb(:,:) ! b,d 1-RDM beta
      double precision, pointer :: bufxbda(:)
      double precision, pointer :: xbda(:,:) ! b,d 1-RDM alpha
      double precision, pointer :: bufxbdb(:)
      double precision, pointer :: xbdb(:,:) ! b,d 1-RDM beta
c
      double precision, pointer :: buffnac(:)
      double precision, pointer :: fnac(:,:)
      double precision, pointer :: fn1aca(:,:) ! a,c 1-electron natural orbital Fock alpha
      double precision, pointer :: fn1ica(:,:) ! i,c 1-electron natural orbital Fock alpha
      double precision, pointer :: fn2ica(:,:) ! i,c 2-electron natural orbital Fock alpha
      double precision, pointer :: fn2ida(:,:) ! i,d 2-electron natural orbital Fock alpha
      double precision, pointer :: buffnad(:)
      double precision, pointer :: fc1aca(:,:) ! a,c 1-electron correlation function Fock alpha
      double precision, pointer :: buffnbc(:)
      double precision, pointer :: fnbc(:,:)
      double precision, pointer :: fn1acb(:,:) ! a,c 1-electron natural orbital Fock beta
      double precision, pointer :: fn1icb(:,:) ! i,c 1-electron natural orbital Fock beta
      double precision, pointer :: fn2icb(:,:) ! i,c 2-electron natural orbital Fock beta
      double precision, pointer :: fn2idb(:,:) ! i,d 2-electron natural orbital Fock beta
      double precision, pointer :: buffnbd(:)
      double precision, pointer :: fc1acb(:,:) ! a,c 1-electron correlation function Fock beta
      double precision, pointer :: buffna(:)
      double precision, pointer :: fna(:,:)
      double precision, pointer :: buffnb(:)
      double precision, pointer :: fnb(:,:)
c
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
      allocate(ocxa(1:nmo),ocxb(1:nmo))
      allocate(buf(max2e),buf2(max2e),scr(mem2))
      allocate(bufaa(maxsbf*nmo),bufba(maxsbf*nmo))
      allocate(bufca(maxsbf*nmo),bufda(maxsbf*nmo))
      allocate(bufab(maxsbf*nmo),bufbb(maxsbf*nmo))
      allocate(bufcb(maxsbf*nmo),bufdb(maxsbf*nmo))
      allocate(bufija(nia*nja),bufijb(nib*njb))
      allocate(bufaca(maxsbf*maxsbf),bufacb(maxsbf*maxsbf))
      allocate(bufbda(maxsbf*maxsbf),bufbdb(maxsbf*maxsbf))
      allocate(bufxbda(maxsbf*maxsbf),bufxbdb(maxsbf*maxsbf))
c
      allocate(buf3aac(maxsbf*maxsbf*max(nia,nja)))
      allocate(buf3aad(maxsbf*maxsbf*max(nia,nja)))
      allocate(buf3abc(maxsbf*maxsbf*max(nia,nja)))
      allocate(buf3abd(maxsbf*maxsbf*max(nia,nja)))
      allocate(buf3bac(maxsbf*maxsbf*max(nib,njb)))
      allocate(buf3bad(maxsbf*maxsbf*max(nib,njb)))
      allocate(buf3bbc(maxsbf*maxsbf*max(nib,njb)))
      allocate(buf3bbd(maxsbf*maxsbf*max(nib,njb)))
c
      allocate(buffnac(max(nia,nja,maxsbf)*nmo))
      allocate(buffnad(max(nia,nja,maxsbf)*nmo))
      allocate(buffnbc(max(nib,njb,maxsbf)*nmo))
      allocate(buffnbd(max(nib,njb,maxsbf)*nmo))
      allocate(buffna(max(nia,nja)*nmo))
      allocate(buffnb(max(nib,njb)*nmo))
c
c     Calculate occupation numbers
c     ----------------------------
c
c     Calculate alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,occa)
      call wfn1_calc_ocx(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,ocxa)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,occb)
      call wfn1_calc_ocx(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,ocxb)
c
c     Calculate alpha-alpha-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                    g_2aa)
c
c     Calculate beta-beta-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
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
c
c             calculate the 2-electron integrals
c
              call int_2e4c(wfn1_inst%basis, ld, lb,
     &                      wfn1_inst%basis, lc, la,
     &                      mem2, scr, max2e, eri)
              den = 0.0d0
c
c             calculate the 1RDM terms
c
              fn1aca(ialo:iahi,iclo:ichi)=>buffnad
              fn1acb(ialo:iahi,iclo:ichi)=>buffnbd
              fn1aca = 0.0d0
              fn1acb = 0.0d0
              daca(ialo:iahi,iclo:ichi)=>bufaca
              dacb(ialo:iahi,iclo:ichi)=>bufacb
              dbda(iblo:ibhi,idlo:idhi)=>bufbda
              dbdb(iblo:ibhi,idlo:idhi)=>bufbdb
              xbda(iblo:ibhi,idlo:idhi)=>bufxbda
              xbdb(iblo:ibhi,idlo:idhi)=>bufxbdb
              daca = 0.0d0
              dacb = 0.0d0
              dbda = 0.0d0
              dbdb = 0.0d0
              xbda = 0.0d0
              xbdb = 0.0d0
c
c             Construct the 1RDM and 2RDM density matrix terms
c
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ii = 1, nmo
                    dbda(ib,id) = dbda(ib,id)
     &                          + nbia(ib,ii)*ndia(id,ii)*occa(ii)
                    dbdb(ib,id) = dbdb(ib,id)
     &                          + nbib(ib,ii)*ndib(id,ii)*occb(ii)
                    xbda(ib,id) = xbda(ib,id)
     &                          + nbia(ib,ii)*ndia(id,ii)*ocxa(ii)
                    xbdb(ib,id) = xbdb(ib,id)
     &                          + nbib(ib,ii)*ndib(id,ii)*ocxb(ii)
                  enddo
                enddo
              enddo
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
c             Derivatives wrt the natural orbitals
c
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ic = iclo, ichi
                    do ia = ialo, iahi
                      fn1aca(ia,ic) = fn1aca(ia,ic)
     &                              + eri(ia,ic,ib,id)*dbdb(ib,id)
                      fn1acb(ia,ic) = fn1acb(ia,ic)
     &                              + eri(ia,ic,ib,id)*dbda(ib,id)
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
              fn1ica(1:nmo,iclo:ichi)=>buffnac
              fn1icb(1:nmo,iclo:ichi)=>buffnbc
              fn1ica = 0.0d0
              fn1icb = 0.0d0
              do ic = iclo, ichi
                do ia = ialo, iahi
                  do ii = 1, nmo
                    fn1ica(ii,ic) = fn1ica(ii,ic)
     &                            + naia(ia,ii)*occa(ii)*fn1aca(ia,ic)
                    fn1icb(ii,ic) = fn1icb(ii,ic)
     &                            + naib(ia,ii)*occb(ii)*fn1acb(ia,ic)
                  enddo
                enddo
              enddo
              do jj = 1, nmo, nja
                fna(1:nmo,jj:min(nmo,jj+nja-1))=>buffna
                fna = 0.0d0
                do kk = jj, min(nmo,jj+nja-1)
                  do ic = iclo, ichi
                    do ii = 1, nmo
                      fna(ii,kk) = fna(ii,kk) 
     &                           + 1.0d0*fn1ica(ii,ic)*ncia(ic,kk)
                    enddo ! ii
                  enddo ! ic
                enddo ! kk
                call ga_acc(wfn1_deriv%g_dna,1,nmo,jj,min(nmo,jj+nja-1),
     &                      fna,nmo,1.0d0)
              enddo ! jj
              do jj = 1, nmo, njb
                fnb(1:nmo,jj:min(nmo,jj+njb-1))=>buffnb
                fnb = 0.0d0
                do kk = jj, min(nmo,jj+njb-1)
                  do ic = iclo, ichi
                    do ii = 1, nmo
                      fnb(ii,kk) = fnb(ii,kk) 
     &                           + 1.0d0*fn1icb(ii,ic)*ncib(ic,kk)
                    enddo ! ii
                  enddo ! ic
                enddo ! kk
                call ga_acc(wfn1_deriv%g_dnb,1,nmo,jj,min(nmo,jj+njb-1),
     &                      fnb,nmo,1.0d0)
              enddo
c
c             Derivatives wrt the correlation functions
c
              fc1aca(ialo:iahi,iclo:ichi)=>buffnad
              fc1acb(ialo:iahi,iclo:ichi)=>buffnbd
              fc1aca = 0.0d0
              fc1acb = 0.0d0
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ic = iclo, ichi
                    do ia = ialo, iahi
                      fc1aca(ia,ic) = fc1aca(ia,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.5d0*(dbda(ib,id)+xbda(ib,id))
     &                                 +dbdb(ib,id))
                      fc1acb(ia,ic) = fc1acb(ia,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.5d0*(dbdb(ib,id)+xbdb(ib,id))
     &                                 +dbda(ib,id))
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
              do ii = 1, nmo
                do ic = iclo, ichi
                  do ia = ialo, iahi
                    v_dca(ii) = v_dca(ii)
     &                        + 1.0d0*naia(ia,ii)*
     &                          fc1aca(ia,ic)*ncia(ic,ii)
                    v_dcb(ii) = v_dcb(ii)
     &                        + 1.0d0*naib(ia,ii)*
     &                          fc1acb(ia,ic)*ncib(ic,ii)
                  enddo
                enddo
              enddo
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
                d3aac(ialo:iahi,iclo:ichi,iialo:iiahi)=>buf3aac
                d3aad(ialo:iahi,idlo:idhi,iialo:iiahi)=>buf3aad
                d3abc(iblo:ibhi,iclo:ichi,iialo:iiahi)=>buf3abc
                d3abd(iblo:ibhi,idlo:idhi,iialo:iiahi)=>buf3abd
                d3bac(ialo:iahi,iclo:ichi,iiblo:iibhi)=>buf3bac
                d3bad(ialo:iahi,idlo:idhi,iiblo:iibhi)=>buf3bad
                d3bbc(iblo:ibhi,iclo:ichi,iiblo:iibhi)=>buf3bbc
                d3bbd(iblo:ibhi,idlo:idhi,iiblo:iibhi)=>buf3bbd
                d3aac = 0.0d0
                d3aad = 0.0d0
                d3abc = 0.0d0
                d3abd = 0.0d0
                d3bac = 0.0d0
                d3bad = 0.0d0
                d3bbc = 0.0d0
                d3bbd = 0.0d0
                fn2ica(iialo:iiahi,iclo:ichi)=>buffnac
                fn2ida(iialo:iiahi,idlo:idhi)=>buffnad
                fn2icb(iiblo:iibhi,iclo:ichi)=>buffnbc
                fn2idb(iiblo:iibhi,idlo:idhi)=>buffnbd
                fn2ica = 0.0d0
                fn2ida = 0.0d0
                fn2icb = 0.0d0
                fn2idb = 0.0d0
                fna(iialo:iiahi,1:nmo)=>buffna
                fnb(iiblo:iibhi,1:nmo)=>buffnb
                fna = 0.0d0
                fnb = 0.0d0
                do id = idlo, idhi
                  do ib = iblo, ibhi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3abd(ib,id,ii) = d3abd(ib,id,ii) 
     &                                  + nbia(ib,jj)*ndia(id,jj)*
     &                                    occija(ii,jj)
                      enddo
                    enddo
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bbd(ib,id,ii) = d3bbd(ib,id,ii) 
     &                                  + nbib(ib,jj)*ndib(id,jj)*
     &                                    occijb(ii,jj)
                      enddo
                    enddo
                  enddo
                  do ia = ialo, iahi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3aad(ia,id,ii) = d3aad(ia,id,ii) 
     &                                  + naia(ia,jj)*ndia(id,jj)*
     &                                    occija(ii,jj)
                      enddo
                    enddo
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bad(ia,id,ii) = d3bad(ia,id,ii) 
     &                                  + naib(ib,jj)*ndib(id,jj)*
     &                                    occijb(ii,jj)
                      enddo
                    enddo
                  enddo
                enddo
                do ic = iclo, ichi
                  do ib = iblo, ibhi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3abc(ib,ic,ii) = d3abc(ib,ic,ii) 
     &                                  + nbia(ib,jj)*ncia(ic,jj)*
     &                                    occija(ii,jj)
                      enddo
                    enddo
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bbc(ib,ic,ii) = d3bbc(ib,ic,ii) 
     &                                  + nbib(ib,jj)*ncib(ic,jj)*
     &                                    occijb(ii,jj)
                      enddo
                    enddo
                  enddo
                  do ia = ialo, iahi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3aac(ia,ic,ii) = d3aac(ia,ic,ii) 
     &                                  + naia(ia,jj)*ncia(ic,jj)*
     &                                    occija(ii,jj)
                      enddo
                    enddo
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bac(ia,ic,ii) = d3bac(ia,ic,ii) 
     &                                  + naib(ib,jj)*ncib(ic,jj)*
     &                                    occijb(ii,jj)
                      enddo
                    enddo
                  enddo
                enddo
                do id = idlo, idhi
                  do ib = iblo, ibhi
                    do ic = iclo, ichi
                      do ia = ialo, iahi
                        do ii = iialo, iiahi
                          den(ia,ic,ib,id) = den(ia,ic,ib,id) 
     &                                     + naia(ia,ii)*ncia(ic,ii)*
     &                                       d3abd(ib,id,ii)*0.5d0
     &                                     + nbia(ib,ii)*ndia(id,ii)*
     &                                       d3aac(ia,ic,ii)*0.5d0
     &                                     - nbia(ib,ii)*ncia(ic,ii)*
     &                                       d3aad(ia,id,ii)*0.5d0
     &                                     - naia(ia,ii)*ndia(id,ii)*
     &                                       d3abc(ib,ic,ii)*0.5d0
                        enddo ! ii
                        do ii = iiblo, iibhi
                          den(ia,ic,ib,id) = den(ia,ic,ib,id) 
     &                                     + naib(ia,ii)*ncib(ic,ii)*
     &                                       d3bbd(ib,id,ii)*0.5d0
     &                                     + nbib(ib,ii)*ndib(id,ii)*
     &                                       d3bac(ia,ic,ii)*0.5d0
     &                                     - nbib(ib,ii)*ncib(ic,ii)*
     &                                       d3bad(ia,id,ii)*0.5d0
     &                                     - naib(ia,ii)*ndib(id,ii)*
     &                                       d3bbc(ib,ic,ii)*0.5d0
                        enddo ! ii
                        do ii = iialo, iiahi
                          fn2ica(ii,ic) = fn2ica(ii,ic) 
     &                                  + naia(ia,ii)*d3abd(ib,id,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
     &                                  - nbia(ib,ii)*d3aad(ia,id,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
                          fn2ida(ii,id) = fn2ida(ii,id) 
     &                                  - naia(ia,ii)*d3abc(ib,ic,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
     &                                  + nbia(ib,ii)*d3aac(ia,ic,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
                        enddo ! ii
                        do ii = iiblo, iibhi
                          fn2icb(ii,ic) = fn2icb(ii,ic) 
     &                                  + naib(ia,ii)*d3bbd(ib,id,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
     &                                  - nbib(ib,ii)*d3bad(ia,id,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
                          fn2idb(ii,id) = fn2idb(ii,id) 
     &                                  - naib(ia,ii)*d3bbc(ib,ic,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
     &                                  + nbib(ib,ii)*d3bac(ia,ic,ii)*
     &                                    eri(ia,ic,ib,id)*0.5d0
                        enddo ! ii
                      enddo ! ia
                    enddo ! ic
                  enddo ! ib
                enddo ! id
                do jj = 1, nmo
                  do ii = iialo, iiahi
                    do ic = iclo, ichi
                      fna(ii,jj) = fna(ii,jj) 
     &                           + 1.0d0*fn2ica(ii,ic)*ncia(ic,jj)
                    enddo
                    do id = idlo, idhi
                      fna(ii,jj) = fna(ii,jj) 
     &                           + 1.0d0*fn2ida(ii,id)*ndia(id,jj)
                    enddo
                  enddo
                enddo
                do jj = 1, nmo
                  do ii = iiblo, iibhi
                    do ic = iclo, ichi
                      fnb(ii,jj) = fnb(ii,jj)
     &                           + 1.0d0*fn2icb(ii,ic)*ncib(ic,jj)
                    enddo
                    do id = idlo, idhi
                      fnb(ii,jj) = fnb(ii,jj)
     &                           + 1.0d0*fn2idb(ii,id)*ndib(id,jj)
                    enddo
                  enddo
                enddo
                call ga_acc(wfn1_deriv%g_dna,iialo,iiahi,1,nmo,
     &                      fna,nia,1.0d0)
                call ga_acc(wfn1_deriv%g_dnb,iiblo,iibhi,1,nmo,
     &                      fnb,nib,1.0d0)
              enddo ! jproc
c
c             Add the 2-electron energy
c
              E_2el = E_2el + 0.5d0*ddot(na*nb*nc*nd,eri,1,den,1)
c
            enddo ! la
          enddo ! lb
        enddo ! lc
      enddo ! ld
c
c     Clean up memory
c     ---------------
c
      deallocate(buf3aac, buf3bac)
      deallocate(buf3aad, buf3bad)
      deallocate(buf3abc, buf3bbc)
      deallocate(buf3abd, buf3bbd)
      deallocate(bufxbda, bufxbdb)
      deallocate(bufbda,  bufbdb)
      deallocate(bufaca,  bufacb)
      deallocate(bufija,  bufijb)
      deallocate(bufcb,   bufdb)
      deallocate(bufab,   bufbb)
      deallocate(bufca,   bufda)
      deallocate(bufaa,   bufba)
      deallocate(buf,     buf2, scr)
      deallocate(ocxa,    ocxb)
      deallocate(occa,    occb)
c
c     Global sums
c     -----------
c
      call ga_sync()
      call ga_dgop(msg_null,E_2el,1,'+')
c
      end subroutine wfn1_e_g_2el
c
c-----------------------------------------------------------------------
