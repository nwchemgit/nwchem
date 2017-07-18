c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_2el(wfn1_param,wfn1_inst,
     &                        wfn1_wave,wfn1_deriv,
     &                        v_dca,v_dcb,g_2aa,g_2bb,
     &                        g_x2aa,g_x2bb,E_2el,E_2el_hf)
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
      integer, intent(inout) :: g_x2aa  ! nmo x nmo matrix for the
                                        ! exchange term derivative of
                                        ! 2-electron density matrix
                                        ! aa-block
      integer, intent(inout) :: g_x2bb  ! nmo x nmo matrix for the
                                        ! exchange term derivative of
                                        ! 2-electron density matrix
                                        ! bb-block
c
      double precision, intent(out) :: E_2el    ! the two-electron energy
      double precision, intent(out) :: E_2el_hf ! the two-electron energy
                                                ! of conventional Hartree-Fock
c
c     Local
c
      double precision, external :: ddot
c
      double precision, allocatable :: occa(:)  ! a-occupation numbers
      double precision, allocatable :: occb(:)  ! b-occupation numbers
c
c     double precision, allocatable :: ocxa(:)  ! a-average exchange 
c                                               ! occupation numbers
c                                               ! of 2RDM
c     double precision, allocatable :: ocxb(:)  ! b-average exchange
c                                               ! occupation numbers
c                                               ! of 2RDM
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
      double precision, pointer :: buf_hf(:)
      double precision, pointer :: den_hf(:,:,:,:)
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
      double precision, pointer :: bufx3aad(:)
      double precision, pointer :: dx3aad(:,:,:)
      double precision, pointer :: bufx3abd(:)
      double precision, pointer :: dx3abd(:,:,:)
      double precision, pointer :: bufx3bad(:)
      double precision, pointer :: dx3bad(:,:,:)
      double precision, pointer :: bufx3bbd(:)
      double precision, pointer :: dx3bbd(:,:,:)
c
      double precision, pointer :: bufx3aac(:)
      double precision, pointer :: dx3aac(:,:,:)
      double precision, pointer :: bufx3abc(:)
      double precision, pointer :: dx3abc(:,:,:)
      double precision, pointer :: bufx3bac(:)
      double precision, pointer :: dx3bac(:,:,:)
      double precision, pointer :: bufx3bbc(:)
      double precision, pointer :: dx3bbc(:,:,:)
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
      double precision, pointer :: bufxija(:)
      double precision, pointer :: ocxija(:,:)
      double precision, pointer :: bufxijb(:)
      double precision, pointer :: ocxijb(:,:)
c
      double precision, pointer :: bufaca(:)
      double precision, pointer :: daca(:,:) ! a,c 1-RDM alpha
      double precision, pointer :: bufacb(:)
      double precision, pointer :: dacb(:,:) ! a,c 1-RDM beta
      double precision, pointer :: bufbda(:)
      double precision, pointer :: dbda(:,:) ! b,d 1-RDM alpha
      double precision, pointer :: bufbdb(:)
      double precision, pointer :: dbdb(:,:) ! b,d 1-RDM beta
c
      double precision, pointer :: bufada(:)
      double precision, pointer :: dada(:,:) ! a,d 1-RDM alpha
      double precision, pointer :: bufadb(:)
      double precision, pointer :: dadb(:,:) ! a,d 1-RDM beta
      double precision, pointer :: bufbca(:)
      double precision, pointer :: dbca(:,:) ! b,c 1-RDM alpha
      double precision, pointer :: bufbcb(:)
      double precision, pointer :: dbcb(:,:) ! b,c 1-RDM beta
c
c     double precision, pointer :: bufxbda(:)
c     double precision, pointer :: xbda(:,:) ! b,d 1-RDM alpha
c     double precision, pointer :: bufxbdb(:)
c     double precision, pointer :: xbdb(:,:) ! b,d 1-RDM beta
c
      double precision, pointer :: buffnaa(:)
      double precision, pointer :: buffnab(:)
      double precision, pointer :: buffnac(:)
      double precision, pointer :: buffnad(:)
      double precision, pointer :: fnac(:,:)
      double precision, pointer :: fn1aca(:,:) ! a,c 1-electron natural
                                               !     orbital Fock alpha
      double precision, pointer :: fn1ica(:,:) ! i,c 1-electron natural
                                               !     orbital Fock alpha
      double precision, pointer :: fn2ica(:,:) ! i,c 2-electron natural
                                               !     orbital Fock alpha
      double precision, pointer :: fn2ida(:,:) ! i,d 2-electron natural
                                               !     orbital Fock alpha
      double precision, pointer :: fc1aca(:,:) ! a,c 1-electron
                                               !     correlation function
                                               !     Fock alpha
      double precision, pointer :: fc1ada(:,:) ! a,d 1-electron
                                               !     correlation function
                                               !     Fock alpha
      double precision, pointer :: fc1bca(:,:) ! b,c 1-electron
                                               !     correlation function
                                               !     Fock alpha
      double precision, pointer :: fc1bda(:,:) ! b,d 1-electron
                                               !     correlation function
                                               !     Fock alpha
      double precision, pointer :: buffnba(:)
      double precision, pointer :: buffnbb(:)
      double precision, pointer :: buffnbc(:)
      double precision, pointer :: buffnbd(:)
      double precision, pointer :: fnbc(:,:)
      double precision, pointer :: fn1acb(:,:) ! a,c 1-electron natural
                                               !     orbital Fock beta
      double precision, pointer :: fn1icb(:,:) ! i,c 1-electron natural
                                               !     orbital Fock beta
      double precision, pointer :: fn2icb(:,:) ! i,c 2-electron natural
                                               !     orbital Fock beta
      double precision, pointer :: fn2idb(:,:) ! i,d 2-electron natural
                                               !     orbital Fock beta
      double precision, pointer :: fc1acb(:,:) ! a,c 1-electron
                                               !     correlation function
                                               !     Fock beta
      double precision, pointer :: fc1adb(:,:) ! a,d 1-electron
                                               !     correlation function
                                               !     Fock beta
      double precision, pointer :: fc1bcb(:,:) ! b,c 1-electron
                                               !     correlation function
                                               !     Fock beta
      double precision, pointer :: fc1bdb(:,:) ! b,d 1-electron
                                               !     correlation function
                                               !     Fock beta
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
      E_2el_hf = 0.0d0
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
c     allocate(ocxa(1:nmo),ocxb(1:nmo))
      allocate(buf(max2e),buf2(max2e),buf_hf(max2e),scr(mem2))
      allocate(bufaa(maxsbf*nmo),bufba(maxsbf*nmo))
      allocate(bufca(maxsbf*nmo),bufda(maxsbf*nmo))
      allocate(bufab(maxsbf*nmo),bufbb(maxsbf*nmo))
      allocate(bufcb(maxsbf*nmo),bufdb(maxsbf*nmo))
      allocate(bufija(nia*nja),bufijb(nib*njb))
      allocate(bufxija(nia*nja),bufxijb(nib*njb))
      allocate(bufaca(maxsbf*maxsbf),bufacb(maxsbf*maxsbf))
      allocate(bufada(maxsbf*maxsbf),bufadb(maxsbf*maxsbf))
      allocate(bufbca(maxsbf*maxsbf),bufbcb(maxsbf*maxsbf))
      allocate(bufbda(maxsbf*maxsbf),bufbdb(maxsbf*maxsbf))
c     allocate(bufxbda(maxsbf*maxsbf),bufxbdb(maxsbf*maxsbf))
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
      allocate(bufx3aac(maxsbf*maxsbf*max(nia,nja)))
      allocate(bufx3aad(maxsbf*maxsbf*max(nia,nja)))
      allocate(bufx3abc(maxsbf*maxsbf*max(nia,nja)))
      allocate(bufx3abd(maxsbf*maxsbf*max(nia,nja)))
      allocate(bufx3bac(maxsbf*maxsbf*max(nib,njb)))
      allocate(bufx3bad(maxsbf*maxsbf*max(nib,njb)))
      allocate(bufx3bbc(maxsbf*maxsbf*max(nib,njb)))
      allocate(bufx3bbd(maxsbf*maxsbf*max(nib,njb)))
c
      allocate(buffnaa(max(nia,nja,maxsbf)*nmo))
      allocate(buffnab(max(nia,nja,maxsbf)*nmo))
      allocate(buffnac(max(nia,nja,maxsbf)*nmo))
      allocate(buffnad(max(nia,nja,maxsbf)*nmo))
      allocate(buffnba(max(nib,njb,maxsbf)*nmo))
      allocate(buffnbb(max(nib,njb,maxsbf)*nmo))
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
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
c     call wfn1_calc_ocx(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,ocxa)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
c     call wfn1_calc_ocx(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,ocxb)
c
c     Calculate alpha-alpha-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                    wfn1_wave%dcnta,g_2aa)
c
c     Calculate beta-beta-occupation numbers
c
      call wfn1_calc_occ2(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                    wfn1_wave%dcntb,g_2bb)
c
c     Calculate exchange alpha-alpha-occupation number derivative
c
      call wfn1_calc_docx2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                     wfn1_wave%dcnta,g_x2aa)
c
c     Calculate exchange beta-beta-occupation number derivative
c
      call wfn1_calc_docx2(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                     wfn1_wave%dcntb,g_x2bb)
c
      call ga_sync()
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
              den_hf(ialo:iahi,iclo:ichi,iblo:ibhi,idlo:idhi)=>buf_hf
c
c             calculate the 2-electron integrals
c
              call int_2e4c(wfn1_inst%basis, ld, lb,
     &                      wfn1_inst%basis, lc, la,
     &                      mem2, scr, max2e, eri)
              den = 0.0d0
              den_hf = 0.0d0
c
c             calculate the 1RDM terms
c
              fn1aca(ialo:iahi,iclo:ichi)=>buffnad
              fn1acb(ialo:iahi,iclo:ichi)=>buffnbd
              fn1aca = 0.0d0
              fn1acb = 0.0d0
              daca(ialo:iahi,iclo:ichi)=>bufaca
              dacb(ialo:iahi,iclo:ichi)=>bufacb
              dada(ialo:iahi,idlo:idhi)=>bufada
              dadb(ialo:iahi,idlo:idhi)=>bufadb
              dbca(iblo:ibhi,iclo:ichi)=>bufbca
              dbcb(iblo:ibhi,iclo:ichi)=>bufbcb
              dbda(iblo:ibhi,idlo:idhi)=>bufbda
              dbdb(iblo:ibhi,idlo:idhi)=>bufbdb
c             xbda(iblo:ibhi,idlo:idhi)=>bufxbda
c             xbdb(iblo:ibhi,idlo:idhi)=>bufxbdb
              daca = 0.0d0
              dacb = 0.0d0
              dada = 0.0d0
              dadb = 0.0d0
              dbca = 0.0d0
              dbcb = 0.0d0
              dbda = 0.0d0
              dbdb = 0.0d0
c             xbda = 0.0d0
c             xbdb = 0.0d0
c
c             Construct the 1RDM and 2RDM density matrix terms
c
              do ic = iclo, ichi
                do ib = iblo, ibhi
                  do ii = 1, nmo
                    dbca(ib,ic) = dbca(ib,ic)
     &                          + nbia(ib,ii)*ncia(ic,ii)*occa(ii)
                    dbcb(ib,ic) = dbcb(ib,ic)
     &                          + nbib(ib,ii)*ncib(ic,ii)*occb(ii)
                  enddo
                enddo
              enddo
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ii = 1, nmo
                    dbda(ib,id) = dbda(ib,id)
     &                          + nbia(ib,ii)*ndia(id,ii)*occa(ii)
                    dbdb(ib,id) = dbdb(ib,id)
     &                          + nbib(ib,ii)*ndib(id,ii)*occb(ii)
c                   xbda(ib,id) = xbda(ib,id)
c    &                          + nbia(ib,ii)*ndia(id,ii)*ocxa(ii)
c                   xbdb(ib,id) = xbdb(ib,id)
c    &                          + nbib(ib,ii)*ndib(id,ii)*ocxb(ii)
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
                do ia = ialo, iahi
                  do ii = 1, nmo
                    dada(ia,id) = dada(ia,id)
     &                          + naia(ia,ii)*ndia(id,ii)*occa(ii)
                    dadb(ia,id) = dadb(ia,id)
     &                          + naib(ia,ii)*ndib(id,ii)*occb(ii)
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
                      den_hf(ia,ic,ib,id) = den_hf(ia,ic,ib,id) 
     &                                    + daca(ia,ic)*dbda(ib,id)
     &                                    + daca(ia,ic)*dbdb(ib,id)
     &                                    + dacb(ia,ic)*dbda(ib,id)
     &                                    + dacb(ia,ic)*dbdb(ib,id)
     &                                    - dada(ia,id)*dbca(ib,ic)
     &                                    - dadb(ia,id)*dbcb(ib,ic)
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
              fc1aca(ialo:iahi,iclo:ichi)=>buffnac
              fc1ada(ialo:iahi,idlo:idhi)=>buffnad
              fc1bca(iblo:ibhi,iclo:ichi)=>buffnaa
              fc1bda(iblo:ibhi,idlo:idhi)=>buffnab
              fc1acb(ialo:iahi,iclo:ichi)=>buffnbc
              fc1adb(ialo:iahi,idlo:idhi)=>buffnbd
              fc1bcb(iblo:ibhi,iclo:ichi)=>buffnba
              fc1bdb(iblo:ibhi,idlo:idhi)=>buffnbb
              fc1aca = 0.0d0
              fc1ada = 0.0d0
              fc1bca = 0.0d0
              fc1bda = 0.0d0
              fc1acb = 0.0d0
              fc1adb = 0.0d0
              fc1bcb = 0.0d0
              fc1bdb = 0.0d0
              do id = idlo, idhi
                do ib = iblo, ibhi
                  do ic = iclo, ichi
                    do ia = ialo, iahi
c
c                     Fock matrix Coulomb terms a-a and a-b
c
                      fc1aca(ia,ic) = fc1aca(ia,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.50d0*dbda(ib,id)+
     &                                 0.50d0*dbdb(ib,id))
                      fc1ada(ia,id) = fc1ada(ia,id)
     &                              + eri(ia,ic,ib,id)*
     &                                0.50d0*dbca(ib,ic)
                      fc1bca(ib,ic) = fc1bca(ib,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                0.50d0*dada(ia,id)
                      fc1bda(ib,id) = fc1bda(ib,id)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.50d0*daca(ia,ic)+
     &                                 0.50d0*dacb(ia,ic))
c
c                     Fock matrix Coulomb terms b-b and b-a
c
                      fc1acb(ia,ic) = fc1acb(ia,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.5d0*dbdb(ib,id)+
     &                                 0.5d0*dbda(ib,id))
                      fc1adb(ia,id) = fc1adb(ia,id)
     &                              + eri(ia,ic,ib,id)*
     &                                0.5d0*dbcb(ib,ic)
                      fc1bcb(ib,ic) = fc1bcb(ib,ic)
     &                              + eri(ia,ic,ib,id)*
     &                                0.5d0*dadb(ia,id)
                      fc1bdb(ib,id) = fc1bdb(ib,id)
     &                              + eri(ia,ic,ib,id)*
     &                                (0.5d0*dacb(ia,ic)+
     &                                 0.5d0*daca(ia,ic))
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
                do id = idlo, idhi
                  do ia = ialo, iahi
                    v_dca(ii) = v_dca(ii)
     &                        - 1.0d0*naia(ia,ii)*
     &                          fc1ada(ia,id)*ndia(id,ii)
                    v_dcb(ii) = v_dcb(ii)
     &                        - 1.0d0*naib(ia,ii)*
     &                          fc1adb(ia,id)*ndib(id,ii)
                  enddo
                enddo
                do ic = iclo, ichi
                  do ib = iblo, ibhi
                    v_dca(ii) = v_dca(ii)
     &                        - 1.0d0*nbia(ib,ii)*
     &                          fc1bca(ib,ic)*ncia(ic,ii)
                    v_dcb(ii) = v_dcb(ii)
     &                        - 1.0d0*nbib(ib,ii)*
     &                          fc1bcb(ib,ic)*ncib(ic,ii)
                  enddo
                enddo
                do id = idlo, idhi
                  do ib = iblo, ibhi
                    v_dca(ii) = v_dca(ii)
     &                        + 1.0d0*nbia(ib,ii)*
     &                          fc1bda(ib,id)*ndia(id,ii)
                    v_dcb(ii) = v_dcb(ii)
     &                        + 1.0d0*nbib(ib,ii)*
     &                          fc1bdb(ib,id)*ndib(id,ii)
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
                ocxija(iialo:iiahi,ijalo:ijahi)=>bufxija
                ocxijb(iiblo:iibhi,ijblo:ijbhi)=>bufxijb
                call ga_get(g_2aa, iialo,iiahi,ijalo,ijahi,occija,nia)
                call ga_get(g_x2aa,iialo,iiahi,ijalo,ijahi,ocxija,nia)
                call ga_get(g_2bb, iiblo,iibhi,ijblo,ijbhi,occijb,nib)
                call ga_get(g_x2bb,iiblo,iibhi,ijblo,ijbhi,ocxijb,nib)
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
                dx3aac(ialo:iahi,iclo:ichi,iialo:iiahi)=>bufx3aac
                dx3aad(ialo:iahi,idlo:idhi,iialo:iiahi)=>bufx3aad
                dx3abc(iblo:ibhi,iclo:ichi,iialo:iiahi)=>bufx3abc
                dx3abd(iblo:ibhi,idlo:idhi,iialo:iiahi)=>bufx3abd
                dx3bac(ialo:iahi,iclo:ichi,iiblo:iibhi)=>bufx3bac
                dx3bad(ialo:iahi,idlo:idhi,iiblo:iibhi)=>bufx3bad
                dx3bbc(iblo:ibhi,iclo:ichi,iiblo:iibhi)=>bufx3bbc
                dx3bbd(iblo:ibhi,idlo:idhi,iiblo:iibhi)=>bufx3bbd
                dx3aac = 0.0d0
                dx3aad = 0.0d0
                dx3abc = 0.0d0
                dx3abd = 0.0d0
                dx3bac = 0.0d0
                dx3bad = 0.0d0
                dx3bbc = 0.0d0
                dx3bbd = 0.0d0
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
                        d3abd(ib,id,ii)  = d3abd(ib,id,ii) 
     &                                   + nbia(ib,jj)*ndia(id,jj)*
     &                                     occija(ii,jj)
                        dx3abd(ib,id,ii) = dx3abd(ib,id,ii) 
     &                                   + nbia(ib,jj)*ndia(id,jj)*
     &                                     ocxija(ii,jj)
                      enddo ! jj
                    enddo ! ii
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bbd(ib,id,ii)  = d3bbd(ib,id,ii) 
     &                                   + nbib(ib,jj)*ndib(id,jj)*
     &                                     occijb(ii,jj)
                        dx3bbd(ib,id,ii) = dx3bbd(ib,id,ii) 
     &                                   + nbib(ib,jj)*ndib(id,jj)*
     &                                     ocxijb(ii,jj)
                      enddo ! jj
                    enddo ! ii
                  enddo ! ib
                  do ia = ialo, iahi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3aad(ia,id,ii)  = d3aad(ia,id,ii) 
     &                                   + naia(ia,jj)*ndia(id,jj)*
     &                                     occija(ii,jj)
                        dx3aad(ia,id,ii) = dx3aad(ia,id,ii) 
     &                                   + naia(ia,jj)*ndia(id,jj)*
     &                                     ocxija(ii,jj)
                      enddo ! jj
                    enddo ! ii
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bad(ia,id,ii)  = d3bad(ia,id,ii) 
     &                                   + naib(ia,jj)*ndib(id,jj)*
     &                                     occijb(ii,jj)
                        dx3bad(ia,id,ii) = dx3bad(ia,id,ii) 
     &                                   + naib(ia,jj)*ndib(id,jj)*
     &                                     ocxijb(ii,jj)
                      enddo ! jj
                    enddo ! ii
                  enddo ! ia
                enddo ! id
                do ic = iclo, ichi
                  do ib = iblo, ibhi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3abc(ib,ic,ii)  = d3abc(ib,ic,ii) 
     &                                   + nbia(ib,jj)*ncia(ic,jj)*
     &                                     occija(ii,jj)
                        dx3abc(ib,ic,ii) = dx3abc(ib,ic,ii) 
     &                                   + nbia(ib,jj)*ncia(ic,jj)*
     &                                     ocxija(ii,jj)
                      enddo ! jj
                    enddo ! ii
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bbc(ib,ic,ii)  = d3bbc(ib,ic,ii) 
     &                                   + nbib(ib,jj)*ncib(ic,jj)*
     &                                     occijb(ii,jj)
                        dx3bbc(ib,ic,ii) = dx3bbc(ib,ic,ii) 
     &                                   + nbib(ib,jj)*ncib(ic,jj)*
     &                                     ocxijb(ii,jj)
                      enddo ! jj
                    enddo ! ii
                  enddo ! ib
                  do ia = ialo, iahi
                    do ii = iialo, iiahi
                      do jj = ijalo, ijahi
                        d3aac(ia,ic,ii)  = d3aac(ia,ic,ii) 
     &                                   + naia(ia,jj)*ncia(ic,jj)*
     &                                     occija(ii,jj)
                        dx3aac(ia,ic,ii) = dx3aac(ia,ic,ii) 
     &                                   + naia(ia,jj)*ncia(ic,jj)*
     &                                     ocxija(ii,jj)
                      enddo ! jj
                    enddo ! ii
                    do ii = iiblo, iibhi
                      do jj = ijblo, ijbhi
                        d3bac(ia,ic,ii)  = d3bac(ia,ic,ii) 
     &                                   + naib(ia,jj)*ncib(ic,jj)*
     &                                     occijb(ii,jj)
                        dx3bac(ia,ic,ii) = dx3bac(ia,ic,ii) 
     &                                   + naib(ia,jj)*ncib(ic,jj)*
     &                                     ocxijb(ii,jj)
                      enddo ! jj
                    enddo ! ii
                  enddo ! ia
                enddo ! ic
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
                        do ii = iialo, iiahi
                          v_dca(ii) = v_dca(ii)
     &                              - 0.5d0*eri(ia,ic,ib,id)*
     &                                (
     &                                 dx3aac(ia,ic,ii)*nbia(ib,ii)*
     &                                 ndia(id,ii)
     &                                +dx3abd(ib,id,ii)*naia(ia,ii)*
     &                                 ncia(ic,ii)
     &                                -dx3aad(ia,id,ii)*nbia(ib,ii)*
     &                                 ncia(ic,ii)
     &                                -dx3abc(ib,ic,ii)*naia(ia,ii)*
     &                                 ndia(id,ii)
     &                                )
     
                        enddo ! ii
                        do ii = iiblo, iibhi
                          v_dcb(ii) = v_dcb(ii)
     &                              - 0.5d0*eri(ia,ic,ib,id)*
     &                                (
     &                                 dx3bac(ia,ic,ii)*nbib(ib,ii)*
     &                                 ndib(id,ii)
     &                                +dx3bbd(ib,id,ii)*naib(ia,ii)*
     &                                 ncib(ic,ii)
     &                                -dx3bad(ia,id,ii)*nbib(ib,ii)*
     &                                 ncib(ic,ii)
     &                                -dx3bbc(ib,ic,ii)*naib(ia,ii)*
     &                                 ndib(id,ii)
     &                                )
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
              E_2el_hf = E_2el_hf+0.5d0*ddot(na*nb*nc*nd,eri,1,den_hf,1)
c
            enddo ! la
          enddo ! lb
        enddo ! lc
      enddo ! ld
c
c     Clean up memory
c     ---------------
c
      deallocate(buffna,   buffnb)
      deallocate(buffnad,  buffnbd)
      deallocate(buffnac,  buffnbc)
      deallocate(buffnab,  buffnbb)
      deallocate(buffnaa,  buffnba)
c
      deallocate(bufx3aac, bufx3bac)
      deallocate(bufx3aad, bufx3bad)
      deallocate(bufx3abc, bufx3bbc)
      deallocate(bufx3abd, bufx3bbd)
c
      deallocate(buf3aac,  buf3bac)
      deallocate(buf3aad,  buf3bad)
      deallocate(buf3abc,  buf3bbc)
      deallocate(buf3abd,  buf3bbd)
c     deallocate(bufxbda,  bufxbdb)
c
      deallocate(bufbda,   bufbdb)
      deallocate(bufbca,   bufbcb)
      deallocate(bufada,   bufadb)
      deallocate(bufaca,   bufacb)
c
      deallocate(bufxija,  bufxijb)
      deallocate(bufija,   bufijb)
c
      deallocate(bufcb,    bufdb)
      deallocate(bufab,    bufbb)
      deallocate(bufca,    bufda)
      deallocate(bufaa,    bufba)
c
      deallocate(buf,      buf2, buf_hf, scr)
c     deallocate(ocxa,     ocxb)
      deallocate(occa,     occb)
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
