c-----------------------------------------------------------------------
c
!> \ingroup wfn1
!> @{
!>
!> \file wfn1_e_g_2el_hf.F
!> The conventional Hartree-Fock 2-electron contributions
!>
!> @}
!>
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Evaluate the conventional Hartree-Fock 2-electron contributions
!>
!> This subroutine evaluates the 2-electron contributions of the
!> conventional Hartree-Fock formalism to the energy, the
!> natural orbital Fock matrices, and the correlation function Fock
!> matrices.
!>
!> Within the Hartree-Fock formalism the Coulomb terms can be
!> constructed from regular 1RDMs. The exchange terms require the
!> explicit evaluation of the 2RDM diagonal. The latter form 
!> \f$O(N^2)\f$ data objects and are therefore distributed over 
!> multiple nodes. 
!>
!> Note that naive evaluation of the 2RDM would have a complexity
!> of \f$O(N^6)\f$ as is evident from the expression
!> \f{eqnarray*}{
!>   \Gamma_{abcd} &=& -\sum_{i,j}N_{ai}N^*_{di}N_{bj}N^*_{cj}d_{ij}
!> \f}
!> where all indices \f$a, b, c, d, i\f$ and \f$j\f$ run over all
!> basis functions. This complexity can be reduced to \f$O(N^5)\f$
!> by breaking the equations down into two steps
!> \f{eqnarray*}{
!>   \gamma_{adj} &=& -\sum_{i}N_{ai}N^*_{di}d_{ij} \\\\
!>   \Gamma_{abcd} &=& \sum_{j}\gamma_{adj}N_{bj}N^*_{cj}
!> \f}
!> Note that this approach is still communication intensive as typically
!> the loops over \f$a\f$ and \f$d\f$ would be partitioned into blocks,
!> and for every block of \f$\gamma_{adj}\f$ all \f$O(N^2)\f$ elements
!> of the diagonal \f$d_{ij}\f$ have to be accessed. Hence formally
!> \f$O(N^4)\f$ data needs to be communicated.
!>
!> Ultimately, the effective complexities can be reduced using, for
!> example, the Cauchy-Schwarz inequality [1-4]. For large molecules
!> only \f$O(N^2)\f$ of the \f$a, b, c, d\f$ indeces, only \f$O(N)\f$
!> of the \f$a, d\f$ indeces would survive. This implies that both the
!> computational as well as the communication complexity are reduced to
!> \f$O(N^3)\f$. At present the Cauchy-Schwarz inequality is not used.
!>
!> On the memory management aspects there are costs associated with
!> allocating memory on the heap as in the Fortran ALLOCATE statement.
!> To avoid these overheads, in particular in the inner loops, Fortran
!> POINTERs are used. This means that the actual memory is associated
!> with buffers allocated outside the loops. The pointers are set 
!> inside the loops to provide convenient array views on the memory. 
!> For clarity pairs of variables are used. The ALLOCATABLE arrays
!> having names BUF_<name> are the buffers that hold the actual memory.
!> Associated with each such array is a POINTER <name> that points to
!> BUF_<name> and provides an appropriate number of dimensions and
!> corresponding index limits.
!>
!> In addition this subroutine inherently deals with \f$\alpha\$ and
!> \f$\beta\f$-spin variables. Variables related to \f$\alpha\f$-spin
!> electrons are identified by names with a trailing "_a" whereas 
!> \f$\beta\f$-electron variables are identified by names with a
!> trailing "_b".
!>
!> Another naming convention is that variables include the remaining
!> indeces. For example \f$\gamma_{adj}\f$ given above would have a
!> name that lists the remaining indeces \f$a, d, j\f$ and the order
!> of these indeces corresponds to the order of the dimensions of the
!> 3-dimensional tensor.
!>
!> This subroutine computes the 2-electron energy term in two 
!> different ways. Once in accordance to the N-representable
!> Hartree-Fock formalism [5], and once, for comparison, according to
!> the conventional Hartree-Fock method [6-8].
!>
!> ### References ###
!>
!> [1] A.-L. Cauchy, "Sur les formules qui résultent de l'emploi du
!>     signe > ou <, et sur les moyennes entre plusieurs quantités",
!>     In Book: Cours d'analyse de l'École Royale Polytechnique, (1821).
!>     Reprinted in Cambridge Library Collection - Mathematics,
!>     Cambridge University Press (2009), pp. 438-459, DOI:
!>     <a href="https://dx.doi.org/10.1017/CBO9780511693328.017>
!>     10.1017/CBO9780511693328.017</a>.
!>
!> [2] V. Bouniakowsky, "Sur quelques inegalités concernant les
!>     intégrales aux différences finies", Memoires de L'Academie
!>     Imperiale des Sciences de St.-Petersbourg, VII, 1 (1859) Nr. 9,
!>     pp. 1-18, URL:
!>     <a href="https://books.google.com/books?id=wZhFAAAAcAAJ>
!>     books.google.com/books?id=wZhFAAAAcAAJ</a>.
!>
!> [3] H. A. Schwarz, "Über ein Flächen kleinsten Flächeninhalts
!>     betreffendes Problem der Variationsrechnung", Acta Societatis
!>     Scientiarum Fennicae, XV: 318 (1885), URL:
!>     <a href="https://books.google.com/books?id=BXtIAQAAMAAJ">
!>     https://books.google.com/books?id=BXtIAQAAMAAJ</a>.
!>
!> [4] J. M. Steele, "The Cauchy-Schwarz Master Class", The Mathematical
!>     Association of America (2004), ISBN: 
!>     <a href="https://isbnsearch.org/isbn/052154677X">
!>     0-521-54677-X</a>.
!>
!> [5] H. J. J. van Dam, "An N-representable Hartree-Fock approach for
!>     non-idempotent density matrixes", in preparation (08/31/2017).
!>
!> [6] C. C. J. Roothaan, "New developments in molecular orbital
!>     theory", Reviews of Modern Physics <b>23</b> (1951) pp. 69-89,
!>     DOI: <a href="https://dx.doi.org/10.1103/RevModPhys.23.69">
!>     10.1103/RevModPhys.23.69</a>.
!>
!> [7] G. G. Hall, "The Molecular Orbital Theory of Chemical Valency.
!>     VIII. A Method of Calculating Ionization Potentials", 
!>     Proceedings of the Royal Society A <b>205</b> (1951) pp. 541-552,
!>     DOI: <a href="https://dx.doi.org/10.1098/rspa.1951.0048">
!>     10.1098/rspa.1951.0048</a>.
!>
!> [8] J. A. Pople, R. K. Nesbet, "Self-consistent orbitals for
!>     radicals", The Journal of Chemical Physics <b>22</b> (1954)
!>     pp. 571-572,  DOI:
!>     <a href="https://dx.doi.org/10.1063/1.1740120">
!>     10.1063/1.1740120</a>.
!>
!> ### TO DO ###
!>
!> 1. Work out whether a given processor needs to execute anything
!>    for a given index of ld, lc, lb, or la and skip empty loops.
!>    At the moment a simple task counter is checked in the inner most
!>    loop. Skipping loop bodies in outer loops can bring further 
!>    efficiency improvements.
!>
!> 2. It seems like in some places consecutive jj loops could be merged
!>    with the potential or reducing the order of intermediate tensors
!>    by 1. If this is actually possible then this could save memory
!>    as well as reduce memory access both of which would be beneficial.
!>
      subroutine wfn1_e_g_2el_hf(wfn1_param,wfn1_inst,wfn1_wave,
     &                           wfn1_deriv,v_dc_a,v_dc_b,
     &                           g_x_aa,g_x_bb,g_dx_aa,g_dx_bb,
     &                           E_2el,E_2el_chf)
      implicit none
c
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param ! the job parameters
      type(wfn1_prminst), intent(in) :: wfn1_inst  ! the job instances
      type(wfn1_wfn),     intent(in) :: wfn1_wave  ! the wavefunction
c
      type(wfn1_drv),  intent(inout) :: wfn1_deriv ! the derivatives
c
c     Correlation function Fock matrices in the natural orbital basis
c     alpha-, and beta-spin respectively
c
      double precision, intent(inout) :: v_dc_a(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dc_b(1:wfn1_param%nmo)
c
c     This implementation assumes that g_x_aa, g_x_bb, g_dx_aa,
c     and g_dx_bb all have the same data distribution. If the actual
c     data distributions are different then this will not affect the
c     results as these differences will be resolved by the Global
c     Arrays. However, this may generate additional communication and
c     hence reduce the efficiency.
c
c     Nmo x nmo global arrays for the 2RDM exchange diagonals.
c     See exchange terms of Eq.(34) of [5].
c
      integer, intent(in) :: g_x_aa ! exchange alpha-alpha block
      integer, intent(in) :: g_x_bb ! exchange beta-beta block
c
c     Nmo x nmo global arrays for the 2RDM exchange diagonals
c     derivatives. See the diagonal terms of Eq.(119) of [5].
c
      integer, intent(in) :: g_dx_aa ! exchange alpha-alpha block
      integer, intent(in) :: g_dx_bb ! exchange beta-beta block
c
c     The Hartree-Fock 2-electron energies
c
c     The 2-electron energy [5]:
c
      double precision, intent(out) :: E_2el ! The 2-electron energy [5]
c
c     The conventional Hartree-Fock 2-electron energy [6-8]:
c
      double precision, intent(out) :: E_2el_chf 
c
c     Local variables:
c
      integer :: nproc ! the number of processors
      integer :: iproc ! the rank of this processor
      integer :: jproc ! a counter over processors
      integer :: itask ! the current task number
c
      integer :: nmo   ! the number of MOs
      integer :: ne_a  ! the number of alpha-electrons
      integer :: ne_b  ! the number of beta-electrons
c
      integer :: ii, jj, kk, ll ! natural orbital labels
      integer :: ia, ib, ic, id ! atomic orbital labels
      integer :: la, lb, lc, ld ! atomic orbital shell labels
      integer :: na, nb, nc, nd ! number of shell functions
c
      integer :: ni, nj ! number of natural orbitals
c
      integer :: ia_lo, ia_hi ! a-index atomic orbital limits
      integer :: ib_lo, ib_hi ! b-index atomic orbital limits
      integer :: ic_lo, ic_hi ! c-index atomic orbital limits
      integer :: id_lo, id_hi ! d-index atomic orbital limits
c
      integer :: ii_lo, ii_hi ! i-index natural orbital limits
      integer :: jj_lo, jj_hi ! j-index natural orbital limits
c
      integer :: max2e  ! quartet maximum no. 2-el integrals
      integer :: mem2   ! amount of scratch space needed
      integer :: maxsbf ! maximum no. of basis functions in a shell
      integer :: nshell ! the number of shells
c
      double precision, allocatable :: occ_a(:) ! a-occupation numbers
      double precision, allocatable :: occ_b(:) ! b-occupation numbers
c
c     Data structure for the electron repulsion integrals, as well as
c     2RDM blocks for the density matrix of [5] and the conventional
c     Hartree-Fock 2RDM.
c
      double precision, allocatable         :: scr(:)
      double precision, allocatable, target :: buf_eri(:)
      double precision, pointer             :: eri(:,:,:,:)
      double precision, allocatable, target :: buf_rdm2(:)
      double precision, pointer             :: rdm2(:,:,:,:)
      double precision, allocatable, target :: buf_rdm2_chf(:)
      double precision, pointer             :: rdm2_chf(:,:,:,:)
c
c     Intermediates in the 2RDM construction
c
      double precision, allocatable, target :: buf_rdm2_adj_a(:)
      double precision, pointer             :: rdm2_adj_a(:,:,:)
      double precision, allocatable, target :: buf_rdm2_adj_b(:)
      double precision, pointer             :: rdm2_adj_b(:,:,:)
c
c     Intermediates in the derivative of the 2RDM construction
c
      double precision, allocatable, target :: buf_rdm2_d_adj_a(:)
      double precision, pointer             :: rdm2_d_adj_a(:,:,:)
      double precision, allocatable, target :: buf_rdm2_d_adj_b(:)
      double precision, pointer             :: rdm2_d_adj_b(:,:,:)
c
c     Data structures for the 1RDM blocks for the Coulomb terms of 
c     the 2RDM of [5] and the conventional Hartree-Fock 2RDM.
c
c     1RDM factors for Coulomb terms
c
      double precision, allocatable, target :: buf_rdm1_ac_a(:)
      double precision, pointer             :: rdm1_ac_a(:,:)
      double precision, allocatable, target :: buf_rdm1_ac_b(:)
      double precision, pointer             :: rdm1_ac_b(:,:)
      double precision, allocatable, target :: buf_rdm1_bd_a(:)
      double precision, pointer             :: rdm1_bd_a(:,:)
      double precision, allocatable, target :: buf_rdm1_bd_b(:)
      double precision, pointer             :: rdm1_bd_b(:,:)
c
c     1RDM factors for exchange terms of the conventional Hartree-Fock
c
      double precision, allocatable, target :: buf_rdm1_ad_a(:)
      double precision, pointer             :: rdm1_ad_a(:,:)
      double precision, allocatable, target :: buf_rdm1_ad_b(:)
      double precision, pointer             :: rdm1_ad_b(:,:)
      double precision, allocatable, target :: buf_rdm1_bc_a(:)
      double precision, pointer             :: rdm1_bc_a(:,:)
      double precision, allocatable, target :: buf_rdm1_bc_b(:)
      double precision, pointer             :: rdm1_bc_b(:,:)
c
c     Fock matrix blocks. 
c
c     For the Fock matrices for the natural orbitals
c
      double precision, allocatable, target :: buf_fn_ac_a(:)
      double precision, pointer             :: fn_ac_a(:,:)
      double precision, allocatable, target :: buf_fn_ac_b(:)
      double precision, pointer             :: fn_ac_b(:,:)
      double precision, allocatable, target :: buf_fn_bcj_a(:)
      double precision, pointer             :: fn_bcj_a(:,:,:)
      double precision, allocatable, target :: buf_fn_bcj_b(:)
      double precision, pointer             :: fn_bcj_b(:,:,:)
c
c     - Same as above but partially transformed
c
      double precision, allocatable, target :: buf_fn_ic_a(:)
      double precision, pointer             :: fn_ic_a(:,:)
      double precision, allocatable, target :: buf_fn_ic_b(:)
      double precision, pointer             :: fn_ic_b(:,:)
      double precision, allocatable, target :: buf_fn_bj_a(:)
      double precision, pointer             :: fn_bj_a(:,:)
      double precision, allocatable, target :: buf_fn_bj_b(:)
      double precision, pointer             :: fn_bj_b(:,:)
c
c     - Same as above but fully transformed
c
      double precision, allocatable, target :: buf_fn_ij_a(:)
      double precision, pointer             :: fn_ij_a(:,:)
      double precision, allocatable, target :: buf_fn_ij_b(:)
      double precision, pointer             :: fn_ij_b(:,:)
c
c     For the Fock matrices for the correlation functions
c
      double precision, allocatable, target :: buf_fc_bcj_a(:)
      double precision, pointer             :: fc_bcj_a(:,:,:)
      double precision, allocatable, target :: buf_fc_bcj_b(:)
      double precision, pointer             :: fc_bcj_b(:,:,:)
c
c     Orbital coefficients
c
      double precision, allocatable, target :: buf_N_di_a(:)
      double precision, pointer             :: N_di_a(:,:)
      double precision, allocatable, target :: buf_N_di_b(:)
      double precision, pointer             :: N_di_b(:,:)
      double precision, allocatable, target :: buf_N_ci_a(:)
      double precision, pointer             :: N_ci_a(:,:)
      double precision, allocatable, target :: buf_N_ci_b(:)
      double precision, pointer             :: N_ci_b(:,:)
      double precision, allocatable, target :: buf_N_bi_a(:)
      double precision, pointer             :: N_bi_a(:,:)
      double precision, allocatable, target :: buf_N_bi_b(:)
      double precision, pointer             :: N_bi_b(:,:)
      double precision, allocatable, target :: buf_N_ai_a(:)
      double precision, pointer             :: N_ai_a(:,:)
      double precision, allocatable, target :: buf_N_ai_b(:)
      double precision, pointer             :: N_ai_b(:,:)
c
c     Occupation numbers of the 2RDM exchange terms
c
      double precision, allocatable, target :: buf_x_aa(:)
      double precision, pointer             :: x_aa(:,:)
      double precision, allocatable, target :: buf_x_bb(:)
      double precision, pointer             :: x_bb(:,:)
      double precision, allocatable, target :: buf_dx_aa(:)
      double precision, pointer             :: dx_aa(:,:)
      double precision, allocatable, target :: buf_dx_bb(:)
      double precision, pointer             :: dx_bb(:,:)
c
c     Functions:
c
      double precision, external :: ddot
c
c     Code:
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      itask = -1
c
      E_2el     = 0.0d0
      E_2el_chf = 0.0d0
c
      nmo  = wfn1_param%nmo
      ne_a = wfn1_param%nea
      ne_b = wfn1_param%neb
c
      if (.not. bas_numcont(wfn1_inst%basis, nshell)) call errquit
     &     ('wfn1_e_g_2el: bas_numcont failed for basis',
     &      wfn1_inst%basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(wfn1_inst%basis, maxsbf)) call errquit
     &     ('wfn1_e_g_2el: bas_nbf_cn_max failed for basis',
     &      wfn1_inst%basis, BASIS_ERR)
c
      call ga_distribution(g_x_aa,0,ii_lo,ii_hi,jj_lo,jj_hi)
      ni = ii_hi-ii_lo+1
      nj = jj_hi-jj_lo+1
      call int_mem_2e4c(max2e,mem2)
c
      allocate(occ_a(1:nmo),occ_b(1:nmo))
      allocate(scr(mem2))
      allocate(buf_eri(max2e),buf_rdm2(max2e),buf_rdm2_chf(max2e))
c
      allocate(buf_rdm2_adj_a(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_adj_b(maxsbf*maxsbf*nj))
c
      allocate(buf_rdm2_d_adj_a(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_d_adj_b(maxsbf*maxsbf*nj))
c
      allocate(buf_N_ai_a(maxsbf*nmo),buf_N_ai_b(maxsbf*nmo))
      allocate(buf_N_bi_a(maxsbf*nmo),buf_N_bi_b(maxsbf*nmo))
      allocate(buf_N_ci_a(maxsbf*nmo),buf_N_ci_b(maxsbf*nmo))
      allocate(buf_N_di_a(maxsbf*nmo),buf_N_di_b(maxsbf*nmo))
c
      allocate(buf_rdm1_ac_a(maxsbf*maxsbf))
      allocate(buf_rdm1_ac_b(maxsbf*maxsbf))
      allocate(buf_rdm1_bd_a(maxsbf*maxsbf))
      allocate(buf_rdm1_bd_b(maxsbf*maxsbf))
c
      allocate(buf_rdm1_ad_a(maxsbf*maxsbf))
      allocate(buf_rdm1_ad_b(maxsbf*maxsbf))
      allocate(buf_rdm1_bc_a(maxsbf*maxsbf))
      allocate(buf_rdm1_bc_b(maxsbf*maxsbf))
c
      allocate(buf_x_aa(ni*nj), buf_x_bb(ni*nj))
      allocate(buf_dx_aa(ni*nj),buf_dx_bb(ni*nj))
c
      allocate(buf_fn_ac_a(maxsbf*maxsbf),buf_fn_ac_b(maxsbf*maxsbf))
      allocate(buf_fn_bcj_a(maxsbf*maxsbf*nj))
      allocate(buf_fn_bcj_b(maxsbf*maxsbf*nj))
      allocate(buf_fn_ic_a(nmo*maxsbf),buf_fn_ic_b(nmo*maxsbf))
      allocate(buf_fn_bj_a(maxsbf*nj),buf_fn_bj_b(maxsbf*nj))
      allocate(buf_fn_ij_a(nmo*nj),buf_fn_ij_b(nmo*nj))
c
      allocate(buf_fc_bcj_a(maxsbf*maxsbf*nj))
      allocate(buf_fc_bcj_b(maxsbf*maxsbf*nj))
c
c     Calculate occupation numbers
c     ----------------------------
c
c     Calculate alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occ_a)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occ_b)
c
c     Calculate alpha-alpha-occupation numbers
c
      call wfn1_calc_ocx(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,g_x_aa)
c
c     Calculate beta-beta-occupation numbers
c
      call wfn1_calc_ocx(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,g_x_bb)
c
c     Calculate exchange alpha-alpha-occupation number derivative
c
      call wfn1_calc_docx(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                    wfn1_wave%dcnta,g_dx_aa)
c
c     Calculate exchange beta-beta-occupation number derivative
c
      call wfn1_calc_docx(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                    wfn1_wave%dcntb,g_dx_bb)
cDEBUG
      write(*,*)"HVD: occ_a: ",occ_a
      write(*,*)"HVD: occ_b: ",occ_b
      write(*,*)"HVD: ocx:"
      call util_flush(6)
      call ga_print(g_x_aa)
      call ga_print(g_x_bb)
      call util_flush(6)
cDEBUG

c
      call ga_sync()
c
      do ld = 1, nshell
        if (.not. bas_cn2bfr(wfn1_inst%basis, ld, id_lo, id_hi))
     &      call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                   wfn1_inst%basis, BASIS_ERR)
c
        nd = id_hi-id_lo+1
        N_di_a(id_lo:id_hi,1:nmo) => buf_N_di_a
        N_di_b(id_lo:id_hi,1:nmo) => buf_N_di_b
        call ga_get(wfn1_wave%g_na,id_lo,id_hi,1,nmo,N_di_a,nd)
        call ga_get(wfn1_wave%g_nb,id_lo,id_hi,1,nmo,N_di_b,nd)
c
        do lb = 1, nshell
          if (.not. bas_cn2bfr(wfn1_inst%basis, lb, ib_lo, ib_hi))
     &        call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                     wfn1_inst%basis, BASIS_ERR)
c
          nb = ib_hi-ib_lo+1
          N_bi_a(ib_lo:ib_hi,1:nmo) => buf_N_bi_a
          N_bi_b(ib_lo:ib_hi,1:nmo) => buf_N_bi_b
          call ga_get(wfn1_wave%g_na,ib_lo,ib_hi,1,nmo,N_bi_a,nb)
          call ga_get(wfn1_wave%g_nb,ib_lo,ib_hi,1,nmo,N_bi_b,nb)
c
c         Construct the bd 1RDM blocks
c
          rdm1_bd_a(ib_lo:ib_hi,id_lo:id_hi) => buf_rdm1_bd_a
          rdm1_bd_b(ib_lo:ib_hi,id_lo:id_hi) => buf_rdm1_bd_b
          rdm1_bd_a = 0.0d0
          rdm1_bd_b = 0.0d0
          do id = id_lo, id_hi
            do ib = ib_lo, ib_hi
              do ii = 1, nmo
                rdm1_bd_a(ib,id) = rdm1_bd_a(ib,id)
     &                           + N_bi_a(ib,ii)*N_di_a(id,ii)*occ_a(ii)
                rdm1_bd_b(ib,id) = rdm1_bd_b(ib,id)
     &                           + N_bi_b(ib,ii)*N_di_b(id,ii)*occ_b(ii)
              enddo ! ii
            enddo ! ib
          enddo ! id
c
          do lc = 1, nshell 
            if (.not. bas_cn2bfr(wfn1_inst%basis, lc, ic_lo, ic_hi))
     &          call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                       wfn1_inst%basis, BASIS_ERR)
c
            nc = ic_hi-ic_lo+1
            N_ci_a(ic_lo:ic_hi,1:nmo) => buf_N_ci_a
            N_ci_b(ic_lo:ic_hi,1:nmo) => buf_N_ci_b
            call ga_get(wfn1_wave%g_na,ic_lo,ic_hi,1,nmo,N_ci_a,nc)
            call ga_get(wfn1_wave%g_nb,ic_lo,ic_hi,1,nmo,N_ci_b,nc)
c
c           Construct the bc 1RDM blocks
c
            rdm1_bc_a(ib_lo:ib_hi,ic_lo:ic_hi) => buf_rdm1_bc_a
            rdm1_bc_b(ib_lo:ib_hi,ic_lo:ic_hi) => buf_rdm1_bc_b
            rdm1_bc_a = 0.0d0
            rdm1_bc_b = 0.0d0
            do ic = ic_lo, ic_hi
              do ib = ib_lo, ib_hi
                do ii = 1, nmo
                  rdm1_bc_a(ib,ic) = rdm1_bc_a(ib,ic)
     &                             + N_bi_a(ib,ii)*N_ci_a(ic,ii)*
     &                               occ_a(ii)
                  rdm1_bc_b(ib,ic) = rdm1_bc_b(ib,ic)
     &                             + N_bi_b(ib,ii)*N_ci_b(ic,ii)*
     &                               occ_b(ii)
                enddo ! ii
              enddo ! ib
            enddo ! ic
c
            do la = 1, nshell
              itask = itask + 1
              if (mod(itask,nproc).ne.iproc) cycle
              if (.not. bas_cn2bfr(wfn1_inst%basis, la, ia_lo, ia_hi))
     &            call errquit('wfn1_e_g_2el: bas_cn2bfr',
     &                         wfn1_inst%basis, BASIS_ERR)
c
              na = ia_hi-ia_lo+1
              N_ai_a(ia_lo:ia_hi,1:nmo) => buf_N_ai_a
              N_ai_b(ia_lo:ia_hi,1:nmo) => buf_N_ai_b
              call ga_get(wfn1_wave%g_na,ia_lo,ia_hi,1,nmo,N_ai_a,na)
              call ga_get(wfn1_wave%g_nb,ia_lo,ia_hi,1,nmo,N_ai_b,na)
c
c             Construct the ac 1RDM blocks
c
              rdm1_ac_a(ia_lo:ia_hi,ic_lo:ic_hi) => buf_rdm1_ac_a
              rdm1_ac_b(ia_lo:ia_hi,ic_lo:ic_hi) => buf_rdm1_ac_b
              rdm1_ac_a = 0.0d0
              rdm1_ac_b = 0.0d0
              do ic = ic_lo, ic_hi
                do ia = ia_lo, ia_hi
                  do ii = 1, nmo
                    rdm1_ac_a(ia,ic) = rdm1_ac_a(ia,ic)
     &                               + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                                 occ_a(ii)
                    rdm1_ac_b(ia,ic) = rdm1_ac_b(ia,ic)
     &                               + N_ai_b(ia,ii)*N_ci_b(ic,ii)*
     &                                 occ_b(ii)
                  enddo ! ii
                enddo ! ia
              enddo ! ic
c
c             Construct the ad 1RDM blocks (conventional Hartree-
c             Fock exchange)
c
              rdm1_ad_a(ia_lo:ia_hi,id_lo:id_hi) => buf_rdm1_ad_a
              rdm1_ad_b(ia_lo:ia_hi,id_lo:id_hi) => buf_rdm1_ad_b
              rdm1_ad_a = 0.0d0
              rdm1_ad_b = 0.0d0
              do id = id_lo, id_hi
                do ia = ia_lo, ia_hi
                  do ii = 1, nmo
                    rdm1_ad_a(ia,id) = rdm1_ad_a(ia,id)
     &                               + N_ai_a(ia,ii)*N_di_a(id,ii)*
     &                                 occ_a(ii)
                    rdm1_ad_b(ia,id) = rdm1_ad_b(ia,id)
     &                               + N_ai_b(ia,ii)*N_di_b(id,ii)*
     &                                 occ_b(ii)
                  enddo ! ii
                enddo ! ia
              enddo ! id
c
c             Set up eri, and rdm2 buffers
c
              eri(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_eri
              rdm2(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_rdm2
              rdm2_chf(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_rdm2_chf
c
c             Calculate the 2-electron integrals
c
              call int_2e4c(wfn1_inst%basis, ld, lb,
     &                      wfn1_inst%basis, lc, la,
     &                      mem2, scr, max2e, eri)
c
              rdm2     = 0.0d0
              rdm2_chf = 0.0d0
c
c             Add the 2RDM terms that can be constructed from
c             1RDM factors.
c
              do id = id_lo, id_hi
                do ib = ib_lo, ib_hi
                  do ic = ic_lo, ic_hi
                    do ia = ia_lo, ia_hi
                      rdm2(ia,ic,ib,id) = rdm2(ia,ic,ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_a(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_a(ib,id)
                      rdm2_chf(ia,ic,ib,id) = rdm2_chf(ia,ic,ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_a(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_a(ib,id)
     &                  - rdm1_ad_a(ia,id)*rdm1_bc_a(ib,ic)
     &                  - rdm1_ad_b(ia,id)*rdm1_bc_b(ib,ic)
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
c
c             Construct contributions to the natural orbital Fock
c             matrices.
c
              fn_ac_a(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fn_ac_a
              fn_ac_b(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fn_ac_b
              fn_ac_a = 0.0d0
              fn_ac_b = 0.0d0
c
              fn_ic_a(1:nmo,ic_lo:ic_hi) => buf_fn_ic_a
              fn_ic_b(1:nmo,ic_lo:ic_hi) => buf_fn_ic_b
              fn_ic_a = 0.0d0
              fn_ic_b = 0.0d0
c
              do id = id_lo, id_hi
                do ib = ib_lo, ib_hi
                  do ic = ic_lo, ic_hi
                    do ia = ia_lo, ia_hi
                      fn_ac_a(ia,ic) = fn_ac_a(ia,ic)
     &                  + eri(ia,ic,ib,id)*
     &                    (rdm1_bd_a(ib,id)+rdm1_bd_b(ib,id))
                      fn_ac_b(ia,ic) = fn_ac_b(ia,ic)
     &                  + eri(ia,ic,ib,id)*
     &                    (rdm1_bd_a(ib,id)+rdm1_bd_b(ib,id))
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
c
              do ic = ic_lo, ic_hi
                do ia = ia_lo, ia_hi
                  do ii = 1, nmo
                    fn_ic_a(ii,ic) = fn_ic_a(ii,ic)
     &                + fn_ac_a(ia,ic)*occ_a(ii)*N_ai_a(ia,ii)
                    fn_ic_b(ii,ic) = fn_ic_b(ii,ic)
     &                + fn_ac_b(ia,ic)*occ_b(ii)*N_ai_b(ia,ii)
                  enddo ! ii
                enddo ! ia
              enddo ! ic
c
              do jj = 1, nmo, nj
                fn_ij_a(1:nmo,jj:min(nmo,jj+nj-1)) => buf_fn_ij_a
                fn_ij_b(1:nmo,jj:min(nmo,jj+nj-1)) => buf_fn_ij_b
                fn_ij_a = 0.0d0
                fn_ij_b = 0.0d0
                do kk = jj, min(nmo,jj+nj-1)
                  do ic = ic_lo, ic_hi
                    do ii = 1, nmo
                      fn_ij_a(ii,kk) = fn_ij_a(ii,kk)
     &                  + fn_ic_a(ii,ic)*N_ci_a(ic,kk)
                      fn_ij_b(ii,kk) = fn_ij_b(ii,kk)
     &                  + fn_ic_b(ii,ic)*N_ci_b(ic,kk)
                    enddo ! ii
                  enddo ! ic
                enddo ! kk
                call ga_acc(wfn1_deriv%g_dna,1,nmo,jj,min(nmo,jj+nj-1),
     &                      fn_ij_a,nmo,1.0d0)
                call ga_acc(wfn1_deriv%g_dnb,1,nmo,jj,min(nmo,jj+nj-1),
     &                      fn_ij_b,nmo,1.0d0)
              enddo ! jj
c
c             Reuse the fn_ac_* data to construct the Coulomb terms
c             to the correlation function Fock matrices.
c
              do ii = 1, nmo
                do ic = ic_lo, ic_hi
                  do ia = ia_lo, ia_hi
                    v_dc_a(ii) = v_dc_a(ii)
     &                + fn_ac_a(ia,ic)*N_ai_a(ia,ii)*N_ci_a(ic,ii)
                    v_dc_b(ii) = v_dc_b(ii)
     &                + fn_ac_b(ia,ic)*N_ai_b(ia,ii)*N_ci_b(ic,ii)
                  enddo ! ia
                enddo ! ic
              enddo ! ii
c
c             Add the 2RDM terms that come from the exchange terms
c             of [5].
c
              do jproc = 0, nproc-1
                call ga_distribution(g_x_aa,jproc,ii_lo,ii_hi,
     &                                            jj_lo,jj_hi)
                ni = ii_hi-ii_lo+1
c
                x_aa(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_x_aa
                x_bb(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_x_bb
                dx_aa(ii_lo:ii_hi,jj_lo:jj_hi) => buf_dx_aa
                dx_bb(ii_lo:ii_hi,jj_lo:jj_hi) => buf_dx_bb
c
                call ga_get(g_x_aa,ii_lo,ii_hi,jj_lo,jj_hi,x_aa,ni)
                call ga_get(g_x_bb,ii_lo,ii_hi,jj_lo,jj_hi,x_bb,ni)
                call ga_get(g_dx_aa,ii_lo,ii_hi,jj_lo,jj_hi,dx_aa,ni)
                call ga_get(g_dx_bb,ii_lo,ii_hi,jj_lo,jj_hi,dx_bb,ni)
c
                rdm2_adj_a(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_a
                rdm2_adj_b(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_b
c
                rdm2_adj_a = 0.0d0
                rdm2_adj_b = 0.0d0
c
                rdm2_d_adj_a(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_adj_a
                rdm2_d_adj_b(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_adj_b
c
                rdm2_d_adj_a = 0.0d0
                rdm2_d_adj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ii = ii_lo, ii_hi
                    do id = id_lo, id_hi
                      do ia = ia_lo, ia_hi
                        rdm2_adj_a(ia,id,jj) = rdm2_adj_a(ia,id,jj)
     &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*x_aa(ii,jj)
                        rdm2_adj_b(ia,id,jj) = rdm2_adj_b(ia,id,jj)
     &                    + N_ai_b(ia,ii)*N_di_b(id,ii)*x_bb(ii,jj)
c
                        rdm2_d_adj_a(ia,id,jj) = rdm2_d_adj_a(ia,id,jj)
     &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*dx_aa(ii,jj)
                        rdm2_d_adj_b(ia,id,jj) = rdm2_d_adj_b(ia,id,jj)
     &                    + N_ai_b(ia,ii)*N_di_b(id,ii)*dx_bb(ii,jj)
                      enddo ! ia
                    enddo ! id
                  enddo ! ii
                enddo ! jj
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
c
c                         Note: sign is absorbed in x_aa and x_bb
c
                          rdm2(ia,ic,ib,id) = rdm2(ia,ic,ib,id)
     &                      + rdm2_adj_a(ia,id,jj)*N_bi_a(ib,jj)*
     &                        N_ci_a(ic,jj)
     &                      + rdm2_adj_b(ia,id,jj)*N_bi_b(ib,jj)*
     &                        N_ci_b(ic,jj)
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
c               Calculate the natural orbital exchange Fock matrix
c               contributions
c
                fn_bcj_a(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fn_bcj_a
                fn_bcj_b(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fn_bcj_b
                fn_bcj_a = 0.0d0
                fn_bcj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
                          fn_bcj_a(ib,ic,jj) = fn_bcj_a(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_adj_a(ia,id,jj)
                          fn_bcj_b(ib,ic,jj) = fn_bcj_b(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_adj_b(ia,id,jj)
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
                fn_bj_a(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fn_bj_a
                fn_bj_b(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fn_bj_b
                fn_bj_a = 0.0d0
                fn_bj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ic = ic_lo, ic_hi
                    do ib = ib_lo, ib_hi
                      fn_bj_a(ib,jj) = fn_bj_a(ib,jj)
     &                  + fn_bcj_a(ib,ic,jj)*N_ci_a(ic,jj)
                      fn_bj_b(ib,jj) = fn_bj_b(ib,jj)
     &                  + fn_bcj_b(ib,ic,jj)*N_ci_b(ic,jj)
                    enddo ! ib
                  enddo ! ic
                enddo ! jj
c
                fn_ij_a(1:nmo,jj_lo:jj_hi) => buf_fn_ij_a
                fn_ij_b(1:nmo,jj_lo:jj_hi) => buf_fn_ij_b
                fn_ij_a = 0.0d0
                fn_ij_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ii = 1, nmo
                    do ib = ib_lo, ib_hi
                      fn_ij_a(ii,jj) = fn_ij_a(ii,jj)
     &                  + fn_bj_a(ib,jj)*N_bi_a(ib,ii)
                      fn_ij_b(ii,jj) = fn_ij_b(ii,jj)
     &                  + fn_bj_b(ib,jj)*N_bi_b(ib,ii)
                    enddo ! ib
                  enddo ! ii
                enddo ! jj
c
                call ga_acc(wfn1_deriv%g_dna,1,nmo,jj_lo,jj_hi,
     &                      fn_ij_a,nmo,1.0d0)
                call ga_acc(wfn1_deriv%g_dnb,1,nmo,jj_lo,jj_hi,
     &                      fn_ij_b,nmo,1.0d0)
c
c               Calculate the correlation function exchange Fock matrix
c               contributions
c
                fc_bcj_a(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fc_bcj_a
                fc_bcj_b(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fc_bcj_b
                fc_bcj_a = 0.0d0
                fc_bcj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
                          fc_bcj_a(ib,ic,jj) = fc_bcj_a(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_d_adj_a(ia,id,jj)
                          fc_bcj_b(ib,ic,jj) = fc_bcj_b(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_d_adj_b(ia,id,jj)
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
                do jj = jj_lo, jj_hi
                  do ic = ic_lo, ic_hi
                    do ib = ib_lo, ib_hi
                      v_dc_a(jj) = v_dc_a(jj)
     &                  + fc_bcj_a(ib,ic,jj)*N_bi_a(ib,jj)*N_ci_a(ic,jj)
                      v_dc_b(jj) = v_dc_b(jj)
     &                  + fc_bcj_b(ib,ic,jj)*N_bi_b(ib,jj)*N_ci_b(ic,jj)
                    enddo ! ib
                  enddo ! ic
                enddo ! jj
c
              enddo ! jproc
c
c             Add the 2-electron energy terms
c
              E_2el     = E_2el     + ddot(na*nb*nc*nd,eri,1,rdm2,1)
              E_2el_chf = E_2el_chf + ddot(na*nb*nc*nd,eri,1,rdm2_chf,1)
c
            enddo ! la
          enddo ! lc
        enddo ! lb
      enddo ! ld
c
      E_2el     = 0.5d0*E_2el
      E_2el_chf = 0.5d0*E_2el_chf
c
      deallocate(buf_fc_bcj_a,buf_fc_bcj_b)
c
      deallocate(buf_fn_ic_a,buf_fn_ic_b)
      deallocate(buf_fn_bcj_a,buf_fn_bcj_b)
      deallocate(buf_fn_ac_a,buf_fn_ac_b)
c
      deallocate(buf_dx_aa,buf_dx_bb)
      deallocate(buf_x_aa, buf_x_bb)
c
      deallocate(buf_rdm1_bc_a,buf_rdm1_bc_b)
      deallocate(buf_rdm1_ad_a,buf_rdm1_ad_b)
c
      deallocate(buf_rdm1_bd_a,buf_rdm1_bd_b)
      deallocate(buf_rdm1_ac_a,buf_rdm1_ac_b)
c
      deallocate(buf_N_di_a,buf_N_di_b)
      deallocate(buf_N_ci_a,buf_N_ci_b)
      deallocate(buf_N_bi_a,buf_N_bi_b)
      deallocate(buf_N_ai_a,buf_N_ai_b)
c
      deallocate(buf_rdm2_d_adj_a,buf_rdm2_d_adj_b)
c
      deallocate(buf_rdm2_adj_a,buf_rdm2_adj_b)
c
      deallocate(buf_eri,buf_rdm2,buf_rdm2_chf)
      deallocate(scr)
      deallocate(occ_a,occ_b)
c
      end
!>
!> @}
c
c-----------------------------------------------------------------------
