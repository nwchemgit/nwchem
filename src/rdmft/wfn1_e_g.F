c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                    wfn1_deriv,wfn1_energy,g_bfmo,g_momo,
     &                    g_momo2,g_momo3,g_momo4,v_dca,v_dcb)
      implicit none
C>
C> \brief Calculate the energy and the derivative of the energy wrt
C> the wavefunction parameters
C>
C> This routine calculates the various energy contributions as well
C> as various Fock matrices that will provide search directions for
C> the energy minimization. Thus the different Fock matrices will be
C> calculated as:
C> - The Fock matrix for the natural orbitals will be returned in the
C>   natural orbital basis (nmo x nmo)
C> - The Fock matrix for the correlation functions will be returned in
C>   the correlation function basis (nmo x nmo)
C>
#include "errquit.fh"
#include "global.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_en.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
      type(wfn1_drv),     intent(inout) :: wfn1_deriv
      type(wfn1_en),      intent(out)   :: wfn1_energy
c
      integer, intent(in) :: g_bfmo  ! nbf x nmo scratch GA
      integer, intent(in) :: g_momo  ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo2 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo3 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo4 ! nmo x nmo scratch GA
c
      double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
      double precision, allocatable :: buf(:,:)
      double precision, allocatable :: bufji(:,:)
c
      integer ilo, ihi ! range limits
      integer jlo, jhi ! range limits
      integer ii       ! counter
      integer jj       ! counter
      integer nmo      ! local copy of the number of MOs
c
      integer nproc
      integer iproc
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
c     Zero all relevant arrays so we can just add contributions
c
c     call ga_zero(g_bfmo)
c     call ga_zero(g_momo)
      call ga_zero(wfn1_deriv%g_dna)
      call ga_zero(wfn1_deriv%g_dnb)
      call ga_zero(wfn1_deriv%g_dca)
      call ga_zero(wfn1_deriv%g_dcb)
c
      wfn1_energy%E_tot = 0.0d0
      wfn1_energy%E_1el = 0.0d0
      wfn1_energy%E_2el = 0.0d0
      wfn1_energy%E_nuc = wfn1_param%Enuc
      wfn1_energy%E_ent = 0.0d0
c
c     Arrange memory for the diagonal correlation function Fock
c     matrices in natural orbital basis
c
c     allocate(v_dca(1:wfn1_param%nmo))
c     allocate(v_dcb(1:wfn1_param%nmo))
c
      v_dca = 0.0d0
      v_dcb = 0.0d0
c
c     Evaluate the various Fock matrix terms
c     --------------------------------------
c
c     The 1-electron terms
c
      call wfn1_e_g_1el(wfn1_param,wfn1_mat,wfn1_wave,wfn1_deriv,
     &                  v_dca,v_dcb,g_bfmo,g_momo,wfn1_energy%E_1el)
c
c     The 2-electron terms
c
      if (wfn1_param%corr_expr == "wfn1") then
        call wfn1_e_g_2el(wfn1_param,wfn1_inst,wfn1_wave,wfn1_deriv,
     &                    v_dca,v_dcb,g_momo,g_momo2,g_momo3,g_momo4,
     &                    wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
      else if (wfn1_param%corr_expr == "wfn1s") then
        call wfn1_e_g_2el_s(wfn1_param,wfn1_inst,wfn1_wave,wfn1_deriv,
     &                      v_dca,v_dcb,g_momo,g_momo2,g_momo3,g_momo4,
     &                      wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
      else
        call errquit("invalid correlation energy option",0,UERR)
      endif
c
c     The entropy term
c
      call wfn1_e_g_ent(wfn1_param,wfn1_wave,v_dca,v_dcb,
     &                  wfn1_energy%E_ent)
c
c     Transform the correlation function Fock matrices from natural
c     orbital to correlation function basis
c
      call ga_copy(wfn1_wave%g_ca,g_momo)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      call ga_get(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          buf(ii,jj) = buf(ii,jj)*v_dca(ii)
        enddo
      enddo
      nmo = wfn1_param%nmo
      call ga_put(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_dgemm('T','N',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,g_momo,
     &              0.0d0,wfn1_deriv%g_dca)
c
      call ga_copy(wfn1_wave%g_cb,g_momo)
      call ga_get(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          buf(ii,jj) = buf(ii,jj)*v_dcb(ii)
        enddo
      enddo
      nmo = wfn1_param%nmo
      call ga_put(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_dgemm('T','N',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,g_momo,
     &              0.0d0,wfn1_deriv%g_dcb)
      deallocate(buf)
c
c     Clean up memory
c
c     deallocate(v_dcb,v_dca)
c
c     Make derivative matrices skew-symmetric
c
      call ga_distribution(wfn1_deriv%g_dna,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dna,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dna,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = -0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            buf(ii,jj) =  0.0d0
          else if (ii.gt.jj) then
            buf(ii,jj) = +0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dna,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      call ga_distribution(wfn1_deriv%g_dnb,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dnb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dnb,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = -0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            buf(ii,jj) =  0.0d0
          else if (ii.gt.jj) then
            buf(ii,jj) = +0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dnb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      call ga_distribution(wfn1_deriv%g_dca,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dca,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = -0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            buf(ii,jj) =  0.0d0
          else if (ii.gt.jj) then
            buf(ii,jj) = +0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      call ga_distribution(wfn1_deriv%g_dcb,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dcb,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = -0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            buf(ii,jj) =  0.0d0
          else if (ii.gt.jj) then
            buf(ii,jj) = +0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      wfn1_energy%E_tot = wfn1_energy%E_1el + wfn1_energy%E_2el
     &                  + wfn1_energy%E_nuc + wfn1_energy%E_ent
c
      call ga_sync()
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_1el(wfn1_param,wfn1_mat,wfn1_wave,wfn1_deriv,
     &                        v_dca,v_dcb,g_bfmo,g_momo,E_1el)
      implicit none
C>
C> \brief Calculate the 1-electron contributions to the energy,
C> the natural orbital and the correlation function Fock matrices
C>
C> The 1-electron term to the Fock matrices will be calculated as:
C> - The natural orbital term will calculated in the natural orbital
C>   basis
C> - The correlation function Fock matrices will be calculated in the
C>   natural orbital basis because the matrix is diagonal in that 
C>   basis
C> Regarding the correlation function Fock matrices these will be
C> calculated in the natural orbital basis in this routine. Once the
C> full diagonal elements have calculated (including the 2-electron
C> terms) then the matrix will be transformed to the correlation
C> function basis. Doing this transformation prematurely just
C> generates unnecesary work.
C> The algorithm is rather simple:
C> 1. Transform the 1-electron energy matrix from the basis functions
C>    to the natural orbital basis
C> 2. Select the diagonal of the result of 1. 
C> 3. Multiply the elements from 2. with the occupation numbers and sum
C>    to get the 1-electron energy.
C> 4. Scale the result from 1. with symmetrized occupation number 
C>    factors to get the natural orbital Fock matrix 1-electron 
C>    contribution term.
C> 5. Add the data from 2. to the diagonal of the correlation function
C>    Fock matrix in the natural orbital basis.
C> 
#include "wfn1_matcnst.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
      type(wfn1_drv),     intent(inout) :: wfn1_deriv
c
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer, intent(inout) :: g_bfmo ! nbf x nmo scratch GA
      integer, intent(inout) :: g_momo ! nmo x nmo scratch GA
c
      double precision, intent(out) :: E_1el    ! the 1-electron energy
c
      double precision, allocatable :: occa(:)  ! the alpha-occupation 
                                                ! numbers
      double precision, allocatable :: occb(:)  ! the beta-occupation 
                                                ! numbers
      double precision, allocatable :: e1a(:)   ! the alpha-1-electron
                                                ! energies
      double precision, allocatable :: e1b(:)   ! the beta-1-electron
                                                ! energies
      double precision, allocatable :: buf(:,:) ! a local buffer
c
c     local
c
      integer iclo, ichi, jclo, jchi ! related to the correlation 
                                     ! functions
      integer imlo, imhi, jmlo, jmhi ! related to the mo x mo 
                                     ! matrix
      integer ii, jj                 ! counters
      integer nmo                    ! number of MOs
      integer nea                    ! number of alpha electrons
      integer neb                    ! number of beta  electrons
      integer iproc, nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      nmo   = wfn1_param%nmo
      nea   = wfn1_param%nea
      neb   = wfn1_param%neb
c
      allocate(occa(1:nmo),occb(1:nmo))
      allocate(e1a(1:nmo),e1b(1:nmo))
c
      occa = 0.0d0
      occb = 0.0d0
      e1a  = 0.0d0
      e1b  = 0.0d0
c
c     Do the ALPHA spin channel
c     -------------------------
c
c     Calculate the alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
c
c     Do step 1.
c
      call ga_dgemm('n','n',wfn1_param%nbf,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_mat%g_tv,wfn1_wave%g_na,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',wfn1_param%nmo,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_wave%g_na,g_bfmo,
     &              0.0d0,g_momo)
c
c     Do step 2.
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1a(ii),1)
      enddo
      call ga_dgop(MT_DBL,e1a,wfn1_param%nmo,"+")
c
c     Do step 3.
c
      do ii = 1, wfn1_param%nmo
        E_1el = E_1el + occa(ii)*e1a(ii)
      enddo
c
c     Do step 4.
c
      allocate(buf(imlo:imhi,jmlo:jmhi))
      call ga_get(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
      do jj = jmlo, jmhi
        do ii = imlo, imhi
          buf(ii,jj) = buf(ii,jj)*occa(ii)
        enddo
      enddo
      call ga_put(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
      deallocate(buf)
      call ga_add(1.0d0,wfn1_deriv%g_dna,
     &            1.0d0,g_momo,wfn1_deriv%g_dna)
c
c     Do step 5.
c
      v_dca = v_dca + e1a
c
c     Do the BETA spin channel
c     ------------------------
c
c     Calculate the beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
c
c     Do step 1.
c
      call ga_dgemm('n','n',wfn1_param%nbf,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_mat%g_tv,wfn1_wave%g_nb,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',wfn1_param%nmo,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_wave%g_nb,g_bfmo,
     &              0.0d0,g_momo)
c
c     Do step 2.
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1b(ii),1)
      enddo
      call ga_dgop(MT_DBL,e1b,wfn1_param%nmo,"+")
c
c     Do step 3.
c
      do ii = 1, wfn1_param%nmo
        E_1el = E_1el + occb(ii)*e1b(ii)
      enddo
c
c     Do step 4.
c
      allocate(buf(imlo:imhi,jmlo:jmhi))
      call ga_get(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
      do jj = jmlo, jmhi
        do ii = imlo, imhi
          buf(ii,jj) = buf(ii,jj)*occb(ii)
        enddo
      enddo
      call ga_put(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
      deallocate(buf)
      call ga_add(1.0d0,wfn1_deriv%g_dnb,
     &            1.0d0,g_momo,wfn1_deriv%g_dnb)
c
c     Do step 5.
c
      v_dcb = v_dcb + e1b
c
c     Done the work, now tidy up
c
      deallocate(e1a,e1b)
      deallocate(occa,occb)
      call ga_sync()
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_ent(wfn1_param,wfn1_wave,v_dca,v_dcb,E_ent)
      implicit none
C>
C> \brief Calculate the entropy derived energy contribution
C>
C> The entropy derived contribution is given by
C> \f{eqnarray}{
C>   E_{\mathrm{ent}} &=& -T_{\mathrm{bath}}S \\\\
C>   S &=& -\sum_i \rho_i\ln(\rho_i) + (1-\rho_i)ln(1-\rho_i)
C> \f}
C> where \f$T_{\mathrm{bath}}\f$ is the bath temperature, \f$S\f$ is
C> the entropy, and \f$\rho_i\f$ is the occupation number of the
C> \f$i\f$-th natural orbital.
C>
C> The derivative of the entropy with respect to the occupation numbers
C> is given by
C> \f{eqnarray}{
C>   \frac{\partial S}{\partial\rho_i}
C>   &=& -\ln(\rho_i)+ln(1-\rho_i)
C> \f}
C> Obviously if \f$\rho_i = 1\f$ or \f$\rho_i = 0\f$ both the entropy
C> and its derivative are undefined. Hence we need to introduce a 
C> cut-off \f$\rho_c\f$. 
C>
C> If \f$\rho_i \lt \rho_c\f$ then 
C> \f{eqnarray}{
C>   S(\rho_i) &=& 0 \\\\
C>   \frac{\partial S}{\partial\rho_i} 
C>   &=& -\ln(\rho_c)
C> \f}
C> If \f$\rho_i \gt 1-\rho_c\f$ then 
C> \f{eqnarray}{
C>   S(\rho_i) &=& 0 \\\\
C>   \frac{\partial S}{\partial\rho_i} 
C>   &=& \ln(1-\rho_c)
C> \f}
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
c
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
      double precision, intent(inout) :: E_ent
c
      double precision, allocatable :: occa(:)  ! the alpha-occupation
                                                ! numbers
      double precision, allocatable :: occb(:)  ! the beta-occupation
                                                ! numbers
      double precision, allocatable :: e1a(:)   ! the alpha-1-electron
                                                ! energies
      double precision, allocatable :: e1b(:)   ! the beta-1-electron
                                                ! energies
      double precision, allocatable :: buf(:,:) ! a local buffer
c
c     local
c
      double precision :: dr, d1r    ! entropy derivatives
      double precision :: rho_c      ! cut-off for occupation numbers
c
      integer iclo, ichi, jclo, jchi ! related to the correlation
                                     ! functions
      integer ii, jj                 ! counters
      integer iproc                  ! the current processor rank
      integer nproc                  ! the total number of processor
                                     ! ranks
c
      if (wfn1_param%Tbath.le.0.0d0) return
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      allocate(occa(1:wfn1_param%nmo),occb(1:wfn1_param%nmo))
      allocate(e1a(1:wfn1_param%nmo),e1b(1:wfn1_param%nmo))
c
      rho_c = 1.0d-10
      occa  = 0.0d0
      occb  = 0.0d0
      e1a   = 0.0d0
      e1b   = 0.0d0
c
c     Do the ALPHA spin channel
c     -------------------------
c
c     Calculate the alpha-occupation numbers
c
      call ga_distribution(wfn1_wave%g_ca,iproc,iclo,ichi,jclo,jchi)
      allocate(buf(iclo:ichi,jclo:jchi))
      call ga_get(wfn1_wave%g_ca,iclo,ichi,jclo,jchi,buf,ichi-iclo+1)
      do jj = 1, wfn1_param%nea
        if(wfn1_wave%icnfa(jj).ge.jclo.and.wfn1_wave%icnfa(jj).le.jchi)
     &  then
          do ii = iclo, ichi
            occa(ii) = occa(ii)
     &               + buf(ii,wfn1_wave%icnfa(jj))
     &               * buf(ii,wfn1_wave%icnfa(jj))
     &               * wfn1_wave%dcnta(jj)
          enddo
        endif
      enddo
      deallocate(buf)
      call ga_dgop(MT_DBL,occa,wfn1_param%nmo,"+")
c
c     Calculate the entropy
c
      do ii = 1, wfn1_param%nmo
        if      (occa(ii).lt.rho_c) then
          dr      = log(rho_c)
          e1a(ii) = e1a(ii) + wfn1_param%Tbath*dr
        else if (occa(ii).gt.(1.0d0-rho_c)) then
          d1r     = log(rho_c)
          e1a(ii) = e1a(ii) + wfn1_param%Tbath*(-d1r)
        else
          dr      = log(occa(ii))
          d1r     = log(1.0d0-occa(ii))
          E_ent   = E_ent + wfn1_param%Tbath
     &            * (occa(ii)*dr+(1.0d0-occa(ii))*d1r)
          e1a(ii) = e1a(ii) + wfn1_param%Tbath*(dr-d1r)
        endif
      enddo
c
c     Do the BETA spin channel
c     ------------------------
c
c     Calculate the beta-occupation numbers
c
      call ga_distribution(wfn1_wave%g_cb,iproc,iclo,ichi,jclo,jchi)
      allocate(buf(iclo:ichi,jclo:jchi))
      call ga_get(wfn1_wave%g_cb,iclo,ichi,jclo,jchi,buf,ichi-iclo+1)
      do jj = 1, wfn1_param%neb
        if(wfn1_wave%icnfb(jj).ge.jclo.and.wfn1_wave%icnfb(jj).le.jchi)
     &  then
          do ii = iclo, ichi
            occb(ii) = occb(ii)
     &               + buf(ii,wfn1_wave%icnfb(jj))
     &               * buf(ii,wfn1_wave%icnfb(jj))
     &               * wfn1_wave%dcntb(jj)
          enddo
        endif
      enddo
      deallocate(buf)
      call ga_dgop(MT_DBL,occb,wfn1_param%nmo,"+")
c
c     Calculate the entropy
c
      do ii = 1, wfn1_param%nmo
        if      (occb(ii).lt.rho_c) then
          dr      = log(rho_c)
          e1b(ii) = e1b(ii) + wfn1_param%Tbath*dr
        else if (occb(ii).gt.(1.0d0-rho_c)) then
          d1r     = log(rho_c)
          e1b(ii) = e1b(ii) + wfn1_param%Tbath*(-d1r)
        else
          dr      = log(occb(ii))
          d1r     = log(1.0d0-occb(ii))
          E_ent   = E_ent + wfn1_param%Tbath
     &            * (occb(ii)*dr+(1.0d0-occb(ii))*d1r)
          e1b(ii) = e1b(ii) + wfn1_param%Tbath*(dr-d1r)
        endif
      enddo
c
c     Store the results
c     -----------------
c
      v_dca = v_dca + e1a
      v_dcb = v_dcb + e1b
c
c     Clean up memory
c
      deallocate(occa,occb,e1a,e1b)
c
      end
c
c-----------------------------------------------------------------------
c
  
