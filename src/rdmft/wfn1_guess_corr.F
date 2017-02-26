c-----------------------------------------------------------------------
c
      subroutine wfn1_guess_corr(wfn1_param,wfn1_wave)
      implicit none
C>
C> \brief Guess the correlation functions 
C>
C> It is essential that the correlation functions have no zero 
C> coefficients as those coefficients will essentially be trapped
C> at that value. 
C>
C> Hence we start by setting for Nelectron correlation functions
C> all Nelectron lowest coefficients to 1. The remaining coefficients
C> will be set to a small value (e.g. 1.0d-3). For the correlation
C> functions beyond Nelectron we do the opposite. Then using basic
C> trigoniometric functions the sign of some of the major coefficient
C> is flipped to ensure a linearly independent set of vectors is 
C> created. Subsequently the set of correlation functions is ortho-
C> normalized to obtain the final guess set.
C>
C> In the ortho-normalization it is important that the order of the 
C> vectors is preserved so that we can sensibly choose which orbitals
C> to occupy. This means that diagonalization based approaches are
C> essentially out of the question. Modified Gramm-Schmidt on the
C> other hand would work.
C>
C> see also wfn1_nxt_init_occ
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "global.fh"
      type(wfn1_prm), intent(in)    :: wfn1_param
      type(wfn1_wfn), intent(inout) :: wfn1_wave
c
      integer nproc   ! the number of ranks
      integer iproc   ! the number of the current rank
      integer ii      ! counter
      integer jj      ! counter
      integer ilo,ihi ! loop limits
      integer jlo,jhi ! loop limits
c
      integer nea     ! local number of alpha electrons
      integer neb     ! local number of alpha electrons
      integer nmo     ! local number of molecular orbitals
      integer ncorea  ! local number of alpha core orbitals
      integer ncoreb  ! local number of beta core orbitals
      integer nvirta  ! local number of alpha virtual orbitals
      integer nvirtb  ! local number of beta virtual orbitals
c
      double precision, allocatable :: buf(:,:)
      double precision, allocatable :: eig(:)
      double precision :: pi
      double precision :: delta
      double precision :: x
      double precision :: dnrm
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      pi = acos(-1.0d0)
      nea = wfn1_param%nea
      neb = wfn1_param%neb
      nmo = wfn1_param%nmo
      ncorea = wfn1_param%ncorea
      ncoreb = wfn1_param%ncoreb
      nvirta = wfn1_param%nvirta
      nvirtb = wfn1_param%nvirtb
      delta = 0.0d0
c
c     Setup the alpha correlation functions
c
      call ga_distribution(wfn1_wave%g_ca,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
c
c     Put initial values in
c
      buf = 0.0d0
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (jj.le.ncorea) then
            if (ii.eq.jj) then
              buf(ii,jj) = 1.0d0
            endif
          else if (jj.le.nea) then
            if (ii.le.ncorea) then
              buf(ii,jj) = 0.0d0
            else if (ii.le.nea) then
              buf(ii,jj) = 1.0d0
            else if (ii.le.nmo-nvirta) then
              buf(ii,jj) = 1.0d-3
            else
              buf(ii,jj) = 0.0d0
            endif
          else if (jj.le.nmo-nvirta) then
            if (ii.le.ncorea) then
              buf(ii,jj) = 0.0d0
            else if (ii.le.nea) then
              buf(ii,jj) = 1.0d-3
            else if (ii.le.nmo-nvirta) then
              buf(ii,jj) = 1.0d0
            else
              buf(ii,jj) = 0.0d0
            endif
          else
            if (ii.eq.jj) then
              buf(ii,jj) = 1.0d0
            endif
          endif
        enddo
      enddo
c
c     Flip signs to make linearly independent functions
c
      if (nea.gt.0) then
        delta = pi/nea
      endif
      do jj = jlo, jhi
        if (jj.le.nea) then
          do ii = ilo, ihi
            if (ii.le.nea) then
              x = (ii-0.5d0)*jj*delta
              if (sin(x).lt.0.0d0) then
                buf(ii,jj) = -buf(ii,jj)
              endif
            endif
          enddo
        endif
      enddo
      if (nmo-nea.gt.0) then
        delta = pi/(nmo-nea)
      endif
      do jj = jlo, jhi
        if (jj.gt.nea) then
          do ii = ilo, ihi
            if (ii.gt.nea) then
              x = (ii-nea-0.5d0)*(jj-nea)*delta
              if (sin(x).lt.0.0d0) then
                buf(ii,jj) = -buf(ii,jj)
              endif
            endif
          enddo
        endif
      enddo
c
c     Normalize vectors
c
      do jj = jlo, jhi
        if (jj.le.ncorea) then
          dnrm = 1.0d0 ! already normalized
        else if (jj.le.nea) then
          dnrm = max(0,nea-ncorea)*(1.0d0**2)
     &         + max(0,nmo-nea-nvirta)*(1.0d-3**2)
        else if (jj.le.nmo-nvirta) then
          dnrm = max(0,nea-ncorea)*(1.0d-3**2)
     &         + max(0,nmo-nea-nvirta)*(1.0d0**2)
        else
          dnrm = 1.0d0 ! already normalized
        endif
        dnrm = 1.0d0/sqrt(dnrm)
        do ii = ilo, ihi
          buf(ii,jj) = dnrm*buf(ii,jj)
        enddo
      enddo
      call ga_put(wfn1_wave%g_ca,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(buf)
c
c     Orthonormalize the vectors
c
      call ga_orthog(wfn1_wave%g_ca,ii,.false.)
c
c     Setup the beta correlation functions
c
      call ga_distribution(wfn1_wave%g_cb,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
c
c     Put initial values in
c
      buf = 0.0d0
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (jj.le.ncoreb) then
            if (ii.eq.jj) then
              buf(ii,jj) = 1.0d0
            endif
          else if (jj.le.neb) then
            if (ii.le.ncoreb) then
              buf(ii,jj) = 0.0d0
            else if (ii.le.neb) then
              buf(ii,jj) = 1.0d0
            else if (ii.le.nmo-nvirtb) then
              buf(ii,jj) = 1.0d-3
            else
              buf(ii,jj) = 0.0d0
            endif
          else if (jj.le.nmo-nvirtb) then
            if (ii.le.ncoreb) then
              buf(ii,jj) = 0.0d0
            else if (ii.le.neb) then
              buf(ii,jj) = 1.0d-3
            else if (ii.le.nmo-nvirtb) then
              buf(ii,jj) = 1.0d0
            else
              buf(ii,jj) = 0.0d0
            endif
          else
            if (ii.eq.jj) then
              buf(ii,jj) = 1.0d0
            endif
          endif
        enddo
      enddo
c
c     Flip signs to make linearly independent functions
c
      if (neb.gt.0) then
        delta = pi/neb
      endif
      do jj = jlo, jhi
        if (jj.le.neb) then
          do ii = ilo, ihi
            if (ii.le.neb) then
              x = (ii-0.5d0)*jj*delta
              if (sin(x).lt.0.0d0) then
                buf(ii,jj) = -buf(ii,jj)
              endif
            endif
          enddo
        endif
      enddo
      if (nmo-neb.gt.0) then
        delta = pi/(nmo-neb)
      endif
      do jj = jlo, jhi
        if (jj.gt.neb) then
          do ii = ilo, ihi
            if (ii.gt.neb) then
              x = (ii-neb-0.5d0)*(jj-neb)*delta
              if (sin(x).lt.0.0d0) then
                buf(ii,jj) = -buf(ii,jj)
              endif
            endif
          enddo
        endif
      enddo
c
c     Normalize vectors
c
      do jj = jlo, jhi
        if (jj.le.ncoreb) then
          dnrm = 1.0d0 ! already normalized
        else if (jj.le.neb) then
          dnrm = max(0,neb-ncoreb)*(1.0d0**2)
     &         + max(0,nmo-neb-nvirtb)*(1.0d-3**2)
        else if (jj.le.nmo-nvirtb) then
          dnrm = max(0,neb-ncoreb)*(1.0d-3**2)
     &         + max(0,nmo-neb-nvirtb)*(1.0d0**2)
        else
          dnrm = 1.0d0 ! already normalized
        endif
        dnrm = 1.0d0/sqrt(dnrm)
        do ii = ilo, ihi
          buf(ii,jj) = dnrm*buf(ii,jj)
        enddo
      enddo
      call ga_put(wfn1_wave%g_cb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(buf)
c
c     Orthonormalize the vectors
c
      call ga_orthog(wfn1_wave%g_cb,ii,.false.)
c
      end
c
c-----------------------------------------------------------------------
