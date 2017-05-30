c-----------------------------------------------------------------------
c
      function wfn1_trace(wfn1_drva,wfn1_drvb)
      implicit none
C>
C> \brief Compute the trace of the gradients
C>
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "wfn1_drv.fh"
#include "wfn1_error.fh"
c
      type(wfn1_drv), intent(in) :: wfn1_drva
      type(wfn1_drv), intent(in) :: wfn1_drvb
      type(wfn1_error)           :: wfn1_trace
c
      integer ilo, ihi, jlo, jhi
      integer ii, jj
      integer itypa, nrowa, ncola
      integer itypb, nrowb, ncolb
c
      integer iproc
      integer nproc
c
      type(wfn1_error) :: error
c
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      call ga_inquire(wfn1_drva%g_dna,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dna,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dna",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dna",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dna",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dnb,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dnb,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dnb",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dnb",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dnb",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dca,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dca,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dca",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dca",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dca",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dcb,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dcb,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dcb",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dcb",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dcb",0,UERR)
c
      call ga_distribution(wfn1_drva%g_dna,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dna,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dna,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          error%errna = error%errna + bufa(ii,jj)*bufb(ii,jj)
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dnb,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dnb,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dnb,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          error%errnb = error%errnb + bufa(ii,jj)*bufb(ii,jj)
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dca,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dca,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dca,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          error%errca = error%errca + bufa(ii,jj)*bufb(ii,jj)
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dcb,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dcb,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dcb,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          error%errcb = error%errcb + bufa(ii,jj)*bufb(ii,jj)
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      error%total = error%errna + error%errnb 
     &            + error%errca + error%errcb
c
      call ga_dgop(MT_DBL,error,5,"+")
      wfn1_trace = error
c
      return
      end function wfn1_trace
c
c-----------------------------------------------------------------------
c
      function wfn1_etrace(wfn1_param,wfn1_wave,wfn1_drva,wfn1_drvb)
      implicit none
C>
C> \brief Compute the trace of the energy gradients
C>
C> This version computes the change of the energy while taking into
C> account the constrains on the energy optimization. I.e. the fact that
C> the energy does not change if the density matrix does not change.
C>
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_error.fh"
c
      type(wfn1_prm), intent(in) :: wfn1_param
      type(wfn1_wfn), intent(in) :: wfn1_wave
      type(wfn1_drv), intent(in) :: wfn1_drva
      type(wfn1_drv), intent(in) :: wfn1_drvb
      type(wfn1_error)           :: wfn1_etrace
c
      integer ilo, ihi, jlo, jhi
      integer ii, jj
      integer itypa, nrowa, ncola
      integer itypb, nrowb, ncolb
      integer nmo
      integer nea, neb
c
      integer iproc
      integer nproc
c
      type(wfn1_error) :: error
c
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
c
      double precision, allocatable :: occ_na(:)
      double precision, allocatable :: occ_nb(:)
      double precision, allocatable :: occ_ca(:)
      double precision, allocatable :: occ_cb(:)
c
      double precision, parameter :: tol = 1.0d-5
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      nmo = wfn1_param%nmo
      nea = wfn1_param%nea
      neb = wfn1_param%neb
c
      call ga_inquire(wfn1_drva%g_dna,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dna,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dna",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dna",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dna",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dnb,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dnb,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dnb",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dnb",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dnb",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dca,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dca,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dca",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dca",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dca",0,UERR)
c
      call ga_inquire(wfn1_drva%g_dcb,itypa,nrowa,ncola)
      call ga_inquire(wfn1_drvb%g_dcb,itypb,nrowb,ncolb)
      if (itypa.ne.itypb) call errquit("wrong type g_dcb",0,UERR)
      if (nrowa.ne.nrowb) call errquit("wrong nrow g_dcb",0,UERR)
      if (ncola.ne.ncolb) call errquit("wrong ncol g_dcb",0,UERR)
c
      allocate(occ_na(1:nmo))
      allocate(occ_nb(1:nmo))
      allocate(occ_ca(1:nmo))
      allocate(occ_cb(1:nmo))
c
      occ_ca = 0.0d0
      occ_cb = 0.0d0
      do ii = 1, nea
        occ_ca(wfn1_wave%icnfa(ii)) = wfn1_wave%dcnta(ii)
      enddo
      do ii = 1, neb
        occ_cb(wfn1_wave%icnfb(ii)) = wfn1_wave%dcntb(ii)
      enddo
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occ_na)
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occ_nb)
c
      call ga_distribution(wfn1_drva%g_dna,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dna,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dna,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (abs(occ_na(ii)-occ_na(jj)).gt.tol) then
            error%errna = error%errna + bufa(ii,jj)*bufb(ii,jj)
          endif
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dnb,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dnb,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dnb,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (abs(occ_nb(ii)-occ_nb(jj)).gt.tol) then
            error%errnb = error%errnb + bufa(ii,jj)*bufb(ii,jj)
          endif
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dca,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dca,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dca,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (abs(occ_ca(ii)-occ_ca(jj)).gt.tol) then
            error%errca = error%errca + bufa(ii,jj)*bufb(ii,jj)
          endif
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_distribution(wfn1_drva%g_dcb,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dcb,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dcb,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if (abs(occ_cb(ii)-occ_cb(jj)).gt.tol) then
            error%errcb = error%errcb + bufa(ii,jj)*bufb(ii,jj)
          endif
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      deallocate(occ_cb)
      deallocate(occ_ca)
      deallocate(occ_nb)
      deallocate(occ_na)
c
      error%total = error%errna + error%errnb 
     &            + error%errca + error%errcb
c
      call ga_dgop(MT_DBL,error,5,"+")
      wfn1_etrace = error
c
      return
      end function wfn1_etrace
c
c-----------------------------------------------------------------------
