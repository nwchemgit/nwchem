c-----------------------------------------------------------------------
c
      double precision function wfn1_trace(wfn1_drva,wfn1_drvb)
      implicit none
C>
C> \brief Compute the trace of the gradients
C>
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "wfn1_drv.fh"
c
      type(wfn1_drv), intent(in) :: wfn1_drva
      type(wfn1_drv), intent(in) :: wfn1_drvb
c
      integer ilo, ihi, jlo, jhi
      integer ii, jj
      integer itypa, nrowa, ncola
      integer itypb, nrowb, ncolb
c
      integer iproc
      integer nproc
c
      double precision :: trace
c
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
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
      trace = 0.0d0
      call ga_distribution(wfn1_drva%g_dna,iproc,ilo,ihi,jlo,jhi)
      allocate(bufa(ilo:ihi,jlo:jhi))
      allocate(bufb(ilo:ihi,jlo:jhi))
      call ga_get(wfn1_drva%g_dna,ilo,ihi,jlo,jhi,bufa,ihi-ilo+1)
      call ga_get(wfn1_drvb%g_dna,ilo,ihi,jlo,jhi,bufb,ihi-ilo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          trace = trace + bufa(ii,jj)*bufb(ii,jj)
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
          trace = trace + bufa(ii,jj)*bufb(ii,jj)
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
          trace = trace + bufa(ii,jj)*bufb(ii,jj)
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
          trace = trace + bufa(ii,jj)*bufb(ii,jj)
        enddo
      enddo
      deallocate(bufb)
      deallocate(bufa)
c
      call ga_dgop(MT_DBL,trace,1,"+")
      wfn1_trace = trace
c
      return
      end function wfn1_trace
c
c-----------------------------------------------------------------------
