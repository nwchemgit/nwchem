      subroutine gw_analytic_wmn(pars,wmn,imo,isp,nri)
      implicit none
#include "gw.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: imo, isp, wmn, nri

      integer oo,ov,vv,erim,nocc
      integer g_a
      integer ilo,ihi,jlo,jhi,iocc
      integer nmo,npoles

      integer,external :: ga_create, ga_destroy

      nmo = pars%nmo
      nocc = pars%nocc(isp)
      oo = pars%g_erioo(isp)
      ov = pars%g_eriov(isp)
      vv = pars%g_erivv(isp)
      erim = pars%g_erim(isp)
      npoles = pars%npoles(isp)

      if (.not.ga_create(mt_dbl,nri,nmo,'temp',nri,0,g_a))
     $  call errquit('gw_analytic_wmn: cannot create array',0,GA_ERR)

      if (imo.le.nocc) then
        ilo = (imo-1)*nocc + 1
        ihi = imo*nocc
        call ga_copy_patch('n',oo,1,nri,ilo,ihi,g_a,1,nri,1,nocc)
        ilo = (imo-1)*(nmo-nocc) + 1
        ihi = imo*(nmo-nocc)
        call ga_copy_patch('n',ov,1,nri,ilo,ihi,g_a,1,nri,nocc+1,nmo)
      else
        do iocc=1,nocc
          ilo = (iocc-1)*(nmo-nocc) + imo - nocc
          call ga_copy_patch('n',ov,1,nri,ilo,ilo,g_a,1,nri,iocc,iocc)
        enddo
        ilo = (imo-nocc-1)*(nmo-nocc) + 1
        ihi = (imo-nocc)*(nmo-nocc)
        call ga_copy_patch('n',vv,1,nri,ilo,ihi,g_a,1,nri,nocc+1,nmo)
      endif

      call ga_dgemm('t','n',nmo,npoles,nri,dsqrt(2d0),g_a,erim,0d0,wmn)
      call ga_elem_multiply(wmn,wmn,wmn)

      if (.not.ga_destroy(g_a))
     $  call errquit('gw_analytic_wmn: can''t destroy array',0,GA_ERR)

      end subroutine
