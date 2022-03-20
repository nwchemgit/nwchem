      subroutine gw_analytic_wmn(pars,wmn,imo,isp,nri)
c
c     Complete Equation (36). On entry, the global array ERIm has the
c     ERIs modified with the sum of the (X+Y) eigenvectors. Here, we
c     do the reamining sum over the fitting basis functions.
c
      implicit none
#include "gw.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      type(gw_params_t)  :: pars
      integer,intent(in) :: imo, isp, wmn, nri

      integer oo,ov,vv,erim,nocc
      integer g_a
      integer ilo,ihi,jlo,jhi,iocc
      integer nmo,npoles,nvir
      integer idx,kmo

      logical,external :: ga_create, ga_destroy

      double precision :: scr(nri)
      double precision :: factor

      nmo = pars%nmo
      nocc = pars%nocc(isp)
      nvir = pars%nvir(isp)
      oo = pars%g_erioo(isp)
      ov = pars%g_eriov(isp)
      vv = pars%g_erivv(isp)
      erim = pars%g_erim                     ! Spin-independent
      npoles = sum(pars%npoles(1:pars%ipol)) ! Total number of poles

      if (.not.ga_create(mt_dbl,nri,nmo,'temp',nri,0,g_a))
     $  call errquit('gw_analytic_wmn: cannot create array',0,GA_ERR)

      !
      ! Gather all ERI vectors for all MO pairs involving the target
      ! orbital
      !
      if (imo.le.nocc) then

        do kmo=1,nocc
          idx = (imo-1)*nocc+kmo-1
          if (pars%me.ne.mod(idx,pars%nprocs)) cycle
          idx = (idx/pars%nprocs) + pars%oolo(isp)
          call ga_get(oo,1,nri,idx,idx,scr,nri)
          call ga_put(g_a,1,nri,kmo,kmo,scr,nri)
        enddo
        do kmo=nocc+1,nmo
          idx = (imo-1)*nvir+kmo-nocc-1
          if (pars%me.ne.mod(idx,pars%nprocs)) cycle
          idx = (idx/pars%nprocs) + pars%ovlo(isp)
          call ga_get(ov,1,nri,idx,idx,scr,nri)
          call ga_put(g_a,1,nri,kmo,kmo,scr,nri)
        enddo
      else
        do kmo=1,nocc
          idx = (kmo-1)*nvir+imo-nocc-1
          if (pars%me.ne.mod(idx,pars%nprocs)) cycle
          idx = (idx/pars%nprocs) + pars%ovlo(isp)
          call ga_get(ov,1,nri,idx,idx,scr,nri)
          call ga_put(g_a,1,nri,kmo,kmo,scr,nri)
        enddo
        do kmo=1,nvir
          idx = (imo-nocc-1)*nvir+kmo-1
          if (pars%me.ne.mod(idx,pars%nprocs)) cycle
          idx = (idx/pars%nprocs) + pars%vvlo(isp)
          call ga_get(vv,1,nri,idx,idx,scr,nri)
          call ga_put(g_a,1,nri,kmo+nocc,kmo+nocc,scr,nri)
        enddo
      endif

      call ga_sync()
      factor = 1d0
      if (pars%ipol.eq.1) factor = dsqrt(2d0)

      !
      ! Contract the ERIs with the modified ERIm
      !
      call ga_dgemm('t','n',nmo,npoles,nri,factor,g_a,erim,0d0,wmn)

      !
      ! We will need the square of these quantities for Equation (37)
      !
      call ga_elem_multiply(wmn,wmn,wmn)

      if (.not.ga_destroy(g_a))
     $  call errquit('gw_analytic_wmn: can''t destroy array',0,GA_ERR)


      end subroutine
