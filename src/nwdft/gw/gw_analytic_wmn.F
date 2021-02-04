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
      integer nmo,npoles,nvir
      integer idx,kmo

      logical,external :: ga_create, ga_destroy

      double precision :: scr(nri)
#ifdef USE_OPENMP
      integer iMaxThreads
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif      

      nmo = pars%nmo
      nocc = pars%nocc(isp)
      nvir = pars%nvir(isp)
      oo = pars%g_erioo(isp)
      ov = pars%g_eriov(isp)
      vv = pars%g_erivv(isp)
      erim = pars%g_erim(isp)
      npoles = pars%npoles(isp)

      if (.not.ga_create(mt_dbl,nri,nmo,'temp',nri,0,g_a))
     $  call errquit('gw_analytic_wmn: cannot create array',0,GA_ERR)

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

      call ga_dgemm('t','n',nmo,npoles,nri,dsqrt(2d0),g_a,erim,0d0,wmn)
      call ga_elem_multiply(wmn,wmn,wmn)

      if (.not.ga_destroy(g_a))
     $  call errquit('gw_analytic_wmn: can''t destroy array',0,GA_ERR)

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      end subroutine
