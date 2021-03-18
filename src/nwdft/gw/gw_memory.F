      subroutine gw_memory(pars)
c
c     Based on tddft_memory
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
   
      double precision ga_dbl, ma_dbl, mo_eris, ao_eris, ma_mb, ga_mb
c
      integer ga_mem, ma_mem, totmopairs, isp
      integer mynpoles, maxnpoles, mynmo, maxqp
      double precision,parameter :: mbytes = 1024d0**2 
      double precision,parameter :: gbytes = 1024d0**3
      double precision dbl2mb

      dbl2mb = dble(ma_sizeof(mt_dbl,1,mt_byte))/mbytes
c
      ga_dbl = dble(ga_memory_avail())/dble(ma_sizeof(mt_dbl,1,mt_byte))
      ma_dbl = dble(ma_inquire_avail(mt_dbl))
      ma_mb = ma_dbl*dbl2mb
      ga_mb = ga_dbl*dbl2mb

      call ga_dgop(mt_dbl,ga_mb,1,"+")

      if (pars%me.eq.0) then
        write(luout,*)
        call util_print_centered(luout,'Memory Information',25,.true.)
        if (ga_uses_ma()) then
          write(luout,9000) ga_mb
        else
          write(luout,9010) ga_mb,ma_mb
        endif
      endif
c
 9000 format(10x,'Available GA+MA space size is ',F12.1,' MBs')
 9010 format(10x,'Available GA space size is    ',F12.1,' MBs',/,
     $       10x,'Available MA per process is   ',F12.1,' MBs')
 9020 format(/,10x,'Peak memory in Initialization')     
 9030 format(10x,'Peak total GA memory is       ',F12.1,' MBs',/,
     $       10x,'Peak MA memory per process is ',F12.1,' MBs')
 9040 format(/,10x,'Peak memory in Analytic Calculation')     
 9050 format(/,10x,'Peak memory in CDGW Calculation')     

      ! ERI computation and transformation
      if (pars%me.eq.0) write(luout,9020)
      totmopairs = 0
      do isp=1,pars%ipol
        totmopairs = totmopairs + pars%nocc(isp)**2 + 
     $                            pars%nocc(isp)*pars%nvir(isp) +
     $                            pars%nvir(isp)*pars%nvqp(isp) 
      enddo
      ga_mem = nbf_cd * (nbf_ao*(nbf_ao+1)/2) + pars%ipol*nbf_ao**2
      ma_mem = nbf_ao*pars%nmo + totmopairs*nbf_cd
      ga_mem = max(ga_mem,totmopairs*nbf_cd)
      ga_mem = max(ga_mem,2*nbf_cd**2)
      if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,ma_mem*dbl2mb
c
      mynpoles = maxval(pars%npoles(1:pars%ipol))/pars%nprocs + 1
      maxnpoles = maxval(pars%npoles(1:pars%ipol))
      mynmo = pars%nmo/pars%nprocs+maxval(pars%nocc(1:pars%ipol))+1
      maxqp = maxval(pars%nqp(1:pars%ipol))

      if (pars%cdgw) then
        if (pars%me.eq.0) write(luout,9050) 
        ga_mem = totmopairs*nbf_cd
        ma_mem = 10*pars%ipol*pars%nmo + 
     $    pars%ipol*pars%ngl*maxqp*mynmo +
     $    pars%ipol*maxqp*mynmo +
     $    pars%ipol*mynpoles +
     $    nbf_cd**2 + 2*nbf_cd*mynmo
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK

        ga_mem = ga_mem + nbf_cd**2
#else
        ga_mem = ga_mem + 2*nbf_cd**2
#endif
#endif
        if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,ma_mem*dbl2mb

      else

        if (pars%me.eq.0) write(luout,9040) 
        ga_mem = (totmopairs + sum(pars%npoles(1:pars%ipol)))*nbf_cd
        ga_mem = ga_mem + sum(pars%npoles(1:pars%ipol))**2

        ma_mem = 10*pars%ipol*pars%nmo +
     $           pars%ipol*maxnpoles +
     $           pars%ipol*mynpoles
        if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,ma_mem*dbl2mb

      endif
c
      return
      end subroutine
