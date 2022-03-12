      subroutine gw_memory(pars)
c
c     Obtain an estimate of the memory usage of the algorithm. This is
c     only a rough guide and likely gives overestimated values. Memory
c     allocated by MKL, SCALAPACK, ELPA is not all taken into account.
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
#ifdef SCALAPACK
#include "scaleig.fh"
#endif
#include "gw.fh"
c
      type(gw_params_t) :: pars
c
   
      double precision ga_dbl, ma_dbl, mo_eris, ao_eris, ma_mb, ga_mb
c
      integer ga_mem, ma_mem, totmopairs, isp
      integer stack_mem, heap_mem
      integer mynpoles, maxnpoles, mynmo, maxqp
      integer heap_buildiw,stack_buildiw,ga_buildiw,extramem_buildiw
      integer heap_buildrn,stack_buildrn,ga_buildrn,extramem_buildrn
      integer stack_scaleig,extramem_scaleig,extramem,nri
      double precision,parameter :: mbytes = 1024d0**2 
      double precision,parameter :: gbytes = 1024d0**3
      double precision dbl2mb

      nri = nbf_cd
      dbl2mb = dble(ma_sizeof(mt_dbl,1,mt_byte))/mbytes
c
      ga_dbl = dble(ga_memory_avail())/dble(ma_sizeof(mt_dbl,1,mt_byte))
      ma_dbl = dble(ma_inquire_avail(mt_dbl))
      ma_mb = ma_dbl*dbl2mb
      ga_mb = ga_dbl*dbl2mb
      stack_mem = ma_inquire_stack(mt_dbl)
      heap_mem = ma_inquire_heap(mt_dbl)

      call ga_dgop(1038,ga_mb,1,"min")


      if (pars%me.eq.0) then
        write(luout,*)
        call util_print_centered(luout,'Memory Information',25,.true.)
        if (ga_uses_ma()) then
          write(luout,9000) ga_mb
        else
          write(luout,9010) ga_mb,stack_mem*dbl2mb,heap_mem*dbl2mb
        endif
      endif
c
 9000 format(10x,'Available GA+MA space size is  ',F12.1,' MBs')
 9010 format(10x,'Available GA per process is    ',F12.1,' MBs',/,
     $       10x,'Available stack per process is ',F12.1,' MBs',/,
     $       10x,'Available heap per process is  ',F12.1,' MBs')
 9020 format(/,10x,'Peak memory in Initialization')     
 9030 format(14x,'Peak GA per process is         ',F12.1,' MBs',/,
     $       14x,'Peak stack per process is      ',F12.1,' MBs',/,
     $       14x,'Peak heap per process is       ',F12.1,' MBs',/,
     $       14x,'Memory not in MA per process   ',F12.1,' MBs')
 9040 format(/,10x,'Peak memory in AFF Calculation') 
 9050 format(/,10x,'Peak memory in CDGW Calculation')     

      !
      ! ERI computation and transformation
      !
      if (pars%me.eq.0) write(luout,9020)
      totmopairs = 0
      do isp=1,pars%ipol
        totmopairs = totmopairs + pars%nocc(isp)**2 + 
     $                            pars%nocc(isp)*pars%nvir(isp) +
     $                            pars%nvir(isp)*pars%nvqp(isp) 
      enddo

      ga_mem = nbf_cd**2 +           ! 2c-eri
     $         pars%ipol*nbf_ao**2 + ! movecs
     $         totmopairs*nbf_cd +   ! 3c-eri
     $         (nbf_ao-pars%nocc(1))**2

      heap_mem = totmopairs*nbf_cd/pars%nprocs + ! 3c-eri
     $         pars%nstates(1)*nbf_cd*nbf_ao_mxnbf_cn/pars%nprocs + ! MObuf
     $         nbf_ao_mxnbf_cn**2*nbf_cd/pars%nprocs +              ! Buf
     $         nbf_ao*pars%nmo + ! movecs
     $         nbf_cd**2         ! 2c-eri
      stack_mem = 2*(nbf_ao-pars%nocc(1))**2 + nri/pars%nprocs
      ma_mem = heap_mem + stack_mem
      extramem = 0

      ga_mem = ga_mem/pars%nprocs
      
      if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,
     $    stack_mem*dbl2mb,heap_mem*dbl2mb,extramem*dbl2mb
c
      mynpoles = maxval(pars%mynpoles(1:pars%ipol))
      maxnpoles = maxval(pars%npoles(1:pars%ipol))
      mynmo = pars%nmo/pars%nprocs+maxval(pars%nocc(1:pars%ipol))+1
      maxqp = maxval(pars%nqp(1:pars%ipol))

      !
      ! CD-GW algorithm
      !
      if (pars%cdgw) then
        if (pars%me.eq.0) write(luout,9050) 

        ga_mem = totmopairs*nbf_cd   ! 3c-MO ERIs

        stack_mem = pars%ipol*pars%nmo + ! Sigma
     $              pars%ipol*pars%nmo + ! Sigma_x
     $              pars%ipol*pars%nmo*3 ! Eref, Emod, Enew

        heap_mem = pars%ipol*mynmo*maxqp*pars%ngl + ! iW
     $             pars%ipol*(mynmo+1)*maxqp   ! iWidx

        if (pars%evgw.or.pars%evgw0) then
          stack_mem = stack_mem + pars%ipol*pars%nmo*pars%eviter*2
        endif

        stack_mem =  stack_mem + pars%ipol*mynpoles ! Eia
        ma_mem = stack_mem + heap_mem

        ga_mem = ga_mem + nri*pars%nmo ! Solutions

        stack_scaleig = 0
        extramem_scaleig = 0
#ifdef SCALAPACK
        if (yev.eq.'v') then
          stack_scaleig = nri**2
        elseif (yev.eq.'x') then
          stack_scaleig = 2*nri**2
        elseif (yev.eq.'d') then
          stack_scaleig = 2*nri**2
#ifdef USE_ELPA
        elseif (yev.eq.'e') then
          extramem_scaleig = nri**2
#endif
        endif
#endif 
        !
        ! *** BUILDIW ***
        !
        extramem_buildiw = maxnpoles*ipol ! factor
        heap_buildiw = (nri*(nri+1))/2 + ! Pi             
     $                  nri*mynmo      + ! Solutions
#ifdef GWCMPLX
#else
     $                 (nri+(nri+1))/2 + ! VX
     $                  nri*mynmo        ! tmp 
        ga_buildiw    = 3*nri**2 + nri
        stack_buildiw = nri
#endif

        !
        ! *** BUILDRN ***
        !
        extramem_buildrn = mynpoles*pars%ipol*2      + ! factor, dfactor
     $                      nri*10                     ! Minres
        stack_buildrn = nri                            ! Pi
        heap_buildrn  = 0
        ga_buildrn = 3*nri**2 + nri

        ga_mem = ga_mem + max(ga_buildiw,ga_buildrn)
        stack_mem = stack_mem + max(stack_buildiw,stack_buildrn) +
     $              stack_scaleig
        heap_mem = heap_mem + max(heap_buildiw,heap_buildrn)
        extramem = max(extramem_buildiw,extramem_buildrn) +
     $             extramem_scaleig
        ma_mem = stack_mem + heap_mem
        ga_mem = ga_mem/pars%nprocs

        if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,
     $    stack_mem*dbl2mb,heap_mem*dbl2mb,extramem*dbl2mb


      !
      !  Spectral decomposition algorithm
      !
      else

        if (pars%me.eq.0) write(luout,9040) 
        ga_mem = (totmopairs + sum(pars%npoles(1:pars%ipol)))*nbf_cd
        ga_mem = ga_mem + sum(pars%npoles(1:pars%ipol))**2

        heap_mem = 10*pars%ipol*pars%nmo +
     $           pars%ipol*maxnpoles +
     $           pars%ipol*mynpoles
        ga_mem = ga_mem/pars%nprocs
        stack_mem = sum(pars%npoles(1:pars%ipol))**2
        extramem = stack_mem
     
        if (pars%me.eq.0) write(luout,9030) ga_mem*dbl2mb,
     $                   stack_mem*dbl2mb,heap_mem*dbl2mb,
     $                   extramem*dbl2mb

      endif
c
      return
      end subroutine
