      subroutine tddft_memory(rtdb,nroots,maxvecs,maxiter,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,nov,tda,algorithm)
c
c $Id: tddft_memory.F,v 1.4 2002-04-28 19:09:13 sohirata Exp $
c
c Itemize all the outstanding memory usage and determine
c how many trial and product vectors we can handle with the
c available memory and alter maxvecs if necessary.
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer rtdb             ! Runtime database
      integer nbf_ao           ! Number of AO basis functions
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta electrons
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nroots           ! Number of roots sought
      integer maxvecs          ! Maximum number of trial vectors
      integer maxiter          ! Maximum number of iterations
      integer ga_dble          ! GA memory size in doubles
      integer ma_dble          ! MA memory size in doubles (stack)
      integer nov(2)           ! Number of occupied virtual pairs
      logical tda              ! True if Tamm-Dancoff approximation
      integer algorithm        ! =1 incore-multiple fock build; =2 incore-single fock build
                               ! =3 disk-multiple fock build;   =4 disk-single fock build
c
      integer i,j,k,n,m
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c -----------------------------
c Retrieve the memory size info
c -----------------------------
c
      ga_dble=ga_memory_avail()/ma_sizeof(mt_dbl,1,mt_byte)
      ma_dble=ma_inquire_avail(mt_dbl)
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'Memory Information',20,.true.)
        if (ga_uses_ma()) then
          write(LuOut,9000) ga_dble
        else
          write(LuOut,9010) ga_dble,ma_dble
        endif
      endif
 9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
 9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
     1       10x,'Available MA space size is    ',i10,' doubles')
c
c -------------------------------------
c Detemine the length of a trial vector
c Occ*Virt
c -------------------------------------
c
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
      enddo
      if (ipol.eq.1) nov(2)=0
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        if (ipol.eq.1) then
          write(LuOut,9100) nov(1)
        else
          write(LuOut,9100) nov(1),nov(2)
        endif
      endif
 9100 format(10x,'Length of a trial vector is     ',2i8)
c
c ----------------------
c Estimate peak MA usage
c ----------------------
c
      n=0
c
c (A+B) (A-B) and their eigenvector matrices
c
      n=n+2*(nroots*30)**2
      if (.not.tda) n=n+2*(nroots*30)**2
c
c Work arrays in tddft_diagon
c
      n=n+3*(nroots*30)**2
      n=n+10*(nroots*30)
c
c Memory buffer (10 MB) for DRA
c
      if ((algorithm.eq.3).or.(algorithm.eq.4))
     1  n=n+int(1.0d7/dfloat(ma_sizeof(mt_dbl,1,mt_byte)))
c
c ------------------------------------------------------
c Estimate peak GA usage besides g_trials and g_products
c ------------------------------------------------------
c
      m=0
c
c MO vectors
c
      m=m+nbf_ao*nbf_ao*ipol
c
c Trial and product density matrices
c 
      m=m+2*nbf_ao*nbf_ao*ipol
      if (.not.tda) m=m+nbf_ao*nbf_ao*ipol
c
c Residual vectors
c
      m=m+(nov(1)+nov(2))*nroots
      if (.not.tda) m=m+(nov(1)+nov(2))*nroots
c
c Symmetrized density matrix and NGA trial and product vectors
c
      if (algorithm.eq.1) then
        m=m+4*nbf_ao*nbf_ao*ipol*(2*nroots)
        if (ipol.eq.2) m=m+4*nbf_ao*nbf_ao*ipol*(2*nroots)
      else if (algorithm.eq.2) then
        m=m+3*nbf_ao*nbf_ao*ipol
        if (ipol.eq.2) m=m+4*nbf_ao*nbf_ao*ipol
      else if (algorithm.eq.3) then
        m=m+4*nbf_ao*nbf_ao*ipol*(2*nroots)
        if (ipol.eq.2) m=m+4*nbf_ao*nbf_ao*ipol*(2*nroots)
        m=m+2*(nov(1)+nov(2))
      else if (algorithm.eq.4) then
        m=m+3*nbf_ao*nbf_ao*ipol
        if (ipol.eq.2) m=m+4*nbf_ao*nbf_ao*ipol
        m=m+2*(nov(1)+nov(2))
      endif
c
c Check GA and MA remaining space at this point
c
      if (n.gt.ma_dble) call errquit
     1  ('tddft_memory: MA space too small',ma_dble)
      if (ga_uses_ma()) then
        ga_dble=ga_dble-m-n
      else
        ga_dble=ga_dble-m
      endif
      if (ga_dble.le.0) call errquit
     1  ('tddft_memory: GA space too small',ga_dble)
c
c Compute GA usage per trial vector
c
      j=0
      do i=1,ipol
        j=j+nov(i)
      enddo
      ! TDA: 2 for g_trials,g_prod_p, and an extra for ddot operation
      ! RPA: 3 for g_trials,g_prod_p,g_prod_m, and an extra for ddot operation
      if (tda) then
        j=j*(2+1)
      else
        j=j*(3+1)
      endif
c
c ------------------------------------------------------------------------------
c Algorithms 1 & 2: allocate 90% of remaining memory for trial & product vectors
c ------------------------------------------------------------------------------
c
      if ((algorithm.eq.1).or.(algorithm.eq.2)) then
        k=int(0.9d0*dfloat(ga_dble)/dfloat(j))
        if (k.lt.maxvecs) then
          if (nodezero.and.util_print('memory information',
     1      print_default)) write(LuOut,9110) maxvecs,k
          maxvecs=k
          if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1      call errquit('tddft_memory: failed writing to rtdb',0)
        endif
        m=m+maxvecs*j
c
c -------------------------------------------------------------------------------
c Algorithms 3 & 4: 10 GB of disk space is assumed (see dra_init in tddft_energy)
c -------------------------------------------------------------------------------
c
      else
        k=int(1.0d10/dfloat(ma_sizeof(mt_dbl,1,mt_byte))/dfloat(j))
        if (k.lt.maxvecs) then
          if (nodezero.and.util_print('memory information',
     1      print_default)) write(LuOut,9110) maxvecs,k
          maxvecs=k
          if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1      call errquit('tddft_memory: failed writing to rtdb',0)
        endif
      endif
      if (maxvecs.lt.nroots) call errquit
     1  ('Too small a subspace',maxvecs)
 9110 format(10x,'Max subspace size reduced from  ',i8,' to ',i8)
c
c ----------------------------
c Print estimated memory usage
c ----------------------------
c
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        write(LuOut,9020) m,n
        if ((algorithm.eq.3).or.(algorithm.eq.4))
     1    write(LuOut,9030) maxvecs*j
      endif
 9020 format(10x,'Estimated peak GA usage is    ',i10,' doubles',/,
     1       10x,'Estimated peak MA usage is    ',i10,' doubles')
 9030 format(10x,'Estimated peak DRA usage is   ',i10,' doubles')
      call util_flush(LuOut)
c
c ------
c Return
c ------
c
      return
      end
