      subroutine tddft_memory(rtdb,nroots,maxvecs,maxiter,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,nov,tda,algorithm)
c
c $Id: tddft_memory.F,v 1.3 2002-04-10 19:17:19 sohirata Exp $
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
      integer algorithm        ! =1 multiple fock build; =2: conventional
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
c
c ----------------------
c Estimate peak MA usage
c ----------------------
c
      n=0
c
c MO eigenvalues, and occupancies
c
      n=n+nbf_ao*ipol
      n=n+nbf_ao*ipol
c
c (A+B) (A-B) eigenvalues, last eigenvalues, 2 residual vectors
c
      n=n+3*(nroots*30*2)
      do i=1,ipol
        n=n+2*((nmo(i)-nfv(i))-(nocc(i)-nfc(i)))*(nocc(i)-nfc(i))*nroots
      enddo
c
c Work arrays in tddft_transform
c
      n=n+3*nbf_ao**2
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
c (A+B) (A-B) and their eigenvector matrices
c
      m=m+2*(nroots*30)**2
      if (.not.tda) m=m+2*(nroots*30)**2
c
c Symmetrized density matrix and NGA trial and product vectors
c
      m=m+nbf_ao*nbf_ao*ipol
      if (algorithm.eq.1) then
        m=m+4*nbf_ao*nbf_ao*ipol*(2*nroots)
      else
        m=m+2*nbf_ao*nbf_ao*ipol
      endif
c
c ------------------------------------------------------------
c Allocate 90% of remaining memory for trial & product vectors
c ------------------------------------------------------------
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
      k=int(0.9d0*dfloat(ga_dble)/dfloat(j))
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        if (ipol.eq.1) then
          write(LuOut,9100) nov(1)
        else
          write(LuOut,9100) nov(1),nov(2)
        endif
      endif
      if (k.lt.maxvecs) then
        if (nodezero.and.util_print('memory information',
     1    print_default)) write(LuOut,9110) maxvecs,k
        maxvecs=k
        if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1    call errquit('tddft_memory: failed writing to rtdb',0)
      endif
      if (maxvecs.lt.nroots) call errquit
     1  ('Too small a subspace',maxvecs)
 9100 format(10x,'Length of a trial vector is     ',2i8)
 9110 format(10x,'Max subspace size reduced from  ',i8,' to ',i8)
c
c ----------------------------
c Print estimated memory usage
c ----------------------------
c
      m=m+maxvecs*j
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        write(LuOut,9020) m,n
      endif
 9020 format(10x,'Estimated peak GA usage is    ',i10,' doubles',/,
     1       10x,'Estimated peak MA usage is    ',i10,' doubles')
      call util_flush(LuOut)
c
c ------
c Return
c ------
c
      return
      end
