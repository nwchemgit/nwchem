      subroutine tddft_memory(rtdb,nroots,maxvecs,maxiter,
     1  ipol,nbf_ao,nocc,nmo,nfc,nfv,nov,tda,algorithm)
c
c $Id$
c
c Itemize all the outstanding memory usage and determine
c how many trial and product vectors we can handle with the
c available memory and alter maxvecs if necessary.
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
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
      double precision ga_dble ! GA memory size in doubles
      double precision ma_dble ! MA memory size in doubles (stack)
      integer nov(2)           ! Number of occupied virtual pairs
      logical tda              ! True if Tamm-Dancoff approximation
      integer algorithm        ! =1 incore-multiple fock build; =2 incore-single fock build
                               ! =3 disk-multiple fock build;   =4 disk-single fock build
c
      integer i,j,k
      double precision n,m
      integer mxtrialf
      parameter(mxtrialf=10)
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c -----------------------------
c Retrieve the memory size info
c -----------------------------
c
      ga_dble=dble(ga_memory_avail())
     1  /dble(ma_sizeof(mt_dbl,1,mt_byte))
      ma_dble=dble(ma_inquire_avail(mt_dbl))
      call ga_dgop(msg_tdmem,ga_dble,1,'+')
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'Memory Information',20,.true.)
        if (ga_uses_ma()) then
          write(LuOut,9000) int(ga_dble)
        else
          write(LuOut,9010) int(ga_dble),int(ma_dble)
c 
c ====================================================================
c Warning! when ga_uses_ma()=.false., "shmmax" may become the limiting
c          memory size, rather than ga_memory_avail.  In that case, 
c          despite the fact that ga_memory_avail gives a memory size
c          sufficient for the calculation, an out-of-memory error can 
c          occur.  Check shmmax by "cat /proc/sys/kernel/shmmax".
c ====================================================================
c
        endif
      endif
 9000 format(10x,'Available GA+MA space size is ',i14,' doubles')
 9010 format(10x,'Available GA space size is    ',i14,' doubles',/,
     1       10x,'Available MA space size is    ',i14,' doubles')
c
c Determine the length of the trial vector: Occ*Virt
      call tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
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
c ========================================================
c Choose the best algorithm (if not specified by the user)
c ========================================================
c
      if (algorithm.eq.0) then
c peak MA usage
        n=0.0d0
        n=n+dble(2*(nroots*mxtrialf)**2)
        if (.not.tda) n=n+dble(2*(nroots*mxtrialf)**2)
        n=n+dble(3*(nroots*mxtrialf)**2)
        n=n+dble(10*(nroots*mxtrialf))
c peak GA usage except for g_trials and g_products and NGA
        m=0.0d0
        m=m+dble(nbf_ao*nbf_ao*ipol)
        m=m+dble(2*nbf_ao*nbf_ao*ipol)
        if (.not.tda) m=m+dble(nbf_ao*nbf_ao*ipol)
        if (ga_uses_ma()) then
          m=ga_dble-m-n
        else
          m=ga_dble-m
        endif
        j=0
        do i=1,ipol
          j=j+nov(i)
        enddo
        if (tda) then
          j=j*(2+1)
        else
          j=j*(3+1)
        endif
c algorithm decision
        if (m.gt.dble(j*(nroots*mxtrialf)
     1           +((2*nbf_ao*ipol)**2)*(2*nroots))) then
          algorithm = 1
          if (nodezero.and.util_print('memory information',
     1      print_default))
     2      write(LuOut,9140) 'Incore multiple tensor contraction' 
        else if (m.gt.dble(j*(nroots*mxtrialf)
     1           +((2*nbf_ao*ipol)**2))) then
          algorithm = 2
          if (nodezero.and.util_print('memory information',
     1      print_default))
     2      write(LuOut,9140) 'Incore single tensor contraction' 
        else if (m.gt.1.0d7/dble(ma_sizeof(mt_dbl,1,mt_byte))
     1           +dble(((2*nbf_ao*ipol)**2)*(2*nroots))) then
          algorithm = 3
          if (nodezero.and.util_print('memory information',
     1      print_default))
     2      write(LuOut,9140) 'Disk-based multiple tensor contraction' 
        else if (m.gt.1.0d7/dble(ma_sizeof(mt_dbl,1,mt_byte))
     1           +dble(((2*nbf_ao*ipol)**2))) then
          algorithm = 4
          if (nodezero.and.util_print('memory information',
     1      print_default))
     2      write(LuOut,9140) 'Disk-based single tensor contraction' 
        else
          call errquit
     1    ('tddft_memory: no algorithm for the memory size',0, MEM_ERR)
        endif
      endif
 9140 format(3x,'       Algorithm : ',a)
c
c ----------------------
c Estimate peak MA usage
c ----------------------
c
      n=0.0d0
c
c (A+B) (A-B) and their eigenvector matrices
c
      n=n+dble(2*(nroots*mxtrialf)**2)
      if (.not.tda) n=n+dble(2*(nroots*mxtrialf)**2)
c
c Work arrays in tddft_diagon
c
      n=n+dble(3*(nroots*mxtrialf)**2)
      n=n+dble(10*(nroots*mxtrialf))
c
c Memory buffer (10 MB) for DRA
c
      if ((algorithm.eq.3).or.(algorithm.eq.4))
     1  n=n+1.0d7/dble(ma_sizeof(mt_dbl,1,mt_byte))
c
c ------------------------------------------------------
c Estimate peak GA usage besides g_trials and g_products
c ------------------------------------------------------
c
      m=0.0d0
c
c MO vectors
c
      m=m+dble(nbf_ao*nbf_ao*ipol)
c
c Trial and product density matrices
c 
      m=m+dble(2*nbf_ao*nbf_ao*ipol)
      if (.not.tda) m=m+dble(nbf_ao*nbf_ao*ipol)
c
c Residual vectors
c
      m=m+dble((nov(1)+nov(2))*nroots)
      if (.not.tda) m=m+dble((nov(1)+nov(2))*nroots)
c
c Symmetrized density matrix and NGA trial and product vectors
c
      if (algorithm.eq.1) then
        m=m+dble(4*nbf_ao*nbf_ao*ipol*(2*nroots))
        if (ipol.eq.2) m=m+dble(4*nbf_ao*nbf_ao*ipol*(2*nroots))
      else if (algorithm.eq.2) then
        m=m+dble(3*nbf_ao*nbf_ao*ipol)
        if (ipol.eq.2) m=m+dble(4*nbf_ao*nbf_ao*ipol)
      else if (algorithm.eq.3) then
        m=m+dble(4*nbf_ao*nbf_ao*ipol*(2*nroots))
        if (ipol.eq.2) m=m+dble(4*nbf_ao*nbf_ao*ipol*(2*nroots))
        m=m+dble(2*(nov(1)+nov(2)))
      else if (algorithm.eq.4) then
        m=m+dble(3*nbf_ao*nbf_ao*ipol)
        if (ipol.eq.2) m=m+dble(4*nbf_ao*nbf_ao*ipol)
        m=m+dble(2*(nov(1)+nov(2)))
      endif
c
c Check GA and MA remaining space at this point
c
      if (n.gt.ma_dble) call errquit
     1  ('tddft_memory: MA space too small',
     . int(n-ma_dble), MA_ERR)
      if (ga_uses_ma()) then
        ga_dble=ga_dble-m-n
      else
        ga_dble=ga_dble-m
      endif
      if (ga_dble.le.0.0d0) call errquit
     1  ('tddft_memory: GA space too small',int(ga_dble), GA_ERR)
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
        k=int(0.9d0*ga_dble/dble(j))
        if (k.lt.maxvecs) then
          if (nodezero.and.util_print('memory information',
     1      print_default)) write(LuOut,9110) maxvecs,k
          maxvecs=k
          if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1      call errquit('tddft_memory: failed writing to rtdb',0,
     &       RTDB_ERR)
        endif
        m=m+dble(maxvecs*j)
c
c -------------------------------------------------------------------------------
c Algorithms 3 & 4: 10 GB of disk space is assumed (see dra_init in tddft_energy)
c -------------------------------------------------------------------------------
c
      else
        k=int(1.0d10/dble(ma_sizeof(mt_dbl,1,mt_byte))/dble(j))
        if (k.lt.maxvecs) then
          if (nodezero.and.util_print('memory information',
     1      print_default)) write(LuOut,9110) maxvecs,k
          maxvecs=k
          if (.not.rtdb_put(rtdb,'tddft:maxvecs',mt_int,1,maxvecs))
     1      call errquit('tddft_memory: failed writing to rtdb',0,
     &       RTDB_ERR)
        endif
      endif
      if (tda) then
        if (maxvecs.lt.nroots+1) call errquit
     1    ('Too small a subspace',maxvecs, MEM_ERR)
      else
        if (maxvecs.lt.2*nroots+2) call errquit
     1    ('Too small a subspace',maxvecs, MEM_ERR)
      endif
 9110 format(10x,'Max subspace size reduced from  ',i8,' to ',i8)
c
c ----------------------------
c Print estimated memory usage
c ----------------------------
c
      if (nodezero.and.util_print('memory information',
     1  print_default)) then
        write(LuOut,9020) int(m),int(n)
        if ((algorithm.eq.3).or.(algorithm.eq.4))
     1    write(LuOut,9030) maxvecs*j
      endif
 9020 format(10x,'Estimated peak GA usage is    ',i14,' doubles',/,
     1       10x,'Estimated peak MA usage is    ',i14,' doubles')
 9030 format(10x,'Estimated peak DRA usage is   ',i14,' doubles')
      if (nodezero) call util_flush(LuOut)
c
c ------
c Return
c ------
c
      return
      end
c
c     ---------------------------------------
c     Determine the length of a trial vector
c     ---------------------------------------
c
      subroutine tddft_lentrialvec(ipol,nocc,nmo,nfc,nfv,nov)
c
      implicit none
c
#include "stdio.fh"
c
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer nocc(2)          ! Number of alpha, beta electrons
      integer nmo(2)           ! Number of alpha, beta electrons
      integer nfc(2),nfv(2)    ! Number of frozen cores/virtuals
      integer nov(2)           ! Number of occupied virtual pairs
c
      integer i
c
      do i=1,ipol
        nov(i)=(nmo(i)-nfv(i)-nocc(i))*(nocc(i)-nfc(i))
c        write(luout,*) i,nov(i),nocc(i)
      enddo
      if (ipol.eq.1) nov(2)=0
c
      return
      end
