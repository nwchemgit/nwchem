c
c
c
c
c
       subroutine moints_mp2_eval(nocc,nvir,orbe,g_k,emp2,epair)
       implicit none
#include "global.fh"
#include "mafdecls.fh"
       integer nocc, nvir
       double precision orbe(nocc+nvir), emp2, epair(*)
       integer g_k
       integer k_exch,ld_exch
       integer i,j,ij,vstart,nnt,nbf,ndim
       integer ilo,ihi,jlo,jhi,my_id
       integer g_epair
       double precision xx,denom
       logical status
       double precision moint_epair_eval
       
       nnt = (nocc*(nocc+1))/2
       if (.not.ga_create(MT_DBL,nnt,1,'epair',1,1,g_epair)) call
     $      errquit('moints_mp2_eval: cannot allocate',0)
       nbf = nocc + nvir
       my_id = ga_nodeid()
       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)
       if (ilo.ne.1) call
     $   errquit('moints_mp2_eval: wrond distrib for exch',0)
c
c Set to handle different exchange matrix sizes
c
       if (ihi.eq.(nbf*nbf)) then
         ndim = nocc+nvir
         vstart = nocc
       else if (ihi.eq.(nvir*nvir)) then
         ndim = nvir
         vstart = 0
       else
         call errquit('moints_mp2_eval: odd exchange matrix size',0)
       endif
c
c Loop over occupied orbital pairs
c
       do i=1,nocc
         do j=1,i
           ij = (i*(i-1))/2 + j
           if ((ij.ge.jlo).and.(ij.le.jhi)) then
             denom = orbe(i) + orbe(j)
             call ga_access(g_k,1,ihi,ij,ij,k_exch,ld_exch)
             xx = moint_epair_eval(ndim,vstart,nocc,nvir,
     $                             dbl_mb(k_exch),orbe,denom)
             call ga_release(g_k,1,ihi,ij,ij)
             if (i.eq.j) xx = xx/2.d0
             call ga_put(g_epair,ij,ij,1,1,xx,1)
           endif
         enddo
       enddo
c
c Get global pair-energies into local array
c
       call ga_sync()
       emp2 = 0.d0
       do ij=1,nnt
         call ga_get(g_epair,ij,ij,1,1,xx,1)
         epair(ij) = xx
         emp2 = emp2 + xx
       enddo
c
c Clean and print
c
       if (.not.ga_destroy(g_epair)) call
     $      errquit('moints_mp2_eval_a: cannot destroy global',0)
       if (ga_nodeid().eq.0) then
         write(6,*)
         call util_print_centered(6,'MP2 and Pair Energies',40,.true.)
         write(6,901) emp2
 901     format(/,15x,'MP2 correction:',f15.10,/)
         write(6,902) ((epair((i*(i-1))/2+j),i,j,j=1,i),i=1,nocc)
 902     format(4(f10.5,1x,'(',2i3,')',2x))
       endif
         
       return
       end





       double precision function moint_epair_eval( ndim, vstart, nocc,
     $                                          nvir, xk, orbe, denom )
       implicit none
       integer ndim, vstart, nocc, nvir
       double precision xk(ndim,ndim), orbe(nocc+nvir), denom
       double precision xx, yy
       integer k,l,kk,ll

       xx = 0.d0
#ifdef DEBUG
       write(6,*)
#endif
       do kk=1,nvir
         do ll=1,nvir
           k = vstart + kk
           l = vstart + ll
           yy = xk(l,k)*(4.d0*xk(l,k) - 2.d0*xk(k,l))
           xx = xx + yy/(denom - orbe(nocc+kk) - orbe(nocc+ll))
#ifdef DEBUG
           write(6,933) k,l,xk(l,k),xk(k,l),yy,xx
 933       format(2i5,2f12.6,5x,2f12.6)
#endif
         enddo
       enddo
       moint_epair_eval = xx
       return
       end







