c
c
c
c
c
       subroutine moints_mp2_eval(nocc,nvir,orbe,g_k,emp2,epair)
C$Id: moints_mp2.F,v 1.18 1996-01-26 23:49:42 rg240 Exp $
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
       double precision moints_epair_eval

       nnt = (nocc*(nocc+1))/2
       if (.not.ga_create(MT_DBL,nnt,1,'epair',1,1,g_epair)) call
     $      errquit('moints_mp2_eval: cannot allocate',0)
       nbf = nocc + nvir
       my_id = ga_nodeid()
       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)

#ifdef DEBUG
       if (ga_nodeid().eq.0) write(6,944)
       do i=0,ga_nnodes()-1
         if (i.eq.my_id) then
           write(6,923) i,ilo,ihi,jlo,jhi
 923       format('moints_mp2_eval: distribution','[',i5,']',3x,
     $            2(2i5,3x))
           call util_flush(6)
	 endif
         call ga_sync()
       enddo
       if (ga_nodeid().eq.0) write(6,944)
 944   format(/,'===================================================',/)
#endif
c
c Check if this node has any data
c
       if ((ilo.ne.1).and.(ihi.ne.-1)) call 
     $ errquit('moints_mp2_eval: wrong distribution for exchange',my_id)
c
c Set to handle different exchange matrix sizes
c
       if (ihi.eq.(nbf*nbf)) then
         ndim = nocc+nvir
         vstart = nocc
       else if ((ihi.eq.(nvir*nvir)).or.(ihi.eq.(nvir*nbf))) then
         ndim = nvir
         vstart = 0
       else if (ihi.ne.-1) then
         call errquit('moints_mp2_eval: odd exchange matrix size',ihi)
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
             xx = moints_epair_eval(ndim,vstart,nocc,nvir,
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
C         write(6,*)
C         call util_print_centered(6,'MP2 and Pair Energies',40,.true.)
         write(6,901) emp2
 901     format(//,15x,'MP2 correction:',f15.10,/)
C         write(6,902) ((epair((i*(i-1))/2+j),i,j,j=1,i),i=1,nocc)
C 902     format(4(f13.8,1x,'(',2i3,')',2x))
       endif
         
       return
       end








       subroutine moints_e2pair( nocc, ostart, occ_lo, occ_hi,
     $                           nvir, nmo2, nbf, g_k, orbe, g_epair )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
       integer nocc, ostart, occ_lo, occ_hi, nvir, nmo2
       double precision orbe(*)
       integer g_k, g_epair
       integer k_exch,ld_exch
       integer ii,i,jj,j,ij,vstart,nnt,nbf,ndim,ij0,ofroz
       integer ilo,ihi,jlo,jhi,my_id
       integer ooff
       double precision xx,denom
       double precision moints_epair_eval
c
c
       ofroz = ostart - 1
       ooff = ((occ_lo-ofroz)*(occ_lo-ofroz-1))/2
       nnt = ((occ_hi-ofroz)*(occ_hi-ofroz+1))/2 - ooff
c
c Check if this node has any data
c
       my_id = ga_nodeid()
       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)
       if ((ilo.ne.1).and.(ihi.ne.-1)) call 
     $ errquit('moints_mp2_eval: wrong distribution for exchange',my_id)
c
c Set to handle different exchange matrix sizes
c Usually K matrix is (nvir,nvir) but could be (nbf,nbf) 
c
       ndim = nvir
       vstart = 0
       if ((ihi.eq.(nbf*nbf)).and.(nmo2.eq.nbf)) then
         ndim = nbf
         vstart = nocc
       else if ((ihi.ne.(nvir*nvir)).and.(ihi.ne.(nvir*nbf)).and.
     $          (ihi.ne.(nbf*nbf))) then
         call errquit('moints_mp2_eval: odd exchange matrix size',ihi)
       endif
c
c Loop over occupied orbital pairs
c
       do i=occ_lo,occ_hi
         ii = i - ofroz
         do j=ostart,i
           jj = j - ofroz
           ij0 = (ii*(ii-1))/2 + jj
           ij = ij0 - ooff
           if ((ij.ge.jlo).and.(ij.le.jhi)) then
             denom = orbe(i) + orbe(j)
             call ga_access(g_k,1,ihi,ij,ij,k_exch,ld_exch)
             xx = moints_epair_eval(ndim,vstart,nocc,nvir,
     $                             dbl_mb(k_exch),orbe,denom)
             call ga_release(g_k,1,ihi,ij,ij)
             if (i.eq.j) xx = xx/2.d0
             call ga_put(g_epair,1,1,ij0,ij0,xx,1)
           endif
         enddo
       enddo
       call ga_sync()

c$$$c
c$$$c *DEBUG*
c$$$c Print out current set of pair energies
c$$$c
c$$$       iii = 0
c$$$       if (ga_nodeid().eq.0) then
c$$$         do i=occ_lo,occ_hi
c$$$  	   ii = i - ofroz
c$$$           do j=ostart,i
c$$$             jj = j - ofroz
c$$$             ij0 = (ii*(ii-1))/2 + jj
c$$$             call ga_get(g_epair,1,1,ij0,ij0,xx,1)
c$$$             iii = iii + 1
c$$$             crap(iii) = xx
c$$$           enddo
c$$$         enddo
c$$$         write(6,883) (crap(jjj),jjj=1,iii)
c$$$ 883     format(8f10.5)
c$$$       endif
         
       return
       end



       double precision function moints_epair_eval( ndim, vstart, nocc,
     $                                          nvir, xk, orbe, denom )
       implicit none
       integer ndim, vstart, nocc, nvir
       double precision xk(ndim,ndim), orbe(nocc+nvir), denom
       double precision xx, yy
       integer k,l,kk,ll

       xx = 0.d0
       do kk=1,nvir
         do ll=1,nvir
           k = vstart + kk
           l = vstart + ll
           yy = xk(l,k)*(4.d0*xk(l,k) - 2.d0*xk(k,l))
           xx = xx + yy/(denom - orbe(nocc+kk) - orbe(nocc+ll))
         enddo
       enddo
       moints_epair_eval = xx
       return
       end
