c
c
c     This routine seems to be dead, and it's call to moints_epair_eval
c     is incorrect.  What else in this file is dead?
c
c
c$$$       subroutine moints_mp2_eval(nocc,nvir,orbe,g_k,emp2,epair)
c$$$C$Id$
c$$$       implicit none
c$$$#include "errquit.fh"
c$$$#include "global.fh"
c$$$#include "mafdecls.fh"
c$$$       integer nocc, nvir
c$$$       double precision orbe(nocc+nvir), emp2, epair(*)
c$$$       integer g_k
c$$$       integer k_exch,ld_exch
c$$$       integer i,j,ij,vstart,nnt,nbf,ndim
c$$$       integer ilo,ihi,jlo,jhi,my_id
c$$$       integer g_epair
c$$$       double precision xx,denom
c$$$       double precision moints_epair_eval
c$$$
c$$$       nnt = (nocc*(nocc+1))/2
*ga:1:0
c$$$       if (.not.ga_create(MT_DBL,nnt,1,'epair',0,0,g_epair)) call
c$$$     $      errquit('moints_mp2_eval: cannot allocate',0)
c$$$       nbf = nocc + nvir
c$$$       my_id = ga_nodeid()
c$$$       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)
c$$$
c$$$#ifdef DEBUG
c$$$       if (ga_nodeid().eq.0) write(6,944)
c$$$       do i=0,ga_nnodes()-1
c$$$         if (i.eq.my_id) then
c$$$           write(6,923) i,ilo,ihi,jlo,jhi
c$$$ 923       format('moints_mp2_eval: distribution','[',i5,']',3x,
c$$$     $            2(2i5,3x))
c$$$           call util_flush(6)
c$$$	 endif
c$$$         call ga_sync()
c$$$       enddo
c$$$       if (ga_nodeid().eq.0) write(6,944)
c$$$ 944   format(/,'===================================================',/)
c$$$#endif
c$$$c
c$$$c Check if this node has any data
c$$$c
c$$$       if ((ilo.ne.1).and.(ihi.ne.-1)) call 
c$$$     $ errquit('moints_mp2_eval: wrong distribution for exchange',my_id)
c$$$c
c$$$c Set to handle different exchange matrix sizes
c$$$c
c$$$       if (ihi.eq.(nbf*nbf)) then
c$$$         ndim = nocc+nvir
c$$$         vstart = nocc
c$$$       else if ((ihi.eq.(nvir*nvir)).or.(ihi.eq.(nvir*nbf))) then
c$$$         ndim = nvir
c$$$         vstart = 0
c$$$       else if (ihi.ne.-1) then
c$$$         call errquit('moints_mp2_eval: odd exchange matrix size',ihi)
c$$$       endif
c$$$c
c$$$c Loop over occupied orbital pairs
c$$$c
c$$$       do i=1,nocc
c$$$         do j=1,i
c$$$           ij = (i*(i-1))/2 + j
c$$$           if ((ij.ge.jlo).and.(ij.le.jhi)) then
c$$$             denom = orbe(i) + orbe(j)
c$$$             call ga_access(g_k,1,ihi,ij,ij,k_exch,ld_exch)
c$$$             xx = moints_epair_eval(ndim,vstart,nvir,
c$$$     $                             dbl_mb(k_exch),orbe(nocc+1),
c$$$     $                             denom)
c$$$             call ga_release(g_k,1,ihi,ij,ij)
c$$$             if (i.eq.j) xx = xx/2.d0
c$$$             call ga_put(g_epair,ij,ij,1,1,xx,1)
c$$$           endif
c$$$         enddo
c$$$       enddo
c$$$c
c$$$c Get global pair-energies into local array
c$$$c
c$$$       call ga_sync()
c$$$       emp2 = 0.d0
c$$$       do ij=1,nnt
c$$$         call ga_get(g_epair,ij,ij,1,1,xx,1)
c$$$         epair(ij) = xx
c$$$         emp2 = emp2 + xx
c$$$       enddo
c$$$c
c$$$c Clean and print
c$$$c
c$$$       if (.not.ga_destroy(g_epair)) call
c$$$     $      errquit('moints_mp2_eval_a: cannot destroy global',0)
c$$$       if (ga_nodeid().eq.0) then
c$$$C         write(6,*)
c$$$C         call util_print_centered(6,'MP2 and Pair Energies',40,.true.)
c$$$         write(6,901) emp2
c$$$ 901     format(//,15x,'MP2 correction:',f15.10,/)
c$$$C         write(6,902) ((epair((i*(i-1))/2+j),i,j,j=1,i),i=1,nocc)
c$$$C 902     format(4(f13.8,1x,'(',2i3,')',2x))
c$$$       endif
c$$$         
c$$$       return
c$$$       end
c$$$





c
c  Evaluates MP2 pair energies from given
c  exchange operator:
c
c            pq               p,q = arbitrary indices
c           K                 i,j = correlated occ. orbitals for this pass
c            ij
c
c  Obviously, the superscript range mo2_lo:mo2_hi must contain the
c  the correlated virtual range, vlo:vhi, viz,
c
c      mo2_lo <= vlo <= vhi <= mo2_hi
c  

       subroutine moints_e2pair(rtdb, nocc, ostart, occ_lo, occ_hi, 
     $                          mo2_lo, mo2_hi, vlo, vhi, g_k, 
     $                          orbe, g_epair)
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
       integer rtdb                                  ! [input]  run-time database handle
       integer nocc                                  ! [input]  number of occupied orbitals (including frozen)
       integer ostart                                ! [input]  first correlated occupied index
       integer occ_lo, occ_hi                        ! [input]  index range of correlated occupied for this pass 
       integer mo2_lo, mo2_hi                        ! [input]  index range of operator superscript
       integer vlo, vhi                              ! [input]  index range of correlated virtuals
       double precision orbe(*)                      ! [input]  orbital energies
       integer g_k                                   ! [input]  exchange operator
       integer g_epair                               ! [output] triangular array of pair energies
c
       integer k_exch,ld_exch
       integer ii,i,jj,j,ij,nnt,ij0,ofroz
       integer ilo,ihi,jlo,jhi,my_id
       integer ooff,nmo2
       double precision xx,denom
       double precision moints_epair_eval
       external moints_epair_eval
       double precision c1,c2       ! constants for moints_epair_eval (because of SCS-MP2)
       double precision fss, fos    ! Spin Component Scaled (SCS) MP2
       logical scs                  ! SCS-MP2 flag
#ifdef BAD_GACCESS
      integer l_local
#endif
c$$$       INTEGER III,JJJ
c$$$       DOUBLE PRECISION CRAP(1000)
c

c     flag and scaling factors for SCS-MP2

       if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $   fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
       if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $   fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
       if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $   scs = .false.       ! flag for SCS-MP2
         if (scs ) then
           c1 = 2.0d0*(fss+fos)
           c2 = 2.0d0*fss
         else
           c1 = 4.0d0
           c2 = 2.0d0
         endif
c
       ofroz = ostart - 1
       ooff = ((occ_lo-ofroz)*(occ_lo-ofroz-1))/2
       nnt = ((occ_hi-ofroz)*(occ_hi-ofroz+1))/2 - ooff
       nmo2 = mo2_hi - mo2_lo + 1
c
c Check if this node has any data
c
       my_id = ga_nodeid()
       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)
       if ((ilo.ne.1).and.(ihi.ne.-1)) call 
     $   errquit('moints_e2pair: wrong distrib for exchange',my_id,
     &       INPUT_ERR)
       if (ihi.ne.-1) then
         if (ihi.lt.(nmo2*nmo2)) call 
     $     errquit('moints_e2pair: exchange dim & param mismatch',0,
     &       INPUT_ERR)
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
#ifdef BAD_GACCESS
             ld_exch=ihi
             if(.not.ma_push_get(MT_DBL,ld_exch,
     $            'scratch buffmp2', l_local, k_exch)) call
     $            errquit('mointsmp2: pushget failed',0,0)
             call ga_get(g_k,1,ihi,ij,ij,dbl_mb(k_exch),ld_exch)
#else
             call ga_access(g_k,1,ihi,ij,ij,k_exch,ld_exch)
#endif
             xx = moints_epair_eval( vlo, vhi, mo2_lo, mo2_hi, 
     $                                 dbl_mb(k_exch), orbe(nocc+1), 
     $                                 denom, c1, c2)
#ifdef BAD_GACCESS
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('mointsmp2: popstack failed',0,0)
#else
             call ga_release(g_k,1,ihi,ij,ij)
#endif
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




       double precision 
     $    function moints_epair_eval(vlo, vhi, mo2_lo, mo2_hi,
     $                                 xk, virorbe, denom, c1, c2 )
       implicit none
#include "errquit.fh"
#include "global.fh"
       integer vlo, vhi, mo2_lo, mo2_hi
       double precision xk(mo2_lo:mo2_hi,mo2_lo:mo2_hi)
       double precision virorbe(vlo:vhi), denom
       double precision xx, yy
       double precision xxxx
       integer a,b
       double precision c1,c2

       xx = 0.d0
       do a=vlo,vhi
         do b=vlo,vhi
           yy = xk(a,b)*(c1*xk(a,b) - c2*xk(b,a))
           xxxx = denom  - virorbe(b) - virorbe(a)
           if (abs(xxxx).lt.1.d-10)
     $       call errquit('moints_epair_eval: zero denominator',0,
     &       INPUT_ERR)
           xx = xx + yy/xxxx
         enddo
       enddo
       moints_epair_eval = xx
       return
       end









