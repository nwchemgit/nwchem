c
c
c This is the interface to use if the exchange operator
c matrices are (nvir*nvir) size. I.e. only unique integrals.
c
c
       subroutine moints_mp2_eval(nbf,nocc,nvir,orbe,
     $                              g_k,emp2,epair)
       implicit none
#include "global.fh"
#include "mafdecls.fh"
       integer nbf, nocc, nvir
       double precision orbe(nbf), emp2, epair(*)
       integer g_k
       integer k_exch, ld_exch
       integer i,j,ij,ilo,ihi,jlo,jhi,n2n2
       integer my_id
       double precision xx,denom
       double precision moint_epair_eval

       n2n2 = nvir*nvir
       my_id = ga_nodeid()
       call ga_distribution(g_k,my_id,ilo,ihi,jlo,jhi)
       if ((ilo.ne.1).or.(ihi.ne.n2n2)) call
     $      errquit('moints_mp2_eval_a: distrib. mismatch',0)
       do i=1,nocc
         do j=1,i
           ij = (i*(i-1))/2 + j
           if ((ij.ge.jlo).and.(ij.le.jhi)) then
             call ga_access(g_k,1,n2n2,ij,ij,k_exch,ld_exch)
             denom = orbe(i) + orbe(j)
             xx = moint_epair_eval(nocc,nvir,dbl_mb(k_exch),
     $                               orbe,denom)
             call ga_release(g_k,1,n2n2,ij,ij)
             if (i.eq.j) xx = xx/2.d0
             emp2 = emp2 + xx
             epair(ij) = xx
           endif
         enddo
       enddo

       call util_print_centered(6,'MP2 and Pair Energies',40,.true.)
       write(6,901) emp2
 901   format(//,10x,'MP2 correction:',f15.10,/)
       write(6,902) ((epair((i*(i-1))/2+j),i,j,j=1,i),i=1,nocc)
 902   format(4(f10.4,1x,'[',2i3,']',2x))
         
       return
       end







c
c Covenience routine to access the exchange matrix easily.
c
       double precision function moint_epair_eval( nocc, nvir, 
     $                                             xk, orbe, denom )
       implicit none
       integer nocc,nvir
       double precision xk(nvir,nvir)
       double precision orbe(nocc+nvir), denom
       double precision xx, yy
       integer k,l

       xx = 0.d0
       do k=1,nvir
         do l=1,nvir
           yy = xk(l,k)*(4.d0*xk(l,k) - 2.d0*xk(k,l))
           xx = xx + yy/(denom - orbe(nocc+k) - orbe(nocc+l))
         enddo
       enddo
       moint_epair_eval = xx
       return
       end










c======================================================================
c
c These routines should never have to be used (only for debugging)
c
c======================================================================
c
c This is the interface to use if the exchange operator
c matrices are (nbf*nbf) size.
c
c Should really use ga_access to use the global K-operator
c matrices directly
c
       subroutine moints_mp2_eval_a(nbf,noper,orbe,g_k,emp2,epair)
       implicit none
#include "global.fh"
#include "mafdecls.fh"
       integer nbf, noper
       double precision orbe(nbf), emp2, epair(*)
       integer g_k
       integer l_scr,k_scr
       integer i,j,ij,vstart
       double precision xx,denom
       logical status
       double precision moint_epair_eval_a

       status = ma_push_get(MT_DBL,(nbf*nbf),'N2 scratch',l_scr,k_scr)
       if (.not.(status)) call errquit('mp2: cannot allocate',0)
       vstart = noper + 1
       do i=1,noper
         do j=1,i
           ij = (i*(i-1))/2 + j
           denom = orbe(i) + orbe(j)
           call ga_get(g_k,1,nbf*nbf,ij,ij,dbl_mb(k_scr),nbf)
           xx = moint_epair_eval_a(nbf,vstart,dbl_mb(k_scr),orbe,denom)
           if (i.eq.j) xx = xx/2.d0
           emp2 = emp2 + xx
           epair(ij) = xx
         enddo
       enddo
       status = ma_pop_stack(l_scr)
       if (.not.(status)) call errquit('mp2: pop stack failed',0)

       call util_print_centered(6,'MP2 and Pair Energies',40,.true.)
       write(6,901) emp2
 901   format(//,10x,'MP2 correction:',f15.10,/)
       write(6,902) ((epair((i*(i-1))/2+j),i,j,j=1,i),i=1,noper)
 902   format(4(f10.4,1x,'[',2i3,']',2x))
         
       return
       end





       double precision function moint_epair_eval_a( nbf, vstart,
     $                                             xk, orbe, denom )
       implicit none
       integer nbf,vstart
       double precision xk(nbf,nbf)
       double precision orbe(nbf), denom
       double precision xx, yy
       integer k,l

       xx = 0.d0
       do k=vstart,nbf
         do l=vstart,nbf
           yy = xk(l,k)*(4.d0*xk(l,k) - 2.d0*xk(k,l))
           xx = xx + yy/(denom - orbe(k) - orbe(l))
         enddo
       enddo
       moint_epair_eval_a = xx
       return
       end



