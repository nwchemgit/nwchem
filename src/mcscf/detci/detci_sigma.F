
c
c
c  Sigma vector Alpha - Alpha contribution
c
c
      subroutine detci_sigmaaa( norb, nsym, nela, nelb, nstra, nstrb,
     $                          nexa, nexb, nekla, neklb,
     $                          osym, ijmap, exa, exb,
     $                          ataba, atabb, ntij, nsblk,
     $                          h, g, f, g_civec, g_sigma )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetcistats.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer nexa                            ! [input] Alpha excitations
      integer nexb                            ! [input] Beta excitations
      integer nekla                           ! [input] Maximum non-zero alpha strings for E_kl
      integer neklb                           ! [input] Maximum non-zero beta strings for E_kl
      integer osym(norb)                      ! [input] Orbital irreps     
      integer ijmap(norb,norb)                ! [input] Map of (i,j) -> ij (symmtery-blocked)
      integer ataba(norb,nela,nsym)           ! [input] Alpha arc weights
      integer atabb(norb,nelb,nsym)           ! [input] Beta arc weights
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer ntij                            ! [input] Symmtery-blocked triangular sum
      integer nsblk                           ! [input] Blocking factor to increase parallelism
      double precision h(ntij)                ! [input] One-electron hamiltonian
      double precision g(ntij,ntij)           ! [input] ERI's
      double precision f(nstrb,nsblk)         ! [scratch] Scratch space
      integer g_civec                         ! [input] CI-vector          
      integer g_sigma                         ! [input/output] Sigma vector
c
      integer Msg_detci_sum
      parameter(Msg_detci_sum=1059)
c
c
c
      integer istr, kstr, jstr
      integer istrlo, istrhi, iistr
      integer msblk, iscnt
      integer iex, kex
      integer rlo, rhi, cilo, cihi
      integer k_ci, k_sig, ldc, lds
      integer myid, numnodes
      integer i, j, k, l
      integer kl, ij, kj, lj
      integer kphase, jphase
      integer isym, ksym, jsym, iksym
      double precision xx, yy, tx, tff, tdot, tgop
      double precision h1(detci_maxorb*detci_maxorb)
c$$$      double precision sdot
c$$$      external sdot

c
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"
c
c
      tff = 0.d0
      tdot = 0.d0
      tgop = 0.d0
      if (.not.ga_compare_distr( g_civec, g_sigma ))
     $  call errquit('detci_sigmaaa: CI vectors do not match',0)
      myid = ga_nodeid()
      numnodes = ga_nnodes()
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_sigmaaa: wrong distrib for CI vector',0)
c
c   Precompute 2e
c
      do k=1,norb
        do l=1,k
          xx = 0.d0
          do j=1,norb
            kj = ijmap(k,j)
            lj = ijmap(l,j)
            xx = xx + g(kj,lj)
          enddo
          kl = ijmap(k,l)
          h1(kl) = xx*0.5d0
        enddo
      enddo
c
c
c
      do istrlo=1,nstrb,nsblk
        istrhi = min((istrlo + nsblk - 1),nstrb)
        msblk = istrhi - istrlo + 1
        tx = util_cpusec()
        call dfill((msblk*nstrb),0.d0,f,1)
        do istr=istrlo,istrhi
          iistr = istr - istrlo + 1
          iscnt = (iistr - 1)*nexb
          do iex=1,nexb
            if (mod((iscnt+iex-1),numnodes).eq.myid) then
              kstr   = exb(1,iex,istr)
              kl     = exb(3,iex,istr)
              kphase = exb(4,iex,istr)
              f(kstr,iistr) = f(kstr,iistr) + kphase*(h(kl) - h1(kl))
c     
              do kex=1,nexb
                jstr   = exb(1,kex,kstr)
                ij     = exb(3,kex,kstr)
                jphase = exb(4,kex,kstr)*kphase
                f(jstr,iistr) = f(jstr,iistr) + 0.5d0*jphase*g(ij,kl)
              enddo
            endif
          enddo
        enddo
        tff = tff + util_cpusec() - tx
        tx = util_cpusec()
        call ga_dgop(Msg_detci_sum, f, (nstrb*msblk), '+' )
        tgop = tgop + util_cpusec() - tx
c
c Data parallel here....
c
        tx = util_cpusec()
        do i=1,nstra
          if ((i.ge.cilo).and.(i.le.cihi)) then
            call ga_access( g_civec, rlo, rhi, i, i, k_ci, ldc )
            call ga_access( g_sigma, rlo, rhi, i, i, k_sig, lds )
            do istr=istrlo,istrhi
              iistr = istr - istrlo + 1
              xx = ddot( nstrb, f(1,iistr), 1, dbl_mb(k_ci), 1 )
              dbl_mb(k_sig+istr-1) = dbl_mb(k_sig+istr-1) + xx
            enddo
            call ga_release( g_civec, rlo, rhi, i, i )
            call ga_release_update( g_sigma, rlo, rhi, i, i )
          endif
        enddo
        tdot = tdot + util_cpusec() - tx
      enddo
c
c
c  Collect stats
c
      detci_aaff_etime = detci_aaff_etime + tff
      detci_aadot_etime = detci_aadot_etime + tdot
      detci_aagop_etime = detci_aagop_etime + tgop
      return
      end












c
c
c  Sigma vector Alpha-Beta contribution 
c
c
      subroutine detci_sigmaab( norb, nsym, nela, nelb, nstra, nstrb,
     $                          nexa, nexb, nekla, neklb,
     $                          osym, ijmap, exa, exb,
     $                          ataba, atabb, ntij, g,
     $                          vrhs, vlhs, vphase, f,
     $                          cprime, sprime, g_civec, g_sigma )
      implicit none
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetcistats.fh"
      integer norb                            ! [input] Orbitals
      integer nsym                            ! [input] Irreps
      integer nela                            ! [input] Alpha electrons
      integer nelb                            ! [input] Beta electrons
      integer nstra                           ! [input] Alpha strings
      integer nstrb                           ! [input] Beta strings
      integer nexa                            ! [input] Alpha excitations
      integer nexb                            ! [input] Beta excitations
      integer nekla                           ! [input] Max non-zero alpha strings for E_kl
      integer neklb                           ! [input] Max non-zero beta strings for E_kl
      integer osym(norb)                      ! [input] Orbital irreps     
      integer ijmap(norb,norb)                ! [input] Map of (i,j) -> ij (symm-blocked)
      integer ataba(norb,nela,nsym)           ! [input] Alpha arc weights
      integer atabb(norb,nelb,nsym)           ! [input] Beta arc weights
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer ntij                            ! [input] Symmtery-blocked triangular sum
      double precision g(ntij,ntij)           ! [input] ERI's
      integer vrhs(nekla)                     ! [scratch] Array of RHS Strings for E_kl
      integer vlhs(*)                         ! [scratch] Array of LHS Strings for E_kl
      integer vphase(nekla)                   ! [scratch] Array of Phases for E_kl
      double precision f(nexb)                ! [scratch] Scratch space
      double precision cprime(nstrb,nekla)    ! [scratch] Gathered CI-vector
      double precision sprime(nstrb,nekla)    ! [scratch] Gathered sigma vector
      integer g_civec                         ! [input] CI-vector          
      integer g_sigma                         ! [input/output] Sigma vector
c
c
      integer i, ij, kl, ii, iii
      integer ak, al, iph
      integer ksym
      integer jstr, kstr
      integer ne_kl, nidx, iex
      integer relv(detci_maxelec)
      integer lelv(detci_maxelec)
      integer oidx(detci_maxorb)
      integer ip(detci_maxelec)
      integer ploop, numnodes, next, myid
      integer nsblk, isblk, slo, shi, sseg
      double precision phase
      double precision tstr, tgath, tdotab, tscat, tsync, tx
      integer nxtval
      external nxtval
c
c
c  Initialize parallel stuff
c
      tstr = 0.d0
      tgath = 0.d0
      tdotab = 0.d0
      tscat = 0.d0
      ploop = -1
      numnodes = ga_nnodes()
      myid = ga_nodeid()
c
c  Block over strings for finer granularity
c
      nsblk = (20*numnodes)/(norb*norb)
      nsblk = max(nsblk,1)
      NSBLK = 2
      if (numnodes.eq.1) nsblk = 1
      sseg = (nekla/nsblk)
      if (mod(nekla,nsblk).ne.0) sseg = sseg + 1
c
c  Loop over all excitation operators
c
c            t
c      E  = a a
c       kl   k l
c
      next = nxtval(numnodes)
      do ak=1,norb
        do al=1,norb
          shi = 0
          do isblk=1,nsblk
            slo = shi + 1
            shi = min((slo + sseg - 1),nekla)
            ploop = ploop + 1
            if (ploop.eq.next) then
              tx = util_cpusec()
              kl = ijmap(ak,al)
c
c Vector of orbital indices except create/annih indices
c Initialize pointer vector
c
              nidx = 0
              do i=1,norb
                if ((i.ne.ak).and.(i.ne.al)) then
                  nidx = nidx + 1
                  oidx(nidx) = i
                endif
              enddo
              do i=1,nela-1
                ip(i) = i
              enddo
c
c Loop through all strings for nidx and (nela-1)
c Accept strings in the block range slo:shi
c Insert orbital index k and l to create 
c LHS and RHS strings where
c
c        |LHS> = E  |RHS>
c                 kl
c
c Push indices into gather/scatter arrays
c
c     
              if (nela.gt.1) then
                ne_kl = 0
                iii = 0
 101            continue
                iii = iii + 1
                if ((iii.ge.slo).and.(iii.le.shi)) then
                  ne_kl = ne_kl + 1
                  iph = 1
                  call detci_ptr2elv( norb, nela, (nela-1), nidx, ip,
     $                                oidx, al, relv, iph )
                  call detci_ptr2elv( norb, nela, (nela-1), nidx, ip,
     $                                oidx, ak, lelv, iph )
                  vphase(ne_kl) = iph
                  vrhs(ne_kl) = detci_elv2str( norb, nela, nsym, osym,
     $                                         ataba, relv)
                  vlhs(ne_kl) = detci_elv2str( norb, nela, nsym, osym,
     $                                         ataba, lelv)
                endif
                if (detci_getnextelv(nidx,(nela-1),ip)) goto 101
              else
                ne_kl = 1
                vrhs(i) = al
                vlhs(i) = ak
                vphase(i) = 1
              endif
              tstr = tstr + util_cpusec() - tx
              tx = util_cpusec()
c
c End loop over possible strings
c
c Gather in CI blocks
c
              call detci_cigather( nstrb, nstra, ne_kl, g_civec, vlhs,
     $                             vphase, sprime )
              call dfill((nstrb*ne_kl), 0.d0, cprime, 1 )
              call transpose( nstrb, ne_kl, sprime, cprime )
              call dfill((nstrb*ne_kl),0.d0,sprime,1)
              tgath = tgath + util_cpusec() - tx
              tx = util_cpusec()
c
c Loop over all beta strings
c
              do jstr=1,nstrb
                call dfill( nstrb, 0.d0, f, 1 )
                do iex=1,nexb
                  vlhs(iex) = exb(1,iex,jstr)
                  ij        = exb(3,iex,jstr)
                  phase     = exb(4,iex,jstr)
                  f(iex)    = phase*g(ij,kl)
                enddo
                call detci_dotabx(jstr, ne_kl, nstrb, nexb, f, vlhs,
     $                            cprime, sprime )
              enddo
              tdotab = tdotab + util_cpusec() - tx
              tx = util_cpusec()
c
c Scatter accumulate result into sigma vector
c
              call transpose( ne_kl, nstrb, sprime, cprime )
              call detci_ciscatter( nstrb, nstra, ne_kl, cprime,
     $                              vrhs, g_sigma )
              tscat = tscat + util_cpusec() - tx
c
c End parallel task
c
              next = nxtval(numnodes)
            endif
          enddo
        enddo
      enddo
      tx = util_cpusec()
      next = nxtval(-numnodes)
      tsync = util_cpusec() - tx
      detci_abstr_etime   = detci_abstr_etime + tstr
      detci_abgath_etime  = detci_abgath_etime + tgath
      detci_abdotab_etime = detci_abdotab_etime + tdotab
      detci_abscat_etime  = detci_abscat_etime + tscat
      detci_absync_etime  = detci_absync_etime + tsync
      return
      end
