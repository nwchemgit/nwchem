c
c
c  Sigma vector Alpha - Alpha contribution
c
c
      subroutine detci_sigmaaa( norb, nsym, nela, nelb, nstra, nstrb,
     $                          nexa, nexb, nekla, neklb,
     $                          osym, mul, ijmap, exa, exb,
     $                          ataba, atabb, ntij,
     $                          h, g, f, g_civec, g_sigma )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
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
      integer mul(nsym,nsym)                  ! [input] Direct product table
      integer ijmap(norb,norb)                ! [input] Map of (i,j) -> ij (symmtery-blocked)
      integer ataba(norb,nela,nsym)           ! [input] Alpha arc weights
      integer atabb(norb,nelb,nsym)           ! [input] Beta arc weights
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer ntij                            ! [input] Symmtery-blocked triangular sum
      double precision h(ntij)                ! [input] One-electron hamiltonian
      double precision g(ntij,ntij)           ! [input] ERI's
      double precision f(nstrb)               ! [scratch] Scratch space
      integer g_civec                         ! [input] CI-vector          
      integer g_sigma                         ! [input/output] Sigma vector
c
c
c
      integer istr, kstr, jstr
      integer iex, kex
      integer rlo, rhi, cilo, cihi
      integer k_ci, k_sig, ldc, lds
      integer myid
      integer i, j, k, l
      integer kl, ij, kj, lj
      integer kphase, jphase
      integer isym, ksym, jsym, iksym
      integer trii(detci_maxorb)
      double precision xx, yy
      double precision ddot
      external ddot
c
c
c
      do i=1,norb
        trii(i) = (i*(i-1))/2
      enddo
c
c
      if (.not.ga_compare_distr( g_civec, g_sigma ))
     $  call errquit('detci_sigmaaa: CI vectors do not match',0)
      myid = ga_nodeid()
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_sigmaaa: wrong distrib for CI vector',0)
c
c
c
      isym = 1
      do istr=1,nstrb
        if ((isym.lt.nsym).and.(istr.gt.ataba(1,1,isym+1)))
     $    isym = isym + 1
        call dfill(nstrb,0.d0,f,1)
        do iex=1,nexb
          kstr   = exb(1,iex,istr)
          ksym   = exb(2,iex,istr)
          kl     = exb(3,iex,istr)
          kphase = exb(4,iex,istr)
          k      = exb(5,iex,istr)
          l      = exb(6,iex,istr)
          iksym = mul(isym,ksym)
C          if (iksym.eq.1) then
            xx = 0.d0
            do j=1,norb
              kj = ijmap(k,j)
              lj = ijmap(l,j)
              xx = xx + g(kj,lj)
            enddo
            f(kstr) = f(kstr) + kphase*(h(kl) - 0.5d0*xx)
C          endif
c
          do kex=1,nexb
            jstr   = exb(1,kex,kstr)
            jsym   = exb(2,kex,kstr)
            ij     = exb(3,kex,kstr)
            jphase = exb(4,kex,kstr)*kphase
            f(jstr) = f(jstr) + 0.5d0*jphase*g(ij,kl)
          enddo
        enddo
c
c Data parallel here....
c
        do i=1,nstra
          if ((i.ge.cilo).and.(i.le.cihi)) then
            call ga_access( g_civec, rlo, rhi, i, i, k_ci, ldc )
            call ga_access( g_sigma, rlo, rhi, i, i, k_sig, lds )
            xx = ddot(nstrb,f,1,dbl_mb(k_ci),1)
            dbl_mb(k_sig+istr-1) = dbl_mb(k_sig+istr-1) + xx
            call ga_release( g_civec, rlo, rhi, i, i )
            call ga_release_update( g_sigma, rlo, rhi, i, i )
          endif
        enddo
c
c            
      enddo
      return
      end






c
c
c  Sigma vector Beta-Beta contribution
c
c
      subroutine detci_sigmabb( norb, nsym, nela, nelb, nstra, nstrb,
     $                          nexa, nexb, nekla, neklb,
     $                          osym, mul, ijmap, exa, exb,
     $                          ataba, atabb, ntij,
     $                          h, g, f, g_civec, g_sigma )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
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
      integer mul(nsym,nsym)                  ! [input] Direct product table
      integer ijmap(norb,norb)                ! [input] Map of (i,j) -> ij (symmtery-blocked)
      integer ataba(norb,nela,nsym)           ! [input] Alpha arc weights
      integer atabb(norb,nelb,nsym)           ! [input] Beta arc weights
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer ntij                            ! [input] Symmtery-blocked triangular sum
      double precision h(ntij)                ! [input] One-electron hamiltonian
      double precision g(ntij,ntij)           ! [input] ERI's
      double precision f(nstra)               ! [scratch] Scratch space
      integer g_civec                         ! [input] CI-vector          
      integer g_sigma                         ! [input/output] Sigma vector
c
c
c
      integer istr, kstr, jstr
      integer iex, kex
      integer rlo, rhi, cilo, cihi
      integer k_ci, k_sig, ldc, lds
      integer myid
      integer i, j, k, l
      integer kl, ij, kj, lj
      integer kphase, jphase
      integer isym, ksym, jsym, iksym
      integer trii(detci_maxorb)
      double precision xx,yy
      double precision tmp(10000)
c
      if (10000 .lt. nstrb) 
     $     call errquit('sigmabb: hard dimension failed',nstrb)
c
c
c
      do i=1,norb
        trii(i) = (i*(i-1))/2
      enddo
c
c
      if (.not.ga_compare_distr( g_civec, g_sigma ))
     $  call errquit('detci_sigmaaa: CI vectors do not match',0)
      myid = ga_nodeid()
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_sigmabb: wrong distrib for CI vector',0)
c
c
c
      isym = 1
      do istr=1,nstra
        if ((isym.lt.nsym).and.(istr.gt.ataba(1,1,isym+1)))
     $    isym = isym + 1
        call dfill(nstra,0.d0,f,1)
        do iex=1,nexa
          kstr   = exa(1,iex,istr)
          ksym   = exa(2,iex,istr)
          kl     = exa(3,iex,istr)
          kphase = exa(4,iex,istr)
          k      = exa(5,iex,istr)
          l      = exa(6,iex,istr)
          iksym = mul(isym,ksym)
C          if (iksym.eq.1) then
            xx = 0.d0
            do j=1,norb
              kj = ijmap(k,j)
              lj = ijmap(l,j)
              xx = xx + g(kj,lj)
            enddo
            f(kstr) = f(kstr) + kphase*(h(kl) - 0.5d0*xx)
C          endif
c
c
          do kex=1,nexa
            jstr   = exa(1,kex,kstr)
            jsym   = exa(2,kex,kstr)
            ij     = exa(3,kex,kstr)
            jphase = exa(4,kex,kstr)*kphase
            i      = exa(5,kex,kstr)
            j      = exa(6,kex,kstr)
            f(jstr) = f(jstr) + 0.5d0*jphase*g(ij,kl)
          enddo
        enddo
c
c  Parallel section here...
c
c$$$        do i=1,nstra
c$$$          if ((i.ge.cilo).and.(i.le.cihi)) then
c$$$            call ga_access( g_civec, rlo, rhi, i, i, k_ci, ldc )
c$$$            call ga_access( g_sigma, rlo, rhi, i, i, k_sig, lds )
c$$$            xx = ddot(nstrb,f,1,dbl_mb(k_ci),1)
c$$$            dbl_mb(k_sig+istr-1) = dbl_mb(k_sig+istr-1) + xx
c$$$            call ga_release( g_civec, rlo, rhi, i, i )
c$$$            call ga_release_update( g_sigma, rlo, rhi, i, i )
c$$$          endif
c$$$        enddo

        do j=1,nstra
           if (f(j).ne.0.0d0 .and. (j.ge.cilo).and.(j.le.cihi)) then
              call ga_access( g_civec, rlo, rhi, j, j, k_ci, ldc )
              do i=1,nstrb
c$$$  xx = 0.d0
c$$$  call ga_get(g_civec,i,i,j,j,yy,1)
                 tmp(i) = f(j)*dbl_mb(k_ci+i-1)
              enddo
              call ga_acc(g_sigma,1,nstrb,istr,istr,tmp,1,1.d0)
              call ga_release( g_civec, rlo, rhi, j, j )
           endif
        enddo
c
c
c            
      enddo
      return
      end









c
c
c  Sigma vector Alpha-Beta contribution 
c
c
      subroutine detci_sigmaab( norb, nsym, nela, nelb, nstra, nstrb,
     $                          nexa, nexb, nekla, neklb,
     $                          osym, mul, ijmap, exa, exb,
     $                          ataba, atabb, ntij, g,
     $                          vrhs, vlhs, vphase, f,
     $                          cprime, sprime, g_civec, g_sigma )
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
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
      integer mul(nsym,nsym)                  ! [input] Direct product table
      integer ijmap(norb,norb)                ! [input] Map of (i,j) -> ij (symm-blocked)
      integer ataba(norb,nela,nsym)           ! [input] Alpha arc weights
      integer atabb(norb,nelb,nsym)           ! [input] Beta arc weights
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer ntij                            ! [input] Symmtery-blocked triangular sum
      double precision g(ntij,ntij)           ! [input] ERI's
      integer vrhs(nekla)                     ! [scratch] Array of RHS Strings for E_kl
      integer vlhs(nekla)                     ! [scratch] Array of LHS Strings for E_kl
      integer vphase(nekla)                   ! [scratch] Array of Phases for E_kl
      double precision f(nstrb)               ! [scratch] Scratch space
      double precision cprime(nstrb,nekla)    ! [scratch] Gathered CI-vector
      double precision sprime(nstrb,nekla)    ! [scratch] Gathered sigma vector
      integer g_civec                         ! [input] CI-vector          
      integer g_sigma                         ! [input/output] Sigma vector
c$$$      double precision civec(nstrb,nstra)     ! [input] CI-vector
c$$$      double precision sigma(nstrb,nstra)     ! [input/output] Sigma vector
c
c
      integer i, ij, kl
      integer ak, al, iph
      integer ksym
      integer jstr, kstr
      integer ne_kl, nidx, iex
      integer relv(detci_maxelec)
      integer lelv(detci_maxelec)
      integer oidx(detci_maxorb)
      integer ip(detci_maxelec)
      integer ploop, numnodes, next, myid
      double precision phase
      integer nxtval
      external nxtval
c
c$$$      INTEGER IREV
c$$$      INTEGER REV_RHS(1000)
c
c Initialize parallel stuff
c
      ploop = -1
      numnodes = ga_nnodes()
      myid = ga_nodeid()
      next = nxtval(numnodes)
c
c
c Loop over all excitation operators
c
c            t
c      E  = a a
c       kl   k l
c
      do ak=1,norb
        do al=1,norb
          ploop = ploop + 1
          if (ploop.eq.next) then
            kl = ijmap(ak,al)
            ne_kl = 0
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
c Insert orbital index k and l to create 
c LHS and RHS strings where
c
c        |LHS> = E  |RHS>
c                 kl

c$$$          CALL IFILL(NSTRA,0,REV_RHS,1)
          
            if (nela.gt.1) then
 101          continue
              ne_kl = ne_kl + 1
              iph = 1
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           al,relv,iph)
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           ak,lelv,iph)
c
c Put index into gather/scatter arrays
c

              vphase(ne_kl) = iph
              vrhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,mul,
     $                                    ataba,relv)
              vlhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,mul,
     $                                    ataba,lelv)

c$$$          IREV =  DETCI_ELV2STR(NORB,NELA,NSYM,OSYM,MUL,ATABA,RELV)
c$$$          REV_RHS(IREV) = NE_KL

              if (detci_getnextelv(nidx,(nela-1),ip)) goto 101
            else
              ne_kl = 1
              vrhs(i) = al
              vlhs(i) = ak
              vphase(i) = 1
            endif

C          WRITE(6,'(20I4)') (VRHS(I),I=1,NE_KL)
C          WRITE(6,'(20I4)') (REV_RHS(I),I=1,NSTRA)

c
c End loop over possible strings
c
c Gather in CI and sigma blocks
c
            call detci_cigather( nstrb, nstra, ne_kl, g_civec, vlhs,
     $                       vphase, cprime )
            call dfill((nstrb*ne_kl),0.d0,sprime,1)
c
c Loop over all beta strings
c
            do jstr=1,nstrb
              call dfill(nstrb,0.d0,f,1)
              do iex=1,nexb
                kstr  = exb(1,iex,jstr)
                ksym  = exb(2,iex,jstr)
                ij    = exb(3,iex,jstr)
                phase = exb(4,iex,jstr)
                f(kstr) = f(kstr) + phase*g(ij,kl)

c$$$              I     = EXB(5,IEX,JSTR)
c$$$              J     = EXB(6,IEX,JSTR)
c$$$              IF ((JSTR.EQ.1).AND.(G(IJ,KL).GT.0.D0))
c$$$     $            WRITE(6,'(A,4I5,2F12.6)') '===',
c$$$     $            AK,AL,I,J,(PHASE*G(IJ,KL)),F(KSTR)

              enddo
              call detci_dotab(jstr,ne_kl,nstrb,f,cprime,sprime)

c$$$            IF ((JSTR.EQ.1).AND.(REV_RHS(1).NE.0))
c$$$     $            WRITE(6,'(A,20X,2F12.6)')
c$$$     $            '>>>',SPRIME(1,REV_RHS(1))
c$$$
            enddo
c
c Scatter accumulate result into sigma vector
c
            call detci_ciscatter(nstrb,nstra,ne_kl,sprime,vrhs,g_sigma)
c
c End parallel task
c
            next = nxtval(numnodes)
          endif
        enddo
      enddo
      next = nxtval(-numnodes)
      return
      end






























