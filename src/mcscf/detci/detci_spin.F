c
c  CI Spin-adaption routine
c  Project off contaminants by Lowdin projection operator
c  S quantum number is determined by the Sz component
c  (S = M = (nela - nelb)/2)
c
c
      subroutine detci_spadpt( norb, nsym, nela, nelb, nstra, nstrb,
     $                         osym, ataba, atabb, nexa, nexb,
     $                         exa, exb, g_civec, g_pvec )
*
* $Id$
*
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nela, nelb
      integer nstra, nstrb
      integer osym(norb)
      integer ataba(0:norb,0:nela,nsym)
      integer atabb(0:norb,0:nelb,nsym)
      integer nexa
      integer nexb
      integer exa(6,nexa,nstra)
      integer exb(6,nexb,nstrb)
      integer g_civec
      integer g_pvec
c$$$      double precision civec(nstrb,nstra)
c$$$      double precision pvec(nstrb,nstra)
c
      integer nel, ssmax
      integer ss, tt
      logical nextodd
c
c  Maximum open shell electrons for
c  nel electrons in norb orbitals
c  ==> maximum S quantum number that can be 
c  obtained from these determinants
c
      nel = nela + nelb
      ssmax = 2*min(nel,norb) - nel
c
c  Desired spin state is obtained from Sz 
c        
      ss = nela - nelb
c
c  Now loop over all spin quantum numbers, T.ne.S
c 
      nextodd = .true.
      do tt=ssmax,0,-2
        if (tt.ne.ss) then
          if (nextodd) then
            nextodd = .false.
            call detci_lowdinproj( norb, nsym, nela, nelb, nstra, nstrb,
     $                             osym, ataba, atabb, nexa, nexb,
     $                             exa, exb, ss, tt, g_civec, g_pvec )
          else
            nextodd = .true.
            call detci_lowdinproj( norb, nsym, nela, nelb, nstra, nstrb,
     $                             osym, ataba, atabb, nexa, nexb,
     $                             exa, exb, ss, tt, g_pvec, g_civec )
c$$$          PRINT*,' TT = ',TT
c$$$          CALL DETCI_PRINT( G_CIVEC, 1.D-3)
          endif
        endif
      if (.not.nextodd) call ga_copy( g_pvec, g_civec )
      enddo
      return
      end







      subroutine detci_lowdinproj( norb, nsym, nela, nelb, nstra, nstrb,
     $                             osym, ataba, atabb, nexa, nexb,
     $                             exa, exb, k2, s2, g_civec, g_pvec )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nela, nelb
      integer nstra, nstrb
      integer osym(norb)
      integer ataba(0:norb,0:nela,nsym)
      integer atabb(0:norb,0:nelb,nsym)
      integer nexa
      integer nexb
      integer exa(6,nexa,nstra)
      integer exb(6,nexb,nstrb)
      integer k2
      integer s2
      integer g_civec
      integer g_pvec
c$$$      double precision civec(nstrb,nstra)
c$$$      double precision pvec(nstrb,nstra)
c
c
c
      integer i, istra, jstra, istrb, jstrb, iex
      integer iib, jjb
      integer eij(detci_maxorb,detci_maxorb)
      integer pij(detci_maxorb,detci_maxorb)
      integer occa(detci_maxorb),occb(detci_maxorb)
      integer elva(detci_maxelec), elvb(detci_maxelec)
      integer alo, ahi, blo, bhi, l_xxci, k_xxci
      integer next, node, numnode
      double precision fac, s, k, xxci, xx, fc
c
      s = 0.5d0*s2
      k = 0.5d0*k2
      fac = 0.25d0*(nela - nelb)**2 + 0.5d0*(nela + nelb) - s*(s+1)
      fc  = 1.d0/(k*(k+1) - s*(s+1))
      fac = fac*fc
      node = ga_nodeid()
      call ga_zero(g_pvec)
c
c
c  Lowdin projection operator:
c
c          2
c         S  - s(s + 1)
c      -------------------
c      k(k + 1) - s(s + 1)
c
c
c       2
c  Use S  = S S   + S (S  + 1)
c            + -     z  z
c
c           1        2   1            b  a
c         = -(n - n )  + -(n + n ) - E  E
c           4  a   b     2  a   b     ij ji
c
c
c                                b  a
c  Apply permutation operator  -E  E    for all i,j open shell
c                                ij ji                
c  Loop over CI vector
c
      call ga_distribution( g_civec, node, blo, bhi, alo, ahi )
      if (.not.ma_push_get(MT_DBL, (bhi-blo+1), 'detci:lowdin',
     $                     l_xxci, k_xxci))
     $    call errquit('detci: cannot allocate lowdin',0, MA_ERR)
      do istra = alo, ahi
         call ifill((detci_maxorb*detci_maxorb),0,eij,1)
         call ifill((detci_maxorb*detci_maxorb),0,pij,1)
         do iex=1,nexa
           eij(exa(6,iex,istra),exa(5,iex,istra)) = exa(1,iex,istra)
           pij(exa(6,iex,istra),exa(5,iex,istra)) = exa(4,iex,istra)
         enddo
         call ga_get(g_civec,blo,bhi,istra,istra,dbl_mb(k_xxci),1)
         call ga_acc(g_pvec,blo,bhi,istra,istra,dbl_mb(k_xxci),1,fac)
         do istrb = blo, bhi
            do iex=1,nexb
              iib = exb(5,iex,istrb)
              jjb = exb(6,iex,istrb)
              if (eij(iib,jjb).ne.0) then
                jstrb = exb(1,iex,istrb)
                jstra = eij(iib,jjb)
                xx = -dbl_mb(k_xxci+istrb-blo)*pij(iib,jjb)
     &               *exb(4,iex,istrb)*fc
                if (dabs(xx).gt.1.0d-12)
     &             call ga_acc(g_pvec,jstrb,jstrb,jstra,jstra,xx,1,1.d0)
              endif
            enddo
         enddo
      enddo
      if (.not.ma_pop_stack(l_xxci))
     $   call errquit('cannot pop stack detci:lowdin',0, MA_ERR)
c
c  End loop over CI vector
c
      return
      end






c
c           2
c  Operate S  on CI vector
c
      subroutine detci_ssq( norb, nsym, nela, nelb, nstra, nstrb,
     $                      osym, ataba, atabb, nexa, nexb,
     $                      exa, exb, civec, pvec )
      implicit none
#include "detciP.fh"
#include "detci.fh"
      integer norb
      integer nsym
      integer nela, nelb
      integer nstra, nstrb
      integer osym(norb)
      integer ataba(0:norb,0:nela,nsym)
      integer atabb(0:norb,0:nelb,nsym)
      integer nexa
      integer nexb
      integer exa(6,nexa,nstra)
      integer exb(6,nexb,nstrb)
      double precision civec(nstrb,nstra)
      double precision pvec(nstrb,nstra)
c
c
c
      integer i, istra, jstra, istrb, jstrb, iex
      integer iib, jjb
      integer eij(detci_maxorb,detci_maxorb)
      integer pij(detci_maxorb,detci_maxorb)
c$$$      integer occa(detci_maxorb),occb(detci_maxorb)
      integer elva(detci_maxelec), elvb(detci_maxelec)
      integer ict
      double precision fac
c
c
c
      ict = 0
      call dfill((nstra*nstrb),0.d0,pvec,1)
c
c                                b  a
c  Apply permutation operator  -E  E    all i,j open shell
c                                ij ji                
c  loop over CI vector
c
      do i=1,nela
        elva(i) = i
      enddo
 101  continue
      istra = detci_elv2str( norb, nela, nsym, osym, ataba, elva )
c$$$      call detci_elv2occ( norb, nelb, elva, occa )
      call ifill((detci_maxorb*detci_maxorb),0,eij,1)
      call ifill((detci_maxorb*detci_maxorb),0,pij,1)
      do iex=1,nexa
        eij(exa(6,iex,istra),exa(5,iex,istra)) = exa(1,iex,istra)
        pij(exa(6,iex,istra),exa(5,iex,istra)) = exa(4,iex,istra)
      enddo

      do i=1,nelb
        elvb(i) = i
      enddo
 102  continue
      istrb = detci_elv2str( norb, nelb, nsym, osym, atabb, elvb )
c$$$      call detci_elv2occ( norb, nelb, elvb, occb )
c
c
c
      do iex=1,nexb
        iib = exb(5,iex,istrb)
        jjb = exb(6,iex,istrb)
        if (eij(iib,jjb).ne.0) then
          jstrb = exb(1,iex,istrb)
          jstra = eij(iib,jjb)
          pvec(jstrb,jstra) = pvec(jstrb,jstra) -
     $                civec(istrb,istra)*pij(iib,jjb)*exb(4,iex,istrb)
        endif
      enddo
c
c
c
      if (detci_getnextelv( norb, nelb, elvb )) goto 102
      if (detci_getnextelv( norb, nela, elva )) goto 101
c
c End loop over CI vector
c                            1         2    1
c Now do trivial components  -(n  - n )  +  -(n - n )
c                            4  a    b      2  a   b
c
      fac = 0.25d0*(nela - nelb)**2 + 0.5d0*(nela + nelb)
      call daxpy((nstra*nstrb),fac,civec,1,pvec,1)
c
c
c
      return
      end







c
c  Insert random errors into g_a
c  to simulate spin-contamination
c
      subroutine detci_random_error( g_a )
      implicit none
#include "global.fh"
      integer g_a
      double precision xx, yy
      integer type, nrow, ncol, i, j
      double precision util_random
      external util_random

      xx = sqrt(ga_ddot(g_a,g_a))*0.001d0
      call ga_inquire(g_a,type,nrow,ncol)
      do i=1,nrow
        do j=1,ncol
          yy = util_random(0)*xx
          call ga_acc(g_a,i,i,j,j,yy,1,1.d0)
        enddo
      enddo
      return
      end

          
          
      
      

