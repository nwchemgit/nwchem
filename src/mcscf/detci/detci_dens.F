      subroutine detci_onepdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                         nexa, nexb, osym, mul, exa, exb, map,
     $                         g_civec, onepdm )
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
      integer osym(norb)                      ! [input] Orbital irreps     
      integer mul(nsym,nsym)                  ! [input] Direct product table
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer map(norb)                       ! [input] Orbital index map
c      double precision civec(nstrb,nstra)     ! [input] CI-vector
      integer g_civec
      double precision onepdm(norb,norb)      ! [output] 1-pdm 
c
      integer ia, ja, ib, jb
      integer ii, jj
      integer iex, ph
      double precision xx, yy, zz
c
c
c
      call dfill((norb*norb),0.d0,onepdm,1)
      do ia=1,nstra
        do iex=1,nexa
          ja = exa(1,iex,ia)
          ph = exa(4,iex,ia)
          ii = map(exa(5,iex,ia))
          jj = map(exa(6,iex,ia))
          xx = 0.d0
          do ib=1,nstrb
            call ga_get(g_civec,ib,ib,ia,ia,yy,1)
            call ga_get(g_civec,ib,ib,ja,ja,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(ib,ia)*civec(ib,ja)
          enddo
          onepdm(ii,jj) = onepdm(ii,jj) + ph*xx
        enddo
      enddo
c
c
c
      do ib=1,nstrb
        do iex=1,nexb
          jb = exb(1,iex,ib)
          ph = exb(4,iex,ib)
          ii = map(exb(5,iex,ib))
          jj = map(exb(6,iex,ib))
          xx = 0.d0
          do ia=1,nstra
            call ga_get(g_civec,ib,ib,ia,ia,yy,1)
            call ga_get(g_civec,jb,jb,ia,ia,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(ib,ia)*civec(jb,ia)
          enddo
          onepdm(ii,jj) = onepdm(ii,jj) + ph*xx
        enddo
      enddo

      return
      end














c
c
c  Two-particle density matrix 
c
c
c              1
c  Gamma     = - <0| E  E   - delta  E   |0>
c       ijkl   2      ij kl        jk il
c
c
c
c
      subroutine detci_twopdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                       nexa, nexb, osym, mul, exa, exb, map,
     $                       g_civec, twopdm )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "detci.fh"
      integer norb                                    ! [input] Orbitals
      integer nsym                                    ! [input] Irreps
      integer nela                                    ! [input] Alpha electrons
      integer nelb                                    ! [input] Beta electrons
      integer nstra                                   ! [input] Alpha strings
      integer nstrb                                   ! [input] Beta strings
      integer nexa                                    ! [input] Alpha excitations
      integer nexb                                    ! [input] Beta excitations
      integer osym(norb)                              ! [input] Orbital irreps     
      integer mul(nsym,nsym)                          ! [input] Direct product table
      integer exa(6,nexa,nstra)                       ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)                       ! [input] Beta excitation lookup table
      integer map(norb)                               ! [input] Orbital index map
      integer g_civec                                 ! [input] CI-vector
c$$$      double precision civec(nstrb,nstra)             ! [input] CI-vector
      double precision twopdm(norb,norb,norb,norb)    ! [output] 2-pdm 
c
c     
      integer ia, ib, ja, jb, ka, kb
      integer iex, jex, ij, kl, i, j, l, k, jj
      double precision xx, yy, zz, p1, p2
      integer l_t, k_t, l_s, k_s
      double precision ddot
      external ddot
c
c
c
c   a  a
c  E  E      contribution
c   ij kl
c
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0)
      call dfill((norb*norb*norb*norb),0.d0,twopdm,1)
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          ij = exa(3,iex,ia)
          p1 = exa(4,iex,ia)
          j  = map(exa(5,iex,ia))
          i  = map(exa(6,iex,ia))
          do jex=1,nexa
            ka = exa(1,jex,ja)
            kl = exa(3,jex,ja)
            p2 = exa(4,jex,ja) * p1 * 0.5d0
            l  = map(exa(5,jex,ja))
            k  = map(exa(6,jex,ja))
            xx = 0.d0
            call ga_get(g_civec,1,nstrb,ka,ka,dbl_mb(k_s),1)
            xx = ddot(nstrb,dbl_mb(k_s),1,dbl_mb(k_t),1) * p2
c$$$            do ib=1,nstrb
c$$$              call ga_get(g_civec,ib,ib,ka,ka,yy,1)
c$$$              call ga_get(g_civec,ib,ib,ia,ia,zz,1)
c$$$              xx = xx + yy*zz
c$$$              xx = xx + civec(ib,ka)*civec(ib,ia)
c$$$            enddo
c$$$            xx = xx * p2
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
          enddo
        enddo
      enddo
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
c
c   b  b
c  E  E      contribution
c   ij kl
c
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0)
      do ib=1,nstrb
        call ga_get(g_civec,ib,ib,1,nstra,dbl_mb(k_t),1)
        do iex=1,nexb
          jb = exb(1,iex,ib)
          ij = exb(3,iex,ib)
          p1 = exb(4,iex,ib)
          j  = map(exb(5,iex,ib))
          i  = map(exb(6,iex,ib))
          do jex=1,nexb
            kb = exb(1,jex,jb)
            kl = exb(3,jex,jb)
            p2 = exb(4,jex,jb) * p1 * 0.5d0
            l  = map(exb(5,jex,jb))
            k  = map(exb(6,jex,jb))
            xx = 0.d0
            call ga_get(g_civec,kb,kb,1,nstra,dbl_mb(k_s),1)
            xx = ddot(nstra,dbl_mb(k_s),1,dbl_mb(k_t),1) * p2
c$$$            do ia=1,nstra
c$$$              call ga_get(g_civec,ib,ib,ia,ia,yy,1)
c$$$              call ga_get(g_civec,kb,kb,ia,ia,zz,1)
c$$$              xx = xx + yy*zz
c$$$              xx = xx + civec(ib,ia)*civec(kb,ia)
c$$$            enddo
c$$$            xx = xx * p2
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
          enddo
        enddo
      enddo
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
c
c   a  b
c  E  E      contribution
c   ij kl
c
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      do ia=1,nstra
        do ib=1,nstrb
          call ga_get(g_civec,ib,ib,ia,ia,yy,1)
c$$$          yy = civec(ib,ia)
          do iex=1,nexa
            ja = exa(1,iex,ia)
            ij = exa(3,iex,ia)
            p1 = exa(4,iex,ia)
            j  = map(exa(5,iex,ia))
            i  = map(exa(6,iex,ia))
            call ga_get(g_civec,1,nstrb,ja,ja,dbl_mb(k_t),1)
            do jex=1,nexb
              jb = exb(1,jex,ib)
              p2 = exb(4,jex,ib) * p1 * 0.5d0
              l  = map(exb(5,jex,ib))
              k  = map(exb(6,jex,ib))
              xx = p2 * yy * dbl_mb(k_t+jb-1)
c$$$              call ga_get(g_civec,jb,jb,ja,ja,zz,1)
c$$$              xx = zz*yy*p2
c$$$              xx = civec(jb,ja)*yy*p2
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
              twopdm(l,k,j,i) = twopdm(l,k,j,i) + xx
            enddo
          enddo
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
c
c
c  This contribution already taken care of above
c
c   b  a
c  E  E      contribution
c   ij kl
c
c$$$      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
c$$$     $     call errquit('cannot allocate temp',0)
c$$$      do ia=1,nstra
c$$$        do ib=1,nstrb
c$$$          call ga_get(g_civec,ib,ib,ia,ia,yy,1)
c$$$          do iex=1,nexa
c$$$            ja = exa(1,iex,ia)
c$$$            kl = exa(3,iex,ia)
c$$$            p1 = exa(4,iex,ia)
c$$$            k  = map(exa(5,iex,ia))
c$$$            l  = map(exa(6,iex,ia))
c$$$            call ga_get(g_civec,1,nstrb,ja,ja,dbl_mb(k_t),1)
c$$$            do jex=1,nexb
c$$$              jb = exb(1,jex,ib)
c$$$              p2 = exb(4,jex,ib) * p1 * 0.5d0
c$$$              i  = map(exb(5,jex,ib))
c$$$              j  = map(exb(6,jex,ib))
c$$$              xx = p2 * yy * dbl_mb(k_t+jb-1)
c$$$c$$$              call ga_get(g_civec,jb,jb,ja,ja,zz,1)
c$$$c$$$              xx = zz*yy*p2
c$$$              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      if (.not.ma_pop_stack(l_t))
c$$$     $   call errquit('cannot pop stack',0)
c
c
c   a  
c  E  delta      contribution
c   il     jk
c
c
      do ia=1,nstra
        do iex=1,nexa
          ja = exa(1,iex,ia)
          p1 = exa(4,iex,ia)
          i  = map(exa(5,iex,ia))
          l  = map(exa(6,iex,ia))
          xx = 0.d0
          do ib=1,nstrb
            call ga_get(g_civec,ib,ib,ia,ia,yy,1)
            call ga_get(g_civec,ib,ib,ja,ja,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(ib,ia)*civec(ib,ja)
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
            jj = map(j)
            twopdm(i,jj,jj,l) = twopdm(i,jj,jj,l) + p2
          enddo
        enddo
      enddo
c
c   b
c  E  delta      contribution
c   il     jk
c
c
      do ib=1,nstrb
        do iex=1,nexb
          jb = exb(1,iex,ib)
          p1 = exb(4,iex,ib)
          i  = map(exb(5,iex,ib))
          l  = map(exb(6,iex,ib))
          xx = 0.d0
          do ia=1,nstra
            call ga_get(g_civec,jb,jb,ia,ia,yy,1)
            call ga_get(g_civec,ib,ib,ia,ia,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(jb,ia)*civec(ib,ia)
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
            jj = map(j)
            twopdm(i,jj,jj,l) = twopdm(i,jj,jj,l) + p2
          enddo
        enddo
      enddo

 834  FORMAT(A2,5X,4I5,F5.1,2E20.10)
 832  FORMAT(A2,20X,2I5,5X,2I2,'<->',2I2,2X,3E20.10)

          
      return
      end








