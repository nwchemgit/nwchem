      subroutine detci_onepdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                         nexa, nexb, osym, mul, exa, exb,
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
          ii = exa(5,iex,ia)
          jj = exa(6,iex,ia)
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
          ii = exb(5,iex,ib)
          jj = exb(6,iex,ib)
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
     $                       nexa, nexb, osym, mul, exa, exb,
     $                       g_civec, twopdm )
      implicit none
#include "global.fh"
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
      integer g_civec                                 ! [input] CI-vector
c$$$      double precision civec(nstrb,nstra)             ! [input] CI-vector
      double precision twopdm(norb,norb,norb,norb)    ! [output] 2-pdm 
c
c     
      integer ia, ib, ja, jb, ka, kb
      integer iex, jex, ij, kl, i, j, l, k
      double precision xx, yy, zz, p1, p2
c
c
c
c   a  a
c  E  E      contribution
c   ij kl
c
      call dfill((norb*norb*norb*norb),0.d0,twopdm,1)
      do ia=1,nstra
        do iex=1,nexa
          ja = exa(1,iex,ia)
          ij = exa(3,iex,ia)
          p1 = exa(4,iex,ia)
          j  = exa(5,iex,ia)
          i  = exa(6,iex,ia)
          do jex=1,nexa
            ka = exa(1,jex,ja)
            kl = exa(3,jex,ja)
            p2 = exa(4,jex,ja) * p1 * 0.5d0
            l  = exa(5,jex,ja)
            k  = exa(6,jex,ja)
            xx = 0.d0
            do ib=1,nstrb
              call ga_get(g_civec,ib,ib,ka,ka,yy,1)
              call ga_get(g_civec,ib,ib,ia,ia,zz,1)
              xx = xx + yy*zz
c$$$              xx = xx + civec(ib,ka)*civec(ib,ia)
c$$$              IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $          PRINT*,':: aa',I,I,K,K,' IstrA, JstrA, IstrB, JstrB ',
c$$$     $                 ia, ka, ib, ib, ' C1 C2 ',yy,zz
            enddo
            xx = xx * p2
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(K.EQ.3).AND.(L.EQ.1))
c$$$     $        WRITE(6,833) 'AA',I,J,K,L,XX,TWOPDM(I,J,K,L)
          enddo
        enddo
      enddo
 833  FORMAT(A2,5X,4I5,2E20.10)
c
c   b  b
c  E  E      contribution
c   ij kl
c
      do ib=1,nstrb
        do iex=1,nexb
          jb = exb(1,iex,ib)
          ij = exb(3,iex,ib)
          p1 = exb(4,iex,ib)
          j  = exb(5,iex,ib)
          i  = exb(6,iex,ib)
          do jex=1,nexb
            kb = exb(1,jex,jb)
            kl = exb(3,jex,jb)
            p2 = exb(4,jex,jb) * p1 * 0.5d0
            l  = exb(5,jex,jb)
            k  = exb(6,jex,jb)
            xx = 0.d0
            do ia=1,nstra
              call ga_get(g_civec,ib,ib,ia,ia,yy,1)
              call ga_get(g_civec,kb,kb,ia,ia,zz,1)
              xx = xx + yy*zz
c$$$              xx = xx + civec(ib,ia)*civec(kb,ia)
c$$$              IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $          PRINT*,':: bb',I,I,K,K,' IstrA, JstrA, IstrB, JstrB ',
c$$$     $                 ia, ia, ib, kb, ' C1 C2 ',yy,zz
            enddo
            xx = xx * p2
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
c$$$            IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $        PRINT*,'>> bb',I,I,K,K,XX,TWOPDM(I,J,K,L)
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(K.EQ.3).AND.(L.EQ.1))
c$$$     $        WRITE(6,833) 'BB',I,J,K,L,XX,TWOPDM(I,J,K,L)
          enddo
        enddo
      enddo
c
c   a  b
c  E  E      contribution
c   ij kl
c
      do ia=1,nstra
        do ib=1,nstrb
          call ga_get(g_civec,ib,ib,ia,ia,yy,1)
c$$$          yy = civec(ib,ia)
          do iex=1,nexa
            ja = exa(1,iex,ia)
            ij = exa(3,iex,ia)
            p1 = exa(4,iex,ia)
            j  = exa(5,iex,ia)
            i  = exa(6,iex,ia)
            do jex=1,nexb
              jb = exb(1,jex,ib)
              kl = exb(3,jex,ib)
              p2 = exb(4,jex,ib) * p1 * 0.5d0
              l  = exb(5,jex,ib)
              k  = exb(6,jex,ib)
              call ga_get(g_civec,jb,jb,ja,ja,zz,1)
              xx = zz*yy*p2
c$$$              IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $          PRINT*,':: ab',I,I,K,K,' IstrA, JstrA, IstrB, JstrB ',
c$$$     $                 ia, ja, ib, jb, ' C1 C2 ',yy,zz
c$$$              xx = civec(jb,ja)*yy*p2
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
c$$$              IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $             PRINT*,'>> ab',I,I,K,K,XX,TWOPDM(I,J,K,L)
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(K.EQ.3).AND.(L.EQ.1))
c$$$     $        WRITE(6,833) 'AB',I,J,K,L,XX,TWOPDM(I,J,K,L)
            enddo
          enddo
        enddo
      enddo
c
c   b  a
c  E  E      contribution
c   ij kl
c
      do ia=1,nstra
        do ib=1,nstrb
          call ga_get(g_civec,ib,ib,ia,ia,yy,1)
          do jex=1,nexb
            jb = exb(1,jex,ib)
            ij = exb(3,jex,ib)
            p1 = exb(4,jex,ib)
            i  = exb(5,jex,ib)
            j  = exb(6,jex,ib)
            do iex=1,nexa
              ja = exa(1,iex,ia)
              kl = exa(3,iex,ia)
              p2 = exa(4,iex,ia) * p1 * 0.5d0
              k  = exa(5,iex,ia)
              l  = exa(6,iex,ia)
              call ga_get(g_civec,jb,jb,ja,ja,zz,1)
              xx = zz*yy*p2
c$$$              IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $          PRINT*,':: ba',I,I,K,K,' IstrA, JstrA, IstrB, JstrB ',
c$$$     $                 ia, ja, ib, jb, ' C1 C2 ',yy,zz
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
c$$$            IF ((I.EQ.J).AND.(K.EQ.L).AND.(I.EQ.K))
c$$$     $        PRINT*,'>> ba',I,I,K,K,XX,TWOPDM(I,J,K,L)
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(K.EQ.3).AND.(L.EQ.1))
c$$$     $        WRITE(6,833) 'BA',I,J,K,L,XX,TWOPDM(I,J,K,L)
            enddo
          enddo
        enddo
      enddo
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
          i  = exa(5,iex,ia)
          l  = exa(6,iex,ia)
          xx = 0.d0
          do ib=1,nstrb
            call ga_get(g_civec,ib,ib,ia,ia,yy,1)
            call ga_get(g_civec,ib,ib,ja,ja,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(ib,ia)*civec(ib,ja)
c$$$            IF ((I.EQ.4).AND.(L.EQ.3))
c$$$     $        WRITE(6,832) 'A ',I,L,IA,IB,JA,IB,YY,ZZ,XX
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(L.EQ.3))
c$$$     $        WRITE(6,834) 'A ',I,J,J,L,P1,P2,TWOPDM(I,J,J,L)
            twopdm(i,j,j,l) = twopdm(i,j,j,l) + p2
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
          i  = exb(5,iex,ib)
          l  = exb(6,iex,ib)
          xx = 0.d0
          do ia=1,nstra
            call ga_get(g_civec,jb,jb,ia,ia,yy,1)
            call ga_get(g_civec,ib,ib,ia,ia,zz,1)
            xx = xx + yy*zz
c$$$            xx = xx + civec(jb,ia)*civec(ib,ia)
c$$$            IF ((I.EQ.4).AND.(L.EQ.3))
c$$$     $        WRITE(6,832) 'B ',I,L,IA,IB,IA,JB,YY,ZZ,XX
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
c$$$            IF ((I.EQ.4).AND.(J.EQ.1).AND.(L.EQ.3))
c$$$     $        WRITE(6,834) 'B ',I,J,J,L,P1,P2,TWOPDM(I,J,J,L)
            twopdm(i,j,j,l) = twopdm(i,j,j,l) + p2
          enddo
        enddo
      enddo

 834  FORMAT(A2,5X,4I5,F5.1,2E20.10)
 832  FORMAT(A2,20X,2I5,5X,2I2,'<->',2I2,2X,3E20.10)

c$$$      WRITE(6,*)
c$$$      WRITE(6,'(A)') '*** TWO-PDM AFTER AB ***'
c$$$      DO I=1,NORB
c$$$        DO J=1,NORB
c$$$          WRITE(6,'(4I5,F12.6)') I,J,J,I,TWOPDM(I,J,J,I)
c$$$        ENDDO
c$$$      ENDDO
c$$$      WRITE(6,*)

c$$$      WRITE(6,*)
c$$$      WRITE(6,'(A)') '*** TWO-PDM ***'
c$$$      DO I=1,NORB
c$$$        DO J=1,NORB
c$$$          DO K=1,NORB
c$$$            DO L=1,NORB
c$$$              WRITE(6,'(4I5,F12.6)') I,J,K,L,TWOPDM(I,J,K,L)
c$$$            ENDDO
c$$$          ENDDO
c$$$        ENDDO
c$$$      ENDDO
c$$$      WRITE(6,*)
          
      return
      end








