      subroutine detci_onepdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                         nexa, nexb, osym, exa, exb, map,
     $                         g_civec, g_civect, onepdm )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
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
      integer exa(6,nexa,nstra)               ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)               ! [input] Beta excitation lookup table
      integer map(norb)                       ! [input] Orbital index map
      integer g_civec                         ! [input] CI-vector
      integer g_civect                        ! [input] CI-vector transpose
      double precision onepdm(norb,norb)      ! [output] 1-pdm 
c
      integer ia, ja, ib, jb
      integer myid, rlo, rhi, cilo, cihi, dilo, dihi
      integer l_t, k_t, l_s, k_s, lds
      integer g_dtmp
      integer ii, jj, iex, ph
      double precision xx, yy, zz
      double precision t1, t2
      DOUBLE PRECISION CRAP(1000)
c
c
c
      myid = ga_nodeid()
      if (.not.(ga_create(MT_DBL, norb, norb, 'd', norb, 1, g_dtmp )))
     $   call errquit('detci_twopdm: cannot create global',0)
      call ga_zero(g_dtmp)
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_onepdm: wrong distrib for CI vector',0)
      call ga_distribution( g_civect, myid, rlo, rhi, dilo, dihi )
      if (((dilo.ne.0).and.(dihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstra)))
     $  call errquit('detci_onepdm: wrong distrib for trp(CI) vec',0)
c
      T1 = UTIL_CPUSEC()
#ifdef DETCI_DENSSERIAL
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
          enddo
          onepdm(ii,jj) = onepdm(ii,jj) + ph*xx
        enddo
      enddo
#else
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      do ia=1,nstra
        call ga_get( g_civec, 1, nstrb, ia, ia, dbl_mb(k_t), 1 )
        do iex=1,nexa
          ja = exa(1,iex,ia)
          if ((ja.ge.cilo).and.(ja.le.cihi)) then
            ph = exa(4,iex,ia)
            ii = map(exa(5,iex,ia))
            jj = map(exa(6,iex,ia))
            call ga_access(g_civec, 1, nstrb, ja, ja, k_s, lds )
            xx = ddot(nstrb, dbl_mb(k_t), 1, dbl_mb(k_s), 1) * ph
            call ga_release( g_civec, 1, nstrb, ja, ja )
            call ga_acc( g_dtmp, ii, ii, jj, jj, xx, 1, 1.d0 )
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
      call ga_sync()
#endif
      T1 = UTIL_CPUSEC() - T1
c
c
c
      T2 = UTIL_CPUSEC()
#ifdef DETCI_DENSSERIAL
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
          enddo
          onepdm(ii,jj) = onepdm(ii,jj) + ph*xx
        enddo
      enddo
#else
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      do ib=1,nstrb
        call ga_get( g_civect, 1, nstra, ib, ib, dbl_mb(k_t), 1 )
        do iex=1,nexb
          jb = exb(1,iex,ib)
          if ((jb.ge.dilo).and.(jb.le.dihi)) then
            ph = exb(4,iex,ib)
            ii = map(exb(5,iex,ib))
            jj = map(exb(6,iex,ib))
            call ga_access( g_civect, 1, nstra, jb, jb, k_s, lds )
            xx = ddot( nstra, dbl_mb(k_t), 1, dbl_mb(k_s), 1 ) * ph
            call ga_release( g_civect, 1, nstra, jb, jb )
            call ga_acc( g_dtmp, ii, ii, jj, jj, xx, 1, 1.d0 )
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
      call ga_sync()
      call ga_get(g_dtmp, 1, norb, 1, norb, onepdm, norb )
#endif
      T2 = UTIL_CPUSEC() - T2
      if (ga_nodeid().eq.0) then
        write(6,911) t1,t2
 911    format('1-pdm timings',2(/,10x,f12.3))
      endif
c
c
c      
      if (.not.ga_destroy(g_dtmp))
     $  call errquit('detci_onepdm:cannot destroy tmp',0)
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
      subroutine detci_twopdm( norb, nsym, nela, nelb, nstra, nstrb,
     $                       nexa, nexb, osym, exa, exb, map,
     $                       g_civec, g_civect, twopdm )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
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
      integer exa(6,nexa,nstra)                       ! [input] Alpha excitation lookup table
      integer exb(6,nexb,nstrb)                       ! [input] Beta excitation lookup table
      integer map(norb)                               ! [input] Orbital index map
      integer g_civec                                 ! [input] CI-vector
      integer g_civect                                ! [input] CI-vector transposed
      double precision twopdm(norb,norb,norb,norb)    ! [output] 2-pdm 
c
c     
      integer ia, ib, ja, jb, ka, kb
      integer iex, jex, ij, kl, ji, lk, jl
      integer i, j, l, k, jj
      double precision xx, yy, zz, p1, p2
      integer l_t, k_t, l_s, k_s, lds
      integer g_dentmp
      integer myid, nn, rlo, rhi, cilo, cihi, dilo, dihi
      double precision t1, t2, t3, t4, t5
      DOUBLE PRECISION CRAP(10000)
C      double precision sdot
C      external sdot
c
c
      myid = ga_nodeid()
      nn = norb*norb
      if (.not.(ga_create(MT_DBL, nn, nn, 'g dens', nn, 1, g_dentmp )))
     $   call errquit('detci_twopdm: cannot create global',0)
      call ga_zero(g_dentmp)
      call dfill((norb*norb*norb*norb),0.d0,twopdm,1)
      call ga_distribution( g_civect, myid, rlo, rhi, dilo, dihi )
      if (((dilo.ne.0).and.(dihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstra)))
     $  call errquit('detci_sigmaaa: wrong distrib for trp(CI)',0)
      call ga_distribution( g_civec, myid, rlo, rhi, cilo, cihi )
      if (((cilo.ne.0).and.(cihi.ne.-1)).and.
     $    ((rlo.ne.1).or.(rhi.ne.nstrb)))
     $  call errquit('detci_sigmaaa: wrong distrib for CI vector',0)
c
c
c
c   a  a
c  E  E      contribution
c   ij kl
c
      T1 = UTIL_CPUSEC()
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0)
#ifdef DETCI_DENSSERIAL      
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          ij = exa(3,iex,ia)
          p1 = exa(4,iex,ia)
          j  = map(exa(5,iex,ia))
          i  = map(exa(6,iex,ia))
          ij = (j-1)*norb + i
          do jex=1,nexa
            ka = exa(1,jex,ja)
            kl = exa(3,jex,ja)
            p2 = exa(4,jex,ja) * p1 * 0.5d0
            l  = map(exa(5,jex,ja))
            k  = map(exa(6,jex,ja))
            call ga_get(g_civec,1,nstrb,ka,ka,dbl_mb(k_s),1)           ! Serial 
            xx = ddot(nstrb,dbl_mb(k_s),1,dbl_mb(k_t),1) * p2
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
          enddo
        enddo
      enddo
#else
      call ga_sync()
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          p1 = exa(4,iex,ia)
          j  = map(exa(5,iex,ia))
          i  = map(exa(6,iex,ia))
          ij = (j-1)*norb + i
          do jex=1,nexa
            ka = exa(1,jex,ja)
            if ((ka.ge.cilo).and.(ka.le.cihi)) then
              p2 = exa(4,jex,ja) * p1 * 0.5d0
              l  = map(exa(5,jex,ja))
              k  = map(exa(6,jex,ja))
              kl = (l-1)*norb + k
              call ga_access(g_civec, 1, nstrb, ka, ka, k_s, lds )     ! Parallel
              xx = ddot(nstrb, dbl_mb(k_s), 1, dbl_mb(k_t), 1 ) * p2
              call ga_release( g_civec, 1, nstrb, ka, ka )
              call ga_acc( g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
            endif
          enddo
        enddo
      enddo
#endif
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
      T1 = UTIL_CPUSEC() - T1
c
c   b  b
c  E  E      contribution
c   ij kl
c
      T2 = UTIL_CPUSEC()
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp2',l_s,k_s))
     $     call errquit('cannot allocate temp',0)
#ifdef DETCI_DENSSERIAL
      do ib=1,nstrb                                                  ! Serial
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
            twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
          enddo
        enddo
      enddo
#else
      call ga_sync()
      do ib=1,nstrb                                                  ! Parallel
        call ga_get(g_civect,1,nstra,ib,ib,dbl_mb(k_t),1)
        do iex=1,nexb
          jb = exb(1,iex,ib)
          p1 = exb(4,iex,ib)
          j  = map(exb(5,iex,ib))
          i  = map(exb(6,iex,ib))
          ij = (j-1)*norb + i
          do jex=1,nexb
            kb = exb(1,jex,jb)
            if ((kb.ge.dilo).and.(kb.le.dihi)) then
              p2 = exb(4,jex,jb) * p1 * 0.5d0
              l  = map(exb(5,jex,jb))
              k  = map(exb(6,jex,jb))
              kl = (l-1)*norb + k
              call ga_access(g_civect, 1, nstra, kb, kb, k_s, lds )
              xx = ddot(nstra, dbl_mb(k_s), 1, dbl_mb(k_t), 1 ) * p2
              call ga_release( g_civect, 1, nstra, kb, kb )
              call ga_acc( g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
            endif
          enddo
        enddo
      enddo
#endif
      if (.not.ma_pop_stack(l_s))
     $   call errquit('cannot pop stack',0)
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
      T2 = UTIL_CPUSEC() - T2
c
c   a  b         b  a 
c  E  E    and  E  E     contribution
c   ij kl        ij kl
c
      T3 = UTIL_CPUSEC()
#ifndef DETCI_DENSAB      
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
#ifdef DETCI_DENSSERIAL
      do ia=1,nstra
        do ib=1,nstrb
          call ga_get(g_civec,ib,ib,ia,ia,yy,1)
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
              twopdm(i,j,k,l) = twopdm(i,j,k,l) + xx
              twopdm(l,k,j,i) = twopdm(l,k,j,i) + xx
            enddo
          enddo
        enddo
      enddo
#else
c$$$      CALL DFILL((NORB*NORB*NORB*NORB),0.d0,CRAP,1)
      call ga_sync()
      do ia=1,nstra                                                ! Parallel
        do iex=1,nexa
          ja = exa(1,iex,ia)
          p1 = exa(4,iex,ia)
          j  = map(exa(5,iex,ia))
          i  = map(exa(6,iex,ia))
          ij = (j-1)*norb + i
          ji = (i-1)*norb + j
C          call ga_get(g_civec,1,nstrb,ja,ja,dbl_mb(k_t),1)
          do ib=1,nstrb
            if ((ib.ge.dilo).and.(ib.le.dihi)) then
              call ga_access(g_civect, 1, nstra, ib, ib, k_s, lds )
              yy = dbl_mb(k_s+ia-1)
              call ga_release(g_civect, 1, nstra, ib, ib )
              do jex=1,nexb
                jb = exb(1,jex,ib)
                call ga_get(g_civec, jb, jb, ja, ja, zz, 1 )
                p2 = exb(4,jex,ib) * p1 * 0.5d0
                l  = map(exb(5,jex,ib))
                k  = map(exb(6,jex,ib))
                kl = (l-1)*norb + k
                lk = (k-1)*norb + l
C                xx = p2 * yy * dbl_mb(k_t+jb-1)
                xx = p2 * yy * zz
                call ga_acc(g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
                call ga_acc(g_dentmp, lk, lk, ji, ji, xx, 1, 1.d0 )
c$$$                CALL XXXPUT(NORB,I,J,K,L,XX,CRAP)
c$$$                CALL XXXPUT(NORB,L,K,J,I,XX,CRAP)
              enddo
            endif
          enddo
        enddo
      enddo
c$$$      CALL GA_SYNC()
c$$$      IF (GA_NODEID().EQ.0) CALL MOINTS_MATPRINT(NN,NN,CRAP)
#endif
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
#endif
      T3 = UTIL_CPUSEC() - T3
c
c   a  
c  E  delta      contribution
c   il     jk
c
c
      T4 = UTIL_CPUSEC()
#ifdef DETCI_DENSSERIAL
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
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
            jj = map(j)
            twopdm(i,jj,jj,l) = twopdm(i,jj,jj,l) + p2
          enddo
        enddo
      enddo
#else
      call ga_sync()
      if (.not.ma_push_get(MT_DBL,nstrb,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      do ia=1,nstra
        call ga_get(g_civec,1,nstrb,ia,ia,dbl_mb(k_t),1)
        do iex=1,nexa
          ja = exa(1,iex,ia)
          if ((ja.ge.cilo).and.(ja.le.cihi)) then
            p1 = exa(4,iex,ia)
            i  = map(exa(5,iex,ia))
            l  = map(exa(6,iex,ia))
            call ga_access(g_civec, 1, nstrb, ja, ja, k_s, lds )
            xx = ddot( nstrb, dbl_mb(k_t), 1, dbl_mb(k_s), 1 )
            call ga_release(g_civec, 1, nstrb, ja, ja, k_s, lds )
            yy = -0.5d0*xx*p1
            do j=1,norb
              ij = (map(j)-1)*norb + i
              jl = (l-1)*norb + map(j)
              call ga_acc(g_dentmp, ij, ij, jl, jl, yy, 1, 1.d0 )
            enddo
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
#endif
      T4 = UTIL_CPUSEC() - T4
c
c   b
c  E  delta      contribution
c   il     jk
c
c
      T5 = UTIL_CPUSEC()
#ifdef DETCI_DENSSERIAL
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
          enddo
          p2 = -0.5d0*xx*p1
          do j=1,norb
            jj = map(j)
            twopdm(i,jj,jj,l) = twopdm(i,jj,jj,l) + p2
          enddo
        enddo
      enddo
#else
      call ga_sync()
      if (.not.ma_push_get(MT_DBL,nstra,'dens:tmp1',l_t,k_t))
     $     call errquit('cannot allocate temp',0)
      do ib=1,nstrb
        call ga_get(g_civect, 1, nstra, ib, ib, dbl_mb(k_t), 1 )
        do iex=1,nexb
          jb = exb(1,iex,ib)
          if ((jb.ge.dilo).and.(jb.le.dihi)) then
            p1 = exb(4,iex,ib)
            i  = map(exb(5,iex,ib))
            l  = map(exb(6,iex,ib))
            call ga_access(g_civect, 1, nstra, jb, jb, k_s, lds )
            xx = ddot(nstra, dbl_mb(k_t), 1, dbl_mb(k_s), 1 )
            call ga_release( g_civect, 1, nstra, jb, jb )
            yy = -0.5d0*xx*p1
            do j=1,norb
              ij = (map(j)-1)*norb + i
              jl = (l-1)*norb + map(j)
              call ga_acc(g_dentmp, ij, ij, jl, jl, yy, 1, 1.d0 )
            enddo
          endif
        enddo
      enddo
      if (.not.ma_pop_stack(l_t))
     $   call errquit('cannot pop stack',0)
      call ga_sync()
      call ga_get(g_dentmp,1,nn,1,nn,twopdm,nn)
#endif
      T5 = UTIL_CPUSEC() - T5
      IF (GA_NODEID().EQ.0) THEN
        WRITE(6,881) T1,T2,T3,T4,T5
 881    FORMAT('2-pdm timings:',5(/,10x,f12.3))
      ENDIF
c
c  Clean up
c
      if (.not.ga_destroy(g_dentmp))
     $  call errquit('detci_twopdm: cannot destroy global dens',0)
          
      return
      end






















#ifdef DETCI_DENSAB
      subroutine detci_twopdm_ab( norb, nsym, nela, nelb, nstra, nstrb,
     $                            nexa, nexb, nekla, neklb,
     $                            osym, map, exa, exb, ataba, atabb, 
     $                            vrhs, vlhs, vphase,
     $                            cprime, sprime, g_civec, twopdm )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
      integer norb                                     ! [input] Orbitals                             
      integer nsym                                     ! [input] Irreps                               
      integer nela                                     ! [input] Alpha electrons                      
      integer nelb                                     ! [input] Beta electrons                       
      integer nstra                                    ! [input] Alpha strings                        
      integer nstrb                                    ! [input] Beta strings                         
      integer nexa                                     ! [input] Alpha excitations                    
      integer nexb                                     ! [input] Beta excitations                     
      integer nekla                                    ! [input] Max non-zero alpha strings for E_kl  
      integer neklb                                    ! [input] Max non-zero beta strings for E_kl   
      integer map(norb)                                ! [input] External remapping of indices        
      integer osym(norb)                               ! [input] Orbital irreps                       
      integer ataba(norb,nela,nsym)                    ! [input] Alpha arc weights                    
      integer atabb(norb,nelb,nsym)                    ! [input] Beta arc weights                     
      integer exa(6,nexa,nstra)                        ! [input] Alpha excitation lookup table        
      integer exb(6,nexb,nstrb)                        ! [input] Beta excitation lookup table         
      integer vrhs(nekla)                              ! [scratch] Array of RHS Strings for E_kl      
      integer vlhs(nekla)                              ! [scratch] Array of LHS Strings for E_kl      
      integer vphase(nekla)                            ! [scratch] Array of Phases for E_kl           
      double precision cprime(nstrb,nekla)             ! [scratch] Gathered CI-vector                 
      double precision sprime(nstrb,nekla)             ! [scratch] Gathered sigma vector              
      integer g_civec                                  ! [input] CI-vector                            
      double precision twopdm(norb,norb,norb,norb)     ! [output] Two-particle density 
c
      integer i, j, k, l
      integer ii, jj, kk, ll
      integer ij, ji, kl, lk
      integer ai, aj
      integer iph
      integer jstr, kstr
      integer ne_kl, nidx, iex
      integer relv(detci_maxelec)
      integer lelv(detci_maxelec)
      integer oidx(detci_maxorb)
      integer ip(detci_maxelec)
      integer g_dentmp
      integer ploop, numnodes, next, myid, nn
      double precision phase, xx
      double precision tx
      integer nxtval
      external nxtval
c      
c
c Initialize parallel stuff
c
      tx = util_cpusec()
      nn = norb*norb
      if (.not.(ga_create(MT_DBL, nn, nn, 'g dens', nn, 1, g_dentmp )))
     $   call errquit('detci_twopdm_ab: cannot create global',0)
      call ga_zero(g_dentmp)
      if (ga_nodeid().eq.0)
     $   call ga_put( g_dentmp, 1, nn, 1, nn, twopdm, nn )
      call ga_sync()
c
      ploop = -1
      numnodes = ga_nnodes()
      myid = ga_nodeid()
      next = nxtval(numnodes)
c
c
c Parallel loop over all alpha excitation operators
c
c            t
c      E  = a a
c       ij   i j
c
      do ai=1,norb
        do aj=1,norb
          ploop = ploop + 1
          if (ploop.eq.next) then
            ij = (map(aj)-1)*norb + map(ai)
            ji = (map(ai)-1)*norb + map(aj)
            ii = map(ai)
            jj = map(aj)
            ne_kl = 0
c
c Vector of orbital indices except create/annih indices
c Initialize pointer vector
c
            nidx = 0
            do i=1,norb
              if ((i.ne.ai).and.(i.ne.aj)) then
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
c
c Put indices into gather/scatter arrays
c
            if (nela.gt.1) then
 101          continue
              ne_kl = ne_kl + 1
              iph = 1
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           aj,relv,iph)
              call detci_ptr2elv(norb,nela,(nela-1),nidx,ip,oidx,
     $                           ai,lelv,iph)
              vphase(ne_kl) = iph
              vrhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,
     $                                    ataba,relv)
              vlhs(ne_kl) = detci_elv2str(norb,nela,nsym,osym,
     $                                    ataba,lelv)
              if (detci_getnextelv(nidx,(nela-1),ip)) goto 101
            else
              ne_kl = 1
              vrhs(i) = aj
              vlhs(i) = ai
              vphase(i) = 1
            endif
c
c Gather in left and right CI columns
c
            call detci_cigather( nstrb, nstra, ne_kl, g_civec, vlhs,
     $                           vphase, cprime )
            call detci_cigather0( nstrb, nstra, ne_kl, g_civec, vrhs,
     $                            sprime )
c
c Loop over all beta strings
c
            do jstr=1,nstrb
              do iex=1,nexb
                kstr  = exb(1,iex,jstr)
                phase = exb(4,iex,jstr)
                l     = map(exb(5,iex,jstr))
                k     = map(exb(6,iex,jstr))
                kk    = k
                ll    = l
                kl    = (l-1)*norb + k
                lk    = (k-1)*norb + l
                xx = ddot(ne_kl, cprime(jstr,1), nstrb,
     $                    sprime(kstr,1), nstrb ) * phase * 0.5d0
                call ga_acc(g_dentmp, ij, ij, kl, kl, xx, 1, 1.d0 )
                call ga_acc(g_dentmp, lk, lk, ji, ji, xx, 1, 1.d0 )
              enddo
            enddo
c
c End parallel task
c
            next = nxtval(numnodes)
          endif
        enddo
      enddo
      next = nxtval(-numnodes)
c$$$      CALL GA_SYNC()
c$$$      NN = NORB*NORB
c$$$      IF (GA_NODEID().EQ.0) CALL MOINTS_MATPRINT(NN,NN,TWOPDM)
c
c
c
      call ga_get(g_dentmp, 1, nn, 1, nn, twopdm, nn )
      if (.not.ga_destroy(g_dentmp))
     $  call errquit('detci_twopdm_ab: cannot destroy global dens',0)
      tx = util_cpusec() - tx
      if (ga_nodeid().eq.0) write(6,881) tx
 881  format('2-pdm-ab timings:',/,10x,f12.3)
c
c
c
      return
      end
#endif









      SUBROUTINE XXXPUT(N,I,J,K,L,XX,YY)
      IMPLICIT NONE
      INTEGER N,I,J,K,L
      DOUBLE PRECISION YY(N,N,N,N)
      DOUBLE PRECISION XX
      YY(I,J,K,L) = YY(I,J,K,L) + XX
      RETURN
      END

      
