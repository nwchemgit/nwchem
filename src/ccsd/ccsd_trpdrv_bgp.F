      subroutine ccsd_trpdrv_bgp(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,einv,
     &     eccsdt,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi,
     $     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     $     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id: ccsd_trpdrv.F,v 2.10 2004-12-16 23:02:58 edo Exp $
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
c     
      double precision t1(*)
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision eorb(*),einv(nvir,nvir)
      double precision Tij(*), Tkj(*), Tia(*), Tka(*)
      double precision Xia(*), Xka(*)
      double precision Jia(*), Jka(*), Jij(*), Jkj(*)
      double precision Kia(*), Kka(*), Kij(*), Kkj(*)
      double precision Dja(*), Djka(*), Djia(*)
c
      double precision emp4,emp5,denom
      double precision emp4a,emp5a,emp4b,emp5b
      double precision eaijk,eorbaux,eccsdt
#ifdef DEBUG_PRINT
      double precision tt0,tt1,trp_time(26)
      integer tt
#endif
      integer jump,dszvv
      integer g_objo,g_objv,g_coul,g_exch
      integer ncor,nocc,nvir
      integer iprt,inode,next,nxtask,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk
      integer start, end, len
      external nxtask
c
c     
      nodes = ga_nnodes()
      iam = ga_nodeid()
C     
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
      inode=-1
      next=nxtask(nodes, 1)
c     
#ifdef DEBUG_PRINT
      do tt = 1, 26
        trp_time(tt) = 0.0d0
      enddo
#endif
c
      print*,'BlueGene/P CCSD(T) code'
c
      do klo = 1, nocc, kchunk
         akold=0
         khi = min(nocc, klo+kchunk-1)
         do a=oseg_lo,oseg_hi
            av=a-ncor-nocc
            do j=1,nocc
               inode=inode+1
               if (inode.eq.next)then
c     
c     Get Dja = Dci,ja for given j, a, all ci
c     
                  start = 1 + (j-1)*lnov
                  len   = lnov
                  end   = start + len - 1
#ifdef DEBUG_PRINT
                  tt0 = ga_wtime()
#endif
                  call ga_get(g_objv,start,end,av,av,Dja,len)
#ifdef DEBUG_PRINT
                  tt1 = ga_wtime()
                  trp_time(1) = trp_time(1) + (tt1-tt0)
#endif
c     
c     Get Tkj = T(b,c,k,j) for given j, klo<=k<=khi, all bc
c     
                  start = (klo-1)*lnvv + 1
                  len   = (khi-klo+1)*lnvv
                  end   = start + len - 1
#ifdef DEBUG_PRINT
                  tt0 = ga_wtime()
#endif
                  call ga_get(g_objo,start,end,j,j,Tkj,len)
#ifdef DEBUG_PRINT
                  tt1 = ga_wtime()
                  trp_time(2) = trp_time(2) + (tt1-tt0)
#endif
c     
c     Get Jkj = J(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
#ifdef DEBUG_PRINT
                  tt0 = ga_wtime()
#endif
                  call ga_get(g_objo,start,end,j,j,Jkj,len)
#ifdef DEBUG_PRINT
                  tt1 = ga_wtime()
                  trp_time(3) = trp_time(3) + (tt1-tt0)
#endif
c     
c     Get Kkj = K(c,l,k,j) for given j, klo<=k<=khi, all cl
c     
                  start = lnovv + lnoov + (klo-1)*lnov + 1
                  len   = (khi-klo+1)*lnov
                  end   = start + len - 1
#ifdef DEBUG_PRINT
                  tt0 = ga_wtime()
#endif
                  call ga_get(g_objo,start,end,j,j,Kkj,len)
#ifdef DEBUG_PRINT
                  tt1 = ga_wtime()
                  trp_time(4) = trp_time(4) + (tt1-tt0)
#endif
c     
                  if (akold .ne. a) then
                     akold = a
c     
c     Get Jka = J(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_coul,1,lnvv,start,end,Jka,lnvv)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(5) = trp_time(5) + (tt1-tt0)
#endif
c     
c     Get Kka = K(b,c,k,a) for given a, klo<=k<=khi, all bc
c     
                     start = (a-oseg_lo)*nocc + klo
                     len   = (khi-klo+1)
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_exch,1,lnvv,start,end,Kka,lnvv)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(6) = trp_time(6) + (tt1-tt0)
#endif
c     
c     Get Tka = Tbl,ka for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,Tka,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(7) = trp_time(7) + (tt1-tt0)
#endif
c     
c     Get Xka = Tal,kb for given a, klo<=k<=khi, all bl
c     
                     start = 1 + lnoov + lnoov + (klo-1)*lnov
                     len   = (khi-klo+1)*lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,Xka,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(8) = trp_time(8) + (tt1-tt0)
#endif
                  endif
c     
c     Get Djka = Dcj,ka for given j, a, klo<=k<=khi, all c
c     
                  do k = klo, khi
                     start = 1 + (j-1)*nvir + (k-1)*lnov
                     len   = nvir
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,
     $                    Djka(1+(k-klo)*nvir),len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(9) = trp_time(9) + (tt1-tt0)
#endif
                  enddo
c     
                  do i=1,nocc
c     
c     Get Tij = T(b,c,i,j) for given j, i, all bc
c     
                     start = (i-1)*lnvv + 1
                     len   = lnvv
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objo,start,end,j,j,Tij,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(10) = trp_time(10) + (tt1-tt0)
#endif
c     
c     Get Jij = J(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objo,start,end,j,j,Jij,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(11) = trp_time(11) + (tt1-tt0)
#endif
c     
c     Get Kij = K(c,l,i,j) for given j, i, all cl
c     
                     start = lnovv + lnoov + (i-1)*lnov + 1
                     len   = lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objo,start,end,j,j,Kij,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(12) = trp_time(12) + (tt1-tt0)
#endif
c     
c     Get Jia = J(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_coul,1,lnvv,start,end,Jia,lnvv)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(13) = trp_time(13) + (tt1-tt0)
#endif
c     
c     Get Kia = K(b,c,i,a) for given a, i, all bc
c     
                     start = (a-oseg_lo)*nocc + i
                     len   = 1
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_exch,1,lnvv,start,end,Kia,lnvv)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(14) = trp_time(14) + (tt1-tt0)
#endif
c     
c     Get Dia = Dcj,ia for given j, i, a, all c
c     
                     start = 1 + (j-1)*nvir + (i-1)*lnov
                     len   = nvir
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,Djia,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(15) = trp_time(15) + (tt1-tt0)
#endif
c     
c     Get Tia = Tbl,ia for given a, i, all bl
c     
                     start = 1 + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,Tia,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(16) = trp_time(16) + (tt1-tt0)
#endif
c     
c     Get Xia = Tal,ib for given a, i, all bl
c     
                     start = 1 + lnoov + lnoov + (i-1)*lnov
                     len   = lnov
                     end   = start + len - 1
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objv,start,end,av,av,Xia,len)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     trp_time(17) = trp_time(17) + (tt1-tt0)
#endif
c     
                     do k=klo,min(khi,i)
c
                        emp4a = 0.0d0
                        emp5a = 0.0d0
                        emp4b = 0.0d0
                        emp5b = 0.0d0
c
                        dszvv=8*lnvv ! assuming 8 bytes per double...
                        call f_memzero(f1n,dszvv)
                        call f_memzero(f2n,dszvv)
                        call f_memzero(f3n,dszvv)
                        call f_memzero(f4n,dszvv)
                        call f_memzero(f1t,dszvv)
                        call f_memzero(f2t,dszvv)
                        call f_memzero(f3t,dszvv)
                        call f_memzero(f4t,dszvv)
c     
c     sum(d) (Jia, Kia)bd * Tkj,cd -> Fbc
c     
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
!                        call ccsd_dovvv_bgp(Jia, Kia,Tkj(1+(k-klo)*lnvv),f1n,f2n,f3n,f4n,nocc,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,Jia,
     1                       nvir,Tkj(1+(k-klo)*lnvv),nvir,1.0d0,
     2                       f1n,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,Kia,
     1                       nvir,Tkj(1+(k-klo)*lnvv),nvir,1.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,Jia,
     1                       nvir,Tkj(1+(k-klo)*lnvv),nvir,1.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,Kia,
     1                       nvir,Tkj(1+(k-klo)*lnvv),nvir,1.0d0,
     2                       f4n,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(18) = trp_time(18) + (tt1-tt0)
#endif
c     
c     sum(d) (Jka, Kka)bd * Tij,cd -> Fbc
c     
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
!                        call ccsd_dovvv_bgp(Jka(1+(k-klo)*lnvv),Kka(1+(k-klo)*lnvv),Tij,f1t,f2t,f3t,f4t,nocc,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,1.0d0,
     2                       f1t,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,1.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,1.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,1.0d0,
     2                       f4t,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(19) = trp_time(19) + (tt1-tt0)
#endif
c     
c     sum(l) (Jij, Kij)cl  * Tkl,ab -> Fbc
c     
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
!                        call ccsd_doooo_bgp(Jkj(1+(k-klo)*lnov),Kkj(1+(k-klo)*lnov),Tia,Xia,f1n,f2n,f3n,f4n,nocc,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f4n,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(20) = trp_time(20) + (tt1-tt0)
#endif
c     
c     sum(l) (Jkj, Kkj)cl  * Tli,ba -> Fbc
c     
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
!                        call ccsd_doooo_bgp(Jij, Kij,Tka(1+(k-klo)*lnov),Xka(1+(k-klo)*lnov),f1t,f2t,f3t,f4t,nocc,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(21) = trp_time(21) + (tt1-tt0)
#endif
c
                        if (iprt.gt.50)then
                           call prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
     $                          f4t, nvir)
                        end if
c
                        eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)
     $                       -eorb(a)
c
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
                        jump = ncor+nocc
                        do b=1,nvir
                          eorbaux=eorb(jump+b)
                          do c=1,nvir
                            einv(c,b)=1.0d0/(eaijk-eorbaux-eorb(jump+c))
                            einv(b,c)=einv(c,b)
                          enddo
                        enddo
c
!                        call ccsd_tengy4_bgp(f1n,f1t,f2n,f2t,
!     1                                       f3n,f3t,f4n,f4t,
!     2                                       einv,emp4a,nvir)
                        call ccsd_tengy4ab_bgp(f1n,f1t,f2n,f2t,
     1                                         f3n,f3t,f4n,f4t,
     2                                         einv,emp4a,emp4b,nvir)
!                        call ccsd_tengy5_bgp(f1n,f1t,f2n,f2t,
!     1                                       f3n,f3t,f4n,f4t,
!     2                                       Dja(1+(i-1)*nvir),
!     3                                       Djia,
!     4                                       t1((k-1)*nvir+1),
!     5                                       einv,emp5a,nvir)
                        call ccsd_tengy5ab_bgp(f1n,f1t,f2n,f2t,
     1                                         f3n,f3t,f4n,f4t,
     2                                         Dja(1+(i-1)*nvir),    ! A
     3                                         Djia,                 ! A
     4                                         t1((k-1)*nvir+1),     ! A
     5                                         Dja(1+(k-1)*nvir),    ! B
     6                                         Djka(1+(k-klo)*nvir), ! B 
     7                                         t1((i-1)*nvir+1),     ! B
     8                                         einv,emp5a,emp5b,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(22) = trp_time(22) + (tt1-tt0)
#endif
c     
!                        if (i.ne.k) then
!#ifdef DEBUG_PRINT
!                           tt0 = ga_wtime()
!#endif
!                           call ccsd_tengy4_bgp(f1t,f1n,f2t,f2n,
!     1                                          f3t,f3n,f4t,f4n,
!     2                                          einv,emp4b,nvir)
!                           call ccsd_tengy5_bgp(f1t,f1n,f2t,f2n,
!     1                                          f3t,f3n,f4t,f4n,
!     2                                          Dja(1+(k-1)*nvir),
!     3                                          Djka(1+(k-klo)*nvir), 
!     4                                          t1((i-1)*nvir+1),
!     5                                          einv,emp5b,nvir)
!#ifdef DEBUG_PRINT
!                           tt1 = ga_wtime()
!                           trp_time(23) = trp_time(23) + (tt1-tt0)
!#endif
c     
!                        else
                        if (i.eq.k) then
                           emp4b = 0.0d0
                           emp5b = 0.0d0
                        end if ! (i.eq.k)
!                        end if ! (i.ne.k)
c
                        emp4 = emp4 + emp4a + emp4b
                        emp5 = emp5 + emp5a + emp5b
c
                     end do
                  end do
                  if (iprt.gt.50)then
                     write(6,1234)iam,a,j,emp4,emp5
 1234                format(' iam aijk',3i5,2e15.5)
                  end if
                  next=nxtask(nodes, 1)
               end if
            end do
         end do
      end do
c     
#ifdef DEBUG_PRINT
c
      do tt = 1, 17
        trp_time(24) = trp_time(24) + trp_time(tt)
      enddo
c
      do tt = 18, 21
        trp_time(25) = trp_time(25) + trp_time(tt)
      enddo
c
      do tt = 22, 23
        trp_time(26) = trp_time(26) + trp_time(tt)
      enddo
c
      call ga_sync()
      if (IAm.eq.0) write(6,87)
   87 format(2x,'node',6x,'ga_get',9x,'dgemm',10x,'tengy')
      call ga_sync()
      write(6,88) IAm,trp_time(24),trp_time(25),trp_time(26)
   88 format(i7,3e15.5)
c
#endif
c
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c     
      end

      subroutine ccsd_tengy4_bgp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1                           einv,emp4,nvir)
      implicit none
      integer nvir,b,c
      double precision emp4,denom
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision einv(nvir,nvir)
c
      do b=1,nvir
        do c=1,nvir
          denom=einv(c,b)
          emp4=emp4+denom*
     &         ( f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b) )*
     &         ( f1t(b,c)-2.0d0*f2t(b,c)-2.0d0*f3t(b,c)+f4t(b,c) )
          emp4=emp4-denom*
     &         (       f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b) )*
     &         ( 2.0d0*f1t(b,c)-f2t(b,c)-f3t(b,c)+2.0d0*f4t(b,c) )
          emp4=emp4+3.0d0*denom*
     &         ( f1n(b,c)*(f1n(b,c)+f3n(c,b)+2.0d0*f4t(c,b))+
     &           f2n(b,c)*f2t(c,b) + f3n(b,c)*f4t(b,c) )
        enddo
      enddo
c
      return
      end

      subroutine ccsd_tengy5_bgp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1                           dintc,dintx,t1v,einv,emp5,nvir)
      implicit none
      integer nvir,b,c
      double precision emp5,denom
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision dintc(nvir),dintx(nvir),t1v(nvir)
      double precision einv(nvir,nvir)
c
      do b=1,nvir
        do c=1,nvir
          denom=einv(c,b)
          emp5=emp5+denom*t1v(b)*dintx(c)*
     &        (
     &                 f1t(b,c)+f2n(b,c)+f4n(c,b)
     &         +4.0d0*(f3n(b,c)+f4t(b,c)+f1n(c,b))
     &         -2.0d0*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                 f1n(b,c)+f2t(b,c)+f3n(c,b))
     &        )
          emp5=emp5+denom*t1v(b)*dintc(c)*
     &        (
     &                 f1n(b,c)+f4n(b,c)+f1t(c,b)
     &         -2.0d0*(f2n(b,c)+f3n(b,c)+f2t(c,b))
     &        )
        enddo
      enddo
c
      return
      end



      subroutine ccsd_tengy4ab_bgp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1                             einv,emp4a,emp4b,nvir)
      implicit none
      integer nvir,b,c
      double precision emp4a,emp4b,denom
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision einv(nvir,nvir)
c
      do b=1,nvir
        do c=1,nvir
          denom=einv(c,b)
c
          emp4a=emp4a+denom*
     &         ( f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b) )*
     &         ( f1t(b,c)-2.0d0*f2t(b,c)-2.0d0*f3t(b,c)+f4t(b,c) )
          emp4a=emp4a+denom*
     &         (       f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b) )*
     &         (-2.0d0*f1t(b,c)+f2t(b,c)+f3t(b,c)-2.0d0*f4t(b,c) )
          emp4a=emp4a+3.0d0*denom*f1n(b,c)*
     &         ( f1n(b,c)+f3n(c,b)+2.0d0*f4t(c,b) )
          emp4a=emp4a+3.0d0*denom*f2n(b,c)*f2t(c,b)
          emp4a=emp4a+3.0d0*denom*f3n(b,c)*f4t(b,c)
c
          emp4b=emp4b+denom*
     &         ( f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b) )*
     &         ( f1n(b,c)-2.0d0*f2n(b,c)-2.0d0*f3n(b,c)+f4n(b,c) )
          emp4b=emp4b+denom*
     &         (       f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b) )*
     &         (-2.0d0*f1n(b,c)+f2n(b,c)+f3n(b,c)-2.0d0*f4n(b,c) )
          emp4b=emp4b+3.0d0*denom*f1t(b,c)*
     &         ( f1t(b,c)+f3t(c,b)+2.0d0*f4n(c,b) )
          emp4b=emp4b+3.0d0*denom*f2t(b,c)*f2n(c,b)
          emp4b=emp4b+3.0d0*denom*f3t(b,c)*f4n(b,c)
c
        enddo
      enddo
c
      return
      end



      subroutine ccsd_tengy5ab_bgp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1                             dintcA,dintxA,t1vA,
     2                             dintcB,dintxB,t1vB,
     3                             einv,emp5a,emp5b,nvir)
      implicit none
      integer nvir,b,c
      double precision emp5a,emp5b,denom
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision dintcA(nvir),dintxA(nvir),t1vA(nvir)
      double precision dintcB(nvir),dintxB(nvir),t1vB(nvir)
      double precision einv(nvir,nvir)
c
      do b=1,nvir
        do c=1,nvir
          denom=einv(c,b)
c
          emp5a=emp5a+denom*t1vA(b)*dintxA(c)*
     &         (
     &                  f1t(b,c)+f2n(b,c)+f4n(c,b)
     &          +4.0d0*(f3n(b,c)+f4t(b,c)+f1n(c,b))
     &          -2.0d0*(f3t(b,c)+f4n(b,c)+f2n(c,b))
     &          -2.0d0*(f1n(b,c)+f2t(b,c)+f3n(c,b))
     &         )
c
          emp5b=emp5b+denom*t1vB(b)*dintcB(c)*
     &         (
     &                  f1t(b,c)+f4t(b,c)+f1n(c,b)
     &          -2.0d0*(f2t(b,c)+f3t(b,c)+f2n(c,b))
     &         )
c
          emp5a=emp5a+denom*t1vA(b)*dintcA(c)*
     &         (
     &                  f1n(b,c)+f4n(b,c)+f1t(c,b)
     &          -2.0d0*(f2n(b,c)+f3n(b,c)+f2t(c,b))
     &         )
c
          emp5b=emp5b+denom*t1vB(b)*dintxB(c)*
     &         (
     &                  f1n(b,c)+f2t(b,c)+f4t(c,b)
     &          +4.0d0*(f3t(b,c)+f4n(b,c)+f1t(c,b))
     &          -2.0d0*(f3n(b,c)+f4t(b,c)+f2t(c,b))
     &          -2.0d0*(f1t(b,c)+f2n(b,c)+f3t(c,b))
     &         )
c
        enddo
      enddo
c
      return
      end

