      subroutine ccsd_trpdrv_bgp2(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi,
     $     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     $     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id$
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
c     
      double precision t1(*),
     &     f1n(*),f1t(*),f2n(*),
     &     f2t(*),f3n(*),f3t(*),f4n(*),f4t(*),eorb(*),
     &     emp4,emp5
      double precision Tij(*), Tkj(*), Tia(*), Tka(*), Xia(*), Xka(*), 
     $     Jia(*), Jka(*), Kia(*), Kka(*),
     $     Jij(*), Jkj(*), Kij(*), Kkj(*), Dja(*), Djka(*), Djia(*)

      integer g_objo,g_objv,ncor,nocc,nvir,iprt,g_coul,
     &     g_exch,oseg_lo,oseg_hi
c     
      double precision eaijk
      integer a,b,c,i,j,k,akold,av,inode,len,next,nxtask
      external nxtask
C     
      Integer Nodes, IAm
c     
      integer klo, khi, start, end
      integer kchunk
c 
#ifdef DEBUG_PRINT
      integer tt
      double precision tt0,tt1,trp_time(26)
#endif
c
      double precision zip
      data zip/0.0d00/
c     
c     apr      call ga_print(g_coul)
c     apr      call ga_print(g_exch)
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
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
                        call dfill(lnvv,zip,f1n,1)
                        call dfill(lnvv,zip,f1t,1)
                        call dfill(lnvv,zip,f2n,1)
                        call dfill(lnvv,zip,f2t,1)
                        call dfill(lnvv,zip,f3n,1)
                        call dfill(lnvv,zip,f3t,1)
                        call dfill(lnvv,zip,f4n,1)
                        call dfill(lnvv,zip,f4t,1)
c     
c     sum(d) (Jia, Kia)bd * Tkj,cd -> Fbc
c     
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
                        call ccsd_dovvv(Jia, Kia,
     $                       Tkj(1+(k-klo)*lnvv),f1n,f2n,f3n,f4n,
     $                       nocc,nvir)
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
                        call ccsd_dovvv(Jka(1+(k-klo)*lnvv),
     $                       Kka(1+(k-klo)*lnvv),
     $                       Tij,f1t,f2t,f3t,f4t,nocc,nvir)
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
                        call ccsd_doooo(Jkj(1+(k-klo)*lnov),
     $                       Kkj(1+(k-klo)*lnov),
     $                       Tia,Xia,
     $                       f1n,f2n,
     $                       f3n,f4n,nocc,nvir)
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
                        call ccsd_doooo(Jij, Kij,
     $                       Tka(1+(k-klo)*lnov),Xka(1+(k-klo)*lnov),
     $                       f1t,f2t,
     $                       f3t,f4t,nocc,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(21) = trp_time(21) + (tt1-tt0)
#endif
c
c                        if (iprt.gt.50)then
c                           call prtfmat(f1n,f1t,f2n,f2t,f3n,f3t,f4n,
c     $                          f4t, nvir)
c                        end if
c
                        eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)
     $                       -eorb(a)
c
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
                        call ccsd_tengy_bgp2(f1n,f1t,f2n,f2t,
     $                                       f3n,f3t,f4n,f4t,
     $                                       Dja(1+(i-1)*nvir),
     $                                       Djia,
     $                                       t1((k-1)*nvir+1),
     $                                       eorb(nocc+ncor+1),
     $                                       eaijk,emp4,emp5,
     $                                       ncor,nocc,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        trp_time(22) = trp_time(22) + (tt1-tt0)
#endif
c     
                        if (i.ne.k)then
#ifdef DEBUG_PRINT
                           tt0 = ga_wtime()
#endif
                           call ccsd_tengy_bgp2(f1t,f1n,f2t,f2n,
     $                                          f3t,f3n,f4t,f4n,
     $                                          Dja(1+(k-1)*nvir),
     $                                          Djka(1+(k-klo)*nvir), 
     $                                          t1((i-1)*nvir+1),
     $                                          eorb(nocc+ncor+1),
     $                                          eaijk,emp4,emp5,
     $                                          ncor,nocc,nvir)
#ifdef DEBUG_PRINT
                           tt1 = ga_wtime()
                           trp_time(23) = trp_time(23) + (tt1-tt0)
#endif
c     
                        end if
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
        !write(6,97) IAm,tt,trp_time(tt)
        trp_time(24) = trp_time(24) + trp_time(tt)
      enddo
      call util_flush(6)
   97 format('node ',i5,': ga_get timer(',i2,') = ',1e15.5)
c
      do tt = 18, 21
        !write(6,98) IAm,tt,trp_time(tt)
        trp_time(25) = trp_time(25) + trp_time(tt)
      enddo
      call util_flush(6)
   98 format('node ',i5,': dgemm timer(',i2,') = ',1e15.5)
c
      do tt = 22, 23
        !write(6,99) IAm,tt,trp_time(tt)
        trp_time(26) = trp_time(26) + trp_time(tt)
      enddo
      call util_flush(6)
   99 format('node ',i5,': tengy timer(',i2,') = ',1e15.5)
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
