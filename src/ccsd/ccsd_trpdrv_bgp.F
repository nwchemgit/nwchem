      subroutine ccsd_trpdrv_bgp(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,einv,
     &     eccsdt,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,emp4,emp5,
     &     oseg_lo,oseg_hi,
     &     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
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
      double precision tt0,tt1
      double precision ga_time(2)
      double precision tengy_time(2)
      double precision dgemm_time
      double precision einv_time
#endif
      integer jump!,dszvv
      integer g_objo,g_objv,g_coul,g_exch
      integer ncor,nocc,nvir
      integer inode,next,nxtask,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk
      external nxtask
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
c
c      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
!       if (occsdps) then
!          call pstat_on(ps_trpdrv)
!       else
!          call qenter('trpdrv',0)
!       endif
      inode=-1
      next=nxtask(nodes, 1)
#ifdef DEBUG_PRINT
      ga_time(1) = 0.0d0
      ga_time(2) = 0.0d0
      tengy_time(1) = 0.0d0
      tengy_time(2) = 0.0d0
      dgemm_time = 0.0d0
      einv_time = 0.0d0
#endif
      do klo = 1, nocc, kchunk
         akold=0
         khi = min(nocc, klo+kchunk-1)
         do a=oseg_lo,oseg_hi
            av=a-ncor-nocc
            do j=1,nocc
               inode=inode+1
               if (inode.eq.next)then
#ifdef DEBUG_PRINT
                  tt0 = ga_wtime()
#endif
                  call ga_get(g_objv,1+(j-1)*lnov,j*lnov,av,av,Dja,lnov)
                  call ga_get(g_objo,(klo-1)*lnvv+1,khi*lnvv,j,j,Tkj,
     &                 (khi-klo+1)*lnvv)
                  call ga_get(g_objo,lnovv+(klo-1)*lnov+1,
     &                 lnovv+khi*lnov,j,j,Jkj,(khi-klo+1)*lnov)
                  call ga_get(g_objo,lnovv+lnoov+(klo-1)*lnov+1,
     &                 lnovv+lnoov+khi*lnov,j,j,Kkj,(khi-klo+1)*lnov)
                  if (akold .ne. a) then
                     akold = a
                     call ga_get(g_coul,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Jka,lnvv)
                     call ga_get(g_exch,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Kka,lnvv)
                     call ga_get(g_objv,1+lnoov+(klo-1)*lnov,
     &                    lnoov+khi*lnov,av,av,Tka,(khi-klo+1)*lnov)
                     call ga_get(g_objv,1+2*lnoov+(klo-1)*lnov,
     &                    2*lnoov+khi*lnov,av,av,Xka,(khi-klo+1)*lnov)
                  endif
                  do k = klo, khi
                     call ga_get(g_objv,1+(j-1)*nvir+(k-1)*lnov,
     &                    j*nvir+(k-1)*lnov,av,av,
     &                    Djka(1+(k-klo)*nvir),nvir)
                  enddo
#ifdef DEBUG_PRINT
                  tt1 = ga_wtime()
                  ga_time(1) = ga_time(1) + (tt1-tt0)
#endif
                  do i=1,nocc
#ifdef DEBUG_PRINT
                     tt0 = ga_wtime()
#endif
                     call ga_get(g_objo,(i-1)*lnvv+1,i*lnvv,j,j,Tij,
     &                    lnvv)
                     call ga_get(g_objo,lnovv+(i-1)*lnov+1,
     &                    lnovv+i*lnov,j,j,Jij,lnov)
                     call ga_get(g_objo,lnovv+lnoov+(i-1)*lnov+1,
     &                    lnovv+lnoov+i*lnov,j,j,Kij,lnov)
                     call ga_get(g_coul,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Jia,lnvv)
                     call ga_get(g_exch,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Kia,lnvv)
                     call ga_get(g_objv,1+(j-1)*nvir+(i-1)*lnov,
     &                    j*nvir+(i-1)*lnov,av,av,Djia,nvir)
                     call ga_get(g_objv,1+lnoov+(i-1)*lnov,
     &                    lnoov+i*lnov,av,av,Tia,lnov)
                     call ga_get(g_objv,1+2*lnoov+(i-1)*lnov,
     &                    2*lnoov+i*lnov,av,av,Xia,lnov)
#ifdef DEBUG_PRINT
                     tt1 = ga_wtime()
                     ga_time(2) = ga_time(2) + (tt1-tt0)
#endif
                     do k=klo,min(khi,i)
!                        emp4a = 0.0d0
!                        emp5a = 0.0d0
!                        emp4b = 0.0d0
!                        emp5b = 0.0d0
!                         dszvv=8*lnvv ! assuming 8 bytes per double...
!                         call f_memzero(f1n,dszvv)
!                         call f_memzero(f2n,dszvv)
!                         call f_memzero(f3n,dszvv)
!                         call f_memzero(f4n,dszvv)
!                         call f_memzero(f1t,dszvv)
!                         call f_memzero(f2t,dszvv)
!                         call f_memzero(f3t,dszvv)
!                         call f_memzero(f4t,dszvv)
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
                        jump = 1+(k-klo)*lnvv
!                        call ccsd_dovvv_bgp(Jia, Kia,Tkj(1+(k-klo)*lnvv),f1n,f2n,f3n,f4n,nocc,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,Jia,
     1                       nvir,Tkj(jump),nvir,0.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,Jia,
     1                       nvir,Tkj(jump),nvir,0.0d0,
     2                       f3n,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,Kia,
     1                       nvir,Tkj(jump),nvir,0.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,Kia,
     1                       nvir,Tkj(jump),nvir,0.0d0,
     2                       f4n,nvir)
!                        call ccsd_dovvv_bgp(Jka(1+(k-klo)*lnvv),Kka(1+(k-klo)*lnvv),Tij,f1t,f2t,f3t,f4t,nocc,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jka(jump),nvir,Tij,nvir,0.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(jump),nvir,Tij,nvir,0.0d0,
     2                       f3t,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(jump),nvir,Tij,nvir,0.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(jump),nvir,Tij,nvir,0.0d0,
     2                       f4t,nvir)
                        jump = 1+(k-klo)*lnov
!                        call ccsd_doooo_bgp(Jkj(1+(k-klo)*lnov),Kkj(1+(k-klo)*lnov),Tia,Xia,f1n,f2n,f3n,f4n,nocc,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(jump),nocc,1.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(jump),nocc,1.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(jump),nocc,1.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Jkj(jump),nocc,1.0d0,
     2                       f4n,nvir)
!                        call ccsd_doooo_bgp(Jij, Kij,Tka(1+(k-klo)*lnov),Xka(1+(k-klo)*lnov),f1t,f2t,f3t,f4t,nocc,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(jump),nvir,Kij,nocc,1.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(jump),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(jump),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(jump),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        dgemm_time = dgemm_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
                        eaijk=eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k)
     $                       -eorb(a)
                        jump = ncor+nocc
                        do b=1,nvir
                          eorbaux=eorb(jump+b)-eaijk
                          do c=1,nvir
                            einv(c,b)=-1.0d0/(eorbaux+eorb(jump+c))
                          enddo
                        enddo
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        einv_time = einv_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
!                        call ccsd_tengy4_bgp(f1n,f1t,f2n,f2t,
!     1                                       f3n,f3t,f4n,f4t,
!     2                                       einv,emp4a,nvir)
!                        call ccsd_tengy4ab_bgp(f1n,f1t,f2n,f2t,
!     1                                         f3n,f3t,f4n,f4t,
!     2                                         einv,emp4a,emp4b,nvir)
!                        call ccsd_tengy5_bgp(f1n,f1t,f2n,f2t,
!     1                                       f3n,f3t,f4n,f4t,
!     2                                       Dja(1+(i-1)*nvir),
!     3                                       Djia,
!     4                                       t1((k-1)*nvir+1),
!     5                                       einv,emp5a,nvir)
!                        call ccsd_tengy5ab_bgp(f1n,f1t,f2n,f2t,
!     1                                         f3n,f3t,f4n,f4t,
!     2                                         Dja(1+(i-1)*nvir),    ! A
!     3                                         Djia,                 ! A
!     4                                         t1((k-1)*nvir+1),     ! A
!     5                                         Dja(1+(k-1)*nvir),    ! B
!     6                                         Djka(1+(k-klo)*nvir), ! B 
!     7                                         t1((i-1)*nvir+1),     ! B
!     8                                         einv,emp5a,emp5b,nvir)
                        call ccsd_tengy_bgp(f1n,f1t,f2n,f2t,
     1                                      f3n,f3t,f4n,f4t,
     2                                      Dja(1+(i-1)*nvir),
     3                                      Djia,
     4                                      t1((k-1)*nvir+1),
     5                                      einv,emp4,emp5,
     6                                      nvir)
#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        tengy_time(1) = tengy_time(1) + (tt1-tt0)
#endif
                        if (i.ne.k) then
#ifdef DEBUG_PRINT
                           tt0 = ga_wtime()
#endif
!                           call ccsd_tengy4_bgp(f1t,f1n,f2t,f2n,
!     1                                          f3t,f3n,f4t,f4n,
!     2                                          einv,emp4b,nvir)
!                           call ccsd_tengy5_bgp(f1t,f1n,f2t,f2n,
!     1                                          f3t,f3n,f4t,f4n,
!     2                                          Dja(1+(k-1)*nvir),
!     3                                          Djka(1+(k-klo)*nvir), 
!     4                                          t1((i-1)*nvir+1),
!     5                                          einv,emp5b,nvir)
                           call ccsd_tengy_bgp(f1t,f1n,f2t,f2n,
     1                                         f3t,f3n,f4t,f4n,
     2                                         Dja(1+(k-1)*nvir),
     3                                         Djka(1+(k-klo)*nvir),
     4                                         t1((i-1)*nvir+1),
     5                                         einv,emp4,emp5,
     6                                         nvir)
#ifdef DEBUG_PRINT
                           tt1 = ga_wtime()
                           tengy_time(2) = tengy_time(2) + (tt1-tt0)
#endif
!                        else
!                        if (i.eq.k) then
!                           emp4b = 0.0d0
!                           emp5b = 0.0d0
!                        end if ! (i.eq.k)
                        end if ! (i.ne.k)
!                        emp4 = emp4 + emp4a + emp4b
!                        emp5 = emp5 + emp5a + emp5b
                     end do
                  end do
                  next=nxtask(nodes, 1)
               end if
            end do
         end do
      end do
#ifdef DEBUG_PRINT
      call ga_sync()
      if (IAm.eq.0) write(6,87)
   87 format(2x,'node',6x,'ga_get',9x,'dgemm',10x,'einv'10x,'tengy')
      call ga_sync()
      write(6,88) IAm,ga_time(1)+ga_time(2),dgemm_time,einv_time,
     &                tengy_time(1)+tengy_time(2)
   88 format(i7,4e15.5)
#endif
      next=nxtask(-nodes, 1)
!       call ga_sync()
!       if (occsdps) then
!          call pstat_off(ps_trpdrv)
!       else
!          call qexit('trpdrv',0)
!       endif
c
      end


      subroutine ccsd_tengy_bgp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     1                          dintc,dintx,t1v,einv,
     2                          emp4,emp5,nvir)
      implicit none
      integer nvir,b,c
      double precision emp4,emp5,denom
      double precision two,three,four
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision dintc(nvir),dintx(nvir),t1v(nvir)
      double precision einv(nvir,nvir)
      data two/2.0d00/,three/3.0d00/,four/4.0d00/
c
      do b=1,nvir
        do c=1,nvir
          denom=einv(c,b)
c
          emp4=emp4+denom*
     &         (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &         (f1t(b,c)-two*f2t(b,c)-two*f3t(b,c)+f4t(b,c))
c
          emp4=emp4-denom*
     &         (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &         (two*f1t(b,c)-f2t(b,c)-f3t(b,c)+two*f4t(b,c))
c
          emp4=emp4+three*denom*(
     &         f1n(b,c)*(f1n(b,c)+f3n(c,b)+two*f4t(c,b))+
     &         f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
c
          emp5=emp5+denom*t1v(b)*dintx(c)*
     &        (     f1t(b,c)+f2n(b,c)+f4n(c,b)
     &         -two*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &              f1n(b,c)+f2t(b,c)+f3n(c,b))
     &         +four*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
          emp5=emp5+denom*t1v(b)*dintc(c)*
     &        (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &          -two*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
c
        enddo
      enddo
c
      return
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

