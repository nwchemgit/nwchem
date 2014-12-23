      subroutine ccsd_trpdrv_omp(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,einv,
     &     eccsdt,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,emp4,emp5,
     &     oseg_lo,oseg_hi,
     &     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id$
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
c     
      integer ncor,nocc,nvir
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
c      integer jump,jump1,jump2
      integer g_objo,g_objv,g_coul,g_exch
      integer inode,next,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk
      integer nxtask
      external nxtask
#ifdef _OPENMP
      integer  omp_get_thread_num
      external omp_get_thread_num
      integer  omp_get_num_threads
      external omp_get_num_threads
      integer  omp_get_max_threads
      external omp_get_max_threads
#endif
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
#ifdef _OPENMP
      if (iam.eq.0) write(6,99) omp_get_max_threads()
   99 format(2x,'Using ',i2,' OpenMP threads in CCSD(T)')
#endif
c
c      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
c
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
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
#ifdef DEBUG_PRINT
                        tt0 = ga_wtime()
#endif
                        if (occsdps) then
                           call pstat_on(ps_doxxx)
                        else
                           call qenter('doxxx',0)
                        endif
c                        jump1 = 1+(k-klo)*lnvv
c                        jump2 = 1+(k-klo)*lnov
!$omp parallel sections
!$omp& firstprivate(nocc,nvir,lnov,lnvv,k,klo)
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f1n,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f3n,nvir)
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f2n,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f4n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f4n,nvir)
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f1t,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f4t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
!$omp end parallel sections
                        if (occsdps) then
                           call pstat_off(ps_doxxx)
                        else
                           call qexit('doxxx',0)
                        endif

#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        dgemm_time = dgemm_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )
c                        jump = ncor+nocc
!$omp parallel do
!$omp& schedule(static)
!$omp& shared(einv,eorb)
!$omp& private(b,c)
!$omp& firstprivate(eaijk)
                        do b=1,nvir
                          do c=1,nvir
                            einv(c,b)=-1.0d0/( eorb(ncor+nocc+b)
     &                                        +eorb(ncor+nocc+c)+eaijk )
                          enddo
                        enddo
!$omp end parallel do

#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        einv_time = einv_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
                        call ccsd_tengy_omp(f1n,f1t,f2n,f2t,
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
                           call ccsd_tengy_omp(f1t,f1n,f2t,f2n,
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

                        end if ! (i.ne.k)
                     end do
                  end do
c                  if (iprt.gt.50)then
c                     write(6,1234)iam,a,j,emp4,emp5
c 1234                format(' iam aijk',3i5,2e15.5)
c                  end if
                  next=nxtask(nodes, 1)
               end if
            end do
            if(ga_nodeid().eq.0) then
               write(6,4321) ' ccsd(t): done ',
     A              a-(ncor+nocc)+((klo-1)/kchunk)*nvir,
     O              ' out of ',(nocc/kchunk)*nvir,
     O              ' progress: ',
     O              ((a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*100d0)/
     D              ((nocc/kchunk)*nvir),
     P              '%'
 4321          format(a,i8,a,i8,a,f6.1,a1)
            endif
         end do
      end do
#ifdef DEBUG_PRINT
      call ga_sync()
      if (iam.eq.0) write(6,87)
   87 format(2x,'node',6x,'ga_get',9x,'dgemm',10x,'einv'10x,'tengy')
      call ga_sync()
      write(6,88) iam,ga_time(1)+ga_time(2),dgemm_time,einv_time,
     &                tengy_time(1)+tengy_time(2)
   88 format(i7,4e15.5)
#endif
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c
      end
