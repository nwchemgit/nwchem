      subroutine ccsd_trpdrv_omp(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
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
      double precision eorb(*)
      double precision Tij(*), Tkj(*), Tia(*), Tka(*)
      double precision Xia(*), Xka(*)
      double precision Jia(*), Jka(*), Jij(*), Jkj(*)
      double precision Kia(*), Kka(*), Kij(*), Kkj(*)
      double precision Dja(*), Djka(*), Djia(*)
! used to make inline threaded tengy correct - for now
      double precision dintc1(nvir),dintx1(nvir),t1v1(nvir)
      double precision dintc2(nvir),dintx2(nvir),t1v2(nvir)
c
      double precision emp4,emp5,denom
      double precision emp4i,emp5i,emp4k,emp5k
      double precision eaijk
      integer g_objo,g_objv,g_coul,g_exch
      integer inode,next,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk
      integer nxtask
      external nxtask
c
c  Dependencies (global array, local array, handle): 
c
c  These are waited on first
c
c      g_objv, Dja,  nbh_objv1
c      g_objv, Djka(1+(k-klo)*nvir), nbh_objv4(k)
c      g_objv, Djia, nbh_objv5
c
c  These are waited on later
c
c      g_objv, Tka,  nbh_objv2
c      g_objv, Xka,  nbh_objv3
c      g_objv, Tia,  nbh_objv6
c      g_objv, Xia,  nbh_objv7
c      g_objo, Tkj,  nbh_objo1
c      g_objo, Jkj,  nbh_objo2
c      g_objo, Kkj,  nbh_objo3
c      g_objo, Tij,  nbh_objo4
c      g_objo, Jij,  nbh_objo5
c      g_objo, Kij,  nbh_objo6
c      g_exch, Kka,  nbh_exch1
c      g_exch, Kia,  nbh_exch2
c      g_coul, Jka,  nbh_coul1
c      g_coul, Jia,  nbh_coul2
c
c  non-blocking handles
c
       integer nbh_objv1,nbh_objv2,nbh_objv3
       integer nbh_objv5,nbh_objv6,nbh_objv7
       integer nbh_objv4(nocc)
c
       integer nbh_objo1,nbh_objo2,nbh_objo3
       integer nbh_objo4,nbh_objo5,nbh_objo6
c
       integer nbh_exch1,nbh_exch2,nbh_coul1,nbh_coul2
c
#ifdef _OPENMP
      integer  omp_get_thread_num
      external omp_get_thread_num
      integer  omp_get_num_threads
      external omp_get_num_threads
      integer  omp_get_max_threads
      external omp_get_max_threads
      if (ga_nodeid().eq.0) write(6,99) omp_get_max_threads()
   99 format(2x,'Using ',i2,' OpenMP threads in CCSD(T)')
#endif
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
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
      do klo = 1, nocc, kchunk
         akold=0
         khi = min(nocc, klo+kchunk-1)
         do a=oseg_lo,oseg_hi
            av=a-ncor-nocc
            do j=1,nocc
               inode=inode+1
               if (inode.eq.next)then

                  call ga_nbget(g_objv,1+(j-1)*lnov,j*lnov,av,av,Dja,
     &                          lnov,nbh_objv1)
                  do k = klo, khi
                     call ga_nbget(g_objv,1+(j-1)*nvir+(k-1)*lnov,
     &                    j*nvir+(k-1)*lnov,av,av,
     &                    Djka(1+(k-klo)*nvir),nvir,nbh_objv4(k))
                  enddo
                  call ga_nbget(g_objo,(klo-1)*lnvv+1,khi*lnvv,j,j,Tkj,
     &                          (khi-klo+1)*lnvv,nbh_objo1)
                  call ga_nbget(g_objo,lnovv+(klo-1)*lnov+1,
     &                          lnovv+khi*lnov,j,j,Jkj,
     &                          (khi-klo+1)*lnov,nbh_objo2)
                  call ga_nbget(g_objo,lnovv+lnoov+(klo-1)*lnov+1,
     &                          lnovv+lnoov+khi*lnov,j,j,Kkj,
     &                          (khi-klo+1)*lnov,nbh_objo3)
                  if (akold .ne. a) then
                     akold = a
                     call ga_nbget(g_coul,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Jka,lnvv,nbh_coul1)
                     call ga_nbget(g_exch,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Kka,lnvv,nbh_exch1)
                     call ga_nbget(g_objv,1+lnoov+(klo-1)*lnov,
     &                    lnoov+khi*lnov,av,av,Tka,(khi-klo+1)*lnov,
     &                    nbh_objv2)
                     call ga_nbget(g_objv,1+2*lnoov+(klo-1)*lnov,
     &                    2*lnoov+khi*lnov,av,av,Xka,(khi-klo+1)*lnov,
     &                    nbh_objv3)
                  endif

                  do i=1,nocc

                     call ga_nbget(g_objv,1+(j-1)*nvir+(i-1)*lnov,
     &                    j*nvir+(i-1)*lnov,av,av,Djia,nvir,nbh_objv5)
                     call ga_nbget(g_objo,(i-1)*lnvv+1,i*lnvv,j,j,Tij,
     &                    lnvv,nbh_objo4)
                     call ga_nbget(g_objo,lnovv+(i-1)*lnov+1,
     &                    lnovv+i*lnov,j,j,Jij,lnov,nbh_objo5)
                     call ga_nbget(g_objo,lnovv+lnoov+(i-1)*lnov+1,
     &                    lnovv+lnoov+i*lnov,j,j,Kij,lnov,nbh_objo6)
                     call ga_nbget(g_coul,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Jia,lnvv,nbh_coul2)
                     call ga_nbget(g_exch,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Kia,lnvv,nbh_exch2)
                     call ga_nbget(g_objv,1+lnoov+(i-1)*lnov,
     &                    lnoov+i*lnov,av,av,Tia,lnov,nbh_objv6)
                     call ga_nbget(g_objv,1+2*lnoov+(i-1)*lnov,
     &                    2*lnoov+i*lnov,av,av,Xia,lnov,nbh_objv7)

                     call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     call ga_nbwait(nbh_objv1) ! Dja
                     call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     call dcopy(nvir,Djia,1,dintx1,1)

                     do k=klo,min(khi,i)

                        call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1,1)
                        call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2,1)
                        call ga_nbwait(nbh_objv4(k)) ! Djka
                        call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2,1)
                        emp4i = 0.0d0
                        emp5i = 0.0d0
                        emp4k = 0.0d0
                        emp5k = 0.0d0
                        if (occsdps) then
                           call pstat_on(ps_doxxx)
                        else
                           call qenter('doxxx',0)
                        endif
!
!  These are the input dependencies for the DGEMM calls below.  
!  We wait on all of them here because GA is not even remotely thread-safe.
!
                        call ga_nbwait(nbh_objv2)
                        call ga_nbwait(nbh_objv3)
                        call ga_nbwait(nbh_objv6)
                        call ga_nbwait(nbh_objv7)
                        call ga_nbwait(nbh_objo1)
                        call ga_nbwait(nbh_objo2)
                        call ga_nbwait(nbh_objo3)
                        call ga_nbwait(nbh_objo4)
                        call ga_nbwait(nbh_objo5)
                        call ga_nbwait(nbh_objo6)
                        call ga_nbwait(nbh_exch1)
                        call ga_nbwait(nbh_exch2)
                        call ga_nbwait(nbh_coul1)
                        call ga_nbwait(nbh_coul2)

!$omp parallel
!$omp& shared(eorb)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t)
!$omp& shared(t1v1,dintc1,dintx1)
!$omp& shared(t1v2,dintc2,dintx2)
!$omp& private(b,c,eaijk,denom)
!$omp& firstprivate(ncor,nocc,nvir,lnov,lnvv,i,j,k,klo)

!$omp sections
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f1n,nvir)
!$omp section
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f2n,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f3n,nvir)
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
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
!$omp section
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f4t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
!$omp end sections

!$omp master
                        if (occsdps) then
                           call pstat_off(ps_doxxx)
                        else
                           call qexit('doxxx',0)
                        endif
                        if (occsdps) then
                           call pstat_on(ps_tengy)
                        else
                           call qenter('tengy',0)
                        endif
!$omp end master

                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )
!$omp do collapse(2)
!$omp& schedule(static)
!$omp& reduction(+:emp5i,emp4i)
!$omp& reduction(+:emp5k,emp4k)
               do b=1,nvir
                 do c=1,nvir
                   denom=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
                   emp4i=emp4i+denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
                   emp4i=emp4i-denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
                   emp4i=emp4i+3*denom*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp4k=emp4k+denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b))*
     &                  (f1n(b,c)-2*f2n(b,c)-2*f3n(b,c)+f4n(b,c))
                   emp4k=emp4k-denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b))*
     &                  (2*f1n(b,c)-f2n(b,c)-f3n(b,c)+2*f4n(b,c))
                   emp4k=emp4k+3*denom*(
     &                  f1t(b,c)*(f1t(b,c)+f3t(c,b)+2*f4n(c,b))+
     &                  f2t(b,c)*f2n(c,b)+f3t(b,c)*f4n(b,c))
                   emp5i=emp5i+denom*t1v1(b)*dintx1(c)*
     &                 (    f1t(b,c)+f2n(b,c)+f4n(c,b)
     &                  -2*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                      f1n(b,c)+f2t(b,c)+f3n(c,b))
     &                  +4*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
                   emp5i=emp5i+denom*t1v1(b)*dintc1(c)*
     &                 (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &                   -2*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
                   emp5k=emp5k+denom*t1v2(b)*dintx2(c)*
     &                 (    f1n(b,c)+f2t(b,c)+f4t(c,b)
     &                  -2*(f3n(b,c)+f4t(b,c)+f2t(c,b)+
     &                      f1t(b,c)+f2n(b,c)+f3t(c,b))
     &                  +4*(f3t(b,c)+f4n(b,c)+f1t(c,b)))
                   emp5k=emp5k+denom*t1v2(b)*dintc2(c)*
     &                 (     f1t(b,c)+f4t(b,c)+f1n(c,b)
     &                   -2*(f2t(b,c)+f3t(b,c)+f2n(c,b)))
                 enddo
               enddo
!$omp end do
!$omp end parallel
                         if (occsdps) then
                            call pstat_off(ps_tengy)
                         else
                            call qexit('tengy',0)
                         endif

                         emp4 = emp4 + emp4i
                         emp5 = emp5 + emp5i
                         if (i.ne.k) then
                             emp4 = emp4 + emp4k
                             emp5 = emp5 + emp5k
                         end if ! (i.ne.k)
                     end do    ! k
                  end do       ! i
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
      call ga_sync()
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c
      end
