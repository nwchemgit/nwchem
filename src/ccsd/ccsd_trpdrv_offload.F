#define ALLOC alloc_if(.true.)  free_if(.false.)
#define FREE  alloc_if(.false.) free_if(.true.)
#define REUSE alloc_if(.false.) free_if(.false.)

      subroutine ccsd_trpdrv_offload(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,emp4,emp5,
     &     oseg_lo,oseg_hi,
     &     kchunk, Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
C     $Id: ccsd_trpdrv_offload.F 26674 2015-01-08 14:36:59Z jhammond $
      implicit none
c
cdir$ ATTRIBUTES OFFLOAD : mic :: dgemm
cdir$ ATTRIBUTES OFFLOAD : mic :: lnov
cdir$ ATTRIBUTES OFFLOAD : mic :: lnvv
cdir$ ATTRIBUTES OFFLOAD : mic :: omp_set_num_threads
cdir$ ATTRIBUTES OFFLOAD : mic :: omp_set_nested
cdir$ ATTRIBUTES OFFLOAD : mic :: mkl_set_num_threads
cdir$ ATTRIBUTES OFFLOAD : mic :: mkl_set_dynamic
cdir$ ATTRIBUTES OFFLOAD : mic :: kmp_set_defaults

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

      double precision Tij(nvir*nvir), Tia(nvir*nocc)
      double precision Xia(nvir*nocc)
      double precision Jia(nvir*nvir), Jij(nvir*nocc)
      double precision Kia(nvir*nvir), Kij(nvir*nocc)
      double precision Djia(nvir)

      double precision dintc1(nvir), dintx1(nvir)
      double precision t1v2(nvir)
      double precision t1v1(nvir*kchunk)
      double precision dintc2(nvir*kchunk)
      double precision dintx2(nvir*kchunk)
      double precision eorb(*)
      double precision Tkj(*), Tka(*)
      double precision Xka(*)
      double precision Jka(*), Jkj(*)
      double precision Kka(*), Kkj(*)
      double precision Dja(*), Djka(*)
! used to make inline threaded tengy correct - for now
c

      double precision emp4,emp5
      double precision emp4z,emp5z,denomz
      double precision emp4iz,emp5iz,emp4kz,emp5kz
      double precision eaijk,ea
      integer g_objo,g_objv,g_coul,g_exch
      integer inode,next,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk,nocc2
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
       integer n_nvir,nc_no1,n_nn,nv_no,klnn,klno,nv_nk
       integer kklnn,kklno
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
CDIR$ ASSUME_ALIGNED f1t: 64
CDIR$ ASSUME_ALIGNED f2n: 64
CDIR$ ASSUME_ALIGNED f2t: 64
CDIR$ ASSUME_ALIGNED f3n: 64
CDIR$ ASSUME_ALIGNED f3t: 64
CDIR$ ASSUME_ALIGNED f4n: 64
CDIR$ ASSUME_ALIGNED f4t: 64
CDIR$ ASSUME_ALIGNED dintc1: 64
CDIR$ ASSUME_ALIGNED dintx1: 64
CDIR$ ASSUME_ALIGNED t1v1: 64
CDIR$ ASSUME_ALIGNED dintc2: 64
CDIR$ ASSUME_ALIGNED dintx2: 64
CDIR$ ASSUME_ALIGNED t1v2: 64
CDIR$ ASSUME_ALIGNED Tij: 64
CDIR$ ASSUME_ALIGNED Tia: 64
CDIR$ ASSUME_ALIGNED Xia: 64
CDIR$ ASSUME_ALIGNED Jia: 64
CDIR$ ASSUME_ALIGNED Jij: 64
CDIR$ ASSUME_ALIGNED Kia: 64
CDIR$ ASSUME_ALIGNED Kij: 64

c
         call omp_set_nested(.true.)
!$omp parallel
!$omp& shared(eorb,kchunk,ea)
!$omp& shared(Tkj,Tka,Xka,Jka,Kka,Jkj,Kkj)
!$omp& shared(ncor,nocc,nvir)
!$omp& shared(emp4z,emp5z)
!$omp& private(n_nvir,n_nn,nv_no,klnn,klno,kklnn,kklno)
!$omp& private(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t)
!$omp& private(Tij,Tia,Xia,Jia,Jij,Kia,Kij)
!$omp& private(Djia,dintc1,dintx1,t1v2,t1v1,dintc2,dintx2)
!$omp& private(i,j,a,klo,k,eaijk,khi)
!$omp single
      nodes = ga_nnodes()
      iam = ga_nodeid()
c
c      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
c
      emp4z=0.0d0
      emp5z=0.0d0
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
      inode=-1
      nocc2=nocc/2
      next=nxtask(nodes, 1)
      n_nvir = nvir*nvir
      nc_no1 = ncor+nocc+1
      nv_no = nvir*nocc
      nv_nk	= nvir*kchunk
      n_nn = ncor+nocc+nvir
      klnn = ((kchunk-1)*lnvv)+n_nvir
      klno = ((kchunk-1)*lnov)+nv_no
!dir$ offload begin target(mic)
         call omp_set_num_threads(64)
         call omp_set_nested(.true.)
         call mkl_set_num_threads(4)
         call mkl_set_dynamic(.false.)
!dir$ end offload

!dir$ offload_transfer target(mic)
     I   nocopy(f1n:length(n_nvir) ALLOC)
     I   nocopy(f1t:length(n_nvir) ALLOC)
     I   nocopy(f2n:length(n_nvir) ALLOC)
     I   nocopy(f2t:length(n_nvir) ALLOC)
     I   nocopy(f3n:length(n_nvir) ALLOC)
     I   nocopy(f3t:length(n_nvir) ALLOC)
     I   nocopy(f4n:length(n_nvir) ALLOC)
     I   nocopy(f4t:length(n_nvir) ALLOC)
     N   in(eorb(ncor+1:n_nn) : ALLOC)
     I   nocopy(Jia:length(n_nvir) ALLOC)
     I   nocopy(Kia:length(n_nvir) ALLOC)
     I   nocopy(Tia:length(nv_no) ALLOC)
     I   nocopy(Xia:length(nv_no) ALLOC)
     I   nocopy(Tij:length(n_nvir) ALLOC)
     I   nocopy(Kij:length(nv_no) ALLOC)
     I   nocopy(Jij:length(nv_no) ALLOC)
     N   nocopy(t1v2:length(nvir) ALLOC)
     N   nocopy(dintc1:length(nvir) ALLOC)
     N   nocopy(dintx1:length(nvir) ALLOC)
     I   nocopy(t1v1:length(nv_nk) ALLOC)
     I   nocopy(dintc2:length(nv_nk) ALLOC)
     I   nocopy(dintx2:length(nv_nk) ALLOC)
     I   nocopy(Tkj:length(klnn)  ALLOC)
     I   nocopy(Kkj:length(klno)  ALLOC)
     I   nocopy(Jkj:length(klno)  ALLOC)
     I   nocopy(Jka:length(klnn)  ALLOC)
     I   nocopy(Tka:length(klno)  ALLOC)
     I   nocopy(Kka:length(klnn)  ALLOC)
     I   nocopy(Xka:length(klno)  ALLOC)
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

      kklnn=((khi-klo)*lnvv)+n_nvir
      kklno=((khi-klo)*lnov)+nv_no
      ea=eorb(a) 
!$omp task  shared( t1,eorb,kchunk,av,
     &  Tkj,Tka,Xka,Jka,Kka,Jkj,Kkj,Dja,Djka,
     &  ncor,nocc,nvir,emp4,emp5,oseg_lo,
     &  g_objo,g_objv,g_coul,g_exch,
     &  nbh_objv1,nbh_objv2,nbh_objv3,nbh_objo1,nbh_objo2,
     &  nbh_objo3,nbh_exch1,nbh_coul1,nbh_objv4)
     &   firstprivate(klo,khi,eaijk)
       call   ccsd_iloop_host(t1,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,emp4,emp5,oseg_lo,
     &     kchunk,Tkj,Tka,Xka,Jka,Kka,Jkj,Kkj, Dja, Djka,
     &     j,a,klo,khi,av,eaijk,
     &     nbh_objv1,nbh_objv2,nbh_objv3,nbh_objo1,nbh_objo2,
     &	  nbh_objo3,nbh_exch1,nbh_coul1,nbh_objv4)
!$omp end task
 !dir$ offload_transfer target(mic)
     I   in(Tkj:length(kklnn)  REUSE)
     I   in(Kkj:length(kklno)  REUSE)
     I   in(Jkj:length(kklno)  REUSE)
     I   in(Jka:length(kklnn)  REUSE)
     I   in(Tka:length(kklno)  REUSE)
     I   in(Kka:length(kklnn)  REUSE)
     I   in(Xka:length(kklno)  REUSE)
                do i=1,nocc2
!$omp  critical
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
!$omp end critical
                     call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     call ga_nbwait(nbh_objv1) ! Dja
                     call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     call dcopy(nvir,Djia,1,dintx1,1)

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
                     do k=klo,min(khi,i)
            call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1((k-klo)*nvir+1),1)
         call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2((k-klo)*nvir+1),1)
                            call ga_nbwait(nbh_objv4(k)) ! Djka	
      call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2((k-klo)*nvir+1),1)	
                     end do

!dir$ offload begin target(mic)  inout(emp4z,emp5z)
     I   in(emp4iz,emp4kz,emp5iz,emp5kz)
     I   in(nvir,nocc,klo,lnvv,lnov,ncor)
     I   in(eaijk,b,c,denomz,k,ea,i,j,khi)
     I   nocopy(Tkj:length(kklnn)  REUSE)
     I   nocopy(Kkj:length(kklno) REUSE)
     I   nocopy(Jkj:length(kklno)  REUSE)
     I   nocopy(Jka:length(kklnn)  REUSE)
     I   nocopy(Tka:length(kklno)  REUSE)
     I   nocopy(Kka:length(kklnn)  REUSE)
     I   nocopy(Xka:length(kklno) REUSE)
     I   nocopy(f1n:length(n_nvir) REUSE)
     I   nocopy(f1t:length(n_nvir) REUSE)
     I   nocopy(f2n:length(n_nvir) REUSE)
     I   nocopy(f2t:length(n_nvir) REUSE)
     I   nocopy(f3n:length(n_nvir) REUSE)
     I   nocopy(f3t:length(n_nvir) REUSE)
     I   nocopy(f4n:length(n_nvir) REUSE)
     I   nocopy(f4t:length(n_nvir) REUSE)
     N   nocopy(eorb(ncor+1:n_nn) : REUSE)
     I   in(t1v1:length(nv_nk) REUSE)
     I   in(dintc2:length(nv_nk) REUSE)
     I   in(dintx2:length(nv_nk) REUSE)
     I   in(Jia:length(n_nvir) REUSE)
     I   in(Kia:length(n_nvir) REUSE)
     I   in(Tia:length(nv_no) REUSE)
     I   in(Xia:length(nv_no) REUSE)
     I   in(Tij:length(n_nvir) REUSE)
     I   in(Kij:length(nv_no) REUSE)
     I   in(Jij:length(nv_no) REUSE)
     N   in(t1v2:length(nvir) REUSE)
     N   in(dintc1:length(nvir) REUSE)
     N   in(dintx1:length(nvir) REUSE)

                     do k=klo,min(khi,i)
                        emp4iz = 0.0d0
                        emp5iz = 0.0d0
                        emp4kz = 0.0d0
                        emp5kz = 0.0d0

                        eaijk=ea - (eorb(ncor+i)
     &                           +eorb(ncor+j)
     &                           +eorb(ncor+k) )
!$omp parallel
!$omp& shared(eorb,eaijk)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t)
!$omp& shared(t1v1,dintc1,dintx1)
!$omp& shared(t1v2,dintc2,dintx2)
!$omp& private(b,c,denomz)
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




!$omp do collapse(2)
!$omp& schedule(static)
!$omp& reduction(+:emp5iz,emp4iz)
!$omp& reduction(+:emp5kz,emp4kz)
               do b=1,nvir
                 do c=1,nvir
                   denomz=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
                   emp4iz=emp4iz+denomz*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
                   emp4iz=emp4iz-denomz*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
                   emp4iz=emp4iz+3*denomz*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp4kz=emp4kz+denomz*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b))*
     &                  (f1n(b,c)-2*f2n(b,c)-2*f3n(b,c)+f4n(b,c))
                   emp4kz=emp4kz-denomz*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b))*
     &                  (2*f1n(b,c)-f2n(b,c)-f3n(b,c)+2*f4n(b,c))
                   emp4kz=emp4kz+3*denomz*(
     &                  f1t(b,c)*(f1t(b,c)+f3t(c,b)+2*f4n(c,b))+
     &                  f2t(b,c)*f2n(c,b)+f3t(b,c)*f4n(b,c))
                   emp5iz=emp5iz+denomz*t1v1((k-klo)*nvir+b)*dintx1(c)*
     &                 (    f1t(b,c)+f2n(b,c)+f4n(c,b)
     &                  -2*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                      f1n(b,c)+f2t(b,c)+f3n(c,b))
     &                  +4*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
                   emp5iz=emp5iz+denomz*t1v1((k-klo)*nvir+b)*dintc1(c)*
     &                 (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &                   -2*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
                   emp5kz=emp5kz+denomz*t1v2(b)*dintx2((k-klo)*nvir+c)*
     &                 (    f1n(b,c)+f2t(b,c)+f4t(c,b)
     &                  -2*(f3n(b,c)+f4t(b,c)+f2t(c,b)+
     &                      f1t(b,c)+f2n(b,c)+f3t(c,b))
     &                  +4*(f3t(b,c)+f4n(b,c)+f1t(c,b)))
                   emp5kz=emp5kz+denomz*t1v2(b)*dintc2((k-klo)*nvir+c)*
     &                 (     f1t(b,c)+f4t(b,c)+f1n(c,b)
     &                   -2*(f2t(b,c)+f3t(b,c)+f2n(c,b)))
                 enddo
               enddo
!$omp end do
!$omp end parallel

                         emp4z = emp4z + emp4iz
                         emp5z = emp5z + emp5iz
                         if (i.ne.k) then
                             emp4z = emp4z + emp4kz
                             emp5z = emp5z + emp5kz
                         end if ! (i.ne.k)
                     end do    ! k
!dir$ end offload
                  end do       ! i
!$omp   taskwait
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
                         emp4 =emp4 + emp4z
                         emp5 = emp5 + emp5z
      call ga_sync()
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
c
!dir$ offload_transfer target(mic)
     I   nocopy(f1n:length(n_nvir) FREE)
     I   nocopy(f1t:length(n_nvir) FREE)
     I   nocopy(f2n:length(n_nvir) FREE)
     I   nocopy(f2t:length(n_nvir) FREE)
     I   nocopy(f3n:length(n_nvir) FREE)
     I   nocopy(f3t:length(n_nvir) FREE)
     I   nocopy(f4n:length(n_nvir) FREE)
     I   nocopy(f4t:length(n_nvir) FREE)
     N   nocopy(eorb(ncor+1:n_nn) : FREE)
     I   nocopy(Jia:length(n_nvir) FREE)
     I   nocopy(Kia:length(n_nvir) FREE)
     I   nocopy(Tia:length(nv_no) FREE)
     I   nocopy(Xia:length(nv_no) FREE)
     I   nocopy(Tij:length(n_nvir) FREE)
     I   nocopy(Kij:length(nv_no) FREE)
     I   nocopy(Jij:length(nv_no) FREE)
     N   nocopy(t1v2:length(nvir) FREE)
     N   nocopy(dintc1:length(nvir) FREE)
     N   nocopy(dintx1:length(nvir) FREE)
     I   nocopy(t1v1:length(nv_nk) FREE)
     I   nocopy(dintc2:length(nv_nk) FREE)
     I   nocopy(dintx2:length(nv_nk) FREE)
     I   nocopy(Tkj:length(klnn) FREE)
     I   nocopy(Kkj:length(klno) FREE)
     I   nocopy(Jkj:length(klno) FREE)
     I   nocopy(Jka:length(klnn) FREE)
     I   nocopy(Tka:length(klno) FREE)
     I   nocopy(Kka:length(klnn) FREE)
     I   nocopy(Xka:length(klno) FREE)

!$omp end single
!$omp end parallel
      end



      subroutine ccsd_iloop_host(t1,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,emp4,emp5,oseg_lo,
     &     kchunk,Tkj,Tka,Xka,Jka,Kka,Jkj,Kkj,Dja,Djka,
     &     j,a,klo,khi,av,eaijk,
     &     nbh_objv1,nbh_objv2,nbh_objv3,nbh_objo1,nbh_objo2,
     &	  nbh_objo3,nbh_exch1,nbh_coul1,nbh_objv4)
C     $Id: ccsd_trpdrv_omp.F 26674 2015-01-08 14:36:59Z jhammond $
      implicit none
c
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
      integer ncor,nocc,nvir
      double precision t1(*)
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision eorb(*)

      double precision Tij(nvir*nvir), Tia(nvir*nocc)
      double precision Xia(nvir*nocc)
      double precision Jia(nvir*nvir), Jij(nvir*nocc)
      double precision Kia(nvir*nvir), Kij(nvir*nocc)
      double precision Djia(nvir)

      double precision Tkj(*), Tka(*)
      double precision Xka(*)
      double precision Jka(*), Jkj(*)
      double precision Kka(*), Kkj(*)
      double precision Dja(*), Djka(*)

      double precision dintc1(nvir),dintx1(nvir),t1v1(nvir)
      double precision dintc2(nvir),dintx2(nvir),t1v2(nvir)
c
      double precision emp4,emp5,denom
      double precision emp4i,emp5i,emp4k,emp5k
      double precision eaijk
      integer g_objo,g_objv,g_coul,g_exch
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,av
      integer klo, khi,nocc2,kchunk



       integer nbh_objv1,nbh_objv2,nbh_objv3
       integer nbh_objv5,nbh_objv6,nbh_objv7
       integer nbh_objv4(nocc)
c
       integer nbh_objo1,nbh_objo2,nbh_objo3
       integer nbh_objo4,nbh_objo5,nbh_objo6
c
       integer nbh_exch1,nbh_exch2,nbh_coul1,nbh_coul2

            nocc2=nocc/2
                  do i=nocc2+1,nocc
!$omp critical
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
!$omp end critical
                     call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     call ga_nbwait(nbh_objv1) ! Dja
                     call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     call dcopy(nvir,Djia,1,dintx1,1)


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

                     do k=klo,min(khi,i)

                        call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1,1)
                        call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2,1)
                        call ga_nbwait(nbh_objv4(k)) ! Djka
                      call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2,1)

                        emp4i = 0.0d0
                        emp5i = 0.0d0
                        emp4k = 0.0d0
                        emp5k = 0.0d0


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

                         emp4 = emp4 + emp4i
                         emp5 = emp5 + emp5i
                         if (i.ne.k) then
                             emp4 = emp4 + emp4k
                             emp5 = emp5 + emp5k
                         end if ! (i.ne.k)
                     end do    ! k
                  end do       ! i
        end
