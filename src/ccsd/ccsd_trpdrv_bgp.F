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
      integer jump,jump1,jump2
      integer g_objo,g_objv,g_coul,g_exch
      integer inode,next,nodes,iam
      integer oseg_lo,oseg_hi
      integer a,b,c,i,j,k,akold,av
      integer klo, khi, kchunk
      integer nxtask
      external nxtask
ceric      integer omp_get_max_threads
ceric      external omp_get_max_threads
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
ceric      if (iam.eq.0) write(6,99) omp_get_max_threads()
ceric   99 format(2x,'Using ',i2,' OpenMP threads in CCSD(T)')
c
c      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
c
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

!$omp sections
!$omp& firstprivate(nvir)
!$omp& private(jump1,jump2)
!$omp section
                        jump1 = 1+(k-klo)*lnvv
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(jump1),nvir,0.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(jump1),nvir,0.0d0,
     2                       f3n,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(jump1),nvir,0.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(jump1),nvir,0.0d0,
     2                       f4n,nvir)
!$omp section
                        jump1 = 1+(k-klo)*lnvv
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jka(jump1),nvir,Tij,nvir,0.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(jump1),nvir,Tij,nvir,0.0d0,
     2                       f3t,nvir)
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(jump1),nvir,Tij,nvir,0.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(jump1),nvir,Tij,nvir,0.0d0,
     2                       f4t,nvir)
!$omp section
                        jump2 = 1+(k-klo)*lnov
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(jump2),nocc,1.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(jump2),nocc,1.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(jump2),nocc,1.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Jkj(jump2),nocc,1.0d0,
     2                       f4n,nvir)
!$omp section
                        jump2 = 1+(k-klo)*lnov
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(jump2),nvir,Kij,nocc,1.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(jump2),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(jump2),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(jump2),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
!$omp end sections

#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        dgemm_time = dgemm_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )
                        jump = ncor+nocc
!$omp parallel do
!$omp& schedule(static)
!$omp& shared(einv,eorb)
!$omp& private(b,c)
!$omp& firstprivate(eaijk,jump)
                        do b=1,nvir
                          do c=1,nvir
                            einv(c,b)=-1.0d0/( eorb(jump+b)
     &                                        +eorb(jump+c)+eaijk )
                          enddo
                        enddo
!$omp end parallel do

#ifdef DEBUG_PRINT
                        tt1 = ga_wtime()
                        einv_time = einv_time + (tt1-tt0)
                        tt0 = ga_wtime()
#endif
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

                        end if ! (i.ne.k)
                     end do
                  end do
                  next=nxtask(nodes, 1)
               end if
            end do
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
c
      end

#ifdef NO_VITALI_TENGY
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
c
          emp5=emp5+denom*t1v(b)*dintc(c)*
     &        (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &          -two*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
c
        enddo
      enddo
c
      return
      end
#else
!
!    Argonne Leadership Computing Facility
!    BlueGene/P version
!    Written by Vitali Morozov <morozov@anl.gov>
!    Updated 20091022
!    -O5 -qipa=noinline -qsmp=omp
!
!    modified from the original (use of einv)
!
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
      double precision e1, e2, e3, e4
      double precision t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
      double precision z0, z1, z3, z5, z9, s0
c
!$omp parallel do
!$omp& private(b,c,e1,e2,e3,e4,z0,z1,z3,z5,z9,s0,
!$omp&         t0,t1,t2,t3,t4,t5,t6,t7,t8,t9)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t,dintx,dintc,t1v,einv)
!$omp& reduction(+:emp5,emp4) schedule(static)
      do b = 1, nvir
         do c = 1, nvir

            s0 = einv(c,b)

            e1 = f1n(c,b) + f4n(c,b) - 2d0*( f2n(c,b) + f3n(c,b) )
            e2 = e1 + 3d0*f1n(c,b)

            e3 = f1t(c,b) -2d0 * f2t(c,b)
            e4 = e3 + f2n(c,b) - 2d0*( f3t(c,b) + f4n(c,b) + f1n(c,b)
     1                                - 2d0*( f3n(c,b) + f4t(c,b) ) )

            emp5 = emp5 + ( t1v(c)*( e1*dintc(b) + e4*dintx(b) )
     1                    + t1v(b)*( e2*dintx(c) + e3*dintc(c) ) ) * s0

            t1 = f1t(c,b) + f3n(c,b)
            t2 = f1n(c,b) + f2t(c,b) + f4n(c,b)
            t9 = f1t(c,b) + f4t(c,b)
            t0 = f2t(c,b) + f3t(c,b)
            t3 = t9 - 2d0*t0
            t4 = t1 + f2n(c,b)
            t5 = 2d0*t9 - t0
            t6 = f3n(c,b) + 2d0*f4t(c,b)
            t7 = f3n(c,b) * f4t(c,b)
            t8 = t1*t3 - f1n(c,b)*t5 +3d0*( f1n(c,b)*f1n(c,b) + t7 )

            z0 = f2t(b,c) + f3t(b,c)
            z9 = f1t(b,c) + f4t(b,c)
            z1 = f1t(b,c) + f3n(b,c)
            z3 = z9 - 2d0*z0
            z5 = 2d0*z9 - z0

            emp4 = emp4 + ( t8 + t2*z3 - t4*z5 
     1           + 3d0*( t6*f1n(b,c) + f2t(c,b)*f2n(b,c) ) ) * s0
         enddo
      enddo
c
      return
      end
#endif


#ifdef DEAD_CODE
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
#endif
