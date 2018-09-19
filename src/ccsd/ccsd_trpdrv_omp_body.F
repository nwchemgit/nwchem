      subroutine ccsd_trpdrv_omp_body(
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     ncor,nocc,nvir, emp4,emp5,
     &     a,i,j,k,klo,
     &     Tij, Tkj, Tia, Tka, Xia, Xka,
     &     Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj,
     &     dintc1,dintx1,t1v1,
     &     dintc2,dintx2,t1v2)
      implicit none
      double precision, intent(inout) :: emp4,emp5
      integer, intent(in) :: ncor,nocc,nvir
      integer, intent(in) :: a,i,j,k, klo
      double precision, intent(in) :: f1n(nvir,nvir),f1t(nvir,nvir)
      double precision, intent(in) :: f2n(nvir,nvir),f2t(nvir,nvir)
      double precision, intent(in) :: f3n(nvir,nvir),f3t(nvir,nvir)
      double precision, intent(in) :: f4n(nvir,nvir),f4t(nvir,nvir)
      double precision, intent(in) :: eorb(*)
      double precision, intent(in) :: Tij(*), Tkj(*), Tia(*), Tka(*)
      double precision, intent(in) :: Xia(*), Xka(*)
      double precision, intent(in) :: Jia(*), Jka(*), Jij(*), Jkj(*)
      double precision, intent(in) :: Kia(*), Kka(*), Kij(*), Kkj(*)
      double precision, intent(in) :: dintc1(nvir),dintx1(nvir)
      double precision, intent(in) :: dintc2(nvir),dintx2(nvir)
      double precision, intent(in) :: t1v1(nvir),t1v2(nvir)
      double precision :: emp4i,emp5i,emp4k,emp5k
      double precision :: eaijk,denom
      integer :: lnov,lnvv
      ! chunking is the loop blocking size in the loop nest
      ! formerly associated with the tengy routine.
      ! we have not explored this paramater space but 32 is
      ! optimal for TLB blocking in matrix transpose on most
      ! architectures (especially x86).
      integer, parameter :: chunking = 32
      integer :: b,c,bb,cc
      lnov=nocc*nvir
      lnvv=nvir*nvir
      emp4i = 0.0d0
      emp5i = 0.0d0
      emp4k = 0.0d0
      emp5k = 0.0d0
!$omp parallel
!$omp& shared(eorb)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t)
!$omp& shared(t1v1,dintc1,dintx1)
!$omp& shared(t1v2,dintc2,dintx2)
!$omp& private(eaijk,denom)
!$omp& firstprivate(ncor,nocc,nvir,lnov,lnvv,i,j,k,klo)
!
! Performance Note:
!
! By definition, the following does not scale to more than 8 threads
! unless nested parallelism (i.e. inside of DGEMM) is used.
! It may be prudent to write a manually threaded wrapper for the
! cases where single-threaded BLAS is used.
!
!$omp sections
!$omp section
      call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1     Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2     f1n,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2     f1n,nvir)
!$omp section
      call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1     Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2     f2n,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2     f2n,nvir)
!$omp section
      call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1     Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2     f3n,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2     f3n,nvir)
!$omp section
      call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1     Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2     f4n,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Xia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2     f4n,nvir)
!$omp section
      call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1     Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2     f1t,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Tka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2     f1t,nvir)
!$omp section
      call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1     Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2     f2t,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2     f2t,nvir)
!$omp section
      call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1     Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2     f3t,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2     f3t,nvir)
!$omp section
      call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1     Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2     f4t,nvir)
      call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1     Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2     f4t,nvir)
!$omp end sections

      eaijk = eorb(a) - ( eorb(ncor+i)+eorb(ncor+j)+eorb(ncor+k) )

!$omp do collapse(2)
!$omp& schedule(static)
!$omp& reduction(+:emp5i,emp4i)
!$omp& reduction(+:emp5k,emp4k)
! WARNING: Do not add IVDEP here.  Code will be incorrect.
           do bb=1,nvir,chunking
             do cc=1,nvir,chunking
               do b=bb,min(bb+chunking-1,nvir)
                 do c=cc,min(cc+chunking-1,nvir)
                   denom=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
! fusing emp[45][ki] accumulates may help vectorization...
                   emp4i=emp4i+denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
     &                        -denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
     &                        +3*denom*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp4k=emp4k+denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b))*
     &                  (f1n(b,c)-2*f2n(b,c)-2*f3n(b,c)+f4n(b,c))
     &                        -denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b))*
     &                  (2*f1n(b,c)-f2n(b,c)-f3n(b,c)+2*f4n(b,c))
     &                        +3*denom*(
     &                  f1t(b,c)*(f1t(b,c)+f3t(c,b)+2*f4n(c,b))+
     &                  f2t(b,c)*f2n(c,b)+f3t(b,c)*f4n(b,c))
                   emp5i=emp5i+denom*t1v1(b)*dintx1(c)*
     &                 (    f1t(b,c)+f2n(b,c)+f4n(c,b)
     &                  -2*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                      f1n(b,c)+f2t(b,c)+f3n(c,b))
     &                  +4*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
     &                        +denom*t1v1(b)*dintc1(c)*
     &                 (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &                   -2*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
                   emp5k=emp5k+denom*t1v2(b)*dintx2(c)*
     &                 (    f1n(b,c)+f2t(b,c)+f4t(c,b)
     &                  -2*(f3n(b,c)+f4t(b,c)+f2t(c,b)+
     &                      f1t(b,c)+f2n(b,c)+f3t(c,b))
     &                  +4*(f3t(b,c)+f4n(b,c)+f1t(c,b)))
     &                        +denom*t1v2(b)*dintc2(c)*
     &                 (     f1t(b,c)+f4t(b,c)+f1n(c,b)
     &                   -2*(f2t(b,c)+f3t(b,c)+f2n(c,b)))
                 enddo
               enddo
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
      end
