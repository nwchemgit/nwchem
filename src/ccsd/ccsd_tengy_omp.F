#ifndef BGP
      subroutine ccsd_tengy_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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
!$omp parallel do collapse(2)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t,dintx,dintc,t1v,einv)
!$omp& reduction(+:emp5,emp4) schedule(static)
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
!$omp end parallel do
c
      return
      end
#else /* BGP */
!
!    Argonne Leadership Computing Facility
!    BlueGene/P version
!    Written by Vitali Morozov <morozov@anl.gov>
!    Updated 20091022
!    -O5 -qipa=noinline -qsmp=omp
!
!    modified from the original (use of einv)
!
      subroutine ccsd_tengy_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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
#endif /* BGP */


#ifdef DEAD_CODE
      subroutine ccsd_tengy4_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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

      subroutine ccsd_tengy5_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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



      subroutine ccsd_tengy4ab_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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



      subroutine ccsd_tengy5ab_omp(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
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
