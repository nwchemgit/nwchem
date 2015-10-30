!
! written by Jeff Hammond, Intel Labs, June 2014
! see https://github.com/jeffhammond/nwchem-tce-triples-kernels/blob/master/src/new_sort4.f for history
!
! There is currently no difference between the OpenMP and non-OpenMP versions
! except for the directives, but preprocessing at fine granularity leads to ugly
! code because Fortran does not have anything like C99/C++11 Pragma.
!
#ifdef USE_OPENMP
      subroutine tce_sort_4(u,s,a,b,c,d,i,j,k,l,f)
      implicit none
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision s(a*b*c*d)
      double precision u(a*b*c*d)
      double precision f
      integer iax,ibx
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      if (abs(f-1.0).lt.1.e-12) then
      if (l.eq.4) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j3 = 1,c
          id(1) = j1
          id(2) = j2
          id(3) = j3
          !id(4) = j4 ! 1234 are only used once and we know l=4, hence ijk must not be 4
          iax = d*(j3-1+c*(j2-1+b*(j1-1)))
          ibx = d*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j4 = 1,d
           ia = j4+iax
           ib = j4+ibx
           s(ib) = u(ia)
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.3) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j4 = 1,d
          id(1) = j1
          id(2) = j2
          !id(3) = j3
          id(4) = j4
          iax = -1+c*(j2-1+b*(j1-1))
          ibx = c*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j3 = 1,c
           ia = j4+d*(j3+iax)
           ib = j3+ibx
           s(ib) = u(ia)
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.2) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j3 = 1,c
         do j4 = 1,d
          id(1) = j1
          !id(2) = j2
          id(3) = j3
          id(4) = j4
          iax = -1+b*(j1-1)
          ibx = b*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j2 = 1,b
           ia = j4+d*(j3-1+c*(j2+iax))
           ib = j2+ibx
           s(ib) = u(ia)
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else ! l.eq.1
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          !id(1) = j1
          id(2) = j2
          id(3) = j3
          id(4) = j4
          ibx = a*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j1 = 1,a
           ia = j4+d*(j3-1+c*(j2-1+b*(j1-1)))
           ib = j1+ibx
           s(ib) = u(ia)
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      endif
      else ! general case: f=/=1.0
      if (l.eq.4) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j3 = 1,c
          id(1) = j1
          id(2) = j2
          id(3) = j3
          !id(4) = j4 ! 1234 are only used once and we know l=4, hence ijk must not be 4
          iax = d*(j3-1+c*(j2-1+b*(j1-1)))
          ibx = d*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j4 = 1,d
           ia = j4+iax
           ib = j4+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.3) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j4 = 1,d
          id(1) = j1
          id(2) = j2
          !id(3) = j3
          id(4) = j4
          iax = -1+c*(j2-1+b*(j1-1))
          ibx = c*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j3 = 1,c
           ia = j4+d*(j3+iax)
           ib = j3+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.2) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j3 = 1,c
         do j4 = 1,d
          id(1) = j1
          !id(2) = j2
          id(3) = j3
          id(4) = j4
          iax = -1+b*(j1-1)
          ibx = b*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j2 = 1,b
           ia = j4+d*(j3-1+c*(j2+iax))
           ib = j2+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else ! l.eq.1
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          !id(1) = j1
          id(2) = j2
          id(3) = j3
          id(4) = j4
          ibx = a*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j1 = 1,a
           ia = j4+d*(j3-1+c*(j2-1+b*(j1-1)))
           ib = j1+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      endif
      endif
      return
      end

      subroutine tce_sortacc_4(u,s,a,b,c,d,i,j,k,l,f)
      implicit none
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision s(a*b*c*d)
      double precision u(a*b*c*d)
      double precision f
      integer iax,ibx
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      if (l.eq.4) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j3 = 1,c
          id(1) = j1
          id(2) = j2
          id(3) = j3
          !id(4) = j4 ! 1234 are only used once and we know l=4, hence ijk must not be 4
          iax = d*(j3-1+c*(j2-1+b*(j1-1)))
          ibx = d*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j4 = 1,d
           ia = j4+iax
           ib = j4+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.3) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j2 = 1,b
         do j4 = 1,d
          id(1) = j1
          id(2) = j2
          !id(3) = j3
          id(4) = j4
          iax = -1+c*(j2-1+b*(j1-1))
          ibx = c*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j3 = 1,c
           ia = j4+d*(j3+iax)
           ib = j3+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else if (l.eq.2) then
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j1 = 1,a
        do j3 = 1,c
         do j4 = 1,d
          id(1) = j1
          !id(2) = j2
          id(3) = j3
          id(4) = j4
          iax = -1+b*(j1-1)
          ibx = b*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j2 = 1,b
           ia = j4+d*(j3-1+c*(j2+iax))
           ib = j2+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      else ! l.eq.1
!$omp  parallel do collapse(3)
!$omp& private(j1,j2,j3,j4,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,i,j,k,l,f,jd)
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          !id(1) = j1
          id(2) = j2
          id(3) = j3
          id(4) = j4
          ibx = a*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
!dec$ vector always nontemporal
!$omp simd
          do j1 = 1,a
           ia = j4+d*(j3-1+c*(j2-1+b*(j1-1)))
           ib = j1+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do
      endif
      return
      end
#else
      subroutine tce_sort_4(u,s,a,b,c,d,i,j,k,l,f)
      implicit none
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision s(a*b*c*d)
      double precision u(a*b*c*d)
      double precision f
      integer iax,ibx
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      if (l.eq.4) then
       do j1 = 1,a
        do j2 = 1,b
         do j3 = 1,c
          id(1) = j1
          id(2) = j2
          id(3) = j3
          !id(4) = j4 ! 1234 are only used once and we know l=4, hence ijk must not be 4
          iax = d*(j3-1+c*(j2-1+b*(j1-1)))
          ibx = d*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j4 = 1,d
           ia = j4+iax
           ib = j4+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else if (l.eq.3) then
       do j1 = 1,a
        do j2 = 1,b
         do j4 = 1,d
          id(1) = j1
          id(2) = j2
          !id(3) = j3
          id(4) = j4
          iax = -1+c*(j2-1+b*(j1-1))
          ibx = c*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j3 = 1,c
           ia = j4+d*(j3+iax)
           ib = j3+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else if (l.eq.2) then
       do j1 = 1,a
        do j3 = 1,c
         do j4 = 1,d
          id(1) = j1
          !id(2) = j2
          id(3) = j3
          id(4) = j4
          iax = -1+b*(j1-1)
          ibx = b*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j2 = 1,b
           ia = j4+d*(j3-1+c*(j2+iax))
           ib = j2+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else ! l.eq.1
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          !id(1) = j1
          id(2) = j2
          id(3) = j3
          id(4) = j4
          ibx = a*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j1 = 1,a
           ia = j4+d*(j3-1+c*(j2-1+b*(j1-1)))
           ib = j1+ibx
           s(ib) = u(ia) * f
          enddo
         enddo
        enddo
       enddo
      endif
      return
      end

      subroutine tce_sortacc_4(u,s,a,b,c,d,i,j,k,l,f)
      implicit none
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision s(a*b*c*d)
      double precision u(a*b*c*d)
      double precision f
      integer iax,ibx
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      if (l.eq.4) then
       do j1 = 1,a
        do j2 = 1,b
         do j3 = 1,c
          id(1) = j1
          id(2) = j2
          id(3) = j3
          !id(4) = j4 ! 1234 are only used once and we know l=4, hence ijk must not be 4
          iax = d*(j3-1+c*(j2-1+b*(j1-1)))
          ibx = d*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j4 = 1,d
           ia = j4+iax
           ib = j4+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else if (l.eq.3) then
       do j1 = 1,a
        do j2 = 1,b
         do j4 = 1,d
          id(1) = j1
          id(2) = j2
          !id(3) = j3
          id(4) = j4
          iax = -1+c*(j2-1+b*(j1-1))
          ibx = c*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j3 = 1,c
           ia = j4+d*(j3+iax)
           ib = j3+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else if (l.eq.2) then
       do j1 = 1,a
        do j3 = 1,c
         do j4 = 1,d
          id(1) = j1
          !id(2) = j2
          id(3) = j3
          id(4) = j4
          iax = -1+b*(j1-1)
          ibx = b*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j2 = 1,b
           ia = j4+d*(j3-1+c*(j2+iax))
           ib = j2+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
      else ! l.eq.1
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          !id(1) = j1
          id(2) = j2
          id(3) = j3
          id(4) = j4
          ibx = a*(id(k)-1+jd(k)*(id(j)-1+jd(j)*(id(i)-1)))
          do j1 = 1,a
           ia = j4+d*(j3-1+c*(j2-1+b*(j1-1)))
           ib = j1+ibx
           s(ib) = s(ib) + u(ia) * f
          enddo
         enddo
        enddo
       enddo
      endif
      return
      end
#endif
