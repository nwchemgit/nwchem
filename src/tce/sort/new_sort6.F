!
! written by Jeff Hammond, Intel Labs, April 2014
! see https://github.com/jeffhammond/nwchem-tce-triples-kernels/blob/master/src/new_sort6.f for history
!
! There is currently no difference between the OpenMP and non-OpenMP versions
! except for the directives, but preprocessing at fine granularity leads to ugly
! code because Fortran does not have anything like C99/C++11 Pragma.
!
#ifdef USE_OPENMP
      subroutine tce_sort_6(u,s,a,b,c,d,e,f,i,j,k,l,m,n,fac)
      implicit none
      integer a,b,c,d,e,f
      integer i,j,k,l,m,n
      double precision s(a*b*c*d*e*f)
      double precision u(a*b*c*d*e*f)
      double precision fac
      integer id(6),jd(6)
      integer ia,ib
      integer j1,j2,j3,j4,j5,j6
      integer jdm,jdl,jdk,jdj
      integer iax,ibx
      ia = 0
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      if (n.eq.6) then
      jdm=jd(m)
      jdl=jd(l)
      jdk=jd(k)
      jdj=jd(j)
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j5 = 1,e
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           !id(6) = j6 ! 123456 are only used once and we know n=6, hence ijklm must not be 6
           iax = f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
           ibx = f*(id(m)-1+jdm*(id(l)-1+jdl*(id(k)-1+jdk
     &                *(id(j)-1+jdj*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j6 = 1,f
            ia = j6+iax
            ib = j6+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.5) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           !id(5) = j5
           id(6) = j6
           iax = -1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1))))
           ibx = e*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j5 = 1,e
            ia = j6+f*(j5+iax)
            ib = j5+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.4) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           !id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+d*(j3-1+c*(j2-1+b*(j1-1)))
           ibx = d*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j4 = 1,d
            ia = j6+f*(j5-1+e*(j4+iax))
            ib = j4+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.3) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           !id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+c*(j2-1+b*(j1-1))
           ibx = c*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j3 = 1,c
            ia = j6+f*(j5-1+e*(j4-1+d*(j3+iax)))
            ib = j3+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.2) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           !id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+b*(j1-1)
           ibx = b*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j2 = 1,b
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2+iax))))
            ib = j2+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else ! if (n.eq.1) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j2 = 1,b
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           !id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           ibx = a*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j1 = 1,a
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
            ib = j1+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      return
      end

      subroutine tce_sortacc_6(u,s,a,b,c,d,e,f,i,j,k,l,m,n,fac)
      implicit none
      integer a,b,c,d,e,f
      integer i,j,k,l,m,n
      double precision s(a*b*c*d*e*f)
      double precision u(a*b*c*d*e*f)
      double precision fac
      integer id(6),jd(6)
      integer ia,ib
      integer j1,j2,j3,j4,j5,j6
      integer jdm,jdl,jdk,jdj
      integer iax,ibx
      ia = 0
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      if (n.eq.6) then
      jdm=jd(m)
      jdl=jd(l)
      jdk=jd(k)
      jdj=jd(j)
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j5 = 1,e
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           !id(6) = j6 ! 123456 are only used once and we know n=6, hence ijklm must not be 6
           iax = f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
           ibx = f*(id(m)-1+jdm*(id(l)-1+jdl*(id(k)-1+jdk
     &                *(id(j)-1+jdj*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j6 = 1,f
            ia = j6+iax
            ib = j6+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.5) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           !id(5) = j5
           id(6) = j6
           iax = -1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1))))
           ibx = e*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j5 = 1,e
            ia = j6+f*(j5+iax)
            ib = j5+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.4) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           !id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+d*(j3-1+c*(j2-1+b*(j1-1)))
           ibx = d*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j4 = 1,d
            ia = j6+f*(j5-1+e*(j4+iax))
            ib = j4+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.3) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j2 = 1,b
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           !id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+c*(j2-1+b*(j1-1))
           ibx = c*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j3 = 1,c
            ia = j6+f*(j5-1+e*(j4-1+d*(j3+iax)))
            ib = j3+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.2) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j1 = 1,a
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           !id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+b*(j1-1)
           ibx = b*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j2 = 1,b
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2+iax))))
            ib = j2+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else ! if (n.eq.1) then
!$omp  parallel do collapse(5)
!$omp& private(j1,j2,j3,j4,j5,j6,iax,ibx,ia,ib,id) shared(s,u)
!$omp& firstprivate(a,b,c,d,e,f,i,j,k,l,fac,jd,jdm,jdl,jdk,jdj)
      do j2 = 1,b
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           !id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           ibx = a*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
!dec$ vector always nontemporal
!$omp simd
           do j1 = 1,a
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
            ib = j1+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      return
      end
#else
      subroutine tce_sort_6(u,s,a,b,c,d,e,f,i,j,k,l,m,n,fac)
      implicit none
      integer a,b,c,d,e,f
      integer i,j,k,l,m,n
      double precision s(a*b*c*d*e*f)
      double precision u(a*b*c*d*e*f)
      double precision fac
      integer id(6),jd(6)
      integer ia,ib
      integer j1,j2,j3,j4,j5,j6
      integer jdm,jdl,jdk,jdj
      integer iax,ibx
      ia = 0
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      if (n.eq.6) then
      jdm=jd(m)
      jdl=jd(l)
      jdk=jd(k)
      jdj=jd(j)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j5 = 1,e
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           !id(6) = j6 ! 123456 are only used once and we know n=6, hence ijklm must not be 6
           iax = f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
           ibx = f*(id(m)-1+jdm*(id(l)-1+jdl*(id(k)-1+jdk
     &                *(id(j)-1+jdj*(id(i)-1)))))
           do j6 = 1,f
            ia = j6+iax
            ib = j6+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.5) then
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           !id(5) = j5
           id(6) = j6
           iax = -1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1))))
           ibx = e*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j5 = 1,e
            ia = j6+f*(j5+iax)
            ib = j5+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.4) then
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           !id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+d*(j3-1+c*(j2-1+b*(j1-1)))
           ibx = d*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j4 = 1,d
            ia = j6+f*(j5-1+e*(j4+iax))
            ib = j4+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.3) then
      do j1 = 1,a
       do j2 = 1,b
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           !id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+c*(j2-1+b*(j1-1))
           ibx = c*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j3 = 1,c
            ia = j6+f*(j5-1+e*(j4-1+d*(j3+iax)))
            ib = j3+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.2) then
      do j1 = 1,a
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           !id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+b*(j1-1)
           ibx = b*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j2 = 1,b
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2+iax))))
            ib = j2+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else ! if (n.eq.1) then
      do j2 = 1,b
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           !id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           ibx = a*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j1 = 1,a
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
            ib = j1+ibx
            s(ib) = u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      return
      end

      subroutine tce_sortacc_6(u,s,a,b,c,d,e,f,i,j,k,l,m,n,fac)
      implicit none
      integer a,b,c,d,e,f
      integer i,j,k,l,m,n
      double precision s(a*b*c*d*e*f)
      double precision u(a*b*c*d*e*f)
      double precision fac
      integer id(6),jd(6)
      integer ia,ib
      integer j1,j2,j3,j4,j5,j6
      integer jdm,jdl,jdk,jdj
      integer iax,ibx
      ia = 0
      jd(1) = a
      jd(2) = b
      jd(3) = c
      jd(4) = d
      jd(5) = e
      jd(6) = f
      if (n.eq.6) then
      jdm=jd(m)
      jdl=jd(l)
      jdk=jd(k)
      jdj=jd(j)
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j5 = 1,e
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           !id(6) = j6 ! 123456 are only used once and we know n=6, hence ijklm must not be 6
           iax = f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
           ibx = f*(id(m)-1+jdm*(id(l)-1+jdl*(id(k)-1+jdk
     &                *(id(j)-1+jdj*(id(i)-1)))))
           do j6 = 1,f
            ia = j6+iax
            ib = j6+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.5) then
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j4 = 1,d
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           !id(5) = j5
           id(6) = j6
           iax = -1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1))))
           ibx = e*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j5 = 1,e
            ia = j6+f*(j5+iax)
            ib = j5+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.4) then
      do j1 = 1,a
       do j2 = 1,b
        do j3 = 1,c
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           id(3) = j3
           !id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+d*(j3-1+c*(j2-1+b*(j1-1)))
           ibx = d*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j4 = 1,d
            ia = j6+f*(j5-1+e*(j4+iax))
            ib = j4+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.3) then
      do j1 = 1,a
       do j2 = 1,b
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           id(2) = j2
           !id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+c*(j2-1+b*(j1-1))
           ibx = c*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j3 = 1,c
            ia = j6+f*(j5-1+e*(j4-1+d*(j3+iax)))
            ib = j3+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else if (n.eq.2) then
      do j1 = 1,a
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           id(1) = j1
           !id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           iax = -1+b*(j1-1)
           ibx = b*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j2 = 1,b
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2+iax))))
            ib = j2+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      else ! if (n.eq.1) then
      do j2 = 1,b
       do j3 = 1,c
        do j4 = 1,d
         do j5 = 1,e
          do j6 = 1,f
           !id(1) = j1
           id(2) = j2
           id(3) = j3
           id(4) = j4
           id(5) = j5
           id(6) = j6
           ibx = a*(id(m)-1+jd(m)*(id(l)-1+jd(l)*(id(k)-1
     &           +jd(k)*(id(j)-1+jd(j)*(id(i)-1)))))
           do j1 = 1,a
            ia = j6+f*(j5-1+e*(j4-1+d*(j3-1+c*(j2-1+b*(j1-1)))))
            ib = j1+ibx
            s(ib) = s(ib) + u(ia) * fac
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      return
      end
#endif
