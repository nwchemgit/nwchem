
      subroutine tce_sort_4(unsorted,sorted,a,b,c,d,
     1                      i,j,k,l,factor)
      implicit none
#include "util.fh"
#include "tce_main.fh"
      integer version
      integer a,b,c,d
      integer i,j,k,l
      integer id(4),jd(4),ia,ib,j1,j2,j3,j4
      double precision sorted(a*b*c*d)
      double precision unsorted(a*b*c*d)
      double precision factor
c
      version = 1000*i+100*j+10*k+l
c
c transpose_1234.log:        The best loop order is:        2134
c transpose_1243.log:        The best loop order is:        1243
c transpose_1324.log:        The best loop order is:        1234
c transpose_1342.log:        The best loop order is:        1342
c transpose_1423.log:        The best loop order is:        1423
c transpose_1432.log:        The best loop order is:        1342
c
      if (version.eq.1234) then
        call transpose_1234_loop_2134(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.1243) then
        call transpose_1243_loop_1243(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.1324) then
        call transpose_1324_loop_1234(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.1342) then
        call transpose_1342_loop_1342(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.1423) then
        call transpose_1423_loop_1423(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.1432) then
        call transpose_1432_loop_1342(unsorted,sorted,a,b,c,d,factor)
c
c transpose_2134.log:        The best loop order is:        2134
c transpose_2143.log:        The best loop order is:        2143
c transpose_2314.log:        The best loop order is:        2134
c transpose_2341.log:        The best loop order is:        2341
c transpose_2413.log:        The best loop order is:        2413
c transpose_2431.log:        The best loop order is:        2341
c
      elseif (version.eq.2134) then
        call transpose_2134_loop_2134(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.2143) then
        call transpose_2143_loop_2143(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.2314) then
        call transpose_2314_loop_2134(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.2341) then
        call transpose_2341_loop_2341(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.2413) then
        call transpose_2413_loop_2413(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.2431) then
        call transpose_2431_loop_2341(unsorted,sorted,a,b,c,d,factor)
c
c transpose_3124.log:        The best loop order is:        1234
c transpose_3142.log:        The best loop order is:        1342
c transpose_3214.log:        The best loop order is:        2134
c transpose_3241.log:        The best loop order is:        2341
c transpose_3412.log:        The best loop order is:        1342
c transpose_3421.log:        The best loop order is:        2341
c
      elseif (version.eq.3124) then
        call transpose_3124_loop_1234(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.3142) then
        call transpose_3142_loop_1342(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.3214) then
        call transpose_3214_loop_2134(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.3241) then
        call transpose_3241_loop_2341(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.3412) then
        call transpose_3412_loop_1342(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.3421) then
        call transpose_3421_loop_2341(unsorted,sorted,a,b,c,d,factor)
c
c transpose_4123.log:        The best loop order is:        1423
c transpose_4132.log:        The best loop order is:        1342
c transpose_4213.log:        The best loop order is:        2413
c transpose_4231.log:        The best loop order is:        2341
c transpose_4312.log:        The best loop order is:        1342
c transpose_4321.log:        The best loop order is:        2341
c
      elseif (version.eq.4123) then
        call transpose_4123_loop_1423(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.4132) then
        call transpose_4132_loop_1342(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.4213) then
        call transpose_4213_loop_2413(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.4231) then
        call transpose_4231_loop_2341(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.4312) then
        call transpose_4312_loop_1342(unsorted,sorted,a,b,c,d,factor)
      elseif (version.eq.4321) then
        call transpose_4321_loop_2341(unsorted,sorted,a,b,c,d,factor)
c
      else
        print*,'something is wrong...'
      endif

      return
      end
        subroutine transpose_4321_loop_3241(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer xdim1,xdim4,rdim1,rdim4
        integer old offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
        rdim1=mod(dim1,4)
        rdim4=mod(dim4,4)
        xdim1=dim1-rdim1
        xdim4=dim4-rdim4
!DEC$ ivdep
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ vector always
!DEC$ loop count min(24), max(40), avg(32)
        do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
         do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ vector always
          do j4 = 1,xdim4,4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ vector always
           do j1 = 1,xdim1,4
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(1+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(1+j1-1))))
            sorted(2+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(2+j1-1))))
            sorted(3+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(3+j1-1))))

            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(1+j4-1)))) = factor*
     &    unsorted(1+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(1+j1+dim1*(j2-1+dim2*(j3-1+dim3*(1+j4-1)))) = factor*
     &    unsorted(1+j4+dim4*(j3-1+dim3*(j2-1+dim2*(1+j1-1))))
            sorted(2+j1+dim1*(j2-1+dim2*(j3-1+dim3*(1+j4-1)))) = factor*
     &    unsorted(1+j4+dim4*(j3-1+dim3*(j2-1+dim2*(2+j1-1))))
            sorted(3+j1+dim1*(j2-1+dim2*(j3-1+dim3*(1+j4-1)))) = factor*
     &    unsorted(1+j4+dim4*(j3-1+dim3*(j2-1+dim2*(3+j1-1))))

            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(2+j4-1)))) = factor*
     &    unsorted(2+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(1+j1+dim1*(j2-1+dim2*(j3-1+dim3*(2+j4-1)))) = factor*
     &    unsorted(2+j4+dim4*(j3-1+dim3*(j2-1+dim2*(1+j1-1))))
            sorted(2+j1+dim1*(j2-1+dim2*(j3-1+dim3*(2+j4-1)))) = factor*
     &    unsorted(2+j4+dim4*(j3-1+dim3*(j2-1+dim2*(2+j1-1))))
            sorted(3+j1+dim1*(j2-1+dim2*(j3-1+dim3*(2+j4-1)))) = factor*
     &    unsorted(2+j4+dim4*(j3-1+dim3*(j2-1+dim2*(3+j1-1))))

            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(3+j4-1)))) = factor*
     &    unsorted(3+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(1+j1+dim1*(j2-1+dim2*(j3-1+dim3*(3+j4-1)))) = factor*
     &    unsorted(3+j4+dim4*(j3-1+dim3*(j2-1+dim2*(1+j1-1))))
            sorted(2+j1+dim1*(j2-1+dim2*(j3-1+dim3*(3+j4-1)))) = factor*
     &    unsorted(3+j4+dim4*(j3-1+dim3*(j2-1+dim2*(2+j1-1))))
            sorted(3+j1+dim1*(j2-1+dim2*(j3-1+dim3*(3+j4-1)))) = factor*
     &    unsorted(3+j4+dim4*(j3-1+dim3*(j2-1+dim2*(3+j1-1))))
           enddo
!DEC$ loop count min(0), max(4), avg(2)
!DEC$ vector always
           do j1 = xdim1+1,dim1,1
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(1+j4-1)))) = factor*
     &    unsorted(1+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(2+j4-1)))) = factor*
     &    unsorted(2+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(3+j4-1)))) = factor*
     &    unsorted(3+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
           enddo
          enddo
!DEC$ loop count min(0), max(4), avg(2)
!DEC$ vector always
          do j4 = xdim4+1,dim4,1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ vector always
           do j1 = 1,xdim1,4
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
            sorted(1+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(1+j1-1))))
            sorted(2+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(2+j1-1))))
            sorted(3+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(3+j1-1))))
           enddo
          enddo
!DEC$ loop count min(0), max(4), avg(2)
!DEC$ vector always
          do j4 = xdim4+1,dim4,1
!DEC$ loop count min(0), max(4), avg(2)
!DEC$ vector always
           do j1 = xdim1+1,dim1,1
            sorted(0+j1+dim1*(j2-1+dim2*(j3-1+dim3*(0+j4-1)))) = factor*
     &    unsorted(0+j4+dim4*(j3-1+dim3*(j2-1+dim2*(0+j1-1))))
           enddo
          enddo
         enddo
        enddo
        return
        end
c
c
c  ADD AUTO-GENERATED CODE HERE
c
c
        subroutine transpose_1234_loop_2134(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_1243_loop_1243(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j4-1+dim4*(j2-1+dim2*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_1324_loop_1234(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j2-1+dim2*(j3-1+dim3*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_1342_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j4-1+dim4*(j3-1+dim3*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_1423_loop_1423(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j2-1+dim2*(j4-1+dim4*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_1432_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j3-1+dim3*(j4-1+dim4*(j1-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2134_loop_2134(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j3-1+dim3*(j1-1+dim1*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2143_loop_2143(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j4-1+dim4*(j1-1+dim1*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2314_loop_2134(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j1-1+dim1*(j3-1+dim3*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2341_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j4-1+dim4*(j3-1+dim3*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2413_loop_2413(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j1-1+dim1*(j4-1+dim4*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_2431_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j3-1+dim3*(j4-1+dim4*(j2-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3124_loop_1234(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j2-1+dim2*(j1-1+dim1*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3142_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j4-1+dim4*(j1-1+dim1*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3214_loop_2134(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j4 = 1,dim4
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j4+dim4*(j1-1+dim1*(j2-1+dim2*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3241_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j4-1+dim4*(j2-1+dim2*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3412_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j1-1+dim1*(j4-1+dim4*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_3421_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j2-1+dim2*(j4-1+dim4*(j3-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4123_loop_1423(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j2-1+dim2*(j1-1+dim1*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4132_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j3-1+dim3*(j1-1+dim1*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4213_loop_2413(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j3 = 1,dim3
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j3+dim3*(j1-1+dim1*(j2-1+dim2*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4231_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j3-1+dim3*(j2-1+dim2*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4312_loop_1342(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j1 = 1,dim1
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j2 = 1,dim2
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j2+dim2*(j1-1+dim1*(j3-1+dim3*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
        subroutine transpose_4321_loop_2341(unsorted,sorted,
     &                           dim1,dim2,dim3,dim4,factor)
        implicit none
        integer dim1,dim2,dim3,dim4
        integer old_offset,new_offset
        integer j1,j2,j3,j4
        double precision sorted(dim1*dim2*dim3*dim4)
        double precision unsorted(dim1*dim2*dim3*dim4)
        double precision factor
!DEC$ prefetch sorted
!DEC$ prefetch unsorted
!DEC$ ivdep
!DEC$ loop count min(24), max(40), avg(32)
        do j2 = 1,dim2
!DEC$ loop count min(24), max(40), avg(32)
         do j3 = 1,dim3
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
          do j4 = 1,dim4
!DEC$ loop count min(24), max(40), avg(32)
!DEC$ unroll(8)
!DEC$ vector always
           do j1 = 1,dim1
            old_offset = j4+dim4*(j3-1+dim3*(j2-1+dim2*(j1-1)))
            new_offset = j1+dim1*(j2-1+dim2*(j3-1+dim3*(j4-1)))
            sorted(new_offset) = unsorted(old_offset) * factor
           enddo
          enddo
         enddo
        enddo
        return
        end
c $Id$
