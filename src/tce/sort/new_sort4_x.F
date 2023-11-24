      subroutine tce_sort_4_x(u,s,a,b,c,d,i,j,k,l,f)
      implicit none
      integer, intent(in) :: a,b,c,d
      integer, intent(in) :: i,j,k,l
      double precision, intent(out) :: s(a*b*c*d)
      double precision, intent(in) :: u(a*b*c*d)
      double precision, intent(in) :: f
      integer :: version
      version = 1000*i+100*j+10*k+l
      if (version.eq.1234) then
        call transpose_1234_x(u,s,a,b,c,d,f)
      elseif (version.eq.1243) then
        call transpose_1243_x(u,s,a,b,c,d,f)
      elseif (version.eq.1324) then
        call transpose_1324_x(u,s,a,b,c,d,f)
      elseif (version.eq.1342) then
        call transpose_1342_x(u,s,a,b,c,d,f)
      elseif (version.eq.1423) then
        call transpose_1423_x(u,s,a,b,c,d,f)
      elseif (version.eq.1432) then
        call transpose_1432_x(u,s,a,b,c,d,f)
      elseif (version.eq.2134) then
        call transpose_2134_x(u,s,a,b,c,d,f)
      elseif (version.eq.2143) then
        call transpose_2143_x(u,s,a,b,c,d,f)
      elseif (version.eq.2314) then
        call transpose_2314_x(u,s,a,b,c,d,f)
      elseif (version.eq.2341) then
        call transpose_2341_x(u,s,a,b,c,d,f)
      elseif (version.eq.2413) then
        call transpose_2413_x(u,s,a,b,c,d,f)
      elseif (version.eq.2431) then
        call transpose_2431_x(u,s,a,b,c,d,f)
      elseif (version.eq.3124) then
        call transpose_3124_x(u,s,a,b,c,d,f)
      elseif (version.eq.3142) then
        call transpose_3142_x(u,s,a,b,c,d,f)
      elseif (version.eq.3214) then
        call transpose_3214_x(u,s,a,b,c,d,f)
      elseif (version.eq.3241) then
        call transpose_3241_x(u,s,a,b,c,d,f)
      elseif (version.eq.3412) then
        call transpose_3412_x(u,s,a,b,c,d,f)
      elseif (version.eq.3421) then
        call transpose_3421_x(u,s,a,b,c,d,f)
      elseif (version.eq.4123) then
        call transpose_4123_x(u,s,a,b,c,d,f)
      elseif (version.eq.4132) then
        call transpose_4132_x(u,s,a,b,c,d,f)
      elseif (version.eq.4213) then
        call transpose_4213_x(u,s,a,b,c,d,f)
      elseif (version.eq.4231) then
        call transpose_4231_x(u,s,a,b,c,d,f)
      elseif (version.eq.4312) then
        call transpose_4312_x(u,s,a,b,c,d,f)
      elseif (version.eq.4321) then
        call transpose_4321_x(u,s,a,b,c,d,f)
      else
        print*,'something is wrong...'
      endif
      end
        subroutine transpose_1234_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j1 = 1:d1, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1243_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            !s(j3+d3*(j4-1+d4*(j2-1+d2*(j1-1)))) = f*u(j4+d4*(j3-1+d3*(j2-1+d2*(j1-1))))
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j4-1+d4*(j2-1+d2*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1324_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j2-1+d2*(j3-1+d3*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1342_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j4-1+d4*(j3-1+d3*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1423_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j2-1+d2*(j4-1+d4*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_1432_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j4 = 1:d4, j2 = 1:d2)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j3-1+d3*(j4-1+d4*(j1-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2134_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j3-1+d3*(j1-1+d1*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2143_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j4-1+d4*(j1-1+d1*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2314_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j1 = 1:d1, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j1-1+d1*(j3-1+d3*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2341_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j4-1+d4*(j3-1+d3*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2413_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j1 = 1:d1, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j1-1+d1*(j4-1+d4*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_2431_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j3-1+d3*(j4-1+d4*(j2-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3124_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j2-1+d2*(j1-1+d1*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3142_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j4-1+d4*(j1-1+d1*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3214_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j1 = 1:d1, j3 = 1:d3, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j4+d4*(j1-1+d1*(j2-1+d2*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3241_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j4-1+d4*(j2-1+d2*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3412_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j2 = 1:d2, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j1-1+d1*(j4-1+d4*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_3421_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j2-1+d2*(j4-1+d4*(j3-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4123_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j2 = 1:d2, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j2-1+d2*(j1-1+d1*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4132_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j1 = 1:d1, j3 = 1:d3, j4 = 1:d4, j2 = 1:d2)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j3-1+d3*(j1-1+d1*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4213_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j1 = 1:d1, j4 = 1:d4, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j3+d3*(j1-1+d1*(j2-1+d2*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4231_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j4 = 1:d4, j1 = 1:d1, j3 = 1:d3)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j3-1+d3*(j2-1+d2*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4312_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j3 = 1:d3, j1 = 1:d1, j4 = 1:d4, j2 = 1:d2)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j2+d2*(j1-1+d1*(j3-1+d3*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
        subroutine transpose_4321_x(u,s,d1,d2,d3,d4,f)
        implicit none
        integer, intent(in) :: d1,d2,d3,d4
        integer :: old_offset,new_offset
        integer :: j1,j2,j3,j4
        double precision, intent(out) :: s(d1*d2*d3*d4)
        double precision, intent(in) :: u(d1*d2*d3*d4)
        double precision, intent(in) :: f
        do concurrent(j2 = 1:d2, j3 = 1:d3, j1 = 1:d1, j4 = 1:d4)
            old_offset = j4+d4*(j3-1+d3*(j2-1+d2*(j1-1)))
            new_offset = j1+d1*(j2-1+d2*(j3-1+d3*(j4-1)))
            s(new_offset) = u(old_offset) * f
        enddo
        end
