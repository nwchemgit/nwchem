*     ******************************
*     *                            *
*     *    Calculate_psi_spin2     *
*     *                            *
*     ******************************

      subroutine Calculate_psi_spin2(ispin,ne,nfft3d,psi,Sab)
      implicit none
      integer ispin,ne(2)
      integer nfft3d
      complex*16 psi(nfft3d,*)
      real*8 Sab

*     **** local variables ****
      integer i,j
      real*8  Sexact,SS
       
      integer taskid,MASTER
      parameter (MASTER=0)


      call Parallel_taskid(taskid)

*     ****************************
*     ***** calculate <S^2> ******
*     ****************************

      if (ispin .eq. 2) then
         Sexact = (dble(ne(1)-ne(2))/2.0d0)
     >          * ((dble(ne(1)-ne(2))/2.0d0)+1.0d0)
         Sab = 0.0d0
         do i=1,ne(1)
         do j=1,ne(2)
*           *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
            call Pack_cc_dot(1,psi(1,i),psi(1,j+ne(1)),SS)
            Sab = Sab + SS*SS
         end do
         end do
         Sab = Sexact + ne(2) - Sab
      else
         Sexact = 0.0d0
         Sab    = 0.0d0
      end if

      if (taskid.eq.MASTER) then
         WRITE(*,'(//A)') '=== Spin Contamination ==='
         write(*,*)
         write(*,*) "<Sexact^2> = ", Sexact
         write(*,*) "<S^2>      = ", Sab
         write(*,*)
      end if

      return
      end
