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

#include "util.fh"
#include "mafdecls.fh"
#include "frac_occ.fh"

*     **** local variables ****
      logical oprint
      integer i,j,ms
      real*8  Sexact,SS,en(2)
       
      integer taskid,MASTER
      parameter (MASTER=0)

*     **** external functions ****
      logical  control_print
      external control_print


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     ****************************
*     ***** calculate <S^2> ******
*     ****************************

      if (ispin .eq. 2) then

         en(1)     = ne(1)
         en(ispin) = ne(ispin)
         if (fractional) then
            do ms=1,ispin
            en(ms) =0.0
            do i=1,ne(ms)
              j = fweight(1)+(i-1)+(ms-1)*ne(1)
              en(ms) = en(ms) + dbl_mb(j)
            end do
            end do
         end if

         Sexact = ((en(1)-en(2))/2.0d0)
     >          * (((en(1)-en(2))/2.0d0)+1.0d0)
         Sab = 0.0d0
         do i=1,ne(1)
         do j=1,ne(2)
*           *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
            call Pack_cc_dot(1,psi(1,i),psi(1,j+ne(1)),SS)
            Sab = Sab + SS*SS
         end do
         end do
         Sab = Sexact + en(2) - Sab
      else
         Sexact = 0.0d0
         Sab    = 0.0d0
      end if

      if (oprint) then
         WRITE(*,'(//A)') '=== Spin Contamination ==='
         write(*,*)
         write(*,*) "<Sexact^2> = ", Sexact
         write(*,*) "<S^2>      = ", Sab
         write(*,*)
      end if

      return
      end
