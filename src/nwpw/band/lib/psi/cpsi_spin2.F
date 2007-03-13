*     ******************************
*     *                            *
*     *    Calculate_cpsi_spin2    *
*     *                            *
*     ******************************

      subroutine Calculate_cpsi_spin2(nbrill,ispin,ne,nfft3d,psi,Sab)
      implicit none
      integer nbrill,ispin,ne(2)
      integer nfft3d
      complex*16 psi(nfft3d,*)
      real*8 Sab

#include "util.fh"
#include "mafdecls.fh"

*     **** local variables ****
      logical oprint
      integer i,j,nb
      integer ii,jj
      real*8  Sexact,SS,en(2),w,ssum
       
      integer taskid,MASTER
      parameter (MASTER=0)

*     **** external functions ****
      logical  control_print
      real*8   brillioun_weight
      external control_print
      external brillioun_weight


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     ****************************
*     ***** calculate <S^2> ******
*     ****************************

      if (ispin .eq. 2) then

         en(1)     = ne(1)
         en(ispin) = ne(ispin)

         w = 0.0d0
         do nb=1,nbrill
            w = w + brillioun_weight(nb)
         end do
         en(1)     = w*ne(1)
         en(ispin) = w*ne(ispin)

         Sexact = ((en(1)-en(2))/2.0d0)
     >          * (((en(1)-en(2))/2.0d0)+1.0d0)
         Sab = 0.0d0
         do nb=1,nbrill
            w = brillioun_weight(nb)
	    ssum=0.0d0
            do i=1,ne(1)
            do j=1,ne(2)
               ii = i      +(nb-1)*(ne(1)+ne(2))
               jj = j+ne(1)+(nb-1)*(ne(1)+ne(2))
 
*              *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
               call Cram_cc_dot(nb,psi(1,ii),psi(1,jj),SS)
               ssum = ssum + SS*SS
c               write(*,*)i,j,SS,(SS*SS)
            end do
            end do
c            write(*,*)"nb= ",nb," ssum =",ssum
c            write(*,*)"delS = ",(ne(2)-ssum)
c            call flush(6)
            Sab=Sab+w*(ne(2)-ssum)
         end do
c         write(*,*)"Sab = ",Sab 
         Sab = Sexact + Sab
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
