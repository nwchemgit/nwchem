*     ******************************
*     *                            *
*     *    Calculate_cpsi_spin2    *
*     *                            *
*     ******************************
c  Calculate the spin contamination for the band case
c  This is defined as (see Szabo and Ostlund auf seite 107)
c     DeltaS= N(down) - sum (i,j) dconjg(S(i,j))*(Si,j)
C  Note that S(i,j) is _complex_!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calculate_cpsi_spin2(nbrillq,ispin,ne,nfft3d,
     >                                psi_tag,Sab)
      implicit none
      integer nbrillq,ispin,ne(2)
      integer nfft3d
      integer psi_tag
      real*8 Sab

#include "util.fh"
#include "bafdecls.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical oprint
      integer i,j,nb
      integer psii_shift,psij_shift
      real*8  Sexact,SS,w,ssum
      complex*16 sij

*     **** external functions ****
      logical  control_print
      integer  cpsi_data_get_ptr
      real*8   brillioun_weight
      external control_print
      external cpsi_data_get_ptr
      external brillioun_weight


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)


*     ****************************
*     ***** calculate <S^2> ******
*     ****************************

      if (ispin .eq. 2) then
         Sexact = ((ne(1)-ne(2))/2.0d0)
     >          * (((ne(1)-ne(2))/2.0d0)+1.0d0)
         Sab = 0.0d0
         do nb=1,nbrillq
            w = brillioun_weight(nb)
	    ssum=0.0d0
            do i=1,ne(1)
            do j=1,ne(2)
               psii_shift = cpsi_data_get_ptr(psi_tag,nb,i)
               psij_shift = cpsi_data_get_ptr(psi_tag,nb,j+ne(1))
 
*              *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
               call Cram_cc_zdot(nb,dbl_mb(psii_shift),
     >                              dbl_mb(psij_shift),sij)
	       SS=dble(dconjg(sij)*sij)
               ssum = ssum + SS
            end do
            end do
            Sab=Sab+w*(dble(ne(2))-ssum)
         end do
         call K1dB_SumAll(Sab)
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
c $Id$
