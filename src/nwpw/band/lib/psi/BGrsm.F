* $Id$
*


*     ***********************************************************
*     *								*
*     *   		 BGrsm (Bundled Grassmann) Library	*	
*     *	   	      (Pack and D3dB implemenation)		*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*     The BGrsm (Bundled Grassmann) library is to be used for handling
* three kinds of data structures.  The first data structure,
* denoted by "f", is a set of n Cram "c" complex vectors.
* The second data structure, denoted by "h", is a set of n
* C3dB "c" complex vectors.   The third data structure denoted by
* "r" is a C3dB "r" real vector.
*
*  Uses: C3dB(1,..) data structure, blas3
*  Assumes: that C3dB(1,...) data structure is initialized
*
* 

*     ****************************
*     *                		 *
*     *        BGrsm_ff_Copy  	 *
*     *                		 *
*     ****************************
*     C = A 
      subroutine BGrsm_ff_Copy(A,C)
      implicit none
      complex*16 A(*)
      complex*16 C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*npack1
      call dcopy(nsize,A,1,C,1)
      return
      end



*     ****************************
*     *                          *
*     *     BGrsm_ff_Copy_tag    *
*     *                          *
*     ****************************
*     C = A
      subroutine BGrsm_ff_Copy_tag(A_tag,C_tag)
      implicit none
      integer A_tag,C_tag
      call cpsi_data_copyall(A_tag,C_tag)
      return
      end




*     ****************************
*     *                		 *
*     *    	BGrsm_f_Zero  	 *
*     *                		 *
*     ****************************
*     A = 0 
      subroutine BGrsm_f_Zero(A)
      implicit none
      complex*16 A(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*npack1
      call dcopy(nsize,0.0d0,0,A,1)
      return
      end 



*     ****************************
*     *                          *
*     *      BGrsm_f_Zero_tag    *
*     *                          *
*     ****************************
*     A = 0
      subroutine BGrsm_f_Zero_tag(A_tag)
      implicit none
      integer A_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,nsize,a_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call cpsi_data_update(A_tag)
      nsize = 2*(neq(1)+neq(2))*npack1
      do nb=1,nbrillq
         a_shift = cpsi_data_get_chnk(A_tag,nb)
         call dcopy(nsize,0.0d0,0,dbl_mb(a_shift),1)
      end do
      call cpsi_data_noupdate(A_tag)
      return
      end

*     ****************************
*     *                		 *
*     *    	BGrsm_hh_Copy  	 *
*     *                		 *
*     ****************************
*
*     C = A 

      subroutine BGrsm_hh_Copy(A,C)
      implicit none
      complex*16 A(*)
      complex*16 C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*nfft3d
      call dcopy(nsize,A,1,C,1)
      return
      end

*     ****************************
*     *                		 *
*     *      BGrsm_fff_Sum   	 *
*     *                		 *
*     ****************************
*
*     C = A + B

      subroutine BGrsm_fff_Sum(A,B,C)
      implicit none
      real*8     A(*)
      real*8     B(*)
      real*8     C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer i,nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*npack1
      do i=1,nsize
         C(i) = A(i) + B(i)
      end do

      return
      end


*     ****************************
*     *                          *
*     *      BGrsm_fff_Sum_tag   *
*     *                          *
*     ****************************
*     C = A + B
      subroutine BGrsm_fff_sum_sub(nsize,A,B,C)
      implicit none
      integer nsize
      real*8 A(*),B(*),C(*)
      integer i
      do i=1,nsize
         C(i) = A(i) + B(i)
      end do
      return
      end
      subroutine BGrsm_fff_Sum_tag(A_tag,B_tag,C_tag)
      implicit none
      integer    A_tag
      integer    B_tag
      integer    C_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,nsize
      integer A_shift,B_shift,C_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call cpsi_data_update(C_tag)
      nsize = 2*(neq(1)+neq(2))*npack1
      do nb=1,nbrillq
         A_shift  = cpsi_data_get_chnk(A_tag,nb)
         B_shift  = cpsi_data_get_chnk(B_tag,nb)
         C_shift  = cpsi_data_get_chnk(C_tag,nb)
         call BGrsm_fff_sum_sub(nsize,dbl_mb(A_shift),
     >                                dbl_mb(B_shift),
     >                                dbl_mb(C_shift))
      end do
      call cpsi_data_noupdate(C_tag)
      return
      end


       
*     ****************************
*     *                		 *
*     *      BGrsm_fff_Sub   	 *
*     *                		 *
*     ****************************
*     C = A - B
      subroutine BGrsm_fff_Sub(A,B,C)
      implicit none
      real*8     A(*)
      real*8     B(*)
      real*8     C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer i,nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*npack1
      do i=1,nsize
         C(i) = A(i) - B(i)
      end do
      return
      end

*     ****************************
*     *                          *
*     *      BGrsm_fff_Sub_tag   *
*     *                          *
*     ****************************
*
*     C = A - B
      subroutine BGrsm_fff_sub_sub(nsize,A,B,C)
      implicit none
      integer nsize
      real*8 A(*),B(*),C(*)
      integer i
      do i=1,nsize
         C(i) = A(i) - B(i)
      end do
      return
      end
      subroutine BGrsm_fff_Sub_tag(A_tag,B_tag,C_tag)
      implicit none
      integer    A_tag
      integer    B_tag
      integer    C_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,nsize
      integer A_shift,B_shift,C_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      nsize = 2*(neq(1)+neq(2))*npack1
      call cpsi_data_update(C_tag)
      do nb=1,nbrillq
         A_shift  = cpsi_data_get_chnk(A_tag,nb)
         B_shift  = cpsi_data_get_chnk(B_tag,nb)
         C_shift  = cpsi_data_get_chnk(C_tag,nb)
         call BGrsm_fff_sub_sub(nsize,dbl_mb(A_shift),
     >                                dbl_mb(B_shift),
     >                                dbl_mb(C_shift))
      end do
      call cpsi_data_noupdate(C_tag)
      return
      end



*     ****************************
*     *                		 *
*     *    	BGrsm_hhh_Sum  	 *
*     *                		 *
*     ****************************
*
*     C = A + B

      subroutine BGrsm_hhh_Sum(A,B,C)
      implicit none
      real*8     A(*)
      real*8     B(*)
      real*8     C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer i,nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*nfft3d
      do i=1,nsize
         C(i) = A(i) - B(i)
      end do
 
      return
      end      

       
*     ****************************
*     *                		 *
*     *       BGrsm_hr_SumSqr 	 *
*     *                		 *
*     ****************************
*
*     C(1:nfft3d,ms) = Sum(ispin,ne,nbrill) A(1:nfft3d,n,nbrill)**2*brillioun_weight()
*     uses - C3dB_cr_aSqrpy

      subroutine  BGrsm_hr_SumSqr(A,C)
      implicit none

#include "cpsi_common.fh"

      complex*16 A(nfft3d,(neq(1)+neq(2)),nbrillq)
      real*8     C(nfft3d,2)


*     **** local variables ****
      integer nbq,nb,ms,n,n1(2),n2(2)
      real*8  weight

*     **** local functions ****
      real*8   brillioun_weight
      external brillioun_weight

      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)

      call C3dB_r_Zero(1,C)
      if (ispin.gt.1) call C3dB_r_Zero(1,C(1,ispin))

      do nbq=1,nbrillq
        nb = nbq
        weight = brillioun_weight(nb)
        do ms=1,ispin
        do n=n1(ms),n2(ms)
           call C3dB_cr_aSqrpy(1,weight,A(1,n,nb),C(1,ms))
        end do
        end do
      end do
      return
      end                   
       

*     ****************************
*     *                          *
*     *       BGrsm_hr_aSqrpy    *
*     *                          *
*     ****************************
*
*     C(1:nfft3d,ms) = Sum(ispin,ne) A(1:nfft3d,n)**2*weight(n)
*     uses - C3dB_cr_aSqrpy

      subroutine  BGrsm_hr_aSqrpy(nfft3d,ne,weight,A,C)
      implicit none
      integer nfft3d,ne
      complex*16 A(nfft3d,ne)
      real*8     weight(*)
      real*8     C(nfft3d)

*     **** local variables ****
      integer n
      real*8  scal2

*     ***** external functions *****
      real*8   lattice_omega
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()

      do n=1,ne
        call C3dB_cr_aSqrpy(1,scal2*weight(n),A(1,n),C)
      end do
      return
      end



*     ****************************
*     *                		 *
*     *     BGrsm_ff_dScale 	 *
*     *                		 *
*     ****************************
*     C = scale*A
      subroutine BGrsm_ff_dScale(scale,A,C)
      implicit none
      real*8     scale
      real*8     A(*)
      real*8     C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer i,nsize
      nsize = 2*nbrillq*(neq(1)+neq(2))*npack1
      do i=1,nsize
         C(i) = scale*A(i) 
      end do
      return
      end


*     ****************************
*     *                          *
*     *      BGrsm_f_dScale_tag  *
*     *                          *
*     ****************************
*     A = alpha*A
      subroutine BGrsm_f_dScale_tag(scale,A_tag)
      implicit none
      real*8 scale
      integer    A_tag

#include "bafdecls.fh"
#include "cpsi_common.fh"

*     **** local variables ****
      integer nb,nsize
      integer A_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      nsize = 2*(neq(1)+neq(2))*npack1
      call cpsi_data_update(A_tag)
      do nb=1,nbrillq
         A_shift  = cpsi_data_get_chnk(A_tag,nb)
         call dscal(nsize,scale,dbl_mb(A_shift),1)
      end do
      call cpsi_data_noupdate(A_tag)
      return
      end



*     ****************************
*     *                		 *
*     *      BGrsm_hh_dScale 	 *
*     *                		 *
*     ****************************
*
*     C = scale*A

      subroutine BGrsm_hh_dScale(scale,A,C)
      implicit none
      real*8 scale
      real*8     A(*)
      real*8     C(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer i,nsize

      nsize = 2*nbrillq*(neq(1)+neq(2))*nfft3d
      do i=1,nsize
         C(i) = scale*A(i) 
      end do
      return
      end

*     ****************************
*     *                		 *
*     *    	BGrsm_fh_fftb 	 *
*     *                		 *
*     ****************************
*
*     B = FFT(-1)[A]
*     uses - C3dB_cr_fft3b,Cram_c_unpack

      subroutine BGrsm_fh_fftb(A,B)
      implicit none
      complex*16 A(*)
      complex*16 B(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer k,nbq,nb,shiftA,shiftB

      shiftA = 1
      shiftB = 1
      do nbq=1,nbrillq
         nb = nbq
         do k=1,(neq(1)+neq(2))
            call Cram_c_Copy(nb,A(shiftA),B(shiftB))
            call Cram_c_unpack(nb,B(shiftB))
            !call C3dB_cr_fft3b(1,B(shiftB))
            call C3dB_cr_pfft3b(1,nb,B(shiftB))
            shiftA = shiftA + npack1
            shiftB = shiftB + nfft3d
         end do
      end do
      return
      end



*     ****************************
*     *                		 *
*     *      BGrsm_hf_fftf 	 *
*     *                		 *
*     ****************************
*
*     A = FFT[A]
*     uses - C3dB_rc_fft3f,Cram_c_pack

      subroutine BGrsm_hg_fftf(B,A)
      implicit none
      complex*16 B(*)
      complex*16 A(*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer k,nbq,nb,shiftA,shiftB

      shiftA = 1
      shiftB = 1
      do nbq=1,nbrillq
         nb = nbq
         do k=1,(neq(1)+neq(2))
           !call C3dB_rc_fft3f(1,B(shiftB))
           call C3dB_rc_pfft3f(1,nb,B(shiftB))
           call Cram_c_pack(nb,B(shiftB))
           call Cram_c_Copy(nb,B(shiftB),A(shiftA))
           shiftA = shiftA + npack1
           shiftB = shiftB + nfft3d
         end do
      end do
      return
      end

*     ****************************
*     *                		 *
*     *     BGrsm_f_Spins        *
*     *                		 *
*     ****************************

*   Returns the  Spin Amplitudes for eigenvalue printout.
*  !!! this routine needs to be checked!!!

c      logical function BGrsm_f_Spins(nbb,Y,S)
      subroutine BGrsm_f_Spins(nbb,Y,S)
      implicit none
      integer nbb
      complex*16 Y(*)
      real*8     S(*)

#include "cpsi_common.fh"
#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      integer k,shft,nb,nb1,nb2,indx,indx2,sindx
      real*8 sx,sy,sz,wght,ssum,s12,s1,s2
      complex*16 z1

*     **** external ****
      real*8   brillioun_weight
      external brillioun_weight

*     **** simple error checking - make sure we are doing a spin-orbit calculation ***
      if (.not.spin_orbit) then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER)
     >      write(luout,*) 
     >      'Warning: BGrsm_f_Spins called when spin-orbit not set'
         return
      end if

      if (nbb.eq.0) then
        nb1=1
        nb2=nbrillq
      else
        nb1=nbb
        nb2=nbb 
      end if
cccccccccccccccccccccccccccccccccccccc
c Should this be done separately for
c each BZ point? Now we are averaging
c over k-points...
cccccccccccccccccccccccccccccccccccccc
      sindx=1
      do k=1,ne(1)
         S(sindx)=0.0d0
         S(sindx+1)=0.0d0
         S(sindx+2)=0.0d0
         sindx=sindx+3
      end do
      shft=ne(1)*npack1
      do nb=nb1,nb2
        wght=brillioun_weight(nb)
        if (nb1.eq.nb2) wght=1.0d0
        indx=((nb-nb1)*shft*2)+1
        sindx=1
        do k=1,ne(1)
           indx2=indx+shft
           call Cram_cc_dot(nb,Y(indx),Y(indx),s1)
           call Cram_cc_dot(nb,Y(indx2),Y(indx2),s2)
           call Cram_cc_zdot(nb,Y(indx),Y(indx2),z1)
           indx=indx+npack1
           sz=s1-s2
           sx=dble(z1)
           sy=dimag(z1)
           S(sindx+2)=S(sindx+2)+sz*wght
           S(sindx+1)=S(sindx+1)+sy*wght
           S(sindx)=S(sindx)+sx*wght
           sindx=sindx+3		
        end do
      end do
      return
      end
