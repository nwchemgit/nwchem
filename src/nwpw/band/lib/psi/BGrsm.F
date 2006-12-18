* $Id: BGrsm.F,v 1.3 2006-12-18 05:36:00 d3p708 Exp $
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
*
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
*     *                		 *
*     *    	BGrsm_f_Zero  	 *
*     *                		 *
*     ****************************
*
*     C = A 

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
*     *                		 *
*     *      BGrsm_fff_Sub   	 *
*     *                		 *
*     ****************************
*
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
*
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
            call C3dB_cr_fft3b(1,B(shiftB))
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
           call C3dB_rc_fft3f(1,B(shiftB))
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
*     *      Returns the  	 *
*     *  Spin Amplitudes         *
*     ****************************
*
*     
*     
      logical function BGrsm_f_Spins(nb,Y,S)
      implicit none
      integer nb
      complex*16 Y(*)
      real*8 S(3,*)

#include "cpsi_common.fh"

*     **** local variables ****
      integer k,ms,shft,nbq,nb1,nb2,indx
      real*8 s1,s2,s3,s4,wght

*     **** external ****
      real*8 zdotc,brillioun_weight
      external zdotc,brillioun_weight

      if (.not.spin_orbit) then
         BGrsm_f_Spins=.false.
         return
      else
         BGrsm_f_Spins=.true.
      end if
      if (nb.eq.0) then
        nb1=1
        nb2=nbrillq
      else
        nb1=nb
        nb2=nb 
      end if
      do k=1,ne(1)
        S(1,k)=0.0d0
        S(2,k)=0.0d0
        S(3,k)=0.0d0
      end do
      shft=ne(1)*npack1
      do nbq=nb1,nb2
        wght=brillioun_weight(nbq)
        if (nb1.eq.nb2) wght=1.0d0
        indx=(nbq-1)*shft*2+1
        do k=1,ne(1)
           s1=zdotc(npack1,Y(indx),1,Y(indx),1)
           s2=zdotc(npack1,Y(indx+shft),1,Y(indx+shft),1)
           s3=zdotc(npack1,Y(indx),1,Y(indx+shft),1)
           s4=zdotc(npack1,Y(indx+shft),1,Y(indx),1)
           indx=indx+npack1
           S(1,k)=S(1,k)+(s1-s2)*wght
           S(2,k)=S(2,k)+s3*wght
           S(3,k)=S(3,k)+s4*wght		
        end do
      end do
 900  continue
      return
      end

