*     ***********************************************************
*     *								*
*     *   		 Stfl (Stiefel) Library		*	
*     *	   	      (Pack and D3dB implemenation)		*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*     The Stfl (Stiefel) library is to be used for handling
* four kinds of data structures.  The first data structure,
* denoted by "g", is a set of n D3dB "c" complex vectors.
* The second data structure, denoted by "h", is a set of n
* D3dB "r" real vectors.  The third data structure, denoted by
* "m", is a (n X n) real matrix.  The fourth data structure denoted
* by "q" is a (n X n) complex matrix.
*
*  Uses: D3dB(1,..) data structure, blas3
*  Assumes: that D3dB(1,...) data structure is initialized
*
*  Future: This data structure will replace D3dB
* 

*     ****************************
*     *                		 *
*     *    	Stfl_gg_Copy  	 *
*     *                		 *
*     ****************************
*
*     C = A 
*     uses - Pack_c_Copy

      subroutine Stfl_gg_Copy(nfft3d,n,nbrill,A,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,n
         call Cram_c_Copy(nb,A(1,i,nb),C(1,i,nb))
      end do
      end do
 
      return
      end


*     ****************************
*     *                		 *
*     *    	Stfl_hh_Copy  	 *
*     *                		 *
*     ****************************
*
*     C = A 
*     uses - C3dB_c_Copy

      subroutine Stfl_hh_Copy(nfft3d,n,nbrill,A,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,n
         call C3dB_c_Copy(1,A(1,i,nb),C(1,i,nb))
      end do
      end do
      return
      end

*     ****************************
*     *                		 *
*     *    	Stfl_ggg_Sum   	 *
*     *                		 *
*     ****************************
*
*     C = A + B
*     uses - Cram_cc_Sum

      subroutine Stfl_ggg_Sum(nfft3d,n,nbrill,A,B,C)
      implicit none
      integer    nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 B(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,n
          call Cram_cc_Sum(nb,A(1,i,nb),B(1,i,nb),C(1,i,nb))
      end do
      end do
 
      return
      end
       
*     ****************************
*     *                		 *
*     *    	Stfl_ggg_Sub   	 *
*     *                		 *
*     ****************************
*
*     C = A - B
*     uses - Pack_cc_Sub

      subroutine Stfl_ggg_Sub(nfft3d,n,nbrill,A,B,C)
      implicit none
      integer    nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 B(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,n
         call Cram_cc_Sub(nb,A(1,i,nb),B(1,i,nb),C(1,i,nb))
      end do
      end do
      return
      end
       
       
*     ****************************
*     *                		 *
*     *    	Stfl_hhh_Sum   	 *
*     *                		 *
*     ****************************
*
*     C = A + B
*     uses - C3dB_cc_Sum

      subroutine Stfl_hhh_Sum(nfft3d,n,nbrill,A,B,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 B(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,n
         call C3dB_cc_Sum(1,A(1,i,nb),B(1,i,nb),C(1,i,nb))
      end do
      end do
 
      return
      end
       
*     ****************************
*     *                		 *
*     *    	Stfl_gmg_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*MM
*     uses - Pack_c_SMul, Pack_cc_daxpy

      subroutine Stfl_gmg_Mul(nfft3d,n,nbrill,A,MM,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      real*8     MM(n,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer j,k,nb

      do nb=1,nbrill
      do k=1,n
         call Cram_c_SMul(nb,MM(1,k,nb),A(1,1,nb),C(1,k,nb))
         do j=2,n
           call Cram_cc_daxpy(nb,MM(j,k,nb),A(1,j,nb),C(1,k,nb))
         end do

      end do
      end do
      return
      end


*     ****************************
*     *                		 *
*     *    	Stfl_gqg_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*Q
*     uses - Cram_c_ZMul, Cram_cc_zaxpy

      subroutine Stfl_gqg_Mul(nfft3d,n,nbrill,A,Q,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 Q(n,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer j,k,nb

      do nb=1,nbrill
      do k=1,n
         call Cram_c_ZMul(nb,Q(1,k,nb),A(1,1,nb),C(1,k,nb))
         do j=2,n
           call Cram_cc_zaxpy(nb,Q(j,k,nb),A(1,j,nb),C(1,k,nb))
         end do
      end do
      end do
      return
      end

*     ****************************
*     *                		 *
*     *    	Stfl_hmh_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*M
*     uses - C3dB_c_SMul, C3dB_cc_daxpy

      subroutine Stfl_hmh_Mul(nfft3d,n,nbrill,A,MM,C)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      real*8     MM(n,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer j,k,nb

      do nb=1,nbrill
      do k=1,n
         call C3dB_c_SMul(1,MM(1,k,nb),A(1,1,nb),C(1,k,nb))
         do j=2,n
           call C3dB_cc_daxpy(1,MM(j,k,nb),A(1,j,nb),C(1,k,nb))
         end do
      end do
      end do
      return
      end


*     ****************************
*     *                		 *
*     *    	Stfl_gg_dScale 	 *
*     *                		 *
*     ****************************
*
*     C = scale*A
*     uses - Cram_c_SMul

      subroutine Stfl_gg_dScale(nfft3d,n,nbrill,scale,A,C)
      implicit none
      integer nfft3d,n,nbrill
      real*8     scale
      complex*16 A(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,n
        call Cram_c_SMul(nb,scale,A(1,k,nb),C(1,k,nb))
      end do
      end do
      return
      end


*     ****************************
*     *                		 *
*     *    	Stfl_hh_dScale 	 *
*     *                		 *
*     ****************************
*
*     C = scale*A
*     uses - C3dB_c_SMul

      subroutine Stfl_hh_dScale(nfft3d,n,nbrill,scale,A,C)
      implicit none
      integer nfft3d,n,nbrill
      real*8 scale
      complex*16 A(nfft3d,n,nbrill)
      complex*16 C(nfft3d,n,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,n
        call C3dB_c_SMul(1,scale,A(1,k,nb),C(1,k,nb))
      end do
      end do
      return
      end



*     ****************************
*     *                		 *
*     *    	Stfl_gh_fftb 	 *
*     *                		 *
*     ****************************
*
*     A = FFT(-1)[A]
*     uses - C3dB_cr_fft3b,Cram_c_unpack

      subroutine Stfl_gh_fftb(nfft3d,n,nbrill,A)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,n
        call Cram_c_unpack(nb,A(1,k,nb))
        call C3dB_cr_fft3b(1,A(1,k,nb))
      end do
      end do
      return
      end



*     ****************************
*     *                		 *
*     *    	Stfl_hg_fftf 	 *
*     *                		 *
*     ****************************
*
*     A = FFT[A]
*     uses - C3dB_rc_fft3f,Cram_c_pack

      subroutine Stfl_hg_fftf(nfft3d,n,nbrill,A)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,n
        call C3dB_rc_fft3f(1,A(1,k,nb))
        call Cram_c_pack(nb,A(1,k,nb))
      end do
      end do

      return
      end

*     ****************************
*     *                		 *
*     *    	Stfl_ggq_dot 	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Stfl_ggq_dot(nfft3d,n,nbrill,A,B,MM)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 B(nfft3d,n,nbrill)
      complex*16 MM(n,n,nbrill)

*     **** local variables ****
      integer j,k,np,nb
      complex*16 sum

      call Parallel_np(np)
      call dcopy(2*n*n*nbrill,0.0d0,0,MM,1)
      do nb=1,nbrill
      do k=1,n
      do j=1,n
       call Cram_cc_izdot(nb,A(1,j,nb),B(1,k,nb),sum) 
       MM(j,k,nb) = sum
      end do
      end do
      end do
      if (np.gt.1) call C3dB_Vector_Sumall(2*n*n*nbrill,MM)
      return
      end

*     ****************************
*     *                		 *
*     *    	Stfl_gg_trace 	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Stfl_gg_trace(nfft3d,n,nbrill,A,B,MM)
      implicit none
      integer nfft3d,n,nbrill
      complex*16 A(nfft3d,n,nbrill)
      complex*16 B(nfft3d,n,nbrill)
      real*8     MM(nbrill)

*     **** local variables ****
      integer j,nb,np
      real*8 sum

      call Parallel_np(np)

      do nb=1,nbrill
         MM(nb) = 0.0d0
         do j=1,n
           call Cram_cc_idot(1,A(1,j,nbrill),B(1,j,nbrill),sum)
           MM(nb) = MM(nb) + sum
         end do
      end do
      if (np.gt.1) call C3dB_Vector_Sumall(nbrill,MM)
      return
      end


*     ****************************
*     *                		 *
*     *    Stfl_g_project	 *
*     *                		 *
*     ****************************
*

      subroutine Stfl_g_project(nfft3d,n,nbrill,Z,Y,YtZ)
      implicit none
      integer	 nfft3d,n,nbrill
      complex*16 Z(nfft3d,n,nbrill)
      complex*16 Y(nfft3d,n,nbrill)
      complex*16 YtZ(n,n,nbrill)

*     ***** local variables *****
      integer i,j,nb,np
      complex*16  sum

      call Parallel_np(np)
      do nb=1,nbrill
      do j=1,n
      do i=1,n
         call Cram_cc_izdot(nb,Y(1,i,nbrill),Z(1,j,nbrill),sum)
         YtZ(i,j,nb) = sum
      end do
      end do
      end do
      if (np.gt.1) call C3dB_Vector_Sumall(2*nbrill*n*n,YtZ)
      return
      end











