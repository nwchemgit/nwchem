*
* $Id: Grassman.F,v 1.10 2003-10-20 16:46:18 bylaska Exp $
*

* $Log: not supported by cvs2svn $
* Revision 1.9  2003/05/07 21:52:26  bylaska
* ...EJB
*
* Revision 1.8  2003/03/09 00:09:12  bylaska
* Added Grsm_ggm_sym_dot subroutine for more efficient calculation of
* product that result in symmetry matrices...EJB
*

*     ***********************************************************
*     *								*
*     *   		 Grsm (Grassman) Library		*	
*     *	   	      (Pack and D3dB implemenation)		*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*    The Grsm (Grassman) library is to be used for handling
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
*     *    	Grsm_gg_Copy  	 *
*     *                		 *
*     ****************************
*
*     C = A 
*     uses - Pack_c_Copy

      subroutine Grsm_gg_Copy(nfft3d,n,A,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer i

c     do i=1,n
c        call Pack_c_Copy(1,A(1,i),C(1,i))
c     end do
      call dcopy(2*n*nfft3d,A,1,C,1)
 
      return
      end


*     ****************************
*     *                		 *
*     *    	Grsm_hh_Copy  	 *
*     *                		 *
*     ****************************
*
*     C = A 
*     uses - D3dB_r_Copy

      subroutine Grsm_hh_Copy(nfft3d,n,A,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer i

c     do i=1,n
c        call D3dB_r_Copy(1,A(1,i),C(1,i))
c     end do
      call dcopy(2*n*nfft3d,A,1,C,1)
 
      return
      end

*     ****************************
*     *                	    	 *
*     *    	Grsm_gg_daxpy  	     *
*     *                	     	 *
*     ****************************
*
*     C = C + alpha*A
*     uses - Pack_cc_Sum

      subroutine Grsm_gg_daxpy(nfft3d,n,alpha,A,C)
      implicit none
      integer    nfft3d,n
      real*8     alpha
      complex*16 A(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****

      call daxpy(2*n*nfft3d,alpha,A,1,C,1)
c     do i=1,n
c         call Pack_cc_daxpy(1,alpha,A(1,i),C(1,i))
c     end do
 
      return
      end

       
       

*     ****************************
*     *                		 *
*     *    	Grsm_ggg_Sum   	 *
*     *                		 *
*     ****************************
*
*     C = A + B
*     uses - Pack_cc_Sum

      subroutine Grsm_ggg_Sum(nfft3d,n,A,B,C)
      implicit none
      integer    nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer i

      do i=1,n
c         call D3dB_cc_Sum(1,A(1,i),B(1,i),C(1,i))
          call Pack_cc_Sum(1,A(1,i),B(1,i),C(1,i))
      end do
 
      return
      end
       

*     ****************************
*     *                		 *
*     *    	Grsm_ggg_Sub   	 *
*     *                		 *
*     ****************************
*
*     C = A - B
*     uses - Pack_cc_Sub

      subroutine Grsm_ggg_Sub(nfft3d,n,A,B,C)
      implicit none
      integer    nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer i

      do i=1,n
c         call D3dB_cc_Sub(1,A(1,i),B(1,i),C(1,i))
         call Pack_cc_Sub(1,A(1,i),B(1,i),C(1,i))
      end do
 
      return
      end
       
       
*     ****************************
*     *                		 *
*     *    	Grsm_hhh_Sum   	 *
*     *                		 *
*     ****************************
*
*     C = A + B
*     uses - D3dB_rr_Sum

      subroutine Grsm_hhh_Sum(nfft3d,n,A,B,C)
      implicit none
      integer nfft3d,n
      real*8 A(nfft3d,n)
      real*8 B(nfft3d,n)
      real*8 C(nfft3d,n)

*     **** local variables ****
      integer i

      do i=1,n
         call D3dB_rr_Sum(1,A(1,i),B(1,i),C(1,i))
      end do
 
      return
      end
       

*     ****************************
*     *                		 *
*     *    	Grsm_gmg_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*MM
*     uses - Pack_c_SMul, Pack_cc_daxpy

      subroutine Grsm_gmg_Mul(nfft3d,n,A,MM,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      real*8     MM(n,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer j,k

      do k=1,n
c         call D3dB_c_SMul(1,MM(1,k),A(1,1),C(1,k))
c         do j=2,n
c           call D3dB_cc_daxpy(1,MM(j,k),A(1,j),C(1,k))
c         end do
         call Pack_c_SMul(1,MM(1,k),A(1,1),C(1,k))
         do j=2,n
           call Pack_cc_daxpy(1,MM(j,k),A(1,j),C(1,k))
         end do

      end do

      return
      end

*     ****************************
*     *                		 *
*     *    	Grsm_gmg_daxpy   	 *
*     *                		 *
*     ****************************
*
*     C = C + A*MM
*     uses - Pack_c_SMul, Pack_cc_daxpy

      subroutine Grsm_gmg_daxpy(nfft3d,n,A,MM,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      real*8     MM(n,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer j,k

      do k=1,n
         do j=1,n
           call Pack_cc_daxpy(1,MM(j,k),A(1,j),C(1,k))
         end do
      end do

      return
      end



*     ****************************
*     *                		 *
*     *    	Grsm_gqg_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*Q
*     uses - Pack_c_SMul, Pack_cc_zaxpy

      subroutine Grsm_gqg_Mul(nfft3d,n,A,Q,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 Q(n,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer j,k


      do k=1,n
c         call D3dB_c_ZMul(1,Q(1,k),A(1,1),C(1,k))
c         do j=2,n
c           call D3dB_cc_zaxpy(1,Q(j,k),A(1,j),C(1,k))
c         end do
         call Pack_c_ZMul(1,Q(1,k),A(1,1),C(1,k))
         do j=2,n
           call Pack_cc_zaxpy(1,Q(j,k),A(1,j),C(1,k))
         end do
      end do

      return
      end

*     ****************************
*     *                		 *
*     *    	Grsm_hmh_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*M
*     uses - D3dB_r_SMul, D3dB_rr_daxpy

      subroutine Grsm_hmh_Mul(nfft3d,n,A,MM,C)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      real*8     MM(n,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer j,k


      do k=1,n
         call D3dB_r_SMul(1,MM(1,k),A(1,1),C(1,k))
         do j=2,n
           call D3dB_rr_daxpy(1,MM(j,k),A(1,j),C(1,k))
         end do
      end do

      return
      end


*     ****************************
*     *                		 *
*     *    	Grsm_gg_dScale 	 *
*     *                		 *
*     ****************************
*
*     C = scale*A
*     uses - Pack_c_SMul

      subroutine Grsm_gg_dScale(nfft3d,n,scale,A,C)
      implicit none
      integer nfft3d,n
      real*8     scale
      complex*16 A(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer k

      do k=1,n
c        call D3dB_c_SMul(1,scale,A(1,k),C(1,k))
        call Pack_c_SMul(1,scale,A(1,k),C(1,k))

      end do

      return
      end



*     ****************************
*     *                		 *
*     *    	Grsm_hh_dScale 	 *
*     *                		 *
*     ****************************
*
*     C = scale*A
*     uses - D3dB_r_SMul

      subroutine Grsm_hh_dScale(nfft3d,n,scale,A,C)
      implicit none
      integer nfft3d,n
      real*8 scale
      complex*16 A(nfft3d,n)
      complex*16 C(nfft3d,n)

*     **** local variables ****
      integer k


      do k=1,n
        call D3dB_r_SMul(1,scale,A(1,k),C(1,k))
      end do

      return
      end



*     ****************************
*     *                		 *
*     *    	Grsm_gh_fftb 	 *
*     *                		 *
*     ****************************
*
*     A = FFT(-1)[A]
*     uses - D3dB_cr_fft3b,Pack_c_unpack

      subroutine Grsm_gh_fftb(nfft3d,n,A)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k


      do k=1,n
        call Pack_c_unpack(1,A(1,k))
        call D3dB_cr_fft3b(1,A(1,k))
      end do

      return
      end



*     ****************************
*     *                		 *
*     *    	Grsm_hg_fftf 	 *
*     *                		 *
*     ****************************
*
*     A = FFT[A]
*     uses - D3dB_rc_fft3f,Pack_c_pack

      subroutine Grsm_hg_fftf(nfft3d,n,A)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k


      do k=1,n
        call D3dB_rc_fft3f(1,A(1,k))
        call Pack_c_pack(1,A(1,k))
      end do

      return
      end


*     ****************************
*     *                		 *
*     *    	Grsm_g_mask 	 *
*     *                		 *
*     ****************************
*
*     uses - mask_C

      subroutine Grsm_g_mask(nfft3d,n,nb,A)
      implicit none
      integer nfft3d,n,nb
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k


      do k=1,n
        call mask_C(nb,A(1,k))
      end do

      return
      end

*     ****************************
*     *                		 *
*     *    	Grsm_ggm_dot 	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Grsm_ggm_dot(nfft3d,n,A,B,MM)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      real*8     MM(n,n)

*     **** local variables ****
      integer j,k,np
      real*8 sum

      call Parallel_np(np)
      call dcopy(n*n,0.0d0,0,MM,1)
      do k=1,n
      do j=1,n
       call Pack_cc_idot(1,A(1,j),B(1,k),sum) 
       MM(j,k) = sum
      end do
      end do
      if (np.gt.1) call D3dB_Vector_Sumall(n*n,MM)

      return
      end


*     ****************************
*     *                          *
*     *     Grsm_ggm_sym_dot     *
*     *                          *
*     ****************************
*
      subroutine Grsm_ggm_sym_dot(nfft3d,n,A,B,MM)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      real*8     MM(n,n)

*     **** local variables ****
      integer j,k,np
      real*8 sum

      call Parallel_np(np)
      call dcopy(n*n,0.0d0,0,MM,1)
      do k=1,n
         call Pack_cc_idot(1,A(1,k),B(1,k),sum)
         MM(k,k) = sum
        do j=k+1,n
         call Pack_cc_idot(1,A(1,j),B(1,k),sum)
         MM(j,k) = sum
         MM(k,j) = sum
        end do
      end do
      if (np.gt.1) call D3dB_Vector_Sumall(n*n,MM)

      return
      end


*     ****************************
*     *                		 *
*     *    	Grsm_ggm_dot2 	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Grsm_ggm_dot2(nfft3d,m,n,A,B,MM)
      implicit none
      integer nfft3d,m,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      real*8     MM(m,m)

*     **** local variables ****
      integer j,k,np
      real*8 sum

      call Parallel_np(np)
      call dcopy(m*m,0.0d0,0,MM,1)
      do k=1,n
      do j=1,n
c       call D3dB_cc_idot(1,A(1,j),B(1,k),sum)
       call Pack_cc_idot(1,A(1,j),B(1,k),sum) 
       MM(j,k) = sum
      end do
      end do
      if (np.gt.1) call D3dB_Vector_Sumall(m*m,MM)

      return
      end


*     ****************************
*     *                		 *
*     *    	Grsm_gg_trace 	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Grsm_gg_trace(nfft3d,n,A,B,MM)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      real*8     MM

*     **** local variables ****
      integer j,np
      real*8 sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


*     **** external functions ****
      real*8   ddot
      external ddot


      call Parallel_np(np)
      MM = 0.0d0
      do j=1,n
c        call Pack_cc_idot(1,A(1,j),B(1,j),sum)
c        MM = MM + sum
        MM = MM 
     >     + ddot(2*nida(1),A(1,j),1,B(1,j),1)
     >     + 2*ddot(2*nidb(1),A(nida(1)+1,j),1,B(nida(1)+1,j),1)
      end do
      if (np.gt.1) call D3dB_Sumall(MM)

      return
      end




*     ****************************
*     *                		 *
*     *    Grsm_h_Zero_Ends	 *
*     *                		 *
*     ****************************
*
*     uses - 

      subroutine Grsm_h_Zero_Ends(nfft3d,n,A)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k

      do k=1,n
        call D3dB_r_Zero_Ends(1,A(1,k))
      end do

      return
      end


*     ****************************
*     *                		 *
*     *    Grsm_g_project	 *
*     *                		 *
*     ****************************
*

      subroutine Grsm_g_project(nfft3d,n,Z,Y,YtZ)
      implicit none
      integer	 nfft3d,n
      complex*16 Z(nfft3d,n)
      complex*16 Y(nfft3d,n)
      real*8     YtZ(n,n)

*     ***** local variables *****
      integer i,j,np
      real*8  sum

      call Parallel_np(np)
      do j=1,n
         do i=1,n
c            call D3dB_cc_idot(1,Y(1,i),Z(1,j),sum)
            call Pack_cc_idot(1,Y(1,i),Z(1,j),sum)
            YtZ(i,j) = sum
         end do
      end do
      if (np.gt.1) call D3dB_Vector_Sumall(n*n,YtZ)


      return
      end


*     **************************************
*     *                                    *
*     *              Grsm_g_QR             *
*     *                                    *
*     **************************************

*     This routine calculates the QR decomposition
* 
*       A-->QR
*
*     Entry - npack,ne,Q
*     Exit  - Q,R

      subroutine Grsm_g_QR(npack,ne,Q,R)
      implicit none
      integer npack,ne
      complex*16 Q(npack,ne)
      real*8     R(ne,ne)

*     **** local variables ****
      integer j,k
      real*8  w

      call dcopy(ne*ne,0.0d0,0,R,1)
      do k=1,ne
         call Pack_cc_dot(1,Q(1,k),Q(1,k),w)
         R(k,k) = dsqrt(w)
         w = 1.0d0/R(k,k)
         call Pack_c_SMul(1,w,Q(1,k),Q(1,k))

         do j=k+1,ne
            call Pack_cc_dot(1,Q(1,k),Q(1,j),w)
            R(k,j) = w
            w = -w
            call Pack_cc_daxpy(1,w,Q(1,k),Q(1,j))
         end do
      end do

      return
      end 









