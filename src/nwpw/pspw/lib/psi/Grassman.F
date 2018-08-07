*
* $Id$
*

* $Log: not supported by cvs2svn $
* Revision 1.19  2009/02/07 03:50:56  bylaska
* Bassi Vectorization Fix...EJB
*
* Revision 1.18  2007/03/29 00:51:45  bylaska
* queue FFT added to HFX and more rigid msgtyping added to queueFFT
* ....EJB
*
* Revision 1.17  2006/10/07 00:30:07  bylaska
* Initial implementation of 2d processor grid parallelization in pspw.  Currently works with:
*
* task pspw steepest_descent
* task pspw energy            (only minimizer 1, minimizer 2?, other minimizers not yet implemented)
*
* Currently only works with USE_MPIF option, straight tcgmsg only partially implemented.  Car-Parrinello, HFX, SIC, and various analysis codes are also not yet ported.
*
*
* The number of processors along the orbital dimension, np_orbital, is entered as follows, e.g.:
*
* nwpw
*    np_orbital  2
* end
*
* The number of processors along the grid dimension, np_grid, is currently defined using np_orbital as
*
* np_grid = np/np_orbital
*
* where np is the total number of processors.
*
* ...EJB
*
* Revision 1.16  2005/07/12 21:26:33  bylaska
* Intitial implementation of pipelined parallel FFT implemented.
*
* ....EJB
*
* Revision 1.15  2005/06/24 19:24:18  bylaska
* various FFT updates, and symmetry point paths added to band_structure...EJB
*
* Revision 1.14  2004/02/22 18:56:38  bylaska
* Input and io changes for simulation_cell
* Grsm_ggm2_sym_dot routine used in paw_overlap_matrix_gen
* ...EBJ
*
* Revision 1.13  2003/10/31 22:16:38  bylaska
* bug fixes ... EJB
*
* Revision 1.12  2003/10/26 22:24:04  bylaska
* performance updates??...EJB
*
* Revision 1.11  2003/10/23 01:37:55  bylaska
* Added Grsm_ggm2_dot subroutine ....EJB
*
* Revision 1.10  2003/10/20 16:46:18  bylaska
* ..EJB
*
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
      !call dcopy(2*n*nfft3d,A,1,C,1)
      call Parallel_shared_vector_copy(.true.,2*n*nfft3d,A,C)
 
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
*     *                          *
*     *         Grsm_ggg_Sum2    *
*     *                          *
*     ****************************
*
*     B = A + B
*     uses - Pack_cc_Sum2

      subroutine Grsm_ggg_Sum2(nfft3d,n,A,B)
      implicit none
      integer    nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)

*     **** local variables ****
      integer i

      do i=1,n
c         call D3dB_cc_Sum(1,A(1,i),B(1,i),C(1,i))
          call Pack_cc_Sum2(1,A(1,i),B(1,i))
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
*     *                          *
*     *         Grsm_gg_dScale1  *
*     *                          *
*     ****************************
*
*     C = scale*A
*     uses - Pack_c_SMul1

      subroutine Grsm_gg_dScale1(nfft3d,n,scale,A)
      implicit none
      integer nfft3d,n
      real*8     scale
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k

      do k=1,n
        call Pack_c_SMul1(1,scale,A(1,k))
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
      logical done
      integer k,indx1,indx2


*     **** external functions ****
      logical  D3dB_cr_pfft3_queue_filled
      external D3dB_cr_pfft3_queue_filled
c      integer tid,Parallel_threadid
c      external    Parallel_threadid

c      tid = Parallel_threadid()

c      do k=1,n
c        call Pack_c_unpack(1,A(1,k))
c        !call D3dB_cr_fft3b(1,A(1,k))
c        !call D3dB_cr_pfft3b(1,1,A(1,k))
c      end do
c      call D3dB_cr_mpfft3b(1,1,n,A)

      indx1=1
      indx2=1
      done = .false.
      do while (.not.done)
        if (indx1.le.n) then
           call D3dB_cr_pfft3b_queuein(1,A(1,indx1))
           indx1 = indx1 + 1
        end if

        if ((D3dB_cr_pfft3_queue_filled()).or.(indx1.gt.n)) then
           call D3dB_cr_pfft3b_queueout(1,A(1,indx2))
           indx2 = indx2 + 1
        end if
        done = ((indx1.gt.n).and.(indx2.gt.n))
      end do

      return
      end



*     ****************************
*     *                		 *
*     *    	Grsm_hg_fftf0 	 *
*     *                		 *
*     ****************************
*
*     A = FFT[A]
*     uses - D3dB_rc_fft3f,Pack_c_mpack

      subroutine Grsm_hg_fftf0(nfft3d,n,A)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)

*     **** local variables ****
      integer k


      !call D3dB_nrc_fft3f(1,n,A)
      do k=1,n
        call D3dB_rc_fft3f(1,A(1,k))
      end do
      call Pack_c_mpack(1,n,A)

      return
      end

*     ****************************
*     *                          *
*     *         Grsm_hg_fftf     *
*     *                          *
*     ****************************
*
*     A = FFT[A]
*     uses - D3dB_rc_pfft3_queue_filled,D3dB_rc_pfft3f_queuein,D3dB_rc_pfft3f_queueout

      subroutine Grsm_hg_fftf(nfft3d,n,A)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)

*     **** local variables ****
      logical done
      integer k,indx1,indx2

*     **** external functions ****
      logical  D3dB_rc_pfft3_queue_filled
      external D3dB_rc_pfft3_queue_filled

      indx1=1
      indx2=1
      done = .false.
      do while (.not.done)
        if (indx1.le.n) then
           call D3dB_rc_pfft3f_queuein(1,A(1,indx1))
           indx1 = indx1 + 1
        end if

        if ((D3dB_rc_pfft3_queue_filled()).or.(indx1.gt.n)) then
           call D3dB_rc_pfft3f_queueout(1,A(1,indx2))
           indx2 = indx2 + 1
        end if
        done = ((indx1.gt.n).and.(indx2.gt.n))
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
      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


      call Parallel_np(np)
      !call dcopy(n*n,0.0d0,0,MM,1)
      call Parallel_shared_vector_zero(.true.,n*n,MM)
      do k=1,n
      do j=1,n
       call Pack_cc_idot(1,A(1,j),B(1,k),e1) 
!$OMP MASTER
       MM(j,k) = e1
!$OMP END MASTER
      end do
      end do
!$OMP BARRIER
      if (np.gt.1) call D3dB_Vector_SumAll(n*n,MM)

      return
      end



*     ****************************
*     *                          *
*     *        Grsm_ggm2_dot     *
*     *                          *
*     ****************************
*
*     uses -

      subroutine Grsm_ggm2_dot(nfft3d,n,A,B,MM,m)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)

      integer m
      real*8     MM(m,m)

*     **** local variables ****
      integer j,k,np
      real*8 sum

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


      call Parallel_np(np)
      !call dcopy(m*m,0.0d0,0,MM,1)
      call Parallel_shared_vector_zero(.true.,m*m,MM)
      do k=1,n
      do j=1,n
       call Pack_cc_idot(1,A(1,j),B(1,k),e1)
!$OMP MASTER
       MM(j,k) = e1
!$OMP END MASTER
      end do
      end do
!$OMP BARRIER
      if (np.gt.1) call D3dB_Vector_Sumall(m*m,MM)

      return
      end

*     ********************************
*     *                              *
*     *        Grsm_ggm2_sym_dot     *
*     *                              *
*     ********************************
*
*     uses -

      subroutine Grsm_ggm2_sym_dot(nfft3d,n,A,B,MM,m)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)

      integer m
      real*8     MM(m,m)

*     **** local variables ****
      integer j,k,np
      real*8 sum

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2

      call Parallel_np(np)
      !call dcopy(m*m,0.0d0,0,MM,1)
      call Parallel_shared_vector_zero(.true.,m*m,MM)
      do k=1,n
         call Pack_cc_idot(1,A(1,k),B(1,k),e1)
!$OMP MASTER
         MM(k,k) = e1
!$OMP END MASTER
        do j=k+1,n
         call Pack_cc_idot(1,A(1,j),B(1,k),e1)
!$OMP MASTER
         MM(j,k) = e1
         MM(k,j) = e1
!$OMP END MASTER
        end do
      end do
!$OMP BARRIER
      if (np.gt.1) call D3dB_Vector_Sumall(m*m,MM)


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

      real*8 e1,e2
      common /eenergy_tmp_common/ e1,e2


      call Parallel_np(np)
      call dcopy(n*n,0.0d0,0,MM,1)
      do k=1,n
         call Pack_cc_idot(1,A(1,k),B(1,k),e1)
!$OMP MASTER
         MM(k,k) = e1
!$OMP END MASTER
        do j=k+1,n
         call Pack_cc_idot(1,A(1,j),B(1,k),e1)
!$OMP MASTER
         MM(j,k) = e1
         MM(k,j) = e1
!$OMP END MASTER
        end do
      end do
!$OMP BARRIER
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
      integer j,np,nida,nidb
      !real*8 sum

      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb


c*     **** common block for pack ****
c      integer nida(0:3),nidb(0:3)
c      integer pack(2,0:3)
c      common / pack_blk / pack,nida,nidb


*     **** external functions ****
      real*8   ddot
      external ddot

      call Parallel2d_np_i(np)
      call Pack_nida(1,nida)
      call Pack_nidb(1,nidb)
!$OMP MASTER
      sa = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:sa)
      do j=1,n
c        call Pack_cc_idot(1,A(1,j),B(1,j),sum)
c        MM = MM + sum
        sa = sa 
     >     + ddot(2*nida,A(1,j),1,B(1,j),1)
     >     + 2*ddot(2*nidb,A(nida+1,j),1,B(nida+1,j),1)
      end do
!$OMP END DO

!$OMP MASTER
      MM = sa
!$OMP END MASTER
!$OMP BARRIER
c      if (np.gt.1) call D3dB_SumAll(MM)

      if (np.gt.1) call D3dB_SumAll(MM)

      return
      end


*     ****************************
*     *                          *
*     *        Grsm_gg_itrace    *
*     *                          *
*     ****************************
* 
*     uses - 

      subroutine Grsm_gg_itrace(nfft3d,n,A,B,MM)
      implicit none
      integer nfft3d,n
      complex*16 A(nfft3d,n)
      complex*16 B(nfft3d,n)
      real*8     MM

*     **** local variables ****
      integer j,np,nida,nidb
      real*8 sum

      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb


c*     **** common block for pack ****
c      integer nida(0:3),nidb(0:3)
c      integer pack(2,0:3)
c      common / pack_blk / pack,nida,nidb


*     **** external functions ****
      real*8   ddot
      external ddot

      call Pack_nida(1,nida)
      call Pack_nidb(1,nidb)

c      call Parallel2d_np_i(np)
!$OMP MASTER
      sa = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:sa)
      do j=1,n
c        call Pack_cc_idot(1,A(1,j),B(1,j),sum)
c        MM = MM + sum
        sa = sa
     >     + ddot(2*nida,A(1,j),1,B(1,j),1)
     >     + 2*ddot(2*nidb,A(nida+1,j),1,B(nida+1,j),1)
      end do
!$OMP END DO

!$OMP MASTER
      MM = sa
!$OMP END MASTER
!$OMP BARRIER
c      if (np.gt.1) call D3dB_SumAll(MM)

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

      call Parallel2d_np_i(np)
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









