*
* $Id: Stiefel.F,v 1.6 2002-05-03 19:44:01 bylaska Exp $
*


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
*     uses - Cram_c_Copy

      subroutine Stfl_gg_Copy(ispin,ne,nfft3d,neall,nbrill,A,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
         call Cram_c_Copy(nb,A(1,i,nb),C(1,i,nb))
      end do
      end do
      return
      end

*     ****************************
*     *                		 *
*     *    	Stfl_g_Zero  	 *
*     *                		 *
*     ****************************
*
*     C = A 
*     uses - Cram_c_Zero

      subroutine Stfl_g_Zero(ispin,ne,nfft3d,neall,nbrill,A)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
         call Cram_c_Zero(nb,A(1,i,nb))
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

      subroutine Stfl_hh_Copy(ispin,ne,nfft3d,neall,nbrill,A,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
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

      subroutine Stfl_ggg_Sum(ispin,ne,nfft3d,neall,nbrill,A,B,C)
      implicit none
      integer    ispin,ne(2)
      integer    nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
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
*     uses - Cram_cc_Sub

      subroutine Stfl_ggg_Sub(ispin,ne,nfft3d,neall,nbrill,A,B,C)
      implicit none
      integer    ispin,ne(2)
      integer    nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
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

      subroutine Stfl_hhh_Sum(ispin,ne,nfft3d,neall,nbrill,A,B,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer i,nb

      do nb=1,nbrill
      do i=1,(ne(1)+ne(2))
         call C3dB_cc_Sum(1,A(1,i,nb),B(1,i,nb),C(1,i,nb))
      end do
      end do
 
      return
      end      

       
*     ****************************
*     *                		 *
*     *       Stfl_hr_SumSqr 	 *
*     *                		 *
*     ****************************
*
*     C(1:nfft3d,ms) = Sum(ispin,ne,nbrill) A(1:nfft3d,n,nbrill)**2*brillioun_weight()
*     uses - C3dB_cr_aSqrpy

      subroutine  Stfl_hr_SumSqr(ispin,ne,nfft3d,neall,nbrill,A,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      real*8     C(nfft3d,2)

*     **** local variables ****
      integer nb,ms,n,n1(2),n2(2)
      real*8  weight

*     **** local functions ****
      real*8   brillioun_weight
      external brillioun_weight

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      call C3dB_r_Zero(1,C)
      if (ispin.gt.1) call C3dB_r_Zero(1,C(1,ispin))

      do nb=1,nbrill
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
*     *                		 *
*     *    	Stfl_gqg_Mul   	 *
*     *                		 *
*     ****************************
*
*     C = A*Q
*     uses - Cram_c_ZMul, Cram_cc_zaxpy

      subroutine Stfl_gqg_Mul(ispin,ne,nfft3d,neall,nbrill,n,A,Q,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill,n
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 Q(n,n,2,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer j,k,nb,ms,jj,kk
      integer n1(2),n2(2)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      call Stfl_g_Zero(ispin,ne,nfft3d,neall,nbrill,C)
      do nb=1,nbrill
      do ms=1,ispin
      do k=n1(ms),n2(ms)
         kk = k-n1(ms)+1
         do j=n1(ms),n2(ms)
           jj = j-n1(ms)+1
           call Cram_cc_zaxpy(nb,Q(jj,kk,ms,nb),A(1,j,nb),C(1,k,nb))
         end do
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

      subroutine Stfl_gg_dScale(ispin,ne,nfft3d,neall,nbrill,scale,A,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      real*8     scale
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,(ne(1)+ne(2))
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

      subroutine Stfl_hh_dScale(ispin,ne,nfft3d,neall,nbrill,scale,A,C)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      real*8 scale
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 C(nfft3d,neall,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,(ne(1)+ne(2))
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
*     B = FFT(-1)[A]
*     uses - C3dB_cr_fft3b,Cram_c_unpack

      subroutine Stfl_gh_fftb(ispin,ne,npack1,nfft3d,neall,nbrill,A,B)
      implicit none
      integer ispin,ne(2)
      integer npack1,nfft3d,neall,nbrill
      complex*16 A(npack1,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,(ne(1)+ne(2))
        call Cram_c_Copy(nb,A(1,k,nb),B(1,k,nb))
        call Cram_c_unpack(nb,B(1,k,nb))
        call C3dB_cr_fft3b(1,B(1,k,nb))
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

      subroutine Stfl_hg_fftf(ispin,ne,nfft3d,neall,nbrill,A)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)

*     **** local variables ****
      integer k,nb

      do nb=1,nbrill
      do k=1,(ne(1)+ne(2))
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

      subroutine Stfl_ggq_dot(ispin,ne,nfft3d,neall,nbrill,n,A,B,Q)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill,n
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      complex*16 Q(n,n,2,nbrill)

*     **** local variables ****
      integer j,k,np,nb,jj,kk
      integer ms,n1(2),n2(2)
      complex*16 sum

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      call Parallel_np(np)
      call dcopy(4*n*n*nbrill,0.0d0,0,Q,1)

      do nb=1,nbrill
      do ms=1,ispin
      do k=n1(ms),n2(ms)
         kk = k-n1(ms)+1
         do j=n1(ms),n2(ms)
           jj = j-n1(ms)+1
           call Cram_cc_izdot(nb,A(1,j,nb),B(1,k,nb),sum) 
           Q(jj,kk,ms,nb) = sum
         end do
      end do
      end do
      end do

      if (np.gt.1) call C3dB_Vector_SumAll(4*n*n*nbrill,Q)
      return
      end

*     ********************************
*     *                	    	     *
*     *    	Stfl_gg_traceall     *
*     *                		     *
*     ********************************
      subroutine Stfl_gg_traceall(ispin,ne,nfft3d,neall,nbrill,A,B,sum)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      real*8 sum

#include "mafdecls.fh"

*     **** local variables ****
      logical value
      integer nb,ms,indx
      integer tmp1(2)

*     **** external functions ****
      real*8   brillioun_weight
      external brillioun_weight

      value =  MA_push_get(mt_dbl,(2*nbrill),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)

      call Stfl_gg_trace(ispin,ne,nfft3d,neall,nbrill,
     >                   A,B,dbl_mb(tmp1(1)))

      sum = 0.0d0
      do nb=1,nbrill
      do ms=1,ispin
         indx = ms + (nb-1)*2
         sum = sum + brillioun_weight(nb)*dbl_mb(tmp1(1)+indx-1)
      end do
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack memory',0)
      return
      end

*     ********************************
*     *                	    	     *
*     *    	Stfl_gg_trace        *
*     *                		     *
*     ********************************
*
*     uses - 

      subroutine Stfl_gg_trace(ispin,ne,nfft3d,neall,nbrill,A,B,MM)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      real*8     MM(2,nbrill)

*     **** local variables ****
      integer ms,j,nb,np,n1(2),n2(2)
      real*8 sum

      call Parallel_np(np)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      call dcopy(2*nbrill,0.0d0,0,MM,1)
      do nb=1,nbrill
      do ms=1,ispin
         do j=n1(ms),n2(ms)
            call Cram_cc_idot(nb,A(1,j,nbrill),B(1,j,nbrill),sum)
            MM(ms,nb) = MM(ms,nb) + sum
         end do
      end do
      end do

      if (np.gt.1) call C3dB_Vector_SumAll(2*nbrill,MM)
      return
      end

*     ***********************************
*     *                                 *
*     *            Stfl_q_diag          *
*     *                                 *
*     ***********************************
      subroutine Stfl_q_diag(ispin,ne,neall,nbrill,n,S,V)
      implicit none
      integer ispin,ne(2)
      integer neall,nbrill,n
      real*8     S(neall,nbrill)
      complex*16 V(n,n,2,nbrill)

#include "mafdecls.fh"

*     **** local variables ****
      logical value
      integer nb,ms
      integer rwork(2),work(2),lwork,ierr

      lwork = 2*ne(1)*ne(1)
      value = MA_push_get(mt_dbl,lwork,
     >                    'rwork',rwork(2),rwork(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,lwork,
     >                    'work',work(2),work(1))
      if (.not. value)
     >   call errquit('Stfl_q_diag:out of stack memory',0)

      do nb=1,nbrill
      do ms=1,ispin
         call ZHEEV('V','L',
     >              ne(ms),
     >              V(1,1,ms,nb),
     >              n,
     >              S(1+(ms-1)*n,nb),
     >              dcpl_mb(work(1)),
     >              lwork,
     >              dbl_mb(rwork(1)),
     >              ierr)
      end do
      end do

      value = MA_pop_stack(work(2))
      value = value.and.MA_pop_stack(rwork(2))
      if (.not. value)
     > call errquit('Stfl_q_diag:error popping stack',0)
      return
      end

*     ***********************************
*     *                                 *
*     *            Stfl_SVD             *
*     *                                 *
*     ***********************************

      subroutine Stfl_SVD(ispin,ne,nfft3d,neall,nbrill,n,A,U,S,V)
      implicit none
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill,n
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 U(nfft3d,neall,nbrill)
      real*8     S(neall,nbrill)
      complex*16 V(n,n,2,nbrill)

#include "mafdecls.fh"

*     **** local variables ****
      logical value
      integer i,nb
      integer tmp2(2)


*     ***** generate V and Sigma^2 ****
      call Stfl_ggq_dot(ispin,ne,nfft3d,neall,nbrill,n,A,A,V)
      call Stfl_q_diag(ispin,ne,neall,nbrill,n,S,V)


*     **** generate U*Sigma ****
      call Stfl_gqg_Mul(ispin,ne,nfft3d,neall,nbrill,n,A,V,U)


*     **** normalize U*sigma, and sqrt(S^2) ****
      value =  MA_push_get(mt_dbl,(neall),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

      do nb=1,nbrill

*        **** normalize U*sigma ****
         do i=1,neall
            call Cram_cc_idot(nb,U(1,i,nb),U(1,i,nb),
     >                        dbl_mb(tmp2(1)+i-1))
         end do
         call C3dB_Vector_SumAll(neall,dbl_mb(tmp2(1)))
      
         do i=1,neall
            dbl_mb(tmp2(1)+i-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+i-1))
   
            call Cram_c_SMul(nb,dbl_mb(tmp2(1)+i-1),
     >                       U(1,i,nb),U(1,i,nb))
         end do

*        **** sqrt(S*S) ****
         do i=1,neall
           if (S(i,nb).lt.0.0d0) S(i,nb) = dabs(S(i,nb))
           S(i,nb) = dsqrt(S(i,nb))
         end do

      end do
      value = MA_pop_stack(tmp2(2))
      if (.not. value) call errquit('error popping stack memory',0)
       
      return
      end


