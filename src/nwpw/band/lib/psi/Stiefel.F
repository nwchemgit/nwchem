*
* $Id$
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
ccccccccccccccc
      call Pneb_ff_copy(0,0,nfft3d,A,C)
      return
cccccccccccc
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

      call Pneb_fff_sum(0,0,nfft3d,A,B,C)
      return
c      do nb=1,nbrill
c      do i=1,(ne(1)+ne(2))
c          call Cram_cc_Sum(nb,A(1,i,nb),B(1,i,nb),C(1,i,nb))
c      end do
c      end do
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

      call Pneb_fff_sub(0,0,nfft3d,A,B,C)
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

      call Pneb_fff_sum(0,0,nfft3d,A,B,C)
 
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
*     *                          *
*     *       Stfl_hr_aSqrpy     *
*     *                          *
*     ****************************
*
*     C(1:nfft3d,ms) = Sum(ispin,ne) A(1:nfft3d,n)**2*weight(n)
*     uses - C3dB_cr_aSqrpy

      subroutine  Stfl_hr_aSqrpy(nfft3d,ne,weight,A,C)
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

      call Pneb_fwf_multiply(0,0,A,Q,nfft3d,C)
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

      call Pneb_ff_scale(0,0,nfft3d,scale,A,C)
      return
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

      call Pneb_ff_scale(0,0,nfft3d,scale,A,C)
      return
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

      call Pneb_ffw_Multiply(0,0,A,B,nfft3d,Q)
      return
      end

*     ********************************
*     *                	    	     *
*     *    	Stfl_gg_traceall     *
*     *                		     *
*     ********************************
      subroutine Stfl_gg_traceall(ispin,ne,nfft3d,neall,nbrill,A,B,sum)
      implicit none
#include "errquit.fh"
      integer ispin,ne(2)
      integer nfft3d,neall,nbrill
      complex*16 A(nfft3d,neall,nbrill)
      complex*16 B(nfft3d,neall,nbrill)
      real*8 sum

      call Pneb_ff_TraceAll(0,0,A,B,nfft3d,sum)
      return
      end

c*     ********************************
c*     *                	    	     *
c*     *    	Stfl_gg_trace        *
c*     *                		     *
c*     ********************************
c*
c*     uses - 
c
c      subroutine Stfl_gg_trace(ispin,ne,nfft3d,neall,nbrill,A,B,MM)
c      implicit none
c      integer ispin,ne(2)
c      integer nfft3d,neall,nbrill
c      complex*16 A(nfft3d,neall,nbrill)
c      complex*16 B(nfft3d,neall,nbrill)
c      real*8     MM(2,nbrill)
c
c      call Pneb_ff_Trace(0,0,A,B,nfft3d,MM)
c      return
c      end

c*     ***********************************
c*     *                                 *
c*     *            Stfl_q_diag          *
c*     *                                 *
c*     ***********************************
c      subroutine Stfl_q_diag(ispin,ne,neall,nbrill,n,S,V)
c      implicit none
c#include "errquit.fh"
c      integer ispin,ne(2)
c      integer neall,nbrill,n
c      real*8     S(neall,nbrill)
c      complex*16 V(n,n,2,nbrill)
c
c#include "bafdecls.fh"
c
c*     **** local variables ****
c      logical value
c      integer nb,ms
c      integer rwork(2),work(2),lwork,ierr
c      integer nbb,mb
c
c      mb=0
c      nbb=0
c      if (nbrill.eq.0) nbb=1
c      call Pneb_w_diag(mb,nbb,n,S,V)
c      return
c      end

c*     ***********************************
c*     *                                 *
c*     *            Stfl_SVD             *
c*     *                                 *
c*     ***********************************
c
c      subroutine Stfl_SVD(ispin,ne,nfft3d,neall,nbrill,n,A,U,S,V)
c      implicit none
c#include "errquit.fh"
c      integer ispin,ne(2)
c      integer nfft3d,neall,nbrill,n
c      complex*16 A(nfft3d,neall,nbrill)
c      complex*16 U(nfft3d,neall,nbrill)
c      real*8     S(neall,nbrill)
c      complex*16 V(n,n,2,nbrill)
c
c#include "bafdecls.fh"
c
c*     **** local variables ****
c      logical value,Pneb_is_init
c      integer i,nb
c      integer tmp2(2)
c      external Pneb_is_init
c
c      if (.not.Pneb_is_init()) then
c        write(0,*)"Pneb is not initialized!"
c      end if
c
c      call Pneb_SVD(0,0,nfft3d,A,U,S,V)
c      return
c
c*     ***** generate V and Sigma^2 ****
c      call Stfl_ggq_dot(ispin,ne,nfft3d,neall,nbrill,n,A,A,V)
c      call Stfl_q_diag(ispin,ne,neall,nbrill,n,S,V)
c
c
c*     **** generate U*Sigma ****
c      call Stfl_gqg_Mul(ispin,ne,nfft3d,neall,nbrill,n,A,V,U)
c
c
c*     **** normalize U*sigma, and sqrt(S^2) ****
c      value =  BA_push_get(mt_dbl,(neall),'tmp2',tmp2(2),tmp2(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c      do nb=1,nbrill
c
c*        **** normalize U*sigma ****
c         do i=1,neall
c            call Cram_cc_idot(nb,U(1,i,nb),U(1,i,nb),
c     >                        dbl_mb(tmp2(1)+i-1))
c         end do
c         call C3dB_Vector_SumAll(neall,dbl_mb(tmp2(1)))
c      
c         do i=1,neall
c            dbl_mb(tmp2(1)+i-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+i-1))
c   
c            call Cram_c_SMul(nb,dbl_mb(tmp2(1)+i-1),
c     >                       U(1,i,nb),U(1,i,nb))
c         end do
c
c*        **** sqrt(S*S) ****
c         do i=1,neall
c           if (S(i,nb).lt.0.0d0) S(i,nb) = dabs(S(i,nb))
c           S(i,nb) = dsqrt(S(i,nb))
c         end do
c
c      end do
c      value = BA_pop_stack(tmp2(2))
c      if (.not. value) call errquit('error popping stack memory',0,
c     &       MA_ERR)
c       
c      return
c      end


