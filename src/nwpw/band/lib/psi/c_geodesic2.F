*
* $Id$
*

*     ***********************************
*     *					*
*     *	       c_geodesic2_init		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine c_geodesic2_init()
      implicit none

#include "bafdecls.fh"
#include "c_geodesic2.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer npack1,neall,nemax,nbrillq

*     **** external functions **** 
      integer  cpsi_ne,cpsi_neq,cpsi_nbrillq,Pneb_w_size,Pneb_4w_size
      external cpsi_ne,cpsi_neq,cpsi_nbrillq,Pneb_w_size,Pneb_4w_size
      integer  cpsi_data_alloc
      external cpsi_data_alloc

      call Cram_max_npack(npack1)
      neall  = cpsi_neq(1)+cpsi_neq(2)
      nemax  = cpsi_ne(1)+cpsi_ne(2)
      nbrillq = cpsi_nbrillq()
      
      Hold_tag  = cpsi_data_alloc(nbrillq,neall,2*npack1)
      Q_tag     = cpsi_data_alloc(nbrillq,neall,2*npack1)
      R_tag     = cpsi_data_alloc(nbrillq,1,2*Pneb_w_size(0,0))
      A_tag     = cpsi_data_alloc(nbrillq,1,2*Pneb_w_size(0,0))
      U_tag     = cpsi_data_alloc(nbrillq,1,2*Pneb_4w_size(0,0))
      S_tag     = cpsi_data_alloc(nbrillq,2*nemax,1)

      return
      end

*     ***********************************
*     *					*
*     *		c_geodesic2_finalize	*
*     *					*
*     ***********************************
*
*     Uses - c_geodesic2 common block
*

      subroutine c_geodesic2_finalize()
      implicit none

#include "bafdecls.fh"
#include "c_geodesic2.fh"
#include "errquit.fh"

*     **** local variables ****


      call cpsi_data_dealloc(Hold_tag)
      call cpsi_data_dealloc(Q_tag)
      call cpsi_data_dealloc(R_tag)
      call cpsi_data_dealloc(A_tag)
      call cpsi_data_dealloc(U_tag)
      call cpsi_data_dealloc(S_tag)
      return
      end



*     ***********************************
*     *					*
*     *	      c_geodesic2_start		*
*     *					*
*     ***********************************
*
*     This routine determines the pxp matrices R and YA, and
* the orthogonal nxp matrix Q.   Q and R are determined from 
* the QR decomposition of the projected direction (I-YY^t)H, and
* YH is defined as the Lagrange Multiplier pxp matrix Y^tH.  
*
*     Uses - c_geodesic2 common block
*
      subroutine c_geodesic2_start(Y_tag,H_tag,max_sigma,dE)
      implicit none
      integer Y_tag,H_tag
      real*8     max_sigma,dE

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,neall,nbrillq
      integer ms,n,ispin,ne(2),nb
      integer shift,shift2,i,j
      integer T_tag

*     **** external functions **** 
      integer  Pneb_4w_size,Pneb_w_size
      external Pneb_4w_size,Pneb_w_size
      integer  cpsi_ispin,cpsi_neq,cpsi_nbrillq,cpsi_ne
      external cpsi_ispin,cpsi_neq,cpsi_nbrillq,cpsi_ne
      integer  cpsi_data_push_stack,cpsi_data_get_allptr
      external cpsi_data_push_stack,cpsi_data_get_allptr
      !real*8   c_electron_eorbit_noocc
      !external c_electron_eorbit_noocc
      real*8   c_electron_eorbit
      external c_electron_eorbit
      integer  cpsi_data_nsize
      external cpsi_data_nsize
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk

      call nwpw_timing_start(10)
      call Cram_max_npack(npack1)
      neall = cpsi_neq(1) + cpsi_neq(2)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      nbrillq = cpsi_nbrillq()


*     **** allocate tmp space ****
      T_tag     = cpsi_data_push_stack(nbrillq,1,2*Pneb_4w_size(0,0))

*     **** Hold <-- H ****
      call Parallel_shared_vector_copy(.true.,
     >                    nbrillq*2*npack1*neall,
     >                    dbl_mb(cpsi_data_get_allptr(H_tag)),
     >                    dbl_mb(cpsi_data_get_allptr(Hold_tag)))


*     **** calculate A=<Y|H>??? ****
      call Pneb_ffw_Multiply(0,0,
     >                       dbl_mb(cpsi_data_get_allptr(Y_tag)),
     >                       dbl_mb(cpsi_data_get_allptr(H_tag)),npack1,
     >                       dbl_mb(cpsi_data_get_allptr(A_tag)))


*     **** calculate Q=(I-YYt)H - should not be necessary but just in case ****
      call Pneb_fwf_Multiply(0,0,dcmplx(1.0d0,0.0d0),
     >                       dbl_mb(cpsi_data_get_allptr(Y_tag)),npack1,
     >                       dbl_mb(cpsi_data_get_allptr(A_tag)),
     >                       dcmplx(0.0d0,0.0d0),
     >                       dbl_mb(cpsi_data_get_allptr(Q_tag)))

      call DAXPY_OMP(nbrillq*2*npack1*neall,(-1.0d0),
     >                dbl_mb(cpsi_data_get_allptr(H_tag)),1,
     >                dbl_mb(cpsi_data_get_allptr(Q_tag)),1)
      call DSCAL_OMP(nbrillq*2*npack1*neall,(-1.0d0),
     >               dbl_mb(cpsi_data_get_allptr(Q_tag)),1)

*     **** calculate QR using Modified Gram-Schmidt ****
      !call Pneb_fw_QR(0,0,dcpl_mb(Q(1)),npack1,dcpl_mb(R(1)))
      call Pneb_fw_QR(0,0,
     >                dbl_mb(cpsi_data_get_allptr(Q_tag)),npack1,
     >                dbl_mb(cpsi_data_get_allptr(R_tag)))


*     **** generate T from A and R ****
*       -     -
*  T = |A, -R^t|
*      |R,  0  |
*       -     -
      call Pneb_AR_to_4w(0,0,
     >                   dbl_mb(cpsi_data_get_allptr(A_tag)),
     >                   dbl_mb(cpsi_data_get_allptr(R_tag)),
     >                   dbl_mb(cpsi_data_get_allptr(T_tag)))

      
*     **** Factor T--> V,W,and S ****
      call Pneb_4w_FactorSkew(0,0,
     >                        dbl_mb(cpsi_data_get_allptr(T_tag)),
     >                        dbl_mb(cpsi_data_get_allptr(U_tag)),
     >                        dbl_mb(cpsi_data_get_allptr(S_tag)),
     >                        max_sigma)

*     **** calculate dE ****
      !dE = 2.0d0*c_electron_eorbit_noocc(H)
      dE = 2.0d0*c_electron_eorbit(H_tag)

      call Pneb_ffw_Multiply(0,0,
     >                       dbl_mb(cpsi_data_get_allptr(Y_tag)),
     >                       dbl_mb(cpsi_data_get_allptr(Y_tag)),npack1,
     >                       dbl_mb(cpsi_data_get_allptr(R_tag)))
      call Pneb_ffw_Multiply(0,0,
     >                       dbl_mb(cpsi_data_get_allptr(Q_tag)),
     >                       dbl_mb(cpsi_data_get_allptr(Q_tag)),npack1,
     >                       dbl_mb(cpsi_data_get_allptr(T_tag)))

c      do nb=1,nbrillq
c       write(*,*) "nb=",nb
c       call write_UU("A=",ne(1),
c     >     dbl_mb(cpsi_data_get_allptr(A_tag)+(nb-1)*2*ne(1)*ne(1)))
c       call write_UU("YtY=",ne(1),
c     >     dbl_mb(cpsi_data_get_allptr(R_tag)+(nb-1)*2*ne(1)*ne(1)))
c       call write_UU("QtQ=",ne(1),
c     >     dbl_mb(cpsi_data_get_allptr(T_tag)+(nb-1)*2*ne(1)*ne(1)))
c       write(*,*) "Sigma=",
c     >  (dbl_mb(cpsi_data_get_allptr(S_tag)+(nb-1)*ne(1)+i-1),i=1,ne(1))
c       write(*,*)
c      end do


*     **** deallocate tmp space ****     
      call cpsi_data_pop_stack(T_tag)

      call nwpw_timing_end(10)
      return
      end

      subroutine write_UU(aa,N,U)
      implicit none
      character*(*) aa
      integer N
      complex*16 U(N,N)

      integer i,j

      write(*,*) "matrix:",aa,":"
      do i=1,N
         write(*,*) (U(i,j), j=1,N)
      end do
      return
      end 


*     ***********************************
*     *					*
*     *	      c_geodesic2_generate_T	*
*     *					*
*     ***********************************
*
*     This routine determines T.  T is defined
* to be a 2nx2n skew symmetric matrix.
*
*       -     -
*  T = |A, -R^t|
*      |R,  0  |
*       -     -
*
*     Entry - n: dimension of matrices A and R
*             A: an nxn skew symmetric matrix
*             R: an nxn matrix
*     Exit - T: a 2nx2n skew symmetric matrix

      subroutine c_geodesic2_generate_T(n,A,R,T)
      implicit none
      integer n
      complex*16 A(n,n)
      complex*16 R(n,n)
      complex*16 T(2*n,2*n)

*     **** local variables ****
      integer i,j

      call Parallel_shared_vector_zero(.true.,8*n*n,T)
      
*     **** copy A to upper-left of T ****
      do j=1,n
      do i=1,n
         T(i,j) = A(i,j)
      end do
      end do

*     **** copy R to lower-left of T ****
      do j=1,n
      do i=1,n
         T(i+n,j) = R(i,j)
      end do
      end do

*     **** copy -R^t to upper-right of T ****
      do j=1,n
      do i=1,n
         T(i,j+n) = -dconjg(R(j,i))
      end do
      end do

      return
      end 


*     ***********************************
*     *					*
*     *		c_geodesic2_get		*
*     *					*
*     ***********************************
*
*     This routine calculates 
*
*   Ynew = Yold*M(t) + Q*N(t)
*
*   where
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*
      subroutine c_geodesic2_get(t,Yold_tag,Ynew_tag)
      implicit none
      real*8     t
      integer Yold_tag
      integer Ynew_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "c_geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,nbrillq
      integer MM_tag,NN_tag

*     **** external functions ****
      integer  cpsi_nbrillq,cpsi_ne
      external cpsi_nbrillq,cpsi_ne
      integer  cpsi_data_push_stack,cpsi_data_get_allptr
      external cpsi_data_push_stack,cpsi_data_get_allptr
      integer  Pneb_w_size
      external Pneb_w_size

      call Cram_max_npack(npack1)
      nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      MM_tag  = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,0))
      NN_tag  = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,0))

      call c_geodesic2_get_MandN(t,MM_tag,NN_tag)

      call Pneb_fwf_Multiply(0,0,dcmplx(1.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Yold_tag)),npack1,
     >             dbl_mb(cpsi_data_get_allptr(MM_tag)),
     >             dcmplx(0.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Ynew_tag)))

      call Pneb_fwf_Multiply(0,0,dcmplx(1.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Q_tag)),npack1,
     >             dbl_mb(cpsi_data_get_allptr(NN_tag)),
     >             dcmplx(1.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Ynew_tag)))

*     **** deallocate tmp space ****     
      call cpsi_data_pop_stack(NN_tag)
      call cpsi_data_pop_stack(MM_tag)
      return
      end


*     ***********************************
*     *					*
*     *		c_geodesic2_transport	*
*     *					*
*     ***********************************
*
*     This routine calculates 
*
*   Hnew = Hold*M(t)    + Yold*R^t*N(t)
*
*   where
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*
      subroutine c_geodesic2_transport(t,Yold_tag,Hnew_tag)
      implicit none
      real*8     t
      integer Yold_tag
      integer Hnew_tag

#include "bafdecls.fh"
#include "c_geodesic2.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer npack1,nbrillq
      integer MM_tag,NN_tag,TT_tag

*     **** external functions ****
      integer  cpsi_nbrillq
      external cpsi_nbrillq
      integer  cpsi_data_push_stack,cpsi_data_get_allptr
      external cpsi_data_push_stack,cpsi_data_get_allptr
      integer  Pneb_4w_size,Pneb_w_size
      external Pneb_4w_size,Pneb_w_size

      call Cram_max_npack(npack1)
      nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      MM_tag  = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,0))
      NN_tag  = cpsi_data_push_stack(nbrillq,1,2*Pneb_w_size(0,0))
      TT_tag  = cpsi_data_push_stack(nbrillq,1,2*Pneb_4w_size(0,0))


      call c_geodesic2_get_MandN(t,MM_tag,NN_tag)

*     **** TT(t) = -R^t*NN(t) ****
      call Pneb_www_Multiply2(0,0,dcmplx(-1.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(R_tag)),
     >             dbl_mb(cpsi_data_get_allptr(NN_tag)),
     >             dcmplx(0.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(TT_tag)))


*     *** Hnew <-- Hold*M(t) + Yold*TT(t) ***
      call Pneb_fwf_Multiply(0,0,dcmplx(1.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Hold_tag)),npack1,
     >             dbl_mb(cpsi_data_get_allptr(MM_tag)),
     >             dcmplx(0.0d0,0.0d0),
     >             dbl_mb(cpsi_data_get_allptr(Hnew_tag)))

*     **** deallocate tmp space ****    
      call cpsi_data_pop_stack(TT_tag)
      call cpsi_data_pop_stack(NN_tag)
      call cpsi_data_pop_stack(MM_tag)

      return
      end

      subroutine c_geodesic2_checkMN(name,n,MM,NN)
      implicit none
      character*(*) name
      integer n
      real*8 MM(n,n),NN(n,n)

*     **** local variables ****
      integer i,j,k
      real*8 tmp(n,n),sum

      do j=1,n
      do i=1,n
        sum = 0.0d0
        do k=1,n
          sum = sum + MM(k,i)*NN(k,j)
        end do
        tmp(i,j) = sum
      end do
      end do

      write(*,*)
      write(*,*) "checkMN:",name
      do i=1,n
        write(*,'(22F8.3)') (tmp(i,j),j=1,n)
      end do
      write(*,*)
      return
      end

*     ***********************************
*     *					*
*     *	      c_geodesic2_get_MandN	*
*     *					*
*     ***********************************
*
*   This routine returns 
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*   where 
* 
*      T =  U*Sigma*U^H, with U=(V+iW)
*
*      is a skew matrix that is decomposed into V,W,and Sigma
*      
      subroutine c_geodesic2_get_MandN(t,M_tag,N_tag)
      implicit none
      real*8 t
      integer M_tag,N_tag

#include "bafdecls.fh"
#include "c_geodesic2.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer nbrillq
      integer AB_tag,RR_tag

*     **** external functions ****
      integer  cpsi_nbrillq,cpsi_ne
      external cpsi_nbrillq,cpsi_ne
      integer  cpsi_data_push_stack,cpsi_data_get_allptr
      external cpsi_data_push_stack,cpsi_data_get_allptr
      integer  Pneb_4w_size
      external Pneb_4w_size

      nbrillq = cpsi_nbrillq()

*     **** allocate tmp space ****
      AB_tag = cpsi_data_push_stack(nbrillq,1,2*Pneb_4w_size(0,0))
      RR_tag = cpsi_data_push_stack(nbrillq,1,2*Pneb_4w_size(0,0))

      call Pneb_4w_RotateSkew(0,0,t,
     >                        dbl_mb(cpsi_data_get_allptr(U_tag)),
     >                        dbl_mb(cpsi_data_get_allptr(S_tag)),
     >                        dbl_mb(cpsi_data_get_allptr(AB_tag)),
     >                        dbl_mb(cpsi_data_get_allptr(RR_tag)))

      call Pneb_4w_to_MN(0,0,
     >                   dbl_mb(cpsi_data_get_allptr(RR_tag)),
     >                   dbl_mb(cpsi_data_get_allptr(M_tag)),
     >                   dbl_mb(cpsi_data_get_allptr(N_tag)))


*     **** deallocate tmp space ****     
      call cpsi_data_pop_stack(RR_tag)
      call cpsi_data_pop_stack(AB_tag)

      return
      end
