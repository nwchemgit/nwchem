*
* $Id: geodesic2.F,v 1.4 2003-06-30 01:23:11 bylaska Exp $
*

*     ***********************************
*     *					*
*     *		geodesic2_init		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic2_init()
      implicit none

#include "mafdecls.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ne(2)

*     **** external functions **** 
      integer  psi_ne
      external psi_ne

      call Pack_npack(1,npack1)
      nemax = psi_ne(1)+psi_ne(2)
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)

      value = MA_alloc_get(mt_dcpl,npack1*nemax,
     >                     'Hold',Hold(2),Hold(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*nemax,
     >                     'Q',Q(2),Q(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                     'R',R(2),R(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                     'A',A(2),A(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                     'V',V(2),V(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,4*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                     'W',W(2),W(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nemax,
     >                     'S',S(2),S(1))
      if (.not. value) call errquit('geodesic2_init:allocating heap',0)

      return
      end

*     ***********************************
*     *					*
*     *		geodesic2_finalize	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic2_finalize()
      implicit none

#include "mafdecls.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value

      value =           MA_free_heap(S(2))   
      value = value.and.MA_free_heap(W(2))
      value = value.and.MA_free_heap(V(2))
      value = value.and.MA_free_heap(A(2))
      value = value.and.MA_free_heap(R(2))
      value = value.and.MA_free_heap(Q(2))
      value = value.and.MA_free_heap(Hold(2))
      if (.not. value) 
     >  call errquit('geodesic2_finalize:freeing heap memory',0)

      return
      end



*     ***********************************
*     *					*
*     *		geodesic2_start		*
*     *					*
*     ***********************************
*
*     This routine determines the pxp matrices R and YA, and
* the orthogonal nxp matrix Q.   Q and R are determined from 
* the QR decomposition of the projected direction (I-YY^t)H, and
* YH is defined as the Lagrange Multiplier pxp matrix Y^tH.  
*
*     Uses - geodesic2 common block
*
      subroutine geodesic2_start(Y,H,max_sigma,dE)
      implicit none
      complex*16 Y(*)
      complex*16 H(*)
      real*8     max_sigma,dE

#include "mafdecls.fh"
#include "geodesic2.fh"


*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ms,n,ispin,ne(2)
      integer shift,shift2,i,j
      integer T(2)

*     **** external functions **** 
      integer  psi_ispin,psi_ne
      real*8   electron_eorbit
      external psi_ispin,psi_ne
      external electron_eorbit

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)
      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)


*     **** allocate tmp space ****
      value = MA_push_get(mt_dbl,4*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                              'T',T(2),T(1))
      if (.not. value) 
     >   call errquit('geodesic2_start: pushing stack',0)


*     **** Hold <-- H ****
      call Grsm_gg_Copy(npack1,nemax,H,dcpl_mb(Hold(1)))

*     **** calculate A=<Y|H> ****
      do ms=1,ispin
c        call Grsm_ggm_dot(npack1,ne(ms),
c    >                     Y(1+(ms-1)*ne(1)*npack1),
c    >                     H(1+(ms-1)*ne(1)*npack1),
c    >                     dbl_mb(A(1)+(ms-1)*ne(1)*ne(1)))
         call Pack_ccm_dot(1,ne(ms),
     >                     Y(1+(ms-1)*ne(1)*npack1),
     >                     H(1+(ms-1)*ne(1)*npack1),
     >                     dbl_mb(A(1)+(ms-1)*ne(1)*ne(1)))
      end do


*     **** calculate Q=(I-YYt)H - should be necessary but just in case ****
      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
c        call Grsm_ggm_dot(npack1,n,                  ! R=Yt*H
c    >                     Y(1+(ms-1)*ne(1)*npack1),
c    >                     H(1+(ms-1)*ne(1)*npack1),
c    >                     dbl_mb(R(1)+shift))
         call Pack_ccm_dot(1,n,                  ! R=Yt*H
     >                     Y(1+(ms-1)*ne(1)*npack1),
     >                     H(1+(ms-1)*ne(1)*npack1),
     >                     dbl_mb(R(1)+shift))

         call Grsm_gmg_Mul(npack1,n,                  ! Q=Y*R
     >                     Y(1+(ms-1)*ne(1)*npack1),
     >                     dbl_mb(R(1)+shift),
     >                     dcpl_mb(Q(1)+(ms-1)*ne(1)*npack1))
      end do
      call Grsm_ggg_Sub(npack1,(ne(1)+ne(2)),         ! Q=H-Q = H-Y*R =H-Y*Yt*H
     >                  H,
     >                  dcpl_mb(Q(1)),
     >                  dcpl_mb(Q(1)))

*     **** calculate QR using Modified Gram-Schmidt ****
      do ms=1,ispin
         shift  = (ms-1)*ne(1)*ne(1)
         shift2 = (ms-1)*ne(1)*npack1
         call Grsm_g_QR(npack1,ne(ms),
     >                  dcpl_mb(Q(1)+shift2),
     >                  dbl_mb(R(1) +shift))
      end do



*     **** generate T ****
      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         call geodesic2_generate_T(ne(ms),
     >                             dbl_mb(A(1)+shift),
     >                             dbl_mb(R(1)+shift),
     >                             dbl_mb(T(1)+4*shift))
      end do

*     **** Factor T--> V,W,and S ****
      do ms=1,ispin
         shift  = (ms-1)*4*ne(1)*ne(1)
         shift2 = (ms-1)*2*ne(1)
         n = 2*ne(ms)
         call Factor_Skew(n,
     >                    dbl_mb(T(1)+shift),
     >                    dbl_mb(V(1)+shift),
     >                    dbl_mb(W(1)+shift),
     >                    dbl_mb(S(1)+shift2))
      end do
      

*     **** calculate dE ****
      dE = 2.0d0*electron_eorbit(H)

    
*     **** deallocate tmp space ****     
      value = MA_pop_stack(T(2))
      if (.not. value) 
     > call errquit('geodesic2_start:popping stack',1)

      
      call nwpw_timing_end(10)

      return
      end

*     ***********************************
*     *					*
*     *	      geodesic2_generate_T	*
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

      subroutine geodesic2_generate_T(n,A,R,T)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 R(n,n)
      real*8 T(2*n,2*n)

*     **** local variables ****
      integer i,j

      call dcopy(4*n*n,0.0d0,0,T,1)
      
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
         T(i,j+n) = -R(j,i)
      end do
      end do

      return
      end 


*     ***********************************
*     *					*
*     *		geodesic2_get		*
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
      subroutine geodesic2_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "mafdecls.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,ispin,ne(2)
      integer ms,n,j,k,shift1,shift2
      integer MM(2),NN(2)

*     **** external functions ****
      integer  psi_ne,psi_ispin
      external psi_ne,psi_ispin


      call Pack_npack(1,npack1)
      ne(1) = psi_ne(1) 
      ne(2) = psi_ne(2)
      ispin = psi_ispin()

*     **** allocate tmp space ****
      value = MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),'MM',
     >                    MM(2),MM(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),'NN',
     >                    NN(2),NN(1))
      if (.not. value) 
     >   call errquit('geodesic2_get: pushing stack',0)

      call geodesic2_get_MandN(t,dbl_mb(MM(1)),dbl_mb(NN(1)))


      do ms=1,ispin
        n = ne(ms)
        if (n.eq.0) go to 100  !*** ferromagnetic case ***
        shift1 = (ms-1)*ne(1)
        shift2 = (ms-1)*ne(1)*ne(1)


        call DGEMM('N','N',2*npack1,n,n,
     >              (1.0d0),
     >              Yold(1+shift1*npack1),2*npack1,
     >              dbl_mb(MM(1)+shift2),n,
     >              (0.0d0),
     >              Ynew(1+shift1*npack1),2*npack1)

        call DGEMM('N','N',2*npack1,n,n,
     >              (1.0d0),
     >              dcpl_mb(Q(1)+shift1*npack1),2*npack1,
     >              dbl_mb(NN(1)+shift2),n,
     >              (1.0d0),
     >              Ynew(1+shift1*npack1),2*npack1)
c       do k=1,n
c         call dcopy(2*npack1,0.0d0,0,Ynew(1+(k-1+shift1)*npack1),1)
c         do j=1,n
c           call Pack_cc_daxpy(1,
c    >                         dbl_mb(MM(1)+j-1+(k-1)*n+shift2),
c    >                         Yold(1+(j-1+shift1)*npack1),
c    >                         Ynew(1+(k-1+shift1)*npack1))
c           call Pack_cc_daxpy(1,
c    >                         dbl_mb(NN(1)+j-1+(k-1)*n+shift2),
c    >                         dcpl_mb(Q(1)+(j-1+shift1)*npack1),
c    >                         Ynew(1+(k-1+shift1)*npack1))
c         end do
c       end do

 100    continue

      end do

*     **** deallocate tmp space ****     
      value =           MA_pop_stack(NN(2))
      value = value.and.MA_pop_stack(MM(2))
      if (.not. value) 
     > call errquit('geodesic2_get:popping stack',1)

      return
      end


*     ***********************************
*     *					*
*     *		geodesic2_transport	*
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
      subroutine geodesic2_transport(t,Yold,Hnew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Hnew(*)

#include "mafdecls.fh"
#include "geodesic2.fh"


*     **** local variables ****
      logical value
      integer npack1,ispin,ne(2)
      integer ms,n,i,j,k,shift,shift1,shift2
      integer MM(2),NN(2),TT(2)
      real*8  sum

*     **** external functions ****
      integer  psi_ne,psi_ispin
      external psi_ne,psi_ispin


      call Pack_npack(1,npack1)
      ne(1) = psi_ne(1) 
      ne(2) = psi_ne(2)
      ispin = psi_ispin()

*     **** allocate tmp space ****
      value = MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),'MM',
     >                    MM(2),MM(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),'NN',
     >                    NN(2),NN(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),'TT',
     >                    TT(2),TT(1))
      if (.not. value) 
     >   call errquit('geodesic2_transport: pushing stack',0)

      call geodesic2_get_MandN(t,dbl_mb(MM(1)),dbl_mb(NN(1)))

*     **** TT(t) = -R^t*NN(t) ****
      do ms=1,ispin
        n=ne(ms)
        shift2 = (ms-1)*ne(1)*ne(1)

        call DGEMM('T','N',n,n,n,
     >              (-1.0d0),
     >              dbl_mb(R(1) +shift2),n,
     >              dbl_mb(NN(1)+shift2),n,
     >              (0.0d0),
     >              dbl_mb(TT(1)+shift2),n)
c       do j=1,n
c       do i=1,n
c          sum =0.0d0
c          shift = (i-1) + (j-1)*n + (ms-1)*ne(1)*ne(1)
c          do k=1,n
c             shift1 = (k-1) + (i-1)*n + (ms-1)*ne(1)*ne(1)
c             shift2 = (k-1) + (j-1)*n + (ms-1)*ne(1)*ne(1)
c             sum = dbl_mb(R(1)+shift1)*dbl_mb(NN(1)+shift2)
c          end do
c          dbl_mb(TT(1)+shift) = -sum
c       end do
c       end do

      end do

*     *** Hnew <-- Hold*M(t) + Yold*TT(t) ***
      do ms=1,ispin
        n = ne(ms)
        if (n.eq.0) go to 100  !*** ferromagnetic case ***
        shift1 = (ms-1)*ne(1)
        shift2 = (ms-1)*ne(1)*ne(1)

        call DGEMM('N','N',2*npack1,n,n,
     >              (1.0d0),
     >              dcpl_mb(Hold(1)+shift1*npack1),2*npack1,
     >              dbl_mb(MM(1)+shift2),n,
     >              (0.0d0),
     >              Hnew(1+shift1*npack1),2*npack1)

        call DGEMM('N','N',2*npack1,n,n,
     >              (1.0d0),
     >              Yold(1+shift1*npack1),2*npack1,
     >              dbl_mb(TT(1)+shift2),n,
     >              (1.0d0),
     >              Hnew(1+shift1*npack1),2*npack1)
c       do k=1,n         
c         call dcopy(2*npack1,0.0d0,0,Hnew(1+(k-1+shift1)*npack1),1)
c         do j=1,n
c           call Pack_cc_daxpy(1,
c    >                         dbl_mb(MM(1)+j-1+(k-1)*n+shift2),
c    >                         dcpl_mb(Hold(1)+(j-1+shift1)*npack1),
c    >                         Hnew(1+(k-1+shift1)*npack1))
c           call Pack_cc_daxpy(1,
c    >                         dbl_mb(TT(1)+j-1+(k-1)*n+shift2),
c    >                         Yold(1+(k-1+shift1)*npack1),
c    >                         Hnew(1+(k-1+shift1)*npack1))
c         end do
c       end do

 100    continue

      end do

*     **** deallocate tmp space ****    
      value =           MA_pop_stack(TT(2)) 
      value = value.and.MA_pop_stack(NN(2))
      value = value.and.MA_pop_stack(MM(2))
      if (.not. value) 
     > call errquit('geodesic2_transport:popping stack',1)


      return
      end

      subroutine geodesic2_checkMN(name,n,MM,NN)
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
*     *	      geodesic2_get_MandN	*
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
      subroutine geodesic2_get_MandN(t,M,N)
      implicit none
      real*8 t
      real*8 M(*),N(*)

#include "mafdecls.fh"
#include "geodesic2.fh"


*     **** local variables ****
      logical value
      integer ispin,ne(2)
      integer ms,i,j,nn,shift,shift2,shift3
      integer RR(2),tmpv(2),tmpw(2)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne

      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)

*     **** allocate tmp space ****
      value =           MA_push_get(mt_dbl,4*(ne(1)*ne(1)+ne(2)*ne(2)),
     >                              'RR',RR(2),RR(1))
      value = value.and.MA_push_get(mt_dbl,4*(ne(1)*ne(1)),
     >                              'tmpv',tmpv(2),tmpv(1))
      value = value.and.MA_push_get(mt_dbl,4*(ne(1)*ne(1)),
     >                              'tmpw',tmpw(2),tmpw(1))
      if (.not. value) 
     >   call errquit('geodesic2_get_MandN: pushing stack',0)

      do ms=1,ispin
         shift  = (ms-1)*4*ne(1)*ne(1)
         shift2 = (ms-1)*2*ne(1)
         nn = ne(ms)
         call Rotate_Skew(2*nn,
     >                    dbl_mb(V(1)+shift),
     >                    dbl_mb(W(1)+shift),
     >                    dbl_mb(S(1)+shift2),
     >                    t,
     >                    dbl_mb(RR(1)+shift),  
     >                    dbl_mb(tmpv(1)),dbl_mb(tmpw(1)))               
         do j=1,nn
         do i=1,nn
            shift  = (i-1)    +(j-1)*nn   + (ms-1)*ne(1)*ne(1)
            shift2 = (i-1)    +(j-1)*2*nn + (ms-1)*4*ne(1)*ne(1)
            shift3 = (i-1+nn) +(j-1)*2*nn + (ms-1)*4*ne(1)*ne(1)
            M(1+shift) = dbl_mb(RR(1)+shift2)
            N(1+shift) = dbl_mb(RR(1)+shift3)
         end do
         end do
      end do


*     **** deallocate tmp space ****     
      value =           MA_pop_stack(tmpw(2))
      value = value.and.MA_pop_stack(tmpv(2))
      value = value.and.MA_pop_stack(RR(2))
      if (.not. value) 
     > call errquit('geodesic2_get_MandN:popping stack',1)

      return
      end
