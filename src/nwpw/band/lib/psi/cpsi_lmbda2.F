*
* $Id$
*
*     ****************************************************
*     *                                                  *
*     *                 cpsi_lmbda2                      *
*     *                                                  *
*     ****************************************************
      subroutine cpsi_lmbda2(ms,nb,npack1,
     >                     psi1,psi2,
     >                     dte,fweight,
     >                     tmp,ierr)
      implicit none
      integer ms,nb,npack1
      complex*16 psi1(*)
      complex*16 psi2(*)
      real*8     dte,fweight(*)
      complex*16 tmp(*)
      integer    ierr

#include "stdio.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
  
      logical value,failed
      integer it
      integer nn,ne(2)
      integer st1,st2
      integer A,B,C,U,D,Ba,Bs,fnm
      complex*16 zero,one,dtec

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      call Parallel_taskid(taskid)


      ierr = 0
      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
      dtec = dcmplx(dte,0.0d0)
*     **** set the matrices dimensions ****
      call Pneb_ne(ne)
      nn   = Pneb_w_size(0,nb)

      A    = 0*nn + 1
      B    = 1*nn + 1
      C    = 2*nn + 1
      Ba   = 3*nn + 1
      Bs   = 4*nn + 1
      fnm  = 5*nn + 1
      st1  = 6*nn + 1
      D    = 7*nn + 1

      U    = Bs
      st2  = B

      !call zcopy(8*nn,zero,0,tmp,1)
      call Parallel_shared_vector_zero(.true.,16*nn,tmp)
     
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::

        call Pneb_ffw_hermit_Multiply(ms,nb,psi2,psi2,npack1,tmp(A))
        call Pneb_ffw_Multiply(ms,nb,psi1,psi2,npack1,tmp(B))
        call Pneb_ffw_hermit_Multiply(ms,nb,psi1,psi1,npack1,tmp(C))

c        write(*,*) "A="
c        call prt_cmat(ne(ms),tmp(A))
c        write(*,*) "B="
c        call prt_cmat(ne(ms),tmp(B))
c        write(*,*) "C="
c        call prt_cmat(ne(ms),tmp(C))

        call cpsi_gen_Ba_Bs(ms,nb,nn,tmp(B),tmp(Bs),tmp(Ba))

        call cpsi_gen_UD(ms,nb,tmp(Bs),tmp(D))

        call cpsi_gen_X(ms,nb,nn,tmp(st1),tmp(st2),
     >                     tmp(A),tmp(Ba),tmp(C),
     >                     tmp(U),tmp(D),tmp(fnm),
     >                     fweight,
     >                     failed)

        if (failed) then
           ierr=10
           if (taskid.eq.MASTER) then
             write(luout,*) "Warning: Lagrange Mult generation failed."
             write(luout,*) " +Try using a smaller time step"
             write(luout,*) " +Gram-Schmidt performed,nb,ms=",nb,ms
           end if
           call Pneb_orthoCheckMake(.true.,ms,nb,npack1,psi2)
        else
           call Pneb_fwf_Multiply(ms,nb,
     >                       dtec,
     >                       psi1,npack1,
     >                       tmp(st1),
     >                       one,
     >                       psi2)
        end if
     
      return
      end


*     ***********************************
*     *                                 *
*     *        cpsi_gen_Ba_Bs           *
*     *                                 *
*     ***********************************
      subroutine cpsi_gen_Ba_Bs(ms,nb,nn,B,Bs,Ba)
      implicit none
      integer ms,nb,nn
      complex*16 B(*),Bs(*),Ba(*)

      call Pneb_w_copy_dagger(ms,nb,B,Ba)

      !call dcopy(nn,B,1,Bs,1)
      call Parallel_shared_vector_copy(.true.,2*nn,B,Bs)

      call DAXPY_OMP(2*nn,1.0d0,Ba,1,Bs,1)
      call DSCAL_OMP(2*nn,0.5d0,Bs,1)

      call DAXPY_OMP(2*nn,-1.0d0,B,1,Ba,1)
      call DSCAL_OMP(2*nn,-0.5d0,Ba,1)

      return
      end

*     ***********************************
*     *                                 *
*     *        cpsi_gen_UD              *
*     *                                 *
*     ***********************************
      subroutine cpsi_gen_UD(ms,nb,Bs,D)
      implicit none
      integer ms,nb
      complex*16 Bs(*)
      real*8     D(*)
      integer i

      !call DSYEV('V','U',n,Bs,n_max, D,Work,2*n_max*n_max,ierr)
      call Pneb_w_diag(ms,nb,D,Bs)

      return
      end



*     ***********************************
*     *                                 *
*     *        cpsi_gen_X               *
*     *                                 *
*     ***********************************
      subroutine cpsi_gen_X(ms,nb,nn,
     >                     X1,tmp,
     >                     A,Ba,C,
     >                     U,D,fnm,
     >                     fweight,
     >                     failed)

      implicit none
      integer ms,nb,nn
      complex*16 X1(*)
      complex*16 tmp(*)
      complex*16 A(*)
      complex*16 Ba(*)
      complex*16 C(*)
      complex*16 U(*)

      real*8 D(*)
      complex*16 fnm(*)
      real*8 fweight(*)
      logical failed

      !**** local variables ****
      integer itrlmd
      real*8  convg
      parameter (itrlmd=120, convg=1.0d-15)

      integer it,i,ne(2)
      real*8  adiff

*     **** external functions ****
      real*8   Pneb_w_dmax
      external Pneb_w_dmax

      call Pneb_ne(ne)

      !**** A = I-A ***
       call DSCAL_OMP(2*nn,(-1.0d0),A,1)
       call Pneb_w_eye(ms,nb,fnm,1.0d0)
       call DAXPY_OMP(2*nn,1.0d0,fnm,1,A,1)

      !*** fnm = I-A ****
      call Parallel_shared_vector_copy(.true.,2*nn,A,fnm)

      !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
      call cpsi_fnm_to_X(ms,nb,fnm,U,D,fweight,tmp)
      call Parallel_shared_vector_copy(.true.,2*nn,fnm,X1)

c       write(*,*) "X1="
c       call prt_cmat(ne(ms),X1)
c       call Pneb_www_Multiply2(ms,nb,dcmplx(1.0d0,0.0d0),X1,X1,
c     >                               dcmplx(0.0d0,0.0d0),tmp)
c       write(*,*) "X1'*X1="
c       call prt_cmat(ne(ms),tmp)

      it     = 0
      failed = .true.
      do while (failed .and. (it.lt.itrlmd))
        it = it + 1

        !*** fnm = X'*C*X ***
        call Pneb_www_Multiply1(ms,nb,dcmplx(1.0d0,0.0d0),C,X1,
     >                                dcmplx(0.0d0,0.0d0),tmp)
        call Pneb_www_Multiply1(ms,nb,dcmplx(1.0d0,0.0d0),X1,tmp,
     >                                dcmplx(0.0d0,0.0d0),fnm)

        !*** fnm = Ba*X - X'*C*X ***
        call Pneb_www_Multiply1(ms,nb,dcmplx(1.0d0,0.0d0),Ba,X1,
     >                                dcmplx(-1.0d0,0.0d0),fnm)

        !*** fnm = Ba*X - X'*Ba - X*C*X ***
        call Pneb_www_Multiply1(ms,nb,dcmplx(-1.0d0,0.0d0),X1,Ba,
     >                                dcmplx(1.0d0,0.0d0),fnm)

        !*** fnm = I-A + Ba*X - X'*Ba - X'*C*X ***
        call DAXPY_OMP(2*nn,1.0d0,A,1,fnm,1)


        !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
        call cpsi_fnm_to_X(ms,nb,fnm,U,D,fweight,tmp)

c       write(*,*) "X1=, it=",it
c       call prt_cmat(ne(ms),fnm)
c       call Pneb_www_Multiply2(ms,nb,dcmplx(1.0d0,0.0d0),fnm,fnm,
c     >                               dcmplx(0.0d0,0.0d0),tmp)
c       write(*,*) "X1'*X1=, it=",it
c       call prt_cmat(ne(ms),tmp)

        !call DMSUB(n_max,n,X1,fnm,tmp)
        !adiff = tmp(idamax(n_max*n,tmp,1))
        !call dcopy(n_max*n,fnm,1,X1,1)
        call Parallel_shared_vector_copy(.true.,2*nn,X1,tmp)
        call DAXPY_OMP(2*nn,-1.0d0,fnm,1,tmp,1)
        adiff = Pneb_w_dmax(ms,nb,tmp)
        call Parallel_shared_vector_copy(.true.,2*nn,fnm,X1)
        !write(*,*) "it,adiff=",it,adiff

        if (adiff.lt.convg) failed = .false.
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          cpsi_fnm_to_X          *
*     *                                 *
*     ***********************************
      subroutine cpsi_fnm_to_X(ms,nb,fnm,U,D,fweight,tmp)
      implicit none
      integer ms,nb
      complex*16 fnm(*)
      complex*16 U(*)
      real*8 D(*)
      real*8 fweight(*)
      complex*16 tmp(*)

      integer i


      !**** fnm = Ut*fnm*U ***
c      call DGEMM('N','N',n,n,n,1.0d0,
c     >           fnm,n_max,
c     >           U,n_max,
c     >           0.0d0,
c     >           tmp,n_max)
c      call DGEMM('T','N',n,n,n,1.0d0,
c     >           U,n_max,
c     >           tmp,n_max,
c     >           0.0d0,
c     >           fnm,n_max)
      call Pneb_www_Multiply1(ms,nb,dcmplx(1.0d0,0.0d0),fnm,U,
     >                              dcmplx(0.0d0,0.0d0),tmp)
      call Pneb_www_Multiply2(ms,nb,dcmplx(1.0d0,0.0d0),U,tmp,
     >                              dcmplx(0.0d0,0.0d0),fnm)

      !**** fnm = (Ut*fnm*U)_nm/(d_n+d_m) ***
c      do j=1,n
c      do i=1,n
c        d2 = D(i,1)+D(j,1)
c        fnm(i,j) = (fnm(i,j)/d2)
c      end do
c      end do
       call Pneb_w_HmldivideDplusD(ms,nb,fnm,D)

      !**** fnm = X = U*{(Ut*fnm*U)_nm/(d_n+d_m)}*Ut ***
c      call DGEMM('N','N',n,n,n,1.0d0,
c     >           U,n_max,
c     >           fnm,n_max,
c     >           0.0d0,
c     >           tmp,n_max)
c      call DGEMM('N','T',n,n,n,1.0d0,
c     >           tmp,n_max,
c     >           U,n_max,
c     >           0.0d0,
c     >           fnm,n_max)
      call Pneb_www_Multiply1(ms,nb,dcmplx(1.0d0,0.0d0),U,fnm,
     >                              dcmplx(0.0d0,0.0d0),tmp)
      call Pneb_www_Multiply3(ms,nb,dcmplx(1.0d0,0.0d0),tmp,U,
     >                              dcmplx(0.0d0,0.0d0),fnm)


c      do j=1,n
c      do i=1,n
c        fnm(i,j) = fnm(i,j)*(2.0d0*fweight(i)/(fweight(i)+fweight(j)))
c      end do
c      end do
      call Pneb_w_Hmlfweightscale(ms,nb,fnm,fweight)

      return
      end

      subroutine prt_cmat(n,C)
      implicit none
      integer n
      complex*16 C(n,n)

      integer i,j

      do i=1,n
         write(*,'(20("(",F8.4,E10.2,")"))') (C(i,j), j=1,n)
      end do

      return
      end

