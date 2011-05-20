*
* $Id$
*


*     ***********************************
*     *                                 *
*     *          paw_psi_MakeOrtho      *
*     *                                 *
*     ***********************************

*   This routine orthonormalizes the orbitals using a modified
* Gram-Schmidt algorithm.
*
      subroutine paw_psi_MakeOrtho(npack,ne,psi)
      implicit none
      integer npack,ne
      double complex psi(npack,ne)

*     **** local variables ****
      integer j,k
      real*8  w

      do k=1,ne
         call paw_overlap_matrix_gen(1,1,psi(1,k),psi(1,k),w)
         w = 1.0d0/dsqrt(w)
c         call Pack_c_SMul(1,w,psi(1,k),psi(1,k))
         call Pack_c_SMul1(1,w,psi(1,k))

         do j=k+1,ne
            call paw_overlap_matrix_gen(1,1,psi(1,k),psi(1,j),w)
            w = -w
            call Pack_cc_daxpy(1,w,psi(1,k),psi(1,j))
         end do
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          paw_psi_CheckOrtho     *
*     *                                 *
*     ***********************************

*   This routine return true if the orbitals are
* orthonormal.

      real*8 function paw_psi_CheckOrtho(npack,ne,psi)
      implicit none
      integer npack,ne
      double complex psi(npack,ne)


*     **** local variables ****
      integer j,k
      real*8  w,error

      error = 0.0d0
      do k=1,ne
        call paw_overlap_matrix_gen(1,1,psi(1,k),psi(1,k),w)
        error = error + dabs(1.0d0-w)

        do j=k+1,ne
           call paw_overlap_matrix_gen(1,1,psi(1,j),psi(1,k),w)
           error = error + dabs(w)
        end do
      end do

      paw_psi_CheckOrtho = error
      return
      end


*     ***********************************
*     *                                 *
*     *          paw_psi_CheckOrtho2    *
*     *                                 *
*     ***********************************

*   This routine return true if the orbitals are
* orthonormal.

      subroutine paw_psi_CheckOrtho2(npack,ne,psi)
      implicit none
      integer npack,ne
      double complex psi(npack,ne)


*     **** local variables ****
      integer j,k
      real*8  w,error

      
      write(*,*)
      do k=1,ne
        call paw_overlap_matrix_gen(1,1,psi(1,k),psi(1,k),w)
        write(*,*) "CheckOrtho2:",k,k,w

        do j=k+1,ne
           call paw_overlap_matrix_gen(1,1,psi(1,j),psi(1,k),w)
        write(*,*) "CheckOrtho2:",j,k,w
        end do
      end do
      write(*,*)

      return
      end



*     ***********************************
*     *                                 *
*     *          paw_psi_lmbda          *
*     *                                 *
*     ***********************************

      subroutine paw_psi_lmbda(ispin,ne,nemax,npack1,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr)

      implicit none
      integer ispin,ne(2),nemax,npack1
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dte
      real*8     lmbda(*)
      real*8     tmp(*)
      integer    ierr

      integer MASTER
      parameter (MASTER=0)

*     **** local variables ****
      logical failed
      integer taskid
      integer n1(2),n2(2)
      integer i,j,ii,jj,ms
      integer n,nn,index
      integer st1,st2
      integer A,B,C,U,D,Ba,Bs,fnm
      integer sl(2)
      real*8  alpha



      call nwpw_timing_start(3)

      call Parallel_taskid(taskid)

      n    = ne(1)
      nn   = n**2
      
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

      call dcopy(8*nn,0.0d0,0,tmp,1)

      sl(1)  = 0*nn + 1
      sl(2)  = 1*nn + 1
      call dcopy(2*nn,0.0d0,0,lmbda,1)

      n1(1)=1
      n2(1)=ne(1)
      n1(2)=ne(1)+1
      n2(2)=ne(1)+ne(2)
      

      do ms=1,ispin
        IF(ne(ms).le.0) go to 640


*       ***** compute the overlap matrices ****
        call paw_overlap_sym_matrix_gen(n,ne(ms),
     >                          psi2(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(A))
        call paw_overlap_matrix_gen(n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(B))
        call paw_overlap_sym_matrix_gen(n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi1(1,n1(ms)),
     >                          tmp(C))


        call paw_psi_gen_Ba_Bs(n,ne(ms),tmp(B),tmp(Bs),tmp(Ba))
        call paw_psi_gen_UD(n,ne(ms),tmp(Bs),tmp(D),lmbda)


        call paw_psi_gen_X(n,ne(ms),tmp(st1),tmp(st2),
     >                     tmp(A),tmp(Ba),tmp(C),
     >                     tmp(U),tmp(D),tmp(fnm),
     >                     failed)

        if (failed) then 
          if (taskid.eq.MASTER) then
            write(6,*) 
     >     'Warning: Lagrange Multiplier generation failed.'
            write(6,*) '        +Try using a smaller time step'
            write(6,*) '        +Gram-Schmidt being performed, spin:',ms
          end if
          call paw_psi_MakeOrtho(npack1,ne(ms),psi2(1,n1(ms)))
        else
          call dcopy(n*ne(ms),tmp(st1),1,lmbda(sl(ms)),1)
          call dscal(n*ne(ms),(1.0d0/dte),lmbda(sl(ms)),1)


*         ****  correction due to the constraint ****
          call dgemm('N','N',2*npack1,ne(ms),ne(ms),
     >              (1.0d0),
     >              psi1(1,n1(ms)),2*npack1,
     >              tmp(st1),n,
     >              (1.0d0),
     >              psi2(1,n1(ms)),2*npack1)

        end if
  640   continue
      end do !*ms*
      call nwpw_timing_end(3)

      return
      end



*     ***********************************
*     *                                 *
*     *        paw_psi_gen_Ba_Bs        *
*     *                                 *
*     ***********************************
      subroutine paw_psi_gen_Ba_Bs(n_max,n,B,Bs,Ba)
      implicit none
      integer n_max,n
      real*8 B(n_max,n)
      real*8 Bs(n_max,n)
      real*8 Ba(n_max,n)

      !*** local variables ***
      integer i,j

      do i=1,n
      do j=1,n
         Bs(i,j) = 0.5d0*(B(i,j)+B(j,i))
         Ba(i,j) = 0.5d0*(B(i,j)-B(j,i))
      end do
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *        paw_psi_gen_UD           *
*     *                                 *
*     ***********************************
      subroutine paw_psi_gen_UD(n_max,n,Bs,D,work)
      implicit none
      integer n_max,n
      real*8 Bs(n_max,n)
      real*8 D(n_max,n)
      real*8 Work(n_max,n)

      !*** local variables ***
      integer ierr

      !call eigen(n_max,n,Bs,D,D(1,2))
      call dsyev('V','U',n,Bs,n_max, D,Work,2*n_max*n_max,ierr)
      return
      end




*     ***********************************
*     *                                 *
*     *        paw_psi_gen_X            *
*     *                                 *
*     ***********************************
      subroutine paw_psi_gen_X(n_max,n,
     >                     X1,tmp,
     >                     A,Ba,C,
     >                     U,D,fnm,
     >                     failed)
     > 
      implicit none
      integer n_max,n
      real*8 X1(n_max,n)
      real*8 tmp(*)
      real*8 A(n_max,n)
      real*8 Ba(n_max,n)
      real*8 C(n_max,n)
      real*8 U(n_max,n)
      real*8 D(n_max,n)
      real*8 fnm(n_max,n)
      logical failed

      !**** local variables ****
      integer itrlmd
      real*8  convg
      parameter (itrlmd=40, convg=1.0d-15)

      integer i,it
      real*8  adiff

      !**** external functions ****
      integer  idamax
      external idamax


      !**** A = I-A ***
      call dscal(n_max*n,(-1.0d0),A,1)
      do i=1,n
         A(i,i) = A(i,i) + 1.0d0
      end do

      !*** fnm = I-A ****
      call dcopy(n_max*n,A,1,fnm,1)

      !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
      call paw_psi_fnm_to_X(n_max,n,fnm,U,D,tmp)
      call dcopy(n_max*n,fnm,1,X1,1)


      it     = 0
      failed = .true.
      do while (failed .and. (it.lt.itrlmd))
        it = it + 1

        !*** fnm = X*C*X ***
        call DMMUL(n_max,n,C,X1,tmp)
        call DMMUL(n_max,n,X1,tmp,fnm)


        !*** fnm = Ba*X - X*C*X ***
        call DMMUL(n_max,n,Ba,X1,tmp)
        call DMSUB(n_max,n,tmp,fnm,fnm)


        !*** fnm = Ba*X - X*Ba - X*C*X ***
        call DMMUL(n_max,n,X1,Ba,tmp)
        call DMSUB(n_max,n,fnm,tmp,fnm)

        !*** fnm = I-A + Ba*X - X*Ba - X*C*X ***
        call DMADD(n_max,n,fnm,A,fnm)


        !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
        call paw_psi_fnm_to_X(n_max,n,fnm,U,D,tmp)

        call DMSUB(n_max,n,fnm,X1,tmp)
        adiff = tmp(idamax(n_max*n,tmp,1))
        call dcopy(n_max*n,fnm,1,X1,1)

        adiff = dabs(adiff)

        if (adiff.lt.convg) failed = .false.
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *        paw_psi_fnm_to_X         *
*     *                                 *
*     ***********************************
      subroutine paw_psi_fnm_to_X(n_max,n,
     >                            fnm,U,D,tmp)
      implicit none
      integer n_max,n
      real*8 fnm(n_max,n)
      real*8 U(n_max,n)
      real*8 D(n_max,n)
      real*8 tmp(n_max,n)

      !**** local variables ****
      integer i,j
      real*8  d2


      !**** fnm = Ut*fnm*U ***
      call dgemm('N','N',n,n,n,1.0d0,
     >           fnm,n_max,
     >           U,n_max,
     >           0.0d0,
     >           tmp,n_max)
      call dgemm('T','N',n,n,n,1.0d0,
     >           U,n_max,
     >           tmp,n_max,
     >           0.0d0,
     >           fnm,n_max)


      !**** fnm = (Ut*fnm*U)_nm/(d_n+d_m) ***
      do j=1,n
      do i=1,n
        d2 = D(i,1)+D(j,1)
        fnm(i,j) = fnm(i,j)/d2
      end do
      end do

      !**** fnm = X = U*{(Ut*fnm*U)_nm/(d_n+d_m)}*Ut ***
      call dgemm('N','N',n,n,n,1.0d0,
     >           U,n_max,
     >           fnm,n_max,
     >           0.0d0,
     >           tmp,n_max)
      call dgemm('N','T',n,n,n,1.0d0,
     >           tmp,n_max,
     >           U,n_max,
     >           0.0d0,
     >           fnm,n_max)

      return
      end 
