
*     ***********************************
*     *                                 *
*     *          psi_lmbda2          	*
*     *                                 *
*     ***********************************

      subroutine psi_lmbda2(ispin,ne,nemax,npack1,
     >                     psi1,psi2,
     >                     dte,fweight,
     >                     lmbda,tmp,ierr)
      implicit none
      integer ispin,ne(2),nemax,npack1
      complex*16 psi1(npack1,nemax)
      complex*16 psi2(npack1,nemax)
      real*8     dte
      real*8     lmbda(*)
      real*8     fweight(*)
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
        call Pack_ccm_sym_dot2(1,n,ne(ms),
     >                          psi2(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(A))
        call Pack_ccmn_dot(1,n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(B))
        call Pack_ccm_sym_dot2(1,n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi1(1,n1(ms)),
     >                          tmp(C))


        call psi_gen_Ba_Bs(n,ne(ms),tmp(B),tmp(Bs),tmp(Ba))
        call psi_gen_UD(n,ne(ms),tmp(Bs),tmp(D),lmbda)


        call psi_gen_X(n,ne(ms),tmp(st1),tmp(st2),
     >                     tmp(A),tmp(Ba),tmp(C),
     >                     tmp(U),tmp(D),tmp(fnm),
     >                     fweight(1+(ms-1)*ne(1)),
     >                     failed)

        if (failed) then
          if (taskid.eq.MASTER) then
            write(6,*)
     >     'Warning: Lagrange Multiplier generation failed.'
            write(6,*) '        +Try using a smaller time step'
            write(6,*) '        +Gram-Schmidt being performed, spin:',ms
          end if
          call Dneall_f_ortho(ms,psi2,npack1)
c          call Grsm_g_MakeOrtho(npack1,ne(ms),psi2(1,n1(ms)))
        else
          call dcopy(n*ne(ms),tmp(st1),1,lmbda(sl(ms)),1)
          call dscal(n*ne(ms),(1.0d0/dte),lmbda(sl(ms)),1)
c         do j=1,ne(ms)
c         do i=1,ne(ms)
c          lmbda(sl(ms)+(i-1)+(j-1)*n)
c    >      =lmbda(sl(ms)+(i-1)+(j-1)*n)*fweight(j)
c         end do
c         end do

*         ****  correction due to the constraint ****
          call DGEMM('N','N',2*npack1,ne(ms),ne(ms),
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
*     *        psi_gen_Ba_Bs   	     	*
*     *                                 *
*     ***********************************
      subroutine psi_gen_Ba_Bs(n_max,n,B,Bs,Ba)
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
*     *        psi_gen_UD           	*
*     *                                 *
*     ***********************************
      subroutine psi_gen_UD(n_max,n,Bs,D,work)
      implicit none
      integer n_max,n
      real*8 Bs(n_max,n)
      real*8 D(n_max,n)
      real*8 Work(n_max,n)

      !*** local variables ***
      integer ierr

      !call eigen(n_max,n,Bs,D,D(1,2))
      call DSYEV('V','U',n,Bs,n_max, D,Work,2*n_max*n_max,ierr)
      return
      end




*     ***********************************
*     *                                 *
*     *        psi_gen_X            	*
*     *                                 *
*     ***********************************
      subroutine psi_gen_X(n_max,n,
     >                     X1,tmp,
     >                     A,Ba,C,
     >                     U,D,fnm,
     >                     fweight,
     >                     failed)
     
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
      real*8 fweight(n)
      logical failed

      !**** local variables ****
      integer itrlmd
      real*8  convg
      parameter (itrlmd=20, convg=1.0d-15)

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
      call psi_fnm_to_X(n_max,n,fnm,U,D,fweight,tmp)
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
        call psi_fnm_to_X(n_max,n,fnm,U,D,fweight,tmp)

        call DMSUB(n_max,n,X1,fnm,tmp)
        adiff = tmp(idamax(n_max*n,tmp,1))
        call dcopy(n_max*n,fnm,1,X1,1)


        if (adiff.lt.convg) failed = .false.
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *        psi_fnm_to_X         *
*     *                                 *
*     ***********************************
      subroutine psi_fnm_to_X(n_max,n,fnm,U,D,fweight,tmp)
      implicit none
      integer n_max,n
      real*8 fnm(n_max,n)
      real*8 U(n_max,n)
      real*8 D(n_max,n)
      real*8 fweight(n)
      real*8 tmp(n_max,n)

      !**** local variables ****
      integer i,j
      real*8  d2


      !**** fnm = Ut*fnm*U ***
      call DGEMM('N','N',n,n,n,1.0d0,
     >           fnm,n_max,
     >           U,n_max,
     >           0.0d0,
     >           tmp,n_max)
      call DGEMM('T','N',n,n,n,1.0d0,
     >           U,n_max,
     >           tmp,n_max,
     >           0.0d0,
     >           fnm,n_max)


      !**** fnm = (Ut*fnm*U)_nm/(d_n+d_m) ***
      do j=1,n
      do i=1,n
        d2 = D(i,1)+D(j,1)
        fnm(i,j) = (fnm(i,j)/d2)
      end do
      end do

      !**** fnm = X = U*{(Ut*fnm*U)_nm/(d_n+d_m)}*Ut ***
      call DGEMM('N','N',n,n,n,1.0d0,
     >           U,n_max,
     >           fnm,n_max,
     >           0.0d0,
     >           tmp,n_max)
      call DGEMM('N','T',n,n,n,1.0d0,
     >           tmp,n_max,
     >           U,n_max,
     >           0.0d0,
     >           fnm,n_max)

      do j=1,n
      do i=1,n
        fnm(i,j) = fnm(i,j)*(2.0d0*fweight(i)/(fweight(i)+fweight(j)))
      end do
      end do

      return
      end





