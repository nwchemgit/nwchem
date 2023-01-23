c $Id: psi_lmbda_paw.f 21177 2011-10-10 17:09:43Z bylaska $

*     ***********************************
*     *                                 *
*     *          psi_lmbda_paw       	*
*     *                                 *
*     ***********************************

      subroutine psi_lmbda_paw(ispin,ne,nemax,npack1,
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
      integer ms,nn
      integer st1,st2
      integer A,B,C,U,D,Ba,Bs,fnm
      integer i,j

      call nwpw_timing_start(3)

      ierr = 0
      call Parallel_taskid(taskid)

      call Dneall_m_size(1,nn)
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

      !call dcopy(8*nn,0.0d0,0,tmp,1)
      call Parallel_shared_vector_zero(.true.,8*nn,tmp)

      do ms=1,ispin
        IF(ne(ms).le.0) go to 640

*       ***** compute the overlap matrices ****
        call Dneall_ffm_sym_Multiply(ms,psi2,psi2,npack1,tmp(A))
        call Dneall_ffm_Multiply(ms,psi1,psi2,npack1,tmp(B))
        call Dneall_ffm_sym_Multiply(ms,psi1,psi1,npack1,tmp(C))

*       *** add the extra paw overlap part ***
        call psp_add_paw_extra_overlap1(ms,psi2,tmp(A))
        call psp_add_paw_extra_overlap2(ms,psi1,psi2,tmp(B))
        call psp_add_paw_extra_overlap1(ms,psi1,tmp(C))

        call psi_gen_Ba_Bs(ms,nn,tmp(B),tmp(Bs),tmp(Ba))

        call psi_gen_UD(ms,tmp(Bs),tmp(D))

        call psi_gen_X_paw(ms,nn,tmp(st1),tmp(st2),
     >                     tmp(A),tmp(Ba),tmp(C),
     >                     tmp(U),tmp(D),tmp(fnm),
     >                     failed)

        if (failed) then
          if (taskid.eq.MASTER) then
            write(*,*)
     >     'Warning: Lagrange Multiplier generation failed.'
            write(*,*) '        +Try using a smaller time step'
            write(*,*) '        +Gram-Schmidt being performed, spin:',ms
          end if
c          call Dneall_f_GramSchmidt(ms,psi2,npack1)
           call Dneall_f_Sortho(ms,psi2,psi1,npack1)

          ierr = 1
        else
          call Dneall_fmf_Multiply(ms,
     >                          psi1,npack1,
     >                          tmp(st1), 1.0d0,
     >                          psi2,1.0d0)
          call dscal_omp(nn,(1.0d0/dte),tmp(st1),1)
          call Dneall_mm_Expand(ms,tmp(st1),lmbda)
        end if
   
  640   continue
      end do !*ms*

      call nwpw_timing_end(3)

      return
      end



*     ***********************************
*     *                                 *
*     *        psi_gen_X_paw            *
*     *                                 *
*     ***********************************
      subroutine psi_gen_X_paw(ms,nn,
     >                     X1,tmp,
     >                     A,Ba,C,
     >                     U,D,fnm,
     >                     failed)
     
      implicit none
      integer ms,nn
      real*8 X1(*)
      real*8 tmp(*)
      real*8 A(*)
      real*8 Ba(*)
      real*8 C(*)
      real*8 U(*)
      real*8 D(*)
      real*8 fnm(*)
      logical failed

      !**** local variables ****
      integer itrlmd
      real*8  convg
      parameter (itrlmd=120, convg=1.0d-15)

      integer it
      real*8  adiff
      common /psi_gen_X_paw_tmp/adiff

*     **** external functions ****
      real*8   Dneall_m_dmax
      external Dneall_m_dmax

      !**** A = I-A ***
       call dscal_omp(nn,(-1.0d0),A,1)
       call Dneall_m_eye(ms,fnm,1.0d0)
       call daxpy_omp(nn,1.0d0,fnm,1,A,1)

      !*** fnm = I-A ****
      !call dcopy(nn,A,1,fnm,1)
      call Parallel_shared_vector_copy(.true.,nn,A,fnm)

      !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
      call psi_fnm_to_X_paw(ms,fnm,U,D,tmp)
      !call dcopy(nn,fnm,1,X1,1)
      call Parallel_shared_vector_copy(.true.,nn,fnm,X1)


      it     = 0
      failed = .true.
      do while (failed .and. (it.lt.itrlmd))
        it = it + 1

        !*** fnm = X*C*X ***
        call Dneall_mmm_Multiply(ms,C, X1,  1.0d0,tmp,0.0d0)
        call Dneall_mmm_Multiply(ms,X1,tmp, 1.0d0,fnm,0.0d0)


        !*** fnm = Ba*X - X*C*X ***
        call Dneall_mmm_Multiply(ms,Ba,X1,1.0d0,fnm,-1.0d0)


        !*** fnm = Ba*X - X*Ba - X*C*X ***
        call Dneall_mmm_Multiply(ms,X1,Ba,-1.0d0,fnm,1.0d0)

        !*** fnm = I-A + Ba*X - X*Ba - X*C*X ***
        call daxpy_omp(nn,1.0d0,A,1,fnm,1)


        !*** solve U*D*Ut*X + X*U*D*Ut = fnm for X ***
        call psi_fnm_to_X_paw(ms,fnm,U,D,tmp)

        !call DMSUB(n_max,n,X1,fnm,tmp)
        !adiff = tmp(idamax(n_max*n,tmp,1))
        !call dcopy(n_max*n,fnm,1,X1,1)
        !call dcopy(nn,X1,1,tmp,1)
        call Parallel_shared_vector_copy(.true.,nn,X1,tmp)
        call daxpy_omp(nn,-1.0d0,fnm,1,tmp,1)
!$OMP MASTER
        adiff = Dneall_m_dmax(ms,tmp)
!$OMP END MASTER
        !call dcopy(nn,fnm,1,X1,1)
        call Parallel_shared_vector_copy(.true.,nn,fnm,X1)

        if (adiff.lt.convg) failed = .false.
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *        psi_fnm_to_X_paw         *
*     *                                 *
*     ***********************************
      subroutine psi_fnm_to_X_paw(ms,fnm,U,D,tmp)
      implicit none
      integer ms
      real*8 fnm(*)
      real*8 U(*)
      real*8 D(*)
      real*8 tmp(*)


      !**** fnm = Ut*fnm*U ***
      call Dneall_mmm_Multiply(ms,fnm,U,1.0d0,tmp,0.0d0)
      call Dneall_mmm_Multiply2(ms,U,tmp,fnm)

      !**** fnm = (Ut*fnm*U)_nm/(d_n+d_m) ***
       call Dneall_m_HmldivideDplusD(ms,fnm,D)

      !**** fnm = X = U*{(Ut*fnm*U)_nm/(d_n+d_m)}*Ut ***
      call Dneall_mmm_Multiply(ms,U,fnm,1.0d0,tmp,0.0d0)
      call Dneall_mmm_Multiply3(ms,tmp,U,fnm)

      return
      end


