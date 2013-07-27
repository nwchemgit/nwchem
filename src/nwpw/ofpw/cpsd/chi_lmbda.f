*     ***************************
*     *				*
*     *		chi_lmbda	*
*     *				*
*     ***************************

      subroutine chi_lmbda(ispin,npack1,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,ierr)

      implicit none
      integer ispin,npack1
      complex*16 psi1(npack1,ispin)
      complex*16 psi2(npack1,ispin)
      real*8     dte
      real*8     lmbda(*)
      integer	 ierr

*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd,idamax
      real*8  convg
      parameter (itrlmd=20, convg=1.0d-15)

*     **** local variables ****
      integer ms,it
      real*8  s11,s12,s21,s22,st1,st2,sa1,sa0
      real*8  sl(2)
      real*8  sum,adiff,alpha

      s11  = 0.0d0
      s12  = 0.0d0
      s21  = 0.0d0
      s22  = 0.0d0
      sa0  = 0.0d0
      sa1  = 0.0d0
      st1  = 0.0d0
      st2  = 0.0d0

      sl(1)  = 0.0d0
      sl(2)  = 0.0d0
      call dcopy(ispin,0.0d0,0,lmbda,1)

      
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      DO 640 ms=1,ispin

*       ***** compute the overlap matrices ****
        call Pack_cc_dot(1,psi2(1,ms),psi2(1,ms),s22)
        call Pack_cc_dot(1,psi2(1,ms),psi1(1,ms),s21)
        call Pack_cc_dot(1,psi1(1,ms),psi1(1,ms),s11)


*       ***** scale the overlap matrices ****
        s22=(1.0d0-s22)*0.5d0/dte
        s21=(1.0d0-s21)*0.5d0
        s12=s21
        s11= (-s11)*0.5d0*dte

        sa0 = s22

        do it=1,itrlmd
          sa1 = s22

          sa1 = sa1 + s21*sa0 + sa0*s12

          st1 = s11*sa0
          sa1 = sa1 + sa0*st1

          st1 = sa1-sa0
          adiff = st1 
          if (adiff.lt.convg) go to 630
          sa0 = sa1

        end do
        ierr=10
        WRITE(6,*) 'ierr=10',adiff
C       return
  630   continue
         lmbda(ms) = sa1
  640 continue

*     **** correction due to the constraint ****
      do ms=1,ispin
         alpha = dte*lmbda(ms)
         call Pack_cc_daxpy(1,alpha,psi1(1,ms),psi2(1,ms))
      end do

      return
      end

c $Id$
