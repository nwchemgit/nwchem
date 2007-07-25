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


*     **** local variables ****
      integer i,j,ii,jj,ms,it,k
      integer index,indext,n,nn
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      integer sl(2)
      real*8  sum,adiff,alpha
*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd,idamax
      real*8  convg
      parameter (itrlmd=20, convg=1.0d-15)



      call current_second(tim1)

      n    = 1
      nn   = n**2
      
      s11  = 0*nn + 1
      s12  = 1*nn + 1
      s21  = 2*nn + 1
      s22  = 3*nn + 1
      sa0  = 4*nn + 1
      sa1  = 5*nn + 1
      st1  = 6*nn + 1
      st2  = 7*nn + 1

      call dcopy(8*nn,0.0d0,0,tmp,1)

      sl(1)  = 0*nn + 1
      sl(2)  = 1*nn + 1
      call dcopy(2*nn,0.0d0,0,lmbda,1)

      
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      DO 640 ms=1,ispin
        IF(ne(ms).le.0) GO TO 640

*       ***** compute the overlap matrices ****
        call Pack_cc_dot(1,psi2(1,ms),psi2(1,ms),s22)
        call Pack_cc_dot(1,psi2(1,ms),psi1(1,ms),s21)
        call Pack_cc_dot(1,psi1(1,ms),psi2(1,ms),s12)
        call Pack_cc_dot(1,psi1(1,ms),psi1(1,ms),s11)

*       ***** scale the overlap matrices ****
        s22=(1.0d0-s22)*0.5d0/dte
        s21=(1.0d0-s21)*0.5d0
        s12=(1.0d0-s12)*0.5d0
        s11= (-s11)*0.5d0*dte

        call dcopy(nn,tmp(s22),1,tmp(sa0),1)
        sa0 = s22

        do it=1,itrlmd
          sa1 = s22
          st1 = s21*sa0
          st2 = sa0*s12

          st1 = st1+st2
          sa1 = st1+sa1

          st1 = s11*sa0
          st2 = sa0*st1

          sa1 = st2+sa1
          st1 = sa1-sa0
          adiff = st1 
          if (adiff.lt.convg) go to 630
          sa0 = sa1

c          CALL dcopy(nn,tmp(s22),1,tmp(sa1),1)
c          CALL DMMUL(n,ne(MS), tmp(s21), tmp(sa0), tmp(st1))
c          CALL DMMUL(n,ne(MS), tmp(sa0), tmp(s12), tmp(st2))
c
c          CALL DMADD(n,ne(MS), tmp(st1), tmp(st2), tmp(st1))
c          CALL DMADD(n,ne(MS), tmp(st1), tmp(sa1), tmp(sa1))
c
c          CALL DMMUL(n,ne(MS), tmp(s11), tmp(sa0), tmp(st1))
c          CALL DMMUL(n,ne(MS), tmp(sa0), tmp(st1), tmp(st2))
c
c          CALL DMADD(n,ne(MS), tmp(st2), tmp(sa1), tmp(sa1))
c          CALL DMSUB(n,ne(MS), tmp(sa1), tmp(sa0), tmp(st1))
c          adiff=tmp(st1 - 1 + (idamax(n*ne(ms),tmp(st1),1)))

c          if(adiff.lt.convg) GO TO 630
c          call dcopy(n*ne(ms),tmp(sa1),1,tmp(sa0),1)
        end do
        ierr=10
        WRITE(6,*) 'ierr=10',adiff
C       return
  630   continue
c        call dcopy(n*ne(ms),tmp(sa1),1,lmbda(sl(ms)),1)
         lmbda(ms) = sa1
  640 continue

*     **** correction due to the constraint ****
      do ms=1,ispin
         alpha = dte*lmbda(ms)
         call Pack_cc_daxpy(1,alpha,psi1(1,ms),psi2(1,ms))
      end do


      return
      end

