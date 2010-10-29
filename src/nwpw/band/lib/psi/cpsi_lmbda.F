*
* $Id$
*
*     ****************************************************
*     *                                                  *
*     *                 cpsi_lmbda                       *
*     *                                                  *
*     ****************************************************
      subroutine cpsi_lmbda(ms,nb,npack1,
     >                     psi1,psi2,
     >                     dte,
     >                     tmp,ierr)
      implicit none
      integer ms,nb,npack1
      complex*16 psi1(*)
      complex*16 psi2(*)
      real*8     dte
      complex*16 tmp(*)
      integer    ierr


*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
  
      logical value,gramschmidt
      integer it
      integer nn
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      real*8 adiff
      complex*16 zero,one,dtec
c      complex*16 zadiff
*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd
      real*8  convg
      parameter (itrlmd=20, convg=1.d-10)

*     **** external functions ****
      integer  Pneb_w_size
      real*8   Pneb_w_dmax
      external Pneb_w_size
      external Pneb_w_dmax

      call Parallel_taskid(taskid)

      gramschmidt = .false.
      ierr = 0
      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
      dtec = dcmplx(dte,0.0d0)
*     **** set the matrices dimensions ****
      nn   = Pneb_w_size(0,nb)
      
      s11  = 0*nn + 1
      s12  = 1*nn + 1
      s21  = 2*nn + 1
      s22  = 3*nn + 1
      sa0  = 4*nn + 1
      sa1  = 5*nn + 1
      st1  = 6*nn + 1
      st2  = 7*nn + 1

      call zcopy(8*nn,zero,0,tmp,1)
     
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::

        call Pneb_ffw_hermit_Multiply(ms,nb,psi2,psi2,npack1,tmp(s22))
        call Pneb_ffw_hermit_Multiply(ms,nb,psi2,psi1,npack1,tmp(s21))
        call Pneb_ffw_hermit_Multiply(ms,nb,psi2,psi1,npack1,tmp(s12))
        call Pneb_ffw_hermit_Multiply(ms,nb,psi1,psi1,npack1,tmp(s11))

        call Pneb_w_scale_s22(ms,nb,dte,tmp(s22))
        call Pneb_w_scale_s21(ms,nb,dte,tmp(s21))
        call Pneb_w_scale_s21(ms,nb,dte,tmp(s12))
        call Pneb_w_scale_s11(ms,nb,dte,tmp(s11))


        call zcopy(nn,tmp(s22),1,tmp(sa0),1)
        do it=1,itrlmd
          call zcopy(nn,tmp(s22),1,tmp(sa1),1)

          call Pneb_www_Multiply1(ms,nb,
     >                            one,
     >                            tmp(s21),tmp(sa0),
     >                            one,
     >                            tmp(sa1))

          call Pneb_www_Multiply2(ms,nb,
     .                            one,
     >                            tmp(sa0),tmp(s12),
     >                            one,
     >                            tmp(sa1))

          call Pneb_www_Multiply3(ms,nb,
     >                            one,
     >                            tmp(s11),tmp(sa0),
     >                            zero,
     >                            tmp(st1))
          call Pneb_www_Multiply1(ms,nb,
     >                            one,
     >                            tmp(sa0),tmp(st1),
     >                            one,
     >                            tmp(sa1))
          call zcopy(nn,tmp(sa1),1,tmp(st1),1)
          call zaxpy(nn,-one,tmp(sa0),1,tmp(st1),1)

          adiff = Pneb_w_dmax(ms,nb,tmp(st1))
          if(adiff.lt.convg) GO TO 630
          call zcopy(nn,tmp(sa1),1,tmp(sa0),1)
        end do

        ierr=10
        if (taskid.eq.MASTER) then
          WRITE(6,*)
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          WRITE(6,*) '        +Try using a smaller time step'
          WRITE(6,*) '        +Gram-Schmidt being performed, spin:',ms
        end if
        call Pneb_orthoCheckMake(.true.,ms,nb,npack1,psi2)
        gramschmidt = .true.

  630   continue

         if (.not.gramschmidt)
     >     call Pneb_fwf_Multiply(ms,nb,
     >                       dtec,
     >                       psi1,npack1,
     >                       tmp(sa1),
     >                       one,
     >                       psi2)

      return
      end

