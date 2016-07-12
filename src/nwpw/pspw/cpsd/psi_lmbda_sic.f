*
* $Id$
*

      subroutine psi_lmbda_sic(ispin,ne,nemaxq,npack1,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr)

      implicit none
      integer ispin,ne(2),nemaxq,npack1
      complex*16 psi1(npack1,nemaxq)
      complex*16 psi2(npack1,nemaxq)
      real*8     dte
      real*8     lmbda(*)
      real*8     tmp(*)
      integer    ierr

*     **** parameters ****
      integer MASTER
      parameter (MASTER=0)

*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd
      real*8  convg
      parameter (itrlmd=150, convg=1.0d-12)

*     **** local variables ****
      logical notgram
      integer taskid
      integer ms,it,nn
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      real*8  adiff
c      real*8  tmp1(1000000)

*     **** external functions ****
      real*8   Dneall_m_dmax
      external Dneall_m_dmax


      call nwpw_timing_start(3)

      call Dneall_m_size(1,nn)
      s11  = 0*nn + 1
      s12  = 1*nn + 1
      s21  = 2*nn + 1
      s22  = 3*nn + 1
      sa0  = 4*nn + 1
      sa1  = 5*nn + 1
      st1  = 6*nn + 1
      st2  = 7*nn + 1
      call dcopy(8*nn,0.0d0,0,tmp,1)

*     ::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      DO 640 ms=1,ispin
        notgram = .true.
        IF(ne(ms).le.0) GO TO 640

        call Dneall_ffm_sym_Multiply(ms,psi2,psi2,npack1,tmp(s22))
        call Dneall_ffm_sym_Multiply(ms,psi2,psi1,npack1,tmp(s21))
        call Dneall_ffm_sym_Multiply(ms,psi1,psi2,npack1,tmp(s12))
        call Dneall_ffm_sym_Multiply(ms,psi1,psi1,npack1,tmp(s11))

****  Begin  ADDED by Kiril *****
        call Dneall_m_Kiril_BTransform(ms,tmp(12),tmp(s21))
c        do ii=1,ne(ms)
c           do jj=1,ne(ms)
c             ii1 = -1+ii+jj*(ii-1)
c             tmp1(ii1+1)=0.5d0*(tmp(s12+ii1)+tmp(s21+ii1))
c
c           enddo
c        enddo
c
c        do ii=1,ne(ms)
c           do jj=1,ii-1
c             ii1 = -1+ii+jj*(ii-1)
c             tmp(s12+ii1)=tmp1(ii1+1)
c             tmp(s21+ii1)=tmp1(ii1+1)
c           enddo
c        enddo
****  End    ADDED by Kiril *****

*       ***** scale the overlap matrices ****
        call Dneall_m_scale_s22(ms,dte,tmp(s22))
        call Dneall_m_scale_s21(ms,dte,tmp(s21))
        call Dneall_m_scale_s21(ms,dte,tmp(s12))
        call Dneall_m_scale_s11(ms,dte,tmp(s11))
          
        call dcopy(nn,tmp(s22),1,tmp(sa0),1)

        do it=1,itrlmd
          CALL dcopy(nn,tmp(s22),1,tmp(sa1),1)

          call Dneall_mmm_Multiply(ms,
     >                              tmp(s21),tmp(sa0),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(s12),1.0d0,
     >                              tmp(sa1),1.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(s11),tmp(sa0),1.0d0,
     >                              tmp(st1),0.0d0)
          call Dneall_mmm_Multiply(ms,
     >                              tmp(sa0),tmp(st1),1.0d0,
     >                              tmp(sa1),1.0d0)
          call dcopy(nn,tmp(sa1),1,tmp(st1),1)
          call daxpy(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)

          adiff = Dneall_m_dmax(ms,tmp(st1))
          if(adiff.lt.convg) GO TO 630
          call dcopy(nn,tmp(sa1),1,tmp(sa0),1)
        end do

        ierr=10
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
          write(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          write(6,*) '        +Try using a smaller time step'
          write(6,*) '        +Gram-Schmidt being performed, spin:',ms
        end if
c        call Dneall_f_ortho(ms,psi2,npack1)
        call Dneall_f_GramSchmidt(ms,psi2,npack1)
        notgram = .false.

  630   continue

*       :::::::::::::::::  correction due to the constraint  :::::::::::::::::
        if (notgram)
     >     call Dneall_fmf_Multiply(ms,
     >                             psi1,npack1,
     >                             tmp(sa1), dte,
     >                             psi2,1.0d0)
        call Dneall_mm_Expand(ms,tmp(sa1),lmbda)

  640 continue

c*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
c      call Dneall_fmf_Multiply(0,
c     >                          psi1,npack1,
c     >                          lmbda, dte,
c     >                          psi2,1.0d0)

      call nwpw_timing_end(3)
      return
      end

