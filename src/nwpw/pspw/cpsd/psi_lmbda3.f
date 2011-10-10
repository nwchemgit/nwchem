*
* $Id$
*

      subroutine psi_lmbda3(ispin,ne,nemax,npack1,
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
      integer	 ierr

      integer MASTER
      parameter (MASTER=0)

*     **** local variables ****
      integer taskid
      integer n1(2),n2(2)
      integer i,j,ii,jj,ms,it,ii1
      integer index,indext,n,nn
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      integer sl(2)
      real*8  adiff,alpha,tmp1(1000000)
*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd,idamax
      real*8  convg
      parameter (itrlmd=200, convg=1.0d-15)


      call nwpw_timing_start(3)

      call Parallel_taskid(taskid)

      n    = ne(1)
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

      n1(1)=1
      n2(1)=ne(1)
      n1(2)=ne(1)+1
      n2(2)=ne(1)+ne(2)
      
*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      DO 640 ms=1,ispin
        IF(ne(ms).le.0) GO TO 640

*       ***** compute the overlap matrices ****
        call Pack_ccm_sym_dot2(1,n,ne(ms),
     >                          psi2(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(s22))
        call Pack_ccmn_dot(1,n,ne(ms),
     >                          psi2(1,n1(ms)),
     >                          psi1(1,n1(ms)),
     >                          tmp(s21))
        call Pack_ccmn_dot(1,n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi2(1,n1(ms)),
     >                          tmp(s12))
        call Pack_ccm_sym_dot2(1,n,ne(ms),
     >                          psi1(1,n1(ms)),
     >                          psi1(1,n1(ms)),
     >                          tmp(s11))


*       ***** scale the overlap matrices ****
        do i=1,ne(ms)
          index = (i-1) + (i-1)*n
      
          tmp(s22+index)=(1.0d0-tmp(s22+index))*0.5d0/dte
          tmp(s21+index)=(1.0d0-tmp(s21+index))*0.5d0
          tmp(s12+index)=(1.0d0-tmp(s12+index))*0.5d0
          tmp(s11+index)= (-tmp(s11+index))*0.5d0*dte
          
           do j=i+1,ne(ms)
             index  = (i-1) + (j-1)*n
             indext = (j-1) + (i-1)*n
   
             tmp(s22+index)= (-tmp(s22+index))*0.5d0/dte
             tmp(s21+index)= (-tmp(s21+index))*0.5d0
             tmp(s12+index)= (-tmp(s12+index))*0.5d0
             tmp(s11+index)= (-tmp(s11+index))*0.5d0*dte

             tmp(s22+indext)=tmp(s22+index)
             tmp(s21+indext)=tmp(s12+index)
             tmp(s12+indext)=tmp(s21+index)
             tmp(s11+indext)=tmp(s11+index)
          end do
        end do

        call dcopy(nn,tmp(s22),1,tmp(sa0),1)

        do it=1,itrlmd
          CALL dcopy(nn,tmp(s22),1,tmp(sa1),1)
c         CALL DMMUL(n,ne(MS), tmp(s21), tmp(sa0), tmp(st1))
c         CALL DMMUL(n,ne(MS), tmp(sa0), tmp(s12), tmp(st2))
c         CALL DMADD(n,ne(MS), tmp(st1), tmp(st2), tmp(st1))
c         CALL DMADD(n,ne(MS), tmp(st1), tmp(sa1), tmp(sa1))

c         CALL DMMUL(n,ne(MS), tmp(s11), tmp(sa0), tmp(st1))
c         CALL DMMUL(n,ne(MS), tmp(sa0), tmp(st1), tmp(st2))
c         CALL DMADD(n,ne(MS), tmp(st2), tmp(sa1), tmp(sa1))
c         CALL DMSUB(n,ne(MS), tmp(sa1), tmp(sa0), tmp(st1))

          call DGEMM('N','N',ne(ms),ne(ms),ne(ms),
     >                (1.0d0),
     >                tmp(s21),n,
     >                tmp(sa0),n,
     >                (1.0d0),
     >                tmp(sa1),n)
          call DGEMM('N','N',ne(ms),ne(ms),ne(ms),
     >                (1.0d0),
     >                tmp(sa0),n,
     >                tmp(s12),n,
     >                (1.0d0),
     >                tmp(sa1),n)
         
          call DGEMM('N','N',ne(ms),ne(ms),ne(ms),
     >                (1.0d0),
     >                tmp(s11),n,
     >                tmp(sa0),n,
     >                (0.0d0),
     >                tmp(st1),n)
         
          call DGEMM('N','N',ne(ms),ne(ms),ne(ms),
     >                (1.0d0),
     >                tmp(sa0),n,
     >                tmp(st1),n,
     >                (1.0d0),
     >                tmp(sa1),n)
          CALL dcopy(nn,tmp(sa1),1,tmp(st1),1)
          call daxpy(nn,(-1.0d0),tmp(sa0),1,tmp(st1),1)

          adiff=tmp(st1 - 1 + (idamax(n*ne(ms),tmp(st1),1)))
          if(adiff.lt.convg) GO TO 630
          call dcopy(n*ne(ms),tmp(sa1),1,tmp(sa0),1)
        end do

        ierr=10
        if (taskid.eq.MASTER) then
          WRITE(6,*) 
     >     'Warning: Lagrange Multiplier tolerance too high:',adiff
          WRITE(6,*) '        +Try using a smaller time step'
          WRITE(6,*) '        +Gram-Schmidt being performed, spin:',ms
        end if
c        call Dneall_f_ortho(ms,psi2,npack1)
        call Dneall_f_GramSchmidt(ms,psi2,npack1)
c        call Grsm_g_MakeOrtho(npack1,ne(ms),psi2(1,n1(ms)))

C       return
  630   continue
        call dcopy(n*ne(ms),tmp(sa1),1,lmbda(sl(ms)),1)
  640 continue

*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      do ms=1,ispin
       IF(ne(ms).le.0) GO TO 650
       call DGEMM('N','N',2*npack1,ne(ms),ne(ms),
     >              dte,
     >              psi1(1,n1(ms)),2*npack1,
     >              lmbda(sl(ms)),n,
     >              (1.0d0),
     >              psi2(1,n1(ms)),2*npack1)
c        do ii=n1(ms),n2(ms)
c           i=ii-n1(ms)+1
c
c           do jj=n1(ms),n2(ms)
c              j=jj-n1(ms)+1
c
c              index = (i-1) + (j-1)*n
c              alpha = dte*lmbda(sl(ms) + index)
c
cc             do k=1,nfft3d
cc                psi2(k,ii) = psi2(k,ii)  + alpha*psi1(k,jj)
cc             end do
c              call Pack_cc_daxpy(1,alpha,psi1(1,jj),psi2(1,ii))
c           end do
c        end do

  650  continue
      end do

      call nwpw_timing_end(3)

      return
      end

