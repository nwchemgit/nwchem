*
* $Id: cpsi_lmbda.F,v 1.3 2002-01-15 17:27:20 edo Exp $
*

      subroutine cpsi_lmbda(ispin,ne,nfft3d,neall,nbrill,n,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr)

      implicit none

      integer ispin,ne(2),nfft3d,neall,nbrill,n
      complex*16 psi1(nfft3d,neall,nbrill)
      complex*16 psi2(nfft3d,neall,nbrill)
      real*8     dte
      complex*16 lmbda(n,n,2,nbrill)
      complex*16 tmp(*)
      integer	 ierr

c     complex*16 q22(n,n,2,nbrill)
c     complex*16 q21(n,n,2,nbrill)
c     complex*16 q12(n,n,2,nbrill)
c     complex*16 q11(n,n,2,nbrill)

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical ortho
      integer n1(2),n2(2)
      integer nb
      integer i,j,ii,jj,ms,it
      integer index,indext,nn
      integer s11,s12,s21,s22,st1,st2,sa1,sa0
      integer sl(2)
      integer iadiff
      real*8 adiff
      complex*16 alpha,zero,one
      complex*16 sum
*     ::::  iteration limit and tolerence for non-liner equations  ::::
      integer itrlmd
      real*8  convg
      parameter (itrlmd=20, convg=1.d-10)

*     **** external functions ****
      integer  izamax
      external izamax


      call Parallel_taskid(taskid)

      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
*     **** set the matrices dimensions ****
      nn   = n**2
      
      s11  = 0*nn + 1
      s12  = 1*nn + 1
      s21  = 2*nn + 1
      s22  = 3*nn + 1
      sa0  = 4*nn + 1
      sa1  = 5*nn + 1
      st1  = 6*nn + 1
      st2  = 7*nn + 1

      

      sl(1)  = 0*nn + 1
      sl(2)  = 1*nn + 1
      

      n1(1)=1
      n2(1)=ne(1)
      n1(2)=ne(1)+1
      n2(2)=ne(1)+ne(2)
      
      
      call zcopy(2*nn*nbrill,0.0d0,0,lmbda,1)
      call zcopy(8*nn,0.0d0,0,tmp,1)
     

*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
c*     **** calculated Overlap matrices ****
c      call Stfl_ggq_dot(ispin,ne,nffft3d,neall,nbrill,n
c     >                  psi2,psi2,q22)
c      call Stfl_ggq_dot(ispin,ne,nffft3d,neall,nbrill,n
c     >                  psi1,psi2,q12)
c      call Stfl_ggq_dot(ispin,ne,nffft3d,neall,nbrill,n
c     >                  psi1,psi2,q11)

      ortho = .true.
      do 650 nb=1,nbrill
      do 640 ms=1,ispin
        IF(ne(ms).le.0) GO TO 640
        do II=n1(ms),n2(ms)
          I   = II - n1(ms) + 1
          index = (i-1) + (i-1)*n
          call Cram_cc_zdot(nb,psi2(1,ii,nb),psi2(1,ii,nb),sum)
          tmp(s22+index)=(1.0d0-sum)*0.5d0/dte

          call Cram_cc_zdot(nb,psi2(1,ii,nb),psi1(1,ii,nb),sum)
          tmp(s21+index)=(1.0d0-sum)*0.5d0
          tmp(s12+index)=dconjg(tmp(s21+index))

          call Cram_cc_zdot(nb,psi1(1,ii,nb),psi1(1,ii,nb),sum)
          tmp(s11+index)= (-sum)*0.5d0*dte


           do JJ=II+1,N2(MS)
             J   = JJ - n1(ms) + 1
             index  = (i-1) + (j-1)*n
             indext = (j-1) + (i-1)*n

             call Cram_cc_zdot(nb,psi2(1,ii,nb),psi2(1,jj,nb),sum)
             tmp(s22+index)=(-sum)*0.5d0/dte

             call Cram_cc_zdot(nb,psi2(1,ii,nb),psi1(1,jj,nb),sum)
             tmp(s21+index)= (-sum)*0.5d0

             call Cram_cc_zdot(nb,psi1(1,ii,nb),psi2(1,jj,nb),sum)
             tmp(s12+index)= (-sum)*0.5d0

             call Cram_cc_zdot(nb,psi1(1,ii,nb),psi1(1,jj,nb),sum)
             tmp(s11+index)= (-sum)*0.5d0*dte

             tmp(s22+indext)=dconjg(tmp(s22+index))
             tmp(s21+indext)=dconjg(tmp(s12+index))
             tmp(s12+indext)=dconjg(tmp(s21+index))
             tmp(s11+indext)=dconjg(tmp(s11+index))
          end do
        end do

        call zcopy(nn,tmp(s22),1,tmp(sa0),1)

        do it=1,itrlmd
          call zcopy(nn,tmp(s22),1,tmp(sa1),1)

          call ZGEMM('N','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(s21),n,
     >                tmp(sa0),n,
     >                zero,
     >                tmp(st1),n)

          call ZGEMM('C','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(sa0),n,
     >                tmp(s12),n,
     >                zero,
     >                tmp(st2),n)

          call ZAXPY(nn,one,tmp(st2),1,tmp(st1),1)

          call ZAXPY(nn,one,
     >               tmp(st1),1,
     >               tmp(sa1),1)

          call ZGEMM('N','C',ne(MS),ne(MS),ne(MS),one,
     >                tmp(s11),n,
     >                tmp(sa0),n,
     >                zero,
     >                tmp(st1),n)

          call ZGEMM('N','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(sa0),n,
     >                tmp(st1),n,
     >                zero,
     >                tmp(st2),n)

          call ZAXPY(nn,one,
     >                tmp(st2),1,
     >                tmp(sa1),1)

          call zcopy(nn,tmp(sa1),1,tmp(st1),1)
          call ZAXPY(nn,-one,
     >                tmp(sa0),1,
     >                tmp(st1),1)

          iadiff = st1 - 1 + izamax(n*ne(ms),tmp(st1),1)
cedo          write(0,*) 8*nn,iadiff,'tmp iad ',tmp(iadiff)
          adiff=  dble(tmp(iadiff))**2
     >         + dimag(tmp(iadiff))**2

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
        ortho = .false.

  630   continue
        call zcopy(nn,tmp(sa1),1,lmbda(1,1,ms,nb),1)
        call zcopy(8*nn,0.0d0,0,tmp,1)

  640 continue
  650 continue



*:::::::::::::::::  correction due to the constraint  :::::::::::::::::
      do nb=1,nbrill
      do ms=1,ispin

        do ii=n1(ms),n2(ms)
           i   = ii - n1(ms) + 1

           do jj=n1(ms),n2(ms)
              j   = jj - n1(ms) + 1

              lmbda(i,j,ms,nb) = dconjg(lmbda(i,j,ms,nb))
              alpha            = dte*lmbda(i,j,ms,nb)

              call Cram_cc_zaxpy(nb,alpha,psi1(1,jj,nb),psi2(1,ii,nb))
           end do
        end do
      end do
      end do

*     **** Ortho CheckMake ****
      if (.not.ortho) then
        call Stfl_g_OrthoCheckMake(ispin,ne,nfft3d,neall,nbrill,psi2)
      end if

      return
      end

