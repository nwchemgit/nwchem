*
* $Id: cpsi_lmbda.F,v 1.1 2001-12-14 20:02:44 bylaska Exp $
*

      subroutine cpsi_lmbda(ispin,ne,nemax,nfft3d,nbrill,
     >                     psi1,psi2,
     >                     dte,
     >                     lmbda,tmp,ierr)

      implicit none

      integer ispin,ne(2),nemax,nfft3d,nbrill
      complex*16 psi1(nfft3d,nemax)
      complex*16 psi2(nfft3d,nemax)
      real*8     dte
      complex*16 lmbda(nemax,nemax,2,nbrill)
      complex*16 tmp(*)
      integer	 ierr


*     **** local variables ****
      integer n1(2),n2(2)
      integer neband1,nb
      integer i,j,ii,jj,iii,jjj,ms,it,k
      integer index,indext,n,nn
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


      one  = dcmplx(1.0d0,0.0d0)
      zero = dcmplx(0.0d0,0.0d0)
*     **** set the matrices dimensions ****
      n    = nemax
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
      neband1 = ne(1)+ne(2)
      
      
      call zcopy(nemax*nemax*2*nbrill,0.0d0,0,lmbda,1)
      call zcopy(8*nn,0.0d0,0,tmp,1)
     

*::::::::::::::::::::::  Lagrangian multipliers  ::::::::::::::::::::::
      do 650 nb=1,nbrill
      do 640 ms=1,ispin
        IF(ne(ms).le.0) GO TO 640
        do II=n1(ms),n2(ms)
          I   = II - n1(ms) + 1
          iii = ii + (nb-1)*neband1
          index = (i-1) + (i-1)*n
          call Cram_cc_zdot(nb,psi2(1,iii),psi2(1,iii),sum)
          tmp(s22+index)=(1.0d0-sum)*0.5d0/dte

          call Cram_cc_zdot(nb,psi2(1,iii),psi1(1,iii),sum)
          tmp(s21+index)=(1.0d0-sum)*0.5d0
          tmp(s12+index)=dconjg(tmp(s21+index))

          call Cram_cc_zdot(nb,psi1(1,iii),psi1(1,iii),sum)
          tmp(s11+index)= (-sum)*0.5d0*dte


           do JJ=II+1,N2(MS)
             J   = JJ - n1(ms) + 1
             jjj = jj + (nb-1)*neband1
             index  = (i-1) + (j-1)*n
             indext = (j-1) + (i-1)*n

             call Cram_cc_zdot(nb,psi2(1,iii),psi2(1,jjj),sum)
             tmp(s22+index)=(-sum)*0.5d0/dte

             call Cram_cc_zdot(nb,psi2(1,iii),psi1(1,jjj),sum)
             tmp(s21+index)= (-sum)*0.5d0

             call Cram_cc_zdot(nb,psi1(1,iii),psi2(1,jjj),sum)
             tmp(s12+index)= (-sum)*0.5d0

             call Cram_cc_zdot(nb,psi1(1,iii),psi1(1,jjj),sum)
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

*         call ZGEMUL(tmp(s21),n,'N',
*    >                tmp(sa0),n,'N',
*    >                tmp(st1),n, ne(MS),ne(MS),ne(MS))
          call ZGEMM('N','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(s21),n,
     >                tmp(sa0),n,
     >                zero,
     >                tmp(st1),n)


*         call ZGEMUL(tmp(sa0),n,'C',
*    >                tmp(s12),n,'N',
*    >                tmp(st2),n, ne(MS),ne(MS),ne(MS))
          call ZGEMM('C','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(sa0),n,
     >                tmp(s12),n,
     >                zero,
     >                tmp(st2),n)

*         call ZGEADD(tmp(st1),n,'N',
*    >                tmp(st2),n,'N',
*    >                tmp(st1),n, ne(MS),ne(MS))
          call ZAXPY(nn,one,tmp(st2),1,tmp(st1),1)

*         call ZGEADD(tmp(st1),n,'N',
*    >                tmp(sa1),n,'N',
*    >                tmp(sa1),n, ne(MS),ne(MS))
          call ZAXPY(nn,one,
     >               tmp(st1),1,
     >               tmp(sa1),1)

*         call ZGEMUL(tmp(s11),n,'N',
*    >                tmp(sa0),n,'C',
*    >                tmp(st1),n, ne(MS),ne(MS),ne(MS))
          call ZGEMM('N','C',ne(MS),ne(MS),ne(MS),one,
     >                tmp(s11),n,
     >                tmp(sa0),n,
     >                zero,
     >                tmp(st1),n)

*         call ZGEMUL(tmp(sa0),n,'N',
*    >                tmp(st1),n,'N',
*    >                tmp(st2),n, ne(MS),ne(MS),ne(MS))
          call ZGEMM('N','N',ne(MS),ne(MS),ne(MS),one,
     >                tmp(sa0),n,
     >                tmp(st1),n,
     >                zero,
     >                tmp(st2),n)

*         call ZGEADD(tmp(st2),n,'N',
*    >                tmp(sa1),n,'N',
*    >                tmp(sa1),n, ne(MS),ne(MS))
          call ZAXPY(nn,one,
     >                tmp(st2),1,
     >                tmp(sa1),1)


*         call ZGESUB(tmp(sa1),n,'N',
*    >                tmp(sa0),n,'N',
*    >                tmp(st1),n, ne(MS),ne(MS))
          call zcopy(nn,tmp(sa1),1,tmp(st1),1)
          call ZAXPY(nn,-one,
     >                tmp(sa0),1,
     >                tmp(st1),1)


          iadiff = st1 - 1 + izamax(n*ne(ms),tmp(st1),1)
          adiff=  dble(tmp(iadiff))**2
     >         + dimag(tmp(iadiff))**2

          if(adiff.lt.convg) GO TO 630
          call zcopy(nn,tmp(sa1),1,tmp(sa0),1)
        end do
        ierr=10
        WRITE(6,*) 'ierr=10',adiff
C       return

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
           iii = ii + (nb-1)*neband1

           do jj=n1(ms),n2(ms)
              j   = jj - n1(ms) + 1
              jjj = jj + (nb-1)*neband1

              lmbda(i,j,ms,nb) = dconjg(lmbda(i,j,ms,nb))
              alpha            = dte*lmbda(i,j,ms,nb)

              call Cram_cc_zaxpy(nb,alpha,psi1(1,jjj),psi2(1,iii))
           end do
        end do
      end do
      end do

      return
      end

