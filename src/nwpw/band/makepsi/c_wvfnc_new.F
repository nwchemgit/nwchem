*
* $Id: c_wvfnc_new.F,v 1.2 2001-12-15 03:11:01 bylaska Exp $
*
      subroutine c_wvfnc_new(wavefunction_filename,
     >             			version,
     >						ngrid,
     >						unita,
     >						ispin,
     >						ne,
     >                      nbrillioun)
      implicit none
      character*30 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2),nbrillioun

#include "mafdecls.fh"
 
*     **** local variables ****
      logical value
      integer MASTER,taskid
      parameter (MASTER=0)

      integer NMAX
      integer filling(2)

      NMAX = ne(1)+ne(2)

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        value = MA_push_get(mt_int,8*NMAX,
     >                    'filling',filling(2),filling(1))
        if (.not. value) call errquit('out of stack memory',0)


        call sub_c_wvfnc_new(NMAX,int_mb(filling(1)),
     >					wavefunction_filename,
     >             		version,
     >					ngrid,
     >					unita,
     >					ispin,
     >					ne,
     >                  nbrillioun)
        write(*,*) "random planewave guess, initial psi:",
     >             wavefunction_filename
        write(*,*) "-   spin, nalpha, nbeta:",ispin,ne
        value = MA_pop_stack(filling(2))
        if (.not. value) call errquit('popping stack memory',0)
      end if
      call ga_sync()

      return
      end


      subroutine sub_c_wvfnc_new(NMAX,filling,
     >                      wavefunction_filename,
     >             			version,
     >						ngrid,
     >						unita,
     >						ispin,
     >						ne,
     >                      nbrillioun)

      implicit none
      integer NMAX
      integer filling(4,NMAX,2)
      character*30 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2),nbrillioun

#include "mafdecls.fh"
 
      logical value
      integer ierr
      character*255 full_filename

      integer nfft1,nfft2,nfft3,nfft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l,l1,l2
      integer i,j,k,ms,n,nb

      double precision p,scale
      double complex cc,cx,sx,zx

*     **** external functions ****
      double precision gcdot,util_random
      external         gcdot,util_random
 

      call getfilling(ne(1),filling)
      if (ispin.eq.2) call getfilling(ne(2),filling(1,1,2))

*     ***** constants *****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      nfft3d=(nfft1)*nfft2*nfft3
      inc2c = nfft1
      inc3c =inc2c*nfft2

      scale=1.0d0/dsqrt(2.0d0)
      zx=(1.0d0,0.0d0)
      sx=(0.0d0,1.0d0)*scale
      cx=(1.0d0,0.0d0)*scale


*     ***** allocate wavefunction memory ****
      value = MA_push_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull_hndl,cfull_indx)
      if (.not. value) call errquit('out of stack memory',0)

*     **** modularize the filling ****
      do ms=1,ispin
        do n=1,ne(ms)
           i = filling(1,n,ms)
           j = filling(2,n,ms)
           k = filling(3,n,ms)
           filling(1,n,ms) = mod(i+inc2c,inc2c)
           filling(2,n,ms) = mod(j+nfft2,nfft2)
           filling(3,n,ms) = mod(k+nfft3,nfft3)
        end do
      end do


 
*:::::::::::::::::  write wavefunction in CPMDV3 format ::::::::::::::
      call util_file_name_noprefix(wavefunction_filename,.false.,
     >                                          .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      call openfile(6,full_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrillioun,1)
 
      do nb=1,nbrillioun
        call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(cfull_indx),1)

        do 150 ms=1,ispin
 
        if (ne(ms).LE.0) GO TO 150
        l1= inc3c*filling(3,1,ms)
     >    + inc2c*filling(2,1,ms)
     >    +       filling(1,1,ms)
        if (filling(4,1,ms).lt.0) cc=sx
        if (filling(4,1,ms).eq.0) cc=zx
        if (filling(4,1,ms).gt.0) cc=cx
        dcpl_mb(cfull_indx+l1) = cc

        call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d)

        dcpl_mb(cfull_indx+l1)=(0.0d0,0.0d0)
 
        DO 130 N=2,NE(ms)
          l1= inc3c*filling(3,n,ms)
     >      + inc2c*filling(2,n,ms)
     >      +       filling(1,n,ms)
          if (filling(4,n,ms).lt.0) cc=sx
          if (filling(4,n,ms).eq.0) cc=zx
          if (filling(4,n,ms).gt.0) cc=cx
          L2=L1
          dcpl_mb(cfull_indx+l1)=cc
          if (filling(1,n,ms).eq.0) then
            l2 = inc3c*mod(nfft3-filling(3,n,ms),nfft3)
     >         + inc2c*mod(nfft2-filling(2,n,ms),nfft2)
     >         +       filling(1,n,ms)
            dcpl_mb(cfull_indx+l2)=dconjg(cc)
          end if

          IF((ABS(filling(4,n,ms)).gt.1)) THEN
            do 125 K=2,NFFT3D-1
c             CFULL(K)=CFULL(K)+(0.5d0-RAND())/DSQRT(dble(NFFT3D))
              dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >                   + (0.5d0-util_random(9))/dsqrt(dble(nfft3d))
  125       continue
            P=gcdot(nfft1,nfft2,nfft3,
     >               dcpl_mb(cfull_indx),
     >               dcpl_mb(cfull_indx))
            P=1.0d0/dsqrt(P)
            call dscal(2*nfft3d,P,dcpl_mb(cfull_indx),1)
          ENDIF
          call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d)

          call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(cfull_indx),1)
          dcpl_mb(cfull_indx+l1)=(0.0d0,0.0d0)
          dcpl_mb(cfull_indx+l2)=(0.0d0,0.0d0)
  130   CONTINUE
 

  150   CONTINUE
 
      end do
      call closefile(6)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = MA_pop_stack(cfull_hndl)
      if (.not. value) call errquit('popping stack memory',0)
      
      return 
      end

      double precision FUNCTION GCDOT(NFFT1,NFFT2,NFFT3,A,B)
      implicit none
      integer nfft1,nfft2,nfft3
      complex*16 A(*),B(*)
*     =============================================================
*     INNER PRODUCT OF COMPLEX VECTORS IN THE RECIPROCAL SPACE
*     ASSUMING GIVEN VEXTORS HAVE AN EVEN PARITY,
*         A(-G) * CONJG(B(-G)) = CONJG(A(G)) * B(G)
*     VECTORS A AND B MUST BE IN THE COMPRESSED ORDER.
*     =============================================================
      double precision sum
      integer k1,k2,k3,inc2,inc3,index

      inc2 = nfft1
      inc3 = inc2*nfft2

      sum  = 0.0d0
      do k3=1,nfft3
      do k2=1,nfft2
      do k1=1,nfft1
         index = (k3-1)*inc3 + (k2-1)*inc2 + k1
         sum = sum + dconjg(A(index))*B(index)
      end do
      end do
      end do

      GCDOT=sum
      return
      end

 

