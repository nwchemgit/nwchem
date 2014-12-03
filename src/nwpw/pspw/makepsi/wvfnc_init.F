*
* $Id$
*


      logical function wvfnc_init(rtdb)
      implicit none
      integer rtdb
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"

      logical value
      integer version
 
      integer ierr

      integer NMAX
      parameter (NMAX=501)
      integer filling(4,NMAX,2)
      integer ne(2),ispin
      logical  restricted

      character*50 cell_name,wavefunction_filename
      character*50 rtdb_unita,rtdb_ngrid
      character*255 full_filename


      integer ngrid(3)
      integer nfft1,nfft2,nfft3,nfft3d,n2ft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l,l1,l2
      integer i,j,k,ms,n,occupation

      double precision unita(3,3),p,scale
      double complex cc,cx,sx,zx

*     **** external functions ****
      double precision gcdotc,util_random
      external         gcdotc,util_random
 
      p = util_random(5291999)  !*** initialize the random sequence ***
      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'wvfnc:wavefunction_filename',
     >                  1,wavefunction_filename)

      value = btdb_get(rtdb,'wvfnc:nup',  mt_int,1,ne(1))
      value = btdb_get(rtdb,'wvfnc:ndown',mt_int,1,ne(2))
      value = btdb_get(rtdb,'wvfnc:up_filling',
     >                 mt_int,(4*ne(1)+1),filling(1,1,1))
      value = btdb_get(rtdb,'wvfnc:down_filling',
     >                 mt_int,(4*ne(2)+1),filling(1,1,2))
      value = btdb_get(rtdb,'wvfnc:restricted', mt_log,1,restricted)
      ispin=2
      if (restricted) ispin=1


*     **** get cell information ****
      value = btdb_cget(rtdb,'wvfnc:cell_name',
     >                  1,cell_name)

      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      rtdb_ngrid = cell_name(1:l)//':ngrid'

      value = btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      value = btdb_get(rtdb,rtdb_ngrid,mt_int,3,ngrid)
     


*     ***** constants *****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      nfft3d=(nfft1/2+1)*nfft2*nfft3
      n2ft3d=2*nfft3d
      inc2c = nfft1/2+1
      inc3c =inc2c*nfft2

      scale=1.0d0/dsqrt(2.0d0)
      zx=(1.0d0,0.0d0)
      sx=(0.0d0,1.0d0)*scale
      cx=(1.0d0,0.0d0)*scale



*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull_hndl,cfull_indx)

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
      occupation = -1
      call iwrite(6,occupation,1)
 
      call dcopy(n2ft3d,0.0d0,0,dcpl_mb(cfull_indx),1)

      do 150 ms=1,ispin
 
        if (ne(ms).LE.0) GO TO 150
c       L1=INC3C*K3(1)+INC2C*K2(1)+K1(1)+1
c       CFULL(L1)=CC(1)
        l1= inc3c*filling(3,1,ms)
     >    + inc2c*filling(2,1,ms)
     >    +       filling(1,1,ms)
        if (filling(4,1,ms).lt.0) cc=sx
        if (filling(4,1,ms).eq.0) cc=zx
        if (filling(4,1,ms).gt.0) cc=cx
        dcpl_mb(cfull_indx+l1) = cc

        write(6,1011) ms,1,filling(1,1,ms),
     >                     filling(2,1,ms),
     >                     filling(3,1,ms),
     >                     filling(4,1,ms)
        call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)

        dcpl_mb(cfull_indx+l1)=(0.0d0,0.0d0)
 
        DO 130 N=2,NE(ms)
c         L1=INC3C*K3(N)+INC2C*K2(N)+K1(N)+1
          l1= inc3c*filling(3,n,ms)
     >      + inc2c*filling(2,n,ms)
     >      +       filling(1,n,ms)
          if (filling(4,n,ms).lt.0) cc=sx
          if (filling(4,n,ms).eq.0) cc=zx
          if (filling(4,n,ms).gt.0) cc=cx
          L2=L1
          dcpl_mb(cfull_indx+l1)=cc
          if (filling(1,n,ms).eq.0) then
c           L2=INC3C*MOD(M3-K3(N),M3)+INC2C*MOD(M2-K2(N),M2)+K1(N)+1
c           CFULL(L2)=DCONJG(CC(N))
            l2 = inc3c*mod(nfft3-filling(3,n,ms),nfft3)
     >         + inc2c*mod(nfft2-filling(2,n,ms),nfft2)
     >         +       filling(1,n,ms)
            dcpl_mb(cfull_indx+l2)=dconjg(cc)
          end if

          IF((ABS(filling(4,n,ms)).gt.1)) THEN
            do 125 K=1,NFFT3D-1
c             CFULL(K)=CFULL(K)+(0.5d0-RAND())/DSQRT(dble(NFFT3D))
              dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >                   + (0.5d0-util_random(0))/dsqrt(dble(nfft3d))
  125       continue
            call gctimereverse(nfft1,nfft2,nfft3,dcpl_mb(cfull_indx))
            P=gcdotc(nfft1,nfft2,nfft3,
     >               dcpl_mb(cfull_indx),
     >               dcpl_mb(cfull_indx))
            write(luout,*) ' Random noise is added to this psi. P=',P
            P=1.0d0/dsqrt(P)
            call dscal(n2ft3d,P,dcpl_mb(cfull_indx),1)
          ENDIF


 1011     FORMAT(' writing <spin=',I2,',psi=',I3,'>    -     (',
     >           3I3,' : ',I3,')')
          WRITE(luout,1011) ms,n,filling(1,n,ms),
     >                       filling(2,n,ms),
     >                       filling(3,n,ms),
     >                       filling(4,n,ms)
          call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)

          dcpl_mb(cfull_indx+l1)=(0.0d0,0.0d0)
          dcpl_mb(cfull_indx+l2)=(0.0d0,0.0d0)
  130   CONTINUE
 

  150 CONTINUE
 
      call closefile(6)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = BA_free_heap(cfull_hndl)
      IF(IERR.EQ.0) THEN
        WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      ELSE
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
        value = .false.
      ENDIF
      call nwpw_message(4)
      
      wvfnc_init = value
      return 
      end
 


      double precision FUNCTION GCDOTC(NFFT1,NFFT2,NFFT3,A,B)
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
     
      sum  = 0.0d0
      inc2 = nfft1/2+1
      inc3 = inc2*nfft2

      do k3=1,nfft3
      do k2=1,nfft2
         do k1=2,nfft1/2+1
            index = (k3-1)*inc3 + (k2-1)*inc2 + k1
            sum = sum + dconjg(A(index))*B(index)
         end do 
      end do 
      end do 
      sum = 2.0d0*sum
     
      do k3=1,nfft3
      do k2=1,nfft2
         index = (k3-1)*inc3 + (k2-1)*inc2 + 1
         sum = sum + dconjg(A(index))*B(index)
      end do 
      end do 
      
      GCDOTC=sum
      return
      end 




      subroutine gctimereverse(nfft1,nfft2,nfft3,A)
      implicit none
      integer nfft1,nfft2,nfft3
      complex*16 A(*)

*     **** local variables ****
      integer i3,j3,k3
      integer i2,j2,k2
      integer index1,index2

      A(1) = dcmplx(dble(A(1)),0.0d0)
      do k3=1,(nfft3/2-1)
         i3 =  k3
         j3 = -k3
         if (i3.lt.0) i3 = i3 + nfft3
         if (j3.lt.0) j3 = j3 + nfft3
         i3 = i3 + 1
         j3 = j3 + 1
         index1 = (i3-1)*(nfft1/2+1)*nfft2 + 1
         index2 = (j3-1)*(nfft1/2+1)*nfft2 + 1

         A(index2) = dconjg(A(index1))

      end do
      index1 = (nfft3/2)*(nfft1/2+1)*nfft2 + 1
      A(index1) = dcmplx(0.0d0,0.0d0)

      do k3=(-nfft3/2+1),(nfft3/2-1)
      do k2=1,(nfft2/2-1)
         i2 =  k2
         i3 =  k3
         j2 = -k2
         j3 = -k3
         if (i2.lt.0) i2 = i2 + nfft2
         if (i3.lt.0) i3 = i3 + nfft3
         if (j2.lt.0) j2 = j2 + nfft2
         if (j3.lt.0) j3 = j3 + nfft3
         i2 = i2 + 1
         i3 = i3 + 1
         j2 = j2 + 1
         j3 = j3 + 1
         index1 = (i3-1)*(nfft1/2+1)*nfft2 + (i2-1)*(nfft1/2+1) + 1
         index2 = (j3-1)*(nfft1/2+1)*nfft2 + (j2-1)*(nfft1/2+1) + 1

         A(index2) = dconjg(A(index1))

      end do
      end do
      k3 = nfft3/2
      do k2=1,(nfft2/2)
         i2 =  k2
         i3 =  k3
         if (i2.lt.0) i2 = i2 + nfft2
         if (i3.lt.0) i3 = i3 + nfft3
         i2 = i2 + 1
         i3 = i3 + 1
         index1 = (i3-1)*(nfft1/2+1)*nfft2 + (i2-1)*(nfft1/2+1) + 1

         A(index1) = dcmplx(0.0d0,0.0d0)
      end do
      k2 = nfft2/2
      do k3=(-nfft3/2+1),(nfft3/2)
         i2 =  k2
         i3 =  k3
         if (i2.lt.0) i2 = i2 + nfft2
         if (i3.lt.0) i3 = i3 + nfft3
         i2 = i2 + 1
         i3 = i3 + 1
         index1 = (i3-1)*(nfft1/2+1)*nfft2 + (i2-1)*(nfft1/2+1) + 1

         A(index1) = dcmplx(0.0d0,0.0d0)
      end do


      return
      end
      

