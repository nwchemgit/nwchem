*
* $Id$
*

      subroutine wvfnc_new(wavefunction_filename,
     >             			version,
     >						ngrid,
     >						unita,
     >						ispin,
     >						nein,oddelcfill)
      implicit none
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,nein(2)
      logical      oddelcfill

#include "bafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
 
*     **** local variables ****
      logical value,fractional
      integer MASTER,taskid
      parameter (MASTER=0)

      integer NMAX
      integer filling(2)
      integer fractional_orbitals(2),ne(2)

*     **** external functions ****
      logical  control_fractional
      integer  control_fractional_orbitals
      external control_fractional
      external control_fractional_orbitals


      ne(1) = nein(1)
      ne(2) = nein(2)
      fractional = control_fractional()
      if (fractional) then
         fractional_orbitals(1) = control_fractional_orbitals(1)
         ne(1) = nein(1) + fractional_orbitals(1)
         if (ispin.eq.2) then
           fractional_orbitals(2) = control_fractional_orbitals(2)
           ne(2) = nein(2) + fractional_orbitals(2)
         end if
      end if

      NMAX = ne(1)+ne(2)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        value = BA_push_get(mt_int,8*NMAX,
     >                    'filling',filling(2),filling(1))
        if (.not. value) 
     >    call errquit('wvfnc_new:out of stack memory',0,MA_ERR)


        call sub_wvfnc_new(NMAX,int_mb(filling(1)),
     >                     wavefunction_filename,
     >                     version,
     >                     ngrid,
     >                     unita,
     >                     ispin,
     >                     ne,
     >                     fractional,
     >                     fractional_orbitals,oddelcfill)
        write(luout,*) "random planewave guess, initial psi:",
     >             wavefunction_filename
        write(luout,*) "-   spin, nalpha, nbeta:",ispin,ne
        value = BA_pop_stack(filling(2))
        if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      end if
      call ga_sync()

      return
      end


      subroutine sub_wvfnc_new(NMAX,filling,
     >                         wavefunction_filename,
     >                         version,
     >                         ngrid,
     >                         unita,
     >                         ispin,
     >                         ne,
     >                         fractional,
     >                         frac_orb,oddelcfill)

      implicit none
      integer NMAX
      integer filling(4,NMAX,2)
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)
      logical      fractional
      integer      frac_orb(2)
      logical      oddelcfill

#include "bafdecls.fh"
#include "errquit.fh"
 
      logical value
      integer ierr
      character*255 full_filename

      integer nfft1,nfft2,nfft3,nfft3d,n2ft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l,l1,l2
      integer i,j,k,ms,n,occupation

      double precision p,scale,fx,fy,fz
      double complex cc,cx,sx,zx,zc,rx,ry

*     **** external functions ****
      double precision gcdotc,util_random
      external         gcdotc,util_random
 

      p = util_random(5291999) !*** initialize the random sequence ****

      call getfilling(.true.,ne(1),filling)
      if (ispin.eq.2) call getfilling(.true.,ne(2),filling(1,1,2))

*     ***** constants *****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      nfft3d=(nfft1/2+1)*nfft2*nfft3
      n2ft3d=2*nfft3d
      inc2c = nfft1/2+1
      inc3c =inc2c*nfft2

      scale=1.0d0/dsqrt(2.0d0)
      zx=dcmplx(1.0d0,0.0d0)
      sx=dcmplx(0.0d0,1.0d0)*scale
      cx=dcmplx(1.0d0,0.0d0)*scale


*     ***** allocate wavefunction memory ****
      value = BA_push_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull_hndl,cfull_indx)
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

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
      if (fractional) then
         occupation = ispin
      else
         occupation = -1
      end if
      call iwrite(6,occupation,1)
 
      call dcopy(n2ft3d,0.0d0,0,dcpl_mb(cfull_indx),1)

      if (version.eq.9) then

         do ms=1,ispin
            cc = zx
            dcpl_mb(cfull_indx) = cc
            do k=2,nfft3d
               dcpl_mb(cfull_indx+k-1) = dcpl_mb(cfull_indx+k-1)
     >                                 + dcmplx((0.5d0-util_random(0)),
     >                                          (0.5d0-util_random(0)))
     >                                  /dsqrt(dble(nfft3d))
            end do
            call gctimereverse(nfft1,nfft2,nfft3,dcpl_mb(cfull_indx))
            P=gcdotc(nfft1,nfft2,nfft3,
     >           dcpl_mb(cfull_indx),
     >           dcpl_mb(cfull_indx))
            P=1.0d0/dsqrt(P)
            call dscal(n2ft3d,P,dcpl_mb(cfull_indx),1)
            call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)
         end do

      else

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
        do k=1,nfft3d
          dcpl_mb(cfull_indx+k-1) = dcpl_mb(cfull_indx+k-1)
     >           + dcmplx((0.5d0-util_random(0)),
     >                    (0.5d0-util_random(0)))
     >              /dsqrt(dble(nfft3d))
c         dcpl_mb(cfull_indx+k-1) = dcpl_mb(cfull_indx+k-1)
c    >           + (0.5d0-util_random(0))
        end do
        zc = dcpl_mb(cfull_indx)
        dcpl_mb(cfull_indx) = dcmplx(dble(zc),0.0d0)
        call gctimereverse(nfft1,nfft2,nfft3,dcpl_mb(cfull_indx))
        P=gcdotc(nfft1,nfft2,nfft3,
     >           dcpl_mb(cfull_indx),
     >           dcpl_mb(cfull_indx))
        P=1.0d0/dsqrt(P)
        call dscal(n2ft3d,P,dcpl_mb(cfull_indx),1)

        call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)

        dcpl_mb(cfull_indx+l1)=dcmplx(0.0d0,0.0d0)
c        call dcopy(n2ft3d,0.0d0,0,dcpl_mb(cfull_indx),1)
 

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
            do 125 K=0,NFFT3D-1
              dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >           + dcmplx((0.5d0-util_random(0)),
     >                    (0.5d0-util_random(0)))
     >              /dsqrt(dble(nfft3d))
c             dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
c    >                   + (0.5d0-util_random(0))
  125       continue
            zc = dcpl_mb(cfull_indx)
            dcpl_mb(cfull_indx) = dcmplx(dble(zc),0.0d0)
            call gctimereverse(nfft1,nfft2,nfft3,dcpl_mb(cfull_indx))
            P=gcdotc(nfft1,nfft2,nfft3,
     >               dcpl_mb(cfull_indx),
     >               dcpl_mb(cfull_indx))
            P=1.0d0/dsqrt(P)
            call dscal(n2ft3d,P,dcpl_mb(cfull_indx),1)
          ENDIF


          call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)

          dcpl_mb(cfull_indx+l1)=dcmplx(0.0d0,0.0d0)
          dcpl_mb(cfull_indx+l2)=dcmplx(0.0d0,0.0d0)
c          call dcopy(n2ft3d,0.0d0,0,dcpl_mb(cfull_indx),1)
  130   CONTINUE
 

  150 CONTINUE
      end if

c     **** add occupation ****
      if (fractional) then
         fx = 1.0d0
         fy = 0.0d0
         fz = 0.5d0
         do ms=1,ispin
           do n=1,ne(ms)
             if (oddelcfill.and.(n.eq.frac_orb(ms))) then
                call dwrite(6,fz,1)
             else if (n.le.frac_orb(ms)) then
                call dwrite(6,fy,1)
             else
                call dwrite(6,fx,1)
             end if
           end do
         end do
      end if

 
      call closefile(6)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = BA_pop_stack(cfull_hndl)
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      
      return 
      end
 

