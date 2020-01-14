*
* $Id$
*

*     **********************************************
*     *                                            *
*     *              c_wvfnc_new                   *
*     *                                            *
*     **********************************************
      subroutine c_wvfnc_new(wavefunction_filename,
     >                       version,
     >                       ngrid,
     >                       unita,
     >                       spin_orbit,
     >                       ispin,
     >                       nein,
     >                       nbrillioun,oddelcfill)
      implicit none
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      logical      spin_orbit
      integer      ispin,nein(2),nbrillioun
      logical      oddelcfill

#include "bafdecls.fh"
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
        value = BA_push_get(mt_int,18*(NMAX+1),
     >                    'filling',filling(2),filling(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)


        call sub_c_wvfnc_new(NMAX,int_mb(filling(1)),
     >                       wavefunction_filename,
     >                       version,
     >                       ngrid,
     >                       unita,
     >                       spin_orbit,
     >                       ispin,
     >                       ne,
     >                       nbrillioun,
     >                       fractional,
     >                       fractional_orbitals,oddelcfill)
        write(*,*) "random planewave guess, initial psi:",
     >             wavefunction_filename
        if (spin_orbit) then
           write(*,'(A,2I5)') "-   spin_orbit, nalpha=nbeta:",
     >                        ne(1),ne(2)
        else
           write(*,'(A,I2,2I5)') "-   spin, nalpha, nbeta:",
     >                           ispin,ne(1),ne(2)
        end if
        call flush(6)
        value = BA_pop_stack(filling(2))
        if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      end if
      call ga_sync()
      return
      end

*     **********************************************
*     *                                            *
*     *              sub_c_wvfnc_new               *
*     *                                            *
*     **********************************************

      subroutine sub_c_wvfnc_new(NMAX,filling,
     >                      wavefunction_filename,
     >                      version,
     >                      ngrid,
     >                      unita,
     >                      spin_orbit,
     >                      ispin,
     >                      ne,
     >                      nbrillioun,
     >                      fractional,frac_orb,oddelcfill)
      implicit none
      integer NMAX
      integer filling(4,(NMAX+1),2)
      character*50 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      logical      spin_orbit
      integer      ispin,ne(2),nbrillioun
      logical      fractional
      integer      frac_orb(2)
      logical      oddelcfill

#include "bafdecls.fh"
#include "errquit.fh"
 
      logical value
      integer ierr
      character*255 full_filename

      integer nfft1,nfft2,nfft3,nfft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l,l1,l2
      integer i,j,k,ms,n,nb,ispin0,neh,occupation

      double precision p,scale,rx,ry,rz
      double complex cc,cx,sx,zx,zc

*     **** external functions ****
      double precision gcdot,util_random
      external         gcdot,util_random
 

      p = util_random(5291999) !*** initialize the random sequence ***
      call getfilling(.false.,ne(1),filling)
      if (ispin.eq.2) call getfilling(.false.,ne(2),filling(1,1,2))

      neh = (ne(1)+mod(ne(1),2))/2

*     ***** constants *****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      nfft3d=(nfft1)*nfft2*nfft3
      inc2c = nfft1
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

      if (spin_orbit) then
         ispin0 = 3
      else
         ispin0 = ispin
      end if
      call iwrite(6,ispin0,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrillioun,1)

      if (fractional) then
         occupation = ispin0
      else
         occupation = -1
      end if
      call iwrite(6,occupation,1)
 
 
      do nb=1,nbrillioun
        call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(cfull_indx),1)

        do 150 ms=1,ispin
 
        if (ne(ms).LE.0) GO TO 150

        n = 1
        if ((spin_orbit).and.
     >      ( ((ms.eq.1).and.(n.gt.neh)).or.
     >        ((ms.eq.2).and.(n.le.neh))   )) then
          call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d) !** should be writing zeros **
        else
          l1= inc3c*filling(3,n,ms)
     >      + inc2c*filling(2,n,ms)
     >      +       filling(1,n,ms)
          if (filling(4,n,ms).lt.0) cc=sx
          if (filling(4,n,ms).eq.0) cc=zx
          if (filling(4,n,ms).gt.0) cc=cx
          dcpl_mb(cfull_indx+l1) = cc
          do k=0,nfft3d-1
            rx = (0.5d0-util_random(0))/dsqrt(dble(nfft3d))
            ry = (0.5d0-util_random(0))/dsqrt(dble(nfft3d))
            dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >                            + dcmplx(rx,ry)
          end do
          zc = dcpl_mb(cfull_indx)
          dcpl_mb(cfull_indx) = dcmplx(dble(zc),0.0d0)
          if (dble(zc).gt.1.0d0) 
     >       dcpl_mb(cfull_indx) = dcmplx(1.0d0,0.0d0)
          p=gcdot(nfft1,nfft2,nfft3,
     >            dcpl_mb(cfull_indx),
     >            dcpl_mb(cfull_indx))
          p=1.0d0/dsqrt(p)
          call dscal(2*nfft3d,p,dcpl_mb(cfull_indx),1)

          call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d)
          dcpl_mb(cfull_indx+l1)=dcmplx(0.0d0,0.0d0)
        end if
 
        DO 130 n=2,ne(ms)
          if ((spin_orbit).and.
     >        ( ((ms.eq.1).and.(n.gt.neh)).or.
     >          ((ms.eq.2).and.(n.le.neh))   )) then
            call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d) !** should be writing zeros **
          else
            l1= inc3c*filling(3,n,ms)
     >        + inc2c*filling(2,n,ms)
     >        +       filling(1,n,ms)
            l2 = inc3c*mod(2*nfft3-filling(3,n,ms),nfft3)
     >         + inc2c*mod(2*nfft2-filling(2,n,ms),nfft2)
     >         +       filling(1,n,ms)
            if (filling(4,n,ms).lt.0) cc=sx
            if (filling(4,n,ms).eq.0) cc=zx
            if (filling(4,n,ms).gt.0) cc=cx
            dcpl_mb(cfull_indx+l1)=cc
            dcpl_mb(cfull_indx+l2)=dconjg(cc)


c           IF((ABS(filling(4,n,ms)).gt.1)) THEN
              do 125 k=0,nfft3d-1
                rx = (0.5d0-util_random(0))/dsqrt(dble(nfft3d))
                ry = (0.5d0-util_random(0))/dsqrt(dble(nfft3d))
                dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >                              + dcmplx(rx,ry)
  125         continue
              zc = dcpl_mb(cfull_indx)
              dcpl_mb(cfull_indx) = dcmplx(dble(zc),0.0d0)
              if (dble(zc).gt.1.0d0) 
     >           dcpl_mb(cfull_indx) = dcmplx(1.0d0,0.0d0)
              p=gcdot(nfft1,nfft2,nfft3,
     >                 dcpl_mb(cfull_indx),
     >                 dcpl_mb(cfull_indx))
              p=1.0d0/dsqrt(p)
              call dscal(2*nfft3d,p,dcpl_mb(cfull_indx),1)
c           END IF
            call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d)
            call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(cfull_indx),1)
            dcpl_mb(cfull_indx+l1)=dcmplx(0.0d0,0.0d0)
            dcpl_mb(cfull_indx+l2)=dcmplx(0.0d0,0.0d0)
          end if
  130   CONTINUE
 

  150   CONTINUE
 
      end do

c     **** add occupation ****
      if (fractional) then
      do nb=1,nbrillioun
         rx = 1.0d0
         ry = 0.0d0
         rz = 0.5d0
         do ms=1,ispin
           do n=1,ne(ms)
             if (oddelcfill.and.(n.eq.(1+ne(ms)-frac_orb(ms)))) then
                call dwrite(6,rz,1)
             else if (n.le.(ne(ms)-frac_orb(ms))) then
                call dwrite(6,rx,1)
             else
                call dwrite(6,ry,1)
             end if
           end do
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

      double precision function gcdot(NFFT1,NFFT2,NFFT3,A,B)
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

      !inc2 = nfft1
      !inc3 = inc2*nfft2

      sum  = 0.0d0
      index = 1
      do k3=1,nfft3
      do k2=1,nfft2
      do k1=1,nfft1
         !index = (k3-1)*inc3 + (k2-1)*inc2 + k1
         sum = sum + dconjg(A(index))*B(index)
         index = index + 1
      end do
      end do
      end do

      gcdot=sum
      return
      end

 

