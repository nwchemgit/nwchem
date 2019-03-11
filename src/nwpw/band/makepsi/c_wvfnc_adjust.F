*
* $Id$
*

*     ***********************************************
*     *                                             *
*     *              c_wvfnc_adjust                 *
*     *                                             *
*     ***********************************************

      subroutine c_wvfnc_adjust(wavefunction_filename,
     >                         ispin,nein,nbrillioun)
      implicit none
      character*50 wavefunction_filename
      integer      ispin,nein(2),nbrillioun

#include "bafdecls.fh"
#include "errquit.fh"
 
*     **** local variables ****
      logical value,fractional
      integer MASTER,taskid
      parameter (MASTER=0)

      integer NMAX
      integer filling(2)
      integer fractional_orbitals(2),ne(2)
      character*255 new_filename,old_filename,emo_filename

*     **** external functions ****
      logical  control_fractional
      integer  control_fractional_orbitals
      external control_fractional
      external control_fractional_orbitals
      character*50 control_input_epsi
      external     control_input_epsi

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
     >    call errquit('wvfnc_adjust:out of stack memory',0,MA_ERR)

        call util_file_name_noprefix('wvfnc_adjust',
     >                    .false.,
     >                    .false.,
     >                     old_filename)
        call util_file_name_noprefix(wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    new_filename)
        call util_file_copy(new_filename,old_filename)

        !*** fetch the emovecs filename ***
        call util_file_name_noprefix(control_input_epsi(),
     >                    .false.,
     >                    .false.,
     >                    emo_filename)

        !*** adjust wavefunctions ***
        call sub_c_wvfnc_adjust(NMAX,int_mb(filling(1)),
     >                     new_filename,
     >                     old_filename,
     >                     emo_filename,
     >                     ispin,
     >                     ne,
     >                     nbrillioun,
     >                     fractional,
     >                     fractional_orbitals)

        !*** remove temporary wvfnc_adjust file ***
        call util_file_unlink(old_filename)


        write(*,*) "wavefunction adjust, new psi:",
     >             wavefunction_filename
        write(*,*) "-   spin, nalpha, nbeta, nbrillioun:",
     >             ispin,ne,nbrillioun
        value = BA_pop_stack(filling(2))
        if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      end if
      call ga_sync()

      return
      end


      subroutine sub_c_wvfnc_adjust(NMAX,filling,
     >                         new_filename,
     >                         old_filename,
     >                         emo_filename,
     >                         ispin,
     >                         ne,
     >                         nbrillioun,
     >                         fractional,
     >                         frac_orb)
      implicit none
      integer NMAX
      integer filling(4,NMAX,2)
      character*255 new_filename
      character*255 old_filename
      character*255 emo_filename
      integer      ispin,ne(2),nbrillioun
      logical      fractional
      integer      frac_orb(2)

#include "bafdecls.fh"
#include "errquit.fh"
 
      logical value,emo_found,emo_used
      character*255 full_filename

      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)

      integer nfft1,nfft2,nfft3,nfft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l,l1,l2
      integer i,j,k,ms,n,occupation,ispint,nb
      integer n0,ms0,n0max,ispin0,ne0(2),nbrillioun0
      integer nx,msx,nxmax,ispinx,nex(2),nbrilliounx

      double precision p,scale
      double complex cc,cx,sx,zx,zc,rx,ry

*     **** external functions ****
      double precision gcdot,util_random
      external         gcdot,util_random
 

      p = util_random(5291999) !*** initialize the random sequence ****

      call getfilling(.true.,ne(1),filling)
      if (ispin.eq.2) call getfilling(.true.,ne(2),filling(1,1,2))


      scale=1.0d0/dsqrt(2.0d0)
      zx=(1.0d0,0.0d0)
      sx=(0.0d0,1.0d0)*scale
      cx=(1.0d0,0.0d0)*scale

 
*     **** write wavefunction in CPMDV3 format ****
      l = index(old_filename,' ') - 1
      call openfile(5,old_filename,l,'r',l)
      call iread(5,version,1)
      call iread(5,ngrid,3)
      call dread(5,unita,9)
      call iread(5,ispin0,1)
      call iread(5,ne0,2)
      call iread(5,occupation,1)
      call iread(5,nbrillioun0,1)

      l = index(new_filename,' ') - 1
      call openfile(6,new_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrillioun,1)
      if (fractional) then
         occupation = ispin
      else
         occupation = -1
      end if
      call iwrite(6,occupation,1)


*     ***** constants *****
      nfft1=ngrid(1)
      nfft2=ngrid(2)
      nfft3=ngrid(3)
      nfft3d=nfft1*nfft2*nfft3
      inc2c = nfft1
      inc3c =inc2c*nfft2

*     ***** allocate wavefunction memory ****
      value = BA_push_get(mt_dcpl,2*nfft3d,
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

*     **** try to read emo_filename ***
      nex(1) = 0
      nex(2) = 0
      emo_found = .false.
      emo_used  = .false.
      inquire(file=emo_filename,exist=emo_found)
      if (emo_found) then
         l = index(emo_filename,' ') - 1
         call openfile(3,emo_filename,l,'r',l)
         call iread(3,version,1)
         call iread(3,ngrid,3)
         call dread(3,unita,9)
         call iread(3,ispinx,1)
         call iread(3,nex,2)
         call iread(3,nbrilliounx,1)
         call iread(3,occupation,1)
      end if

      ms0 = 1
      msx = 1
      do nb=1,nbrillioun
      do ms=1,ispin
         n0 = 1
         nx = 1
         n0max = ne0(ms0)
         nxmax = nex(ms0)
         do n=1,ne(ms)

            !*** read from old filename ***
            if (n.le.n0max) then
               call dread(5,dcpl_mb(cfull_indx),2*nfft3d)
               n0 = n0 + 1

            !*** read from emo_filename ***
            else if (n.le.(n0max+nxmax)) then
               call dread(3,dcpl_mb(cfull_indx),2*nfft3d)
               nx = nx + 1
               emo_used  = .true.
               call sub_c_wvfnc_project_out(nfft3d,
     >                                  dcpl_mb(cfull_indx),
     >                                  dcpl_mb(cfull_indx+nfft3d),
     >                                  old_filename)
               P=gcdot(nfft1,nfft2,nfft3,
     >                 dcpl_mb(cfull_indx),
     >                 dcpl_mb(cfull_indx))
               P=1.0d0/dsqrt(P)
               call dscal(2*nfft3d,P,dcpl_mb(cfull_indx),1)

            !*** generate new random wavefunction ***
            else
               call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(cfull_indx),1)
               l1= inc3c*filling(3,n,ms)
     >           + inc2c*filling(2,n,ms)
     >           +       filling(1,n,ms)
               if (filling(4,n,ms).lt.0) cc=sx
               if (filling(4,n,ms).eq.0) cc=zx
               if (filling(4,n,ms).gt.0) cc=cx
                l2=l1
                dcpl_mb(cfull_indx+l1)=cc
                if (filling(1,n,ms).eq.0) then
                  l2 = inc3c*mod(nfft3-filling(3,n,ms),nfft3)
     >               + inc2c*mod(nfft2-filling(2,n,ms),nfft2)
     >               +       filling(1,n,ms)
                  dcpl_mb(cfull_indx+l2)=dconjg(cc)
                end if
c                if((ABS(filling(4,n,ms)).gt.1)) then
                  do 125 k=0,nfft3d-1
                    dcpl_mb(cfull_indx+k) = dcpl_mb(cfull_indx+k)
     >                 + dcmplx((0.5d0-util_random(0)),
     >                          (0.5d0-util_random(0)))
     >                    /dsqrt(dble(nfft3d))
  125             continue
                  zc = dcpl_mb(cfull_indx)
                  dcpl_mb(cfull_indx) = dcmplx(dble(zc),0.0d0)
                  if (dble(zc).gt.1.0d0)
     >               dcpl_mb(cfull_indx) = dcmplx(1.0d0,0.0d0)
                  call sub_c_wvfnc_project_out(nfft3d,
     >                                  dcpl_mb(cfull_indx),
     >                                  dcpl_mb(cfull_indx+nfft3d),
     >                                  old_filename)

                  P=gcdot(nfft1,nfft2,nfft3,
     >                     dcpl_mb(cfull_indx),
     >                     dcpl_mb(cfull_indx))
                  P=1.0d0/dsqrt(P)
                  call dscal(2*nfft3d,P,dcpl_mb(cfull_indx),1)
c                end if
            end if

            call dwrite(6,dcpl_mb(cfull_indx),2*nfft3d)
            !n0 = n0 + 1
         end do
         
         ms0 = ms0 + 1
         msx = msx + 1

*        **** rewind the wavefunction read ****
         if (ms0.gt.ispin0) then
            ms0 = 1
            call closefile(5)
            l = index(old_filename,' ') - 1
            call openfile(5,old_filename,l,'r',l)
            call iread(5,version,1)
            call iread(5,ngrid,3)
            call dread(5,unita,9)
            call iread(5,ispin0,1)
            call iread(5,ne0,2)
            call iread(5,nbrillioun0,2)
            call iread(5,occupation,1)
         end if

*        **** read remaining wvfunctions in spin ****
         if ((ispin0.eq.2).and.(ispin.eq.2).and.(n0.le.n0max)) then
           do n=n0,n0max
             call dread(5,dcpl_mb(cfull_indx),2*nfft3d)
           end do
         end if

*        **** read remaining emo_filename wavefunctions in spin ***
         if (emo_found) then
            if (msx.gt.ispinx) then
               msx = 1
               call closefile(3)
               l = index(emo_filename,' ') - 1
               call openfile(3,emo_filename,l,'r',l)
               call iread(3,version,1)
               call iread(3,ngrid,3)
               call dread(3,unita,9)
               call iread(3,ispinx,1)
               call iread(3,nex,2)
               call iread(3,nbrilliounx,1)
               call iread(3,occupation,1)
            end if

*           **** read remaining wvfunctions in spin ****
            if ((ispinx.eq.2).and.(ispin.eq.2).and.(nx.le.nxmax)) then
              do n=nx,nxmax
                call dread(3,dcpl_mb(cfull_indx),2*nfft3d)
              end do
            end if
         end if

      end do
      end do

c     **** add occupation - don't use previous start from scratch ****
      if (fractional) then
      do nb=1,nbrillioun
         rx = 1.0d0
         ry = 0.0d0
         do ms=1,ispin
           do n=1,ne(ms)
             if (n.le.ne0(ms)) then
                call dwrite(6,rx,1)
             else
                call dwrite(6,ry,1)
             end if
           end do
         end do
      end do
      end if
 
      if (emo_found) then
         call closefile(3)
         if (emo_used) call util_file_unlink(emo_filename) !*** remove emo_filename ***
      end if
      call closefile(5)
      call closefile(6)
 
      value = BA_pop_stack(cfull_hndl)
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)
      
      return 
      end
 
      subroutine sub_c_wvfnc_project_out(nfft3d,psi,tpsi,old_filename)
      implicit none
      integer nfft3d
      complex*16 psi(*)
      complex*16 tpsi(*)
      character*255 old_filename

      integer l,version,ngrid(3)
      integer i,occupation
      integer ispin0,ne0(2),nbrillioun0

      double precision p,unita(9)

*     **** external functions ****
      double precision gcdot
      external         gcdot

      l = index(old_filename,' ') - 1
      call openfile(4,old_filename,l,'r',l)
      call iread(4,version,1)
      call iread(4,ngrid,3)
      call dread(4,unita,9)
      call iread(4,ispin0,1)
      call iread(4,ne0,2)
      call iread(4,nbrillioun0,1)
      call iread(4,occupation,1)
      do i=1,ne0(1)
         call dread(4,tpsi,2*nfft3d)
         p = gcdot(ngrid(1),ngrid(2),ngrid(3),tpsi,psi)
         call daxpy(2*nfft3d,-p,tpsi,1,psi,1)
      end do
      call closefile(4)

      return
      end

