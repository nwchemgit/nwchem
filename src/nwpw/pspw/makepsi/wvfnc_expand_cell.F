*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   wvfnc_expand_cell		*
*     *						*
*     *******************************************

      logical function wvfnc_expand_cell(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin,dne(2)

      character*50 new_wavefunction_filename,filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2


      integer ngrid(3)
      integer dngrid(3)
      integer cell_expand(3),i,j,k
      integer cfull(2),dcfull(2)
      integer nfft3d,n2ft3d
      integer dnfft3d,dn2ft3d
      integer ms,n,l,occupation


      double precision unita(3,3)

      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'xpndr:old_wavefunction_filename',
     >                  1,old_wavefunction_filename)
      value = btdb_cget(rtdb,'xpndr:new_wavefunction_filename',
     >                  1,new_wavefunction_filename)

      value = btdb_get(rtdb,'nwpw:cell_expand',mt_int,3,cell_expand)

 
      call util_file_name_noprefix(old_wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)

      l = index(full_filename,' ') - 1
      call openfile(5,full_filename,l,'r',l)
      call iread(5,version,1)
      call iread(5,ngrid,3)
      call dread(5,unita,9)
      call iread(5,ispin,1)
      call iread(5,ne,2)
      call iread(5,occupation,1)
 
      dngrid(1) = cell_expand(1)*ngrid(1)
      dngrid(2) = cell_expand(2)*ngrid(2)
      dngrid(3) = cell_expand(3)*ngrid(3)
      dne(1) = cell_expand(1)*cell_expand(2)*cell_expand(3)*ne(1)
      dne(2) = cell_expand(1)*cell_expand(2)*cell_expand(3)*ne(2)
      unita(1,1) = unita(1,1)*cell_expand(1)
      unita(2,1) = unita(2,1)*cell_expand(1)
      unita(3,1) = unita(3,1)*cell_expand(1)
      unita(1,2) = unita(1,2)*cell_expand(2)
      unita(2,2) = unita(2,2)*cell_expand(2)
      unita(3,2) = unita(3,2)*cell_expand(2)
      unita(1,3) = unita(1,3)*cell_expand(3)
      unita(2,3) = unita(2,3)*cell_expand(3)
      unita(3,3) = unita(3,3)*cell_expand(3)

      l = index(new_wavefunction_filename,' ') - 1
      filename = new_wavefunction_filename(1:l)//".expander"
      call util_file_name(filename,
     >                    .true.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      call openfile(6,full_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,dngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,dne,2)
      call iwrite(6,occupation,1)

       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
      dnfft3d = (dngrid(1)/2+1)*dngrid(2)*dngrid(3)
       n2ft3d = 2* nfft3d
      dn2ft3d = 2*dnfft3d

      write(*,109) old_wavefunction_filename
      write(*,110) new_wavefunction_filename
      write(*,111) ngrid(1), ngrid(2), ngrid(3),
     >            dngrid(1),dngrid(2),dngrid(3)
  109 format(' old_filename: ',A)
  110 format(' new_filename: ',A)
  111 format(' converting  : ',I3,'x',I3,'x',I3,' --> ', 
     >                     I3,'x',I3,'x',I3)
 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))

      value = BA_alloc_get(mt_dcpl,dnfft3d,
     >                     'dcfull',dcfull(2),dcfull(1))

      do ms=1,ispin
        do n=1,ne(ms)
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)

          write(*,'(A,I5,A,I2)') "converting .... psi:", n," spin:",ms
          do k=0,cell_expand(3)-1
          do j=0,cell_expand(2)-1
          do i=0,cell_expand(1)-1
          write(*,*) "           ....    :", i,j,k
          call wvfnc_expand_cell_convert(i,j,k,
     >                               ngrid,dcpl_mb(cfull(1)),
     >                               dngrid,dcpl_mb(dcfull(1)))

          call dwrite(6,dcpl_mb(dcfull(1)),dn2ft3d)
          end do
          end do
          end do

        end do
      end do
      call closefile(5)
      call closefile(6)

c     *** copy wvfnc_expander to new_wavefunction_filename ****
      call util_file_name_noprefix(new_wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)
      call util_file_copy(full_filename,full_filename2)
      call util_file_unlink(full_filename)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = BA_free_heap(cfull(2))
      value = BA_free_heap(dcfull(2))
      !IF(IERR.EQ.0) THEN
      !  WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      !ELSE
      !  WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      !  value = .false.
      !ENDIF
      !call nwpw_message(4)
      
      wvfnc_expand_cell = value
      return 
      end
 

*     ***************************************************
*     *							*
*     *	 	   wvfnc_expand_cell_convert   		*
*     *							*
*     ***************************************************

      subroutine wvfnc_expand_cell_convert(ishift,jshift,kshift,
     >                              ngrid,psi1,dngrid,psi2)
      implicit none
      integer ishift,jshift,kshift
      integer ngrid(3)
      complex*16 psi1(*)
      integer dngrid(3)
      complex*16 psi2(*)

*     **** local variables ****
      integer nfft3d,dnfft3d,n2ft3d,dn2ft3d
      integer inc2,inc3,dinc2,dinc3
      integer nxh,nyh,nzh
      integer i,j,k
      integer i1,j1,k1
      integer i2,j2,k2
      integer fi,fj,fk
      integer indx,dindx
      real*8  pi,gi,gj,gk
      complex*16 bw1,bw2,bw3
      complex*16 cw1,cw2,cw3
       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
      dnfft3d = (dngrid(1)/2+1)*dngrid(2)*dngrid(3)
       n2ft3d = 2* nfft3d
      dn2ft3d = 2*dnfft3d
       inc2 = ( ngrid(1)/2+1)
      dinc2 = (dngrid(1)/2+1)
       inc3 = ( ngrid(1)/2+1)* ngrid(2)
      dinc3 = (dngrid(1)/2+1)*dngrid(2)

      nxh = ngrid(1)/2
      nyh = ngrid(2)/2
      nzh = ngrid(3)/2
      fi  = dngrid(1)/ngrid(1)
      fj  = dngrid(2)/ngrid(2)
      fk  = dngrid(3)/ngrid(3)
      gi  = 1.0d0/dble(fi)
      gj  = 1.0d0/dble(fj)
      gk  = 1.0d0/dble(fk)

      pi = 4.0d0*datan(1.0d0)
      cw1=dcmplx(dcos(0.5d0*ishift*pi),-dsin(0.5d0*ishift*pi))
      cw2=dcmplx(dcos(0.5d0*jshift*pi),-dsin(0.5d0*jshift*pi))
      cw3=dcmplx(dcos(0.5d0*kshift*pi),-dsin(0.5d0*kshift*pi))
c     write(*,*) "cws:",cw1,cw2,cw3
      cw1=dcmplx(dcos(gi*pi*ishift),-dsin(gi*pi*ishift))
      cw2=dcmplx(dcos(gj*pi*jshift),-dsin(gj*pi*jshift))
      cw3=dcmplx(dcos(gk*pi*kshift),-dsin(gk*pi*kshift))

      call dcopy(dn2ft3d,0.0d0,0,psi2,1)
      do k=-nzh+1,nzh-1
      do j=-nyh+1,nyh-1
      do i=0,nxh-1

         i1 = i
         j1 = j
         k1 = k
         !if (i1.ge.0) i2 = fi*i1 + ishift
         !if (i1.lt.0) i2 = fi*i1 - ishift
         !if (j1.ge.0) j2 = fj*j1 + jshift
         !if (j1.lt.0) j2 = fj*j1 - jshift
         !if (k1.ge.0) k2 = fk*k1 + kshift
         !if (k1.lt.0) k2 = fk*k1 - kshift
         i2 = fi*i1 + ishift
         j2 = fj*j1 + jshift
         k2 = fk*k1 + kshift
         



         if (i1 .lt. 0) i1 = i1 + ngrid(1)
         if (j1 .lt. 0) j1 = j1 + ngrid(2)
         if (k1 .lt. 0) k1 = k1 + ngrid(3)

         if (i2 .lt. 0) i2 = i2 + dngrid(1)
         if (j2 .lt. 0) j2 = j2 + dngrid(2)
         if (k2 .lt. 0) k2 = k2 + dngrid(3)

         indx   = (k1)*inc3  +(j1)*inc2  + i1+1
         dindx  = (k2)*dinc3 +(j2)*dinc2 + i2+1

       
         psi2(dindx) = psi1(indx)*cw1*cw2*cw3


c        i1 = i
c        j1 = -j
c        k1 = -k
c        i2 = fi*i1
c        j2 = fj*j1
c        k2 = fk*k1

c        if (i1 .lt. 0) i1 = i1 + ngrid(1)
c        if (j1 .lt. 0) j1 = j1 + ngrid(2)
c        if (k1 .lt. 0) k1 = k1 + ngrid(3)

c        if (i2 .lt. 0) i2 = i2 + dngrid(1)
c        if (j2 .lt. 0) j2 = j2 + dngrid(2)
c        if (k2 .lt. 0) k2 = k2 + dngrid(3)

c        indx   = (k1)*inc3  +(j1)*inc2  + i1+1
c        dindx  = (k2)*dinc3 +(j2)*dinc2 + i2+1

c        psi2(dindx) = psi1(indx)


      end do
      end do
      end do

      return 
      end
 


