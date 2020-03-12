*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   wvfnc_expander		*
*     *						*
*     *******************************************

      logical function wvfnc_expander(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin

      character*50 new_wavefunction_filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2


      integer ngrid(3)
      integer dngrid(3)
      integer cfull(2),dcfull(2),occ1(2)
      integer nfft3d,n2ft3d
      integer dnfft3d,dn2ft3d
      integer ms,n,l,occupation
      double precision unita(3,3)

      logical  control_print
      external control_print

      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'xpndr:old_wavefunction_filename',
     >                  1,old_wavefunction_filename)
      value = btdb_cget(rtdb,'xpndr:new_wavefunction_filename',
     >                  1,new_wavefunction_filename)

      value = btdb_get(rtdb,'xpndr:ngrid',mt_int,3,dngrid)

 
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
 
      call util_file_name('wvfnc_expander',
     >                    .true.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      call openfile(6,full_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,dngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,occupation,1)

       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
      dnfft3d = (dngrid(1)/2+1)*dngrid(2)*dngrid(3)
       n2ft3d = 2* nfft3d
      dn2ft3d = 2*dnfft3d
      if (control_print(print_medium)) then
      write(luout,109) old_wavefunction_filename
      write(luout,110) new_wavefunction_filename
      write(luout,111) ngrid(1), ngrid(2), ngrid(3),
     >            dngrid(1),dngrid(2),dngrid(3)
      end if
  109 format(' old_filename: ',A)
  110 format(' new_filename: ',A)
  111 format(' converting  : ',I3,'x',I3,'x',I3,' --> ', 
     >                     I3,'x',I3,'x',I3)
 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))

      value = BA_alloc_get(mt_dcpl,dnfft3d,
     >                     'dcfull',dcfull(2),dcfull(1))
      value = BA_alloc_get(mt_dbl,ne(1)+ne(2),
     >                     'occ1',occ1(2),occ1(1))

      do ms=1,ispin
        do n=1,ne(ms)
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)
          if (control_print(print_medium))
     >      write(luout,'(A,I5,A,I2)') 
     >      "converting .... psi:", n," spin:",ms
          call wvfnc_expander_convert(ngrid,dcpl_mb(cfull(1)),
     >                               dngrid,dcpl_mb(dcfull(1)))

          call dwrite(6,dcpl_mb(dcfull(1)),dn2ft3d)

        end do
      end do
      if (occupation.gt.0) then
          call dread(5,dbl_mb(occ1(1)),(ne(1)+ne(2)))
          call dwrite(6,dbl_mb(occ1(1)),(ne(1)+ne(2)))
      end if
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
      value = BA_free_heap(occ1(2))
      !IF(IERR.EQ.0) THEN
      !  WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      !ELSE
      !  WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      !  value = .false.
      !ENDIF
      !call nwpw_message(4)
      
      wvfnc_expander = value
      return 
      end
 


*     ***************************************************
*     *							*
*     *	 	   wvfnc_expander_convert		*
*     *							*
*     ***************************************************

      subroutine wvfnc_expander_convert(ngrid,psi1,dngrid,psi2)
      implicit none
      integer ngrid(3)
      complex*16 psi1(*)
      integer dngrid(3)
      complex*16 psi2(*)

*     **** local variables ****
      integer nfft3d,dnfft3d,n2ft3d,dn2ft3d
      integer inc2,inc3,dinc2,dinc3
      integer i,j,k,i2,j2,k2,n1,n2,n3
      integer jdiff,kdiff,indx,dindx
      logical jreverse,kreverse
       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
      dnfft3d = (dngrid(1)/2+1)*dngrid(2)*dngrid(3)
       n2ft3d = 2* nfft3d
      dn2ft3d = 2*dnfft3d
       inc2 = ( ngrid(1)/2+1)
      dinc2 = (dngrid(1)/2+1)
       inc3 = ( ngrid(1)/2+1)* ngrid(2)
      dinc3 = (dngrid(1)/2+1)*dngrid(2)


      n1 = ngrid(1)
      n2 = ngrid(2)
      n3 = ngrid(3)
      if (n1.gt.dngrid(1)) n1 = dngrid(1)
      if (n2.gt.dngrid(2)) n2 = dngrid(2)
      if (n3.gt.dngrid(3)) n3 = dngrid(3)

      jdiff = dngrid(2) - ngrid(2)
      kdiff = dngrid(3) - ngrid(3)
      jreverse = (jdiff.lt.0)
      kreverse = (kdiff.lt.0)
      if (jreverse) jdiff = -jdiff
      if (kreverse) kdiff = -kdiff

      call dcopy(dn2ft3d,0.0d0,0,psi2,1)
      do k=0,n3-1
      do j=0,n2-1
      do i=0,(n1/2)
         indx  = i
         dindx = i

         if (k.lt. (n3/2)) then
           k2 = k
         else
           k2 = kdiff + k
         end if

         if (j.lt. (n2/2)) then
           j2 = j
         else
           j2 = jdiff + j
         end if
         if (jreverse) then
            indx  = indx  + j2*inc2
            dindx = dindx + j *dinc2
         else
            indx  = indx  + j *inc2
            dindx = dindx + j2*dinc2
         end if
         if (kreverse) then
            indx  = indx  + k2*inc3
            dindx = dindx + k *dinc3
         else
            indx  = indx  + k *inc3
            dindx = dindx + k2*dinc3
         end if

         psi2(dindx+1) = psi1(indx+1)
      end do
      end do
      end do

      return 
      end
 

