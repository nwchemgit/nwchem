*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   chi_expander 		*
*     *						*
*     *******************************************

      logical function chi_expander(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin

      character*50 new_wavefunction_filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2


      integer ngrid(3)
      integer dngrid(3)
      integer nfft3d,n2ft3d
      integer dnfft3d,dn2ft3d
      integer ms,n,l
      integer cfull(2),dcfull(2)


      double precision unita(3,3)

      value = .false.
      version = 9

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
      value = BA_push_get(mt_dcpl,nfft3d,
     >                    'cfull',cfull(2),cfull(1))

      value = BA_push_get(mt_dcpl,dnfft3d,
     >                    'dcfull',dcfull(2),dcfull(1))

      do ms=1,ispin
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)

          write(*,*) "converting .... spin:",ms
          call wvfnc_expander_convert(ngrid,dcpl_mb(cfull(1)),
     >                               dngrid,dcpl_mb(dcfull(1)))

          call dwrite(6,dcpl_mb(dcfull(1)),dn2ft3d)

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
 
 9999 value =           BA_pop_stack(dcfull(2))
      value = value.and.BA_pop_stack(cfull(2))
      !IF(IERR.EQ.0) THEN
      !  WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      !ELSE
      !  WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      !  value = .false.
      !ENDIF
      !call nwpw_message(4)
      
      chi_expander = value
      return 
      end

