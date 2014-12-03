*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   wvfnc_convert_restricted     *
*     *						*
*     *******************************************

      logical function wvfnc_convert_restricted(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin,nne(2),iispin

      character*50 new_wavefunction_filename,filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2,full_filename3


      integer ngrid(3)
      integer dngrid(3)
      integer cfull(2),dcfull(2)
      integer nfft3d,n2ft3d
      integer dnfft3d,dn2ft3d
      integer ms,n,l,ma_type,occupation


      double precision unita(3,3)


      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'ionize:old_wavefunction_filename',
     >                  1,old_wavefunction_filename)
      value = btdb_cget(rtdb,'ionize:new_wavefunction_filename',
     >                  1,new_wavefunction_filename)

 
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
 
      iispin = 2
      nne(1) = ne(1)
      nne(2) = ne(1)
      l = index(new_wavefunction_filename,' ') - 1
      filename = new_wavefunction_filename(1:l)//".xpander"
      call util_file_name(filename,
     >                    .true.,
     >                    .false.,
     >                    full_filename2)
      l = index(full_filename2,' ') - 1
      call openfile(6,full_filename2,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,iispin,1)
      call iwrite(6,nne,2)
      call iwrite(6,occupation,1)

       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
       n2ft3d = 2* nfft3d

      write(*,109) old_wavefunction_filename
      write(*,110) new_wavefunction_filename
      write(*,111) 
  109 format(' old_filename: ',A)
  110 format(' new_filename: ',A)
  111 format(' converting  :  restricted  --> unrestricted')
 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))


      do ms=1,1
        do n=1,ne(ms)
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)

          write(*,*) "copying .... psi:", n," spin:",1
          call dwrite(6,dcpl_mb(cfull(1)),n2ft3d)

        end do
      end do
      call closefile(5)


      l = index(full_filename,' ') - 1
      call openfile(5,full_filename,l,'r',l)
      call iread(5,version,1)
      call iread(5,ngrid,3)
      call dread(5,unita,9)
      call iread(5,ispin,1)
      call iread(5,ne,2)
      do ms=1,1
        do n=1,ne(ms)
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)

          write(*,*) "copying .... psi:", n," spin:",2
          call dwrite(6,dcpl_mb(cfull(1)),n2ft3d)

        end do
      end do
      call closefile(5)
      call closefile(6)

c     *** copy wvfnc_expander to new_wavefunction_filename ****
      call util_file_name_noprefix(new_wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename3)
      call util_file_copy(full_filename2,full_filename3)
      call util_file_unlink(full_filename2)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = BA_free_heap(cfull(2))
      !IF(IERR.EQ.0) THEN
      !  WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      !ELSE
      !  WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      !  value = .false.
      !ENDIF
      !call nwpw_message(4)
      
      wvfnc_convert_restricted = value
      return 
      end
