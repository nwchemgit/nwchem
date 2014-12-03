*
* $Id$
*

      logical function wvfnc_reformat(rtdb)
      implicit none
      integer rtdb
#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin

      character*50 new_wavefunction_filename,filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2
      character*50 cell_name,rtdb_unita

      integer ngrid(3)
      integer cfull(2)
      integer nfft3d,n2ft3d
      integer ms,n,l,occupation

      double precision unita(3,3)
      double precision unita2(3,3)

      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'cgsd:output_wavefunction_filename',
     >                  1,old_wavefunction_filename)
      value = value.and.
     >        btdb_cget(rtdb,'cgsd:input_wavefunction_filename',
     >                  1,new_wavefunction_filename)

      if (.not. value)
     >   call errquit('wvfnc_reformat: error reading filenames ',0,
     &       RTDB_ERR)

*     **** cell: stuff ****
      value = value.and.
     >        btdb_cget(rtdb,'cgsd:cell_name',1,cell_name)
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      value = value.and.
     >        btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita2)

      if (.not. value)
     >   call errquit('wvfnc_reformat: error reading cell ',0, RTDB_ERR)


 
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
 
      l = index(new_wavefunction_filename,' ') - 1
      filename = new_wavefunction_filename(1:l)//".reform"
      call util_file_name(filename,
     >                    .true.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      call openfile(6,full_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita2,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,occupation,1)
       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
       n2ft3d = 2* nfft3d

      write(*,*)
      write(*,*) '  -------------------------  '
      write(*,*) ' | adjusting wavefunctions | '
      write(*,*) '  -------------------------  '
      call nwpw_message(1)
      write(*,*)

      write(*,109) old_wavefunction_filename
      write(*,110) new_wavefunction_filename
      write(*,111) 
      write(*,112)  unita(1,1), unita(2,1), unita(3,1),
     >             unita2(1,1),unita2(2,1),unita2(3,1)
      write(*,113)  unita(1,2), unita(2,2), unita(3,2),
     >             unita2(1,2),unita2(2,2),unita2(3,2)
      write(*,112)  unita(1,3), unita(2,3), unita(3,3),
     >             unita2(1,3),unita2(2,3),unita2(3,3)
  109 format(' old filename: ',A)
  110 format(' new filename: ',A)

  111 format(2x,'unita:',28x,'unita2:')
  112 format(3F7.3, 13x, 3F7.3)
  113 format(3F7.3, 5x,"-->",5x, 3F7.3)
 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))


      do ms=1,ispin
        do n=1,ne(ms)
          write(*,*) "copying .... psi:", n," spin:",ms
          call  dread(5,dcpl_mb(cfull(1)),n2ft3d)
          call dwrite(6,dcpl_mb(cfull(1)),n2ft3d)
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
      IF(IERR.EQ.0) THEN
        WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      ELSE
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
        value = .false.
      ENDIF
      call nwpw_message(4)
      write(*,*)
      
      wvfnc_reformat = value
      return 
      end
 

