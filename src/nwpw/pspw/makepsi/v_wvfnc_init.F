*
* $Id$
*

      logical function v_wvfnc_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

      logical value
      integer version,occupation
 
      integer ierr

      integer NMAX
      parameter (NMAX=100)
      integer ne(2),ispin
      logical  restricted

      character*50 cell_name,v_wavefunction_filename
      character*50 rtdb_unita,rtdb_ngrid
      character*255 full_filename


      integer ngrid(3)
      integer nfft1,nfft2,nfft3,nfft3d,n2ft3d
      integer inc2c,inc3c
      integer cfull_indx,cfull_hndl,l
      integer ms,n

      double precision unita(3,3)

 
      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'v_wvfnc:v_wavefunction_filename',
     >                  1,v_wavefunction_filename)

      value = btdb_get(rtdb,'v_wvfnc:nup',  mt_int,1,ne(1))
      value = btdb_get(rtdb,'v_wvfnc:ndown',mt_int,1,ne(2))

      value = btdb_get(rtdb,'v_wvfnc:restricted', mt_log,1,restricted)
      ispin=2
      if (restricted) ispin=1


*     **** get cell information ****
      value = btdb_cget(rtdb,'v_wvfnc:cell_name',
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


*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull_hndl,cfull_indx)


      

 
*:::::::::::::::::  write wavefunction in CPMDV3 format ::::::::::::::
      call util_file_name_noprefix(v_wavefunction_filename,
     >                    .false.,
     >                    .false.,
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

        call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)

        do N=2,NE(ms)
          call dwrite(6,dcpl_mb(cfull_indx),n2ft3d)
        end do  
 

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
      
      v_wvfnc_init = value
      return 
      end
 

