*
* $Id$
*

      subroutine v_wvfnc_new(v_wavefunction_filename,
     >                      version,
     >                      ngrid,
     >                      unita,
     >                      ispin,
     >                      ne)
      implicit none
      character*50 v_wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)

#include "bafdecls.fh"
 
*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value
      character*255 full_filename
      integer nfft1,nfft2,nfft3,nfft3d,n2ft3d
      integer inc2c,inc3c,occupation
      integer cfull_indx,cfull_hndl,l
      integer ms,n,ierr

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then 
      value = .false.

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

      ierr=0
      GO TO 9999
 
 9110 ierr=10
      GO TO 9999
 9111 ierr=11
      GO TO 9999
 
 9999 value = BA_free_heap(cfull_hndl)

      write(*,*) 
     >  "Generated initial wavefunctions velocities:",
     >  v_wavefunction_filename
      end if
      call ga_sync()
      
      return 
      end
 

