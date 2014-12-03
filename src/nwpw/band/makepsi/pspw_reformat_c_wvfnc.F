*
* $Id$
*

      logical function  pspw_reformat_c_wvfnc(nbrill)
      implicit none
      integer nbrill

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
 
*     **** local variables ****
      logical value,oprint,hprint
      character*50 wavefunction_filename,filename
      character*255 full_filename,full_filename2
      character*255 full_filename3
      integer version,ierr
      integer ne(2),ispin,neall
      integer ngrid(3)
      integer cfull(2),cfull_out(2)
      integer nfft3d,n2ft3d,nfft3d_out
      integer n,l,nb,occupation

      double precision unita(3,3)

*     ***** local functions ****
      logical  control_print
      external control_print

      character*50 control_input_psi
      external     control_input_psi


      value = .false.
      hprint = control_print(print_high)
      oprint = control_print(print_medium)



      version = 3  ! or 4


      wavefunction_filename = control_input_psi() 
      call util_file_name_noprefix(wavefunction_filename,
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
 
      version = 5
      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".junk"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)

      l = index(full_filename2,' ') - 1
      call openfile(6,full_filename2,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrill,1)
      call iwrite(6,occupation,1)
       
       nfft3d = ( ngrid(1)/2+1)*ngrid(2)*ngrid(3)
       n2ft3d = 2* nfft3d
       nfft3d_out =  ngrid(1)*ngrid(2)*ngrid(3)

      if (oprint) then
      write(*,*)
      write(*,*) 
     > '+ converting wavefunctions from pspw format to band format'

      end if

      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".pspw"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename3)

      if (oprint) then
      write(*,109) full_filename
      write(*,110) full_filename3
      end if
  109 format(' filename converted    : ',A)
  110 format(' pspw filename moved to: ',A)

      neall = ne(1) + ne(2) 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d_out,
     >                     'cfull_out',cfull_out(2),cfull_out(1))
      if (.not.value) 
     > call errquit("pspw_reformat_c_wvfnc:push stack",0, MA_ERR)

       
      do nb=1,nbrill
        do n=1,neall

          if (hprint) then
          write(*,111) "converting ... psi:",n," nb:",nb
          end if
  111     format(A,I4,A,I2,A,I4)

          call  dread(5,dcpl_mb(cfull(1)),n2ft3d)
          call  convert_real_to_complex(ngrid(1),ngrid(2),ngrid(3),
     >                                  dcpl_mb(cfull(1)),
     >                                  dcpl_mb(cfull_out(1)))
          call dwrite(6,dcpl_mb(cfull_out(1)),2*nfft3d_out)
        end do
      end do
      call closefile(5)
      call closefile(6)
      value = BA_free_heap(cfull(2))
      value = BA_free_heap(cfull_out(2))

c     *** copy .junk wvfncr to wavefunction_filename ****
      call util_file_copy(full_filename,full_filename3)
      call util_file_copy(full_filename2,full_filename)
      call util_file_unlink(full_filename2)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 IF(IERR.ne.0) THEN
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
        value = .false.
      ENDIF
      
      pspw_reformat_c_wvfnc = value
      return 
      end



      logical function  pspw_reformat_c_ewvfnc(nbrill)
      implicit none
      integer nbrill

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
 
*     **** local variables ****
      logical value,oprint,hprint
      character*50 wavefunction_filename,filename
      character*255 full_filename,full_filename2
      character*255 full_filename3
      integer version,ierr
      integer ne(2),ispin,neall
      integer ngrid(3)
      integer cfull(2),cfull_out(2)
      integer nfft3d,n2ft3d,nfft3d_out
      integer n,l,nb,occupation

      double precision unita(3,3)

*     ***** local functions ****
      logical  control_print
      external control_print

      character*50 control_input_epsi
      external     control_input_epsi

      value = .false.
      hprint = control_print(print_high)
      oprint = control_print(print_medium)

      version = 3  ! or 4

      wavefunction_filename = control_input_epsi() 
      call util_file_name_noprefix(wavefunction_filename,
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
 
      version = 5
      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".junk"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)

      l = index(full_filename2,' ') - 1
      call openfile(6,full_filename2,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrill,1)
      call iwrite(6,occupation,1)
       
       nfft3d = ( ngrid(1)/2+1)*ngrid(2)*ngrid(3)
       n2ft3d = 2* nfft3d
       nfft3d_out =  ngrid(1)*ngrid(2)*ngrid(3)

      if (oprint) then
      write(*,*)
      write(*,*) 
     > '+ converting wavefunctions from pspw format to band format'

      end if

      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".pspw"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename3)

      if (oprint) then
      write(*,109) full_filename
      write(*,110) full_filename3
      end if
  109 format(' filename converted    : ',A)
  110 format(' pspw filename moved to: ',A)

      neall = ne(1) + ne(2) 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d_out,
     >                     'cfull_out',cfull_out(2),cfull_out(1))
      if (.not.value) 
     > call errquit("pspw_reformat_c_ewvfnc:push stack",0, MA_ERR)

       
      do nb=1,nbrill
        do n=1,neall

          if (hprint) then
          write(*,111) "converting ... psi:",n," nb:",nb
          end if
  111     format(A,I4,A,I2,A,I4)

          call  dread(5,dcpl_mb(cfull(1)),n2ft3d)
          call  convert_real_to_complex(ngrid(1),ngrid(2),ngrid(3),
     >                                  dcpl_mb(cfull(1)),
     >                                  dcpl_mb(cfull_out(1)))
          call dwrite(6,dcpl_mb(cfull_out(1)),2*nfft3d_out)
        end do
      end do
      call closefile(5)
      call closefile(6)
      value = BA_free_heap(cfull(2))
      value = BA_free_heap(cfull_out(2))

c     *** copy .junk wvfncr to wavefunction_filename ****
      call util_file_copy(full_filename,full_filename3)
      call util_file_copy(full_filename2,full_filename)
      call util_file_unlink(full_filename2)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 IF(IERR.ne.0) THEN
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
        value = .false.
      ENDIF
      
      pspw_reformat_c_ewvfnc = value
      return 
      end



      logical function  pspw_reformat_c_v_wvfnc(nbrill)
      implicit none
      integer nbrill

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"

*     **** local variables ****
      logical value,oprint,hprint
      character*50 wavefunction_filename,filename
      character*255 full_filename,full_filename2
      character*255 full_filename3
      integer version,ierr
      integer ne(2),ispin,neall
      integer ngrid(3)
      integer cfull(2),cfull_out(2)
      integer nfft3d,n2ft3d,nfft3d_out
      integer n,l,nb,occupation

      double precision unita(3,3)

*     ***** local functions ****
      logical  control_print
      external control_print

      character*50 control_input_v_psi
      external     control_input_v_psi

      value = .false.
      hprint = control_print(print_high)
      oprint = control_print(print_medium)

      version = 3  ! or 4


      wavefunction_filename = control_input_v_psi()
      call util_file_name_noprefix(wavefunction_filename,
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

      version = 5
      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".junk"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)

      l = index(full_filename2,' ') - 1
      call openfile(6,full_filename2,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,ne,2)
      call iwrite(6,nbrill,1)
      call iwrite(6,occupation,1)

       nfft3d = ( ngrid(1)/2+1)*ngrid(2)*ngrid(3)
       n2ft3d = 2* nfft3d
       nfft3d_out =  ngrid(1)*ngrid(2)*ngrid(3)


      if (oprint) then
      write(*,*)
      write(*,*)
     > '+ converting wavefunctions from pspw format to band format'

      end if

      l = index(wavefunction_filename,' ') - 1
      filename = wavefunction_filename(1:l)//".pspw"
      call util_file_name_noprefix(filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename3)

      if (oprint) then
      write(*,109) wavefunction_filename
      write(*,110) filename
      end if
  109 format(' filename converted    : ',A)
  110 format(' pspw filename moved to: ',A)

      neall = ne(1) + ne(2)
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,nfft3d_out,
     >                     'cfull_out',cfull_out(2),cfull_out(1))
      if (.not.value)
     > call errquit("pspw_reformat_c_ewvfnc:push stack",0, MA_ERR)

      do nb=1,nbrill
        do n=1,neall

          if (hprint) then
          write(*,111) "converting ... psi:",n," nb:",nb
          end if
  111     format(A,I4,A,I2,A,I4)

          call  dread(5,dcpl_mb(cfull(1)),n2ft3d)
          call  convert_real_to_complex(ngrid(1),ngrid(2),ngrid(3),
     >                                  dcpl_mb(cfull(1)),
     >                                  dcpl_mb(cfull_out(1)))
          call dwrite(6,dcpl_mb(cfull_out(1)),2*nfft3d_out)
        end do
      end do
      call closefile(5)
      call closefile(6)
      value = BA_free_heap(cfull(2))
      value = BA_free_heap(cfull_out(2))

c     *** copy . junk wvfncr to wavefunction_filename ****
      call util_file_copy(full_filename,full_filename3)
      call util_file_copy(full_filename2,full_filename)
      call util_file_unlink(full_filename2)
      IERR=0
      GO TO 9999

 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999

 9999 IF(IERR.ne.0) THEN
        WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
        value = .false.
      ENDIF

      pspw_reformat_c_v_wvfnc = value
      return
      end








      subroutine convert_real_to_complex(nx,ny,nz,A,B)
      implicit none
      integer nx,ny,nz
      complex*16 A(nx/2+1,ny,nz)
      complex*16 B(nx,ny,nz)

*     **** local variables ****
      integer nxh,nyh,nzh
      integer k1,k2,k3
      integer i,j,k
      integer ii,jj,kk

      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      call dcopy(2*nx*ny*nz,0.0d0,0,B,1)

      do k3 = -nzh+1,nzh-1
      do k2 = -nyh+1,nyh-1
      do k1 = 0,nxh-1
         i = k1
         j = k2
         k = k3
         if (i .lt. 0) i = i + nx
         if (j .lt. 0) j = j + ny
         if (k .lt. 0) k = k + nz
         B(i+1,j+1,k+1) = A(i+1,j+1,k+1)

         if (k1.ne.0) then
           ii = -k1
           jj = -k2
           kk = -k3
           if (ii .lt. 0) ii = ii + nx
           if (jj .lt. 0) jj = jj + ny
           if (kk .lt. 0) kk = kk + nz
           B(ii+1,jj+1,kk+1) = dconjg(A(i+1,j+1,k+1))
         end if
      end do
      end do
      end do
      return
      end 


 
 

