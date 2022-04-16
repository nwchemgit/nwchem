*
* $Id: dipole_motion.F 22501 2012-05-20 06:00:35Z d3y133 $
*

*     ***************************************
*     *                                     *
*     *            dipole_motion_init       *
*     *                                     *
*     ***************************************
      subroutine dipole_motion_init(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      logical       found,found_bak
      character*50  filename
      character*255 full_filename,full_bak
      integer l1,l2
      real*8 time,dipole(3)

      real*8 time_shift
      common /dipole_motion_common/ time_shift

*     **** external functions ****
      logical  control_dipole_motion
      external control_dipole_motion


      if (control_dipole_motion()) then

        if (.not.btdb_cget(rtdb,'nwpw:dipole_motion_filename',
     >                      1,filename))
     >     call util_file_prefix('dipole_motion',filename)

        call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

        call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER)  then

*        **** check for backup file ****
         call util_file_name_noprefix('DIPOLE99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*) 
            write(*,*) "DIPOLE99-bak exists:"
            l1=index(full_bak,' ')
            l2=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:l2),
     >                 " to ",full_filename(1:l2)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if

         inquire(file=full_filename,exist=found)
         if (found) then

*         **** make a new backup file ***
          call util_file_copy(full_filename,full_bak)

          time_shift = 0.0d0
          open(unit=36,file=full_filename,form='formatted',status='old')
            do while (found)
               read(36,*,ERR=30,END=30) time_shift,dipole
            end do
  30        continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
           backspace 36
#endif
         else
          time_shift = 0.0d0
          open(unit=36,file=full_filename,form='formatted')
         end if
        end if

      end if
      return
      end


*     ***************************************
*     *                                     *
*     *            dipole_motion_end        *
*     *                                     *
*     ***************************************
      subroutine dipole_motion_end()
      implicit none

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)
      character*255 full_bak

*     **** external functions ****
      logical  control_dipole_motion
      external control_dipole_motion

      if (control_dipole_motion()) then

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
           close(unit=36)

*           **** remove backup file ***
            call util_file_name_noprefix('DIPOLE99-bak',.false.,
     >                                .false.,
     >                                full_bak)
            call util_file_unlink(full_bak)
        end if

      end if
      return
      end


*     ***************************************
*     *                                     *
*     *            dipole_motion_write      *
*     *                                     *
*     ***************************************

      subroutine dipole_motion_write(periodic,time,
     >                               ispin,ne,neq,npack1,nfft3d,dn,psi)
      implicit none
      logical periodic
      real*8  time
      integer ispin,ne(2),neq(2),npack1,nfft3d
      real*8  dn(*)
      complex*16 psi(*)

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer i
      real*8 dipole(3)

      real*8 time_shift
      common /dipole_motion_common/ time_shift

*     **** external functions ****
      logical  control_dipole_motion
      external control_dipole_motion

      if (control_dipole_motion()) then

        if (periodic) then
           call Calculate_Resta_Dipole(.false.,ispin,ne,neq,
     >                                 npack1,nfft3d,psi,
     >                                 dipole)
        else
           call Calculate_Molecular_Dipole(ispin,ne,2*nfft3d,dn,dipole)
        end if

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
           write(36,110) time+time_shift,dipole
           call util_flush(36)
        end if
  110 format(e14.6,9e14.6)

      end if

      return
      end


