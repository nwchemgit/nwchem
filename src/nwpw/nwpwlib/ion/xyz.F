
*     ***************************
*     *				*
*     *		xyz_init	*
*     *				*
*     ***************************
      subroutine xyz_init()
*
* $Id$
*
      implicit none 

      integer   MASTER
      parameter (MASTER=0)

      logical found,found_bak
      integer taskid,l1,l2
      character*50 filename
      character*255 full_filename,full_bak

*     **** external functions ****
      character*50 control_xyz
      external     control_xyz


      call Parallel_taskid(taskid)

*     **** produce XYZ FILE ****
      if (taskid.eq.MASTER) then

         filename = control_xyz()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)

*        **** check for backup file ***
         call util_file_name_noprefix('XYZ99-bak',.false.,
     >                                  .false.,
     >                                  full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*)
            write(*,*) "XYZ99-bak exists:"
            l1=index(full_bak,' ')
            l2=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:l2),
     >                 " to ",full_filename(1:l2)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if

         inquire(file=full_filename,exist=found)
         if (found) then

*           **** make a new backup file ***
            call util_file_copy(full_filename,full_bak)

           open(unit=18,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
             read(18,*,end=100) 
           end do
  100      continue
#if defined(FUJITSU) || defined(PSCALE) || defined(SOLARIS) || defined(__crayx1) || defined(GCC46)
           backspace 18
#endif
         else
           open(unit=18,file=full_filename,form='formatted',
     >           status='new')
         end if
      end if


      return
      end


*     ***************************
*     *				*
*     *		xyz_end 	*
*     *				*
*     ***************************
      subroutine xyz_end()
      implicit none

      integer   MASTER
      parameter (MASTER=0)

      integer taskid
      character*255 full_bak

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         close(unit=18)

*        **** remove backup file ***
         call util_file_name_noprefix('XYZ99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)
      end if

      return
      end


*     ***************************
*     *				*
*     *		xyz_write	*
*     *				*
*     ***************************
      subroutine xyz_write()
      implicit none

      real*8    AACONV
      parameter (AACONV=0.529177d0)
      integer   MASTER
      parameter (MASTER=0)

      integer taskid,i

*     **** external functions ***
      character*2 ion_aname_nocarat
      integer     ion_nion
      real*8      ion_rion, ion_vion
      external    ion_aname_nocarat
      external    ion_nion
      external    ion_rion, ion_vion
      integer ierr

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(18,110) ion_nion()
         do i=1,ion_nion()
            write(18,111) ion_aname_nocarat(i),
     >                  ion_rion(1,i)*AACONV,
     >                  ion_rion(2,i)*AACONV,
     >                  ion_rion(3,i)*AACONV,
     >                  ion_vion(1,i)*AACONV,
     >                  ion_vion(2,i)*AACONV,
     >                  ion_vion(3,i)*AACONV
         end do
         call util_flush(18)
      end if
  110 format(I5/)
  111 format(A2,5x, 6e12.4)

      return
      end

*     ***************************
*     *                         *
*     *  xyz_novelecity_write   *
*     *                         *
*     ***************************
      subroutine xyz_novelocity_write()
      implicit none

      real*8    AACONV
      parameter (AACONV=0.529177d0)
      integer   MASTER
      parameter (MASTER=0)

      integer taskid,i

*     **** external functions ***
      character*2 ion_aname_nocarat
      integer     ion_nion
      real*8      ion_rion, ion_vion
      external    ion_aname_nocarat
      external    ion_nion
      external    ion_rion, ion_vion
      integer ierr

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         write(18,110) ion_nion()
         do i=1,ion_nion()
            write(18,111) ion_aname_nocarat(i),
     >                  ion_rion(1,i)*AACONV,
     >                  ion_rion(2,i)*AACONV,
     >                  ion_rion(3,i)*AACONV
         end do
         call util_flush(18)
      end if
  110 format(I5/)
  111 format(A2,5x, 6e12.4)

      return
      end



