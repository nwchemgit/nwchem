
*     ***************************
*     *				*
*     *		xyz_init	*
*     *				*
*     ***************************
      subroutine xyz_init()
*
* $Id: xyz.F,v 1.3 2001-11-30 23:52:03 bylaska Exp $
*
      implicit none 

      integer   MASTER
      parameter (MASTER=0)

      logical found
      integer taskid
      character dum
      character*30 filename
      character*255 full_filename

*     **** external functions ****
      character*30 control_xyz
      external     control_xyz


      call Parallel_taskid(taskid)

*     **** produce XYZ FILE ****
      if (taskid.eq.MASTER) then
         filename = control_xyz()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         inquire(file=full_filename,exist=found)

         if (found) then
           open(unit=18,file=full_filename,form='formatted',
     >          status='old')
           do while (found)
           read(18,*,end=100) dum
           end do
  100      continue
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

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         close(unit=18)
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
      character*2 ion_aname
      integer     ion_nion, pspw_qmmm_nion
      real*8      ion_rion, ion_vion
      external    ion_aname
      external    ion_nion, pspw_qmmm_nion
      external    ion_rion, ion_vion

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(18,110) (ion_nion()+pspw_qmmm_nion())
         do i=1,ion_nion()
            write(18,111) ion_aname(i),
     >                  ion_rion(1,i)*AACONV,
     >                  ion_rion(2,i)*AACONV,
     >                  ion_rion(3,i)*AACONV,
     >                  ion_vion(1,i)*AACONV,
     >                  ion_vion(2,i)*AACONV,
     >                  ion_vion(3,i)*AACONV
         end do
         call pspw_qmmm_PrintXYZ(18)
         call util_flush(18)
      end if
  110 format(I5/)
  111 format(A2,5x, 6e12.4)

      return
      end



