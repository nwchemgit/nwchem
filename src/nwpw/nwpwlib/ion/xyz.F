
*     ***************************
*     *				*
*     *		xyz_init	*
*     *				*
*     ***************************
      subroutine xyz_init()
*
* $Id: xyz.F,v 1.1 2001-08-30 17:56:11 bylaska Exp $
*
      implicit none 

      integer   MASTER
      parameter (MASTER=0)

      integer taskid,l
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
         l = index(filename,' ') - 1
         open(unit=18,file=full_filename,form='formatted')
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
      character*2 psp_atom
      integer     ion_nion, ion_katm,Waterpsp_nwater
      real*8      ion_rion, ion_vion
      external    psp_atom
      external    ion_nion, ion_katm,Waterpsp_nwater
      external    ion_rion, ion_vion

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(18,110) (ion_nion()+3*Waterpsp_nwater())
         do i=1,ion_nion()
            write(18,111) psp_atom(ion_katm(i)),
     >                  ion_rion(1,i)*AACONV,
     >                  ion_rion(2,i)*AACONV,
     >                  ion_rion(3,i)*AACONV,
     >                  ion_vion(1,i)*AACONV,
     >                  ion_vion(2,i)*AACONV,
     >                  ion_vion(3,i)*AACONV
         end do
         call Waterpsp_PrintXYZ(18)
         call util_flush(18)
      end if
  110 format(I5/)
  111 format(A2,5x, 6e12.4)

      return
      end



