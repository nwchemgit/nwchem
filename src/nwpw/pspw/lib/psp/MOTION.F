
*     ***************************
*     *				*
*     *		MOTION_init	*
*     *				*
*     ***************************
      subroutine MOTION_init(rtdb)
*
* $Id: MOTION.F,v 1.1 2002-01-14 21:06:28 bylaska Exp $
*
      implicit none 
      integer   rtdb

#include "rtdb.fh"


      integer   MASTER
      parameter (MASTER=0)

      logical value,found
      integer taskid
      character*30 filename
      character*255 full_filename


*     **** external functions ***
      real*8   lattice_omega
      integer  ion_nion,control_it_out,pspw_qmmm_nion
      external lattice_omega
      external ion_nion,control_it_out,pspw_qmmm_nion


      call Parallel_taskid(taskid)

      value = rtdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                    full_filename)

*     **** produce MOTION FILE ****
      if (taskid.eq.MASTER) then
         inquire(file=full_filename,exist=found)

         open(unit=19,file=full_filename,form='formatted')

         write(19,*) control_it_out(),(ion_nion()+pspw_qmmm_nion()),
     >               lattice_omega()
      end if

      return
      end


*     ***************************
*     *				*
*     *		MOTION_end 	*
*     *				*
*     ***************************
      subroutine MOTION_end()
      implicit none

      integer   MASTER
      parameter (MASTER=0)

      integer taskid

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         close(unit=19)
      end if

      return
      end

*     ***************************
*     *							*
*     *		MOTION_write		*
*     *							*
*     ***************************
      subroutine MOTION_write(time)
      implicit none
      real*8 time

      integer   MASTER
      parameter (MASTER=0)

      integer taskid,i

*     **** external functions ***
      logical     pspw_qmmm_found
      integer     ion_nion,pspw_qmmm_nion
      real*8      ion_rion, ion_vion
      real*8      pspw_qmmm_rion,pspw_qmmm_vion
      external    pspw_qmmm_found
      external    ion_nion,pspw_qmmm_nion
      external    ion_rion, ion_vion
      external    pspw_qmmm_rion,pspw_qmmm_vion

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(19,110) time
         do i=1,ion_nion()
            write(19,111) ion_rion(1,i),
     >                    ion_rion(2,i),
     >                    ion_rion(3,i),
     >                    ion_vion(1,i),
     >                    ion_vion(2,i),
     >                    ion_vion(3,i)
         end do
         if (pspw_qmmm_found()) then
           do i=1,pspw_qmmm_nion()
              write(19,111) pspw_qmmm_rion(1,i),
     >                      pspw_qmmm_rion(2,i),
     >                      pspw_qmmm_rion(3,i),
     >                      pspw_qmmm_vion(1,i),
     >                      pspw_qmmm_vion(2,i),
     >                      pspw_qmmm_vion(3,i)

           end do
         end if
      end if
  110 format(e14.6)
  111 format(6e14.6)

      return
      end



