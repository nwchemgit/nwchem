*
* $Id: MOTION.F,v 1.6 2004-05-16 05:38:07 edo Exp $

*  Note - The format of the MOTION file was changed on 5-8-2002.
*
* Format of MOTION file:
*
*  time nion omega
*   x(1) y(1) z(1) vx(1) vy(1) vz(1)
*   x(2) y(2) z(2) vx(2) vy(2) vz(2)
*   x(3) y(3) z(3) vx(3) vy(3) vz(3)
*   ...
*   x(nion) y(nion) z(nion) vx(nion) vy(nion) vz(nion)
*  time2 nion omega
*   ...
*

*     ***************************
*     *				*
*     *	      MOTION_init	*
*     *				*
*     ***************************
*
*  This routine initializes the MOTION file, which
* is used to keep track of ion positions and velocities.
*

      subroutine MOTION_init(rtdb)
      implicit none 
      integer   rtdb

#include "rtdb.fh"


      integer   MASTER
      parameter (MASTER=0)

      logical value,found
      integer taskid
      real*8 ch_tmp
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

*        **** MOTION FILE already exists - parse to EOF ****
         if (found) then         
           open(unit=19,file=full_filename,form='formatted',
     >          status='old')
           do while(.true.)
             read(19,*,ERR=30,END=30) ch_tmp
           end do
 30        continue
#if defined(FUJITSU_SOLARIS) || defined(PSCALE) || defined(__crayx1)
           backspace 19
#endif

*        **** MOTION FILE does not exist ****
         else
            open(unit=19,file=full_filename,form='formatted',
     >           status='new')
         end if

c        ! taken out on 5-8-2002
c         open(unit=19,file=full_filename,form='formatted')
c         write(19,*) control_it_out(),(ion_nion()+pspw_qmmm_nion()),
c     >               lattice_omega()
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
*     *				*
*     *		MOTION_write    *
*     *				*
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
      real*8      ion_rion, ion_vion,lattice_omega
      real*8      pspw_qmmm_rion,pspw_qmmm_vion
      external    pspw_qmmm_found
      external    ion_nion,pspw_qmmm_nion
      external    ion_rion, ion_vion,lattice_omega
      external    pspw_qmmm_rion,pspw_qmmm_vion

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(19,110) time,
     >                 (ion_nion()+pspw_qmmm_nion()),
     >                 lattice_omega()
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
         call util_flush(19)
      end if
  110 format(e14.6,I6,e14.6)
  111 format(6e14.6)

      return
      end



