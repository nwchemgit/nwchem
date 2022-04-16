*
* $Id$

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

#include "btdb.fh"


      integer   MASTER
      parameter (MASTER=0)

      logical value,found,found_bak
      integer taskid,l1,l2
      real*8 ch_tmp
      character*50 filename
      character*255 full_filename,full_bak


*     **** external functions ***
      real*8   lattice_omega
      integer  ion_nion,control_it_out,control_code
      external lattice_omega
      external ion_nion,control_it_out,control_code


      call Parallel_taskid(taskid)

      if ((control_code().eq.1).or.(control_code().eq.13)) then
         if (.not.btdb_cget(rtdb,'cpsd:ion_motion_filename',1,filename))
     >      call util_file_prefix('ion_motion',filename)
      else
         if (.not.btdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >      call util_file_prefix('ion_motion',filename)
      end if

      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                              full_filename)
      

*     **** produce MOTION FILE ****
      if (taskid.eq.MASTER) then

*        **** check for backup file ****
         call util_file_name_noprefix('MOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*)
            write(*,*) "MOTION99-bak exists:"
            l1=index(full_bak,' ')
            l2=index(full_filename,' ')
            write(*,*) "   Copying ",full_bak(1:l2),
     >                 " to ",full_filename(1:l2)
            write(*,*)
            call util_file_copy(full_bak,full_filename)
         end if


*        **** MOTION FILE already exists - parse to EOF ****
         inquire(file=full_filename,exist=found)
         if (found) then         

*          **** make a new backup file ***
           call util_file_copy(full_filename,full_bak)

           open(unit=19,file=full_filename,form='formatted',
     >          status='old')
           do while(.true.)
             read(19,*,ERR=30,END=30) ch_tmp
           end do
 30        continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1)||defined(GCC46)
           backspace 19
#endif

*        **** MOTION FILE does not exist ****
         else
            open(unit=19,file=full_filename,form='formatted',
     >           status='new')
         end if

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
      character*255 full_bak

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         close(unit=19)

*        **** remove backup file ***
         call util_file_name_noprefix('MOTION99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         call util_file_unlink(full_bak)
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

      integer   MASTER,taskid
      parameter (MASTER=0)

      integer i

*     **** external functions ***
      integer     ion_nion
      external    ion_nion
      real*8      ion_rion,ion_vion,lattice_omega,lattice_unita
      external    ion_rion,ion_vion,lattice_omega,lattice_unita
      character*2 ion_aname_nocarat
      external    ion_aname_nocarat
      character*4 ion_aname
      external    ion_aname

      call Parallel_taskid(taskid)
    
      if (taskid.eq.MASTER) then
         write(19,110) time,
     >                 ion_nion(),
     >                 lattice_omega(),
     >     lattice_unita(1,1),lattice_unita(2,1),lattice_unita(3,1),
     >     lattice_unita(1,2),lattice_unita(2,2),lattice_unita(3,2),
     >     lattice_unita(1,3),lattice_unita(2,3),lattice_unita(3,3)
         do i=1,ion_nion()
            write(19,111) i,ion_aname_nocarat(i),ion_aname(i),
     >                    ion_rion(1,i),
     >                    ion_rion(2,i),
     >                    ion_rion(3,i),
     >                    ion_vion(1,i),
     >                    ion_vion(2,i),
     >                    ion_vion(3,i)
         end do
         call util_flush(19)
      end if
  110 format(e14.6,I6,e14.6,9e14.6)
  111 format(i6,a3,a5,6e14.6)

      return
      end

