*
* $Id: fei_output.F,v 1.4 2006-12-27 19:47:52 bylaska Exp $
*
      subroutine fei_init(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      logical       found
      character*30  filename
      character*255 full_filename
      integer nion,ii,iii
      real*8 E,x,y,z,fx,fy,fz
      real*8 a1x,a1y,a1z
      real*8 a2x,a2y,a2z
      real*8 a3x,a3y,a3z

*     **** external functions ****
      logical  control_Fei
      external control_Fei


      if (control_Fei()) then

        if (.not.rtdb_cget(rtdb,'cpmd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)

        call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

        call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER)  then
         inquire(file=full_filename,exist=found)
         if (found) then
          open(unit=35,file=full_filename,form='formatted',status='old')
            do while (found)
               read(35,*,ERR=30,END=30) nion
               read(35,*,ERR=30,END=30) E
               read(35,*,ERR=30,END=30) a1x,a1y,a1z
               read(35,*,ERR=30,END=30) a2x,a2y,a2z
               read(35,*,ERR=30,END=30) a3x,a3y,a3z
               do ii=1,nion
                  read(35,*,ERR=30,END=30) iii,x,y,z,fx,fy,fz
               end do
            end do
  30        continue
#if defined(FUJITSU_SOLARIS) || defined(PSCALE) || defined(__crayx1)
           backspace 35
#endif
         else
          open(unit=35,file=full_filename,form='formatted')
         end if
        end if

      end if
      return
      end

      subroutine fei_end()
      implicit none

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

*     **** external functions ****
      logical  control_Fei
      external control_Fei

      if (control_Fei()) then

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) close(unit=35)

      end if
      return
      end


      subroutine fei_output(E,ion_fion)
      implicit none
      real*8 E
      real*8  ion_fion(3,*)

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer ii,nion

*     **** external functions ****
      logical  control_Fei
      integer  ion_nion
      real*8   lattice_unita,ion_rion
      external control_Fei
      external ion_nion
      external lattice_unita,ion_rion

      if (control_Fei()) then

        call Parallel_taskid(taskid)

        if (taskid.eq.MASTER) then
        nion = ion_nion()
        write(35,100) nion
        write(35,101) E
        write(35,102) lattice_unita(1,1),
     >                lattice_unita(2,1),
     >                lattice_unita(3,1),
     >                lattice_unita(1,2),
     >                lattice_unita(2,2),
     >                lattice_unita(3,2),
     >                lattice_unita(1,3),
     >                lattice_unita(2,3),
     >                lattice_unita(3,3)
        do ii=1,nion
          write(35,103) ii,
     >                  ion_rion(1,ii),ion_rion(2,ii),ion_rion(3,ii),
     >                  ion_fion(1,ii),ion_fion(2,ii),ion_fion(3,ii)
        end do
        call util_flush(35)
        end if

  100   format(I8)
  101   format(F20.6)
  102   format(3F10.5/3F10.5/3F10.5)
  103   format(I5,6F10.5)

      end if

      return
      end
