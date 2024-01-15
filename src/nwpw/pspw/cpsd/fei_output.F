*
* $Id$
*
      subroutine fei_init(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      logical       found,found_bak
      character*50  filename
      character*2 tag2
      character*4 tag4
      character*255 full_filename,full_bak
      integer nion,ii,iii,l1,l2
      real*8 E,x,y,z,fx,fy,fz
      real*8 a1x,a1y,a1z
      real*8 a2x,a2y,a2z
      real*8 a3x,a3y,a3z

*     **** external functions ****
      logical  control_Fei
      external control_Fei
      integer  control_code
      external control_code


      if (control_Fei()) then

      if ((control_code().eq.1).or.(control_code().eq.13)) then
           if (.not.btdb_cget(rtdb,'cpsd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)
        else
           if (.not.btdb_cget(rtdb,'cpmd:fei_filename',1,filename))
     >     call util_file_prefix('fei',filename)
        end if
        call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

        call util_file_name_noprefix(filename,.false.,
     >                               .false.,
     >                               full_filename)

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER)  then

*        **** check for backup file ****
         call util_file_name_noprefix('FEI99-bak',.false.,
     >                                .false.,
     >                                full_bak)
         inquire(file=full_bak,exist=found_bak)
         if (found_bak) then
            write(*,*) 
            write(*,*) "FEI99-bak exists:"
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

          open(unit=35,file=full_filename,form='formatted',status='old')
            do while (found)
               read(35,*,ERR=30,END=30) nion
               read(35,*,ERR=30,END=30) E
               read(35,*,ERR=30,END=30) a1x,a1y,a1z
               read(35,*,ERR=30,END=30) a2x,a2y,a2z
               read(35,*,ERR=30,END=30) a3x,a3y,a3z
               do ii=1,nion
                  read(35,*,ERR=30,END=30) iii,tag2,tag4,x,y,z,fx,fy,fz
               end do
            end do
  30        continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
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
      character*255 full_bak

*     **** external functions ****
      logical  control_Fei
      external control_Fei

      if (control_Fei()) then

        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
           close(unit=35)

*           **** remove backup file ***
            call util_file_name_noprefix('FEI99-bak',.false.,
     >                                .false.,
     >                                full_bak)
            call util_file_unlink(full_bak)
        end if

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
      character*2 ion_aname_nocarat
      character*4 ion_aname
      external control_Fei
      external ion_nion
      external lattice_unita,ion_rion
      external   ion_aname_nocarat
      external   ion_aname

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
          write(35,103) ii,ion_aname_nocarat(ii),ion_aname(ii),
     >                  ion_rion(1,ii),ion_rion(2,ii),ion_rion(3,ii),
     >                  ion_fion(1,ii),ion_fion(2,ii),ion_fion(3,ii)
        end do
        call util_flush(35)
        end if

  100   format(I8)
  101   format(F20.6)
  102   format(3F10.5/3F10.5/3F10.5)
  103   format(I5,A3,A5,6F10.5)

      end if

      return
      end
