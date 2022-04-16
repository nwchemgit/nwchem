c     $Id$
*     *********************************
*     *                               *
*     *          CIF_write            *
*     *                               *
*     *********************************
      subroutine CIF_write(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"

      integer taskid,MASTER
      parameter (MASTER=0)

      logical       shift_cell,found
      integer       flen,l
      character*50 filename
      character*255 full_filename
      character ch_tmp
      logical  control_print
      external control_print

      if (btdb_cget(rtdb,'nwpw:cif_filename',1,filename)) then
         flen  = index(filename,' ') - 1
         filename = filename(1:flen)//'.cif'
     
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)

         if (.not.btdb_get(rtdb,'nwpw:cif_shift_cell',
     >                     mt_log,1,shift_cell)) 
     >      shift_cell = .false.

         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
            l = index(full_filename,' ') -1
            if (control_print(print_medium)) then
              write(luout,*) 
              write(luout,*) 
     >        'Writing Crystallographic Information File:', 
     >                   full_filename(1:l)
              if (shift_cell) then
                write(luout,*) 
     >          '  - cell shifted (nwpw:cif_shift_cell .true.)'
              else
                write(luout,*) 
     >          '  - cell not shifted (nwpw:cif_shift_cell .false.)'
              end if
            endif

            inquire(file=full_filename,exist=found)

*           **** CIF FILE already exists - parse to EOF ****
            if (found) then
              open(unit=20,file=full_filename,form='formatted',
     >             status='old')
              do while(.true.)
                read(20,*,ERR=30,END=30) ch_tmp
              end do
 30           continue
#if defined(FUJITSU) || defined(PSCALE) || defined(__crayx1) || defined(GCC46)
              backspace 20
#endif

*           **** CIF FILE does not exist ****
            else
              open(unit=20,file=full_filename,form='formatted')
            end if

            call ion_Print_CIF(20,shift_cell)
            close(20)
         end if
      end if


      return
      end

