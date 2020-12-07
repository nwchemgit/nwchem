
*     *****************************
*     *                           *
*     *     pspw_director_init    *
*     *                           *
*     *****************************
      subroutine pspw_director_init(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "util.fh"
#include "stdio.fh"

*     **** local variables ****
      logical   lprint
      integer   MASTER,taskid,l
      parameter (MASTER=0)

      character*255 filename,full_filename

      logical use_director
      common /pspw_director_common/ use_director

*     **** external functions ****
      logical  control_use_director,control_print
      external control_use_director,control_print


      use_director = control_use_director()

      if (use_director) then
         call Parallel_taskid(taskid)
         if (.not.btdb_cget(rtdb,'nwpw:director_filename',
     >                      1,filename))
     >      call util_file_prefix('director',filename)
            call util_file_name_noprefix(filename,.false.,.false.,
     >                                   full_filename)

          l = index(full_filename,' ') -1
          lprint= ((taskid.eq.MASTER).and.control_print(print_medium))
          if (lprint) write(luout,1210) full_filename(1:l)

      end if

 1210 FORMAT(/' pspw_director filename:',A)
      return
      end 

*     ***************************
*     *                         *
*     *     pspw_director       *
*     *                         *
*     ***************************
*
*  This routine initializes the director file, which
* is used to keep track of ion positions and velocities.

      subroutine pspw_director(rtdb)
      implicit none
      integer   rtdb

#include "btdb.fh"
#include "inp.fh"
#include "stdio.fh"

      integer   MASTER,taskid
      parameter (MASTER=0)

      logical found
      character*50 filename
      character*255 full_filename,line


      logical use_director
      common /pspw_director_common/ use_director

*     **** external functions ***
      real*8   lattice_omega
      integer  ion_nion,control_it_out
      external lattice_omega
      external ion_nion,control_it_out

      if (use_director) then

         call Parallel_taskid(taskid)

         if (.not.btdb_cget(rtdb,'nwpw:director_filename',1,filename))
     >       call util_file_prefix('director',filename)

         call util_file_name_noprefix(filename,.false.,.false.,
     >                                full_filename)

         call pspw_reset_cmd_director()

         if (taskid.eq.MASTER) then

            inquire(file=full_filename,exist=found)
            if (found) then
               open(unit=83,file=full_filename,
     >              form='formatted',status='old')
               do while (found)
                  read(83,'(A)',ERR=30,END=30) line
                  if (inp_strlen(line).gt.9) 
     >               call pspw_add_cmd_director(line)
               end do
   30          continue 
               close(83)

               call util_file_unlink(full_filename)

            end if
         end if

         call pspw_brdcst_cmd_director()
         call pspw_run_cmd_director()

      end if

      return
      end


*     *************************************
*     *                                   *
*     *     pspw_reset_cmd_director       *
*     *                                   *
*     *************************************
      subroutine pspw_reset_cmd_director()
      implicit none

*     **** pspw_director_lines common block ****
      character*255 lines(10)
      integer nlines
      common /pspw_director_lines/ lines,nlines

      nlines = 0
      return
      end 


*     *************************************
*     *                                   *
*     *      pspw_add_cmd_director        *
*     *                                   *
*     *************************************
      subroutine pspw_add_cmd_director(line)
      implicit none
      character*(*) line

*     **** pspw_director_lines common block ****
      character*255 lines(10)
      integer nlines
      common /pspw_director_lines/ lines,nlines

      nlines = nlines + 1
      lines(nlines) = line
      return
      end


*     *************************************
*     *                                   *
*     *      pspw_brdcst_cmd_director     *
*     *                                   *
*     *************************************
      subroutine pspw_brdcst_cmd_director()
      implicit none

*     **** local variables ****
      integer   MASTER,taskid
      parameter (MASTER=0)

      integer i

*     **** pspw_director_lines common block ****
      character*255 lines(10)
      integer nlines
      common /pspw_director_lines/ lines,nlines

      call Parallel_Brdcst_ivalue(MASTER,nlines)
      do i=1,nlines
         call pspw_director_brdcst_string(lines(i))
      end do
      return
      end




*     *************************************
*     *                                   *
*     *      pspw_run_cmd_director        *
*     *                                   *
*     *************************************
      subroutine pspw_run_cmd_director()
      implicit none

#include "bafdecls.fh"
#include "inp.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical found
      integer indx,indx2,ll1,ll2,i,ind,ii
      integer r0_ptr,r1_ptr,r2_ptr,v2_ptr,nion

      character*255 filename,cmd

*     **** pspw_director_lines common block ****
      character*255 lines(10)
      integer nlines
      common /pspw_director_lines/ lines,nlines

*     **** external functions ****
      integer  ion_nion,ion_rion_indx_ptr
      external ion_nion,ion_rion_indx_ptr
      real*8   ion_rion2,ion_rion,ion_vion
      external ion_rion2,ion_rion,ion_vion
      
      call Parallel_taskid(taskid)

      do i=1,nlines
         filename = ''
         cmd = ''

         filename = lines(i)
         indx     = index(filename," ")
         filename = filename(indx+1:)

         if (inp_strlen(filename).gt.3) then
            do while (filename(1:1).eq.' ') 
               filename = filename(2:)
            end do
      
            indx2 = index(filename," ")
            cmd      = filename(indx2+1:)
            filename = filename(:indx2-1)
            ll1 = inp_strlen(filename)
            ll2 = inp_strlen(cmd)

            if (inp_contains(.false.,"writestatus",lines(i),ind)) then

              if (taskid.eq.MASTER) then
                 write(*,*) "add_cmd_director: writestatus "
     >                      //filename(1:ll1)

                 open(unit=82,file=filename(1:ll1),
     >                action='write',position='append')
                 write(82,'(A)') cmd(1:ll2)
                 close(82)
              end if

            else if (inp_contains(.false.,"loadmovecs",
     >               lines(i),ind)) then
               if (taskid.eq.MASTER) then
                  write(*,*) "add_cmd_director: loadmovecs "
     >                       //filename(1:ll1)
               end if
               call psi_tmp_read_full_filename(filename)

            else if (inp_contains(.false.,"savemovecs",
     >               lines(i),ind)) then
               if (taskid.eq.MASTER) then
                  write(*,*) "add_cmd_director: savemovecs "
     >                       //filename(1:ll1)
               end if
               call psi_tmp_write_full_filename(filename)

            else if (inp_contains(.false.,"loadgeometry",
     >                            lines(i),ind)) then
               r0_ptr = ion_rion_indx_ptr(0)
               r1_ptr = ion_rion_indx_ptr(1)
               r2_ptr = ion_rion_indx_ptr(2)
               v2_ptr = ion_rion_indx_ptr(3)
               nion   = ion_nion()
               if (taskid.eq.MASTER) then
                  write(*,*) "add_cmd_director: loadgeometry "
     >                       //filename(1:ll1)
                  inquire(file=filename(1:ll1),exist=found)
                  if (found) then
                     open(unit=83,file=filename(1:ll1),
     >                    action='read',status='old')
                     do ii=1,nion
                        read(83,*) dbl_mb(r2_ptr+3*(ii-1)),
     >                             dbl_mb(r2_ptr+3*(ii-1)+1),
     >                             dbl_mb(r2_ptr+3*(ii-1)+2),
     >                             dbl_mb(r1_ptr+3*(ii-1)),
     >                             dbl_mb(r1_ptr+3*(ii-1)+1),
     >                             dbl_mb(r1_ptr+3*(ii-1)+2),
     >                             dbl_mb(r0_ptr+3*(ii-1)),
     >                             dbl_mb(r0_ptr+3*(ii-1)+1),
     >                             dbl_mb(r0_ptr+3*(ii-1)+2),
     >                             dbl_mb(v2_ptr+3*(ii-1)),
     >                             dbl_mb(v2_ptr+3*(ii-1)+1),
     >                             dbl_mb(v2_ptr+3*(ii-1)+2) 
                     end do
                     close(83)
                     call util_file_unlink(filename(1:ll1))
                  end if
               end if
               call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r2_ptr))
               call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r1_ptr))
               call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(r0_ptr))
               call Parallel_Brdcst_values(MASTER,3*nion,dbl_mb(v2_ptr))

               call Nose_zero_thermostats()
            

            else if (inp_contains(.false.,"savegeometry",
     >               lines(i),ind)) then
               r0_ptr = ion_rion_indx_ptr(0)
               r1_ptr = ion_rion_indx_ptr(1)
               r2_ptr = ion_rion_indx_ptr(2)
               v2_ptr = ion_rion_indx_ptr(3)
               nion   = ion_nion()
               if (taskid.eq.MASTER) then
                  write(*,*) "add_cmd_director: savegeometry "
     >                       //filename(1:ll1)
                  call util_file_unlink(filename(1:ll1))
                  open(unit=82,file=filename(1:ll1),
     >                 action='write',status='new')
                  do ii=1,nion
                     write(82,'(12E24.15)') dbl_mb(r2_ptr+3*(ii-1)),
     >                                     dbl_mb(r2_ptr+3*(ii-1)+1),
     >                                     dbl_mb(r2_ptr+3*(ii-1)+2),
     >                                     dbl_mb(r1_ptr+3*(ii-1)),
     >                                     dbl_mb(r1_ptr+3*(ii-1)+1),
     >                                     dbl_mb(r1_ptr+3*(ii-1)+2),
     >                                     dbl_mb(r0_ptr+3*(ii-1)),
     >                                     dbl_mb(r0_ptr+3*(ii-1)+1),
     >                                     dbl_mb(r0_ptr+3*(ii-1)+2),
     >                                     dbl_mb(v2_ptr+3*(ii-1)),
     >                                     dbl_mb(v2_ptr+3*(ii-1)+1),
     >                                     dbl_mb(v2_ptr+3*(ii-1)+2) 
                  end do
                  close(82)
               end if

            end if
         end if

      end do

      return
      end


*     *************************************
*     *                                   *
*     *    pspw_director_brdcst_string    *
*     *                                   *
*     *************************************
      subroutine pspw_director_brdcst_string(mystring)
      implicit none
      character*(*) mystring

#include "inp.fh"

      integer   MASTER,taskid
      parameter (MASTER=0)

      integer istring(255),ilen,i
      character*255 tmpstring

      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
         ilen = inp_strlen(mystring)
         do i=1,ilen
            istring(i) = ichar(mystring(i:i))
         end do
      end if
      call Parallel_Brdcst_ivalue(MASTER,ilen)
      call Parallel_Brdcst_ivalues(MASTER,ilen,istring)
      tmpstring = ''
      do i=1,ilen
         tmpstring(i:i) = char(istring(i))
      end do
      mystring = tmpstring

      return
      end

