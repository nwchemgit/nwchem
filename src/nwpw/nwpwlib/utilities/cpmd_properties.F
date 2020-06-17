c
c     $Id$
c

*     *********************************
*     *                               *
*     *          cpmd_properties      *
*     *                               *
*     *********************************
      subroutine cpmd_properties(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value,doit,oprint,mdxs,fei,dipole_motion
      character*80 filename
      character*255 full_filename,full_filename2

*     
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      if (.not.btdb_get(rtdb,'nwpw:cpmd_properties',mt_log,1,doit))
     >   doit = .false.

      if (oprint) then
         write(luout,*) 
         write(luout,*) "================================"
         write(luout,*) "==== CPMD Property Analysis ===="
         write(luout,*) "================================"
         if (doit) then
            write(luout,*) 
            write(luout,*) "CPMD property analysis is on. ",
     >         "To turn analysis off set cpmd_properties, i.e."
            write(luout,*) "nwpw"
            write(luout,*) "   cpmd_properties off"
            write(luout,*) "end"
            write(luout,*)
         else
            write(luout,*)
            write(luout,*) "CPMD property analysis is off. ",
     >         "To turn analysis on set cpmd_properties, i.e."
            write(luout,*) "nwpw"
            write(luout,*) "   cpmd_properties on"
            write(luout,*) "end"
            write(luout,*)
         end if
      end if

      if (doit) then

      if (.not.btdb_cget(rtdb,'cpmd:ion_motion_filename',1,filename))
     >    call util_file_prefix('ion_motion',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

      if (.not.btdb_cget(rtdb,'cpmd:unfolded_xyz_filename',1,filename))
     >    call util_file_prefix('unfolded.xyz',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename2)

      call pair_dist(full_filename)
      call pair_distab(full_filename)
      call velocity_autocorrelation(rtdb,full_filename)
      call velocity_spectrum(rtdb,full_filename)
      call generate_unfolded_xyz(full_filename,full_filename2)


     

*     **** perform shake force averaging if Fei set *****
      if (.not.btdb_get(rtdb,'nwpw:fei',mt_log,1,fei))
     >   fei = .false.

      if (fei) then
         call shake_analysis(rtdb)
         call FixIon_analysis(rtdb)
      end if

*     **** MD XS analysis ****
      if (.not.btdb_get(rtdb,'nwpw:md_xs',mt_log,1,mdxs)) mdxs=.false.
      if (mdxs) call md_xs_analysis(rtdb,full_filename)


*     **** perform dipole correlation and power spectrum if dipole_motion set *****
      if (.not.btdb_get(rtdb,'nwpw:dipole_motion',mt_log,
     >                  1,dipole_motion))
     >   dipole_motion = .false.

      if (dipole_motion) then
       if(.not.btdb_cget(rtdb,'nwpw:dipole_motion_filename',1,filename))
     >   call util_file_prefix('dipole_motion',filename)
       call util_file_name_noprefix(filename,.false.,
     >                              .false.,
     >                              full_filename)

         call dipole_motion_autocorrelation(rtdb,full_filename)
      end if


      end if

      return
      end

*     *********************************
*     *                               *
*     *      task_aimd_properties     *
*     *                               *
*     *********************************
*
*   This routine is a stand-alone task that calls cpmd_properties.
*
      logical function task_aimd_properties(rtdb)
      implicit none
      integer rtdb

#include "btdb.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,dummy,hasbak,hasxyz,mdxs,mdxsonly
      real*8 unita(3,3)
      character*80 filename,filenamexyz,filenamebak
      character*255 full_filename,full_filenamexyz

*     **** external functions ****
      logical  control_read
      external control_read
      integer  control_np_orbital
      external control_np_orbital

      value = .true.

      !call Parallel_Init()
      !value = value.and.control_read(999,rtdb)
      !call Parallel2d_Init(control_np_orbital())

*     **** save nwpw:cpmd_properties and change it .true. ****
      if (.not.btdb_get(rtdb,'nwpw:cpmd_properties',mt_log,1,dummy))
     >   dummy = .false.
      value = value.and.
     >        btdb_put(rtdb,'nwpw:cpmd_properties',mt_log,1,.true.)

*     **** if nwpw_aimd_properties_xyz exist, then                         ****
*     **** convert xyz filename to MOTION filename (prefix.ion_motion_tmp) ****
      hasxyz = .false.
      if (btdb_cget(rtdb,'nwpw:aimd_properties_xyz',1,filenamexyz)) then
         hasxyz = .true.
         if (.not.btdb_get(rtdb,'cell_default:unita',
     >                     mt_dbl,9,unita)) then
            unita(1,1) = 1.0d0 !*** unita=Identity will be interpreted by md_xs_analysis as no unit cell
            unita(2,1) = 0.0d0
            unita(3,1) = 0.0d0
            unita(1,2) = 0.0d0
            unita(2,2) = 1.0d0
            unita(3,2) = 0.0d0
            unita(1,3) = 0.0d0
            unita(2,3) = 0.0d0
            unita(3,3) = 1.0d0
         end if
         call util_file_name_noprefix(filenamexyz,.false.,
     >                                .false.,
     >                                full_filenamexyz)
         call util_file_prefix('ion_motion_tmp',filename)
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                                full_filename)
         call nwpw_xyz_to_motion(unita,
     >                           full_filenamexyz,
     >                           full_filename)
         
*        **** save nwpw:ion_motion_filename and change it to '.ion_motion_tmp' ****
         hasbak = .true.
         if (.not.btdb_cget(rtdb,'cpmd:ion_motion_filename',
     >                      1,filenamebak))
     >      hasbak = .false.
         value = value.and.
     >           btdb_cput(rtdb,'cpmd:ion_motion_filename',1,filename)

      end if

      call Parallel_Init()
      value = value.and.control_read(999,rtdb)
      call Parallel2d_Init(control_np_orbital())
  
*     **** MD XS analysis only ****
      if (.not.btdb_get(rtdb,'nwpw:cpmd_properties:md_xs_only',
     >                  mt_log,1,mdxsonly)) mdxsonly=.false.
      if (.not.btdb_get(rtdb,'nwpw:md_xs',mt_log,1,mdxs)) mdxs=.false.
      if (mdxsonly) then
         if (mdxs) call md_xs_analysis(rtdb,filename)

*     **** do all cpmd analysisl ****
      else
         call cpmd_properties(rtdb)
      end if

      call ga_sync()
      call Parallel2d_Finalize()
      call Parallel_Finalize()

*     **** reset nwpw:ion_motion_filename to its saved value ****
      if (hasxyz) then
         if (hasbak) then
            value = value.and.
     >              btdb_cput(rtdb,'cpmd:ion_motion_filename',
     >                         1,filenamebak)

         else
            value = value.and.
     >              rtdb_delete(rtdb,'cpmd:ion_motion_filename')
         end if
      end if

*     **** reset nwpw:cpmd_properties to its saved value ****
      value = value.and.
     >        btdb_put(rtdb,'nwpw:cpmd_properties',mt_log,1,dummy)

      task_aimd_properties = value
      return
      end 


*     *********************************
*     *                               *
*     *     nwpw_xyz_to_motion        *
*     *                               *
*     *********************************
      subroutine nwpw_xyz_to_motion(unita,xyzfilename,motionfilename)
      implicit none
      real*8        unita(3,3)
      character*(*) xyzfilename
      character*(*) motionfilename

#include "bafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"

c     **** parameters ****
      real*8    AACONV
      parameter (AACONV=0.529177d0)
      integer MASTER,taskid
      parameter (MASTER=0)

*     **** local variables ****
      logical found_motion
      integer i,natoms,nframes,t
      real*8  time1,volume
      real*8  x,y,z,v1,v2,v3
      character*2 tag2

*     *** synchronization ***
      call ga_sync()

      !call Parallel_taskid(taskid)
      taskid = ga_nodeid()

      if (taskid.eq.MASTER) then

      inquire(file=motionfilename,exist=found_motion)
      if (found_motion) then

         write(luout,*)
         write(luout,*) 
     >     '===== Generating MOTION file from xyz trajectory file ===='
         write(luout,*) 
         write(luout,*) "xyz filename         : ",xyzfilename
         write(luout,*) "MOTION filename      : ",motionfilename
         write(luout,*) "MOTION filename already exists. Not copying."
         write(luout,*) "Delete this file if you want to regenerate it."
         write(luout,*)

      else

      volume = unita(1,1)*(unita(2,2)*unita(3,3)-unita(3,2)*unita(2,3))
     >       + unita(2,1)*(unita(3,2)*unita(1,3)-unita(1,2)*unita(3,3))
     >       + unita(3,1)*(unita(1,2)*unita(2,3)-unita(2,2)*unita(1,3))
      volume = dabs(volume)


*     *** get the number of frames from MOTION ***
      call xyz_nframes(xyzfilename,nframes,natoms)

      write(luout,*)
      write(luout,*) 
     =  '===== Generating MOTION file from xyz trajectory file ===='
      write(luout,*) 
      write(luout,*) "nframes = ",nframes
      write(luout,*) "natoms  = ",natoms
      write(luout,*) 
      write(luout,*) "xyz filename         : ",xyzfilename
      write(luout,*) "MOTION filename      : ",motionfilename


      open(unit=18,file=xyzfilename,form='formatted')
      open(unit=19,file=motionfilename,form='formatted')

      time1 = 0.0d0
      do t=1,nframes
        
         read(18,*,ERR=30,END=30) natoms
         read(18,*,ERR=30,END=30) 
         write(19,110) time1,natoms,volume,unita

         do i=1,natoms
            read(18,*,ERR=30,END=30) tag2,x,y,z,v1,v2,v3
            write(19,111) i,tag2,tag2,
     >                    x/AACONV, y/AACONV, z/AACONV,
     >                    v1/AACONV,v2/AACONV,v3/AACONV
         end do
         time1 = time1 + 50.0d0
         call util_flush(19)

      end do
  30  continue
      close(unit=18)
      close(unit=19)

      end if
      end if

*     *** synchronization ***
      call ga_sync()

      return
  110 format(e14.6,I6,e14.6,9e14.6)
  111 format(i6,a3,a5,6e14.6)

      end

