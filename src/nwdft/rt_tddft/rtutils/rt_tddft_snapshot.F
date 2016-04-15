      integer function rt_tddft_find_itreference (params)
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_find_itreference: "
      double precision, parameter :: thresh = 1d-9

C     == Variables ==
      integer it, itreference
      double precision tt
      logical found

      found = .false.
      itreference = -99

      do it = 1, params%nt
         tt = params%tmin + (it-1)*params%dt
         
         if ( abs(tt - params%viz_treference) .lt. thresh) then
            if (found) then
               call errquit (pname//"multiple values found",0,0)
            endif
            
            itreference = it
            found = .true.
         endif
      enddo

      if (.not. found) then
         call errquit (pname//"failed to find value",0,0)
      endif

      call ga_sync ()

      if (itreference < 1)
     $     call errquit (pname//"bad itreference",0,0)

      rt_tddft_find_itreference = itreference
      
      end function



      logical function rt_tddft_at_snap (params, tt)
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_at_snap: "
      double precision, parameter :: thresh = 1d-9

C     == Variables ==
      logical at_snap


C
C     if visualization is active, only make snapshots between the
C     specified range.
C
      if ((tt.le.params%viz_tend).and.(tt.ge.params%viz_tstart)) then
         at_snap = .true.
      else
         at_snap = .false.
      endif

      
C
C     also if this is the "reference" density matrix, make a snapshot
C
      if ( abs(tt - params%viz_treference) .lt. thresh) then
         at_snap = .true.
      endif

      at_snap = at_snap .and. params%viz_active
      
      rt_tddft_at_snap = at_snap

      end function



C====================================================================
C     
C     Generate a filename for density matrix dump (snapshot) with
C     with time step index.
C
      subroutine rt_tddft_snapshot_fname (params, it, tt, fout)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: it
      double precision, intent(in)  :: tt

      
C     == Outputs ==
      character(50), intent(out)    :: fout  !hardcoded str len; fix if changing fname format

      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_snapshot_fname: "

      if (params%nt .gt. 99999999)
     $     call errquit (pname//"nt too big; fix formatting", 0, 0)
      
C      write (fout, "(a,i0.10,a,e12.6,a)")
C     $        "ptot_ao_re.", it, ".", tt, "-au"

      write (fout, "(a,i0.10)")
     $        "ptot_ao_re.", it
      
      end subroutine


C====================================================================
      subroutine rt_tddft_snapshot_fname_cube (params, it, tt, fout)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: it
      double precision, intent(in)  :: tt

C     == Outputs ==
      character(50), intent(out)    :: fout  !note hardcoded from fname call

      
C     == Parameters ==
      character(*) ,parameter :: pname="rt_tddft_snapshot_fname_cube: "

      if (params%nt .gt. 99999999)
     $     call errquit (pname//"nt too big; fix formatting", 0, 0)
      
      if (params%viz_subgs) then
         write (fout, "(a,i0.10,a)")
     $        "density_subgs.", it, ".cube"
      else
         write (fout, "(a,i0.10,a)")
     $        "density.", it, ".cube"
      endif
      
      end subroutine



C====================================================================
C
C     Dumps the supplied real part of the dens mat in AO basis to file;
C     we later read in this snapshot to make density plots, etc.
C
      subroutine rt_tddft_snapshot_save (params, it, tt, g_densre_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      integer, intent(in)           :: g_densre_ao  !re part of *total* dens mat in AO basis

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_snapshot_save: "


C     == Variables ==
      character(255) fname   
      character(50) fname_tail  !note hardcoded str length; update if changing fname format
      double precision elapsed

      if (params%prof) call prof_start (elapsed)

      call rt_tddft_snapshot_fname (params, it, tt, fname_tail)
      call util_file_name (fname_tail, .false., .false., fname)

      if (.not. dmat_io_dump (g_densre_ao, trim(fname)))
     $     call errquit (pname//"failed to dump densao_re",0,0)

      if (params%prof) call prof_end (elapsed, "Saving snapshot")
      end subroutine


C====================================================================
C
C     Convert all stored density matrix snapshots (real part in AO
C     basis) to charge densities .cube files using dplot.  Note: dplot
C     will clean up after itself, which includes shutting down the
C     integrals, etc.  Call this at the very end.
C
      subroutine rt_tddft_snapshot_dplot (params, subgs)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      logical, intent(in)           :: subgs  !if true, subtract ground state density


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_snapshot_dplot: "


C     == Variables ==
      character(255) fname_dens, fname_dens_gs, fname_cube
      character(50) fname_tail  !note hardcoded str length; update if changing fname format
      character(50) fname_tail_cube
      integer it, itreference
      double precision tt
      integer isnap
      integer me
      integer num_basis,bas_handles(2)
      double precision elapsed
      character(64) title


C     == External ==
      logical, external :: dplot
C      logical, external :: rt_tddft_atpoint
      logical, external :: rt_tddft_at_snap
      integer, external :: rt_tddft_find_itreference


      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid ()

      if (me.eq.0) then
         call util_print_centered (luout,
     $        "Post-processing of density matrix snapshots", 20, .true.)

         write (luout, *) ""
         write (luout, *) ""
      endif


C     (determine index of density matrix we want to subtract)
      itreference = rt_tddft_find_itreference (params)
      
C
C     If desired, set up dplot to subtract ground state: rho = rho(t) -
C     rho(t=0):
C
      if (subgs) then
         
C         call rt_tddft_snapshot_fname(params,1, params%tmin, fname_tail)  !this is hardcoded for t=0

         call rt_tddft_snapshot_fname (params, itreference,
     $        params%viz_treference, fname_tail)
         
         call util_file_name(fname_tail, .false.,.false., fname_dens_gs)
         
         if (.not.rtdb_cput(params%rtdb,"dplot:File_Mat2",
     $        1,fname_dens_gs))
     $        call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)
         call rt_tddft_print_notice (
     $        "Subtracting reference density matrix")
      else
         call rt_tddft_print_notice (
     $        "Not subtracting reference density matrix")
      endif

c     int cleanup to avoid clashes in dplot
            call int_terminate

C
C     Loop over all times and convert each stored dens mat snapshot to
C     cube.
C
      isnap = 0
      it = 1
      do while (it .lt. params%nt)
         tt = params%tmin + (it-1)*params%dt

C         if (rt_tddft_atpoint(it, params%nt, params%nsnapshots)) then
         if (rt_tddft_at_snap (params, tt)) then

            isnap = isnap + 1

            if (me.eq.0) then
               write (luout, "(a,i0,a,f12.3,a)")
     $              "Postprocessing snapshot ",
     $              isnap, ", t = ", tt, " au ..."
            endif

C     (title of the slice is the time in atomic units)
            write (title, "(a,1es14.7,a,1es14.7,a)")
     $           "t = ", tt, " au = ", tt*au2fs, " fs"

            call rt_tddft_snapshot_fname (params, it, tt, fname_tail)
            call util_file_name (fname_tail,.false.,.false., fname_dens)

            call rt_tddft_snapshot_fname_cube (params, it, tt,
     $           fname_tail_cube)
            
C     changed to match new dplot api (no leading ./perm)
c$$$            call util_file_name (trim(fname_tail_cube),
c$$$     $        .false., .false., fname_cube)
            fname_cube = trim(fname_tail_cube)

C     (load dplot params in rtdb)
            if (.not.rtdb_cput(params%rtdb,'dplot:Title',1,title))
     &           call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)

            if (.not.rtdb_cput(params%rtdb,"dplot:File_Out" ,
     $           1,fname_cube))
     $           call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)
            
            if (.not.rtdb_cput(params%rtdb,"dplot:File_Mat1",
     $           1,fname_dens))
     $           call errquit(pname//"Write failed to rtdb",0,RTDB_ERR)

C     (call dplot to convert)
            if (.not. dplot (params%rtdb))
     $           call errquit (pname//"dplot call failed", 0, 0)
         endif

         it = it + 1
      enddo

c     int cleanup to avoid clashes in dplot
      bas_handles(1)=params%ao_bas_han
      num_basis=1
      if(params%iVcoul_opt.eq.1) then
         num_basis=2
         bas_handles(2)=params%cd_bas_han
      endif
      call int_init(params%rtdb,num_basis,bas_handles)

      if (params%prof) call prof_end (elapsed,
     $     "Converting snapshots to cube files")
      end subroutine
c $Id$
