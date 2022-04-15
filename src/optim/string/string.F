c $Id$

C> \ingroup task
C> @{
!
! Zero temperature string method (ZTS) to compute minimum energy paths (MEPs)
! The method works faster if the initial end points are both minima
! This code allows one to input a guess at the "middle" point also
! See:
! http://cims.nyu.edu/~eve2/string.htm (existed in 2011)
! W. E, W. Ren, and E. Vanden-Eijnden, Phys. Rev. B: 66, 052301 (2002).
! W. E, W. Ren, and E. Vanden-Eijnden, J. Chem. Phys.: 126, 164103 (2007).
! http://www.math.princeton.edu/string/index.html (existed in 2011)
!
!   It should be noted that the end points will go downhill because of 
!   this evolution in "time".  The non-end points also slip downhill, but 
!   the reparameterization forces them to still form a curve between
!   the two end points.  With that in mind, they should eventually
!   form a minimum energy path.
!

! The RTDB parameters that control the code are:
! string:tol   : The convergence criteria (based on RMS coordinate change)
! string:stepsize  : The "Time step" small=slow, big=unstable
! string:nbeads  : The number of points on the string more=>slow & stable
! string:maxit     : The number of time steps
! string:interpol  : Linear or spline or other.  Not very important.
! string:freeze1  : freeze endpoint 1
! string:freezen  : freeze endpoint N
! string:linopt   : Turn on/off linear optimization
! string:projection1  : Turn off/on projection of gradient
! string:pathguess  : which method for guessing path from 1 to TS to N


! This code requires at two or three geometry blocks to be defined
! "geometry"  Endpoint #1
! "endgeom"   Endpoint #2
! "midgeom"   Optional transition state
! At this point the code does not verify that they geometries are the same atoms
! in the same order, but they need to be for the code to work right.
! The main routine can handle having all geometries specified, but the 
! driver does support reading them in yet.

      logical function task_string(rtdb)
      implicit none

! External stuff
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"


! Input parameters that do not get changed
      integer rtdb   ! The runtime database

! Local variables for memory blocks
      integer l_dg, k_dg,l_energy, k_energy, l_dcoords, k_dcoords
      integer l_dVdx, k_dVdx,l_coords_old, k_coords_old
      integer l_wrk, k_wrk,l_coords, k_coords
      integer l_tmp_list, k_tmp_list
      integer l_mass,k_mass

! Local variables
      integer maxit,nbeads,natoms,middle,geom,geom2,size,i,interpol
      real*8 tol,stepsize,t,t1,t2,t3,rms1,rms2
      logical converged, freeze1, freezeN, newchain,hasmiddle,impose
      logical oprint,value
      character*255 geom_name
      integer geom_len,nbeads_exist
      integer lold_coords,kold_coords,j1,j2,shift,shift1,shift2,m
      integer string_algorithm
      integer ks,km,ke,ifit(2),wfit(2),nfit

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name


      oprint = (ga_nodeid().eq.0)

! Set convergence criteria
      if (.not. rtdb_get(rtdb,'string:tol',mt_dbl,1,tol))
     $      tol = 0.0001d0

! Set the time step size
      if (.not. rtdb_get(rtdb,'string:stepsize',mt_dbl,1,stepsize))
     $      stepsize = 0.1d0

! Number of points on the path
      if (.not. rtdb_get(rtdb,'string:nbeads',mt_int,1,nbeads))
     $      nbeads=10

! Number of steps
      if (.not. rtdb_get(rtdb,'string:maxit',mt_int,1,maxit))
     $      maxit = 25

! They type of interpolation
      if (.not. rtdb_get(rtdb,'string:interpol',mt_int,1,interpol))
     $      interpol = 1

! nhist=m
      if (.not. rtdb_get(rtdb,'string:m',mt_int,1,m)) m = 3

! algorithm
      if(.not.rtdb_get(rtdb,'string:algorithm',mt_int,1,
     >                 string_algorithm))
     >     string_algorithm = 3

! The general geometry structure
      if (.not. geom_create(geom, 'geometry'))
     &   call errquit('task_string: geom_create?', 911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &   call errquit('task_string: no geometry ', geom, RTDB_ERR)

! Get number of atoms
      if (.not. geom_ncent(geom,natoms))
     $     call errquit('task_string: natoms?',geom, GEOM_ERR)

! Get memory
      if(.not.ma_push_get(mt_dbl,nbeads,'dg',l_dg, k_dg))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbeads,'energy',l_energy, k_energy))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbeads,'dcoord',l_dcoords, k_dcoords))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      size = 3*natoms*nbeads
      if(.not.ma_push_get(mt_dbl,size,'dVdx',l_dVdx, k_dVdx))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'coords_old',l_coords_old,
     $     k_coords_old))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'wrk',l_wrk, k_wrk))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'coords',l_coords, k_coords))
     $     call errquit('task_string: memory', size, MA_ERR)

      size = 3*natoms*nbeads*(2*m+1)
      if(.not.ma_push_get(mt_dbl,size,'tmp_list',
     $                    l_tmp_list, k_tmp_list))
     $     call errquit('task_string: memory', size, MA_ERR)

      size = 3*natoms
      if(.not.ma_push_get(mt_dbl,size,'mass',l_mass, k_mass))
     $     call errquit('task_string: memory', size, MA_ERR)
      if (.not.geom_masses_get(geom,natoms,dbl_mb(k_mass)))
     $     call errquit('task_string: masses?',geom,GEOM_ERR)
      call dscal(natoms,1822.89d0,dbl_mb(k_mass),1)


! Warning: No way to read in all nbead geomtries

      newchain = .false.
      if (.not.rtdb_get(rtdb,'bead_list:new',mt_log,1,newchain))
     >  newchain = .true.

*     **** initialize bead list ****
      if (newchain) then

!        *** Get the coordinates of geometry 1 ***
         if (.not. geom_cart_coords_get(geom, dbl_mb(k_coords)))
     $        call errquit('task_string: geom?',0, GEOM_ERR)

!        *** Get the coordinates of last geometry ***
         size = k_coords + (nbeads-1)*3*natoms
         ke = size
         ks = k_coords
         km = -1
         if (.not. geom_create(geom2, 'geom2'))
     &      call errquit('task_string: geom_create?', 921, GEOM_ERR)
         if (.not. geom_rtdb_load(rtdb, geom2, 'endgeom'))
     &      call errquit('task_string: no endgeom ', geom2, RTDB_ERR)
         if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     $           call errquit('task_string: geom?',0, GEOM_ERR)

!        *** Middle geometry -- if it exists ***
         if (.not.rtdb_get(rtdb,'string:hasmiddle',mt_log,1,hasmiddle))
     >      hasmiddle = .false.
         middle = 0
         if (hasmiddle) then
            if (.not. geom_rtdb_load(rtdb, geom2, 'midgeom')) then
               middle = 0
               hasmiddle = .false.
            else
              middle = (nbeads/2)+1
              size = k_coords + (middle-1)*3*natoms
              km = size
              if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     $              call errquit('task_string: geom?',0, GEOM_ERR)
            end if
         end if

         if (.not. geom_destroy(geom2))
     &      call errquit('task_string: geom_destroy?', geom2, GEOM_ERR)

c        **** impose ****
         if (.not.rtdb_get(rtdb,'string:impose',mt_log,1,impose))
     >      impose = .false.

         if (impose) then
            value= MA_push_get(mt_int,(2*natoms),'ifit',ifit(2),ifit(1))
            value = value.and.
     >               MA_push_get(mt_dbl,(natoms),'wfit',wfit(2),wfit(1))
            if(.not.value) 
     >        call errquit('string_initial_path failed',3,MA_ERR)

            if (hasmiddle) then
               call neb_impose(natoms,dbl_mb(ks),dbl_mb(km),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)

               if (oprint) then
                  write(luout,*)
     >            " - Imposing midgeom geometry onto",
     >            " starting (geometry) geometry"
                  write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
                  write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
               end if

               call neb_impose(natoms,dbl_mb(km),dbl_mb(ke),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)
               if (oprint) then
                  write(luout,*)
     >            " - Imposing endgeom geometry onto midgeom geometry"
                  write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
                  write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
               end if
            else
               call neb_impose(natoms,dbl_mb(ks),dbl_mb(ke),
     >                   nfit,int_mb(ifit(1)),dbl_mb(wfit(1)),rms1,rms2)
               if (oprint) then
                  write(luout,*)
     >            " - Imposing endgeom geometry",
     >            " onto starting (geometry) geometry"
                  write(luout,'(A,F10.6)') "    + initial rmsq = ",rms1
                  write(luout,'(A,F10.6)') "    + imposed rmsq = ",rms2
               end if
            end if

            value =           MA_pop_stack(wfit(2))
            value = value.and.MA_pop_stack(ifit(2))
            if(.not.value) 
     >         call errquit('string_initial_path failed',4,MA_ERR)
         end if

         


*     **** restart - read in old bead list ****
      else
         middle = -1

         if (.not.rtdb_get(rtdb,'bead:size',
     >                     mt_int,1,nbeads_exist))
     >     call errquit('getting string:nbeads_exist failed',4,RTDB_ERR)

         if (oprint) 
     >   write(luout,*) "nbeads,nbeads_exist=",nbeads,nbeads_exist

         if (.not. geom_create(geom2, 'geom2'))
     &      call errquit('task_string: geom_create?', 921, GEOM_ERR)

*        **** resize the bead list ****
         if (nbeads_exist.ne.nbeads) then

            !size = 3*natoms*nbeads*(2*m+1)
            size = 3*natoms*nbeads
            if(.not.ma_push_get(mt_dbl,size,'old_coords',
     >          lold_coords, kold_coords))
     >          call errquit('task_string: memory', size, MA_ERR)

            do i=1,nbeads_exist
               size = kold_coords + (i-1)*3*natoms
               geom_name   = 'bead'//bead_index_name(i)//':geom'
               geom_len = inp_strlen(geom_name)
               if (.not.geom_rtdb_load(rtdb,geom2,
     >             geom_name(1:geom_len)))
     >            call errquit('task_string: geom2?',0, GEOM_ERR)
               if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     >            call errquit('task_string: geom2?',0, GEOM_ERR)
            end do
            do i=1,nbeads
              t = (i-1)/dble(nbeads-1)

              j1 = t*(nbeads_exist-1) + 1
              j2 = j1+1
              t1 = (j1-1)/dble(nbeads_exist-1)
              t2 = (j2-1)/dble(nbeads_exist-1)
              t3 = (t-t1)/(t2-t1)

              if (j2.gt.nbeads_exist) then
                 t3 = 0.0d0
                 j2=nbeads_exist
              end if

              shift  = (i-1)*3*natoms
              shift1 = (j1-1)*3*natoms
              shift2 = (j2-1)*3*natoms

              call dcopy(3*natoms,dbl_mb(kold_coords+shift1),1,
     >                   dbl_mb(k_coords+shift),1)
              call dscal(3*natoms,(1.0d0-t3),dbl_mb(k_coords+shift),1)
              call daxpy(3*natoms,t3,dbl_mb(kold_coords+shift2),1,
     >                   dbl_mb(k_coords+shift),1)
           end do
            if (.not.MA_pop_stack(lold_coords))
     >          call errquit('task_string: memory', size, MA_ERR)

*        **** read in the bead list ****
         else
            do i=1,nbeads_exist
               size = k_coords + (i-1)*3*natoms
               geom_name   = 'bead'//bead_index_name(i)//':geom'
               geom_len = inp_strlen(geom_name)
               if (oprint) 
     >            write(luout,*) "loading geom2=",geom_name(1:geom_len)
               if (.not.geom_rtdb_load(rtdb,geom2,
     >             geom_name(1:geom_len)))
     >            call errquit('task_string: geom2?',0, GEOM_ERR)
               if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     >            call errquit('task_string: geom2?',0, GEOM_ERR)
            end do
         end if

         if (.not. geom_destroy(geom2))
     &      call errquit('task_string: geom2_destroy?',geom2,GEOM_ERR)

      end if

! Do we freeze bead 1 or N?
      if (.not. rtdb_get(rtdb, 'string:freeze1', mt_log,1,freeze1))
     &        freeze1 = .false.
      if (.not. rtdb_get(rtdb, 'string:freezen', mt_log,1,freezeN))
     &        freezeN = .false.


! Call the code
      if (.not. rtdb_put(rtdb,'opt:string',mt_log,1,.true.))
     $     call errquit('task_string: rtdb_put failed',0,
     &       RTDB_ERR)
      call zts_meps(maxit,nbeads,tol,stepsize,string_algorithm,
     &    natoms,middle,
     &    dbl_mb(k_coords),geom,rtdb,dbl_mb(k_dg),dbl_mb(k_energy),
     &    dbl_mb(k_dVdx),dbl_mb(k_dcoords),dbl_mb(k_coords_old),
     &    dbl_mb(k_wrk),
     &    m, dbl_mb(k_tmp_list),
     &       dbl_mb(k_tmp_list+3*natoms*nbeads*m),
     &       dbl_mb(k_tmp_list+3*natoms*nbeads*2*m),
     &       dbl_mb(k_mass),
     &    interpol,freeze1,freezeN,converged)
      if (.not. rtdb_put(rtdb,'opt:string',mt_log,1,.false.))
     $     call errquit('task_string: rtdb_put failed',0,
     &       RTDB_ERR)


!     *** Set the coordinates of beads to new geometry objects ***
      if (.not.rtdb_put(rtdb,'bead:size',mt_int,1,nbeads))
     >   call errquit('setting bead:size failed',4,RTDB_ERR)

      do i=1,nbeads
         size = k_coords + (i-1)*3*natoms
         geom_name   = 'bead'//bead_index_name(i)//':geom'
         geom_len = inp_strlen(geom_name)
         if (.not. geom_cart_coords_set(geom, dbl_mb(size)))
     $      call errquit('task_string: geom?',0, GEOM_ERR)
         if (.not.geom_rtdb_store(rtdb,geom,geom_name(1:geom_len)))
     $      call errquit('task_string: geom?',0, GEOM_ERR)
      end do

      if (.not.rtdb_put(rtdb,'bead_list:new',mt_log,1,.false.))
     > call errquit('setting bead_list:new failed',4,RTDB_ERR)


! Return memory - Being lazy, just pop is all in one fell swoop
      if(.not. ma_chop_stack(l_dg))
     &   call errquit('task_string: ma?',99,MA_ERR)

      if (.not. geom_destroy(geom))
     &   call errquit('task_string: geom_destroy?', geom, GEOM_ERR)

! Done
      !task_string = converged
      task_string = .true.
      return
      end
C> @}

*     ****************************************************
*     *                                                  *
*     *                 zts_meps                         *
*     *                                                  *
*     ****************************************************
      subroutine zts_meps(maxit,nbeads,tol,stepsize,string_algorithm,
     &    natoms,middle,
     &    coords,geom,rtdb,dg,energy,dVdx,dcoords,coords_old,wrk,
     &    m,coords_list,dVdx_list,HdVdx,masses,
     &    interpol,freeze1,freezeN,converged)
      implicit none
      
! External stuff
      real*8 zts_distance
      external zts_distance
      integer ga_nodeid
      external ga_nodeid
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"

! Input parameters that do not get changed
      integer maxit  ! Maximum number of iterations (200?)
      integer nbeads  ! Number of images or beads along the path (3 or more)
      real*8 tol ! Tolerance for convergence  (1.0d-7)
      real*8 stepsize ! The "time" step size (0.0001)
      integer string_algorithm
      integer natoms ! The number of atoms
      integer middle ! The "middle" coordinates bead number - 0 if none, -1 all
                     !      of the coordinates are filled
      integer rtdb   ! The runtime database
      integer interpol ! The style of interpolation
      logical freeze1, freezeN ! Freeze ends of string

! Input/output
      real*8  coords(3,natoms,nbeads)
      integer geom   ! A geometry record, should match molecule of interest
      logical converged  ! Did we converge
      
! Variables used locally, but not passed in or out
      real*8 dg(nbeads) ! 0,.5,1 for three, 0,1/3,2/3,1 for four
      real*8 energy(nbeads)
      real*8 dVdx(3,natoms,nbeads)
      real*8 dcoords(nbeads)
      real*8 coords_old(3,natoms,nbeads)
      real*8 wrk(3,natoms,nbeads)

      integer m
      real*8 coords_list(*)
      real*8 dVdx_list(*)
      real*8 HdVdx(3,natoms,nbeads)
      real*8 masses(*)

! Local variables
      integer i, j, k, nstep,print_shift,print_count,ii,ng,shift,itm
      real*8 tmp, tmp2,t
      logical badgeom,oprint,stalled
      real*8 g1(3),g2(3),g3(3)
      integer iistart, iiend
      real*8 hess, e0, e1, a2, dsgrad,sum,sum2,e3
      logical linopt
      logical projection1,found,status,finishedstep
      integer pathguess,manyjj,goodjj,jj,algorithm
      real*8  sum0,sum0_old,alpha,stepsize0,stepsize1,stepsizem
      real*8  gmax,grms,xmax,xrms
      real*8  gmax0,grms0,xmax0,xrms0

      character     ch_tmp
      character*255 filename,full_filename,perm_name,movecs_name
      character*255 filename2,full_filename2
      integer movecslen,permlen

      character*7 bead_index_name
      external    bead_index_name

      double precision ddot
      external ddot

! Setup the problem
! "Secret" options with defaults that generally do not need changed
      if (.not. rtdb_get(rtdb, 'string:linopt', mt_log,1,linopt))
     &        linopt = .false.
      if (.not.rtdb_get(rtdb,'string:projection1',mt_log,1,projection1))
     &        projection1 = .false.
      if (.not. rtdb_get(rtdb,'string:pathguess',mt_int,1,pathguess))
     $      pathguess = 2

      if (.not.rtdb_get(rtdb,'string:print_shift',mt_int,1,print_shift))
     >   print_shift = 0

      ng = 3*natoms*nbeads
      oprint = (ga_nodeid() .eq. 0)

      do i = 1, nbeads
        dg(i) = dble(i-1)/dble(nbeads-1)
      enddo
      
      if (middle .eq. 0) then  ! Flexible linear interpolation to get TS
        print_count = 0
        call zts_guessall(natoms,nbeads,coords,geom)
      else if (middle .gt. 1 .and. middle .lt. nbeads) then ! TS guess
        print_count = 0
        call zts_guess(natoms,nbeads,middle,coords,pathguess)
      else if (middle .eq. -1) then  ! We know exactly what we are doing
        ! Do nothing.  All coordinates were setup on entry.
         if (.not.rtdb_get(rtdb,'string:print_count',
     >                     mt_int,1,print_count))
     >      print_count = 0
      else  ! This makes no sense
        call errquit('zts_meps: fatal error', 0, 0)
      endif

*     **** set .string_epath filename ****s
      call util_file_prefix('string_epath',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

! Now we do the work
      if (ga_nodeid() .eq. 0) then
        write(luout,'(a)')'@zts'
        write(luout,'(a)')'@zts String method.'
        write(luout,'(a,f9.5)')'@zts Temperature          = ', 0.0d0
        write(luout,'(a,f9.5)')'@zts Covergence Tolerance = ', tol
        write(luout,'(a,f9.5)')'@zts Step Size            = ', stepsize
        write(luout,'(a,i9)')  '@zts Maximum Time Steps   = ', maxit
        write(luout,'(a,i9)')  '@zts Number of replicas   = ', nbeads
        write(luout,'(a,i9)')  '@zts Number of histories  = ', m
        if (string_algorithm.eq.0) then
           write(luout,'(a)') '@zts algorithm            = 0'//
     >                        ' (QN Fixed Point)'
        else if (string_algorithm.eq.1) then
           write(luout,'(a)') '@zts algorithm            = 1'//
     >                        ' (Damped Verlet)'
        else 
           write(luout,'(a)') '@zts algorithm            = 3'//
     >                        ' (QN Fixed Point - Damped Verlet)'
        end if
        write(luout,'(a,i9)')  '@zts String Interpolator  = ', interpol
        if (freeze1) then
          write(luout,'(a)')   '@zts First Replica        = frozen'
        else
          write(luout,'(a)')   '@zts First Replica        = moves'
        endif
        if (freezeN) then
          write(luout,'(a)')   '@zts Last Replica         = frozen'
        else
          write(luout,'(a)')   '@zts Last Replica         = moves'
        endif
        write(luout,'(a)')'@zts'
        write(luout,'(a)')
     &    '@zts  Step       xrms       xmax        E start       '
     &     //'E middle          E end          E max      E average'
     &     //' Walltime'
      endif
      if (maxit .lt. 1) then ! Just printing out initial geometries
         do i = 1, nbeads
           energy(i) = 0.0d0
         enddo
      endif

!     *** Look for atoms that are too close***
      do i = 1, nbeads
         call zts_Robinson(natoms,geom,coords(1,1,i),badgeom)
      end do
!     ***  Project out net motion ***
      do i = 2, nbeads
        call zts_min_motion(natoms,coords(1,1,i),coords(1,1,i-1),geom)
      end do

      !*** Calculate the energy and forces for each bead ***
      call zts_runall(rtdb,geom,natoms,nbeads,0,
     >                      freeze1,freezeN,projection1,
     >                      coords,energy,dVdx)
      sum0 = ddot(ng,dVdx,1,dVdx,1)
      sum0_old = sum0
      manyjj = 0
      goodjj = 0
      stepsize0 = 1.0d0
      stepsize1 = stepsize
      stepsizem = stepsize/16.0d0
      
      algorithm = string_algorithm
      gmax = 0.0d0
      grms = 0.0d0
      if (string_algorithm.ge.2) then
         gmax = 0.0d0
         do i=1,nbeads
            do j=1,natoms
               do k=1,3
                  if (dabs(dVdx(k,j,i)).gt.gmax) gmax=dabs(dVdx(k,j,i))
               end do
            end do
         end do
         if (gmax.gt.0.5d0) then
            algorithm = 1
            stepsize = stepsize1
         else
            algorithm  = 0
            stepsize = stepsize0
         end if
      end if

      call dcopy(ng,coords,1,coords_old,1)
      itm = 0
      do nstep = 1, maxit
       if (oprint) write(luout,*) "string: iteration #",nstep

       !**** lmbfgs update ****
       if (algorithm.eq.0) then
         if (oprint) write(*,*) "string: Fixed Point step"
         itm = itm+1

         call dcopy(3*natoms*nbeads,coords,1,coords_old,1)

         if (itm.le.m) then
             shift = (itm-1)*ng+1
             call dcopy(ng,coords_old,1,coords_list(shift),1)
             call dcopy(ng,dVdx,1,dVdx_list(shift),1)
             call neb_lmbfgs(ng,itm,
     >                       coords_list,
     >                       dVdx_list,
     >                       HdVdx)
          else
             do ii=1,m-1
                shift = (ii-1)*ng+1
                call dcopy(ng,coords_list(shift+ng),1,
     >                        coords_list(shift),   1)
                call dcopy(ng,dVdx_list(shift+ng),1,
     >                        dVdx_list(shift),   1)
              end do
             shift = (m-1)*ng+1
             call dcopy(ng,coords_old,1,coords_list(shift),1)
             call dcopy(ng,dVdx,1,dVdx_list(shift),1)
             call neb_lmbfgs(ng,m,
     >                       coords_list,
     >                       dVdx_list,
     >                       HdVdx)
          end if
          sum  = ddot(ng,dVdx,1,HdVdx,1) /dble(nbeads)
          sum2 = ddot(ng,HdVdx,1,HdVdx,1)/dble(nbeads)
          if (oprint) write(luout,*) "string: <g|s>=",sum
          if (oprint) write(luout,*) "string: <s|s>=",sum2
          if((sum.le.0.0d0).or.(sum .gt.1.0d0).or.(sum2.gt.10.0d0)) then
             call dcopy(ng,dVdx,1,HdVdx,1)
             itm = 0
             if (oprint) write(luout,*) "string: s=g"
          end if

          finishedstep = .false.
          alpha = stepsize
          jj = 0
          do while ((.not.finishedstep).and.(jj.lt.2))
             jj = jj + 1

             sum  = ddot(ng,dVdx,1,HdVdx,1) 
             if ((sum.le.0).or.(sum0.gt.sum0_old)) then
                call dcopy(ng,dVdx,1,HdVdx,1)
                !itm = 0
                if (oprint) write(luout,*) "string: s=g"
             end if
             call string_move(alpha,interpol,geom,natoms,nbeads,
     >                      dg,wrk,dcoords,coords_old,HdVdx,coords)

             call zts_runall(rtdb,geom,natoms,nbeads,nstep,
     >                      freeze1,freezeN,projection1,
     >                      coords,energy,dVdx)
             sum0 = ddot(ng,dVdx,1,dVdx,1)

             finishedstep = (sum0.le.sum0_old).or.
     >                      (stepsize.le.stepsizem)
             if (.not.finishedstep) then
                alpha = 0.5d0*alpha
                manyjj = manyjj + 1
                goodjj = 0
                if (manyjj.gt.2) then
                   stepsize = 0.5d0*stepsize
                   if (stepsize.lt.stepsizem) stepsize = stepsizem
                   manyjj = 0
                   if (oprint) 
     >             write(luout,*) "string: reducing stepsize=",stepsize
                end if
             end if
             if ((sum0.le.sum0_old).and.(jj.eq.1)) then
                goodjj = goodjj + 1
                if (goodjj.gt.10) then
                   stepsize = 2.0d0*stepsize
                   if (stepsize.gt.stepsize0) stepsize = stepsize0
                   goodjj = 0
                end if
             end if
             if (oprint)
     >          write(luout,*) "string: sum0,sum0_old=",
     >                     sum0,sum0_old,jj,finishedstep,alpha,itm

          end do
          sum0_old = sum0

       !**** damped verlet update ****
        else
         if (oprint) write(*,*) "string: Damped Verlet step, stepsize=",
     >                          stepsize
          call string_verlet_move(stepsize,interpol,geom,natoms,nbeads,
     >                    dg,wrk,masses,dcoords,coords_old,dVdx,coords)

          call zts_runall(rtdb,geom,natoms,nbeads,nstep,
     >                    freeze1,freezeN,projection1,
     >                    coords,energy,dVdx)

        end if


         gmax0 = gmax
         grms0 = grms
         xrms0 = xrms
         xmax0 = xmax

         !**** Check convergence ****
         gmax = 0.0d0
         grms = 0.0d0
         do i=1,nbeads
         do j=1,natoms
         do k=1,3
            grms = grms + dVdx(k,j,i)**2
            if (dabs(dVdx(k,j,i)).gt.gmax) gmax = dabs(dVdx(k,j,i))
         end do
         end do
         end do
         grms = dsqrt(grms/dble(ng))
         xrms = 0.0d0
         xmax = 0.0d0
         do i = 1, nbeads
           tmp = zts_distance(natoms,coords(1,1,i),coords_old(1,1,i))
           xrms = xrms + tmp
           if (tmp.gt.xmax) xmax = tmp
         enddo
         xrms = xrms/dble(nbeads*natoms)

         if (oprint) write(luout,*) "string: gmax,grms,xrms,xmax=",
     >                              gmax,grms,xrms,xmax

         stalled = ((gmax.ge.gmax0).or.(grms.ge.grms0)).and.(nstep.gt.1)

         !*** neb algorithm 3 - switch between fixed point and damped Verlet ***
         if (string_algorithm.ge.3) then
      
            !*** switch to fixed point if Xmax less than 1.0 ***
            if ((gmax.lt.0.5d0).and.
     >          (gmax.lt.gmax0).and.
     >          (grms.lt.grms0).and.
     >          (algorithm.ne.0)    ) then
               algorithm = 0
               itm = 0
               if (stepsize.gt.stepsize0) stepsize = stepsize0
               if (oprint) write(luout,*) "string: switching to fixed ",
     >                                    " point, stepsize=",stepsize
            end if

            !*** if stalled switch to damped Verlet ***
            if ((stalled).and.(algorithm.ne.1)) then
               algorithm = 1
               call dcopy(ng,coords,1,coords_old,1)
               if (stepsize.gt.stepsize1) stepsize = stepsize1
               if (oprint) write(luout,*) "string: switching to damped",
     >                                    " Verlet, stepsize=",stepsize
            end if

         end if

         if (ga_nodeid().eq.0) then
           tmp2 = energy(1)
           do i = 2, nbeads
             if (energy(i).gt.tmp2) tmp2=energy(i)
           enddo
           e3 = 0.0d0
           do i=1,nbeads
              e3 = e3 + energy(i)
           end do
           e3 = e3/dble(nbeads)
           write(luout,
     &    '(a5,i5,1x,f10.6,1x,f10.6,1x,f14.7,1x,f14.7,1x,f14.7,1x,f14.7,
     &       1x,F14.7,f9.1)'
     &            )'@zts ',nstep,xrms,xmax,energy(1),
     &                            energy(1+nbeads/2),
     &                            energy(nbeads),tmp2,e3,util_wallsec()


*           ***************************************
*           ***** print out .string_epath file ****
*           ***************************************
            inquire(file=full_filename,exist=found)
*           **** FILE already exists - parse to EOF ****
            if (found) then
              open(unit=19,file=full_filename,form='formatted',
     >             status='old')
              do while(.true.)
                read(19,*,ERR=30,END=30) ch_tmp
              end do
 30           continue
#if defined(FUJITSU) || defined(PSCALE) || defined(SOLARIS) || defined(__crayx1) || defined(GCC46)
              backspace 19
#endif
            write(19,*) " "
*           **** .string_epath FILE does not exist ****
            else
              open(unit=19,file=full_filename,form='formatted')
            end if
            write(19,*)
     > "#-------------------------------------------------------"
            write(19,*)  "# String Path iteration  = ",nstep
            write(19,*)  "# algorithm              = ",string_algorithm
            write(19,*)  "# String Interpolator    = ",interpol
            write(19,*)  "# nbeads                 = ",nbeads
            write(19,*)  "# step Size              = ",stepsize
            write(19,*)  "# convergence            = ",xrms,xmax
            write(19,*)  "# convergence tol (Xrms,Xmax)=",tol,tol
            write(19,*)
     > "#-------------------------------------------------------"
            write(luout,*) "string: Path Energy #",nstep
            do i=1,nbeads
               t = (i-1)/dble(nbeads-1)
               write(19,*) t,energy(i)
               write(luout,*) "string: ",i,energy(i)
            end do
            close(19)

         endif
         print_count = print_count + 1
         if (print_shift.gt.0) then
            if (mod(print_count,print_shift).eq.0) then
               call util_file_prefix(
     >            'stringpath'//bead_index_name(print_count)//'.xyz',
     >            filename2)
               call util_file_name_noprefix(filename2,.false.,
     >                             .false.,
     >                             full_filename2)

               if (ga_nodeid() .eq. 0) then
                   open(unit=23,file=full_filename2,form='formatted')
                   do i=1,nbeads
                      if (.not.geom_cart_coords_set(geom,coords(1,1,i)))
     >                 call errquit('zts_meps:set geometry',i,GEOM_ERR)
                      if (.not.geom_print_xyz(geom,23))
     >                call errquit('zts_meps:print geometry',i,GEOM_ERR)
                   end do
                   close(23)
               end if
            end if
         end if
         if (.not.rtdb_put(rtdb,'string:print_count',
     >                     mt_int,1,print_count))
     >      call errquit('setting print:print_count failed',4,RTDB_ERR)

         if ((xrms.lt.tol).and.(xmax.lt.tol)) goto 747
      end do


 747  continue
      if (nstep .ge. maxit) then
         if (oprint) 
     >   write(luout,'(a)') 
     >     '@zts The string calculation failed to converge'
         converged = .false.
      else
         if (oprint)
     &   write(luout,'(a)') '@zts The string calculation converged'
         converged = .true.
      endif


*     **** write out xyz trajectory ****
      call util_file_prefix('string_final.xyz',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
      if (ga_nodeid() .eq. 0) then
      open(unit=19,file=full_filename,form='formatted')

       write(luout,*)
       do i = 1, nbeads
         write(luout,'(a16,i5,a21,f20.12)') '@zts Bead number = ', i ,
     &       ' Potential Energy = ', energy(i)
         if ( .not. geom_cart_coords_set( geom, coords(1,1,i) ) )
     &        call errquit('zts_meps: set geometry',i,GEOM_ERR)
         if ( .not. geom_print_xyz(geom,6))
     &        call errquit('zts_meps: print geometry',i,GEOM_ERR)
         if ( .not. geom_print_xyz(geom,19))
     &        call errquit('zts_meps: print geometry',i,GEOM_ERR)
       enddo
       close(19)
       end if


*     **** write out final path energies ****
      call util_file_prefix('string_final_epath',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
      if (ga_nodeid() .eq. 0) then
         open(unit=19,file=full_filename,form='formatted')
         write(19,*)
     > "#-------------------------------------------------------"
         write(19,*)  "# String Path iteration  = ",nstep
         write(19,*)  "# algorithm              = ",string_algorithm
         write(19,*)  "# String Interpolator    = ",interpol
         write(19,*)  "# nbeads                 = ",nbeads
         write(19,*)  "# step Size              = ",stepsize
         write(19,*)  "# convergence            = ",tmp
         write(19,*)  "# convergence tol (Xrms) = ",tol
         write(19,*)
     > "#-------------------------------------------------------"
         do i=1,nbeads
            t = (i-1)/dble(nbeads-1)
            write(19,*) t,energy(i)
         end do
         close(19)
      endif

      if (ga_nodeid() .eq. 0) write(luout,*)
      
      end
      
      subroutine zts_linear_fit(nvar,nbeads,Y,Xin,wrk)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses linear interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! Work storage
      real*8 wrk(nvar,nbeads)
      ! local variables
      integer ivar, ibead, point
      real*8 X, XinRight1, XinLeft1, frac

      point = 1
      do ibead = 2, nbeads-1
        ! Find where this point is in the data
        X = dble(ibead-1)/dble(nbeads-1)  ! runs from zero to one
 11     if (Xin(point+1) .le. X) then
          point = point + 1
          goto 11
        endif
        if (point .eq. nbeads) point = nbeads-1 ! paranoid
        ! Solve linear interpolation function for data
        XinLeft1 = Xin(point)
        XinRight1 = Xin(point+1)
        frac = (X-XinLeft1)/(XinRight1-XinLeft1)
        do ivar = 1, nvar
          wrk(ivar,ibead) = (1.0d0-frac)*Y(ivar,point)+
     &                      frac*Y(ivar,point+1)
        enddo 
      enddo
      do ibead = 2, nbeads-1 ! end points do not change
        do ivar = 1, nvar
          Y(ivar,ibead) = wrk(ivar,ibead)
        enddo
      enddo
      
      return
      end
      
      subroutine zts_spline_fit(nvar,nbeads,Y,Xin,dg)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses spline interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      real*8 dg(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! local variables
      integer ivar, ibead
      real*8 X
      
      integer MAX_N
      PARAMETER (MAX_N=750)
      real*8 Y_tmp(MAX_N),fit(MAX_N)

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)

      do ivar = 1, nvar  ! loop over all the coordinates
        do ibead = 1, nbeads
          Y_tmp(ibead) = Y(ivar,ibead)
        enddo
        call zts_spline_Setup(Xin,Y_tmp,nbeads,fit) ! (x,y,n,fit)

        do ibead = 2, nbeads-1
          X = dg(ibead)
          call zts_spline(Xin,Y_tmp,nbeads,fit,X,Y(ivar,ibead)) ! (x,y,n,fit,xin,yout)
        enddo
      enddo
      return
      end
      
      subroutine zts_spline_Setup(x,y,n,fit)
      ! Input
      integer n
      real*8 x(n),y(n)
      ! Output
      real*8 fit(n)
      ! Local
      integer MAX_N
      PARAMETER (MAX_N=750)
      integer i
      real*8 tmp1,tmp2,u(MAX_N)
! This routine solves the tri-diagonal problem of a natural spline
      fit(1)=0.0d0
      u(1)=0.0d0
      
      if (n .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)

      do i=2,n-1
        tmp1=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        tmp2=tmp1*fit(i-1)+2
        fit(i)=(tmp1-1.0d0)/tmp2
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-tmp1*u(i-1))/tmp2
      enddo

      fit(n)=0.0d0
      do i=n-1,1,-1
        fit(i)=fit(i)*fit(i+1)+u(i)
      enddo
      return
      end

      subroutine zts_spline(x,y,n,fit,xin,yout)
      integer n
      real*8 xin,yout,x(n),fit(n),y(n)
      integer try,high,low
      real*8 a,b,h
      low=1
      high=n
      do while (high-low .gt. 1)
        try=(high+low)/2
        if(x(try).gt.xin)then
          high=try
        else
          low=try
        endif
      enddo
      h=x(high)-x(low)
      if (h.eq.0.0d0) then
        yout=(y(low)+y(high))/2.0d0 ! That is not good
      else
        a=(x(high)-xin)/h
        b=(xin-x(low))/h
        yout=a*y(low)+b*y(high)+((a**3-a)*fit(low)+(b**3-b)*fit(high))*
     &           (h**2)/6.0d0
      endif
      return
      end

*     ****************************************
*     *                                      *
*     *           zts_min_motion             *
*     *                                      *
*     ****************************************
      subroutine zts_min_motion(natoms,new ,old, geom)
      implicit none
#include "errquit.fh"
      integer natoms, geom
      real*8 new(3,natoms), old(3,natoms)

      real*8 improve
      integer isystype
      logical geom_systype_get
      external geom_systype_get

      if (.not. geom_systype_get(geom, isystype))
     $     call errquit('zts_min_motion: isystype?',0, GEOM_ERR)
      
      improve = 1.0d10
      do while (improve .gt. 0.001d0)
         improve = 0.0d0
         call zts_xyz(1,natoms,new,old,improve)
         call zts_xyz(2,natoms,new,old,improve)
         call zts_xyz(3,natoms,new,old,improve)
         if (isystype .eq. 0) then
           call zts_xyz(4,natoms,new,old,improve)
           call zts_xyz(5,natoms,new,old,improve)
           call zts_xyz(6,natoms,new,old,improve)
         else if (isystype .eq. 1) then
           call zts_xyz(1,natoms,new,old,improve) ! Does not change X
         endif
      enddo

      return
      end

      subroutine zts_xyz(which,natoms,new,old,improve)
      implicit none
      real*8 zts_distance
      external zts_distance
      integer which
      integer natoms
      real*8 new(3,natoms), old(3,natoms)
      real*8 improve
      
      real*8 tmp(3,natoms)
      real*8  left,  right,  try, middle
      real*8 xleft, xright, xtry
      real*8 initial, final

      initial = zts_distance(natoms,old,new)
      xleft   = -0.1d0
      xright  =  0.1d0
      middle  = initial

 99   continue  ! Walk unto we bracket the minimum
      call zts_copy_coords(natoms,new,tmp,which,xleft)
      left = zts_distance(natoms,old,tmp)
      if (left .lt. middle) then
         middle = left
         call zts_copy_coords(natoms,tmp,new,0,0.0d0)
         goto 99
      endif
      call zts_copy_coords(natoms,new,tmp,which,xright)
      right = zts_distance(natoms,old,tmp)
      if (right .lt. middle) then
         middle = right
         call zts_copy_coords(natoms,tmp,new,0,0.0d0)
         goto 99
      endif

! Now we have a middle that is better than left and right - go hunting
! xmiddle is always zero, xleft is negative, xright is positive
! We don't try to be smart with this linear search, since
! it is fast and non-linear

      do while ( (xright - xleft) .gt. 0.001d0)
        xtry = xright*0.5d0
        call zts_copy_coords(natoms,new,tmp,which,xtry)
        try = zts_distance(natoms,old,tmp)
        if (try .ge. middle) then
          right  = try
          xright = xtry
        else
          xleft  = -xtry
          left   = middle
          xright = xright-xtry
          right  = right
          middle = try
          call zts_copy_coords(natoms,tmp,new,0,0.0d0)
        endif
        xtry = xleft*0.5d0
        call zts_copy_coords(natoms,new,tmp,which,xtry)
        try = zts_distance(natoms,old,tmp)
        if (try .ge. middle) then
          left  = try
          xleft = xtry
        else
          xright  = -xtry
          right   = middle
          xleft   = xleft-xtry
          left    = left
          middle  = try
          call zts_copy_coords(natoms,tmp,new,0,0.0d0)
        endif
      enddo

      final = zts_distance(natoms,old,new)
      improve = improve + (initial - final)
      return
      end

! Make this a function so that we can easily mass weight it
! or otherwise play with it
      real*8 function zts_distance(n,c1,c2)
      implicit none
#include "errquit.fh"
      integer n
      real*8 c1(3,n),c2(3,n)
      integer j
      real*8 tmp

! Warning:  Does not deal with periodic wrapping yet

      tmp = 0.0d0
      do j = 1, n
           tmp = tmp + sqrt(
     &          (c1(1,j) - c2(1,j))**2 +
     &          (c1(2,j) - c2(2,j))**2 +
     &          (c1(3,j) - c2(3,j))**2)
      enddo
      zts_distance = tmp
      return
      end
      
! By abstracting this out, we can treat angles and distances the same above

      subroutine zts_copy_coords(n,in,out,which,shift)
      implicit none
      integer n
      real*8 in(3,n),out(3,n)
      integer which
      real*8 shift
      
      integer i,j,k
      real*8 tmp1, tmp2
      
      
      do i = 1, n
          out(1,i) = in(1,i)
          out(2,i) = in(2,i)
          out(3,i) = in(3,i)
      enddo
      if (which .eq. 0 .or. shift .eq. 0.0d0) then
!       Nothing to do
      else if (which .ge. 1 .and. which .le. 3) then
        do i = 1, n
          out(which,i) = out(which,i) + shift
        enddo
      else if (which .ge. 4 .and. which .le. 6) then ! Angles
        j = which - 3
        k = which - 2
        if (k .eq. 4) k = 1
        do i = 1, n
           tmp1 = out(j,i)*cos(shift) - out(k,i)*sin(shift)
           tmp2 = out(j,i)*sin(shift) + out(k,i)*cos(shift)
           out(j,i) = tmp1
           out(k,i) = tmp2
        enddo
      else
        call errquit('zts_copy_coords: bad which',0,0)
      endif
      return
      end

! If atoms are too close, then fix it
! Subroutine is named after the phrase "Danger, Will Robinson"
! Warning, no support for periodic images yet
       subroutine zts_Robinson(natoms,geom,coords,stillbad)
       implicit none

! Global includes
#include "nwc_const.fh"
       logical geom_cart_get_charges
       external geom_cart_get_charges
! Input only
       integer natoms
       integer geom
! Output only
       logical stillbad
! Input/Output
       real*8 coords(3,natoms)
! Local variables
       integer max_bad
       parameter (max_bad = 250)
       real*8 bad(4,max_bad)
       real*8 typei, typej
       integer is_h
       real*8 min_h_h, min_x_h, min_x_x, min_dist
       integer trys
       integer num_bad
       real*8 siner, coser, cosed, sined 
       real*8 r
       integer i,j,try,ij
       real*8 charges(nw_max_atom)

! Warning: These assume Angstom units and non-fractional
! Otherwise they are way too hard a constraint
       min_h_h = 0.7d0
       min_x_h = 0.9d0
       min_x_x = 1.0d0

       trys = 10

       if (.not. geom_cart_get_charges(geom,natoms,charges))
     1        call errquit('zts_Robinson:charges',i, 0)

       do try = 1, trys
        stillbad=.false.
       
        num_bad=0

!       Find the bad ones

        do i=1,natoms
         typei = charges(i)
         do j=1,i-1
           typej = charges(j)
           r = sqrt((coords(1,j) - coords(1,i))**2
     1          +   (coords(2,j) - coords(2,i))**2
     2          +   (coords(3,j) - coords(3,i))**2)
           is_h = 0
           if (typei .le. 2) is_h = is_h + 1
           if (typej .le. 2) is_h = is_h + 1

           if (((r .lt. min_h_h) .and. (is_h .eq. 2)) .or.
     1         ((r .lt. min_x_h) .and. (is_h .eq. 1)) .or.
     2         ((r .lt. min_x_x) .and. (is_h .eq. 0))) then
               num_bad=num_bad+1
               if(num_bad.gt.max_bad) then
                  call errquit('zts_Robinson:too many bad',num_bad, 0)
               endif
               bad(1,num_bad)=i
               bad(2,num_bad)=j
               bad(3,num_bad)=typei
               bad(4,num_bad)=typej
           endif
         enddo
        enddo

!       Fix the bad ones

        do ij=1,num_bad
          stillbad = .true.
          typei = bad(3,ij)
          typej = bad(4,ij)
          if (typej .lt. typei) then  ! Move lighter atom
            i = nint(bad(2,ij))
            j = nint(bad(1,ij))
          else
            i = nint(bad(1,ij))
            j = nint(bad(2,ij))
          endif
          is_h = 0
          if (typei .le. 2) is_h = is_h + 1
          if (typej .le. 2) is_h = is_h + 1

          if (is_h .eq. 2) then
            min_dist = min_h_h
          else if (is_h .eq. 1) then
            min_dist = min_x_h
          else
            min_dist = min_x_x
          endif
          if      (try .le. 2) then
             min_dist = min_dist / 3.0d0
          else if (try. le. 4) then
             min_dist = min_dist / 2.0d0
          else if (try .le. 5) then
             min_dist = min_dist / 1.5d0
          endif
          ! What we have now
          r = sqrt((coords(1,i)-coords(1,j))**2 + 
     1             (coords(2,i)-coords(2,j))**2 +
     2             (coords(3,i)-coords(3,j))**2)
          if (r .lt. 0.0001d0) then
           ! different amounts in each direction to make sure symmetry is broken
           if (coords(3,i) .gt. coords(3,j)) then
             coords(3,i) = coords(3,i) + 0.001d0
           else
             coords(3,i) = coords(3,i) - 0.001d0
           endif
           if (coords(2,i) .gt. coords(2,j)) then
             coords(2,i) = coords(2,i) + 0.0001d0
           else
             coords(2,i) = coords(2,i) - 0.0001d0
           endif
           if (coords(1,i) .gt. coords(1,j)) then
             coords(1,i) = coords(1,i) + 0.00001d0
           else
             coords(1,i) = coords(1,i) - 0.00001d0
           endif
          else if (r .lt. min_dist) then
           siner = (coords(3,i)-coords(3,j))/r
           coser = sqrt(1.0d0-siner*siner)
           if(coser .lt. 1d-8) then ! Force 0*0/0 to be zero
             cosed = 0.0d0
             sined = 0.0d0
           else
             cosed = (coords(1,i)-coords(1,j))/(r*coser)
             sined = (coords(2,i)-coords(2,j))/(r*coser)
           endif
           ! What we get instead
           coords(3,i) = min_dist*siner + coords(3,j)
           coords(2,i) = (min_dist * coser) * sined + coords(2,j)
           coords(1,i) = (min_dist * coser) * cosed + coords(1,j)
          endif
        enddo
        if (.not. stillbad) goto 395
       enddo ! "try" loop

 395   continue

       return
       end


      subroutine zts_akima_fit(nvar,nbeads,Y,Xin,dg)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses Akima spline interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      real*8 dg(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! local variables
      integer ivar, ibead
      real*8 X
      
      integer MAX_N
      PARAMETER (MAX_N=750)
      real*8 Y_tmp(MAX_N)

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)
! WARNING: Would be faster, if we did a "setup" like natural spline code
      do ivar = 1, nvar  ! loop over all the coordinates
        do ibead = 1, nbeads
          Y_tmp(ibead) = Y(ivar,ibead)
        enddo
        do ibead = 2, nbeads-1
          X = dg(ibead)
          call zts_akima_spline(nbeads,Xin,Y_tmp,X,Y(ivar,ibead))
        enddo
      enddo
      return
      end
      

!
! Akima cubic spline interpolation - This style avoids oscillations
! See: Akima, Hiroshi
! A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures
! J. ACM 17(4)
! October 1970
! 589-602
! doi = 10.1145/321607.321609
!

      subroutine zts_akima_spline(ndim,x,y,xin,yout)
      implicit none
! Input
      integer ndim
      real*8 x(ndim)
      real*8 y(ndim)
      real*8 xin
! Output
      real*8 yout
! Local
      integer low, high, try, i
      real*8 dydx(5),ddydx(4),t(2)
      real*8 tmp, h, dx

! First we find where we are
      low=1
      high=ndim
      do while (high-low .gt. 1)
        try=(high+low)/2
        if(x(try).gt.xin)then
          high=try
        else
          low=try
        endif
      enddo
      try = high ! Value to the right, not left

      if(try.eq.2) then
          do i=3,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(2) = 2.0d0*dydx(3) - dydx(4)
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
      else if(try.eq.3) then
          do i=2,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
      else if(try.eq.ndim) then
          do i=1,3
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(4) = 2.0d0*dydx(3) - dydx(2)
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
      else if(try.eq.ndim-1) then
          do i=1,4
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
      else
          do i=1,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
      end if
      do i=1,4
          ddydx(i) = abs(dydx(i+1) - dydx(i))
      end do
      do i=1,2
          tmp = ddydx(i+2) + ddydx(i)
          if(tmp .eq. 0) then
            t(i) = 0.0d0
          else
            t(i) = (ddydx(i+2)*dydx(i+1)+ddydx(i)*dydx(i+2))/tmp
          end if
      end do
      h = x(try)-x(try-1)
      dx = xin - x(try-1)
      yout = y(try-1)
     1      + (dx)   *(t(1))
     2      + (dx**2)*(3.0d0*dydx(3)-2.0d0*t(1)-t(2))/h
     3      + (dx**3)*(t(1)+t(2)-2.0d0*dydx(3))      /(h**2)
    
      return
      end



      subroutine zts_guess(natoms,nbeads,middle,coords,guess)
      implicit none
! Input
      integer natoms, nbeads, middle,guess
! Input/Output
      real*8 coords(3,natoms,nbeads)
! Local variables
      integer i,j,k
      real*8 dgg
      real*8 a,b,c
      real*8 pi
      parameter (pi=3.141592653589793238462643383279d0)

!     Simple linear interpolation

      if (guess .eq. 1) then
        do i = 2, middle-1
         dgg = dble(i-1)/dble(middle-1)
         do j = 1, natoms
           do k = 1, 3
            coords(k,j,i) = coords(k,j,1)*(1-dgg) +
     &                      coords(k,j,middle)*dgg
           enddo
         enddo
        enddo
        do i = middle+1, nbeads-1
         dgg = dble(i-middle)/dble(nbeads-middle)
         do j = 1, natoms
           do k = 1, 3
            coords(k,j,i) = coords(k,j,middle)*(1-dgg) +
     &                      coords(k,j,nbeads)*dgg
           enddo
         enddo
        enddo
      else if (guess .eq. 2) then ! Quadratic in x,y,z
        do j = 1, natoms
          do k = 1, 3
            dgg = dble(middle-1)/dble(nbeads-1)
            a=coords(k,j,1)
            b=((coords(k,j,middle)-coords(k,j,1))-
     &         (coords(k,j,nbeads)-coords(k,j,1))*dgg**2) /
     &        (dgg-dgg**2) 
            c=(coords(k,j,nbeads)-coords(k,j,1)) - b
            do i = 2, nbeads-1
              dgg = dble(i-1)/dble(nbeads-1)
              coords(k,j,i) = a+b*dgg+c*dgg**2
            enddo
          enddo
        enddo
      else if (guess .eq. 3) then ! Quadratic both parts with dx=0 at ends
        ! this just turns out to be linear with more points need ends
        do i = 2, middle-1
         dgg = dble(i-1)/dble(middle-1)
         do j = 1, natoms
           do k = 1, 3
            a = coords(k,j,1)
            b = coords(k,j,middle) - coords(k,j,1)
            coords(k,j,i) = a + b*(3*dgg**2 + 2*dgg**3)
           enddo
         enddo
        enddo
        do i = middle+1, nbeads-1
         dgg = dble(i-middle)/dble(nbeads-middle)
         do j = 1, natoms
           do k = 1, 3
            a = coords(k,j,middle)
            b = coords(k,j,nbeads) - coords(k,j,middle)
            coords(k,j,i) = a + b*(3*dgg**2 + 2*dgg**3)
           enddo
         enddo
        enddo
      else if (guess .eq. 4) then ! f=a+b*x+c*sin(pi*x)
        do j = 1, natoms
          do k = 1, 3
            dgg = dble(middle-1)/dble(nbeads-1)
            a=coords(k,j,1)
            b=coords(k,j,nbeads) - coords(k,j,1)
            c=((coords(k,j,middle)-coords(k,j,1)) - b*dgg)/sin(pi*dgg)
            do i = 2, nbeads-1
              dgg = dble(i-1)/dble(nbeads-1)
              coords(k,j,i) = a+b*dgg+c*sin(pi*dgg)
            enddo
          enddo
        enddo
      else
        call errquit('zts_guess: Bad Guess type', guess, 0)
      endif

      return
      end

      subroutine zts_guessall(natoms,nbeads,coords,geom)
      implicit none
!     Input
      integer natoms, nbeads, geom
!     Input/Output
      real*8 coords(3,natoms,nbeads)
!     Local variabls
      integer i, ix, iix, k, j
      integer MAX_N
      PARAMETER (MAX_N=750)
      logical done(MAX_N)
      logical badgeom
      real*8 dgg

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)
      done(1) = .true.
      do i = 2, nbeads - 1
        done(i) = .false.
      enddo
      done(nbeads) = .true.
 625  continue ! Loop start
        ix = -1
        iix= -1
        do i = 1, nbeads ! Find un-done region
          if (ix .eq. -1 .and. .not. done(i)) then
            ix = i-1
          endif
          if (ix .ne. -1 .and. done(i)) then
            iix = i
            goto 626
          endif
        enddo
        goto 627 ! Nothing found, we are done
 626    continue
        i = (iix+ix)/2
        dgg = dble(i-ix)/dble(iix-ix)
        do j = 1, natoms
          do k = 1, 3
            coords(k,j,i) = coords(k,j,ix )*(1-dgg) +
     &                      coords(k,j,iix)*dgg
          enddo
        enddo
        call zts_Robinson(natoms,geom,coords(1,1,i),badgeom)
        done(i) = .true.
        goto 625
      ! End of loop
 627  continue

        return
        end

      subroutine zts_guessall_periodic(natoms,nbeads,coords,geom,amat)
      implicit none
!     Input
      integer natoms, nbeads, geom
!     Input/Output
      real*8 coords(3,natoms,nbeads)
      real*8 amat(3,3)
!     Local variabls
      integer i1,i2,i3
      integer i, ix, iix, k, j
      integer MAX_N
      PARAMETER (MAX_N=750)
      logical done(MAX_N)
      logical badgeom
      real*8 dgg
      real*8 dx,dy,dz,bmat(3,3),vol
      real*8 f1,f2,f3


      !*** generate bmat ***
      bmat(1,1) = amat(2,2)*amat(3,3) - amat(3,2)*amat(2,3)
      bmat(2,1) = amat(3,2)*amat(1,3) - amat(1,2)*amat(3,3)
      bmat(3,1) = amat(1,2)*amat(2,3) - amat(2,2)*amat(1,3)

      bmat(1,2) = amat(2,3)*amat(3,1) - amat(3,3)*amat(2,1)
      bmat(2,2) = amat(3,3)*amat(1,1) - amat(1,3)*amat(3,1)
      bmat(3,2) = amat(1,3)*amat(2,1) - amat(2,3)*amat(2,1)

      bmat(1,3) = amat(2,1)*amat(3,2) - amat(3,1)*amat(2,2) 
      bmat(2,3) = amat(3,1)*amat(1,2) - amat(1,1)*amat(3,2) 
      bmat(3,3) = amat(1,1)*amat(2,2) - amat(2,1)*amat(1,2) 
      vol=amat(1,1)*bmat(1,1)+amat(2,1)*bmat(2,1)+amat(3,1)*bmat(3,1)
      call dscal(9,1.0d0/vol,bmat,1)


      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)
      done(1) = .true.
      do i = 2, nbeads - 1
        done(i) = .false.
      enddo
      done(nbeads) = .true.
 625  continue ! Loop start
        ix = -1
        iix= -1
        do i = 1, nbeads ! Find un-done region
          if (ix .eq. -1 .and. .not. done(i)) then
            ix = i-1
          endif
          if (ix .ne. -1 .and. done(i)) then
            iix = i
            goto 626
          endif
        enddo
        goto 627 ! Nothing found, we are done
 626    continue
        i = (iix+ix)/2
        dgg = dble(i-ix)/dble(iix-ix)
       
        do j = 1, natoms
          dx = coords(1,j,iix) - coords(1,j,ix )
          dy = coords(2,j,iix) - coords(2,j,ix )
          dz = coords(3,j,iix) - coords(3,j,ix )
          f1 = dx*bmat(1,1) + dy*bmat(2,1) + dz*bmat(3,1)
          f2 = dx*bmat(1,2) + dy*bmat(2,2) + dz*bmat(3,2)
          f3 = dx*bmat(1,3) + dy*bmat(2,3) + dz*bmat(3,3)
          f1 = f1 - DNINT(f1)
          f2 = f2 - DNINT(f2)
          f3 = f3 - DNINT(f3)
          dx = coords(1,j,ix ) 
     >       + dgg*(amat(1,1)*f1+amat(1,2)*f2+amat(1,3)*f3)
          dy = coords(2,j,ix ) 
     >       + dgg*(amat(2,1)*f1+amat(2,2)*f2+amat(2,3)*f3)
          dz = coords(3,j,ix ) 
     >       + dgg*(amat(3,1)*f1+amat(3,2)*f2+amat(3,3)*f3)

          coords(1,j,i) = dx
          coords(2,j,i) = dy
          coords(3,j,i) = dz

        enddo
        !call zts_Robinson(natoms,geom,coords(1,1,i),badgeom)
        done(i) = .true.
        goto 625
      ! End of loop
 627  continue

        return
        end







*     ****************************************
*     *                                      *
*     *             zts_runall               *
*     *                                      *
*     ****************************************
      subroutine zts_runall(rtdb,geom,natoms,nbeads,nstep,
     >                      freeze1,freezeN,projection1,
     >                      coords,energy,dVdx)
      implicit none
      integer rtdb,geom
      integer natoms,nbeads,nstep
      logical freeze1,freezeN,projection1,tmpVar
      real*8 coords(3,natoms,nbeads)
      real*8 energy(nbeads)
      real*8 dVdx(3,natoms,nbeads)

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "util.fh"
#include "inp.fh"

*     **** local variables ****
      logical skip1,skip2,status,oprint,ignore
      integer i,j,k
      integer permlen,movecslen
      character*255 movecs_name,perm_name
      real*8 g1(3),g2(3),g3(3),tmp

*     **** external functions ****
      logical  task_gradient
      external task_gradient
      character*7 bead_index_name
      external    bead_index_name

      oprint = (ga_nodeid().eq.0)

      !*** Calculate the energy and forces for each bead ***
      call dcopy(3*natoms*nbeads,0.0d0,0,dVdx,1)
      do i=1,nbeads

         if ( .not. geom_cart_coords_set( geom, coords(1,1,i) ) )
     >      call errquit('zts_meps: set geometry',0,GEOM_ERR)
         call sym_geom_project(geom, 1.0d-6)

         if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     >     call errquit('zts_meps: store geometry',0,GEOM_ERR)

         !*** Skip the gradient for frozen stuff, unless first step - need E ***
         skip1 = freeze1.and.(i.eq.1).and.(nstep.gt.0)
         skip2 = freezeN.and.(i.eq.nbeads).and.(nstep.gt.0)
         if (.not.(skip1.or.skip2)) then

*           *** copy bead movecs to taskmovecs ***
            call util_file_prefix('movecs',perm_name)
            movecs_name = 'bead'//bead_index_name(i)//'.movecs'
            call util_file_name_resolve(perm_name, .false.)
            call util_file_name_resolve(movecs_name, .false.)
            permlen = inp_strlen(perm_name)
            movecslen = inp_strlen(movecs_name)
            if (ga_nodeid().eq.0) then
              inquire(file=movecs_name,exist=status)
              if (status) then
                 call util_file_copy(movecs_name(1:movecslen),
     >                               perm_name(1:permlen))
              end if
            end if

*           *** run gradient task ***
            if (.not.rtdb_get(rtdb,"task:ignore",mt_log,1,ignore))
     >         ignore = .false.
            if(.not.rtdb_put(rtdb,"scf:converged",mt_log,1,.false.))
     >        call errquit('scf:converged put',0,0)
            if(.not.rtdb_put(rtdb,"dft:converged",mt_log,1,.false.))
     >        call errquit('dft:converged put',0,0)

            if (oprint) write(*,*) "string: running bead ",i
            call ga_sync()
*            if ( .not. (task_gradient(rtdb).or.ignore))
            tmpVar = task_gradient(rtdb)
            if ( .not. (tmpVar.or.ignore))
     >         call errquit('zts_meps: gradient failed',0, 0)


*           *** copy taskmovecs to bead movecs ****
            if (ga_nodeid().eq.0) then
               inquire(file=perm_name,exist=status)
               if (status) then
                  call util_file_copy(perm_name(1:permlen),
     >                                movecs_name(1:movecslen))
               end if
            end if

            if (.not. rtdb_get(rtdb,'task:energy',mt_dbl,1,energy(i)))
     &         call errquit('zts_meps: failed getting energy',0,
     &         RTDB_ERR)
            if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, 3*natoms,
     &          dVdx(1,1,i))) call errquit(
     &          'zts_meps: failed getting gradient',0,RTDB_ERR)
            if (oprint) 
     >        write(*,'(A,I3,A,F18.6)') 
     >        "string: finished bead ",i," energy=",energy(i)
         end if
      end do

      if (freeze1) then
        do j = 1, natoms
          do k = 1, 3
            dVdx(k,j,1) = 0.0d0
          end do
        end do
      end if

      if (freezeN) then
        do j = 1, natoms
          do k = 1, 3
            dVdx(k,j,nbeads) = 0.0d0
          end do
        end do
      end if

      !**** Projection of gradient onto perpedicular to path? ****
      if (projection1) then
        do i = 2,nbeads - 1
          do j = 1,natoms
            do k = 1, 3
              g1(k) = coords(k,j,i) - coords(k,j,i+1)
              g2(k) = coords(k,j,i) - coords(k,j,i-1)
            enddo
            call cross_product(g1,g2,g3) ! Vector perpendicular to plane
            tmp = 0.0d0
            do k = 1, 3
              tmp = tmp + g3(k)**2
            enddo
            if (tmp.gt.1.0d-16) then
              tmp = 1.0d0/dsqrt(tmp)
              do k = 1, 3
                g3(k) = g3(k)*tmp ! Normalize
              end do
              ! Projection
              call cross_product(dVdx(1,j,i),g3,g2)
              call cross_product(g3,g2,dVdx(1,j,i))
            endif
          enddo
        enddo
      endif

      return
      end


*     ************************************************
*     *                                              *
*     *               string_move                    *
*     *                                              *
*     ************************************************
      subroutine string_move(alpha,interpol,geom,natoms,nbeads,
     >                       dg,wrk,dcoords,coords1,HdVdx,coords2)
      implicit none
      real*8 alpha
      integer interpol,geom
      integer natoms,nbeads
      real*8 dg(*),wrk(*)
      real*8 dcoords(nbeads)
      real*8 coords1(3,natoms,nbeads)
      real*8 HdVdx(3,natoms,nbeads)
      real*8 coords2(3,natoms,nbeads)

*     **** local variables ****
      integer i,j,k
      logical badgeom

*     **** external functions ****
      real*8   zts_distance
      external zts_distance

      do i=1,nbeads
        do j=1,natoms
          do k=1,3
            coords2(k,j,i) = coords1(k,j,i) - alpha*HdVdx(k,j,i)
          end do
        end do
      end do

!     *** Project out net motion - Can rotate/translater frozen N bead ***
      do i = 2, nbeads
         call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
     >                        geom)
      end do

      !*** reparameterize - This is where the magic happens ****
      dcoords(1) = 0.0d0
      do i=2,nbeads
        dcoords(i) = 
     &     zts_distance(natoms,coords2(1,1,i),coords2(1,1,i-1))
      end do
      do i=2,nbeads
        dcoords(i) = dcoords(i) + dcoords(i-1) ! Now relative to bead 1
      end  do
      do i=1,nbeads
        dcoords(i) = dcoords(i)/dcoords(nbeads) ! Runs from zero to one
      end do
      do i=2,nbeads-1
        if (dcoords(i).le.dcoords(i-1)) then  ! paranoid
            dcoords(i) = (99.0d0*dcoords(i-1)+dcoords(i+1))/100.0d0
        end if
      end do
      if (dcoords(nbeads) .le. dcoords(nbeads-1)) then  ! paranoid
        dcoords(nbeads-1) = (dcoords(nbeads-2)+dcoords(nbeads))/2.0d0
      end if

      ! We now interpolate new points along the string
      ! This adds more points near regions of big dcoords
      ! This really a normalized arc length on a spline
      if (interpol .eq. 1) then
        call zts_linear_fit(natoms*3,nbeads,coords2,dcoords,wrk)
      else if (interpol .eq. 2) then
        call zts_spline_fit(natoms*3,nbeads,coords2,dcoords,dg)
      else if (interpol .eq. 3) then
        call zts_akima_fit(natoms*3,nbeads,coords2,dcoords,dg)
      else
         call errquit('zts_meps: fatal error',0,0)
      end if

!     Project out net motion - Can rotate/translater frozen N bead
      do i=2,nbeads
        call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
     >                      geom)
      end do

c!     Look for atoms that are too close
c      do i=1,nbeads
c        call zts_Robinson(natoms,geom,coords2(1,1,i),badgeom)
c      end do
c
c!      Project out net motion
c      do i=2,nbeads
c        call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
c     >                      geom)
c      end do
      return 
      end

      


*     ************************************************
*     *                                              *
*     *               string_verlet_move             *
*     *                                              *
*     ************************************************
      subroutine string_verlet_move(dt,interpol,geom,natoms,nbeads,
     >                       dg,wrk,mass,dcoords,coords1,dVdx,coords2)
      implicit none
      real*8 dt
      integer interpol,geom
      integer natoms,nbeads
      real*8 dg(*),wrk(*)
      real*8 mass(*)
      real*8 dcoords(nbeads)
      real*8 coords1(3,natoms,nbeads)
      real*8 dVdx(3,natoms,nbeads)
      real*8 coords2(3,natoms,nbeads)

*     **** local variables ****
      integer i,j,k
      logical badgeom
      real*8 v

*     **** external functions ****
      real*8   zts_distance
      external zts_distance

*     *** c1 <- 2*c1 - c0 + dti*g          ***
*     ***    <- c1 + ((c1-c0)/t)*t + dti*g ***
      do i=1,nbeads
        do j=1,natoms
          do k=1,3
            v = coords2(k,j,i) - coords1(k,j,i)
            if (v*dVdx(k,j,i).gt.0.0d0) v = 0.0d0
            coords1(k,j,i)= coords2(k,j,i)
            coords2(k,j,i)= coords2(k,j,i)+v-(dt*dt/mass(j))*dVdx(k,j,i)
          end do
        end do
      end do

!     *** Project out net motion - Can rotate/translater frozen N bead ***
      do i = 2, nbeads
         call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
     >                        geom)
      end do

      !*** reparameterize - This is where the magic happens ****
      dcoords(1) = 0.0d0
      do i=2,nbeads
        dcoords(i) = 
     &     zts_distance(natoms,coords2(1,1,i),coords2(1,1,i-1))
      end do
      do i=2,nbeads
        dcoords(i) = dcoords(i) + dcoords(i-1) ! Now relative to bead 1
      end  do
      do i=1,nbeads
        dcoords(i) = dcoords(i)/dcoords(nbeads) ! Runs from zero to one
      end do
      do i=2,nbeads-1
        if (dcoords(i).le.dcoords(i-1)) then  ! paranoid
            dcoords(i) = (99.0d0*dcoords(i-1)+dcoords(i+1))/100.0d0
        end if
      end do
      if (dcoords(nbeads) .le. dcoords(nbeads-1)) then  ! paranoid
        dcoords(nbeads-1) = (dcoords(nbeads-2)+dcoords(nbeads))/2.0d0
      end if

      ! We now interpolate new points along the string
      ! This adds more points near regions of big dcoords
      ! This really a normalized arc length on a spline
      if (interpol .eq. 1) then
        call zts_linear_fit(natoms*3,nbeads,coords2,dcoords,wrk)
      else if (interpol .eq. 2) then
        call zts_spline_fit(natoms*3,nbeads,coords2,dcoords,dg)
      else if (interpol .eq. 3) then
        call zts_akima_fit(natoms*3,nbeads,coords2,dcoords,dg)
      else
         call errquit('zts_meps: fatal error',0,0)
      end if

!     Project out net motion - Can rotate/translater frozen N bead
      do i=2,nbeads
        call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
     >                      geom)
      end do

c!     Look for atoms that are too close
c      do i=1,nbeads
c        call zts_Robinson(natoms,geom,coords2(1,1,i),badgeom)
c      end do
c
c!      Project out net motion
c      do i=2,nbeads
c        call zts_min_motion(natoms,coords2(1,1,i),coords2(1,1,i-1),
c     >                      geom)
c      end do
      return 
      end

      


