c
c     initialize coordinates and velocities
c
      subroutine qmd_coords_vels(rtdb,geom,nat,istep,ekin,rstfilename,
     $       rstfileunit,r,v,m,lact,nactive,restart,do_namd,state,
     $       nstates,coup,zden)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "global.fh"
#include "geom.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer rtdb               ! rtdb
      integer geom               ! geom
      integer istep              ! step count
      integer nat                ! number of atoms
      integer nactive            ! number of active atoms
      integer rstfileunit        ! restart file unit number
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision ekin      ! kinetic energy
      logical lact(nat)          ! active atoms
      logical restart            ! restart flag
      character(*) rstfilename   ! restart file name
c
c     NAMD variables
      logical do_namd            ! NAMD flag
      integer state              ! current state
      integer nstates            ! number of states
c     THESE ARE NOT ALLOCATED IF do_namd IS FALSE
      double precision coup(nstates,nstates) ! coupling matrix
      double complex zden(nstates,nstates)   ! density matrix
c
c     Local variables
      integer i, j
      integer logsize
      double precision gscale
      logical rst_exists
      
c
      character*30 pname
c
c     preliminaries
      pname='qmd_coords_vels: '
c
      logsize=MA_sizeof(MT_LOG,1,MT_BYTE)
c
c     Check for restart file is restart flag is true
      rst_exists=.false.
      if (restart) then
        if (ga_nodeid().eq.0)
     $    inquire(file=rstfilename,exist=rst_exists)
        call ga_brdcst(Msg_Vec_Stat+MSGINT,rst_exists,logsize,0)
        if (.not.rst_exists) then
          if (ga_nodeid().eq.0) then
            write(luout,*) 'No restart file found'
            write(luout,*) 'Using default values'
          end if
        end if
      end if
c
      restart=rst_exists
c
      if (restart) then
        if (ga_nodeid().eq.0) then
          write(luout,*) 'Restart file found'
          write(luout,*) 'Reading positions, velocities, KE,'
          write(luout,*) 'and last time step from file.'
        end if
c
c       read velocities from file
        call qmd_read_rst(rstfilename,rstfileunit,nat,v,r,ekin,
     $         istep,state,nstates,do_namd,zden,coup)
       
      else
c
c       pull coordinates and velocites from geometry object
        if (.not.geom_coords_vels_get(geom,r,v))
     $    call errquit(pname//'geom_coords_vels_get',0,GEOM_ERR)
c
c       determine if velocities were present in input file by looking
c        for the largest element of velocities array
        if (maxval(dabs(v)).lt.1.d-8) then
         if (ga_nodeid().eq.0) write(*,*) 'no velocities found in input'
c
c        randomize velocities for the given temperature
         call qmd_randomize(nat,m,v)
        end if
c
c       starting from beginning
        istep=1
      end if ! restart
c
c     remove velocities from inactive atoms
      call qmd_freeze_inactive(nat,m,v,lact)
c
c     remove global translations and rotations
c     kinetic energy is the same on exit as was on entry
c     don't do it if we have inactive atoms
      if (nat.eq.nactive) call qmd_fix_momentum(nat,m,v,r)
c
c     persist updated coordinates in geom
      if (.not.geom_cart_coords_set(geom,r))
     $   call errquit(pname//'geom_cart_coords_set',0,GEOM_ERR)
c
c     persist updated velocity in geom
      if (.not.geom_vel_set(geom,v))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
c     store updated geometry
      if(.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &   call errquit(pname//'geom_rtdb_store',0,GEOM_ERR)
c
      return
      end
