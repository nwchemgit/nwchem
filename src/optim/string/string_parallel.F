*     ****************************************
*     *                                      *
*     *             zts_runall_parallel      *
*     *                                      *
*     ****************************************
!>
!> @brief Run all beads concurrently
!>
!> Rather than looping over all beads this subroutine runs
!> all beads at the same time. To do so the GA processor
!> group is split into as many subgroups as there are beads.
!> In every subgroup an appropriate copy of the RTDB with
!> the right geometry and the movecs are created. Subquently
!> every subgroup calls the gradient module on its bead.
!>
!> @param rtdb        global runtime database file handle
!> @param geom        global geometry handle
!> @param natoms      the number of atoms in the geometry
!> @param nbeads      the number of beads in the string
!> @param nstep       the number of the current step
!> @param freeze1     shall the structure of the 1st bead be frozen?
!> @param freezeN     shall the structure of the last bead be frozen?
!> @param projection1 shall the gradient be projected on the
!>                    perpendicular of the path?
!> @param coords      the atomic coordinates of all beads
!> @param energy      the total energies for all beads
!> @param dVdx        the atomic gradients of all beads
!>
      subroutine zts_runall_parallel(
     >                      rtdb,geom,natoms,nbeads,nstep,
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
#include "stdio.fh"
#include "bgj_common.fh"

*     **** local variables ****
      logical skip1,skip2,status,oprint,ignore
      integer i,j,k
      integer permlen,movecslen
      character*255 movecs_name,perm_name
      real*8 g1(3),g2(3),g3(3),tmp
      integer nprocs_per_bead
      integer nprocs_remainder
      integer nprocs_for_bead(nbeads)
      character*255 global_rtdb_fname
      integer glbl_len
      character*255 bead_rtdb_fname
      character*255 file_prefix
      character*255 bead_out_fname
      integer bead_len, prefix_len
      character*255 movecs_out
      character*255 bead_movecs_out
      integer bead_movecs_len
      integer global_pgroup, global_iproc, global_nproc
      integer bead_pgroup, bead_iproc, bead_nproc
      integer ibead_lo, ibead_hi, mybead, imiddle_bead
      integer nactive_beads
      integer ilow, ihigh, ii
      integer np, ip, iplo, iphi, iproc, nproc
      integer nmygrpproc
      character*132 line
      integer ifile, ierr
      real*8  tmp_energy
      real*8, allocatable :: tmp_dVdx(:,:)
      integer, allocatable :: iproc_to_bead(:)
      integer, allocatable :: mygrpprocs(:)

*     **** external functions ****
      logical  task_gradient
      external task_gradient
      character*7 bead_index_name
      external    bead_index_name
#ifndef GANXTVAL
      logical, external :: util_statldb
#endif

      oprint = (ga_nodeid().eq.0)
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      allocate(iproc_to_bead(0:nproc-1))
      allocate(mygrpprocs(0:nproc-1))
      allocate(tmp_dVdx(1:3,1:natoms))
      if (nstep.eq.0) then
#ifndef GANXTVAL
        !*** Turn dynamics load balancing off, never turn back on!
        if (util_statldb(.true.,rtdb)) then
          if (iproc.eq.0) 
     &      write(luout,*)' Switched to static load balancing'
        endif
#endif
        !*** Do not use disk for DFT scratch files, never turn back on!
        if(.not.rtdb_put(rtdb,'dft:largenode',mt_log,1,.true.))
     &    call errquit('zts_runall_parallel: rtdb_put nodisk failed',
     &                 1,RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'dft:noio', mt_int, 1, 1))
     &    call errquit('zts_runall_parallel: rtdb_put noio failed',
     &                 1,RTDB_ERR)
        if(.not.rtdb_put(rtdb,'dft:xcreplicated',mt_log,1,.false.))
     &    call errquit('zts_runall_parallel: rtdb_put xcreplicated '
     &                 //'fail',1,RTDB_ERR)
        !*** Switch to direct to avoid IO for integrals
        call fock_force_direct(rtdb)
        if (.not. rtdb_put(rtdb, 'int2e:filesize', MT_INT, 1, -1))
     &    call errquit("zts_runall_parallel: rtdb put "
     &                 //"int2e:filesize failed",0,RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'int2e:memsize', MT_INT, 1, -1))
     &    call errquit("zts_runall_parallel: rtdb put int2e:memsize "
     &                 //"failed",0,RTDB_ERR)
      endif

      !*** Calculate the energy and forces for each bead ***
      call dcopy(3*natoms*nbeads,0.0d0,0,dVdx,1)
      !*** Skip the gradient for frozen stuff, unless first step - need E ***
      !*** Node 0 needs to keep first and last energy for global sum      ***
      skip1 = freeze1.and.(nstep.gt.0)
      skip2 = freezeN.and.(nstep.gt.0)
      ibead_lo = 1
      ibead_hi = nbeads
      if (iproc.eq.0) then
        if (skip1) then
          ibead_lo = 2
        endif
        if (skip2) then
          ibead_hi = nbeads-1
        endif
      endif
      do i=ibead_lo, ibead_hi
        energy(i) = 0.0d0
      enddo
      nactive_beads=nbeads
      if (skip1) nactive_beads=nactive_beads-1
      if (skip2) nactive_beads=nactive_beads-1
      nprocs_per_bead=nproc/nactive_beads
      nprocs_remainder=mod(nproc,nactive_beads)
      !*** All processors need to agree on who works on what
      ibead_lo = 1
      ibead_hi = nbeads
      if (skip1) then
        ibead_lo = 2
      endif
      if (skip2) then
        ibead_hi = nbeads-1
      endif
      nprocs_for_bead=0
      do i = ibead_lo, ibead_hi
        nprocs_for_bead(i) = nprocs_per_bead
      enddo
      if (nprocs_remainder.gt.0) then
        !*** transition state SCF typically converges slower
        !*** add remainder procs to the middle of the string
        imiddle_bead = nbeads/2
        if (mod(nbeads,2).ne.0) imiddle_bead = imiddle_bead + 1
        ilow  = imiddle_bead - nprocs_remainder/2
        ihigh = imiddle_bead + nprocs_remainder/2
        if (mod(nprocs_remainder,2).eq.0) ilow = ilow + 1
        do i = ilow, ihigh
          nprocs_for_bead(i) = nprocs_for_bead(i) + 1
        enddo
      endif
      !*** Check we did this right
      np = 0
      do i=1,nbeads
        np = np + nprocs_for_bead(i)
      enddo
      if (np.ne.nproc) then
        call errquit('zts_runall_parallel: nprocs_for_bead '//
     &               'miscalculated',np,UERR)
      endif
      !*** Build table for global proc-id to bead to work on
      ip = 0
      do i=1,nbeads
        iplo = ip
        iphi = ip + nprocs_for_bead(i)-1
        do ii = iplo, iphi
          iproc_to_bead(ii) = i
        enddo
        ip = ip + nprocs_for_bead(i)
      enddo
      mybead = iproc_to_bead(iproc)
      nmygrpproc = 0
      do ip = 0, nproc-1
        if (mybead.eq.iproc_to_bead(ip)) then
          mygrpprocs(nmygrpproc) = ip
          nmygrpproc = nmygrpproc + 1
        endif
      enddo
      !*** Switching to the bead local group
      global_pgroup = ga_pgroup_get_default()
      global_iproc  = iproc
      global_nproc  = nproc
      bead_pgroup   = ga_pgroup_create(mygrpprocs,nmygrpproc)
      bead_iproc    = ga_pgroup_nodeid(bead_pgroup)
      bead_nproc    = ga_pgroup_nnodes(bead_pgroup)
      !*** Each bead will need to have its own RTDB with its own
      !*** geometry and movecs file
      !flush(rtdb)
      if (.not. rtdb_getfname(rtdb,global_rtdb_fname)) 
     &  call errquit('zts_runall_parallel: rtdb_getfname failed',
     &               0,RTDB_ERR)
      if (.not.rtdb_close(rtdb,'keep')) then
        call errquit("zts_runall_parallel: closing global "//
     &                "rtdb failed",0,RTDB_ERR)
      endif
      if (nstep.eq.0) then
        !*** We need to replicate the RTDB for every bead
        do i = 1, nbeads
          glbl_len = inp_strlen(global_rtdb_fname)
          bead_rtdb_fname = global_rtdb_fname(1:(glbl_len-3))
          bead_rtdb_fname = bead_rtdb_fname(1:(glbl_len-3))//
     &                      bead_index_name(i)//'.db'
          bead_len = inp_strlen(bead_rtdb_fname)
          if (global_iproc.eq.0) then
            inquire(file=global_rtdb_fname,exist=status)
            if (status) then
              call util_file_copy(global_rtdb_fname(1:glbl_len),
     &                            bead_rtdb_fname(1:bead_len))
            endif
          endif
        enddo
      endif
      call ga_sync()
      !*** Switch process group
      call ga_pgroup_set_default(bead_pgroup)
      bead_iproc = ga_nodeid()
      bead_nproc = ga_nnodes()
      !*** Redirect standard output to additional files
      if (global_iproc.ne.0) then
        if (bead_iproc.eq.0) then
          glbl_len = inp_strlen(global_rtdb_fname)
          bead_out_fname = global_rtdb_fname(1:(glbl_len-3))
          bead_out_fname = bead_out_fname(1:(glbl_len-3))//
     &                     bead_index_name(mybead)//'.out'
#ifdef USE_SUBGROUPS
          LuOut = 6
#endif
          open(UNIT=LuOut,FILE=bead_out_fname)
        endif
      endif
      glbl_len = inp_strlen(global_rtdb_fname)
      file_prefix = global_rtdb_fname(1:(glbl_len-3))//
     &              bead_index_name(mybead)
      prefix_len = inp_strlen(file_prefix)
      bead_rtdb_fname = file_prefix(1:prefix_len)//'.db'
      bead_len = inp_strlen(bead_rtdb_fname)
      if (.not.rtdb_open(bead_rtdb_fname(1:bead_len), 'old', rtdb))
     &  call errquit('zts_runall_parallel: opening bead rtdb '//
     &               'file failed',0,RTDB_ERR)
      bgj_rtdb = rtdb
      call util_file_prefix_set(file_prefix)
      if (.not.rtdb_put(rtdb,"scf:converged",mt_log,1,.false.))
     >  call errquit('zts_runall_parallel: scf:converged put',0,0)
      if (.not.rtdb_put(rtdb,"dft:converged",mt_log,1,.false.))
     >  call errquit('zts_runall_parallel: dft:converged put',0,0)

      if (nstep.eq.0) then
        !*** We to update the file prefix on the bead RTDB file
        call util_file_info_rtdb(rtdb)
        !*** We need to reset the output movecs file
        if (rtdb_cget(rtdb,'scf:output vectors',1,movecs_out)) then
          bead_movecs_len = inp_strlen(movecs_out)-7
          bead_movecs_out = movecs_out(1:bead_movecs_len)//
     >                      bead_index_name(mybead)//'.movecs'
        else
          bead_movecs_out = bead_rtdb_fname(1:bead_len-3)//
     >                      '.movecs'
        endif
        if (.not.rtdb_cput(rtdb,'scf:output vectors',1,
     >                     bead_movecs_out)) then
          call errquit('zts_runall_parallel: scf:output vectors put',
     &                 0,0)
        endif
        if (rtdb_cget(rtdb,'dft:output vectors',1,movecs_out)) then
          bead_movecs_len = inp_strlen(movecs_out)-7
          bead_movecs_out = movecs_out(1:bead_movecs_len)//
     >                      bead_index_name(mybead)//'.movecs'
        else
          bead_movecs_out = bead_rtdb_fname(1:bead_len-3)//
     >                      '.movecs'
        endif
        if (.not.rtdb_cput(rtdb,'dft:output vectors',1,
     &                     bead_movecs_out)) then
          call errquit('zts_runall_parallel: dft:output vectors put',
     &                 0,0)
        endif
        !*** We want only one bead to print the basis set information
        !*** and suppress that output for all other beads
        if (mybead.ne.1) then
          if(.not.rtdb_put(rtdb,"basisprint:ao basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:ao basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:cd basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:cd basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:xc basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:xc basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:ecp basis",mt_log,1,
     &                     .false.))
     &      call errquit("zts_runall_parallel: 'basisprint:ecp basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:ri-mp2 basis",mt_log,1,
     &                     .false.))
     &      call errquit("zts_runall_parallel: "
     &                 //"'basisprint:ri-mp2 basis'"
     &                 //" put",0,RTDB_ERR)
        endif
      else if (nstep.eq.1) then
        !*** We need to reset the input movecs file
        if (rtdb_cget(rtdb,'scf:output vectors',1,movecs_out)) then
          if (.not.rtdb_cput(rtdb,'scf:input vectors',1,
     >                       movecs_out)) then
            call errquit('zts_runall_parallel: scf:input vectors put',
     &                   0,0)
          endif
        endif
        if (rtdb_cget(rtdb,'dft:output vectors',1,movecs_out)) then
          if (.not.rtdb_cput(rtdb,'dft:input vectors',1,
     >                       movecs_out)) then
            call errquit('zts_runall_parallel: dft:input vectors put',
     &                   0,0)
          endif
        endif 
        !*** Suppress the basis set output now also for the first bead
        !*** (we have already printed this information in step 0).
        if (mybead.eq.1) then
          if(.not.rtdb_put(rtdb,"basisprint:ao basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:ao basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:cd basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:cd basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:xc basis",mt_log,1,.false.))
     &      call errquit("zts_runall_parallel: 'basisprint:xc basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:ecp basis",mt_log,1,
     &                     .false.))
     &      call errquit("zts_runall_parallel: 'basisprint:ecp basis'"
     &                 //" put",0,RTDB_ERR)
          if(.not.rtdb_put(rtdb,"basisprint:ri-mp2 basis",mt_log,1,
     &                     .false.))
     &      call errquit("zts_runall_parallel: "
     &                 //"'basisprint:ri-mp2 basis'"
     &                 //" put",0,RTDB_ERR)
        endif
      endif
      !*** Update the geometry 
      if ( .not. geom_cart_coords_set( geom, coords(1,1,mybead) ) )
     >   call errquit('zts_runall_parallel: set geometry',0,GEOM_ERR)
      call sym_geom_project(geom, 1.0d-6)
      if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     >  call errquit('zts_runall_parallel: store geometry',0,GEOM_ERR)
      !*** Run the SCF gradient calculation
      if (.not.rtdb_get(rtdb,"task:ignore",mt_log,1,ignore))
     >   ignore = .false.
      if(.not.rtdb_put(rtdb,"scf:converged",mt_log,1,.false.))
     >  call errquit('zts_runall_parallel: scf:converged put',0,0)
      if(.not.rtdb_put(rtdb,"dft:converged",mt_log,1,.false.))
     >  call errquit('zts_runall_parallel: dft:converged put',0,0)

      if (bead_iproc.eq.0) write(LuOut,*)"string: running bead ",mybead
      call ga_sync()
      tmpVar = task_gradient(rtdb)
      if (.not.(tmpVar.or.ignore))
     &  call errquit("zts_runall_parallel: gradient failed",0,UERR)

      if (bead_iproc.eq.0) then
        if (.not. rtdb_get(rtdb,'task:energy',mt_dbl,1,energy(mybead)))
     &      call errquit('zts_runall_parallel: failed getting energy',0,
     &      RTDB_ERR)
        write(LuOut,'(A,I3,A,F18.6)')
     &        " string: finished bead ",mybead," energy=",energy(mybead)
        if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, 3*natoms,
     &      dVdx(1,1,mybead))) call errquit(
     &      'zts_runall_parallel: failed getting gradient',0,RTDB_ERR)
      else
        if (.not. rtdb_get(rtdb,'task:energy',mt_dbl,1,tmp_energy))
     &      call errquit('zts_runall_parallel: failed getting energy',0,
     &      RTDB_ERR)
        if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, 3*natoms,
     &      tmp_dVdx)) call errquit(
     &      'zts_runall_parallel: failed getting gradient',0,RTDB_ERR)
      endif

      !*** Need to go from processor group to global
      if (.not.rtdb_close(rtdb,'keep'))
     &  call errquit('zts_runall_parallel close failed for pgroup',
     &               0,RTDB_ERR)
      if (global_iproc.ne.0) then
        if (bead_iproc.eq.0) then
          close(UNIT=LuOut)
#ifdef USE_SUBGROUPS
          LuOut = 6
#endif
        endif
      endif
      call ga_sync()
      call ga_pgroup_set_default(global_pgroup)
      if (.not.ga_pgroup_destroy(bead_pgroup))
     &  call errquit('zts_runall_parallel: destroying bead pgroup '//
     &               'failed',0,GA_ERR)
      if (.not.rtdb_open(global_rtdb_fname(1:glbl_len), 'old', rtdb))
     &  call errquit('zts_runall_parallel: opening global rtdb file '//
     &               'failed',0,RTDB_ERR)
      bgj_rtdb = rtdb
      !*** Go back to the global file prefix
      if(.not. rtdb_cget(rtdb,'file_prefix',1,file_prefix))
     &     call errquit('zts_runall_parallel: rtdb get file_prefix '
     &                  //'failed',0,0)
      call util_file_prefix_set(file_prefix)
      !*** Echo stdout from other processor groups to the global stdout
      if (ga_nodeid().eq.0) then
        if (.not.util_get_io_unit(ifile)) then
          call errquit('zts_runall_parallel: out of free file units',
     &                 0,UERR)
        endif
        do ii = ibead_lo+1, ibead_hi
          glbl_len = inp_strlen(global_rtdb_fname)
          bead_out_fname = global_rtdb_fname(1:(glbl_len-3))
          bead_out_fname = bead_out_fname(1:(glbl_len-3))//
     &                   bead_index_name(ii)//'.out'
          open(UNIT=ifile,FILE=bead_out_fname)
          do while(.true.)
            line=""
            read(ifile,'(a132)',ERR=1000,IOSTAT=ierr)line
            write(luout,'(a)')line(1:inp_strlen(line))
          enddo
 1000     close(UNIT=ifile)
        enddo
      endif

      call ga_dgop(202312,energy,nbeads,'+')
      call ga_dgop(202313,dVdx,3*natoms*nbeads,'+')
      
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
      deallocate(tmp_dVdx)
      deallocate(iproc_to_bead)
      deallocate(mygrpprocs)

      return
      end
