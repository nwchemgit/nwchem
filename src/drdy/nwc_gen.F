C $Id$
      Subroutine nwc_gen (x,v,dx,fmat,hess,xs,natom,n3,n3tm,ndbar,op)
C
C  Run NWChem task to get potential and first and second derivatives
C   for geometry x.
C
      Implicit None
#include "errquit.fh"
#include "nwc_drdyP.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "printlevels.fh"
#include "rtdb.fh"
#include "stdio.fh"
      Integer natom, n3, n3tm, ndbar, op
      Double Precision v
      Double Precision x(n3),dx(n3),fmat(n3tm,n3),
     &      hess(ndbar),xs(n3), velocities(3,n3)
*
      Integer ndima
      Parameter (ndima=3)
      Double Precision amat(ndima,ndima),
     &      rr(ndima,ndima),rrs(ndima,ndima),
     *   rri(ndima,ndima)
*
      Integer geom
      Integer i,j,jj
      double precision threquiv
      character*16 groupname
      character*255 dummy_file
      Integer ncenter_B4_autosym, ncenter, nata
      Integer current_print_level,nops
      logical oautosym
      logical status_rtdb, status_ignore
*     logical copy_sym, copy_c1, oautosym
*

      call util_print_get_level(current_print_level)
      call util_print_set_level(print_none)
*
      if (op.eq.DRDY_CODE_HESSIAN) then
*
* ... delete finite difference files that may exist.
        call util_file_name('hess',  .false., .false.,dummy_file)
        call util_file_unlink(dummy_file)
        call util_file_name('fd_ddipole',.false.,.false.,dummy_file)
        call util_file_unlink(dummy_file)
      endif
*
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('nwc_gen: geom_create failed',911,
     &       GEOM_ERR)
      if (.not.geom_set_user_units(geom,'a.u.'))
     &      call errquit('nwc_gen: geom_set_user_units failed',911,
     &       GEOM_ERR)
      if (.not.geom_cart_set(geom,natom,atomic_labels,x,atomic_charge))
     &      call errquit('nwc_gen: geom_cart_set failed',911,
     &       GEOM_ERR)
      if (.not.geom_masses_set(geom,natom,nwcmass))
     &      call errquit('nwc_gen:geom_masses_set failed',911,
     &       GEOM_ERR)
      if (.not.geom_atomct_set(geom,natom,atomic_atomct))
     &      call errquit('nwc_gen:geom_atomct_set failed',911,
     &       GEOM_ERR)
*      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
*     &      call errquit('nwc_gen: geom_rtdb_store failed',911,
*     &       RTDB_ERR)
c
c Now figure out if we are going to use autosym and if so, get the
c value for the threshold and take care of the symmetry
c
      status_rtdb   = rtdb_parallel(.false.)

      if (.not. rtdb_get(my_rtdb,'drdy:autosym',mt_log,1,oautosym))
     &  call errquit('nwc_gen: problem getting autosym from rtdb',555,
     &       RTDB_ERR)
      if (oautosym) then
        if (.not.rtdb_get(my_rtdb,'drdy:threquiv',mt_dbl,1,threquiv))
     &   call errquit('nwc_gen:problem getting threquiv from rtdb',555,
     &       RTDB_ERR)
        ncenter_B4_autosym = natom
        ncenter = natom
        call dcopy(n3,x,1,dx,1)   ! temporary use of dx 
        call dcopy(ncenter,atomic_charge,1,copy_charge,1)
        call dcopy(ncenter,atomic_atomct,1,copy_atomct,1)
        do i = 1,ncenter
          copy_labels(i) = atomic_labels(i)
        enddo
        call dcopy(3*n3,0.d0,0,velocities,1)
        call geom_auto_sym(my_rtdb,geom,dx,
     &      copy_charge,copy_labels,copy_atomct,ncenter,
     &      threquiv,groupname,velocities)
        if (op.ne.DRDY_CODE_SPENERGY)
     &    write(luout,10) groupname
        if (groupname(1:2).ne."C1") then
        if (geom_group_set(geom,groupname)) then
          if (.not.geom_cart_set(geom,ncenter,copy_labels,dx,
     &         copy_charge))
     &      call errquit('nwc_gen: geom_cart_set failed', 0,
     &       GEOM_ERR)
          if (ncenter_B4_autosym .ne. ncenter) call errquit
     &       ('nwc_gen: autosym bug : number of atoms wrong',
     &         ncenter, INPUT_ERR)
        else
          write(luout,*) ' autosym detected unknown group ',
     &                      groupname
          call errquit('nwc_gen: autosym: invalid group',0, INPUT_ERR)
        endif
c
c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c
        nata = ncenter
        call sym_nwc(geom,my_rtdb,nata,.false.,1.0d00,threquiv,nops,
     %       .false.)
c
c     Check that if we used autosym that we ended up with the
c     same no. of atoms ... if we don't then autosym and nwchemsym
c     don't agree on the orientation of point group elements
c
        if (.not. geom_ncent(geom,ncenter)) call errquit
     $       ('nwc_gen: geom_cent?',0, GEOM_ERR)
        if (ncenter_B4_autosym .ne. ncenter) call errquit
     $       ('nwc_gen: autosym bug : too many atoms',ncenter,
     &       INPUT_ERR)
c
c     Force exact symetry on the coordinates
c
        call sym_geom_project(geom, threquiv)
        endif
      endif
c
c store as default geometry
c
      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
     &      call errquit('nwc_gen: geom_rtdb_store failed',911,
     &       RTDB_ERR)
      if (.not.geom_destroy(geom))
     &      call errquit('nwc_gen: geom_destroy failed',911, GEOM_ERR)
      call util_print_set_level(current_print_level)
*      write(luout,*)' x after  ',x
c
c Actually do the deed
c
      status_ignore   = rtdb_parallel(status_rtdb)
      status_rtdb   = rtdb_parallel(.false.)

      if ((op.ge.DRDY_CODE_SPENERGY).and.
     &    (op.le.DRDY_CODE_HESSIAN)) then
        call drdy_synch(op,'nwc_gen')
      else
        write(luout,*) 'Trying to calculate an unknown op code'
        call errquit('nwc_gen: Unknown op code',op, INPUT_ERR)
      endif
      status_ignore   = rtdb_parallel(status_rtdb)
      status_rtdb   = rtdb_parallel(.false.)
c
c Get the appropriate values to hand back
c
      call drdy_nwc_get_energy(my_rtdb,v)
      if (op.ge.DRDY_CODE_GRADIENT) then
         call drdy_nwc_get_coords(my_rtdb,xs)
         call drdy_nwc_get_gradient(my_rtdb,dx)
      endif
      if (op.eq.DRDY_CODE_HESSIAN) 
     &   call drdy_nwc_get_hessian(my_rtdb,hess)
*      
C
C  Get transformation matrix (rotation matrix amat) from the gaussian
C     standard orientation (xs) to the original orientation (x)
      if (op.ge.DRDY_CODE_GRADIENT) then
        call drdy_rotmat(x,xs,amat,rr,rrs,rri,natom,ndima)
        call drdy_rotg(dx,amat,rr,natom,n3)
      endif
C
      if (op.eq.DRDY_CODE_HESSIAN) then
        jj = 0
        do i = 1,n3
          do j = 1,i
            jj = jj + 1
            fmat(i,j) = hess(jj)
            fmat(j,i) = hess(jj)
          enddo
        enddo
C  Transform hessian matrix from gaussian standard orientation to the
C     original orientation
        call drdy_rotf(fmat,amat,rr,natom,n3,n3tm)
      endif
*
      status_ignore = rtdb_parallel(status_rtdb)
      return
10    Format (1x,'Symmetry group is ',a3)
1000  Format (1x,79a1)
1001  Format (' ')
1002  Format (1x,2i2)
1003  Format (1x,5a1,1p3e20.10)
      End
