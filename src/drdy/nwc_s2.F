      Subroutine nwc_s2 (x,v,dx,fmat,hess,xs,amass,natom,n3,n3tm,ndbar)
C
C  Run NWChem task to get potential and first and second derivatives
C   for geometry x.
C
      Implicit None
#include "nwc_drdyP.fh"
#include "geom.fh"
#include "printlevels.fh"
#include "stdio.fh"
      Integer natom, n3, n3tm, ndbar
      Double Precision v
      Double Precision x(n3),dx(n3),fmat(n3tm,n3),
     &      hess(ndbar),xs(n3),amass(n3)
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
      Integer current_print_level
      logical copy_sym, copy_c1
C 
C Transform to cartesian coordinates
*      write(6,*)' x before ',x
      Call drdy_trans (2,n3,amass,x,dx)
*
      call util_print_get_level(current_print_level)
      call util_print_set_level(print_none)
*
* ... delete finite difference files that may exist.
      call util_file_name('hess',  .false., .false.,dummy_file)
      call util_file_unlink(dummy_file)
      call util_file_name('fd_ddipole',  .false., .false.,dummy_file)
      call util_file_unlink(dummy_file)
*
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('nwc_s2: geom_create failed',911)
      if (.not.geom_set_user_units(geom,'a.u.'))
     &      call errquit('nwc_s2: geom_set_user_units failed',911)
      if (.not.geom_cart_set(geom,natom,atomic_labels,x,atomic_charge))
     &      call errquit('nwc_s2: geom_cart_set failed',911)
*      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
*     &      call errquit('nwc_s2: geom_rtdb_store failed',911)
      ncenter_B4_autosym = natom
      ncenter = natom
      threquiv = 0.01d00
      call dcopy(n3,x,1,dx,1)   ! temporary use of dx 
      call dcopy(ncenter,atomic_charge,1,copy_charge,1)
      do i = 1,ncenter
        copy_labels(i) = atomic_labels(i)
      enddo
      call geom_auto_sym(my_rtdb,geom,dx,
     &      copy_charge,copy_labels,ncenter,
     &      threquiv,groupname)
      copy_sym = .false.
      copy_c1 = .false.
      if (geom_group_set(geom,groupname)) then
        copy_c1 = ncenter_B4_autosym .ne. ncenter
        if (.not.copy_c1) copy_sym = .true.
      else
        copy_c1 = .true.
      endif
      if (copy_c1.or.copy_sym) then
        if      (copy_c1) then
          write(luout,*)' autosym failed using original coordinates'
          if (.not.geom_cart_set(geom,ncenter,atomic_labels,x,
     &          atomic_charge))
     &          call errquit('nwc_s2: geom_cart_set failed',911)
        else if (copy_sym) then
          if (.not.geom_cart_set(geom,ncenter,copy_labels,dx,
     &          copy_charge))
     &          call errquit('nwc_s2: geom_cart_set failed',911)
          nata = ncenter
          call sym_nwc(geom,my_rtdb,nata,.false.,1.0d00,threquiv)
          if (.not.geom_ncent(geom,ncenter)) call errquit
     &          ('nwc_s2:failed geom_ncent',911)
          if (ncenter_B4_autosym.ne.ncenter) call errquit
     &          ('nwc_s2: sym_nwc failed ',911)
          call sym_geom_project(geom, threquiv)
*          write(luout,*)' autosym worked'
        endif
      else
        call errquit('nwc_s2: logic error in copy_{c1|sym}',911)
      endif
      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
     &      call errquit('nwc_s2: geom_rtdb_store failed',911)
      if (.not.geom_destroy(geom))
     &      call errquit('nwc_s2: geom_destroy failed',911)
      call util_print_set_level(current_print_level)
*      write(6,*)' x after  ',x
      call drdy_synch(DRDY_CODE_HESSIAN,'nwc_s2')
      call drdy_nwc_get_energy(my_rtdb,v)
      call drdy_nwc_get_coords(my_rtdb,xs)
      call drdy_nwc_get_gradient(my_rtdb,dx)
      call drdy_nwc_get_hessian(my_rtdb,hess)
*      
C
C  Get transformation matrix (rotation matrix amat) from the gaussian
C     standard orientation (xs) to the original orientation (x)
      call drdy_rotmat(x,xs,amat,rr,rrs,rri,natom,n3,ndima)
      call drdy_rotg(dx,amat,rr,natom,n3,ndima)
C
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
      call drdy_rotf(fmat,amat,rr,natom,n3,ndima,n3tm)
C  Transform coordinates and gradient to mass-weight cartesians
      Call drdy_trans (1,n3,amass,x,dx)
C  Transform hessian matrix to mass-weighted coordinates
      Call drdy_tranf (1,n3tm,n3,amass,fmat)
*
      return
1000  Format (1x,79a1)
1001  Format (' ')
1002  Format (1x,2i2)
1003  Format (1x,5a1,1p3e20.10)
600   format (1x,' construct gaussian data file:')
      End
