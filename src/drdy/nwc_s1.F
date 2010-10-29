C $Id$
      Subroutine nwc_s1 (x,v,dx,xs,amass,natom,n3)
C THIS ROUTINE IS NOW DEAD.  LOOK AT NWC_GEN INSTEAD.
C
C  Run NWChem task to get potential and first derivatives
C   for geometry x.
C
      Implicit None
#include "errquit.fh"
#include "nwc_drdyP.fh"
#include "geom.fh"
#include "stdio.fh"
#include "printlevels.fh"
      Integer natom, n3
      Double Precision v
      Double Precision x(n3),dx(n3),xs(n3),amass(n3)
*
      Integer ndima
      Parameter (ndima=3)
      Double Precision amat(ndima,ndima),
     &      rr(ndima,ndima),rrs(ndima,ndima),
     *   rri(ndima,ndima)
*
      Integer geom
      double precision threquiv
      character*16 groupname
      Integer i, ncenter_B4_autosym, ncenter, nata
      Integer current_print_level
      logical copy_sym, copy_c1
*
      call util_print_get_level(current_print_level)
      call util_print_set_level(print_none)
*      
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('nwc_s1: geom_create failed',911, GEOM_ERR)
      if (.not.geom_set_user_units(geom,'a.u.'))
     &      call errquit('nwc_s1: geom_set_user_units failed',911,
     &       GEOM_ERR)
      if (.not.geom_cart_set(geom,natom,atomic_labels,
     &      x,atomic_charge))
     &      call errquit('nwc_s1: geom_cart_set failed',911,
     &       GEOM_ERR)
*      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
*     &      call errquit('nwc_s1: geom_rtdb_store failed',911,
*     &       RTDB_ERR)
      ncenter_B4_autosym = natom
      ncenter = natom
      threquiv = 0.01d00
      call dcopy(n3,x,1,dx,1)  ! temporary use of dx 
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
     &          call errquit('nwc_s1: geom_cart_set failed',911,
     &       GEOM_ERR)
        else if (copy_sym) then
          if (.not.geom_cart_set(geom,ncenter,copy_labels,dx,
     &          copy_charge))
     &          call errquit('nwc_s1: geom_cart_set failed',911,
     &       GEOM_ERR)
          nata = ncenter
          call sym_nwc(geom,my_rtdb,nata,.false.,1.0d00,threquiv)
          if (.not.geom_ncent(geom,ncenter)) call errquit
     &          ('nwc_s1:failed geom_ncent',911, GEOM_ERR)
          if (ncenter_B4_autosym.ne.ncenter) call errquit
     &          ('nwc_s1: sym_nwc failed ',911, GEOM_ERR)
          call sym_geom_project(geom, threquiv)
*          write(luout,*)' autosym worked'
        endif
      else
        call errquit('nwc_s1: logic error in copy_{c1|sym}',911,
     &       INPUT_ERR)
      endif
      if (.not.geom_rtdb_store(my_rtdb,geom,'geometry'))
     &      call errquit('nwc_s1: geom_rtdb_store failed',911,
     &       RTDB_ERR)
      if (.not.geom_destroy(geom))
     &      call errquit('nwc_s1: geom_destroy failed',911,
     &       GEOM_ERR)
      call util_print_set_level(current_print_level)
*      write(luout,*)' x after  ',x
      call drdy_synch(DRDY_CODE_GRADIENT,'nwc_s1')
      call drdy_nwc_get_energy(my_rtdb,v)
      call drdy_nwc_get_coords(my_rtdb,xs)
      call drdy_nwc_get_gradient(my_rtdb,dx)
*      
C
C  Get transformation matrix (rotation matrix amat) from the gaussian
C     standard orientation (xs) to the original orientation (x)
      call drdy_rotmat(x,xs,amat,rr,rrs,rri,natom,n3,ndima)
      call drdy_rotg(dx,amat,rr,natom,n3,ndima)
*
      return
1000  Format (1x,79a1)
1001  Format (' ')
1002  Format (1x,2i2)
1003  Format (1x,5a1,1p3e20.10)
600   format (1x,' construct gaussian data file:')
      End
