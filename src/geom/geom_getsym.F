C> \ingroup geom
C> @{
      subroutine geom_getsym0(rtdb,geom,oautosym,
     ,     ncenter,isys,
     ,     name,tags,
     ,     threquiv,scale,coords,charge,velocities,atomct)
C     $Id$
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "geom.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      character*(*) name,tags(*)
      logical oautosym
      integer rtdb,geom
      integer isys
      integer ncenter
      double precision scale,coords(3,*),charge(*),
     ,     velocities(*),threquiv,atomct(*)
c
      logical oprint,oprint_sym
      integer ncenter_before_auto
      integer nata              ! num. of unique atoms
      integer itries
      logical g_loaded
      integer k_coordsym,l_coordsym,nops
      character*16 groupname
c
      oprint=ga_nodeid().eq.0
      oprint_sym=util_print('autosym symmetry information',
     ,     print_high).and.oprint
      itries=0
      g_loaded=.false.
 2020    continue
      if(oautosym) then
c     
c     ----- detect symmetry -----
c     
         if(.not.g_loaded) then
         if (.not. geom_cart_get(geom,ncenter,tags,coords,charge))
     $           call errquit('geom_getsym0: geom_cart_get failed', 0,
     &       GEOM_ERR)
            if (.not. geom_vel_get(geom, velocities))
     $        call errquit('geom_getsym0: geom_vel_get failed', 0,
     &           GEOM_ERR)
            endif
            g_loaded=.true.
            ncenter_before_auto = ncenter
            call geom_auto_sym(rtdb,geom,coords,charge,tags,atomct,
     $           ncenter,threquiv,groupname,velocities)
            if (geom_group_set(geom,groupname)) then
               if (.not.geom_cart_set(geom,ncenter,tags,coords,charge))
     $             call errquit('geom_getsym0: geom_cart_set failed', 0,
     &       GEOM_ERR)
               if (.not. geom_vel_set(geom, velocities))
     $              call errquit('geom_getsym0: geom_vel_set failed', 0,
     &       GEOM_ERR)
               if (ncenter_before_auto .ne. ncenter) call errquit
     $             ('geom_getsym0: autosym bug : number of atoms wrong',
     $              ncenter, GEOM_ERR)
            else
               write(LuOut,*) ' autosym detected unknown group ', 
     &                        groupname
               call errquit('geom_getsym0: autosym: invalid group',0,
     &       INPUT_ERR)
            endif
         endif

c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c
         nata=ncenter
         call sym_nwc(geom,rtdb,nata,oprint_sym,scale,threquiv,nops)
c     
c     Check that if we used autosym that we ended up with the 
c     same no. of atoms ... if we don't then autosym and nwchemsym
c     don't agree on the orientation of point group elements
c     
         if (oautosym) then
            if (.not. geom_ncent(geom,ncenter)) call errquit
     $           ('geom_getsym0: geom_cent?',0, GEOM_ERR)
            if (ncenter_before_auto .gt. ncenter) call errquit
     $           ('geom_getsym0: autosym bug : too few atoms',ncenter,
     &       GEOM_ERR)
#if 0            
            if (ncenter_before_auto .ne. ncenter) call errquit
     $           ('geom_getsym0: autosym bug : too many atoms',ncenter,
     &       GEOM_ERR)
#else
            if (ncenter_before_auto .ne. ncenter) then
               threquiv=threquiv*0.1d0
               itries=itries+1
            if (itries .gt. 4) call errquit
     $           ('geom_getsym0: too many tries',itries,
     &       GEOM_ERR)
            ncenter=ncenter_before_auto
            if (oprint)
     W        write(luout,*) ' decreasing threquiv to ',threquiv
               goto 2020
            endif
#endif
         endif
c     
         if (oprint_sym.and.(isys.eq.0)) call sym_print_all(geom,
     $        .true.,.true.,.true.,.true.,.true.)
c
         return
         end
C>
C> \brief Establish symmetry related transformations and unique atom
C> lists
C>
C> Establishes the symmetry related transformation and the unique atom
C> list. The molecule might be re-oriented along the moments of 
C> inertia. The data in the geometry is updated as a result.
C>
      subroutine geom_getsym(rtdb,geom,name)
      implicit none
#include "nwc_const.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "inp.fh"
      integer rtdb !< [Input] The RTDB handle
      integer geom !< [Input] The geometry handle
c
      double precision threquiv,scale
      character*(*) name !< [Input] The name of the geometry
      logical geom_get_user_scale,oautosym,s
      external geom_get_user_scale
      character*256 tmp
      integer k
c
      threquiv=1d-2 ! hardwired for now
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit('geom_getsym: failed getting user scale',0,
     &       GEOM_ERR)
      oautosym=.true.
c
      call geom_getsym0(rtdb,geom,oautosym,
     ,     ncenter(geom),isystype(geom),
     ,     name,tags(1,geom),
     ,     threquiv,scale,coords(1,1,geom),
     ,     charge(1,geom),velocities(1,1,geom),geom_atomct(1,geom))
c
c     store stuff to rtdb
c
      tmp = 'geometry:'//trans(geom)(1:lent(geom))
      k = inp_strlen(tmp)+1
      s = .true.
      tmp(k:) = ' ' 
      tmp(k:) = ':no. unique centers'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, ncenter_unique(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group number'
      s = s .and. rtdb_put(rtdb, tmp, mt_int, 1, group_number(geom))
      tmp(k:) = ' ' 
      tmp(k:) = ':group name'
      s = s .and. rtdb_cput(rtdb, tmp, 1, group_name(geom))
      return
      end
C> @}
      subroutine geom_findrot(nat,oldx,newx,rot)
      implicit none
#include "global.fh"
#include "util.fh"
#include "util_params.fh"
#include "stdio.fh"
      integer nat
      double precision oldx(3,*)
      double precision newx(3,*)
      double precision rot(3,3) ! [out]
c
      double precision c(4,4),v(4,4),q(4),d(4)
      double precision xmid,ymid,zmid,norm,det
      integer i
      double precision xxyx, xxyy, xxyz,
     *     xyyx, xyyy, xyyz, xzyx, xzyy, xzyz
      double precision x1,y1,z1,x2,y2,z2
      double precision xrott,yrott,zrott
      logical oprint, oprint_info
c
      oprint = util_print('mo guess', print_default)
      oprint_info = (ga_nodeid() .eq. 0 .and. oprint)
c
c     ----------------------------------------------------------
c     calculate the 3x3 rotation matrix which transforms the old
c     geometry into the active geometry.  Literature reference:
c     S. K. Kearsley, "On the Orthogonal Transformation Used for
c     Structural Comparisons", Acta Crystallographica Section A,
c     45, 208-210 (1989). 
c     ----------------------------------------------------------
c
c     Find the centroid of the old and active geometries
c     and translate it to the origin.
c 
        xmid = 0.0d0
        ymid = 0.0d0
        zmid = 0.0d0
        norm = 0.0d0
        do i = 1, nat
           xmid = xmid + oldx(1,i)
           ymid = ymid + oldx(2,i)
           zmid = zmid + oldx(3,i)
           norm = norm + 1.0d0
        end do
        xmid = xmid / norm
        ymid = ymid / norm
        zmid = zmid / norm
        do i = 1, nat
           oldx(1,i) = oldx(1,i) - xmid
           oldx(2,i) = oldx(2,i) - ymid
           oldx(3,i) = oldx(3,i) - zmid
        end do
c
c     Repeat for the first structure.
c     This centroid position is saved.
c
        xmid = 0.0d0
        ymid = 0.0d0
        zmid = 0.0d0
        norm = 0.0d0
        do i = 1, nat
           xmid = xmid + newx(1,i)
           ymid = ymid + newx(2,i)
           zmid = zmid + newx(3,i)
           norm = norm + 1.0d0
        end do
        xmid = xmid / norm
        ymid = ymid / norm
        zmid = zmid / norm
        do i = 1, nat
           newx(1,i) = newx(1,i) - xmid
           newx(2,i) = newx(2,i) - ymid
           newx(3,i) = newx(3,i) - zmid
        end do
c
c     Use the quaternion - based method to achieve the best fit
c     superposition of the two sets of coordinates.
c
        xxyx = 0.0d0
        xxyy = 0.0d0
        xxyz = 0.0d0
        xyyx = 0.0d0
        xyyy = 0.0d0
        xyyz = 0.0d0
        xzyx = 0.0d0
        xzyy = 0.0d0
        xzyz = 0.0d0
        do i = 1, nat
           x1 = newx(1,i)
           y1 = newx(2,i)
           z1 = newx(3,i)
           x2 = oldx(1,i)
           y2 = oldx(2,i)
           z2 = oldx(3,i)
           xxyx = xxyx + x1*x2
           xxyy = xxyy + y1*x2
           xxyz = xxyz + z1*x2
           xyyx = xyyx + x1*y2
           xyyy = xyyy + y1*y2
           xyyz = xyyz + z1*y2
           xzyx = xzyx + x1*z2
           xzyy = xzyy + y1*z2
           xzyz = xzyz + z1*z2
        end do
        c(1,1) = xxyx + xyyy + xzyz
        c(1,2) = xzyy - xyyz
        c(2,2) = xxyx - xyyy - xzyz
        c(1,3) = xxyz - xzyx
        c(2,3) = xxyy + xyyx
        c(3,3) = xyyy - xzyz - xxyx
        c(1,4) = xyyx - xxyy
        c(2,4) = xzyx + xxyz
        c(3,4) = xyyz + xzyy
        c(4,4) = xzyz - xxyx - xyyy
        c(2,1) = c(1,2)
        c(3,1) = c(1,3)
        c(3,2) = c(2,3)
        c(4,1) = c(1,4)
        c(4,2) = c(2,4)
        c(4,3) = c(3,4)
c
c     diagonalize the quadratic form matrix
c
      call dcopy(16,c,1,v,1)
      call util_jacobi (4,v,4,d)
c
        q(1) = v(1,4)
        q(2) = v(2,4)
        q(3) = v(3,4)
        q(4) = v(4,4)
c
c     Assemble the rotation matrix that superimposes the molecules.
c
        rot(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
        rot(1,2) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
        rot(1,3) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
        rot(2,1) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
        rot(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
        rot(2,3) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
        rot(3,1) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
        rot(3,2) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
        rot(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)
c
      if (oprint_info)   then
        call util_print_centered(luout,'Rotation matrix',15,.true.)
          do i = 1,3
            write(luout,4) rot(1,i), rot(2,i), rot(3,i)
          enddo
          write(luout,*)
          call util_flush(luout)
      endif
   4  format(3F12.6)
c
c print out the rotated coordinates for debugging purposes
c
      do i = 1, nat
         x2 = oldx(1,i)
         y2 = oldx(2,i)
         z2 = oldx(3,i)
         xrott = x2*rot(1,1) + y2*rot(2,1) + z2*rot(3,1)
         yrott = x2*rot(1,2) + y2*rot(2,2) + z2*rot(3,2)
         zrott = x2*rot(1,3) + y2*rot(2,3) + z2*rot(3,3)
         x2 = xrott + xmid
         y2 = yrott + ymid
         z2 = zrott + zmid
         x2 = x2 * cau2ang
         y2 = y2 * cau2ang
         z2 = z2 * cau2ang
         if (ga_nodeid().eq.0) then
         write(6,4) x2,y2,z2
         endif
      end do
c
c      call ga_sync()
c
c calculate the determinant of the rotation matrix
c
      det = rot(1,1)*(rot(2,2)*rot(3,3)-rot(3,2)*rot(2,3))
     &     -rot(1,2)*(rot(2,1)*rot(3,3)-rot(2,3)*rot(3,1))
     &     +rot(1,3)*(rot(2,1)*rot(3,2)-rot(2,2)*rot(3,1))
c
      if (ga_nodeid().eq.0) write(6,*)'det(rot) ',det
      return
      end
      subroutine geom_m3inv(r,rinv)
      implicit none
      double precision, intent(in) :: r(3,3)
      double precision, intent(out) :: rinv(3,3)
c
      double precision detinv
c
      detinv = 1/(r(1,1)*r(2,2)*r(3,3) - r(1,1)*r(2,3)*r(3,2) -
     *     r(1,2)*r(2,1)*r(3,3) + r(1,2)*r(2,3)*r(3,1) +
     *     r(1,3)*r(2,1)*r(3,2) - r(1,3)*r(2,2)*r(3,1))

      rinv(1,1) =  detinv * (r(2,2)*r(3,3) - r(2,3)*r(3,2))
      rinv(2,1) = -detinv * (r(2,1)*r(3,3) - r(2,3)*r(3,1))
      rinv(3,1) =  detinv * (r(2,1)*r(3,2) - r(2,2)*r(3,1))
      rinv(1,2) = -detinv * (r(1,2)*r(3,3) - r(1,3)*r(3,2))
      rinv(2,2) =  detinv * (r(1,1)*r(3,3) - r(1,3)*r(3,1))
      rinv(3,2) = -detinv * (r(1,1)*r(3,2) - r(1,2)*r(3,1))
      rinv(1,3) =  detinv * (r(1,2)*r(2,3) - r(1,3)*r(2,2))
      rinv(2,3) = -detinv * (r(1,1)*r(2,3) - r(1,3)*r(2,1))
      rinv(3,3) =  detinv * (r(1,1)*r(2,2) - r(1,2)*r(2,1))

      return
      end
