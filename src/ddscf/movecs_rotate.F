      subroutine movecs_rotate(rtdb, geom, basis, iset, g_vecs)
C$Id: movecs_rotate.F,v 1.1 2002-12-03 17:16:48 d3k958 Exp $
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
c
c     Rotate the MO's calculated for a particular geometry to the current geometry
c     using the symmetry routine sym_bas_op.
c
      integer rtdb, geom, basis, g_vecs(2)
      integer iset ! [in] 1 = alpha 2 =beta
      integer file_nbf, nbf, i, j, nmo, nsets
      character*80 geomname, basisname, filename, title
      character*80 rotate(2), name(100) !max_geom=100 
      character*20 scftype
      character*16 atom_tags
      double precision cold(3,1000), cact(3,1000) ! max # atoms
      double precision x1(1000), y1(1000), z1(1000)
      double precision x2(1000), y2(1000), z2(1000)
      double precision c(4,4), v(4,4), rot(3,3), q(4), d(4)
      double precision xxyx,xxyy,xxyz
      double precision xyyx,xyyy,xyyz
      double precision xzyx,xzyy,xzyz
      double precision xmid, ymid, zmid, norm
      double precision xrot, yrot, zrot
      double precision zan(1000) ! max # atoms
      integer l_occ, k_occ, l_evals, k_evals, g_c, g_v
      integer ngeom, len, leng, oldgeom, nat
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical int_normalize
      external int_normalize
      logical oprint, oprint_info
c
      oprint = util_print('mo guess', print_default)
      oprint_info = (ga_nodeid() .eq. 0 .and. oprint)
      
      if (oprint_info) then
         write(6,*)
         call util_print_centered(6,'Orbital rotation guess',15,
     $        .true.)
         call util_flush(6)
      endif

      if (.not. rtdb_cget(rtdb, 'vectors:rotate', 2, rotate))
     $     call errquit('movecs_rotate: no file/geometry', 0)
      geomname = rotate(1)
      filename = rotate(2)
c
c     -------------------------------------------------
c     get cartesian coordinates for the active geometry
c     -------------------------------------------------
       if (.not.geom_cart_get(geom,nat,atom_tags,cact,zan))
     $     call errquit('movecs_rotate: geom_cart_get', 0)
c
c     ---------------------
c     get the old geometry 
c     ---------------------
c
      if (.not.rtdb_get(rtdb, 'geometry:ngeom', mt_int, 1, ngeom))
     $     call errquit('movecs_rotate: cant get ngeom_rtdb', 0)
c
      if (.not.rtdb_cget(rtdb, 'geometry:names', ngeom, name))
     $     call errquit('movecs_rotate: cant get ngeom_rtdb', 0)
c
      leng = inp_strlen(geomname)
      do i = 1, ngeom
        len = inp_strlen(name(i))
        if (name(i)(1:len) .eq. geomname(1:leng))
     $    oldgeom = i
      end do
c
       if (.not.geom_rtdb_load(rtdb, oldgeom, geomname))
     $     call errquit('movecs_rotate: cant load geometry', 0)
c
       if (.not.geom_cart_get(oldgeom,nat,atom_tags,cold,zan))
     $     call errquit('movecs_rotate: geom_cart_get', 0)
c
       do i = 1,nat
         x1(i) = cold(1,i)
         y1(i) = cold(2,i)
         z1(i) = cold(3,i)
       enddo
c
       do i = 1,nat
         x2(i) = cact(1,i)
         y2(i) = cact(2,i)
         z2(i) = cact(3,i)
       enddo
c
c     ----------------------------------------------------------
c     calculate the 3x3 rotation matrix which transforms the old
c     geometry into the active geometry
c     ----------------------------------------------------------
c
c     find the centroid of the old and active geometries
c     and translate it to the origin

      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nat
         xmid = xmid + x2(i)
         ymid = ymid + y2(i)
         zmid = zmid + z2(i)
         norm = norm + 1.0d0
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, nat
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do
c
c     now repeat for the first structure, note
c     that this centroid position gets returned
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nat
         xmid = xmid + x1(i)
         ymid = ymid + y1(i)
         zmid = zmid + z1(i)
         norm = norm + 1.0d0
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, nat
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do
c
c  now use the quaternion-based method to achieve the best fit
c  superposition of the two sets of coordinates

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
         xxyx = xxyx + x1(i)*x2(i)
         xxyy = xxyy + y1(i)*x2(i)
         xxyz = xxyz + z1(i)*x2(i)
         xyyx = xyyx + x1(i)*y2(i)
         xyyy = xyyy + y1(i)*y2(i)
         xyyz = xyyz + z1(i)*y2(i)
         xzyx = xzyx + x1(i)*z2(i)
         xzyy = xzyy + y1(i)*z2(i)
         xzyz = xzyz + z1(i)*z2(i)
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
c
      if (.not. ga_create(mt_dbl, 4,4, 'c',1,1,g_c))
     $     call errquit('movecs_rotate: could not create g_c',0)
c 
      if (.not. ga_create(mt_dbl, 4,4, 'v',1,1,g_v))
     $     call errquit('movecs_rotate: could not create g_v',0)
c 
      call ga_zero(g_c)
      call ga_zero(g_v)
c
      do i = 1,4
      do j = 1,4
        call ga_put(g_c,i,i,j,j,c(i,j),1)
      enddo
      enddo
c
c      call ga_print (g_c)
c
      call ga_symmetrize(g_c)
      call ga_diag_std(g_c, g_v, d)
c
c      call ga_print (g_v)
c
      do i = 1,4
      do j = 1,4
        call ga_get(g_v,i,i,j,j,v(i,j),1)
      enddo
      enddo
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble rotation matrix that superimposes the molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
c      do i = 1, nat
c         xrot = x2(i)*rot(1,1) + y2(i)*rot(1,2) + z2(i)*rot(1,3)
c         yrot = x2(i)*rot(2,1) + y2(i)*rot(2,2) + z2(i)*rot(2,3)
c         zrot = x2(i)*rot(3,1) + y2(i)*rot(3,2) + z2(i)*rot(3,3)
c         x2(i) = xrot + xmid
c         y2(i) = yrot + ymid
c         z2(i) = zrot + zmid
c      print *,x2(i)/1.889725989d0
c     $       ,y2(i)/1.889725989d0
c     $       ,z2(i)/1.889725989d0
c      end do
c
      if (.not.ga_destroy(g_c)) 
     $    call errquit('movecs_rotate: could not destroy g_c',0)
c
      if (.not.ga_destroy(g_v)) 
     $    call errquit('movecs_rotate: could not destroy g_v',0)
c
c     ---------------------------------------------
c     read in vectors calculated with old geometry
c     ---------------------------------------------
c
      call util_file_name_resolve(filename, .false.)
c
      if (.not. movecs_read_header(filename, title,
     $     basisname, scftype, file_nbf, nsets, nmo, 2))
     $     call errquit
     $     ('movecs_rotate: failed to read old header', 0)
c
      if (.not. (bas_numbf(basis, nbf)))
     $     call errquit('movecs rotate: bas_numbf?', 0)
c
      if (nbf .lt. file_nbf) then
          call errquit
     $     ('movecs_rotate: nbf on file inconsistent with basis', 0)
      endif
c
      if (.not. ma_push_get(mt_dbl, nbf, 'occ',l_occ, k_occ))
     $     call errquit('movecs_rotate: ma occ', nbf)
c
      if (.not. ma_push_get(mt_dbl, nbf, 'evals',l_evals, k_evals))
     $     call errquit('movecs_rotate: ma occ', nbf)
c
      if (.not. movecs_read(filename, iset, dbl_mb(k_occ),
     $     dbl_mb(k_evals), g_vecs(iset))) call errquit
     $     ('movecs_project: failed reading vectors', 0)
c
c      call ga_print (g_vecs(iset))
c




c
      if (.not. ma_pop_stack(l_evals)) call errquit
     $     ('movecs_rotate: ma pop', 0)
      if (.not. ma_pop_stack(l_occ)) call errquit
     $     ('movecs_rotate: ma pop', 0)


      end
