      subroutine movecs_rotate(rtdb, geom, basis, nset, g_vecs,
     .     evals, foccs)
C$Id: movecs_rotate.F,v 1.8 2002-12-17 02:39:29 edo Exp $
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
c
c     Rotate the MO's calculated at an old geometry to the active geometry.
c     Syntax: vectors input rotate geomname filename
c
      integer rtdb, geom, basis, g_vecs(2)
      double precision evals(*),foccs(*)
c
      integer iset, nset ! [in] 1 = rhf, 2 = uhf (2 sets_
      integer file_nbf, nbf, i, j, nmo, nsets
      character*80 geomname, basisname, filename, title
      character*80 rotate(2) !max_geom=100 
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
      integer l_occ, k_occ,  g_c, g_v
      integer l_a, k_a, l_b, k_b
      integer oldgeom, nat,iptr
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      logical int_normalize
      external int_normalize
      logical oprint, oprint_info
c
      oprint = util_print('mo guess', print_default)
      oprint_info = (ga_nodeid() .eq. 0 .and. oprint)
c      
      if (.not. rtdb_cget(rtdb, 'vectors:rotate', 2, rotate))
     $     call errquit('movecs_rotate: no file/geometry', 0)
      geomname = rotate(1)
      filename = rotate(2)
c
      if (oprint_info) then
         write(6,*)
         call util_print_centered(6,'Orbital rotation guess',15,
     $        .true.)
         write(6,1)
     $        filename(1:inp_strlen(filename)),
     $        geomname(1:inp_strlen(geomname))
 1       format(/
     $        ' input vectors filename    = ',a/
     $        ' input vectors geometry    = ',a/)
         call util_flush(6)
      endif

c     -------------------------------------------------
c     get cartesian coordinates for the active geometry
c     -------------------------------------------------
       if (.not.geom_cart_get(geom,nat,atom_tags,cact,zan))
     $     call errquit('movecs_rotate: geom_cart_get', 0)
c
c     ----------------------------------------------
c     get cartesian coordinates for the old geometry
c     ----------------------------------------------
c
      if (.not. geom_create(oldgeom, geomname))
     $     call errquit('movecs_rotate: cant create oldgeom', 0)
c
      if (.not.geom_rtdb_load(rtdb, oldgeom, geomname))
     $     call errquit('movecs_rotate: cant load geometry', 0)
c
      if (.not.geom_cart_get(oldgeom,nat,atom_tags,cold,zan))
     $     call errquit('movecs_rotate: geom_cart_get', 0)
c
      if (.not. geom_destroy(oldgeom))
     $     call errquit('movecs_rotate: cant destroy oldgeom', 0)
c
       if (ga_nodeid().eq.0) then
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
c
c     ----------------------------------------------------------
c     calculate the 3x3 rotation matrix which transforms the old
c     geometry into the active geometry.  Literature reference:
c     S. K. Kearsley, "On the Orthogonal Transformation Used for
c     Structural Comparisons", Acta Crystallographica Section A,
c     45, 208-210 (1989).  Adapted from Tinker routine quatfit.f
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
c     Repeat for the first structure.
c     This centroid position is saved.
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
c     Use the quaternion-based method to achieve the best fit
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
      endif   ! ga_nodeid().eq.0
      call ga_sync()
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
      if (ga_nodeid().eq.0) then
        do i = 1,4
        do j = 1,4
          call ga_put(g_c,i,i,j,j,c(i,j),1)
        enddo
        enddo
      endif
      call ga_sync()
c
      call ga_symmetrize(g_c)
      call ga_diag_std(g_c, g_v, d)
c
      if (ga_nodeid().eq.0) then
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
c     Assemble the rotation matrix that superimposes the molecules.
c
        rot(1,1) = q(1)**2.0d0 + q(2)**2.0d0 - q(3)**2.0d0 - q(4)**2.0d0
        rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
        rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
        rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
        rot(2,2) = q(1)**2.0d0 - q(2)**2.0d0 + q(3)**2.0d0 - q(4)**2.0d0
        rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
        rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
        rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
        rot(3,3) = q(1)**2.0d0 - q(2)**2.0d0 - q(3)**2.0d0 + q(4)**2.0d0
c
      endif  ! ga_nodeid().eq.0
      call ga_sync()
c
   4  format(3F12.6)
      if (oprint_info)   then
        call util_print_centered(6,'Rotation matrix',15,.true.)
          do i = 1,3
            write(6,4) rot(i,1), rot(i,2), rot(i,3)
          enddo
      endif
      call ga_sync()
c
c  for error checking:
c  -------------------
c     rotate second molecule to best fit with first molecule
c
c      if (ga_nodeid().eq.0) then
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
c      endif
c      call ga_sync()
c
      if (.not.ga_destroy(g_c)) 
     $    call errquit('movecs_rotate: could not destroy g_c',0)
c
      if (.not.ga_destroy(g_v)) 
     $    call errquit('movecs_rotate: could not destroy g_v',0)
c
c     ---------------------------------------------
c     Read in vectors calculated at the old geometry.
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
      do iset = 1,nset
         iptr=(iset-1)*nbf+1
      if (.not. movecs_read(filename, iset, foccs(iptr),
     $     evals(iptr), g_vecs(iset))) call errquit
     $     ('movecs_project: failed reading vectors', 0) 
      enddo

c     ---------------------------------------------------
c     apply rotation matrix to MO vectors, one at a time
c     ---------------------------------------------------
c
      if (.not. ma_push_get(mt_dbl, nbf, 'vectors in', l_a, k_a))
     $     call errquit('sym_movecs_apply_rot: ma?', nbf)

      if (.not. ma_push_get(mt_dbl, nbf, 'vectors out', l_b, k_b))
     $     call errquit('sym_movecs_apply_rot: ma?', nbf)
c
      call ga_sync()
c
      if (ga_nodeid().eq.0) then
c
        do iset = 1,nset
c
        do j = 1,nmo
         do i = 1,nbf
           dbl_mb(k_a+i-1)=0.0d0
           dbl_mb(k_b+i-1)=0.0d0
         enddo
         call ga_get(g_vecs(iset),1,nbf,j,j,dbl_mb(k_a),1)
         call sym_movecs_apply_rot(basis, rot, dbl_mb(k_a), dbl_mb(k_b))
         call ga_put(g_vecs(iset),1,nbf,j,j,dbl_mb(k_b),1)
c         write(6,*)
c         write(6,*)'mo# ',j
c         do i = 1,nbf
c         write (6,11) dbl_mb(k_a+i-1),  dbl_mb(k_b+i-1)
c         enddo
        enddo
c
        enddo   !iset
c
      endif
      call ga_sync()
c
      if (.not. ma_pop_stack(l_b)) call errquit
     $     ('movecs_rotate: ma pop', 0)
      if (.not. ma_pop_stack(l_a)) call errquit
     $     ('movecs_rotate: ma pop', 0)
c
 11   format(2f12.6)
      return 
      end
c
c
      subroutine sym_movecs_apply_rot(basis, matrix, v, t)
C$Id: movecs_rotate.F,v 1.8 2002-12-17 02:39:29 edo Exp $
      implicit none
#include "geom.fh"
#include "sym.fh"
#include "bas.fh"
#include "mafdecls.fh"
      integer basis                  ! [input]
      double precision matrix(3,3)   ! [input]
      double precision v(*)          ! [input]
      double precision t(*)          ! [output]
c
c     Apply the rotation matrix to the vector of basis function
c     coefficients in v(*), returning the result in t(*).
c     Adapted from ../symmetry/sym_mo_ap_of.F
c
      integer nbf, nbf_per_sh
      integer geom
      integer maxf, maxr, ang_max
      integer l_r, k_r
      integer i
      integer ncent, icent, jcent, ish, jsh, ilo, ihi, jlo, jhi
      integer ishlo, ishhi, jshlo, jshhi
      integer type, ngen, nprim, sphcart
      integer opinv
      logical do_spherical
c
c     Get basis and geom info
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_movecs_apply_rot: no geometry from basis?', basis)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('sym_movecs_apply_rot: bas_numbf?', basis)
      if (.not. bas_high_angular(basis, ang_max))
     $     call errquit('sym_movecs_apply_rot: bas_high_angular', basis)
      if (.not. geom_ncent(geom, ncent)) call errquit
     $     ('sym_bas_irreps: geom_ncent ?', geom)
c     
c     Get the transformation matrices for shells
c
      do_spherical = bas_is_spherical(basis)
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
      if (.not. ma_push_get(MT_DBL, maxr, 'r', l_r, k_r))
     $     call errquit('sym_movecs_apply_rot: ma failed', maxr)
c
      call sym_bas_rot
     &    (geom, matrix, dbl_mb(k_r), maxf, ang_max, do_spherical)
c      call sym_inv_op(geom, op, opinv)
c
**      call dfill(nbf, 0.0d0, t, 1) ! Now defined to add in to result
c
      do icent = 1, ncent
         jcent = icent
c         jcent = sym_center_map(geom, icent, opinv)
         if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $        call errquit('sym_movecs_apply_rot: bas_ce2cnr',0)
         if (.not. bas_ce2cnr(basis, jcent, jshlo, jshhi))
     $        call errquit('sym_movecs_apply_rot: bas_ce2cnr',0)
c
         do ish = ishlo, ishhi
            jsh = jshlo + ish - ishlo ! OP ish -> jsh
c
            if (.not. bas_continfo(basis, ish, type, nprim, ngen,
     $           sphcart)) call errquit
     $           ('sym_movecs_apply_rot: bas_continfo', ish)
            if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $           call errquit('sym_movecs_apply_rot: bas_cn2bfr', ish)
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('sym_movecs_apply_rot: bas_cn2bfr', jsh)
c
            if (type .eq. -1) then
*               call sym_movecs_app_rot2(dbl_mb(k_r), maxf, ang_max,
*     $              0, v(ilo), 1, 1, t(jlo))
               t(jlo) = t(jlo) + v(ilo)
               call sym_movecs_app_rot2(dbl_mb(k_r), maxf, ang_max,
     $              1, v(ilo+1), 3, 1, t(jlo+1))
            else if (type .eq. 0) then
               do i = 0, ngen-1
                  t(jlo+i) = t(jlo+i) + v(ilo+i)  ! S functions are easy
               enddo
            else
               nbf_per_sh = (ihi-ilo+1)/ngen
               call sym_movecs_app_rot2(dbl_mb(k_r), maxf, ang_max,
     $              type, v(ilo), nbf_per_sh, ngen, t(jlo))
            endif
         enddo
      enddo
c
c*      write(6,*) ' sym_movecs_apply_rot: ', op
c      do i = 1, nbf
c         write(6,11) v(i), t(i)
c 11      format(2f12.6)
c      enddo
c
      if (.not. ma_pop_stack(l_r)) call errquit
     $     ('sym_movecs_apply_rot: ma?', 0)
c
      end
      subroutine sym_movecs_app_rot2(r, maxf, ang_max,
     $              type, v, nbf_per_sh, ngen, t)
      implicit none
      integer maxf, ang_max
      double precision r(maxf, maxf, 0:ang_max)
      integer type
      integer nbf_per_sh
      integer ngen
      double precision v(nbf_per_sh, ngen)
      double precision t(nbf_per_sh, ngen)
c
      integer i, j, k
      double precision s
c
c     Apply symmetry operation to MO coefficients of a generally
c     contracted shell
c
      do k = 1, ngen
         do i = 1, nbf_per_sh
            s = v(i,k)
            if (s .ne. 0.0d0) then
               do j = 1, nbf_per_sh
                  t(j,k) = t(j,k) + r(j,i,type)*s
               enddo
            endif
         enddo
      enddo
c
      end
      subroutine sym_bas_rot(geom, matrix, r, maxf, ang_max,
     &    otransform)
C$Id: movecs_rotate.F,v 1.8 2002-12-17 02:39:29 edo Exp $
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c     
      integer geom         ! [input]  geometry handle
      integer op           ! [input]  Desired operator
      double precision matrix(3,3) ! [input]  Desired operator
      integer maxf         ! [input]  Leading dimension of r
      integer ang_max      ! [input]  Max. ang. momentum of shell
      logical otransform   ! [input]  transform to spherical basis
c     
c     Return the transformation matrices for basis functions
c     under the specified rotation matrix.  Adapted from
c     ../symmetry/syb_bas_op.F
c     
c     NOTE that the identity operation is NOT included
c     
c     NOTE that only cartesian shells are supported
c     
c     Let X(I,L) be the I'th function in the shell with angular
c     momentum L.  The application of a symmetry operator will
c     map shell X into an equivalent shell on a possibly different
c     center and will also mix up the components of the shell
c     according to
c     
c     R op X(I,L) = sum(J) X(J,L)*R(J,I,L)
c     
c     A function X(J,L) is just a product of the cartesian
c     coordinates raised to some power.  This product may be
c     expressed as
c     
c     X(I,L) = X(IM,L-1) * X(M,1)
c     
c     X(M,1) is of course just a component of a P function and
c     the indices IM and M are determined so that the product
c     yields X(I,L).
c     
c     Thus, the action of R may be expressed in terms of lower
c     angular momenta shells
c     
c     R op X(I,L) = R op X(IM,L-1) * X(M,1)
c     
c     = sum(JK) sum(K) X(JK,L-1) * X(K,1) * 
c     R(JK,IM,L-1) * R(K,M,1)
c     
c     = sum(J) X(J,L)*R(J,I,L) (by definition)
c     
c     We also have from the transformation of coordinates
c     (note that for consistency coordinates are written as
c     a 3 component column vector and basis functions in a
c     shell as a row vector)
c     
c     (R op r)(i) = sum(j) O(i,j) r(j)
c     
c     or
c     
c     R op r = O * r
c     
c     where O is the transformation matrix for coordinates
c     under operation R.
c     
c     And we also have the following result (obtained by
c     considering transformations of the coordinates instead
c     of the functions)
c     
c     R op f(r) = f(R^-1 op r)
c     
c     Considering how the P (L=1) functions transform
c     (note a p-shell is the row vector X(I,1) = (rT)(I))
c     
c     R op X(r;I,1) = X(R^-1 op r; I, 1)
c     = (O(R^-1)*r)T
c     = rT * O(R^-1)T
c     = rT * O(R)
c     = X(r;J,1) * R(J,I,1)
c     
c     Hence, R(J,I,1) = O(J,I) (but note that they are
c     applied from opposite sides).  Also, note that any
c     translational components of R are irrelevant.
c     
c     Finally, R(1,1,1) = 1 (s functions).
c     
      double precision  ! [output] basis function transformation under operation op
     &    R(1:maxf,1:maxf,0:ang_max) 
      integer II, JJ, KK, LL, MM, IIMM, JJKK, i, j,  nLL, nLLm1
      integer ijk(3), pqr(3), lmn(3)
c     
c     Ensure cartesian basis function polynomials are defined
c     
      call defnxyz(ang_max)
c     
*     call util_print_centered(6,
*     $     'Transformation of basis shells under symmetry operators', 
*     $     30, .true.)
c     
c     Loop through operators in the group
c     
*      write(6,1) op
*    1    format(/' Symmetry operator ', i3)
      call dfill(maxf*maxf*(ang_max+1), 0.0d0, R, 1)
c     
c     Loop thru angular momenta
c     
      do LL = 0, ang_max
         nLL = ((LL+1)*(LL+2))/2
         if (LL .eq. 0) then
c     
c     S functions are easy
c     
            R(1,1,0) = 1.0d0
c     
         else if (LL .eq. 1) then
c     
c     P functions transform using the same matrix that
c     transforms coordinates without the translational
c     components ... except of course that the matrix
c     is applied from the other side.
c     
c            call sym_get_cart_op(geom, op, matrix)
            do i = 1, 3
               do j = 1, 3
                  R(j,i,1) = matrix(i,j)
               enddo
            enddo
c     
         else
c     
c     Higher functions are built recursively
c     
            nLLm1 = ((LL)*(LL+1))/2
c     
c     Loop thru cartesian functions with sum of powers LL
c     and determine the powers of x, y, z
c     
            do II = 1, nLL
               call getNxyz(LL, II, ijk)
c     
c     Express this function as a product of functions in (LL-1)*(1) and
c     determine the indices of these functions (IIMM and MM)
c     
               do MM = 1, 3
                  if (ijk(MM) .gt. 0) then
                     ijk(MM) = ijk(MM)-1
                     call matchNxyz(LL-1, IIMM, ijk)
                     goto 10
                  endif
               enddo
               call errquit('sym_bas_rot: should never get here ', 0)
 10            continue
c     
c     Loop thru LL-1 functions in the result index
c     
               do JJKK = 1, nLLm1
                  call getNxyz(LL-1, JJKK, pqr)
c     
c     Loop thru P functions and add onto the LL-1 functions
c     to get a new LL function
c     
                  do KK = 1, 3
                     lmn(1) = pqr(1)
                     lmn(2) = pqr(2)
                     lmn(3) = pqr(3)
                     lmn(KK) = lmn(KK) + 1
                     call matchNxyz(LL, JJ, lmn)
c     
                     R(JJ,II,LL) = R(JJ,II,LL) + 
     $                    R(JJKK, IIMM, LL-1) * R(KK, MM, 1)
c     
                  enddo
               enddo
            enddo
         endif
c     
*         write(6,2) LL
* 2       format(/'      Angular momentum ', i2)
*         call output(R(1,1,LL), 1, nLL, 1, nLL, maxf, maxf, 1)
      enddo
      if (.not.otransform) return
c
      call sym_bas_rot_spherical(maxf,ang_max,r)
c
      end
      subroutine sym_bas_rot_spherical(maxf,ang_max,r)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer maxf    ! [input] leading dimension(s) of R
      integer ang_max ! [input] trailing dimension of R
      double precision R(1:maxf,1:maxf,0:ang_max)  ! [input/output] cartesian/spherical
c::local
      double precision sum(36,36)  ! intermediate matrix (handles j functions)
      double precision val
*     double precision scale
      integer lval
      integer l2s
      integer s1, s2, c1, c2, r1, r2
      double precision sover(6,6,0:2)
c::functions
      logical spcart_init
      external spcart_init
c::statement function ----- start
      integer iic,iis,iil,lval_max
      double precision Dtrans, Dinvtrans, Dscale
      integer sindex
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
      Dinvtrans(iis,iic,iil) =
     &    dbl_mb((int_mb(k_invsp2c_lindx+iil))+
     &    (iic-1)*(2*iil+1)+iis+iil)
      sindex(iis,lval) = lval+iis+1
c--(2*lm+1)*l + l+i + (lm-l):m=2,a=5
c--(2*lm+1)*l + i + lm      :m=2,a=3
c--lm*(2*l+1) + l + i       :m=2,a=3
      Dscale(iis,lval,lval_max) =
     &    dbl_mb(k_cart_norm_scale + lval_max*(2*lval+1)+lval+iis)
c::statement function ----- end
c
      if (.not.spcart_init(ang_max,.true.,.false.))
     &    call errquit
     &    ('sym_bas_rot_spherical: could not initialize spcart',911)
c
      call dfill ((6*6*3),0.0d00,sover,1)
      sover(1,1,0) = 1.0d00
      sover(1,1,1) = 1.0d00
      sover(2,2,1) = 1.0d00
      sover(3,3,1) = 1.0d00
      sover(1,1,2) = 1.0d00
      sover(2,2,2) = 1.0d00/3.0d00
      sover(3,3,2) = 1.0d00/3.0d00
      sover(4,4,2) = 1.0d00
      sover(5,5,2) = 1.0d00/3.0d00
      sover(6,6,2) = 1.0d00
      sover(1,4,2) = 1.0d00/3.0d00
      sover(1,6,2) = 1.0d00/3.0d00
      sover(4,6,2) = 1.0d00/3.0d00
      sover(4,1,2) = 1.0d00/3.0d00
      sover(6,1,2) = 1.0d00/3.0d00
      sover(6,4,2) = 1.0d00/3.0d00
c
      do lval = 2,ang_max         ! skip s and p transformations
        l2s = (lval+1)*(lval+2)/2
*        write(6,*)' cartesian matrix for angular momentum ',lval
*        call output(R(1,1,lval),1,l2s,1,l2s,maxf,maxf,1)
        call dfill((36*36),0.0d0,sum,1)
        do c2 = 1,l2s
          do s1 = -lval,lval
            val = 0.0d00
            do c1 = 1,l2s
              val = val +
     &            dinvtrans(s1,c1,lval)*R(c1,c2,lval)
            enddo
            sum(sindex(s1,lval),c2) = val
          enddo
        enddo
* now sum is sum(spherical,cartesian,implied lval)
*        write(6,*)' intermediate sum(sph,car)'
*        call output(sum,1,(2*lval+1),1,l2s,36,36,1)
        do s1 = -lval,lval
          do s2 = -lval,lval
            val = 0.0d00
            do c2 = 1,l2s
              val = val +
     &            dtrans(c2,s2,lval)*sum(sindex(s1,lval),c2)
            enddo
            r1 = sindex(s1,lval)
            r2 = sindex(s2,lval)
            R(r1,r2,lval) = val
          enddo
        enddo
*        write(6,*)' spherical matrix for angular momentum ',lval
*        call output(R(1,1,lval),1,(2*lval+1),1,(2*lval+1),maxf,maxf,1)
*        do s1 = 2*lval+2,maxf
*          do s2 = 2*lval+2,maxf
*            R(s1,s2,lval) = 1234.0d00
*          enddo
*        enddo
      enddo
c      
      end
