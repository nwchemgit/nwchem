C
C  $Id: hess_dipder.F,v 1.1 2002-01-19 01:05:02 windus Exp $
C
      subroutine hess_dipder(geom,basis,dder,g_dens,n3xyz,ncent)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
c  This routine takes the dipole derivative integrals and condenses
c  them with the total density to get an initial contribution into
c  the dipole derivative tensor.  The density passed in is assumed
c  to be the total density.
c
      integer geom            ! [input] geometry handle
      integer basis           ! [input] basis handle
      integer ncent           ! [input] number of atomic centers
      integer n3xyz           ! [input] 3 * ncent
      double precision dder(3, n3xyz) ! [input/output] dipole derivative tensor
      integer g_dens          ! [input] total density
c
      integer i, iat, index, nshell, nbf_max, nint, nintout
      integer l_scr, k_scr, l_buf, k_buf, l_dens, k_dens
      integer lmax, maxbuf, maxscr, maxd
      integer me, nproc, ishell, jshell, ilo, ihi, jlo, jhi
      integer ijshell, idim, jdim, idatom(3)
      double precision center(3)
c
      logical geom_center_of_mass 
      external geom_center_of_mass
c
      lmax = 1
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('hess_dipder: bad basis', 555)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hess_dipder: bas_numcont failed for basis', basis)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &     ('hess_dipder: bas_nbf_cn_max failed',555)
      if (.not. geom_center_of_mass(geom, center)) call errquit
     *     ('hess_dipder: geom_center_of_mass failed',555)
c
c  Get nuclear contribution
c
      do i = 1, 3
        do iat = 1, ncent
          index = (iat-1)*3 + i
          dder(i,index) = dder(i,index) + charge(iat,geom)
        enddo
      enddo
c
c  Get integral buffers
c
      maxd = nbf_max**2
      call int_mem_dipole(maxbuf,maxscr,basis,basis,lmax)
      maxbuf = 9*maxbuf  ! need xyz derivatives of dipole
      maxscr = max(100000,maxscr)
      if (.not. ma_push_get(mt_dbl,maxbuf,'intd buffer',l_buf,k_buf))
     *  call errquit('hess_dipder:could not get intd buffer',maxbuf)
      if (.not. ma_push_get(mt_dbl,maxscr,'intd scratch',
     *  l_scr, k_scr))
     *  call errquit('hess_dipder: scratch alloc failed',maxscr)
      if(.not. ma_push_get(mt_dbl, maxd, 'mult:dens', l_dens, k_dens))
     &     call errquit('hess_dipder: cannot allocate dens', maxd)
c
c  Start the computation
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &        ('hess_dipder: bas_cn2bfr failed for basis', basis)
         idim = ihi - ilo + 1
 
         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('hess_dipder: bas_cn2bfr', basis)
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c     Get the density patch, make the integrals and contract
c
               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $              dbl_mb(k_dens), idim)
c
               call intd_mpolel(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              maxbuf, dbl_mb(k_buf),nintout,idatom)
c
               call contract_ddip(dbl_mb(k_dens), idim, jdim,
     &              idatom, dbl_mb(k_buf), maxbuf, dder, n3xyz)
c
c              ind = 0
c              do i = 1, idim
c                 ioff = (i-1)*jdim + k_dens - 1
c                 do m = 1, ((l+1)*(l+2))/2
c                    moff = (l*(l+1)*(l+2))/6 + k_m - 1 + m
c                    do j = 1, jdim
c                       dbl_mb(moff) = dbl_mb(moff) -
c    $                       dbl_mb(j+ioff)*dbl_mb(ind)
c                       ind = ind + 1
c                    end do
c     enddo
c              end do
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c  Get rid of extra memory
c
      if (.not. ma_chop_stack(l_buf)) call errquit
     *  ('hess_dipder: problem chopping stack',555)
c
      return
      end
c
      subroutine contract_ddip(dens, idim, jdim, idatom, 
     *        buf, lbuf, dder, n3xyz)
      implicit none
c
      integer idim            ! [input] i dim of dens
      integer jdim            ! [input] j dim of dens
      integer idatom(3)       ! [input] derivative wrt atom
      integer lbuf            ! [input] length of int. 
      integer n3xyz           ! [input] 3 times no. of atoms
      double precision dens(idim, jdim) ! [input] density
      double precision buf(lbuf)  ! [input] integrals
      double precision dder(3, n3xyz) ! [output]
c
      integer index, at, i, atind, ip1, ip2
      double precision xdip, ydip, zdip, zero
      data zero /0.0d0/
c
      index = 1
      do at = 1, 2                ! over each center
        if (idatom(at).ge.1) then
          do i = 1, 3             ! x, y, z derivative
            xdip = zero
            ydip = zero
            zdip = zero
            atind = (idatom(at) - 1) * 3 + 1
c
            do ip1 = 1, idim      ! over basis functions
              do ip2 = 1, jdim
                xdip = xdip + dens(ip1,ip2) * buf(index)
                index = index + 1
                ydip = ydip + dens(ip1,ip2) * buf(index)
                index = index + 1
                zdip = zdip + dens(ip1,ip2) * buf(index)
                index = index + 1
              enddo
            enddo
c
            dder(i,atind) = dder(i,atind) + xdip
            atind = atind + 1
            dder(i,atind) = dder(i,atind) + ydip
            atind = atind + 1
            dder(i,atind) = dder(i,atind) + zdip
            atind = atind + 1
          enddo
        endif
      enddo
c
      return
      end
c
      subroutine hess_sxdip(g_sx, dder, g_dip, nocc, ndens, 
     *        n3xyz, num)
      implicit none
c
#include "global.fh"
c
c  This routine combines the overlap derivatives and the dipole
c  integrals in the MO basis and puts the result into the dipole
c  derivative tensor.
c
      integer g_sx           ! [input] Sx in the MO basis
      integer g_dip          ! [input] dipole integrals in MO basis
      integer ndens          ! [input] no. of densities
      integer nocc(ndens)    ! [input] no. of occ.
      integer n3xyz          ! [input] 3 times the no. of atoms
      integer num            ! [input] no. of basis functions
      double precision dder(3, n3xyz) ! [input/output] dipole der. tensor
c
      integer ndns, offset1, offset2, ind, i
      integer alo(3), ahi(3), blo(3), bhi(3)
c
      alo(2) = 1
      alo(3) = 1
      blo(2) = 1
      blo(3) = 1
      do ndns = 1, ndens
        ahi(2) = nocc(ndns)
        ahi(3) = nocc(ndns)
        bhi(2) = nocc(ndns)
        bhi(3) = nocc(ndns)
        offset1 = (ndns - 1) * n3xyz
        offset2 = (ndns - 1) * 3
        do ind = 1, n3xyz
          alo(1) = offset1 + ind
          ahi(1) = offset1 + ind
          do i = 1, 3
            blo(1) = offset2 + i
            bhi(1) = offset2 + i
            dder(i, ind) = dder(i, ind) +
     *                     nga_ddot_patch(g_sx, 'N', alo, ahi,
     *                                    g_dip, 'N', blo, bhi)
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine hess_resdip(g_res, dder, g_dip, nocc, nvirt, ndens,
     *        n3xyz, num)
      implicit none
c
#include "global.fh"
c
c  This routine combines the responses and the dipole
c  integrals in the MO basis and puts the result into the dipole
c  derivative tensor.  Note that the responses are 2D packed
c  with the alpha and beta together.
c
      integer g_res          ! [input] responses in the MO basis
      integer g_dip          ! [input] dipole integrals in MO basis
      integer ndens          ! [input] no. of densities
      integer nocc(ndens)    ! [input] no. of occ.
      integer nvirt(ndens)   ! [input] no. of virt.
      integer n3xyz          ! [input] 3 times the no. of atoms
      integer num            ! [input] no. of basis functions
      double precision dder(3, n3xyz) ! [input/output] dipole der. tensor
c
      integer ndns, offset, ind, i
      integer alo(2), ahi(2), blo(3), bhi(3)
c
      ahi(1) = 0    ! initialize
      do ndns = 1, ndens
        alo(1) = 1
        if (ndns.ne.1) alo(1) = alo(1) + nocc(ndns-1)*nvirt(ndns-1)
        ahi(1) = ahi(1) + nocc(ndns)*nvirt(ndns)
        blo(2) = 1
        bhi(2) = nocc(ndns)
        blo(3) = nocc(ndns) + 1
        bhi(3) = num
        offset = (ndns - 1) * 3
        do ind = 1, n3xyz
          alo(2) = ind
          ahi(2) = ind
          do i = 1, 3
            blo(1) = offset + i
            bhi(1) = offset + i
            dder(i, ind) = dder(i, ind) +
     *                     nga_ddot_patch(g_res, 'N', alo, ahi,
     *                                    g_dip, 'T', blo, bhi)
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine hess_dipint(rtdb,geom,basis,g_dip)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
c
c  This routine takes the dipole integrals and stores them in
c  3 matrices.  Note that int_init and int_terminate are called
c  inside this routine.
c
      integer rtdb            ! [input] rtdb handle
      integer geom            ! [input] geometry handle
      integer basis           ! [input] basis handle
      integer g_dip           ! [output] dipole integrals
c
      integer nshell, nbf_max
      integer l_scr, k_scr, l_buf, k_buf, l_dens, k_dens
      integer lmax, maxbuf, maxscr, maxd
      integer me, nproc, ishell, jshell, ilo, ihi, jlo, jhi
      integer ijshell, idim, jdim, nint, offset, i
      integer lo(3), hi(3), ld(2)
      double precision center(3)
c
      lmax = 1
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('hess_dipder: bad basis', 555)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hess_dipder: bas_numcont failed for basis', basis)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &     ('hess_dipder: bas_nbf_cn_max failed',555)
      if (.not. geom_center_of_mass(geom, center)) call errquit
     *     ('hess_dipder: geom_center_of_mass failed',555)
c
c  Get integral buffers
c
      call int_init(rtdb,1,basis)
      maxd = nbf_max**2
      call int_mem_dipole(maxbuf,maxscr,basis,basis,lmax)
      maxbuf = 3*maxbuf  ! need xyz derivatives
      maxscr = max(100000,maxscr)
      if (.not. ma_push_get(mt_dbl,maxbuf,'intd buffer',l_buf,k_buf))
     *  call errquit('hess_dipder:could not get intd buffer',maxbuf)
      if (.not. ma_push_get(mt_dbl,maxscr,'intd scratch',
     *  l_scr, k_scr))
     *  call errquit('hess_dipder: scratch alloc failed',maxscr)
      if(.not. ma_push_get(mt_dbl, maxd, 'mult:dens', l_dens, k_dens))
     &     call errquit('hess_dipder: cannot allocate dens', maxd)
c
c  Start the computation
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &        ('hess_dipder: bas_cn2bfr failed for basis', basis)
         idim = ihi - ilo + 1
 
         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('hess_dipder: bas_cn2bfr', basis)
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c  Get the integrals
c
               call int_mpole(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              maxbuf, dbl_mb(k_buf))
c
               lo(2) = ilo
               hi(2) = ihi
               lo(3) = jlo
               hi(3) = jhi
               ld(1) = 1
               ld(2) = idim
c
               offset = nint
               do i = 1, 3
                  lo(1) = i
                  hi(1) = i
                  call nga_acc(g_dip, lo, hi, dbl_mb(k_buf+offset),
     &                         ld, 1.0d0)
                  offset = offset + nint
               end do
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c  Get rid of extra memory
c
      if (.not. ma_chop_stack(l_buf)) call errquit
     *  ('hess_dipder: problem chopping stack',555)
c
c  Terminate integral api
c
      call int_terminate()
c
      return
      end
