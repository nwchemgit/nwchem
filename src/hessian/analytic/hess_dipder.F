C
C  $Id$
C
      subroutine hess_dipder(geom,basis,dder,g_dens,n3xyz,ncent)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "msgids.fh"
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
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      lmax = 1
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('hess_dipder: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hess_dipder: bas_numcont failed for basis', basis,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &     ('hess_dipder: bas_nbf_cn_max failed',555, BASIS_ERR)
      if (.not. geom_center_of_mass(geom, center)) call errquit
     *     ('hess_dipder: geom_center_of_mass failed',555, GEOM_ERR)
c     center (1) = 0.0d00
c     center (2) = 0.0d00
c     center (3) = 0.0d00
c
c  Get integral buffers
c
      maxd = nbf_max**2
      call int_mem_dipole(maxbuf,maxscr,basis,basis,lmax)
      maxbuf = 9*maxbuf  ! need xyz derivatives of dipole
      maxscr = max(100000,maxscr)
      if (.not. ma_push_get(mt_dbl,maxbuf,'intd buffer',l_buf,k_buf))
     *  call errquit('hess_dipder:could not get intd buffer',maxbuf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'intd scratch',
     *  l_scr, k_scr))
     *  call errquit('hess_dipder: scratch alloc failed',maxscr,
     &       MA_ERR)
      if(.not. ma_push_get(mt_dbl, maxd, 'mult:dens', l_dens, k_dens))
     &     call errquit('hess_dipder: cannot allocate dens', maxd,
     &       MA_ERR)
c
c  Get nuclear contribution
c
      if (me.eq.0) then
        do i = 1, 3
          do iat = 1, ncent
            index = (iat-1)*3 + i
            dder(i,index) = dder(i,index) + charge(iat,geom)
          enddo
        enddo
      endif
c
c  Start the computation
c
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &        ('hess_dipder: bas_cn2bfr failed for basis', basis,
     &       BASIS_ERR)
         idim = ihi - ilo + 1
 
         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('hess_dipder: bas_cn2bfr', basis,
     &       BASIS_ERR)
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
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c  Get rid of extra memory
c
      if (.not. ma_chop_stack(l_buf)) call errquit
     *  ('hess_dipder: problem chopping stack',555,
     &       MA_ERR)
c
c  Make sure everyone has the full contribution
c
      call ga_dgop(msg_hess_exy,dder,n3xyz*3,'+')
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
      double precision dder(3, 3, *) ! [output]
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
            atind = idatom(at)
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
            dder(1,i,atind) = dder(1,i,atind) - xdip
            dder(2,i,atind) = dder(2,i,atind) - ydip
            dder(3,i,atind) = dder(3,i,atind) - zdip
          enddo
        endif
      enddo
c
      return
      end
c
      subroutine hess_sxdip(g_sx, dder, g_dip, nocc, ndens, 
     *        n3xyz, nbf)
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
      integer nbf            ! [input] no. of basis functions
      double precision dder(3, 3, *) ! [input/output] dipole der. tensor
c
      integer ndns, offset1, offset2, ind, i, j, k, nxyz
      integer alo(3), ahi(3), blo(3), bhi(3)
      double precision factor
c
      if (ndens.eq.1) then
        factor = 2.0d00
      else
        factor = 1.0d00
      endif
      nxyz = n3xyz/3
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
        ind = 1
        do k = 1, nxyz
        do j = 1, 3
          alo(1) = offset1 + ind
          ahi(1) = offset1 + ind
          do i = 1, 3
            blo(1) = offset2 + i
            bhi(1) = offset2 + i
            dder(i, j, k) = dder(i, j, k) + factor *
     *                     nga_ddot_patch(g_sx, 'N', alo, ahi,
     *                                    g_dip, 'N', blo, bhi)
          enddo
          ind = ind + 1
        enddo
        enddo
      enddo
c
      return
      end
c
      subroutine hess_resdip(g_res, dder, g_dip, nocc, nvirt, ndens,
     *        n3xyz, nmo)
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
      integer nmo            ! [input] no. of MOs (alpha or beta)
      double precision dder(3, 3, *) ! [input/output] dipole der. tensor
c
      integer ndns, offset, ind, i, j, k, nxyz
      integer alo(2), ahi(2), blo(3), bhi(3)
      double precision factor
c
      if (ndens.eq.1) then
        factor = 4.0d00
      else
        factor = 2.0d00
      endif
      nxyz = n3xyz/3
      ahi(1) = 0    ! initialize
      do ndns = 1, ndens
        alo(1) = 1
        if (ndns.ne.1) alo(1) = alo(1) + nocc(ndns-1)*nvirt(ndns-1)
        ahi(1) = ahi(1) + nocc(ndns)*nvirt(ndns)
        blo(2) = 1
        bhi(2) = nocc(ndns)
        blo(3) = nocc(ndns) + 1
        bhi(3) = nmo
        offset = (ndns - 1) * 3
        ind = 1
        do k = 1, nxyz
        do j = 1, 3
          alo(2) = ind
          ahi(2) = ind
          do i = 1, 3
            blo(1) = offset + i
            bhi(1) = offset + i
            dder(i, j, k) = dder(i, j, k) + factor *
     *                     nga_ddot_patch(g_res, 'N', alo, ahi,
     *                                    g_dip, 'T', blo, bhi)
          enddo
          ind = ind + 1
        enddo
        enddo
      enddo
c
      return
      end
c
      subroutine hess_dipint(rtdb,geom,basis,g_dip,ndens,nbf)
      implicit none
#include "errquit.fh"
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
      integer ndens           ! [input] number of density matrices
      integer nbf             ! [input] number of basis functions
c
      integer nshell, nbf_max
      integer l_scr, k_scr, l_buf, k_buf
      integer lmax, maxbuf, maxscr, maxd
      integer me, nproc, ishell, jshell, ilo, ihi, jlo, jhi
      integer ijshell, idim, jdim, nint, numint, offset, i
      integer lo(3), hi(3), ld(2)
      integer blo(3), bhi(3)
      double precision center(3)
c
      lmax = 1
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('hess_dipder: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('hess_dipder: bas_numcont failed for basis', basis,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &     ('hess_dipder: bas_nbf_cn_max failed',555, BASIS_ERR)
      if (.not. geom_center_of_mass(geom, center)) call errquit
     *     ('hess_dipder: geom_center_of_mass failed',555, BASIS_ERR)
c     center (1) = 0.0d00
c     center (2) = 0.0d00
c     center (3) = 0.0d00
c
c  Get integral buffers
c
      call int_init(rtdb,1,basis)
      maxd = nbf_max**2
      call int_mem_dipole(maxbuf,maxscr,basis,basis,lmax)
      maxbuf = 3*maxbuf  ! need xyz derivatives
      maxscr = max(100000,maxscr)
      maxscr = max(maxscr, maxbuf) ! this allows it to be used for temp
      if (.not. ma_push_get(mt_dbl,maxbuf,'intd buffer',l_buf,k_buf))
     *  call errquit('hess_dipder:could not get intd buffer',maxbuf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'intd scratch',
     *  l_scr, k_scr))
     *  call errquit('hess_dipder: scratch alloc failed',maxscr,
     &       MA_ERR)
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
     &        ('hess_dipder: bas_cn2bfr failed for basis', basis,
     &       BASIS_ERR)
         idim = ihi - ilo + 1
 
         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('hess_dipder: bas_cn2bfr', basis,
     &       BASIS_ERR)
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c  Get the integrals
c
               call int_mpolel(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              maxbuf, dbl_mb(k_buf),numint)
c
               call hess_dip_in_g(g_dip, dbl_mb(k_buf),
     &              dbl_mb(k_scr), ilo, ihi, jlo, jhi, idim, jdim)
c
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c  Get rid of extra memory
c
      if (.not. ma_chop_stack(l_buf)) call errquit
     *  ('hess_dipder: problem chopping stack',555, MA_ERR)
c
c  Terminate integral api
c
      call int_terminate()
c
c  Now copy the integrals to the beta block if UHF
c
      if (ndens.eq.2) then
        lo(1) = 1
        hi(1) = 3
        lo(2) = 1
        hi(2) = nbf
        lo(3) = 1
        hi(3) = nbf
        blo(1) = 4
        bhi(1) = 6
        blo(2) = 1
        bhi(2) = nbf
        blo(3) = 1
        bhi(3) = nbf
        call nga_copy_patch('N', g_dip, lo, hi,
     *      g_dip, blo, bhi)
      endif
c
      return
      end
c
      subroutine hess_dip_in_g(g_dip, buf, scr,
     *   ilo, ihi, jlo, jhi, idim, jdim)
      implicit none
c
      integer g_dip          ! [output] dipole matrices
      integer ilo            ! [input] low i index
      integer ihi            ! [input] high i index
      integer jlo            ! [input] low j index
      integer jhi            ! [input] high j index
      integer idim           ! [input] i dimension
      integer jdim           ! [input] j dimension
      double precision buf(3,jdim,idim)   ! [input] integrals
      double precision scr(idim,jdim)     ! [input] scratch space
c
      integer i, j, xyz, nint
      integer lo(3), hi(3), ld(2)
c
      nint = idim * jdim
c
      ld(1) = 1
      ld(2) = idim
      lo(2) = ilo
      hi(2) = ihi
      lo(3) = jlo
      hi(3) = jhi
      do xyz = 1, 3
        lo(1) = xyz
        hi(1) = xyz
        call dfill(nint,0.0d0,scr,1)
        do j = 1, jdim ! gather proper integrals
          do i = 1, idim
            scr(i,j) = scr(i,j) + buf(xyz,j,i)
          enddo
        enddo
        call nga_put(g_dip,lo,hi,scr,ld)
      enddo
      return
      end
