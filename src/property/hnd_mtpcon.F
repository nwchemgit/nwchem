      subroutine hnd_mtpcon(basis,geom,g_dens,mtpval,mpole,centerl)
c
c $Id$
c
c     This routine calculates the electronic contribution of a given
c     multipole for a given density at the grid points defined in points
c     and adds the nuclear and point charge contributions
c
c     It returns an array ((mpole+1)*(mpole+2)/2) which holds all 
c     (mpole+1)*(mpole+2)/2 components at the mpole coordinate
c     Ordering is defined in NWints/int (defNxyz, getNxyz, matchNxyz)
c
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
c
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer mpole    ! [input] multipole order 
      double precision centerl(3) ! [input] coordinates for coord mpole
      double precision mtpval(*)  ! [output] mpole values
c
      integer ishell, jshell, ijshell, nshell, nbf_max, me, nproc
      integer ilo, ihi, jlo, jhi, idim, jdim, nint
      integer l_dens, k_dens, l_scr, k_scr, l_buf, k_buf
      integer maxbuf, maxscr, ncomp, mtpint
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     ----- calculate buffer and scratch space -----
c           buffer = (lmax*(lmax+1)/2)^2 * (max(nder*3,1) 
c           scratch = see hnd_elfder wrapper routine
c
      call int_init_1emtpl(maxbuf,maxscr,basis,mpole)
c
      if (.not. bas_geom(basis, geom)) call errquit
     $   ('hnd_elfcon: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $   ('hnd_elfcon: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &   ('hnd_elfcon: bas_nbf_cn_max failed',555, BASIS_ERR)
c
      if (.not. ma_push_get(mt_dbl,nbf_max*nbf_max,'dens patch',l_dens,
     &    k_dens)) call errquit('hnd_elfcon: ma 1 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'scratch',l_scr,k_scr))
     &    call errquit('hnd_elfcon: ma 2 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxbuf,'int buf',l_buf,k_buf))
     &    call errquit('hnd_elfcon: ma 3 failed',911,MA_ERR)
c
c     Zero mtpval result array
c
      ncomp = (mpole+1)*(mpole+2)/2
      call dcopy(ncomp,0.0d0,0,mtpval,1)
c
c     ----- calculate multipole integral component(s) at all points -----
c
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi)) call
     &            errquit('hnd_elfcon: bas_cn2bfr',basis,BASIS_ERR)
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c     Get the density patch, make the integrals and contract
c
               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $                     dbl_mb(k_dens), idim)
c

               call dcopy(maxbuf,0.0d0,0,dbl_mb(k_buf),1)
               call int_mpolel(basis,ishell,basis,jshell,mpole,centerl,
     &                         maxscr,dbl_mb(k_scr),maxbuf,
     &                         dbl_mb(k_buf),mtpint)
               call multi_reduce2(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           mtpval,idim,jdim,ncomp)
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c     Collect components from all the nodes for all points
c
      call ga_sync()
      call ga_dgop(msg_efgs_col,mtpval,ncomp,'+')
c
c     Clean up MA data blocks
c
      if (.not.ma_pop_stack(l_buf)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_buf failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_scr)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_scr failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_dens)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_dens failed',911,MA_ERR)
      return
      end
c
      subroutine multi_reduce2(mblock,block,rblock,idim,jdim,nblock)
c
      implicit none
      integer idim,jdim, nblock
      double precision mblock(nblock,jdim,idim), block(idim,jdim)
      double precision rblock(nblock)
c
      integer iblock,i,j
c
      do iblock = 1, nblock
         do i = 1, idim
            do j = 1, jdim
               rblock(iblock)=rblock(iblock)+mblock(iblock,j,i)*
     &                                       block(i,j)
            enddo
         enddo
      enddo
c
      return
      end
