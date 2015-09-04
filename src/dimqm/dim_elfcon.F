      subroutine dim_elfcon(basis,geom,g_dens,points,npt,elfval,nder)
c
c $Id: hnd_elfcon.F 19707 2010-10-29 17:59:36Z d3y133 $
c
c     This routine calculates the electronic contribution of the
c     electronic integral defined by nder for a given density at 
c     the grid points defined in points. 
c
c     It returns an array (max(nder*3,1),npts) which holds all 
c     max(nder*3,1) components for each grid point
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
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      double precision points(3,npt) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
c
      integer ishell, jshell, ijshell, nshell, nbf_max, me, nproc
      integer ilo, ihi, jlo, jhi, idim, jdim, nint
      integer l_dens, k_dens, l_scr, k_scr, l_buf, k_buf
      integer maxbuf, maxscr, i
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     ----- calculate buffer and scratch space -----
c           buffer = (lmax*(lmax+1)/2)^2 * (max(nder*3,1) * ngridpoints
c           scratch = see hnd_elfder wrapper routine
c
      call int_init_1eelec(maxbuf,maxscr,basis,nder,npt)
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
c     Zero elfval result array
c
      call dcopy(max(nder*3,1)*npt,0.0d0,0,elfval,1)
c
c     ----- calculate electronic integral component(s) at all points -----
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
               call dim_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt)
c
               call dim_multi_reduce(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           elfval,idim,jdim,npt*(max(nder*3,1)))
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c
c     Collect components from all the nodes for all points
c
      call ga_sync()
      call ga_dgop(msg_efgs_col,elfval,npt*(max(nder*3,1)),'+')
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
      subroutine dim_multi_reduce(mblock,block,rblock,idim,jdim,nblock)
c
      implicit none
#include "stdio.fh"
      integer idim,jdim,nblock
      double precision mblock(idim,jdim,nblock), block(idim,jdim)
      double precision rblock(nblock)
c
      integer iblock,i,j
c
c      write(luout,*) "integral patch:"
c      write(luout,*) mblock
      do iblock = 1, nblock
         do i = 1, idim
            do j = 1, jdim
               rblock(iblock)=rblock(iblock)+mblock(i,j,iblock)*
     &                                       block(i,j)
            enddo
         enddo
      enddo
c
      return
      end
