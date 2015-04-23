      subroutine hnd_elfcon(basis,geom,g_dens,points,npt,elfval,nder)
c
c $Id$
c
c     This routine calculates the electronic contribution of the
c     electronic integral defined by nder for a given density at 
c     the grid points defined in points. 
c
c     It returns an array (max(nder*3,1),npts) which holds all 
c     max(nder*3,1) components for each grid point
c
      implicit none
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      double precision points(*) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
      call hnd_elfcon_0(basis,geom,g_dens,points,npt,elfval,nder,
     D     .false.)
      return
      end 
      subroutine hnd_elfcon_schw(basis,geom,g_dens,points,npt,elfval,
     N     nder)
c
c     modified hnd_elfcon: Schwarz screening is used
c     more details below
c
      implicit none
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      double precision points(*) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
      call hnd_elfcon_0(basis,geom,g_dens,points,npt,elfval,nder,
     D     .true.)
      return
      end 
      subroutine hnd_elfcon_0(basis,geom,g_dens,points,npt,elfval,nder,
     D     doschwarz)
c
c     kernel routine for hnd_elfcon
c     can do schwarz screening if needed
c     Careful: since it initialize and cleans up schwarz, could conflict
c     previous schwarz initializations
#define M12 1
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
#include "schwarz.fh"
c
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      double precision points(3,npt) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
      logical doschwarz ! [input] use and init schwarz screening
c
      integer ishell, jshell, ijshell, nshell, nbf_max, me, nproc
      integer ilo, ihi, jlo, jhi, idim, jdim, nint
      integer l_dens, k_dens, l_scr, k_scr, l_buf, k_buf
      integer maxbuf, maxscr, i
      integer nxtask, task_size,next
      external nxtask
      double precision thresh,pmat_max,dabsmax,schw_ij
      external dabsmax
      parameter (thresh=1d-12)
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
      if (doschwarz) call schwarz_init(geom, basis)
      ijshell = 0
      task_size = 1
      next = nxtask(nproc,task_size)+1
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
         idim = ihi - ilo + 1

#ifndef M12
         do jshell = 1, nshell
#else
         do jshell=1,ishell
#endif
            if(doschwarz) then
               schw_ij=schwarz_shell(ishell,jshell)
            else
               schw_ij=1d99
            endif
            if(schw_ij.gt.thresh*1d-1) then
            ijshell = ijshell + 1
            if (ijshell.eq.next) then
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

               if(doschwarz) then
                  pmat_max = dabsmax(idim*jdim,dbl_mb(k_dens))
               else
                  pmat_max=1d0
               endif
c
           if(schw_ij*pmat_max.gt.thresh) then
c
               call int_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt)
#ifdef M12
               if(ishell.ne.jshell)
     D              call dscal(idim*jdim,2d0,dbl_mb(k_dens),1)
#endif
!DEC$ NOINLINE
               call multi_reduce(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           elfval,idim,jdim,npt*(max(nder*3,1)))
               endif
               next = nxtask(nproc,task_size)+1
               endif
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
      if(doschwarz) call schwarz_tidy()
c     Collect components from all the nodes for all points
c
      next = nxtask(-nproc,task_size)

      call ga_mask_sync(.false.,.true.)
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
      subroutine multi_reduce(mblock,block,rblock,idim,jdim,nblock)
c
      implicit none
      integer idim,jdim,nblock
      double precision mblock(idim,jdim,nblock), block(idim,jdim)
      double precision rblock(nblock) ! [out]
c
      integer iblock,i,j
         do j = 1, jdim
      do iblock = 1, nblock
!DEC$ LOOP COUNT MAX=40
            do i = 1, idim
               rblock(iblock)=rblock(iblock)+mblock(i,j,iblock)*
     &                                       block(i,j)
            enddo
         enddo
      enddo
c
      return
      end
