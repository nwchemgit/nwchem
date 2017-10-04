#define DEBUG 1
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
     D     .false.,.true.)
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
     D     .true.,.true.)
      return
      end
cc AJL/Begin/FDE 
      subroutine hnd_elfcon_fde(basis,geom,g_dens,points,npt,elfval,
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
      double precision elfval(*)     ! [output] efg values for each
coord
      call hnd_elfcon_0(basis,geom,g_dens,points,npt,elfval,nder,
     D     .false.,.false.)
      return
      end
cc AJL/End
      subroutine hnd_elfcon_0(basis,geom,g_dens,points,npt,elfval,nder,
     D     doschwarz,do_parallel)
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
#include "util.fh"
c
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      double precision points(3,npt) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
      logical doschwarz ! [input] use and init schwarz screening
cc AJL/Begin/FDE
      logical do_parallel ! [input] for calculating FDE potential we turn
                          ! parallel off as we are already working with 
                          ! nxtask for grid points
cc AJL/End
c
      integer ishell, jshell, ijshell, nshell, nbf_max, me, nproc
      integer ilo, ihi, jlo, jhi, idim, jdim, nint
      integer l_dens, k_dens, l_scr, k_scr, l_buf, k_buf
      integer maxbuf, maxscr, i
      integer nxtask, task_size,next
      integer npasses,i_pass,avail_ma,iptr,npt_pass,mxbynpt,
     M     maxbuf_pass
      external nxtask
      double precision thresh,pmat_max,dabsmax,schw_ij
      external dabsmax
      logical oprint
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
c
check if enough MA mem is avail, if not then mult. passes
c
      avail_ma = ma_inquire_avail(mt_dbl)*7/10
c maxbuf is eq to maxscr, therefore 2d0 factor
      npasses=maxbuf*2d0/avail_ma 
 123  npasses=npasses+1
      oprint=npasses.gt.1.and.ga_nodeid().eq.0
      if(oprint)  then
         write(6,321) ' hnd_elfcon: passes      = ',npasses
      endif
c     check if last loop ok
      npt_pass=npt/npasses
c     value of npt_pass at last loop
      npt_pass=npt-(npasses-1)*npt_pass
      mxbynpt=maxbuf/npt
      maxbuf_pass=mxbynpt*npt_pass
      if(2*maxbuf_pass.gt.avail_ma) goto 123
      if (do_parallel) call ga_igop(msg_efgs_col-1,npasses,1,'max')
      npt_pass=npt/npasses
      mxbynpt=maxbuf/npt
      if(oprint)  then
         write(6,321) ' hnd_elfcon: npt         = ',npt
         write(6,321) ' hnd_elfcon: npt_pass    = ',npt_pass
#ifdef DEBUG
         write(6,321) ' hnd_elfcon: mxbynpt     = ',mxbynpt
         write(6,321) ' hnd_elfcon: maxbuf      = ',maxbuf
         write(6,321) ' hnd_elfcon: maxbuf_pass = ',
     A        mxbynpt*npt_pass
#endif
 321     format(a,i10)
      endif
      if(do_parallel) call ga_sync()
c
c     Zero elfval result array
c
      call dcopy(max(nder*3,1)*npt,0.0d0,0,elfval,1)
c     2ndcall to init_1eelec with new npt_pass
      call int_init_1eelec(maxbuf,maxscr,basis,nder,npt_pass)
c
c     alloc scr and buf for ints
c
      maxscr=max(maxbuf_pass,maxscr)
      if(oprint)  then
         write(6,321) ' hnd_elfcon: maxbuf_pass2= ',maxbuf_pass
      endif
      if (.not. ma_push_get(mt_dbl,maxscr,'scratch',l_scr,k_scr))
     &    call errquit('hnd_elfcon: ma 2 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxbuf_pass,'int buf',
     L        l_buf,k_buf))
     &    call errquit('hnd_elfcon: ma 3 failed',maxbuf_pass,MA_ERR)
      if (doschwarz) call schwarz_init(geom, basis)
      do i_pass=1,npasses
      iptr=(i_pass-1)*npt_pass+1
      npt_pass=npt/npasses
      if(i_pass.eq.npasses) npt_pass=npt-iptr+1

      maxbuf_pass=mxbynpt*npt_pass
      maxscr=maxbuf_pass
      if(oprint)  then
         write(6,421) ' hnd_elfcon: pass = ',i_pass,
     E        ' of ',npasses,100d0*i_pass/npasses,
     E        '% done at ',util_wallsec(),'s.'
 421     format(a,i10,a,i10,f6.1,a,f16.2,a)
#ifdef DEBUG
         write(6,321) '   hnd_elfcon: iptr        = ',iptr
         write(6,321) '   hnd_elfcon: npt_pass    = ',npt_pass
         write(6,321) '   hnd_elfcon: maxbuf_pass = ',maxbuf_pass
#endif
      endif
      call dcopy(maxbuf_pass,0.0d0,0,dbl_mb(k_buf),1)
c
c     ----- calculate electronic integral component(s) at all points -----
c
      ijshell = 0
      task_size = 1
cc AJL/Begin/FDE
      if (do_parallel) then
        next = nxtask(nproc,task_size)+1
      else
        next = 1
      endif
cc AJL/End
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
     &                         nder,points(1,iptr),npt_pass)
#ifdef M12
               if(ishell.ne.jshell)
     D              call dscal(idim*jdim,2d0,dbl_mb(k_dens),1)
#endif
!DEC$ NOINLINE
               call multi_reduce(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           elfval(iptr),idim,jdim,
     N              npt_pass*(max(nder*3,1)))
               endif
cc AJL/Begin/FDE
               if (do_parallel) then
                 next = nxtask(nproc,task_size)+1
               else 
                 next = next + 1
               endif
            end if
cc AJL/End
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
      if (do_parallel) next=nxtask(-nproc,task_size)
      enddo ! i_pass
      if(doschwarz) call schwarz_tidy()
      if (.not.ma_chop_stack(l_scr)) then
       call ma_summarize_allocated_blocks()
         call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_scr failed',1,MA_ERR)
      endif
c     Collect components from all the nodes for all points
c
cc AJL/Begin/FDE
      if (do_parallel) then
        call ga_dgop(msg_efgs_col,elfval,npt*(max(nder*3,1)),'+')
      end if
cc AJL/End
c
c     Clean up MA data blocks
c
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
