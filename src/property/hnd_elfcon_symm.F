      subroutine hnd_elfcon_symm(basis,geom,g_dens,points,npt,
     &                           elfval,nder)
c
c $Id: hnd_elfcon.F,v 1.4 2004-12-20 22:41:07 bert Exp $
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
#include "zora.fh" ! Use g_AtNr,nlist 
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
      integer l_munu,k_munu,n_munu,n_munu1
      integer l_dens1,k_dens1
      integer l_AtNr,k_AtNr
      ! g_munuV6 contains unique munu-EFG
      ! in sequence: xx,yy,zz,xy,xz,yz
      ! each chunk is nbf*nbf/2 in size
      integer ipolmunu
      integer g_munuV6,g_munu_rhoS,
     &        g_zora_scale_munu(2)
      double precision val
      integer count,j,k,nbf,ndir,indx,indx1,ii

      integer count2,indx2

      double precision ac,ac1           ! Added by FA
      external multi_reduce     ! defined in hnd_eflcon.F
      external multi_reduce_off ! for off-diag elements
      common/munu_EFG/g_munuV6,g_munu_rhoS,
     &                g_zora_scale_munu,ipolmunu ! munu-matrices
c     Note.- g_munu_rhoS is created in dft_zora_rhos.F
c            g_munuV6    is created in hnd_elfcon_symm.F
c            g_munu_rhoS is already available but not used yet.
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
      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('hnd_elfcon: bas_numbf failed',555, BASIS_ERR)
c
      if (.not. ma_push_get(mt_dbl,nbf_max*nbf_max,'dens patch',l_dens,
     &    k_dens)) call errquit('hnd_elfcon: ma 1 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'scratch',l_scr,k_scr))
     &    call errquit('hnd_elfcon: ma 2 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxbuf,'int buf',l_buf,k_buf))
     &    call errquit('hnd_elfcon: ma 3 failed',911,MA_ERR)
c ------- create munu-scratch-array--- START
c  Allocate memory for l_AtNr,k_AtNr 
c  npt=Nr of atoms
      if (.not.ma_alloc_get(mt_dbl,npt,'AtNr',
     &    l_AtNr,k_AtNr))
     &    call errquit('hnd_elfcon_symm: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,npt,dbl_mb(k_AtNr),1)

c     nlist, Nr of selected atoms defined in zora.fh
      ndir=6 ! Nr. directions xx,yy,zz,xy,xz,yz
      n_munu=nbf*nbf*ndir*nlist
c     write(*,*) "nbf=",nbf," n_munu=",n_munu,
c    &           " nlist=",nlist
       if (.not. ma_push_get(mt_dbl,n_munu,'munu',l_munu,k_munu))
     &    call errquit('hnd_elfcon: ma-munu 3 failed',911,MA_ERR)
c    Reset arr dbl_mb(k_munu+count)
       call dcopy(n_munu,0.0d0,0,dbl_mb(k_munu),1)
c ------- create munu-scratch-array--- END
c
c     Zero elfval result array
c
      call dcopy(max(nder*3,1)*npt,0.0d0,0,elfval,1)
c
c     ----- calculate electronic integral component(s) at all points -----
c
c =============== Part 1: main diagonal elements ============ START
      ijshell = 0
      do ishell = 1, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)

c         write(*,17) ishell,ilo,ihi
c 17      format('(ishell,ilo,ihi)=(',i8,',',i8,',',i8,')')

         idim = ihi - ilo + 1
         jshell=ishell
         ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
               jhi=ihi
               jlo=ilo
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c     Get the density patch, make the integrals and contract
c
               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $                     dbl_mb(k_dens), idim)
c
               call int_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt) 
               call multi_reduce(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           elfval,idim,jdim,
     &                           npt*(max(nder*3,1)))
c ================= TEST-copy-munu === START
               call copy_munu_mat(dbl_mb(k_buf),
     &                            ndir,npt,
     &                            dbl_mb(k_AtNr),nlist,
     &                            idim,jdim,npt*(max(nder*3,1)),
     &                            ilo,ihi,jlo,jhi,
     &                            dbl_mb(k_munu),nbf)
c ================= TEST-copy-munu === END
            end if  ! mod parallel loop
      end do    ! ishell
c =============== Part 1: main diagonal elements ============ END
c =============== Part 2: Off-main diagonal elements ======== START
      ijshell = 0
      do ishell = 2, nshell
c
c     get basis info
c
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, ishell-1
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info
c
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi)) call
     &            errquit('hnd_elfcon: bas_cn2bfr',basis,BASIS_ERR)
               jdim = jhi - jlo + 1
               nint = idim * jdim

c         write(*,18) ishell,ilo,ihi,
c     &               jshell,jlo,jhi
c 18      format('off:(ishell,ilo,ihi)=(',i8,',',i8,',',i8,')',
c     &             ' (jshell,jlo,jhi)=(',i8,',',i8,',',i8,')')
c
c     Get the density patch, make the integrals and contract
c
               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $                     dbl_mb(k_dens), idim)
c
               call int_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt)
c
               call multi_reduce_off(dbl_mb(k_buf),dbl_mb(k_dens),
     &                               elfval,idim,jdim,
     &                               npt*(max(nder*3,1)))
c ================= TEST-copy-munu === START
               call copy_munu_mat(dbl_mb(k_buf),
     &                            ndir,npt,
     &                            dbl_mb(k_AtNr),nlist,
     &                            idim,jdim,npt*(max(nder*3,1)),
     &                            ilo,ihi,jlo,jhi,
     &                            dbl_mb(k_munu),nbf)
c ================= TEST-copy-munu === END
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c     Collect components from all the nodes for all points
c
      call ga_sync()
      call ga_dgop(msg_efgs_col,elfval,npt*(max(nder*3,1)),'+')
      call ga_dgop(msg_efgs_col,dbl_mb(k_munu),n_munu,'+')

c ------ Collect from nodes munu matrix ---- START
c ---- PRINTOUT selected munu matrix ---- START
        n_munu1=nbf*(nbf+1)/2*ndir*nlist
c       write(*,*) "n_munu1=",n_munu1
c --------- create g_munuV6 -------- START
         if (.not. ga_create(mt_dbl,1,n_munu1,
     &                       'hnd_elfcon_symm: g_munu',
     $                       0,0,g_munuV6))
     $       call errquit('hnd_elfcon_symm:', 0,
     &                    GA_ERR)
        call ga_zero(g_munuV6)
c --------- create g_munuV6 -------- END
c +++++++++ ONLY FOR TEST get dens ++ START
       maxbuf=nbf*nbf
       if (.not. ma_push_get(mt_dbl,maxbuf,'dens1',l_dens1,k_dens1))
     &    call errquit('hnd_elfcon: ma 3 failed',911,MA_ERR)
c        write(*,*) "----- g_dens-check ---- START"
c        call ga_print(g_dens)
c        write(*,*) "----- g_dens-check ---- END"
       call ga_get(g_dens, 1, nbf, 1, nbf,
     $             dbl_mb(k_dens1), nbf)
c +++++++++ ONLY FOR TEST get dens ++ END

       count=1
       do ii=1,nlist
        do k=1,ndir
          indx1=nbf*nbf*ndir*(ii-1)+
     &          nbf*nbf*(k-1)
c +++++++++++ CALC-EFGs +++++ to check +++ START 
c        ac=0.0d0
c         do i=1,nbf
c          do j=1,nbf
c           indx=indx1+nbf*(i-1)+j-1
c           indx2=nbf*(i-1)+j-1
c           ac=ac+dbl_mb(k_munu+indx)*dbl_mb(k_dens1+indx2)      
c          enddo
c         enddo
c      write(*,22) ii,k,ac
c 22   format('EFG-(',i8,',',i8,')=',f15.8)
c +++++++++++ CALC-EFGs +++++ to check +++ END
          ac1=0.0d0
          do i=1,nbf
           indx=indx1+nbf*(i-1)+i-1
c ======= ONLY-TEST==START
c           indx2=nbf*(i-1)+i-1
c           ac1=ac1+dbl_mb(k_munu+indx)*dbl_mb(k_dens1+indx2)
c ======= ONLY-TEST==END

c           write(*,23) ii,k,i,i,dbl_mb(k_munu+indx)
c 23        format('DIAG:arr_munu(',i8,',',
c     &            i8,',',i8,',',i8,')=',f15.8)
           val=dbl_mb(k_munu+indx)
           call ga_fill_patch(g_munuV6,1,1,count,count,val)
           count=count+1
          enddo ! end-loop-i
          do i=2,nbf
           do j=1,i-1
            indx=indx1+nbf*(i-1)+j-1
c ======= ONLY-TEST==START
c           indx2=nbf*(i-1)+j-1
c           ac1=ac1+2.0d0*dbl_mb(k_munu+indx)*dbl_mb(k_dens1+indx2)
c ======= ONLY-TEST==END

c            write(*,22) ii,k,i,j,dbl_mb(k_munu+indx)
c 22         format('OFF:arr_munu(',i8,',',i8,',',
c     &             i8,',',i8,')=',f15.8)
            val=dbl_mb(k_munu+indx)
            call ga_fill_patch(g_munuV6,1,1,count,count,val)
            count=count+1
           enddo ! end-loop-j
          enddo ! end-loop-i

c      write(*,27) ii,k,ac1
c 27   format('EFG-symm-(',i8,',',i8,')=',f15.8)

        enddo ! end-loop-k
       enddo ! end-loop-ii
c        write(*,*) "---------- g_munuV6------------- START"
c        call ga_print(g_munuV6)
c        write(*,*) "---------- g_munuV6------------- END"

c       stop

c ---- PRINTOUT selected munu matrix ---- END
c =============== Part 2: Off-main diagonal elements ======== END
c
        if (.not.ma_pop_stack(l_dens1)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_dens1 failed',911,MA_ERR)


        if (.not.ma_pop_stack(l_munu)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_munu failed',911,MA_ERR)
c     write(*,*) "AFTER-cleaning-1"
c     Clean up MA data blocks
c
      if (.not.ma_pop_stack(l_buf)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_buf failed',911,MA_ERR)
c     write(*,*) "AFTER-cleaning-2"
      if (.not.ma_pop_stack(l_scr)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_scr failed',911,MA_ERR)
c     write(*,*) "AFTER-cleaning-3"
       if (.not.ma_pop_stack(l_dens)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_dens failed',911,MA_ERR)
c     write(*,*) "AFTER-cleaning-4"
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('dft_zora_utils: ma_free_heap l_AtNr',0, MA_ERR)

      return
      end

      subroutine multi_reduce_off(mblock,block,rblock,idim,jdim,nblock)
c
      implicit none
      integer idim,jdim,nblock
      double precision mblock(idim,jdim,nblock), block(idim,jdim)
      double precision rblock(nblock)
      double precision factor
      integer iblock,i,j

      factor=2.0d0
c
      do iblock = 1, nblock
         do i = 1, idim
            do j = 1, jdim
               rblock(iblock)=rblock(iblock)+mblock(i,j,iblock)*
     &                                       block(i,j)*factor
            enddo
         enddo
      enddo
c
      return
      end

      subroutine copy_munu_mat(mblock,
     &                         ndir,nat,
     &                         atm_list,nlist,
     &                         idim,jdim,nblock,
     &                         ilo,ihi,jlo,jhi,
     &                         arr_munu,nbf)
cc     dir_slc, 1 -> xx, 2 -> yy, 3 -> zz
c               4 -> xy, 5 -> xz, 6 -> yz
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
      
      integer atm_slc,dir_slc,nbf
      integer ilo,ihi,jlo,jhi
      integer iblock_slc,ndir
      integer idim,jdim,nblock
      integer nlist
      double precision mblock(idim,jdim,nblock)
      double precision arr_munu(nbf,nbf,ndir,nlist)
      integer i,j,k,ii,jj
      integer iat,nat
      double precision atm_list(nat) 
c
      do iat=1,nlist
       atm_slc=atm_list(iat) 
       do dir_slc=1,ndir
        iblock_slc=ndir*(atm_slc-1)+dir_slc
c        write(*,110) atm_slc,dir_slc,iblock_slc
c 110    format('(atm_slc,dir_slc,iblock_slc)=(',
c     &          i8,',',i8,',',i8,')')
       ii=1
         do i = ilo, ihi
           jj=1
           do j = jlo, jhi
            arr_munu(i,j,dir_slc,iat)=
     &           mblock(ii,jj,iblock_slc) ! gives upper-triang
            arr_munu(j,i,dir_slc,iat)=
     &           mblock(ii,jj,iblock_slc) ! gives lower-triang
            jj=jj+1
           enddo ! end-loop-j
           ii=ii+1
         enddo ! end-loop-i
       enddo ! end-loop-dir_scl
      enddo ! end-loop-iat
      return
      end


