      subroutine hnd_elfcon_symm(
     &         basis,     !  in: basis handle
     &         geom,      !  in: geom  handle
     &         g_dens,    !  in: electron density
     &         points,    !  in: (x,y,z) list of centers 
     &         npt,       !  in: total number of centers
     &         elfval,    ! out: EFG value list
     &         nder,      !  in: =2 means second derivative
     &         efgfile)   !  in: = 0 NOT NLMO/NBO =1 DO NLMO/NBO
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
c     Description: This routine is more efficient than
c                  hnd_elfcon() because it uses symmetry advantage
c                  and also, it generates
c                  g_munuV6, matrix neccessary in NLMO analysis
c     Author: Fredy W. Aquino
c     Date  : 06-22-11
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "zora.fh" 
      integer basis    ! [input] basis set
      integer geom     ! [input] geometry
      integer g_dens   ! [input] GA with density
      integer npt      ! [input] number of coord points
      integer nder     ! [input] electronic integral type
      integer efgfile  ! in: = 0 NOT NLMO/NBO =1 DO NLMO/NBO
      integer ishell, jshell, ijshell, nshell, nbf_max, me, nproc
      integer ilo, ihi, jlo, jhi, idim, jdim, nint
      integer l_dens, k_dens, 
     &        l_scr , k_scr, 
     &        l_buf , k_buf
      integer l_munu,k_munu,n_munu,n_munu1
      integer maxbuf, maxscr, i
      integer l_dens1,k_dens1
      integer count,j,k,nbf,ndir,indx,indx1,ii
      integer count2,indx2
      double precision points(3,npt) ! [input] coordinates for points
      double precision elfval(*)     ! [output] efg values for each coord
      double precision val
      double precision ac,ac1  
      external multi_reduce,    ! defined in hnd_eflcon.F
     &         multi_reduce_off,! for off-diag elements
     &         get_munuV6

      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     ----- calculate buffer and scratch space -----
c           buffer = (lmax*(lmax+1)/2)^2 * (max(nder*3,1) * ngridpoints
c           scratch = see hnd_elfder wrapper routine

      call int_init_1eelec(maxbuf,maxscr,basis,nder,npt)
      if (.not. bas_geom(basis, geom)) call errquit
     $   ('hnd_elfcon: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $   ('hnd_elfcon: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &   ('hnd_elfcon: bas_nbf_cn_max failed',555, BASIS_ERR)
      if (.not. bas_numbf(basis,nbf)) call errquit
     &   ('hnd_elfcon: bas_numbf failed',555, BASIS_ERR)

      if (.not. ma_push_get(mt_dbl,nbf_max*nbf_max,'dens patch',l_dens,
     &    k_dens)) call errquit('hnd_elfcon: ma 1 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxscr,'scratch',l_scr,k_scr))
     &    call errquit('hnd_elfcon: ma 2 failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxbuf,'int buf',l_buf,k_buf))
     &    call errquit('hnd_elfcon_symm: ma 3 failed',911,MA_ERR)
c ------- create munu-scratch-array--- START
      if (efgfile.eq.1) then ! Do NLMO/NBO analysis
       ndir=6 ! Nr. directions xx,yy,zz,xy,xz,yz
       n_munu=nbf*nbf*ndir*npt
       if (.not. ma_push_get(mt_dbl,n_munu,'munu',l_munu,k_munu))
     &    call errquit('hnd_elfcon: ma-munu 3 failed',911,MA_ERR)
c    Reset arr dbl_mb(k_munu+count)
       call dcopy(n_munu,0.0d0,0,dbl_mb(k_munu),1)
      endif
c ------- create munu-scratch-array--- END
c
c     Zero elfval result array

      call dcopy(max(nder*3,1)*npt,0.0d0,0,elfval,1)
c
c     ----- calculate electronic integral component(s) at all points -----
c
c =============== Part 1: main diagonal elements ============ START
      ijshell = 0
      do ishell = 1, nshell

c     get basis info

         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
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

               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $                     dbl_mb(k_dens), idim)

               call int_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt) 
               call multi_reduce(dbl_mb(k_buf),dbl_mb(k_dens),
     &                           elfval,idim,jdim,
     &                           npt*(max(nder*3,1)))
               if (efgfile.eq.1) then ! Do NLMO/NBO analysis
                call copy_munu_mat(dbl_mb(k_buf),
     &                             ndir,npt,
     &                             idim,jdim,npt*(max(nder*3,1)),
     &                             ilo,ihi,jlo,jhi,
     &                             dbl_mb(k_munu),nbf)
               endif
            end if  ! mod parallel loop
      end do    ! ishell
c =============== Part 1: main diagonal elements ============ END
c =============== Part 2: Off-main diagonal elements ======== START
      ijshell = 0
      do ishell = 2, nshell
c
c     get basis info

         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('hnd_elfcon: bas_cn2bfr failed for basis',basis,BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, ishell-1
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
c
c     get basis info

               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi)) call
     &            errquit('hnd_elfcon: bas_cn2bfr',basis,BASIS_ERR)
               jdim = jhi - jlo + 1
               nint = idim * jdim
c
c     Get the density patch, make the integrals and contract

               call ga_get(g_dens, ilo, ihi, jlo, jhi,
     $                     dbl_mb(k_dens), idim)

               call int_1eelec(basis,ishell,basis,jshell,maxscr,
     &                         dbl_mb(k_scr),nint,dbl_mb(k_buf),
     &                         nder,points,npt)

               call multi_reduce_off(dbl_mb(k_buf),
     &                               dbl_mb(k_dens),
     &                               elfval, ! out: EFG val
     &                               idim,jdim,
     &                               npt*(max(nder*3,1)))
               if (efgfile.eq.1) then ! Do NLMO/NBO analysis
                call copy_munu_mat(dbl_mb(k_buf),
     &                             ndir,npt,
     &                             idim,jdim,npt*(max(nder*3,1)),
     &                             ilo,ihi,jlo,jhi,
     &                             dbl_mb(k_munu),nbf)
               endif
            end if  ! mod parallel loop
         end do   ! jshell
      end do    ! ishell
c     Collect components from all the nodes for all points
      call ga_sync()
      call ga_dgop(msg_efgs_col,elfval,npt*(max(nder*3,1)),'+')
      if (efgfile.eq.1) then
       call ga_dgop(msg_efgs_col,dbl_mb(k_munu),n_munu,'+')
c =============== Part 2: Off-main diagonal elements ======== END
       call get_munuV6(dbl_mb(k_munu), ! in: array with matrices
     &                 n_munu,         ! in: size(vmunu)
     &                 ndir,           ! in: =6 = xx,yy,zz,xy,xz,yz
     &                 npt,            ! in:
     &                 nbf,            ! in: nr. basis functions
     &                 lbl_nlmo)       ! in: defined in zora.fh
        if (.not.ma_pop_stack(l_munu)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_munu failed',911,MA_ERR)
      endif
c     Clean up MA data blocks

      if (.not.ma_pop_stack(l_buf)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_buf failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_scr)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_scr failed',911,MA_ERR)
       if (.not.ma_pop_stack(l_dens)) call errquit
     &   ('hnd_elfcon, ma_pop_stack of l_dens failed',911,MA_ERR)
      return
      end

      subroutine multi_reduce_off(mblock,block,rblock,idim,jdim,nblock)

      implicit none
      integer idim,jdim,nblock
      double precision mblock(idim,jdim,nblock), block(idim,jdim)
      double precision rblock(nblock)
      double precision factor
      integer iblock,i,j

      factor=2.0d0
      do iblock = 1, nblock
         do i = 1, idim
            do j = 1, jdim
               rblock(iblock)=rblock(iblock)+mblock(i,j,iblock)*
     &                                       block(i,j)*factor
            enddo
         enddo
      enddo
      return
      end

      subroutine copy_munu_mat(mblock,
     &                         ndir,nat,
     &                         idim,jdim,nblock,
     &                         ilo,ihi,jlo,jhi,
     &                         arr_munu, ! out :
     &                         nbf)      ! in  : nr of basis functions
cc     dir_slc, 1 -> xx, 2 -> yy, 3 -> zz
c               4 -> xy, 5 -> xz, 6 -> yz

      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"     
      integer dir_slc,nbf
      integer ilo,ihi,jlo,jhi
      integer iblock_slc,ndir
      integer idim,jdim,nblock
      integer iat,nat
      double precision mblock(idim,jdim,nblock)
      double precision arr_munu(nbf,nbf,ndir,nat)
      integer i,j,k,ii,jj

      do iat=1,nat
       do dir_slc=1,ndir
        iblock_slc=ndir*(iat-1)+dir_slc
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

      subroutine get_munuV6(vmunu,   ! in: array with matrices
     &                      n_munu,  ! in: size(vmunu)
     &                      ndir,    ! in: =6 = xx,yy,zz,xy,xz,yz
     &                      nlist,   ! in:
     &                      nbf,     ! in: nr. basis functions
     &                      lbl_nlmo)! in: defined in zora.fh
      ! g_munuV6 contains unique munu-EFG
      ! in sequence: xx,yy,zz,xy,xz,yz
      ! each chunk is nbf*(nbf+1)/2 in size
      implicit none
#include "nwc_const.fh" 
#include "errquit.fh" 
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
c ------------- added definitions -------- START
       character*(*) lbl_nlmo
       character*255 zorafilename
       logical dft_zoraEFGZ4_NLMOAnalysis_write
       external dft_zoraEFGZ4_NLMOAnalysis_write,
     &          util_file_name
c ------------- added definitions -------- END
       integer ii,i,j,k,indx,indx1
       integer ndir,nlist,nbf,count,n_munu,n_munu1
       double precision val,vmunu(n_munu)
       integer g_munuV6,g_munu_rhoS,g_dens,
     &         g_zora_scale_munu(2),ipolmunu
c ------ Collect from nodes munu matrix ---- START
c ---- PRINTOUT selected munu matrix ---- START
       n_munu1=nbf*(nbf+1)/2*ndir*nlist
c --------- create g_munuV6 -------- START
         if (.not. ga_create(mt_dbl,1,n_munu1,
     &        'hnd_elfcon_symm: g_munu',0,0,g_munuV6))
     $    call errquit('hnd_elfcon_symm:',0,GA_ERR)
        call ga_zero(g_munuV6)
c --------- create g_munuV6 -------- END
       count=1
       do ii=1,nlist
        do k=1,ndir
          indx1=nbf*nbf*ndir*(ii-1)+
     &          nbf*nbf*(k-1)
          do i=1,nbf
           indx=indx1+nbf*(i-1)+i-1
           val=vmunu(indx+1)
           call ga_fill_patch(g_munuV6,1,1,count,count,val)
           count=count+1
          enddo ! end-loop-i
          do i=2,nbf
           do j=1,i-1
            indx=indx1+nbf*(i-1)+j-1
            val=vmunu(indx+1)
            call ga_fill_patch(g_munuV6,1,1,count,count,val)
            count=count+1
           enddo ! end-loop-j
          enddo ! end-loop-i
        enddo ! end-loop-k
       enddo ! end-loop-ii

c       if (ga_nodeid().eq.0)
c     &   write(*,*) '--------- g_munuV6 ---------- START'
c        call ga_print(g_munuV6)
c        if (ga_nodeid().eq.0)
c     &   write(*,*) '--------- g_munuV6 ---------- END'
        
c ---------> Write NMLO analysis data: 3 of 3 ----- START
        call util_file_name(lbl_nlmo,.false.,.false.,zorafilename)
        if (.not.dft_zoraEFGZ4_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz
     &              nlist, ! in: list of selected atoms 
     &                  3, ! in: writing order =1,2,3
     &           ipolmunu, ! in: write for ndada=1
     &  g_zora_scale_munu, ! in: write for ndada=1
     &        g_munu_rhoS, ! in: write for ndata=2
     &             g_dens, ! in: write for ndata=2
     &           g_munuV6))! in: write for ndata=3
     &     call errquit('get_rhoS: dft_zoraNLMO_write failed',
     &                  0,DISK_ERR)
c ---------> Write NMLO analysis data: 3 of 3 ----- END
        if (.not. ga_destroy(g_munuV6)) call errquit(
     &    'dft_zora_utils: ga_destroy failed ',0, GA_ERR)
      return
      end
