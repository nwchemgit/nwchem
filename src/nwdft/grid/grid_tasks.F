      subroutine grid_tasks(rtdb, coord, qsh, nqlen,
     &                      aqbuf, npts_per_ang,nq_tot)
c
c$Id: grid_tasks.F,v 1.18 2000-03-18 02:14:59 edo Exp $
c
      implicit none
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      integer rtdb, qsh, nqlen, npts_per_ang, nq_tot
      integer lrij, irij, lrq, irq, lp, ip, lexpo, iexpo
      integer ndbl
      double precision coord, aqbuf, rad,rad_sh,
     ,     rmax,rmin,rad_write
      dimension coord(3,ncenters)
      dimension qsh(3,nqlen)
      dimension aqbuf(4,npts_per_ang)
c
      logical oreadfile, owritefile
      integer filesize, memsize
      character*255 gridfilename
      character*4 whichf
c
      integer iqsh, irsh, ictr, iang, npts_iqsh
      integer nq,nproc,me,g_buck,g_buckind,
     ,     k_qat,l_qat,nqsh_at,iii,i_a,ncu,
     ,     nx,ny,nz,fnxyz,nxyz,nqtotc,nq_cube,ncubenz

      logical trunc_neigh
      logical grid_file_open
      external grid_file_open
      logical grid_file_close
      external grid_file_close
      logical grid_test_mem
      external grid_test_mem
c
c     input   
c        qsh(1,iqsh) = irsh | where within nrad
c        qsh(2,iqsh) = ictr | current center
c        qsh(3,iqsh) = iang | angular quad for given irsh radial
c
c     Open grid pts file
c
      call grid_packing_info
c
c     temp set gridfilename to gridpts
c
c     Grid_pts caching/file information (from scf_get_info)
c
      if (.not. rtdb_get(rtdb, 'grid:filesize',
     $     MT_INT, 1, filesize)) filesize = 0
      write(0,*) ' filesize ',filesize
      if (.not. rtdb_get(rtdb, 'grid:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
      write(0,*) ' memsize ',memsize
c
c     The opening routine will put the .pid on the grid_pts filename
c     (hence even tho' parallel file, open as sequential)
c
      if (.not. rtdb_cget(rtdb, 'grid:filename', 1, gridfilename)) 
     $   call util_file_name('gridpts',.true.,.false.,gridfilename)
      if (.not. rtdb_cput(rtdb, 'grid:filename', 1, gridfilename)) 
     $   call errquit('grid_tasks: trouble putting gridfilename', 0)
c
      oreadfile = .false.
      owritefile= .false.
c
      if (filesize.gt.0 .or. memsize.gt.0) then
         if (.not. grid_file_open(gridfilename, memsize, filesize))
     $        call errquit('grid_tasks: opening gridpts?', 0)
      endif
c
      owritefile= grid_test_mem (rtdb, geom, nq_tot)
c
c
c     get partitioning function
c

      if(lssw) do_interp = xc_interp(whichf)
c
c     Allocate some scratch space for construction of the spatial weights
c
      ndbl = (ncenters*(ncenters+1))/2
      if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &   call errquit('  in grid_tasks: cannot allocate Rij',0)
C      call a_dist(coord, Dbl_MB(iRij), ncenters)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'rq',lrq,irq))
     &   call errquit('grid_tasks: cannot allocate rq',0)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'p',lp,ip))
     &   call errquit('grid_tasks: cannot allocate rq',0)
      if (.not.MA_Push_get(mt_dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_tasks: cannot allocate rq',4)
      if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &   call errquit('grid_tasks: cannot allocate qat',5)
c      
      nproc=ga_nnodes()
      me=ga_nodeid()
c
      if (.not.rtdb_get(rtdb, 'dft:trunc_neigh', mt_log, 1,
     &   trunc_neigh))trunc_neigh=.true.

      call grid_buckinit(g_buck,g_buckind)

      nqtotc=0
      ncubenz=0
      do i_a =1,ncenters
c
c     get shells centers on i_a
c
         nqsh_at=0
         do iqsh=1,nqshells
            if(qsh(2,iqsh).eq.i_a) then
               nqsh_at=nqsh_at+1
               int_mb(k_qat+nqsh_at-1)=iqsh
            endif
         enddo
c
c     get no of segment
c
         rad=dble(rad_cutoff(1,iatype(i_a)))
         nxyz=nint((rad*2d0-l_rect)/(l_rect*2d0))
c
c     loop over cubes centers on i_a
c
         do ncu=1,ncubes_type(iatype(i_a))
            fnxyz=int_mb(k_nxyz+(iatype(i_a)-1)*ncubesmx+ncu-1)
cedo            write(0,*) ' ncu ',ncu,nxyz
            call n2nxyz(nxyz,fnxyz,nx,ny,nz)
            
!            write(0,*) ' nxyz ',fnxyz,nx,ny,nz
            call grid_rmaxmin(nx,ny,nz,rmax,rmin,l_rect)
            call dfill(4*npts_per_ang, 0.0d0, aqbuf, 1)
            rad_write=rmax
            nq_cube=0
         do 100 iii = 1, nqsh_at
            iqsh=int_mb(k_qat+iii-1)          
c

            irsh = qsh(1,iqsh)         
            rad_sh=rpts(irsh,iatype(i_a))
c
c     check if radius between rmin and rmax
c            
            if(rad_sh.ge.rmin.and.rad_sh.le.rmax) then
cedo               rad_write=max(rad_write,rad_sh)
               ictr = qsh(2,iqsh)
               iang = qsh(3,iqsh)
               nq = npts_per_ang
c     
c     Generate all angular quadrature points for this radial shell
c     
c     
            call grid_qpts(rad_sh,coord(1,ictr), ictr, irsh, iang, 
     ,              aqbuf, nq, 
     ,           nx,ny,nz)
cedo            if(nq.ne.0) write(0,*) ' nq ',nq
c     
c     modify the weight to include the spatial weight (and symmetry factor)
c     
            call grid_setspac(rtdb, coord, aqbuf,
     ,           dbl_mb(irij), 
     &           dbl_mb(irq), dbl_mb(ip), ictr, nq,
     ,           dbl_mb(iexpo), rad_sh, trunc_neigh)
            npts_iqsh = nq
            
            
c     
c     since only generating symmetry unique quadrature 
c     scale weights with symmetry factor if desired
c     
            if(npts_iqsh.ne.0) then
               call grid_sym(ictr, aqbuf, npts_iqsh)
c            
c     now pack this aqbuf into the write buffer (write and update task counter as needed)
c     
c               write(0,*) ' nq in cube ',npts_iqsh,rad_sh
               call grid_file_write(npts_iqsh, ictr, rad_write, aqbuf)
               
               nqtotc=nqtotc+npts_iqsh
               nq_cube=nq_cube+npts_iqsh
            endif
         endif
 100     continue
c         if(nq_cube.ne.0)write(0,*) ' nq_cube ',nq_cube,rad_write
         if(nq_cube.ne.0) ncubenz=ncubenz+1
         enddo ! end loop over cubes
      enddo !end loop over atoms

      write(0,*) ' NQTOTC NQ/cube',nqtotc,nqtotc/ncubenz
      call grid_buckend(g_buck,g_buckind)
!      if(.true.) call errquit(' out ',1948)
c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      call grid_file_write(nq, n_tot_tasks, rad, aqbuf)
c      write(6,*)' n_tot_tasks = ', n_tot_tasks
c
c     clean up
c
      if (.not. MA_Pop_Stack(l_qat))
     &   call errquit('grid_tasks: pop stack failed.',1)
      if (.not. MA_Pop_Stack(lexpo))
     &   call errquit('grid_tasks: pop stack failed.',1)
      if (.not. MA_Pop_Stack(lp))
     &   call errquit('grid_tasks: pop stack failed.',2)
      if (.not. MA_Pop_Stack(lrq))
     &   call errquit('grid_tasks: pop stack failed.',3)
      if (.not. MA_Pop_Stack(lrij))
     &   call errquit('grid_tasks: pop stack failed.',4)
      return
      end 
