      subroutine grid_tasks(rtdb, xyz, qsh, nqlen,
     &                      aqbuf, npts_per_ang,nq_tot)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "sym.fh"
#include "geom.fh"
c
      integer rtdb, qsh, nqlen, npts_per_ang, nq_tot
      integer lrij, irij, lrq, irq, lp, ip, lexpo, iexpo
      integer ndbl
      double precision xyz, aqbuf, rad,rad_sh,
     ,     rmax,rmin,rad_write,symfact,enuc
      dimension xyz(3,ncenters)
      dimension qsh(3,nqlen)
      dimension aqbuf(4,npts_per_ang)
c
      logical oreadfile, owritefile
      integer filesize, memsize
      character*255 gridfilename
c
      integer iqsh, irsh, ictr, iang, npts_iqsh
      integer nq,nproc,me,
     ,     k_qat,l_qat,nqsh_at,iii,i_a,ncu,
     ,     fnxyz,nxyz,nqtotc,nq_cube,ncubenz,iaia
      logical grid_file_open
      external grid_file_open
      logical grid_test_mem
      external grid_test_mem
      integer nxtask,n1,n2,icount
      external nxtask
      integer iang_done
      integer nops ! [out]
      integer k_symops(2),l_symops(2) ! [out]  symops handles
      integer k_gridscr,l_gridscr ! [out]  grid handles
      integer k_iscr,l_iscr ! [out] iscr handles
      integer nq_left
      integer eaf_size_in_mb
      integer grid_nxyz
      external grid_nxyz
      logical  grid_file_write, grid_written
      external grid_file_write
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
      if (.not. rtdb_get(rtdb, 'grid:memsize',
     $     MT_INT, 1, memsize)) memsize = 0
      if (.not. rtdb_get(rtdb, 'grid:eaf_size_in_mb',
     $     MT_INT, 1, eaf_size_in_mb)) eaf_size_in_mb = 0
c
c     The opening routine will put the .pid on the grid_pts filename
c     (hence even tho' parallel file, open as sequential)
c
      call util_file_name('gridpts',.true.,.false.,gridfilename)
c
      oreadfile = .false.
      owritefile= .false.
c
      if (filesize.gt.0 .or. memsize.gt.0) then
         if(.not.geom_nuc_rep_energy(geom, Enuc)) 
     .        call errquit(' gridquad: enucrep failed ',0, GEOM_ERR)
         if (.not. grid_file_open(enuc,
     .        gridfilename, memsize, filesize,eaf_size_in_mb))
     $        call errquit('grid_tasks: opening gridpts?', 0, DISK_ERR)
      endif
c
      owritefile= grid_test_mem (rtdb,  nq_tot)
c
c
c     get partitioning function
c

c
c     Allocate some scratch space for construction of the spatial weights
c
      ndbl = (ncenters*(ncenters+1))/2
      if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &   call errquit('  in grid_tasks: cannot allocate Rij',0, MA_ERR)
C      call a_dist(xyz, Dbl_MB(iRij), ncenters)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'rq',lrq,irq))
     &   call errquit('grid_tasks: cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,ncenters*npts_per_ang,'p',lp,ip))
     &   call errquit('grid_tasks: cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_tasks: cannot allocate rq',4, MA_ERR)
      if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &   call errquit('grid_tasks: cannot allocate qat',5, MA_ERR)
c      
      nproc=ga_nnodes()
      me=ga_nodeid()
c
      call dft_tstamp(' Before buckinit.')

      call grid_buckinit(rtdb)

      nqtotc=0
      ncubenz=0
      call dft_tstamp(' Before atom loop')
c
c     parallel loop init
c
      icount = 1
      n1 = 0
      n2 = nxtask(nproc,icount)
      call grid_sympts0(geom,nops,k_symops,l_symops,
     ,     k_gridscr,l_gridscr,k_iscr,l_iscr,nqmax)
      do 111 i_a =1,ncenters
         if(oskel) then
            if (.not.sym_atom(geom, i_a, symfact))  goto 111
         endif
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
         iaia=iatype(i_a)
         rad=dble(rad_cutoff(1,iaia))
         nxyz=grid_nxyz(rad,lcube(iaia))
c
c     loop over cubes centers on i_a
c
         do ncu=1,ncubes_type(iatype(i_a))
            if (n1.eq.n2 )then
c     parallelization here ... might need file exchange when done
c     for load balancing
               fnxyz=int_mb(k_nxyz+(iatype(i_a)-1)*ncubesmx+ncu-1)
               fnxyz=int_mb(k_nxyz+(iatype(i_a)-1)*ncubesmx+ncu-1)
               call grid_rmaxmin(fnxyz,rmax,rmin,lcube(iaia),nxyz)
               nq_cube=0
               do 100 iii = 1, nqsh_at
                  iqsh=int_mb(k_qat+iii-1)          
                  irsh = qsh(1,iqsh)         
                  rad_sh=rpts(irsh,iatype(i_a))
c     
c     check if radius between rmin and rmax
c     
                  if(rad_sh.ge.rmin.and.rad_sh.le.rmax) then
!                     write(6,4321) iii,rad_sh
! 4321                format( ' Q ',I3,' RAD ',f12.6)
                     ictr = qsh(2,iqsh)
                     iang = qsh(3,iqsh)
                     nq = npts_per_ang
c     
c     Generate all angular quadrature points for this radial shell
c     and keep only those inside cube
c     
                     call grid_sympts(.false.,
     .                    rad_sh,xyz(1,ictr), ictr, 
     .                    irsh, iang, iang_done,nq_left,
     ,                    aqbuf, nq, fnxyz,nxyz,
     ,                    dbl_mb(k_gridscr),nops,dbl_mb(k_symops(1)),
     .                    dbl_mb(k_symops(2)),int_mb(k_iscr))
c     
c     modify the weight to include the spatial weight (and symmetry factor)
c     
                     if(ncenters.gt.1) then
                        call grid_setspac( xyz, aqbuf,
     ,                       dbl_mb(irij), 
     &                       dbl_mb(irq), dbl_mb(ip), ictr, nq,
     ,                       dbl_mb(iexpo), rad_sh)
                     endif
                     npts_iqsh = nq
c     
c     since only generating symmetry unique quadrature 
c     scale weights with symmetry factor if desired
c     
                     if(npts_iqsh.ne.0) then
                        if(oskel) call dscal(
     /                       npts_iqsh, symfact, aqbuf(4,1), 4)
c                        if(oskel) call grid_sym(ictr,aqbuf,npts_iqsh)
c     
c     now pack this aqbuf into the write buffer 
c     (write and update task counter as needed)
c     
                        grid_written = grid_file_write(npts_iqsh, ictr, 
     ,                       rad_write, aqbuf,  ncu)
                        
                        nqtotc=nqtotc+npts_iqsh
                        nq_cube=nq_cube+npts_iqsh
!                        write(6,1234) nx,ny,nz,iii,npts_iqsh,nqtotc
! 1234                   format( ' N ',3I3,' Q ',I3,' NP ',2I5)
                     endif
                  endif
 100           continue ! end loop over angular shells
               if(nq_cube.ne.0) ncubenz=ncubenz+1
               n1 = n1 + 1
               n2 = nxtask(nproc,icount)
            else
               n1 = n1 + 1
            endif
         enddo                  ! end loop over cubes
 111  continue ! endd loop over atoms
      call grid_sympts2(l_symops,l_gridscr)
      n1 = nxtask(-nproc,icount)
      call dft_tstamp(' done atom loop')

      write(LuOut,*) ' me ',me,
     ,     ' NQTOTC NQ/cube',nqtotc,nqtotc/ncubenz
      call grid_buckend()
!      if(.true.) call errquit(' out ',1948)
c
c     Write last buffer of grid pts and close file     
c
      nq = 1
      n_tot_tasks = 0 ! sets nctr to 0; writes last buf and returns #tasks
      grid_written = grid_file_write(nq, n_tot_tasks, rad, aqbuf,0)
c
c     clean up
c
      if (.not. MA_Pop_Stack(l_qat))
     &   call errquit('grid_tasks: pop stack failed.',1, MA_ERR)
      if (.not. MA_Pop_Stack(lexpo))
     &   call errquit('grid_tasks: pop stack failed.',1, MA_ERR)
      if (.not. MA_Pop_Stack(lp))
     &   call errquit('grid_tasks: pop stack failed.',2, MA_ERR)
      if (.not. MA_Pop_Stack(lrq))
     &   call errquit('grid_tasks: pop stack failed.',3, MA_ERR)
      if (.not. MA_Pop_Stack(lrij))
     &   call errquit('grid_tasks: pop stack failed.',4, MA_ERR)
      return
      end 
