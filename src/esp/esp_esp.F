      subroutine esp_esp(ene,occ,dns,x,q,grid,val)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
#include "schwarz.fh"
#include "stdio.fh"
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      real*8 ene(2*nbf),occ(2*nbf),dns(mbf,mbf)
      real*8 x(3,natoms),q(natoms)
      real*8 grid(3,mxgrid),val(mxgrid)
c
      integer i,j,ish,jsh,ix
      integer ilo,ihi,jlo,jhi
      real*8 dist,vt
      integer imin,nval,jshi
      real*8 dmin,fact
      character*10 today,now
c
      integer maxbuf,maxscr,avail
      integer batch,no_batches,ngrid_batches
      integer iptr
      logical ldummy(1)
c
c     get the density matrix, occupation and energies
c
      call esp_denmat(occ,ene)
c
      if (lscrn.ne.0) call schwarz_init(igeom,ibasis)
c
c     get electrostatic potential on the grid points
c
      do i=1,ngrid
         val(i)=0.0d0
      end do
c
c     estimate mem. req.
c
      call int_init_1eelec(maxbuf,maxscr,ibasis,0,ngrid)
      avail=ma_inquire_avail(mt_dbl)*8/10
      no_batches=nint(dble(maxbuf+maxscr)/dble(avail))
      no_batches=max(no_batches,1)
      ngrid_batches=ngrid/no_batches
#ifdef DEBUG
      if(me.eq.0) then
      write(6,*) ' avail ',avail
      write(6,*) ' no_batches ',no_batches
      write(6,*) ' ngrid_batches ',ngrid_batches
      endif
#endif
      iptr=1
      do batch=1,no_batches
         if(batch.eq.no_batches) ngrid_batches=ngrid-iptr+1
         call hnd_elfcon_schw(ibasis,igeom,lg_d,
     C        grid(1,iptr),ngrid_batches,val(iptr),0)
         iptr=iptr+ngrid_batches
      enddo
      if(me.eq.0) then
      do 8 i=1,ngrid
         imin=0
         val(i)=-val(i)
         vt=val(i)
         do 9 ix=1,natoms
            dist=sqrt((grid(1,i)-x(1,ix))*(grid(1,i)-x(1,ix))+
     +          (grid(2,i)-x(2,ix))*(grid(2,i)-x(2,ix))+
     +          (grid(3,i)-x(3,ix))*(grid(3,i)-x(3,ix)))
            if(imin.eq.0) then
               dmin=dist
               imin=ix
            endif
            if(dmin.gt.dist) then
               dmin=dist
               imin=ix
            endif
            val(i)=val(i)+q(ix)/dist
    9    continue
    8 continue
      endif
c
      if(lscrn.ne.0) call schwarz_tidy()
c
      if(np.gt.1) then
        call ga_brdcst(mre_002,val,ngrid*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
c     deallocate memory for the integrals
c
      if(.not.ga_destroy(lg_d))
     + call md_abort('esp_denmat: ga_destroy lg_d failed',me)
c
c     integral termination
c
      call int_terminate()
c
      if(me.eq.0) then
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
      rewind(lfngrd)
      write(lfngrd,'(i10,f20.10)') ngrid,charge
      do i=1,ngrid
        write(lfngrd,'(4f20.10)') (grid(j,i),j=1,3),val(i)
      enddo
      close(unit=lfngrd)
      endif
c
      return
      end
