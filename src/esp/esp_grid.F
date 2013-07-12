      subroutine esp_grid(x,q,g,v)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(3,natoms),q(natoms),g(3,mxgrid),v(mxgrid)
      real*8 d,dmin,rcut2,rad2
      integer i,ix,iy,iz,ia,j
c
      if(newgrd.eq.0) then
        if(me.eq.0) then
          call util_file_name('grid',.false.,.false.,grdfil)
          open(unit=lfngrd,file=grdfil,form='formatted',status='old')
          rewind(lfngrd)
          read(lfngrd,'(i10,f20.10)') ngrid,charge
          do 10 i=1,ngrid
            read(lfngrd,'(4f20.10)') (g(j,i),j=1,3),v(i)
   10     continue
          close(unit=lfngrd)
        endif
        call ga_brdcst(mre_013,v,ngrid*ma_sizeof(mt_dbl,1,mt_byte),0)
      else
c
c       determine actual number of grid points
c
        if(me.eq.0) then
          rcut2=rcut*rcut
          ngrid=1
          do 1 iz=0,mgrid(3)+1
            do 2 iy=0,mgrid(2)+1
              do 3 ix=0,mgrid(1)+1
                g(1,ngrid)=xmin(1)-rcut+dble(ix)*spac
                g(2,ngrid)=xmin(2)-rcut+dble(iy)*spac
                g(3,ngrid)=xmin(3)-rcut+dble(iz)*spac
                dmin=rcut2
                ia=1
                do 4 i=1,natoms
                  rad2=(cnm2au*radius(int(q(i))))**2
                  d=(g(1,ngrid)-x(1,i))*(g(1,ngrid)-x(1,i))
     +             +(g(2,ngrid)-x(2,i))*(g(2,ngrid)-x(2,i))
     +             +(g(3,ngrid)-x(3,i))*(g(3,ngrid)-x(3,i))
                  if(rad2-d.gt.small) goto 3
                  if(dmin-d.gt.small) dmin=d
    4           continue
                if(rcut2-dmin.gt.small) then
cx                write(*,'(i5,7f12.6)') ngrid,(x(i,1),i=1,3),
cx     +                   (g(i,ngrid),i=1,3),sqrt(dmin)
                  ngrid=ngrid+1
                  if(ngrid.gt.mxgrid)
     +               call md_abort('Not enough memory available',0)
                endif
    3         continue
    2       continue
    1     continue
          ngrid=ngrid-1
        endif
      endif
c
c     broadcast grid to other nodes
c
      if(np.gt.0) then
        call ga_brdcst(mre_001,ngrid,ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(mre_002,g,3*ngrid*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
c     determine gridpoints handled by current node
c
      if(np.gt.1) then
        ngfr=me*(ngrid/np)+1
        ngto=(me+1)*(ngrid/np)
        if(me.eq.np-1) ngto=ngrid
      else
        ngfr=1
        ngto=ngrid
      endif
c
      return
      end
