      subroutine hnd_xyzgia(gijkl,ng,nimax,njmax,nkmax,nlmax)
      implicit double precision (a-h,o-z)
      dimension gijkl(ng,3,3,nlmax,nkmax,njmax,nimax)
c
c     ----- shift down for x1,y1,z1 -giao- integrals -----
c
      do ni=1,nimax-1
         do nj=1,njmax
            do nk=1,nkmax
               do nl=1,nlmax
                  do ixyz=1,3
                     do ig=1,ng
      gijkl(ig,ixyz,2,nl,nk,nj,ni)=gijkl(ig,ixyz,1,nl,nk,nj,ni+1)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     ----- shift down for x2,y2,z2 -giao- integrals -----
c
      do ni=1,nimax
         do nj=1,njmax
            do nk=1,nkmax-1
               do nl=1,nlmax
                  do ixyz=1,3
                     do ig=1,ng
      gijkl(ig,ixyz,3,nl,nk,nj,ni)=gijkl(ig,ixyz,1,nl,nk+1,nj,ni)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c $Id$
