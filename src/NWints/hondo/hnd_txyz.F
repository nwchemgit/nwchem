c
c $Id$
c
      subroutine hnd_txyz(xt,yt,zt,xs,ys,zs,ni,nj,aj)
      implicit none
c
      double precision xt, yt, zt, xs, ys, zs, aj
      integer ni, nj, i, j
      dimension xt(ni,*),yt(ni,*),zt(ni,*)
      dimension xs(ni,*),ys(ni,*),zs(ni,*)
c
      do 10 i=1,ni
      xt(i,1)=(xs(i,1  )            -xs(i,3  )*(aj+aj))*aj
      yt(i,1)=(ys(i,1  )            -ys(i,3  )*(aj+aj))*aj
      zt(i,1)=(zs(i,1  )            -zs(i,3  )*(aj+aj))*aj
   10 continue
c
      if(nj.eq.1) return
c
      do 20 i=1,ni
      xt(i,2)=(xs(i,2  )*dble(2+2-1)-xs(i,4  )*(aj+aj))*aj
      yt(i,2)=(ys(i,2  )*dble(2+2-1)-ys(i,4  )*(aj+aj))*aj
      zt(i,2)=(zs(i,2  )*dble(2+2-1)-zs(i,4  )*(aj+aj))*aj
   20 continue
c
      if(nj.eq.2) return
c
      do 30 j=3,nj
      do 30 i=1,ni
      xt(i,j)=(xs(i,j  )*dble(j+j-1)-xs(i,j+2)*(aj+aj))*aj
     1        -xs(i,j-2)*dble(((j-1)*(j-2))/2)
      yt(i,j)=(ys(i,j  )*dble(j+j-1)-ys(i,j+2)*(aj+aj))*aj
     1        -ys(i,j-2)*dble(((j-1)*(j-2))/2)
      zt(i,j)=(zs(i,j  )*dble(j+j-1)-zs(i,j+2)*(aj+aj))*aj
     1        -zs(i,j-2)*dble(((j-1)*(j-2))/2)
   30 continue
c
      return
      end
