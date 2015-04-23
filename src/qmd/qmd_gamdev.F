c
      double precision function qmd_gamdev(nn,idum)
c
      implicit none
c
      integer idum,nn,j
c
      double precision am,e,s,v1,v2,x,y
c
      external qmd_ran1
      double precision qmd_ran1
c
      if(nn.lt.6)then
        x=1.d0
        do j=1,nn
          x=x*qmd_ran1(idum)
        end do
        x=-log(x)
      else
1       v1=2.d0*qmd_ran1(idum)-1.d0
        v2=2.d0*qmd_ran1(idum)-1.d0
       if(v1**2+v2**2.gt.1.d0)goto 1
        y=v2/v1
        am=nn-1
        s=sqrt(2.d0*am+1.d0)
        x=s*y+am
       if(x.le.0.d0)goto 1
        e=(1.d0+y**2)*exp(am*log(x/am)-s*y)
       if(qmd_ran1(idum).gt.e)goto 1
      end if
      qmd_gamdev=x
      return
      end
