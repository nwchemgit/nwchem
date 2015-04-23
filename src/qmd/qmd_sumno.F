c
      double precision function qmd_sumno(nn,idum)
c
      implicit none
c
      integer :: nn, idum
c
      external qmd_gasdev
      double precision qmd_gasdev
      external qmd_gamdev
      double precision qmd_gamdev
c
c returns the sum of n independent gaussian noises squared
c (i.e. equivalent to summing the square of the return values of nn
c calls to gasdev)
c      
      if(nn==0) then
        qmd_sumno=0.d0
      else if(nn==1) then
        qmd_sumno=qmd_gasdev(idum)**2
      else if(modulo(nn,2)==0) then
        qmd_sumno=2.d0*qmd_gamdev(nn/2,idum)
      else
        qmd_sumno=2.d0*qmd_gamdev((nn-1)/2,idum) + qmd_gasdev(idum)**2
      end if
      return
      end
