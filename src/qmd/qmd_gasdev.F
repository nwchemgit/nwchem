c
      double precision function qmd_gasdev(idum)
c
      implicit none
c
      integer idum,iset
c
      double precision fac,gset,rsq,v1,v2
c
      external qmd_ran1
      double precision qmd_ran1
c
      save iset,gset
      data  iset/0/
c
      if (iset.eq.0) then
1       v1=2.d0*qmd_ran1(idum)-1.
        v2=2.d0*qmd_ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=dsqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        qmd_gasdev=v2*fac
        iset=1
      else
        qmd_gasdev=gset
        iset=0
      endif
c
      return
      end
