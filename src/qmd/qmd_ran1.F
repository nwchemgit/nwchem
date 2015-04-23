c
      double precision function qmd_ran1(idum)
c
      implicit none
c
      integer idum,IA,IM,IQ,IR,NTAB,NDIV
      double precision AM,EPS,RNMX
c
      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     & NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c
      integer j,k,iv(NTAB),iy
      save iv,iy
      data iv /NTAB*0/, iy /0/
c
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      qmd_ran1=min(AM*iy,RNMX)
c
      return
      end
