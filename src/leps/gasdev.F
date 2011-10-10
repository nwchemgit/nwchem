C     program gasdev1
C     double precision gks_gasdev,gasdev,seed
C      seed = -1234
C      write(*,*)  'before gasdev'
C      gasdev=gks_gasdev(seed)
C      write(*,*)  'Afer gasdev'
C     end
C    rundom number generators from gaussian distribution, 
C    from Numerical Recipies
C    it is initialized with -IDUM

      double precision FUNCTION GKS_GASDEV(IDUM)
      save
c     implicit undefined (a-z)
      integer idum,iset
      double precision v1,v2,r,fac,gset,gks_ran1
      DATA ISET /0/
c      write(*,*) 'inside GKS_GASDEV'
      IF(ISET.EQ.0) THEN
   10   V1=2.0E+00*GKS_RAN1(IDUM)-1.0E+00
        V2=2.0E+00*GKS_RAN1(IDUM)-1.0E+00
        R=V1**2+V2**2
        IF(R.GE.1.0E+00) GO TO 10
        FAC=SQRT(-2.0E+00*LOG(R)/R)
        GSET      =V1*FAC
        GKS_GASDEV=V2*FAC
        ISET=1
      ELSE
        GKS_GASDEV=GSET
        ISET=0
      ENDIF
c      write(*,*) idum,iset,gks_gasdev
C9999 format(' ---gasdev--- = ',2i15,f10.6)
      RETURN
      END
ccccccccc/ccccccccc/ccccccccc/ccccccccc/ccccccccc/ccccccccc/ccccccccc/ccccccccc/
      double precision function gks_ran1(idum)
c
c     N.R. FUNCTION gks_ran2(idum)
c
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     1   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     2   NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      gks_ran1=min(AM*iy,RNMX)
      return
      END
ccccccc

c $Id$
