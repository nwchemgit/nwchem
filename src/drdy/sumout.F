C $Id$
      Subroutine drdy_sumout(ns)
C
C Sort RPH interpolation information to give points along MEP in 
C  increasing value of s (iflg#0, final sort and print out summary)
C  
      Implicit none
      Integer ns
C
#include "drdyP.fh"
C
      Integer i,is,nfreq,ntot
C
      write (fu6,600)
      write (fu6,601) sr,vr*ckcal
      write (fu6,602) sp,vp*ckcal
      if (lgs(1).ne.0) write (fu6,603) 0.0d00,vspsv*ckcal,
     *   (xr(i,5)*cang/amass(i),i=1,n3)
      do is = 1,ns
         write (fu6,604) sgrid(is),vgrid(is)*ckcal,
     *   (xgrid(i,is)*cang/amass(i),i=1,n3)
      enddo
      ntot = nf(1)+nf(2)
      write (fu6,610)
      if (ntot.eq.0) then
         write (fu6,611) sr
      else 
         write (fu6,611) sr,(freqr(i,1),i=1,ntot)
      endif
      ntot = nf(3)+nf(4)
      if (ntot.eq.0) then
         write (fu6,613) sp
      else
         write (fu6,613) sp,(freqr(i,3),i=1,ntot)
      endif
      nfreq = nf(5)
      if (lgs(1).ne.0) write (fu6,615) 0.0d00,
     *   (freqr(i,5),i=1,nfreq+1)
      write (fu6,*) ' '
      do is = 1,ns
         write (fu6,616) sgrid(is),
     *    (freqg(i,is),i=1,nfreq)
      enddo
      return
600   Format (//,' Summary of MEP data',/,5x,'R = reactants, ',
     * 'P = products, * = saddle point or starting geometry',//,
     * '  Cartesian coordinates, x',//,t9,'s (au)',3x,'V (kcal/mol)',
     * 8x,'x (Angstroms)')
601   Format (' R',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
602   Format (' P',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
603   Format (' *',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
604   Format (t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
610   Format (//,'  Frequencies',//,
     *  t9,'s (au)',T37,'frequencies (cm^-1)')
611   Format (' R',t4,0pf12.5,(T35,0p9f9.2))
613   Format (' P',t4,0pf12.5,(T35,0p9f9.2))
615   Format (' *',t4,0pf12.5,(T35,0p9f9.2))
616   Format (t4,0pf12.5,(T35,0p9f9.2))
      End
