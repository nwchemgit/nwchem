      subroutine sinfo(
     $     ish, li, ci, ei, ngi, xi,
     $     jsh, lj, cj, ej, ngj, xj,
     $     ksh, lk, ck, ek, ngk, xk,
     $     lsh, ll, cl, el, ngl, xl,
     $     jtype)
c $Id$
      implicit none
c     
c     Takes info about centers and shells, precomputes some axis
c     transformation stuff and stuffs data into internal common blocks
c     
c     ish - shell no.,used only to determine if contraciton data is the 
c           same as in the immediately preceding call
c     li  - the highest angular momenteum of shell (s=0, p=1)
c     ci  - s and p contraction coefficients (unused p-coeffs must be zero)
c     ei  - exponents of primtives
c     ngi - no. of primitives
c     xi  - cartesian coords of atom
c     jtype - 1-6 specifies 0000,0001,0011,0101,0111,1111
c
c     Saving of data between calls is temporarily disabled ...
c     need to redo checking based on coordinates and unique
c     shell labels.
c
      integer ish, jsh, ksh, lsh, ngi, ngj, ngk, ngl
      integer li, lj, lk, ll
      double precision ci(ngi,2), cj(ngj,2), ck(ngk,2), cl(ngl,2)
      double precision ei(ngi), ej(ngj), ek(ngk), el(ngl)
      double precision xi(3), xj(3), xk(3), xl(3)
      integer jtype
c     
#include "maxc.fh"
#include "miscg.fh"
#include "shllfo.fh"
#include "spgeom.fh"
#include "qgeom.fh"
#include "const.fh"
#include "pgeom.fh"
#include "ginf.fh"
#include "picon.fh"
#include "errquit.fh"
c     
      double precision sixty, done, dzero, pt5, pt7, pt9, pt0001, tenm12
c     
      integer i, j, k, l, ind
      integer ijatom, klatom, ijshel, klshel
      logical oijat, oklat
      double precision qq, qqtest, qqqq, qqq, eab, gbeab
      double precision abx, aby, abz, cdx, cdy, cdz
      double precision ppq1, ppq2, ppq3, pmq1, pmq2, pmq3
      double precision p2, temp, p3333, sinp, p3131, csai, cpai
c
      integer ijshlo, ijato, klshlo, klato
      integer*4 jtype4
      save ijshlo, ijato, klshlo, klato ! STATIC variables
c     
      data ijshlo, ijato, klshlo, klato /4*0/
      data sixty/60.0d0/
      data dzero/0.0d0/,pt5/0.5d0/,pt7/0.7d0/,pt9/0.9d0/,done/1.0d0/
*
*rjh pt0001 was 1d-4 but this resulted in discontinuities in the
*rjh energy of order 1d-6 for displacements of order less than 1d-3 bohr
*rjh for certain orientations of molecules.  The smaller value pushes
*rjh the discontinuity to displacements of order 1d-5.
*
      data pt0001/1.0d-06/,tenm12/1.0d-12/
c     
c     Determine if atom/shell information has changed ... if so save info
c     
*      if ((max(iat,jat,kat,lat) .gt. 1000) .or.
*     $     (max(ish,jsh,ksh,lsh) .gt. 10000)) call errquit
*     $     ('sinfo: Wow, this calculation is really big!', 0)
c     
*      ijatom = iat*1000  + jat
*      klatom = kat*1000  + lat
*      ijshel = ish*10000 + jsh
*      klshel = ksh*10000 + lsh
      jtype4=jtype
      ijatom = 1
      klatom = 1
      ijshel = 1
      klshel = 1
c     
      if(ijshel .ne. ijshlo) then
         nga = ngi
         ngb = ngj
         la  = abs(li)
         lb  = abs(lj)
         ngangb = nga*ngb
         mab = la+lb-1
c
c     If sp shell then ci(*,1) = s coeffs, ci(*,2) = p coeffs
c     If s  shell then ci(*,1) = s coeffs, no p coeffs
c     If p  shell then ci(*,1) = p coeffs, no s coeffs
c
         if (li .eq. -1) then
            do i = 1, nga
               cmaxa(i) = max(abs(ci(i,1)),abs(ci(i,2)))
               ag(i)    = ei(i)
               csa(i)   = ci(i,1)
               cpa(i)   = ci(i,2)
            enddo
         else if (li .eq. 0) then
            do i = 1, nga
               cmaxa(i) = abs(ci(i,1))
               ag(i)    = ei(i)
               csa(i)   = ci(i,1)
               cpa(i)   = 0.0d0
            enddo
         else
            do i = 1, nga
               cmaxa(i) = abs(ci(i,1))
               ag(i)    = ei(i)
               csa(i)   = 0.0d0
               cpa(i)   = ci(i,1)
            enddo
         endif
         if (lj .eq. -1) then
            do j = 1, ngb
               cmaxb(j) = max(abs(cj(j,1)),abs(cj(j,2)))
               bg(j)    = ej(j)
               csb(j)   = cj(j,1)
               cpb(j)   = cj(j,2)
            enddo
         else if (lj .eq. 0) then
            do j = 1, ngb
               cmaxb(j) = abs(cj(j,1))
               bg(j)    = ej(j)
               csb(j)   = cj(j,1)
               cpb(j)   = 0.0d0
            enddo
         else
            do j = 1, ngb
               cmaxb(j) = abs(cj(j,1))
               bg(j)    = ej(j)
               csb(j)   = 0.0d0
               cpb(j)   = cj(j,1)
            enddo
         endif
      endif
c     
      if(klshel .ne. klshlo) then
*     klshlo = klshel
         ngc = ngk
         ngd = ngl
         lc  = abs(lk)
         ld  = abs(ll)
         if (lk .eq. -1) then
            do k = 1, ngc
               cmaxc(k) = max(abs(ck(k,1)),abs(ck(k,2)))
               cg(k)    = ek(k)
               csc(k)   = ck(k,1)
               cpc(k)   = ck(k,2)
            enddo
         else if (lk .eq. 0) then
            do k = 1, ngc
               cmaxc(k) = abs(ck(k,1))
               cg(k)    = ek(k)
               csc(k)   = ck(k,1)
               cpc(k)   = 0.0d0
            enddo
         else
            do k = 1, ngc
               cmaxc(k) = abs(ck(k,1))
               cg(k)    = ek(k)
               csc(k)   = 0.0d0
               cpc(k)   = ck(k,1)
            enddo
         endif
         if (ll .eq. -1) then
            do l = 1, ngd
               cmaxd(l) = max(abs(cl(l,1)),abs(cl(l,2)))
               dg(l)    = el(l)
               csd(l)   = cl(l,1)
               cpd(l)   = cl(l,2)
            enddo
         else if (ll .eq. 0) then
            do l = 1, ngd
               cmaxd(l) = abs(cl(l,1))
               dg(l)    = el(l)
               csd(l)   = cl(l,1)
               cpd(l)   = 0.0d0
            enddo
         else
            do l = 1, ngd
               cmaxd(l) = abs(cl(l,1))
               dg(l)    = el(l)
               csd(l)   = 0.0d0
               cpd(l)   = cl(l,1)
            enddo
         endif
      endif
c     
c     fill common misc
c     
      mcd = lc+ld-1
c     
c     Construct coordinate dependent info
c
      oijat = ijatom.eq.ijato
      oklat = klatom.eq.klato
      if(.not.(oijat.and.oklat)) then
         if( .not. oijat) then
            ax = xi(1)
            ay = xi(2)
            az = xi(3)
            bx = xj(1)
            by = xj(2)
            bz = xj(3)
            abx = bx-ax
            aby = by-ay
            abz = bz-az
            rabsq = abx**2+aby**2+abz**2
            rab = dsqrt(rabsq)
            if (rab .ne. 0.0d0) then
               p31 = abx/rab
               p32 = aby/rab
               p33 = abz/rab
            else
               p31 = dzero
               p32 = dzero
               p33 = done
            endif
*     ijato=ijatom
         endif
         if( .not. oklat ) then
            cx = xk(1)
            cy = xk(2)
            cz = xk(3)
            dx = xl(1)
            dy = xl(2)
            dz = xl(3)
            cdx = dx-cx
            cdy = dy-cy
            cdz = dz-cz
            rcdsq = cdx**2+cdy**2+cdz**2
            rcd = dsqrt(rcdsq)
            if (rcd .ne. 0.0d0) then
               q31 = cdx/rcd
               q32 = cdy/rcd
               q33 = cdz/rcd
            else
               q31 = dzero
               q32 = dzero
               q33 = done
            endif
*     klato=klatom
         endif
         cosg = p31*q31+p32*q32+p33*q33
         cosg = dmin1(done,cosg)
         cosg = dmax1(-done,cosg)
         sing = dsqrt(done-cosg*cosg)
         if (  dabs(cosg)-pt9) 300,300,240
 240     ppq1 = p31+q31
         ppq2 = p32+q32
         ppq3 = p33+q33
         pmq1 = p31-q31
         pmq2 = p32-q32
         pmq3 = p33-q33
         p21 = pmq2*ppq3-ppq2*pmq3
         p22 = pmq3*ppq1-ppq3*pmq1
         p23 = pmq1*ppq2-ppq1*pmq2
         p2 = dsqrt(p21*p21+p22*p22+p23*p23)
         sing = pt5*p2
         if (sing-tenm12) 280,260,260
 260     temp = done/p2
         p21 = p21*temp
         p22 = p22*temp
         p23 = p23*temp
         go to 380
 280     if (  dabs(p31)-pt7) 360,360,320
 300     p21 = (p32*q33-p33*q32)/sing
         p22 = (p33*q31-p31*q33)/sing
         p23 = (p31*q32-p32*q31)/sing
         go to 380
 320     p3333 = p33*p33
         p3333 = dmin1(done,p3333)
         sinp = dsqrt(done-p3333)
         p21 = p32/sinp
         p22 = -p31/sinp
         p23 = dzero
         go to 380
 360     p3131 = p31*p31
         p3131 = dmin1(done,p3131)
         sinp = dsqrt(done-p3131)
         p21 = dzero
         p22 = p33/sinp
         p23 = -p32/sinp
 380     q21 = p21
         q22 = p22
         q23 = p23
         p11 = p22*p33-p23*p32
         p12 = p23*p31-p21*p33
         p13 = p21*p32-p22*p31
         q11 = q22*q33-q23*q32
         q12 = q23*q31-q21*q33
         q13 = q21*q32-q22*q31
         acx = (cx-ax)*p11+(cy-ay)*p12+(cz-az)*p13
         acy = (cx-ax)*p21+(cy-ay)*p22+(cz-az)*p23
         if (  dabs(acy)-pt0001) 400,400,420
 400     acy = dzero
 420     continue
         acz = (cx-ax)*p31+(cy-ay)*p32+(cz-az)*p33
         acy2 = acy*acy
      endif
      if( ijshel.ne.ijshlo ) then
         ind=0
*     ijshlo=ijshel
         do 38000 i = 1,nga
            ga = ag(i)
            csai = csa(i)
            cpai = cpa(i)
            do 36000 j = 1,ngb
               ind = ind+1
               gb = bg(j)
               gab = ga+gb
               gp(ind) = gab
               eab = done/gab
               ep(ind) = eab
               gbeab = gb*eab
               app(ind) = gbeab*rab
               bpp(ind) = app(ind)-rab
               qqq = ga*gbeab*rabsq
               if (qqq-sixty) 18000,18000,10000
10000          ismlp(ind) = 2
               dp00p(ind) = dzero
#if 1
               jtype4=jtype-3
               if (jtype4) 36000,36000,12000
#else
               if (jtype-3) 36000,36000,12000
#endif
12000          dp01p(ind) = dzero
               conp(ind) = dzero
#if 1
               jtype4=jtype-5
               if (jtype4) 14000,14000,16000
#else
               if (jtype-5) 14000,14000,16000
#endif
14000          bpp(ind) = bpp(ind)*gab
               go to 36000
16000          dp10p(ind) = dzero
               dp11p(ind) = dzero
               go to 36000
18000          qq =  dexp(-qqq)*eab
               qqtest = cmaxa(i)*cmaxb(j)*qq
               if (qqtest-error1) 22000,22000,20000
20000          ismlp(ind) = 0
               go to 28000
22000          if (qqtest-error2) 26000,26000,24000
24000          ismlp(ind) = 1
               go to 28000
26000          ismlp(ind) = 2
28000          qqqq = pito52*qq
               dp00p(ind) = qqqq*csai*csb(j)
#if 1
               jtype4=jtype-3
               if (jtype4) 36000,36000,30000
#else
               if (jtype-3) 36000,36000,30000
#endif
30000          dp01p(ind) = qqqq*csai*cpb(j)
#if 1
               jtype4=jtype-5
               if (jtype4) 32000,32000,34000
#else
               if (jtype-5) 32000,32000,34000
#endif
32000          conp(ind) = dp01p(ind)*eab
               dp00p(ind) = dp00p(ind)*gab/dp01p(ind)
               bpp(ind) = bpp(ind)*gab
               go to 36000
34000          dp10p(ind) = qqqq*cpai*csb(j)
               dp11p(ind) = qqqq*cpai*cpb(j)
               conp(ind) = dp11p(ind)
               dp00p(ind) = dp00p(ind)/dp11p(ind)
               dp01p(ind) = dp01p(ind)/dp11p(ind)
               dp10p(ind) = dp10p(ind)/dp11p(ind)
36000       continue
38000    continue
      endif
      return
      end
