      Subroutine hf1mkr_fast(Axyz,Bxyz,Cxyz,zan,ncenters,
     &                  alpha,Pxyz,RS,PC,ff,R,
     &                  R0,R0C,IJK,NPP,Lp,Lp3,CENTER)
c $Id$
c
c  Acknowledgment of AI-Assisted Development
c
c  Portions of this source code were developed with the assistance of
c  Claude (Anthropic), an AI assistant, specifically Claude Opus 5
c  (accessed via API/Claude.ai, 2026).
c     Centers are processed in batches of nc so that the vector length
c     seen by igamma, hfmkr and the local loops is MVEC = NPP*nc rather
c     than NPP.  This matters because NPP is frequently 1, in which case
c     the original code called igamma and hfmkr ncenters times on
c     vectors of length one.
c
c     ---------------------------------------------------------------
c     SCRATCH REQUIREMENT CHANGE.  The caller must allocate
c
c        R  :  MVMAXA*(Lp+1)*Lp3   doubles   (was NPP*(Lp+1)*Lp3)
c        PC :  MVMAXA*3            doubles   (was NPP*3)
c        ff :  MVMAXA*3            doubles   (was NPP*3)
c
c     where MVMAXA = max(MVMAX,NPP_max) and MVMAX is the parameter set
c     below.  Pxyz, RS, R0, R0C, IJK are unchanged.
c
c     ff only needs 2 columns; the third has always been unused, both
c     here and in the original code.  It is kept at 3 so that the
c     caller's existing "3" sizing expression needs no edit.
c
c     To disable batching entirely set MVMAX=1: nc becomes 1, MVEC
c     becomes NPP, and the routine reduces to the original behaviour.
c     ---------------------------------------------------------------
c
      Implicit none
c::passed
      Integer ncenters                  ! number of nuclei
      Integer NPP                       ! number of primitive products
      Integer Lp                        ! angular momentum sum La+Lb+MXD
      Integer Lp3                       ! number of angular momentum functions
      Double precision Axyz(3),Bxyz(3)  ! Cartesian coordinates of centers
      Double precision Cxyz(3,ncenters) ! Nuclear cartesian coordinates
      Double precision zan(ncenters)    ! Nuclear charges
      Double precision alpha(2,NPP)     ! basis function exponents copy
      Double precision Pxyz(3,NPP)      ! coordinates of product center
      Double precision RS(NPP)          ! scale factor
      Double precision PC(*)            ! (MVEC,3)        P to nucleus
      Double precision ff(*)            ! (MVEC,3)        recursion factors
      Double precision R(*)             ! (MVEC,0:Lp,Lp3) R000j / RIKJj
      Double precision R0(NPP,Lp3)      ! R integrals
      Double precision R0C(ncenters,NPP,Lp3) ! R integrals per center
      integer IJK(0:Lp,0:Lp,0:Lp)       ! index array
      Logical CENTER
c::local
      integer ic,icb,nc,nk,nv,j,k,m,n,i,off,ldj
      integer MVEC
      Double precision a,b,apb,f1,f2,PCx,PCy,PCz
      Double precision zk,acc
      Double precision Ax1,Ax2,Ax3,Bx1,Bx2,Bx3,Cx,Cy,Cz
      Double precision PI, PI4          ! pi, 4/pi
      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)
#include "hf1.fh"
c
c Define the auxiliary function integrals necessary to compute the nuclear
c attraction integrals (NAIs). These integrals are scaled by an appropriate
c factor, RS, defined as
c
c         / a + b \ 1/2
c   RS = | ------- |
c         \  PI/4 /
c
c The scale factor for the Hermite expansion coefficients is assumed to be
c
c         /   PI  \ 3/2     /   a b   __2 \
c   ES = | ------- |    EXP| - -----  AB   |
c         \ a + b /         \  a + b      /
c
c Therefore,
c
c            2 PI        /   a b   __2 \
c   ES RS = -------  EXP| - -----  AB   |
c            a + b       \  a + b      /
c
c******************************************************************************
c
c --- batch size and leading dimension -------------------------------------
      nc = MVMAX/NPP
      if (nc .lt. 1)        nc = 1
      if (nc .gt. ncenters) nc = ncenters
      MVEC = NPP*nc
c
c     stride between successive j planes of R is MVEC;
c     stride between successive n planes is MVEC*(Lp+1)
      ldj = MVEC*(Lp+1)
c
c --- zero the accumulator -------------------------------------------------
c     n outer / m inner: stride-1 on R0(m,n).  The original nest had
c     these inverted and did not vectorize.
      do 101 n = 1,Lp3
       do 100 m = 1,NPP
        R0(m,n) = 0.D0
  100  continue
  101 continue
c
c --- cache A and B coordinates -------------------------------------------
      Ax1 = Axyz(1)
      Ax2 = Axyz(2)
      Ax3 = Axyz(3)
      Bx1 = Bxyz(1)
      Bx2 = Bxyz(2)
      Bx3 = Bxyz(3)

      do 110 m = 1,NPP

c Define the center "P".

       a = alpha(1,m)
       b = alpha(2,m)
       apb = a + b
       f1 = a/apb
       f2 = b/apb
       Pxyz(1,m) = f1*Ax1 + f2*Bx1
       Pxyz(2,m) = f1*Ax2 + f2*Bx2
       Pxyz(3,m) = f1*Ax3 + f2*Bx3

c Define the scaling factor.

       RS(m) = sqrt(apb*PI4)


  110 continue
c
c==========================================================================
c --- loop over batches of centers ----------------------------------------
c==========================================================================
      do 200 icb = 1,ncenters,nc
       nk = min(nc,ncenters-icb+1)     ! centers in this batch
       nv = NPP*nk                     ! lanes actually in use
c
c ----- PC, ff and the incomplete gamma argument --------------------------
c       lane index i = (k-1)*NPP + m  for center k, primitive m
       do 130 k = 1,nk
        ic  = icb + k - 1
        off = (k-1)*NPP
        Cx  = Cxyz(1,ic)
        Cy  = Cxyz(2,ic)
        Cz  = Cxyz(3,ic)
        do 120 m = 1,NPP
         i   = off + m
         apb = alpha(1,m) + alpha(2,m)
         PCx = Pxyz(1,m) - Cx
         PCy = Pxyz(2,m) - Cy
         PCz = Pxyz(3,m) - Cz
         PC(i)          = PCx
         PC(MVEC+i)     = PCy
         PC(2*MVEC+i)   = PCz
         ff(i)          = RS(m)
         ff(MVEC+i)     = -2.D0*apb
c        apb*(PC.PC) is identical to -0.5*ff(m,2)*(PC.PC)
         R(i)           = apb*(PCx*PCx + PCy*PCy + PCz*PCz)
  120   continue
  130  continue
c
c ----- pad the unused tail ----------------------------------------------
c       igamma and hfmkr are called with length MVEC, not nv, so that
c       their internal declarations PQ(NAB,3) and R(NAB,0:Lp,Lp3) use
c       the same leading dimension MVEC as the manual indexing above.
c       Lanes nv+1..MVEC are set so they produce finite values that are
c       subsequently never accumulated.  T=0 gives F_j(0)=1/(2j+1),
c       and ff=0 then zeroes those lanes in loop 145.
       do 135 i = nv+1,MVEC
        PC(i)          = 0.D0
        PC(MVEC+i)     = 0.D0
        PC(2*MVEC+i)   = 0.D0
        ff(i)          = 0.D0
        ff(MVEC+i)     = 0.D0
        R(i)           = 0.D0
  135  continue

c Evaluate the incomplete gamma function.

       call igamma(R,MVEC,Lp)
c
c Define the initial auxiliary functions (i.e., R000j, j=0,Lp).
c       flat loop of length MVEC; vectorizes even when NPP=1.
c       The Lp=0 special case is unnecessary here: the j loop runs once.
       do 145 j = 0,Lp
        off = j*MVEC
        do 140 i = 1,MVEC
         R(off+i) = ff(i)*R(off+i)
         ff(i)    = ff(i)*ff(MVEC+i)
  140   continue
  145  continue
c
c ----- recursively build the remaining auxiliary functions (RIJKj) ------
       call hfmkr(R,IJK,PC,MVEC,Lp,Lp3)
c
c ----- transfer to R0 ---------------------------------------------------
c       R(i,0,n) is at flat offset (n-1)*ldj + i
       if (CENTER) then
c
        do 160 n = 1,Lp3
         off = (n-1)*ldj
         do 155 k = 1,nk
          ic = icb + k - 1
          zk = zan(ic)
          i  = off + (k-1)*NPP
          do 150 m = 1,NPP
           R0C(ic,m,n) = -zk*R(i+m)
           R0(m,n) = R0(m,n) + R0C(ic,m,n)
  150     continue
  155    continue
  160   continue
c
       else if (NPP .eq. 1) then
c
c       lane index equals the center index, so R(off+k) is stride-1 in k
c       and the center sum becomes a contiguous reduction of length nk
        do 175 n = 1,Lp3
         off = (n-1)*ldj
         acc = 0.D0
         do 170 k = 1,nk
          acc = acc + zan(icb+k-1)*R(off+k)
  170    continue
         R0(1,n) = R0(1,n) - acc
  175   continue
c
       else
c
        do 190 n = 1,Lp3
         off = (n-1)*ldj
         do 185 k = 1,nk
          zk = zan(icb+k-1)
          i  = off + (k-1)*NPP
          do 180 m = 1,NPP
           R0(m,n) = R0(m,n) - zk*R(i+m)
  180     continue
  185    continue
  190   continue
c
       end if
c
  200 continue
c
      end
