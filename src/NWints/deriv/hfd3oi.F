*........................................#define DEBUG_D3O
      subroutine hfd3ois(
     &    Axyz,Aprims,Acoef,NPA,NCA,La,
     &    Bxyz,Bprims,Bcoef,NPB,NCB,Lb,
     &    Cxyz,Cprims,Ccoef,NPC,NCC,Lc,
     &    d3OI,Nintin,
     &    DryRun,W0,maxW0)
* $Id$
      implicit none
#include "stdio.fh"
*
*  This routine computes the 3 center overlap derivative integrals
*       for the contracted gaussian function information passed in
*  There are 3 component vectors d/dA d/dB and d/dC in the result 
*  matrix.
*

*:: passed
      integer La  ! [input] angular momentum of A function
      integer NPA ! [input] number of primitive gausians on A
      integer NCA ! [input] number of general contractions on B
      double precision Axyz(3)        ! [input] coordinates of center A
      double precision Aprims(NPA)    ! [input] primitive exponents on A
      double precision Acoef(NPA,NCA) ! [input] contraction coefficients on A 

      integer Lb  ! [input] angular momentum of B function
      integer NPB ! [input] number of primitive gausians on B
      integer NCB ! [input] number of general contractions on B
      double precision Bxyz(3) ! [input] coordinates of center B
      double precision Bprims(NPB)    ! [input] primitive exponents on B
      double precision Bcoef(NPB,NCB) ! [input] contraction coefficients on B 

      integer Lc  ! [input] angular momentum of C function
      integer NPC ! [input] number of primitive gausians on C
      integer NCC ! [input] number of general contractions on C
      double precision Cxyz(3) ! [input] coordinates of center C
      double precision Cprims(NPC)    ! [input] primitive exponents on C
      double precision Ccoef(NPC,NCC) ! [input] contraction coefficients on C 
      
      double precision d3OI(*) ! [output] 3-center derivative overlap integrals
      integer Nintin           ! [input] leading dimension of d3OI
      
      logical DryRun  ! [input]
      integer maxW0   ! if (!DryRun) [input]  size of scratch array
*                     ! if (DryRun)  [output] estimated size 
      double precision W0(maxW0) ! [scratch] array.

*:: local
      integer i_top   ! [W0 pointer] next pointer value in W0
      integer i_used  ! [W0 pointer] current memory used on stack of W0
      integer i_save  ! [W0 pointer] pointer set to mark reuse point on W0 stack
      integer i_alpha ! [W0 pointer] exponents and Kabc 
      integer i_G     ! [W0 pointer] Coordinates of global center 
*................................... This is exponent dependent
      integer i_E     ! [W0 pointer] Hermite expansion coefficients
      integer i_Gt    ! [W0 pointer] Coordinates of global center offset by recusrion center
      integer i_ABC2i ! [W0 pointer] 1/2(a+b+c) used for recursion
      integer i_exabc ! [W0 pointer] exponents for performance
      integer MaxMem  ! maximum memory determined with DryRun=True call
      integer Nabc    ! Number of exponent "triplets" after screening
      integer Lap, Lbp, Lcp ! incremented L values (for derivs)
      integer Labc    ! Total derivative angular momentum
      integer sizeE   ! Size of E hermite expansion coefficients
      integer nint    ! block of integral size
*:: exceptions
* right now general contractions are out      
      if ((NCA.gt.1).or.(NCB.gt.1).or.(NCC.gt.1)) then
        write(luout,*)' hfd3oi: General contractions not supported'
        stop ' hfd3oi: fatal error '
      endif
#if defined(DEBUG_D3O)
      write(6,*)' hfd3ois debug'
      write(6,*)' La, Lb, Lc ',La, Lb, Lc
      write(6,*)' Nintin ',Nintin
      write(6,*)' maxW0  ',maxW0
#endif
      call dfill(maxW0,0.0d00,W0,1)

* allocate memory on W0 stack
      i_alpha = 1                               ! needed always
      i_top   = i_alpha + (NPA*NPB*NPC)*4
      i_used  = i_top - 1
      i_save  = i_top ! all pointers needed

      if (i_used.gt.maxW0) then
        write(luout,*)' hfd3oi: insufficient scratch space '
        write(luout,*)' hfd3oi: scratch required  :',i_used
        write(luout,*)' hfd3oi: scratch allocated :',maxW0
        stop ' hfd3oi: fatal error '
      endif
      
      MaxMem = i_used  ! take care of compiler warnings

      if (DryRun) then
        MaxMem = i_used
        Nabc = NPA*NPB*NPC
      else
** Define the prefactor of the charge distributions 
** Determine the primitive exponents used after screening
        call hf1set3(
     &      Axyz,Aprims,Acoef,NPA,
     &      Bxyz,Bprims,Bcoef,NPB,
     &      Cxyz,Cprims,Ccoef,NPC,
     &      W0(i_alpha),Nabc)
      endif

** Define the global center of the charge distribution
*** shifted angular momentum      
      Lap  = La + 1   ! incremented for integral derivatives
      Lbp  = Lb + 1   ! incremented for integral derivatives
      Lcp  = Lc + 1   ! incremented for integral derivatives
      Labc = Lap + Lbp + Lcp ! maximum total angular momentum value 

*** determine size of expansion coefficient array E
      sizeE = Nabc*3         ! exponent range (1:Nabc) and xyz range (1:3)
      sizeE = sizeE*(Labc+1) ! total angular momentum range (0:Labc)
      sizeE = sizeE*(Lap+1)  ! angular momentum range for A (0:Lap)
      sizeE = sizeE*(Lbp+1)  ! angular momentum range for B (0:Lbp)      
      sizeE = sizeE*(Lcp+1)  ! angular momentum range for C (0:Lcp)
***: allocate off W0 stack
      i_E    = i_top          ! i_E needed always
      i_G    = i_E + sizeE    ! i_G needed only to compute E coeffs
      i_top  = i_G + Nabc*3
      i_save = i_G     ! i_alpha and i_E are needed later
      i_used = i_top - 1

      if (i_used.gt.maxW0) then
        write(luout,*)' hfd3oi: insufficient scratch space '
        write(luout,*)' hfd3oi: scratch required  :',i_used
        write(luout,*)' hfd3oi: scratch allocated :',maxW0
        stop ' hfd3oi: fatal error '
      endif

      if (DryRun) then
        MaxMem = max( MaxMem, i_used)
      else
        call hfctr3(Axyz,Bxyz,Cxyz,W0(i_alpha),W0(i_G),Nabc)
      endif

** Determine the Hermite expansion coefficients (E)      
*** compute memory pointer for W0 stack

      i_Gt     = i_top            ! i_Gt needed only to compute E
      i_ABC2I  = i_Gt    + Nabc*3 ! i_ABC2I needed only to compute E
      i_top    = i_ABC2I + Nabc*3 
      i_used   = i_top - 1 

      if (i_used.gt.maxW0) then
        write(luout,*)' hfd3oi: insufficient scratch space '
        write(luout,*)' hfd3oi: scratch required  :',i_used
        write(luout,*)' hfd3oi: scratch allocated :',maxW0
        stop ' hfd3oi: fatal error '
      endif

      if (DryRun) then
        MaxMem = max( MaxMem, i_used)
      else
#if defined(DEBUG_D3O)
        write(6,*)' i_alpha = ', i_alpha
        write(6,*)' i_G     = ', i_G
        write(6,*)' i_Gt    = ', i_Gt
        write(6,*)' i_ABC2I = ', i_ABC2I
        write(6,*)' i_E     = ', i_E
#endif
        call hf1mke3(Axyz,Bxyz,Cxyz,
     &      W0(i_alpha),
     &      W0(i_G),
     &      W0(i_Gt),
     &      W0(i_ABC2I),
     &      W0(i_E),
     &      Nabc,Lap,Lbp,Lcp)

      endif
c
** set up pointers to copy exponents to transformed array (for performance)
*compute-F-Separately:      Fsize = Nabc*3
*compute-F-Separately:      Fsize = Fsize*(La+1)
*compute-F-Separately:      Fsize = Fsize*(Lb+1)
*compute-F-Separately:      Fsize = Fsize*(Lc+1)
*compute-F-Separately:      Fsize = Fsize*2
*compute-F-Separately:      i_F     = isave
*compute-F-Separately:      i_exabc = i_F + Fsize
      i_exabc = i_save
      i_top   = i_exabc + Nabc*3
      i_used  = i_top - 1
      if (i_used.gt.maxW0) then
        write(luout,*)' hfd3oi: insufficient scratch space '
        write(luout,*)' hfd3oi: scratch required  :',i_used
        write(luout,*)' hfd3oi: scratch allocated :',maxW0
        stop ' hfd3oi: fatal error '
      endif

c Return the maximum amount of scratch space required by DryRun.

      if( DryRun )then
        MaxMem = max( MaxMem, i_used)
        maxW0 = MaxMem
        return
      end if

** copy exponents to performance array
      call hfdabc_exc(W0(i_alpha),W0(i_exabc),Nabc)
*compute-F-Separately:** compute F coefficients
*compute-F-Separately:      call hfmkF3(Nabc,La,Lb,Lc,Labc,W0(i_exabc),W0(i_F),W0(i_E))
** Compute the 3-center overlap derivative integrals

      nint =      (La+1)*(La+2)/2
      nint = nint*(Lb+1)*(Lb+2)/2
      nint = nint*(Lc+1)*(Lc+2)/2
      call hfdabc(W0(i_E),W0(i_exabc),d3oi,nint,Nabc,La,Lb,Lc,Labc)

      end
      subroutine hfdabc(E,ex,d3oi,nint,Nabc,La,Lb,Lc,Labc)
      implicit none
c
      integer La
      integer Lb
      integer Lc
      integer Labc
      integer Nabc
      integer nint
      
      double precision d3oi(nint,9)   ! 9 = d/dA,d/dB,d/dC

      double precision E(Nabc,3,0:Labc,0:(La+1),0:(Lb+1),0:(Lc+1))

      double precision ex(Nabc,3)  ! 1=A, 2=B, 3=C
c
      integer Nxyz(3)
      integer La2, Lb2, Lc2
      integer n
      integer m
      integer ica, icb, icc
      integer ia, ja, ka
      integer ib, jb, kb
      integer ic, jc, kc
      double precision sum1, sum2, sum3
      double precision a, b
      double precision E1, E2, E3
      double precision F1, F2, F3
c
      La2 = (La+1)*(La+2)/2
      Lb2 = (Lb+1)*(Lb+2)/2
      Lc2 = (Lc+1)*(Lc+2)/2
c
      n = 0
c
* loop over A center shell components
      do ica = 1,La2

        call getNxyz(La,ica,Nxyz)
        Ia = Nxyz(1)
        Ja = Nxyz(2)
        Ka = Nxyz(3)

        do icb = 1,Lb2

          call getNxyz(Lb,icb,Nxyz)
          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)

          do icc = 1,Lc2

            call getNxyz(Lc,icc,Nxyz)
            Ic = Nxyz(1)
            Jc = Nxyz(2)
            Kc = Nxyz(3)

            n = n + 1

            sum1 = 0.0d00
            sum2 = 0.0d00
            sum3 = 0.0d00
            do m = 1, Nabc
              E1 = E(m,1,0,Ia,Ib,Ic)
              E2 = E(m,2,0,Ja,Jb,Jc)
              E3 = E(m,3,0,Ka,Kb,Kc)
              a = ex(m,1)
              F1 = 2.0d00*a*E(m,1,0,(Ia+1),Ib,Ic) 
              if (Ia.gt.0) F1 = F1 -  dble(Ia)*E(m,1,0,(Ia-1),Ib,Ic)
              F2 = 2.0d00*a*E(m,2,0,(Ja+1),Jb,Jc) 
              if (Ja.gt.0) F2 = F2 -  dble(Ja)*E(m,2,0,(Ja-1),Jb,Jc)
              F3 = 2.0d00*a*E(m,3,0,(Ka+1),Kb,Kc) 
              if (Ka.gt.0) F3 = F3 -  dble(Ka)*E(m,3,0,(Ka-1),Kb,Kc)
              sum1 = sum1 + F1*E2*E3
              sum2 = sum2 + E1*F2*E3
              sum3 = sum3 + E1*E2*F3
            enddo
            d3oi(n,1) = sum1
            d3oi(n,2) = sum2
            d3oi(n,3) = sum3

            sum1 = 0.0d00
            sum2 = 0.0d00
            sum3 = 0.0d00
            do m = 1, Nabc
              E1 = E(m,1,0,Ia,Ib,Ic)
              E2 = E(m,2,0,Ja,Jb,Jc)
              E3 = E(m,3,0,Ka,Kb,Kc)
              b = ex(m,2)
              F1 = 2.0d00*b*E(m,1,0,Ia,(Ib+1),Ic) 
              if (Ib.gt.0) F1 = F1 - dble(Ib)*E(m,1,0,Ia,(Ib-1),Ic)
              F2 = 2.0d00*b*E(m,2,0,Ja,(Jb+1),Jc) 
              if (Jb.gt.0) F2 = F2 - dble(Jb)*E(m,2,0,Ja,(Jb-1),Jc)
              F3 = 2.0d00*b*E(m,3,0,Ka,(Kb+1),Kc) 
              if (Kb.gt.0) F3 = F3 - dble(Kb)*E(m,3,0,Ka,(Kb-1),Kc)
              sum1 = sum1 + F1*E2*E3
              sum2 = sum2 + E1*F2*E3
              sum3 = sum3 + E1*E2*F3
            enddo
            d3oi(n,4) = sum1
            d3oi(n,5) = sum2
            d3oi(n,6) = sum3

          enddo
        enddo
      enddo
c
* use translational invariance to compute d/dc = -d/da - d/db
c
*(1) d/dc = d/da
      call dcopy((3*nint),d3oi(1,1),1,d3oi(1,7),1)
*(2) d/dc = -d/da
      call dscal((3*nint),-1.0d00,d3oi(1,7),1)
*(3) d/dc = -d/da - d/db                  
      call daxpy((3*nint),-1.0d00,d3oi(1,4),1,d3oi(1,7),1)
      end
      subroutine hfdabc_exc(alpha,exabc,Nabc)
      implicit none
      integer Nabc
      double precision alpha(4,Nabc)
      double precision exabc(Nabc,3)
c
      integer i

      do i = 1, Nabc
        exabc(i,1) = alpha(1,i)
        exabc(i,2) = alpha(2,i)
        exabc(i,3) = alpha(3,i)
      enddo
      end
      subroutine hfmkF3(Nabc,La,Lb,Lc,Labc,exabc,F,E)
      implicit none
c
      integer Nabc
      integer Labc
      integer La
      integer Lb
      integer Lc
      double precision exabc(Nabc,3)
      double precision E(Nabc,3,0:Labc,0:(La+1),0:(Lb+1),0:(Lc+1))
      double precision F(Nabc,3,0:La,0:Lb,0:Lc,2)
c index description for F
c 1  [m]      Nabc
c 1  [xyz]    3
c 0  [ia]     La
c 0  [ib]     Lb
c 0  [ic]     Lc
c 1  [center] 2  (only centers 1 and 2 are computed center 3 by translational invariance
c
      integer ia, ib, ic
      integer m, ixyz
c
c for center A 
c
** ia = 0
      do ib = 0, Lb
        do ic = 0, Lc
          do ixyz = 1,3
            do m = 1,Nabc
              F(m,ixyz,0,ib,ic,1) =
     &            2.0d00*exabc(m,1)*E(m,ixyz,0,1,ib,ic)
            enddo
          enddo
        enddo
      enddo
** ia = 1 -> La
      do ia = 1, La
        do ib = 0, Lb
          do ic = 0, Lc
            do ixyz = 1,3
              do m = 1,Nabc
              F(m,ixyz,ia,ib,ic,1) =
     &              2.0d00*exabc(m,1)*E(m,ixyz,0,(ia+1),ib,ic) -
     &              dble(ia)*E(m,ixyz,0,(ia-1),ib,ic)
              enddo
            enddo
          enddo
        enddo
      enddo

c
c for center B 
c
** ib = 0
      do ia = 0, La
        do ic = 0, Lc
          do ixyz = 1,3
            do m = 1,Nabc
              F(m,ixyz,ia,0,ic,2) =
     &            2.0d00*exabc(m,2)*E(m,ixyz,0,ia,1,ic)
            enddo
          enddo
        enddo
      enddo
** ib = 1 -> Lb
      do ia = 0, La
        do ib = 1, Lb
          do ic = 0, Lc
            do ixyz = 1,3
              do m = 1,Nabc
                F(m,ixyz,ia,ib,ic,2) =
     &              2.0d00*exabc(m,2)*E(m,ixyz,0,ia,(ib+1),ic) -
     &              dble(ib)*E(m,ixyz,0,ia,(ib-1),ic)
              enddo
            enddo
          enddo
        enddo
      enddo
c
c for center C
c
*. . . .  These "F" values should not be needed since d/dc is computed
*         by translational invariance
*not-needed:** ic = 0
*not-needed:      do ia = 0, Lb
*not-needed:        do ib = 0, Lc
*not-needed:          do ixyz = 1,3
*not-needed:            do m = 1,Nabc
*not-needed:              F(m,ixyz,ia,ib,0,3) =
*not-needed:     &            2.0d00*exabc(m,3)*E(m,ixyz,0,ia,ib,1)
*not-needed:            enddo
*not-needed:          enddo
*not-needed:        enddo
*not-needed:      enddo
*not-needed:** ic = 1 -> Lc
*not-needed:      do ia = 0, La
*not-needed:        do ib = 0, Lb
*not-needed:          do ic = 1, Lc
*not-needed:            do ixyz = 1,3
*not-needed:              do m = 1,Nabc
*not-needed:              F(m,ixyz,ia,ib,ic,3) =
*not-needed:     &              2.0d00*exabc(m,3)*E(m,ixyz,0,ia,ib,(ic+1)) -
*not-needed:     &              dble(ic)*E(m,ixyz,0,ia,ib,(ic-1))
*not-needed:              enddo
*not-needed:            enddo
*not-needed:          enddo
*not-needed:        enddo
*not-needed:      enddo
c
      end
