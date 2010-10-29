      subroutine hf3pot(
     &    Axyz,Aprims,Acoef,NPA,NCA,La,
     &    Bxyz,Bprims,Bcoef,NPB,NCB,Lb,
     &    Cxyz,Cprims,Ccoef,NPC,NCC,Lc,
     &    b3nai,b3nai_sz,nint,
     &    DryRun,Scr,lscr)
c
c $Id$
c
      implicit none
c:: includes
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c:: passed 
      integer La  ! [input] angular momentum for center a
      integer Lb  ! [input] angular momentum for center b
      integer Lc  ! [input] angular momentum for center c
      integer NPA ! [input] number of primitives on center a
      integer NPB ! [input] number of primitives on center b
      integer NPC ! [input] number of primitives on center c
      integer NCA ! [input] number of general contractions on center a
      integer NCB ! [input] number of general contractions on center b
      integer NCC ! [input] number of general contractions on center c
      double precision Axyz(3) ! [input] coordinates for center a
      double precision Bxyz(3) ! [input] coordinates for center b
      double precision Cxyz(3) ! [input] coordinates for center c
      double precision Aprims(NPA) ! [input] primitive exponents on center a
      double precision Bprims(NPB) ! [input] primitive exponents on center b
      double precision Cprims(NPC) ! [input] primitive exponents on center c
      double precision Acoef(NPA,NCA) ! [input] primitive coefficients on center a
      double precision Bcoef(NPB,NCB) ! [input] primitive coefficients on center b
      double precision Ccoef(NPC,NCC) ! [input] primitive coefficients on center c
      integer b3nai_sz                 ! [input] input size of 3 center nai buffer
      double precision b3nai(b3nai_sz) ! [output] 3 center nai buffer
      integer nint                     ! [output] actual number of integrals computed
      logical DryRun                   ! [input] run routine only to test memory size needed
      integer lscr                     ! [input] size of scratch array
      double precision scr(lscr)       ! [scratch] scratch array for computations
c::local 
      logical gencon ! general contraction flag
      Integer Nabc   ! product of primitives
      Integer maxmem ! memory returned in lscr for dryrun
      Integer Lg     ! total angular momentum
      Integer Lg3    ! l3 of total angular momentum
c:: scratch pointers 
      Integer i_top   ! top of scratch
      Integer i_save  ! point at which scratch can be reused
      Integer i_alpha ! scr ptr for exponents and Hermite expansion coef prefactor
      Integer i_E     ! Hermite expansion coefficients
      Integer i_G1    ! center G coordinates (first use;through making of E)
      Integer i_G2    ! center G coordinates (second use; making of R)
      Integer i_Gt    ! Temp array for G coorinates
      Integer i_ABC2I ! Temp array for prefactors of E
      Integer i_R0    ! the auxillary hermite integrals required for sum
      Integer i_ijk   ! pointer for auxillary hermite integrals
      Integer i_RS    ! prefactors for auxillary hermite integrals
      Integer i_GC    ! center for auxillary hermite integral recursion
      Integer i_r     ! the auxillary hermite integrals required for recursion
      Integer i_ff    ! prefactor array for forming recusion set of R's 

#if defined(INTDEBUG)
      write(luout,*)' hf3pot debug'
      if (.not.dryrun) then
        write(luout,*)'b3nai_sz',b3nai_sz
        write(luout,*)'dryrun',dryrun
        write(luout,*)'lscr',lscr
        call hf_print_set(1)
        call hf_print('hf3pot a',Axyz,Aprims,Acoef,NPA,NCA,La)
        call hf_print('hf3pot b',Bxyz,Bprims,Bcoef,NPB,NCB,Lb)
        call hf_print('hf3pot c',Cxyz,Cprims,Ccoef,NPC,NCC,Lc)
        call hf_print_set(0)
      endif
#endif
c
c Compute 3-center nuclear attraction integrals for ECPs
c
      gencon = nca.gt.1.or.ncb.gt.1.or.ncc.gt.1
      if (gencon) call errquit
     &    (' hf3pot not ready for general contractions ',911,
     &       INT_ERR)
c
      i_alpha = 1
      i_top   = i_alpha + (NPA*NPB*NPC)*4 - 1

#if defined(INTDEBUG)
      write(luout,*)'i_alpha :',i_alpha
      write(luout,*)'i_top   :',i_top
      write(luout,*)npa,npb,npc,(NPA*NPB*NPC),(NPA*NPB*NPC)*4
#endif
      if (i_top.gt.lscr) then
        write(luout,*)' hf3pot: insufficient scratch space '
        write(luout,*)'       : needed    :',i_top
        write(luout,*)'       : allocated :',lscr
        call errquit('hf3pot: fatal error ',1, INT_ERR)
      endif
c
c.. determine actual product pairs to be kept. 
c   The following line is to take care of compiler warnings.
c
      maxmem = i_top
      if (DryRun) then
        maxmem = i_top
        Nabc = NPA*NPB*NPC
      else
        call hf1set3(
     &      Axyz,Aprims,Acoef,NPA,
     &      Bxyz,Bprims,Bcoef,NPB,
     &      Cxyz,Cprims,Ccoef,NPC,
     &      scr(i_alpha),Nabc)
      endif
#if defined(INTDEBUG)
      write(luout,*)
      write(luout,*)' nabc = ',nabc
      write(luout,*)
#endif

c Define the center of the charge distribution.
c A+B -> P; P + C -> G

      Lg = La+Lb+Lc ! total angular momentum

c
      i_E    = i_top + 1
      i_G1   = i_E   + Nabc*3*(Lg+1)*(La+1)*(Lb+1)*(Lc+1)
      i_save = i_G1  ! i_alpha and I_E saved
      i_top  = i_G1  + Nabc*3 - 1

#if defined(INTDEBUG)
      write(luout,*)
      write(luout,*)'i_alpha :',i_alpha
      write(luout,*)'i_E     :',i_E
      write(luout,*)'i_G1    :',i_G1
      write(luout,*)'i_save  :',i_save
      write(luout,*)'i_top   :',i_top
      write(luout,*)
#endif

      if (i_top.gt.lscr) then
        write(luout,*)' hf3pot: insufficient scratch space '
        write(luout,*)'       : needed    :',i_top
        write(luout,*)'       : allocated :',lscr

        write(luout,*)'i_alpha :',i_alpha
        write(luout,*)'i_E     :',i_E
        write(luout,*)'i_G1    :',i_G1
        call errquit('hf3pot: fatal error ',2, INT_ERR)
      endif
c
      if (DryRun) then
        maxmem = max(maxmem, i_top)
      else
        call hfctr3(Axyz,Bxyz,Cxyz,scr(i_alpha),scr(i_G1),Nabc)
      endif
c
c Define the Hermite linear expansion coefficients
c
      i_Gt    = i_top   + 1
      i_ABC2I = i_Gt    + Nabc*3
      i_top   = i_ABC2I + Nabc*3 - 1

#if defined(INTDEBUG)
      write(luout,*)
      write(luout,*)'i_save  :',i_save
      write(luout,*)'i_alpha :',i_alpha
      write(luout,*)'i_E     :',i_E
      write(luout,*)'i_G1    :',i_G1
      write(luout,*)'i_Gt    :',i_Gt
      write(luout,*)'i_ABC2I :',i_ABC2I
      write(luout,*)'i_top   :',i_top
      write(luout,*)
#endif 

      if (i_top.gt.lscr) then
        write(luout,*)' hf3pot: insufficient scratch space '
        write(luout,*)'       : needed    :',i_top
        write(luout,*)'       : allocated :',lscr

        write(luout,*)'i_alpha :',i_alpha
        write(luout,*)'i_E     :',i_E
        write(luout,*)'i_G1    :',i_G1
        write(luout,*)'i_Gt    :',i_Gt
        write(luout,*)'i_ABC2I :',i_ABC2I
        call errquit('hf3pot: fatal error ',3, INT_ERR)
      endif
c
      if (DryRun) then
        maxmem = max(maxmem, i_top)
      else
        call hf1mke3(Axyz, Bxyz, Cxyz,
     &      scr(i_alpha), scr(i_G1), scr(i_Gt), scr(i_ABC2I),
     &      scr(i_E), Nabc, La, Lb, Lc)
      endif
c
      Lg3   = ((Lg+1)*(Lg+2)*(Lg+3))/6

#if defined(INTDEBUG)
      write(luout,*)
      write(luout,*)'lg,lg3',lg,lg3
      write(luout,*)
#endif

      i_R0  = i_save 
      i_IJK = i_R0   + Nabc*Lg3 
      i_RS  = i_IJK  + (Lg+1)*(Lg+1)*(Lg+1)
      i_GC  = i_RS   + Nabc    
      i_G2  = i_GC   + Nabc*3  
      i_ff  = i_G2   + Nabc*3  
      i_R   = i_ff   + 2*Nabc  
      i_top = i_R    + Nabc*(Lg+1)*Lg3 - 1

#if defined(INTDEBUG)
      write(luout,*)
      write(luout,*)'i_save  :',i_save
      write(luout,*)'i_alpha :',i_alpha
      write(luout,*)'i_E     :',i_E
      write(luout,*)'i_R0    :',i_R0
      write(luout,*)'i_IJK   :',i_IJK
      write(luout,*)'i_RS    :',i_RS
      write(luout,*)'i_GC    :',i_GC
      write(luout,*)'i_G2    :',i_G2
      write(luout,*)'i_ff    :',i_ff
      write(luout,*)'i_R     :',i_R
      write(luout,*)'i_top   :',i_top
      write(luout,*)
#endif

      if (i_top.gt.lscr) then
        write(luout,*)' hf3pot: insufficient scratch space '
        write(luout,*)'       : needed    :',i_top
        write(luout,*)'       : allocated :',lscr

        write(luout,*)'i_alpha :',i_alpha
        write(luout,*)'i_E     :',i_E
        write(luout,*)'i_R0    :',i_R0
        write(luout,*)'i_IJK   :',i_IJK
        write(luout,*)'i_RS    :',i_RS
        write(luout,*)'i_GC    :',i_GC
        write(luout,*)'i_G2    :',i_G2
        write(luout,*)'i_ff    :',i_ff
        write(luout,*)'i_R     :',i_R
        call errquit('hf3pot: fatal error ',3, INT_ERR)
      endif
c
      if (DryRun) then
        maxmem = max(maxmem, i_top)
      else
        call hf3mkr(Axyz, Bxyz, Cxyz, scr(i_alpha), scr(i_G2),
     &      scr(i_RS), scr(i_GC), scr(i_ff), scr(i_R), scr(i_R0),
     &      scr(i_IJK), Nabc, Lg, Lg3)
      endif
c
      if (DryRun) then
        lscr = maxmem
        return
      endif
c
      call hf3PEabc(b3nai, scr(i_E), scr(i_R0), scr(i_IJK),
     &    Nabc, La, Lb, Lc, Lg, Lg3, nint, b3nai_sz)
c
      end
