      subroutine int_1edk_ga(basis_in, g_in, type, oskel)
c $Id: int_1edk_ga.F,v 1.7 1999-11-03 00:35:33 bert Exp $
      implicit none
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "global.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "apiP.fh"
#include "bas.fh"
#include "sym.fh"
#include "basP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "bgj.fh"
c     
c     Compute the Douglas-Kroll integrals and place then in a given global array.
c
c     There are three forms of Douglas-Kroll that are or can be used:
c     1. Free-particle projection operators, similar to Foldy-Wouthuysen 
c        transformation.
c     2. External potential projection operators without pxVp terms 
c        (vn2 approximation).
c     3. External potential projection operators, i.e  complete Douglas-Kroll.
c
c     Method:
c     A. Determine uncontracted basis (usebasis) from input basis (basis_in)
c     B. Create p^2/2 matrix and diagonalize to get real to momentum space 
c        transformation.
c     C. Create diagonal matrices (store as vectors) of momentum operators in 
c        momentum space.
c     D. Compute all terms of E, V and p.Vp terms that are needed for form 1 
c        (see above).
c     E. Compute all combined terms that are needed for form 2 (see above) 
c        and add to DKg_a.
c     F. Compute the three W terms with pxVp required in form 3 (see above) 
c        and add them to the DKg_a array
c     G. Transform DKg_a to real space
c     H. Transform DKg_a from uncontracted basis to contracted basis
c     
c     arguments
c     
      integer basis_in          ! [input] basis sets for which integrals are computed
      integer g_in              ! [output] GA handle to array of final integrals
      integer type              ! [input] type of integral, 1=kinetic, 2=potential
c     
c     local variables
c     
      integer form_type         ! form of Douglas-Kroll (see first section header)
      integer nbf, nbf_in       ! # of basis functions for usebasis and basis_in
      integer usebasis          ! Uncontracted basis to be used for computation of integrals
      integer Tg_a              ! GA handle kinetic energy integral matrix T 
      integer Ug_a              ! GA handle transformation matrix r -> p space 
      integer Vg_a              ! GA handle potential energy integral matrix V
      integer Vdg_a             ! GA handle p.Vp and later pxVP integral matrix
      integer DKg_a             ! GA handle Douglas-Kroll integrals in uncontracted basis
      integer tempg_a           ! Temporary g_a matrix #1 handle
      integer tmpg_a            ! Temporary g_a matrix #2 handle
      integer Mixg_a            ! GA handle for transformation between contr and uncont
      integer p2m_a, p2Im_a     ! MA handle and index p^2 vector 
      integer Em_a, EIm_a       ! MA handle and index E vector
      integer Am_a, AIm_a       ! MA handle and index A vector
      integer Km_a, KIm_a       ! MA handle and index K vector
      integer tempm_a, tempIm_a ! Temporary m_a vector handle and index
      integer iloop             
      integer basis_array(2)
c
      double precision c,c2  ! Speed of light c (in au) plus c^2 
c
      logical oskel
c     
      c=137.0359895d0
      c2=c*c
c
      call ga_sync()
c
c     Create an uncontracted basis from basis_in
c
      call bas_create_uncont(basis_in,'D-K basis',usebasis)
      call int_normalize(usebasis)
c
c     Re-initialize integrals so that pVp integral memory requirement is set
c
      call int_terminate
      basis_array(1) = basis_in
      basis_array(2) = usebasis
      call int_init(bgj_get_rtdb_handle(),2,basis_array)
c
c     Define the dimensions of the arrays
c
      if (.not. bas_numbf(usebasis,nbf))
     $   call errquit('int_1edk_ga: bas_numbf failed',0)
c
c     Create DKg_a of nbf*nbf to store integrals of uncontracted basis
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA D-K',1,1,DKg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA D-K',0)
c
c     Create a temp g_a of nbf*nbf as temp/help matrix array
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 1',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 2',1,1,tmpg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0)
c
c     Create transformation matrix r-space to g-space, U
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Umatrix',1,1,Ug_a))
     $     call errquit('int_1edk_ga: ga_create failed Umatrix',0)
      call ga_zero(Ug_a)
c
c     Get kinetic energy matrix p^2/2 in uncontracted basis, diagonalize
c     Store transformation matrix and vector of eigenvalues
c
c     Construct, fill Tg_a with kinetic energy matrix integrals
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Tmatrix',1,1,Tg_a))
     $     call errquit('int_1edk_ga: ga_create failed Tmatrix',0)
      call ga_zero(Tg_a)
      call int_1e_oldga(usebasis,usebasis,Tg_a,'kinetic',oskel)
c
c     Diagonalize T in the orthogonal basis, using the overlap S as metric to get XU
c
c     Put eigenvalues in p2g_a and multiply by 2.0 to get p^2,
c     destroy T because we don't need that one anymore
c
      call ga_zero(tempg_a)
      call int_1e_oldga(usebasis,usebasis,tempg_a,'overlap',oskel)
      if (.not. ma_push_get(MT_DBL,nbf,'psquared',p2m_a,p2Im_a))
     $     call errquit('int_1edk_ga: ma_push_get failed psquared',0)
      call ga_diag(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
c     call ga_diag_std(Tg_a,Ug_a,dbl_mb(p2Im_a))
      if (.not. ga_destroy(Tg_a)) 
     $     call errquit('int_1edk_ga: ga_destroy Tg_a',0)
      call dscal(nbf,2.0d0,dbl_mb(p2Im_a),1)
c
c     Construct vector E from eigenvalues: E = Sqrt(p^2 c^2 + c^4)
c
      if (.not. ma_push_get(MT_DBL,nbf,'Evector',Em_a,EIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Evector',0)
      do 10 iloop = 1, nbf
         dbl_mb(EIm_a+iloop-1) = c*Sqrt(dbl_mb(p2Im_a+iloop-1)+c2)
   10 continue
c
c     Construct vector A from E: A = Sqrt((E+c^2)/2E) and
c               vector K from E: K = Sqrt(c/(E+c^2))
c
      if (.not. ma_push_get(MT_DBL,nbf,'Avector',Am_a,AIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Avector',0)
      if (.not. ma_push_get(MT_DBL,nbf,'Kvector',Km_a,KIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Kvector',0)
      do 20 iloop = 1, nbf
         dbl_mb(KIm_a+iloop-1) = c/(dbl_mb(EIm_a+iloop-1)+c2)  
         dbl_mb(AIm_a+iloop-1) = Sqrt((dbl_mb(EIm_a+iloop-1)+c2)/
     $                                (2.0d0*dbl_mb(EIm_a+iloop-1)))
   20 continue
c
c     Clear given global array and put in relativistic kinetic energy E on
c     diagonal (E matrix diagonal in momentum space)
c
c     We have to subtract mc^2 to get the non-relativistic type integrals
c
      call ga_zero(DKg_a)
      call ga_copy_todiag(DKg_a,dbl_mb(EIm_a))
      call ga_add_todiag(DKg_a,-c2)
c
c     We have our kinetic energy integrals
c 
      if (type .eq. 1) goto 100
c
c     If we just need the potential energy integrals we have to clear Dkg_a
c
      if (type .eq. 2) call ga_zero(DKg_a)
c
c     Get potential energy matrix V in uncontracted basis, add AVA term to global array:
c
c     <phi|AVA|phi> = AU^tVUA
c
c     First create U^tVU in Vg_a
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Vmatrix',1,1,Vg_a))
     $     call errquit('int_1edk_ga: ga_create failed Vmatrix',0)
      call ga_zero(Vg_a)
      call int_1e_oldga(usebasis,usebasis,Vg_a,'potential',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Vg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Vg_a)
c
c     Do A Vg_a A, remember A is diagonal stored as vector and store in DKg_a
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,Dkg_a)
c
c     Get double derivative potential energy matrix p.Vp in uncontracted basis, 
c     add ARVRA term to global array:
c
c     <phi|ARVRA|phi> = AKU^t(p.Vp)UKA
c
c     First create U^t(p.Vp)U in Vdg_a
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Vdmatrix',1,1,Vdg_a))
     $     call errquit('int_1edk_ga: ga_create failed Vdmatrix',0)
c
c     get p.Vp integrals and put them in Vdg_a
c
      call ga_zero(Vdg_a)
      call int_1e_oldga(usebasis,usebasis,Vdg_a,'pVp',oskel)
c
c     Do U^t Vdg_a U to transform to momentum space
c
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Vdg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Vdg_a)
c
c     Do A K Vdg_a K A, remember A and K are diagonal stored as vector and add to DKg_a
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,DKg_a)
c
c     At this point we have the basic integrals for the free-particle projection
c     operators. The remaining terms come from the external potential. Skip and 
c     construct final integrals if form_type=1 (see header).
c 
      if (form_type .eq. 1) goto 100
c
c     Scale U^tVU (in Vg_a) and U^tp.VpU (in Vdg_a) with 1/(Ep+Ep') for W terms
c
      call ga_scale_divEpp(Vg_a,dbl_mb(EIm_a))
      call ga_scale_divEpp(Vdg_a,dbl_mb(EIm_a))
c
c     Compute - W(1)E(p)W(1) term and add to global array:
c
c     - <phi|WEW|phi> = - UAK Vdg_a KAE(p)A          Vg_a  AU^t     
c
c                       - UA  Vg_a  AE(p)AK          Vdg_a KAU^t 
c
c                       + UA  Vg_a  Ap^2K^2E(p)A     Vg_a  AU^t
c
c                       + UAK Vdg_a KAE(p)K^-2p^-2AK Vdg_a KAU^t
c
c                       - UAK pxVp  KAE(p)K^-2p^-2AK pxVp  KAU^t
c
c     (Last term done at the end, reusing p.Vp space)
c
c     Use an intermediate array to generate combined scaling factors
c
      if (.not. ma_push_get(MT_DBL,nbf,'MA temp',tempm_a,tempIm_a))
     $     call errquit('D-K integrals: ma_push_get failed MA temp',0)
c
c     Term 1: - UAK Vdg_a KAE(p)A Vg_a AU^t
c     
      do 30 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(EIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   30 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     Term 2: - UA Vg_a AE(p)AK Vdg_a KAU^t
c     (Same intermediate array as Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rh(tmpg_a,dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     Term 3: + UA Vg_a Ap^2K^2E(p)A Vg_a AU^t
c
      do 33 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(p2Im_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(EIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   33 continue
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmpg_a,DKg_a)
c
c     Term 4: + UAK Vdg_a KAE(p)K^-2p^-2AK Vdg_a KAU^t
c
      do 34 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(EIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   34 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmpg_a,DKg_a)
c
c     Compute - 1/2 E(p)W(1)^2 term and add to global array:
c
c     - <phi|EW^2|phi> = - UE(p)AK Vdg_a KAA          Vg_a  AU^t
c                     
c                        - UE(p)A  Vg_a  AAK          Vdg_a KAU^t
c      
c                        + UE(p)A  Vg_a  Ap^2K^2A     Vg_a  AU^t
c
c                        + UE(p)AK Vdg_a KAK^-2p^-2AK Vdg_a KAU^t
c
c                        - UE(p)AK pxVp  KAK^-2p^-2AK pxVp  KAU^t
c
c     (Last term done at the end, reusing p.Vp space)
c
c     Term 1: - UE(p)AK Vdg_a KAA Vg_a AU^t
c
      do 40 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   40 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     Term 2: - UE(p)A Vg_a AAK Vdg_a KAU^t
c     (Intermediate array is the same as in Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     Term 3: + UE(p)A Vg_a Ap^2K^2A Vg_a AU^t
c
      do 41 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(p2Im_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   41 continue
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     Term 4: + UE(p)AK Vdg_a KAK^-2p^-2AK Vdg_a KAU^t
c
      do 42 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   42 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     Compute - 1/2 W(1)^2E(p) term and add to global array:
c
c     - <phi|EW^2|phi> = - UAK Vdg_a KAA          Vg_a  AE(p)U^t
c                      
c                        - UA  Vg_a  AAK          Vdg_a KAE(p)U^t
c      
c                        + UA  Vg_a  Ap^2K^2A     Vg_a  AE(p)U^t
c
c                        + UAK Vdg_a KAK^-2p^-2AK Vdg_a KAE(p)U^t
c
c                        - UAK pxVp  KAK^-2p^-2AK pxVp  KAE(p)U^t
c
c     (Last term done at the end, reusing p.Vp space)
c
c     Term 1: - UAK Vdg_a KAA Vg_a AE(p)U^t
c
      do 50 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   50 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     Term 2: - UA Vg_a AAK Vdg_a KAE(p)U^t
c     (Intermediate array is the same as in Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a) 
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     Term 3: + UA Vg_a Ap^2K^2A Vg_a AE(p)U^t
c
      do 51 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(p2Im_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
   51 continue
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     Term 4: + UAK Vdg_a KAK^-2p^-2AK Vdg_a KAE(p)U^t
c
      do 52 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   52 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     At this point we have the basic integrals for the external potential projection
c     operators without the pxVp type terms. These remaining pxVp terms come from 
c     are needed for the full Douglas-Kroll set of integrals. If we neglect these 
c     terms we get an incomplete Douglas-Kroll transformation a la Haberlen&Rosch 
c     (CPL, 199, 491, 1992) denoted vn2.
c     Skip and construct final integrals if form_type=2 (see header).
c 
      if (form_type .eq. 2) then
         if (.not. ma_pop_stack(tempm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack tempm_a',0)
         goto 100
      endif
c
c     Unfortunately for the full Douglas-Kroll we do need for each term 
c     integrals of the type (pxVp)D(pxVp). We don't need the Vdg_a array 
c     anymore so we will reuse that for to get the pxVp integrals and to 
c     construct the final three terms.
c
c     NOTE: If we neglect these terms we get an incomplete Douglas-Kroll
c           transformation a la Haberlen&Rosch (CPL, 199, 491, 1992)
c
c
c     Get pxVp (in Vdg_a) and transform to U^tVU and scale with 1/(Ep+Ep') for W terms
c
      call ga_zero(Vdg_a)
      call int_1e_oldga(usebasis,usebasis,Vdg_a,'pxVp',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Vdg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Vdg_a)
      call ga_scale_divEpp(Vdg_a,dbl_mb(EIm_a))
c
c     W(1)E(p)W(1): - UAK     pxVp  KAE(p)K^-2p^-2AK pxVp  KAU^t
c
      do 60 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(EIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   60 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     E(p)W(1)^2:   - UE(p)AK pxVp  KAK^-2p^-2AK     pxVp  KAU^t
c
      do 61 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   61 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     W(1)^2E(p):   - UAK     pxVp  KAK^-2p^-2AK     pxVp  KAE(p)U^t
c
      do 62 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   62 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     At this point we have the integrals we want in momentum space. We have to 
c     backtransform to real space and, as we work in a primitive basis, contract
c     them back into the contracted basis.
c
      if (.not. ma_pop_stack(tempm_a))
     $    call errquit('int_1edk_ga: ma_pop_stack tempm_a',0)
  100 continue
c
c     Get overlap and multiply with Ug_a to get SXU
c
      call ga_zero(tempg_a)
      call ga_zero(tmpg_a)
      call int_1e_oldga(usebasis, usebasis, tempg_a, 'overlap', .false.)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,
     $              tmpg_a)
c
c     Transform final integrals back to real space integrals U V U^t
c
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,DKg_a,0.0d0,
     $              tempg_a)
      call ga_dgemm('n','t',nbf,nbf,nbf,1.0d0,tempg_a,tmpg_a,0.0d0,
     $              DKg_a)
c
c     Free MA's and GA's that are not needed anymore at this time
c
      if (.not. ma_pop_stack(Km_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Km_a',0)
      if (.not. ma_pop_stack(Am_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Am_a',0)
      if (.not. ma_pop_stack(Em_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Em_a',0)
      if (.not. ma_pop_stack(p2m_a))
     $     call errquit('int_1edk_ga: ma_pop_stack p2m_a',0)
      if (type .gt. 1) then
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0)
      endif
      if (.not. ga_destroy(Ug_a))
     $     call errquit('int_1edk_ga: ga_destroy Ug_a',0)
      if (.not. ga_destroy(tmpg_a))
     $     call errquit('int_1edk_ga: ga_destroy tmpg_a',0)
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0)
c
c     We now have our integrals in a primitive basis in DKg_a. If the original
c     basis is contracted we have to transform to this basis. Get transformation
c     matrix, combine with DKg_a and put final integrals in g_in
c
c     Approach: 
c
c     For j = 1, ncontr
c         g_in(i,j) = sum(k,l) S^-1(i,k) Scontr-uncon(k,l) DKg_in(l,j)
c     
c     S_-1 is inverse overlap in the uncontracted basis
c     Scontr_uncon is the mixed basis overlap
c
      if (.not. bas_numbf(basis_in,nbf_in))
     $   call errquit('int_1edk_ga: bas_numbf failed',0)
      if (.not. ga_create(MT_DBL,nbf,nbf_in,'Mixed overlap',1,1,
     $     Mixg_a)) call errquit('int_1edk_ga: create mixed overlap',0)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0)
      if (.not. ga_create(MT_DBL,nbf,nbf_in,'GA temp',1,1,tmpg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA tmp',0)
c
c     Create mixed overlap between basis_in and usebasis
c
      call ga_zero(Mixg_a)
      call int_1e_oldga(usebasis,basis_in,Mixg_a,'overlap',.false.)
c
c     Create overlap matrix usebasis, invert to S^-1 and multiply with Mixg_a
c
      call ga_zero(tempg_a)
      call int_1e_oldga(usebasis,usebasis,tempg_a,'overlap',.false.)
      call ga_matpow(tempg_a, -1.0d0, 0.0d0)
      call ga_dgemm('n','n',nbf,nbf_in,nbf,1.0d0,tempg_a,Mixg_a,
     $              0.0d0,tmpg_a)
      call ga_copy(tmpg_a,Mixg_a)
c
c     Create integrals in contracted basis g_in via Mixg_a x DKg_a Mixg_a^T
c
      if (.not. ga_destroy(tmpg_a))
     $     call errquit('int_1edk_ga: ga_destroy tmpg_a',0)
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0)
      if (.not. ga_create(MT_DBL,nbf_in,nbf,'GA temp',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0)
      call ga_zero(g_in)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf_in,nbf,nbf,1.0d0,Mixg_a,
     $              DKg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf_in,nbf_in,nbf,1.0d0,tempg_a,Mixg_a,
     $              0.0d0,g_in)
c
c     Destroy transformation matrices Mixg_a, tempg_a and DKg_a
c
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0)
      if (.not. ga_destroy(Mixg_a))
     $     call errquit('int_1edk_ga: ga_destroy Mixg_a',0)
      if (.not. ga_destroy(DKg_a))
     $     call errquit('int_1edk_ga: ga_destroy DKg_a',0)
c
c     We don't need the uncontracted basis anymore, destroy it also
c
      if (.not. bas_destroy(usebasis))
     $   call errquit('int_1edk_ga: bas_destroy usebasis failed',0)
c
      return
c
      end
