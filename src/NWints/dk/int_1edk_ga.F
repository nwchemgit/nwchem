c $Id: int_1edk_ga.F,v 1.23 1999-12-30 23:06:35 bert Exp $
c:tex-\subsection{int\_1edk\_ga}
c:tex-This routine computes the Douglas-Kroll(-Hess) integrals for  
c:tex-a given basis set and stores them in a global array. 
c:tex-
c:tex-\noindent Author: W. A. de Jong
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1edk_ga(basis_in, g_in, type, kskel)
c:tex-\end{verbatim}
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
#include "rel_consts.fh"
#include "bgj.fh"
c     
c     Compute the Douglas-Kroll integrals and ADD them to a given global array.
c
c     There are three forms of Douglas-Kroll that are or can be used:
c     1. Free-particle projection operators, similar to Foldy-Wouthuysen 
c        transformation.
c     2. External potential projection operators without pxVp terms 
c        (vn2 approximation). This is also Bernd Hess' implementation.
c     3. External potential projection operators, i.e  complete Douglas-Kroll.
c
c     Method:
c     A. Determine uncontracted basis (dk_fitbas) from input basis (basis_in)
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
c:tex-{\it Argument list:}
c:tex-\begin{verbatim}
      integer basis_in          ! [input] basis sets for which integrals are computed
      integer g_in              ! [output] GA handle to array of final integrals
      integer type              ! [input] type of int, 1=kinetic, 2=potential, 3=both
      logical kskel             ! [input] Skeleton integrals, yes or no?
c:tex-\end{verbatim}
c     
c     local variables
c     
      integer nbf, nbf_in       ! # of basis functions for dk_fitbas and basis_in
      integer Tg_a              ! GA handle kinetic energy integral matrix T 
      integer Ug_a              ! GA handle transformation matrix r -> p space 
      integer Vg_a              ! GA handle potential energy integral matrix V
      integer Vdg_a             ! GA handle p.Vp and later pxVP integral matrix
      integer DKg_a             ! GA handle Douglas-Kroll integrals in uncontracted basis
      integer X_ga, Y_ga, Z_ga  ! GA arrays for th epxVp integrals (if needed)
      integer tempg_a           ! Temporary g_a matrix #1 handle
      integer tmpg_a            ! Temporary g_a matrix #2 handle
      integer Mixg_a            ! GA handle for transformation between contr and uncont
      integer p2m_a, p2Im_a     ! MA handle and index p^2 vector 
      integer Em_a, EIm_a       ! MA handle and index E vector
      integer Am_a, AIm_a       ! MA handle and index A vector
      integer Km_a, KIm_a       ! MA handle and index K vector
      integer tempm_a, tempIm_a ! Temporary m_a vector handle and index
      integer iloop             
      logical oskel            
c
c     cau = 137.0360000d0         ! Speed of light by Bernd Hess
c
c     At this point we set oskel to false. We cannot use oskel because of the
c     transformations we have to do with the full integral matrices.
c
      oskel = .false.
c
c     Normalize the uncontracted fitting basis 
c
      call int_normalize(dk_fitbas)
c
c     Define the dimensions of the arrays
c
      if (.not. bas_numbf(dk_fitbas,nbf))
     $   call errquit('int_1edk_ga: bas_numbf failed',0)
c
c     Create DKg_a of nbf*nbf to store integrals of uncontracted basis
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA D-K',1,1,DKg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA D-K',0)
      call ga_zero(DKg_a)
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
      call int_1e_oldga(dk_fitbas,dk_fitbas,Tg_a,'kinetic',oskel)
c
c     Diagonalize T in the orthogonal basis, using the overlap S as metric to get XU
c
c     Put eigenvalues in p2g_a and multiply by 2.0 to get p^2,
c     destroy T because we don't need that one anymore
c
      call ga_zero(tempg_a)
      call int_1e_oldga(dk_fitbas,dk_fitbas,tempg_a,'overlap',oskel)
      if (.not. ma_push_get(MT_DBL,nbf,'psquared',p2m_a,p2Im_a))
     $     call errquit('int_1edk_ga: ma_push_get failed psquared',0)
#ifdef PARALLEL_DIAG
      call ga_diag(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
#else
      call ga_diag_seq(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
#endif

      if (.not. ga_destroy(Tg_a)) 
     $     call errquit('int_1edk_ga: ga_destroy Tg_a',0)
      call dscal(nbf,2.0d0,dbl_mb(p2Im_a),1)
c
c     Construct vector E from eigenvalues: E = Sqrt(p^2 c^2 + c^4)
c                                            = c*Sqrt(p^2 + c^2)
c
      if (.not. ma_push_get(MT_DBL,nbf,'Evector',Em_a,EIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Evector',0)
      do 10 iloop = 1, nbf
         dbl_mb(EIm_a+iloop-1) = cau*Sqrt(dbl_mb(p2Im_a+iloop-1)+csq)
   10 continue
c
c     Construct vector A from E: A = Sqrt((E+c^2)/2E) and
c               vector K from E: K = c/(E+c^2)
c
      if (.not. ma_push_get(MT_DBL,nbf,'Avector',Am_a,AIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Avector',0)
      if (.not. ma_push_get(MT_DBL,nbf,'Kvector',Km_a,KIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Kvector',0)
      do 20 iloop = 1, nbf
         dbl_mb(KIm_a+iloop-1) = cau/(dbl_mb(EIm_a+iloop-1)+csq)  
         dbl_mb(AIm_a+iloop-1) = Sqrt((dbl_mb(EIm_a+iloop-1)+csq)/
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
      call ga_add_todiag(DKg_a,-csq)
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
      call int_1e_oldga(dk_fitbas,dk_fitbas,Vg_a,'potential',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Vg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Vg_a)
c
c     Do A Vg_a A, remember A is diagonal stored as vector and store in DKg_a
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,DKg_a)
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
      call int_1e_oldga(dk_fitbas,dk_fitbas,Vdg_a,'pVp',oskel)
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
c     construct final integrals if dk_form_type=1 (see header).
c 
      if (dk_form_type .eq. 1) then     
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0)
         goto 100
      endif
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
      call ga_scale_lh(tempg_a,dbl_mb(KIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     Term 2: - UA Vg_a AE(p)AK Vdg_a KAU^t
c     (Same intermediate array as Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_scale_rh(tempg_a,dbl_mb(KIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a)
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
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
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
      if (.not. ga_destroy(Vdg_a))
     $     call errquit('int_1edk_ga: ga_destroy Vdg_a',0)
      if (.not. ga_destroy(Vg_a))
     $     call errquit('int_1edk_ga: ga_destroy Vg_a',0)
c
c     At this point we have the basic integrals for the external potential projection
c     operators without the pxVp type terms. These remaining pxVp terms come from 
c     are needed for the full Douglas-Kroll set of integrals. If we neglect these 
c     terms we get an incomplete Douglas-Kroll transformation a la Haberlen&Rosch 
c     (CPL, 199, 491, 1992) denoted vn2.
c     Skip and construct final integrals if dk_form_type=2 (see header).
c 
      if (dk_form_type .eq. 2) then
         if (.not. ma_pop_stack(tempm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack tempm_a',0)
         goto 100
      endif
c
c     Unfortunately for the full Douglas-Kroll we do need for each term 
c     integrals of the type (pxVp)D(pxVp). We now need three GA's to 
c     store the three components of the integrals and to construct the 
c     final three terms.
c
c     NOTE: If we neglect these terms we get an incomplete Douglas-Kroll
c           transformation a la Haberlen&Rosch (CPL, 199, 491, 1992)
c
c
c     Get pxVp X_ga, Y_ga and Z_ga, and transform to U^tVU and scale 
c     with 1/(Ep+Ep') for W terms
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp X',1,1,X_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA X_ga',0)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Y',1,1,Y_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA Y_ga',0)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Z',1,1,Z_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA Z_ga',0)
c
      call int_pxvp_ga(dk_fitbas,X_ga,Y_ga,Z_ga,oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,X_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,X_ga)
      call ga_scale_divEpp(X_ga,dbl_mb(EIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Y_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Y_ga)
      call ga_scale_divEpp(Y_ga,dbl_mb(EIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Z_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Z_ga)
      call ga_scale_divEpp(Z_ga,dbl_mb(EIm_a))
c
c     W(1)E(p)W(1): - UAK     pxVp  KAE(p)K^-2p^-2AK pxVp  KAU^t
c
c     For each component of pxVp (i.e. x, y, z)
c
      do 60 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(EIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   60 continue
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmpg_a)
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,1.0d0,
     $              tmpg_a)
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,1.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     E(p)W(1)^2:   - UE(p)AK pxVp  KAK^-2p^-2AK     pxVp  KAU^t
c
c     For each component of pxVp (i.e. x, y, z)
c
      do 61 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   61 continue
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmpg_a)
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,1.0d0,
     $              tmpg_a)
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,1.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     W(1)^2E(p):   - UAK     pxVp  KAK^-2p^-2AK     pxVp  KAE(p)U^t
c
c     For each component of pxVp (i.e. x, y, z)
c
      do 62 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
   62 continue
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmpg_a)
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,1.0d0,
     $              tmpg_a)
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,1.0d0,
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
      if (.not. ga_destroy(Z_ga))
     $    call errquit('int_1edk_ga: ga_destroy Z_ga',0)
      if (.not. ga_destroy(Y_ga))
     $    call errquit('int_1edk_ga: ga_destroy Y_ga',0)
      if (.not. ga_destroy(X_ga))
     $    call errquit('int_1edk_ga: ga_destroy X_ga',0)
      if (.not. ma_pop_stack(tempm_a))
     $    call errquit('int_1edk_ga: ma_pop_stack tempm_a',0)
  100 continue
c
c     Get overlap and multiply with Ug_a to get SXU
c
      call ga_zero(tempg_a)
      call ga_zero(tmpg_a)
      call int_1e_oldga(dk_fitbas, dk_fitbas, tempg_a,'overlap',.false.)
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
c     Create mixed overlap between basis_in and dk_fitbas
c
      call ga_zero(Mixg_a)
      call int_1e_oldga(dk_fitbas,basis_in,Mixg_a,'overlap',.false.)
c
c     Create overlap matrix dk_fitbas, invert to S^-1 and multiply with Mixg_a
c
      call ga_zero(tempg_a)
      call int_1e_oldga(dk_fitbas,dk_fitbas,tempg_a,'overlap',.false.)
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
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf_in,nbf,nbf,1.0d0,Mixg_a,
     $              DKg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf_in,nbf_in,nbf,1.0d0,tempg_a,Mixg_a,
     $              1.0d0,g_in)
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
      return
c
      end
