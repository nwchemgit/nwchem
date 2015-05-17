c $Id$
C> \ingroup nwint
C> @{
C>
C> \brief Compute Douglas-Kroll(-Hess) integrals and store them
C> in a global array
C>    
C> Compute the Douglas-Kroll integrals and ADD them to a given global array.
C>
C> There are three forms of Douglas-Kroll that are or can be used:
C> * 1. Free-particle projection operators, similar to Foldy-Wouthuysen 
C>   transformation.
C> * 2. External potential projection operators without pxVp terms 
C>   (vn2 approximation). This is also Bernd Hess' implementation.
C> * 3. External potential projection operators, i.e  complete Douglas-Kroll.
C>
C> Method:
C> * A. Determine uncontracted basis (dk_fitbas) from input basis (basis_in)
C> * B. Create p^2/2 matrix and diagonalize to get real to momentum space 
C>   transformation.
C> * C. Create diagonal matrices (store as vectors) of momentum operators in 
C>   momentum space.
C> * D. Compute all terms of E, V and p.Vp terms that are needed for form 1 
C>   (see above).
C> * E. Compute all combined terms that are needed for form 2 (see above) 
C>   and add to DKg_a.
C> * F. Compute the three W terms with pxVp required in form 3 (see above) 
C>   and add them to the DKg_a array
C> * G. Transform DKg_a to real space
C> * H. Transform DKg_a from uncontracted basis to contracted basis
C>
C> Author: W. A. de Jong
C>
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
#include "errquit.fh"
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
      integer basis_in          !< [Input] basis sets for which integrals are computed
      integer g_in              !< [Output] GA handle to array of final integrals
      integer type              !< [Input] type of int,
                                !< * 1=kinetic,
                                !< * 2=potential,
                                !< * 3=both
      logical kskel             !< [Input] Skeleton integrals, yes or no?
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
      integer X_ga, Y_ga, Z_ga  ! GA arrays for the pxVp integrals (if needed)
      integer tempg_a           ! Temporary g_a matrix #1 handle
      integer tmpg_a            ! Temporary g_a matrix #2 handle
      integer Mixg_a            ! GA handle for transformation between contr and uncont
      integer p2m_a, p2Im_a     ! MA handle and index p^2 vector 
      integer Em_a, EIm_a       ! MA handle and index E vector
      integer Am_a, AIm_a       ! MA handle and index A vector
      integer Km_a, KIm_a       ! MA handle and index K vector
      integer KEm_a, KEIm_a     ! MA handle and index E-mc^2 vector
      integer tempm_a, tempIm_a ! Temporary m_a vector handle and index
      integer iloop             
      logical oskel            
      integer Ev1g_a            ! GA handle first-order even term in uncontracted basis
      integer RE1Rg_a           ! GA handle R*E1*R
      integer tmp2g_a           ! Temporary g_a matrix #3 handle for DK3
      integer XpxVp_ga          ! Unmodified pxVp for DK3
      integer YpxVp_ga          ! Unmodified pxVp for DK3
      integer ZpxVp_ga          ! Unmodified pxVp for DK3
      integer tmp3g_a           ! Temporary g_a matrix #3 handle for
                                ! full DK3
      integer info
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
c
c     cau = 137.0360000d0         ! Speed of light by Bernd Hess
c
c     At this point we set oskel to false. We cannot use oskel because of the
c     transformations we have to do with the full integral matrices.
c
      oskel = .false.
c
c     Define the dimensions of the arrays
c
      if (.not. bas_numbf(dk_fitbas,nbf))
     $   call errquit('int_1edk_ga: bas_numbf failed',0, BASIS_ERR)
c
c     Create DKg_a of nbf*nbf to store integrals of uncontracted basis
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA D-K',1,1,DKg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA D-K',0,
     &       GA_ERR)
      call ga_zero(DKg_a)
c
c     Create a temp g_a of nbf*nbf as temp/help matrix array
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 1',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 2',1,1,tmpg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
c
c     Create transformation matrix r-space to g-space, U
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Umatrix',1,1,Ug_a))
     $     call errquit('int_1edk_ga: ga_create failed Umatrix',0,
     &       GA_ERR)
      call ga_zero(Ug_a)
c
c     Get kinetic energy matrix p^2/2 in uncontracted basis, diagonalize
c     Store transformation matrix and vector of eigenvalues
c
c     Construct, fill Tg_a with kinetic energy matrix integrals
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Tmatrix',1,1,Tg_a))
     $     call errquit('int_1edk_ga: ga_create failed Tmatrix',0,
     &       GA_ERR)
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
     $     call errquit('int_1edk_ga: ma_push_get failed psquared',0,
     &       MA_ERR)
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
c     pdsygv is unstable since it uses the inaccurate pdsyevx
cold         call ga_pdsygv(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
c     S=chol*chol_t
      info= ga_cholesky('L',tempg_a)
c     inv(a)
      info = ga_llt_i('L',tempg_a,-1)
c     U=chol(inv(a))_t
      info = ga_cholesky('L',tempg_a,-1)

c     f_prime=U_t*f*inv(chol)*U
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,tempg_a,Tg_a,0.0d0,tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,tempg_a,
     C     0.0d0,tg_a)
c     diag f_prime c_prime
      call ga_pdsyev(tg_a, tmpg_a, dbl_mb(p2im_a), 0)
c     c=inv(chol)_t*c_prime
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,tmpg_a,
     Z     0.0d0,ug_a)
#else
      call ga_diag(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
#endif
#else
      call ga_diag_seq(Tg_a,tempg_a,Ug_a,dbl_mb(p2Im_a))
#endif
c
      if (.not. ga_destroy(Tg_a)) 
     $     call errquit('int_1edk_ga: ga_destroy Tg_a',0,
     &       GA_ERR)
      call dscal(nbf,2.0d0,dbl_mb(p2Im_a),1)
c
c     Construct vector E from eigenvalues: E = Sqrt(p^2 c^2 + c^4)
c                                            = c*Sqrt(p^2 + c^2)
c
      if (.not. ma_push_get(MT_DBL,nbf,'Evector',Em_a,EIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Evector',0,
     &       MA_ERR)
      do 10 iloop = 1, nbf
         dbl_mb(EIm_a+iloop-1) = cau*Sqrt(dbl_mb(p2Im_a+iloop-1)+csq)
   10 continue
c
c     Construct vector A from E: A = Sqrt((E+c^2)/2E) and
c               vector K from E: K = c/(E+c^2)
c
      if (.not. ma_push_get(MT_DBL,nbf,'Avector',Am_a,AIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Avector',0,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL,nbf,'Kvector',Km_a,KIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed Kvector',0,
     &       MA_ERR)
      if (dk_form_type .ne. 0) then   ! Douglas-Kroll     
      do 20 iloop = 1, nbf
         dbl_mb(KIm_a+iloop-1) = cau/(dbl_mb(EIm_a+iloop-1)+csq)  
         dbl_mb(AIm_a+iloop-1) = Sqrt((dbl_mb(EIm_a+iloop-1)+csq)/
     $                                (2.0d0*dbl_mb(EIm_a+iloop-1)))
   20 continue
      else  
      do iloop = 1, nbf
         dbl_mb(KIm_a+iloop-1) = cau/(dbl_mb(EIm_a+iloop-1)+csq)  
         dbl_mb(AIm_a+iloop-1) = (Sqrt(1.0d0+csq*dbl_mb(p2Im_a+iloop-1)/
     $     ((dbl_mb(EIm_a+iloop-1)+csq)*(dbl_mb(EIm_a+iloop-1)+csq))))/
     $                                  (dbl_mb(EIm_a+iloop-1)+csq)
      enddo    
      endif
c
c     Clear given global array and put in relativistic kinetic energy E on
c     diagonal (E matrix diagonal in momentum space)
c
c     We have to subtract mc^2 to get the non-relativistic type integrals
c
      call ga_zero(DKg_a)
      if (.not. ma_push_get(MT_DBL,nbf,'Kinetic',KEm_a,KEIm_a))
     $     call errquit('int_1edk_ga: ma_push_get failed KEvector',0,
     &       MA_ERR)
      do 21 iloop = 1, nbf
         dbl_mb(KEIm_a+iloop-1) = (dbl_mb(p2Im_a+iloop-1)*csq)/
     $                            (dbl_mb(EIm_a+iloop-1)+csq)
   21 continue
      call ga_copy_todiag(DKg_a,dbl_mb(KEIm_a))
      if (.not. ma_pop_stack(KEm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack KEm_a',0,
     &       MA_ERR)
c
c     We have our kinetic energy integrals
c 
      if (type .eq. 1) goto 100
c
c     If we just need the potential energy integrals we have to clear Dkg_a
c
      if (type .eq. 2) call ga_zero(DKg_a)
c
c     construct RE1R for DK3
c
      if (dk_form_type .ge. 4) then   ! higher-order DK
      if (.not. ga_create(MT_DBL,nbf,nbf,'RE1Rmatrix',1,1,RE1Rg_a))
     $     call errquit('int_1edk_ga: ga_create failed RE1Rmatrix',0,
     &       GA_ERR)
      call ga_zero(RE1Rg_a)
      endif
c
c     Get potential energy matrix V in uncontracted basis, add AVA term to global array:
c
c     <phi|AVA|phi> = AU^tVUA
c
c     First create U^tVU in Vg_a
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Vmatrix',1,1,Vg_a))
     $     call errquit('int_1edk_ga: ga_create failed Vmatrix',0,
     &       GA_ERR)
      call ga_zero(Vg_a)
      call int_1e_oldga(dk_fitbas,dk_fitbas,Vg_a,'potential',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Vg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Vg_a)
      call ga_symmetrize(Vg_a)
c
c     Do A Vg_a A, remember A is diagonal stored as vector and store in DKg_a
c
      call ga_copy(Vg_a,tempg_a)
      if (dk_form_type .ne. 0) then   ! Douglas-Kroll     
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,DKg_a)
      if (dk_form_type .ge. 4) then   ! higher-order DK
c
c     Do R^2 (A Vg_a A) R^2 = K^2 p^2 (A Vg_a A) p^2 K^2, 
c     and store in RE1Rg_a for DK3
c
      call ga_scale_rlh(tempg_a,dbl_mb(p2Im_a),dbl_mb(p2Im_a))
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_dadd(1.0d0,RE1Rg_a,1.0d0,tempg_a,RE1Rg_a)
      endif
      else   ! RESC
      if (.not. ma_push_get(MT_DBL,nbf,'MA temp',tempm_a,tempIm_a))
     $     call errquit('D-K integrals: ma_push_get failed MA temp',0,
     &       MA_ERR)
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = Sqrt(dbl_mb(KIm_a+iloop-1))  
      enddo   
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = 1.0d0/dbl_mb(AIm_a+iloop-1)  
      enddo   
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,cau*2.0d0,tempg_a,DKg_a)   ! multiply 2c
      endif
c
c     Get double derivative potential energy matrix p.Vp in uncontracted basis, 
c     add ARVRA term to global array:
c
c     <phi|ARVRA|phi> = AKU^t(p.Vp)UKA
c
c     First create U^t(p.Vp)U in Vdg_a
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'Vdmatrix',1,1,Vdg_a))
     $     call errquit('int_1edk_ga: ga_create failed Vdmatrix',0,
     &       GA_ERR)
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
      call ga_symmetrize(Vdg_a)
c
c     Do A K Vdg_a K A, remember A and K are diagonal stored as vector and add to DKg_a
c
      call ga_copy(Vdg_a,tempg_a)
      if (dk_form_type .ne. 0) then   ! Douglas-Kroll     
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,DKg_a)
      if (dk_form_type .ge. 4) then   ! higher-order DK
c
c     add A K Vdg_a K A to RE1Rg_a for DK3
c
      call ga_dadd(1.0d0,RE1Rg_a,1.0d0,tempg_a,RE1Rg_a)
      endif
      else   ! RESC
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = 1.0d0/dbl_mb(AIm_a+iloop-1)  
      enddo   
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tempg_a,DKg_a)
      endif
c
c     copy current DKg_a (i.e.,1st-order even term) into Ev1g_a for DK3
c
      if (dk_form_type .ge. 4) then   ! higher-order DK
      if (.not. ga_create(MT_DBL,nbf,nbf,'Ev1matrix',1,1,Ev1g_a))
     $     call errquit('int_1edk_ga: ga_create failed Ev1matrix',0,
     &       GA_ERR)
      call ga_zero(Ev1g_a)
      call ga_copy(DKg_a,Ev1g_a)
      endif
c
c     At this point we have the basic integrals for the free-particle projection
c     operators. The remaining terms come from the external potential. Skip and 
c     construct final integrals if dk_form_type=1 (see header).
c 
      if (dk_form_type .eq. 1) then     
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0,
     &       GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
         goto 100
      endif
      if (dk_form_type .eq. 0) then   ! RESC
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0, GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
         if (.not. ma_pop_stack(tempm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack tempm_a',0, GA_ERR)
         goto 100
      endif
c
c     Scale U^tVU (in Vg_a) and U^tp.VpU (in Vdg_a) with 1/(Ep+Ep') for W terms
c
      call ga_scale_divEpp(Vg_a,dbl_mb(EIm_a))
      call ga_symmetrize(Vg_a)
      call ga_scale_divEpp(Vdg_a,dbl_mb(EIm_a))
      call ga_symmetrize(Vdg_a)
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
     $     call errquit('D-K integrals: ma_push_get failed MA temp',0,
     &       MA_ERR)
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
c
c     At this point we have the basic integrals for the external potential projection
c     operators without the pxVp type terms. These remaining pxVp terms come from 
c     are needed for the full Douglas-Kroll set of integrals. If we neglect these 
c     terms we get an incomplete Douglas-Kroll transformation a la Haberlen&Rosch 
c     (CPL, 199, 491, 1992) denoted vn2.
c     Skip and construct final integrals if dk_form_type=2 (see header).
c 
      if (dk_form_type .eq. 2) then
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0, GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
         if (.not. ma_pop_stack(tempm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack tempm_a',0, MA_ERR)
         goto 100
      endif
      if (dk_form_type .eq. 3) goto 200   ! complete DK2
c
c     Compute - W(1)E(1)W(1) term and add to global array:
c
c     - <phi|WEW|phi> = - UAKVdg_aKA R^-2 RE(1)R AVg_aAU^t     
c
c                       - UAVg_aA RE(1)R R^-2 AKVdg_aKAU^t 
c
c                       + UAVg_aA RE(1)R AVg_aAU^t
c
c                       + UAKVdg_aKA R^-2 RE(1)R R^-2 AKVdg_aKAU^t
c
c                       + UAK(pxV'p)KA R^-2 AK(pxVp)KA AVg_aAU^t (from term 1)
c
c                       + UAVg_aA AK(pxVp)KA R^-2 K(pxV'p)KAU^t (from term 2)
c
c                       - UAK(pxV'p)KA R^-2 RE(1)R R^-2 AK(pxV'p)KAU^t (from term 4)
c
c                       - UAK(pxV'p)KA R^-2 AK(pxVp)KA R^-2 AKVdg_aKAU^t (from term 4)
c
c                       - UAKVdg_aKA R^-2 AK(pxVp)KA R^-2 AK(pxV'p)KAU^t (from term 4)
c
c     (Last 5 terms done at the end, reusing p.Vp space)
c
c     Term 1: - UAKVdg_aKA R^-2 RE(1)R AVg_aAU^t
c           = - UKAVdg_aAK^-1p^-2 RE(1)R AVg_aAU^t
c     
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 3',1,1,tmp2g_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(KIm_a))
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) 
     $   = 1.0d0/dbl_mb(KIm_a+iloop-1)/dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,RE1Rg_a,0.0d0,
     $              tmpg_a)
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmp2g_a,DKg_a)
c
c     Term 2: - UAVg_aA RE(1)R R^-2 AKVdg_aKAU^t
c           = - UAVg_aA RE(1)R p^-2K^-1 AVdg_aAKU^t
c     
      call ga_zero(tempg_a)
      call ga_transpose(tmp2g_a,tempg_a)
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tempg_a,DKg_a)
c
c     Term 3: + UAVg_aA RE(1)R AVg_aAU^t  
c     
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,RE1Rg_a,0.0d0,
     $              tmpg_a)
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmp2g_a,DKg_a)
c
c     Term 4: + UAKVdg_aKA R^-2 RE(1)R R^-2 AKVdg_aKAU^t
c           = + UKAVdg_aA K^-1p^-2 RE(1)R p^-2K^-1 AVdg_aAKU^t
c     
      call ga_copy(RE1Rg_a,tempg_a)
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)/dbl_mb(KIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vdg_a,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmpg_a,DKg_a)
c
c     Compute + 1/2 E(1)W(1)^2 term and add to global array:
c
c     + <phi|EW^2|phi> = + UE(1) AKVdg_aKA AVg_aAU^t
c                     
c                        + UE(1) AVg_aA AKVdg_aKA U^t
c      
c                        - UE(1) AVg_aA p^2K^2 AVg_aAU^t
c
c                        - UE(1) AKVdg_aKA K^-2p^-2 AKVdg_aKAU^t
c
c                        - UAK(pxVp)KA AK(pxV'p)KA AVg_aAU^t (from term 1)
c                     
c                        - UAK(pxVp)KA AVg_aA AK(pxV'p)KA U^t (from term 2)
c      
c                        + UE(1) AK(pxV'p)KA K^-2p^-2 AK(pxV'p)KAU^t (from term 4)
c
c                        + UAK(pxVp)KA AK(pxV'p)KA K^-2p^-2 AKVdg_aKAU^t (from term 4)
c
c                        + UAK(pxVp)KA AKVdg_aKA K^-2p^-2 AK(pxV'p)KAU^t (from term 4)
c
c     (Last 5 terms done at the end, reusing p.Vp space)
c
c     Term 1: + UE(1)AK Vdg_a KAA Vg_a AU^t
c
      do 340 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
  340 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Ev1g_a,tmpg_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
c
c     Term 2: + UE(1)A Vg_a AAK Vdg_a KAU^t
c     (Intermediate array is the same as in Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Ev1g_a,tmpg_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
c
c     Term 3: - UE(1)A Vg_a Ap^2K^2A Vg_a AU^t
c
      do 341 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(p2Im_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
  341 continue
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Ev1g_a,tmpg_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tempg_a,DKg_a)
c
c     Term 4: - UE(1)AK Vdg_a KAK^-2p^-2AK Vdg_a KAU^t
c
      do 342 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
  342 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Ev1g_a,tmpg_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tempg_a,DKg_a)
c
c     Compute + 1/2 W(1)^2E(1) term and add to global array:
c
c     + <phi|W^2E|phi> = + UAKVdg_aKA AVg_aA E(1)U^t
c                      
c                        + UAVg_aA AKVdg_aKA E(1)U^t
c      
c                        - UAVg_aA p^2K^2 AVg_aA E(1)U^t
c
c                        - UAKVdg_aKA K^-2p^-2 AKVdg_aKA E(1)U^t
c                      
c                        - UAVg_aA AK(pxV'p)KA AK(pxVp)KAU^t (from term2)
c
c                        - UAK(pxV'p)KA AVg_aA AK(pxVp)KAU^t (from term1)
c                      
c                        + UAK(pxV'p)KA K^-2p^-2 AK(pxV'p)KA E(1)U^t (from term 4)
c                      
c                        + UAKVdg_aKA K^-2p^-2 AK(pxV'p)KA AK(pxVp)KAU^t (from term 4)
c
c                        + UAK(pxV'p)KA K^-2p^-2 AKVdg_aKA AK(pxVp)KAU^t (from term 4)
c
c     (Last 5 terms done at the end, reusing p.Vp space)
c
c     Term 1: + UAK Vdg_a KAA Vg_a AE(1)U^t
c
      do 350 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
  350 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,Ev1g_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
c
c     Term 2: + UA Vg_a AAK Vdg_a KAE(1)U^t
c     (Intermediate array is the same as in Term 1)
c
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Vg_a,tempg_a,0.0d0,
     $              tmpg_a) 
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,Ev1g_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
c
c     Term 3: - UA Vg_a Ap^2K^2A Vg_a AE(1)U^t
c
      do 351 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(p2Im_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(KIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)
  351 continue
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,Ev1g_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tempg_a,DKg_a)
c
c     Term 4: - UAK Vdg_a KAK^-2p^-2AK Vdg_a KAE(1)U^t
c
      do 352 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = dbl_mb(AIm_a+iloop-1)*
     $                              dbl_mb(AIm_a+iloop-1)/
     $                              dbl_mb(p2Im_a+iloop-1)
  352 continue
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Vdg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,Ev1g_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tempg_a,DKg_a)
c
      if (dk_form_type .eq. 4) then   ! incomplete DK3
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0, GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
         if (.not. ma_pop_stack(tempm_a))
     $       call errquit('int_1edk_ga: ma_pop_stack tempm_a',0, MA_ERR)
         if (.not. ga_destroy(Ev1g_a))
     $        call errquit('int_1edk_ga: ga_destroy W1sqg_a',0, GA_ERR)
         if (.not. ga_destroy(RE1Rg_a))
     $        call errquit('int_1edk_ga: ga_destroy RE1Rg_a',0, GA_ERR)
         if (.not. ga_destroy(tmp2g_a))
     $        call errquit('int_1edk_ga: ga_destroy tmp2g_a',0, GA_ERR)
         goto 100
      endif
  200 continue
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
     $     call errquit('int_1edk_ga: ga_create failed GA X_ga',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Y',1,1,Y_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA Y_ga',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Z',1,1,Z_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA Z_ga',0,
     &       GA_ERR)
      if (dk_form_type .eq. 5) then   ! complete DK3
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp X',1,1,XpxVp_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA XpxVp_ga',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Y',1,1,YpxVp_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA YpxVp_ga',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp Z',1,1,ZpxVp_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA ZpxVp_ga',0,
     &       GA_ERR)
      endif
c
      call int_pxvp_ga(dk_fitbas,X_ga,Y_ga,Z_ga,oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,X_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,X_ga)
      if (dk_form_type .ge. 5) call ga_copy(X_ga,XpxVp_ga)   ! store unmodified pxVp into XpxVp for DK3
      call ga_scale_divEpp(X_ga,dbl_mb(EIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Y_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Y_ga)
      if (dk_form_type .ge. 5) call ga_copy(Y_ga,YpxVp_ga)   ! store unmodified pxVp into XpxVp for DK3
      call ga_scale_divEpp(Y_ga,dbl_mb(EIm_a))
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,Z_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,Z_ga)
      if (dk_form_type .ge. 5) call ga_copy(Z_ga,ZpxVp_ga)   ! store unmodified pxVp into XpxVp for DK3
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
      if (dk_form_type .eq. 3) then   ! complete DK2
         if (.not. ga_destroy(Z_ga))
     $        call errquit('int_1edk_ga: ga_destroy Z_ga',0, GA_ERR)
         if (.not. ga_destroy(Y_ga))
     $        call errquit('int_1edk_ga: ga_destroy Y_ga',0, GA_ERR)
         if (.not. ga_destroy(X_ga))
     $        call errquit('int_1edk_ga: ga_destroy X_ga',0, GA_ERR)
         if (.not. ma_pop_stack(tempm_a))
     $        call errquit('int_1edk_ga: ma_pop_stack tempm_a',0,
     &       MA_ERR)
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0, GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
         goto 100
      endif
c
c     W(1)E(1)W(1): cross-product term 7
c
c     - UAK(pxV'p)KA R^-2 RE(1)R R^-2 AK(pxV'p)KAU^t
c   = - UKA(pxV'p) AK^-1p^-2 RE(1)R p^-2K^-1A (pxV'p)AKU^t
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp 4',1,1,tmp3g_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)/dbl_mb(KIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(RE1Rg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,X_ga,0.0d0,
     $              tmp3g_a)
      call ga_copy(tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Y_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Y_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Z_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Z_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     W(1)E(1)W(1): cross-product terms 8 + 9
c
c     - UAK(pxV'p)KA R^-2 AK(pxVp)KA R^-2 AKVdg_aKAU^t
c     - UAKVdg_aKA   R^-2 AK(pxVp)KA R^-2 AK(pxV'p)KAU^t
c   = - UKA(pxV'p) A^2p^-2(pxVp)p^-2A^2 Vdg_aAKU^t
c     - UKAVdg_a   A^2p^-2(pxVp)p^-2A^2 (pxV'p)AKU^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(XpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(YpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Y_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(ZpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Z_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
      call ga_copy(Vdg_a,tempg_a)
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp3g_a,tempg_a,0.0d0,
     $              tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
      call ga_copy(XpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(YpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(ZpxVp_ga,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
      call ga_copy(Vdg_a,tempg_a)
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,tmp3g_a,0.0d0,
     $              tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-1.0d0,tmpg_a,DKg_a)
c
c     W(1)E(1)W(1): cross-product terms 5 + 6
c
c     + UAK(pxV'p)KA R^-2 AK(pxVp)KA AVg_aAU^t
c     + UAVg_aA AK(pxVp)KA R^-2 K(pxV'p)KAU^t
c   = + UAK(pxV'p)p^-2A^2(pxVp)AK AVg_aAU^t
c     + UAVg_aA KA(pxVp)A^2p^-2(pxV'p)KAU^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(XpxVp_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(YpxVp_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Y_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(ZpxVp_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Z_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_scale_rlh(tmp3g_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmp3g_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp3g_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmpg_a,DKg_a)
c
      call ga_copy(XpxVp_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(YpxVp_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(ZpxVp_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_scale_rlh(tmp3g_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmp3g_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,tmp3g_a,0.0d0,
     $              tmpg_a)
      call ga_dadd(1.0d0,DKg_a,1.0d0,tmpg_a,DKg_a)
c
c     E(1)W(1)W(1) + W(1)W(1)E(1): cross-product term 7
c
c     + UE(1) AK(pxV'p)KA K^-2p^-2 AK(pxV'p)KAU^t
c     + UAK(pxV'p)KA K^-2p^-2 AK(pxV'p)KA E(1)U^t
c   = + UE(1) AK(pxV'p)A^2p^-2(pxV'p)KAU^t
c     + UAK(pxV'p)A^2p^-2(pxV'p)KA E(1)U^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,X_ga,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmpg_a)
c
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Y_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp2g_a,tmpg_a)
c
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Z_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp2g_a,tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Ev1g_a,tmpg_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,Ev1g_a,0.0d0,
     $              tempg_a)
      call ga_dadd(1.0d0,DKg_a,0.5d0,tempg_a,DKg_a)
c
c     E(1)W(1)W(1) + W(1)W(1)E(1): cross-product term 8
c
c    + UAK(pxVp)KA AK(pxV'p)KA K^-2p^-2 AKVdg_aKAU^t
c     + UAKVdg_aKA K^-2p^-2 AK(pxV'p)KA AK(pxVp)KAU^t
c   = + UAK(pxVp)K^2A^2(pxV'p)A^2p^-2Vdg_aKAU^t
c     + UAKVdg_ap^-2A^2(pxV'p)A^2K^2(pxVp)KAU^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)*dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(X_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,XpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,YpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,ZpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp3g_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)*dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,XpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,YpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,ZpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,tmp3g_a,0.0d0,
     $              tmpg_a)
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     E(1)W(1)W(1) + W(1)W(1)E(1): cross-product term 9
c
c     + UAK(pxVp)KA AKVdg_aKA K^-2p^-2 AK(pxV'p)KAU^t
c     + UAK(pxV'p)KA K^-2p^-2 AKVdg_aKA AK(pxVp)KAU^t
c   = + UAK(pxVp)K^2A^2Vdg_aA^2p^-2(pxV'p)KAU^t
c     + UAK(pxV'p)p^-2A^2Vdg_aA^2K^2(pxVp)KAU^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)*dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,XpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,X_ga,0.0d0,
     $              tmp3g_a)
      call ga_copy(tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,YpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Y_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,ZpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Z_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)*dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(Vdg_a,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    /dbl_mb(p2Im_a+iloop-1)
      enddo
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,XpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,YpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Y_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,ZpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Z_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,0.5d0,tmpg_a,DKg_a)
c
c     E(1)W(1)W(1) + W(1)W(1)E(1): cross-product term 5
c
c     - UAK(pxVp)KA AK(pxV'p)KA AVg_aAU^t
c     - UAVg_aA AK(pxV'p)KA AK(pxVp)KAU^t
c   = - UAK(pxVp)K^2A^2(pxV'p)KA AVg_aAU^t
c     - UAVg_aA AK(pxV'p)K^2A^2(pxVp)KAU^t
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)*dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(X_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,XpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,YpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_lh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,ZpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_scale_rlh(tmp3g_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmp3g_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp3g_a,tempg_a,0.0d0,
     $              tmpg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
      call ga_copy(X_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,XpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_copy(tmp2g_a,tmp3g_a)
c
      call ga_copy(Y_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,YpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_copy(Z_ga,tempg_a)
      call ga_scale_rh(tempg_a,dbl_mb(tempIm_a))
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,ZpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_dadd(1.0d0,tmp3g_a,1.0d0,tmp2g_a,tmp3g_a)
c
      call ga_scale_rlh(tmp3g_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmp3g_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_zero(tmpg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,tmp3g_a,0.0d0,
     $              tmpg_a)
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
c     E(1)W(1)W(1) + W(1)W(1)E(1): cross-product term 6
c
c     - UAK(pxVp)KA AVg_aA AK(pxV'p)KA U^t 
c     - UAK(pxV'p)KA AVg_aA AK(pxVp)KAU^t 
c   = - UAK(pxVp)  KA^2Vg_aA^2K (pxV'p)KA U^t 
c     - UAK(pxV'p) KA^2Vg_aA^2K (pxVp)KAU^t 
c
      do iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1)
     $   = dbl_mb(AIm_a+iloop-1)*dbl_mb(AIm_a+iloop-1)
     $    *dbl_mb(KIm_a+iloop-1)
      enddo
      call ga_copy(Vg_a,tempg_a)
      call ga_scale_rlh(tempg_a,dbl_mb(tempIm_a),dbl_mb(tempIm_a))
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,XpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,X_ga,0.0d0,
     $              tmp3g_a)
      call ga_copy(tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,YpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Y_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,ZpxVp_ga,tempg_a,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmp2g_a,Z_ga,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,XpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,YpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Y_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_zero(tmp2g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,ZpxVp_ga,0.0d0,
     $              tmp2g_a)
      call ga_zero(tmp3g_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,Z_ga,tmp2g_a,0.0d0,
     $              tmp3g_a)
      call ga_dadd(1.0d0,tmpg_a,1.0d0,tmp3g_a,tmpg_a)
c
      call ga_scale_rlh(tmpg_a,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(tmpg_a,dbl_mb(AIm_a),dbl_mb(AIm_a))
      call ga_dadd(1.0d0,DKg_a,-0.5d0,tmpg_a,DKg_a)
c
      if (dk_form_type .eq. 5) then   ! complete DK3
         if (.not. ga_destroy(Z_ga))
     $        call errquit('int_1edk_ga: ga_destroy Z_ga',0, GA_ERR)
         if (.not. ga_destroy(Y_ga))
     $        call errquit('int_1edk_ga: ga_destroy Y_ga',0, GA_ERR)
         if (.not. ga_destroy(X_ga))
     $        call errquit('int_1edk_ga: ga_destroy X_ga',0, GA_ERR)
         if (.not. ma_pop_stack(tempm_a))
     $        call errquit('int_1edk_ga: ma_pop_stack tempm_a',0,
     &       MA_ERR)
         if (.not. ga_destroy(Ev1g_a))
     $        call errquit('int_1edk_ga: ga_destroy Ev1g_a',0, GA_ERR)
         if (.not. ga_destroy(RE1Rg_a))
     $        call errquit('int_1edk_ga: ga_destroy RE1Rg_a',0, GA_ERR)
         if (.not. ga_destroy(tmp2g_a))
     $        call errquit('int_1edk_ga: ga_destroy tmp2g_a',0, GA_ERR)
         if (.not. ga_destroy(XpxVp_ga))
     $        call errquit('int_1edk_ga: ga_destroy XpxVp_ga',0, GA_ERR)
         if (.not. ga_destroy(YpxVp_ga))
     $        call errquit('int_1edk_ga: ga_destroy YpxVp_ga',0, GA_ERR)
         if (.not. ga_destroy(ZpxVp_ga))
     $        call errquit('int_1edk_ga: ga_destroy ZpxVp_ga',0, GA_ERR)
         if (.not. ga_destroy(tmp3g_a))
     $        call errquit('int_1edk_ga: ga_destroy tmp3g_a',0, GA_ERR)
         if (.not. ga_destroy(Vdg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vdg_a',0, GA_ERR)
         if (.not. ga_destroy(Vg_a))
     $        call errquit('int_1edk_ga: ga_destroy Vg_a',0, GA_ERR)
      endif
  100 continue
c
c     Symmetrize final set of integrals
c
      call ga_symmetrize(DKg_a)
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
      call ga_symmetrize(DKg_a)
c
c     Free MA's and GA's that are not needed anymore at this time
c
      if (.not. ma_pop_stack(Km_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Km_a',0, MA_ERR)
      if (.not. ma_pop_stack(Am_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Am_a',0, MA_ERR)
      if (.not. ma_pop_stack(Em_a))
     $     call errquit('int_1edk_ga: ma_pop_stack Em_a',0, MA_ERR)
      if (.not. ma_pop_stack(p2m_a))
     $     call errquit('int_1edk_ga: ma_pop_stack p2m_a',0, MA_ERR)
      if (.not. ga_destroy(Ug_a))
     $     call errquit('int_1edk_ga: ga_destroy Ug_a',0, GA_ERR)
      if (.not. ga_destroy(tmpg_a))
     $     call errquit('int_1edk_ga: ga_destroy tmpg_a',0, GA_ERR)
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0, GA_ERR)
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
     $   call errquit('int_1edk_ga: bas_numbf failed',0, BASIS_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf_in,'Mixed overlap',1,1,
     $     Mixg_a)) call errquit('int_1edk_ga: create mixed overlap',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA temp',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL,nbf,nbf_in,'GA temp',1,1,tmpg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA tmp',0,
     &       GA_ERR)
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
     $     call errquit('int_1edk_ga: ga_destroy tmpg_a',0, GA_ERR)
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0, GA_ERR)
      if (.not. ga_create(MT_DBL,nbf_in,nbf,'GA temp',1,1,tempg_a))
     $     call errquit('int_1edk_ga: ga_create failed GA temp',0,
     &       GA_ERR)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf_in,nbf,nbf,1.0d0,Mixg_a,
     $              DKg_a,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf_in,nbf_in,nbf,1.0d0,tempg_a,Mixg_a,
     $              1.0d0,g_in)
      call ga_symmetrize(g_in)
c
c     Destroy transformation matrices Mixg_a, tempg_a and DKg_a
c
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0, GA_ERR)
      if (.not. ga_destroy(Mixg_a))
     $     call errquit('int_1edk_ga: ga_destroy Mixg_a',0, GA_ERR)
      if (.not. ga_destroy(DKg_a))
     $     call errquit('int_1edk_ga: ga_destroy DKg_a',0, GA_ERR)
c
      return
c
      end
C>
C> @}
