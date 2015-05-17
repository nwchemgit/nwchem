c $Id$
C> \ingroup nwint
C> @{
C>
C> Compute the one-electron spin-orbit Douglas-Kroll(-Hess) integrals
C> and store them in a global array
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
c:tex-This routine computes the one-electron spin-orbit
c:tex-Douglas-Kroll(-Hess) integrals for  
c:tex-a given basis set and stores them in a global array. 
c:tex-
c:tex-\noindent Author: W. A. de Jong
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1edk_so_ga(basis_in, g_in, type, kskel)
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
      integer g_in(3)           !< [Output] GA handle to array of final integrals
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
      integer V_ga              ! GA handle potential energy integral matrix V
      integer DKg_a             ! GA handle Douglas-Kroll integrals in uncontracted basis
      integer X_ga(3)           ! GA arrays for the pxVp integrals 
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
      integer C_ga(3)           ! Collected final spin-orbit integrals
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
      if(info.ne.0) call errquit('int1edk: cholesky ',0,0)
c     inv(a)
      info = ga_llt_i('L',tempg_a,-1)
      if(info.ne.0) call errquit('int1edk: cholesky ',0,0)
c     U=chol(inv(a))_t
      info = ga_cholesky('L',tempg_a,-1)
      if(info.ne.0) call errquit('int1edk: cholesky ',0,0)

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
      do 20 iloop = 1, nbf
         dbl_mb(KIm_a+iloop-1) = cau/(dbl_mb(EIm_a+iloop-1)+csq)  
         dbl_mb(AIm_a+iloop-1) = Sqrt((dbl_mb(EIm_a+iloop-1)+csq)/
     $                                (2.0d0*dbl_mb(EIm_a+iloop-1)))
   20 continue
c
c     Start spin-orbit integral assembly
c
      call ga_zero(g_in(1))
      call ga_zero(g_in(2))
      call ga_zero(g_in(3))
c
c     Get collection matrices C_ga(1), (2), and (3)
c
      do iloop = 1, 3
        if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp X',1,1,C_ga(iloop)))
     $      call errquit('int_1edk_ga: ga_create failed GA C_ga',iloop,
     &       GA_ERR)
      enddo
c
c     Get pxVp matrices X_ga, Y_ga and Z_ga, and transform to U^tVU 
c
      do iloop = 1, 3
        if (.not. ga_create(MT_DBL,nbf,nbf,'GA pxVp X',1,1,X_ga(iloop)))
     $      call errquit('int_1edk_ga: ga_create failed GA X_ga',iloop,
     &       GA_ERR)
      enddo
c
      call int_pxvp_ga(dk_fitbas,X_ga(1),X_ga(2),X_ga(3),oskel)
      do iloop = 1, 3
         call ga_zero(tempg_a)
         call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,X_ga(iloop),
     $                 0.0d0,tempg_a)
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,
     $                 X_ga(iloop))
      enddo
c
c     <phi|AR<pxVp>RA|phi> = AKU^t(<pxVp>)UKA
c
c     Do A K pxVp K A, remember A and K are diagonal stored as vector and add to DKg_a
c
      do iloop = 1, 3
         call ga_scale_rlh(X_ga(iloop),dbl_mb(KIm_a),dbl_mb(KIm_a))
         call ga_scale_rlh(X_ga(iloop),dbl_mb(AIm_a),dbl_mb(AIm_a))
      enddo
c
c     Copy terms to collection matrices
c
      do iloop = 1, 3
         call ga_copy(X_ga(iloop),C_ga(iloop))
      enddo
c
c     At this point we have the basic integrals for the free-particle projection
c     operators. The remaining terms come from the external potential. Skip and
c     construct final integrals if dk_form_type=1 (see header).
c
      if (dk_form_type .eq. 1) goto 1002
c
c     Start on the W terms
c
c     The pxVp integrals now will be scaled with 1/(Ep+Ep') for W terms
c
      do iloop = 1, 3
         call ga_scale_divEpp(X_ga(iloop),dbl_mb(EIm_a))
      enddo
c
c     First terms with just V:
c     -1/2   A K pxVp K A   A     V    A E              from W^2E
c     -1/2   A     V    A   A K pxVp K A E              from W^2E
c     -1/2 E A K pxVp K A   A     V    A                from EW^2
c     -1/2 E A     V    A   A K pxVp K A                from EW^2
c     -1     A K pxVp K A E A     V    A                from WEW 
c     -1     A     V    A E A K pxVp K A                from WEW 
c
c     Get the V integrals, do A U^t Vdg_a U A and scale with 1/(Ep+Ep')
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA Vtype Z',1,1,V_ga))
     $     call errquit('int_1edk_ga: ga_create failed GA V_ga',0,
     &       GA_ERR)
      call ga_zero(V_ga)
      call int_1e_oldga(dk_fitbas,dk_fitbas,V_ga,'potential',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,V_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,V_ga)
      call ga_scale_divEpp(V_ga,dbl_mb(EIm_a))
      call ga_scale_rlh(V_ga,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
c     Loop over the three components
c     First pxVp.V, then v.pxVp, and terms with center E 
c
      do iloop = 1, 3
         call ga_zero(tmpg_a)
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(iloop),V_ga,0.0d0,
     $              tmpg_a)
         call ga_copy(tmpg_a,tempg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
         call ga_dadd(1.0d0,C_ga(iloop),-0.5d0,tmpg_a,C_ga(iloop))
         call ga_dadd(1.0d0,C_ga(iloop),-0.5d0,tempg_a,C_ga(iloop))
c
         call ga_zero(tmpg_a)
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,V_ga,X_ga(iloop),0.0d0,
     $                 tmpg_a)
         call ga_copy(tmpg_a,tempg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
         call ga_dadd(1.0d0,C_ga(iloop),-0.5d0,tmpg_a,C_ga(iloop))
         call ga_dadd(1.0d0,C_ga(iloop),-0.5d0,tempg_a,C_ga(iloop))
 
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(iloop),tmpg_a,
     $                 0.0d0,tempg_a)
         call ga_dadd(1.0d0,C_ga(iloop),-1.0d0,tempg_a,C_ga(iloop))
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,X_ga(iloop),
     $                 0.0d0,tempg_a)
         call ga_dadd(1.0d0,C_ga(iloop),-1.0d0,tempg_a,C_ga(iloop))
      enddo
c
c     Next terms with p.Vp:
c     +1/2   A K pxVp K A K^-2    p^-2 A K p.Vp K A E   from W^2E
c     +1/2   A K p.Vp K A K^-2    p^-2 A K pxVp K A E   from W^2E
c     +1/2 E A K pxVp K A K^-2    p^-2 A K p.Vp K A     from EW^2
c     +1/2 E A K p.Vp K A K^-2    p^-2 A K pxVp K A     from EW^2
c     +1     A K pxVp K A K^-2E(p)p^-2 A K p.Vp K A     from WEW 
c     +1     A K p.Vp K A K^-2E(p)p^-2 A K pxVp K A     from WEW 
c
c     Get the p.Vp integrals, do A K U^t Vdg_a U K A and scale with 1/(Ep+Ep')
c
      call ga_zero(V_ga)
      call int_1e_oldga(dk_fitbas,dk_fitbas,V_ga,'pVp',oskel)
      call ga_zero(tempg_a)
      call ga_dgemm('t','n',nbf,nbf,nbf,1.0d0,Ug_a,V_ga,0.0d0,tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,V_ga)
      call ga_scale_divEpp(V_ga,dbl_mb(EIm_a))
      call ga_scale_rlh(V_ga,dbl_mb(KIm_a),dbl_mb(KIm_a))
      call ga_scale_rlh(V_ga,dbl_mb(AIm_a),dbl_mb(AIm_a))
c
c     Loop over the three components
c     First pxVp.p.Vp, then p.Vp.pxVp, and terms with center E 
c
c     We need 1/(k^2p^2), put them in array first
c
      if (.not. ma_push_get(MT_DBL,nbf,'MA temp',tempm_a,tempIm_a))
     $     call errquit('D-K integrals: ma_push_get failed MA temp',0,
     &       MA_ERR)
      do 34 iloop = 1, nbf
         dbl_mb(tempIm_a+iloop-1) = 1.0d0/                 
     $                              dbl_mb(p2Im_a+iloop-1)/
     $                              dbl_mb(KIm_a+iloop-1)/
     $                              dbl_mb(KIm_a+iloop-1)
   34 continue
c
      do iloop = 1, 3
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(iloop),tmpg_a,
     $                 0.0d0,tempg_a)
         call ga_copy(tempg_a,tmpg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
         call ga_dadd(1.0d0,C_ga(iloop),+0.5d0,tmpg_a,C_ga(iloop))
         call ga_dadd(1.0d0,C_ga(iloop),+0.5d0,tempg_a,C_ga(iloop))
c
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(tempIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,X_ga(iloop),
     $                 0.0d0,tempg_a)
         call ga_copy(tempg_a,tmpg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
         call ga_dadd(1.0d0,C_ga(iloop),+0.5d0,tmpg_a,C_ga(iloop))
         call ga_dadd(1.0d0,C_ga(iloop),+0.5d0,tempg_a,C_ga(iloop))
c
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(iloop),tmpg_a,
     $                 0.0d0,tempg_a)
         call ga_dadd(1.0d0,C_ga(iloop),+1.0d0,tempg_a,C_ga(iloop))
         call ga_zero(tempg_a)
         call ga_copy(V_ga,tmpg_a)
         call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
         call ga_scale_rh(tmpg_a,dbl_mb(tempIm_a))
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,X_ga(iloop),
     $                 0.0d0,tempg_a)
         call ga_dadd(1.0d0,C_ga(iloop),+1.0d0,tempg_a,C_ga(iloop))
      enddo
c
c     At this point we probably have the set of integrals that would be similar 
c     to those for DKH, without the real pxVp terms, hence pxVp x pxVp.
c     If we neglect these remaining terms we get an incomplete Douglas-Kroll 
c     transformation a la Haberlen&Rosch (CPL, 199, 491, 1992) denoted vn2.
c     However, this splitting is not clear to me, especially because the previous 
c     pxVp terms seem to cancel each other and the final pxVp x pxVp term is 
c     significant for the core spin-orbit splitting.
c
c     Hence, following DKH type skipping of integrals will be taken out for now.
c     if (dk_form_type .eq. 2) goto 1001
c
c     Finally terms with pxVp x pxVp:
c     -1/2   A K pxVp K A K^-2    p^-2 A K pxVp K A E   from W^2E
c     -1/2 E A K pxVp K A K^-2    p^-2 A K pxVp K A     from EW^2
c     -1     A K pxVp K A K^-2E(p)p^-2 A K pxVp K A     from WEW 
c
      call ga_zero(tempg_a) 
      call ga_zero(tmpg_a)
      call ga_zero(V_ga)
      call ga_copy(X_ga(3),V_ga)
      call ga_scale_lh(V_ga,dbl_mb(tempIm_a))
      call ga_copy(X_ga(2),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(2),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(3),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_copy(tempg_a,tmpg_a)
      call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,C_ga(1),-0.5d0,tmpg_a,C_ga(1))
      call ga_dadd(1.0d0,C_ga(1),-0.5d0,tempg_a,C_ga(1))
      call ga_zero(tempg_a) 
      call ga_scale_lh(V_ga,dbl_mb(EIm_a))
      call ga_copy(X_ga(2),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(2),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(3),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_dadd(1.0d0,C_ga(1),-1.0d0,tempg_a,C_ga(1))
c
      call ga_zero(tempg_a) 
      call ga_zero(tmpg_a)
      call ga_zero(V_ga)
      call ga_copy(X_ga(1),V_ga)
      call ga_scale_lh(V_ga,dbl_mb(tempIm_a))
      call ga_copy(X_ga(3),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(3),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(1),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_copy(tempg_a,tmpg_a)
      call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,C_ga(2),-0.5d0,tmpg_a,C_ga(2))
      call ga_dadd(1.0d0,C_ga(2),-0.5d0,tempg_a,C_ga(2))
      call ga_zero(tempg_a) 
      call ga_scale_lh(V_ga,dbl_mb(EIm_a))
      call ga_copy(X_ga(3),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(3),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(1),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_dadd(1.0d0,C_ga(2),-1.0d0,tempg_a,C_ga(2))
c
      call ga_zero(tempg_a) 
      call ga_zero(tmpg_a)
      call ga_zero(V_ga)
      call ga_copy(X_ga(2),V_ga)
      call ga_scale_lh(V_ga,dbl_mb(tempIm_a))
      call ga_copy(X_ga(1),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(1),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(2),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_copy(tempg_a,tmpg_a)
      call ga_scale_lh(tempg_a,dbl_mb(EIm_a))
      call ga_scale_rh(tmpg_a,dbl_mb(EIm_a))
      call ga_dadd(1.0d0,C_ga(3),-0.5d0,tmpg_a,C_ga(3))
      call ga_dadd(1.0d0,C_ga(3),-0.5d0,tempg_a,C_ga(3))
      call ga_zero(tempg_a) 
      call ga_scale_lh(V_ga,dbl_mb(EIm_a))
      call ga_copy(X_ga(1),tmpg_a)
      call ga_scale_lh(tmpg_a,dbl_mb(tempIm_a))
      call ga_scale_lh(tmpg_a,dbl_mb(EIm_a))
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,X_ga(1),V_ga,0.0d0,
     $                 tempg_a)
      call ga_dgemm('n','n',nbf,nbf,nbf,-1.0d0,X_ga(2),tmpg_a,1.0d0,
     $                 tempg_a)
      call ga_dadd(1.0d0,C_ga(3),-1.0d0,tempg_a,C_ga(3))
c
 1001 if (.not. ma_pop_stack(tempm_a))
     $     call errquit('int_1edk_ga: ma_pop_stack tempm_a',0, MA_ERR)
      if (.not. ga_destroy(V_ga))
     $     call errquit('int_1edk_ga: ga_destroy V_ga',0, GA_ERR)
c
c     Get overlap and multiply with Ug_a to get SXU
c
 1002 call ga_zero(tempg_a)
      call ga_zero(tmpg_a)
      call int_1e_oldga(dk_fitbas, dk_fitbas, tempg_a,'overlap',.false.)
      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tempg_a,Ug_a,0.0d0,
     $              tmpg_a)
c
c     Transform final integrals back to real space integrals U V U^t
c
      do iloop = 1, 3
         call ga_zero(tempg_a)
         call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,tmpg_a,C_ga(iloop),
     $                 0.0d0,tempg_a)
         call ga_dgemm('n','t',nbf,nbf,nbf,1.0d0,tempg_a,tmpg_a,0.0d0,
     $                 C_ga(iloop))
      enddo
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
      do iloop = 1, 3
         if (.not. ga_destroy(X_ga(iloop)))
     $        call errquit('int_1edk_ga: ga_destroy X_ga',iloop, GA_ERR)
      enddo
c
c     We now have our integrals in a primitive basis in X_ga, Y_ga and Z_ga. If the original
c     basis is contracted we have to transform to this basis. Get transformation
c     matrix, combine with primitive integrals and put final integrals in g_in(1-3)
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
      do iloop = 1, 3
         call ga_zero(tempg_a)
         call ga_dgemm('t','n',nbf_in,nbf,nbf,1.0d0,Mixg_a,
     $                  C_ga(iloop),0.0d0,tempg_a)
          call ga_dgemm('n','n',nbf_in,nbf_in,nbf,1.0d0,tempg_a,Mixg_a,
     $                  1.0d0,g_in(3-iloop+1))
      enddo
c
c     Destroy transformation matrices Mixg_a, tempg_a and DKg_a
c
      if (.not. ga_destroy(tempg_a))
     $     call errquit('int_1edk_ga: ga_destroy tempg_a',0, GA_ERR)
      if (.not. ga_destroy(Mixg_a))
     $     call errquit('int_1edk_ga: ga_destroy Mixg_a',0, GA_ERR)
      do iloop = 1, 3
         if (.not. ga_destroy(C_ga(iloop)))
     $        call errquit('int_1edk_ga: ga_destroy C_ga',iloop, GA_ERR)
      enddo
c
      return
c
      end
C>
C> @}
