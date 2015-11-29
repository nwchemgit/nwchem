C> \ingroup wfn1
C> @{
C>
C> \brief Orthogonalize \f$\Psi_I\f$ on the set of orthogonal states
C> \f$\Psi_J,\; J=1,I-1\f$
C>
C> Assume we have an orthonormal set of N-electron wavefunctions 
C> \f$\left\{\Psi_J\right\}\f$ where \f$1\le J \lt I\f$. 
C> Assume also that \f$\Psi_I\f$ is normalized. We want to
C> orthonormalize \f$\Psi_I\f$ onto that set. Because of the way the
C> superposition principle works for WFN1 style wavefunctions it is
C> most efficient to use Gramm-Schmidt (even though that is less 
C> accurate than modified Gramm-Schmidt). 
C>
C> The approach then is to:
C>
C> * Calculate \f$\forall_{J<I}, c_{IJ} = \left\langle\Psi_I\right|\left.\Psi_J\right\rangle\f$.
C>   The new wavefunction then is \f$\Psi'_I = A\left(\Psi_I-\sum_{J=1}^{I-1}c_{IJ}\Psi_J\right)\f$
C>   where \f$A=1/\sqrt{1+\sum_J(c_{IJ}^2-2c_{IJ})}\f$.
C>
C> * Calculate the 1-electron density matrix
C>   \f$D = \left|\Psi'_I\right\rangle\left\langle\Psi'_I\right|\f$.
C>
C> * Diagonalize \f$D\f$ to obtain new natural orbitals 
C>   \f$\left\{^I\phi^{D\sigma}\right\}\f$ and associated occupation
C>   numbers.
C>
C> * Update the occupation functions to match the new occupation
C>   numbers.
C>
C> In particular the requirement to construct and diagonalize a 
C> density matrix makes the simultaneous orthonormalization on all
C> lower states attractive (as opposed to a modified Gramm-Schmidt).
C>
      subroutine wfn1_orthogonalize(nst,nbf,nea,neb,ov,r_noa,r_nob,
     &                              r_ofa,r_ofb,ovla,ovlb,cij,
     &                              docca,doccb,da,db)
      implicit none
c
      integer nst !< [Input] The number of states (\f$=I\f$)
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$-electrons
      integer neb !< [Input] The number of \f$\beta\f$-electrons
c
      double precision ov(nbf,nbf) !< [Input] The overlap matrix
      double precision r_noa(nbf,nbf,nst) !< [In/Output] The
      !< \f$\alpha\f$ natural orbitals
      double precision r_nob(nbf,nbf,nst) !< [In/Output] The
      !< \f$\beta\f$ natural orbitals
      double precision r_ofa(nbf,nbf,nst) !< [In/Output] The
      !< \f$\alpha\f$ occupation functions
      double precision r_ofb(nbf,nbf,nst) !< [In/Output] The
      !< \f$\beta\f$ occupation functions
      double precision ovla(nea,nea) !< [Scratch] The \f$\alpha\f$
      !< overlaps
      double precision ovlb(neb,neb) !< [Scratch] The \f$\beta\f$
      !< overlaps
      double precision cij(nst) !< [Scratch] The state overlaps
      double precision docca(nbf) !< [Scratch] The \f$\alpha\f$
      !< occupation numbers
      double precision doccb(nbf) !< [Scratch] The \f$\beta\f$
      !< occupation numbers
      double precision da(nbf,nbf) !< [Scratch] The \f$\alpha\f$
      !< density matrix
      double precision db(nbf,nbf) !< [Scratch] The \f$\beta\f$
      !< density matrix
c
      integer ii !< Counter over natural orbitals
      integer ir !< Counter over occupation functions
      integer ia !< Counter over basis functions
      integer ib !< Counter over basis functions
c
      integer GENEIG !< Parameter for Generalized Eigenvalue Problem
      parameter(GENEIG = 1)
c
c     First query LAPACK DSYGV for memory requirements
c
      lwork = -1
      info  =  0
      call dsygv(GENEIG,'V','U',nbf,da,nbf,db,nbf,docca,ww,lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_orthogonalize: a: dsygv failed ",info,UERR)
      endif
      lwork = ww+0.2d0
      if (.not.ma_push_get(MT_DBL,lwork,"work dsygv",l_w,k_w)) then
        call errquit("wfn1_orthogonalize: alloc work failed",
     +               ma_sizeof(MT_DBL,lwork,MT_BYTE),MA_ERR)
      endif
c
c     Calculate the overlaps and the normalization factor
c
      cij(nst) = 1.0d0
      do jst = 1, nst-1
        call wfn1_ovlp_2_states(nea,nbf,r_ofa(1,1,jst),r_noa(1,1,jst),
     +                          r_ofa(1,1,nst),r_noa(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovla)
        call wfn1_ovlp_2_states(neb,nbf,r_ofb(1,1,jst),r_nob(1,1,jst),
     +                          r_ofb(1,1,nst),r_nob(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovlb)
        Sa_b = wfn1_overlap_bo(nea,ovla,dbl_mb(k_m1),dbl_mb(k_m2))
        Sb_b = wfn1_overlap_bo(neb,ovlb,dbl_mb(k_m1),dbl_mb(k_m2))
        cij(jst) = Sa_b*Sb_b
      enddo
      aa = 1.0d0
      do jst = 1, nst-1
        aa = aa + cij(jst)**2-2.0d0*cij(jst)
      enddo
c
c     Construct the density matrix of the orthonormalize state
c
c     Do diagonal terms
c
      call wfn1_add_dmat(nbf,nea,neb,1.0d0,
     +                   r_noa(1,1,nst),r_nob(1,1,nst),
     +                   r_ofa(1,1,nst),r_ofb(1,1,nst),
     +                   0.0d0,da,db,docca,doccb)
      do jst = 1, nst-1
        call wfn1_add_dmat(nbf,nea,neb,cij(jst)**2,
     +                     r_noa(1,1,jst),r_nob(1,1,jst),
     +                     r_ofa(1,1,jst),r_ofb(1,1,jst),
     +                     1.0d0,da,db,docca,doccb)
      enddo
c
c     Do cross terms
c
      do jst = 1, nst-1
c
c       Do <Y_J|Y_N>
c
        call wfn1_ovlp_2_states(nea,nbf,r_ofa(1,1,jst),r_noa(1,1,jst),
     +                          r_ofa(1,1,nst),r_noa(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovla)
        call wfn1_ovlp_2_states(neb,nbf,r_ofb(1,1,jst),r_nob(1,1,jst),
     +                          r_ofb(1,1,nst),r_nob(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovlb)
        call wfn1_overlap_ldr(nea,ovla,lpa,rpa,lma,rma,va,ipa,tmp)
        call wfn1_overlap_ldr(neb,ovlb,lpb,rpb,lmb,rmb,vb,ipb,tmp)
        call wfn1_LR2d(nbf,nea,lpa,rpa,lma,rma,
     +                 r_ofa(1,1,jst),r_ofa(1,1,nst),
     +                 r_ofa2(1,1,jst),r_ofa2(1,1,nst),tmp)
        call wfn1_LR2d(nbf,neb,lpb,rpb,lmb,rmb,
     +                 r_ofb(1,1,jst),r_ofb(1,1,nst),
     +                 r_ofb2(1,1,jst),r_ofb2(1,1,nst),tmp)
        call wfn1_add_tdmat(nbf,nea,neb,-cij(jst)*ipa*ipb,va,vb,
     +                 r_noa(1,1,jst),r_nob(1,1,jst),
     +                 r_ofa2(1,1,jst),r_ofb2(1,1,jst),
     +                 r_noa(1,1,nst),r_nob(1,1,nst),
     +                 r_ofa2(1,1,nst),r_ofb2(1,1,nst),1.0d0,
     +                 da,db,docca,doccb)
c
c       Do <Y_N|Y_J>
c
        call wfn1_ovlp_2_states(nea,nbf,r_ofa(1,1,nst),r_noa(1,1,nst),
     +                          r_ofa(1,1,jst),r_noa(1,1,jst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovla)
        call wfn1_ovlp_2_states(neb,nbf,r_ofb(1,1,nst),r_nob(1,1,nst),
     +                          r_ofb(1,1,jst),r_nob(1,1,jst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovlb)
        call wfn1_overlap_ldr(nea,ovla,lpa,rpa,lma,rma,va,ipa,tmp)
        call wfn1_overlap_ldr(neb,ovlb,lpb,rpb,lmb,rmb,vb,ipb,tmp)
        call wfn1_LR2d(nbf,nea,lpa,rpa,lma,rma,
     +                 r_ofa(1,1,nst),r_ofa(1,1,jst),
     +                 r_ofa2(1,1,nst),r_ofa2(1,1,jst),tmp)
        call wfn1_LR2d(nbf,neb,lpb,rpb,lmb,rmb,
     +                 r_ofb(1,1,nst),r_ofb(1,1,jst),
     +                 r_ofb2(1,1,nst),r_ofb2(1,1,jst),tmp)
        call wfn1_add_tdmat(nbf,nea,neb,-cij(jst)*ipa*ipb,va,vb,
     +                 r_noa(1,1,nst),r_nob(1,1,nst),
     +                 r_ofa2(1,1,nst),r_ofb2(1,1,nst),
     +                 r_noa(1,1,jst),r_nob(1,1,jst),
     +                 r_ofa2(1,1,jst),r_ofb2(1,1,jst),1.0d0,
     +                 da,db,docca,doccb)
c
      enddo
c
c     Normalize the new density matrices
c
      call dscal(nbf*nbf,1.0d0/aa,da,1)
      call dscal(nbf*nbf,1.0d0/aa,db,1)
c
c     Calculate new natural orbitals
c
      call dcopy(nbf*nbf,da,1,r_noa(1,1,nst),1)
      call dcopy(nbf*nbf,ov,1,da,1)
      call dsygv(GENEIG,'V','U',nbf,r_noa(1,1,nst),nbf,da,nbf,docca,
     +           dbl_mb(k_w),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_orthogonalize: b: dsygv failed ",info,UERR)
      endif
      call dcopy(nbf*nbf,db,1,r_nob(1,1,nst),1)
      call dcopy(nbf*nbf,ov,1,db,1)
      call dsygv(GENEIG,'V','U',nbf,r_nob(1,1,nst),nbf,db,nbf,doccb,
     +           dbl_mb(k_w),lwork,info)
      if (info.ne.0) then
        call errquit("wfn1_orthogonalize: c: dsygv failed ",info,UERR)
      endif
c
c     Update the occupation to re-establish the wavefunction
c
      call wfn1_match_occ(.false.,nbf,nea,docca,r_ofa(1,1,nst),m1,m2,
     +                    m3,m4,v1,.false.)
      call wfn1_match_occ(.false.,nbf,neb,doccb,r_ofb(1,1,nst),m1,m2,
     +                    m3,m4,v1,.false.)
c
      end
C>
C> @}
