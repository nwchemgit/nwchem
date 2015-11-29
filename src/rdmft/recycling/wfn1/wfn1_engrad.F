C> \ingroup wfn1
C> @{
C>
C> \brief Evaluate the energies and optionally gradients for a number
C> of states
C>
C> Given the current values of the occupation functions, natural
C> orbitals and the Lagrangians of a number of "states" evaluate
C> the energy of each state. Here we define a state as a set of values
C> of the variables that belong together. Different sets may represent
C> different electronic states, or different points on the same
C> electronic state. Nevertheless, each set defines a different 
C> independent energy expression we refer to as a state. 
C>
C> Optionally, this routine also evaluates the gradients of these
C> states. However, for most optimization schemes it is expected that
C> an energy and a gradient will be requested for a single state or
C> that only the energy is requested for multiple states.
C>
      subroutine wfn1_engrad(rtdb,geom,nst,nbf,nea,neb,h1,eri,erix,ov,
     +                       sfac,ofac,r_noa,r_nob,r_ofa,r_ofb,
     +                       ehfa,ehfb,ewfa,ewfb,enta,entb,temperature,
     +                       power,nperma,npermb,perma,permb,
     +                       signa,signb,ovla,ovlb,
     +                       En,serr,oerr,oprint)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer rtdb !< [Input] The RTDB handle
      integer geom !< [Input] The geometry handle
      integer nbf  !< [Input] The number of basis functions
      integer nst  !< [Input] The number of states for which
                   !< energies and gradients should be evaluated
      integer nea  !< [Input] The number of \f$\alpha\f$-spin
                   !< electrons for each state
      integer neb  !< [Input] The number of \f$\beta\f$-spin
                   !< electrons for each state
      integer nperma !< [Input] The number alpha electron permutations
      integer npermb !< [Input] The number beta electron permutations
      double precision h1(nbf,nbf) !< [Input] The 1-electron integrals
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                            !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                             !< integrals for exchange
      double precision ov(nbf,nbf) !< [Input] The overlap integrals
      double precision sfac !< [Input] The scale factor for the spin
                            !< penality function
      double precision ofac !< [Input] The scale factor for the orbital
                            !< order penality function
      double precision r_noa(nbf,nbf,nst) !< [Input] The 
      !< \f$\alpha\f$-spin natural orbitals. This is a 3D array with the
      !< dimensions:
      !< - 1. The basis functions
      !< - 2. The natural orbitals
      !< - 3. The states
      double precision r_nob(nbf,nbf,nst) !< [Input] The 
      !< \f$\beta\f$-spin natural orbitals. This is a 3D array with the
      !< same dimensions as `r_noa`.
      double precision r_ofa(nbf,nbf,nst) !< [Input] The
      !< \f$\alpha\f$-spin occupation functions. This is a 3D array with
      !< the dimension:
      !< - 1. The natural orbitals
      !< - 2. The occupation functions
      !< - 3. The states
      double precision r_ofb(nbf,nbf,nst) !< [Input] The 
      !< \f$\beta\f$-spin occupation functions. This is a 3D array with
      !< the same dimensions as `r_ofa`.
c
      double precision ehfa(nbf) !< [Input] The Hartree-Fock alpha
      !< orbital energies.
      double precision ehfb(nbf) !< [Input] The Hartree-Fock beta
      !< orbital energies.
      double precision ewfa(nbf) !< [Output] The WFN1 alpha
      !< orbital energies.
      double precision ewfb(nbf) !< [Output] The WFN1 beta
      !< orbital energies.
      double precision enta(nbf) !< [Output] The alpha orbital
      !< entropies.
      double precision entb(nbf) !< [Output] The beta orbital entropies.
      double precision temperature !< [Input] The electron temperature
      double precision power !< [Input] The matrix power
c
      integer perma(nea,nperma) !< [Input] The alpha electron
                                !< permutations
      integer permb(neb,npermb) !< [Input] The beta electron
                                !< permutations
      double precision signa(nperma) !< [Input] The alpha permutation
                                     !< signs
      double precision signb(npermb) !< [Input] The beta permutation
                                     !< signs
      double precision ovla(nea,nea) !< [Scratch] The alpha overlaps
      double precision ovlb(neb,neb) !< [Scratch] The beta overlaps
c
      double precision En !< [Output] The total energies of all 
                          !< states
      double precision serr !< [Output] The \f$ S^2 \f$ error
      double precision oerr(nst) !< [Output] The orbital order error
c
      logical oprint !< [Input] If .true. print some energy data
c
C     Local Variables
C
      double precision E1el
      double precision E2elc
      double precision E2elcsr
      double precision E2elx
      double precision E2elxsr
      double precision ES2err
      double precision Eorder
      double precision Exc
      double precision Enuc
      double precision Sa,    Sb
      double precision Sa_p,  Sb_p
      double precision Sa_b,  Sb_b
      double precision Sa_e,  Sb_e
      double precision TSa, TSb
      double precision dSa, dSb, dSt
      double precision s2t, s2a
      double precision emin
cDEBUG
c     double precision E1enta(nbf), E1entb(nbf)
c     double precision E2Centa(nbf), E2Centb(nbf)
c     double precision E2Xenta(nbf), E2Xentb(nbf)
cDEBUG
c
      integer net
      integer len
      integer l_m1, k_m1
      integer l_m2, k_m2
      integer l_m3, k_m3
      integer l_m4, k_m4
      integer l_m5, k_m5
      integer l_v1, k_v1
      integer l_v2, k_v2
      integer l_v3, k_v3
      integer l_a , k_a 
      integer l_b , k_b
      integer ii, jj, imin, ist, ipar1, ipar2
C
C     Functions
C
      logical          geom_nuc_rep_energy
      external         geom_nuc_rep_energy
      double precision wfn1_s2, wfn1_overlap_exp, wfn1_overlap_p
      external         wfn1_s2, wfn1_overlap_exp, wfn1_overlap_p
      double precision wfn1_overlap_bo
      external         wfn1_overlap_bo
      double precision xc_hfexch !< The exact exchange weight in DFT
      external         xc_hfexch 
cDEBUG
      logical stat
cDEBUG
C
C     Code
C
      En      = 0.0d0
      E1el    = 0.0d0
      E2elc   = 0.0d0
      E2elcsr = 0.0d0
      E2elx   = 0.0d0
      E2elxsr = 0.0d0
      Eorder  = 0.0d0
cDEBUG
c     E1enta   = 0.0d0
c     E1entb   = 0.0d0
c     E2Centa  = 0.0d0
c     E2Centb  = 0.0d0
c     E2Xenta  = 0.0d0
c     E2Xentb  = 0.0d0
cDEBUG
      if (.not.geom_nuc_rep_energy(geom,Enuc))
     &  call errquit("wfn1_energy: Enuc failed",0,GEOM_ERR)
c
      len = nbf*nbf
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"M1",l_m1,k_m1))
     +   call errquit("wfn1_engrad: could not allocate M1",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"M2",l_m2,k_m2))
     +   call errquit("wfn1_engrad: could not allocate M2",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"M3",l_m3,k_m3))
     +   call errquit("wfn1_engrad: could not allocate M3",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"M4",l_m4,k_m4))
     +   call errquit("wfn1_engrad: could not allocate M4",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"M5",l_m5,k_m5))
     +   call errquit("wfn1_engrad: could not allocate M5",0,MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"V1",l_v1,k_v1))
     +   call errquit("wfn1_engrad: could not allocate V1",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"V2",l_v2,k_v2))
     +   call errquit("wfn1_engrad: could not allocate V2",0,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"V3",l_v3,k_v3))
     +   call errquit("wfn1_engrad: could not allocate V3",0,MA_ERR)
cDEBUG
c     if (.not.ma_push_get(MT_DBL,nbf,"VA",l_a,k_a))
c    +   call errquit("wfn1_engrad: could not allocate VA",0,MA_ERR)
c     if (.not.ma_push_get(MT_DBL,nbf,"VB",l_b,k_b))
c    +   call errquit("wfn1_engrad: could not allocate VB",0,MA_ERR)
cDEBUG
c
      s2t = (0.5d0*(nea-neb))*(0.5d0*(nea-neb)+1.0d0)
      s2a = wfn1_s2(nbf,nea,neb,r_ofa,r_ofb,r_noa,r_nob,ov,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_v1))
cDEBUG
c     write(*,*)'*** HVD A: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      serr = (s2a-s2t)**1
c     ES2err = sfac*(serr*serr)
      ES2err = sfac*serr
c
      call wfn1_1el(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +              r_ofa(1,1,nst),r_ofb(1,1,nst),
     +              h1,E1el,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD B: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_c(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +              r_ofa(1,1,nst),r_ofb(1,1,nst),
     +              eri,E2elc,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD C: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_hfx(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +              r_ofa(1,1,nst),r_ofb(1,1,nst),ov,
     +              erix,E2elx,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2),power)
cDEBUG
c     write(*,*)'*** HVD D: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_csr(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +              r_ofa(1,1,nst),r_ofb(1,1,nst),
     +              eri,E2elcsr,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD E: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_csr(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +              r_ofa(1,1,nst),r_ofb(1,1,nst),
     +              erix,E2elxsr,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD F: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_entropy(.true.,nbf,nea,r_ofa(1,1,nst),enta,
     +              dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD G: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_entropy(.true.,nbf,neb,r_ofb(1,1,nst),entb,
     +              dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD H: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_xc(rtdb,nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +             r_ofa(1,1,nst),r_ofb(1,1,nst),Exc,dbl_mb(k_m1),
     +             dbl_mb(k_m2))
cDEBUG
c     write(*,*)'*** HVD I: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
c
c     Calculate the orbital energies
c
      call dfill(nbf,0.0d0,ehfa,1)
      call dfill(nbf,0.0d0,ehfb,1)
      call dfill(nbf,0.0d0,ewfa,1)
      call dfill(nbf,0.0d0,ewfb,1)
      call wfn1_1el_orb(nbf,nea,r_noa(1,1,nst),r_ofa(1,1,nst),h1,ehfa,
     +                  dbl_mb(k_m1),dbl_mb(k_v1))
cDEBUG
c     write(*,*)'*** HVD J: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_1el_orb(nbf,neb,r_nob(1,1,nst),r_ofb(1,1,nst),h1,ehfb,
     +                  dbl_mb(k_m1),dbl_mb(k_v1))
cDEBUG
c     write(*,*)'*** HVD K: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_C_orb(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +                  r_ofa(1,1,nst),r_ofb(1,1,nst),eri,
     +                  ewfa,ewfb,
     +                  dbl_mb(k_m1),dbl_mb(k_m2),
     +                  dbl_mb(k_m3),dbl_mb(k_m4),
     +                  dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD L: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      call wfn1_2el_hfX_orb(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
     +                  r_ofa(1,1,nst),r_ofb(1,1,nst),erix,
     +                  ewfa,ewfb,
     +                  dbl_mb(k_m1),dbl_mb(k_m2),
     +                  dbl_mb(k_m3),dbl_mb(k_m4),
     +                  dbl_mb(k_v1),dbl_mb(k_v2))
cDEBUG
c     write(*,*)'*** HVD M: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      do ii = 1, nbf
        ehfa(ii) = ehfa(ii) + ewfa(ii)
        ehfb(ii) = ehfb(ii) + ewfb(ii)
      enddo
c     call wfn1_2el_Cslf_orb(nbf,nea,neb,r_noa(1,1,nst),r_nob(1,1,nst),
c    +                  r_ofa(1,1,nst),r_ofb(1,1,nst),eri,
c    +                  ewfa,ewfb,
c    +                  dbl_mb(k_m1),dbl_mb(k_m2),
c    +                  dbl_mb(k_m3),dbl_mb(k_m4),
c    +                  dbl_mb(k_v1),dbl_mb(k_v2))
c     call wfn1_2el_hfXslf_orb(nbf,nea,neb,
c    +                  r_noa(1,1,nst),r_nob(1,1,nst),
c    +                  r_ofa(1,1,nst),r_ofb(1,1,nst),erix,
c    +                  ewfa,ewfb,
c    +                  dbl_mb(k_m1),dbl_mb(k_m2),
c    +                  dbl_mb(k_m3),dbl_mb(k_m4),
c    +                  dbl_mb(k_v1),dbl_mb(k_v2))
c
      E2elcsr = -E2elcsr
      E2elx = E2elx*xc_hfexch()
c     En = Enuc+E1el+E2elc+E2elx+E2elcsr+E2elxsr
      En = Enuc+E1el+E2elc+E2elx+ES2err
      En = Enuc+E1el+E2elc+E2elx+Exc
      dSa = 0.0d0
      dSb = 0.0d0
      do ii = 1, nbf
        dSa = dSa + enta(ii)
        dSb = dSb + entb(ii)
c       ewfa(ii) = ewfa(ii)/9.48d0
c       ewfb(ii) = ewfb(ii)/9.48d0
cInf    ewfa(ii) = ewfa(ii)/30.0d0 ! Used with entropy separation 
cInf    ewfb(ii) = ewfb(ii)/30.0d0 ! approach. Too small for average
                                   ! entropy approach
cOne    ewfa(ii) = ewfa(ii)/26.0d0 ! Still too large a factor but the
cOne    ewfb(ii) = ewfb(ii)/26.0d0 ! electrons get smeared too much
                                   ! already. The infinite degeneracy
                                   ! entropy model is clearly better.
cHACK
        ewfa(ii) = ewfa(ii)/100.0d0
        ewfb(ii) = ewfb(ii)/100.0d0
      enddo
      dSt = dSa+dSb
      net = nea+neb
      TSa = 0.0d0
      TSb = 0.0d0
c     dSa = dSa/3.0d0
c     dSb = dSb/3.0d0
      do ii = 1, nbf
        if (oprint) then
c         write(*,'(" alpha ent: ",i4,4f14.8)')ii,
c    +      (ewfa(ii)-ehfa(ii))*enta(ii)/(3.0d0/2.0d0),
c    +      (ewfa(ii)-ehfa(ii))*enta(ii),
c    +      enta(ii),ehfa(ii),ewfa(ii)
c         write(*,'(" alpha ent: ",i4,3f14.8)')ii,
c    +      enta(ii),ehfa(ii),ewfa(ii)
c         write(*,'(" alpha ent: ",i4,5f14.8)')ii,
c    +      enta(ii),ewfa(ii),enta(ii)*ewfa(ii),ehfa(ii),
c    +                        ehfa(ii)+enta(ii)*ewfa(ii)
          if (nea.gt.0) then
c           write(*,'(" alpha ent: ",i4,4f14.8)')ii,
c    +        enta(ii),ewfa(ii),ewfa(ii)*dSa/nea,
c    +                          ehfa(ii)+ewfa(ii)*dSa/nea
            write(*,'(" alpha ent: ",i4,4f14.8)')ii,
     +        enta(ii),ewfa(ii),ewfa(ii)*dSt/net,
     +                          ehfa(ii)+ewfa(ii)*dSt/net
          else
            write(*,'(" alpha ent: ",i4,4f14.8)')ii,
     +        enta(ii),ewfa(ii),0.0d0,
     +                          ehfa(ii)+0.0d0
          endif
        endif
      enddo
      if (oprint) then
        write(*,*)
      endif
      do ii = 1, nbf
c       TSa = TSa + (ewfa(ii)-ehfa(ii))*enta(ii)*(3.0d0/2.0d0)
c       TSa = TSa + max(0.0d0,temperature)*enta(ii)
c       TSa = TSa + ewfa(ii)*enta(ii)
c       TSa = TSa + ewfa(ii)*(dSa+dSb)
c       TSa = TSa + ewfa(ii)*(dSa+dSb)/(nea+neb)
        if (nea.gt.0) then
c         TSa = TSa + ewfa(ii)*dSa/nea
          TSa = TSa + ewfa(ii)*dSt/net
        endif
      enddo
      do ii = 1, nbf
        if (oprint) then
c         write(*,'(" beta  ent: ",i4,4f14.8)')ii,
c    +      (ewfb(ii)-ehfb(ii))*entb(ii)/(3.0d0/2.0d0),
c    +      (ewfb(ii)-ehfb(ii))*entb(ii),
c    +      entb(ii),ehfb(ii),ewfb(ii)
c         write(*,'(" beta  ent: ",i4,3f14.8)')ii,
c    +      entb(ii),ehfb(ii),ewfb(ii)
c         write(*,'(" beta  ent: ",i4,5f14.8)')ii,
c    +      entb(ii),ewfb(ii),entb(ii)*ewfb(ii),ehfb(ii),
c    +                        ehfb(ii)+entb(ii)*ewfb(ii)
          if (neb.gt.0) then
c           write(*,'(" beta  ent: ",i4,4f14.8)')ii,
c    +        entb(ii),ewfb(ii),ewfb(ii)*dSb/neb,
c    +                          ehfb(ii)+ewfb(ii)*dSb/neb
            write(*,'(" beta  ent: ",i4,4f14.8)')ii,
     +        entb(ii),ewfb(ii),ewfb(ii)*dSt/net,
     +                          ehfb(ii)+ewfb(ii)*dSt/net
          else
            write(*,'(" beta  ent: ",i4,4f14.8)')ii,
     +        entb(ii),ewfb(ii),0.0d0,
     +                          ehfb(ii)+0.0d0
          endif
        endif
      enddo
      if (oprint) then
        write(*,*)
      endif
      do ii = 1, nbf
c       TSb = TSb + (ewfb(ii)-ehfb(ii))*entb(ii)*(3.0d0/2.0d0)
c       TSb = TSb + max(0.0d0,temperature)*entb(ii)
c       TSb = TSb + ewfb(ii)*entb(ii)
c       TSb = TSb + ewfb(ii)*(dSa+dSb)
c       TSb = TSb + ewfb(ii)*(dSa+dSb)/(nea+neb)
        if (neb.gt.0) then
c         TSb = TSb + ewfb(ii)*dSb/neb
          TSb = TSb + ewfb(ii)*dSt/net
        endif
      enddo
      if (temperature.gt.0.0d0) then
        En = En + TSa + TSb
      endif
c
c     Calculate the overlap penalty
c
      do ist = 1, nst
        oerr(ist) = 0.0d0
      enddo
      do ist = 1, nst-1
        call wfn1_ovlp_2_states(nea,nbf,r_ofa(1,1,ist),r_noa(1,1,ist),
     +                          r_ofa(1,1,nst),r_noa(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovla)
cDEBUG
c     write(*,*)'*** HVD N: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
        call wfn1_ovlp_2_states(neb,nbf,r_ofb(1,1,ist),r_nob(1,1,ist),
     +                          r_ofb(1,1,nst),r_nob(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovlb)
cDEBUG
c     write(*,*)'*** HVD O: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
c       sa = abs(wfn1_overlap_exp(nea,nperma,ovla,perma,signa))
c       sb = abs(wfn1_overlap_exp(neb,npermb,ovlb,permb,signb))
c       Sa_e = wfn1_overlap_exp(nea,nperma,ovla,perma,signa)
c       Sb_e = wfn1_overlap_exp(neb,npermb,ovlb,permb,signb)
c       Sa_p = wfn1_overlap_p(nea,ovla,dbl_mb(k_m1),dbl_mb(k_m2))
c       Sb_p = wfn1_overlap_p(neb,ovlb,dbl_mb(k_m1),dbl_mb(k_m2))
        Sa_e = wfn1_overlap_bo(nea,ovla,dbl_mb(k_m1),dbl_mb(k_m2))
        Sb_e = wfn1_overlap_bo(neb,ovlb,dbl_mb(k_m1),dbl_mb(k_m2))
cDEBUG
c     write(*,*)'*** HVD P: wfn1_engrad'
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
cDEBUG
c       write(*,*)'*** wfn1_overlap_ldr_tst: alpha'
cDEBUG
c       call wfn1_overlap_ldr_tst(nea,ovla,dbl_mb(k_m1),dbl_mb(k_m2),
c    +                            dbl_mb(k_m3),dbl_mb(k_m4),
c    +                            dbl_mb(k_v1),ipar1,dbl_mb(k_m5))
cDEBUG
c       write(*,*)'*** wfn1_overlap_ldr_tst: beta'
cDEBUG
c       call wfn1_overlap_ldr_tst(neb,ovlb,dbl_mb(k_m1),dbl_mb(k_m2),
c    +                            dbl_mb(k_m3),dbl_mb(k_m4),
c    +                            dbl_mb(k_v2),ipar2,dbl_mb(k_m5))
cDEBUG
c       write(*,*)'*** wfn1_overlap_ldr_tst: done'
cDEBUG
        Sa_b = 1.0d0*ipar1
        Sb_b = 1.0d0*ipar2
        do ii = 1, nea
          Sa_b = Sa_b*dbl_mb(k_v1+ii-1)
        enddo
        do ii = 1, neb
          Sb_b = Sb_b*dbl_mb(k_v2+ii-1)
        enddo
        
        Sa   = sa_e
        Sb   = sb_e
c       if (abs(sa-sa_p).gt.1.0d-7) then
c         write(*,*)ist,nst,':engrad: sa,sa_p = ',sa,sa_p
c       endif
c       if (abs(sb-sb_p).gt.1.0d-7) then
c         write(*,*)ist,nst,':engrad: sb,sb_p = ',sb,sb_p
c       endif
c       if (abs(sa-sa_b).gt.1.0d-7) then
c         write(*,*)ist,nst,':engrad: sa,sa_b = ',sa,sa_b
c       endif
c       if (abs(sb-sb_b).gt.1.0d-7) then
c         write(*,*)ist,nst,':engrad: sb,sb_b = ',sb,sb_b
c       endif
cDEBUG
c       if (oprint) then
c         write(*,*)'*** wfn1_engrad: sa =',ist,sa
c         call hess_hssout(ovla,nea,nea,nea)
c         write(*,*)'*** wfn1_engrad: sb =',ist,sb
c         call hess_hssout(ovlb,neb,neb,neb)
c         call hess_hssout(dbl_mb(k_m1),neb,neb,neb)
c       endif
cDEBUG
c       en = en + (sa*sb*1.0d0)**2
        oerr(ist) = oerr(ist) + (sa*sb)**2
c       oerr(ist) = oerr(ist) + (sa*sb)
      enddo
c
c     call wfn1_order(nbf,nea,neb,ehfa,ehfb,ewfa,ewfb,ofac,
c    +                Eorder,oerr)
c     En = En + Eorder
c
      if (oprint) then
        call wfn1_ovlp_2_states(nea,nbf,r_ofa(1,1,nst),r_noa(1,1,nst),
     +                          r_ofa(1,1,nst),r_noa(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovla)
        call wfn1_ovlp_2_states(neb,nbf,r_ofb(1,1,nst),r_nob(1,1,nst),
     +                          r_ofb(1,1,nst),r_nob(1,1,nst),ov,
     +                          dbl_mb(k_m1),dbl_mb(k_m2),ovlb)
        sa = wfn1_overlap_exp(nea,nperma,ovla,perma,signa)
        sb = wfn1_overlap_exp(neb,npermb,ovlb,permb,signb)
        write(*,*)'The 1-electron energy             = ',E1el
        write(*,*)'The Coulomb energy                = ',E2elc
        write(*,*)'The HF exchange energy            = ',E2elx
        write(*,*)'The DFT xc energy                 = ',Exc
c       write(*,*)'The spin penalty energy           = ',ES2err,s2t,s2a
        write(*,*)'The spin penalty energy           = ',ES2err,sfac
        write(*,*)'The alpha entropy                 = ',TSa
        write(*,*)'The beta entropy                  = ',TSb
        write(*,*)'The Coulomb self repulsion energy = ',E2elcsr
        write(*,*)'The self exchange energy          = ',E2elxsr
        write(*,*)'The wavefunction overlap          = ',oerr
        write(*,*)'The wavefunction overlap          = ',sa,sb
        write(*,*)'The total wfn1 energy             = ',En
        write(*,*)
      endif
c
      if (.not.ma_pop_stack(l_v3))
     +   call errquit("wfn1_engrad: could not deallocate V3",0,MA_ERR)
      if (.not.ma_pop_stack(l_v2))
     +   call errquit("wfn1_engrad: could not deallocate V2",0,MA_ERR)
      if (.not.ma_pop_stack(l_v1))
     +   call errquit("wfn1_engrad: could not deallocate V1",0,MA_ERR)
c
      if (.not.ma_pop_stack(l_m5))
     +   call errquit("wfn1_engrad: could not deallocate M5",0,MA_ERR)
      if (.not.ma_pop_stack(l_m4))
     +   call errquit("wfn1_engrad: could not deallocate M4",0,MA_ERR)
      if (.not.ma_pop_stack(l_m3))
     +   call errquit("wfn1_engrad: could not deallocate M3",0,MA_ERR)
      if (.not.ma_pop_stack(l_m2))
     +   call errquit("wfn1_engrad: could not deallocate M2",0,MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     +   call errquit("wfn1_engrad: could not deallocate M1",0,MA_ERR)
c
      end
C> @}

