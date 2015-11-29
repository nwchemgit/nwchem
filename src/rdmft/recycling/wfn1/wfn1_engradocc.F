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
      subroutine wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                          factor,ofac,r_noa,r_nob,r_ofa,r_ofb,
     +                          ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                          occa,occb,En,deviation,oerr,oprint)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer rtdb !< [Input] The RTDB handle
      integer nbf !< [Input] The number of basis functions
      integer nst !< [Input] The number of states for which
                  !< energies and gradients should be evaluated
      integer nea !< [Input] The number of \f$\alpha\f$-spin
                  !< electrons for each state
      integer neb !< [Input] The number of \f$\beta\f$-spin
                  !< electrons for each state
      double precision h1(nbf,nbf) !< [Input] The 1-electron integrals
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                            !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                             !< integrals for exchange
      double precision ov(nbf,nbf) !< [Input] The overlap integrals
      double precision factor !< [Input] The scale factor for the 
                              !< occupation penality function
      double precision ofac   !< [Input] The scale factor for the 
                              !< orbital order penality function
      double precision r_noa(nbf,nbf) !< [Input] The 
      !< \f$\alpha\f$-spin natural orbitals. This is a 2D array with the
      !< dimensions:
      !< - 1. The basis functions
      !< - 2. The natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] The 
      !< \f$\beta\f$-spin natural orbitals. This is a 2D array with the
      !< same dimensions as `r_noa`.
      double precision r_ofa(nbf,nbf) !< [Input] The
      !< \f$\alpha\f$-spin occupation functions. This is a 2D array with
      !< the dimension:
      !< - 1. The natural orbitals
      !< - 2. The occupation functions
      double precision r_ofb(nbf,nbf) !< [Input] The 
      !< \f$\beta\f$-spin occupation functions. This is a 2D array with
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
      double precision enta(nbf) !< [Output] The alpha orbital entropies
      double precision entb(nbf) !< [Output] The beta orbital entropies
      double precision power !< [Input] The matrix power
c
      double precision occa(nbf) !< [Input] The alpha electron 
      !< occupation numbers
      double precision occb(nbf) !< [Input] The beta electron 
      !< occupation numbers
c
      double precision En !< [Output] The total energies of all 
                          !< states
      double precision deviation !< [Output] The occupation error
      double precision oerr !< [Output] The orbital ordering error
c
      logical oprint !< [Input] If .true. print stuff.
c
C     Local Variables
C
      double precision E1el
      double precision E2elc
      double precision E2elcsr
      double precision E2elx
      double precision E2elxsr
      double precision Eocc
      double precision Eorder
      double precision ES2err
      double precision Exc
      double precision Sa, Sb
      double precision s2t, s2a
c
      double precision E1a(nbf),  E1b(nbf)
      double precision E2Ca(nbf), E2Cb(nbf)
      double precision E2Xa(nbf), E2Xb(nbf)
c
      integer len
      integer l_m1, k_m1
      integer l_m2, k_m2
      integer l_m3, k_m3
      integer l_m4, k_m4
      integer l_m5, k_m5
      integer l_v1, k_v1
      integer l_v2, k_v2
      integer l_v3, k_v3
      integer ii
C
C     Functions
C
c     double precision wfn1_s2
c     external         wfn1_s2
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
      Eocc    = 0.0d0
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
c
c     s2t = (0.5d0*(nea-neb))*(0.5d0*(nea-neb)+1.0d0)
c     s2a = wfn1_s2(nbf,nea,neb,r_ofa,r_ofb,r_noa,r_nob,ov,
c    +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_v1))
c     serr = (s2a-s2t)**1
c     ES2err = sfac*serr
      ES2err = 0.0d0
c
      call wfn1_1el(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,
     +              h1,E1el,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_v1),dbl_mb(k_v2))
      call wfn1_2el_c(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,
     +              eri,E2elc,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
      call wfn1_2el_hfx(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,ov,
     +              erix,E2elx,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2),power)
      call wfn1_2el_csr(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,
     +              eri,E2elcsr,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
      call wfn1_2el_csr(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,
     +              erix,E2elxsr,
     +              dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3),dbl_mb(k_m4),
     +              dbl_mb(k_v1),dbl_mb(k_v2))
      call wfn1_occ(nbf,nea,neb,r_ofa,r_ofb,occa,occb,factor,Eocc,
     +              deviation,dbl_mb(k_v1))
      call wfn1_entropy(.true.,nbf,nea,r_ofa,enta,dbl_mb(k_v1),
     +              dbl_mb(k_v2))
      call wfn1_entropy(.true.,nbf,neb,r_ofb,entb,dbl_mb(k_v1),
     +              dbl_mb(k_v2))
      call wfn1_xc(rtdb,nbf,nea,neb,r_noa,r_nob,
     +             r_ofa,r_ofb,Exc,dbl_mb(k_m1),
     +             dbl_mb(k_m2))
      call dfill(nbf,0.0d0,ewfa,1)
      call dfill(nbf,0.0d0,ewfb,1)
      call dfill(nbf,0.0d0,E1a,1)
      call dfill(nbf,0.0d0,E1b,1)
      call dfill(nbf,0.0d0,E2Ca,1)
      call dfill(nbf,0.0d0,E2Cb,1)
      call dfill(nbf,0.0d0,E2Xa,1)
      call dfill(nbf,0.0d0,E2Xb,1)
      call wfn1_1el_orb(nbf,nea,r_noa,r_ofa,h1,E1a,dbl_mb(k_m1),
     +                  dbl_mb(k_v1))
      call wfn1_1el_orb(nbf,neb,r_nob,r_ofb,h1,E1b,dbl_mb(k_m1),
     +                  dbl_mb(k_v1))
      call wfn1_2el_C_orb(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,eri,
     +                  E2Ca,E2Cb,
     +                  dbl_mb(k_m1),dbl_mb(k_m2),
     +                  dbl_mb(k_m3),dbl_mb(k_m4),
     +                  dbl_mb(k_v1),dbl_mb(k_v2))
      call wfn1_2el_hfX_orb(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,erix,
     +                  E2Xa,E2Xb,
     +                  dbl_mb(k_m1),dbl_mb(k_m2),
     +                  dbl_mb(k_m3),dbl_mb(k_m4),
     +                  dbl_mb(k_v1),dbl_mb(k_v2))
      do ii = 1, nbf
        ewfa(ii) = E1a(ii) + E2Ca(ii) + E2Xa(ii)
        ewfb(ii) = E1b(ii) + E2Cb(ii) + E2Xb(ii)
        if (oprint) then
          write(*,'(" orb: ",i3,8f22.8)')ii,ewfa(ii),ewfb(ii),
     +          E1a(ii),E1b(ii),E2Ca(ii),E2Cb(ii),E2Xa(ii),E2Xb(ii)
        endif
      enddo
c     call wfn1_l1(do_grad,nbf,nst,l_occa,l_occb,r_noa,r_nob,
c    +             r_l1a,r_l1b,ov,El1,r_dnoa,r_dnob,r_dl1a,r_dl1b,
c    +             dbl_mb(k_m1),dbl_mb(k_m2),dbl_mb(k_m3))
c     call wfn1_l3(do_grad,nbf,nst,nea,neb,r_ofa,r_ofb,r_l3a,r_l3b,
c    +             El3,r_dofa,r_dofb,r_dl3a,r_dl3b,dbl_mb(k_m1))
c
      E2elcsr = -E2elcsr
c     En = E1el+E2elc+E2elx+E2elcsr+E2elxsr+Eocc
c     En = E1el+E2elc+E2elx+Eocc
      En = E1el+E2elc+E2elx
      Sa = 0.0d0
      Sb = 0.0d0
c     do ii = 1, nbf
c       if (oprint) then
c         write(*,'(" alpha ent: ",i4,4f14.8)')ii,
c    +      (ewfa(ii)-ehfa(ii))*enta(ii)/(3.0d0/2.0d0),
c    +      enta(ii),ehfa(ii),ewfa(ii)
c       endif
c     enddo
      if (oprint) then
        write(*,*)
      endif
      do ii = 1, nea
        Sa = Sa + (ewfa(ii)-ehfa(ii))*enta(ii)*(3.0d0/2.0d0)
      enddo
c     do ii = 1, nbf
c       if (oprint) then
c         write(*,'(" beta  ent: ",i4,4f14.8)')ii,
c    +      (ewfb(ii)-ehfb(ii))*entb(ii)/(3.0d0/2.0d0),
c    +      entb(ii),ehfb(ii),ewfb(ii)
c       endif
c     enddo
      if (oprint) then
        write(*,*)
      endif
      do ii = 1, neb
        Sb = Sb + (ewfb(ii)-ehfb(ii))*entb(ii)*(3.0d0/2.0d0)
      enddo
      En = En + Sa + Sb
c
c     Add penalty function to ensure that ewf is always larger than
c     ehf.
c
      Eorder = 0.0d0
      call wfn1_order(nbf,nea,neb,ehfa,ehfb,ewfa,ewfb,ofac,
     +                Eorder,oerr)
c     do ii = 1, nea
c       Eorder = Eorder
c    +  + factor*1.0d6*(abs(ewfa(ii)-ehfa(ii))-(ewfa(ii)-ehfa(ii)))
c     enddo
c     do ii = 1, neb
c       Eorder = Eorder
c    +  + factor*1.0d6*(abs(ewfb(ii)-ehfb(ii))-(ewfb(ii)-ehfb(ii)))
c     enddo
c     En = En + Eorder
cDEBUG
c       write(*,*)'The 1-electron energy             = ',E1el
c       write(*,*)'The Coulomb energy                = ',E2elc
c       write(*,*)'The HF exchange energy            = ',E2elx
c       write(*,*)'The Coulomb self repulsion energy = ',E2elcsr
c       write(*,*)'The self exchange energy          = ',E2elxsr
c       write(*,*)'The occupation number term        = ',Eocc
c       write(*,*)'The orbital order penalty         = ',Eorder
c       write(*,*)'The total wfn1 energy             = ',En
c       write(*,*)
cDEBUG
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

