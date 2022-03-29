*----------------------------------------------------------------------*
      subroutine occ_orbgrad(i_hybrid,ivcc,
     &                       ccvec1,ccvec2,vec1,vec2,
     &                       ittact,iooexcc,iooexc,
     &                       nooexc,namp,nspin,
     &                       xkapnrm,xkresnrm,xlresnrm,
     &                       lu_ccamp,lu_lamp,lu_sig,
     &                       lukappa,lu_ogrd,lu_odia,
     &                       lu_1den,lu_2den,
     &                       lu1int_o,lu2int_o,
     &                       luintm1,luintm2,luref)
*----------------------------------------------------------------------*
*
*     purpose: driver routine to calculate orbital gradient
*
*----------------------------------------------------------------------*
      include "wrkspc.inc"
      include "orbinp.inc"
      include "glbbas.inc"
      include "cgas.inc"
      include "cands.inc"
      include "cc_exc.inc"
      include "oper.inc"
      include "cstate.inc"
      include "ctcc.inc"

      integer, parameter ::
     &     ntest = 50

      integer, intent(in) ::
     &     namp,
     &     i_hybrid,ivcc,
     &     lu_ccamp, lu_lamp, lu_sig,
     &     lukappa, luintm1, luintm2, luref,
     &     lu1int_o,lu2int_o,
     &     iooexcc(*), iooexc(*), nooexc(nspin), ittact(*)
      real(8), intent(inout) ::
     &     ccvec1(*), ccvec2(*), vec1(*), vec2(*)

      character(8) :: cctype

      logical ::
     &     l_ahap
      real(8), external ::
     &     inprod, inprdd

      call atim(cpu0,wall0)

      if (ntest.ge.10) then
        write(6,*) 'Welcome to occ_orbgrad'
        write(6,*) '======================'
        write(6,*) ' namp = ',namp
        write(6,*) ' i_hybrid = ',i_hybrid
        write(6,*) ' lukappa, lu_ogrd = ',lukappa, lu_ogrd
        write(6,*) ' lu_ccamp, lu_lamp = ',lu_ccamp, lu_lamp
        write(6,*) ' luintm1, luintm2, luref = ',luintm1, luintm2, luref
      end if

      lblk = -1
      cctype(1:6) = 'GEN_CC'
      idum = 0
      icspc = ietspc
      icsm  = irefsm
      isspc = ietspc
      issm  = irefsm
      mx_term = 100
      icc_exc = 1
      xconv = 1.0d-20

      nooexc_tot = nooexc(1)
      if (nspin.eq.2) nooexc_tot = nooexc_tot + nooexc(2) 

      call memman(idum,idum,'MARK  ',idum,'ORBGRD')

      lu_lambda = iopen_nus('CCLAMBDA')
      if (luintm1.lt.0.or.luintm2.lt.0) then
        lusc1 = iopen_nus('CCDENSCR1')
        lusc2 = iopen_nus('CCDENSCR2')
        lusc3 = iopen_nus('CCDENSCR3')
        lusc4 = iopen_nus('CCDENSCR4')
      else
        lusc3 = iopen_nus('CCDENSCR3')
      end if

* calculate 1 and 2 density
      ! get right state
      if (luintm1.gt.0) then
        lurst = luintm1
      else
        lurst = -luintm1
        call vec_from_disc(ccvec1,namp,1,lblk,lu_ccamp)
        call expt_ref2(luref,lurst,lusc1,lusc2,lusc3,xconv,mx_term,
     &               ccvec1,dum,vec1,vec2,namp,cctype,0)
      end if

      if (ivcc.eq.1) then
        ! normalize state
        xnorm = sqrt(inprdd(vec1,vec1,lurst,lurst,1,lblk))
        call sclvcd(lurst,lusc3,1d0/xnorm,vec1,1,lblk)
        lurst = lusc3
      end if

      if (ivcc.eq.0) then
        ! get left state, if necessary
        ! get refstate in ITSPC:
        call expciv(irefsm,ietspc,luref,itspc,lusc3,
     &            lblk,lusc4,1,0,idc,0)
        if (luintm2.gt.0) then
          ! add |ref> to obtain lambda
          lulst = lu_lambda
          call vecsmd(vec1,vec2,1d0,1d0,lusc3,luintm2,lu_lambda,1,lblk)
        else
          lulst = lu_lambda
          ! sum(mu) L(mu) tau(mu) |HF>
          call vec_from_disc(ccvec2,namp,1,lblk,lu_lamp)
          icspc = itspc
          isspc = itspc
          call sig_gcc(vec1,vec2,lusc3,lusc2,ccvec2)
          ! add |hf>
          call vecsmd(vec1,vec2,1d0,1d0,lusc3,lusc2,lusc1,1,lblk)
          ! and multiply with exp(-t\dag)
          if (luintm1.gt.0) 
     &        call vec_from_disc(ccvec1,namp,1,lblk,lu_ccamp)
          call scalve(ccvec1,-1d0,namp)
          ! conjugate amplitudes (reorder)  and operators
          call conj_ccamp(ccvec1,1,ccvec2)
          call conj_t
          ! exp(-t\dag), result on lulst
          call expt_ref2(lusc1,lulst,lusc2,lusc3,lusc4,xconv,mx_term,
     &       ccvec2,dum,vec1,vec2,namp,cctype,0)
          call conj_t
        end if
      else
        ! just get a copy of the right state
        lulst = lu_lambda
        call copvcd(lurst,lulst,vec1,1,lblk)
      end if

      ntbsq = ntoob*ntoob
      lrho1 = nspin*ntbsq
      lrho2 = nspin*ntbsq*(ntbsq+1)/2 + (nspin-1)*ntbsq**2

      iden = 2
      isepab = 0
      if (nspin.eq.2) isepab = 1
c      call densi2(iden,work(krho1),work(krho2),
c     &     vec1,vec2,lulst,lurst,exps2,
c     &     0,work(ksrho1))
      if (ivcc.eq.0) then
        icspc = ietspc
        isspc = itspc
      else
        ! well, no games with VCC; it always needs the full (FCI) space
        icspc = ietspc
        isspc = ietspc
      end if
c dbg
      print *,'call to densi2_ab'
      call util_flush(6)
c dbg
      call densi2_ab(iden,work(krho1),work(krho2),
     &     vec1,vec2,lulst,lurst,exps2,
     &     0,work(ksrho1),isepab)
c dbg
      print *,'after densi2_ab'
      call util_flush(6)
c dbg

      ! symmetrize one-body density matrix
      do ispin = 1, nspin
        ioff = (ispin-1)*ntbsq
        call sym_blmat(work(krho1+ioff),1,ntoob)
      end do
      ! symmetrize two-body density matrix
      do ispc = 1, nspin*2-1
        ioff = (ispc-1)*ntbsq*(ntbsq+1)/2
        imod = 0
        if (ispc.eq.3) imod = 1

        call sym_2dens(work(krho2+ioff),ntoob,imod,0.01d0)
      end do

      if (lu_1den.gt.0) 
     &     call vec_to_disc(work(krho1),lrho1,1,lblk,lu_1den)
      if (lu_2den.gt.0) then
        call vec_to_disc(work(krho2),lrho2,1,lblk,lu_2den)
      end if

* calculate eff. fock matrix
      icc_exc = 0
      call memman(klfoo,nspin*ntoob**2,'ADDL  ',2,'FEFF  ')
      call fock_mat_ab(work(klfoo),2,nspin)
c      call fock_mat(work(klfoo),2)
* get orbital gradient
      call memman(kle1,nooexc_tot,'ADDL  ',2,'E1    ')

      ! Hybrid: the p/h gradient is already on disc, read it
      if (i_hybrid.eq.1)
     &     call vec_from_disc(work(kle1),nooexc_tot,1,lblk,lu_ogrd)

* orbital gradient types:
* 2 :  general orbital gradient at kappa != 0
* 1 :  CEP (current expansion point) orbital gradient
* 3 :  numerical general orbital gradient
      iogtyp = 1
      
      if (iogtyp.eq.1) then
        do ispin = 1, nspin
          ifoo = (ispin-1)*ntoob**2
          ie1  = (ispin-1)*nooexc(1)
          imode = 0
          if (i_hybrid.eq.1) imode=-1
          call f_to_e1(work(klfoo+ifoo),work(kle1+ie1),
     &         1,iooexcc(2*ie1+1),
     &         nooexc(ispin),imode)
          if (ntest.ge.100) then
            if (i_hybrid.eq.1)
     &           write(6,*) 'After adding inactive/active rotations:'
            if (nspin.eq.1)
     &           write(6,*) 'The orbital gradient:'
            if (nspin.eq.2.and.ispin.eq.1)
     &           write(6,*) 'The orbital gradient (alpha):'
            if (nspin.eq.2.and.ispin.eq.2)
     &           write(6,*) 'The orbital gradient (beta):'
            call wrt_excvec(work(kle1+ie1),
     &           iooexcc(2*ie1+1),nooexc(ispin))
          end if
        end do

        if (nspin.eq.2) then
          iexc2 = nooexc(1)
          do iexc = 1, nooexc(1)
            iorb = iooexcc((iexc-1)*2+1) 
            jorb = iooexcc((iexc-1)*2+2)
            itp = itpfto(iorb)
            jtp = itpfto(jorb)
            l_ahap = i_iadx(itp).eq.2.and.i_iadx(jtp).eq.2.and.
     &           ihpvgas(itp).ne.ihpvgas(jtp)
            if (.not.l_ahap) then
              do
                iexc2 = iexc2+1
                if (iexc2.gt.nooexc_tot) stop 'oha!'
                iorb2 = iooexcc((iexc2-1)*2+1)
                jorb2 = iooexcc((iexc2-1)*2+2)
                if (iorb.eq.iorb2.and.jorb.eq.jorb2) then
                  avg = 0.5d0*(work(kle1-1+iexc)+work(kle1-1+iexc2))
                  work(kle1-1+iexc) = avg
                  work(kle1-1+iexc2)= avg
                  exit
                end if
              end do
            end if
          end do
          if (ntest.ge.100) then
            if (i_hybrid.eq.1)
     &           write(6,*) 'After averaging non-ap/ah rotations:'
            do ispin = 1, nspin
              ie1  = (ispin-1)*nooexc(1)
              if (nspin.eq.1)
     &             write(6,*) 'The orbital gradient:'
              if (nspin.eq.2.and.ispin.eq.1)
     &             write(6,*) 'The orbital gradient (alpha):'
              if (nspin.eq.2.and.ispin.eq.2)
     &             write(6,*) 'The orbital gradient (beta):'
              call wrt_excvec(work(kle1+ie1),
     &             iooexcc(2*ie1+1),nooexc(ispin))
            end do
          end if

        end if

      else if (iogtyp.eq.2) then
        stop 'do your work'
      else if (iogtyp.eq.3) then
cedo disabled since routine is missing     
        write(6,*) 'orbgrd_num2 missing '
        stop 'orbgrd_num2'
c        call orbgrd_num2(0,work(kle1),lukappa,
c     &       lu1int_o,lu2int_o,iooexcc,nooexc)
      else
        write(6,*) 'unknown iogtyp (',iogtyp,')'
        stop 'occ_orbgrad'
      end if
      xkresnrm = sqrt(inprod(work(kle1),work(kle1),nooexc_tot))

      if (ntest.ge.10) then
        do ispin = 1, nspin
          ie1  = (ispin-1)*nooexc(1)
          if (nspin.eq.1)
     &         write(6,*) 'The orbital gradient:'
          if (nspin.eq.2.and.ispin.eq.1)
     &         write(6,*) 'The orbital gradient (alpha):'
          if (nspin.eq.2.and.ispin.eq.2)
     &         write(6,*) 'The orbital gradient (beta):'
          call wrt_excvec(work(kle1+ie1),
     &         iooexcc(2*ie1+1),nooexc(ispin))
        end do
      end if

      call vec_to_disc(work(kle1),nooexc_tot,1,lblk,lu_ogrd)

      call relunit(lu_lambda,'delete')

      if (luintm1.lt.0.or.luintm2.lt.0) then
        call relunit(lusc1,'delete')
        call relunit(lusc2,'delete')
        call relunit(lusc3,'delete')
        call relunit(lusc4,'delete')      
      else
        call relunit(lusc3,'delete')
      end if

      ! add orbital-gradient contributions to residual of
      ! left-hand vector
      if (i_hybrid.eq.1) then
        call memman(ko1c,ntoob,'ADDL  ',1,'OGRD C')
        call memman(ko1a,ntoob,'ADDL  ',1,'OGRD A')

        call vec_from_disc(ccvec1,namp,1,lblk,lu_sig)

        do ispin = 1, nspin
          call compress_t1s(ccvec1,work(klfoo),1,work(ko1c),work(ko1a),
     &     work(klsobex),work(klsox_to_ox),
     &     work(klcobex_tp),work(klibsobex),
     &     nspobex_tp,ispin,-2)
        end do

        xlresnrm = sqrt(inprod(ccvec1,ccvec1,namp))

        call vec_to_disc(ccvec1,namp,1,lblk,lu_sig)

      end if

      iaprhess=1
c      iaprhess=0
      if (lu_odia.le.0) iaprhess=0
      if (iaprhess.eq.1) then
CCCCC
c          call memman(kle3,nooexc_tot,'ADDL  ',2,'NT KAP')
c          call memman(kle4,nooexc_tot,'ADDL  ',2,'NT KAP')
CCCCC
        imode = 0
        do ispin = 1, nspin
          ifoo = (ispin-1)*ntoob**2
          ie1  = (ispin-1)*nooexc(1)
          itt  = (ispin-1)*ngas**2

          if (nspin.eq.1) ispc = 0
          if (nspin.eq.2) ispc = ispin
          call diag_orbhes(work(kle1+ie1),work(klfoo+ifoo),
     &         iooexc(1+ifoo),nooexc(ispin),1,ittact(1+itt),ispc)
CCCCC
c          xinc = 0.00001d0
c          luoinc = iopen_uus()
c          lumodu = iopen_uus()
c          
c          idx = 0
c          if (ispin.eq.2) idx = nooexc(1)
c          do iexc = 1, nooexc(ispin)
c            idx = idx+1
c            ! increment kappa
c            work(kle3:kle3-1+nooexc_tot) = 0d0
c            work(kle3-1+idx) = xinc
c            call vec_to_disc(work(kle3),nooexc_tot,1,-1,luoinc)
c            ! get transf.-mat. from current kappa
c            call kap2u(3,-1000,lukappa,lumodu,iooexcc,nooexc,nspin)
c            ! update transf.mat.
c            call kap2u(-3,luoinc,-1000,lumodu,iooexcc,nooexc,nspin)
c            ! get new integrals with mod. transf.mat.
c            call tra_kappa(-1,lumodu,iooexcc,nooexc,nspin,
c     &                      1,lu1int_o,lu2int_o)
c            ! get Fock
c            call fock_mat_ab(work(klfoo),2,nspin)
c            ! get E1(+)
c            do jspin = 1, nspin
c              jfoo = (jspin-1)*ntoob**2
c              je1  = (jspin-1)*nooexc(1)
c              imode = 0
c              call f_to_e1(work(klfoo+jfoo),work(kle3+je1),
c     &             1,iooexcc(2*ie1+1),
c     &             nooexc(jspin),imode)
c            end do
c
c            ! increment kappa -
c            work(kle4:kle4-1+nooexc_tot) = 0d0
c            work(kle4-1+idx) = -xinc
c            call vec_to_disc(work(kle4),nooexc_tot,1,-1,luoinc)
c            ! get transf.-mat. from current kappa
c            call kap2u(3,-1000,lukappa,lumodu,iooexcc,nooexc,nspin)
c            ! update transf.mat.
c            call kap2u(-3,luoinc,-1000,lumodu,iooexcc,nooexc,nspin)
c            ! get new integrals with mod. transf.mat.
c            call tra_kappa(-1,lumodu,iooexcc,nooexc,nspin,
c     &                      1,lu1int_o,lu2int_o)
c            ! get Fock
c            call fock_mat_ab(work(klfoo),2,nspin)
c            ! get E1(-)
c            do jspin = 1, nspin
c              jfoo = (jspin-1)*ntoob**2
c              je1  = (jspin-1)*nooexc(1)
c              imode = 0
c              call f_to_e1(work(klfoo+jfoo),work(kle4+je1),
c     &             1,iooexcc(2*je1+1),
c     &             nooexc(jspin),imode)
c            end do
c
c            iorb = iooexcc((idx-1)*2+1) 
c            jorb = iooexcc((idx-1)*2+2)
c
c            xdia = (work(kle3-1+idx)-work(kle4-1+idx))/(2d0*xinc)
c            print *,' i,j: ',iorb,jorb
c            print *,' numerical  ',xdia
c            print *,' analytical ',work(kle1-1+idx)
c
c            work(kle1-1+idx) = xdia
c          end do
c
c          ! restore law and order:
c          call kap2u(3,-1000,lukappa,lumodu,iooexcc,nooexc,nspin)
c          ! update transf.mat.
c          call tra_kappa(-1,lumodu,iooexcc,nooexc,nspin,
c     &                      1,lu1int_o,lu2int_o)
c
c          call relunit(luoinc,'delete')
c          call relunit(lumodu,'delete')
CCCCC
          if (ntest.ge.100) then
            write(6,*) 'The diagonal orbital Hessian (unmodified):'
            if (ispin.eq.1.and.nspin.eq.2) write(6,*) '   alpha part'
            if (ispin.eq.2.and.nspin.eq.2) write(6,*) '   beta part'
            call wrt_excvec(work(kle1+ie1),
     &                      iooexcc(1+2*ie1),nooexc(ispin))
          end if
        end do

        if (nspin.eq.2) then
          iexc2 = nooexc(1)
          do iexc = 1, nooexc(1)
            iorb = iooexcc((iexc-1)*2+1) 
            jorb = iooexcc((iexc-1)*2+2)
            itp = itpfto(iorb)
            jtp = itpfto(jorb)
            l_ahap = i_iadx(itp).eq.2.and.i_iadx(jtp).eq.2.and.
     &           ihpvgas(itp).ne.ihpvgas(jtp)
            if (.not.l_ahap) then
              do
                iexc2 = iexc2+1
                if (iexc2.gt.nooexc_tot) stop 'oha2!'
                iorb2 = iooexcc((iexc2-1)*2+1)
                jorb2 = iooexcc((iexc2-1)*2+2)
                if (iorb.eq.iorb2.and.jorb.eq.jorb2) then
                  avg = 0.5d0*(work(kle1-1+iexc)+work(kle1-1+iexc2))
                  work(kle1-1+iexc) = avg
                  work(kle1-1+iexc2)= avg
                  exit
                end if
              end do
            end if
          end do
        end if

        xmin = 1d15
        do iexc = 1, nooexc_tot
c          work(kle1-1+iexc) = work(kle1-1+iexc)
c          if (abs(work(kle1-1+iexc)).le.1d-4) work(kle1-1+iexc)=1d+13
          if (work(kle1-1+iexc).lt.-0.2d0)
     &         work(kle1-1+iexc)=-work(kle1-1+iexc)
          xmin = min(xmin,work(kle1-1+iexc))
        end do
        xsh = 0d0
        xshm = 2d-2
        if (xmin.lt.xshm) then
          xsh = xshm - xmin!2d-1 - xmin
          work(kle1:kle1-1+nooexc_tot) =
     &         work(kle1:kle1-1+nooexc_tot)+xsh
        end if
        if (ntest.ge.50) then
          write(6,*) 'The diagonal orbital Hessian (modified):'
          write(6,*) ' shift was: ',xsh
          do ispin = 1, nspin
            ie1  = (ispin-1)*nooexc(1)
            if (ispin.eq.1.and.nspin.eq.2) write(6,*) '   alpha part'
            if (ispin.eq.2.and.nspin.eq.2) write(6,*) '   beta part'
            call wrt_excvec(work(kle1+ie1),
     &                      iooexcc(1+2*ie1),nooexc(ispin))
          end do
        end if
        call vec_to_disc(work(kle1),nooexc_tot,1,lblk,lu_odia)
      end if

      call memman(idum,idum,'FLUSM  ',idum,'ORBGRD')

      call atim(cpu,wall)

      call prtim(6,'time for orbital gradient',cpu-cpu0,wall-wall0)

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine mdfy_hss(xhss,ittact,iooexcc,itpfto,
     &                    ihpvgas_ab,ld_hpv,
     &                    nooexc,ngas,nspin)
*----------------------------------------------------------------------*
*     a diagonal hessian approximated as 
*
*      H(ij,ij) = F_ii - F_jj
*
*     is input. Some empirical rules are used to modify i.pt.
*     inactive/active rotations
*
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nooexc(2),ngas,nspin,ld_hpv,
     &     ittact(ngas,ngas,nspin),iooexcc(2,*),itpfto(*),
     &     ihpvgas_ab(ld_hpv,nspin)
      real(8), intent(inout) ::
     &     xhss(*)

      integer, parameter ::
     &     ntest = 00

      integer ::
     &     nooexc_tot,
     &     ispin, iexc, idx, iorb, jorb, itp, jtp

      nooexc_tot = nooexc(1)
      if (nspin.eq.2) nooexc_tot = nooexc_tot + nooexc(2)

      idx = 0
      do ispin = 1, nspin
c        iooff = (ispin-1)*nooexc(1)
        do iexc = 1, nooexc(ispin)
          idx = idx + 1
          iorb = iooexcc(1,idx)
          jorb = iooexcc(2,idx)
          itp = itpfto(iorb)
          jtp = itpfto(jorb)
          if (ihpvgas_ab(itp,ispin).eq.
     &        ihpvgas_ab(jtp,ispin)    ) then
            xhss(idx) = 0.025d0 * xhss(idx)
          end if
          if (xhss(idx).lt.0d0) xhss(idx)=0.1d0
        
        end do
      end do

      if (ntest.ge.100) then
        write(6,*) 'Modified approx. orbital hessian:'
        idx = 0
        do ispin = 1, nspin
          if (nspin.eq.2.and.ispin.eq.1) write(6,*) 'alpha:'
          if (nspin.eq.2.and.ispin.eq.2) write(6,*) 'beta:'
          do iexc = 1, nooexc(ispin)
            idx = idx+1
            write(6,'(x,2i5,g20.10)')
     &           iooexcc(1,idx),
     &           iooexcc(2,idx),
     &           xhss(idx)
          end do
        end do
      end if

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine new_fock(fock,rho1,diag_h,nooexc,iooexcc)
*----------------------------------------------------------------------*
      
      include "wrkspc.inc"
      include "glbbas.inc"
      include "cintfo.inc"
      include "orbinp.inc"
      include "lucinp.inc"

      integer, parameter ::
     &     ntest = 00

      real(8), intent(inout) ::
     &     fock(*), rho1(*), diag_h(*)
      integer, intent(in) ::
     &     iooexcc(2,*), nooexc

      integer ::
     &     ioff(nsmob)

      call swapve(work(krho1),rho1,ntoob*ntoob)
      call copvec(work(kint1),fock,nint1)
      call fifam(fock)

      if (ntest.ge.100) then
        write(6,*) 'new fock matrix:'
        call aprblm2(fock,ntoobs,ntoobs,nsmob,1)
      end if

      ! precalculate offsets of symmetry blocks
      idx = 0
      do ism = 1, nsmob
        ioff(ism) = idx
        idx = idx + (ntoobs(ism)+1)*ntoobs(ism)/2
      end do

      do iexc = 1, nooexc
        idx = iooexcc(1,iexc)
        jdx = iooexcc(2,iexc)
        ! convert type to symmetry ordering:
        ism = ismfto(idx)
        jsm = ismfto(jdx)
        ii = ireots(idx) - ibso(ism) + 1
        jj = ireots(jdx) - ibso(jsm) + 1
        ! diagonal element addresses in *upper* triangles
        iid = ioff(ism) + (ii+1)*ii/2
        jjd = ioff(jsm) + (jj+1)*jj/2
        
          print *,'iexc, idx, jdx, ii, jj: ',iexc, idx, jdx, ii, jj
          print *,'   ism, jsm, ndimi, ndimj: ', ism, jsm, ibso(ism)
          print *,'   iid, jjd: ', iid, jjd
          print *,'             ',fock(iid) , fock(jjd)

        diag_h(iexc) =  4d0 * (fock(iid) - fock(jjd))
      end do
      if (ntest.ge.50) then
        write(6,*) 'diagonal Orbital Hessian:'
        call wrtmat(diag_h,nooexc,1,nooexc,1)
      end if

      call swapve(work(krho1),rho1,ntoob*ntoob)

      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine omg2ogrd(luomg,lu_ogrd,ccvec1,
     &             iooexcc,nooexc,n_cc_amp,nspin,
     &             xkresnrm,xomgnrm)
*----------------------------------------------------------------------*
*
*     rearrage single-excitation part of Omega as orbital gradient
*     for Brueckner CC
*
*----------------------------------------------------------------------*
      include "wrkspc.inc"
      include "orbinp.inc"
      include "lucinp.inc"
      include "ctcc.inc"

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     luomg, lu_ogrd,
     &     nooexc(2), n_cc_amp, iooexcc(*), nspin
      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp)
      real(8), intent(out) ::
     &     xkresnrm,xomgnrm

      real(8), external ::
     &     inprod

      if (ntest.ge.10) then
        write(6,*) 'OMG2OGRD'
        write(6,*) '========'
        write(6,*) ' luomg, lu_ogrd: ',luomg, lu_ogrd
        write(6,*) ' nooexc, n_cc_amp: ',nooexc, n_cc_amp, nspin
      end if

      lblk = -1
      idum = 0
      call memman(idum,idum,'MARK  ',idum,'OM2OGR')

      nooexc_tot = nooexc(1)
      if (nspin.eq.2) nooexc_tot = nooexc_tot + nooexc(2)

      ! get workspace
      logrd = 0
      do ism = 1, nsmob
        logrd = logrd + ntoobs(ism)*ntoobs(ism)
      end do

      call memman(kogrd,logrd,'ADDL  ',2,'OGRD  ')
      call memman(kogrdc,nooexc_tot,'ADDL  ',2,'OGRDC ')
      call memman(ko1c,ntoob,'ADDL  ',1,'OGRD C')
      call memman(ko1a,ntoob,'ADDL  ',1,'OGRD A')

      ! read from disc
      call vec_from_disc(ccvec1,n_cc_amp,1,lblk,luomg)
      
      ! and sort singles into new array
      do ispin = 1, nspin
c      iab = 1 ! prelim

        call expand_t1s_new(ccvec1,work(kogrd),1,work(ko1c),work(ko1a),
     &     work(klsobex),work(klsox_to_ox),
     &     work(klcobex_tp),work(klibsobex),
     &     nspobex_tp,ispin)

        ! compress orbital gradient
        ioff = (ispin-1)*nooexc(1)*2
        koff = (ispin-1)*nooexc(1)
        call comprs_kap(work(kogrd),work(kogrdc+koff),
     &       iooexcc(ioff+1),nooexc(ispin))

      end do

      ! save modified omega to disc 
      call zero_t1(ccvec1)
      xomgnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))

      call vec_to_disc(ccvec1,n_cc_amp,1,lblk,luomg)

      ! scale it with factor
      if (nspin.eq.1) fac = 4d0
      if (nspin.eq.2) fac = 2d0
      call scalve(work(kogrdc),fac,nooexc_tot)

      xkresnrm = sqrt(inprod(work(kogrdc),work(kogrdc),nooexc_tot))

      ! save orbital gradient to disc
      call vec_to_disc(work(kogrdc),nooexc_tot,1,lblk,lu_ogrd)

      if (ntest.ge.100) then
        write(6,*) ' Brueckner orbital gradient'
        write(6,*) ' =========================='
        do ispin = 1, nspin
          ioff = (ispin-1)*nooexc(1)*2
          koff = (ispin-1)*nooexc(1)
          if (ispin.eq.1.and.nspin.eq.2) write(6,*) '  alpha'
          if (ispin.eq.2.and.nspin.eq.2) write(6,*) '  beta'
          call wrt_excvec(work(kogrdc+koff),iooexcc(ioff+1),
     &                    nooexc(ispin))
        end do
      end if

      idum = 0
      call memman(idum,idum,'FLUSM ',idum,'OM2OGR')

      return
      end

      subroutine occ_scnd_num(i_obcc,
     &           ccvec1,ccvec2,vec1,vec2,
     &           ittact,iooexcc,iooexc,
     &           nooexc,n_cc_amp,nspin,
     &           lu1into,lu2into,luref,
     &           luccamp,lu_lamp,lukappa,
     &           lutrvec,lutrv_l,lutrv_o,
     &           lursig,lusig_l,lusig_o,
     &           iactt,iactl,iacto)

      include 'implicit.inc'

      integer, intent(in) ::
     &     nooexc(nspin)
      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)

      real(8) ::
     &     xinc(2)
      integer ::
     &     luomg(2), lugrd(2), luogrd(2)

      real(8), external ::
     &     inprdd

      ninc = 2
      delta = 1d-5
      xinc(1) = delta
      xinc(2) = -delta

      lumodt = iopen_nus('NUMMODT')
      lumodl = iopen_nus('NUMMODL')
      lumodo = iopen_nus('NUMMODO')

      luomg(1) = iopen_nus('OMGINC1')
      luomg(2) = iopen_nus('OMGINC2')
      lugrd(1) = iopen_nus('GRDINC1')
      lugrd(2) = iopen_nus('GRDINC2')
      luogrd(1) = iopen_nus('OGRDINC1')
      luogrd(2) = iopen_nus('OGRDINC2')

      luintm1 = iopen_nus('NUMCCINTM1')
      luintm2 = iopen_nus('NUMCCINTM2')
      luintm3 = iopen_nus('NUMCCINTM3')

      lurhs = iopen_nus('NUM_RHS')

      write(6,*) ' NUMERICAL HESSIAN FOR OCC:'

      xnt2 = inprdd(ccvec1,ccvec1,lutrvec,lutrvec,1,-1)
      xnl2 = inprdd(ccvec1,ccvec1,lutrv_l,lutrv_l,1,-1)
      xno2 = inprdd(ccvec1,ccvec1,lutrv_o,lutrv_o,1,-1)
      write(6,*) ' NORM OF INPUT VECTOR: ',sqrt(xnt2+xnl2+xno2)
      write(6,*) ' NORM OF COMPONENTS:   ',
     &     sqrt(xnt2),sqrt(xnl2),sqrt(xno2)

      do iinc = 1, ninc
        write(6,*) ' OUTPUT FOR XINC = ',XINC(IINC),' FOLLOWS:'
        ! modify parameters by increment
        call vecsmd(ccvec1,ccvec2,1d0,xinc(iinc),
     &       luccamp,lutrvec,lumodt,1,-1)
        call vecsmd(ccvec1,ccvec2,1d0,xinc(iinc),
     &       lu_lamp,lutrv_l,lumodl,1,-1)
c V A:
c        call vecsmd(ccvec1,ccvec2,1d0,xinc(iinc),
c     &       lukappa,lutrv_o,lumodo,1,-1)
c
c V B:
        ! get transf.-mat. from current kappa
        call kap2u(3,-1000,lukappa,lumodo,iooexcc,nooexc,nspin)
        call sclvcd(lutrv_o,luintm1,xinc(iinc),ccvec1,1,-1)
        ! update transf.mat.
        call kap2u(-3,luintm1,-1000,lumodo,iooexcc,nooexc,nspin)
c

        ! get new orbitals with mod. transf.mat.
        call tra_kappa(-1,lumodo,iooexcc,nooexc,nspin,
     &                 1,lu1into,lu2into)

        ! call vector function
        call cc_vec_fnc2(ccvec1,ccvec2,ecc,eccl,
     &                   xampnrm,xomgnrm,xdum,
     &                   vec1,vec2,1,'GEN_CC',
     &                   xdum,
     &                   lumodt,luomg(iinc),lumodl,
     &                   luintm1,luintm2,luintm3)

        ! call rhs and jacobian lh trafo
        call zero_ord_rhs(ccvec1,vec1,vec2,lumodt,lurhs,luintm1)
        lr_switch = 1
        iadd_rhs = 1
        itex_sm = 1
        call jac_t_vec2(lr_switch,iadd_rhs,0,1,1,
     &           ccvec1,ccvec2,vec1,vec2,
     &           n_cc_amp,n_cc_amp,
     &           ecc1,xlampnrm,xlresnrm,
     &           lumodt,luomg(iinc),lumodl,lugrd(iinc),lurhs,
     &           luintm1,luintm2,luintm3)

        ! call orbital gradient
        if (iacto.eq.0) then
          namp = sum(nooexc(1:nspin))
          ccvec1(1:namp) = 0d0
          call vec_to_disc(ccvec1,namp,1,-1,luogrd(iinc))
        else
          if (i_obcc.eq.1) then
            call omg2ogrd(luomg(iinc),luogrd(iinc),ccvec1,
     &           iooexcc,nooexc,n_cc_amp,nspin,
     &           xkresnrm,xomgnrm)          
          end if
          ivcc = 0
          call occ_orbgrad(i_obcc,ivcc,
     &         ccvec1,ccvec2,vec1,vec2,
     &         ittact,iooexcc,iooexc,
     &         nooexc,n_cc_amp,nspin,
     &         xkapnrm,xkresnrm,xlresnrm,
     &         lumodt,lumodl,lugrd(iinc),
     &         ludum,luogrd(iinc),-1,
     &         -1,-1,
     &         lu1into,lu2into,
     &         luintm1,luintm3,luref)
        end if

      end do

      ! calculate numerical sigma-vectors
      fac = 1d0/(2d0*delta)
      call vecsmd(ccvec1,ccvec2,fac,-fac,luomg(1),luomg(2),lursig,1,-1)
      call vecsmd(ccvec1,ccvec2,fac,-fac,lugrd(1),lugrd(2),lusig_l,1,-1)
c TEST --- no t/l contribution
c          do ii = 1, 10
c            print *,'T,L set to 0D0 !!!'
c          end do
c          ccvec1(1:n_cc_amp)=0d0
c          call vec_to_disc(ccvec1,n_cc_amp,1,-1,lursig)
c          call vec_to_disc(ccvec1,n_cc_amp,1,-1,lusig_l)
c TEST
      call vecsmd(ccvec1,ccvec2,fac,-fac,
     &     luogrd(1),luogrd(2),lusig_o,1,-1)
c TEST --- no kappa contribution
c          if (nspin.eq.2) stop 'not possible'
c          do ii = 1, 10
c            print *,'kappa set to 0D0 !!!'
c          end do
c          ccvec1(1:nooexc(1))=0d0
c          call vec_to_disc(ccvec1,nooexc(1),1,-1,lusig_o)
c TEST

      xnt2 = inprdd(ccvec1,ccvec1,lursig,lursig,1,-1)
      xnl2 = inprdd(ccvec1,ccvec1,lusig_l,lusig_l,1,-1)
      xno2 = inprdd(ccvec1,ccvec1,lusig_o,lusig_o,1,-1)
      write(6,*) ' NORM OF OUTPUT VECTOR: ',sqrt(xnt2+xnl2+xno2)
      write(6,*) ' NORM OF COMPONENTS:   ',
     &     sqrt(xnt2),sqrt(xnl2),sqrt(xno2)

      ! delete files
      call relunit(lumodt,'delete')
      call relunit(lumodl,'delete')
      call relunit(lumodo,'delete')

      call relunit(luomg(1),'delete')
      call relunit(luomg(2),'delete')
      call relunit(lugrd(1),'delete')
      call relunit(lugrd(2),'delete')
      call relunit(luogrd(1),'delete')
      call relunit(luogrd(2),'delete')

      call relunit(luintm1,'delete')
      call relunit(luintm2,'delete')
      call relunit(luintm3,'delete')

      call relunit(lurhs,'delete')

      ! restore the old orbitals
      call tra_kappa(lukappa,-1,iooexcc,nooexc,nspin,
     &               1,lu1into,lu2into)

      end
c $Id$
