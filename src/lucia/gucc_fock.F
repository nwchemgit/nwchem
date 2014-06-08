************************************************************************
      subroutine gucc_fock(irefspc,itrefspc,
     &                     ccvec1,ccvec2,ccvec3,ccvec4,
     &                     civec1,civec2,c2vec,
     &                     n_cc_typ,i_cc_typ,
     &                     namp_cc_typ,ioff_cc_typ,
     &                     iopsym,mxb_ci,
     &                     luc,luamp,luomg,ludia,
     &                     luec,luhc)
************************************************************************
*
* purpose : driver for tests of fock-space GUCC
*
*  ak, early 2004
*
************************************************************************
*
* units:
*   luc   = definition of reference function
*   luamp = amplitude vectors (also output for most recent vector)
*   luampold = scratch containing old vectors from previous iterations
*           (on input it may also be a first trial vector)
*   luomg = error vectors
*   ludia = diagonal preconditioner
*   luec  = scratch for exp(G)|ref>
*   luhc  = scratch for H exp(G)|ref>

* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cecore.inc'
      include 'gtbce.inc'
      include 'cintfo.inc'
      include 'glbbas.inc'
      include 'opti.inc'
* constants
      integer, parameter ::
     &     ntest = 5

* arrays
      integer ::
     &     ioff_cc_typ(n_cc_typ), namp_cc_typ(n_cc_typ),
     &     i_cc_typ(4*ngas,n_cc_typ)
      real*8 ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp)
* local
      logical ::
     &     calc_Omg, calc_gradE, tstgrad, tst_hss, comm_ops,
     &     do_eag, do_foo, do_hss
      character*8 cctype
      integer ::
     &     ictp(n_cc_typ)
* external functions
      real*8, external :: inprod, inprdd

* ============================================================
* initialize : restart, set coefs to zero
* ============================================================

      call atim(cpu0,wall0)

      nprint = 1
      lblk = -1

      if (ntest.ge.5) then
        write(6,*) '======================='
        write(6,*) 'entered gucc_fock with:'
        write(6,*) '======================='
        write(6,*) ' iopsym = ',iopsym
      end if

* unit inits
      lusc1 = iopen_nus('GTBSCR1')
      lusc2 = iopen_nus('GTBSCR2')
      lusc3 = iopen_nus('GTBSCR3')
      lusc4 = iopen_nus('GTBSCR4')
      lusc5 = iopen_nus('GTBSCR5')
c      lusc6 = iopen_nus('GTBSCR6')
c      lusc7 = iopen_nus('GTBSCR7')
c      lusc8 = iopen_nus('GTBSCR8')
c      lusc9 = iopen_nus('GTBSCR9')

      call memman(idum,idum,'MARK  ',idum,'FOCKSP')

! be sure to have the correct h1 integrals:
      call copvec(work(kint1o),work(kint1),nint1)

      npairs = 4*ntoob**2
      nheff = npairs*(npairs+1)/2
      call memman(kheff,nheff,'ADDL  ',2,'H_EFF ')
      call memman(kheff2,npairs*npairs,'ADDL  ',2,'H_EFF ')
      nwmat = npairs*npairs
      call memman(kwmat,nwmat,'ADDL  ',2,'W_MAT ')
      call memman(keig,npairs,'ADDL  ',2,'EIGEN ')

      npairs = (2*ntoob-1)*(2*ntoob)/2
      nel = nactel

      print *,' NEL = ', nel

      call diag_effH(work(kheff),work(kheff2),
     &     work(kwmat),work(keig),nel)

      call memchk2('a effH')

c     call logwmat(work(kwmat))

      call cleanvec(work(kwmat),npairs**2)

      call rearr_g(work(kwmat),ccvec1,
     &     n_cc_typ,i_cc_typ,ioff_cc_typ)

      igtbmod = 0
      iopsym = 0
      
      ! important:
      call scalve(ccvec1,-1d0,n_cc_amp)
      call gtbce_E(igtbmod,elen1,variance,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
      h1 = elen1+ecore

      ! effective Hamiltonian:
      call rearr_g(work(kheff2),ccvec3,
     &     n_cc_typ,i_cc_typ,ioff_cc_typ)

      call memchk2('a r g')

      ! recalc energy with effective hamilt.
      isigden = 1
      call sigden_cc(civec1,civec2,luc,lusc1,ccvec3,isigden)      
      hf1 = ecore - inprdd(civec1,civec2,luc,lusc1,1,-1)
      print *,'TEST: E(HF) = ',hf1


      ! create a 2-fold excited determinant
      ccvec2(ioff_cc_typ(25)+2) = 1d0
      isigden = 1
      call sigden_cc(civec1,civec2,luc,lusc1,ccvec2,isigden)

      call gtbce_E(igtbmod,elen2,variance,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               lusc1,lusc2,lusc3,lusc4,lusc5)
                    !ref   eGref HeGref

      h2 = elen2+ecore
      ! <ex|e(-G)He(G)|0>:
      h12 = inprdd(civec1,civec2,lusc2,luhc,1,-1)
      h21 = inprdd(civec1,civec2,lusc3,luec,1,-1)

      print *,' :: ',h1 , h21
      print *,' :: ',h12, h2

      ! Heff|expG 1> on lusc4
      isigden = 1
      call sigden_cc(civec1,civec2,luec,lusc4,ccvec3,isigden)      
      h1t = ecore - inprdd(civec1,civec2,luec,lusc4,1,-1)
      ! Heff|expG 2> on lusc5
      isigden = 1
      call sigden_cc(civec1,civec2,lusc2,lusc5,ccvec3,isigden)      
      h2t = ecore - inprdd(civec1,civec2,lusc2,lusc5,1,-1)

      h12t = -inprdd(civec1,civec2,lusc2,lusc4,1,-1)
      h21t = -inprdd(civec1,civec2,luec,lusc5,1,-1)

      print *,'TEST:'
      print *,' :: ',h1t , h21t
      print *,' :: ',h12t, h2t


      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')
      call relunit(lusc3,'delete')
      call relunit(lusc4,'delete')
      call relunit(lusc5,'delete')
c      call relunit(lusc6,'delete')
c      call relunit(lusc7,'delete')
c      call relunit(lusc8,'delete')
c      call relunit(lusc9,'delete')

      call memman(idum,idum,'FLUSM ',idum,'FOCKSP')

      return

      end

************************************************************************

      subroutine diag_effH(h_eff,h_eff_out,wmat,eig,nel)
************************************************************************
*
*     set up 
*           H(pq,rs) = (pq|rs) + 1/(N-1)(h(pq)dlt(rs)+h(rs)dlt(pq))
*
*     and diagonalize in the sense:
*
*       W(tu,pq)H(pq,rs)W(rs,vw) = e(tu)dlt(tu,vw)
*
************************************************************************

      include 'wrkspc.inc'
      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cecore.inc'
      include 'gtbce.inc'
      include 'opti.inc'
      include 'cintfo.inc'
      include 'glbbas.inc'
      include 'oper.inc'

      dimension h_eff(*), wmat(*), eig(*), h_eff_out(*)

      ! set up h_eff

      npairs = 2*ntoob*(2*ntoob-1)/2
      h_eff(1:npairs*(npairs+1)/2) = 0d0
      xel = dble(nel)

      ipr = 0

      i_unrorb = 0
      ispcas = 0
      i_use_simtrh = 0

      do mp = 1, 2
c old:
c      do ip = 1, ntoob
c        do mr = 1, mp
c
c new:
      do mr = 1, mp
        do ip = 1, ntoob
c
          irmx = ntoob
          if (mr.eq.mp) irmx = ip-1
        do ir = 1, irmx
          ipr = ipr + 1
!          ips = (imp-1)*2*ntoob + ims

          iqs = 0

          do mq = 1, 2
c old:
c          do iq = 1, ntoob
c            do ms = 1, mq
c
c new:
          do ms = 1, mq
            do iq = 1, ntoob
c
              ismx = ntoob
              if (ms.eq.mq) ismx = iq-1
            do is = 1, ismx
              iqs = iqs + 1
!              iqr = (imq-1)*2*ntoob + imr
 
              if (ipr.lt.iqs) cycle ! hermitian operator

              idxpqrs = ipr*(ipr-1)/2 + iqs

c              write(6,*) 'ip,iq,ir,is: ',ip,iq,ir,is
c              write(6,*) 'mp,mq,mr,ms: ',mp,mq,mr,ms
cc
c              write(6,*) 'ipr,iqs,idxpqrs: ', ipr,iqs,idxpqrs
cc
              if (mp.eq.mq.and.mr.eq.ms) then
c
c                write(6,*) ' before    : ', h_eff(idxpqrs)
c
                h_eff(idxpqrs) = h_eff(idxpqrs) + gtijkl(ip,iq,ir,is)
c                
c                write(6,*) ' +2-el int : ', h_eff(idxpqrs)
c                  
                if (ip.eq.iq)
     &               h_eff(idxpqrs) = h_eff(idxpqrs) +
     &                                  1d0/(xel-1d0)*geth1i_2(ir,is)

                if (ir.eq.is)
     &               h_eff(idxpqrs) = h_eff(idxpqrs) +
     &                                  1d0/(xel-1d0)*geth1i_2(ip,iq)

c                write(6,*) ' +1-el int : ', h_eff(idxpqrs)
              end if
              if (mp.eq.ms.and.mq.eq.mr) then
c
c                write(6,*) ' before    : ', h_eff(idxpqrs)
c
                h_eff(idxpqrs) = h_eff(idxpqrs) - gtijkl(ip,is,ir,iq)
                
c                write(6,*) ' +2-el int : ', h_eff(idxpqrs)
                  
                if (ip.eq.is)
     &               h_eff(idxpqrs) = h_eff(idxpqrs) -
     &                                  1d0/(xel-1d0)*geth1i_2(ir,iq)

                if (ir.eq.iq)
     &               h_eff(idxpqrs) = h_eff(idxpqrs) -
     &                                  1d0/(xel-1d0)*geth1i_2(ip,is)

c                write(6,*) ' +1-el int : ', h_eff(idxpqrs)

              end if

              idx2pqrs = (ipr-1)*npairs+iqs
              idx2qpsr = (iqs-1)*npairs+ipr

              h_eff_out(idx2pqrs) = h_eff(idxpqrs)
              h_eff_out(idx2qpsr) = h_eff(idxpqrs)

            end do
            end do
          end do
          end do
        end do
        end do
      end do
      end do

      write(6,*) 'the hamiltonian:'
c      write(6,*) ' norm^2 = ',ddot(npairs**2,h_eff,1,h_eff,1)
      
      call prtrlt(h_eff,npairs)

      call copdia(h_eff,eig,npairs,1)

      write(6,*) 'the diagonal of the hamiltonian:'

      call wrtmat(eig,npairs,1,npairs,1)

      ! call diagonalizer

c      call eigen(h_eff,wmat,npairs,0,1)
      call eigen(h_eff,wmat,npairs,0,0)

      call copdia(h_eff,eig,npairs,1)

      write(6,*) 'the eigenvalues of the hamiltonian:'

      call wrtmat(eig,npairs,1,npairs,1)

      call memman(keigr,npairs,'ADDL  ',2,'EIG RE')
      call memman(keigi,npairs,'ADDL  ',2,'EIG IM')

      call dumsort(eig,npairs,npairs,work(keigr),work(keigi))
      call reo_vec(work(keigr),eig,npairs,work(keigi),1)
      call copvec(work(keigi),eig,npairs)

      write(6,*) 'the eigenvalues of the hamiltonian (sorted):'

      call wrtmat(eig,npairs,1,npairs,1)

      write(6,*) 'E_nuc = ',ecore

      ener = ecore
      do ii = 1, nel*(nel-1)/2
        ener = ener + eig(ii)
        print *,ii,eig(ii),ener
      end do

      write(6,*) 'the W-matrix: ',sqrt(ddot(npairs**2,wmat,1,wmat,1))

      call wrtmat2(wmat,npairs,npairs,npairs,npairs)

c      return

c      call memman(keigr,npairs,'ADDL  ',2,'EIG RE')
c      call memman(keigi,npairs,'ADDL  ',2,'EIG IM')
      call memman(klneig,npairs,'ADDL  ',2,'EIG IM')
      call memman(keigv,npairs**2,'ADDL  ',2,'EIGVEC')
      call memman(keigvr,npairs**2,'ADDL  ',2,'EIGVEC')
      call memman(keigvi,npairs**2,'ADDL  ',2,'EIGVEC')
      call memman(keigvr2,npairs**2,'ADDL  ',2,'EIGVEC')
      call memman(keigvi2,npairs**2,'ADDL  ',2,'EIGVEC')
      call memman(kgmatr,npairs**2,'ADDL  ',2,'GMAT R')
      call memman(kgmati,npairs**2,'ADDL  ',2,'GMAT I')
      call memman(kiscr,npairs,'ADDL  ',1,'INTSCR')
      call memman(kscr,npairs,'ADDL  ',2,'RELSCR')

      call logumat(npairs,work(kgmatr),wmat,
     &     work(keigv),work(keigvr),work(keigvi))
c
c      call copvec(wmat,work(keigvr),npairs**2)
c
c      call rg(npairs,npairs,work(keigvr),
c     &     work(keigr),work(keigi),1,work(keigv),
c     &     work(kiscr),work(kscr),ierr)
c      call nrmvec(npairs,work(keigv),work(keigi))
c
c      print *,'eigenvalues of W-matrix:'
c
c      do ii = 1, npairs
c        eigr = work(keigr-1+ii)
c        eigi = work(keigi-1+ii)
c        rad  = eigr*eigr+eigi*eigi
c        costau = eigr
c        tau = acos(eigr)*sign(1d0,eigi)
c        work(klneig-1+ii) = tau
c        print '(i4,2(2x,e20.10),2(2x,f15.10))',ii,eigr,eigi,rad,
c     &       tau/3.1415926535897932d0*180d0
c      end do
c
cc      print *,'the eigenvectors'
cc      call wrtmat2(work(keigv),npairs,npairs,npairs,npairs)
c
c      ! sort components into real and imaginary matrices
c
c      ii = 0
c      do while(ii.lt.npairs)
c        ii = ii+1
c        eigr = work(keigr-1+ii)
c        eigi = work(keigi-1+ii)
c        if (eigi.eq.0d0) then
c          call copvec(work(keigv+(ii-1)*npairs),
c     &                work(keigvr+(ii-1)*npairs),npairs)
c          call setvec(work(keigvi+(ii-1)*npairs),0d0,npairs)
c        else
c          call copvec(work(keigv+(ii-1)*npairs),
c     &                work(keigvr+(ii-1)*npairs),npairs)
c          call copvec(work(keigv+(ii-1)*npairs),
c     &                work(keigvr+(ii)*npairs),npairs)
c          call copvec(work(keigv+(ii)*npairs),
c     &                work(keigvi+(ii-1)*npairs),npairs)
c          call copvec(work(keigv+(ii)*npairs),
c     &                work(keigvi+(ii)*npairs),npairs)
c          call scalve(work(keigvi+(ii)*npairs),-1d0,npairs)
c          ii = ii+1  ! additional increment
c        end if
c
c      end do
c
cc      print *,'the eigenvectors (re):'
cc      call wrtmat2(work(keigvr),npairs,npairs,npairs,npairs)
cc      print *,'the eigenvectors (im):'
cc      call wrtmat2(work(keigvi),npairs,npairs,npairs,npairs)
c
cc      ! test the trafo matrices:
cc
cc      ! A^T A - B^T B:
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &     1d0,work(keigvr),npairs,
cc     &         work(keigvr),npairs,
cc     &     0d0,work(keigvr2),npairs)
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &     1d0,work(keigvi),npairs,
cc     &         work(keigvi),npairs,
cc     &     1d0,work(keigvr2),npairs)
cc
cc      ! A^T B + B^T A:
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &     1d0,work(keigvr),npairs,
cc     &         work(keigvi),npairs,
cc     &     0d0,work(keigvi2),npairs)
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &    -1d0,work(keigvi),npairs,
cc     &         work(keigvr),npairs,
cc     &     1d0,work(keigvi2),npairs)
cc
cc      print *,'U^T U, real part:'
cc      call wrtmat2(work(keigvr2),npairs,npairs,npairs,npairs)
cc
cc      print *,'U^T U, imaginary part:'
cc      call wrtmat2(work(keigvi2),npairs,npairs,npairs,npairs)
cc
cc
cc      ! W A
cc      call dgemm('n','n',npairs,npairs,npairs,
cc     &           1d0,wmat,npairs,
cc     &               work(keigvr),npairs,
cc     &           0d0,work(keigvr2),npairs)
cc
cc      ! W B
cc      call dgemm('n','n',npairs,npairs,npairs,
cc     &           1d0,wmat,npairs,
cc     &               work(keigvi),npairs,
cc     &           0d0,work(keigvi2),npairs)
cc
cc      ! Re D = A^T W A - B^T W B
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &           1d0,work(keigvi),npairs,
cc     &               work(keigvi2),npairs,
cc     &           0d0,work(kgmatr),npairs)
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &           1d0,work(keigvr),npairs,
cc     &               work(keigvr2),npairs,
cc     &           1d0,work(kgmatr),npairs)
cc
cc      ! Im D = A^T W B + B^T W A
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &           1d0,work(keigvr),npairs,
cc     &               work(keigvi2),npairs,
cc     &           0d0,work(kgmati),npairs)
cc      call dgemm('t','n',npairs,npairs,npairs,
cc     &          -1d0,work(keigvi),npairs,
cc     &               work(keigvr2),npairs,
cc     &           1d0,work(kgmati),npairs)
cc
cc      print *,'U^T W U, real part:'
cc      xnrm = sqrt(ddot(npairs**2,work(kgmatr),1,work(kgmatr),1 ))
cc      print *,'norm = ',xnrm
cc
cc      call wrtmat2(work(kgmatr),npairs,npairs,npairs,npairs)
cc
cc      print *,'U^T W U, imaginary part:'
cc      call wrtmat2(work(kgmati),npairs,npairs,npairs,npairs)
c
c      ! D A^+
c      do ii = 1, npairs
c        tau = work(klneig-1+ii) 
c        do jj = 1, npairs
c          ijadr = (jj-1)*npairs+ii-1  ! uahhh
c          jiadr = (ii-1)*npairs+jj-1  ! uahhh
c          work(keigvr2+ijadr) = tau * work(keigvr+jiadr)
c        end do
c      end do
c
c      ! D B^+
c      do ii = 1, npairs
c        tau = work(klneig-1+ii) 
c        do jj = 1, npairs
c          ijadr = (jj-1)*npairs+ii-1
c          jiadr = (ii-1)*npairs+jj-1
c          work(keigvi2+ijadr) = tau * work(keigvi+jiadr)
c        end do
c      end do
c
c      ! Re G = - B D A^+ + A D B^+
c      call dgemm('n','n',npairs,npairs,npairs,
c     &           -1d0,work(keigvi),npairs,
c     &                work(keigvr2),npairs,
c     &            0d0,work(kgmatr),npairs)
c      call dgemm('n','n',npairs,npairs,npairs,
c     &           +1d0,work(keigvr),npairs,
c     &                work(keigvi2),npairs,
c     &            1d0,work(kgmatr),npairs)
c
c      ! Im G = A D A^+ + B D B^+
c      call dgemm('n','n',npairs,npairs,npairs,
c     &            1d0,work(keigvr),npairs,
c     &                work(keigvr2),npairs,
c     &            0d0,work(kgmati),npairs)
c      call dgemm('n','n',npairs,npairs,npairs,
c     &            1d0,work(keigvi),npairs,
c     &                work(keigvi2),npairs,
c     &            1d0,work(kgmati),npairs)
c
c      print *,'G, real part:'
c
c      call wrtmat2(work(kgmatr),npairs,npairs,npairs,npairs)
c
c      print *,'G, imaginary part:'
c      xnrm = sqrt(ddot(npairs**2,work(kgmatr),1,work(kgmatr),1 ))
c      print *,'norm = ',xnrm
c
c      if (xnrm.gt.1d-10)
c     &     call wrtmat2(work(kgmati),npairs,npairs,npairs,npairs)

      call expmat(work(kgmati),work(kgmatr),work(keigv),work(keigvr),
     &     npairs,1d-20)

      print *,'W from exp(G):'
      call wrtmat2(work(kgmati),npairs,npairs,npairs,npairs)
**

      call copvec(work(kgmatr),wmat,npairs**2)

c      stop 'end of test'

      end

      subroutine expmat(expx,xmat,xscr1,xscr2,ndim,thrsh)

*     calculate exp(X), returned on expx, of (ndim,ndim)-matrix X,
*     input on xmat, by Taylor expansion (threshold thrsh)
*     xscr is a scratch matrix of the same dimensions as xmat, expx

      implicit none

      integer, parameter :: ntest = 100, maxn = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     thrsh
      real(8), intent(inout) ::
     &     expx(ndim,ndim), xmat(ndim,ndim),
     &     xscr1(ndim,ndim),  xscr2(ndim,ndim)

      logical ::
     &     conv
      integer ::
     &     n, ndim2, ii
      real(8) ::
     &     xnrm, fac

      real(8), external ::
     &     inprod

      expx(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)
      xscr2(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)

      do ii = 1, ndim
        expx(ii,ii) = expx(ii,ii) + 1d0
      end do

      ndim2 = ndim*ndim
      n = 1
      conv = .false.

      do while (.not.conv)
        n = n+1
        if (n.gt.maxn) exit

        fac = 1d0/dble(n)

        ! Xscr = 1/N Xscr * X 
        call matml7(xscr1,xscr2,xmat,
     &              ndim,ndim,
     &              ndim,ndim,
     &              ndim,ndim,0d0,fac,0)

        xnrm = sqrt(inprod(xscr1,xscr1,ndim2))
        if (xnrm.lt.thrsh) conv = .true.

        if (ntest.ge.10)
     &       write(6,*) ' N = ',n,'  |1/N! X^N| = ',xnrm

        call vecsum(expx,expx,xscr1,1d0,1d0,ndim2)

        call copvec(xscr1,xscr2,ndim2)

      end do

      if (.not.conv) then
        write(6,*) ' Taylor expansion of exp(X) did not converge!'
        stop 'expmat'
      end if

      return
      end
c
      subroutine rearr_g(gmat_in,gmat_out,ntss_tp,itss_tp,ibtss_tp)

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 1000

* input
      real(8), intent(inout) ::
     &     gmat_in(2*ntoob*(2*ntoob-1)/2,2*ntoob*(2*ntoob-1)/2),
     &     gmat_out(*)

c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)


      npairs = 2*ntoob*(2*ntoob-1)/2

      if (ntest.ge.1000) then
        write(6,*) ' input amplitudes: '
        call wrtmat2(gmat_in,npairs,npairs,npairs,npairs)
        write(6,*) 'ibtss_tp:'
        call iwrtma(ibtss_tp,1,ntss_tp,1,ntss_tp)
      end if
        

        ! loop over types
        do itss = 1, ntss_tp
          idx = ibtss_tp(itss) - 1

          ! identify two-particle excitations:
          nel_ca = ielsum(itss_tp(1,1,itss),ngas)
          nel_cb = ielsum(itss_tp(1,2,itss),ngas)
          nel_aa = ielsum(itss_tp(1,3,itss),ngas)
          nel_ab = ielsum(itss_tp(1,4,itss),ngas)
          nc = nel_ca + nel_cb
          na = nel_aa + nel_ab
          if (na.ne.2) cycle
          ! transform occupations to groups
          call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
          call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
          call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
          call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)
          
          if (mscomb_cc.ne.0) then
            call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &           itss_tp(1,3,itss),itss_tp(1,4,itss),
     &           ngas,idiag)
          else
            idiag = 0
          end if
        
          ! loop over symmetry blocks
          ism = 1 ! totally symmetric operators
          do ism_c = 1, nsmst
            ism_a = multd2h(ism,ism_c)
            do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      ! idx_c contains:
                      ! p1, p2, h1, h2
                      ! idx_s contains
                      ! mp1, mp2, mh1, mh2
                      
                      print *,'p, r: ',idx_c(1),idx_c(2)
                      print *,'      ',idx_s(1),idx_s(2)
                      print *,'q, s: ',idx_c(3),idx_c(4)
                      print *,'      ',idx_s(3),idx_s(4)

                      idxpr = igtidx(idx_c(1),idx_c(2),
     &                              idx_s(1),idx_s(2),ntoob)
                      idxqs = igtidx(idx_c(3),idx_c(4),
     &                              idx_s(3),idx_s(4),ntoob)

                      print *,' >>> ',idx,' ---> ',idxpr,idxqs 

                      gmat_out(idx) = gmat_in(idxpr,idxqs)

c testing coverage:
                      if (gmat_in(idxpr,idxqs).ne.1d100) then
                        gmat_in(idxpr,idxqs) = 1d100
                      else
                        gmat_in(idxpr,idxqs) =
     &                       gmat_in(idxpr,idxqs)+1d100
                      end if

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
            end do ! ism_ca
          end do ! ism_c 

        end do ! itss

      if (ntest.ge.1000) then
        write(6,*) 'coverage:'
        call wrtmat2(gmat_in,npairs,npairs,npairs,npairs)
        write(6,*) ' output amplitudes: '
        call wrt_cc_vec2(gmat_out,6,'GEN_CC')
      end if

      return
      end

      integer function igtidx(ip,ir,mp,mr,ndim)

      implicit none

      integer, intent(in) ::
     &     ip, ir,
     &     mp, mr, ndim

      integer ::
     &     ioff, idx1, idx2

      ! get ipr:
      if (mp.eq.1.and.mr.eq.1) then
        ioff = 0
        idx1 = max(ip,ir)
        idx2 = min(ip,ir)
      else if (mp.eq.1.and.mr.eq.2) then
        ioff = ndim*(ndim-1)/2 ! NB: no diagonal in a/a block
        idx1 = ir
        idx2 = ip
      else if (mp.eq.2.and.mr.eq.1) then
        ioff = ndim*(ndim-1)/2
        idx1 = ip
        idx2 = ir
      else if (mp.eq.2.and.mr.eq.2) then
        ioff = ndim*(ndim-1)/2 + ndim*ndim
        idx1 = max(ip,ir)
        idx2 = min(ip,ir)
      else
        stop 'mp, mr ?'
      end if

      print *,'ioff,idx1,idx2: ',ioff,idx1,idx2

      if (mp.eq.mr) then
        if (ip.eq.ir) stop 'ip, ir ?'
        if (idx1.eq.1) stop 'max(ip,ir) ?'
        igtidx = ioff + (idx1-1)*(idx1-2)/2+idx2
      else
        igtidx = ioff + (idx1-1)*ndim+idx2
      end if

      print *,' >> gtidx = ',igtidx

      return

      end 

      subroutine cleanvec(vec,ndim)

      implicit none

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     vec(ndim)

      integer ::
     &     ii
      real(8) ::
     &     thr

      thr = 1000d0*epsilon(1d0)
      
      print *,'thr = ',thr
      print *,'ndim = ',ndim

      do ii = 1, ndim
        if (abs(vec(ii)).lt.thr) vec(ii)=0d0
      end do

      return

      end
c $Id$
