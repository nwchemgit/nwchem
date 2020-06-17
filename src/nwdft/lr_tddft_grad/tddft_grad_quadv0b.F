      subroutine tddft_grad_quadv0b(ictr, ldew, ipol, grad, 
     1  nbf_ao, mbf_ao, natoms, nq, chi, delchi, 
     2  heschi, Amat, Cmat, ccdel, 
     3  ibf, iniz, ifin, cetobfr,
     4  PPP, F, zzz, hh, g_dens, force, 
     5  oactive, nmat, calctype)
c
c     $Id$
c
c This routine constructs the gradient of XC-functional terms with
c respect to nuclear coordinates.  Ideas for this routine are borrowed
c from both dftg_quadv0b and xc_rhogen (when what = 1).
c
c NOTE!  The Amat and Cmat defined here are not the same as the Amat
c and Cmat from xc_eval_fnl
c 
      Implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "dftpara.fh"
c     
      logical ldew              ! true if weight derivatives are included [input]
      integer natoms
      logical grad, oactive(natoms)
c
      integer ictr
      integer ipol, nbf_ao, mbf_ao,  nq
      integer nmat, imat
      integer ibf(mbf_ao), iniz(natoms), ifin(natoms)
      integer cetobfr(2,*)
      integer calctype
c
c Scratch array for product of XC functional, basis functions, and 
c gradients of basis functions.
      double precision zzz(*)
c Scratch array for a density matrix from g_dens.
      double precision PPP(*)
c Used for storing part of PPP. 
      double precision F(*)
c GA Stuff
      integer g_dens(2)
c     
c     scratch array for hessian product
c     
      double precision hh(nq,3,*)
c     
c Basis functions and their gradients
c     
      double precision chi(nq,mbf_ao), delchi(nq,3,mbf_ao),
     1                 heschi(nq,6,mbf_ao)
c     
c XC-functional part 
c     
      double precision Amat(nq,ipol,nmat), Cmat(nq,3,ipol,nmat)
c     
c Force vector
c     
      double precision force(3,natoms)
c
c Local stuff
c
      logical ldo
c
      integer iat, n, ispin, mu,  nu, mu1, zzzptr
      integer ifinia, inizia, nnia, nbfia, ifirst, ilast
      integer nbhand
c
      double precision C1, C2, C3, z1, z2, z3, Am
      double precision out(3)
      double precision ccdel(nq,mbf_ao,*)
c
      if (grad) then
        do imat = 1, nmat
          do ispin = 1, ipol
            do nu = 1, mbf_ao
              do n = 1, nq
c delchi => 1 = x, 2 = y, 3 = z
                ccdel(n,nu,ispin) = 
     1            Cmat(n,1,ispin,imat)*delchi(n,1,nu) +
     2            Cmat(n,2,ispin,imat)*delchi(n,2,nu) +
     3            Cmat(n,3,ispin,imat)*delchi(n,3,nu)
              enddo ! n
            enddo ! nu
          enddo ! ispin
        enddo ! imat
      endif
      do iat = 1, natoms
        ldo = .true.
c ictr comes from grid_file_read, and it determines the atom center
c the grid points belong to.  We only do the gradients when the
c grid points don't belong to the current atom defined by iat when
c quadrature weight derivatives are needed.  The reason for doing
c this is clarified at the end of the code (employing Newton's third
c law to prevent double counting).
        if (ldew) ldo = ldo.and.(iat.ne.ictr)
        if (ldo) then
          if (iniz(iat).ne.0) then

c Here we select parts of the density matrix to grab.
            inizia = iniz(iat)
            ifinia = ifin(iat)
            ifirst = cetobfr(1,iat)
            ilast  = cetobfr(2,iat)
            nbfia = ilast-ifirst+1
            nnia = ifinia-inizia+1
c
c            write(6,*) 'inizia', inizia
c            write(6,*) 'ifinia', ifinia
c            write(6,*) 'ifirst', ifirst
c            write(6,*) 'ilast', ilast
c            write(6,*) 'nbfia', nbfia
c            write(6,*) 'nnia', nnia
c
c out stores each gradient
            out(1)=0.0d0
            out(2)=0.0d0
            out(3)=0.0d0
c Here, the code is modified because g_dens is dimensioned differently
c for response calculations.
            do imat = 1, nmat
              do ispin = 1, ipol
c We need both the ground state density matrix and the perturbed
c density matrix to get the contribution from the XC-functional to the
c TDDFT gradients.

c Use this for the ground state density matrix (for the gradient of
c the ground state density).  Note that g_dens(2) and g_dens(4) store
c the ground state density matrices in TDDFT gradients.
                if (calctype.eq.1) then
                  call ga_get(g_dens(2*imat+2*(ispin-1)), 
     1              ifirst, ilast, 1, nbf_ao, PPP, nbfia)
c Use this for the perturbed density matrix (for the gradient of the
c perturbed density).  Note that g_dens(1) and g_dens(3) store the
c penturbed density matrices in TDDFT gradients.
                else if (calctype.eq.2) then
                  call ga_get(g_dens(imat+2*(ispin-1)), 
     1              ifirst, ilast, 1, nbf_ao, PPP, nbfia)
                endif
                do mu = inizia, ifinia
c mu1 is used for indexing the scratch array hh (which starts at 1, 
c not inizia) 
                  mu1 = mu - inizia + 1
c First term of the XC-functional contribution to the ES gradients 
c (similar to the related term in the GS gradients, see Eq. 9 
c in jcp_98_5612).  This involves the XC-functional contribution
c from the partial derivative with respect to the the density.
                  do n = 1, nq
                    Am = Amat(n,ispin,imat)
c This is the x, y, z contribution to the gradient.
                    hh(n,1,mu1) =  Am*delchi(n,1,mu)
                    hh(n,2,mu1) =  Am*delchi(n,2,mu)
                    hh(n,3,mu1) =  Am*delchi(n,3,mu)
                  enddo
c Here we do the GGA only part, including the partial derivative of 
c the XC-functional with respect to the density gradient.
                  if (grad) then
                    do n = 1, nq
c C1 => dgaa, C2 => dgab, C3 => dgbb
                      C1 = Cmat(n,1,ispin,imat)
                      C2 = Cmat(n,2,ispin,imat)
                      C3 = Cmat(n,3,ispin,imat)
c Takes the LDA part and adds in the GGA parts, multiplied with the 
c Hessian of the basis function.
c heschi: 1 => dxdx    2 => dxdy    3 => dxdz 
c         4 => dydy    5 => dydz    6 => dzdz
                      hh(n,1,mu1) = hh(n,1,mu1) + C1*heschi(n,1,mu) +
     1                                            C2*heschi(n,2,mu) +
     2                                            C3*heschi(n,3,mu)
                      hh(n,2,mu1) = hh(n,2,mu1) + C1*heschi(n,2,mu) +
     1                                            C2*heschi(n,4,mu) +
     2                                            C3*heschi(n,5,mu)
                      hh(n,3,mu1) = hh(n,3,mu1) + C1*heschi(n,3,mu) +
     1                                            C2*heschi(n,5,mu) +
     2                                            C3*heschi(n,6,mu)
                    enddo ! n
                  endif ! grad
                enddo ! mu
c zzz => chi(nu)*Amat*delchi(mu) + chi(nu)*Cmat*heschi(mu)
c Left multiply by the basis functions.
                call dgemm('t', 'n', 3*nnia, mbf_ao, nq, 1.0d0, hh, nq,
     1               chi, nq, 0.0d0, zzz, 3*nnia)
c For a GGA, left multiply by the gradient of the basis functions.
                if (grad) then
c zzz => delchi(nu)*Cmat*delchi(mu) + zzz
c     => delchi(nu)*Cmat*delchi(mu) + chi(nu)*Amat*delchi(mu)
c      + chi(nu)*Cmat*heschi(mu) 
c 
                  call dgemm('t', 'n', 3*nnia, mbf_ao, nq, 1.0d0,
     1                 delchi(1,1,inizia), nq,
     2                 ccdel(1,1,ispin), nq, 1.0d0, zzz, 3*nnia)
                endif
c Gather elements from the density matrix in PPP into F.
                call gat_mat(F, PPP,  nbfia, nbf_ao,  nnia, mbf_ao,
     1            ifirst, 1, ibf(inizia), ibf)
c     zzz(3,nnia,mbf)*F(nnia,mbf)
                if (ipol.eq.1) then
                  call dgemv('n', 3, nnia*mbf_ao, -2.0d0, zzz, 3,
     1              F, 1, 1.0d0, out, 1)
                else
                  if (calctype.eq.1) then
                    call dgemv('n', 3, nnia*mbf_ao, -2.0d0, zzz, 3,
     1                F, 1, 1.0d0, out, 1)
                  else if (calctype.eq.2) then
                    call dgemv('n', 3, nnia*mbf_ao, -4.0d0, zzz, 3,
     1                F, 1, 1.0d0, out, 1)
                  endif
                endif
              enddo ! ispin
c iat => any atom that isn't the atom the grid points belong to
c This condition is for the gradient with respect to the nuclear 
c coordinates, without doing gradients of quadrature weights (i.e.
c doing the gradients with respect to nuclear coordinates explicitly).
              if (oactive(iat)) then
                force(1,iat) = force(1,iat) + out(1)
                force(2,iat) = force(2,iat) + out(2)
                force(3,iat) = force(3,iat) + out(3)
              endif
c ictr => the atom that the grid points belong to
c This condition is for the gradient with respect to the nuclear 
c coordinates, when we instead do the gradient with respect to the
c electronic coordinates (due to the definition of a Gaussian basis
c function, the sign becomes negative since the coordinates are
c always written Rnuc - rel, Rnuc:{Xnuc,Ynuc,Znuc} and 
c rel:{xel,yel,zel}).
              if (ldew.and.oactive(ictr)) then
                force(1,ictr) = force(1,ictr) - out(1)
                force(2,ictr) = force(2,ictr) - out(2)
                force(3,ictr) = force(3,ictr) - out(3)
              endif
            enddo ! imat
          endif ! iniz(iat).neq.0
        endif ! ldo
      enddo ! iat
c
      return 
      end 
