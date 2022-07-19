      subroutine dftg_quadv0_b_adft(ictr,ldew,basis,ipol,grad,
     $          kslap,nbf,mbf,nat,nq,chi,delchi,heschi,d3chi,
     $          Amat,Cmat,Lmat,scr,ibf,iniz,ifin,cetobfr,
     $          force,oactive,cd_coef)

      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "dftpara.fh"
      logical ldew
      logical grad
      logical kslap
      integer ictr
      integer basis
      integer ipol
      integer nbf
      integer mbf
      integer nat
      integer nq
      logical oactive(nat)
      integer ibf(mbf), iniz(nat), ifin(nat)
      integer cetobfr(2,*)
      double precision chi(nq,mbf),delchi(nq,3,mbf)
      double precision heschi(nq,6,mbf),d3chi(nq,10,mbf)
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
      double precision Lmat(nq,ipol)
      double precision force(3,nat)
      double precision scr(nq,10)
      double precision cd_coef(nbf,ipol)
      double precision gradient(3)
      double precision coeff(nbf)
      double precision ddot
      external ddot

      logical ldo
      integer iat, n, ispin, mu, ixyz
      integer ifinia, inizia, nnia, nbfia, ifirst, ilast
      
      do iat=1,nat
        ldo = .true.
        if (ldew) ldo = ldo .and. (iat.ne.ictr)
        if (.not.ldo) cycle
        ldo = (oactive(iat).or.(ldew .and. oactive(ictr)))
        if (.not.ldo) cycle
        if (iniz(iat).eq.0) cycle

        inizia = iniz(iat)
        ifinia = ifin(iat)
        ifirst = cetobfr(1,iat)
        ilast  = cetobfr(2,iat)
        nbfia  = ilast - ifirst + 1
        nnia   = ifinia - inizia + 1

        gradient(:) = 0d0

        do ispin=1,ipol

          do mu=inizia,ifinia
            coeff(mu-inizia+1) = cd_coef(ibf(mu),ispin)
          enddo

          call dgemv('n',3*nq,nnia,1d0,delchi(1,1,inizia),3*nq,coeff,1,
     $                0d0,scr,1)     
          call dgemv('t',nq,3,1d0,scr,nq,amat(1,ispin),1,1d0,
     $                gradient,1)

          if (grad) then
            call dgemv('n',6*nq,nnia,1d0,heschi(1,1,inizia),6*nq,coeff,
     $                  1,0d0,scr,1)
            gradient(1) = gradient(1) + 
     $           ddot(3*nq,scr,1,cmat(1,1,ispin),1)
            gradient(2) = gradient(2) +
     $           ddot(nq,scr(1,2),1,cmat(1,1,ispin),1)       
            gradient(2) = gradient(2) +
     $           ddot(2*nq,scr(1,4),1,cmat(1,2,ispin),1)       
            gradient(3) = gradient(3) +
     $           ddot(nq,scr(1,3),1,cmat(1,1,ispin),1)
            gradient(3) = gradient(3) +
     $           ddot(2*nq,scr(1,5),1,cmat(1,2,ispin),1)
          endif

          if (kslap) then
            call dgemv('n',10*nq,nnia,1d0,d3chi(1,1,inizia),10*nq,
     $                  coeff,1,0d0,scr,1)       
            scr(:,1) = scr(:,1) + scr(:,4) + scr(:,6)
            gradient(1) = gradient(1) + ddot(nq,scr,1,lmat(1,ispin),1)
            scr(:,1) = scr(:,2) + scr(:,7) + scr(:,9)
            gradient(2) = gradient(2) + ddot(nq,scr,1,lmat(1,ispin),1)
            scr(:,1) = scr(:,3) + scr(:,8) + scr(:,10)
            gradient(3) = gradient(3) + ddot(nq,scr,1,lmat(1,ispin),1)
          endif
        enddo

        if (oactive(iat)) then
          force(:,iat) = force(:,iat)-gradient(:)
        endif
        if (ldew.and.oactive(ictr)) then
          force(:,ictr) = force(:,ictr)+gradient(:)
        endif
      enddo

      return
      end
