      Subroutine xc_tabcd(TOLL,
     &     Tmat,TTmat,Amat,Bmat,Cmat,Dmat,chi,delchi,
     &     ipol,nq,nbf,max_at_bf,max_at_bf2,GRAD,basis,
     &     natoms,iniz,ifin,g_vxc,ibf)

C$Id: xc_tabcd.F,v 1.5 1996-05-04 21:08:39 og845 Exp $
      implicit none
#include "bas.fh"
      Logical GRAD
      integer basis
      
      integer max_at_bf ! [input]
      integer max_at_bf2 ! [input]
      integer ipol  ! [input]
      integer nq    ! [input]
      integer nbf    ! [input]
      integer natoms ! [input]
      double precision TOLL
      integer g_vxc(*)   ! GA handle [input]

c-->

      double precision Tmat(max_at_bf2),TTmat(max_at_bf2)

c--> Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)

c--> Sampling Matrices for [Products of] Basis Functions & Gradients

      double precision Bmat(nq,max_at_bf)
      double precision Dmat(3,nq,max_at_bf)
      integer iniz(natoms),ifin(natoms)

c--> Basis Functions & Gradients

      double precision chi(nq,nbf),delchi(3,nq,nbf)

      integer ibf(nbf)
c
      double precision A_MAX,FUNC_MAX,AXF_MAX,ATOM_MAX,
     ,     B_MAX,FUNC_MAXJ
      integer iat,inizia,ifinia,nbfia,nnia,ifirst,ilast
      integer jat,inizja,ifinja,nbfja,nnja,jfirst,jlast
      integer ii,mu,mu1
      integer n
      double precision chi1
      double precision dabsmax
      external dabsmax
c
c Compute the matrix product for the XC potential and energy:
c
c              T = transpose(A*B) + transpose(C*D)
c
c
c******************************************************************************

c
      A_MAX=dabsmax(nq*ipol,Amat)
      FUNC_MAX=dabsmax(nq*nbf,chi)
      AXF_MAX=A_MAX*FUNC_MAX
      do 430 iat = 1,natoms
        if(iniz(iat).ne.0) then
          if (.not.bas_ce2bfr(basis,iat,ifirst,ilast))
     &         call errquit('Exiting in tabcd',43)

          inizia=iniz(iat)
          ifinia=ifin(iat)
          nnia=ifinia-inizia+1
          ATOM_MAX=dabsmax(nnia*nq,chi(1,inizia))
          if((AXF_MAX*ATOM_MAX).lt.TOLL) goto 430
          nbfia=ilast-ifirst+1
          do 405 ii=1,ipol
            do mu=1,nnia
              mu1=mu+inizia-1
              if(GRAD) then
                do n=1,nq
                  chi1=chi(n,mu1)
                  Bmat(n,mu) = Amat(ii,n)*chi1 +
     +                 delchi(1,n,mu1)*Cmat(ii,1,n) +
     +                 delchi(2,n,mu1)*Cmat(ii,2,n) +
     +                 delchi(3,n,mu1)*Cmat(ii,3,n)
                  Dmat(1,n,mu) = Cmat(ii,1,n)*chi1
                  Dmat(2,n,mu) = Cmat(ii,2,n)*chi1
                  Dmat(3,n,mu) = Cmat(ii,3,n)*chi1
                enddo
              else

cdir$ UNROLL 3
                do n=1,nq
                  Bmat(n,mu) = chi(n,mu1)*Amat(ii,n)
                enddo
              endif
            enddo
          B_MAX=dabsmax(nnia*nq,Bmat)

            do 168 jat=1,iat
              if(iniz(jat).ne.0) then
                inizja=iniz(jat)
                ifinja=ifin(jat)
                  FUNC_MAXJ=dabsmax(
     &            (ifinja-inizja+1)*nq,chi(1,inizja))
                  if(B_MAX*FUNC_MAXJ.lt.TOLL) goto 168
                if (.not.bas_ce2bfr(basis,jat,jfirst,jlast))
     &               call errquit('Exiting in tabcd',44)
                nbfja=jlast-jfirst+1
                nnja=ifinja-inizja+1
                call dgemm('t','n',nnia,nnja,nq,1.d0,Bmat,nq,
     &               chi(1,inizja),nq,0.d0,Tmat,nnia)
                if (GRAD) call dgemm(
     &               't','n',nnia,nnja,3*nq,1.d0,Dmat,3*nq,
     &               delchi(1,1,inizja),3*nq,1.d0,Tmat,nnia)
                call scat_mat(TTmat,Tmat,nbfia,nbfja,nnia,nnja,
     &               ifirst,jfirst,ibf(inizia),ibf(inizja))
                call upd_atom_block(g_vxc(ii), basis,
     $               iat, jat, TTmat)
                if(iat.ne.jat) then
                  call transp_mat(TTmat,Tmat,nbfia,nbfja)
                  call upd_atom_block(g_vxc(ii), basis,
     $                 jat, iat, Tmat)
                endif

              endif
  168       continue
  405     continue

        endif

  430 continue


      return
      end


