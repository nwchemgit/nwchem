      Subroutine xc_rhogen(TOLL,
     &     basis,iga_dens,
     &     max_at_bf,natoms,ipol,nq,nbf,mbf,
     &     ldelley,GRAD,
     &     F,Pmat,ff,zz,
     &     chi,delchi,ibf,iniz,ifin,
     &     rho,delrho,gamma,
     &     rhoat)
c
c$Id: xc_rhogen.F,v 1.5 1995-10-17 08:55:36 d3g681 Exp $
c
      Implicit double precision (a-h,o-z)
      Implicit integer (i-n)

#include "bas.fh"
      Logical GRAD,ldelley
      Integer basis

    
c--> Density Matrix
      Dimension iga_dens(*)
      Dimension Pmat(max_at_bf*max_at_bf),F(max_at_bf*max_at_bf)
      Dimension ibf(mbf),iniz(natoms),ifin(natoms)

c--> Basis Functions & Gradients

      Dimension chi(nq,mbf),delchi(3,nq,mbf)

c--> Sampling Matrices for [Products of] Basis Functions & Gradients

      Dimension zz(3,nq,max_at_bf)
      Dimension ff(nq,max_at_bf)

c--> Charge Density & Its Cube Root

      Dimension rho(ipol*(ipol+1)/2,nq),rhoat(nq,natoms)

c--> Charge Density Gradient

      Dimension delrho(ipol,3,nq)

c--> Charge Density Gradient Dot Products

      Dimension gamma(ipol*(ipol+1)/2,nq)
c
c Evaluate the charge density and its gradient at each of the sampling 
c points and then compute its cube root.
c
c***************************************************************************
c
c     tolerance setting
c
      rho_tol0=TOLL*.1d0
      rho_tol1=0.5D0*rho_tol0

      i0=ipol-1

c     Loop over products of AO basis functions to define the matrix bands.

      FUNC_MAX=dabsmax(nq*mbf,chi)

      do 230 iat = 1,natoms
        inizia=iniz(iat)
        if(inizia.ne.0) then
          if (.not.bas_ce2bfr(basis,iat,ifirst,ilast))
     &         call errquit('Exiting in rho_gen',13)
          nbfia=ilast-ifirst+1
          nnia=ifin(iat)-inizia+1

          FUNC_MAXI=dabsmax(nnia*nq,chi(1,inizia))
          FUNC_MAX2=FUNC_MAXI*FUNC_MAX

          do jat=1,iat
            if(iniz(jat).ne.0) then
              inizja=iniz(jat)
            if (.not.bas_ce2bfr(basis,jat,jfirst,jlast))
     &           call errquit('Exiting in rho_gen',14)
              do 215 ii=1,ipol

                call get_atom_block(iga_dens(ii), basis,
     $               iat, jat, Pmat, idim, jdim)
                
                P_MAX=dabsmax(idim*jdim,Pmat)
                if((P_MAX*FUNC_MAX2).lt.rho_tol0) goto 215

                nbfja=jlast-jfirst+1
                nnja=ifin(jat)-inizja+1

                call gat_mat(F,Pmat,nbfia,nbfja,nnia,nnja,
     &               ifirst,jfirst,ibf(inizia),ibf(inizja))

                call dgemm('n','n',nq,nnja,nnia,1.d0,chi(1,inizia),nq,
     &               F,nnia,0.d0,ff,nq)

                 if(iat.ne.jat) 
     &               call dscal(nq*nnja,2.d0,ff,1)
                 if(ldelley) then
                   if(iat.ne.jat) then
                     do mu=inizja,ifin(jat)
                       mu1=mu-inizja+1
                       do n=1,nq
                         aaa= chi(n,mu)*ff(n,mu1)
                         rho(i0+ii,n) = rho(i0+ii,n) + aaa
                         rhoat(n,iat)=rhoat(n,iat)+ aaa
                         rhoat(n,jat)=rhoat(n,jat)+ aaa
                       enddo
                     enddo
                   else
                     do mu=inizja,ifin(jat)
                       mu1=mu-inizja+1
                       do n=1,nq
                         aaa= chi(n,mu)*ff(n,mu1)
                         rho(i0+ii,n) = rho(i0+ii,n) + aaa
                         rhoat(n,iat)=rhoat(n,iat)+ aaa
                       enddo
                     enddo
                   endif
                 else
                   do mu=inizja,ifin(jat)
                     mu1=mu-inizja+1
cdir$ UNROLL 3
                     do n=1,nq
                       rho(i0+ii,n) = rho(i0+ii,n) + 
     +                      chi(n,mu)*ff(n,mu1)
                     enddo
                   enddo
                 endif

                 if(grad)then
                   call dgemm('n','n',3*nq,nnja,nnia,1.d0,
     &                  delchi(1,1,inizia),3*nq,
     &                  F,nnia,0.d0,zz,3*nq)
                 if(iat.ne.jat) 
     &               call dscal(3*nq*nnja,2.d0,zz,1)
                 do mu=inizja,ifin(jat)
                   mu1=mu-inizja+1
                   do n=1,nq
                     delrho(ii,1,n) = delrho(ii,1,n) + 
     +                    delchi(1,n,mu)*ff(n,mu1) +
     +                     chi(n,mu)*zz(1,n,mu1)
                     delrho(ii,2,n) = delrho(ii,2,n) + 
     +                    delchi(2,n,mu)*ff(n,mu1) +
     +                     chi(n,mu)*zz(2,n,mu1)
                     delrho(ii,3,n) = delrho(ii,3,n) + 
     +                    delchi(3,n,mu)*ff(n,mu1) +
     +                     chi(n,mu)*zz(3,n,mu1)
                   enddo
                 enddo
               endif
                
  215         continue
            endif
          enddo

        endif

  230 continue
      if(ipol.eq.2) then
        call cut_tol(nq,rho(3,1),3,rho_tol1)
        call cut_tol(nq,rho(2,1),3,rho_tol1)
        do 240 n = 1,nq

          rho(1,n) = rho(2,n) + rho(3,n)

  240   continue
        if(GRAD) then
          do 250 n = 1,nq
            
            gamma(1,n) = delrho(1,1,n)*delrho(1,1,n) +
     &           delrho(1,2,n)*delrho(1,2,n) +
     &           delrho(1,3,n)*delrho(1,3,n)
            
            gamma(2,n) = delrho(2,1,n)*delrho(2,1,n) +
     &           delrho(2,2,n)*delrho(2,2,n) +
     &           delrho(2,3,n)*delrho(2,3,n)
            
            gamma(3,n) = delrho(1,1,n)*delrho(2,1,n) +
     &           delrho(1,2,n)*delrho(2,2,n) +
     &           delrho(1,3,n)*delrho(2,3,n)
            
  250     continue
        endif
      else
        call cut_tol(nq,rho,1,rho_tol0)
        if(GRAD) then

          do  n = 1,nq
            gamma(1,n) = delrho(1,1,n)*delrho(1,1,n)
          enddo
          do  n = 1,nq
            gamma(1,n) = gamma(1,n)+delrho(1,2,n)*delrho(1,2,n)
          enddo
          do  n = 1,nq
            gamma(1,n) = gamma(1,n)+delrho(1,3,n)*delrho(1,3,n)
          enddo
        endif
      endif
      if(ldelley) 
     &     call cut_tol(nq*natoms,rhoat,1,rho_tol0)
      return
      end
      
      Subroutine cut_tol(n,vec,inc,tol)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      dimension vec(n)
      i=1
cdir$ UNROLL 3
      do ii=1,n
        if(vec(i).lt.tol) vec(i)=tol
        i=i+inc
      enddo
      return
      end
