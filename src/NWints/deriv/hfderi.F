      Subroutine hfderi(Ep,Eq,R0,IJK,ERI,E3,sum,
     &                  NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &                  MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &                  canAB,canCD,canPQ)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
#include "sh_order.fh"

      Logical canAB,canCD,canPQ

c--> Hermite Linear Expansion Coefficients

      Dimension Ep(3,NPP,0:MXD,0:(La+Lb),0:La,0:Lb)
      Dimension Eq(3,NPQ,0:MXD,0:(Lc+Ld),0:Lc,0:Ld)

c--> Auxiliary Function Integrals & Index

      Dimension R0((NPQ*NPP),*),IJK(0:Lr,0:Lr,0:Lr)

c--> ERI

      Dimension ERI((NPQ*NPP),Nint,(Nder1*Nder2*Nder3))

c--> Derivative Indices

      Dimension inder1(6,Nder1),inder2(3,Nder2),inder3(3,Nder3)

c--> Scratch Space

      Dimension E3(*),sum((NPQ*NPP),(Nder1*Nder2),*)
c
c Compute gradient electron repulsion integrals (ERI).
c
c     Formula:
c
c               __
c               \     Ic,Id;n10   Jc,Jd;n11   Kc,Kd;n12
c     ERI  =    /   Ex          Ey          Ez          SUM
c               --    Iq          Jq          Kq           Iq,Jq,Kq
c            Iq,Jq,Kq
c
c                            __
c                           \     Lq   Ia,Ib;n7   Ja,Jb;n8   Ka,Kb;n9
c         SUM          =    /  (-1)  Ex         Ey         Ez         R
c            Iq,Jq,Kq       --         Ip         Jp         Kp        Ir,Jr,Kr
c                        Ip,Jp,Kp
c
c                                Ir = (Ip+n1) + (Iq+n4)
c                        where   Jr = (Jp+n2) + (Jq+n5)
c                                Kr = (Kp+n3) + (Kq+n6)
c
c                         and    Lq = (Iq+n4) + (Jq+n5) + (Kq+n6)
c
c Currently, the 12 gradient integrals can be grouped into 3 sets. The 
c integrals in each set can be evaluated with a single call to this routine
c The derivative indices specify the order of differentiation for each of 
c the 12 coordinates:
c
c     Set 1, inder1(i,n), i=1-3:   d/dPx, d/dPy, d/dPz;       n1,  n2,  n3
c                         i=4-6:   d/dQx, d/dQy, d/dQz;       n4,  n5,  n6
c     Set 2, inder2(i,n), i=1-3:   d/dRx, d/dRy, d/dRz;       n7,  n8,  n9
c     Set 3, inder3(i,n), i=1-3:   d/dSx, d/dSy, d/dSz;      n10, n11, n12
c
c N.B. This routine cannot simultaneously evaluate integrals from
c      different sets!
c
c******************************************************************************
      double precision zot
      integer la6,lb6,lc6,ld6
      integer nder12
      logical lneg
      la6=(la*(la+1)*(la+2))/6
      lb6=(lb*(lb+1)*(lb+2))/6
      lc6=(lc*(lc+1)*(lc+2))/6
      ld6=(ld*(ld+1)*(ld+2))/6
      nder12=nder1*nder2
      if(nder12.ne.1.and.nder12.ne.3) call
     C     errquit(' hfderi: nder12.ne.1.and.ne.3 ',nder12,0)
    
c

c Initialize the block of ERIs.

      call dcopy(nder12*Nder3*Nint*NPQ*NPP,0d0,0,ERI,1)

c General case:  [ab|cd]

c Define the number of shell components on each center.

       La2 = ((La+1)*(La+2))/2
       Lb2 = ((Lb+1)*(Lb+2))/2
       Lc2 = ((Lc+1)*(Lc+2))/2
       Ld2 = ((Ld+1)*(Ld+2))/2

       Lqmax  = Lc + Ld
       Lqmax3 = ((Lqmax+1)*(Lqmax+2)*(Lqmax+3))/6

c Loop over the components of the "a" and "b" shells.

       nn = 0
       md_limit = Ld2
       mb_limit = Lb2
       
       do 500 ma = 1,La2

        nc =  la6 + ma
        ia = Ixyz(1,nc)
        ja = Ixyz(2,nc)
        ka = Ixyz(3,nc)

        if( canAB )  mb_limit = ma

        do 490 mb = 1,mb_limit

           if( canPQ )then

            if( canAB )then
             ij = (ma*(ma-1))/2 + mb
            else
             ij = (ma-1)*Lb2 + mb
            end if


           end if


           nc =  lb6 + mb
           ib = Ixyz(1,nc)
           jb = Ixyz(2,nc)
           kb = Ixyz(3,nc)

c Sum across (Ip,Jp,Kp) for each value of (Iq,Jq,Kq).

         call dcopy((NPQ*NPP)*nder12*Lqmax3,0d0,0,sum,1)

         do 425 Ip = 0,Ia+Ib
         do 425 Jp = 0,Ja+Jb
         do 425 Kp = 0,Ka+Kb

          do 424 id2 = 1,Nder2

           n7 = inder2(1,id2)
           n8 = inder2(2,id2)
           n9 = inder2(3,id2)

c Define the product of the Hermite expansions coefficients for 
c charge distribution "P".

           do 410 mp = 1,NPP
            E3(mp) = Ep(1,mp,n7,Ip,Ia,Ib)*
     &               Ep(2,mp,n8,Jp,Ja,Jb)*
     &               Ep(3,mp,n9,Kp,Ka,Kb)
  410      continue

           do 420 Iq = 0,Lqmax
           do 420 Jq = 0,Lqmax-Iq
           do 420 Kq = 0,Lqmax-Iq-Jq

            nq = IJK(Iq,Jq,Kq)

             id4 = (id2-1)*Nder1 

            do 419 id1 = 1,Nder1

             n1 = inder1(1,id1)
             n2 = inder1(2,id1)
             n3 = inder1(3,id1)

             n4 = inder1(4,id1)
             n5 = inder1(5,id1)
             n6 = inder1(6,id1)

             Ir = (Ip+Iq)+n1+n4
             Jr = (Jp+Jq)+n2+n5
             Kr = (Kp+Kq)+n3+n6

             nr = IJK(Ir,Jr,Kr)

             id4 = id4 + 1

c Include the factor of (-1)**Lq

             Lq = (Iq+n4) + (Jq+n5) + (Kq+n6)

             lneg=( and(Lq,1).eq.1 )

              mr = 0
              if(npq.eq.1) then
                 if(lneg) then
                    do  mp = 1,NPP
                       sum(mp,id4,nq) = sum(mp,id4,nq) -e3(mp)*R0(mp,nr)
                    enddo
                 else
                    do  mp = 1,NPP
                       sum(mp,id4,nq) = sum(mp,id4,nq) +e3(mp)*R0(mp,nr)
                    enddo
                 endif
              else
              do  mp = 1,NPP
                 if(lneg) then
                    zot=-e3(mp)
                 else
                    zot=e3(mp)
                 endif
              do mq = mr+1,mr+NPQ

               sum(mq,id4,nq) = sum(mq,id4,nq) +zot*R0(mq,nr)

               enddo
               mr = mr + NPQ
            enddo
            endif


  419       continue

  420      continue

  424     continue

  425    continue

c Loop over the components of the "c" and "d" shells.

         do 480 mc = 1,Lc2

            nc =  lc6 + mc
            ic = Ixyz(1,nc)
            jc = Ixyz(2,nc)
            kc = Ixyz(3,nc)

          if( canCD ) md_limit = mc

          do 470 md = 1,md_limit

           if( canPQ )then

            if( canAB )then
             kl = (mc*(mc-1))/2 + md
            else
             kl = (mc-1)*Ld2 + md
            end if

            if( kl.gt.ij ) go to 480

           end if

             nc =  ld6 + md
             id = Ixyz(1,nc)
             jd = Ixyz(2,nc)
             kd = Ixyz(3,nc)

           nn = nn + 1

c Sum across (Iq,Jq,Kq).

           do 450 Iq = 0,Ic+Id
           do 450 Jq = 0,Jc+Jd
           do 450 Kq = 0,Kc+Kd

            nq = IJK(Iq,Jq,Kq)

            id5 = 0
            if(npq.eq.1) then
               do  id3 = 1,Nder3

                  n10 = inder3(1,id3)
                  n11 = inder3(2,id3)
                  n12 = inder3(3,id3)
                  zot=
     *                       Eq(1,1,n10,Iq,Ic,Id)*
     &                       Eq(2,1,n11,Jq,Jc,Jd)*
     &                       Eq(3,1,n12,Kq,Kc,Kd)
                  if(nder12.eq.1) then
                     id5 = id5 + 1
                     
                     do  mp = 1,NPP

                        ERI(mp,nn,id5) = ERI(mp,nn,id5) + 
     +                       sum(mp,1,nq)*zot
                        
                     enddo
                  else
                     do  mp = 1,NPP
                        ERI(mp,nn,id5+1) = ERI(mp,nn,id5+1) + 
     +                       sum(mp,1,nq)*zot
                        ERI(mp,nn,id5+2) = ERI(mp,nn,id5+2) + 
     +                       sum(mp,2,nq)*zot
                        ERI(mp,nn,id5+3) = ERI(mp,nn,id5+3) + 
     +                       sum(mp,3,nq)*zot
                     enddo
                     id5 = id5 + 3
                  endif
               enddo
            else
            do 449 id3 = 1,Nder3

             n10 = inder3(1,id3)
             n11 = inder3(2,id3)
             n12 = inder3(3,id3)

c Define the product of the Hermite expansions coefficients for 
c charge distribution "Q".

             if(NPQ.eq.1) then
                zot = Eq(1,1,n10,Iq,Ic,Id)*
     &               Eq(2,1,n11,Jq,Jc,Jd)*
     &               Eq(3,1,n12,Kq,Kc,Kd)
             else
             do mq = 1,NPQ
                E3(mq) = Eq(1,mq,n10,Iq,Ic,Id)*
     &               Eq(2,mq,n11,Jq,Jc,Jd)*
     &               Eq(3,mq,n12,Kq,Kc,Kd)
             enddo
             endif

             if(NPP.eq.1) then
                if(NPQ.eq.1) then
                   if(nder12.eq.1) then
                      ERI(1,nn,id5) = ERI(1,nn,id5) + 
     +                     zot*sum(1,1,nq)
                      id5=id5+1
                   else
                      ERI(1,nn,id5+1) = ERI(1,nn,id5) + 
     +                     zot*sum(1,1,nq)
                      ERI(1,nn,id5+2) = ERI(1,nn,id5) + 
     +                     zot*sum(1,2,nq)
                      ERI(1,nn,id5+3) = ERI(1,nn,id5) + 
     +                     zot*sum(1,3,nq)
                      id5 = id5 + 3
                   endif
                else
                   do id4 = 1,nder12
                      id5 = id5 + 1
                      do  mq = 1,NPQ
                         ERI(mq,nn,id5) = ERI(mq,nn,id5) + 
     +                        E3(mq)*sum(mq,id4,nq)
                      enddo
                   enddo
                endif
             else
                mr = 0
                if(nder12.eq.1) then
              id5 = id5 + 1
                   do mp = 1,NPP
                      do  mq = 1,NPQ
                         mr = mr + 1

                         ERI(mr,nn,id5) = ERI(mr,nn,id5) + 
     +                        E3(mq)*sum(mr,1,nq)

                    enddo
                 enddo
                else

              do 440 mp = 1,NPP
                 do 440 mq = 1,NPQ
                    mr = mr + 1

                    ERI(mr,nn,id5+1) = ERI(mr,nn,id5+1) + 
     +                   E3(mq)*sum(mr,1,nq)
                    ERI(mr,nn,id5+2) = ERI(mr,nn,id5+2) + 
     +                   E3(mq)*sum(mr,2,nq)
                    ERI(mr,nn,id5+3) = ERI(mr,nn,id5+3) + 
     +                   E3(mq)*sum(mr,3,nq)

 440             continue
              id5 = id5 + 3
                 
          endif
             endif

  449       continue

            endif
  450      continue

  470     continue
  480    continue
  490   continue
  500  continue

      end
