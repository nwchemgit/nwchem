      Subroutine hferi(Ep,Eq,R0,IJK,ERI,E3,sum,MPP,MPQ,
     &                 NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,MXD,
     &                 canAB,canCD,canPQ)
c $Id$

      Implicit none
#include "sh_order.fh"
      integer mpp,mpq,npp,npq,nint
      integer la,lb,lc,ld,lr,mxd
      integer ld2,lqmax,lqmax3
      integer iq,jq,nn,ia,ja,ka,mb_limit
      integer ij,mb,ib,jb,kb,ip,jp,kp,mp,ir,jr,lq,nq
      integer kr,nr,mq,mc,ic,jc,kc,md_limit,kl,md,id,jd,kd
      integer lb2,lc2,kq,ma,nc,la6,lb6,lc6,ld6
      integer nq1,iq1,iqfin,iqlast
      Logical canAB,canCD,canPQ

c--> Hermite Linear Expansion Coefficients

      Double precision Ep(3,NPP,0:MXD,0:(La+Lb),0:La,0:Lb)
      Double precision Eq(3,NPQ,0:MXD,0:(Lc+Ld),0:Lc,0:Ld)

c--> Auxiliary Function Integrals & Index

      Double precision R0(MPQ,MPP,*)
      Integer IJK(0:Lr,0:Lr,0:Lr)

c--> ERI

      Double precision ERI(Nint)

c--> Scratch Space

      Double precision E3(*),sum(MPQ,*),erinn,zot,
     ,     zot1,zot2
      Integer Nxyz(3),nqmx
c
c Compute electron repulsion integrals (ERI).
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
c N.B.  For simple ERI (i.e., no derivative integrals)  n[1-12] = 0!
c
c******************************************************************************
*rak:      integer num_pg, num_tot, num_qg
*rak:      save num_pg, num_tot, num_qg
*rak:      data num_pg /0/
*rak:      data num_qg /0/
*rak:      data num_tot /0/
*rak:
*rak:      if(MPP.ge.MPQ) then
*rak:        num_pg = num_pg + 1
*rak:      else
*rak:        num_qg = num_qg + 1
*rak:      endif
*rak:      num_tot = num_tot + 1
*rak:      if (num_tot.lt.3200.or.mod(num_tot,1000).eq.0) then
*rak:        write(6,*)' num_pg  ',num_pg
*rak:        write(6,*)' num_qg  ',num_qg
*rak:        write(6,*)' num_tot ',num_tot
*rak:      endif
c
c General case:  [ab|cd]

c Define the number of shell components on each center.

      Lb2 = ((Lb+1)*(Lb+2))/2
      Lc2 = ((Lc+1)*(Lc+2))/2
      Ld2 = ((Ld+1)*(Ld+2))/2
      la6=(la*(la+1)*(la+2))/6
      lb6=(lb*(lb+1)*(lb+2))/6
      lc6=(lc*(lc+1)*(lc+2))/6
      ld6=(ld*(ld+1)*(ld+2))/6
      md_limit = Ld2
       
c Initialize the block of ERIs.

      Lqmax  = Lc + Ld
      Lqmax3 = ((Lqmax+1)*(Lqmax+2)*(Lqmax+3))/6

c Loop over the components of the "A" and "B" shells.
      nqmx=0
      do Iq = 0,Lqmax
         do Jq = 0,Lqmax-Iq
!DEC$ LOOP COUNT MAX=10, MIN=1
            do Kq = 0,Lqmax-Iq-Jq
               nqmx = max(IJK(Iq,Jq,Kq),nqmx)
            enddo
         enddo
      enddo
      nn = 0
      do ma = 1,((La+1)*(La+2))/2

       nc =  la6 + ma
       ia = Ixyz(1,nc)
       ja = Ixyz(2,nc)
       ka = Ixyz(3,nc)
         

        if( canAB )then
           mb_limit = ma
           ij = (ma*(ma-1))/2 
        else
           mb_limit = Lb2
           ij = (ma-1)*Lb2 
        end if

        do mb = 1,mb_limit
           ij=ij+1


       nc =  lb6 + mb
       ib = Ixyz(1,nc)
       jb = Ixyz(2,nc)
       kb = Ixyz(3,nc)

c Sum across (Ip,Jp,Kp) for each value of (Iq,Jq,Kq).
          call dcopy(mpq*nqmx,0d0,0,sum,1)

          do Ip = 0,Ia+Ib
             do Jp = 0,Ja+Jb
                do Kp = 0,Ka+Kb
                   if(MPP.eq.1) then
                      E3(1) = Ep(1,1,0,Ip,Ia,Ib)*
     &                     Ep(2,1,0,Jp,Ja,Jb)*
     &                     Ep(3,1,0,Kp,Ka,Kb)
                   else
                      
c Define the product of the Hermite expansions coefficients for 
c overlap distribution "P".
!DEC$ LOOP COUNT MAX=30, MIN=1
#ifdef PGI_NOSIMD
cpgi$l novect
#endif
                      do mp = 1,MPP
                         E3(mp) = Ep(1,mp,0,Ip,Ia,Ib)*
     &                        Ep(2,mp,0,Jp,Ja,Jb)*
     &                        Ep(3,mp,0,Kp,Ka,Kb)
                      end do
                   endif

                do Iq = 0,Lqmax
                   Ir = Ip + Iq
                  do Jq = 0,Lqmax-Iq
                     Jr = Jp + Jq
                      Lq = Iq + Jq -1
                    do Kq = 0,Lqmax-Iq-Jq
                       
                      nq = IJK(Iq,Jq,Kq)

                      Kr = Kp + Kq

                      nr = IJK(Ir,Jr,Kr)

c Include the factor of (-1)**(Iq+Jq+Kq).
                      Lq=Lq+1
#if defined(GCC4) || defined(PGLINUX)
                      if(IAND(Lq,1).eq.1)then
#else
                      if(AND(Lq,1).eq.1)then
#endif
                        if(MPQ.eq.1) then
!DEC$ LOOP COUNT MAX=30, MIN=1
                           do mp = 1,MPP
                              sum(1,nq) = sum(1,nq)-E3(mp)*R0(1,mp,nr)
                           end do
                        else
                        do mp = 1,MPP
                           zot=-E3(mp)
!DEC$ LOOP COUNT MAX=30, MIN=1
                          do mq = 1,MPQ
                            sum(mq,nq) = sum(mq,nq)+zot*R0(mq,mp,nr)
                          end do
                        end do
                        endif
                      else
                        if(MPQ.eq.1) then
!DEC$ LOOP COUNT MAX=30, MIN=1
                           do mp = 1,MPP
                              sum(1,nq) = sum(1,nq)+E3(mp)*R0(1,mp,nr)
                           end do
                        else
                        do mp = 1,MPP
                           zot=E3(mp)
!DEC$ LOOP COUNT MAX=30, MIN=1
                          do mq = 1,MPQ
                            sum(mq,nq) = sum(mq,nq)+zot*R0(mq,mp,nr)
                          end do
                        end do
                        endif
                      end if

                    end do
                  end do
                end do

              end do
            end do
          end do

c Loop over the components of the "C" and "D" shells.

          do mc = 1,Lc2

             nc =  lc6 + mc
             ic = Ixyz(1,nc)
             jc = Ixyz(2,nc)
             kc = Ixyz(3,nc)

            if( canCD ) md_limit = mc


            if( canAB )then
               kl = (mc*(mc-1))/2 
            else
               kl = (mc-1)*Ld2 
            end if
            do md = 1,md_limit

              if( canPQ )then
                 kl=kl+1
                 if( kl.gt.ij ) go to 480
              end if

             nc =  ld6 + md
             id = Ixyz(1,nc)
             jd = Ixyz(2,nc)
             kd = Ixyz(3,nc)

              nn = nn + 1

c Sum across (Iq,Jq,Kq).
c Define the product of the Hermite expansion coefficients for 
c overlap distribution "Q" and calculate eri
              erinn=0d0
                
              if(MPQ.eq.1) then
                 do Iq = 0,Ic+Id
                    do Jq = 0,Jc+Jd
!DEC$ LOOP COUNT MAX=30, MIN=1
                       do Kq = 0,Kc+Kd
                          nq = IJK(Iq,Jq,Kq)
                          ERInn = ERInn + 
     &                         Eq(1,1,0,Iq,Ic,Id)*
     &                         Eq(2,1,0,Jq,Jc,Jd)*
     &                         Eq(3,1,0,Kq,Kc,Kd)*sum(1,nq)
                       enddo
                    enddo
                 enddo
              else
                 iqfin=0
                 iqlast=ic+id
                 if(iqlast.gt.0) then
                    iqfin=-1
#if defined(GCC4) || defined(PGLINUX)
                    if(iand(iqlast,1).eq.0) then
#else
                    if(and(iqlast,1).eq.0) then
#endif
                       iqfin=ic+id
                       iqlast=iqlast-2
                    endif
                    do Iq = 0,iqlast,2
                       iq1=iq+1
                       do Jq = 0,Jc+Jd
                          do Kq = 0,Kc+Kd
                             nq = IJK(Iq,Jq,Kq)
                             nq1 = IJK(Iq1,Jq,Kq)
!DEC$ LOOP COUNT MAX=30, MIN=1
                             do mq = 1,MPQ
                                ERInn = ERInn + 
     &                               (Eq(1,mq,0,Iq,Ic,Id)*sum(mq,nq)+
     &                               Eq(1,mq,0,Iq1,Ic,Id)*sum(mq,nq1))*
     &                               Eq(2,mq,0,Jq,Jc,Jd)*
     &                               Eq(3,mq,0,Kq,Kc,Kd)
                             end do
                          end do
                       enddo
                    enddo
                 endif
                 if(iqfin.ne.-1) then
                    do Jq = 0,Jc+Jd
                       do Kq = 0,Kc+Kd
                          nq = IJK(Iqfin,Jq,Kq)
!DEC$ LOOP COUNT MAX=30, MIN=1
                          do mq = 1,MPQ
                             ERInn = ERInn + 
     &                            Eq(1,mq,0,Iqfin,Ic,Id)*
     &                            Eq(2,mq,0,Jq,Jc,Jd)*
     &                            Eq(3,mq,0,Kq,Kc,Kd)*sum(mq,nq)
                          end do
                       end do
                    enddo
                 endif
              endif
              eri(nn)=erinn

            end do              ! md
  480       continue
          end do                ! mc
        end do                  ! mb
      end do                    ! ma

      end
