c $Id$
      Subroutine hferi_gen(Ep,Eq,PairP,PairQ,R0,IJK,ERI,E3,t1,t2,t3,t4,
     &    MPP,MPQ,NPP,NPQ,La,Lb,Lc,Ld,La2,Lb2,Lc2,Ld2,Lqmax,Lqmax3,Lr,
     &    Acoefs,Bcoefs,Ccoefs,Dcoefs,NPA,NPB,NPC,NPD,NCA,NCB,NCC,NCD,
     &    MXD,canAB,canCD,canPQ)

      Implicit none
#include "errquit.fh"

      Logical canAB,canCD,canPQ

c--> Number of pairs in pass, total number

      Integer MPP,MPQ, NPP,NPQ

c--> Angular momenta and number of cartesian components

      Integer La,Lb,Lc,Ld, La2,Lb2,Lc2,Ld2, Lqmax,Lqmax3,Lr

c--> Derivative order

      Integer MXD

c--> Number of primitives, number of contracted functions

      Integer NPA,NPB,NPC,NPD, NCA,NCB,NCC,NCD

c--> Hermite Linear Expansion Coefficients

      Double precision Ep(3,NPP,0:MXD,0:(La+Lb),0:La,0:Lb)
      Double precision Eq(3,NPQ,0:MXD,0:(Lc+Ld),0:Lc,0:Ld)

c--> Pair index arrays

      Integer PairP(2,NPP)
      Integer PairQ(2,NPQ)

c--> Auxiliary Function Integrals & Index

      Double precision R0(MPQ,MPP,*)
      Integer IJK(0:Lr,0:Lr,0:Lr)

c--> Contraction coefficients

      Double precision Acoefs(NPA,NCA)
      Double precision Bcoefs(NPB,NCB)
      Double precision Ccoefs(NPC,NCC)
      Double precision Dcoefs(NPD,NCD)

c--> ERI

      Double precision ERI(Ld2,NCD,Lc2,NCC,Lb2*NCB,La2*NCA)

c--> Scratch Space

      Double precision E3(*)
      Double precision t1(MPQ,Lqmax3,NPB,NCA)  ! 1st 1/4 trans
      Double precision t2(MPQ,Lqmax3,MPP)      ! 2nd 1/4 trans, primitives
      Double precision t3(MPQ)                 ! 3rd 1/4 trans
      Double precision t4(NPD)                 ! 4th 1/4 trans
      Integer Nxyz(3)

c--> Local variables

      Integer Lq,ninti
      Integer ma,mb,mc,md, mp,mq, nq,nr
      Integer Ia,Ja,Ka, Ib,Jb,Kb, Ic,Jc,Kc, Id,Jd,Kd
      Integer Ip,Jp,Kp, Iq,Jq,Kq, Ir,Jr,Kr
      Integer ipa,ipb,ipc,ipd, ica,icb,icc,icd
      Integer ie,icab,NCAB
      Integer ind_ca,ind_cb
      Double precision ca,cb

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
c General case:  [ab|cd]

      if (canAB .or. canCD .or. canPQ) 
     &    call errquit ('hferi_gen:can''t do canonical integrals',911,
     &       INT_ERR)

c Initialize the block of ERIs.

      NCAB = NCA*NCB
      ninti = La2*Lb2*Lc2*Ld2*NCA*NCB*NCC*NCD
      call dfill (ninti,0.0d0,eri,1)

c Loop over the components of the "A" and "B" shells.

      do ma = 1,La2

        call getNxyz(La,ma,Nxyz)
        Ia = Nxyz(1)
        Ja = Nxyz(2)
        Ka = Nxyz(3)

        do mb = 1,Lb2

          call getNxyz(Lb,mb,Nxyz)
          Ib = Nxyz(1)
          Jb = Nxyz(2)
          Kb = Nxyz(3)

c Sum across (Ip,Jp,Kp) for each value of (Iq,Jq,Kq).

          call dfill(MPP*MPQ*Lqmax3,0.0d00,t2,1)

          do Ip = 0,Ia+Ib
            do Jp = 0,Ja+Jb
              do Kp = 0,Ka+Kb

c Define the product of the Hermite expansion coefficients for 
c overlap distribution "P".

                do mp = 1,MPP
                  E3(mp) = Ep(1,mp,0,Ip,Ia,Ib)
     &                    *Ep(2,mp,0,Jp,Ja,Jb)
     &                    *Ep(3,mp,0,Kp,Ka,Kb)
                end do

                do Iq = 0,Lqmax
                  do Jq = 0,Lqmax-Iq
                    do Kq = 0,Lqmax-Iq-Jq

                      nq = IJK(Iq,Jq,Kq)
                      Ir = Ip + Iq
                      Jr = Jp + Jq
                      Kr = Kp + Kq
                      nr = IJK(Ir,Jr,Kr)

c Include the factor of (-1)**(Iq+Jq+Kq).
c Sum over Hermite functions

                      Lq = Iq + Jq + Kq
                      if( mod(Lq,2).eq.0 )then
                        do mp = 1,MPP
                          do mq = 1,MPQ
                            t2(mq,nq,mp) = t2(mq,nq,mp) 
     &                          + E3(mp)*R0(mq,mp,nr)
                          end do
                        end do
                      else
                        do mp = 1,MPP
                          do mq = 1,MPQ
                            t2(mq,nq,mp) = t2(mq,nq,mp) 
     &                          - E3(mp)*R0(mq,mp,nr)
                          end do
                        end do
                      end if

                    end do
                  end do
                end do

              end do
            end do
          end do

c Contract over shell a

          call dfill(MPQ*Lqmax3*NPB*NCA,0.0d0,t1,1)
          do ica = 1,NCA
            do mp = 1,MPP
              ipa = PairP(1,mp)
              ipb = PairP(2,mp)
              ca = Acoefs(ipa,ica)
              do nq = 1,Lqmax3
                do mq = 1,MPQ
                  t1(mq,nq,ipb,ica) = t1(mq,nq,ipb,ica) + 
     &                t2(mq,nq,mp)*ca
                end do
              end do
            end do
          end do

c Contract over shell b

          call dfill(MPQ*Lqmax3*NCAB,0.0d0,t2,1)
          icab = 0
          do ica = 1,NCA
            do icb = 1,NCB
              icab = icab+1
              do ipb = 1,NPB
                cb = Bcoefs(ipb,icb)
                do nq = 1,Lqmax3
                  do mq = 1,MPQ
                    t2(mq,nq,icab) = t2(mq,nq,icab)
     &                  +t1(mq,nq,ipb,ica)*cb
                  end do
                end do
              end do
            end do
          end do

            
c Loop over the components of the "C" and "D" shells.

          do mc = 1,Lc2

            call getNxyz(Lc,mc,Nxyz)
            Ic = Nxyz(1)
            Jc = Nxyz(2)
            Kc = Nxyz(3)

            do md = 1,Ld2

              call getNxyz(Ld,md,Nxyz)
              Id = Nxyz(1)
              Jd = Nxyz(2)
              Kd = Nxyz(3)

c Define the product of the Hermite expansion coefficients for 
c overlap distribution "Q".

              ie = 0
              do Iq = 0,Ic+Id
                do Jq = 0,Jc+Jd
                  do Kq = 0,Kc+Kd
                    do mq = 1,MPQ
                      ie = ie+1
                      E3(ie) = Eq(1,mq,0,Iq,Ic,Id)
     &                        *Eq(2,mq,0,Jq,Jc,Jd)
     &                        *Eq(3,mq,0,Kq,Kc,Kd)
                    end do
                  end do
                end do
              end do

c Sum across (Iq,Jq,Kq).

              icab = 0
              ind_ca = ma
              do ica = 1,NCA
                ind_cb = mb
                do icb = 1,NCB
                  call dfill (MPQ,0.0d0,t3,1)
                  icab = icab+1
                  ie = 0
                  do Iq = 0,Ic+Id
                    do Jq = 0,Jc+Jd
                      do Kq = 0,Kc+Kd
                        
                        nq = IJK(Iq,Jq,Kq)

c Sum over Hermite functions to give integrals

                        do mq = 1,MPQ
                          ie = ie+1
                          t3(mq) = t3(mq) + E3(ie)*t2(mq,nq,icab)
                        end do

                      end do
                    end do
                  end do

c Contract over C shell

                  do icc = 1,NCC
                    call dfill (NPD,0.0d0,t4,1)
                    do mq = 1,MPQ
                      ipc = PairQ(1,mq)
                      ipd = PairQ(2,mq)
                      t4(ipd) = t4(ipd) + t3(mq)*Ccoefs(ipc,icc)
                    end do

c Contract over D shell

                    do icd = 1,NCD
                      do ipd = 1,NPD
                        ERI(md,icd,mc,icc,ind_cb,ind_ca) =
     &                      ERI(md,icd,mc,icc,ind_cb,ind_ca) +
     &                      t4(ipd)*Dcoefs(ipd,icd)
                      end do
                    end do
                  end do

                  ind_cb = ind_cb+Lb2
                end do ! icb
                ind_ca = ind_ca+La2
              end do ! ica

            end do
          end do

        end do
      end do

      end
