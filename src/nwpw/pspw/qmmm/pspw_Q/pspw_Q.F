c
c $Id$
c


*     **********************************
*     *	                               *
*     *          Q_Electrostatic	   *
*     *                                *
*     **********************************

*   This routine calculates the Electrostatic energy between two centers 
*
      real*8 function Q_Electrostatic(r1,q1,r2,q2)
      implicit none
      real*8 r1(3),q1
      real*8 r2(3),q2

*     **** local variables ****
      real*8 xo,yo,zo,d
     

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d = dsqrt(xo*xo + yo*yo + zo*zo)

      Q_Electrostatic = q1*q2/d
      return 
      end

*     **********************************
*     *	                               *
*     *      Q_Electrostatic_Force 	   *
*     *                                *
*     **********************************

*   This routine calculates the Electrostaticn force between two centers
*
      subroutine Q_Electrostatic_Force(r1,q1,f1,r2,q2,f2)
      implicit none
      real*8 r1(3),q1
      real*8 f1(3)
      real*8 r2(3),q2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,d0,d2
      real*8 der
      real*8 fx,fy,fz
     

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)
  
      d2 = (xo*xo + yo*yo + zo*zo)
      d0 = dsqrt(d2)

      der = -q1*q2/d2
      fx = -(xo/d0)*der
      fy = -(yo/d0)*der
      fz = -(zo/d0)*der  

      f2(1) = f2(1) + fx
      f2(2) = f2(2) + fy
      f2(3) = f2(3) + fz

      f1(1) = f1(1) - fx
      f1(2) = f1(2) - fy
      f1(3) = f1(3) - fz
      return 
      end




*     **********************************
*     *                                *
*     *      Q_Electrostatic_SmForce   *
*     *                                *
*     **********************************

*   This routine calculates the Electrostaticn force between two centers
*
      subroutine Q_Electrostatic_SmForce(Sm,r1,q1,f1,r2,q2,f2)
      implicit none
      real*8 Sm
      real*8 r1(3),q1
      real*8 f1(3)
      real*8 r2(3),q2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,d0,d2
      real*8 der
      real*8 fx,fy,fz


      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)

      d2 = (xo*xo + yo*yo + zo*zo)
      d0 = dsqrt(d2)

      der = -Sm*q1*q2/d2
      fx = -(xo/d0)*der
      fy = -(yo/d0)*der
      fz = -(zo/d0)*der

      f2(1) = f2(1) + fx
      f2(2) = f2(2) + fy
      f2(3) = f2(3) + fz

      f1(1) = f1(1) - fx
      f1(2) = f1(2) - fy
      f1(3) = f1(3) - fz
      return
      end



c     *************************************
c     *                                   *
c     *           pspw_Q_E                *
c     *                                   *
c     *************************************
      real*8 function pspw_Q_E(nion,nion_qm,katm,psp_zv,amass,
     >                          nfrag,indx_frag_start,size_frag,
     >                          kfrag,
     >                          switch_Rin,switch_Rout,
     >                          self_interaction,
     >                          nshl3d,rcell,
     >                          rion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      real*8  psp_zv(*),amass(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*) 
      integer kfrag(*)
      real*8  switch_Rin(*),switch_Rout(*)
      logical self_interaction(*)
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*)


*     **** local variables ****
      integer dutask,taskid,np
      integer w1,a,ks1,kk1,n1
      integer w2,b,ks2,kk2,n2
      real*8  Rin1,Rout1,rw1_cm(3)
      real*8  Rin2,Rout2,rw2_cm(3)
      real*8  Rin,Rout,E

*     **** external functions ****
      real*8   pspw_Q_E_frag_frag,Q_Electrostatic_self
      external pspw_Q_E_frag_frag,Q_Electrostatic_self

      call nwpw_timing_start(41) 
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      dutask = 0

      E = 0.0d0

      do w1=1,nfrag-1
         if (dutask.eq.taskid) then
         ks1 = indx_frag_start(w1)
         n1  = size_frag(w1)
         Rin1  = switch_Rin(kfrag(w1))
         Rout1 = switch_Rout(kfrag(w1))
         call pspw_Q_cm(n1,ks1,amass,rion,rw1_cm)

         do w2=w1+1,nfrag
            ks2 = indx_frag_start(w2)
            n2  = size_frag(w2)
            Rin2  = switch_Rin(kfrag(w2))
            Rout2 = switch_Rout(kfrag(w2))
            call pspw_Q_cm(n2,ks2,amass,rion,rw2_cm)

            Rin  = 0.5d0*(Rin1 +Rin2)
            Rout = 0.5d0*(Rout1+Rout2)
            E = E + pspw_Q_E_frag_frag(n1,ks1,rw1_cm,
     >                                 n2,ks2,rw2_cm,
     >                                 Rin,Rout,
     >                                 katm,psp_zv,rion,
     >                                 nshl3d,rcell)
         end do
         end if
         dutask = mod(dutask+1,np)
      end do
     

c     **** take out MM/MM Coulomb self energy ****
      do w1 = 1,nfrag
         if (.not.self_interaction(kfrag(w1))) then
         if (dutask.eq.taskid) then
         ks1 = indx_frag_start(w1)
         do a=1,size_frag(w1)
         do b=a+1,size_frag(w1)
            kk1 = ks1+a-1
            kk2 = ks1+b-1
            E = E - Q_Electrostatic_self(rion(1,kk1),psp_zv(katm(kk1)),
     >                                   rion(1,kk2),psp_zv(katm(kk2)))
         end do
         end do
         end if
         dutask = mod(dutask+1,np)
         end if
      end do
      if (np.gt.1) call Parallel_SumAll(E)

      call nwpw_timing_end(41) 
      pspw_Q_E = E
      return
      end

c     *********************************************
c     *                                           *
c     *              pspw_Q_E_frag_frag           *
c     *                                           *
c     *********************************************

      real*8 function pspw_Q_E_frag_frag(n1,ks1,rw1_cm,
     >                                   n2,ks2,rw2_cm,
     >                                   Rin,Rout,
     >                                   katm,psp_zv,rion,
     >                                   nshl3d,rcell)
      implicit none
      integer n1,ks1
      real*8  rw1_cm(3)
      integer n2,ks2
      real*8  rw2_cm(3)
      real*8  Rin,Rout
      integer katm(*)
      real*8  psp_zv(*)
      real*8  rion(3,*)
      integer nshl3d
      real*8  rcell(nshl3d,3)


*     **** local variables ****
      integer l,a,b,kk1,kk2
      real*8  x,y,z,rcm,Sm
      real*8  r1(3),q1
      real*8  r2(3),q2
      real*8  E

*     **** external functions ****
      real*8   pspw_Q_Switching,Q_Electrostatic
      external pspw_Q_Switching,Q_Electrostatic

      E = 0.0d0
      do l=1,nshl3d
         x = rw1_cm(1) - rw2_cm(1) + rcell(l,1)
         y = rw1_cm(2) - rw2_cm(2) + rcell(l,2)
         z = rw1_cm(3) - rw2_cm(3) + rcell(l,3)
         rcm = dsqrt(x**2 + y**2 + z**2)
         if (rcm.lt.Rout) then
            Sm = pspw_Q_Switching(Rin,Rout,rcm)-1.0d0
            kk1 = ks1
            do a = 1,n1
               r1(1) = rion(1,kk1) + rcell(l,1)
               r1(2) = rion(2,kk1) + rcell(l,2)
               r1(3) = rion(3,kk1) + rcell(l,3)
               q1    = psp_zv(katm(kk1))

               kk2 = ks2
               do b = 1,n2
                  r2(1) = rion(1,kk2) 
                  r2(2) = rion(2,kk2)
                  r2(3) = rion(3,kk2)
                  q2    = psp_zv(katm(kk2))
                  E = E + Sm*Q_Electrostatic(r1,q1,r2,q2)
                  kk2 = kk2 + 1
               end do
               kk1 = kk1 + 1
            end do
         end if
      end do

      pspw_Q_E_frag_frag = E
      return
      end


      

*     **********************************
*     *                                *
*     *      Q_Electrostatic_self      *
*     *                                *
*     **********************************

*   This routine calculates the Electrostatic energy between two centers
*
      real*8 function Q_Electrostatic_self(r1,q1,r2,q2)
      implicit none
      real*8 r1(3),q1
      real*8 r2(3),q2

*     **** local variables ****
      real*8 xo,yo,zo,d

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)

      d = dsqrt(xo*xo + yo*yo + zo*zo)

      Q_Electrostatic_self = q1*q2/d
      return
      end




*     ***************************************
*     *	                                    *
*     *      Q_Electrostatic_mForce_self    *
*     *                                     *
*     ***************************************

*   This routine calculates minus the Electrostatic force between two centers
*
      subroutine Q_Electrostatic_mForce_self(r1,q1,f1,r2,q2,f2)
      implicit none
      real*8 r1(3),q1
      real*8 f1(3)
      real*8 r2(3),q2
      real*8 f2(3)

*     **** local variables ****
      real*8 xo,yo,zo,d0,d2
      real*8 der
      real*8 fx,fy,fz
     

      xo = r2(1)-r1(1)
      yo = r2(2)-r1(2)
      zo = r2(3)-r1(3)

      d2 = (xo*xo + yo*yo + zo*zo)
      d0 = dsqrt(d2)

      der = -q1*q2/d2
      fx = -(xo/d0)*der
      fy = -(yo/d0)*der
      fz = -(zo/d0)*der  

      f2(1) = f2(1) - fx
      f2(2) = f2(2) - fy
      f2(3) = f2(3) - fz

      f1(1) = f1(1) + fx
      f1(2) = f1(2) + fy
      f1(3) = f1(3) + fz
      return 
      end




c     *************************************
c     *                                   *
c     *           pspw_Q_fion             *
c     *                                   *
c     *************************************
      subroutine pspw_Q_fion(nion,nion_qm,katm,psp_zv,amass,
     >                       nfrag,indx_frag_start,size_frag,
     >                       kfrag,
     >                       switch_Rin,switch_Rout,
     >                       self_interaction,
     >                       nshl3d,rcell,
     >                       rion,fion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      real*8  psp_zv(*),amass(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*) 
      integer kfrag(*)
      real*8  switch_Rin(*),switch_Rout(*)
      logical self_interaction(*)
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*),fion(3,*)


*     **** local variables ****
      integer dutask,taskid,np
      integer w1,a,ks1,kk1,n1
      integer w2,b,ks2,kk2,n2
      real*8  Rin1,Rout1,rw1_cm(3),m1
      real*8  Rin2,Rout2,rw2_cm(3),m2
      real*8  Rin,Rout

      call nwpw_timing_start(41) 
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      dutask = 0

!$OMP DO
      do w1=1,nfrag-1
         if (mod(w1-1,np).eq.taskid) then
         ks1 = indx_frag_start(w1)
         n1  = size_frag(w1)
         Rin1  = switch_Rin(kfrag(w1))
         Rout1 = switch_Rout(kfrag(w1))
         call pspw_Q_cmm(n1,ks1,amass,rion,rw1_cm,m1)

         do w2=w1+1,nfrag
            ks2 = indx_frag_start(w2)
            n2  = size_frag(w2)
            Rin2  = switch_Rin(kfrag(w2))
            Rout2 = switch_Rout(kfrag(w2))
            call pspw_Q_cmm(n2,ks2,amass,rion,rw2_cm,m2)

            Rin  = 0.5d0*(Rin1 +Rin2)
            Rout = 0.5d0*(Rout1+Rout2)
            call pspw_Q_fion_frag_frag(n1,ks1,rw1_cm,m1,
     >                                 n2,ks2,rw2_cm,m2,
     >                                 Rin,Rout,
     >                                 katm,psp_zv,amass,rion,fion,
     >                                 nshl3d,rcell)
         end do
         end if
         !dutask = mod(dutask+1,np)
      end do
!$OMP END DO

c     **** take out MM/MM Coulomb self energy ****
!$OMP DO
      do w1 = 1,nfrag
         if (.not.self_interaction(kfrag(w1))) then
         if (mod(w1-1,np).eq.taskid) then
         ks1 = indx_frag_start(w1)
         do a=1,size_frag(w1)
         do b=a+1,size_frag(w1)
            kk1 = ks1+a-1
            kk2 = ks1+b-1
            call Q_Electrostatic_mForce_self(
     >               rion(1,kk1),psp_zv(katm(kk1)),fion(1,kk1),
     >               rion(1,kk2),psp_zv(katm(kk2)),fion(1,kk2))
         end do
         end do
         end if
         !dutask = mod(dutask+1,np)
         end if
      end do
!$OMP END DO

      call nwpw_timing_end(41) 
      return
      end

c     *********************************************
c     *                                           *
c     *              pspw_Q_fion_frag_frag           *
c     *                                           *
c     *********************************************

      subroutine pspw_Q_fion_frag_frag(n1,ks1,rw1_cm,m1,
     >                                 n2,ks2,rw2_cm,m2,
     >                                 Rin,Rout,
     >                                 katm,psp_zv,amass,rion,fion,
     >                                 nshl3d,rcell)
      implicit none
      integer n1,ks1
      real*8  rw1_cm(3),m1
      integer n2,ks2
      real*8  rw2_cm(3),m2
      real*8  Rin,Rout
      integer katm(*)
      real*8  psp_zv(*)
      real*8  amass(*)
      real*8  rion(3,*)
      real*8  fion(3,*)
      integer nshl3d
      real*8  rcell(nshl3d,3)


*     **** local variables ****
      integer l,a,b,kk1,kk2
      real*8  x,y,z,rcm,Sm,dSm
      real*8  r1(3),q1
      real*8  r2(3),q2
      real*8  E

*     **** external functions ****
      real*8   Q_Electrostatic
      external Q_Electrostatic


      do l=1,nshl3d
         x = rw1_cm(1) - rw2_cm(1) + rcell(l,1)
         y = rw1_cm(2) - rw2_cm(2) + rcell(l,2)
         z = rw1_cm(3) - rw2_cm(3) + rcell(l,3)
         rcm = dsqrt(x**2 + y**2 + z**2)
         if (rcm.lt.Rout) then
            call pspw_Q_dSwitching(Rin,Rout,rcm,Sm,dSm)
            Sm = Sm - 1.0d0

*           *** calculate E, and -Sm*grad(E) ****
            E = 0.0d0
            kk1 = ks1
            do a = 1,n1
               r1(1) = rion(1,kk1) + rcell(l,1)
               r1(2) = rion(2,kk1) + rcell(l,2)
               r1(3) = rion(3,kk1) + rcell(l,3)
               q1    = psp_zv(katm(kk1))

               kk2 = ks2
               do b = 1,n2
                  r2(1) = rion(1,kk2) 
                  r2(2) = rion(2,kk2)
                  r2(3) = rion(3,kk2)
                  q2    = psp_zv(katm(kk2))
                  E = E + Q_Electrostatic(r1,q1,r2,q2)

                  call Q_Electrostatic_SmForce(Sm,r1,q1,fion(1,kk1),
     >                                            r2,q2,fion(1,kk2))
                  kk2 = kk2 + 1
               end do
               kk1 = kk1 + 1
            end do

*           *** calculate -E*grad(Sm) ****
            kk1 = ks1
            do a=1,n1
               fion(1,kk1) = fion(1,kk1) - E*dSm*(x/rcm)*amass(kk1)/m1
               fion(2,kk1) = fion(2,kk1) - E*dSm*(y/rcm)*amass(kk1)/m1
               fion(3,kk1) = fion(3,kk1) - E*dSm*(z/rcm)*amass(kk1)/m1
               kk1 = kk1 + 1
            end do
            kk2 = ks2
            do b=1,n2
               fion(1,kk2) = fion(1,kk2) + E*dSm*(x/rcm)*amass(kk2)/m2
               fion(2,kk2) = fion(2,kk2) + E*dSm*(y/rcm)*amass(kk2)/m2
               fion(3,kk2) = fion(3,kk2) + E*dSm*(z/rcm)*amass(kk2)/m2
               kk2 = kk2 + 1
            end do

         end if
      end do

      return
      end


*     ******************************
*     *                            *
*     *       pspw_Q_Switching     *
*     *                            *
*     ******************************
      real*8 function pspw_Q_Switching(Rin,Rout,r)
      implicit none
      real*8 Rin,Rout,r

*     **** ST2 parameters ****
c      real*8 Rl,Ru,c1,c2
c      parameter (Rl=2.0160d0/0.529177d0)
c      parameter (Ru=3.1287d0/0.529177d0)
c      parameter (c1 = (Ru-Rl)**3)
c      parameter (c2 = (3.0d0*Ru-Rl))
       real*8 c1,c2,s

      if (r.le.Rin) then
         s = 0.0d0
      else if (r.ge.Rout) then
         s = 1.0d0
      else
         c1 = (Rout-Rin)**3
         c2 = 3.0d0*Rout-Rin
         s = ((r-Rin)**2) * (c2-2.0d0*r)/c1
      end if

      pspw_Q_Switching = s
      return
      end


*     ******************************
*     *                            *
*     *       pspw_Q_dSwitching     *
*     *                            *
*     ******************************
      subroutine pspw_Q_dSwitching(Rin,Rout,r,S,dS)
      implicit none
      real*8 Rin,Rout,r,S,dS

*     **** ST2 parameters ****
c      real*8 Rl,Ru,c1,c2
c      parameter (Rl=2.0160d0/0.529177d0)
c      parameter (Ru=3.1287d0/0.529177d0)
c      parameter (c1 = (Ru-Rl)**3)
c      parameter (c2 = (3.0d0*Ru-Rl))
       real*8 c1,c2

      if (r.le.Rin) then
         S  = 0.0d0
         dS = 0.0d0
      else if (r.ge.Rout) then
         S  = 1.0d0
         dS = 0.0d0
      else
         c1 = (Rout-Rin)**3
         c2 = 3.0d0*Rout-Rin
         S  =       ((r-Rin)**2) * (c2-2.0d0*r)/c1
         
         dS = 2.0d0*(r-Rin)    * ((c2-2.0d0*r)/c1)
     >      - 2.0d0*((r-Rin)**2) /c1
      end if

      return
      end

*     ******************************
*     *                            *
*     *       pspw_Q_cm            *
*     *                            *
*     ******************************
      subroutine pspw_Q_cm(n,ks,amass,rion,rcm)
      implicit none
      integer n,ks
      real*8 amass(*),rion(3,*),rcm(3)

*     **** local variables ****
      integer k,kk
      real*8  m

      m      = 0.0d0
      rcm(1) = 0.0d0
      rcm(2) = 0.0d0
      rcm(3) = 0.0d0
      kk = ks
      do k=1,n
         rcm(1) = rcm(1) + amass(kk)*rion(1,kk)
         rcm(2) = rcm(2) + amass(kk)*rion(2,kk)
         rcm(3) = rcm(3) + amass(kk)*rion(3,kk)
         m = m + amass(kk)
         kk = kk + 1
      end do
      rcm(1) = rcm(1)/m
      rcm(2) = rcm(2)/m
      rcm(3) = rcm(3)/m

      return
      end

      
*     ******************************
*     *                            *
*     *       pspw_Q_cmm           *
*     *                            *
*     ******************************
      subroutine pspw_Q_cmm(n,ks,amass,rion,rcm,m)
      implicit none
      integer n,ks
      real*8 amass(*),rion(3,*),rcm(3)
      real*8  m

*     **** local variables ****
      integer k,kk

      m      = 0.0d0
      rcm(1) = 0.0d0
      rcm(2) = 0.0d0
      rcm(3) = 0.0d0
      kk = ks
      do k=1,n
         rcm(1) = rcm(1) + amass(kk)*rion(1,kk)
         rcm(2) = rcm(2) + amass(kk)*rion(2,kk)
         rcm(3) = rcm(3) + amass(kk)*rion(3,kk)
         m = m + amass(kk)
         kk = kk + 1
      end do
      rcm(1) = rcm(1)/m
      rcm(2) = rcm(2)/m
      rcm(3) = rcm(3)/m

      return
      end


