c
c $Id$
c

*     *************************
*     *                       *
*     *     pspw_LJ_init      *
*     *                       *
*     *************************
      subroutine pspw_LJ_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "pspw_LJ.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"

*     **** local variables
      logical value,oprint
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j,l
      integer nkatm,nkatm_qm
      real*8  sigma_epsilon(2),ei,si,ej,sj
      character*50 rtdbname

*     **** external functions ****
      logical     control_print
      character*4 ion_atom
      integer     ion_nkatm,ion_nkatm_qm
      integer     ion_katm
      external    control_print
      external    ion_atom
      external    ion_nkatm,ion_nkatm_qm
      external    ion_katm
      

      nkatm    = ion_nkatm()
      nkatm_qm = ion_nkatm_qm()

*     **** allocate LJ parameters ****
      value = btdb_parallel(.true.)
      value = BA_alloc_get(mt_dbl,nkatm,
     >                 'sigma_atom',sigma_atom(2),sigma_atom(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nkatm,
     >                 'epsilon_atom',epsilon_atom(2),epsilon_atom(1))

      value = value.and.
     >        BA_alloc_get(mt_dbl,nkatm*nkatm,
     >                 'sigma_all',sigma_all(2),sigma_all(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nkatm*nkatm,
     >                 'epsilon_all',epsilon_all(2),epsilon_all(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)


      
*     **** read in LJ parameters ***
      do i=1,nkatm
         rtdbname = 'pspw_LJ_param_ion:'//ion_atom(i)
         j = index(rtdbname,' ')-1
         value = btdb_get(rtdb,rtdbname(1:j),mt_dbl,2,sigma_epsilon)
         if (.not. value) call dcopy(2,0.0d0,0,sigma_epsilon,1)

        dbl_mb(sigma_atom(1)+i-1)  =sigma_epsilon(1)/(0.529177d0)
        dbl_mb(epsilon_atom(1)+i-1)=sigma_epsilon(2)/(23.06d0*27.2116d0)
      end do


*     **** generate pairwise LJ parameters ****
      do j=1,nkatm
         ej = dbl_mb(epsilon_atom(1)+j-1)
         sj = dbl_mb(sigma_atom(1)  +j-1) 
         do i=1,nkatm
            ei = dbl_mb(epsilon_atom(1)+i-1)
            si = dbl_mb(sigma_atom(1)  +i-1) 

            rtdbname='pspw_LJ_param_ion_ion:'//ion_atom(i)//ion_atom(j)
            value = btdb_get(rtdb,rtdbname,mt_dbl,2,sigma_epsilon)
            if (.not. value) then
               dbl_mb(sigma_all(1)  +(j-1)*nkatm+(i-1)) = 0.5d0*(si+sj)
               dbl_mb(epsilon_all(1)+(j-1)*nkatm+(i-1)) = dsqrt(ei*ej)
            else
               dbl_mb(sigma_all(1)+(j-1)*nkatm+(i-1)) 
     >         = sigma_epsilon(1)/0.529177d0
               dbl_mb(epsilon_all(1)  +(j-1)*nkatm+(i-1)) 
     >         = sigma_epsilon(2)/(23.06d0*27.2116d0)
            end if
         end do
      end do

*     **** check for LJ include_qm option ****
      if (.not.btdb_get(rtdb,'pspw_LJ_include_qm',mt_log,1,include_qm))
     >   include_qm = .false.
      
      
      call Parallel_taskid(taskid)
      oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
      if (oprint) then
          write(luout,*)
          write(luout,*) 'Ion LJ Parameters:'
          if (include_qm)
     >       write(luout,*) '- including QM/QM LJ interactions'
          do i=1,nkatm
            write(luout,'(A4,4x,A12,E10.3,A9,A12,E10.3,A9)') 
     >       ion_atom(i),
     >       'sigma:',  
     >       dbl_mb(sigma_atom(1)+i-1)*(0.529177d0), 
     >       'Angstrom',
     >       'epsilon:',
     >       dbl_mb(epsilon_atom(1)+i-1)*(23.06d0*27.2116d0),
     >       'kcal/mol'
          end do
          write(luout,*)
          if (include_qm) then
             do j=1,nkatm_qm
             do i=j,nkatm_qm
               if (dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1).gt.1e-11)
     >          write(luout,'(A4,1x,A4,4x,A7,E10.3,A9,A12,E10.3,A9)') 
     >          ion_atom(i),ion_atom(j),
     >          'sigma:',  
     >          dbl_mb(sigma_all(1)+(j-1)*nkatm+i-1)*(0.529177d0), 
     >          'Angstrom',
     >          'epsilon:',
     >       dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1)*(23.06d0*27.2116d0),
     >          'kcal/mol'
             end do
             end do
          end if
          do j=1,nkatm_qm
          do i=nkatm_qm+1,nkatm
            if (dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1).gt.1e-11)
     >       write(luout,'(A4,1x,A4,4x,A7,E10.3,A9,A12,E10.3,A9)') 
     >       ion_atom(i),ion_atom(j),
     >       'sigma:',  
     >       dbl_mb(sigma_all(1)+(j-1)*nkatm+i-1)*(0.529177d0), 
     >       'Angstrom',
     >       'epsilon:',
     >       dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1)*(23.06d0*27.2116d0),
     >       'kcal/mol'
          end do
          end do
          do j=nkatm_qm+1,nkatm
          do i=j,nkatm
            if (dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1).gt.1e-11)
     >       write(luout,'(A4,1x,A4,4x,A7,E10.3,A9,A12,E10.3,A9)') 
     >       ion_atom(i),ion_atom(j),
     >       'sigma:',  
     >       dbl_mb(sigma_all(1)+(j-1)*nkatm+i-1)*(0.529177d0), 
     >       'Angstrom',
     >       'epsilon:',
     >       dbl_mb(epsilon_all(1)+(j-1)*nkatm+i-1)*(23.06d0*27.2116d0),
     >       'kcal/mol'
          end do
          end do
          write(luout,*)
      end if

      return
      end

*     *************************
*     *                       *
*     *     pspw_LJ_end	      *
*     *                       *
*     *************************
      subroutine pspw_LJ_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_LJ.fh"
#include "errquit.fh"

      logical value

      value =           BA_free_heap(sigma_atom(2))
      value = value.and.BA_free_heap(epsilon_atom(2))
      value = value.and.BA_free_heap(sigma_all(2))
      value = value.and.BA_free_heap(epsilon_all(2))
      if (.not.value) call errquit('pspw_LJ_end: error BA_free_heap',
     >                             0,MA_ERR)
      return
      end


*     ***************************
*     *				*
*     *    pspw_LJ_param_Ion	*
*     *				*
*     ***************************
      subroutine pspw_LJ_param_Ion(i,s,e)
      implicit none
      integer i
      real*8  s,e

#include "bafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ka = ion_katm(i)
      s  = dbl_mb(sigma_all(1)  + ka-1)
      e  = dbl_mb(epsilon_all(1)+ ka-1)
      return
      end


c     *************************************
c     *                                   *
c     *           pspw_LJ_E               *
c     *                                   *
c     *************************************
      real*8 function pspw_LJ_E(nion,nion_qm,katm,
     >                          nfrag,indx_frag_start,size_frag,kfrag,
     >                          self_interaction,lmbda,
     >                          nshl3d,rcell,
     >                          rion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*)
      integer kfrag(*)
      logical self_interaction(*)
      real*8  lmbda
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*)

#include "bafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,jj,nkatm
      integer w1,a,k1,kk1,n1
      integer w2,b,k2,kk2
      real*8  e1,s1,e2,s2,E

*     **** external functions ****
      integer  ion_nkatm
      real*8   pspw_VLJ_E_periodic,pspw_VLJ_E_periodic_self
      real*8   pspw_VLJ_E_onecell,pspw_VLJ_E_periodic_image
      external ion_nkatm
      external pspw_VLJ_E_periodic,pspw_VLJ_E_periodic_self
      external pspw_VLJ_E_onecell,pspw_VLJ_E_periodic_image
      
      call nwpw_timing_start(40)
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nkatm = ion_nkatm()

      E = 0.0d0
      dutask = 0

c     **** QM/QM LJ energy ****
      if (include_qm) then
         do j = 1,nion_qm-1
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=j+1,nion_qm
               ii = katm(i)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
                  s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
                  E = E + pspw_VLJ_E_periodic(e1,s1,rion(1,i),rion(1,j),
     >                                        nshl3d,rcell)
               end if
            end do
            end if
            dutask = mod(dutask+1,np)
         end do
      end if

c     **** QM/MM LJ energy ****
      do j = nion_qm+1,nion
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=1,nion_qm
            ii = katm(i)
            e1 = lmbda*dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
            if (e1.gt.1.0d-11) then
               s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
               E = E + pspw_VLJ_E_periodic(e1,s1,rion(1,i),rion(1,j),
     >                                     nshl3d,rcell)
            end if
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

c     **** MM/MM LJ 1 cell energy ****
      do w1 = 1,nfrag-1
      if (dutask.eq.taskid) then
      do w2 = w1+1,nfrag
         k1 = indx_frag_start(w1)
         k2 = indx_frag_start(w2)
         kk1 = k1
         do a=1,size_frag(w1)
            kk2 = k2
            do b=1,size_frag(w2)
               ii = katm(kk1)
               jj = katm(kk2)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
                  s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
                  E = E + pspw_VLJ_E_onecell(e1,s1,
     >                                       rion(1,kk1),rion(1,kk2))
               end if
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if
      dutask = mod(dutask+1,np)
      end do

c     **** MM/MM LJ self energy ****
      do w1=1,nfrag
         if (self_interaction(kfrag(w1))) then
         if (dutask.eq.taskid) then
           k1 = indx_frag_start(w1)
           n1 = size_frag(w1)
           kk1 = k1
           do a=1,n1-1
             kk2 = kk1 + 1
             do b=a+1,n1
               ii = katm(kk1)
               jj = katm(kk2)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
                  s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
                  E = E + pspw_VLJ_E_onecell(e1,s1,
     >                                       rion(1,kk1),rion(1,kk2))
               end if
               kk2 = kk2 + 1
             end do
             kk1 = kk1 + 1
           end do
         end if
         dutask = mod(dutask+1,np)
         end if
      end do


c     **** MM/MM LJ self image energy ****
      if (nshl3d.gt.1) then

      do j = nion_qm+1,nion
       if (dutask.eq.taskid) then
         jj = katm(j)
         e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+jj-1)
         s1 = dbl_mb(sigma_all(1)  +(jj-1)*nkatm+jj-1)
         if ((e1.gt.1.0d-11)) then
            E = E + pspw_VLJ_E_periodic_image(e1,s1,nshl3d,rcell)
         end if
       end if
       dutask = mod(dutask+1,np)
      end do


c     **** MM/MM LJ image energy ****
      do j = (nion_qm+1),(nion-1)
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=j+1,nion
            ii = katm(i)
            e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
            if (e1.gt.1.0d-11) then
               s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
               E  = E + pspw_VLJ_E_periodic_self(e1,s1,
     >                                          rion(1,i),rion(1,j),
     >                                          nshl3d,rcell)
            end if
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

      end if

      if (np.gt.1) call Parallel_SumAll(E)
      call nwpw_timing_end(40)

      pspw_LJ_E = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_onecell           *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_onecell(epsilon12,sigma12,r1,r2)
      implicit none
      real*8 epsilon12,sigma12
      real*8 r1(3)
      real*8 r2(3)
      

*     **** local variables ****
      real*8  dx,dy,dz,r
      real*8  u,u6,u12,E

      E         = 0.0d0
      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      
      r = dsqrt(dx**2 + dy**2 + dz**2)
      u = (sigma12/r)
      u6  = u**6
      u12 = u6**2
      E = E + 4.0d0*epsilon12*(u12-u6)

      pspw_VLJ_E_onecell = E
      return
      end


c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_periodic          *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_periodic(epsilon12,sigma12,r1,r2,
     >                                    nshl3d,rcell)
      implicit none
      real*8  epsilon12,sigma12
      real*8  r1(3)
      real*8  r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  u,u6,u12,E

      E         = 0.0d0
      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         E = E + 4.0d0*epsilon12*(u12-u6)
      end do

      pspw_VLJ_E_periodic = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_periodic_self     *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_periodic_self(epsilon12,sigma12,r1,r2,
     >                                         nshl3d,rcell)
      implicit none
      real*8 epsilon12,sigma12
      real*8 r1(3)
      real*8 r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  u,u6,u12,E

      E         = 0.0d0
      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         E = E + 4.0d0*epsilon12*(u12-u6)
      end do

      pspw_VLJ_E_periodic_self = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_periodic_image    *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_periodic_image(e1,s1,nshl3d,rcell)
      implicit none
      real*8 e1,s1
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  x,y,z,r
      real*8  u,u6,u12,E

      E         = 0.0d0
      do l=2,nshl3d
         x = rcell(l,1)
         y = rcell(l,2)
         z = rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (s1/r)
         u6  = u**6
         u12 = u6**2
c         E = E + 4.0d0*epsilon12*(u12-u6)
         E = E + 2.0d0*e1*(u12-u6)
      end do

      pspw_VLJ_E_periodic_image = E
      return
      end


c     *************************************
c     *                                   *
c     *           pspw_LJ_fion            *
c     *                                   *
c     *************************************
      subroutine pspw_LJ_fion(nion,nion_qm,katm,
     >                        nfrag,indx_frag_start,size_frag,
     >                        kfrag,
     >                        self_interaction,lmbda,
     >                        nshl3d,rcell,
     >                        rion,fion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*) 
      integer kfrag(*)
      logical self_interaction(*)
      real*8  lmbda
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*)
      real*8  fion(3,*)

#include "bafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,jj,nkatm
      integer w1,a,k1,kk1,n1
      integer w2,b,k2,kk2
      real*8  e1,s1,e2,s2

*     **** external functions ****
      integer  ion_nkatm
      external ion_nkatm

      call nwpw_timing_start(40)
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nkatm = ion_nkatm()
      dutask = 0

c     **** QM/QM LJ force ****
      if (include_qm) then
         do j = 1,nion_qm-1
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=j+1,nion_qm
               ii = katm(i)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
               s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
               call pspw_VLJ_fion_periodic(e1,s1,
     >                                     rion(1,i),fion(1,i),
     >                                     rion(1,j),fion(1,j),
     >                                     nshl3d,rcell)
               end if
            end do
            end if
            dutask = mod(dutask+1,np)
         end do
      end if

c     **** QM/MM LJ energy ****
      do j = nion_qm+1,nion
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=1,nion_qm
            ii = katm(i)
            e1 = lmbda*dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
            if (e1.gt.1.0d-11) then
            s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
            call pspw_VLJ_fion_periodic(e1,s1,
     >                                  rion(1,i),fion(1,i),
     >                                  rion(1,j),fion(1,j),
     >                                  nshl3d,rcell)
            end if
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

c     **** MM/MM LJ 1 cell energy ****
      do w1 = 1,nfrag-1
      if (dutask.eq.taskid) then
      do w2 = w1+1,nfrag
         k1 = indx_frag_start(w1)
         k2 = indx_frag_start(w2)
         kk1 = k1
         do a=1,size_frag(w1)
            kk2 = k2
            do b=1,size_frag(w2)
               ii = katm(kk1)
               jj = katm(kk2)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
               s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
               call pspw_VLJ_fion_onecell(e1,s1,
     >                                    rion(1,kk1),fion(1,kk1),
     >                                    rion(1,kk2),fion(1,kk2))
               end if
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if
      dutask = mod(dutask+1,np)
      end do

c     **** MM/MM LJ self energy ****
      do w1=1,nfrag
         if (self_interaction(kfrag(w1))) then
         if (dutask.eq.taskid) then
           k1 = indx_frag_start(w1)
           n1 = size_frag(w1)
           kk1 = k1
           do a=1,n1-1
             kk2 = kk1 + 1
             do b=a+1,n1
               ii = katm(kk1)
               jj = katm(kk2)
               e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
               if (e1.gt.1.0d-11) then
               s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
               call pspw_VLJ_fion_onecell(e1,s1,
     >                                    rion(1,kk1),fion(1,kk1),
     >                                    rion(1,kk2),fion(1,kk2))
               end if
               kk2 = kk2 + 1
             end do
             kk1 = kk1 + 1
           end do
         end if
         dutask = mod(dutask+1,np)
         end if
      end do

      if (nshl3d.gt.1) then

c     **** MM/MM LJ self image energy - no force ****
c     **** MM/MM LJ image energy ****
      do j = (nion_qm+1),(nion-1)
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=j+1,nion
            ii = katm(i)
            e1 = dbl_mb(epsilon_all(1)+(jj-1)*nkatm+ii-1)
            if (e1.gt.1.0d-11) then
            s1 = dbl_mb(sigma_all(1)+(jj-1)*nkatm+ii-1)
            call pspw_VLJ_fion_periodic_self(e1,s1,
     >                                       rion(1,i),fion(1,i),
     >                                       rion(1,j),fion(1,j),
     >                                       nshl3d,rcell)
            end if
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

      end if


      call nwpw_timing_end(40)
      return
      end


c     *********************************************
c     *                                           *
c     *              pspw_VLJ_fion_periodic       *
c     *                                           *
c     *********************************************

      subroutine pspw_VLJ_fion_periodic(epsilon12,sigma12,
     >                                  r1,f1,r2,f2,
     >                                  nshl3d,rcell)
      implicit none
      real*8 epsilon12,sigma12
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  u,u6,u12,dVLJ

      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         dVLJ = -(4.0d0*epsilon12/r)*(12.0d0*u12-6.0d0*u6)
         f1(1) = f1(1) - (x/r)*dVLJ
         f1(2) = f1(2) - (y/r)*dVLJ
         f1(3) = f1(3) - (z/r)*dVLJ
         f2(1) = f2(1) + (x/r)*dVLJ
         f2(2) = f2(2) + (y/r)*dVLJ
         f2(3) = f2(3) + (z/r)*dVLJ
      end do

      return
      end



c     *********************************************
c     *                                           *
c     *         pspw_VLJ_fion_periodic_self       *
c     *                                           *
c     *********************************************

      subroutine pspw_VLJ_fion_periodic_self(epsilon12,sigma12,
     >                                       r1,f1,r2,f2,
     >                                       nshl3d,rcell)
      implicit none
      real*8  epsilon12,sigma12
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  u,u6,u12,dVLJ

      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         dVLJ = -(4.0d0*epsilon12/r)*(12.0d0*u12-6.0d0*u6)
         f1(1) = f1(1) - (x/r)*dVLJ
         f1(2) = f1(2) - (y/r)*dVLJ
         f1(3) = f1(3) - (z/r)*dVLJ
         f2(1) = f2(1) + (x/r)*dVLJ
         f2(2) = f2(2) + (y/r)*dVLJ
         f2(3) = f2(3) + (z/r)*dVLJ
      end do

      return
      end



c     *********************************************
c     *                                           *
c     *             pspw_VLJ_fion_onecell         *
c     *                                           *
c     *********************************************

      subroutine pspw_VLJ_fion_onecell(epsilon12,sigma12,r1,f1,r2,f2)
      implicit none
      real*8 sigma12,epsilon12
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)

*     **** local variables ****
      real*8  x,y,z,r
      real*8  u,u6,u12,dVLJ

      !sigma12   = 0.5d0*(s1+s2)
      !epsilon12 = dsqrt(e1*e2)
      x = r1(1) - r2(1)
      y = r1(2) - r2(2)
      z = r1(3) - r2(3)
      r = dsqrt(x**2 + y**2 + z**2)
      u = (sigma12/r)
      u6  = u**6
      u12 = u6**2
      dVLJ = -(4.0d0*epsilon12/r)*(12.0d0*u12-6.0d0*u6)
      f1(1) = f1(1) - (x/r)*dVLJ
      f1(2) = f1(2) - (y/r)*dVLJ
      f1(3) = f1(3) - (z/r)*dVLJ
      f2(1) = f2(1) + (x/r)*dVLJ
      f2(2) = f2(2) + (y/r)*dVLJ
      f2(3) = f2(3) + (z/r)*dVLJ

      return
      end

