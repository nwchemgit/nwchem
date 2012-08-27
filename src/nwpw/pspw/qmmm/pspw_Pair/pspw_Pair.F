c
c $Id: pspw_Pair.F 22780 2012-08-26 19:56:12Z bylaska $
c

*     *************************
*     *                       *
*     *     pspw_Pair_init    *
*     *                       *
*     *************************
      subroutine pspw_Pair_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "pspw_Pair.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j,l,k,nn,npair
      integer nkatm,nkatm_qm
      character*50 rtdbname

*     **** external functions ****
      character*4 ion_atom
      character*7 c_index_name
      integer     ion_nkatm,ion_nkatm_qm
      integer     ion_katm
      external    ion_atom
      external    c_index_name
      external    ion_nkatm,ion_nkatm_qm
      external    ion_katm
      
      nkatm    = ion_nkatm()
      nkatm_qm = ion_nkatm_qm()

*     **** allocate Pair parameters ****
      value = rtdb_parallel(.true.)

      include_qm = .false.
      npairs = 0
      do j=1,nkatm
         do i=1,nkatm
           rtdbname='pspw_Pair_ion_ion_npair:'//ion_atom(i)//ion_atom(j)
            if (rtdb_get(rtdb,rtdbname,mt_int,1,npair)) then
                npairs = npairs + npair
                if ((i.le.nkatm_qm).and.(j.le.nkatm_qm)) 
     >             include_qm = .true.
            end if
         end do
      end do

      if (npairs.gt.0) then
         value = value.and.
     >        MA_alloc_get(mt_int,nkatm*nkatm,
     >                 'npair_all',npair_all(2),npair_all(1))
         value = value.and.
     >        MA_alloc_get(mt_int,nkatm*nkatm,
     >                 'pair_start',pair_start(2),pair_start(1))

         call icopy(nkatm*nkatm,0,0,int_mb(npair_all(1)),1)
         call icopy(nkatm*nkatm,0,0,int_mb(pair_start(1)),1)

         value = value.and.
     >        MA_alloc_get(mt_int,npairs,
     >                 'type_all',type_all(2),type_all(1))
         value = value.and.
     >        MA_alloc_get(mt_dbl,3*npairs,
     >                 'param_all',param_all(2),param_all(1))
         if (.not. value) 
     >   call errquit('pspw_Pair_init:out of heap memory',0, MA_ERR)

      

*        **** read in Pair potential parameters ****
         nn = 0
         do j=1,nkatm
            do i=1,nkatm
               rtdbname='pspw_Pair_ion_ion_npair:'
     >                  //ion_atom(i)//ion_atom(j)
               if (.not.rtdb_get(rtdb,rtdbname,mt_int,1,npair))
     >            npair = 0

               int_mb(npair_all(1)+(j-1)*nkatm+(i-1))  = npair
               int_mb(pair_start(1)+(j-1)*nkatm+(i-1)) = nn
               do k=1,npair
                  rtdbname = 'pspw_Pair_ion_ion_type:'
     >                     //c_index_name(k)
     >                     //ion_atom(i)//ion_atom(j)
                  value = value.and.
     >                    rtdb_get(rtdb,rtdbname,mt_int,1,
     >                             int_mb(type_all(1)+nn))

                  rtdbname = 'pspw_Pair_ion_ion_param:'
     >                     //c_index_name(k)
     >                     //ion_atom(i)//ion_atom(j)
                  value = value.and.
     >                    rtdb_get(rtdb,rtdbname,mt_dbl,3,
     >                             dbl_mb(param_all(1)+3*nn))
                  nn = nn + 1
               end do

            end do
         end do

         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
             write(luout,*)
             write(luout,*) 'Pair Ion-Ion Parameters (units=a.u.):'
             if (include_qm)
     >          write(luout,*) '- including QM/QM Pair interactions'
             do j=1,nkatm
             do i=j,nkatm
               npair = int_mb(npair_all(1)+(j-1)*nkatm+(i-1))
               nn    = int_mb(pair_start(1)+(j-1)*nkatm+(i-1))
               if (npair.gt.0) then
                  write(luout,'(A4,1x,A4)') ion_atom(i),ion_atom(j)
                  do k=1,npair
                     if (int_mb(type_all(1)+nn+k-1).eq.1) then
                        write(luout,'(4x,A30,E10.3,A3,E10.3,A3,E10.3)')
     >                    '- Potential = A*(rij-B)**C, A:',  
     >                    dbl_mb(param_all(1)+3*(nn+k-1)),
     >                    ' B:',  
     >                    dbl_mb(param_all(1)+3*(nn+k-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+3*(nn+k-1)+2)
                     else
                        write(luout,'(4x,A30,E10.3,A3,E10.3)') 
     >                    '- Potential = A*exp(rij/B), A:',  
     >                    dbl_mb(param_all(1)+3*(nn+k-1)),
     >                    ' B:',  
     >                    dbl_mb(param_all(1)+3*(nn+k-1)+1)
                     end if
                  end do
               endif
             end do
             end do
             write(luout,*)
         end if

      end if
      return
      end

*     *************************
*     *                       *
*     *     pspw_Pair_end     *
*     *                       *
*     *************************
      subroutine pspw_Pair_end()
      implicit none

#include "mafdecls.fh"
#include "pspw_Pair.fh"
#include "errquit.fh"

      logical value

      if (npairs.gt.0) then
         value =           MA_free_heap(npair_all(2))
         value = value.and.MA_free_heap(pair_start(2))
         value = value.and.MA_free_heap(type_all(2))
         value = value.and.MA_free_heap(param_all(2))
         if (.not.value) 
     >      call errquit('pspw_Pair_end: error MA_free_heap',
     >                   0,MA_ERR)
      end if
      return
      end



c     *************************************
c     *                                   *
c     *           pspw_Pair_E             *
c     *                                   *
c     *************************************
      real*8 function pspw_Pair_E(nion,nion_qm,katm,
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

#include "mafdecls.fh"
#include "pspw_Pair.fh"

*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,jj,nkatm
      integer w1,a,k1,kk1,n1,npair,istart,k
      integer w2,b,k2,kk2
      real*8  E

*     **** external functions ****
      integer  ion_nkatm
      real*8   pspw_VPair_E_periodic,pspw_VPair_E_periodic_self
      real*8   pspw_VPair_E_onecell,pspw_VPair_E_periodic_image
      external ion_nkatm
      external pspw_VPair_E_periodic,pspw_VPair_E_periodic_self
      external pspw_VPair_E_onecell,pspw_VPair_E_periodic_image
      
      call nwpw_timing_start(40)
      E = 0.0d0

      if (npairs.gt.0) then
      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nkatm = ion_nkatm()

      dutask = 0

c     **** QM/QM VPower energy ****
      if (include_qm) then
         do j = 1,nion_qm-1
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=j+1,nion_qm
               ii = katm(i)
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0)
     >            E = E + pspw_VPair_E_periodic(npair,
     >                           int_mb(type_all(1) +istart),
     >                           dbl_mb(param_all(1)+3*istart),
     >                           rion(1,i),rion(1,j),
     >                           nshl3d,rcell)
            end do
            end if
            dutask = mod(dutask+1,np)
         end do
      end if

c     **** QM/MM VPower energy ****
      do j = nion_qm+1,nion
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=1,nion_qm
            ii = katm(i)
            npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
            istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
            if (npair.gt.0)
     >         E = E + pspw_VPair_E_periodic(npair,
     >                        int_mb(type_all(1) +istart),
     >                        dbl_mb(param_all(1)+3*istart),
     >                        rion(1,i),rion(1,j),
     >                        nshl3d,rcell)
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

c     **** MM/MM VPower 1 cell energy ****
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
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0)
     >            E = E + pspw_VPair_E_onecell(npair,
     >                        int_mb(type_all(1) +istart),
     >                        dbl_mb(param_all(1)+3*istart),
     >                        rion(1,i),rion(1,j))
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if
      dutask = mod(dutask+1,np)
      end do

c     **** MM/MM VPower self energy ****
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
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0)
     >            E = E + pspw_VPair_E_onecell(npair,
     >                        int_mb(type_all(1) +istart),
     >                        dbl_mb(param_all(1)+3*istart),
     >                        rion(1,i),rion(1,j))
               kk2 = kk2 + 1
             end do
             kk1 = kk1 + 1
           end do
         end if
         dutask = mod(dutask+1,np)
         end if
      end do

c     **** MM/MM VPair self image energy ****
      if (nshl3d.gt.1) then

      do j = nion_qm+1,nion
       if (dutask.eq.taskid) then
         jj = katm(j)
         npair  = int_mb(npair_all(1) +(jj-1)*nkatm+jj-1)
         istart = int_mb(pair_start(1)+(jj-1)*nkatm+jj-1)
         if (npair.gt.0)
     >      E = E + pspw_VPair_E_periodic_image(npair,
     >                  int_mb(type_all(1) +istart),
     >                  dbl_mb(param_all(1)+3*istart),
     >                  nshl3d,rcell)

       end if
       dutask = mod(dutask+1,np)
      end do


c     **** MM/MM VPair image energy ****
      do j = (nion_qm+1),(nion-1)
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=j+1,nion
            ii = katm(i)
            npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
            istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
            if (npair.gt.0)
     >         E = E + pspw_VPair_E_periodic_self(npair,
     >                  int_mb(type_all(1) +istart),
     >                  dbl_mb(param_all(1)+3*istart),
     >                  rion(1,i),rion(1,j),
     >                  nshl3d,rcell)
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

      end if

      if (np.gt.1) call Parallel_SumAll(E)

      end if
      call nwpw_timing_end(40)

      pspw_Pair_E = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VPair_E_onecell        *
c     *                                           *
c     *********************************************

      real*8 function pspw_VPair_E_onecell(n,t,p,r1,r2)
      implicit none
      integer n,t(*)
      real*8 p(3,*)
      real*8 r1(3)
      real*8 r2(3)

*     **** local variables ****
      integer k
      real*8  dx,dy,dz,r
      real*8  E

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      r = dsqrt(dx**2 + dy**2 + dz**2)

      E = 0.0d0
      do k=1,n
         if (t(k).eq.1) then
             E = E + p(k,1)*( (r-p(2,k))**p(3,k) )
         else if (t(k).eq.2) then
             E = E + p(k,1)*dexp( r/p(2,k) )
         endif
      end do

      pspw_VPair_E_onecell = E
      return
      end


c     *********************************************
c     *                                           *
c     *              pspw_VPair_E_periodic        *
c     *                                           *
c     *********************************************

      real*8 function pspw_VPair_E_periodic(n,t,p,r1,r2,
     >                                    nshl3d,rcell)
      implicit none
      integer n,t(*)
      real*8  p(3,*)
      real*8  r1(3)
      real*8  r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E

      E         = 0.0d0
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
                E = E + p(k,1)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.2) then
                E = E + p(k,1)*dexp( r/p(2,k) )
            endif
         end do
      end do
      pspw_VPair_E_periodic = E
      return
      end


c     *********************************************
c     *                                           *
c     *         pspw_VPair_E_periodic_self        *
c     *                                           *
c     *********************************************
      real*8 function pspw_VPair_E_periodic_self(n,t,p,r1,r2,
     >                                         nshl3d,rcell)
      implicit none
      integer n,t(*)
      real*8 p(3,*)
      real*8 r1(3)
      real*8 r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E

      E         = 0.0d0
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
                E = E + p(k,1)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.2) then
                E = E + p(k,1)*dexp( r/p(2,k) )
            endif
         end do
      end do

      pspw_VPair_E_periodic_self = E
      return
      end

c     *********************************************
c     *                                           *
c     *           pspw_VPair_E_periodic_image    *
c     *                                           *
c     *********************************************

      real*8 function pspw_VPair_E_periodic_image(n,t,p,nshl3d,rcell)
      implicit none
      integer n,t(*)
      real*8 p(3,*)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  x,y,z,r
      real*8  E

      E = 0.0d0
      do l=2,nshl3d
         x = rcell(l,1)
         y = rcell(l,2)
         z = rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
                E = E + p(k,1)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.2) then
                E = E + p(k,1)*dexp( r/p(2,k) )
            endif
         end do
      end do
      pspw_VPair_E_periodic_image = E
      return
      end




c     *************************************
c     *                                   *
c     *           pspw_Pair_fion          *
c     *                                   *
c     *************************************
      subroutine pspw_Pair_fion(nion,nion_qm,katm,
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

#include "mafdecls.fh"
#include "pspw_Pair.fh"

*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,jj,nkatm,npair,istart
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

c     **** QM/QM Pair force ****
      if (include_qm) then
         do j = 1,nion_qm-1
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=j+1,nion_qm
               ii = katm(i)
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0) 
     >            call pspw_VPair_fion_periodic(npair,
     >                             int_mb(type_all(1)+istart),
     >                             dbl_mb(param_all(1)+3*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell)
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
            npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
            istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
            if (npair.gt.0) 
     >         call pspw_VPair_fion_periodic(npair,
     >                             int_mb(type_all(1) +istart),
     >                             dbl_mb(param_all(1)+3*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell)
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
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0) 
     >            call pspw_VPair_fion_onecell(npair,
     >                             int_mb(type_all(1) +istart),
     >                             dbl_mb(param_all(1)+3*istart),
     >                             rion(1,kk1),fion(1,kk1),
     >                             rion(1,kk2),fion(1,kk2))
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if
      dutask = mod(dutask+1,np)
      end do

c     **** MM/MM Pair self energy ****
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
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0) 
     >            call pspw_VPair_fion_onecell(npair,
     >                             int_mb(type_all(1) +istart),
     >                             dbl_mb(param_all(1)+3*istart),
     >                             rion(1,kk1),fion(1,kk1),
     >                             rion(1,kk2),fion(1,kk2))
               kk2 = kk2 + 1
             end do
             kk1 = kk1 + 1
           end do
         end if
         dutask = mod(dutask+1,np)
         end if
      end do

      if (nshl3d.gt.1) then

c     **** MM/MM Pair self image energy - no force ****
c     **** MM/MM Pair image energy ****
      do j = (nion_qm+1),(nion-1)
         if (dutask.eq.taskid) then
         jj = katm(j)
         do i=j+1,nion
            ii = katm(i)
            npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
            istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
            if (npair.gt.0) 
     >         call pspw_VPair_fion_periodic_self(npair,
     >                             int_mb(type_all(1) +istart),
     >                             dbl_mb(param_all(1)+3*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell)
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
c     *           pspw_VPair_fion_periodic        *
c     *                                           *
c     *********************************************
      subroutine pspw_VPair_fion_periodic(n,t,p,
     >                                  r1,f1,r2,f2,
     >                                  nshl3d,rcell)
      implicit none
      integer n,t(*)
      real*8 p(3,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
               dVPair =  p(k,1)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
            else if (t(k).eq.2) then
               dVPair =  p(k,1)/p(2,k) * dexp( r/p(2,k) )
            else
               dVPair = 0.0d0
            endif
            f1(1) = f1(1) - (x/r)*dVPair
            f1(2) = f1(2) - (y/r)*dVPair
            f1(3) = f1(3) - (z/r)*dVPair
            f2(1) = f2(1) + (x/r)*dVPair
            f2(2) = f2(2) + (y/r)*dVPair
            f2(3) = f2(3) + (z/r)*dVPair
         end do
      end do

      return
      end



c     *********************************************
c     *                                           *
c     *      pspw_VPair_fion_periodic_self       *
c     *                                           *
c     *********************************************
      subroutine pspw_VPair_fion_periodic_self(n,t,p,
     >                                       r1,f1,r2,f2,
     >                                       nshl3d,rcell)
      implicit none
      integer n,t(*)
      real*8  p(3,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
               dVPair = p(k,1)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
            else if (t(k).eq.2) then
               dVPair = p(k,1)/p(2,k) * dexp( r/p(2,k) )
            else
               dVPair = 0.0d0
            endif
            f1(1) = f1(1) - (x/r)*dVPair
            f1(2) = f1(2) - (y/r)*dVPair
            f1(3) = f1(3) - (z/r)*dVPair
            f2(1) = f2(1) + (x/r)*dVPair
            f2(2) = f2(2) + (y/r)*dVPair
            f2(3) = f2(3) + (z/r)*dVPair
         end do
      end do

      return
      end


c     *********************************************
c     *                                           *
c     *          pspw_VPair_fion_onecell          *
c     *                                           *
c     *********************************************
      subroutine pspw_VPair_fion_onecell(n,t,p,r1,f1,r2,f2)
      implicit none
      integer n,t(*)
      real*8 p(3,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)

*     **** local variables ****
      integer k
      real*8  x,y,z,r
      real*8  dVPair

      x = r1(1) - r2(1)
      y = r1(2) - r2(2)
      z = r1(3) - r2(3)
      r = dsqrt(x**2 + y**2 + z**2)
      do k=1,n
         if (t(k).eq.1) then
            dVPair =  p(k,1)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
         else if (t(k).eq.2) then
            dVPair =  p(k,1)/p(2,k) * dexp( r/p(2,k) )
         else
            dVpair = 0.0d0
         endif
         f1(1) = f1(1) - (x/r)*dVPair
         f1(2) = f1(2) - (y/r)*dVPair
         f1(3) = f1(3) - (z/r)*dVPair
         f2(1) = f2(1) + (x/r)*dVPair
         f2(2) = f2(2) + (y/r)*dVPair
         f2(3) = f2(3) + (z/r)*dVPair
      end do

      return
      end

