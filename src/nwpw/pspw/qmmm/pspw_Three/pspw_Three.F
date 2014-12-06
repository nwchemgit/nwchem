c
c $Id$
c

*     *************************
*     *                       *
*     *     pspw_Three_init   *
*     *                       *
*     *************************
      subroutine pspw_Three_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "pspw_Three.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j,l,k,nn,ntriple,shift
      integer nkatm,nkatm_qm,nkatm2
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
      nkatm2 = nkatm*nkatm

      value = btdb_parallel(.true.)

*     **** determine the ntriples ****
      include_qm = .false.
      ntriples = 0
      do k=1,nkatm
         do j=1,nkatm
            do i=1,nkatm

*          **** check for standard triple potentials ****
           rtdbname='pspw_Three_ion_ion_ion_ntriple:'
     >              //ion_atom(i)//ion_atom(j)//ion_atom(k)
            if (btdb_get(rtdb,rtdbname,mt_int,1,ntriple)) then
                ntriples = ntriples + ntriple
                if ((i.le.nkatm_qm).and.
     >              (j.le.nkatm_qm).and.
     >              (k.le.nkatm_qm)) 
     >             include_qm = .true.
            end if
         end do
      end do

      if (ntriples.gt.0) then

*        **** allocate Three parameters ****
         value = BA_alloc_get(mt_int,nkatm*nkatm*nkatm,
     >                 'ntriple_all',ntriple_all(2),ntriple_all(1))
         value = value.and.
     >        BA_alloc_get(mt_int,nkatm*nkatm*nkatm,
     >                 'triple_start',triple_start(2),triple_start(1))
         if (.not. value) 
     >   call errquit('pspw_Three_init:out of heap memory',1,MA_ERR)
         call icopy(nkatm*nkatm*nkatm,0,0,int_mb(ntriple_all(1)),1)
         call icopy(nkatm*nkatm*nkatm,0,0,int_mb(triple_start(1)),1)

         value = value.and.
     >        BA_alloc_get(mt_int,ntriples,
     >                 'type_all',type_all(2),type_all(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,4*ntriples,
     >                 'param_all',param_all(2),param_all(1))
         if (.not. value) 
     >   call errquit('pspw_Three_init:out of heap memory',2,MA_ERR)

*        **** Generate Three potential parameters ****
         nn = 0
         do k=1,nkatm
          do j=1,nkatm
            do i=1,nkatm

*              **** read and add standard pair potentials ****
               rtdbname='pspw_Pair_ion_ion_ion_ntriple:'
     >                  //ion_atom(i)//ion_atom(j)//ion_atom(k)
               if (.not.btdb_get(rtdb,rtdbname,mt_int,1,ntriple))
     >            ntriple = 0

               shift = (k-1)*nkatm2+(j-1)*nkatm+(i-1)
               int_mb(ntriple_all(1)+shift)  = npair
               int_mb(triple_start(1)+shift) = nn
               do l=1,ntriple
                  rtdbname = 'pspw_Three_ion_ion_ion_type:'
     >                     //c_index_name(l)
     >                     //ion_atom(i)//ion_atom(j)//ion_atom(k)
                  value = value.and.
     >                    btdb_get(rtdb,rtdbname,mt_int,1,
     >                             int_mb(type_all(1)+nn))

                  rtdbname = 'pspw_Three_ion_ion_ion_param:'
     >                     //c_index_name(l)
     >                     //ion_atom(i)//ion_atom(j)//ion_atom(k)
                  value = value.and.
     >                    btdb_get(rtdb,rtdbname,mt_dbl,4,
     >                             dbl_mb(param_all(1)+4*nn))
                  nn = nn + 1
               end do

            end do
          end do
         end do

*        **** write out Three potential data ****
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER) then
             write(luout,*)
     >       'Three-Body Ion-Ion-Ion Parameters (units=a.u.):'
             if (include_qm)
     >          write(luout,*) 
     >          '- including QM/QM Three-Body interactions'
             do i=1,nkatm
             do j=i,nkatm
             do k=j,nkatm
               shift = (k-1)*nkatm2+(j-1)*nkatm+(i-1))
               ntriple = int_mb(ntriple_all(1)+shift)
               nn      = int_mb(triple_start(1)+shift)
               if (ntriple.gt.0) then
                  write(luout,'(A4,1x,A4,1x,A4)') 
     >            ion_atom(i),ion_atom(j),ion_atom(k)

                  do l=1,ntriple
                     if (int_mb(type_all(1)+nn+l-1).eq.1) then

                write(luout,'(4x,A49,E14.6,A5,E14.6,A3,E14.6,A3,E14.6)')
     >              '- Potential=A*exp(-rij/rho)-C/rij**6-D/rij**8, A:',  
     >                    dbl_mb(param_all(1)+4*(nn+l-1)),
     >                    ' rho:',  
     >                    dbl_mb(param_all(1)+4*(nn+l-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+4*(nn+l-1)+2),
     >                    ' D:',  
     >                    dbl_mb(param_all(1)+4*(nn+l-1)+3)

                     end if
                  end do
               endif
             end do
             end do
             end do
             write(luout,*)
         end if

      end if !*** ntriples.gt.0 ****

      return
      end

*     *************************
*     *                       *
*     *     pspw_Three_end    *
*     *                       *
*     *************************
      subroutine pspw_Three_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_Three.fh"
#include "errquit.fh"

      logical value

      if (ntriples.gt.0) then
         value =           BA_free_heap(ntriple_all(2))
         value = value.and.BA_free_heap(triple_start(2))
         value = value.and.BA_free_heap(type_all(2))
         value = value.and.BA_free_heap(param_all(2))
         if (.not.value) 
     >      call errquit('pspw_Three_end: error BA_free_heap',
     >                   0,MA_ERR)
      end if
      return
      end

c     *************************************
c     *                                   *
c     *       pspw_gen_Neighbor_List      *
c     *                                   *
c     *************************************
      subroutine pspw_gen_Neighbor_List(nion,rion,
     >               nshl3d,rcell,
     >               r2max,nlist_max,nlist)
      implicit none
      integer nion
      real*8 rion(3,*)
      integer nshl3d
      real*8 rcell(3,*)
      real*8 rmax
      integer nlist_max,nlist(nlist_max,*)

      ijkmax(1) = 0
      ijkmax(2) = 0
      ijkmax(3) = 0
      ijkmin(1) = 99999
      ijkmin(2) = 99999
      ijkmin(2) = 99999
      imax = 0
      do k=1,nion
         ijk(1,k) = nint(rion(1,k)/rmax)
         ijk(2,k) = nint(rion(2,k)/rmax)
         ijk(3,k) = nint(rion(3,k)/rmax)

         if (ijk(1,k).gt.ijkmax(1)) ijkmax(1) = ijk(1,k)
         if (ijk(2,k).gt.ijkmax(2)) ijkmax(2) = ijk(2,k)
         if (ijk(3,k).gt.ijkmax(3)) ijkmax(3) = ijk(3,k)

         if (ijk(1,k).lt.ijkmin(1)) ijkmin(1) = ijk(1,k)
         if (ijk(2,k).lt.ijkmin(2)) ijkmin(2) = ijk(2,k)
         if (ijk(3,k).lt.ijkmin(3)) ijkmin(3) = ijk(3,k)
      end do
      nijk(1) = ijkmax(1)-ijkmin(1)+1
      nijk(2) = ijkmax(2)-ijkmin(2)+1
      nijk(3) = ijkmax(3)-ijkmin(3)+1

      do k=ijkmin(3),ijkmax(3)
         do j=ijkmin(2),ijkmax(2)
            do i=ijkmin(1),ijkmax(1)
            end do
         end do
      end do

      return
      end



c     *************************************
c     *                                   *
c     *           pspw_Three_E            *
c     *                                   *
c     *************************************
      real*8 function pspw_Three_E(nion,nion_qm,katm,
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
#include "pspw_Three.fh"

*     **** local variables ****
      integer dutask,taskid,np
      integer i,j,ii,jj,nkatm,nkatm2
      integer w1,a,k1,kk1,n1,npair,istart,k
      integer w2,b,k2,kk2
      real*8  E

*     **** external functions ****
      integer  ion_nkatm
      real*8   pspw_VThree_E_periodic,pspw_VThree_E_periodic_self
      real*8   pspw_VThree_E_onecell,pspw_VThree_E_periodic_image
      external ion_nkatm
      external pspw_VThree_E_periodic,pspw_VThree_E_periodic_self
      external pspw_VThree_E_onecell,pspw_VThree_E_periodic_image
      
      call nwpw_timing_start(40)
      E = 0.0d0

      if (ntriple.gt.0) then

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      nkatm = ion_nkatm()
      nkatm2 = nkatm*nkatm

      dutask = 0

c     **** create neighbor list ****
      do k=1,nion
         nlist = 0
         do j=1,nion
            r2 = (rion(1,j)-rion(1,k))**2
               + (rion(2,j)-rion(2,k))**2
               + (rion(3,j)-rion(3,k))**2
            if (r2.le.r2_neighbor) then
               int_mb(list_neighbor(1)+(k-1)*nion+nlist) = j
               nlist = nlist + 1
            end if
         end do
         int_mb(nlist_neighbor(1)+k-1) = nlist
      end do


c     **** QM/QM VThree energy ****
      if (include_qm) then
         do k = 1,nion_qm-1
           if (dutask.eq.taskid) then
              kk = katm(k)

              nlist = int_mb(nlist_neighbor(1)+k-1)
              do r=1,nlist
                 j=int_mb(list_neighbor(1)+(k-1)*nion+r)
                 jj = katm(j)
                 do s=1,nlist
                    i=int_mb(list(1)+(k-1)*nion+s)
                    ii = katm(i)

                     shift = (kk-1)*nkatm2 +(jj-1)*nkatm+ii-1
                     ntriple  = int_mb(ntriple_all(1)+shift)
                     istart = int_mb(pair_start(1)+shift)
                     if (ntriple.gt.0)
     >                  E = E + pspw_VThree_E_periodic(ntriple,
     >                           int_mb(type_all(1) +istart),
     >                           dbl_mb(param_all(1)+4*istart),
     >                           rion(1,i),rion(1,j),rion(1,k)
     >                           nshl3d,rcell)
                 end do
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
     >                        dbl_mb(param_all(1)+4*istart),
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
     >                        dbl_mb(param_all(1)+4*istart),
     >                        rion(1,kk1),rion(1,kk2))
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
     >                        dbl_mb(param_all(1)+4*istart),
     >                        rion(1,kk1),rion(1,kk2))
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
     >                  dbl_mb(param_all(1)+4*istart),
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
     >                  dbl_mb(param_all(1)+4*istart),
     >                  rion(1,i),rion(1,j),
     >                  nshl3d,rcell)
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

      end if

      if (np.gt.1) call Parallel_SumAll(E)

      end if !*** npairs.gt.0 ***

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
      real*8 p(4,*)
      real*8 r1(3)
      real*8 r2(3)

*     **** local variables ****
      integer k
      real*8  dx,dy,dz,r
      real*8  E,u,u6,u12

      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      r = dsqrt(dx**2 + dy**2 + dz**2)

      E = 0.0d0
      do k=1,n
         if (t(k).eq.1) then
            u   = (p(2,k)/r)
            u6  = u**6
            u12 = u6**2
            E = E + 4.0d0*p(1,k)*(u12-u6)
         else if (t(k).eq.2) then
            E = E + p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
         else if (t(k).eq.3) then
            E = E + p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6 - p(4,k)/r**8
         else if (t(k).eq.4) then
            E = E + p(1,k)*( (r-p(2,k))**p(3,k) )
         else if (t(k).eq.5) then
            E = E + p(1,k)*dexp(-r/p(2,k))
         else if (t(k).eq.6) then
            E = E + p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))**2
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
      real*8  p(4,*)
      real*8  r1(3)
      real*8  r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E,u,u6,u12

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
               u   = (p(2,k)/r)
               u6  = u**6
               u12 = u6**2
               E = E + 4.0d0*p(1,k)*(u12-u6)
            else if (t(k).eq.2) then
               E = E + p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
            else if (t(k).eq.3) then
               E = E+ p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6 - p(4,k)/r**8
            else if (t(k).eq.4) then
                E = E + p(1,k)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.5) then
                E = E + p(1,k)*dexp(-r/p(2,k))
            else if (t(k).eq.6) then
                E = E + p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))**2
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
      real*8 p(4,*)
      real*8 r1(3)
      real*8 r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E,u,u6,u12

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
               u   = (p(2,k)/r)
               u6  = u**6
               u12 = u6**2
               E = E + 4.0d0*p(1,k)*(u12-u6)
            else if (t(k).eq.2) then
               E = E + p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
            else if (t(k).eq.3) then
               E = E+ p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6 - p(4,k)/r**8
            else if (t(k).eq.4) then
               E = E + p(1,k)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.5) then
               E = E + p(1,k)*dexp(-r/p(2,k))
            else if (t(k).eq.6) then
                E = E + p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))**2
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
      real*8 p(4,*)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  x,y,z,r
      real*8  E,u,u6,u12

      E = 0.0d0
      do l=2,nshl3d
         x = rcell(l,1)
         y = rcell(l,2)
         z = rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         do k=1,n
            if (t(k).eq.1) then
               u   = (p(2,k)/r)
               u6  = u**6
               u12 = u6**2
               !E = E + 4.0d0*p(1,k)*(u12-u6)
               E = E + 2.0d0*p(1,k)*(u12-u6)
            else if (t(k).eq.2) then
               !E = E + p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
               E = E + 0.5d0*p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
            else if (t(k).eq.3) then
               E=E+0.5d0*p(1,k)*dexp(-r/p(2,k))-p(3,k)/r**6-p(4,k)/r**8
            else if (t(k).eq.4) then
               !E = E +p(1,k)*( (r-p(2,k))**p(3,k) )
               E = E + 0.5d0*p(1,k)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.5) then
               !E = E + p(1,k)*dexp( r/p(2,k) )
               E = E + 0.5d0*p(1,k)*dexp(-r/p(2,k))
            else if (t(k).eq.6) then
                E = E+0.50d0*p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))**2
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

#include "bafdecls.fh"
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

      if (npairs.gt.0) then

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
     >                             dbl_mb(param_all(1)+4*istart),
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
     >                             dbl_mb(param_all(1)+4*istart),
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
     >                             dbl_mb(param_all(1)+4*istart),
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
     >                             dbl_mb(param_all(1)+4*istart),
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
     >                             dbl_mb(param_all(1)+4*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell)
         end do
         end if
         dutask = mod(dutask+1,np)
      end do

      end if !*** nshl3d.gt.1 ***

      end if !*** npairs.gt.0 ***


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
      real*8 p(4,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

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
               u = (p(2,k)/r)
               u6  = u**6
               u12 = u6**2
               dVPair = -(4.0d0*p(1,k)/r)*(12.0d0*u12-6.0d0*u6)
            else if (t(k).eq.2) then
               dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
            else if (t(k).eq.3) then
               dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
     >                                               +8.0d0*p(4,k)/r**9
            else if (t(k).eq.4) then
               dVPair =  p(1,k)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
            else if (t(k).eq.5) then
               dVPair = -p(1,k)/p(2,k)*dexp(-r/p(2,k))
            else if (t(k).eq.6) then
               dVpair = 2.0d0*p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))
     >                       *p(2,k)*dexp(-p(2,k)*(r-p(3,k)))
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
      real*8  p(4,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

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
               u = (p(2,k)/r)
               u6  = u**6
               u12 = u6**2
               dVPair = -(4.0d0*p(1,k)/r)*(12.0d0*u12-6.0d0*u6)
            else if (t(k).eq.2) then
               dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
            else if (t(k).eq.3) then
               dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
     >                                               +8.0d0*p(4,k)/r**9
            else if (t(k).eq.4) then
               dVPair = p(1,k)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
            else if (t(k).eq.5) then
               dVPair = -p(1,k)/p(2,k)*dexp(-r/p(2,k))
            else if (t(k).eq.6) then
               dVpair = 2.0d0*p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))
     >                       *p(2,k)*dexp(-p(2,k)*(r-p(3,k)))
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
      real*8 p(4,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)

*     **** local variables ****
      integer k
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

      x = r1(1) - r2(1)
      y = r1(2) - r2(2)
      z = r1(3) - r2(3)
      r = dsqrt(x**2 + y**2 + z**2)
      do k=1,n
         if (t(k).eq.1) then
            u = (p(2,k)/r)
            u6  = u**6
            u12 = u6**2
            dVPair = -(4.0d0*p(1,k)/r)*(12.0d0*u12-6.0d0*u6)
         else if (t(k).eq.2) then
            dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
         else if (t(k).eq.3) then
            dVPair= -p(1,k)/p(2,k)*dexp(-r/p(2,k))+6.0d0*p(3,k)/r**7
     >                                            +8.0d0*p(4,k)/r**9
         else if (t(k).eq.4) then
            dVPair =  p(1,k)*p(3,k)*( (r-p(2,k))**(p(3,k)-1.0d0) )
         else if (t(k).eq.5) then
            dVPair = -p(1,k)/p(2,k)*dexp(-r/p(2,k))
         else if (t(k).eq.6) then
            dVpair = 2.0d0*p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))
     >                    *p(2,k)*dexp(-p(2,k)*(r-p(3,k)))
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

