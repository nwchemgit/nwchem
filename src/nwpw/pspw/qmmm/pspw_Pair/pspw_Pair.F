c
c $Id$
c

*     *************************
*     *                       *
*     *     pspw_Pair_init    *
*     *                       *
*     *************************
      subroutine pspw_Pair_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "pspw_Pair.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"

*     **** local variables
      logical value,oprint
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      logical include_lj_qm
      integer i,j,l,k,nn,npair
      integer nkatm,nkatm_qm
      integer sigma_atom(2),epsilon_atom(2)
      real*8 ei,ej,e,si,sj,s,A,B,s6,sigma_epsilon(2)
      character*50 rtdbname

*     **** external functions ****
      logical     control_print
      character*4 ion_atom
      character*7 c_index_name
      integer     ion_nkatm,ion_nkatm_qm
      integer     ion_katm
      external    control_print
      external    ion_atom
      external    c_index_name
      external    ion_nkatm,ion_nkatm_qm
      external    ion_katm
      
      nkatm    = ion_nkatm()
      nkatm_qm = ion_nkatm_qm()

      value = btdb_parallel(.true.)

*     **** allocate and read in LJ atom parameters ***
      value = BA_push_get(mt_dbl,nkatm,
     >                 'sigma_atom',sigma_atom(2),sigma_atom(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nkatm,
     >                 'epsilon_atom',epsilon_atom(2),epsilon_atom(1))
      if (.not. value) 
     >   call errquit('pspw_Pair_init:out of stack memory',0,MA_ERR)

      do i=1,nkatm
         rtdbname = 'pspw_LJ_param_ion:'//ion_atom(i)
         k = index(rtdbname,' ')-1
         if (.not.btdb_get(rtdb,rtdbname(1:k),mt_dbl,2,sigma_epsilon))
     >      call dcopy(2,0.0d0,0,sigma_epsilon,1)

        dbl_mb(sigma_atom(1)+i-1)  =sigma_epsilon(1)/(0.529177d0)
        dbl_mb(epsilon_atom(1)+i-1)=sigma_epsilon(2)/(23.06d0*27.2116d0)
      end do

*     **** check for LJ include_qm option ****
      if (.not.btdb_get(rtdb,'pspw_LJ_include_qm',
     >                  mt_log,1,include_lj_qm))
     >   include_lj_qm = .false.


*     **** determine the npairs ****
      include_qm = .false.
      npairs = 0
      do j=1,nkatm
         do i=1,nkatm

*           **** check for LJ atom potentials ***
            if ( ((i.gt.nkatm_qm).or.(j.gt.nkatm_qm)) .or.
     >           (include_lj_qm))                      then
               ei = dbl_mb(epsilon_atom(1)+i-1)
               ej = dbl_mb(epsilon_atom(1)+j-1)
               e = dsqrt(ei*ej)
               if (dabs(e).gt.1.0d-12) then
                 npairs = npairs + 1
                 if ((i.le.nkatm_qm).and.(j.le.nkatm_qm)) 
     >             include_qm = .true.
               end if
            end if

*          **** check for standard pair potentials ****
           rtdbname='pspw_Pair_ion_ion_npair:'//ion_atom(i)//ion_atom(j)
            if (btdb_get(rtdb,rtdbname,mt_int,1,npair)) then
                npairs = npairs + npair
                if ((i.le.nkatm_qm).and.(j.le.nkatm_qm)) 
     >             include_qm = .true.
            end if
         end do
      end do

      if (npairs.gt.0) then

*        **** allocate Pair parameters ****
         value = BA_alloc_get(mt_int,nkatm*nkatm,
     >                 'npair_all',npair_all(2),npair_all(1))
         value = value.and.
     >        BA_alloc_get(mt_int,nkatm*nkatm,
     >                 'pair_start',pair_start(2),pair_start(1))
         if (.not. value) 
     >   call errquit('pspw_Pair_init:out of heap memory',1, MA_ERR)
         call icopy(nkatm*nkatm,0,0,int_mb(npair_all(1)),1)
         call icopy(nkatm*nkatm,0,0,int_mb(pair_start(1)),1)

         value = value.and.
     >        BA_alloc_get(mt_int,npairs,
     >                 'type_all',type_all(2),type_all(1))
         value = value.and.
     >        BA_alloc_get(mt_dbl,6*npairs,
     >                 'param_all',param_all(2),param_all(1))
         if (.not. value) 
     >   call errquit('pspw_Pair_init:out of heap memory',2, MA_ERR)

*        **** Generate Pair potential parameters ****
         nn = 0
         do j=1,nkatm
            do i=1,nkatm

*              **** read and add standard pair potentials ****
               rtdbname='pspw_Pair_ion_ion_npair:'
     >                  //ion_atom(i)//ion_atom(j)
               if (.not.btdb_get(rtdb,rtdbname,mt_int,1,npair))
     >            npair = 0

               int_mb(npair_all(1)+(j-1)*nkatm+(i-1))  = npair
               int_mb(pair_start(1)+(j-1)*nkatm+(i-1)) = nn
               do k=1,npair
                  rtdbname = 'pspw_Pair_ion_ion_type:'
     >                     //c_index_name(k)
     >                     //ion_atom(i)//ion_atom(j)
                  value = value.and.
     >                    btdb_get(rtdb,rtdbname,mt_int,1,
     >                             int_mb(type_all(1)+nn))

                  rtdbname = 'pspw_Pair_ion_ion_param:'
     >                     //c_index_name(k)
     >                     //ion_atom(i)//ion_atom(j)
                  value = value.and.
     >                    btdb_get(rtdb,rtdbname,mt_dbl,6,
     >                             dbl_mb(param_all(1)+6*nn))
                  nn = nn + 1
               end do


*              **** add LJ potentials from LJ atom parameters ***
               if ( ((i.gt.nkatm_qm).or.(j.gt.nkatm_qm)) .or.
     >              (include_lj_qm))                      then
                  ei = dbl_mb(epsilon_atom(1)+i-1)
                  ej = dbl_mb(epsilon_atom(1)+j-1)
                  si = dbl_mb(sigma_atom(1)+i-1)
                  sj = dbl_mb(sigma_atom(1)+j-1)
                  e = dsqrt(ei*ej)
                  s = 0.5d0*(si+sj)
                  if (dabs(e).gt.1.0d-12) then

                     int_mb(npair_all(1)+(j-1)*nkatm+(i-1))  = npair + 1
                     int_mb(type_all(1)+nn) = 1
                     dbl_mb(param_all(1)+6*nn)   = e
                     dbl_mb(param_all(1)+6*nn+1) = s
                     nn = nn + 1

                  end if
               end if

            end do
         end do

*        **** write out Pair potential data ****
         call Parallel_taskid(taskid)
         oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
         if (oprint) then
             write(luout,*)
             write(luout,*) 'Ion LJ Atom Parameters:'
             if (include_lj_qm)
     >          write(luout,*) '- including QM/QM LJ atom interactions'
             do i=1,nkatm
               write(luout,'(A4,4x,A12,E10.3,A9,A12,E10.3,A9)')
     >          ion_atom(i),
     >          'sigma:',
     >          dbl_mb(sigma_atom(1)+i-1)*(0.529177d0),
     >          'Angstrom',
     >          'epsilon:',
     >          dbl_mb(epsilon_atom(1)+i-1)*(23.06d0*27.2116d0),
     >          'kcal/mol'
             end do
             write(luout,*)
             write(luout,*) 'Pair Ion-Ion Parameters (units=a.u.):'
             if (include_qm)
     >          write(luout,*) '- including QM/QM Pair interactions'
             do i=1,nkatm
             do j=i,nkatm
               npair = int_mb(npair_all(1)+(j-1)*nkatm+(i-1))
               nn    = int_mb(pair_start(1)+(j-1)*nkatm+(i-1))
               if (npair.gt.0) then
                  write(luout,'(A4,1x,A4)') ion_atom(i),ion_atom(j)
                  do k=1,npair
                     if (int_mb(type_all(1)+nn+k-1).eq.1) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6)') 
     >              '- Potential=4*e*((s/rij)**12-(s/rij)**6),      e:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' s:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1)
                     else if (int_mb(type_all(1)+nn+k-1).eq.2) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*exp(-rij/rho) - C/rij**6,        A:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' rho:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)
                     else if (int_mb(type_all(1)+nn+k-1).eq.3) then
                write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*exp(-rij/rho)-C/rij**6-D/rij**8, A:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' rho:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2),
     >                    ' D:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+3)
                     else if (int_mb(type_all(1)+nn+k-1).eq.4) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*(rij-B)**C,                      A:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' B:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)
                     else if (int_mb(type_all(1)+nn+k-1).eq.5) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6)') 
     >              '- Potential=A*exp(-rij/rho),                   A:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' rho:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1)
                     else if (int_mb(type_all(1)+nn+k-1).eq.6) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=De*(1-exp(-a*(rij-re)) )**2,      De:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' a:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' re:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)
                     else if (int_mb(type_all(1)+nn+k-1).eq.7) then
                write(luout,'(4x,A,A,A,
     >                        /48x,A5,E14.6,A5,E14.6,A5,E14.6,A5,E14.6,
     >                        /48x,A5,E14.6,A5,E14.6)')
     >                    '- Potential=A*exp(-rij/rho)',
     >                    '-fdamp(6,Cd;rij)*C/rij**6',
     >                    '-fdamp(8,Dd;rij)*D/rij**8,',
     >                    ' A:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' rho:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2),
     >                    ' D:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+3),
     >                    ' Cd:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+4),
     >                    ' Dd:',  
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+5)
                     else if (int_mb(type_all(1)+nn+k-1).eq.8) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6)')
     >              '- Potential=A*exp(-(rij/rho)**2),              A:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' rho:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1)
                     else if (int_mb(type_all(1)+nn+k-1).eq.9) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*exp(-(r/B)**C)                   A:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' B:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)

                     else if (int_mb(type_all(1)+nn+k-1).eq.10) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*exp(-(r/B)**C)/r                 A:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' B:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)

                     else if (int_mb(type_all(1)+nn+k-1).eq.11) then
                        write(luout,'(4x,A49,E14.6,A5,E14.6,A5,E14.6)')
     >              '- Potential=A*(1-tanh((r-B)/C))                A:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)),
     >                    ' B:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+1),
     >                    ' C:',
     >                    dbl_mb(param_all(1)+6*(nn+k-1)+2)

                     end if
                  end do
               endif
             end do
             end do
             write(luout,*)
         end if

      end if !*** npairs.gt.0 ****

*     **** deallocate memory for LJ atom data ****
      value =           BA_pop_stack(epsilon_atom(2))
      value = value.and.BA_pop_stack(sigma_atom(2))
      if (.not. value) 
     >   call errquit('pspw_Pair_init:popping stack',3,MA_ERR)

      return
      end

*     *************************
*     *                       *
*     *     pspw_Pair_end     *
*     *                       *
*     *************************
      subroutine pspw_Pair_end()
      implicit none

#include "bafdecls.fh"
#include "pspw_Pair.fh"
#include "errquit.fh"

      logical value

      if (npairs.gt.0) then
         value =           BA_free_heap(npair_all(2))
         value = value.and.BA_free_heap(pair_start(2))
         value = value.and.BA_free_heap(type_all(2))
         value = value.and.BA_free_heap(param_all(2))
         if (.not.value) 
     >      call errquit('pspw_Pair_end: error BA_free_heap',
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

#include "bafdecls.fh"
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
     >                           dbl_mb(param_all(1)+6*istart),
     >                           rion(1,i),rion(1,j),
     >                           nshl3d,rcell)
            end do
            end if
            dutask = mod(dutask+1,np)
         end do
      end if

c     **** QM/MM VPower energy ****
      if (nion_qm.gt.0) then
         do j = nion_qm+1,nion
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=1,nion_qm
               ii = katm(i)
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0)
     >            E = E + lmbda*pspw_VPair_E_periodic(npair,
     >                        int_mb(type_all(1) +istart),
     >                        dbl_mb(param_all(1)+6*istart),
     >                        rion(1,i),rion(1,j),
     >                        nshl3d,rcell)
            end do
            end if
            dutask = mod(dutask+1,np)
         end do
      end if

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
     >                        dbl_mb(param_all(1)+6*istart),
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
     >                        dbl_mb(param_all(1)+6*istart),
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
     >                  dbl_mb(param_all(1)+6*istart),
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
     >                  dbl_mb(param_all(1)+6*istart),
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


c     *************************************
c     *                                   *
c     *           pspw_Pair_Emix          *
c     *                                   *
c     *************************************
      real*8 function pspw_Pair_Emix(nion,nion_qm,katm,
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

c        **** QM/MM VPower energy ****
         do j = nion_qm+1,nion
            if (dutask.eq.taskid) then
            jj = katm(j)
            do i=1,nion_qm
               ii = katm(i)
               npair  = int_mb(npair_all(1) +(jj-1)*nkatm+ii-1)
               istart = int_mb(pair_start(1)+(jj-1)*nkatm+ii-1)
               if (npair.gt.0)
     >            E = E + pspw_VPair_E_periodic(npair,
     >                        int_mb(type_all(1) +istart),
     >                        dbl_mb(param_all(1)+6*istart),
     >                        rion(1,i),rion(1,j),
     >                        nshl3d,rcell)
            end do
            end if
            dutask = mod(dutask+1,np)
         end do

         if (np.gt.1) call Parallel_SumAll(E)

      end if !*** npairs.gt.0 ***

      call nwpw_timing_end(40)

      pspw_Pair_Emix = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VPair_fdamp             *
c     *                                           *
c     *********************************************
c
c    This routine calculates the dispersion damping function
c
c    f(nd,b;r) = 1 - c*exp(-b*r) * Sum(k=0,nd) (b*r)**k/k!

      real*8 function pspw_VPair_fdamp(nd,b,c,r)
      implicit none
      integer nd
      real*8 b,c,r

*     **** local variables ****
      integer k
      real*8 f,f0,num,denom

      f0    = 1.0d0
      num   = 1.0d0
      denom = 1.0d0
      do k=1,nd
          num   = num * (b*r)
          denom = denom * k
          f0 = f0 + num/denom
      end do
      f = 1.0d0 - c*dexp(-b*r)*f0

      pspw_VPair_fdamp = f
      return
      end

c     *********************************************
c     *                                           *
c     *              pspw_VPair_dfdamp            *
c     *                                           *
c     *********************************************
c
c    This routine calculates the derivative of the dispersion damping function
c
c    f(nd,b;r) = 1 - c*exp(-b*r) * Sum(k=0,nd) (b*r)**k/k!
c              = 1 - c*exp(-b*r) * (1 + Sum(k=1,nd) (b*r)**k/k!)
c    df(nd,b;r) = c*b*exp(-b*r) * Sum(k=0,nd) (b*r)**k/k!
c               -   c*exp(-b*r) * Sum(k=1,nd) (b*r)**k/k! * (k/r)
c
      real*8 function pspw_VPair_dfdamp(nd,b,c,r)
      implicit none
      integer nd
      real*8 b,c,r

*     **** local variables ****
      integer k
      real*8 df,f0,f1,num,denom

      f0    = 1.0d0
      f1    = 0.0d0
      num   = 1.0d0
      denom = 1.0d0
      do k=1,nd
          num   = num * (b*r)
          denom = denom * k
          f0 = f0 + (num/denom)
          f1 = f1 + (num/denom) * (k/r)
      end do
      df = c*dexp(-b*r)*(b*f0 - f1)

      pspw_VPair_dfdamp = df
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VPair_E_onecell         *
c     *                                           *
c     *********************************************

      real*8 function pspw_VPair_E_onecell(n,t,p,r1,r2)
      implicit none
      integer n,t(*)
      real*8 p(6,*)
      real*8 r1(3)
      real*8 r2(3)

*     **** local variables ****
      integer k
      real*8  dx,dy,dz,r
      real*8  E,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp
      external pspw_VPair_fdamp

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
         else if (t(k).eq.7) then
            E = E + p(1,k)*dexp(-r/p(2,k)) 
     >            - p(3,k)/r**6 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >            - p(4,k)/r**8 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
         else if (t(k).eq.8) then
            E = E + p(1,k)*dexp(-(r/p(2,k))**2)
         else if (t(k).eq.9) then
            E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))
         else if (t(k).eq.10) then
            E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))/r
         else if (t(k).eq.11) then
            E = E + p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k)))
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
      real*8  p(6,*)
      real*8  r1(3)
      real*8  r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp
      external pspw_VPair_fdamp

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
            else if (t(k).eq.7) then
               E = E + p(1,k)*dexp(-r/p(2,k)) 
     >            - p(3,k)/r**6 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >            - p(4,k)/r**8 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
            else if (t(k).eq.8) then
                E = E + p(1,k)*dexp(-(r/p(2,k))**2)
            else if (t(k).eq.9) then
                E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))
            else if (t(k).eq.10) then
                E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))/r
            else if (t(k).eq.11) then
               E = E + p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k)))
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
      real*8 p(6,*)
      real*8 r1(3)
      real*8 r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  E,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp
      external pspw_VPair_fdamp

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
            else if (t(k).eq.7) then
               E = E + p(1,k)*dexp(-r/p(2,k)) 
     >            - p(3,k)/r**6 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >            - p(4,k)/r**8 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
            else if (t(k).eq.8) then
               E = E + p(1,k)*dexp(-(r/p(2,k))**2)
            else if (t(k).eq.9) then
               E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))
            else if (t(k).eq.10) then
               E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))/r
            else if (t(k).eq.11) then
               E = E + p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k)))
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
      real*8 p(6,*)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  x,y,z,r
      real*8  E,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp
      external pspw_VPair_fdamp

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
               E = E + 4.0d0*p(1,k)*(u12-u6)
            else if (t(k).eq.2) then
               E = E+ p(1,k)*dexp(-r/p(2,k)) - p(3,k)/r**6
            else if (t(k).eq.3) then
               E=E+p(1,k)*dexp(-r/p(2,k))-p(3,k)/r**6-p(4,k)/r**8
            else if (t(k).eq.4) then
               E = E +p(1,k)*( (r-p(2,k))**p(3,k) )
            else if (t(k).eq.5) then
               E = E + p(1,k)*dexp(r/p(2,k))
            else if (t(k).eq.6) then
                E = E+p(1,k)*(1.0d0-dexp(-p(2,k)*(r-p(3,k))))**2
            else if (t(k).eq.7) then
               E = E + p(1,k)*dexp(-r/p(2,k)) 
     >            - p(3,k)/r**6 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >            - p(4,k)/r**8 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
            else if (t(k).eq.8) then
               E = E + p(1,k)*dexp((r/p(2,k))**2)
            else if (t(k).eq.9) then
               E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))
            else if (t(k).eq.10) then
               E = E + p(1,k)*dexp(-(r/p(2,k))**p(3,k))/r
            else if (t(k).eq.11) then
               E = E + p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k)))
            endif
         end do
      end do
      pspw_VPair_E_periodic_image = 0.5d0*E
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
     >                             dbl_mb(param_all(1)+6*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell,1.0d0)
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
     >                             dbl_mb(param_all(1)+6*istart),
     >                             rion(1,i),fion(1,i),
     >                             rion(1,j),fion(1,j),
     >                             nshl3d,rcell,lmbda)
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
     >                             dbl_mb(param_all(1)+6*istart),
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
     >                             dbl_mb(param_all(1)+6*istart),
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
     >                             dbl_mb(param_all(1)+6*istart),
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
     >                                  nshl3d,rcell,lmbda)
      implicit none
      integer n,t(*)
      real*8 p(6,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8 lmbda

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp,pspw_Vpair_dfdamp
      external pspw_VPair_fdamp,pspw_Vpair_dfdamp

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
            else if (t(k).eq.7) then
               dVpair = -p(1,k)/p(2,k)*dexp(-r/p(2,k)) 
     >         + 6.0d0*p(3,k)/r**7 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >         + 8.0d0*p(4,k)/r**9 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
     >         - p(3,k)/r**6 * pspw_VPair_dfdamp(6,p(5,k),1.0d0,r)
     >         - p(4,k)/r**8 * pspw_VPair_dfdamp(8,p(6,k),1.0d0,r)
            else if (t(k).eq.8) then
               dVPair = -(2.0d0*r*p(1,k)/p(2,k)**2)*dexp(-(r/p(2,k))**2)
            else if (t(k).eq.9) then
               dVPair = -(p(1,k)*p(3,k)*(r/p(2,k))**p(3,k))
     >                   *dexp(-(r/p(2,k))**p(3,k))/r
            else if (t(k).eq.10) then
               dVPair = -p(1,k)*dexp(-(r/p(2,k))**p(3,k))
     >                         *(p(3,k)*(r/p(2,k))**p(3,k)+1.0d0)
     >                         /(r**2)
            else if (t(k).eq.11) then
               dVPair = -p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k))**2)
     >                  /p(3,k)
            else
               dVPair = 0.0d0
            endif
            f1(1) = f1(1) - (x/r)*dVPair*lmbda
            f1(2) = f1(2) - (y/r)*dVPair*lmbda
            f1(3) = f1(3) - (z/r)*dVPair*lmbda
            f2(1) = f2(1) + (x/r)*dVPair*lmbda
            f2(2) = f2(2) + (y/r)*dVPair*lmbda
            f2(3) = f2(3) + (z/r)*dVPair*lmbda
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
      real*8  p(6,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l,k
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp,pspw_Vpair_dfdamp
      external pspw_VPair_fdamp,pspw_Vpair_dfdamp

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
            else if (t(k).eq.7) then
               dVpair = -p(1,k)/p(2,k)*dexp(-r/p(2,k)) 
     >         + 6.0d0*p(3,k)/r**7 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >         + 8.0d0*p(4,k)/r**9 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
     >         - p(3,k)/r**6 * pspw_VPair_dfdamp(6,p(5,k),1.0d0,r)
     >         - p(4,k)/r**8 * pspw_VPair_dfdamp(8,p(6,k),1.0d0,r)
            else if (t(k).eq.8) then
               dVPair = -(2.0d0*r*p(1,k)/p(2,k)**2)*dexp(-(r/p(2,k))**2)
            else if (t(k).eq.9) then
               dVPair = -(p(1,k)*p(3,k)*(r/p(2,k))**p(3,k))
     >                   *dexp(-(r/p(2,k))**p(3,k))/r
            else if (t(k).eq.10) then
               dVPair = -p(1,k)*dexp(-(r/p(2,k))**p(3,k))
     >                         *(p(3,k)*(r/p(2,k))**p(3,k)+1.0d0)
     >                         /(r**2)
            else if (t(k).eq.11) then
               dVPair = -p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k))**2)
     >                  /p(3,k)
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
      real*8 p(6,*)
      real*8 r1(3),f1(3)
      real*8 r2(3),f2(3)

*     **** local variables ****
      integer k
      real*8  x,y,z,r
      real*8  dVPair,u,u6,u12

*     **** external functions ****
      real*8   pspw_VPair_fdamp,pspw_Vpair_dfdamp
      external pspw_VPair_fdamp,pspw_Vpair_dfdamp

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
         else if (t(k).eq.7) then
            dVpair = -p(1,k)/p(2,k)*dexp(-r/p(2,k)) 
     >      + 6.0d0*p(3,k)/r**7 * pspw_VPair_fdamp(6,p(5,k),1.0d0,r)
     >      + 8.0d0*p(4,k)/r**9 * pspw_VPair_fdamp(8,p(6,k),1.0d0,r)
     >      - p(3,k)/r**6 * pspw_VPair_dfdamp(6,p(5,k),1.0d0,r)
     >      - p(4,k)/r**8 * pspw_VPair_dfdamp(8,p(6,k),1.0d0,r)
         else if (t(k).eq.8) then
            dVPair = -(2.0d0*r*p(1,k)/p(2,k)**2)*dexp(-(r/p(2,k))**2)
         else if (t(k).eq.9) then
            dVPair = -(p(1,k)*p(3,k)*(r/p(2,k))**p(3,k))
     >                *dexp(-(r/p(2,k))**p(3,k))/r
         else if (t(k).eq.10) then
            dVPair = -p(1,k)*dexp(-(r/p(2,k))**p(3,k))
     >                      *(p(3,k)*(r/p(2,k))**p(3,k)+1.0d0)
     >                      /(r**2)
         else if (t(k).eq.11) then
            dVPair = -p(1,k)*(1.0d0-dtanh((r-p(2,k))/p(3,k))**2)
     >               /p(3,k)
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

