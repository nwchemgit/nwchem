*     *************************************
*     *                                   *
*     *          kp_structure             *
*     *                                   *
*     *************************************
      subroutine kp_structure(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"


*     **** parallel variables ****
      integer  taskid,taskid_k,np,np_i,np_j,np_k
      integer  MASTER
      parameter(MASTER=0)

      integer nbrillall,k,isize,NN,i,j,n1,n2,n3,ii,i1,i2,i3
      integer ne(2),ne_excited(2),neall(2)
      integer itet(3,4),teigs(2,4),eigs(2),hml(2),vnm(2)
      integer ix,iy,iz,ipiv(4),ierr
      real*8  ktet(3,4),kstet(3,4),kstetm(3,4),dktet(3,4),rx,ry,rz
      real*8  ks(3),kv(3),kold(3),dk(3),dk2,pathlength,dist,e1,e2,EV
      real*8  vandermonde(4,4),A(4,4),f(4)
      character*255 full_filename
 

*     **** external functions ****
      logical  Dneall_mne_allocate_block,Dneall_m_free
      logical  psi_spin_orbit,control_spin_orbit,control_print
      integer  psi_ispin,psi_ne,control_num_kvectors_structure
      integer  control_ksvector_index
      integer  control_excited_ne
      real*8   lattice_unitg
      external Dneall_mne_allocate_block,Dneall_m_free
      external psi_spin_orbit,control_spin_orbit,control_print
      external psi_ispin,psi_ne,control_num_kvectors_structure
      external control_ksvector_index
      external control_excited_ne
      external lattice_unitg


      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (psi_ispin().eq.2) ne_excited(2) = control_excited_ne(2)
      neall(1) = ne(1) + ne_excited(1)
      neall(2) = ne(2) + ne_excited(2)
      NN=neall(1)-neall(2)
      EV=27.2116d0

      nbrillall = control_num_kvectors_structure()

      if (nbrillall.gt.1) then

      call Parallel_taskid(taskid)

*     **** allocate eigs ****
      isize = neall(1)+neall(2)
      if (.not.BA_push_get(mt_dbl,isize,'eigs',eigs(2),eigs(1)))
     >  call errquit(
     > 'kp_structure:out of stack memory',0,MA_ERR)
      if (.not.Dneall_mne_allocate_block(0,1,neall,hml))
     >  call errquit(
     > 'kb_structure:out of stack memory',1,MA_ERR)
      if (.not.Dneall_mne_allocate_block(0,3,neall,vnm))
     >  call errquit(
     > 'kb_structure:out of stack memory',1,MA_ERR)



*     **** produce eigenvalue band file(s) ****
      if (psi_ispin().eq.1) then
        call util_file_name('restricted_pspw_kp',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
      else
        call util_file_name('alpha_pspw_kp',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=58,file=full_filename,form='formatted')
        end if
        call util_file_name('beta_pspw_kp',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
        if (taskid.eq.MASTER) then
         open(unit=59,file=full_filename,form='formatted')
        end if
      end if

      if (taskid.eq.MASTER) then
         write(luout,*)
         write(luout,*) "pspw kp band structure calculation" 
         write(luout,*)
      end if
      call epsi_generate_kb_vnm(dbl_mb(vnm(1)))

      do k=1,nbrillall
         call control_ksvector_structure(k,ks)
         call control_kvector_structure(k,kv)

         if (taskid.eq.MASTER) then
            if (k.eq.1) then
               pathlength = 0.0d0
            else
               dist =        (kv(1)-kold(1))**2
               dist = dist + (kv(2)-kold(2))**2
               dist = dist + (kv(3)-kold(3))**2
               dist = dsqrt(dist)
               pathlength = pathlength + dist
            end if
            kold(1) = kv(1)
            kold(2) = kv(2)
            kold(3) = kv(3)
            write(luout,1508) k,pathlength,ks,kv
            write(luout,1500)
         end if

      
         call kb_vnmtoH(kv,neall,dbl_mb(vnm(1)),dbl_mb(hml(1)))
         call Dneall_mne_diagonalize(0,neall,
     >                              dbl_mb(hml(1)),dbl_mb(eigs(1)),
     >                              .true.)
        

         if (taskid.eq.MASTER) then
            do i=0,NN-1
               e1 = dbl_mb(eigs(1)+neall(1)-i-1)
               write(luout,1510) e1,e1*EV
            end do
            do i=0,neall(2)-1
               e1 = dbl_mb(eigs(1)+neall(1)-i-NN-1)
               e2 = dbl_mb(eigs(1)+neall(1)+neall(2)-i-1)
               write(luout,1510)  e1,e1*EV,e2,e2*EV
            end do
            write(58,'(1000E14.6)') pathlength,
     >          (dbl_mb(eigs(1)+i-1),i=1,neall(1))
            if (psi_ispin().eq.2) 
     >      write(59,'(1000E14.6)') pathlength,
     >             (dbl_mb(eigs(1)+i-1),i=1,neall(2))
         end if

      end do

      if (taskid.eq.MASTER) then
        close(58)
        if (psi_ispin().eq.2) close(59)
      end if

      if (.not.Dneall_m_free(vnm))
     >  call errquit(
     > 'kp_structure:popping stack memory',1,MA_ERR)
      if (.not.Dneall_m_free(hml))
     >  call errquit(
     > 'kp_structure:popping stack memory',1,MA_ERR)
      if (.not.BA_pop_stack(eigs(2))) 
     >  call errquit(
     > 'kp_structure:popping stack memory',0,MA_ERR)



      end if

 1500 FORMAT(/' kp orbital energies:')
 1508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'pathlength=',f10.6,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
c 1509 FORMAT(/' Brillouin zone tet-point: ',i3,i5
c     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
c     >       /'          =<',3f8.3,'>')
 1510 FORMAT(2(E18.7,' (',F8.3,'eV)'))


      return
      end 

      subroutine kb_vnmtoH(kv,neall,vnm,H)
      implicit none
      integer neall(2)
      real*8 kv(3),vnm(*),H(*)

*     **** local variables ****
      integer ms,i,j,vshift
      real*8  f,w,dk2

*     **** external functions ****
      integer  psi_ispin,psi_ne,Dneall_mne_size
      real*8   Dneall_mne_value,psi_eigenvalue,psi_virtual
      external psi_ispin,psi_ne,Dneall_mne_size
      external Dneall_mne_value,psi_eigenvalue,psi_virtual

      dk2   = kv(1)*kv(1) + kv(2)*kv(2) + kv(3)*kv(3)
      vshift = Dneall_mne_size(0,neall)
      do ms=1,psi_ispin()
         do j=1,neall(ms)
         do i=j,neall(ms)
            w =    -kv(1)*Dneall_mne_value(0,neall,ms,i,j,vnm)
            w = w - kv(2)*Dneall_mne_value(0,neall,ms,i,j,
     >                                         vnm(1+vshift))
            w = w - kv(3)*Dneall_mne_value(0,neall,ms,i,j,
     >                                         vnm(1+vshift+vshift))
            if (i.eq.j) then
               if (i.le.psi_ne(ms)) then
                  f = psi_eigenvalue(ms,i)
               else
                  f = psi_virtual(ms,i-psi_ne(ms))
               end if
               w = w + (f + 0.5d0*dk2)
               call Dneall_mne_set_value(w,0,neall,ms,i,j,H)
            else
               call Dneall_mne_set_value(w,0,neall,ms,i,j,H)
               call Dneall_mne_set_value(w,0,neall,ms,j,i,H)
            end if
         end do
         end do
      end do

      return
      end

c $Id$
