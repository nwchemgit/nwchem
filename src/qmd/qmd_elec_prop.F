c
c     Electronic propagation
c
      subroutine qmd_elec_prop(rtdb,nat,zden,en1,en2,coup1,coup2,
     $                         nstates,v,mass,KE,PE,curr_state,dtn,
     $                         nstep,dte,do_tdks,deco)
c
      implicit none
c
#include "qmd_common.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util_params.fh"
c
      integer rtdb
      integer nat                ! number of atoms
      integer curr_state
      integer nstates
      integer nstep
      double complex zden(nstates,nstates)
      double precision en1(nstates)
      double precision en2(nstates)
      double precision coup1(nstates,nstates)
      double precision coup2(nstates,nstates)
      double precision dtn, dte
      double precision mass(nat) ! masses
      double precision v(3,nat)  ! velocities
      double precision KE, PE
      logical do_tdks
      logical deco
c
      integer instate
      integer elec_steps
      integer l_hi, k_hi
      integer l_hp, k_hp
      integer l_zham, k_zham
      integer l_ztden1, k_ztden1
      integer l_ztden2, k_ztden2
      integer l_zrk4s, k_zrk4s
      integer i,j
      double precision new_KE
      double precision en_gap
      double precision thresh
c
      character*30 pname
c
c     preliminaries
      pname='qmd_elec_prop: '
c
c     THRESHOLD HARD CODED
      thresh=1.d-8
c
      if (nstep.eq.1) then
c       if first step, we only have only coupling, need to wait
c       electronic integration will not start until ts=2
        coup1(:,:)=coup2(:,:)
        en1(:)=en2(:)
      else
c       take care of some initial things
        call qmd_kinetic(nat,mass,v,KE)
        instate=curr_state
        elec_steps=nint(dtn/dte)
c
c       we don't do the electronic propagation on nstep.eq.1
c       since we don't have enough values of the couplings/energies
c       to extrapolate/interpolate so we double up on nstep.eq.2
        if (nstep.eq.2) elec_steps=2.d0*elec_steps
c
c       allocate local memory
        if (.not.ma_push_get(mt_dbl,nstates,'hopint',l_hi,k_hi))
     $     call errquit(pname//'ma_push_get hopint',nstates,MA_ERR)
        if (.not.ma_push_get(mt_dbl,nstates,'hopprob',l_hp,k_hp))
     $     call errquit(pname//'ma_push_get hopprob',nstates,MA_ERR)
        if (.not.ma_push_get(mt_dcpl,nstates**2,'zham',l_zham,k_zham))
     $     call errquit(pname//'ma_push_get zham',nstates**2,MA_ERR)
        if (.not.ma_push_get(mt_dcpl,nstates**2,'ztden1',
     $                        l_ztden1,k_ztden1))
     $     call errquit(pname//'ma_push_get ztden1',nstates**2,MA_ERR)
        if (.not.ma_push_get(mt_dcpl,nstates**2,'ztden2',
     $                        l_ztden2,k_ztden2))
     $     call errquit(pname//'ma_push_get ztden2',nstates**2,MA_ERR)
        if (.not.ma_push_get(mt_dcpl,nstates**2,'zrk4s',
     $                        l_zrk4s,k_zrk4s))
     $     call errquit(pname//'ma_push_get zrk4s',nstates**2,MA_ERR)
c
c       propagate
        call qmd_elec_rk4(nstates,curr_state,elec_steps,dte,zden,en1,
     $                    en2,coup1,coup2,dbl_mb(k_hi),dbl_mb(k_hp),
     $                    dcpl_mb(k_zham),dcpl_mb(k_ztden1),
     $                    dcpl_mb(k_ztden2),dcpl_mb(k_zrk4s),dtn,nstep)
c
c       Check if energy of current state is lower than the energy
c        of state below it or higher than state above it
c        (i.e. negative excitation energy)
        if (curr_state.gt.1) then
          if (en2(curr_state).lt.en2(curr_state-1))
     $      curr_state=curr_state-1
        end if
        if (curr_state.lt.nstates) then
          if (en2(curr_state).gt.en2(curr_state+1))
     $      curr_state=curr_state+1
        end if
c
c       decoherence 
        if (deco) call qmd_deco(nstates,curr_state,dtn,KE,en2,
     $                          zden)
c
        if (curr_state.ne.instate) then
c        Don't have coupling vector, so using isotropic rescale
CCCc       Hammes-Schiffer and Tully, J. Chem. Phys. 101, 4657 (1994)
CCC         alpha=0.d0
CCC         beta=0.d0
CCC         do i=1,nat
CCC          alpha=alpha+0.5d0*sum(g(:,i)*g(:,i))/(mass(i)*camu2au)
CCC          beta=beta+sum(v(:,i)*g(:,i))
CCC         end do
CCC         en_gap=state_en(instate,2)-state_en(curr_state,2)
CCC         discrim=beta**2+(4.d0*alpha*en_gap)
CCC         if (discrim.ge.0.d0) then
CCC          if (beta.lt.0.d0) then
CCC           gamm=(beta+dsqrt(discrim))/(2.d0*alpha)
CCC          else
CCC           gamm=(beta-dsqrt(discrim))/(2.d0*alpha)
CCC          end if
CCC          do i=1,nat
CCC           v(:,i)=v(:,i)-((gamm*g(:,i))/(mass(i)*camu2au))
CCC          end do
          en_gap=en2(curr_state)-en2(instate)
          if ((KE-en_gap).gt.thresh) then
             PE=en2(curr_state)
             new_KE=KE-en_gap
             call qmd_rescale(nat,mass,v,KE,new_KE)
             if (.not.do_tdks) then
               if (curr_state.eq.1) then
                 if (.not.rtdb_cput(rtdb,'task:theory',1,'dft'))
     $             call errquit(pname//'change theory fail',0,RTDB_ERR)
               else 
                 if (.not.rtdb_cput(rtdb,'task:theory',1,'tddft'))
     $             call errquit(pname//'change theory fail',0,RTDB_ERR)
                 if (.not.rtdb_put(rtdb,'tddft_grad:isinglet_roots',
     $                             mt_int,1,curr_state-1)) 
     $             call errquit(pname//'failed to write iroot',0,
     $                          RTDB_ERR)
                 if (.not.rtdb_put(rtdb,'tddft_grad:iroots',mt_int,1,
     $                             curr_state-1))
     $             call errquit(pname//'failed to write iroot',0,
     $                          RTDB_ERR)
                 if (.not.rtdb_put(rtdb,'tddft:target',mt_int,1,
     $                             curr_state-1))
     $             call errquit(pname//'failed to write targ',0,
     $                          RTDB_ERR)
               end if
             end if
          else
            curr_state=instate
            v(:,:)=-v(:,:)
            PE=en2(curr_state)
          end if
        end if
c
        coup1(:,:)=coup2(:,:)
        en1(:)=en2(:)
c
c       free memory
        if(.not.ma_chop_stack(l_hi))
     +   call errquit(pname//'ma_chop_stack hopint',0,MA_ERR)
        if(.not.ma_verify_allocator_stuff())
     +   call errquit(pname//'ma_verify_allocator_stuff',0,MA_ERR)
      end if
c
      call ga_sync()
c
      return
      end
c
c     RK4 integration
c
      subroutine qmd_elec_rk4(nstates,curr_state,elec_steps,dte,zden,
     $                        en1,en2,coup1,coup2,hopint,hopprob,zham,
     $                        ztden1,ztden2,zrk4s,dtn,nstep)
c
      implicit none
c
#include "qmd_common.fh"
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer nstates
      integer curr_state
      integer elec_steps
      integer nstep
      double precision dte
      double precision dtn
      double precision en1(nstates)
      double precision en2(nstates)
      double precision coup1(nstates,nstates)
      double precision coup2(nstates,nstates)
      double precision hopint(nstates)
      double precision hopprob(nstates)
      double complex zden(nstates,nstates)
      double complex zham(nstates,nstates)
      double complex ztden1(nstates,nstates)
      double complex ztden2(nstates,nstates)
      double complex zrk4s(nstates,nstates)
c
      integer istep
      integer irk4
      integer is, js
      integer i, j
      double precision pfac1(4)
      double precision pfac2(4)
      double precision en_gap
      double precision alpha, beta
      double precision gamm, discrim
      double precision time, zeta
      double precision hopsum
      double precision trace, thresh
      double complex zscr1
      double complex zi1, z1
      parameter (zi1=(0.d0,1.d0))
      parameter (z1=(1.d0,0.d0))
c
      external qmd_ran1
      double precision qmd_ran1
c
      character*30 pname
c
c     preliminaries
      pname='qmd_elec_rk4: '
c
c     THRESHOLD HARD CODED
      thresh=1.d-8
c
      pfac1(1)=dte*0.d0
      pfac1(2)=dte*5.d-1
      pfac1(3)=dte*5.d-1
      pfac1(4)=dte*1.d0
      pfac2(1)=dte/6.d0
      pfac2(2)=dte/3.d0
      pfac2(3)=dte/3.d0
      pfac2(4)=dte/6.d0
      hopint(:)=(-2.d0)*dble(zden(:,curr_state))*coup1(:,curr_state)
c     Propagate density matrix
c      dP/dt = -i[H,P]
c      dP/dt = -i(HP-PH)
c      dP/dt = +i(PH-HP)
      do istep=0,elec_steps-1
c       Construct Hamiltonian
c       H=V-id; V->energies, d->couplings
        zham=(0.d0,0.d0)
        do i=1,nstates
          if (nstep.eq.2) then
c           need extra integration
c           extrapolate and interpolate energies at (t-dt) and (t)
            zham(i,i)=en1(i)+((dble(istep)-5.d-1*
     $                dble(elec_steps))/(5.d-1*dble(elec_steps)))*
     $                (en2(i)-en1(i))
          else
c           interpolate energies at (t-dt) and t
            zham(i,i)=en1(i)+(dble(istep)/dble(elec_steps))*
     $                (en2(i)-en1(i))
          end if
          do j=1,nstates
            if (i.ne.j) then
              if (nstep.eq.2) then
c               need extra integration
c               extrapolate and interpolate couplings at (t-dt) and (t)
                zham(j,i)=zham(j,i)-zi1*(coup1(j,i)+((dble(istep)-5.d-1*
     $                    dble(elec_steps))/(5.d-1*dble(elec_steps)))*
     $                    (coup2(j,i)-coup1(j,i)))
              else
c               interpolate couplings at (t-dt) and t
                zham(j,i)=zham(j,i)-zi1*(coup1(j,i)+((dble(istep)/
     $                    dble(elec_steps))*(coup2(j,i)-coup1(j,i))))
              end if
            end if
          end do
        end do
        zrk4s(:,:)=(0.d0,0.d0)
        ztden2(:,:)=zden(:,:)
c       4th order Runge-Kutta
c        P(t+1)=P(t)+(dt/6)*(k1+2*k2+2*k3+k4)
c         k1=i{P1*H-H*P1}, P1=P(t)
c         k2=i{P2*H-H*P2}, P2=P(t)+(dt/2)*k1
c         k3=i{P3*H-H*P3}, P3=P(t)+(dt/2)*k2
c         k4=i{P4*H-H*P4}, P4=P(t)+dt*k3
        do irk4=1,4
         ztden1(:,:)=zden(:,:)+pfac1(irk4)*zrk4s(:,:)
         zrk4s(:,:)=(0.d0,1.d0)*(matmul(ztden1(:,:),zham(:,:))-
     $                 matmul(zham(:,:),ztden1(:,:)))
         ztden2(:,:)=ztden2(:,:)+pfac2(irk4)*zrk4s(:,:)
        end do
        zden(:,:)=ztden2(:,:)
c       Composite Simpson's rule integration
        if ((mod(istep,2).ne.0).and.(istep.ne.elec_steps-1)) then
         hopint(:)=hopint(:)-4.d0*(2.d0*dble(zden(:,curr_state)*
     $             dble(zi1*zham(:,curr_state))))
        else if ((mod(istep,2).eq.0).and.(istep.ne.elec_steps-1)) then
         hopint(:)=hopint(:)-2.d0*(2.d0*dble(zden(:,curr_state)*
     $             dble(zi1*zham(:,curr_state))))
        else
         hopint(:)=hopint(:)-2.d0*dble(zden(:,curr_state)*
     $             dble(zi1*zham(:,curr_state)))
        end if
        if (nstep.eq.2) then
          time=dble(nstep-2)*dtn+dble(istep)*dte
        else
          time=dble(nstep-1)*dtn+dble(istep)*dte
        end if
C        if (ga_nodeid().eq.0) then
C          write(luout,'((A4,f15.5),$)')'EP',time
C          do i=1,nstates
C            write(luout,'((f18.8),$)') abs(zden(i,i))
C          end do
C          write(luout,*)
C        end if
      end do
c
c     Checks
      trace=0.d0
      do i=1,nstates
        trace=trace+abs(zden(i,i))
      end do
      if (abs(trace-1.d0).gt.thresh)
     +   call errquit(pname//'Tr[P]=/=1',0,CALC_ERR)
c      ztden1(:,:)=(0.d0,0.d0)
c      ztden1=matmul(zden,zden)-zden
c      if (maxval(abs(ztden1)).gt.thresh)
c     +   call errquit(pname//'[P^2-P]=/=0',0,CALC_ERR)
c
      hopint(:)=hopint(:)*(dte/3.d0)
      hopint(curr_state)=0.d0
      hopprob(:)=max(0.d0,hopint(:)/abs(zden(curr_state,curr_state)))
      zeta=qmd_ran1(idum)
      do i=1,nstates
        hopsum=0.d0
        do j=1,i
          hopsum=hopsum+hopprob(j)
        end do
        if ((hopsum-zeta).gt.thresh) then
          curr_state=i
          exit
        end if
      end do
c
      return
      end
c
c     Decoherence
c       Granucci & Persico J. Chem. Phys. 126, 134114 (2007)
c
      subroutine qmd_deco(nstates,cs,dtn,KE,en,zden)
c
      implicit none
c
#include "qmd_common.fh"
#include "stdio.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer nstates
      integer cs
      double precision dtn
      double precision KE
      double precision en(nstates)
      double complex zden(nstates,nstates)
c
      integer j, k
      double precision tau1
      double precision tau2
      double precision trP
c
      character*30 pname
c
c     preliminaries
      pname='qmd_deco: '
c
      trP=0.d0
      do j=1,nstates
        if (j.ne.cs) then
          tau1=exp(-dtn/((1.d0/abs(en(j)-en(cs)))*(1.d0+(1.d-1/KE))))
          zden(j,j)=zden(j,j)*tau1*tau1
          trP=trP+abs(zden(j,j))
        end if
      end do
      do j=1,nstates
        do k=1,nstates
          if (j.ne.k) then
            if (j.eq.cs) then
              tau1=exp(-dtn/((1.d0/abs(en(k)-en(cs)))*
     $             (1.d0+(1.d-1/KE))))
              tau2=sqrt((1.d0-trP)/abs(zden(cs,cs)))
            else if (k.eq.cs) then
              tau1=exp(-dtn/((1.d0/abs(en(j)-en(cs)))*
     $             (1.d0+(1.d-1/KE))))
              tau2=sqrt((1.d0-trP)/abs(zden(cs,cs)))
            else
              tau1=exp(-dtn/((1.d0/abs(en(j)-en(cs)))*
     $             (1.d0+(1.d-1/KE))))
              tau2=exp(-dtn/((1.d0/abs(en(k)-en(cs)))*
     $             (1.d0+(1.d-1/KE))))
            end if
            zden(j,k)=zden(j,k)*tau1*tau2
          end if
        end do
      end do
      zden(cs,cs)=(1.d0,0.d0)-trP
c
      return
      end
