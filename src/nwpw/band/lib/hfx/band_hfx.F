*
* $Id$
*

*     *************************
*     *                       *
*     *     band_init_HFX     *
*     *                       *
*     *************************
      subroutine band_init_HFX(rtdb,nbrill0,ispin0,ne)
      implicit none
      integer rtdb
      integer nbrill0,ispin0
      integer ne(2)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "band_hfx.fh"

*     **** local variables ****
      logical value
      integer ma_type
      integer n1,n2,n3,mapping,ms,neq(2)

*     **** external functions ****
      integer  control_version,control_mapping,Pneb_nbrillq
      external control_version,control_mapping,Pneb_nbrillq
      integer  cpsi_data_alloc
      external cpsi_data_alloc

      nbrillioun = nbrill0
      nbrillq = Pneb_nbrillq()
      ispin = ispin0
      norbs(1) = 0
      norbs(2) = 0
      ehfx = 0.0d0
      hfx_on = .false.
      call C3dB_nfft3d(1,nfft3d)

      if (.not.btdb_get(rtdb,'band:HFX',mt_log,1,hfx_on))
     >   hfx_on = .false.


*     **** get the number of HFX orbitals ****
      if (hfx_on) then
         norbs(1) = ne(1)
         norbs(2) = ne(2)
         neall = ne(1)+ne(2)

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_screening_radius',
     >                      mt_dbl,1,rcut)) 
     >       rcut = 8.0d0

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_screening_power',
     >                      mt_dbl,1,pp)) 
     >       pp = 8.0d0

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_screening_type',
     >                      mt_int,1,flag)) 
     >       flag = 0

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_relax',
     >                      mt_log,1,relaxed)) 
     >       relaxed = .true.

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_solver_type',
     >                      mt_int,1,solver_type)) then

            solver_type = 1
         end if

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_parameter',
     >                       mt_dbl,1,HFX_parameter))
     >       HFX_parameter = 1.0d0

         if (.not. btdb_get(rtdb,
     >                      'band:HFX_print_orbital_contribution',
     >                       mt_log,1,orb_contribution))
     >       orb_contribution = .false.
 

*        **** initialize coulomb_screened ****
         call c_coulomb_screened_init(flag,rcut,pp)

*      
*        **** initialize orb_contribution ****
c         do ms=1,ispin
c           value = BA_alloc_get(mt_dbl,norbs(ms),
c     >                'ehfx_orb',ehfx_orb(2,ms),ehfx_orb(1,ms))
c           if (.not. value)   
c     >       call errquit('band_init_HFX: out of heap memory',1, MA_ERR)
c         end do

      end if


c     **** define extra psi and Hpsi  ****
      call Parallel3d_np_k(npk)
      call Parallel3d_taskid_k(taskid_k)
      npkrot  = npk/2
      npkeven = (mod(npk,2).eq.0)
      if (.not.npkeven) npkrot = npkrot+1
      replicated = (npk.gt.1)
      if (hfx_on.and.replicated) then
         nbrillq_max = nbrillq
         call Parallel_IMaxAll(nbrillq_max)
         nrsize = 2*nbrillq_max*neall*nfft3d
         psi_rep1_tag  = cpsi_data_alloc(nbrillq_max,neall,2*nfft3d)
         psi_rep2_tag  = cpsi_data_alloc(nbrillq_max,neall,2*nfft3d)
         Hpsi_rep1_tag = cpsi_data_alloc(nbrillq_max,neall,2*nfft3d)
         Hpsi_rep2_tag = cpsi_data_alloc(nbrillq_max,neall,2*nfft3d)
         value = BA_alloc_get(mt_dbl,4*nbrillq_max,'kw_rep1',
     >                        kw_rep1(2),kw_rep1(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,4*nbrillq_max,'kw_rep2',
     >                        kw_rep2(2),kw_rep2(1))
         value = value.and.
     >           BA_alloc_get(mt_int,npk,'nbrillq_rep',
     >                        nbrillq_rep(2),nbrillq_rep(1))
         if (.not. value)   
     >      call errquit('band_init_HFX: out of heap memory',1, MA_ERR)

         do n1=1,npk
            int_mb(nbrillq_rep(1)+n1-1) = 0
         end do
         int_mb(nbrillq_rep(1)+taskid_k) = nbrillq
         call Parallela_Vector_ISumAll(3,npk,int_mb(nbrillq_rep(1)))
      end if

   
      return
      end


*     *************************
*     *                       *
*     *     band_end_HFX      *
*     *                       *
*     *************************
      subroutine band_end_HFX()
      implicit none

#include "bafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer i,ms

*     **** external functions ****
      integer  control_version
      external control_version

      if ((norbs(1)+norbs(2)).gt.0) then

*       **** print out orbital contributions ****
c        if (orb_contribution) then
c           call Parallel_taskid(taskid)
c           if (taskid.eq.MASTER) then
c              write(6,487)
c              write(6,488)
c              do ms=1,ispin
c              do i=1,norbs(ms)
c                write(6,489) 
c     >            ms,i,
c     >            dbl_mb(ehfx_orb(1,ms)+i-1)
c              end do
c              end do
c           end if
c  487   format(//,'== Orbital Contributions to HFX ==')
c  488   format(/1x,'orbital',15x,
c     >         'HF_Exchange') 	
c  489   format(1x,i3,i7,2x,e18.6)
c        end if


*       **** deallocate memory ****
        value = .true.
c        do ms=1,ispin
c          value = value.and.BA_free_heap(ehfx_orb(2,ms))
c        end do

*        **** end coulomb_screened ****
         call c_coulomb_screened_end()

*        **** deallocate replicated space if necessary ****
         if (replicated) then
            call cpsi_data_dealloc(psi_rep1_tag)
            call cpsi_data_dealloc(psi_rep2_tag)
            call cpsi_data_dealloc(Hpsi_rep1_tag)
            call cpsi_data_dealloc(Hpsi_rep2_tag)
            value =           BA_free_heap(kw_rep1(2))
            value = value.and.BA_free_heap(kw_rep2(2))
            value = value.and.BA_free_heap(nbrillq_rep(2))
            if (.not. value)   
     >         call errquit('band_end_HFX: free heap memory',1,MA_ERR)
         end if

      end if

      return
      end

*     *************************
*     *                       *
*     *     band_print_HFX    *
*     *                       *
*     *************************
      subroutine band_print_HFX(unit)
      implicit none
      integer unit

#include "bafdecls.fh"
#include "band_hfx.fh"

*     **** local variables ****
      integer i,ms
      real*8   control_attenuation
      external control_attenuation

      if (hfx_on) then
        if (relaxed) then
          write(unit,1001)
        else
          write(unit,1002)
        end if
        write(unit,1006)
        if (rcut.ge.0.0d0) write(unit,1008) rcut
        if (rcut.ge.0.0d0) write(unit,1009) pp
        if (rcut.ge.0.0d0) write(unit,1011) flag
        if ((rcut.ge.0.0d0).and.(flag.eq.2)) 
     >     write(unit,1012) control_attenuation()
        if (hfx_parameter.ne.1.0d0) write(unit,1010) hfx_parameter
        write(unit,*)
      end if

      return
 1001 FORMAT(6x,"- HFX relaxed")
 1002 FORMAT(6x,"- HFX unrelaxed")
 1003 FORMAT(6x,"- HFX restricted orbitals :",10I5)
 1004 FORMAT(6x,"- HFX alpha orbitals:",10I5)
 1005 FORMAT(6x,"- HFX beta orbitals :",10I5)

 1006 FORMAT(6x,"- HFX screened coulomb solver")
 1007 FORMAT(6x,"- HFX free-space coulomb solver")
 1008 FORMAT(6x,"- HFX screening radius(band:HFX_screening_radius):",
     >       E10.3)
 1009 FORMAT(6x,"- HFX screening power (band:HFX_screening_power) :",
     >       E10.3)
 1010 FORMAT(6x,"- HFX scaling parameter (band:HFX_parameter)     :",
     >       E10.3)
 1011 FORMAT(6x,"- HFX screening type (band:HFX_screening_type)   :",
     >       I2)
 1012 FORMAT(6x,"- attenuation parameter (nwpw:attenuation)       :",
     >       E10.3)
      end



*     ****************************
*     *                    	 *
*     *     band_potential_HFX   *
*     *                          *
*     ****************************
      subroutine band_potential_HFX(ispin0,psi_r_tag,Hpsi_r_tag)
      implicit none
      integer    ispin0
      integer    psi_r_tag
      integer    Hpsi_r_tag

#include "bafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

      integer nb1,nb2
      integer istart,iend,jstart,jend,imodn,imodtask
      integer ms,l,q,n,indx1,indx2,Levels,neq(2)
      integer nbrillq_rep1,r1,r2
      integer request2(4),request3(4),request4(4)
      real*8  kv1(3),kv2(3),w2

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_get_allptr
      external cpsi_data_get_chnk,cpsi_data_get_allptr
      real*8   brillioun_weight,brillioun_k
      external brillioun_weight,brillioun_k
c      integer  Butter_Levels,Dneall_na_ptr
c      external Butter_Levels,Dneall_na_ptr

      call nwpw_timing_start(33)
      ehfx = 0.0d0
      phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

c        **** reduceall algorithm ****
         if (replicated) then

*           **** initialize replicated data ****
            call dcopy(4*nbrillq_max,0.0d0,0,dbl_mb(kw_rep1(1)),1)
            call dcopy(4*nbrillq_max,0.0d0,0,dbl_mb(kw_rep2(1)),1)
            do nb1=1,nbrillq
               dbl_mb(kw_rep1(1)+4*(nb1-1))   = brillioun_k(1,nb1)
               dbl_mb(kw_rep1(1)+4*(nb1-1)+1) = brillioun_k(2,nb1)
               dbl_mb(kw_rep1(1)+4*(nb1-1)+2) = brillioun_k(3,nb1)
               dbl_mb(kw_rep1(1)+4*(nb1-1)+3) = brillioun_weight(nb1)
            end do
            r2 = taskid_k
            nbrillq_rep1 = nbrillq
            call dcopy(nrsize,0.0d0,0,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),1)
            call dcopy(nrsize,0.0d0,0,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep2_tag)),1)
            call dcopy(nrsize,0.0d0,0,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_rep1_tag)),1)
            call dcopy(nrsize,0.0d0,0,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_rep2_tag)),1)

            call dcopy(nbrillq*2*nfft3d*neall,
     >                 dbl_mb(cpsi_data_get_allptr(psi_r_tag)),1,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),1)
            call dcopy(nbrillq*2*nfft3d*neall,
     >                 0.0d0,0,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_r_tag)),1)


            call Parallela_start_rotate(3,1,12,
     >                                 dbl_mb(kw_rep1(1)),4*nbrillq_max,
     >                                 dbl_mb(kw_rep2(1)),4*nbrillq_max,
     >                                 request2)
            call Parallela_start_rotate(3,1,13,
     >           dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),nrsize,
     >           dbl_mb(cpsi_data_get_allptr(psi_rep2_tag)),nrsize,
     >           request3)

            ehfx   = 0.0d0
            phfx   = 0.0d0
            kv1(1) = 0.0d0
            kv1(2) = 0.0d0
            kv1(3) = 0.0d0
            kv2(1) = 0.0d0
            kv2(2) = 0.0d0
            kv2(3) = 0.0d0
            call c_coulomb_screened_set(kv1,kv2)
            do nb1=1,nbrillq
               call band_potential_HFX_sub(brillioun_weight(nb1),
     >                  dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb1)),
     >                  dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)))
            end do
            do nb1=1,nbrillq-1
               do nb2=nb1+1,nbrillq
                  kv1(1) = brillioun_k(1,nb1)
                  kv1(2) = brillioun_k(2,nb1)
                  kv1(3) = brillioun_k(3,nb1)
                  kv2(1) = brillioun_k(1,nb2)
                  kv2(2) = brillioun_k(2,nb2)
                  kv2(3) = brillioun_k(3,nb2)
                  call c_coulomb_screened_set(kv1,kv2)
                  call band_potential_HFX_sub2(brillioun_weight(nb1),
     >                                         brillioun_weight(nb2),
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb2)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb2)))
               end do
            end do

            do r1=1,npkrot-1
               r2 = mod(r2-1+npk,npk)
               nbrillq_rep1 = int_mb(nbrillq_rep(1)+r2)

               call Parallela_end_rotate(request2)
               call dcopy(4*nbrillq_max,
     >                    dbl_mb(kw_rep2(1)),1,dbl_mb(kw_rep1(1)),1)

               call Parallela_end_rotate(request3)
               call dcopy(nrsize,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep2_tag)),1,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),1)

               if (r1.gt.1) then
                  call Parallela_end_rotate(request4)
                  call dcopy(nrsize,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_rep2_tag)),1,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_rep1_tag)),1)
               end if

               if ((r1.lt.(npkrot-1)).or.(npkeven)) then
                  call Parallela_start_rotate(3,1,12,
     >                                 dbl_mb(kw_rep1(1)),4*nbrillq_max,
     >                                 dbl_mb(kw_rep2(1)),4*nbrillq_max,
     >                                 request2)
                  call Parallela_start_rotate(3,1,13,
     >              dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),nrsize,
     >              dbl_mb(cpsi_data_get_allptr(psi_rep2_tag)),nrsize,
     >              request3)
               end if

               do nb1=1,nbrillq
                  do nb2=1,nbrillq_rep1
                     kv1(1) = brillioun_k(1,nb1)
                     kv1(2) = brillioun_k(2,nb1)
                     kv1(3) = brillioun_k(3,nb1)
                     kv2(1) = dbl_mb(kw_rep1(1)+4*(nb2-1))
                     kv2(2) = dbl_mb(kw_rep1(1)+4*(nb2-1)+1)
                     kv2(3) = dbl_mb(kw_rep1(1)+4*(nb2-1)+2)
                     w2     = dbl_mb(kw_rep1(1)+4*(nb2-1)+3)

                     call c_coulomb_screened_set(kv1,kv2)
                     call band_potential_HFX_sub2(
     >                    brillioun_weight(nb1),w2,
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag,nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(psi_rep1_tag,nb2)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_rep1_tag,nb2)))
                  end do
               end do

               if (r1.lt.(npkrot-1)) then
                  call Parallela_start_rotate(3,1,14,
     >              dbl_mb(cpsi_data_get_allptr(Hpsi_rep1_tag)),nrsize,
     >              dbl_mb(cpsi_data_get_allptr(Hpsi_rep2_tag)),nrsize,
     >              request4)
               else
                  call Parallela_start_rotate(3,-r1,14,
     >              dbl_mb(cpsi_data_get_allptr(Hpsi_rep1_tag)),nrsize,
     >              dbl_mb(cpsi_data_get_allptr(Hpsi_rep2_tag)),nrsize,
     >              request4)
               end if

            end do
            if (npkeven) then
               r2 = mod(r2-1+npk,npk)
               nbrillq_rep1 = int_mb(nbrillq_rep(1)+r2)

               call Parallela_end_rotate(request2)
               call dcopy(4*nbrillq_max,
     >                    dbl_mb(kw_rep2(1)),1,dbl_mb(kw_rep1(1)),1)

               call Parallela_end_rotate(request3)
               call dcopy(nrsize,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep2_tag)),1,
     >                 dbl_mb(cpsi_data_get_allptr(psi_rep1_tag)),1)

               do nb1=1,nbrillq
                  do nb2=1,nbrillq_rep1
                     kv1(1) = brillioun_k(1,nb1)
                     kv1(2) = brillioun_k(2,nb1)
                     kv1(3) = brillioun_k(3,nb1)
                     kv2(1) = dbl_mb(kw_rep1(1)+4*(nb2-1))
                     kv2(2) = dbl_mb(kw_rep1(1)+4*(nb2-1)+1)
                     kv2(3) = dbl_mb(kw_rep1(1)+4*(nb2-1)+2)
                     w2     = dbl_mb(kw_rep1(1)+4*(nb2-1)+3)

                     call c_coulomb_screened_set(kv1,kv2)
                     call band_potential_HFX_sub3(
     >                    brillioun_weight(nb1),w2,
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag,nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(psi_rep1_tag,nb2)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)))
                  end do
               end do
            end if

            call Parallela_end_rotate(request4)
            call daxpy(nbrillq*2*nfft3d*neall,1.0d0,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_rep2_tag)),1,
     >                 dbl_mb(cpsi_data_get_allptr(Hpsi_r_tag)),1)

            call Parallel_SumAll(ehfx)
            call Parallel_SumAll(phfx)

         else
            ehfx   = 0.0d0
            phfx   = 0.0d0
            kv1(1) = 0.0d0
            kv1(2) = 0.0d0
            kv1(3) = 0.0d0
            kv2(1) = 0.0d0
            kv2(2) = 0.0d0
            kv2(3) = 0.0d0
            call c_coulomb_screened_set(kv1,kv2)
            do nb1=1,nbrillq
            
               call dcopy(2*nfft3d*neall,0.0d0,0,
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)),1)
               call band_potential_HFX_sub(brillioun_weight(nb1),
     >                  dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb1)),
     >                  dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)))
            end do
            do nb1=1,nbrillq-1
               do nb2=nb1+1,nbrillq
                  kv1(1) = brillioun_k(1,nb1)
                  kv1(2) = brillioun_k(2,nb1)
                  kv1(3) = brillioun_k(3,nb1)
                  kv2(1) = brillioun_k(1,nb2)
                  kv2(2) = brillioun_k(2,nb2)
                  kv2(3) = brillioun_k(3,nb2)
           
                  call c_coulomb_screened_set(kv1,kv2)
                  call band_potential_HFX_sub2(brillioun_weight(nb1),
     >                                         brillioun_weight(nb2),
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(psi_r_tag, nb2)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb1)),
     >                    dbl_mb(cpsi_data_get_chnk(Hpsi_r_tag,nb2)))
               end do
            end do
            call Parallel_SumAll(ehfx)
            call Parallel_SumAll(phfx)
         end if


      end if
      call nwpw_timing_end(33)

      return
      end


*     *************************
*     *                       *
*     *     band_energy_HFX   *
*     *                       *
*     *************************
      subroutine band_energy_HFX(ispin0,psi_r_tag,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      integer psi_r_tag
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

      integer q,n,indx1,indx2

      call nwpw_timing_start(33)
c     **** calculate HFX energy  ****
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
         call band_energy_HFX_sub(ispin0,psi_r_tag,ehfx_out,phfx_out)

c     **** nothing to do ****
      else
         ehfx_out = ehfx
         phfx_out = phfx
      end if
      call nwpw_timing_end(33)

      return
      end

*     ****************************
*     *                          *
*     *  band_potential_HFX_orb  *
*     *                          *
*     ****************************
      subroutine band_potential_HFX_orb(nb1,ms1,psi_r_tag,orb_r,Horb_r)
      implicit none
      integer nb1,ms1,psi_r_tag
      complex*16 orb_r(*),Horb_r(*)

#include "bafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      integer nbq2
      real*8 kv1(3),kv2(3)

*     **** external functions ****
      integer  cpsi_data_get_chnk
      external cpsi_data_get_chnk
      real*8   brillioun_weight,brillioun_k,brillioun_all_k
      external brillioun_weight,brillioun_k,brillioun_all_k

      if ((norbs(ms1).ne.0).and.relaxed) then
         call dcopy(2*nfft3d,0.0d0,0,Horb_r,1)
         kv1(1) = brillioun_all_k(1,nb1)
         kv1(2) = brillioun_all_k(2,nb1)
         kv1(3) = brillioun_all_k(3,nb1)
         do nbq2=1,nbrillq
            kv2(1) = brillioun_k(1,nbq2)
            kv2(2) = brillioun_k(2,nbq2)
            kv2(3) = brillioun_k(3,nbq2)
            call c_coulomb_screened_set(kv1,kv2)
            call band_potential_HFX_orb_sub(ms1,brillioun_weight(nbq2),
     >               dbl_mb(cpsi_data_get_chnk(psi_r_tag,nbq2)),
     >               orb_r,Horb_r)
         end do
         call K1dB_Vector_SumAll(2*nfft3d,Horb_r)

      end if
      return
      end



*     *************************
*     *                       *
*     *     band_HFX          *
*     *                       *
*     *************************
      logical function band_HFX()
      implicit none

#include "band_hfx.fh"

      band_HFX= hfx_on
      return
      end

*     *************************
*     *                       *
*     *   band_HFX_relaxed    *
*     *                       *
*     *************************
      logical function band_HFX_relaxed()
      implicit none

#include "bafdecls.fh"
#include "band_hfx.fh"

      band_hfx_relaxed = relaxed
      return
      end





c***************** sub/replicated routines *****************************

*     ********************************
*     *                    	     *
*     *     band_potential_HFX_sub   *
*     *                              *
*     ********************************
*
* calculates the exact exchange potentials and k=l
*
      subroutine band_potential_HFX_sub(weight,psi_r,Hpsi_r)
      implicit none

#include "band_hfx.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      real*8     weight
      complex*16 psi_r(nfft3d,*)
      complex*16 Hpsi_r(nfft3d,*)

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,ms
      integer dn(2),vij(2),tmp1(2),tmp2(2),index1,index2
      integer psir1,psir2,hpsir1,hpsir2
      integer i1,j1,k1,NN
      integer i2,j2,k2
      integer i3,j3,k3
      integer kv1,kv2
      real*8  scal1,scal2,dv,eh,ph,weight2

*     **** external functions ****
      real*8   lattice_omega,c_icoulomb_screened_e
      logical   C3dB_rc_pfft3_queue_filled,C3dB_cr_pfft3_queue_filled
      external lattice_omega,c_icoulomb_screened_e
      external  C3dB_rc_pfft3_queue_Cilled,C3dB_cr_pfft3_queue_filled

      !ehfx = 0.0d0
      !phfx = 0.0d0
      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

         call C3dB_nx(1,n1)
         call C3dB_ny(1,n2)
         call C3dB_nz(1,n3)
         value = BA_push_get(mt_dcpl,nfft3d,'dn_hfx',dn(2),dn(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'vij_hfx',vij(2),vij(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp1_hfx',tmp1(2),tmp1(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp2_hfx',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('band_potential_HFX:out of stack memory',0,
     >      MA_ERR)

 
         weight2 = weight*weight
         scal1 = 1.0d0/dble(n1*n2*n3)
         scal2 = 1.0d0/lattice_omega()
         dv = scal1/scal2
 
*        ***** screened coulomb solver ****
           do ms=1,ispin
c              call dcopy(norbs(ms),0.0d0,0,dbl_mb(ehfx_orb(1,ms)),1)
              NN = norbs(ms)*(norbs(ms)+1)/2
              i1 = 1
              j1 = 1
              k1 = 1
              i2 = 1
              j2 = 1
              k2 = 1
              i3 = 1
              j3 = 1
              k3 = 1
              done = .false.
              do while (.not.done)

                 if (k1.le.NN) then

*                   **** generate dnij for Vij  ****
                    call C3dB_bb_Mul(1,psi_r(1,j1),
     >                                 psi_r(1,i1),dcpl_mb(dn(1)))
                    call C3dB_b_SMul1(1,scal2*scal1,dcpl_mb(dn(1)))
                    call C3dB_rc_pfft3f_queuein(0,dcpl_mb(dn(1)))
 
                    k1 = k1 + 1
                    j1 = j1 + 1
                    if (j1.gt.i1) then
                       j1 = 1
                       i1 = i1 + 1
                    end if
                 end if
 
                 if (     ((C3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >               .and.(k2.le.NN)) then

                       call C3dB_rc_pfft3f_queueout(0,dcpl_mb(dn(1)))

*                      **** generate Vcoul ****
                       eh = c_icoulomb_screened_e(dcpl_mb(dn(1)))
                       call c_coulomb_screened_v(dcpl_mb(dn(1)),
     >                                           dcpl_mb(vij(1)))

*                      **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
                       eh = eh*hfx_parameter*weight2
                       if (ispin.eq.1) eh = eh + eh
                       ph = 2.0d0*eh
                       ehfx = ehfx - eh
                       phfx = phfx - ph
c                       dbl_mb(ehfx_orb(1,ms)+i2-1) 
c     >                  = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
                       if (i2.ne.j2) then
                          ehfx = ehfx - eh
                          phfx = phfx - ph
c                          dbl_mb(ehfx_orb(1,ms)+i2-1) 
c     >                     = dbl_mb(ehfx_orb(1,ms)+i2-1) - eh
                       end if

                       call C3dB_cr_pfft3b_queuein(0,dcpl_mb(vij(1)))

                       k2 = k2 + 1
                       j2 = j2 + 1
                       if (j2.gt.i2) then
                          j2 = 1
                          i2 = i2 + 1
                       end if
                 end if

                 if ((C3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

                    call C3dB_cr_pfft3b_queueout(0,dcpl_mb(vij(1)))

*                   **** apply hfx_parameter ****
                    call C3dB_b_SMul1(1,hfx_parameter*weight,
     >                                dcpl_mb(vij(1)))

*                   **** generate (Vij)*psi_r ***
                    call C3dB_bb_ncMul(1,dcpl_mb(vij(1)),
     >                                 psi_r(1,j3),
     >                                 dcpl_mb(tmp1(1)))

*                   **** add -(Vij)*psi_r to Hpsi_r ***
                    call C3dB_bb_Sub2(1,dcpl_mb(tmp1(1)),Hpsi_r(1,i3))


                    !**** include off diagonal terms ****
                    if (i3.ne.j3) then
*                      **** generate (Vij)*psi_r ***
                       call C3dB_bb_Mul(1,dcpl_mb(vij(1)),
     >                                 psi_r(1,i3),
     >                                 dcpl_mb(tmp2(1)))

*                      **** add -(Vij)*psi_r to Hpsi_r ***
                       call C3dB_bb_Sub2(1,dcpl_mb(tmp2(1)),
     >                                   Hpsi_r(1,j3))
                    end if

                    k3 = k3 + 1
                    j3 = j3 + 1
                    if (j3.gt.i3) then
                       j3 = 1
                       i3 = i3 + 1
                    end if
                 end if

              done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
           end do !**** while ****

c           call Parallel_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
       end do !**** ms *****

 
         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(vij(2))
         value = value.and.BA_pop_stack(dn(2))
         if (.not. value) 
     >      call errquit('band_potential_HFX:popping stack memory',0,
     >      MA_ERR)

      end if

      return
      end


*     ********************************
*     *                    	     *
*     *     band_potential_HFX_sub2  *
*     *                              *
*     ********************************
*
* calculates the exact exchange potentials and k!=l
*
      subroutine band_potential_HFX_sub2(weight1,weight2,
     >                                   psi1_r,psi2_r,
     >                                   Hpsi1_r,Hpsi2_r)
      implicit none

#include "band_hfx.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      real*8     weight1,weight2
      complex*16  psi1_r(nfft3d,*), psi2_r(nfft3d,*)
      complex*16 Hpsi1_r(nfft3d,*),Hpsi2_r(nfft3d,*)

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,ms
      integer dn12(2),v12(2),tmp1(2),tmp2(2),index1,index2
      integer psir1,psir2,hpsir1,hpsir2
      integer i1,j1,k1,N,NN
      integer i2,j2,k2
      integer i3,j3,k3
      integer kv1,kv2
      real*8  scal1,scal2,dv,eh,ph,weight12

*     **** external functions ****
      real*8   lattice_omega,c_icoulomb_screened_e
      logical   C3dB_rc_pfft3_queue_filled,C3dB_cr_pfft3_queue_filled
      external lattice_omega,c_icoulomb_screened_e
      external  C3dB_rc_pfft3_queue_Cilled,C3dB_cr_pfft3_queue_filled

      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

         call C3dB_nx(1,n1)
         call C3dB_ny(1,n2)
         call C3dB_nz(1,n3)
         value = BA_push_get(mt_dcpl,nfft3d,'dn12_hfx',dn12(2),dn12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'v12_hfx',v12(2),v12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp1_hfx',tmp1(2),tmp1(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp2_hfx',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('band_potential_HFX:out of stack memory',0,
     >      MA_ERR)

 
         weight12 = weight1*weight2
         scal1 = 1.0d0/dble(n1*n2*n3)
         scal2 = 1.0d0/lattice_omega()
         dv = scal1/scal2
 
*        ***** screened coulomb solver ****
           do ms=1,ispin
              N  = norbs(ms)
              NN = N*N
              i1 = 1
              j1 = 1
              k1 = 1
              i2 = 1
              j2 = 1
              k2 = 1
              i3 = 1
              j3 = 1
              k3 = 1
              done = .false.
              do while (.not.done)

                 if (k1.le.NN) then

*                   **** generate dn12 for v12 ****
                    call C3dB_bb_Mul(1,psi1_r(1,j1),
     >                                 psi2_r(1,i1),
     >                                 dcpl_mb(dn12(1)))
                    call C3dB_b_SMul1(1,scal2*scal1,dcpl_mb(dn12(1)))
                    call C3dB_rc_pfft3f_queuein(0,dcpl_mb(dn12(1)))
 
                    k1 = k1 + 1
                    j1 = j1 + 1
                    if (j1.gt.N) then
                       j1 = 1
                       i1 = i1 + 1
                    end if
                 end if
 
                 if (     ((C3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >               .and.(k2.le.NN)) then

                       call C3dB_rc_pfft3f_queueout(0,dcpl_mb(dn12(1)))

*                      **** generate V12 ****
                       eh=2.0d0*c_icoulomb_screened_e(dcpl_mb(dn12(1)))
                       call c_coulomb_screened_v(dcpl_mb(dn12(1)),
     >                                           dcpl_mb(v12(1)))

                       call C3dB_cr_pfft3b_queuein(0,dcpl_mb(v12(1)))

*                      **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
                       eh = eh*hfx_parameter*weight12
                       if (ispin.eq.1) eh = eh + eh
                       ph = 2.0d0*eh
                       ehfx = ehfx - eh
                       phfx = phfx - ph

                       k2 = k2 + 1
                       j2 = j2 + 1
                       if (j2.gt.N) then
                          j2 = 1
                          i2 = i2 + 1
                       end if
                 end if

                 if ((C3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

                    call C3dB_cr_pfft3b_queueout(0,dcpl_mb(v12(1)))

*                   **** generate conjg(V12)*psi1_r and add -conjg(V12)*psi1_r to Hpsi2_r ****
                    call C3dB_bb_ncMul(1,dcpl_mb(v12(1)),psi1_r(1,j3),
     >                                 dcpl_mb(tmp1(1)))
*                   **** apply hfx_parameter ****
                    call C3dB_b_SMul1(1,hfx_parameter*weight1,
     >                                dcpl_mb(tmp1(1)))
                    call C3dB_bb_Sub2(1,dcpl_mb(tmp1(1)),Hpsi2_r(1,i3))


*                   **** generate (V12)*psi2_r and add -(V12)*psi2_r to Hpsi1_r ****
                    call C3dB_bb_Mul(1,dcpl_mb(v12(1)),psi2_r(1,i3),
     >                                 dcpl_mb(tmp2(1)))
*                   **** apply hfx_parameter ****
                    call C3dB_b_SMul1(1,hfx_parameter*weight2,
     >                                dcpl_mb(tmp2(1)))
                    call C3dB_bb_Sub2(1,dcpl_mb(tmp2(1)),Hpsi1_r(1,j3))


                    k3 = k3 + 1
                    j3 = j3 + 1
                    if (j3.gt.N) then
                       j3 = 1
                       i3 = i3 + 1
                    end if
                 end if

              done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
           end do !**** while ****

       end do !**** ms *****

 
         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(v12(2))
         value = value.and.BA_pop_stack(dn12(2))
         if (.not. value) 
     >   call errquit('band_potential_HFX_sub2:popping stack memory',
     >                0,MA_ERR)

      end if

      return
      end

*     ********************************
*     *                    	     *
*     *     band_potential_HFX_sub3  *
*     *                              *
*     ********************************
*
* calculates the exact exchange potentials and k!=l
*
      subroutine band_potential_HFX_sub3(weight1,weight2,
     >                                   psi1_r,psi2_r,
     >                                   Hpsi1_r)
      implicit none

#include "band_hfx.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      real*8     weight1,weight2
      complex*16  psi1_r(nfft3d,*), psi2_r(nfft3d,*)
      complex*16 Hpsi1_r(nfft3d,*)

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,ms
      integer dn12(2),v12(2),tmp2(2)
      integer i1,j1,k1,N,NN
      integer i2,j2,k2
      integer i3,j3,k3
      integer kv1,kv2
      real*8  scal1,scal2,dv,eh,ph,weight12

*     **** external functions ****
      real*8   lattice_omega,c_icoulomb_screened_e
      logical   C3dB_rc_pfft3_queue_filled,C3dB_cr_pfft3_queue_filled
      external lattice_omega,c_icoulomb_screened_e
      external  C3dB_rc_pfft3_queue_Cilled,C3dB_cr_pfft3_queue_filled

      if (((norbs(1)+norbs(2)).ne.0).and.relaxed) then

         call C3dB_nx(1,n1)
         call C3dB_ny(1,n2)
         call C3dB_nz(1,n3)
         value = BA_push_get(mt_dcpl,nfft3d,'dn12_hfx',dn12(2),dn12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'v12_hfx',v12(2),v12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp2_hfx',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('band_potential_HFX:out of stack memory',0,
     >      MA_ERR)

 
         weight12 = weight1*weight2
         scal1 = 1.0d0/dble(n1*n2*n3)
         scal2 = 1.0d0/lattice_omega()
         dv = scal1/scal2
 
*        ***** screened coulomb solver ****
           do ms=1,ispin
              N  = norbs(ms)
              NN = N*N
              i1 = 1
              j1 = 1
              k1 = 1
              i2 = 1
              j2 = 1
              k2 = 1
              i3 = 1
              j3 = 1
              k3 = 1
              done = .false.
              do while (.not.done)

                 if (k1.le.NN) then

*                   **** generate dn12 for v12 ****
                    call C3dB_bb_Mul(1,psi1_r(1,j1),
     >                                 psi2_r(1,i1),
     >                                 dcpl_mb(dn12(1)))
                    call C3dB_b_SMul1(1,scal2*scal1,dcpl_mb(dn12(1)))
                    call C3dB_rc_pfft3f_queuein(0,dcpl_mb(dn12(1)))
 
                    k1 = k1 + 1
                    j1 = j1 + 1
                    if (j1.gt.N) then
                       j1 = 1
                       i1 = i1 + 1
                    end if
                 end if
 
                 if (     ((C3dB_rc_pfft3_queue_filled()).or.(k1.gt.NN))
     >               .and.(k2.le.NN)) then

                       call C3dB_rc_pfft3f_queueout(0,dcpl_mb(dn12(1)))

*                      **** generate V12 ****
                       eh = c_icoulomb_screened_e(dcpl_mb(dn12(1)))
                       call c_coulomb_screened_v(dcpl_mb(dn12(1)),
     >                                           dcpl_mb(v12(1)))

                       call C3dB_cr_pfft3b_queuein(0,dcpl_mb(v12(1)))

*                      **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
                       eh = eh*hfx_parameter*weight12
                       if (ispin.eq.1) eh = eh + eh
                       ph = 2.0d0*eh
                       ehfx = ehfx - eh
                       phfx = phfx - ph

                       k2 = k2 + 1
                       j2 = j2 + 1
                       if (j2.gt.N) then
                          j2 = 1
                          i2 = i2 + 1
                       end if
                 end if

                 if ((C3dB_cr_pfft3_queue_filled()).or.(k2.gt.NN)) then

                    call C3dB_cr_pfft3b_queueout(0,dcpl_mb(v12(1)))

*                   **** generate (V12)*psi2_r and add -(V12)*psi2_r to Hpsi1_r ****
                    call C3dB_bb_Mul(1,dcpl_mb(v12(1)),psi2_r(1,i3),
     >                                 dcpl_mb(tmp2(1)))
*                   **** apply hfx_parameter ****
                    call C3dB_b_SMul1(1,hfx_parameter*weight2,
     >                                dcpl_mb(tmp2(1)))
                    call C3dB_bb_Sub2(1,dcpl_mb(tmp2(1)),Hpsi1_r(1,j3))

                    k3 = k3 + 1
                    j3 = j3 + 1
                    if (j3.gt.N) then
                       j3 = 1
                       i3 = i3 + 1
                    end if
                 end if

              done = (k1.gt.NN).and.(k2.gt.NN).and.(k3.gt.NN)
           end do !**** while ****

       end do !**** ms *****

       value =           BA_pop_stack(tmp2(2))
       value = value.and.BA_pop_stack(v12(2))
       value = value.and.BA_pop_stack(dn12(2))
       if (.not. value) 
     > call errquit('band_potential_HFX_sub3:popping stack memory',
     >              0,MA_ERR)

      end if

      return
      end






*     *****************************
*     *                           *
*     *     band_energy_HFX_sub   *
*     *                           *
*     *****************************
      subroutine band_energy_HFX_sub(ispin0,psi_r_tag,ehfx_out,phfx_out)
      implicit none
      integer ispin0
      integer psi_r_tag
      real*8 ehfx_out
      real*8 phfx_out

#include "bafdecls.fh"
#include "band_hfx.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,j,n1,n2,n3,ms,k1
      integer dn(2),tmp1(2),index1,index2
      real*8  scal1,scal2,dv,eh,ph

*     **** external functions ****
      real*8   lattice_omega,coulomb_screened_e
      external lattice_omega,coulomb_screened_e

      
      if (((norbs(1)+norbs(2)).ne.0).and.(.not.relaxed)) then
        ehfx = 0.0d0
        phfx = 0.0d0

        call C3dB_nx(1,n1)
        call C3dB_ny(1,n2)
        call C3dB_nz(1,n3)
        value = BA_push_get(mt_dcpl,(2*nfft3d),'dn_hfx',dn(2),dn(1))
        value = value.and.
     >          BA_push_get(mt_dcpl,(nfft3d),'tmp1_hfx',tmp1(2),tmp1(1))
        if (.not. value) 
     >     call errquit('band_energy_HFX_sub:out of stack memory',0,
     >       MA_ERR)

        scal1 = 1.0d0/dble(n1*n2*n3)
        scal2 = 1.0d0/lattice_omega()
        dv = scal1/scal2

        k1 = 1
c        do ms=1,ispin
c        do i=1,norbs(ms)
c         dbl_mb(ehfx_orb(1,ms)+i-1) = 0.0d0
c         do j=1,i
c
c            if (mod(k1,npj).eq.taskid_j) then
c              index1 = (i-1)*n2ft3d+1
c              index2 = (j-1)*n2ft3d+1 
c
c*             **** generate dnij ****
c              call D3dB_rr_Mul(1,psi_r(index1),psi_r(index2),
c     >                         dbl_mb(dn(1)))
cc              call D3dB_r_SMul(1,scal2,dbl_mb(dn(1)),dbl_mb(dn(1)))
c              call D3dB_r_SMul1(1,scal2,dbl_mb(dn(1)))
c              call D3dB_r_Zero_Ends(1,dbl_mb(dn(1)))
c
c*             ***** screened coulomb solver ****
c              if (solver_type.eq.1) then
c
c*               **** generate dng ****
cc                call D3dB_r_SMul(1,scal1,dbl_mb(dn(1)),
cc     >                                   dbl_mb(dn(1)))
c                call D3dB_r_SMul1(1,scal1,dbl_mb(dn(1)))
c                call D3dB_rc_pfft3f(1,0,dbl_mb(dn(1)))
c                call Pack_c_pack(0,dbl_mb(dn(1)))
c
c*               **** get Ecoul energy ****
c                eh = coulomb_screened_e(dbl_mb(dn(1)))
c       
c*             ***** free-space coulomb solver ****
c              else
c                 call coulomb2_v(dbl_mb(dn(1)),dbl_mb(tmp1(1)))
c                 call D3dB_rr_dot(1,dbl_mb(dn(1)),dbl_mb(tmp1(1)),eh)
c                 eh = 0.5d0*eh*dv
c              end if
c
c*             **** apply hfx_parameter, double eh for restricted, and calculcate ph ****
c              eh = eh*hfx_parameter
c              if (ispin0.eq.1) eh = eh + eh
c              ph = 2.0d0*eh
c
c              ehfx = ehfx - eh
c              phfx = phfx - ph
c              dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)-eh
c
c              !**** include off diagonal terms ****
c              if (i.ne.j) then
c                 ehfx = ehfx - eh
c                 phfx = phfx - ph
c                 dbl_mb(ehfx_orb(1,ms)+i-1) = dbl_mb(ehfx_orb(1,ms)+i-1)
c     >                                      - eh
c              end if
c
c            end if
c            k1 = k1 + 1
c
c         end do
c        end do
c        call D1dB_Vector_SumAll(norbs(ms),dbl_mb(ehfx_orb(1,ms)))
c        end do

        value =           BA_pop_stack(tmp1(2))
        value = value.and.BA_pop_stack(dn(2))
        if (.not. value) 
     >     call errquit('band_energy_HFX_sub:popping stack memory',0,
     >       MA_ERR)

c        call D1dB_SumAll(ehfx)
c        call D1dB_SumAll(phfx)
      end if

      ehfx_out = ehfx
      phfx_out = phfx

      return
      end



      subroutine band_potential_HFX_orb_sub(ms,weight1,
     >                                      psi1_r,orb_r,Horb_r)
      implicit none

#include "band_hfx.fh"
#include "bafdecls.fh"
#include "errquit.fh"

      integer ms
      real*8     weight1
      complex*16  psi1_r(nfft3d,*)
      complex*16 orb_r(nfft3d),Horb_r(nfft3d)

*     **** local variables ****
      logical value,done
      integer i,j,n1,n2,n3,msshift
      integer dn12(2),v12(2),tmp1(2)
      integer j1,j2,j3,k1,k2,k3,N
      real*8  scal1,scal2

*     **** external functions ****
      real*8   lattice_omega,c_icoulomb_screened_e
      logical   C3dB_rc_pfft3_queue_filled,C3dB_cr_pfft3_queue_filled
      external lattice_omega,c_icoulomb_screened_e
      external  C3dB_rc_pfft3_queue_Cilled,C3dB_cr_pfft3_queue_filled

      if ((norbs(ms).gt.0).and.relaxed) then

         call C3dB_nx(1,n1)
         call C3dB_ny(1,n2)
         call C3dB_nz(1,n3)
         value = BA_push_get(mt_dcpl,nfft3d,'dn12_hfx',dn12(2),dn12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'v12_hfx',v12(2),v12(1))
         value = value.and.
     >           BA_push_get(mt_dcpl,nfft3d,'tmp1_hfx',tmp1(2),tmp1(1))
         if (.not. value) 
     >      call errquit('band_potential_HFX:out of stack memory',0,
     >      MA_ERR)
 
         scal1 = 1.0d0/dble(n1*n2*n3)
         scal2 = 1.0d0/lattice_omega()
 
*        ***** screened coulomb solver ****
         msshift = (ms-1)*norbs(1)
         N  = norbs(ms)
         j1 = 1
         j2 = 1
         j3 = 1
         k1 = 1
         k2 = 1
         k3 = 1
         done = .false.
         do while (.not.done)

            if (k1.le.N) then

*              **** generate dn12 for v12 ****
               call C3dB_bb_Mul(1,psi1_r(1,j1+msshift),
     >                            orb_r,
     >                            dcpl_mb(dn12(1)))
               call C3dB_b_SMul1(1,scal2*scal1,dcpl_mb(dn12(1)))
               call C3dB_rc_pfft3f_queuein(0,dcpl_mb(dn12(1)))
               k1 = k1 + 1
               j1 = j1 + 1
            end if
 
            if (     ((C3dB_rc_pfft3_queue_filled()).or.(k1.gt.N))
     >          .and.(k2.le.N)) then

               call C3dB_rc_pfft3f_queueout(0,dcpl_mb(dn12(1)))

*              **** generate V12 ****
               call c_coulomb_screened_v(dcpl_mb(dn12(1)),
     >                                   dcpl_mb(v12(1)))
               call C3dB_cr_pfft3b_queuein(0,dcpl_mb(v12(1)))

               k2 = k2 + 1
               j2 = j2 + 1
            end if

            if ((C3dB_cr_pfft3_queue_filled()).or.(k2.gt.N)) then

               call C3dB_cr_pfft3b_queueout(0,dcpl_mb(v12(1)))

*              **** generate conjg(V12)*psi1_r and add -conjg(V12)*psi1_r to Horb_r ****
               call C3dB_bb_ncMul(1,dcpl_mb(v12(1)),
     >                            psi1_r(1,j3+msshift),
     >                            dcpl_mb(tmp1(1)))

*              **** apply hfx_parameter ****
               call C3dB_b_SMul1(1,hfx_parameter*weight1,
     >                                dcpl_mb(tmp1(1)))
               call C3dB_bb_Sub2(1,dcpl_mb(tmp1(1)),Horb_r)

               k3 = k3 + 1
               j3 = j3 + 1
            end if

            done = (k1.gt.N).and.(k2.gt.N).and.(k3.gt.N)
         end do !**** while ****
 
         value =           BA_pop_stack(tmp1(2))
         value = value.and.BA_pop_stack(v12(2))
         value = value.and.BA_pop_stack(dn12(2))
         if (.not. value) 
     >   call errquit('band_potential_HFX_orb_sub:popping stack memory',
     >                0,MA_ERR)

      end if

      return
      end

