c
c $Id$
c

*     **********************************
*     *	                               *
*     *         pspw_2qmmm_init        *
*     *                                *
*     **********************************

      subroutine pspw_2qmmm_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
#include "pspw_2qmmm.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,oprint
      integer i,j,ii,jj,kk,ll,mm,ia,ja,ni,nj
      integer frag_size,shake_ptr,indx_ptr
      real*8  volume,f,pi
      character*80 rtdb_name
      integer natmp,nbtmp,nabsize,nindxsize,ndsqsize
      integer nbindxsize,nKr0size
      integer naindxsize,nQr0size
      integer ncbindxsize,ncKr0size
      integer nmbindxsize,nmKr0size
      integer ndindxsize,ndKr0size
      integer nddindxsize,nddKr0size
      integer ncnindxsize0,ncnindxsize1,ncnindxsize2,ncnKr0size,n1,n2
      integer nbbindxsize0,nbbindxsize1,nbbKr0size


      integer ncons_hb
      character*255 namespace,tag

*     **** external functions ****
      logical     control_print
      integer     ion_nion,ion_nion_qm,ion_nion_mm,control_version
      integer     ion_rion_ptr,inp_strlen
      character*7 c_index_name
      real*8      lattice_unita,pspw_qmmm_dihedral_value
      external    control_print
      external    ion_nion,ion_nion_qm,ion_nion_mm,control_version
      external    ion_rion_ptr,inp_strlen
      external    c_index_name
      external    lattice_unita,pspw_qmmm_dihedral_value


*     **** reading number of frags ****
      rtdb_name = 'pspw_2qmmm_nfrag'
      if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nfrag)) nfrag = 0



cccc      qmmm_found  =  (ion_nion_mm().gt.0).and.(nfrag.gt.0)
      qmmm_found  =  (nfrag.gt.0)
      shake_found = .false.
      lmbda_flag  = .false.
      qmmm_lmbda  = 1.0d0

*     **** don't do anything if there are no MM atoms ****
      if (qmmm_found) then


*        **** reading qmmm_lmbda ****
         rtdb_name = 'pspw_2qmmm_lmbda'
         lmbda_flag = .false.
         qmmm_lmbda=1.0d0
         if (btdb_get(rtdb,rtdb_name,mt_dbl,1,qmmm_lmbda))
     >      lmbda_flag=.true.


*        **** reading number of kinds of frags ****
         rtdb_name = 'pspw_2qmmm_nkfrag'
         value = btdb_get(rtdb,rtdb_name,mt_int,1,nkfrag)
         if (.not. value)
     >   call errquit(
     >   'pspw_2qmmm_init:failed reading nkfrag',0,RTDB_ERR)

         value = BA_alloc_get(mt_dbl,3*ion_nion(),'fion2tmp',
     >                        fion2tmp(2),fion2tmp(1))
         value = BA_alloc_get(mt_int,nfrag,'index_nonfrag_start',
     >                        indx_frag_start(2),indx_frag_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nfrag,'size_nonfrag',
     >                        size_frag(2),size_frag(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nfrag,'nonkfrag',kfrag(2),kfrag(1))

         value = value.and.
     >           BA_alloc_get(mt_log,nkfrag,'nonself_interaction',
     >                        self_interaction(2),self_interaction(1))
         value = value.and.
     >           BA_alloc_get(mt_log,nkfrag,'nonincell_frag',
     >                        incell_frag(2),incell_frag(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nkfrag,'nonswitch_Rin',
     >                        switch_Rin(2),switch_Rin(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nkfrag,'nonswitch_Rout',
     >                        switch_Rout(2),switch_Rout(1))

*        **** shake variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnshake',
     >                        nshake(2),nshake(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnab_shake_start',
     >                        nab_shake_start(2),nab_shake_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnindx_shake_start',
     >                        nindx_shake_start(2),nindx_shake_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonndsq_shake_start',
     >                        ndsq_shake_start(2),ndsq_shake_start(1))

*        **** bond variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnbond',
     >                        nbond(2),nbond(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnindx_bond_start',
     >                        nindx_bond_start(2),nindx_bond_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnKr0_bond_start',
     >                        nKr0_bond_start(2),nKr0_bond_start(1))

*        **** angle variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnangle',
     >                        nangle(2),nangle(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnindx_angle_start',
     >                        nindx_angle_start(2),nindx_angle_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnKr0_angle_start',
     >                        nKr0_angle_start(2),nKr0_angle_start(1))

*        **** cbond variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonncbond',
     >                        ncbond(2),ncbond(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnindx_cbond_start',
     >                        nindx_cbond_start(2),nindx_cbond_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnKr0_cbond_start',
     >                        nKr0_cbond_start(2),nKr0_cbond_start(1))

*        **** mbond variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnmbond',
     >                        nmbond(2),nmbond(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnindx_mbond_start',
     >                        nindx_mbond_start(2),nindx_mbond_start(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnKr0_mbond_start',
     >                        nKr0_mbond_start(2),nKr0_mbond_start(1))


*        **** dihedral spring variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonndihedral_spring',
     >                        ndihedral_spring(2),ndihedral_spring(1))
         value = value.and.
     >     BA_alloc_get(mt_int,nkfrag,'nonnindx_dihedral_spring_start',
     >                        nindx_dihedral_spring_start(2),
     >                        nindx_dihedral_spring_start(1))
         value = value.and.
     >      BA_alloc_get(mt_int,nkfrag,'nonnKr0_dihedral_spring_start',
     >                      nKr0_dihedral_spring_start(2),
     >                      nKr0_dihedral_spring_start(1))

*        **** dihedral variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonndihedral',
     >                        ndihedral(2),ndihedral(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx_dihedral_start',
     >                        nindx_dihedral_start(2),
     >                        nindx_dihedral_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnKr0_dihedral_start',
     >                      nKr0_dihedral_start(2),
     >                      nKr0_dihedral_start(1))

*        **** coord variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonncoord',
     >                        ncoord(2),ncoord(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx0_coord_start',
     >                        nindx0_coord_start(2),
     >                        nindx0_coord_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx1_coord_start',
     >                        nindx1_coord_start(2),
     >                        nindx1_coord_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx2_coord_start',
     >                        nindx2_coord_start(2),
     >                        nindx2_coord_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnKr0_coord_start',
     >                      nKr0_coord_start(2),
     >                      nKr0_coord_start(1))


*        **** bondings variables ****
         value = value.and.
     >           BA_alloc_get(mt_int,nkfrag,'nonnbondings',
     >                        nbondings(2),nbondings(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx0_bondings_start',
     >                        nindx0_bondings_start(2),
     >                        nindx0_bondings_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnindx1_bondings_start',
     >                        nindx1_bondings_start(2),
     >                        nindx1_bondings_start(1))
         value = value.and.
     >         BA_alloc_get(mt_int,nkfrag,'nonnKr0_bondings_start',
     >                      nKr0_bondings_start(2),
     >                      nKr0_bondings_start(1))

         if (.not.value) 
     >    call errquit('cannot allocate heap memory for 2qmmm',0,
     >       MA_ERR)


*        **** set index_start and size_frag****
         ii = 0
         jj = 0
         nabsize   = 0
         nindxsize = 0
         ndsqsize  = 0
         nbindxsize = 0
         nKr0size   = 0
         naindxsize = 0
         nQr0size   = 0
         ncbindxsize = 0
         ncKr0size   = 0
         nmbindxsize = 0
         nmKr0size   = 0
         ndindxsize = 0
         ndKr0size   = 0
         nddindxsize = 0
         nddKr0size   = 0
         ncnindxsize0 = 0
         ncnindxsize1 = 0
         ncnindxsize2 = 0
         ncnKr0size   = 0
         nbbindxsize0 = 0
         nbbindxsize1 = 0
         nbbKr0size   = 0

         do ia=1,nkfrag

            rtdb_name = 'pspw_2qmmm_frag_size:'//c_index_name(ia)
            value = btdb_get(rtdb,rtdb_name,mt_int,1,frag_size)
            if (.not. value) frag_size = 1

            rtdb_name='pspw_2qmmm_frag_nindex_start:'//c_index_name(ia)
            value = btdb_get(rtdb,rtdb_name,mt_int,1,ni)
            if (.not. value) ni = 1
            
            do i=1,ni
               int_mb(size_frag(1)+ii+i-1) = frag_size
               int_mb(kfrag(1)+ii+i-1)     = ia
            end do

            rtdb_name = 'pspw_2qmmm_frag_index_start:'//c_index_name(ia)
            value = btdb_get(rtdb,rtdb_name,mt_int,
     >                       ni,int_mb(indx_frag_start(1)+ii))
            if (.not. value) int_mb(indx_frag_start(1)+ii) = 1
            ii = ii + ni

            rtdb_name = 'pspw_2qmmm_frag_self_interaction:'
     >                  //c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_log,1,
     >                       log_mb(self_interaction(1)+ia-1)))
     >          log_mb(self_interaction(1)+ia-1) = .false.

            rtdb_name = 'pspw_2qmmm_frag_incell:'
     >                  //c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_log,1,
     >                       log_mb(incell_frag(1)+ia-1)))
     >          log_mb(incell_frag(1)+ia-1) = .true.

            rtdb_name = 'pspw_2qmmm_frag_switch_Rin:'//c_index_name(ia)
            value = btdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rin(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rin(1)+ia-1) =(2.0160d0/0.529177d0)

            rtdb_name = 'pspw_2qmmm_frag_switch_Rout:'//c_index_name(ia)
            value = btdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rout(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rout(1)+ia-1) =(3.1287d0/0.529177d0)


*           **** shake stuff  - get nshake  ****
            int_mb(nab_shake_start(1)+ia-1) = nabsize + 1
            rtdb_name = 'pspw_2qmmm_frag_nshk:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nshake(1)+ia-1))) 
     >         int_mb(nshake(1)+ia-1) = 0
            nabsize = nabsize + int_mb(nshake(1)+ia-1)


*           **** shake stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_shake_start(1)+ia-1) = nindxsize + 1
            int_mb(ndsq_shake_start(1) +ia-1) = ndsqsize  + 1

            do ja=1,int_mb(nshake(1)+ia-1)
               rtdb_name = 'pspw_2qmmm_frag_na:'
     >                     //c_index_name(ia)//c_index_name(ja)
               if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,natmp))
     >            natmp = 0
               nindxsize = nindxsize + natmp

               rtdb_name = 'pspw_2qmmm_frag_nb:'
     >                     //c_index_name(ia)//c_index_name(ja)
               if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nbtmp))
     >            nbtmp = 0
               ndsqsize = ndsqsize + nbtmp
            end do

*           **** bond stuff  - get nbond  ****
            rtdb_name = 'pspw_2qmmm_frag_nbond:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nbond(1)+ia-1))) 
     >         int_mb(nbond(1)+ia-1) = 0

*           **** bond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_bond_start(1)+ia-1) = nbindxsize + 1
            int_mb(nKr0_bond_start(1) +ia-1) = nKr0size  + 1
            nbindxsize = nbindxsize + 2*int_mb(nbond(1)+ia-1)
            nKr0size   = nKr0size   + 2*int_mb(nbond(1)+ia-1)


*           **** angle stuff  - get nangle  ****
            rtdb_name = 'pspw_2qmmm_frag_nangle:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nangle(1)+ia-1)))
     >         int_mb(nangle(1)+ia-1) = 0

*           **** angle stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_angle_start(1)+ia-1) = naindxsize + 1
            int_mb(nKr0_angle_start(1) +ia-1) = nQr0size  + 1
            naindxsize = naindxsize + 3*int_mb(nangle(1)+ia-1)
            nQr0size   = nQr0size   + 2*int_mb(nangle(1)+ia-1)

*           **** cbond stuff  - get ncbond  ****
            rtdb_name = 'pspw_2qmmm_frag_ncbond:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ncbond(1)+ia-1))) 
     >         int_mb(ncbond(1)+ia-1) = 0

*           **** cbond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_cbond_start(1)+ia-1) = ncbindxsize + 1
            int_mb(nKr0_cbond_start(1) +ia-1) = ncKr0size  + 1
            ncbindxsize = ncbindxsize + 4*int_mb(ncbond(1)+ia-1)
            ncKr0size   = ncKr0size   + 3*int_mb(ncbond(1)+ia-1)

*           **** mbond stuff  - get nmbond  ****
            rtdb_name = 'pspw_2qmmm_frag_nmbond:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nmbond(1)+ia-1))) 
     >         int_mb(nmbond(1)+ia-1) = 0

*           **** mbond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_mbond_start(1)+ia-1) = nmbindxsize + 1
            int_mb(nKr0_mbond_start(1) +ia-1) = nmKr0size  + 1
            nmbindxsize = nmbindxsize + 3*int_mb(nmbond(1)+ia-1)
            nmKr0size   = nmKr0size   + 6*int_mb(nmbond(1)+ia-1)

*           **** dihedral_spring stuff  - get ndihedral_spring  ****
            rtdb_name =
     >        'pspw_2qmmm_frag_ndihedral_spring:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ndihedral_spring(1)+ia-1)))
     >         int_mb(ndihedral_spring(1)+ia-1) = 0

*           **** dihedral_spring stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_dihedral_spring_start(1)+ia-1) = ndindxsize + 1
            int_mb(nKr0_dihedral_spring_start(1) +ia-1) = ndKr0size  + 1
            ndindxsize = ndindxsize + 4*int_mb(ndihedral_spring(1)+ia-1)
            ndKr0size  = ndKr0size  + 2*int_mb(ndihedral_spring(1)+ia-1)

*           **** dihedral stuff  - get ndihedral  ****
            rtdb_name =
     >        'pspw_2qmmm_frag_ndihedral:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ndihedral(1)+ia-1)))
     >         int_mb(ndihedral(1)+ia-1) = 0

*           **** dihedral stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_dihedral_start(1)+ia-1) = nddindxsize + 1
            int_mb(nKr0_dihedral_start(1) +ia-1) = nddKr0size  + 1
            nddindxsize = nddindxsize + 4*int_mb(ndihedral(1)+ia-1)
            nddKr0size  = nddKr0size  + 3*int_mb(ndihedral(1)+ia-1)


*           **** coordination number stuff  - get ncoord  ****
            rtdb_name =
     >        'pspw_2qmmm_frag_ncoord:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ncoord(1)+ia-1)))
     >         int_mb(ncoord(1)+ia-1) = 0

*           **** coordination number stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx0_coord_start(1)+ia-1) = ncnindxsize0 + 1
            int_mb(nindx1_coord_start(1)+ia-1) = ncnindxsize1 + 1
            int_mb(nindx2_coord_start(1)+ia-1) = ncnindxsize2 + 1
            int_mb(nKr0_coord_start(1) +ia-1)  = ncnKr0size  + 1
            do i=1,int_mb(ncoord(1)+ia-1)
               rtdb_name = 'pspw_2qmmm_frag_coord_indxsize1:'
     >                     //c_index_name(ia)//c_index_name(i)
               if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n1)) n1 =0
               rtdb_name = 'pspw_2qmmm_frag_coord_indxsize2:'
     >                     //c_index_name(ia)//c_index_name(i)
               if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n2)) n2 =0
               ncnindxsize0 = ncnindxsize0 + 1
               ncnindxsize1 = ncnindxsize1 + n1
               ncnindxsize2 = ncnindxsize2 + n2
            end do
            ncnKr0size  = ncnKr0size  + 6*int_mb(ncoord(1)+ia-1)


*           **** bondings stuff  - get nbondings ****
            rtdb_name =
     >        'pspw_2qmmm_frag_nbondings:'//c_index_name(ia)
            if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nbondings(1)+ia-1)))
     >         int_mb(nbondings(1)+ia-1) = 0

*           **** bondings stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx0_bondings_start(1)+ia-1) = nbbindxsize0 + 1
            int_mb(nindx1_bondings_start(1)+ia-1) = nbbindxsize1 + 1
            int_mb(nKr0_bondings_start(1) +ia-1) = nbbKr0size  + 1
            do i=1,int_mb(nbondings(1)+ia-1)
               rtdb_name = 'pspw_2qmmm_frag_indxsize1_bondings:'
     >                     //c_index_name(ia)//c_index_name(i)
               if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n1)) n1 =0
               nbbindxsize0 = nbbindxsize0 + 1
               nbbindxsize1 = nbbindxsize1 + 2*n1
               nbbKr0size  = nbbKr0size    + n1 + 2
            end do

         end do


*        *** more shake ****
         shake_found = (nabsize.gt.0)
         if (shake_found) then

            value = BA_alloc_get(mt_int,nabsize,'nonna',na(2),na(1))
            value = value.and.
     >              BA_alloc_get(mt_int,nabsize,'nonnb',nb(2),nb(1))
            value = value.and.
     >              BA_alloc_get(mt_int,nindxsize,'nonindx_shake',
     >                        indx_shake(2),indx_shake(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,ndsqsize,'nondsq_shake',
     >                        dsq_shake(2),dsq_shake(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm shake',0,MA_ERR)

            nabsize   = 0
            nindxsize = 0
            ndsqsize  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nshake(1)+ia-1)
                  rtdb_name = 'pspw_2qmmm_frag_na:'
     >                        //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,natmp))
     >               natmp = 0
                  int_mb(na(1)+nabsize) = natmp
                  rtdb_name = 'pspw_2qmmm_frag_nb:'
     >                        //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,nbtmp))
     >               nbtmp = 0
                  int_mb(nb(1)+nabsize) = nbtmp
                  nabsize = nabsize + 1

                  rtdb_name = 'pspw_2qmmm_frag_indx_shake:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,natmp,
     >                int_mb(indx_shake(1)+nindxsize)))
     >               int_mb(indx_shake(1)+nindxsize) = 0
                  nindxsize = nindxsize + natmp

                  rtdb_name = 'pspw_2qmmm_frag_dsq_shake:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,nbtmp,
     >                dbl_mb(dsq_shake(1)+ndsqsize)))
     >               dbl_mb(dsq_shake(1)+ndsqsize) = 0.0d0
                  ndsqsize = ndsqsize + nbtmp

               end do
            end do
            call pspw_2qmmm_ion_add_constraint()
         end if

*        *** more bond ****
         bond_found = (nbindxsize.gt.0)
         if (bond_found) then
            value = BA_alloc_get(mt_int,nbindxsize,'nonindx_bond',
     >                        indx_bond(2),indx_bond(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nKr0size,'nonKr0_bond',
     >                        Kr0_bond(2),Kr0_bond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm bond',0,MA_ERR)

            nbindxsize = 0
            nKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nbond(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_bond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,2,
     >                int_mb(indx_bond(1)+nbindxsize)))
     >               int_mb(indx_bond(1)+nbindxsize) = 0
                  nbindxsize = nbindxsize + 2

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_bond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                dbl_mb(Kr0_bond(1)+nKr0size)))
     >               dbl_mb(Kr0_bond(1)+nKr0size) = 0.0d0
                  nKr0size = nKr0size + 2

               end do
            end do

         end if



*        *** more angle ****
         angle_found = (naindxsize.gt.0)
         if (angle_found) then
            value = BA_alloc_get(mt_int,naindxsize,'nonindx_angle',
     >                        indx_angle(2),indx_angle(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nQr0size,'nonKr0_angle',
     >                        Kr0_angle(2),Kr0_angle(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm angle',0,MA_ERR)

            naindxsize = 0
            nQr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nangle(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_angle:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,3,
     >                int_mb(indx_angle(1)+naindxsize)))
     >               int_mb(indx_angle(1)+naindxsize) = 0
                  naindxsize = naindxsize + 3

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_angle:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                dbl_mb(Kr0_angle(1)+nQr0size)))
     >               dbl_mb(Kr0_angle(1)+nQr0size) = 0.0d0
                  nQr0size = nQr0size + 2

               end do
            end do

         end if

*        *** more cbond ****
         cbond_found = (ncbindxsize.gt.0)
         if (cbond_found) then
            value = BA_alloc_get(mt_int,ncbindxsize,'nonindx_cbond',
     >                        indx_cbond(2),indx_cbond(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,ncKr0size,'nonKr0_cbond',
     >                        Kr0_cbond(2),Kr0_cbond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm cbond',0,MA_ERR)

            ncbindxsize = 0
            ncKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ncbond(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_cbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,4,
     >                int_mb(indx_cbond(1)+ncbindxsize)))
     >               int_mb(indx_cbond(1)+ncbindxsize) = 0
                  ncbindxsize = ncbindxsize + 4

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_cbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,3,
     >                dbl_mb(Kr0_cbond(1)+ncKr0size)))
     >               dbl_mb(Kr0_cbond(1)+ncKr0size) = 0.0d0
                  ncKr0size = ncKr0size + 3

               end do
            end do

         end if

*        *** more mbond ****
         mbond_found = (nmbindxsize.gt.0)
         if (mbond_found) then
            value = BA_alloc_get(mt_int,nmbindxsize,'nonindx_mbond',
     >                        indx_mbond(2),indx_mbond(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nmKr0size,'nonKr0_mbond',
     >                        Kr0_mbond(2),Kr0_mbond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm mbond',0,MA_ERR)

            nmbindxsize = 0
            nmKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nmbond(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_mbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,3,
     >                int_mb(indx_mbond(1)+nmbindxsize)))
     >               int_mb(indx_mbond(1)+nmbindxsize) = 0
                  nmbindxsize = nmbindxsize + 3

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_mbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,6,
     >                dbl_mb(Kr0_mbond(1)+nmKr0size)))
     >               dbl_mb(Kr0_mbond(1)+nmKr0size) = 0.0d0
                  nmKr0size = nmKr0size + 6

               end do
            end do

         end if

*        *** more dihedral_spring ****
         dihedral_spring_found = (ndindxsize.gt.0)
         if (dihedral_spring_found) then
            value = BA_alloc_get(mt_int,ndindxsize,
     >              'nonindx_dihedral_spring',
     >               indx_dihedral_spring(2),
     >               indx_dihedral_spring(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,ndKr0size,
     >                           'nonKr0_dihedral_spring',
     >                            Kr0_dihedral_spring(2),
     >                            Kr0_dihedral_spring(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm dihedral_spring',
     >       0,MA_ERR)

            ndindxsize = 0
            ndKr0size  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ndihedral_spring(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_dihedral_spring:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,4,
     >               int_mb(indx_dihedral_spring(1)+ndindxsize)))
     >               int_mb(indx_dihedral_spring(1)+ndindxsize) = 0
                  ndindxsize = ndindxsize + 4

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_dihedral_spring:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,2,
     >               dbl_mb(Kr0_dihedral_spring(1)+ndKr0size))) then
                     dbl_mb(Kr0_dihedral_spring(1)+ndKr0size) = 0.001d0
                     dbl_mb(Kr0_dihedral_spring(1)+ndKr0size+1) 
     >                 = pspw_qmmm_dihedral_value(
     >                    int_mb(indx_dihedral_spring(1)+ndindxsize-4),
     >                    dbl_mb(ion_rion_ptr()))
                  end if
                  if (dbl_mb(Kr0_dihedral_spring(1)+ndKr0size+1)
     >               .gt.(10.0d0)) then
                     dbl_mb(Kr0_dihedral_spring(1)+ndKr0size+1) 
     >                 = pspw_qmmm_dihedral_value(
     >                    int_mb(indx_dihedral_spring(1)+ndindxsize-4),
     >                    dbl_mb(ion_rion_ptr()))
                  end if
                  ndKr0size = ndKr0size + 2

               end do
            end do
         end if

*        *** more dihedral ****
         dihedral_found = (nddindxsize.gt.0)
         if (dihedral_found) then
            value = BA_alloc_get(mt_int,nddindxsize,
     >              'nonindx_dihedral',
     >               indx_dihedral(2),
     >               indx_dihedral(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nddKr0size,
     >                           'nonKr0_dihedral',
     >                            Kr0_dihedral(2),
     >                            Kr0_dihedral(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm dihedral',
     >       0,MA_ERR)

            nddindxsize = 0
            nddKr0size  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ndihedral(1)+ia-1)

                  rtdb_name = 'pspw_2qmmm_frag_indx_dihedral:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,4,
     >               int_mb(indx_dihedral(1)+nddindxsize)))
     >               int_mb(indx_dihedral(1)+nddindxsize) = 0
                  nddindxsize = nddindxsize + 4

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_dihedral:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,3,
     >               dbl_mb(Kr0_dihedral(1)+nddKr0size)))
     >               dbl_mb(Kr0_dihedral(1)+nddKr0size) = 0.0d0
                  nddKr0size = nddKr0size + 3

               end do
            end do
         end if

*        *** more coordination number ****
         coord_found = (ncnindxsize1.gt.0)
         if (coord_found) then
            value = BA_alloc_get(mt_int,ncnindxsize0,
     >              'nonnsize1_coord',
     >               nsize1_coord(2),
     >               nsize1_coord(1))
            value = value.and.BA_alloc_get(mt_int,ncnindxsize0,
     >              'nonnsize2_coord',
     >               nsize2_coord(2),
     >               nsize2_coord(1))
            value = value.and.BA_alloc_get(mt_int,ncnindxsize1,
     >              'nonindx1_coord',
     >               indx1_coord(2),
     >               indx1_coord(1))
            value = value.and.BA_alloc_get(mt_int,ncnindxsize2,
     >              'nonindx2_coord',
     >               indx2_coord(2),
     >               indx2_coord(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,ncnKr0size,
     >                           'nonKr0_coord',
     >                            Kr0_coord(2),
     >                            Kr0_coord(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm coord',
     >       0,MA_ERR)


            ncnindxsize0 = 0
            ncnindxsize1 = 0
            ncnindxsize2 = 0
            ncnKr0size  = 0 
            do ia=1,nkfrag
               do ja=1,int_mb(ncoord(1)+ia-1)
                  rtdb_name = 'pspw_2qmmm_frag_coord_indxsize1:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n1)) n1=0
                  rtdb_name = 'pspw_2qmmm_frag_coord_indxsize2:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n2)) n2=0
                  int_mb(nsize1_coord(1)+ncnindxsize0) = n1
                  int_mb(nsize2_coord(1)+ncnindxsize0) = n2
                  ncnindxsize0 = ncnindxsize0 + 1

                  rtdb_name = 'pspw_2qmmm_frag_coord_indx1:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,n1,
     >               int_mb(indx1_coord(1)+ncnindxsize1)))
     >               int_mb(indx1_coord(1)+ncnindxsize1) = 0

                  rtdb_name = 'pspw_2qmmm_frag_coord_indx2:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,n2,
     >               int_mb(indx2_coord(1)+ncnindxsize2)))
     >               int_mb(indx2_coord(1)+ncnindxsize2) = 0
                  ncnindxsize1 = ncnindxsize1 + n1
                  ncnindxsize2 = ncnindxsize2 + n2

                  rtdb_name = 'pspw_2qmmm_frag_coord_Kr0:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,6,
     >               dbl_mb(Kr0_coord(1)+ncnKr0size)))
     >               dbl_mb(Kr0_coord(1)+ncnKr0size) = 0.0d0
                  ncnKr0size = ncnKr0size + 6
               end do
            end do
         end if

*        *** more bondings ****
         bondings_found = (nbbindxsize1.gt.0)
         if (bondings_found) then
            value = BA_alloc_get(mt_int,nbbindxsize0,
     >              'nonnsize1_bondings',
     >               nsize1_bondings(2),
     >               nsize1_bondings(1))
            value = value.and.BA_alloc_get(mt_int,nbbindxsize1,
     >              'nonindx1_bondings',
     >               indx1_bondings(2),
     >               indx1_bondings(1))
            value = value.and.
     >              BA_alloc_get(mt_dbl,nbbKr0size,
     >                           'nonKr0_bondings',
     >                            Kr0_bondings(2),
     >                            Kr0_bondings(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for 2qmmm bondings',
     >       0,MA_ERR)

            nbbindxsize0 = 0
            nbbindxsize1 = 0
            nbbKr0size  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nbondings(1)+ia-1)
                  rtdb_name = 'pspw_2qmmm_frag_indxsize1_bondings:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,n1)) n1=0
                  int_mb(nsize1_bondings(1)+nbbindxsize0) = n1
                  nbbindxsize0 = nbbindxsize0 + 1

                  rtdb_name = 'pspw_2qmmm_frag_indx1_bondings:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_int,2*n1,
     >               int_mb(indx1_bondings(1)+nbbindxsize1)))
     >               int_mb(indx1_bondings(1)+nbbindxsize1) = 0
                  nbbindxsize1 = nbbindxsize1 + 2*n1

                  rtdb_name = 'pspw_2qmmm_frag_Kr0_bondings:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.btdb_get(rtdb,rtdb_name,mt_dbl,n1+2,
     >               dbl_mb(Kr0_bondings(1)+nbbKr0size)))
     >               dbl_mb(Kr0_bondings(1)+nbbKr0size) = 0.0d0
                  nbbKr0size = nbbKr0size + n1 + 2
               end do
            end do
         end if


*       ***** write out header info ****
        call Parallel_taskid(taskid)
        oprint= ((taskid.eq.MASTER).and.control_print(print_medium))
        if (oprint) then
           write(luout,*)
           write(luout,'(1x,A)') "Auxiliary Potential Parameters:"
           write(luout,'(2x,"fragments start (index_start) =",10I5)')
     >         (int_mb(indx_frag_start(1)+mm-1),mm=1,nfrag)
           write(luout,'(2x,"fragments kind                =",10I5)')
     >         (int_mb(kfrag(1)+mm-1),mm=1,nfrag)

           if (lmbda_flag)
     >        write(luout,'(2x,A,E19.10)')
     >        "auxilliary_lambda coupling parameter  = ",qmmm_lmbda

           do ia=1,nkfrag
             write(luout,'(/2x,A,I5)') "non-fragment: ",ia
             write(luout,'(2x,A,I5)') " - non-fragment size = ",
     >                         int_mb(size_frag(1)+ia-1)
             if (log_mb(self_interaction(1)+ia-1)) then
                write(luout,*) "  - ionic non-fragment"
             else
                write(luout,*) "  - covalent non-fragment"
             end if
             if (log_mb(incell_frag(1)+ia-1)) then
                write(luout,*) "  - incell non-fragment"
             else
                write(luout,*) "  - not incell non-fragment"
             end if
             write(luout,*) "  - Switching Parameters = ",
     >                  dbl_mb(switch_Rin(1)+ia-1),
     >                  dbl_mb(switch_Rout(1)+ia-1)

             nj = int_mb(nshake(1)+ia-1)
             i  = int_mb(nab_shake_start(1)+ia-1)
             ii = int_mb(nindx_shake_start(1)+ia-1)
             jj = int_mb(ndsq_shake_start(1)+ia-1)
             if (nj.gt.0) then
             do ja=1,nj
                natmp = int_mb(na(1)+i-1+ja-1)
                nbtmp = int_mb(nb(1)+i-1+ja-1)
                indx_ptr  = indx_shake(1)+ii-1
                shake_ptr = dsq_shake(1) +jj-1
                if (natmp.eq.nbtmp) then
                   write(luout,'(A,5(I4,F8.3))')   
     >                "   - shake = ",
     >                (int_mb(indx_ptr+j-1),
     >                 dsqrt(dbl_mb(shake_ptr+j-1)),
     >                 j=1,nbtmp)
             else
                   write(luout,'(A,5(I4,F8.3))')   
     >                "   - shake = ",
     >                (int_mb(indx_ptr+j-1),
     >                 dsqrt(dbl_mb(shake_ptr+j-1)),
     >                 j=1,nbtmp),
     >                int_mb(indx_ptr+natmp-1)
             end if
               ii = ii + natmp
               jj = jj + nbtmp
             end do
             end if

             ii = int_mb(nindx_bond_start(1)+ia-1)
             jj = int_mb(nKr0_bond_start(1) +ia-1)
             do ja=1,int_mb(nbond(1)+ia-1)
                indx_ptr  = indx_bond(1)+ii-1
                shake_ptr = Kr0_bond(1) +jj-1

             write(luout,'(A,2I4,4x,A,F11.6,4x,A,F11.6)') 
     >       "   - Bond Spring Parameters =  ",
     >                  (int_mb(indx_ptr+j-1),j=1,2),
     >                  'K=',
     >                  dbl_mb(shake_ptr),
     >                  're=',
     >                  dbl_mb(shake_ptr+1)
               ii = ii + 2
               jj = jj + 2
             end do

             ii = int_mb(nindx_angle_start(1)+ia-1)
             jj = int_mb(nKr0_angle_start(1) +ia-1)
             do ja=1,int_mb(nangle(1)+ia-1)
                indx_ptr  = indx_angle(1)+ii-1
                shake_ptr = Kr0_angle(1) +jj-1

             write(luout,'(A,3I4,4x,A,F11.6,4x,A,F11.6)') 
     >       "   - Angle Spring Parameters = ",
     >                  (int_mb(indx_ptr+j-1),j=1,3),
     >                  'Ktheta=',
     >                  dbl_mb(shake_ptr),
     >                  'thetae=',
     >                  dbl_mb(shake_ptr+1)
               ii = ii + 3
               jj = jj + 2
             end do

             ii = int_mb(nindx_cbond_start(1)+ia-1)
             jj = int_mb(nKr0_cbond_start(1) +ia-1)
             do ja=1,int_mb(ncbond(1)+ia-1)
                indx_ptr  = indx_cbond(1)+ii-1
                shake_ptr = Kr0_cbond(1) +jj-1

             write(luout,'(A,A,2I4,4x,A,F11.6, 
     >             /35x, A,2I4,4x,A,F11.6,4x,A,F11.6)') 
     >       "   - Cross Bond Spring Parameters: ",
     >                  'indexes1=',
     >                  (int_mb(indx_ptr+j-1),j=1,2),
     >                  'r1e=',
     >                  dbl_mb(shake_ptr),
     >                  'indexes2=',
     >                  (int_mb(indx_ptr+j-1),j=3,4),
     >                  'r2e=',
     >                  dbl_mb(shake_ptr+1),
     >                  'K12=',
     >                  dbl_mb(shake_ptr+2)
               ii = ii + 4
               jj = jj + 3
             end do

             ii = int_mb(nindx_mbond_start(1)+ia-1)
             jj = int_mb(nKr0_mbond_start(1) +ia-1)
             do ja=1,int_mb(nmbond(1)+ia-1)
                indx_ptr  = indx_mbond(1)+ii-1
                shake_ptr = Kr0_mbond(1) +jj-1

                kk = int_mb(indx_ptr+2)

                if (kk.eq.1) then
                        write(luout,'(3x,A49,E14.6,A5,E14.6,/5x,A,2I4)')
     >              '- Potential=4*e*((s/rij)**12-(s/rij)**6),      e:',
     >                    dbl_mb(shake_ptr),
     >                    ' s:',
     >                    dbl_mb(shake_ptr+1),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.2) then
                        write(luout,'(3x,A49,E14.6,A5,E14.6,A5,E14.6,
     >                                /5x,A,2I4)')
     >              '- Potential=A*exp(-rij/rho) - C/rij**6,        A:',
     >                    dbl_mb(shake_ptr),
     >                    ' rho:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' C:',
     >                    dbl_mb(shake_ptr+2),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.3) then
                write(luout,'(3x,A49,E14.6,A5,E14.6,A5,E14.6,A5,E14.6,
     >                        /5x,A,2I4)')
     >              '- Potential=A*exp(-rij/rho)-C/rij**6-D/rij**8, A:',
     >                    dbl_mb(shake_ptr),
     >                    ' rho:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' C:',
     >                    dbl_mb(shake_ptr+2),
     >                    ' D:',
     >                    dbl_mb(shake_ptr+3),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.4) then
                        write(luout,'(3x,A49,E14.6,A5,E14.6,A5,E14.6,
     >                                /5x,A,2I4)')
     >              '- Potential=A*(rij-B)**C,                      A:',
     >                    dbl_mb(shake_ptr),
     >                    ' B:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' C:',
     >                    dbl_mb(shake_ptr+2),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.5) then
                        write(luout,'(3x,A49,E14.6,A5,E14.6,/5x,A,2I4)')
     >              '- Potential=A*exp(-rij/rho),                   A:',
     >                    dbl_mb(shake_ptr),
     >                    ' rho:',
     >                    dbl_mb(shake_ptr+1),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.6) then
c             write(luout,'(A,A,2I4,3(4x,A,F11.6))')
c     >       "   - Morse Parameters: ",
c     >                  "indexes=",
c     >                  (int_mb(indx_ptr+j-1),j=1,2),
c     >                  "De=", dbl_mb(shake_ptr),
c     >                  "rho=",dbl_mb(shake_ptr+1),
c     >                  "re=", dbl_mb(shake_ptr+2)
                       write(luout,'(3x,A49,E14.6,A5,E14.6,A5,E14.6,
     >                              /5x,A,2I4)')
     >              '- Potential=De*(1-exp(-a*(rij-re)) )**2,      De:',
     >                    dbl_mb(shake_ptr),
     >                    ' a:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' re:',
     >                    dbl_mb(shake_ptr+2),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.7) then
                write(luout,'(3x,A,A,A,
     >                        /48x,A5,E14.6,A5,E14.6,A5,E14.6,A5,E14.6,
     >                        /48x,A5,E14.6,A5,E14.6,/5x,A,2I4)')
     >                    '- Potential=A*exp(-rij/rho)',
     >                    '-fdamp(6,Cd;rij)*C/rij**6',
     >                    '-fdamp(8,Dd;rij)*D/rij**8,',
     >                    ' A:',
     >                    dbl_mb(shake_ptr),
     >                    ' rho:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' C:',
     >                    dbl_mb(shake_ptr+2),
     >                    ' D:',
     >                    dbl_mb(shake_ptr+3),
     >                    ' Cd:',
     >                    dbl_mb(shake_ptr+4),
     >                    ' Dd:',
     >                    dbl_mb(shake_ptr+5),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.8) then
                      write(luout,'(3x,A49,E14.6,A5,E14.6,/5x,A,2I4)')
     >              '- Potential=A*exp(-(rij/rho)**2),              A:',
     >                    dbl_mb(shake_ptr),
     >                    ' rho:',
     >                    dbl_mb(shake_ptr+1) ,
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                else if (kk.eq.9) then
                        write(luout,'(3x,A49,E14.6,A5,E14.6,A5,E14.6,
     >                                /5x,A,2I4)')
     >              '- Potential=A*exp(-(rij/B)**C),                A:',
     >                    dbl_mb(shake_ptr),
     >                    ' B:',
     >                    dbl_mb(shake_ptr+1),
     >                    ' C:',
     >                    dbl_mb(shake_ptr+2),
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2)
                end if

               ii = ii + 3
               jj = jj + 6
             end do

             pi = 4.0*datan(1.0d0)
             ii = int_mb(nindx_dihedral_spring_start(1)+ia-1)
             jj = int_mb(nKr0_dihedral_spring_start(1) +ia-1)
             do ja=1,int_mb(ndihedral_spring(1)+ia-1)
                indx_ptr  = indx_dihedral_spring(1)+ii-1
                shake_ptr = Kr0_dihedral_spring(1) +jj-1

             write(luout,'(A,A,4I4,2(4x,A,F11.6))')
     >       "   - Dihedral Spring Parameters: ",
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,4),
     >                  "Kphi=",dbl_mb(shake_ptr),
     >                  "phie=",dbl_mb(shake_ptr+1)*180.0d0/pi
               ii = ii + 4
               jj = jj + 2
             end do

             pi = 4.0*datan(1.0d0)
             ii = int_mb(nindx_dihedral_start(1)+ia-1)
             jj = int_mb(nKr0_dihedral_start(1) +ia-1)
             do ja=1,int_mb(ndihedral(1)+ia-1)
                indx_ptr  = indx_dihedral(1)+ii-1+ja-1
                shake_ptr = Kr0_dihedral(1) +jj-1+ja-1

             write(luout,'(A,A,4I4,3(4x,A,F11.6))')
     >       "   - Dihedral Parameters: ",
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,4),
     >                  "Vphi=",dbl_mb(shake_ptr),
     >                  "nphi=",dbl_mb(shake_ptr+1),
     >                  "phie=",dbl_mb(shake_ptr+2)*180.0d0/pi
               ii = ii + 4
               jj = jj + 3
             end do

             ll = int_mb(nindx0_coord_start(1)+ia-1)
             ii = int_mb(nindx1_coord_start(1)+ia-1)
             jj = int_mb(nindx2_coord_start(1)+ia-1)
             kk = int_mb(nKr0_coord_start(1)  +ia-1)
             do ja=1,int_mb(ncoord(1)+ia-1)
                i  = indx1_coord(1)+ii-1
                j  = indx2_coord(1)+jj-1
                ni = int_mb(nsize1_coord(1)+ll-1)
                nj = int_mb(nsize2_coord(1)+ll-1)
                shake_ptr = Kr0_coord(1)+kk-1

                if (dbl_mb(shake_ptr+5).gt.0) then
                write(luout,'(A,/6x,A,F11.6,4x,A,F11.6,
     >                      /6x,A,F11.6,4x,A,F11.6)')
     >       "   - Coordination Number Parameters (Sprik): ",
     >                  "n=    ",dbl_mb(shake_ptr),
     >                  "r0=   ",dbl_mb(shake_ptr+2),
     >                  "K=    ",dbl_mb(shake_ptr+3),
     >                  "cn0=  ",dbl_mb(shake_ptr+4)
               else
                write(luout,'(A,/6x,A,F11.6,4x,A,F11.6,4x,A,F11.6,
     >                      /6x,A,F11.6,4x,A,F11.6)')
     >       "   - Coordination Number Parameters (LJ): ",
     >                  "n=    ",dbl_mb(shake_ptr),
     >                  "m=    ",dbl_mb(shake_ptr+1),
     >                  "r0=   ",dbl_mb(shake_ptr+2),
     >                  "K=    ",dbl_mb(shake_ptr+3),
     >                  "cn0=  ",dbl_mb(shake_ptr+4)
               end if
               write(luout,'(6x,"Index1 =",10I5)')
     >         (int_mb(i+mm-1),mm=1,ni)
               write(luout,'(6x,"Index2 =",10I5)')
     >         (int_mb(j+mm-1),mm=1,nj)

               ll = ll + 1
               ii = ii + ni
               jj = jj + nj
               kk = kk + 6
             end do

             ll = int_mb(nindx0_bondings_start(1)+ia-1)
             ii = int_mb(nindx1_bondings_start(1)+ia-1)
             kk = int_mb(nKr0_bondings_start(1)  +ia-1)
             do ja=1,int_mb(nbondings(1)+ia-1)
                i  = indx1_bondings(1)+ii-1
                ni = int_mb(nsize1_bondings(1)+ll-1)
                shake_ptr = Kr0_bondings(1)+kk-1

                write(luout,'(A,/6x,A,F11.6,4x,A,F11.6,4x,A,I5)')
     >                "   - Bondings Parameters: ",
     >                "-- K0=",dbl_mb(shake_ptr),
     >                "gamma0=",dbl_mb(shake_ptr+1),
     >                "nb=",ni
                write(luout,'(6x,A)') "-- coefficient index1 index2" 
                do mm=1,ni
                   write(luout,'(6x,F14.6,I7,I7)') 
     >                dbl_mb(shake_ptr+1+mm),
     >                int_mb(i+2*mm-2),
     >                int_mb(i+2*mm-1)
                end do
                ii = ii + 2*ni
                kk = kk + ni+2
                ll = ll + 1
             end do


             
           end do
 
        end if


      end if

      return
      end


*     **********************************
*     *	                               *
*     *         pspw_2qmmm_end 	       *
*     *                                *
*     **********************************

      subroutine pspw_2qmmm_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_2qmmm.fh"

      logical value
      integer  control_version
      external control_version

      if (qmmm_found) then
        value =           BA_free_heap(fion2tmp(2))
        value = value.and.BA_free_heap(indx_frag_start(2))
        value = value.and.BA_free_heap(size_frag(2))
        value = value.and.BA_free_heap(kfrag(2))
        value = value.and.BA_free_heap(self_interaction(2))
        value = value.and.BA_free_heap(incell_frag(2))
        value = value.and.BA_free_heap(switch_Rin(2))
        value = value.and.BA_free_heap(switch_Rout(2))
        value = value.and.BA_free_heap(nshake(2))
        value = value.and.BA_free_heap(nab_shake_start(2))
        value = value.and.BA_free_heap(nindx_shake_start(2))
        value = value.and.BA_free_heap(ndsq_shake_start(2))
        if (shake_found) then
          value = value.and.BA_free_heap(na(2))
          value = value.and.BA_free_heap(nb(2))
          value = value.and.BA_free_heap(indx_shake(2))
          value = value.and.BA_free_heap(dsq_shake(2))
        end if

        value = value.and.BA_free_heap(nbond(2))
        value = value.and.BA_free_heap(nindx_bond_start(2))
        value = value.and.BA_free_heap(nKr0_bond_start(2))
        if (bond_found) then
          value = value.and.BA_free_heap(indx_bond(2))
          value = value.and.BA_free_heap(Kr0_bond(2))
        end if

        value = value.and.BA_free_heap(nangle(2))
        value = value.and.BA_free_heap(nindx_angle_start(2))
        value = value.and.BA_free_heap(nKr0_angle_start(2))
        if (angle_found) then
          value = value.and.BA_free_heap(indx_angle(2))
          value = value.and.BA_free_heap(Kr0_angle(2))
        end if

        value = value.and.BA_free_heap(ncbond(2))
        value = value.and.BA_free_heap(nindx_cbond_start(2))
        value = value.and.BA_free_heap(nKr0_cbond_start(2))
        if (cbond_found) then
          value = value.and.BA_free_heap(indx_cbond(2))
          value = value.and.BA_free_heap(Kr0_cbond(2))
        end if

        value = value.and.BA_free_heap(nmbond(2))
        value = value.and.BA_free_heap(nindx_mbond_start(2))
        value = value.and.BA_free_heap(nKr0_mbond_start(2))
        if (mbond_found) then
          value = value.and.BA_free_heap(indx_mbond(2))
          value = value.and.BA_free_heap(Kr0_mbond(2))
        end if

        value = value.and.BA_free_heap(ndihedral_spring(2))
        value = value.and.BA_free_heap(nindx_dihedral_spring_start(2))
        value = value.and.BA_free_heap(nKr0_dihedral_spring_start(2))
        if (dihedral_spring_found) then
          value = value.and.BA_free_heap(indx_dihedral_spring(2))
          value = value.and.BA_free_heap(Kr0_dihedral_spring(2))
        end if

        value = value.and.BA_free_heap(ndihedral(2))
        value = value.and.BA_free_heap(nindx_dihedral_start(2))
        value = value.and.BA_free_heap(nKr0_dihedral_start(2))
        if (dihedral_found) then
          value = value.and.BA_free_heap(indx_dihedral(2))
          value = value.and.BA_free_heap(Kr0_dihedral(2))
        end if

        value = value.and.BA_free_heap(ncoord(2))
        value = value.and.BA_free_heap(nindx0_coord_start(2))
        value = value.and.BA_free_heap(nindx1_coord_start(2))
        value = value.and.BA_free_heap(nindx2_coord_start(2))
        value = value.and.BA_free_heap(nKr0_coord_start(2))
        if (coord_found) then
          value = value.and.BA_free_heap(nsize1_coord(2))
          value = value.and.BA_free_heap(nsize2_coord(2))
          value = value.and.BA_free_heap(indx1_coord(2))
          value = value.and.BA_free_heap(indx2_coord(2))
          value = value.and.BA_free_heap(Kr0_coord(2))
        end if

        value = value.and.BA_free_heap(nbondings(2))
        value = value.and.BA_free_heap(nindx0_bondings_start(2))
        value = value.and.BA_free_heap(nindx1_bondings_start(2))
        value = value.and.BA_free_heap(nKr0_bondings_start(2))
        if (bondings_found) then
          value = value.and.BA_free_heap(nsize1_bondings(2))
          value = value.and.BA_free_heap(indx1_bondings(2))
          value = value.and.BA_free_heap(Kr0_bondings(2))
        end if

        if (.not.value) call errquit('cannot free heap memory',0,
     &       MA_ERR)

      end if

      return
      end

*     **********************************
*     *	                               *
*     *         pspw_2qmmm_found        *
*     *                                *
*     **********************************
      logical function pspw_2qmmm_found()
      implicit none

#include "pspw_2qmmm.fh"

      pspw_2qmmm_found = qmmm_found
      return
      end


*     **********************************
*     *                                *
*     *    pspw_2qmmm_shake_found      *
*     *                                *
*     **********************************
      logical function pspw_2qmmm_shake_found()
      implicit none

#include "pspw_2qmmm.fh"

      pspw_2qmmm_shake_found = shake_found
      return
      end



*     **********************************
*     *                                *
*     *        pspw_2qmmm_spring_E     *
*     *                                *
*     **********************************
      real*8 function pspw_2qmmm_spring_E()
      implicit none

#include "pspw_2qmmm.fh"

      !**** external functions ****
      real*8   pspw_2qmmm_spring_E_sub
      external pspw_2qmmm_spring_E_sub

      pspw_2qmmm_spring_E = pspw_2qmmm_spring_E_sub(qmmm_lmbda)
      return
      end 

*     **********************************
*     *                                *
*     *        pspw_2qmmm_Emix         *
*     *                                *
*     **********************************
      real*8 function pspw_2qmmm_Emix()
      implicit none

      !**** external functions ****
      real*8   pspw_2qmmm_spring_E_sub
      external pspw_2qmmm_spring_E_sub

      pspw_2qmmm_Emix = pspw_2qmmm_spring_E_sub(1.0d0)
      return
      end 


*     **********************************
*     *                                *
*     *        pspw_2qmmm_spring_Esub  *
*     *                                *
*     **********************************
      real*8 function pspw_2qmmm_spring_E_sub(lmbda)
      implicit none
      real*8 lmbda

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_2qmmm.fh"

c     **** local variables ****
      integer dutask,taskid,np
      integer w1,wk1,ks1,nbs,nas,ncbs,nbm,ii,jj,kk,ll,nds,ndds,ncns,nbbs
      real*8  E

c     **** external functions ***
      integer  ion_rion_ptr
      external ion_rion_ptr
      real*8   pspw_qmmm_spring_bond_frag
      external pspw_qmmm_spring_bond_frag
      real*8   pspw_qmmm_spring_angle_frag
      external pspw_qmmm_spring_angle_frag
      real*8   pspw_qmmm_spring_cbond_frag
      external pspw_qmmm_spring_cbond_frag
      real*8   pspw_qmmm_morse_bond_frag
      external pspw_qmmm_morse_bond_frag
      real*8   pspw_qmmm_dihedral_spring_frag
      external pspw_qmmm_dihedral_spring_frag
      real*8   pspw_qmmm_dihedral_frag
      external pspw_qmmm_dihedral_frag
      real*8   pspw_qmmm_spring_coordnum_frag
      external pspw_qmmm_spring_coordnum_frag
      real*8   pspw_qmmm_spring_bondings_frag
      external pspw_qmmm_spring_bondings_frag


      call Parallel_np(np)
      call Parallel_taskid(taskid)
      dutask = 0

      E = 0.0d0
      do w1=1,nfrag
         if (dutask.eq.taskid) then
         ks1  = int_mb(indx_frag_start(1)+w1-1)
         wk1  = int_mb(kfrag(1)+w1-1)
         nbs  = int_mb(nbond(1)+wk1-1)
         nas  = int_mb(nangle(1)+wk1-1)
         ncbs = int_mb(ncbond(1)+wk1-1)
         nbm  = int_mb(nmbond(1)+wk1-1)
         nds   = int_mb(ndihedral_spring(1)+wk1-1)
         ndds  = int_mb(ndihedral(1)+wk1-1)
         ncns  = int_mb(ncoord(1)+wk1-1)
         nbbs  = int_mb(nbondings(1)+wk1-1)
         if (nbs.gt.0) then
            jj  = int_mb(nindx_bond_start(1)+wk1-1)
            kk  = int_mb(nKr0_bond_start(1) +wk1-1)
            E = E + pspw_qmmm_spring_bond_frag(nbs,
     >                                int_mb(indx_bond(1)+jj-1),
     >                                dbl_mb(Kr0_bond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (nas.gt.0) then
            jj  = int_mb(nindx_angle_start(1)+wk1-1)
            kk  = int_mb(nKr0_angle_start(1) +wk1-1)
            E = E + pspw_qmmm_spring_angle_frag(nas,
     >                                int_mb(indx_angle(1)+jj-1),
     >                                dbl_mb(Kr0_angle(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (ncbs.gt.0) then
            jj  = int_mb(nindx_cbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_cbond_start(1) +wk1-1)
            E = E + pspw_qmmm_spring_cbond_frag(ncbs,
     >                                int_mb(indx_cbond(1)+jj-1),
     >                                dbl_mb(Kr0_cbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (nbm.gt.0) then
            jj  = int_mb(nindx_mbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_mbond_start(1) +wk1-1)
            E = E + pspw_qmmm_morse_bond_frag(nbm,
     >                                int_mb(indx_mbond(1)+jj-1),
     >                                dbl_mb(Kr0_mbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (nds.gt.0) then
            jj  = int_mb(nindx_dihedral_spring_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_spring_start(1) +wk1-1)
            E = E + pspw_qmmm_dihedral_spring_frag(nds,
     >                        int_mb(indx_dihedral_spring(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral_spring(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (ndds.gt.0) then
            jj  = int_mb(nindx_dihedral_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_start(1) +wk1-1)
            E = E + pspw_qmmm_dihedral_frag(ndds,
     >                        int_mb(indx_dihedral(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (ncns.gt.0) then
            ii  = int_mb(nindx0_coord_start(1)+wk1-1)
            jj  = int_mb(nindx1_coord_start(1)+wk1-1)
            kk  = int_mb(nindx2_coord_start(1)+wk1-1)
            ll  = int_mb(nKr0_coord_start(1)  +wk1-1)
            E = E + pspw_qmmm_spring_coordnum_frag(ncns,
     >                        int_mb(nsize1_coord(1)+ii-1),
     >                        int_mb(indx1_coord(1) +jj-1),
     >                        int_mb(nsize2_coord(1)+ii-1),
     >                        int_mb(indx2_coord(1) +kk-1),
     >                        dbl_mb(Kr0_coord(1) +ll-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if 
         if (nbbs.gt.0) then
            ii  = int_mb(nindx0_bondings_start(1)+wk1-1)
            jj  = int_mb(nindx1_bondings_start(1)+wk1-1)
            kk  = int_mb(nKr0_bondings_start(1) +wk1-1)
            E = E + pspw_qmmm_spring_bondings_frag(nbbs,
     >                        int_mb(nsize1_bondings(1)+ii-1),
     >                        int_mb(indx1_bondings(1)+jj-1),
     >                        dbl_mb(Kr0_bondings(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if

         end if
         dutask = mod(dutask+1,np)
        
      end do
      if (np.gt.1) call Parallel_SumAll(E)

      pspw_2qmmm_spring_E_sub = E*lmbda
      return
      end


*     **********************************
*     *                                *
*     *      pspw_2qmmm_spring_fion    *
*     *                                *
*     **********************************

      subroutine pspw_2qmmm_spring_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_2qmmm.fh"

c     **** local variables ****
      integer dutask,taskid,np
      integer w1,wk1,ks1,nbs,nas,ncbs,nbm,ii,jj,kk,ll,nds,ndds,ncns,nbbs
      integer nion3

c     **** external functions ***
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      dutask = 0

      nion3 = 3*ion_nion()
      call Parallel_shared_vector_zero(.true.,nion3,dbl_mb(fion2tmp(1)))
!$OMP DO
      do w1=1,nfrag
         if (mod(w1-1,np).eq.taskid) then
         ks1  = int_mb(indx_frag_start(1)+w1-1)
         wk1  = int_mb(kfrag(1)+w1-1)
         nbs  = int_mb(nbond(1)+wk1-1)
         nas  = int_mb(nangle(1)+wk1-1)
         ncbs = int_mb(ncbond(1)+wk1-1)
         nbm  = int_mb(nmbond(1)+wk1-1)
         nds  = int_mb(ndihedral_spring(1)+wk1-1)
         ndds  = int_mb(ndihedral(1)+wk1-1)
         ncns  = int_mb(ncoord(1)+wk1-1)
         nbbs  = int_mb(nbondings(1)+wk1-1)
         if (nbs.gt.0) then
            jj  = int_mb(nindx_bond_start(1)+wk1-1)
            kk  = int_mb(nKr0_bond_start(1) +wk1-1)
            call pspw_qmmm_spring_bond_frag_fion(nbs,
     >                                int_mb(indx_bond(1)+jj-1),
     >                                dbl_mb(Kr0_bond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (nas.gt.0) then
            jj  = int_mb(nindx_angle_start(1)+wk1-1)
            kk  = int_mb(nKr0_angle_start(1) +wk1-1)
            call pspw_qmmm_spring_angle_frag_fion(nas,
     >                                int_mb(indx_angle(1)+jj-1),
     >                                dbl_mb(Kr0_angle(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (ncbs.gt.0) then
            jj  = int_mb(nindx_cbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_cbond_start(1) +wk1-1)
            call pspw_qmmm_spring_cbond_frag_fion(ncbs,
     >                                int_mb(indx_cbond(1)+jj-1),
     >                                dbl_mb(Kr0_cbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (nbm.gt.0) then
            jj  = int_mb(nindx_mbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_mbond_start(1) +wk1-1)
            call pspw_qmmm_morse_bond_frag_fion(nbm,
     >                                int_mb(indx_mbond(1)+jj-1),
     >                                dbl_mb(Kr0_mbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (nds.gt.0) then
            jj  = int_mb(nindx_dihedral_spring_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_spring_start(1) +wk1-1)
            call pspw_qmmm_dihedral_spring_frag_fion(nds,
     >                        int_mb(indx_dihedral_spring(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral_spring(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (ndds.gt.0) then
            jj  = int_mb(nindx_dihedral_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_start(1) +wk1-1)
            call pspw_qmmm_dihedral_frag_fion(ndds,
     >                        int_mb(indx_dihedral(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (ncns.gt.0) then
            ii  = int_mb(nindx0_coord_start(1)+wk1-1)
            jj  = int_mb(nindx1_coord_start(1)+wk1-1)
            kk  = int_mb(nindx2_coord_start(1)+wk1-1)
            ll  = int_mb(nKr0_coord_start(1) +wk1-1)
            call pspw_qmmm_spring_coordnum_frag_fion(ncns,
     >                        int_mb(nsize1_coord(1)+ii-1),
     >                        int_mb(indx1_coord(1) +jj-1),
     >                        int_mb(nsize2_coord(1)+ii-1),
     >                        int_mb(indx2_coord(1) +kk-1),
     >                        dbl_mb(Kr0_coord(1)   +ll-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if
         if (nbbs.gt.0) then
            ii  = int_mb(nindx0_bondings_start(1)+wk1-1)
            jj  = int_mb(nindx1_bondings_start(1)+wk1-1)
            kk  = int_mb(nKr0_bondings_start(1) +wk1-1)
            call pspw_qmmm_spring_bondings_frag_fion(nbbs,
     >                        int_mb(nsize1_bondings(1)+ii-1),
     >                        int_mb(indx1_bondings(1)+jj-1),
     >                        dbl_mb(Kr0_bondings(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        dbl_mb(fion2tmp(1)+3*(ks1-1)))
         end if

         end if
         !dutask = mod(dutask+1,np)
      end do
!$OMP END DO
      call DAXPY_OMP(nion3,qmmm_lmbda,dbl_mb(fion2tmp(1)),1,fion,1)

      return
      end




*     ******************************
*     *                            *
*     *       pspw_2qmmm_shake     *
*     *                            *
*     ******************************
*   shakes the 2fragments 
*
      subroutine pspw_2qmmm_shake(r2,r1)
      implicit none
      real*8 r2(3,*),r1(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_2qmmm.fh"

*     **** local variables ****
      integer maxit
      double precision tol
      parameter (maxit=1000,tol=1.0d-4)

      integer w1,wk1,ks1
      integer nshk,ii,jj,kk

*     **** external functions ****
      integer  ion_amass_ptr
      external ion_amass_ptr

!$OMP DO
      do w1=1,nfrag
         wk1  = int_mb(kfrag(1)+w1-1)
         nshk = int_mb(nshake(1)+wk1-1)
         if (nshk.gt.0) then
            ks1 = int_mb(indx_frag_start(1)+w1-1)
            ii  = int_mb(nab_shake_start(1)  +wk1-1)
            jj  = int_mb(nindx_shake_start(1)+wk1-1)
            kk  = int_mb(ndsq_shake_start(1) +wk1-1)
            call pspw_qmmm_shake_frag(nshk,
     >                                int_mb(na(1)+ii-1),
     >                                int_mb(nb(1)+ii-1),
     >                                int_mb(indx_shake(1)+jj-1),
     >                                dbl_mb(dsq_shake(1) +kk-1),
     >                                tol,maxit,
     >                                dbl_mb(ion_amass_ptr()+ks1-1),
     >                                r2(1,ks1),r1(1,ks1))
     >
         end if
      end do
!$OMP END DO
  

      return
      end

*     ***************************************
*     *                                     *
*     *      pspw_2qmmm_ion_add_constraint  *
*     *                                     *
*     ***************************************
      subroutine pspw_2qmmm_ion_add_constraint()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "pspw_2qmmm.fh"

*     **** local variables ****
      integer w1,wk1,ks1,ns
      integer nshk,ii,jj,kk

      do w1=1,nfrag
         wk1  = int_mb(kfrag(1)+w1-1)
         nshk = int_mb(nshake(1)+wk1-1)
         if (nshk.gt.0) then
            ks1 = int_mb(indx_frag_start(1)+w1-1)
            ii  = int_mb(nab_shake_start(1)  +wk1-1)
            jj  = int_mb(nindx_shake_start(1)+wk1-1)
            kk  = int_mb(ndsq_shake_start(1) +wk1-1)
            do ns=1,nshk
               call ion_add_constraint(int_mb(nb(1)+ii-1+ns-1))
            end do
         end if
      end do
      return
      end



*     **********************************
*     *                                *
*     *        pspw_2qmmm_lambda       *
*     *                                *
*     **********************************
      real*8 function pspw_2qmmm_lambda()
      implicit none

#include "pspw_qmmm.fh"

      pspw_2qmmm_lambda = qmmm_lmbda
      return
      end

*     **********************************
*     *                                *
*     *     pspw_2qmmm_lambda_flag     *
*     *                                *
*     **********************************
      logical function pspw_2qmmm_lambda_flag()
      implicit none

#include "pspw_qmmm.fh"

      pspw_2qmmm_lambda_flag = lmbda_flag
      return
      end


