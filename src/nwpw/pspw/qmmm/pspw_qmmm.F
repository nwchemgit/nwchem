c
c $Id: pspw_qmmm.F,v 1.35 2009-03-26 00:19:01 bylaska Exp $
c

***********************************************************
*                                                         *
*                    PSPW-QMMM  module                    *
*                                                         *
* Top-level interface to Model_potential and Model_charge *
* modules                                                 *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*    -- developed by Eric J. Bylaska on October 18,2001   *
*                                                         *
***********************************************************
*
*
*  

*     **********************************
*     *	                               *
*     *          pspw_qmmm_init        *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer i,j,ii,jj,ia,ja,ni,nj,frag_size,shake_ptr,indx_ptr
      real*8  volume
      character*80 rtdb_name
      integer natmp,nbtmp,nabsize,nindxsize,ndsqsize
      integer nbindxsize,nKr0size
      integer naindxsize,nQr0size
      integer ncbindxsize,ncKr0size
      integer nmbindxsize,nmKr0size
      integer ndindxsize,ndKr0size
      integer nddindxsize,nddKr0size

*     **** external functions ****
      integer     ion_nion,ion_nion_qm,ion_nion_mm,control_version
      integer     ewald_nshl3d,ewald_rcell_ptr,ion_rion_ptr
      character*7 c_index_name
      real*8      lattice_unita
      external    ion_nion,ion_nion_qm,ion_nion_mm,control_version
      external    ewald_nshl3d,ewald_rcell_ptr,ion_rion_ptr
      external    c_index_name
      external    lattice_unita

      rtdb_name = 'pspw_qmmm_auxon'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_log,1,value)) value = .false.

      qmmm_found     = (ion_nion_mm().gt.0).or. value 
      auxiliary_only = (ion_nion_mm().le.0).and.value
      shake_found    = .false.
      nfrag  = 0
      nkfrag = 0
      

*     **** don't do anything if there are no MM atoms ****
      if (qmmm_found) then

*        **** set the nshl3d and rcell ****
         if (.not.auxiliary_only) then
         if (control_version().eq.3) then
            nshl3d   = ewald_nshl3d()
            rcell(1) = ewald_rcell_ptr()
            periodic = .true. 
         else
            periodic = .false.
            nshl3d = 1
            if (.not.MA_alloc_get(mt_dbl,3,'rcell',rcell(2),rcell(1))) 
     >      call errquit(
     >      'pspw_qmmm_init:failed allocating rcell from heap',0,MA_ERR)
            dbl_mb(rcell(1))   = 0.0d0
            dbl_mb(rcell(1)+1) = 0.0d0
            dbl_mb(rcell(1)+2) = 0.0d0
         end if
         end if

*        **** reading number of frags ****
         rtdb_name = 'pspw_qmmm_nfrag'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,nfrag)) nfrag = 0

*        **** reading number of kinds of frags ****
         rtdb_name = 'pspw_qmmm_nkfrag'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,nkfrag)) nkfrag = 0

         if (nfrag.eq.0) go to 999

         value = MA_alloc_get(mt_int,nfrag,'index_frag_start',
     >                        indx_frag_start(2),indx_frag_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nfrag,'size_frag',
     >                        size_frag(2),size_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nfrag,'kfrag',kfrag(2),kfrag(1))

         value = value.and.
     >           MA_alloc_get(mt_log,nkfrag,'self_interaction',
     >                        self_interaction(2),self_interaction(1))
         value = value.and.
     >           MA_alloc_get(mt_log,nkfrag,'incell_frag',
     >                        incell_frag(2),incell_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,nkfrag,'switch_Rin',
     >                        switch_Rin(2),switch_Rin(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,nkfrag,'switch_Rout',
     >                        switch_Rout(2),switch_Rout(1))

*        **** shake variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nshake',
     >                        nshake(2),nshake(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nab_shake_start',
     >                        nab_shake_start(2),nab_shake_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nindx_shake_start',
     >                        nindx_shake_start(2),nindx_shake_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'ndsq_shake_start',
     >                        ndsq_shake_start(2),ndsq_shake_start(1))

*        **** bond variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nbond',
     >                        nbond(2),nbond(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nindx_bond_start',
     >                        nindx_bond_start(2),nindx_bond_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nKr0_bond_start',
     >                        nKr0_bond_start(2),nKr0_bond_start(1))

*        **** angle variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nangle',
     >                        nangle(2),nangle(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nindx_angle_start',
     >                        nindx_angle_start(2),nindx_angle_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nKr0_angle_start',
     >                        nKr0_angle_start(2),nKr0_angle_start(1))

*        **** cbond variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'ncbond',
     >                        ncbond(2),ncbond(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nindx_cbond_start',
     >                        nindx_cbond_start(2),nindx_cbond_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nKr0_cbond_start',
     >                        nKr0_cbond_start(2),nKr0_cbond_start(1))

*        **** mbond variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nmbond',
     >                        nmbond(2),nmbond(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nindx_mbond_start',
     >                        nindx_mbond_start(2),nindx_mbond_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nKr0_mbond_start',
     >                        nKr0_mbond_start(2),nKr0_mbond_start(1))

*        **** dihedral spring variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'ndihedral_spring',
     >                        ndihedral_spring(2),ndihedral_spring(1))
         value = value.and.
     >         MA_alloc_get(mt_int,nkfrag,'nindx_dihedral_spring_start',
     >                        nindx_dihedral_spring_start(2),
     >                        nindx_dihedral_spring_start(1))
         value = value.and.
     >         MA_alloc_get(mt_int,nkfrag,'nKr0_dihedral_spring_start',
     >                      nKr0_dihedral_spring_start(2),
     >                      nKr0_dihedral_spring_start(1))

*        **** dihedral variables ****
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'ndihedral',
     >                        ndihedral(2),ndihedral(1))
         value = value.and.
     >         MA_alloc_get(mt_int,nkfrag,'nindx_dihedral_start',
     >                        nindx_dihedral_start(2),
     >                        nindx_dihedral_start(1))
         value = value.and.
     >         MA_alloc_get(mt_int,nkfrag,'nKr0_dihedral_start',
     >                      nKr0_dihedral_start(2),
     >                      nKr0_dihedral_start(1))


         if (.not.value) 
     >    call errquit('cannot allocate heap memory for qmmm',0,
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

         do ia=1,nkfrag

            rtdb_name = 'pspw_qmmm_frag_size:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,frag_size)
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading frag_size',0,RTDB_ERR)

            rtdb_name = 'pspw_qmmm_frag_nindex_start:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,ni)
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading nindx_start',0,RTDB_ERR)
            
            do i=1,ni
               int_mb(size_frag(1)+ii+i-1) = frag_size
               int_mb(kfrag(1)+ii+i-1)     = ia
            end do

            rtdb_name = 'pspw_qmmm_frag_index_start:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,
     >                       ni,int_mb(indx_frag_start(1)+ii))
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading indx_start',0,RTDB_ERR)
            ii = ii + ni

            rtdb_name = 'pspw_qmmm_frag_self_interaction:'
     >                  //c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_log,1,
     >                       log_mb(self_interaction(1)+ia-1)))
     >          log_mb(self_interaction(1)+ia-1) = .false.

            rtdb_name = 'pspw_qmmm_frag_incell:'
     >                  //c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_log,1,
     >                       log_mb(incell_frag(1)+ia-1)))
     >          log_mb(incell_frag(1)+ia-1) = .true.

            rtdb_name = 'pspw_qmmm_frag_switch_Rin:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rin(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rin(1)+ia-1) =(2.0160d0/0.529177d0)

            rtdb_name = 'pspw_qmmm_frag_switch_Rout:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rout(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rout(1)+ia-1) =(3.1287d0/0.529177d0)


*           **** shake stuff  - get nshake  ****
            int_mb(nab_shake_start(1)+ia-1) = nabsize + 1
            rtdb_name = 'pspw_qmmm_frag_nshk:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nshake(1)+ia-1))) 
     >         int_mb(nshake(1)+ia-1) = 0
            nabsize = nabsize + int_mb(nshake(1)+ia-1)


*           **** shake stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_shake_start(1)+ia-1) = nindxsize + 1
            int_mb(ndsq_shake_start(1) +ia-1) = ndsqsize  + 1

            do ja=1,int_mb(nshake(1)+ia-1)
               rtdb_name = 'pspw_qmmm_frag_na:'
     >                     //c_index_name(ia)//c_index_name(ja)
               if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,natmp))
     >            natmp = 0
               nindxsize = nindxsize + natmp

               rtdb_name = 'pspw_qmmm_frag_nb:'
     >                     //c_index_name(ia)//c_index_name(ja)
               if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,nbtmp))
     >            nbtmp = 0
               ndsqsize = ndsqsize + nbtmp
            end do

*           **** bond stuff  - get nbond  ****
            rtdb_name = 'pspw_qmmm_frag_nbond:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nbond(1)+ia-1))) 
     >         int_mb(nbond(1)+ia-1) = 0

*           **** bond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_bond_start(1)+ia-1) = nbindxsize + 1
            int_mb(nKr0_bond_start(1) +ia-1) = nKr0size  + 1
            nbindxsize = nbindxsize + 2*int_mb(nbond(1)+ia-1)
            nKr0size   = nKr0size   + 2*int_mb(nbond(1)+ia-1)


*           **** angle stuff  - get nangle  ****
            rtdb_name = 'pspw_qmmm_frag_nangle:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nangle(1)+ia-1)))
     >         int_mb(nangle(1)+ia-1) = 0

*           **** angle stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_angle_start(1)+ia-1) = naindxsize + 1
            int_mb(nKr0_angle_start(1) +ia-1) = nQr0size  + 1
            naindxsize = naindxsize + 3*int_mb(nangle(1)+ia-1)
            nQr0size   = nQr0size   + 2*int_mb(nangle(1)+ia-1)

*           **** cbond stuff  - get ncbond  ****
            rtdb_name = 'pspw_qmmm_frag_ncbond:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ncbond(1)+ia-1))) 
     >         int_mb(ncbond(1)+ia-1) = 0

*           **** cbond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_cbond_start(1)+ia-1) = ncbindxsize + 1
            int_mb(nKr0_cbond_start(1) +ia-1) = ncKr0size  + 1
            ncbindxsize = ncbindxsize + 4*int_mb(ncbond(1)+ia-1)
            ncKr0size   = ncKr0size   + 3*int_mb(ncbond(1)+ia-1)

*           **** mbond stuff  - get nmbond  ****
            rtdb_name = 'pspw_qmmm_frag_nmbond:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(nmbond(1)+ia-1))) 
     >         int_mb(nmbond(1)+ia-1) = 0

*           **** mbond stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_mbond_start(1)+ia-1) = nmbindxsize + 1
            int_mb(nKr0_mbond_start(1) +ia-1) = nmKr0size  + 1
            nmbindxsize = nmbindxsize + 2*int_mb(nmbond(1)+ia-1)
            nmKr0size   = nmKr0size   + 3*int_mb(nmbond(1)+ia-1)


*           **** dihedral_spring stuff  - get ndihedral_spring  ****
            rtdb_name = 
     >        'pspw_qmmm_frag_ndihedral_spring:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ndihedral_spring(1)+ia-1)))
     >         int_mb(ndihedral_spring(1)+ia-1) = 0

*           **** dihedral_spring stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_dihedral_spring_start(1)+ia-1) = ndindxsize + 1
            int_mb(nKr0_dihedral_spring_start(1) +ia-1) = ndKr0size  + 1
            ndindxsize = ndindxsize + 4*int_mb(ndihedral_spring(1)+ia-1)
            ndKr0size  = ndKr0size  + 2*int_mb(ndihedral_spring(1)+ia-1)


*           **** dihedral stuff  - get ndihedral  ****
            rtdb_name =  
     >        'pspw_qmmm_frag_ndihedral:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                        int_mb(ndihedral(1)+ia-1)))
     >         int_mb(ndihedral(1)+ia-1) = 0

*           **** dihedral stuff  - determine nindxsize and ndsqsize  ****
            int_mb(nindx_dihedral_start(1)+ia-1) = nddindxsize + 1
            int_mb(nKr0_dihedral_start(1) +ia-1) = nddKr0size  + 1
            nddindxsize = nddindxsize + 4*int_mb(ndihedral(1)+ia-1)
            nddKr0size  = nddKr0size  + 3*int_mb(ndihedral(1)+ia-1)

         end do


*        *** more shake ****
         shake_found = (nabsize.gt.0)
         if (shake_found) then

            value = MA_alloc_get(mt_int,nabsize,'na',na(2),na(1))
            value = value.and.
     >              MA_alloc_get(mt_int,nabsize,'nb',nb(2),nb(1))
            value = value.and.
     >              MA_alloc_get(mt_int,nindxsize,'indx_shake',
     >                        indx_shake(2),indx_shake(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,ndsqsize,'dsq_shake',
     >                        dsq_shake(2),dsq_shake(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm shake',0,MA_ERR)

            nabsize   = 0
            nindxsize = 0
            ndsqsize  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nshake(1)+ia-1)
                  rtdb_name = 'pspw_qmmm_frag_na:'
     >                        //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,natmp))
     >               natmp = 0
                  int_mb(na(1)+nabsize) = natmp
                  rtdb_name = 'pspw_qmmm_frag_nb:'
     >                        //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,nbtmp))
     >               nbtmp = 0
                  int_mb(nb(1)+nabsize) = nbtmp
                  nabsize = nabsize + 1

                  rtdb_name = 'pspw_qmmm_frag_indx_shake:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,natmp,
     >                int_mb(indx_shake(1)+nindxsize)))
     >               int_mb(indx_shake(1)+nindxsize) = 0
                  nindxsize = nindxsize + natmp

                  rtdb_name = 'pspw_qmmm_frag_dsq_shake:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,nbtmp,
     >                dbl_mb(dsq_shake(1)+ndsqsize)))
     >               dbl_mb(dsq_shake(1)+ndsqsize) = 0.0d0
                  ndsqsize = ndsqsize + nbtmp

               end do
            end do
         end if

*        *** more bond ****
         bond_found = (nbindxsize.gt.0)
         if (bond_found) then
            value = MA_alloc_get(mt_int,nbindxsize,'indx_bond',
     >                        indx_bond(2),indx_bond(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,nKr0size,'Kr0_bond',
     >                        Kr0_bond(2),Kr0_bond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm bond',0,MA_ERR)

            nbindxsize = 0
            nKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nbond(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_bond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,2,
     >                int_mb(indx_bond(1)+nbindxsize)))
     >               int_mb(indx_bond(1)+nbindxsize) = 0
                  nbindxsize = nbindxsize + 2

                  rtdb_name = 'pspw_qmmm_frag_Kr0_bond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                dbl_mb(Kr0_bond(1)+nKr0size)))
     >               dbl_mb(Kr0_bond(1)+nKr0size) = 0.0d0
                  nKr0size = nKr0size + 2

               end do
            end do

         end if



*        *** more angle ****
         angle_found = (naindxsize.gt.0)
         if (angle_found) then
            value = MA_alloc_get(mt_int,naindxsize,'indx_angle',
     >                        indx_angle(2),indx_angle(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,nQr0size,'Kr0_angle',
     >                        Kr0_angle(2),Kr0_angle(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm angle',0,MA_ERR)

            naindxsize = 0
            nQr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nangle(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_angle:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,3,
     >                int_mb(indx_angle(1)+naindxsize)))
     >               int_mb(indx_angle(1)+naindxsize) = 0
                  naindxsize = naindxsize + 3

                  rtdb_name = 'pspw_qmmm_frag_Kr0_angle:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,2,
     >                dbl_mb(Kr0_angle(1)+nQr0size)))
     >               dbl_mb(Kr0_angle(1)+nQr0size) = 0.0d0
                  nQr0size = nQr0size + 2

               end do
            end do

         end if

*        *** more cbond ****
         cbond_found = (ncbindxsize.gt.0)
         if (cbond_found) then
            value = MA_alloc_get(mt_int,ncbindxsize,'indx_cbond',
     >                        indx_cbond(2),indx_cbond(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,ncKr0size,'Kr0_cbond',
     >                        Kr0_cbond(2),Kr0_cbond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm cbond',0,MA_ERR)

            ncbindxsize = 0
            ncKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ncbond(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_cbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,4,
     >                int_mb(indx_cbond(1)+ncbindxsize)))
     >               int_mb(indx_cbond(1)+ncbindxsize) = 0
                  ncbindxsize = ncbindxsize + 4

                  rtdb_name = 'pspw_qmmm_frag_Kr0_cbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,3,
     >                dbl_mb(Kr0_cbond(1)+ncKr0size)))
     >               dbl_mb(Kr0_cbond(1)+ncKr0size) = 0.0d0
                  ncKr0size = ncKr0size + 3

               end do
            end do

         end if

*        *** more mbond ****
         mbond_found = (nmbindxsize.gt.0)
         if (mbond_found) then
            value = MA_alloc_get(mt_int,nmbindxsize,'indx_mbond',
     >                        indx_mbond(2),indx_mbond(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,nmKr0size,'Kr0_mbond',
     >                        Kr0_mbond(2),Kr0_mbond(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm mbond',0,MA_ERR)

            nmbindxsize = 0
            nmKr0size   = 0
            do ia=1,nkfrag
               do ja=1,int_mb(nmbond(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_mbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,2,
     >                int_mb(indx_mbond(1)+nmbindxsize)))
     >               int_mb(indx_mbond(1)+nmbindxsize) = 0
                  nmbindxsize = nmbindxsize + 2

                  rtdb_name = 'pspw_qmmm_frag_Kr0_mbond:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,3,
     >                dbl_mb(Kr0_mbond(1)+nmKr0size)))
     >               dbl_mb(Kr0_mbond(1)+nmKr0size) = 0.0d0
                  nmKr0size = nmKr0size + 3

               end do
            end do

         end if

*        *** more dihedral_spring ****
         dihedral_spring_found = (ndindxsize.gt.0)
         if (dihedral_spring_found) then
            value = MA_alloc_get(mt_int,ndindxsize,
     >              'indx_dihedral_spring',
     >               indx_dihedral_spring(2),
     >               indx_dihedral_spring(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,ndKr0size,
     >                           'Kr0_dihedral_spring',
     >                            Kr0_dihedral_spring(2),
     >                            Kr0_dihedral_spring(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory for qmmm dihedral_spring',
     >       0,MA_ERR)

            ndindxsize = 0
            ndKr0size  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ndihedral_spring(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_dihedral_spring:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,2,
     >               int_mb(indx_dihedral_spring(1)+ndindxsize)))
     >               int_mb(indx_dihedral_spring(1)+ndindxsize) = 0
                  ndindxsize = ndindxsize + 4

                  rtdb_name = 'pspw_qmmm_frag_Kr0_dihedral_spring:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,3,
     >               dbl_mb(Kr0_dihedral_spring(1)+ndKr0size))) 
     >               dbl_mb(Kr0_dihedral_spring(1)+ndKr0size) = 0.0d0
                  ndKr0size = ndKr0size + 2

               end do
            end do
         end if

*        *** more dihedral ****
         dihedral_found = (nddindxsize.gt.0)
         if (dihedral_found) then
            value = MA_alloc_get(mt_int,nddindxsize,
     >              'indx_dihedral',
     >               indx_dihedral(2),
     >               indx_dihedral(1))
            value = value.and.
     >              MA_alloc_get(mt_dbl,nddKr0size,
     >                           'Kr0_dihedral',
     >                            Kr0_dihedral(2),
     >                            Kr0_dihedral(1))
            if (.not.value) call errquit( 
     >       'cannot allocate heap memory for qmmm dihedral',
     >       0,MA_ERR)

            nddindxsize = 0
            nddKr0size  = 0
            do ia=1,nkfrag
               do ja=1,int_mb(ndihedral(1)+ia-1)

                  rtdb_name = 'pspw_qmmm_frag_indx_dihedral:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_int,2,
     >               int_mb(indx_dihedral(1)+nddindxsize)))
     >               int_mb(indx_dihedral(1)+nddindxsize) = 0
                  nddindxsize = nddindxsize + 4

                  rtdb_name = 'pspw_qmmm_frag_Kr0_dihedral:'
     >                       //c_index_name(ia)//c_index_name(ja)
                  if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,3,
     >               dbl_mb(Kr0_dihedral(1)+nddKr0size)))      
     >               dbl_mb(Kr0_dihedral(1)+nddKr0size) = 0.0d0
                  nddKr0size = nddKr0size + 3

               end do
            end do
         end if


*       ***** initialize Model, LJ, Pol, VIB, and CAV ****
        if (.not.auxiliary_only) call pspw_LJ_init(rtdb)

*       ***** make sure mm fragments are in cell ****
        if (periodic) then
           rtdb_name = 'nwpw:qmmm_fragcell'
           if (.not.rtdb_get(rtdb,rtdb_name,mt_log,1,value)) 
     >        value = .true.
           if (value) then
              call Parallel_taskid(taskid)
              if (taskid.eq.MASTER) then
                 write(*,*) "Making sure that fragments are connnected"
                 write(*,*) " - to turn this option off: ",
     >                      "set nwpw:qmmm_fragcell .false."
              end if
              call pspw_qmmm_fragcell(dbl_mb(ion_rion_ptr()))
           end if
           call pspw_qmmm_incell1(dbl_mb(ion_rion_ptr()))
        end if


*       ***** write out header info ****
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
           do ia=1,nkfrag
             write(*,*) "fragment ",ia
             write(*,'(1x,A,I4)') "  - fragment size = ",
     >                         int_mb(size_frag(1)+ia-1)
             if (log_mb(self_interaction(1)+ia-1)) then
                write(*,*) "  - ionic fragment"
             else
                write(*,*) "  - covalent fragment"
             end if
             if (log_mb(incell_frag(1)+ia-1)) then
                write(*,*) "  - incell fragment"
             else
                write(*,*) "  - not incell fragment"
             end if
             write(*,*) "  - Switching Parameters = ",
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
                   write(*,'(A,5(I4,F8.3))')   
     >                "   - shake = ",
     >                (int_mb(indx_ptr+j-1),
     >                 dsqrt(dbl_mb(shake_ptr+j-1)),
     >                 j=1,nbtmp)
             else
                   write(*,'(A,5(I4,F8.3))')   
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

             write(*,'(A,2I4,4x,A,F11.6,4x,A,F11.6)') 
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

             write(*,'(A,3I4,4x,A,F11.6,4x,A,F11.6)') 
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

             write(*,'(A,A,2I4,4x,A,F11.6, 
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

             write(*,'(A,A,2I4,3(4x,A,F11.6))') 
     >       "   - Morse Parameters: ",
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,2),
     >                  "De=", dbl_mb(shake_ptr),
     >                  "rho=",dbl_mb(shake_ptr+1),
     >                  "re=", dbl_mb(shake_ptr+2)
               ii = ii + 2
               jj = jj + 3
             end do

             ii = int_mb(nindx_dihedral_spring_start(1)+ia-1)
             jj = int_mb(nKr0_dihedral_spring_start(1) +ia-1)
             do ja=1,int_mb(ndihedral_spring(1)+ia-1)
                indx_ptr  = indx_dihedral_spring(1)+ii-1
                shake_ptr = Kr0_dihedral_spring(1) +jj-1

             write(*,'(A,A,4I4,2(4x,A,F11.6))')
     >       "   - Dihedral Spring Parameters: ",
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,4),
     >                  "Kphi=", dbl_mb(shake_ptr),
     >                  "phie=",dbl_mb(shake_ptr+1)
               ii = ii + 4
               jj = jj + 2
             end do

             ii = int_mb(nindx_dihedral_start(1)+ia-1)
             jj = int_mb(nKr0_dihedral_start(1) +ia-1)
             do ja=1,int_mb(ndihedral(1)+ia-1)
                indx_ptr  = indx_dihedral(1)+ii-1
                shake_ptr = Kr0_dihedral(1) +jj-1

             write(*,'(A,A,4I4,3(4x,A,F11.6))')
     >       "   - Dihedral Parameters: ",
     >                  "indexes=",
     >                  (int_mb(indx_ptr+j-1),j=1,4),
     >                  "Vphi=", dbl_mb(shake_ptr),
     >                  "nphi=",dbl_mb(shake_ptr+1),
     >                  "phie=", dbl_mb(shake_ptr+2)
               ii = ii + 4
               jj = jj + 3
             end do

           end do
        end if

 999    continue
        call pspw_2qmmm_init(rtdb)
      end if

      return
      end


*     **********************************
*     *	                               *
*     *         pspw_qmmm_end 	       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

      logical value
      integer  control_version
      external control_version
      logical  pspw_2qmmm_found
      external pspw_2qmmm_found

      if (qmmm_found) then
        if (nfrag.gt.0) then
        value =           MA_free_heap(indx_frag_start(2))
        value = value.and.MA_free_heap(size_frag(2))
        value = value.and.MA_free_heap(kfrag(2))
        value = value.and.MA_free_heap(self_interaction(2))
        value = value.and.MA_free_heap(incell_frag(2))
        value = value.and.MA_free_heap(switch_Rin(2))
        value = value.and.MA_free_heap(switch_Rout(2))
        value = value.and.MA_free_heap(nshake(2))
        value = value.and.MA_free_heap(nab_shake_start(2))
        value = value.and.MA_free_heap(nindx_shake_start(2))
        value = value.and.MA_free_heap(ndsq_shake_start(2))
        if (shake_found) then
          value = value.and.MA_free_heap(na(2))
          value = value.and.MA_free_heap(nb(2))
          value = value.and.MA_free_heap(indx_shake(2))
          value = value.and.MA_free_heap(dsq_shake(2))
        end if

        value = value.and.MA_free_heap(nbond(2))
        value = value.and.MA_free_heap(nindx_bond_start(2))
        value = value.and.MA_free_heap(nKr0_bond_start(2))
        if (bond_found) then
          value = value.and.MA_free_heap(indx_bond(2))
          value = value.and.MA_free_heap(Kr0_bond(2))
        end if

        value = value.and.MA_free_heap(nangle(2))
        value = value.and.MA_free_heap(nindx_angle_start(2))
        value = value.and.MA_free_heap(nKr0_angle_start(2))
        if (angle_found) then
          value = value.and.MA_free_heap(indx_angle(2))
          value = value.and.MA_free_heap(Kr0_angle(2))
        end if

        value = value.and.MA_free_heap(ncbond(2))
        value = value.and.MA_free_heap(nindx_cbond_start(2))
        value = value.and.MA_free_heap(nKr0_cbond_start(2))
        if (cbond_found) then
          value = value.and.MA_free_heap(indx_cbond(2))
          value = value.and.MA_free_heap(Kr0_cbond(2))
        end if

        value = value.and.MA_free_heap(nmbond(2))
        value = value.and.MA_free_heap(nindx_mbond_start(2))
        value = value.and.MA_free_heap(nKr0_mbond_start(2))
        if (mbond_found) then
          value = value.and.MA_free_heap(indx_mbond(2))
          value = value.and.MA_free_heap(Kr0_mbond(2))
        end if

        value = value.and.MA_free_heap(ndihedral_spring(2))
        value = value.and.MA_free_heap(nindx_dihedral_spring_start(2))
        value = value.and.MA_free_heap(nKr0_dihedral_spring_start(2))
        if (dihedral_spring_found) then
          value = value.and.MA_free_heap(indx_dihedral_spring(2))
          value = value.and.MA_free_heap(Kr0_dihedral_spring(2))
        end if

        value = value.and.MA_free_heap(ndihedral(2))
        value = value.and.MA_free_heap(nindx_dihedral_start(2))
        value = value.and.MA_free_heap(nKr0_dihedral_start(2))
        if (dihedral_found) then
          value = value.and.MA_free_heap(indx_dihedral(2))
          value = value.and.MA_free_heap(Kr0_dihedral(2))
        end if



        if (.not.auxiliary_only) then
        if (control_version().ne.3) then
        value = value.and.MA_free_heap(rcell(2))
        end if
        end if

        if (.not.value) call errquit('cannot free heap memory',0,
     &       MA_ERR)

        end if

        if (pspw_2qmmm_found())  call pspw_2qmmm_end()
        if (.not.auxiliary_only) call pspw_LJ_end()
      end if

      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_found        *
*     *                                *
*     **********************************
      logical function pspw_qmmm_found()
      implicit none

#include "pspw_qmmm.fh"

      pspw_qmmm_found = qmmm_found
      return
      end


*     **********************************
*     *                                *
*     *    pspw_qmmm_shake_found       *
*     *                                *
*     **********************************
      logical function pspw_qmmm_shake_found()
      implicit none

#include "pspw_qmmm.fh"

      pspw_qmmm_shake_found = shake_found
      return
      end


*     **********************************
*     *	                               *
*     *         pspw_qmmm_LJ_E 	       *
*     *                                *
*     **********************************

      real*8 function pspw_qmmm_LJ_E()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      real*8   pspw_LJ_E
      external pspw_LJ_E

      if (auxiliary_only) then
         pspw_qmmm_LJ_E = 0.0d0
      else
         pspw_qmmm_LJ_E = pspw_LJ_E(ion_nion(),
     >                              ion_nion_qm(),
     >                              int_mb(ion_katm_ptr()),
     >                              nfrag,
     >                              int_mb(indx_frag_start(1)),
     >                              int_mb(size_frag(1)),
     >                              int_mb(kfrag(1)),
     >                              log_mb(self_interaction(1)),
     >                              nshl3d,dbl_mb(rcell(1)),
     >                              dbl_mb(ion_rion_ptr()))
      end if
      return
      end



*     **********************************
*     *                                *
*     *         pspw_qmmm_Q_E          *
*     *                                *
*     **********************************
      real*8 function pspw_qmmm_Q_E()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  psp_zv_ptr,ion_amass_ptr
      external psp_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      if (auxiliary_only) then
         pspw_qmmm_Q_E = 0.0d0
      else
         pspw_qmmm_Q_E = pspw_Q_E(ion_nion(),
     >                            ion_nion_qm(),
     >                            int_mb(ion_katm_ptr()),
     >                            dbl_mb(psp_zv_ptr()),
     >                            dbl_mb(ion_amass_ptr()),
     >                            nfrag,
     >                            int_mb(indx_frag_start(1)),
     >                            int_mb(size_frag(1)),
     >                            int_mb(kfrag(1)),
     >                            dbl_mb(switch_Rin(1)),
     >                            dbl_mb(switch_Rout(1)),
     >                            log_mb(self_interaction(1)),
     >                            nshl3d,dbl_mb(rcell(1)),
     >                            dbl_mb(ion_rion_ptr()))
      end if
      return
      end


*     **********************************
*     *                                *
*     *         pspw_qmmm_mmq_Q_E          *
*     *                                *
*     **********************************

      real*8 function pspw_qmmm_mmq_Q_E()
      implicit none  

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"
         
c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  mmq_zv_ptr,ion_amass_ptr
      external mmq_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      if (auxiliary_only) then
         pspw_qmmm_mmq_Q_E = 0.0d0
      else
         pspw_qmmm_mmq_Q_E = pspw_Q_E(ion_nion(),
     >                            ion_nion_qm(),
     >                            int_mb(ion_katm_ptr()),
     >                            dbl_mb(mmq_zv_ptr()),
     >                            dbl_mb(ion_amass_ptr()),
     >                            nfrag,
     >                            int_mb(indx_frag_start(1)),
     >                            int_mb(size_frag(1)),
     >                            int_mb(kfrag(1)),
     >                            dbl_mb(switch_Rin(1)),
     >                            dbl_mb(switch_Rout(1)),
     >                            log_mb(self_interaction(1)),
     >                            nshl3d,dbl_mb(rcell(1)),
     >                            dbl_mb(ion_rion_ptr()))
      end if
      return
      end




*     **********************************
*     *                                *
*     *        pspw_qmmm_spring_E      *
*     *                                *
*     **********************************

      real*8 function pspw_qmmm_spring_E()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** local variables ****
      integer dutask,taskid,np
      integer w1,wk1,ks1,nbs,nas,ncbs,nbm,jj,kk,nds,ndds
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
      real*8   pspw_2qmmm_spring_E
      external pspw_2qmmm_spring_E
      logical  pspw_2qmmm_found
      external pspw_2qmmm_found

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
            E = E + pspw_qmmm_dihedral_spring_frag(nbm,
     >                        int_mb(indx_dihedral_spring(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral_spring(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         if (ndds.gt.0) then
            jj  = int_mb(nindx_dihedral_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_start(1) +wk1-1)
            E = E + pspw_qmmm_dihedral_frag(nbm,
     >                        int_mb(indx_dihedral(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)))
         end if
         end if
         dutask = mod(dutask+1,np)
        
      end do
      if (np.gt.1) call Parallel_SumAll(E)

c     **** non-fragments ****
      if (pspw_2qmmm_found()) E = E + pspw_2qmmm_spring_E()

      pspw_qmmm_spring_E = E
      return
      end


*     **********************************
*     *                                *
*     *      pspw_qmmm_spring_fion     *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_spring_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** local variables ****
      integer dutask,taskid,np
      integer w1,wk1,ks1,nbs,nas,ncbs,nbm,jj,kk,nds,ndds

c     **** external functions ***
      integer  ion_rion_ptr
      external ion_rion_ptr
      logical  pspw_2qmmm_found
      external pspw_2qmmm_found

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      dutask = 0

      do w1=1,nfrag
         if (dutask.eq.taskid) then
         ks1  = int_mb(indx_frag_start(1)+w1-1)
         wk1  = int_mb(kfrag(1)+w1-1)
         nbs  = int_mb(nbond(1)+wk1-1)
         nas  = int_mb(nangle(1)+wk1-1)
         ncbs = int_mb(ncbond(1)+wk1-1)
         nbm  = int_mb(nmbond(1)+wk1-1)
         nds  = int_mb(ndihedral_spring(1)+wk1-1)
         ndds  = int_mb(ndihedral(1)+wk1-1)
         if (nbs.gt.0) then
            jj  = int_mb(nindx_bond_start(1)+wk1-1)
            kk  = int_mb(nKr0_bond_start(1) +wk1-1)
            call pspw_qmmm_spring_bond_frag_fion(nbs,
     >                                int_mb(indx_bond(1)+jj-1),
     >                                dbl_mb(Kr0_bond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                fion(1,ks1))
         end if
         if (nas.gt.0) then
            jj  = int_mb(nindx_angle_start(1)+wk1-1)
            kk  = int_mb(nKr0_angle_start(1) +wk1-1)
            call pspw_qmmm_spring_angle_frag_fion(nas,
     >                                int_mb(indx_angle(1)+jj-1),
     >                                dbl_mb(Kr0_angle(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                fion(1,ks1))
         end if
         if (ncbs.gt.0) then
            jj  = int_mb(nindx_cbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_cbond_start(1) +wk1-1)
            call pspw_qmmm_spring_cbond_frag_fion(ncbs,
     >                                int_mb(indx_cbond(1)+jj-1),
     >                                dbl_mb(Kr0_cbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                fion(1,ks1))
         end if
         if (nbm.gt.0) then
            jj  = int_mb(nindx_mbond_start(1)+wk1-1)
            kk  = int_mb(nKr0_mbond_start(1) +wk1-1)
            call pspw_qmmm_morse_bond_frag_fion(nbm,
     >                                int_mb(indx_mbond(1)+jj-1),
     >                                dbl_mb(Kr0_mbond(1) +kk-1),
     >                                dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                                fion(1,ks1))
         end if
         if (nds.gt.0) then
            jj  = int_mb(nindx_dihedral_spring_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_spring_start(1) +wk1-1)
            call pspw_qmmm_dihedral_spring_frag_fion(nds,
     >                        int_mb(indx_dihedral_spring(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral_spring(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        fion(1,ks1))
         end if
         if (ndds.gt.0) then
            jj  = int_mb(nindx_dihedral_start(1)+wk1-1)
            kk  = int_mb(nKr0_dihedral_start(1) +wk1-1)
            call pspw_qmmm_dihedral_frag_fion(ndds,
     >                        int_mb(indx_dihedral(1)+jj-1),
     >                        dbl_mb(Kr0_dihedral(1) +kk-1),
     >                        dbl_mb(ion_rion_ptr()+3*(ks1-1)),
     >                        fion(1,ks1))
         end if


         end if
         dutask = mod(dutask+1,np)
      end do

c     *** non-fragments ***
      if (pspw_2qmmm_found()) call pspw_2qmmm_spring_fion(fion)
      return
      end



*     **********************************
*     *                                *
*     *         pspw_qmmm_LJ_fion      *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_LJ_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr

      if (.not.auxiliary_only) then
         call  pspw_LJ_fion(ion_nion(),
     >                      ion_nion_qm(),
     >                      int_mb(ion_katm_ptr()),
     >                      nfrag,
     >                      int_mb(indx_frag_start(1)),
     >                      int_mb(size_frag(1)),
     >                      int_mb(kfrag(1)),
     >                      log_mb(self_interaction(1)),
     >                      nshl3d,dbl_mb(rcell(1)),
     >                      dbl_mb(ion_rion_ptr()),fion)
      end if
      return
      end



*     **********************************
*     *                                *
*     *         pspw_qmmm_Q_fion       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_Q_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  psp_zv_ptr,ion_amass_ptr
      external psp_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      if (.not.auxiliary_only) then
         call pspw_Q_fion(ion_nion(),
     >                    ion_nion_qm(),
     >                    int_mb(ion_katm_ptr()),
     >                    dbl_mb(psp_zv_ptr()),
     >                    dbl_mb(ion_amass_ptr()),
     >                    nfrag,
     >                    int_mb(indx_frag_start(1)),
     >                    int_mb(size_frag(1)),
     >                    int_mb(kfrag(1)),
     >                    dbl_mb(switch_Rin(1)),
     >                    dbl_mb(switch_Rout(1)),
     >                    log_mb(self_interaction(1)),
     >                    nshl3d,dbl_mb(rcell(1)),
     >                    dbl_mb(ion_rion_ptr()),fion)
      end  if
      return
      end

*     **********************************
*     *                                *
*     *         pspw_qmmm_mmq_Q_fion   *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_mmq_Q_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  mmq_zv_ptr,ion_amass_ptr
      external mmq_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      if (.not.auxiliary_only) then
         call pspw_Q_fion(ion_nion(),
     >                    ion_nion_qm(),
     >                    int_mb(ion_katm_ptr()),
     >                    dbl_mb(mmq_zv_ptr()),
     >                    dbl_mb(ion_amass_ptr()),
     >                    nfrag,
     >                    int_mb(indx_frag_start(1)),
     >                    int_mb(size_frag(1)),
     >                    int_mb(kfrag(1)),
     >                    dbl_mb(switch_Rin(1)),
     >                    dbl_mb(switch_Rout(1)),
     >                    log_mb(self_interaction(1)),
     >                    nshl3d,dbl_mb(rcell(1)),
     >                    dbl_mb(ion_rion_ptr()),fion)
      end  if
      return
      end




*     **********************************
*     *                                *
*     *         pspw_qmmm_fion         *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*      **** local variables ****
      integer ftmp(2),nion3,np

*     *** external functions ****
      integer  ion_nion
      external ion_nion
      call Parallel_np(np)

      if (np.eq.1) then
         if (.not.auxiliary_only) then
            call pspw_qmmm_LJ_fion(fion)
            call pspw_qmmm_Q_fion(fion)
         end if
         call pspw_qmmm_spring_fion(fion)
      else
         nion3 = 3*ion_nion()
         if (.not.MA_push_get(mt_dbl,nion3,'ftmp',ftmp(2),ftmp(1)))
     >  call errquit('pspw_qmmm_fion:cannot push stack memory',0,MA_ERR)

         call dcopy(nion3,0.0,0,dbl_mb(ftmp(1)),1)
         if (.not.auxiliary_only) then
            call pspw_qmmm_LJ_fion(dbl_mb(ftmp(1)))
            call pspw_qmmm_Q_fion(dbl_mb(ftmp(1)))
         end if
         call pspw_qmmm_spring_fion(dbl_mb(ftmp(1)))
         call Parallel_Vector_SumAll(nion3,dbl_mb(ftmp(1)))

         call daxpy(nion3,1.0d0,dbl_mb(ftmp(1)),1,fion,1)


         if (.not.MA_pop_stack(ftmp(2)))
     >   call errquit('pspw_qmmm_fion:cannot pop stack memory',0,MA_ERR)
      end if
      return
      end


*     **********************************
*     *                                *
*     *         pspw_qmmm_mmq_fion     *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_mmq_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*      **** local variables ****
      integer ftmp(2),nion3,np

*     *** external functions ****
      integer  ion_nion
      external ion_nion
      call Parallel_np(np)

      if (np.eq.1) then
         if (.not.auxiliary_only) then
            call pspw_qmmm_LJ_fion(fion)
            call pspw_qmmm_mmq_Q_fion(fion)
         end if
         call pspw_qmmm_spring_fion(fion)
      else
         nion3 = 3*ion_nion()
         if (.not.MA_push_get(mt_dbl,nion3,'ftmp',ftmp(2),ftmp(1)))
     >  call errquit('pspw_qmmm_fion:cannot push stack memory',0,MA_ERR)

         call dcopy(nion3,0.0,0,dbl_mb(ftmp(1)),1)
         if (.not.auxiliary_only) then
            call pspw_qmmm_LJ_fion(dbl_mb(ftmp(1)))
            call pspw_qmmm_mmq_Q_fion(dbl_mb(ftmp(1)))
         end if
         call pspw_qmmm_spring_fion(dbl_mb(ftmp(1)))
         call Parallel_Vector_SumAll(nion3,dbl_mb(ftmp(1)))

         call daxpy(nion3,1.0d0,dbl_mb(ftmp(1)),1,fion,1)

         if (.not.MA_pop_stack(ftmp(2)))
     >   call errquit('pspw_qmmm_fion:cannot pop stack memory',0,MA_ERR)
      end if
      return
      end


*     ******************************
*     *                            *
*     *       pspw_qmmm_frag_cm    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_frag_cm(size_frag,kfrag_s,rion,rcm,mtotal)
      implicit none
      integer size_frag,kfrag_s
      real*8 rion(3,*),rcm(3),mtotal

#include "mafdecls.fh"

*     **** local variables ****
      integer k,kk,ptr
      real*8  m

*     **** external functions ****
      integer  ion_amass_ptr
      external ion_amass_ptr

      ptr = ion_amass_ptr()

      mtotal = 0.0d0
      rcm(1) = 0.0d0
      rcm(2) = 0.0d0
      rcm(3) = 0.0d0
      kk = kfrag_s
      do k=1,size_frag
         m = dbl_mb(ptr+kk-1)
         rcm(1) = rcm(1) + m*rion(1,kk)
         rcm(2) = rcm(2) + m*rion(2,kk)
         rcm(3) = rcm(3) + m*rion(3,kk)
         mtotal = mtotal + m
         kk = kk + 1
      end do
      rcm(1) = rcm(1)/mtotal
      rcm(2) = rcm(2)/mtotal
      rcm(3) = rcm(3)/mtotal

      return
      end

*     ******************************
*     *                            *
*     *       pspw_qmmm_fragcell   *
*     *                            *
*     ******************************
*
      subroutine pspw_qmmm_fragcell(r1)
      implicit none
      real*8 r1(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1

      do w1=1,nfrag
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call lattice_fragcell(n1,r1(1,ks1))
      end do
      return
      end




*     ******************************
*     *                            *
*     *       pspw_qmmm_incell1    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell1(r1)
      implicit none
      real*8 r1(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1,ia
      real*8  rcm(3),mtotal

      do w1=1,nfrag
c        ia = int_mb(kfrag(1)+w1-1)
c        if (log_mb(incell_frag(1)+ia-1)) then
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell1_frag(rcm,n1,r1(1,ks1))
c        end if
      end do
      return
      end 

*     ******************************
*     *                            *
*     *       pspw_qmmm_incell2    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell2(r1,r2)
      implicit none
      real*8 r1(3,*),r2(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1,ia
      real*8  rcm(3),mtotal

      do w1=1,nfrag
c        ia = int_mb(kfrag(1)+w1-1)
c        if (log_mb(incell_frag(1)+ia-1)) then
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell2_frag(rcm,n1,r1(1,ks1),r2(1,ks1))
c        end if
      end do
      return
      end


*     ******************************
*     *                            *
*     *       pspw_qmmm_incell3    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell3(r1,r2,r3)
      implicit none
      real*8 r1(3,*),r2(3,*),r3(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1,ia
      real*8  rcm(3),mtotal


      do w1=1,nfrag
c        ia = int_mb(kfrag(1)+w1-1)
c        if (log_mb(incell_frag(1)+ia-1)) then
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell3_frag(rcm,n1,r1(1,ks1),r2(1,ks1),r3(1,ks1))
c        end if
      end do

      return
      end



*     ******************************
*     *                            *
*     *       pspw_qmmm_shake      *
*     *                            *
*     ******************************
*   shakes the fragments 
*
      subroutine pspw_qmmm_shake(r2,r1)
      implicit none
      real*8 r2(3,*),r1(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer maxit
      double precision tol
      parameter (maxit=1000,tol=1.0d-4)

      integer w1,wk1,ks1
      integer nshk,ii,jj,kk

*     **** external functions ****
      integer  ion_amass_ptr
      external ion_amass_ptr
      logical  pspw_2qmmm_found
      external pspw_2qmmm_found

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
  
c     **** non-frags ****
      if (pspw_2qmmm_found()) call pspw_2qmmm_shake(r2,r1)

      return
      end

*     ******************************
*     *                            *
*     *   pspw_qmmm_shake_frag     *
*     *                            *
*     ******************************
*   shakes a fragment 
*
      subroutine pspw_qmmm_shake_frag(nshk,na,nb,indx,dsq,
     >                                tol,maxit,
     >                                mass,r2,r1)
      implicit none
      integer nshk, na(*),nb(*)
      integer indx(*)
      real*8  dsq(*)
      real*8  tol
      integer maxit
      real*8  mass(*)
      real*8 r2(*),r1(*)

*     **** local variables ****
      integer ii,jj,sk

      ii = 1
      jj = 1
      do sk=1,nshk
         call shake_chain3(na(sk),indx(ii),nb(sk),
     >                     tol,maxit,
     >                     dsq(jj),
     >                     mass,r2,r1)
         ii = ii + na(sk)
         jj = jj + nb(sk)
      end do
      return
      end

