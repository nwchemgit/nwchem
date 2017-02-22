      subroutine mp2_make_cphf_rhs(scftype,transform,basis,nmo,nbf,
     $     noa,nva,nob,nvb,irs_a,irs_b,map_a,map_b,
     $     num_oa,num_ob,num_va,num_vb,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,
     $     sym_lo_ob,sym_hi_ob,sym_lo_vb,sym_hi_vb,
     $     g_vecs_a,g_vecs_b,g_lai_a,g_lai_b,g_rhs)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      character*(*) scftype
      logical transform
      integer irs_a(*),irs_b(*),map_a(*),map_b(*)
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer num_oa(0:*),num_va(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer nmo,nbf,noa,nva,nob,nvb
      integer basis
      character*256 cphf_rhs
      integer g_vecs_a,g_vecs_b,g_lai_a,g_lai_b
      integer g_rhs
      logical file_write_ga
      external file_write_ga
      integer g_lai_a_tr,g_lai_b_tr
c     
      call cphf_fname('cphf_rhs',cphf_rhs)
*ga:1:0
      if(.not.ga_create(mt_dbl,nva,noa,'temp',0,0,g_lai_a_tr))
     $     call errquit('could not alloc tran array',1, GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl,nvb,nob,'temp',0,0,g_lai_b_tr))
     $     call errquit('could not alloc tran array',1, GA_ERR)
      call ga_transpose(g_lai_a,g_lai_a_tr)
      call ga_transpose(g_lai_b,g_lai_b_tr)
      call cphf_rhs_reshape(scftype, 
     $     g_rhs,g_lai_a_tr,g_lai_b_tr,noa,nva,
     $     nob,nvb,nmo,map_a,map_b)
      if(.not.file_write_ga(cphf_rhs,g_rhs))
     $     call errquit('could not write cphf_rhs',1, DISK_ERR)
      if(.not.ga_destroy(g_lai_a_tr))call errquit(
     $     'could not destroy lai_a_tr handle',1, GA_ERR)
      if(.not.ga_destroy(g_lai_b_tr))call errquit(
     $     'could not destroy lai_b_tr handle',1, GA_ERR)
      end
c     
      subroutine cphf_rhs_reshape(scftype,g_rhs,g_lai_a,g_lai_b,noa,
     $     nva,nob,nvb,nmo,map_a,map_b)
      implicit none
#include "global.fh"
      character*(*) scftype
      integer g_rhs,g_lai_a,g_lai_b
      integer noa,nva,nob,nvb,nmo
      integer ioff
      integer map_a(*), map_b(*)
c     
      integer i, a, ii, aa, ia
      double precision xx
c     
*     write(6,*) ' map_a ', (map_a(i),i=1,noa+nva)
*     write(6,*) ' map_b ', (map_b(i),i=1,nob+nvb)
c     
      ioff = 0
      do i = 1+ga_nodeid(),noa,ga_nnodes() ! Loop thru SCF order
         do a = 1,nva
            ii = map_a(i)
            aa = map_a(a+noa)-noa
            call ga_get(g_lai_a,aa,aa,ii,ii,xx,1)
            ia = ioff+a+(i-1)*nva
            call ga_put(g_rhs,ia,ia,1,1,xx,1)
         end do
      end do
c     
c$$$  call ga_copy_patch('n',g_lai_a,1,nva,1,noa,
c$$$  $     g_rhs,1,nva*noa,1,1)
      if (scftype .eq. 'UHF') then
         ioff=nva*noa
         do i = 1+ga_nodeid(),nob,ga_nnodes() ! Loop thru SCF order
            do a = 1,nvb
               ii = map_b(i)
               aa = map_b(a+nob)-nob
               call ga_get(g_lai_b,aa,aa,ii,ii,xx,1)
               ia = ioff+a+(i-1)*nvb
               call ga_put(g_rhs,ia,ia,1,1,xx,1)
            end do
         end do
c$$$  call ga_copy_patch('n',g_lai_b,1,nvb,1,nob,
c$$$  $        g_rhs,ioff+1,ioff+(nvb*nob),1,1)
      else
         call ga_dscal(g_rhs, 2.0d0)
      end if
c     
      end
