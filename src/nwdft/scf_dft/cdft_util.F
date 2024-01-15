      subroutine cdft_init(iatms1,iatms2,fatms1,fatms2,g_constr,nconstr,
     $                     ipop,dtype,constr, geom, ncenters,ao_bas_han,
     $                     nbf_ao, g_s, g_shalf)
      implicit none
      integer nconstr, ipop
      integer iatms1(nconstr),fatms1(nconstr)
      integer iatms2(nconstr),fatms2(nconstr)
      integer dtype(nconstr),g_constr(nconstr)
      double precision constr(nconstr)
      integer geom, ncenters, ao_bas_han, nbf_ao 
      integer g_s, g_shalf
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "geom.fh"
      integer i, j
      integer g_x, g_t
      double precision coor(3), nuc_charge, charge
      character*16 t
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

c convert charge to number of electrons
      do i = 1, nconstr
c
c       Check the atom centers on the constraints
c
        if (iatms1(i).lt.1) then
          if (ga_nodeid().eq.0) then
            write(LuOut,*)'cdft_init: constraint no.',nconstr,
     &                    ' lower atom number is ',iatms1(i),
     &                    ' but should be at least 1'
            call util_flush(LuOut)
          endif
          call errquit('cdft_init: lower center no. out of range',
     &                 iatms1(i),INPUT_ERR)
        endif
        if (iatms2(i).ne.0.and.iatms2(i).lt.1) then
          if (ga_nodeid().eq.0) then
            write(LuOut,*)'cdft_init: constraint no.',nconstr,
     &                    ' lower atom number is ',iatms2(i),
     &                    ' but should be at least 1'
            call util_flush(LuOut)
          endif
          call errquit('cdft_init: lower center no. out of range',
     &                 iatms2(i),INPUT_ERR)
        endif
        if (fatms1(i).gt.ncenters) then
          if (ga_nodeid().eq.0) then
            write(LuOut,*)'cdft_init: constraint no.',nconstr,
     &                    ' upper atom number is ',fatms1(i),
     &                    ' but should be at most ',ncenters
            call util_flush(LuOut)
          endif
          call errquit('cdft_init: upper center no. out of range',
     &                 fatms1(i),INPUT_ERR)
        endif
        if (iatms2(i).ne.0.and.fatms2(i).gt.ncenters) then
          if (ga_nodeid().eq.0) then
            write(LuOut,*)'cdft_init: constraint no.',nconstr,
     &                    ' upper atom number is ',fatms2(i),
     &                    ' but should be at most ',ncenters
            call util_flush(LuOut)
          endif
          call errquit('cdft_init: upper center no. out of range',
     &                 fatms2(i),INPUT_ERR)
        endif
c
c       Done checking
c
        if(dtype(i).eq.1) then
          charge = 0d0
          do j = iatms1(i), fatms1(i)
            if(.not.geom_cent_get(geom, j, t, coor, nuc_charge))
     $        call errquit('cdft_init: geom_cent_get fail', 0, GA_ERR)
            charge = charge + nuc_charge
          enddo
          if(iatms2(i).ne.0) then
            do j = iatms2(i), fatms2(i)
              if(.not.geom_cent_get(geom, j, t, coor, nuc_charge))
     $          call errquit('cdft_init: geom_cent_get fail', 0, GA_ERR)
              charge = charge - nuc_charge
            enddo
          endif
          constr(i) = charge - constr(i)
        endif
      enddo

      if(ipop.eq.1) return   !!!!  grid sa, done in ../grid/grid_quadv0.F

      g_t = ga_create_atom_blocked(geom, AO_bas_han, 'temp')

      if(ipop.eq.2) then    !!!!  mulliken
       do i = 1, nconstr
        call mull_sa(ao_bas_han, nbf_ao, g_s, g_constr(i), 
     $                iatms1(i),fatms1(i))
        if(iatms2(i).ne.0) then
          call mull_sa(ao_bas_han, nbf_ao, g_s, g_t, 
     $                iatms2(i),fatms2(i))
          call ga_dadd(1.d0, g_constr(i), -1.d0, g_t, g_constr(i))
        endif
       enddo
      elseif(ipop.eq.3) then    !!!!  lowdin
       do i = 1, nconstr
        call trans_ol(ao_bas_han, nbf_ao,g_constr(i),g_shalf,g_shalf,
     $                 iatms1(i),fatms1(i))
        if(iatms2(i).ne.0) then
          call trans_ol(ao_bas_han, nbf_ao,g_t,g_shalf,g_shalf,
     $                 iatms2(i),fatms2(i))
          call ga_dadd(1.d0, g_constr(i), -1.d0, g_t, g_constr(i))
        endif
       enddo
      elseif(ipop.eq.4) then    !!!!  projected lowdin
       do i = 1, nconstr
        g_x = ga_create_atom_blocked(geom, AO_bas_han, 'x-AO ovl')
        call davidson_pop(geom, ncenters, ao_bas_han, nbf_ao,
     $                       g_s, g_shalf, g_x)
        call trans_ol(ao_bas_han,nbf_ao,g_constr(i),g_shalf,g_x,
     $                 iatms1(i),fatms1(i))
        if(iatms2(i).ne.0) then
          call trans_ol(ao_bas_han, nbf_ao,g_t,g_shalf,g_shalf,
     $                 iatms2(i),fatms2(i))
          call ga_dadd(1.d0, g_constr(i), -1.d0, g_t, g_constr(i))
        endif
        if (.not. ga_destroy(g_x)) call errquit
     &     ('cdft_init: could not destroy g_x', 0, GA_ERR)
       enddo
      else
        call errquit('choice of pop not implemented', 0, UNKNOWN_ERR)
      endif

      if (.not. ga_destroy(g_t)) call errquit
     &     ('cdft_init: could not destroy g_t', 0, GA_ERR)
      return
      end


      subroutine mull_sa(bas_han, nbf, g_s, g_sa, iatom,fatom)
      implicit none
      integer bas_han, nbf, iatom, fatom  
      integer g_s, g_sa
#include "bas.fh"
#include "global.fh"
      logical status
      integer i, j, ifirst, ilast
      double precision stmp(nbf,nbf)

      status = bas_ce2bfr(bas_han, iatom, ifirst, i)
      status = bas_ce2bfr(bas_han, fatom, i, ilast)

      call ga_get(g_s,1,nbf,1,nbf,stmp,nbf)

      do i = 1, nbf
         if(i.le.ilast .and. i.ge.ifirst) then
           do j = 1, nbf
             if(j.gt.ilast .or. j.lt.ifirst) then
               stmp(i,j)=0.5d0*stmp(i,j)
             endif
           enddo
         else
           do j = 1, nbf
             if(j.le.ilast .and. j.ge.ifirst) then
               stmp(i,j)=0.5d0*stmp(i,j)
             else
               stmp(i,j)=0.d0
             endif
           enddo
         endif
      enddo

      call ga_put(g_sa,1,nbf,1,nbf,stmp,nbf)

      return
      end


      subroutine trans_ol(bas_han, nbf, g_p, g_l, g_r, iatom,fatom)
      implicit none

      integer bas_han, nbf, iatom, fatom  
      integer g_p, g_l, g_r  

#include "bas.fh"
#include "global.fh"

      logical status
      integer i, ifirst, ilast

      status = bas_ce2bfr(bas_han, iatom, ifirst, i)
      status = bas_ce2bfr(bas_han, fatom, i, ilast)

      call ga_matmul_patch('N', 'N', 1d0, 0d0,
     $                     g_l, 1, nbf, ifirst, ilast,
     $                     g_r, ifirst, ilast, 1, nbf,
     $                     g_p, 1, nbf, 1, nbf)

      return
      end


      subroutine lowd_pop( geom, basis, iga_dens,iga_s, g_tmp )

C$Id$
      Implicit none
#include "errquit.fh"
      integer geom,basis
      integer iga_s              ! overlap   GA handle
      integer iga_dens           ! dens. mat GA handle
      integer iga_ps             ! product   GA handle
      integer g_tmp             

      integer natoms,nshells
      integer lPSmat,iPSmat,lqatom,iqatom,lqshell,iqshell
      integer iatom,ilo,ihi,nbf,max_at_bf2
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical status

#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "cscfps.fh"
#include "mafdecls.fh"

      if (oscfps) call pstat_on(ps_mull)
      

c*****************************************************************************

c     
      if(.not.geom_ncent(geom, natoms))
     &     call errquit(' exiting in mull_pop',0, GEOM_ERR)
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit(' exiting in mull_pop',1, BASIS_ERR)
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in mull_op',1, BASIS_ERR)
      max_at_bf2 = 0
      do iatom = 1, natoms
        if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $       call errquit('mul_pop: bas_ce2bfr failed', iatom,
     &       BASIS_ERR)
        max_at_bf2 = max(max_at_bf2, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf2*max_at_bf2

      if(.not.MA_Push_Get(mt_dbl,max_at_bf2,'PS',lPSmat,iPSmat))
     &     call errquit(' exiting in mull_pop: insuff stack',21, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,natoms,'q atom',lqatom,iqatom))
     &     call errquit(' exiting in mull_pop: insuff stack',22, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,nshells,'q shell',lqshell,iqshell))
     &     call errquit(' exiting in mull_pop: insuff stack',3, MA_ERR)

      iga_PS=ga_create_atom_blocked(geom, basis, 'PS product')

      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_dens,iga_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_s,g_tmp,0.d0,iga_PS)
      call mull_calc(basis,natoms, nshells,max_at_bf2,iga_PS,
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat), 'total')
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33, MA_ERR)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34, MA_ERR)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35, MA_ERR)
      status= ga_destroy(iga_PS)
c
      if (oscfps) call pstat_off(ps_mull)
c
      return
      end


      subroutine get_mat12(geom, bas_han, nbf, g_m, g_m12p, g_m12m)
      implicit none
      integer geom, bas_han, nbf, g_m, g_m12p, g_m12m
#include "geom.fh"
#include "global.fh"
#include "errquit.fh"
      double precision evals(nbf)
      integer g_vecs, g_tmp
      integer i
      double precision toll_s, zero
      parameter (toll_s = 1.d-16, zero = 0.d0)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      do i = 1, nbf
         evals(i)=zero
      enddo

      g_vecs = ga_create_atom_blocked(geom, bas_han,
     &                                 'tmp eig vecs')
      call ga_zero(g_vecs)
#if defined(PARALLEL_DIAG)
      call ga_diag_std(g_m, g_vecs, evals)
#else
      call ga_diag_std_seq(g_m, g_vecs, evals)
#endif

      g_tmp = ga_create_atom_blocked(geom, bas_han, 'tmp')
      call diis_bld12(toll_s, evals, g_vecs, g_m12m, g_tmp, 2)
      call diis_bld12(toll_s, evals, g_vecs, g_m12p, g_tmp, 3)

      if(.not. ga_destroy(g_tmp)) call errquit
     $   ('get_mat12: could not destroy g_tmp', 0, GA_ERR)
      if(.not. ga_destroy(g_vecs)) call errquit
     $   ('get_mat12: could not destroy g_vecs', 0, GA_ERR)

      return
      end


      subroutine get_blocked_s(natoms, bas_han, nbf, g_s, g_sblock)
      implicit none

      integer natoms, bas_han, nbf
      integer g_s, g_sblock
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      logical status
      integer iatom, ifirst, ilast, i, j, ij, ji
      integer i_s, l_s
      double precision zero
      parameter (zero = 0.d0)

      if(.not.ma_push_get(mt_dbl,nbf*nbf,'sblock',l_s,i_s))
     &     call errquit(' get_blocked_s: insuff stack',21, MA_ERR)

      call ga_print(g_s)
      call ga_get(g_s,1,nbf,1,nbf,dbl_mb(i_s),nbf)

      do iatom = 1, natoms
         status = bas_ce2bfr(bas_han, iatom, ifirst, ilast)
         do j = ifirst, ilast
            do i = ilast+1, nbf
            ij = (j-1)*nbf+i
            ji = (i-1)*nbf+j
            dbl_mb(i_s+ij-1)=zero
            dbl_mb(i_s+ji-1)=zero
            enddo
         enddo
      enddo

      call ga_put(g_sblock,1,nbf,1,nbf,dbl_mb(i_s),nbf)
      call ga_print(g_sblock)

      if(.not.ma_pop_stack(l_s))
     &     call errquit(' exiting in get_blocked_s',35, MA_ERR)

      return
      end
      

      subroutine davidson_pop(geom, natoms, bas_han, nbf,g_s,g_sx,g_x)
      implicit none
      integer geom, natoms, bas_han, nbf, g_s, g_sx, g_x
#include "geom.fh"
#include "global.fh"
#include "errquit.fh"
      integer g_w, g_r, g_tmp1, g_tmp2
      double precision one, zero
      parameter (one = 1.d0, zero = 0.d0)
      logical status

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      g_w = ga_create_atom_blocked(geom, bas_han,'sw block')
      g_r = ga_create_atom_blocked(geom, bas_han,'sr block')
      g_tmp1 = ga_create_atom_blocked(geom, bas_han,'tmp 1')
      g_tmp2 = ga_create_atom_blocked(geom, bas_han,'tmp 2')

      call get_blocked_s(natoms, bas_han, nbf, g_s, g_w)

c  store g_sblock12p in g_tmp1, g_sblock12m in g_w
      call get_mat12(geom, bas_han, nbf, g_w, g_tmp1, g_w)

c  g_r = g_w * g_s * g_w
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_w,g_s,zero,g_tmp2)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_tmp2,g_w,zero,g_r)

c  store g_r12p in g_r, g_r12m in g_tmp2
      call get_mat12(geom, bas_han, nbf, g_r, g_r, g_tmp2)

c g_x = g_r12p * g_sblock12p
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_r,g_tmp1,zero,g_x)
      call ga_print(g_x)

c g_sx = g_s * g_sblock12m * g_r12m
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_s,g_w,zero,g_r)
      call ga_dgemm('n','n',nbf,nbf,nbf,one,g_r,g_tmp2,zero,g_sx)

      status= ga_destroy(g_w)
      status= ga_destroy(g_r)
      status= ga_destroy(g_tmp1)
      status= ga_destroy(g_tmp2)
      if(.not.status) call errquit
     $   ('davidson pop: could not destroy global arrays', 0, GA_ERR)

      return
      end


      subroutine grid_make_sa(rtdb)
      implicit none
#include "rtdb.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
      integer nbf_ao2
      integer nconstr, g_constr(max_constr)
      integer iatms1(max_constr), fatms1(max_constr) 
      integer iatms2(max_constr), fatms2(max_constr)
      integer isa, lsa, isa2, lsa2, i

      nbf_ao2=nbf_ao*nbf_ao

      if (.not. rtdb_get(rtdb, 'dft:nconstr', mt_int, 1, nconstr))
     $  call errquit('grid_make_sa get nconstr failed',0, RTDB_ERR)

      if (.not. (
     $          rtdb_get(rtdb, 'dft:iatms1', mt_int, nconstr, iatms1(1))
     $    .and. rtdb_get(rtdb, 'dft:fatms1', mt_int, nconstr, fatms1(1))
     $    .and. rtdb_get(rtdb, 'dft:iatms2', mt_int, nconstr, iatms2(1))
     $    .and. rtdb_get(rtdb, 'dft:fatms2', mt_int, nconstr, fatms2(1))
     $    ) )
     &  call errquit('grid_make_sa: atm not found', 0, RTDB_ERR)
      if (.not. rtdb_get(rtdb,'dft:gconstr',mt_int,nconstr,g_constr(1)))
     $  call errquit('grid_make_sa get g_constr failed',0, RTDB_ERR)

      if (.not.ma_push_get(mt_dbl,nbf_ao2,'sa',lsa,isa))
     &   call errquit('frag_dens: cannot allocate sa',0, MA_ERR)
      do i = 1, nconstr
       call dfill(nbf_ao2,0.d0,dbl_mb(isa),1)
       call frag_dens(iatms1(i),fatms1(i),dbl_mb(isa))
       if(iatms2(i).ne.0) then
         if (.not.ma_push_get(mt_dbl,nbf_ao2,'sa2',lsa2,isa2))
     &    call errquit('frag_dens: cannot allocate sa2',0, MA_ERR)
         call dfill(nbf_ao2,0.d0,dbl_mb(isa2),1)
         call frag_dens(iatms2(i),fatms2(i),dbl_mb(isa2))
         call daxpy(nbf_ao2,-1.d0,dbl_mb(isa2),1,dbl_mb(isa),1)
         if(.not.MA_pop_stack(lsa2))
     &      call errquit("grid_make_sa: pop sa2", 100, MA_ERR)
       endif
       if(ga_nodeid().eq.0)
     G  call ga_put(g_constr(i),1,nbf_ao,1,nbf_ao,dbl_mb(isa),nbf_ao)
       call ga_sync()
      enddo
      if(.not.MA_pop_stack(lsa))
     &     call errquit("grid_make_sa: pop sa", 100, MA_ERR)
      return
      end


      subroutine frag_dens(iatm,fatm,sa)
      implicit none
c Qin / partly copied grid/grid_cckohn.F
c    / to find out the total density on the grid points centered from iatm to fatm
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "cdft.fh"
#include "cgridfile.fh"
#include "grid_cube.fh"
#include "geom.fh"
      integer iatm, fatm
      double precision sa(nbf_ao,nbf_ao)

      integer iqsh, istep, nxyz, ncontrset
      integer ixyz, lxyz, icharge, lcharge, itags, ltags
      integer lrqbuf,irqbuf,lqxyz,iqxyz,lqwght,iqwght, nqpts,ncube,
     &     ictr_buf

      double precision rad
      integer lbas_cset_info, ibas_cset_info,
     &     lbas_cent_info, ibas_cent_info,
     &     ldocset, idocset,
     &     l_rchi_atom,i_rchi_atom,
     &     l_rq,i_rq,lniz, iniz,
     &     lchi_ao, ichi_ao 
      integer i
      logical grid_file_rewind
      external grid_file_rewind
c
c     ao basis set info used by xc_eval_basis
c
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('frag_dens:bas_numcont',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*ncenters, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('frag_dens: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('frag_dens: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), ncenters)
c
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('frag_dens: cannot allocate ccdocset',
     .     ncontrset, MA_ERR)
      do i=1,ncontrset
         log_mb(idocset+i-1)=.true.
      enddo
c
      if(.not.MA_push_get(MT_int, ncenters, 'iniz',
     &     lniz, iniz))
     &     call errquit("frag_dens:iniz",0, MA_ERR)
      do i= 1, ncenters
         int_mb(iniz+i-1)=1
      enddo
c
      nxyz = 3*ncenters
      if (.not.MA_push_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('frag_dens: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,ncenters,'charge',lcharge,icharge))
     &   call errquit('frag_dens: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,ncenters*16,'tags',ltags,itags))
     &   call errquit('frag_dens: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('frag_dens: geom_cart_get failed',74, GEOM_ERR)

      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('frag_dens: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('frag_dens: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4,
     &     'quad pts buffer', lrqbuf, irqbuf))
     &     call errquit('frag_dens: quad buffer', 3, MA_ERR)

      if (.not. grid_file_rewind())
     $     call errquit('frag_dens: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
c     loop over records in the grid file.
c
      do 200 iqsh = 1, n_rec_in_file
c
c     Define the current range of radial shells and integration center.
c
         call grid_file_read(n_per_rec, nqpts, ictr_buf,
     &        rad,dbl_mb(irqbuf),nsubb)

         if(ictr_buf.lt.iatm .or. ictr_buf.gt.fatm) goto 200

         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts, UNKNOWN_ERR)
         if(nqpts.eq.0) goto 200
         istep=0
c
c     loop over the cubes?
c
         do  ncube=1,nsubb
c
c     temp put buf into currently used arrays qxyz and qwght
c
            call grid_repack(dbl_mb(irqbuf), dbl_mb(iqxyz),
     &           dbl_mb(iqwght), nqpts, rad,istep)
            if(nqpts.ne.0) then
c
c     compute the basis functions over the grid
c
               if(.not.MA_Push_get(MT_dbl, ncenters, 'rchi_atom',
     &              l_rchi_atom,i_rchi_atom))
     &              call errquit("frag_dens:rchi_atom",0, MA_ERR)
               if(.not.MA_Push_get(MT_dbl, nqpts*ncenters, 'rq',
     &              l_rq,i_rq))
     &              call errquit("frag_dens:rq",0, MA_ERR)
               if (.not.MA_Push_Get(mt_dbl, nqpts*nbf_ao,
     &              'chi_ao', lchi_ao, ichi_ao))
     &              call errquit('frag_dens: chi_ao',0, MA_ERR)
c
               call qdist(dbl_mb(i_rchi_atom), dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters)
               call xc_eval_basis(ao_bas_han, 0, dbl_mb(ichi_ao),
     &              0d0, 0d0, 0d0, dbl_mb(i_rq),
     &              dbl_mb(iqxyz), dbl_mb(ixyz), nqpts, ncenters,
     &              int_mb(iniz), log_mb(idocset),
     &              int_mb(ibas_cent_info), int_mb(ibas_cset_info))
               call acc_sa(dbl_mb(ichi_ao), dbl_mb(iqwght), sa,
     $                      nbf_ao, nqpts)

               if(.not.MA_pop_stack(lchi_ao))
     &              call errquit("frag_dens: pop chi_ao", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rq))
     &              call errquit("frag_dens: pop rq", 100, MA_ERR)
               if(.not.MA_pop_stack(l_rchi_atom))
     &             call errquit("frag_dens: pop rchi_atom",100,MA_ERR)
            endif
         enddo
 200  continue
      if(.not.MA_pop_stack(lrqbuf))
     &     call errquit("frag_dens: pop rqbuf", 100, MA_ERR)
      if(.not.MA_pop_stack(lqwght))
     &     call errquit("frag_dens: pop qwght", 100, MA_ERR)
      if(.not.MA_pop_stack(lqxyz))
     &     call errquit("frag_dens: pop qxyz", 100, MA_ERR)
      if(.not.MA_pop_stack(ltags))
     &     call errquit("frag_dens: pop tags", 100, MA_ERR)
      if(.not.MA_pop_stack(lcharge))
     &     call errquit("frag_dens: pop charge", 100, MA_ERR)
      if(.not.MA_pop_stack(lxyz))
     &     call errquit("frag_dens: pop xyz", 100, MA_ERR)
      if(.not.MA_pop_stack(lniz))
     &     call errquit("frag_dens: pop niz", 100, MA_ERR)
      if(.not.MA_pop_stack(ldocset))
     &     call errquit("frag_dens: pop docset", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cset_info))
     &     call errquit("frag_dens: pop bas_cset_info", 100, MA_ERR)
      if(.not.MA_pop_stack(lbas_cent_info))
     &     call errquit("frag_dens: pop bas_cent_info", 100, MA_ERR)

      return
      end


      subroutine acc_sa(chi_ao,qwght,sa,nbf,nqpts)
      implicit none
      integer nbf, nqpts
      double precision chi_ao(nqpts,nbf),qwght(nqpts),sa(nbf,nbf)
      integer i, j, k
      do i = 1, nbf
        do j = 1, nbf
          do k = 1, nqpts
            sa(i,j) = sa(i,j) + chi_ao(k,i)*qwght(k)*chi_ao(k,j)
          enddo
        enddo
      enddo
      return
      end


      subroutine cdft_deriv(nconstr, g_constr,
     $                    g_movecs, g_dens, evals,
     $                    noc, nmo, grad, hess,
     $                    g_tmp, fddl)
      implicit none
      integer nconstr 
c      double precision constr(nconstr)
      integer g_constr(nconstr), g_movecs, g_dens, g_tmp
      double precision evals(*)
      integer noc, nmo
      double precision grad(nconstr),hess(nconstr,nconstr)
      double precision fddl
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      integer lgtia, igtia, ngtia, iigtia, jjgtia
      integer ioc, iuoc, nuoc, i, j
      double precision temp
      double precision zero, onem, one
      double precision delta_e,eps
      parameter(zero = 0.d0, one = 1.d0, onem = -one, eps=1d-9)

      integer g_t
      logical status

      nuoc = nmo - noc
      do i = 1, nconstr
        grad(i) = grad(i) + ga_ddot(g_dens,g_constr(i)) 
      enddo

      ngtia=noc*nuoc*nconstr
      if(.not.MA_Push_Get(MT_Dbl,ngtia,'gtia',lgtia,igtia))
     $       call errquit('derivatives: cannot allocate gtia)',0,MA_ERR)
c
      call dfill(ngtia, 0.0d0, dbl_mb(igtia), 1)
c
      status=ga_create(mt_dbl, nmo, nmo, 'ttmp', 0, nmo, g_t)
      do i = 1, nconstr
        call ga_dgemm('t', 'n', nmo, nmo,
     $      nmo, one, g_movecs, g_constr(i), zero, g_tmp)
        call ga_dgemm('n', 'n', nmo, nmo,
     $      nmo, one, g_tmp, g_movecs, zero, g_t)
        do ioc = 1+ga_nodeid(), noc,ga_nnodes()
#if 0
           do iuoc = noc+1, nmo
              call ga_get(g_t,ioc,ioc,iuoc,iuoc,temp,1)
              j = (iuoc-noc) + (ioc-1)*nuoc
     $                + (i-1)*noc*nuoc
              dbl_mb(igtia+j-1) = temp
              enddo
#else
              j = igtia+
     P              (ioc-1)*nuoc
     $                + (i-1)*noc*nuoc
              call ga_get(g_t,ioc,ioc,noc+1,nmo,dbl_mb(j),1)
#endif
        enddo
      enddo
      status=ga_destroy(g_t)
      call ga_dgop(2015,dbl_mb(igtia),ngtia,'+')

      do i = 1, nconstr
        do ioc = 1, noc
           do 10 iuoc = noc+1, nmo
              iigtia = (iuoc-noc) + (ioc-1)*nuoc
     $                       +(i-1)*noc*nuoc
              delta_e=evals(ioc)-evals(iuoc)
              if(abs(delta_e).lt.eps) then
                 if(ga_nodeid().eq.0) then
                    write(6,*) ' degenerate evals ',ioc,iuoc,
     C                   evals(ioc),evals(iuoc),dbl_mb(igtia+iigtia-1)
                 endif
                 temp=0d0
              else
                 temp=dbl_mb(igtia+iigtia-1)/(evals(ioc)-evals(iuoc))
              endif
              do j = 1, nconstr
                jjgtia = (iuoc-noc) + (ioc-1)*nuoc
     $                       +(j-1)*noc*nuoc
                hess(i,j) = hess(i,j) + fddl*
     $                      temp*dbl_mb(igtia+jjgtia-1)
              enddo
10         continue
        enddo
      enddo

      if(.not.ma_pop_stack(lgtia))
     $   call errquit('cdft_deriv: cannot pop stack',0,MA_ERR)

      return
      end


      subroutine find_la(lambda, dl, ddl, do_ddl, dlold, x1, x2,
     $                    lbrak,check)
      implicit none

      double precision lambda, dl, ddl, dlold, x1, x2
      logical do_ddl, lbrak, check

      double precision zero, golden
      parameter(zero=0.d0, golden=1.6d0)

      check = .false.

      if(dlold.eq.zero) then
         lbrak = .true.
         x1 = zero
         dlold = dl
         x2 = -dl/ddl
         if(abs(x2).gt.0.1d0) x2 = sign(0.1d0,x2)
         lambda = x2
         do_ddl = .false.
         return
      elseif(lbrak .and. dlold*dl.gt.zero) then
         if(abs(dlold).lt.abs(dl)) then
           lambda=x2
           x2 = x1 + golden*(x1-x2)
           x1=lambda
           dlold = dl
         else
           x2 = x2 + golden*(x2-x1)
         endif
         lambda = x2
         do_ddl = .false.
         return
      elseif(lbrak .and. dlold*dl.lt.zero) then
         lambda = (x1+x2)/2.d0
         lbrak = .false.
         do_ddl = .true.
         return
      endif

      if(dl*dlold.lt.zero) then
         x2 = lambda
      else
         x1 = lambda
         dlold = dl
      endif

c!! Here try Newton step.  If it's got out of [x1, x2], then bisection.

      lambda = lambda - dl/ddl
      if((lambda-x1)*(lambda-x2).gt.zero) 
     $    lambda = (x1+x2)/2.d0
      if(abs(x1-x2).lt.1d-10) then
         check = .true.
      endif

      return
      end


      subroutine cdft_newt(nconstr, multipliers, grad, hess, grad_old,
     $                      p, check,u,v,w,work)
      implicit none
      integer nconstr
      double precision multipliers(nconstr), grad(nconstr),
     $                 hess(nconstr, nconstr), grad_old(nconstr)
      double precision p(nconstr),u(nconstr,nconstr),
     $                 w(nconstr),v(nconstr,nconstr), 
     $                 work(5*nconstr), wmin, s
      logical backtrack, check
      integer i, j, info
      double precision x0, x1, back, y1, a, b, c
      double precision thres
      parameter (thres=1d-8)

      check = .false.
      backtrack = .false.

      x0=0d0
      x1=0d0
      do i = 1, nconstr
         x0=x0+grad(i)*grad(i)
         x1=x1+grad_old(i)*grad_old(i)
      enddo
      if((x0-x1).gt.1d-8) backtrack = .true.

      if(backtrack) then
        x0=0d0
        x1=0d0
        y1=0d0
        do i = 1, nconstr
           x0=x0+grad_old(i)*p(i)
           x1=x1+grad(i)*p(i)
           do j = 1 , nconstr
             y1=y1+p(i)*hess(i,j)*p(j)
           enddo
        enddo
        a=x0-x1+y1
        b=2d0*(x1-x0)-y1
        c=b*b-4d0*a*x0
        if(c.lt.0) then
          back = x0/(x0-x1)
        else
          back = (-b-sqrt(c))/(2d0*a)
        endif
        if(back.gt.5d-1) back = 5d-1
        if(back.lt.1d-1) back = 1d-1
        do i = 1, nconstr
          multipliers(i) = multipliers(i) - (1.d0-back)*p(i)
          p(i) = p(i)*back
          if(dabs(p(i)).lt.1d-10) check = .true.
        enddo
        return
      endif

      do i = 1, nconstr
        grad_old(i) = grad(i)
        p(i) = -grad(i)
      enddo

c use svd in case there are redundant constraints
      do i = 1, nconstr
        p(i) = 0.d0
      enddo
      call dgesvd('A','A',nconstr,nconstr,hess,nconstr,
     $            w,u,nconstr,v,nconstr,work,5*nconstr,info)
c
      wmin = w(1) * thres
      do 60 j=1, nconstr
        if(w(j).le.wmin)GO TO 60
        s = 0.0d0
        do 40 i=1,nconstr
           s = s - u(i,j) * grad(i)  ! minus grad is used
   40   continue
        s = s/w(j)
        do 50 i=1,nconstr
           p(i)=p(i)+s*v(j,i)
   50   continue
   60 continue

      do i = 1, nconstr
        multipliers(i) = multipliers(i) + p(i)
      enddo

      return
      end


      subroutine grad_cdens_lowd(rtdb, g_dens, g_cdens)
      implicit none
      integer rtdb
      integer g_dens(2)              ! input
      integer g_cdens             ! input/output

#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "util.fh"
#include "tcgmsg.fh"
#include "schwarz.fh"
#include "msgids.fh"
c#include "dftps.fh"

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_s, g_svecs, g_shalf
      integer g_t1, g_t2, g_t3
      integer isvals, lsvals
      integer nconstr, ic
      integer iatms1(max_constr), iatms2(max_constr)
      integer fatms1(max_constr), fatms2(max_constr)
      double precision multipliers(max_constr)
      integer dtype(max_constr)

      logical status
      integer i, j, ifirst(2), ilast(2)
      double precision fij
      double precision toll_s
      parameter (toll_s=1d-10)

      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call ga_zero(g_s)
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals,
     &   isvals))
     &   call errquit('grad_cdft: cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
      call int_1e_ga(AO_bas_han, AO_bas_han, g_s, 'overlap', oskel)
      if (oskel)call sym_symmetrize(geom, AO_bas_han, .false., g_s)
c
c     Diag Overlap Matrix
c
      g_svecs = ga_create_atom_blocked(geom, AO_bas_han,
     &                                 'AO ovl eig vecs')
      call ga_zero(g_svecs)
#if defined(PARALLEL_DIAG)
      call ga_diag_std(g_s, g_svecs, Dbl_MB(isvals))
#else
      call ga_diag_std_seq(g_s, g_svecs, Dbl_MB(isvals))
#endif
c
      g_shalf = ga_create_atom_blocked(geom, AO_bas_han, 'shalf')
      call ga_zero(g_shalf)
      g_t1 = ga_create_atom_blocked(geom, AO_bas_han, 'temp1')
      g_t2 = ga_create_atom_blocked(geom, AO_bas_han, 'temp2')
      g_t3 = ga_create_atom_blocked(geom, AO_bas_han, 'temp3')
      call diis_bld12(toll_s, dbl_mb(isvals), g_svecs, g_shalf,
     &    g_t1,  3)           !!!! get S**1/2 in g_shalf

       if (.not. rtdb_get(rtdb, 'dft:nconstr', mt_int, 1, nconstr))
     &   call errquit('grid_cden_lowd: get nconstr failed', 0, RTDB_ERR)
      if (.not. (
     $          rtdb_get(rtdb,'dft:iatms1', mt_int, nconstr,iatms1(1))
     $    .and. rtdb_get(rtdb,'dft:iatms2', mt_int, nconstr,iatms2(1))
     $    .and. rtdb_get(rtdb,'dft:fatms1', mt_int, nconstr,fatms1(1))
     $    .and. rtdb_get(rtdb,'dft:fatms2', mt_int, nconstr,fatms2(1))
     $    ))
     & call errquit('grad_cden_lowd: atms failed', 2700, RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'dft:multipliers', mt_dbl, nconstr,
     $     multipliers(1)))
     & call errquit('grad_cden_lowd: get multipliers failed',0,RTDB_ERR)
      if (.not. rtdb_get(rtdb,'dft:dtype', mt_int, nconstr, dtype(1)))
     &   call errquit('grad_cden_lowd: get dtype', 3000, RTDB_ERR)

      do ic = 1, nconstr

      call ga_copy(g_svecs, g_t2)
      call ga_copy(g_shalf, g_t3)

      status = bas_ce2bfr(ao_bas_han, iatms1(ic), ifirst(1), i)
      status = bas_ce2bfr(ao_bas_han, fatms1(ic), i, ilast(1))
      if(iatms2(ic).ne.0) then
       status = bas_ce2bfr(ao_bas_han, iatms2(ic), ifirst(2), i)
       status = bas_ce2bfr(ao_bas_han, fatms2(ic), i, ilast(2))
      endif

      multipliers(ic) = -multipliers(ic) * 2d0

c build  t1_li = \sum_v P_lv * c_vi
      if(dtype(ic).eq.1) then
       call ga_dgemm('n', 'n', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     $                      g_dens(1), g_t2, 0d0, g_t1)
      else
       call ga_dgemm('n', 'n', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     $                      g_dens(2), g_t2, 0d0, g_t1)
      endif

c build  s_lj = \sum_u s^1/2_lu * c_uj
      call ga_matmul_patch('t', 'n', 1.d0, 0.d0,
     $                      g_t3, 1, nbf_ao, ifirst(1), ilast(1),
     $                      g_t2, ifirst(1), ilast(1), 1, nbf_ao,
     $                      g_s, 1, nbf_ao, 1, nbf_ao)
      if(iatms2(ic).ne.0) then
        call ga_matmul_patch('t', 'n', -1.d0, 1.d0,
     $                      g_t3, 1, nbf_ao, ifirst(2), ilast(2),
     $                      g_t2, ifirst(2), ilast(2), 1, nbf_ao,
     $                      g_s, 1, nbf_ao, 1, nbf_ao)
      endif

c build  t3_ij (x_ij) = \sum_l s_lj * t1_li
      call ga_dgemm('t', 'n', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     $                      g_t1, g_s, 0d0, g_t3)
c build f_ij*x_ij
      do i = 1, nbf_ao
        do j = 1, nbf_ao
          call ga_get(g_t3,i,i,j,j,fij,1)
          fij = fij
     $          /(sqrt(dbl_mb(isvals+i-1))+sqrt(dbl_mb(isvals+j-1)))
          call ga_put(g_s,i,i,j,j,fij,1)
        enddo
      enddo
      call ga_dadd(multipliers(ic), g_s, 0.d0, g_t3, g_t3)
c
c build s_mj (g_mj)
      call ga_dgemm('n', 'n', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     $                      g_t2, g_t3, 0d0, g_s)

c finally build t_mn
      call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao, 1.d0,
     $                      g_s, g_t2, 0d0, g_t1)

      call ga_symmetrize(g_t1)

      call ga_add(1.d0,g_cdens,1d0,g_t1,g_cdens)

      enddo

      if (.not. ga_destroy(g_t1)) call errquit
     &   ('grad_cdft: could not destroy g_t1', 1, GA_ERR)
      if (.not. ga_destroy(g_t2)) call errquit
     &   ('grad_cdft: could not destroy g_t2', 1, GA_ERR)
      if (.not. ga_destroy(g_t3)) call errquit
     &   ('grad_cdft: could not destroy g_t3', 1, GA_ERR)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('grad_cdft: could not destroy g_s', 1, GA_ERR)
      if (.not. ga_destroy(g_shalf)) call errquit
     &   ('grad_cdft: could not destroy g_shalf', 1, GA_ERR)
      if (.not. ga_destroy(g_svecs)) call errquit
     &   ('grad_cdft: could not destroy g_svecs', 1, GA_ERR)

      if (.not.ma_pop_stack(lsvals))
     &     call errquit('grad_cdft: cannot pop stack', 12, MA_ERR)

      return
      end
