      Subroutine fukui(geom, nbf_ao, AO_bas_han, g_movecs,
     &                 noc, k_eval,  ipol, oskel, tol2e)
c     
c     $Id: fukui.F,v 1.5 1999-10-13 03:43:53 d3j191 Exp $
c     
      implicit none
      integer g_dens_HOMO(2), g_dens_LUMO(2), g_orb, g_dens_ss, ik,
     &        g_movecs(2), noc(2), isp, k_eval(2), g_s, geom, nbf_ao,
     &        ipol, AO_bas_han, me, l_temp_vec, i_temp_vec
      double precision eig_lumo(2), eig_homo(2), jfac(4), kfac(4),
     &                 ecoul_orb, gap_spin, tol2e,
     &                 mu_n_mas, mu_n_men, mu_n_cer,
     &                 mu_s_mas, mu_s_men, mu_s_cer
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oprint_fukui, oskel
      me = ga_nodeid()
      oprint_fukui = util_print('Fukui information', print_high)
      g_orb = ga_create_atom_blocked(geom, AO_bas_han, 'ga_orb')
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'temp vec',
     &   l_temp_vec, i_temp_vec))
     &   call errquit('fukui: cannot allocate temp vec',0)
      do isp = 1, 2
         g_dens_HOMO(isp) = ga_create_atom_blocked(geom,
     &                      AO_bas_han, 'ga_dens_orb_homo')
         g_dens_LUMO(isp) = ga_create_atom_blocked(geom,
     &                      AO_bas_han, 'ga_dens_orb_lumo')
         call ga_zero(g_dens_HOMO(isp))
         call ga_zero(g_dens_LUMO(isp))
      end do
      g_dens_ss = ga_create_atom_blocked(geom, AO_bas_han,
     &            'ga_dens_orb_ss')
      do isp = 1, ipol
         do ik = noc(isp), (noc(isp)+1)
            call ga_zero(g_orb)
            call ga_get(g_movecs(isp), 1, nbf_ao, ik, ik,
     &                  Dbl_MB(i_temp_vec), nbf_ao)
            call ga_put(g_orb, 1, nbf_ao, ik, ik, 
     &                  Dbl_MB(i_temp_vec), nbf_ao)
            if (ik.eq.noc(isp)) then
               eig_homo(isp) = dbl_mb(k_eval(isp) + 
     &                         noc(isp) - 1)
               call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao,
     &                       1.0d00, g_orb, g_orb, 0.d00, 
     &                       g_dens_HOMO(isp))
            else
               eig_lumo(isp) = dbl_mb(k_eval(isp) + noc(isp))
               call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao,
     &                       1.0d00, g_orb, g_orb, 0.d00,
     &                       g_dens_LUMO(isp))
            end if
         end do
      end do
c
      if (.not.ma_pop_stack(l_temp_vec))
     &   call errquit('fukui: cannot pop stack',0)
c      
      if (ipol.eq.1) then
         call ga_copy(g_dens_HOMO(1), g_dens_HOMO(2))
         call ga_copy(g_dens_LUMO(1), g_dens_LUMO(2))
         eig_homo(2) = eig_homo(1)
         eig_lumo(2) = eig_lumo(1)
      end if
c
      if (me.eq.0.and.oprint_fukui)
     &    call dft_header(' Reactivity Parameters ')
c
      mu_n_mas = 0.5d00*(eig_lumo(1) + eig_lumo(2))
      mu_n_men = 0.5d00*(eig_homo(1) + eig_homo(2))
      mu_n_cer = 0.5*(mu_n_mas + mu_n_men)
      mu_s_mas = 0.5d00*(eig_lumo(1) - eig_homo(2))
      mu_s_men = 0.5d00*(eig_homo(1) - eig_lumo(2))
      mu_s_cer = 0.5*(mu_s_mas + mu_s_men)
c
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call ga_zero(g_s)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_LUMO(2), 0.5d00,
     &            g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header('  Condensed Fukui function [fnn(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_HOMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fnn(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(-0.5d00, g_dens_LUMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fsn(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(-0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_HOMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fsn(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fss(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(-0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fns(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(1), 0.5d00,
     &             g_dens_LUMO(2), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fss(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(1), -0.5d00,
     &             g_dens_LUMO(2), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fns(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s, 'total')
c
      mu_n_men = 0.5d00*(eig_homo(1) + eig_homo(2))
      mu_n_cer = 0.5*(mu_n_mas + mu_n_men)
      mu_s_mas = 0.5d00*(eig_lumo(1) - eig_homo(2))
      mu_s_men = 0.5d00*(eig_homo(1) - eig_lumo(2))
      mu_s_cer = 0.5*(mu_s_mas + mu_s_men)
      if (me.eq.0.and.oprint_fukui)then
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,*) '    mu_n(+)    mu_n(-)    mu_n(0)'
         write(LuOut,'(3f11.4)')  mu_n_mas, mu_n_men, mu_n_cer
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,*) '    mu_s(+)    mu_s(-)    mu_s(0)'
         write(LuOut,'(3f11.4)')  mu_s_mas, mu_s_men,mu_s_cer
         write(LuOut,*) ' ------------------------------------'
      endif
      call ga_zero(g_dens_ss)
      call ga_add(-1.0d00,g_dens_HOMO(2),1.0d00,
     &             g_dens_LUMO(1),g_dens_ss)
      kfac(1) = 0.d00
      jfac(1) = 1.0d0
      jfac(2) = 1.0d0
      kfac(2) = 0d0
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_ss, g_orb)
      ecoul_orb = 0.5d0*ga_ddot(g_dens_ss,g_orb)
      gap_spin = eig_lumo(1)-eig_homo(2) + ecoul_orb
      if (me.eq.0.and.oprint_fukui)then
        write(LuOut,*) ' ------------------------------------'
        write(LuOut,*) ' Vertical Singlet-Triplet splitting energy'
        write(LuOut,*) '  delta_eig    integral        GAP'
        write(LuOut,'(3f11.4)') eig_lumo(1)-eig_homo(2),ecoul_orb,
     &                          gap_spin  
        write(LuOut,*) ' ------------------------------------'
      end if
      if (.not. ga_destroy(g_orb)) call errquit
     &   ('fukui: could not destroy g_orb', 0)
      if (.not. ga_destroy(g_dens_HOMO(1))) call errquit
     &   ('fukui: could not destroy g_dens_HOMO', 0)
      if (.not. ga_destroy(g_dens_LUMO(1))) call errquit
     &   ('fukui: could not destroy g_dens_LUMO', 0)
      if (.not. ga_destroy(g_dens_HOMO(2))) call errquit
     &   ('fukui: could not destroy g_dens_HOMO_bet', 0)
      if (.not. ga_destroy(g_dens_LUMO(2))) call errquit
     &   ('fukui: could not destroy g_dens_LUMO_bet', 0)
      if (.not. ga_destroy(g_dens_ss)) call errquit
     &   ('fukui: could not destroy g_dens_ss', 0)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('fukui: could not destroy g_s', 0)
      return
      end

      Subroutine mull_pop_fuk( geom, basis, iga_dens,iga_s, string)

C$Id: fukui.F,v 1.5 1999-10-13 03:43:53 d3j191 Exp $
      Implicit none
      integer geom,basis
      integer iga_s              ! overlap   GA handle
      integer iga_dens           ! dens. mat GA handle
      integer iga_ps             ! product   GA handle
      character*(*) string      ! Identifier for ecce output
                                ! total/open/alpha/beta/mp2/ccsd/...

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
     &     call errquit(' exiting in mull_pop',0)
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit(' exiting in mull_pop',1)
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in mull_op',1)
      max_at_bf2 = 0
      do iatom = 1, natoms
        if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $       call errquit('mul_pop: bas_ce2bfr failed', iatom)
        max_at_bf2 = max(max_at_bf2, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf2*max_at_bf2

      if(.not.MA_Push_Get(mt_dbl,max_at_bf2,'PS',lPSmat,iPSmat))
     &     call errquit(' exiting in mull_pop: insuff stack',21)
      if(.not.MA_Push_Get(mt_dbl,natoms,'q atom',lqatom,iqatom))
     &     call errquit(' exiting in mull_pop: insuff stack',22)
      if(.not.MA_Push_Get(mt_dbl,nshells,'q shell',lqshell,iqshell))
     &     call errquit(' exiting in mull_pop: insuff stack',3)

      iga_PS=ga_create_atom_blocked(geom, basis, 'PS product')

      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_dens,iga_s,0.d0,iga_PS)
      call mull_calc_fuk(basis,natoms, nshells,max_at_bf2,iga_PS,
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat), string)
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35)
      status= ga_destroy(iga_PS)
c
      if (oscfps) call pstat_off(ps_mull)
c
      return
      end
c
c
c
      Subroutine mull_calc_fuk(basis, natoms, nshells,max_at_bf2,iga_PS,
     &     qatom,qshell,PSmat, string)

      Implicit none
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer basis
      integer natoms,nshells
      integer iga_PS             ! product   GA handle
      integer ifirst,ilast,nbf_at,max_at_bf2
      integer ish1,ish2,ish,nn,iat,mu
      integer me,nproc, geom
      double precision psmu, coord(3), qnuc
      double precision qatom(natoms),qshell(nshells),PSmat(max_at_bf2)
      character*2 symbol
      character*16 tag
      character*32 element
      character*(*) string
      integer atn
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)

      if (.not. bas_geom(basis, geom)) call errquit
     $     ('mull_pop: bas_geom failed',basis)

      do iat=me+1,natoms,nproc
        if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &       call errquit(' exiting in mull_pop',4)
        call get_atom_block(iga_PS, basis,
     $       iat, iat, PSmat, nbf_at, nbf_at)
        mu=0
        do ish=ish1,ish2
          if (.not. bas_cn2bfr(basis,ish,ifirst,ilast))
     &         call errquit(' exiting in mull_pop.',5)
          do nn=ifirst,ilast
            mu=mu+1
            psmu=PSmat((mu-1)*nbf_at+mu)
            qshell(ish)=qshell(ish)+psmu
          enddo
          qatom(iat)=qatom(iat)+qshell(ish)
        enddo
      enddo
      call ga_sync
      call ga_dgop(Msg_Mull1,qatom,natoms,'+')
      call ga_dgop(Msg_Mull2,qshell,nshells,'+')
      if(me.eq.0) then
        write(LuOut,1)
    1   format(/'    Atom       Condensed Fukui ')
        write(luout,11)
 11     format( ' -----------   ----------------')
        do iat=1,natoms
          if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &         call errquit(' exiting in mull_pop',4)
c
          if (.not. geom_cent_get(geom, iat, tag, coord, qnuc)) 
     $         call errquit('mull_pop: geom_cent_tag failed',0)
c
          if (.not. geom_tag_to_element(tag, symbol, element, atn))
     $         symbol = 'X'
          if (ish2.ge.ish1) then
             write(LuOut,2) iat,symbol,nint(qnuc),qatom(iat) 
 2           format(1x,i4,1x,a2,i4,1x,f10.4)
          endif
        enddo
      endif
c
      call ga_sync

      return
      end
