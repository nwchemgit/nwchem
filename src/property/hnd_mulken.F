      subroutine hnd_mulken(rtdb,basis,geom)
* $Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "util.fh"
      integer   rtdb, basis, geom
      integer   ga_create_atom_blocked
      external  ga_create_atom_blocked
      logical status
      integer g_s, g_spin, g_p, g_ps, atn
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      integer bf_offset,sh_offset
      integer natoms,nshells,nbf,lmax
      integer l_tags,k_tags,l_aogr,k_aogr,l_atov,k_atov,
     ,     l_atsh,k_atsh,l_scr,k_scr
      integer iseg,idens,i,irow,m,iatom,jatom,jlo,jhi,ii,ilo,ihi,ij
      integer ij_offset,jshell,jtype,l_b,k_b,l_v,k_v,l_f,k_f
      integer idum
      double precision dsum,dums,temp
      external dsum
      character*3 scftyp
      character symbol*2, symbolj*2, tag*16, element*20, thedens*10
      character*1 types(8)
      logical prt_ovlp_pop     ! if .true. print the overlap populations
      logical prt_ao_pop       ! if .true. print the AO gross overlap
                               ! populations
      logical prt_atm_ovlp_pop ! if .true. print the overlap populations
                               ! condensed to atoms
      data types/'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K'/
c
      if (ga_nodeid().eq.0) write(luout,2100) 
c
c     Decide what we need to print
c
      prt_ovlp_pop     = ((ga_nodeid().eq.0).and.
     &                   util_print("overlap population",print_debug))
      prt_ao_pop       = ((ga_nodeid().eq.0).and.
     &                   util_print("ao gross population",print_debug))
      prt_atm_ovlp_pop = ((ga_nodeid().eq.0).and.
     &                   util_print("atomic overlap population",
     &                   print_debug))
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Determine dimensions:
c     natoms : number of atoms
c     nshells: number of shells
c     nbf    : number of basis functions
c     lmax   : highest angular momentum
c
      if(.not.geom_ncent(geom, natoms))
     &   call errquit('hnd_mulken: geom_ncent failed',0, GEOM_ERR)
      if(.not. bas_numcont(basis,nshells))
     &   call errquit('hnd_mulken: geom_numcont failed',1, BASIS_ERR)
      if(.not. bas_numbf(basis,nbf))
     &   call errquit('hnd_mulken: geom_numbf failed',1, BASIS_ERR)
      if (.not.bas_high_angular(basis,lmax))
     &   call errquit('hnd_mulken: lmax error',lmax, INT_ERR)
c
c     Node zero collects data in MA's
c
      if (ga_nodeid().eq.0) then
         if(.not. ma_push_get(mt_byte,nbf*16,'bf tags',l_tags,k_tags))
     $    call errquit('hnd_mulken, alloc tags failed',0, MA_ERR)
c
c     Node zero keeps a list of atom/basis function tags in format
c     Atom#  Atomtag  Function (Char*16)
c
         call bas_vec_info(basis,byte_mb(k_tags))
c
         if(.not.ma_push_get(mt_dbl,nbf*ndens,'aogross',l_aogr,k_aogr))
     &      call errquit('hnd_mulken, alloc aogr failed',0, MA_ERR)
         call dcopy(nbf*ndens, 0.0d0, 0, dbl_mb(k_aogr), 1)
         if(.not.ma_push_get(mt_dbl,natoms*natoms,'atgross',l_atov,
     &      k_atov)) call errquit('hnd_mulken, atov failed',0, MA_ERR)
         if(.not.ma_push_get(mt_dbl,natoms*(lmax+2)*ndens,'atshell',
     &      l_atsh,k_atsh)) call errquit('hnd_mulken, atsh failed',0,
     &      MA_ERR)
         call dcopy(natoms*(lmax+2)*ndens, 0.0d0, 0, dbl_mb(k_atsh), 1)
         if(.not.ma_push_get(mt_dbl,nbf,'scratch',l_scr,k_scr))
     &     call errquit('hnd_mulken, alloc scr failed',911, MA_ERR)
      endif
c
c     Global array for overlap
c
      g_s = ga_create_atom_blocked(geom, basis,'hnd_mulken s')
      g_p = ga_create_atom_blocked(geom, basis,'hnd_mulken p')
c
c     Loop over densities. For RHF only one density, otherwise
c     there will be three, alpha, beta and total
c
      do idens = 1, ndens
         if (ndens.gt.1) then
            if (idens.eq.1) thedens = 'Alpha Spin'
            if (idens.eq.2) thedens = 'Beta Spin'
            if (idens.eq.3) thedens = 'Total'
         else
            thedens = 'Total'
         endif
c
c     Calculate overlap population
c
         call ga_zero(g_s)
         call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
         call ga_elem_multiply(g_s,g_dens(idens),g_s)
c
         if (ga_nodeid().gt.0) goto 1000
c
         bf_offset = k_aogr + nbf * (idens-1)
         sh_offset = k_atsh + natoms * (lmax+2) * (idens-1)
c
c     Print overlap population in blocks of seven at a time
c
         if (prt_ovlp_pop) then
           write(luout,2120) thedens
           do iseg = 1, nbf, 7
              write (luout,2128) (i,i=iseg,min(iseg+6,nbf))
              do irow = 1, nbf
                 call ga_get(g_s,irow,irow,iseg,min(iseg+6,nbf),
     &                       dbl_mb(k_scr),1)
                 write (luout,2148) irow,
     &                       (byte_mb(k_tags+(irow-1)*16+m),m=0,15),
     &                       (dbl_mb(k_scr+i),i=0,min(nbf-iseg,6))
              enddo
           enddo
         endif
c
c     Calculate gross population for ao's and atoms
c     For each row also print the gross population
c
         if (prt_ao_pop) write(luout,2298) thedens
         call dcopy(natoms*natoms, 0.0d0, 0, dbl_mb(k_atov), 1)
         do irow = 1, nbf
c
            if(.not.bas_bf2ce(basis,irow,iatom)) call
     &         errquit('hnd_mulken, bas_bf2ce failed',911, BASIS_ERR)
c
c     Get a row of the overlap population and calculate
c     the gross population for ao's and atoms
c
            call ga_get(g_s,irow,irow,1,nbf,dbl_mb(k_scr),1)
c
c     Use row to calculate the total gross population in ao's
c     (i.e. collapse whole row to single number)
c
            dbl_mb(bf_offset+irow-1) = dsum(nbf,dbl_mb(k_scr),1)
c
c     Print total gross population in ao's 
c
            if (prt_ao_pop)
     &        write(luout,2297) irow,
     &                          (byte_mb(k_tags+(irow-1)*16+m),m=0,15),
     &                          dbl_mb(bf_offset+irow-1)
c
c     Use row to calculate the total overlap population in atoms
c
            do jatom = 1, natoms
               if(.not.bas_ce2bfr(basis,jatom,jlo,jhi)) call
     &            errquit('hnd_mulken, bas_bf2ce failed',911, BASIS_ERR)
               ij_offset = k_atov+(jatom-1)*natoms+iatom-1
               dbl_mb(ij_offset)=dbl_mb(ij_offset) + 
     &                           dsum(jhi-jlo+1,dbl_mb(k_scr+jlo-1),1)
            enddo
         enddo
c
c     Print overlap population condensed to atoms
c
         if (prt_atm_ovlp_pop) then
           write(luout,2396) thedens
           do iatom = 1, natoms, 7
              write (luout,2128) (i,i=iatom,min(iatom+6,natoms))
              do jatom = 1, natoms
                 ij_offset = k_atov+(jatom-1)*natoms+iatom-1
                 if (.not. geom_cent_tag(geom, jatom, tag))
     $               call errquit('bas_vec_info: tag', jatom, BASIS_ERR)
                 status = geom_tag_to_element(tag, symbol, element, atn)
                 write (luout,2397) jatom,symbol,(dbl_mb(ij_offset+i),
     &                                        i=0,min(natoms-iatom,6))
              enddo
           enddo
         endif
c
c     Compress from ao's to shells
c     Use total gross population in ao's array dbl_mb(bf_offset)
c
         do jshell = 1, nshells
            if (.not. bas_cn2bfr(basis, jshell, jlo, jhi)) call errquit
     &         ('hnd_mulken: bas_cn2bfr failed ',basis,BASIS_ERR)
            if (.not.bas_cn2ce(basis,jshell,jatom)) call
     &         errquit('hnd_mulken, bas_bf2ce failed',911, BASIS_ERR)
            if (.not.bas_continfo(basis,jshell,jtype,idum,idum,idum)) 
     &         call errquit('hnd_mulken, continfo fail',911, BASIS_ERR)
            ij_offset = sh_offset+(jatom-1)*(lmax+2)+jtype
            dbl_mb(ij_offset)=dbl_mb(ij_offset) + 
     &                        dsum(jhi-jlo+1,dbl_mb(bf_offset+jlo-1),1)
         enddo
c
c     Print shell populations per atom
c
         write(luout,2499) thedens,(types(i),i=1,lmax+1)
         write(luout,2497)
         do iatom = 1, natoms
            ij_offset = sh_offset+(iatom-1)*(lmax+2)
            if (.not. geom_cent_tag(geom, iatom, tag))
     $          call errquit('bas_vec_info: tag', iatom, BASIS_ERR)
            status = geom_tag_to_element(tag, symbol, element, atn)
            write(luout,2498) iatom,symbol,(dbl_mb(ij_offset+i-1),
     &                                      i=1,lmax+1)
         enddo
         call util_flush(luout)
c
c     Calculate total gross population on atoms
c     Use shell population in dbl_mb(sh_offset)
c     Print total gross population on atoms from dbl_mb(sh_offset)
c
         write(luout,2595) thedens
         do jatom = 1, natoms
            ij_offset = sh_offset+(jatom-1)*(lmax+2)
            if (.not. geom_cent_tag(geom, jatom, tag))
     $          call errquit('bas_vec_info: tag', jatom, BASIS_ERR)
            status = geom_tag_to_element(tag, symbol, element, atn)
            dbl_mb(ij_offset+lmax+1) = dsum(lmax+1,dbl_mb(ij_offset),1)
            write(luout,2592) jatom,symbol,dble(atn),
     &                        dbl_mb(ij_offset+lmax+1)
         enddo
c
c     Calculate orbital and atomic spin densities
c
         if (idens.eq.3) then
            if (prt_ao_pop) then
              write(luout,2694)
              do irow = 1, nbf
                 bf_offset = k_aogr + nbf 
                 write(luout,2297) irow,
     &                 (byte_mb(k_tags+(irow-1)*16+m),m=0,15),
     &                 (dbl_mb(k_aogr+irow-1)-dbl_mb(bf_offset+irow-1))
              enddo
            endif
            write(luout,2693)
            do irow = 1, natoms
               sh_offset = k_atsh+(irow-1)*(lmax+2)
               ij_offset = k_atsh+natoms*(lmax+2)+(irow-1)*(lmax+2)
               if (.not. geom_cent_tag(geom, irow, tag))
     $             call errquit('bas_vec_info: tag', irow, BASIS_ERR)
               status = geom_tag_to_element(tag, symbol, element, atn)
               write(luout,2592) irow,symbol,dble(atn),
     &               (dbl_mb(sh_offset+lmax+1)-dbl_mb(ij_offset+lmax+1))
            enddo
         endif
c
c     At this point all the atoms come back together
c
 1000 call ga_sync()
c
      enddo ! for a density block
c
c     Node zero cleans up MA's
c
      if (ga_nodeid().eq.0) then
         if(.not.ma_pop_stack(l_scr)) call errquit('hnd_mulken, 
     &      ma_pop_stack of scr failed',911,MA_ERR)
         if(.not.ma_pop_stack(l_atsh)) call errquit('hnd_mulken, 
     &      ma_pop_stack of atsh failed',911,MA_ERR)
         if(.not.ma_pop_stack(l_atov)) call errquit('hnd_mulken, 
     &      ma_pop_stack of atov failed',911,MA_ERR)
         if(.not.ma_pop_stack(l_aogr)) call errquit('hnd_mulken, 
     &      ma_pop_stack of aogr failed',911,MA_ERR)
         if(.not.ma_pop_stack(l_tags)) call errquit('hnd_mulken, 
     &      ma_pop_stack of tags failed',911,MA_ERR)
      endif
c
c     ----- bond index and valency analysis -----
c
      if (ndens.eq.3) then
         g_spin = ga_create_atom_blocked(geom, basis,'hnd_mulken spin')
         call ga_add(1.0d0,g_dens(1),-1.0d0,g_dens(2),g_spin)
         g_ps = ga_create_atom_blocked(geom, basis,'hnd_mulken ps')
      endif
c
c     ----- closed and open shell systems population analysis -----
c
c     Original reference for NDO systems:
c     D.R. Armstrong, P.G. Perkins and J.J.P. Stewart
c     J.C.S. Dalton Trans. 838 (1973)
c     Extension  first principles: I. Mayer, Chem. Phys. Lett.
c     97,270(1983).
c     Number of free electrons: O.G. Stradella, H.O. Villar, E.A. Castro
c     Theor. Chim. Acta(Berl.) 70, 67(1986)
c
c       d * s = p (non - hermitean)
c
      call ga_zero(g_s)
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
      call ga_dgemm('N','T',nbf,nbf,nbf,1.0d0,g_dens(ndens),g_s,
     &              0.0d0,g_p)
      if (ndens.eq.3) call ga_dgemm('N','T',nbf,nbf,nbf,1.0d0,g_spin,
     &                              g_s,0.0d0,g_ps)
c
c     ----- -b-, valency , and free electron -----
c
      if (.not.ma_push_get(mt_dbl,natoms*natoms,'b mat',l_b,k_b)) 
     &   call errquit('hnd_mulken, b mat failed',0, MA_ERR)
      call dcopy(natoms*natoms, 0.0d0, 0, dbl_mb(k_b), 1)
      if(.not.ma_push_get(mt_dbl,natoms,'v vec',l_v,k_v)) 
     &   call errquit('hnd_mulken, v vec failed',0, MA_ERR)
      call dcopy(natoms, 0.0d0, 0, dbl_mb(k_v), 1)
      if(.not.ma_push_get(mt_dbl,natoms,'f vec',l_f,k_f)) 
     &   call errquit('hnd_mulken, f vec failed',0, MA_ERR)
      call dcopy(natoms, 0.0d0, 0, dbl_mb(k_f), 1)
c
      do iatom = 1, natoms
         ii = (iatom-1)*natoms+iatom-1
         if (.not.bas_ce2bfr(basis,iatom,ilo,ihi)) call
     &      errquit('hnd_mulken, bas_bf2ce failed',911, BASIS_ERR)
         do jatom = 1, natoms
            ij = (iatom-1)*natoms+jatom-1
            if (.not.bas_ce2bfr(basis,jatom,jlo,jhi)) call
     &         errquit('hnd_mulken, bas_bf2ce failed',911, BASIS_ERR)
            dbl_mb(k_b+ij) = ga_ddot_patch(g_p,'N',jlo,jhi,ilo,ihi,
     &                                     g_p,'T',ilo,ihi,jlo,jhi) 
            if(ndens.eq.3) dbl_mb(k_b+ij) = dbl_mb(k_b+ij) + 
     &                           ga_ddot_patch(g_ps,'N',jlo,jhi,ilo,ihi,
     &                                         g_ps,'T',ilo,ihi,jlo,jhi)
         enddo
         dums=0.0d0
         if(ndens.eq.3) dums=ga_ddot_patch(g_ps,'N',ilo,ihi,ilo,ihi,
     &                                     g_ps,'T',ilo,ihi,ilo,ihi)
         dbl_mb(k_v+iatom-1)=-ga_ddot_patch(g_p,'N',ilo,ihi,ilo,ihi,
     &                                      g_p,'T',ilo,ihi,ilo,ihi)
         do i = ilo,ihi
            call ga_get(g_p,i,i,i,i,temp,1)
            dbl_mb(k_v+iatom-1)=dbl_mb(k_v+iatom-1)+2.0d0*temp
         enddo
         dbl_mb(k_f+iatom-1)=(dbl_mb(k_v+iatom-1)-dums+
     &                        dbl_mb(k_b+ii))/2.0d0
         do jatom=1,natoms
            if(iatom.ne.jatom) then
               ij = (iatom-1)*natoms+jatom-1
               dbl_mb(k_f+iatom-1)=dbl_mb(k_f+iatom-1)-dbl_mb(k_b+ij)
            endif
         enddo
         dbl_mb(k_b+ii)=0.0d0
      enddo
c
c     Only node 0 prints, rest skips this section
c
      if (ga_nodeid().gt.0) goto 2000
c
c     ----- print all bond indices -----
c
      write(luout,9999)
      do jatom = 1, natoms, 8
         do iatom = 1, natoms
            write(luout,9998) (iatom,i,dbl_mb(k_b+(iatom-1)*natoms+i-1),
     &                         i=jatom,min(jatom+7,natoms))
         enddo
      enddo
c
c     ----- now print only the largest one -----
c
      write(luout,9993)
      do iatom = 1, natoms
         do jatom = iatom+1, natoms
            ij = (iatom-1)*natoms+jatom-1
            if (dbl_mb(k_b+ij).gt.0.1d0) then
               if (.not. geom_cent_tag(geom, iatom, tag))
     $            call errquit('bas_vec_info: tag', iatom, BASIS_ERR)
               status = geom_tag_to_element(tag, symbol, element, atn)
               if (.not. geom_cent_tag(geom, jatom, tag))
     $            call errquit('bas_vec_info: tag', jatom, BASIS_ERR)
               status = geom_tag_to_element(tag, symbolj, element, atn)
               write(luout,9994) iatom,symbol,jatom,symbolj,
     &                           dbl_mb(k_b+ij)
            endif
         enddo
      enddo
c
c     sum bond indices for a given atom
c     v(a) = n(a) + sum of b(a,b) is the mulliken charge
c     spin(a) = sum of b(a,b)-v(a)
c     All directly calculated when printing...
c
c     ----- print analysis -----
c
      do jatom = 1, natoms
          ii = (jatom-1)*natoms+jatom-1
          dums = 0.0d0
          do iatom = 1, natoms
            ij = (iatom-1)*natoms+jatom-1
            dums = dums + dbl_mb(k_b+ij)
          enddo
          dbl_mb(k_b+ii) = dums
      enddo
c
      write(luout,9997)
      do iatom = 1, natoms
         ii = (iatom-1)*natoms+iatom-1
         if (.not. geom_cent_tag(geom, iatom, tag))
     $         call errquit('bas_vec_info: tag', iatom, BASIS_ERR)
         status = geom_tag_to_element(tag, symbol, element, atn)
         write(luout,9996) iatom,symbol,dbl_mb(k_v+iatom-1),
     &                     dbl_mb(k_f+iatom-1),dbl_mb(k_b+ii),
     &                     (dbl_mb(k_b+ii)+dbl_mb(k_f+iatom-1)),
     &                     (dbl_mb(k_v+iatom-1)-dbl_mb(k_b+ii))
      enddo

 2000 call ga_sync()
c
c     Clean up memory
c
      if(.not.ma_pop_stack(l_f)) call errquit('hnd_mulken, 
     &   ma_pop_stack of l_f failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_v)) call errquit('hnd_mulken, 
     &   ma_pop_stack of l_v failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_b)) call errquit('hnd_mulken, 
     &   ma_pop_stack of l_b failed',911,MA_ERR)
      if (ndens.eq.3) then
         if(.not.ga_destroy(g_ps)) call
     &     errquit('hnd_mulken, ga_destroy failed',911, BASIS_ERR)
         if(.not.ga_destroy(g_spin)) call
     &     errquit('hnd_mulken, ga_destroy failed',911, BASIS_ERR)
      endif
      if(.not.ga_destroy(g_p)) call
     &         errquit('hnd_mulken, ga_destroy failed',911, BASIS_ERR)
      if(.not.ga_destroy(g_s)) call
     &         errquit('hnd_mulken, ga_destroy failed',911, BASIS_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('hnd_mulken: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
      call ga_sync()
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 2100 format(//,10x,28(1h-),/,10x,'Mulliken population analysis',/,
     & 10X,28(1h-))
 2120 format(/,10x,'----- ',a10,' overlap population -----')
 2128 format(/,23x,7(6x,i3,6x),/)
 2148 format(i5,2x,16a,7f15.10)
 2298 format(/,10x,'----- ',a10,' gross population in ao -----',/)
 2297 format(i5,2x,16a,f11.5)
 2396 format(/,10x,'----- ',a10,
     &       ' overlap population condensed to atoms -----',/)
 2397 format(i5,2x,a2,14x,7f15.10)
 2499 format(/,1x,a10,' S,P,D,... shell population',/,1x,32(1h-),
     &       /,'    Atom  ',8(8x,a1))
 2497 format(1x,86(1h-))
 2498 format(1x,i5,1x,a2,2x,8f10.5)
 2595 format(/,10x,'----- ',a10,' gross population on atoms ----',/)
 2592 format(10x,i5,2x,a2,f6.1,f12.5)
 2694 format(/,10x,'----- Spin population in ao -----',/)
 2693 format(/,10x,'----- Atomic spin population -----',/)
 9999 format(/,10x,'----- Bond indices -----')
 9998 format(8(1x,i2,'-',i2,f10.5))
 9997 format(/,62x,'  Free electrons',5x,'        Valency',/,
     1         28x,'  Number of   ',3x,'   Sum of    ',
     2          6x,'+ Bond indices ',5x,'   - Bond indices    ',/,
     3         13x,'   Valency',5x,'Free electrons',4x,'Bond indices',
     4          5x,'=Mulliken charge',5x,'= Net spin population')
 9996 format(1x,i5,2x,a2,3x,f10.5,5x,f10.5,7x,f10.5,8x,f10.5,
     1      14x,f10.5)
 9995 format(/)
 9994 format(5(i3,1x,a4,'-',i3,1x,a4,f10.5))
 9993 format(/,1x,'Large bond indices',/,1x,18(1h-))
      end
