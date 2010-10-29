      subroutine datoms(geom, basis, hatom, g_dens, oprin,
     +     pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +     cvec, copn, smin, qmin, transf, cc, nbb, oprint_ener,
     $     odens, eatom,spec)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
      integer geom, basis
      logical odens,spec        ! If true make the global density
      double precision eatom    ! Return sum of atomic energies
c     
*rak: wrapper routine now.  doatoms_do does the actual work
c     
      integer nbb               ! No. of primtives used to dimension arrays
      integer g_dens
      real *8 hatom(*)
      real *8 pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      real *8 h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      real *8 cvec(*), copn(*), smin(nbb,*), qmin(nbb,*),transf(*),cc(*)
      logical oprin, oprint_ener
c
      logical status_ma
      integer nat_dim
      integer h_odone, k_odone
      integer h_tags,  k_tags
c
      if (.not. geom_ncent(geom, nat_dim)) call errquit
     $     ('guess atom dens: geom_ncent failed', geom, GEOM_ERR)
      if (.not. ma_push_get
     &    (MT_LOG,nat_dim,'atomscf:doatomd:odone',h_odone,k_odone))
     &    call errquit('ma get atomscf:doatomd:odone failed',911,
     &       MA_ERR)
      if (.not. ma_push_get
     &    (MT_BYTE,16*nat_dim,'atomscf:doatomd:tags',h_tags,k_tags))
     &    call errquit('ma get atomscf:doatomd:tags failed',911,
     &       MA_ERR)

      call datoms_do(geom, basis, hatom, g_dens, oprin,
     &     pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     &     cvec, copn, smin, qmin, transf, cc, nbb, oprint_ener,
     &     odens, eatom, byte_mb(k_tags), log_mb(k_odone), nat_dim,spec)
      status_ma = .true.
      status_ma = status_ma .and. ma_pop_stack(h_tags)
      status_ma = status_ma .and. ma_pop_stack(h_odone)
      if (.not.status_ma) call errquit
     &    ('doatoms: fatal error popping tags/odone stack',911, MEM_ERR)
      end
      subroutine datoms_do(geom, basis, hatom, g_dens, oprin,
     +     pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +     cvec, copn, smin, qmin, transf, cc, nbb, oprint_ener,
     $     odens, eatom, tags, odone, nat_dim,spec)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cguess.fh"
#include "msgids.fh"
#include "stdio.fh"
      integer geom, basis
      logical odens,spec        ! If true make the global density
      double precision eatom    ! Return sum of atomic energies
c     
c...  subroutine to coordinate atom-scf calls and d-matrix gathering
c...  for atomic startup
c...  h,t  : full h-matrix and t-matrix to supply to atom
c...  d    : full density matrix as return parameter
c...  
c...  **note** data is transferred to atom directly via common/cguess/
c     
      integer nbb               ! No. of primtives used to dimension arrays
      integer g_dens
      real *8 hatom(*)
      real *8 pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      real *8 h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      real *8 cvec(*), copn(*), smin(nbb,*), qmin(nbb,*),transf(*),cc(*)
      logical oprin, oprint_ener
c     
      integer nat_dim
      logical odone(nat_dim)
      character*16 tags(nat_dim)
c     
      double precision charge, coords(3)
      character*16 tag, element
      character*2 symbol
      integer atn
      integer natoms
      integer shift
      integer icnlo, icnhi, ibflo, ibfhi, jbfhi, jbflo, jat
c     
      integer ic(4,nb),iiloc(nb,4)
c     
      logical omyjob            ! True if this process does this atom
      integer nproc, me, icount ! Parallel processing info
c     
      integer ispdf, iat, j
      integer isymax
      real *8 toteng, znps
c     
      integer maxtyp
      parameter (maxtyp = 6)
c
      logical ospherical
c     
      toteng = 0.0d0
c     
c     Determine tags and which centers are identifiable atoms
c     with both positive charge and basis funcitons.
c     Odone keeps track of which centers have been done ... mark
c     centers that do not seem to be ataoms as done already.
c     
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('guess atom dens: geom_ncent failed', geom, GEOM_ERR)
      if (natoms .gt. nat_dim) call errquit
     $     ('guess atom dens: mxnat too small ', natoms, UNKNOWN_ERR)
c
      do iat = 1 , natoms
         if (.not. geom_cent_get(geom, iat, tag, coords, charge))
     $        call errquit('guess atom dens: geom_cent_get?', iat,
     &       GEOM_ERR)
         if (.not. bas_ce2cnr(basis, iat, icnlo, icnhi))
     $        call errquit('guess atom dens: bas_ce2cnr?', iat,
     &       BASIS_ERR)
         odone(iat) =  .not.
     $        (geom_tag_to_element(tag, symbol, element, atn)
     $        .and.  charge.ge.1.0d0 .and. icnhi.gt.0)
         tags(iat)  = tag
      enddo
c     
c     Loop thru atoms doing the calculation
c     
      me = ga_nodeid()
      nproc = ga_nnodes()
      icount = 0
c     
      do iat = 1, natoms
c     
         if (.not. odone(iat)) then
            if (.not. geom_cent_get(geom, iat, tag, coords, charge))
     $           call errquit('guess atom dens: geom_cent_get?', iat,
     &       GEOM_ERR)
c     
c     Decide if this atom is my task
c     
            omyjob = mod(icount,nproc) .eq. me
            icount = icount + 1
            if (.not. omyjob) then
               do jat = iat, natoms
                  if (tags(jat) .eq. tag) odone(jat) = .true.
               enddo
               goto 1010
            endif
c     
*       write(6,*) me,' doing ', iat, ' ', tag, element
c     
            call atscf_bas_gather(geom, basis, iat, ic, iiloc,
     $           isymax, hatom, znps, element, tag, ospherical)
c     
            call atomd(oprin,luout,znps,ic,isymax,hatom,
     +           pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +           cvec, copn, smin, qmin, transf, cc, nbb,
     +           element, tag,spec)
c
*ZIJING      call atom_eig_vec(znps,iat,cc,nsym,ncsh,nosh,nccup,nbc)
c     
c     loop thru atoms with the same tag and add in energy and density
c     
            if (.not. bas_ce2bfr(basis, iat, ibflo, ibfhi)) 
     $           call errquit('guess atom dens: bas_ce2bfr', iat,
     &       BASIS_ERR)
            do jat = iat, natoms
               if ((.not. odone(jat)) .and. tags(iat).eq.tags(jat)) then
                  toteng = toteng + energ
                  odone(jat) = .true.
c
                  if (odens) then
                     if (.not. bas_ce2bfr(basis, jat, jbflo, jbfhi))
     $                    call errquit('guess atom dens:bas_ce2bfr',jat,
     &       BASIS_ERR)
c     
c     shift adjusts iiloc() from atom iat to atom jat
c     
                     shift = jbflo - ibflo
                     do ispdf = 1, 4
                        do j = 1, nbc(ispdf)
                           iiloc(j,ispdf) = iiloc(j,ispdf) + shift
                        enddo
                     enddo
c     
                     call creded(g_dens, dt, iiloc, nbb, ospherical)
c     
c     undo the shifting ready for next atom ... clearer than fudging shift
c     
                     do ispdf = 1, 4
                        do j = 1, nbc(ispdf)
                           iiloc(j,ispdf) = iiloc(j,ispdf) - shift
                        enddo
                     enddo
                  endif
               endif
            enddo
         endif
c     
 1010    continue               ! Jump here if not this node task
      enddo
c     
      call ga_sync()
c     
      call ga_dgop(Msg_AtSCF, toteng, 1, '+')
      if (ga_nodeid().eq.0 .and. oprint_ener) then
         write (luout,6010) toteng
         call util_flush(luout)
      endif
c
      eatom = toteng
c
      return
 6010 format(' Sum of atomic energies: ',f20.8)
      end
      subroutine atscf_bas_gather(geom, basis, iat, ic, iiloc,
     $     isymax, hatom, znps, element, tag, ospherical)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "bas.fh"
#include "cguess.fh"
#include "geom.fh"
      integer geom              ! [input]
      integer basis             ! [input]
      integer iat               ! [input]
      integer ic(4,nb)          ! [input]
      integer iiloc(nb,4)       ! [output]
      integer isymax            ! [output]
      double precision hatom(*) ! [output]
      double precision znps     ! [output]
      character*16 tag, element ! [output]
      logical ospherical        ! [output]
c     
c     Gather all of the basis/atomic information necessary in order
c     to call atomd() to perform an atomic SCF calculation on atom
c     number iat.
c     
      integer nshell
      double precision charge, coords(3)
      character*2 symbol
      integer atn
      integer type, nprim, ngen, sphcart, ispdf
      integer nbf_per_shell, ibf, igen, i, j, k, kh, kkzc
      double precision ee, fac
      integer icnlo, icnhi, ibflo, ibfhi, icont
      integer iisch(nb,4),iicnt(nb,4)
      logical status
      double precision expnt(nw_max_prim)
      double precision coeff(nw_max_prim*nw_max_gen_con)
      double precision pi32
c*debug-s
*      integer di,dj
*      integer mycalls
*      save mycalls
*      data mycalls /0/
c*debug-e
      data pi32/5.56832799683170d0/
c*debug-s
*      mycalls = mycalls + 1
*      write(6,*)'bas_gat:i:',mycalls,'  geom  : ',geom
*      write(6,*)'bas_gat:i:',mycalls,'  basis : ',basis
*      write(6,*)'bas_gat:i:',mycalls,'  iat   : ',iat
*      write(6,*)'bas_gat:i:',mycalls,'  ic    : ',
*     &    ((ic(di,dj),di=1,4),dj=1,nb)
c*debug-s
c     
c     gather shell / symmetry info
c     
      if (.not. bas_ce2cnr(basis, iat, icnlo, icnhi))
     $     call errquit('guess atom dens: bas_ce2cnr?', iat,
     &       BASIS_ERR)
      do i = 1 , 4
         nbc(i) = 0
      enddo
c     
c     nbc  # shells / symmetry
c     iisch  contains index of shell
c     iiloc  contains position of starting ao of shell in "real" world
c     translate to 1 (s)
c     
      nshell = 0
c     
c     loop thru general contractions on atom
c     
      do icont = icnlo, icnhi
         if (.not. bas_continfo(basis, icont, type, nprim, ngen,
     $        sphcart))
     $        call errquit('guess atom dens: bas_continfo', icont,
     &       BASIS_ERR)
*debug:              write(6,*)' basis handle    =',basis
*debug:              write(6,*)' for contraction =',icont
*debug:              write(6,*)' type            =',type
*debug:              write(6,*)' nprim           =',nprim
*debug:              write(6,*)' ngen            =',ngen
*debug:              write(6,*)' sphcart         =',sphcart
         if (.not. bas_cn2bfr(basis, icont, ibflo, ibfhi))
     $        call errquit('guess atom dens: bas_cn2bfr', icont,
     &       BASIS_ERR)
         if (nprim.gt.nw_max_prim .or. ngen.gt.nw_max_gen_con) 
     $       call errquit ('guess atom dens: nprim/ngen too big', 0,
     &       UNKNOWN_ERR)
c     
         if (type .le. 3) then  ! Only s/p/d/f functions
            nbf_per_shell = (ibfhi - ibflo + 1) / ngen
            if (type .eq. -1) nbf_per_shell = 1
            ibf = ibflo
c     
c     loop thru shells inside contractions
c     
            do igen = 1, ngen
               nshell = nshell + 1
c     
               if (type .eq. -1) then
                  ispdf = igen  ! 1/2 = s/p
               else
                  ispdf = type + 1
               endif
c     
               nbc(ispdf) = nbc(ispdf) + 1
               iicnt(nbc(ispdf),ispdf) = icont
               iisch(nbc(ispdf),ispdf) = igen
               iiloc(nbc(ispdf),ispdf) = ibf
c     
               ibf = ibf + nbf_per_shell
            enddo
         endif
      enddo
      ospherical = sphcart .eq. 1
c..   
c..   we gathered symmetry/shell info ; now get the real thing
c..   
      kkzc = 0
      kh = 0
      isymax = 0
      do ispdf = 1 , 4
c..   nbas = total # primitives for this symmetry
         nbas(ispdf) = 0
         if (nbc(ispdf).gt.0) isymax = ispdf
         do j = 1 , nbc(ispdf)
c     
c     get contraction coeffs and exponents
c     
            icont = iicnt(j,ispdf)
            igen  = iisch(j,ispdf)
            if (.not. bas_get_coeff(basis, icont, coeff))
     $           call errquit('guess atom dens: get_coeff', icont,
     &       BASIS_ERR)
            if (.not. bas_get_exponent(basis, icont, expnt))
     $           call errquit('guess atom dens: get_expnt', icont,
     &       BASIS_ERR)
            if (.not. bas_continfo(basis, icont, type, nprim,
     $           ngen, sphcart))
     $           call errquit('guess atom dens: continfo', icont,
     &       BASIS_ERR)
            
c..   ic = # number of primitives /contracted /symmetry
            ic(ispdf,j) = nprim
            nbas(ispdf) = nbas(ispdf) + nprim
c..   gather the primitives
            if(kkzc+nprim.gt.nb) call errquit(
     C           ' cguess.fh: increase nb to',kkzc+nprim,0)
            do k = 1, nprim
               kkzc = kkzc + 1
               zeta(kkzc) = expnt(k)
               cont(kkzc) = coeff(k + nprim*(igen-1))
c...  get contraction coeffs as we are used to
               ee = 2*zeta(kkzc)
               fac = pi32/(ee*sqrt(ee))
               if (ispdf.eq.2) then
                  fac = 0.5d0*fac/ee
               else if (ispdf.eq.3) then
*     fac = (1.0d0/3.0d0)  *  0.75d0*fac/(ee*ee)
                  fac = 0.75d0*fac/(ee*ee)
               else if (ispdf.eq.4) then
                  fac = 1.875d0*fac/(ee**3)
               end if
               cont(kkzc) = cont(kkzc)*sqrt(fac)
               ucont(kkzc) = 1.0d0/sqrt(fac)
            enddo
c...  in the pseudopotential case, we would be involved in
c...  gathering  the h integrals for the contracted aos
c...  so that they are added in at the right time (in atomd)
c...  use proper offset to use pure d or f functions (ioffhp)
c...  only the comments remain ...
            do k = 1 , j
               kh = kh + 1
               hatom(kh) = 0.0d0
            enddo
c...  
         enddo
      enddo
c..   
c..   all prepared call  atomd
c..   zeta,cont,nbas,nbc,nbas,ic,zn are passed via cguess
c..   energ and the density matrix dt are received via cguess
c..   note zn is the real nuclear charge / znps is the effective charg
c..   
      status = geom_cent_get(geom, iat, tag, coords, charge)
      status = geom_tag_to_element(tag, symbol, element, atn)
      zn = dble(atn)
      znps = charge
c
c call ecp integrals to fill hatom properly
c
      call atscf_ecp(geom,basis,tag,hatom,nbas,
     &    zeta,cont,ucont,kkzc)
*debug-s
*      write(6,*)'bas_gat:o:',mycalls,'  ic    : ',
*     &    ((ic(di,dj),di=1,4),dj=1,nb)
*      write(6,*)'bas_gat:o:',mycalls,'  iiloc  : ',
*     &    ((iiloc(di,dj),di=1,nb),dj=1,4)
*      write(6,*)'bas_gat:o:',mycalls,'  isymax : ',isymax
*      write(6,*)'bas_gat:o:',mycalls,'  znps   : ',znps
*      write(6,*)'bas_gat:o:',mycalls,'  tag    : ',tag
*      write(6,*)'bas_gat:o:',mycalls,'  element: ',element
*      write(6,*)'bas_gat:o:',mycalls,'  ospherical: ',ospherical
*debug-e
c     
      end
