      subroutine datoms(geom, basis, hatom, g_dens, oprin,
     + pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     + cvec, copn, smin, qmin, transf, cc, nbb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "cguess.fh"
#include "msgids.fh"
      integer geom, basis
c
c...   subroutine to coordinate atom-scf calls and d-matrix gathering
c...   for atomic startup
c...   h,t  : full h-matrix and t-matrix to supply to atom
c...   d    : full density matrix as return parameter
c...
c...   **note** data is transferred to atom directly via common/cguess/
c
      integer nbb               ! No. of primtives used to dimension arrays
      integer g_dens
      real *8 hatom(*)
      real *8 pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      real *8 h(*), dc(*), dos(*), dt(*), dold(*), ss(*)
      real *8 cvec(*), copn(*), smin(nbb,*), qmin(nbb,*),transf(*),cc(*)
      logical oprin, status
      integer iwrite
c
c     These should be dynamically allocated
c
      integer mxnat, mxnprim, mxngen
      parameter (mxnat = 500, mxnprim = 20, mxngen = 10)
      logical odone(mxnat)
      character*16 tags(mxnat)
      double precision expnt(mxnprim)
      double precision coeff(mxnprim*mxngen)
c
      double precision charge, coords(3)
      character*16 tag, element
      character*2 symbol
      integer atn
      integer natoms
      integer nbf_per_shell, ibf, igen, shift
      integer icnlo, icnhi, ibflo, ibfhi, icont, jbfhi, jbflo, jat
c
      integer ic(4,nb),iiloc(nb,4),iisch(nb,4),iicnt(nb,4)
c
      logical omyjob            ! True if this process does this atom
      integer nproc, me, icount ! Parallel processing info
c
      integer i, ispdf, iat, j, k, nshell
      integer type, nprim, ngen, sphcart
      integer kkzc, kh, isymax
      real *8 pi32, toteng, ee, fac, znps
c
c..    we need xy for d and xyz for f in gathering h-ints
c..    so set proper offsets in ioffhp (see do 140)
c
      integer ioffhp(4)
      data ioffhp/0,0,3,9/
      data iwrite/6/
c
      integer maxtyp
      parameter (maxtyp = 6)
      integer minf(maxtyp),maxf(maxtyp)
      data minf  / 1, 2,  5, 11, 21, 1 /
      data maxf  / 1, 4, 10, 20, 35, 4 /
      data pi32/5.56832799683170d0/
c
      toteng = 0.0d0
c
c     Determine tags and which centers are identifiable atoms
c     with both positive charge and basis funcitons.
c     Odone keeps track of which centers have been done ... mark
c     centers that do not seem to be ataoms as done already.
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('guess atom dens: geom_ncent failed', geom)
      if (natoms .gt. mxnat) call errquit
     $     ('guess atom dens: mxnat too small ', natoms)
      do iat = 1 , natoms
         if (.not. geom_cent_get(geom, iat, tag, coords, charge))
     $        call errquit('guess atom dens: geom_cent_get?', iat)
         if (.not. bas_ce2cnr(basis, iat, icnlo, icnhi))
     $        call errquit('guess atom dens: bas_ce2cnr?', iat)
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
     $           call errquit('guess atom dens: geom_cent_get?', iat)
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
*            write(6,*) me,' doing ', iat, ' ', tag
c
c...  gather  shell / symmetry info
c
            if (.not. bas_ce2cnr(basis, iat, icnlo, icnhi))
     $           call errquit('guess atom dens: bas_ce2cnr?', iat)
            do i = 1 , 4
               nbc(i) = 0
            enddo
c
c     nbc  # shell's / symmetry
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
     $              sphcart))
     $              call errquit('guess atom dens: bas_continfo', icont)
               if (.not. bas_cn2bfr(basis, icont, ibflo, ibfhi))
     $              call errquit('guess atom dens: bas_cn2bfr', icont)
               if (nprim.gt.mxnprim .or. ngen.gt.mxngen) call errquit
     $              ('guess atom dens: nprim/ngen too big', 0)
               if (sphcart .ne. 0) call errquit
     $              ('guess atoms dens: cannot do sphericals', 0)
c
               if (type .le. 3) then ! Only s/p/d/f functions
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
                        ispdf = igen ! 1/2 = s/p
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
     $                 call errquit('guess atom dens: get_coeff', icont)
                  if (.not. bas_get_exponent(basis, icont, expnt))
     $                 call errquit('guess atom dens: get_expnt', icont)
                  if (.not. bas_continfo(basis, icont, type, nprim,
     $                 ngen, sphcart))
     $                 call errquit('guess atom dens: continfo', icont)
                  
c..   ic = # number of primitives /contracted /symmetry
                  ic(ispdf,j) = nprim
                  nbas(ispdf) = nbas(ispdf) + nprim
c..   gather the primitives
                  do k = 1, nprim
                     kkzc = kkzc + 1
                     zeta(kkzc) = expnt(k)
                     cont(kkzc) = coeff(k + nprim*(igen-1))
c...  get contraction coeff's as we are used to
                     ee = 2*zeta(kkzc)
                     fac = pi32/(ee*sqrt(ee))
                     if (ispdf.eq.2) then
                        fac = 0.5d0*fac/ee
                     else if (ispdf.eq.3) then
                        fac = (1.0d0/3.0d0)  *  0.75d0*fac/(ee*ee)
                     else if (ispdf.eq.4) then
                        fac = 1.875d0*fac/(ee**3)
                     end if
                     cont(kkzc) = cont(kkzc)*sqrt(fac)
                  enddo
c...  in the pseudopotential case, we would be involved in
c...  gathering  the h integrals for the contracted ao's
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
            zn = charge
            znps = charge
c     
            status = geom_cent_get(geom, iat, tag, coords, charge)
            status = geom_tag_to_element(tag, symbol, element, atn)
            call atomd(oprin,iwrite,znps,ic,isymax,hatom,
     +           pcap, qcap, fc, fo, s, u, t, h, dc, dos, dt, dold, ss,
     +           cvec, copn, smin, qmin, transf, cc, nbb,
     +           element, tag)
c..   
c
c     loop thru atoms with the same tag and add in energy and density
c
            if (.not. bas_ce2bfr(basis, iat, ibflo, ibfhi)) 
     $           call errquit('guess atom dens: bas_ce2bfr', iat)
            do jat = iat, natoms
               if ((.not. odone(jat)) .and. tags(iat).eq.tags(jat)) then
                  if (.not. bas_ce2bfr(basis, jat, jbflo, jbfhi))
     $                 call errquit('guess atom dens: bas_ce2bfr', jat)
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
                  odone(jat) = .true.
                  toteng = toteng + energ
                  call creded(g_dens, dt, iiloc, nbb)
c
c     undo the shifting ready for next atom ... clearer than fudging shift
c
                  do ispdf = 1, 4
                     do j = 1, nbc(ispdf)
                        iiloc(j,ispdf) = iiloc(j,ispdf) - shift
                     enddo
                  enddo
               endif
            enddo
         endif
c
 1010    continue               ! Jump here if not this node's task
      enddo
c
      call ga_sync()
c
      call ga_dgop(Msg_AtSCF, toteng, 1, '+')
      if (ga_nodeid().eq.0) then
         call util_print_centered(6,
     $        'Superposition of Atomic Density Guess', 20, .true.)
         write (iwrite,6010) toteng
         call util_flush(iwrite)
      endif
c
      return
 6010 format(' Sum of atomic energies: ',f20.8)
      end
