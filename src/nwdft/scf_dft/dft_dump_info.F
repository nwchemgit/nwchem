      subroutine dft_dump_info(me)
c
      implicit none
c
      integer i, me
c      
#include "cdft.fh"
c
      write(*,*)' me = ',me
c
c     Logicals
c
      write(*,*)' LOGICALS '
      write(*,*)' nodamping = ', nodamping
      write(*,*)' nolevelshifting = ', nolevelshifting
      write(*,*)' nodiis = ', nodiis
      write(*,*)' ltight = ', ltight
      write(*,*)' lmaxov = ', lmaxov
      write(*,*)' incore = ', incore
      write(*,*)' ldelley = ', ldelley
      write(*,*)' store_wght = ', store_wght
      write(*,*)' XCFIT = ', XCFIT
      write(*,*)' CDFIT = ', CDFIT
      write(*,*)' levelshift = ', levelshift
      write(*,*)' damp = ', damp
      write(*,*)' diis = ', diis
      write(*,*)' direct = ', direct
      write(*,*)' oskel = ', oskel
      write(*,*)' oadapt = ', oadapt
      write(*,*)' lcfac(i) = ', (lcfac(i),i=1,numfunc)
      write(*,*)' nlcfac(i) = ', (nlcfac(i),i=1,numfunc)
      write(*,*)' lxfac(i) = ', (lxfac(i),i=1,numfunc)
      write(*,*)' nlxfac(i) = ', (nlxfac(i),i=1,numfunc)
      write(*,*)' xccomb(i) = ', (xccomb(i),i=1,numfunc)
c
c     Integers
c
      write(*,*)' INTEGERS '
      write(*,*)' iAOacc = ', iAOacc
      write(*,*)' iCDacc = ', iXCacc
      write(*,*)' iXCacc = ', iXCacc
      write(*,*)' itol2e = ', itol2e
      write(*,*)' imull = ', imull
      write(*,*)' iangquad = ', iangquad
      write(*,*)' iterations = ', iterations
      write(*,*)' iguess = ', iguess
      write(*,*)' mult = ', mult
      write(*,*)' ipol = ', ipol
      write(*,*)' nfock = ', nfock
      write(*,*)' ncydp = ', ncydp
      write(*,*)' ncysh = ', ncysh
      write(*,*)' ncyds = ', ncyds
      write(*,*)' noio = ', noio
      write(*,*)' ndamp = ', ndamp
      write(*,*)' idecomp = ', idecomp
      write(*,*)' nagrid = ', nagrid
      write(*,*)' nquad_task = ', nquad_task
      write(*,*)' nrinc = ', nrinc
      write(*,*)' nrmax = ', nrmax
      write(*,*)' geom = ', geom
      write(*,*)' ncenters = ', ncenters
      write(*,*)' nbf = ', nbf
      write(*,*)' nradpts = ', nradpts
      write(*,*)' nang_leb_pts = ', nang_leb_pts
      write(*,*)' nbf_ao = ', nbf_ao
      write(*,*)' nbf_cd = ', nbf_cd
      write(*,*)' nbf_xc = ', nbf_xc
      write(*,*)' MXqshells = ', MXqshells
      write(*,*)' MXrad = ', MXrad
      write(*,*)' MXagrid = ', MXagrid
      write(*,*)' MXtheta = ', MXtheta
      write(*,*)' MXphi = ', MXphi
      write(*,*)' IAC_MAX = ', IAC_MAX
      write(*,*)' nqshells = ', nqshells
      write(*,*)' iqrot = ', iqrot
      write(*,*)' ntypes = ', ntypes
      write(*,*)' nquad_task = ', nquad_task
      write(*,*)' dft_ntags_bsmx = ', dft_ntags_bsmx
      write(*,*)' dft_max_cent = ', dft_max_cent
      write(*,*)' AO_bas_han = ', AO_bas_han
      write(*,*)' CD_bas_han = ', CD_bas_han
      write(*,*)' XC_bas_han = ', XC_bas_han
      write(*,*)' nrad: ',(nrad(i), i = 1, dft_ntags_bsmx)
      write(*,*)' nang: ',(nang(i), i = 1, dft_ntags_bsmx)
c
c     Doubles
c
      write(*,*)' DOUBLES '
      write(*,*)' rlshift = ', rlshift
      write(*,*)' dampon = ', dampon
      write(*,*)' dampoff = ', dampoff
      write(*,*)' levlon = ', levlon
      write(*,*)' levloff = ', levloff
      write(*,*)' diison = ', diison
      write(*,*)' diisoff = ', diisoff
      write(*,*)' rcharge = ', rcharge
      write(*,*)' tol_rho = ', tol_rho
      write(*,*)' cfac = ', (cfac(i),i=1,numfunc)
      write(*,*)' xfac =', (xfac(i),i=1,numfunc)
      write(*,*)' hl_tol = ', hl_tol
      write(*,*)' e_conv = ', e_conv
      write(*,*)' d_conv = ', d_conv
      write(*,*)' g_conv = ', g_conv
c
c     Characters
c
      write(*,*)' CHARACTERS '
      write(*,*)' title = ', title
      write(*,*)' gridtype = ', gridtype
      write(*,*)' xname = ', xname
      write(*,*)' cname = ', cname
c
      return
      end
c $Id$
