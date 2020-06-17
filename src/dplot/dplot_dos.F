      subroutine dplot_dosinit(rtdb,filename)
C $Id$
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer rtdb ![in]
      character*(*) filename ! [in]
      double precision traslx,sigma ! [in]
      integer nstep,change ! [in]
      double precision e1,e2 ![in]
      double precision scalex,scaley
      double precision vde
      logical do_gaussbr
c
      logical movecs_read,movecs_read_header,int_normalize
      external movecs_read,movecs_read_header,int_normalize
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*255 title,basis_name,scftype
      integer nbf 
      integer geom,basis
      integer nsets             ! Returns no. of functions in each set
      integer ldnmo             ! Inputs size of nmo
      integer nmo(2)        
      integer iset
      integer g_vecs(2)
      integer l_eval,k_eval,l_focc,k_focc,
     K     k_proj,l_proj,k_c,l_c,k_sc,l_sc,k_fproj,l_fproj
      integer l_dos,k_dos
      double precision span,trasly,emin,emax
      character*256 dfile
      integer nproj,angproj
      logical doproj
c
c     geom 1st
      if (.not.geom_create(geom,'geometry'))
     &     call errquit('dplotdos: geom_create failed',0,0)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &     call errquit('dplotdos: geom_rtdb_load failed',0,0)
c     basis 2nd
      if (.not.bas_create(basis,'ao basis'))
     &   call errquit('dplotdos: bas_create failed',0,0)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &   call errquit('dplotdos: bas_rtdb_load failed',0,0)
      if (.not.int_normalize(rtdb,basis))
     &   call errquit('dplotdos: int_normalize failed',0,0)
      if (.not.bas_numbf(basis,nbf))
     &   call errquit('dplotdos: bas_numbf failed',0,0)
c     get evals and evecs
      ldnmo=2
      call util_file_name_resolve(filename,.false.)
      if(.not.movecs_read_header(filename(1:inp_strlen(filename)), 
     A     title, basis_name, 
     $     scftype, nbf, nsets, nmo, ldnmo))
     E     call errquit(filename(1:inp_strlen(filename))//
     F     ': mreadh failed ',0,0)
      g_vecs(1) = ga_create_atom_blocked(geom, basis, 'Scr')
      if(nsets.eq.2) 
     G     g_vecs(2) = ga_create_atom_blocked(geom, basis, 'scr2')
c
c     alloc focc,eval,evecs
c
      if(.not.MA_Push_Get(MT_dbl,nsets*nbf,'evals',l_eval,k_eval))
     A     call errquit(' pushget failed ',0,0)
      if(.not.MA_Push_Get(MT_dbl,nsets*nbf,'foccs',l_focc,k_focc))
     A     call errquit(' pushget failed ',0,0)
      do iset=1,nsets
         if(.not.movecs_read(filename(1:inp_strlen(filename)), 
     I        iset, 
     F        dbl_mb(k_focc+(iset-1)*nbf), 
     E        dbl_mb(k_eval+(iset-1)*nbf), 
     G        g_vecs(iset))) call
     A        errquit(' movecsread failed ',0,0)
      enddo
      change=0
      e1=0d0
      e2=0d0
      traslx=0d0
      if (.not.rtdb_get(rtdb,'dplot:dos_nstep',mt_int,
     &   1,nstep)) nstep=1000
      if (.not.rtdb_get(rtdb,'dplot:dos_sigma',mt_dbl,
     &   1,sigma)) sigma=0.025d0
      if (.not.rtdb_get(rtdb,'dplot:dos_span',mt_dbl,
     &   1,span)) span=5d0
      if (.not.rtdb_get(rtdb,'dplot:dos_emin',mt_dbl,
     &   1,emin)) emin=0d0
      if (.not.rtdb_get(rtdb,'dplot:dos_emax',mt_dbl,
     &   1,emax)) emax=0d0
      if (.not.rtdb_get(rtdb,'dplot:dos_doproj',mt_log,
     &   1,doproj)) doproj=.false.
      if (.not.rtdb_get(rtdb,'dplot:dos_trasly',mt_dbl,
     &   1,trasly)) trasly=0d0
      if (.not.rtdb_get(rtdb,'dplot:dos_traslx',mt_dbl,
     &   1,traslx)) traslx=0d0
      if (.not.rtdb_get(rtdb,'dplot:dos_scalex',mt_dbl,
     &   1,scalex)) scalex=1d0
      if (.not.rtdb_get(rtdb,'dplot:dos_scaley',mt_dbl,
     &   1,scaley)) scaley=1d0
      if (.not.rtdb_get(rtdb,'dplot:dos_vde',mt_dbl,
     &   1,vde)) vde=-999d0
      if (.not.rtdb_get(rtdb,'dplot:do_gaussbr',mt_log,
     &   1,do_gaussbr)) do_gaussbr=.false.
      If (.not.rtdb_cget(rtdb,'dplot:File_Out',1,dfile))
     &   dfile='fort.16'
      call util_file_name_resolve(dfile,.false.)
      If (ga_nodeid().eq.0 ) Then
         open(unit=16,file=dfile(1:inp_strlen(dfile)),
     &        status='unknown', access='sequential',form='formatted',
     &           err=111)
      End If
      if(ga_nodeid().eq.0) write(0,*) ' DOS output goes to file ',
     A     dfile(1:inp_strlen(dfile))
      if(ga_nodeid().eq.0) write(0,*) ' dopro ',doproj
      k_proj=0
      k_fproj=0
      if(doproj) then
         if (.not.rtdb_get(rtdb,'dplot:angproj',mt_int,
     &        1,angproj)) angproj=-1
         if(ga_nodeid().eq.0) write(0,*) ' angproj ',angproj
         if(angproj.ne.-1) then
            call int_init(rtdb,1,basis)
         if(.not.MA_Push_Get(MT_int,nbf,'aproj',l_proj,k_proj))
     A     call errquit(' pushget failed ',0,0)
c
c     set all projections true
c         
         do iset=1,nbf
            int_mb(k_proj+iset-1)=1
         enddo
            call dplot_getaproj(geom,basis,nbf,angproj,
     A        nproj,int_mb(k_proj))
         if(.not.MA_Push_Get(MT_dbl,nbf,'fproj',l_fproj,k_fproj))
     A     call errquit(' pushget failed ',0,0)
         if(.not.MA_Push_Get(MT_dbl,nbf,'ccc',l_c,k_c))
     A     call errquit(' pushget failed ',0,0)
         if(.not.MA_Push_Get(MT_dbl,nbf,'sscc',l_sc,k_sc))
     A     call errquit(' pushget failed ',0,0)
c     neeed to generalize to uhf
         call dplot_pmull(geom,basis,g_vecs(1),nbf,
     p     dbl_mb(k_fproj),nproj,int_mb(k_proj),
     c     dbl_mb(k_c),dbl_mb(k_sc))
            call int_terminate
         if(.not.ma_chop_stack(l_c)) call errquit
     D     ('dplotdos: machopstack failed ',0,0)
         else
            doproj=.false.
         endif
      endif
      if(.not.MA_Push_Get(MT_dbl,2*nstep,'evals',l_dos,k_dos))
     A     call errquit(' pushget failed ',0,0)
      call dplot_dos(nbf,nmo,span,traslx,nstep,sigma,
     C     scalex,scaley,vde,do_gaussbr,dbl_mb(k_dos),
     E     emin,emax,
     P     doproj,nproj,int_mb(k_proj),dbl_mb(k_fproj),trasly,
     C     change,e1,e2,nsets,
     A     dbl_mb(k_focc),dbl_mb(k_eval))
      if(.not.ma_chop_stack(l_eval)) call errquit
     D     ('dplotdos: machopstack failed ',0,0)
      if(.not.ga_destroy(g_vecs(1))) call errquit
     A ('dplotdos: gadestr failed ',0,0)

      if(nsets.eq.2) then
      if(.not.ga_destroy(g_vecs(2))) call errquit
     A ('dplotdos: gadestr failed ',0,0)
      endif
      if(.not.geom_destroy(geom)) call errquit
     A ('dplotdos: geomdestr failed ',0,0)
      if(.not.bas_destroy(basis)) call errquit
     A ('dplotdos: basdestr failed ',0,0)
      return
 111  call errquit(' failed to open file '//dfile(1:inp_strlen(dfile)),
     C   0,0)
      end
      subroutine dplot_dos(nbf,nmo,span,traslx,nstep,sigma,
     S     scalex,scaley,vde,do_gaussbr,dos,
     E     emin,emax,
     P     doproj,nproj,shproj,fproj,trasly,
     C     change,e1,e2,ipol,
     A     foccs,evals)
      implicit none
#include "global.fh"
#include "stdio.fh"
      integer nbf ![in]
      integer nmo(*) ![in]
      logical do_gaussbr ! [in]
      double precision traslx,sigma,trasly ! [in]
      double precision scalex,scaley,vde
      double precision span ![in]
      integer nstep,change ! [in]
      integer ipol ![in]
      double precision e1,e2 ![in]
      logical doproj
      integer nproj,shproj(*)
      double precision fproj(*)
      double precision evals(nbf,2),
     * foccs(nbf,2),dos(nstep,*)
      double precision emin,emax ! [in]
c
      integer i1,i2 ![in] 1,nmo?
      integer nocc(2)
      integer spin,i
      double precision eps
      parameter(eps=1d-12)
      double precision focc,eval
      integer ien
      double precision de,ei
      double precision sqpi,factpi
      double precision au2ev,efermi,fval
      double precision tinteg,pinteg,adde,renrm
      double precision nelec 
      integer stepsz,ien_0,ien_1
      parameter(au2ev=27.21d0)
c
c     traslx=-[E(N)-E(N-1)-e_HOMO]
c
      sigma=sigma/au2ev*0.5d0
      factpi=0.5d0/acos(-1d0)
      sqpi=sigma*sqrt(acos(-1d0))
      nelec=0d0
      efermi=0d0
      do spin=1,ipol
c
         i1=1
         i2=nmo(spin)
        do i=i1,i2
          if(foccs(i,spin).ge.eps) nelec=nelec+foccs(i,spin)
          if(i.lt.i2) then
             if(foccs(i,spin).gt.eps.and.foccs(i+1,spin).lt.eps) then
                 nocc(spin)=i
                 efermi=evals(i,spin)
            if(ga_nodeid().eq.0) write(luout,*)' homo eval no. ',
     N                nocc(spin),' energy (eV)=',efermi*au2ev

                endif
c     HOMO
          endif
       enddo
      enddo
c shift to match vde
      if(abs(traslx).lt.1d-9.and.vde.ne.-999d0) then
c     e_homo+traslx=vde
         traslx=-efermi*au2ev-vde
      endif
      if(abs(emin).lt.1d-5.and.abs(emax).lt.1d-5) then
c     use span to define energy range
         emin=efermi-span/au2ev
         emax=efermi+span/au2ev
      else
c     user defined emin emax in eV
         emin=(emin-traslx)/au2ev
         emax=(emax-traslx)/au2ev
      endif
      if(ga_nodeid().eq.0) write(luout,*) ' au extrm ',emin,emax
      if(ga_nodeid().eq.0) write(luout,*) ' ev extrm ',emin*au2ev,
     E     emax*au2ev
      if(ga_nodeid().eq.0) write(luout,*) ' sigma ',sigma,
     E sigma*au2ev*2d0
      if(ga_nodeid().eq.0) write(luout,*) ' scalex ',scalex
      if(ga_nodeid().eq.0) write(luout,*) ' vde ',vde
      if(ga_nodeid().eq.0) write(luout,*) ' traslx trasly ',traslx,
     E     trasly,change
      if(change.eq.1) then
         emax=(-e1+traslx)/au2ev
         emin=(-e2+traslx)/au2ev
         de=(emax-emin)/(nstep-1)
      else
         de=(emax-emin)/(nstep-1)
      endif
      tinteg=0d0
      pinteg=0d0
      do spin=1,ipol
        do ien=1,nstep
          dos(ien,spin)=0d0
        enddo
        do i=i1,i2
           eval=evals(i,spin)
           focc=foccs(i,spin)
          ei=emin
          fval=1d0
          if(doproj) fval=fproj(i)
          do ien=0,nstep-1
c
c gaussian
c
             if(do_gaussbr) then
                adde=exp(-((eval-ei)/sigma)**2)/sigma/sqrt(acos(-1d0))
c
c lorentzian broadening
c
             else
                adde=sigma/(sigma*sigma+(eval-ei)**2)*factpi
             endif
           if(foccs(i,spin).gt.eps) then
              tinteg=tinteg+adde
              pinteg=pinteg+adde*fval
           endif
            dos(ien+1,spin)=dos(ien+1,spin)+fval*adde
            ei=ei+de
          enddo
        enddo
      enddo
      ei=emin
      if(ga_nodeid().eq.0) then
      renrm=dble(nelec)/tinteg
      if(ga_nodeid().eq.0) then
         write(luout,'(A,F14.2)')' no. electrons in energy range ',nelec
         write(luout,'(A,F14.2)')' integrated PDOS ',pinteg*renrm
      endif
      if(scalex.lt.0d0) then
         stepsz=-1
         ien_0=nstep
         ien_1=1
      else
         stepsz=1
         ien_0=1
         ien_1=nstep
      endif
      do ien=ien_0,ien_1,stepsz
c        write(6,11) ei,dos(ien+1,1)*fac,dos(ien+1,2)*fac
c        write(6,11) ei*au2ev+trasl,(dos(ien+1,1)+dos(ien+1,2))*fac
        ei=emin+(ien-1)*de
         if(ipol.eq.1) then
            write(16,11) 
     X           scalex*(ei*au2ev+traslx),
     Y           scaley*(dos(ien+1,1)*renrm+trasly)
         else
            write(16,11) 
     X           scalex*(ei*au2ev+traslx),
     Y           scaley*(dos(ien+1,1)*renrm+trasly),
     Y           scaley*(dos(ien+1,2)*renrm+trasly)
         endif
   11   format(3(2x,1pe14.6))
cold        ei=ei+de
      enddo
      endif
      call ga_sync()
      
      end
      subroutine dplot_pmull(geom,basis,g_evecs,nbf,
     p     projsh,nproj,shproj,
     c     c,sc)
      implicit none
#include "global.fh"
      integer geom
      integer basis,g_evecs
      double precision projsh(*)
      integer nbf
      integer nproj
      integer shproj(*)
      double precision c(nbf),sc(nbf)
c
      integer i,j
      integer mu
c
c     forms p_shell^i (i-th MO)=
c     \sum_mu part of shell \sum_nu C_mu_i C_nu_i S_mu_nu
c     = \sum'_mu C_mu_i (\sum_nu S_mu_nu C_nu_i)
c     = \sum'_mu C_mu_i SC_mu_i
c
      integer g_s,g_sc
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c     get mulliken
c     get S
      g_s = ga_create_atom_blocked(geom, basis, 'AO ovl')
      call ga_zero(g_s)
      g_sc = ga_create_atom_blocked(geom, basis, 'AO ovl')
      call int_1e_ga(basis, basis, g_s, 'overlap', .false.)
c     get SC
      call ga_dgemm('n', 'n', nbf, nbf, nbf,1d0,g_s,
     $     g_evecs, 0d0, g_sc)
      call dcopy(nbf,0d0,0,projsh,1)
      do i=1+ga_nodeid(),nbf,ga_nnodes()
         if(shproj(i).ne.0) then
            mu=i
c     get  mu-th row
         call ga_get(g_evecs, mu,mu,1, nbf, c,1)
         call ga_get(g_sc, mu,mu,1, nbf, sc,1)
         do j=1,nbf
            projsh(j)=projsh(j)+c(j)*sc(j)
         enddo
         endif
      enddo
      call ga_sync()
      call ga_dgop(1984,projsh,nbf,'+')
      if(.not.ga_destroy(g_s)) call errquit
     A ('dplotdos: gadestr failed ',0,0)
      if(.not.ga_destroy(g_sc)) call errquit
     A ('dplotdos: gadestr failed ',0,0)
      return
      end
      subroutine dplot_getaproj(geom,basis,nbf,angproj,
     A        nproj,proj)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"
      integer geom
      integer basis
      integer nbf
      integer angproj
      integer nproj
      integer proj(*)
c     
      integer natoms,iat,ic1,ic2,icset,ifirst,ilast,
     ,     l,nprim,dum,ifun
c     
      if(.not.geom_ncent(geom,natoms)) call errquit(
     g     ' geomncent ',0,0)
      nproj=0
      do iat=1,natoms
        if (.not.bas_ce2cnr(basis,iat,ic1,ic2))
     &       call errquit('Exiting in xc_signf.',11, BASIS_ERR)
       do icset = ic1,ic2
c
c       get info about current contraction set
c
        if (.not. bas_cn2bfr(basis, icset, ifirst, ilast))
     &       call errquit('Exiting in grid_signf.',4, BASIS_ERR)
            if (.not. bas_continfo(basis, icset,  l ,nprim,
     &                             dum, dum))
     &       call errquit('Exiting in grid_signf.',4, BASIS_ERR)

            if(l.ne.angproj) then
               do ifun=ifirst,ilast
                  nproj=nproj+1
                  proj(ifun)=0
               enddo
            endif
         enddo
      enddo
c
      return
      end
