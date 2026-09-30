c     DFTD4 interface developed with assistance from OpenAI Codex.
c     Implementation and numerical validation were reviewed and tested by MPH.
c     Validatation performed via tests varying processor count, S66 benchmarks
c     against native DFTD4 command line interface with multiple functionals
c     optimization + frequency tests against model systems (see dftd4 tests
c     in QA/tests) against references (inquire MPH for details). 
c     
      subroutine xc_dftd4(rtdb,geom,value,gradient,what)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
      integer rtdb,geom
      double precision value,gradient(*)
      character*(*) what
      integer n,nreal,ltag,itag,lxyz,ixyz,lz,iz,lq,iq,lg,ig,lm,im
      integer i,j,mode,status,version,mlen
      double precision charge,energy,d4param(6)
      character*128 method
      character*256 message
      logical disp,xc_chkdispauto
      logical geom_ncent,geom_cart_get2
      external xc_chkdispauto
      external geom_ncent,geom_cart_get2

      disp=.false.
      if (.not.rtdb_get(rtdb,'dft:disp',mt_log,1,disp)) disp=.false.
      if (disp.or.xc_chkdispauto()) call errquit(
     & 'DFTD4 cannot be combined with DFT DISP/automatic dispersion',
     & 0,INPUT_ERR)
      call xc_dftd4_method(rtdb,method,mlen)

      if (.not.geom_ncent(geom,n))
     & call errquit('xc_dftd4: geom_ncent failed',0,GEOM_ERR)
      if (.not.MA_push_get(MT_Byte,n*16,'d4 tags',ltag,itag))
     & call errquit('xc_dftd4: tags allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,3*n,'d4 xyz',lxyz,ixyz))
     & call errquit('xc_dftd4: xyz allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Int,n,'d4 numbers',lz,iz))
     & call errquit('xc_dftd4: numbers allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,n,'d4 charges',lq,iq))
     & call errquit('xc_dftd4: charges allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Int,n,'d4 atom map',lm,im))
     & call errquit('xc_dftd4: map allocation failed',0,MA_ERR)
      if (.not.geom_cart_get2(geom,n,Byte_MB(itag),Dbl_MB(ixyz),
     & Dbl_MB(iq),Int_MB(iz)))
     & call errquit('xc_dftd4: geom_cart_get2 failed',0,GEOM_ERR)
      nreal=0
      do i=1,n
         if (.not.((Byte_MB(itag+16*(i-1)).eq.'b'.or.
     &               Byte_MB(itag+16*(i-1)).eq.'B').and.
     &              (Byte_MB(itag+16*(i-1)+1).eq.'q'.or.
     &               Byte_MB(itag+16*(i-1)+1).eq.'Q'))) then
            nreal=nreal+1
            Int_MB(im+nreal-1)=i
            Int_MB(iz+nreal-1)=Int_MB(iz+i-1)
            do j=1,3
               Dbl_MB(ixyz+3*(nreal-1)+j-1)=
     &              Dbl_MB(ixyz+3*(i-1)+j-1)
            enddo
         endif
      enddo
      charge=0.0d0
      if (.not.rtdb_get(rtdb,'charge',MT_DBL,1,charge)) charge=0.0d0

      mode=0
      if (what.eq.'forces') mode=1
      if (.not.MA_push_get(MT_Dbl,3*nreal,'d4 derivative',lg,ig))
     & call errquit('xc_dftd4: derivative allocation failed',0,MA_ERR)
      call dfill(3*nreal,0.0d0,Dbl_MB(ig),1)
      energy=0.0d0
      status=0
      version=0
      message=' '
      if (ga_nodeid().eq.0) then
         call nw_dftd4_eval(nreal,Int_MB(iz),Dbl_MB(ixyz),charge,
     &        method,mlen,mode,energy,Dbl_MB(ig),d4param,version,status,
     &        message)
         if (status.ne.0) then
            write(LuOut,'(A,A)') ' DFTD4 error: ',message
            call errquit('DFTD4 evaluation failed',status,
     &                   UNKNOWN_ERR)
         endif
      endif
      call ga_dgop(msg_xc_coef,energy,1,'+')
      if (mode.eq.1)
     & call ga_dgop(msg_xc_coef,Dbl_MB(ig),3*nreal,'+')

      if (what.eq.'energy') then
         value=value+energy
         if (ga_nodeid().eq.0) then
            write(LuOut,100) method(1:mlen),d4param,energy
 100        format(/,
     &       '      DFTD4 Dispersion Correction',/,
     &       '      ---------------------------',/,
     &       '          parameters for functional: ',A,/,
     &       '              s6=',F10.6,'  s8=',F10.6,
     &       '  s9=',F10.6,/,
     &       '              a1=',F10.6,'  a2=',F10.6,
     &       '  alp=',F10.6,/,
     &       '          dispersion energy        : ',F22.12,/)
         endif
      elseif (what.eq.'forces') then
         if (ga_nodeid().eq.0) then
            write(LuOut,105)
            do i=1,nreal
               write(LuOut,110) Int_MB(im+i-1),
     &              (Dbl_MB(ig+3*(i-1)+j-1),j=1,3)
            enddo
 105        format(/,18X,'DFTD4 DISPERSION GRADIENT',/,
     &       5X,'atom',10X,'x',14X,'y',14X,'z')
 110        format(5X,I4,3(1X,F18.12))
            write(LuOut,*)
         endif
         do i=1,nreal
            j=Int_MB(im+i-1)
            call daxpy(3,1.0d0,Dbl_MB(ig+3*(i-1)),1,
     &                 gradient(1+3*(j-1)),1)
         enddo
      else
         call errquit('xc_dftd4: invalid request',0,UNKNOWN_ERR)
      endif

      if (.not.MA_chop_stack(ltag))
     & call errquit('xc_dftd4: MA cleanup failed',0,MA_ERR)
      return
      end

      subroutine xc_dftd4_to_hessian(rtdb)
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
      integer rtdb,geom,n,nreal,ltag,itag,lxyz,ixyz,lz,iz,lq,iq,lm,im
      integer lh,ih,lt,it,lu,iu,n3,n3real,nhess,nhesst,i,j,k,l,ib,jb
      integer status,version
      integer mode,mlen
      double precision charge,energy,d4param(6)
      character*128 method
      character*256 message
      character*(nw_max_path_len) filehess
      logical disp,xc_chkdispauto
      external xc_chkdispauto

      disp=.false.
      if (.not.rtdb_get(rtdb,'dft:disp',mt_log,1,disp)) disp=.false.
      if (disp.or.xc_chkdispauto()) call errquit(
     & 'DFTD4 cannot be combined with DFT DISP/automatic dispersion',
     & 0,INPUT_ERR)
      call xc_dftd4_method(rtdb,method,mlen)
      if (.not.rtdb_cget(rtdb,'geometry',1,message)) message='geometry'
      if (.not.geom_create(geom,message))
     & call errquit('xc_dftd4_hess: geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,message))
     & call errquit('xc_dftd4_hess: geometry load failed',0,GEOM_ERR)
      if (.not.geom_ncent(geom,n))
     & call errquit('xc_dftd4_hess: geom_ncent failed',0,GEOM_ERR)
      n3=3*n
      nhess=n3*n3
      nhesst=n3*(n3+1)/2
      if (.not.MA_push_get(MT_Byte,n*16,'d4h tags',ltag,itag))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,n3,'d4h xyz',lxyz,ixyz))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Int,n,'d4h numbers',lz,iz))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,n,'d4h charges',lq,iq))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Int,n,'d4h atom map',lm,im))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,nhess,'d4 hessian',lh,ih))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,nhess,'old hessian',lt,it))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.MA_push_get(MT_Dbl,nhesst,'tri hessian',lu,iu))
     & call errquit('xc_dftd4_hess: allocation failed',0,MA_ERR)
      if (.not.geom_cart_get2(geom,n,Byte_MB(itag),Dbl_MB(ixyz),
     & Dbl_MB(iq),Int_MB(iz)))
     & call errquit('xc_dftd4_hess: geometry read failed',0,GEOM_ERR)
      nreal=0
      do i=1,n
         if (.not.((Byte_MB(itag+16*(i-1)).eq.'b'.or.
     &               Byte_MB(itag+16*(i-1)).eq.'B').and.
     &              (Byte_MB(itag+16*(i-1)+1).eq.'q'.or.
     &               Byte_MB(itag+16*(i-1)+1).eq.'Q'))) then
            nreal=nreal+1
            Int_MB(im+nreal-1)=i
            Int_MB(iz+nreal-1)=Int_MB(iz+i-1)
            do j=1,3
               Dbl_MB(ixyz+3*(nreal-1)+j-1)=
     &              Dbl_MB(ixyz+3*(i-1)+j-1)
            enddo
         endif
      enddo
      n3real=3*nreal
      charge=0.0d0
      if (.not.rtdb_get(rtdb,'charge',MT_DBL,1,charge)) charge=0.0d0
      call dfill(nhess,0.0d0,Dbl_MB(ih),1)
      mode=2
      energy=0.0d0
      status=0
      message=' '
      if (ga_nodeid().eq.0) then
         call nw_dftd4_eval(nreal,Int_MB(iz),Dbl_MB(ixyz),charge,
     &        method,mlen,
     &        mode,energy,Dbl_MB(ih),d4param,version,status,message)
         if (status.ne.0) then
            write(LuOut,'(A,A)') ' DFTD4 Hessian error: ',message
            call errquit('DFTD4 Hessian failed',status,
     &                   UNKNOWN_ERR)
         endif

         call util_file_name('hess',.false.,.false.,filehess)
         open(unit=69,file=filehess,form='formatted',status='old',
     &        err=900)
         call dfill(nhess,0.0d0,Dbl_MB(it),1)
         do i=1,nhesst
            read(69,*,err=901,end=902) Dbl_MB(iu+i-1)
         enddo
         close(69)
         call trin2squa(Dbl_MB(it),Dbl_MB(iu),n3)
         do i=1,nreal
            do j=1,3
               k=3*(Int_MB(im+i-1)-1)+j
               do l=1,nreal
                  do jb=1,3
                     ib=3*(Int_MB(im+l-1)-1)+jb
                     Dbl_MB(it+k-1+n3*(ib-1))=
     &                Dbl_MB(it+k-1+n3*(ib-1))+
     &                Dbl_MB(ih+3*(i-1)+j-1+
     &                       n3real*(3*(l-1)+jb-1))
                  enddo
               enddo
            enddo
         enddo
         call stpr_wrt_fd_from_sq(Dbl_MB(it),n3,filehess)
         write(LuOut,'(A,A)') ' DFTD4 contribution added to Hessian: ',
     &                        filehess
      endif
      if (.not.MA_chop_stack(ltag))
     & call errquit('xc_dftd4_hess: MA cleanup failed',0,MA_ERR)
      if (.not.geom_destroy(geom))
     & call errquit('xc_dftd4_hess: geom_destroy failed',0,GEOM_ERR)
      call ga_sync
      return
 900  call errquit('xc_dftd4_hess: cannot open Hessian',0,DISK_ERR)
 901  call errquit('xc_dftd4_hess: cannot read Hessian',0,DISK_ERR)
 902  call errquit('xc_dftd4_hess: short Hessian file',0,DISK_ERR)
      end

      subroutine xc_dftd4_method(rtdb,method,mlen)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cdft.fh"
      integer rtdb,mlen,i,nlibxc,libxc_status,native_status
      character*(*) method
      character*128 name1,name2,spec
      method=' '

#ifdef USE_LIBXC
      if (libxcon) then
         call nwchem_libxc_dftd4_names(name1,name2,nlibxc,
     &                                 libxc_status)
         if (libxc_status.gt.0)
     &      call errquit('DFTD4 cannot use scaled LibXC terms',
     &                   libxc_status,INPUT_ERR)
         if (nlibxc.lt.1.or.nlibxc.gt.2)
     &      call errquit('DFTD4 needs one XC functional or an X/C pair',
     &                   nlibxc,INPUT_ERR)
         call xc_dftd4_lower(name1)
         if (nlibxc.eq.1) then
            method=trim(name1)
         else
            call xc_dftd4_lower(name2)
            call xc_dftd4_libxc_pair(name1,name2,method)
         endif
c        DFTD4 4.2.0 lacks these exact LibXC aliases.
         if (trim(method).eq.'hyb_mgga_xc_mn12_sx') method='mn12sx'
         if (trim(method).eq.'hyb_gga_xc_b2plyp') method='b2plyp'
         if (trim(method).eq.'hyb_gga_xc_b2gpplyp')
     &      method='b2gpplyp'
         mlen=len_trim(method)
         return
      endif
#endif

c     Named NWChem combinations whose input spelling differs from the
c     canonical DFTD4 parameter name are mapped explicitly.
      if (xccomb(2)) then
         method='b3lyp'
      elseif (xccomb(4)) then
         method='b97'
      elseif (xccomb(10)) then
         method='pbe0'
      elseif (xccomb(18)) then
         method='tpssh'
      elseif (xccomb(21)) then
         method='mpw1b95'
      elseif (xccomb(22)) then
         method='mpwb1k'
      elseif (xccomb(23)) then
         method='b1b95'
      elseif (xccomb(25)) then
         method='pw6b95'
      elseif (xccomb(33)) then
         method='m06l'
      elseif (xccomb(35)) then
         method='m06'
      elseif (xccomb(37)) then
         method='b97d'
      elseif (xccomb(55)) then
         method='bhlyp'
      elseif (xccomb(56)) then
         method='b3p86'
      elseif (xccomb(57)) then
         method='b3pw91'
      elseif (xccomb(58)) then
         method='pbe'
      elseif (xccomb(59)) then
         method='hse03'
      elseif (xccomb(66)) then
         method='scan'
      elseif (xccomb(71)) then
         method='rscan'
      elseif (xccomb(73)) then
         method='r2scan'
      elseif (xccomb(76)) then
         method='wb97'
      elseif (xccomb(77)) then
         method='wb97x-2008'
      elseif (xccomb(81)) then
         method='mn12sx'
      elseif (xccomb(86)) then
         method='r2scan0'
      else
c        Native NWChem X/C pairs use the original specification so
c        omitted and explicit unit coefficients are equivalent.
         spec=' '
         if (.not.rtdb_cget(rtdb,'dft:xc_spec',1,spec))
     &      call errquit('DFTD4 could not identify the XC functional',
     &                   0,INPUT_ERR)
         call xc_dftd4_lower(spec)
         call xc_dftd4_native_pair(spec,method,native_status)
         if (native_status.eq.2)
     &      call errquit('DFTD4 cannot use scaled native X/C terms',
     &                   0,INPUT_ERR)
         if (native_status.lt.0)
     &      call errquit('DFTD4 cannot map this native X/C pair',
     &                   0,INPUT_ERR)
         if (native_status.eq.1) then
            mlen=len_trim(method)
            return
         endif
c        A single NWChem XC keyword can be passed directly. DFTD4
c        performs
c        the final parameter lookup and reports unsupported functionals.
         method=adjustl(spec)
         do i=1,len_trim(method)
            if (method(i:i).eq.' '.or.method(i:i).eq.char(9))
     &         call errquit(
     &         'DFTD4 needs a named XC functional or LibXC X/C pair',
     &         0,INPUT_ERR)
         enddo
      endif
      mlen=len_trim(method)
      return
      end

      subroutine xc_dftd4_native_pair(spec,method,status)
      implicit none
      character*(*) spec,method
      character*64 word(4),name(2)
      integer status,nword,nname,i,j,ios
      double precision coeff(2),value,eps
      logical xc_dftd4_native_names
      external xc_dftd4_native_names
      parameter (eps=1.0d-12)

      method=' '
      do i=1,4
         word(i)=' '
      enddo
      nword=0
      i=1
 10   if (i.gt.len_trim(spec)) goto 30
      if (spec(i:i).eq.' '.or.spec(i:i).eq.char(9)) then
         i=i+1
         goto 10
      endif
      if (nword.eq.4) then
         status=0
         return
      endif
      nword=nword+1
      j=i
 20   if (j.le.len_trim(spec)) then
         if (spec(j:j).ne.' '.and.spec(j:j).ne.char(9)) then
            j=j+1
            goto 20
         endif
      endif
      word(nword)=spec(i:j-1)
      i=j+1
      goto 10

 30   nname=0
      coeff(1)=1.0d0
      coeff(2)=1.0d0
      do i=1,nword
         read(word(i),*,iostat=ios) value
         if (ios.eq.0) then
            if (nname.lt.1.or.nname.gt.2) then
               status=0
               return
            endif
            coeff(nname)=value
         else
            nname=nname+1
            if (nname.gt.2) then
               status=0
               return
            endif
            name(nname)=word(i)
         endif
      enddo
      status=0
      if (nname.ne.2) return
      if (abs(coeff(1)-1.0d0).gt.eps.or.
     &    abs(coeff(2)-1.0d0).gt.eps) then
         status=2
         return
      endif

      status=1
      if (xc_dftd4_native_names(name(1),name(2),'becke88','lyp')) then
         method='blyp'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'becke88','perdew86')) then
         method='bp86'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'becke88','perdew91')) then
         method='bpw'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'becke88','cpbe96')) then
         method='bpbe'
      elseif (xc_dftd4_native_names(name(1),name(2),'mpw91','lyp')) then
         method='mpwlyp'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'mpw91','perdew91')) then
         method='mpwpw91'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xpbe96','cpbe96')) then
         method='pbe'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'rpbe','cpbe96')) then
         method='rpbe'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'revpbe','cpbe96')) then
         method='revpbe'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xperdew91','perdew91')) then
         method='pw91'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xperdew91','perdew86')) then
         method='pwp'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xperdew86','cpbe96')) then
         method='pw86pbe'
      elseif (xc_dftd4_native_names(name(1),name(2),'optx','lyp')) then
         method='olyp'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'optx','cpbe96')) then
         method='opbe'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xtpss03','ctpss03')) then
         method='tpss'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xscan','cscan')) then
         method='scan'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xr2scan','cr2scan')) then
         method='r2scan'
      elseif (xc_dftd4_native_names(name(1),name(2),'xm06','cm06')) then
         method='m06'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'xm06-l','cm06-l')) then
         method='m06l'
      elseif (xc_dftd4_native_names(name(1),name(2),
     &                              'gill96','lyp')) then
         method='glyp'
      else
         status=-1
      endif
      return
      end

      logical function xc_dftd4_native_names(name1,name2,xname,cname)
      implicit none
      character*(*) name1,name2,xname,cname
      xc_dftd4_native_names=(
     & (trim(name1).eq.trim(xname).and.trim(name2).eq.trim(cname)).or.
     & (trim(name2).eq.trim(xname).and.trim(name1).eq.trim(cname)))
      return
      end

      subroutine xc_dftd4_libxc_pair(name1,name2,method)
      implicit none
      character*(*) name1,name2,method
      logical xc_dftd4_pair
      external xc_dftd4_pair
      method=trim(name1)//':'//trim(name2)
      if (xc_dftd4_pair(name1,name2,'gga_x_am05','gga_c_am05'))
     &   method='am05'
      if (xc_dftd4_pair(name1,name2,'gga_x_b88','gga_c_lyp'))
     &   method='blyp'
      if (xc_dftd4_pair(name1,name2,'gga_x_b88','gga_c_pbe'))
     &   method='bpbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_b88','gga_c_p86'))
     &   method='bp86'
      if (xc_dftd4_pair(name1,name2,'gga_x_b88','gga_c_pw91'))
     &   method='bpw'
      if (xc_dftd4_pair(name1,name2,'gga_x_mpw91','gga_c_lyp'))
     &   method='mpwlyp'
      if (xc_dftd4_pair(name1,name2,'gga_x_mpw91','gga_c_pw91'))
     &   method='mpwpw91'
      if (xc_dftd4_pair(name1,name2,'gga_x_optx','gga_c_lyp'))
     &   method='olyp'
      if (xc_dftd4_pair(name1,name2,'gga_x_optx','gga_c_pbe'))
     &   method='opbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_pbe','gga_c_pbe'))
     &   method='pbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_rpbe','gga_c_pbe'))
     &   method='rpbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_pbe_r','gga_c_pbe'))
     &   method='revpbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_pbe_sol',
     &                               'gga_c_pbe_sol')) method='pbesol'
      if (xc_dftd4_pair(name1,name2,'gga_x_pw86','gga_c_pbe'))
     &   method='pw86pbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_rpw86','gga_c_pbe'))
     &   method='rpw86pbe'
      if (xc_dftd4_pair(name1,name2,'gga_x_pw91','gga_c_pw91'))
     &   method='pw91'
      if (xc_dftd4_pair(name1,name2,'gga_x_pw91','gga_c_p86'))
     &   method='pwp'
      if (xc_dftd4_pair(name1,name2,'mgga_x_tpss','mgga_c_tpss'))
     &   method='tpss'
      if (xc_dftd4_pair(name1,name2,'mgga_x_revtpss',
     &                               'mgga_c_revtpss')) method='revtpss'
      if (xc_dftd4_pair(name1,name2,'mgga_x_scan','mgga_c_scan'))
     &   method='scan'
      if (xc_dftd4_pair(name1,name2,'mgga_x_rscan','mgga_c_rscan'))
     &   method='rscan'
      if (xc_dftd4_pair(name1,name2,'mgga_x_r2scan',
     &                               'mgga_c_r2scan')) method='r2scan'
      if (xc_dftd4_pair(name1,name2,'mgga_x_m06','mgga_c_m06'))
     &   method='m06'
      if (xc_dftd4_pair(name1,name2,'mgga_x_m06_l','mgga_c_m06_l'))
     &   method='m06l'
      if (xc_dftd4_pair(name1,name2,'gga_x_g96','gga_c_lyp'))
     &   method='glyp'
      return
      end

      logical function xc_dftd4_pair(name1,name2,want1,want2)
      implicit none
      character*(*) name1,name2,want1,want2
      xc_dftd4_pair=(trim(name1).eq.trim(want1).and.
     &                trim(name2).eq.trim(want2)).or.
     &               (trim(name1).eq.trim(want2).and.
     &                trim(name2).eq.trim(want1))
      return
      end

      subroutine xc_dftd4_lower(text)
      implicit none
      character*(*) text
      integer i,k
      do i=1,len_trim(text)
         k=ichar(text(i:i))
         if (k.ge.ichar('A').and.k.le.ichar('Z'))
     &      text(i:i)=char(k+ichar('a')-ichar('A'))
      enddo
      return
      end
