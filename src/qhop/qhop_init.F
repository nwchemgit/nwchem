      subroutine qhop_init(lfntop,filtop,lfnout,me)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer lfntop,lfnout,me
      integer i,numbi,numbd,idata(9)
      character*255 filtop
      character*4 string
c
      if(me.eq.0) then
c
c     reading in the qhop par files
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(lfntop)
   99 continue
      read(lfntop,1001,end=9999,err=9999) string
 1001 format(a4)
      if(string(1:4).ne.'qhop') goto 99
c
      read(lfntop,1011) mxseq,mxar,maxpar,maxseq
 1011 format(2i5,i10,i5)
c
      read(lfntop,1025) (shbit(i),i=1,5)
 1025 format(5i5)
c
      endif
c
      idata(1)=mxseq
      idata(2)=mxar
      idata(3)=maxpar
      idata(4)=maxseq
      do 1 i=1,5
      idata(4+i)=shbit(i)
    1 continue
      numbi=ma_sizeof(mt_int,1,mt_byte)
      call ga_brdcst(hop_01,idata,9*numbi,0)
      mxseq=idata(1)
      mxar=idata(2)
      maxpar=idata(3)
      maxseq=idata(4)
      do 2 i=1,5
      shbit(i)=idata(4+i)
    2 continue
c
c     allocate memory for parameter arrays
c
      if(.not.ma_push_get(mt_int,2*maxseq,'iptseq',l_ptseq,i_ptseq))
     + call md_abort('Failed to allocate iptseq',0)
      if(.not.ma_push_get(mt_int,mxseq,'iptarat',l_ptarat,i_ptarat))
     + call md_abort('Failed to allocate iptarat',0)
      if(.not.ma_push_get(mt_int,4*mxar,'iarat',l_iarat,i_iarat))
     + call md_abort('Failed to allocate iarat',0)
      if(.not.ma_push_get(mt_dbl,2*mxar,'racs',l_racs,i_racs))
     + call md_abort('Failed to allocate racs',0)
      if(.not.ma_push_get(mt_int,maxpar,'iptpar',l_ptpar,i_ptpar))
     + call md_abort('Failed to allocate iptpar',0)
      if(.not.ma_push_get(mt_dbl,31*maxpar,'par',l_par,i_par))
     + call md_abort('Failed to allocate par',0)
c
      call qhop_init2(lfntop,lfnout,me,
     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
     + dbl_mb(i_par+28*maxpar))
c
      close(unit=lfntop)
c
c     broadcast parameters to all processors
c     --------------------------------------
c
      numbi=ma_sizeof(mt_int,1,mt_byte)
      numbd=ma_sizeof(mt_dbl,1,mt_byte)
c
      call ga_brdcst(hop_02,int_mb(i_ptseq),2*maxseq*numbi,0)
      call ga_brdcst(hop_03,int_mb(i_ptarat),mxseq*numbi,0)
      call ga_brdcst(hop_04,int_mb(i_iarat),4*mxar*numbi,0)
      call ga_brdcst(hop_05,int_mb(i_ptpar),maxpar*numbi,0)
      call ga_brdcst(hop_06,dbl_mb(i_racs),2*mxar*numbd,0)
      call ga_brdcst(hop_07,dbl_mb(i_par),31*maxpar*numbd,0)
c
c      call qhop_test(lfnout,me,
c     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
c     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
c     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
c     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
c     + dbl_mb(i_par+28*maxpar))
c
      return
 9999 call md_abort('Error reading qhop data from topology file',me)
      return
      end
      subroutine qhop_init2(lfntop,lfnout,me,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     input variables
c     ---------------
c
c     matm            : dimension of atom arrays
c     natm            : number of atoms
c     mseq            : dimension of residue arrays
c     nseq            : number of residues
c     nparms          : number of parameter sets
c     mparms          : number of parameter sets plus 1
c     cseq            : residue names
c     catm(1,1:natm)  : atom names
c     catm(k,1:natm)  : atom types in parameter set k-1
c     latm(5,1:natm)  : residue number
c     latm(10,1:natm) : + heavy atom attached to protonatable hydrogen
c                       - protonatable hydrogen, relative index to heavy atom
c     lseq(5,1:nseq)  : number of protonation states of residue
c     lseq(6,1:nseq)  : current protonation state of residue
c
      integer i,j,ipar,numpar,larat,numseq,nseq,lfntop,lfnout,me
c
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
c
      integer seqin(2),hvin(2),prtin(2)	 
      character*10 dummy
c
      integer numbi,numbd,numb
      character*4 string
c
      if(me.eq.0) then
c
      rewind(lfntop)
   99 continue
      read(lfntop,1000,end=9999,err=9999) string
 1000 format(a4)
      if(string(1:4).ne.'qhop') goto 99
c
c     reading in the qhop par files
c
      read(lfntop,1011) numseq,larat,numpar,nseq
 1011 format(2i5,i10,i5)
c
      read(lfntop,1025) (shbit(i),i=1,5)
 1025 format(5i5)
c
      do 1 i=1,nseq
      read(lfntop,1010) j,dummy,(ptseq(i,j),j=1,2)
 1010 format(i5,1x,a10,2i5)
    1 continue
c
      do 2 i=1,numseq
      read(lfntop,1020) j,dummy,ptarat(i)
 1020 format(i5,1x,a10,1x,i5)
    2 continue
c
      do 3 i=1,larat
      read(lfntop,1021) j,(arat(i,j),j=1,4), deq(i),aneq(i)  
 1021 format(i5,4i4,2f8.4)
    3 continue
c
      do 4 ipar=1,numpar 
      read(lfntop,1001)ptpar(ipar)
 1001 format(i10)
c       
      read(lfntop,1003) (e12fxy(ipar,i),i=1,3)
 1003 format(3f12.6)
      read(lfntop,1004) (defe12(ipar,i),i=1,7)
 1004 format(3f12.6,f8.4,2f12.6,f12.8)
      read(lfntop,1005) (zpef(ipar,i),i=1,6)
 1005 format(6f12.6)
      read(lfntop,1006) (tdsgl(ipar,i),i=1,6)
 1006 format(6f12.6)
      read(lfntop,1007) (tunnel(ipar,i),i=1,7)
 1007 format(f8.2,f8.3,f10.5,e12.3,f10.5,2e12.3)
      read(lfntop,1008) (tstval(ipar,i),i=1,5)
 1008 format(f10.2,2f10.4,f10.6,f10.4) 
      write(lfnout,1009)ptpar(ipar)
 1009 format(i10)
      write(lfnout,1003)(e12fxy(ipar,i),i=1,3)     	
      write(lfnout,1004)(defe12(ipar,i),i=1,7)
      write(lfnout,1005)(zpef(ipar,i),i=1,6)
      write(lfnout,1006)(tdsgl(ipar,i),i=1,6)
      write(lfnout,1007)(tunnel(ipar,i),i=1,7)
      write(lfnout,1008)(tstval(ipar,i),i=1,5)
    4 continue
c
      endif
c
      return
  112 call md_abort('err opening input file',0)
      return
  113 call md_abort( 'err opening output file',0)
      return
 9999 call md_abort('Error reading qhop data from topology file',me)
      return
      end
      subroutine qhop_test(lfnout,me,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy)
c
      implicit none
c
#include "qhop_common.fh"
c
      integer lfnout,me
c
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
c
      integer i
c
c      write(*,'(i5,5f12.6)') me,(tunnel(i,1),i=1,min(5,maxpar))
c
      return
      end
