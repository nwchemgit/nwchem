      subroutine qhop_init(lfnhop,filhop,lfnout,me)
c
c $Id: qhop_init.F,v 1.1 2003-07-17 19:36:51 d3j191 Exp $
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
      integer i,j,ipar,numpar,larat,numseq,nseq,lfnhop,lfnout,me
      integer seqin(2),hvin(2),prtin(2)	 
      character*10 dummy
      character*255 filhop
c
      integer numbi,numbd,numb
c
      if(me.eq.0) then
c
      open(unit=lfnhop,file=filhop(1:index(filhop,' ')-1),err=112)
c
c     reading in the qhop par files
c
      read(lfnhop,1011) numseq,larat,numpar,nseq
 1011 format(2i5,i10,i5)
c
      read(lfnhop,1025) (shbit(i),i=1,5)
 1025 format(5i5)
c
      do 1 i=1,nseq
      read(lfnhop,1010) j,dummy,(ptseq(i,j),j=1,2)
 1010 format(i5,1x,a10,2i5)
    1 continue
c
      do 2 i=1,numseq
      read(lfnhop,1020) j,dummy,ptarat(i)
 1020 format(i5,1x,a10,1x,i5)
    2 continue
c
      do 3 i=1,larat
      read(lfnhop,1021) j,(arat(i,j),j=1,4), deq(i),aneq(i)  
 1021 format(i5,4i4,2f8.4)
    3 continue
c
      do 4 ipar=1,numpar 
      read(lfnhop,1001)ptpar(ipar)
 1001 format(i10)
c       
      read(lfnhop,1003) (e12fxy(ipar,i),i=1,3)
 1003 format(3f12.6)
      read(lfnhop,1004) (defe12(ipar,i),i=1,7)
 1004 format(3f12.6,f8.4,2f12.6,f12.8)
      read(lfnhop,1005) (zpef(ipar,i),i=1,6)
 1005 format(6f12.6)
      read(lfnhop,1006) (tdsgl(ipar,i),i=1,6)
 1006 format(6f12.6)
      read(lfnhop,1007) (tunnel(ipar,i),i=1,7)
 1007 format(f8.2,f8.3,f10.5,e12.3,f10.5,2e12.3)
      read(lfnhop,1008) (tstval(ipar,i),i=1,5)
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
      close(lfnhop)
c
      endif
c
c     broadcast parameters to all processors
c     --------------------------------------
c
      numbi=ma_sizeof(mt_int,1,mt_byte)
      numbd=ma_sizeof(mt_dbl,1,mt_byte)
c
      numb=numbd*(31*maxpar+2*mxar)+numbi*(3*maxseq+maxpar+4*mxar+5)
c
      call ga_brdcst(hop_01,tunnel,numb,0)
c
      return
  112 call errquit('err opening input file',0)
      return
  113 call errquit( 'err opening output file',0)
      return
      end
