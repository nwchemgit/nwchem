      subroutine argos_wtrest(lfn)
c
      implicit none
c
#include "argos.fh"
c
      integer lfn
c
      if(me.ne.0) return
c
      write(lfn,1000)
 1000 format('restart input')
      write(lfn,1001) ntype,mdtype
 1001 format(11i7)
      write(lfn,1001) nfpair,nflong
      write(lfn,1001) lwtype,lstype,nfrest,keepr
      write(lfn,1002) irun,iequi,idacq,mrun,mequi,mdacq,ldacq
 1002 format(11i7)
      write(lfn,1003) stime,tstep
 1003 format(2f12.6)
      write(lfn,1004) rshort,rlong
 1004 format(2f12.6)
      write(lfn,1005) mshitw,tlwsha
      write(lfn,1005) mshits,tlssha
 1005 format(i7,f12.6)
      write(lfn,1006) ipscal,prsext,prsrlx,compr,ipopt
 1006 format(i5,e12.5,f12.6,e12.5,i5)
      write(lfn,1007) itscal,tmpext1,tmprlx,tmsrlx,tmpext2,
     + tann1,tann2
 1007 format(i5,6f12.6)
      tmpext=tmpext1
      write(lfn,1008) nfgaus,ivopt,tgauss,iseed
 1008 format(2i7,f12.6,i12)
      write(lfn,1009) nfoutp,nfstat,nfprop,nfnewf,ibatch
 1009 format(11i7)
      write(lfn,1009) ibinar,iformt
      write(lfn,1009) nfcoor,nfscoo,nfvelo,nfsvel,nfforc,nfsfor
      write(lfn,1009) nffree
      write(lfn,1009) nfcntr,nfslow
      write(lfn,1010) nfrdf,numrdf,ngc,ngr,ngl,ngrww,ngrsw,ngrss
 1010 format(8i7)
      write(lfn,1011) rrdf,drdf
 1011 format(2f12.6)
      write(lfn,1012) numdis,lendis
 1012 format(2i7)
      write(lfn,1013) numhis,lenhis
 1013 format(11i7)
      write(lfn,1014) nfdip,ndip,rdip
 1014 format(2i7,f12.6)
      write(lfn,1015) nfkirk,nkirk,rkirk
 1015 format(2i7,f12.6)
c
      return
      end
c $Id$
