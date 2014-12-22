      subroutine ccsd_tengy(f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
     &                      dintc,dintx,t1v,eorb,eaijk,
     $                      emp4,emp5,ncor,nocc,nvir)
      implicit none
#include "ccsdps.fh"
      integer nvir,ncor,nocc
      integer b,c
      double precision emp4,emp5,denom,eaijk
      double precision one,two,three,four
      double precision f1n(nvir,nvir),f1t(nvir,nvir)
      double precision f2n(nvir,nvir),f2t(nvir,nvir)
      double precision f3n(nvir,nvir),f3t(nvir,nvir)
      double precision f4n(nvir,nvir),f4t(nvir,nvir)
      double precision dintc(nvir),dintx(nvir),t1v(nvir)
      double precision eorb(*)
      data one/1.0d00/,two/2.0d00/,three/3.0d00/,four/4.0d00/
c
      if (occsdps) then
         call pstat_on(ps_tengy)
      else
         call qenter('tengy',0)
      endif
c
      do b=1,nvir
        do c=1,nvir
          denom=one/(eaijk-eorb(ncor+nocc+b)-eorb(ncor+nocc+c))
          emp4=emp4+denom*
     &         (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &         (f1t(b,c)-two*f2t(b,c)-two*f3t(b,c)+f4t(b,c))
c
          emp4=emp4-denom*
     &         (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &         (two*f1t(b,c)-f2t(b,c)-f3t(b,c)+two*f4t(b,c))
c
          emp4=emp4+three*denom*(
     &         f1n(b,c)*(f1n(b,c)+f3n(c,b)+two*f4t(c,b))+
     &         f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
c
          emp5=emp5+denom*t1v(b)*dintx(c)*
     &        (     f1t(b,c)+f2n(b,c)+f4n(c,b)
     &         -two*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &              f1n(b,c)+f2t(b,c)+f3n(c,b))
     &         +four*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
          emp5=emp5+denom*t1v(b)*dintc(c)*
     &        (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &          -two*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
c
        enddo
      enddo
c
      if (occsdps) then
         call pstat_off(ps_tengy)
      else
         call qexit('tengy',0)
      endif
c
      return
      end
