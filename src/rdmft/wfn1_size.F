c-----------------------------------------------------------------------
c
      integer function wfn1_ibd_proc(ndouble,nproc)
      implicit none
C>
C> \brief Return the number of bytes per processor needed to store the
C> given number of doubles
C>
C> This function assumes that the memory requirements are equally 
C> distributed. Hence the number may be inaccurate for highly 
C> inhomogenous distributions.
C>
#include "mafdecls.fh"
c
      integer, intent(in) :: ndouble ! the number of doubles
      integer, intent(in) :: nproc   ! the number of processors
c
      wfn1_ibd_proc = ma_sizeof(MT_BYTE,1,MT_DBL)*
     &                (ndouble+nproc-1)/nproc
c
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function wfn1_ibi_proc(ninteger,nproc)
      implicit none
C>
C> \brief Return the number of bytes per processor needed to store the
C> given number of integers
C>
C> This function assumes that the memory requirements are equally 
C> distributed. Hence the number may be inaccurate for highly 
C> inhomogenous distributions.
C>
#include "mafdecls.fh"
c
      integer, intent(in) :: ninteger ! the number of integers
      integer, intent(in) :: nproc    ! the number of processors
c
      wfn1_ibi_proc = ma_sizeof(MT_BYTE,1,MT_INT)*
     &                (ninteger+nproc-1)/nproc
c
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function wfn1_ibd(ndouble)
      implicit none
C>
C> \brief Return the number of bytes needed to store the given number of
C> doubles
C>
#include "mafdecls.fh"
c
      integer, intent(in) :: ndouble ! the number of doubles
c
      wfn1_ibd = ma_sizeof(MT_BYTE,1,MT_DBL)*ndouble
c
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function wfn1_ibi(ninteger)
      implicit none
C>
C> \brief Return the number of bytes needed to store the given number of
C> integers
C>
#include "mafdecls.fh"
c
      integer, intent(in) :: ninteger
c
      wfn1_ibi = ma_sizeof(MT_BYTE,1,MT_DBL)*ninteger
c
      return
      end
c
c-----------------------------------------------------------------------
