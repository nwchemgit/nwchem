      subroutine dec_fpe
*
* $Id: dec_fpe.f,v 1.2 1997-10-31 20:45:31 d3e129 Exp $
*
      implicit integer*4 (a-z)
      include '/usr/include/for_fpe_flags.f'
      data first /1/
c
c     Disable underflow traps on the Alpha.  C6H6 in 6-31g 
c     fails inside texas if this is not done.  Don't know why.
c
      if (first.gt.1) return
      first = first + 1
c
      mask = for_get_fpe()
c$$$      write(6,*) ' FPE_M_TRAP_UND ' , iand(FPE_M_TRAP_UND,mask)
c$$$      write(6,*) ' FPE_M_TRAP_OVF ' , iand(FPE_M_TRAP_OVF,mask)
c$$$      write(6,*) ' FPE_M_TRAP_DIV0 ' , iand(FPE_M_TRAP_DIV0,mask)
c$$$      write(6,*) ' FPE_M_TRAP_INV ' , iand(FPE_M_TRAP_INV,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_00 ' , iand(FPE_M_RESERVED_00,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_01 ' , iand(FPE_M_RESERVED_01,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_02 ' , iand(FPE_M_RESERVED_02,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_03 ' , iand(FPE_M_RESERVED_03,mask)
c$$$      write(6,*) ' FPE_M_MSG_OVF ' , iand(FPE_M_MSG_OVF,mask)
c$$$      write(6,*) ' FPE_M_MSG_UND ' , iand(FPE_M_MSG_UND,mask)
c$$$      write(6,*) ' FPE_M_MSG_DIV0 ' , iand(FPE_M_MSG_DIV0,mask)
c$$$      write(6,*) ' FPE_M_MSG_INV ' , iand(FPE_M_MSG_INV,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_04 ' , iand(FPE_M_RESERVED_04,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_05 ' , iand(FPE_M_RESERVED_05,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_06 ' , iand(FPE_M_RESERVED_06,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_07 ' , iand(FPE_M_RESERVED_07,mask)
c$$$      write(6,*) ' FPE_M_ABRUPT_UND ' , iand(FPE_M_ABRUPT_UND,mask)
c$$$      write(6,*) ' FPE_M_ABRUPT_OVF ' , iand(FPE_M_ABRUPT_OVF,mask)
c$$$      write(6,*) ' FPE_M_ABRUPT_DIV0 ' , iand(FPE_M_ABRUPT_DIV0,mask)
c$$$      write(6,*) ' FPE_M_ABRUPT_INV ' , iand(FPE_M_ABRUPT_INV,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_08 ' , iand(FPE_M_RESERVED_08,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_09 ' , iand(FPE_M_RESERVED_09,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_10 ' , iand(FPE_M_RESERVED_10,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_11 ' , iand(FPE_M_RESERVED_11,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_12 ' , iand(FPE_M_RESERVED_12,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_13 ' , iand(FPE_M_RESERVED_13,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_14 ' , iand(FPE_M_RESERVED_14,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_15 ' , iand(FPE_M_RESERVED_15,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_16 ' , iand(FPE_M_RESERVED_16,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_17 ' , iand(FPE_M_RESERVED_17,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_18 ' , iand(FPE_M_RESERVED_18,mask)
c$$$      write(6,*) ' FPE_M_RESERVED_19 ' , iand(FPE_M_RESERVED_19,mask)
c
      mask = 0
      mask = ior(mask, FPE_M_TRAP_OVF)
      mask = ior(mask, FPE_M_TRAP_DIV0)
      mask = ior(mask, FPE_M_TRAP_INV)
      mask = ior(mask, FPE_M_ABRUPT_UND)
c
      mask = for_set_fpe(mask)
c
      end
