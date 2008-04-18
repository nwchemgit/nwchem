      integer natms,atmtype,typmol

      real*8  ccc,vvv,fff
      real*8  typmass,typchge

      character typname*8

      dimension atmtype(mxatms)
      dimension ccc(mxatms,3),vvv(mxatms,3),fff(mxatms,3)
      dimension typmass(mxtype),typchge(mxtype)
      dimension typname(mxtype),typmol(mxtype)

      common/cm_atom/atmtype,ccc,vvv,fff,natms
      common/cm_type/typmass,typchge,typname,typmol
