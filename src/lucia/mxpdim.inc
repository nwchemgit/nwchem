* contains all PARAMETERS defining LUCIA
      PARAMETER(MXPLNC = 102)
      PARAMETER (MXPIRR = 8)
      PARAMETER (MXPOBS = 8 )
      PARAMETER (MXPR4T = 10 )
      PARAMETER (MXPR0T = 1 )
      PARAMETER (MXPORB = 500)
      PARAMETER (MXPXOT = 9)
      PARAMETER (MXPXST = 100 )
      PARAMETER (MXPSHL = 100)
      PARAMETER (MXPL = 20 )
      PARAMETER (MXPXT = 25)
      PARAMETER (MXPICI = 30)
      PARAMETER (MXPSTT = 30000)
      PARAMETER (MXPCSM = 8)
      PARAMETER (MXPCTP = 30)
      PARAMETER (MXCNSM = 8)
      PARAMETER (MTYP = 30)
      PARAMETER (MXP_NINTP_FOR_EX=50)
      PARAMETER (MXP2EIARR = 10)
      PARAMETER (MXPORBSETS = 10)
      PARAMETER (MXPOPORB = 50)
      PARAMETER (MXPNSCRFIL = 10)
*
      PARAMETER (MXPNEL = 1000)
*
*. Note : MXPNGAS = MXPR4T+6 !!
*. Required in order to handle GAS and RAS within /LUCINP/   
      PARAMETER (MXPNGAS = 16)
      PARAMETER(MXPNSMST = 8)
*. Largest allowed division of space for perturbation operator
      PARAMETER (MXPPTSPC=20)
*. Largest cycle of spaces in multispace method
      PARAMETER(MXPCYC = 100)
*
      PARAMETER(MXPTSOB = 60) 
*. Largest number of spin creation or annihilations in CC operators
      PARAMETER(MXPLCCOP = 20)
*. For the product wave-function business 
*. Max number of orbitals per orbitalsubspace
      PARAMETER(MXP_NOBSBSPC_PRDWF = 20)
*. Max number of wavefunction per orbital subspace
      PARAMETER(MXP_NWF_SUBSPC_PRDWF=20)
*. Allowed number of fragments of molecule - for VB project
      PARAMETER (MXPFRAG = 10)
*. Allowed number of supersymmetry-irreps- derived from 
*. MXPL assuming Dinfh
      PARAMETER(MXP_NIRREP_SUPSYM = 2*(2*MXPL+1))
*. Allowed number of supersymmetry-irreps- derived from MXPL assuming Dinfh
      PARAMETER(MXP_NSUPSYM = 2*(MXPL + 1 + MXPL*(MXPL+1)))
*. Allowed number of VB general spaces
      PARAMETER(MXPNVBGNPSC = 10)
*. Dimension of work 
      PARAMETER (MXPWRD = 1)
*
      
