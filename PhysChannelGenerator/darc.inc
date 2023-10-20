C-----------------------------------------------------------------------
C   darc.inc
C   This file is INCLUDEd in the DARC modules.
C   It is used to set dimensions. 
C   The following table indicates where the dimensions are used.
C   This version is for use with dstg0,1,2 only, then use dto3 etc.
C
C           dstg0  dstg1  dstg2
C
C   MXCH                    *  
C   MXI1             *      *
C   MXI2             *      *
C   MXNB             *      *  
C   MXNC                    *
C   MXNK             *      *  
C   MXNL                    *  
C   MXNM                    *
C   MXNP      *      *     
C   MXNW      *             *
C   MXP1      *      *
C
C-----------------------------------------------------------------------
C   MXCH   channels
C   MXI1   radial integrals
C   MXI2   marks position of radial integrals
C   MXNB   continuum orbitals per KAPPA value
C   MXNC   relativistic configurations
C   MXNK   KAPPA values
C   MXNL   target levels or correlation functions
C   MXNM   MCP angular coefficients
C   MXNP   points in radial mesh
C   MXNW   relativistic orbitals
C   MXP1   mesh points used to store orbitals
C-----------------------------------------------------------------------
      INTEGER MXCH,MXI1,MXI2,MXNB,MXNC
      INTEGER MXNK,MXNL,MXNM,MXNP,MXNW,MXP1
C     jinrui 20141224
      INTEGER NIRX
      PARAMETER (NIRX=20)
      PARAMETER (MXCH = 1000)
      PARAMETER (MXI1 = 40000000)
      PARAMETER (MXI2 = 2900)
      PARAMETER (MXNB = 80)
      PARAMETER (MXNC = 80000)
      PARAMETER (MXNK = 100)
      PARAMETER (MXNL =21000)
      PARAMETER (MXNM = 80000000)
      PARAMETER (MXNP = 5000)
      PARAMETER (MXNW = 500)
      PARAMETER (MXP1 = 900000)
C-----------------------------------------------------------------------
