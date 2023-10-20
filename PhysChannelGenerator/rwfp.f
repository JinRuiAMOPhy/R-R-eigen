      IMPLICIT REAL*8(A-H,O-Z)
	character fname*20
      PARAMETER(NMX=1000,NMWF=12000,jm=35,MAXWF=350)
      dimension WMU(MAXWF),IPHZ(MAXWF)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      COMMON /NP/ RNK,ANK,VNK,NK
c gaoxiang 2012.10.6 to read in the structure factor to form boundary condition
      common /gaobnd/ skli
c gaoxiang 2012.10.6
      DATA ALF/137.03602D0/
      open(5,file='rwf.in',status='old')
      open(1,file='radata.dat',status='unknown')
      open(2,file='rwf.out',form='formatted')
      open(22,file='WF.DAT',form='unformatted')
c gaoxiang 2012.10.6
	open(13,file='eband.txt',form='formatted',status='unknown')
C jinrui  2014.04.04
      open(10,file='grid.out',status='unknown')
	open(11,file='grid.in',status='unknown')
      READ (5,5000) ZN,DELL,RMAX
      READ (5,5001) NWF,NWFO,NDEBUG,NMIN,NFACT
 5001 FORMAT(8I10)
      IF (DELL.EQ.0.D0) DELL=1.D-14
      NSTEP=50
      Z=ZN
      XION=Z-1.D0
      RNUC=0.D0
      NB=10
      D=2.
      R=10.D0/D**NB
      DO 1 I=1,NB
      R=R*D
      RH(I)=R
    1 VRX(I)=0.D0
      V0=0.D0
      RFAC=NFACT
      IF (NFACT.EQ.0) RFAC=1.D0
      IF (ZN.NE.0.D0) GO TO 6992
       LOG=1
      CALL RDEN(LOG,NB)
      LEIZN=ZN
      HP=1.D0/H
      Z=ZN
      IF(RNUC.GT.0.D0) Z=0.
      DO 2 I=1,NB
    2 VRX(I)=Z-ZN*VRX(I)
 6992 ZATOM=ZN
      XKOOP=0.D0
      IF(NDEBUG.LT.3) GOTO 4
      WRITE (6,6000) NB
 6000 FORMAT(10I10)
      WRITE (6,6001) (VRX(I),I=1,NB)
 6001 FORMAT(1P5E15.7)
      WRITE (6,6001) (RH(I),I=1,NB)
      WRITE (6,6001) V0
    4 Z=ZN
CC------P.H.Z.-----------------
CC  FOR CORRECTING THE QUANTUM DEFECT OF CONTINUUM
      DO 111,I=1,MAXWF
      IPHZ(I)=0
      WMU(I)=0.D0
 111  CONTINUE
      JPHZ=0
      FPHZL=0.D0
      FPHZJ=0.D0
CC---------------------
c**************HAN.050610****************
	NNNHAN=0
	RRRHAN=0.0
	RDELN=0.0
	RDEL1=0.0
      RNTHANNNNN=0.0
	HANNNNN=0.0 
c**************HAN.050610****************

c gaoxiang 2012.10.3 reset the initia value of orbital energy
      esave=-(ZN)**2
c gaoxiang 2012.10.3
      DO 100 IW=1,NWF
c gaoxiang 2012.10.6 to read in the structure factor to form boundary condition
c      READ (5,5000) FN,FL,FJ,E,RFA,SCALE
      READ (5,5000) FN,FL,FJ,E,RFA,SCALE,skli
5002  format(6F10.4,E20.10)      
c gaoxiang 2012.10.6
CC------P.H.Z--------------------------
      IF(FPHZL.NE.FL.OR.FPHZJ.NE.FJ) JPHZ=JPHZ+1
      IPHZ(JPHZ)=IPHZ(JPHZ)+1
      FPHZJ=FJ
      FPHZL=FL
CC---------------------------
      IF(SCALE.EQ.0.D0) SCALE=50.D0
 5000 FORMAT(8F10.4)
      HMIN=96.+(1.-ZATOM)*0.79
      MMIN=HMIN
      MMIN=(MMIN/8+1)*8
      MMIN=MAX0(MMIN,32)
      NMIN=MAX0(MMIN,NMIN)
      HMIN=DFLOAT(NMIN)
      ZN=ZATOM
      RUNIT=1.D0/(ZATOM*SCALE)
      Z=ZN
      IF(RNUC.GT.0.D0)  Z=0.D0
c      IF(E.LT.0.D0) GOTO 3700
c      FACT=RFA
C      IF(RFA.EQ.0.D0) FACT=10.
c
c---- Xiao-Min changes FACT  from 10 to 10, Sept. 23, 1994
c
c      IF(RFA.EQ.0.D0) FACT=10.
c      FACL=dfloat(NFACT)
c      IF(NFACT.EQ.0) FACL=16.D0
c      P1=DSQRT(E*(2.D0+E/ALF**2))
c      RENG=1.D0/P1
c      DP=RENG/FACL
c      RENG=RENG*FACT
c      RN=RENG
c      ZI=XION+1.D0
c      IF(ZI.EQ.ZATOM) GOTO 3605
c      RN=RH(NB)
c      ZSC=VRX(NB)
c      DO 3600 I=1,NB
c      II=NB-I+1
c      IF(VRX(II).LT.ZSC) GOTO 3605
c 3600 RN=RH(II)
c 3605 RN=DMAX1(RN,RENG)
c      RN=DMAX1(RN,RMAX)
c      HFT=DLOG(RN*ZATOM*SCALE)
c      NHP=1.D0/DLOG(RN/(RN-DP))+1.
c      HP=DFLOAT(NHP)
c      HP=DMAX1(HP,HMIN)
c      N=HFT*HP
c      IF(N.LT.NMWF) GO TO 7000
c      HSCALE=DFLOAT(N)/DFLOAT(NMWF)
c      WRITE(6,6005) N,HP
c      HP=HP/HSCALE
c      HP=DMAX1(HP,32.D0)
c      GO TO 7000
c 3700 CONTINUE
c      FNN=-(XION+1.)**2/E
c      FNN=DSQRT(FNN)
c      FNN=DMAX1(FNN,1.D0)
c      FNN=DMIN1(FNN,FN)
c gaoxiang 2012.10.3
      fnn=fn
c gaoxiang 2012.10.3
      FACT=1.D0+0.5D0*DMAX1(0.D0,(FNN-3.D0))
      FACT=FACT*(1.D0+0.5D0*DMAX1(0.D0,(FL-1.D0)))
      FACT=FACT*(15.D0+12.D0*DMAX1(1.D0,XION))
      RFACT=RFAC
      IF(RFA.GT.0.D0) RFACT=RFA
      RFT=RFACT*FACT
      IF (FN.EQ.1) RFT=2.5D0*RFT
      RN=RFT/(XION+1.D0)
      HFT=DLOG(RN*ZATOM*SCALE)
      HP=HMIN
      N=HP*HFT
      IF (FN.EQ.1.D0) GOTO 7000
      N=N+(FN-1.)*NSTEP
      RN=RN*FN*FN
c gaoxiang 2012.10.3 redefine the value of rn
      rn=DMIN1(RN,RMAX)
c gaoxiang 2012.10.3 redefine the value of rn
      HFT=HFT+2.D0*DLOG(FN)
      NHP=N/HFT
      IF (NHP.GE.NMIN) GO TO 3800
      NHP=NMIN
      N=NHP*HFT
 3800 HP=DFLOAT(NHP)
      IF(N.LT.NMWF) GO TO 7000
      HSCALE=N/NMWF*(1.D0+200.D0/E)
c gaoxiang 2012.10.3 may have problem of hscale!
      WRITE(6,6005) N,HP
      HP=HP/HSCALE
 6005 FORMAT(' **************  N=',I5,'  H=',E10.3,' **************')
      HP=DMAX1(HP,32.D0)
 7000 H=1.D0/HP
      N=(N/4+1)*4
      N=MIN0(NMWF,N)
c      N=4000
      IF (NDEBUG.GE.1) WRITE (6,6002) N,RN,HP,H
 6002 FORMAT(I5,1P9E13.6)
      EPS=RN*2
      D=DEXP(H)
      R=RN/D**N
c      IF(E.LT.0.D0) GO TO 7100
c      RTEST=R*D
c      RTEST1=RN*(1.D0-1.D0/D)
c      IF(RTEST.LE.RUNIT) GO TO 7100
c      IF(RTEST.LE.DP) GO TO 7100
c      IF(RTEST1.LE.DP) GO TO 7100
c      WRITE(6,7101) RTEST,RUNIT,DP
c 7101 FORMAT(' Caution !!! '/,' First step ',d12.4,' should be less'
c     x ' than ',d12.4,' and ',d12.4)
c      WRITE(6,7102) RTEST1,DP
c 7102 FORMAT(' Caution !!! '/,' Last step ',d12.4,' should be less'
c     x ' than ',d12.4)
 7100 DO 5 I=1,N
      R=R*D
      RAD(I)=R
      A(I)=0.
      B(I)=0.
    5 CONTINUE
      NK=2
      DO 24 I=1,N
      R=RAD(I)
   10 IF(R.GE.RH(NB)) GOTO 20
      IF (R.LE.RH(NK)) GO TO 15
      NK=NK+1
      GO TO 10
   15 NKK=NK+1
      IF (NKK.GT.NB) NKK=NB-1
      VR(I)=CUBINT(R,NKK,RH,VRX,NB,NDEBUG)
      GO TO 24
   20 VR(I)=VRX(NB)
   24 CONTINUE
      IF (NDEBUG.LT.2) GO TO 25
      WRITE (6,6001) (RAD(I),I=1,N)
      WRITE (6,6001) (VR(I),I=1,N)
   25 CONTINUE
      EMAX=0
      EMIN=0
c      IF (E.GE.0.) GO TO 30
c gaoxiang 2012.10.3 reset the initia value of orbital energy
      bcond=1
c	print*, 'bcond=',bcond,rn,rmax
      if (rn.ne.rmax) then
        EMAX=-(XION+1.D0)**2/FN**2
        EMAX=EMAX/3.
        EMIN=-(ZN/FN)**2
	  bcond=0
	else
	  emin=esave
        emin=-(zn)**2
        emax=(fn*3.14159265/rn)**2
      endif
   30 CONTINUE
      CALL DIFFER(DELL,EMAX,EMIN)
c gaoxiang 2012.10.3 save the orbital energy for estimating the next one
      esave=e
c gaoxiang 2012.10.3
      XION1=XION+1
      S=2.D0*(FL-FJ)
c gaoxiang 2012.10.6 note the s above is just opposite as defined but is the same as the one in subroutine MATCH!
      DO 35 I=1,N
  35  B(I)=S*B(I)
c   40 IF(E.GT.0.D0) GOTO 51
c      PHI=RQD(E,FN,FJ,XION1)
ccccccc------P.H.Z.---------------
c       WMU(IW)=PHI
c-----------------------------
CCCC   DETERMINATION OF RCUT FOR BOUND STATES           LIU  11/4/91
c      EPSS=2.0D-36 * DSQRT(-2.D0*E)
c      DO 45 I=1,N
c      DNST=A(N-I+1)*A(N-I+1)+B(N-I+1)*B(N-I+1)
c      IF(DNST.LT.EPSS) GOTO 45
c      AAA=DFLOAT(1-I)
c      RN=RN*DEXP(AAA/HP)
c      N=N-I+1
c      GOTO 50
c   45 CONTINUE
ccccc---------P.H.Z.---------------------
c   51 continue
c      WMU(IW)=PHI
c      PHI=CON(WMU,IPHZ(JPHZ),IW,PHI)
c      WMU(IW)=PHI
CC------------------------------------
c   50 WRITE(6,6004) IW,NWFO
c 6004 FORMAT(' *****',I5,' *****',88X,I5)
CCCC                                                   11/4/91
c      WRITE (6,6002) N,HP,RN,E,FN,FL,FJ,PHI,XION1,ZATOM
c      WRITE (2,*) N,HP,RN,E,FN,FL,FJ,PHI,XION1,ZATOM
       write(22)N,HP,RN,E,FN,FL,FJ,PHI,XION1,ZATOM
c	WRITE (2,*) '************************'
      WRITE (6,6002) N,HP,RN,E,FN,FL,FJ,XION1,ZATOM
c gaoxiang 2012.10.6
      write(13,'(4f12.6)') fn,fl,fj,e
      WRITE (2,*) N,HP,RN,E,FN,FL,FJ,phi,XION1,ZATOM
	WRITE (2,*) '************************'
	DO I=1,N
C      WRITE (2,*) (A(I),I=1,N),(B(I),I=1,N)
      WRITE (2,*)RAD(I),A(I),B(I)
	ENDDO
      WRITE (22) (A(I),I=1,N),(B(I),I=1,N)

      IF(NDEBUG.LE.1) GOTO 101
      WRITE (6,6003) (A(I),I=1,N)
      WRITE (6,6003) (B(I),I=1,N)

101   continue

      IF(NDEBUG.eq.0) then
      ifn=nint(fn)
      ifl=nint(fl)
      write(fname,145)ifn,ifl,fj
      write(6,*)fname
      open(unit=888,file=fname,status='unknown')
      write(888,*)(RAD(I),A(I),B(I),i=1,n)
  145  format ('state',i3,i2,f4.1)
      close(888)
	endif
! 6003 FORMAT (1P10E12.4)
C*****************************HANXY.050616************************
	WRITE (10,*) '************************************************'
	write(10,'(A6,F5.1,A6,F5.1,A6,F5.1)')
     & '    N=',FN,'    L=',FL,'    J=',FJ
	DELTN=RAD(N)-RAD(N-1)
	RMAXHAN=RAD(N)
	HMIDD=RMAXHAN/(RMAXHAN-DELTN)
	HANN=DLOG(HMIDD)
	NNHAN=N
C******************DEFINE THE Number OF THE GRID***************
  191	NNHAN=NNHAN+100
  	RNTHANN=RMAXHAN/((RMAXHAN/(RMAXHAN-DELTN))**(NNHAN-1)-1)
	DELT1=RNTHANN*((RMAXHAN/(RMAXHAN-DELTN))-1)
	IF (DELT1.GT.1.03D-8) GOTO 191
C**************************************************************
	RMXGRASP=RNTHANN*(EXP((NNHAN-1)*HANN)-1)
c
	write(10,'(A7,A14,E14.6)')'GRASP: ','RMAX=',RMXGRASP
	write(10,'(A7,A14,E14.6)')'GRASP: ','R(N)-R(N-1)=',
     *   RNTHANN*exp((NNHAN-1)*HANN)*(1-exp(-HANN))
	RMINGRASP=RNTHANN*(exp(HANN)-1)
	write(10,'(A7,A14,E14.6)')'GRASP: ','R(2)-R(1)=',RMINGRASP
c
	write(10,'(A7,A14,E14.6)')'RWF: ','RMAX=',RAD(N)
	write(10,'(A7,A14,E14.6)')'RWF: ','R(N)-R(N-1)=',
     *   RAD(N)-RAD(N-1)
	write(10,'(A7,A14,E14.6)')'RWF: ','R(2)-R(1)=',RAD(2)-RAD(1)
c
	IF (RMINGRASP.GT.(RAD(2)-RAD(1))) THEN
	write(10,*)'THE INITIAL STEP OF GRID IS TOO BIG'
	END IF 
c
	write(10,'(A10,A12,I5)')'FOR GRASP:',' N = ',NNHAN
	write(10,'(A22,E14.6)')'             RNT = ',RNTHANN
	write(10,'(A22,E14.6)')'               H = ', HANN 


  100 CONTINUE
c gaoxiang 2012.10.6
      close(13)

c********************HAN.050610**************************
      IF(RRRHAN.LT.RMXGRASP)THEN
	NNNHAN=NNHAN
	RRRHAN=RMXGRASP
	RDELN=RNTHANN*exp((NNHAN-1)*HANN)*(1-exp(-HANN))
	RDEL1=RMINGRASP
      RNTHANNNNN=RNTHANN
	HANNNNN=HANN 
	END IF
c	DO I=1,N
C      WRITE (2,*) (A(I),I=1,N),(B(I),I=1,N)
c      WRITE (2,*)RAD(I),A(I),B(I)
c	ENDDO
      IF(NDEBUG.LE.1) GOTO 1001
      WRITE (6,6003) (A(I),I=1,N)
      WRITE (6,6003) (B(I),I=1,N)
 6003 FORMAT (1P10E12.4)
 1001 CONTINUE
      write(10,*)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
	write(10,*)'THE FINAL RESULTS'
 	write(10,'(A7,A14,E14.6)')'GRASP: ','RMAX=',RRRHAN
	write(10,'(A7,A14,E14.6)')'GRASP: ','R(N)-R(N-1)=',
     *   RDELN
	write(10,'(A7,A14,E14.6)')'GRASP: ','R(2)-R(1)=',RDEL1    
	write(10,'(A10,A12,I5)')'FOR GRASP:',' N = ',NNNHAN
	write(10,'(A22,E14.6)')'             RNT = ',RNTHANNNNN
	write(10,'(A22,E14.6)')'               H = ', HANNNNN 
	write(11,'(A10,A12,I5)')'FOR GRASP:',' N = ',NNNHAN
	write(11,'(A22,E14.6)')'             RNT = ',RNTHANNNNN
	write(11,'(A22,E14.6)')'               H = ', HANNNNN 

      STOP
      END
      DOUBLE PRECISION FUNCTION ADLINT(Y,N,H,X)
C     SIMPSON INTEGRATION
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(  N),X(  N)
      HT=H/3.D0
      SUM=0.D0
      AA=0.D0
      I=0
      NN=N/2
      DO 10 J=1,NN
      I=I+1
      AB=Y(I)*4.D0
      I=I+1
      AC=Y(I)
      SUM=SUM+(AA+AB+AC)*HT
      X(J)=SUM
   10 AA=AC
      ADLINT=SUM
      RETURN
      END
      DOUBLE PRECISION FUNCTION RQD(E,XN,FJ,XION)
      IMPLICIT REAL*8(A-H,O-Z)
      DATA C/137.03602D0/
      XK= FJ+0.5D0
      GAM=DSQRT(XK*XK-(XION/C)**2)
      RW=E/(C*C)
      RQ1=1.D0+RW
      RN=RW/2.D0
      RQ2=DSQRT(1.D0+RN)
      RQ=XION*RQ1/RQ2
      QD=DSQRT(-0.5D0/E)
      U=QD*RQ+XK-GAM
      RQD=XN-U
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE DIFFER(DELL,EMAX,EMIN)
C  CONTROL SUBROUTINE FOR THE INTEGRATION OF THE DIFFERENTIAL EQUATIONS.
C  FINDS EIGENVALUES. NORMALIZES ORBITAL FUNCTIONS.
C  CALLED BY HEX AND BIND.
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DIMENSION ZZX(NMWF)
      DATA C/137.03602D0/
      IF(FN.GT.0.D0)
     XDELK=0.01/FN**3
c      BCOND=0.
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      H3=H/3.0D0
    1 CALL START1
      CALL DIFF1 (KI,FN1)
      IF(NDEBUG.GE.5)
     1WRITE (6,6000) EMIN,E,EMAX,FN,FN1,FL,FJ
 6000 FORMAT (3D16.8,4F10.4)
C  FN1 = NUMBER OF NODES + L + 1. FN1 SHOULD EQUAL FN.
c      IF (E.GE.0.D0) GO TO 100
      IF (FN-FN1) 2,4,13
C  TOO MANY NODES.
    2 IF (E.LT.EMAX) EMAX=E
      E=0.5D0*(E+EMIN)
      DL1=DABS(EMAX-EMIN)/(DABS(E)+1.D0)
      IF (DL1.GE.DELK) GO TO 1
      WRITE (6,15) FN1,FN,FL,FJ,EMIN,E,EMAX
   15 FORMAT (16H0TOO MANY NODES.,4F4.1,3D16.8)
c gaoxiang 2012.10.3
      if(emin.lt.0.d0) then
        EMIN=1.2D0*EMIN
	else
	  emin=0.8*emin
	endif
c gaoxiang 2012.10.3
      GO TO 1
C  CORRECT NUMBER OF NODES.
    4 AKI=A(KI)
      BKI=B(KI)
      CALL START2 (KI,KJ)
      CALL DIFF2 (KI,KJ)
      AIK=A(KI)
      BIK=B(KI)
      RA=AKI/AIK
      RB=BKI/BIK
      DO 5 K=KI,KJ
      A(K)=A(K)*RA
    5 B(K)=B(K)*RA
      DG=DEXP(H*G)
      RG=RN**G/DG**N
      DO 6 K=1,KJ
      RG=RG*DG
      A(K)=A(K)*RG
    6 B(K)=B(K)*RG
      IF (KJ.EQ.N) GO TO 8
      RJ=RAD(KJ)
      DO 7 K=KJ,N
      R=RAD(K)
c      W=DEXP(-DSQRT(-2.0D0*E)*(R-RJ))
c gaoxiang 2012.10.3
c      appp = DSQRT(-2.0D0*E)*(R-RJ)
      appp = DSQRT(dabs(-2.0D0*E))*(R-RJ)
      w = 0.d0
      if(appp.lt.150.d0) w = dexp(-appp)
      A(K)=A(KJ)*W
    7 B(K)=B(KJ)*W
    8 DO 9 K=1,N
      R=RAD(K)
    9 Y(K)=(A(K)*A(K)+B(K)*B(K))*R
      W=ADLINT(Y,N,H,ZZX)
      W=W+Y(1)/(G+G+1.0D0)
      DE=CS*A(KI)*B(KI)*(1.0D0-RB/RA)/W
      DL1=DABS(EMAX-EMIN)/(DABS(E)+1.D0)
      DL=DABS(DE/E)
      IF(NDEBUG.GE.5) WRITE(6,6000) E,DE,DL1
      IF ((DL.GT.DELK) .AND. (DL1.LT.DELK)) GO TO 12
      IF (DE.GT.0.0) EMIN=E
      IF (DE.LT.0.0) EMAX=E
      E=E+DE
      IF (DL.LT.DELL) GO TO 10
      IF (E.GT.EMAX) E=0.5D0*(E-DE+EMAX)
      IF (E.LT.EMIN) E=0.5D0*(E-DE+EMIN)
      GO TO 1
   10 X=1.D0/DSQRT(W)
      DO 11 K=1,N
      A(K)=A(K)*X
   11 B(K)=B(K)*X
      RETURN
   12 WRITE (6,16) FN,FL,FJ,EMIN,E,EMAX,DE
   16 FORMAT (16H0 DE TOU  LARGE.,3F4.1,4D15.8)
      GO TO 10
C  TOO FEW NODES.
   13 IF (E.GT.EMIN) EMIN=E
      E=0.5D0*(E+EMAX)
      DL1=DABS(EMAX-EMIN)/(DABS(E)+1.D0)
      IF (DL1.LT.DE LK) GO TO 14
      GO TO 1
   14 WRITE (6,17) FN 1,FN,FL,FJ,EMIN,E,EMAX
   17 FORMAT(' TOO FEW NODES',3X,4F4.1,3D16.8)
c gaoxiang 2012.10.3
      if (emax.lt.0.d0) then
        EMAX=0.8*EMAX
c        EMAX=DMAX1(1.0D0,10.D0/RN**2)
      else
	  emax=1.2*emax
	endif
c gaoxiang 2012.10.3
      E=0.5D0*(E+EMAX)
      GO TO 1
  100 DG=DEXP(H*G)
      RG=RN**G/DG**N
      DO 105 K=1,N
      RG=RG*DG
      A(K)=A(K)*RG
  105 B(K)=B(K)*RG
      CALL MATCH
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE DIFF1(KI,FN1)
C  INTEGRATES OUTWARD TO CLASSICAL TURNING POINT. COUNTS N0DES IN THE
C  MAJOR COMPONENT OF THE RADIAL WAVE FUNCTION.
C  CALL ED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      H3=H/3.0D0
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      Q11=-G+S*(FJ+0.5D0)
      Q22=-G-S*(FJ+0.5D0)
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=DA(4)
      DB(3)=DB(4)
      FN1=1.0D0+FL
      M=N-10
      KM=N-1
c      IF (E.GE.0.) GO TO 2
      DO 1 K=11,M
      KM=N-K
      R=RAD(KM)
      IF (E*R+Z-VR(KM).GT.0.0) GO TO 2
    1 CONTINUE
    2 KI=KM+1
      DO 3 K=5,KI
      R=RAD(K)
      RP21=(E*R-VR(K)+Z)/CS
      RP12=-RP21-2.0D0*CS*R
      A(K)=A(K-4)+8.0D0*(DA(3)-0.5D0*DA(2)+DA(1))
      B(K)=B(K-4)+8.0D0*(DB(3)-0.5D0*DB(2)+DB(1))
      DA(4)=H3*(Q11*A(K)+RP12*B(K))
      DB(4)=H3*(Q22*B(K)+RP21*A(K))
      A(K)=A(K-1)+1.125D0*DA(4)+2.375D0*DA(3)-0.625D0*DA(2)+.125D0*DA(1)
      B(K)=B(K-1)+1.125D0*DB(4)+2.375D0*DB(3)-0.625D0*DB(2)+.125D0*DB(1)
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=H3*(Q11*A(K)+RP12*B(K))
      DB(3)=H3*(Q22*B(K)+RP21*A(K))
   3  FN1=FN1+DABS(DSIGN(1.0D0,A(K))-DSIGN(1.0D0,A(K-1)))/2.D0
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE DIFF2(KI,KJ)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
C   INTEGRATES INWARD TO THE CLASSICAL TURNING POINT
C   CALLDE BY DIFFER
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      H3=-H/3.0D0
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      Q11=-G+S*(FJ+0.5D0)
      Q22=-G-S*(FJ+0.5D0)
      K=KJ-3
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=DA(4)
      DB(3)=DB(4)
    1 K=K-1
      R=RAD(K)
      RP21=(E*R-VR(K)+Z)/CS
      RP12=-RP21-2.0D0*CS*R
      A(K)=A(K+4)+8.0D0*(DA(3)-0.5D0*DA(2)+DA(1))
      B(K)=B(K+4)+8.0D0*(DB(3)-0.5D0*DB(2)+DB(1))
      DA(4)=H3*(Q11*A(K)+RP12*B(K))
      DB(4)=H3*(Q22*B(K)+RP21*A(K))
      A(K)=A(K+1)+1.125D0*DA(4)+2.375D0*DA(3)-0.625D0*DA(2)+.125D0*DA(1)
      B(K)=B(K+1)+1.125D0*DB(4)+2.375D0*DB(3)-0.625D0*DB(2)+.125D0*DB(1)
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DB(2)=DB(3)
      DA(3)=H3*(Q11*A(K)+RP12*B(K))
      DB(3)=H3*(Q22*B(K)+RP21*A(K))
      IF (K.GT.KI) GO TO 1
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE START1
C  CALLED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      H3=H/3.0D0
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      GP=G+S*(FJ+0.5D0)
      GM=G-S*(FJ+0.5D0)
      P21=(E-V0)/CS
      P12=-P21-2.0D0*CS
      Q11=-GM
      Q21=Z/CS
      Q12=-Q21
      Q22=-GP
      A0(1)=DSQRT(S*GP)
      B0(1)=DSQRT(-S*GM)
      FI=0.0
      DO 1 I=2,5
      FI=FI+1.0D0
      W1=P12*B0(I-1)/FI
      W2=P21*A0(I-1)/FI
      A0(I)=W1+(Q11*W1+Q12*W2)/(FI+2.0D0*G)
    1 B0(I)=W2+(Q21*W1+Q22*W2)/(FI+2.0D0*G)
      DO 2 K=1,5
      R=RAD(K)
      A(K)=A0(1)+R*(A0(2)+R*(A0(3)+R*(A0(4)+R*A0(5))))
      B(K)=B0(1)+R*(B0(2)+R*(B0(3)+R*(B0(4)+R*B0(5))))
      DA(K)=H3*(Q11*A(K)+(R*P12+Q12)*B(K))
    2 DB(K)=H3*(Q22*B(K)+(R*P21+Q21)*A(K))
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE START2(KI,KJ)
C  STARTS INWARD INTEGRATION. OUTER BOUNDARY CONDITION IS SET HERE.
C  CALLED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,EEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3Y(NMWF)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
c gaoxiang 2012.10.6 to read in the structure factor to form boundary condition
      common /gaobnd/ skli
c gaoxiang 2012.10.6
      DATA C/137.03602D0/
      H3=-H/3.0D0
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      Q22=-G-S*(FJ+0.5D0)
      Q11=-G+S*(FJ+0.5D0)
c gaoxiang 2012.10.6
      DO 1 K=KI,N
      R=RAD(K)
c      IF ((E*R-VR(K)+Z)*R+EPS) 2,1,1
    1  KJ=K
      IF (BCOND.NE.0.) GO TO 3
C  ISOLATED ATOM BOUNDARY CONDITIONS.
    2 RZ=(Z-VR(KJ))/R
      RL=(FJ+0.5D0)/R
      RK=-S*(FJ+0.5D0)/R
      P=-2.0D0*(E+RZ)+RL*RL
    8 FORMAT (' THE TURNING POINT IS BEYOND RN. A LARGER VALUE OF RN MAY
     1 BE NEEDED.')
      IF (P.LT.0.0) WRITE (6,8)
      P=DSQRT(P)
      APA=-P+0.5D0*(RZ-RL*RL)/(R*P**2)
      A(KJ)=1.0D0
      B(KJ)=CS*(APA+RK)/(-2.0D0*C*C-E-RZ)
      GO TO 4
C  WIGNER-SEITZ BOUNDARY CONDITIONS.
c    3 IF (BCOND.NE.1.) GO TO 23
c      B(KJ)=DMOD(FL,2.D0)
c      A(KJ)=1.0D0-B(KJ)
c      GO TO 40
c   23 IF(BCOND.NE.2.)GO TO 33
c          A(KJ)=DMOD(FL,2.D0)
c      B(KJ)=1.-A(KJ)
c      GO TO 40
c   33 IF(BCOND.NE.3.) GO TO 35
c      A(KJ)=1.
c      B(KJ)=1.
c      GO TO 40
c   35 A(KJ)=1.
c      B(KJ)=-1.
    3 continue
c gaoxiang 2012.10.6 to read in the structure factor to form boundary condition
c 2012.10.6 gaoxiang using the non-relativistic approximation first!
c      B(KJ)=1.d0
c      dle=(fl*skli+2.d0*(fl+1.d0)*(2.d0*fl+1.d0))
c     */(skli-2.d0*(2.d0*fl+1.d0))
c	rkappa=-s*(fj+0.5d0)
c      A(KJ)=-2.d0*C*rn/(dle+rkappa+1.d0)
c	print*, 'rn=',rn,c,fn,fl,fj,rkappa,dle,skli,a(kj)
c   40 CONTINUE
      B(kj)=1.d0
	A(kj)=0.d0

c gaoxiang 2012.10.6
    4 DO 5 L=1,4
      K=KJ-L
      A(K)=A(KJ)
    5 B(K)=B(KJ)
      DO 7 I=1,4
      K=KJ+1
      DO 6 L=1,5
      K=K-1
      R=RAD(K)
      RP21=(E*R-VR(K)+Z)/CS
      RP12=-RP21-2.0D0*CS*R
      DA(L)=H3*(Q11*A(K)+RP12*B(K))
    6 DB(L)=H3*(Q22*B(K)+RP21*A(K))
      A(KJ-1)=A(KJ)+(251.0D0*DA(1)+646.0D0*DA(2)-264.0D0*DA(3)
     1+106.0D0*DA(4)-19.0*DA(5))/240.0D0
      B(KJ-1)=B(KJ)+(251.0D0*DB(1)+646.0D0*DB(2)-264.0D0*DB(3)
     1+106.0D0*DB(4)-19.0*DB(5))/240.0D0
      A(KJ-2)=A(KJ)+(116.0D0*DA(1)+496.0D0*DA(2)+96.0D0*DA(3)
     1+16.0D0*DA(4)-4.0D0*DA(5))/120.0D0
      B(KJ-2)=B(KJ)+(116.0D0*DB(1)+496.0D0*DB(2)+96.0D0*DB(3)
     1+16.0D0*DB(4)-4.0D0*DB(5))/120.0D0
      A(KJ-3)=A(KJ)+(81.0D0*DA(1)+306.0D0*DA(2)+216.0D0*DA(3)
     1+126.0D0*DA(4)-9.0D0*DA(5))/80.0D0
      B(KJ-3)=B(KJ)+(81.0D0*DB(1)+306.0D0*DB(2)+216.0D0*DB(3)
     1+126.0D0*DB(4)-9.0D0*DB(5))/80.0D0
      A(KJ-4)=A(KJ)+(56.0D0*DA(1)+256.0D0*DA(2)+96.0D0*DA(3)
     1+256.0D0*DA(4)+56.0D0*DA(5))/60.0D0
    7 B(KJ-4)=B(KJ)+(56.0D0*DB(1)+256.0D0*DB(2)+96.0D0*DB(3)
     1+256.0D0*DB(4)+56.0D0*DB(5))/60.0D0
      RETURN
      END
C SEGMENT CMMN1+NP,DATA1
      SUBROUTINE RDEN(LOG,IREAD)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(JM),XL(JM),XJ(JM),XE(JM),XZ(JM),BCOND,XNNN,
     2  FN,FL,E,FJ,VR(NMWF),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3 Y(NMWF)
      COMMON/DATA1/
     1 A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5)
     2,RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DIMENSION ALP(10),CDEN(NMX),XRD(NMX)
      EQUIVALENCE (XRD(1),RH(1))
      READ (LOG,1000) ALP
 1000 FORMAT (10A8)
      WRITE (6,6000) ALP
 6000 FORMAT (' *****  INPUT FOR TRIAL POTENTIAL ***** ',/1X,10A8)
      READ (LOG,1001) NP,JP,RNP,HP,ZNP,XIP,PHP,EPP,DEP,DELP,ALP(1),I
 1001 FORMAT (2I5,F6.3,F4.0,2F5.0,4F10.7,2X,A5,I3)
      IF(NP.GE.NMX) STOP ' STOP IN RDEN BECAUSE NP IS TOO LARGER ! '
      IF (NDEBUG.GE.1)
     XWRITE(6,1001) NP,JP,RNP,HP,ZNP,XIP,PHP,EPP,DEP,DELP,ALP(1),I
      READ (LOG,1002) NCYC,ALP(1),I
 1002 FORMAT (I5,67X,A5,I3)
      IF (NDEBUG.GE.1)
     XWRITE (6,1002) NCYC,ALP(1),I
      READ (LOG,1003) KPTP,XALPP,XLATP,XNNP,FEP,RNUP,ALP(1),I
 1003 FORMAT (I10,F10.7,9X,F5.1,F10.5,F5.1,D10.3,13X,A5,I3)
      IF (NDEBUG.GE.1)
     XWRITE (6,1003) KPTP,XALPP,XLATP,XNNP,FEP,RNUP,ALP(1),I
      DO 10 L=1,JP
      READ (LOG,1004) XNP,XLP,XJP,EIGP,OCCP,SIGP,BCNP,ALP(1),I
      IF (NDEBUG.GE.1)
     XWRITE(6,1004) XNP,XLP,XJP,EIGP,OCCP,SIGP,BCNP,ALP(1),I
   10 CONTINUE
 1004 FORMAT(3F5.1,F14.7,F7.4,5X,F9.4,E12.4,10X,A5,I3)
      DO 20 L=1,NP,5
      JL=L+4
      READ (LOG,1005) (CDEN(K),K=L,JL),ALP(1),I
 1005 FORMAT (5E14.7,2X,A5,I3)
      IF (NDEBUG.GE.1)
     XWRITE (6,6001) (CDEN(K),K=L,JL),ALP(1),I
   20 CONTINUE
 6001 FORMAT (1X,1P5E14.7,2X,A5,I3)
      RHO=CDEN(1)
      ZN=ZNP
      XION=XIP
      RNUC=RNUP
      H=1.D0/HP
      IREAD=NP
   60 DO 70 L=1,NP,5
      JL=L+4
      READ (LOG,1005) (CDEN(K),K=L,JL),ALP(1),I
      IF (NDEBUG.GE.1)
     XWRITE (6,6001) (CDEN(K),K=L,JL),ALP(1),I
   70 CONTINUE
      Z=ZN
      IF(RNUC.GT.0.0D0) Z=0.D0
      DO 72 L=1,NP
   72 VRX(L)=    CDEN(L)
      HP=1.D0/HP
      D=DEXP(HP)
      R=RNP/D**NP
      DO 78 L=1,NP
      R=R*D
   78 XRD(L)=R
      A1=XRD(1)*RHO/3.D0
      B1=RHO/2.D0
      Q1=Z-VRX(1)*ZN
      X=XRD(1)/(RNUC+1.D-40)
      Z=ZN
      RVN=0.D0
      IF (RNUC.EQ.0.D0) GO TO 80
      RVN=-Z*X*(1.5D0-0.5D0*X*X)
   80 CONTINUE
      V0=(Q1-RVN-A1+XRD(1)*B1)/XRD(1)
      IF (RNUC.EQ.0.D0) RETURN
      V0=V0-1.5D0*Z/(RNUC+1.D-40)
      RETURN
      END
      DOUBLE PRECISION FUNCTION CUBINT(XJ,N,X,Y,NX,NDEBUG)
C  CUBIC LAGRANGE INTERPOLATI0N TO FIND Y(XI). N IS A LOCATION IN THE
C  TABLE OF X NEAR XI. UXUALLY XI IS BETWEEN X(N-1) AND X(N)
C  CALLED BY POT.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NX), Y(NX)
      IF (N.GE.NX) GO TO 1
      XI=XJ
      AMM=(XI-X(N-1))*(XI-X(N))*(XI-X(N+1))
      AM=(XI-X(N-2))*(XI-X(N))*(XI-X(N+1))
      A=(XI-X(N-2))*(XI-X(N-1))*(XI-X(N+1))
      AP=(XI-X(N-2))*(XI-X(N-1))*(XI-X(N))
      XI=X(N-2)
      BMM=(XI-X(N-1))*(XI-X(N))*(XI-X(N+1))
      XI=X(N-1)
      IF (BMM.EQ.0.0) GO TO 1
      BM=(XI-X(N-2))*(XI-X(N))*(XI-X(N+1))
      XI=X(N)
      IF (BM.EQ.0.0) GO TO 1
      B=(XI-X(N-2))*(XI-X(N-1))*(XI-X(N+1))
      XI=X(N+1)
      IF (B.EQ.0.0) GO TO 1
      BP=(XI-X(N-2))*(XI-X(N-1))*(XI-X(N))
      IF (BP.EQ.0.0) GO TO 1
      CUBINT=Y(N-2)*AMM/BMM+Y(N-1)*AM/BM+Y(N)*A/B+Y(N+1)*AP/BP
      RETURN
    1 CUBINT=Y(NX)
      IF(NDEBUG.LT.4) RETURN
      WRITE(6,2) N,X(N-2),X(N-1),X(N)
    2 FORMAT(' TROUBLE IN CUBINT. SET TO Y(N).',I5,4D16.8)
      RETURN
      END
C SEGMENT CMMN1,DATA1
      SUBROUTINE CRF(KAP,W,Z0,RO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C     ENERGY NORMALIZATION:  ENERGY IN NATURAL UNIT.
C
C
      COMMON/RK/RK(3,4),Y(4),YC(4),DY(4),H,VH,X,EP,EM,DKAP,V,AZ
      COMMON /CWF/G(2),F(2)
C     P*R .GE. 10.D0
      DATA PI/3.14159265358979D0/
      DATA ALPHA/137.0360/
      CON=2.*ALPHA
      DKAP=KAP
      EM=-W/CON
      EP=CON-EM
      P=DSQRT(W*(1.+W/CON**2))
      FAC=W/(CON*P)
      GNUP=Z0/P
      GNU=GNUP*(1.+2.*W/CON**2)
      AZ=Z0/ALPHA
      GAM=DSQRT(KAP**2-AZ**2)
      CVR=1.D0
      CVI=0.D0
      CWR=1.D0
      CWI=0.D0
      CTVR=1.D0
      CTVI=0.D0
      CTWR=1.D0
      CTWI=0.D0
      C3R=GAM
      C3I=-GNU
      C4R=-GAM
      C4I=-GNU
      C1R=C3R-1.D0
      C1I=C3I
      C2R=C4R-1.D0
      C2I=C4I
      RIR=0.D0
      RII=-1.D0
      RR=0.1D0*((GAM+10)**2+GNU**2)/P
      NP=(RR-RO)*P
      IF(NP.LE.0) NP=1
      RR=RO+NP/P
      TPR=2.*P*RR
      DO 11 J=1,100
      TPRJ=TPR*DFLOAT(J)
      C1R=C1R+1.D0
      C2R=C2R+1.D0
      C3R=C3R+1.D0
      C4R=C4R+1.D0
      CALL CPRD(CTVR,CTVI,RIR,RII,DR,DI)
      CALL CPRD(DR,DI,C1R,C1I,CTVR,CTVI)
      CALL CPRD(CTVR,CTVI,C2R,C2I,DR,DI)
      CTVR=DR/TPRJ
      CTVI=DI/TPRJ
      CVR=CVR+CTVR
      CVI=CVI+CTVI
      CALL CPRD(CTWR,CTWI,RIR,RII,DR,DI)
      CALL CPRD(DR,DI,C3R,C3I,CTWR,CTWI)
      CALL CPRD(CTWR,CTWI,C4R,C4I,DR,DI)
      CTWR=DR/TPRJ
      CTWI=DI/TPRJ
      CWR=CWR+CTWR
      CWI=CWI+CTWI
      T1=DSQRT(CTVR*CTVR+CTVI*CTVI)
      T2=DSQRT(CTWR*CTWR+CTWI*CTWI)
      TEST=DMAX1(T1,T2)
      D1=DSQRT(CVR*CVR+CVI*CVI)
      D2=DSQRT(CWR*CWR+CWI*CWI)
      DEST=DMIN1(D1,D2)
      TEST=TEST/DEST
      IF(TEST.LE.1.D-6) GO TO 12
   11 CONTINUE
      WRITE(6,6001) TEST,DEST
 6001 FORMAT(' ***** CONVERGENCE IN FOLLOWING C WFS ',1P2E15.7)
   12 GN=-GNUP
      CALL CPRD(GN,DKAP,CWR,CWI,DR,DI)
      CWR=DR/TPR
      CWI=DI/TPR
      GR=CVR+CWR
      GI=CVI+CWI
      FR=CVR-CWR
      FI=CVI-CWI
      ARG=0.5D0*TPR+GNU*DLOG(TPR)
      XDK=DABS(DKAP)
      T1=DATAN(GNUP/XDK)
      IF(DKAP.GT.0.D0) T1=-PI-T1
      T2=DATAN(GNU/GAM)
      GN=-GNU
      DCP=GAMLGI(GAM,GN)
      DELC=T1-T2-PI*GAM+2.D0*DCP
      DELC=0.5D0*DELC
      ARG=ARG+DELC
      EC=DCOS(ARG)
      ES=DSIN(ARG)
      CALL CPRD(EC,ES,GR,GI,DR,DI)
      G(1)=DR
      G(2)= DI
      CALL CPRD(EC,ES,FR,FI,DR,DI)
      F(1)=-DI*FAC
      F(2)= DR*FAC
      XN=G(1)*F(2)-G(2)*F(1)
      XN=XN*PI
      XN=DSQRT(DABS(XN))
      DO 20  I=1,2
      G(I)=G(I)/XN
   20 F(I)=F(I)/XN
      N8=16
      H=-1.D0/(DFLOAT(N8)*P)
      VH=0.5D0*H
      X=RR
      ND= DFLOAT(N8)*NP
      Y(1)=G(1)
      Y(2)=F(1)
      Y(3)=G(2)
      Y(4)=F(2)
      V=-AZ/X
      CALL SLOPE
      DO 30 I=1,ND
      CALL RUNGE
      XN=Y(1)*Y(4)-Y(2)*Y(3)
   30 CONTINUE
      G(1)=Y(1)
      G(2)=Y(3)
      F(1)=Y(2)
      F(2)=Y(4)
      RETURN
      END
      SUBROUTINE CPRD(AR,AI,BR,BI,PR,PI)
      IMPLICIT REAL*8(A-H,O-Z)
      PR=AR*BR-AI*BI
      PI=AR*BI+AI*BR
      RETURN
      END
      DOUBLE PRECISION FUNCTION GAMLGI(R,AI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(4)
      DATA C/12.D0,-360.D0,1260.D0,-1680.D0/
      UR=R+10.D0
      UI=AI
      CALL ICOM(UR,UI,VR,VI)
      CALL CPRD(VR,VI,VR,VI,WR,WI)
      XR=0.5D0*DLOG(UR*UR+UI*UI)
      XI=DATAN2(UI,UR)
      U5=UR-0.5D0
      CALL CPRD(U5,UI,XR,XI,SR,SI)
      SI=SI-UI
      DO 5 K=1,4
      T=VI/C(K)
      CALL CPRD(VR,VI,WR,WI,VR,VI)
    5 SI=SI+T
      UR=R
      DO 10 K=1,10
      XI=DATAN2(UI,UR)
      SI=SI-XI
   10 UR=UR+1.D0
   15 GAMLGI=SI
      RETURN
      END
      SUBROUTINE ICOM(XR,XI,YR,YI)
      IMPLICIT REAL*8 (A-H,O-Z)
      D=XR*XR+XI*XI
      YR=XR/D
      YI=-XI/D
      RETURN
      END
C SEGMENT DATA1
      SUBROUTINE RUNGE
C     4-TH ORDER RUNGE-KUTTA METHOD
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/RK/RK(3,4),Y(4),YC(4),DY(4),H,VH,X,EP,EM,CAP,V,Z
C     H STEP SIZE   VH=0.5D0*H
C     X POSO ITION
C     EP   M+E   137*2+KE/137
C     EM   M-E
C     V   POTENTIAL
C     Z   ZI/137
      DO 1 J=1,4
    1 YC(J)=Y(J)
      DO 2 J=1,4
      RK(1,J)=H*DY(J)
    2 Y(J)=YC(J)+0.5D0*RK(1,J)
      X=X+VH
      V=-Z/X
      CALL SLOPE
      DO 3 J=1,4
      RK(2,J)=H*DY(J)
    3 Y(J)=YC(J)+0.5D0*RK(2,J)
      CALL SLOPE
      DO 4 J=1,4
      RK(3,J)=H*DY(J)
    4 Y(J)=YC(J)+RK(3,J)
      X=X+VH
      V=-Z/X
      CALL SLOPE
      DO 5 J=1,4
    5 Y(J)=YC(J)+(RK(1,J)+2.D0*(RK(2,J)+RK(3,J))+H*DY(J))/6.D0
      RETURN
      END
C SEGMENT DATA1
      SUBROUTINE SLOPE
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/RK/RK(3,4),Y(4),YC(4),DY(4),H,VH,X,EP,EM,CAP,V,Z
      EPP=EP-V
      EMM=EM+V
      CAPX=CAP/X
      DY(1)=EPP*Y(2)-CAPX*Y(1)
      DY(2)=EMM*Y(1)+CAPX*Y(2)
      DY(3)=EPP*Y(4)-CAPX*Y(3)
      DY(4)=EMM*Y(3)+CAPX*Y(4)
      RETURN
      END
C SEGMENT CMMN1,DATA1
      SUBROUTINE MATCH
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,JM=35)
      COMMON/DATA1/
     3A(NMWF),B(NMWF),RAD(NMWF),VRX(NMX),XKOOP,DA(5),DB(5),A0(5),B0(5),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
c      COMMON/DATA1/A(NMWF),B(NMWF),RAD(NMWF),DY(863),KP,N,J,NC,NDEBUG,NP
      COMMON/CMMN1/D1(14),XION,PHI,D2(185),FL,E,FJ,D3(4012)
      COMMON /CWF /G(2),F(2)
      DATA PI/3.141592654/
      ZI=XION+1.D0
      RR=RAD(N)
      S=2.D0*(FL-FJ)
      KAP=S*(FJ+0.5D0)
      ER=2.D0*E
      CALL CRF(KAP,ER,ZI,RR)
      DOM=G(1)*F(2)-F(1)*G(2)
      XC=(A(N)*F(2)-S*B(N)*G(2))/DOM
      XS= (A(N)*F(1)-S*B(N)*G(1))/DOM
      PHI=DATAN2(XS,XC)/PI
      IF(PHI.LT.0.D0) PHI=PHI+2.D0
      XN=DSQRT((XS*XS+XC*XC)*137.036D0)
      XN=1.D0/XN
   80 DO 100 I=1,N
      A(I)=A(I)*XN
  100 B(I)=B(I)*XN
      IF(DABS(DOM-0.3183099D0).LE.0.001D0) RETURN
      WRITE (6,6001) DOM,XC,XS,C,XN
 6001 FORMAT(' *****',1P8E12.4)
      RETURN
      END
cc---------------P.H.Z.------------------
      DOUBLE PRECISION FUNCTION CON(WMU,NW,IW,PHI)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMX=1000,NMWF=12000,jm=35,MAXWF=350)
      dimension WMU(MAXWF)
      INTEGER NW,IW
      CON=PHI
      IF(NW.LE.1) RETURN
      NN=IW-1
      X=WMU(NN)-PHI
      IF (X.EQ.0.D0) RETURN
      XX=DABS(X)
      S=X/XX
      N=XX+0.2D0
      N=(N/2)*2
      CON=PHI+S*DFLOAT(N)
      RETURN
      END

