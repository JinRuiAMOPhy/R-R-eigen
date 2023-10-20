      IMPLICIT REAL*8(A-H,O-Z)
C
C   ---- This program is modified by Xiao-Min Tong at Arpil 25, 1994
C       (1). extend the dimension of R from 441 to NMX = 1000
C       (2). parameterlized the jm = 35
C

      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   VINDEN(jm),DA(5),DV(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      COMMON /NP/ RNK,ANK,VNK,NK
      DIMENSION ERR(5),RHO(NMX),UR(NMX)
      DATA C/137.03602D0/
C
C --- file allocate:  1 : input try potential;  15 : input card contral
C                     9 : rscf potential output. 50 :vscf potential output
C
      open(1,file='try.pot',status='unknown')
      open(15,file='rscf.in',status='old')
      open(4,form='unformatted',status='scratch')
      open(9,file='radata.dat',status='unknown')
      open(50,file='vscf.dat',status='unknown')
      RHO(:) = ZERO
      NCASE=0
    1 CONVRG=0.5D0
      NCASE=NCASE+1
      ET4=0.5D0
      XNNN=1.0D0
      ET3=0.0
      CALL READIN(IREAD,BCON)
      D=DEXP(H)
      R=RN/D**N

      DO 2 K=1,N
      R=R*D
      RAD(K)=R
      A(K)=0.0
      B(K)=0.0
      DENS(K)=0.0
      RHNLJ(K)=0.0
      VR(K)=0.0
      VRX(K)=0.0
      UR(K)=0.
    2 Y(K)=0.0
      NUM1=0.0
      IF(NCASE.EQ.1) CALL START0
      ERR(1)=0.1D21
      ERR(2)=ERR(1)
      ERR(3)=ERR(1)
      ERR(4)=ERR(1)
      ERR(5)=ERR(1)
      CONVRG=ERR(5)
   21 FORMAT (1H1,10A8)
      WRITE (6,21) AAME
      NCYCLE=NC1
      PH=PHI
C  TAPE 4 IS FOR TEMPORARY STORAGE OF THE RADIAL WAVE FUNCTIONS.
      IF (IREAD.EQ.0) GO TO 3
      LOG=1
      IF((ZN-XION).NE.1.D0) CALL RDEN(LOG,IREAD)
      IF(IREAD.GE.5) GO TO 8
      VR(1)=0.
    3 REWIND 4
    6 XKOOP=0.0
      ERR(1)=ERR(2)
      ERR(2)=ERR(3)
      ERR(3)=ERR(4)
      ERR(4)=ERR(5)
      IF((ZN-XION).NE.1.D0) CALL POT (RH,VR)
c      IF(NDEBUG.GE.2)
c      write(6,*)'zn-xion',zn-xion,zn,xion
c      WRITE(6,6001) (VR(I),I=1,N)
c      write(6,*)'rh'
c      WRITE(6,6001) (RH(I),I=1,N)
c	write(6,*)'here is vscf',n
c	do i=1,n
c      WRITE(6,6001) RH(I),VR(I)
c	enddo

 6001 FORMAT(1P10E12.4)
      R1=RNK
      ERR(5)=CONVRG
      IF(NDEBUG.NE.0)
     1WRITE (6,6000) ERR
 6000 FORMAT(5E16.8)
      ERRM=DMIN1(ERR(1),ERR(2),ERR(3),ERR(4),ERR(5))
      IF (PHI.LT.0.0) GO TO 8
      IF ((NC1-NCYCLE).LT.10) GO TO 8
      IF (ERRM.NE.ERR(1)) GO TO 8
      WRITE (6,22)
   22 FORMAT (//' ERROR IN POTENTIAL IS NO LONGER DECREASING. ITERATION'
     1' TERMINATED ALTHOUGH CONVERGENCE CRITERION IS NOT MET.'//)
      CONVRG=1.0D0
    8 IF (CONVRG.LT.DELRVR) CONVRG=1.0D0
      SUMMA=0.0
      PVT=0.0
      DO 9 K=1,N
    9 DENS(K)=0.0
      REWIND 4
      DELM=0.01D0
      IF(DABS(ET3-ET4).LT.0.1D0) DELM=0.001D0
      IF(DABS(ET3-ET4).LT.0.001D0)DELM=DEL
      DELL=DMIN1(DELM,DEL*DSQRT(CONVRG/DELRVR))
      IF (CONVRG.EQ.1.0D0) DELL=DEL
C  THE USE OF A CONVERGENCE CRITERION, DELL, GREATER THAN DEL MAY LEAD
C  TO INSTABILITY. REPLACE IT WITH DELL=DEL IF NECESSARY
      DO 13 I=1,J
      BCOND=0.0
      IF (XN(I).LT.0.) BCOND=BCON
      FN=DABS(XN(I))
      FL=XL(I)
      FJ=XJ(I)
      E=XE(I)
      IF(NDEBUG.GE.4)
     1WRITE(6,6001) FN,FL,FJ,XE(I)
      CALL DIFFER(DELL)
      IF(NDEBUG.GE.3)
     1WRITE(6,6001) FN,FL,FJ,E
      DO 10 K=1,N
   10 RHNLJ(K)=A(K)**2+B(K)**2
      SUMMA=SUMMA+XZ(I)*E
      IF (CONVRG.EQ.0.0) GO TO 11
      WRITE (4) (A(K),K=1,N)
      WRITE (4) (B(K),K=1,N)
   11 XE(I)=E
      DO 12 K=1,N
   12 DENS(K )=DENS(K)+XZ(I)*RHNLJ(K)
   13 CONTINUE
      IF((ZN-XION).EQ.1.D0) CONVRG=1.D0
      IF((ZN-XION).EQ.1.D0) GOTO 200
C  IF PHI IS NEGATIVE THE CHARGE DENSITY AVERAGING RATIO IS ABS(PHI).
C  IF PHI IS POSITIVE THE RATIO STARTS AT 1.0 AND DECREASES T0 THI.
      IF (PHI.GT.0.0) GO TO 200
      PH=DABS(PHI)
   16 PS=1.0D0-PH
      DO 17 K=1,N
   17 RH(K)=PS*RH(K)+PH*DENS(K)
      GO TO 300
  200 CALL FEEDBK(N,NUM1,RHO,UR,RH,DENS,A,B,ERROR,Y)
  300 CONTINUE
      ET4=ET3
      ET3=SUMMA-EV3+EV-EE2
      IF (RNUC.GT.0.0) ET3=ET3+EV4
      BNK = 0.d0 ! rui jin
      WRITE (6,23) RNK,ANK,BNK,Q,EE2,EV,EV3,SUMMA,ET3
   23 FORMAT(' RNK=',F8.5,'ANK=',F5.1,'BNK=',F11.5,'Q=',F10.5,'EE2=',
     1F10.4,'EV=',F11.4,' EV3=',F11.4,' SUMMA=',F11.4,'ET3=',F11.4)
   18 CONTINUE
      IF (CONVRG.EQ.1.0D0) CALL OUTPUT
      IF (CONVRG.EQ.1.0D0) CALL PUNCHC(4)
      IF (CONVRG.EQ.1.0D0) GO TO 1
      NCYCLE=NCYCLE-1
      IF (NCYCLE.NE.0) GO TO 3
      WRITE (6,19) NC1
   19 FORMAT(//' PROBLEM NOT CONVERGED AFTER',I4,'CYCLES'//)
      CONVRG=1.0D0
      GO TO 18
      END
C SEGMENT CMMN1+NP
      SUBROUTINE RDEN(LOG,IREAD)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DIMENSION ALP(10),CDEN(nmx),XRD(NMX)
      READ (LOG,1000) ALP
      WRITE (6,6000) ALP
      READ (LOG,1001) NP,JP,RNP,HP,ZNP,XIP,PHP,EPP,DEP,DELP,ALP(1),I
      WRITE(6,1001) NP,JP,RNP,HP,ZNP,XIP,PHP,EPP,DEP,DELP,ALP(1),I
      READ (LOG,1002) NCYC,ALP(1),I
      WRITE (6,1002) NCYC,ALP(1),I
      READ (LOG,1003) KPTP,XALPP,XLATP,XNNP,FEP,RNUP,ALP(1),I
      WRITE (6,1003) KPTP,XALPP,XLATP,XNNP,FEP,RNUP,ALP(1),I
      DO 10 L=1,JP
      READ (LOG,1004) XNP,XLP,XJP,EIGP,OCCP,SIGP,BCNP,ALP(1),I
      IF(L.GT.J) GO TO 10
      IF (IREAD.EQ.2) GO TO 5
      IF(IREAD.EQ.3) GO TO 5
      IF (IREAD.EQ.6) GO TO 5
      GO TO 10
    5 CONTINUE
      XN(L)=XNP
      IF(BCOND.NE.0.) XN(L)=-DABS(XN(L))
      XJ(L)=XJP
      XL(L)=XLP
      XE(L)=EIGP
      IF((2.*XJ(L)+1.).LT.XZ(L)) GO TO 200
   10 WRITE(6,1004) XNP,XLP,XJP,EIGP,OCCP,SIGP,BCNP,ALP(1),I
 1004 FORMAT(3F5.1,F14.7,F7.4,5X,F9.4,E12.4,10X,A5,I3)
      DO 20 L=1,NP,5
      JL=L+4
      READ (LOG,1005) (CDEN(K),K=L,JL),ALP(1),I
 1005 FORMAT (5E14.7,2X,A5,I3)
   20 WRITE (6,6001) (CDEN(K),K=L,JL),ALP(1),I
 6001 FORMAT (1X,1P5E14.7,2X,A5,I3)
      GE=(ZN-XION)/(ZNP-XIP)
      IF (RNP.NE.RN) GO TO 30
      HPP=1./H
      IF (HP.NE.HPP) GO TO 30
      DO 25 L=1,NP
   25 RH(L)=CDEN(L)*GE
      GO TO 60
   30 HP=1./HP
      D=DEXP(HP)
      R=RNP/D**NP
      DO 35 L=1,NP
      R=R*D
   35 XRD(L)=R
      NK=2
      DO 50 L=1,N
      R=RAD(L)
   40 IF(R.LE.XRD(NK)) GO TO 45
      NK=NK+1
      GO TO 40
   45 NKK=NK+1
   50 RH(L)=CUBINT(R,NKK,XRD,CDEN,NP)*GE
      DO 52 L=1,N
      IF(RH(L).GE.0.) GO TO 52
      RH(L)=-RH(L)
   52 CONTINUE
      HP=1./HP
 6002 FORMAT  ('  ***** INTERPOLATE SUCCESSFULLY ***** ')
      DO 55 L=1,N
   55 A(L)=RH(L)*RAD(L)
      W=ADLINT(A,N,H)
      GE=(ZN-XION)/W
      DO 58 L=1,N
   58 RH(L)=GE*RH(L)
      WRITE(6,6003) W
      WRITE(6,6003) (RH(I),I=1,N)
      WRITE (6,6002)
 6003 FORMAT(1P10E12.4)
   60 DO 70 L=1,NP,5
      JL=L+4
      READ (LOG,1005) (CDEN(K),K=L,JL),ALP(1),I
   70 WRITE (6,6001) (CDEN(K),K=L,JL),ALP(1),I
      Z=ZN
      IF(RNUC.GT.0.0) Z=0
      IF (RNP.NE.RN) GO TO 75
      HPP=1./H
      IF (HP.NE.HPP) GO TO 75
      DO 72 L=1,NP
   72 VR(L)=Z-ZN*CDEN(L)
      RETURN
   75 HP=1./HP
      D=DEXP(HP)
      R=RNP/D**NP
      DO 78 L=1,NP
      R=R*D
   78 XRD(L)=R
      NK=2
      DO 80 L=1,N
      R=RAD(L)
   90 IF(R.LE.XRD(NK)) GO TO 100
      NK=NK+1
      GO TO 90
  100 NKK=NK+1
      VR(L)=CUBINT(R,NKK,XRD,CDEN,NP)
      VR(L)=Z-ZN*VR(L)
   80 CONTINUE
      WRITE (6,6003) (VR(L),L=1,N)
      WRITE(6,6002)
      RETURN
  200 WRITE (6,6004)
 1000 FORMAT (10A8)
 6000 FORMAT (' *****  INPUT FOR TRIAL POTENTIAL ***** ',/1X,10A8)
 1001 FORMAT (2I5,F6.3,F4.0,2F5.0,4F10.7,2X,A5,I3)
 1002 FORMAT (I5,67X,A5,I3)
 1003 FORMAT (I10,F10.7,9X,F5.1,F10.5,F5.1,D10.3,13X,A5,I3)
 6004 FORMAT('  INPUT IS NOT CONSISTENT WITH TAPE DATA; NEED SUBSHELL IN
     XPUT')
      STOP
      END
C SEGMENT CMMN1+NP
      SUBROUTINE READIN(IREAD,BCON)
      IMPLICIT REAL*8(A-H,O-Z),integer*4 (i-n)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PN(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      ilog = 15
    1 SKIP=0.0
   11 FORMAT (10A8)
      READ (ilog,11,END=100) AAME
   12 FORMAT (2I5,F6.3,F4.3,2F5.0,4F10.7)
      READ (ilog,12) N,J,RN,H,ZN,XION,PHI,EPS,DEL,DELRVR
      IF ((J.LE.0) .OR. (ZN.LE.0.0)) STOP
      IF (N.EQ.0) N=421
      IF (RN.EQ.0.0) RN=60.0D0
      IF (H.EQ.0.0) H=32.0D0
      IF (PHI.EQ.0.0) PHI=0.3D0
      IF (EPS.LT.RN) EPS=75.0D0
      IF (DEL.EQ.0.0) DEL=0.000005D0
      IF (DELRVR.EQ.0.0) DELRVR=0.00001D0
C  THE ABOVE VALUES OF N, RN, H, ETC. SUFFICE FOR MOST ATOMS AND IONS.
      READ(ilog,13) NC1,IREAD,NDEBUG,NPRINT
   13 FORMAT(16I5)
      IF (NC1.EQ.0) NC1=30
   14 FORMAT (I10,F10.7,10X,2F5.1,D10.3)
      READ (ilog,14) KPOT,XALPHA,XLATTR,FEX,RNUC
   18 FORMAT (51X,'   TABLE OF INITIAL DATA'/1H ,10A8/1H ,2I5,F6.2,3F4.0
     1  ,2F8.3,2F10.7/1H ,I5,10X,2I5/1H,I10,F10.7,10X,2F5.1,D10.3)
      WRITE (6,18) AAME, N,J,RN,H,ZN,XION,PHI,EPS,DEL,
     1   DELRVR, NC1,     NDEBUG,NPRINT, KPOT,XALPHA,XLATTR,FEX,RNUC
      BCOND=0.
      H=1.0D0/H
      Z=ZN
      IF (RNUC.EQ.0.0) GO TO 3
C  IF RNUC IS SET TO A VALUE GREATER THAN 1.0 IT IS ASSUMED THAT THIS IS
C  THE NUCLEAR MASS NUMBER, A, AND THE NUCLEAR RADIUS IS COMPUTED FROM
C  THE FORMULA R=1.1*A**(1/3) FERMIS.
      IF (RNUC.GE.1.0D0) RNUC=0.0000208D0*RNUC**(1.0D0/3.0D0)
      Z=0.0
      R=RN*DEXP(-(N-3)*H)
      IF (RNUC.GE.R) GO TO 3
      WRITE (6,22)
   22 FORMAT (' LESS THAN THREE MESH POINTS INSIDE NUCLEUS.')
      SKIP=1.0D0
    3 IF (KPOT.EQ.1) GO TO 2
      IF (KPOT.EQ.4) GO TO 2
      IF (KPOT.NE.2) GO TO 6
      GO TO 5
    2 READ (ilog,20) XALPHA,XNNN,XLATTR
   20 FORMAT(3F10.5)
      IF (XALPHA.LT.0.0D0) XALPHA=2.0D0/3.0D0
      GO TO 4
    5 XALPHA=2.0D0/3.0D0
    6 CONTINUE
   15 FORMAT (' XALPHA NOT SPECIFIED. XALPHA MADE 1.0')
      IF (XALPHA.EQ.0.0) WRITE (6,15)
      IF (XALPHA.EQ.0.0) XALPHA=1.0D0
    4 SUM=0.0D0
   16 FORMAT (3F5.1,F14.7,F18.10,20X,F4.0)
      IF(IREAD.EQ.2) GO TO 200
      IF(IREAD.EQ.6) GO TO 200
      IF(IREAD.EQ.3) GO TO 200
      bcond=0.d0
      DO 7 I=1,J
      READ (ilog,16) XN(I),XL(I),XJ(I),XE(I),XZ(I)
      WRITE (6,16) XN(I),XL(I),XJ(I),XE(I),XZ(I)
      IF ((XN(I)-XL(I)-1.0D0).LT.0.0) GO TO 8
      IF ((DABS(XL(I)-XJ(I))-0.5D0).NE.0.0) GO TO 8
      IF ((2.0D0*XJ(I)+1.0D0).LT.XZ(I)) GO TO 8
      IF (IREAD.EQ.4) XN(I)=-DABS(XN(I))
    7 CONTINUE
      IF (IREAD.EQ.4) READ (5,2000) BCOND
      GO TO 210
  200 CONTINUE
      READ (ilog,2000) (XZ(I),I=1,J)
 2000 FORMAT(10F8.4)
      IF(IREAD.EQ.3) READ(ilog,2000) BCOND
  210 CONTINUE
      DO 205 I=1,J
  205 SUM=SUM+XZ(I)
      BCON=BCOND
      IF(DABS(SUM+XION-ZN).GT.1.0D-7) GO TO 9
      IF(RN.LT.EPS) GOTO 10
      WRITE(6,6000)
 6000 FORMAT('ERROR ON RN OR EPS ')
      SKIP=1.0
      GO TO 10
    8 WRITE (6,18) AAME
   17 FORMAT(' ERROR ON EIGENVALUE CARD',I5)
      WRITE (6,17) I
      SKIP=1.0D0
      GO TO 10
    9 WRITE (6,18) AAME
   19 FORMAT (' TOTAL NUMBER OF ELECTRONS IN ERROR')
      WRITE (6,19)
      SKIP=1.0D0
   10 IF (SKIP.EQ.1.0D0) GO TO 1
      RETURN
  100 STOP
      END
      FUNCTION ZEFF(E,FN,FK)
      IMPLICIT REAL*8 (A-H,O-Z),INTEGER*4 (I-N)
      DATA O/137.03602/
      G=1.+E/O**2
      G=1./G**2-1.
      A=1.+G
      B=FN-FK
      C=DSQRT(G)
      D=FN*(FN-2.D0*FK)
      F=DSQRT(B*B-A*D)
      Z=(C*B+C*F)/A
      ZEFF=Z*O
      RETURN
      END
C SEGMENT CMMN1+NP
      SUBROUTINE DIFF1(KI,FN1)
C  INTEGRATES OUTWARD TO CLASSICAL TURNING POINT. COUNTS N0DES IN THE
C  MAJOR COMPONENT OF THE RADIAL WAVE FUNCTION.
C  CALL ED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
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
      DO 1 K=11,M
      KM=N-K
      R=RAD(KM)
      IF (E*R+Z-VR(KM).GT.0.0) GO TO 2
    1 CONTINUE
    2 KI=KM+1
      DO 3 K=5,KI
      R=RAD(K)
      RP21=(E*R-VR(K)-VRX(K)+Z)/CS
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
C SEGMENT CMMN1+NP
      SUBROUTINE DIFFER(DELL)
C  CONTROL SUBROUTINE FOR THE INTEGRATION OF THE DIFFERENTIAL EQUATIONS.
C  FINDS EIGENVALUES. NORMALIZES ORBITAL FUNCTIONS.
C  CALLED BY HEX AND BIND.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   RINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      EMIN=-(ZN/FN)**2
      EMAX=-((XION+XLATTR)/FN)**2/2.5
      DELK=0.01
      IF(BCOND.NE.0.) EMAX=DMAX1(1.D0,10./RN**2)
      IF (E.GT.EMAX.OR.E.LT.EMIN) E=0.5D0*(EMAX+EMIN)
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
      IF (FN-FN1) 2,4,13
C  TOO MANY NODES.
    2 IF (E.LT.EMAX) EMAX=E
      E=0.5D0*(E+EMIN)
      DL1=DABS(EMAX-EMIN)/(DABS(E)+1.D0)
      IF (DL1.GE.DELK) GO TO 1
      WRITE (6,15) FN1,FN,FL,FJ,EMIN,E,EMAX
   15 FORMAT (16H0TOO MANY NODES.,4F4.1,3D16.8)
      EMIN=1.2D0*EMIN
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
      appp = DSQRT(-2.0D0*E)*(R-RJ)
      w = 0.d0
      if(appp.lt.150.d0) w = dexp(-appp)
c      W=DEXP(-DSQRT(-2.0D0*E)*(R-RJ))
      A(K)=A(KJ)*W
    7 B(K)=B(KJ)*W
    8 DO 9 K=1,N
      R=RAD(K)
    9 Y(K)=(A(K)*A(K)+B(K)*B(K))*R
      W=ADLINT(Y,N,H)
      W=W+Y(1)/(G+G+1.0D0)
      DE=CS*A(KI)*B(KI)*(1.0D0-RB/RA)/W
      DL1=DABS(EMAX-EMIN)/(DABS(E)+1.D0)
      DL=DABS(DE/E)
      IF ((DL.GT.DELK) .AND. (DL1.LT.DELK)) GO TO 12
      IF (DL1.LT.DELL) GO TO 12
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
   17 FORMAT (16H0TOO FEW  NODES.,4F4.1,3D16.8)
      IF (BCOND.NE.0.) STOP
      BCOND=1.
      WRITE (6,6001) BCOND
 6001 FORMAT (' ***** RESET B.C.=1. *****',F10.5)
      EMAX=DMAX1(1.0D0,10.D0/RN**2)
      E=0.5D0*(E+EMAX)
      GO TO 1
      END
C SEGMENT CMMN1+NP
      SUBROUTINE START2(KI,KJ)
C  STARTS INWARD INTEGRATION. OUTER BOUNDARY CONDITION IS SET HERE.
C  CALLED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,EEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      H3=-H/3.0D0
      S=2.0D0*(FJ-FL)
      CS=C*S
      G=DSQRT((FJ+0.5D0)**2-(Z/CS)**2)
      Q22=-G-S*(FJ+0.5D0)
      Q11=-G+S*(FJ+0.5D0)
      DO 1 K=KI,N
      R=RAD(K)
      IF ((E*R-VR(K)+Z)*R+EPS) 2,1,1
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
    3 IF (BCOND.NE.1.) GO TO 23
      B(KJ)=DMOD(FL,2.D0)
      A(KJ)=1.0D0-B(KJ)
      GO TO 40
   23 IF(BCOND.NE.2.)GO TO 33
          A(KJ)=DMOD(FL,2.D0)
      B(KJ)=1.-A(KJ)
      GO TO 40
   33 IF(BCOND.NE.3.) GO TO 35
      A(KJ)=1.
      B(KJ)=1.
      GO TO 40
   35 A(KJ)=1.
      B(KJ)=-1.
   40 CONTINUE
    4 DO 5 L=1,4
      K=KJ-L
      A(K)=A(KJ)
    5 B(K)=B(KJ)
      DO 7 I=1,4
      K=KJ+1
      DO 6 L=1,5
      K=K-1
      R=RAD(K)
      RP21=(E*R-VR(K)-VRX(K)+Z)/CS
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
C SEGMENT CMMN1+NP
      SUBROUTINE START1
C  CALLED BY DIFFER.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4 BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
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
C SEGMENT CMMN1+NP
      SUBROUTINE OUTPUT
C  OUTPUT PRINTS DATA FROM THE CALCULATION. IF NPRINT = 1 THE ORBITAL
C  FUNCTIONS ARE READ FROM TAPE 4 AND PRINTED.
C  CALLED BY HEX.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      character*4 alj
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DIMENSION ALJ(50)
      DATA             ALJ/6H 1S1/2,6H 2S1/2,6H 2P1/2,6H 2P3/2,
     1   6H 3S1/2,6H 3P1/2,6H 3P3/2,6H 3D3/2,6H 3D5/2,6H 4S1/2,
     2   6H 4P1/2,6H 4P3/2,6H 4D3/2,6H 4D5/2,6H 4F5/2,6H 4F7/2,
     3   6H 5S1/2,6H 5P1/2,6H 5P3/2,6H 5D3/2,6H 5D5/2,6H 5F5/2,
     4   6H 5F7/2,6H 5G7/2,6H 5G9/2,6H 6S1/2,6H 6P1/2,6H 6P3/2,
     5   6H 6D3/2,6H 6D5/2,6H 6F5/2,6H 6F7/2,6H 6G7/2,6H 6G9/2,
     6   6H 6H9/2,6H6H11/2,6H 7S1/2,6H 7P1/2,6H 7P7/2,6H 7D3/2,
     7   6H 7D5/2,6H 7F5/2,6H 7F7/2,6H 7G7/2,6H 7G9/2,6H 7H9/2,
     8   6H7H11/2,6H7I11/2,6H7I17/2,6H 8S1/2/
C     DATA P4/0.0795774715459477D0/
C
      

      REWIND 4
      WRITE (6,11) AAME
   11 FORMAT (1H1,30X,10A8,//)
      WRITE (6,13) XALPHA,XLATTR,XNNN,FEX,RNUC
   13 FORMAT(' EXCHANGE COEFFICIENT=',F10.7,5X,' TAIL CORRECTION=',
     1F3.0,5X,' XNNN=',F7.5,2X,' FEX=',F4.1,5X,' NUCLEAR RADIUS=',D10.3)
      IF (FEX.EQ.1.0D0) WRITE (6,21)
   21 FORMAT (' FREE ELECTRON EXCHANGE POTENTIAL')
      IF (KPOT.EQ.1) GO TO 1
      IF (KPOT.EQ.2) GO TO 2
      IF (KPOT.EQ.4) GO TO 31
      IF (FEX.LE.1.0D0) WRITE (6,12)
   12 FORMAT (1H0,23X,7HORBITAL,4X,1HN,4X,1HL,4X,1HJ,5X,9HELECTRONS,5X,
     110HEIGENVALUE,5X,2HPV,/)
      IF (FEX.EQ.2.0D0) WRITE (6,24)
   24 FORMAT (1H0,14X,7HORBITAL,4X,1HN,4X,1HL,4X,1HJ,5X,9HELECTRONS,5X,
     110HEIGENVALUE,5X,2HPV,10X,18HH-F BINDING ENERGY,/)
      IF (FEX.EQ.3.0D0) WRITE (6,25)
   25 FORMAT (1H0,14X,7HORBITAL,4X,1HN,4X,1HL,4X,1HJ,5X,9HELECTRONS,5X,
     110HEIGENVALUE,6X,11HE(N)-E(N-1),6X,11HE(N)-E(N+1),/)
      IFEX=FEX+1
      GO TO 3
    1 WRITE (6,26) R1
   26 FORMAT(1H0,40X,' NEW POTENTIAL USED. R1=',F8.5)
      IFEX=4
      WRITE (6,25)
      GO TO 3
   31 WRITE (6,32)
   32 FORMAT(1H0,40X,'ROSEN-LINDGREN POTENTIAL USED.')
      WRITE (6,25)
      GO TO 3
    2 WRITE (6,27)
   27 FORMAT (1H0,45X,' KOHN-SHAM POTENTIAL USED. ')
      IFEX=4
      WRITE (6,25)
      PVT=0.
    3 DO 7 I=1,J
      READ (4) (A(K),K=1,N)
      READ (4) (B(K),K=1,N)
      XN(I)=DABS(XN(I))
      FN=XN(I)
      FK=XJ(I)+0.5
      NI=XN(I)*(XN(I)-2.0D0)+XL(I)+XJ(I)+1.501D0
      IF(BCOND.NE.0.) GO TO 500
      PV(I)=ZEFF(XE(I),FN,FK)
      PV(I)=ZN-PV(I)
      GO TO 600
  500 E=XE(I)
      S=2.*(XL(I)-XJ(I))
      C=S*137.036
      RHNI=A(N)**2+B(N)**2
      PV(I)=(RHNI*(E*RN+Z-VR(N)-VRX(N))+2.*(C*B(N))**2*RN-2.*FK*C*A(N)*
     1B(N))/RN**3
      RHNI=RHNI/RN**3
      PV(I)=PV(I)-0.158837*XALPHA*RHNI*(RN**XNNN*RH(N))**.333333333333D0
      PV(I)=294.0*PV(I)
      PVT=PVT+PV(I)*XZ(I)
  600 CONTINUE
      KM=1
      WFM=0.
      DO 900 K=1,N
      WF=A(K)*A(K)
      IF(WF.LT.WFM) GO TO 900
      WFM=WF
      KM=K
  900 CONTINUE
      BINDEN(I)=RAD(KM)
    5 WRITE (6,28) ALJ(NI),XN(I),XL(I),XJ(I),XZ(I),XE(I),PV(I),BINDEN(I)
   28 FORMAT (15X,A6,2X,3F5.1,F10.1,1PD19.7,D15.7,D19.7)
    7 CONTINUE
      IF(BCOND.NE.0.) WRITE(6,6001) PVT,BCOND
 6001 FORMAT(' ***** PRESSURE(MBAR)=',1PE15.7,'  FOR B.C.(',0PF4.1,')')
      EKIN=SUMMA-EV3
      IF (RNUC.EQ.0.0) EKIN=EKIN-EV4
      EE2=-EE2
      WRITE (6,15) ZN,Q,SUMMA,ET3,EKIN,EV4,EV,EE2
   15 FORMAT (1H0,39X,15HNUCLEAR CHARGE=,F12.6/28X,27HINTEGRAL OF CHARGE
     1 DENSITY=,F12.6/25X,30HSUM OF THE ENERGY EIGENVALUES=,F12.4/42X,1
     23HTOTAL ENERGY=,F12.4/40X,15HKINETIC ENERGY=,F12.4/34X,22HE-N POT
     3ENTIAL ENERGY=,F12.4/34X,21HE-E POTENTIAL ENERGY=,F12.4/39X,16HEXC
     4HANGE ENERGY=,F12.4)
      EE2=-EE2
      WRITE (6,11) AAME
      WRITE (6,23) RN,N,H
   23 FORMAT(12X,41H TABLE OF RADII. RN=RMAX*EXP(N-NMAX)*H).,' RMAX=',
     1F8.4,' NMAX=',I5,' H=',F10.8/)
      WRITE (6,17) (RAD(K),K=1,N)
c     WRITE (50,17) (RAD(K),K=1,N)
   17 FORMAT (1H 1P8D14.7)
      DO 9 K=1,N
    9 VRX(K)=(Z-VR(K))/ZN
      WRITE (6,11) AAME
      WRITE (6,16)
c     WRITE (50,16)
   16 FORMAT (53X,'RH0(R)'/)
      WRITE (6,17) (RH(K),K=1,N)
      WRITE (6,11) AAME
      WRITE (6,22)
c     WRITE (50,17) (RH(K),K=1,N)
c     WRITE (50,11) AAME
c     WRITE (50,22)
   22 FORMAT (47X,'POTENTIAL SCREENING FACTORS'/)
      WRITE (6,17) (VRX(K),K=1,N)
c     WRITE (50,17) (VRX(K),K=1,N)
      llp=2
	write(50,*)'here is vscf',zn
      write(50,*)'   rad           vscf   '
	do i=1,n
      WRITE(50,'(2e14.5)') rad(i),-VRX(I)*ZN/rad(i)
 	enddo

      IF (NPRINT.NE.1) RETURN
      REWIND 4
      REWIND 3
      DO 10 K=1,J
      READ (4) (A(I),I=1,N)
      READ (4) (B(I),I=1,N)
      NI=XN(K)*(XN(K)-2.0D0)+XL(K)+XJ(K)+1.501D0
      WRITE (6,20) AAME
   20 FORMAT (1H1,10A8/)
      WRITE (6,18) ALJ(NI)
c      WRITE (3,18) ALJ(NI)
   18 FORMAT(1H ,20X,'R TIMES MAJOR COMPONENT OF THE ',A6,' ORBITAL')
c      EEEE=XE(K)*0.0272D0
c      WRITE(3,41) XN(K),XL(K),XJ(K),EEEE,XN(K),BCOND
   41 FORMAT(1X,3F5.1,F14.7,F7.4,F4.0)
      DO 44 I=1,N
      IF(ABS(A(I)).LT.1.D-50) A(I)=0.D0
      IF(ABS(B(I)).LT.1.D-50) B(I)=0.D0
 44   CONTINUE
      HP=1.D0/H
c     WRITE(10) N,HP,RN,XE(K),XN(K),XL(K),XJ(K),PHI,XION,Z
c     WRITE(10) (A(I),I=1,N),(B(I),I=1,N)
c      WRITE(3,42) A
c      WRITE(3,42) B
 42   FORMAT(5E14.7)
      WRITE (6,17) (A(I),I=1,N)
      WRITE (6,20) AAME
      WRITE (6,19) ALJ(NI)
   19 FORMAT(1H ,20X,'R TIMES MINOR COMPONENT OF THE ',A6,' ORBITAL')
      WRITE (6,17) (B(I),I=1,N)
   10 CONTINUE
c      CLOSE(10)
      RETURN
      END
C SEGMENT CMMN1+NP
      SUBROUTINE POT(RHO,UR)
C  POT COMPUTES THE SELF-CONSISTENT FIELD P0TENTIAL FUNCTION FR0N THE
C  IS COMPUTED ITIS COMPARED WITH THE OLD ONE AND CONVRG IS SET EQUAL
CC  CHARGE DENSITY AND SEVERAL T0TAL ENERGIES. AS THE POTENTIAL FUNCTION
C  TO THE MAXIMUM DIFFERENCE BETWEEN THE TWO. POT IS ALSO USED IN THE
C  COMPUTATION OF ELECTRON BINDING ENERGIES BY BIND.
C  CALLED BY HEX AND PIND.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      COMMON /NP/ RNK,ANK,BNK,NK
      DIMENSION RHO(NMX), UR(NMX),C(NMX), DC(5)
C
      DATA PI/3.141592653589793D0/
      CONST=1.5D0/PI*QBRT4(0.75D0*PI)
      CONSTT=CONST*XALPHA
      S2=2.0D0
      Z=ZN
      IF (RNUC .EQ.0.0) ZN=0.0
      NK=0
      RNK=RN
      ANK=Z-XLATTR-XION-XKOOP
      RU=ANK
      IF (RU.LT.0.0) RU=0.0
      S3=S2+1.0D0
      ENUC=0.0
      IF (RNUC.GT.0.0) RU=RU-ZN
      H3=H/3.0D0
      R=RAD(1)
      A(1)=RHO(1)*R/S3
      B(1)=RHO(1)/S2
      DB(1)=H3*RHO(1)
      DA(1)=DB(1)*R
      S5=S2+3.0D0
      C(1)=RHO(1)*R**3/S5
      DC(1)=DA(1)*R**2
      R=RAD(2)
      A(2)=RHO(2)*R/S3
      B(2)=RHO(2)/S2
      DB(2)=H3*RHO(2)
      DA(2)=DB(2)*R
      NC=0
      C(2)=RHO(2)*R**3/S5
      DC(2)=DA(2)*R**2
      print *,'rho1',rho(1),rho(2), R, H, ANK, RNK, BNK
      print *, 'a1b1b2a2b2c2', A(1), A(2), B(1),B(2),B(1),C(2)
      DO 1 K=3,N
      R=RAD(K)
      DB(3)=H3*RHO(K)
      DA(3)=DB(3)*R
      DC(3)=DA(3)*R**2
      A(K)=A(K-2)+DA(3)+4.0D0*DA(2)+DA(1)
      B(K)=B(K-2)+DB(3)+4.0D0*DB(2)+DB(1)
      C(K)=C(K-2)+DC(3)+4.0D0*DC(2)+DC(1)
      IF (A(K).LT.ANK) NK=K
      IF (R.LE.RNUC) NC=K
      DA(1)=DA(2)
      DB(1)=DB(2)
      DA(2)=DA(3)
      DC(1)=DC(2)
      DC(2)=DC(3)
    1 DB(2)=DB(3)
      print *, (B(K),k=1,5)
      NC=NC+1
      NK=MIN0(NK+1,N)
      QA=A(N)
      ZB=-Z*B(N)
      EZ=0.0
      IF (RNUC.EQ.0.0) GO TO 2
      ANC=CUBINT(RNUC,NC,RAD,A,N)
      BNC=CUBINT(RNUC,NC,RAD,B,N)
      CNC=CUBINT(RNUC,NC,RAD,C,N)
      EZ=ZB-Z*(1.5D0*ANC/RNUC-BNC-0.5D0*CNC/RNUC**3)
      ZB=EZ
    2 IF (KPOT.NE.1) GO TO 4
      RNK=0.0
      BNK=0.0
      DX=0.0
      IF (NK.LT.3) GO TO 5
      RNK=CUBINT(ANK,NK,A,RAD,N)
      BNK=CUBINT(ANK,NK,A,B,N)
      DX=DLOG(RAD(NK)/RNK)
      GO TO 5
    4 NK=N
      RNK=RN
      BNK=B(N)
      DX=0.0
    5 DO 6 K=1,NK
      R=RAD(K)
      IF(RHO(K).LT.0)
     *WRITE(6,6348) K,RHO(K),R
 6348 FORMAT(5X,I4,1P1E14.6,8X,1P1D12.6)
    6 Y(K)=-CONSTT*QBRT4(RHO(K)*R**XNNN)
      VNK=0.0
      V0=BNK-VNK+Y(1)/RAD(1)-1.5D0*ZN/(RNUC+1.D-40)
      print *, 'v0,bnk,rad1',v0,bnk,rad(1)      
      BN=B(N)
      R0=0.0
      EPSLON=0.0
      DO 17 K=1,N
      R=RAD(K)
      IF (R.GE.RNUC) GO TO 8
      X=R/RNUC
      RVN=-Z*X*(1.5D0-0.5D0*X**2)
      GO TO 9
    8 RVN=-ZN
    9 IF (KPOT.NE.1) GO TO 11
      IF (K.LT.NK) GO TO 12
      GO TO 13
   11 IF (R0.GT.0.0) GO TO 13
   12 RV=A(K)+R*(BNK-B(K)-VNK)+Y(K)+RVN
      IF( RV.GE.RU) GO TO 114
      GO TO 14
  114 R0=1.
   13 RV=RU
   14 CONTINUE
      IF (XKOOP.NE.0.0) GO TO 16
      RX=RV
      IF(RX.LE.0.) RX=1.
      ERROR=DABS(RV-UR(K))/RX
      IF (ERROR.LE.EPSLON) GO TO 16
      EPSLON=ERROR
      NERROR=K
   16  UR(K)=RV
      RHNLJ(K)=RHO(K)*UR(K)
      Y(K)=0.5D0*CONST*QBRT4(RHO(K)*R**XNNN)*RHO(K)
      A(K)=A(K)*RHO(K)
   17 CONTINUE
      ZN=Z
      IF (RNUC.GT.0.0) Z=0.0
      EY=0.0
      ED=0.0
      ER=0.0
      IF (NK.LT.10) GO TO 19
      EY=ADLINT(Y,NK,H)+Y(1)/3.0D0
      ED=ADLINT(A,NK,H)+A(1)/5.0D0
      ER=ADLINT(RHNLJ,NK,H)+RHNLJ(1)/3.0D0
      IF (NK.GE.N) GO TO 19
      EY=EY- ENDCOR(Y,DX,NK,H)
      ED=ED- ENDCOR(A,DX,NK,H)
      ER=ER- ENDCOR(RHNLJ,DX,NK,H)-(ZN-Z)*(BN-BNK)
   19 IF (XKOOP.NE.0.0) RETURN
      EE2=EY
      EV=ED
      EV3=ER
      Q=QA
      EV4=ZB
      CONVRG=DABS(EPSLON)
      IF(CONVRG.EQ.1.) CONVRG=1.01
      WRITE (6,20) EPSLON,NERROR
   20 FORMAT(' MAXIMUM ERROR IN V(R) IS',D15.8,'AT THE',I5,'TH POINT')
      RETURN
      END
C SEGMENT CMMN1+NP
      SUBROUTINE PUNCHC(K)
C  PUNCHC PUNCHES AN OUTPUT DECK OF THE SAME FORM AS THE INPUT DECK
C  CALLED BY HEX.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DIMENSION CDEN(nmx)

      ilog=9
      WRITE (ilog,8) AAME
    8 FORMAT (10A8)
      G=1.0D0/H
      N0=1
      WRITE (ilog,9) N,J,RN,G,ZN,XION,PHI,EPS,DEL,DELRVR,AAME(1),N0
    9 FORMAT (2I5,F6.3,F4.0,2F5.0,4F10.7,2X,A5,I3)
      N0=2
      WRITE (ilog,10) NC1,AAME(1),N0
   10 FORMAT (I5,67X,A5,I3)
      N0=3
      WRITE (ilog,11) KPOT,XALPHA,XLATTR,XNNN,FEX,RNUC,AAME(1),N0
   11 FORMAT (I10,F10.7,9X,F5.1,F10.5,F5.1,D10.3,13X,A5,I3)
      DO 2 I=1,J
      BCOND=0.0
      IF (XN(I).LT.0.0) BCOND=1.0D0
      XN(I)=DABS(XN(I))
      N0=N0+1
      WRITE(ilog,12) XN(I),XL(I),XJ(I),XE(I),XZ(I),PV(I),BINDEN(I),
     1   AAME(1),N0
    2 CONTINUE
   12 FORMAT (3F5.1,F14.7,F7.4,5X,F9.4,E12.4,10X,A5,I3)
      FN=0.0
      FL=0.0
      FJ=0.0
      E=0.00
      N5=(N/5)*5+5
      DO 20 I=1,N5
   20 CDEN(I)=0.
      DO 25 I=1,N
   25 CDEN(I)=RH(I)
      DO 3 I=1,N,5
      L=I+4
      N0=N0+1
    3 WRITE (ilog,14) (CDEN(I1),I1=I,L),AAME(1),N0
   14 FORMAT (1P5D14.7,2X,A5,I3)
      DO 30 I=1,N
   30 CDEN(I)=VRX(I)
      DO 33 I=1,N,5
      L=I+4
      N0=N0+1
   33 WRITE (ilog,14) (CDEN(I1),I1=I,L),AAME(1),N0

c   jinrui for orbital energy reading 2015.05.28       
      open(43,file='atomlib.out')
      write(43,*)int(zn)
      write(43,*)int(xion),J,0.d0
      do i = 1, j
c      write(*,*) i , j , xn(i),xl(i),xj(i)
      xkap = ( xl(i) - xj(i)) * ( 2.0 * xj(i) + 1.0) 
      kap = int ( xkap )
      nxn = int(xn(i))
      WRITE(43,19) nxn ,kap,int(XZ(i)),XE(I)
      enddo
19    format(1X,3I5,E16.8)

c      
      RETURN
      END
C SEGMENT CMMN1+NP
      SUBROUTINE START0
C  COMPUTES AN INITIAL RADIAL CHARGE DISTRIBUTION
C     CALLED BY HEX.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XALTTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT

      ALPHA=0.3058D0*(ZN-XION)**(1.0D0/3.0D0)
      BETA=5.8632D0*DMAX1(0.0D0,ZN-XION-2.0D0)*ALPHA
      GAMMA=DSQRT(4.0D0+3.2D0*XION)
      G2=GAMMA**3
      IF ((ZN-XION).LT.2.0D0) G2=0.5D0*(ZN-XION)*G2
      DO 1 K=1,N
      R=RAD(K)
      X=ALPHA*R
    1 RH(K)=BETA*DSQRT(X)*DEXP(-3.0D0*X)+G2*R**2*DEXP(-GAMMA*R)
    2 FORMAT(8D16.8)
      print "(3E32.22)", RH(1), RH(2),RH(3)
      IF (NDEBUG.NE.0) WRITE (6,2) (RH(I),I=1,N)
      RETURN
      END
      SUBROUTINE FEEDBK(N1,NUM,RHP1,DEP1,RHP2,DEP2,STE,STORE,ERROR,AVF)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      DIMENSION STORE(NMX),STE(NMX),RHP1(NMX),RHP2(NMX),AVF(NMX)
     1                 ,ERROR(5),DEP1(NMX),DEP2(NMX)
      DO 1 I=1,N1
      STE(I)=0.0
      STORE(I)=0.0
    1 AVF(I)=0.0
      ASUM=0.0
      ALPH=0.0
      NUM=NUM+1
      IF(NUM-2) 2,2,4
    2 DO 3 I=1,N1
    3 STE(I)=0.5
      GO TO 15
    4 IF(ERROR(5).GE.ERROR(4))GO TO 2
      DO 12 I=1,N1
      AVF(I)=DEP2(I)-DEP1(I)-RHP2(I)+RHP1(I)
      IF(AVF(I).EQ.0.0.AND.RHP2(I).EQ.0.0) GO TO 5
      IF(AVF(I).NE.0.0.AND.RHP2(I).EQ.0.0)  GO TO 6
      IF(DABS(AVF(I)/RHP2(I))-0.001)5,5,6
    5 CONTINUE
      ALPH=0.5
      GO TO 11
    6 ALPH=(DEP2(I)-DEP1(I))/AVF(I)
      IF(ALPH)7,10,8
    7 ALPH=0.0
      GO TO 11
    8 IF(0.5-ALPH) 9,10,10
    9 ALPH=0.5
   10 CONTINUE
   11 STORE(I)=ALPH
   12 CONTINUE
      STORE(1)=0.5
      STE(2)=STORE(2)
      ASUM=STORE(1)+STORE(2)+STORE(3)+STORE(4)+STORE(5)
      K=N1-3
      DO 13 I=3,K
      STE(I)=ASUM*0.2
   13 ASUM=ASUM-STORE(I-2)+STORE(I+3)
      STE(K+1)=STE(K)
      STE(K+2)=STE(K)
      STE(K+3)=STE(K)
   15 CONTINUE
      DO 16 I=1,N1
      RHP1(I)=RHP2(I)
      DEP1(I)=DEP2(I)
   16 CONTINUE
      DO 14 I=1,N1
      STE(I)=DABS(STE(I))
      RHP2(I)=DEP2(I)+STE(I)*(RHP2(I)-DEP2(I))
      IF(RHP2(I).LT.0)
     *WRITE(6,398) RHP2(I),I,STE(I)
   14 CONTINUE
  398 FORMAT('IN THE FEEDBK  RHP2(I) IS'/D14.6,5X,I3,2X,'STE'/D10.7)
      RETURN
      END
      DOUBLE PRECISION FUNCTION ENDCOR(Y,DX,N,H)
C  QUADRATIC END CORRECTION FOR INTEGRATION.
C  CALLED BY POT.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      DIMENSION Y(NMX)
      D=DX/H
      ENDCOR =DX*(Y(N)-0.25D0*(Y(N+1)-Y(N-1))*D+
     1      0.166666666666667D0*(Y(N+1)-2.0D0*Y(N)+Y(N-1))*D*D)
      RETURN
      END
      DOUBLE PRECISION FUNCTION CUBINT(XJ,N,X,Y,NX)
C  CUBIC LAGRANGE INTERPOLATI0N TO FIND Y(XI). N IS A LOCATION IN THE
C  TABLE OF X NEAR XI. UXUALLY XI IS BETWEEN X(N-1) AND X(N)
C  CALLED BY POT.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
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
      WRITE (6,2) N,X(N-2),X(N-1),X(N),X(N+1)
    2 FORMAT(' TROUBLE IN CUBINT. SET TO Y(N).',I5,4D16.8)
      RETURN
      END
      DOUBLE PRECISION FUNCTION ADLINT(Y,N,H)
C  EULER-MACLAURIN INTEGRATION FORMULA. N MUST EXCEED 10.
C  CALLED BY POT,DIFFER,AND BIND.
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
      DIMENSION Y(NMX)
      DATA CA,CB,CC,CD,CE/0.329861111D0,1.320833333D0,0.766666667D0,
     11.101388889D0,0.981250000D0/
      SUM=CA*(Y(1)+Y(N))+CB*(Y(2)+Y(N-1))+CC*(Y(3)+Y(N-2))+CD*(Y(4)+Y(N-
     13))+CE*(Y(5)+Y(N-4))
      M=N-5
      DO 1 K=6,M
    1 SUM=SUM+Y(K)
      ADLINT=SUM*H
      RETURN
      END
      DOUBLE PRECISION FUNCTION QBRT4(A)
      IMPLICIT REAL*8(A-H,O-Z)
      QBRT4=A**(0.3333333333333333D0)
      RETURN
      END
C SEGMENT CMMN1+NP
      SUBROUTINE DIFF2(KI,KJ)
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (NMX=2000,jm=50)
C   INTEGRATES INWARD TO THE CLASSICAL TURNING POINT
C   CALLDE BY DIFFER
      COMMON /CMMN1/ AAME(10),RN,H,Z,ZN,XION,PHI,EPS,DEL,DELRVR,XALPHA,
     1   XLATTR,FEX,RNUC,XN(jm),XL(jm),XJ(jm),XE(jm),XZ(jm),BCOND,XNNN,
     2   FN,FL,E,FJ,VR(NMX),Q,EV,ED,EV3,ER,EV4,EZ,EE2,EY,V0,ET3,SUMMA,
     3   PVT,CONVRG,PV(jm),A(NMX),B(NMX),RHNLJ(NMX),DENS(NMX),Y(NMX),
     4   BINDEN(jm),DA(5),DB(5),A0(5),B0(5),VRX(NMX),XKOOP,R1,RAD(NMX),
     5   RH(NMX),KPOT,N,J,NC1,NDEBUG,NPRINT
      DATA C/137.03602D0/
      do 1000 i=kj+1,n
      a(i) = 0.d0
      b(i) = 0.d0
 1000 continue
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
      RP21=(E*R-VR(K)-VRX(K)+Z)/CS
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
