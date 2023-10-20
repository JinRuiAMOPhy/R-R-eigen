      PROGRAM DRIVER
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'darc.inc'
      CHARACTER AELEM*2,ASTAT*30(MXCH)
      DIMENSION ITJPO(MXCH),ISPAR(MXCH),KC(MXCH),ETAR(MXCH),
     +QDT(MXCH),NSRT(MXCH),XJ(mxch),lcf(mxch)
      COMMON /XQDT/QDT,NKC,NSRT
      COMMON /XORBS/XEORB(MXCH),IQDT(MXCH),IEORB
      common /elemt/ion,AELEM
      character (len = 10) :: how_to_calc_orb
      open(1,file='analyz.in')
      print *,"element "
      read(1,*)AELEM
      print *, 'charge state'
      read(1,*)ion
      write(*,*)'input the     total 2J+1      for N+1 system'
      READ(1,*)J2P1
      write(*,*)'input the     PARITY    for N+1 system   odd -1
     +   even   1'
      READ(1,*)NPTY
      write(*,*)'input number of    CONTINUUM ORBITAL KAPPA'
      READ(1,*)NKC
      READ(1,*)how_to_calc_orb
      DO IKC=1,NKC
        II=MOD(IKC,2)
        IF(II.EQ.1)THEN
          KC(IKC)=-(IKC/2+1)
        ELSEIF(II.EQ.0)THEN
          KC(IKC)=(IKC/2)
        ENDIF
      ENDDO
      call findK(kc(1),kcmin,nkc)
      call findK(kc(nkc),kcmax,nkc)
      !write(*,*)(KC(I),I=1,NKC)
      !write(*,*)'MAX KC   MIN KC', KCMAX,KCMIN
      ASTAT(:) = '' 
      call nist(NAST,ASTAT,ETAR,xj,ISPAR,lcf)
      do iast=1,nast
        itjpo(iast)=int(xj(iast)*2.0)
      enddo
c      !write(*,*)(ASTAT(iast),iast=1,nast)
c      write(*,*)'reading nist done'
      close(1)
      IF((index(how_to_calc_orb, 'read') /= 0) .or. 
     & (index(how_to_calc_orb, 'QDT') /= 0) )THEN
        OPEN(2,FILE='QDT.IN')
        DO IKC=1,NKC
        READ(2,*)IQDT(IKC),NSRT(IKC),QDT(IKC)
        ENDDO
        CLOSE(2)
      ELSE if (index(how_to_calc_orb, 'calc') /= 0) then
        CALL RSCFORB(NKC)
      ENDIF
c      write(*,*)'reading QDT done'

      xion=dble(ion)
      !write(*,*)'kcmax,kcmin',kcmax,kcmin,nast
      CALL CONT(ASTAT,ETAR,ISPAR,ITJPO,J2P1,KCMAX,KCMIN,NAST,NPTY,
     +XION,lcf,how_to_calc_orb)

      END

      SUBROUTINE A2J(BB,JJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER BB*4
      BB=ADJUSTR(BB)
      IONE=0
      IDEC=0
      JJ=0
      IF(BB(4:4).EQ.'2'.AND.BB(3:3).EQ.'/')THEN

        IF(BB(1:1).NE.' ')IDEC=ICHAR(BB(1:1))-ICHAR('0')
        IONE=ICHAR(BB(2:2))-ICHAR('0')
        JJ=IDEC*10+IONE
      ELSE
        IF(BB(3:3).NE.' ')IDEC=ICHAR(BB(3:3))-ICHAR('0')
        IONE=ICHAR(BB(4:4))-ICHAR('0')
        JJ=IDEC*10+IONE
        JJ=2*JJ
      ENDIF
      END
      SUBROUTINE A2PTY(BB,JJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER BB*4
      BB=ADJUSTR(BB)
      IF(BB(4:4).EQ.'-')THEN
        JJ=-1
      ELSEIF(BB(4:4).EQ.'+')THEN
        JJ=1
      ENDIF
      END
      
      SUBROUTINE CONT(ASTAT,ETAR,ISPAR,ITJPO,J2P1,KCMAX,KCMIN,NAST,
     +NPTY,QTAR,lcf,how_to_calc_orb)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      INCLUDE 'darc.inc'
C
C  Statement functions
C
      LOGICAL ITRG
C
C  Parameter variables
C
      DOUBLE PRECISION EPS10
      PARAMETER (EPS10=1.D-10)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION QTAR
      PARAMETER (ZERO=0.D0)
C
C  Argument variables
C
C
      CHARACTER*2 NHC(MXCH)
     
      INTEGER ISPAR(MXCH)
      INTEGER ITJPO(MXCH)
      INTEGER ITAR
      INTEGER Iwrite
      INTEGER J2P1
      INTEGER K2P(MXCH)
      INTEGER KCMAX
      INTEGER KCMIN
      INTEGER NAKK(MXCH)
      INTEGER NAST
      INTEGER NCHAN
      INTEGER NPTY
      INTEGER NTARG(MXCH)
C
C  Local variables
C
      character ASTAT*30(MXCH)
      integer lcf(mxch),INVINDX(mxch)
c JINRUI 2016.11.24
      DOUBLE PRECISION XEORB
      INTEGER IEORB,IQDT
      CHARACTER*2 LAB(22)
      DOUBLE PRECISION AMAX
      DOUBLE PRECISION DMU,EORB
      DOUBLE PRECISION ETAR(MXCH)
      INTEGER I,IA,IPAR,ITEST
      INTEGER ITESTX,J,J1,J2
      INTEGER J3,JAG,JCHANX,JLEV
      INTEGER JT,K,K1,KAPA
      INTEGER KAPAJ,KJ,KXMAX,KXMIN
      INTEGER LAG,MK,NCHAN1
      INTEGER NJ,NL
      INTEGER JINTAR(MXCH),JINTJ(MXCH),JINSPAR(MXCH),
     +JINK(MXCH),JINNJ(MXCH),JINNL(MXCH),INDX(MXCH),
     +NC(MXCH),NCTAR(MXCH,MXCH)
      DOUBLE PRECISION ECONT(MXCH),DMIU(MXCH)      
!     rui jin
      character(len = 10 ) :: how_to_calc_orb
      character(len = 4) :: xo, xt, ALS
      character(len = 30 ) :: ATMP
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C   This function tests the triangular condition.
C   The arguments are 2J.
      ITRG(J1,J2,J3) = ABS(J1-J2).LE.J3.AND.J3.LE.J1+J2.AND.MOD(J1+J2+J3!
     +,2).EQ.0
C
      COMMON /XORBS/XEORB(MXCH),IQDT(MXCH),IEORB
      DATA LAB/'s ','p-','p ','d-','d ','f-','f ','g-','g ','h-','h ',  !
     +   'i-','i ','j-','j ','k-','k ','l-','l ','m-','m ','  '/

      interface !from_astate
        character(len=4) function from_astate(As)
        character(len=30) :: As
        end function from_astate
      end interface 
C-----------------------------------------------------------------------
      iwrite=999
      open(iwrite,file='TEST.OUT')
      OPEN(777,FILE='Channel.csv')
C-----------------------------------------------------------------------
C
C   Determine the channels allowed in this symmetry.
C
C   NCHAN    number of channels
C   NTARG(I) target level for channel I
C   K2P(I)   kappa-value for channel I
C   NAKK(I)  K-value for channel I
C
C   K2P  : -1, 1,-2, 2,-3, 3, ...
C   NAKK :  1, 2, 3, 4, 5, 6, ...
C
      NCHAN = 0
      NCHAN1 = 0
      KXMIN = -1
      KXMAX = -1
C
      DO J = 1,NAST
      !write(*,*)'nast',j
C
C   For each target level J determine the 2J value and parity.
C
cccccc  jinrui      
        IA=J      
        JLEV = ITJPO(IA)
C
        IF (ISPAR(IA).LT.0) THEN
          ITEST = 1
        ELSE
          ITEST = 0
        ENDIF
C
C   Examine the kappa values.
C
C   J2P1 is 2J+1 for the continnum state
C   JLEV is 2J for the target state
C   JAG  is 2J for the continuum orbital
C
        K1 = J2P1+JLEV
        DO K = 1,K1
C
          IF (MOD(K,2).EQ.1) THEN
            KAPA = -(K+1)/2
          ELSE
            KAPA = K/2
          ENDIF
C
          MK = ABS(KAPA)
          JAG = MK+MK-1
C
C  Test the triangular condition.
C
          IF (ITRG(JLEV,JAG,J2P1-1)) THEN
C
            LAG = MK+(SIGN(1,KAPA)-1)/2
            ITESTX = ITEST+LAG
C
            IF (MOD(ITESTX,2).NE.0) THEN
              IPAR = -1
            ELSE
              IPAR = 1
            ENDIF
C
C  Test the parity condition.
C
            IF (IPAR.EQ.NPTY) THEN
C
              IF (KXMIN.LT.0) THEN
                KXMIN = K
                KXMAX = K
              ELSE
                IF (K.GT.KXMAX) THEN
                  KXMAX = K
                ELSE
                  IF (K.LT.KXMIN) THEN
                    KXMIN = K
                  ENDIF
                ENDIF
              ENDIF
              write(999,*)'K',K,'KAPA',KAPA
C
              IF (K.GT.KCMAX .OR. K.LT.KCMIN) THEN
C
C
                NCHAN1 = NCHAN1+1
                write (Iwrite,3080) J,KAPA
C
              ELSE
C
                NCHAN = NCHAN+1
C
c                IF (NCHAN.LE.IDMTST(1)) THEN
                  NTARG(NCHAN) = J
                  K2P(NCHAN) = KAPA
                  NAKK(NCHAN) = K
                  IF (K.GT.21) THEN
                    NHC(NCHAN) = LAB(22)
                  ELSE
                    NHC(NCHAN) = LAB(K)
                  ENDIF
c                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
      IF (NCHAN1.GT.0) write (Iwrite,3010) NCHAN1
      write (Iwrite,3020) KXMIN,KXMAX
c     jinrui 20141224
C
      IF (NCHAN.EQ.0) THEN
        write (Iwrite,3070)
        RETURN
      ENDIF
C
C   Write out the data defining the channels.
C
      JCHANX = J2P1-1
c     jinrui 
      write(777,'(4a8)')'Nchan','2*J_TOT','NPTY'
      write(777,'(4I8)')NCHAN,JCHANX,NPTY
      write(*,'(4a8)')'Nchan','2*J_TOT','NPTY'
      write(*,'(4I8)')NCHAN,JCHANX,NPTY
      
C
      IF (MOD(JCHANX,2).EQ.0) THEN
        JCHANX = JCHANX/2
        IF (NPTY.EQ.1) THEN
          write (Iwrite,3030) NCHAN,JCHANX
        ELSE
          write (Iwrite,3040) NCHAN,JCHANX
        ENDIF
      ELSE
        IF (NPTY.EQ.1) THEN
          write (Iwrite,3090) NCHAN,JCHANX
        ELSE
          write (Iwrite,3100) NCHAN,JCHANX
        ENDIF
      ENDIF
C
      write (Iwrite,3110)
      write(777,3130)
      write(*,*)'nchan',nchan
      DO I = 1,NCHAN
        K = K2P(I)
        MK = ABS(K)
        NJ = MK+MK-1
        NL = MK+(SIGN(1,K)-1)/2
        ITAR=NTARG(I)
        JINTAR(I)=ITAR
        JINTJ(I)=ITJPO(ITAR)
        JINSPAR(I)=ISPAR(ITAR)
        JINK(I)=K
        JINNJ(I)=NJ
        JINNL(I)=NL
        NC(ITAR)=NC(ITAR)+1
        NCTAR(ITAR,NC(ITAR))=I
        CALL FINDQDT(K,NL,QTAR,DMU,EORB,how_to_calc_orb)
        ECONT(I)=ETAR(ITAR)+EORB
        DMIU(I)=DMU
        write (Iwrite,3050) I,ITAR,ITJPO(ITAR)-1,
     +  ISPAR(ITAR),NHC(I),K,NJ,NL,ETAR(ITAR)+EORB
c      jinrui 20141224
        write(777,3121)
     +I,ITAR,ITJPO(ITAR),ISPAR(ITAR),ETAR(ITAR),NHC(I),K,NJ,NL,
     +ETAR(ITAR),EORB,DMU        
      ENDDO
      write(*,*)'channel count done '
      OPEN(666,FILE='Channel-order.csv')
      write(666,3130)
      CALL POP(ECONT,NCHAN,INDX)
      CALL INV(INDX,INVINDX,NCHAN)
 3120 FORMAT(4(I4,','),F15.4,',',A,',',A,',',A8,',',4(I8,','),F15.4,',')
      DO I=1,NCHAN
        J=INDX(I)
        write(xt, "(F4.1)") float(JINTJ(J))/2.0
        write(xo, "(F4.1)") float(JINNJ(J))/2.0
        ATMP = ASTAT(JINTAR(J))
        ALS = from_astate(ATMP)
        print "(A,1x,A,2i3)", trim(adjustl(ATMP)), trim(ALS), 
     &   JINTJ(J),
     &   JINSPAR(J)
        write(666,3120)I,JINTAR(J),JINTJ(J),JINSPAR(J),
     +  ETAR(JINTAR(J)),trim(adjustl(ATMP)),
     &trim(adjustl(ALS))//
     &trim(adjustl(xt))//
     &trim(NHC(J)(1:1))//
     &trim(adjustl(xo)),
     &  NHC(J),JINK(J),
     +  JINNJ(J),JINNL(J),J,ECONT(I)
      ENDDO
      CALL TAB(ECONT,JINTAR,JINTJ,JINSPAR,ETAR,NHC,JINK,
     +JINNJ,JINNL,DMIU,NCHAN,INDX,NAST,NC,NCTAR,
     +ITJPO,ISPAR,ASTAT,lcf,QTAR,INDX)
C-----------------------------------------------------------------------
 3121 FORMAT(4(I4,','),F15.4,',',A8,',',3(I8,','),3(F15.4,','))      
 3130 FORMAT('CHANNEL,TARG,2JT,PT,ETARG,ORB,KAPc,2Jc,lc,id,e,MIU')
 3010 FORMAT (/1X,' The number of channels excluded was ',I5)
 3020 FORMAT (/1X,' The K-value ranges from ',I4,' to ',I4)
 3030 FORMAT (/1X,I5,' channels generated with  J = ',I2,               !
     +' and even parity')
 3040 FORMAT (/1X,I5,' channels generated with  J = ',I2,               !
     +' and odd parity')
 3050 FORMAT (1X,4I4,8X,A2,2I4,'/2',I4,F15.4)
 3060 FORMAT (/1X,I5,' continuum CSFs generated')
 3070 FORMAT (/' ***********************************************'/      !
     +' ***          WARNING from CONT              ***'/               !
     +' *** No channels generated for this symmetry ***'/               !
     +' ***********************************************')
 3080 FORMAT (/                                                         !
     +' ***********************************************************'/   !
     +' ***              WARNING from CONT                      ***'/   !
     +' *** Channel with level ',I4,' and kappa ',I4,' not included ***'!
     +/' ***********************************************************')
 3090 FORMAT (/1X,I5,' channels generated with  J = ',I2,               !
     +'/2 and even parity')
 3100 FORMAT (/1X,I5,' channels generated with  J = ',I2,               !
     +'/2 and odd parity')
 3110 FORMAT (/10X,'target'/10X,'state ',10X,'K',3X,'J',5X,'L'/)
      END
      character(len=4) function from_astate(astate) 
      implicit none
      character(len = 30), intent(in) :: astate
      integer :: ibeg, iend
      ibeg = index(astate, '_')
      iend = index(astate(ibeg+1:len(astate)), '_') -1 + ibeg
      if(astate(iend:iend) == '*') then
        from_astate = astate(ibeg+1:iend-1)
      else 
        from_astate = astate(ibeg+1:iend)
      end if
      !print *, astate,ibeg, iend, from_astate
      return
      end function from_astate

      subroutine pop(f,n,indx)
      implicit real*8 (a-h,o-z)
      parameter (mm=1000)
      dimension f(mm),indx(mm)
      do i=1,n
        indx(i)=i
      enddo
      do 100 i=1,n
      do 100 j=i,n
        if(f(i).gt.f(j))then
          tmp=f(j)
          ik=indx(j)
          f(j)=f(i)
          indx(j)=indx(i)
          f(i)=tmp
          indx(i)=ik

        endif
100   continue
      end subroutine
      
      SUBROUTINE TAB(ECONT,JINTAR,JINTJ,JINSPAR,ETAR,NHC,JINK,
     +JINNJ,JINNL,DMIU,NCHAN,INVINDX,NAST,NC,NCTAR,
     +ITJPO,ISPAR,ASTAT,lcf,QTAR,INDX)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'darc.inc'
      DIMENSION ITJPO(MXCH),ISPAR(MXCH),KC(MXCH),ETAR(MXCH),
     +QDT(MXCH),NSRT(MXCH)
      INTEGER JINTAR(MXCH),JINTJ(MXCH),JINSPAR(MXCH),INDX(MXCH),
     +JINK(MXCH),JINNJ(MXCH),JINNL(MXCH),INVINDX(MXCH),
     +NC(MXCH),NCTAR(MXCH,MXCH),IPOS(MXCH,MXCH),lcf(mxch)
      DOUBLE PRECISION ECONT(MXCH),DMIU(MXCH)
      CHARACTER*2 NHC(MXCH)
      CHARACTER ASTAT*30(MXCH)
      COMMON /XQDT/QDT,NKC,NSRT
      COMMON /XORBS/XEORB(MXCH),IQDT(MXCH),IEORB
      OPEN(888,FILE='channel-open.csv')
1     FORMAT(5(','),<NCHAN>(F10.4,','))      
      write(888,1)(ECONT(I),I=1,NCHAN)
      write(888,1)(ECONT(I)/QTAR**2,I=1,NCHAN)
2     FORMAT(I4,',',A<lentarg>,',',2(I4,','),F10.4)
3     FORMAT(I4,',',A<lentarg>,',',2(I4,','),F10.4)
!4     FORMAT(<4+IPOS(I,I2)>(','),A)
4     FORMAT(<IPOS(I,I2)>(','),A)
      
      DO I=1,NCHAN
        J=invindx(I)
        JTA=JINTAR(J)
        INC=NC(JTA)
c        write(*,*)'jta',jta,'inc',inc
        DO I2=1,INC
          ICH=NCTAR(JTA,I2)
          IF(ICH.EQ.J)then
            IPOS(JTA,I2)=I
c            write(*,*)'   ich',ich,'j',j,'ipos',IPOS(JTA,I2)
          endif
        ENDDO
      ENDDO

      DO I=1,NAST
        INC=NC(I)
        lentarg=lcf(i)
        IF(INC.GT.0)THEN
cc         write(*,*)'TARG CHAN  POS'
c          DO I2=1,INC
c            write(*,'(3I5,3x,A)')I,I2,IPOS(I,I2),NHC(NCTAR(I,I2))
c          ENDDO
          write(888,3, advance = 'no')
     &       I,ASTAT(I),ITJPO(I),ISPAR(I),ETAR(I)
          DO I2=1,INC
            write(888,4)NHC(NCTAR(I,I2))
          ENDDO
        ELSE
          write(888,2)I,ASTAT(I),ITJPO(I),ISPAR(I),ETAR(I)
        ENDIF
      ENDDO
      CLOSE(888)
      END
       subroutine nist(ix,Acon,E,xj,np,lcf)
       implicit real*8 (a-h,o-z)
       parameter (MXL=1000)
       character*2 AELEM
       character ACONFLS*50, J*5, AHEAD*1,ACON*30(MXL)
       dimension lcf(mxl),E(mxl),mainl(mxl),XJ(mxl),NP(mxl)
       common /state/ACONFLS, J
       common /length/lenls,len3,lenconf
       common /elemt/ion,AELEM
       open(2,file='test.out')
       call FDLEV(ACON,E,lcf,ix,mainl,XJ,NP)
       open(333,file='target.in')
       write(333,*)ix
       do i=1,ix
         write(333,33)ACON(i),E(i),xj(i),np(i)
       enddo
       close(333)
33     format(A<lcf(i)>,E20.10,F7.4,i3)       
       end

       subroutine FDLEV(ACON,E,lcf,ix,mainl,XJ,NP)
       implicit real*8 (a-h,o-z)
       parameter (MXL=1000)
       character*2 ELEM(9),AELEM
       character AFILE*10,ATMP*20
       character ACON*30(MXL)
       dimension E(mxl),lcf(mxl),mainl(mxl),XJ(mxl),NP(mxl)
       logical EX
       common /elemt/ion,AELEM
       common /state/ACONFLS, J
       write(ATMP,*)AELEM
       LEN1=len_trim(adjustl(ATMP))
       AFILE(1:LEN1)=adjustl(ATMP)
       LEN2=LEN1+1
       AFILE(LEN1+1:LEN2)='.'
       LEN1=LEN2
       write(ATMP,'(i3)')ion
       LEN2=len1+len_trim(adjustl(ATMP))
       AFILE(LEN1+1:LEN2)=adjustl(ATMP)
       AFILE(LEN2+1:LEN2+5)='.NIST'
       inquire(file=AFILE,exist=EX)
       if(EX)then
       open(1,file=AFILE)
       call READLINE(ACON,E,lcf,ix,mainl,XJ,NP)
       endif
       end
       
       
       subroutine READLINE(ACON,E,lcf,ix,mainl,XJ,NP)
       implicit real*8 (a-h,o-z)
       parameter (MXL=1000)
       character*300 TERM
       character AHEAD*1,ACON*30(MXL), ATMP * 30 
       dimension ispec(6),inum(2),icpl(2),isml(2)
       dimension E(mxl),lcf(mxl),mainl(mxl),XJ(mxl),NP(mxl)
       data ISPEC/32,63,40,41,45,124/       
C---------------------------------------------------------------
C             SPACE, ?  (  )   -  |
C---------------------------------------------------------------       
       common /state/ACONFLS, J
c       common /length/lenls,len3,lenconf
      
11     format(A300)         
       iline = 0   
       ix=0
       mainl=1
12     read(1,11,err=77,end=77)TERM
       iline = iline +1
       AHEAD=TERM(1:1) 
       ihead=ich(AHEAD)
c       write(*,*)AHEAD,ihead, iline
       if(ihead.eq.45)goto 12
       if(ihead.eq.8)goto 12
       if(ihead.eq.9)goto 12
       if(ihead.eq.7.or.ihead.eq.32)then
         call READCOMP(term,ix,E,ACON,lcf,mainl,XJ,NP)
         goto 12
c   reading complete lines
       else
         write(*,*)'There might be something wrong with NIST.IN'
       endif
77     write(*,*)'read done!'     
       end
      subroutine readcomp(term,ix,E,ACON,lcf,mainl,XJ,NP)
      implicit real*8 (a-h,o-z)
      parameter (MXL=1000)
      character*300 TERM
      character ACONFLS*50, J*5, AHEAD*1,Acon*30(mxl),AP*5, ATMP * 30
      dimension ispec(6),inum(2),icpl(2),isml(2),XJ(mxl),NP(mxl)
      dimension E(mxl),lcf(mxl),MAINL(MXL)
      data ISPEC/32,63,40,41,45,124/
C---------------------------------------------------------------
C             SPACE, ?  (  )   -  |
C---------------------------------------------------------------
      common /state/ACONFLS, J
      common /length/lenls,len3,lenconf
      imak=0  
      LEN1=1
      LEN2=0
      len3=0
      i=1
      lend=0
      E1=0.0
      idec=0
      idot =0
      AMTP =''
      do while(imak.le.3.and.i.le.300)
        AHEAD=TERM(I:I)
        if(imak.lt.2)then
          if(AHEAD.EQ.'|')then
            imak=imak+1
            goto 44
          elseif(AHEAD.NE.'|'.and.AHEAD.NE.' ')then
            len2=len2+1
            ACONFLS(LEN2:LEN2)=AHEAD
            goto 44
          elseif(len2.gt.1.and.AHEAD.eq.' '.and.lend.lt.1)then
            len2=len2+1
            ACONFLS(LEN2:LEN2)='_'
            lend=1
            goto 44
          endif
        elseif(imak.lt.3)then
          if(AHEAD.EQ.'|')then
            imak=imak+1
            goto 44
          elseif(AHEAD.NE.'|'.and.AHEAD.NE.' ')then
            len3=len3+1
            J(LEN3:LEN3)=AHEAD
            goto 44
          endif  
        elseif(len3.ge.1)then
          if(ahead.eq.'('.or.ahead.eq.'?'.or.ahead.eq.'+')then
            goto 55
          elseif(AHEAD.ne.' '.and.AHEAD.ne. '|'.and.idot.lt.1
     +    .and.AHEAD.ne.'.')then
            E1=E1*10.0+dble(ICHAR(AHEAD)-ichar('0'))
          elseif(AHEAD.eq.'.')then
            idot =idot+1
          elseif(idot.gt.0.and.AHEAD.ne.'|'.and.AHEAD.ne.'
     +    '.and.AHEAD.ne.'.')then
            idec=idec+1
            E1=E1+dble(ICHAR(AHEAD)-ichar('0'))/10.0**idec
          elseif(AHEAD.eq.'|')then
            imak =imak+1
            goto 44
          endif
        endif        
   44   i=i+1
      enddo
55    if(len3.ge.1)then
        ix=ix+1
        if(len2.gt.0)then
          lenls=len2
        else
          mainl(ix)=0
        endif      
        aconfls=adjustl(aconfls)
        j=adjustl(j)
        lenls1=lenls+1
        aconfls(lenls1:lenls1)='_'
        if(aconfls(lenls1-1:lenls1-1).eq.'*')then
          NP(ix)=-1
        else
          NP(ix)=1
        endif
        lenconf=lenls1+len3
        aconfls(lenls1+1:lenconf)=J        
C        write(*,*)J,'len3',len3
        XJ(ix)=XA2J(J(1:len3),len3)
        acon(ix)(1:lenconf)=ACONFLS(1:lenconf)
        print '(A)', trim(acon(ix))
        E(ix)=E1
        lcf(ix)=lenconf
      endif
      end
      function ICH(A)
      implicit real*8 (a-h,o-z)
      parameter(isec = 9)
      character*1 A
      dimension ispec(6),inum(2),icpl(2),isml(2)
      data ISPEC/32,63,40,41,45,124/
C--------------------------------------------------------------
C           SPACE,     ?       (        )        -        |     
C---------------------------------------------------------------       
      data INUM/49,57/,ICPL/65,90/,ISML/97,122/
C--------------------------------------------------------------
C           NUMBER     CAPITAL LETTER  SMALL LETTER
      ia=ichar(A)
      do i=1,6
        i2=ispec(i)
        if(A.eq.char(I2))then
          ICH=I2  
        endif
      enddo
      if(IA.GE.INUM(1).and.IA.LE.INUM(2))then
        ICH=7 
      else if(IA.GE.ICPL(1).and.IA.LE.ICPL(2))then
        ICH=8
      else if(IA.GE.ISML(1).and.IA.LE.ISML(2))then
        ICH=9
      endif
      end
      function XA2J(BB,L)
      implicit real*8 (a-h,o-z)
      parameter(isec = 9)
      character AA*5,A*1(5),BB*5
      dimension ispec(6),inum(2),icpl(2),isml(2)
      AA=adjustl(BB(1:L))
      if(L.eq.1)then
         A(1)=AA(1:1)
         XA2J=dble(ichar(A(1))-ichar('0'))
      elseif(L.eq.2)then
         A(1)=AA(1:1)
         A(2)=AA(2:2)
         XA2J=dble(ichar(A(1))-ichar('0'))*10.0+
     +   dble(ichar(A(2))-ichar('0'))
      elseif(L.eq.3)then
         A(1)=AA(1:1)
         XA2J=DBLE(ichar(A(1))-ichar('0'))/2.0
      elseif(L.eq.4)then
         A(1)=AA(1:1)
         A(2)=AA(2:2)
         XA2J=DBLE(ichar(A(1))-ichar('0'))*10.0+
     +   DBLE(ichar(A(2))-ichar('0'))
         XA2J=XA2J/2.0
      endif
      end
      SUBROUTINE INV(INDX,INVINDX,NCHAN)
      implicit real*8(a-h,o-z)
      INCLUDE 'darc.inc'
      dimension INDX(mxch),INVINDX(MXCH)
      do i=1,nchan
      J=INDX(I)
      INVINDX(J)=I        
      enddo
      !write(*,'(20I3)')(indx(i),i=1,nchan)
      !write(*,'(20I3)')(invindx(i),i=1,nchan)
      
      end
      subroutine findK(kapa,k,nkc)
      implicit real*8(a-h,o-z)
      INCLUDE 'darc.inc'
      dimension kc(mxch)
      DO IKC=1,NKC
        II=MOD(IKC,2)
        IF(II.EQ.1)THEN
          KC(IKC)=-(IKC/2+1)
        ELSEIF(II.EQ.0)THEN
          KC(IKC)=(IKC/2)
        ENDIF
      ENDDO
      do k=1,nkc
      if(kapa.eq.kc(K))then
        return
      endif
      enddo      
      end
      subroutine makerscf
c
c   input file: laputa.in, configurations
c  output file: laputa.out, database of transition energy ,probobility
c               and variance
c 
c-----------------------------------------------------------------------
c ==================== begin program ======================================
 
      implicit real*8(a-h,o-z)
      parameter (nmx=2000,jm=200,mxline=1000)
      parameter (mxion=95,ntwomx=10000)
      parameter (mxphoto=5000,mxnk=50)

c ... set the maximum number of:
c        nmx, radius mesh points;  jm, bound electron orbitors;
c        mxion, ionization degree for an element (i.e. nuclear charge )

      character symb*2,symbn1*2
      character*4 wrd(70),dsym(70)*4
      character aac1*4,aac2*2,aac3*1
      character symfig*160

      dimension azc(jm),nncc(jm),kkcc(jm),nnzz(jm)
      dimension dd0(mxline),dds(mxline,jm),wwd(mxline,jm)
     *         ,ffo(mxline)
      dimension nayi(100),kayi(100),nayf(100),kayf(100)
c gaoxiang 2015.12.24 add to resolve bf output
      dimension spetbf(mxphoto),ffkap(mxphoto),spetbff(mxion,jm,mxphoto)
      dimension use1(mxphoto),use2(mxphoto),use3(mxphoto),
     *          use4(mxphoto),use5(mxphoto),fracsp(10)

      dimension xx(mxnk),yy(mxnk),w(mxnk),
     *          c(mxnk),iop(2),tab(3)
      dimension nncn1(jm),azcn1(jm),symbn1(jm)
      dimension occpp(27,20)

      common /tranar/ norbb,narr(70),karr(70),nppi(mxline),kppi(mxline),
     *   nppf(mxline),kppf(mxline),nartol,npii(mxline),kpii(mxline),
     *   npff(mxline),kpff(mxline),ntrm
c .........................................................................

      common /cmmn1/ aame(10),rn,h,z,zn,xion,phi,eps,del,delrvr,xalpha,
     1   xlattr,fex,rnuc,xn(jm),xl(jm),xj(jm),xe(jm),xz(jm),bcond,xnnn,
     2   fn,fl,e,fj,vr(nmx),q,ev,ed,ev3,er,ev4,ez,ee2,ey,v0,et3,summa,
     3   pvt,convrg,pv(jm),a(nmx),b(nmx),rhnlj(nmx),dens(nmx),y(nmx),
     4   vinden(jm),da(5),dv(5),a0(5),b0(5),vrx(nmx),xkoop,r1,rad(nmx),
     5   rh(nmx),kpot,n,j,nc1,ndbg(5),nprint
      common /nctrn/itype
      common /ipd /pip,vipd

c ..........................................................................

      common /card / nnc(jm),nzc(jm),llc(jm),kkc(jm),fjc(jm),symb(jm)
c ..........................................................................
      common /plsma/ znuc,amass,temper,rou,cDoppl
      common /eeee/ ephin,ephmx,epht(mxphoto),spetbb(mxphoto)
      common /mmmm/ mesh
      common /utauta/dd0,dds,wwd,ffo
      common /bnfbnf/ffcc(jm,50),eebf(jm,50),ecore(jm)
      common /utautt/norbt

      data pi,evcm,aut,avgdat/3.1415926d0,8056.479d0,27.21d0,6.02d+23/
      data wrd /' 1s ',' 2s ',' 2p ',' 2p+',' 3s ',' 3p ',' 3p+',' 3d ',
     * ' 3d+',' 4s ',' 4p ',' 4p+',' 4d ',' 4d+',' 4f ',' 4f+',
     * ' 5s ',' 5p ',' 5p+',' 5d ',' 5d+',' 5f ',' 5f+',' 5g ',' 5g+',
     * ' 6s ',' 6p ',' 6p+',' 6d ',' 6d+',' 6f ',' 6f+',' 6g ',' 6g+',
     * ' 7s ',' 7p ',' 7p+',' 7d ',' 7d+',' 7f ',' 7f+',' 7g ',' 7g+',
     * ' 8s ',' 8p ',' 8p+',' 8d ',' 8d+',' 8f ',' 8f+',' 8g ',' 8g+',
     * ' 9s ',' 9p ',' 9p+',' 9d ',' 9d+',' 9f ',' 9f+',' 9g ',' 9g+',
     * '10s ','10p ','10p+','10d ','10d+','10f ','10f+','10g ','10g+'/

      open(5,file='laputa.in',status='old')
      open(15,file='rscf.in',status='unknown',form='formatted')
      itype=1
c      read(5,*) einev,emxev,mesh,istrg
      ephmx=emxev/aut
      ephin=einev/aut
      eh=(ephmx-ephin)/dfloat(mesh)
      e0=ephin
      do i=1,mesh
         e0=e0+eh
         epht(i)=e0
         spetbb(i)=0.0
         spetbf(i)=0.0
      end do

c gaoxinag 2016.3.24 ground bf

      nr=421
      rnr=60.0d0
      hr=32.0d0
      phir=0.3d0
      epsr=75.0d0
      delr=0.000005d0
      delrrr=0.00001d0
      nc1=30
      iread=0
      do i=1,5
        ndbg(i)=0
      enddo
      nprint=0
      kpot=3
      xalpha=1.d0
      xlattr=1.d0
      fex=0.d0
      rnuc=0.d0

      znr=26.d0
      nzion=0
c      amass=55.845
c      iread=iiread

      nnlog = 15
c gaoxiang 2016.3.24 ground bf
      read(5,11) aame
      read(5,*) znr, nzion, norbt, amass
      read(5,9321) (nnc(i),symb(i),azc(i),i=1,norbt)
   11 format (10a8)
   16 format (3f5.1,f14.7,f7.4,20x,f4.0)
 9321 format(70(i2,a2,f6.4))

      rewind nnlog
       xzz=nzion
       write (nnlog,19) aame, nr,norbt,rnr,hr,znr,xzz,phir,epsr,delr,
     1   delrrr, nc1,iread,(ndbg(ind),ind=1,5),nprint,
     1   kpot,xalpha,xlattr,fex,rnuc


   19 format (10a8/2i5,f6.2,f4.0,2f5.0,4f10.7
     1  /8i5/i10,f10.7,10x,2f5.1,d10.3)

      do in=1,norbt
        call orbitor(llc(in),kkc(in),fjc(in),azc(in),symb(in))
        xnnc=nnc(in)
        xllc=llc(in)
        xfjc=fjc(in)
        xnzc=azc(in)
        xeec=-(znr/xnnc)**2/2.d0
        write (nnlog,16) xnnc,xllc,xfjc,xeec,xnzc
      enddo
      close(nnlog)
      end

      subroutine orbitor(lc,kc,fc,fnz,sym)
      implicit real*8(a-h,o-z)

c this subroutine determin kap,l,j,for symbol of electron
c
      character sym*2
      if(sym.eq.'s ') then
      lc=0
      kc=-1
      fc=0.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration s'
      return
      end if
c
      if(sym.eq.'p ') then
      lc=1
      kc=1
      fc=0.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration p'
      return
      end if
c
      if(sym.eq.'p+') then
      lc=1
      kc=-2
      fc=1.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration p+'
      return
      end if
c
      if(sym.eq.'d ') then
      lc=2
      kc=2
      fc=1.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration d'
      return
      end if
c
      if(sym.eq.'d+') then
      lc=2
      kc=-3
      fc=2.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration d+'
      return
      end if
c
      if(sym.eq.'f ') then
      lc=3
      kc=3
      fc=2.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration f'
      return
      end if
c
      if(sym.eq.'f+') then
      lc=3
      kc=-4
      fc=3.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration f+'
      return
      end if
c
      if(sym.eq.'g ') then
      lc=4
      kc=4
      fc=3.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration g'
      return
      end if
c
      if(sym.eq.'g+') then
      lc=4
      kc=-5
      fc=4.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration g+'
      return
      end if
c
      if(sym.eq.'h ') then
      lc=5
      kc=5
      fc=4.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration h'
      return
      end if
c
      if(sym.eq.'h+') then
      lc=5
      kc=-6
      fc=5.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration h+'
      return
      end if
c
      if(sym.eq.'i ') then
      lc=6
      kc=6
      fc=5.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration i'
      return
      end if
c
      if(sym.eq.'i+') then
      lc=6
      kc=-7
      fc=6.5
      nwin=2.*fc+1
      if(fnz.gt.nwin) stop '! error in configuration i+'
      return
      end if
      end
c
      double precision function cnm(n)
      implicit real*8 (a-h,o-z)
      if(n.eq.0) then
      cnm=1.0
      return
      end if
      nn=n
      sum=1.0
      do 10 i=1,nn
      sum=sum*i
   10 continue
      cnm=sum
      return
      end
c
      subroutine kaptran(k,kap)
      implicit real*8(a-h,o-z)
      dimension kap(3)
      if(k.eq.0) stop 'error in kaptran '
      if(k.gt.0) then
      kap(1)=k-1
      kap(2)=-k
      kap(3)=k+1
      else
      kap(1)=k+1
      kap(2)=-k
      kap(3)=k-1
      end if
      return
      end
      SUBROUTINE FINDQDT(K,NL,QTAR,DMU,EORB, how_to_calc_orb)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'darc.inc'
      DIMENSION QDT(MXCH),NSRT(MXCH)
      COMMON /XQDT/QDT,NKC,NSRT 
      COMMON /XORBS/XEORB(MXCH),IQDT(MXCH),IEORB
      character (len= 10) how_to_calc_orb
      IF((index(how_to_calc_orb, 'read') /= 0) .or.
     & (index(how_to_calc_orb, 'QDT') /= 0) )THEN
        DO I=1,NKC
          IF(K.EQ.IQDT(I))then
            DMU=QDT(I)
            XN=DBLE(NSRT(I))
            CYCLE
          ENDIF
        ENDDO
        EORB1=-QTAR**2/(XN-DMU)**2
        IF(NL.NE.0)THEN
          EORB2=-QTAR**2/DBLE(NL)**2
          EORB=MAX(EORB1,EORB2)
        ELSE
          EORB=EORB1
        ENDIF
      ELSE if (index(how_to_calc_orb, 'calc') /= 0) then
        DO I=1,NKC
          IF(K.EQ.IQDT(I))then
            EORB=XEORB(I)                
            CYCLE
          ENDIF
        ENDDO
      ENDIF
      END
      SUBROUTINE RSCFORB(NKC)
      USE IFPORT
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'darc.inc'
      CHARACTER RSCFP*50,RWFP*50
      DIMENSION XL(MXCH),XJ(MXCH),XNN(MXCH)
      COMMON /XORBS/XEORB(MXCH),IQDT(MXCH),IEORB
      DATA RSCFP/'rscfp>lk'/
      DATA RWFP/'rwfp>lk'/
      CALL MAKERSCF
      write(*,*)'make rscf done'
      III=SYSTEM(RSCFP)
      write(*,*)'rscfp done'
      CALL MAKERWFP
      III=SYSTEM(RWFP)
      OPEN(33,FILE='eband.txt')
      DO I=1,NKC
      READ(33,*)XNN(I),XL(I),XJ(I),XEORB(I)
      IQDT(I)=LJ2KAP(XL(I),XJ(I))
      ENDDO
      END
      FUNCTION LJ2KAP(XL,XJ)
      IMPLICIT REAL*8(A-H,O-Z)
      I=INT(2.0*XJ+1.0)
      ISGN=INT(2.0*XL-2.0*XJ)
      LJ2KAP=I*ISGN/2
      END
      SUBROUTINE MAKERWFP
c
c   input file: laputa.in, configurations
c  output file: laputa.out, database of transition energy ,probobility
c               and variance
c 
c-----------------------------------------------------------------------
c ==================== begin program ======================================
 
      implicit real*8(a-h,o-z)
      parameter (nmx=2000,jm=200,mxline=1000)
      parameter (mxion=95,ntwomx=10000)
      parameter (mxphoto=5000,mxnk=50)

c ... set the maximum number of:
c        nmx, radius mesh points;  jm, bound electron orbitors;
c        mxion, ionization degree for an element (i.e. nuclear charge )

      character symb*2,symbn1*2
      character*4 wrd(70),dsym(70)*4
      character aac1*4,aac2*2,aac3*1
      character symfig*160

      dimension azc(jm),nncc(jm),kkcc(jm),nnzz(jm)
      dimension dd0(mxline),dds(mxline,jm),wwd(mxline,jm)
     *         ,ffo(mxline)
      dimension nayi(100),kayi(100),nayf(100),kayf(100)
c gaoxiang 2015.12.24 add to resolve bf output
      dimension spetbf(mxphoto),ffkap(mxphoto),spetbff(mxion,jm,mxphoto)
      dimension use1(mxphoto),use2(mxphoto),use3(mxphoto),
     *          use4(mxphoto),use5(mxphoto),fracsp(10)

      dimension xx(mxnk),yy(mxnk),w(mxnk),
     *          c(mxnk),iop(2),tab(3)
      dimension nncn1(jm),azcn1(jm),symbn1(jm)
      dimension occpp(27,20)

      common /tranar/ norbb,narr(70),karr(70),nppi(mxline),kppi(mxline),
     *   nppf(mxline),kppf(mxline),nartol,npii(mxline),kpii(mxline),
     *   npff(mxline),kpff(mxline),ntrm
c .........................................................................

      common /cmmn1/ aame(10),rn,h,z,zn,xion,phi,eps,del,delrvr,xalpha,
     1   xlattr,fex,rnuc,xn(jm),xl(jm),xj(jm),xe(jm),xz(jm),bcond,xnnn,
     2   fn,fl,e,fj,vr(nmx),q,ev,ed,ev3,er,ev4,ez,ee2,ey,v0,et3,summa,
     3   pvt,convrg,pv(jm),a(nmx),b(nmx),rhnlj(nmx),dens(nmx),y(nmx),
     4   vinden(jm),da(5),dv(5),a0(5),b0(5),vrx(nmx),xkoop,r1,rad(nmx),
     5   rh(nmx),kpot,n,j,nc1,ndbg(5),nprint
      common /nctrn/itype
      common /ipd /pip,vipd

c ..........................................................................

      common /card / nnc(jm),nzc(jm),llc(jm),kkc(jm),fjc(jm),symb(jm)
c ..........................................................................
      common /plsma/ znuc,amass,temper,rou,cDoppl
      common /eeee/ ephin,ephmx,epht(mxphoto),spetbb(mxphoto)
      common /mmmm/ mesh
      common /utauta/dd0,dds,wwd,ffo
      common /bnfbnf/ffcc(jm,50),eebf(jm,50),ecore(jm)
      common /utautt/norbt

      data pi,evcm,aut,avgdat/3.1415926d0,8056.479d0,27.21d0,6.02d+23/
      data wrd /' 1s ',' 2s ',' 2p ',' 2p+',' 3s ',' 3p ',' 3p+',' 3d ',
     * ' 3d+',' 4s ',' 4p ',' 4p+',' 4d ',' 4d+',' 4f ',' 4f+',
     * ' 5s ',' 5p ',' 5p+',' 5d ',' 5d+',' 5f ',' 5f+',' 5g ',' 5g+',
     * ' 6s ',' 6p ',' 6p+',' 6d ',' 6d+',' 6f ',' 6f+',' 6g ',' 6g+',
     * ' 7s ',' 7p ',' 7p+',' 7d ',' 7d+',' 7f ',' 7f+',' 7g ',' 7g+',
     * ' 8s ',' 8p ',' 8p+',' 8d ',' 8d+',' 8f ',' 8f+',' 8g ',' 8g+',
     * ' 9s ',' 9p ',' 9p+',' 9d ',' 9d+',' 9f ',' 9f+',' 9g ',' 9g+',
     * '10s ','10p ','10p+','10d ','10d+','10f ','10f+','10g ','10g+'/

      open(5,file='lrwfp.in',status='old')
      open(15,file='rwf.in',status='unknown',form='formatted')
      itype=1
c      read(5,*) einev,emxev,mesh,istrg
      ephmx=emxev/aut
      ephin=einev/aut
      eh=(ephmx-ephin)/dfloat(mesh)
	e0=ephin
	do i=1,mesh
	  e0=e0+eh
	  epht(i)=e0
	  spetbb(i)=0.0
	  spetbf(i)=0.0
	end do

c gaoxinag 2016.3.24 ground bf

      nr=421
      rnr=60.0d0
      hr=32.0d0
      phir=0.3d0
      epsr=75.0d0
      delr=0.000005d0
      delrrr=0.00001d0
      nc1=30
      iread=0
      do i=1,5
        ndbg(i)=0
      enddo
      nprint=0
      kpot=3
      xalpha=1.d0
      xlattr=1.d0
      fex=0.d0
      rnuc=0.d0

      znr=26.d0
      nzion=0
c      amass=55.845
c      iread=iiread

      nnlog = 15
c gaoxiang 2016.3.24 ground bf
      write(*,*)'makerwfp'
      read(5,11) aame
      read(5,*) znr, nzion, norbt, amass
      read(5,9321) (nnc(i),symb(i),azc(i),i=1,norbt)
   11 format (10a8)
   16 format (3f5.1,f14.7,f7.4,20x,f4.0)
 9321 format(70(i2,a2,f6.4))

      rewind nnlog
      xzz=nzion
      dell=0.d0 
      RMAX1=120.d0
      NWF=norbt
      NWF0=0
      NDEBUG=0
      NMIN=0
      NFACT=0
      zxx=0.d0
      write(nnlog,5000)zxx,dell,RMAX1
      write(nnlog,5001)NWF,NWF0,NDEBUG,NMIN,NFACT
5000  FORMAT(8F10.4)      
5001  FORMAT(8I10)      
      RFA=0.0
      SCALE=0.0
      skli=-0.1E11
      do in=1,norbt
        call orbitor(llc(in),kkc(in),fjc(in),azc(in),symb(in))
        xnnc=nnc(in)
        xllc=llc(in)
        xfjc=fjc(in)
        xnzc=azc(in)
        xeec=-(znr/xnnc)**2/2.d0
        write (nnlog,5002)xnnc,xllc,xfjc,xeec,RWA,scale,skli
        write (*,5002)xnnc,xllc,xfjc,xeec,RWA,scale,skli
      enddo
5002  FORMAT(6F10.4,E20.10)
      close(nnlog)
      end
