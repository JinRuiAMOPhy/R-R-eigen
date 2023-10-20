      module numeric
      use constants
      contains 
C**TRANSFORM BETWEEN ORTHOGONAL MATRIX AND EULER'S ANGLE*************
      SUBROUTINE UANGCON(U,A,R,N,IND,DS,DC,ndimen,iq)
C  |---------------------------------------------|
C  |  NOTICE: VARIABLE DIMENSION IS USED: Ndimen |
C  |          THE CORRESPONDING ELEMENTS SHOULD  |
C  |          HAVE THE EXACTLY DIMENSION BETWEEN |
C  |          MAIN AND SUBROUTINE PROGRAM        |
C  |---------------------------------------------|
C    A: EULER'S ANGLE WITH SUBSCRIPTS (I,J);
C  IND: MATRIX-->ANGLE: IND=0; ANGLE-->MATRIX: IND=1
C    K: NUMBER OF ANGLE;
C    N: DIMENSION OF MATRIX;
C    R: SHOULD BE UNIT MATRIX AFTER TERMINATING THE PROCESS OF MATRIX-->ANGLE
C    U: ORTHOGONAL MATRIX
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(Ndimen,Ndimen),R(Ndimen,Ndimen),
     1 A(Ndimen,Ndimen),DS(Ndimen,Ndimen),DC(Ndimen,Ndimen)
      common  /E/EEN,ee1(ms)

      IF (IND.EQ.1) GOTO 40
      DO 5 I=1,N
      DO 5 J=1,N
c      write(6,'(5e12.7)')r(i,j)
5     R(I,J)=U(I,J)

      DO 30 J=N,2,-1
      T=1.d0
C  CALCULATE A(1,J)...A(J-1,J) OF ORTHOGONAL MATRIX "DO 10"
      DO 10 I=1,J-1
      DS(I,J)=R(J,I)/T
      if(dabs(DS(I,J)).gt.1.d0) then
       print*,ds(i,j),r(j,i),t
       if(dabs(DS(I,J))-1.d0.le.0.5d0) then
        tmp=ds(i,j)/dabs(ds(i,j))
        ds(i,j)=tmp
       endif
      endif
      A(I,J)=DASIN(DS(I,J))
      DC(I,J)=DCOS(A(I,J))
      T=T*DC(I,J)

C--IN CASE SOME ANGLE A(I,J)=90 DEGREE--|
        IF(DABS(DC(I,J)).LT.1D-10.AND.I.LT.J-1) THEN
        WRITE(*,*) 'A(',I,J,')-90 DEGREE<1D-10'
        DO 8 IA=I+1,J-1
        A(IA,J)=0
        DS(IA,J)=0
8       DC(IA,J)=1
        GOTO 12
        ENDIF
10      CONTINUE
12      CONTINUE

C--IN CASE THE DIAGONAL ELEMENT OF THE ORTHOGONAL MATRIX IS NEGATIVE--|
        IF(R(J,J).LT.0) THEN
        APIA=DASIN(0.5D0)*6
        A(J-1,J)=APIA-A(J-1,J)
        IF(A(J-1,J).GT.APIA) A(J-1,J)=A(J-1,J)-APIA*2.D0
        DC(J-1,J)=-DC(J-1,J)
        ENDIF

C  ELIMINATE A(1,J)...A(J-1,J) FROM ORTHOGONAL MATRIX "DO 20"
      DO 30 I=J-1,1,-1
c       write
        DO 20 M=1,N
      RMI=R(M,I)
        RMJ=R(M,J)
        R(M,I)=RMI*DC(I,J)-RMJ*DS(I,J)
20    R(M,J)=RMI*DS(I,J)+RMJ*DC(I,J)
30    CONTINUE
      DO 35 I=1,N
c      write(6,*)i,r(i,i)
      write(26,*)i,r(i,i)
      IF(dabs(R(I,I)).LT.0.9) THEN
      write(6,*)i,r(i,i)
        WRITE(*,*) 'ERROR: U-matrix in file:mchc.out NO ORTHOGONAL at'
        WRITE(*,*) 'energy point E=',EEN
cq        STOP
        iq = 1
      ENDIF
35    CONTINUE
      GOTO 80

40    CALL UNITMAT(U,N,ndimen)
      DO 60 J=2,N
      DO 60 I=1,J-1
      DC(I,J)=DCOS(A(I,J))
      DS(I,J)=DSIN(A(I,J))
      DO 50 M=1,N
        UMI=U(M,I)
        UMJ=U(M,J)
      U(M,I)=UMI*DC(I,J)+UMJ*DS(I,J)
      U(M,J)=-UMI*DS(I,J)+UMJ*DC(I,J)
      RMI=R(M,I)
        RMJ=R(M,J)
        R(M,I)=RMI*DC(I,J)+RMJ*DS(I,J)
50    R(M,J)=-RMI*DS(I,J)+RMJ*DC(I,J)
60    CONTINUE
80    RETURN
      END

C**UNIT MATRIX WITH Nth DIMENSION************************************
      SUBROUTINE UNITMAT(UM,N,ndimen)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION UM(Ndimen,Ndimen)
      DO 20 I=1,N
      DO 10 J=1,N
10    UM(I,J)=0.d0
20    UM(I,I)=1.d0
      RETURN
      END

C**SEQUENCING MIUalfa************************************************
      SUBROUTINE SEQUENC(U,DMIU,DMIU0,N)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ndimen=ms)
      !parameter (me=10000,ms=500,pai=3.1415926)
      DIMENSION U(NDIMEN,NDIMEN),DMIU(NDIMEN),ID(NDIMEN),DMIU0(NDIMEN),
     +DY(NDIMEN)
      common  /E/EEN,ee1(ms)
      COMMON /DU/DU(NDIMEN,NDIMEN)
      COMMON /DM/DD(NDIMEN)
      IF(N.gt.NDIMEN) STOP 'ERROR: NDIMEN in program should be change'

C--SEQUENCING ACCORDING ABOVE U MATRIX----|
      DO 5 I=1,N
      UM=0.
      DO 5 K=1,N
      UX=0.
      DO 3 J=1,N
3     UX=UX+U(J,K)*DU(J,I)
      IF(ABS(UX).GT.ABS(UM)) THEN
        UM=UX
        ID(I)=K
      ENDIF
5     CONTINUE

C--MARKING FROM LARGEST TO SMALLEST----|
c      DO 20 I=1,N
c      D=-10000.
c      DO 10 J=1,N
c      IF(D.LT.DMIU(J)) THEN
c        D=DMIU(J)
c        ID(I)=J
c      ENDIF
c10    CONTINUE
c      DMIU(ID(I))=-10000.
c20    CONTINUE

      DO 30 I=1,N
30    DY(I)=DMIU0(ID(I))
      DO 40 I=1,N
40    DMIU0(I)=DY(I)

      DO 50 I=1,N
      IF((DD(I)-DMIU(ID(I))).GT.0.5) THEN
        DMIU(ID(I))=DMIU(ID(I))+1.
      ELSEIF((DD(I)-DMIU(ID(I))).LT.-0.5) THEN
        DMIU(ID(I))=DMIU(ID(I))-1.
      ENDIF
50    DD(I)=DMIU(ID(I))
      DO 60 I=1,N
60    DMIU(I)=DD(I)

      NNT=0
61    NT=0

C--COLUMN----|
      DO 66 I=1,N
      DL=0.
      DO 63 J=1,N
63    DL=DL+DU(J,I)*U(J,ID(I))
      IF(DL.LT.0) THEN
        DO 65 K=1,N
65      U(K,ID(I))=-U(K,ID(I))
        NT=1
        NNT=NNT+1
      ENDIF
66    CONTINUE

C--ROW----|
      DO 69 I=1,N
      DH=0.
      DO 67 J=1,N
67    DH=DH+DU(I,J)*U(I,ID(J))
      IF(DH.LT.0) THEN
        DO 68 K=1,N
68      U(I,K)=-U(I,K)
        NT=1
        NNT=NNT+1
      ENDIF
69    CONTINUE

      IF(NT.EQ.1) THEN
        IF(NNT.GT.300) THEN
          WRITE(*,*) 'U-matrix data in file:mchc.out ERROR at'
          WRITE(*,*) 'energy point E=',EEN
          STOP
        ENDIF
        GOTO 61
      ENDIF

      DO 75 I=1,N
      DO 75 J=1,N
75    DU(I,J)=U(I,ID(J))
      DO 80 I=1,N
      DO 80 J=1,N
80    U(I,J)=DU(I,J)
      RETURN
      END
      SUBROUTINE MINV(U,N,W)
CC
C     ************************************************
C     SUBROUTINE FOR MARTRIX CALCULATION
C     N:THE ROW NUMBER OF MARTRIX.
C     W:THE VALUE OF THE MARTRIX.
C     ************************************************
      PARAMETER(ndimen=ms)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NDIMEN,NDIMEN)
      IF(N.GT.NDIMEN) STOP 'ERROR: NDIMEN IN SUB-PROGRAM TOO SMALL'
      w=1.
      l=1
      i=1
5     if(i.eq.n.or.w.eq.0) then
      else
        if (DABS(u(i,i)).gt.1D-30) then
        else
          j=i+1
6         if (j.gt.n.or.l.eq.0) then
          else
            if (Dabs(u(j,i)).le.1D-30) then
            else
              do 7 k=i,n
              t=u(i,k)
                u(i,k)=u(j,k)
                u(j,k)=t
7               continue
                w=-w
                l=0
             end if

             j=j+1
             goto 6
           end if
           if(l.ne.0) w=0
         end if
         if(w.eq.0) then
         else
              do 8 m=i+1,n
              x=u(m,i)/u(i,i)
              do 8 l=i+1,n
              u(m,l)=u(m,l)-x*u(i,l)
8              continue
         end if
         i=i+1
         goto 5
      end if
      if(w.eq.0) then
      else
        do 9 i=1,n
9        w=w*u(i,i)
      endif
c      write(*,*)'w',w
      return
      end

      subroutine xmatcor(u,v,a,n)
      implicit real*8 (a-h,o-z)
      !parameter (ms=500)
      dimension u(ms,ms),a(ms,ms),v(ms,ms)

      do i=1,n
        do j=1,n
         a(i,j)=0.d0
         do k=1,n
          a(i,j)=a(i,j)+u(k,i)*v(k,j)
         enddo
        enddo
      enddo
      return
      end

c*********************************************************************
      subroutine pivot(u,u2,mord,msgn,ns,ms)
      implicit real*8 (a-h,o-z)
      dimension u(ms,ms),mord(ms),msgn(ms),u2(ms,ms),mp(ms)
      do is=1,ns
      do js=1,ns
      mord(is)=is
      mp(is)=is
      msgn(is)=1
      u2(is,js)=u(is,js)
      enddo
      enddo

      do is=1,ns
      tmp=u2(is,is)
      do ks=is,ns
      if(dabs(u2(is,ks)) > dabs(tmp))then
        tmp=u2(is,ks)
        u2(is,ks)=u2(is,is)
        u2(is,is)=tmp
        mptmp=mp(is)
        mp(is)=mp(ks)
        mp(ks)=mptmp
        do ls=1,ns
          if(ls /= is)then
           tmp2=u2(ls,is)
           u2(ls,is)=u2(ls,ks)
           u2(ls,ks)=tmp2
          endif
        enddo
      endif
      enddo
      if(u2(is,is) < ZERO)then
        msgn(is) = -1
        do ls=1,ns
        u2(ls,is)=-u2(ls,is)
        enddo
      endif
      enddo
      do is=1,ns
c	mord(mp(is))=is
c	zdl order
        mord(is)=mp(is)
      enddo
      return
      end subroutine pivot

!--------- Angular momentum coupling ---
c********************************************************************
      real*8 function s9j(a1,a2,a3,b1,b2,b3,c1,c2,c3)
      implicit real*8(a-h,o-z)
c     calculate 9j symbol   a1  a2  a3
c                           b1  b2  b3
c                           c1  c2  c3
      jmax1=dmin1(a1+c3,c2+b1)*2.d0
      jmin1=dmax1(dabs(a1-c3),dabs(c2-b1))*2.d0
      jmax=dmin1(a2+b3,b1+c2)*2.d0
      jmin=dmax1(dabs(a2-b3),dabs(b1-c2))*2.d0
      jmax=min0(jmax1,jmax)
      jmin=max0(jmin1,jmin)
      s9j=0.d0
      do 10 j=jmin,jmax,2
      xj=dfloat(j)/2.d0
      phase=(-1.d0)**j
   10 s9j=s9j+phase*(2.d0*xj+1.d0)*s6j(a1,b1,c1,c2,c3,xj)
     1*s6j(a2,b2,c2,b1,xj,b3)
     2*s6j(a3,b3,c3,xj,a1,a2)
      return
      end

      subroutine sixj  (ia,ib,ic,id,ie,ig,sax)
          implicit real*8(a-h,o-z)
      common/fat/pfat(176)
c     note: do not use this subroutine]
c     use function s6j and racah
      l1= iabs(ib-ie)
      if(ia-l1) 4000,12,12
   12 if(ia-ib-ie) 13,13,4000
   13 l1= iabs(ic-ig)
      if(ia-l1) 4000,14,14
   14 if(ia-ic-ig) 15,15,4000
   15 l1= iabs(id-ic)
      if(ie-l1) 4000,16,16
   16 if(ie-id-ic) 17,17,4000
   17 l1= iabs(id-ig)
      if(ib-l1) 4000,18,18
   18 if(ib-id-ig) 19,19,4000
   19 l1=ia+ib+ie-((ia+ib+ie)/2)*2
      if(l1) 4000,20,4000
   20 l1=ia+ic+ig-((ia+ic+ig)/2)*2
      if(l1) 4000,21,4000
   21 l1=ic+id+ie-((ic+id+ie)/2)*2
      if(l1) 4000,22,4000
   22 l1=ib+id+ig-((ib+id+ig)/2)*2
      if(l1) 4000,23,4000
 4000 sax=0.
      return
 5000 sax=0.
      write(6,6005)
 6005 format('  n factorial exists,n>175 ')
      stop
   23 continue
c   mlst now calculate sax as it not zero.
c    lsi g page 99 of edmonds book.
      i1=(ia+ib-ie)/2+1
      i2=(ia-ib+ie)/2+1
      i3=(-ia+ib+ie)/2+1
      i4=(ia+ib+ie+2)/2+1
      d11=pfat(i1)
      d12=pfat(i2)
      d13=pfat(i3)-pfat(i4)
      d1=(d11+d12+d13)/2.d0
      i1=(ia+ic-ig)/2+1
      i2=(ia-ic+ig)/2+1
      i3=(-ia+ic+ig)/2+1
      i4=(ia+ic+ig+2)/2+1
      d11=pfat(i1)
      d12=pfat(i2)
      d13=pfat(i3)-pfat(i4)
      d2=(d11+d12+d13)/2.d0
      i1=(id+ib-ig)/2+1
      i2=(id-ib+ig)/2+1
      i3=(-id+ib+ig)/2+1
      i4=(id+ib+ig+2)/2+1
      d11=pfat(i1)
      d12=pfat(i2)
      d13=pfat(i3)-pfat(i4)
      d3=(d11+d12+d13)/2.d0
      i1=(id+ic-ie)/2+1
      i2=(id-ic+ie)/2+1
      i3=(-id+ic+ie)/2+1
      i4=(id+ic+ie+2)/2+1
      d11=pfat(i1)
      d12=pfat(i2)
      d13=pfat(i3)-pfat(i4)
      d4=(d11+d12+d13)/2.d0
      dd=d1+d2+d3+d4
      ss=0.
      k1=ia+ib+ie
      k2=ia+ic+ig
      k3=id+ib+ig
      k4=id+ic+ie
      k5=ia+ib+id+ic
      k6=ib+ie+ic+ig
      k7=ie+ia+ig+id
      kmax=min0(k5,k6,k7)
      kmin=max0(k1,k2,k3,k4)
      item=max0(kmin,k5,k6,k7)
      l1=kmin/2-1
      sign1=(-1.d0)**l1
      k=kmin
   25 continue
      sign1=-sign1
      in=(k+2)/2+1
      i1=(k-k1)/2+1
      i2=(k-k2)/2+1
      i3=(k-k3)/2+1
      i4=(k-k4)/2+1
      i5=(k5-k)/2+1
      i6=(k6-k)/2+1
      i7=(k7-k)/2+1
      t=pfat(in)-pfat(i1)-pfat(i2)-pfat(i3)-pfat(i4)-pfat(i5)-pfat(i6)-
     1  pfat(i7)+dd
      if(dabs(t).ge.175) go to 5000
      t=sign1*dexp(t)
      ss=ss+t
      k=k+2
      if(k-kmax) 25,25,24
   24 continue
      sax=ss
      return
      end
      real*8 function s6j(aa,ab,ae,ad,ac,ag)
c     calculate 6j symbol aa,ab,ae
c                         ad,ac,ag
      implicit real*8(a-h,o-z)
      data ii/0/
      if(ii.ne.0) go to 1
      ii=1
      dum=cgc(0d0,0d0,0d0,0d0,0d0,0d0)
    1 ia=2*aa
      ib=2*ab
      ic=2*ac
      id=2*ad
      ie=2*ae
      ig=2*ag
      call sixj(ia,ib,ic,id,ie,ig,s6s)
      s6j=s6s
      return
      end
      real*8 function  cgc(e1,e2,e3,q1,q2,q3)
      implicit real*8(a-h,o-z)
      common/fat/p(176)
      data ii /0/
      if(ii.eq.0) go to 3
    5 emax=e1+e2
      emin=dabs(e1-e2)
      if((e3.gt.emax).or.(e3.lt.emin)) go to 1
      if(q3.ne.(q1+q2)) go to 1
      j1=     (q3-e1+e2)
      j2=     (e3-e1+e2)
      j3=     (e3+q3)
      j1=max0(0,j1)
      j2=min0(j2,j3)
      b=0.
      j1=j1+1
      j2=j2+1
      do 2 jj=j1,j2
      j=jj-1
      is=j+(e2+q2)
      s=1.
      if(2*(is/2).ne.is) s=-1.
      i1=     (e3+e2+q1)-j+1
      i2=     (e1-q1)+j+1
      i3=     (e3+e2-e1)-j+1
      i4=     (e3+q3)-j+1
      i5=j+1
      i6=     (e1-e2-q3)+j+1
      t=p(i1)-p(i3)+p(i2)-p(i4)-p(i5)-p(i6)
      t=dexp(t)
    2 b=b+s*t
      i1=     (e3+e1-e2)+1
      i2=     (e3+e1+e2)+2
      a=p(i1)-p(i2)
      i1=     (e3-e1+e2)+1
      i2=     (e1-q1)+1
      a=a+p(i1)-p(i2)
      i1=     (-e3+e1+e2)+1
      i2=     (e1+q1)+1
      a=a+p(i1)-p(i2)
      i1=     (e3+q3)+1
      i2=     (e2-q2)+1
      a=a+p(i1)-p(i2)
      i1=     (e3-q3)+1
      i2=     (e2+q2)+1
      a=a+p(i1)-p(i2)
      a=a/2.
      b=b*dexp(a)
      a=(2.*e3+1.)
      cgc=dsqrt(a)*b
      return
    1 cgc=0.
      return
    3 p(1)=0.
      p(2)=0.
      do 4 i=3,176
    4 p(i)=p(i-1)+dlog(dfloat(i-1))
      ii=1
      go to 5
      end

      function geometryjjls9j(ns, jjls, LSJ, UgSign, chan_targcf) 
     &  result(Ug)
      use constants
      implicit none
      integer, intent(in) :: ns
      real(long), dimension(ms, 6), intent(in) :: jjls
      real(long), dimension(ms, 3), intent(in) :: LSJ
      real(long), dimension(:), intent(in) :: UgSign
      integer, dimension(:), intent(in) :: chan_targcf
      real(long), dimension(ms, ms) :: Ug
      integer :: is, js, iblck, ib, nh, istart, is2, js2
      integer, dimension(ms) :: imark
      real(long) :: a1, a2, a3, b1,b2, b3, c1, c2, c3
      real(long) :: c9j, xj3, xs2
      real(long), dimension(ms) :: xl1, xl2, xl3, xs1, xs3, xj1, xj2
      xs1(:) = jjls(:,1)
      xl1(:) = jjls(:,2)
      xj1(:) = jjls(:,3)
      xs2 = jjls(1,4)
      xl2(:) = jjls(:,5)
      xj2(:) = jjls(:,6)
      xs3(:) = LSJ(:,1)
      xl3(:) = LSJ(:,2)
      xj3 = LSJ(1,3)
      loop_targ: do is=1,ns
        if(is.eq.1)then
         iblck=1
        elseif(chan_targcf(is).ne.chan_targcf(is-1))then
         iblck=iblck+1
        endif
        imark(iblck)=is
      end do loop_targ

      do 200 ib=1,iblck
      if(ib.eq.1)then
      nh=imark(ib)
      istart=0
      else
      nh=imark(ib)-imark(ib-1)
      istart=imark(ib-1)
      endif
      do 200 is2=1,nh
      is=is2+istart
      do 200 js2=1,nh
      js=js2+istart
      if(xl1(is).eq.xl1(js).and.xl2(is).eq.xl2(js)
     &  .and.xs1(is).eq.xs1(js))then
      a1=xl1(is)
      a2=xs1(is)
      a3=xj1(is)
      b1=xl2(is)
      b2=xs2
      b3=xj2(is)
      c1=xl3(js)
      c2=xs3(js)
      c3=xj3
      c9j = s9j(a1,b1,c1,a2,b2,c2,a3,b3,c3)
      c9j = dsqrt((2*A3+1)*(2*B3+1)*(2*C1+1)*(2*C2+1))*C9J
      else
      c9j=0.d0
      endif
      ug(is,js)=c9j*UgSign(js)
200   continue
      end
      end module numeric
