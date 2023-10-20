      module assign_eigenchann 
      use constants
      use stdio
      implicit real*8 (a-h,o-z)
      !parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      parameter (ms1=500,ryd=13.6058)
      logical exham
      dimension xmiu(me,ms),dal(me,ms,ms),e(me),u(me,ms,ms),
     +dav(me,ms,ms),dli(me,ms,ms),iord(me,ms),xmiu2(me,ms),
     +dal2(me,ms,ms),u2(me,ms,ms),dav2(me,ms,ms),iord2(me,ms),
     +ie2(me),nsgn(ms),isgnin(ms),mords(ms),msgns(ms),s(ms,ms),
     +s2(ms,ms),xtmp(me,ms),utmp(me,ms,ms),dltmp(me,ms,ms),
     +dvtmp(me,ms,ms)
      dimension ug(ms,ms)
      data ie2/me*0/
      data ii/0/,xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /control2/iconnect,idescend
      common /ham/ham(ms,ms)
      common /exham/exham
      common /E2/ee1(ms1)

      contains 
      
      subroutine step2()

      call openfile
      xmiu=0.0
      dal=0.0
      dli=0.0
      ug=0.0
      ie=0
      istg=0
      call geometryjjls9j(ug)
      call putumiuda(ie,istg,ug,xmiu,dal,dli)
      ie=1
      if(exham)then
      write(129,*)'geometry U from eigen'
      call readham(ug)
      call putumiuda(ie,istg,ug,xmiu,dal,dli)
      endif
      istg=1
      call stringput
      call manualswap(ie2,iord)
      open(345,file='test')
      do 100 ie=1,ne
c     read in the sign output from smooth-pivot.
      read(134,'(50I10)')iee,(isgnin(i),i=1,ns)
      call readumiuda(xmiu,u,dal,dav,dli,ie,e)
      if(ie.eq.1)then
      call preput
      call smat(ug,u,ie,s)
      write(345,*)'s'
      do i=1,ns
      write(345,'(100E15.5)')(s(i,j),j=1,ns)
      enddo
      if(idescend.eq.0)then
      call pivot(s,s2,mords,msgns,ns,ms)
      else
      do i=1,ns
      msgns(i)=1
      mords(i)=i
      enddo      
      endif
      write(345,'(100(10x,I5))')(mords(i),i=1,ns)
      write(345,*)'s2'
      do i=1,ns
      write(345,'(100E15.5)')(s2(i,j),j=1,ns)
      enddo
      endif
c      write(*,*)'ie',ie
      call swap2(xmiu,xmiu2,u,u2,dal,dav,dal2,dav2,mords,msgns,ie)
c      write(*,*)'ie',ie
      write(345,*)'u2'
      do i=1,ns
      write(345,'(100E15.5)')(u2(ie,i,j),j=1,ns)
      enddo
      etot=e(ie)*xqp1**2
      if(ie2(ie).ne.0)then
      call swapumiuda(xtmp,xmiu2,utmp,u2,dltmp,dvtmp,
     +dal2,dav2,ie,iord,msgns)
      do i=1,ns
      xmiu2(ie,i)=xtmp(ie,i)
      do m1=1,mfin
      dal2(ie,i,m1)=dltmp(ie,i,m1)
      dav2(ie,i,m1)=dvtmp(ie,i,m1)
      enddo
      do j=1,ns
      u2(ie,i,j)=utmp(ie,i,j)
      enddo
      enddo
 
      endif
      write(345,*)'u2'
      do i=1,ns
      write(345,'(100E15.5)')(u2(ie,i,j),j=1,ns)
      enddo
      call putout(istg,ie,etot,xmiu2,u2,dal2,dav2,dli,ug,e)
      if(iconnect.eq.ie)call descend(u2,ie)
100   continue
      call putvmat(ug,u2,ie,etot)
      write(*,*)'putvmat done'
      call closefile

      end subroutine step2

      subroutine descend(u,ie)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926)
      dimension u(me,ms,ms)
      common /control/irad,mxe,nchop,ns,ne,mfin
      write(139,*)ns,ie
      do i=1,ns
      write(139,'(50f10.5)')(u(ie,i,j),j=1,ns)
      enddo
      close(139)
      end
      subroutine closefile
      close(123)
      close(124)
      close(125)
      close(126)
      close(127)
      close(128)
      end
      subroutine preput
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      character afile*72(ms,4),A1*4
      common /control/irad,mxe,nchop,ns,ne,mfin,mlast
      common /E2/ee1(ms1)
      do m1=1,mlast
      write(A1,'(I3)')m1
      A1=adjustl(A1)
      len=len_trim(a1)
      write(afile(m1,1),128)a1
      write(afile(m1,2),129)a1
      write(afile(m1,3),130)a1
      write(afile(m1,4),131)a1
      ifile=1000+m1
      jfile=1500+m1
      kfile=2000+m1
      mfile=2500+m1
128   format('DalphaL2-smooth.',A<len>,'.out')
129   format('DalphaV2-smooth.',A<len>,'.out')
130   format('Da2.',A<len>,'.out')
131   format('Da2.',A<len>,'.plot')
      open(ifile,file=afile(m1,1))
      write(ifile,'(2a20)')'E','DalphaL'
      open(jfile,file=afile(m1,2))
      write(jfile,'(2a20)')'E','DalphaV'
      open(kfile,file=afile(m1,3))
      open(mfile,file=afile(m1,4))
      enddo
      end
     
      subroutine putout(istg,ie,etot,xmiu2,u2,dal2,dav2,dli,ug,e)
      implicit real*8 (a-h,o-z)
      integer a
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      character ale(ms),aln(ms),alnp1(ms)
      character afile*72,A1*4
     1,name*13,ang*13(ms)
      dimension e(me),dli(me,ms,ms),iord(me,ms),xmiu2(me,ms),
     +dal2(me,ms,ms),u2(me,ms,ms),dav2(me,ms,ms),
     +utmp(ms,ms),iord2(me,ms),
     +ie2(me),nsgn(ms),isgnin(ms)
      dimension 
     +plus(ms),
     +ug(ms,ms),vmat(ms,ms),s(ms,ms),u3(ms,ms),a(ms)
      data ie2/me*0/
      data ii/0/,xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin,mlast
      common /E2/ee1(ms1)
      common /filename/afile(ms)
      if(irad.eq.1)then
      do m1=1,mfin
      ifile=1000+m1
      jfile=1500+m1
      kfile=2000+m1
      mfile=2500+m1
      write(126)(dal2(ie,ia,m1),ia=1,ns)
      write(126)(dav2(ie,ia,m1),ia=1,ns)
      write(126)(dli(ie,ia,m1),ia=1,ns)
      write(ifile,100)etot,(dal2(ie,ia,m1),ia=1,ns)
      write(jfile,100)etot,(dav2(ie,ia,m1),ia=1,ns)
     
      if(a(m1).eq.0)then
      write(mfile,*)mlast,m1,ee1(m1)
      a(m1)=1
      endif
      write(kfile,100)etot,(dal2(ie,ia,m1),ia=1,ns)
      write(mfile,100)etot,(dal2(ie,ia,m1),ia=1,ns)
      enddo
      endif

100   format(100E20.10)

101   format(A3,I3,100E20.10)      
188   format('E = ',F13.5,2x,'Ryd.',F13.5,2x,'eV',i4,'th energy points')
      write(129,188)etot,(e(ie)-EGRD)*xqp1**2*RYD,ie
107   format(16x,<ns>(2x,i3,5x))      
108   format(16x,<ns>(2x,i3,a,a,'_',i1,a),10x,A3)
      nprint = ns
      write(124)mxe,nchop
      write(124)e(ie)
      write(124)(xmiu2(ie,is)+plus(is),is=1,ns)
      call mat3copy2(ie,u2,u3)
      call putumiuda(ie,istg,u3,xmiu2,dal2,dli)
      if(irad.eq.1)then
      write(137,'(50E15.6)')etot,(dli(ie,is,1),is=1,ns)
      
      endif
      call putsmat(ug,u2,ie,etot,e)
      call mat3copy2(ie,u2,u3)
      call genmiuang(ie,u3,xmiu2,etot)
      end
      subroutine putsmat(ug,u2,ie,etot,e)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      character aeigchan*20(ms),aionchan*20(ms)
      dimension ug(ms,ms),u2(me,ms,ms),e(me),s(ms,ms),u3(ms,ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /channels/aeigchan,aionchan
      common /E2/ee1(ms1)
      call mat3copy2(ie,u2,u3)
      call xmatcor(ug,u3,s,ns)
      write(999,188)etot,(e(ie)-EGRD)*xqp1**2*RYD,ie
188   format('E = ',F13.5,2x,'Ryd.',F13.5,2x,'eV',i4,'th energy points')
      write(999,107)(is,is=1,ns)
107   format(11x,<ns>(2x,i3,5x))      
      do is=1,ns
      write(999,120)is,aeigchan(is),(s(is,js),js=1,ns)
120   format(i2,1x,A6,' |',<ns>(f10.5))
      enddo
      write(999,*)
      end
      subroutine smat(ug,u2,ie,s)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      character aeigchan*20(ms),aionchan*20(ms)
      dimension ug(ms,ms),u2(me,ms,ms),e(me),s(ms,ms),u3(ms,ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /channels/aeigchan,aionchan
      common /E2/ee1(ms1)
      call mat3copy2(ie,u2,u3)
      call xmatcor(ug,u3,s,ns)
      end
      subroutine
     +swapumiuda(xmiu2,xmiu,u2,u,dal2,dav2,dal,dav,ie,iord,isgnin)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      character ale(ms),aln(ms),alnp1(ms)
     1,name*13,ang*13(ms)
      dimension xmiu(me,ms),dal(me,ms,ms),e(me),u(me,ms,ms),
     +dav(me,ms,ms),dli(me,ms,ms),iord(me,ms),xmiu2(me,ms),
     +dal2(me,ms,ms),u2(me,ms,ms),dav2(me,ms,ms),
     +utmp(ms,ms),iord2(me,ms),
     +ie2(me),nsgn(ms),i2(ms),isgnin(ms)
      dimension sgnUG(ms),plus(ms),ugb(ms,ms)
      data ie2/me*0/
      data ii/0/,xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /E2/ee1(ms1)
      do 1 is=1,ns 
      xmiu2(ie,iord(ie,is))=xmiu(ie,is)
1     enddo
      do 3 js=1,ns
      do 2 is=1,ns
      u(ie,is,js)=u(ie,is,js)*isgnin(js)
2     enddo
      if(irad.eq.1)then
      do m1=1,mfin
      dal(ie,js,m1)=dal(ie,js,m1)*isgnin(js)
      dav(ie,js,m1)=dav(ie,js,m1)*isgnin(js)
      write(*,*)m1,mfin,ie,dal(ie,js,m1)
      enddo
      endif
3     enddo
      do 5 js=1,ns
      do 4 is=1,ns
      u2(ie,is,iord(ie,js))=u(ie,is,js)
4     enddo
      if(irad.eq.1)then
      do m1=1,mfin
      dal2(ie,iord(ie,js),m1)=dal(ie,js,m1)
      dav2(ie,iord(ie,js),m1)=dav(ie,js,m1)
      enddo
      endif
5     enddo
      do 7 is=1,ns
      if(u2(ie,is,is).lt.0.d0)then
      do 6 ks=1,ns
      u2(ie,ks,is)=-u2(ie,ks,is)
6     enddo
      if(irad.eq.1)then
      do m1=1,mfin
      dal2(ie,is,m1)=-dal2(ie,is,m1)
      dav2(ie,is,m1)=-dav2(ie,is,m1)
      enddo
      endif
      nsgn(is)=-1
      endif
7     enddo
      end
      subroutine copyumiuda(xmiu2,xmiu,u2,u,dal2,dav2,dal,dav,ie)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      dimension xmiu(me,ms),dal(me,ms,ms),u(me,ms,ms),
     +dav(me,ms,ms),dli(me,ms,ms),xmiu2(me,ms),
     +dal2(me,ms,ms),u2(me,ms,ms),dav2(me,ms,ms)
      dimension plus(ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /E2/ee1(ms1)
      do 1 is=1,ns 
      xmiu2(ie,is)=xmiu(ie,is)
1     enddo
      do 3 js=1,ns
      do 2 is=1,ns
      u2(ie,is,js)=u(ie,is,js)
2     enddo
      if(irad.eq.1)then
      do m1=1,mfin
      dal2(ie,js,m1)=dal(ie,js,m1)
      dav2(ie,js,m1)=dav(ie,js,m1)
      enddo
      endif
3     enddo
      end

      subroutine readumiuda(xmiu,u,dal,dav,dli,ie,e)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      dimension xmiu(me,ms),dal(me,ms,ms),e(me),u(me,ms,ms),
     +dav(me,ms,ms),dli(me,ms,ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin,mlast
      common /E2/ee1(ms1)
      read(123)mxe,nchop
      read(123)e(ie)
      een=e(ie)
      read(123)(xmiu(ie,is),is=1,ns)
      do is=1,ns
      read(123)(u(ie,is,js),js=1,ns)
      enddo
      if(irad.eq.1)then
      if(ie.eq.1)read(125)mlast
      read(125)etot,wdet,mfin,(ee1(m1),m1=1,mfin)
      write(*,*)'bbb',etot,wdet,mfin,(ee1(m1),m1=1,mfin)
      write(126)etot,wdet,mfin,(ee1(m1),m1=1,mfin)
      do m1=1,mfin
      read(125)(dal(ie,ialfa,m1),ialfa=1,ns)
      read(125)(dav(ie,ialfa,m1),ialfa=1,ns)
      read(125)(dli(ie,ialfa,m1),ialfa=1,ns)      
      enddo
      endif          
      end
      subroutine manualswap(ie2,iord)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      dimension iord(me,ms),iord2(me,ms),ie2(me),nsgn(ms),i1(me),i2(me)
      data ii/0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      read(127,*)nsec
      do isec=1,nsec
      read(127,*)i1(isec),i2(isec)
      read(127,*)(iord2(isec,is),is=1,ns)
      enddo
c      write(*,*)(iord2(1,is),is=1,ns)
      do isec=1,nsec
      do ii=i1(isec),i2(isec)
      ie2(ii)=1
      do is=1,ns
      iord(ii,is)=iord2(isec,is)
      enddo
      enddo
      enddo
      end
      subroutine anglab(nang)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058,mm=ms*ms)
      character anum1*3,anum2*3,ang*13(mm)
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /angles/ang
      nang=0
      do 10 ithet=1,Ns-1
      do 10 jthet=ithet+1,Ns
      write(anum1,'(i3)')ithet
      anum1=adjustl(anum1)
      lenthet1=len_trim(anum1)
      write(anum2,'(i3)')jthet
      anum2=adjustl(anum2)
      lenthet2=len_trim(anum2)          
      nang=nang+1
      write(ang(nang),789)anum1,'->',anum2      
   10 enddo
789   format(a<lenthet1>,a2,a<lenthet2>)        
      end
      subroutine stringput
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058,mm=ms*ms)
      character ale(ms),aln(ms),alnp1(ms)
     1,name*13,ang*13(mm)
      character aeigchan*20(ms),aionchan*20(ms)
      dimension xmiu(me,ms),dal(me,ms,ms),e(me),u(me,ms,ms),
     +dav(me,ms,ms),dli(me,ms,ms),iord(me,ms),A(ms,ms),xmiu2(me,ms),
     +dal2(me,ms,ms),R(ms,ms),u2(me,ms,ms),dav2(me,ms,ms),
     +utmp(ms,ms),DS(ms,ms),dc(ms,ms),iord2(me,ms),
     +ie2(me),nsgn(ms),isgnin(ms)
      dimension xl1(ms),xs1(ms),XJ1(ms),XL2(ms),XJ2(ms),
     +xl3(ms),xs3(ms),sgnUG(ms),plus(ms),ugb(ms,ms),
     +ug(ms,ms),vmat(ms,ms),s(ms,ms)
      data ie2/me*0/
      data ii/0/,xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /angles/ang
      common /channels/aeigchan,aionchan
      call anglab(nang)
      write(128,*)'E'
      write(128,*)'Ryd'
      if(irad.eq.1)then
      write(135,*)'E'
      write(136,*)'E'
      write(135,*)'Ryd'        
      write(136,*)'Ryd'        
      endif
109   format('- ',<ns>A,<ns*(ns-1)/2>A13)
110   format('- ',<ns>(7x,a))      
      if(ns.gt.1)then
c      write(*,*)(aeigchan(js),js=1,ns),(ang(jth),jth=1,nang)
      write(128,109)(aeigchan(js),js=1,ns),(ang(jth),jth=1,nang)
      else
      write(128,110)(aeigchan(js),js=1,ns)
      endif
      if(irad.eq.1)then
      write(136,110)(aeigchan(js),js=1,ns)
      write(135,110)(aeigchan(js),js=1,ns)
      endif
      end
      subroutine genmiuang(ie,u2,xmiu2,etot)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926)
      dimension A(ms,ms),xmiu2(me,ms),R(ms,ms),u2(ms,ms),
     +DS(ms,ms),dc(ms,ms),plus(ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      iq=0
c      write(*,*)'ie=',ie
      CALL UANGCON(U2,A,R,Ns,0,DS,DC,ms,iq)
      WRITE(128,'(500(1PE16.7))')etot,
     +(xmiu2(ie,is)+plus(is),is=1,ns),
     +((A(I,J)/pai,J=I+1,Ns),I=1,Ns-1)
      WRITE(133,'(I3,500(1PE16.7))') Ns,een,
     +(xmiu2(ie,is)+plus(is),is=1,ns)
     +,((A(I,J),J=I+1,Ns),I=1,Ns-1)
      end
      subroutine putvmat(ug,u,ie,etot)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      dimension ug(ms,ms),vmat(me,ms,ms),u(me,ms,ms)
      character aeigchan*20(ms),aionchan*20(ms),name*13
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /channels/aeigchan,aionchan
c	write up the vmat colomn in each file unit named as number//chan
c     into the output channel idisk=222
      idisk=222
      do 204 i=1,ns
      if(i.eq.1)then
      write(name,'(i2,a)')i,'st.Eigchan'
      elseif(i.eq.2)then
      write(name,'(i2,a)')i,'nd.Eigchan'
      elseif(i.eq.3)then
      write(name,'(i2,a)')i,'rd.Eigchan'
      else
      write(name,'(i2,a)')i,'th.Eigchan'
      endif
      open(idisk,file=name)
      write(idisk,'(100a15)')'E',(name,k=2,ns)
      write(idisk,'(a15)')'Ryd'
110   format('Ryd.',<ns>A)
      write(idisk,110)(aeigchan(js),js=1,ns)
      do 203 ie=1,ne
      do 202 is=1,ns
      vmat(ie,is,i)=0.d0
      do 202 ks=1,ns
      vmat(ie,is,i)=vmat(ie,is,i)+ug(ks,is)*u(ie,ks,i)
202   continue
      write(idisk,'(100E15.5)')etot,(vmat(ie,is,i),is=1,ns)
203   continue
      close(idisk)
204   continue
      end
      subroutine readham(ug)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      character ale(ms),aln(ms),alnp1(ms)
      dimension xl1(ms),xs1(ms),XJ1(ms),XL2(ms),XJ2(ms),
     +xl3(ms),xs3(ms),sgnUG(ms),ugb(ms,ms),vmat(ms,ms),s(ms,ms)
      dimension plus(ms),ug(ms,ms),u3(ms,ms)
      dimension ihamsign(ms),ihamord(ms),ug2(ms,ms)
      data xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      read(888,*)
      read(888,*)
      do 888 i=1,ns
      read(888,899)(ug(i,j),j=1,ns)
 888  enddo
899   format(16x,<ns>(f10.5))
      read(444,*)(ihamord(i),i=1,ns)
      read(444,*)(ihamsign(i),i=1,ns)
      do 777 i=1,ns
      do 666 j=1,ns
      k=ihamord(j)
      ug2(i,k)=ug(i,j)
 666  enddo
      do 999 j=1,ns
      ug2(i,j)=dble(ihamsign(j))*ug2(i,j)
 999  enddo
 777  enddo
      do 222 i=1,ns
      do 222 j=1,ns
      ug(i,j)=ug2(i,j)
  222 enddo
      end
      subroutine geometryjjls9j(ug)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      character ale(ms),aln(ms),alnp1(ms)
      dimension xl1(ms),xs1(ms),XJ1(ms),XL2(ms),XJ2(ms),imark(ms),
     +xl3(ms),xs3(ms),sgnUG(ms),ugb(ms,ms),vmat(ms,ms),s(ms,ms)
      dimension plus(ms),ug(ms,ms),u3(ms,ms),itarcf(ms)
      data xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      ug=0.0
      write(*,*)'input number of channel, J3'
      read(111,*)ns2,xj3
      if(ns2.ne.ns)stop ' wrong dimension ns'
      write(*,*)'input the LS coupling scheme in LJ 
     1inereasing order'	
      do 11 is=1,ns
      read(111,*)XL1(is),XS1(IS),XJ1(IS),XL2(is),XJ2(IS)
      ale(is) = OrbAngMom(XL2(is))
      aln(is) = TargAngMom(XL1(is))
11    continue
      write(*,*)'input the LS coupling scheme in jj 
     1inereasing order'
      do 12 is=1,ns
      read(111,*)xl3(is),xs3(is),itarcf(is)
      alnp1(is) = OrbAngMom(xl3(is))
      if(is.eq.1)then
       iblck=1
      elseif(itarcf(is).ne.itarcf(is-1))then
       iblck=iblck+1
      endif
      imark(iblck)=is
12    continue
      call achannels(xs1,xs3,aln,ale,alnp1,xj1,xj2)
      read(114,*) xqp1
      read(114,*)(sgnUG(is),is=1,ns)
      read(114,*)(plus(is),is=1,ns)
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
     1.and.xs1(is).eq.xs1(js))then
      a1=xl1(is)
      a2=xs1(is)
      a3=xj1(is)
      b1=xl2(is)
      b2=xs2
      b3=xj2(is)
      c1=xl3(js)
      c2=xs3(js)
      c3=xj3
      c9j=s9j(a1,b1,c1,a2,b2,c2,a3,b3,c3)
      c9j=dsqrt((2*A3+1)*(2*B3+1)*(2*C1+1)*(2*C2+1))*C9J
      else
      c9j=0.d0
      endif
      ug(is,js)=c9j*sgnug(js)
200   continue
      write(129,'(a)')'geometrical Uo Matrix:  <jj|LS> '
      end
      subroutine putumiuda(ie,istg,u,xmiu,dal,dli)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,ms1=500,pai=3.1415926,ryd=13.6058)
      dimension xmiu(me,ms),dal(me,ms,ms),e(me),u(ms,ms),
     +dli(me,ms,ms),ub(ms,ms)
      character aeigchan*20(ms),aionchan*20(ms)
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /channels/aeigchan,aionchan
      common /E2/ee1(ms1)
      write(129,107)(is,is=1,ns)
107   format(16x,<ns>(2x,i3,5x))      
108   format(16x,<ns>(A6,4x))
109   format(16x,<10*ns>A1)
      write(129,108)(aeigchan(js),js=1,ns)
      write(129,109)('-',js=1,10*ns)
      if(istg.ne.0)then
      write(129,110)(xmiu(ie,is),is=1,ns)
      endif
      if(istg.ne.0.and.irad.eq.1)then
      write(129,111)(dal(ie,is,1),is=1,ns)
      write(*,111)(dal(ie,is,1),is=1,ns)
      endif
      if(istg.ne.0)then
      write(129,109)('-',js=1,10*ns)
      endif
110   format('QD',14x,<ns>F10.5)
111   format('Da',14x,<ns>F10.5)
      do is=1,ns
      if(istg.ne.0.and.irad.eq.1)then
      write(129,112)is,aionchan(is),(u(is,js),js=1,ns),dli(ie,is,1)
      else
      write(129,113)is,aionchan(is),(u(is,js),js=1,ns)
      endif
      enddo
112   format(i2,1x,A11,1x,'|',<ns>(f10.5),' | ',f20.10)
113   format(i2,1x,A11,1x,'|',<ns>(f10.5))
      call matcopy(ub,u)
      call minv(ub,ns,detug)
      if(istg.eq.0)then
      write(129,'(a,f10.5)')'|Ug|=',detug
      else
      write(129,'(a,f10.5)')'|U|=',detug
      endif
      end
      subroutine matcopy(utmp,u2)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926)
      dimension u2(ms,ms),utmp(ms,ms)
      common /control/irad,mxe,nchop,ns,ne,mfin
      do is=1,ns
      do js=1,ns
      utmp(is,js)=u2(is,js)
      enddo
      enddo
      end
      subroutine mat3copy2(ie,u,u2)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926)
      dimension u2(ms,ms),u(me,ms,ms)
      common /control/irad,mxe,nchop,ns,ne,mfin
      do is=1,ns
      do js=1,ns
      u2(is,js)=u(ie,is,js)
      enddo
      enddo
      end
      subroutine achannels(xs1,xs3,aln,ale,alnp1,xj1,xj2)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      character ale(ms),aln(ms),alnp1(ms),aeigchan*20(ms),
     +aionchan*20(ms)
      dimension xs1(ms),xs3(ms),is1(ms),is3(ms),xj1(ms),xj2(ms)
      dimension plus(ms),ug(ms,ms),u3(ms,ms)
      data ii/0/,xs2/0.5d0/
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /channels/aeigchan,aionchan
      do js=1,ns
      is1(js)=int(2.0*xs1(js)+1.0)
      is3(js)=int(2.0*xs3(js)+1.0)
      enddo
      do is=1,ns
      write(aeigchan(is),100)is1(is),aLn(is),ale(is),is3(is),alnp1(is)
      aeigchan(is)=adjustl(aeigchan(is))
      write(aionchan(is),101)is1(is),aln(is),xj1(is),ale(is),xj2(is)
      aionchan(is)=adjustl(aionchan(is))
      enddo
100   format(i3,a,a,'_',i1,a)
101   format('(',i1,a,f3.1,')',a,f3.1)
      
      end
      subroutine openfile
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926,ryd=13.6058)
      logical exham
      common /E/EEN,EGRD,xqp1,plus
      common /control/irad,mxe,nchop,ns,ne,mfin,mlast
      common /control2/iconnect,idescend
      common /ham/ham(ms,ms)
      common /exham/exham
      open(134,file='sign.out')
      open(111,file='jjls.in',status='unknown')
      open(777,file='groundlev')
      read(777,*)EGRD
      open(999,file='correla.out')
      write(*,*)'photonic process  irad=1;'
      write(*,*)'collision process irad=0'
      read(*,*) irad
      write(*,*)irad
      if(irad.eq.1)then
      open(125,file='Dalpha-unf-smooth',form='unformatted',status='old')
      open(126,file='Dalpha-unf-eigen-smooth',form='unformatted')
c      open(135,file='DalphaL-smooth.out')
c      open(136,file='DalphaV-smooth.out')
c      open(138,file='Da.out')
      open(137,file='Di.out')
      endif
      open(114,file='order',status='unknown')
      open(123,file='ESUM-unf-smooth',form='unformatted',status='old')
      open(124,file='ESUM-unf-smooth-2',form='unformatted')
      open(127,file='smooth.in',status='old')
      open(128,file='miu-smooth-2.test')
      open(129,file='smooth-ug.out')
      open(133,file='miuang.out')
      open(139,file='descend2.out')
      read(123)mxe,nchop
      rewind(123)
      ns=nchop
      ne=mxe
      write(*,*)'iconnect =0, the last one'
      read(*,*)iconnect
      if(iconnect.eq.0)iconnect=mxe
      write(*,*)'keep the order descended from last energy section  1'
      write(*,*)'reorder according to geometry  0'
      read(*,*)idescend
c      write(*,*)'The number of initial state of PI process?'
c      read(*,*)ni
      inquire(file='HAM-sum.csv',exist=exham)
      if(exham)then
      open(888,file='HAM-sum.csv')
      open(444,file='Ham.order')
      endif
      end

      subroutine 
     +swap2(xmiu,xmiu2,u,upivot,dal,dav,dal2,dav2,mords,msgns,ie)
      implicit real*8 (a-h,o-z)
      parameter (me=500,ms=200,pai=3.1415926)
      dimension mords(ms),msgns(ms)
      dimension xmiu(me,ms),xmiu2(me,ms),u(me,ms,ms),upivot(me,ms,ms),
     +dav(me,ms,ms),dal2(me,ms,ms),dav2(me,ms,ms),dal(me,ms,ms)
      common /control/irad,mxe,nchop,ns,ne,mfin
      common /control2/iconnect,idescend
      do is=1,ns
      do js=1,ns
      upivot(ie,is,js)=u(ie,is,mords(js))
      enddo
      do js=1,ns
      upivot(ie,is,js)=upivot(ie,is,js)*msgns(js)
      enddo
      enddo
      do js=1,ns 
      xmiu2(ie,js)=xmiu(ie,mords(js))
      enddo

      if(irad.eq.1)then
      do m1=1,mfin
      do js=1,ns      
      dal2(ie,js,m1)=dal(ie,mords(js),m1)
      dav2(ie,js,m1)=dav(ie,mords(js),m1)
      enddo
      do 7 is=1,ns
      dal2(ie,is,m1)=msgns(is)*dal2(ie,is,m1)
      dav2(ie,is,m1)=msgns(is)*dav2(ie,is,m1)
7     enddo
      enddo
      endif
      end
       
      end module assign_eigenchann
