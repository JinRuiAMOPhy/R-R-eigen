      program MAKERWFP
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


