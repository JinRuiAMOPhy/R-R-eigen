      module correlate_pivot
      use constants
      use structs
      use init_final
      implicit none 
      integer :: ie, i, is, ialfa, js
! Caution namespace of these will be global !!
      type(control) :: CTRL
      type(SDmat) :: SD
      type(SDmat) :: SD_CP
      type(SDmat) :: SD_inher
      type(fname) :: FN
      real(long), dimension(ms,ms) :: s, spivot, Ustand, U2Ug
      real(long), dimension(me) :: E
      integer, dimension(me,ms) :: iord, imsgn
      integer, dimension(ms) :: mord, msgn
! global data
      contains 

      subroutine step1(DIR)
      implicit none
      real(long), dimension(ms,ms) :: ut
      character(len = longlong_str), intent(in) :: DIR
      real(long) :: det
      integer :: i, j

      call initiate(DIR, FN, CTRL, SD, SD_inher)
      loop_E : do ie=1,SD%mxe
        write(FN%unit_dbg,'(A,i4)')'# Energy point',ie
        call readumiuda(ie,e,SD)
        call contrljjorder(SD%U)
        if(ie == 1)then
          if(CTRL%yes_inherit)then
            call xmatcor(SD_inher%U,SD%U,s,SD%nchop)
            write(FN%unit_dbg,"(A)")
     &      "# Pivot correlation with inherit matrix(partially)"
            call pivot(s,spivot,mord,msgn,SD%nchop,ms)
          else if (index(CTRL%which_corr_u_mat, 'ug') /= 0) then
            write(FN%unit_dbg,"(A)")"# Pivot correlation with UG"
            call xmatcor(SD%Ug, SD%U, U2Ug, SD%nchop)
            call pivot(U2Ug,spivot,mord,msgn,SD%nchop,ms)
            write(FN%unit_dbg,"(A)")"# U2Ug"
            do i = 1, SD%nchop
             call print_vector_comment(U2Ug(i,:), SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
            end do 
          else if (index(CTRL%which_corr_u_mat, 'first') /= 0) then
            write(FN%unit_dbg,"(A)")"# Pivot the first point"
            call pivot(SD%U,SD_CP%U,mord,msgn,SD%nchop,ms)
          else if (index(CTRL%which_corr_u_mat, 'prev') /= 0) then
            write(FN%unit_dbg,"(A)")"# Pivot the first point"
            call pivot(SD%U,SD_CP%U,mord,msgn,SD%nchop,ms)
          endif !yes_inherit
          if(CTRL%yes_print == 1)then
            write(FN%unit_dbg,"(A)")"# Correlation matrix"
            do i = 1, SD%nchop
             call print_vector_comment(s(i,:), SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
            end do 
          end if 
          iord(ie,:)=mord(:)
          imsgn(ie,:)=msgn(:)
          call swapUmiu(SD%U,SD_CP%U,SD%miu,SD_CP%miu,msgn,mord)
          if(CTRL%yes_print) then
            write(FN%unit_dbg,"(A)")"# after swapUmiu"
            do i = 1, SD%nchop
             call print_vector_comment(SD_CP%U(i,:), SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
            end do 
            write(FN%unit_dbg,"(A)")"# swap order (alpha)"
            call print_int_vector_comment(mord, SD%nchop, 10, 
     &       "(i10)", FN%unit_dbg)
             call print_int_vector_comment(msgn, SD%nchop, 10, 
     &       "(i10)", FN%unit_dbg)
          end if !yes_print
        else if (ie > 1 ) then
          call xmatcor(Ustand,SD%U,s,SD%nchop)
          if(CTRL%yes_print) then
              write(FN%unit_dbg,"(A)")"# Correlation with Standard U"
              do i = 1, SD%nchop
               call print_vector_comment(s(i,:), SD%nchop, 10, 
     &         "(f10.4)", FN%unit_dbg)
              end do 
          end if 
          call pivot(s,spivot,mord,msgn,SD%nchop,ms)
          call swapUmiu(SD%U,SD_CP%U,SD%miu,SD_CP%miu,msgn,mord)
          iord(ie,:)=mord(:)
          imsgn(ie,:)=msgn(:)
        endif
        if (CTRL%yes_skip .and. SD%Excit > CTRL%skip_beg .and. 
     &     SD%Excit < CTRL%skip_end) cycle
        if(CTRL%yes_rad) 
     &    call swapda(SD%DL,SD%DV,SD_CP%DL,SD_CP%DV,mord,msgn)
        write(FN%unit_dbg,'(A)')
     &   '# => order and sign'
        call print_int_vector_comment(mord, SD%nchop, 10, 
     &       "(i10)", FN%unit_dbg)
        call print_int_vector_comment(msgn, SD%nchop, 10, 
     &       "(i10)", FN%unit_dbg)
        if(CTRL%yes_man_swp) then
          write(FN%unit_dbg, "(A)", advance = 'no')
     &     '# Manually swapped after '
          write(FN%unit_dbg, "(A)", advance = 'yes')
     &     ' sort the correlation matrix'
          call print_int_vector_comment(CTRL%man_swp(ie,:), 
     &    SD%nchop, 10, "(i10)", FN%unit_dbg)
          call print_int_vector_comment(CTRL%man_sgn(ie,:), 
     &    SD%nchop, 10, "(i10)", FN%unit_dbg)
          call swapUmiu(SD_CP%U,SD_CP%U,SD_CP%miu,SD_CP%miu,
     &      CTRL%man_sgn(ie,:),CTRL%man_swp(ie,:))
          if(CTRL%yes_rad)
     &      call swapda(SD_CP%DL,SD_CP%DV,SD_CP%DL,
     &SD_CP%DV,CTRL%man_swp(ie,:),CTRL%man_sgn(ie,:))
        end if
        if(ie == 1) call preput
        call put(ie,SD_CP%U,SD_CP%DL,SD_CP%DV,SD%Di, 
     &     SD%Excit,SD_CP%miu,msgn,mord)
        call genmiuang(SD_CP%U,SD_CP%miu)
        ut(:,:) = SD_CP%U(:,:) 
        call minv(ut, SD%nchop, det)
        write(FN%unit_dbg, "('# |U| = ', f10.5)") det
        call xmatcor(SD_CP%U, SD%Ug, U2Ug, SD%nchop)
        write(FN%unit_dbg, "(A)") '# Correlation with UG'
        do i = 1, SD%nchop
          call print_vector_comment(U2Ug(i,:), SD%nchop, 10,
     &         "(f10.4)", FN%unit_dbg)
        end do 
        if (index(CTRL%which_corr_u_mat, 'ug') /= 0) then
          Ustand(:,:) = SD%Ug(:,:)
        else if(index(CTRL%which_corr_u_mat, 'prev') /= 0) then
          Ustand(:,:) = SD_CP%U(:,:)
        else if(ie == 1) then
          Ustand(:,:) = SD_CP%U(:,:)
        end if

        if(ie == CTRL%herit_p) call leave_herit(SD_CP%U)
      end do loop_E
        write(FN%unit_dbg,"('# ',A)") 'Eigenchannel labels'
        write(FN%unit_dbg,"('# ')",advance = 'no') 
        do i = 1, SD%nchop
          write(FN%unit_dbg,"('(',A,')',A,1x)",advance = 'no') 
     &      trim(SD%chan_ion_str(i)),trim(SD%chan_eig_str(i))
        end do 
        write(FN%unit_dbg,"('# ',A)") ''
        write(FN%unit_dbg,"('# ',A)") 'Eigenchannel labels (LaTex)'
        write(FN%unit_dbg,"('# ')",advance = 'no') 
        do i = 1, SD%nchop
          write(FN%unit_dbg,"('(',A,')',A,1x)",advance = 'no') 
     &      trim(ionchan_str_Latex(SD%chan_ion_str(i))),
     &      trim(LS_str_Latex(SD%chan_eig_str(i)))
        end do 
      end subroutine step1

      subroutine leave_herit(Ustand)
      use constants
      use structs
      use constants
      use structs
      implicit none
      integer :: isort, jn
      real(long), dimension(:,:) :: Ustand

      write(FN%unit_herit,*) SD%nchop
      do isort=1,SD%nchop
        write(FN%unit_herit,'(50f10.5)')(Ustand(isort,jn),jn=1,SD%nchop)
      enddo
      end subroutine leave_herit

      subroutine swapUmiu(u,U_after,xmiu,xmiu2,msgn,mord)
      use constants
      use structs
      implicit none
      integer :: is, js
      integer, dimension(:) :: mord, msgn
      real(long), dimension(:) :: xmiu, xmiu2
      real(long), dimension(:,:) :: u, u_after
        do is=1,SD%nchop
        do js=1,SD%nchop
          U_after(is,js)=u(is,mord(js))
        enddo
        do js=1,SD%nchop
          U_after(is,js)=U_after(is,js)*msgn(js)
        enddo
        enddo
        do js=1,SD%nchop
          xmiu2(js)=xmiu(mord(js))
        enddo
      end subroutine swapUmiu

      subroutine swapda(dal,dav,dal2,dav2,iord,imsgn)
      use constants
      use structs
      implicit none 
      integer :: m1, is
      integer, dimension(:) :: iord, imsgn
      real(long), dimension(:,:) :: dal, dav, dal2, dav2
      if(.not.CTRL%yes_rad) return 
      do m1=1,SD%mfin
        do is=1,SD%nchop
        dal2(is,m1)=dal(iord(is),m1)
        dav2(is,m1)=dav(iord(is),m1)
        enddo        
        do is=1,SD%nchop
        dal2(is,m1)=dal2(is,m1)*imsgn(is)
        dav2(is,m1)=dav2(is,m1)*imsgn(is)
        enddo       
      enddo 
      return 
      end subroutine swapda


      subroutine genmiuang(u2,xmiu2)
      use constants
      use structs
      implicit none 
      integer, parameter :: iq = 0 ! matrix to angle 
      integer :: is, i, j
      real(long), dimension(:), intent(in) :: XMIU2
      real(long), dimension(:,:), intent(in) :: U2
      real(long), dimension(ms,ms) :: A, R, DS, DC

      CALL UANGCON(U2,A,R,SD%nchop,0,DS,DC,ms,iq)
      WRITE(FN%unit_S_out_for,'(I3,4000(1PE16.7))')SD%nchop,SD%Excit,
     &(xmiu2(is)+CTRL%xplus(is),is=1,SD%nchop), 
     &((A(I,J),J=I+1,SD%nchop),I=1,SD%nchop-1)
      end subroutine genmiuang


      subroutine preput
      use constants
      use structs
      implicit none 
      integer :: m1, len, ifile, jfile
      character(len=3) ::A1
      character(len = longlong_str) :: fnameV, fnameL
      
      do m1=1,SD%mlast
        write(A1,'(I3)')m1
        A1=adjustl(A1)
        len=len_trim(a1)
        write(fnameV,"(A,'DalphaV-smooth.',A,'.out')")
     &   trim(adjustl(CTRL%path)),trim(adjustl(a1))
        write(fnameL,"(A,'DalphaL-smooth.',A,'.out')")
     &   trim(adjustl(CTRL%path)),trim(adjustl(a1))
        ifile=1000+m1
        jfile=1500+m1
        open(ifile,file=fnameL)
        !write(ifile,'(2a20)')'E','DalphaL'
        open(jfile,file=fnameV)
        !write(jfile,'(2a20)')'E','DalphaV'
      enddo
      end
      subroutine put(ie,u2,dal2,dav2,dli,Excit,xmiu2,msgn,mord)
      use constants
      use structs
      implicit none
      integer, intent(in) :: ie
      integer :: m1, ifile, jfile, ialfa, js, is
      real(long) :: Excit
      real(long), dimension(:,:), intent(in) :: dli,dal2,
     &dav2, u2
      real(long), dimension(:), intent(in) :: xmiu2
      integer, dimension(ms) :: mord, msgn
      if(CTRL%yes_rad)then
        do m1=1,SD%mfin
          ifile=1000+m1
          jfile=1500+m1
          write(FN%unit_D_out_unf)(dal2(ialfa,m1),ialfa=1,SD%nchop)
          write(FN%unit_D_out_unf)(dav2(ialfa,m1),ialfa=1,SD%nchop)
          write(FN%unit_D_out_unf)(dli(ialfa,m1),ialfa=1,SD%nchop)
          write(ifile,'(100E20.10)')
     &    Excit,(dal2(ialfa,m1),ialfa=1,SD%nchop)
          write(jfile,'(100E20.10)')
     &    Excit,(dav2(ialfa,m1),ialfa=1,SD%nchop)
        enddo
      endif
      write(FN%unit_S_out_unf)(SD_CP%miu(js)+CTRL%xplus(js),
     + js=1,SD%nchop)
      do is=1,SD%nchop
        write(FN%unit_S_out_unf)(u2(is,js),js=1,SD%nchop)
      enddo

      write(FN%unit_dbg,'(A)') '#'
      write(FN%unit_dbg,'(12x,2x,1000(4x,a5,1x))') 
     &  (SD%chan_eig_str(i), i = 1, SD%nchop)
      write(FN%unit_dbg,'("  miu",7x,50f10.5)')
     & (SD_CP%miu(js)+CTRL%xplus(js),
     & js=1,SD%nchop) 
      if(CTRL%yes_rad)then
      do m1 = 1, SD%mlast
       write(FN%unit_dbg,'("  Da",8x,50f10.5)')(SD_CP%DL(js,m1),
     & js=1,SD%nchop)
      end do 
      endif
      write(FN%unit_dbg,'(1000a)') ('-',i=1, SD%nchop * 10 + 12) 
      do is=1,SD%nchop
      write(FN%unit_dbg, "(A10,1x,'|')", advance = 'no')
     &  SD%chan_ion_str(is)
      call print_vector(SD_CP%U(is,:), SD%nchop, 10, 
     &       "(f10.5)", FN%unit_dbg)
      enddo
      end subroutine put

      subroutine readumiuda(ie,e,SD)
      use constants
      use structs
      implicit none 
      type(SDmat), intent(inout) :: SD
      integer, intent(in) :: ie
      integer :: mxe, nchop
      integer :: is, js, m1, ialfa, mt
      real(long), dimension(me), intent(out) :: e
      real(long) :: etot, wdet, etr
      read(FN%unit_S_in)mxe,nchop
      read(FN%unit_S_in)e(ie)
      SD%Excit=e(ie)
      read(FN%unit_S_in)(SD%miu(is),is=1,SD%nchop)
      do is=1,SD%nchop
        read(FN%unit_S_in)(SD%U(is,js),js=1,SD%nchop)
      enddo
      write(FN%unit_S_out_unf)mxe,nchop
      write(FN%unit_S_out_unf)SD%Excit
      if(CTRL%yes_print) then
        write(FN%unit_dbg,'("#",200A)')('*',is=1,20)
        write(FN%unit_dbg,"('# Read in at E = ', f10.5,"//
     &   " 'Ryd./Z^2(IP1=0)')")
     &   SD%Excit
        write(FN%unit_dbg,"(A)")'# miu'
        call print_vector_comment(SD%miu, SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
        write(FN%unit_dbg,"(A)")'# U'
        do is=1,SD%nchop
          call print_vector_comment(SD%U(is,:), SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
        enddo
      end if

      if(.not.CTRL%yes_rad) return 

      if(ie == 1)then
        read(FN%unit_D_in)SD%mlast ! # of <i| in <I|D|f> data blocks
        write(FN%unit_D_out_unf)SD%mlast
      endif
      ! the mfin's block of <mfin|D|f>
      read(FN%unit_D_in)etot,wdet,SD%mfin,(SD%Eistat(m1),m1=1,SD%mfin)
      write(FN%unit_D_out_unf)etot,wdet,SD%mfin,
     & (SD%Eistat(m1),m1=1,SD%mfin)

      do m1=1,SD%mfin
        read(FN%unit_D_in)etr,mt, (SD%DL(ialfa,m1),ialfa=1,SD%nchop)
        read(FN%unit_D_in)(SD%DV(ialfa,m1),ialfa=1,SD%nchop)
        read(FN%unit_D_in)(SD%Di(ialfa,m1),ialfa=1,SD%nchop)
        if(CTRL%yes_print) then
          write(FN%unit_dbg,"('# <b',i3,'|D(L)|f> readin')")m1
          call print_vector_comment(SD%DL(:,m1), SD%nchop, 10, 
     &       "(f10.4)", FN%unit_dbg)
        end if
      enddo
      end subroutine readumiuda

      subroutine contrljjorder(u)
      use constants
      use structs
      implicit none
      integer :: is, js
      real(long), dimension(:,:), intent(inout) :: U
      real(long), dimension(ms,ms) :: utmp
      if(index(CTRL%JJordTyp,'rev') /= 0)then
        write(*,*) ' Cauiton JJ order reversed for the '
        write(*,*)'conneting process!'
        do is=1,SD%nchop
          do js=1,SD%nchop
            utmp(is,js)=u(SD%nchop-is+1,js)
          enddo
        enddo
      else if(index(CTRL%JJordTyp,'readin') /= 0)then
        do is=1,SD%nchop
          do js=1,SD%nchop
            utmp(CTRL%jjord(is),js)=u(is,js)
          enddo
        enddo
      else if(index(CTRL%JJordTyp,'no') /= 0)then
        return
      endif

      if(CTRL%yes_print) then
        write(FN%unit_dbg,"('# ',A)")'U_ialfa jj reordered'
        do is=1,SD%nchop
         call print_vector_comment(utmp(is,:), SD%nchop, 10, 
     &   "(f10.4)", FN%unit_dbg)
        enddo
      end if 
      u(:,:) = utmp(:,:)
      return
      end subroutine contrljjorder

      end module correlate_pivot
