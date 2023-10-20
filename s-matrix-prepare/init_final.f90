module init_final
   use constants
   use structs
   use darray
   use stdio
   use numeric
   implicit none
contains
   subroutine from_str_get_dbg(CTRL)
      implicit none
      type(control), intent(inout) :: CTRL
      if(index(CTRL%dbgstr, 'all') /= 0 .or. index(CTRL%dbgstr, 'yes') /= 0) then
         CTRL%yes_print= .true.
      end if
      return
   end subroutine from_str_get_dbg

   subroutine read_control(FN, CTRL, SD, SD_inher, keys)
      implicit none
      type(fname), intent(inout) :: FN
      type(control), intent(inout) :: CTRL
      type(SDmat), intent(inout) :: SD
      type(SDmat), intent(inout) :: SD_inher
      integer, parameter :: funit = 15
      integer :: i, ncol, j, ntmp, ie, it, it2, iseg, nseg, ibegin, iend
      logical :: file_ex
      character(len = short_str), dimension(:), allocatable, intent(out) :: keys
      character(len = longlong_str), dimension(:), allocatable :: vals
      character(len = longlong_str) :: skip_str, JJord_str, plusOne_str, &
         ug_sign_str, chan_ion_target_cf_str, chan_ion_str, chan_eig_str, &
         manual_swap, manual_sgn
      character(len = short_str), dimension(ms) :: tmp_chan_eig_str
      character(len = long_str), dimension(me) :: man_swp, man_sgn
      real(long), dimension(ms) :: UgSign
      real(long), dimension(ms,ms) :: Ut
      real(long) :: det
      integer, dimension(ms) :: chan_targcf, mord,msgn
      logical :: yes_ug_swap

      CTRL%herit_p = -1
      CTRL%JJordTyp = ''
      plusOne_str = ''
      JJord_str = ''
      skip_str = ''
      chan_ion_str = ''
      chan_eig_str = ''
      ug_sign_str = ''
      chan_ion_target_cf_str = ''
      yes_ug_swap = .False.
      CTRL%which_corr_u_mat = 'first'
      manual_swap = ''
      manual_sgn = ''
      FN%S_in = trim(CTRL%path)//trim(FN%S_in)
      FN%dbg = trim(CTRL%path)//trim(FN%dbg)
      SD%JT = -1
      inquire(file = FN%S_in, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%S_in)
         stop
      end if
      
      open(unit =FN%unit_dbg, file = FN%dbg, action = 'write', &
           form = 'formatted')
      open(unit = FN%unit_S_in, file = FN%S_in, action = 'read', &
           form = 'unformatted')
      read(FN%unit_S_in) SD%mxe, SD%nchop
      CTRL%mxe = SD%mxe
      CTRL%nchop = SD%nchop
      rewind(FN%unit_S_in)
      write(FN%unit_dbg,"('# Energy points', i3, ' openchan ',i3)") SD%mxe, SD%nchop

      FN%smooth_ctrl= trim(CTRL%path)//trim(FN%smooth_ctrl)
      inquire(file = FN%smooth_ctrl, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%smooth_ctrl)
         stop
      end if
      open(FN%unit_smooth_ctrl, file = FN%smooth_ctrl, action = 'read', &
           form = 'formatted')
      call read_key_val_array(FN%unit_smooth_ctrl, keys, vals)

      do i = 1, size(keys)
         write(FN%unit_dbg,"('# ', A, 2x, A)")trim(keys(i)), trim(vals(i))
         if(keys(i) == 'jjorder') then
            JJord_str = vals(i)
         else if(keys(i) == 'inherit') then
            if(index(vals(i), 'yes') == 0) then 
              CTRL%yes_inherit = .True.
              FN%inherit = 'inherit.txt'
            else if (index(vals(i), 'no') == 0) then
              CTRL%yes_inherit = .False.
            else 
              CTRL%yes_inherit = .True.
              FN%inherit = vals(i)
            end if 
         else if(keys(i) == 'plus_princ') then
            plusOne_str = vals(i)
         else if(keys(i) == 'chan_ion') then
            chan_ion_str = vals(i)
            call get_string_array_space(chan_ion_str, SD%chan_ion_str)
            if(size(SD%chan_ion_str) < SD%nchop) then
               stop '# too few ion chan are provided'
            else if(size(SD%chan_ion_str) > SD%nchop) then
               write(FN%unit_dbg,"('# ', A, 2x, i3,1x,A)") &
                 '# too many ion chan are provided, only the first',SD%nchop, &
                 'ones will be used'
            end if
         else if(keys(i) == 'chan_eig') then
            chan_eig_str = vals(i)
         else if(keys(i) == 'ug_sign') then
            ug_sign_str = vals(i)
         else if(keys(i) == 'chan_ion_target_cf') then
            chan_ion_target_cf_str = vals(i)
         else if(keys(i) == 'manual_swap') then
            manual_swap = vals(i)
         else if(keys(i) == 'manual_sign') then
            manual_sgn = vals(i)
         else if(keys(i) == 'herit_point') then
              read(vals(i), *) CTRL%herit_p
         else if(keys(i) == 'correlate_point') then
            read(vals(i), *) CTRL%which_corr_u_mat
         else if(keys(i) == 'skip') then
            if(index(vals(i),'no') ==0) then
               skip_str = vals(i)
               CTRL%yes_skip = .True.
            else if(index(vals(i),'no') /=0) then
               CTRL%yes_skip = .False.
            end if 
         else if(keys(i) == 'radiation') then
            if(index(vals(i),'yes') /=0) CTRL%yes_rad= .True.
            if(index(vals(i),'no') /=0) CTRL%yes_rad = .False.
         else if(keys(i) == 'debug') then
            read(vals(i), *) CTRL%dbgstr
            call from_str_get_dbg(CTRL)
         else if(keys(i) == 'jt' .or. keys(i) == 'j_tot') then
            read(vals(i), *) SD%JT
         end if
      end do 
      if(SD%JT == -1 .and. chan_eig_str == '') &
        stop 'neither JT nor eigen-channels are provided'
      if(index(JJord_str, 'rev') /= 0) then
        CTRL%JJordTyp = 'rev'
      else if(index(JJord_str, 'no') /= 0) then
        CTRL%JJordTyp = 'no'
      else
        CTRL%JJordTyp = 'readin'
        !print *, JJord_str
        ncol = from_line_get_ncol(JJord_str)
        if(ncol /= SD%nchop ) &
           stop 'wrong dimension of jjord'
        !allocate(T%IPs(ncol))
        !print*, ncol
        read(vals(i), *)(CTRL%jjord(j), j = 1, ncol)
         !print*, CTRL%jjord(1:ncol)
      end if
      if(skip_str /= ''.or. index(skip_str, 'no') /=0 ) then
        ncol = from_line_get_ncol(skip_str)
        read(skip_str, "(2f)")CTRL%skip_beg, CTRL%skip_end
      end if
      if(chan_ion_str == '') stop 'jjls notations are missing'
      !print *, 'chan_ion_str:', trim(chan_ion_str)
      do j = 1, SD%nchop
        SD%jjls(j, :) = read_jjls(SD%chan_ion_str(j))
        !print "(2A,6f5.2)", trim(SD%chan_ion_str(j)), '->', SD%jjls(j, :)
      end do 

      if(chan_ion_target_cf_str /= '') then
        read(chan_ion_target_cf_str, *) (SD%chan_targcf(i), i = 1, SD%nchop)
      else 
        SD%chan_targcf(:) = (/(1, i = 1, SD%nchop) /)
      end if 

      if(chan_eig_str /= '') then
        call get_string_array_space(chan_eig_str, SD%chan_eig_str)
        do j = 1, SD%nchop
          SD%LSJ(j, :) = read_LSJ(SD%chan_eig_str(j))
        end do 
      else
        call gen_LSJ_chan(SD)
      end if

      if(chan_eig_str == '') then
        call gen_eig_str(SD)
      end if

      if(ug_sign_str /= '') then
        read(ug_sign_str, *) (UgSign(i), i = 1, SD%nchop)
      else 
        UgSign(:) = (/(1, i = 1, SD%nchop) /)
      end if

      SD%UG(:,:) = geometryjjls9j(SD%nchop, SD%jjls, SD%LSJ, UgSign, SD%chan_targcf)
      call pivot(SD%UG,SD%UG,mord,msgn,SD%nchop,ms)

      write(FN%unit_dbg,'(12x,2x,1000(4x,a5,1x))') &
     &  (SD%chan_eig_str(i), i = 1, SD%nchop)
      write(FN%unit_dbg,'(1000a)') ('-',i=1, SD%nchop * 10 + 12) 
      do i=1,SD%nchop
        write(FN%unit_dbg, "(A10,1x,'|')", advance = 'no') &
     &  SD%chan_ion_str(i)
        call print_vector(SD%Ug(i,:), SD%nchop, 10, &
     &       "(f10.5)", FN%unit_dbg)
      enddo
      ut(:,:) = SD%Ug(:,:)
      call minv(ut, SD%nchop, det)
      write(FN%unit_dbg, "('# |Ug| = ', f10.5)") det

      do i = 1, SD%nchop
        if(mord(i) /= i) yes_ug_swap = .True.
      end do 
      if(yes_ug_swap) then
        write(FN%unit_dbg, "(A)") "#--- Ug was swapped automatically ---"
        call print_int_vector_comment(mord,  SD%nchop, 10, &
     &  "(i4)", FN%unit_dbg)
        call print_int_vector_comment(msgn,  SD%nchop, 10, &
     &  "(i4)", FN%unit_dbg)
        call swap_str_array(SD%chan_eig_str, tmp_chan_eig_str, mord, SD%nchop)
        SD%chan_eig_str(1:SD%nchop) = tmp_chan_eig_str(1:SD%nchop)
      end if
      if(index(plusOne_str, 'yes') /= 0) then
        read(plusOne_str,*)(CTRL%xplus(i),i=1,SD%nchop)
        !print*, CTRL%xplus(1:SD%nchop)
      else
        do i =1,SD%nchop
          CTRL%xplus(i)=0.d0
        enddo
      endif
      if(CTRL%herit_p == -1) then
        CTRL%herit_p = SD%mxe
      end if 

      call get_string_array_sep(manual_swap, man_swp, ';', nseg)
      if(nseg > 0 ) then
        CTRL%yes_man_swp = .True.
        write(FN%unit_dbg, "('# manual swap order in step2')")
      else 
        CTRL%yes_man_swp = .False.
      end if 
      do iseg = 1, nseg
        it = index(man_swp(iseg), '-')
        read(man_swp(iseg)(1:it-1), *) ibegin
        it2 = index(man_swp(iseg), ':')
        read(man_swp(iseg)(it+1:it2-1), *) iend
        write(FN%unit_dbg, "('# ',A)") trim(man_swp(iseg)) 
        do ie = 1, SD%mxe
          if (ie >= ibegin .and. ie <= iend) then
            read(man_swp(iseg)(it2+1:long_str), *) (CTRL%man_swp(ie, i), i = 1, SD%nchop)
          else
            do i = 1, SD%nchop
              CTRL%man_swp(ie, i) = i
            end do
          end if
        end do 
      end do 

      call get_string_array_sep(manual_sgn, man_sgn, ';', nseg)
      if( nseg > 0 ) &
      &  write(FN%unit_dbg, "('# manual swap order in step2')")
      do iseg = 1, nseg
        it = index(man_sgn(iseg), '-')
        read(man_sgn(iseg)(1:it-1), *) ibegin
        it2 = index(man_sgn(iseg), ':')
        read(man_sgn(iseg)(it+1:it2-1), *) iend
        write(FN%unit_dbg, "('# ',A)") trim(man_sgn(iseg)) 
        do ie = 1, SD%mxe
          if (ie >= ibegin .and. ie <= iend) then
            read(man_sgn(iseg)(it2+1:long_str), *) (CTRL%man_sgn(ie, i), i = 1, SD%nchop)
          else
            do i = 1, SD%nchop
              CTRL%man_sgn(ie, i) = 1
            end do
          end if
        end do 
      end do 
      deallocate(vals)

      close(FN%unit_smooth_ctrl)

      if(CTRL%yes_rad) then
        FN%D_in = trim(CTRL%path)//trim(FN%D_in)
        open(FN%unit_D_in,file=FN%D_in,form='unformatted', &
       &action='read')
      end if
      FN%S_out_unf = trim(CTRL%path)//trim(FN%S_out_unf)
      open(FN%unit_S_out_unf,file=FN%S_out_unf,form='unformatted', &
     &action='write')
      FN%S_out_for= trim(CTRL%path)//trim(FN%S_out_for)
      open(FN%unit_S_out_for,file=FN%S_out_for,form='formatted', &
     &action='write')
      FN%D_out_unf = trim(CTRL%path)//trim(FN%D_out_unf)
      open(FN%unit_D_out_unf,file=FN%D_out_unf,form='unformatted', &
     &action='write')
      FN%D_out_for = trim(CTRL%path)//trim(FN%D_out_for)
      open(FN%unit_D_out_for,file=FN%D_out_for,form='formatted', &
     &action='write')
      FN%herit = trim(CTRL%path)//trim(FN%herit)
      open(FN%unit_herit,file=FN%herit,form='formatted', &
     &action='write')
      FN%U_sign= trim(CTRL%path)//trim(FN%U_sign)
      open(FN%unit_u_sign,file=FN%U_sign, form='formatted', &
     &action='write')
      
      if(CTRL%yes_inherit) call inherit(FN, CTRL, SD_inher)
      print*,"smooth.inp done"

      return
   !   call comment_str_among_char('INPUT ERROR', '*', stderr)
   end subroutine read_control 
   subroutine inherit(FN, CTRL, SD_inher)
     implicit none
     type(fname), intent(inout) :: FN
     type(control), intent(inout) :: CTRL
     type(SDmat), intent(out) :: SD_inher
     integer :: i, j
     inquire(file=FN%inherit,exist=CTRL%yes_inherit)
     SD_inher%U = ZERO
     if(.not.CTRL%yes_inherit)then
       write(*,*)'not inherit'
       return
     else
       FN%inherit = trim(CTRL%path)//trim(FN%inherit)
       open(unit = FN%unit_inherit,file=FN%inherit,form='formatted', &
         &action='read')
       read(FN%unit_inherit,*)SD_inher%nchop
       write(FN%unit_dbg,*)'inheritage matrix'
       write(FN%unit_dbg,*)SD_inher%nchop
       if(SD_inher%nchop >= CTRL%nchop) stop 'wrong inherit U dimension, check plz'
       do i=1,SD_inher%nchop
       read(FN%unit_inherit,'(50f10.5)')(SD_inher%U(i,j),j=1,SD_inher%nchop)
       write(FN%unit_dbg,'(50f10.5)')(SD_inher%U(i,j),j=1,SD_inher%nchop)
       enddo
       do i=SD_inher%nchop+1,CTRL%nchop
       SD_inher%U(i,i)=1.0
       enddo
       close(FN%unit_inherit)
     endif
   end subroutine inherit

   subroutine initiate(DIR, FN, CTRL, SD, SD_inher)
      implicit none
      type(fname), intent(inout) :: FN
      character(len = longlong_str), intent(in) :: DIR
      type(control), intent(inout) :: CTRL
      type(SDmat), intent(inout) :: SD
      type(SDmat), intent(inout) :: SD_inher
      character(len = short_str), dimension(:), allocatable :: keys
      character(len = short_str) :: sfmtstr 
      integer :: i, funit_dbg
      CTRL%path = DIR
      CTRL%dbgstr = ''
      write(FN%unit_dbg, "('#----- Parameters -------- ')") 
      call read_control(FN, CTRL, SD, SD_inher, keys)
      if(.not.check_item_in_list(keys,'radiation')) then
         write(FN%unit_dbg, "('radiation')", advance = 'no')
         write(FN%unit_dbg, "(' (internal) = ')", advance = 'no')
         write(FN%unit_dbg, "(A)") CTRL%dbgstr
      else
         write(FN%unit_dbg, *)
      end if
      !write(funit_dbg, "('# IP(in cm-1):')", advance = 'yes')
      !sfmtstr = '(E20.10)'
      !call print_vector_comment(T%IPs, T%n_IP, print_col_l_fmt, &
      !     sfmtstr, funit_dbg)
      !write(funit_dbg, "('# IP index for each chann:')", advance = 'yes')
      !sfmtstr = '(i4)'
      !call print_int_vector_comment(T%IP_seq, T%nchan, print_col_s_fmt, &
      !     sfmtstr, funit_dbg)
      write(FN%unit_dbg, "('#----- Parameters -------- ')") 
      deallocate(keys)

      return 
   end subroutine initiate

   subroutine gen_eig_str(SD) 
      implicit none
      type(SDmat), intent(inout) :: SD
      integer :: i
      real(long) :: S, L
      character(len=4) :: AS
      character(len=1) :: AL
      do i = 1, SD%nchop
        S = SD%LSJ(i,1); L = SD%LSJ(i,2)
        write(AS, "(i4)") int(TWO*S+ONE)
        AL = TargAngMom(L)
        write(SD%chan_eig_str(i) , "(A,A)") trim(adjustl(AS)), trim(adjustl(AL))
      end do 
      return
   end subroutine gen_eig_str

   subroutine gen_LSJ_chan(SD) 
      implicit none
      type(SDmat), intent(inout) :: SD
      integer :: i, n
      real(long), dimension(ms,6) :: LS
      real(long) :: s2, l2, t2, St, Lt, Jt, sc, lc, jc
      LS = ZERO
      sc = HALF
      n = 0
      print *, 'SD%JT', SD%JT
      do i = 1, SD%nchop
        St = SD%jjls(i, 1); Lt = SD%jjls(i, 2); Jt = SD%jjls(i, 3)
        lc = SD%jjls(i, 5); jc = SD%jjls(i, 6)
        t2 = SD%chan_targcf(i)
        !print '(A,8f5.1)', trim(SD%chan_ion_str(i)), st, lt, jt, lc,jc 
        !print '(A,8f5.1)','s1+s2,s1-s2',St + sc,abs(St - sc)
        !print '(A,8f5.1)','L1+L2,L1-L2',Lt + lc,abs(Lt - lc)
        do s2 = abs(St - sc), St + sc, ONE
          do l2 = abs(Lt - lc), Lt + lc, ONE
            if((SD%JT > l2 + s2) .or. (SD%JT < abs(l2 - s2))) cycle
            if((SD%JT > Jt + jc) .or. (SD%JT < abs(Jt - jc))) cycle
            if(.not. this_term_exist_in_2D_array(LS, St,Lt,lc, S2, L2, t2, n)) then
              n = n + 1
              LS(n,1) = t2; LS(n,2) = St
              LS(n,3) = Lt; LS(n,4) = lc; LS(n, 5) = S2; LS(n,6) = L2
         !     print '(A, 8f5.1)', 'y', s2, l2
            else
         !     print '(A, 8f5.1)', 'n', s2, l2
            end if
          end do 
        end do
      end do
      do i = 1, n
        print "(i3,1x,'->',8f5.1)", i , LS(i,5:6)
        SD%LSJ(i,1) = LS(i,5); SD%LSJ(i,2) = LS(i,6)
        SD%LSJ(i,3) = SD%JT
      end do 
      return
   end subroutine gen_LSJ_chan
   function this_term_exist_in_2D_array(LS,St, Lt, lc, s2,l2,t2, n) result(yes) 
     implicit none
     integer :: i
     integer, intent(in) :: n
     real(long), dimension(ms,6), intent(in) :: LS
     real(long), intent(in) ::St, Lt, lc, s2, l2, t2
     logical :: yes
     yes = .False.
     if(n == 0) then
       yes = .False.
       return
     end if
     do i = 1, n
       if((LS(i,1) == t2) .and. (LS(i,2) == St) .and. (LS(i,3) == Lt) .and. &
          (LS(i,4) == lc) .and. (LS(i,5) == s2) .and. (LS(i,6) == l2)) then
         yes = .True.
         return
       else
         yes = .False.
       end if
     end do 
     return
   end function this_term_exist_in_2D_array

   subroutine finalize(FN)
      type(fname), intent(inout) :: FN
      !if(allocated(S%d_mu)) deallocate(S%d_mu)
      !close(funit_dbg)
     close(FN%unit_u_sign); close(FN%unit_U_order); close(FN%unit_herit)
     close(FN%unit_inherit); close(FN%unit_D_out_for); close(FN%unit_D_out_unf)
     close(FN%unit_dbg); close(FN%unit_S_in)
     close(FN%unit_S_out_unf); close(FN%unit_S_out_for); close(FN%unit_D_in)
   end subroutine finalize
end module init_final
