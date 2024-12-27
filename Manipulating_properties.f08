module changing_parameters
  use system_parameters
  use transfrm_info
  use conversion_routines
  use Energy_module
  use gradient_module
  
  implicit none
  
  real*8 :: ppS=0.005d0 !ppI=pulling or pushing Stepsize
  
contains
  
  subroutine update_node_cnnctd()
    implicit none
    
    
    
  end subroutine update_node_cnnctd
  
  
  subroutine update_double_node()
    implicit none

    
  end subroutine update_double_node
  
  
  subroutine update_count_this_dn()
    
    
  end subroutine update_count_this_dn
  
  
  
  subroutine change_k_spr_of_one_spr(spr_nmbr,hm)
    implicit none
    integer :: spr_nmbr
    real*8  :: hm
    
    k_spr(spr_nmbr) = (1.d0 + hm) * k_spr(spr_nmbr) 
    
  end subroutine change_k_spr_of_one_spr
  
  
  subroutine change_k_spr_of_mltipl_spr(which_sprs,spr_cnt,hm)
    implicit none
    integer :: spr_cnt
    integer :: which_sprs(1:spr_cnt)
    real*8  :: hm(1:spr_cnt)
    integer :: i
    
    do i = 1,spr_cnt
       k_spr(which_sprs(i)) = (1.d0 + hm(i)) * k_spr(which_sprs(i))
    enddo
    
  end subroutine change_k_spr_of_mltipl_spr
  
  
  
  subroutine change_k_spr_of_a_typ_of_spr(which_typ,typ_cnt,hm)
    implicit none
    integer :: typ_cnt
    integer :: which_typ
    real*8  :: hm(1:typ_cnt)
    integer :: i,count
    
    count = 1
    
    do i = 1,N_spr
       if (typ_spr(i) .eq. which_typ) then
          k_spr(i) = (1.d0 + hm(count)) * k_spr(i)
          count = count + 1
       endif
    enddo
    
  end subroutine change_k_spr_of_a_typ_of_spr
  
  subroutine switch_off_all_sprs_except_some(spr_cnt,which_sprs)
    implicit none
    integer :: spr_cnt
    integer :: which_sprs(1:spr_cnt)
    real*8  :: k_tmp(1:N_spr)
    
    integer :: i
    
    k_tmp = 0.0d0
    k_tmp = k_spr
    k_spr = 0.0d0
    
    do i = 1,spr_cnt
       k_spr(which_sprs(i)) = k_tmp(which_sprs(i))
    enddo
    
  end subroutine switch_off_all_sprs_except_some
  
  !****************************************************************!
  
  subroutine change_l0_of_one_spr(spr_nmbr,hm)
    implicit none
    integer :: spr_nmbr
    real*8  :: hm
    
    l0(spr_nmbr) = (1.d0 + hm) * l0(spr_nmbr) 
    
  end subroutine change_l0_of_one_spr
  
  
  subroutine change_l0_of_mltipl_spr(which_sprs,spr_cnt,hm)
    implicit none
    integer :: spr_cnt
    integer :: which_sprs(1:spr_cnt)
    real*8  :: hm(1:spr_cnt)
    integer :: i
    
    do i = 1,spr_cnt
       l0(which_sprs(i)) = (1.d0 + hm(i)) * l0(which_sprs(i))
    enddo
    
  end subroutine change_l0_of_mltipl_spr
  
  
  subroutine change_l0_of_a_typ(which_typ,typ_cnt,hm)
    implicit none
    integer :: typ_cnt
    integer :: which_typ
    real*8  :: hm(1:typ_cnt)
    integer :: i,count
    
    count = 1
    
    do i = 1,N_spr
       if (typ_spr(i) .eq. which_typ) then
          l0(i) = (1.d0 + hm(count)) * l0(i)
          count = count + 1
       endif
    enddo
    
    
  end subroutine change_l0_of_a_typ
  
  
  !***************************************************************!
  
  
  subroutine update_spr_typ(spr_nmbr,new_typ)
    implicit none
    integer :: spr_nmbr
    integer :: new_typ
    
    typ_spr(spr_nmbr) = new_typ 
    
  end subroutine update_spr_typ
  
  
  !****************************************************************!

  subroutine change_k_area_of_a_cell(cell_no,hm)
    implicit none
    integer :: cell_no
    real*8  :: hm
    
    k_area(cell_no) = (1.d0 + hm) * k_area(cell_no)
    
  end subroutine change_k_area_of_a_cell
  
  
  subroutine change_k_area_of_mltipl_cell(which_cells,cell_cnt,hm)
    implicit none
    integer :: cell_cnt
    integer :: which_cells(1:cell_cnt)
    real*8  :: hm(1:cell_cnt)
    integer :: i
    
    do i = 1,cell_cnt
       k_area(which_cells(i)) = (1.d0 + hm(i)) * k_area(which_cells(i)) 
    enddo
    
  end subroutine change_k_area_of_mltipl_cell
  
  
  subroutine change_k_area_of_a_typ_of_cell(which_typ,typ_cnt,hm)
    implicit none
    integer :: typ_cnt
    integer :: which_typ
    real*8  :: hm(1:typ_cnt)
    integer :: i,count
    
    count = 1
    
    do i = 1,N_cell
       if (typ_area(i) .eq. which_typ) then
          k_area(i) = (1.d0 + hm(count)) * k_area(i)
          count = count + 1
       endif
    enddo
    
  end subroutine change_k_area_of_a_typ_of_cell
  
  
  subroutine switch_off_all_areas_except_some(area_cnt,which_areas)
    implicit none
    integer :: area_cnt
    integer :: which_areas(1:area_cnt)
    real*8  :: k_tmp(1:N_cell)
    
    integer :: i
    
    k_tmp  = 0.0d0
    k_tmp  = k_area
    k_area = 0.0d0
    
    do i = 1,area_cnt
       k_area(which_areas(i)) = k_tmp(which_areas(i))
    enddo
    
    
  end subroutine switch_off_all_areas_except_some
  
  
  !********************************************************************!
  
  subroutine change_A0_of_one_cell(cell_no,hm)
    implicit none
    integer :: cell_no
    real*8  :: hm !hm = how much

    A0(cell_no) = (1+hm)*A0(cell_no)
    
  end subroutine change_A0_of_one_cell
  
  
  subroutine change_A0_of_mltipl_cell(which_cells,cell_cnt,hm)
    implicit none
    integer :: cell_cnt
    integer :: which_cells(1:cell_cnt)
    real*8  :: hm(1:cell_cnt)
    integer :: i
    
    do i = 1,cell_cnt
       A0(which_cells(i)) = (1+hm(i)) * A0(which_cells(i))
    enddo
    
  end subroutine change_A0_of_mltipl_cell
  
  subroutine change_A0_of_a_type_of_cell(which_typ,typ_cnt,hm)
    implicit none
    integer :: which_typ
    integer :: typ_cnt
    real*8  :: hm(1:typ_cnt)
    integer :: i,count
    
    count = 1
    
    do i = 1,N_cell
       if (typ_area(i) .eq. which_typ) then
          A0(i) = (1.d0 + hm(count)) * A0(i)
          count = count + 1
       endif
    enddo
    
  end subroutine change_A0_of_a_type_of_cell
  
  subroutine change_A0_ofAPair_frombottom(pairNo,hmuch)
    implicit none
    integer, intent(in) :: pairNo
    real*8,  intent(in) :: hmuch
    
    real*8  :: hm(1:2)
    integer :: reducing_No(1:2)
    integer :: cellNo
    integer :: i
    
    reducing_No(1) = 2*PairNo-2
    reducing_No(2) = 2*PairNo-1
    
    hm(1:2) = hmuch
    
    do i = 1,2
       cellNo = N_cell - reducing_No(i)
       !write(*,*) cellNo,"cellNo"
       !if (i==2) stop
       
       A0(cellNo) = (1.0d0+hm(i))*A0(cellNo)
       
    enddo
    
  end subroutine change_A0_ofAPair_frombottom
  
  subroutine get_k_sprs_of_rectangular_shape(cell_no)
    
    implicit none
    integer :: cell_no
    integer :: cell_typ
    
    integer :: spr1,spr2,spr3,spr4
    real*8  :: numertr,dnumertr 
    real*8  :: res
    
    cell_typ = typ_area(cell_no)
    !write(*,*) cell_typ,"cell typ"
    
    if (cell_typ.ne.1 .AND. cell_typ.ne.4) then !1=lft_rght cell, 4=clft_rght_bot
       write(*,*) "Not Rectanguler"
       stop
    endif
    
    if (cell_no.eq.1 .or. cell_no.eq.(nsl+1)) then
       write(*,*) "Not Comapatible for boundary cell yet"
       stop
    endif
    
    
    !spr1,spr2 are vertical ones, spr3,spr4 are the horizontal ones
    
    if (cell_typ .eq. 1) then
       spr1 = area_spr(cell_no,1)
       spr2 = area_spr(cell_no,4)
       spr3 = area_spr(cell_no,2)
       spr4 = area_spr(cell_no,3)
    elseif (cell_typ .eq. 4) then
       spr1 = area_spr(cell_no,2)
       spr2 = area_spr(cell_no,3)
       spr3 = area_spr(cell_no,1)
       spr4 = area_spr(cell_no,4)
    endif
    
    !VERTICAL FORCE
    
    numertr  = k_area(cell_no)*(A0(cell_no)-A(cell_no))*l(spr3) 
    write(*,*) k_area(cell_no) , A0(cell_no), A(cell_no),l(spr3), numertr,"k_area,A0,A,l(spr3),num"
    
    dnumertr = denumerator(spr1,spr2)
    
    write(*,*) numertr,dnumertr,"num,denum"
    res = numertr / dnumertr
    
    k_spr(spr1) = res
    k_spr(spr2) = res
    
    !HORIZONTAL FORCE
    
    numertr  = k_area(cell_no) * (A0(cell_no) - A(cell_no)) * l(spr1)
    write(*,*) k_area(cell_no) , A0(cell_no), A(cell_no),l(spr1),"k_area,A0,A,l(spr1)"
    
    dnumertr = denumerator(spr3,spr4)
    
    write(*,*) numertr,dnumertr,"num,denum"
    res = 2.0*numertr / dnumertr
    
    k_spr(spr3) = res
    k_spr(spr4) = res
    
  contains
    
    real*8 function denumerator(sprA,sprB)
      implicit none
      integer :: sprA,sprB
      
      denumerator = l(sprA) - l0(sprA) + l(sprB) - l0(sprB)
      
      write(*,*) sprA,l(sprA),l0(sprA),"sprA,l,l0"
      write(*,*) sprB,l(sprB),l0(sprB),"sprB,l,l0"
      write(*,*) denumerator,"dnumertr"
      
      if (abs(denumerator-0.00) .lt. 1e-15) then
         write(*,*) "Denumertetor cant be zero,sbrtn:get_k_sprs_rect_shape,flnm:Manipulating_prop"
         stop
      endif
      
    end function denumerator
    
  end subroutine get_k_sprs_of_rectangular_shape
  
  
  
  subroutine get_spr_typ_cnt(typ_no,typ_cnt)
    implicit none
    integer :: typ_no
    integer :: typ_cnt
    integer :: i
    
    typ_cnt = 0
    
    do i = 1,N_spr
       if (typ_spr(i) .eq. typ_no) then
          typ_cnt = typ_cnt + 1
       endif
    enddo
    
  end subroutine get_spr_typ_cnt
  
  
  subroutine get_area_typ_cnt(which_typ,typ_cnt)
    implicit none
    integer :: which_typ
    integer :: typ_cnt
    integer :: i
    
    typ_cnt = 0
    
    do i = 1,N_cell
       if (typ_area(i) .eq. which_typ) then
          typ_cnt = typ_cnt + 1
       endif
    enddo
    
  end subroutine get_area_typ_cnt
  
  
  subroutine change_fctr_of_bot4_cells
    implicit none
    integer :: cell_cnt
    integer, allocatable :: which_cells(:)
    integer :: i
    
    cell_cnt = 4
    
    allocate(which_cells(1:cell_cnt))
    
    which_cells(1) = N_cell - 3
    which_cells(2) = N_cell - 2
    which_cells(3) = N_cell - 1
    which_cells(4) = N_cell
    
    write(*,*) which_cells(1:4),"which_cells"
    
    fctr(which_cells(1)) = 1.05d0 
    fctr(which_cells(2)) = 1.05d0 
    fctr(which_cells(3)) = 1.05d0 
    fctr(which_cells(4)) = 1.05d0
    
    do i = 1,cell_cnt
       A0(which_cells(i)) = fctr(which_cells(i)) * A0(which_cells(i))
       write(*,*) A0(which_cells(i)),which_cells(i),"insd change_fctr_of_bot4_cells"
    enddo
    
  end subroutine change_fctr_of_bot4_cells
  
  
  subroutine shortening_the_pulling_sprs(hmuch)
    
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: pulling_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 3) then
             pulling_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) pulling_spr(1:2),"pulling sprs"
    
    write(*,*) l0(pulling_spr(1)),l0(pulling_spr(2)),"l0 pulling spr bfr shortening"
    
    do i = 1,2
       spr_no = pulling_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(pulling_spr(1)),l0(pulling_spr(2)),"l0 pulling spr aft shortening"
    
  end subroutine shortening_the_pulling_sprs
  
  
  
  subroutine shortening_the_pulley_sprs(hmuch)
    
    implicit none
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: pulley_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr
       if (i.eq.2) cell_no = crght_top_cell_nmbr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 6) then
             pulley_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) pulley_spr(1:2),"pulley sprs"
    write(*,*) l0(pulley_spr(1)),l0(pulley_spr(2)),"l0 pulley spr bfr shortening"
    
    do i = 1,2
       spr_no = pulley_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(pulley_spr(1)),l0(pulley_spr(2)),"l0 pulley spr aft shortening"
    
    
  end subroutine shortening_the_pulley_sprs
  
  
  subroutine shortening_SecndMidSpr_frmbot(hmuch)
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: middle_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)

    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2*(ncclb-1)
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2*(nccrb-1)
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 3) then
             middle_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) middle_spr(1:2),"middle sprs"
    
    write(*,*) l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr bfr shortening"
       
    do i = 1,2
       spr_no = middle_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr aft shortening"
    
    
  end subroutine shortening_SecndMidSpr_frmbot
  
  
  
  subroutine shortening_ThrdMidSpr_frmbot(hmuch)
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: middle_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2*(ncclb-2)
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2*(nccrb-2)
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 3) then
             middle_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) middle_spr(1:2),"middle sprs"
    
    write(*,*) l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr bfr shortening"
       
    do i = 1,2
       spr_no = middle_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr aft shortening"
    
    
  end subroutine shortening_ThrdMidSpr_frmbot
  
  
  subroutine shortening_SecndSideSpr_frmbot(hmuch)
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: side_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2*(ncclb-1)
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2*(nccrb-1)
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 4) then
             side_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) side_spr(1:2),"side sprs"
    
    write(*,*) l0(side_spr(1)),l0(side_spr(2)),"l0 side spr bfr shortening"
       
    do i = 1,2
       spr_no = side_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(side_spr(1)),l0(side_spr(2)),"l0 side spr aft shortening"
    
    
  end subroutine shortening_SecndSideSpr_frmbot
  
  
  subroutine shortening_ThrdSideSpr_frmbot(hmuch)
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: side_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2*(ncclb-2)
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2*(nccrb-2)
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 4) then
             side_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) side_spr(1:2),"side sprs"
    
    write(*,*) l0(side_spr(1)),l0(side_spr(2)),"l0 side spr bfr shortening"
       
    do i = 1,2
       spr_no = side_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(side_spr(1)),l0(side_spr(2)),"l0 side spr aft shortening"
    
  end subroutine shortening_ThrdSideSpr_frmbot
  
  subroutine shortening_ThrdConnSpr_frmbot(hmuch)
    implicit none
    real*8,intent(in) :: hmuch

    integer :: cell_no,spr_no
    integer :: conn_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr + 2*(ncclb-2)
       if (i.eq.2) cell_no = crght_top_cell_nmbr + 2*(nccrb-2)
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (j==jmax) then
             conn_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) conn_spr(1:2),"conn sprs"
    
    write(*,*) l0(conn_spr(1)),l0(conn_spr(2)),"l0 conn spr bfr shortening"
       
    do i = 1,2
       spr_no = conn_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(conn_spr(1)),l0(conn_spr(2)),"l0 conn spr aft shortening"
    
    
  end subroutine shortening_ThrdConnSpr_frmbot
  
  
  subroutine shortening_ApicalSpr_nxtToPulleySpr(hmuch)
    implicit none
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: ApiclSpr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = ncl
       if (i.eq.2) cell_no = ncl + ncr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          write(*,*) spr_no,":spr_no:"
          
          if (typ_spr(spr_no) .eq. 3) then
             ApiclSpr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) ApiclSpr(1:2),"ApicalSpr"
    write(*,*) l0(ApiclSpr(1:2)),"l0 ApiclSpr bfr shortening"
    
    do i = 1,2
       spr_no = ApiclSpr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(ApiclSpr(1:2)),"l0 ApiclSpr aft shortening"
    
    
  end subroutine shortening_ApicalSpr_nxtToPulleySpr
  
  subroutine shortening_Diagonal_Sprs(hmuch)
    
    implicit none
    real*8,intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: diag_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr
       if (i.eq.2) cell_no = crght_top_cell_nmbr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 7) then
             diag_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) diag_spr(1:2),"diagonal sprs"
    
    write(*,*) l0(diag_spr(1)),l0(diag_spr(2)),"l0 diagonal spr bfr shortening"
    
    do i = 1,2
       spr_no = diag_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(diag_spr(1)),l0(diag_spr(2)),"l0 diagonal spr aft shortening"
    
    
    
  end subroutine shortening_Diagonal_Sprs
  
  
  subroutine shortening_basalFurrow_pair2(hmuch)
    implicit none
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: basal_furrowP2(1:2)
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = clft_top_cell_nmbr
       if (i.eq.2) cell_no = crght_top_cell_nmbr
       
       write(*,*) cell_no,i,"cell_no,i"

       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          write(*,*) spr_no,":spr_no:"
          
          if (typ_spr(spr_no) .eq. 8) then
             basal_furrowP2(i) = spr_no
             write(*,*) spr_no,"spr_no"
          endif
          
       enddo
       
    enddo

    write(*,*) basal_furrowP2(1:2),"basal furrow pair2 sprs"
    write(*,*) l0(basal_furrowP2(1)),l0(basal_furrowP2(2)),"l0 pulley spr bfr shortening"

    do i = 1,2
       spr_no = basal_furrowP2(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(basal_furrowP2(1)),l0(basal_furrowP2(2)),"l0 pulley spr aft shortening"
    
    
  end subroutine shortening_basalFurrow_pair2
  
  
  subroutine shortening_VertSprBfrPulley(hmuch)
    implicit none
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: VertSpr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = ncl
       if (i.eq.2) cell_no = ncl + ncr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          write(*,*) spr_no,":spr_no:"
          
          if (j==jmax) then
             VertSpr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) VertSpr(1:2),"Vertical Spr"
    write(*,*) l0(VertSpr(1:2)),"l0 VertSpr bfr shortening"
    
    
    do i = 1,2
       spr_no = VertSpr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(VertSpr(1:2)),"l0 VertSpr aft shortening"
    
  end subroutine shortening_VertSprBfrPulley
  
  subroutine shortening_BasalSprBfrPulley(hmuch)
    implicit none
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: BasalSpr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i.eq.1) cell_no = ncl
       if (i.eq.2) cell_no = ncl + ncr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       if (i.eq.1) jmax = area_spr(cell_no,0)
       if (i.eq.2) jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          write(*,*) spr_no,":spr_no:"
          
          if (typ_spr(spr_no).eq.4) then
             BasalSpr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) BasalSpr(1:2),"Basal Spr"
    write(*,*) l0(BasalSpr(1:2)),"l0 BasalSpr bfr shortening"
    
    
    do i = 1,2
       spr_no = BasalSpr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*) l0(BasalSpr(1:2)),"l0 BasalSpr aft shortening"
    
    
  end subroutine shortening_BasalSprBfrPulley
  
  
  subroutine shortening_MidSpr_frmtop(serial,hmuch)
    
    implicit none
    integer,intent(in) :: serial
    real*8, intent(in) :: hmuch
    
    integer :: cell_no,spr_no
    integer :: middle_spr(1:2)
    
    integer :: i,j,jmax
    real*8  :: hm(1:2)
    
    jmax = -10
    hm(1:2) = hmuch
    
    do i = 1,2 !2 is not for dmnsn, 1 for lft,1 for rght
       
       if (i==1) cell_no = clft_top_cell_nmbr  + 2*serial
       if (i==2) cell_no = crght_top_cell_nmbr + 2*serial
       
       write(*,*) cell_no,i,"cell_no,i"
       
       jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) .eq. 3) then
             middle_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) middle_spr(1:2),"middle sprs"
    
    write(*,*)l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr bfr shortening"
       
    do i = 1,2
       spr_no = middle_spr(i)
       
       l0(spr_no) = (1.0d0-hm(i)) * l0(spr_no)
       
    enddo
    
    write(*,*)l0(middle_spr(1)),l0(middle_spr(2)),"l0 middle spr aft shortening"
    
  end subroutine shortening_MidSpr_frmtop
  
  subroutine get_the_MidSpr_frmtop(serial,middle_spr)
    implicit none
    integer, intent(in) :: serial
    integer, intent(inout) :: middle_spr(1:2)

    integer :: i,j,jmax
    integer :: cell_no,spr_no
    
    jmax = -10
    
    do i = 1,2
       
       if (i==1) cell_no = clft_top_cell_nmbr  + 2*serial
       if (i==2) cell_no = crght_top_cell_nmbr + 2*serial
       
       write(*,*) cell_no,i,"cell_no,i"
       
       jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) == 3) then
             middle_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) middle_spr(1:2),"middle sprs"
    
  end subroutine get_the_MidSpr_frmtop
  
  
  subroutine get_the_DiagSpr(diag_spr)
    implicit none
    integer, intent(inout) :: diag_spr(1:2)
    
    integer :: i,j,jmax
    integer :: cell_no,spr_no
    
    
    jmax = -10
    
    do i = 1,2
       
       if (i==1) cell_no = clft_top_cell_nmbr
       if (i==2) cell_no = crght_top_cell_nmbr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) == 7) then
             diag_spr(i) = spr_no
          endif
          
       enddo
       
    enddo
    
    write(*,*) diag_spr(1:2),"diagonal sprs"
    
  end subroutine get_the_DiagSpr
  
  subroutine get_the_ApicalSpr_bfrPulley(apcl_spr)
    implicit none
    integer, intent(out) :: apcl_spr(1:2)
    
    integer :: i,j,jmax
    integer :: cell_no,spr_no
    
    
    jmax = -10
    
    do i = 1,2
       
       if (i==1) cell_no = ncl
       if (i==2) cell_no = ncl + ncr
       
       write(*,*) cell_no,i,"cell_no,i"
       
       jmax = area_spr(cell_no,0)
       
       do j = 1,jmax
          
          spr_no = area_spr(cell_no,j)
          
          if (typ_spr(spr_no) == 3) then
             apcl_spr(i) = spr_no
          endif
          
       enddo
       
    enddo

    open(unit=125,file='Apcl_spr.dat')
    write(125,fmt=*) apcl_spr(1:2),"apical sprs bfr pulley"
    close(125)
    
  end subroutine get_the_ApicalSpr_bfrPulley
  
  subroutine get_the_ApclSpr(apcl_spr)
    implicit none
    integer, intent(out) :: apcl_spr(1:N_cell)
    
    integer :: Top_cell,spr_nm
    integer :: i,j,jmax
    
    open(unit=120,file='Apical_spr.dat')
    
    Top_cell = ncl+ncr

    apcl_spr = -1
    
    do i = 1,N_cell
       jmax = area_spr(i,0)
       
       do j = 1,jmax
          spr_nm = area_spr(i,j)
          
          if (i.le.Top_cell) then
             if (typ_spr(spr_nm)==1.or.typ_spr(spr_nm)==3) then
                apcl_spr(i) = spr_nm
             endif
                
          elseif (i.gt.Top_cell.and.i.le.(Top_cell+2)) then
             if (typ_spr(spr_nm)==6) then
                apcl_spr(i) = spr_nm
             endif
             
          else
             if (typ_spr(spr_nm)==3) then
                apcl_spr(i) = spr_nm
             endif
          endif
          
       enddo
       write(120,fmt=*) apcl_spr(i),i,"apcl_spr,cell"
    enddo
    
    close(120)
    
  end subroutine get_the_ApclSpr
  
  subroutine shortening_ApclSpr(apcl_spr,hmuch)
    implicit none
    
    integer,intent(in) :: apcl_spr(1:N_cell)
    real*8 ,intent(in) :: hmuch
    
    real*8  :: hm(1:N_cell)
    integer :: i,spr_nm
    
    open(unit=141,file='ApclSpr_BfrAft.dat',position="append")
    write(141,fmt=*) "l0(apcl_spr),apcl_spr(i),i"
    
    do i = 1,N_cell
       
       if (apcl_spr(i).ne.(-1)) then
          
          write(141,fmt=*) l0(apcl_spr(i)),apcl_spr(i),i,"bfr"
          
          spr_nm = apcl_spr(i)
          hm(i)  = hmuch
          l0(spr_nm) = (1.0d0-hm(i))*l0(spr_nm)
          
          write(141,fmt=*) l0(apcl_spr(i)),apcl_spr(i),i,"aft"
          
       endif
       
    enddo
    
    close(141)
    
  end subroutine shortening_ApclSpr
  
  
  subroutine get_the_BaslSpr(basl_spr)
    implicit none
    integer, intent(out) :: basl_spr(1:N_cell)
    
    integer :: Top_cell,spr_nm
    integer :: i,j,jmax
    
    open(unit=120,file='Basl_spr.dat')
    
    Top_cell = ncl+ncr
    
    do i = 1,N_cell
       jmax = area_spr(i,0)
       
       do j = 1,jmax
          spr_nm = area_spr(i,j)
          
          if (i.le.Top_cell) then
             if (typ_spr(spr_nm)==2.or.typ_spr(spr_nm)==4) then
                basl_spr(i) = spr_nm
             endif
                
          elseif (i.gt.Top_cell.and.i.le.(Top_cell+2)) then
             if (typ_spr(spr_nm)==8) then
                basl_spr(i) = spr_nm
             endif
             
          else
             if (typ_spr(spr_nm)==4) then
                basl_spr(i) = spr_nm
             endif
          endif
       
       enddo
       write(120,fmt=*) basl_spr(i),i,"basl_spr,cell"
    enddo
    
    close(120)
    
  end subroutine get_the_BaslSpr
  
  subroutine get_the_LtrlSpr(ltrl_spr)
    implicit none
    integer, intent(out) :: ltrl_spr(1:N_cell)
    
    integer :: Top_cell,spr_nm
    integer :: i,j,jmax
    
    open(unit=120,file='ltrl_spr.dat')
    
    Top_cell = ncl+ncr
    
    do i = 1,N_cell
       jmax = area_spr(i,0)
       
       do j = 1,jmax
          spr_nm = area_spr(i,j)
          
          if (i.le.Top_cell) then
             if (typ_spr(spr_nm)==5) then
                ltrl_spr(i) = spr_nm
             endif
                
          elseif (i.gt.Top_cell.and.i.le.(Top_cell+2)) then
             if (typ_spr(spr_nm)==7) then
                ltrl_spr(i) = spr_nm
             endif
             
          else
             if (typ_spr(spr_nm)==5) then
                ltrl_spr(i) = spr_nm
             endif
          endif
       
       enddo
       write(120,fmt=*) ltrl_spr(i),i,"ltrl_spr,cell"
    enddo
    
    close(120)
    
  end subroutine get_the_LtrlSpr
  
  subroutine shortening_sprPair(Pair,hmuch)
    implicit none
    integer, intent(in) :: Pair(1:2)
    real*8 , intent(in) :: hmuch
    
    real*8  :: hm(1:2)
    integer :: i,spr_nm
    
    write(*,*) Pair(1:2),"spr pair"
    write(*,*) l0(Pair(1)),l0(Pair(2)),"l0s bfr shortn"
    
    do i = 1,2
       spr_nm = Pair(i)
       hm(i)  = hmuch
       l0(spr_nm) = (1.0d0-hm(i))*l0(spr_nm)
    enddo
    
    write(*,*) l0(Pair(1)),l0(Pair(2)),"l0s aft shortn"
    
  end subroutine shortening_sprPair
  
  subroutine shortening_singleSpr(spr_nm,hmuch)
    implicit none
    integer, intent(in) :: spr_nm
    real*8 , intent(in) :: hmuch
    
   
    write(*,*) l0(spr_nm),spr_nm,"l0 bfr shortn"
    
    l0(spr_nm) = (1.0d0-hmuch)*l0(spr_nm)
    
    write(*,*) l0(spr_nm),spr_nm,"l0s aft shortn"
    
  end subroutine shortening_singleSpr
  
  subroutine switch_off_alpha_beta_excpt_bottom4_cells
    implicit none
    integer :: N_actvCellOptm,N_semiActvCellOptm
    integer,allocatable :: actv_cellOptm(:),semi_actvCellOptm(:)
    logical :: actv
    
    integer :: spr_nm,area_nm
    integer :: i,j
    
    spr_nm = -10 ; area_nm = -10
    
    N_actvCellOptm     = 4
    N_semiActvCellOptm = 2
    
    allocate(actv_cellOptm(1:N_actvCellOptm))
    allocate(semi_actvCellOptm(1:N_semiActvCellOptm))
    
    actv = .False.
    
    actv_cellOptm(1) = N_cell - 3
    actv_cellOptm(2) = N_cell - 2
    actv_cellOptm(3) = N_cell - 1
    actv_cellOptm(4) = N_cell
    
    semi_actvCellOptm(1) = N_cell - 5
    semi_actvCellOptm(2) = N_cell - 4
    
    do i = 1,N_cell
       actv = .False.
       
       do j = 1,N_actvCellOptm
          if (i.eq.actv_cellOptm(j)) then
             actv = .True.
             exit
          endif
       enddo
       
       if (actv.eqv..False.) then
          beta(i) = 0.0d0
          area_nm = i
       endif
       
       if (area_nm.ne.semi_actvCellOptm(1) .AND. area_nm.ne.semi_actvCellOptm(2)) then
          
          do j = 1,max_spr_area
             spr_nm = area_spr(area_nm,j)
          
             if (spr_nm .ne. (-1)) alpha(spr_nm) = 0.0d0 
          
          enddo
          
       elseif (area_nm.eq.semi_actvCellOptm(1) .or. area_nm.eq.semi_actvCellOptm(2))then
          
          do j = 1,max_spr_area
             spr_nm = area_spr(area_nm,j)
             !write(*,*)
             if (typ_spr(spr_nm) .ne. 5) then
                if (spr_nm .ne. (-1)) alpha(spr_nm) = 0.0d0 
             endif
             
          enddo
          
       endif
       
       
    enddo
    
    do i = 1,N_cell
       write(*,*) beta(i),i,"alpha of cell_nm"
    enddo
    
    
    do i = 1,N_spr
       write(*,*) alpha(i),i,"beta of spr_nm"
    enddo
    
    !stop
    
  end subroutine switch_off_alpha_beta_excpt_bottom4_cells
  
  
  subroutine switchOff_optimizer_for_bottom4_CellSpr
    implicit none
    
    integer :: N_actvCellOptm
    integer,allocatable :: actv_cellOptm(:)
    
    integer :: spr_nm,area_nm
    integer :: nspr_area !no of spr cnn to that area
    integer :: i,j
    
    spr_nm = -10 ; area_nm = -10
    
    N_actvCellOptm = 4
    
    allocate(actv_cellOptm(1:N_actvCellOptm))
    
    actv_cellOptm(1) = N_cell - 3
    actv_cellOptm(2) = N_cell - 2
    actv_cellOptm(3) = N_cell - 1
    actv_cellOptm(4) = N_cell
    
    do i = 1,N_actvCellOptm
       area_nm = actv_cellOptm(i)
       optmCell(area_nm) = 0
       
       nspr_area = area_spr(area_nm,0)
       
       do j = 1,nspr_area
          spr_nm = area_spr(area_nm,j)
          optmSpr(spr_nm) = 0
       enddo
       
    enddo
    
    
    do i = 1,N_cell
       write(*,*) optmCell(i),i,"optmCell of cell_nm"
    enddo
    
    
    do i = 1,N_spr
       write(*,*) optmSpr(i),i,"optmSpr of spr_nm"
    enddo
    
    
  end subroutine switchOff_optimizer_for_bottom4_CellSpr
  
  
  subroutine switchOff_optimizer_for_bot6_Cell_and_bot4Spr
    implicit none
    
    integer :: N_actvCellOptm,N_semiactvCellOptm
    integer,allocatable :: actv_cellOptm(:),semiactv_cellOptm(:)
    
    integer :: spr_nm,area_nm
    integer :: nspr_area !no of spr cnn to that area
    integer :: i,j
    
    spr_nm = -10 ; area_nm = -10
    
    N_actvCellOptm = 4 ; N_semiactvCellOptm = 2
    
    allocate(actv_cellOptm(1:N_actvCellOptm),semiactv_cellOptm(1:N_semiactvCellOptm))
    
    actv_cellOptm(1) = N_cell - 3
    actv_cellOptm(2) = N_cell - 2
    actv_cellOptm(3) = N_cell - 1
    actv_cellOptm(4) = N_cell
    
    do i = 1,N_actvCellOptm
       area_nm = actv_cellOptm(i)
       optmCell(area_nm) = 0
       
       nspr_area = area_spr(area_nm,0)
       
       do j = 1,nspr_area
          spr_nm = area_spr(area_nm,j)
          optmSpr(spr_nm) = 0
       enddo
       
    enddo
    
    semiactv_CellOptm(1) = N_cell-5
    semiactv_CellOptm(2) = N_cell-4
    
    do i = 1,N_semiactvCellOptm
       area_nm = semiactv_cellOptm(i)
       optmCell(area_nm) = 0
    enddo
    
    do i = 1,N_cell
       write(*,*) optmCell(i),i,"optmCell of cell_nm"
    enddo
    
    
    do i = 1,N_spr
       write(*,*) optmSpr(i),i,"optmSpr of spr_nm"
    enddo
    
  end subroutine switchOff_optimizer_for_bot6_Cell_and_bot4Spr
  
  
  
  subroutine switchOff_optimizer_for_allA0
    implicit none
    integer :: spr_nm,area_nm
    integer :: nspr_area !no of spr cnn to that area
    integer :: i,j
    
    spr_nm = -10 ; area_nm = -10
    
    
    do i = 1,N_cell
       optmCell(i) = 0
       
       !nspr_area = area_spr(area_nm,0)
       
       !do j = 1,nspr_area
        !  spr_nm = area_spr(area_nm,j)
         ! optmSpr(spr_nm) = 0
       !enddo
       
    enddo
    
    
    do i = 1,N_cell
       write(*,*) optmCell(i),i,"optmCell of cell_nm"
    enddo

    
    do i = 1,N_spr
       write(*,*) optmSpr(i),i,"optmSpr of spr_nm"
    enddo
    
    
  end subroutine switchOff_optimizer_for_allA0
  
  
  
  
  
  subroutine switch_on_alphabeta_ofFirstPair_ofbottomCell
    implicit none
    
    integer :: N_actvCellOptm
    integer,allocatable :: actv_cellOptm(:)
    
    integer :: spr_nm,area_nm
    integer :: nspr_area !no of spr cnn to that area
    integer :: i,j
    
    spr_nm = -10 ; area_nm = -10
    
    N_actvCellOptm = 2
    
    allocate(actv_cellOptm(1:N_actvCellOptm))
    
    actv_cellOptm(1) = N_cell - 5
    actv_cellOptm(2) = N_cell - 4
        
    
    do i = 1,N_actvCellOptm
       area_nm = actv_cellOptm(i)
       beta(area_nm) = 1.0d0
       
       nspr_area = area_spr(area_nm,0)
       
       do j = 1,nspr_area
          spr_nm = area_spr(area_nm,j)
          alpha(spr_nm) = 1.0d0
       enddo
       
    enddo

    
    do i = 1,N_cell
       write(*,*) beta(i),i,"alpha of cell_nm"
    enddo
    
    
    do i = 1,N_spr
       write(*,*) alpha(i),i,"beta of spr_nm"
    enddo
    
    
  end subroutine switch_on_alphabeta_ofFirstPair_ofbottomCell
  
  
  subroutine switch_off_the_PenF
    implicit none
    integer :: N_actvCellOptm,N_semiActvCellOptm
    integer,allocatable :: actv_cellOptm(:),semi_actvCellOptm(:)
    logical :: actv
    
    integer :: spr_nm,area_nm
    integer :: i,j
    
    alpha(1:N_spr) = 0.0d0
    beta(1:N_cell) = 0.0d0
    
    do i = 1,N_cell
       write(*,*) beta(i),i,"alpha of cell_nm"
    enddo
    
    
    do i = 1,N_spr
       write(*,*) alpha(i),i,"beta of spr_nm"
    enddo
    
    !stop
    
  end subroutine switch_off_the_PenF
  
  
  subroutine put_tension_in_horzntl_spr_lft_rght
    implicit none
    
    integer :: i
    integer :: cnt
    real*8  :: hm
    
    cnt = 0
    hm = 0.10d0
    
    do i = (cnt+1),(cnt+nsl)
       if (typ_spr(i).eq.3 .or. typ_spr(i).eq.4) then
          l0(i) = (1-hm)*l0(i)
       endif
    enddo
    
    write(*,*) l0(1:nsl),"l0_lft"
    
    cnt = nsl
    
    do i = (cnt+1),(cnt+nsr)
       if (typ_spr(i).eq.3 .or. typ_spr(i).eq.4) then
          l0(i) = (1-hm)*l0(i)
       endif
    enddo
    
    write(*,*) l0((nsl+1):(nsl+nsr)),"l0_rght"
    !stop
  end subroutine put_tension_in_horzntl_spr_lft_rght
  
  
  subroutine successiveChangein_f
    implicit none
    
    integer :: cellNo_of_fmax!it serves no purpose in current alg
    real*8  :: fmax
    real*8  :: reltv_flft(1:N_lftCells)
    real*8  :: reltv_frght(1:N_rghtCells)
    real*8  :: reltv_fcntrl(1:N_cntrlCells)
    real*8  :: reltv_fendRgn(1:N_endRgnCells)
    real*8  :: reltv_fval
    integer :: cell_no
    integer :: RgnIndctr
    
    integer :: i,j,imax,jmax
    real*8  :: fval
    
    if (stageNo==1 .and. stageType==1) then
       imax = 3
    elseif (stageNo==4 .and. stageType==1) then
       imax = 3
    elseif (stageNo==4 .and. stageType==2) then
       imax = 2
    endif
    
    
    do i = 1,imax
       
       if (i==1) then
          jmax = N_lftCells
          cellNo_of_fmax = lftSideCell(N_lftCells)
          RgnIndctr = 1
          call get_fmax_And_reltvf(fmax,reltv_flft,N_lftCells,RgnIndctr)
          
       elseif (i==2) then
          jmax = N_rghtCells
          cellNo_of_fmax = rghtSideCell(N_rghtCells)
          RgnIndctr = 2
          call get_fmax_And_reltvf(fmax,reltv_frght,N_rghtCells,RgnIndctr)
          
       elseif (i==3) then
          
          if (stageNo==1 .and. stageType==1) then
             jmax = N_cntrlCells
             cellNo_of_fmax = cntrlCell(N_cntrlCells)
             RgnIndctr = 3
             call get_fmax_And_reltvf(fmax,reltv_fcntrl,N_cntrlCells,RgnIndctr)
             
          elseif (stageNo==4 .and. stageType==1) then
             jmax = N_endRgnCells
             cellNo_of_fmax = endRgnCell(N_endRgnCells)
             write(*,*) endRgnCell(1:N_endRgnCells),"ERC"
             RgnIndctr = 3
             call get_fmax_And_reltvf(fmax,reltv_fendRgn,N_endRgnCells,RgnIndctr)
             
          endif
          
       endif
       
       
       do j = 1,jmax
          
          if (i==1) then
             cell_no = lftSideCell(j)
             reltv_fval = reltv_flft(j)
             
          elseif (i==2) then
             cell_no = rghtSideCell(j)
             reltv_fval = reltv_frght(j)
             
          elseif (i==3) then
             
             if (stageNo==1 .and. stageType==1) then
                cell_no = cntrlCell(j)
                reltv_fval = reltv_fcntrl(j)
                
             elseif (stageNo==4 .and. stageType==1) then
                cell_no = endRgnCell(j)
                reltv_fval = reltv_fendRgn(j)
                
             endif
             
          endif
          
          call get_f(cell_no,fmax,cellNo_of_fmax,reltv_fval,fval)
          !if (cell_no.ne.17) fctr(cell_no) = fval-0.25d0
          !if (cell_no==17) fctr(cell_no) = fval
          fctr(cell_no) = fval
          write(*,*) fctr(cell_no),fval,cell_no,"fctr,fval,cellno"
          write(*,*) " "
          
       enddo
       
    enddo    
    !stop
    
  end subroutine successiveChangein_f
  
  
  subroutine get_fmax_And_reltvf(fmax,reltv_f,N_cellRgn,RgnIndctr)
    implicit none
    integer, intent(in) :: N_cellRgn !lft/rght its 8,for cntrl/end 1
    integer, intent(in) :: RgnIndctr !lft=1,rght=2,cntrl/end=3
    
    real*8, intent(out)  :: fmax
    real*8, intent(out)  :: reltv_f(1:N_cellRgn)
    
    fmax = 1.20d0
    
    if (stageNo==1 .and. stageType==1) then
       
       if (RgnIndctr.le.2) then
          reltv_f(1) = 0.80d0 !0.50d0 
          reltv_f(2) = 0.75d0 !0.50d0
          reltv_f(3) = 0.70d0 !0.50d0
          reltv_f(4) = 0.65d0 !0.50d0
          reltv_f(5) = 0.60d0 !0.50d0
          reltv_f(6) = 0.55d0 !0.50d0
          reltv_f(7) = 0.50d0 !0.50d0
          reltv_f(8) = 0.50d0 !0.50d0
          
       elseif (RgnIndctr==3) then
          reltv_f(1) = 3.25d0
       endif
       
    elseif (stageNo==4) then
       
       if (strctNo==1) then
          
          if (RgnIndctr.le.2) then
             reltv_f(1) = 0.75d0
             reltv_f(2) = 0.83d0
             reltv_f(3) = 0.75d0
             reltv_f(4) = 0.50d0
             reltv_f(5) = 1.05d0
             reltv_f(6) = 1.70d0
             reltv_f(7) = 1.15d0
             reltv_f(8) = 0.90d0
             
           elseif (RgnIndctr==3) then !end cell for stg4,typ1
             reltv_f(1) = 0.85d0
          endif
          
       elseif (strctNo==2) then
          
          if (RgnIndctr.le.2) then
             reltv_f(1) = 1.15d0
             reltv_f(2) = 1.15d0
             reltv_f(3) = 1.05d0
             reltv_f(4) = 1.05d0
             reltv_f(5) = 1.80d0
             reltv_f(6) = 1.20d0
             reltv_f(7) = 1.13d0
             reltv_f(8) = 1.05d0
          elseif (RgnIndctr==3) then
             reltv_f(1) = 1.00d0
          endif
          
       endif
       
       !reltv_f = 2.0d0*reltv_f
       
    endif
    
  end subroutine get_fmax_And_reltvf
  
  
  subroutine get_f(cell_no,fmax,CellNoof_fmax,reltv_fval,fval)
    implicit none
    integer :: cell_no
    integer :: CellNoof_fmax
    real*8  :: fmax
    real*8  :: reltv_fval,fval
    
    real*8  :: srVar1 !sr=SubRoutine
    
    write(*,*) CellNoof_fmax,"Cell No of fmax"
    write(*,*) reltv_fval,fmax,A0(CellNoof_fmax),"relfval,fmax,A0-fmax"
    write(*,*) A0(cell_no),"A0-cellno"
    
    !srVar1  = (1.d0-reltv_fval)*((1.d0-fmax)*A0(CellNoof_fmax))**2
    !write(*,*) srVar1,"v1"
    
    !fval  = 1.0d0+sqrt(srVar1/(A0(cell_no))**2)
    !write(*,*) fval,"v2"
    
    fval = 1.0d0 + (reltv_fval)*(fmax-1.0d0)
    
  end subroutine get_f
  
  subroutine changeA0_using_f(fval)
    implicit none
    real*8, intent(in) :: fval(1:N_cell)
    integer :: i
    
    do i = 1,N_cell
       write(*,*) fval(i),i,"fval"
    enddo
    
    do i = 1,N_cell
       A0(i) = fval(i)*A0(i)
       write(*,*) A0(i),i,"A0_val"
    enddo
    
    !if (strctNo==1) stop
    
  end subroutine changeA0_using_f
  
  
  subroutine get_init_and_trgt_Vals_fromTwoPairCells()
    implicit none
    
    
    
  end subroutine get_init_and_trgt_Vals_fromTwoPairCells
  
  
  subroutine get_trnsfrng_cells(trns_frmCell,trns_toCell)
    implicit none
    !integer :: cnt
    integer, intent(inout) :: trns_frmCell(1:2)
    integer, intent(inout) :: trns_toCell(1:2)
    
    !cnt = ncl + ncr + ncclt + nccrt
    
    !do i = 1,ncclb
    
    !enddo
    trns_frmCell(1) = 15 ; trns_frmCell(2) = 16
    trns_toCell(1)  = 13 ; trns_toCell(2)  = 14
    
  end subroutine get_trnsfrng_cells
  
  
  
  subroutine smallChnge_and_equilibrate(trns_frmCell,trns_tocell,max_prcnt)
    implicit none
    
    integer, intent(in) :: trns_frmcell(1:2)
    integer, intent(in) :: trns_tocell(1:2)
    real*8,  intent(in) :: max_prcnt
    
    integer :: trns_frmspr,trns_tospr
    integer :: i,j,m
    integer :: jmax,jmax1,jmax2,mmax
    
    integer :: area1,area2,spr_type,spr1,spr2
    integer :: area_nm,spr_nm
    
    real*8  :: incr,prcnt
    real*8  :: A0_strt(1:2),A0_fnsh(1:2),A0_range(1:2)
    real*8  :: l0_strt(1:2,1:3),l0_fnsh(1:2,1:3),l0_range(1:2,1:3)
    
    integer :: funcChoice,N_variabls,iter
    real*8  :: ftol,fret,fpp
    
    do i = 1,2
       area1 = trns_toCell(i)
       area2 = trns_frmCell(i)
       
       A0_strt(i) = A0(area1)
       A0_fnsh(i) = A0(area2)
       
       A0_range(i) = A0_fnsh(i)-A0_strt(i)
       
       jmax1 = area_spr(area1,0)
       jmax2 = area_spr(area2,0)
       
       if (jmax1.ne.jmax2) then
          write(*,*) "spr nos must be equal in both cell"
          stop
       endif
       
       jmax = jmax1
       
       do j = 1,jmax-1
          spr_type = typ_spr(j)
          call get_sametypSprfromTwoArea(area1,area2,spr_type,spr1,spr2)
          
          trns_toSpr  = spr1
          trns_frmSpr = spr2
          
          l0_strt(i,j)  = l0(trns_toSpr)
          l0_fnsh(i,j)  = l0(trns_frmSpr)
          
          l0_range(i,j)  = l0_fnsh(i,j) - l0_strt(i,j)
       enddo
       
    enddo
    
    
    incr  = 0.10d0
    mmax  = int(max_prcnt/incr)
    
    do m = 1,mmax
       prcnt = m*incr
       
       do i = 1,2
          area_nm = trns_toCell(i)
          
          A0(area_nm) = A0_strt(i) + prcnt*A0_range(i)
          
          jmax = area_spr(area_nm,0)
          
          do j = 1,jmax-1
             spr_nm     = area_spr(area_nm,j)          
             l0(spr_nm) = l0_strt(i,j) + prcnt*l0_range(i,j)
          enddo
          
       enddo
       
       
       call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       
       funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
       call frprmn(coordntes,grd_mv,iter,funcChoice,N_variabls,ftol,fret,fpp)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)
       
       funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06    
       call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_variabls,ftol,fret,fpp)
       
    enddo
    
    
  end subroutine smallChnge_and_equilibrate
  
  
  
  
  subroutine TrgtPropChng_and_equilibrate(trns_frmCell,trns_tocell,max_prcnt)
    implicit none
    
    integer, intent(in) :: trns_frmcell(1:2)
    integer, intent(in) :: trns_tocell(1:2)
    real*8,  intent(in) :: max_prcnt
    
    integer :: trns_frmspr,trns_tospr
    integer :: i,j,m
    integer :: jmax,jmax1,jmax2,mmax
    
    integer :: area1,area2,spr_type,spr1,spr2
    integer :: area_nm,spr_nm
    
    real*8  :: incr,prcnt
    real*8  :: At_strt(1:2),At_trgt(1:2),At_rnge(1:2)
    real*8  :: Lt_strt(1:2,1:3),Lt_trgt(1:2,1:3),Lt_rnge(1:2,1:3)
    
    integer :: funcChoice,N_variabls,iter
    real*8  :: ftol,fret,fpp
    
    do i = 1,2
       area1 = trns_toCell(i)
       area2 = trns_frmCell(i)
       
       At_strt(i) = At(area1)
       At_trgt(i) = At(area2)
       
       At_rnge(i) = At_trgt(i)-At_strt(i)
       
       jmax1 = area_spr(area1,0)
       jmax2 = area_spr(area2,0)
       
       if (jmax1.ne.jmax2) then
          write(*,*) "spr nos must be equal in both cell"
          stop
       endif
       
       jmax = jmax1
       
       do j = 1,jmax-1
          spr_type = typ_spr(j)
          call get_sametypSprfromTwoArea(area1,area2,spr_type,spr1,spr2)
          
          trns_toSpr  = spr1
          trns_frmSpr = spr2
          
          Lt_strt(i,j)  = Lt(trns_toSpr)
          Lt_trgt(i,j)  = Lt(trns_frmSpr)
          
          Lt_rnge(i,j)  = Lt_trgt(i,j) - Lt_strt(i,j)
       enddo
       
    enddo
    
    
    incr  = 0.10d0
    mmax  = int(max_prcnt/incr)
    
    do m = 1,mmax
       prcnt = m*incr
       
       do i = 1,2
          area_nm = trns_toCell(i)
          
          At(area_nm) = At_strt(i) + prcnt*At_rnge(i)
          
          jmax = area_spr(area_nm,0)
          
          do j = 1,jmax-1
             spr_nm     = area_spr(area_nm,j)          
             Lt(spr_nm) = Lt_strt(i,j) + prcnt*Lt_rnge(i,j)
          enddo
          
       enddo
       
       
       if (m==1) call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       
       funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
       call frprmn(coordntes,grd_mv,iter,funcChoice,N_variabls,ftol,fret,fpp)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)
       
       funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06    
       call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_variabls,ftol,fret,fpp)
       coordntes(1:N_mvCoordnte) = coordntes_xy(1:N_mvCoordnte)
       grd_mv(1:N_mvCoordnte)    = grdmv_xy(1:N_mvCoordnte)
       call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)
       
    enddo
    
    
  end subroutine TrgtPropChng_and_equilibrate
  
  
  
  
  subroutine get_sametypSprfromTwoArea(area1,area2,spr_type,spr1,spr2)
    implicit none
    integer, intent(in)  :: area1,area2
    integer, intent(in)  :: spr_type
    integer, intent(out) :: spr1,spr2
    
    integer :: i,j,jmax
    integer :: area_nm,spr_nm
    
    do i = 1,2
       
       if(i==1) area_nm = area1
       if(i==2) area_nm = area2
       
       jmax = area_spr(area_nm,0)
       
       do j = 1,jmax-1
          spr_nm = area_spr(area_nm,j)
          
          if (typ_spr(j)==spr_type) then
             if (i==1) spr1 = spr_nm
             if (i==2) spr2 = spr_nm
          endif
       enddo
       
    enddo
    
  end subroutine get_sametypSprfromTwoArea
  
  subroutine find_theCommon_spr(area1,area2,spr_type)
    implicit none
    integer :: area1,area2
    integer :: spr_type
    integer :: j,jmax
    
    
    
    
  end subroutine find_theCommon_spr
  
  
  
  subroutine fctr_and_l0s_frstPairCntrlBot_cells(max_prcnt)
    implicit none
    real*8,intent(in) :: max_prcnt
    integer :: cell_cnt
    
    integer, allocatable :: which_cells(:)
    integer, allocatable :: A0_ini(:)
    
    real*8  :: incr,hm
    integer :: i,m
    integer :: Nstp
    integer :: funcChoice,N_variabls,iter
    real*8  :: ftol,fret,fpp
    
    cell_cnt = 2
    
    allocate(which_cells(1:cell_cnt),A0_ini(1:cell_cnt))
    
    which_cells(1) = N_cell - 5
    which_cells(2) = N_cell - 4
    
    write(*,*) which_cells(1:2),"which_cells"
    
    incr = 0.05d0
    Nstp = int(max_prcnt/incr)
    
    A0_ini(1) = A0(which_cells(1))
    A0_ini(2) = A0(which_cells(2))
    
    
    do m = 1,Nstp
       
       hm = m*incr
       
       fctr(which_cells(1)) = (1.0d0 + hm)
       fctr(which_cells(2)) = (1.0d0 + hm)
       
       !hm = hm*2.0d0
       
       do i = 1,cell_cnt
          A0(which_cells(i)) = fctr(which_cells(i)) * A0_ini(i)
          write(*,*) A0(which_cells(i)),which_cells(i),"insd change_fctr_of_bot4_cells"          
          call shorten_sprs_of_an_area(which_cells(i),hm)
       enddo
       
       funcChoice = 2 ; N_variabls = N_mvCoordnte_withl0A0 ; ftol=1.d-03
       call frprmn(coordntes,grd_mv,iter,funcChoice,N_variabls,ftol,fret,fpp)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       grdmv_xy(1:N_mvCoordnte)     = grd_mv(1:N_mvCoordnte)
       
       funcChoice = 1 ; N_variabls = N_mvCoordnte ; ftol=1.d-06
       call frprmn(coordntes_xy,grdmv_xy,iter,funcChoice,N_variabls,ftol,fret,fpp)
       
    enddo
    
  end subroutine fctr_and_l0s_frstPairCntrlBot_cells
  
  
  subroutine shorten_sprs_of_an_area(area_nm,hm)
    implicit none
    integer, intent(in) :: area_nm
    real*8 , intent(in) :: hm
    
    integer :: i,j,imax
    integer :: spr_nm
    
    imax = area_spr(area_nm,0)
    
    do i = 1,imax
       
       spr_nm = area_spr(area_nm,i)
       
       if (typ_spr(spr_nm)==3 .or. typ_spr(spr_nm)==4) then
          l0(spr_nm) = (1.0d0-(hm))*l0(spr_nm)
       else
          l0(spr_nm) = (1.0d0-(2.0d0*hm))*l0(spr_nm)
       endif
       
    enddo  
    
  end subroutine shorten_sprs_of_an_area
  
  
  
  subroutine keep_Tension_Const_inBoundarySpr(Tensions,N_bndrySpr)
    implicit none
    
    integer,intent(in) :: N_bndrySpr
    real*8 ,intent(in) :: Tensions(1:N_bndrySpr)
    
    integer :: i
    
  end subroutine keep_Tension_Const_inBoundarySpr
  
  
  subroutine decrease_kspr_and_karea_except_first3_lftrghtCells(hm)
    
    implicit none
    real*8,intent(in)  :: hm
    integer :: i
    
    
    do i = 1,N_cell
       if (i.le.3) then
          continue
       elseif ((i-ncl).gt.0 .AND.(i-ncl).le.3 ) then
          continue
       else
          write(*,*) i,"cell_nm"
          k_area(i) = (1.0d0-hm)*k_area(i)
       endif
    enddo
    
    do i = 1,N_spr
       if (i.le.(nsl-3)) then
          continue
          
       !elseif(i.gt.(nsl-3) .and. i.le.nsl)
        !  write(*,*) i,"spr_nm"
         ! k_spr(i) = (1.0d0-hm)*k_spr(i)
          
       elseif ((i-nsl).gt.0 .AND.(i-nsl).le.9 ) then
          continue
       else
          write(*,*) i,"spr_nm"
          k_spr(i) = (1.0d0-hm)*k_spr(i)
       endif
    enddo
    
    !stop
  end subroutine decrease_kspr_and_karea_except_first3_lftrghtCells
  
  
  subroutine store_and_reload_data(dumcoordntes_xy,duml0,dumA0)
    implicit none
    real*8, intent(inout) :: dumcoordntes_xy(1:N_mvCoordnte)
    real*8 :: coordntes_xyStore(1:N_mvCoordnte)
    real*8, intent(inout) :: duml0(1:N_spr)
    real*8 :: l0Store(1:N_spr)
    real*8, intent(inout) :: dumA0(1:N_cell)
    real*8 :: A0Store(1:N_cell)
    
    call store_data
    call reload_data
    
    
  contains
    
    subroutine store_data
      implicit none
      
      coordntes_xyStore = dumcoordntes_xy
      l0Store = duml0
      A0Store = dumA0
      
    end subroutine store_data
    
    
    subroutine reload_data
      
      dumcoordntes_xy = coordntes_xyStore
      duml0 = l0Store
      dumA0 = A0Store
      
    end subroutine reload_data
    
  end subroutine store_and_reload_data
  
  
  subroutine get_CgXNode_withStep(step_no)
    implicit none
    integer, intent(in) :: step_no
    integer :: i,j,jmax
    integer :: cnt_strt,cnt_node
    
    do i = 1,2
       
       if (i==1) then
          jmax=nvsl
          cnt_strt = 0
       elseif (i==2) then
          jmax=nvsr
          cnt_strt = 2*nvsl
       endif
       
       !write(*,*) jmax,"jmax"
       
       do j = 1,jmax
          
          if (j.ne.1) then
             continue
          else
             cnt_node = cnt_strt + (2*j-1)
             
             if (i==1) then !set to pulling now, to push make it i=2,nxt i=1
                CgXNode(cnt_node)   = (step_no-1)*ppS
                CgXNode(cnt_node+1) = (step_no-1)*ppS
                
             elseif (i==2) then
                CgXNode(cnt_node)   = (step_no-1)*(-ppS)
                CgXNode(cnt_node+1) = (step_no-1)*(-ppS)
             endif
             
          endif
          
       enddo
       
    enddo
    
    open(unit=121,file="CgXNode_value.dat")
    
    do i = 1,N_node
       write(unit=121,fmt=*) "CgXNode=",CgXNode(i),"for i=",i
       !if (step_no==2) write(*,*) "CgXNode=",CgXNode(i),"for Node=",i
    enddo
    
    close(121)
    !if (step_no==2) stop
    
  end subroutine get_CgXNode_withStep
  
  
  subroutine increase_bottomSurface_CgXNode(step_no)
    implicit none
    integer, intent(in) :: step_no
    integer :: i,j,jmax
    integer :: cnt_strt,cnt_node
    
    do i = 1,2
       
       if (i==1) then
          jmax=nvsl
          cnt_strt = 0
       elseif (i==2) then
          jmax=nvsr
          cnt_strt = 2*nvsl
       endif
       
       
       do j = 1,jmax
          
          if (j.ne.1) then
             continue
          else
             cnt_node = cnt_strt + (2*j-1)
             
             if (i==1) then
                CgXNode(cnt_node+1) = (step_no-1)*(0.020d0)
                
             elseif (i==2) then
                CgXNode(cnt_node+1) = (step_no-1)*(-0.020d0)
                
             endif
             
          endif
          
       enddo
       
    enddo
    
    do i = 1,N_node
       !write(unit=121,fmt=*) "CgXNode=",CgXNode(i),"for i=",i
       if (step_no==10) write(*,*) "CgXNode=",CgXNode(i),"for Node=",i
    enddo
    
    !if (step_no==10) stop
    
  end subroutine increase_bottomSurface_CgXNode
  
  
  
  subroutine updt_kspr_with_l0
    implicit none
    integer :: i
    
    open(unit=90,file='updt_Prop.dat',position='append')
    
    do i = 1,N_spr
       write(90,fmt=*) k_spr(i),l0(i),i,"k_spr,l0,spr_no (bfr)"
       k_spr(i) = 1.00d0/l0(i)
       write(90,fmt=*) k_spr(i),i,"k_spr,spr_no (aft)"
    enddo
    
    write(90,fmt=*) " "  
    close(90)
    
  end subroutine updt_kspr_with_l0
  
  subroutine updt_karea_with_A0
    implicit none
    integer :: i
    
    open(unit=90,file='updt_Prop.dat',position='append')
    
    do i = 1,N_cell
       write(90,fmt=*) k_area(i),A0(i),i,"k_area,A0,cell_no (bfr)"
       k_area(i) = 1.00d0/A0(i)
       write(90,fmt=*) k_area(i),i,"k_area,cell_no (aft)"
    enddo
    
    write(90,fmt=*) " "
    close(90)
    
  end subroutine updt_karea_with_A0
  
  
  subroutine get_thrshold_angleInfo(dum_coordntesxy,thrshAnglInfo,thrshAngl)
    implicit none
    real*8, intent(in)  :: dum_coordntesxy(1:N_mvCoordnte)
    integer,intent(out) :: thrshAnglInfo(1:N_thrshAngl,1:2)
    real*8, intent(out) :: thrshAngl(1:N_thrshAngl)
    
    real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer :: i,j,jstrt,jfnsh
    integer :: cnt
    integer :: node_nm,area_nm,area_serial
    
    real*8  :: CN(1:2),NN1(1:2),NN2(1:2)
    real*8  :: Phi
    
    call coordntes_to_nodes(dum_coordntesxy,dum_Nodes)
    
    open(unit=176,file='thrshAnglInfo.dat')
    
    cnt = 1
    
    do i = 1,2
       if (i==1) then
          jstrt = 1
          jfnsh = ncl
       elseif (i==2) then
          jstrt = ncl+1
          jfnsh = ncl+ncr
       endif
       
       do j = jstrt,jfnsh
          area_nm = j
          node_nm = area_node(area_nm,1)
          
          thrshAnglInfo(cnt,1) = node_nm
          thrshAnglInfo(cnt,2) = area_nm
          
          write(176,fmt=*) thrshAnglInfo(cnt,1:2),cnt,"thrshAngleInfo,cnt"
          
          call get_area_serial(node_nm,area_nm,area_serial)
          call get_CN_NN1_NN2_nodes(node_nm,area_serial,dum_Nodes,CN,NN1,NN2)
          call get_Phi(CN,NN1,NN2,Phi)
          
          thrshAngl(cnt) = Phi
          
          if(i==2 .and.j==jfnsh) then
             continue
          else
             cnt = cnt+1
          endif
          
       enddo
       
    enddo
    
    write(176,fmt=*) cnt,N_thrshAngl,"cnt,N_thrshAngl"
    close(176)
    
    if (cnt.ne.N_thrshAngl) then
       write(*,*) "thrshAngle count is not matching"
       write(*,*) "sb:get_thrshold_angles,fl:Manipualting_prop"
       stop
    endif
    
  end subroutine get_thrshold_angleInfo
  
  
  subroutine get_curr_thrshldAngls(dum_coordntesxy,thrshAnglInfo,curr_thrshAngl)
    implicit none
    real*8, intent(in)  :: dum_coordntesxy(1:N_mvCoordnte)
    integer,intent(in)  :: thrshAnglInfo(1:N_thrshAngl,1:2)
    real*8, intent(out) :: curr_thrshAngl(1:N_thrshAngl)
    
    real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
    integer :: i,j,jstrt,jfnsh
    integer :: cnt
    integer :: node_nm,area_nm,area_serial
    
    real*8  :: CN(1:2),NN1(1:2),NN2(1:2)
    real*8  :: Phi
    
    open(unit=190,file='ThrshAngle.dat')
    call coordntes_to_nodes(dum_coordntesxy,dum_Nodes)
    
    cnt = 1
    
    do i = 1,2
       if (i==1) then
          jstrt = 1
          jfnsh = ncl
       elseif (i==2) then
          jstrt = ncl+1
          jfnsh = ncl+ncr
       endif
       
       do j = jstrt,jfnsh
          node_nm = thrshAnglInfo(cnt,1)
          area_nm = thrshAnglInfo(cnt,2)
          
          call get_area_serial(node_nm,area_nm,area_serial)
          write(190,fmt=*) node_nm,area_nm,area_serial,"node_nm,area_nm,area_serial"
          call get_CN_NN1_NN2_nodes(node_nm,area_serial,dum_Nodes,CN,NN1,NN2)
          call get_Phi(CN,NN1,NN2,Phi)
          
          curr_thrshAngl(cnt) = Phi
          
          cnt = cnt+1
       enddo
          
    enddo
    
    close(190)
    
  end subroutine get_curr_thrshldAngls
  
  subroutine cmpr_thrshAngls(curr_thrshAngl,thrshAngl,thrsh_lgcl)
    implicit none
    real*8, intent(in)  :: curr_thrshAngl(1:N_thrshAngl)
    real*8, intent(in)  :: thrshAngl(1:N_thrshAngl)
    logical,intent(out) :: thrsh_lgcl
    
    integer :: i,j,jstrt,jfnsh
    
    open(unit=146,file="cmpr_thrshAngls.dat")
    write(146,fmt=*) "curr_thrshAngl,preset_thrshAngl,Angle No,logical"
    
    thrsh_lgcl = .False.
    
    do i = 1,N_thrshAngl
       !if (i.le.4) then
          if (curr_thrshAngl(i) .lt. thrshAngl(i)) then
             thrsh_lgcl = .True.
             write(146,fmt=*) curr_thrshAngl(i),thrshAngl(i),i,thrsh_lgcl
             exit
          endif
          
          write(146,fmt=*) curr_thrshAngl(i),thrshAngl(i),i,thrsh_lgcl
       !else
        !  continue
       !endif
    enddo
    
  end subroutine cmpr_thrshAngls
  
  
  subroutine symmetricize_prop
    implicit none
    integer :: i
    integer :: Hlf_Nspr
    
    integer,allocatable :: lftSpr(:),rghtSpr(:)
    
    open(unit=198,file='symmetricCell.dat',position='append')
       
    Hlf_Nspr = N_spr/2
    allocate(lftSpr(1:Hlf_Nspr),rghtSpr(1:Hlf_Nspr))
    
    do i = 1,Hlf_Ncell
       write(198,fmt=*) lftSideCell(i),rghtSideCell(i),i
       k_area(lftSideCell(i)) = k_area(rghtSideCell(i))
       A0(lftSideCell(i))     = A0((rghtSideCell(i)))
       write(198,fmt=*) A0(lftSideCell(i)),A0(rghtSideCell(i)),k_area(lftSideCell(i)),k_area(rghtSideCell(i))
    enddo
    
    call lftSide_and_rghtSide_SprPairs(lftSpr,rghtSpr,Hlf_Nspr)
    
    do i = 1,Hlf_Nspr
       write(198,fmt=*) lftSpr(i),rghtSpr(i),i
       k_spr(lftSpr(i)) = k_spr(rghtSpr(i))
       l0(lftSpr(i))    = l0(rghtSpr(i))
       write(198,fmt=*) l0(lftSpr(i)),l0(rghtSpr(i)),k_spr(lftSpr(i)),k_spr(rghtSpr(i))
    enddo
    
    close(198)
    
  end subroutine symmetricize_prop
  
  subroutine lftSide_and_rghtSide_SprPairs(lftSpr,rghtSpr,Hlf_Nspr)
    implicit none
    integer, intent(in)  :: Hlf_Nspr
    integer, intent(out) :: lftSpr(1:Hlf_Nspr),rghtSpr(1:Hlf_Nspr)
    
    integer :: i,cnt
    integer :: lftCell,rghtCell
    integer :: sprL1,sprL2,sprL3
    integer :: sprR1,sprR2,sprR3
    
    cnt = 1
    
    do i = 1,Hlf_Ncell
       lftCell  = lftSideCell(i)
       rghtCell = rghtSideCell(i)
       
       sprL1 = 3*lftCell-2 ; sprR1 = 3*rghtCell-2
       sprL2 = 3*lftCell-1 ; sprR2 = 3*rghtCell-1
       sprL3 = 3*lftCell   ; sprR3 = 3*rghtCell
       
       lftSpr(cnt)   = sprL1 ; rghtSpr(cnt)   = sprR1  
       lftSpr(cnt+1) = sprL2 ; rghtSpr(cnt+1) = sprR2
       lftSpr(cnt+2) = sprL3 ; rghtSpr(cnt+2) = sprR3
       
       cnt = cnt+3
       
    enddo
    
    if ((cnt-1).ne.Hlf_Nspr) then
       stop
       write(*,*) "cnt is not equal to Half_Nspr"
    endif
    
  end subroutine lftSide_and_rghtSide_SprPairs
  
  subroutine symmetricize_nodes
    implicit none
    integer :: Hlf_Nnode
    integer :: i
    integer, allocatable :: lftSNode(:),rghtSNode(:)
    
    Hlf_Nnode = N_node/2
    
  end subroutine symmetricize_nodes  
  
  ! tbd = TO BE DETERMINED
  
  subroutine split_spr_srialy_with_ratio_tbd(SprNm,Choice,ks_splt,l0_splt,l_splt) 
    implicit none
    integer, intent(in)  :: SprNm,Choice
    real*8 , intent(out) :: ks_splt(1:2),l0_splt(1:2),l_splt(1:2)
    
    real*8  :: nodeV1(1:N_dmnsn),nodeV2(1:N_dmnsn),nodeV3(1:N_dmnsn)
    integer :: node1,node2,node3
    integer :: no_of_split
    real*8  :: rat_mltp = 1000.00d0
    real*8  :: len1,len2
    real*8  :: ratio_of_len
    real*8  :: ratio1,ratio2
    integer :: N_sprV,N_nodeV
    
    real*8,  allocatable :: ksV(:),l0V(:),lV(:)
    real*8,  allocatable :: nxyV(:,:)
    
    interface
       real*8 function distance(Point1,Point2,N_dim)
         implicit none
         real*8  :: Point1(1:N_dim),Point2(1:N_dim)
         integer :: N_dim
       end function distance
    end interface
    
    if (Choice == 0) then
       continue
    elseif (Choice == 1) then
       continue
    elseif (Choice==2) then
       
       N_nodeV = N_nodeSES
       allocate(nxyV(1:N_nodeV,1:N_dmnsn))
       nxyV    = -1.0d20
       
       N_sprV  = N_sprSES
       allocate(ksV(1:N_sprV),l0V(1:N_sprV),lV(1:N_sprV))
       ksV=-1.0d20 ; l0V=-1.0d20 ; lV=-1.0d20
       
       ksV(1:N_sprV) = k_sprStrES(1:N_sprV)
       l0V(1:N_sprV) = l0_StrES(1:N_sprV)
       lV(1:N_sprV)  = l_StrES(1:N_sprV)
       
       nxyV(1:N_nodeV,1:N_dmnsn) = node_xyStrES(1:N_nodeV,1:N_dmnsn)
       
    endif
    
    if (spr_node(SprNm,0) == 2) write(*,*) "Ratio needs to as an Input"
    if (spr_node(SprNm,0) == 2) stop
    
    no_of_split = spr_node(SprNm,0)-1
    
    if (no_of_split .gt. 2) write(*,*) "not ready for generalized splitting"
    if (no_of_split .gt. 2) stop
    
    node1 = spr_node(SprNm,1)
    node2 = spr_node(SprNm,2)
    node3 = spr_node(SprNm,3)
    
    nodeV1(1:N_dmnsn) = nxyV(node1,1:N_dmnsn)
    nodeV2(1:N_dmnsn) = nxyV(node2,1:N_dmnsn)
    nodeV3(1:N_dmnsn) = nxyV(node3,1:N_dmnsn)

    write(*,*) nodeV1(1:2),nodeV2(1:2),nodeV3(1:2),"nodeVs"
    
    len1 = distance(nodeV1,nodeV2,N_dmnsn)
    len2 = distance(nodeV2,nodeV3,N_dmnsn)
    
    ratio_of_len = len1/len2
    ratio1       = ratio_of_len*rat_mltp
    ratio2       = rat_mltp

    write(*,*) ksV(SprNm),l0V(SprNm),lV(SprNm),"ks-l0-l bfr_splt"
    write(*,*) len1,len2,ratio_of_len,ratio1,ratio2,"prp val is spr splt"
    
    ks_splt(1) = ((ratio1+ratio2)/ratio1)*(ksV(SprNm))
    ks_splt(2) = ((ratio1+ratio2)/ratio2)*(ksV(SprNm))
    
    l0_splt(1) = (ratio1/(ratio1+ratio2))*(l0V(SprNm)) 
    l0_splt(2) = (ratio2/(ratio1+ratio2))*(l0V(SprNm))
    
    l_splt(1) = (ratio1/(ratio1+ratio2))*(lV(SprNm)) 
    l_splt(2) = (ratio2/(ratio1+ratio2))*(lV(SprNm))
    
    write(*,*) ks_splt(1:2),l0_splt(1:2),l_splt(1:2),"l_aft_splt"
    
  end subroutine split_spr_srialy_with_ratio_tbd
  
  
  subroutine split_spr_srialy_equal_len(SprNm,Choice,ks_splt,l0_splt,l_splt) 
    implicit none
    integer, intent(in)  :: SprNm,Choice
    real*8 , intent(out) :: ks_splt(1:2),l0_splt(1:2),l_splt(1:2)
    
    real*8               :: ratio1,ratio2
    integer              :: N_sprV
    real*8,  allocatable :: ksV(:),l0V(:),lV(:)
    real*8,  allocatable :: nxyV(:,:)
    
    
    if (Choice == 0) then
       continue
    elseif (Choice == 1) then
       continue
    elseif (Choice==2) then
       
       N_sprV  = N_sprSES
       allocate(ksV(1:N_sprV),l0V(1:N_sprV),lV(1:N_sprV))
       ksV=-1.0d20 ; l0V=-1.0d20 ; lV=-1.0d20
       
       ksV(1:N_sprV) = k_sprStrES(1:N_sprV)
       l0V(1:N_sprV) = l0_StrES(1:N_sprV)
       lV(1:N_sprV)  = l_StrES(1:N_sprV)
       
    endif
    
    
    ratio1 = 1.00d0 ; ratio2 = 1.00d0
    
    ks_splt(1) = ((ratio1+ratio2)/ratio1)*(ksV(SprNm))
    ks_splt(2) = ((ratio1+ratio2)/ratio2)*(ksV(SprNm))
    
    write(*,*) ks_splt(1:2),"ks_aft_splt"
    
    l0_splt(1) = (ratio1/(ratio1+ratio2))*(l0V(SprNm)) 
    l0_splt(2) = (ratio2/(ratio1+ratio2))*(l0V(SprNm))
    
    write(*,*) l0_splt(1:2),"l0_aft_splt"
    
    l_splt(1) = (ratio1/(ratio1+ratio2))*(lV(SprNm)) 
    l_splt(2) = (ratio2/(ratio1+ratio2))*(lV(SprNm))
    
    write(*,*) l_splt(1:2),"l_aft_splt"
    
  end subroutine split_spr_srialy_equal_len

  
  subroutine cmbine_spr_prlelly_cnnctd(NsprToPr,ks_bp,l0_bp,l_bp,ks_cm,l0_cm,l_cm)
    implicit none
    integer, intent(in)  :: NsprToPr
    real*8 , intent(in)  :: ks_bp(1:NsprToPr),l0_bp(1:NsprToPr),l_bp(1:NsprToPr)
    real*8 , intent(out) :: ks_cm,l0_cm,l_cm
    
    integer :: i,j
    real*8  :: l0_sum
    
    ks_cm  = 0.0d0 ; l0_cm = 0.0d0 ; l_cm = 0.0d0
    l0_sum = 0.0d0
    
    do i = 1,NsprToPr
       ks_cm  = ks_cm  + ks_bp(i)
       l0_sum = l0_sum + l0_bp(i)
    enddo
    
    
    l0_cm = l0_sum/NsprToPr
    l_cm  = l_bp(1)
    
    write(*,*) ks_bp(1:2),"ks_bp"
    write(*,*) l0_bp(1:2),"l0_bp"
    write(*,*) l_bp(1:2),"l_bp"
    
  end subroutine cmbine_spr_prlelly_cnnctd
  
  
  subroutine cmbine_spr_serially_cnnctd(NsprToCm,ks_bc,l0_bc,l_bc,ks_ac,l0_ac,l_ac)
    implicit none
    integer, intent(in)  :: NsprToCm
    real*8 , intent(in)  :: ks_bc(1:NsprToCm),l0_bc(1:NsprToCm),l_bc(1:NsprToCm)! bc = bfr cmbining
    real*8 , intent(out) :: ks_ac,l0_ac,l_ac ! ac = aft cmbining
    
    integer :: i,j
    real*8  :: ks_acInv
    
    ks_acInv=0.0d0 ; l0_ac=0.0d0 ; l_ac=0.0d0
    
    do i = 1,NsprToCm
       ks_acInv = ks_acInv + (1.0d0/ks_bc(i))
       l0_ac    = l0_ac + l0_bc(i)
       l_ac     = l_ac  + l_bc(i)
    enddo
    
    ks_ac = 1.0d0/ks_acInv
    
    write(*,*) ks_ac,l0_ac,l_ac,"ks-l0-l (ac)"
    
  end subroutine cmbine_spr_serially_cnnctd
  
  
end module changing_parameters

