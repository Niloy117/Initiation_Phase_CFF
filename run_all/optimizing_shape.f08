module stiff_bottom4cell_holding_it_in_prev_pos
  use system_parameters
  use transfrm_info
  use PenF_module
  use PenF_minimizer_mod
  implicit none
  
  integer :: N_activeSpr,N_activePair
  integer,allocatable :: activeSpr(:)

  integer :: pairA_spr(1:2)
  integer :: pairB_spr(1:2)
  integer :: pairC_spr(1:2)
  integer :: pairD_spr(1:2)
  integer :: pairE_spr(1:2)
  integer :: pairF_spr(1:2)
  integer :: pairG_spr(1:2)
  
contains

  subroutine minimize_activePairSprs
    implicit none
    integer :: hm
    integer :: spr_or_area
    integer :: pairSpr(1:2)
    logical :: ext_flg

    real*8  :: PenF_prv,PenF_new,PenF_diff
    real*8  :: EPS

    integer :: i
    integer :: lp_cnt

    EPS = 1e-02
    lp_cnt = 1
    
    ext_flg=.False.
    hm = 2 
    spr_or_area = 1

    call get_activePairSprs

    write(*,*) A0(9:12),"A0 insd minimize_activePairSprs"
    
    PenF_prv = 1.0d0
    !PenF_prv = PenF(node_xy)
    write(*,*) PenF_prv,"PenF_prv"

    
    do while(.not. ext_flg)

       do i = 1,N_activePair
          if (i==1) pairSpr = pairA_spr
          if (i==2) pairSpr = pairB_spr
          if (i==3) pairSpr = pairC_spr
          if (i==4) pairSpr = pairE_spr
          if (i==5) pairSpr = pairD_spr
          if (i==6) pairSpr = pairF_spr
          if (i==7) pairSpr = pairG_spr

          !if (i.gt.5) continue 
          if (i==8) then
             write(*,*) "Exceed max pair,sbrtn: minimize_activePairSprs,flnm:optimizing_shape"
             stop
          endif


          call minimize_PenF_wrt_l0s(hm,pairSpr,spr_or_area)
          !if (i.ge.5) continue
       enddo

       
       PenF_new = PenF(node_xy)
       write(*,*) PenF_prv,PenF_new,"PenF_prv,PenF_new"
       
       PenF_diff = PenF_prv-PenF_new
       
       lp_cnt = lp_cnt + 1
       write(*,*) lp_cnt,"lp_cnt"
       
       if (PenF_diff .lt. 0) then
          write(*,*) "Bug in minimization, sbrtn: minimize_activePairSprs,flnm:optimizing_shape.f08"
          !PenF_prv = PenF_prv
          !stop
       elseif (PenF_diff .lt. EPS) then
          ext_flg = .True.
       else
          PenF_prv = PenF_new
       endif

    enddo


  end subroutine minimize_activePairSprs
  
  

  
  subroutine get_activePairSprs
    implicit none

    call get_N_activespr
    call pairwise_list_of_active_spr_bot4_cell

    
  end subroutine get_activePairSprs


  
  subroutine get_N_activespr
    implicit none

    N_activeSpr = nsecclb*2 + nseccrb*2 + 2 

    if(mod(N_activeSpr,2) .ne. 0) then
       write(*,*) "N_activeSpr cant be odd here,sbrtn:get_active_spr,flnm:optimizing_shape"
       stop

    else
       N_activePair = N_activeSpr/2
    endif

    
  end subroutine get_N_activespr

  subroutine pairwise_list_of_active_spr_bot4_cell
    implicit none
    integer :: pair_spr(1:2)
    integer :: area1,area2
    integer :: half_pair,type
    
    integer :: i
    
    pair_spr = 0

    half_pair = N_activePair/2
    
    do i = 1,N_activePair
       
       if (i .le. half_pair) then
          area1 = N_cell-1 ; area2 = N_cell

          if (i==1) type = 4
          if (i==2) type = 5
          if (i==3) type = 3

       elseif (i .gt. half_pair) then
          
          if ( i .ne. N_activePair) then
             area1 = N_cell-3 ; area2 = N_cell - 2
          else
             area1 = N_cell-5 ; area2 = N_cell - 4
          endif

          if ((i-half_pair) == 1) type = 4
          if ((i-half_pair) == 2) type = 5 
          if ((i-half_pair) == 3) type = 3
          if ((i-half_pair) == 4) type = 5
          
       endif

       call get_SprPair(area1,area2,type,pair_spr)
       
       if (i==1) pairA_spr = pair_spr
       if (i==2) pairB_spr = pair_spr
       if (i==3) pairC_spr = pair_spr
       if (i==4) pairD_spr = pair_spr
       if (i==5) pairE_spr = pair_spr
       if (i==6) pairF_spr = pair_spr
       if (i==7) pairG_spr = pair_spr

       if (i.gt.7) then
          write(*,*) "Exceded the expected pair, sbrtn: pairwise_list,flnm: optimizing_shape"
          stop
       endif
       
       
    enddo

    write(*,*) pairA_spr,"A 32,35"
    write(*,*) pairB_spr,"B 33,36"
    write(*,*) pairC_spr,"C 31,34" 
    write(*,*) pairD_spr,"D 26,29"
    write(*,*) pairE_spr,"E 27,30"
    write(*,*) pairF_spr,"F 25,28"
    write(*,*) pairG_spr,"G 21,24" 
    
    
  end subroutine pairwise_list_of_active_spr_bot4_cell

  
  subroutine get_SprPair(area1,area2,type,spr_pair)
    integer :: area1,area2
    integer :: type
    integer :: spr_pair(1:2)

    integer :: i
    integer :: imax1,imax2,imax
    integer :: spr_nm
    
    imax1 = area_spr(area1,0)
    imax2 = area_spr(area2,0)
    
    imax = imax1 + imax2
    !write(*,*) i
    write(*,*) area1,area2,"area1,area2"
    do i = 1,imax
       if (i.le.imax1) then
          spr_nm = area_spr(area1,i)
          
          if (typ_spr(spr_nm).eq. type) then
             spr_pair(1) = spr_nm
          endif
          
       else
          spr_nm = area_spr(area2,i-imax1)
          
          if (typ_spr(spr_nm) .eq. type) then
             spr_pair(2) = spr_nm
          endif
          
       endif
    enddo

  end subroutine get_SprPair

  
end module stiff_bottom4cell_holding_it_in_prev_pos


