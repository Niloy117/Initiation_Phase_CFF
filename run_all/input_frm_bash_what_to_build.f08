module sys_building_info 
  implicit none
  
  integer :: stageNo,stageType
  
  integer :: max_blck=9,actv_blck! depends on stageNo,stageType
  
  integer,allocatable :: wntbb(:) !wntbb = what needs to be built
  
  logical :: lft_toBe_built, rght_toBe_built, clft_top_toBe_built
  logical :: additnl_node_toBe_built, crght_top_toBe_built
  logical :: clft_bot_toBe_built,crght_bot_toBe_built
  logical :: midlCell_tobebuilt,endCell_tobebuilt
  
  integer,allocatable :: whabb(:) 
  
  logical :: lft_alrdy_built, rght_alrdy_built, clft_top_alrdy_built
  logical :: additnl_node_alrdy_built, crght_top_alrdy_built
  logical :: clft_bot_alrdy_built,crght_bot_alrdy_built
  logical :: midlCell_alrdy_built,endCell_alrdy_built
  
  integer,allocatable :: dependency_array(:,:)
  
contains
  
  subroutine get_stage_info
    implicit none
    
    if (stageNo==1) then
       if (stageType==1) actv_blck = 3 !left, right and central
       if (stageType==2) actv_blck = 2 !left and right
       
    elseif (stageNo==2) then 
       if (stageType==1) continue
       if (stageType==2) continue
       
    elseif (stageNo==3) then 
       if (stageType==1) continue
       if (stageType==2) continue
       
    elseif (stageNo==4) then
       if (stageType==1) actv_blck = 7 !adding an end cell
       if (stageType==2) actv_blck = 6
    endif 

    call alloc_init_and_read_wntbb_array
    call alloc_and_init_whabb_array
    
  end subroutine get_stage_info
  
  subroutine alloc_init_and_read_wntbb_array
    implicit none
    integer :: i
    integer :: count_blcks
    
    allocate(wntbb(0:max_blck))
    
    wntbb = 0
    
    lft_toBe_built          = .False.
    rght_toBe_built         = .False.
    clft_top_toBe_built     = .False.
    additnl_node_toBe_built = .False.
    crght_top_toBe_built    = .False.
    
    clft_bot_toBe_built     = .False.
    crght_bot_toBe_built    = .False.
    midlCell_tobebuilt      = .False.
    endCell_tobebuilt       = .False.
    
    call get_wntbb
    
    write(*,*) wntbb,"wntbb"
    
    call get_wntbbLgcls
    
    count_blcks = 0
    
    do i = 1,max_blck
       if (wntbb(i) .ne. 0) then
          count_blcks = count_blcks + 1
       else
          continue
       endif
    enddo
    write(*,*) count_blcks,"count_blcks_after_loop"
    
    if (count_blcks .ne. wntbb(0)) then
       write(*,*) "Number of blocks to be built are not matching,"
       write(*,*) "flnm:input_frm_bash_what_to_build, sbrtn: alloc_init_and_read_wntbb_array"
       stop
    endif
    
  end subroutine alloc_init_and_read_wntbb_array
  
  subroutine get_wntbb
    implicit none

    if (stageNo==1) then
       
       if (stageType==1) then
          wntbb(0) = actv_blck
          wntbb(1) = 1
          wntbb(2) = 2
          wntbb(8) = 8 !middle cell
       elseif (stageType==2) then
          wntbb(0) = actv_blck
          wntbb(1) = 1
          wntbb(2) = 2
       endif

    elseif (stageNo==2) then
       write(*,*) "Fix stage 2(1L),flnm:input_frm_bash_what_to_build"
       stop      
       
    elseif (stageNo==3) then
       write(*,*) "Fix stage 3(2E),flnm:input_frm_bash_what_to_build"
       stop
       
    elseif (stageNo==4) then
       
       if (stageType==1) then
          wntbb(0) = actv_blck
          wntbb(1) = 1
          wntbb(2) = 2
          wntbb(3) = 3
          wntbb(5) = 5
          wntbb(6) = 6
          wntbb(7) = 7
          wntbb(9) = 9 !end cell         
       elseif (stageType==2) then
          wntbb(0) = actv_blck
          wntbb(1) = 1
          wntbb(2) = 2
          wntbb(3) = 3
          wntbb(5) = 5
          wntbb(6) = 6
          wntbb(7) = 7
       endif
       
    endif
    
  end subroutine get_wntbb

  subroutine get_wntbbLgcls
    implicit none
    integer :: i
    
    open(unit=11,file='check_whats_going_on.dat')
    write(11,fmt=*) "list of blocks that has to be built are :"
    
    do i = 1,max_blck
       
       if (wntbb(i) .ne. 0) then

          if (wntbb(i) == 1) then
             lft_toBe_built  = .True.
             write(11,fmt=*) "Left Block"
          elseif (wntbb(i) == 2) then
             rght_toBe_built = .True.
             write(11,fmt=*) "Right Block"
          elseif (wntbb(i) == 3) then
             clft_top_toBe_built    = .True.
             write(11,fmt=*) "Central Left Top Block"
          elseif (wntbb(i) == 4) then
             additnl_node_toBe_built = .True.
             write(11,fmt=*) "Additional Node Block"
          elseif (wntbb(i) == 5) then
             crght_top_toBe_built   = .True.
             write(11,fmt=*) "Central Right Top Block"
             
          elseif (wntbb(i) == 6) then
             clft_bot_toBe_built    = .True.
             write(11,fmt=*) "Central Left Bottom Block"
          elseif (wntbb(i) == 7) then
             crght_bot_toBe_built    = .True.
             write(11,fmt=*) "Central Right Bottom Block"
          elseif (wntbb(i) == 8) then
             midlCell_tobebuilt = .True.
             write(11,fmt=*) "Middle Cell Block"
          elseif (wntbb(i) == 9) then
             endCell_tobebuilt = .True.
             write(11,fmt=*) "End Cell Block"
             
          elseif (wntbb(i).gt.max_blck) then
             write(*,*) "MORE THAN MAX BLOCK,PROGRAM TERMINATED,"
             write(*,*) "flnm:input_frm_bash_what_to_build, sbrtn: alloc_init_and_read_wntbb_array"
             stop
          endif
       endif
       
    enddo
    
    close(11)
    
  end subroutine get_wntbbLgcls
  
  subroutine alloc_and_init_whabb_array
    implicit none
    
    allocate(whabb(1:max_blck))
    
    whabb = 0
    
    lft_alrdy_built          = .False.
    rght_alrdy_built         = .False.
    clft_top_alrdy_built     = .False.
    additnl_node_alrdy_built = .False.
    crght_top_alrdy_built    = .False.
    
    clft_bot_alrdy_built     = .False.
    crght_bot_alrdy_built    = .False.
    midlCell_alrdy_built     = .False.
    endCell_alrdy_built      = .False.
    
  end subroutine alloc_and_init_whabb_array
  
  
end module sys_building_info



