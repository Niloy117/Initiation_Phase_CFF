module generating_the_shape
  
  use sys_building_info
  use system_parameters
  use diff_stage
  use strct_variables
  use node_variables
  use spring_variables
  use area_variables
  use grvtnl_variables
  use bend_variables
  use SRyp_variables
  use curve_variables
  use end_info
  use moving_coordnte_variables
  
  implicit none
  !integer, allocatable :: dependency_array(:,:)
  
contains
  
  subroutine generate_blck_params
    implicit none
    integer :: i
    integer :: lp_end
    
    lp_end = wntbb(0)
    
    !write(*,*) lp_end,"lp_end"
    
    !call alloc_init_and_build_dependency_array (not needed right now, will see if needed later)

    
    call get_AspectRatio_LaLb
    
    do i = 1,max_blck
       !write(*,*) wntbb(i),"wntbb(i)"
       if (wntbb(i) .eq. 0) then
          continue
       else
          call generate_indvd_blck_params(wntbb(i))
       endif
       
    enddo
    
    !now will come the test dependent routines
    
    !call lft_rght_common_params
    call global_params
    call wall_cell_lngth_and_other_params
    
  end subroutine generate_blck_params


  
  subroutine generate_indvd_blck_params(blck_id)
    implicit none
    integer :: blck_id
    
    if (blck_id == 1) then
       call lft_ladder_params
    elseif (blck_id == 2) then
       call rght_ladder_params
    elseif (blck_id == 3) then
       call clft_top_params
    elseif (blck_id == 4) then
       call additnl_node_params
    elseif (blck_id == 5) then
       call crght_top_params
    elseif (blck_id == 6) then
       call clft_bot_params
    elseif (blck_id == 7) then
       call crght_bot_params
    elseif (blck_id == 8) then
       call cntrlCell_params
    elseif (blck_id == 9) then
       call endCell_params
    endif
    
  end subroutine generate_indvd_blck_params
  
  
  
  
  subroutine generate_blcks
    implicit none
    integer :: i
    integer :: lp_end
   
    !integer :: depending_blck(1:max_blck)
    
   
    call allocate_and_initialize_node_variables
    call allocate_and_initialize_spring_variables
    call allocate_and_initialize_area_variables
    call allocate_and_initialize_grvVars
    call allocate_and_initialize_bend_variables
    
    call allocate_and_initialize_end_node_variables
    call allocate_and_initialize_spring_end_variables
    call allocate_and_initialize_area_end_variables
    
    lp_end = wntbb(0)
    
    do i = 1,max_blck
       !call get_depending_info(wntbb(i),depending_blck)
       if (wntbb(i) == 0) then
          continue
       else
          call generate_indvd_blck(wntbb(i))
          whabb(i) = wntbb(i)
       endif
       call update_lgicls_of_whabb_array      
    enddo
    
    if (stageNo==1 .or. stageNo==2) then
       call get_InitiatorCell
       call shifting_origin_at_InitiatorCell
       
    elseif (stageNo==3) then
       write(*,*) "BRING MY ATTENTION in generating_sys"
       stop
    elseif (stageNo==4) then
       call get_InitiatorCell
       call shifting_origin_at_Pulley_Node
    endif
    
    call get_N_optmSpr_and_optmSpr
    call get_N_optmCell_and_optmCell
        
    call get_SideWiseCellNo
    call get_select_xy
    
    call print_NodeVars
    call print_spring_variables
    call print_area_variables
    
  end subroutine generate_blcks

  
  subroutine generate_indvd_blck(blck_id)
    implicit none
    integer :: blck_id!,prev_blck_id
    integer :: cnnctng_node
    !integer :: depending_blck(1:max_blck)
    
    if (blck_id == 1) then
       call lft_Nodes_Types_and_Ends
       call lft_spring_Properties
       call lft_area_Properties
       !call lft_blck_end
       
    elseif (blck_id == 2) then
       call rght_Nodes_Types_and_Ends
       call rght_spring_Properties
       call rght_area_Properties
       !call rght_blck_end
       
    elseif (blck_id == 3) then
       call clft_top_Nodes_Types_and_Ends
       call clft_top_spring_Properties
       call clft_top_area_Properties
       
    elseif (blck_id == 4) then
       !call get_cnnctng_node(prev_blck_id,cnnctng_node)
       !call additnl_Nodes_Types_and_Ends(cnnctng_node)
       !call additnl_spring_Properties
       continue
       
    elseif (blck_id == 5) then
       call crght_top_Nodes_Types_and_Ends
       call crght_top_spring_Properties
       call crght_top_area_Properties
       
    elseif (blck_id == 6) then
       call clft_bot_Nodes_Types_and_Ends
       call clft_bot_spring_Properties
       call clft_bot_area_Properties
       
    elseif (blck_id == 7) then
       call crght_bot_Nodes_Types_and_Ends
       call crght_bot_spring_Properties
       call crght_bot_area_Properties
       
    elseif (blck_id == 8) then
       call cntrlCell_Nodes_Types_and_Ends
       call cntrlCell_spring_Properties
       call cntrlCell_area_Properties

    elseif (blck_id == 9) then
       call endRgnCell_Nodes_Types_and_Ends
       call endRgnCell_spring_Properties
       call endRgnCell_area_Properties
    else
       continue
    endif

  end subroutine generate_indvd_blck


  subroutine destroy_blcks
    implicit none
    
    call deallocate_node_variables
    call deallocate_spring_variables
    call deallocate_area_variables
    call deallocate_grvVars
    call deallocate_bend_variables
    call deallocate_curve_variables
    
    call deallocate_end_node_variables
    call deallocate_spring_end_variables
    call deallocate_area_end_variables
    
  end subroutine destroy_blcks


  
  subroutine get_depending_info(blck_id,depending_blck)
    implicit none
    
    integer :: blck_id
    integer :: depending_blck(1:max_blck)
    
    integer :: lp_end
    integer :: i
    
    depending_blck = -1
    
    if (dependency_array(blck_id,0)==0 .or. dependency_array(blck_id,0)==(-1)) then
       continue
    else
       
       lp_end = dependency_array(blck_id,0)
       
       do i = 1,lp_end
          depending_blck(i) = dependency_array(blck_id,i)
          !call get_blck_ends(depending_blck)
       enddo
    endif
    
  end subroutine get_depending_info
  
  subroutine get_cnnctng_node(prev_blck_id,cnnctng_node)
    implicit none
    integer, intent(in)  :: prev_blck_id
    integer, intent(out) :: cnnctng_node

    if (prev_blck_id == 1) then
       cnnctng_node = lft_endNode(2)
    elseif (prev_blck_id == 2) then
       cnnctng_node = rght_endNode(2)
    elseif (prev_blck_id == 3) then
       cnnctng_node = clft_top_endNode(2)
    endif !!if needed will add more prev_blck_id
    
  end subroutine get_cnnctng_node

  subroutine get_cnnctng_area(prev_blck_id,cnnctng_area)
    implicit none

    integer, intent(in)  :: prev_blck_id
    integer, intent(out) :: cnnctng_area

    if (prev_blck_id == 1) then
       cnnctng_area = 1
    elseif (prev_blck_id == 2) then
       cnnctng_area = 2
    endif
    
  end subroutine get_cnnctng_area
  
  subroutine update_lgicls_of_whabb_array
    implicit none
       
    if (whabb(1) .ne. 0) then
       lft_alrdy_built          = .True.
    elseif (whabb(2) .ne. 0) then
       rght_alrdy_built         = .True.
    elseif (whabb(3) .ne. 0) then
       clft_top_alrdy_built    = .True.
    elseif (whabb(4) .ne. 0) then
       additnl_node_alrdy_built = .True.
    elseif (whabb(5) .ne. 0) then
       crght_top_alrdy_built   = .True.
    elseif (whabb(6) .ne. 0) then
       clft_bot_alrdy_built = .True.
    elseif (whabb(7) .ne. 0) then
       crght_bot_alrdy_built   = .True.  
    endif
    
  end subroutine update_lgicls_of_whabb_array









  subroutine alloc_init_and_build_dependency_array
    implicit none
    
    allocate(dependency_array(1:max_blck,0:max_blck))
    
    dependency_array = -1
    
    call get_the_dependency_array
    
  contains

    subroutine get_the_dependency_array
      implicit none
      integer :: i
      
      do i = 1,max_blck
         if (i==1) then !lft_blck
            dependency_array(i,0) = 0 !dependcy_array(i,0) tells on how many blcks it depends on 
         elseif (i==2) then !rght_blck
            dependency_array(i,0) = 0
         elseif (i==3) then !
            dependency_array(i,0) = 1
            dependency_array(i,1) = 1
         elseif (i==4) then
            dependency_array(i,0) = 2
            dependency_array(i,1) = 3
            dependency_array(i,2) = 4
         elseif (i==5) then
            dependency_array(i,0) = 2
            dependency_array(i,1) = 2
            dependency_array(i,2) = 3
         else
            continue
         endif
      enddo
      
    end subroutine get_the_dependency_array
    
  end subroutine alloc_init_and_build_dependency_array






  
end module generating_the_shape



module perturbing_the_shape
  !will try to write more generalized way of perturbing later, after having a good idea of how the testing will be done
  use generating_the_shape !maintain hierarchy(after generating the shape, I will use perturbing_the_shape)
  implicit none
  
contains
  
  subroutine change_coordnte_of_a_node(node_nmbr,dx,dy)
    implicit none
    integer :: node_nmbr,other_node !if any
    real*8  :: dx
    real*8  :: dy
    
    node_xy(node_nmbr,1) = node_xy(node_nmbr,1) + dx
    node_xy(node_nmbr,2) = node_xy(node_nmbr,2) + dy

    if(node_cnnctd(node_nmbr).ne.0 .and. count_this_dn(node_nmbr)==1) then
       other_node = node_cnnctd(node_nmbr)
       node_xy(other_node,1:2) = node_xy(node_nmbr,1:2)
    endif
    
  end subroutine change_coordnte_of_a_node
  
  subroutine change_coordnte_of_multiple_nodes(N_nodestoBechanged,which_nodes,dx,dy)
    
    implicit none
    integer :: N_nodestoBechanged
    integer :: which_nodes(1:N_nodestoBechanged)
    real*8  :: dx(1:N_nodestoBechanged)
    real*8  :: dy(1:N_nodestoBechanged)
    
    integer :: i
    integer :: node_nmbr,other_node
    
    dx = 0.0d0
    dy = 0.0d0

    do i = 1,N_nodestoBechanged
       node_nmbr = which_nodes(i)

       node_xy(node_nmbr,1) = node_xy(node_nmbr,1) + dx(i)
       node_xy(node_nmbr,2) = node_xy(node_nmbr,2) + dy(i)

       if(node_cnnctd(node_nmbr).ne.0.and.count_this_dn(node_nmbr).ne.(-1)) then
          other_node = node_cnnctd(node_nmbr)
          node_xy(other_node,1:2) = node_xy(node_nmbr,1:2)
       endif
    
    enddo
    
  end subroutine change_coordnte_of_multiple_nodes

  subroutine change_coordnte_of_a_type_of_node(type_of_node,dx,dy)
    implicit none
    integer :: type_of_node
    integer :: i
    real*8  :: dx,dy
    
    do i = 1,N_node
       if (node_typ(i) == type_of_node) then
          node_xy(i,1) = node_xy(i,1) + dx
          node_xy(i,2) = node_xy(i,2) + dy
       else
          continue
       endif
    enddo
    
  end subroutine change_coordnte_of_a_type_of_node

end module perturbing_the_shape


