module moving_coordnte_info
  use system_parameters
  implicit none
  real*8 , allocatable :: coordntes(:),coordntesStr(:),coordntesStrNI(:)
  real*8 , allocatable :: coordntesStrAN(:),coordntesStrES(:),coordntesStrNVFR(:)
  
  real*8 , allocatable :: coordntes_xy(:),coordntes_xyStr(:),coordntes_xyStrNI(:)
  real*8 , allocatable :: coordntes_xyStrAN(:),coordntes_xyStrES(:),coordntes_xyStrNVFR(:)
  
  integer, allocatable :: coordnte_of_which_node(:),coordnte_of_which_nodeStr(:),coordnte_of_which_nodeStrNI(:)
  integer, allocatable :: coordnte_of_which_nodeStrAN(:),coordnte_of_which_nodeStrES(:),coordnte_of_which_nodeStrNVFR(:)
  
  integer, allocatable :: x_or_y(:),x_or_yStr(:),x_or_yStrNI(:),x_or_yStrAN(:),x_or_yStrES(:),x_or_yStrNVFR(:)
  
end module moving_coordnte_info

module moving_coordnte_variables
  use cell_info
  use system_parameters
  use moving_coordnte_info
  
contains
  
  subroutine get_all_moving_coordnte_variables
    implicit none
    
    call get_N_mvCoordnte
    call get_N_mvCoordnte_withl0A0
    write(*,*) N_mvCoordnte,N_mvCoordnte_withl0A0,"N_mv and with l0A0" ! ; stop 'aft N_mv'
    call alloc_and_init_moving_coordnte_variables
    call moving_coordnte_values
    
  end subroutine get_all_moving_coordnte_variables
  
  subroutine get_all_moving_coordnte_variables_wo_StrVars
    implicit none
    
    call get_N_mvCoordnte
    call get_N_mvCoordnte_withl0A0
    write(*,*) N_mvCoordnte,N_mvCoordnte_withl0A0,"N_mv and with l0A0"
    !stop
    call alloc_and_init_moving_coordnte_variables_wo_StrVars
    call moving_coordnte_values
    
  end subroutine get_all_moving_coordnte_variables_wo_StrVars
  
  subroutine get_N_mvCoordnte
    implicit none
    integer :: i
    integer :: count_coordnte
    
    count_coordnte = 0
    
    do i = 1,N_node
       
       if (node_typ(i) .eq. 0) then !fixed node
          continue
          
       elseif (node_typ(i) .ne. 0) then
          
          if (node_cnnctd(i) .eq. 0) then
             if (node_typ(i) .eq. 1) then !free node
                count_coordnte = count_coordnte + 2
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                count_coordnte = count_coordnte + 1
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             if (count_this_dn(i) .eq. 1) then
                if (node_typ(i) .eq. 1) then !free node
                   count_coordnte = count_coordnte + 2
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   count_coordnte = count_coordnte + 1
                endif
             elseif (count_this_dn(i) .ne. 1) then
                continue
             endif
             
          endif
       endif
       
    enddo
    
    N_mvCoordnte = count_coordnte
    write(*,*) N_mvCoordnte,"N_mvCoordnte"
    
  end subroutine get_N_mvCoordnte
  
  
  subroutine get_N_mvCoordnte_withl0A0
    implicit none
    
    call get_N_mvCoordnte
    
    N_optmSpr  = 0
    N_optmCell = 0
    
    !write(*,*) l0_variatn_lgcl,A0_variatn_lgcl,"l0_variatn_lgcl,A0_variatn_lgcl"
    !write(*,*) N_mvCoordnte,"N_mvCoordnte"
    !write(*,*) N_mvCoordnte_withl0A0,"N_mvCoordnte_withl0A0"
    
    if((l0_variatn_lgcl.eqv..True.).AND.(A0_variatn_lgcl.eqv..False.)) then
       
       call get_optmSpr_count(N_optmSpr)
       N_mvCoordnte_withl0A0=N_mvCoordnte+N_optmSpr
       
    elseif((l0_variatn_lgcl.eqv..False.).AND.(A0_variatn_lgcl.eqv..True.)) then
       
       call get_optmCell_count(N_optmCell)
       N_mvCoordnte_withl0A0=N_mvCoordnte+N_optmCell
       
    elseif((l0_variatn_lgcl.eqv..True.).AND.(A0_variatn_lgcl.eqv..True.)) then
       
       call get_optmSpr_count(N_optmSpr)
       call get_optmCell_count(N_optmCell)
       N_mvCoordnte_withl0A0=N_mvCoordnte+N_optmSpr+N_optmCell
       
    elseif((l0_variatn_lgcl.eqv..False.).AND.(A0_variatn_lgcl.eqv..False.)) then
       !write(*,*) "Entering 4"
       write(*,*) "Both spring and area logical are False"
       N_mvCoordnte_withl0A0=N_mvCoordnte
    endif
    
    !write(*,*) N_mvCoordnte,"N_mvCoordnte"
    !write(*,*) N_mvCoordnte_withl0A0,"N_mvCoordnte_withl0A0"
    !stop
    
  end subroutine get_N_mvCoordnte_withl0A0
  
  subroutine get_optmSpr_count(lam_optmSpr)
    implicit none
    integer, intent(inout) :: lam_optmSpr
    integer :: i

    lam_optmSpr = 0
    
    do i = 1,N_spr
       if (optmSpr(i).ne.0 .AND. optmSpr(i).ne.1) then
          write(*,*) "optmSpr is neither 0 nor 1"
          stop
       else
          lam_optmSpr = lam_optmSpr + optmSpr(i)
       endif
    enddo
    
  end subroutine get_optmSpr_count

  subroutine get_optmCell_count(lam_optmCell)
    implicit none
    integer, intent(inout) :: lam_optmCell
    integer :: i

    lam_optmCell = 0
    
    do i = 1,N_cell
       if (optmCell(i).ne.0 .AND. optmCell(i).ne.1) then
          write(*,*) "optmCell is neither 0 nor 1"
          stop
       else
          lam_optmCell = lam_optmCell + optmCell(i)
       endif
    enddo
    
  end subroutine get_optmCell_count
  
  subroutine alloc_and_init_moving_coordnte_variables
    implicit none
    
    !write(*,*) coordntes
    
    allocate(coordntes(1:N_mvCoordnte_withl0A0))
    allocate(coordntesStr(1:N_mvCoordnte_withl0A0))
    !allocate(coordntesStrNVFR(1:N_mvCoordnte_withl0A0))
    
    allocate(coordntes_xy(1:N_mvCoordnte),coordntes_xyStr(1:N_mvCoordnte))
    !allocate(coordntes_xyStrNVFR(1:N_mvCoordnte))
    
    allocate(coordnte_of_which_node(1:N_mvCoordnte))
    allocate(coordnte_of_which_nodeStr(1:N_mvCoordnte))
    !allocate(coordnte_of_which_nodeStrNVFR(1:N_mvCoordnte))
    
    allocate(x_or_y(1:N_mvCoordnte),x_or_yStr(1:N_mvCoordnte))
    !allocate(x_or_yStrNVFR(1:N_mvCoordnte))
    
    coordntes              = -1.d30 ; coordntesStr              = -1.d30
    coordntes_xy           = -1.d30 ; coordntes_xyStr           = -1.d30
    coordnte_of_which_node = -1     ; coordnte_of_which_nodeStr = -1 
    x_or_y                 = -1     ; x_or_yStr                 = -1
    
    !coordntesStrNVFR              = -1.d30
    !coordntes_xyStrNVFR           = -1.d30
    !coordnte_of_which_nodeStrNVFR = -1 
    !x_or_yStrNVFR                 = -1
    
    
  end subroutine alloc_and_init_moving_coordnte_variables
  
  subroutine alloc_and_init_moving_coordnte_variables_wo_StrVars
    implicit none
    
    !write(*,*) coordntes
    
    allocate(coordntes(1:N_mvCoordnte_withl0A0))
    allocate(coordntes_xy(1:N_mvCoordnte))
    allocate(coordnte_of_which_node(1:N_mvCoordnte))
    allocate(x_or_y(1:N_mvCoordnte))
    
    coordntes              = -1.d30 
    coordntes_xy           = -1.d30 
    coordnte_of_which_node = -1
    x_or_y                 = -1
    
  end subroutine alloc_and_init_moving_coordnte_variables_wo_StrVars
  
  subroutine moving_coordnte_values
    implicit none
    integer :: i,j,jmax
    integer :: count_coordnte
    
    count_coordnte = 1
    
    do i = 1,N_node
       
       if (node_typ(i) .eq. 0) then
          continue
       elseif (node_typ(i) .ne. 0) then
          
          if (node_cnnctd(i) .eq. 0) then
             
             if (node_typ(i) .eq. 1) then !free node
                coordntes(count_coordnte:(count_coordnte+1)) = node_xy(i,1:N_dmnsn)
                coordnte_of_which_node(count_coordnte:(count_coordnte+1)) = i   
                x_or_y(count_coordnte)   = 1
                x_or_y(count_coordnte+1) = 2
                
                count_coordnte = count_coordnte + 2
                
             elseif (node_typ(i) .eq. 2) then !x free, y fixed
                coordntes(count_coordnte) = node_xy(i,1) !y fixed
                coordnte_of_which_node(count_coordnte) = i
                x_or_y(count_coordnte)   = 1
                
                count_coordnte = count_coordnte + 1
             endif
             
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) .eq. 1) then
                
                if (node_typ(i) .eq. 1) then !free node
                   
                   coordntes(count_coordnte:(count_coordnte+1)) = node_xy(i,1:N_dmnsn)
                   coordnte_of_which_node(count_coordnte:(count_coordnte+1)) = i
                   
                   x_or_y(count_coordnte)   = 1
                   x_or_y(count_coordnte+1) = 2
                   count_coordnte = count_coordnte + 2
                   
                elseif (node_typ(i) .eq. 2) then !x free, y fixed
                   coordntes(count_coordnte) = node_xy(i,1) !y fixed
                   coordnte_of_which_node(count_coordnte) = i
                   x_or_y(count_coordnte)   = 1
                   
                   count_coordnte = count_coordnte + 1
                endif
                
             elseif (count_this_dn(i) .ne. 1) then
                continue
             endif
             
          endif
          
       endif
       
    enddo
    
    !write(*,*) coordntes,"coordntes"
    !write(*,*) count_coordnte,"CC"
    
    count_coordnte = count_coordnte - 1
    
    if (count_coordnte.ne.N_mvCoordnte) then
       write(*,*) "flnm:conv_for_moving_sys,sb:moving_coor_value"
    endif
    
    write(*,*) N_mvCoordnte,N_mvCoordnte_withl0A0,"N_mv and with l0A0"
    
    do i=1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       if (i==1) then
          if(l0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       elseif (i==2) then
          if(A0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       endif
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (optmSpr(j)==1) then
                coordntes(count_coordnte+1) = l0(j)
                count_coordnte = count_coordnte + 1
             elseif (optmSpr(j)==0) then
                continue
             elseif (optmSpr(j).ne.0 .AND. optmSpr(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif
             
          elseif (i==2) then
             
             if (optmCell(j)==1) then
                coordntes(count_coordnte+1) = A0(j)
                count_coordnte = count_coordnte + 1
             elseif (optmCell(j)==0) then
                continue
             elseif (optmCell(j).ne.0 .AND. optmCell(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif
             
          endif
          
       enddo
    enddo
    
    if (count_coordnte.ne.N_mvCoordnte_withl0A0) then
       write(*,*) "count_coor =",count_coordnte,"N_mvCoordnte_withl0A0 =",N_mvCoordnte_withl0A0
       write(*,*) "count_coor not equal to N_mvCoordnte_withl0A0"
       write(*,*) "flnm:info_mod,sbrtn:nodesl0A0_to_coordntes"
       stop
    endif
    
    coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
    
  end subroutine moving_coordnte_values
  
  
  subroutine deallocate_moving_coordnte_variables
    implicit none
    
    deallocate(coordntes,coordntesStr)
    deallocate(coordntes_xy,coordntes_xyStr)
    deallocate(coordnte_of_which_node,coordnte_of_which_nodeStr)
    deallocate(x_or_y,x_or_yStr)
    
  end subroutine deallocate_moving_coordnte_variables
  
  subroutine deallocate_moving_coordnte_variables_wo_StrVars
    implicit none
    
    deallocate(coordntes)
    deallocate(coordntes_xy)
    deallocate(coordnte_of_which_node)
    deallocate(x_or_y)
    
  end subroutine deallocate_moving_coordnte_variables_wo_StrVars
  
end module moving_coordnte_variables



module conversion_routines
  use cell_info
  use system_parameters
  use moving_coordnte_info
  use moving_coordnte_variables
  
  implicit none
  
contains
  
  subroutine nodes_to_coordntes(nodes_conv,coordntes_conv)
    implicit none
    
    real*8 :: nodes_conv(1:N_node,1:N_dmnsn)
    real*8 :: coordntes_conv(1:N_mvCoordnte)
    
    integer :: i,count_coordnte
    
    count_coordnte = 1
    
    do i = 1,N_node
       
       if (node_typ(i) == 0) then !fixed node
          continue
       elseif (node_typ(i) .ne. 0) then
          
          if (node_cnnctd(i) == 0) then
             
             if (node_typ(i) == 1) then !free node
                coordntes_conv(count_coordnte:(count_coordnte+1)) = nodes_conv(i,1:N_dmnsn)
                count_coordnte = count_coordnte + 2
                
             elseif (node_typ(i) == 2) then !x free, y fixed
                coordntes_conv(count_coordnte) = nodes_conv(i,1) !nodes_conv(i,2) if x fixed, y free
                count_coordnte = count_coordnte + 1
             endif
             
          elseif (node_cnnctd(i) .ne. 0) then
             
             if (count_this_dn(i) == 1) then
                
                if (node_typ(i) == 1) then !free node
                   coordntes_conv(count_coordnte:(count_coordnte+1)) = nodes_conv(i,1:N_dmnsn)
                   count_coordnte = count_coordnte + 2
                   
                elseif (node_typ(i) == 2) then !x free, y fixed
                   coordntes_conv(count_coordnte) = nodes_conv(i,1)
                   count_coordnte = count_coordnte + 1
                endif
                
             elseif (count_this_dn(i) .ne. 1) then
                continue
             endif
             
          endif
       endif
       
    enddo
    
    
  end subroutine nodes_to_coordntes

  
  subroutine nodesl0A0_to_coordntes(nodes_conv,l0_conv,A0_conv,coordntes_conv)
    implicit none   
    real*8  :: nodes_conv(1:N_node,1:N_dmnsn)
    real*8  :: l0_conv(1:N_spr),A0_conv(1:N_cell)
    real*8  :: coordntes_conv(1:N_mvCoordnte_withl0A0)

    real*8  :: xy(1:N_mvCoordnte)
    integer :: cnt
    integer :: i,j,jmax
    
    call nodes_to_coordntes(nodes_conv,xy)

    coordntes_conv(1:N_mvCoordnte) = xy(1:N_mvCoordnte)
    cnt = N_mvCoordnte
    
    !coordntes_conv((cnt+1):(cnt+N_spr)) = l0(1:N_spr) 
    !cnt = cnt + N_spr

    !coordntes_conv((cnt+1):(cnt+N_cell)) = A0(1:N_cell)
    !cnt = cnt + N_cell

    do i=1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       if (i==1) then
          if(l0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       elseif (i==2) then
          if(A0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       endif
       
       do j = 1,jmax

          if (i==1) then

             if (optmSpr(j)==1) then
                coordntes_conv(cnt+1) = l0(j)
                cnt = cnt + 1
             elseif (optmSpr(j)==0) then
                continue
             elseif (optmSpr(j).ne.0 .AND. optmSpr(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif
             
          elseif (i==2) then

             if (optmCell(j)==1) then
                coordntes_conv(cnt+1) = A0(j)
                cnt = cnt + 1
             elseif (optmCell(j)==0) then
                continue
             elseif (optmCell(j).ne.0 .AND. optmCell(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif

          endif
         
       enddo
    enddo

    if (cnt.ne.N_mvCoordnte_withl0A0) then
       write(*,*) "cnt not equal to N_mvCoordnte_withl0A0"
       write(*,*) "flnm:info_mod,sbrtn:nodesl0A0_to_coordntes"
       stop
    endif
    
  end subroutine nodesl0A0_to_coordntes

  
  subroutine coordntes_to_nodes(coordntes_conv,nodes_conv)
    use system_parameters
    
    implicit none
    integer :: i
    real*8  :: coordntes_conv(1:N_mvCoordnte)
    real*8  :: nodes_conv(1:N_node,1:N_dmnsn)
    integer :: node_no,xy_indicator
    
    nodes_conv = node_xy
    
    !do i=1,10
    !   write(*,*) coordntes_conv(i),"bfr dte"
    !   write(*,*) node_xy(i,1:2),nodes_conv(i,1:2),"bfr the"
    !enddo
    !write(*,*) node_xy(4,1:2),"node_xy 4xy"
    !write(*,*) nodes_conv(9,1:2),"cnv9"
    !write(*,*) nodes_conv(10,1:2),"cnv10"
    
    !open(unit=102,file="node_xy_in_coor_to_nodes.dat")
    !write(unit=102,fmt=*) node_xy
    
    !write(*,*) N_mvCoordnte,"inside_coordntes_to_nodes_sbrtn"
    
    do i = 1,N_mvCoordnte
       node_no       = coordnte_of_which_node(i)
       xy_indicator  = x_or_y(i)
       
       nodes_conv(node_no,xy_indicator) = coordntes_conv(i)
    enddo
    
    !write(*,*) N_doubleNode,"Ndn"
    !write(*,*) double_node,"dn"
    
    do i=1,N_doubleNode
       nodes_conv(double_node(i,2),1:2) = nodes_conv(double_node(i,1),1:2)
    enddo
    
    !write(*,*) " "
    !do i=1,10
    !   write(*,*) coordntes_conv(i),"aft dte"
    !   write(*,*) node_xy(i,1:2),nodes_conv(i,1:2),"aft the"
    !enddo
    
  end subroutine coordntes_to_nodes
  

  subroutine coordntes_to_nodesl0A0(coordntes_conv,nodes_conv,l0_conv,A0_conv)
    implicit none
    real*8  :: coordntes_conv(1:N_mvCoordnte_withl0A0)
    real*8  :: nodes_conv(1:N_node,1:N_dmnsn)
    real*8  :: l0_conv(1:N_spr),A0_conv(1:N_cell)
    
    real*8  :: xy(1:N_mvCoordnte)
    integer :: i,j,jmax
    integer :: lam_l0variatn,lam_A0variatn
    integer :: cnt
    
    
    xy(1:N_mvCoordnte) = coordntes_conv(1:N_mvCoordnte) 
    call coordntes_to_nodes(xy,nodes_conv)
    cnt = N_mvCoordnte
    
    lam_l0variatn = 0 ; lam_A0variatn = 0
    
    !write(*,*) l0_variatn_lgcl,A0_variatn_lgcl,"l0 varitn lgcl,A0 variatn lgcl"

    do i=1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       if (i==1) then
          if(l0_variatn_lgcl.eqv..False.) then
             l0_conv = l0
             cycle
          endif
       elseif (i==2) then
          if(A0_variatn_lgcl.eqv..False.) then
             A0_conv = A0
             cycle
          endif
       endif

       do j = 1,jmax

          if (i==1) then

             if (optmSpr(j)==1) then                
                l0_conv(j) = coordntes_conv(cnt+1)
                cnt = cnt + 1
             elseif (optmSpr(j)==0) then
                l0_conv(j) = l0(j)
             elseif (optmSpr(j).ne.0 .AND. optmSpr(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif
             
          elseif (i==2) then

             if (optmCell(j)==1) then
                A0_conv(j) = coordntes_conv(cnt+1)
                cnt = cnt + 1
             elseif (optmCell(j)==0) then
                A0_conv(j) = A0(j)
             elseif (optmCell(j).ne.0 .AND. optmCell(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif

          endif
         
       enddo
    enddo

    if (cnt.ne.N_mvCoordnte_withl0A0) then
       write(*,*) "cnt not equal to N_mvCoordnte_withl0A0"
       write(*,*) "flnm:info_mod,sbrtn:nodesl0A0_to_coordntes"
       stop
    endif

    
  end subroutine coordntes_to_nodesl0A0

  subroutine coordntesxyl0A0_to_coordntes(coordntesxy_conv,l0_conv,A0_conv,coordntes_conv)
    implicit none
    real*8 :: coordntesxy_conv(1:N_mvCoordnte)
    real*8 :: l0_conv(1:N_spr)
    real*8 :: A0_conv(1:N_cell)
    real*8 :: coordntes_conv(1:N_mvCoordnte_withl0A0)
    
    integer :: i,j,jmax
    integer :: cnt
    
    coordntes_conv(1:N_mvCoordnte) = coordntesxy_conv(1:N_mvCoordnte)
    cnt = N_mvCoordnte
    
    do i=1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       if (i==1) then
          if(l0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       elseif (i==2) then
          if(A0_variatn_lgcl.eqv..False.) then
             cycle
          endif
       endif
       
       do j = 1,jmax

          if (i==1) then

             if (optmSpr(j)==1) then                
                coordntes_conv(cnt+1) = l0_conv(j)
                cnt = cnt + 1
             elseif (optmSpr(j)==0) then
                continue
             elseif (optmSpr(j).ne.0 .AND. optmSpr(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif
             
          elseif (i==2) then

             if (optmCell(j)==1) then
                coordntes_conv(cnt+1) = A0_conv(j)
                cnt = cnt + 1
             elseif (optmCell(j)==0) then
                continue
             elseif (optmCell(j).ne.0 .AND. optmCell(j).ne.1) then
                write(*,*) "Neither 0 nor 1,info_mod,nodesl0A0_to_coor"
                stop
             endif

          endif
         
       enddo
    enddo

  end subroutine coordntesxyl0A0_to_coordntes
  
  subroutine trim_1D_array(A,max,cut)
    implicit none
    integer :: max,cut
    real*8, allocatable, intent(inout) :: A(:)
    real*8  :: B(1:cut)
    
    
    if (max .lt. cut) then
       write(*,*) "inside trim_1D_array:max is less than cut"
       stop
       
    else
       B(1:cut) =  A(1:cut)
       
       deallocate(A)
       allocate(A(1:cut))
       
       A(1:cut) = B(1:cut)
    endif
    
  end subroutine trim_1D_array
  
  
  subroutine write_coordntes_as_nodes(dum_coordntes)
    implicit none
    real*8  :: dum_coordntes(1:N_mvCoordnte)
    real*8  :: nodes_tmp(1:N_node,1:N_dmnsn)
    integer :: i
    
    nodes_tmp = -1e30
    call coordntes_to_nodes(dum_coordntes,nodes_tmp)
    
    do i = 1,N_node
       if (node_typ(i) .eq. 0) then
          write(*,*) " "
          
       elseif (node_typ(i) .eq. 1) then
          write(*,*) nodes_tmp(i,1:2), i
       elseif (node_typ(i) .eq. 2) then
          write(*,*) nodes_tmp(i,1:2), i
       endif

    enddo

    
  end subroutine write_coordntes_as_nodes


  subroutine change_node_type_updt_moving_vars(hm,which,what)
    implicit none
    integer :: hm
    integer :: which(1:hm)
    integer :: what(1:hm)

    integer :: i

    write(*,*) coordntes,"coordntes"
    call coordntes_to_nodesl0A0(coordntes,node_xy,l0,A0)

    do i = 1,N_node
       write(*,*) node_xy(i,1:2)
    enddo

    do i = 1,N_spr
       write(*,*) l0(i),i,"l0_spr"
    enddo

    do i = 1,N_cell
       write(*,*) A0(i),i,"A0_cell"
    enddo
    
    do i = 1,hm
       node_typ(which(i)) = what(i)
    enddo

    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    
  end subroutine change_node_type_updt_moving_vars


  subroutine save_config_and_use_it_as_input
    implicit none

    real*8 :: A0_inpt(1:N_cell)
    real*8 :: l0_inpt(1:N_spr)
    real*8 :: nodeXY_inpt(1:N_node,1:N_dmnsn)

    integer :: i,j,jmax
    integer :: N_item

    jmax = -10
    N_item = 3
    
    do i = 1,N_item
       if (i.eq.1) jmax = N_cell
       if (i.eq.2) jmax = N_spr
       if (i.eq.3) jmax = N_node
       
       do j = 1,jmax
          if (i.eq.1) then
             A0_inpt(j) = A0(j)
          elseif (i.eq.2) then
             l0_inpt(j) = l0(j)
          elseif (i.eq.3) then
             nodeXY_inpt(j,1:N_dmnsn) = node_xy(j,1:N_dmnsn)
          endif
       enddo
       
    enddo


    do i = 1,N_item
       if (i.eq.1) jmax = N_cell
       if (i.eq.2) jmax = N_spr
       if (i.eq.3) jmax = N_node
       
       do j = 1,jmax
          if (i.eq.1) then
             A0(j) = A0_inpt(j)
          elseif (i.eq.2) then
             l0(j) = l0_inpt(j)
          elseif (i.eq.3) then
             node_xy(j,1:N_dmnsn) = nodeXY_inpt(j,1:N_dmnsn)
          endif
       enddo
       
    enddo
    
  end subroutine save_config_and_use_it_as_input

  subroutine print_coordnteVars
    implicit none

    integer :: i

    open(unit=172,file='coordnteVars.dat')
    
    do i = 1,N_mvCoordnte
       write(172,fmt=*) coordntes_xy(i),coordnte_of_which_node(i),x_or_y(i)
    enddo
    
    close(172)
    
  end subroutine print_coordnteVars
  
  
end module conversion_routines


