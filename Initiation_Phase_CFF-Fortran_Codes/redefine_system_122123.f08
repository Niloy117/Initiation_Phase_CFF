
module redefining_system_module
  use system_parameters
  use transfrm_info
  use mltiple_calls_together
  use switch_and_unswitch_models
  
  implicit none
  real*8,allocatable :: nodeXY_Saved(:,:)
  
  integer :: sprDem !spring to be demolished
  integer :: sprDemForCmbning
  real*8  :: ydis_frmVitlnMem
  real*8  :: ydis_aftCmbingSpr
  integer :: VF_BotmSurfcNodes(1:2),VF_UpprSurfcNodes(1:2)
  integer :: N_apclSideSgmntd
  integer :: INperCellT1,INperCellT2
  integer :: NumNonCrtclCell,frstCrtclCell
  
  integer              :: sprVFtopLayr
  integer, allocatable :: sprApclCntrl(:),sprApclNC1L(:),sprApclNC1R(:)
  integer, allocatable :: sprApclNC2L(:),sprApclNC2R(:)
  integer, allocatable :: sprApclNCL(:,:),sprApclNCR(:,:)
  
contains
  
  subroutine allocate_rdfnmod_vars
    implicit none
    
    allocate(nodeXY_Saved(1:N_node,1:N_dmnsn))
    
  end subroutine allocate_rdfnmod_vars
  
  subroutine deallocate_rdfnmod_vars
    implicit none
    
    if (CyclNo.gt.1) then
       deallocate(nodeXY_Saved)
    endif
    
  end subroutine deallocate_rdfnmod_vars
  
  subroutine redefine_system
    implicit none
    
    call deallocate_moving_coordnte_variables
    call deallocate_transfrm_variables
    call deallocate_phiInfo
    call deallocate_curve_variables
    
    call redefining_system_parameters_stage4
    call redefining_nodeTypes_Ends_stage4
    call get_nodeXY_coordntesXYbfrRdfn
    
    call redefining_spring_types_and_ends_stage4
    call redefining_area_types_stage4
    call get_SideWiseCellNo
    
    call redefining_CgXNode_stage4
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call get_all_moving_coordnte_variables
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call get_all_the_transfrms
    
  end subroutine redefine_system
  
  subroutine redefine_system_for_changing_only_node_typ(node_typ_chng_case)
    implicit none
    integer, intent(in) :: node_typ_chng_case
    
    call deallocate_moving_coordnte_variables_wo_StrVars
    call change_node_typ(node_typ_chng_case)
    
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
  end subroutine redefine_system_for_changing_only_node_typ
  
  subroutine change_node_typ(node_typ_chng_case)
    implicit none
    integer, intent(in) :: node_typ_chng_case
    integer             :: node1,node2,node3,node4
    integer             :: i,j
    integer             :: nodeLv,nodeRv
    
    if (node_typ_chng_case == 1) then
       
       node1 = (Hlf_Ncell+1)*2 - 1
       node2 = 2*(Hlf_Ncell+1)*2 - 1
       
       write(*,*) node1,node2,"node1-2"
       write(*,*) node_typ(node1),node_typ(node2),"curr_node_typ"
       
       node_typ(node1) = 0
       node_typ(node2) = 0
    
       write(*,*) node_typ(node1),node_typ(node2),"changed_node_typ"
       
    elseif (node_typ_chng_case == 2) then
       
       if (CellsMeet==0) then
          
          node1 = 2*(Hlf_Ncell+1)*2 + 1
          node2 = 2*(Hlf_Ncell+1)*2 + 2
          
          write(*,*) node1,node2,"node1-2"
          write(*,*) node_typ(node1),node_typ(node2),"curr_node_typ"
          
          node_typ(node1) = 1
          
       elseif (CellsMeet.gt.0) then
          write(*,*) "CellsMeet must be greater than 0"
       endif
       
    elseif (node_typ_chng_case == 3) then
       
       node1 = (Hlf_Ncell+1)*2 - 1
       node2 = 2*(Hlf_Ncell+1)*2 - 1
       node3 = node1-2
       node4 = node2-2
       
       write(*,*) node1,node2,node3,node4,"node1-2-3-4"
       write(*,*) node_typ(node1),node_typ(node2),node_typ(node3),node_typ(node4),"curr_node_typ"
       
       node_typ(node1) = 1 ; node_typ(node2) = 1
       node_typ(node3) = 1 ; node_typ(node4) = 1
       
       write(*,*) node_typ(node1),node_typ(node2),node_typ(node3),node_typ(node4),"changed_node_typ"
       
    elseif (node_typ_chng_case == 4) then
       
       do i = 1,addedNCPair
          nodeLv           = (Hlf_Ncell+1)*2 - 1 - 2*(i-1)
          nodeRv           = 2*(Hlf_Ncell+1)*2 - 1 - 2*(i-1)
          node_typ(nodeLv) = 1
          node_typ(nodeRv) = 1
       enddo
       
    endif
    
  end subroutine change_node_typ
  
  subroutine redefining_system_parameters_stage4
    implicit none

    open(unit=167,file='redefine_system.dat')
    
    strctNo = 2
    write(167,fmt=*) strctNo, "strct No"
    
    !lft_ladder_params
    ncl  = ncl - 1
    nsl  = nsecl * ncl
    nvsl = ncl + 1
    
    lft_node = 2*nvsl
    
    write(167,fmt=*) ncl,nsl,nvsl,lft_node,"ncl,nsl,nvsl,lft_Node"
    
    !rght_ladder_params
    ncr  = ncr - 1
    nsr  = nsecr * ncr
    nvsr = ncr + 1
    
    rght_node     = 2*nvsr
    
    write(167,fmt=*) ncr,nsr,nvsr,rght_node,"ncl,nsl,nvsl,lft_Node"
    
    !clft_top_params & clft_rght_params
    
    call clft_top_params
    call crght_top_params

    write(167,fmt=*) nsclt,nvsclt,clft_top_node,clft_top_spr,clft_top_end_spr,&
         clft_top_cell_nmbr,"clft_top_params"
    write(167,fmt=*) nscrt,nvscrt,crght_top_node,crght_top_spr,crght_top_end_spr,&
         crght_top_cell_nmbr,"crght_top_params"
    
    
    !clft_bot_params
    ncclb   = ncclb + 1
    nsclb   = nsecclb * ncclb
    nvsclb  = ncclb
    
    clft_bot_node      = 2*ncclb
    clft_bot_end_node  = lft_node + rght_node + clft_top_node &
         + additnl_node + crght_top_node + clft_bot_node
    
    clft_bot_spr       = nsclb
    clft_bot_end_spr   = nsl + nsr + clft_top_spr + additnl_spr &
         + crght_top_spr + clft_bot_spr
    
    first_spr_of_clft_bot = nsl+ nsr+ clft_top_spr+ crght_top_spr+ 1 
    clft_bot_cell_nmbr    = ncl + ncr + ncclt + nccrt + 2*ncclb - 1 !2 for symmetry
    write(167,fmt=*) ncclb,nsclb,nvsclb,clft_bot_node,clft_bot_end_node,&
         clft_bot_spr,clft_bot_end_spr,first_spr_of_clft_bot, clft_bot_cell_nmbr,"clft bot params"

    !crght_bot_params
    nccrb = nccrb + 1
    nscrb   = nseccrb * nccrb
    nvscrb  = nccrb

    crght_bot_node     = 2*nccrb
    crght_bot_end_node = lft_node + rght_node + clft_top_node &
         + additnl_node + crght_top_node + clft_bot_node + crght_bot_node
    
    crght_bot_spr     = nscrb
    crght_bot_end_spr = nsl + nsr + clft_top_spr + additnl_spr &
         + crght_top_spr + clft_bot_spr + crght_bot_spr

    first_spr_of_crght_bot = nsl+ nsr+ clft_top_spr+ crght_top_spr+ nsecclb+ 1  
    crght_bot_cell_nmbr    = ncl + ncr + ncclt + nccrt + ncclb + nccrb

    write(167,fmt=*) nccrb,nscrb,nvscrb,crght_bot_node,crght_bot_end_node,&
         crght_bot_spr,crght_bot_end_spr,first_spr_of_crght_bot, crght_bot_cell_nmbr,"crght bot params"
    
    call global_params !checking as redefining wont change these params
    
    write(167,fmt=*) N_node,N_spr,N_cell,Hlf_Ncell,"N_node,N_spr,N_cell,Hlf_Ncell"
    write(167,fmt=*) N_lftCells,N_rghtCells,"N lft and rght Cells"

    close(167)
    
  end subroutine redefining_system_parameters_stage4
  
  
  subroutine redefining_nodeTypes_Ends_stage4
    implicit none

    count_nodes = 0
    
    do i = 1,nvsl

       if (i .eq. 1) then
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 2
          
       elseif (i.eq.2) then
          node_typ(count_nodes+1) = 2 
          node_typ(count_nodes+2) = 1
          
       else
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 1
          
       endif
       
       if (i .eq. nvsl) then
          lft_endNode(0) = 2
          lft_endNode(1) = count_nodes + 1 
          lft_endNode(2) = count_nodes + 2
       endif
       
       count_nodes = count_nodes + 2
    enddo
    
    do i = 1,nvsr
       
       if (i==1) then
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 2
          
       elseif (i==2) then
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 1
       else
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 1
       endif
       
       if (i .eq. nvsr) then
          rght_endNode(0) = 2
          rght_endNode(1) = count_nodes + 1
          rght_endNode(2) = count_nodes + 2
       endif
       
       count_nodes = count_nodes + 2
    enddo
    
    write(*,*) count_nodes,"count_nodes rdfn1"
    
    do i = 1,clft_top_node
       
       if(i .eq. 1) then
          node_typ(count_nodes+1)    = 0
       elseif (i .eq. 2) then  
          node_typ(count_nodes+1)    = 1
       elseif (i .eq. 3) then
          node_typ(count_nodes+1)    = 1  
       endif
       
       if (i.eq.clft_top_node) then
          clft_top_endNode(0) = 3
          clft_top_endNode(1) = count_nodes - 1 
          clft_top_endNode(2) = count_nodes
          clft_top_endNode(3) = count_nodes + 1
          !write(*,*) "Entered"
       endif
       
       count_nodes = count_nodes + 1
    enddo
    
    write(*,*) clft_top_endNode(0:3),"clft"
    !stop
    
    do i = 1,crght_top_node
       
       if (i .eq. 1) then   !Pulley
          node_typ(count_nodes+1) = 0
       elseif (i .eq. 2) then
          node_typ(count_nodes+1) = 1
       else
          node_typ(count_nodes+1) = 1
       endif
       
       if (i.eq.crght_top_node) then
          crght_top_endNode(0) = 3
          crght_top_endNode(1) = count_nodes - 1 
          crght_top_endNode(2) = count_nodes
          crght_top_endNode(3) = count_nodes + 1
       endif
       
       count_nodes = count_nodes + 1
    enddo

    node_typ((count_nodes+1):N_node) = 1
    
    !clft_bot_endNode and crght_bot_endNode wont change
    
  end subroutine redefining_nodeTypes_Ends_stage4
  
  
  subroutine redefining_spring_types_and_ends_stage4
    implicit none
    integer :: i
    integer :: count_spr
    
    count_spr = 0
    
    do i = 1,nsl
       count_spr = count_spr + 1
       
       if (i .eq. 1) typ_spr(count_spr) = 1
       if (i .eq. 2) typ_spr(count_spr) = 2
       
       if (i.gt.2) then
          if (mod(i,3) .eq. 0) then
             typ_spr(count_spr) = 5
          elseif (mod(i,3) .eq. 1) then
             typ_spr(count_spr) = 3 
          else
             typ_spr(count_spr) = 4
          endif
       endif
       
       
       if (i.eq.nsl) then
          lft_endSpring(0) = 1
          lft_endSpring(1) = count_spr
       endif
    enddo
    
    do i = 1,nsr
       count_spr = count_spr + 1
       
       if (i.eq.1) typ_spr(count_spr) = 1
       if (i.eq.2) typ_spr(count_spr) = 2
       
       if (i.gt.2) then
          if (mod(i,3) .eq. 0) then
             typ_spr(count_spr) = 5
          elseif (mod(i,3) .eq. 1) then
             typ_spr(count_spr) = 3
          else
             typ_spr(count_spr) = 4
          endif
       endif

       if (i.eq.nsr) then
          rght_endSpring(0) = 1
          rght_endSpring(1) = count_spr
       endif
          
    enddo

    
    do i = 1, clft_top_spr
       count_spr = count_spr+1

       if (i.eq.1) then
          typ_spr(count_spr) = 6
          
       elseif (i .eq. 2) then 
          typ_spr(count_spr) = 7 
          
       elseif (i .eq. 3) then
          typ_spr(count_spr) = 8
          
       endif
       
       if (i.eq.1) then
          clft_top_endSpring(1) = count_spr
       elseif (i.eq.clft_top_spr) then
          clft_top_endSpring(2) = count_spr
       endif
       
    enddo
    
    
    do i = 1, crght_top_spr
       count_spr = count_spr + 1
       
       if (i .eq. 1) then
          typ_spr(count_spr) = 6
          
       elseif (i.eq.2) then
          typ_spr(count_spr) = 7
          
       elseif (i.eq.3) then
          typ_spr(count_spr) = 8
          
       endif

       if (i .eq. 1) then
          crght_top_endSpring(1) = count_spr
       elseif (i .eq. crght_top_spr) then
          crght_top_endSpring(2) = count_spr
       endif
       
    enddo
    
    
    do i = 1, nsclb
      
       if (i.ne.1 .and. mod(i,nsecclb).eq.1) then !nsecclb=3
          count_spr = count_spr + 4
       else
          count_spr = count_spr + 1
       endif
       
       if (mod(i,3).eq.1) then
          typ_spr(count_spr) = 3
       elseif (mod(i,3) .eq. 2) then 
          typ_spr(count_spr) = 4
       elseif (mod(i,3) .eq. 0) then
          typ_spr(count_spr) = 5
       endif
       
       
       if (i==nsclb) then
          clft_bot_endSpring(1) = count_spr
       endif
       
    enddo
    
    count_spr = nsl + nsr + nsclt + nscrt + nsecclb
    
    do i = 1,nscrb 
       
       if (i.ne.1 .and. mod(i,nseccrb).eq.1) then !nseccrb=3
          count_spr = count_spr + 4
       else
          count_spr = count_spr + 1
       endif
       
       if (mod(i,3).eq.1) then
          typ_spr(count_spr) = 3   
       elseif (mod(i,3) .eq. 2) then 
          typ_spr(count_spr) = 4
       elseif (mod(i,3) .eq. 0) then
          typ_spr(count_spr) = 5
       endif
       
       if (i==nscrb) then
          crght_bot_endSpring(1) = count_spr
       endif
       
    enddo
    
  end subroutine redefining_spring_types_and_ends_stage4
  
  
  subroutine redefining_area_types_stage4
    implicit none
    integer :: j
    integer :: count_area
    
    count_area = 0

    do j = 1,ncl
       count_area = count_area + 1
       typ_area(count_area) = 1
    enddo
    
    do j = 1,ncr
       count_area = count_area + 1
       typ_area(count_area) = 1
    enddo

    do j = 1,ncclt
       count_area = count_area + 1
       typ_area(count_area) = 2 
    enddo

    do j = 1,nccrt
       count_area = count_area + 1
       typ_area(count_area) = 2 
    enddo
    

    count_area = ncl + ncr + ncclt + nccrt
    
    do j = 1,ncclb
       
       if (j.eq.1) then
          count_area = count_area + 1
          typ_area(count_area) = 3
          
       elseif (j.eq.2) then
          count_area = count_area + 2
          typ_area(count_area) = 4
          
       else
          count_area = count_area + 2
          typ_area(count_area) = 5
          
       endif
       
    enddo

    count_area = ncl + ncr + ncclt + nccrt + 1
    
    do j = 1,nccrb
       
       if (j.eq.1) then
          count_area = count_area+1
          typ_area(count_area) = 3
          
       elseif (j.eq.2) then
          count_area = count_area+2
          typ_area(count_area) = 4
          
       else
          count_area = count_area+2
          typ_area(count_area) = 5
          
       endif
       
    enddo
    
  end subroutine redefining_area_types_stage4

  
  subroutine redefining_CgXNode_stage4
    implicit none
    integer :: nodes_Rdfn(1:N_node)
    integer :: i
    integer :: node_nm
    real*8  :: CgXNode_prv(1:N_node)

    CgXNode_prv = CgXNode
    
    call get_rdfnNodes(nodes_Rdfn)

    open(unit=34,file='CgX_NodeRdfn.dat')

    
    do i = 1,N_node
       node_nm = nodes_Rdfn(i)
       CgXNode(i) = CgXNode_prv(node_nm)
       write(unit=34,fmt=*) CgXNode(i),i,node_nm,"CgXNode,i,node_nm"
    enddo

    close(34)
    
  end subroutine redefining_CgXNode_stage4
  
  
  subroutine get_nodeXY_coordntesXYbfrRdfn
    implicit none
    integer :: nodes_Rdfn(1:N_node)
    integer :: i,j
    integer :: node_nm
    
    call get_rdfnNodes(nodes_Rdfn)
    
    open(unit=211,file='check_nodeXY_aftRdfn.dat')
    open(unit=212,file='node_xyBfrRdfn.dat')

    do i = 1,N_node
       node_nm = nodes_Rdfn(i)
       
       !write(*,*) nodeXY_Saved(node_nm,1:2),i,node_typ(i),"node_xy,node,typ"
       
       do j = 1,N_dmnsn
          node_xy(i,j) = nodeXY_Saved(node_nm,j)
          nodeXY_strctTN(3,i,j) = nodeXY_Saved(node_nm,j)
       enddo

       
       write(211,fmt=*) node_xy(i,1:2),i,node_typ(i),nodes_Rdfn(i),"node_xy,node,node,typ,rd"
       write(212,fmt=*) nodeXY_Saved(i,1:2),i,"nodeXYsaved"
       
    enddo
    
    !call nodes_to_coordntes(node_xy,coordntes_xy)
    
    !write(*,*) coordntes_xy,"coordntes_xy"

    close(211)
    close(212)
    
  end subroutine get_nodeXY_coordntesXYbfrRdfn

  
  subroutine get_rdfnNodes1(nodes_Rdfn)
    implicit none
    integer,intent(inout) :: nodes_Rdfn(1:N_node)
    integer :: i
    
    open(unit=305,file='nodesRdfn.dat')
    
    do i = 1,N_node
       if (i.le.8) then
          nodes_Rdfn(i) = i
          
       elseif(i.gt.8 .and. i.le.24) then
          
          if (i==9)  nodes_Rdfn(i) = 11
          if (i==10) nodes_Rdfn(i) = 12
          if (i==11) nodes_Rdfn(i) = 13
          if (i==12) nodes_Rdfn(i) = 14
          if (i==13) nodes_Rdfn(i) = 15
          if (i==14) nodes_Rdfn(i) = 16
          if (i==15) nodes_Rdfn(i) = 17
          if (i==16) nodes_Rdfn(i) = 18
          
          if (i==17) nodes_Rdfn(i) = 21
          if (i==18) nodes_Rdfn(i) = 9
          if (i==19) nodes_Rdfn(i) = 10
          
          if (i==20) nodes_Rdfn(i) = 24
          if (i==21) nodes_Rdfn(i) = 19
          if (i==22) nodes_Rdfn(i) = 20
          
          if (i==23) nodes_Rdfn(i) = 22
          if (i==24) nodes_Rdfn(i) = 23
          
       elseif (i.gt.24) then
          nodes_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd node is =",nodes_Rdfn(i), "for node old =",i
       
    enddo

    close(305)
    
    
  end subroutine get_rdfnNodes1
  
  subroutine get_rdfnNodes(nodes_Rdfn)
    implicit none
    integer, intent(inout) :: nodes_Rdfn(1:N_node)
    
    integer :: i
    integer :: CurrNvsl,PrevNvsl
    integer :: CurrNvsr,PrevNvsr
    integer :: CurrNnl,PrevNnl
    integer :: CurrNnr,PrevNnr
    integer :: cnlr,PrevCnlr
    integer :: PrevNnCntrlTop
    
    open(unit=305,file='nodesRdfnRwrt.dat')

    CurrNvsl = nvsl ; PrevNvsl = nvsl+1
    CurrNvsr = nvsr ; PrevNvsr = nvsr+1

    CurrNnl = 2*CurrNvsl ; PrevNnl = 2*PrevNvsl 
    CurrNnr = 2*CurrNvsr ; PrevNnr = 2*PrevNvsr

    cnlr = CurrNnl+CurrNnr
    PrevCnlr = PrevNnl + PrevNnr
    
    PrevNnCntrlTop = 2*PrevNvsl + 2*PrevNvsr + clft_top_node + crght_top_node
    
    do i = 1,N_node
       
       if (i.le.CurrNnl) then
          nodes_Rdfn(i) = i 
       elseif (i.gt.CurrNnl .and. i.le.cnlr) then
          nodes_Rdfn(i) = i+2 !2 is dueto 1 vs is chopped off  
       elseif (i.gt.cnlr .and. i.le.PrevNnCntrlTop) then
          
          if (i==(cnlr+1)) nodes_Rdfn(i) = PrevCnlr+1
          if (i==(cnlr+2)) nodes_Rdfn(i) = PrevNnl-1
          if (i==(cnlr+3)) nodes_Rdfn(i) = PrevNnl
          
          if (i==(cnlr+4)) nodes_Rdfn(i) = PrevCnlr+clft_top_node+1
          if (i==(cnlr+5)) nodes_Rdfn(i) = PrevCnlr-1
          if (i==(cnlr+6)) nodes_Rdfn(i) = PrevCnlr
          
          if (i==(PrevNnCntrlTop-3)) nodes_Rdfn(i) = PrevCnlr+2
          if (i==(PrevNnCntrlTop-2)) nodes_Rdfn(i) = PrevCnlr+clft_top_node
          if (i==(PrevNnCntrlTop-1)) nodes_Rdfn(i) = i
          if (i==(PrevNnCntrlTop-0)) nodes_Rdfn(i) = i
          
       elseif (i.gt.PrevNnCntrlTop) then
          nodes_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd node is =",nodes_Rdfn(i), "for node old =",i
    enddo

    close(305)

  end subroutine get_rdfnNodes
  
  
  subroutine save_nodeXY
    implicit none
    integer :: i,j
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          nodeXY_Saved(i,j) = node_xy(i,j)
       enddo
    enddo
    
  end subroutine save_nodeXY
  
  
  subroutine demolish_spring_and_redefine_system(sprNo)
    implicit none
    integer, intent(in) :: sprNo
    
    sprDem = sprNo
    
    call store_SysVars_and_Arrays_wwo_dmlsh_spr ! 0 call inside
    call get_Sysvars_dmlsh_spr ! 0 call inside
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr ! 4 calls inside
    
    call get_NodeVars_wwo_dmlsh_spr
    call get_SprVars_dmlsh_spr
    call get_Cgvars_wwo_dmlsh_spr
    call get_kphi_and_nodePhiTyp_wwo_dmlsh_spr
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_DS_trnsfrmVars
    call readjst_all_trnsfrms
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) "FINALLY CAME HERE IC"
    
  end subroutine demolish_spring_and_redefine_system
  
  subroutine adding_node_inIC_ApclMem_and_redefine_system
    implicit none
    integer :: cnt
    
    write(*,*) N_mvCoordnte,"N_mvCoordnte Bfr"
    
    call store_SysVars_and_Arrays_adding_node_inIC_ApclMem
    call get_Sysvars_addingNode_inIC_ApclMem
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr ! 4 calls inside
    
    call get_NodeVars_adding_node_inIC_ApclMem
    call get_SprVars_adding_node_inIC_ApclMem
    call get_CgVars_adding_node_inIC_ApclMem
    call get_kphi_and_nodePhiTyp_adding_node_inIC_ApclMem
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_adding_node_inIC_ApclMem_trnsfrmVars
    call readjst_all_trnsfrms_inIC_ApclMem
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) "REDEFINING inIC_Apcl"
    
  end subroutine adding_node_inIC_ApclMem_and_redefine_system
  
  
  subroutine combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM2 !CM2=CellsMeet2
    implicit none
    real*8 :: E,Ea,Es,Eg
    real*8 :: Es1,Ea1,Eg1,Eb1

    !           READ THIS FIRST ABOUT THE HASHED BLACK
    !                          ||
    !                          ||
    !                         \__/
    ! ###################################################################
    ! The following deallocate and allocate routines are NEEDED ONLY if
    ! you have multiple REDEFINE_systems in the same simulation procedures
    
    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    call deallocate_and_reallocate_transfrmStr_variables
    call deallocate_all_gradient_variables_StrVars
    call get_all_gradient_variables_StrVars
    
    ! ###################################################################
    
    call store_SysVars_and_Arrays_combining_cntrl_and_neighApcl_sprs
    call get_Sysvars_combining_cntrl_and_neighApcl_sprs
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr !4CallsInside
    
    call get_NodeVars_combining_cntrl_and_neighApcl_sprs
    call get_SprVars_combining_cntrl_and_neighApcl_sprs
    call get_CgVars_combining_cntrl_and_neighApcl_sprs
    call get_kphi_and_nodePhiTyp_combining_cntrl_and_neighApcl_sprs
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_combining_cntrl_and_neighApcl_sprs
    call readjst_all_trnsfrms_frm_reading_cellsmeet_file
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) node_xy(21,1:2),"21_1";write(*,*) node_xy(138,1:2),"138_1"
    
    call Energy_CHECK_withPrinting
    call Find_Analytical_and_Numerical_Mismatch
    !write(*,*) "stopped aft mismatch chk" ; stop
    
  end subroutine combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM2
  
  
  subroutine combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM1
    implicit none
    real*8  :: E,Ea,Es,Eg
    real*8  :: Ea1=0.0d0,Es1=0.0d0,Eg1=0.0d0,Eb1=0.0d0
    integer :: i,j,jmax
    
    ! ###################################################################
    ! The following deallocate and allocate routines are NEEDED ONLY if
    ! you have multiple REDEFINE_systems in the same simulation procedures
    
    write(*,*) Es1,Ea1,Eg1,Eb1,"Es1-Ea1-Eg1-Eb1 bfr calling even"
    write(*,*) " "
    
    do i = 1,4
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(*,*) k_spr(j), l0(j),l(j),j,"sprprop  bfr calling"
          if (i==2) write(*,*) k_area(j),A0(j),A(j),j,"cellprop bfr calling"
          if (i==3) write(*,*) CgXNode(j),CgYNode(j),k_phi(j,1:4),j,"CgX/Y,kphi prop bfr calling"
          if (i==4) write(*,*) node_xy(j,1),node_xy(j,2),j,"node_xy(j,1:2)"
       enddo
       
       write(*,*) " "
    enddo
    
    Es1 = spr_E(node_xy,l0)
    Ea1 = area_E(node_xy,A0)
    Eg1 = grvtnl_E(node_xy)
    Eb1 = bend_E(node_xy)
    
    write(*,*) Es1,Ea1,Eg1,Eb1,"Es1-Ea1-Eg1-Eb1 bfr restrcting"
    
    Es=0.0d0 ; Ea=0.0d0 ; Eg=0.0d0
    
    do i = 1,N_spr
       if (i==1) write(*,*) Es,"Es begin"
       Es = Es + 0.50d0*k_spr(i)*(l(i)-l0(i))**2
       write(*,*) i,k_spr(i),l(i),l0(i),Es,"Es bfr restrcting"
    enddo
    
    do i = 1,N_cell
       Ea = Ea + 0.50d0*k_area(i)*(A(i)-A0(i))**2
       write(*,*) i,k_area(i),A(i),A0(i),Ea,"Ea bfr restrcting"
    enddo
    
    do i = 1,N_node
       Eg = Eg + CgXNode(i)*(node_xy(i,1)) + CgYNode(i)*(node_xy(i,2))
       write(*,*) i,CgXNode(i),CgYNode(i),node_xy(i,1:2),Eg,"Eg bfr restrcting"
    enddo 
    
    
    write(*,*) N_node,N_spr,N_cell,"var 0"
    write(*,*) max_node_area,"max_node_area 0"
    write(*,*) modelID,CellsMeet,"model_ID-CellsMeet 0"
    
    
    
    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    call deallocate_and_reallocate_transfrmStr_variables
    call deallocate_all_gradient_variables_StrVars
    call get_all_gradient_variables_StrVars
    
    ! ###########################################################################################
    
    call store_SysVars_and_Arrays_combining_cntrl_and_neighApcl_sprs
    call get_Sysvars_combining_cntrl_and_neighApcl_sprs
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr ! 4 calls inside
    
    call get_NodeVars_combining_cntrl_and_neighApcl_sprs
    call get_SprVars_combining_cntrl_and_neighApcl_sprs
    call get_CgVars_combining_cntrl_and_neighApcl_sprs
    call get_kphi_and_nodePhiTyp_combining_cntrl_and_neighApcl_sprs
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_combining_cntrl_and_neighApcl_sprs

    write(*,*) N_node,N_spr,N_cell,"var"
    write(*,*) max_node_area,"max_node_area"
    write(*,*) modelID,CellsMeet,"model_ID-CellsMeet"
    
    call readjst_all_trnsfrms_frm_reading_cellsmeet_file
    
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    
    !write(*,*) node_xy(21,1:2),"21_1"
    !write(*,*) node_xy(138,1:2),"138_1"
    
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E1"
    
    
    Es=0.0d0 ; Ea=0.0d0 ; Eg=0.0d0
    
    do i = 1,N_spr
       if (i==1) write(*,*) Es,"Es begin"
       Es = Es + 0.50d0*k_spr(i)*(l(i)-l0(i))**2
       write(*,*) i,k_spr(i),l(i),l0(i),Es,"Es aft E1"
    enddo
    
    do i = 1,N_cell
       Ea = Ea + 0.50d0*k_area(i)*(A(i)-A0(i))**2
       write(*,*) i,k_area(i),A(i),A0(i),Ea,"Ea aft E1"
    enddo
    
    do i = 1,N_node
       Eg = Eg + CgXNode(i)*(node_xy(i,1)) + CgYNode(i)*(node_xy(i,2))
       write(*,*) i,CgXNode(i),CgYNode(i),node_xy(i,1:2),Eg,"Eg aft E1"
    enddo 
    
    Es1 = spr_E(node_xy,l0)
    Ea1 = area_E(node_xy,A0)
    Eg1 = grvtnl_E(node_xy)
    Eb1 = bend_E(node_xy)
    
    write(*,*) Es1,Ea1,Eg1,Eb1,"Es1-Ea1-Eg1-Eb1 bfr E2"
    
    
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E2"
    
    call Find_Analytical_and_Numerical_Mismatch
    !write(*,*) "stopped"
    !stop 'stopped aft mismatch chk'
    
  end subroutine combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM1
  
  subroutine create_two_pulleyNodeSystem_intheMiddleOfApclside
    implicit none

    
    
  end subroutine create_two_pulleyNodeSystem_intheMiddleOfApclside
  
  subroutine redefine_system_wo_demolishSpr
    implicit none
    
    call deallocate_and_reallocate_StrVars
    call store_SysVars_and_Arrays_wwo_dmlsh_spr
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr
    
    call get_NodeVars_wwo_dmlsh_spr
    call get_SprVars_wo_demlsh_spr
    call get_Cgvars_wwo_dmlsh_spr
    call get_kphi_and_nodePhiTyp_wwo_dmlsh_spr
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_DS_trnsfrmVars
    call readjst_all_trnsfrms_wo_dmlsh
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) "FINALLY CAME HERE WO"
    
  end subroutine redefine_system_wo_demolishSpr
  
  
  subroutine redefine_system_wo_demolishSpr_NI_system
    implicit none
    real*8  :: E,Ea,Es,Eg
    real*8  :: Ea1=0.0d0,Es1=0.0d0,Eg1=0.0d0,Eb1=0.0d0
    integer :: i,j,jmax
    
    call deallocate_and_reallocate_StrVars ; write(*,*) "IF prblm arise, arrays arent allocated"
    call store_SysVars_and_Arrays_wwo_dmlsh_spr
    call deallocate_and_reallocate_arrays_wwo_dmlsh_spr
    
    call get_NodeVars_wwo_dmlsh_spr_NIsystem
    call get_SprVars_wo_demlsh_spr
    call get_Cgvars_wwo_dmlsh_spr
    call get_kphi_and_nodePhiTyp_wwo_dmlsh_spr
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call alloc_and_init_Arry_forContinued_NI_trnsfrm
    allctn_Continued_NI_trnsfrmArray = 1 !IMPRTNT_LINE 
    
    call store_NI_trnsfrms_forContinuedTrnsfrm
    call deallocate_and_reallocate_DS_trnsfrmVars_NI_system
    call readjst_all_trnsfrms_frm_reading_cellsmeet_file
    
    call store_all_gradient_variables_NIcontinued
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E1"
    
    
    Es=0.0d0 ; Ea=0.0d0 ; Eg=0.0d0
    
    do i = 1,N_spr
       if (i==1) write(*,*) Es,"Es begin"
       Es = Es + 0.50d0*k_spr(i)*(l(i)-l0(i))**2
       write(*,*) i,k_spr(i),l(i),l0(i),Es,"Es aft E1"
    enddo
    
    do i = 1,N_cell
       Ea = Ea + 0.50d0*k_area(i)*(A(i)-A0(i))**2
       write(*,*) i,k_area(i),A(i),A0(i),Ea,"Ea aft E1"
    enddo
    
    do i = 1,N_node
       Eg = Eg + CgXNode(i)*(node_xy(i,1)) + CgYNode(i)*(node_xy(i,2))
       write(*,*) i,CgXNode(i),CgYNode(i),node_xy(i,1:2),Eg,"Eg aft E1"
    enddo 
    
    Es1 = spr_E(node_xy,l0)
    Ea1 = area_E(node_xy,A0)
    Eg1 = grvtnl_E(node_xy)
    Eb1 = bend_E(node_xy)
    
    write(*,*) Es1,Ea1,Eg1,Eb1,"Es1-Ea1-Eg1-Eb1 bfr E2"
    
    
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E2"
    
    call Find_Analytical_and_Numerical_Mismatch
    !write(*,*) "stopped"
    !stop 'stopped aft mismatch chk'
    write(*,*) "SLEEPING FOR 1 secs"
    call sleep(1)
    
  end subroutine redefine_system_wo_demolishSpr_NI_system
  
  
  subroutine deallocate_and_reallocate_StrVars
    implicit none
    
    deallocate(node_xyStr,node_typStr,count_this_dnStr)
    deallocate(node_cnnctdStr,double_nodeStr)
    
    deallocate(typ_sprStr,k_sprStr,l0_Str,l_Str)
    deallocate(nodePhi_typStr,k_phiStr,CgXNode_Str,CgYNode_Str)
    deallocate(coordntes_xyStr)
    
    allocate(node_xyStr(1:N_node,1:N_dmnsn))
    allocate(node_typStr(1:N_node),count_this_dnStr(1:N_node))
    allocate(node_cnnctdStr(1:N_node),double_nodeStr(1:N_node,1:N_dmnsn))
    
    allocate(typ_sprStr(1:N_spr),k_sprStr(1:N_spr))
    allocate(l0_Str(1:N_spr),l_str(1:N_spr))
    
    allocate(nodePhi_typStr(1:N_node),k_phiStr(1:N_node,1:max_Phi_node))
    allocate(CgXNode_Str(1:N_node),CgYNode_Str(1:N_node))
    allocate(coordntes_xyStr(1:N_mvCoordnte))
    
    node_xyStr       = -1.0d20  ; node_typStr    = -1
    count_this_dnStr = -1       ; node_cnnctdStr =  0
    double_nodeStr   = -1

    typ_sprStr  = -1      ; k_sprStr    = -1.0d20
    l0_Str      = -1.0d25 ; l_Str       = -1.0d25
    
    nodePhi_typStr  = -10     ; k_phiStr    = -1.0d30
    CgXNode_Str     = -1.0d30 ; CgYNode_Str = -1.0d30
    coordntes_xyStr = -1.0d30
    
  end subroutine deallocate_and_reallocate_StrVars
  
  subroutine store_SysVars_and_Arrays_wwo_dmlsh_spr
    implicit none
    
    node_xyStr       = node_xy
    node_typStr      = node_typ   
    node_cnnctdStr   = node_cnnctd 
    double_nodeStr   = double_node
    count_this_dnStr = count_this_dn
    
    typ_sprStr = typ_spr
    k_sprStr   = k_spr   
    l0_Str     = l0 
    l_Str      = l
    
    nodePhi_typStr = nodePhi_typ
    k_phiStr       = k_phi
    CgXNode_Str    = CgXNode
    CgYNode_Str    = CgYNode
    
    N_nodeS = N_node
    N_sprS  = N_spr
    N_phiS  = N_phi
    
    coordntes_xyStr = coordntes_xy
    
  end subroutine store_SysVars_and_Arrays_wwo_dmlsh_spr
  
  subroutine store_SysVars_and_Arrays_adding_node_inIC_ApclMem
    implicit none ! AN = Adding Node
    
    node_xyStrAN       = node_xy
    node_typStrAN      = node_typ   
    node_cnnctdStrAN   = node_cnnctd 
    double_nodeStrAN   = double_node
    count_this_dnStrAN = count_this_dn
    
    typ_sprStrAN = typ_spr
    k_sprStrAN   = k_spr   
    l0_StrAN     = l0 
    l_StrAN      = l
    
    nodePhi_typStrAN = nodePhi_typ
    k_phiStrAN       = k_phi
    CgXNode_StrAN     = CgXNode
    
    N_nodeSAN = N_node
    N_sprSAN  = N_spr
    N_phiSAN  = N_phi
    
    coordntes_xyStrAN = coordntes_xy
    
  end subroutine store_SysVars_and_Arrays_adding_node_inIC_ApclMem
  
  
  subroutine store_SysVars_and_Arrays_combining_cntrl_and_neighApcl_sprs
    implicit none
    
    node_xyStrES       = node_xy
    node_typStrES      = node_typ   
    node_cnnctdStrES   = node_cnnctd 
    double_nodeStrES   = double_node
    count_this_dnStrES = count_this_dn
    
    typ_sprStrES = typ_spr
    k_sprStrES   = k_spr
    l0_StrES     = l0 
    l_StrES      = l
    
    nodePhi_typStrES = nodePhi_typ
    k_phiStrES       = k_phi
    CgXNode_StrES    = CgXNode
    CgYNode_StrES    = CgYNode
    
    N_nodeSES = N_node
    N_sprSES  = N_spr
    N_phiSES  = N_phi
    
    coordntes_xyStrES = coordntes_xy
    
    write(*,*) modelID,"modelID in store"
    
  end subroutine store_SysVars_and_Arrays_combining_cntrl_and_neighApcl_sprs
  
  subroutine get_Sysvars_dmlsh_spr
    implicit none
    
    if (stageNo==1 .and. stageType==1) then
       N_node = N_node+1 !only from S1T1 to S2T1
       N_spr  = N_spr-1  !only from S1T1 to S2T1
    else
       write(*,*) "fl:redefine_sys,sb:get_Sysvars_dmlsh_spr"
    endif
    
  end subroutine get_Sysvars_dmlsh_spr
  
  subroutine get_Sysvars_addingNode_inIC_ApclMem
    implicit none
    N_node = N_node+2
  end subroutine get_Sysvars_addingNode_inIC_ApclMem
  
  subroutine get_Sysvars_combining_cntrl_and_neighApcl_sprs
    implicit none
    
    N_node = N_node-1
    N_spr  = N_spr -1  
    
  end subroutine get_Sysvars_combining_cntrl_and_neighApcl_sprs
  
  subroutine deallocate_and_reallocate_arrays_wwo_dmlsh_spr !wwo=with or without
    implicit none
    
    deallocate(node_xy,node_typ,node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(nodePhi_typ,k_phi,CgXNode,CgYNode)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_arrays_wwo_dmlsh_spr
  
  subroutine deallocate_and_reallocate_arrays_forAdding_VF_region
    implicit none
    
    deallocate(node_xy,node_typ,node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(k_area,A0,A)
    deallocate(nodePhi_typ,k_phi,CgXNode,CgYNode)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_area_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_arrays_forAdding_VF_region
  
  subroutine deallocate_and_reallocate_arrays_for_crtclDesgnd_apclSurfc
    implicit none
    
    deallocate(node_xy,node_typ,node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(k_area,A0,A)
    deallocate(nodePhi_typ,k_phi)
    deallocate(CgXNode,CgYNode)
    deallocate(AmpApclYP,AlphaApclYP,epsApclYP,powrApclYP,activtnFctrApcl)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_area_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    call allocate_and_initialize_SRypVars_wo_StrVars
    
  end subroutine deallocate_and_reallocate_arrays_for_crtclDesgnd_apclSurfc
  
  subroutine get_NodeVars_wwo_dmlsh_spr
    implicit none
    integer :: i
    integer :: lft_nodeAct,rght_nodeAct 
    integer :: lim1,lim2,lim3,lim4
    
    if (dmlshDecsn==1) then
       node_xy(1:N_nodeS,1:N_dmnsn) = node_xyStr(1:N_nodeS,1:N_dmnsn)
       write(*,*) origin(1:2) , "origin"
       node_xy(N_nodeS+1,1:N_dmnsn) = origin(1:2)
       
    elseif (dmlshDecsn==0) then
       node_xy(1:N_nodeS,1:N_dmnsn) = node_xyStr(1:N_nodeS,1:N_dmnsn)
    endif
    
    lft_nodeAct  = lft_node-2*CellsMeet  !I am not changing presuming nothing changes
    rght_nodeAct = rght_node-2*CellsMeet !no need of dmlshDecsn
    
    lim1 = lft_nodeAct
    lim2 = lft_node
    lim3 = lft_node+rght_nodeAct
    lim4 = lft_node+rght_node
    
    write(*,*) lft_nodeAct,lft_node,rght_nodeAct,rght_node,"nodes"
    write(*,*) dmlshDecsn,lim1,lim2,lim3,lim4,"lims"
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          node_typ(i) = node_typStr(i)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          node_typ(i) = 1
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          node_typ(i) = node_typStr(i)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          node_typ(i) = 1
          
       elseif (i.gt.lim4) then
          node_typ(i) = 0
          
          if (i.ne.N_node) write(*,*) i,N_node,"fl:redefine_sys,sb:rdf_nodeTyp"
          if (i.ne.N_node) stop
          
       endif   
       
    enddo
    
    lft_endNode(0)=2 ; rght_endNode(0)=2
    
    lft_endNode(1)  = lft_endNode(1)-2
    lft_endNode(2)  = lft_endNode(2)-2
    rght_endNode(1) = rght_endNode(1)-2
    rght_endNode(2) = rght_endNode(2)-2
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_wwo_dmlsh_spr
  
  subroutine get_NodeVars_VF_region
    implicit none ! IC=InitiatorCell;NC=NeighbouringCell
    
    if (N_nodeSNVFR .ne. N_node) stop 'N_node must be = N_nodeNVFR'
    
    node_xy(1:N_node,1:N_dmnsn) = node_xyStrNVFR(1:N_nodeSNVFR,1:N_dmnsn)
    node_typ(1:N_node)          = node_typStrNVFR(1:N_nodeSNVFR)
    node_cnnctd(1:N_node)       = node_cnnctdStrNVFR(1:N_nodeSNVFR)
    double_node(1:N_node,1:2)   = double_nodeStrNVFR(1:N_nodeSNVFR,1:2)
    count_this_dn(1:N_node)     = count_this_dnStrNVFR(1:N_nodeSNVFR)
    
    call print_NodeVars
    
  end subroutine get_NodeVars_VF_region
  
  subroutine get_NodeVars_crtclDesgnd_apclSurfc
    implicit none
    integer :: lim1,lim2,lim3,lim4
    integer :: INperCellNonCrtcl,INperCellCrtcl
    integer :: strtOrg,fnshOrg,strtStr,fnshStr
    integer :: i,j
    integer :: antORpost,crtclCellNum
    
    INperCellNonCrtcl = NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl     ;write(*,*)INperCellNonCrtcl,"INnonCrtc"
    INperCellCrtcl    = NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl;write(*,*)INperCellCrtcl,"INCrtcl"
    
    lim1 = numRegulrNode
    lim2 = numRegulrNode + (Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellNonCrtcl) 
    
    strtOrg=1 ; fnshOrg=lim1 ; strtStr=1 ; fnshStr=lim1
    call copy_nodexy_nodetyp_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
    strtOrg=lim1+1 ; fnshOrg=lim2 ; strtStr=lim1+1 ; fnshStr=lim2
    call copy_nodexy_nodetyp_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
    
    if (NCP_CrtclApSrfc.gt.0) then
       antORpost = -1   
       do i = 1,NCP_CrtclApSrfc
          crtclCellNum = Hlf_Ncell-(i-1)
          call find_nodeVars_for_crtclCells(antORpost,crtclCellNum)
       enddo
    endif
    
    strtOrg=numRegulrNode+(Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellNonCrtcl)+&
         (NCP_CrtclApSrfc*INperCellCrtcl)+1
    fnshOrg=strtOrg-1+(Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellNonCrtcl)
    strtStr=numRegulrNode+(Hlf_Ncell)*(INperCellNonCrtcl)+1
    fnshStr=strtStr-1+(Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellNonCrtcl)
    
    call copy_nodexy_nodetyp_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
    
    if (NCP_CrtclApSrfc.gt.0) then
       antORpost = +1   
       do i = 1,NCP_CrtclApSrfc
          crtclCellNum = Hlf_Ncell+Hlf_Ncell-(i-1)
          call find_nodeVars_for_crtclCells(antORpost,crtclCellNum)
       enddo
    endif
    
    if (VF_regionModelled==0) crtclCellNum=N_cell
    if (VF_regionModelled==1) crtclCellNum=N_cell-1
    
    call find_nodeVars_for_initiatrCells(crtclCellNum)
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_crtclDesgnd_apclSurfc
  
  subroutine copy_nodexy_nodetyp_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
    implicit none
    integer, intent(in) :: strtOrg,fnshOrg
    integer, intent(in) :: strtStr,fnshStr
    
    node_xy(strtOrg:fnshOrg,1:N_dmnsn) = node_xyStr(strtStr:fnshStr,1:N_dmnsn) 
    node_typ(strtOrg:fnshOrg)          = node_typStr(strtStr:fnshStr)
    
  end subroutine copy_nodexy_nodetyp_btwnTwoSys
  
  subroutine find_nodeVars_for_crtclCells(antORpost,crtclCellNum)
    implicit none
    integer, intent(in)  :: antORpost
    integer, intent(in)  :: crtclCellNum
    integer, allocatable :: nodesInAS(:)
    integer              :: NumNodesInAS
    integer              :: ApBsLt,n1,n2,currCrtclCellIndx
    real*8               :: fn(1:N_dmnsn),sn(1:N_dmnsn)
    integer              :: cntNodeWoCrtclCell,cntNodeCurr
    integer              :: BsLtINstrtPrv,BsLtINstrtCurr,BsLtINfnshPrv,BsLtINfnshCurr
    integer              :: i,j,imax,jmax
    
    NumNodesInAS= 2+NAEC_Apcl                            ! 2 for trmnl nodes
    allocate(nodesInAS(1:NumNodesInAS))  ; nodesInAS = -1 ; ApBsLt = 1
    call get_nodesInASide(crtclCellNum,ApBsLt,NumNodesInAS,nodesInAS)
    
    INperCellT1   = (NAEC_Apcl     +NAEC_Bsal+NAEC_Ltrl)
    INperCellT2   = (NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl)
    NumNonCrtclCell  = Hlf_Ncell-NCP_CrtclApSrfc
    frstCrtclCell = NumNonCrtclCell+1
    
    if (antORpost==-1) currCrtclCellIndx = crtclCellNum-frstCrtclCell
    if (antORpost==+1) currCrtclCellIndx = crtclCellNum-Hlf_Ncell-frstCrtclCell
    
    if (antORpost==-1) then
       cntNodeCurr        = (numRegulrNode)+(NumNonCrtclCell*INperCellT1)+(currCrtclCellIndx*INperCellT2)+1
    elseif (antORpost==+1) then
       cntNodeCurr = (numRegulrNode)+(2*NumNonCrtclCell*INperCellT1)+(NCP_CrtclApSrfc*INperCellT2) &
            +(currCrtclCellIndx*INperCellT2)+1
    endif
    
    cntNodeWoCrtclCell = (numRegulrNode)+(crtclCellNum-1)*(INperCellT1)
       
    write(*,*) antORpost,cntNodeCurr,crtclCellNum,INperCellT1,INperCellT2,frstCrtclCell,&
         currCrtclCellIndx,cntNodeWoCrtclCell,"wr 001"
    
    do i = 1,(NumNodesInAS-1)
       n1 = nodesInAS(i) ; n2 = nodesInAS(i+1)
       fn(1:2) = node_xyStr(n1,1:2) ; sn(1:2) = node_xyStr(n2,1:2)
       write(*,*) n1,n2,fn(1:2),sn(1:2),"n1,n2,fn,sn"
       
       if (i.ne.(NumNodesInAS-1)) then
          node_xy((cntNodeCurr),1:2)            = (fn(1:2)+sn(1:2))*0.5000d0
          node_xy((cntNodeCurr+1),1:2)          = sn(1:2)
          node_typ(cntNodeCurr:(cntNodeCurr+1)) = node_typStr(n2) 
          write(*,*)cntNodeCurr,node_xy((cntNodeCurr),1:2),(cntNodeCurr+1),node_xy((cntNodeCurr+1),1:2),"node1"
          write(*,*)cntNodeCurr,node_typ(cntNodeCurr),(cntNodeCurr+1),node_typ(cntNodeCurr+1),"node Tp1"
          cntNodeCurr = cntNodeCurr+2 ; write(*,*) cntNodeCurr,"CNC 1"
          
       elseif (i==(NumNodesInAS-1)) then
          node_xy((cntNodeCurr),1:2) = (fn(1:2)+sn(1:2))*0.5000d0
          node_typ(cntNodeCurr)      = node_typStr(n1)
          write(*,*) cntNodeCurr,node_xy((cntNodeCurr),1:2),"node2"
          write(*,*) cntNodeCurr,node_typ(cntNodeCurr),"node Tp2"
          cntNodeCurr = cntNodeCurr+1 ; write(*,*) cntNodeCurr,"CNC 2"
       endif
    enddo
    
    BsLtINstrtPrv  = cntNodeWoCrtclCell+NAEC_Apcl+1
    BsLtINfnshPrv  = cntNodeWoCrtclCell+NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl
    BsLtINstrtCurr = cntNodeCurr
    BsLtINfnshCurr = cntNodeCurr+NAEC_Bsal+NAEC_Ltrl-1
    
    write(*,*) BsLtINstrtPrv,BsLtINfnshPrv,BsLtINstrtCurr,BsLtINfnshCurr,"BsLtIN"
    call copy_nodexy_nodetyp_btwnTwoSys(BsLtINstrtCurr,BsLtINfnshCurr,BsLtINstrtPrv,BsLtINfnshPrv)
    write(*,*) node_typ(BsLtINstrtCurr:BsLtINfnshCurr),BsLtINstrtCurr,BsLtINfnshCurr,"01"
    
  end subroutine find_nodeVars_for_crtclCells
  
  
  subroutine find_nodeVars_for_initiatrCells(crtclCellNum)
    implicit none
    integer, intent(in)  :: crtclCellNum
    integer, allocatable :: nodesInAS(:)
    integer              :: NumNodesInAS
    integer              :: INperCellT1,INperCellT2
    integer              :: NumNonCrtclCell,frstCrtclCell
    integer              :: ApBsLt,n1,n2
    real*8               :: fn(1:N_dmnsn),sn(1:N_dmnsn)
    integer              :: cntNodeWoCrtclCell,cntNodeCurr
    integer              :: BsLtINstrtPrv,BsLtINstrtCurr,BsLtINfnshPrv,BsLtINfnshCurr
    integer              :: i,j,imax,jmax
    
    NumNodesInAS= 2+NAEC_Apcl                            ! 2 for trmnl nodes
    allocate(nodesInAS(1:NumNodesInAS))  ; nodesInAS = -1 ; ApBsLt = 1
    call get_nodesInASideOfIC(crtclCellNum,ApBsLt,NumNodesInAS,nodesInAS)
    
    INperCellT1       = (NAEC_Apcl     +NAEC_Bsal+NAEC_Ltrl)
    INperCellT2       = (NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl)
    NumNonCrtclCell      = Hlf_Ncell-NCP_CrtclApSrfc
    
    cntNodeCurr       = (numRegulrNode)+2*(NumNonCrtclCell*INperCellT1)+&
         2*(NCP_CrtclApSrfc*INperCellT2)+1
    cntNodeWoCrtclCell = (numRegulrNode)+(crtclCellNum-1)*(INperCellT1)
    
    write(*,*) cntNodeCurr,crtclCellNum,INperCellT1,INperCellT2,cntNodeWoCrtclCell, "wr002"
    
    do i = 1,(NumNodesInAS-1)
       n1 = nodesInAS(i) ; n2 = nodesInAS(i+1)
       fn(1:2) = node_xyStr(n1,1:2) ; sn(1:2) = node_xyStr(n2,1:2)
       write(*,*) n1,n2,fn(1:2),sn(1:2),"n1,n2,fn,sn"
       
       if (i.ne.(NumNodesInAS-1)) then
          node_xy((cntNodeCurr),1:2)   = (fn(1:2)+sn(1:2))*0.5000d0
          node_xy((cntNodeCurr+1),1:2) = sn(1:2)
          node_typ(cntNodeCurr:(cntNodeCurr+1)) = node_typStr(n2) 
          write(*,*)cntNodeCurr,node_xy((cntNodeCurr),1:2),cntNodeCurr+1,node_xy((cntNodeCurr+1),1:2),"node1IC"
          write(*,*)cntNodeCurr,node_typ(cntNodeCurr),(cntNodeCurr+1),node_typ(cntNodeCurr+1),"node Tp1 IC"
          cntNodeCurr = cntNodeCurr+2 ; write(*,*) cntNodeCurr,"CNC 1 IC"
       elseif (i==(NumNodesInAS-1)) then
          node_xy((cntNodeCurr),1:2) = (fn(1:2)+sn(1:2))*0.5000d0
          node_typ(cntNodeCurr)      = node_typStr(n1)
          write(*,*) cntNodeCurr,node_xy((cntNodeCurr),1:2),"node2IC"
          write(*,*) cntNodeCurr,node_typ(cntNodeCurr),"node Tp2 IC"
          cntNodeCurr = cntNodeCurr+1 ; write(*,*) cntNodeCurr,"CNC 2 IC"
       endif
    enddo
    
    BsLtINstrtPrv  = cntNodeWoCrtclCell+NAEC_Apcl+1
    BsLtINfnshPrv  = cntNodeWoCrtclCell+NAEC_Apcl+NAEC_Bsal
    BsLtINstrtCurr = cntNodeCurr
    BsLtINfnshCurr = cntNodeCurr+NAEC_Bsal-1
    
    write(*,*) BsLtINstrtPrv,BsLtINfnshPrv,BsLtINstrtCurr,BsLtINfnshCurr,"BsLtIN" 
    call copy_nodexy_nodetyp_btwnTwoSys(BsLtINstrtCurr,BsLtINfnshCurr,BsLtINstrtPrv,BsLtINfnshPrv)
    write(*,*) node_typ(BsLtINstrtCurr:BsLtINfnshCurr),BsLtINstrtCurr,BsLtINfnshCurr,"02"
    
  end subroutine find_nodeVars_for_initiatrCells
  
  
  subroutine get_NodeVars_woAddingNodesInICNC_apclSide_VF_region
    implicit none
    continue
  end subroutine get_NodeVars_woAddingNodesInICNC_apclSide_VF_region
  
  subroutine get_NodeVars_wwo_dmlsh_spr_NIsystem
    implicit none
    integer :: i
    integer :: lft_node_Cnnctd_To_VM,rght_node_Cnnctd_To_VM 
    integer :: TN_end
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    
    if (dmlshDecsn==1) then
       if (CellsMeet .ne. 0) stop 'dmlshDecsn 1 only should be t CM = 0'
       
       TN_end = 2*(Hlf_Ncell+1)*2 ! 24+24=48
       
       node_xy(1:TN_end,1:N_dmnsn)          = node_xyStr(1:TN_end,1:N_dmnsn)   !1-48
       node_xy(TN_end+1,1:N_dmnsn)          = origin(1:2)                      !49
       node_xy((TN_end+2):N_node,1:N_dmnsn) = node_xyStr((TN_end+1):N_nodeS,1:N_dmnsn)!50~Rest
       
    elseif (dmlshDecsn==0) then
       node_xy(1:N_nodeS,1:N_dmnsn) = node_xyStr(1:N_nodeS,1:N_dmnsn)
    endif
    
    lft_node_Cnnctd_To_VM  = (Hlf_Ncell+1)*2   - 2*CellsMeet  
    rght_node_Cnnctd_To_VM = 2*(Hlf_Ncell+1)*2 - 2*CellsMeet 
    
    lim1 = lft_node_Cnnctd_To_VM
    lim2 = (Hlf_Ncell+1)*2
    lim3 = rght_node_Cnnctd_To_VM
    lim4 = 2*(Hlf_Ncell+1)*2
    
    if (CellsMeet==0)   lim5 = lim4 
    if (CellsMeet.gt.0) lim5 = lim4+1
    
    lim6 = N_node
    
    write(*,*) lft_node_Cnnctd_To_VM,rght_node_Cnnctd_To_VM,dmlshDecsn,"params"
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims"
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          node_typ(i) = node_typStr(i)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          node_typ(i) = 1
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          node_typ(i) = node_typStr(i)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          node_typ(i) = 1
          
       elseif (i.gt.lim4 .and. i.le.lim5) then
          node_typ(i) = 0
          
       elseif (i.gt.lim5 .and. i.le.lim6) then
          if (CellsMeet == 0) node_typ(i) = node_typStr(i-1)
          if (CellsMeet.gt.0) node_typ(i) = node_typStr(i)
       endif
       
    enddo
    
    lft_endNode(0)  = 2   ; rght_endNode(0) = 2
    lft_endNode(1)  = ((Hlf_Ncell+1)*2-1)   - 2*CellsMeet
    lft_endNode(2)  = ((Hlf_Ncell+1)*2-0)   - 2*CellsMeet
    rght_endNode(1) = (2*(Hlf_Ncell+1)*2-1) - 2*CellsMeet
    rght_endNode(2) = (2*(Hlf_Ncell+1)*2-0) - 2*CellsMeet
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_wwo_dmlsh_spr_NIsystem
  
  
  subroutine get_NodeVars_adding_node_inIC_ApclMem
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: i,j
    integer :: lftsideTopN,rghtsideTopN
    
    lim1 = (Hlf_Ncell)*(2)    ! =22
    lim2 = lim1+2             ! =24
    lim3 = lim2+(Hlf_Ncell*2) ! =46
    lim4 = lim3+2             ! =48
    lim5 = lim4+2             ! =50
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,"lims inside get_NodeVars_adding_node_inIC_ApclMem"
    
    ydis_frmVitlnMem = 0.20d0
    write(*,*) ydis_frmVitlnMem,"ydis"
    
    open(unit=382,file='addingNodeIC_nodeVarschk.dat')
    
    do i = 1,N_node
       
       if (i.le.lim1) then ! (1-22)
          
          node_xy(i,1:N_dmnsn) = node_xyStrAN(i,1:N_dmnsn)
          node_typ(i)          = node_typStrAN(i)
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then ! (23-24)
          
          if ((i-lim1)==1) then ! 23
             node_xy(i,1) = node_xyStrAN(i,1)
             node_xy(i,2) = node_xyStrAN(i,2) - ydis_frmVitlnMem
             node_typ(i)  = 1 ! free
             
          elseif ((i-lim1)==2) then ! 24
             node_xy(i,1:N_dmnsn) = node_xyStrAN(i,1:N_dmnsn)
             node_typ(i)          = node_typStrAN(i)
          endif
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then ! (25-46)
          node_xy(i,1:N_dmnsn) = node_xyStrAN(i,1:N_dmnsn)
          node_typ(i)          = node_typStrAN(i)
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then ! (47-48)
          
          if ((i-lim3)==1) then ! 47
             node_xy(i,1) = node_xyStrAN(i,1)
             node_xy(i,2) = node_xyStrAN(i,2) - ydis_frmVitlnMem
             node_typ(i)  = 1 ! free
             
          elseif ((i-lim3)==2) then ! 48
             node_xy(i,1:N_dmnsn) = node_xyStrAN(i,1:N_dmnsn)
             node_typ(i)          = node_typStrAN(i)
          endif
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then ! (49-50)
          
          if ((i-lim4)==1) then
             
             lftsideTopN          = (Hlf_Ncell)*2 + 1
             node_xy(i,1:N_dmnsn) = node_xyStrAN(lftsideTopN,1:N_dmnsn)
             node_typ(i)          = node_typStrAN(lftsideTopN)
             
          elseif ((i-lim4)==2) then
             
             rghtsideTopN         = (Hlf_Ncell+1)*2 + (Hlf_Ncell)*2 + 1 
             node_xy(i,1:N_dmnsn) = node_xyStrAN(rghtsideTopN,1:N_dmnsn)
             node_typ(i)          = node_typStrAN(rghtsideTopN)
          endif
          
       endif
       
       write(382,*) i,node_xy(i,1:N_dmnsn),node_typ(i)
       
    enddo
    
    close(382)
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_adding_node_inIC_ApclMem
  
  
  subroutine get_NodeVars_combining_cntrl_and_neighApcl_sprs
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7,lim8
    integer :: i,j
    integer :: tbdn1,tbdn2 ! to be double node 1 and 2 
    
    ydis_aftCmbingSpr = 0.05d0
    
    if (CellsMeet==1) then
       
       lim1 = (Hlf_Ncell-1)*2 + 2 ! =22
       lim2 = lim1 + 2            ! =24
       
       lim3 = lim2 + lim1         ! =46
       lim4 = lim3 + 2            ! =48
       
       lim5 = lim4 + 1
       lim6 = N_node
       
       tbdn1 = lim1+1 ; tbdn2 = lim3+1
       write(*,*) tbdn1,tbdn2,"tbdn1 and tbdn2"
       write(*,*) origin(1:2),"origin_Val"
       
       do i = 1,N_node
          
          if (i.le.lim1) then ! (1~22)
             
             node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
             node_typ(i)          = node_typStrES(i)
             
          elseif ((i.gt.lim1) .and. (i.le.lim2)) then ! (23~24)
             
             if ((i-lim1)==1) then ! 23
                node_xy(i,1) = origin(1)
                node_xy(i,2) = (node_xyStrES(tbdn1,2) + node_xyStrES(tbdn2,2))*(0.50d0)
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim1)==2) then ! 24
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
             
          elseif ((i.gt.lim2) .and. (i.le.lim3)) then ! (25~46)
             
             node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
             node_typ(i)          = node_typStrES(i)
             
          elseif ((i.gt.lim3) .and. (i.le.lim4)) then 
             
             if ((i-lim3)==1) then ! 47
                node_xy(i,1) = origin(1)
                node_xy(i,2) = node_xy(tbdn1,2) 
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim3)==2) then ! 48
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
             
          elseif ((i.gt.lim4) .and. (i.le.lim5)) then
             
             node_xy(i,1:N_dmnsn) = origin(1:N_dmnsn)
             node_typ(i)          = 0 ! fixed node
             
          elseif ((i.gt.lim5) .and. (i.le.lim6)) then
              node_xy(i,1:N_dmnsn) = node_xyStrES(i+1,1:N_dmnsn)
              node_typ(i)          = node_typStrES(i+1)
          endif
          
       enddo
       
    elseif (CellsMeet==2) then
       
       lim1 = (Hlf_Ncell-2)*2 + 2 ! =20
       lim2 = lim1 + 2            ! =22
       lim3 = lim2 + 2            ! =24
       
       lim4 = lim3 + lim1         ! =44
       lim5 = lim4 + 2            ! =46
       lim6 = lim5 + 2            ! =48
       
       lim7  = lim6 + 1           ! =49
       lim8  = N_node             ! =139
       
       tbdn1 = lim2+1 ; tbdn2 = lim5+1
       write(*,*) tbdn1,tbdn2,"tbdn1 and tbdn2"
       write(*,*) origin(1:2),"origin_Val"
       
       do i = 1,N_node
       
          if (i.le.lim1) then ! (1~20)
             
             node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
             node_typ(i)          = node_typStrES(i)
             
          elseif ((i.gt.lim1) .and. (i.le.lim2)) then ! (21~22)
             
             if ((i-lim1)==1) then ! 21 
                node_xy(i,1) = origin(1)
                node_xy(i,2) = origin(2) - ydis_aftCmbingSpr
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim1)==2) then ! 22
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
          
          elseif ((i.gt.lim2) .and. (i.le.lim3)) then ! (23~24)
             
             if ((i-lim2)==1) then ! 23
                node_xy(i,1) = origin(1)
                node_xy(i,2) = (node_xyStrES(tbdn1,2) + node_xyStrES(tbdn2,2))*(0.50d0)
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim2)==2) then ! 24
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
          
          elseif ((i.gt.lim3) .and. (i.le.lim4)) then ! (25~44)
             
             node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
             node_typ(i)          = node_typStrES(i)
             
          elseif ((i.gt.lim4) .and. (i.le.lim5)) then ! (45~46)
             
             if ((i-lim4)==1) then ! 45 
                node_xy(i,1) = origin(1)
                node_xy(i,2) = origin(2) - ydis_aftCmbingSpr
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim4)==2) then ! 46
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
             
             
          elseif ((i.gt.lim5) .and. (i.le.lim6)) then
             
             
             if ((i-lim5)==1) then ! 47
                node_xy(i,1) = origin(1)
                node_xy(i,2) = node_xy(tbdn1,2) 
                node_typ(i)  = 1 ! free
                
             elseif ((i-lim5)==2) then ! 48
                node_xy(i,1:N_dmnsn) = node_xyStrES(i,1:N_dmnsn)
                node_typ(i)          = node_typStrES(i)
             endif
             
          elseif ((i.gt.lim6) .and. (i.le.lim7)) then ! 49
             
             node_xy(i,1:N_dmnsn) = origin(1:N_dmnsn)
             node_typ(i)          = 0 ! fixed node
             
          elseif ((i.gt.lim7) .and. (i.le.lim8)) then ! (50~Rest)
             
             node_xy(i,1:N_dmnsn) = node_xyStrES(i+1,1:N_dmnsn)
             node_typ(i)          = node_typStrES(i+1)
          endif
          
       enddo
       
    endif
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_combining_cntrl_and_neighApcl_sprs
  
  subroutine get_SprVars_dmlsh_spr
    implicit none
    integer :: i,count_spr
    
    count_spr=0
    
    do i = 1,N_sprS
       
       if (i.lt.sprDem) then
          
          count_spr          = count_spr+1
          typ_spr(count_spr) = typ_sprStr(i)
          k_spr(count_spr)   = k_sprStr(i)
          l0(count_spr)      = l0_Str(i)
          l(count_spr)       = l_Str(i)
          
       elseif (i==SprDem) then
          continue
          
       elseif (i.gt.SprDem) then
          
          count_spr          = count_spr+1
          typ_spr(count_spr) = typ_sprStr(i)
          k_spr(count_spr)   = k_sprStr(i)
          l0(count_spr)      = l0_Str(i)
          l(count_spr)       = l_Str(i)
          
       endif
       
       if (count_spr==(lft_endSpring(1)-2) .or. count_spr==(rght_endSpring(1)-2)) then
          write(*,*) lft_endSpring(1),rght_endSpring(1)
          typ_spr(count_spr) = 6
       endif
       
       !write(*,*) typ_spr(count_spr),k_spr(count_spr),l0(count_spr),count_spr,N_sprS,"sprProp"
    enddo
    
    write(*,*) "Adjust l of invaginating sprs, adding length (dist_to_pulley+dwnwrds dist)"
    call sleep(1)
    
  end subroutine get_SprVars_dmlsh_spr
  
  subroutine get_SprVars_wo_demlsh_spr
    implicit none
    integer :: meetingApclSprlft,meetingApclSprRght
    integer :: PrVmeetingApclSprlft,PrVmeetingApclSprRght
    integer :: nsprsInACell
    
    if (CellsMeet.lt.1) stop 'CM must be gt 1 in the routine, sb: get_SprVars_wo_demlsh_spr'
    
    typ_spr(1:N_sprS) = typ_sprStr(1:N_sprS)
    k_spr(1:N_sprS)   = k_sprStr(1:N_sprS)
    l0(1:N_sprS)      = l0_Str(1:N_sprS)
    l(1:N_sprS)       = l_Str(1:N_sprS)
    
    if (modelID==1) nsprsInACell = nsecl
    if (modelID==2) nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    if (modelID == 1) then
       
       meetingApclSprlft     = lft_endSpring(1)  -2-(CellsMeet-0)*(nsprsInACell)
       meetingApclSprRght    = rght_endSpring(1) -2-(CellsMeet-0)*(nsprsInACell)
       PrVmeetingApclSprlft  = lft_endSpring(1)  -2-(CellsMeet-1)*(nsprsInACell)
       PrVmeetingApclSprrght = rght_endSpring(1) -2-(CellsMeet-1)*(nsprsInACell)
       
       typ_spr(meetingApclSprlft)     = 6
       typ_spr(meetingApclSprRght)    = 6
       typ_spr(PrVmeetingApclSprlft)  = 3
       typ_spr(PrVmeetingApclSprRght) = 3
       
       write(*,*) meetingApclSprlft,meetingApclSprRght,"meetCell -m1"
       write(*,*) PrVmeetingApclSprlft,PrVmeetingApclSprRght,"PrVmeetCell -m1"
       
    elseif (modelID == 2) then
       
       meetingApclSprlft    = (Hlf_Ncell*nsprsInACell)   -2-(CellsMeet-0)*(nsprsInACell)
       meetingApclSprRght   = (2*Hlf_Ncell*nsprsInACell) -2-(CellsMeet-0)*(nsprsInACell)
       PrVmeetingApclSprlft = (Hlf_Ncell*nsprsInACell)   -2-(CellsMeet-1)*(nsprsInACell)
       PrVmeetingApclSprrght= (2*Hlf_Ncell*nsprsInACell) -2-(CellsMeet-1)*(nsprsInACell)
       
       typ_spr(meetingApclSprlft)     = 6
       typ_spr(meetingApclSprRght)    = 6
       typ_spr(PrVmeetingApclSprlft)  = 3
       typ_spr(PrVmeetingApclSprRght) = 3

       write(*,*) meetingApclSprlft,meetingApclSprRght,"meetCell - m2"
       write(*,*) PrVmeetingApclSprlft,PrVmeetingApclSprRght,"PrVmeetCell -m2"
       
    endif
    
    call print_sprVars_wo_Lt_alpha_optmSpr
    
  end subroutine get_SprVars_wo_demlsh_spr
  
  subroutine get_SprVars_adding_node_inIC_ApclMem
    implicit none
    
    ! nothing_changes in this case
    
    write(*,*) N_spr,N_sprSAN
    if (N_spr.ne.N_sprSAN) stop 'N_spr =/ N_sprSAN'
    
    typ_spr(1:N_sprSAN) = typ_sprStrAN(1:N_sprSAN)
    k_spr(1:N_sprSAN)   = k_sprStrAN(1:N_sprSAN)
    l0(1:N_sprSAN)      = l0_StrAN(1:N_sprSAN)
    l(1:N_sprSAN)       = l_StrAN(1:N_sprSAN)
    
  end subroutine get_SprVars_adding_node_inIC_ApclMem
  
  subroutine get_SprVars_for_VF_region
    implicit none
    integer :: bndryApclSpr
    integer :: sprInApclSide
    
    write(*,*) N_spr,N_sprSNVFR,"N_spr not= N_sprNVFR"
    
    typ_spr(1:N_sprSNVFR) = typ_sprStrNVFR(1:N_sprSNVFR)
    k_spr(1:N_sprSNVFR)   = k_sprStrNVFR(1:N_sprSNVFR)
    l0(1:N_sprSNVFR)      = l0_StrNVFR(1:N_sprSNVFR)
    l(1:N_sprSNVFR)       = l_StrNVFR(1:N_sprSNVFR)
    
    !typ_spr(N_sprSNVFR+1) = 10 !connects two to be meeting  Cells Dirctly
    
    !sprInApclSide = NAEC_Apcl+1 ; write(*,*)sprInApclSide,NAEC_Apcl,"sAp"
    !bndryApclSpr  = 1
    
    !l0(N_sprSNVFR+1)      = real(NAEC_Apcl+1)*(l0(bndryApclSpr)+l0(bndryApclSpr)+l0(bndryApclSpr))
    !k_spr(N_sprSNVFR+1)   = 1.00d0/l0(N_sprSNVFR+1) 
    
    !write(*,*) (N_sprSNVFR+1),k_spr(N_sprSNVFR+1),l0(N_sprSNVFR+1),"NsVF"
    
  end subroutine get_SprVars_for_VF_region
  
  subroutine get_SprVars_crtclDesgnd_apclSurfc
    implicit none
    integer :: i,j,imax,jmax
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: strtOrg,fnshOrg,strtStr,fnshStr
    integer :: cntCurrSys
    real*8  :: diffks,diffl0,diffl
    
    write(*,*) N_sprS,"N_sprS in prv"
    
    lim1 = (Hlf_Ncell-NCP_CrtclApSrfc)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    lim2 = (Hlf_Ncell)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    lim3 = (Hlf_Ncell+Hlf_Ncell-NCP_CrtclApSrfc)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    lim4 = (Hlf_Ncell+Hlf_Ncell)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    lim5 = (Hlf_Ncell+Hlf_Ncell)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3) + (NAEC_Apcl+NAEC_Bsal+2)
    write(*,*) lim1,lim2,lim3,lim4,lim5,"lim in spr apclSurfc"
    !lim1 = 90 ; lim2 = 99 ; lim3 = 189 ; lim4 = 198 ; lim5 = 204
    
    cntCurrSys = 1
    
    do i = 1,N_sprS
       
       if (i.le.lim1) then
          
          strtOrg=cntCurrSys ; fnshOrg= cntCurrSys ; strtStr=i ; fnshStr=i ! strt and end same for 1 element
          call copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
          cntCurrSys = cntCurrSys+1 
          
       elseif ((i.gt.lim1).and.(i.le.lim2)) then
          
          if ((i-lim1).le.(NAEC_Apcl+1)) then ! Apcl Spr
             call convert_spring_into_pieces_btwnTwoSys(i,cntCurrSys)
          elseif ((i-lim1).gt.(NAEC_Apcl+1)) then ! Bsal/Ltrl Spr
             strtOrg=cntCurrSys ; fnshOrg= cntCurrSys ; strtStr=i ; fnshStr=i ! strt and end same for 1 element
             call copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
             cntCurrSys = cntCurrSys+1 
          endif
          
       elseif ((i.gt.lim2).and.(i.le.lim3)) then
          strtOrg=cntCurrSys ; fnshOrg= cntCurrSys ; strtStr=i ; fnshStr=i ! strt and end same for 1 element
          call copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
          cntCurrSys = cntCurrSys+1
       elseif ((i.gt.lim3).and.(i.le.lim4)) then
          
          if ((i-lim3).le.(NAEC_Apcl+1)) then ! Apcl Spr
             call convert_spring_into_pieces_btwnTwoSys(i,cntCurrSys)
          elseif ((i-lim3).gt.(NAEC_Apcl+1)) then ! Bsal/Ltrl Spr
             strtOrg=cntCurrSys ; fnshOrg= cntCurrSys ; strtStr=i ; fnshStr=i ! strt and end same for 1 element
             call copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
             cntCurrSys = cntCurrSys+1
          endif
          
       elseif ((i.gt.lim4).and.(i.le.lim5)) then 
          if ((i-lim4).le.(NAEC_Apcl+1)) then ! Apcl Spr
             call convert_spring_into_pieces_btwnTwoSys(i,cntCurrSys)
          elseif ((i-lim4).gt.(NAEC_Apcl+1)) then ! Bsal/Ltrl Spr
             strtOrg=cntCurrSys ; fnshOrg= cntCurrSys ; strtStr=i ; fnshStr=i ! strt and end same for 1 element
             call copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
             cntCurrSys = cntCurrSys+1
          endif
       endif
       
    enddo
    
    call print_sprVars_crtclDesgnd_apclSurfc(lim1,lim2,lim3,lim4,lim5)
    
  end subroutine get_SprVars_crtclDesgnd_apclSurfc
  
  subroutine copy_sprPrps_btwnTwoSys(strtOrg,fnshOrg,strtStr,fnshStr)
    implicit none
    integer, intent(in) :: strtOrg,fnshOrg,strtStr,fnshStr
    
    write(*,*) strtOrg,fnshOrg,strtStr,fnshStr,"strtfnsh org-str"
    
    typ_spr(strtOrg:fnshOrg) = typ_sprStr(strtStr:fnshStr)
    k_spr(strtOrg:fnshOrg)   = k_sprStr(strtStr:fnshStr)
    l0(strtOrg:fnshOrg)      = l0_Str(strtStr:fnshStr) 
    l(strtOrg:fnshOrg)       = l_Str(strtStr:fnshStr)
    
  end subroutine copy_sprPrps_btwnTwoSys
  
  subroutine convert_spring_into_pieces_btwnTwoSys(CnvrtingSpr,cntCurrSys)
    implicit none
    integer, intent(in)  :: CnvrtingSpr
    integer, intent(out) :: cntCurrSys  
    integer              :: j,jmax
    real*8               :: addApNodesPerSgmntdApSprReal
    
    jmax                         = addApNodesPerSgmntdApSpr+1
    addApNodesPerSgmntdApSprReal = addApNodesPerSgmntdApSpr*(1.0000d0)  
    
    do j = 1,jmax
       typ_spr(cntCurrSys) = typ_sprStr(CnvrtingSpr)
       k_spr(cntCurrSys)   = k_sprStr(CnvrtingSpr)*(addApNodesPerSgmntdApSprReal+1.0000d0)
       l0(cntCurrSys)      = l0_Str(CnvrtingSpr)/(addApNodesPerSgmntdApSprReal+1.0000d0)
       l(cntCurrSys)       = l_Str(CnvrtingSpr)/(addApNodesPerSgmntdApSprReal+1.0000d0)
       cntCurrSys          = cntCurrSys+1
    enddo
    
  end subroutine convert_spring_into_pieces_btwnTwoSys
  
  subroutine get_SprVars_combining_cntrl_and_neighApcl_sprs
    implicit none
    integer :: i,count_spr
    integer :: lim1,lim2,lim3,lim4
    integer :: nsprsInACell
    integer :: unChngdCellL,unChngdCellR
    real*8  :: kval(1:2),lval(1:2),l0val(1:2)
    
    if (modelID==1) then
       write(*,*) "modelID should be 2 instead of",modelID
       stop
    endif
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    unChngdCellL = Hlf_Ncell-1 ; unChngdCellR = (unChngdCellL)+(Hlf_Ncell-1)
    count_spr    = 0
    
    write(*,*) nsprsInACell,unChngdCellL,unChngdCellR,"get_SprVars in combining_cntrl_and_nei"
    write(*,*) count_spr,"count_spr"
    
    sprDemForCmbning = (N_cell-1)*(nsprsInACell)+1; write(*,*) sprDemForCmbning,"sprDemForCmb"
    
    call get_spr_prp_aft_cmbining_spr(kval,lval,l0val)
    
    lim1 = (unChngdCellL*nsprsInACell)
    lim2 = lim1 + nsprsInACell
    lim3 = lim2 + lim1
    lim4 = lim3 + nsprsInACell
    
    write(*,*) N_sprSES,"SES"
    
    do i = 1,N_sprSES
       
       if (i.le.lim1) then
          
          count_spr          = count_spr+1
          typ_spr(count_spr) = typ_sprStrES(i)
          k_spr(count_spr)   = k_sprStrES(i)
          l0(count_spr)      = l0_StrES(i)
          l(count_spr)       = l_StrES(i)
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if ((i-lim1)==1) then
             
             count_spr          = count_spr+1
             typ_spr(count_spr) = typ_sprStrES(i)
             k_spr(count_spr)   = kval(1)
             l0(count_spr)      = l0val(1)
             l(count_spr)       = lval(1)
             
          elseif ((i-lim1).gt.1) then
             
             count_spr          = count_spr+1
             typ_spr(count_spr) = typ_sprStrES(i)
             k_spr(count_spr)   = k_sprStrES(i)
             l0(count_spr)      = l0_StrES(i)
             l(count_spr)       = l_StrES(i)
             
          endif
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          count_spr          = count_spr+1
          typ_spr(count_spr) = typ_sprStrES(i)
          k_spr(count_spr)   = k_sprStrES(i)
          l0(count_spr)      = l0_StrES(i)
          l(count_spr)       = l_StrES(i)
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          if ((i-lim3)==1) then
             
             count_spr          = count_spr+1
             typ_spr(count_spr) = typ_sprStrES(i)
             k_spr(count_spr)   = kval(2)
             l0(count_spr)      = l0val(2)
             l(count_spr)       = lval(2)
             
          elseif ((i-lim3).gt.1) then
             
             count_spr          = count_spr+1
             typ_spr(count_spr) = typ_sprStrES(i)
             k_spr(count_spr)   = k_sprStrES(i)
             l0(count_spr)      = l0_StrES(i)
             l(count_spr)       = l_StrES(i)
             
          endif
          
       elseif (i==SprDemForCmbning) then
          continue
          
       elseif ((i.gt.SprDemForCmbning) .and. (i.le.N_sprSES)) then
          
          count_spr          = count_spr+1
          typ_spr(count_spr) = typ_sprStrES(i)
          k_spr(count_spr)   = k_sprStrES(i)
          l0(count_spr)      = l0_StrES(i)
          l(count_spr)       = l_StrES(i)
          
       endif
       
    enddo
    
    call print_sprVars_wo_Lt_alpha_optmSpr
    
  end subroutine get_SprVars_combining_cntrl_and_neighApcl_sprs
  
  
  subroutine get_spr_prp_aft_cmbining_spr(kval,lval,l0val)
    implicit none
    real*8, intent(out) :: kval(1:2),lval(1:2),l0val(1:2)
    
    integer :: nsprsInACell,caseV
    integer :: sprL,sprC,sprR
    
    real*8  :: ks_L,ks_C,ks_R
    real*8  :: l0_L,l0_C,l0_R
    real*8  :: l_L ,l_C ,l_R
    
    real*8  :: ks_C1,ks_C2,ks_P,ks_P1,ks_P2
    real*8  :: l0_C1,l0_C2,l0_P,l0_P1,l0_P2
    real*8  :: l_C1 ,l_C2 ,l_P ,l_P1 ,l_P2
    
    real*8  :: ks_spltL(1:2),l0_spltL(1:2),l_spltL(1:2)
    real*8  :: ks_spltC(1:2),l0_spltC(1:2),l_spltC(1:2)
    real*8  :: ks_spltR(1:2),l0_spltR(1:2),l_spltR(1:2)
    
    integer             :: NsprToPr,NsprToCm
    real*8, allocatable :: ks_bcm1(:),l0_bcm1(:),l_bcm1(:)
    real*8, allocatable :: ks_bcm2(:),l0_bcm2(:),l_bcm2(:)
    real*8, allocatable :: ks_bcm3(:),l0_bcm3(:),l_bcm3(:)
    real*8, allocatable :: ks_bcm4(:),l0_bcm4(:),l_bcm4(:)
    
    real*8  :: ks_acm1,ks_acm2,ks_acm3,l0_acm1,l0_acm2,l0_acm3,l_acm1,l_acm2,l_acm3
    real*8  :: ks_acm4,l0_acm4,l_acm4
    integer :: choice
    
    open(unit=389,file='get_spr_cmbining.dat')
    
    nsprsInACell = (NAEC_Apcl+1+NAEC_Bsal+1+NAEC_Ltrl+1)
    
    sprL = (Hlf_Ncell-1)*(nsprsInACell) + 1
    sprC = (N_cell-1)*(nsprsInACell)    + 1
    sprR = (N_cell-2)*(nsprsInACell)    + 1
    
    write(*,*) sprL,sprC,sprR,"sprS bfr Cmbining"
    
    if (CellsMeet == 1) then
       
       write(389,*) k_sprStrES(sprL),k_sprStrES(sprC),k_sprStrES(sprR),"k_spr strES"
       write(389,*) l0_strES(sprL),l0_strES(sprC),l0_strES(sprR),"l0 strES"
       write(389,*) l_strES(sprL),l_strES(sprC),l_strES(sprR),"l strES"
       
       choice = 2
       call split_spr_srialy_with_ratio_tbd(sprL,choice,ks_spltL,l0_spltL,l_spltL) 
       call split_spr_srialy_equal_len(sprC,choice,ks_spltC,l0_spltC,l_spltC)
       call split_spr_srialy_with_ratio_tbd(sprR,choice,ks_spltR,l0_spltR,l_spltR)
       
       write(389,*) ks_spltL(1:2),ks_spltC(1:2),ks_spltR(1:2),"ks_splt L-C-R"
       write(389,*) l0_spltL(1:2),l0_spltC(1:2),l0_spltR(1:2),"l0_splt L-C-R"
       write(389,*) l_spltL(1:2), l_spltC(1:2), l_spltR(1:2), "l_splt L-C-R"
       
       NsprToPr = 2
       allocate(ks_bcm1(1:NsprToPr),l0_bcm1(1:NsprToPr),l_bcm1(1:NsprToPr))
       allocate(ks_bcm2(1:NsprToPr),l0_bcm2(1:NsprToPr),l_bcm2(1:NsprToPr))
       
       ks_bcm1=-1.0d30 ; l0_bcm1=-1.0d30 ; l_bcm1=-1.0d30
       ks_bcm2=-1.0d30 ; l0_bcm2=-1.0d30 ; l_bcm2=-1.0d30
       
       NsprToCm = 2
       caseV    = 1
       call get_ks_prop_bfr_comb(NsprToCm,caseV,ks_spltL,ks_spltC,ks_spltR,ks_bcm1)
       call get_l0_prop_bfr_comb(NsprToCm,caseV,l0_spltL,l0_spltC,l0_spltR,l0_bcm1)
       call get_l_prop_bfr_comb(NsprToCm, caseV,l_spltL ,l_spltC ,l_spltR ,l_bcm1 )
       
       write(389,*) " "
       write(389,*) ks_bcm1(1:2),l0_bcm1(1:2),l_bcm1(1:2),"bcm1"
       
       caseV = 2
       call get_ks_prop_bfr_comb(NsprToCm,caseV,ks_spltL,ks_spltC,ks_spltR,ks_bcm2)
       call get_l0_prop_bfr_comb(NsprToCm,caseV,l0_spltL,l0_spltC,l0_spltR,l0_bcm2)
       call get_l_prop_bfr_comb(NsprToCm, caseV,l_spltL ,l_spltC ,l_spltR ,l_bcm2 )
       
       write(389,*) ks_bcm2(1:2),l0_bcm2(1:2),l_bcm2(1:2),"bcm2"
       
       call cmbine_spr_prlelly_cnnctd(NsprToCm,ks_bcm1,l0_bcm1,l_bcm1,ks_acm1,l0_acm1,l_acm1)
       call cmbine_spr_prlelly_cnnctd(NsprToCm,ks_bcm2,l0_bcm2,l_bcm2,ks_acm2,l0_acm2,l_acm2)
       
       write(389,*) " "
       write(389,*) ks_acm1,l0_acm1,l_acm1,"acm1"
       write(389,*) ks_acm2,l0_acm2,l_acm2,"acm2"
       
       NsprToCm   = 2
       allocate(ks_bcm3(1:NsprToCm),l0_bcm3(1:NsprToCm),l_bcm3(1:NsprToCm))
       
       ks_bcm3(1) = ks_spltL(1) ; ks_bcm3(2) = ks_acm1
       l0_bcm3(1) = l0_spltL(1) ; l0_bcm3(2) = l0_acm1
       l_bcm3(1)  = l_spltL(1)  ; l_bcm3(2)  = l_acm1
       write(389,*) ks_bcm3(1:2),l0_bcm3(1:2),l_bcm3(1:2),"bcm3"
       
       call cmbine_spr_serially_cnnctd(NsprToCm,ks_bcm3,l0_bcm3,l_bcm3,ks_acm3,l0_acm3,l_acm3)
       write(389,*) ks_acm3,l0_acm3,l_acm3,"acm3"
       
       NsprToCm   = 2
       allocate(ks_bcm4(1:NsprToCm),l0_bcm4(1:NsprToCm),l_bcm4(1:NsprToCm))
       
       ks_bcm4(1) = ks_spltR(1) ; ks_bcm4(2) = ks_acm2
       l0_bcm4(1) = l0_spltR(1) ; l0_bcm4(2) = l0_acm2
       l_bcm4(1)  = l_spltR(1)  ; l_bcm4(2)  = l_acm2
       write(389,*) ks_bcm4(1:2),l0_bcm4(1:2),l_bcm4(1:2),"bcm4"
       
       call cmbine_spr_serially_cnnctd(NsprToCm,ks_bcm4,l0_bcm4,l_bcm4,ks_acm4,l0_acm4,l_acm4)
       write(389,*) ks_acm4,l0_acm4,l_acm4,"acm4"
       
       kval(1)  = ks_acm3 ; kval(2)  = ks_acm4
       l0val(1) = l0_acm3 ; l0val(2) = l0_acm4
       lval(1)  = l_acm3  ; lval(2)  = l_acm4
       
       write(389,*) kval(1:2), "kval"
       write(389,*) l0val(1:2),"l0val"
       write(389,*) lval(1:2), "lval"
       
    elseif (CellsMeet == 2) then
       
       ks_L = k_sprStrES(sprL) ; ks_C = k_sprStrES(sprC) ; ks_R = k_sprStrES(sprR)
       l0_L = l0_StrES(sprL)   ; l0_C = l0_StrES(sprC)   ; l0_R = l0_StrES(sprR) 
       l_L  = l_StrES(sprL)    ; l_C  = l_StrES(sprC)    ; l_R  = l_StrES(sprR)
       
       write(*,*) ks_L,ks_C,ks_R,"ks bfr L,C,R"
       write(*,*) l0_L,l0_C,l0_R,"l0 bfr L,C,R"
       write(*,*) l_L ,l_C ,l_R, "l  bfr L,C,R"
       
       ! Splitting Central Spring into Two
       
       ks_C1 = (ks_C)*(2.0d0) ; ks_C2 = ks_C1
       l0_C1 = (l0_C)/(2.0d0) ; l0_C2 = l0_C1
       l_C1  = (l_C) /(2.0d0) ; l_C2  = l_C1
       
       write(*,*) ks_C1,ks_C2,"ks_C's aft split"
       write(*,*) l0_C1,l0_C2,"l0_C's aft split"
       write(*,*) l_C1 ,l_C2 ,"l_C's aft split"
       
       ! combinging four springs
       
       ks_P = ks_L + ks_C1 + ks_C2 + ks_R
       l0_P = (1.0d0/4.0d0) * (l0_L + l0_C1 + l0_C2 + l0_R)
       l_P  = (1.0d0/4.0d0) * (l_L  + l_C1  + l_C2  + l_R )
       
       write(*,*) ks_P,l0_P,l_P,"ks-l0-l P's aft combining"
       
       ! splitting single spring parallelly into two
       
       ks_P1 = ks_P/2.0d0 ; ks_P2 = ks_P1
       l0_P1 = l0_P       ; l0_P2 = l0_P1
       l_P1  = l_P        ; l_P2  = l_P1
       
       write(*,*) ks_P1,ks_P2,"ks_P1 & ks_P2"
       write(*,*) l0_P1,l0_P2,"l0_P1 & l0_P2"
       write(*,*) l_P1 ,l_P2 ,"l_P1  & l_P2 "
       
       kval(1)  = ks_P1 ; kval(2)  = ks_P2
       l0val(1) = l0_P1 ; l0val(2) = l0_P2
       lval(1)  = l_P1  ; lval(2)  = l_P2
       
       write(*,*) kval(1:2),"kval"
       write(*,*) l0val(1:2),"l0val"
       write(*,*) lval(1:2),"lval"
       
    endif
    
    close(389)
    
  end subroutine get_spr_prp_aft_cmbining_spr
  
  subroutine get_ks_prop_bfr_comb(NsprToCm,caseV,ks_spltL,ks_spltC,ks_spltR,ks_bp)
    implicit none
    integer, intent(in)  :: NsprToCm,caseV
    real*8 , intent(in)  :: ks_spltL(1:2),ks_spltC(1:2),ks_spltR(1:2)
    real*8 , intent(out) :: ks_bp(1:NsprToCm)
    
    if (caseV==1) then
       ks_bp(1) = ks_spltL(2) ; ks_bp(2) = ks_spltC(1)
       write(*,*) ks_bp(1:2),"ks_bp1"
    elseif (caseV==2) then   
       ks_bp(1)=ks_spltR(2) ; ks_bp(2)=ks_spltC(2)
       write(*,*) ks_bp(1:2),"ks_bp2"
    endif
    
  end subroutine get_ks_prop_bfr_comb
  
  subroutine get_l0_prop_bfr_comb(NsprToCm,caseV,l0_spltL,l0_spltC,l0_spltR,l0_bp)
    implicit none
    integer, intent(in)  :: NsprToCm,caseV
    real*8 , intent(in)  :: l0_spltL(1:2),l0_spltC(1:2),l0_spltR(1:2)
    real*8 , intent(out) :: l0_bp(1:NsprToCm)
    
    if (caseV==1) then
       l0_bp(1) = l0_spltL(2) ; l0_bp(2) = l0_spltC(1)
       write(*,*) l0_bp(1:2),"l0_bp1"
    elseif (caseV==2) then   
       l0_bp(1)=l0_spltR(2) ; l0_bp(2)=l0_spltC(2)
       write(*,*) l0_bp(1:2),"l0_bp2"
    endif
    
  end subroutine get_l0_prop_bfr_comb
  
  subroutine get_l_prop_bfr_comb(NsprToCm,caseV,l_spltL,l_spltC,l_spltR,l_bp)
    implicit none
    integer, intent(in)  :: NsprToCm,caseV
    real*8 , intent(in)  :: l_spltL(1:2),l_spltC(1:2),l_spltR(1:2)
    real*8 , intent(out) :: l_bp(1:NsprToCm)
    
    if (caseV==1) then
       l_bp(1) = l_spltL(2) ; l_bp(2) = l_spltC(1)
       write(*,*) l_bp(1:2),"l_bp1"
    elseif (caseV==2) then   
       l_bp(1)=l_spltR(2) ; l_bp(2)=l_spltC(2)
       write(*,*) l_bp(1:2),"l_bp2"
    endif
    
  end subroutine get_l_prop_bfr_comb
  
  subroutine get_CellVars_VF_region
    implicit none
    real*8 :: VF_area
    
    write(*,*) N_cell,N_cellSNVFR,"N_cell in VF"
    
    k_area(1:N_cellSNVFR) = k_areaStrNVFR(1:N_cellSNVFR)
    A0(1:N_cellSNVFR)     = A0_StrNVFR(1:N_cellSNVFR)
    A(1:N_cellSNVFR)      = A_StrNVFR(1:N_cellSNVFR)
    
    !call vitelline_fluid_region(VF_area)
    VF_area         = 0.0000d0
    
    A(N_cell)       = VF_area
    A0(N_cell)      = A(N_cell)
    k_area(N_cell)  = k_area(1) !(1.0000d0)/A0(1)
    
    write(*,*) k_area(N_cell),A0(N_cell),A(N_cell),"VF area prop"
    
  end subroutine get_CellVars_VF_region
  
  subroutine get_CellVars_crtclDesgnd_apclSurfc
    implicit none
    write(*,*) N_cell,N_cellS,"N_cell~N_cellS"
    if (N_cell .ne. N_cellS) stop 'inconsistent N_cell'
    
    k_area(1:N_cell) = k_areaStr(1:N_cell)
    A0(1:N_cell)     = A0_Str(1:N_cell)
    A(1:N_cell)      = A_Str(1:N_cell)
    
  end subroutine get_CellVars_crtclDesgnd_apclSurfc
  
  subroutine get_Cgvars_wwo_dmlsh_spr
    implicit none
    
    CgXNode(1:N_node)  = 0.00d0 ; CgYNode(1:N_node) = 0.00d0
    CgXNode(1:N_nodeS) = CgXNode_Str(1:N_nodeS)
    CgYNode(1:N_nodeS) = CgYNode_Str(1:N_nodeS)
    
  end subroutine get_Cgvars_wwo_dmlsh_spr
  
  subroutine get_CgVars_adding_node_inIC_ApclMem
    implicit none
    
    CgXNode(1:N_node) = 0.00d0 ; CgYNode(1:N_node) = 0.00d0
    CgXNode(1:N_nodeSAN) = CgXNode_StrAN(1:N_nodeSAN)
    
    !CgYNode(1:N_nodeSAN) = CgYNode_StrAN(1:N_nodeSAN)
    
  end subroutine get_CgVars_adding_node_inIC_ApclMem
  
  subroutine get_CgVars_VF_region
    implicit none
    
    CgXNode(1:N_nodeSNVFR) = CgXNode_StrNVFR(1:N_nodeSNVFR)
    CgYNode(1:N_nodeSNVFR) = CgYNode_StrNVFR(1:N_nodeSNVFR)
    
  end subroutine get_CgVars_VF_region
  
  subroutine get_CgVars_crtclDesgnd_apclSurfc
    implicit none
    integer :: i,j,imax,jmax
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: cntCurrNode,cntPrvNode
    real*8  :: TINY=1.0d-15
    
    write(*,*) INperCellT1,INperCellT2,NumNonCrtclCell,frstCrtclCell,"info At the begin"
    
    lim1 = numRegulrNode
    lim2 = lim1 + NumNonCrtclCell*INperCellT1
    lim3 = lim1 + Hlf_Ncell*INperCellT1
    lim4 = lim3 + NumNonCrtclCell*INperCellT1
    lim5 = lim3 + Hlf_Ncell*INperCellT1
    lim6 = lim5 + (NAEC_Apcl+NAEC_Bsal+2)
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims in CgVars_crtclDesgnd"
    cntCurrNode = 1
    
    do i = 1,N_nodeS
       
       if (i.le.lim1) then
          CgXNode(cntCurrNode)  = CgXNode_Str(i)
          CgYNode(cntCurrNode)  = CgYNode_Str(i)
          cntCurrNode           = cntCurrNode+1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          CgXNode(cntCurrNode)  = CgXNode_Str(i)
          CgYNode(cntCurrNode)  = CgYNode_Str(i)
          cntCurrNode           = cntCurrNode+1
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          if ((i-lim2).le.(NAEC_Apcl)) then ! Apcl Node
             
             if ((i-lim2)==1) then
                cntPrvNode = i
                call get_CgXnYn_distrbtn_forInsrtdNodes(cntPrvNode,cntCurrNode)
                write(*,*) cntPrvNode,cntCurrNode,"cnt Prv Curr in i =",i
             elseif ((i-lim2).gt.1) then
                write(*,*) 'cases being considered in first apical node'
                cycle
             endif
             
          elseif ((i-lim2).gt.(NAEC_Apcl)) then
             CgXNode(cntCurrNode)  = CgXNode_Str(i)
             CgYNode(cntCurrNode)  = CgYNode_Str(i)
             cntCurrNode           = cntCurrNode+1
          endif
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          CgXNode(cntCurrNode)  = CgXNode_Str(i)
          CgYNode(cntCurrNode)  = CgYNode_Str(i)
          cntCurrNode           = cntCurrNode+1
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          if ((i-lim4).le.(NAEC_Apcl)) then ! Apcl Node
             
             if ((i-lim4)==1) then
                cntPrvNode = i
                call get_CgXnYn_distrbtn_forInsrtdNodes(cntPrvNode,cntCurrNode)
                write(*,*) cntPrvNode,cntCurrNode,"cnt Prv Curr in i =",i
             elseif ((i-lim4).gt.1) then
                write(*,*) 'cases being considered in first apical node'
                cycle
             endif
             
          elseif ((i-lim4).gt.(NAEC_Apcl)) then
             CgXNode(cntCurrNode)  = CgXNode_Str(i)
             CgYNode(cntCurrNode)  = CgYNode_Str(i)
             cntCurrNode           = cntCurrNode+1
          endif
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          
          if ((i-lim5).le.(NAEC_Apcl)) then ! Apcl Node
             
             if ((i-lim5)==1) then
                cntPrvNode = i
                call get_CgXnYn_distrbtn_forInsrtdNodes(cntPrvNode,cntCurrNode)
                write(*,*) cntPrvNode,cntCurrNode,"cnt Prv Curr in i =",i
             elseif ((i-lim5).gt.1) then
                write(*,*) 'cases being considered in first apical node'
                cycle
             endif
             
          elseif ((i-lim5).gt.(NAEC_Apcl)) then
             CgXNode(cntCurrNode)  = CgXNode_Str(i)
             CgYNode(cntCurrNode)  = CgYNode_Str(i)
             cntCurrNode           = cntCurrNode+1
          endif
          
       endif
       
    enddo
    
  end subroutine get_CgVars_crtclDesgnd_apclSurfc
  
  subroutine get_CgXnYn_distrbtn_forInsrtdNodes(cntPrvNode,cntCurrNode)
    implicit none
    integer, intent(inout) :: cntPrvNode,cntCurrNode
    integer                :: cntPrvNodeSave,cntCurrNodeSave
    real*8                 :: TINYval=1.0d-15
    integer                :: i,j
    real*8                 :: sumOfCgY
    
    cntPrvNodeSave=cntPrvNode; cntCurrNodeSave=cntCurrNode; write(*,*)cntPrvNodeSave,cntCurrNodeSave,"cntChk 1"
    
    if (abs(CgXNode_Str(cntPrvNode)).gt.TINYval) then
       write(*,*) CgXNode_Str(cntPrvNode),cntPrvNode,"CgXNode val"
       stop 'CgXNode must not be gt Zero in InsrtdNode'
    endif
    
    CgXNode(cntCurrNode:(cntCurrNode+NAEC_ApclCrtcl-1)) = 0.0000d0
    write(*,*) cntCurrNode,"cntCurrNode val"
    
    if (CgYNode_Str(cntPrvNode).le.TINYval) then
       do i = 1,NAEC_ApclCrtcl
          CgYNode(cntCurrNode) = 0.0000d0
          cntCurrNode          = cntCurrNode+1
       enddo
       
    elseif (CgYNode_Str(cntPrvNode).gt.TINYval) then
       sumOfCgY = 0.00d0
       do i = 1,NAEC_Apcl
          sumOfCgY  = sumOfCgY + CgYNode_Str(cntPrvNode+i-1)
       enddo
       
       write(*,*) sumOfCgY,"sumOfCgY"
       
       do i = 1,NAEC_ApclCrtcl
          CgYNode(cntCurrNode) = (sumOfCgY)/real(NAEC_ApclCrtcl*1.0000d0) 
          cntCurrNode          = cntCurrNode+1
       enddo
       
    endif
    
    write(*,*) cntPrvNode,cntPrvNodeSave,cntCurrNode,"cntChk 2"
    if (cntPrvNode.ne.cntPrvNodeSave) stop 'cnt prv should not change'
    
  end subroutine get_CgXnYn_distrbtn_forInsrtdNodes
  
  
  subroutine get_CgVars_combining_cntrl_and_neighApcl_sprs
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: i,j
    integer :: cnnctdNodeVal
    
    lim1 = (Hlf_Ncell+1)*4 + 1 ! =49
    lim2 = N_node              ! =139
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          
          if (node_cnnctd(i) == 0) then
             CgXNode(i) = CgXNode_StrES(i)
             CgYNode(i) = CgYNode_StrES(i)
             
          elseif (node_cnnctd(i) .ne. 0) then
             cnnctdNodeVal = node_cnnctd(i)
             write(*,*) cnnctdNodeVal,"cnnNodeVal"
             
             CgXNode(i)    = CgXNode_StrES(i) + CgXNode_StrES(cnnctdNodeVal) 
             CgYNode(i)    = CgYNode_StrES(i) + CgYNode_StrES(cnnctdNodeVal)
             
          endif
          
       elseif ((i.gt.lim1).and.(i.le.lim2)) then ! inserted nodes
          
          CgXNode(i) = CgXNode_StrES(i+1)
          CgYNode(i) = CgYNode_StrES(i+1)
          
       endif
       
    enddo
    
  end subroutine get_CgVars_combining_cntrl_and_neighApcl_sprs
  
  subroutine get_kphi_and_nodePhiTyp_wwo_dmlsh_spr
    implicit none
    real*8  :: k_phiVal
    
    k_phiVal = k_phiStr(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal !that's because no bending
    
    nodePhi_typ(1:N_node) = 1
    
  end subroutine get_kphi_and_nodePhiTyp_wwo_dmlsh_spr
  
  subroutine get_kphi_and_nodePhiTyp_adding_node_inIC_ApclMem
    implicit none
    real*8 :: k_phiVal
    
    k_phiVal                       = k_phiStrAN(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal !that's because no bending
    nodePhi_typ(1:N_node)          = 1
    
  end subroutine get_kphi_and_nodePhiTyp_adding_node_inIC_ApclMem
  
  subroutine get_kphi_and_nodePhiTyp_VF_region
    implicit none
    real*8 :: k_phiVal
    
    write(*,*) N_node,N_nodeSNVFR,"N_node must be equal to N_nodeSNVFR"
    
    k_phiVal                       = k_phiStrNVFR(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal
    nodePhi_typ(1:N_node)          = 1
    
  end subroutine get_kphi_and_nodePhiTyp_VF_region

  subroutine get_kphi_and_nodePhiTyp_crtclDesgnd_apclSurfc
    implicit none
    real*8 :: k_phiVal
    
    write(*,*) N_node,N_nodeS,"N_node must be equal to N_nodeSNVFR"
    
    k_phiVal                       = k_phiStr(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal
    nodePhi_typ(1:N_node)          = 1
    
  end subroutine get_kphi_and_nodePhiTyp_crtclDesgnd_apclSurfc
  
  
  subroutine get_kphi_and_nodePhiTyp_combining_cntrl_and_neighApcl_sprs
    implicit none
    real*8 :: k_phiVal
    
    k_phiVal                       = k_phiStrES(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal !that's because no bending
    nodePhi_typ(1:N_node)          = 1
    
  end subroutine get_kphi_and_nodePhiTyp_combining_cntrl_and_neighApcl_sprs
  
  
  subroutine get_SRyp_props_crtclDesgnd_apclSurfc
    implicit none
    
    call get_SRYP_props ; call get_modified_ApmApclYP_for_AddedinsrtdNodes
    call readjst_get_actvtn_fctrs_crtclDesgnd_apclSurfc
    
  end subroutine get_SRyp_props_crtclDesgnd_apclSurfc
  
  subroutine get_modified_ApmApclYP_for_AddedinsrtdNodes
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: cntCurrApNode,cntPrvApNode
    integer :: i,j
    
    write(*,*) numApclNodes,"NumApclNodes"
    
    lim1 = (Hlf_Ncell-NCP_CrtclApSrfc)*(NAEC_Apcl+1)
    lim2 = lim1 + (NAEC_Apcl+1) + 1
    lim3 = lim2 + (Hlf_Ncell-NCP_CrtclApSrfc)*(NAEC_Apcl+1)
    lim4 = lim3 + (NAEC_Apcl+1) + 1
    lim5 = lim4 + (NAEC_Apcl)
    
    cntCurrApNode = 1
    
    do i = 1,numApclNodes
       
       if (i.le.lim1) then
          cntCurrApNode = cntCurrApNode+1
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if (((i-lim1)==1) .or. (i==lim2)) then
             cntCurrApNode = cntCurrApNode+1
          else   
             if ((i-lim1)==2) then
                cntPrvApNode = i
                call get_AmpApcl_distrbtn_forInsrtdNodes(cntPrvApNode,cntCurrApNode)
                write(*,*) cntPrvApNode,cntCurrApNode,"cnt Prv Curr in i =",i 
             elseif ((i-lim1).gt.2) then
                continue
             endif
          endif
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          cntCurrApNode = cntCurrApNode+1
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          if (((i-lim3)==1) .or. (i==lim4)) then
             cntCurrApNode = cntCurrApNode+1
          else
             if ((i-lim3)==2) then
                cntPrvApNode = i
                call get_AmpApcl_distrbtn_forInsrtdNodes(cntPrvApNode,cntCurrApNode)
                write(*,*) cntPrvApNode,cntCurrApNode,"cnt Prv Curr in i =",i 
             elseif ((i-lim3).gt.2) then
                continue
             endif
          endif
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          if ((i-lim4)==1) then
             cntPrvApNode = i
             call get_AmpApcl_distrbtn_forInsrtdNodes(cntPrvApNode,cntCurrApNode)
             write(*,*) cntPrvApNode,cntCurrApNode,"cnt Prv Curr in i =",i 
          elseif ((i-lim4).gt.1) then
             continue
          endif       
       endif
       
    enddo
    
    
  end subroutine get_modified_ApmApclYP_for_AddedinsrtdNodes
  
  subroutine get_AmpApcl_distrbtn_forInsrtdNodes(cntPrvApNode,cntCurrApNode)
    implicit none
    integer, intent(inout) :: cntPrvApNode,cntCurrApNode
    integer                :: cntPrvApNodeSave,cntCurrApNodeSave
    real*8                 :: TINYval=1.0d-15
    integer                :: i,j
    real*8                 :: sumOfAmpApclYP
    
    cntPrvApNodeSave=cntPrvApNode; cntCurrApNodeSave=cntCurrApNode; write(*,*)cntPrvApNodeSave,cntCurrApNodeSave,"cntChk 1"
    
    sumOfAmpApclYP = 0.00d0
    do i = 1,NAEC_Apcl
       sumOfAmpApclYP  = sumOfAmpApclYP + AmpApclYP_Str(cntPrvApNode+i-1)
    enddo
    
    write(*,*) sumOfAmpApclYP,"sumOfAmpApclYP"
    
    do i = 1,NAEC_ApclCrtcl
       AmpApclYP(cntCurrApNode) = (sumOfAmpApclYP)/real(NAEC_ApclCrtcl) 
       cntCurrApNode            = cntCurrApNode+1
    enddo
    
    write(*,*) cntPrvApNode,cntPrvApNodeSave,cntCurrApNode,"cntChk 2"
    if (cntPrvApNode.ne.cntPrvApNodeSave) stop 'cnt prv should not change'
      
  end subroutine get_AmpApcl_distrbtn_forInsrtdNodes
  
  subroutine readjst_get_actvtn_fctrs_crtclDesgnd_apclSurfc()
    implicit none
    integer :: NumActvCell
    integer :: NNApclL,NNApclR           !NNApclL=Number of Nodes in Apical (Left/Right)
    integer :: InactvNodesL,InactvNodesR 
    integer :: i,j,lim1,lim2,lim3,lim4,lim5
    
    NumActvCell   = addedNCPair
    
    NNApclL       = (Hlf_Ncell-NCP_CrtclApSrfc)*(1+NAEC_Apcl) + (NCP_CrtclApSrfc)*(1+NAEC_ApclCrtcl) + 1
    NNApclR       = NNApclL
    
    if (NCP_CrtclApSrfc.gt.NumActvCell) then
       InactvNodesL  = NNApclL - (NumActvCell)*(1+NAEC_ApclCrtcl)
    elseif (NCP_CrtclApSrfc.le.NumActvCell) then
       InactvNodesL  = NNApclL - (NumActvCell-NCP_CrtclApSrfc)*(1+NAEC_Apcl)-(NCP_CrtclApSrfc)*(1+NAEC_ApclCrtcl)
    endif
    
    InactvNodesR  = InactvNodesL
    
    lim1 = InactvNodesL
    lim2 = NNApclL
    lim3 = NNApclL + InactvNodesR
    lim4 = NNApclL + NNApclR
    lim5 = numApclNodes
    
    do i = 1,numApclNodes  
       if (i.le.lim1) then
          activtnFctrApcl(i) = 0
       elseif (i.gt.lim1 .and. i.le.lim2) then
          activtnFctrApcl(i) = 1
       elseif (i.gt.lim2 .and. i.le.lim3) then
          activtnFctrApcl(i) = 0
       elseif (i.gt.lim3 .and. i.le.lim4) then
          activtnFctrApcl(i) = 1
       elseif (i.gt.lim4 .and. i.le.lim5) then
          activtnFctrApcl(i) = 1   
       endif
       write(*,*) activtnFctrApcl(i),i,apclNodes(i),"chkActv"
    enddo
    
  end subroutine readjst_get_actvtn_fctrs_crtclDesgnd_apclSurfc
  
  
  
  
  subroutine store_DS_trnsfrms
    implicit none
    !Renaming purpose, DS=Demolished Spring,TN=Terminal Node
    
    call store_TN_trnsfrms
    
  end subroutine store_DS_trnsfrms
  
  subroutine deallocate_and_reallocate_DS_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = max_node_sprS
    max_area_spr  = max_area_sprS
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS 
    max_node_area = max_node_areaS 
    max_spr_area  = max_spr_areaS
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
  end subroutine deallocate_and_reallocate_DS_trnsfrmVars
  
  
  subroutine deallocate_and_reallocate_DS_trnsfrmVars_NI_system
    implicit none
    integer :: i,j
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = max_node_sprS2
    max_area_spr  = max_area_sprS2
    max_spr_node  = max_spr_nodeS2
    max_area_node = max_area_nodeS2 
    max_node_area = max_node_areaS2 
    max_spr_area  = max_spr_areaS2
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
    write(*,*) max_node_spr,max_node_sprS2,"mns"
    
    open(unit=433,file='spr_node_dealloc_and_realloc_chk.dat')
    do i = 1,N_spr
       write(433,*) spr_node(i,0:max_node_spr),"spr-node"
    enddo
    close(433)
    
  end subroutine deallocate_and_reallocate_DS_trnsfrmVars_NI_system
  
  
  subroutine deallocate_and_reallocate_adding_node_inIC_ApclMem_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = max_node_sprS+1
    max_area_spr  = max_area_sprS
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS 
    max_node_area = max_node_areaS+1
    max_spr_area  = max_spr_areaS
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
  end subroutine deallocate_and_reallocate_adding_node_inIC_ApclMem_trnsfrmVars
  
  
  subroutine deallocate_and_reallocate_combining_cntrl_and_neighApcl_sprs
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = max_node_sprS-1
    max_area_spr  = max_area_sprS
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS 
    max_node_area = max_node_areaS-1
    max_spr_area  = max_spr_areaS
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
  end subroutine deallocate_and_reallocate_combining_cntrl_and_neighApcl_sprs
  
  subroutine deallocate_and_reallocate_VF_region_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    write(*,*) max_node_sprS,max_area_sprS,max_spr_nodeS,max_area_nodeS,max_node_areaS,max_spr_areaS,"max_N/A/S_N/A/S"
    
    max_node_spr  = max_node_sprS
    max_area_spr  = max_area_sprS
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS 
    max_node_area = max(max_node_areaS,(2*(1+addedNCPair)+(1+(2*addedNCPair))*(NAEC_Apcl)))
    max_spr_area  = max(max_spr_areaS,(1*(NAEC_Apcl+1)+(2*addedNCPair)*(NAEC_Apcl+1)))
    
    write(*,*) max_node_area,max_spr_area,"max_node_area/max_spr_area"
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
  end subroutine deallocate_and_reallocate_VF_region_trnsfrmVars
  
  
  subroutine deallocate_and_reallocate_crtclDesgnd_apclSurfc_trnsfrmVars
    implicit none
    integer :: max_node_areaT1,max_node_areaT2,max_node_areaT3
    integer :: max_spr_areaT1,max_spr_areaT2,max_spr_areaT3
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    write(*,*) max_node_sprS,max_area_sprS,max_spr_nodeS,max_area_nodeS,max_node_areaS,max_spr_areaS,"max_N/A/S_N/A/S"
    
    max_node_spr  = max_node_sprS
    max_area_spr  = max_area_sprS
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS
    
    call get_max_node_and_spr_inAn_areaOfPrvSys(max_node_areaS,max_spr_areaS)
    write(*,*) max_node_areaS,max_spr_areaS,"max_node_areaS,max_spr_areaS value using routine" 
    
    max_node_areaT1 = max_node_areaS
    max_node_areaT2 = 5+(NAEC_ApclCrtcl+NAEC_Bsal+2*NAEC_Ltrl)
    max_spr_areaT1  = max_spr_areaS
    max_spr_areaT2  = (NAEC_ApclCrtcl+1+NAEC_Bsal+1+2*(NAEC_Ltrl+1))
    
    if (NCP_CrtclApSrfc.gt.addedNCPair) then
       max_node_areaT3 = 2*(1+addedNCPair) + (1+(2*addedNCPair))*(NAEC_Apcl) + &
            ((addedNCPair*2)+1)*(NAEC_ApclCrtcl-NAEC_Apcl) 
       max_spr_areaT3 = 1*(NAEC_Apcl+1) + (2*addedNCPair)*(NAEC_Apcl+1) + &
            ((2*addedNCPair)+1)*(NAEC_ApclCrtcl-NAEC_Apcl)
       
    elseif (NCP_CrtclApSrfc.le.addedNCPair) then
       max_node_areaT3 = 2*(1+addedNCPair) + (1+(2*addedNCPair))*(NAEC_Apcl) + &
            ((NCP_CrtclApSrfc*2)+1)*(NAEC_ApclCrtcl-NAEC_Apcl)
       max_spr_areaT3 = 1*(NAEC_Apcl+1) + (2*addedNCPair)*(NAEC_Apcl+1) + &
         ((2*NCP_CrtclApSrfc)+1)*(NAEC_ApclCrtcl-NAEC_Apcl)
    endif
    
    max_node_area = max(max_node_areaT1,max_node_areaT2,max_node_areaT3)
    max_spr_area  = max(max_spr_areaT1,max_spr_areaT2,max_spr_areaT3)
    
    write(*,*) max_node_area,max_spr_area,"max_node_area/max_spr_area aft comprsn"
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))    
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    spr_node  = -1 ; spr_area  = -1 
    node_spr  = -1 ; node_area = -1
    area_node = -1 ; area_spr  = -1
    
    Nlist     = -1
    
  end subroutine deallocate_and_reallocate_crtclDesgnd_apclSurfc_trnsfrmVars
  
  subroutine get_max_node_and_spr_inAn_areaOfPrvSys(max_node_areaV,max_spr_areaV)
    implicit none
    integer, intent(out) :: max_node_areaV,max_spr_areaV
    integer              :: i,j
    integer              :: max_node_areaVCurr,max_spr_areaVCurr
    
    max_node_areaV=0     ; max_spr_areaV=0
    max_node_areaVCurr=0 ; max_spr_areaVCurr=0
    
    do i = 1,N_cellS
       max_node_areaVCurr = area_nodeS(i,0)
       max_spr_areaVCurr  = area_sprS(i,0)
       
       if (max_node_areaVCurr.gt.max_node_areaV) then
          max_node_areaV = max_node_areaVCurr
       endif
       
       if (max_spr_areaVCurr.gt.max_spr_areaV) then
          max_spr_areaV = max_spr_areaVCurr
       endif
    enddo
    
  end subroutine get_max_node_and_spr_inAn_areaOfPrvSys
  
  subroutine readjst_all_trnsfrms
    implicit none
    
    call readjst_node_to_other_trnsfrms
    call readjst_spr_to_other_trnsfrms
    call readjst_area_to_other_trnsfrms
    call get_Nlist
    call print_all_trnsfrms
    
  end subroutine readjst_all_trnsfrms
  
  subroutine readjst_node_to_other_trnsfrms
    implicit none
    integer :: i
    integer :: lft_nodeAct,rght_nodeAct 
    integer :: lim1,lim2,lim3,lim4
    integer :: dn11,dn12
    
    lft_nodeAct  = lft_node-2
    rght_nodeAct = rght_node-2
    
    lim1 = lft_nodeAct
    lim2 = lft_node
    lim3 = lft_node+rght_nodeAct
    lim4 = lft_node+rght_node
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:2),i,"nodes,i"
    enddo
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          
          if ((i-lim1)==1) then
             
             dn11 = double_node(1,1)
             dn12 = double_node(1,2)
             
             write(*,*) i,dn11,dn12,"i,dn11,dn12, if i.ne.dn11, stop next line"
             !if (i.ne.dn11) stop
             
             node_spr(i,0) = 4
             node_spr(i,1) = node_sprS(dn11,1)
             node_spr(i,2) = node_sprS(dn11,2)
             node_spr(i,3) = node_sprS(dn12,1)
             node_spr(i,4) = node_sprS(dn12,2)
             
             node_area(i,0) = 3
             node_area(i,1) = node_areaS(dn11,1)
             node_area(i,2) = node_areaS(dn12,1)
             node_area(i,3) = node_areaS(dn11,2)
             
          elseif ((i-lim1)==2) then
             node_spr(i,0:2) = node_sprS(i,0:2)
             node_spr(i,3)   = node_sprS(i,3)-1 !49=50-1
             
             node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
             
          endif 
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          
          if ((i-lim3)==1) then
             write(*,*) i,dn11,dn12,"i,dn11,dn12, if i.ne.dn12, stop next line"
             !if (i.ne.dn12) stop
             
             node_spr(i,0:max_spr_node)   = node_spr(dn11,0:max_spr_node)
             node_area(i,0:max_area_node) = node_area(dn11,0:max_area_node)
             
          elseif ((i-lim3)==2) then
             node_spr(i,0:2) = node_sprS(i,0:2)
             node_spr(i,3)   = node_sprS(i,3)-1 !49=50-1
             
             node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
             
          endif
          
          
       elseif (i.gt.lim4) then
          
          if (i.ne.N_node) write(*,*) i,N_node,"fl:redefine_sys,sb:readjst_all"
          if (i.ne.N_node) stop
          
          node_spr(i,0)   = 4
          node_spr(i,1:2) = 3*ncl-2 !3*8-2 
          node_spr(i,3:4) = 3*(ncl+ncr)-2 !3*16-2
          
          node_area(i,0) = 2
          node_area(i,1) = ncl
          node_area(i,2) = ncl+ncr
          
       endif   
       
    enddo
    
  end subroutine readjst_node_to_other_trnsfrms
  
  
  subroutine readjst_spr_to_other_trnsfrms
    implicit none
    integer :: i
    integer :: nslAct,nsrAct 
    integer :: lim1,lim2,lim3,lim4
    
    nslAct = nsl-nsecl !I am not changing presuming nothing changes
    nsrAct = nsr-nsecr
    
    lim1 = nslAct
    lim2 = nsl
    lim3 = nsl+nsrAct
    lim4 = nsl+nsr
    
    do i = 1,N_spr
       
       if (i.le.lim1) then
          spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
          spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then

          if ((i-lim1)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_nodeS+1 !origin
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim1)==2) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr) 
             
          elseif ((i-lim1)==3) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          endif
             
       elseif (i.gt.lim2 .and. i.le.lim3) then
          spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
          spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          
          if ((i-lim3)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_nodeS+1 !origin
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim3)==2) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim3)==3) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          endif
          
       elseif (i.gt.lim4) then 
          spr_node(i,0:max_node_spr) = spr_nodeS(i+1,0:max_node_spr) !49 & 50
          spr_area(i,0:max_area_spr) = spr_areaS(i+1,0:max_area_spr) !49 & 50
          
       endif
       
    enddo
    
  end subroutine readjst_spr_to_other_trnsfrms
  
  
  subroutine readjst_area_to_other_trnsfrms
    implicit none
    integer :: i
    integer :: nclAct,ncrAct 
    integer :: lim1,lim2,lim3,lim4
    
    nclAct = ncl-1
    ncrAct = ncr-1
    
    lim1 = nclAct
    lim2 = ncl
    lim3 = ncl+ncrAct
    lim4 = ncl+ncr
    write(*,*) lim1,lim2,lim3,lim4,"lims"
    !call sleep(10)
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          area_node(i,0:max_node_area) = area_nodeS(i,0:max_node_area)
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          area_node(i,0)   = 5
          area_node(i,1:4) = area_nodeS(i,1:4)
          area_node(i,5)   = N_nodeS+1

          area_spr(i,0:max_spr_area) = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          area_node(i,0:max_node_area) = area_nodeS(i,0:max_node_area)
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          area_node(i,0)   = 5
          area_node(i,1)   = area_nodeS(i,1)
          area_node(i,2)   = N_nodeS+1
          area_node(i,3:5) = area_nodeS(i,2:4)
          
          area_spr(i,0:max_spr_area)  = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim4) then
          area_node(i,0)   = 3
          area_node(i,1:3) = area_nodeS(i,2:4)
          
          area_spr(i,0) = 3
          area_spr(i,1) = area_sprS(i,1)
          area_spr(i,2) = area_sprS(i,2)
          area_spr(i,3) = area_sprS(i,4)
          
          if (i.ne.N_cell) write(*,*) i,N_cell,"i=/N_cell"
          if (i.ne.N_cell) stop
          
       endif
       
    enddo
    
  end subroutine readjst_area_to_other_trnsfrms
  
  
  subroutine readjst_all_trnsfrms_wo_dmlsh
    implicit none
    
    call readjst_node_to_other_trnsfrms_wo_dmlsh
    call readjst_spr_to_other_trnsfrms_wo_dmlsh
    call readjst_area_to_other_trnsfrms_wo_dmlsh
    call get_Nlist
    call print_all_trnsfrms
    
  end subroutine readjst_all_trnsfrms_wo_dmlsh
  
  subroutine readjst_node_to_other_trnsfrms_wo_dmlsh
    implicit none
    integer :: i
    integer :: lft_nodeAct,lft_nodePrvAct
    integer :: rght_nodeAct,rght_nodePrvAct
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: dn11,dn12

    open(unit=142,file='readjst_lims.dat',position='append')
    
    lft_nodeAct  = lft_node-CellsMeet*2  ; lft_nodePrvAct  = lft_node-(CellsMeet-1)*2
    rght_nodeAct = rght_node-CellsMeet*2 ; rght_nodePrvAct = rght_node-(CellsMeet-1)*2
    
    lim1 = lft_nodeAct !14 to 12
    lim2 = lft_nodePrvAct !16 to 14
    lim3 = lft_node !18
    lim4 = lft_node+rght_nodeAct
    lim5 = lft_node+rght_nodePrvAct
    lim6 = lft_node+rght_node
    
    
    write(142,*)lft_nodeAct,lft_nodePrvAct,"lft"
    write(142,*)rght_nodeAct,rght_nodePrvAct,"rght"
    write(142,*)lft_node,rght_node,"LRN"
    
    write(142,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims"
    write(142,*) CellsMeet,"CM "
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims"
    write(*,*) CellsMeet,"CM "
    
    close(142)
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          
          if ((i-lim1)==1) then
             dn11 = double_node(1,1)
             dn12 = double_node(1,2)
             write(*,*) i,dn11,dn12,"i,dn11,dn12, if i.ne.dn11, stop next line"
             if (i.ne.dn11) stop
             
             node_spr(i,0) = 6
             node_spr(i,1) = node_sprS(dn11,1)
             node_spr(i,2) = node_sprS(dn11,2)
             node_spr(i,3) = node_sprS(dn12,1)
             node_spr(i,4) = node_sprS(dn12,2)
             node_spr(i,5) = node_sprS(dn11,3)
             node_spr(i,6) = node_sprS(dn12,3)
             
             node_area(i,0) = 4
             node_area(i,1) = node_areaS(dn11,1)
             node_area(i,2) = node_areaS(dn12,1)
             node_area(i,3) = node_areaS(dn11,2)
             node_area(i,4) = node_areaS(dn12,2)
             
          elseif ((i-lim1)==2) then
             node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
             node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
             
          endif
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim4 .and. i.le.lim5) then
          
          if ((i-lim4)==1) then
             write(*,*) i,dn11,dn12,"i,dn11,dn12, if i.ne.dn12, stop next line"
             if (i.ne.dn12) stop
             
             node_spr(i,0:max_spr_node)   = node_spr(dn11,0:max_spr_node)
             node_area(i,0:max_area_node) = node_area(dn11,0:max_area_node)
          elseif ((i-lim4)==2) then
             node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
             node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
             
          endif
          
       elseif (i.gt.lim5 .and. i.le.lim6) then
          node_spr(i,0:max_spr_node)   = node_sprS(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaS(i,0:max_area_node)
          
       elseif (i.gt.lim6) then
          node_spr(i,0) = node_sprS(i,0)
          node_spr(i,1) = node_sprS(i,1)-nsecl
          node_spr(i,2) = node_sprS(i,2)-nsecl
          node_spr(i,3) = node_sprS(i,3)-nsecr
          node_spr(i,4) = node_sprS(i,4)-nsecr
          
          node_area(i,0) = node_areaS(i,0)
          node_area(i,1) = node_areaS(i,1)-1
          node_area(i,2) = node_areaS(i,2)-1
       endif
       
    enddo
    
  end subroutine readjst_node_to_other_trnsfrms_wo_dmlsh
  
  subroutine readjst_spr_to_other_trnsfrms_wo_dmlsh
    implicit none
    integer :: i
    integer :: nslAct,nslPrvAct
    integer :: nsrAct,nsrPrvAct 
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    
    nslAct = nsl-CellsMeet*nsecl ; nslPrvAct = nsl-(CellsMeet-1)*nsecl
    nsrAct = nsr-CellsMeet*nsecr ; nsrPrvAct = nsr-(CellsMeet-1)*nsecr
    
    lim1 = nslAct !18 to 15
    lim2 = nslPrvAct !21 to 18
    lim3 = nsl !24
    lim4 = nsl+nsrAct
    lim5 = nsl+nsrPrvAct
    lim6 = nsl+nsr
    
    
    do i = 1,N_spr
       
       if (i.le.lim1) then !1-15
          spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
          spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then !16-18
          
          if ((i-lim1)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_node !37
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim1)==2 .or. (i-lim1)==3) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          endif
          
       elseif (i.gt.lim2 .and. i.le.lim3) then !19-24
          
          if ((i-lim2)==1) then
             spr_node(i,0) = 2
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = spr_nodeS(i,3)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim2).gt.1) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          endif
          
       elseif (i.gt.lim3 .and. i.le.lim4) then !25-39
          spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
          spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim4 .and. i.le.lim5) then !40-42
          
          if ((i-lim4)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_node
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim4)==2 .or. (i-lim4)==3) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          endif
          
       elseif (i.gt.lim5 .and. i.le.lim6) then !43-48
          
          if ((i-lim5)==1) then
             spr_node(i,0) = 2
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = spr_nodeS(i,3)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim5).gt.1) then
             spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          endif
          
       elseif (i.gt.lim6) then
          spr_node(i,0:max_node_spr) = spr_nodeS(i,0:max_node_spr)
          spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
       endif
       
    enddo
    
  end subroutine readjst_spr_to_other_trnsfrms_wo_dmlsh
  
  
  subroutine readjst_area_to_other_trnsfrms_wo_dmlsh
    implicit none
    integer :: i
    integer :: nclAct,nclPrvAct
    integer :: ncrAct,ncrPrvAct 
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    
    nclAct = ncl-CellsMeet ; nclPrvAct = ncl-(CellsMeet-1)
    ncrAct = ncr-CellsMeet ; ncrPrvAct = ncr-(CellsMeet-1)
    
    lim1 = nclAct !5
    lim2 = nclPrvAct !6
    lim3 = ncl !8
    lim4 = ncl+ncrAct
    lim5 = ncl+ncrPrvAct
    lim6 = ncl+ncr
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims in Area_to_other_adjst"
    
    do i = 1,N_cell
       
       if (i.le.lim1) then !1-5
          area_node(i,0:max_node_area) = area_nodeS(i,0:max_node_area)
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then !6 
          area_node(i,0) = 5
          area_node(i,1:4) = area_nodeS(i,1:4)
          area_node(i,5)   = N_node
          
          area_spr(i,0:max_spr_area) = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim2 .and. i.le.lim3) then !7-8
          area_node(i,0)   = 4
          area_node(i,1:4) = area_nodeS(i,1:4)
          
          area_spr(i,0:max_spr_area) = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then !9-13
          area_node(i,0:max_node_area) = area_nodeS(i,0:max_node_area)
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim4 .and. i.le.lim5) then !14
          area_node(i,0)   = 5
          area_node(i,1)   = area_nodeS(i,1)
          area_node(i,2)   = N_node
          area_node(i,3:5) = area_nodeS(i,2:4)
          
          area_spr(i,0:max_spr_area) = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim5 .and. i.le.lim6) then !15-16

          if ((i-lim5)==1) then
             area_node(i,0)   = 4
             area_node(i,1)   = area_nodeS(i,1)
             area_node(i,2:4) = area_nodeS(i,3:5)
          elseif ((i-lim5).gt.1) then
             area_node(i,0)   = 4
             area_node(i,1:4) = area_nodeS(i,1:4)
          endif
          
          area_spr(i,0:max_spr_area)  = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim6) then
          area_node(i,0:max_node_area) = area_nodeS(i,0:max_node_area)
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       endif
       
    enddo
    
    
  end subroutine readjst_area_to_other_trnsfrms_wo_dmlsh
  
  
  subroutine readjst_all_trnsfrms_inIC_ApclMem
    implicit none
    
    call readjst_node_to_other_trnsfrms_inIC_ApclMem
    call readjst_spr_to_other_trnsfrms_inIC_ApclMem
    call readjst_area_to_other_trnsfrms_inIC_ApclMem
    call print_all_trnsfrms ; write(*,*) "MADE UPTO THIS"
    call get_Nlist
    call print_all_trnsfrms
    
  end subroutine readjst_all_trnsfrms_inIC_ApclMem
  
  subroutine readjst_node_to_other_trnsfrms_inIC_ApclMem
    implicit none
    integer :: i,j,imax,jmax
    
    write(*,*) N_nodeS,N_node,"N_nodes bfr aft"
    write(*,*) node_spr(N_nodeS+1,0:max_spr_node),"node_spr for N_node+1"
    write(*,*) node_spr(N_nodeS+2,0:max_spr_node),"node_spr for N_node+2"
    
    write(*,*) node_area(N_nodeS+1,0:max_area_node),"node_areaFor N_node+1"
    write(*,*) node_area(N_nodeS+2,0:max_area_node),"node_areaFor N_node+2"
    
    node_spr(1:N_nodeS,0:max_spr_node) = node_sprS(1:N_nodeS,0:max_spr_node)
    
    node_spr(N_nodeS+1,0)   = 4
    node_spr(N_nodeS+1,1:2) = 3*ncl-2
    node_spr(N_nodeS+1,3:4) = 3*(ncl+ncr)+1
    
    node_spr(N_nodeS+2,0)   = 4
    node_spr(N_nodeS+2,1:2) = 3*(ncl+ncr)-2
    node_spr(N_nodeS+2,3:4) = 3*(ncl+ncr)+1
    
    node_area(1:N_nodeS,0:max_area_node) = node_areaS(1:N_nodeS,0:max_area_node)
    
    node_area(N_nodeS+1,0) = 2
    node_area(N_nodeS+1,1) = ncl
    node_area(N_nodeS+1,2) = N_cell
    
    node_area(N_nodeS+2,0) = 2
    node_area(N_nodeS+2,1) = ncl+ncr
    node_area(N_nodeS+2,2) = N_cell
    
    open(unit=443,file='readjst_NTO_inIC_ApclMem.dat')
    
    imax = 2
    jmax = N_node
    
    do i = 1,imax
       
       jmax = N_node
       
       do j = 1,jmax
          if (i==1) write(443,*) j,node_spr(j,0:max_spr_node)
          if (i==2) write(443,*) j,node_area(j,0:max_area_node)
       enddo
       write(443,*) " "
    enddo
    
    close(443)
    
  end subroutine readjst_node_to_other_trnsfrms_inIC_ApclMem
  
  
  subroutine readjst_spr_to_other_trnsfrms_inIC_ApclMem
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: i,j,imax,jmax
    
    lim1 = (Hlf_Ncell-1)*3         ! unchanged left [1 ~ 10*3=30]
    lim2 = lim1 + 3                ! total left spr [31 ~ 33]
    lim3 = lim2 + (Hlf_Ncell-1)*3  ! unchaned right [34 ~ 63]
    lim4 = lim3 + 3                ! total right [64 ~ 66]
    lim5 = lim4 + 2                ! two central cells
    
    write(*,*) lim1,lim2,lim3,lim4,lim5, "lims in "
    write(*,*) N_spr,N_sprS,"N_spr"
    
    !IMPORTANT Explanation: In this system I need the max_node_spr=4 which was =3 in previous cases, that is why
    !whenever I was copying spr_node(i,0:~) = spr_nodeS(i,0:~), I was writing (max_node_spr-1) in the place of ~
    
    do i = 1,N_spr
       
       if (i.le.lim1) then
          spr_node(i,0:max_node_spr-1) = spr_nodeS(i,0:max_node_spr-1)
          spr_area(i,0:max_area_spr)   = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          
          if ((i-lim1)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_nodeS+1
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim1) .gt. 1) then
             spr_node(i,0:max_node_spr-1) = spr_nodeS(i,0:max_node_spr-1)
             spr_area(i,0:max_area_spr)   = spr_areaS(i,0:max_area_spr) 
          endif
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          spr_node(i,0:max_node_spr-1) = spr_nodeS(i,0:max_node_spr-1)
          spr_area(i,0:max_area_spr)   = spr_areaS(i,0:max_area_spr)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          
          if ((i-lim3)==1) then
             spr_node(i,0) = 3
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_nodeS+2
             spr_node(i,3) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim3) .gt. 1) then
             spr_node(i,0:max_node_spr-1) = spr_nodeS(i,0:max_node_spr-1)
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
          endif
          
       elseif (i.gt.lim4 .and. i.le.lim5) then ! central springs 
          
          if ((i-lim4)==1) then
             spr_node(i,0) = 4
             spr_node(i,1) = spr_nodeS(i,1)
             spr_node(i,2) = N_nodeS+1
             spr_node(i,3) = N_nodeS+2
             spr_node(i,4) = spr_nodeS(i,2)
             
             spr_area(i,0:max_area_spr) = spr_areaS(i,0:max_area_spr)
             
          elseif ((i-lim4)==2) then
             spr_node(i,0:max_node_spr-1) = spr_nodeS(i,0:max_node_spr-1)
             spr_area(i,0:max_area_spr)   = spr_areaS(i,0:max_area_spr)
          endif
          
       endif
       
    enddo

    
    open(unit=445,file='readjst_STO_inIC_ApclMem.dat')
    
    imax = 2
    jmax = N_spr
    
    do i = 1,imax
       
       jmax = N_spr
       
       do j = 1,jmax
          if (i==1) write(445,*) j,spr_node(j,0:max_node_spr)
          if (i==2) write(445,*) j,spr_area(j,0:max_area_spr)
       enddo
       write(445,*) " "
    enddo
    
    close(445)
    
  end subroutine readjst_spr_to_other_trnsfrms_inIC_ApclMem
  
  subroutine readjst_area_to_other_trnsfrms_inIC_ApclMem
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: nodesInAreaBfr,nodesInAreaAft
    integer :: i,j,imax,jmax
    
    lim1 = Hlf_Ncell-1
    lim2 = Hlf_Ncell
    lim3 = lim2+(Hlf_Ncell-1)
    lim4 = lim3+1
    lim5 = N_cell

    write(*,*) lim1,lim2,lim3,lim4,lim5,"lims in readjst_ATO_trnsfr_inIC_ApclMem"
    call sleep(1)
    
    !IMPORTANT Explntn: In this systm I need the max_node_area=6 which was =5 in previous cases, that is why
    !whenever I was copying area_node(i,0:~) = area_nodeS(i,0:~), I was writing (max_node_area-1) in place of ~
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          area_node(i,0:max_node_area-1) = area_nodeS(i,0:max_node_area-1)
          area_spr(i,0:max_spr_area)     = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim1 .and. i.le.lim2) then
          
          nodesInAreaBfr = area_nodeS(i,0)
          nodesInAreaAft = nodesInAreaBfr+1
          
          area_node(i,0) = nodesInAreaAft
          area_node(i,1:nodesInAreaBfr) = area_nodeS(i,1:nodesInAreaBfr)
          area_node(i,nodesInAreaAft)   = N_nodeS+1
          
          area_spr(i,0:max_spr_area) = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim2 .and. i.le.lim3) then
          
          area_node(i,0:max_node_area-1) = area_nodeS(i,0:max_node_area-1)
          area_spr(i,0:max_spr_area)     = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim3 .and. i.le.lim4) then
          
          nodesInAreaBfr = area_nodeS(i,0)
          nodesInAreaAft = nodesInAreaBfr+1
          
          area_node(i,0)                = nodesInAreaAft
          area_node(i,1)                = area_nodeS(i,1)
          area_node(i,2)                = N_nodeS+2
          area_node(i,3:nodesInAreaAft) = area_nodeS(i,2:nodesInAreaBfr)
          
          area_spr(i,0:max_spr_area)   = area_sprS(i,0:max_spr_area)
          
       elseif (i.gt.lim4 .and. i.le.lim5) then
          
          nodesInAreaBfr = area_nodeS(i,0)
          nodesInAreaAft = nodesInAreaBfr+2
          
          area_node(i,0)                = nodesInAreaAft
          area_node(i,1)                = area_nodeS(i,1)
          area_node(i,2)                = N_nodeS+2
          area_node(i,3)                = N_nodeS+1
          area_node(i,4:nodesInAreaAft) = area_nodeS(i,2:nodesInAreaBfr)
          
          area_spr(i,0:max_spr_area)    = area_sprS(i,0:max_spr_area)
          
       endif
       
    enddo
    
    open(unit=446,file='readjst_ATO_inIC_ApclMem.dat')
    
    imax = 2
    jmax = N_cell
    
    do i = 1,imax
       
       jmax = N_cell
       
       do j = 1,jmax
          if (i==1) write(446,*) j,area_node(j,0:max_node_area)
          if (i==2) write(446,*) j,area_spr(j,0:max_spr_area)
       enddo
       write(446,*) " "
    enddo
    
    close(446)
    
  end subroutine readjst_area_to_other_trnsfrms_inIC_ApclMem
  
  
  
  subroutine readjst_all_trnsfrms_VF_region
    implicit none
    
    call get_necessary_variabls_for_Editingtrnsfrm
    call readjst_node_to_other_trnsfrmsVF_region
    call readjst_spr_to_other_trnsfrmsVF_region
    call readjst_area_to_other_trnsfrmsVF_region
    call print_all_trnsfrms ; write(*,*) "MADE UPTO THIS VF"
    call get_Nlist
    call print_all_trnsfrms
    
  end subroutine readjst_all_trnsfrms_VF_region
  
  
  
  subroutine get_necessary_variabls_for_Editingtrnsfrm
    implicit none
    
    Hlf_Ncell = (N_cell-1)/2 ; write(*,*) Hlf_Ncell,"Hlf_Ncell-trnsfrmVF"
    
    VF_BotmSurfcNodes(1) =   (Hlf_Ncell+1)*2 - 1
    VF_BotmSurfcNodes(2) = 2*(Hlf_Ncell+1)*2 - 1
    
    VF_UpprSurfcNodes(1) =   ((Hlf_Ncell+1)-addedNCPair)*2 - 1
    VF_UpprSurfcNodes(2) = (2*(Hlf_Ncell+1)-addedNCPair)*2 - 1
    
    write(*,*) VF_BotmSurfcNodes(1:2),"VFb"
    write(*,*) VF_UpprSurfcNodes(1:2),"VFu"
    
    numINperCell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
    
  end subroutine get_necessary_variabls_for_Editingtrnsfrm
  
  subroutine readjst_node_to_other_trnsfrmsVF_region
    implicit none
    
    if (NI_incldAS_woPP==0) call readjst_node_to_other_trnsfrmsVF_region_methd1
    if (NI_incldAS_woPP==1) call readjst_node_to_other_trnsfrmsVF_region_methd2
    
  end subroutine readjst_node_to_other_trnsfrmsVF_region
  
  subroutine readjst_node_to_other_trnsfrmsVF_region_methd1
    implicit none
    write(*,*) N_nodeSNVFR,N_node,max_spr_node,max_spr_nodeS,"[[ 1"
    
    node_spr(1:N_nodeSNVFR,0:max_spr_node)=node_sprS(1:N_nodeSNVFR,0:max_spr_node)
    node_area(1:N_nodeSNVFR,0:max_area_node)=node_areaS(1:N_nodeSNVFR,0:max_area_node)
    
    node_spr(VF_UpprSurfcNodes(1),0) = 4
    node_spr(VF_UpprSurfcNodes(1),4) = N_spr
    node_spr(VF_UpprSurfcNodes(2),0) = 4
    node_spr(VF_UpprSurfcNodes(2),4) = N_spr
    
    write(*,*) node_spr(VF_UpprSurfcNodes(1),0:max_spr_node),"NS VFu1 m1"
    write(*,*) node_spr(VF_UpprSurfcNodes(2),0:max_spr_node),"NS VFu2 m1"
    
    node_area(VF_UpprSurfcNodes(1),0) = 3 ; node_area(VF_UpprSurfcNodes(1),3) = N_cell
    node_area(VF_UpprSurfcNodes(2),0) = 3 ; node_area(VF_UpprSurfcNodes(2),3) = N_cell
    node_area(VF_BotmSurfcNodes(1),0) = 3 ; node_area(VF_BotmSurfcNodes(1),3) = N_cell
    node_area(VF_BotmSurfcNodes(2),0) = 3 ; node_area(VF_BotmSurfcNodes(2),3) = N_cell
    
    write(*,*)node_area(VF_UpprSurfcNodes(1),0:max_area_node),"NAVFu1 m1"
    write(*,*)node_area(VF_UpprSurfcNodes(2),0:max_area_node),"NAVFu2 m1"
    write(*,*)node_area(VF_BotmSurfcNodes(1),0:max_area_node),"NAVFb1 m1"
    write(*,*)node_area(VF_BotmSurfcNodes(2),0:max_area_node),"NAVFb2 m1"
    
  end subroutine readjst_node_to_other_trnsfrmsVF_region_methd1
  
  
  subroutine readjst_node_to_other_trnsfrmsVF_region_methd2
    implicit none
    integer :: j,jmax,k,kmax
    integer :: nodeNm,jmaxHlf
    integer :: lftNode,rghtNode
    
    write(*,*) N_nodeSNVFR,N_node,max_spr_node,max_spr_nodeS,"[[ 2"
    
    node_spr(1:N_nodeSNVFR,0:max_spr_node)=node_sprS(1:N_nodeSNVFR,0:max_spr_node)
    node_area(1:N_nodeSNVFR,0:max_area_node)=node_areaS(1:N_nodeSNVFR,0:max_area_node)
    
    write(*,*) node_spr(VF_UpprSurfcNodes(1),0:max_spr_node),"NS VFu1 m2"
    write(*,*) node_spr(VF_UpprSurfcNodes(2),0:max_spr_node),"NS VFu2 m2"
    
    jmax = 1+addedNCPair ; write(*,*) jmax,addedNCPair,"jmax-addedNCPair"
    
    do j = 1,jmax
       
       lftNode  = (2*Hlf_Ncell)+1 - (j-1)*2
       rghtNode = (2*(Hlf_Ncell+1))*2-1 - (j-1)*2
       
       node_area(lftNode,0)  = 3 ; node_area(lftNode,3)  = N_cell
       node_area(rghtNode,0) = 3 ; node_area(rghtNode,3) = N_cell
       
       write(*,*) node_area(lftNode,0:max_area_node),lftNode,j,"VF node-areaL"
       write(*,*) node_area(rghtNode,0:max_area_node),rghtNode,j,"VF node-areaR"
       
    enddo
    
    jmax        = 1+2*addedNCPair
    kmax        = NAEC_Apcl
    jmaxHlf     = jmax/2
    
    do j = 1,jmax 
       kmax = NAEC_Apcl
       
       do k = 1,kmax
          
          if (j.le.jmaxHlf)     nodeNm = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-jmaxHlf+(j-1))*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k      ! 48 + (11-jmaxHlf+(j-1))*6 + k 
          if (j == (jmaxHlf+1)) nodeNm = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell*2)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k
          if (j.gt.(jmaxHlf+1)) nodeNm = 2*(Hlf_Ncell+1)*2+((Hlf_Ncell*2)-(j-jmaxHlf-1))*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k  ! 48 + (22-(j-2-1))*6 + k   
          
          node_area(nodeNm,0) = 2
          node_area(nodeNm,2) = N_cell
          write(*,*) node_area(nodeNm,0:max_area_node),nodeNm,"NAVF IN m2"
          
       enddo
       
    enddo
    
  end subroutine readjst_node_to_other_trnsfrmsVF_region_methd2
  
  subroutine readjst_spr_to_other_trnsfrmsVF_region
    implicit none
    
    !if (addedNCpair == 2) call readjst_spr_to_other_trnsfrmsVF_region_methd1
    !if (addedNCpair.ne.2) call readjst_spr_to_other_trnsfrmsVF_region_methd2
    call readjst_spr_to_other_trnsfrmsVF_region_methd2
    
  end subroutine readjst_spr_to_other_trnsfrmsVF_region
  
  
  subroutine readjst_spr_to_other_trnsfrmsVF_region_methd1
    implicit none
    integer :: nsprsInACell
    integer :: i,j
    
    write(*,*) "modelID = ",modelID
    if (modelID == 1) write(*,*) "not for modelID=1"
    write(*,*) N_sprSNVFR,N_spr,"N_sprSNVFR,N_spr"
    
    spr_node(1:N_sprSNVFR,0:max_node_spr) = spr_nodeS(1:N_sprSNVFR,0:max_node_spr)
    spr_area(1:N_sprSNVFR,0:max_area_spr) = spr_areaS(1:N_sprSNVFR,0:max_area_spr)
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    write(*,*)Hlf_Ncell,nsprsInACell,"Hlf_Ncell-nsprsInACell in VF spr"
    
    allocate(sprApclCntrl(1:(NAEC_Apcl+1))) ; sprApclCntrl = -1
    allocate(sprApclNC1L(1:(NAEC_Apcl+1)))  ; sprApclNC1L  = -1
    allocate(sprApclNC1R(1:(NAEC_Apcl+1)))  ; sprApclNC1R  = -1
    allocate(sprApclNC2L(1:(NAEC_Apcl+1)))  ; sprApclNC2L  = -1
    allocate(sprApclNC2R(1:(NAEC_Apcl+1)))  ; sprApclNC2R  = -1
    
    do i = 1,(NAEC_Apcl+1) 
       sprApclCntrl(i) = (2*Hlf_Ncell)  *(nsprsInACell) + i
       sprApclNC1L(i)  = (Hlf_Ncell-1)  *(nsprsInACell) + i
       sprApclNC1R(i)  = (2*Hlf_Ncell-1)*(nsprsInACell) + i
       sprApclNC2L(i)  = (Hlf_Ncell-2)  *(nsprsInACell) + i
       sprApclNC2R(i)  = (2*Hlf_Ncell-2)*(nsprsInACell) + i
    enddo
    
    write(*,*) sprApclCntrl(1:(NAEC_Apcl+1)),"sprApclCntrl"
    write(*,*) sprApclNC1L(1:(NAEC_Apcl+1)), "sprApclNC1L"
    write(*,*) sprApclNC1R(1:(NAEC_Apcl+1)), "sprApclNC1R"
    write(*,*) sprApclNC2L(1:(NAEC_Apcl+1)), "sprApclNC2L"
    write(*,*) sprApclNC2R(1:(NAEC_Apcl+1)), "sprApclNC2R"
    
    do i = 1,(NAEC_Apcl+1)
       
       spr_area(sprApclCntrl(i),0) = 2
       spr_area(sprApclCntrl(i),1) = N_cell-1
       spr_area(sprApclCntrl(i),2) = N_cell
       write(*,*) spr_area(sprApclCntrl(i),0:max_area_spr),i,"SA VFarea CC"
       
       spr_area(sprApclNC1L(i),0)  = 2
       spr_area(sprApclNC1L(i),1)  = Hlf_Ncell
       spr_area(sprApclNC1L(i),2)  = N_cell
       write(*,*) spr_area(sprApclNC1L(i),0:max_area_spr),i,"SA VFarea NC1L"
       
       spr_area(sprApclNC1R(i),0)  = 2
       spr_area(sprApclNC1R(i),1)  = 2*Hlf_Ncell
       spr_area(sprApclNC1R(i),2)  = N_cell
       write(*,*) spr_area(sprApclNC1R(i),0:max_area_spr),i,"SA VFarea NC1R"
       
       spr_area(sprApclNC2L(i),0)  = 2
       spr_area(sprApclNC2L(i),1)  = Hlf_Ncell-1
       spr_area(sprApclNC2L(i),2)  = N_cell
       write(*,*) spr_area(sprApclNC1L(i),0:max_area_spr),i,"SA VFarea NC1L"
       
       spr_area(sprApclNC2R(i),0)  = 2
       spr_area(sprApclNC2R(i),1)  = 2*Hlf_Ncell-1
       spr_area(sprApclNC2R(i),2)  = N_cell
       write(*,*) spr_area(sprApclNC1R(i),0:max_area_spr),i,"SA VFarea NC1R"
       
    enddo
    
    !spr_area(sprVFtopLayr,0) = 1
    !spr_area(sprVFtopLayr,1) = N_cell
    !write(*,*) spr_area(sprVFtopLayr,0:max_area_spr),"SA VFarea TopLayr"
    
  end subroutine readjst_spr_to_other_trnsfrmsVF_region_methd1
  
  subroutine readjst_spr_to_other_trnsfrmsVF_region_methd2
    implicit none
    integer :: nsprsInACell
    integer :: i,j
    
    write(*,*) "modelID = ",modelID
    if (modelID == 1) write(*,*) "not for modelID=1"
    write(*,*) N_sprSNVFR,N_spr,"N_sprSNVFR,N_spr"
    
    spr_node(1:N_sprSNVFR,0:max_node_spr) = spr_nodeS(1:N_sprSNVFR,0:max_node_spr)
    spr_area(1:N_sprSNVFR,0:max_area_spr) = spr_areaS(1:N_sprSNVFR,0:max_area_spr)
    
    nsprsInACell=(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3) ;write(*,*)Hlf_Ncell,nsprsInACell,"Hlf_Ncell-nsprInC VF spr"
    
    allocate(sprApclCntrl(1:(NAEC_Apcl+1)))             ; sprApclCntrl = -1
    allocate(sprApclNCL(1:addedNCPair,1:(NAEC_Apcl+1))) ; sprApclNCL   = -1
    allocate(sprApclNCR(1:addedNCPair,1:(NAEC_Apcl+1))) ; sprApclNCR   = -1  
    
    do i = 0,addedNCPair  
       if (i==0) then
          do j = 1,(NAEC_Apcl+1)
              sprApclCntrl(j) = (2*Hlf_Ncell)  *(nsprsInACell) + j
           enddo
           write(*,*) sprApclCntrl(1:(NAEC_Apcl+1)),"sprApclCntrl"
       else
          do j = 1,(NAEC_Apcl+1)
             sprApclNCL(i,j) = (Hlf_Ncell-i)  *(nsprsInACell) + j
             sprApclNCR(i,j) = (2*Hlf_Ncell-i)*(nsprsInACell) + j
          enddo
          write(*,*) sprApclNCL(i,1:(NAEC_Apcl+1)),i,"sprApclNCL"
          write(*,*) sprApclNCR(i,1:(NAEC_Apcl+1)),i,"sprApclNCR"
       endif
    enddo
    
    do i = 0,addedNCPair  
       if (i==0) then
          do j = 1,(NAEC_Apcl+1)
             spr_area(sprApclCntrl(j),0) = 2
             spr_area(sprApclCntrl(j),1) = N_cell-1
             spr_area(sprApclCntrl(j),2) = N_cell
             write(*,*) spr_area(sprApclCntrl(j),0:max_area_spr),j,"SA VFarea CC"
           enddo
       else
          do j = 1,(NAEC_Apcl+1)
             spr_area(sprApclNCL(i,j),0) = 2
             spr_area(sprApclNCL(i,j),1) = Hlf_Ncell-(i-1)
             spr_area(sprApclNCL(i,j),2) = N_cell
             
             spr_area(sprApclNCR(i,j),0) = 2
             spr_area(sprApclNCR(i,j),1) = (2*Hlf_Ncell)-(i-1)
             spr_area(sprApclNCR(i,j),2) = N_cell
             
             write(*,*) spr_area(sprApclNCL(i,j),0:max_area_spr),i,j,"SA VFarea NCL"
             write(*,*) spr_area(sprApclNCR(i,j),0:max_area_spr),i,j,"SA VFarea NCR"
          enddo
       endif
    enddo
    
  end subroutine readjst_spr_to_other_trnsfrmsVF_region_methd2
  
  subroutine readjst_area_to_other_trnsfrmsVF_region
    implicit none
    
    if (NI_incldAS_woPP == 0) call readjst_area_to_other_trnsfrmsVF_region_methd1
    if (NI_incldAS_woPP == 1) call readjst_area_to_other_trnsfrmsVF_region_methd3
    
  end subroutine readjst_area_to_other_trnsfrmsVF_region
  
  subroutine readjst_area_to_other_trnsfrmsVF_region_methd1
    implicit none
    integer :: sprAS,sprBS,sprLS !(spr in Apcl/Basal/Ltrl) side
    
    sprAS = NAEC_Apcl+1 ; sprBS = NAEC_Bsal+1 ; sprLS = NAEC_Ltrl+1
    
    area_node(1:N_cellSNVFR,0:max_node_area) = area_nodeS(1:(N_cell-1),0:max_node_area)
    area_spr(1:N_cellSNVFR,0:max_spr_area)   = area_sprS(1:N_cellSNVFR,0:max_spr_area)
    
    area_node(N_cell,0) = 4
    area_node(N_cell,1) = VF_UpprSurfcNodes(1) 
    area_node(N_cell,2) = VF_BotmSurfcNodes(1)
    area_node(N_cell,3) = VF_BotmSurfcNodes(2)
    area_node(N_cell,4) = VF_UpprSurfcNodes(2)
    
    write(*,*) area_node(N_cell,0:max_node_area),"AN VFcell 1"
    
    area_spr(N_cell,0)                         = 1 + 3*sprAS
    area_spr(N_cell,1:sprAS)                   = sprApclNC1L(1:sprAS)
    area_spr(N_cell,(sprAS+1))                 = sprVFtopLayr
    area_spr(N_cell,(sprAS+2):(2*sprAS+1))     = sprApclCntrl(1:sprAS)
    area_spr(N_cell,(2*sprAS+2):(3*sprAS+1))   = sprApclNC1R(1:sprAS)
    
    write(*,*) area_spr(N_cell,0:max_spr_area),"AS VFcell 1"
    
  end subroutine readjst_area_to_other_trnsfrmsVF_region_methd1
  
  subroutine readjst_area_to_other_trnsfrmsVF_region_methd2
    implicit none
    integer :: k,kmax
    integer :: currIndx
    integer :: sprAS,sprBS,sprLS !(spr in Apcl/Basal/Ltrl) side
    
    sprAS = NAEC_Apcl+1 ; sprBS = NAEC_Bsal+1 ; sprLS = NAEC_Ltrl+1 
    
    area_node(1:N_cellSNVFR,0:max_node_areaS) = area_nodeS(1:(N_cell-1),0:max_node_areaS)
    area_spr(1:N_cellSNVFR,0:max_spr_areaS)   = area_sprS(1:N_cellSNVFR,0:max_spr_areaS)
    
    kmax                = NAEC_Apcl
    
    area_node(N_cell,0) = 2+(NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    
    area_node(N_cell,1) = VF_UpprSurfcNodes(1)
    currIndx            = 1
    
    do k = 1,kmax
       area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-2)*(numINperCell)+k
       currIndx                     = currIndx+1
    enddo
    
    area_node(N_cell,currIndx+1) = VF_UpprSurfcNodes(1) + 2
    currIndx                     = currIndx+1
    
    do k = 1,kmax
       area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-1)*(numINperCell)+k
       currIndx                     = currIndx+1
    enddo
    
    area_node(N_cell,currIndx+1) = VF_BotmSurfcNodes(1)
    currIndx                     = currIndx+1
    
    do k = 1,kmax
       area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell*2)*(numINperCell)+k
       currIndx                     = currIndx+1
    enddo
    
    area_node(N_cell,currIndx+1) = VF_BotmSurfcNodes(2)
    currIndx                     = currIndx+1
    
    do k = 1,kmax
       area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2 + ((Hlf_Ncell*2)-1)*(numINperCell) + (NAEC_Apcl)-(k-1)
       currIndx                     = currIndx+1
    enddo
    
    area_node(N_cell,currIndx+1) = VF_BotmSurfcNodes(2)-2
    currIndx                     = currIndx+1
    
    do k = 1,kmax
       area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2 + ((Hlf_Ncell*2)-2)*(numINperCell) + (NAEC_Apcl)-(k-1)
       currIndx                     = currIndx+1
    enddo
    
    area_node(N_cell,currIndx+1) = VF_UpprSurfcNodes(2)
    
    write(*,*) area_node(N_cell,0:max_node_area),"AN VFcell 2"
    
    
    !area_spr(N_cell,0)                         = 1 + 1*sprAS + (2*1)*sprAS !(previous case)
    !area_spr(N_cell,1:sprAS)                   = sprApclNC1L(1:sprAS)
    !area_spr(N_cell,(sprAS+1))                 = sprVFtopLayr
    !area_spr(N_cell,(sprAS+2):(2*sprAS+1))     = sprApclCntrl(1:sprAS)
    !area_spr(N_cell,(2*sprAS+2):(3*sprAS+1))   = sprApclNC1R(1:sprAS)
    
    area_spr(N_cell,0)                       = 0 + 1*sprAS + (2*2)*sprAS
    
    area_spr(N_cell,(0*sprAS+1):(1*sprAS))   = sprApclNC2L(1:sprAS)
    area_spr(N_cell,(1*sprAS+1):(2*sprAS))   = sprApclNC1L(1:sprAS)
    area_spr(N_cell,(2*sprAS+1):(3*sprAS))   = sprApclCntrl(1:sprAS)
    area_spr(N_cell,(3*sprAS+1):(4*sprAS))   = sprApclNC1R(1:sprAS)
    area_spr(N_cell,(4*sprAS+1):(5*sprAS))   = sprApclNC2R(1:sprAS)
    
    write(*,*) area_spr(N_cell,0:max_spr_area),"AS VFcell 2"
    
  end subroutine readjst_area_to_other_trnsfrmsVF_region_methd2
  
  
  subroutine readjst_area_to_other_trnsfrmsVF_region_methd3
    implicit none
    integer :: k,kmax,i,imax,imaxHlf
    integer :: currIndx
    integer :: sprAS,sprBS,sprLS !(spr in Apcl/Basal/Ltrl) side
    integer :: num_sgmntA,num_sgmntB
    integer :: hlf_sgmntA,hlf_sgmntB
    integer :: cnt_sgmntA,cnt_sgmntB
    integer :: diffB1,diffB2
    
    sprAS = NAEC_Apcl+1 ; sprBS = NAEC_Bsal+1 ; sprLS = NAEC_Ltrl+1 
    
    area_node(1:N_cellSNVFR,0:max_node_areaS) = area_nodeS(1:(N_cell-1),0:max_node_areaS)
    area_spr(1:N_cellSNVFR,0:max_spr_areaS)   = area_sprS(1:N_cellSNVFR,0:max_spr_areaS)
    
    num_sgmntA = 2 + (addedNCPair*2)     ; num_sgmntB = 1 + (addedNCPair*2)
    hlf_sgmntA = num_sgmntA/2            ; hlf_sgmntB = num_sgmntB/2
    write(*,*) num_sgmntA,hlf_sgmntA,num_sgmntB,hlf_sgmntB,"num_sgmnt and half"
    
    imax = num_sgmntA + num_sgmntB ; write(*,*) imax,"imax"
    kmax = NAEC_Apcl
    
    area_node(N_cell,0) = 2*(1+addedNCPair) + (1+(2*addedNCPair))*(NAEC_Apcl);write(*,*)area_node(N_cell,0),"AN"
    
    currIndx=0 ; cnt_sgmntA=0 ; cnt_sgmntB=0
    
    do i = 1,imax
       
       if (mod(i,2)==1) then
          cnt_sgmntA = cnt_sgmntA+1
          
          if (cnt_sgmntA.le.hlf_sgmntA) then
             
             area_node(N_cell,currIndx+1) = VF_UpprSurfcNodes(1) + 2*(cnt_sgmntA-1)
             currIndx = currIndx+1
             
          elseif ((cnt_sgmntA.gt.hlf_sgmntA).and.(cnt_sgmntA.le.num_sgmntA)) then
             area_node(N_cell,currIndx+1) = VF_BotmSurfcNodes(2) - 2*(cnt_sgmntA-hlf_sgmntA-1)
             currIndx = currIndx+1
             
          endif
             
       elseif (mod(i,2)==0) then
          cnt_sgmntB = cnt_sgmntB+1
          
          if (cnt_sgmntB .le. hlf_sgmntB) then
             
             diffB1 = hlf_sgmntB - (cnt_sgmntB-1)
             
             do k = 1,kmax
                area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-diffB1)*(numINperCell)+k
                currIndx                     = currIndx+1
             enddo
             
          elseif (cnt_sgmntB == (hlf_sgmntB+1)) then
             
             do k = 1,kmax
                area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell*2)*(numINperCell)+k
                currIndx                     = currIndx+1
             enddo
             
          elseif ((cnt_sgmntB.gt.(hlf_sgmntB+1)) .and. (cnt_sgmntB.le.num_sgmntB)) then
             
             diffB2 = cnt_sgmntB - (hlf_sgmntB+1)
             
             do k = 1,kmax
                area_node(N_cell,currIndx+1) = 2*(Hlf_Ncell+1)*2+((Hlf_Ncell*2)-diffB2)*(numINperCell)+&
                     (NAEC_Apcl)-(k-1)
                currIndx                     = currIndx+1
             enddo
             
          endif
          
       endif
       
    enddo

    do i = 1,area_node(N_cell,0)
       write(*,*) i,area_node(N_cell,i),"ani"
    enddo

    write(*,*) " "
    write(*,*) area_node(N_cell,0:max_node_area),"anii"
    
    area_spr(N_cell,0) = 0 + 1*sprAS + (2*addedNCPair)*sprAS
    
    imax    = 1 + (2*addedNCPair)
    imaxHlf = imax/2 ; write(*,*) imaxHlf,"imaxHlf"
    
    do i = 1,imax
       
       if (i.le.imaxHlf) then
          area_spr(N_cell,((i-1)*sprAS+1):(i*sprAS)) = sprApclNCL((addedNCPair-(i-1)),1:sprAS)
          
       elseif (i==(imaxHlf+1)) then
          area_spr(N_cell,((i-1)*sprAS+1):(i*sprAS)) = sprApclCntrl(1:sprAS)
          
       elseif ((i.gt.(imaxHlf+1)) .and. (i.le.imax)) then
          area_spr(N_cell,((i-1)*sprAS+1):(i*sprAS)) = sprApclNCR((i-(imaxHlf+1)),1:sprAS)
          
       endif
       
    enddo
    
    write(*,*) area_spr(N_cell,0:max_spr_area),"AS i"
    
  end subroutine readjst_area_to_other_trnsfrmsVF_region_methd3
  
  
  subroutine readjst_all_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    
    call readjst_node_to_other_trnsfrms_crtclDesgnd_apclSurfc
    call readjst_spr_to_other_trnsfrms_crtclDesgnd_apclSurfc
    call readjst_area_to_other_trnsfrms_crtclDesgnd_apclSurfc
    call print_all_trnsfrms ; write(*,*) "MADE UPTO THIS crtclDesgnd_apclSurfc"
    call get_Nlist
    call print_all_trnsfrms
    
  end subroutine readjst_all_trnsfrms_crtclDesgnd_apclSurfc
  
  
  subroutine readjst_node_to_other_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    integer :: cntCurrNode,cntPrvNode
    integer :: cntCurrSpr ,cntPrvSpr
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: i,j,imax,jmax
    integer :: sprNmPrv,sprNmCurr
    
    cntCurrNode=1 ; cntPrvNode=1
    cntCurrSpr =1 ; cntPrvSpr =1
    
    !write(*,*) N_nodeS,N_node,"insd sb:readjst_node_to_other_trnsfrms_crtclDesgnd_apclSurfc"
    
    lim1 = 2*(Hlf_Ncell+1)*2
    lim2 = lim1 + (Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1)
    lim3 = lim2 + (NCP_CrtclApSrfc)*(INperCellT1)
    lim4 = lim3 + (Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1)
    lim5 = lim4 + (NCP_CrtclApSrfc)*(INperCellT1)
    lim6 = lim5 + (NAEC_Apcl+NAEC_Bsal)
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lim 1-6"
    
    do i = 1,N_nodeS
       
       if (i.le.lim1) then
          node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
          jmax                                    = node_sprS(i,0)
          node_spr(cntCurrNode,0)                 = jmax
          
          do j = 1,jmax
             sprNmPrv  = node_sprS(i,j)
             call get_currSprNm_crtclDesgnd_apclSurfc(sprNmPrv,sprNmCurr)
             node_spr(cntCurrNode,j) = sprNmCurr
          enddo
          
          cntCurrNode = cntCurrNode+1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
          jmax                                    = node_sprS(i,0)
          node_spr(cntCurrNode,0)                 = jmax
          node_spr(cntCurrNode,1:jmax)            = node_sprS(i,1:jmax)
          cntCurrNode                             = cntCurrNode+1
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          if ((i-lim2).le.NAEC_Apcl) then
             
             if ((i-lim2)==1) then
                
                do j = 1,NAEC_ApclCrtcl
                   node_spr(cntCurrNode,0)   = 2
                   
                   if (j==1) then
                      node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
                      node_spr(cntCurrNode,1:2)               = node_sprS(i,1:2)
                      cntCurrNode                             = cntCurrNode+1
                   elseif (j.gt.1) then
                      node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
                      node_spr(cntCurrNode,1:2)               = node_sprS(i,1:2) + (j-1)
                      cntCurrNode                             = cntCurrNode+1
                   endif
                   
                enddo
                
             elseif ((i-lim2).gt.1) then
                continue
             endif
             
          elseif ((i-lim2).gt.NAEC_Apcl) then
             node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
             jmax                                    = node_sprS(i,0)
             node_spr(cntCurrNode,0)                 = jmax
             node_spr(cntCurrNode,1:jmax)            = node_sprS(i,1:jmax) + (addApNodesPerApSpr)
             cntCurrNode                             = cntCurrNode+1
          endif
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
          jmax                                    = node_sprS(i,0)
          node_spr(cntCurrNode,0)                 = jmax
          node_spr(cntCurrNode,1:jmax)            = node_sprS(i,1:jmax) + (addApNodesPerApSpr)
          cntCurrNode                             = cntCurrNode+1
          
        
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          if ((i-lim4).le.NAEC_Apcl) then
             
             if ((i-lim4)==1) then
             
                do j = 1,NAEC_ApclCrtcl
                   node_spr(cntCurrNode,0) = 2
                   
                   if (j==1) then
                      node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
                      node_spr(cntCurrNode,1:2)               = node_sprS(i,1:2) + (addApNodesPerApSpr)
                      cntCurrNode                             = cntCurrNode+1
                   elseif (j.gt.1) then
                      node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
                      node_spr(cntCurrNode,1:2)               = node_sprS(i,1:2) + (addApNodesPerApSpr) + (j-1)
                      cntCurrNode                             = cntCurrNode+1
                   endif
                
                enddo
                
             elseif ((i-lim4).gt.1) then
                continue
             endif
             
          elseif ((i-lim4).gt.NAEC_Apcl) then
             node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
             jmax                                    = node_sprS(i,0)
             node_spr(cntCurrNode,0)                 = jmax
             node_spr(cntCurrNode,1:jmax)            = node_sprS(i,1:jmax) + 2*(addApNodesPerApSpr)
             cntCurrNode                             = cntCurrNode+1
          endif
          
       elseif  ((i.gt.lim5) .and. (i.le.lim6)) then
          
          if ((i-lim5).le.NAEC_Apcl) then
             
             if ((i-lim5)==1) then
                
                do j = 1,NAEC_ApclCrtcl
                   node_spr(cntCurrNode,0)   = 2
                   
                   if (j==1) then
                      node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
                      node_spr(cntCurrNode,1:2)               = node_sprS(i,1:2) + 2*(addApNodesPerApSpr)
                      cntCurrNode                             = cntCurrNode+1
                   elseif (j.gt.1) then
                      node_area(cntCurrNode,0:max_area_node) = node_areaS(i,0:max_area_node)
                      node_spr(cntCurrNode,1:2)              = node_sprS(i,1:2) + 2*(addApNodesPerApSpr)+ (j-1)
                      cntCurrNode                            = cntCurrNode+1
                   endif
                
                enddo
                
             elseif ((i-lim5).gt.1) then
                continue
             endif
             
          elseif ((i-lim5).gt.NAEC_Apcl) then
             node_area(cntCurrNode,0:max_area_nodeS) = node_areaS(i,0:max_area_nodeS)
             jmax                                    = node_sprS(i,0)
             node_spr(cntCurrNode,0)                 = jmax
             node_spr(cntCurrNode,1:jmax)            = node_sprS(i,1:jmax) + 3*(addApNodesPerApSpr)
             cntCurrNode                             = cntCurrNode+1
          endif
          
       endif
       
    enddo
    
  end subroutine readjst_node_to_other_trnsfrms_crtclDesgnd_apclSurfc
  
  
  subroutine get_currSprNm_crtclDesgnd_apclSurfc(sprNmPrv,sprNmCurr)
    implicit none
    integer, intent(in)  :: sprNmPrv
    integer, intent(out) :: sprNmCurr
    integer              :: firstMrgin,secndMrgin,thirdMrgin,forthMrgin
    
    firstMrgin = (Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1+3) + 1
    secndMrgin = (Hlf_Ncell+Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1+3) + 1
    thirdMrgin = (2*Hlf_Ncell)*(INperCellT1+3) + 1
    forthMrgin = (2*Hlf_Ncell)*(INperCellT1+3) + (NAEC_Apcl+1+NAEC_Bsal+1)
    
    !write(*,*) firstMrgin,secndMrgin,thirdMrgin,forthMrgin,"margins in sb:get_currSprNm_crtclDesgnd_apclSurfc"
    
    if (sprNmPrv.le.firstMrgin) then
       sprNmCurr = sprNmPrv
    elseif ((sprNmPrv.gt.firstMrgin) .and. (sprNmPrv.le.secndMrgin)) then
       sprNmCurr = sprNmPrv + addApNodesPerApSpr
    elseif ((sprNmPrv.gt.secndMrgin) .and. (sprNmPrv.le.thirdMrgin)) then
       sprNmCurr = sprNmPrv + 2*addApNodesPerApSpr
    elseif ((sprNmPrv.gt.thirdMrgin) .and. (sprNmPrv.le.forthMrgin)) then
       sprNmCurr = sprNmPrv + 3*addApNodesPerApSpr   
    endif
    
    !write(*,*) sprNmCurr,sprNmPrv,"sprnm curr prv"
    
  end subroutine get_currSprNm_crtclDesgnd_apclSurfc
  
  
  subroutine readjst_spr_to_other_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: cntCurrSpr
    integer :: i,j,jmax
    integer :: nodeNmFrmPrvSys
    integer :: cntAddNodeMltplr
    
    lim1 = (Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1+3)
    lim2 = (Hlf_Ncell)*(INperCellT1+3)
    lim3 = (2*Hlf_Ncell-NCP_CrtclApSrfc)*(INperCellT1+3)
    lim4 = (2*Hlf_Ncell)*(INperCellT1+3)
    lim5 = (2*Hlf_Ncell)*(INperCellT1+3) + (NAEC_Apcl+NAEC_Bsal+2)
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,"lims in sb: readjst_spr_to_other_trnsfrms_crtclDesgnd_apclSurfc"
    
    cntCurrSpr = 1
    
    do i = 1,N_sprS
       
       if (i.le.lim1) then
          spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(i,0:max_area_sprS)
          spr_node(cntCurrSpr,0:max_node_sprS) = spr_nodeS(i,0:max_node_sprS)
          cntCurrSpr                           = cntCurrSpr + 1
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          cntAddNodeMltplr = 1
          call get_sprOthr_readjst_crtclDesgnd_apclSurfc(i,lim1,cntAddNodeMltplr,cntCurrSpr)
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(i,0:max_area_sprS)
          spr_node(cntCurrSpr,0)               = spr_nodeS(i,0)
          
          do j = 1,spr_nodeS(i,0)
             nodeNmFrmPrvSys = spr_nodeS(i,j)
             if (nodeNmFrmPrvSys.le.numRegulrNode) then
                spr_node(cntCurrSpr,j) = spr_nodeS(i,j)
             elseif (nodeNmFrmPrvSys.gt.numRegulrNode) then
                spr_node(cntCurrSpr,j) = spr_nodeS(i,j) + (cntAddNodeMltplr*addApNodesPerApSpr)
             endif
          enddo
          cntCurrSpr = cntCurrSpr + 1
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          cntAddNodeMltplr = 2
          call get_sprOthr_readjst_crtclDesgnd_apclSurfc(i,lim3,cntAddNodeMltplr,cntCurrSpr)
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          cntAddNodeMltplr = 3
          call get_sprOthr_readjst_crtclDesgnd_apclSurfc(i,lim4,cntAddNodeMltplr,cntCurrSpr)
       endif
       
    enddo
    
  end subroutine readjst_spr_to_other_trnsfrms_crtclDesgnd_apclSurfc
  
  subroutine get_sprOthr_readjst_crtclDesgnd_apclSurfc(iV,limV,cntAddNodeMltplr,cntCurrSpr)
    implicit none
    integer, intent(in)  :: iV,limV,cntAddNodeMltplr
    integer, intent(out) :: cntCurrSpr
    
    integer              :: j,nodeNmFrmPrvSys
    
    if ((iV-limV).le.(NAEC_Apcl+1)) then
       
       if ((iV-limV)==1) then
          
          spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(iV,0:max_area_sprS) 
          spr_node(cntCurrSpr,0)               = 2
          spr_node(cntCurrSpr,1)               = spr_nodeS(iV,1)
          spr_node(cntCurrSpr,2)               = spr_nodeS(iV,2)+ (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
          cntCurrSpr                           = cntCurrSpr+1
          
       elseif ((iV-limV)==2) then
          
          do j = 1,(NAEC_ApclCrtcl-1)
             spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(iV,0:max_area_sprS)
             spr_node(cntCurrSpr,0)               = 2
             spr_node(cntCurrSpr,1:2)             = spr_nodeS(iV,1:2) + (j-1) + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
             cntCurrSpr                           = cntCurrSpr+1
          enddo
          
       elseif ((iV-limV)==(NAEC_Apcl+1)) then
          spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(iV,0:max_area_sprS)
          spr_node(cntCurrSpr,0)               = 2
          spr_node(cntCurrSpr,1)               = spr_nodeS(iV,1) + (cntAddNodeMltplr*addApNodesPerApSpr)
          spr_node(cntCurrSpr,2)               = spr_nodeS(iV,2)
          cntCurrSpr                           = cntCurrSpr+1
       endif
       
    elseif (((iV-limV).gt.(NAEC_Apcl+1)) .and. ((iV-limV).le.(INperCellT1+3))) then
       
       spr_area(cntCurrSpr,0:max_area_sprS) = spr_areaS(iV,0:max_area_sprS)
       spr_node(cntCurrSpr,0)               = 2
       
       do j = 1,spr_node(cntCurrSpr,0)
          nodeNmFrmPrvSys  = spr_nodeS(iV,j)   
          if (nodeNmFrmPrvSys.le.numRegulrNode) then
             spr_node(cntCurrSpr,j) = nodeNmFrmPrvSys
          elseif (nodeNmFrmPrvSys.gt.numRegulrNode) then
             spr_node(cntCurrSpr,j) = nodeNmFrmPrvSys + (cntAddNodeMltplr*addApNodesPerApSpr)
          endif
          
       enddo
       
       cntCurrSpr = cntCurrSpr+1
    endif
    
  end subroutine get_sprOthr_readjst_crtclDesgnd_apclSurfc
  
  
  subroutine readjst_area_to_other_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    
    call readjst_area_to_node_trnsfrms_crtclDesgnd_apclSurfc
    call readjst_area_to_spr_trnsfrms_crtclDesgnd_apclSurfc
    
  end subroutine readjst_area_to_other_trnsfrms_crtclDesgnd_apclSurfc
  
  subroutine readjst_area_to_node_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: i,j,jmax,jlm1,jlm2,jlm3,jlm4,k,kmax
    integer :: cntCurrArea
    integer :: nodeFrmPrv,nodeCurr
    integer :: cntAddNodeMltplr
    
    cntCurrArea = 1
    
    lim1 = Hlf_Ncell - NCP_CrtclApSrfc
    lim2 = Hlf_Ncell
    lim3 = Hlf_Ncell + Hlf_Ncell - NCP_CrtclApSrfc
    lim4 = Hlf_Ncell + Hlf_Ncell
    lim5 = Hlf_Ncell + Hlf_Ncell + 1
    lim6 = Hlf_Ncell + Hlf_Ncell + 1 + 1 ! VF_area
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,"lims in sb:readjst_area_to_node_trnsfrms_crtclDesgnd_apclSurfc"
    write(*,*) N_cellS,"_cellS" 
    do i = 1,N_cellS
       
       if (i.le.lim1) then
          area_node(cntCurrArea,0:max_node_areaS) = area_nodeS(i,0:max_node_areaS)
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then ! Anterior Side Cell
          
          area_node(cntCurrArea,0) = 4 + (2*NAEC_Ltrl) + (NAEC_Bsal) + (NAEC_ApclCrtcl)
          cntAddNodeMltplr         = (i-lim1)
          
          jmax = area_node(cntCurrArea,0)
          jlm1 = 1+(NAEC_Ltrl)+1
          jlm2 = jlm1+(NAEC_Bsal)+1+(NAEC_Ltrl)+1
          jlm3 = jlm2+(NAEC_Apcl+1)
          !write(*,*) jlm1,jlm2,jlm3,jmax,"jlm"
          
          do j = 1,jlm3
             
             if (j.le.jlm1) then
                area_node(cntCurrArea,j) = area_nodeS(i,j)
                
             elseif ((j.gt.jlm1) .and. (j.le.jlm2)) then
                
                nodeFrmPrv = area_nodeS(i,j)
                if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
                if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr*addApNodesPerApSpr)
                area_node(cntCurrArea,j) = nodeCurr
                !write(*,*) area_node(cntCurrArea,j),area_nodeS(i,j),i,j,"cc"
                
             elseif ((j.gt.jlm2) .and. (i.le.jlm3)) then
                
                if ((j-jlm2)==1) then
                   nodeFrmPrv = area_nodeS(i,j)
                   kmax       = NAEC_ApclCrtcl
                   do k = 1,kmax
                      area_node(cntCurrArea,j+k-1) = nodeFrmPrv+ (cntAddNodeMltplr*addApNodesPerApSpr)-(k-1)
                   enddo
                elseif ((j-jlm2).gt.1) then
                   continue
                endif
                
             endif
             
          enddo
          
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then

          area_node(cntCurrArea,0) = area_nodeS(i,0)
          jmax                     = area_node(cntCurrArea,0)
          
          do j = 1,jmax
             nodeFrmPrv = area_nodeS(i,j) !; write(*,*) nodeFrmPrv,j,"npv"
             if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
             if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr*addApNodesPerApSpr)
              area_node(i,j) = nodeCurr
          enddo
          
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          area_node(cntCurrArea,0) = 4 + (2*NAEC_Ltrl) + (NAEC_Bsal) + (NAEC_ApclCrtcl)
          cntAddNodeMltplr         = NCP_CrtclApSrfc + (i-lim3)
          
          jmax = area_node(cntCurrArea,0)
          jlm1 = 1+(NAEC_Apcl)+1
          jlm2 = jlm1+(NAEC_Ltrl)+1+(NAEC_Bsal)+1
          jlm3 = jlm2+(NAEC_Ltrl)
          !write(*,*) jlm1,jlm2,jlm3,jmax,"jlm"
          
          do j = 1,jlm3
             
             if (j.le.jlm1) then
                
                if (j==1) then
                   area_node(cntCurrArea,j) = area_nodeS(i,j)
                   
                elseif (j==2) then
                   nodeFrmPrv = area_nodeS(i,j)
                   kmax       = NAEC_ApclCrtcl
                   do k = 1,kmax
                      area_node(cntCurrArea,j+k-1) = nodeFrmPrv + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)+(k-1)
                      !write(*,*) (j+k-1),area_node(cntCurrArea,j+k-1),nodeFrmPrv,addApNodesPerApSpr,(cntAddNodeMltplr-1),j,k,"aa"
                   enddo
                   
                elseif (j==jlm1) then
                   area_node(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_nodeS(i,j)
                   
                else
                   continue
                endif
                
             elseif ((j.gt.jlm1) .and. (j.le.jlm2)) then
               
                nodeFrmPrv = area_nodeS(i,j)
                if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
                if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr*addApNodesPerApSpr)
                area_node(cntCurrArea,j+NAEC_ApclCrtcl-2) = nodeCurr
                
                !write(*,*) nodeFrmPrv,j,(j+NAEC_ApclCrtcl-2),nodeCurr,"jlm12"
                
             elseif ((j.gt.jlm2) .and. (j.le.jlm3)) then 
                area_node(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_nodeS(cntCurrArea,j) + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
             endif
             
          enddo
          
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          area_node(cntCurrArea,0) = 4 + (2*NAEC_Ltrl) + (NAEC_Bsal) + (NAEC_ApclCrtcl)
          cntAddNodeMltplr         = cntAddNodeMltplr  + 1
          
          jmax = area_node(cntCurrArea,0)
          jlm1 = 1+(NAEC_Apcl)+1
          jlm2 = jlm1+(NAEC_Ltrl)+1
          jlm3 = jlm2+(NAEC_Bsal)+1
          jlm4 = jlm3+(NAEC_Ltrl)
          
          !write(*,*) jlm1,jlm2,jlm3,jlm4,jmax,"jlm"
          
          do j = 1,jlm4
             
             if (j.le.jlm1) then
                
                nodeFrmPrv = area_nodeS(i,j)
                if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
                if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr-1-NCP_CrtclApSrfc)*(addApNodesPerApSpr)
                area_node(cntCurrArea,j) = nodeCurr
                
             elseif ((j.gt.jlm1) .and. (j.le.jlm2)) then
                
                nodeFrmPrv = area_nodeS(i,j)
                if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
                if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr)*(addApNodesPerApSpr)
                area_node(cntCurrArea,j) = nodeCurr
                
             elseif ((j.gt.jlm2) .and. (j.le.jlm3)) then
                
                nodeFrmPrv = area_nodeS(i,j)
                if (nodeFrmPrv.le.numRegulrNode) nodeCurr = nodeFrmPrv
                if (nodeFrmPrv.gt.numRegulrNode) nodeCurr = nodeFrmPrv + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
                area_node(cntCurrArea,j) = nodeCurr
                
             elseif ((j.gt.jlm3) .and. (j.le.jlm4)) then
                
                if ((j-jlm3)==1) then
                   nodeFrmPrv = area_nodeS(i,j)
                   kmax       = NAEC_ApclCrtcl
                   do k = 1,kmax
                      area_node(cntCurrArea,j+k-1) = nodeFrmPrv + (cntAddNodeMltplr*addApNodesPerApSpr)-(k-1)
                   enddo
                   
                else
                   continue
                endif
                
             endif
             
          enddo
          
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          call readjst_VF_area_node_with_moreAddedNode
       endif
       
    enddo  
    
  end subroutine readjst_area_to_node_trnsfrms_crtclDesgnd_apclSurfc
  
  subroutine readjst_VF_area_node_with_moreAddedNode
    implicit none
    integer :: numNodePrvSys,lftNodeVF,rghtNodeVF
    integer :: i,j,unDistrbdLftNode,unDistrbdRghtNode
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: cntMltpL,cntMltpR,cntCurrNode
    integer :: incrsOfNI
    
    numNodePrvSys = area_nodeS(N_cell,0)
    lftNodeVF     = (numNodePrvSys-(NAEC_Apcl+2))/2
    rghtNodeVF    = lftNodeVF
    
    unDistrbdLftNode  = max((lftNodeVF-NCP_CrtclApSrfc*(NAEC_Apcl+1)),0)
    unDistrbdRghtNode = unDistrbdLftNode
    
    write(*,*) numNodePrvSys,lftNodeVF,rghtNodeVF,"c1"
    write(*,*) unDistrbdLftNode,unDistrbdRghtNode,"c2"
    
    lim1 = unDistrbdLftNode
    lim2 = lftNodeVF
    lim3 = lftNodeVF+NAEC_Apcl+2
    lim4 = lftNodeVF+NAEC_Apcl+2+rghtNodeVF-unDIstrbdRghtNode
    lim5 = lftNodeVF+NAEC_Apcl+2+rghtNodeVF

    write(*,*) lim1,lim2,lim3,lim4,lim5,"LIM"
    area_node(N_cell,0) = area_nodeS(N_cell,0) + ((NCP_CrtclApSrfc*2)+1)*(NAEC_ApclCrtcl-NAEC_Apcl) 
    write(*,*) area_node(N_cell,0),"area_node VF max"
    
    cntMltpL    = 1 ; cntMltpR = 1
    cntCurrNode = 1
    write(*,*) area_nodeS(N_cell,1:22)
    do i = 1,numNodePrvSys
       
       if (i.le.lim1) then
          area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
          !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a1"
          cntCurrNode = cntCurrNode + 1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if (mod((i-lim1),(NAEC_Apcl+1))==1) then
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
             write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a2"
             cntCurrNode = cntCurrNode+1
             
          elseif (mod((i-lim1),(NAEC_Apcl+1))==2) then
             
             do j = 1,NAEC_ApclCrtcl
                area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i) + (j-1) + (cntMltpL-1)*(addApNodesPerApSpr)
                !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a3"
                cntCurrNode = cntCurrNode + 1
             enddo
             cntMltpL = cntMltpL + 1
          else
             continue
          endif
          
       elseif ((i.gt.lim2).and.(i.le.lim3)) then
          
          if (mod((i-lim2),(NAEC_Apcl+2))==1) then
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
             !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a4"
             cntCurrNode = cntCurrNode+1
             
          elseif (mod((i-lim2),(NAEC_Apcl+2))==2) then
             
             do j = 1,NAEC_ApclCrtcl
                area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i) + (j-1) + (2*NCP_CrtclApSrfc)*(addApNodesPerApSpr)
                !write(*,*) area_node(N_cell,cntCurrNode),area_nodeS(N_cell,i),i,cntCurrNode,"a5"
                cntCurrNode = cntCurrNode + 1
             enddo
             
          elseif(mod((i-lim2),(NAEC_Apcl+2))==0) then
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
             !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a6"
             cntCurrNode = cntCurrNode+1
             
          else
             continue
          endif
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          if (mod((i-lim3),(NAEC_Apcl+1))==1) then
             
             do j = 1,NAEC_ApclCrtcl
                area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i) - (j-1) + (2*NCP_CrtclApSrfc-(cntMltpR-1))*(addApNodesPerApSpr)
                !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a7"
                cntCurrNode = cntCurrNode + 1
             enddo
             cntMltpR = cntMltpR + 1
             
          elseif (mod((i-lim3),(NAEC_Apcl+1))==0) then
             
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
             !write(*,*) area_node(N_cell,cntCurrNode),i,cntCurrNode,"a8"
             cntCurrNode = cntCurrNode+1
             
          else
             continue
          endif
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          if (mod((i-lim4),(NAEC_Apcl+1))==0) then
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i)
             !write(*,*) area_node(N_cell,cntCurrNode),area_nodeS(N_cell,i),i,cntCurrNode,"a9"
             cntCurrNode = cntCurrNode + 1
             
          elseif(mod((i-lim4),(NAEC_Apcl+1)).ne.0) then
             area_node(N_cell,cntCurrNode) = area_nodeS(N_cell,i) + (NCP_CrtclApSrfc*addApNodesPerApSpr)
             !write(*,*) area_node(N_cell,cntCurrNode),area_nodeS(N_cell,i),i,cntCurrNode,"a10"
             cntCurrNode = cntCurrNode + 1
          endif
          
       endif
       
    enddo
    
  end subroutine readjst_VF_area_node_with_moreAddedNode
  
  
  subroutine readjst_area_to_spr_trnsfrms_crtclDesgnd_apclSurfc
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: i,j,jmax,jlm1,jlm2,jlm3,jlm4,k,kmax
    integer :: cntCurrArea
    integer :: nodeFrmPrv,nodeCurr
    integer :: cntAddNodeMltplr
    integer :: maxCnnSprInArea
    
    cntCurrArea      = 1
    cntAddNodeMltplr = 0
    
    lim1 = Hlf_Ncell - NCP_CrtclApSrfc
    lim2 = Hlf_Ncell
    lim3 = Hlf_Ncell + Hlf_Ncell - NCP_CrtclApSrfc
    lim4 = Hlf_Ncell + Hlf_Ncell
    lim5 = Hlf_Ncell + Hlf_Ncell + 1
    lim6 = Hlf_Ncell + Hlf_Ncell + 1 + 1
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims in sb:readjst_area_to_spr_trnsfrms_crtclDesgnd_apclSurfc"
    
    do i = 1,N_cellS

        if (i.le.lim1) then
           area_spr(cntCurrArea,0:max_spr_areaS) = area_sprS(i,0:max_spr_areaS)
           cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then ! Anterior Side Distrbd Cell
          
          area_spr(cntCurrArea,0) = (NAEC_ApclCrtcl+1) + (NAEC_Bsal+1) + 2*(NAEC_Ltrl+1)
          jmax                    = area_sprS(i,0)
          
          jlm1                    = NAEC_Ltrl+1
          jlm2                    = jlm1 + NAEC_Apcl+1
          jlm3                    = jlm2 + NAEC_Bsal+1+NAEC_Ltrl+1
          
          cntAddNodeMltplr        = cntAddNodeMltplr+1
          write(*,*) jlm1,jlm2,jlm3,jmax,"jlm"
          
          do j = 1,jmax
             
             if (j.le.jlm1) then
                area_spr(cntCurrArea,j) = area_sprS(i,j)
                write(*,*) area_spr(cntCurrArea,j),area_sprS(i,j),i,j,"sA1"
                
             elseif ((j.gt.jlm1).and.(j.le.jlm2)) then   
                
                if ((j-jlm1)==1) then
                   kmax = NAEC_ApclCrtcl+1
                   do k = 1,kmax
                      area_spr(cntCurrArea,j+k-1) = area_sprS(i,j) + (k-1)
                      write(*,*) area_spr(cntCurrArea,j+k-1),area_sprS(i,j),i,j+k-1,k,"sA2"
                   enddo
                elseif ((j-jlm1).gt.1) then
                   continue
                endif
                
             elseif ((j.gt.jlm2).and.(j.le.jlm3)) then  
                area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_sprS(i,j) + (cntAddNodeMltplr*addApNodesPerApSpr)
                write(*,*) area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2),area_sprS(i,j),i,j+NAEC_ApclCrtcl-2,"sA3"
             endif
             
          enddo
          
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          area_spr(cntCurrArea,0) = area_sprS(i,0) 
          maxCnnSprInArea         = area_spr(cntCurrArea,0)
          !write(*,*) maxCnnSprInArea,"mx"
          
          area_spr(cntCurrArea,1:maxCnnSprInArea) = area_sprS(i,1:maxCnnSprInArea)+(cntAddNodeMltplr*addApNodesPerApSpr)
          !write(*,*) area_spr(cntCurrArea,1:maxCnnSprInArea),"mx2"
          !write(*,*) area_spr(cntCurrArea,(maxCnnSprInArea+1):(maxCnnSprInArea+3)),"mx3"
          cntCurrArea = cntCurrArea+1
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then ! Posterior Side Distrbd Cell
          
          area_spr(cntCurrArea,0) = (NAEC_Apcl+1) + (NAEC_Bsal+1) + 2*(NAEC_Ltrl+1)
          jmax                    = area_sprS(i,0)
          
          jlm1                    = NAEC_Ltrl+1
          jlm2                    = jlm1 + NAEC_Apcl+1
          jlm3                    = jlm2 + NAEC_Bsal+1+NAEC_Ltrl+1
          
          cntAddNodeMltplr        = cntAddNodeMltplr+1
          
          do j = 1,jmax
             
             if (j.le.jlm1) then
                area_spr(cntCurrArea,j) = area_sprS(i,j) + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
                write(*,*) area_spr(cntCurrArea,j),area_sprS(i,j),i,cntCurrArea,j,"sP1"
                
             elseif ((j.gt.jlm1).and.(j.le.jlm2)) then  
                
                if ((j-jlm1)==1) then
                   kmax = NAEC_ApclCrtcl+1
                   do k = 1,kmax
                      area_spr(cntCurrArea,j+k-1) = area_sprS(i,j) + (k-1) + (cntAddNodeMltplr-1)*(addApNodesPerApSpr)
                      write(*,*) area_spr(cntCurrArea,j+k-1),area_sprS(i,j),i,j+k-1,k,"sP2"
                   enddo
                elseif ((j-jlm1).gt.1) then
                   continue
                endif
                
             elseif ((j.gt.jlm2).and.(j.le.jlm3)) then 
                area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_sprS(i,j) + (cntAddNodeMltplr*addApNodesPerApSpr)
                write(*,*) area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2),area_sprS(i,j),i,j+NAEC_ApclCrtcl-2,"sP3"
             endif
             
          enddo
          
          cntCurrArea = cntCurrArea+1

       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          
          area_spr(cntCurrArea,0) = (NAEC_Ltrl+1) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
          jmax                    = area_spr(cntCurrArea,0)
          
          jlm1                    = NAEC_Ltrl+1
          jlm2                    = jlm1+NAEC_Apcl+1
          jlm3                    = jlm2+NAEC_Bsal+1
          jlm4                    = jlm3+NAEC_Ltrl+1
          
          cntAddNodeMltplr        = cntAddNodeMltplr+1
          
          do j = 1,jmax
             
             if (j.le.jlm1) then
                area_spr(cntCurrArea,j) = area_sprS(i,j) + (NCP_CrtclApSrfc)*(cntAddNodeMltplr)
                
             elseif ((j.gt.jlm1).and.(j.le.jlm2))  then
                
                if ((j-jlm1)==1) then
                   kmax = NAEC_ApclCrtcl+1
                   do k = 1,kmax
                      area_spr(cntCurrArea,j+k-1) = area_sprS(i,j) + (k-1) + (2*NCP_CrtclApSrfc)*(addApNodesPerApSpr)
                   enddo
                elseif ((j-jlm1).gt.1) then
                   continue
                endif
                
             elseif ((j.gt.jlm2).and.(j.le.jlm3)) then  
                area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_sprS(i,j) + (2*NCP_CrtclApSrfc+1)*(addApNodesPerApSpr)
             elseif ((j.gt.jlm3).and.(j.le.jlm4)) then 
                area_spr(cntCurrArea,j+NAEC_ApclCrtcl-2) = area_sprS(i,j) + (2*NCP_CrtclApSrfc)*(addApNodesPerApSpr)
             endif
             
          enddo
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          call readjst_VF_area_spr_with_moreAddedNode
       endif
      
    enddo
    
  end subroutine readjst_area_to_spr_trnsfrms_crtclDesgnd_apclSurfc
  
  
  subroutine readjst_VF_area_spr_with_moreAddedNode
    implicit none
    integer :: numSprPrvSys,lftSprVF,rghtSprVF
    integer :: i,j,unDistrbdLftSpr,unDistrbdRghtSpr
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: cntMltpL,cntMltpR,cntCurrSpr
    integer :: incrsOfNI
    
    numSprPrvSys = area_nodeS(N_cell,0)
    lftSprVF     = (numSprPrvSys-(NAEC_Apcl+1))/2
    rghtSprVF    = lftSprVF
    
    unDistrbdLftSpr  = max((lftSprVF-NCP_CrtclApSrfc*(NAEC_Apcl+1)),0)
    unDistrbdRghtSpr = unDistrbdLftSpr
    
    write(*,*) numSprPrvSys,lftSprVF,rghtSprVF,"c1"
    write(*,*) unDistrbdLftSpr,unDistrbdRghtSpr,"c2"
    
    lim1 = unDistrbdLftSpr
    lim2 = lftSprVF
    lim3 = lftSprVF+NAEC_Apcl+1
    lim4 = lftSprVF+NAEC_Apcl+1+rghtSprVF-unDIstrbdRghtSpr
    lim5 = lftSprVF+NAEC_Apcl+1+rghtSprVF
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,"LIM SPR"
    area_spr(N_cell,0) = area_sprS(N_cell,0) + ((NCP_CrtclApSrfc*2)+1)*(NAEC_ApclCrtcl-NAEC_Apcl) 
    
    cntMltpL    = 1 ; cntMltpR = 1
    cntCurrSpr = 1
    write(*,*) area_nodeS(N_cell,1:21),"AS 1-21"
    
    do i = 1,numSprPrvSys
       
       if (i.le.lim1) then
          
          area_spr(N_cell,cntCurrSpr) = area_sprS(N_cell,i)
          !write(*,*) area_spr(N_cell,cntCurrSpr),i,cntCurrSpr,"a1"
          cntCurrSpr = cntCurrSpr + 1
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if ((i-lim1)==1) then
             
             do j = 1,(NAEC_ApclCrtcl+1)
                area_spr(N_cell,cntCurrSpr) = area_sprS(N_cell,i) + (j-1) + (cntMltpL-1)*(addApNodesPerApSpr)
                !write(*,*) area_spr(N_cell,cntCurrSpr),i,cntCurrSpr,"a3"
                cntCurrSpr = cntCurrSpr + 1
             enddo
             cntMltpL = cntMltpL + 1
          else
             continue
          endif
          
       elseif ((i.gt.lim2).and.(i.le.lim3)) then
          
        
             
          if ((i-lim2)==1) then
             
             do j = 1,NAEC_ApclCrtcl+1
                area_spr(N_cell,cntCurrSpr) = area_sprS(N_cell,i) + (j-1) + (2*NCP_CrtclApSrfc)*(addApNodesPerApSpr)
                !write(*,*) area_spr(N_cell,cntCurrSpr),area_sprS(N_cell,i),i,cntCurrSpr,"a5"
                cntCurrSpr = cntCurrSpr + 1
             enddo
             
          else
             continue
          endif
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          if ((i-lim3)==1) then
             
             do j = 1,NAEC_ApclCrtcl+1
                area_spr(N_cell,cntCurrSpr) = area_sprS(N_cell,i) + (j-1) + (cntMltpR)*(addApNodesPerApSpr)
                !write(*,*) area_spr(N_cell,cntCurrSpr),i,cntCurrSpr,"a7"
                cntCurrSpr = cntCurrSpr + 1
             enddo
             cntMltpR = cntMltpR + 1
             
          else
             continue
          endif
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
             
             area_spr(N_cell,cntCurrSpr) = area_sprS(N_cell,i) + (NCP_CrtclApSrfc*addApNodesPerApSpr)
             !write(*,*) area_spr(N_cell,cntCurrSpr),area_sprS(N_cell,i),i,cntCurrSpr,"a10"
             cntCurrSpr = cntCurrSpr + 1
          
       endif
       
    enddo
    
  end subroutine readjst_VF_area_spr_with_moreAddedNode
  
  
  subroutine readjst_all_trnsfrms_frm_reading_cellsmeet_file
    implicit none
    
    call readjst_node_to_other_trnsfrms_frm_reading_cellsmeet_file
    call readjst_spr_to_other_trnsfrms_frm_reading_cellsmeet_file
    call readjst_area_to_other_trnsfrms_frm_reading_cellsmeet_file
    
    call get_Nlist
    write(*,*) Nlist(1,1,1:2),"1"
    call print_all_trnsfrms
    write(*,*) Nlist(1,1,1:2),"2"
    
    write(*,*) "IT works upTo readjst_file"
    !stop
    
  end subroutine readjst_all_trnsfrms_frm_reading_cellsmeet_file
  
  subroutine readjst_node_to_other_trnsfrms_frm_reading_cellsmeet_file
    implicit none
    integer :: i,j
    integer :: nodeV
    
    if (modelID==1) then

       if (CellsMeet==0) then
          open(unit=271,file='node_spr_cellsMeetsTN_eq_0.dat')
          open(unit=272,file='node_area_cellsMeetsTN_eq_0.dat')
          open(unit=371,file='node_other_cellsMeetsTN_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=271,file='node_spr_cellsMeetsTN_eq_1.dat')
          open(unit=272,file='node_area_cellsMeetsTN_eq_1.dat')
          open(unit=371,file='node_other_cellsMeetsTN_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=271,file='node_spr_cellsMeetsTN_eq_2.dat')
          open(unit=272,file='node_area_cellsMeetsTN_eq_2.dat')
          open(unit=371,file='node_other_cellsMeetsTN_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=271,file='node_spr_cellsMeetsTN_eq_3.dat')
          open(unit=272,file='node_area_cellsMeetsTN_eq_3.dat')
          open(unit=371,file='node_other_cellsMeetsTN_eq_3.dat')
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3 (TN)"
          stop
       endif
       
    elseif (modelID==2) then
       
       if (CellsMeet==0) then
          open(unit=271,file='node_spr_cellsMeetsNI_eq_0.dat')
          open(unit=272,file='node_area_cellsMeetsNI_eq_0.dat')
          open(unit=371,file='node_other_cellsMeetsNI_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=271,file='node_spr_cellsMeetsNI_eq_1.dat')
          open(unit=272,file='node_area_cellsMeetsNI_eq_1.dat')
          open(unit=371,file='node_other_cellsMeetsNI_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=271,file='node_spr_cellsMeetsNI_eq_2.dat')
          open(unit=272,file='node_area_cellsMeetsNI_eq_2.dat')
          open(unit=371,file='node_other_cellsMeetsNI_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=271,file='node_spr_cellsMeetsNI_eq_3.dat')
          open(unit=272,file='node_area_cellsMeetsNI_eq_3.dat')
          open(unit=371,file='node_other_cellsMeetsNI_eq_3.dat')
          
       elseif (CellsMeet==4) then
          open(unit=271,file='node_spr_cellsMeetsNI_eq_4.dat')
          open(unit=272,file='node_area_cellsMeetsNI_eq_4.dat')
          open(unit=371,file='node_other_cellsMeetsNI_eq_4.dat')
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3-4 (NI) "
          stop
       endif
       
    endif
    
    
    do i = 1,2
       
       do j = 1,N_node
          
          if (i==1) then
             read(271,*)  nodeV,node_spr(nodeV,0:max_spr_node) 
             write(371,*) nodeV,node_spr(nodeV,0:max_spr_node)
             
          elseif (i==2) then
             read(272,*)  nodeV,node_area(nodeV,0:max_area_node)
             write(371,*) nodeV,node_area(nodeV,0:max_area_node)
          endif
          
       enddo
       
       write(371,*) " "
       
    enddo
    
    close(271)
    close(272)
    close(371)
    
  end subroutine readjst_node_to_other_trnsfrms_frm_reading_cellsmeet_file
  
  
  subroutine readjst_spr_to_other_trnsfrms_frm_reading_cellsmeet_file 
    implicit none
    integer :: i,j
    integer :: sprV
    
    if (modelID==1) then
       
       if (CellsMeet==0) then
          open(unit=273,file='spr_node_cellsMeetsTN_eq_0.dat')
          open(unit=274,file='spr_area_cellsMeetsTN_eq_0.dat')
          open(unit=373,file='spr_other_cellsMeetsTN_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=273,file='spr_node_cellsMeetsTN_eq_1.dat')
          open(unit=274,file='spr_area_cellsMeetsTN_eq_1.dat')
          open(unit=373,file='spr_other_cellsMeetsTN_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=273,file='spr_node_cellsMeetsTN_eq_2.dat')
          open(unit=274,file='spr_area_cellsMeetsTN_eq_2.dat')
          open(unit=373,file='spr_other_cellsMeetsTN_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=273,file='spr_node_cellsMeetsTN_eq_3.dat')
          open(unit=274,file='spr_area_cellsMeetsTN_eq_3.dat')
          open(unit=373,file='spr_other_cellsMeetsTN_eq_3.dat')   
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3 (spr-adjsmnt-TN)"
          stop
       endif
       
    elseif (modelID==2) then
       
       if (CellsMeet==0) then
          open(unit=273,file='spr_node_cellsMeetsNI_eq_0.dat')
          open(unit=274,file='spr_area_cellsMeetsNI_eq_0.dat')
          open(unit=373,file='spr_other_cellsMeetsNI_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=273,file='spr_node_cellsMeetsNI_eq_1.dat')
          open(unit=274,file='spr_area_cellsMeetsNI_eq_1.dat')
          open(unit=373,file='spr_other_cellsMeetsNI_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=273,file='spr_node_cellsMeetsNI_eq_2.dat')
          open(unit=274,file='spr_area_cellsMeetsNI_eq_2.dat')
          open(unit=373,file='spr_other_cellsMeetsNI_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=273,file='spr_node_cellsMeetsNI_eq_3.dat')
          open(unit=274,file='spr_area_cellsMeetsNI_eq_3.dat')
          open(unit=373,file='spr_other_cellsMeetsNI_eq_3.dat')   
          
       elseif (CellsMeet==4) then
          open(unit=273,file='spr_node_cellsMeetsNI_eq_4.dat')
          open(unit=274,file='spr_area_cellsMeetsNI_eq_4.dat')
          open(unit=373,file='spr_other_cellsMeetsNI_eq_4.dat')   
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3-4 (spr-adjsmnt-NI)"
          stop
       endif
       
    endif
    
    write(*,*) N_spr,"nspr val in read"
    
    
    do i = 1,2
       
       do j = 1,N_spr
          
          if (i==1) then
             read(273,*)  sprV,spr_node(sprV,0:max_node_spr)
             write(373,*) sprV,spr_node(sprV,0:max_node_spr)
          elseif (i==2) then
             read(274,*)  sprV,spr_area(sprV,0:max_area_spr)
             write(373,*) sprV,spr_area(sprV,0:max_area_spr)
          endif
          
       enddo
       write(373,*) " "
    enddo
    
    close(273)
    close(274)
    close(373)
    
  end subroutine readjst_spr_to_other_trnsfrms_frm_reading_cellsmeet_file
  
  subroutine readjst_area_to_other_trnsfrms_frm_reading_cellsmeet_file
    implicit none
    integer :: i,j
    integer :: areaV
    
    
    if (modelID==1) then
       
       if (CellsMeet==0) then
          open(unit=275,file='area_node_cellsMeetsTN_eq_0.dat')
          open(unit=276,file='area_spr_cellsMeetsTN_eq_0.dat')
          open(unit=375,file='area_other_cellsMeetsTN_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=275,file='area_node_cellsMeetsTN_eq_1.dat')
          open(unit=276,file='area_spr_cellsMeetsTN_eq_1.dat')
          open(unit=375,file='area_other_cellsMeetsTN_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=275,file='area_node_cellsMeetsTN_eq_2.dat')
          open(unit=276,file='area_spr_cellsMeetsTN_eq_2.dat')
          open(unit=375,file='area_other_cellsMeetsTN_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=275,file='area_node_cellsMeetsTN_eq_3.dat')
          open(unit=276,file='area_spr_cellsMeetsTN_eq_3.dat')
          open(unit=375,file='area_other_cellsMeetsTN_eq_3.dat')
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3 (area-adjstmnt TN)"
          stop
       endif
       
    elseif (modelID==2) then
       
       if (CellsMeet==0) then
          open(unit=275,file='area_node_cellsMeetsNI_eq_0.dat')
          open(unit=276,file='area_spr_cellsMeetsNI_eq_0.dat')
          open(unit=375,file='area_other_cellsMeetsNI_eq_0.dat')
          
       elseif (CellsMeet==1) then
          open(unit=275,file='area_node_cellsMeetsNI_eq_1.dat')
          open(unit=276,file='area_spr_cellsMeetsNI_eq_1.dat')
          open(unit=375,file='area_other_cellsMeetsNI_eq_1.dat')
          
       elseif (CellsMeet==2) then
          open(unit=275,file='area_node_cellsMeetsNI_eq_2.dat')
          open(unit=276,file='area_spr_cellsMeetsNI_eq_2.dat')
          open(unit=375,file='area_other_cellsMeetsNI_eq_2.dat')
          
       elseif (CellsMeet==3) then
          open(unit=275,file='area_node_cellsMeetsNI_eq_3.dat')
          open(unit=276,file='area_spr_cellsMeetsNI_eq_3.dat')
          open(unit=375,file='area_other_cellsMeetsNI_eq_3.dat')
          
       elseif (CellsMeet==4) then
          open(unit=275,file='area_node_cellsMeetsNI_eq_4.dat')
          open(unit=276,file='area_spr_cellsMeetsNI_eq_4.dat')
          open(unit=375,file='area_other_cellsMeetsNI_eq_4.dat')
          
          
       else
          write(*,*) "CellsMeet not equal to anyone 0-1-2-3-4 (area-adjstmnt NI)"
          stop
       endif
       
    endif
    
    do i = 1,2
       
       do j = 1,N_cell
          
          if (i==1) then
             read(275,*)  areaV,area_node(areaV,0:max_node_area)
             write(375,*) areaV,area_node(areaV,0:max_node_area)
          elseif (i==2) then
             read(276,*)  areaV,area_spr(areaV,0:max_spr_area)
             write(375,*) areaV,area_spr(areaV,0:max_spr_area)
          endif
          
       enddo
       write(375,*) " "
    enddo
    
    close(275)
    close(276)
    close(375)
    
  end subroutine readjst_area_to_other_trnsfrms_frm_reading_cellsmeet_file
  
  subroutine redefining_system_parameters_stage1_type1
    implicit none
    
    open(unit=167,file='redefine_system_S1T1.dat')
    
    !lft_ladder_params
    ncl  = ncl - 1
    nsl  = nsecl * ncl
    nvsl = ncl + 1
    
    lft_node = 2*nvsl
    
    write(167,*) ncl,nsl,nvsl,lft_node,"ncl,nsl,nvsl,lft_Node"
    
    !rght_ladder_params
    ncr  = ncr - 1
    nsr  = nsecr * ncr
    nvsr = ncr + 1
    
    rght_node = 2*nvsr
    
    write(167,*) ncr,nsr,nvsr,rght_node,"ncl,nsl,nvsl,lft_Node"
    
    
    
    write(167,*) N_node,N_spr,N_cell,Hlf_Ncell,"N_node,N_spr,N_cell,Hlf_Ncell"
    write(167,*) N_lftCells,N_rghtCells,"N lft and rght Cells"

    close(167)
    
  end subroutine redefining_system_parameters_stage1_type1
  
  
  subroutine add_cellFor_VitellineFluid_region_and_redefine_system
    ! VF=VitellineFluid
    implicit none
    logical  :: lgcl_dealloc=.False.
    integer  :: i,j
    
    !addedNCPair = 3
    
    call getDecsn_of_dealloctn_reallocatn_of_movingVars(lgcl_dealloc)
    write(*,*) "lgcl_dealloc =",lgcl_dealloc
    
    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    call deallocate_and_reallocate_transfrmStr_variables
    call deallocate_all_gradient_variables_StrVars
    call get_all_gradient_variables_StrVars
    
    call store_SysVars_and_Arrays_bfrAddingCell_for_VF_region
    call get_Sysvars_aftAddingCell_for_VF_region_WO_PrtlApclSgmntn
    call deallocate_and_reallocate_arrays_forAdding_VF_region
    
    call get_NodeVars_VF_region
    call get_SprVars_for_VF_region
    call get_CellVars_VF_region
    call get_CgVars_VF_region
    call get_kphi_and_nodePhiTyp_VF_region
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_DS_trnsfrms
    call deallocate_and_reallocate_VF_region_trnsfrmVars
    call readjst_all_trnsfrms_VF_region
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) "AFT adding VF"
    
  end subroutine add_cellFor_VitellineFluid_region_and_redefine_system
  
  
  subroutine store_SysVars_and_Arrays_bfrAddingCell_for_VF_region
    implicit none
    
    node_xyStrNVFR       = node_xy
    node_typStrNVFR      = node_typ   
    node_cnnctdStrNVFR   = node_cnnctd 
    double_nodeStrNVFR   = double_node
    count_this_dnStrNVFR = count_this_dn
    
    typ_sprStrNVFR       = typ_spr
    k_sprStrNVFR         = k_spr
    l0_StrNVFR           = l0 
    l_StrNVFR            = l
    
    k_areaStrNVFR        = k_area
    A0_StrNVFR           = A0
    A_StrNVFR            = A
    
    nodePhi_typStrNVFR   = nodePhi_typ
    k_phiStrNVFR         = k_phi
    
    CgXNode_StrNVFR      = CgXNode
    CgYNode_StrNVFR      = CgYNode
    
    N_nodeSNVFR          = N_node
    N_sprSNVFR           = N_spr
    N_phiSNVFR           = N_phi
    N_cellSNVFR          = N_cell
    
    coordntes_xyStrNVFR  = coordntes_xy
    
  end subroutine store_SysVars_and_Arrays_bfrAddingCell_for_VF_region
  
  subroutine get_Sysvars_aftAddingCell_for_VF_region_WO_PrtlApclSgmntn
    implicit none
    
    N_cell = N_cell + 1 ; write(*,*)N_cell,"N_cell aft adding VF"
    N_spr  = N_spr  + 0 ; write(*,*)N_spr, "N_spr  aft adding VF"
    N_node = N_node + 0 ; write(*,*)N_node,"N_node  aft adding VF"    
    
  end subroutine get_Sysvars_aftAddingCell_for_VF_region_WO_PrtlApclSgmntn
  
  subroutine get_Sysvars_aftAddingCell_for_VF_region_WT_PrtlApclSgmntn
    implicit none
    integer :: N_apclSprCreated
    
    N_apclSideSgmntd = 3
    
    N_cell = N_cell + 1
    N_spr  = N_spr  + 1
    N_node = N_node     
    
  end subroutine get_Sysvars_aftAddingCell_for_VF_region_WT_PrtlApclSgmntn
  
  
  subroutine getDecsn_of_dealloctn_reallocatn_of_movingVars(lgcl_dealloc)
    implicit none
    logical, intent(out) :: lgcl_dealloc
    integer              :: N_mvCoordnteVal
    
    N_mvCoordnteVal = N_mvCoordnte
    call get_N_mvCoordnte
    
    write(*,*) N_mvCoordnteVal,N_mvCoordnte,"N_mvCoor bfr and aft"
    
    lgcl_dealloc = .True.
    
    if (N_mvCoordnteVal == N_mvCoordnte) lgcl_dealloc=.False.
    if (N_mvCoordnteVal.ne.N_mvCoordnte) lgcl_dealloc=.True.
    
  end subroutine getDecsn_of_dealloctn_reallocatn_of_movingVars
  
  
  
  subroutine change_num_of_InsrtdNodes_inApclSide_and_redefine_system()
    implicit none
    
    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    call deallocate_and_reallocate_transfrmStr_variables
    call deallocate_all_gradient_variables_StrVars
    call get_all_gradient_variables_StrVars
    
    call store_SysVars_and_Arrays_for_crtclDesgnd_apclSurfc
    call get_Sysvars_for_crtclDesgnd_apclSurfc
    call deallocate_and_reallocate_arrays_for_crtclDesgnd_apclSurfc
    
    call get_NodeVars_crtclDesgnd_apclSurfc
    call get_SprVars_crtclDesgnd_apclSurfc
    call get_CellVars_crtclDesgnd_apclSurfc
    call get_CgVars_crtclDesgnd_apclSurfc
    call get_kphi_and_nodePhiTyp_crtclDesgnd_apclSurfc
    call get_SRyp_props_crtclDesgnd_apclSurfc
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_anyTypPrv_trnsfrms
    call deallocate_and_reallocate_crtclDesgnd_apclSurfc_trnsfrmVars
    call readjst_all_trnsfrms_crtclDesgnd_apclSurfc
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    write(*,*) "AFT crtclDesgnd_apclSurfc "
    
  end subroutine change_num_of_InsrtdNodes_inApclSide_and_redefine_system
  
  subroutine  store_SysVars_and_Arrays_for_crtclDesgnd_apclSurfc
    implicit none
    
    call resolve_allocataion_problems_with_system_variabls
    
    node_xyStr          = node_xy
    node_typStr         = node_typ   
    node_cnnctdStr      = node_cnnctd 
    double_nodeStr      = double_node
    count_this_dnStr    = count_this_dn
    
    typ_sprStr          = typ_spr
    k_sprStr            = k_spr
    l0_Str              = l0 
    l_Str               = l
    
    k_areaStr           = k_area
    A0_Str              = A0
    A_Str               = A
    
    nodePhi_typStr      = nodePhi_typ
    k_phiStr            = k_phi
    
    CgXNode_Str         = CgXNode
    CgYNode_Str         = CgYNode
    
    AmpApclYP_Str       = AmpApclYP
    AlphaApclYP_Str     = AlphaApclYP
    epsApclYP_Str       = epsApclYP
    powrApclYP_Str      = powrApclYP
    activtnFctrApcl_Str = activtnFctrApcl
    
    N_nodeS             = N_node
    N_sprS              = N_spr
    N_phiS              = N_phi
    N_cellS             = N_cell
    
    numApclNodesS       = numApclNodes
    numBsalNodesS       = numBsalNodes
    numLtrlNodesS       = numLtrlNodes
    
    coordntes_xyStr     = coordntes_xy
    
  end subroutine store_SysVars_and_Arrays_for_crtclDesgnd_apclSurfc
  
  subroutine get_Sysvars_for_crtclDesgnd_apclSurfc
    implicit none
    
    N_cell = N_cell + 0                                                  ; write(*,*) N_cell,"N_cell crtclDsnd"
    N_spr  = N_spr  + ((2*NCP_CrtclApSrfc)+1)*(NAEC_ApclCrtcl-NAEC_Apcl) ; write(*,*) N_spr,"N_spr crtclDsnd"
    N_node = N_node + ((2*NCP_CrtclApSrfc)+1)*(NAEC_ApclCrtcl-NAEC_Apcl) ; write(*,*) N_node,"N_node crtclDsnd"
    
    call get_positional_Nodes_with_diff_types_of_InsrtdNodes ! apcl_bsal_ltrl_Nodes
    
    ! the Case here is = change_num_of_InsrtdNodes_inApclSide
    
  end subroutine get_Sysvars_for_crtclDesgnd_apclSurfc
  
  subroutine resolve_allocataion_problems_with_system_variabls
    implicit none
    
    call resolve_allocation_problems_with_node_variabls
    call resolve_allocation_problems_with_spr_variabls
    call resolve_allocation_problems_with_area_variabls
    call resolve_allocation_problems_with_bend_variabls
    call resolve_allocation_problems_with_grvtnl_variabls
    call resolve_allocation_problems_with_SRyp_variabls
    
  end subroutine resolve_allocataion_problems_with_system_variabls
  
  subroutine resolve_allocation_problems_with_node_variabls
    implicit none
    
    if (.not.allocated(node_xyStr)) then
       allocate(node_xyStr(1:N_node,1:N_dmnsn)) ; node_xyStr = -1.0d20
    else
       deallocate(node_xyStr)
       allocate(node_xyStr(1:N_node,1:N_dmnsn)) ; node_xyStr = -1.0d20
    endif
    
    if (.not.allocated(node_typStr)) then
       allocate(node_typStr(1:N_node)) ; node_typStr = -1
    else
       deallocate(node_typStr) ; allocate(node_typStr(1:N_node)) ;  node_typStr = -1
    endif
    
    if (.not.allocated(node_cnnctd)) then
       allocate(node_cnnctd(1:N_node)) ; node_cnnctd = -1
    else
       deallocate(node_cnnctd) ; allocate(node_cnnctd(1:N_node)) ;  node_cnnctd = -1
    endif
    
    if (.not.allocated(double_nodeStr)) then
       allocate(double_nodeStr(1:N_node,1:2)) ; double_nodeStr = -1
    else
       deallocate(double_nodeStr); allocate(double_nodeStr(1:N_node,1:2)) ; double_nodeStr = -1
    endif
    
    if (.not.allocated(count_this_dnStr)) then
       allocate(count_this_dnStr(1:N_node)) ; count_this_dnStr = -1
    else
       deallocate(count_this_dnStr) ; allocate(count_this_dnStr(1:N_node)) ;  count_this_dnStr = -1
    endif
    
  end subroutine resolve_allocation_problems_with_node_variabls
  
  subroutine resolve_allocation_problems_with_spr_variabls
    implicit none
    
    if (.not.allocated(typ_sprStr)) then
       allocate(typ_sprStr(1:N_spr)) ; typ_sprStr = -1
    else
       deallocate(typ_sprStr) ; allocate(typ_sprStr(1:N_spr)) ; typ_sprStr = -1
    endif
    
    if (.not.allocated(k_sprStr)) then
       allocate(k_sprStr(1:N_spr)) ; k_sprStr = -1.0d20
    else
       deallocate(k_sprStr) ; allocate(k_sprStr(1:N_spr)) ; k_sprStr = -1.0d20
    endif
    
    if (.not.allocated(l0_Str)) then
       allocate(l0_Str(1:N_spr)) ; l0_Str = -1.0d25
    else
       deallocate(l0_Str) ; allocate(l0_Str(1:N_spr)) ; l0_Str = -1.0d25
    endif
    
    if (.not.allocated(l_Str)) then
       allocate(l_Str(1:N_spr)) ; l_Str = -1.0d25
    else
       deallocate(l_Str) ; allocate(l_Str(1:N_spr)) ; l_Str = -1.0d25
    endif
    
  end subroutine resolve_allocation_problems_with_spr_variabls
  
  
  subroutine resolve_allocation_problems_with_area_variabls
    implicit none
    
    if (.not.allocated(k_areaStr)) then
       allocate(k_areaStr(1:N_cell)) ; k_areaStr = -1.0d20
    else
       deallocate(k_areaStr) ; allocate(k_areaStr(1:N_cell)) ; k_areaStr = -1.0d20
    endif
    
    if (.not.allocated(A0_Str)) then
       allocate(A0_Str(1:N_cell)) ; A0_Str = -1.0d25
    else
       deallocate(A0_Str) ; allocate(A0_Str(1:N_cell)) ; A0_Str = -1.0d25
    endif
    
    if (.not.allocated(A_Str)) then
       allocate(A_Str(1:N_cell)) ; A_Str = -1.0d25
    else
       deallocate(A_Str) ; allocate(A_Str(1:N_cell)) ; A_Str = -1.0d25
    endif
    
  end subroutine resolve_allocation_problems_with_area_variabls
  
  subroutine resolve_allocation_problems_with_bend_variabls
    implicit none
    
    if (.not.allocated(nodePhi_typStr)) then
       allocate(nodePhi_typStr(1:N_node)) ; nodePhi_typStr = -10
    else
       deallocate(nodePhi_typStr) ; allocate(nodePhi_typStr(1:N_node)) ;  nodePhi_typStr = -10
    endif
    
    if (.not.allocated(k_phiStr)) then
       allocate(k_phiStr(1:N_node,1:2)) ; k_phiStr = -1.0d30
    else
       deallocate(k_phiStr); allocate(k_phiStr(1:N_node,1:2)) ; k_phiStr = -1.0d30
    endif
    
  end subroutine resolve_allocation_problems_with_bend_variabls
  
  subroutine resolve_allocation_problems_with_grvtnl_variabls
    implicit none
    
    if (.not.allocated(CgXNode_Str)) then
       allocate(CgXNode_Str(1:N_node)) ; CgXNode_Str = -10
    else
       deallocate(CgXNode_Str) ; allocate(CgXNode_Str(1:N_node)) ;  CgXNode_Str = -10
    endif
    
    if (.not.allocated(CgYNode_Str)) then
       allocate(CgYNode_Str(1:N_node)) ; CgYNode_Str = -10
    else
       deallocate(CgYNode_Str) ; allocate(CgYNode_Str(1:N_node)) ;  CgYNode_Str = -10
    endif
    
  end subroutine resolve_allocation_problems_with_grvtnl_variabls
  
  subroutine resolve_allocation_problems_with_SRyp_variabls
    implicit none
    
    if (.not.allocated(AmpApclYP_Str)) then
       allocate(AmpApclYP_Str(1:numApclNodes)) ; AmpApclYP_Str(1:numApclNodes) = -1.0d30 
    else
       deallocate(AmpApclYP_Str) ; allocate(AmpApclYP_Str(1:numApclNodes)) ; AmpApclYP_Str(1:numApclNodes) = -1.0d30
    endif
    
    if (.not.allocated(AlphaApclYP_Str)) then
       allocate(AlphaApclYP_Str(1:numApclNodes)) ; AlphaApclYP_Str(1:numApclNodes) = -1.0d30 
    else
       deallocate(AlphaApclYP_Str);allocate(AlphaApclYP_Str(1:numApclNodes));AlphaApclYP_Str(1:numApclNodes)=-1.0d30
    endif
    
    if (.not.allocated(epsApclYP_Str)) then
       allocate(epsApclYP_Str(1:numApclNodes)) ; epsApclYP_Str(1:numApclNodes) = -1.0d30 
    else
       deallocate(epsApclYP_Str) ; allocate(epsApclYP_Str(1:numApclNodes)) ; epsApclYP_Str(1:numApclNodes) = -1.0d30
    endif
    
    if (.not.allocated(powrApclYP_Str)) then
       allocate(powrApclYP_Str(1:numApclNodes)) ; powrApclYP_Str(1:numApclNodes) = -1.0d30 
    else
       deallocate(powrApclYP_Str); allocate(powrApclYP_Str(1:numApclNodes)); powrApclYP_Str(1:numApclNodes) = -1.0d30
    endif
    
    if (.not.allocated(activtnFctrApcl_Str)) then
       allocate(activtnFctrApcl_Str(1:numApclNodes)) ; activtnFctrApcl_Str(1:numApclNodes) = -1 
    else
       deallocate(activtnFctrApcl_Str) ; allocate(activtnFctrApcl_Str(1:numApclNodes))
       activtnFctrApcl_Str(1:numApclNodes) = -1
    endif
    
  end subroutine resolve_allocation_problems_with_SRyp_variabls
  
  subroutine add_nodes_InApicalSides_and_redefine_system
    implicit none
    logical :: lgcl_dealloc=.False.
    integer :: i,j
    
    call getDecsn_of_dealloctn_reallocatn_of_movingVars(lgcl_dealloc);write(*,*) "lgcl_dealloc =",lgcl_dealloc
    call deallocate_moving_coordnte_variables
    call get_all_moving_coordnte_variables
    call deallocate_and_reallocate_transfrmStr_variables
    call deallocate_all_gradient_variables_StrVars
    call get_all_gradient_variables_StrVars
    
    !call store_SysVars_and_Arrays_bfr_AddingNode_inApclSide
    !call get_Sysvars_bfr_AddingNode_inApclSide
    !call deallocate_and_reallocate_arrays_bfr_AddingNode_inApclSide
    
    
  end subroutine add_nodes_InApicalSides_and_redefine_system
  
  subroutine store_SysVars_and_Arrays_bfr_AddingNode_inApclSide
    implicit none
    
    
  end subroutine store_SysVars_and_Arrays_bfr_AddingNode_inApclSide
  
  
  subroutine print_sprVars_crtclDesgnd_apclSurfc(lim1,lim2,lim3,lim4,lim5)
    implicit none
    integer, intent(in) :: lim1,lim2,lim3,lim4,lim5
    integer             :: i,j,cntCurrSys
    real*8              :: diffks,diffl0,diffl
    
    cntCurrSys = 1
    
    do i = 1,N_sprS
       
       if (i.le.lim1) then
          diffks = k_sprStr(i)-k_spr(cntCurrSys)
          diffl0 = l0_Str(i)-l0(cntCurrSys)
          diffl  = l_Str(i) - l(cntCurrSys)
          
          write(*,*) k_sprStr(i),k_spr(cntCurrSys),diffks,l0_Str(i),l0(cntCurrSys),diffl0,&
               l_Str(i),l(cntCurrSys),diffl,i,cntCurrSys,"lim1"
          cntCurrSys = cntCurrSys+1
          
       elseif ((i.gt.lim1).and.(i.le.lim2)) then   
          
          if ((i-lim1).le.(NAEC_Apcl+1)) then
             
             do j = 1,(addApNodesPerSgmntdApSpr+1)
                if (j==1)   write(*,*) k_sprStr(i),k_spr(cntCurrSys),l0_Str(i),l0(cntCurrSys),&
                     l_Str(i),l(cntCurrSys),i,cntCurrSys,"lim2"
                if (j.gt.1) write(*,*) k_spr(cntCurrSys),l0(cntCurrSys),l(cntCurrSys),i,cntCurrSys,"lim2"
                cntCurrSys = cntCurrSys+1
             enddo
             
          elseif ((i-lim1).gt.(NAEC_Apcl+1)) then
             
             diffks = k_sprStr(i)-k_spr(cntCurrSys)
             diffl0 = l0_Str(i)-l0(cntCurrSys)
             diffl  = l_Str(i) - l(cntCurrSys)
             
             write(*,*) k_sprStr(i),k_spr(cntCurrSys),diffks,l0_Str(i),l0(cntCurrSys),diffl0,&
                  l_Str(i),l(cntCurrSys),diffl,i,cntCurrSys,"lim2"
             cntCurrSys = cntCurrSys+1 
          endif
          
       elseif ((i.gt.lim2).and.(i.le.lim3)) then
          
          diffks = k_sprStr(i)-k_spr(cntCurrSys)
          diffl0 = l0_Str(i)-l0(cntCurrSys)
          diffl  = l_Str(i) - l(cntCurrSys)
          
          write(*,*) k_sprStr(i),k_spr(cntCurrSys),diffks,l0_Str(i),l0(cntCurrSys),diffl0,&
               l_Str(i),l(cntCurrSys),diffl,i,cntCurrSys,"lim3"
          cntCurrSys = cntCurrSys+1
          
       elseif ((i.gt.lim3).and.(i.le.lim4)) then
          
          if ((i-lim3).le.(NAEC_Apcl+1)) then
             
             do j = 1,(addApNodesPerSgmntdApSpr+1)
                if (j==1)   write(*,*) k_sprStr(i),k_spr(cntCurrSys),l0_Str(i),l0(cntCurrSys),&
                     l_Str(i),l(cntCurrSys),i,cntCurrSys,"lim4"
                if (j.gt.1) write(*,*) k_spr(cntCurrSys),l0(cntCurrSys),l(cntCurrSys),i,cntCurrSys,"lim4"
                cntCurrSys = cntCurrSys+1
             enddo
             
          elseif ((i-lim3).gt.(NAEC_Apcl+1)) then
             
             diffks = k_sprStr(i)-k_spr(cntCurrSys)
             diffl0 = l0_Str(i)-l0(cntCurrSys)
             diffl  = l_Str(i) - l(cntCurrSys)
             
             write(*,*) k_sprStr(i),k_spr(cntCurrSys),diffks,l0_Str(i),l0(cntCurrSys),diffl0,&
                  l_Str(i),l(cntCurrSys),diffl,i,cntCurrSys,"lim4"
             cntCurrSys = cntCurrSys+1 
          endif
          
          
       elseif ((i.gt.lim4).and.(i.le.lim5)) then
          
          if ((i-lim4).le.(NAEC_Apcl+1)) then
             
             do j = 1,(addApNodesPerSgmntdApSpr+1)
                if (j==1)   write(*,*) k_sprStr(i),k_spr(cntCurrSys),l0_Str(i),l0(cntCurrSys),&
                     l_Str(i),l(cntCurrSys),i,cntCurrSys,"lim5"
                if (j.gt.1) write(*,*) k_spr(cntCurrSys),l0(cntCurrSys),l(cntCurrSys),i,cntCurrSys,"lim5"
                cntCurrSys = cntCurrSys+1
             enddo
             
          elseif ((i-lim4).gt.(NAEC_Apcl+1)) then
             
             diffks = k_sprStr(i)-k_spr(cntCurrSys)
             diffl0 = l0_Str(i)-l0(cntCurrSys)
             diffl  = l_Str(i) - l(cntCurrSys)
             
             write(*,*) k_sprStr(i),k_spr(cntCurrSys),diffks,l0_Str(i),l0(cntCurrSys),diffl0,&
                  l_Str(i),l(cntCurrSys),diffl,i,cntCurrSys,"lim5"
             cntCurrSys = cntCurrSys+1 
          endif
       endif
       
    enddo
    
  end subroutine print_sprVars_crtclDesgnd_apclSurfc
  
  subroutine Energy_CHECK_withPrinting
    implicit none
    real*8 :: E,Ea,Es,Eg
    real*8 :: Es1,Ea1,Eg1,Eb1
    
    E = Energy(node_xy,l0,A0) ; write(*,*) E,"E1"
    
    Es=0.0d0 ; Ea=0.0d0 ; Eg=0.0d0
    
    do i = 1,N_spr
       Es = Es + (0.5000d0)*(k_spr(i)*(l(i)-l0(i))**2)
       write(*,*) i,k_spr(i),l(i),l0(i),Es,"Es"
    enddo
    
    do i = 1,N_cell
       Ea = Ea + (0.5000d0)*(k_area(i)*(A(i)-A0(i))**2)
       write(*,*) i,k_area(i),A(i),A0(i),Ea,"Ea"
    enddo
    
    do i = 1,N_node
       Eg = Eg + CgXNode(i)*(node_xy(i,1)) + CgYNode(i)*(node_xy(i,2))
       write(*,*) i,CgXNode(i),CgYNode(i),node_xy(i,1:2),Eg,"Eg"
    enddo 
    
    Es1 = spr_E(node_xy,l0)
    Ea1 = area_E(node_xy,A0)
    Eg1 = grvtnl_E(node_xy)
    Eb1 = bend_E(node_xy)
    
    write(*,*) Es1,Ea1,Eg1,Eb1,"Es1-Ea1-Eg1-Eb1" ; call sleep(1)
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E2"
    
  end subroutine Energy_CHECK_withPrinting
  
  
end module redefining_system_module
