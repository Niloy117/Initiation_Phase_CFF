module vitelline_fluid_calctn_mod
  use system_parameters
  use cell_info
  use transfrm_info
  use triangle_routines
  
  implicit none
  
contains
  
  real*8 function vitelline_fluid_region(area_nm,dum_Nodes)
    implicit none
    integer, intent(in) :: area_nm
    real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8              :: VF_area
    
    if (NI_incldAS_woPP == 0) call vitelline_fluid_region_methd1(VF_area,dum_Nodes)
    if (NI_incldAS_woPP == 1) call vitelline_fluid_region_methd5(VF_area,dum_Nodes)
    
    vitelline_fluid_region = VF_area
    
  end function vitelline_fluid_region
  
  subroutine vitelline_fluid_region_methd1(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: node1,node2,node3,node4
    real*8              :: N1(1:3),N2(1:3),N3(1:3),N4(1:3)
    real*8              :: sgmntd_triangl1,sgmntd_triangl2
    
    VF_area = 0.00d0 ; sgmntd_triangl1 = 0.00d0 ; sgmntd_triangl2 = 0.00d0
    
    node1 = ((N_cell-1)/2)*2 - 1 ; node2 = node1 + 2
    node3 = (N_cell*2)-1         ; node4 = node3 - 2
    
    write(*,*) node1,node2,node3,node4,N_cell,"nodes"
    write(*,*) dum_Nodes(node1,1:2),dum_Nodes(node1,1:2),"aaa"
    N1(1:3) = 0.0000d0 ; N2(1:3) = 0.0000d0 ; N3(1:3) = 0.0000d0 ; N4(1:3) = 0.0000d0
    N1(1:2) = dum_Nodes(node1,1:2) ; N2(1:2) = dum_Nodes(node2,1:2)
    N3(1:2) = dum_Nodes(node3,1:2) ; N4(1:2) = dum_Nodes(node4,1:2)
    
    sgmntd_triangl1 = calc_triangle(N1,N2,N3) ; write(*,*) sgmntd_triangl1,"st1"
    sgmntd_triangl2 = calc_triangle(N1,N3,N4) ; write(*,*) sgmntd_triangl2,"st2"
    
    VF_area = sgmntd_triangl1 + sgmntd_triangl2 ; write(*,*) VF_area,"vf"
    
  end subroutine vitelline_fluid_region_methd1
  
  
  subroutine vitelline_fluid_region_methd2(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    real*8              :: N1(1:3),N2(1:3),N3(1:3)
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node2Num,node3Num
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: MiddleOfEndNodes(1:2)
    real*8              :: sgmntd_triangl
    
    
    NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    MiddleOfEndNodes(1)       = 0.5000d0 * (EndNodesOFVFarea(1,1) + EndNodesOFVFarea(2,1)) 
    MiddleOfEndNodes(2)       = 0.5000d0 * (EndNodesOFVFarea(1,2) + EndNodesOFVFarea(2,2)) 
    
    !write(*,*) EndNodesOFVFarea(1,1:2),"End1"
    !write(*,*) EndNodesOFVFarea(2,1:2),"End2"
    !write(*,*) MiddleOfEndNodes(1:2),"Middle"
    
    N1(1:3) = 0.0000d0 ; N2(1:3) = 0.0000d0 ; N3(1:3) = 0.0000d0
    N1(1:2) = MiddleOfEndNodes(1:2)
    
    VF_area = 0.0000d0
    
    do i = 1,(NumNodesInVFarea-1)
       node2Num = i ; node3Num = i+1
       
       N2(1:2) = NodesOfVF(node2num,1:2)
       N3(1:2) = NodesOfVF(node3num,1:2)
       
       sgmntd_triangl = 0.00d0
       sgmntd_triangl = calc_triangle(N1,N2,N3)  !; write(*,*) sgmntd_triangl,i,"st i"    
       VF_area        = VF_area + sgmntd_triangl !; write(*,*) VF_area,i,"VF_area"
       
    enddo
    
  end subroutine vitelline_fluid_region_methd2
  
  subroutine vitelline_fluid_region_methd3(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node1Num,node2Num,nodeNum
    real*8              :: x1,y1,x2,y2,xval
    real*8              :: base,height
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: vrtically_sgmntd_area
    real*8              :: yTINYplus=+0.0000000000001d0
    real*8              :: yTINYminus=-0.0000000000001d0
    
    NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    VF_area = 0.0000d0
    
    do i = 1,(NumNodesInVFarea-1)
       
       node1Num = area_node(N_cell,i)
       node2Num = area_node(N_cell,i+1) !; write(*,*) node1Num,node2Num,"node1-2 Num"
       
       x1 = dum_Nodes(node1Num,1) ; x2 = dum_Nodes(node2Num,1) !; write(*,*) x1,x2,"x1-x2"
       y1 = dum_Nodes(node1Num,2) ; y2 = dum_Nodes(node2Num,2) !; write(*,*) y1,y2,"y1-y2"
       
       
       if (y1.gt.yTINYplus) then
          
          if (y2.gt.yTINYplus) then
             vrtically_sgmntd_area = 0.0000d0
             
          elseif (y2.lt.yTINYminus) then
             
             call find_the_xval_on_yZero(x1,y1,x2,y2,xVal)
             base = abs(x2-xval) ; height = abs(y2)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             vrtically_sgmntd_area = 0.0000d0
          endif
          
          
       elseif (y1.lt.yTINYminus) then
          
          if (y2.gt.yTINYplus) then
             
             call find_the_xval_on_yZero(x1,y1,x2,y2,xVal)
             base = abs(x1-xval) ; height = abs(y1)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
             
          elseif (y2.lt.yTINYminus) then
             vrtically_sgmntd_area = 0.5000d0 * (abs(y1)+abs(y2)) * abs(x2-x1)
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             base = abs(x2-x1) ; height = abs(y1)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height) 
             
          endif
          
          
       elseif ((y1.gt.yTINYminus) .and. (y1.lt.yTINYplus)) then
          
          if (y2.gt.yTINYplus) then
             
             vrtically_sgmntd_area = 0.0000d0
             
          elseif (y2.lt.yTINYminus) then
             base = abs(x2-x1) ; height = abs(y2)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height) 
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             
             vrtically_sgmntd_area = 0.0000d0
             
          endif
          
          
          
       endif
       
       VF_area = VF_area + vrtically_sgmntd_area !; write(*,*) VF_area,i,"VF_area"
       
    enddo
    
  end subroutine vitelline_fluid_region_methd3
  
  
  subroutine vitelline_fluid_region_methd3_write(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8,  intent(in) :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node1Num,node2Num,nodeNum
    real*8              :: x1,y1,x2,y2,xval
    real*8              :: base,height
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: vrtically_sgmntd_area
    real*8              :: yTINYplus=+0.0000000000001d0
    real*8              :: yTINYminus=-0.0000000000001d0
    real*8              :: VF_areaBfr
    
    open(unit=142,file='vitelline_fluid_region_methd3_write1.dat')
    
    do i = 1,area_node(N_cell,0)
       nodeNum = area_node(N_cell,i) ; write(142,*) nodeNum,i,"NodeNum" 
       write(142,*) dum_Nodes(nodeNum,1:2),nodeNum,"node_xyV"
    enddo
    write(142,*) " "
    
    NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    write(142,*) NumNodesInVFarea,area_node(N_cell,0),"NumNodeInVF equiv test"
    
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    VF_area = 0.0000d0
      
    do i = 1,(NumNodesInVFarea-1)
       
       node1Num = area_node(N_cell,i)
       node2Num = area_node(N_cell,i+1) ; write(142,*) node1Num,node2Num,"node1-2 Num"
       
       x1 = dum_Nodes(node1Num,1) ; x2 = dum_Nodes(node2Num,1) ; write(142,*) x1,x2,"x1-x2"
       y1 = dum_Nodes(node1Num,2) ; y2 = dum_Nodes(node2Num,2) ; write(142,*) y1,y2,"y1-y2" ; write(142,*) " "
       
       
       if (y1.gt.yTINYplus) then
            
          if (y2.gt.yTINYplus) then
             vrtically_sgmntd_area = 0.0000d0
             write(142,*) y1,y2,x1,x2,i,"cond 11"
             
          elseif (y2.lt.yTINYminus) then
             
             call find_the_xval_on_yZero_write(x1,y1,x2,y2,xVal)
             base = abs(x2-xval) ; height = abs(y2)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
             write(142,*) y1,y2,x1,x2,i,"cond 12"
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             vrtically_sgmntd_area = 0.0000d0
             write(142,*) y1,y2,x1,x2,i,"cond 13"
          endif
          
          
       elseif (y1.lt.yTINYminus) then
          
          if (y2.gt.yTINYplus) then
             call find_the_xval_on_yZero_write(x1,y1,x2,y2,xVal)
             base = abs(x1-xval) ; height = abs(y1)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
             write(142,*) y1,y2,x1,x2,i,"cond 21"
             
          elseif (y2.lt.yTINYminus) then
             vrtically_sgmntd_area = 0.5000d0 * (abs(y1)+abs(y2)) * abs(x2-x1)
             write(142,*) y1,y2,x1,x2,i,"cond 22"
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             base = abs(x2-x1) ; height = abs(y1)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height) 
             write(142,*) y1,y2,x1,x2,i,"cond 23"
          endif
          
          
       elseif ((y1.gt.yTINYminus) .and. (y1.lt.yTINYplus)) then
          
          if (y2.gt.yTINYplus) then
             vrtically_sgmntd_area = 0.0000d0
             write(142,*) y1,y2,x1,x2,i,"cond 31"
             
          elseif (y2.lt.yTINYminus) then
             base = abs(x2-x1) ; height = abs(y2)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height) 
             write(142,*) y1,y2,x1,x2,i,"cond 32"
             
          elseif ((y2.lt.yTINYplus) .and. (y2.gt.yTINYminus)) then
             vrtically_sgmntd_area = 0.0000d0
             write(142,*) y1,y2,x1,x2,i,"cond 33"
          endif
          
       endif
       
       VF_areaBfr = VF_area
       VF_area    = VF_area + vrtically_sgmntd_area
       write(142,*) VF_areaBfr,VF_area,i,"VF_area"
       
    enddo
    
  end subroutine vitelline_fluid_region_methd3_write
  
  
  subroutine vitelline_fluid_region_methd5(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8, intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node1Num,node2Num,nodeNum
    real*8              :: x1,y1,x2,y2,xval
    real*8              :: base,sumOfPLlines
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: vrtically_sgmntd_area
    
    !NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    
    if (VF_regionModelled==1 .and. NI_AS_wtCrtclSurfc==0) then
       NumNodesInVFarea = 2*(1+addedNCPair)+(1+(2*addedNCPair))*(NAEC_Apcl)
    elseif (VF_regionModelled==1 .and. NI_AS_wtCrtclSurfc==1) then
       NumNodesInVFarea = 2*(1+addedNCPair) + (1+(2*addedNCPair))*(NAEC_Apcl) + &
            ((NCP_CrtclApSrfc*2)+1)*(NAEC_ApclCrtcl-NAEC_Apcl)
    endif
    
    !write(*,*) NumNodesInVFarea,addedNCPair,NCP_CrtclApSrfc,NAEC_ApclCrtcl,NAEC_Apcl,"check it"
    !write(*,*) area_node(N_cell,0),"aft chk it"
    
    NumNodesInVFarea = area_node(N_cell,0)
    
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    VF_area = 0.0000d0
    
    do i = 1,(NumNodesInVFarea-1)
       node1Num = area_node(N_cell,i)
       node2Num = area_node(N_cell,i+1)
       x1       = dum_Nodes(node1Num,1) ; x2 = dum_Nodes(node2Num,1)
       y1       = dum_Nodes(node1Num,2) ; y2 = dum_Nodes(node2Num,2)
       base     = (x2-x1) ; sumOfPLlines = -(y1+y2)
       
       vrtically_sgmntd_area = 0.5000d0*(base)*(sumOfPLlines)      
       VF_area               = VF_area + vrtically_sgmntd_area
    enddo
    
  end subroutine vitelline_fluid_region_methd5
  
  
  
  subroutine vitelline_fluid_region_methd5_write(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8, intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node1Num,node2Num,nodeNum
    real*8              :: x1,y1,x2,y2,xval
    real*8              :: base,sumOfPLlines
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: vrtically_sgmntd_area
    
    
    open(unit=143,file='vitelline_fluid_region_methd5_write.dat')
    
    do i = 1,area_node(N_cell,0)
       nodeNum = area_node(N_cell,i) ; write(143,*) nodeNum,i,"NodeNum" 
       write(143,*) dum_Nodes(nodeNum,1:2),nodeNum,"node_xyV"
    enddo
    write(143,*) " "
    
    
    NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    VF_area = 0.0000d0
    
    do i = 1,(NumNodesInVFarea-1)
       
       node1Num = area_node(N_cell,i)
       node2Num = area_node(N_cell,i+1) ; write(143,*) node1Num,node2Num,"node1-2 Num"
       
       x1 = dum_Nodes(node1Num,1) ; x2 = dum_Nodes(node2Num,1) ; write(143,*) x1,x2,"x1-x2"
       y1 = dum_Nodes(node1Num,2) ; y2 = dum_Nodes(node2Num,2) ; write(143,*) y1,y2,"y1-y2"
       
       base = (x2-x1) ; sumOfPLlines = -(y1+y2) ; write(143,*) base,sumOfPLlines,"base-SumOfPLlines"
       vrtically_sgmntd_area = 0.5000d0*(base)*(sumOfPLlines) 
             
       VF_area = VF_area + vrtically_sgmntd_area ; write(143,*) vrtically_sgmntd_area,VF_area,i,"VS,VF_area"
       write(143,*) " "
    enddo
    
    close(143)
    
  end subroutine vitelline_fluid_region_methd5_write
  
  
  subroutine vitelline_fluid_region_methd4(VF_area,dum_Nodes)
    implicit none
    real*8, intent(out) :: VF_area
    real*8, intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    
    integer             :: NumNodesInVFarea,NodeNm,cntNode
    integer             :: k,kmax,i
    integer             :: node1Num,node2Num
    real*8              :: x1,y1,x2,y2,xval
    real*8              :: base,height
    real*8, allocatable :: NodesOfVF(:,:)
    real*8              :: EndNodesOFVFarea(1:2,1:2)
    real*8              :: vrtically_sgmntd_area
    real*8              :: yTINY=0.0000000000001d0
    
    integer             :: nodeNum
      
    open(unit=142,file='vitelline_fluid_region_methd4.dat')
    
    do i = 1,area_node(N_cell,0)
       nodeNum = area_node(N_cell,i) ; write(142,*) nodeNum,i,"NodeNum" 
       write(*,*) dum_nodes(nodeNum,1:2),nodeNum,"node_xyV"
    enddo
    
    NumNodesInVFarea = 2 + (NAEC_Apcl) + 4*(1+(NAEC_Apcl))
    allocate(NodesOfVF(1:NumNodesInVFarea,1:2)) ; NodesOfVF = -1.0d30
    write(142,*) NumNodesInVFarea,area_node(i,0),"NumNodeInVF equiv test"
    
    call nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    
    VF_area = 0.0000d0
    
    do i = 1,(NumNodesInVFarea-1)
         
       node1Num = area_node(N_cell,i)
       node2Num = area_node(N_cell,i+1) ; write(142,*) node1Num,node2Num,"node1-2 Num"
       
       x1 = dum_Nodes(node1Num,1) ; x2 = dum_Nodes(node2Num,1) ; write(142,*) x1,x2,"x1-x2"
       y1 = dum_Nodes(node1Num,2) ; y2 = dum_Nodes(node2Num,2) ; write(142,*) y1,y2,"y1-y2"
       
       if (y1.gt.yTINY) then
          
          if (y2.lt.yTINY) then
             call find_the_xval_on_yZero(x1,y1,x2,y2,xVal)
             base = abs(x2-xval) ; height = abs(y2)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
             
             write(142,*) base,height,vrtically_sgmntd_area,"BHV"
             
          elseif (y2.gt.yTINY) then
             vrtically_sgmntd_area = 0.0000d0
          endif
          
       elseif (y1.lt.yTINY) then
          
          if (y2.lt.yTINY) then
             vrtically_sgmntd_area = 0.5000d0 * (abs(y1)+abs(y2)) * abs(x2-x1)
          elseif (y2.gt.yTINY) then
             call find_the_xval_on_yZero(x1,y1,x2,y2,xVal)
             base = abs(x1-xval) ; height = abs(y1)
             vrtically_sgmntd_area = 0.5000d0*(base)*(height)
          endif
          
       endif
       
       VF_area = VF_area + vrtically_sgmntd_area ; write(142,*) VF_area,i,"VF_area"
       
    enddo
    
    close(142)
    
  end subroutine vitelline_fluid_region_methd4
  
  
  subroutine nodes_of_VF(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    implicit none
    integer, intent(in)  :: NumNodesINVFarea
    real*8,  intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 , intent(out) :: nodesOfVF(1:NumNodesINVFarea,1:2)
    real*8 , intent(out) :: EndNodesOFVFarea(1:2,1:2)
    
    if (area_node(N_cell,0).le.0) then
       call nodes_of_VF_manualApproach(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    elseif (area_node(N_cell,0).gt.0) then
       call nodes_of_VF_automtdApproach(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    endif
    
  end subroutine nodes_of_VF
  
  subroutine nodes_of_VF_manualApproach(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    implicit none
    integer, intent(in)  :: NumNodesINVFarea
    real*8,  intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 , intent(out) :: nodesOfVF(1:NumNodesINVFarea,1:2)
    real*8 , intent(out) :: EndNodesOFVFarea(1:2,1:2)
    
    integer              :: NodeNm,cntNode
    integer              :: k,kmax
    
    write(*,*) "HAVE TO ADJUST IT for more added Node Case"
    
    cntNode                = 1
    NodeNm                 = ((N_cell-1)/2)*2 - 1  ; write(*,*) NodeNm,"NodeNm p1"
    NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
    cntNode                = cntNode+1
    
    EndNodesOFVFarea(1,1:2)= dum_Nodes(NodeNm,1:2) 
    
    kmax = NAEC_Apcl+1
    
    do k = 1,kmax
       NodeNm                 = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k
       write(*,*) NodeNm,"NodeNm p2"
       NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
       cntNode                = cntNode+1
    enddo
    
    NodeNm                 = ((N_cell-1)/2)*2 + 1  ; write(*,*) NodeNm,"NodeNm p3"
    NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
    cntNode                = cntNode+1
    
    do k = 1,kmax
       NodeNm                 = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell*2)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k
       write(*,*) NodeNm,"NodeNm p4"
       NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
       cntNode                = cntNode+1
    enddo
    
    NodeNm                 = (N_cell*2)-1        ; write(*,*) NodeNm,"NodeNm p5"  
    NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
    cntNode                = cntNode+1
    
    do k = 1,kmax
       NodeNm                 = 2*(Hlf_Ncell+1)*2+((Hlf_Ncell*2)-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+k
       write(*,*) NodeNm,"NodeNm p6"
       NodesOfVF(cntNode,1:2) = dum_Nodes(NodeNm,1:2)
       cntNode                = cntNode+1
    enddo
    
    NodeNm                    = (N_cell*2)-3        ; write(*,*) NodeNm,"NodeNm p7"
    NodesOfVF(cntNode,1:2)    = dum_Nodes(NodeNm,1:2)
    
    EndNodesOFVFarea(2,1:2)   = dum_Nodes(NodeNm,1:2)
    
  end subroutine nodes_of_VF_manualApproach
  
  
  subroutine nodes_of_VF_automtdApproach(NumNodesINVFarea,dum_Nodes,NodesOfVF,EndNodesOFVFarea)
    implicit none
    integer, intent(in)  :: NumNodesINVFarea
    real*8,  intent(in)  :: dum_Nodes(1:N_node,1:N_dmnsn)
    real*8 , intent(out) :: nodesOfVF(1:NumNodesINVFarea,1:2)
    real*8 , intent(out) :: EndNodesOFVFarea(1:2,1:2)
    
    integer :: i,j
    integer :: NodeNm
    
    if (area_node(N_cell,0).ne.NumNodesINVFarea) then
       write(*,*) area_node(N_cell,0),NumNodesINVFarea,"area_node(N_cell,0),NumNodesINVFarea"
       stop 'NumNodesINVFarea=/area_node(N_cell,0)'
    endif
    
    do i = 1,NumNodesINVFarea
       NodeNm = area_node(N_cell,i)
       if (i==1)                EndNodesOFVFarea(1,1:2)  = dum_Nodes(NodeNm,1:2)
       if (i==NumNodesINVFarea) EndNodesOFVFarea(2,1:2)  = dum_Nodes(NodeNm,1:2)
       nodesOfVF(i,1:2) = dum_Nodes(NodeNm,1:2)
    enddo
    
  end subroutine nodes_of_VF_automtdApproach
  
  subroutine find_the_xval_on_yZero(x1,y1,x2,y2,xVal)
    implicit none
    real*8, intent(in)  :: x1,y1,x2,y2
    real*8, intent(out) :: xVal
    real*8              :: slopeV,cval_fn,cval_sn,diffInCval
    
    slopeV  = (y2-y1)/(x2-x1)
    cval_fn = y1 - (slopeV)*x1
    cval_sn = y2 - (slopeV)*x2
    
    diffInCval =  (cval_fn-cval_sn)/(cval_fn*100.0000d0)
    xVal       = -(cval_fn)/(slopeV)
    
  end subroutine find_the_xval_on_yZero
  
  subroutine find_the_xval_on_yZero_write(x1,y1,x2,y2,xVal)
    implicit none
    real*8, intent(in)  :: x1,y1,x2,y2
    real*8, intent(out) :: xVal
    real*8              :: slopeV,cval_fn,cval_sn,diffInCval
    
    slopeV  = (y2-y1)/(x2-x1)
    cval_fn = y1 - (slopeV)*x1
    cval_sn = y2 - (slopeV)*x2
    
    diffInCval =  (cval_fn-cval_sn)/(cval_fn*100.0000d0)
    xVal       = -(cval_fn)/(slopeV)
    
    write(142,*) slopeV,cval_fn,cval_sn,diffInCval,xval,"slope-cval"
    
  end subroutine find_the_xval_on_yZero_write
  
  
end module vitelline_fluid_calctn_mod
