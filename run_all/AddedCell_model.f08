module Adding_cells
  use nodeInsrtd_cntrlStates
  implicit none
  integer :: N_addedCellL,N_addedCellR
  integer :: N_addedCell
  integer :: N_addedSprL,N_addedSprR
  integer :: N_addedSpr
  
  integer :: N_TNsAdded,Hlf_NtnsAdded ! TN = Terminal Node
  integer :: N_INsAdded ! IN = Inserted Node
  integer :: Hlf_TNsNAC,Hlf_INsNAC
  integer :: N_TNsNAC,N_TNsAC
  integer :: N_lstTN_NAC,Hlf_NtnNAC,N_lstTN_AC,Hlf_NtnAC
  integer :: last_INLftAC
  
  integer :: nsprs_InACell,nsprs_InEndCell
  
  real*8 , allocatable :: node_xyNAC(:,:)
  integer, allocatable :: node_typNAC(:),node_cnnctdNAC(:)
  integer, allocatable :: double_nodeNAC(:,:),count_this_dnNAC(:)
  
  integer, allocatable :: typ_sprNAC(:)
  real*8 , allocatable :: k_sprNAC(:) ,l0_NAC(:),l_NAC(:)
  real*8 , allocatable :: k_areaNAC(:),A0_NAC(:),A_NAC(:)
  
  real*8 , allocatable :: nodePhi_typNAC(:),k_phiNAC(:,:)
  real*8 , allocatable :: CgXNodeNAC(:),CgYNodeNAC(:)
  
  integer :: N_nodeNAC,N_sprNAC,N_cellNAC,Hlf_NsprNAC,Hlf_NcellNAC
  integer :: N_mvCoordnteNAC,N_mvCoordnte_withl0A0NAC
  
  real*8,  allocatable :: coordntesNAC(:),coordntes_xyNAC(:)
  integer, allocatable :: coordnte_of_which_nodeNAC(:),x_or_yNAC(:)
  
  integer :: allctn_cnt=0
  
  integer :: max_node_sprNAC ,max_node_areaNAC
  integer :: max_spr_nodeNAC ,max_spr_areaNAC
  integer :: max_area_nodeNAC,max_area_sprNAC
  
  integer :: allctn_cntTrns=0,allctn_cntTrnsBfrMet=0
  
  integer, allocatable :: node_sprNAC(:,:)  , node_areaNAC(:,:)
  integer, allocatable :: spr_nodeNAC(:,:)  , spr_areaNAC(:,:)
  integer, allocatable :: area_nodeNAC(:,:) , area_sprNAC(:,:)
  real*8 , allocatable :: NlistNAC(:,:,:)
  
  real*8, allocatable  :: gradientNAC(:,:),gradient_anaNAC(:,:),gradient_numNAC(:,:)
  real*8, allocatable  :: grd_mvNAC(:),grdmv_xyNAC(:)
  
  integer :: NcellAddedVrsn
  integer :: N_addedCycl,NcellPropToBeRead,NsprPropToBeRead
  
  real*8, allocatable :: ksCS1(:),ksCS2(:),ksCS3(:),ksCS4(:)
  real*8, allocatable :: l0CS1(:),l0CS2(:),l0CS3(:),l0CS4(:)
  real*8, allocatable :: kaCS1(:),kaCS2(:),kaCS3(:),kaCS4(:)
  real*8, allocatable :: A0CS1(:),A0CS2(:),A0CS3(:),A0CS4(:)
  real*8, allocatable :: CgXCS1(:),CgXCS2(:),CgXCS3(:),CgXCS4(:)
  
  integer :: FirstNonAddedCellL,FirstNonAddedCellR

  integer, allocatable :: node_sprBfrMet(:,:),node_areaBfrMet(:,:)
  integer, allocatable :: spr_nodeBfrMet(:,:),spr_areaBfrMet(:,:)
  integer, allocatable :: area_nodeBfrMet(:,:),area_sprBfrMet(:,:)
  integer, allocatable :: NlistBfrMet(:,:,:)
  
  integer             :: N_CS
  integer             :: allctnCycl_cnt=0
  real*8, allocatable :: ka_cycl(:,:,:),A0_cycl(:,:,:),ks_cycl(:,:,:),l0_cycl(:,:,:)
  
contains
  
  subroutine switchto_ACmodel_and_mayOrmaynt_switchback
    implicit none
    real*8  :: E,Es,Esp,Ea,Eap
    integer :: i,j,jmax,nodeNm
    integer :: numOfCellsToAdd
    
    call get_gradient(node_xy,l0,A0,gradient)
    call Find_Analytical_and_Numerical_Mismatch
    
    numOfCellsToAdd=4
    call switchto_AC_model(numOfCellsToAdd) !AC=Added Cell model, NAC=NonAddedCell model

    AddedCellModelInitiation = 1
    E                        = Energy(node_xy,l0,A0) ; write(*,*) E,"E in sb: switchto_ACmodel_and_mayOrmaynt"
    
    !call Equilibrate_AddedCell_model
    
    FirstNonAddedCellL = 1 + N_addedCellL
    FirstNonAddedCellR = 1 + Hlf_Ncell + N_addedCellR
    call make_AddedCellIdenticaltoFirstNonAddedCell
    
    !stop
    
    call PropValueReadForAddedCell
    
    if (AC_NACswitchback==1) then 
       !call deallocate_repetitive_arrays_ACmodel
       !call switchback_to_NAC_model
       write(*,*) "switcback_to_NAC not written yet"
       stop
    elseif (AC_NACswitchback==0) then
       continue
    endif
    
  end subroutine switchto_ACmodel_and_mayOrmaynt_switchback
  
  
  subroutine switchto_ACmodel_prptycopyFrmFirstIC_and_Equilibrate(numOfCellsToAdd)
    implicit none
    integer, intent(in) :: numOfCellsToAdd
    
    real*8  :: E,Es,Esp,Ea,Eap
    integer :: i,j,jmax,nodeNm
    
    call get_gradient(node_xy,l0,A0,gradient)
    call Find_Analytical_and_Numerical_Mismatch
    
    call switchto_AC_model(numOfCellsToAdd) !AC=Added Cell model, NAC=NonAddedCell model
    
    AddedCellModelInitiation = 1
    E                        = Energy(node_xy,l0,A0) ; write(*,*) E,"E in sb: switchto_ACmodel_and_mayOrmaynt"
    
    FirstNonAddedCellL = 1 + N_addedCellL
    FirstNonAddedCellR = 1 + Hlf_Ncell + N_addedCellR
    call make_AddedCellIdenticaltoFirstNonAddedCell
    
    call Equilibrate_AddedCell_model
    
  end subroutine switchto_ACmodel_prptycopyFrmFirstIC_and_Equilibrate
  
  subroutine switchto_AC_model(numOfCellsToAdd) !AC=Added Cell method
    implicit none
    integer, intent(in) :: numOfCellsToAdd
    integer :: i
    
    if (modelID==1) stop 'AC models are written for only NI system in switchto_AC_model'
    
    call store_NAC_SysVars_and_Arrays
    call get_AC_SysVars(numOfCellsToAdd)
    call deallocate_and_reallocate_NAC_Arrays
    
    call get_NodeVars_ForACmodel
    call get_AreaVars_ForACmodel
    call get_SprVars_ForACmodel
    call get_kphi_and_nodePhiTyp_ForACmodel
    call get_CgVars_ForACmodel
    
    call store_all_moving_coordnte_variables_NAC
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_NAC_trnsfrms
    call deallocate_and_reallocate_AC_trnsfrmVars
    
    do i = 1,N_lstTN_NAC
       write(*,*) node_sprNAC(i,0:max_spr_node),i,"node_sprNAC 1"
    enddo
    
    open(unit=112,file='TrnsfrmDataChk.dat')
    
    write(112,*) N_node,N_spr,N_cell,Hlf_Ncell
    write(112,*) N_nodeNAC,N_sprNAC,N_cellNAC
    write(112,*) Hlf_NsprNAC,Hlf_NcellNAC
    
    write(112,*) N_TNsNAC,N_TNsAC
    write(112,*) N_lstTN_NAC,Hlf_NtnNAC,N_lstTN_AC,Hlf_NtnAC
    write(112,*) last_INLftAC
    
    write(112,*) N_addedCellL,N_addedCellR
    write(112,*) nsprs_InACell
    
    write(112,*) NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl
    write(112,*) max_spr_node,max_area_node
    write(112,*) max_node_spr,max_area_spr
    write(112,*) max_node_area,max_spr_area
    
    close(112)
    
    
    call update_all_the_transfrms
    call print_all_NAC_trnsfrms
    call print_all_AC_trnsfrms_automated
    
    call updt_A0l0_usingTrnsfrms
    
    call store_gradient_system_variablesNAC
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    
  end subroutine switchto_AC_model
  
  
  subroutine store_NAC_SysVars_and_Arrays
    implicit none
    
    if (allctn_cnt==0) then
       
       allocate(node_xyNAC(1:N_node,1:N_dmnsn),node_typNAC(1:N_node))
       allocate(node_cnnctdNAC(1:N_node),double_nodeNAC(1:N_node,1:2),count_this_dnNAC(1:N_node))
       allocate(typ_sprNAC(1:N_spr),k_sprNAC(1:N_spr),l0_NAC(1:N_spr),l_NAC(1:N_spr))
       allocate(k_areaNAC(1:N_cell),A0_NAC(1:N_cell),A_NAC(1:N_cell))
       
       allocate(nodePhi_typNAC(1:N_node),k_phiNAC(1:N_node,1:max_Phi_node))
       allocate(CgXNodeNAC(1:N_node),CgYNodeNAC(1:N_node))
       
       allocate(gradientNAC(1:N_node,1:N_dmnsn),gradient_anaNAC(1:N_node,1:N_dmnsn))
       allocate(gradient_numNAC(1:N_node,1:N_dmnsn))
       allocate(grd_mvNAC(1:N_mvCoordnte_withl0A0),grdmv_xyNAC(1:N_mvCoordnte))
       
       allctn_cnt = 1
       
    else
       continue
    endif
    
    node_xyNAC(1:N_node,1:N_dmnsn)  = node_xy(1:N_node,1:N_dmnsn)
    node_typNAC(1:N_node)           = node_typ(1:N_node)
    node_cnnctdNAC(1:N_node)        = node_cnnctd(1:N_node)
    double_nodeNAC(1:N_node,1:2)    = double_node(1:N_node,1:2)
    count_this_dnNAC(1:N_node)      = count_this_dn(1:N_node)
    
    typ_sprNAC(1:N_spr) = typ_spr(1:N_spr)
    k_sprNAC(1:N_spr)   = k_spr(1:N_spr)
    l0_NAC(1:N_spr)     = l0(1:N_spr)
    l_NAC(1:N_spr)      = l(1:N_spr)
    
    k_areaNAC(1:N_cell) = k_area(1:N_cell)
    A0_NAC(1:N_cell)    = A0(1:N_cell)
    A_NAC(1:N_cell)     = A(1:N_cell)
    
    nodePhi_typNAC(1:N_node)          = nodePhi_typ(1:N_node)
    k_phiNAC(1:N_node,1:max_Phi_node) = k_phi(1:N_node,1:max_Phi_node)
    CgXNodeNAC(1:N_node)              = CgXNode(1:N_node)
    CgYNodeNAC(1:N_node)              = CgYNode(1:N_node)
    
    N_nodeNAC       = N_node
    N_sprNAC        = N_spr
    N_cellNAC       = N_cell
    
    Hlf_NsprNAC     = (N_spr-NAEC_Ltrl-1)/2
    Hlf_NcellNAC    = N_cell/2
    
  end subroutine store_NAC_SysVars_and_Arrays
  
  subroutine get_AC_SysVars(numOfCellsToAdd)
    implicit none
    integer, intent(in) :: numOfCellsToAdd
    integer             :: sprAdded
    
    N_addedCellL   = numOfCellsToAdd
    N_addedCellR   = numOfCellsToAdd
    
    N_addedCell    = N_addedCellL + N_addedCellR
    N_cell         = N_cell + N_addedCell
    Hlf_Ncell      = N_cell/2
    NcellAddedVrsn = N_cell
    
    nsprs_InACell   = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    nsprs_InEndCell = (NAEC_Ltrl+1)
    
    N_addedSprL   = (N_addedCellL)*(nsprs_InACell)
    N_addedSprR   = (N_addedCellR)*(nsprs_InACell)
    N_addedSpr    = N_addedSprL + N_addedSprR
    
    if (CellsMeet == 0) N_lstTN_NAC = 2*(Hlf_NcellNAC+1)*2     ! N_nodeS
    if (CellsMeet.gt.0) N_lstTN_NAC = 2*(Hlf_NcellNAC+1)*2 + 1
    
    N_TNsNAC      = N_lstTN_NAC
    N_lstTN_AC    = (N_lstTN_NAC) + (N_addedCellL+N_addedCellR)*2
    
    N_TNsAdded   = N_addedCell*2
    N_INsAdded   = (NAEC_apcl+NAEC_bsal+NAEC_ltrl)*N_addedCell
    N_node       = N_node + N_TNsAdded + N_INsAdded 
    
    sprAdded     = ((NAEC_apcl+1) + (NAEC_bsal+1) + (NAEC_ltrl+1))*N_addedCell
    N_spr        = N_spr + sprAdded
    
    last_INLftAC = N_lstTN_AC + ((N_node-2-N_lstTN_AC)/2)
    
    write(*,*) N_addedCellL,N_addedCellR,N_addedCell,N_cell,Hlf_Ncell,NcellAddedVrsn,"print vars1"
    write(*,*) nsprs_InACell,nsprs_InEndCell,N_addedSprL,N_addedSprR,N_addedSpr,"print vars2"
    write(*,*) N_lstTN_NAC,Hlf_NcellNAC,N_TNsNAC,"print vars3"
    write(*,*) N_lstTN_AC,N_TNsAdded,N_INsAdded,N_node,sprAdded,N_spr,last_INLftAC,"print vars4"
    write(*,*) N_nodeS,"N_nodeS"
    
  end subroutine get_AC_SysVars
  
  subroutine deallocate_and_reallocate_NAC_Arrays
    implicit none
    
    deallocate(node_xy,node_typ)
    deallocate(node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(k_area,A0,A)
    
    deallocate(nodePhi_typ,k_phi)
    deallocate(CgXNode,CgYNode)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_area_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_NAC_Arrays
  
  subroutine get_NodeVars_ForACmodel
    implicit none
    real*8  :: prvLTopxy(1:2),prvRTopxy(1:2)
    real*8  :: prvLBotxy(1:2),prvRBotxy(1:2)
    integer :: prvLTopNT,prvRTopNT !NT=NodeType
    integer :: prvLBotNT,prvRBotNT
    
    integer :: ACL_cnt,ACR_cnt
    integer :: cnt_INsL,cnt_INsR
    integer :: i,j,jmax,m,m1,m2
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: cellNm,LftOrRght
    integer :: apcl_End(1:2),bsal_End(1:2),ltrl_End(1:2)
    integer :: NAEC_val
    integer :: N1,N2
    real*8  :: node1(1:2),node2(1:2)
    integer :: End_addedIN
    integer :: nodeNm
    integer :: st1,st2,end1,end2
    
    real*8, allocatable :: IntrmdNodes(:,:)
    
    open(unit=110,file='get_NodeVars_AC.dat')
    
    Hlf_NtnsAdded = (N_TNsAdded/2)   !L/R both equal that why
    Hlf_TNsNAC    = (N_lstTN_NAC)/2  !N_nodeS/2
    Hlf_INsNAC    = (N_nodeNAC-N_lstTN_NAC)/2
    
    write(110,*) Hlf_NtnsAdded,Hlf_TNsNAC,Hlf_INsNAC,"Hlf_NtnsAdded,Hlf_TNsNAC,Hlf_INsNAC"
    write(*,*) Hlf_NtnsAdded,Hlf_TNsNAC,Hlf_INsNAC,"Hlf_NtnsAdded,Hlf_TNsNAC,Hlf_INsNAC"
    
    N_TNsAC = (N_lstTN_NAC) + (N_TNsAdded) !Num of TNs in AddedCell
    write(110,*) N_lstTN_NAC,N_TNsAdded,N_TNsAC,"N_lstTN_NAC,N_TNsAdded,N_TNsAC"
    
    
    prvLTopxy(1:2) = node_xyNAC(1,1:N_dmnsn)
    prvLBotxy(1:2) = node_xyNAC(2,1:N_dmnsn)
    prvRTopxy(1:2) = node_xyNAC((Hlf_TNsNAC+1),1:N_dmnsn)
    prvRBotxy(1:2) = node_xyNAC((Hlf_TNsNAC+2),1:N_dmnsn)
    
    prvLTopNT = node_typNAC(1)
    prvLBotNT = node_typNAC(2)
    prvRTopNT = node_typNAC(Hlf_TNsNAC+1)
    prvRBotNT = node_typNAC(Hlf_TNsNAC+2)
    
    write(110,*) prvLTopxy(1:2),prvLTopNT,"LTop"
    write(110,*) prvLBotxy(1:2),prvLBotNT,"LBot"
    write(110,*) prvRTopxy(1:2),prvRTopNT,"RTop"
    write(110,*) prvRBotxy(1:2),prvRBotNT,"LBot"
    
    ACL_cnt = 0 ; ACR_cnt = 0
    
    lim1 = Hlf_NtnsAdded
    lim2 = lim1 + Hlf_TNsNAC
    lim3 = lim2 + Hlf_NtnsAdded
    lim4 = lim3 + Hlf_TNsNAC
    lim5 = lim4 + 1
    
    write(110,*) lim1,lim2,lim3,lim4,lim5,"lims"
    
    do i = 1,N_TNsAC
       
       if (i.le.lim1) then
          
          if (mod(i,2)==1) then
             node_xy(i,1) = prvLTopxy(1) - (N_addedCellL-ACL_cnt)*(La1) 
             node_xy(i,2) = prvLTopxy(2)
             
             if ((i/2)==0)   node_typ(i)  = prvLTopNT
             if ((i/2).gt.0) node_typ(i)  = 2
             
          elseif (mod(i,2)==0) then
             node_xy(i,1) = prvLBotxy(1) - (N_addedCellL-ACL_cnt)*(La1)
             node_xy(i,2) = prvLBotxy(2)

             if ((i/2)==1)   node_typ(i)  = prvLBotNT
             if ((i/2).gt.1) node_typ(i)  = 1
             
             ACL_cnt = ACL_cnt+1
          endif
          
          write(110,*) i,node_xy(i,1:2),node_typ(i),"i"
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          node_xy(i,1:2) = node_xyNAC((i-Hlf_NtnsAdded),1:2)
          node_typ(i)    = node_typNAC(i-Hlf_NtnsAdded)
          
          write(110,*) i,node_xy(i,1:2),node_typ(i),"i"
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          if (mod(i,2)==1) then
             node_xy(i,1) = prvRTopxy(1) + (N_addedCellR-ACR_cnt)*(La1) 
             node_xy(i,2) = prvRTopxy(2)
             
             if (((i-lim2)/2)==0)   node_typ(i)  = prvRTopNT
             if (((i-lim2)/2).gt.0) node_typ(i)  = prvRTopNT
             
          elseif (mod(i,2)==0) then
             node_xy(i,1) = prvRBotxy(1) + (N_addedCellR-ACR_cnt)*(La1)
             node_xy(i,2) = prvRBotxy(2)
             
             if (((i-lim2)/2)==1)   node_typ(i)  = prvRBotNT
             if (((i-lim2)/2).gt.1) node_typ(i)  = 1
             
             ACR_cnt = ACR_cnt+1
          endif
          
          write(110,*) i,node_xy(i,1:2),node_typ(i),"i"
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          node_xy(i,1:2) = node_xyNAC((i-N_TNsAdded),1:2)
          node_typ(i)    = node_typNAC(i-N_TNsAdded)
          
          write(110,*) i,node_xy(i,1:2),node_typ(i),"i"
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          node_xy(i,1:2) = node_xyNAC((i-N_TNsAdded),1:2)
          node_typ(i)    = node_typNAC(i-N_TNsAdded)

          write(110,*) i,node_xy(i,1:2),node_typ(i),"i"
          
       endif
       
    enddo
    
    cnt_INsL = 1 ; cnt_INsR = 1

    write(*,*) NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,"NAEC_val"

    do i = 1,2
       
       if (i==1) jmax = N_addedCellL
       if (i==2) jmax = N_addedCellR
       
       do j = 1,jmax
          
          if (i==1) then
             cellNm = j ; LftOrRght = 1
             call get_TN_valsFor_INsXY(cellNm,LftorRght,apcl_End,bsal_End,ltrl_End)
             
          elseif (i==2) then
             cellNm = Hlf_Ncell + j ; LftOrRght = 2
             call get_TN_valsFor_INsXY(cellNm,LftOrRght,apcl_End,bsal_End,ltrl_End)
          endif
          
          do m = 1,3 ! 3 for apcl,bsal,ltrl
             
             if (m==1) then
                NAEC_val = NAEC_apcl
                N1 = apcl_End(1)
                N2 = apcl_End(2)
                
             elseif (m==2) then
                NAEC_val = NAEC_bsal
                N1 = bsal_End(1)
                N2 = bsal_End(2)
                
             elseif (m==3) then
                NAEC_val = NAEC_ltrl
                N1 = ltrl_End(1)
                N2 = ltrl_End(2)
             endif
             
             if (NAEC_val.gt.0) then
                node1(1:2) = node_xy(N1,1:2)
                node2(1:2) = node_xy(N2,1:2)
                
                allocate(IntrmdNodes(1:NAEC_val,1:2)) ; IntrmdNodes = -1.0d30
                call get_allIntrmdNodes(node1,node2,NAEC_val,IntrmdNodes)
                
                do m1 = 1,NAEC_val
                   if (i==1) nodeNm = N_TNsAC + cnt_INsL
                   if (i==2) nodeNm = N_TNsAC + Hlf_INsNAC-1 + (cnt_INsL-1) + cnt_INsR
                   
                   node_xy(nodeNm,1:2) = IntrmdNodes(m1,1:2)
                   node_typ(nodeNm)    = 1
                   
                   if (i==1) cnt_INsL = cnt_INsL+1
                   if (i==2) cnt_INsR = cnt_INsR+1
                   
                enddo
                
                deallocate(IntrmdNodes)
             endif
             
          enddo
          
       enddo
       
       End_addedIN = nodeNm
       write(110,*) End_addedIN,"EndAddedIN"
       write(*,*) End_addedIN,Hlf_INsNAC,"End_addedIN"
       
       if (i==1) then
          node_xy((End_addedIN+1):(End_addedIN+Hlf_INsNAC-1),1:2) = node_xyNAC((N_lstTN_NAC+1):(N_lstTN_NAC+Hlf_INsNAC-1),1:2)
          node_typ(End_addedIN+1:End_addedIN+Hlf_INsNAC-1)    = 1

          do m2 = (End_addedIN+1),(End_addedIN+Hlf_INsNAC-1)
             write(110,*) m2,node_xy(m2,1:2),node_typ(m2),"m2"
          enddo
          
       elseif (i==2) then
          
          node_xy((End_addedIN+1):(End_addedIN+Hlf_INsNAC-1),1:2) = node_xyNAC((N_lstTN_NAC+Hlf_INsNAC-1+1):(N_nodeNAC-2),1:2)
          node_typ(End_addedIN+1:End_addedIN+Hlf_INsNAC-1)    = 1

          do m2 = (End_addedIN+1),(End_addedIN+Hlf_INsNAC-1)
             write(110,*) m2,node_xy(m2,1:2),node_typ(m2),"m2"
          enddo
          
       endif
       
       
       if (i==2) then ! This one is for Inserted nodes at the last Cell
          st1 = N_node-NAEC_Ltrl+1    ; end1 = N_node
          st2 = N_nodeNAC-NAEC_Ltrl+1 ; end2 = N_nodeNAC
          write(*,*) st1,end1,st2,end2,"st-end"
          
          node_xy(st1:end1,1:2) = node_xyNAC(st2:end2,1:2)
          node_typ((N_node-NAEC_Ltrl+1):N_node) = 1

          do m2 = st1,st2
             write(110,*) node_xy(m2,1:2),node_typ(m2),m2,"m2"
          enddo
          
       endif
       
    enddo
    
    open(unit=68,file='nodesInSepFl.dat')
    
    do i = 1,N_node
       !write(68,*) i,node_xy(i,1:2),node_typ(i),"nodes"
       write(68,*) node_xy(i,1:2)
    enddo
    
    close(68)
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
    close(110)
    
  end subroutine get_NodeVars_ForACmodel
  
  subroutine get_TN_valsFor_INsXY(cellNm,LftorRght,apcl_End,bsal_End,ltrl_End)
    implicit none
    integer, intent(in)  :: cellNm,LftorRght
    integer, intent(out) :: apcl_End(1:2),bsal_End(1:2),ltrl_End(1:2)
    
    if (LftorRght == 1) then
       apcl_End(1) = 1 + (cellNm-1)*2 ; apcl_End(2) = 3 + (cellNm-1)*2 
       bsal_End(1) = 2 + (cellNm-1)*2 ; bsal_End(2) = 4 + (cellNm-1)*2
       ltrl_End(1) = 3 + (cellNm-1)*2 ; ltrl_End(2) = 4 + (cellNm-1)*2
       
    elseif (LftorRght == 2) then
       apcl_End(1) = 1 + (cellNm)*2   ; apcl_End(2) = 3 + (cellNm)*2
       bsal_End(1) = 2 + (cellNm)*2   ; bsal_End(2) = 4 + (cellNm)*2
       ltrl_End(1) = 3 + (cellNm)*2   ; ltrl_End(2) = 4 + (cellNm)*2 
    endif
    
  end subroutine get_TN_valsFor_INsXY
  
  
  subroutine get_AreaVars_ForACmodel
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    
    lim1 = N_addedCellL
    lim2 = lim1 + Hlf_NcellNAC
    lim3 = lim2 + N_addedCellR
    lim4 = lim3 + Hlf_NcellNAC
    lim5 = N_cell
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          k_area(i) = -1.0d30
          A0(i)     = -1.0d30
          A(i)      = -1.0d25
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          k_area(i) = k_areaNAC(i-N_addedCellL)
          A0(i)     = A0_NAC(i-N_addedCellL)
          A(i)      = A_NAC(i-N_addedCellL)
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          k_area(i) = -1.0d30
          A0(i)     = -1.0d30
          A(i)      = -1.0d25
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          k_area(i) = k_areaNAC(i-N_addedCell)
          A0(i)     = A0_NAC(i-N_addedCell)
          A(i)      = A_NAC(i-N_addedCell)
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          k_area(i) = k_areaNAC(i-N_addedCell)
          A0(i)     = A0_NAC(i-N_addedCell)
          A(i)      = A_NAC(i-N_addedCell)
       endif
       
    enddo
    
    
  end subroutine get_AreaVars_ForACmodel
  
  
  subroutine get_SprVars_ForACmodel
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: i,j,jmax
    integer :: sprNm,sprNmofNAC
    integer :: sprOfFirstCellLNAC,sprOfSecndCellLNAC
    integer :: sprOfFirstCellRNAC,sprOfSecndCellRNAC
    integer :: sprNmRdctn
    
    lim1 = N_addedCellL
    lim2 = lim1 + Hlf_NcellNAC
    lim3 = lim2 + N_addedCellR
    lim4 = lim3 + Hlf_NcellNAC
    lim5 = N_cell
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,"lims in get_SprVars_ForACmodel"
    
    do i = 1,N_cell
       
       if (i.ne.N_cell) jmax = nsprs_InACell
       if (i == N_cell) jmax = nsprs_InEndCell
       
       do j = 1,jmax
          
          sprNm = (i-1)*nsprs_InACell + j
          
          if (i.le.lim1) then
             k_spr(sprNm) = -1.0d30
             l0(sprNm)    = -1.0d30
             l(sprNm)     = -1.0d25
             
             if (i==1) then
                sprOfFirstCellLNAC = j
                typ_spr(sprNm)     = typ_sprNAC(sprNm)
                
             elseif (i.gt.1) then
                sprOfSecndCellLNAC = 1*(nsprs_InACell) + j
                typ_spr(sprNm)     = typ_sprNAC(sprOfSecndCellLNAC)
             endif
             
          elseif ((i.gt.lim1) .and. (i.le.lim2)) then
             
             sprNmRdctn   = N_addedCellL*nsprs_InACell
             
             k_spr(sprNm) = k_sprNAC(sprNm-sprNmRdctn)
             l0(sprNm)    = l0_NAC(sprNm-sprNmRdctn)
             l(sprNm)     = l_NAC(sprNm-sprNmRdctn)
             
             if ((i-lim1)==1) then
                sprOfSecndCellLNAC = 1*(nsprs_InACell) + j
                typ_spr(sprNm)    = typ_sprNAC(sprOfSecndCellLNAC)
                
             elseif ((i-lim1).gt.1) then
                sprNmofNAC     = (i-N_addedCellL-1) + j
                typ_spr(sprNm) = typ_sprNAC(sprNmofNAC) 
             endif
          
          elseif ((i.gt.lim2) .and. (i.le.lim3)) then
             k_spr(sprNm) = -1.0d30
             l0(sprNm)    = -1.0d30
             l(sprNm)     = -1.0d25
             
             if ((i-lim2)==1) then
                sprOfFirstCellRNAC = (Hlf_NcellNAC)*(nsprs_InACell) + j
                typ_spr(sprNm)     = typ_sprNAC(sprOfFirstCellRNAC)
                
             elseif ((i-lim2).gt.1) then
                sprOfSecndCellRNAC = (Hlf_NcellNAC+1)*(nsprs_InACell) + j
                typ_spr(sprNm)     = typ_sprNAC(sprOfSecndCellRNAC)
             endif
             
          elseif ((i.gt.lim3) .and. (i.le.lim4)) then
             
             sprNmRdctn   = N_addedCell*nsprs_InACell
             
             k_spr(sprNm) = k_sprNAC(sprNm-sprNmRdctn)
             l0(sprNm)    = l0_NAC(sprNm-sprNmRdctn)
             l(sprNm)     = l_NAC(sprNm-sprNmRdctn)
             
             if ((i-lim3)==1) then
                sprOfSecndCellRNAC = (Hlf_NcellNAC+1)*(nsprs_InACell) + j
                typ_spr(sprNm)    = typ_sprNAC(sprOfSecndCellRNAC)
                
             elseif ((i-lim3).gt.1) then
                sprNmofNAC     = (i-N_addedCell-1)*(nsprs_InACell) + j
                typ_spr(sprNm) = typ_sprNAC(sprNmofNAC) 
             endif
             
          elseif ((i.gt.lim4) .and. (i.le.lim5)) then
             
             sprNmRdctn   = N_addedCell*nsprs_InACell
             
             k_spr(sprNm) = k_sprNAC(sprNm-sprNmRdctn)
             l0(sprNm)    = l0_NAC(sprNm-sprNmRdctn)
             l(sprNm)     = l_NAC(sprNm-sprNmRdctn)
             
             sprNmofNAC     = (i-N_addedCell-1)*(nsprs_InACell) + j
             typ_spr(sprNm) = typ_sprNAC(sprNmofNAC)
             
          endif
          
          write(*,*) sprNm,k_spr(sprNm),l0(sprNm),l(sprNm),"sprNm-ks-l0-l"
       enddo
       
    enddo
    
  end subroutine get_SprVars_ForACmodel
  
  subroutine get_kphi_and_nodePhiTyp_ForACmodel
    implicit none
    integer :: i
    real*8  :: k_phiVal
    
    k_phiVal = k_phiNAC(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal
    
    do i = 1,N_node
       if (i.le.N_TNsAC) then
          nodePhi_typ(i) = 1
       elseif ((i.gt.N_TNsAC) .and. (i.le.N_node)) then
          nodePhi_typ(i) = 2
       endif
    enddo
    
    do i=1,N_node
       if (nodePhi_typ(i)==2) then
          k_phi(i,1:max_Phi_node) = 0.00d0
       endif
    enddo
    
  end subroutine get_kphi_and_nodePhiTyp_ForACmodel
  
  
  subroutine get_CgVars_ForACmodel
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5,lim6
    integer :: i
    
    lim1 = Hlf_NtnsAdded
    lim2 = lim1 + Hlf_TNsNAC
    lim3 = lim2 + Hlf_NtnsAdded
    lim4 = lim3 + Hlf_TNsNAC
    lim5 = lim4 + 1
    lim6 = N_node

    open(unit=55,file='CgNodeAC.dat')
    
    write(55,*) lim1,lim2,lim3,lim4,lim5,lim6,"lims Cg"
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          CgXNode(i) = 0.0d0
          CgYNode(i) = 0.0d0
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          CgXNode(i) = CgXNodeNAC(i-Hlf_NtnsAdded)
          CgYNode(i) = CgYNodeNAC(i-Hlf_NtnsAdded)
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          CgXNode(i) = 0.0d0
          CgYNode(i) = 0.0d0
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          CgXNode(i) = CgXNodeNAC(i-N_TNsAdded)
          CgYNode(i) = CgYNodeNAC(i-N_TNsAdded)
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          CgXNode(i) = CgXNodeNAC(i-N_TNsAdded)
          CgYNode(i) = CgYNodeNAC(i-N_TNsAdded)
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          CgXNode(i) = 0.0d0
          CgYNode(i) = 0.0d0
       endif
       
       write(55,*) CgXNode(i),CgYNode(i),i,"CgXNode-YNode,i"
    enddo
    
    close(55)
    
  end subroutine get_CgVars_ForACmodel
  
  subroutine store_all_moving_coordnte_variables_NAC
    implicit none

    N_mvCoordnteNAC           = N_mvCoordnte
    N_mvCoordnte_withl0A0NAC  = N_mvCoordnte_withl0A0
    
    coordntesNAC    = coordntes
    coordntes_xyNAC = coordntes_xy
    
    coordnte_of_which_nodeNAC = coordnte_of_which_node
    x_or_yNAC                 = x_or_y 
    
    
  end subroutine store_all_moving_coordnte_variables_NAC
  
  
  subroutine store_NAC_trnsfrms
    implicit none

    max_node_sprNAC  = max_node_spr
    max_node_areaNAC = max_node_area
    max_spr_nodeNAC  = max_spr_node
    max_spr_areaNAC  = max_spr_area
    max_area_nodeNAC = max_area_node
    max_area_sprNAC  = max_area_spr
    
    if (allctn_cntTrns==0) then
       allocate(node_sprNAC(1:N_nodeNAC,0:max_spr_nodeNAC) ,node_areaNAC(1:N_nodeNAC,0:max_area_nodeNAC))
       allocate(spr_nodeNAC(1:N_sprNAC,0:max_node_sprNAC) ,spr_areaNAC(1:N_sprNAC,0:max_area_sprNAC))
       allocate(area_nodeNAC(1:N_cellNAC,0:max_node_areaNAC),area_sprNAC(1:N_cellNAC,0:max_spr_areaNAC))
       allocate(NlistNAC(1:N_nodeNAC,1:max_area_nodeNAC,1:N_dmnsn))
       allctn_cntTrns = 1
    endif
    
    node_sprNAC  = node_spr
    node_areaNAC = node_area   
    spr_nodeNAC  = spr_node
    spr_areaNAC  = spr_area
    area_nodeNAC = area_node
    area_sprNAC  = area_spr
    
    NlistNAC = Nlist
    
    call print_all_NAC_trnsfrms
    
  end subroutine store_NAC_trnsfrms
  
  subroutine deallocate_and_reallocate_AC_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = max_node_sprNAC !did not change in AC and NAC
    max_node_area = max_node_areaNAC
    max_spr_node  = max_spr_nodeNAC
    max_spr_area  = max_spr_areaNAC
    max_area_node = max_area_nodeNAC
    max_area_spr  = max_area_sprNAC
    
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
    
  end subroutine deallocate_and_reallocate_AC_trnsfrmVars
  
  
  subroutine update_all_the_transfrms
    implicit none
    
    !call update_node_spr_readVrsn
    !call update_node_area_readVrsn
    !call update_spr_node_readVrsn
    !call update_spr_area_readVrsn
    !call update_area_node_readVrsn
    !call update_area_spr_readVrsn
    
    call read_node_spr_FrmNAC
    call read_node_area_FrmNAC
    call update_node_other_AddedCell
    
    call read_spr_node_FrmNAC
    call read_spr_area_FrmNAC
    call update_spr_other_AddedCell
    
    call read_area_node_FrmNAC
    call read_area_spr_FrmNAC
    call update_area_other_AddedCell
    
    call get_Nlist
    
  end subroutine update_all_the_transfrms
  
  subroutine read_node_spr_FrmNAC
    implicit none
    integer :: i,node_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='node_spr_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='node_spr_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='node_spr_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='node_spr_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='node_spr_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='node_spr_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='node_spr_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='node_spr_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='node_spr_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='node_spr_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_node_spr_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='node_sprWOtxtNAC.dat')
    endif
    
    do i = 1,N_nodeNAC
       read(67,*) node_nm,node_sprNAC(i,0:max_spr_node)
    enddo
    
    close(67)
    
  end subroutine read_node_spr_FrmNAC
  
  subroutine read_node_area_FrmNAC
    implicit none
    integer :: i,node_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='node_area_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='node_area_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='node_area_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='node_area_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='node_area_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='node_area_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='node_area_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='node_area_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='node_area_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='node_area_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_node_area_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='node_areaWOtxtNAC.dat')
    endif
        
    do i = 1,N_nodeNAC
       read(67,*) node_nm,node_areaNAC(i,0:max_area_node)
    enddo
    
    close(67)
    
  end subroutine read_node_area_FrmNAC
  
  subroutine read_spr_node_FrmNAC
    implicit none
    integer :: i,spr_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='spr_node_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='spr_node_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='spr_node_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='spr_node_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='spr_node_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='spr_node_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='spr_node_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='spr_node_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='spr_node_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='spr_node_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_spr_node_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='spr_nodeWOtxtNAC.dat')
    endif
    
    
    do i = 1,N_sprNAC
       read(67,*) spr_nm,spr_nodeNAC(i,0:max_node_spr)
    enddo
    
    close(67)
    
  end subroutine read_spr_node_FrmNAC
  
  subroutine read_spr_area_FrmNAC
    implicit none
    integer :: i,spr_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='spr_area_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='spr_area_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='spr_area_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='spr_area_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='spr_area_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='spr_area_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='spr_area_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='spr_area_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='spr_area_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='spr_area_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_spr_area_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='spr_areaWOtxtNAC.dat')
    endif
    
    open(unit=67,file='spr_areaWOtxtNAC.dat')
    
    do i = 1,N_sprNAC
       read(67,*) spr_nm,spr_areaNAC(i,0:max_area_spr)
    enddo
    
    close(67)
    
  end subroutine read_spr_area_FrmNAC
  
  
  subroutine read_area_node_FrmNAC
    implicit none
    integer :: i,area_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='area_node_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='area_node_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='area_node_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='area_node_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='area_node_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='area_node_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='area_node_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='area_node_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='area_node_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='area_node_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_area_node_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='area_nodeWOtxtNAC.dat')
    endif
    
    
    do i = 1,N_cellNAC
       read(67,*) area_nm,area_nodeNAC(i,0:max_node_area)
    enddo
    
    close(67)
    
  end subroutine read_area_node_FrmNAC
  
  
  subroutine read_area_spr_FrmNAC
    implicit none
    integer :: i,area_nm
    
    if (modelID==1) then   
       if (CellsMeet==0) open(unit=67,file='area_spr_TNCM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='area_spr_TNCM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='area_spr_TNCM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='area_spr_TNCM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='area_spr_TNCM4_WOtxtNAC.dat')
       
    elseif (modelID==2) then
       if (CellsMeet==0) open(unit=67,file='area_spr_NICM0_WOtxtNAC.dat')
       if (CellsMeet==1) open(unit=67,file='area_spr_NICM1_WOtxtNAC.dat')
       if (CellsMeet==2) open(unit=67,file='area_spr_NICM2_WOtxtNAC.dat')
       if (CellsMeet==3) open(unit=67,file='area_spr_NICM3_WOtxtNAC.dat')
       if (CellsMeet==4) open(unit=67,file='area_spr_NICM4_WOtxtNAC.dat')
       
    else
       write(*,*) 'openning a previous file, chk if that okay, sb: read_area_spr_FrmNAC'
       write(*,*) CellsMeet,modelID,"cellsmeet and modelID"
       open(unit=67,file='area_sprWOtxtNAC.dat')
    endif
    
    
    do i = 1,N_cellNAC
       read(67,*) area_nm,area_sprNAC(i,0:max_spr_area)
    enddo
    
    close(67)
    
  end subroutine read_area_spr_FrmNAC
  
  
  subroutine update_node_other_AddedCell
    implicit none
    integer :: dnStrtsL,dnStrtsR
    integer :: lastNodeLbfrDN,lastNodeRbfrDN
    integer :: i,j,jStr,jEnd,jmax,jmaxS,jmaxA
    integer :: prvOddNode,prvEvnNode
    integer :: node_equiv,modVal
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7
    integer :: cnt_NI,cnt_NIA
    integer :: currSprStrt,prvSprStrt,First_spr
    integer :: spr_cnctdAct,sprInNAC
    integer :: node_cnctdAct,areaInNAC,area_cnt
    
    Hlf_NtnNAC = N_lstTN_NAC/2 ; Hlf_NtnAC = N_lstTN_AC/2
    
    write(*,*) Hlf_NtnNAC,Hlf_NtnAC,"Hlf_NtnNAC,Hlf_NtnAC"
    write(*,*) N_lstTN_NAC,N_lstTN_AC,"N_lstTN"
    
    
    do i = 1,N_lstTN_NAC
       write(*,*) node_sprNAC(i,0:max_spr_node),i,"node_sprNAC"
    enddo
    
    do i = 1,2
       
       if (i==1) then
          jStr = 1              ; jEnd = Hlf_NtnNAC
       elseif (i==2) then
          jStr = (Hlf_NtnNAC+1) ; jEnd = N_lstTN_NAC
       endif
       
       write(*,*) jStr,jEnd,"jStr,jEnd"
       
       do j = jStr,jEnd
          
          spr_cnctdAct = node_sprNAC(j,0) ; write(*,*) spr_cnctdAct,j,max_spr_node,"spr_cnctdAct,j"
          
          if (spr_cnctdAct == max_spr_node) then
             if (i==1) dnStrtsL = j
             if (i==2) dnStrtsR = j
             exit
          endif
          
       enddo
       
    enddo
    
    lastNodeLbfrDN = dnStrtsL-1 ! last node left  before Double Node
    lastNodeRbfrDN = dnStrtsR-1 ! last node right before Double Node
    
    write(*,*) lastNodeLbfrDN,lastNodeRbfrDN,dnStrtsL,dnStrtsR,"lastNodeLbfrDN,lastNodeRbfrDN"
    
    lim1 = lastNodeLbfrDN
    lim2 = lim1 + (N_addedCellL*2)
    lim3 = Hlf_NtnAC
    lim4 = lastNodeRbfrDN + (N_addedCellL*2)
    lim5 = lim4 + (N_addedCellR*2)
    lim6 = N_TNsAC
    lim7 = N_node
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,lim7,"lims in update_node_other_AddedCell"
    
    cnt_NI   = 1
    cnt_NIA  = 1
    area_cnt = 0
    
    do i = 1,N_node
       
       if (i.le.lim1) then
          node_spr(i,0:max_spr_node)   = node_sprNAC(i,0:max_spr_node)
          node_area(i,0:max_area_node) = node_areaNAC(i,0:max_area_node)
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if (mod(i,2)==1) then
             prvOddNode = i-2
             
             jmaxS          = node_spr(prvOddNode,0)
             node_spr(i,0)  = node_spr(prvOddNode,0)
             
             do j = 1,jmaxS
                node_spr(i,j) = node_spr(prvOddNode,j) + nsprs_InACell 
             enddo
             
             jmaxA           = node_area(prvOddNode,0)
             node_area(i,0)  = node_area(prvOddNode,0)
             
             do j = 1,jmaxA
                node_area(i,j) = node_area(prvOddNode,j) + 1 
             enddo
             
          elseif (mod(i,2)==0) then
             prvEvnNode = i-2
             
             jmaxS         = node_spr(prvEvnNode,0)
             node_spr(i,0) = node_spr(prvEvnNode,0)
             
             do j = 1,jmaxS
                node_spr(i,j) = node_spr(prvEvnNode,j) + nsprs_InACell
             enddo
             
             jmaxA          = node_area(prvEvnNode,0)
             node_area(i,0) = node_area(prvEvnNode,0) 
             
             do j = 1,jmaxA
                node_area(i,j) = node_area(prvEvnNode,j) + 1
             enddo
             
          endif
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          node_equiv = i-(N_addedCellL*2)
          
          node_spr(i,0) = node_sprNAC(node_equiv,0)
          jmaxS         = node_sprNAC(node_equiv,0)
          
          do j = 1,jmaxS
             sprInNAC = node_sprNAC(node_equiv,j)
             
             if (sprInNAC .le. Hlf_NsprNAC) then
                node_spr(i,j) = node_sprNAC(node_equiv,j) + (N_addedCellL*nsprs_InACell)
             elseif (sprInNAc .gt. Hlf_NsprNAC) then
                node_spr(i,j) = node_sprNAC(node_equiv,j) + (N_addedCellL+N_addedCellR)*nsprs_InACell
             endif
          enddo
          
          node_area(i,0)   = node_areaNAC(node_equiv,0)
          jmaxA            = node_areaNAC(node_equiv,0)
          
          do j = 1,jmaxA
             areaInNAC = node_areaNAC(node_equiv,j)
             
             if (areaInNAC .le. Hlf_NcellNAC) node_area(i,j) =node_areaNAC(node_equiv,j)+(N_addedCellL)
             if (areaInNAC .gt. Hlf_NcellNAC) node_area(i,j) =node_areaNAC(node_equiv,j)+(N_addedCellL+N_addedCellR)
             
          enddo
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          node_equiv = i - (N_addedCellL*2)
          
          node_spr(i,0) = node_sprNAC(node_equiv,0)
          jmaxS         = node_sprNAC(node_equiv,0)
          
          do j = 1,jmaxS
             node_spr(i,j) = node_sprNAC(node_equiv,j) + (N_addedCellL*nsprs_InACell)
          enddo
          
          node_area(i,0)   = node_areaNAC(node_equiv,0)
          jmaxA            = node_areaNAC(node_equiv,0)
          
          do j = 1,jmaxA
             node_area(i,j) = node_areaNAC(node_equiv,j) + N_addedCellL
          enddo
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          if (mod(i,2)==1) then
             prvOddNode = i-2
             
             jmaxS         = node_spr(prvOddNode,0)
             node_spr(i,0) = node_spr(prvOddNode,0)
             
             do j = 1,jmaxS
                node_spr(i,j) = node_spr(prvOddNode,j) + nsprs_InACell 
             enddo
             
             jmaxA          = node_area(prvOddNode,0)
             node_area(i,0) = node_area(prvOddNode,0)
             
             do j = 1,jmaxA
                node_area(i,j) = node_area(prvOddNode,j) + 1
             enddo
             
             
          elseif (mod(i,2)==0) then
             prvEvnNode = i-2
             
             jmaxS         = node_spr(prvEvnNode,0)
             node_spr(i,0) = node_spr(prvEvnNode,0)
             
             do j = 1,jmaxS
                node_spr(i,j) = node_spr(prvEvnNode,j) + nsprs_InACell
             enddo
             
             
             jmaxA          = node_area(prvEvnNode,0)
             node_area(i,0) = node_area(prvEvnNode,0)
             
             do j = 1,jmaxA
                node_area(i,j) = node_area(prvEvnNode,j) + 1
             enddo
             
          endif
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          
          node_equiv = i-(N_addedCellL+N_addedCellR)*2
          
          node_spr(i,0) = node_sprNAC(node_equiv,0)
          jmaxS         = node_sprNAC(node_equiv,0)
          
          do j = 1,jmaxS
             sprInNAC      = node_sprNAC(node_equiv,j)
             
             if (sprInNAC .le. Hlf_NsprNAC) then
                node_spr(i,j) = node_sprNAC(node_equiv,j) + (N_addedCellL*nsprs_InACell)
                
             elseif (sprInNAC .gt. Hlf_NsprNAC) then
                node_spr(i,j) = node_sprNAC(node_equiv,j) + (N_addedCellL+N_addedCellR)*nsprs_InACell 
             endif
             
          enddo
          
          node_area(i,0) = node_areaNAC(node_equiv,0)
          jmaxA          = node_areaNAC(node_equiv,0)
          
          do j = 1,jmaxA
             areaInNAC    = node_areaNAC(node_equiv,j)
             if (areaInNAC .le. Hlf_NcellNAC) node_area(i,j) = node_areaNAC(node_equiv,j)+ (N_addedCellL)
             if (areaInNAC .gt. Hlf_NcellNAC) node_area(i,j) = node_areaNAC(node_equiv,j)+(N_addedCellL+N_addedCellR)   
          enddo
          
       elseif ((i.gt.lim6) .and. (i.le.lim7)) then
          
          if ((i-lim6) == 1) prvSprStrt = -1
          if ((i-lim6).gt.1) prvSprStrt = node_spr((i-1),1)
          
          if (mod(cnt_NI,2) == 1) then
             
             if (mod(cnt_NI,4)==1) currSprStrt = prvSprStrt + 3
             if (mod(cnt_NI,4)==3) currSprStrt = prvSprStrt + 2
             
             if (i==(N_node-1)) then
                
                modVal = mod(cnt_NI,4)
                
                if (modVal.ne.1) then
                   write(*,*) "modVal should be 1",modVal,cnt_NI,N_node
                   stop
                endif
                
                currSprStrt = prvSprStrt + 2 ! exception
                
             endif
             
             cnt_NI = cnt_NI+1
             
          elseif (mod(cnt_NI,2) == 0) then
             currSprStrt = prvSprStrt + 1
             cnt_NI      = cnt_NI+1
          endif
          
          node_spr(i,0) = 2 
          node_spr(i,1) = currSprStrt
          node_spr(i,2) = currSprStrt+1
          
          if (mod(cnt_NIA,2) == 1) then
             
             if (mod(cnt_NIA,4)==1) then
                area_cnt = area_cnt+1
                
                node_area(i,0) = 1
                node_area(i,1) = area_cnt
                
             elseif (mod(cnt_NIA,4)==3) then
                node_area(i,0) = 2
                node_area(i,1) = area_cnt
                node_area(i,2) = area_cnt+1
             endif
             
             cnt_NIA = cnt_NIA+1
             
          elseif (mod(cnt_NIA,2) == 0) then
             
             node_area(i,0:max_area_node) = node_area((i-1),0:max_area_node)
             cnt_NIA = cnt_NIA+1
          endif
          
       endif
       
       if (i==last_INLftAC) then
          node_area((last_INLftAC-1),2) = N_cell
          node_area(last_INLftAC,2)     = N_cell
       endif
       
    enddo
    
    
    open(unit=89,file='node_sprWOtxtAC_General.dat')
    open(unit=87,file='node_areaWOtxtAC_General.dat')
    
    do i = 1,N_node
       write(89,*) i,node_spr(i,0:max_spr_node)
       write(87,*) i,node_area(i,0:max_area_node)
    enddo
    
    close(89)
    close(87)
    
  end subroutine update_node_other_AddedCell
  
  subroutine update_spr_other_AddedCell
    implicit none
    integer :: lastLTSprLbfrMeet,lastLTSprRbfrMeet
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7
    integer :: i,j,jmax,jStr,jEnd,jmaxN,jmaxA
    integer :: node_cnctdAct,node_nm,area_nm
    integer :: spr_equiv,node_equiv,area_equiv
    
    do i = 1,2
       
       if (i==1) then
          jStr = 1*(nsprs_InACell)                ; jEnd = (Hlf_NcellNAC)*(nsprs_InACell)
       elseif (i==2) then
          jStr = (Hlf_NcellNAC+1)*(nsprs_InACell) ; jEnd = (N_cellNAC)*(nsprs_InACell)
       endif
       
       write(*,*) jStr,jEnd,"jStr,jEnd"
       
       do j = jStr,jEnd
          
          node_cnctdAct = spr_nodeNAC(j,0)
          !write(*,*) node_cnctdAct,"node_cnctdAct"
          
          if (node_cnctdAct == max_node_spr) then
             
             if (i==1) then
                lastLTSprLbfrMeet = j-1
                exit
             elseif (i==2) then
                LastLTSprRbfrMeet = j-1
                exit
             endif
             
          endif
          
       enddo
    enddo
    
    write(*,*) lastLTSprLbfrMeet,lastLTSprRbfrMeet,"lastLT"
    
    lim1 = lastLTSprLbfrMeet
    lim2 = lim1 + (N_addedCellL)*(nsprs_InACell)
    lim3 = (Hlf_Ncell)*(nsprs_InACell)
    lim4 = lastLTSprRbfrMeet + (N_addedCellL)*(nsprs_InACell)
    lim5 = lim4 + (N_addedCellR)*(nsprs_InACell)
    lim6 = (N_cell-1)*(nsprs_InACell)
    lim7 = lim6 + (NAEC_Ltrl+1)
    
    write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,lim7,"limsSpr"
    
    do i = 1,N_spr
       
       if (i.le.lim1) then
          
          spr_equiv = i
          
          jmaxN         = spr_nodeNAC(spr_equiv,0)
          jmaxA         = spr_areaNAC(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_nodeNAC(spr_equiv,j)
             
             if (node_equiv .le. N_lstTN_NAC) node_nm = node_equiv
             if (node_equiv .gt. N_lstTN_NAC) node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2
             
             spr_node(i,j) = node_nm
             
          enddo

          do j = 1,jmaxA
             area_equiv    = spr_areaNAC(spr_equiv,j)
             area_nm       = area_equiv
             spr_area(i,j) = area_nm
          enddo
          
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          spr_equiv = i - nsprs_InACell
          
          jmaxN         = spr_node(spr_equiv,0)
          jmaxA         = spr_area(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_node(spr_equiv,j)
             
             if (node_equiv .le. N_lstTN_AC) node_nm = node_equiv + 2
             if (node_equiv .gt. N_lstTN_AC) node_nm = node_equiv + (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)

             spr_node(i,j) = node_nm
             
          enddo
          
          do j = 1,jmaxA
             area_equiv    = spr_area(spr_equiv,j)
             area_nm       = area_equiv + 1
             spr_area(i,j) = area_nm
          enddo
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          spr_equiv = i - (nsprs_InACell)*(N_addedCellL)
          
          jmaxN         = spr_nodeNAC(spr_equiv,0)
          jmaxA         = spr_areaNAC(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_nodeNAC(spr_equiv,j)
             
             if (node_equiv .lt. N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL)*2
                
             elseif (node_equiv == N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2
                
             elseif (node_equiv.gt.N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2 &
                     + (N_addedCellL)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             spr_node(i,j) = node_nm
             
          enddo
          
          do j = 1,jmaxA
             area_equiv = spr_areaNAC(spr_equiv,j)

             if (area_equiv .ne. N_cellNAC) area_nm = area_equiv + N_addedCellL
             if (area_equiv == N_cellNAC)   area_nm = area_equiv + N_addedCellL +  N_addedCellR
             
             spr_area(i,j) = area_nm
          enddo
          
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          spr_equiv = i - (N_addedCellL)*(nsprs_InACell)
          
          jmaxN         = spr_nodeNAC(spr_equiv,0)
          jmaxA         = spr_areaNAC(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_nodeNAC(spr_equiv,j)
             
             if (node_equiv .le. N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL)*2
                
             elseif (node_equiv .gt. N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2 &
                     + (N_addedCellL)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             spr_node(i,j) = node_nm
             
          enddo
          
          do j = 1,jmaxA
             area_equiv = spr_areaNAC(spr_equiv,j)
             area_nm    = area_equiv + N_addedCellL
             spr_area(i,j) = area_nm
          enddo
          
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          spr_equiv = i - nsprs_InACell
          
          jmaxN         = spr_node(spr_equiv,0)
          jmaxA         = spr_area(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_node(spr_equiv,j)
             
             if (node_equiv .le. N_lstTN_AC) node_nm = node_equiv + 2
             if (node_equiv .gt. N_lstTN_AC) node_nm = node_equiv + (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             
             spr_node(i,j) = node_nm
             
          enddo
          
          do j = 1,jmaxA
             area_equiv    = spr_area(spr_equiv,j)
             area_nm       = area_equiv + 1
             spr_area(i,j) = area_nm
          enddo
          
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          
          spr_equiv = i - (nsprs_InACell)*(N_addedCellL+N_addedCellR)
          
          jmaxN         = spr_nodeNAC(spr_equiv,0)
          jmaxA         = spr_areaNAC(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_nodeNAC(spr_equiv,j)
             
             if (node_equiv .le. N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2
                
             elseif (node_equiv.gt.N_lstTN_NAC) then
                node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2 &
                     + (N_addedCellL+N_addedCellR)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             spr_node(i,j) = node_nm
             
          enddo
          
          
          do j = 1,jmaxA
             area_equiv = spr_areaNAC(spr_equiv,j)
             area_nm    = area_equiv + N_addedCellL + N_addedCellR

             spr_area(i,j) = area_nm
          enddo
          
          
       elseif ((i.gt.lim6) .and. (i.le.lim7)) then
          
          spr_equiv = i - (nsprs_InACell)*(N_addedCellL+N_addedCellR)

          jmaxN         = spr_nodeNAC(spr_equiv,0)
          jmaxA         = spr_areaNAC(spr_equiv,0)
          
          spr_node(i,0) = jmaxN
          spr_area(i,0) = jmaxA
          
          do j = 1,jmaxN
             
             node_equiv = spr_nodeNAC(spr_equiv,j)

             if (node_equiv .le. N_lstTN_NAC) then
                
                if ((i-lim6) == 1)   node_nm = node_equiv + (N_addedCellL*2)
                if ((i-lim6) .gt. 1) node_nm = node_equiv + (N_addedCellL+N_addedCellR)*2
                
             elseif (node_equiv .gt. N_lstTN_NAC) then
                node_nm = node_equiv +  + (N_addedCellL+N_addedCellR)*2 &
                     + (N_addedCellL+N_addedCellR)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             spr_node(i,j) = node_nm
             
          enddo
          
          do j = 1,jmaxA
             area_equiv = spr_areaNAC(spr_equiv,j)
             area_nm    = area_equiv + N_addedCellL + N_addedCellR
             spr_area(i,j) = area_nm
          enddo
          
       endif
       
    enddo
    
    open(unit=89,file='spr_nodeWOtxtAC_General.dat')
    open(unit=87,file='spr_areaWOtxtAC_General.dat')
    
    do i = 1,N_spr
       write(89,*) i,spr_node(i,0:max_node_spr)
       write(87,*) i,spr_area(i,0:max_area_spr)
    enddo
    
    close(89)
    close(87)
    
  end subroutine update_spr_other_AddedCell
  
  
  subroutine update_area_other_AddedCell
    implicit none
    integer :: lastAreaLbfrMeet,lastAreaRbfrMeet
    integer :: i,j,jStr,jEnd,jmax,jmaxN,Hlf_jmaxN,jmaxS
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7
    integer :: node_nm,spr_nm
    integer :: node_nmPrv,node_nmEquiv,node_nmCur
    integer :: prvArea,area_Equiv
    integer :: node_cnctdAct
    integer :: prvAreaSpr,curAreaSpr
    integer :: cur_spr,spr_equiv 
    
    do i = 1,2
       
       if (i==1) then
          jStr = 1 ; jEnd = Hlf_NcellNAC
       elseif (i==2) then
          jStr = (Hlf_NcellNAC+1) ; jEnd = N_cellNAC
       endif
       
       !write(*,*) jStr,jEnd,"jStr,jEnd"
       
       do j = jStr,jEnd
          
          node_cnctdAct = area_nodeNAC(j,0)
          !write(*,*) node_cnctdAct,"node_cnctdAct"
          
          if (node_cnctdAct == max_node_area) then
             if (i==1) then
                lastAreaLbfrMeet = j-1
                exit
             elseif (i==2) then
                lastAreaRbfrMeet = j-1
                exit
             endif
          endif
          
       enddo
       
    enddo
    
    !write(*,*) lastAreaLbfrMeet,lastAreaRbfrMeet,"lastArea(L/R)"
    
    lim1 = lastAreaLbfrMeet
    lim2 = lim1 + N_addedCellL
    lim3 = Hlf_Ncell
    lim4 = lastAreaRbfrMeet + N_addedCellL
    lim5 = lastAreaRbfrMeet + N_addedCellL + N_addedCellR
    lim6 = N_cell-1
    lim7 = N_cell
    
    !write(*,*) lim1,lim2,lim3,lim4,lim5,lim6,lim7,"lims Area"
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          
          jmaxN          = area_nodeNAC(i,0)
          area_node(i,0) = area_nodeNAC(i,0)
          
          do j = 1,jmaxN
             node_nm = area_nodeNAC(i,j)
             
             if (node_nm .le. N_lstTN_NAC) area_node(i,j) = node_nm
             if (node_nm .gt. N_lstTN_NAC) area_node(i,j) = node_nm + (N_addedCellL+N_addedCellR)*2
          enddo
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          prvArea = i-1
          
          jmaxN          = area_node(prvArea,0)
          area_node(i,0) = area_node(prvArea,0)
          
          do j = 1,jmaxN
             node_nmPrv = area_node(prvArea,j)
             !write(*,*) node_nmPrv,prvArea,"node_nmPrv"
             
             if (node_nmPrv .le. N_lstTN_AC) then
                node_nmCur = node_nmPrv + 2
             elseif (node_nmPrv .gt. N_lstTN_AC) then
                node_nmCur = node_nmPrv + (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             area_node(i,j) = node_nmCur
             
          enddo
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          area_equiv = i-N_addedCellL
          
          jmaxN          = area_nodeNAC(area_equiv,0)
          area_node(i,0) = area_nodeNAC(area_equiv,0)
          
          do j = 1,jmaxN
             
             node_nmEquiv = area_nodeNAC(area_equiv,j)
             
             if (node_nmEquiv .le. N_lstTN_NAC) then
                
                if (node_nmEquiv .ne. N_lstTN_NAC) then
                   area_node(i,j) = (node_nmEquiv) + (N_addedCellL*2)
                elseif (node_nmEquiv == N_lstTN_NAC) then
                   area_node(i,j) = (node_nmEquiv) + (N_addedCellL+N_addedCellR)*2
                endif
                
             elseif (node_nmEquiv .gt. N_lstTN_NAC) then
                area_node(i,j) = (node_nmEquiv) + (N_addedCellL+N_addedCellR)*2 + (N_addedCellL)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
          enddo
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          area_equiv = i-N_addedCellL
          
          jmaxN          = area_nodeNAC(area_equiv,0)
          area_node(i,0) = area_nodeNAC(area_equiv,0)
          
          do j = 1,jmaxN
             node_nmEquiv   = area_nodeNAC(area_equiv,j)
             
             if (node_nmEquiv .le. N_lstTN_NAC) then
                area_node(i,j) = (node_nmEquiv) + (N_addedCellL)*2 
                
             elseif (node_nmEquiv .gt. N_lstTN_NAC) then
                area_node(i,j) = (node_nmEquiv) + (N_addedCellL+N_addedCellR)*2 + (N_addedCellL)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
                
             endif
             
          enddo
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          prvArea = i-1
          
          jmaxN          = area_node(prvArea,0)
          area_node(i,0) = area_node(prvArea,0)
          
          do j = 1,jmaxN
             
             node_nmPrv = area_node(prvArea,j)
             
             if (node_nmPrv .le. N_lstTN_AC) then
                node_nmCur = node_nmPrv + 2
             elseif (node_nmPrv .gt. N_lstTN_AC) then
                node_nmCur = node_nmPrv + (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             area_node(i,j) = node_nmCur
             
          enddo
          
       elseif ((i.gt.lim5) .and. (i.le.lim6)) then
          
          area_equiv = i-(N_addedCellL+N_addedCellR)
          !write(*,*) area_equiv,"AE"
          
          jmaxN          = area_nodeNAC(area_equiv,0)
          area_node(i,0) = area_nodeNAC(area_equiv,0)
          
          do j = 1,jmaxN
             
             node_nmEquiv   = area_nodeNAC(area_equiv,j)
             
             if (node_nmEquiv .le. N_lstTN_NAC) then
                node_nmCur  = (node_nmEquiv) + (N_addedCellL+N_addedCellR)*2
                
             elseif (node_nmEquiv .gt. N_lstTN_NAC) then
                node_nmCur  = (node_nmEquiv) + (N_addedCellL+N_addedCellR)*2 &
                     + (N_addedCellL+N_addedCellR)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
             endif
             
             area_node(i,j) = node_nmCur 
             
          enddo
          
       elseif ((i.gt.lim6) .and. (i.le.lim7)) then
          
          area_equiv = i-(N_addedCellL+N_addedCellR)
          
          jmaxN          = area_nodeNAC(area_equiv,0)
          area_node(i,0) = area_nodeNAC(area_equiv,0)
          
          Hlf_jmaxN = jmaxN/2
          
          !write(*,*) Hlf_jmaxN,"Hlf_jmaxN"
          
          do j = 1,jmaxN
             
             node_nmEquiv = area_nodeNAC(area_equiv,j)
             
             if (j.le.Hlf_jmaxN) then
                
                if (node_nmEquiv .le. N_lstTN_NAC) then
                   node_nmCur = node_nmEquiv + (N_addedCellL)*2
                   
                elseif (node_nmEquiv .gt. N_lstTN_NAC) then
                   node_nmCur = node_nmEquiv + (N_addedCellL+N_addedCellR)*2 + (N_addedCellL)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
                endif
                
             elseif (j.gt.Hlf_jmaxN) then
                
                if (node_nmEquiv .le. N_lstTN_NAC) then
                   node_nmCur = node_nmEquiv + (N_addedCellL+N_addedCellR)*2
                   
                elseif (node_nmEquiv .gt. N_lstTN_NAC) then
                   node_nmCur = node_nmEquiv + (N_addedCellL+N_addedCellR)*2 &
                        + (N_addedCellL+N_addedCellR)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
                   
                endif
                
             endif
             
             area_node(i,j) = node_nmCur
             
          enddo
          
       endif
       
    enddo
    
    lim1 = Hlf_NcellNAC
    lim2 = Hlf_Ncell
    lim3 = lim2+Hlf_NcellNAC
    lim4 = N_cell-1
    lim5 = N_cell
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          area_spr(i,0:max_spr_area) = area_sprNAC(i,0:max_spr_area)
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          prvArea = i-1
          
          jmaxS          = area_spr(prvArea,0)
          area_spr(i,0)  = area_spr(prvArea,0)
          
          do j = 1,jmaxS
             prvAreaSpr = area_spr(prvArea,j)
             curAreaSpr = prvAreaSpr + (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
             
             area_spr(i,j) = curAreaSpr
          enddo
          
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          area_equiv = i-N_addedCellL
          
          jmaxS           = area_sprNAC(area_equiv,0)
          area_spr(i,0)   = area_sprNAC(area_equiv,0)
          
          do j = 1,jmaxS
             area_spr(i,j) = area_sprNAC(area_equiv,j) + (NAEC_Apcl+1+NAEC_Bsal+1+NAEC_Ltrl+1)*(N_addedCellL)
          enddo
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          prvArea = i-1
          
          jmaxS          = area_spr(prvArea,0)
          area_spr(i,0)  = area_spr(prvArea,0)
          
          do j = 1,jmaxS
             prvAreaSpr = area_spr(prvArea,j)
             curAreaSpr = prvAreaSpr + (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
             
             area_spr(i,j) = curAreaSpr
          enddo
          
       elseif ((i.gt.lim4) .and. (i.le.lim5)) then
          
          area_equiv = i-(N_addedCellL+N_addedCellR)
          
          jmaxS           = area_sprNAC(area_equiv,0)
          area_spr(i,0)   = area_sprNAC(area_equiv,0)
            
          do j = 1,jmaxS
             
             spr_equiv = area_sprNAC(area_equiv,j)
             
             !write(*,*) Hlf_NsprNAC,spr_equiv,"spr_equiv"
             
             if (spr_equiv.le.Hlf_NsprNAC) then
                cur_spr = spr_equiv + (NAEC_Apcl+1+NAEC_Bsal+1+NAEC_Ltrl+1)*(N_addedCellL)
             elseif (spr_equiv.gt.Hlf_NsprNAC) then
                cur_spr = spr_equiv + (NAEC_Apcl+1+NAEC_Bsal+1+NAEC_Ltrl+1)*(N_addedCellL+N_addedCellR)
             endif

             !write(*,*) cur_spr,"cur_spr"
             area_spr(i,j) = cur_spr
             
          enddo
          
       endif
       
    enddo
    
    open(unit=89,file='area_nodeWOtxtAC_General.dat')
    open(unit=87,file='area_sprWOtxtAC_General.dat')
    
    do i = 1,N_cell
       write(89,*) i,area_node(i,0:max_node_area)
       write(87,*) i,area_spr(i,0:max_spr_area)
    enddo
    
    close(89)
    close(87)
    
  end subroutine update_area_other_AddedCell
  
  
  subroutine update_node_spr_readVrsn
    implicit none
    integer :: i
    integer :: nodeNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*) "N_cell not equal to AddedVrsn"
    endif
    
    open(unit=221,file='node_sprWOtxtAC.dat')
    open(unit=223,file='node_sprWOtxtACchkData.dat')
    
    do i = 1,N_node
       read(221,*)  nodeNm,node_spr(nodeNm,0:max_spr_node)
       write(223,*) nodeNm,node_spr(nodeNm,0:max_spr_node)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_node_spr_readVrsn
  
  subroutine update_node_area_readVrsn
    implicit none
    integer :: i
    integer :: nodeNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*)"N_cell not equal to AddedVrsn"
       stop
    endif
    
    open(unit=221,file='node_areaWOtxtAC.dat')
    open(unit=223,file='node_areaWOtxtAC.chkData.dat')
    
    write(*,*) N_node,"BFR Genjam"
    write(*,*) max_area_node,"BFR Genjam"
    
    do i = 1,N_node
       read(221,*)  nodeNm,node_area(nodeNm,0:max_area_node)
       write(223,*) nodeNm,node_area(nodeNm,0:max_area_node)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_node_area_readVrsn
  
  subroutine update_spr_node_readVrsn
    implicit none
    integer :: sprNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*) "N_cell not equal to AddedVrsn"
    endif
    
    open(unit=221,file='spr_nodeWOtxtAC.dat')
    open(unit=223,file='spr_nodeWOtxtACchkData.dat')
    
    do i = 1,N_spr
       read(221,*)  sprNm,spr_node(sprNm,0:max_node_spr)
       write(223,*) sprNm,spr_node(sprNm,0:max_node_spr)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_spr_node_readVrsn
  
  subroutine update_spr_area_readVrsn
    implicit none
    integer :: i
    integer :: sprNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*) "N_cell not equal to AddedVrsn"
    endif
    
    open(unit=221,file='spr_areaWOtxtAC.dat')
    open(unit=223,file='spr_areaWOtxtACchkData.dat')
    
    do i = 1,N_spr
       read(221,*)  sprNm,spr_area(sprNm,0:max_area_spr)
       write(223,*) sprNm,spr_area(sprNm,0:max_area_spr)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_spr_area_readVrsn
  
  subroutine update_area_node_readVrsn
    implicit none
    integer :: i
    integer :: areaNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*) "N_cell not equal to AddedVrsn"
    endif
    
    open(unit=221,file='area_nodeWOtxtAC.dat')
    open(unit=223,file='area_nodeWOtxtACchkData.dat')
    
    do i = 1,N_cell
       read(221,*) areaNm,area_node(areaNm,0:max_node_area) 
       read(223,*) areaNm,area_node(areaNm,0:max_node_area)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_area_node_readVrsn
  
  subroutine update_area_spr_readVrsn
    implicit none
    integer :: i
    integer :: areaNm
    
    if (N_cell.ne.NcellAddedVrsn) then
       write(*,*) "N_cell not equal to AddedVrsn"
    endif
    
    open(unit=221,file='area_sprWOtxtAC.dat')
    open(unit=223,file='area_sprWOtxtACchkData.dat')
    
    do i = 1,N_cell
       read(221,*) areaNm,area_spr(areaNm,0:max_spr_area) 
       read(223,*) areaNm,area_spr(areaNm,0:max_spr_area)
    enddo
    
    close(221)
    close(223)
    
  end subroutine update_area_spr_readVrsn

  subroutine update_node_spr_writeVrsn
    implicit none
    
    continue
  end subroutine update_node_spr_writeVrsn
  
  subroutine update_area_node
    implicit none
    integer :: lim1,lim2,lim3,lim4,lim5
    integer :: diffInCellNm,diffInTNs,diffInINs
    integer :: cellNmNAC,cellNmNxtNAC
    integer :: nodeInNxtcellNAC,nodeShdBeInCurrcellNAC,nodeInNAC
    integer :: i,j
    integer :: nodeDiffForTNsEachCell,nodeDiffForINsEachCell
    integer :: nodeInNxtNAC
    
    lim1 = N_addedCellL
    lim2 = lim1 + Hlf_NcellNAC
    lim3 = lim2 + N_addedCellR
    lim4 = lim3 + Hlf_NcellNAC
    lim5 = N_cell
    
    nodeDiffForTNsEachCell = 2
    nodeDiffForINsEachCell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          
          if (i==1) then
             diffInCellNm = i-N_addedCellL-1  ! will take prop from cell 1
             cellNmNAC    = 1
          elseif ((i.gt.1) .and. (i.le.lim1)) then
             diffInCellNm = i-N_addedCellL-2 ! will take prop from cell 2
             cellNmNAC    = 2
          endif
          
          diffInTNs = (diffInCellNm)*(nodeDiffForTNsEachCell)
          diffInINs = (diffInCellNm)*(nodeDiffForINsEachCell)
          
          do j = 0,max_node_area
             
             if (j==0) then
                area_node(i,j) = area_nodeNAC(cellNmNAC,j)
                
             elseif (j.gt.0) then
                nodeInNAC = area_nodeNAC(cellNmNAC,j)
                
                if (nodeInNxtNAC .le. N_lstTN_NAC) then               
                   area_node(i,j) = nodeInNAC - diffInTNs
                   
                elseif (nodeInNxtNAC .gt. N_lstTN_NAC) then
                   area_node(i,j) = nodeInNAC - diffInINs
                endif
                
             endif
             
          enddo
          
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          if ((i-lim1) == 1) then
             
             cellNmNAC    = i-N_addedCellL
             cellNmNxtNAC = cellNmNAC+1 ! For reversing prop from 2 to 1
             
             do j = 0,max_node_area
                
                if (j==0) then
                   area_node(i,j) = area_nodeNAC(cellNmNxtNAC,j)
                   
                elseif (j.gt.0) then
                   nodeInNxtcellNAC = area_nodeNAC(cellNmNxtNAC,j)
                   
                   if (nodeInNxtNAC .le. N_lstTN_NAC) then
                      nodeShdBeInCurrcellNAC = nodeInNxtcellNAC-nodeDiffForTNsEachCell
                      area_node(i,j) = nodeShdBeInCurrcellNAC
                      
                   elseif (nodeInNxtNAC .gt. N_lstTN_NAC) then
                      nodeShdBeInCurrcellNAC = nodeInNxtcellNAC - nodeDiffForINsEachCell
                      area_node(i,j) = nodeShdBeInCurrcellNAC
                   endif
                   
                endif
                
             enddo
             
          else
             cellNmNAC = i-N_addedCellL
             area_node(i,0:max_node_area) = area_nodeNAC(cellNmNAC,0:max_node_area)
          endif
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          if ((i-lim2)==1) then
             diffInCellNm = i-(N_addedCellR)-(Hlf_Ncell+1)  ! will take prop from cell 12
             cellNmNAC    = Hlf_NcellNAC + 1
          elseif (((i-lim2).gt.1) .and. (i.le.lim3)) then
             diffInCellNm = i-(N_addedCellR)-(Hlf_Ncell+2) ! will take prop from cell 13
             cellNmNAC    = Hlf_NcellNAC + 2
          endif
          
          diffInTNs = (diffInCellNm)*(nodeDiffForTNsEachCell)
          diffInINs = (diffInCellNm)*(nodeDiffForINsEachCell)
          
          do j = 0,max_node_area
             
             if (j==0) then
                area_node(i,j) = area_nodeNAC(cellNmNAC,j)
                
             elseif (j.gt.0) then
                nodeInNAC = area_nodeNAC(cellNmNAC,j)
                
                if (nodeInNxtNAC .le. N_lstTN_NAC) then               
                   area_node(i,j) = nodeInNAC - diffInTNs
                   
                elseif (nodeInNxtNAC .gt. N_lstTN_NAC) then
                   area_node(i,j) = nodeInNAC - diffInINs
                endif
                
             endif
             
          enddo
          
          
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          if ((i-lim3) == 1) then
             
             cellNmNAC    = i-(N_addedCellL+N_addedCellR)
             cellNmNxtNAC = cellNmNAC+1 ! For reversing prop from 13 to 12
             
             do j = 0,max_node_area
                
                if (j==0) then
                   area_node(i,j) = area_nodeNAC(cellNmNxtNAC,j)
                   
                elseif (j.gt.0) then
                   nodeInNxtcellNAC = area_nodeNAC(cellNmNxtNAC,j)
                   
                   if (nodeInNxtcellNAC == (-1)) then
                      continue
                   elseif (nodeInNxtcellNAC .ne. (-1)) then
                      
                      if (nodeInNxtcellNAC .le. N_lstTN_NAC) then
                         nodeShdBeInCurrcellNAC = nodeInNxtcellNAC-2
                         area_node(i,j) = nodeShdBeInCurrcellNAC
                         
                      elseif (nodeInNxtcellNAC .gt. N_lstTN_NAC) then
                         nodeShdBeInCurrcellNAC = nodeInNxtcellNAC - (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
                         area_node(i,j) = nodeShdBeInCurrcellNAC      
                      endif
                      
                   endif
                endif
                
             enddo
             
          else
             cellNmNAC = i-(N_addedCellL+N_addedCellR)
             area_node(i,0:max_node_area) = area_nodeNAC(cellNmNAC,0:max_node_area)
          endif
          
       elseif (i.gt.lim4) then
          continue
       endif
       
       
    enddo
    
  end subroutine update_area_node
  
  
  subroutine print_all_NAC_trnsfrms
    implicit none
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 2
    
    open(unit=148,file='NAC_Transform_Variables.dat')
    
    do i = 1,N_itm
       
       if (i==1) write(148,*) "node_spr:"
       if (i==2) write(148,*) "node_area:"
       
       do j = 1,N_nodeNAC
          write(148,*) "node_nm =",j
          
          if (i==1) then
             write(148,*) node_sprNAC(j,0:max_spr_nodeNAC)
          elseif (i==2) then
             write(148,*) node_areaNAC(j,0:max_area_nodeNAC)
          endif
          
       enddo
       
       write(148,*) " "
       
    enddo
    
    
    do i = 1,N_itm
       if (i==1) write(148,fmt=*) "spr_node:"
       if (i==2) write(148,fmt=*) "spr_area:"
       
       do j = 1,N_sprNAC
          write(148,fmt=*) "spr_nm =",j
          
          if (i==1) then
             write(148,fmt=*) spr_nodeNAC(j,0:max_node_sprNAC) 
          elseif (i==2) then
             write(148,fmt=*) spr_areaNAC(j,0:max_area_sprNAC)
          endif
          
       enddo
       write(148,fmt=*) " "
       
    enddo
    
    do i = 1,N_itm
       if (i==1) write(148,fmt=*) "area_node:"
       if (i==2) write(148,fmt=*) "area_spr:"
       
       do j = 1,N_cellNAC
          write(148,fmt=*) "area_nm =",j
          
          if (i==1) then
             write(148,fmt=*) area_nodeNAC(j,0:max_node_areaNAC) 
          elseif (i==2) then
             write(148,fmt=*) area_sprNAC(j,0:max_spr_areaNAC)
          endif
          
       enddo
       write(148,fmt=*) " "
       
    enddo
    
       
    do i = 1,N_nodeNAC
       write(148,fmt=*) "Nlist for  node_nm =",i
       
       jmax = node_areaNAC(i,0)
       
       do j = 1,jmax
          write(148,fmt=*) NlistNAC(i,j,1:2)
       enddo
       
       write(148,fmt=*) " "
    enddo
    
    write(148,fmt=*) N_phi,"N_phi"
    write(148,fmt=*) max_node_sprNAC,max_area_sprNAC,"max_node_spr,max_area_spr"
    write(148,fmt=*) max_spr_nodeNAC,max_area_nodeNAC,"max_spr_node,max_area_node"
    write(148,fmt=*) max_node_areaNAC,max_spr_areaNAC,"max_node_area,max_spr_area"
    
    write(148,fmt=*) N_node,N_nodeNAC,"N_node,N_nodeNAC"
    write(148,fmt=*) N_spr,N_sprNAC,"N_spr,N_sprNAC"
    
    close(148)
    
  end subroutine print_all_NAC_trnsfrms

  
  ! subroutine print_all_NAC_trnsfrms
  !   implicit none
  !   integer :: i,j,jmax
  !   integer :: N_itm
    
  !   N_itm = 2
    
  !   open(unit=148,file='NAC_Transform_Variables.dat')
    
  !   do i = 1,N_itm
       
  !      if (i==1) write(148,*) "node_spr:"
  !      if (i==2) write(148,*) "node_area:"
       
  !      do j = 1,N_nodeNAC
  !         write(148,*) "node_nm =",j
          
  !         if (i==1) then
  !            write(148,*) node_sprNAC(j,0:max_spr_nodeNAC)
  !         elseif (i==2) then
  !            write(148,*) node_areaNAC(j,0:max_area_nodeNAC)
  !         endif
          
  !      enddo
       
  !      write(148,*) " "
       
  !   enddo
    
    
  !   do i = 1,N_itm
  !      if (i==1) write(148,fmt=*) "spr_node:"
  !      if (i==2) write(148,fmt=*) "spr_area:"
       
  !      do j = 1,N_sprNAC
  !         write(148,fmt=*) "spr_nm =",j
          
  !         if (i==1) then
  !            write(148,fmt=*) spr_nodeNAC(j,0:max_node_sprNAC) 
  !         elseif (i==2) then
  !            write(148,fmt=*) spr_areaNAC(j,0:max_area_sprNAC)
  !         endif
          
  !      enddo
  !      write(148,fmt=*) " "
       
  !   enddo
    
  !   do i = 1,N_itm
  !      if (i==1) write(148,fmt=*) "area_node:"
  !      if (i==2) write(148,fmt=*) "area_spr:"
       
  !      do j = 1,N_cellNAC
  !         write(148,fmt=*) "area_nm =",j
          
  !         if (i==1) then
  !            write(148,fmt=*) area_nodeNAC(j,0:max_node_areaNAC) 
  !         elseif (i==2) then
  !            write(148,fmt=*) area_sprNAC(j,0:max_spr_areaNAC)
  !         endif
          
  !      enddo
  !      write(148,fmt=*) " "
       
  !   enddo
    
       
  !   do i = 1,N_nodeNAC
  !      write(148,fmt=*) "Nlist for  node_nm =",i
       
  !      jmax = node_areaNAC(i,0)
       
  !      do j = 1,jmax
  !         write(148,fmt=*) NlistNAC(i,j,1:2)
  !      enddo
       
  !      write(148,fmt=*) " "
  !   enddo
    
  !   write(148,fmt=*) N_phi,"N_phi"
  !   write(148,fmt=*) max_node_sprNAC,max_area_sprNAC,"max_node_spr,max_area_spr"
  !   write(148,fmt=*) max_spr_nodeNAC,max_area_nodeNAC,"max_spr_node,max_area_node"
  !   write(148,fmt=*) max_node_areaNAC,max_spr_areaNAC,"max_node_area,max_spr_area"
    
  !   write(148,fmt=*) N_node,N_nodeNAC,"N_node,N_nodeNAC"
  !   write(148,fmt=*) N_spr,N_sprNAC,"N_spr,N_sprNAC"
    
  !   close(148)
    
  ! end subroutine print_all_NAC_trnsfrms
  
  
  subroutine print_all_AC_trnsfrms_automated
    implicit none
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 2
    
    open(unit=148,file='AC_Transform_Variables.dat')
    
    do i = 1,N_itm
       
       if (i==1) write(148,*) "node_spr:"
       if (i==2) write(148,*) "node_area:"
       
       do j = 1,N_node
          !write(148,*) "node_nm =",j
          
          if (i==1) then
             write(148,*) j,node_spr(j,0:max_spr_node)
          elseif (i==2) then
             write(148,*) j,node_area(j,0:max_area_node)
          endif
          
       enddo
       
       write(148,*) " "
       
    enddo
    
    
    do i = 1,N_itm
       if (i==1) write(148,*) "spr_node:"
       if (i==2) write(148,*) "spr_area:"
       
       do j = 1,N_spr
          !write(148,*) "spr_nm =",j
          
          if (i==1) then
             write(148,*) j,spr_node(j,0:max_node_spr) 
          elseif (i==2) then
             write(148,*) j,spr_area(j,0:max_area_spr)
          endif
          
       enddo
       write(148,*) " "
       
    enddo
    
    do i = 1,N_itm
       if (i==1) write(148,*) "area_node:"
       if (i==2) write(148,*) "area_spr:"
       
       do j = 1,N_cell
          !write(148,*) "area_nm =",j
          
          if (i==1) then
             write(148,*) j,area_node(j,0:max_node_area) 
          elseif (i==2) then
             write(148,*) j,area_spr(j,0:max_spr_area)
          endif
          
       enddo
       write(148,*) " "
       
    enddo
    
       
    do i = 1,N_node
       write(148,*) "Nlist for  node_nm =",i
       
       jmax = node_area(i,0)
       
       do j = 1,jmax
          write(148,*) Nlist(i,j,1:2)
       enddo
       
       write(148,*) " "
    enddo
    
    write(148,*) N_phi,"N_phi"
    write(148,*) max_node_spr,max_area_spr,"max_node_spr,max_area_spr"
    write(148,*) max_spr_node,max_area_node,"max_spr_node,max_area_node"
    write(148,*) max_node_area,max_spr_area,"max_node_area,max_spr_area"
    
    write(148,*) N_node,N_node,"N_node,N_node"
    write(148,*) N_spr,N_spr,"N_spr,N_spr"
    
    close(148)
    
  end subroutine print_all_AC_trnsfrms_automated
  
  subroutine updt_A0l0_usingTrnsfrms
    implicit none
    integer :: lim1,lim2,lim3,lim4
    integer :: Nnode_poly,Nnode_Spr
    integer :: i,j,jmax,m,m1
    integer :: node_nm,spr_nm,area_nm
    
    real*8, allocatable :: Nodes(:,:),Nodes_Spr(:,:)
    
    interface
       
       real*8 function Polygon_Area(Nodes,Nnode_poly)
         use Triangle_routines
         implicit none
         integer, intent(in) :: Nnode_poly
         real*8 , intent(in) :: Nodes(1:Nnode_poly,1:3)
       end function Polygon_Area
       
       real*8 function LengthForArbtrySgmnt(NodesVal,Nnode_SprVal)
         implicit none
         integer, intent(in) :: Nnode_SprVal
         real*8 , intent(in) :: NodesVal(1:Nnode_SprVal,1:3)
       end function LengthForArbtrySgmnt
       
    end interface
    
    
    lim1 = N_addedCellL
    lim2 = Hlf_Ncell
    lim3 = Hlf_Ncell + N_addedCellR
    lim4 = N_cell
    
    write(*,*) lim1,lim2,lim3,lim4,"lims in updt_A0l0_usingTrnsfrms"
    
    do i = 1,N_cell
       
       if (i.le.lim1) then
          
          Nnode_poly = area_node(i,0)
          allocate(Nodes(1:Nnode_Poly,1:3))
          
          do j = 1,Nnode_poly
             node_nm = area_node(i,j)
             Nodes(j,1:N_dmnsn) = node_xy(node_nm,1:N_dmnsn)
          enddo
          
          A(i)      = Polygon_Area(Nodes,Nnode_poly)
          A0(i)     = A(i)
          k_area(i) = (1.0d0)/(A0(i))
          
          deallocate(Nodes)
          
          do m = 1,nsprs_InACell
             spr_nm = (i-1)*nsprs_InACell + m
             
             
             Nnode_Spr = spr_node(spr_nm,0)
             allocate(Nodes_Spr(1:Nnode_Spr,1:3))
             
             do m1 = 1,Nnode_Spr
                node_nm = spr_node(spr_nm,m1)
                Nodes_Spr(m1,1:N_dmnsn) = node_xy(node_nm,1:N_dmnsn)
                if (N_dmnsn==2) Nodes_Spr(m1,3) = 0.0d0
             enddo
             
             l(spr_nm)     = LengthForArbtrySgmnt(Nodes_Spr,Nnode_Spr)
             l0(spr_nm)    = l(spr_nm)
             k_spr(spr_nm) = (1.0d0)/(l0(spr_nm))
             
             deallocate(Nodes_Spr)
             
          enddo
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          continue
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          Nnode_poly = area_node(i,0)
          allocate(Nodes(1:Nnode_Poly,1:3))
          
          do j = 1,Nnode_poly
             node_nm = area_node(i,j)
             Nodes(j,1:N_dmnsn) = node_xy(node_nm,1:N_dmnsn)
          enddo
          
          A(i)      = Polygon_Area(Nodes,Nnode_poly)
          A0(i)     = A(i)
          k_area(i) = (1.0d0)/(A0(i))
          
          deallocate(Nodes)
          
          do m = 1,nsprs_InACell
             spr_nm = (i-1)*nsprs_InACell + m
             
             
             Nnode_Spr = spr_node(spr_nm,0)
             allocate(Nodes_Spr(1:Nnode_Spr,1:3))
             
             do m1 = 1,Nnode_Spr
                node_nm = spr_node(spr_nm,m1)
                Nodes_Spr(m1,1:N_dmnsn) = node_xy(node_nm,1:N_dmnsn)
                if (N_dmnsn==2) Nodes_Spr(m1,3) = 0.0d0
             enddo
             
             l(spr_nm)     = LengthForArbtrySgmnt(Nodes_Spr,Nnode_Spr)
             l0(spr_nm)    = l(spr_nm)
             k_spr(spr_nm) = (1.0d0)/(l0(spr_nm))
             
             deallocate(Nodes_Spr)
             
          enddo
          
       elseif ((i.gt.lim3) .and.(i.le.lim4)) then
          continue
       endif
       
    enddo

    open(unit=56,file='propAftUpdtngA0l0.dat')
    
    do i = 1,3
       
       if (i==1) then
          jmax = N_cell
          write(56,*) "Cell Props ::"
       elseif (i==2) then
          jmax = N_spr
          write(56,*) "Spr Props ::"
       elseif (i==3) then
          jmax = N_node
          write(56,*) "Node Props ::" 
       endif

       
       do j = 1,jmax
          
          if (i==1) then
             write(56,*) j,k_area(j),A(j),A0(j)
             
          elseif (i==2) then
             write(56,*) j,k_spr(j),l(j),l0(j),typ_spr(j)
             
          elseif (i==3) then
             write(56,*) j,CgXNode(j),nodePhi_typ(j)
             
          endif
          
       enddo
       
       write(56,*) " "
    enddo
    
    close(56)
    
  end subroutine updt_A0l0_usingTrnsfrms

  
  subroutine store_gradient_system_variablesNAC
    implicit none
    
    gradientNAC     = gradient
    gradient_anaNAC = gradient_ana
    gradient_numNAC = gradient_num
    grd_mvNAC       = grd_mv
    grdmv_xyNAC     = grdmv_xy
    
  end subroutine store_gradient_system_variablesNAC
  
  
  subroutine PropValueReadForAddedCell
    implicit none
    real*8 , allocatable :: ksV(:),l0V(:),kaV(:),A0V(:),CgXV(:)
    integer              :: CS1,CS2,CS3,CS4
    
    open(unit=32,file='PropValueReadForAddedCell.dat')
    
    allocate(ksCS1(1:N_sprNAC),  ksCS2(1:N_sprNAC),  ksCS3(1:N_sprNAC),  ksCS4(1:N_sprNAC))
    allocate(l0CS1(1:N_sprNAC),  l0CS2(1:N_sprNAC),  l0CS3(1:N_sprNAC),  l0CS4(1:N_sprNAC))
    allocate(kaCS1(1:N_cellNAC), kaCS2(1:N_cellNAC), kaCS3(1:N_cellNAC), kaCS4(1:N_cellNAC))
    allocate(A0CS1(1:N_cellNAC), A0CS2(1:N_cellNAC), A0CS3(1:N_cellNAC), A0CS4(1:N_cellNAC))
    allocate(CgXCS1(1:N_nodeNAC),CgXCS2(1:N_nodeNAC),CgXCS3(1:N_nodeNAC),CgXCS4(1:N_nodeNAC))
    
    allocate(ksV(1:N_sprNAC),l0V(1:N_sprNAC),kaV(1:N_cellNAC),A0V(1:N_cellNAC),CgXV(1:N_nodeNAC))
    
    CS1=13 ; CS2=14 ; CS3=15 ; CS4=16
    
    call readStrctPropsForContinuing_MoreAddedCell(CS1,ksV,l0V,kaV,A0V,CgXV)
    ksCS1(1:N_sprNAC)   = ksV(1:N_sprNAC)   ; l0CS1(1:N_sprNAC)  = l0V(1:N_sprNAC)
    kaCS1(1:N_cellNAC)  = kaV(1:N_cellNAC)  ; A0CS1(1:N_cellNAc) = A0V(1:N_cellNAC)
    CgXCS1(1:N_nodeNAC) = CgXV(1:N_nodeNAC)
    
    call readStrctPropsForContinuing_MoreAddedCell(CS2,ksV,l0V,kaV,A0V,CgXV)
    ksCS2(1:N_sprNAC)   = ksV(1:N_sprNAC)   ; l0CS2(1:N_sprNAC)  = l0V(1:N_sprNAC)
    kaCS2(1:N_cellNAC)  = kaV(1:N_cellNAC)  ; A0CS2(1:N_cellNAc) = A0V(1:N_cellNAC)
    CgXCS2(1:N_nodeNAC) = CgXV(1:N_nodeNAC)
    
    call readStrctPropsForContinuing_MoreAddedCell(CS3,ksV,l0V,kaV,A0V,CgXV)
    ksCS3(1:N_sprNAC)   = ksV(1:N_sprNAC)   ; l0CS3(1:N_sprNAC)  = l0V(1:N_sprNAC)
    kaCS3(1:N_cellNAC)  = kaV(1:N_cellNAC)  ; A0CS3(1:N_cellNAc) = A0V(1:N_cellNAC)
    CgXCS3(1:N_nodeNAC) = CgXV(1:N_nodeNAC)
    
    call readStrctPropsForContinuing_MoreAddedCell(CS4,ksV,l0V,kaV,A0V,CgXV)
    ksCS4(1:N_sprNAC)   = ksV(1:N_sprNAC)   ; l0CS4(1:N_sprNAC)  = l0V(1:N_sprNAC)
    kaCS4(1:N_cellNAC)  = kaV(1:N_cellNAC)  ; A0CS4(1:N_cellNAc) = A0V(1:N_cellNAC)
    CgXCS4(1:N_nodeNAC) = CgXV(1:N_nodeNAC)
    
    close(32)
    
  end subroutine PropValueReadForAddedCell
  
  subroutine get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
    implicit none
    integer, intent(in)  :: N_lims
    integer, intent(in)  :: limVal(1:N_lims)
    
    real*8 , intent(out) :: ka_inpt(1:N_cell),ka_otpt(1:N_cell)
    real*8 , intent(out) :: A0_inpt(1:N_cell),A0_otpt(1:N_cell)
    real*8 , intent(out) :: ks_inpt(1:N_spr) ,ks_otpt(1:N_spr)
    real*8 , intent(out) :: l0_inpt(1:N_spr) ,l0_otpt(1:N_spr)
    
    integer :: cellCntL,cellCntR
    integer :: sprNmCurr,sprNmNAC
    integer :: i,j,jmax
    
    open(unit=56,file='get_InsOutsForAC.dat',position='append')
    
    cellCntL = 1 ; cellCntR = 1 + Hlf_NcellNAC
    
    do i = 1,N_cell
       
       if (i.le.limVal(1)) then

          if (limVal(1)==0) then
             continue
             
          elseif (limVal(1).ne.0) then
             
             ka_inpt(i) = k_area(i) ; ka_otpt(i) = k_area(i) 
             A0_inpt(i) = A0(i)     ; A0_otpt(i) = A0(i)
             
             write(56,*) i,limVal(1),"cellNm-lim1"
             write(56,*) " "
             
             do j = 1,nsprs_InACell
                
                sprNmCurr = (i-1)*nsprs_InACell + j
                write(56,*) sprNmCurr,"sprNmCurr"
                
                ks_inpt(sprNmCurr) = k_spr(sprNmCurr) ; ks_otpt(sprNmCurr) = k_spr(sprNmCurr) 
                l0_inpt(sprNmCurr) = l0(sprNmCurr)    ; l0_otpt(sprNmCurr) = l0(sprNmCurr)
             
             enddo
             
             write(56,*) " "

          endif
             
       elseif ((i.gt.limVal(1)) .and. (i.le.limVal(2))) then
          
          if (hlfCycl==1) then
             ka_inpt(i) = kaCS1(cellCntL) ; ka_otpt(i) = kaCS2(cellCntL)
             A0_inpt(i) = A0CS1(cellCntL) ; A0_otpt(i) = A0CS2(cellCntL)
             
          elseif (hlfCycl==2) then
             ka_inpt(i) = kaCS3(cellCntL) ; ka_otpt(i) = kaCS4(cellCntL)
             A0_inpt(i) = A0CS3(cellCntL) ; A0_otpt(i) = A0CS4(cellCntL)
          endif
          
          write(56,*) i,limVal(1:2),"cellNm-lim1,2"
          write(56,*) " "
          
          do j = 1,nsprs_InACell
             
             sprNmCurr = (i-1)*nsprs_InACell + j
             sprNmNAC  = (cellCntL-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,sprNmNAC,"sprNmCurrNAC"
             
             if (hlfCycl==1) then
                ks_inpt(sprNmCurr) = ksCS1(sprNmNAC) ; ks_otpt(sprNmCurr) = ksCS2(sprNmNAC)
                l0_inpt(sprNmCurr) = l0CS1(sprNmNAC) ; l0_otpt(sprNmCurr) = l0CS2(sprNmNAC)
                
             elseif (hlfCycl==2) then
                ks_inpt(sprNmCurr) = ksCS3(sprNmNAC) ; ks_otpt(sprNmCurr) = ksCS4(sprNmNAC)
                l0_inpt(sprNmCurr) = l0CS3(sprNmNAC) ; l0_otpt(sprNmCurr) = l0CS4(sprNmNAC)
             endif
             
          enddo
          
          cellCntL = cellCntL+1
          write(56,*) " "
          
       elseif ((i.gt.limVal(2)) .and. (i.le.limVal(3))) then
          
          ka_inpt(i) = k_area(i) ; ka_otpt(i) = k_area(i)
          A0_inpt(i) = A0(i)     ; A0_otpt(i) = A0(i)
          
          write(56,*) i,limVal(2:3),"cellNm-lim2,3"
          write(56,*) " "
          
          do j = 1,nsprs_InACell
             sprNmCurr = (i-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,"sprNmCurr"
             
             ks_inpt(sprNmCurr) = k_spr(sprNmCurr) ; ks_otpt(sprNmCurr) = k_spr(sprNmCurr) 
             l0_inpt(sprNmCurr) = l0(sprNmCurr)    ; l0_otpt(sprNmCurr) = l0(sprNmCurr)
             
          enddo
          
          write(56,*) " "
          
       elseif ((i.gt.limVal(3)) .and. (i.le.limVal(4))) then
          
          ka_inpt(i) = k_area(i) ; ka_otpt(i) = k_area(i)
          A0_inpt(i) = A0(i)     ; A0_otpt(i) = A0(i)
          
          write(56,*) i,limVal(3:4),"cellNm-lim3,4"
          write(56,*) " "
          
          do j = 1,nsprs_InACell
             sprNmCurr = (i-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,"sprNmCurr"
             
             ks_inpt(sprNmCurr) = k_spr(sprNmCurr) ; ks_otpt(sprNmCurr) = k_spr(sprNmCurr)
             l0_inpt(sprNmCurr) = l0(sprNmCurr)    ; l0_otpt(sprNmCurr) = l0(sprNmCurr)
             
          enddo
          
          write(56,*) " "
          
       elseif ((i.gt.limVal(4)) .and. (i.le.limVal(5))) then

          if (hlfCycl==1) then
             ka_inpt(i) = kaCS1(cellCntR) ; ka_otpt(i) = kaCS2(cellCntR)
             A0_inpt(i) = A0CS1(cellCntR) ; A0_otpt(i) = A0CS2(cellCntR)
             
          elseif (hlfCycl==2) then
             ka_inpt(i) = kaCS3(cellCntR) ; ka_otpt(i) = kaCS4(cellCntR)
             A0_inpt(i) = A0CS3(cellCntR) ; A0_otpt(i) = A0CS4(cellCntR)
          endif
          
          write(56,*) i,limVal(4:5),"cellNm-lim4,5"
          write(56,*) " "
           
          do j = 1,nsprs_InACell
             
             sprNmCurr = (i-1)*nsprs_InACell + j
             sprNmNAC  = (cellCntR-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,sprNmNAC,"sprNmCurrNAC"
             
             if (hlfCycl==1) then
                ks_inpt(sprNmCurr) = ksCS1(sprNmNAC) ; ks_otpt(sprNmCurr) = ksCS2(sprNmNAC)
                l0_inpt(sprNmCurr) = l0CS1(sprNmNAC) ; l0_otpt(sprNmCurr) = l0CS2(sprNmNAC)
                
             elseif (hlfCycl==2) then
                ks_inpt(sprNmCurr) = ksCS3(sprNmNAC) ; ks_otpt(sprNmCurr) = ksCS4(sprNmNAC)
                l0_inpt(sprNmCurr) = l0CS3(sprNmNAC) ; l0_otpt(sprNmCurr) = l0CS4(sprNmNAC)
             endif
                
          enddo
          
          cellCntR = cellCntR+1
          
          write(56,*) " "
          
       elseif ((i.gt.limVal(5)) .and. (i.le.limVal(6))) then
          
          ka_inpt(i) = k_area(i) ; ka_otpt(i) = k_area(i)
          A0_inpt(i) = A0(i)     ; A0_otpt(i) = A0(i)
          write(56,*) i,limVal(5:6),"cellNm-lim5,6"
          write(56,*) " "
          
          do j = 1,nsprs_InACell
             sprNmCurr = (i-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,"sprNmCurr"
             
             ks_inpt(sprNmCurr) = k_spr(sprNmCurr) ; ks_otpt(sprNmCurr) = k_spr(sprNmCurr)
             l0_inpt(sprNmCurr) = l0(sprNmCurr)    ; l0_otpt(sprNmCurr) = l0(sprNmCurr)
          enddo
          
          write(56,*) " "
          
       elseif ((i.gt.limVal(6)) .and. (i.le.limVal(7))) then
          
          ka_inpt(i) = k_area(i) ; ka_otpt(i) = k_area(i)
          A0_inpt(i) = A0(i)     ; A0_otpt(i) = A0(i)
          write(56,*) i,limVal(6:7),"cellNm-lim6,7"
          write(56,*) " "
          
          do j = 1,nsprs_InEndCell
             sprNmCurr = (i-1)*nsprs_InACell + j
             write(56,*) sprNmCurr,"sprNmCurr"
             
             ks_inpt(sprNmCurr) = k_spr(sprNmCurr) ; ks_otpt(sprNmCurr) = k_spr(sprNmCurr)
             l0_inpt(sprNmCurr) = l0(sprNmCurr)    ; l0_otpt(sprNmCurr) = l0(sprNmCurr)
          enddo
          
       endif
       
    enddo
    
    do i = 1,2
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          if (i==1) write(56,*) ks_inpt(j),ks_otpt(j),l0_inpt(j),l0_otpt(j),j,"spr"
          if (i==2) write(56,*) ka_inpt(j),ka_otpt(j),A0_inpt(j),A0_otpt(j),j,"cell"
       enddo
       
       write(56,*) " "
       
    enddo
    
    
    close(56) 
    
    
  end subroutine get_InOutsForAC

  
  
  subroutine get_InOutsForAC_Cycl(cyclNm,hlfcyclNm,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
    implicit none
    integer, intent(in)  :: cyclNm,hlfcyclNm
    real*8 , intent(out) :: ka_inpt(1:N_cell),ka_otpt(1:N_cell)
    real*8 , intent(out) :: A0_inpt(1:N_cell),A0_otpt(1:N_cell)
    real*8 , intent(out) :: ks_inpt(1:N_spr) ,ks_otpt(1:N_spr)
    real*8 , intent(out) :: l0_inpt(1:N_spr) ,l0_otpt(1:N_spr)

    integer :: i,j,jmax
    integer :: csIn,csOt
    
    if (hlfCyclNm==1) then
       csIn = 1 ; csOt = 2
    elseif (hlfCyclNm==2) then
       csIn = 3 ; csOt = 4
    endif
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             ka_inpt(j) = ka_cycl(cyclNm,csIn,j) ; ka_otpt(j) = ka_cycl(cyclNm,csOt,j)
             A0_inpt(j) = A0_cycl(cyclNm,csIn,j) ; A0_otpt(j) = A0_cycl(cyclNm,csOt,j)
             
          elseif (i==2) then
             ks_inpt(j) = ks_cycl(cyclNm,csIn,j) ; ks_otpt(j) = ks_cycl(cyclNm,csOt,j)
             l0_inpt(j) = l0_cycl(cyclNm,csIn,j) ; l0_otpt(j) = l0_cycl(cyclNm,csOt,j)
          endif
          
       enddo
       
    enddo
    
  end subroutine get_InOutsForAC_Cycl
  
  
  
  !subroutine get_CellNosAndSprNosForAddedCycle(cyclNm,cellNos,SprNos)
    !implicit none
    !integer, intent(in)  :: cyclNm
    !integer, intent(out) :: cellNos(1:NcellPropToBeRead)
    !integer, intent(out) :: sprNos(1:NsprPropToBeRead)
    
    !integer :: i,j,jmax
    !integer :: skipCell,skipSpr
    
    !continue
    
  !end subroutine get_CellNosAndSprNosForAddedCycle
  
  subroutine make_AddedCellIdenticaltoFirstNonAddedCell
    implicit none
    integer :: lim1,lim2,lim3,lim4
    integer :: sprNmInFirstNonAddedCellL,sprNmInFirstNonAddedCellR
    integer :: sprNmInCurrCell
    integer :: i,j,jmax
    
    lim1 = FirstNonAddedCellL-1
    lim2 = Hlf_Ncell
    lim3 = FirstNonAddedCellR
    lim4 = N_cell
    
    write(*,*) lim1,lim2,lim3,lim4,"lims in make_AddedCellIdentical"
    
    do i = 1,N_cell
       
       if (i.le.(lim1)) then
          
          k_area(i) = k_area(FirstNonAddedCellL)
          A0(i)     = A0(FirstNonAddedCellL)
          
          do j = 1,nsprs_InACell
             sprNmInFirstNonAddedCellL = (nsprs_InACell)*(FirstNonAddedCellL-1) + j
             sprNmInCurrCell           = (nsprs_InACell)*(i-1) + j
             
             k_spr(sprNmInCurrCell) = k_spr(sprNmInFirstNonAddedCellL)
             l0(sprNmInCurrCell)    = l0(sprNmInFirstNonAddedCellL)
          enddo
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          continue
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          k_area(i) = k_area(FirstNonAddedCellR)
          A0(i)     = A0(FirstNonAddedCellR)
          
          do j = 1,nsprs_InACell
             sprNmInFirstNonAddedCellR = (nsprs_InACell)*(FirstNonAddedCellR-1) + j
             sprNmInCurrCell           = (nsprs_InACell)*(i-1) + j
             
             k_spr(sprNmInCurrCell) = k_spr(sprNmInFirstNonAddedCellR)
             l0(sprNmInCurrCell)    = l0(sprNmInFirstNonAddedCellR)
          enddo
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          continue
       endif
       
    enddo
    
    call Equilibrate_AddedCell_model
    
    open(unit=89,file='make_AddedCellIdenticaltoFirstNonAddedCell.dat')

    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          if (i==1) write(89,*) k_area(i),A0(i),i
          if (i==2) write(89,*) k_spr(i),l0(i),i
       enddo
       
       write(89,*) " "
    enddo
    
    close(89)
    
  end subroutine make_AddedCellIdenticaltoFirstNonAddedCell
  
  subroutine readStrctPropsForContinuing_MoreAddedCell(strctNum,ksV,l0V,kaV,A0V,CgXV)
    implicit none
    integer, intent(in)  :: strctNum
    real*8 , intent(out) :: ksV(1:N_sprNAC) ,l0V(1:N_sprNAC)
    real*8 , intent(out) :: kaV(1:N_cellNAC),A0V(1:N_cellNAC)
    real*8 , intent(out) :: CgXV(1:N_nodeNAC)
    
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "Not for stage 1, type 1, routine: readStrctPropsForContinuing_MoreAddedCell"
    endif
    
    N_itm = 3
    
    if (modelID==1) flnm='strctPropsAddedCellTN'
    if (modelID==2) flnm='strctPropsAddedCellNI'
    
    
    if (strctNum.le.9) then
       write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
    elseif (strctNum.gt.9 .and. strctNum.le.99) then
       write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
    endif
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=21,file=trim(adjustl(full_flnm)))
    
    do i=1,N_itm
       
       if (i==1) jmax=N_sprNAC
       if (i==2) jmax=N_cellNAC
       if (i==3) jmax=N_nodeNAC
       
       do j = 1,jmax
          
          if (i==1) read(21,*) ksV(j),l0V(j)
          if (i==2) read(21,*) kaV(j),A0V(j)
          if (i==3) read(21,*) CgXV(j)
          
       enddo
       
    enddo
    
    close(21)
    
  end subroutine readStrctPropsForContinuing_MoreAddedCell
  
  
  
  subroutine checkAreaRoutineChk
    implicit none
    integer :: Nnode_Poly
    
    real*8, allocatable  :: Nodes(:,:)
    real*8, allocatable  :: AreaVal(:)
    
    integer :: i,j,jmax
    integer :: areaNm,nodeNm
    
    interface
       real*8 function Polygon_Area(NodesA,Nnode_PolyA)
         implicit none
         integer, intent(in) :: Nnode_PolyA
         real*8 , intent(in) :: NodesA(1:Nnode_PolyA,1:3)
       end function Polygon_Area
    end interface
    
    N_node = 175
    N_cell = 29
    
    max_node_area = 11

    !!use deallocate(node_xy and other arrays) based on where its being called
    
    allocate(node_xy(1:N_node,1:2))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(AreaVal(1:N_cell))
    
    open(unit=34,file='nodesReadAddedCell.dat')
    open(unit=35,file='area_nodeWOtxtAC.dat')
    
    do i = 1,2
       
       if (i==1) jmax = N_node
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          if (i==1) read(34,*) node_xy(j,1:2)
          if (i==2) read(35,*) areaNm,area_node(j,0:max_node_area)
          if (i==1) write(*,*) node_xy(j,1:2)
          if (i==2) write(*,*) areaNm,area_node(j,0:max_node_area)
       enddo
       
    enddo
    
    close(34)
    close(35)
    
    
    do i = 1,N_cell
       Nnode_Poly = area_node(i,0)
       
       allocate(Nodes(1:Nnode_Poly,1:3))

       do j = 1,Nnode_Poly
          nodeNm = area_node(i,j)
          Nodes(j,1:2) = node_xy(nodeNm,1:2)
          
       enddo
       
       AreaVal(i) = Polygon_Area(Nodes,Nnode_Poly)
       write(*,*) AreaVal(i),i,"AreaVal"
       
       deallocate(Nodes)
       
    enddo
    
  end subroutine checkAreaRoutineChk
  
  subroutine updt_trnsfrm_addedCellNI_aftMet
    implicit none
    integer :: i,j,jmax
    integer :: N_itm
    integer :: nodeNm,sprNm,areaNm
    
    call allocate_and_storeBfrMetTrnsfrms
    
    call node_spr_aftMet
    call node_area_aftMet
    call spr_node_aftMet
    call spr_area_aftMet
    call area_node_aftMet
    call area_spr_aftMet
    
    !call get_Nlist
    !call print_all_trnsfrms_AftMet
    
  end subroutine updt_trnsfrm_addedCellNI_aftMet
  
  
  subroutine allocate_and_storeBfrMetTrnsfrms
    implicit none
    
    if (allctn_cntTrnsBfrMet==0) then
       
       allocate(node_sprBfrMet(1:N_node,0:max_spr_node)  , node_areaBfrMet(1:N_node,0:max_area_node))
       allocate(spr_nodeBfrMet(1:N_spr,0:max_node_spr)   , spr_areaBfrMet(1:N_spr,0:max_area_spr))
       allocate(area_nodeBfrMet(1:N_cell,0:max_node_area), area_sprBfrMet(1:N_cell,0:max_spr_area))
       allocate(NlistBfrMet(1:N_node,1:max_area_node,1:2))
       
       allctn_cntTrnsBfrMet = 1
       
    endif
    
    node_sprBfrMet  = node_spr
    node_areaBfrMet = node_area
    
    spr_nodeBfrMet  = spr_node
    spr_areaBfrMet  = spr_area
    
    area_nodeBfrMet = area_node
    area_sprBfrMet  = area_spr
    
    NlistBfrMet = Nlist
    
  end subroutine allocate_and_storeBfrMetTrnsfrms
  
  
  subroutine node_spr_aftMet
    implicit none
    integer :: lastTN
    integer :: prvNode
    integer :: FirstDNsPrv(1:2),FirstDNsNow(1:2)
    integer :: i,imax
    integer :: node1,node2
    
    
    lastTN = N_TNsAC
    
    imax = node_spr(lastTN,0)
    
    do i = 1,imax
       prvNode            = node_sprBfrMet(lastTn,i) 
       node_spr(lastTN,i) = prvNode - nsprs_InACell
    enddo
    
    FirstDNsPrv(1:2) = double_node(1,1:2)
    FirstDNsNow(1)   = double_node(1,1) - 2
    FirstDNsNow(2)   = double_node(1,2) - 2
    
    write(*,*) FirstDNsPrv(1:2),"Fdn prv"
    write(*,*) FirstDNsNow(1:2),"Fdn now"
    
    node1 = FirstDnsNow(1) ; node2 = FirstDNsNow(2)
    node_spr(node1,0) = node_sprBfrMet(node1,0) + node_sprBfrMet(node2,0) 
    
    node_spr(node1,1:2) = node_sprBfrMet(node1,1:2)
    node_spr(node1,3:4) = node_sprBfrMet(node2,1:2)
    node_spr(node1,5)   = node_sprBfrMet(node1,3)
    node_spr(node1,6)   = node_sprBfrMet(node2,3)
    
    node_spr(node2,0:max_spr_node) = node_spr(node1,0:max_spr_node)
    
  end subroutine node_spr_aftMet
  
  
  subroutine node_area_aftMet
    implicit none
    integer :: lastTN
    integer :: prvArea
    integer :: FirstDNsPrv(1:2),FirstDNsNow(1:2)
    integer :: i,imax
    integer :: node1,node2
    
    open(unit=27,file='node_areaTrns.dat')
    
    lastTN = N_TNsAC
    
    imax = node_area(lastTN,0)
    
    do i = 1,imax
       prvArea              = node_areaBfrMet(lastTn,i) 
       node_area(lastTN,i)  = prvArea - 1
    enddo
    
    FirstDNsPrv(1:2) = double_node(1,1:2)
    FirstDNsNow(1)   = double_node(1,1) - 2
    FirstDNsNow(2)   = double_node(1,2) - 2

    write(*,*) FirstDNsPrv(1:2),"FirstDNsPrv"
    write(*,*) FirstDNsNow(1:2),"FirstDNsNow"
    
    node1 = FirstDnsNow(1) ; node2 = FirstDNsNow(2)
    node_area(node1,0) = node_areaBfrMet(node1,0) + node_areaBfrMet(node2,0) 
    
    node_area(node1,1:2) = node_areaBfrMet(node1,1:2)
    node_area(node1,3:4) = node_areaBfrMet(node2,1:2)
    
    node_area(node2,0:max_area_node) = node_area(node1,0:max_area_node)
    
    close(27)
    
  end subroutine node_area_aftMet
  
  subroutine spr_node_aftMet
    implicit none
    integer :: lastTN
    integer :: areaLB,areaRB !areaLeftBfr,areaRightBfr
    integer :: areaLN,areaRN !areaLeftNow,areaRightNow
    integer :: apclLB,apclRB
    integer :: apclLN,apclRN
    
    integer :: i,j,jmax
    integer :: nodeNm
    
    lastTN = N_TNsAC
    
    areaLB = node_areaBfrMet(lastTN,1)    ; areaRB = node_areaBfrMet(lastTN,2)
    apclLB = (areaLB-1)*(nsprs_InACell)+1 ; apclRB = (areaRB-1)*(nsprs_InACell)+1
    
    write(*,*) areaLB,areaRB,"area B"
    write(*,*) apclLB,apclRB,"apcl B"
    
    do i = 1,2
       
       if (i==1) jmax = spr_node(apclLB,0)
       if (i==2) jmax = spr_node(apclRB,0)
       
       do j = 1,jmax
          
          if (i==1) nodeNm = spr_node(apclLB,j)
          if (i==2) nodeNm = spr_node(apclRB,j)
          
          if (nodeNm == lastTN) then
             
             if (i==1) then
                spr_node(apclLB,j)   = spr_nodeBfrMet(apclLB,j+1)
                spr_node(apclLB,j+1) = -1
                
                spr_node(apclLB,0)   = 2 
                exit
                
             elseif (i==2) then
                spr_node(apclRB,j)   = spr_nodeBfrMet(apclRB,j+1)
                spr_node(apclRB,j+1) = -1
                
                spr_node(apclRB,0)   = 2 
                exit  
                
             endif
             
          endif
          
       enddo
       
    enddo
    
    areaLN = areaLB-1                     ; areaRN = areaRB-1
    apclLN = (areaLN-1)*(nsprs_InACell)+1 ; apclRN = (areaRN-1)*(nsprs_InACell)+1
  
    do i = 1,2
       
       if (i==1) jmax = spr_node(apclLN,0)
       if (i==2) jmax = spr_node(apclRN,0)
       
       do j = 1,jmax
          
          if (i==1) nodeNm = spr_node(apclLN,j)
          if (i==2) nodeNm = spr_node(apclRN,j)
          
          if (j==jmax) then
           
             if (i==1) then
                
                spr_node(apclLN,j)   = lastTN
                spr_node(apclLN,j+1) = spr_nodeBfrMet(apclLN,j)
                
                spr_node(apclLN,0)   = 3
                exit
                
             elseif (i==2) then
                
                spr_node(apclRN,j)   = lastTN
                spr_node(apclRN,j+1) = spr_nodeBfrMet(apclRN,j)
                
                spr_node(apclRN,0)   = 3 
                exit  
                
             endif
             
          endif
          
       enddo
       
    enddo
    
    
  end subroutine spr_node_aftMet
  
  
  subroutine spr_area_aftMet
    implicit none
    
    continue !nothing will change
    
  end subroutine spr_area_aftMet
  
  
  subroutine area_node_aftMet
    implicit none
    integer :: lastTN
    integer :: areaLB,areaRB
    integer :: areaLN,areaRN
    
    integer :: nodeNm,areaNm,nodePostn
    integer :: nodePostnOfLastTNLft,nodePostnOfLastTNRgt
    integer :: i,j,jmax
    
    lastTN = N_TNsAC
    
    areaLB = node_areaBfrMet(lastTN,1) ; areaRB = node_areaBfrMet(lastTN,2)
    
    do i = 1,2
       
       if (i==1) areaNm = areaLB
       if (i==2) areaNm = areaRB
       
       jmax = area_nodeBfrMet(areaNm,0)
       
       do j = 1,jmax
          
          nodeNm = area_nodeBfrMet(areaNm,j)
          
          if (nodeNm == lastTN) then
             nodePostn = j
             
             if (i==1) nodePostnOfLastTNLft = nodePostn
             if (i==2) nodePostnOfLastTNRgt = nodePostn
             
             if (nodePostn == jmax) then
                area_node(areaNm,j) = -1
                
             elseif (nodePostn .ne. jmax) then
                area_node(areaNm,j:(jmax-1)) = area_nodeBfrMet(areaNm,(j+1):(jmax))
                area_node(areaNm,jmax)       = -1
                exit
             endif
             
          endif
          
       enddo
       
       area_node(areaNm,0) = area_nodeBfrMet(areaNm,0) - 1
       
    enddo
    
    areaLN = areaLB-1 ; areaRN = areaRB-1
    
    do i = 1,2
       
       if (i==1) areaNm = areaLN
       if (i==2) areaNm = areaRN
       
       jmax = max_node_area
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (j==nodePostnOfLastTNLft) then
                area_node(areaNm,j) = lastTN
                
                if (j == max_node_area) then
                   continue
                elseif (j .ne. max_node_area) then
                   area_node(areaNm,j) = lastTN
                   area_node(areaNm,(j+1):(jmax)) = area_nodeBfrMet(areaNm,(j):(jmax-1))
                   exit
                endif
                
             endif
             
          elseif (i==2) then
             
             if (j==nodePostnOfLastTNRgt) then
                area_node(areaNm,j) = lastTN
                
                if (j == max_node_area) then
                   continue
                elseif (j .ne. max_node_area) then
                   area_node(areaNm,j) = lastTN
                   area_node(areaNm,(j+1):(jmax)) = area_nodeBfrMet(areaNm,(j):(jmax-1))
                   exit
                endif
                
             endif
             
          endif
          
       enddo
       
       area_node(areaNm,0) = area_nodeBfrMet(areaNm,0) + 1
       
    enddo
    
  end subroutine area_node_aftMet
  
  
  subroutine area_spr_aftMet
    implicit none
    
    continue !nothing will change
    
  end subroutine area_spr_aftMet
  
  subroutine print_all_trnsfrms_AftMet
    implicit none
    integer :: i,j,jmax
    integer :: max
    integer :: N_itm
    
    N_itm = 2
    
    open(unit=145,file='TransfrmVarsAftMet.dat')
    
    do i = 1,N_itm
       
       if (i==1) write(145,fmt=*) "node_spr:"
       if (i==2) write(145,fmt=*) "node_area:"
       
       do j = 1,N_node
          
          if (i==1) write(145,fmt=*) j,node_spr(j,0:max_spr_node)
          if (i==2) write(145,fmt=*) j,node_area(j,0:max_area_node)
          
       enddo
       
       write(145,fmt=*) " "
       write(545,fmt=*) " "
       
    enddo
    
    
    do i = 1,N_itm
       
       if (i==1) write(145,*) "spr_node:"
       if (i==2) write(145,*) "spr_area:"
       
       do j = 1,N_spr
          
          if (i==1) write(145,*) j,spr_node(j,0:max_node_spr)   
          if (i==2) write(145,*) j,spr_area(j,0:max_area_spr)
          
       enddo
       
       write(145,*) " "
       
    enddo
    
    do i = 1,N_itm
       
       if (i==1) write(145,fmt=*) "area_node:"
       if (i==2) write(145,fmt=*) "area_spr:"
       
       do j = 1,N_cell
          
          if (i==1) write(145,*) j,area_node(j,0:max_node_area)
          if (i==2) write(145,*) j,area_spr(j,0:max_spr_area)
          
       enddo
       
       write(145,*) " "
       
    enddo

    write(145,*) "Nlist:"
    
    
    do i = 1,N_node
       
       write(145,*) "Nlist for  node_nm =",i
       
       jmax = node_area(i,0)
       
       do j = 1,jmax
          write(145,*) Nlist(i,j,1:2)
       enddo
       
       write(145,*) " "
    enddo
    
    close(145)
    
  end subroutine print_all_trnsfrms_AftMet
  
  
  subroutine get_decision_of_redefining_AddedCellNI(tol_Rdfn,lgcl_rdfn)
    implicit none
    real*8,  intent(in)    :: tol_Rdfn
    logical, intent(inout) :: lgcl_rdfn
    
    logical :: lgcl_NeighSwitch
    integer :: Pulley,LftNeigh,RghtNeigh
    integer :: LastCellAboveLft,LastCellAboveRght
    
    if ((modelID==2) .or.(AddedCellModelInitiation==1)) then
       continue
    else
       write(*,*) "sbrtn: get_decision_of_redefining_AddedCellV not compatible"
       stop
    endif
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    lgcl_rdfn = .False.
    
    open(unit=96,file='OriginDis.dat',position='append')
    
    call find_LastCellofTopLayerAddedCell(LastCellAboveLft,LastCellAboveRght)
    
    LftNeigh  = 2*(LastCellAboveLft)    - 1
    RghtNeigh = 2*(LastCellAboveRght+1) - 1
    write(*,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh NI"
    write(*,*) origin(1:2),"origin for NIadded"    
    
    Origin_dis(1:2) = dist_to_Origin(LftNeigh,RghtNeigh)
    write(*,*) Origin_dis(1:2),"Origin dis"
    write(96,fmt=*) Origin_dis(1:2)
    close(96)
    
    lgcl_NeighSwitch = .False.
    call Origin_Neighbour_switch(LftNeigh,RghtNeigh,lgcl_NeighSwitch)
    
    
    if (Origin_dis(1).lt.tol_Rdfn .and.Origin_dis(2).lt.tol_Rdfn) then
       write(*,*) "Both are less"
       lgcl_rdfn = .True.
       
    elseif (Origin_dis(1).lt.tol_Rdfn.and.Origin_dis(2).gt.tol_Rdfn) then
       write(*,*) "Lft is less"
       
    elseif (Origin_dis(1).gt.tol_Rdfn.and.Origin_dis(2).lt.tol_Rdfn) then
       write(*,*) "Rght is less"
       
    elseif (Origin_dis(1).gt.tol_Rdfn.and.Origin_dis(2).lt.tol_Rdfn) then
       write(*,*) "Both are big"
       
    elseif (lgcl_NeighSwitch .eqv. .True.) then
       lgcl_rdfn = .True.
       write(*,*) "Neighbour of Origin switched"
    endif
    
  end subroutine get_decision_of_redefining_AddedCellNI
  
  
  subroutine find_LastCellofTopLayerAddedCell(LastCellAboveLft,LastCellAboveRght)
    implicit none
    integer, intent(out) :: LastCellAboveLft,LastCellAboveRght
    integer :: LastTN
    
    LastTN = N_TNsAC

    LastCellAboveLft  = node_area(LastTN,1)
    LastCellAboveRght = node_area(LastTN,2)
    
  end subroutine find_LastCellofTopLayerAddedCell
  
  
  subroutine taking_lftrghtOriginNeigh_downwards_AC_and_updtSys(LftNeigh,RghtNeigh)
    implicit none
    integer, intent(inout) :: LftNeigh,RghtNeigh
    
    integer :: LastCellAboveLft,LastCellAboveRght 
    real*8  :: E,ydis_frmPulley
    
    open(unit=76,file='tkingLftRghtdownAndUpdt.dat')
    
    call updt_trnsfrm_addedCellNI_aftMet
    
    ydis_frmPulley = 0.01d0
    
    call find_LastCellofTopLayerAddedCell(LastCellAboveLft,LastCellAboveRght)
    
    LftNeigh  = 2*LastCellAboveLft  + 1
    RghtNeigh = 2*LastCellAboveRght + 3
    
    write(76,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh"
    
    node_xy(LftNeigh,1) = origin(1)
    node_xy(LftNeigh,2)  = origin(2) - ydis_frmPulley
    
    node_xy(RghtNeigh,1) = origin(1)
    node_xy(RghtNeigh,2) = origin(2) - ydis_frmPulley
    
    node_typ(LftNeigh)  = 1
    node_typ(RghtNeigh) = 1
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call get_Nlist
    call print_all_trnsfrms_AftMet
    
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    close(76)

    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E"
    
    call get_gradient(node_xy,l0,A0,gradient)
    call Find_Analytical_and_Numerical_Mismatch
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExprmntAdded,FrameAdded)
    
    
  end subroutine taking_lftrghtOriginNeigh_downwards_AC_and_updtSys

  
  ! subroutine get_IFProps_AtConnectingCycls
  !   implicit none
  !   real*8 :: kaF1(1:N_cell),kaI2(1:N_cell),kaIF(1:N_cell) !IF=InterFace
  !   real*8 :: A0F1(1:N_cell),A0I2(1:N_cell),A0IF(1:N_cell)
  !   real*8 :: ksF1(1:N_spr) ,ksI2(1:N_spr) ,ksIF(1:N_spr)
  !   real*8 :: l0F1(1:N_spr) ,l0I2(1:N_spr) ,l0IF(1:n_spr)
    
  !   integer :: FirstCellL,FirstCellR
  !   integer :: LastCellL ,LastCellR
    
  !   integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7
  !   integer :: N_lims,cyclNm
  !   integer, allocatable :: limVal(:)
    
  !   integer :: csFnsh,csStrt
  !   integer :: cyclBfr,cyclAft
  !   real*8  :: weight
  !   integer :: i,j,jmax
  !   integer :: IF_no
    
  !   real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_cell),A0_inpt(1:N_cell),A0_otpt(1:N_cell)
  !   real*8  :: ks_inpt(1:N_spr) ,ks_otpt(1:N_spr) ,l0_inpt(1:N_spr) ,l0_otpt(1:N_spr)
    
  !   open(unit=292,file='limsValInConnCycl.dat')

  !   weight = 0.50d0
    
  !   allocate(ka_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
  !   allocate(A0_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
  !   allocate(ks_cycl(1:N_addedCycl,1:N_CS,1:N_spr))
  !   allocate(l0_cycl(1:N_addedCycl,1:N_CS,1:N_spr))
    
    
  !   do cyclNm = 1,N_addedCycl
       
  !      FirstCellL = N_addedCellL - (cyclNm-1) 
  !      FirstCellR = Hlf_Ncell + N_addedCellR - (cyclNm-1)
       
  !      LastCellL  = FirstCellL + (Hlf_NcellNAC-2)
  !      LastCellR  = FirstCellR + (Hlf_NcellNAC-2)
       
  !      lim1 = FirstCellL-1
  !      lim2 = LastCellL
  !      lim3 = Hlf_Ncell
  !      lim4 = FirstCellR-1
  !      lim5 = LastCellR
  !      lim6 = N_cell-1
  !      lim7 = N_cell
       
  !      N_lims = 7
  !      allocate(limVal(1:N_lims))
       
  !      do i = 1,N_lims
  !         if (i==1) limVal(i) = lim1
  !         if (i==2) limVal(i) = lim2
  !         if (i==3) limVal(i) = lim3
  !         if (i==4) limVal(i) = lim4
  !         if (i==5) limVal(i) = lim5
  !         if (i==6) limVal(i) = lim6
  !         if (i==7) limVal(i) = lim7
  !      enddo
       
  !      write(292,*) limVal(1:7),"lims"
       
  !      hlfCycl = 1
  !      call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)

  !      ka_cycl(cyclNm,1,1:N_cell) = ka_inpt(1:N_cell)
  !      ka_cycl(cyclNm,2,1:N_cell) = ka_otpt(1:N_cell)
  !      A0_cycl(cyclNm,1,1:N_cell) = A0_inpt(1:N_cell)
  !      A0_cycl(cyclNm,2,1:N_cell) = A0_otpt(1:N_cell)
  !      ks_cycl(cyclNm,1,1:N_spr)  = ks_inpt(1:N_spr)
  !      ks_cycl(cyclNm,2,1:N_spr)  = ks_otpt(1:N_spr)
  !      l0_cycl(cyclNm,1,1:N_spr)  = l0_inpt(1:N_spr)
  !      l0_cycl(cyclNm,2,1:N_spr)  = l0_otpt(1:N_spr)
       
  !      hlfCycl = 2
  !      call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
       
  !      ka_cycl(cyclNm,3,1:N_cell) = ka_inpt(1:N_cell)
  !      ka_cycl(cyclNm,4,1:N_cell) = ka_otpt(1:N_cell)
  !      A0_cycl(cyclNm,3,1:N_cell) = A0_inpt(1:N_cell)
  !      A0_cycl(cyclNm,4,1:N_cell) = A0_otpt(1:N_cell)
  !      ks_cycl(cyclNm,3,1:N_spr)  = ks_inpt(1:N_spr)
  !      ks_cycl(cyclNm,4,1:N_spr)  = ks_otpt(1:N_spr)
  !      l0_cycl(cyclNm,3,1:N_spr)  = l0_inpt(1:N_spr)
  !      l0_cycl(cyclNm,4,1:N_spr)  = l0_otpt(1:N_spr)
       
  !      deallocate(limVal)
       
  !   enddo
    
  !   csFnsh  = 4 ; csStrt  = 1 
    
  !   do IF_no = 1,(N_addedCycl-1)
       
  !      cyclBfr = IF_no ; cyclAft = IF_no+1
       
  !      kaF1(1:N_cell) = ka_cycl(cyclBfr,csFnsh,1:N_cell)
  !      kaI2(1:N_cell) = ka_cycl(cyclAft,csStrt,1:N_cell)
       
  !      A0F1(1:N_cell) = A0_cycl(cyclBfr,csFnsh,1:N_cell)
  !      A0I2(1:N_cell) = A0_cycl(cyclAft,csStrt,1:N_cell)
       
  !      ksF1(1:N_spr)  = ks_cycl(cyclBfr,csFnsh,1:N_spr)
  !      ksI2(1:N_spr)  = ks_cycl(cyclAft,csStrt,1:N_spr)
       
  !      l0F1(1:N_spr)  = l0_cycl(cyclBfr,csFnsh,1:N_spr)
  !      l0I2(1:N_spr)  = l0_cycl(cyclAft,csStrt,1:N_spr)
       
  !      do i = 1,2
          
  !         if (i==1) jmax = N_cell
  !         if (i==2) jmax = N_spr
          
  !         do j = 1,jmax
             
  !            if (i==1) then
  !               kaIF(j) = weight*kaF1(j) + (1.0d0-weight)*kaI2(j)
  !               A0IF(j) = weight*A0F1(j) + (1.0d0-weight)*A0I2(j)
                
  !            elseif (i==2) then
  !               ksIF(j) = weight*ksF1(j) + (1.0d0-weight)*ksI2(j) 
  !               l0IF(j) = weight*l0F1(j) + (1.0d0-weight)*l0I2(j) 
  !            endif
             
  !         enddo
       
  !      enddo
       
  !      ka_cycl(cyclBfr,csFnsh,1:N_cell) = kaIF(1:N_cell)
  !      ka_cycl(cyclAft,csStrt,1:N_cell) = kaIF(1:N_cell)
  !      A0_cycl(cyclBfr,csFnsh,1:N_cell) = A0IF(1:N_cell)
  !      A0_cycl(cyclAft,csStrt,1:N_cell) = A0IF(1:N_cell)
       
  !      ks_cycl(cyclBfr,csFnsh,1:N_spr)  = ksIF(1:N_spr)
  !      ks_cycl(cyclAft,csStrt,1:N_spr)  = ksIF(1:N_spr)
  !      l0_cycl(cyclBfr,csFnsh,1:N_spr)  = l0IF(1:N_spr)
  !      l0_cycl(cyclAft,csStrt,1:N_spr)  = l0IF(1:N_spr)
       
  !   enddo
    
  !   close(292)
    
  ! end subroutine get_IFProps_AtConnectingCycls
  
end module Adding_cells
