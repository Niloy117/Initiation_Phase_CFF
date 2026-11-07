module switch_and_unswitch_models
  use system_parameters
  use transfrm_info
  use conversion_routines
  use mltiple_calls_together
  use storing_changing_restoring_routines
  
  use Wfunc_and_its_derivative
  use changing_parameters
  use storing_changing_restoring_routines
  
  implicit none
  
  integer :: Exprmnt_NI=37,  Frame_NI=1    !67 !change 37 to 30
  integer :: ExprmntAdded=32,FrameAdded=1
  integer :: Exprmnt_NIVF=38,Frame_NIVF=1
  integer :: Exprmnt_Init=40,FramePrdctnInit=1
  
  integer :: Nstp_PT=10
  integer :: NI_saveOrNot
  integer :: prtnVal
  integer :: N_progCycl
  integer :: ProgCyclNo,strctVno,strctVno1,strctVno2
  integer :: NI_alrdySaved
  integer :: allctn_Continued_NI_trnsfrmArray=0
  
  real*8,  allocatable :: node_xyNIsaved(:,:)
  integer, allocatable :: typ_sprNIsaved(:)
  real*8,  allocatable :: k_sprNIsaved(:) ,l0_NIsaved(:),l_NIsaved(:)
  real*8,  allocatable :: k_areaNIsaved(:),A0_NIsaved(:),A_NIsaved(:)
  
  real*8, allocatable  :: A0_stepF(:),ka_stepF(:) !with F=first, from strc 1-3
  real*8, allocatable  :: l0_stepF(:),ks_stepF(:) !if needed will design frm 3-2
  real*8, allocatable  :: A0_I(:),ka_I(:),l0_I(:),ks_I(:)
  real*8, allocatable  :: A0_F(:),ka_F(:),l0_F(:),ks_F(:)
  
  real*8, allocatable  :: ksvNI(:),kavNI(:),l0vNI(:),A0cNI(:),CgXvNI(:),CgYvNI(:)
  real*8, allocatable  :: ksv2NI(:),kav2NI(:),l0v2NI(:),A0v2NI(:),CgXv2NI(:),CgYv2NI(:)
  
contains
  
  subroutine switchto_NI_model_run_and_switchbackto_TN
    implicit none
    real*8  :: E
    real*8  :: A0_str1,A0_str2,A0_str3,A0_str4
    real*8  :: l0_str1,l0_str2,l0_str3
    integer :: i
    real*8  :: Ea_aft,Ea_bfr,Num_grd
    real*8  :: Es,Ea,Eb,Eg
    real*8  :: A1,A2,A3
    
    if (nodeInsrtnStrts==1) continue
    if (nodeInsrtnStrts.ne.1) then
       write(*,*) "stop Reason: nodeInsrtnStrts must be zero to continue"
       write(*,*) "sb:switchto_NI_model_run_and_switchbackto_TN, fl:node_insertion"
       stop
    endif
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0) ; write(*,*) E,"E" ; call get_gradient(node_xy,l0,A0,gradient)
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    
    Frame_NI = Frame_NI+1
    call Find_Analytical_and_Numerical_Mismatch
    
    call sleep(1)
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine switchto_NI_model_run_and_switchbackto_TN
  
  subroutine Equilibrate_AddedCell_model
    implicit none
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExprmntAdded,FrameAdded)
    
    FrameAdded = FrameAdded+1
    write(*,*) (FrameAdded-1),"in FrameAdded"
    
  end subroutine Equilibrate_AddedCell_model
  
  subroutine Equilibrate_only_NI_model
    implicit none
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    
    Frame_NI = Frame_NI+1
    write(*,*) (Frame_NI-1),"in NI"
    
  end subroutine Equilibrate_only_NI_model
  
  subroutine Equilibrate_only_NI_model_withVF_region
    implicit none
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NIVF,Frame_NIVF)
    
    Frame_NIVF = Frame_NIVF+1
    write(*,*) (Frame_NIVF-1),"in NIVF"
    
  end subroutine Equilibrate_only_NI_model_withVF_region
  
  subroutine WO_Equilibrate_onlysaveConfig_withVF_region
    implicit none
    
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NIVF,Frame_NIVF)
    
    Frame_NIVF = Frame_NIVF+1
    write(*,*) (Frame_NIVF-1),"in NIVF"
    
  end subroutine WO_Equilibrate_onlysaveConfig_withVF_region
  
  
  subroutine Equilibrate_only_TN_model(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    
    FrmNo = FrmNo+1
    write(*,*) (FrmNo-1),"in only TN"
    
  end subroutine Equilibrate_only_TN_model
  
  
  subroutine save_only_NI_config
    implicit none
    
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    
    Frame_NI = Frame_NI+1
    write(*,*) (Frame_NI-1),"in NI"
    
  end subroutine save_only_NI_config
  
  subroutine switchto_NI_giveTheProps_run_and_switchbackto_TN()
    implicit none
    
  end subroutine switchto_NI_giveTheProps_run_and_switchbackto_TN
  
  subroutine switch_to_NI_and_Equilibrate
    implicit none
    real*8  :: E
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E in sb:switch_to_NI_and_Equilibrate"
    call get_gradient(node_xy,l0,A0,gradient)
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    
  end subroutine switch_to_NI_and_Equilibrate
  
  subroutine aftrEqulibrate_backto_TN
    implicit none
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
  end subroutine aftrEqulibrate_backto_TN
  
  subroutine switch_to_NI_and_not_Equilibrate
    implicit none
    real*8  :: E
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E in switch_to_NI_and_not_Equilibrate"
    call get_gradient(node_xy,l0,A0,gradient)
    
    !call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    
  end subroutine switch_to_NI_and_not_Equilibrate
  
  subroutine not_Equlibrate_backto_TN
    implicit none
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
  end subroutine not_Equlibrate_backto_TN
  
  
  subroutine switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr(saveFlnm,sgmntdOrNot,WhereTo)
    implicit none
    character(len=100), intent(in) :: saveFlnm
    integer, intent(in) :: sgmntdOrNot,whereTo
    real*8              :: E
    
    if (nodeInsrtnStrts==1) continue
    if (nodeInsrtnStrts.ne.1) then
       write(*,*) "stop Reason: nodeInsrtnStrts must be zero to continue"
       write(*,*) "sb:switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr, fl:node_insertion"
       stop
    endif
    
    if (sgmntdOrNot==0) then
       
       if (whereTo.ne.0) then
          write(*,*) "stop Reason: if sgmntd is 0, whereTo mustbe 0"
          stop
       endif
       
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       write(*,*) E,"E"
       call get_gradient(node_xy,l0,A0,gradient)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       
       Frame_NI = Frame_NI+1
       call Find_Analytical_and_Numerical_Mismatch

       call sleep(1)
       call write_in_given_file(saveFlnm)
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (sgmntdOrNot==1) then
       
       if (whereTo==1) then
          call switch_to_NI_model
          E = Energy(node_xy,l0,A0)
          write(*,*) E,"E"
          call get_gradient(node_xy,l0,A0,gradient)
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          
          Frame_NI = Frame_NI+1
          call Find_Analytical_and_Numerical_Mismatch

          call sleep(1)
          
       elseif (whereTo==2) then
          call write_in_given_file(saveFlnm)
          call deallocate_repetitive_arrays
          call switchback_to_TN_model
       endif
       
    endif
    
  end subroutine switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr
  
  
  
  subroutine switchto_NI_model_run_switchbackto_TN_seprtlyOrtogthrNOFlnm(sgmntdOrNot,WhereTo)
    implicit none
    integer, intent(in) :: sgmntdOrNot,whereTo
    real*8              :: E
    
    if (nodeInsrtnStrts==1) continue
    if (nodeInsrtnStrts.ne.1) then
       write(*,*) "stop Reason: nodeInsrtnStrts must be zero to continue"
       write(*,*) "sb:switchto_NI_model_run_save_switchbackto_TN_seprtlyOrtogthr, fl:node_insertion"
       stop
    endif
    
    if (sgmntdOrNot==0) then
       
       if (whereTo.ne.0) then
          write(*,*) "stop Reason: if sgmntd is 0, whereTo mustbe 0"
          stop
       endif
       
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       write(*,*) E,"E"
       call get_gradient(node_xy,l0,A0,gradient)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       
       Frame_NI = Frame_NI+1
       call Find_Analytical_and_Numerical_Mismatch

       call sleep(1)
       call deallocate_repetitive_arrays
       call switchback_to_TN_model

    elseif (sgmntdOrNot==1) then
       
       if (whereTo==1) then
          call switch_to_NI_model
          E = Energy(node_xy,l0,A0)
          write(*,*) E,"E"
          call get_gradient(node_xy,l0,A0,gradient)
          
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          
          Frame_NI = Frame_NI+1
          call Find_Analytical_and_Numerical_Mismatch

          call sleep(1)
          
       elseif (whereTo==2) then
          call deallocate_repetitive_arrays
          call switchback_to_TN_model
       endif
       
    endif
    
  end subroutine switchto_NI_model_run_switchbackto_TN_seprtlyOrtogthrNOFlnm
  
  
  
  subroutine switch_to_NI_model
    implicit none
    integer :: cnt_lp
    
    modelID   = 2
    NAEC      = 2
    NAEC_Apcl = 0 ; NAEC_Bsal = NAEC ; NAEC_Ltrl = NAEC
    
    call store_TN_SysVars_and_Arrays
    call get_NI_Sysvars
    call deallocate_and_reallocate_TN_Arrays
    
    call get_NodeVars
    call get_SprVars
    call get_kphi_and_nodePhiTyp
    call get_CgXvars
    call get_CgYvars
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    !stop
    call store_TN_trnsfrms
    call deallocate_and_reallocate_for_NI_trnsfrmVars
    call get_all_the_transfrms
    call print_all_trnsfrms
    !stop
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
  end subroutine switch_to_NI_model
  
  
  subroutine switch_to_NI_model_withNodesInAllmembranes
    implicit none
    
    NI_incldAS_woPP = 1
    modelID         = 2
    
    NAEC      = 2
    NAEC_Apcl = NAEC ; NAEC_Bsal = NAEC ; NAEC_Ltrl = NAEC
    
    call store_TN_SysVars_and_Arrays
    call get_NI_Sysvars_withNodesInAllSides
    call deallocate_and_reallocate_TN_Arrays
    
    call get_NodeVars
    call get_SprVars
    call get_kphi_and_nodePhiTyp
    call get_CgXvars
    call get_CgYvars
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_TN_trnsfrms
    call deallocate_and_reallocate_for_NI_trnsfrmVars
    call get_all_the_transfrms
    call print_all_trnsfrms
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    call print_prop_of_system
    
  end subroutine switch_to_NI_model_withNodesInAllmembranes
  
  subroutine switch_to_NIsaved_model
    implicit none
    real*8 :: E
    
    modelID = 2
    NAEC    = 2
    
    call store_TN_SysVars_and_Arrays
    call get_NI_Sysvars
    call deallocate_and_reallocate_TN_Arrays
    
    call get_NodeVars_NIsaved
    call get_SprVars_NIsaved
    call get_AreaVarsNIsaved
    call get_kphi_and_nodePhiTyp
    call get_CgXvars
    
    call store_all_moving_coordnte_variables
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    !stop
    call store_TN_trnsfrms
    call deallocate_and_reallocate_for_NI_trnsfrmVars
    call get_all_the_transfrms
    call print_all_trnsfrms
    !stop
    
    call store_all_gradient_variables
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
    !E = Energy(node_xy,l0,A0)
    !write(*,*) E,"E"
    !call get_gradient(node_xy,l0,A0,gradient)
    
  end subroutine switch_to_NIsaved_model
  
  subroutine allocate_and_initialize_NIsavedVars
    implicit none
    
    allocate(node_xyNIsaved(1:N_node,1:N_dmnsn))
    allocate(typ_sprNIsaved(1:N_spr))
    allocate(k_sprNIsaved(1:N_spr))
    allocate(l0_NIsaved(1:N_spr),l_NIsaved(1:N_spr))
    allocate(k_areaNIsaved(1:N_cell))
    allocate(A0_NIsaved(1:N_cell),A_NIsaved(1:N_cell))
    
    node_xyNIsaved = -1.0d030
    typ_sprNIsaved = -1
    k_sprNIsaved   = -1.0d30
    l0_NIsaved     = -1.0d30
    l_NIsaved      = -1.0d30
    k_areaNIsaved  = -1.0d30
    A0_NIsaved     = -1.0d30
    A_NIsaved      = -1.0d30
    
  end subroutine allocate_and_initialize_NIsavedVars
  
  
  subroutine get_NIsaved_Properties
    implicit none

    node_xyNIsaved(1:N_node,1:N_dmnsn) = node_xy(1:N_node,1:N_dmnsn)
    
    typ_sprNIsaved(1:N_spr) = typ_spr(1:N_spr)
    k_sprNIsaved(1:N_spr)   = k_spr(1:N_spr)
    l0_NIsaved(1:N_spr)     = l0(1:N_spr)
    l(1:N_spr)              = l(1:N_spr)
    k_areaNIsaved(1:N_cell) = k_area(1:N_cell)
    A0_NIsaved(1:N_cell)    = A0(1:N_cell)
    A_NIsaved(1:N_cell)     = A(1:N_cell)
    
  end subroutine get_NIsaved_Properties
  
  subroutine store_TN_SysVars_and_Arrays
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
    
  end subroutine store_TN_SysVars_and_Arrays
  
  
  subroutine get_NI_Sysvars
    implicit none
    
    N_node = N_node + NAEC*N_curve
    N_spr  = N_spr  + NAEC*N_curve 
    
  end subroutine get_NI_Sysvars
  
  
  subroutine get_NI_Sysvars_withNodesInAllSides
    implicit none
    
    call get_N_curve
    
    write(*,*) N_node,N_spr,N_curve,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,N_cell,"inside NI_Sysvars BFR"
    
    N_node = N_node + (NAEC)*(N_curve)
    N_spr  = N_spr  + (NAEC)*(N_curve) 
    
    write(*,*) N_node,N_spr,N_curve,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,N_cell,"inside NI_Sysvars AFT"
    
  end subroutine get_NI_Sysvars_withNodesInAllSides
  
  
  subroutine deallocate_and_reallocate_TN_Arrays
    implicit none
    
    deallocate(node_xy,node_typ,node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(nodePhi_typ,k_phi,CgXNode,CgYNode)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_TN_Arrays
  
  subroutine get_NodeVars
    implicit none
    integer :: N_addedNodes
    integer :: i,j,cnt
    integer :: N1,N2
    real*8  :: node1(1:2),node2(1:2)
    integer :: spr_nm
    
    real*8, allocatable :: IntrmdNodes(:,:)
    
    node_xy(1:N_nodeS,1:N_dmnsn) = node_xyStr(1:N_nodeS,1:N_dmnsn)
    
    N_addedNodes = NAEC*N_curve
    cnt          = N_nodeS
    
    allocate(IntrmdNodes(1:NAEC,1:2))       ; IntrmdNodes     = -1.d20
    allocate(insrtd_vrtxNode(1:N_node,1:2)) ; insrtd_vrtxNode = -1
    
    call get_curve_spr_trnsfrm_and_is_it_a_curve
    
    do i = 1,N_curve
       
       spr_nm = curve_spr(i)
       
       call get_trmnlNode_nums_frmSprNm(spr_nm,N1,N2)
       call get_xy_of_trmnlNodes(spr_nm,N1,N2,node1,node2)
       call get_allIntrmdNodes(node1,node2,NAEC,IntrmdNodes)
       
       if (typ_sprStr(spr_nm)==8) then
          write(*,*) IntrmdNodes(1,1:2),"Intrmd 1" ; write(*,*) IntrmdNodes(2,1:2),"Intrmd 2"
       endif   
       
       do j = 1,NAEC
          cnt = cnt+1
          node_xy(cnt,1:2)       = IntrmdNodes(j,1:2) !; write(*,*) cnt,node_xy(cnt,1:2),"cnt,node_xy"
          insrtd_vrtxNode(cnt,1) = N1
          insrtd_vrtxNode(cnt,2) = N2
       enddo
       
    enddo
    
    !call shift_the_IntrmdNodes_to_help_equilibrium
    
    call get_node_typs_for_insrtdNodes
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars
  
  subroutine get_nodesInASide(CellNum,ApBsLt,NumNodesInASide,nodesInASide)
    implicit none
    integer, intent(in)  :: CellNum,ApBsLt,NumNodesInASide
    integer, intent(out) :: nodesInASide(1:NumNodesInASide)
    
    integer :: i,j
    integer :: INperCell

    INperCell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) ; write(*,*) INperCell,"INperCell"
    
    do i = 1,NumNodesInASide
       
       if (i==1) then
          if (CellNum.le.Hlf_Ncell) then
             if (ApBsLt==1) nodesInASide(i) = (CellNum-1)*2+1
             if (ApBsLt==2) nodesInASide(i) = (CellNum-1)*2+2
             if (ApBsLt==3) nodesInASide(i) = (CellNum-1)*2+3
          elseif ((CellNum.gt.Hlf_Ncell).and.(CellNum.le.(2*Hlf_Ncell))) then
             if (ApBsLt==1) nodesInASide(i) = (CellNum)*2+1
             if (ApBsLt==2) nodesInASide(i) = (CellNum)*2+2
             if (ApBsLt==3) nodesInASide(i) = (CellNum)*2+3
          endif
       elseif (i==NumNodesInASide) then
           if (CellNum.le.Hlf_Ncell) then
             if (ApBsLt==1) nodesInASide(i) = (CellNum)*2+1
             if (ApBsLt==2) nodesInASide(i) = (CellNum)*2+2
             if (ApBsLt==3) nodesInASide(i) = (CellNum)*2+2
          elseif ((CellNum.gt.Hlf_Ncell).and.(CellNum.le.(2*Hlf_Ncell))) then
             if (ApBsLt==1) nodesInASide(i) = (CellNum+1)*2+1
             if (ApBsLt==2) nodesInASide(i) = (CellNum+1)*2+2
             if (ApBsLt==3) nodesInASide(i) = (CellNum+1)*2+2
          endif
       else
          if (ApBsLt==1) nodesInASide(i)= numRegulrNode+(CellNum-1)*(INperCell)+(i-1)  
          if (ApBsLt==2) nodesInASide(i)= numRegulrNode+(CellNum-1)*(INperCell)+(NAEC_Apcl)+(i-1)
          if (ApBsLt==3) nodesInASide(i)= numRegulrNode+(CellNum-1)*(INperCell)+(NAEC_Apcl+NAEC_Bsal)+(i-1)
       endif
       
    enddo
    
  end subroutine get_nodesInASide
  
   subroutine get_nodesInASideOfIC(CellNum,ApBsLt,NumNodesInASide,nodesInASide)
    implicit none
    integer, intent(in)  :: CellNum,ApBsLt,NumNodesInASide
    integer, intent(out) :: nodesInASide(1:NumNodesInASide)
    
    integer :: i,j
    integer :: INperCell
    
    if (ApBsLt==3) stop 'IC doesnt have LT side'
    
    INperCell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) ; write(*,*) INperCell,"INperCell"
    
    do i = 1,NumNodesInASide
       
       if (i==1) then
          if (ApBsLt==1) nodesInASide(i) = (Hlf_Ncell)*2+1
          if (ApBsLt==2) nodesInASide(i) = (Hlf_Ncell)*2+2
       elseif (i==NumNodesInASide) then
          if (ApBsLt==1) nodesInASide(i) = (2*(Hlf_Ncell+1))*2-1
          if (ApBsLt==2) nodesInASide(i) = (2*(Hlf_Ncell+1))*2
       else
          if (ApBsLt==1) nodesInASide(i)= numRegulrNode+(2*Hlf_Ncell)*(INperCell)+(i-1)  
          if (ApBsLt==2) nodesInASide(i)= numRegulrNode+(2*Hlf_Ncell)*(INperCell)+(NAEC_Apcl)+(i-1)
       endif
       
    enddo
    
  end subroutine get_nodesInASideOfIC
  
    
  subroutine get_NodeVars_NIsaved
    implicit none
    integer :: N_addedNodes
    integer :: i,j,cnt
    integer :: N1,N2
    real*8  :: node1(1:2),node2(1:2)
    integer :: spr_nm
    
    real*8, allocatable :: IntrmdNodes(:,:)
    
    node_xy(1:N_node,1:N_dmnsn) = node_xyNIsaved(1:N_node,1:N_dmnsn)
    
    N_addedNodes = NAEC*N_curve
    cnt = N_nodeS
    
    allocate(IntrmdNodes(1:NAEC,1:2))
    IntrmdNodes = -1.d20
    allocate(insrtd_vrtxNode(1:N_node,1:2))
    insrtd_vrtxNode = -1
    
    call get_curve_spr_trnsfrm_and_is_it_a_curve
    
    do i = 1,N_curve
       spr_nm = curve_spr(i)
       
       if (spr_node(spr_nm,0)==2) then
          N1 = spr_node(spr_nm,1)
          N2 = spr_node(spr_nm,2)
       endif
       
       !write(*,*) N1,N2,"N1-N2"
       
       if (typ_sprStr(spr_nm).ne.8) then
          node1(1:2) = node_xy(N1,1:2)
          node2(1:2) = node_xy(N2,1:2)
          
       elseif (typ_sprStr(spr_nm)==8) then
          node1(1:2) = node_xy(N2,1:2)
          node2(1:2) = node_xy(N1,1:2)
          write(*,*) N1,N2,"N1,N2 of typ_spr = 8"
          write(*,*) node1(1:2),"xy of N2"
          write(*,*) node2(1:2),"xy of N1"
       endif
       
       
       call get_allIntrmdNodes(node1,node2,NAEC,IntrmdNodes)
       
       if (typ_sprStr(spr_nm)==8) then
          write(*,*) IntrmdNodes(1,1:2),"Intrmd 1"
          write(*,*) IntrmdNodes(2,1:2),"Intrmd 2"
       endif   
       
       do j = 1,NAEC
          cnt = cnt+1
          !node_xy(cnt,1:2) = IntrmdNodes(j,1:2)
          !write(*,*) cnt,node_xy(cnt,1:2),"cnt,node_xy"
          insrtd_vrtxNode(cnt,1) = N1
          insrtd_vrtxNode(cnt,2) = N2
       enddo
       
    enddo
    
    
    node_typ(1:N_nodeS) = node_typStr(1:N_nodeS)
    node_typ((N_nodeS+1):N_node) = 1
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_NodeVars_NIsaved
  
  subroutine get_trmnlNode_nums_frmSprNm(spr_nm,N1,N2)
    implicit none
    integer, intent(in)  :: spr_nm
    integer, intent(out) :: N1,N2
    integer              :: nodesIntheSpr
    
    nodesIntheSpr = spr_node(spr_nm,0)
    
    if (nodesIntheSpr == 2) then
       N1 = spr_node(spr_nm,1)
       N2 = spr_node(spr_nm,2)  ! ; write(*,*) spr_nm,N1,N2,"N1-N2"
       
    elseif(nodesIntheSpr.gt.2) then
       write(*,*) spr_nm,((spr_nm/3)+1),mod(spr_nm,3),nodesIntheSpr,"sprInfo for nodesIntheSpr>2"
       stop 'not for nodesInthespr > 2'
    endif
    
  end subroutine get_trmnlNode_nums_frmSprNm
  
  subroutine get_xy_of_trmnlNodes(spr_nm,N1,N2,node1,node2)
    implicit none
    integer, intent(in)  :: spr_nm,N1,N2
    real*8 , intent(out) :: node1(1:N_dmnsn),node2(1:N_dmnsn)
    
    !%%%%%% typ_sprStr = 8 is no more relevant now, was previously designed for stage4 cases %%%%%%
    
    if (typ_sprStr(spr_nm).ne.8) then
       node1(1:2) = node_xy(N1,1:2)
       node2(1:2) = node_xy(N2,1:2)
       
    elseif (typ_sprStr(spr_nm)==8) then
       node1(1:2) = node_xy(N2,1:2)
       node2(1:2) = node_xy(N1,1:2)
       write(*,*) N1,N2,"N1,N2 (typ_spr=8)"; write(*,*)node1(1:2),node2(1:2),"xy of N2,N1 for typ=8"
    endif
    
  end subroutine get_xy_of_trmnlNodes
  
  subroutine get_allIntrmdNodes(node1,node2,No_ofNodes,IntrmdNodes)
    implicit none
    integer, intent(in)  :: No_ofNodes 
    real*8,  intent(in)  :: node1(1:2),node2(1:2)
    real*8,  intent(out) :: IntrmdNodes(1:No_ofNodes,1:2)
    
    integer :: i
    real*8  :: m1,m2,dm1,dm2
    real*8  :: nodeI(1:2)
    
    !write(*,*) node1(1:2),"node1"
    !write(*,*) node2(1:2),"node2"
    
    dm1 = 0.0d0 ; dm2 = 0.0d0
    
    do i = 1,No_ofNodes
       m1 = 1.0d0 + dm1
       m2 = real(No_ofNodes) - dm2
       
       call intrmdNode(node1,node2,m1,m2,nodeI)
       IntrmdNodes(i,1:2) = nodeI(1:2)
       
       dm1 = real(i)
       dm2 = real(i)
    enddo
    
  end subroutine get_allIntrmdNodes
  
  subroutine intrmdNode(node1,node2,m1,m2,nodeIntrmd)
    implicit none
    real*8, intent(in)  :: node1(1:2),node2(1:2)
    real*8, intent(in)  :: m1,m2
    real*8, intent(out) :: nodeIntrmd(1:2)
    
    real*8 :: x1,y1,x2,y2,xm,ym
    
    !write(*,*) node1(1:2),"node1"
    !write(*,*) node2(1:2),"node2"
    
    x1 = node1(1) ; y1 = node1(2)
    x2 = node2(1) ; y2 = node2(2)
    
    xm = (m1*x2+m2*x1)/(m1+m2)
    ym = (m1*y2+m2*y1)/(m1+m2)
    
    nodeIntrmd(1) = xm
    nodeIntrmd(2) = ym
    
    !write(*,*) nodeIntrmd(1:2),"nodeIntrmd"
    
  end subroutine intrmdNode
  
  subroutine get_node_typs_for_insrtdNodes()
    implicit none
    integer :: i,istrt,iend,imax,j,jmax,jLp
    integer :: cnt_apcl,cnt_bsal,cnt_ltrl
    integer :: cellNm,cnt_node,INinOneSide
    integer :: apclINstrt,apclINfnsh
    integer :: bsalINstrt,bsalINfnsh
    integer :: ltrlINstrt,ltrlINfnsh
    integer :: cellNmBfrFreeAS
    integer :: nodeStrt,nodeFnsh
    
    node_typ(1:N_nodeS) = node_typStr(1:N_nodeS)
    
    write(*,*) "check Hlf_Ncell value = ",Hlf_Ncell
    if (Hlf_Ncell .ne.11) stop 'inconsistent H;f_Ncell'
    
    cnt_apcl=0 ; cnt_bsal=0 ; cnt_ltrl=0
    INinOneSide = (Hlf_Ncell)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) ; write(*,*) "INinOneSide = ",INinOneSide
    
    if (NAEC_Apcl.ne.0) cnt_apcl=1 
    if (NAEC_Bsal.ne.0) cnt_bsal=1
    if (NAEC_Ltrl.ne.0) cnt_ltrl=1
    
    node_typ(1:N_nodeS) = node_typStr(1:N_nodeS)
    istrt = N_nodeS+1   ; iend = N_node
    
    if (IN_apclSideCaseNo==1) cellNmBfrFreeAS=Hlf_Ncell-addedNCPair
    write(*,*)"IN_apclSideCaseNo = ",IN_apclSideCaseNo
    
    do i = 1,(Hlf_Ncell)
       cellNm = i
       
       nodeStrt = (N_nodeS)+1+(i-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
       nodeFnsh = (N_nodeS)+0+(i-0)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)
       write(*,*) nodeStrt,nodeFnsh,"nodeStrt-nodeFnsh"
       
       apclINstrt = (nodeStrt+0)                   ; apclINfnsh = (nodeStrt+NAEC_Apcl-1)  
       bsalINstrt = (nodeStrt+NAEC_Apcl)           ; bsalINfnsh = (nodeStrt+NAEC_Apcl+NAEC_Bsal-1)
       ltrlINstrt = (nodeStrt+NAEC_Apcl+NAEC_Bsal) ; ltrlINfnsh = (nodeStrt+NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl-1)
       
       write(*,*) apclINstrt,apclINfnsh,"Apcl SF"
       write(*,*) bsalINstrt,bsalINfnsh,"Bsal_SF"
       write(*,*) ltrlINstrt,ltrlINfnsh,"Ltrl_SF"
       write(*,*) " "
       
       
       if (cnt_apcl==1) then
          
          if (cellNm.le.cellNmBfrFreeAS) node_typ(apclINstrt:apclINfnsh) = 2 ! x free, y fixed
          if (cellNm.gt.cellNmBfrFreeAS) node_typ(apclINstrt:apclINfnsh) = 1 ! free
          
       endif
       
       if (cnt_bsal==1) node_typ(bsalINstrt:bsalINfnsh) = 1 ! Basal and Lateral nodes are always 
       if (cnt_ltrl==1) node_typ(ltrlINstrt:ltrlINfnsh) = 1 ! free nodes, that's why 1.
       
    enddo
    
    node_typ((N_nodeS+INinOneSide+1):(N_nodeS+INinOneSide+INinOneSide)) = node_typ((N_nodeS+1):(N_nodeS+INinOneSide)) 
    node_typ((N_nodeS+(2*INinOneSide)+1):N_node) = 1 ! Initiator Cell
    
  end subroutine get_node_typs_for_insrtdNodes
  
  subroutine get_kphi_and_nodePhiTyp
    implicit none
    integer :: i
    real*8  :: k_phiVal
    
    k_phiVal                       = k_phiStr(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal
    
    
    do i = 1,N_node
       if (i.le.N_nodeS) then
          nodePhi_typ(i) = 1
       elseif (i.gt.N_nodeS) then
          nodePhi_typ(i) = 2
       endif
    enddo
    
    do i=1,N_node
       if (nodePhi_typ(i)==2) then
          k_phi(i,1:max_Phi_node) = 0.00d0
       endif
    enddo
    
    open(unit=743,file='kphi_and_nodePhi.dat')
    
    do i = 1,N_node
       write(743,*) nodePhi_typ(i),k_phi(i,1:max_Phi_node),i,"nodePhi,k_phi,i"
    enddo
    
    close(743)
    
  end subroutine get_kphi_and_nodePhiTyp
  
  subroutine get_curve_spr_trnsfrm_and_is_it_a_curve !might be shiftd elsewhere
    implicit none
    
    if (NI_incldAS_woPP == 0) call get_curve_spr_trnsfrm_and_is_it_a_curve_methd1
    if (NI_incldAS_woPP == 1) call get_curve_spr_trnsfrm_and_is_it_a_curve_methd2

    write(*,*) "NI_incldAS_woPP = ",NI_incldAS_woPP
    if ((NI_incldAS_woPP.ne.0) .and. (NI_incldAS_woPP.ne.1)) stop 'incompatible NI_incldAS_woPP'
    
  end subroutine get_curve_spr_trnsfrm_and_is_it_a_curve
  
  subroutine get_curve_spr_trnsfrm_and_is_it_a_curve_methd1
    implicit none
    integer :: i,j,jmax,jEnd
    integer :: spr_nm,typ_nm
    integer :: cnt,cnt_Pspr !Pspr=P  
    
    allocate(curve_spr(1:N_curve))
    allocate(is_it_a_curve(1:N_sprS))
    
    curve_spr     = -1
    is_it_a_curve = 0
    
    cnt  = 0
    jmax = 3
    cnt_Pspr = 0
    
    do i = 1,N_cell  
       !write(*,*) " "
       
       if (i==N_cell) then
          
          if (stageNo==1 .and. stageType==1) then
             if (CellsMeet==0)   jEnd = 2
             if (CellsMeet.gt.0) jEnd = 1
          elseif (stageNo==4 .and. stageType==1) then
             jEnd = 1
          elseif (stageNo==4 .and. stageType==2) then
             jEnd = jmax
          endif
          
       elseif (i.ne.N_cell) then
          jEnd = jmax
       endif
       
       
       if (stageNo==1 .and. stageType==1) then
          !write(*,*) jEnd,"jEnd"
          
          do j = 1,jEnd
             
             spr_nm = (i-1)*jmax + j
             typ_nm = typ_sprStr(spr_nm)
             !write(*,*) spr_nm,typ_nm,"spr,typ"
             
             if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                cnt = cnt+1
                !write(*,*) cnt,spr_nm,"cnt,spr"
                curve_spr(cnt) = spr_nm
                is_it_a_curve(spr_nm) = 1
             endif
             
          enddo
          
       elseif (stageNo==2 .or. stageNo==3) then
          write(*,*) "fl:node_insertion_model,sb:get_curve_spr_trnsfrm"
          stop
          
       elseif (stageNo==4) then
          
          do j = 1,jEnd
             
             spr_nm = (i-1)*jmax + j
             typ_nm = typ_sprStr(spr_nm)
             !write(*,*) spr_nm,typ_nm,"spr_nm,typ_nm"
             
             if (i.le.(ncl+ncr)) then
                
                if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                   cnt = cnt+1
                   !write(*,*) cnt,spr_nm,"cnt,spr"
                   curve_spr(cnt) = spr_nm
                   is_it_a_curve(spr_nm) = 1
                endif
                
             elseif (i.gt.(ncl+ncr) .and. i.le.(ncl+ncr+2)) then
                
                if (typ_nm .ne. 6) then
                   cnt = cnt+1
                   !write(*,*) cnt,spr_nm,"cnt,spr"
                   curve_spr(cnt) = spr_nm
                   is_it_a_curve(spr_nm) = 1
                endif
                
             elseif (i.gt.(ncl+ncr+2)) then
                
                if (typ_nm==4 .or. typ_nm==5) then
                   cnt = cnt+1
                   !write(*,*) cnt,spr_nm,"cnt,spr"
                   curve_spr(cnt) = spr_nm
                   is_it_a_curve(spr_nm) = 1
                endif
                
             endif
             
          enddo
       endif
    enddo
    
    !write(*,*) "cnt =",cnt,"and N_curve =",N_curve
    
    if (cnt.ne.N_curve) then
       write(*,*) "cnt must be equal to N_curve"
       write(*,*) "cnt =",cnt,"and N_curve =",N_curve
       stop
    endif
    
  end subroutine get_curve_spr_trnsfrm_and_is_it_a_curve_methd1
  

  
  
  subroutine get_curve_spr_trnsfrm_and_is_it_a_curve_methd2
    implicit none
    integer :: i,imax,j,jmax,jEnd
    integer :: spr_nm,typ_nm,cnt
    
    write(*,*) "NI_incldAS_woPP =",NI_incldAS_woPP
    
    if (NI_incldAS_woPP == 0) stop 'not Compatible for with Pulley System'
    if (NI_incldAS_woPP == 1) continue
    
    allocate(curve_spr(1:N_curve))
    allocate(is_it_a_curve(1:N_sprS))
    
    curve_spr     = -1
    is_it_a_curve = 0
    
    if (VF_regionModelled==0) imax = N_cell
    if (VF_regionModelled==1) imax = N_cell-1
    
    Hlf_Ncell = imax/2 ; write(*,*) Hlf_Ncell,N_sprS,"Hlf_Ncell and N_sprS"
    cnt       = 0
    jmax      = 3
    
    do i = 1,imax  
       
       if (i.ne.imax) jEnd = jmax
       if (i == imax) jEnd = N_sprS-((Hlf_Ncell*2)*3)
       
       write(*,*) "jEnd =",jEnd,"for cellNo =",i
          
       do j = 1,jEnd
          spr_nm                = (i-1)*jmax + j
          typ_nm                = typ_sprStr(spr_nm)
          cnt                   = cnt+1
          curve_spr(cnt)        = spr_nm
          is_it_a_curve(spr_nm) = 1
          write(*,*) spr_nm,typ_nm,cnt,spr_nm,"sprPrps"
       enddo
       
    enddo
    
    write(*,*) "cnt =",cnt,"and N_curve =",N_curve
    
    if (cnt.ne.N_curve) then
       write(*,*) "cnt must be equal to N_curve"
       write(*,*) "cnt =",cnt,"and N_curve =",N_curve
       stop
    endif
    
  end subroutine get_curve_spr_trnsfrm_and_is_it_a_curve_methd2
  
  
  
  
  
  subroutine shift_the_IntrmdNodes_to_help_equilibrium
    implicit none
    integer :: i
    
    do i = (N_nodeS+1),N_node
       
       if (mod(i,2)==0) then
          
          if (node_xy(i,1).lt.0.0d0) then
             node_xy(i,1) = node_xy(i,1) + x_shft
          elseif (node_xy(i,1).gt.0.0d0) then
             node_xy(i,1) = node_xy(i,1) - x_shft
          endif
          
       endif
       
    enddo
    
    
  end subroutine shift_the_IntrmdNodes_to_help_equilibrium
  
  
  
  subroutine deallocate_and_reallocate_spr_vars
    implicit none
    
    deallocate(typ_spr,k_spr,l0,l)
    
    call allocate_and_initialize_spring_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_spr_vars
  
  subroutine get_SprVars
    implicit none
    integer :: i,j,jmax,k
    integer :: spr_nm
    integer :: cnt,cnt_Pspr
    integer :: cnt_PrvSys
    integer :: spr_Curr
    
    call TN_and_NI_modelSpr_trnsfrmtn
    
    cnt_PrvSys = 1
    cnt  = 0
    jmax = 3
    
    spr_Curr = 0
    
    do i = 1,N_sprS
       !write(*,*) trmnl_intrmdSpr(i,0:(NAEC+1)),i,"trmnl_intrmd,spr"
       !if (i==(N_sprS-1) .or. i==N_sprS) write(*,*) typ_sprStr(i),"typ 49-50"
       
       if (trmnl_intrmdSpr(i,0)==1) then
          !write(*,*) "En1"
          typ_spr(spr_Curr+1)  = typ_sprStr(i)
          k_spr(spr_Curr+1)    = k_sprStr(i)
          l0(spr_Curr+1)       = l0_Str(i)
          l(spr_Curr+1)        = l_Str(i)
          
          !write(*,*)k_spr(spr_curr+1),l0(spr_curr+1),l(spr_curr+1),(spr_curr+1)
          
          spr_Curr = spr_Curr+1
          !write(*,*) i,Spr_Curr,"i,Spr_Curr"
          
       elseif (trmnl_intrmdSpr(i,0)==(NAEC+1)) then
          !write(*,*) "En2"
          
          do k=1,(NAEC+1)
             typ_spr(spr_Curr+k) = typ_sprStr(i)
             k_spr(spr_Curr+k)   = k_sprStr(i)*(NAEC+1)
             l0(spr_Curr+k)      = l0_Str(i)/(NAEC+1)
             l(spr_Curr+k)       = l_Str(i)/(NAEC+1)
          enddo
          
          !write(*,*)k_spr(spr_Curr+1),l0(spr_curr+1),l(spr_curr+1),(spr_curr+1)
          !write(*,*)k_spr(spr_Curr+2),l0(spr_curr+2),l(spr_curr+2),(spr_curr+2)
          !write(*,*) k_spr(spr_Curr+1),k_spr(spr_Curr+2),k_sprStr(i),"k_spr"
          !write(*,*) l0(spr_Curr+1),l0(spr_Curr+2),l0_Str(i),"l0"
          
          spr_Curr = spr_Curr+(NAEC+1) !2=NAEC+1
          !write(*,*) i,Spr_Curr,"i,Spr_Curr"
       endif
       
       !checking if spr_Curr=N_spr
       
       if (i==N_sprS) then        
          if (trmnl_intrmdSpr(i,0)==1) then
             !spr_Curr= spr_Curr-1
             continue
          elseif (trmnl_intrmdSpr(i,0)==(NAEC+1)) then
             !spr_Curr = spr_Curr-2 !2=(NAEC+1)
             continue
          endif
          
          if (spr_Curr.ne.N_spr) then
             write(*,*) "spr_Curr =",spr_Curr,"Nspr =",N_spr
             write(*,*) "spr_Curr is not equal to N_spr"
             stop
          endif
          
       endif
       
    enddo
    
    call print_sprVars_wo_Lt_alpha_optmSpr
    
  end subroutine get_SprVars
  
  
  subroutine get_SprVars_NIsaved
    implicit none
    
    call TN_and_NI_modelSpr_trnsfrmtn
    
    typ_spr(1:N_spr) = typ_sprNIsaved(1:N_spr)
    k_spr(1:N_spr)   = k_sprNIsaved(1:N_spr)
    l0(1:N_spr)      = l0_NIsaved(1:N_spr)
    l(1:N_spr)       = l_NIsaved(1:N_spr)
    
    call print_sprVars_wo_Lt_alpha_optmSpr
    
  end subroutine get_SprVars_NIsaved
  
  subroutine get_AreaVarsNIsaved
    implicit none
    
    k_area(1:N_cell) = k_areaNIsaved(1:N_cell)
    A0(1:N_cell)     = A0_NIsaved(1:N_cell)
    A(1:N_cell)      = A_NIsaved(1:N_cell)
    
  end subroutine get_AreaVarsNIsaved
  
  subroutine get_CgXvars
    implicit none
    
    CgXNode(1:N_node)  = 0.00d0 
    CgXNode(1:N_nodeS) = CgXNode_Str(1:N_nodeS)
    
  end subroutine get_CgXvars
  
  subroutine get_CgYvars
    implicit none
    
    CgYNode(1:N_node)  = 0.00d0 
    CgYNode(1:N_nodeS) = CgYNode_Str(1:N_nodeS)
    
  end subroutine get_CgYvars
  
  subroutine TN_and_NI_modelSpr_trnsfrmtn
    implicit none
    
    if (NI_incldAS_woPP == 0) call TN_and_NI_modelSpr_trnsfrmtn_methd1
    if (NI_incldAS_woPP == 1) call TN_and_NI_modelSpr_trnsfrmtn_methd2
    
  end subroutine TN_and_NI_modelSpr_trnsfrmtn
  
  subroutine TN_and_NI_modelSpr_trnsfrmtn_methd1
    implicit none
    integer :: i,j,k,jmax,jEnd
    integer :: spr_Prv,spr_Curr
    integer :: typ_nm
    
    allocate(trmnl_intrmdSpr(1:N_sprS,0:(NAEC+1)))
    
    jmax = 3
    trmnl_intrmdSpr = -1
    spr_Prv = 0 ; spr_Curr = 0
    
    do i = 1,N_cell
       
       if (i==N_cell) then
          
          if (stageNo==1 .and. stageType==1) then
             if (CellsMeet==0)   jEnd = 2
             if (CellsMeet.gt.0) jEnd = 1
          elseif (stageNo==4 .and. stageType==1) then
             jEnd = 1
          elseif (stageNo==4 .and. stageType==2) then
             jEnd = jmax
          endif
          
       elseif (i.ne.N_cell) then
          jEnd = jmax
       endif
       
       if (stageNo==1 .and. stageType==1) then
          
          do j = 1,jEnd
             
             spr_Prv = (i-1)*jmax + j
             typ_nm = typ_sprStr(spr_Prv)
             
             if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                
                trmnl_intrmdSpr(spr_Prv,0) = NAEC+1
                
                do k = 1,(NAEC+1)
                   trmnl_intrmdSpr(spr_Prv,k) = spr_Curr+k
                enddo
                
                spr_Curr = spr_Curr + (NAEC+1) !2=NAEC+1
                
             else
                trmnl_intrmdSpr(spr_Prv,0) = 1
                trmnl_intrmdSpr(spr_Prv,1) = spr_Curr+1
                spr_Curr = spr_Curr+1
             endif
             
          enddo
          
       elseif (stageNo==2 .or. stageNo==3) then
          write(*,*) "fl:node_insertion,sb: TI_and_NI model spr"
       elseif (stageNo==4) then
          
          do j = 1,jEnd
             spr_Prv = (i-1)*jmax + j
             typ_nm = typ_sprStr(spr_Prv)
             
             if (i.le.(ncl+ncr)) then
                
                if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                   
                   trmnl_intrmdSpr(spr_Prv,0) = NAEC+1
                   
                   do k = 1,(NAEC+1)
                      trmnl_intrmdSpr(spr_Prv,k) = spr_Curr+k
                   enddo
                   
                   spr_Curr = spr_Curr + (NAEC+1) !2=NAEC+1
                   
                else
                   trmnl_intrmdSpr(spr_Prv,0) = 1
                   trmnl_intrmdSpr(spr_Prv,1) = spr_Curr+1
                   spr_Curr = spr_Curr+1
                endif
                
             elseif ((i.gt.(ncl+ncr)) .and. (i.le.(ncl+ncr+2))) then
                
                if (typ_nm.ne.6) then
                   
                   trmnl_intrmdSpr(spr_Prv,0) = NAEC+1
                   
                   do k = 1,(NAEC+1)
                      trmnl_intrmdSpr(spr_Prv,k) = spr_Curr+k
                   enddo
                   
                   spr_Curr = spr_Curr + (NAEC+1) !2=NAEC+1
                   
                else
                   trmnl_intrmdSpr(spr_Prv,0) = 1
                   trmnl_intrmdSpr(spr_Prv,1) = spr_Curr+1
                   spr_Curr = spr_Curr+1
                   
                endif
                
             elseif ((i.gt.(ncl+ncr+2)) .and. (i.le.N_cell)) then
                
                if (typ_nm==4 .or. typ_nm==5) then
                   
                   trmnl_intrmdSpr(spr_Prv,0) = NAEC+1
                   
                   do k = 1,(NAEC+1)
                      trmnl_intrmdSpr(spr_Prv,k) = spr_Curr+k
                   enddo
                   
                   spr_Curr = spr_Curr + (NAEC+1) !2=NAEC+1
                   
                else
                   trmnl_intrmdSpr(spr_Prv,0) = 1
                   trmnl_intrmdSpr(spr_Prv,1) = spr_Curr+1
                   spr_Curr = spr_Curr+1
                endif
                
             endif
             
          enddo
          
       endif
       
    enddo
    
    
    open(unit=197,file='TrmnlIntrmdSpr1.dat')
    
    write(197,fmt=*) "       ","spr_nm","    ","trmnl_intrmdSpr"
    
    do i = 1,N_sprS
       write(197,fmt=*) i,trmnl_intrmdSpr(i,0:(NAEC+1))
    enddo
    
    close(197)
    
  end subroutine TN_and_NI_modelSpr_trnsfrmtn_methd1
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine TN_and_NI_modelSpr_trnsfrmtn_methd2
    implicit none
    integer :: i,j,k,jmax,jEnd
    integer :: spr_Prv,spr_Curr
    integer :: typ_nm,chk_val
    
    allocate(trmnl_intrmdSpr(1:N_sprS,0:(NAEC+1)))
    
    jmax            = 3
    trmnl_intrmdSpr = -1
    spr_Prv         = 0 ; spr_Curr = 0
    
    do i = 1,N_cell
       
       if (i==N_cell) then
          if (CellsMeet==0)   jEnd = 2
          if (CellsMeet.gt.0) jEnd = 1
       elseif (i.ne.N_cell) then
          jEnd = jmax
       endif
       
       do j = 1,jEnd
             
          spr_Prv = (i-1)*jmax + j
          typ_nm = typ_sprStr(spr_Prv)
          
          chk_val=0 ; call typ_spr_chk(typ_nm,chk_val)
          
          if (chk_val == 1) then
             trmnl_intrmdSpr(spr_Prv,0) = NAEC+1
             do k = 1,(NAEC+1)
                trmnl_intrmdSpr(spr_Prv,k) = spr_Curr+k
             enddo
             spr_Curr = spr_Curr + (NAEC+1) !2=NAEC+1
          else
             trmnl_intrmdSpr(spr_Prv,0) = 1
             trmnl_intrmdSpr(spr_Prv,1) = spr_Curr+1
             spr_Curr = spr_Curr+1
          endif
          
       enddo
          
    enddo
    
    
    open(unit=197,file='TrmnlIntrmdSpr2.dat')
    
    write(197,fmt=*) "       ","spr_nm","    ","trmnl_intrmdSpr"
    
    do i = 1,N_sprS
       write(197,fmt=*) i,trmnl_intrmdSpr(i,0:(NAEC+1))
    enddo
    
    close(197)
    
  end subroutine TN_and_NI_modelSpr_trnsfrmtn_methd2
  
  
  subroutine store_all_moving_coordnte_variables
    implicit none
    
    N_mvCoordnteS             = N_mvCoordnte
    N_mvCoordnte_withl0A0_S   = N_mvCoordnte_withl0A0
    
    coordntesStr              = coordntes
    coordntes_xyStr           = coordntes_xy
    
    coordnte_of_which_nodeStr = coordnte_of_which_node
    x_or_yStr                 = x_or_y 
    
  end subroutine store_all_moving_coordnte_variables
  
  subroutine store_all_moving_coordnte_variables_NI
    implicit none
    
    N_mvCoordnteSNI             = N_mvCoordnte
    N_mvCoordnte_withl0A0_SNI   = N_mvCoordnte_withl0A0
    
    coordntesStrNI              = coordntes
    coordntes_xyStrNI           = coordntes_xy
    
    coordnte_of_which_nodeStrNI = coordnte_of_which_node
    x_or_yStrNI                 = x_or_y 
    
  end subroutine store_all_moving_coordnte_variables_NI
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !transfrm arrays edit
  
  subroutine store_TN_trnsfrms
    implicit none
    
    max_node_sprS  = max_node_spr
    max_area_sprS  = max_area_spr
    max_spr_nodeS  = max_spr_node
    max_area_nodeS = max_area_node
    max_node_areaS = max_node_area
    max_spr_areaS  = max_spr_area
    
    node_sprS  = node_spr
    node_areaS = node_area   
    spr_nodeS  = spr_node
    spr_areaS  = spr_area
    area_nodeS = area_node
    area_sprS  = area_spr
    
    NlistS = Nlist
    
    call print_all_stored_trnsfrms
    
  end subroutine store_TN_trnsfrms
  
  subroutine store_anyTypPrv_trnsfrms
    implicit none
    
    max_node_sprS  = max_node_spr
    max_area_sprS  = max_area_spr
    max_spr_nodeS  = max_spr_node
    max_area_nodeS = max_area_node
    max_node_areaS = max_node_area
    max_spr_areaS  = max_spr_area
    
    call resolve_allocataion_problems_with_trnsfrm_arrays
    
    node_sprS  = node_spr
    node_areaS = node_area   
    spr_nodeS  = spr_node
    spr_areaS  = spr_area
    area_nodeS = area_node
    area_sprS  = area_spr
    
    NlistS = Nlist
    
    call print_all_stored_trnsfrms_v2
    
  end subroutine store_anyTypPrv_trnsfrms
  
  subroutine resolve_allocataion_problems_with_trnsfrm_arrays
    implicit none
    
    if (.not.allocated(node_sprS)) then
       allocate(node_sprS(1:N_nodeS,0:max_spr_nodeS)) ; node_sprS = -1
    else
       deallocate(node_sprS)
       allocate(node_sprS(1:N_nodeS,0:max_spr_nodeS)) ; node_sprS = -1
    endif
    
    if (.not.allocated(node_areaS)) then
       allocate(node_areaS(1:N_nodeS,0:max_area_nodeS)) ; node_areaS = -1
    else
       deallocate(node_areaS)
       allocate(node_areaS(1:N_nodeS,0:max_area_nodeS)) ; node_areaS = -1
    endif
    
    
    if (.not.allocated(spr_nodeS)) then
       allocate(spr_nodeS(1:N_sprS,0:max_node_sprS)) ; spr_nodeS = -1
    else
       deallocate(spr_nodeS)
       allocate(spr_nodeS(1:N_sprS,0:max_node_sprS)) ; spr_nodeS = -1
    endif
    
    
    if (.not.allocated(spr_areaS)) then
       allocate(spr_areaS(1:N_sprS,0:max_area_sprS)) ; spr_areaS = -1
    else
       deallocate(spr_areaS)
       allocate(spr_areaS(1:N_sprS,0:max_area_sprS)) ; spr_areaS = -1
    endif
    
    
    if (.not.allocated(area_nodeS)) then
       allocate(area_nodeS(1:N_cellS,0:max_node_areaS)) ; area_nodeS = -1
    else
       deallocate(area_nodeS)
       allocate(area_nodeS(1:N_cellS,0:max_node_areaS)) ; area_nodeS = -1
    endif
    
    
    if (.not.allocated(area_sprS)) then
       allocate(area_sprS(1:N_cellS,0:max_spr_areaS)) ; area_sprS = -1
    else
       deallocate(area_sprS)
       allocate(area_sprS(1:N_cellS,0:max_spr_areaS)) ; area_sprS = -1
    endif
    
    
    if (.not.allocated(NlistS)) then
       allocate(NlistS(1:N_nodeS,1:max_area_nodeS,1:2)) ; NlistS = -1
    else
       deallocate(NlistS)
       allocate(NlistS(1:N_nodeS,1:max_area_nodeS,1:2)) ; NlistS = -1
    endif
    
    
  end subroutine resolve_allocataion_problems_with_trnsfrm_arrays
  
  subroutine store_NI_trnsfrms_forContinuedTrnsfrm
    implicit none
    
    max_node_sprS2  = max_node_spr
    max_area_sprS2  = max_area_spr
    max_spr_nodeS2  = max_spr_node
    max_area_nodeS2 = max_area_node
    max_node_areaS2 = max_node_area
    max_spr_areaS2  = max_spr_area
    
    node_sprS2  = node_spr
    node_areaS2 = node_area   
    spr_nodeS2  = spr_node
    spr_areaS2  = spr_area
    area_nodeS2 = area_node
    area_sprS2  = area_spr
    
    NlistS2 = Nlist
    
    call print_all_stored_trnsfrms_S2
    
  end subroutine store_NI_trnsfrms_forContinuedTrnsfrm
  
  
  subroutine store_NI_trnsfrms
    implicit none
    
    max_node_sprSNI  = max_node_spr
    max_area_sprSNI  = max_area_spr
    max_spr_nodeSNI  = max_spr_node
    max_area_nodeSNI = max_area_node
    max_node_areaSNI = max_node_area
    max_spr_areaSNI  = max_spr_area
    
    node_sprSNI  = node_spr
    node_areaSNI = node_area   
    spr_nodeSNI  = spr_node
    spr_areaSNI  = spr_area
    area_nodeSNI = area_node
    area_sprSNI  = area_spr
    
    NlistSNI = Nlist
    
    call print_all_stored_trnsfrms
    
  end subroutine store_NI_trnsfrms
  
  
  subroutine deallocate_and_reallocate_for_NI_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)

    if (SystemTyp .ne. 1) then
       max_node_spr  = 3 
       max_area_spr  = 2 
       max_spr_node  = 6
       max_area_node = 4 
       max_node_area = 5+(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+NAEC_Ltrl) ! 5+(4 sides)-->MaxNo_curve_in_Each_cell*NAEC
       max_spr_area  = 4+(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+NAEC_Ltrl) ! 4+(4 sides)-->MaxNo_curve_in_Each_cell*NAEC
       
    elseif (SystemTyp == 1) then
       max_node_spr  = 4 ! different than prev blck
       max_area_spr  = 2 
       max_spr_node  = 6
       max_area_node = 4 
       max_node_area = 6+(3*NAEC) ! 5+(3*1) --> MaxNo_curve_in_Each_cell*NAEC !diff from prev blck
       max_spr_area  = 4+(3*NAEC) ! 4+(3*1) --> MaxNo_curve_in_Each_cell*NAEC
    endif
    
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
    
  end subroutine deallocate_and_reallocate_for_NI_trnsfrmVars
  
  
  subroutine deallocate_and_reallocate_for_TN_trnsfrmVars
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    max_node_spr  = 3 
    max_area_spr  = 2 
    max_spr_node  = 6
    max_area_node = 4 
    max_node_area = 5
    max_spr_area  = 4
    
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
    
  end subroutine deallocate_and_reallocate_for_TN_trnsfrmVars
  
  
  subroutine store_all_gradient_variables
    implicit none
    
    call store_gradient_system_variables
    call store_gradient_spr_variables
    call store_gradient_area_variables
    call store_gradient_grvtnl_variables
    call store_gradient_bend_variables
    call store_gradient_SRyp_variables
    
  end subroutine store_all_gradient_variables
  
  subroutine store_all_gradient_variables_NIcontinued
    implicit none
    
    call store_gradient_system_variables_S2
    call store_gradient_spr_variables_S2
    call store_gradient_area_variables_S2
    call store_gradient_grvtnl_variables_S2
    call store_gradient_bend_variables_S2
    
  end subroutine store_all_gradient_variables_NIcontinued
  
  
  subroutine store_all_gradient_variables_NI
    implicit none
    
    call store_gradient_system_variables_NI
    call store_gradient_spr_variables_NI
    call store_gradient_area_variables_NI
    call store_gradient_grvtnl_variables_NI
    call store_gradient_bend_variables_NI
    
  end subroutine store_all_gradient_variables_NI
  
  subroutine switchback_to_TN_model
    implicit none
    
    modelID = 1
    write(*,*) modelID,"modelID"
    !stop
    
    call restore_TN_vars
    call deallocate_and_reallocate_TN_Arrays
    call restore_TN_arrays
    
    call dealloc_and_restore_TN_moving_variables
    
    call restore_TN_trnsfrm_vars
    call dealloc_and_restore_TN_trnsfrm_Arrays
    
    call deallocate_and_restore_TN_gradient_variables
    
  end subroutine switchback_to_TN_model
  
  subroutine restore_TN_vars
    implicit none
    
    N_node = N_nodeS
    N_spr  = N_sprS
    N_phi  = N_phiS
    
  end subroutine restore_TN_vars
  
  subroutine restore_TN_arrays
    implicit none
    
    node_xy       = node_xyStr
    node_typ      = node_typStr   
    node_cnnctd   = node_cnnctdStr 
    double_node   = double_nodeStr
    count_this_dn = count_this_dnStr
    
    typ_spr = typ_sprStr
    k_spr   = k_sprStr   
    l0      = l0_Str 
    l       = l_Str
    
    k_phi   = k_phiStr
    CgXNode = CgXNode_Str
    CgYNode = CgYNode_Str
    
  end subroutine restore_TN_arrays
  
  subroutine dealloc_and_restore_TN_moving_variables
    implicit none
    
    call deallocate_moving_coordnte_variables_wo_StrVars
    
    N_mvCoordnte          = N_mvCoordnteS
    N_mvCoordnte_withl0A0 = N_mvCoordnte_withl0A0_S

    call alloc_and_init_moving_coordnte_variables_wo_StrVars
    call restore_all_moving_coordnte_variables
    
  end subroutine dealloc_and_restore_TN_moving_variables
  
  subroutine restore_all_moving_coordnte_variables
    implicit none
    
    coordntes    = coordntesStr
    coordntes_xy = coordntes_xyStr
    
    coordnte_of_which_node = coordnte_of_which_nodeStr
    x_or_y                 = x_or_yStr 
    
  end subroutine restore_all_moving_coordnte_variables

  
  subroutine restore_TN_trnsfrm_vars
    implicit none
    
    max_node_spr  = max_node_sprS
    max_area_spr  = max_area_sprS
    
    max_spr_node  = max_spr_nodeS
    max_area_node = max_area_nodeS
    
    max_node_area = max_node_areaS
    max_spr_area  = max_spr_areaS
    
  end subroutine restore_TN_trnsfrm_vars
  
  subroutine dealloc_and_restore_TN_trnsfrm_Arrays
    implicit none
    
    deallocate(node_spr,node_area)
    deallocate(spr_node,spr_area)
    deallocate(area_node,area_spr)
    deallocate(Nlist)
    
    allocate(spr_node(1:N_spr,0:max_node_spr))
    allocate(spr_area(1:N_spr,0:max_area_spr))
    allocate(node_spr(1:N_node,0:max_spr_node))
    allocate(node_area(1:N_node,0:max_area_node))
    allocate(area_spr(1:N_cell,0:max_spr_area))
    allocate(area_node(1:N_cell,0:max_node_area))
    allocate(Nlist(1:N_node,1:max_area_node,1:2))
    
    node_spr  = node_sprS
    node_area = node_areaS    
    spr_node  = spr_nodeS
    spr_area  = spr_areaS
    area_node = area_nodeS
    area_spr  = area_sprS
    
    Nlist = NlistS
    
  end subroutine dealloc_and_restore_TN_trnsfrm_Arrays

  subroutine deallocate_and_restore_TN_gradient_variables
    implicit none
    
    call deallocate_all_gradient_variables_wo_StrVars
    call restore_all_gradient_variables
    
  end subroutine deallocate_and_restore_TN_gradient_variables
  
  subroutine restore_all_gradient_variables
    implicit none
    
    call restore_gradient_system_variables
    call restore_gradient_spr_variables
    call restore_gradient_area_variables
    call restore_gradient_grvtnl_variables
    call restore_gradient_bend_variables

  end subroutine restore_all_gradient_variables
    
  subroutine deallocate_repetitive_arrays
    implicit none

    deallocate(curve_spr)
    deallocate(is_it_a_curve)
    deallocate(trmnl_intrmdSpr)
    deallocate(insrtd_vrtxNode)
    
  end subroutine deallocate_repetitive_arrays
  
  ! subroutine check_curve_orientation
  !   implicit none
  !   integer :: spr_nm
  !   real*8  :: cm1,cm2
    
  !   do i = 1,N_curve

  !      node_nm = N_nodeS + ()
  !      spr_nm  = curve_spr(i)
       
  !      if (typ_sprS(spr_nm)==5) then
  !         cm1 = node_xyBfr(node_nm,1)
  !         cm2 = node_xy(node_nm,1)

  !         if (abs(cm1).lt.abs(cm2)) then
  !            call shorten_RelatdSprs_and_equilibrate(spr_nm)
  !         endif
          
  !      endif
       
  !   enddo
    
  ! end subroutine check_curve_orientation
  

   subroutine get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    implicit none
    
    real*8,  intent(in)    :: tol_Rdfn
    logical, intent(inout) :: lgcl_rdfn
    
    logical :: lgcl_NeighSwitch
    integer :: Pulley,LftNeigh,RghtNeigh
    integer :: LastCellAboveLft,LastCellAboveRght
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    lgcl_rdfn = .False.
    
    if (stageNo==1 .and. stageType==1) then
       
       open(unit=96,file='OriginDis.dat',position='append')
       
       if (modelID==1) then
          call get_LftRghtNeigh_ofOrigin(LftNeigh,RghtNeigh)
          write(*,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh"
          
       elseif (modelID==2) then
          call find_LastCellofTopLayer(LastCellAboveLft)
          LastCellAboveRght = Hlf_Ncell + LastCellAboveLft
          
          LftNeigh  = 2*(LastCellAboveLft+1)  - 1
          RghtNeigh = 2*(LastCellAboveRght+2) - 1
          write(*,*) LftNeigh,RghtNeigh,"LftNeigh,RghtNeigh NI"
          write(*,*) origin(1:2),"origin for NI"
       endif
       
          
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
          
       elseif (Origin_dis(1).gt.tol_Rdfn.and.Origin_dis(2).gt.tol_Rdfn) then
          write(*,*) "Both are big"
          
       elseif (lgcl_NeighSwitch .eqv. .True.) then
          lgcl_rdfn = .True.
          write(*,*) "Neighbour of Origin switched"
       endif
       
    elseif (stageNo==4) then
       
       open(unit=98,file='PulleyDis.dat',position='append')
       call get_Pulley_and_Neighbours(Pulley,LftNeigh,RghtNeigh)
       write(*,*) Pulley,LftNeigh,RghtNeigh,"Pulley,LftNeigh,RghtNeigh"
       
       Pulley_dis(1:2) = dist_to_Pulley(Pulley,LftNeigh,RghtNeigh)
       write(*,*) Pulley_dis(1:2),"Pulley dis"
       write(98,fmt=*) Pulley_dis(1:2)
       close(98)
       
       lgcl_NeighSwitch = .False.
       call Pulley_Neighbour_switch(Pulley,LftNeigh,RghtNeigh,lgcl_NeighSwitch)
       
       
       if (Pulley_dis(1).lt.tol_Rdfn .and.Pulley_dis(2).lt.tol_Rdfn) then
          write(*,*) "Both are less"
          lgcl_rdfn = .True.
          
       elseif (Pulley_dis(1).lt.tol_Rdfn .and. Pulley_dis(2).gt.tol_Rdfn) then
          write(*,*) "Lft is less"
          
       elseif (Pulley_dis(1).gt.tol_Rdfn .and. Pulley_dis(2).lt.tol_Rdfn) then
          write(*,*) "Rght is less"
          
       elseif (Pulley_dis(1).gt.tol_Rdfn .and. Pulley_dis(2).gt.tol_Rdfn) then
          write(*,*) "Both are big"
          
       elseif (lgcl_NeighSwitch .eqv. .True.) then
          lgcl_rdfn = .True.
          write(*,*) "Neighbour of Pulley switched"
       endif
       
    endif
    
    
  end subroutine get_decision_of_redefining
  
  
  subroutine get_decision_of_redefining_diffWay(tol_Rdfn,lgcl_rdfn)
    implicit none
    
    real*8,  intent(in)    :: tol_Rdfn
    logical, intent(inout) :: lgcl_rdfn
    
    logical :: lgcl_NeighSwitch
    integer :: LftNeigh,RghtNeigh
    
    LftNeigh  = (Hlf_Ncell+1)*2 - 1 - (CellsMeet)*2
    RghtNeigh = LftNeigh + (Hlf_Ncell+1)*2
    write(*,*) LftNeigh,RghtNeigh,"Lft-Rght Neigh"
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    
    lgcl_rdfn = .False.
    
    Origin_dis(1:2)  = dist_to_Origin(LftNeigh,RghtNeigh)  ; write(*,*) Origin_dis(1:2),"Origin dis Pos1"
    lgcl_NeighSwitch = .False. ; call Origin_Neighbour_switch(LftNeigh,RghtNeigh,lgcl_NeighSwitch)
       
    if ((Origin_dis(1).le.tol_Rdfn).and.(Origin_dis(2).le.tol_Rdfn)) then
       write(*,*) "Both are less" ; lgcl_rdfn = .True.
       
    elseif ((Origin_dis(1).le.tol_Rdfn).and.(Origin_dis(2).gt.tol_Rdfn)) then
       write(*,*) "Lft is less"
       
    elseif ((Origin_dis(1).gt.tol_Rdfn).and.(Origin_dis(2).le.tol_Rdfn)) then
       write(*,*) "Rght is less"
       
    elseif ((Origin_dis(1).gt.tol_Rdfn).and.(Origin_dis(2).gt.tol_Rdfn)) then
       write(*,*) "Both are big"
       
    elseif (lgcl_NeighSwitch .eqv. .True.) then
       lgcl_rdfn = .True. ; write(*,*) "Neighbour of Origin switched"
    endif
    
    
  end subroutine get_decision_of_redefining_diffWay
  
  
  subroutine get_decision_of_tensionPositivity(NAECval,sprL,sprR,tol_Tnsn,lgcl_Tnsn,TnsnVal)
    implicit none
    integer, intent(in)    :: NAECval
    integer, intent(in)    :: sprL(1:(NAECval+1)),sprR(1:(NAECval+1))
    real*8,  intent(in)    :: tol_Tnsn
    logical, intent(inout) :: lgcl_Tnsn
    real*8,  intent(in)    :: TnsnVal(1:N_spr)
    
    integer :: lftSpr,rghtSpr
    real*8  :: TnsnLft,TnsnRght
    integer :: i,j
    
    lgcl_Tnsn = .False.
    
    open(unit=76,file='TnsnPositivity.dat',position='append')
    
    write(76,*) NAECVal,"NAEC"
    write(76,*) sprL(1:(NAECval+1)),"sprL"
    write(76,*) sprR(1:(NAECval+1)),"sprR"
    
    do i = 1,(NAECval+1)
       
       lftSpr   = sprL(i) ; rghtSpr = sprR(i)
       TnsnLft  = TnsnVal(lftSpr)
       TnsnRght = TnsnVal(rghtSpr)
       
       write(76,*) lftSpr,rghtSpr,TnsnLft,TnsnRght,"sprs & Tnsn"
       
       if ( (abs(TnsnLft) .lt. tol_Tnsn) .or. (abs(TnsnLft) .lt. tol_Tnsn) ) then
          lgcl_Tnsn = .True.
          write(76,*) lftSpr,rghtSpr,TnsnLft,TnsnRght,"spr & Tnsn"
          exit
       endif
       
    enddo

    close(76)
    
  end subroutine get_decision_of_tensionPositivity
    
  function dist_to_Pulley(Pulley,LftNeigh,RghtNeigh)
    
    implicit none
    real*8  :: dist_to_Pulley(1:2)
    integer, intent(in) :: Pulley,LftNeigh,RghtNeigh

    real*8 ::  Pulley_xy(1:2),LftNeigh_xy(1:2),RghtNeigh_xy(1:2)
    
    interface
       
       real*8 function distance(Point1,Point2,N_dim)
         implicit none
         integer, intent(in) :: N_dim
         real*8,  intent(in) :: Point1(1:N_dim)
         real*8,  intent(in) :: Point2(1:N_dim)
       end function distance
       
    end interface
    
    Pulley_xy(1:2)     = node_xy(Pulley,1:2)
    LftNeigh_xy(1:2)   = node_xy(LftNeigh,1:2)
    RghtNeigh_xy(1:2)  = node_xy(RghtNeigh,1:2)
    
    dist_to_Pulley(1) = distance(Pulley_xy,LftNeigh_xy,N_dmnsn)
    dist_to_Pulley(2) = distance(Pulley_xy,RghtNeigh_xy,N_dmnsn)
    
  end function dist_to_Pulley
  
  function dist_to_Origin(LftNeigh,RghtNeigh)
    
    implicit none
    real*8              :: dist_to_Origin(1:2)
    integer, intent(in) :: LftNeigh,RghtNeigh
    real*8              ::  Origin_xy(1:2),LftNeigh_xy(1:2),RghtNeigh_xy(1:2)
    
    interface
       
       real*8 function distance(Point1,Point2,N_dim)
         implicit none
         integer, intent(in) :: N_dim
         real*8,  intent(in) :: Point1(1:N_dim)
         real*8,  intent(in) :: Point2(1:N_dim)
       end function distance
       
    end interface
    
    Origin_xy(1:2)     = origin(1:2)
    LftNeigh_xy(1:2)   = node_xy(LftNeigh,1:2)
    RghtNeigh_xy(1:2)  = node_xy(RghtNeigh,1:2)
    
    dist_to_Origin(1) = distance(Origin_xy,LftNeigh_xy,N_dmnsn)
    dist_to_Origin(2) = distance(Origin_xy,RghtNeigh_xy,N_dmnsn)
    
  end function dist_to_Origin

  subroutine find_LastCellofTopLayer(LastCellAbove)
    implicit none
    integer, intent(out) :: LastCelLAbove
    real*8  :: yDis,ZERO
    integer :: i
    
    open(unit=54,file='LastCellTopLayer.dat')
    
    ZERO = 0.00d0
    
    
    do i = 1,(N_nodeS/2),2
       write(54,*) i,"i"
       yDis = node_xy(i,2)
       
       if (yDis.lt.ZERO) then
          LastCellAbove = (i-2)/2
          write(54,*) LastCellAbove,"LastCellAbove"
          exit
       endif
       
    enddo
    
    close(54)
    
  end subroutine find_LastCellofTopLayer
  
  subroutine find_FrstCellofBottomEndToAdjust(FrstCellAdjstBottm)
    implicit none
    integer, intent(out) :: FrstCellAdjstBottm
    
    open(unit=54,file='FrstCellAdjstBottm.dat')
    
    FrstCellAdjstBottm = Hlf_Ncell-1

    write(54,*) FrstCellAdjstBottm,"FrstCellAdjstBotm"
    
    close(54)
    
  end subroutine find_FrstCellofBottomEndToAdjust
    
  subroutine ApprchnCellSprShortening_and_Equilibrate_UntilMetNI(cellL,cellR,hmuch,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellL,celLR
    real*8,  intent(inout) :: hmuch
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: TINYval
    logical :: lgcl_Rdfn,lgcl_Tnsn
    real*8  :: tol_Rdfn,tol_Tnsn
    
    integer :: sprCntEachCell
    integer :: sprCntApcl,sprCntBsal,sprCntLtrl
    integer :: Nstp,sprPair(1:2)
    real*8  :: hm
    integer :: i,j,jmax,m
    integer :: NAECval
    
    integer, allocatable   :: APsprL(:),APsprR(:)
    integer, allocatable   :: BSsprL(:),BSsprR(:)
    integer, allocatable   :: LTsprL(:),LTsprR(:)
    
    integer, allocatable   :: spAP(:,:),spBS(:,:),spLT(:,:)
    real*8 , allocatable   :: Initl0AP(:,:),Initl0BS(:,:),Initl0LT(:,:)
    
    real*8 :: TnsnVal(1:N_spr)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    open(unit=15,file='DiagShortenUntilMet.dat')
    
    TINYval = 1.0d-15
    
    !NAEC_Apcl = 0 ; NAEC_Bsal = NAEC ; NAEC_Ltrl = NAEC
    
    allocate(APsprL(1:(NAEC_Apcl+1)),APsprR(1:(NAEC_Apcl+1)))
    allocate(BSsprL(1:(NAEC_Bsal+1)),BSsprR(1:(NAEC_Bsal+1)))
    allocate(LTsprL(1:(NAEC_Ltrl+1)),LTsprR(1:(NAEC_Ltrl+1)))
    
    allocate(spAP(1:(NAEC_Apcl+1),1:2))
    allocate(spBS(1:(NAEC_Bsal+1),1:2))
    allocate(spLT(1:(NAEC_Ltrl+1),1:2))
    
    allocate(Initl0AP(1:(NAEC_Apcl+1),1:2))
    allocate(Initl0BS(1:(NAEC_Bsal+1),1:2))
    allocate(Initl0LT(1:(NAEC_Ltrl+1),1:2))
    
    sprCntApcl = NAEC_Apcl+1    
    sprCntBsal = NAEC_Bsal+1
    sprCntLtrl = NAEC_Ltrl+1
    
    sprCntEachCell = sprCntApcl + sprCntBsal + sprCntLtrl
    
    do i = 1,3
       
       if (i==1) jmax = (NAEC_Apcl+1)
       if (i==2) jmax = (NAEC_Bsal+1)
       if (i==3) jmax = (NAEC_Ltrl+1)
       
       do j = 1,jmax
          
          if (i==1) then
             APsprL(j) = (sprCntEachCell*(cellL-1)) + j
             APsprR(j) = (sprCntEachCell*(cellR-1)) + j
             write(15,*) APsprL(j),APsprL(j),j,"Apcl Sprs"
             
          elseif (i==2) then
             BSsprL(j) = (sprCntEachCell*(cellL-1)) + (sprCntApcl) + j
             BSsprR(j) = (sprCntEachCell*(cellR-1)) + (sprCntApcl) + j
             write(15,*) BSsprL(j),BSsprL(j),j,"Bsal Sprs"
             
          elseif (i==3) then
             LTsprL(j) = (sprCntEachCell*(cellL-1)) + (sprCntApcl+sprCntBsal) + j
             LTsprR(j) = (sprCntEachCell*(cellR-1)) + (sprCntApcl+sprCntBsal) + j
             write(15,*) LTsprL(j),LTsprL(j),j,"Ltrl Sprs"
          endif
             
       enddo
    enddo
    
    do i = 1,3
       
       if (i==1) jmax = (NAEC_Apcl+1)
       if (i==2) jmax = (NAEC_Bsal+1)
       if (i==3) jmax = (NAEC_Ltrl+1)
       
       do j = 1,jmax
          
          if (i==1) then
             spAP(j,1)     = APsprL(j)     ; spAP(j,2)     = APsprR(j)
             Initl0AP(j,1) = l0(spAP(j,1)) ; Initl0AP(j,2) = l0(spAP(j,2))
             write(15,*) spAP(j,1:2),j,"spAP and j"
             write(15,*) Initl0AP(j,1:2),j,"Initl0AP and j"
             
          elseif (i==2) then
             spBS(j,1)     = BSsprL(j)     ; spBS(j,2)     = BSsprR(j)
             Initl0BS(j,1) = l0(spBS(j,1)) ; Initl0BS(j,2) = l0(spBS(j,2))
             write(15,*) spBS(j,1:2),j,"spBS and j"
             write(15,*) Initl0BS(j,1:2),j,"Initl0BS and j"
             
          elseif (i==3) then
             spLT(j,1)     = LTsprL(j)     ; spLT(j,2)     = LTsprR(j)
             Initl0LT(j,1) = l0(spLT(j,1)) ; Initl0LT(j,2) = l0(spLT(j,2))
             write(15,*) spLT(j,1:2),j,"spLT and j"
             write(15,*) Initl0LT(j,1:2),j,"Initl0LT and j"
          endif
          
       enddo
       
    enddo
    
    ! i=1
    
    ! do
    !    open(unit=17,file='ChkTnsnLT.dat',position='append')
       
       
    !    lgcl_Rdfn = .False. ; tol_Rdfn = 0.05d0
    !    call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn)
       
    !    lgcl_Tnsn=.False. ; tol_Tnsn=0.02d0 ; NAECval=NAEC_Ltrl
    !    TnsnVal(1:N_spr) = TnsnComprsn(node_xy,l0)
       
    !    do m=1,N_spr
    !       write(17,*) k_spr(m),l0(m),TnsnVal(m),m,"bfr"
    !    enddo
    !    write(17,*) " "
       
    !    call get_decision_of_tensionPositivity(NAECval,LTsprL,LTsprR,tol_Tnsn,lgcl_Tnsn,TnsnVal)

    !    do m=1,N_spr
    !       write(17,*) k_spr(m),l0(m),TnsnVal(m),m,"aft"
    !    enddo

    !    write(17,*) " "
       
    !    do m=1,N_node
    !       write(17,*) node_xy(m,1:2),m,"nodes"
    !    enddo
       
    !    close(17)
       
    !    if (lgcl_Rdfn .eqv. .True.) then
    !       hmuch = hm
    !       write(15,*) (i-1),"steps needed"
    !       write(15,*) (hmuch)
    !       exit
          
    !    elseif (lgcl_Tnsn .eqv. .True.) then
    !       hmuch = hm
    !       write(15,*) (i-1),"steps needed Negative"
    !       write(15,*) (hmuch)
    !       exit
          
    !    else
          
    !       do j = 1,(NAEC_Ltrl+1)
    !          l0(spLT(j,1)) = Initl0LT(j,1)
    !          l0(spLT(j,2)) = Initl0LT(j,2)
    !       enddo
          
    !       do j = 1,(NAEC_Ltrl+1)
    !          sprPair(1:2) = spLT(j,1:2)
    !          hm = (i)*(hmuch)
    !          call shortening_sprPair(sprPair,hm)
    !          write(15,*) l0(sprPair(1)),l0(sprPair(2)),j,"l0 AftShrtn insd main sb"
    !       enddo
          
    !       call Equilibrate_system
    !       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    !       FrmNo = FrmNo + 1
    !       write(15,*) (FrmNo-1),"FrmNo"
          
    !       i=i+1
          
    !    endif
       
    ! enddo
    
    ! write(15,*) ExpNo,(FrmNo-1),"ExpNo,FrmNo aft Ltrl"
    ! write(15,*) " "
    
    lgcl_Rdfn=.False. !this line wont be here
    lgcl_Tnsn=.True.  !and this line
    
    if ((lgcl_Rdfn.eqv..False.) .and. (lgcl_Tnsn.eqv..True.)) then
       
       i=1
       
       do
          open(unit=18,file='ChkTnsnAP.dat',position='append')
          
          lgcl_Rdfn = .False. ; tol_Rdfn = 0.05d0
          call get_decision_of_redefining(tol_Rdfn,lgcl_Rdfn)
          
          do m=1,N_spr
             write(18,*) k_spr(m),l0(m),TnsnVal(m),m,"bfr"
          enddo
          write(18,*) " "
          
          lgcl_Tnsn=.False. ; tol_Tnsn=0.02d0 ; NAECval=NAEC_Apcl
          call get_decision_of_tensionPositivity(NAECval,APsprL,APsprR,tol_Tnsn,lgcl_Tnsn,TnsnVal)
          
          do m=1,N_spr
             write(18,*) k_spr(m),l0(m),TnsnVal(m),m,"aft"
          enddo
          
          write(18,*) " "
          
          do m=1,N_node
             write(18,*) node_xy(m,1:2),m,"nodes"
          enddo
          
          close(18)
          
          if (lgcl_Rdfn .eqv. .True.) then
             hmuch = hm
             write(15,*) (i-1),"steps needed"
             write(15,*) (hmuch)
             exit
             
          elseif (lgcl_Tnsn .eqv. .True.) then
             hmuch = hm
             write(15,*) (i-1),"steps needed Negative"
             write(15,*) (hmuch)
             exit
             
          else
             
             do j = 1,(NAEC_Apcl+1)
                l0(spAP(j,1)) = Initl0AP(j,1)
                l0(spAP(j,2)) = Initl0AP(j,2)
             enddo
             
             do j = 1,(NAEC_Apcl+1)
                sprPair(1:2) = spAP(j,1:2)
                hm = (i)*(hmuch)
                call shortening_sprPair(sprPair,hm)
                write(15,*) l0(sprPair(1)),l0(sprPair(2)),j,"l0 AftShrtnInsd main sb"
             enddo
             
             call Equilibrate_system
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
             FrmNo = FrmNo + 1
             
             i=i+1
             
          endif
          
       enddo
       
    endif
    
    if (lgcl_Tnsn .eqv. .True.) then
       write(*,*) "Still not meeting"
    endif
    
    write(15,*) ExpNo,(FrmNo-1),"ExpNo,FrmNo aft Apcl"
    write(15,*) " "
    close(15)
    
  end subroutine ApprchnCellSprShortening_and_Equilibrate_UntilMetNI
  
  subroutine lgcls_and_weighted_arrays(cellL,cellR,sprLgcl,cellLgcl,weightSpr,weightCell)
    implicit none
    integer, intent(in)  :: cellL,cellR
    integer, intent(out) :: sprLgcl(1:N_spr),cellLgcl(1:N_cell)
    real*8 , intent(out) :: weightSpr(1:N_spr),weightCell(1:N_cell)
    
    integer :: sprCntEachCell
    integer :: apSprL(1:(NAEC_Apcl+1)),apSprR(1:(NAEC_Apcl+1))
    integer :: bsSprL(1:(NAEC_Bsal+1)),bsSprR(1:(NAEC_Bsal+1))
    integer :: ltSpr1L(1:(NAEC_Ltrl+1)),ltSpr1R(1:(NAEC_Ltrl+1))
    integer :: ltSpr2L(1:(NAEC_Ltrl+1)),ltSpr2R(1:(NAEC_Ltrl+1))
    integer :: sprL,sprR
    real*8  :: weightSprVal,weightCellVal
    integer :: i,j,jmax
    
    sprCntEachCell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    weightSprVal  = 0.40d0
    !weightCellVal = 1.00d0
    
    cellLgcl(cellL) = 0 ; cellLgcl(cellR) = 0
    
    do i = 1,4
       
       if (i==1) jmax = NAEC_Ltrl+1
       if (i==2) jmax = NAEC_Apcl+1
       if (i==3) jmax = NAEC_Bsal+1
       if (i==4) jmax = NAEC_Ltrl+1 
       
       do j = 1,jmax
          
          if (i==1) then
             ltSpr1L(j) = (cellL-2)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             ltSpr1R(j) = (cellR-2)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             
             sprL = ltSpr1L(j) ; sprR = ltSpr1R(j)
             sprLgcl(sprL) = 1 ; sprLgcl(sprR) = 1
             
             weightSpr(sprL) = weightSprVal
             weightSpr(sprR) = weightSprVal
             
          elseif (i==2) then
             apSprL(j) = (cellL-1)*(sprCntEachCell) + j
             apSprR(j) = (cellR-1)*(sprCntEachCell) + j
             
             sprL = apSprL(j)  ; sprR = apSprR(j)
             sprLgcl(sprL) = 1 ; sprLgcl(sprR) = 1
             
             weightSpr(sprL) = weightSprVal
             weightSpr(sprR) = weightSprVal
             
          elseif (i==3) then
             bsSprL(j) = (cellL-1)*(sprCntEachCell) + (NAEC_Apcl+1) + j
             bsSprR(j) = (cellR-1)*(sprCntEachCell) + (NAEC_Apcl+1) + j
             
             sprL = bsSprL(j)  ; sprR = bsSprR(j)
             sprLgcl(sprL) = 1 ; sprLgcl(sprR) = 1
             
             weightSpr(sprL) = weightSprVal
             weightSpr(sprR) = weightSprVal
             
          elseif (i==4) then
             ltSpr2L(j) = (cellL-1)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             ltSpr2R(j) = (cellR-1)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             
             sprL = ltSpr2L(j) ; sprR = ltSpr2R(j)
             sprLgcl(sprL) = 1 ; sprLgcl(sprR) = 1
             
             weightSpr(sprL) = weightSprVal
             weightSpr(sprR) = weightSprVal
             
          endif
          
       enddo
       
    enddo
    
    open(unit=40,file='lgcl_abd_weighted.dat')
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          if (i==1) write(40,*) cellLgcl(j),weightCell(j),j,"lgcl & weight of cell"
          if (i==2) write(40,*) sprLgcl(j), weightSpr(j) ,j,"lgcl & weight of spr"
       enddo
       
       write(40,*) " "
       
    enddo
    
    close(40)
    
  end subroutine lgcls_and_weighted_arrays
  
  subroutine get_lgclsForCellsMatch(cellLgcl)
    implicit none
    integer, intent(out) :: cellLgcl(1:N_cell)
    integer :: i
    
    !cellLgcl(1:6)   = 1 ; cellLgcl((Hlf_Ncell+1):(Hlf_Ncell+6))   = 1
    !cellLgcl(7)     = 0 ; cellLgcl(Hlf_Ncell+7)                   = 0
    !cellLgcl(8)     = 0 ; cellLgcl(Hlf_Ncell+8)                   = 0
    !cellLgcl(9)     = 1 ; cellLgcl(Hlf_Ncell+9)                   = 1
    !cellLgcl(10:11) = 0 ; cellLgcl((Hlf_Ncell+10):(Hlf_Ncell+11)) = 0

    cellLgcl(1:6)   = 0 ; cellLgcl((Hlf_Ncell+1):(Hlf_Ncell+6))   = 0
    cellLgcl(7)     = 0 ; cellLgcl(Hlf_Ncell+7)                   = 0
    cellLgcl(8)     = 0 ; cellLgcl(Hlf_Ncell+8)                   = 0
    cellLgcl(9)     = 0 ; cellLgcl(Hlf_Ncell+9)                   = 0
    cellLgcl(10)    = 1 ; cellLgcl((Hlf_Ncell+10))                = 1
    cellLgcl(11)    = 0 ; cellLgcl((Hlf_Ncell+11))                = 0
    
    cellLgcl(N_cell) = 0
    
    open(unit=23,file='lgclsForCellsMatch.dat')

    do i = 1,N_cell
       write(23,*) cellLgcl(i),i,"cellLgcl"
    enddo
    
    close(23)
    
  end subroutine get_lgclsForCellsMatch
  
  
  subroutine get_sprCur_and_sprNxt(cellCur,cellNxt,sprsInCell,sprCntEachCell,sprCurs,sprNxts)
    implicit none
    integer, intent(in)  :: cellCur,cellNxt
    integer, intent(in)  :: sprsInCell,sprCntEachCell
    integer, intent(out) :: sprCurs(1:sprsInCell),sprNxts(1:sprsInCell)
    
    integer :: sprCnt
    integer :: i,j,jmax
    integer :: apSprCur(1:(NAEC_Apcl+1)) ,apSprNxt(1:(NAEC_Apcl+1))
    integer :: bsSprCur(1:(NAEC_Bsal+1)) ,bsSprNxt(1:(NAEC_Bsal+1))
    integer :: ltSpr1Cur(1:(NAEC_Ltrl+1)),ltSpr1Nxt(1:(NAEC_Ltrl+1))
    integer :: ltSpr2Cur(1:(NAEC_Ltrl+1)),ltSpr2Nxt(1:(NAEC_Ltrl+1))

    open(unit=21,file='sprCurs_and_sprNxts.dat',position='append')
    write(21,*) cellCur,cellNxt,sprsInCell,"cellCur-Nxt and sprsInCell"
    
    sprCnt = 1
    
    do i = 1,4
       
       if (i==1) jmax = NAEC_Ltrl+1
       if (i==2) jmax = NAEC_Apcl+1
       if (i==3) jmax = NAEC_Bsal+1
       if (i==4) jmax = NAEC_Ltrl+1 
       
       do j = 1,jmax
          
          if (i==1) then
             
             if ((cellCur==1) .or. (cellCur==(Hlf_Ncell+1))) then
                continue
             else
                
                ltSpr1Cur(j) = (cellCur-2)*(sprCntEachCell)+ (NAEC_Apcl+1)+ (NAEC_Bsal+1)+ j
                ltSpr1Nxt(j) = (cellNxt-2)*(sprCntEachCell)+ (NAEC_Apcl+1)+ (NAEC_Bsal+1)+ j
                
                sprCurs(sprCnt) = ltSpr1Cur(j) ; sprNxts(sprCnt) = ltSpr1Nxt(j)
                sprCnt = sprCnt+1
                
             endif
                
          elseif (i==2) then
             
             apSprCur(j) = (cellCur-1)*(sprCntEachCell) + j
             apSprNxt(j) = (cellNxt-1)*(sprCntEachCell) + j
             
             sprCurs(sprCnt) = apSprCur(j)  ; sprNxts(sprCnt) = apSprNxt(j)
             sprCnt = sprCnt+1
             
          elseif (i==3) then
             
             bsSprCur(j) = (cellCur-1)*(sprCntEachCell) + (NAEC_Apcl+1) + j
             bsSprNxt(j) = (cellNxt-1)*(sprCntEachCell) + (NAEC_Apcl+1) + j
             
             sprCurs(sprCnt) = bsSprCur(j)  ; sprNxts(sprCnt) = bsSprNxt(j)
             sprCnt = sprCnt+1
             
          elseif (i==4) then
             
             ltSpr2Cur(j) = (cellCur-1)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             ltSpr2Nxt(j) = (cellNxt-1)*(sprCntEachCell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             
             sprCurs(sprCnt) = ltSpr2Cur(j) ; sprNxts(sprCnt) = ltSpr2Nxt(j)
             sprCnt = sprCnt+1
             
          endif
          
       enddo
       
    enddo
    
    if ((sprCnt-1) .ne. sprsInCell) then
       write(*,*) (sprCnt-1),sprsInCell,"sprCnt & sprsInCell not matching"
       stop
    endif

    do i = 1,sprsInCell
       write(21,*) sprCurs(i),sprNxts(i),(sprNxts(i)-sprCurs(i)),i,"sprCurs-nxt-Diff-i"
    enddo
    
    write(21,*) " "
    
    close(21)
    
  end subroutine get_sprCur_and_sprNxt


  
  subroutine symmetricize_node_for_NIsystem(LftRghtCopy)
    implicit none
    integer, intent(in) :: LftRghtCopy
    
    integer :: rglrNode,LrglrNode,RrglrNode
    integer :: addtNode,LaddtNode,RaddtNode
    integer :: tot_Lnodes,tot_Rnodes
    integer :: hlf_rglrNode,hlf_addtNode
    
    integer :: i,j,jmax,m
    integer :: cntNodeS,cntNode,cnt_nodeSkip
    integer :: lftNode,rghtNode
    real*8  :: x,y
    integer :: nodeNm,reason
    real*8  :: ForceBfr(1:N_node,1:2),ForceAft(1:N_node,1:2)
    real*8  :: ffc
    integer :: NnodesIn
    real*8  :: Ebfr,Eaft
    
    integer, allocatable :: LNodes(:),RNodes(:)
    
    if (modelID==1) then
       write(*,*) "not written for TN system"
       stop
    elseif (NAEC.ne.2) then
       write(*,*) "not written for NAEC not equals 2 system"
       stop
    endif
    
    open(unit=91,file='symmNode.dat')
    open(unit=92,file='grdForce.dat')
    open(unit=94,file='nodeValues.dat')
    
    rglrNode  = (N_nodeS-1)
    LrglrNode = (rglrNode/2)
    RrglrNode = (rglrNode/2)
    
    addtNode  = (N_node-1) - (rglrNode) - 2
    LaddtNode = (addtNode/2)
    RaddtNode = (addtNode/2)
    
    tot_Lnodes = LrglrNode + LaddtNode + 1
    tot_Rnodes = RrglrNode + RaddtNode + 1
    
    write(91,*) tot_Lnodes,tot_Rnodes,"tot_Lnodes,tot_Rnodes"
    
    if (tot_Lnodes.ne.tot_Rnodes) then
       write(*,*) "tot_Lnodes .ne. tot_Rnodes"
       stop
    endif
    
    hlf_rglrNode = (rglrNode/2)
    hlf_addtNode = (addtNode/2)
    
    write(91,*) rglrNode,LrglrNode,RrglrNode,"rglrNode, L/R rglrNode"
    write(91,*) addtNode,LaddtNode,RaddtNode,"addtNode, L/R addtNode"
    write(91,*) hlf_rglrNode,hlf_addtNode,"hlf rglr, hlf addt"
    write(91,*) " "
    
    write(91,*) tot_Lnodes,tot_Rnodes,"tot L/R nodes"
    allocate(LNodes(1:tot_Lnodes),RNodes(1:tot_Rnodes))
    cntNode = 1
    
    
    do i = 1,2
       
       if (i==1) then
          jmax     = hlf_rglrNode
          cntNodeS = 0
          
       elseif (i==2) then
          jmax     = (hlf_addtNode+1)
          cntNodeS = rglrNode+1
       endif
       
       do j = 1,jmax
          
          if (i==1) then
             lftNode  = cntNodeS + j
             rghtNode = hlf_rglrNode + cntNodeS + j
             
             write(91,*) cntNode,"cntNode"
             LNodes(cntNode) = lftNode ; RNodes(cntNode) = rghtNode
             cntNode = cntNode+1
             
          elseif (i==2) then
             
             if (j.ne.jmax) then
                lftNode  = cntNodeS + j
                rghtNode = hlf_addtNode + cntNodeS + j
                
                write(91,*) cntNode,"cntNode"
                LNodes(cntNode) = lftNode ; RNodes(cntNode) = rghtNode
                cntNode = cntNode+1
                
             elseif (j==jmax) then
                lftNode  = N_node-1
                rghtNode = N_node
                
                write(91,*) cntNode,"cntNode"
                LNodes(cntNode) = lftNode ; RNodes(cntNode) = rghtNode
                cntNode = cntNode+1
                
             endif
             
          endif
          
          write(91,*) lftNode,rghtNode,"lft-rght"
       enddo
       
    enddo
    
    do i = 1,tot_Lnodes
       write(91,*) LNodes(i),RNodes(i),i,"LNodes,RNodes,i"
    enddo
    
    filechk = 1
    
    Ebfr = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    ForceBfr = -gradient
    
    do i = 1,N_node
       write(92,*) ForceBfr(i,1:2),i,"ForceBfr"
    enddo
    
    write(92,*) " "
    
    do i = 1,3
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             ffc = k_area(j) * (A0(j)-A(j))
             write(92,*) ffc,k_area(j),A0(j),A(j),j,"ka(A0-A) [Bfr]"
             
             NnodesIn = area_node(j,0)
             write(94,*) area_node(j,1:NnodesIn),j,"nodesInarea [Bfr]"
             
             do m = 1,NnodesIn
                nodeNm = area_node(j,m) 
                write(94,*) node_xy(m,1:2),m,"node_xy [Bfr]"
             enddo
             
          elseif (i==2) then
             ffc = k_spr(j)  * (l0(j)-l(j))
             write(92,*) ffc,k_spr(j),l0(j),l(j),j,"ks(l0-l) [Bfr]"
             
             NnodesIn = spr_node(j,0)
             write(94,*) spr_node(j,1:NnodesIn),j,"nodesInspr [Bfr]"
             
             do m = 1,NnodesIn
                nodeNm = spr_node(j,m)
                write(94,*) node_xy(m,1:2),m,"node_xy [Bfr]"
             enddo
             
          elseif (i==3) then
             ffc = CgXNode(j)
             write(92,*) ffc,j,"CgXNode [Bfr]"
             
          endif
          
       enddo
       
       write(92,*) " "
       write(94,*) " "
       
    enddo
    
    
    do i = 1,tot_Lnodes
       
       lftNode = LNodes(i) ; rghtNode = RNodes(i)
       
       if (LftRghtCopy==1) then
          node_xy(lftNode,1)  = -node_xy(rghtNode,1)
          node_xy(lftNode,2)  = node_xy(rghtNode,2)
          
       elseif (LftRghtCopy==2) then
          node_xy(rghtNode,1) = -node_xy(lftNode,1)
          node_xy(rghtNode,2) = node_xy(lftNode,2)
          
       endif
       
    enddo
    
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:2),i,"i~i"
    enddo
    
    close(91)
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    Eaft = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    ForceAft = -gradient
    
    filechk = 0
    
    write(92,*) " "
    write(94,*) " "
    
    do i = 1,N_node
       write(92,*) ForceAft(i,1:2),i,"ForceAft"
    enddo
    
    
    do i = 1,3
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             ffc = k_area(j) * (A0(j)-A(j))
             write(92,*) ffc,k_area(j),A0(j),A(j),j,"ka(A0-A) [Aft]"
             
             NnodesIn = area_node(j,0)
             write(94,*) area_node(j,1:NnodesIn),j,"nodesInarea [Aft]"
             
             do m = 1,NnodesIn
                nodeNm = area_node(j,m) 
                write(94,*) node_xy(m,1:2),m,"node_xy [Aft]"
             enddo
             
          elseif (i==2) then
             ffc = k_spr(j)  * (l0(j)-l(j))
             write(92,*) ffc,k_spr(j),l0(j),l(j),j,"ks(l0-l) [Aft]"
             
             NnodesIn = spr_node(j,0)
             write(94,*) spr_node(j,1:NnodesIn),j,"nodesInspr [Aft]"
             
             do m = 1,NnodesIn
                nodeNm = spr_node(j,m)
                write(94,*) node_xy(m,1:2),m,"node_xy [Aft]"
             enddo
             
          elseif (i==3) then
             ffc = CgXNode(j)
             write(92,*) ffc,j,"CgXNode [Aft]"
             
          endif
          
       enddo
       
       write(92,*) " "
       write(94,*) " "
       
    enddo
    
    
    close(92)
    close(94)
    
  end subroutine symmetricize_node_for_NIsystem
  
  subroutine read_nodeXY_frmDatFile(LNodes,RNodes,tot_Lnodes,tot_Rnodes)
    implicit none
    integer, intent(in) :: LNodes(1:tot_Lnodes),RNodes(1:tot_Rnodes)
    integer, intent(in) :: tot_Lnodes,tot_Rnodes
    
    integer :: cnt_nodeSkip
    real*8  :: x,y
    integer :: nodeNm
    integer :: lftNode,rghtNode
    integer :: i,j,reason
    
    open(unit=11,file='NI_modelNW006.dat')
    
    i=1
    
    cnt_nodeSkip = 0
    
    do
       read(11,*,IOSTAT=reason) x,y,nodeNm
       
       if ((i+cnt_nodeSkip).ne.nodeNm) then
          
          write(*,*) (i+cnt_nodeSkip),nodeNm,"SKIP"
          
          do j = 1,tot_LNodes
             lftNode = LNodes(j) ; rghtNode = RNodes(j)
             
             if (rghtNode==(i+cnt_nodeSkip)) then
                node_xy(rghtNode,1:2) = node_xy(lftNode,1:2)
                write(*,*) node_xy(lftNode,1:2),"lft node_xy"
                write(*,*) nodeNm,lftNode,rghtNode,"node"
                exit
             endif
             
          enddo
          
          cnt_nodeSkip = cnt_nodeSkip+1
          write(*,*) cnt_nodeSkip,"cnt_nodeSkip val"
          
       endif
       
       write(*,*) reason,"iostat"
       
       if (reason.ne.0) exit
       
       node_xy(nodeNm,1) = x ; node_xy(nodeNm,2) = y
       write(*,*) node_xy(nodeNm,1:2),nodeNm,i,"node_xy,nodeNm,i"
       i = i+1
    enddo
    
    write(*,*) " "
    write(*,*) cnt_nodeSkip,"Skip"
    
    close(11)
    
    
    do i = 1,N_node
       write(*,*) node_xy(i,1:2),i,"i~"
    enddo
    
  end subroutine read_nodeXY_frmDatFile
  
  
  subroutine convrsn_of_TNprop_aft_chnging_NIprop
    implicit none
    character(len=100) :: buffrFlnm
    
    integer :: NsprStr,NcellStr,NnodeStr
    integer :: apclTN,bsalTN,ltrlTN
    integer :: sprInApcl,sprInBsal,sprInLtrl
    integer :: nsprsInACellTN,nsprsInACellNI
    integer :: NAEC_val
    integer :: cntNI,sprNm
    integer :: i,j,jmax
    
    integer, allocatable :: apclSpr(:),bsalSpr(:),ltrlSpr(:)
    real*8 , allocatable :: ksNI(:),l0NI(:)
    real*8, allocatable  :: ksStr(:),kaStr(:)
    real*8, allocatable  :: l0Str(:),A0Str(:)
    
    if (modelID==1) then
       write(*,*) "Not compatible for modelID=1, we need NI model at the beginning"
    elseif (modelID==2) then
       continue
    endif
    
    buffrFlnm='NI_propStorge.dat'
    call save_NIpropInABuffrFile(buffrFlnm)
    NsprStr = N_spr ; NcellStr = N_cell ; NnodeStr = N_node
    
    allocate(ksStr(1:N_spr),kaStr(1:N_cell))
    allocate(l0Str(1:N_spr),A0Str(1:N_cell))
    
    ksStr = k_spr ; kaStr = k_area
    l0Str = l0    ; A0Str = A0
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    sprInApcl = NAEC_apcl+1 ; sprInBsal = NAEC_bsal+1 ; sprInLtrl = NAEC_ltrl+1
    allocate(apclSpr(1:sprInApcl),bsalSpr(1:sprInBsal),ltrlSpr(1:sprInLtrl))
    apclSpr=0 ; bsalSpr=0 ; ltrlSpr=0
    
    nsprsInACellTN = 3
    nsprsInACellNI = (sprInApcl)+(sprInBsal)+(sprInLtrl)
    
    do i = 1,N_cell
       
       if (i.ne.N_cell) jmax = nsprsInACellTN
       
       if (i==N_cell) then
          if (cellsMeet==0)    jmax = 2
          if (cellsMeet.gt.1)  jmax = 1
       endif
       
       if (jmax==3) then
          
          apclTN = (i-1)*(nsprsInACellTN) + 1 
          bsalTN = (i-1)*(nsprsInACellTN) + 2
          ltrlTN = (i-1)*(nsprsInACellTN) + 3
          
       elseif (jmax==2) then
          
          apclTN = (i-1)*(nsprsInACellTN) + 1
          bsalTN = (i-1)*(nsprsInACellTN) + 2
          
       elseif (jmax==1) then
          bsalTN = (i-1)*(nsprsInACellTN) + 1
       endif
       
       do j = 1,jmax
          
          if (j==1) then
             
             if (jmax.ne.1) then
                
                sprNm = apclTN ; NAEC_val = NAEC_apcl
                allocate(ksNI(1:(NAEC_val+1)))

             elseif (jmax==1) then
                
                sprNm = bsalTN ; NAEC_val = NAEC_bsal
                allocate(ksNI(1:(NAEC_val+1)))
             endif
             
          elseif (j==2) then
             sprNm = bsalTN ; NAEC_val = NAEC_bsal
             allocate(ksNI(1:(NAEC_val+1)))
             
          elseif (j==3) then
             sprNm = ltrlTN ; NAEC_val = NAEC_ltrl
             allocate(ksNI(1:(NAEC_val+1)))
          endif
          
          ksNI(1:(NAEC_val+1)) = ksStr(cntNI:(cntNI+NAEC_val+1))
          l0NI(1:(NAEC_val+1)) = l0Str(cntNI:(cntNI+NAEC_val+1))
          
          k_spr(sprNm)         = ksprTNsystem(ksNI,NAEC_val)
          l0(sprNm)            = l0TNsystem(ksNI,NAEC_val)
          
          cntNI = cntNI + (NAEC_val+1)
          
       enddo
       
    enddo
    
    deallocate(ksStr,kaStr)
    deallocate(l0Str,A0Str)
    
  end subroutine convrsn_of_TNprop_aft_chnging_NIprop
  
  subroutine  save_NIpropInABuffrFile(buffrFlnm)
    implicit none
    character(len=100), intent(in) :: buffrFlnm
    
    integer :: i,j,jmax
    integer :: N_itm
    
    open(unit=212,file=trim(adjustl(buffrFlnm)))
    
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) write(212,*) k_spr(j),l0(j),j
          if (i==2) write(212,*) k_area(j),A0(j),j
          if (i==3) write(212,*) CgXNode(j),j
       enddo
       
       write(212,*) " "
       
    enddo
    
    close(212)
    
  end subroutine save_NIpropInABuffrFile
  
  real*8 function ksprTNsystem(ksNI,NAEC_val)
    implicit none
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: ksNI(1:(NAEC_val+1))
    
    integer :: i,j
    real*8  :: reverse_kspr
    
    reverse_kspr = 0.0d0
    
    do i = 1,(NAEC_val+1)
       reverse_kspr = reverse_kspr + (1.0d0)/(ksNI(i))
    enddo
    
    ksprTNsystem = (1.0d0)/(reverse_kspr)
    
  end function ksprTNsystem
  
  real*8 function l0TNsystem(l0NI,NAEC_val)
    implicit none
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: l0NI(1:(NAEC_val+1))
    
    integer :: i
    real*8  :: suml0
    
    suml0 = 0.0d0
    
    do i = 1,(NAEC_val+1)
       suml0 = suml0 + l0NI(i)
    enddo
    
    l0TNsystem = suml0
    
  end function l0TNsystem
  
  subroutine write_in_given_file(flnm)
    implicit none
    character(len=100), intent(in) :: flnm
    integer :: i,j,jmax,N_itm
    
    write(*,*) trim(adjustl(flnm))
    N_itm = 4
    
    open(unit=901,file=trim(adjustl(flnm)))
    
    do i = 1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(901,*) k_spr(j),l0(j)
          if (i==2) write(901,*) k_area(j),A0(j)
          if (i==3) write(901,*) CgXNode(j)
          if (i==4) write(901,*) node_xy(j,1:2)
       enddo
       
       write(901,*) " "
       
    enddo
    
    close(901)
    
  end subroutine write_in_given_file
  
  
  subroutine Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    write(*,*) " "
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    call switchto_NI_model_run_and_switchbackto_TN
    
  end subroutine Equilibrate_bothTN_NImodel
  
  
  subroutine Equilibrate_bothTN_NImodel_with_PreSimltd_DatFile(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer            :: objTyp
    character(len=100) :: flnm
    real*8             :: E
    
    call get_nodeXY_frmNxtFile(ExpNo,FrmNo)
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    call switch_to_NI_model
    call get_nodeXY_frmNxtFile(Exprmnt_NI,Frame_NI)
    call nodes_to_coordntes(node_xy,coordntes_xy)
    E = Energy(node_xy,l0,A0) ; write(*,*) E,"E in PreSimltd_DatFile"
    call print_prop_of_system
    call get_gradient(node_xy,l0,A0,gradient)
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine Equilibrate_bothTN_NImodel_with_PreSimltd_DatFile
    
  subroutine notEquilibrate_saveTN_switchtoNI_notEquilibrate(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    call switch_to_NI_and_not_Equilibrate
    call not_Equlibrate_backto_TN
    
  end subroutine notEquilibrate_saveTN_switchtoNI_notEquilibrate
  
  subroutine check_trnsfrms_bfrAft_SystemChange(BfrOrAft,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: BfrorAft
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,j,imax,jmax
    integer :: nelmntInARow
    integer, allocatable :: rowVal(:)
    
    nelmntInARow = 14
    allocate(rowVal(1:nelmntInARow))
    
    if (BfrOrAft==-1) then
       
       open(unit=323,file='trnsfrm_bfr_TN.dat')
       open(unit=324,file='trnsfrm_bfr_NI.dat')
    elseif (BfrOrAft==+1) then
       open(unit=323,file='trnsfrm_aft_TN.dat')
       open(unit=324,file='trnsfrm_aft_NI.dat')
    endif
    
    imax = 6
    
    do i = 1,imax 
       if (i.le.2)                  jmax=N_node
       if ((i.gt.2) .and. (i.le.4)) jmax=N_spr
       if ((i.gt.4) .and. (i.le.6)) jmax=N_cell
       do j = 1,jmax
          call getting_rowVal(i,j,nelmntInARow,rowVal)
          write(323,*) rowVal(1:nelmntInaRow)
          write(*,*) rowVal(1:nelmntInaRow),"rowVal"
       enddo
    enddo
    
    write(*,*) "ENTEred in the sbrtn1"
    call sleep(2)
    call switch_to_NI_and_Equilibrate
    
    imax = 6
    
    do i = 1,imax    
       if (i.le.2)                  jmax=N_node
       if ((i.gt.2) .and. (i.le.4)) jmax=N_spr
       if ((i.gt.4) .and. (i.le.6)) jmax=N_cell
       do j = 1,jmax
          call getting_rowVal(i,j,nelmntInARow,rowVal)
          write(324,*) rowVal(1:nelmntInaRow)
          write(*,*) rowVal(1:nelmntInaRow),"rowVal"
       enddo
    enddo
    write(*,*) "ENTEred in the sbrtn2"
    
    call aftrEqulibrate_backto_TN
    
    close(323)
    close(324)
    
  end subroutine Check_Trnsfrms_BfrAft_SystemChange

  subroutine getting_rowVal(iVal,jVal,nelmntInARow,rowVal)
    implicit none
    integer, intent(in)  :: iVal,jVal,nelmntInARow
    integer, intent(out) :: rowVal(1:nelmntInARow)
    
    rowVal(1:nelmntInARow) = -1
    rowVal(1) = jVal
    
    if (iVal==1) rowVal(2:(max_spr_node+2))  = node_spr(jVal,0:max_spr_node)
    if (iVal==2) rowVal(2:(max_area_node+2)) = node_area(jVal,0:max_area_node)
    if (iVal==3) rowVal(2:(max_node_spr+2))  = spr_node(jVal,0:max_node_spr)
    if (iVal==4) rowVal(2:(max_area_spr+2))  = spr_area(jVal,0:max_area_spr)
    if (iVal==5) rowVal(2:(max_node_area+2)) = area_node(jVal,0:max_node_area)
    if (iVal==6) rowVal(2:(max_spr_area+2))  = area_spr(jVal,0:max_spr_area)
    
  end subroutine getting_rowVal
  
  subroutine check_grd_at_both_TNandNI
    implicit none
    real*8, allocatable :: ksSave(:),l0Save(:)
    real*8, allocatable :: kaSave(:),A0Save(:)
    real*8, allocatable :: CgXSave(:)
    real*8, allocatable :: grdValA(:,:),grdValN(:,:) ,diff_AorN(:,:)
    
    integer :: AnaOrNumSave
    integer :: chk_com,i
    real*8  :: E
    real*8  :: TOL_diff=2.0e-8
    integer :: mrk1,mrk2
    
    open(unit=82,file='chkTN_GRD.dat',position='append')
    open(unit=83,file='chkNI_GRD.dat',position='append')
    
    allocate(ksSave(1:N_spr) ,l0Save(1:N_spr))
    allocate(kaSave(1:N_cell),A0Save(1:N_cell))
    allocate(CgXSave(1:N_node))
    
    allocate(grdValA(1:N_node,1:N_dmnsn))
    allocate(grdValN(1:N_node,1:N_dmnsn))
    allocate(diff_AorN(1:N_node,1:N_dmnsn))
    
    grdValA = 0.00d0 ; grdValN = 0.0d0 ; diff_AorN = -1.0d30
    
    ksSave  = k_spr  ; l0Save = l0
    kaSave  = k_area ; A0Save = A0
    CgXSave = CgXNode
    
    AnaOrNumSave = Analtcl_or_Numrcl
    
    chk_com = 1
    
    if (chk_com==1) then
       k_area = 0.00d0 ; CgXNode = 0.00d0
    elseif (chk_com==2) then
       k_spr  = 0.0d0  ; CgXNode = 0.00d0
    elseif (chk_com==3) then
       k_spr  = 0.00d0 ; k_area = 0.00d0
    endif
    
    Analtcl_or_Numrcl = 1
    call get_gradient(node_xy,l0,A0,grdValA)
    Analtcl_or_Numrcl = 2
    call get_gradient(node_xy,l0,A0,grdValN)
    
    
    do i = 1,N_node
       diff_AorN(i,1:2) = grdValA(i,1:2)-grdValN(i,1:2)
       
       if (abs(diff_AorN(i,1)) .le. TOL_diff) mrk1=0
       if (abs(diff_AorN(i,1)) .gt. TOL_diff) mrk1=-1
       if (abs(diff_AorN(i,2)) .le. TOL_diff) mrk2=0
       if (abs(diff_AorN(i,2)) .gt. TOL_diff) mrk2=-1
       
       write(82,*) grdValA(i,1:2),grdValN(i,1:2),diff_AorN(i,1:2),i,mrk1,mrk2
    enddo
    
    
    k_spr   = ksSave
    k_area  = kaSave
    CgXNode = CgXSave
    
    Analtcl_or_Numrcl = AnaOrNumSave
    
    deallocate(ksSave,l0Save)
    deallocate(kaSave,A0Save)
    deallocate(CgXSave)
    
    deallocate(grdValA)
    deallocate(grdValN)
    deallocate(diff_AorN)
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E"
    call get_gradient(node_xy,l0,A0,gradient)
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    !Frame_NI = Frame_NI+1
    
    allocate(ksSave(1:N_spr) ,l0Save(1:N_spr))
    allocate(kaSave(1:N_cell),A0Save(1:N_cell))
    allocate(CgXSave(1:N_node))
    
    allocate(grdValA(1:N_node,1:N_dmnsn))
    allocate(grdValN(1:N_node,1:N_dmnsn))
    allocate(diff_AorN(1:N_node,1:N_dmnsn))
    
    grdValA = 0.00d0 ; grdValN = 0.0d0 ; diff_AorN = -1.0d30
    
    ksSave  = k_spr  ; l0Save = l0
    kaSave  = k_area ; A0Save = A0
    CgXSave = CgXNode
    
    AnaOrNumSave = Analtcl_or_Numrcl
    
    chk_com = 1
    
    if (chk_com==1) then
       k_area = 0.00d0 ; CgXNode = 0.00d0
    elseif (chk_com==2) then
       k_spr  = 0.0d0  ; CgXNode = 0.00d0
    elseif (chk_com==3) then
       k_spr  = 0.00d0 ; k_area = 0.00d0
    endif
    
    Analtcl_or_Numrcl = 1
    call get_gradient(node_xy,l0,A0,grdValA)
    Analtcl_or_Numrcl = 2
    call get_gradient(node_xy,l0,A0,grdValN)
    
    do i = 1,N_node
       diff_AorN(i,1:2) = grdValA(i,1:2)-grdValN(i,1:2)
       
       if (abs(diff_AorN(i,1)) .le. TOL_diff) mrk1=0
       if (abs(diff_AorN(i,1)) .gt. TOL_diff) mrk1=-1
       if (abs(diff_AorN(i,2)) .le. TOL_diff) mrk2=0
       if (abs(diff_AorN(i,2)) .gt. TOL_diff) mrk2=-1
       
       write(83,*) grdValA(i,1:2),grdValN(i,1:2),diff_AorN(i,1:2),i,mrk1,mrk2
    enddo
    
    k_spr   = ksSave
    k_area  = kaSave
    CgXNode = CgXSave
    
    Analtcl_or_Numrcl = AnaOrNumSave
    
    deallocate(ksSave,l0Save)
    deallocate(kaSave,A0Save)
    deallocate(CgXSave)
    
    deallocate(grdValA)
    deallocate(grdValN)
    deallocate(diff_AorN)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(82)
    close(83)
    
  end subroutine check_grd_at_both_TNandNI
  
  subroutine get_nodeXY_frmNxtFile(ExpNo,FrmNo)
    implicit none
    integer, intent(in)  :: ExpNo
    integer, intent(out) :: FrmNo
    
    integer            :: objTyp,cntLp,error
    character(len=100) :: flnm
    
    
    objTyp = 1
    call get_fileNames_inside_the_modFile(ExpNo,FrmNo,objTyp,flnm)
    
    write(*,*) trim(adjustl(flnm))
    
    open(unit=249,file=trim(adjustl(flnm)))

    cntLp = 1
    
    do
       read(unit=249,fmt=*,IOSTAT=error) node_xy(cntLp,1:N_dmnsn)
       write(*,*) node_xy(cntLp,1:2),cntLp,"in check"
       cntLp = cntLp+1
       if (error.lt.0) exit
    enddo
    
    close(249)
    
  end subroutine get_nodeXY_frmNxtFile

  subroutine print_prop_of_system
    implicit none
    integer :: i,j,jmax
    real*8  :: Es,Ea,Eg,E_tot
    real*8  :: SumEs=0.0d0,SumEa=0.0d0,SumEg=0.0d0
    
    do i = 1,3
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             Es = 0.50d0*k_spr(j)*(l(j)-l0(j))**2
             write(*,*) k_spr(j),l(j),l0(j),Es,j,"Es"
             SumEs = SumEs+Es
             
          elseif (i==2) then
             
             Ea = 0.50d0*k_area(j)*(A(j)-A0(j))**2
             write(*,*) k_area(j),A(j),A0(j),Ea,j,"Ea"
             SumEa = SumEa+Ea
             
          elseif (i==3) then
             
             Eg = CgXNode(j)*(node_xy(j,1)) + CgYNode(j)*(node_xy(j,2))
             write(*,*) CgXNode(j),CgYNode(j),node_xy(j,1:2),Eg,j,"Eg"
             SumEg = SumEg+Eg
             
          endif
          
       enddo
       
    enddo
    
    E_tot = Es+Ea+Eg
    
    write(*,*) E_tot,"E_tot"
    
    if (E_tot.gt.300.0d0) then
       write(*,*) "E_tot HIGH"
       stop
    endif
    
  end subroutine print_prop_of_system


  subroutine print_prop_of_system_wtSRyp
    implicit none
    integer :: i,j,jmax
    real*8  :: Es,Ea,Eg,Er,E_tot
    real*8  :: SumEs=0.0d0,SumEa=0.0d0,SumEg=0.0d0,SumEr=0.0d0
    
    do i = 1,4
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=numApclNodes
       
       do j = 1,jmax
          
          if (i==1) then
             
             Es = 0.50d0*k_spr(j)*(l(j)-l0(j))**2
             write(*,*) k_spr(j),l(j),l0(j),Es,j,"Es"
             SumEs = SumEs+Es
             
          elseif (i==2) then
             
             Ea = 0.50d0*k_area(j)*(A(j)-A0(j))**2
             write(*,*) k_area(j),A(j),A0(j),Ea,j,"Ea"
             SumEa = SumEa+Ea
             
          elseif (i==3) then
             
             Eg = CgXNode(j)*(node_xy(j,1)) + CgYNode(j)*(node_xy(j,2))
             write(*,*) CgXNode(j),CgYNode(j),node_xy(j,1:2),Eg,j,"Eg"
             SumEg = SumEg+Eg
             
          elseif (i==4) then
             
             Er = (AmpApclYP(j)*exp(AlphaApclYP(j)*(node_xy(apclNodes(j),2)))) / &
                  ((-node_xy(apclNodes(j),2)+epsApclYP(j))**(powrApclYP(j)))
             write(*,*) AmpApclYP(j),AlphaApclYP(j),epsApclYP(j),powrApclYP(j),Er,apclNodes(j),j,"Er"
             SumEr = SumEr+Er
             
          endif
          
       enddo
       
    enddo
    
    E_tot = SumEs+SumEa+SumEg+SumEr
    
    write(*,*) E_tot,"E_tot"
    
    if (E_tot.gt.300.0d0) then
       write(*,*) "E_tot HIGH"
    endif
    
  end subroutine print_prop_of_system_wtSRyp
  
  
  subroutine switch_to_TN_model
    implicit none
    
    modelID             = 1
    switchBackFrmNItoTN = 1
    
    call allocate_variabls_neededToStore_for_NIsys
    
    call store_NI_SysVars_and_Arrays
    call get_TN_Sysvars
    call deallocate_and_reallocate_TN_Arrays
    
    call get_TN_NodeVars
    call get_TN_SprVars
    call get_TN_kphi_and_nodePhiTyp
    call get_TN_CgXvars                ; call get_TN_CgYvars
    
    call store_all_moving_coordnte_variables_NI
    call deallocate_moving_coordnte_variables_wo_StrVars 
    call get_all_moving_coordnte_variables_wo_StrVars
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    call store_NI_trnsfrms
    call deallocate_and_reallocate_for_TN_trnsfrmVars
    call get_all_the_transfrms
    call print_all_trnsfrms
    
    call store_all_gradient_variables_NI
    call deallocate_all_gradient_variables_wo_StrVars
    call get_all_gradient_variables_wo_StrVars
    
  end subroutine switch_to_TN_model
  
  subroutine allocate_variabls_neededToStore_for_NIsys
    implicit none
    
    allocate(node_xyStrNI(1:N_node,1:N_dmnsn))            ; node_xyStrNI       = -1.0d30
    allocate(node_typStrNI(1:N_node))                     ; node_typStrNI      = -1
    allocate(node_cnnctdStrNI(1:N_node))                  ; node_cnnctdStrNI   = -1
    allocate(double_nodeStrNI(1:N_node,1:N_dmnsn))        ; double_nodeStrNI   = -1 
    allocate(count_this_dnStrNI(1:N_node))                ; count_this_dnStrNI = -1
    
    allocate(typ_sprStrNI(1:N_spr))                       ; typ_sprStrNI = -1
    allocate(k_sprStrNI(1:N_spr))                         ; k_sprStrNI   = -1.0d30
    allocate(l0_StrNI(1:N_spr))                           ; l0_StrNI     = -1.0d30
    allocate(l_StrNI(1:N_spr))                            ; l_StrNI      = -1.0d30
    
    allocate(nodePhi_typStrNI(1:N_node))                  ; nodePhi_typStrNI = -1
    allocate(k_phiStrNI(1:N_node,1:max_Phi_node))         ; k_phiStrNI       = -1.0d30
    allocate(CgXNode_StrNI(1:N_node))                     ; CgXNode_StrNI    = -1.0d30 
    allocate(CgYNode_StrNI(1:N_node))                     ; CgYNode_StrNI    = -1.0d30
    
    allocate(coordntesStrNI(1:N_mvCoordnte_withl0A0))     ; coordntesStrNI              = -1.0d30
    allocate(coordntes_xyStrNI(1:N_mvCoordnte))           ; coordntes_xyStrNI           = -1.0d30
    allocate(coordnte_of_which_nodeStrNI(1:N_mvCoordnte)) ; coordnte_of_which_nodeStrNI = -1
    allocate(x_or_yStrNI(1:N_mvCoordnte))                 ; x_or_yStrNI                 = -1
    
    allocate(spr_nodeSNI(1:N_spr,0:max_node_spr))         ; spr_nodeSNI  = -1
    allocate(spr_areaSNI(1:N_spr,0:max_area_spr))         ; spr_areaSNI  = -1 
    allocate(node_sprSNI(1:N_node,0:max_spr_node))        ; node_sprSNI  = -1
    allocate(node_areaSNI(1:N_node,0:max_spr_node))       ; node_areaSNI = -1
    allocate(area_sprSNI(1:N_cell,0:max_spr_area))        ; area_sprSNI  = -1
    allocate(area_nodeSNI(1:N_cell,0:max_node_area))      ; area_nodeSNI = -1
    allocate(NlistSNI(1:N_node,1:max_area_node,1:2))      ; NlistSNI     = -1
    
    allocate(gradientSNI(1:N_node,1:N_dmnsn))             ; gradientSNI     = 10.0d5
    allocate(gradient_anaSNI(1:N_node,1:N_dmnsn))         ; gradient_anaSNI = 10.0d5
    allocate(gradient_numSNI(1:N_node,1:N_dmnsn))         ; gradient_numSNI = 10.0d5    
    allocate(grd_mvSNI(1:N_mvCoordnte_withl0A0))          ; grd_mvSNI       = 10.0d5
    allocate(grdmv_xySNI(1:N_mvCoordnte))                 ; grdmv_xySNI     = 10.0d5
    
    allocate(grd_sprSNI(1:N_node,1:N_dmnsn))              ; grd_sprSNI      = -1.d30
    allocate(delE_delRSNI(1:N_spr))                       ; delE_delRSNI    = -1.d30
    allocate(delR_delNodeSNI(1:N_node,1:max_spr_node,1:2)); delR_delNodeSNI = -1.d30
    
    allocate(grd_areaSNI(1:N_node,1:N_dmnsn))             ; grd_areaSNI     = -1.d30
    allocate(delE_delASNI(1:N_cell))                      ; delE_delASNI    = -1.d30
    allocate(delA_delNodeSNI(1:N_node,1:max_area_node,1:2));delA_delNodeSNI = -1.d30
    
    allocate(grd_grvtnlSNI(1:N_node,1:N_dmnsn))                  ; grd_grvtnlSNI  = -1.0d30
    allocate(grd_bendSNI(1:N_node,1:N_dmnsn))                    ; grd_bendSNI    = -1.0d30
    allocate(delB_delCNSNI(1:N_node,1:max_area_node,1:N_dmnsn))  ; delB_delCNSNI  = -1.0d30
    allocate(delB_delNN1SNI(1:N_node,1:max_area_node,1:N_dmnsn)) ; delB_delNN1SNI = -1.0d30
    allocate(delB_delNN2SNI(1:N_node,1:max_area_node,1:N_dmnsn)) ; delB_delNN2SNI = -1.0d30
    
    
  end subroutine allocate_variabls_neededToStore_for_NIsys
  
  
  subroutine alloc_and_init_Arry_forContinued_NI_trnsfrm
    implicit none
    
    write(*,*)          modelID,"modelID in allocate_and_initialize_transfrm_variables"
    if (modelID .ne. 2) stop 'model ID must be 2 in Here in allocate_and_initialize_transfrm_variables'
    
    if (allctn_Continued_NI_trnsfrmArray==1) then ! ALLOCATION ALREADY HAPPENED
       continue
    elseif  (allctn_Continued_NI_trnsfrmArray==0) then
       
       max_node_sprS2   = 3  
       max_area_sprS2   = 2
       
       max_spr_nodeS2   = 6
       max_area_nodeS2  = 4
       
       max_node_areaS2  = 5 + (NAEC_Bsal) + (NAEC_Ltrl) + (NAEC_Ltrl)
       max_spr_areaS2   = 4 + (NAEC_Bsal) + (NAEC_Ltrl) + (NAEC_Ltrl)
       
       allocate(spr_nodeS2(1:N_spr,0:max_node_spr), spr_areaS2(1:N_spr,0:max_area_spr))
       allocate(node_sprS2(1:N_node,0:max_spr_node),node_areaS2(1:N_node,0:max_spr_node))
       allocate(area_sprS2(1:N_cell,0:max_spr_area),area_nodeS2(1:N_cell,0:max_node_area))
       allocate(NlistS2(1:N_node,1:max_area_node,1:2))
       
       spr_nodeS2  = -1 ; spr_areaS2  = -1
       node_sprS2  = -1 ; node_areaS2 = -1
       area_nodeS2 = -1 ; area_sprS2  = -1 
       NlistS2     = -1
       
       
       allocate(gradientS2(1:N_node,1:N_dmnsn))     ; gradientS2     = 10.0d5
       allocate(gradient_anaS2(1:N_node,1:N_dmnsn)) ; gradient_anaS2 = 10.0d5
       allocate(gradient_numS2(1:N_node,1:N_dmnsn)) ; gradient_numS2 = 10.0d5    
       allocate(grd_mvS2(1:N_mvCoordnte_withl0A0))  ; grd_mvS2       = 10.0d5
       allocate(grdmv_xyS2(1:N_mvCoordnte))         ; grdmv_xyS2     = 10.0d5
       
       allocate(grd_sprS2(1:N_node,1:N_dmnsn))               ; grd_sprS2      = -1.d30
       allocate(delE_delRS2(1:N_spr))                        ; delE_delRS2    = -1.d30
       allocate(delR_delNodeS2(1:N_node,1:max_spr_node,1:2)) ; delR_delNodeS2 = -1.d30
       
       allocate(grd_areaS2(1:N_node,1:N_dmnsn))                     ; grd_areaS2     = -1.d30
       allocate(delE_delAS2(1:N_cell))                              ; delE_delAS2    = -1.d30
       allocate(delA_delNodeS2(1:N_node,1:max_area_node,1:N_dmnsn)) ; delA_delNodeS2 = -1.d30
       
       allocate(grd_grvtnlS2(1:N_node,1:N_dmnsn))                  ; grd_grvtnlS2  = -1.0d30
       allocate(grd_bendS2(1:N_node,1:N_dmnsn))                    ; grd_bendS2    = -1.0d30
       allocate(delB_delCNS2(1:N_node,1:max_area_node,1:N_dmnsn))  ; delB_delCNS2  = -1.0d30
       allocate(delB_delNN1S2(1:N_node,1:max_area_node,1:N_dmnsn)) ; delB_delNN1S2 = -1.0d30
       allocate(delB_delNN2S2(1:N_node,1:max_area_node,1:N_dmnsn)) ; delB_delNN2S2 = -1.0d30
       
    endif
    
  end subroutine alloc_and_init_Arry_forContinued_NI_trnsfrm
  
  
  subroutine store_NI_SysVars_and_Arrays !sub1
    implicit none
    
    node_xyStrNI       = node_xy
    node_typStrNI      = node_typ   
    node_cnnctdStrNI   = node_cnnctd 
    double_nodeStrNI   = double_node
    count_this_dnStrNI = count_this_dn
    
    typ_sprStrNI = typ_spr
    k_sprStrNI   = k_spr   
    l0_StrNI     = l0 
    l_StrNI      = l
    
    nodePhi_typStrNI = nodePhi_typ
    k_phiStrNI       = k_phi
    CgXNode_StrNI    = CgXNode
    CgYNode_StrNI    = CgYNode
    
    N_nodeSNI = N_node
    N_sprSNI  = N_spr
    N_phiS    = N_phiSNI
    
  end subroutine store_NI_SysVars_and_Arrays
  
  subroutine get_TN_Sysvars !sub2
    implicit none
    
    write(*,*) N_node,N_spr,"inside get_TN1"
    
    N_node = N_node - (NAEC_Bsal+NAEC_Ltrl)*(Hlf_Ncell*2) - (NAEC_Bsal)*1
    N_spr  = N_spr  - (NAEC_Bsal+NAEC_Ltrl)*(Hlf_Ncell*2) - (NAEC_Bsal)*1
    
    write(*,*) N_node,N_spr,"inside get_TN2"
    
  end subroutine get_TN_Sysvars
  
  subroutine deallocate_and_reallocate_NI_Arrays !sub3
    implicit none
    
    deallocate(node_xy,node_typ,node_cnnctd,double_node,count_this_dn)
    deallocate(typ_spr,k_spr,l0,l)
    deallocate(nodePhi_typ,k_phi,CgXNode,CgYNode)
    
    call allocate_and_initialize_node_variables_wo_StrVars
    call allocate_and_initialize_spring_variables_wo_StrVars
    call allocate_and_initialize_grvVars_wo_StrVars
    call allocate_and_initialize_bend_variables_wo_StrVars
    
  end subroutine deallocate_and_reallocate_NI_Arrays
  
  
  subroutine get_TN_NodeVars ! sub4
    implicit none
    
    node_xy(1:N_node,1:N_dmnsn) = node_xyStrNI(1:N_node,1:N_dmnsn)
    node_typ(1:N_node)          = node_typStrNI(1:N_node)
    
    call get_list_of_double_nodes_method2
    call nodes_cnnctd_and_count_this_dn
    call print_NodeVars
    
  end subroutine get_TN_NodeVars
  
  
  subroutine get_TN_SprVars !sub5
    implicit none
    integer :: nsprsInACell
    integer :: cellNm,NAEC_val
    integer :: ApclSpNI(1:(NAEC_Apcl+1)),BsalSpNI(1:(NAEC_Bsal+1)),LtrlSpNI(1:(NAEC_Ltrl+1))
    integer :: ApclSpTN,BsalSpTN,LtrlSpTN
    real*8  :: lval=0.0d0
    integer :: sprTN,sprNI
    integer :: i1,i2,i3
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i1 = 1,Hlf_Ncell
       
       cellNm  = i1
       ApclSpTN=(cellNm-1)*3+1 ; BsalSpTN=(cellNm-1)*3+2 ; LtrlSpTN=(cellNm-1)*3+3
       call get_ApclBsalLtrl_OfNIsystm(cellNm,ApclSpNI,BsalSpNI,LtrlSpNI)
       
       do i2 = 1,3
          
          if (i2==1) then
             
             NAEC_val = NAEC_Apcl
             sprTN    = ApclSpTN
             sprNI    = ApclSpNI(1)
             
          elseif (i2==2) then
             
             NAEC_val = NAEC_Bsal
             sprTN    = BsalSpTN
             sprNI    = BsalSpNI(1)
             
          elseif (i2==3) then
             
             NAEC_val = NAEC_Ltrl
             sprTN    = LtrlSpTN
             sprNI    = LtrlSpNI(1)
             
          endif
          
          typ_spr(sprTN) = typ_spr(sprNI)
          k_spr(sprTN)   = k_spr(sprNI)/(NAEC_val+1)
          l0(sprTN)      = l0(sprNI)*(NAEC_val+1)
          
          lval = 0.0d0
          
          do i3 = 1,(NAEC_val+1)
             if (i2==1) lval = lval+l(ApclSpNI(i3))
             if (i2==2) lval = lval+l(BsalSpNI(i3))
             if (i2==3) lval = lval+l(LtrlSpNI(i3))
          enddo
          
          if (i2==1) l(ApclSpTN) = lval
          if (i2==2) l(BsalSpTN) = lval
          if (i2==3) l(LtrlSpTN) = lval
          
       enddo
       
       write(*,*) i1,ApclSpTN,BsalSpTN,LtrlSpTN,"TN spr"
       write(*,*) typ_spr(ApclSpTN),k_spr(ApclSpTN),l0(ApclSpTN),l(ApclSpTN),i1,"Apcl Prop"
       write(*,*) typ_spr(BsalSpTN),k_spr(BsalSpTN),l0(BsalSpTN),l(BsalSpTN),i1,"Bsal Prop"
       write(*,*) typ_spr(LtrlSpTN),k_spr(LtrlSpTN),l0(LtrlSpTN),l(LtrlSpTN),i1,"Ltrl Prop"
       
    enddo
    
  end subroutine get_TN_SprVars
  
  subroutine get_TN_kphi_and_nodePhiTyp
    implicit none
    integer :: i
    real*8  :: k_phiVal
    
    k_phiVal                       = k_phiStrNI(1,1)
    k_phi(1:N_node,1:max_Phi_node) = k_phiVal
    nodePhi_typ(1:N_node)          = 1
    
  end subroutine get_TN_kphi_and_nodePhiTyp
  
  
  subroutine get_TN_CgXvars
    implicit none
    
    CgXNode(1:N_node) = CgXNode_StrNI(1:N_node)
    
  end subroutine get_TN_CgXvars
  
  
  subroutine get_TN_CgYvars
    implicit none
    
    CgYNode(1:N_node) = CgYNode_StrNI(1:N_node)
    
  end subroutine get_TN_CgYvars
  
  
  
  subroutine get_ApclBsalLtrl_OfNIsystm(cellNm,ApclSp,BsalSp,LtrlSp)
    implicit none
    integer, intent(in)  :: cellNm
    integer, intent(out) :: ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1))
    
    integer :: i,j,jmax
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,3
       
       if (i==1) jmax=(NAEC_Apcl+1)
       if (i==2) jmax=(NAEC_Bsal+1)
       if (i==3) jmax=(NAEC_Ltrl+1)
       
       do j = 1,jmax
          if (i==1) ApclSp(j) = (cellNm-1)*(nsprsInACell) + j
          if (i==2) BsalSp(j) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
          if (i==3) LtrlSp(j) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
       enddo
       
    enddo
    
  end subroutine get_ApclBsalLtrl_OfNIsystm
  
  
  subroutine get_ApclBsalLtrl_OfNIsystmCrtclApclSrfc(cellNm,ApclSp,BsalSp,LtrlSp,NAEC_ApclV)
    implicit none
    integer, intent(in)                 :: cellNm
    integer, allocatable, intent(inout) :: ApclSp(:),BsalSp(:),LtrlSp(:)
    integer, intent(out)                :: NAEC_ApclV
    
    if (cellNm.gt.(Hlf_Ncell-NCP_CrtclApSrfc)) NAEC_ApclV = NAEC_ApclCrtcl
    if (cellNm.le.(Hlf_Ncell-NCP_CrtclApSrfc)) NAEC_ApclV = NAEC_Apcl
    
    write(*,*) NAEC_ApclCrtcl,NAEC_Apcl,NAEC_ApclV,"which values is picked"
    
    allocate(ApclSp(1:(NAEC_ApclV+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)))
    ApclSp = -1 ; BsalSp = -1 ; LtrlSp = -1
    
    if (NAEC_ApclV==NAEC_Apcl) then
       call get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase1(cellNm,ApclSp,BsalSp,LtrlSp)
    elseif (NAEC_ApclV==NAEC_ApclCrtcl) then
       call get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase2(cellNm,ApclSp,BsalSp,LtrlSp)
    endif
    
  end subroutine get_ApclBsalLtrl_OfNIsystmCrtclApclSrfc
  
  subroutine get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase1(cellNm,ApclSp,BsalSp,LtrlSp) ! Sprs BEFORE Crtcl ApSrfc being Introduced
    implicit none
    integer, intent(in)  :: cellNm
    integer, intent(out) :: ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)) 
    integer              :: i,j,jmax
    
    call get_ApclBsalLtrl_OfNIsystm(cellNm,ApclSp,BsalSp,LtrlSp) ! just renamed for clarity purpose
    
  end subroutine get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase1
  
  
  subroutine get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase2(cellNm,ApclSp,BsalSp,LtrlSp) ! Sprs AFTER Crtcl ApSrfc being Introduced
    implicit none
    integer, intent(in)  :: cellNm
    integer, intent(out) :: ApclSp(1:(NAEC_ApclCrtcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1))
    integer              :: i,j,jmax
    integer              :: nsprsInACellT1,nsprsInACellT2
    integer              :: sprsBFRtheCell,firstCrtclCell
    
    nsprsInACellT1 = (NAEC_Apcl+1)     +(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    nsprsInACellT2 = (NAEC_ApclCrtcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    Hlf_Ncell      = (N_cell-1)/2
    firstCrtclCell = Hlf_Ncell-NCP_CrtclApSrfc+1
    
    if (cellNm.lt.firstCrtclCell) then
       write(*,*) Hlf_Ncell,N_cell,cellNm,firstCrtclCell,"cellNm val"
       stop 'cellNm must be greater than firstCrtclCell: check prv val statemnt'
    endif
    
    sprsBFRtheCell = (cellNm-firstCrtclCell)*(nsprsInACellT2) + (firstCrtclCell-1)*(nsprsInACellT1)
    write(*,*) cellNm,firstCrtclCell,nsprsInACellT1,nsprsInACellT2,"chkBfr aprsBfrtheCell"
    write(*,*) sprsBFRtheCell,"sprsBFRtheCell"
    
    do i = 1,3
       
       if (i==1) jmax=(NAEC_ApclCrtcl+1)
       if (i==2) jmax=(NAEC_Bsal+1)
       if (i==3) jmax=(NAEC_Ltrl+1)
       
       do j = 1,jmax
          if (i==1) ApclSp(j) = sprsBFRtheCell + j
          if (i==2) BsalSp(j) = sprsBFRtheCell + (NAEC_ApclCrtcl+1) + j
          if (i==3) LtrlSp(j) = sprsBFRtheCell + (NAEC_ApclCrtcl+1) + (NAEC_Bsal+1) + j
       enddo
       
    enddo
    
  end subroutine get_ApclBsalLtrl_OfNIsystm_WT_CrtclSrfcCase2
  
  
  subroutine get_ApclBsal_singl_cell_NIsys_with_CM0(ApclSp,BsalSp)
    implicit none
    integer, intent(out) :: ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1))
    
    integer :: i,j,jmax
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    do i = 1,2
       
       if (i==1) jmax=(NAEC_Apcl+1)
       if (i==2) jmax=(NAEC_Bsal+1)
       
       do j = 1,jmax
          if (i==1) ApclSp(j) = (N_cell-1)*(nsprsInACell) + j
          if (i==2) BsalSp(j) = (N_cell-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
       enddo
       
    enddo
    
  end subroutine get_ApclBsal_singl_cell_NIsys_with_CM0
  
  
  subroutine get_ApclBsal_singl_cell_NIsys_with_CM0_and_VFarea(ApclSp,BsalSp)
    implicit none
    integer, intent(out) :: ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1))
    
    integer :: i,j,jmax
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    do i = 1,2
       
       if (i==1) jmax=(NAEC_Apcl+1)
       if (i==2) jmax=(NAEC_Bsal+1)
       
       do j = 1,jmax
          if (i==1) ApclSp(j) = (N_cell-2)*(nsprsInACell) + j
          if (i==2) BsalSp(j) = (N_cell-2)*(nsprsInACell) + (NAEC_Apcl+1) + j
       enddo
       
    enddo
    
  end subroutine get_ApclBsal_singl_cell_NIsys_with_CM0_and_VFarea
  
  subroutine get_ApclBsal_singl_cell_NIsys_with_CM0_VFarea_and_CrtclSrfc(ApclSp,BsalSp,NAEC_ApclV)
    implicit none
    integer, allocatable, intent(inout) :: ApclSp(:),BsalSp(:)
    integer,              intent(out)   :: NAEC_ApclV
    integer                             :: i,j,jmax
    integer                             :: nsprsInACellT1,nsprsInACellT2
    integer                             :: NsprBFRcrtclApSrfc,NsprAFTcrtclApSrfc
    
    allocate(ApclSp(1:(NAEC_ApclCrtcl+1)),BsalSp(1:(NAEC_Bsal+1)))
    
    nsprsInACellT1 = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    nsprsInACellT2 = (NAEC_ApclCrtcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    NsprBFRcrtclApSrfc = (N_cell-2-2*NCP_CrtclApSrfc)*(nsprsInACellT1)
    NsprAFTcrtclApSrfc = (2*NCP_CrtclApSrfc)*(nsprsInACellT2)
    
    write(*,*) nsprsInACellT1,nsprsInACellT2,"nsprsInACellT1-T2"
    write(*,*) NsprBFRcrtclApSrfc,NsprAFTcrtclApSrfc,"Nspr(BFR/AFT)crtclApSrfc"
    
    do i = 1,2
       
       if (i==1) jmax=(NAEC_ApclCrtcl+1)
       if (i==2) jmax=(NAEC_Bsal+1)
       
       do j = 1,jmax
          if (i==1) ApclSp(j) = NsprBFRcrtclApSrfc + NsprAFTcrtclApSrfc + j
          if (i==2) BsalSp(j) = NsprBFRcrtclApSrfc + NsprAFTcrtclApSrfc + (NAEC_ApclCrtcl+1) + j
       enddo
       
    enddo

    NAEC_ApclV = NAEC_ApclCrtcl
    
  end subroutine get_ApclBsal_singl_cell_NIsys_with_CM0_VFarea_and_CrtclSrfc
  
  subroutine get_Bsal_singl_cell_NIsys_with_CM_gt_0(BsalSp)
    implicit none
    integer, intent(out) :: BsalSp(1:(NAEC_Bsal+1))
    
    integer :: j,jmax
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    jmax=(NAEC_Bsal+1)
    
    do j = 1,jmax
       BsalSp(j) = (N_cell-1)*(nsprsInACell) + j
    enddo
    
  end subroutine get_Bsal_singl_cell_NIsys_with_CM_gt_0
  
  
  
end module switch_and_unswitch_models

