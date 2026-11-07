module node_info
  implicit none
  
  real*8 ,allocatable :: node_xy(:,:),     node_xyStr(:,:),     node_xyStrAN(:,:),     node_xyStrES(:,:)
  integer,allocatable :: node_typ(:),      node_typStr(:),      node_typStrAN(:),      node_typStrES(:)
  integer,allocatable :: node_cnnctd(:),   node_cnnctdStr(:),   node_cnnctdStrAN(:),   node_cnnctdStrES(:)
  integer,allocatable :: double_node(:,:), double_nodeStr(:,:), double_nodeStrAN(:,:), double_nodeStrES(:,:)
  integer,allocatable :: count_this_dn(:), count_this_dnStr(:), count_this_dnStrAN(:), count_this_dnStrES(:)
  
  real*8 ,allocatable :: node_xyStrNI(:,:)    , node_xyStrNVFR(:,:)
  integer,allocatable :: node_typStrNI(:)     , node_typStrNVFR(:)
  integer,allocatable :: node_cnnctdStrNI(:)  , node_cnnctdStrNVFR(:)
  integer,allocatable :: double_nodeStrNI(:,:), double_nodeStrNVFR(:,:)
  integer,allocatable :: count_this_dnStrNI(:), count_this_dnStrNVFR(:)
  
  ! NVFR = Non Vitelline Fluid Region Case
  
end module node_info

module spr_info
  implicit none
  
  integer,allocatable :: typ_spr(:),typ_sprStr(:),typ_sprStrAN(:),typ_sprStrES(:),typ_sprStrNI(:),typ_sprStrNVFR(:)
  real*8, allocatable :: k_spr(:),k_sprStr(:),k_sprStrAN(:),k_sprStrES(:),k_sprStrNI(:),k_sprStrNVFR(:)
  real*8, allocatable :: l0(:)   ,l0_Str(:)  ,l0_StrAN(:)  ,l0_StrES(:)  ,l0_StrNI(:)  ,l0_StrNVFR(:)
  real*8, allocatable :: l(:)    ,l_Str(:)   ,l_StrAN(:)   ,l_StrES(:)   ,l_StrNI(:)   ,l_StrNVFR(:) 
  
  real*8, allocatable :: Lt(:)   !Lt=target length
  real*8, allocatable :: alpha(:)
  integer,allocatable :: optmSpr(:)
  
end module spr_info

module area_info
  implicit none
  
  integer,allocatable :: typ_area(:)
  real*8, allocatable :: k_area(:), k_areaStr(:), k_areaStrNVFR(:)
  real*8, allocatable :: A0(:)    , A0_Str(:)   , A0_StrNVFR(:)   
  real*8, allocatable :: fctr(:)                      ! fctr = fctr A0 to be altered
  integer,allocatable :: area_sides(:)
  real*8, allocatable :: A(:),At(:),A_Str(:),A_StrNVFR(:)      ! At = target area
  real*8, allocatable :: beta(:)
  integer,allocatable :: optmCell(:)
  
  integer,allocatable :: lftSideCell(:),rghtSideCell(:)
  integer,allocatable :: cntrlCell(:),endRgnCell(:)
  
end module area_info

module bend_info
  implicit none
  
  integer, allocatable :: nodePhi_typ(:),      nodePhi_typStr(:),  nodePhi_typStrAN(:)
  integer, allocatable :: nodePhi_typStrES(:), nodePhi_typStrNI(:),nodePhi_typStrNVFR(:)
  real*8,  allocatable :: k_phi(:,:),          k_phiStr(:,:),      k_phiStrAN(:,:)
  real*8,  allocatable :: k_phiStrES(:,:),     k_phiStrNI(:,:),    k_phiStrNVFR(:,:)
  
end module bend_info


module diag_info
  implicit none
  
  real*8, allocatable :: D(:),Dt(:)
  real*8, allocatable :: gammaa(:)
  
end module diag_info

module grvtnl_info
  implicit none
  
  real*8               :: CgX,CgX_max,CgY,CgY_max
  real*8               :: CgX_strt,CgX_fnsh,CgY_strt,CgY_fnsh
  integer              :: N_CgXEval,N_CgYEval
  integer              :: N_region,N_imp_region
  integer              :: div_imp_region
  integer              :: prim_CgXNodes,addtnl_CgXNodes,prim_CgYNodes,addtnl_CgYNodes
  integer, allocatable :: imp_region(:)
  real*8,  allocatable :: CgX_val(:),CgY_val(:)
  integer              :: count_cgX,count_cgY
  real*8               :: region_incr,imp_region_incr
  real*8               :: imp_region_CgX_strt
  integer, allocatable :: grv_Node(:)
  logical              :: grvNode_lgcl
  
  real*8, allocatable  :: CgXNode(:),CgXNode_Str(:),CgXNode_StrAN(:),CgXNode_StrES(:),CgXNode_StrNI(:),CgXNode_StrNVFR(:)
  real*8, allocatable  :: CgYNode(:),CgYNode_Str(:),CgYNode_StrAN(:),CgYNode_StrES(:),CgYNode_StrNI(:),CgYNode_StrNVFR(:)
  
end module grvtnl_info

module SRyp_info
  implicit none
  integer             :: numRegulrNode,numRegulrNodeS
  integer             :: numInsrtdNode,numInsrtdNodeS
  
  integer             :: numApclNodes,numApclNodesS
  integer             :: numBsalNodes,numBsalNodesS
  integer             :: numLtrlNodes,numLtrlNodesS
  
  integer,allocatable :: apclNodes(:),bsalNodes(:),ltrlNodes(:)
  
  real*8, allocatable :: AmpApclYP(:),       AmpApclYP_Str(:)
  real*8, allocatable :: AlphaApclYP(:),     AlphaApclYP_Str(:)
  real*8, allocatable :: epsApclYP(:),       epsApclYP_Str(:)
  real*8, allocatable :: powrApclYP(:),      powrApclYP_Str(:)
  
  integer,allocatable :: activtnFctrApcl(:), activtnFctrApcl_Str(:)
  integer,allocatable :: activtnFctrBsal(:), activtnFctrBsal_Str(:)
  integer,allocatable :: activtnFctrLtrl(:), activtnFctrLtrl_Str(:)
  
end module SRyp_info

module curve_info
  implicit none
  
  real*8, allocatable :: curveX1Y1(:,:)
  real*8, allocatable :: curveX2Y2(:,:)
  real*8, allocatable :: P1(:),P2(:),dP(:)
  real*8, allocatable :: TC(:)
  real*8, allocatable :: R(:)
  real*8, allocatable :: xc_crv(:),yc_crv(:)
  real*8, allocatable :: AnglS(:),AnglE(:) !S=Start,E=End
  
  real*8,  allocatable :: curve_spr(:)
  integer, allocatable :: is_it_a_curve(:)
  
end module curve_info


module strct_info

  implicit none
  real*8, allocatable :: l0_strctTN(:,:),ks_strctTN(:,:)
  real*8, allocatable :: A0_strctTN(:,:),ka_strctTN(:,:)
  real*8, allocatable :: CgX_strctTN(:,:),CgY_strctTN(:,:)
  
  real*8, allocatable :: l0_strctNI(:,:),ks_strctNI(:,:)
  real*8, allocatable :: A0_strctNI(:,:),ka_strctNI(:,:)
  real*8, allocatable :: CgX_strctNI(:,:),CgY_strctNI(:,:)
  
  integer :: N_mvCoordnte1,  N_mvCoordnte2,  N_mvCoordnte3,  N_mvCoordnte4
  integer :: N_mvCoordnte5,  N_mvCoordnte6,  N_mvCoordnte7,  N_mvCoordnte8
  integer :: N_mvCoordnte9,  N_mvCoordnte10, N_mvCoordnte11, N_mvCoordnte12
  integer :: N_mvCoordnte13, N_mvCoordnte14, N_mvCoordnte15, N_mvCoordnte16
  integer :: N_mvCoordnte17, N_mvCoordnte18
  
  integer :: NMCWL0A01,  NMCWL0A02,  NMCWL0A03,  NMCWL0A04
  integer :: NMCWL0A05,  NMCWL0A06,  NMCWL0A07,  NMCWL0A08
  integer :: NMCWL0A09,  NMCWL0A010, NMCWL0A011, NMCWL0A012
  integer :: NMCWL0A013, NMCWL0A014, NMCWL0A015, NMCWL0A016
  integer :: NMCWL0A017, NMCWL0A018
  
  !NMCWL0A01 = N_mvCoordnte_withl0A0 for strct1
  
  real*8, allocatable :: coordntesXY_strct1(:), coordntesXY_strct2(:)
  real*8, allocatable :: coordntesXY_strct3(:), coordntesXY_strct4(:)
  real*8, allocatable :: coordntesXY_strct5(:), coordntesXY_strct6(:)
  real*8, allocatable :: coordntesXY_strct7(:), coordntesXY_strct8(:)
  real*8, allocatable :: coordntesXY_strct9(:), coordntesXY_strct10(:)
  real*8, allocatable :: coordntesXY_strct11(:),coordntesXY_strct12(:)
  real*8, allocatable :: coordntesXY_strct13(:),coordntesXY_strct14(:)
  real*8, allocatable :: coordntesXY_strct15(:),coordntesXY_strct16(:)
  real*8, allocatable :: coordntesXY_strct17(:),coordntesXY_strct18(:)
  
  real*8, allocatable :: coordntes_strct1(:), coordntes_strct2(:)
  real*8, allocatable :: coordntes_strct3(:), coordntes_strct4(:)
  real*8, allocatable :: coordntes_strct5(:), coordntes_strct6(:)
  real*8, allocatable :: coordntes_strct7(:), coordntes_strct8(:)
  real*8, allocatable :: coordntes_strct9(:), coordntes_strct10(:)
  real*8, allocatable :: coordntes_strct11(:),coordntes_strct12(:)
  real*8, allocatable :: coordntes_strct13(:),coordntes_strct14(:)
  real*8, allocatable :: coordntes_strct15(:),coordntes_strct16(:)
  real*8, allocatable :: coordntes_strct17(:),coordntes_strct18(:)
  
  real*8, allocatable :: nodeXY_strctTN(:,:,:),nodeXY_strctNI(:,:,:)
  
end module strct_info


module cell_info
  use node_info
  use spr_info
  use area_info
  use bend_info
  use SRyp_info
  use diag_info
  use grvtnl_info
  use curve_info
  use strct_info
end module cell_info


module node_end_info
  implicit none
  integer :: max_endNode
  
  integer, allocatable :: lft_endNode(:),rght_endNode(:)
  integer, allocatable :: clft_top_endNode(:),additnl_node_endNode(:),crght_top_endNode(:)
  integer, allocatable :: clft_bot_endNode(:),crght_bot_endNode(:)
  
  real*8, allocatable  :: lft_nodeXY(:,:), rght_nodeXY(:,:) !bndry_nodes_xy(not sure if needed for later)
  
end module node_end_info

module spr_end_info
  implicit none
  integer :: max_endSpring
  
  integer, allocatable :: lft_endSpring(:),rght_endSpring(:),cntrlCell_Spring(:)
  integer, allocatable :: clft_top_endSpring(:),additnl_node_endSpring(:)
  integer, allocatable :: crght_top_endSpring(:)
  integer, allocatable :: clft_bot_endSpring(:),crght_bot_endSpring(:)
  integer, allocatable :: endCell_Spring(:)
  
end module spr_end_info

module area_end_info
  implicit none
  integer :: max_endArea
  
  integer, allocatable :: lft_endArea(:), rght_endArea(:)
  integer, allocatable :: clft_top_endArea(:),additnl_node_endArea(:),&
       crght_top_endArea(:)
  integer, allocatable :: clft_bot_endArea(:),crght_bot_endArea(:)
  
end module area_end_info


module end_info
  use node_end_info
  use spr_end_info
  use area_end_info
end module end_info


module non_changing_parameters
  implicit none
  integer,parameter :: N_dmnsn = 2
end module non_changing_parameters


module system_parameters
  use sys_building_info
  use cell_info
  use end_info
  use strct_info
  use non_changing_parameters
  
  implicit none
  
  integer :: modelID,switchBackFrmNItoTN
  integer :: AddedCellModelInitiation
  integer :: N_struct,strctNo
  integer :: NAEC                          ! NAEC=nodes added in each curve
  integer :: NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,NAEC_ApclCrtcl
  integer :: numINperCell                  ! num of Inserted Node Per Cell
  
  real*8  :: AR_lft_rght,AR_cntrl !AR=Aspect Ratio(Lb/La)
  real*8  :: La1,Lb1,La2,Lb2
  real*8  :: bndryCell_xdis,howlong_bndryCell
  integer :: N_SysBlck
  
  logical, allocatable :: build_sys(:)
  
  integer :: sys_order ! =0(lft_frst);=1(rght_frst)
  logical :: lgcl_llp,lgcl_rlp,lgcl_cltp,lgcl_anp,lgcl_crtp,lgcl_clbp,lgcl_crbp
  
  integer :: ncl,ncr,ncc,ncclt,nccrt,ncclb,nccrb,ncer!nccrb=Noofcellcntrlrghtbot
  integer :: nsl,nsr,nsc, nsclt,nscrt,nsclb,nscrb,nser
  integer :: ndl,ndr,ndc, ndclt,ndcrt,ndclb,ndcrb,nder
  
  integer :: nsecl ,nsecr, nsecc, nsecclt ,nseccrt ,nsecclb ,nseccrb, nsecer
  integer :: nseclS,nsecrS,nseccS,nseccltS,nseccrtS,nsecclbS,nseccrbS,nsecerS !S=Store
  
  !nseccrb = no of springs each cell, central right bottom
  integer :: ndecl,ndecr, ndecclt,ndeccrt,ndecclb,ndeccrb
  
  integer :: nvsl,nvsr, nvsclt,nvscrt, nvsclb,nvscrb
  !nvscrb = no of verical sides central right bottom
  
  integer :: N_node,lft_node,rght_node,cntrl_node,endRgn_node,blck_node
  integer :: clft_top_node,crght_top_node,additnl_node,clft_bot_node,crght_bot_node
  
  integer :: lft_end_node,rght_end_node,clft_top_end_node,additnl_end_node
  integer :: crght_top_end_node,clft_bot_end_node,crght_bot_end_node
  integer :: N_spr,blck_spr,clft_top_spr,additnl_spr,crght_top_spr,N_nospr,clft_bot_spr,crght_bot_spr
  
  integer :: lft_end_spr,rght_end_spr,clft_top_end_spr
  integer :: additnl_end_spr,crght_top_end_spr,clft_bot_end_spr,crght_bot_end_spr

  integer :: first_spr_of_clft_bot,first_spr_of_crght_bot
  integer :: N_diag,clft_top_diag,crght_top_diag,clft_bot_diag,crght_bot_diag
  
  integer :: N_phi
  integer :: N_nodeS    ,N_sprS    ,N_phiS
  integer :: N_nodeSNI  ,N_sprSNI  ,N_phiSNI
  integer :: N_nodeSAN  ,N_sprSAN  ,N_phiSAN
  integer :: N_nodeSES  ,N_sprSES  ,N_phiSES
  integer :: N_nodeSNVFR,N_sprSNVFR,N_phiSNVFR,N_cellS,N_cellSNVFR
  
  real*8  :: x_wall,Lx
  real*8  :: unt_len, unt_len_cb, cntrl_gap
  real*8  :: unt_lenX,unt_lenY
  real*8  :: deflc_frm_pulley,deflc_frm_cntrl_top_blck,deflc_for_endRgn
  real*8  :: initial_l0_ofCP
  real*8  :: origin(1:N_dmnsn)
  real*8  :: Pulley_dis(1:N_dmnsn),Origin_dis(1:N_dmnsn)
  
  integer :: N_cell,Hlf_Ncell,blck_cell,cntrl_cell
  integer :: N_lftCells,N_rghtCells,N_cntrlCells,N_endRgnCells
  
  integer :: cntrl_cell_nmbr,clft_top_cell_nmbr,crght_top_cell_nmbr
  integer :: clft_bot_cell_nmbr,crght_bot_cell_nmbr,endRgn_cell_nmbr
  integer :: nnecl,nnecr,nnecc, nnecclt,nneccrt, nnecclb,nneccrb, nnecer 
  !nneccrb = no of nodes each cell central right bottom
  
  integer :: N_doubleNode
  integer :: count_nodes!,count_dn
  integer :: count_spr,count_area,count_diag
  logical :: findn_DoubleNode(1:2) 
  
  integer :: N_mvCoordnte,N_mvCoordnte_withl0A0,N_variabls
  integer :: N_mvCoordnteS,N_mvCoordnte_withl0A0_S
  integer :: N_mvCoordnteSNI,N_mvCoordnte_withl0A0_SNI
  
  logical :: spr_lgcl,area_lgcl,grvtnl_lgcl,bend_lgcl,SRyp_lgcl
  logical :: spr_lgcl_cnstr,area_lgcl_cnstr,grvtnl_lgcl_cnstr,bend_lgcl_cnstr,SRyp_lgcl_cnstr
  
  logical :: l0_variatn_lgcl,A0_variatn_lgcl!if l0,A0 would also be considred as variables
  integer :: select_xy !for grvtnl Enrgy
  integer :: funcChoice
  
  integer :: N_optmSpr,N_optmCell
  integer :: max_Phi_node !maxm_Phi(Angle)_connectd_to_a_Node
  
  integer :: N_thrshAngl
  integer :: N_curve
  real*8  :: x_shft,y_shft
  
  integer :: InitiatorCell
  integer :: filechk

  integer :: CyclNo
  integer :: nclMin,ncrMin
  integer :: N_actvCell,N_actvSpr,N_actvNode

  integer :: EquilAlgrthm
  integer :: grdMismatch
  integer :: dmlshDecsn
  integer :: MeetData,SaveOrReadData
  integer :: N_seq
  integer :: CellsMeet
  integer :: nodeInsrtnStrts,AC_NACswitchback
  integer :: hlfCycl ! this is to find out either in first/second part of cycle 
  
  character(len=100) :: directryNm1,directryNm2,directryFlnm
  
  integer :: whatsDirectryNm1
  integer :: SystemTyp
  integer :: RUN_WTorWO_Force
  integer :: EnergyChk=0
  integer :: printIn_sprE=0
  
  real*8  :: restrctring_xval
  integer :: VF_regionModelled=0,addedNCPair=0
  integer :: NI_incldAS_woPP=0    ! NodeInsrtd inclding Apcl Side without PulleyPnt
  integer :: NI_AS_wtCrtclSurfc=0
  integer :: IN_apclSideCaseNo=0
  integer :: NCP_CrtclApSrfc=0,addApNodesPerApSpr=0,addApNodesPerSgmntdApSpr=0
  
  ! NCP_CrtclApSrfc ~Number of CellPairs with Critical Apical Surface
  ! addApNodesPerApSpr ~ additional apical nodes per apical side for CRTCL surfc
  
contains
  
  subroutine initialize_system_params
    implicit none
    
    modelID                  = 1 ! 01=Terminal node model, 02=Node insertion model
    switchBackFrmNItoTN      = 0
    AddedCellModelInitiation = 0
    
    N_struct  = 3 ; strctNo   = 0
    NAEC      = 0
    NAEC_Apcl = 0 ; NAEC_Bsal = 0 ; NAEC_Ltrl = 0 ; NAEC_ApclCrtcl = 0
    
    
    AR_lft_rght = 10000.0d0
    AR_cntrl    = 10000.0d0
    
    La1 = -1000.0d0 ; Lb1 = -1000.0d0
    La2 = -1000.0d0 ; Lb2 = -1000.0d0
    
    bndryCell_xdis    = -1000.0d0
    howlong_bndryCell = 100.0d0
    
    N_SysBlck = 0
    
    ncl=0; ncr=0; ncc=0; ncclt=0; nccrt=0; ncclb=0; nccrb=0
    nsl=0; nsr=0; nsc=0; nsclt=0; nscrt=0; nsclb=0; nscrb=0
    
    nnecl=0;nnecr=0; nnecc=0;nnecclt=0;nneccrt=0; nnecclb=0;nneccrb=0 ; nnecer=0  
    nsecl=0;nsecr=0; nsecc=0;nsecclt=0;nseccrt=0; nsecclb=0;nseccrb=0 ; nsecer=0 
    
    nvsl=0; nvsr=0; nvsclt=0; nvscrt=0; nvsclb=0; nvscrb=0
        
    blck_node     = 0; clft_top_node  = 0
    additnl_node  = 0; crght_top_node = 0
    clft_bot_node = 0; crght_bot_node  = 0
    
    lft_end_node       = 0; rght_end_node    = 0    
    clft_top_end_node  = 0; additnl_end_node = 0; crght_top_end_node = 0
    clft_bot_end_node  = 0; crght_bot_end_node = 0
    
    clft_top_spr = 0; additnl_spr = 0; crght_top_spr = 0
    clft_bot_spr = 0; crght_bot_spr = 0
    
    first_spr_of_clft_bot = 0 ; first_spr_of_crght_bot = 0
    
    x_wall     = -1.0d30
    unt_len    = -1.0d30
    unt_len_cb = -1.0d30
    unt_lenX   = -1.0d30
    unt_lenY   = -1.0d30
    cntrl_gap  = -1.0d30
    Lx         = -1.0d30
    
    deflc_frm_pulley = -1.0d30
    deflc_frm_cntrl_top_blck = -1.0d30
    deflc_for_endRgn = -1.0d30
    
    N_doubleNode = 0
    
    count_nodes = 0
    !count_dn    = 1
    count_spr   = 0
    count_diag  = 0
    count_area  = 0
    
    findn_DoubleNode(1:2) = .False.
    
    N_node = 0
    N_spr  = 0 ; N_nospr   = 0
    N_cell = 0 ; Hlf_Ncell = 0
    N_lftCells = 0 ; N_rghtCells = 0 ; N_cntrlCells = 0 ; N_endRgnCells = 0
    N_phi = 0
    
    N_doubleNode = 0
    
    N_mvCoordnte          = 0 ; N_mvCoordnteS           = 0
    N_mvCoordnte_withl0A0 = 0 ; N_mvCoordnte_withl0A0_S = 0
    N_variabls            = 0
    
    spr_lgcl    = .False.
    area_lgcl   = .False.
    grvtnl_lgcl = .False.
    bend_lgcl   = .False.
    SRyp_lgcl   = .False.
    
    spr_lgcl_cnstr    = .False.
    area_lgcl_cnstr   = .False.
    grvtnl_lgcl_cnstr = .False.
    bend_lgcl_cnstr   = .False.
    SRyp_lgcl_cnstr   = .False.
    
    l0_variatn_lgcl = .False.
    A0_variatn_lgcl = .False.
    
    select_xy  = 0
    funcChoice = 0
    
    max_Phi_node = 0
    N_thrshAngl  = 0
    N_curve      = 0

    x_shft = 0.01d0 ; y_shft = 0.01d0
    
    Pulley_dis(1:2) = -1000.0d0 ; Origin_dis(1:2) = -1000.0d0
    Origin(1:2)     = -1000.0d0
    
    InitiatorCell = 0
    filechk = 0
    
    CyclNo = 0
    nclMin = 0 ; ncrMin = 0
    N_actvCell = 0 ; N_actvSpr = 0
    
    EquilAlgrthm = 0
    grdMismatch  = 0
    dmlshDecsn   = -1 
    MeetData     = 0 ; SaveOrReadData = 0 
    
    N_seq = 0
    CellsMeet = 0
    hlfCycl = 0
    
    whatsDirectryNm1 = 1
    
    if (whatsDirectryNm1==1) then
       write(directryNm1,*) "/home"
    elseif (whatsDirectryNm1==2) then
       write(directryNm1,*)"/LOCAL_quagga"
    endif
    
    write(directryNm2,*)"/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/"
    write(directryFlnm,*)trim(adjustl(directryNm1))//trim(adjustl(directryNm2))

    
    AC_NACswitchback = 0
    SystemTyp        = 0 ! Normal Case,Bfr adding 2 Node in Cntrl Cell Ap Membrn
    RUN_WTorWO_Force = 0 ! 0 = no force pulling down, 1 = force pulling down
    
  end subroutine initialize_system_params
  
  subroutine reinitialize_system_params
    implicit none
    
    AR_lft_rght = 10000.0d0
    AR_cntrl    = 10000.0d0
    
    La1 = -1000.0d0 ; Lb1 = -1000.0d0
    La2 = -1000.0d0 ; Lb2 = -1000.0d0
    
    bndryCell_xdis    = -1000.0d0
    howlong_bndryCell = 100.0d0
    
    N_SysBlck = 0
    
    ncl=0; ncr=0; ncc=0; ncclt=0; nccrt=0; ncclb=0; nccrb=0
    nsl=0; nsr=0; nsc=0; nsclt=0; nscrt=0; nsclb=0; nscrb=0
    
    nnecl=0; nnecr=0; nnecc=0; nnecclt=0; nneccrt=0; nnecclb=0; nneccrb=0 ; nnecer=0
    nsecl=0; nsecr=0; nsecc=0; nsecclt=0; nseccrt=0; nsecclb=0; nseccrb=0 ; nsecer=0 
    
    nvsl=0; nvsr=0; nvsclt=0; nvscrt=0; nvsclb=0; nvscrb=0
        
    blck_node     = 0; clft_top_node  = 0
    additnl_node  = 0; crght_top_node = 0
    clft_bot_node = 0; crght_bot_node  = 0
    
    lft_end_node       = 0; rght_end_node    = 0    
    clft_top_end_node  = 0; additnl_end_node = 0; crght_top_end_node = 0
    clft_bot_end_node  = 0; crght_bot_end_node = 0
    
    clft_top_spr = 0; additnl_spr = 0; crght_top_spr = 0
    clft_bot_spr = 0; crght_bot_spr = 0
    
    first_spr_of_clft_bot = 0 ; first_spr_of_crght_bot = 0
    
    x_wall     = -1.0d30
    unt_len    = -1.0d30
    unt_len_cb = -1.0d30
    unt_lenX   = -1.0d30
    unt_lenY   = -1.0d30
    cntrl_gap  = -1.0d30
    Lx         = -1.0d30
    
    deflc_frm_pulley = -1.0d30
    deflc_frm_cntrl_top_blck = -1.0d30
    deflc_for_endRgn = -1.0d30
    
    N_doubleNode     = 0
    
    count_nodes = 0
    count_spr   = 0
    count_diag  = 0
    count_area  = 0
    
    findn_DoubleNode(1:2) = .False.

    N_node = 0
    N_spr  = 0 ; N_nospr   = 0
    N_cell = 0 ; Hlf_Ncell = 0
    N_lftCells = 0 ; N_rghtCells = 0 ; N_cntrlCells = 0 ; N_endRgnCells = 0
    N_phi = 0
    
    N_doubleNode = 0
    
    N_mvCoordnte          = 0 ; N_mvCoordnteS           = 0
    N_mvCoordnte_withl0A0 = 0 ; N_mvCoordnte_withl0A0_S = 0
    N_variabls            = 0
    
    select_xy  = 0
    funcChoice = 0
    
    N_thrshAngl = 0
    N_curve     = 0
    
  end subroutine reinitialize_system_params

  
  
  subroutine giving_big_values_to_global_params!Not needed right now
    implicit none
    integer :: sparse_value

    sparse_value = 100
    
    N_node = sparse_value
    N_spr  = sparse_value
    N_cell = sparse_value
    
    N_doubleNode = sparse_value
    N_mvCoordnte = sparse_value
    
  end subroutine giving_big_values_to_global_params
  
  subroutine get_sys_lgcls(sprLgclInput,areaLgclInput,grvLgclInput,&
       bendLgclInput,SRypLgclInput)
    implicit none
    logical :: sprLgclInput,areaLgclInput,grvLgclInput
    logical :: bendLgclInput,SRypLgclInput
    
    spr_lgcl    = sprLgclInput
    area_lgcl   = areaLgclInput
    grvtnl_lgcl = grvLgclInput
    bend_lgcl   = bendLgclInput
    SRyp_lgcl   = SRypLgclInput
    
  end subroutine get_sys_lgcls
  
  subroutine get_sys_lgcls_cnstr(sprLgclInput,areaLgclInput,grvLgclInput,&
       bendLgclInput,SRypLgclInput)
    implicit none
    logical :: sprLgclInput,areaLgclInput,grvLgclInput
    logical :: bendLgclInput,SRypLgclInput
    
    spr_lgcl_cnstr    = sprLgclInput
    area_lgcl_cnstr   = areaLgclInput
    grvtnl_lgcl_cnstr = grvLgclInput
    bend_lgcl_cnstr   = bendLgclInput
    SRyp_lgcl_cnstr   = SRypLgclInput
    
  end subroutine get_sys_lgcls_cnstr
  
  subroutine get_l0A0_varitn_lgcls(l0_lgcl,A0_lgcl)
    implicit none
    logical :: l0_lgcl
    logical :: A0_lgcl

    l0_variatn_lgcl = l0_lgcl
    A0_variatn_lgcl = A0_lgcl
    
  end subroutine get_l0A0_varitn_lgcls

  subroutine get_AspectRatio_LaLb
    implicit none

    AR_lft_rght = 5.00d0
    AR_cntrl    = 2.00d0
    
    La1 = 1.50d0
    Lb1 = AR_lft_rght * La1
    
    La2 = 2.50d0
    Lb2 = AR_cntrl * La2

    howlong_bndryCell = 1.0d0
    
    !if (strctNo==1) then
     !  bndryCell_xdis    = howlong_bndryCell * La1
    !elseif (strctNo==2) then
       !!bndryCell_xdis    = (howlong_bndryCell+1.0d0) * La1
       !bndryCell_xdis    = howlong_bndryCell * La1
    !elseif (strctNo==4) then
     !  bndryCell_xdis    = howlong_bndryCell * La1
    !endif
    
    bndryCell_xdis    = howlong_bndryCell * La1

  end subroutine get_AspectRatio_LaLb
  
  subroutine lft_ladder_params
    implicit none
    
    if (stageNo.ne.4) then
       ncl = 11 !it was 8
       
    elseif (stageNo==4) then
       if (strctNo==1) ncl = 4  !Num of cells lft !_NCL !(prv it was 4)
       if (strctNo==2) ncl = 3
       if (strctNo==4) ncl = 5 !(For cyclic purpose) 
    endif
    
    nsecl = 3  !Num of spr in each cell lft/spr IDs
    ndecl = 2  !Num of diag in each cell lft/diag IDs
    nnecl = 4  !Num of nodes each cell lft/node IDs
    
    nsl  = nsecl * ncl  !nsl = Num of spr left
    !nCRl = nsecl * ncl  !nCRl = Num of crv left
    ndl  = ndecl * ncl  !ndl = Num of diag left
    nvsl = ncl + 1  !nvsl = Num of vrtcl side lft
    
    lft_node = 2*nvsl
     
    !write(*,*) lft_end_node,lft_end_spr,"LEN,LES"
    
  end subroutine lft_ladder_params
  
  
  subroutine tst_lft_ladder_params
    implicit none
    
    if (strctNo==1) ncl = 4
    if (strctNo==2) ncl = 2
    
    nsecl = 3
    ndecl = 2
    nnecl = 4
    
    nsl  = nsecl * ncl
    ndl  = ndecl * ncl
    nvsl = ncl + 1 
    
    lft_node = 2*nvsl
    
    
  end subroutine tst_lft_ladder_params
  
  
  subroutine rght_ladder_params
    implicit none
    
    if (stageNo.ne.4) then
       ncr = 11 ! it was 8
       
    elseif (stageNo==4) then
       if (strctNo==1) ncr = 4 !Num of cells rght !_NCR !(prv it was 4)
       if (strctNo==2) ncr = 3
       if (strctNo==4) ncr = 5 !(For cyclic purpose)
    endif
    
    nsecr = 3
    ndecr = 2
    nnecr = 4
    
    nsr  = nsecr * ncr
    !nCRr = nsecr * ncr
    ndr  = ndecr * ncr
    nvsr = ncr + 1
    
    rght_node = 2*nvsr
      
  end subroutine rght_ladder_params
  
  subroutine cntrlCell_params
    implicit none
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "cntrlCell_params should not be executed"
       stop
    endif
    
    ncc   = 1
    nsecc = 2
    nnecc = 4
    
    nsc = nsecc * ncc
    cntrl_node = 0
    cntrl_cell_nmbr = ncl+ncr+1
    
  end subroutine cntrlCell_params
  
  subroutine clft_top_params
    implicit none
      
    ncclt   = 1 !Num of cells cntrl lft top
    nsecclt = 3 !Num of spr in each cell cntrl lft top
    ndecclt = 5 !Num of dig in each cell cntrl lft top 5c2-5 
    nnecclt = 5 !Num of nodes in each cell cntrl lft top
    
    nsclt   = nsecclt * ncclt !nsclt = Num of spr cntrl lft top
    !nCRclt  = 0 ! Num of crv cntrl lft top 
    ndclt   = ndecclt * ncclt
    
    nvsclt  = 1  !nvscl = Num of vrtcl (littlebit titled may be) side cntrl lft top (not_needed right now)
    
    clft_top_node      = 3*ncclt
    clft_top_end_node  = lft_node + rght_node + clft_top_node
    
    clft_top_spr       = nsclt
    clft_top_diag      = ndclt
    clft_top_end_spr   = nsl + nsr + clft_top_spr
    clft_top_cell_nmbr = ncl + ncr + ncclt
    
  end subroutine clft_top_params
  
  subroutine additnl_node_params
    implicit none
    
    additnl_node = 1
    additnl_end_node = lft_node + rght_node + clft_top_node + additnl_node
    
    additnl_spr  = 1
    additnl_end_spr = nsl + nsr + clft_top_spr + additnl_spr
    
  end subroutine additnl_node_params
  
  subroutine crght_top_params
    implicit none
    
    nccrt   = 1
    nseccrt = 3
    ndeccrt = 5
    nneccrt = 5 
    
    nscrt   = nseccrt * nccrt
    !nCRcrt  = 0 !Num of crv central rght top
    ndcrt   = ndeccrt * nccrt
    nvscrt  = 1 
      
    crght_top_node = 3*nccrt
    !write(*,*) crght_top_node,"inside cntrl_rght_top_params"
    
    crght_top_end_node = lft_node + rght_node + clft_top_node &
         + additnl_node + crght_top_node
      
    crght_top_spr     = nscrt
    crght_top_diag    = ndcrt
    
    crght_top_end_spr = nsl + nsr + clft_top_spr + additnl_spr + crght_top_spr
    
    crght_top_cell_nmbr = ncl + ncr + ncclt + nccrt
    
  end subroutine crght_top_params
  
  
  subroutine clft_bot_params
    implicit none
    
    if (strctNo==1.or.strctNo==4)  then
       ncclb = 3 !Num of cells cntrl lft bottom
    elseif (strctNo==2) then
       ncclb = 4
    endif
    
    nsecclb = 3 !Num of spr in each cell cntrl lft bottom
    ndecclb = 2 !Num of diag in each cell cntrl lft bottom
    nnecclb = 4 !Num of nodes in each cell cntrl lft bottom(in top it is 5)
    
    nsclb   = nsecclb * ncclb !nsclb = Num of spr cntrl lft bottom
    !nCRclb  = nsecclb * ncclb !nCRclb = Num of crv cntrl lft bottom
    ndclb   = ndecclb * ncclb
    nvsclb  = ncclb  !nvsclb = Num of vrtcl (littlebit titled may be) side cntrl lft bottom (not_needed right now)
    
    clft_bot_node      = 2*ncclb
    clft_bot_end_node  = lft_node + rght_node + clft_top_node &
         + additnl_node + crght_top_node + clft_bot_node
    
    clft_bot_spr       = nsclb
    clft_bot_diag      = ndclb
    clft_bot_end_spr   = nsl + nsr + clft_top_spr + additnl_spr &
         + crght_top_spr + clft_bot_spr
    
    first_spr_of_clft_bot = nsl+ nsr+ clft_top_spr+ crght_top_spr+ 1 
    clft_bot_cell_nmbr    = ncl + ncr + ncclt + nccrt + 2*ncclb - 1!2 for smmtry
    
  end subroutine clft_bot_params
  
  
  subroutine crght_bot_params
    implicit none
    
    if (strctNo==1.or.strctNo==4) then
       nccrb = 3
    elseif (strctNo==2) then
       nccrb = 4
    endif
    
    nseccrb = 3
    ndeccrb = 2
    nneccrb = 4 !in top it is 5
    
    nscrb   = nseccrb * nccrb
    !nCRcrb  = nseccrb * nccrb !num of crv cntrl rght bottom
    ndcrb   = ndeccrb * nccrb
    nvscrb  = nccrb 
    
    crght_bot_node = 2*nccrb
    !write(*,*) crght_bot_node,"inside cntrl_rght_bot_params"
    
    crght_bot_end_node = lft_node + rght_node + clft_top_node &
         + additnl_node + crght_top_node + clft_bot_node + crght_bot_node
      
    crght_bot_spr     = nscrb
    crght_bot_diag    = ndcrb
    crght_bot_end_spr = nsl + nsr + clft_top_spr + additnl_spr &
         + crght_top_spr + clft_bot_spr + crght_bot_spr
    
    first_spr_of_crght_bot = nsl+ nsr+ clft_top_spr+ crght_top_spr+ nsecclb+ 1  
    crght_bot_cell_nmbr    = ncl + ncr + ncclt + nccrt + ncclb + nccrb
    
  end subroutine crght_bot_params
  
  subroutine endCell_params
    implicit none
    
    if (stageNo==4 .and. stageType==1) then
       continue
    else
       write(*,*) "endCell_params should not be executed"
       stop
    endif
    
    ncer   = 1 !no of cells end region
    nsecer = 1 !no of sprs each cell end region
    nnecer = 3 !no of nodes each cell end region
    
    nser = ncer*nsecer !no of sprs end region
    endRgn_node = 0 !due to no additnl nodes
    endRgn_cell_nmbr = ncl + ncr + ncclt + nccrt + ncclb + nccrb + 1
    
  end subroutine endCell_params
  
  subroutine lft_rght_common_params
    implicit none
    
    blck_node = lft_node + rght_node  
    blck_spr  = nsl + nsr
    blck_cell = ncl + ncr
    
    write(*,*) blck_node,blck_spr,blck_cell,"inside_lft_rght_common"
    
  end subroutine lft_rght_common_params
  
  
  subroutine global_params
    implicit none
    
    write(*,*)  lft_node,rght_node,clft_top_node,additnl_node,&
         crght_top_node,clft_bot_node,crght_bot_node,endRgn_node,"bfr_N_node calc"
    
    N_node  = lft_node + rght_node + cntrl_node + clft_top_node + additnl_node &
         + crght_top_node + clft_bot_node + crght_bot_node + endRgn_node
    
    N_spr   = nsl + nsr + nsc + clft_top_spr + additnl_spr + crght_top_spr &
         + clft_bot_spr + crght_bot_spr + nser
    N_diag  = ndl + ndr + ndc + clft_top_diag + crght_top_diag &
         + clft_bot_diag + crght_bot_diag
    
    N_nospr = 2
    N_cell  = ncl + ncr + ncc + ncclt + nccrt + ncclb + nccrb + ncer
    Hlf_Ncell = N_cell/2
    
    N_lftCells    = ncl + ncclt + ncclb
    N_rghtCells   = ncr + nccrt + nccrb
    N_cntrlCells  = ncc
    N_endRgnCells = ncer
    
    write(*,*) N_node,N_spr,N_cell,"N_node,N_spr,N_cell"
    write(*,*) N_lftCells,N_rghtCells,N_cntrlCells,N_endRgnCells,"N lft,rght & midl Cells"
    
    N_thrshAngl = ncl+ncr
    !N_curve     = nCRl + nCRr + nCRclt + nCRcrt + nCRclb + nCRcrb
    !stop
    
  end subroutine global_params
  
  
  subroutine wall_cell_lngth_and_other_params
    use sys_building_info
    implicit none
    
    x_wall     = 0.0d0
    unt_len    = 1.0d0
    unt_lenX   = 1.0d0
    unt_lenY   = 1.0d0
    unt_len_cb = 0.50d0
    
    if ((clft_top_toBe_built.eqv..True.) .and. (crght_top_toBe_built.eqv..True.)) then
       cntrl_gap = 2.0d0
    elseif ((clft_top_toBe_built.eqv..True.) .or. (crght_top_toBe_built.eqv..True.))then
       cntrl_gap = 1.0d0
    else                       !both_are_false
       cntrl_gap = 0.0d0
    endif
    
    if (stageNo==1) then
       if (stageType==1) Lx = 2.0d0*bndryCell_xdis+(nvsl-2)*La1+(nvsr-2)*La1+La1
       if (stageType==2) Lx = 2.0d0*bndryCell_xdis+(nvsl-2)*La1+(nvsr-2)*La1
       
    elseif (stageNo==2) then
       write(*,*) "flnm:info_modules.f08,sb:wall_cell_lngth"
    elseif (stageNo==3) then
       write(*,*) "flnm:info_modules.f08,sb:wall_cell_lngth"
       
    elseif (stageNo==4) then
       
       if (stageType==1 .or. stageType==2) then
          
          if (nvsl.gt.0 .and. nvsr.gt.0) then
             Lx = 2.0d0*bndryCell_xdis + (nvsl-3)*La1 + (nvsr-3)*La1 + (cntrl_gap*Lb2)
          elseif (nvsl.gt.0 .and. nvsr.eq.0) then
             Lx = bndryCell_xdis + (nvsl-3)*La1 + (cntrl_gap*Lb2)
          elseif (nvsl.eq.0 .and. nvsr.gt.0) then
             Lx = bndryCell_xdis + (nvsr-3)*La1 + (cntrl_gap*Lb2)
          else
             Lx = (cntrl_gap*Lb2)
          endif
          
       endif
       
    endif
    
    write(*,*) nvsl,nvsr,cntrl_gap,"nvsl,nvsr,cntrl_gap"
    write(*,*) La1,Lx,"La1,Lx"
    write(*,*) bndryCell_xdis,nvsl,nvsr,La1,"~~"
    !stop
    
    deflc_frm_pulley         = 0.20d0 !its used as percent like (0.2*Lb1)
    deflc_frm_cntrl_top_blck = 1.00d0 
    deflc_for_endRgn         = 0.20d0 !its used as absolute
    
  end subroutine wall_cell_lngth_and_other_params
  
  
  
  
  subroutine get_Pulley_and_Neighbours(Pulley,LftNeigh,RghtNeigh)
    implicit none
    integer :: Pulley
    integer :: LftNeigh
    integer :: RghtNeigh
    
    Pulley    = rght_endNode(2) + 1
    LftNeigh  = lft_endNode(1)
    RghtNeigh = rght_endNode(1)
    
  end subroutine get_Pulley_and_Neighbours
  
  
  subroutine Origin_Neighbour_switch(LftNeigh,RghtNeigh,lgcl_NeighSwitch)
    implicit none
    integer, intent(in)    :: LftNeigh,RghtNeigh
    logical, intent(inout) :: lgcl_NeighSwitch
    
    real*8 :: xdis_lft,xdis_rght
    
    lgcl_NeighSwitch = .False.
    
    xdis_lft  = origin(1) - node_xy(LftNeigh,1)
    xdis_rght = origin(1) - node_xy(RghtNeigh,1)
    
    if (xdis_lft.lt.0.0d0 .AND. xdis_rght.gt.0.0d0) then
       lgcl_NeighSwitch = .True.
    elseif (xdis_lft.lt.0.0d0 .AND. xdis_rght.lt.0.0d0) then
       write(*,*) "Both Neighbours in the right"
    elseif (xdis_lft.gt.0.0d0 .AND. xdis_rght.gt.0.0d0) then
       write(*,*) "Both Neighbours in the left"
    else
       continue
    endif
    
  end subroutine Origin_Neighbour_switch
  
  
  subroutine Pulley_Neighbour_switch(Pulley,LftNeigh,RghtNeigh,lgcl_NeighSwitch)
    implicit none
    integer, intent(in)    :: Pulley,LftNeigh,RghtNeigh
    logical, intent(inout) :: lgcl_NeighSwitch
    
    real*8 :: xdis_lft,xdis_rght
    
    lgcl_NeighSwitch = .False.
    
    xdis_lft  = node_xy(Pulley,1) - node_xy(LftNeigh,1)
    xdis_rght = node_xy(Pulley,1) - node_xy(RghtNeigh,1)
    
    if (xdis_lft.lt.0.0d0 .AND. xdis_rght.gt.0.0d0) then
       lgcl_NeighSwitch = .True.
    elseif (xdis_lft.lt.0.0d0 .AND. xdis_rght.lt.0.0d0) then
       write(*,*) "Both Neighbours in the right"
    elseif (xdis_lft.gt.0.0d0 .AND. xdis_rght.gt.0.0d0) then
       write(*,*) "Both Neighbours in the left"
    else
       continue
    endif
    
  end subroutine Pulley_Neighbour_switch
  
  
  
  subroutine get_LftRghtNeigh_ofOrigin(LftNeigh,RghtNeigh)
    implicit none
    integer :: LftNeigh
    integer :: RghtNeigh
    
    LftNeigh  = lft_endNode(1)
    RghtNeigh = rght_endNode(1)

  end subroutine get_LftRghtNeigh_ofOrigin
  
end module system_parameters


module strct_variables
  use sys_building_info
  use non_changing_parameters
  use cell_info
  use system_parameters

  implicit none
  integer :: NsprMax,NcellMax,NnodeMax
  
contains
  
  subroutine allocate_and_initialize_struct_variables
    
    implicit none
    integer :: m,i,j
    integer :: jmax
    integer :: N_elemnt
      
    write(*,*) N_spr,N_cell,"N_spr,N_cell in struct"
    
    NsprMax  = 250 
    NcellMax = 250
    NnodeMax = 250
    
    allocate(l0_strctTN(1:N_struct,1:NsprMax) ,ks_strctTN(1:N_struct,1:NsprMax))
    allocate(A0_strctTN(1:N_struct,1:NcellMax),ka_strctTN(1:N_struct,1:NcellMax))
    allocate(CgX_strctTN(1:N_struct,1:NnodeMax),CgY_strctTN(1:N_struct,1:NnodeMax))
    
    allocate(l0_strctNI(1:N_struct,1:NsprMax) ,ks_strctNI(1:N_struct,1:NsprMax))
    allocate(A0_strctNI(1:N_struct,1:NcellMax),ka_strctNI(1:N_struct,1:NcellMax))
    allocate(CgX_strctNI(1:N_struct,1:NnodeMax),CgY_strctNI(1:N_struct,1:NnodeMax))
    
    allocate(nodeXY_strctTN(1:N_struct,1:NnodeMax,1:N_dmnsn))
    allocate(nodeXY_strctNI(1:N_struct,1:NnodeMax,1:N_dmnsn))
    
    l0_strctTN  = -1.d20 ; ks_strctTN  = -1.d20
    A0_strctTN  = -1.d20 ; ka_strctTN  = -1.d20
    CgX_strctTN = -1.d20 ; CgY_strctTN = -1.d20
    
    l0_strctNI  = -1.d20 ; ks_strctNI  = -1.d20
    A0_strctNI  = -1.d20 ; ka_strctNI  = -1.d20
    CgX_strctNI = -1.d20 ; CgY_strctNI = -1.d20 
    
    N_mvCoordnte1 = 0 ; N_mvCoordnte2 = 0
    NMCWL0A01 = 0     ; NMCWL0A02 = 0 
    nodeXY_strctTN = -1.d20
    
    N_elemnt = 4
    
    do m = 1,N_struct
       
       do i = 1,N_elemnt
          if (i==1) jmax = NsprMax
          if (i==2) jmax = NcellMax
          if (i==3) jmax = NnodeMax
          if (i==4) jmax = NnodeMax
          
          do j = 1,jmax
             !if (i==1) write(*,*) l0_strctTN(m,j),ks_strctTN(m,j),"l0_ks_strctTN"
             !if (i==2) write(*,*) A0_strctTN(m,j),ka_strctTN(m,j),"A0_ka_strctTN"
             !if (i==3) write(*,*) CgX_strctTN(m,j),CgY_strctTN(m,j),"CgX and CgY_strctTN"
             !if (i==4) write(*,*) nodeXY_strctTN(m,j,1:2),"nodeXY_strctTN"
             continue
          enddo
          
       enddo
    enddo
    
  end subroutine allocate_and_initialize_struct_variables
  
end module strct_variables

module node_variables
  use sys_building_info
  use non_changing_parameters
  use cell_info
  use end_info
  use system_parameters
  
  implicit none
  
contains
  
  subroutine allocate_and_initialize_node_variables
    
    implicit none
    
    !write(*,*) N_node,"inside alloc and init node vars"
    
    allocate(node_xy(1:N_node,1:N_dmnsn))
    allocate(node_xyStr(1:N_node,1:N_dmnsn))
    allocate(node_xyStrAN(1:N_node,1:N_dmnsn))    ! AN=Added Node
    allocate(node_xyStrES(1:(N_node+2),1:N_dmnsn))! ES=Equivalent_Spr systm +2 as comes aft AN systm
    allocate(node_xyStrNVFR(1:N_node,1:N_dmnsn))
    
    allocate(node_typ(1:N_node))
    allocate(node_typStr(1:N_node))
    allocate(node_typStrAN(1:N_node))
    allocate(node_typStrES(1:N_node+2))
    allocate(node_typStrNVFR(1:N_node))
    
    allocate(node_cnnctd(1:N_node))
    allocate(node_cnnctdStr(1:N_node))
    allocate(node_cnnctdStrAN(1:N_node))
    allocate(node_cnnctdStrES(1:(N_node+2)))
    allocate(node_cnnctdStrNVFR(1:N_node))
    
    allocate(count_this_dn(1:N_node))
    allocate(count_this_dnStr(1:N_node))
    allocate(count_this_dnStrAN(1:N_node))
    allocate(count_this_dnStrES(1:(N_node+2)))
    allocate(count_this_dnStrNVFR(1:N_node))
    
    allocate(double_node(1:N_node,1:2))
    allocate(double_nodeStr(1:N_node,1:2))
    allocate(double_nodeStrAN(1:N_node,1:2))
    allocate(double_nodeStrES(1:(N_node+2),1:2))
    allocate(double_nodeStrNVFR(1:N_node,1:2))
    
    node_xy  = -1.0d20 ; node_xyStr  = -1.0d20 ; node_xyStrAN  = -1.0d20 ; node_xyStrES  = -1.0d20
    node_typ      = -1 ; node_typStr      = -1 ; node_typStrAN      = -1 ; node_typStrES      = -1
    count_this_dn = -1 ; count_this_dnStr = -1 ; count_this_dnStrAN = -1 ; count_this_dnStrES = -1
    node_cnnctd   =  0 ; node_cnnctdStr   =  0 ; node_cnnctdStrAN   =  0 ; node_cnnctdStrES   =  0
    double_node   = -1 ; double_nodeStr   = -1 ; double_nodeStrAN   = -1 ; double_nodeStrES   = -1
    
    node_xyStrNVFR  = -1.0d20
    node_typStrNVFR      = -1
    count_this_dnStrNVFR = -1
    node_cnnctdStrNVFR   =  0
    double_nodeStrNVFR   = -1
    
  end subroutine allocate_and_initialize_node_variables
  
  subroutine allocate_and_initialize_node_variables_wo_StrVars
    implicit none
    
    allocate(node_xy(1:N_node,1:N_dmnsn))
    allocate(node_typ(1:N_node))
    allocate(node_cnnctd(1:N_node))
    allocate(count_this_dn(1:N_node))
    allocate(double_node(1:N_node,1:N_dmnsn))
    
    node_xy       = -1.0d20
    node_typ      = -1
    count_this_dn = -1
    node_cnnctd   =  0
    double_node   = -1
    
  end subroutine allocate_and_initialize_node_variables_wo_StrVars
  
  subroutine deallocate_node_variables
    implicit none
    
    deallocate(node_xy,node_xyStr,node_xyStrAN)
    deallocate(node_typ,node_typStr,node_typStrAN)
    deallocate(node_cnnctd,node_cnnctdStr,node_cnnctdStrAN)
    deallocate(count_this_dn,count_this_dnStr,count_this_dnStrAN)
    deallocate(double_node,double_nodeStr,double_nodeStrAN)
    
  end subroutine deallocate_node_variables

  
  subroutine allocate_and_initialize_end_node_variables
    implicit none
    
    max_endNode = 3
    
    allocate(lft_endNode(0:max_endNode),rght_endNode(0:max_endNode))
    allocate(clft_top_endNode(0:max_endNode),crght_top_endnode(0:max_endNode))
    allocate(additnl_node_endNode(0:max_endNode))
    allocate(clft_bot_endNode(0:max_endNode),crght_bot_endnode(0:max_endNode))
    
    allocate(lft_nodeXY(1:max_endNode,1:N_dmnsn),rght_nodeXY(1:max_endNode,1:N_dmnsn))

    lft_endNode = 0; rght_endNode = 0;
    clft_top_endNode = 0; crght_top_endNode = 0
    additnl_node_endNode = 0
    clft_bot_endNode = 0; crght_bot_endNode = 0
    
    lft_nodeXY = -1.0d20; rght_nodeXY = -1.0d20
    
  end subroutine allocate_and_initialize_end_node_variables 

  
  subroutine deallocate_end_node_variables
    implicit none
    
    deallocate(lft_endNode,rght_endNode)
    deallocate(clft_top_endNode,crght_top_endnode)
    deallocate(additnl_node_endNode)
    deallocate(clft_bot_endNode,crght_bot_endnode)
    deallocate(lft_nodeXY,rght_nodeXY)
    
  end subroutine deallocate_end_node_variables
  
  
  
  subroutine lft_Nodes_Types_and_Ends
    implicit none
    integer :: i
    
    if (stageNo==1) then
       
       do i = 1,nvsl
          
          if (i==1) then 
             node_xy((count_nodes+1):(count_nodes+2),1) = x_wall 
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
             
          elseif ((i.gt.1).and.(i.le.nvsl)) then
             
             node_xy((count_nodes+1):(count_nodes+2),1) = x_wall + (i-1)*La1     
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
             
          endif
          
          if (i==1) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 2 !1
             
          elseif (i==2) then 
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
             
             lft_nodeXY(1:lft_endNode(0),1:N_dmnsn) = node_xy(lft_endNode(1):lft_endNode(2),1:N_dmnsn)
             
          endif
          
          write(*,*) node_xy(count_nodes+1,1:2),"node_xy_lft"
          write(*,*) node_xy(count_nodes+2,1:2),"node_xy_lft"
          
          count_nodes = count_nodes + 2
          
       enddo

       
       
    elseif (stageNo==2 .or. stageNo==3) then
       write(*,*) "BRING MY ATTENTION in info_mod,lft_Nodes"
       stop
    elseif (stageNo==4) then
       
       do i = 1,nvsl
          !lft_node_strt = count_nodes + 1
          if (i==1) then 
             node_xy((count_nodes+1):(count_nodes+2),1) = x_wall 
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
             
          elseif (i==2) then
             node_xy((count_nodes+1):(count_nodes+2),1) = x_wall + bndryCell_xdis
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
          
          elseif (i.ne.nvsl) then
             
             node_xy((count_nodes+1):(count_nodes+2),1) = x_wall + bndryCell_xdis + (i-2)*La1
          
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
             
          elseif (i.eq.nvsl) then
             node_xy((count_nodes+1),1) = x_wall + bndryCell_xdis + (i-3)*La1 + ((Lb2*unt_len)/2.0d0)
             
             node_xy((count_nodes+2),1) = x_wall + bndryCell_xdis + (i-3)*La1 + ((La1*unt_len_cb)/2.0d0)
             
             node_xy(count_nodes+1,2) = Lb1
             node_xy(count_nodes+2,2) = 0.0d0
             
          endif
          
          if (i==1) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 2
             
          elseif (i==2) then 
             node_typ(count_nodes+1) = 2 
             node_typ(count_nodes+2) = 1 !later it should be 1
             
          else
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 1
             
          endif
          
          if (i .eq. nvsl) then
             lft_endNode(0) = 2
             lft_endNode(1) = count_nodes + 1 
             lft_endNode(2) = count_nodes + 2
             
             lft_nodeXY(1:lft_endNode(0),1:N_dmnsn) = node_xy(lft_endNode(1):lft_endNode(2),1:N_dmnsn)
             
          endif
          
          !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_lft"
          !write(*,*) node_xy(count_nodes+2,1:2),"node_xy_lft"
          
          count_nodes = count_nodes + 2
          
       enddo
       

    endif
    
  end subroutine lft_Nodes_Types_and_Ends
  
  
  
  subroutine tst_lft_Nodes_Types_and_Ends
    implicit none
    integer :: i
    
    do i = 1,nvsl
       
       if (i==1) then 
          node_xy((count_nodes+1):(count_nodes+2),1) = x_wall 
          node_xy(count_nodes+1,2) = Lb1
          node_xy(count_nodes+2,2) = 0.0d0
          
       elseif (i==2) then
          node_xy((count_nodes+1):(count_nodes+2),1) = x_wall + bndryCell_xdis
          node_xy(count_nodes+1,2) = Lb1
          node_xy(count_nodes+2,2) = 0.0d0
          
       else
          
          node_xy((count_nodes+1):(count_nodes+2),1) = x_wall + bndryCell_xdis + (i-2)*La1
          
          node_xy(count_nodes+1,2) = Lb1
          node_xy(count_nodes+2,2) = 0.0d0
          
          
       endif

       if (i==1) then
          node_typ(count_nodes+1) = 2
          node_typ(count_nodes+2) = 2
       else
          node_typ((count_nodes+1):(count_nodes+2)) = 1
       endif

       
       if (i .eq. nvsl) then
          lft_endNode(0) = 2
          lft_endNode(1) = count_nodes + 1 
          lft_endNode(2) = count_nodes + 2
          
          lft_nodeXY(1:lft_endNode(0),1:N_dmnsn) = node_xy(lft_endNode(1):lft_endNode(2),1:N_dmnsn)
          
       endif
       
       count_nodes = count_nodes + 2
       
    enddo   
    
  end subroutine tst_lft_Nodes_Types_and_Ends
  
  
  
  subroutine rght_Nodes_Types_and_Ends
    
    implicit none
    integer :: i
    
    write(*,*) Lx,unt_len,"Lx,unt_len"
    
    
    if (stageNo==1.or.stageNo==2) then
       
       do i = 1,nvsr
          
          if (i==1) then
             node_xy((count_nodes+1):(count_nodes+2),1) = Lx
             node_xy((count_nodes+1),2)                 = Lb1
             node_xy((count_nodes+2),2)                 = 0.0d0
             
          elseif (i.gt.1.and.i.le.nvsr) then
             node_xy((count_nodes+1):(count_nodes+2),1) = Lx - (i-1)*La1
             
             node_xy((count_nodes+1),2) = Lb1
             node_xy((count_nodes+2),2) = 0.0d0
             
          endif
          
          if (i.eq.1) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 2 !1
             
          elseif (i.eq.2) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 1 !later it should be 1
          else
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 1
          endif
          
          if (i .eq. nvsr) then
             rght_endNode(0) = 2
             rght_endNode(1) = count_nodes + 1
             rght_endNode(2) = count_nodes + 2
             
             rght_nodeXY(1:rght_endNode(0),1:N_dmnsn) = node_xy(rght_endNode(1):rght_endNode(2),1:N_dmnsn)
             
          endif
          
          !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_rght"
          !write(*,*) node_xy(count_nodes+2,1:2),"node_xy_rght"
       
          count_nodes = count_nodes + 2
          
       enddo
       
    elseif (stageNo==4) then
       
       do i = 1,nvsr
       
          if (i==1) then
             node_xy((count_nodes+1):(count_nodes+2),1) = Lx
             node_xy((count_nodes+1),2)                 = Lb1
             node_xy((count_nodes+2),2)                 = 0.0d0
             
          elseif (i==2) then
             node_xy((count_nodes+1):(count_nodes+2),1) = Lx - bndryCell_xdis
             
             node_xy((count_nodes+1),2)                 = Lb1
             node_xy((count_nodes+2),2)                 = 0.0d0
             
          elseif (i.ne.nvsr) then
             node_xy((count_nodes+1):(count_nodes+2),1) = Lx - bndryCell_xdis - (i-2)*La1
             
             node_xy((count_nodes+1),2) = Lb1
             node_xy((count_nodes+2),2) = 0.0d0
             
          elseif (i.eq.nvsr) then
             
             node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (i-3)*La1 - (Lb2*unt_len)/2.0d0
             node_xy((count_nodes+2),1) = Lx - bndryCell_xdis - (i-3)*La1 - ((La1*unt_len_cb)/2.0d0)
             
             node_xy((count_nodes+1),2) = Lb1
             node_xy((count_nodes+2),2) = 0.0d0
             
          endif
          
          if (i.eq.1) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 2
             
          elseif (i.eq.2) then
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 1 !later it should be 1
          else
             node_typ(count_nodes+1) = 2
             node_typ(count_nodes+2) = 1
          endif
          
          if (i .eq. nvsr) then
             rght_endNode(0) = 2
             rght_endNode(1) = count_nodes + 1
             rght_endNode(2) = count_nodes + 2
             
             rght_nodeXY(1:rght_endNode(0),1:N_dmnsn) = node_xy(rght_endNode(1):rght_endNode(2),1:N_dmnsn)
             
          endif
          
          !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_rght"
          !write(*,*) node_xy(count_nodes+2,1:2),"node_xy_rght"
       
          count_nodes = count_nodes + 2
          
       enddo
       
    endif
    
  end subroutine rght_Nodes_Types_and_Ends

  
  subroutine cntrlCell_Nodes_Types_and_Ends
    implicit none

    !this subroutine is not needed as no extra addintional node has been added in central
    !cell, still it has been added to maintain similarity
    
    continue
    
    
  end subroutine cntrlCell_Nodes_Types_and_Ends
  
  
  subroutine clft_top_Nodes_Types_and_Ends
    implicit none
    !integer :: depending_blck(1:max_blck)
    integer :: i
    !integer :: start,finish
    
    !start  = count_nodes + 1 !this no is PULLEY
    !finish = count_nodes + cntrl_lft_node
    
    do i = 1,clft_top_node
       
       if (i==1) then
          node_xy((count_nodes+1),1) = bndryCell_xdis + (ncl-2)*La1 + Lb2
          node_xy((count_nodes+1),2) = (Lb1*unt_len)
          
          node_typ(count_nodes+1)    = 0
          
       elseif (i==2) then
          node_xy((count_nodes+1),1) = bndryCell_xdis + (ncl-2)*La1 + Lb2
          node_xy((count_nodes+1),2) = Lb1 * (unt_len - deflc_frm_pulley)
          
          node_typ(count_nodes+1)    = 1
          
       elseif (i==3) then
          node_xy((count_nodes+1),1) = bndryCell_xdis + (ncl-2)*La1 + (unt_len_cb*La1)
          node_xy((count_nodes+1),2) = 0.0d0
          
          node_typ(count_nodes+1)    = 1  
          
       endif
       
       if (i==clft_top_node) then
          clft_top_endNode(0) = 3
          clft_top_endNode(1) = count_nodes - 1 
          clft_top_endNode(2) = count_nodes
          clft_top_endNode(3) = count_nodes + 1 !15
       endif
       
       !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_clft_top"
       
       count_nodes = count_nodes + 1
       
    enddo
    
  end subroutine clft_top_Nodes_Types_and_Ends
  
  
  subroutine additnl_Nodes_Types_and_Ends(cnnctng_node)
    implicit none
    integer, intent(in) :: cnnctng_node
    
    !node_xy((count_nodes+1),1) = (ncl+nccl)*unt_len
    node_xy((count_nodes+1),1) = node_xy(cnnctng_node,1) 
    node_xy((count_nodes+1),2) = node_xy(cnnctng_node,2) - (Lb1*unt_len)

    additnl_node_endNode(0) = 1
    additnl_node_endNode(1) = count_nodes + 1

    node_typ(count_nodes+1)    = 0
    
    !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_additnl_node"

    count_nodes = count_nodes + 1
    
  end subroutine additnl_Nodes_Types_and_Ends
  
  subroutine crght_top_Nodes_Types_and_Ends
    implicit none
    integer :: i
    !integer :: start,finish
    
    !start  = additnl_end_node + 1
    !finish = cntrl_rght_end_node
    
    !cntrl_rght_end_node = finish
    
    do i = 1,crght_top_node
       if (i .eq. 1) then   !Pulley
          node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (ncr-2)*La1 - Lb2
          node_xy((count_nodes+1),2) = (Lb1*unt_len)
          
          node_typ(count_nodes+1)    = 0
          
       elseif (i .eq. 2) then
          node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (ncr-2)*La1 - Lb2
          node_xy((count_nodes+1),2) = Lb1 * (unt_len-deflc_frm_pulley)
          
          node_typ(count_nodes+1)    = 1
          
       else
          node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (ncr-2)*La1 - (unt_len_cb*La1)
          node_xy((count_nodes+1),2) = 0.0d0
          
          node_typ(count_nodes+1)    = 1
          
       endif
       
       if (i.eq.crght_top_node) then
          crght_top_endNode(0) = 3
          crght_top_endNode(1) = count_nodes - 1 
          crght_top_endNode(2) = count_nodes
          crght_top_endNode(3) = count_nodes + 1
       endif
       !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_crght_top"
       
       count_nodes = count_nodes + 1
       
    enddo
    write(*,*) count_nodes,"count_nodes CRT"
    
  end subroutine crght_top_Nodes_Types_and_Ends
  
  
  subroutine clft_bot_Nodes_Types_and_Ends
    implicit none
    !integer :: depending_blck(1:max_blck)
    integer :: i,j
    !integer :: start,finish
    real*8  :: prv_cell_node_xy(1:2,1:2)
        
    !start  = count_nodes + 1 !this no is PULLEY
    !finish = count_nodes + cntrl_lft_node
    
    prv_cell_node_xy = -10000.00d0
    
    do i = 1,ncclb   
       
       do j = 1,2 !j=1 for middle rgn, j=2 for side rgn
          
          if(j == 1) then
             
             if (i.eq.1) then
                node_xy((count_nodes+1),1) = bndryCell_xdis + (ncl-2)*La1 + Lb2
                node_xy((count_nodes+1),2) = node_xy(clft_top_endNode(2),2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
                
                node_typ(count_nodes+1) = 1
                
                prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                
             elseif (i.eq.2) then
                
                if (stageNo==4 .and. stageType==1) then

                   if (i==ncclb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - deflc_for_endRgn
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                      
                   elseif (i.ne.ncclb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                   endif
                   
                elseif (stageNo==4 .and. stageType==2) then
                   
                   node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                   node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
                
                   node_typ(count_nodes+1) = 1
                   prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                   
                endif
                
             elseif (i.ge.3) then
                
                if (stageNo==4 .and. stageType==1) then

                   if (i==ncclb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - deflc_for_endRgn
                      !write(*,*) node_xy((count_nodes+1),1:2),deflc_for_endRgn,"dflc_for_endRgn"
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                      
                   elseif (i.ne.ncclb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - La2
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2) 
                   endif
                   
                elseif (stageNo==4 .and. stageType==2) then 
                   node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                   node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - La2
                   
                   node_typ(count_nodes+1) = 1
                   prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                   
                endif
                
             endif
             
             count_nodes = count_nodes + 1
             
          elseif (j == 2) then
             
             if (i.eq.1) then
                node_xy((count_nodes+1),1) = bndryCell_xdis + (ncl-2)*La1 + (unt_len_cb*La1)
                node_xy((count_nodes+1),2) = - (deflc_frm_pulley*Lb1)/2.0d0 
                
                write(*,*) (count_nodes+1),node_xy((count_nodes+1),2),"insd clft_bot_Nodes_Types_and_Ends"
                
                node_typ(count_nodes+1)    = 1
                
                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
                
             elseif (i.eq.2) then
                node_xy((count_nodes+1),1) = prv_cell_node_xy(2,1)
                node_xy((count_nodes+1),2) = prv_cell_node_xy(2,2) - ((deflc_frm_pulley*Lb1)/2.0d0)
                
                node_typ(count_nodes+1) = 1
                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
                
             else
                node_xy((count_nodes+1),1) = prv_cell_node_xy(2,1)
                node_xy((count_nodes+1),2) = prv_cell_node_xy(2,2) - La2
                
                node_typ(count_nodes+1) = 1
                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
                
             endif
             
             count_nodes = count_nodes + 1
             
          endif
          
          !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_clft_bot"
          !write(*,*) count_nodes,"l1"
       enddo
       
       if (i .ne. ncclb) then
          count_nodes = count_nodes + 2
       endif
       
       !write(*,*) count_nodes,"l2"
       
       if (i .eq. ncclb) then
          clft_bot_endNode(0) = 2
          clft_bot_endNode(1) = count_nodes - 1
          clft_bot_endNode(2) = count_nodes
       endif
       
    enddo
    
    write(*,*) count_nodes,"aft clft_bot_Nodes_Types_and_Ends"
    
  end subroutine clft_bot_Nodes_Types_and_Ends
  
  
  subroutine crght_bot_Nodes_Types_and_Ends
    implicit none
    !integer :: depending_blck(1:max_blck)
    integer :: i,j
    !integer :: start,finish
    real*8  :: prv_cell_node_xy(1:2,1:2)
    !start  = count_nodes + 1 !this no is PULLEY
    !finish = count_nodes + cntrl_lft_node

    prv_cell_node_xy = -10000.00d0
    
    count_nodes = count_nodes - nnecclb*(ncclb-1)
    
    do i = 1,nccrb
       
       do j = 1,2
          
          if(j == 1) then
             
             if (i .eq. 1) then
                node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (ncr-2)*La1 - Lb2
                node_xy((count_nodes+1),2) = node_xy(crght_top_endNode(2),2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
             
                node_typ(count_nodes+1) = 1                
                prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)

             elseif (i.eq.2) then

                if (stageNo==4 .and. stageType==1) then
                   
                   if (i==nccrb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - deflc_for_endRgn
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                      
                   elseif (i.ne.ncclb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
                      
                      node_typ(count_nodes+1) = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                   endif
                   
                elseif (stageNo==4 .and. stageType==2) then
                   node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                   node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - ((deflc_frm_cntrl_top_blck*Lb1)/2.0d0)
                
                   node_typ(count_nodes+1) = 1
                   prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)   
                endif
                
                
             elseif (i.ge.3) then
                
                if (stageNo==4 .and. stageType==1) then

                   if (i==nccrb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - deflc_for_endRgn
                      
                      node_typ(count_nodes+1)    = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)

                   elseif (i.ne.nccrb) then
                      node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                      node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - La2
                      
                      node_typ(count_nodes+1)    = 1
                      prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                   endif
                   
                elseif (stageNo==4 .and. stageType==2) then 
                   node_xy((count_nodes+1),1) = prv_cell_node_xy(1,1)
                   node_xy((count_nodes+1),2) = prv_cell_node_xy(1,2) - La2
                
                   node_typ(count_nodes+1)    = 1
                   prv_cell_node_xy(1,1:2) = node_xy((count_nodes+1),1:2)
                endif
                
             endif
             
             count_nodes = count_nodes + 1
             write(*,*) count_nodes,"count_nodes in CRB nodes_types j1"
             
          elseif (j .eq. 2) then
             
             if (i .eq. 1) then
                node_xy((count_nodes+1),1) = Lx - bndryCell_xdis - (ncr-2)*La1 - (unt_len_cb*La1)
                node_xy((count_nodes+1),2) = - (deflc_frm_pulley*Lb1)/2.0d0
                write(*,*) (count_nodes+1),node_xy((count_nodes+1),2),"insd clft_bot_Nodes_Types_and_Ends"
                node_typ(count_nodes+1)    = 1
                
                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
                
             elseif (i .eq. 2) then
                
                node_xy((count_nodes+1),1) = prv_cell_node_xy(2,1)
                node_xy((count_nodes+1),2) = prv_cell_node_xy(2,2) - ((deflc_frm_pulley*Lb1)/2.0d0)

                node_typ(count_nodes+1)    = 1

                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
                
                
             else
                node_xy((count_nodes+1),1) = prv_cell_node_xy(2,1)
                node_xy((count_nodes+1),2) = prv_cell_node_xy(2,2) - La2

                node_typ(count_nodes+1)    = 1

                prv_cell_node_xy(2,1:2) = node_xy((count_nodes+1),1:2)
             endif

             count_nodes = count_nodes + 1
             write(*,*) count_nodes,"count_nodes in CRB nodes_types j2"
             
          endif

          !write(*,*) j,clft_bot_node,"j,clft_bot_node"
          !if (j.eq.clft_bot_node) then
           !  crght_bot_endNode(0) = 2
            ! crght_bot_endNode(1) = count_nodes
             !crght_bot_endNode(2) = count_nodes + 1
          !endif
          
          !write(*,*) node_xy(count_nodes+1,1:2),"node_xy_clft_bot"
          
       enddo
       
       if (i .ne. nccrb) then
          count_nodes = count_nodes + 2
       elseif (i .eq. nccrb) then
          continue
       endif
       
       if (i .eq. nccrb) then
          crght_bot_endNode(0) = 2
          crght_bot_endNode(1) = count_nodes - 1
          crght_bot_endNode(2) = count_nodes
       endif
       
    enddo
    
  end subroutine crght_bot_Nodes_Types_and_Ends
  
  subroutine endRgnCell_Nodes_Types_and_Ends
    implicit none
    
    !this subroutine is not needed as no extra addintional node has been added in central
    !cell, still it has been added to maintain similarity
    
    continue
    
  end subroutine endRgnCell_Nodes_Types_and_Ends

  subroutine shifting_origin_at_Pulley_Node
    implicit none
    
    integer :: Pulley
    real*8  :: Pulley_xy(1:2)

    integer :: i,j
    
    Pulley = rght_endNode(2) + 1
    !write(*,*) rght_endNode(2),"rghtEN"
    Pulley_xy(1) = node_xy(Pulley,1)
    Pulley_xy(2) = node_xy(Pulley,2)
    
    do i = 1,N_node
       
       do j = 1,N_dmnsn
          node_xy(i,j) = node_xy(i,j) - Pulley_xy(j)
       enddo
       
    enddo
    
  end subroutine shifting_origin_at_Pulley_Node 
  
  subroutine shifting_origin_at_InitiatorCell
    implicit none
    integer :: vrtcl1,vrtcl2
    integer :: node1,node2
    integer :: i,j

    !write(*,*) InitiatorCell,"Init Cell"
    
    if (mod(N_cell,2)==0) then
       origin(1) = node_xy(lft_EndNode(1),1)
       origin(2) = node_xy(lft_EndNode(1),1)
       
    elseif (mod(N_cell,2)==1) then
       vrtcl1 = nvsl
       vrtcl2 = nvsl+nvsr
       
       node1  = 2*vrtcl1-1
       node2  = 2*vrtcl2-1

       !write(*,*) node1,node2,"n1,n2"
       !write(*,*) node_xy(node1,1:2),"nxy1"
       !write(*,*) node_xy(node2,1:2),"nxy2"
       
       origin(1) = (node_xy(node1,1) + node_xy(node2,1))/2.0d0
       origin(2) = (node_xy(node1,2) + node_xy(node2,2))/2.0d0   
       
    endif
    
    do i = 1,N_node
       !write(*,*) node_xy(i,1:2),"before"
       
       do j = 1,N_dmnsn
          node_xy(i,j) = node_xy(i,j) - origin(j)
       enddo
       !write(*,*) node_xy(i,1:2),"aft"
    enddo

    origin(1) = origin(1)-origin(1)
    origin(2) = origin(2)-origin(2)
    
    write(*,*) origin(1:2),"origin"
    
  end subroutine shifting_origin_at_InitiatorCell
  
  
  subroutine print_NodeVars
    implicit none
    integer :: i
    
    open(unit=175,file='NodeVars.dat')!,position='append')
    
    write(175,*) "   ","Node X","   ","Node Y",""
    
    do i = 1,N_node
       write(175,*) node_xy(i,1:2),i
    enddo
    
    write(175,*) " "
    write(175,*) "   ","node_type","   ","node_cnnctd","   ","count_this_dn",""
    
    do i = 1,N_node
       write(175,*) node_typ(i),node_cnnctd(i),count_this_dn(i),i
    enddo

    write(175,*) " "
    write(175,*) "double_node",""
    
    do i = 1,N_node
       write(175,*) double_node(i,1:2)
    enddo
    
    write(175,*) " "
    write(175,*) " "
    
    close(175)
    
  end subroutine print_NodeVars
  

  subroutine print_NodeVars1
    implicit none
    integer :: i
    
    open(unit=175,file='NodeVars1.dat')!,position='append')
    
    write(175,*) "   ","Node X","   ","Node Y",""
    
    do i = 1,N_node
       write(175,*) node_xy(i,1:2),i
    enddo
    
    write(175,*) " "
    write(175,*) "   ","node_type","   ","node_cnnctd","   ","count_this_dn",""
    
    do i = 1,N_node
       write(175,*) node_typ(i),node_cnnctd(i),count_this_dn(i),i
    enddo

    write(175,*) " "
    write(175,*) "double_node",""
    
    do i = 1,N_node
       write(175,*) double_node(i,1:2)
    enddo
    
    write(175,*) " "
    write(175,*) " "
    
    close(175)
    
  end subroutine print_NodeVars1
  
  
  subroutine print_NodeVars2
    implicit none
    integer :: i
    
    open(unit=175,file='NodeVars2.dat')!,position='append')
    
    write(175,*) "   ","Node X","   ","Node Y",""
    
    do i = 1,N_node
       write(175,*) node_xy(i,1:2),i
    enddo
    
    write(175,*) " "
    write(175,*) "   ","node_type","   ","node_cnnctd","   ","count_this_dn",""
    
    do i = 1,N_node
       write(175,*) node_typ(i),node_cnnctd(i),count_this_dn(i),i
    enddo
    
    write(175,*) " "
    write(175,*) "double_node",""
    
    do i = 1,N_node
       write(175,*) double_node(i,1:2)
    enddo
    
    write(175,*) " "
    write(175,*) " "
    
    close(175)
    
  end subroutine print_NodeVars2
  
  
  subroutine nodes_cnnctd_and_count_this_dn
    implicit none
    integer :: i,j

    node_cnnctd    = 0
    count_this_dn  = -1
    
    do i = 1,N_doubleNode
       do j = 1,2
          if(j.eq.1) then
             node_cnnctd(double_node(i,j)) = double_node(i,(j+1))
             count_this_dn(double_node(i,j)) = 1
          else
             node_cnnctd(double_node(i,j)) = double_node(i,(j-1))
          endif
       enddo
    enddo
    
    do j = 1,N_node
       !write(*,*) "Connected nodes and count_this_dn for",j,"=",node_cnnctd(j),count_this_dn(j)
    end do
    
  end subroutine nodes_cnnctd_and_count_this_dn
  
  
  subroutine get_list_of_double_nodes_mthod1
    implicit none
    integer :: i 
    
    double_node = 0
    
    do i = 1,N_doubleNode
       if (i .eq. 1) then
          double_node(i,1) = blck_node + 1
          double_node(i,2) = additnl_end_node + 1
       elseif (i .eq. 2) then
          double_node(i,1) = blck_node + 2
          double_node(i,2) = additnl_end_node + 2
       endif
       !write(*,*) "double_node_number",i,"=",double_node(i,1:2)
    enddo
    
  end subroutine get_list_of_double_nodes_mthod1
  
  subroutine get_list_of_double_nodes_method2
    implicit none
    integer :: i,j,k
    integer :: twice_N_doubleNode
    integer :: count_dn,count_sq
    real*8  :: TINY = 1e-15
    
    integer, allocatable :: double_node_tmp(:,:)
    logical :: lgcl_sq
    
    count_dn = 1
    write(*,*) N_node,"N_node"
    
    do i = 1,N_node
       
       do j = 1,N_node
          if (j.eq.i) then
             continue
          else            
             findn_DoubleNode(1) = .False.
             findn_DoubleNode(2) = .False.
             
             do k = 1,N_dmnsn  
                if(abs(node_xy(i,k)-node_xy(j,k)) .lt. TINY) then
                   findn_DoubleNode(k) = .True.
                   !write(*,*) findn_DoubleNode(k),i,j
                else
                   findn_DoubleNode(k) = .False.
                   !write(*,*) findn_DoubleNode(k),i,j
                endif
             enddo
             
             if ((findn_DoubleNode(1).eqv..True.) .AND. (findn_DoubleNode(2).eqv..True.)) then
                !write(*,*) count_dn,"count_dn"
                double_node(count_dn,1) = i
                double_node(count_dn,2) = j
                count_dn = count_dn + 1
                write(*,*) double_node((count_dn-1),1:2),"double_node"
                !write(*,*) findn_DoubleNode(1),findn_DoubleNode(2),"findn_DoubleNode12"
                write(*,*) count_dn,i,j,"cdn,i,j"
                write(*,*) node_xy(i,1:2),node_xy(j,1:2),i,j,"nodes i,j"
             endif
                         
          endif
       enddo
    enddo
    
    twice_N_doubleNode = count_dn - 1
    N_doubleNode = twice_N_doubleNode/2
    !write(*,*) double_node(1:6,1:2),"double_Nodes_all"
    
    count_sq = 1
    allocate(double_node_tmp(1:N_doubleNode,1:2))
    double_node_tmp = -1
    lgcl_sq = .False.
    
    do i = 1,twice_N_doubleNode
       
       if (count_sq .eq. 1) then
          double_node_tmp(count_sq,1:2) = double_node(1,1:2)
          !write(*,*) double_node_tmp(count_sq,1:2),count_sq,"dn0"
          count_sq = count_sq + 1
       else
          lgcl_sq = .False.
          
          do j = 1,count_sq
             if (double_node(i,1) .eq. double_node_tmp(j,2)) then
                !write(*,*) double_node(i,1),double_node_tmp(j,2),"dn1"
                lgcl_sq = .True.
                exit
             endif
          enddo
          
          if (lgcl_sq .eqv. .False.) then 
             double_node_tmp(count_sq,1:2) = double_node(i,1:2)
             count_sq = count_sq + 1
             !write(*,*) double_node(i,1),double_node_tmp((count_sq-1),1:2),"dn3"
             !write(*,*) count_sq,"count_sq"
          endif
          
       endif
       
    enddo
    
    write(*,*) double_node_tmp,"dn_tmp"
    
    double_node(1:N_doubleNode,1:2) = double_node_tmp
    double_node((N_doubleNode+1):(twice_N_doubleNode),1:2) = -1
    
    open(unit=143,file='doubleNodeList.dat',position='append')
    
    do i = 1,N_doubleNode
       write(143,fmt=*) double_node(i,1:2),i,"dn,dn_order"
    enddo
    
    close(143)
    
  end subroutine get_list_of_double_nodes_method2
  
  subroutine first_or_second_node_of_double_node(node_num,double_node_num,node_positn)
    implicit none
    integer, intent(in)  :: node_num
    integer, intent(out) :: double_node_num,node_positn
    integer :: i,j
    logical :: dn_exists

    dn_exists = .False.
    
    do i = 1,N_doubleNode
       do j = 1,2
          if (node_num .eq. double_node(i,j)) then
             double_node_num = i
             node_positn     = j
             dn_exists = .True.
          endif
       enddo
    enddo

    if (dn_exists .eqv. .True.) then
       continue
    elseif (dn_exists .eqv. .False.) then
       write(*,*) "WARNING: Input Node is not a Double Node",node_num
       write(*,*) "SBRTN: first_or_second_node_if_double_node; info_modules.f08"
       call sleep(1)
       !stop
    endif
    
  end subroutine first_or_second_node_of_double_node
  
  
  subroutine are_they_double_nodes(lgcl_dn,node1,node2)
    !wont be used, this subroutine performs like, if I put two nodes, and anyone of them are first node of any double node, then lgcl will be positive, for example node 13, and node 17 or node 13 and node 18 is given, then lgcl will be true, but it wont serve my purpose
    
    implicit none
    integer, intent(in)  :: node1,node2
    integer :: n1,n2
    logical, intent(out) :: lgcl_dn 

    integer :: n1_tmp
    integer :: i

    n1 = node1
    n2 = node2
    lgcl_dn = .False.
    
    if (n1 .lt. n2) then
       continue
    elseif (n1 .gt. n2) then
       n1_tmp = n1
       n1 = n2
       n2 = n1_tmp
    elseif (n1 .eq. n2) then
       write(*,*) "Both nodes cant be same, sbrtn: are_they_double_nodes"
       stop
    endif
    
    do i = 1,N_doublenode
       if (n1 .eq. double_Node(i,1)) then
          write(*,*) n1,double_Node(i,1),"n1,dn"
          lgcl_dn = .True.
       endif  
    enddo
    
  end subroutine are_they_double_nodes

  subroutine are_they_same_dn(lgcl_sameDn,node1,node2)

    implicit none
    logical, intent(out) :: lgcl_sameDn  
    integer, intent(in)  :: node1,node2
    integer :: i,j
    integer :: dn_no,pstn
    logical :: dn_exists

    lgcl_sameDn = .False.
    dn_exists   = .False.
    ! do i = 1,N_node
    !    write(*,*) double_Node(i,1:2),"double node for i =",i
    ! enddo
    
    do i = 1,N_doubleNode
       do j = 1,2
          if (node1 .eq. double_Node(i,j)) then
             dn_no = i
             pstn  = j
             dn_exists = .True.
             exit
          endif
       enddo
       if (dn_exists .eqv. .True.) exit
    enddo

    if (dn_exists .eqv. .True.) then
       if (pstn .eq. 1) then
          if (node2 .eq. double_Node(dn_no,2)) then
             lgcl_sameDn = .True.
          endif
          
       elseif (pstn .eq. 2) then
          if (node2 .eq. double_Node(dn_no,1)) then
             lgcl_sameDn = .True.
          endif
          
       endif
    endif
    
  end subroutine are_they_same_dn
  
  subroutine is_it_a_dn(node_nm,lgcl_dn)
    !This sbrtn will say if any node is a comp of any double node
    
    implicit none
    integer :: node_nm
    logical :: lgcl_dn

    integer :: i,j

    lgcl_dn = .False.
    
    do i = 1,N_doubleNode

       do j = 1,2 !2 is NOT related to N_dmnsn here
          if (node_nm .eq. double_Node(i,j)) then
             lgcl_dn = .True.
          endif
       enddo
       
    enddo
    
    
  end subroutine is_it_a_dn

  
  subroutine get_the_other_dn(given,other)
    implicit none
    integer, intent(in)  :: given
    integer, intent(out) :: other
    
    integer :: double_node_num,node_positn

    call first_or_second_node_of_double_node(given,double_node_num,node_positn)
    if (node_positn.eq.1) then
       other = double_node(double_node_num,2)
    elseif (node_positn.eq.2) then
       other = double_node(double_node_num,1)
    endif
    
  end subroutine get_the_other_dn

 
  
end module node_variables



module spring_variables  
  use sys_building_info
  use cell_info
  use end_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains
  
  subroutine allocate_and_initialize_spring_variables
    implicit none
    
    !write(*,*) N_spr,"inside alloc and init spr vars"
    
    allocate(typ_spr(1:N_spr),     typ_sprStr(1:N_spr),    typ_sprStrAN(1:N_spr))
    allocate(typ_sprStrES(1:N_spr),typ_sprStrNVFR(1:N_spr))
    
    allocate(k_spr(1:N_spr)     ,k_sprStr(1:N_spr)    ,k_sprStrAN(1:N_spr))
    allocate(k_sprStrES(1:N_spr),k_sprStrNVFR(1:N_spr))
    
    allocate(l0(1:N_spr)      ,l0_Str(1:N_spr)   ,l0_StrAN(1:N_spr))
    allocate(l0_StrES(1:N_spr),l0_StrNVFR(1:N_spr))
    
    allocate(l(1:N_spr)      ,l_Str(1:N_spr)     ,l_StrAN(1:N_spr))
    allocate(l_StrES(1:N_spr),l_StrNVFR(1:N_spr))
    
    allocate(Lt(1:N_spr),alpha(1:N_spr),optmSpr(1:N_spr))
    
    typ_spr = -1      ; typ_sprStr  = -1      ; typ_sprStrAN  = -1      ; typ_sprStrES  = -1 
    k_spr   = -1.0d20 ; k_sprStr    = -1.0d20 ; k_sprStrAN    = -1.0d20 ; k_sprStrES    = -1.0d20 
    l0      = -1.0d25 ; l0_Str      = -1.0d25 ; l0_StrAN      = -1.0d25 ; l0_StrES      = -1.0d25
    l       = -1.0d25 ; l_Str       = -1.0d25 ; l_StrAN       = -1.0d25 ; l_StrES       = -1.0d25
    
    typ_sprStrNVFR  = -1 
    k_sprStrNVFR    = -1.0d20 
    l0_StrNVFR      = -1.0d25
    l_StrNVFR       = -1.0d25
    
    Lt      = -1.0d25
    alpha   = -1.0d20
    optmSpr = -100
    
  end subroutine allocate_and_initialize_spring_variables
  
  subroutine allocate_and_initialize_spring_variables_wo_StrVars
    implicit none
    
    allocate(typ_spr(1:N_spr),k_spr(1:N_spr),l0(1:N_spr),l(1:N_spr))
    
    typ_spr  = -1
    k_spr    = -1.0d20
    l0       = -1.0d25
    l        = -1.0d25
    
  end subroutine allocate_and_initialize_spring_variables_wo_StrVars
  
  subroutine deallocate_spring_variables
    implicit none
    
    deallocate(typ_spr,typ_sprStr,typ_sprStrAN)
    deallocate(k_spr,k_sprStr,k_sprStrAN)
    deallocate(l0,l0_Str,l0_StrAN)
    deallocate(l,l_Str,l_StrAN)
    
    deallocate(Lt,alpha,optmSpr)
    
  end subroutine deallocate_spring_variables
  
  
  subroutine allocate_and_initialize_spring_end_variables
    implicit none
    
    max_endSpring = 2
    
    allocate(lft_endSpring(0:max_endSpring),rght_endSpring(0:max_endSpring))
    allocate(cntrlCell_Spring(0:max_endSpring))
    allocate(clft_top_endSpring(0:max_endSpring),crght_top_endSpring(0:max_endSpring))
    allocate(additnl_node_endSpring(0:max_endSpring))
    allocate(clft_bot_endSpring(0:max_endSpring),crght_bot_endSpring(0:max_endSpring))
    allocate(endCell_Spring(0:max_endSpring))
    
    lft_endSpring = -1; rght_endSpring = -1; cntrlCell_Spring = -1
    clft_top_endSpring = -1; crght_top_endSpring = -1
    additnl_node_endSpring = -1
    clft_bot_endSpring = -1; crght_bot_endSpring = -1
    endCell_Spring = -1
    
  end subroutine allocate_and_initialize_spring_end_variables

  subroutine deallocate_spring_end_variables
    implicit none
    
    deallocate(lft_endSpring,rght_endSpring,cntrlCell_Spring)
    deallocate(clft_top_endSpring,crght_top_endSpring)
    deallocate(additnl_node_endSpring)
    deallocate(clft_bot_endSpring,crght_bot_endSpring)
    deallocate(endCell_Spring)
    
  end subroutine deallocate_spring_end_variables
  
  
  subroutine lft_spring_Properties
    implicit none
    integer :: i
    
    !write(*,*) nsl,"nsl_in_lst"
    
    do i = 1,nsl
       count_spr = count_spr+1
       
       if (i .eq. 1) then
          typ_spr(count_spr) = 1
          l0(count_spr)      = bndryCell_xdis
       elseif (i .eq. 2) then
          typ_spr(count_spr) = 2
          l0(count_spr)      = bndryCell_xdis
          
       elseif (i.gt.2) then
          if (mod(i,3) .eq. 0) then
             typ_spr(count_spr) = 5
             
             if (stageNo==1 .or. stageNo==2) then
                l0(count_spr) = Lb1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-lft_spr"
             elseif  (stageNo==4) then
                if(i.le.(nsl-nsecl)) l0(count_spr) = Lb1
                if(i.gt.(nsl-nsecl)) l0(count_spr) = sqrt(((Lb2-La1*unt_len_cb)/2.0d0)**2 + Lb1**2)
             endif
             
          elseif (mod(i,3) .eq. 1) then
             typ_spr(count_spr) = 3

             if (stageNo==1.or.stageNo==2) then
                l0(count_spr) = La1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-lft_spr"
             elseif (stageNo==4) then
                if(i.le.(nsl-nsecl)) l0(count_spr) = La1             
                if(i.gt.(nsl-nsecl)) l0(count_spr) = Lb2/2.0d0
             endif
             
          else
             typ_spr(count_spr) = 4
             
             if (stageNo==1.or.stageNo==2) then
                l0(count_spr) = La1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-lft_spr"
             elseif (stageNo==4) then
                if(i.le.(nsl-nsecl)) l0(count_spr) = La1
                if(i.gt.(nsl-nsecl)) l0(count_spr) = (La1*unt_len_cb)/2.0d0
             endif
             
          endif
          
       endif

       if (i==1 .or. i==2) then
          !k_spr(count_spr) = (1.0d0)/(howlong_bndryCell)
          k_spr(count_spr) = 1.0d0
       else
          k_spr(count_spr) = 1.0d0
       endif
       
       Lt(count_spr)    = l0(count_spr)
       alpha(count_spr) = 1.0d0
       
       if (i.eq.nsl) then
          lft_endSpring(0) = 1
          lft_endSpring(1) = count_spr
       endif
       
       !write(*,*) typ_spr(i),"type of spring i",i

       if (i==nsl) initial_l0_ofCP = l0(count_spr)
    enddo
    
    !write(*,*) lft_endSpring(0:1),"Les"
    open(unit=87,file='initiall0_ofCP.dat')
    write(87,fmt=*) initial_l0_ofCP,"1"
    close(87)
    
  end subroutine lft_spring_Properties


  
  subroutine tst_lft_spring_Properties
    implicit none
    integer :: i
    
    !write(*,*) nsl,"nsl_in_lst"
    
    do i = 1,nsl
       count_spr = count_spr+1
       
       if (i .eq. 1) then
          typ_spr(count_spr) = 1
          l0(count_spr)      = bndryCell_xdis
       elseif (i .eq. 2) then
          typ_spr(count_spr) = 2
          l0(count_spr)      = bndryCell_xdis
          
       elseif (i.gt.2) then
          
          if (mod(i,3) .eq. 0) then
             typ_spr(count_spr) = 5  
             l0(count_spr) = Lb1
                
          elseif (mod(i,3) .eq. 1) then
             typ_spr(count_spr) = 3   
             l0(count_spr) = La1
           
          else
             typ_spr(count_spr) = 4
             l0(count_spr) = La1

          endif

       endif

       if (i==1 .or. i==2) then
          !k_spr(count_spr) = (1.0d0)/(howlong_bndryCell)
          k_spr(count_spr) = 1.0d0
       else
          k_spr(count_spr) = 1.0d0
       endif
       
       Lt(count_spr)    = l0(count_spr)
       alpha(count_spr) = 1.0d0
       
       if (i.eq.nsl) then
          lft_endSpring(0) = 1
          lft_endSpring(1) = count_spr
       endif
       
       !write(*,*) typ_spr(i),"type of spring i",i
       
    enddo

    !write(*,*) lft_endSpring(0:1),"Les"
    
  end subroutine tst_lft_spring_Properties



  
  subroutine rght_spring_Properties
    implicit none
    integer :: i
    !integer :: start,finish

    !start  = nsl + 1
    !finish = nsl + nsr
    !write(*,*) start,finish,"strt_finsh_in_rght_spr_types"
    
    do i = 1,nsr
       count_spr = count_spr+1
       
       
       if (i.eq.1) then
          typ_spr(count_spr) = 1
          l0(count_spr)      = bndryCell_xdis
          
       elseif (i.eq.2) then
          typ_spr(count_spr) = 2
          l0(count_spr)      = bndryCell_xdis
          
       elseif (i.gt.2) then
          if (mod(i,3) .eq. 0) then
             typ_spr(count_spr) = 5
             
             if (stageNo==1.or.stageNo==2) then
                l0(count_spr) = Lb1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-rght_spr"
             elseif (stageNo==4) then
                if(i.le.(nsr-nsecr)) l0(count_spr) = Lb1
                if(i.gt.(nsr-nsecr)) l0(count_spr) = sqrt(((Lb2-La1*unt_len_cb)/2.0d0)**2 + Lb1**2)
             endif
             
          elseif (mod(i,3) .eq. 1) then
             typ_spr(count_spr) = 3
             
             if (stageNo==1.or.stageNo==2) then
                l0(count_spr) = La1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-rght_spr"
             elseif (stageNo==4) then
                if(i.le.(nsr-nsecr)) l0(count_spr) = La1
                if(i.gt.(nsr-nsecr)) l0(count_spr) = Lb2/2.0d0
             endif
             
          else
             typ_spr(count_spr) = 4
             
             if (stageNo==1.or.stageNo==2) then
                l0(count_spr) = La1
             elseif (stageNo==3) then
                write(*,*) "BRING MY ATTENTION,info-rght_spr"
             elseif (stageNo==4) then
                if(i.le.(nsr-nsecr)) l0(count_spr) = La1
                if(i.gt.(nsr-nsecr)) l0(count_spr) = (La1*unt_len_cb)/2.0d0
             endif
             
          endif
          
       endif
       
       if (i==1 .or. i==2) then
          !k_spr(count_spr) = (1.0d0)/(howlong_bndryCell) !needed for long enough first cell
          k_spr(count_spr) = 1.0d0
       else
          k_spr(count_spr) = 1.0d0
       endif
       
       Lt(count_spr)    = l0(count_spr)
       alpha(count_spr) = 1.0d0
       
       
       !write(*,*) typ_spr(i),"type of spring i",i
       
       
       if (i.eq.nsr) then
          rght_endSpring(0) = 1
          rght_endSpring(1) = count_spr
       endif
       if (i==nsr) initial_l0_ofCP = l0(count_spr)
    enddo
    
    open(unit=87,file='initiall0_ofCP.dat')
    write(87,fmt=*) initial_l0_ofCP,"2"
    close(87)
    
    !write(*,*) rght_endSpring(0:1),"reS"
    
  end subroutine rght_spring_Properties

  
  subroutine cntrlCell_spring_Properties
    implicit none
    integer :: i
    
    if (stageNo.ne.1 .or. stageType.ne.1) then
       write(*,*) "stageNo and stageType variation,fl:info_mod,sb:cntrlCell_spring"
    endif
    
    do i = 1,nsc
       count_spr = count_spr+1
       
       if (i==1) then
          typ_spr(count_spr) = 3
          l0(count_spr)      = La1
       elseif (i==2) then
          typ_spr(count_spr) = 4
          l0(count_spr)      = La1
       endif
       
       k_spr(count_spr) = 1.0d0
       Lt(count_spr)    = l0(count_spr)
       alpha(count_spr) = 1.0d0
       
       if (i==1) cntrlCell_Spring(0) = 2
       cntrlCell_Spring(i) = count_spr
       
    enddo
    
    
  end subroutine cntrlCell_spring_Properties

  
  subroutine clft_top_spring_Properties
    implicit none
    integer :: i
    
    !write(*,*) blck_spr,cntrl_lft_end_spr,"blck_spr,cntrl_lft_spr"

    clft_top_endSpring(0) = 2 !in lft_rght, it was done inside
    
    do i = 1, clft_top_spr
       count_spr = count_spr+1

       if (i.eq.1) then
          typ_spr(count_spr) = 6  ! count_spr = 13,Pulley
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = (Lb2/2.0d0) + (deflc_frm_pulley*Lb1)
          ! write(*,*) l0(count_spr)
          ! write(*,*) Lb2,deflc_frm_pulley,Lb1
          ! stop
          alpha(count_spr)   = 1.0d0
          
       elseif (i .eq. 2) then 
          typ_spr(count_spr) = 7 ! count_spr = 14,Diagonal
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = sqrt((Lb1*(unt_len-deflc_frm_pulley))**2 + &
               (Lb2*unt_len-La1*unt_len_cb)**2)
          
          !write(*,*) unt_len,deflc_frm_pulley,unt_len_cb,"unt"
          !write(*,*) La,Lb,"LaLb"
          !write(*,*) l0(count_spr),"l0"
          
          alpha(count_spr)   = 1.0d0
          
       elseif (i .eq. 3) then
          typ_spr(count_spr) = 8  ! count_spr = 15,Basal_cntrl
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = (La1*unt_len_cb)/2.0d0
          alpha(count_spr)   = 1.0d0
          !write(*,*) count_spr,"basal cntrl"
       endif

       Lt(count_spr)      = l0(count_spr)
       
       if (i.eq.1) then
          clft_top_endSpring(1) = count_spr
       elseif (i.eq.clft_top_spr) then
          clft_top_endSpring(2) = count_spr
       endif
       
    enddo

    !write(*,*) clft_top_endSpring(0:2),"CLTes"
    
  end subroutine clft_top_spring_Properties
  
  subroutine additnl_spring_Properties
    implicit none
    
    count_spr = count_spr + 1
    typ_spr(count_spr) = 9 !addrnl_type
    k_spr(count_spr)   = 1.0d0
    l0(count_spr)      = La2
    Lt(count_spr)      = l0(count_spr)
    alpha(count_spr)   = 1.0d0
    
    additnl_node_endSpring(0) = 1
    additnl_node_endSpring(1) = count_spr

    !write(*,*) additnl_node_endSpring(0:1),"ANes"
    
  end subroutine additnl_spring_Properties
  
  subroutine crght_top_spring_Properties
    implicit none
    integer :: i

    crght_top_endSpring(0) = 2
    
    do i = 1, crght_top_spr
       count_spr = count_spr + 1
       !write(*,*) count_spr,"count_spr"
       
       if (i .eq. 1) then
          typ_spr(count_spr) = 6  ! count_spr = 16,Pulley
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = (Lb2/2.0d0) + (deflc_frm_pulley*Lb1)
          alpha(count_spr)   = 1.0d0
          
       elseif (i.eq.2) then
          typ_spr(count_spr) = 7  ! count_spr = 17, Diagonal
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = sqrt((Lb1*(unt_len-deflc_frm_pulley))**2 + &
               (Lb2*unt_len-La1*unt_len_cb)**2)
          !write(*,*) l0(count_spr),"l0"
          
          alpha(count_spr)   = 1.0d0
          
       elseif (i.eq.3) then
          typ_spr(count_spr) = 8  ! count_spr = 18,Basal_cntrl
          k_spr(count_spr)   = 1.0d0
          l0(count_spr)      = (La1*unt_len_cb)/2.0d0
          alpha(count_spr)   = 1.0d0
          !write(*,*) count_spr,"basal cntrl"
       endif

       Lt(count_spr) = l0(count_spr)
       
       if (i .eq. 1) then
          crght_top_endSpring(1) = count_spr
       elseif (i .eq. crght_top_spr) then
          crght_top_endSpring(2) = count_spr
       endif
       
    enddo
    !write(*,*) crght_top_endSpring(0:2),"CRTes"
    
  end subroutine crght_top_spring_Properties


  
  subroutine clft_bot_spring_Properties
    implicit none
    integer :: i
    
    clft_bot_endSpring(0) = 1
    
    do i = 1,nsclb
      
       if (i.ne.1 .and. mod(i,nsecclb).eq.1) then !nsecclb=3
          count_spr = count_spr + 4
       else
          count_spr = count_spr + 1
       endif
       
       write(*,*) count_spr,"count_spr_in_clft_bot"
       
       if (mod(i,3).eq.1) then !Apcl
          typ_spr(count_spr) = 3  ! count_spr = 19,Apical_NonBoundary
          k_spr(count_spr)   = 1.0d0
          
          if ((i/nsecclb).lt.2) then
             l0(count_spr) = (Lb1*deflc_frm_cntrl_top_blck)/2.0d0
             !write(*,*) l0(count_spr),"l0"
             
          else
             if (stageNo==4 .and. stageType==1) then

                if ((i/nsecclb).lt.(ncclb-1)) then
                   l0(count_spr) = La2
                elseif ((i/nsecclb)==(ncclb-1)) then
                   l0(count_spr) = deflc_for_endRgn
                endif
                
             elseif (stageNo==4 .and. stageType==2) then
                l0(count_spr) = La2
             endif
             
          endif
          
          alpha(count_spr)   = 1.0d0
          
       elseif (mod(i,3) .eq. 2) then !Bsal
          typ_spr(count_spr) = 4 ! count_spr = 20,Basal_NonBoundary 
          k_spr(count_spr)   = 1.0d0
          
          if ((i/nsecclb) .lt. 2) then
             l0(count_spr) = (Lb1*deflc_frm_pulley)/2.0d0 !0.10*Lb
             !write(*,*) l0(count_spr),"l0"
             
          else
             l0(count_spr) = La2
          endif
          
          alpha(count_spr)   = 1.0d0
          
       elseif (mod(i,3) .eq. 0) then
          typ_spr(count_spr) = 5  ! count_spr = 21,Verical
          k_spr(count_spr)   = 1.0d0
          
          if ((i/nsecclb) == 1) then
             l0(count_spr) = sqrt((Lb2*unt_len-La1*unt_len_cb)**2 + &
                  ((Lb1*(deflc_frm_cntrl_top_blck-deflc_frm_pulley))/2.0d0)**2)
             !write(*,*) l0(count_spr),"l0"
             
          else
             if (stageNo==4 .and. stageType==1) then
                
                if ((i/nsecclb).lt.ncclb) then
                   l0(count_spr) = Lb2*unt_len - La1*unt_len_cb

                elseif ((i/nsecclb)==ncclb) then
                   l0(count_spr) = sqrt((Lb2*unt_len-La1*unt_len_cb)**2 + &
                        (La2-deflc_for_endRgn)**2)
                  
                endif
                
             elseif (stageNo==4 .and. stageType==2) then
                l0(count_spr) = Lb2*unt_len - La1*unt_len_cb
             endif
          endif
          
          alpha(count_spr)   = 1.0d0
          
       endif
       
       Lt(count_spr) = l0(count_spr)
       
       write(*,*) k_spr(count_spr),l0(count_spr),typ_spr(count_spr),count_spr,"CLB"
       
       if (i == nsclb) then
          clft_bot_endSpring(1) = count_spr
       endif
       
    enddo
    
    !write(*,*) clft_bot_endSpring(0:1),"CLBes"
    
  end subroutine clft_bot_spring_Properties
  
  
  subroutine crght_bot_spring_Properties
    implicit none
    integer :: i
    
    crght_bot_endSpring(0) = 1 
    count_spr = nsl + nsr + nsclt + nscrt + nsecclb
    
    do i = 1,nscrb 
       
       if (i.ne.1 .and. mod(i,nseccrb).eq.1) then !nseccrb=3
          count_spr = count_spr + 4
       else
          count_spr = count_spr + 1
       endif
       
       write(*,*) count_spr,"count_spr_in_crght_bot"
       
       if (mod(i,3).eq.1) then
          typ_spr(count_spr) = 3
          k_spr(count_spr)   = 1.0d0
          
          if ((i/nseccrb).lt.2) then
             l0(count_spr) = (Lb1*deflc_frm_cntrl_top_blck)/2.0d0
             !write(*,*) l0(count_spr),"l0"
             
          else
             if (stageNo==4 .and. stageType==1) then
                
                if ((i/nseccrb).lt.(nccrb-1)) then
                   l0(count_spr) = La2
                elseif ((i/nseccrb)==(nccrb-1)) then
                   l0(count_spr) = deflc_for_endRgn
                endif
                
             elseif (stageNo==4 .and. stageType==2) then
                l0(count_spr) = La2
             endif
             
          endif
          
          alpha(count_spr)   = 1.0d0
          
       elseif (mod(i,3) .eq. 2) then 
          typ_spr(count_spr) = 4
          k_spr(count_spr)   = 1.0d0
          
          if ((i/nseccrb) .lt. 2) then
             l0(count_spr) = (Lb1*deflc_frm_pulley)/2.0d0 !0.10*Lb
          else
             l0(count_spr) = La2
          endif
          
          alpha(count_spr)   = 1.0d0
           
       elseif (mod(i,3) .eq. 0) then
          typ_spr(count_spr) = 5
          k_spr(count_spr)   = 1.0d0
           
          if ((i/nseccrb) .eq. 1) then
             l0(count_spr) = sqrt((Lb2*unt_len-La1*unt_len_cb)**2 + &
                  ((Lb1*(deflc_frm_cntrl_top_blck-deflc_frm_pulley))/2.0d0)**2)
             !write(*,*) l0(count_spr),"l0"
             
          else
             if (stageNo==4 .and. stageType==1) then
                
                if ((i/nseccrb).lt.nccrb) then
                   l0(count_spr) = Lb2*unt_len - La1*unt_len_cb
                   
                elseif ((i/nseccrb)==nccrb) then
                   l0(count_spr) = sqrt((Lb2*unt_len-La1*unt_len_cb)**2 + &
                        (La2-deflc_for_endRgn)**2)
                   
                endif
                
             elseif (stageNo==4 .and. stageType==2) then
                l0(count_spr) = Lb2*unt_len - La1*unt_len_cb
             endif
             
          endif
          
          alpha(count_spr)   = 1.0d0
          
       endif
       
       Lt(count_spr) = l0(count_spr)
       
       if (i == nscrb) then
          crght_bot_endSpring(1) = count_spr
       endif
       
    enddo
    
    !write(*,*) crght_bot_endSpring(0:2),"CRBes"
    
  end subroutine crght_bot_spring_Properties
  
  subroutine endRgnCell_spring_Properties
    implicit none
    integer :: i

    if (stageNo.ne.4 .or. stageType.ne.1) then
       write(*,*) "stageNo and stageType variation,fl:info_mod,sb:cntrlCell_spring"
    endif
    
    do i = 1,nser
       count_spr = count_spr + 1
       
       if (i==1) then
          typ_spr(count_spr) = 4
          l0(count_spr)      = 2.0d0*(Lb2*unt_len - La1*unt_len_cb)
          
          k_spr(count_spr) = 1.0d0
          Lt(count_spr)    = l0(count_spr)
          alpha(count_spr) = 1.0d0
          
          endCell_Spring(1)  = count_spr
       endif
       
    enddo
    
  end subroutine endRgnCell_spring_Properties
  
  
  subroutine get_N_optmSpr_and_optmSpr
    implicit none   
    integer :: i
    
    N_optmSpr = N_spr
    
    do i = 1,N_spr
       optmSpr(i) = 1
       !write(*,*) optmSpr(i),i,"optmSpr,Spr No"
    enddo
    
  end subroutine get_N_optmSpr_and_optmSpr
  
  subroutine joint_spr_Prop(spr1,spr2,ks_jnt,l0s_jnt)
    implicit none
    integer, intent(in)  :: spr1,spr2
    real*8,  intent(out) :: ks_jnt,l0s_jnt
    
    real :: ks1,ks2,l0s1,l0s2
    
    write(*,*) spr1,spr2,"spr1-2"
    
    ks1  = k_spr(spr1)
    ks2  = k_spr(spr2)
    l0s1 = l0(spr1)
    l0s2 = l0(spr2)
    write(*,*) ks1,ks2,l0s1,l0s2,"ks1-2,l0s1-2"
    
    ks_jnt  = ks1+ks2
    l0s_jnt = (ks1*l0s1 + ks2*l0s2)/(ks1+ks2)
    
    write(*,*) ks_jnt,l0s_jnt,"ks_jnt,l0s_jnt"
     
  end subroutine joint_spr_Prop
   
  subroutine print_spring_variables
    implicit none
    integer :: i
    
    open(unit=145,file='Spring_Variables.dat')
    
    write(145,fmt=*) "typ_spr","k_spr","l0","Lt","alpha","optmSpr","spr_nm"
    
    do i = 1,N_spr
       write(145,fmt=*) typ_spr(i),l0(i),k_spr(i),Lt(i),alpha(i),optmSpr(i),i
    enddo
    
    write(145,fmt=*) ""
    close(145)
    
  end subroutine print_spring_variables
  
  subroutine print_sprVars_wo_Lt_alpha_optmSpr
    implicit none
    integer :: i
    
    open(unit=146,file='Spring_Vars.dat')
    
    write(146,fmt=*) "typ_spr","k_spr","l0","l","spr_nm"
    
    do i = 1,N_spr
       write(146,fmt=*) typ_spr(i),k_spr(i),l0(i),l(i),i
    enddo
    
    write(146,fmt=*) ""
    
    close(146)
    
  end subroutine print_sprVars_wo_Lt_alpha_optmSpr
  
end module spring_variables
 
 
module diagonal_variables
  use sys_building_info
  use cell_info
  use end_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains

  subroutine allocate_and_initialize_diagonal_variables
    implicit none

    allocate(D(1:N_diag),Dt(1:N_diag))
    allocate(gammaa(1:N_diag))

    D      = -1e25
    Dt     = -1e25
    gammaa = -1e20
    
  end subroutine allocate_and_initialize_diagonal_variables

  subroutine deallocate_diagonal_variables
    implicit none
    
    deallocate(D,Dt)
    deallocate(gammaa)
     
  end subroutine deallocate_diagonal_variables

  subroutine lft_diagonal_Properties
    implicit none
    integer :: i
    
    do i = 1,ndl
       count_diag = count_diag+1
       Dt(count_diag)     = sqrt(unt_lenX**2 + unt_lenY**2)
       gammaa(count_diag) = 1.0d0
    enddo
    
  end subroutine lft_diagonal_Properties

  subroutine rght_diagonal_Properties
    implicit none
    integer :: i
    
    do i = 1,ndr
       count_diag = count_diag + 1
       Dt(count_diag)     = sqrt(unt_lenX**2 + unt_lenY**2)
       gammaa(count_diag) = 1.0d0
    enddo
    
  end subroutine rght_diagonal_Properties

  subroutine clft_top_diagnonal_Properties
    implicit none
    integer :: i
    
    do i = 1,clft_top_diag
       count_diag = count_diag + 1
       
       if (i.eq.1) then
          Dt(count_diag) = sqrt(unt_lenX**2 + deflc_frm_Pulley**2)
       elseif (i.eq.2) then
          Dt(count_diag) = sqrt(unt_lenY**2 + unt_len_cb**2)
       elseif (i.eq.3) then
          Dt(count_diag) = sqrt(unt_lenX**2 + unt_lenY**2)
       elseif (i.eq.4) then
          Dt(count_diag) = sqrt(unt_lenX**2 + (unt_lenY-deflc_frm_Pulley)**2)
       elseif (i.eq.5) then
          Dt(count_diag) = sqrt((unt_lenX-unt_len_cb)**2 + unt_lenY**2)
       endif
       
    enddo
    
  end subroutine clft_top_diagnonal_Properties


  subroutine crght_top_diagonal_Properties
    implicit none
    integer :: i
    
    do i = 1,crght_top_diag
       count_diag = count_diag + 1
       
       if (i.eq.1) then
          Dt(count_diag) = sqrt(unt_lenX**2 + deflc_frm_Pulley**2)
       elseif (i.eq.2) then
          Dt(count_diag) = sqrt(unt_lenY**2 + unt_len_cb**2)
       elseif (i.eq.3) then
          Dt(count_diag) = sqrt(unt_lenX**2 + unt_lenY**2)
       elseif (i.eq.4) then
          Dt(count_diag) = sqrt(unt_lenX**2 + (unt_lenY-deflc_frm_Pulley)**2)
       elseif (i.eq.5) then
          Dt(count_diag) = sqrt((unt_lenX-unt_len_cb)**2 + unt_lenY**2)
       endif
       
    enddo

    
  end  subroutine crght_top_diagonal_Properties


  subroutine clft_bot_diagonal_Properties
    implicit none
    

    
  end subroutine clft_bot_diagonal_Properties


  subroutine crght_bot_diagonal_Properties


  end subroutine crght_bot_diagonal_Properties

  subroutine print_diagonal_variables


  end subroutine print_diagonal_variables



end module diagonal_variables

module area_variables
  use sys_building_info
  use cell_info
  use end_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains
  
  subroutine allocate_and_initialize_area_variables
    implicit none
    
    allocate(k_area(1:N_cell),A0(1:N_cell),A(1:N_cell))
    allocate(At(1:N_cell),fctr(1:N_cell),beta(1:N_cell),optmCell(1:N_cell))
    allocate(typ_area(1:N_cell),area_sides(1:N_cell))
    allocate(lftSideCell(1:N_lftCells),rghtSideCell(1:N_rghtCells))
    allocate(cntrlCell(1:N_cntrlCells),endRgnCell(1:N_endRgnCells))
    allocate(k_areaStrNVFR(1:N_cell),A0_StrNVFR(1:N_cell),A_StrNVFR(1:N_cell))
    
    typ_area   = -1
    area_sides = -1
    k_area     = -1.0d20 ; k_areaStrNVFR = -1.0d20  
    A0         = -1.0d25 ; A0_StrNVFR    = -1.0d25
    A          = -1.0d25 ; A_StrNVFR     = -1.0d25
    At         = -1.0d25
    fctr       = -1.0d25
    beta       = -1.0d20
    optmCell   = -100
    
    lftSideCell  = -1
    rghtSideCell = -1
    cntrlCell    = -1
    endRgnCell   = -1
    
  end subroutine allocate_and_initialize_area_variables
  
  
  subroutine allocate_and_initialize_area_variables_wo_StrVars
    implicit none
    
    allocate(k_area(1:N_cell),A0(1:N_cell),A(1:N_cell))
    
    k_area(1:N_cell) = -1.0d20
    A0(1:N_cell)     = -1.0d25
    A(1:N_cell)      = -1.0d25
    
  end subroutine allocate_and_initialize_area_variables_wo_StrVars
  
  subroutine deallocate_area_variables
    implicit none
    
    deallocate(k_area,A,A0,At,fctr,beta,optmCell)
    deallocate(typ_area,area_sides)
    deallocate(lftSideCell,rghtSideCell,cntrlCell,endRgnCell)
    
  end subroutine deallocate_area_variables

  
  subroutine allocate_and_initialize_area_end_variables
    implicit none
    
    max_endArea   = 2
    
    allocate(lft_endArea(1:max_endArea),rght_endArea(1:max_endArea))
    allocate(clft_top_endArea(1:max_endArea),crght_top_endArea(1:max_endArea))
    allocate(additnl_node_endArea(1:max_endArea))
    allocate(clft_bot_endArea(1:max_endArea),crght_bot_endArea(1:max_endArea))
    
    lft_endArea = 0; rght_endArea = 0
    clft_top_endArea = 0; crght_top_endArea = 0
    additnl_node_endArea = 0
    clft_bot_endArea = 0; crght_bot_endArea = 0
    
  end subroutine allocate_and_initialize_area_end_variables


  subroutine deallocate_area_end_variables
    implicit none
    
    deallocate(lft_endArea,rght_endArea)
    deallocate(clft_top_endArea,crght_top_endArea)
    deallocate(additnl_node_endArea)
    deallocate(clft_bot_endArea,crght_bot_endArea)
    
    
  end subroutine deallocate_area_end_variables
  
  
  subroutine lft_area_Properties
    implicit none
    integer :: i
    
    do i = 1,ncl
       count_area = count_area + 1
       typ_area(count_area) = 1
       area_sides(count_area) = 4
       k_area(count_area)   = 1.0d0
       
       if(i.eq.1) then
          A0(count_area) = bndryCell_xdis*Lb1
       elseif(i.ne.1 .and. i.ne.ncl) then
          A0(count_area) = La1*Lb1
          
       elseif(i.eq.ncl) then
          
          if (stageNo==1.or.stageNo==2) then
             A0(count_area) = La1*Lb1
          elseif (stageNo==3) then
             write(*,*) "BRING MY ATTENTION,info_modules,sb:lft_area_Properties"
          elseif (stageNo==4) then
             A0(count_area) = 0.5* ((Lb2/2.0d0) + ((La1*unt_len_cb)/2.0d0))* (Lb1)
          endif
          
       endif
       
       At(count_area)       = A0(count_area)
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0
       
    enddo
    
  end subroutine lft_area_Properties


  subroutine tst_lft_area_Properties
    implicit none
    integer :: i
    
    do i = 1,ncl
       count_area = count_area + 1
       typ_area(count_area) = 1
       area_sides(count_area) = 4
       k_area(count_area)   = 1.0d0
       
       if(i.eq.1) then
          A0(count_area) = bndryCell_xdis*Lb1
       elseif(i.ne.1) then
          A0(count_area) = La1*Lb1          
       endif
       
       At(count_area)       = A0(count_area)
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0
       
    enddo
    
  end subroutine tst_lft_area_Properties
  
  
  
  subroutine rght_area_Properties
    implicit none
    integer :: i
    
    do i = 1,ncr
       count_area = count_area + 1
       typ_area(count_area) = 1
       area_sides(count_area) = 4
       k_area(count_area)   = 1.0d0
       
       if (i.eq.1) then
          A0(count_area) = bndryCell_xdis*Lb1
       elseif (i.ne.1 .and. i.ne.ncr) then
          A0(count_area) = La1*Lb1
       elseif (i.eq.ncr) then

          if (stageNo==1) then
             A0(count_area) = La1*Lb1
          elseif (stageNo==2 .or. stageNo==3) then
             write(*,*) "BRING MY ATTENTION, fl:info,sb:rght_area"
          elseif (stageNo==4) then
             A0(count_area) =0.5*((Lb2/2.0d0) + ((La1*unt_len_cb)/2.0d0))* (Lb1)
          endif
          
       endif
       
       At(count_area)       = A0(count_area)
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0
       
    enddo
    
  end subroutine rght_area_Properties
  
  
  subroutine cntrlCell_area_Properties
    implicit none
    integer :: i
    
    do i = 1,ncc
       
       count_area = count_area + 1
       !write(*,*) count_area,"ca"
       typ_area(count_area)   = 1
       area_sides(count_area) = 4
       k_area(count_area)     = 1.0d0
       
       A0(count_area)       = La1*Lb1
       At(count_area)       = A0(count_area)!not needed, simlrty with othr rtns
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0 !not needed, simlrty with othr rtns
       
    enddo
    
  end subroutine cntrlCell_area_Properties
  
  subroutine clft_top_area_Properties
    implicit none
    integer :: i
    
    do i = 1,ncclt
       count_area = count_area+1
       typ_area(count_area) = 2 !change frm lft_rght
       area_sides(count_area) = 5
       k_area(count_area)   = 1.0d0
       
       A0(count_area)       = 0.5d0*(Lb2*unt_len-La1*unt_len_cb)*Lb1*(unt_len-deflc_frm_pulley) &
            + (La1*unt_len_cb)*Lb1 + (Lb2*unt_len-La1*unt_len_cb)*(Lb1*deflc_frm_pulley) &
            - 0.5* ((Lb2/2.0d0) + ((La1*unt_len_cb)/2.0d0))* (Lb1)
       
       !Triangle+ 2Rectangle(0.5*0.9*0.8 + 0.10*1.0 +0.9*0.2) - A0(ncl)
       !write(*,*) A0(count_area),"A0"
       
       At(count_area)       = A0(count_area)
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0
    enddo
    
  end subroutine clft_top_area_Properties

  
  subroutine crght_top_area_Properties
    implicit none
    integer :: i
    
    do i = 1,nccrt
       count_area = count_area+1
       typ_area(count_area) = 2 !change frm lft_rght
       area_sides(count_area) = 5
       k_area(count_area)   = 1.0d0
       A0(count_area)       = 0.5d0*(Lb2*unt_len-La1*unt_len_cb)*Lb1*(unt_len-deflc_frm_pulley) &
            + (La1*unt_len_cb)*Lb1 + (Lb2*unt_len-La1*unt_len_cb)*(Lb1*deflc_frm_pulley) &
            - 0.5* ((Lb2/2.0d0) + ((La1*unt_len_cb)/2.0d0))* (Lb1)

       !Triangle+ 2Rectangle(0.5*0.9*0.8 + 0.10*1.0 + 0.9*0.2) - A0(ncr)
       write(*,*) A0(count_area),"A0"
       
       At(count_area)       = A0(count_area)
       fctr(count_area)     = 1.0d0
       beta(count_area)     = 1.0d0
    enddo
    
  end subroutine crght_top_area_Properties

  

  subroutine clft_bot_area_Properties
    implicit none
    integer :: i
    real*8  :: prll1,prll2
    
    count_area = ncl + ncr + ncclt + nccrt
    
    do i = 1,ncclb
       
       if (i.eq.1) then
          count_area = count_area+1
          typ_area(count_area) = 3
          k_area(count_area)   = 1.0d0
          
          prll1 = Lb1*(deflc_frm_pulley/2.0d0)
          prll2 = Lb1*(deflc_frm_cntrl_top_blck/2.0d0)

          A0(count_area) = 0.50d0*(prll1+prll2)*&
               (Lb2*unt_len-La1*unt_len_cb)

          !0.5*(0.1+0.5)*0.9 [Trapezoid]
          !write(*,*) A0(count_area),"A0"
          
          
          fctr(count_area)     = 1.0d0
          beta(count_area)     = 1.0d0
          
       elseif (i.eq.2) then
          count_area = count_area+2
          typ_area(count_area) = 4
          k_area(count_area)   = 1.0d0
          
          prll1 = Lb1*(deflc_frm_pulley/2.0d0)
          prll2 = Lb1*(deflc_frm_cntrl_top_blck/2.0d0)

          A0(count_area) = 0.50d0*(prll1+prll2)*&
               (Lb2*unt_len-La1*unt_len_cb)
          ! 0.5*(0.1+0.5)*0.9 [Trapezoid]
          !write(*,*) A0(count_area),"A0"
          
          fctr(count_area)     = 1.0d0
          beta(count_area)     = 1.0d0
          
       else
          count_area = count_area + 2
          k_area(count_area)   = 1.0d0
          
          if (stageNo==4 .and. stageType==1) then

             if (i.ne.(ncclb)) then
                typ_area(count_area) = 5
                A0(count_area) = (Lb2*unt_len-La1*unt_len_cb)*La2
                
             elseif (i==(ncclb)) then
                typ_area(count_area) = 6
                prll1 = La2 
                prll2 = deflc_for_endRgn
                
                A0(count_area) = 0.50*(prll1+prll2)*(Lb2*unt_len-La1*unt_len_cb)
             endif
             
          elseif (stageNo==4 .and. stageType==2) then
             typ_area(count_area) = 5
             A0(count_area)       = (Lb2*unt_len-La1*unt_len_cb)*La2
             !write(*,*) A0(count_area),"A0"
          endif
          
          fctr(count_area)     = 1.0d0
          beta(count_area)     = 1.0d0
          
       endif
       
       At(count_area) = A0(count_area)
       area_sides(count_area) = 4
       
    enddo
    
  end subroutine clft_bot_area_Properties



  subroutine crght_bot_area_Properties
    implicit none
    integer :: i
    real*8  :: prll1,prll2
    
   count_area = ncl + ncr + ncclt + nccrt + 1
   
   do i = 1,nccrb
      
      if (i.eq.1) then
         count_area = count_area+1
         typ_area(count_area) = 3
         k_area(count_area)   = 1.0d0
         
         prll1 = Lb1*(deflc_frm_pulley/2.0d0)
         prll2 = Lb1*(deflc_frm_cntrl_top_blck/2.0d0)
         
         A0(count_area) = 0.50d0*(prll1+prll2)*&
              (Lb2*unt_len-La1*unt_len_cb)
         
         ! 0.5*(0.1+0.5)*0.9 [Trapezoid]
         !write(*,*) A0(count_area),"A0"
         
         
         fctr(count_area)     = 1.0d0
         beta(count_area)     = 1.0d0
         
      elseif (i.eq.2) then
         count_area = count_area+2
         typ_area(count_area) = 4
         k_area(count_area)   = 1.0d0

         prll1 = Lb1*(deflc_frm_pulley/2.0d0)
         prll2 = Lb1*(deflc_frm_cntrl_top_blck/2.0d0)
         
         A0(count_area) = 0.50d0*(prll1+prll2)*&
              (Lb2*unt_len-La1*unt_len_cb)
         
         ! 0.5*(0.1+0.5)*0.9 [Trapezoid]
         !write(*,*) A0(count_area),"A0"
         
         
         fctr(count_area)     = 1.0d0
         beta(count_area)     = 1.0d0
         
      else
         count_area = count_area+2
         k_area(count_area)   = 1.0d0

         if (stageNo==4 .and. stageType==1) then
            
            if (i.ne.(ncclb)) then
               typ_area(count_area) = 5
               A0(count_area) = (Lb2*unt_len-La1*unt_len_cb)*La2
               
            elseif (i==(ncclb)) then
               typ_area(count_area) = 6
               prll1 = La2 
               prll2 = deflc_for_endRgn
               
               A0(count_area) = 0.50*(prll1+prll2)*(Lb2*unt_len-La1*unt_len_cb)
            endif
            
         elseif (stageNo==4 .and. stageType==2) then
            typ_area(count_area) = 5
            A0(count_area)       = (Lb2*unt_len-La1*unt_len_cb)*La2
            !write(*,*) A0(count_area),"A0"
         endif
         
         fctr(count_area)     = 1.0d0
         beta(count_area)     = 1.0d0
         
      endif
      
      At(count_area) = A0(count_area)
      area_sides(count_area) = 4
      
   enddo
   
 end subroutine crght_bot_area_Properties
 
 subroutine endRgnCell_area_Properties
   implicit none
   integer :: i

   do i = 1,ncer
      
      count_area = count_area+1
      !write(*,*) count_area,"ca"
      typ_area(count_area)   = 7
      area_sides(count_area) = 3
      k_area(count_area)     = 1.0d0
      
      A0(count_area) = 0.5* (2.0d0 * (Lb2*unt_len-La1*unt_len_cb)) * (La2-deflc_for_endRgn)
      At(count_area)       = A0(count_area)!not needed, simlrty with othr rtns
      fctr(count_area)     = 1.0d0
      beta(count_area)     = 1.0d0 !not needed, simlrty with othr rtns
      
   enddo
   
 end subroutine endRgnCell_area_Properties
 
 subroutine get_SideWiseCellNo
   implicit none
   integer :: i

   open(unit=156,file='LftRghtSideCell.dat',position='append')
   
   call get_lftsideCellNo
   call get_rghtsideCellNo
   call get_cntrlCellNo
   call get_endRgnCellNo
   
   do i = 1,N_lftCells
      if (i==1) write(156,fmt=*) "LeftSide Cell List :"
      write(156,fmt=*) lftSideCell(i),i,"lftSideCell"
   enddo
   
   do i = 1,N_rghtCells
      if (i==1) write(156,fmt=*) "RightSide Cell List :"
      write(156,fmt=*) rghtSideCell(i),i,"rghtSideCell"
   enddo
   
   do i = 1,N_cntrlCells
      if (i==1) write(156,fmt=*) "Central Cell List :"
      write(156,fmt=*) cntrlCell(i),i,"cntrlCell"
   enddo

   do i = 1,N_endRgnCells
      if (i==1) write(156,fmt=*) "EndRegion Cell List :"
      write(156,fmt=*) endRgnCell(i),i,"endRgnCell"
   enddo
   
   close(156)
   
 end subroutine get_SideWiseCellNo
 
 subroutine get_lftsideCellNo
   implicit none
   integer :: cnt_cell,cnt_lft
   integer :: i
   
   cnt_cell = 0 ; cnt_lft = 0

   do i = 1,ncl
      lftSideCell(cnt_lft+i) = cnt_cell + i
   enddo

   cnt_cell = ncl + ncr
   cnt_lft  = ncl
   
   do i = 1,ncclt
      lftSideCell(cnt_lft+i) = cnt_cell + i
   enddo

   cnt_cell = cnt_cell + ncclt + nccrt
   cnt_lft  = cnt_lft + ncclt
   
   do i = 1,ncclb
      lftSideCell(cnt_lft+i) = cnt_cell + (2*i-1)
   enddo

   cnt_lft = cnt_lft + ncclb
   
   if (cnt_lft .ne. N_lftCells) then
      write(*,*)"flnm:info_modules,sbrtn:get_lftsideCell"
      stop
   endif
   
 end subroutine get_lftsideCellNo
 
 subroutine get_rghtsideCellNo
   implicit none
   integer :: cnt_cell,cnt_rght
   integer :: i
   
   cnt_cell = ncl ; cnt_rght = 0
   
   do i = 1,ncr
      rghtSideCell(cnt_rght+i) = cnt_cell + i
   enddo
   
   cnt_cell = ncl + ncr + ncclt
   cnt_rght = ncr
   
   do i = 1,nccrt
      rghtSideCell(cnt_rght+i) = cnt_cell + i
   enddo
   
   cnt_cell = cnt_cell + nccrt + 1
   cnt_rght = cnt_rght + nccrt
   
   do i = 1,nccrb
      rghtSideCell(cnt_rght+i) = cnt_cell + (2*i-1)
   enddo
   
   cnt_rght = cnt_rght + nccrb
   
   if (cnt_rght .ne. N_rghtCells) then
      write(*,*)"flnm:info_modules,sbrtn:get_rghtsideCell"
      stop
   endif
   
 end subroutine get_rghtsideCellNo
 
 
 subroutine get_cntrlCellNo
   implicit none
   
   if (stageNo==1 .and. stageType==1) then
      cntrlCell(1) = (ncl+ncr) + 1
   else
      write(*,*) "Check get_cntrlCellNo in info_modules.f08"
      call sleep(1)
   endif
   
 end subroutine get_cntrlCellNo
 
 
 subroutine get_endRgnCellNo
   implicit none
   
   if (stageNo==4 .and. stageType==1) then
      endRgnCell(1) = (ncl+ncr) + (ncclt+nccrt) + (ncclb+nccrb) + 1
   else
      write(*,*) "Check get_endRgnCellNo in info_modules.f08"
      call sleep(1)
   endif
   
 end subroutine get_endRgnCellNo
 
 subroutine get_N_optmCell_and_optmCell
   implicit none
   integer :: i
   
   N_optmCell = N_cell
   
   do i = 1,N_cell
      optmCell(i) = 1
      !write(*,*) optmCell(i),i,"optmCell,CellNo"
   enddo
   
 end subroutine get_N_optmCell_and_optmCell
 
 subroutine print_area_variables
   implicit none
   integer :: i

   open(unit=145,file='Area_Variables.dat')
   
   write(145,fmt=*) "area_type","A0","k_area","At","beta","fctr","area_nm"
   
   do i = 1,N_cell
      write(145,fmt=*) typ_area(i),A0(i),k_area(i),At(i),beta(i),fctr(i),i
   enddo

   write(145,fmt=*) ""
   close(145)
   
 end subroutine print_area_variables
 
 subroutine print_area_variables_only_AA0ka
   implicit none
   integer :: i

   open(unit=245,file='Area_Variables_onlyAA0ka.dat')
   
   write(245,fmt=*) "A","A0","k_area","area_nm"
   
   do i = 1,N_cell
      write(245,fmt=*) A(i),A0(i),k_area(i),i
   enddo
   
   write(245,fmt=*) ""
   close(245)
   
 end subroutine print_area_variables_only_AA0ka
 
 
end module area_variables


module bend_variables
  use sys_building_info
  use system_parameters

  implicit none

contains

  subroutine allocate_and_initialize_bend_variables
    implicit none
    
    max_Phi_node = 4

    allocate(nodePhi_typ(1:N_node))
    allocate(nodePhi_typStr(1:N_node),nodePhi_typStrAN(1:N_node),nodePhi_typStrES(1:(N_node+2)))
    allocate(nodePhi_typStrNVFR(1:N_node))
    
    allocate(k_phi(1:N_node,1:max_Phi_node))
    allocate(k_phiStr(1:N_node,1:max_Phi_node))
    allocate(k_phiStrAN(1:N_node,1:max_Phi_node),k_phiStrES(1:N_node,1:max_Phi_node))
    allocate(k_phiStrNVFR(1:N_node,1:max_Phi_node))
    
    nodePhi_typ(1:N_node)        = -10
    nodePhi_typStr(1:N_node)     = -10
    nodePhi_typStrAN(1:N_node)   = -10
    nodePhi_typStrES(1:N_node)   = -10 !unrealistic
    nodePhi_typStrNVFR(1:N_node) = -10
    
    k_phi(1:N_node,1:max_Phi_node)        = -1.0d30
    k_phiStr(1:N_node,1:max_Phi_node)     = -1.0d30
    k_phiStrAN(1:N_node,1:max_Phi_node)   = -1.0d30
    k_phiStrES(1:N_node,1:max_Phi_node)   = -1.0d30
    k_phiStrNVFR(1:N_node,1:max_Phi_node) = -1.0d30
    
  end subroutine allocate_and_initialize_bend_variables
  
  subroutine allocate_and_initialize_bend_variables_wo_StrVars
    implicit none

    allocate(nodePhi_typ(1:N_node))
    allocate(k_phi(1:N_node,1:max_Phi_node))

    nodePhi_typ  = -10 !unrealistic
    
    k_phi(1:N_node,1:max_Phi_node) = -1.0d30
    
  end subroutine allocate_and_initialize_bend_variables_wo_StrVars
  
  subroutine deallocate_bend_variables
    implicit none

    deallocate(nodePhi_typ,nodePhi_typStr,nodePhi_typStrAN)
    deallocate(k_phi,k_phiStr,k_phiStrAN)
    
  end subroutine deallocate_bend_variables
  
  subroutine get_bend_variable_values
    implicit none
    integer :: i
    
    nodePhi_typ(1:N_node) = 1
    k_phi(1:N_node,1:max_Phi_node) = 0.0d0
    
  end subroutine get_bend_variable_values
  
  subroutine print_bend_variables
    implicit none
    integer :: i

    open(unit=145,file='Bend_Variables.dat')

    write(145,fmt=*) "nodePhi_typ","   ","node_nm"
    
    do i = 1,N_node
       write(145,fmt=*) nodePhi_typ(i),i 
    enddo
    
    write(145,fmt=*) " "
    
    write(145,fmt=*) "k_phi","node_nm"
    
    do i = 1,N_node
       write(145,fmt=*) k_phi(i,1:max_Phi_node),i
    enddo
    
    write(145,fmt=*) " "
    
    close(145)
    
  end subroutine print_bend_variables
  
end module bend_variables


module grvtnl_variables
  use sys_building_info
  use cell_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains
  
  subroutine grvtnl_variables_values(N_imp_region_dum,div_imp_region_dum,imp_region_dum)
    
    implicit none
    integer, intent(in) :: N_imp_region_dum
    integer, intent(in) :: div_imp_region_dum
    integer, intent(in) :: imp_region_dum(1:N_imp_region_dum)
    
    !range change (CgX_finish,CgY_finish)
    !intermidiate gap (prim_CgXNodes change)
    
    CgX = 0.0d0 ; CgX_max = 0.0d0 !CgX_max will be valued later()
    CgX_strt = 0.00d0 
    CgX_fnsh = 0.25d0  
    
    !N_CgXEval = -1
    N_region = 5
    
    !N_imp_region = 1
    !allocate(imp_region(1:N_imp_region))
    !imp_region(1) = 1
    
    prim_CgXNodes = N_region + 1
    
    !div_imp_region = 100
    addtnl_CgXNodes = N_imp_region_dum*(div_imp_region_dum-1)
    
    N_CgXEval = prim_CgXNodes + addtnl_CgXNodes ; N_CgYEval = N_CgXEval
    write(*,*) N_CgXEval,N_CgYEval,"N_CgXEval-N_CgYEval"
    
    allocate(CgX_val(1:N_CgXEval),CgY_val(1:N_CgYEval))
    CgX_val = 0.00d0 ; CgY_val = 0.00d0
    
    count_cgX = 1
    
    region_incr = (CgX_fnsh-CgX_strt)/(N_region)
    imp_region_incr = (region_incr)/(div_imp_region_dum)
    
    
  end subroutine grvtnl_variables_values
  
  
  subroutine rgion_dvsn(N_imp_region_dum,div_imp_region_dum,imp_region_dum)
    implicit none
    integer, intent(in) :: N_imp_region_dum
    integer, intent(in) :: div_imp_region_dum
    integer, intent(in) :: imp_region_dum(1:N_imp_region_dum)
    
    integer :: i,j
    
    call grvtnl_variables_values(N_imp_region_dum,div_imp_region_dum,imp_region_dum)
    
    do i = 1,N_region
       CgX_val(count_cgX) = CgX_strt + (i-1)*region_incr
       
       do j = 1,N_imp_region_dum
          
          if (i==imp_region_dum(j)) then             
             imp_region_CgX_strt = CgX_val(count_cgX)
             call imp_rgn_CgXs
             !count_cgX = count_cgX + (div_imp_region-1)
          endif
          
       enddo
       
       if (i.eq.N_region) then
          CgX_val(count_cgX+1) = CgX_fnsh
       endif
       
       count_cgX = count_cgX + 1
       
    enddo
    
    call print_grvtnl_variables
   

  contains
    
    subroutine imp_rgn_CgXs
      
      implicit none
      integer :: i1
      !write(*,*) count_cgX,"inside_imp_regn_Cgs"
      count_cgX = count_cgX + 1
      
      !write(*,*) div_imp_region,"div_imp_region"
      !write(*,*) imp_region_CgX_strt,imp_region_incr,"imp_region_CgX_strt and incr"
      
      do i1 = 1,(div_imp_region-1)
         CgX_val(count_cgX) = imp_region_CgX_strt + i1*imp_region_incr
         
         if (i1.ne.(div_imp_region-1)) then
            count_cgX = count_cgX + 1
         endif
      enddo
      
    end subroutine imp_rgn_CgXs

  end subroutine rgion_dvsn

  
  
  subroutine print_grvtnl_variables
    implicit none
    integer :: i
    
    !do i = 1,N_CgXEval
       !write(*,*) CgX_val(i),CgY_val,"CgX_val,CgY_val"
    !enddo
    
  end subroutine print_grvtnl_variables
  
  
  subroutine allocate_and_initialize_grvVars
    implicit none
    
    allocate(CgXNode(1:N_node),CgXNode_Str(1:N_node),CgXNode_StrAN(1:N_node),CgXNode_StrES(1:(N_node+2)))
    allocate(CgXNode_StrNVFR(1:N_node))
    allocate(CgYNode(1:N_node),CgYNode_Str(1:N_node),CgYNode_StrAN(1:N_node),CgYNode_StrES(1:(N_node+2)))
    allocate(CgYNode_StrNVFR(1:N_node))
    
    CgXNode         = -1.0d30 ; CgYNode         = -1.0d30
    CgXNode_Str     = -1.0d30 ; CgYNode_Str     = -1.0d30
    CgXNode_StrAN   = -1.0d30 ; CgYNode_StrAN   = -1.0d30
    CgXNode_StrES   = -1.0d30 ; CgYNode_StrES   = -1.0d30
    CgXNode_StrNVFR = -1.0d30 ; CgYNode_StrNVFR = -1.0d30
    
  end subroutine allocate_and_initialize_grvVars
  
  subroutine allocate_and_initialize_grvVars_wo_StrVars
    implicit none
    
    allocate(CgXNode(1:N_node),CgYNode(1:N_node))
    CgXNode = -1.0d30 
    CgYNode = -1.0d30
    
  end subroutine allocate_and_initialize_grvVars_wo_StrVars
  
  subroutine get_CgNode
    
    implicit none
    integer :: i,j,jmax
    integer :: cnt_node,cnt_strt

    open(unit=121,file='Cg_Node.dat')
    
    CgXNode(1:N_node) = 0.0d0
    CgYNode(1:N_node) = 0.0d0
    
    do i = 1,2
       
       if (i==1) then
          jmax=nvsl
          cnt_strt = 0
       elseif (i==2) then
          jmax=nvsr
          cnt_strt = 2*nvsl
       endif
       
       write(*,*) jmax,"jmax in CgNode"
       
       do j = 1,jmax
          
          if (j.ne.1) then
             continue
          else
             cnt_node = cnt_strt + (2*j-1)
             
             if (i==1) then
                CgXNode(cnt_node)   = 0.00d0
                CgXNode(cnt_node+1) = 0.00d0
                
             elseif (i==2) then
                CgXNode(cnt_node)   = -0.00d0
                CgXNode(cnt_node+1) = -0.00d0
             endif
             
          endif
          
       enddo
       
    enddo
    
    do i = 1,N_node
       write(unit=121,fmt=*) "CgXNode=",CgXNode(i),"CgYNode=",CgYNode(i),"for i=",i
    enddo
    
    close(121)
    
  end subroutine get_CgNode
  
  subroutine get_select_xy
    implicit none
    
    !select_xy = 1
    !select_xy = 2
    select_xy = 3
    
  end subroutine get_select_xy
  
  subroutine deallocate_grvVars
    implicit none
    
    deallocate(CgXNode,CgXNode_Str,CgXNode_StrAN)
    deallocate(CgYNode,CgYNode_Str,CgYNode_StrAN)
    
  end subroutine deallocate_grvVars
  
  
  subroutine print_grvtnl_vars
    implicit none
    integer :: i
    
    open(unit=145,file='Gravitational_Variables.dat')
    write(145,*) "CgXNode","CgYNode","node_nm"
    
    do i = 1,N_node
       write(145,*) CgXNode(i),CgYNode(i),i
    enddo
    
    write(145,*) ""
    close(145)
    
  end subroutine print_grvtnl_vars
  
end module grvtnl_variables

module SRyp_variables
  use sys_building_info
  use cell_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains
  
  subroutine allocate_and_initialize_SRypVars
    implicit none
    
    write(*,*) numApclNodes,"numApclNodes CHK HERE"
    
    allocate(AmpApclYP(1:numApclNodes))           ; AmpApclYP(1:numApclNodes)           = -1.0d30 
    allocate(AlphaApclYP(1:numApclNodes))         ; AlphaApclYP(1:numApclNodes)         = -1.0d30 
    allocate(epsApclYP(1:numApclNodes))           ; epsApclYP(1:numApclNodes)           = -1.0d30
    allocate(powrApclYP(1:numApclNodes))          ; powrApclYP(1:numApclNodes)          = -1.0d30
    allocate(activtnFctrApcl(1:numApclNodes))     ; activtnFctrApcl(1:numApclNodes)     = -1
    
    allocate(AmpApclYP_Str(1:numApclNodes))       ; AmpApclYP_Str(1:numApclNodes)       = -1.0d30 
    allocate(AlphaApclYP_Str(1:numApclNodes))     ; AlphaApclYP_Str(1:numApclNodes)     = -1.0d30 
    allocate(epsApclYP_Str(1:numApclNodes))       ; epsApclYP_Str(1:numApclNodes)       = -1.0d30
    allocate(powrApclYP_Str(1:numApclNodes))      ; powrApclYP_Str(1:numApclNodes)      = -1.0d30
    allocate(activtnFctrApcl_Str(1:numApclNodes)) ; activtnFctrApcl_Str(1:numApclNodes) = -1
    
  end subroutine allocate_and_initialize_SRypVars
  
  subroutine allocate_and_initialize_SRypVars_wo_StrVars
    implicit none
    
    allocate(AmpApclYP(1:numApclNodes))           ; AmpApclYP(1:numApclNodes)           = -1.0d30 
    allocate(AlphaApclYP(1:numApclNodes))         ; AlphaApclYP(1:numApclNodes)         = -1.0d30 
    allocate(epsApclYP(1:numApclNodes))           ; epsApclYP(1:numApclNodes)           = -1.0d30
    allocate(powrApclYP(1:numApclNodes))          ; powrApclYP(1:numApclNodes)          = -1.0d30
    allocate(activtnFctrApcl(1:numApclNodes))     ; activtnFctrApcl(1:numApclNodes)     = -1
    
  end subroutine allocate_and_initialize_SRypVars_wo_StrVars
  
  subroutine get_SRYP_props
    implicit none
    
    AmpApclYP(1:numApclNodes)       = 0.00010d0
    AlphaApclYP(1:numApclNodes)     = 30.00d0   !250.0000d0
    epsApclYP(1:numApclNodes)       = 0.300d0
    powrApclYP(1:numApclNodes)      = 2.00d0
    activtnFctrApcl(1:numApclNodes) = 0
    
  end subroutine get_SRYP_props
  
  subroutine get_positionalNodes()
    implicit none
    integer :: i,i_fctr,j,jmax,k,kmax
    integer :: numApclInACell,numBsalInACell,numLtrlInACell
    integer :: apclCnt,bsalCnt,ltrlCnt
    integer :: cntInsrtdNode
    
    if (modelID==1)         stop 'Not written for modelID=1'
    if (NI_incldAS_woPP==0) stop 'Not written for IN only in basal and lateral sides' 
    
    if (VF_regionModelled==0) Hlf_Ncell = N_cell/2
    if (VF_regionModelled==1) Hlf_Ncell =(N_cell-1)/2
    
    write(*,*) modelID,NI_incldAS_woPP,VF_regionModelled,Hlf_Ncell,"props get_positionalNodes"
    
    numRegulrNode = 2*((Hlf_Ncell+1)*2)
    numInsrtdNode = N_node - numRegulrNode; write(*,*)numRegulrNode,numInsrtdNode,N_node,"numNodes"
    
    numApclNodes = (Hlf_Ncell+1)*(2) + (Hlf_Ncell*2)*(NAEC_Apcl) + (NAEC_Apcl)
    numBsalNodes = (Hlf_Ncell+1)*(2) + (Hlf_Ncell*2)*(NAEC_Bsal) + (NAEC_Bsal)
    numLtrlNodes = (Hlf_Ncell*2)*(NAEC_Ltrl)
    
    write(*,*) numApclNodes,numBsalNodes,numLtrlNodes,"num of Apcl Bsal and Ltrl Nodes"
    write(*,*) " "
    
    allocate(apclNodes(1:numApclNodes),bsalNodes(1:numBsalNodes),ltrlNodes(1:numLtrlNodes))
    apclNodes=-1 ; bsalNodes=-1 ; ltrlNodes=-1
    
    numApclInACell=1+NAEC_Apcl ; numBsalInACell=1+NAEC_Bsal ; numLtrlInACell=NAEC_Ltrl
    cntInsrtdNode = 0
    apclCnt=0 ; bsalCnt=0 ; ltrlCnt=0
    
    do i = 1,(2*Hlf_Ncell)
       
       if (i.le.Hlf_Ncell) i_fctr=i-1
       if (i.gt.Hlf_Ncell) i_fctr=i-0
       
       if (i.le.Hlf_Ncell) then
          apclCnt = 1+i_fctr*numApclInACell
          bsalCnt = 1+i_fctr*numBsalInACell
       elseif ((i.gt.Hlf_Ncell) .and. (i.le.(2*Hlf_Ncell))) then
          apclCnt = 1+(i_fctr-1)*(numApclInACell)+1
          bsalCnt = 1+(i_fctr-1)*(numBsalInACell)+1
       endif
       
       ltrlCnt = (i-1)*(NAEC_Ltrl)
       
       apclNodes(apclCnt) = i_fctr*2 + 1
       bsalNodes(bsalCnt) = i_fctr*2 + 2
       
       write(*,*) apclCnt,bsalCnt,apclNodes(apclCnt),bsalNodes(bsalCnt),"apcl-bsal Rglr"
       write(*,*) " "
       
       do j = 1,3
          
          if (j==1) kmax=NAEC_Apcl
          if (j==2) kmax=NAEC_Bsal
          if (j==3) kmax=NAEC_Ltrl
          
          do k = 1,kmax
             
             if (j==1) apclNodes(apclCnt+k) = numRegulrNode + cntInsrtdNode + 1
             if (j==2) bsalNodes(bsalCnt+k) = numRegulrNode + cntInsrtdNode + 1
             if (j==3) ltrlNodes(ltrlCnt+k) = numRegulrNode + cntInsrtdNode + 1
             
             cntInsrtdNode = cntInsrtdNode+1
             if (j==1) write(*,*) apclCnt,apclCnt+k,cntInsrtdNode,"cntInsrtdNode A"
             if (j==2) write(*,*) bsalCnt,bsalCnt+k,cntInsrtdNode,"cntInsrtdNode B"
             if (j==3) write(*,*) ltrlCnt,ltrlCnt+k,cntInsrtdNode,"cntInsrtdNode L"
             
          enddo
          
          if (j==1) write(*,*) apclNodes((apclCnt+1):(apclCnt+NAEC_Apcl)),"apcl IN"
          if (j==2) write(*,*) bsalNodes((bsalCnt+1):(bsalCnt+NAEC_Bsal)),"bsal IN"
          if (j==3) write(*,*) ltrlNodes((ltrlCnt+1):(ltrlCnt+NAEC_Ltrl)),"ltrl IN"
          
       enddo
       
       if (i==Hlf_Ncell)   apclNodes(1+(i_fctr+1)*(numApclInACell)) = (i_fctr+1)*2 + 1
       if (i==Hlf_Ncell)   bsalNodes(1+(i_fctr+1)*(numBsalInACell)) = (i_fctr+1)*2 + 2
       if (i==2*Hlf_Ncell) apclNodes(2+(i_fctr+0)*(numApclInACell)) = (i_fctr+1)*2 + 1
       if (i==2*Hlf_Ncell) bsalNodes(2+(i_fctr+0)*(numBsalInACell)) = (i_fctr+1)*2 + 2
       
       if (i==Hlf_Ncell)     write(*,*) (1+(i_fctr+1)*(numApclInACell)),(1+(i_fctr+1)*(numBsalInACell)),&
            apclNodes(1+(i_fctr+1)*(numApclInACell)),bsalNodes(1+(i_fctr+1)*(numBsalInACell))
       if (i==(2*Hlf_Ncell)) write(*,*) (2+(i_fctr+0)*(numApclInACell)),(2+(i_fctr+0)*(numBsalInACell)),&
            apclNodes(2+(i_fctr+0)*(numApclInACell)),bsalNodes(2+(i_fctr+0)*(numBsalInACell))
       
       write(*,*) " "
       
       
       if (i==2*Hlf_Ncell) then ! add for INITIATOR CELL
          
          do j = 1,2
             if (j==1) kmax=NAEC_Apcl
             if (j==2) kmax=NAEC_Bsal
             do k = 1,kmax
                if (j==1) apclNodes(apclCnt+NAEC_Apcl+1+k) = numRegulrNode + cntInsrtdNode + 1
                if (j==2) bsalNodes(bsalCnt+NAEC_Bsal+1+k) = numRegulrNode + cntInsrtdNode + 1
                cntInsrtdNode = cntInsrtdNode+1
                if (j==1) write(*,*) apclCnt,apclCnt+NAEC_Apcl+1+k,cntInsrtdNode,"cntInsrtdNode A-end"
                if (j==2) write(*,*) bsalCnt,bsalCnt+NAEC_Bsal+1+k,cntInsrtdNode,"cntInsrtdNode B-end"
             enddo
             if (j==1)write(*,*) apclNodes((apclCnt+NAEC_Apcl+1+1):(apclCnt+1+2*NAEC_Apcl)),"apIN-end"
             if (j==2)write(*,*) bsalNodes((bsalCnt+NAEC_Bsal+1+1):(bsalCnt+1+2*NAEC_Bsal)),"bsIN-end"
          enddo
          
       endif
       
    enddo
    
    
    
    do i = 1,3
       
       if (i==1) jmax=numApclNodes
       if (i==2) jmax=numBsalNodes
       if (i==3) jmax=numLtrlNodes
       
       do j = 1,jmax
          if (i==1) write(*,*) j,apclNodes(j),"AP"
          if (i==2) write(*,*) j,bsalNodes(j),"BS"
          if (i==3) write(*,*) j,ltrlNodes(j),"LT"
       enddo
       
       write(*,*) " "
    enddo
    
  end subroutine get_positionalNodes
  
  
  subroutine get_positional_Nodes_with_diff_types_of_InsrtdNodes
    implicit none
    integer :: i,imax,i_fctr,j,jmax,k,kmax
    integer :: numApclInACell,numBsalInACell,numLtrlInACell
    integer :: apclCnt,bsalCnt,ltrlCnt
    integer :: apclCntLft, bsalCntLft, ltrlCntLft
    integer :: apclCntRght,bsalCntRght,ltrlCntRght
    integer :: cntInsrtdNode,RglrNodeLft,InsrtdNodeLft,lftNode,rghtNode
    
    if (modelID==1)           stop 'Not written for modelID=1'
    if (NI_incldAS_woPP==0)   stop 'not written for NI_incldAS_woPP=0--NI basal and lateral surface' 
    if (NI_AS_wtCrtclSurfc==0) stop 'not written for NI_AS_wtCrtclSurfc=0--without crtcl surfc'
    
    if (VF_regionModelled==0) Hlf_Ncell = N_cell/2
    if (VF_regionModelled==1) Hlf_Ncell =(N_cell-1)/2
    
    write(*,*) modelID,NI_incldAS_woPP,VF_regionModelled,Hlf_Ncell,"props get_positionalNodes"
    
    numRegulrNode = 2*((Hlf_Ncell+1)*2)
    numInsrtdNode = N_node - numRegulrNode; write(*,*)numRegulrNode,numInsrtdNode,N_node,"numNodes"
    
    numApclNodes  = (Hlf_Ncell+1)*(2) + (Hlf_Ncell-NCP_CrtclApSrfc)*(2)*(NAEC_Apcl) &
         + (NCP_CrtclApSrfc*2)*(NAEC_ApclCrtcl) + (NAEC_ApclCrtcl)
    
    numBsalNodes  = (Hlf_Ncell+1)*(2) + (Hlf_Ncell*2)*(NAEC_Bsal) + (NAEC_Bsal)
    numLtrlNodes  = (Hlf_Ncell*2)*(NAEC_Ltrl)
    
    write(*,*) numApclNodes,numBsalNodes,numLtrlNodes,"num of Apcl Bsal and Ltrl Nodes" ; write(*,*) " "
    
    call chk_allocatn_and_allocate_the_srfcNodes
    apclNodes=-1 ; bsalNodes=-1 ; ltrlNodes=-1
    
    numApclInACell = 1+NAEC_Apcl ; numBsalInACell = 1+NAEC_Bsal ; numLtrlInACell = NAEC_Ltrl
    cntInsrtdNode  = 1
    apclCnt        = 1 ; bsalCnt = 1 ; ltrlCnt = 1
    
    
    do i = 1,Hlf_Ncell
       
       apclNodes(apclCnt) = (i-1)*2 +1 ; apclCnt = apclCnt+1
       bsalNodes(bsalCnt) = (i-1)*2 +2 ; bsalCnt = bsalCnt+1
       
       do j = 1,3
          
          if (j==1) then
             if (i.le.(Hlf_Ncell-NCP_CrtclApSrfc)) kmax = NAEC_Apcl
             if (i.gt.(Hlf_Ncell-NCP_CrtclApSrfc)) kmax = NAEC_ApclCrtcl
          elseif (j==2) then
             kmax = NAEC_Bsal
          elseif (j==3) then
             kmax = NAEC_Ltrl
          endif
          
          do k = 1,kmax
             
             if (j==1) then
                apclNodes(apclCnt) = numRegulrNode + cntInsrtdNode
                write(*,*) apclCnt,cntInsrtdNode,apclNodes(apclCnt),"A insr"
                apclCnt = apclCnt+1 ; cntInsrtdNode = cntInsrtdNode+1
                
             elseif (j==2) then
                bsalNodes(bsalCnt) = numRegulrNode + cntInsrtdNode
                write(*,*) bsalCnt,cntInsrtdNode,bsalNodes(bsalCnt),"B insr"
                bsalCnt = bsalCnt+1 ; cntInsrtdNode = cntInsrtdNode+1
                
             elseif (j==3) then
                ltrlNodes(ltrlCnt) = numRegulrNode + cntInsrtdNode
                write(*,*) ltrlCnt,cntInsrtdNode,ltrlNodes(ltrlCnt),"L insr"
                ltrlCnt = ltrlCnt+1 ; cntInsrtdNode = cntInsrtdNode+1
             endif
             
          enddo
          
       enddo
       
    enddo
    
    apclNodes(apclCnt) = (Hlf_Ncell)*2 + 1 ; bsalNodes(bsalCnt) = (Hlf_Ncell)*2 + 2 
    write(*,*) apclCnt,bsalCnt,apclNodes(apclCnt),bsalNodes(bsalCnt),"A B p1"
    
    apclCntLft    = apclCnt ; bsalCntLft = bsalCnt ; ltrlCntLft = ltrlCnt-1
    apclCnt = apclCnt+1 ; bsalCnt = bsalCnt+1 ; write(*,*) apclCnt,bsalCnt,"A B p2"
    
    RglrNodeLft   = ((Hlf_Ncell+1)*2)
    InsrtdNodeLft = (Hlf_Ncell-NCP_CrtclApSrfc)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) &
         + (NCP_CrtclApSrfc)*(NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl)
    
    write(*,*) apclCntLft,bsalCntLft,ltrlCntLft,RglrNodeLft,InsrtdNodeLft,"cnt Lft"
    
    
    do j = 1,3
       
       if (j==1) kmax = apclCntLft
       if (j==2) kmax = bsalCntLft
       if (j==3) kmax = ltrlCntLft
       
       do k = 1,kmax
          
          if (j==1) lftNode = apclNodes(k)
          if (j==2) lftNode = bsalNodes(k)
          if (j==3) lftNode = ltrlNodes(k)
          
          if (lftNode .le. numRegulrNode) rghtNode = lftNode+RglrNodeLft
          if (lftNode .gt. numRegulrNode) rghtNode = lftNode+InsrtdNodeLft
          
          if (j==1) apclNodes(apclCntLft+k) = rghtNode
          if (j==2) bsalNodes(bsalCntLft+k) = rghtNode
          if (j==3) ltrlNodes(ltrlCntLft+k) = rghtNode

          if (j==1) write(*,*) (apclCntLft+k),lftNode,rghtNode,j,k,"LN"
       enddo
       
    enddo
    
    apclCnt = apclCnt + apclCntLft
    bsalCnt = bsalCnt + bsalCntLft
    ltrlCnt = ltrlCnt + ltrlCntLft
    
    apclNodes(apclCnt) = 2*(Hlf_Ncell+1)*2 - 1 ; bsalNodes(bsalCnt) = 2*(Hlf_Ncell+1)*2 
    write(*,*) apclCnt,bsalCnt,apclNodes(apclCnt),bsalNodes(bsalCnt),"A B p1"
    apclCnt = apclCnt+1 ; bsalCnt = bsalCnt+1 ; write(*,*) apclCnt,bsalCnt,"A B p2"
    
    apclCntRght = 2*apclCntLft ; bsalCntRght = 2*bsalCntLft ; ltrlCntRght = 2*ltrlCntLft
    write(*,*) apclCntRght,bsalCntRght,ltrlCntRght,"apcl-bsal-ltrl Cnt Rght"
    
    do j = 1,2
       if (j==1) kmax = NAEC_ApclCrtcl
       if (j==2) kmax = NAEC_Bsal
       
       do k = 1,kmax
          if (j==1) apclNodes(apclCntRght+k) = (apclCntRght+bsalCntRght+ltrlCntRght) + k
          if (j==2) bsalNodes(bsalCntRght+k) = (apclCntRght+bsalCntRght+ltrlCntRght+NAEC_ApclCrtcl) + k
          if (j==1) write(*,*) j,k,(apclCntRght+k),apclNodes(apclCntRght+k),"InsrCRT j1"
          if (j==2) write(*,*) j,k,(bsalCntRght+k),bsalNodes(bsalCntRght+k),"InsrCRT j2"   
       enddo
       
    enddo

    do i = 1,3
       
       if (i==1) jmax=numApclNodes
       if (i==2) jmax=numBsalNodes
       if (i==3) jmax=numLtrlNodes
       
       do j = 1,jmax
          if (i==1) write(*,*) apclNodes(j),j,"apcl wri"
          if (i==2) write(*,*) bsalNodes(j),j,"bsal wri"
          if (i==3) write(*,*) ltrlNodes(j),j,"ltrl wri"
       enddo
       
       write(*,*) " "
    enddo
       
  end subroutine get_positional_Nodes_with_diff_types_of_InsrtdNodes
  
  
  subroutine chk_allocatn_and_allocate_the_srfcNodes
    implicit none
    
    if (.not.allocated(apclNodes)) then
       allocate(apclNodes(1:numApclNodes))
    else
       deallocate(apclNodes) ; allocate(apclNodes(1:numApclNodes))
    endif
    
    if (.not.allocated(bsalNodes)) then
       allocate(bsalNodes(1:numBsalNodes))
    else
       deallocate(bsalNodes) ; allocate(bsalNodes(1:numBsalNodes))
    endif
    
    if (.not.allocated(ltrlNodes)) then
       allocate(ltrlNodes(1:numLtrlNodes))
    else
       deallocate(ltrlNodes) ; allocate(ltrlNodes(1:numLtrlNodes))
    endif
    
  end subroutine chk_allocatn_and_allocate_the_srfcNodes
  
  subroutine get_actvtn_fctrs_v0()
    implicit none
    integer :: j,jmax,k,kmax
    integer :: nodeNm,jmaxHlf
    integer :: lftNode,rghtNode
    integer :: sumOfAllNAEC
    
    sumofAllNAEC = NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl
    jmax         = 1+addedNCPair           ! this 2*1 reprsnt the added pairs
    
    do j = 1,(jmax-1)
       
       lftNode  = (2*Hlf_Ncell)+1 - (j-1)*2
       rghtNode = (2*(Hlf_Ncell+1))*2-1 - (j-1)*2
       
       activtnFctrApcl(lftNode)  = 1
       activtnFctrApcl(rghtNode) = 1
       
       write(*,*) activtnFctrApcl(lftNode),activtnFctrApcl(rghtNode),"actvtvFctr Rglr"
       
    enddo
    
    jmax        = 1+2*addedNCPair
    kmax        = NAEC_Apcl
    jmaxHlf     = jmax/2
    
    do j = 1,jmax 
       kmax = NAEC_Apcl
       do k = 1,kmax
          if (j.le.jmaxHlf)     nodeNm = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell-jmaxHlf+(j-1))*(sumOfAllNAEC)+k
          if (j == (jmaxHlf+1)) nodeNm = 2*(Hlf_Ncell+1)*2+(Hlf_Ncell*2)*(sumOfAllNAEC)+k
          if (j.gt.(jmaxHlf+1)) nodeNm = 2*(Hlf_Ncell+1)*2+((Hlf_Ncell*2)-(j-jmaxHlf-1))*(sumOfAllNAEC)+k 
          
          activtnFctrApcl(nodeNm) = 1
          write(*,*) activtnFctrApcl(nodeNm),"activtnFctr Insrtd"
       enddo   
    enddo
    
  end subroutine get_actvtn_fctrs_v0
  
  subroutine get_actvtn_fctrs()
    implicit none
    integer :: NumActvCell
    integer :: NNApclL,NNApclR           !NNApclL=Number of Nodes in Apical (Left/Right)
    integer :: InactvNodesL,InactvNodesR 
    integer :: i,j,lim1,lim2,lim3,lim4,lim5
    
    NumActvCell   = addedNCPair
    
    NNApclL       = (Hlf_Ncell)*(1+NAEC_Apcl) + 1
    NNApclR       = NNApclL
    InactvNodesL  = NNApclL - (NumActvCell*(1+NAEC_Apcl))
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
    
  end subroutine get_actvtn_fctrs
  
end module SRyp_variables


module curve_variables
  use sys_building_info
  use cell_info
  use end_info
  use non_changing_parameters
  use system_parameters
  
  implicit none
  
contains
  
  subroutine allocate_and_initialize_curve_variables
    implicit none

    open(unit=109,file='curve_vars.dat')
    
    call get_N_curve
    write(109,*) N_curve,"N_curve"
    
    allocate(curveX1Y1(1:N_curve,1:N_dmnsn))
    allocate(curveX2Y2(1:N_curve,1:N_dmnsn))
    allocate(P1(1:N_curve),P2(1:N_curve),dP(1:N_curve))
    allocate(TC(1:N_curve))
    allocate(R(1:N_curve))

    allocate(xc_crv(1:N_curve),yc_crv(1:N_curve))
    allocate(AnglS(1:N_curve),AnglE(1:N_curve))
    
    curveX1Y1 = -1.d30
    curveX2Y2 = -1.d30
    
    P1 = -1.d30 ; P2 = -1.d30 ; dP = -1.d30
    TC = -1.d30 
    R  = -1.d30

    xc_crv = -1.d30 ; yc_crv = -1.d30
    AnglS  = 500.0d0 ; AnglE = 500.0d0
    
    close(109)
    
  end subroutine allocate_and_initialize_curve_variables
  
  subroutine deallocate_curve_variables
    implicit none
    
    deallocate(curveX1Y1,curveX2Y2)
    deallocate(P1,P2,dP)
    deallocate(TC,R)
    
    deallocate(xc_crv,yc_crv)
    deallocate(AnglS,AnglE)

  end subroutine deallocate_curve_variables
  
  
  subroutine get_N_curve
    implicit none
    
    if (NI_incldAS_woPP == 0) call get_N_curve_methd1
    if (NI_incldAS_woPP == 1) call get_N_curve_methd2
    
  end subroutine get_N_curve
  
  subroutine get_N_curve_methd1
    implicit none
    integer :: i,j,jmax,jEnd
    integer :: spr_nm,typ_nm
    integer :: cnt,cnt_Pspr !Pspr=P  
    
    cnt      = 0
    jmax     = 3
    cnt_Pspr = 0
    
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
             
             spr_nm = (i-1)*jmax + j
             typ_nm = typ_spr(spr_nm)
             
             if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                cnt = cnt+1
             endif
             
          enddo
          
       elseif (stageNo==2 .or. stageNo==3) then
          write(*,*) "fl:info_module,sb:get_N_curve"
          stop
          
       elseif (stageNo==4) then
          
          do j = 1,jEnd
             
             spr_nm = (i-1)*jmax + j
             typ_nm = typ_spr(spr_nm)
             
             if (i.le.(ncl+ncr)) then
                
                if (typ_nm==2 .or. typ_nm==4 .or. typ_nm==5) then
                   cnt = cnt+1
                endif
                
             elseif (i.gt.(ncl+ncr) .and. i.le.(ncl+ncr+2)) then
                
                if (typ_nm .ne. 6) then
                   cnt = cnt+1
                endif
                
             elseif (i.gt.(ncl+ncr+2)) then
                
                if (typ_nm==4 .or. typ_nm==5) then
                   cnt = cnt+1
                endif
                
             endif
             
          enddo
          
       endif
       
    enddo
    
    N_curve = cnt
    write(*,*) N_curve,"N_curve"
    
  end subroutine get_N_curve_methd1
  
  subroutine get_N_curve_methd2
    implicit none
    integer :: cnt_apcl,cnt_bsal,cnt_ltrl,cnt_Sum
    integer :: i,j
    
    cnt_apcl=0 ; cnt_bsal=0 ; cnt_ltrl=0
    
    if (NAEC_Apcl .ne. 0) cnt_apcl = 1
    if (NAEC_Bsal .ne. 0) cnt_bsal = 1
    if (NAEC_Ltrl .ne. 0) cnt_ltrl = 1
    
    cnt_Sum = cnt_apcl+cnt_bsal+cnt_ltrl
    
    if (NI_incldAS_woPP==0) Hlf_Ncell=(N_cell)/2
    if (NI_incldAS_woPP==1) Hlf_Ncell=(N_cell-1)/2
    
    N_curve = 0
    
    if (CellsMeet == 0) N_curve = (Hlf_Ncell)*(cnt_Sum)*2 + 1 + 1 
    if (CellsMeet.gt.0) N_curve = (Hlf_Ncell)*(cnt_Sum)*2 + 1 
    
    write(*,*) NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl,"NAEC_val in N_curve-methd2"
    write(*,*) Hlf_Ncell,N_curve,"N_curve methd2"
    
  end subroutine get_N_curve_methd2
  
  
end module curve_variables

