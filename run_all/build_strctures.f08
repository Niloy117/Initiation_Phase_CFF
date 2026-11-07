module build_struct
  use calls_for_tests
  use redefining_system_module
  use PropTrnsferToNxtCell
  use switch_and_unswitch_models
  use nodeInsrtd_cntrlStates
  use Adding_cells
  use Manipulate_and_equilibrate
  use Control_State_Length_Match
  implicit none
  
  integer, allocatable :: thrshAnglInfo(:,:)
  real*8 , allocatable :: thrshAngl(:)
  real*8 , allocatable :: curr_thrshAngl(:)
  
  logical :: thrsh_lgcl
  integer :: inpt_rngeFlag
  integer :: readOrWriteForRunningCycle
  integer :: arbtry_NstepDcsn=0
  
contains
  
  subroutine build_struct_And_saveProp_S4
    implicit none
    integer :: i
    
    call get_blcks
    call updt_kspr_with_l0
    call updt_karea_with_A0
    !call print_the_system
    !stop
    
    call get_all_the_transfrms
    call print_the_system
    call get_all_moving_coordnte_variables!convrsn_fr_mving_sys.f08
    call get_all_gradient_variables !gradient_calcEnrgy.f08
    call get_all_grdPenF_variables !gradient_calcPenF.f08   
    !stop
    
    if(strctNo==2) then
       N_struct = 3
       call allocate_and_initialize_struct_variables
    endif
    
    call print_the_system
    
    !call test_05
    !call updt_kspr_with_l0
    !call updt_karea_with_A0
    !call saveA0l0_ofStrct(strctNo) !(stageNo,strctNo)
    
    call read_strctProps(strctNo)
    call Strct_to_actual(strctNo)
    
    if (strctNo==1) then
       call allocate_and_initialize_thrshld_vars
       call get_thrshold_angleInfo(coordntes_xy,thrshAnglInfo,thrshAngl)!Manipulating_prop
    endif
    
  end subroutine build_struct_And_saveProp_S4
  
  subroutine build_struct_And_saveProp_S1
    implicit none
    
    if (strctNo==1) then
       call get_blcks
       call updt_kspr_with_l0
       call updt_karea_with_A0
       !call print_the_system
       !stop
       
       call get_all_the_transfrms
       call print_the_system
       call get_all_moving_coordnte_variables!convrsn_fr_mvin
       call get_all_gradient_variables !gradient_calcEnrgy
       call get_all_grdPenF_variables !gradient_calcPenF
       
       call allocate_and_initialize_struct_variables
    endif
    
    if (stageType==1) then
       !call test_S1T1
       call test_S1T1_TargetVersion
       call saveA0l0_ofStrct(strctNo)
       !call read_strctProps(strctNo)
       !call Strct_to_actual(strctNo)
       
    elseif (stageType==2) then
       call test_S1T2
    endif
    
  end subroutine build_struct_And_saveProp_S1
  
  subroutine adjusting_routines_in_trnsfrm_strcts_S1(strct_INPT,strct_OTPT)
    implicit none
    integer, intent(in) :: strct_INPT,strct_OTPT
    
    integer :: LftNeigh=0,RghtNeigh=0
    logical :: lgcl_rdfn
    real*8  :: tol_Rdfn
    integer :: sprNo
    integer :: i,cnting
    real*8  :: Enrgy_frst!,Esf,Eaf,Egf,Ebf,Esindvd
    integer :: cntV
    integer :: printORread
    
    if (strct_INPT==1 .and. strct_OTPT==2) then ! 1
       
       call get_blcks
       call updt_kspr_with_l0
       call updt_karea_with_A0
       
       call get_all_the_transfrms
       call print_the_system
       call get_all_moving_coordnte_variables
       call get_all_gradient_variables
       call get_all_grdPenF_variables
       
       call allocate_and_initialize_struct_variables
       
       Enrgy_frst = Energy(node_xy,l0,A0)
       write(*,*) Enrgy_frst,"Enrgy_frst"
       
       grd_mv = 10.0d5 ; grdmv_xy = 10.0d5       
       
       call nodesl0A0_to_coordntes(node_xy,l0,A0,coordntes)
       coordntes_xy(1:N_mvCoordnte) = coordntes(1:N_mvCoordnte)
       
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
    elseif (strct_INPT==3 .and. strct_OTPT==4) then ! 2
       
       lgcl_rdfn=.False. ; tol_Rdfn=0.05d0 ; call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " ";write(*,*) "HURRAY, THE CELLS MEET 1st Time !!!!!!!";write(*,*) " "
          CellsMeet = 1
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " :( Cells did not meet yet, stop here to adjust (1st)"
          stop
       endif
       
       
       write(*,*) "Befre IT"
       call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       
       sprNo      = cntrlCell_Spring(1)
       write(*,*) sprNo,cntrlCell_Spring(1),"cntrlCell_Spring 1"
       
       dmlshDecsn = 1
       call demolish_spring_and_redefine_system(sprNo)
       
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       call switch_to_NI_model
       call print_trnsfrm_at_diffrntCellsMeet
       printORread=1 ; call printORread_NIorTN_trnsfrms_inFiles(printORread)
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (strct_INPT==4 .and. strct_OTPT==5) then ! 3
       continue
       
    elseif (strct_INPT==6 .and. strct_OTPT==7) then ! 4
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.59 !0.08d0
       write(*,*) "change this Rdfn to 0.05"
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 2nd Time !!!!!!!"
          write(*,*) " "
          CellsMeet = 2
             
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (2nd)"     
          stop
       endif
       
       call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       call switch_to_NI_model
       call print_trnsfrm_at_diffrntCellsMeet
       printORread=1 ; call printORread_NIorTN_trnsfrms_inFiles(printORread)
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (strct_INPT==7 .and. strct_OTPT==8) then ! 5
       continue
       
    elseif (strct_INPT==9 .and. strct_OTPT==10) then ! 6
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.70 !0.11d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 3rd Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=3
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (3rd)"
          stop
       endif
       
       call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       call switch_to_NI_model
       call print_trnsfrm_at_diffrntCellsMeet
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (strct_INPT==10 .and. strct_OTPT==11) then ! 7
       continue
       
    elseif (strct_INPT==12 .and. strct_OTPT==13) then ! 8
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.30d0 !0.18d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " ";write(*,*) "HURRAY, THE CELLS MEET 4th Time !!!!!!!";write(*,*) " "
          CellsMeet=4
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust (4th)"
          stop
       endif
       
       call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       call switch_to_NI_model
       call print_trnsfrm_at_diffrntCellsMeet
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (strct_INPT==13 .and. strct_OTPT==14) then ! 9
       continue
       
    elseif (strct_INPT==15 .and. strct_OTPT==16) then ! 10
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.09d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          write(*,*) " "
          write(*,*) "HURRAY, THE CELLS MEET 5th Time !!!!!!!"
          write(*,*) " "
          
          CellsMeet=5
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          write(*,*) " :( Cells did not meet yet, stop here to adjust"
          stop
       endif
       
       call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       !Frame = Frame+1
       !call switchto_NI_model_run_and_switchbackto_TN
       
       call switch_to_NI_model
       call print_trnsfrm_at_diffrntCellsMeet
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    endif
    
  end subroutine adjusting_routines_in_trnsfrm_strcts_S1
  
  
  subroutine adjusting_routines_for_NI_systems(CellsMeetV,PCSorICS)
    implicit none
    integer, intent(in) :: CellsMeetV
    integer, intent(in) :: PCSorICS

    real*8  :: tol_Rdfn
    logical :: lgcl_Rdfn
    integer :: LftNeigh=0,RghtNeigh=0
    
    write(*,*) CellsMeetV,CellsMeet,"CellsMeetV-CellsMeet"
    if (CellsMeetV .ne. CellsMeet) stop 'inconsistent CellsMeet' 
    
    
    if (PCSorICS == 2) then ! PCS=1 ; ICS=2 [this will be true for ICS]
       
       lgcl_Rdfn=.False. ; tol_Rdfn=0.07d0
       call get_decision_of_redefining_diffWay(tol_Rdfn,lgcl_Rdfn)
       
       if (lgcl_rdfn .eqv. .True.) then
          call taking_lftrghtOriginNeigh_downwards_diffWay(LftNeigh,RghtNeigh)
          CellsMeet = CellsMeet+1 ; write(*,*) CellsMeet,"CM pos1"
       else
          stop 'lgcl_rdfn is False'
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr_NI_system
       call Equilibrate_only_NI_model
       
    endif

  end subroutine adjusting_routines_for_NI_systems
  
  subroutine Property_Trnsfr_ToNextCell
    
    implicit none
    real*8  :: ks_I(1:N_spr) ,ks_F(1:N_spr) ,l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8  :: ka_I(1:N_cell),ka_F(1:N_cell),A0_I(1:N_cell),A0_F(1:N_cell)
    
    integer :: InitFirstPTCell,FirstPTCell,BotmCellPExclusn,ApprchingCell
    
    integer :: cellPTNum,sprPTNum
    integer :: prtn,N_itm
    integer :: i,j,jmax
    
    logical :: lgcl_rdfn
    real*8  :: tol_rdfn
    integer :: LftNeigh,RghtNeigh
    integer :: ProgCycl
    
    integer :: NCP
    real*8  :: hmMltplL(1:3)
    integer :: strctV,strctV1,strctV2
    integer :: strghtOrNI
    
    
    N_progCycl   = 3
    NI_saveOrNot = 1
    
    call allocate_ProgrsnStgProp
    call allocate_ProgrsnStgNIProp
    
    call ProgStgPropTrnsfr
    
    
    do ProgCycl=1,N_progCycl
       
       ProgCyclNo = ProgCycl
       
       !if (ProgCycl==1) then
       nodeInsrtnStrts=1
       !elseif (ProgCycl.ne.1) then
          !nodeInsrtnStrts=0
       !endif
       
       !if (ProgCycl==2) stop
       
       InitFirstPTCell  = 3
       FirstPTCell      = 3-ProgCycl+1
       BotmCellPExclusn = InitFirstPTCell-FirstPTCell+1
       ApprchingCell    = FirstPTCell+4
       
       call getPT_arrays(ProgCycl,FirstPTCell,BotmCellPExclusn,cellPT,sprPT,nodePT)
       call get_SmthCellSprsArray(ApprchingCell)
       
       prtn      = 1 ; prtnVal   = prtn
       strctVno1 = 1 ; strctVno2 = 2
       !call get_PT_InptOtpt(prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       call get_PT_InptOtpt_WO_Intr(prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)!Intrmdte
       call trnsfrm_btwn_strcts_PT(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       strctV = 1 ; strghtOrNI = 1
       call save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0_I,l0_I,ka_I,ks_I)
       
       NCP=5+(ProgCycl-1);hmMltplL(1)=0.0d0;hmMltplL(2)=0.0d0;hmMltplL(3)=-0.60d0
       
       strctV = 2 ; strctVno2 = 2
       call Stepwise_SprShorteningNC_and_Equilibrate_UntilMet(Exprmnt,Frame,NCP,hmMltplL)
       
       A0_F(1:N_cell) = A0(1:N_cell)     ; l0_F(1:N_spr) = l0(1:N_spr)
       ka_F(1:N_cell) = k_area(1:N_cell) ; ks_F(1:N_spr) = k_spr(1:N_spr)
       
       strctV=2 ; strghtOrNI=1
       call save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0_F,l0_F,ka_F,ks_F)
       
       call writeForChkDataProgCycl(ProgCycl,strctV,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       call sleep(1)
       
       if (lgcl_rdfn.eqv..True.) then
          
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*) "THE CELLS MEET 6th Time !!!!!!!"
          if (ProgCycl==2) write(*,*) "THE CELLS MEET 7th Time !!!!!!!"
          if (ProgCycl==3) write(*,*) "THE CELLS MEET 8th Time !!!!!!!"
          
          write(*,*) " "
          
          if (ProgCycl==1) CellsMeet=6
          if (ProgCycl==2) CellsMeet=7
          if (ProgCycl==3) CellsMeet=8
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*)":( Cells didnt meet,stop here for ProgC=1"
          if (ProgCycl==2) write(*,*)":( Cells didnt meet,stop here for ProgC=2"
          if (ProgCycl==3) write(*,*)":( Cells didnt meet,stop here for ProgC=3"
          
          stop
          
       endif
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       
       prtn=2 ; prtnVal = prtn
       strctVno1 = 3 ; strctVno2 = 4
       !call get_PT_InptOtpt(prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       call get_PT_InptOtpt_WO_Intr(Prtn,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       call trnsfrm_btwn_strcts_PT(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       strctV=3 ; strghtOrNI = 1
       call save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0_I,l0_I,ka_I,ks_I)
       strctV=4 ; strghtOrNI = 1
       call save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0_F,l0_F,ka_F,ks_F)
       
       
       strghtOrNI = 2
       call resave_NI_strctV3_ifNeeded(ProgCycl,strghtOrNI)
       
       call writeForChkDataProgCycl(ProgCycl,strctV,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
    enddo
    
    
  end subroutine Property_Trnsfr_ToNextCell
  
  
  subroutine writeForChkDataProgCycl(ProgCycl,strctV,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer,intent(in) :: ProgCycl,strctV
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    integer :: i,j,jmax
    integer :: N_itm
    
    open(unit=93,file='ProgStgDataChk.dat',position='append')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             write(93,*) ks_I(j),l0_I(j),(strctV-1),ks_F(j),l0_F(j),strctV
          elseif (i==2) then
             write(93,*) ka_I(j),A0_I(j),(strctV-1),ka_F(j),A0_F(j),strctV
          elseif (i==3) then
             write(93,*) CgXNode(j)
          endif
          
       enddo
       
       write(93,*) " "
       
    enddo
    
    close(93)
    
  end subroutine writeForChkDataProgCycl
  
  subroutine PrgrsnStgCyclicMovement
    implicit none
    integer :: ProgCycl,strctV
    real*8  :: ksV(1:N_spr),kaV(1:N_cell)
    real*8  :: l0V(1:N_spr),A0V(1:N_cell)
    real*8  :: CgXV(1:N_node)
    
    real*8  :: ks_I(1:N_spr) ,ks_F(1:N_spr) ,l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8  :: ka_I(1:N_cell),ka_F(1:N_cell),A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    integer :: LftNeigh,RghtNeigh
    
    integer :: strghtOrNI
    
    N_progCycl   = 3
    !NI_saveOrNot = 1
    
    do ProgCycl = 1,N_progCycl
       
       if (ProgCycl==1) nodeInsrtnStrts=1
       
       strctV=1 ; strghtOrNI = 1 
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgXV)
       ks_I(1:N_spr)     = ksV(1:N_spr)   ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)  ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgXV(1:N_node) ; CgYNode(1:N_node) = 0.0d0
       
       strctV=2 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgXV)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       call trnsfrm_btwn_strcts_PT(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       
       if (lgcl_rdfn.eqv..True.) then
          
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*) "THE CELLS MEET 6th Time !!!!!!!"
          if (ProgCycl==2) write(*,*) "THE CELLS MEET 7th Time !!!!!!!"
          if (ProgCycl==3) write(*,*) "THE CELLS MEET 8th Time !!!!!!!"
          
          write(*,*) " "
          
          if (ProgCycl==1) CellsMeet=6
          if (ProgCycl==2) CellsMeet=7
          if (ProgCycl==3) CellsMeet=8
          
       elseif (lgcl_rdfn.eqv..False.) then
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*)" :( Cells didnt meet,stop here for ProgC=1"
          if (ProgCycl==2) write(*,*)" :( Cells didnt meet,stop here for ProgC=2"
          if (ProgCycl==3) write(*,*)" :( Cells didnt meet,stop here for ProgC=3"
          
          stop
          
       endif
       
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       
       strctV=3 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgXV)
       ks_I(1:N_spr)     = ksV(1:N_spr)   ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)  ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgXV(1:N_node) ; CgXNode(1:N_node) = 0.0d0
       
       
       strctV=4 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgXV)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       call trnsfrm_btwn_strcts_PT(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       
    enddo
    
    
  end subroutine PrgrsnStgCyclicMovement
  
  
  subroutine PrgrsnStgCyclicMovement_Separately
    implicit none
    integer :: ProgCycl
    
    real*8, allocatable :: ks_I(:),ks_F(:),l0_I(:),l0_F(:)
    real*8, allocatable :: ka_I(:),ka_F(:),A0_I(:),A0_F(:)
    real*8, allocatable :: CgX_V(:)
    
    real*8, allocatable :: ksV(:),l0V(:)
    real*8, allocatable :: kaV(:),A0V(:)
    
    real*8, allocatable :: ks_SV(:),l0_SV(:)
    real*8, allocatable :: ka_SV(:),A0_SV(:)
    
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn,lgcl_meet
    integer :: LftNeigh,RghtNeigh
    
    integer :: strctV,strghtOrNI   
    real*8  :: E
    integer :: i
    
    call switch_to_NI_model
    call allocate_and_initialize_NIsavedVars
    call get_NIsaved_Properties
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    N_progCycl   = 3
    
    do ProgCycl = 1,N_progCycl
       
       if (ProgCycl==1) nodeInsrtnStrts=1
       !if (ProgCycl==2) stop
       
       allocate(ks_I(1:N_spr) ,l0_I(1:N_spr))
       allocate(ks_F(1:N_spr) ,l0_F(1:N_spr))
       allocate(ka_I(1:N_cell),A0_I(1:N_cell))
       allocate(ka_F(1:N_cell),A0_F(1:N_cell))
       allocate(CgX_V(1:N_node))
       
       allocate(ksV(1:N_spr) ,l0V(1:N_spr))
       allocate(kaV(1:N_cell),A0V(1:N_cell))
       
       strctV=1 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_I(1:N_spr)     = ksV(1:N_spr)    ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)   ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgX_V(1:N_node) ; CgYNode(1:N_node) = 0.0d0
       
       strctV=2 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       call trnsfrm_btwn_strct_PT_woNI(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       deallocate(ks_I,l0_I)
       deallocate(ks_F,l0_F)
       deallocate(ka_I,A0_I)
       deallocate(ka_F,A0_F)
       deallocate(CgX_V)
       
       deallocate(ksV,l0V)
       deallocate(kaV,A0V)
       
       call switch_to_NIsaved_model
       
       allocate(ks_I(1:N_spr) ,l0_I(1:N_spr))
       allocate(ks_F(1:N_spr) ,l0_F(1:N_spr))
       allocate(ka_I(1:N_cell),A0_I(1:N_cell))
       allocate(ka_F(1:N_cell),A0_F(1:N_cell))
       allocate(CgX_V(1:N_node))
       
       allocate(ksV(1:N_spr) ,l0V(1:N_spr))
       allocate(kaV(1:N_cell),A0V(1:N_cell))
       
       strctV=1 ; strghtOrNI = 2
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_I(1:N_spr)     = ksV(1:N_spr)    ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)   ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgX_V(1:N_node) ; CgYNode(1:N_node) = 0.0d0
       
       strctV=2 ; strghtOrNI = 2
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       !if (ProgCycl==2) then
          
        !  call Equilibrate_system
         ! call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          !Frame_NI = Frame_NI+1
          
       !endif
       
       call trnsfrm_btwn_strct_PT_woNI(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       call decision_of_meeting_NI(lgcl_meet)
       
       if (lgcl_meet.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
          call get_NIsaved_Properties
          
       elseif (lgcl_meet.eqv..False.) then
          write(*,*) "cells not met for NI when ProgCycl =",ProgCycl
          stop
       endif
       
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
       deallocate(ks_I,l0_I)
       deallocate(ks_F,l0_F)
       deallocate(ka_I,A0_I)
       deallocate(ka_F,A0_F)
       deallocate(CgX_V)
       
       deallocate(ksV,l0V)
       deallocate(kaV,A0V)
       
       
       lgcl_rdfn = .False.
       tol_Rdfn  = 0.05d0
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       if (lgcl_rdfn.eqv..True.) then
          
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*) "THE CELLS MEET 6th Time !!!!!!!"
          if (ProgCycl==2) write(*,*) "THE CELLS MEET 7th Time !!!!!!!"
          if (ProgCycl==3) write(*,*) "THE CELLS MEET 8th Time !!!!!!!"
          
          write(*,*) " "
          
          if (ProgCycl==1) CellsMeet=6
          if (ProgCycl==2) CellsMeet=7
          if (ProgCycl==3) CellsMeet=8
          
       elseif (lgcl_rdfn.eqv..False.) then
          
          write(*,*) " "
          
          if (ProgCycl==1) write(*,*)" :( Cells didnt meet,stop here for ProgC=1"
          if (ProgCycl==2) write(*,*)" :( Cells didnt meet,stop here for ProgC=2"
          if (ProgCycl==3) write(*,*)" :( Cells didnt meet,stop here for ProgC=3"
          
          stop
          
       endif
       
       if (lgcl_rdfn.eqv..True.) then
          call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
       elseif (lgcl_rdfn.eqv..True.) then
          write(*,*) "lgcl is not true,stop"
          stop
       endif
       
       dmlshDecsn = 0
       call redefine_system_wo_demolishSpr
       call Equilibrate_system
       
       
       allocate(ks_I(1:N_spr) ,l0_I(1:N_spr))
       allocate(ks_F(1:N_spr) ,l0_F(1:N_spr))
       allocate(ka_I(1:N_cell),A0_I(1:N_cell))
       allocate(ka_F(1:N_cell),A0_F(1:N_cell))
       allocate(CgX_V(1:N_node))
       
       allocate(ksV(1:N_spr) ,l0V(1:N_spr))
       allocate(kaV(1:N_cell),A0V(1:N_cell))
       
       strctV=3 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_I(1:N_spr)     = ksV(1:N_spr)    ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)   ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgX_V(1:N_node) ; CgYNode(1:N_node) = 0.0d0
       
       
       strctV=4 ; strghtOrNI = 1
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       call trnsfrm_btwn_strct_PT_woNI(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       
       deallocate(ks_I,l0_I)
       deallocate(ks_F,l0_F)
       deallocate(ka_I,A0_I)
       deallocate(ka_F,A0_F)
       deallocate(CgX_V)
       
       deallocate(ksV,l0V)
       deallocate(kaV,A0V)
       
       !call switch_to_NI_model
       call switch_to_NIsaved_model
       
       allocate(ks_I(1:N_spr) ,l0_I(1:N_spr))
       allocate(ks_F(1:N_spr) ,l0_F(1:N_spr))
       allocate(ka_I(1:N_cell),A0_I(1:N_cell))
       allocate(ka_F(1:N_cell),A0_F(1:N_cell))
       allocate(CgX_V(1:N_node))
       
       allocate(ksV(1:N_spr) ,l0V(1:N_spr))
       allocate(kaV(1:N_cell),A0V(1:N_cell))
       
       strctV=3 ; strghtOrNI = 2
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_I(1:N_spr)     = ksV(1:N_spr)    ; l0_I(1:N_spr)     = l0V(1:N_spr)
       ka_I(1:N_cell)    = kaV(1:N_cell)   ; A0_I(1:N_cell)    = A0V(1:N_cell)
       CgXNode(1:N_node) = CgX_V(1:N_node) ; CgYNode(1:N_node) = 0.0d0
       
       
       strctV=4 ; strghtOrNI = 2
       call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgX_V)
       ks_F(1:N_spr)  = ksV(1:N_spr)  ; l0_F(1:N_spr)  = l0V(1:N_spr)
       ka_F(1:N_cell) = kaV(1:N_cell) ; A0_F(1:N_cell) = A0V(1:N_cell)
       
       !open(unit=81,file='nodeVal.dat')
       
       !do i = 1,N_node
        !  write(81,*) node_xy(i,1:2),i,"nodes"
       !enddo
       
       !close(81)
       
       !call get_Nlist
       !E = Energy(node_xy,l0,A0)
       !write(*,*) E,"E"
       !call get_gradient(node_xy,l0,A0,gradient)
       
       
       call trnsfrm_LittleBitReverseStrct3to4(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       stop
       call trnsfrm_btwn_strct_PT_woNI(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
       
       !if (ProgCycl==1) then
        !  call Equilibrate_system
         !call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          !Frame_NI = Frame_NI+1
       !endif
       
       call get_NIsaved_Properties
       
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
       !k_spr(1:N_spr)   = ks_TN(1:N_spr) !here ks_TN stores ks_F values
       !l0(1:N_spr)      = l0_TN(1:N_spr)
       !l(1:N_spr)
       
       !k_area(1:N_cell) = ka_TN(1:N_cell)
       !A0(1:N_cell)     = A0_TN(1:N_cell)
       
       deallocate(ks_I,l0_I)
       deallocate(ks_F,l0_F)
       deallocate(ka_I,A0_I)
       deallocate(ka_F,A0_F)
       deallocate(CgX_V)
       
       deallocate(ksV,l0V)
       deallocate(kaV,A0V)

       if (ProgCycl==1) stop
       
    enddo
    
    
  end subroutine PrgrsnStgCyclicMovement_Separately
  
  
  subroutine RunCyclicProcessForAddedCellNI
    implicit none
    integer :: FirstCellL,FirstCellR
    integer :: LastCellL ,LastCellR
    
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7
    integer :: N_lims
    integer :: sprNm,sprNmNAC
    integer :: cellCnt,cyclNm
    
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_cell)
    real*8  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell)
    real*8  :: ks_inpt(1:N_spr) ,ks_otpt(1:N_spr)
    real*8  :: l0_inpt(1:N_spr) ,l0_otpt(1:N_spr)
    
    integer, allocatable :: limVal(:)
    integer :: LftNeigh,RghtNeigh
    real*8  :: E
    real*8  :: tol_Rdfn
    logical :: lgcl_Rdfn
    integer :: weightType
    integer :: sprNum,areaNum,nodeNum
    integer :: callLoctn,strctToRead
    
    write(*,*) "ENTERED HERE in RUN"
    
    N_addedCycl = 4 ; N_CS = 4 ! check the N_addedCellL in AddedCell_model.f08 script
    readOrWriteForRunningCycle = 2
    arbtry_NstepDcsn=1 ! Very important, change it to 0, if needed 11 frame per half Cycle
    
    
    if (readOrWriteForRunningCycle == 1) then
       
       open(unit=192,file='writeVrsnForRunningCycl.dat')
       
       do cyclNm = 1,N_addedCycl
          
          FirstCellL = N_addedCellL - (cyclNm-1) 
          FirstCellR = Hlf_Ncell + N_addedCellR - (cyclNm-1)
          
          LastCellL  = FirstCellL + (Hlf_NcellNAC-2)
          LastCellR  = FirstCellR + (Hlf_NcellNAC-2)
          
          lim1 = FirstCellL-1
          lim2 = LastCellL
          lim3 = Hlf_Ncell
          lim4 = FirstCellR-1
          lim5 = LastCellR
          lim6 = N_cell-1
          lim7 = N_cell
          
          N_lims = 7
          allocate(limVal(1:N_lims))
          
          do i = 1,N_lims
             if (i==1) limVal(i) = lim1
             if (i==2) limVal(i) = lim2
             if (i==3) limVal(i) = lim3
             if (i==4) limVal(i) = lim4
             if (i==5) limVal(i) = lim5
             if (i==6) limVal(i) = lim6
             if (i==7) limVal(i) = lim7
          enddo
          
          write(192,*) limVal(1:7),"lims"
          
          hlfCycl = 1
          call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          call get_Props_AtConnectingCycls(cyclNm,hlfCycl,limVal,N_lims)
          !call get_InOutsForAC_Cycl(cyclNm,hlfCycl,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          call trnsfrm_btwn_strcts_AddedCellNI(ks_inpt,ks_otpt,l0_inpt,l0_otpt,ka_inpt,ka_otpt,A0_inpt,A0_otpt)
          
          tol_Rdfn = 0.10d0 ; lgcl_Rdfn = .False.
          call get_decision_of_redefining_AddedCellNI(tol_Rdfn,lgcl_Rdfn)
          
          write(*,*) lgcl_Rdfn,"Rdfn"
          
          if (lgcl_Rdfn.eqv..True.) then
             call taking_lftrghtOriginNeigh_downwards_AC_and_updtSys(LftNeigh,RghtNeigh)
          elseif (lgcl_Rdfn.eqv..False.) then
             write(*,*) "lgcl is not true,stop"
             
             call get_gradient(node_xy,l0,A0,gradient)
             call Find_Analytical_and_Numerical_Mismatch
             
             stop
          endif
          
          hlfCycl = 2
          call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          call get_Props_AtConnectingCycls(cyclNm,hlfCycl,limVal,N_lims)
          !call get_InOutsForAC_Cycl(cyclNm,hlfCycl,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          call trnsfrm_btwn_strcts_AddedCellNI(ks_inpt,ks_otpt,l0_inpt,l0_otpt,ka_inpt,ka_otpt,A0_inpt,A0_otpt)
          
          deallocate(limVal)
          
          !if(cyclNm==10) exit
       enddo
       
       do weightType = 1,3
          call write_ConnectingCyclProps(weightType)
          call get_IFProps_AtConnectingCycls(weightType)
       enddo
       
       close(192)
       
    elseif (readOrWriteForRunningCycle == 2) then
       
       open(unit=193,file='readVrsnForRunningCycl.dat')
       
       weightType = 3
       
       allocate(ka_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
       allocate(A0_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
       allocate(ks_cycl(1:N_addedCycl,1:N_CS,1:N_spr))
       allocate(l0_cycl(1:N_addedCycl,1:N_CS,1:N_spr))
       
       call read_Props_AtConnectingCycls(weightType)
       call writeCyclPropsAftRead
       
       write(*,*) "ENTERED HERE ALSO RUN 2"
       
       do cyclNm = 1,N_addedCycl
          
          hlfCycl = 1
          call get_InOutsForAC_Cycl(cyclNm,hlfCycl,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          
          if (cyclNm==1) continue
          call trnsfrm_btwn_strcts_AddedCellNI(ks_inpt,ks_otpt,l0_inpt,l0_otpt,ka_inpt,ka_otpt,A0_inpt,A0_otpt)
          write(*,*) "stopiing here"
          stop
          
          tol_Rdfn = 0.10d0 ; lgcl_Rdfn = .False.
          call get_decision_of_redefining_AddedCellNI(tol_Rdfn,lgcl_Rdfn)
          
          !** WILL RUN THE FORCE PERTURBATION ROUTINE HERE **!
          
          !if (prtrbtnRunLoctn == 2) then
          CyclVal = cyclNm
          callLoctn = 2 ; strctToRead = 2
          call force_prtrbtn_applied_at_boundry_of_CS2(callLoctn,strctToRead,ExprmntAdded,FrameAdded)
          call Extrapolate_Adjust_strctPropWithCond_in_MAE_woSwitchBack(ExprmntAdded,FrameAdded)
          stop
          !endif
          
          !** RER
          
          write(*,*) lgcl_Rdfn,"Rdfn"
          if (lgcl_Rdfn.eqv..True.) then
             call taking_lftrghtOriginNeigh_downwards_AC_and_updtSys(LftNeigh,RghtNeigh)
          elseif (lgcl_Rdfn.eqv..False.) then
             write(*,*) "lgcl is not true,stop"
             
             call get_gradient(node_xy,l0,A0,gradient)
             call Find_Analytical_and_Numerical_Mismatch
             stop
          endif
          
          hlfCycl = 2
          call get_InOutsForAC_Cycl(cyclNm,hlfCycl,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
          call trnsfrm_btwn_strcts_AddedCellNI(ks_inpt,ks_otpt,l0_inpt,l0_otpt,ka_inpt,ka_otpt,A0_inpt,A0_otpt)
          
       enddo
       
       close(193)
       
    endif
    
  end subroutine RunCyclicProcessForAddedCellNI
  
  
  subroutine productionRunForIntiatnPhase(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    real*8  :: ks_I(1:N_spr) , ks_F(1:N_spr)
    real*8  :: l0_I(1:N_spr) , l0_F(1:N_spr)
    real*8  :: ka_I(1:N_cell), ka_F(1:N_cell)
    real*8  :: A0_I(1:N_cell), A0_F(1:N_cell)
    real*8  :: CgX_I(1:N_node),CgX_F(1:N_node)
    real*8  :: CgY_I(1:N_node),CgY_F(1:N_node)
    integer :: PCSfl,PCSval
    integer :: ICSfl,ICSval
    integer :: i,j,N_step
    real*8  :: EnrgyVal,EsV,EaV,EgV,EbV
    
    arbtry_NstepDcsn      = 1
    if (arbtry_NstepDcsn == 0) N_step = 10 
    if (arbtry_NstepDcsn == 1) N_step = 6
    
    if (PCSorICS==1) then
       call read_the_prps_of_High_ks_tst_woNodeTYP(CMv,PCSorICS)
       call nodes_to_coordntes(node_xy,coordntes_xy)
       call turnOffbndryForceToChk
       
    elseif (PCSorICS==2) then
       FramePrdctnInit  = (N_step+1) + (CMv-2)*((N_step+1)*2)
       readNodeFilesTyp = 1 ; call read_config_and_start_simlnFrm_there(Exprmnt_Init,FramePrdctnInit)
       readNodeFilesTyp = 0
       call nodes_to_coordntes(node_xy,coordntes_xy)
       
       FramePrdctnInit = FramePrdctnInit+1 
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_Init,FramePrdctnInit)
       FramePrdctnInit = FramePrdctnInit !; stop 'aft pull sys'
    endif
    
    write(*,*) CMv,PCSorICS,"pp~"
    
    if (PCSorICS==1) then
       PCSfl=CMv ; PCSval=1 ; ICSfl=CMv+1 ; ICSval=2
       call read_the_prps_of_High_ks_tst_PCSorICS(PCSfl,PCSval,ks_I,l0_I,ka_I,A0_I,CgX_I,CgY_I)
       call read_the_prps_of_High_ks_tst_PCSorICS(ICSfl,ICSval,ks_F,l0_F,ka_F,A0_F,CgX_F,CgY_F)
    elseif (PCSorICS==2) then
       ICSfl=CMv ; ICSval=2 ; PCSfl=CMv   ; PCSval=1 
       call read_the_prps_of_High_ks_tst_PCSorICS(ICSfl,ICSval,ks_I,l0_I,ka_I,A0_I,CgX_I,CgY_I)
       call read_the_prps_of_High_ks_tst_PCSorICS(PCSfl,PCSval,ks_F,l0_F,ka_F,A0_F,CgX_F,CgY_F)
    endif
    
    call trnsfrm_btwn_strcts_varying_kskal0A0CgXCgY(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F,&
         CgX_I,CgX_F,CgY_I,CgY_F)
    
  end subroutine productionRunForIntiatnPhase
  
  subroutine get_Props_AtConnectingCycls(cyclNm,hlfCyclNm,limVal,N_lims)
    implicit none
    integer, intent(in) :: cyclNm,hlfCyclNm,N_lims
    integer, intent(in) :: limVal(1:N_lims)
    
    integer :: FirstCellL,FirstCellR
    integer :: LastCellL ,LastCellR
    integer :: i,j,jmax
    
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_cell),A0_inpt(1:N_cell),A0_otpt(1:N_cell)
    real*8  :: ks_inpt(1:N_spr) ,ks_otpt(1:N_spr) ,l0_inpt(1:N_spr) ,l0_otpt(1:N_spr)
    
    open(unit=292,file='limsValInConnCycl.dat')
    
    if (allctnCycl_cnt==0) then
       allocate(ka_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
       allocate(A0_cycl(1:N_addedCycl,1:N_CS,1:N_cell))
       allocate(ks_cycl(1:N_addedCycl,1:N_CS,1:N_spr))
       allocate(l0_cycl(1:N_addedCycl,1:N_CS,1:N_spr))

       allctnCycl_cnt = 1
    endif
    
    
    write(292,*) limVal(1:7),"lims"
    write(292,*) N_lims,"N_lims"
    
    if (hlfCyclNm == 1) then
       call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
       
       ka_cycl(cyclNm,1,1:N_cell) = ka_inpt(1:N_cell)
       ka_cycl(cyclNm,2,1:N_cell) = ka_otpt(1:N_cell)
       A0_cycl(cyclNm,1,1:N_cell) = A0_inpt(1:N_cell)
       A0_cycl(cyclNm,2,1:N_cell) = A0_otpt(1:N_cell)
       ks_cycl(cyclNm,1,1:N_spr)  = ks_inpt(1:N_spr)
       ks_cycl(cyclNm,2,1:N_spr)  = ks_otpt(1:N_spr)
       l0_cycl(cyclNm,1,1:N_spr)  = l0_inpt(1:N_spr)
       l0_cycl(cyclNm,2,1:N_spr)  = l0_otpt(1:N_spr)
       
    elseif (hlfCycl == 2) then
       call get_InOutsForAC(limVal,N_lims,ka_inpt,ka_otpt,A0_inpt,A0_otpt,ks_inpt,ks_otpt,l0_inpt,l0_otpt)
       
       ka_cycl(cyclNm,3,1:N_cell) = ka_inpt(1:N_cell)
       ka_cycl(cyclNm,4,1:N_cell) = ka_otpt(1:N_cell)
       A0_cycl(cyclNm,3,1:N_cell) = A0_inpt(1:N_cell)
       A0_cycl(cyclNm,4,1:N_cell) = A0_otpt(1:N_cell)
       ks_cycl(cyclNm,3,1:N_spr)  = ks_inpt(1:N_spr)
       ks_cycl(cyclNm,4,1:N_spr)  = ks_otpt(1:N_spr)
       l0_cycl(cyclNm,3,1:N_spr)  = l0_inpt(1:N_spr)
       l0_cycl(cyclNm,4,1:N_spr)  = l0_otpt(1:N_spr)
       
    endif
    
    close(292)
    
  end subroutine get_Props_AtConnectingCycls


  subroutine write_ConnectingCyclProps(weightType)
    implicit none
    integer, intent(in) :: weightType
    
    character(len=100)  :: flnmks,flnml0,flnmka,flnmA0
    character(len=100)  :: flnmbr1,flnmbr2
    character(len=100)  :: fullflnm1,fullflnm2

    integer :: i,j,jmax
    integer :: cyclNm,csNo
    
    if (weightType==1) then
       flnmks='ksWT1Cycl'
       flnml0='l0WT1Cycl'
       flnmka='kaWT1Cycl'
       flnmA0='A0WT1Cycl'
       
    elseif (weightType==2) then
       flnmks='ksWT2Cycl'
       flnml0='l0WT2Cycl'
       flnmka='kaWT2Cycl'
       flnmA0='A0WT2Cycl'
       
    elseif (weightType==3) then
       flnmks='ksWT3Cycl'
       flnml0='l0WT3Cycl'
       flnmka='kaWT3Cycl'
       flnmA0='A0WT3Cycl'
    endif
    
    do cyclNm = 1,N_addedCycl
       
       if (cyclNm.le.9) write(flnmbr1,'(i1.1)') cyclNm
       if ((cyclNm.gt.9) .and. (cyclNm.le.99)) write(flnmbr1,'(i2.2)') cyclNm
       
       do csNo = 1,N_CS
          
          write(flnmbr2,'(a,i1.1,a)') 'CS',csNo,'.dat'
          
          do i = 1,2
             
             if (i==1) then
                
                fullflnm1=trim(adjustl(flnmks))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                fullflnm2=trim(adjustl(flnml0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                
                open(unit=46,file=trim(adjustl(fullflnm1)))
                open(unit=47,file=trim(adjustl(fullflnm2)))
                
                jmax = N_spr
                
             elseif (i==2) then
                
                fullflnm1=trim(adjustl(flnmka))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                fullflnm2=trim(adjustl(flnmA0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                
                open(unit=48,file=trim(adjustl(fullflnm1)))
                open(unit=49,file=trim(adjustl(fullflnm2)))
                
                jmax = N_cell
             endif
             
             do j = 1,jmax
                
                if (i==1) then
                   write(46,*) ks_cycl(cyclNm,csNo,j),j
                   write(47,*) l0_cycl(cyclNm,csNo,j),j
                   
                elseif (i==2) then
                   write(48,*) ka_cycl(cyclNm,csNo,j),j
                   write(49,*) A0_cycl(cyclNm,csNo,j),j
                endif
                
             enddo
             
             if (i==1) then
                close(46)
                close(47)
             elseif (i==2) then
                close(48)
                close(49)
             endif
             
          enddo
          
       enddo

       !if (cyclNm==10) exit
    enddo
    
    
  end subroutine write_ConnectingCyclProps
  
  subroutine get_IFProps_AtConnectingCycls(weightType)
    implicit none
    integer, intent(in) :: weightType
    
    real*8 :: kaF1(1:N_cell),kaI2(1:N_cell),kaIF(1:N_cell) !IF=InterFace
    real*8 :: A0F1(1:N_cell),A0I2(1:N_cell),A0IF(1:N_cell)
    real*8 :: ksF1(1:N_spr) ,ksI2(1:N_spr) ,ksIF(1:N_spr)
    real*8 :: l0F1(1:N_spr) ,l0I2(1:N_spr) ,l0IF(1:n_spr)
    
    integer :: csFnsh,csStrt
    integer :: cyclBfr,cyclAft
    real*8  :: weight
    integer :: i,j,jmax
    integer :: IF_no,cyclNm,csNo
    
    character(len=100) :: flnmks,flnml0,flnmka,flnmA0
    character(len=100) :: flnmbr1,flnmbr2,flnmbr3,flnmbr4
    character(len=100) :: fullflnm1,fullflnm2,fullflnm3,fullflnm4
    
    
    if (weightType==1) weight = 0.00d0
    if (weightType==2) weight = 0.50d0
    if (weightType==3) weight = 1.00d0
    
    if (weightType==1) then
       flnmks='ksWT1Cycl'
       flnml0='l0WT1Cycl'
       flnmka='kaWT1Cycl'
       flnmA0='A0WT1Cycl'
       
    elseif (weightType==2) then
       flnmks='ksWT2Cycl'
       flnml0='l0WT2Cycl'
       flnmka='kaWT2Cycl'
       flnmA0='A0WT2Cycl'
       
    elseif (weightType==3) then
       flnmks='ksWT3Cycl'
       flnml0='l0WT3Cycl'
       flnmka='kaWT3Cycl'
       flnmA0='A0WT3Cycl'
    endif

    csFnsh = 4 ; csStrt  = 1 
    
    do IF_no = 1,(N_addedCycl-1)
       
       cyclBfr = IF_no ; cyclAft = IF_no+1
       
       kaF1(1:N_cell) = ka_cycl(cyclBfr,csFnsh,1:N_cell)
       kaI2(1:N_cell) = ka_cycl(cyclAft,csStrt,1:N_cell)
       
       A0F1(1:N_cell) = A0_cycl(cyclBfr,csFnsh,1:N_cell)
       A0I2(1:N_cell) = A0_cycl(cyclAft,csStrt,1:N_cell)
       
       ksF1(1:N_spr)  = ks_cycl(cyclBfr,csFnsh,1:N_spr)
       ksI2(1:N_spr)  = ks_cycl(cyclAft,csStrt,1:N_spr)
       
       l0F1(1:N_spr)  = l0_cycl(cyclBfr,csFnsh,1:N_spr)
       l0I2(1:N_spr)  = l0_cycl(cyclAft,csStrt,1:N_spr)
       
       do i = 1,2
          
          if (i==1) jmax = N_cell
          if (i==2) jmax = N_spr
          
          do j = 1,jmax
             
             if (i==1) then
                kaIF(j) = weight*kaF1(j) + (1.0d0-weight)*kaI2(j)
                A0IF(j) = weight*A0F1(j) + (1.0d0-weight)*A0I2(j)
                
             elseif (i==2) then
                ksIF(j) = weight*ksF1(j) + (1.0d0-weight)*ksI2(j) 
                l0IF(j) = weight*l0F1(j) + (1.0d0-weight)*l0I2(j) 
             endif
             
          enddo
       
       enddo
       
       !ka_cycl(cyclBfr,csFnsh,1:N_cell) = kaIF(1:N_cell)
       !ka_cycl(cyclAft,csStrt,1:N_cell) = kaIF(1:N_cell)
       !A0_cycl(cyclBfr,csFnsh,1:N_cell) = A0IF(1:N_cell)
       !A0_cycl(cyclAft,csStrt,1:N_cell) = A0IF(1:N_cell)
       
       !ks_cycl(cyclBfr,csFnsh,1:N_spr)  = ksIF(1:N_spr)
       !ks_cycl(cyclAft,csStrt,1:N_spr)  = ksIF(1:N_spr)
       !l0_cycl(cyclBfr,csFnsh,1:N_spr)  = l0IF(1:N_spr)
       !l0_cycl(cyclAft,csStrt,1:N_spr)  = l0IF(1:N_spr)

       if (cyclBfr.le.9) then
          write(flnmbr1,'(i1.1)') cyclBfr
          
       elseif ((cyclBfr.gt.9) .and. (cyclBfr.le.99)) then
          write(flnmbr1,'(i2.2)') cyclBfr
       endif
       
       if (csFnsh.le.9) then
          write(flnmbr2,'(a,i1.1,a)') 'CS',csFnsh,'.dat'
          
       elseif ((csFnsh.gt.9) .and. (csFnsh.le.99)) then  
          write(flnmbr2,'(a,i2.2,a)') 'CS',csFnsh,'.dat'
       endif
       
       if (cyclAft.le.9) then
          write(flnmbr3,'(i1.1)') cyclAft
          
       elseif ((cyclAft.gt.9) .and. (cyclAft.le.99)) then
          write(flnmbr3,'(i2.2)') cyclAft
       endif

       if (csStrt.le.9) then
          write(flnmbr4,'(a,i1.1,a)') 'CS',csStrt,'.dat'
          
       elseif ((csStrt.gt.9) .and. (csStrt.le.99)) then
          write(flnmbr4,'(a,i2.2,a)') 'CS',csStrt,'.dat'
       endif
          
          
       do i = 1,2
          
          if (i==1) then
             
             fullflnm1=trim(adjustl(flnmks))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
             fullflnm2=trim(adjustl(flnml0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
             fullflnm3=trim(adjustl(flnmks))//trim(adjustl(flnmbr3))//trim(adjustl(flnmbr4))
             fullflnm4=trim(adjustl(flnml0))//trim(adjustl(flnmbr3))//trim(adjustl(flnmbr4))
             
             open(unit=46,file=trim(adjustl(fullflnm1)))
             open(unit=47,file=trim(adjustl(fullflnm2)))
             open(unit=50,file=trim(adjustl(fullflnm3)))
             open(unit=51,file=trim(adjustl(fullflnm4)))
             
             jmax = N_spr
             
          elseif (i==2) then
             
             fullflnm1=trim(adjustl(flnmka))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
             fullflnm2=trim(adjustl(flnmA0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
             fullflnm3=trim(adjustl(flnmka))//trim(adjustl(flnmbr3))//trim(adjustl(flnmbr4))
             fullflnm4=trim(adjustl(flnmA0))//trim(adjustl(flnmbr3))//trim(adjustl(flnmbr4))
             
             open(unit=48,file=trim(adjustl(fullflnm1)))
             open(unit=49,file=trim(adjustl(fullflnm2)))
             open(unit=52,file=trim(adjustl(fullflnm3)))
             open(unit=53,file=trim(adjustl(fullflnm4)))
             
             jmax = N_cell
          endif
          
          do j = 1,jmax
             if (i==1) then
                write(46,*) ksIF(j),j
                write(47,*) l0IF(j),j
                write(50,*) ksIF(j),j
                write(51,*) l0IF(j),j
                
             elseif (i==2) then
                write(48,*) kaIF(j),j
                write(49,*) A0IF(j),j
                write(52,*) kaIF(j),j
                write(53,*) A0IF(j),j
             endif
             
          enddo
          
          if (i==1) then
             close(46)
             close(47)
             close(50)
             close(51)
          elseif (i==2) then
             close(48)
             close(49)
             close(52)
             close(53)
          endif
          
       enddo
       
       if (IF_no==9) exit
    enddo
    
    
  end subroutine get_IFProps_AtConnectingCycls
  
  
  
  subroutine read_Props_AtConnectingCycls(weightType)
    implicit none
    integer, intent(in) :: weightType
    
    real*8 :: kaF1(1:N_cell),kaI2(1:N_cell),kaIF(1:N_cell) !IF=InterFace
    real*8 :: A0F1(1:N_cell),A0I2(1:N_cell),A0IF(1:N_cell)
    real*8 :: ksF1(1:N_spr) ,ksI2(1:N_spr) ,ksIF(1:N_spr)
    real*8 :: l0F1(1:N_spr) ,l0I2(1:N_spr) ,l0IF(1:n_spr)
    
    integer :: csFnsh,csStrt
    integer :: cyclBfr,cyclAft
    real*8  :: weight
    integer :: i,j,jmax
    integer :: cyclNm,csNo
    integer :: sprNm,areaNm
    
    character(len=100) :: flnmks,flnml0,flnmka,flnmA0
    character(len=100) :: flnmbr1,flnmbr2
    character(len=100) :: fullflnm1,fullflnm2
    
    if (weightType==1) then
       flnmks='ksWT1Cycl'
       flnml0='l0WT1Cycl'
       flnmka='kaWT1Cycl'
       flnmA0='A0WT1Cycl'
       
    elseif (weightType==2) then
       flnmks='ksWT2Cycl'
       flnml0='l0WT2Cycl'
       flnmka='kaWT2Cycl'
       flnmA0='A0WT2Cycl'
       
    elseif (weightType==3) then
       flnmks='ksWT3Cycl'
       flnml0='l0WT3Cycl'
       flnmka='kaWT3Cycl'
       flnmA0='A0WT3Cycl'
    endif
    
    do cyclNm = 1,N_addedCycl
       
       write(flnmbr1,'(i1.1)') cyclNm
       
       do csNo = 1,N_CS
          
          write(flnmbr2,'(a,i1.1,a)') 'CS',csNo,'.dat'
          
          do i = 1,2
             
             if (i==1) then
                
                fullflnm1=trim(adjustl(flnmks))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                fullflnm2=trim(adjustl(flnml0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                
                open(unit=46,file=trim(adjustl(fullflnm1)))
                open(unit=47,file=trim(adjustl(fullflnm2)))
                
                jmax = N_spr
                
             elseif (i==2) then
                
                fullflnm1=trim(adjustl(flnmka))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                fullflnm2=trim(adjustl(flnmA0))//trim(adjustl(flnmbr1))//trim(adjustl(flnmbr2))
                
                open(unit=48,file=trim(adjustl(fullflnm1)))
                open(unit=49,file=trim(adjustl(fullflnm2)))
                
                jmax = N_cell
             endif
             
             do j = 1,jmax
                
                if (i==1) then
                   read(46,*) ks_cycl(cyclNm,csNo,j),sprNm
                   read(47,*) l0_cycl(cyclNm,csNo,j),sprNm
                   
                elseif (i==2) then
                   read(48,*) ka_cycl(cyclNm,csNo,j),areaNm
                   read(49,*) A0_cycl(cyclNm,csNo,j),areaNm
                   
                endif
                
             enddo
             
             if (i==1) then
                close(46)
                close(47)
             elseif (i==2) then
                close(48)
                close(49)
             endif
             
          enddo
          
       enddo
       
    enddo
    
    
  end subroutine read_Props_AtConnectingCycls
  
  
  subroutine writeCyclPropsAftRead
    implicit none
    integer :: cyclNm,csNo
    integer :: i,j,jmax
    
    open(unit=36,file='ksCyclVal.dat')
    open(unit=37,file='l0CyclVal.dat')
    open(unit=38,file='kaCyclVal.dat')
    open(unit=39,file='A0CyclVal.dat')
     
    do cyclNm = 1,N_addedCycl
       
       do csNo = 1,N_CS
             
          do i = 1,2
             
             if (i==1) jmax = N_spr
             if (i==2) jmax = N_cell
             
             do j = 1,jmax
                
                if (i==1) then
                   write(36,*) ks_cycl(cyclNm,csNo,j),j
                   write(37,*) l0_cycl(cyclNm,csNo,j),j
                   
                elseif (i==2) then
                   write(38,*) ka_cycl(cyclNm,csNo,j),j
                   write(39,*) A0_cycl(cyclNm,csNo,j),j
                endif
                
             enddo

             if (i==1) then
                write(36,*) " "
                write(37,*) " "
               
             elseif (i==2) then
                write(38,*) " "
                write(39,*) " "
             endif
             
          enddo
          
       enddo
       
    enddo

    close(36)
    close(37)
    close(38)
    close(39)
    
  end subroutine writeCyclPropsAftRead
  
  
  
  
  subroutine destroy_struct
    implicit none
    
    call destroy_blcks
    call deallocate_transfrm_variables
    call deallocate_phiInfo
    
    call deallocate_moving_coordnte_variables
    call deallocate_all_gradient_variables
    call deallocate_all_grdPenF_variables
    
  end subroutine destroy_struct
  
  
  subroutine saveA0l0_ofStrct(struct_No)
    implicit none
    
    integer, intent(in) :: struct_No
    integer :: i,j,imax

    if (modelID==2) then
       write(*,*) "routine: saveA0l0_ofStrct isn't for modelID=2"
       stop
    endif
    
    if (stageNo==4 .and. stageType==1) then
       
       if (struct_No==1) then
          open(unit= 312,file='InitFinlStrucProp1S4T1.dat')
          open(unit= 313,file='strctProps1S4T1.dat')
       elseif (struct_No==2) then
          open(unit=312,file='InitFinlStrucProp2S4T1.dat')
          open(unit=313,file='strctProps2S4T1.dat')
       elseif (struct_No==3) then
          open(unit=312,file='IntrmdStrucPropS4T1.dat')
          open(unit=313,file='strctProps3S4T1.dat')
       endif
       
    elseif (stageNo==4 .and. stageType==2) then
       
       if (struct_No==1) then
          open(unit= 312,file='InitFinlStrucProp1S4T2.dat')
          open(unit= 313,file='strctProps1S4T2.dat')
       elseif (struct_No==2) then
          open(unit=312,file='InitFinlStrucProp2S4T2.dat')
          open(unit=313,file='strctProps2S4T2.dat')
       elseif (struct_No==3) then
          open(unit=312,file='IntrmdStrucPropS4T2.dat')
          open(unit=313,file='strctProps3S4T2.dat')
       endif
       
    elseif (stageNo==1 .and. stageType==1) then
       
       if (struct_No==1) then
          open(unit= 312,file='InitFinlStrucProp1S1T1.dat')
          open(unit= 313,file='strctProps1S1T1.dat')
       elseif (struct_No==2) then
          open(unit=312,file='InitFinlStrucProp2S1T1.dat')
          open(unit=313,file='strctProps2S1T1.dat')
       elseif (struct_No==3) then
          open(unit=312,file='InitFinlStrucProp3S1T1.dat')
          open(unit=313,file='strctProps3S1T1.dat')
       elseif (struct_No==4) then
          open(unit= 312,file='InitFinlStrucProp4S1T1.dat')
          open(unit= 313,file='strctProps4S1T1.dat')
       elseif (struct_No==5) then
          open(unit= 312,file='InitFinlStrucProp5S1T1.dat')
          open(unit= 313,file='strctProps5S1T1.dat')
       elseif (struct_No==6) then
          open(unit= 312,file='InitFinlStrucProp6S1T1.dat')
          open(unit= 313,file='strctProps6S1T1.dat')
       elseif (struct_No==7) then
          open(unit= 312,file='InitFinlStrucProp7S1T1.dat')
          open(unit= 313,file='strctProps7S1T1.dat')
       elseif (struct_No==8) then
          open(unit= 312,file='InitFinlStrucProp8S1T1.dat')
          open(unit= 313,file='strctProps8S1T1.dat')
       elseif (struct_No==9) then
          open(unit= 312,file='InitFinlStrucProp9S1T1.dat')
          open(unit= 313,file='strctProps9S1T1.dat')
       elseif (struct_No==10) then
          open(unit= 312,file='InitFinlStrucProp10S1T1.dat')
          open(unit= 313,file='strctProps10S1T1.dat')
       elseif (struct_No==11) then
          open(unit= 312,file='InitFinlStrucProp11S1T1.dat')
          open(unit= 313,file='strctProps11S1T1.dat')
       elseif (struct_No==12) then
          open(unit= 312,file='InitFinlStrucProp12S1T1.dat')
          open(unit= 313,file='strctProps12S1T1.dat')
       elseif (struct_No==13) then
          open(unit= 312,file='InitFinlStrucProp13S1T1.dat')
          open(unit= 313,file='strctProps13S1T1.dat')
          
       elseif (struct_No==14) then
          open(unit= 312,file='InitFinlStrucProp14S1T1.dat')
          open(unit= 313,file='strctProps14S1T1.dat')
          
       elseif (struct_No==15) then
          open(unit= 312,file='InitFinlStrucProp15S1T1.dat')
          open(unit= 313,file='strctProps15S1T1.dat')  
          
       elseif (struct_No==16) then
          open(unit= 312,file='InitFinlStrucProp16S1T1.dat')
          open(unit= 313,file='strctProps16S1T1.dat')
          
       elseif (struct_No==17) then
          open(unit= 312,file='InitFinlStrucProp17S1T1.dat')
          open(unit= 313,file='strctProps17S1T1.dat')  
          
       elseif (struct_No==18) then
          open(unit= 312,file='InitFinlStrucProp18S1T1.dat')
          open(unit= 313,file='strctProps18S1T1.dat') 
       endif
       
    endif
    
    write(unit=312,fmt=*) struct_No,"structNo"    
    !call sleep(1)
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
    
    l0_strctTN(struct_No,1:N_spr)  = l0(1:N_spr)
    ks_strctTN(struct_No,1:N_spr)  = k_spr(1:N_spr)
    A0_strctTN(struct_No,1:N_cell) = A0(1:N_cell)
    ka_strctTN(struct_No,1:N_cell) = k_area(1:N_cell)
    CgX_strctTN(struct_No,1:N_node) = CgXNode(1:N_node)
    
    if (struct_No==1) then
       N_mvCoordnte1 = N_mvCoordnte
       allocate(coordntesXY_strct1(1:N_mvCoordnte1))
       coordntesXY_strct1 = coordntes_xy
       
       NMCWL0A01 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct1(1:NMCWL0A01))
       coordntes_strct1 = coordntes
       
    elseif (struct_No==2) then
       N_mvCoordnte2 = N_mvCoordnte
       allocate(coordntesXY_strct2(1:N_mvCoordnte2))
       coordntesXY_strct2 = coordntes_xy
       
       NMCWL0A02 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct2(1:NMCWL0A02))
       coordntes_strct2 = coordntes
       
    elseif (struct_No==3) then
       
       if (stageNo==4) then
          continue
       elseif (stageNo==1 .and. stageType==1) then
          N_mvCoordnte3 = N_mvCoordnte
          allocate(coordntesXY_strct3(1:N_mvCoordnte3))
          coordntesXY_strct3 = coordntes_xy
          
          NMCWL0A03 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct3(1:NMCWL0A03))
          coordntes_strct3 = coordntes
       endif
       
    elseif (struct_No==4) then
       N_mvCoordnte4 = N_mvCoordnte
       allocate(coordntesXY_strct4(1:N_mvCoordnte4))
       coordntesXY_strct4 = coordntes_xy
       
       NMCWL0A04 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct4(1:NMCWL0A04))
       coordntes_strct4 = coordntes
       
    elseif (struct_No==5) then
       N_mvCoordnte5 = N_mvCoordnte
       allocate(coordntesXY_strct5(1:N_mvCoordnte5))
       coordntesXY_strct5 = coordntes_xy
       
       NMCWL0A05 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct5(1:NMCWL0A05))
       coordntes_strct5 = coordntes
       
    elseif (struct_No==6) then
       N_mvCoordnte6 = N_mvCoordnte
       allocate(coordntesXY_strct6(1:N_mvCoordnte6))
       coordntesXY_strct6 = coordntes_xy
       
       NMCWL0A06 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct6(1:NMCWL0A06))
       coordntes_strct6 = coordntes
       
    elseif (struct_No==7) then
       N_mvCoordnte7 = N_mvCoordnte
       allocate(coordntesXY_strct7(1:N_mvCoordnte7))
       coordntesXY_strct7 = coordntes_xy
       
       NMCWL0A07 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct7(1:NMCWL0A07))
       coordntes_strct7 = coordntes

    elseif (struct_No==8) then
       N_mvCoordnte8 = N_mvCoordnte
       allocate(coordntesXY_strct8(1:N_mvCoordnte8))
       coordntesXY_strct8 = coordntes_xy
       
       NMCWL0A08 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct8(1:NMCWL0A08))
       coordntes_strct8 = coordntes
       
    elseif (struct_No==9) then
       N_mvCoordnte9 = N_mvCoordnte
       allocate(coordntesXY_strct9(1:N_mvCoordnte9))
       coordntesXY_strct9 = coordntes_xy
       
       NMCWL0A09 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct9(1:NMCWL0A09))
       coordntes_strct9 = coordntes
       
    elseif (struct_No==10) then
       N_mvCoordnte10 = N_mvCoordnte
       allocate(coordntesXY_strct10(1:N_mvCoordnte10))
       coordntesXY_strct10 = coordntes_xy
       
       NMCWL0A010 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct10(1:NMCWL0A010))
       coordntes_strct10 = coordntes
       
    elseif (struct_No==11) then
       N_mvCoordnte11 = N_mvCoordnte
       allocate(coordntesXY_strct11(1:N_mvCoordnte11))
       coordntesXY_strct11 = coordntes_xy
       
       NMCWL0A011 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct11(1:NMCWL0A011))
       coordntes_strct11 = coordntes
       
    elseif (struct_No==12) then
       N_mvCoordnte12 = N_mvCoordnte
       allocate(coordntesXY_strct12(1:N_mvCoordnte12))
       coordntesXY_strct12 = coordntes_xy
       
       NMCWL0A012 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct12(1:NMCWL0A012))
       coordntes_strct12 = coordntes
       
    elseif (struct_No==13) then
       N_mvCoordnte13 = N_mvCoordnte
       allocate(coordntesXY_strct13(1:N_mvCoordnte13))
       coordntesXY_strct13 = coordntes_xy
       
       NMCWL0A013 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct13(1:NMCWL0A013))
       coordntes_strct13 = coordntes
       
    elseif (struct_No==14) then
       N_mvCoordnte14 = N_mvCoordnte
       allocate(coordntesXY_strct14(1:N_mvCoordnte14))
       coordntesXY_strct14 = coordntes_xy
       
       NMCWL0A014 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct14(1:NMCWL0A014))
       coordntes_strct14 = coordntes

    elseif (struct_No==15) then
       N_mvCoordnte15 = N_mvCoordnte
       allocate(coordntesXY_strct15(1:N_mvCoordnte15))
       coordntesXY_strct15 = coordntes_xy
       
       NMCWL0A015 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct15(1:NMCWL0A015))
       coordntes_strct15 = coordntes
       
    elseif (struct_No==16) then
       N_mvCoordnte16 = N_mvCoordnte
       allocate(coordntesXY_strct16(1:N_mvCoordnte16))
       coordntesXY_strct16 = coordntes_xy
       
       NMCWL0A016 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct16(1:NMCWL0A016))
       coordntes_strct16 = coordntes
       
    elseif (struct_No==17) then
       N_mvCoordnte17 = N_mvCoordnte
       allocate(coordntesXY_strct17(1:N_mvCoordnte17))
       coordntesXY_strct17 = coordntes_xy
       
       NMCWL0A017 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct17(1:NMCWL0A017))
       coordntes_strct17 = coordntes
       
    elseif (struct_No==18) then
       N_mvCoordnte18 = N_mvCoordnte
       allocate(coordntesXY_strct18(1:N_mvCoordnte18))
       coordntesXY_strct18 = coordntes_xy
       
       NMCWL0A018 = N_mvCoordnte_withl0A0
       allocate(coordntes_strct18(1:NMCWL0A018))
       coordntes_strct18 = coordntes
       
    endif
    
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          nodeXY_strctTN(struct_No,i,j) = node_xy(i,j)
       enddo
    enddo
    
    do i = 1,2
       do j = 1,N_spr
          if (i==1) write(312,*) "l0_strctTN=",l0_strctTN(struct_No,j),"spr_nm=",j
          if (i==2) write(312,*) "ks_strctTN=",ks_strctTN(struct_No,j),"spr_nm=",j
          if (i==1) write(313,*) ks_strctTN(struct_No,j),l0_strctTN(struct_No,j) !together once
       enddo
       write(312,*) " "
       write(313,*) " "
    enddo
    
    do i = 1,2
       do j = 1,N_cell
          if (i==1) write(312,fmt=*) "A0_strctTN=",A0_strctTN(struct_No,j),"area_nm=",j
          if (i==2) write(312,fmt=*) "ka_strctTN=",ka_strctTN(struct_No,j),"area_nm=",j
          if (i==1) write(313,fmt=*) ka_strctTN(struct_No,j),A0_strctTN(struct_No,j)
       enddo
       write(312,fmt=*) " "
       write(313,fmt=*) " "
    enddo
    
    do i = 1,N_node
       write(312,fmt=*) "CgX_strct=",CgX_strctTN(struct_No,i),"node_nm=",i
       write(313,fmt=*) CgX_strctTN(struct_No,i)
    enddo
    write(312,fmt=*) " "
    write(313,fmt=*) " "
    
    if (struct_No==1)  imax=N_mvCoordnte1
    if (struct_No==2)  imax=N_mvCoordnte2
    if (struct_No==3)  imax=N_mvCoordnte3
    if (struct_No==4)  imax=N_mvCoordnte4
    if (struct_No==5)  imax=N_mvCoordnte5
    if (struct_No==6)  imax=N_mvCoordnte6
    if (struct_No==7)  imax=N_mvCoordnte7
    if (struct_No==8)  imax=N_mvCoordnte8
    if (struct_No==9)  imax=N_mvCoordnte9
    if (struct_No==10) imax=N_mvCoordnte10
    if (struct_No==11) imax=N_mvCoordnte11
    if (struct_No==12) imax=N_mvCoordnte12
    if (struct_No==13) imax=N_mvCoordnte13
    if (struct_No==14) imax=N_mvCoordnte14
    if (struct_No==15) imax=N_mvCoordnte15
    if (struct_No==16) imax=N_mvCoordnte16
    if (struct_No==17) imax=N_mvCoordnte17
    if (struct_No==18) imax=N_mvCoordnte18
    
    do i = 1,imax

       if (struct_No==1) then
          write(312,*) "coordntes_xy=",coordntesXY_strct1(i),"coordntes_XY1=",i
       elseif (struct_No==2) then
          write(312,*) "coordntes_xy=",coordntesXY_strct2(i),"coordntes_XY2=",i
       elseif (struct_No==3) then
          
          if (stageNo==1 .and. stageType==1) then
             write(312,*)"coordntes_xy=",coordntesXY_strct3(i),"coordntes_XY3=",i
          elseif (stageNo==4) then
             continue
          endif
          
       elseif (struct_No==4) then
          write(312,*) "coordntes_xy=",coordntesXY_strct4(i),"coordntes_XY4=",i
       elseif (struct_No==5) then
          write(312,*) "coordntes_xy=",coordntesXY_strct5(i),"coordntes_XY5=",i
       elseif (struct_No==6) then
          write(312,*) "coordntes_xy=",coordntesXY_strct6(i),"coordntes_XY6=",i
       elseif (struct_No==7) then
          write(312,*) "coordntes_xy=",coordntesXY_strct7(i),"coordntes_XY7=",i
       elseif (struct_No==8) then
          write(312,*) "coordntes_xy=",coordntesXY_strct8(i),"coordntes_XY8=",i
       elseif (struct_No==9) then
          write(312,*) "coordntes_xy=",coordntesXY_strct9(i),"coordntes_XY9=",i
       elseif (struct_No==10) then
          write(312,*) "coordntes_xy=",coordntesXY_strct10(i),"coordntes_XY10=",i
       elseif (struct_No==11) then
          write(312,*) "coordntes_xy=",coordntesXY_strct11(i),"coordntes_XY11=",i
       elseif (struct_No==12) then
          write(312,*) "coordntes_xy=",coordntesXY_strct12(i),"coordntes_XY12=",i
       elseif (struct_No==13) then
          write(312,*) "coordntes_xy=",coordntesXY_strct13(i),"coordntes_XY13=",i
       elseif (struct_No==14) then
          write(312,*) "coordntes_xy=",coordntesXY_strct14(i),"coordntes_XY14=",i
       elseif (struct_No==15) then
          write(312,*) "coordntes_xy=",coordntesXY_strct15(i),"coordntes_XY15=",i
       elseif (struct_No==16) then
          write(312,*) "coordntes_xy=",coordntesXY_strct16(i),"coordntes_XY16=",i
       elseif (struct_No==17) then
          write(312,*) "coordntes_xy=",coordntesXY_strct17(i),"coordntes_XY17=",i
       elseif (struct_No==18) then
          write(312,*) "coordntes_xy=",coordntesXY_strct18(i),"coordntes_XY18=",i
       endif
       
       if (struct_No==1) then
          write(313,*) coordntesXY_strct1(i)
       elseif (struct_No==2) then
          write(313,*) coordntesXY_strct2(i)
          
       elseif (struct_No==3) then
          
          if (stageNo==1 .and. stageType==1) then
             write(313,*) coordntesXY_strct3(i)
          elseif (stageNo==4) then
             continue
          endif
          
       elseif (struct_No==4) then
          write(313,*) coordntesXY_strct4(i)
       elseif (struct_No==5) then
          write(313,*) coordntesXY_strct5(i)
       elseif (struct_No==6) then
          write(313,*) coordntesXY_strct6(i)
       elseif (struct_No==7) then
          write(313,*) coordntesXY_strct7(i)
       elseif (struct_No==8) then
          write(313,*) coordntesXY_strct8(i)
       elseif (struct_No==9) then
          write(313,*) coordntesXY_strct9(i)
       elseif (struct_No==10) then
          write(313,*) coordntesXY_strct10(i)
       elseif (struct_No==11) then
          write(313,*) coordntesXY_strct11(i)
       elseif (struct_No==12) then
          write(313,*) coordntesXY_strct12(i)
       elseif (struct_No==13) then
          write(313,*) coordntesXY_strct13(i)
       elseif (struct_No==14) then
          write(313,*) coordntesXY_strct14(i)
       elseif (struct_No==15) then
          write(313,*) coordntesXY_strct15(i)
       elseif (struct_No==16) then
          write(313,*) coordntesXY_strct16(i)
       elseif (struct_No==17) then
          write(313,*) coordntesXY_strct17(i)
       elseif (struct_No==18) then
          write(313,*) coordntesXY_strct18(i)
       endif
       
    enddo
    
    write(312,fmt=*) " "
    write(313,fmt=*) " "
    
    do i = 1,N_node
       write(312,fmt=*) "node_xy=",nodeXY_strctTN(struct_No,i,1:N_dmnsn),"node_nm",i
       write(313,fmt=*) nodeXY_strctTN(struct_No,i,1:N_dmnsn)
    enddo
    
    close(312)
    close(313)
    
  end subroutine saveA0l0_ofStrct
  
  
  subroutine allocate_and_initialize_thrshld_vars
    implicit none
    
    allocate(thrshAnglInfo(1:N_thrshAngl,1:2))
    allocate(thrshAngl(1:N_thrshAngl))
    allocate(curr_thrshAngl(1:N_thrshAngl))
    
    thrshAnglInfo  = -1
    thrshAngl      = -1.0d-30
    curr_thrshAngl = -1.0d-30
    
    thrsh_lgcl = .False.
    
  end subroutine allocate_and_initialize_thrshld_vars
  
  subroutine deallocate_thrshld_vars
    implicit none
    
    deallocate(thrshAnglInfo,thrshAngl,curr_thrshAngl)
    
  end subroutine deallocate_thrshld_vars
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine get_thirdStrctProp(strct_INPT,strct_TRGT)
    
    implicit none
    integer, intent(in) :: strct_INPT,strct_TRGT
    
    integer :: N_step,N_add
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,m,i1
    
    real*8  :: A0_inpt(1:N_cell),A0_trgt(1:N_cell),A0_rnge(1:N_cell)
    real*8  :: l0_inpt(1:N_spr),l0_trgt(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_inpt(1:N_cell),ka_trgt(1:N_spr),ka_rnge(1:N_cell)
    real*8  :: ks_inpt(1:N_spr),ks_trgt(1:N_spr),ks_rnge(1:N_spr)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    
    integer :: funcChoice,N_variabls,iter
    real*8  :: ftol,fret,fpp
    
    integer :: Pulley,LftNeigh,RghtNeigh
    real*8  :: Pulley_dis(1:2)
    real*8  :: tol_Rdfn
    logical :: lgcl_rdfn
    
    real*8  :: Enrgy_val
    integer :: step_no
    real*8  :: hm
    integer :: pairNo
    
    real*8  :: fracChangeLP4(1:2)
    
    inpt_rngeFlag = 1
    
    open(unit=179,file='get_thirdStrctProp.dat')
    open(unit=180,file='Track_thirdStrctChanges.dat')
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          if (modelID==1) node_xy(i,j) = nodeXY_strctTN(strct_INPT,i,j)
          if (modelID==2) node_xy(i,j) = nodeXY_strctNI(strct_INPT,i,j)
       enddo
    enddo
    
    Exprmnt = 16 ; Frame = 1
    
    N_step   = 10 ; N_add = 5
    stepSize = (1.0d0)/N_step
    
    tol_Rdfn = 0.05d0
    
    call get_inpt_rnge(strct_INPT,strct_TRGT,A0_inpt,A0_trgt,A0_rnge,ka_inpt,&
         ka_trgt,ka_rnge,l0_inpt,l0_trgt,l0_rnge,ks_inpt,ks_trgt,ks_rnge)
    call get_Ift(IftS,IftC)
    
    
    do i = 1,(N_step)
       
       step = i*stepSize
       
       if (i.ne.N_step) call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
       
       
       if (i==1) then
          if (modelID==1) CgXNode(1:N_node) = CgX_strctTN(1,1:N_node)
          if (modelID==2) CgXNode(1:N_node) = CgX_strctNI(1,1:N_node)
          
          k_phi(1:N_node,1:max_Phi_node) = 0.00d0
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       if (i.ne.N_step) then
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          call switchto_NI_model_run_and_switchbackto_TN
          
          write(179,fmt=*) (Frame-1),"frame where intrpltn occurs"
          write(*,*) (Frame-1),"frame where intrpltn occurs"
          
          write(180,*) (Frame-1),"bfr Tilt Frame"
          write(180,*) l0(7),l0(19),l0(10),l0(22),"discrep l0s bfr cmpr"
          call cmpr_and_decide_to_reduce_tilt(Exprmnt,Frame,thrshAngl,curr_thrshAngl,thrshAnglInfo,thrsh_lgcl)
          write(180,*) (Frame-1),"aft reduce Tilt Frame"
          write(180,*) l0(7),l0(19),l0(10),l0(22),"discrep l0s aft cmpr"
          
          call Interpltn_adjstmnt_for_spr(l0_inpt,N_step,i,l0_trgt,l0_rnge)
          
          if(i==(N_step-1)) then

             write(180,*) (Frame-1),"bfr LP4(crit spr,ltrl2 of apprch cell) rlx"
             write(180,*) l0(7),l0(19),l0(10),l0(22),"discrep l0s bfr lp4"
             call StepWiseChangingl0_AP1AP2_andthen_BP3_andthenLP4(Exprmnt,Frame,fracChangeLP4)
             write(180,*) (Frame-1),"aft LP4(crit spr,ltrl2 of apprch cell) rlx"
             write(180,*) l0(7),l0(19),l0(10),l0(22),"discrep l0s aft lp4"

             call Interpltn_adjstmnt_for_spr(l0_inpt,N_step,i,l0_trgt,l0_rnge)
          endif
          
       endif
       
       
       lgcl_rdfn = .False.
       call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
       
       if (lgcl_rdfn .eqv. .True.) then
          strctNo = 3
          call updt_kspr_with_l0
          call updt_karea_with_A0
          call saveA0l0_ofStrct(strctNo)
          
          strctNo = 1
          write(*,*) "got upto here"
          !stop
          exit
          
       elseif (lgcl_rdfn .eqv. .False.) then
          continue
          !CgXNode(1:N_node) = CgX_strctTN/NI(1,1:N_node)
       endif
       
       if(i==(N_step-1) .AND. lgcl_rdfn.eqv..True.) then
          call StepWise_unchangingl0_LP4(Exprmnt,Frame,fracChangeLP4)
          call Interpltn_adjstmnt_for_spr(l0_inpt,N_step,i,l0_trgt,l0_rnge)
       endif
       
    enddo
    
    close(179)
    close(180)
    
    write(*,*) lgcl_rdfn,"lgcl_rdfn"
    
  end subroutine get_thirdStrctProp
  
  
  subroutine TwoStep_trnsfrmtn
    implicit none
    integer :: strct_INPT,strct_OTPT
    integer :: i
    
    open(unit=128,file='l0A0_aftRstrc.dat')
    
    Exprmnt = 17
    !Frame = 1 !this will be on for cycle using no loop 
    
    if (CyclNo==1) then
       Frame = 1
    else
       continue
    endif
    
    inpt_rngeFlag = 2
    
    strct_INPT = 1 ; strct_OTPT = 3
    call trnsfrm_btwn_strcts(strct_INPT,strct_OTPT)
    
    open(unit=124,file='Varying_CgXNode.dat',position='append')
    
    do i=1,N_node
       write(124,fmt=*) CgXNode(i),i
    enddo
    
    close(124)
    !stop
    
    call restructure_system
    
    do i = 1,N_spr
       write(128,fmt=*) l0(i),k_spr(i),i
    enddo
    
    write(128,fmt=*) " "
    
    do i = 1,N_cell
       write(128,fmt=*) A0(i),k_area(i),i
    enddo
    
    close(128)
    
    strct_INPT = 3 ; strct_OTPT = 2    
    call trnsfrm_btwn_strcts(strct_INPT,strct_OTPT)
    
  end subroutine TwoStep_trnsfrmtn
  
  
  subroutine trnsfrm_btwn_strcts(strct_INPT,strct_OTPT)
    implicit none
    integer, intent(in) :: strct_INPT
    integer, intent(in) :: strct_OTPT
    
    real*8  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell)
    real*8  :: l0_inpt(1:N_spr),l0_otpt(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_spr),ka_rnge(1:N_cell)
    real*8  :: ks_inpt(1:N_spr),ks_otpt(1:N_spr),ks_rnge(1:N_spr)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          !node_xy(i,j) = nodeXY_strctTN(strct_INPT,i,j)
          if (modelID==1) node_xy(i,j) = nodeXY_strctTN(strct_INPT,i,j)
          if (modelID==2) node_xy(i,j) = nodeXY_strctNI(strct_INPT,i,j)
       enddo
    enddo
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 10
    stepSize = (1.0d0)/N_step
    
    call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt,&
         ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    call get_Ift(IftS,IftC)
    
    
    open(unit=86,file='trnsfrm_btwn_strct.dat')
    
    N_itm = 9
    
    do i = 1,N_itm
       
       if (i.le.4) jmax=N_spr
       if ((i.gt.4) .and. (i.le.8)) jmax=N_cell
       if ((i.gt.8) .and. (i.le.10)) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctTN(1,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctTN(3,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctTN(1,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctTN(3,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctTN(1,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctTN(3,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctTN(1,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctTN(3,j),j
             if (i==9) write(86,fmt=*) node_xy(j,1:2),nodeXY_strctTN(1,j,1:2)
             
          elseif (modelID==2) then
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctNI(1,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctNI(3,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctNI(1,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctNI(3,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctNI(1,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctNI(3,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctNI(1,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctNI(3,j),j
             if (i==9) write(86,fmt=*) node_xy(j,1:2),nodeXY_strctNI(1,j,1:2)
          endif
          
       enddo
       
       write(86,fmt=*) " "
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
       
       if (i==0) then
          
          if (strct_INPT==1) then
             if (modelID==1) CgXNode(1:N_node) = CgX_strctTN(strct_INPT,1:N_node)
             if (modelID==2) CgXNode(1:N_node) = CgX_strctNI(strct_INPT,1:N_node)
          elseif (strct_INPT==3) then
             if (modelID==1) CgXNode(1:N_node) = CgX_strctTN(strct_OTPT,1:N_node)
             if (modelID==2) CgXNode(1:N_node) = CgX_strctNI(strct_OTPT,1:N_node)
          endif
          
          k_phi(1:N_node,1:max_Phi_node) = 0.00d0
          
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)      
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo    
    
  end subroutine trnsfrm_btwn_strcts
  
  subroutine trnsfrm_btwn_strcts_S1(strct_INPT,strct_OTPT)
    implicit none
    integer, intent(in) :: strct_INPT
    integer, intent(in) :: strct_OTPT
    
    real*8  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell)
    real*8  :: l0_inpt(1:N_spr),l0_otpt(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_spr),ka_rnge(1:N_cell)
    real*8  :: ks_inpt(1:N_spr),ks_otpt(1:N_spr),ks_rnge(1:N_spr)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 10
    stepSize = (1.0d0)/N_step
    
    call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt,&
         ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    call get_Ift(IftS,IftC)
    
    
    open(unit=86,file='trnsfrm_btwn_strct_S1.dat')
    
    N_itm = 9
    
    do i = 1,N_itm
       
       if (i.le.4) jmax=N_spr
       if ((i.gt.4) .and. (i.le.8)) jmax=N_cell
       if ((i.gt.8) .and. (i.le.10)) jmax=N_node
       
       do j = 1,jmax

          if (modelID==1) then
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctTN(strct_INPT,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctTN(strct_OTPT,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctTN(strct_INPT,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctTN(strct_OTPT,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctTN(strct_INPT,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctTN(strct_OTPT,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctTN(strct_INPT,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctTN(strct_OTPT,j),j
             if (i==9) write(86,*) node_xy(j,1:2),nodeXY_strctTN(strct_INPT,j,1:2)
             
          elseif (modelID==2) then
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctNI(strct_INPT,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctNI(strct_OTPT,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctNI(strct_INPT,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctNI(strct_OTPT,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctNI(strct_INPT,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctNI(strct_OTPT,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctNI(strct_INPT,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctNI(strct_OTPT,j),j
             if (i==9) write(86,*) node_xy(j,1:2),nodeXY_strctNI(strct_INPT,j,1:2)
          endif
             
       enddo
       
       write(86,fmt=*) " "
    enddo
    
    close(86)
    !stop
    
    do i = 0,(N_step)
       
       step = i*stepSize
          
       call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
       
       if (i==0) then
          if (modelID==1) CgXNode(1:N_node) = CgX_strctTN(strct_INPT,1:N_node)
          if (modelID==2) CgXNode(1:N_node) = CgX_strctNI(strct_INPT,1:N_node)
          
          k_phi(1:N_node,1:max_Phi_node) = 0.00d0
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
    
  end subroutine trnsfrm_btwn_strcts_S1
  
  
  subroutine trnsfrm_btwn_strcts_S1_withPressAdjstmnt(strct_INPT,strct_OTPT,seqNo)
    implicit none
    integer, intent(in) :: strct_INPT
    integer, intent(in) :: strct_OTPT
    integer, intent(in) :: seqNo
    
    real*8 , allocatable :: A0_inpt(:),A0_otpt(:),A0_rnge(:)
    real*8 , allocatable :: l0_inpt(:),l0_otpt(:),l0_rnge(:)
    real*8 , allocatable :: ka_inpt(:),ka_otpt(:),ka_rnge(:)
    real*8 , allocatable :: ks_inpt(:),ks_otpt(:),ks_rnge(:)
    integer, allocatable :: IftS(:),IftC(:)
    real*8 , allocatable :: TnsnVal(:),PresVal(:) 
    
    integer :: N_step,N_itm
    integer :: cnt
    real*8  :: stepSize,step,E
    integer :: i,j,imax,jmax,i1,m,m1
    integer :: strct1,strct2,strct3,strct4,strctToRead,modelWntd
    integer :: strctFrm,strctTo
    integer :: CS1,CS2,CS3,CS4
    real*8  :: rsFctr
    integer :: howManyPBC,IfBotmMostIncluded
    integer :: howManyPFC,cntrlSt
    integer :: howManyIC !IC=Incoming Cells
    real*8  :: hmuch
    integer :: saveWhom,LastRglrCell
    integer :: copyDirectn,TopOrBottom
    integer :: readingDcsn,cell_no
    integer :: RghtstrtTN
    integer :: cellNm
    
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
    
    !open(unit=87,file='PressAdjstmntNIChk.dat',position='append')
    
    if (seqNo.lt.9) then
       write(*,*)"routine: trnfrm_btwn_strcts_S1_withPressAdjstmt for seqNo>=9"
       stop
    endif
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    if (strct_INPT==13) then
      call set_lower_bndry_to_AvgPostn_and_Equilibrate(strct_INPT,Exprmnt_NI,Frame_NI)
      rsFctr = 1.00d0
      call rescaling_the_system(strct_INPT,rsFctr,Exprmnt_NI,Frame_NI)
    endif
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    allocate(IftS(1:N_spr),IftC(1:N_cell))
    allocate(A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell))
    allocate(ka_inpt(1:N_cell),ka_otpt(1:N_cell),ka_rnge(1:N_cell))
    allocate(l0_inpt(1:N_spr), l0_otpt(1:N_spr) ,l0_rnge(1:N_spr))
    allocate(ks_inpt(1:N_spr), ks_otpt(1:N_spr) ,ks_rnge(1:N_spr))
    
    N_step   = 2
    stepSize = (1.0d0)/N_step
    
    call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt,&
         ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    CgXNode(1:N_node)              = CgX_strctTN(strct_INPT,1:N_node)
    CgYNode(1:N_node)              = 0.0d0
    k_phi(1:N_node,1:max_Phi_node) = 0.00d0
    
    
    if (hlfCycl==1) then
       
       write(*,*) (Frame_NI-1),"Bfr CS1"
       call NI_switch_and_read_strctProps_withAddedCellNI(strct_INPT,A0_inpt,ka_inpt,l0_inpt,ks_inpt)
       saveWhom = 3
       call saveA0l0_ofStrct_AddedCell_TNorNI(strct_INPT,saveWhom)
       write(*,*) (Frame_NI-1),"Aft CS1";write(*,*) (Frame_NI-1),"Bfr CS2"
       
       call NI_switch_and_read_strctProps_withAddedCellNI(strct_OTPT,A0_otpt,ka_otpt,l0_otpt,ks_otpt)
       saveWhom = 3
       call saveA0l0_ofStrct_AddedCell_TNorNI(strct_OTPT,saveWhom)
       write(*,*) (Frame_NI-1),"Aft CS2"
       call print_the_properties
       
    elseif (hlfCycl==2) then
       
       write(*,*) (Frame_NI-1),"Bfr CS3"
       call NI_switch_and_borrow_strctPropsfrom_HlfCycl1(strct_INPT,A0_inpt,ka_inpt,l0_inpt,ks_inpt)
       saveWhom = 3
       call saveA0l0_ofStrct_AddedCell_TNorNI(strct_INPT,saveWhom)
       write(*,*) (Frame_NI-1),"Aft CS3"
       
       write(*,*) (Frame_NI-1),"Bfr CS4"
       call NI_switch_and_read_strctProps_withAddedCellNI(strct_OTPT,A0_otpt,ka_otpt,l0_otpt,ks_otpt)
       saveWhom = 3
       call saveA0l0_ofStrct_AddedCell_TNorNI(strct_OTPT,saveWhom)
       write(*,*) (Frame_NI-1),"Aft CS4"
       
    endif
    
    deallocate(IftS,IftC)
    deallocate(A0_inpt,A0_otpt,A0_rnge)
    deallocate(ka_inpt,ka_otpt,ka_rnge)
    deallocate(l0_inpt,l0_otpt,l0_rnge)
    deallocate(ks_inpt,ks_otpt,ks_rnge)
    
    !close(87)
    
    do m = 1,2 !1 for TN, 1 for NI
       
       modelID = m
       
       allocate(IftS(1:N_spr),IftC(1:N_cell))
       allocate(A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell))
       allocate(ka_inpt(1:N_cell),ka_otpt(1:N_cell),ka_rnge(1:N_cell))
       allocate(l0_inpt(1:N_spr), l0_otpt(1:N_spr) ,l0_rnge(1:N_spr))
       allocate(ks_inpt(1:N_spr), ks_otpt(1:N_spr) ,ks_rnge(1:N_spr))
       
      
       
       call get_Ift(IftS,IftC)
       call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt, &
            ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
       
       A0     = A0_inpt
       k_spr  = ks_inpt
       l0     = l0_inpt
       k_area = ka_inpt
       
       if (m==1) then
          N_step=2  ; stepSize = (1.0d0)/N_step
       elseif (m==2) then
          N_step=2  ; stepSize = (1.0d0)/N_step
       endif
       
       do i = 0,(N_step)
          
          step = i*stepSize
          call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
          
          if (i==0) then
             
             if (modelID==1) then
                
                CgXNode(1:N_node) = CgX_strctTN(strct_INPT,1:N_node)
                CgYNode(1:N_node) = 0.0d0
                node_xy(1:N_node,1:N_dmnsn) = nodeXY_strctTN(strct_INPT,1:N_node,1:N_dmnsn)
                call nodes_to_coordntes(node_xy,coordntes_xy)
                
             elseif (modelID==2) then
                
                CgXNode(1:N_node) = CgX_strctNI(strct_INPT,1:N_node)
                CgYNode(1:N_node) = 0.0d0
                node_xy(1:N_node,1:N_dmnsn) = nodeXY_strctNI(strct_INPT,1:N_node,1:N_dmnsn)
                call nodes_to_coordntes(node_xy,coordntes_xy)
                
                write(*,*) strct_INPT,strct_OTPT,"strct"
                do m1 = 1,N_node
                   write(*,*) node_xy(m1,1:2),m1,"m1 node"
                enddo
                
             endif
             
             k_phi(1:N_node,1:max_Phi_node)   = 0.00d0
             
          endif
          
          call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
          call Equilibrate_system
          
          if(i==0) write(*,*) "check_here HEREEEEEEE"
          
          if (m==1) then
             call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
             Frame = Frame+1
          elseif (m==2) then
             call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
             Frame_NI = Frame_NI+1
             
             !if ((Frame_NI-1)==9) then
             !   write(*,*) "Frame_NI=",Frame_NI-1
             !   write(*,*) "Entered in 9"
             !   call diagonal_tnsn_rls_test
             !   write(*,*) "Stopping in build_strctures, diag_test"
             !   stop
             !endif
          endif
          
          if (m==2) then
             
             if (i==0) then
                
                if (hlfCycl == 1) cntrlSt = 1
                if (hlfCycl == 2) cntrlSt = 3
                
                allocate(TnsnVal(1:N_spr),PresVal(1:N_cell))
                TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
                PresVal(1:N_cell) = Pressure(node_xy,A0)
                call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
                deallocate(TnsnVal,PresVal)
                
             elseif (i==N_step) then
                
                if (hlfCycl == 1) cntrlSt = 2
                if (hlfCycl == 2) cntrlSt = 4
                
                allocate(TnsnVal(1:N_spr),PresVal(1:N_cell))
                TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
                PresVal(1:N_cell) = Pressure(node_xy,A0)
                call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
                deallocate(TnsnVal,PresVal)
                
             endif
             
          endif
          
       enddo
       
       if (m==1) then
          call switch_to_NI_model
          E = Energy(node_xy,l0,A0)
          call get_gradient(node_xy,l0,A0,gradient)
       elseif (m==2) then
          
          if (seqNo == 10) then !this condn is for continuing to AddedCell model,if not clear it
             continue
          else
             call deallocate_repetitive_arrays
             call switchback_to_TN_model
          endif
          
       endif
       
       deallocate(IftS,IftC)
       deallocate(A0_inpt,A0_otpt,A0_rnge)
       deallocate(ka_inpt,ka_otpt,ka_rnge)
       deallocate(l0_inpt,l0_otpt,l0_rnge)
       deallocate(ks_inpt,ks_otpt,ks_rnge)
       
    enddo
    
    if (strct_OTPT==16) then
       continue
    endif
    
    if (seqNo==9) then
       !call diagonal_tnsn_rls_test 
       !write(*,*) "stopping_after_the_test" ; stop
       call boundary_force_test
       write(*,*) "stopping_after_the_boundary_force_test" ; stop
       continue
    endif
    
    if (seqNo==10) then
       cellNm=5 ; call make_control_state_perfect_for_supp_doc(cellNm)
       call print_the_properties
       write(*,*) "stopping Aft Seq=10"
       stop
    endif
    
  end subroutine trnsfrm_btwn_strcts_S1_withPressAdjstmnt
  
  
  subroutine trnsfrm_btwn_strcts_inTwoSteps_S1(strct_INPT,strct_OTPT)
    implicit none
    integer, intent(in) :: strct_INPT
    integer, intent(in) :: strct_OTPT
    
    real*8  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell)
    real*8  :: l0_inpt(1:N_spr),l0_otpt(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_spr),ka_rnge(1:N_cell)
    real*8  :: ks_inpt(1:N_spr),ks_otpt(1:N_spr),ks_rnge(1:N_spr)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    integer :: cellNm
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 2
    stepSize = (1.0d0)/N_step
    
    call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt,&
         ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    call get_Ift(IftS,IftC)
    
    open(unit=86,file='trnsfrm_btwn_strct_S1.dat')
    
    N_itm = 9
    
    do i = 1,N_itm
       
       if (i.le.4) jmax=N_spr
       if ((i.gt.4) .and. (i.le.8)) jmax=N_cell
       if ((i.gt.8) .and. (i.le.10)) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctTN(strct_INPT,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctTN(strct_OTPT,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctTN(strct_INPT,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctTN(strct_OTPT,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctTN(strct_INPT,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctTN(strct_OTPT,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctTN(strct_INPT,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctTN(strct_OTPT,j),j
             if (i==9) write(86,fmt=*) node_xy(j,1:2),nodeXY_strctTN(strct_INPT,j,1:2)
             
          elseif (modelID==2) then
             
             if (i==1) write(86,fmt=*) l0_inpt(j),l0_strctNI(strct_INPT,j),j
             if (i==2) write(86,fmt=*) l0_otpt(j),l0_strctNI(strct_OTPT,j),j
             if (i==3) write(86,fmt=*) ks_inpt(j),ks_strctNI(strct_INPT,j),j
             if (i==4) write(86,fmt=*) ks_otpt(j),ks_strctNI(strct_OTPT,j),j
             if (i==5) write(86,fmt=*) A0_inpt(j),A0_strctNI(strct_INPT,j),j
             if (i==6) write(86,fmt=*) A0_otpt(j),A0_strctNI(strct_OTPT,j),j
             if (i==7) write(86,fmt=*) ka_inpt(j),ka_strctNI(strct_INPT,j),j
             if (i==8) write(86,fmt=*) ka_otpt(j),ka_strctNI(strct_OTPT,j),j
             if (i==9) write(86,*) node_xy(j,1:2),nodeXY_strctNI(strct_INPT,j,1:2)
             
          endif
          
       enddo
       
       write(86,fmt=*) " "
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
       
       if (i==0) then
          if (modelID==1) CgXNode(1:N_node) = CgX_strctTN(strct_INPT,1:N_node)
          if (modelID==2) CgXNode(1:N_node) = CgX_strctNI(strct_INPT,1:N_node)
          
          CgYNode(1:N_node)              = 0.00d0
          k_phi(1:N_node,1:max_Phi_node) = 0.00d0
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       !if (i==0) then
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          !call switchto_NI_model_run_and_switchbackto_TN !(uncomment if you wanna run)
       !else
       !   call Equilibrate_bothTN_NImodel_with_PreSimltd_DatFile(Exprmnt,Frame)
       !endif
       
    enddo
    
  end subroutine trnsfrm_btwn_strcts_inTwoSteps_S1
  
  
  subroutine trnsfrm_btwn_strcts_inTwoSteps_S1_with_condn(hlfCycl,strct_INPT,strct_OTPT)
    implicit none
    integer, intent(in) :: hlfCycl
    integer, intent(in) :: strct_INPT
    integer, intent(in) :: strct_OTPT
    
    real*8  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell)
    real*8  :: l0_inpt(1:N_spr),l0_otpt(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_inpt(1:N_cell),ka_otpt(1:N_spr),ka_rnge(1:N_cell)
    real*8  :: ks_inpt(1:N_spr),ks_otpt(1:N_spr),ks_rnge(1:N_spr)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1,j1,j1max
    
    integer :: demlishingSpr
    real*8  :: initialL,tolrncL
    integer :: cnt_lp
    logical :: lgcl_condn
    real*8  :: diff(1:2),E
    integer :: BfrOrAft
    integer :: chk_condn_needed=0
    integer :: fileNo
    integer :: primNodeS,primNodeE,intrNodeS,intrNodeE
    integer :: CMv,PCSorICS
    integer :: CellsMeetV
    
    real*8, allocatable :: gradA(:,:),gradN(:,:)
    real*8, allocatable :: grdPnt1(:,:),grdPnt2(:,:),grdDiff(:,:)
    real*8, allocatable :: kaSaveV(:),CgXSaveV(:),CgYSaveV(:)
    
    real*8              :: diff_nodeXY(1:150,1:2)=10.0d30
    real*8, allocatable :: ksv1NI(:),kav1NI(:),l0v1NI(:),A0v1NI(:),CgXv1NI(:),CgYv1NI(:)
    real*8, allocatable :: ksDiff(:),kaDiff(:),l0Diff(:),A0Diff(:),CgXDiff(:),CgYDiff(:)
    
    write(*,*) hlfCycl,strct_INPT,strct_OTPT,"hlfCycl-strct_INPT-strct_OTPT"
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 2
    stepSize = (1.0d0)/N_step
    
    call get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,ka_inpt,&
         ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    call get_Ift(IftS,IftC)
    
    open(unit=86,file='trnsfrm_btwn_strct_S1condn.dat')
    
    N_itm = 9
    demlishingSpr = (N_cell-1)*3 + 1
    
    primNodeS=1 ; primNodeE=48 ; intrNodeS=51 ; intrNodeE=140
    
    do i = 1,N_itm
       
       if (i.le.4)                   jmax=N_spr
       if ((i.gt.4) .and. (i.le.8))  jmax=N_cell
       if ((i.gt.8) .and. (i.le.10)) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             
             if (i==1) write(86,*) l0_inpt(j),l0_strctTN(strct_INPT,j),j
             if (i==2) write(86,*) l0_otpt(j),l0_strctTN(strct_OTPT,j),j
             if (i==3) write(86,*) ks_inpt(j),ks_strctTN(strct_INPT,j),j
             if (i==4) write(86,*) ks_otpt(j),ks_strctTN(strct_OTPT,j),j
             if (i==5) write(86,*) A0_inpt(j),A0_strctTN(strct_INPT,j),j
             if (i==6) write(86,*) A0_otpt(j),A0_strctTN(strct_OTPT,j),j
             if (i==7) write(86,*) ka_inpt(j),ka_strctTN(strct_INPT,j),j
             if (i==8) write(86,*) ka_otpt(j),ka_strctTN(strct_OTPT,j),j
             if (i==9) write(86,*) node_xy(j,1:2),nodeXY_strctTN(strct_INPT,j,1:2)
             
          elseif (modelID==2) then
             
             if (i==1) write(86,*) l0_inpt(j),l0_strctNI(strct_INPT,j),j
             if (i==2) write(86,*) l0_otpt(j),l0_strctNI(strct_OTPT,j),j
             if (i==3) write(86,*) ks_inpt(j),ks_strctNI(strct_INPT,j),j
             if (i==4) write(86,*) ks_otpt(j),ks_strctNI(strct_OTPT,j),j
             if (i==5) write(86,*) A0_inpt(j),A0_strctNI(strct_INPT,j),j
             if (i==6) write(86,*) A0_otpt(j),A0_strctNI(strct_OTPT,j),j
             if (i==7) write(86,*) ka_inpt(j),ka_strctNI(strct_INPT,j),j
             if (i==8) write(86,*) ka_otpt(j),ka_strctNI(strct_OTPT,j),j
             if (i==9) write(86,*) node_xy(j,1:2),nodeXY_strctNI(strct_INPT,j,1:2)
             
          endif
          
       enddo
       
       write(86,*) " "
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
       
       if (i==0) then
          
          if (modelID==1) then
             CgXNode(1:N_node) = CgX_strctTN(strct_INPT,1:N_node)
             CgYNode(1:N_node) = CgY_strctTN(strct_INPT,1:N_node)
          elseif (modelID==2) then
             CgXNode(1:N_node) = CgX_strctNI(strct_INPT,1:N_node)
             CgYNode(1:N_node) = CgY_strctTN(strct_INPT,1:N_node)
          endif
          
          k_phi(1:N_node,1:max_Phi_node) = 0.00d0
       endif
       
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_bothTN_NImodel(Exprmnt,Frame)
       
       if (hlfCycl==1) then 
          
          if (i==0) then
             call redevelop_the_initiation_Phase()
             stop 'stopping aft redevelop_the_initiation_Phase'
          endif
          
       endif
       
    enddo
    
    write(*,*) "Stop here: sb:trnsfrm_btwn_strcts_inTwoSteps_S1_with_condn"
    !stop
    
  end subroutine trnsfrm_btwn_strcts_inTwoSteps_S1_with_condn
  
  
  subroutine redevelop_the_initiation_Phase()
    implicit none
    integer :: redevlop_PCSorICS
    
    call taken_from_inside_of_tbsitsS1wc !This routin is for CM0 PCS to CM1ICS/CM2ICS with Force
    !call taken_from_inside_of_tbsitsS1wc_to_modelVF
    
    call redevlop_the_PCS_of_initiationPhase()
    stop 'aft redev'
    
  end subroutine redevelop_the_initiation_Phase
  
  
  subroutine redevlop_the_PCS_of_initiationPhase()
    implicit none
    integer :: CellsMeetV,PCSorICS
    integer :: EXEorREAD,CMv,strctNum
    integer :: i,j,jmax
    integer :: READFRMFILE
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    ! EXPLANATION OF WHAT THE ROUTINE DOES  !!!!
    ! EXEorREAD=1 [goes to manpltns_needtodo_for_cntrlStates & do manipulations-(notNeeded anymore)]
    ! EXEorREAD=2 [read data files, e.g (sprPrps/cellPrps..)_CM2PCS.dat, and continues from there]
    ! READFRMFILE=0 -->  will take us to the get_fix_the_PCS routine and generate the files, e.g.,
    !                     (sprPrps/convSprPrps/cellPrps/nodePrps..)_CM(1/2/3/4)PCS_HighKs.dat
    ! READFRMFILE=1 -->  will take us to the get_fix_the_PCS routine and continue with the already
    !                    generated files
    
    
    EXEorREAD=2 ; CMv=CellsMeet ; PCSorICS=1 ; write(*,*)CellsMeet,CMv,PCSorICS,"11~"   !CM1PCS
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    READFRMFILE=1 ; call get_fix_the_PCS(CMv,PCSorICS,READFRMFILE)
    
    EXEorREAD=2 ; CMv=CellsMeet+1 ; PCSorICS=2 ; write(*,*)CellsMeet,CMv,PCSorICS,"22~" !CM2ICS
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS); write(*,*) "IS IT HERE"
    READFRMFILE=1 ; call get_fix_the_ICS(CMv,PCSorICS,READFRMFILE)
    CellsMeetV=CellsMeet;PCSorICS=2;call adjusting_routines_for_NI_systems(CellsMeetV,PCSorICS)
    !call productionRunForIntiatnPhase(CMv,PCSorICS)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    EXEorREAD=2 ; CMv=CellsMeet ; PCSorICS=1 ; write(*,*)CellsMeet,PCSorICS,"CM-PCSorICS 2"!CM2PCS
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    READFRMFILE=1 ; call get_fix_the_PCS(CMv,PCSorICS,READFRMFILE)
    
    
    EXEorREAD=2 ; CMv=CellsMeet+1 ; PCSorICS=2 ; write(*,*) CellsMeet,PCSorICS,"CM-PCSor(ICS)"
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    READFRMFILE=1 ; call get_fix_the_ICS(CMv,PCSorICS,READFRMFILE)
    CellsMeetV=CellsMeet;PCSorICS=2;call adjusting_routines_for_NI_systems(CellsMeetV,PCSorICS)
    !call productionRunForIntiatnPhase(CMv,PCSorICS)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    EXEorREAD=2; CMv=CellsMeet ; PCSorICS=1 ; write(*,*) CellsMeet,PCSorICS,"CM-PCSorICS 3"!CM3PCS
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    READFRMFILE=1 ; call get_fix_the_PCS(CMv,PCSorICS,READFRMFILE)
    
    EXEorREAD=2 ; CMv=CellsMeet+1 ; PCSorICS=2 ; write(*,*) CellsMeet,PCSorICS,"CM-PCSor(ICS)"
    call manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    READFRMFILE=1 ; call get_fix_the_ICS(CMv,PCSorICS,READFRMFILE)
    CellsMeetV=CellsMeet;PCSorICS=2;call adjusting_routines_for_NI_systems(CellsMeetV,PCSorICS)
    !call productionRunForIntiatnPhase(CMv,PCSorICS)
    
    call replaceProps_however_retain_only_Cg();PCSorICS=1;write(*,*)CellsMeet,PCSorICS,"CM4PC"
    call Equilibrate_only_NI_model;write(*,*)"FrmNI aftReplcePropsRetain_only_Cg=",(Frame_NI-1)
    call modify_prp_frm_CM3PCS_wo_Highks_value_and_change_stpBystp
    
    READFRMFILE=1 ; call get_fix_the_PCS(CMv,PCSorICS,READFRMFILE)
    
  end subroutine redevlop_the_PCS_of_initiationPhase
  
  
  subroutine redevelop_the_ICS_of_initiationPhase()
    implicit none
    
    
  end subroutine redevelop_the_ICS_of_initiationPhase
  
  !subroutine redev_blckwise_design(whichCM,PCSorICS)
  !  implicit none
  !  integer :: readFrmFiles(1:maxNumCM)
    
  !  call get_the_readfrmFiles_val()
  !  call 
    
  !end subroutine redev_blckwise_design
  
  
  subroutine get_fix_the_PCS(CMv,PCSorICS,readFrmFile)
    implicit none
    integer, intent(in) :: CMv,PCSorICS,readFrmFile
    
    real*8  :: ks_str(1:N_spr)  ,l0_str(1:N_spr)
    real*8  :: ka_str(1:N_cell) ,A0_str(1:N_cell)
    real*8  :: CgX_str(1:N_node),CgY_str(1:N_node)
    real*8  :: nodeXY_str(1:N_node,1:N_dmnsn)
    integer :: lpVal
    real*8  :: CgXChng,TIIF_val
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! STORING props for general case
    
    ks_str(1:N_spr)           = k_spr(1:N_spr)
    l0_str(1:N_spr)           = l0(1:N_spr)
    ka_str(1:N_cell)          = k_area(1:N_cell)
    A0_str(1:N_cell)          = A0(1:N_cell)
    CgX_str(1:N_node)         = CgXNode(1:N_node)
    CgY_str(1:N_node)         = CgYNode(1:N_node)
    nodeXY_Str(1:N_node,1:2)  = node_xy(1:N_node,1:2)
    
    if (readFrmFile==0) then
       
       call higher_kspr_tst_with_same_tension_length_and_adjstd_restLength(CMv,PCSorICS)
       call save_the_prps_of_High_ks_tst(CMv,PCSorICS)
       
       write(*,*) "Frame_NI after save =", (Frame_NI-1)
       stop 'aft save'
       
    elseif (readFrmFile==1) then
       
       call read_the_prps_of_High_ks_tst(CMv,PCSorICS)
       call Equilibrate_only_NI_model 
       
       if (CMv==4) then
          
          
          ! call rscale_basedOn_bndryPress()
          
          ! *************************
          ! To APPLY BOUNDARY Tnsn & check if the system can return based on the node position (y23)
          do lpVal = 1,4
             call pullAtBndry_and_toChkIf_the_SystmCanReturn(CMv,lpVal)
          enddo
          ! *************************
          
          !call make_equivalent_TIIF_for_all_CMPCS(CMv)
          !TIIF_val = 0.50d0 ; call make_InptTIIF_for_all_CMPCS_frmZeroTIIF(CMv,TIIF_val)
          
          !call take_equivlnt_TIIF_systm_corrct_IC_nodePstn_and_curvtr(CMv)
          !call rscale_basedOn_bndryPress() ; write(*,*) Frame_NI-1,"aft rscle aft take_eq ... sb"
          
          !TIIF_val=0.50d0 ; call TIIF_adjstmnt_ToASpcfcVal(CMv,TIIF_val)
          !write(*,*) Frame_NI-1,"aft TIIF adjst"
          
          !call read_TIIFfile_props_and_apply_force(CMv)
          
          !call remove_vrtclForce_and_toChkIf_the_PCScanreturn(CMv)
          
          !call apply_forces_at_ApclSrfc_to_create_flatBase
          
          !call take_currnt_dats_and_edit_forces_to_match(CMv)
          
          !call save_the_prps_of_High_ks_tst_buffr(CMv,PCSorICS)
          
          if (CMv==4) then
             stop 'aft CMv=4'
          endif
          
       endif
       
       !call productionRunForIntiatnPhase(CMv,PCSorICS)
       
    endif
    
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! PERTURBATION TESTS HERE
    
    !if (CMv==4 .and. PCSorICS==1) call Y_force_perturbation_to_PCS
    !call displace_pulley_virtually_to_check_force ; stop 'aft displce' 
    
    !if (PCSorICS==1) then
    !   call SensitvtyTst_vrtclForce_onTrmnlNode_ofIC_Apcl
    !endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k_spr(1:N_spr)        = ks_str(1:N_spr)   
    l0(1:N_spr)           = l0_str(1:N_spr)          
    k_area(1:N_cell)      = ka_str(1:N_cell)
    A0(1:N_cell)          = A0_str(1:N_cell)
    CgXNode(1:N_node)     = CgX_str(1:N_node)
    CgYNode(1:N_node)     = CgY_str(1:N_node)   
    node_xy(1:N_node,1:2) = nodeXY_Str(1:N_node,1:2)
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call Equilibrate_only_NI_model
    
  end subroutine get_fix_the_PCS
  
  
  subroutine get_fix_the_ICS(CMv,PCSorICS,readFrmFile)
    implicit none
    integer, intent(in) :: CMv,PCSorICS,readFrmFile
    
    real*8  :: ks_str(1:N_spr)  ,l0_str(1:N_spr)
    real*8  :: ka_str(1:N_cell) ,A0_str(1:N_cell)
    real*8  :: CgX_str(1:N_node),CgY_str(1:N_node)
    real*8  :: nodeXY_str(1:N_node,1:N_dmnsn)
    integer :: InitPCSfl,TrmnlPCSfl,N_step
    real*8  :: ksI(1:N_spr),l0I(1:N_spr),kaI(1:N_cell),A0I(1:N_cell),CgXI(1:N_node),CgYI(1:N_node)
    real*8  :: ksT(1:N_spr),l0T(1:N_spr),kaT(1:N_cell),A0T(1:N_cell),CgXT(1:N_node),CgYT(1:N_node)
    integer :: PCSorICS_val
    
    real*8, allocatable :: InitFlProp(:),TrmnlFlProp(:)

    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! STORING props for general case
    
    ks_str(1:N_spr)           = k_spr(1:N_spr)
    l0_str(1:N_spr)           = l0(1:N_spr)
    ka_str(1:N_cell)          = k_area(1:N_cell)
    A0_str(1:N_cell)          = A0(1:N_cell)
    CgX_str(1:N_node)         = CgXNode(1:N_node)
    CgY_str(1:N_node)         = CgYNode(1:N_node)
    nodeXY_Str(1:N_node,1:2)  = node_xy(1:N_node,1:2)
    
    write(*,*) CMv,PCSorICS,readFrmFile,"must check val"
    
    if (readFrmFile==0) then
    
       InitPCSfl=CMv-1 ; TrmnlPCSfl=CMv ; write(*,*) InitPCSfl,TrmnlPCSfl,CMv,"Init-Trmn-CMv"
       
       PCSorICS_val=1
       call read_the_prps_of_High_ks_tst_PCSorICS(InitPCSfl, PCSorICS_val,ksI,l0I,kaI,A0I,CgXI,CgYI)
       call read_the_prps_of_High_ks_tst_PCSorICS(TrmnlPCSfl,PCSorICS_val,ksT,l0T,kaT,A0T,CgXT,CgYT)
       call read_the_nodePrps_frm_PCS_data(InitPCSfl,PCSorICS_val)
       call print_aftReadfrm_PCSfiles(ksI,ksT,l0I,l0T,kaI,kaT,A0I,A0T,CgXI,CgXT,CgYI,CgYT)
       
       N_step = 10
       call trnsfrm_btwn_PCS_strcts_toFindICS(N_step,ksI,ksT,l0I,l0T,kaI,kaT,A0I,A0T,&
            CgXI,CgXT,CgYI,CgYT)
       
       call save_the_prps_of_High_ks_tst(CMv,PCSorICS)
       write(*,*)"Frame_NI after save =",(Frame_NI-1);write(*,*)CMv,PCSorICS;stop 'aft save in ICS'
       
    elseif (readFrmFile==1) then
       
       call read_the_prps_of_High_ks_tst(CMv,PCSorICS)
       call Equilibrate_only_NI_model
       
    endif
    
    k_spr(1:N_spr)        = ks_str(1:N_spr)   
    l0(1:N_spr)           = l0_str(1:N_spr)          
    k_area(1:N_cell)      = ka_str(1:N_cell)
    A0(1:N_cell)          = A0_str(1:N_cell)
    CgXNode(1:N_node)     = CgX_str(1:N_node)
    CgYNode(1:N_node)     = CgY_str(1:N_node)   
    node_xy(1:N_node,1:2) = nodeXY_Str(1:N_node,1:2)
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call Equilibrate_only_NI_model
    
  end subroutine get_fix_the_ICS
  
  
  subroutine trnsfrm_btwn_PCS_strcts_toFindICS(N_step,ksI,ksT,l0I,l0T,kaI,kaT,A0I,A0T,&
       CgXI,CgXT,CgYI,CgYT)
    
    implicit none
    integer, intent(in)  :: N_step
    real*8 , intent(in)  :: ksI(1:N_spr)  ,ksT(1:N_spr)  ,l0I(1:N_spr)  ,l0T(1:N_spr)
    real*8 , intent(in)  :: kaI(1:N_spr)  ,kaT(1:N_spr)  ,A0I(1:N_spr)  ,A0T(1:N_spr)
    real*8 , intent(in)  :: CgXI(1:N_node),CgXT(1:N_node),CgYI(1:N_node),CgYT(1:N_node)
    
    real*8               :: ksRnge(1:N_spr),  l0Rnge(1:N_spr),   kaRnge(1:N_cell)
    real*8               :: A0Rnge(1:N_cell), CgXRnge(1:N_node), CgYRnge(1:N_node) 
    !real*8               :: RngeCombined(1:2*(N_spr+N_cell+N_node))
    integer              :: IFTspr(1:N_spr),  IFTcell(1:N_cell), IFTnode(1:N_node)
    !integer, allocatable :: IFTcombined(:,:)
    real*8               :: stepSize,stepFrctn,tol_Rdfn
    logical              :: lgcl_Rdfn
    integer              :: Hlf_Nstep
    
    ksRnge(1:N_spr)   = ksT(1:N_spr)   - ksI(1:N_spr)   !; RngeCombined(1:N_spr)         = 
    l0Rnge(1:N_spr)   = l0T(1:N_spr)   - l0I(1:N_spr)   !; RngeCombined(N_spr+1:2*N_spr) =
    kaRnge(1:N_cell)  = kaT(1:N_cell)  - kaI(1:N_cell)  !; RngeCombined(2*N_spr+1:)
    A0Rnge(1:N_cell)  = A0T(1:N_cell)  - A0I(1:N_cell)  !;
    CgXRnge(1:N_node) = CgXT(1:N_node) - CgXI(1:N_node) !;
    CgYRnge(1:N_node) = CgYT(1:N_node) - CgYI(1:N_node) !;
    
    call print_inpt_otpt_rnge(ksI,ksT,ksRnge,l0I,l0T,l0Rnge,kaI,kaT,kaRnge,&
       A0I,A0T,A0Rnge,CgXI,CgXT,CgXRnge,CgYI,CgYT,CgYRnge)
    
    IFTspr(1:N_spr)   = 1 !; IFTcombined(1:N_spr)                              = IFTspr(1:N_spr)
    IFTcell(1:N_cell) = 1 !; IFTcombined((N_spr+1):(N_spr+N_cell))             = IFTcell(1:N_cell)
    IFTnode(1:N_node) = 1 !; IFTcombined((N_spr+N_cell+1):(N_spr+N_cell+N_node)= IFTnode(1:N_node)
    
    stepSize          = (1.0000d0)/(N_step) ; write(*,*) N_step,stepSize,"N_step,stepSize"
    Hlf_Nstep         = N_step/2            ; write(*,*) N_step,Hlf_Nstep,"N_step,Hlf_Nstep"
    
    do i = 0,(N_step)
       stepFrctn = (i)*(stepSize)
       call stepUp_variables_excpt_DiagSpr(i,stepFrctn,Iftspr,IFTcell,IFTnode,&
            A0I,A0Rnge,kaI,kaRnge,l0I,l0Rnge,ksI,ksRnge,CgXI,CgXRnge,CgYI,CgYRnge)
       call Equilibrate_only_NI_model ; write(*,*) (Frame_NI-1),i,"bfr applying bndry F"
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa();write(*,*)Frame_NI-1,i,"aft applng bndryF"
       
       if (i==0) call make_diagonal_spr_tnsn_reduced_toFindICS !; stop 'aft Hlf_Nstep'
       if (i==0) write(*,*) "aft making Tension close to zero"
       
       lgcl_Rdfn=.False. ; tol_Rdfn=0.07d0
       call get_decision_of_redefining_diffWay(tol_Rdfn,lgcl_Rdfn)
       
       if (lgcl_rdfn .eqv. .True.) then
          write(*,*) "lgcl_Rdfn is True with CMv =",CellsMeet+1 ; exit
       endif
       
    enddo
    
  end subroutine trnsfrm_btwn_PCS_strcts_toFindICS
  
  subroutine print_inpt_otpt_rnge(ksI,ksT,ksRnge,l0I,l0T,l0Rnge,kaI,kaT,kaRnge,&
       A0I,A0T,A0Rnge,CgXI,CgXT,CgXRnge,CgYI,CgYT,CgYRnge)
    implicit none
    real*8, intent(in) :: ksI(1:N_spr),  ksT(1:N_spr),  ksRnge(1:N_spr)
    real*8, intent(in) :: l0I(1:N_spr),  l0T(1:N_spr),  l0Rnge(1:N_spr)
    real*8, intent(in) :: kaI(1:N_cell), kaT(1:N_cell), kaRnge(1:N_cell)
    real*8, intent(in) :: A0I(1:N_cell), A0T(1:N_cell), A0Rnge(1:N_cell)
    real*8, intent(in) :: CgXI(1:N_node),CgXT(1:N_node),CgXRnge(1:N_node)
    real*8, intent(in) :: CgYI(1:N_node),CgYT(1:N_node),CgYRnge(1:N_node)
    integer            :: i,j,jmax,numElmntToPrint
    
    open(unit=459,file='print_props_InitFinlRNge.dat')
    
    numElmntToPrint = 6
    
    do i = 1,numElmntToPrint
       
       if (i==1 .or. i==2) jmax=N_spr
       if (i==3 .or. i==4) jmax=N_cell
       if (i==5 .or. i==6) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(459,*) ksT(j), ksI(j), ksRnge(j), j,"ks"
          if (i==2) write(459,*) l0T(j), l0I(j), l0Rnge(j), j,"l0"
          if (i==3) write(459,*) kaT(j), kaI(j), kaRnge(j), j,"ka"
          if (i==4) write(459,*) A0T(j), A0I(j), A0Rnge(j), j,"A0"
          if (i==5) write(459,*) CgXT(j),CgXI(j),CgXRnge(j),j,"CgX"
          if (i==6) write(459,*) CgYT(j),CgYI(j),CgYRnge(j),j,"CgY"
       enddo
       
       write(459,*) " "
       
    enddo
    
    close(459)
    
  end subroutine print_inpt_otpt_rnge
  
  
  subroutine print_aftReadfrm_PCSfiles(ksI,ksT,l0I,l0T,kaI,kaT,A0I,A0T,CgXI,CgXT,CgYI,CgYT)
    implicit none
    real*8, intent(in) :: ksI(1:N_spr),  ksT(1:N_spr)
    real*8, intent(in) :: l0I(1:N_spr),  l0T(1:N_spr)
    real*8, intent(in) :: kaI(1:N_cell), kaT(1:N_cell)
    real*8, intent(in) :: A0I(1:N_cell), A0T(1:N_cell)
    real*8, intent(in) :: CgXI(1:N_node),CgXT(1:N_node)
    real*8, intent(in) :: CgYI(1:N_node),CgYT(1:N_node)
    
    integer            :: i,j,jmax,numElmntToPrint
    
    open(unit=261,file='print_props_aftReadfrmPCS.dat')
    
    numElmntToPrint = 6
    
    do i = 1,numElmntToPrint
       
       if (i==1 .or. i==2) jmax=N_spr
       if (i==3 .or. i==4) jmax=N_cell
       if (i==5 .or. i==6) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(261,*) ksT(j), ksI(j), j,"ks"
          if (i==2) write(261,*) l0T(j), l0I(j), j,"l0"
          if (i==3) write(261,*) kaT(j), kaI(j), j,"ka"
          if (i==4) write(261,*) A0T(j), A0I(j), j,"A0"
          if (i==5) write(261,*) CgXT(j),CgXI(j),j,"CgX"
          if (i==6) write(261,*) CgYT(j),CgYI(j),j,"CgY"
       enddo
       
       write(261,*) " "
       
    enddo
    
    close(261)
    
  end subroutine print_aftReadfrm_PCSfiles
  
  
  
  !subroutine combining_multipl_1DimArrays_into_2DimArray(numArrays,maxNumElmnt)
  !  implicit none
  !  integer, intent(in) :: numArrays,maxNumElmnt
  !  integer, intent(in) :: array1(1:max),a
    
    
  !end subroutine combining_arrays
  
  subroutine trnsfrm_btwn_strcts_PT(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: ks_rnge(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell),A0_rnge(1:N_cell)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    real*8  :: Esv,EaV
    integer :: itrnNo,Progcycl
    integer :: strctV1,strctV2
    integer :: strghtOrNI
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 10
    stepSize = (1.0d0)/N_step
    
    imax          = 2
    NI_alrdySaved = 0
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
          endif
          
       enddo
       
    enddo
    
    call get_Ift(IftS,IftC)
    
    open(unit=86,file='trnsfrm_btwn_strct_PT.dat')
    
    N_itm = 5
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.5)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(86,*) ks_I(j),ks_F(j),j
          if (i==2) write(86,*) l0_I(j),l0_F(j),j
          if (i==3) write(86,*) ka_I(j),ka_F(j),j
          if (i==4) write(86,*) A0_I(j),A0_F(j),j
          if (i==5) write(86,*) CgXNode(j),j
       enddo
       
       write(86,fmt=*) " "
       
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_I,A0_rnge,ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       if (i.gt.0) then
          
          open(unit=232,file="Echeck_PTT.dat")
          
          EsV=0.0d0 ; EaV=0.0d0
          
          do j = 1,N_spr
             EsV = EsV + k_spr(j)*((l(j)-l0(j))**2)
             write(232,*) k_spr(j),l(j),l0(j),EsV,j,"EsV"
          enddo
          
          write(232,*) " "
          
          do j = 1,N_cell
             EaV = EaV + k_area(j)*((A(j)-A0(j))**2)
             write(232,*) k_area(j),A(j),A0(j),EaV,j,"EaV"
          enddo
          
          close(232)
          
       endif
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       
       if (NI_saveOrNot == 0) then
          call switchto_NI_model_run_and_switchbackto_TN
          
       elseif (NI_saveOrNot == 1) then
          
          ProgCycl   = ProgCyclNo
          strctV1    = strctVno1
          strctV2    = strctVno2
          strghtOrNI = 2
          itrnNo     = i
          
          
          if (NI_alrdySaved==0) then
             call switchto_NI_model_save_and_switchbackto_TN(ProgCycl,strctV1,strctV2,strghtOrNI,itrnNo)
             
          elseif (NI_alrdySaved==1) then
             continue
          endif
          
       endif
       
       if (i==N_step) then
          NI_alrdySaved = 0
       endif
       
    enddo
    
  end subroutine trnsfrm_btwn_strcts_PT
  
  
  subroutine trnsfrm_btwn_strct_PT_woNI(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer,intent(in) :: strghtOrNI
    
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: ks_rnge(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell),A0_rnge(1:N_cell)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1,j1,i1max,j1max
    real*8  :: Esv,EaV
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step   = 10
    stepSize = (1.0d0)/N_step
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
          endif
          
       enddo
       
    enddo
    
    call get_Ift(IftS,IftC)
    
    open(unit=86,file='trnsfrm_btwn_strct_PT_woNI.dat')
    
    N_itm = 5
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.5)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(86,*) ks_I(j),ks_F(j),j
          if (i==2) write(86,*) l0_I(j),l0_F(j),j
          if (i==3) write(86,*) ka_I(j),ka_F(j),j
          if (i==4) write(86,*) A0_I(j),A0_F(j),j
          if (i==5) write(86,*) CgXNode(j),j
       enddo
       
       write(86,fmt=*) " "
       
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_I,A0_rnge,ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       if (i.gt.0) then
          
          open(unit=232,file="Echeck_PTT.dat")
          
          EsV=0.0d0 ; EaV=0.0d0
          
          do j = 1,N_spr
             EsV = EsV + k_spr(j)*((l(j)-l0(j))**2)
             write(232,*) k_spr(j),l(j),l0(j),EsV,j,"EsV"
          enddo
          
          write(232,*) " "
          
          do j = 1,N_cell
             EaV = EaV + k_area(j)*((A(j)-A0(j))**2)
             write(232,*) k_area(j),A(j),A0(j),EaV,j,"EaV"
          enddo
          
          close(232)
          
       endif
       
       
       if (strghtOrNI==2) then
          
          if (Frame_NI==45) then
             
             open(unit=14,file='DiscrepencyInProgCycl.dat',position='append')
             
             i1max = 3
             
             do i1 = 1,i1max
                
                if (i1==1) j1max=N_spr
                if (i1==2) j1max=N_cell
                if (i1==3) j1max=N_node
                
                do j1 = 1,j1max
                   if (i1==1) write(14,*) k_spr(j1),l(j1),l0(j1),j1
                   if (i1==2) write(14,*) k_area(j1),A(j1),A0(j1),j1
                   if (i1==3) write(14,*) CgXNode(j1),node_xy(j1,1:2),j1
                enddo
                
                write(14,*) " "
                
             enddo
             
             
             close(14)
             
          endif
       endif
       
       call Equilibrate_system
       
       if (strghtOrNI==1) then
          
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          
       elseif (strghtOrNI==2) then
          
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          
          if (Frame_NI==44) then
             
             open(unit=14,file='DiscrepencyInProgCycl.dat',position='append')
             
             i1max = 3
             
             do i1 = 1,i1max
                
                if (i1==1) j1max=N_spr
                if (i1==2) j1max=N_cell
                if (i1==3) j1max=N_node
                
                do j1 = 1,j1max
                   if (i1==1) write(14,*) k_spr(j1),l(j1),l0(j1),j1
                   if (i1==2) write(14,*) k_area(j1),A(j1),A0(j1),j1
                   if (i1==3) write(14,*) CgXNode(j1),node_xy(j1,1:2),j1
                enddo
                
                write(14,*) " "
                
             enddo
             
             close(14)
             
          endif
          
          !if (Frame_NI==45) stop
          
          Frame_NI = Frame_NI+1
          
       endif
       
    enddo
    
    
  end subroutine trnsfrm_btwn_strct_PT_woNI
  
  
  subroutine trnsfrm_LittleBitReverseStrct3to4(strghtOrNI,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer, intent(in) :: strghtOrNI
    
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: ks_rnge(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell),A0_rnge(1:N_cell)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: N_step,N_stepToRun
    integer :: N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1,j1,i1max,j1max
    real*8  :: Esv,EaV 
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step   = 50 ; N_stepToRun = 15
    stepSize = (-1.0d0)/N_step
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
          endif
          
       enddo
       
    enddo
    
    call get_Ift(IftS,IftC)
    
    open(unit=84,file='trnsfrm_LittleBitS3toS4.dat')
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.5)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(84,*) ks_I(j),ks_F(j),j
          if (i==2) write(84,*) l0_I(j),l0_F(j),j
          if (i==3) write(84,*) ka_I(j),ka_F(j),j
          if (i==4) write(84,*) A0_I(j),A0_F(j),j
          if (i==5) write(84,*) CgXNode(j),j
       enddo
       
       write(86,fmt=*) " "
       
    enddo
    
    close(84)
    
    
    do i = 0,(N_stepToRun)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_I,A0_rnge,ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       call Equilibrate_system
       
       if (strghtOrNI==1) then
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
          Frame = Frame+1
          
       elseif (strghtOrNI==2) then
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          Frame_NI = Frame_NI+1
       endif
       
    enddo
    
  end subroutine trnsfrm_LittleBitReverseStrct3to4
  
  
  subroutine trnsfrm_btwn_strcts_AddedCellNI(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: ks_rnge(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell),A0_rnge(1:N_cell)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    real*8  :: Esv,EaV
    integer :: itrnNo,Progcycl
    integer :: strctV1,strctV2
    integer :: strghtOrNI
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    if (arbtry_NstepDcsn==0) N_step = 10
    if (arbtry_NstepDcsn==1) N_step = 6
    
    stepSize = (1.0d0)/N_step
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
          endif
          
       enddo
       
    enddo
    
    call get_Ift(IftS,IftC)
    
    open(unit=86,file='trnsfrm_btwn_strct_ACNI.dat')
    
    N_itm = 5
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.5)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(86,*) ks_I(j),ks_F(j),j
          if (i==2) write(86,*) l0_I(j),l0_F(j),j
          if (i==3) write(86,*) ka_I(j),ka_F(j),j
          if (i==4) write(86,*) A0_I(j),A0_F(j),j
          if (i==5) write(86,*) CgXNode(j),j
       enddo
       
       write(86,fmt=*) " "
       
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_I,A0_rnge,ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExprmntAdded,FrameAdded)
       FrameAdded = FrameAdded+1
       
    enddo
    
  end subroutine trnsfrm_btwn_strcts_AddedCellNI
  
  
  
  subroutine trnsfrm_btwn_strcts_varying_kskal0A0CgXCgY(ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F,&
       CgX_I,CgX_F,CgY_I,CgY_F)
    implicit none
    real*8, intent(in) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(in) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(in) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(in) :: A0_I(1:N_cell),A0_F(1:N_cell)
    real*8, intent(in) :: CgX_I(1:N_node),CgX_F(1:N_node)
    real*8, intent(in) :: CgY_I(1:N_node),CgY_F(1:N_node)
    
    real*8  :: ks_rnge(1:N_spr),  l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell), A0_rnge(1:N_cell)
    real*8  :: CgX_rnge(1:N_node),CgY_rnge(1:N_node)
    
    integer :: IftS(1:N_spr),IftC(1:N_cell),IftN(1:N_node)
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    real*8  :: Esv,EaV
    integer :: itrnNo,strctV1,strctV2
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    if (arbtry_NstepDcsn==0) N_step = 10
    if (arbtry_NstepDcsn==1) N_step = 6
    
    stepSize = (1.0d0)/N_step
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
             
          elseif (i==3) then
             CgX_rnge(j) = CgX_F(j)-CgX_I(j)
             CgY_rnge(j) = CgY_F(j)-CgY_I(j)
          endif
          
       enddo
       
    enddo
    
    !call get_Ift(IftS,IftC)
    
    IftS(1:N_spr)=1 ; IftC(1:N_cell)=1 ; IftN(1:N_node)=1
    
    open(unit=86,file='trnsfrm_btwn_strct_ACNI.dat')
    
    N_itm = 6
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.6)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(86,*) ks_I(j), ks_F(j),j
          if (i==2) write(86,*) l0_I(j), l0_F(j),j
          if (i==3) write(86,*) ka_I(j), ka_F(j),j
          if (i==4) write(86,*) A0_I(j), A0_F(j),j
          if (i==5) write(86,*) CgX_I(j),CgX_F(j),j
          if (i==6) write(86,*) CgY_I(j),CgY_F(j),j
       enddo
       
       write(86,fmt=*) " "
       
    enddo
    
    close(86)
    
    do i = 0,(N_step)
       
       step = i*stepSize
       
       call stepUp_variables_includingCg(i,step,IftS,IftC,IftN,A0_I,A0_rnge,&
            ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge,CgX_I,CgX_rnge,CgY_I,CgY_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_Init,FramePrdctnInit)
       FramePrdctnInit = FramePrdctnInit+1
       
    enddo
    
  end subroutine trnsfrm_btwn_strcts_varying_kskal0A0CgXCgY
  
  
  subroutine resave_NI_strctV3_ifNeeded(ProgCycl,strghtOrNI)
    implicit none
    integer, intent(in) :: ProgCycl
    integer, intent(in) :: strghtOrNI
    
    real*8, allocatable :: A0_I(:),A0_F(:),l0_I(:),l0_F(:)
    real*8, allocatable :: ka_I(:),ka_F(:),ks_I(:),ks_F(:)
    real*8, allocatable :: CgX_V(:)
    
    integer :: strctV
    real*8  :: E
    integer :: i,j,jmax
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    write(*,*) E,"E"
    call get_gradient(node_xy,l0,A0,gradient)
    
    allocate(A0_I(1:N_cell),A0_F(1:N_cell))
    allocate(l0_I(1:N_spr) ,l0_F(1:N_spr) )
    allocate(ka_I(1:N_cell),ka_F(1:N_cell))
    allocate(ks_I(1:N_spr) ,ks_F(1:N_spr) )
    allocate(CgX_V(1:N_node))
    
    strctV = 3 !; strghtOrNI = 2
    call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ks_I,l0_I,ka_I,A0_I,CgX_V)
    CgXNode(1:N_node) = CgX_V(1:N_node) ; CgYNode(1:N_node) = 0.0d0
    strctV = 4 !; strghtOrNI = 2
    call readProgrsnProps(ProgCycl,strctV,strghtOrNI,ks_F,l0_F,ka_F,A0_F,CgX_V)
    
    open(unit=96,file='compare_revPropChk_bfr.dat')
    open(unit=97,file='compare_revStrct3Prop_bfr.dat')
    
    do i = 1,2
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             write(96,*) k_spr(j),l0(j),j
             write(97,*) ks_I(j),l0_I(j),j
             
          elseif (i==2) then
             write(96,*) k_area(j),A0(j),j
             write(97,*) ka_I(j),A0_I(j),j
          endif
          
       enddo
       
       write(96,*) " "
       write(97,*) " "
       
    enddo
    
    close(96)
    close(97)
    
    call compare_revMovment_and_update(ProgCycl,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    
    strctV=3 !; strghtOrNI = 2
    call save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0_I,l0_I,ka_I,ks_I)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model 
    
  end subroutine resave_NI_strctV3_ifNeeded
  
  
  
  
  subroutine compare_revMovment_and_update(ProgCycl,ks_I,ks_F,l0_I,l0_F,ka_I,ka_F,A0_I,A0_F)
    implicit none
    integer, intent(in)   :: ProgCycl
    real*8, intent(inout) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    real*8, intent(inout) :: l0_I(1:N_spr) ,l0_F(1:N_spr)
    real*8, intent(inout) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8, intent(inout) :: A0_I(1:N_cell),A0_F(1:N_cell)
    
    real*8  :: ks_rnge(1:N_spr),l0_rnge(1:N_spr)
    real*8  :: ka_rnge(1:N_cell),A0_rnge(1:N_cell)
    integer :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: N_step,N_itm
    real*8  :: stepSize,step
    integer :: i,j,imax,jmax,i1
    
    integer :: stepOne
    real*8  :: y0,y1,ydif
    integer :: chkNode
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    N_step = 10 ; stepOne = 1
    stepSize = (1.0d0)/N_step
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             ks_rnge(j) = ks_F(j)-ks_I(j)
             l0_rnge(j) = l0_F(j)-l0_I(j)
             
          elseif (i==2) then
             ka_rnge(j) = ka_F(j)-ka_I(j)
             A0_rnge(j) = A0_F(j)-A0_I(j)
          endif
          
       enddo
       
    enddo
    
    
    call get_Ift(IftS,IftC)
    
    open(unit=196,file='compare_revMovment.dat')
    
    N_itm = 5
    
    do i = 1,N_itm
       
       if (i.le.2) jmax=N_spr
       if ((i.gt.2) .and. (i.le.4)) jmax=N_cell
       if ((i.gt.4) .and. (i.le.5)) jmax=N_node
       
       do j = 1,jmax
          if (i==1) write(196,*) ks_I(j),ks_F(j),j
          if (i==2) write(196,*) l0_I(j),l0_F(j),j
          if (i==3) write(196,*) ka_I(j),ka_F(j),j
          if (i==4) write(196,*) A0_I(j),A0_F(j),j
          if (i==5) write(196,*) CgXNode(j),j
       enddo
       
       write(196,fmt=*) " "
       
    enddo
    
    close(196)
    
    
    do i = 0,StepOne !(N_step)
       
       step = i*stepSize
       
       call stepUp_variables(i,step,IftS,IftC,A0_I,A0_rnge,ka_I,ka_rnge,l0_I,l0_rnge,ks_I,ks_rnge)
       call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       
       Frame_NI = Frame_NI+1
       
       
       open(unit=43,file='DN_chk.dat',position='append')
       write(43,*) double_node(1,1),"dn"
       
       chkNode = 13 - 2*(ProgCycl-1)
       
       if (i==0) y0 = node_xy(chkNode,2)
       if (i==1) y1 = node_xy(chkNode,2)
       
       write(43,*) chkNode,y0,y1,"chkNode-y0-y1"
       close(43)
       
    enddo
    
    ydif = y0-y1
    
    if (ydif .ge. 0.00) then
       continue
    elseif (ydif .lt. 0.00) then
       ks_I(1:N_spr)  = k_spr(1:N_spr)
       l0_I(1:N_spr)  = l0(1:N_spr)
       ka_I(1:N_cell) = k_area(1:N_cell)
       A0_I(1:N_cell) = A0(1:N_cell)
    endif
    
    k_spr(1:N_spr)   = ks_F(1:N_spr)
    l0(1:N_spr)      = l0_F(1:N_spr)
    k_area(1:N_cell) = ka_F(1:N_cell)
    A0(1:N_cell)     = A0_F(1:N_cell)

    open(unit=92,file='compare_revPropChk_aft.dat')
    open(unit=97,file='compare_revStrct3Prop_aft.dat')
    
    do i = 1,2
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             write(92,*) k_spr(j),l0(j),j
             write(97,*) ks_I(j),l0_I(j),j
             
          elseif (i==2) then
             write(92,*) k_area(j),A0(j),j
             write(97,*) ka_I(j),A0_I(j),j
          endif
          
       enddo
       
       write(92,*) " "
       write(97,*) " "
       
    enddo
    
    close(92)
    close(97)
    
  end subroutine compare_revMovment_and_update
  
  
  subroutine restructure_system
    implicit none
    integer :: Pulley,LftNeigh,RghtNeigh
    
    call get_Pulley_and_Neighbours(Pulley,LftNeigh,RghtNeigh)
    call deallocate_rdfnmod_vars
    call allocate_rdfnmod_vars
    call taking_lftrghtPulleyNeigh_downwards(Pulley,LftNeigh,RghtNeigh)
    call save_nodeXY   
    call redefine_system
    call redefine_Prop_strct3
    call redefine_Props
    
    call deallocate_thrshld_vars
    call allocate_and_initialize_thrshld_vars
    call get_thrshold_angleInfo(coordntes_xy,thrshAnglInfo,thrshAngl)
    
  end subroutine restructure_system
  
  subroutine stepUp_variables(stepNo,step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
    implicit none
    integer,intent(in) :: stepNo
    real*8, intent(in) :: step
    real*8, intent(in) :: A0_inpt(1:N_cell),A0_rnge(1:N_cell)
    real*8, intent(in) :: ka_inpt(1:N_cell),ka_rnge(1:N_cell)
    real*8, intent(in) :: l0_inpt(1:N_spr),l0_rnge(1:N_spr)
    real*8, intent(in) :: ks_inpt(1:N_spr),ks_rnge(1:N_spr)
    integer,intent(in) :: IftS(1:N_spr),IftC(1:N_cell) 
    
    real*8  :: stepSpr(1:N_spr),stepCell(1:N_cell)
    integer :: i,j,imax,jmax
    real*8  :: EsV,EaV
    
    
    call get_step(step,stepNo,IftS,IftC,stepSpr,stepCell)
    
    open(unit=130,file='step.dat',position='append')
    write(130,fmt=*) step,"step"
    
    do i=1,N_spr
       write(130,fmt=*) IftS(i),stepSpr(i),i,"IftS,stepSpr,spr_nm"
    enddo
    
    do i=1,N_cell
       write(130,fmt=*) IftC(i),stepCell(i),i,"IftC,stepCell,cell_nm"
    enddo
    close(130)
    
    
    open(unit=131,file='InStepVars.dat',position='append')
    write(131,fmt=*) stepNo,"stepNo"
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          
          if(i==1) then
             A0(j) = A0_inpt(j) + stepCell(j)*A0_rnge(j)
             k_area(j) = ka_inpt(j) + stepCell(j)*ka_rnge(j)
          elseif(i==2) then
             l0(j) = l0_inpt(j) + stepSpr(j)*l0_rnge(j)
             k_spr(j) = ks_inpt(j) + stepSpr(j)*ks_rnge(j)
          endif
          
          if (i==1) write(131,fmt=*) A0(j),k_area(j),j,"A0,k_area,cell_no"
          if (i==2) write(131,fmt=*) l0(j),k_spr(j) ,j,"l0,k_spr, spr_no"
       enddo
       
       write(131,fmt=*) " "
    enddo
    
    close(131)
    
  end subroutine stepUp_variables
  
  
  subroutine stepUp_variables_includingCg(stepNo,step,IftS,IftC,IftN,A0_inpt,A0_rnge,&
       ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge,CgX_inpt,CgX_rnge,CgY_inpt,CgY_rnge)
    implicit none
    integer,intent(in) :: stepNo
    real*8, intent(in) :: step
    integer,intent(in) :: IftS(1:N_spr),IftC(1:N_cell),IftN(1:N_node)
    real*8, intent(in) :: A0_inpt(1:N_cell), A0_rnge(1:N_cell)
    real*8, intent(in) :: ka_inpt(1:N_cell), ka_rnge(1:N_cell)
    real*8, intent(in) :: l0_inpt(1:N_spr),  l0_rnge(1:N_spr)
    real*8, intent(in) :: ks_inpt(1:N_spr),  ks_rnge(1:N_spr)
    real*8, intent(in) :: CgX_inpt(1:N_node),CgX_rnge(1:N_node)
    real*8, intent(in) :: CgY_inpt(1:N_node),CgY_rnge(1:N_node)
    
    real*8  :: stepSpr(1:N_spr),stepCell(1:N_cell),stepNode(1:N_node)
    integer :: i,j,imax,jmax
    real*8  :: EsV,EaV
    
    
    call get_step_includingCg(step,stepNo,IftS,IftC,IftN,stepSpr,stepCell,stepNode)
    
    open(unit=132,file='InStepVarsIncludingCg.dat',position='append')
    write(132,*) stepNo,"stepNo" ; write(132,*) " "
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if(i==1) then
             l0(j)       = l0_inpt(j)  + stepSpr(j)  * l0_rnge(j)
             k_spr(j)    = ks_inpt(j)  + stepSpr(j)  * ks_rnge(j)
             
          elseif(i==2) then
             A0(j)       = A0_inpt(j)  + stepCell(j) * A0_rnge(j)
             k_area(j)   = ka_inpt(j)  + stepCell(j) * ka_rnge(j)
             
          elseif (i==3) then
             CgXNode(j)  = CgX_inpt(j) + stepNode(j) * CgX_rnge(j)
             CgYNode(j)  = CgY_inpt(j) + stepNode(j) * CgY_rnge(j)
          endif
          
          if (i==1) write(132,*) l0(j),     l0_inpt(j), l0_rnge(j),j, "l0Info spr_no"
          if (i==1) write(132,*) k_spr(j),  ks_inpt(j), ks_rnge(j),j, "ksInfo spr no"
          
          if (i==2) write(132,*) A0(j),     A0_inpt(j), A0_rnge(j),j, "A0Info cell_no"
          if (i==2) write(132,*) k_area(j), ka_inpt(j), ka_rnge(j),j, "kaInfo cell no"
          
          if (i==3) write(132,*) CgXNode(j),CgX_inpt(j),CgX_rnge(j),j,"CgXInfo node_no"
          if (i==3) write(132,*) CgYNode(j),CgY_inpt(j),CgY_rnge(j),j,"CgYInfo node_no"
          
       enddo
       
       write(132,*) " "
    enddo
    
    close(132)
    
  end subroutine stepUp_variables_includingCg
  
  
  subroutine stepUp_variables_excpt_DiagSpr(stepNo,step,IftS,IftC,IftN,A0_inpt,A0_rnge,&
       ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge,CgX_inpt,CgX_rnge,CgY_inpt,CgY_rnge)
    implicit none
    integer,intent(in) :: stepNo
    real*8, intent(in) :: step
    integer,intent(in) :: IftS(1:N_spr),IftC(1:N_cell),IftN(1:N_node)
    real*8, intent(in) :: A0_inpt(1:N_cell), A0_rnge(1:N_cell)
    real*8, intent(in) :: ka_inpt(1:N_cell), ka_rnge(1:N_cell)
    real*8, intent(in) :: l0_inpt(1:N_spr),  l0_rnge(1:N_spr)
    real*8, intent(in) :: ks_inpt(1:N_spr),  ks_rnge(1:N_spr)
    real*8, intent(in) :: CgX_inpt(1:N_node),CgX_rnge(1:N_node)
    real*8, intent(in) :: CgY_inpt(1:N_node),CgY_rnge(1:N_node)
    
    real*8  :: stepSpr(1:N_spr),stepCell(1:N_cell),stepNode(1:N_node)
    integer :: diagsprL(1:(NAEC_Ltrl+1)),diagsprR(1:(NAEC_Ltrl+1))
    real*8  :: l0_diagsprLstr(1:(NAEC_Ltrl+1)),l0_diagsprRstr(1:(NAEC_Ltrl+1))
    real*8  :: ks_diagsprLstr(1:(NAEC_Ltrl+1)),ks_diagsprRstr(1:(NAEC_Ltrl+1))
    integer :: i,j,imax,jmax
    real*8  :: EsV,EaV
    
    if (stepNo==0) call find_the_invagniating_and_invaginated_cell
    call get_the_diagSprs_and_strTheDiagPrps(diagsprL,diagsprR,l0_diagsprLstr,l0_diagsprRstr,&
         ks_diagsprLstr,ks_diagsprRstr)
    
    call get_step_includingCg(step,stepNo,IftS,IftC,IftN,stepSpr,stepCell,stepNode)
    
    open(unit=132,file='InStepVarsIncludingCg.dat',position='append')
    write(132,*) stepNo,"stepNo" ; write(132,*) " "
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if(i==1) then
             l0(j)       = l0_inpt(j)  + stepSpr(j)  * l0_rnge(j)
             k_spr(j)    = ks_inpt(j)  + stepSpr(j)  * ks_rnge(j)
             
          elseif(i==2) then
             A0(j)       = A0_inpt(j)  + stepCell(j) * A0_rnge(j)
             k_area(j)   = ka_inpt(j)  + stepCell(j) * ka_rnge(j)
             
          elseif (i==3) then
             CgXNode(j)  = CgX_inpt(j) + stepNode(j) * CgX_rnge(j)
             CgYNode(j)  = CgY_inpt(j) + stepNode(j) * CgY_rnge(j)
          endif
          
          if (i==1) write(132,*) l0(j),     l0_inpt(j), l0_rnge(j),j, "l0Info spr_no"
          if (i==1) write(132,*) k_spr(j),  ks_inpt(j), ks_rnge(j),j, "ksInfo spr no"
          
          if (i==2) write(132,*) A0(j),     A0_inpt(j), A0_rnge(j),j, "A0Info cell_no"
          if (i==2) write(132,*) k_area(j), ka_inpt(j), ka_rnge(j),j, "kaInfo cell no"
          
          if (i==3) write(132,*) CgXNode(j),CgX_inpt(j),CgX_rnge(j),j,"CgXInfo node_no"
          if (i==3) write(132,*) CgYNode(j),CgY_inpt(j),CgY_rnge(j),j,"CgYInfo node_no"
          
       enddo
       
       write(132,*) " "
    enddo
    
    close(132)

    call set_diagSprPrp_back(diagsprL,diagsprR,l0_diagsprLstr,l0_diagsprRstr,&
         ks_diagsprLstr,ks_diagsprRstr)
    
  end subroutine stepUp_variables_excpt_DiagSpr
  
  
  subroutine Slowly_stepUp_variables_and_Equilibrate(Exprmnt,Frame,stepNo,&
       step,IftS,IftC,A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge)
    implicit none
    
    integer, intent(in) :: Exprmnt
    integer, intent(inout) :: Frame
    
    integer, intent(in) :: stepNo
    real*8,  intent(in) :: step
    real*8,  intent(in) :: A0_inpt(1:N_cell),A0_rnge(1:N_cell)
    real*8,  intent(in) :: ka_inpt(1:N_cell),ka_rnge(1:N_cell)
    real*8,  intent(in) :: l0_inpt(1:N_spr),l0_rnge(1:N_spr)
    real*8,  intent(in) :: ks_inpt(1:N_spr),ks_rnge(1:N_spr)
    integer, intent(in) :: IftS(1:N_spr),IftC(1:N_cell) 
    
    real*8  :: stepSprCur(1:N_spr),stepCellCur(1:N_cell)
    real*8  :: stepSprPrv(1:N_spr),stepCellPrv(1:N_cell)
    real*8  :: stepSprRng(1:N_spr),stepCellRng(1:N_cell)
    real*8  :: stepSpr(1:N_spr),stepCell(1:N_cell)
    
    integer :: i,j,imax,jmax,m
    
    real*8 :: PrvStep
    integer :: N_inStep
    
    PrvStep = 0.1d0*(stepNo-1)
    
    call get_step(PrvStep,stepNo,IftS,IftC,stepSprPrv,stepCellPrv)
    call get_step(step,stepNo,IftS,IftC,stepSprCur,stepCellCur)

    stepSprRng  = stepSprCur  - stepSprPrv
    stepCellRng = stepCellCur - stepCellPrv
    
    open(unit=130,file='SlowStep.dat',position='append')
    write(130,fmt=*) step,PrvStep,"step,PrvStep"

    write(130,fmt=*) "IftS","StepSprCur","StepSprPrv","StepSprRng","spr_nm"
    
    do i=1,N_spr
       write(130,fmt=*) IftS(i),stepSprCur(i),stepSprPrv(i),stepSprRng,i
    enddo
    
    write(130,fmt=*) "IftC","StepCellCur","StepCellPrv","StepCellRng","cell_nm" 

    do i=1,N_cell
       write(130,fmt=*) IftC(i),stepCellCur(i),stepCellPrv(i),stepCellRng(i),i
    enddo
    
    close(130)
    
    
    open(unit=131,file='InStepVars.dat',position='append')
    write(131,fmt=*) stepNo,"stepNo"
    
    imax = 2

    N_inStep = 5 
    
    do m = 1,N_inStep
       call get_inStpIncr(stepCellPrv,stepCellRng,N_cell,N_inStep,m,stepCell)
       call get_inStpIncr(stepSprPrv,stepSprRng,N_spr,N_inStep,m,stepSpr)
       
       do i = 1,imax
          
          if (i==1) jmax=N_cell
          if (i==2) jmax=N_spr
          
          do j = 1,jmax
             if(i==1) then
                A0(j) = A0_inpt(j) + stepCell(j)*A0_rnge(j)
                k_area(j) = ka_inpt(j) + stepCell(j)*ka_rnge(j)
             elseif(i==2) then
                l0(j) = l0_inpt(j) + stepSpr(j)*l0_rnge(j)
                k_spr(j) = ks_inpt(j) + stepSpr(j)*ks_rnge(j)
             endif
             
             if (i==1) write(131,fmt=*) A0(j),k_area(j),j,"A0,k_area,cell_no"
             if (i==2) write(131,fmt=*) l0(j),k_spr(j) ,j,"l0,k_spr, spr_no"
          enddo
          
          write(131,fmt=*) " "
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
       Frame = Frame+1
       call switchto_NI_model_run_and_switchbackto_TN
       
    enddo
 
    close(131)
    
    write(*,*) l0(37),l0(40),"l0 37 and 40"
    
  end subroutine Slowly_stepUp_variables_and_Equilibrate
  
  
  subroutine get_inpt_rnge(strct_INPT,strct_OTPT,A0_inpt,A0_otpt,A0_rnge,&
       ka_inpt,ka_otpt,ka_rnge,l0_inpt,l0_otpt,l0_rnge,ks_inpt,ks_otpt,ks_rnge)
    
    implicit none
    integer, intent(in)    :: strct_INPT,strct_OTPT
    
    real*8, intent(inout)  :: A0_inpt(1:N_cell),A0_otpt(1:N_cell),A0_rnge(1:N_cell)
    real*8, intent(inout)  :: l0_inpt(1:N_spr),l0_otpt(1:N_spr),l0_rnge(1:N_spr)
    real*8, intent(inout)  :: ka_inpt(1:N_cell),ka_otpt(1:N_spr),ka_rnge(1:N_cell)
    real*8, intent(inout)  :: ks_inpt(1:N_spr),ks_otpt(1:N_spr),ks_rnge(1:N_spr)
    
    integer :: i,j,imax,jmax
    integer :: cell_InptOtpt(1:N_cell)
    integer :: cell_nm
    integer :: spr_InptOtpt(1:N_spr)
    integer :: spr_nm
    
    
    cell_InptOtpt = 0
    
    if (stageNo==1 .and. stageType==1) then
       do i=1,N_cell
          cell_InptOtpt(i) = i
       enddo
       
    elseif (stageNo==4 .and. stageType==1) then
       if (modelID==1) then
          call get_cellInptOtpt(cell_InptOtpt)
       elseif (modelID==2) then
          write(*,*) "get_cellInptOtpt is not updated for modelID=2,stop1"
          stop
       endif
    elseif (stageNo==4 .and. stageType==2) then
       if (modelID==1) then
          call get_cellInptOtpt(cell_InptOtpt)
       elseif (modelID==2) then
          write(*,*) "get_cellInptOtpt is not updated for modelID=2,stop2"
          stop
       endif
    endif
    
    spr_InptOtpt  = 0 
    
    if (stageNo==1 .and. stageType==1) then
       do i=1,N_spr
          spr_InptOtpt(i) = i
       enddo
       
    elseif (stageNo==4 .and. stageType==1) then
       if (modelID==1) then
          call get_sprInptOtpt(spr_InptOtpt)
       elseif (modelID==2) then
          write(*,*) "get_sprInptOtpt is not updated for modelID=2,stop1"
          stop
       endif
    elseif (stageNo==4 .and. stageType==2) then
       if (modelID==1) then
          call get_sprInptOtpt(spr_InptOtpt)
       elseif (modelID==2) then
          write(*,*) "get_sprInptOtpt is not updated for modelID=2,stop2"
          stop
       endif
    endif
    
    imax = 2
    
    open(unit=301,file='input_output_l0A0.dat')
    
    do i = 1,imax
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       if(i==1) write(unit=301,fmt=*) "A0 inpt-otpt,ka inpt-otpt,cell no"
       if(i==2) write(unit=301,fmt=*) "l0 inpt-otpt,ks inpt-otpt,spr no"
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (modelID==1) then
                A0_inpt(j) = A0_strctTN(strct_INPT,j)
                ka_inpt(j) = ka_strctTN(strct_INPT,j)
                
             elseif (modelID==2) then
                A0_inpt(j) = A0_strctNI(strct_INPT,j)
                ka_inpt(j) = ka_strctNI(strct_INPT,j)
             endif
             
             cell_nm    = cell_InptOtpt(j)

             if (modelID==1) then
                A0_otpt(j) = A0_strctTN(strct_OTPT,cell_nm)
                ka_otpt(j) = ka_strctTN(strct_OTPT,cell_nm)
                
             elseif (modelID==2) then
                A0_otpt(j) = A0_strctNI(strct_OTPT,cell_nm)
                ka_otpt(j) = ka_strctNI(strct_OTPT,cell_nm)
             endif
             
             write(unit=301,fmt=*) A0_inpt(j),A0_otpt(j),ka_inpt(j),ka_otpt(j),j
             
          endif
          
          if (i==2) then
             
             if (modelID==1) then
                l0_inpt(j) = l0_strctTN(strct_INPT,j)
                ks_inpt(j) = ks_strctTN(strct_INPT,j)
                
             elseif (modelID==2) then
                l0_inpt(j) = l0_strctNI(strct_INPT,j)
                ks_inpt(j) = ks_strctNI(strct_INPT,j)
             endif
             
             spr_nm     = spr_InptOtpt(j)

             if (modelID==1) then
                l0_otpt(j) = l0_strctTN(strct_OTPT,spr_nm)
                ks_otpt(j) = ks_strctTN(strct_OTPT,spr_nm)
                
             elseif (modelID==2) then
                l0_otpt(j) = l0_strctNI(strct_OTPT,spr_nm)
                ks_otpt(j) = ks_strctNI(strct_OTPT,spr_nm)
             endif
             
             write(unit=301,fmt=*) l0_inpt(j),l0_otpt(j),ks_inpt(j),ks_otpt(j),j
             
          endif
          
       enddo

       write(unit=301,fmt=*) " "
       
    enddo
    
    
    do i = 1,imax
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             A0_rnge(j) = A0_otpt(j)-A0_inpt(j)
             ka_rnge(j) = ka_otpt(j)-ka_inpt(j)
          elseif (i==2) then
             l0_rnge(j) = l0_otpt(j)-l0_inpt(j)
             ks_rnge(j) = ks_otpt(j)-ks_inpt(j)
          endif
          
       enddo
    enddo
    
    
  end subroutine get_inpt_rnge
  
  subroutine get_Ift(IftS,IftC) !Ift=Interpolation Func Type
    implicit none
    integer, intent(inout) :: IftS(1:N_spr),IftC(1:N_cell)
    integer :: i,j,jmax
    integer :: diag_spr(1:2),apcl_spr(1:2)
    
    do i = 1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          if(i==1) IftS(j) = 1
          if(i==2) IftC(j) = 1
       enddo
       
    enddo
    
    !call get_the_DiagSpr(diag_spr)
    !IftS(diag_spr(1)) = 2 ; IftS(diag_spr(2)) = 2
    !call get_the_ApicalSpr_bfrPulley(apcl_spr)
    !IftS(apcl_spr(1)) = 3 ; IftS(apcl_spr(2)) = 3
    
  end subroutine get_Ift
  
  subroutine get_step(step,stepNo,IftS,IftC,stepSpr,stepCell)
    implicit none
    real*8, intent(in)    :: step
    integer,intent(in)    :: stepNo
    integer,intent(in)    :: IftS(1:N_spr),IftC(1:N_cell)
    real*8, intent(inout) :: stepSpr(1:N_spr),stepCell(1:N_cell)
    
    integer :: i,j,jmax
    real*8  :: hlfPi
    integer :: Ift_nm
    
    hlfPi= (2.0d0)*atan(1.0)
    
    do i = 1,2
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (IftS(j)==1) then
                stepSpr(j) = step
             elseif (IftS(j)==2) then
                stepSpr(j) = 1.0d0-cos(hlfPi*Step)
             elseif (IftS(j)==3) then
                Ift_nm = IftS(j)
                stepSpr(j) = func_StepSpr(stepNo,Ift_nm)
                
             elseif (IftS(j).ne.1 .and. IftS(j).ne.2 .and. IftS(j).ne.3) then
                write(*,*) "IftS should be 1,2 or 3 sbrtn : get_step()"
                stop
             endif
             
          elseif (i==2) then
             
             if (IftC(j)==1) then
                stepCell(j) = step
             elseif (IftC(j)==2) then
                stepCell(j) = 1.0d0-cos(hlfPi*Step)
             elseif (IftC(j).ne.1 .and. IftC(j).ne.2) then
                write(*,*) "IftC should be 1 or 2,sbrtn : get_step()"
                stop
             endif
             
          endif
          
       enddo
       
    enddo
    
  end subroutine get_step
  
  
  subroutine get_step_includingCg(step,stepNo,IftS,IftC,IftN,stepSpr,stepCell,stepNode)
    implicit none
    real*8, intent(in)    :: step
    integer,intent(in)    :: stepNo
    integer,intent(in)    :: IftS(1:N_spr),   IftC(1:N_cell),    IftN(1:N_node)
    real*8, intent(inout) :: stepSpr(1:N_spr),stepCell(1:N_cell),stepNode(1:N_node)
    
    integer :: i,j,jmax
    real*8  :: hlfPi
    integer :: Ift_nm
    
    hlfPi= (2.0d0)*atan(1.0)

    open(unit=131,file='stepIncludingCg.dat',position='append')
    write(131,*) stepNo,step,"stepFrctn"
    
    do i = 1,3
       
       if (i==1) then
          jmax=N_spr
          write(131,*) "Writing Spring Props at different Step"
       elseif (i==2) then
          jmax=N_cell
          write(131,*) "Writing Cell Props at different Step"
       elseif (i==3) then
          jmax=N_node
          write(131,*) "Writing Node Props at different Step"
       endif
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (IftS(j)==1) then
                stepSpr(j) = step
             elseif (IftS(j)==2) then
                stepSpr(j) = 1.0d0-cos(hlfPi*Step)
             elseif (IftS(j)==3) then
                Ift_nm = IftS(j)
                stepSpr(j) = func_StepSpr(stepNo,Ift_nm)
             elseif (IftS(j).ne.1 .and. IftS(j).ne.2 .and. IftS(j).ne.3) then
                stop 'IftS should be 1,2 or 3 sbrtn : get_step()'
             endif
             
             write(131,*) j,IftS(j),stepSpr(j),"IftS-stepSpr" 
             
          elseif (i==2) then
             
             if (IftC(j)==1) then
                stepCell(j) = step
             elseif (IftC(j)==2) then
                stepCell(j) = 1.0d0-cos(hlfPi*Step)
             else
                stop 'IftC should be 1 or 2,sbrtn : get_step()'
             endif
             
             write(131,*) j,IftC(j),stepCell(j),"IftC-stepCell" 
             
          elseif (i==3) then
             
             if (IftN(j)==1) then
                stepNode(j) = step
             else
                stop 'IftN should be 1,sbrtn : get_step()'
             endif
             
             write(131,*) j,IftN(j),stepNode(j),"IftN-stepNode" 
             
          endif
          
       enddo
       
       write(131,*) " "
    enddo

    close(131)
    
  end subroutine get_step_includingCg
  
  
  real*8 function func_StepSpr(stepNo,Ift_nm)
    implicit none
    integer, intent(in) :: stepNo
    integer, intent(in) :: Ift_nm
    
    real*8  :: stepSizeF,stepSizeS !F=Fast;S=Slow
    integer :: N_slowStp,N_fastStp
    real*8  :: slow_prtn,fast_prtn
    
    stepSizeS = 0.00d0
    N_slowStp = 6 
    slow_prtn = stepSizeS*N_slowStp
    
    N_fastStp = 10 - N_slowStp
    fast_prtn = 1.0d0 - slow_prtn 
    stepSizeF = fast_prtn/N_fastStp
    
    if (stepNo.le.N_slowStp) then
       func_StepSpr = stepSizeS*stepNo
    elseif (stepNo.gt.N_slowStp) then
       func_StepSpr = (stepNo-N_slowStp)*stepSizeF
    endif
    
    open(unit=133,file='Step_for_Ift3.dat',position='append')
    write(133,fmt=*) func_StepSpr,"func_StepSpr"
    close(133)
    
  end function func_StepSpr
  
  subroutine get_inStpIncr(stepElmPrv,stepElmRng,N_elmnt,N_inStep,currStep,stepElm)
    implicit none
    integer, intent(in)  :: N_elmnt 
    real*8 , intent(in)  :: stepElmPrv(1:N_elmnt)
    real*8 , intent(in)  :: stepElmRng(1:N_elmnt)
    integer, intent(in)  :: N_inStep,currStep
    real*8 , intent(out) :: stepElm(1:N_elmnt)

    real*8  :: stepSizeV
    integer :: i
    
    stepSizeV = (1.0d0/N_inStep)*currStep

    open(unit=186,file='InstpIncr.dat',position='append')

    write(186,fmt=*) "stepElm","Elm No"
    
    do i = 1,N_elmnt  
       stepElm(i) = stepElmPrv(i) + stepSizeV*stepElmRng(i) 
       write(186,fmt=*) stepElm(i),i
    enddo
    
    close(186)
    
  end subroutine get_inStpIncr
  
  subroutine get_rdfnSprs1(sprs_Rdfn)
    implicit none
    integer, intent(inout) :: sprs_Rdfn(1:N_spr)

    integer :: i
    
    open(unit=305,file='sprsRdfn.dat')

    do i = 1,N_spr
       
       if (i.le.9) then
          sprs_Rdfn(i) = i
          
       elseif (i.gt.9 .and. i.le.30) then
          
          if (i==10) sprs_Rdfn(i) = 13
          if (i==11) sprs_Rdfn(i) = 14
          if (i==12) sprs_Rdfn(i) = 15
          if (i==13) sprs_Rdfn(i) = 16
          if (i==14) sprs_Rdfn(i) = 17
          if (i==15) sprs_Rdfn(i) = 18
          if (i==16) sprs_Rdfn(i) = 19
          if (i==17) sprs_Rdfn(i) = 20
          if (i==18) sprs_Rdfn(i) = 21
          
          if (i==19) sprs_Rdfn(i) = 10
          if (i==20) sprs_Rdfn(i) = 12
          if (i==21) sprs_Rdfn(i) = 11
          
          if (i==22) sprs_Rdfn(i) = 22
          if (i==23) sprs_Rdfn(i) = 24
          if (i==24) sprs_Rdfn(i) = 23
          
          if (i==25) sprs_Rdfn(i) = 25
          if (i==26) sprs_Rdfn(i) = 27
          if (i==27) sprs_Rdfn(i) = 26
          
          if (i==28) sprs_Rdfn(i) = 28
          if (i==29) sprs_Rdfn(i) = 30
          if (i==30) sprs_Rdfn(i) = 29
          

       elseif (i.gt.30) then
          sprs_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd spr is =",sprs_Rdfn(i), "for spr old =",i
       
    enddo

    close(305)
    
  end subroutine get_rdfnSprs1
  
  subroutine get_rdfnSprs(sprs_Rdfn)
    implicit none
    integer, intent(inout) :: sprs_Rdfn(1:N_spr)
    
    integer :: i    
    integer :: PrevNsl,PrevNsr
    integer :: CurrNsl,CurrNsr
    integer :: PrevNsCntrlTop
    integer :: cslr
    
    open(unit=305,file='sprsRdfnRwrt.dat')
    
    CurrNsl = nsl ; PrevNsl = nsl+nsecl*1
    CurrNsr = nsr ; PrevNsr = nsr+nsecr*1 
    
    PrevNsCntrlTop = PrevNsl + PrevNsr + nsclt + nscrt
    cslr = CurrNsl+CurrNsr
    
    write(305,fmt=*) CurrNsl,CurrNsr,PrevNsCntrlTop
    write(305,fmt=*) PrevNsl,PrevNsr
    write(305,fmt=*) nsecl
    
    do i = 1,N_spr
       
       if (i.le.CurrNsl) then
          sprs_Rdfn(i) = i
          
       elseif (i.gt.CurrNsl .and. i.le.cslr) then
          sprs_Rdfn(i) = i+(nsecl*1) !1 dueto 1 cell chopped off
          
       elseif (i.gt.cslr .and. i.le.PrevNsCntrlTop) then          
          if (i==(cslr+1)) sprs_Rdfn(i) = PrevNsl-2
          if (i==(cslr+2)) sprs_Rdfn(i) = PrevNsl
          if (i==(cslr+3)) sprs_Rdfn(i) = PrevNsl-1
          
          if (i==(cslr+4)) sprs_Rdfn(i) = PrevNsl+PrevNsr-2
          if (i==(cslr+5)) sprs_Rdfn(i) = PrevNsl+PrevNsr
          if (i==(cslr+6)) sprs_Rdfn(i) = PrevNsl+PrevNsr-1
          
          if (i==(cslr+7)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt-2
          if (i==(cslr+8)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt
          if (i==(cslr+9)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt-1
          
          if (i==(cslr+10)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt+nscrt-2
          if (i==(cslr+11)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt+nscrt
          if (i==(cslr+12)) sprs_Rdfn(i) = PrevNsl+PrevNsr+nsclt+nscrt-1
          
       elseif (i.gt.PrevNsCntrlTop) then
          sprs_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd spr is =",sprs_Rdfn(i), "for spr old =",i
       
    enddo
    
    close(305)
    
  end subroutine get_rdfnSprs

  subroutine get_rdfnAreas1(areas_Rdfn)
    implicit none
    integer, intent(inout) :: areas_Rdfn(1:N_cell)
    
    integer :: i
    
    open(unit=305,file='areasRdfn.dat')

    do i = 1,N_cell
       
       if (i.le.3) then
          areas_Rdfn(i) = i
          
       elseif (i.gt.3 .and. i.le.7) then
          
          if (i==4)  areas_Rdfn(i) = 5
          if (i==5)  areas_Rdfn(i) = 6
          if (i==6)  areas_Rdfn(i) = 7      
          if (i==7)  areas_Rdfn(i) = 4
     
       elseif (i.gt.7) then
          areas_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd area is =",areas_Rdfn(i), "for area old =",i
       
    enddo

    close(305)

  end subroutine get_rdfnAreas1

  subroutine get_rdfnAreas(areas_Rdfn)
    implicit none
    integer, intent(inout) :: areas_Rdfn(1:N_cell)

    integer :: i
    integer :: PrevNcl,PrevNcr
    integer :: CurrNcl,CurrNcr
    integer :: PrevNcCntrlTop
    integer :: cclr
    
    open(unit=305,file='areasRdfnRwrt.dat')
    
    CurrNcl = ncl ; PrevNcl = ncl+1
    CurrNcr = ncr ; PrevNcr = ncr+1

    PrevNcCntrlTop = PrevNcl+PrevNcr+ncclt+nccrt
    cclr = CurrNcl + CurrNcr
    
    do i = 1,N_cell
       
       if (i.le.CurrNcl) then
          areas_Rdfn(i) = i
       elseif (i.gt.CurrNcl .and. i.le.cclr) then
          areas_Rdfn(i) = i+1 !1 dueto 1 cell chopped off
          
       elseif (i.gt.cclr.and.i.le.PrevNcCntrlTop) then
          if (i==(cclr+1)) areas_Rdfn(i) = PrevNcl
          if (i==(cclr+2)) areas_Rdfn(i) = PrevNcl+PrevNcr
          if (i==(PrevNcCntrlTop-1)) areas_Rdfn(i) = i
          if (i==(PrevNcCntrlTop-0)) areas_Rdfn(i) = i
          
       elseif (i.gt.PrevNcCntrlTop) then
          areas_Rdfn(i) = i
       endif
       
       write(unit=305,fmt=*)"Rdfnd area is =",areas_Rdfn(i), "for area old =",i
       
    enddo

    close(305)
    
  end subroutine get_rdfnAreas
  
  
  subroutine get_cellInptOtpt(cell_InptOtpt)
    implicit none
    integer,intent(inout) :: cell_InptOtpt(1:N_cell)

    integer :: nclPrv,ncrPrv
    
    open(unit=302,file='cellInptOtpt.dat')

    nclPrv = ncl-1 ; ncrPrv = ncr-1
    
    do i = 1,N_cell
       
       if (inpt_rngeFlag == 1) then
          
          if (i.le.(ncl-1)) then
             cell_InptOtpt(i) = i
             
          elseif(i.gt.(ncl-1) .and. i.le.(ncl+ncr-1)) then
             
             if (i==(ncl))   cell_InptOtpt(i) = (nclPrv)+(ncrPrv)+1
             if (i==(ncl+1)) cell_InptOtpt(i) = (nclPrv)+1
             if (i==(ncl+2)) cell_InptOtpt(i) = (nclPrv)+2
             if (i==(ncl+3)) cell_InptOtpt(i) = (nclPrv)+3
             
          elseif (i.gt.(ncl+ncr-1)) then
             cell_InptOtpt(i) = i
          endif
          
          write(unit=302,fmt=*) "output is =",cell_InptOtpt(i), "for cell input =",i

       elseif (inpt_rngeFlag == 2) then
          cell_InptOtpt(i) = i
       else
          write(*,*) "Faulty inpt_rngeFlag value,sb:get_cellInptOtpt,fl:build_s"
       endif
       
    enddo
    
    close(302)
    
  end subroutine get_cellInptOtpt
  

  subroutine get_sprInptOtpt(spr_InptOtpt)
    
    implicit none
    integer, intent(inout) :: spr_InptOtpt(1:N_spr)

    integer :: i 

    open(unit=303,file='sprInptOtpt.dat')

       
    do i = 1,N_spr
       if (inpt_rngeFlag == 1) then !inpt_rngeFlag is the key that will let us escape the non-generality of the routine with nums
          
          if (i.le.9) then
             spr_InptOtpt(i) = i
             
          elseif(i.gt.9 .and. i.le.30) then
             
             if (i==10) spr_InptOtpt(i) = 19
             if (i==11) spr_InptOtpt(i) = 21
             if (i==12) spr_InptOtpt(i) = 20
             if (i==13) spr_InptOtpt(i) = 10
             if (i==14) spr_InptOtpt(i) = 11
             if (i==15) spr_InptOtpt(i) = 12
             if (i==16) spr_InptOtpt(i) = 13
             if (i==17) spr_InptOtpt(i) = 14
             if (i==18) spr_InptOtpt(i) = 15
             if (i==19) spr_InptOtpt(i) = 16
             if (i==20) spr_InptOtpt(i) = 17
             if (i==21) spr_InptOtpt(i) = 18
             if (i==22) spr_InptOtpt(i) = 22
             if (i==23) spr_InptOtpt(i) = 24
             if (i==24) spr_InptOtpt(i) = 23
             if (i==25) spr_InptOtpt(i) = 25
             if (i==26) spr_InptOtpt(i) = 27
             if (i==27) spr_InptOtpt(i) = 26
             if (i==28) spr_InptOtpt(i) = 28
             if (i==29) spr_InptOtpt(i) = 30
             if (i==30) spr_InptOtpt(i) = 29
             
          elseif (i.gt.30) then
             spr_InptOtpt(i) = i
          endif
          
          write(unit=303,fmt=*) "output is =",spr_InptOtpt(i), "for spr input =",i
          
       elseif (inpt_rngeFlag == 2) then
          spr_InptOtpt(i) = i
       else
          write(*,*)"Faulty inpt_rngeFlag value,sb:get_sprInptOtpt,fl:build_s"
       endif
       
    enddo

    close(303)

  end subroutine get_sprInptOtpt
  
  
  subroutine redefine_Prop_strct3
    implicit none
    real*8 :: A0_Prv(1:N_cell),ka_Prv(1:N_cell)
    real*8 :: l0_Prv(1:N_spr) ,ks_Prv(1:N_spr)
    
    integer :: sprs_Rdfn(1:N_spr),areas_Rdfn(1:N_cell) 
    integer :: spr_nm,area_nm
    
    integer :: N_item
    integer :: i,j,jmax

    if (modelID==1) then
       A0_Prv(1:N_cell) = A0_strctTN(3,1:N_cell)   
       ka_Prv(1:N_cell) = ka_strctTN(3,1:N_cell)
       
       l0_Prv(1:N_spr)  = l0_strctTN(3,1:N_spr)
       ks_Prv(1:N_spr)  = ks_strctTN(3,1:N_spr) 
       
    elseif (modelID==2) then
       A0_Prv(1:N_cell) = A0_strctNI(3,1:N_cell)   
       ka_Prv(1:N_cell) = ka_strctNI(3,1:N_cell)
       
       l0_Prv(1:N_spr)  = l0_strctNI(3,1:N_spr)
       ks_Prv(1:N_spr)  = ks_strctNI(3,1:N_spr)
    endif
    
    sprs_Rdfn  = 0
    areas_Rdfn = 0
    
    open(unit=302,file='l0A0_bfrRdfn.dat')
    
    write(302,fmt=*) "l0_strct","ks_strct","spr"
    
    do i = 1,N_spr
       if (modelID==1) write(unit=302,fmt=*) l0_strctTN(3,i),ks_strctTN(3,i),i
       if (modelID==2) write(unit=302,fmt=*) l0_strctNI(3,i),ks_strctNI(3,i),i
    enddo
    
    write(302,fmt=*) " "
    
    write(302,fmt=*) "A0_strct","ka_strct","area"
    
    do i = 1,N_cell
       if (modelID==1) write(unit=302,fmt=*) A0_strctTN(3,i),ka_strctTN(3,i),i
       if (modelID==2) write(unit=302,fmt=*) A0_strctNI(3,i),ka_strctNI(3,i),i
    enddo
    
    close(302)

    if (modelID==1) then
       call get_rdfnSprs(sprs_Rdfn)
       call get_rdfnAreas(areas_Rdfn)
    elseif (modelID==2) then
       write(*,*) "get_rdfnSprs and get_rdfnAreas are not ready fpr modelID=2"
       write(*,*) "stop is in sbrtn : redefine_Prop_strct3, in flnm: build_strctures"
       stop
    endif
       
    N_item = 2
    
    
    do i = 1,N_item
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             area_nm = areas_Rdfn(j)
             
             if (modelID==1) then
                A0_strctTN(3,j) = A0_Prv(area_nm)
                ka_strctTN(3,j) = ka_Prv(area_nm)
                
             elseif (modelID==2) then
                write(*,*) "get_rdfnAreas not ready"
                stop
             endif
                
          elseif (i==2) then
             spr_nm = sprs_Rdfn(j)
             
             if (modelID==1) then
                l0_strctTN(3,j) = l0_Prv(spr_nm)
                ks_strctTN(3,j) = ks_Prv(spr_nm)
                
             elseif (modelID==2) then
                write(*,*) "get_rdfnSprs not ready"
                stop
             endif
             
          endif
          
       enddo
       
    enddo
    
    open(unit=302,file='l0A0_aftRdfn.dat')
    
    write(302,fmt=*) "l0_strct","ks_strct","spr"
    
    do i = 1,N_spr
       if (modelID==1) write(unit=302,fmt=*) l0_strctTN(3,i),ks_strctTN(3,i),i
       if (modelID==2) stop
    enddo
    
    write(unit=302,fmt=*) " "
    
    write(302,fmt=*) "A0_strct","ka_strct","area"
    
    do i = 1,N_cell
       if (modelID==1) write(unit=302,fmt=*) A0_strctTN(3,i),ka_strctTN(3,i),i
       if (modelID==2) stop
    enddo
    
    close(302)
    
  end subroutine redefine_Prop_strct3
  
  subroutine redefine_Props
    implicit none
    real*8 :: A0_Prv(1:N_cell),ka_Prv(1:N_cell)
    real*8 :: A_Prv(1:N_cell),At_Prv(1:N_cell),fctr_Prv(1:N_cell)
    real*8 :: beta_Prv(1:N_cell),optmCell_Prv(1:N_Cell)
    
    real*8 :: l0_Prv(1:N_spr), ks_Prv(1:N_spr)
    real*8 :: l_Prv(1:N_spr) , Lt_Prv(1:N_spr)
    real*8 :: alpha_Prv(1:N_spr),optmSpr_Prv(1:N_spr)
    
    integer :: sprs_Rdfn(1:N_spr),areas_Rdfn(1:N_cell) 
    integer :: spr_nm,area_nm
    
    integer :: N_item
    integer :: i,j,jmax
     
    
    A_Prv  = A
    At_Prv = At
    beta_Prv = beta
    optmCell_Prv = optmCell
    fctr_Prv = fctr
    
    l_Prv  = l
    Lt_Prv = Lt
    alpha_Prv   = alpha 
    optmSpr_Prv = optmSpr
    
    sprs_Rdfn  = 0
    areas_Rdfn = 0
    
    open(unit=303,file='Prprty_bfr.dat')
    
    write(303,fmt=*) "l0_bfr","kspr_bfr","spr"
    
    do i = 1,N_spr
       write(unit=303,fmt=*) l0(i),k_spr(i),i
    enddo
    
    write(unit=303,fmt=*) " "
    
    write(303,fmt=*) "A0_bfr","karea_bfr","area"
    
    do i = 1,N_cell
       write(unit=303,fmt=*) A0(i),k_area(i),i
    enddo
    
    close(303)
    
    call get_rdfnSprs(sprs_Rdfn)
    call get_rdfnAreas(areas_Rdfn)
    
    N_item = 2
    
    
    do i = 1,N_item
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          if (i==1) then
             area_nm     = areas_Rdfn(j)
             A(j)        = A_Prv(area_nm)
             At(j)       = At_Prv(area_nm)
             beta(j)     = beta_Prv(area_nm)
             optmCell(j) = optmCell_Prv(area_nm)
             fctr(j)     = fctr_Prv(area_nm)
             
          elseif (i==2) then
             spr_nm     = sprs_Rdfn(j)
             l(j)       = l_Prv(spr_nm)
             Lt(j)      = Lt_Prv(spr_nm)
             alpha(j)   = alpha_Prv(spr_nm)
             optmSpr(j) = optmSpr_Prv(spr_nm)
             
          endif
          
       enddo
       
    enddo
    
    open(unit=303,file='Prprty_aft.dat')
    
    write(303,fmt=*) "l0_aft","kspr_aft","spr"
    
    do i = 1,N_spr
       write(unit=303,fmt=*) l0(i),k_spr(i),i
    enddo
    
    write(303,fmt=*) " "
    write(303,fmt=*) "A0_bfr","karea_bfr","area"
    
    do i = 1,N_cell
       write(unit=303,fmt=*) A0(i),k_area(i),i
    enddo
    
    close(303)
    
    
  end subroutine redefine_Props
  
  subroutine updt_trgt(A0_inpt,A0_rnge,ka_inpt,ka_rnge,l0_inpt,l0_rnge,ks_inpt,ks_rnge,A0_trgt,ka_trgt,l0_trgt,ks_trgt)
    
    implicit none
    real*8, intent(in) :: A0_inpt(1:N_cell),A0_rnge(1:N_cell)
    real*8, intent(in) :: ka_inpt(1:N_cell),ka_rnge(1:N_cell)
    real*8, intent(in) :: l0_inpt(1:N_spr),l0_rnge(1:N_spr)
    real*8, intent(in) :: ks_inpt(1:N_spr),ks_rnge(1:N_spr)

    real*8, intent(out) :: A0_trgt(1:N_cell),ka_trgt(1:N_spr)
    real*8, intent(out) :: l0_trgt(1:N_spr) ,ks_trgt(1:N_spr)

    integer :: i,imax,j,jmax

    open(unit=137,file='Updtd_Trgt.dat')
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       if (i==1) write(137,fmt=*) "A0_trgt","","ka_trgt","","area"
       if (i==2) write(137,fmt=*) "l0_trgt","","ks_trgt","","spr"
       
       do j = 1,jmax
          
          if(i==1) then
             A0_trgt(j) = A0_inpt(j) + A0_rnge(j)
             ka_trgt(j) = ka_inpt(j) + ka_rnge(j)
          elseif(i==2) then
             l0_trgt(j) = l0_inpt(j) + l0_rnge(j)
             ks_trgt(j) = ks_inpt(j) + ks_rnge(j)
          endif
          
          if (i==1) write(137,fmt=*) A0_trgt(j),ka_trgt(j),j
          if (i==2) write(137,fmt=*) l0_trgt(j),ks_trgt(j),j
          
       enddo
       
       write(137,fmt=*) " "
       
    enddo

    close(137)

  end subroutine updt_trgt
  
  
  subroutine Interpltn_adjstmnt_for_spr(l0_inpt,Nstp,currStp,l0_trgt,l0_rnge)
    implicit none
    real*8, intent(in)    :: l0_trgt(1:N_spr)
    integer,intent(in)    :: Nstp,currStp
    real*8, intent(inout) :: l0_inpt(1:N_spr),l0_rnge(1:N_spr)
    
    real*8  :: Pf,Pc,Pi,Prnge !Pf=Property's final value
    integer :: i
    
    
    do i = 1,N_spr
       if (i==12) then
          open(unit=177,file='IntrpltnChk.dat',position='append')
       endif
       
       Pf = l0_trgt(i)
       Pc = l0(i)
       Pi = l0_inpt(i)
       Prnge = l0_rnge(i)
       
       if (i==12) write(177,fmt=*) Pf,Pc,Pi,Prnge,"bfr updt"
       
       call get_the_updtd_initProp_and_rnge(Pf,Pc,Nstp,currStp,Pi,Prnge)
       
       if (i==12) write(177,fmt=*) Pf,Pc,Pi,Prnge,"aft updt"
       
       l0_inpt(i) = Pi 
       l0_rnge(i) = Prnge
       
       if (i==12) then
          close(177)
       endif
       
    enddo
    
  end subroutine Interpltn_adjstmnt_for_spr

  
  subroutine Interpltn_adjstmnt_for_cell(A0_inpt,Nstp,currStp,A0_trgt,A0_rnge)
    implicit none
    real*8, intent(in)    :: A0_trgt(1:N_spr)
    integer,intent(in)    :: Nstp,currStp
    real*8, intent(inout) :: A0_inpt(1:N_spr),A0_rnge(1:N_spr)
    
    real*8  :: Pf,Pc,Pi,Prnge !Pf=Property's final value
    integer :: i
    
    do i = 1,N_cell
       if (i==25) then
          open(unit=177,file='IntrpltnChk.dat')
       endif
       
       Pf = A0_trgt(i)
       Pc = A0(i)
       Pi = A0_inpt(i)
       Prnge = A0_rnge(i)
       
       if (i==25) write(177,fmt=*) Pf,Pc,Pi,Prnge,"bfr updt"
       
       call get_the_updtd_initProp_and_rnge(Pf,Pc,Nstp,currStp,Pi,Prnge)
       
       if (i==25) write(177,fmt=*) Pf,Pc,Pi,Prnge,"aft updt"
       
       A0_inpt(i) = Pi 
       A0_rnge(i) = Prnge
       
       if (i==25) then
          close(177)
       endif
       
    enddo
    
  end subroutine Interpltn_adjstmnt_for_cell


  
  subroutine get_the_updtd_initProp_and_rnge(Pf,Pc,Nstp,currStp,Pi,Prnge)
    implicit none
    real*8, intent(in)  :: Pf,Pc
    integer,intent(in)  :: Nstp,currStp
    real*8, intent(inout) :: Pi,Prnge
    
    real*8  :: ratio,NstpReal,currStpReal
    
    write(177,fmt=*) " "
    write(177,fmt=*) Pf,Pc,Nstp,currStp,"in"
    
    NstpReal=real(Nstp) ; currStpReal=real(currStp)
    write(177,fmt=*) NstpReal,currStpReal,"NstpReal,currStpReal"
    
    ratio = (NstpReal/(NstpReal-currStpReal))
    write(177,fmt=*) ratio,"ratio"
    
    Pi = Pf - ratio * (Pf-Pc)
    write(177,*) Pi
    Prnge = Pf - Pi
    
    write(177,fmt=*) Pi,Prnge,"out"
    write(177,fmt=*) " "
    
  end subroutine get_the_updtd_initProp_and_rnge
  
  
  subroutine print_the_system
    implicit none
    
    call print_NodeVars
    call print_spring_variables
    call print_area_variables
    call print_grvtnl_vars
    call print_bend_variables
    call print_all_trnsfrms
    call print_coordnteVars
    
  end subroutine print_the_system
  
  subroutine read_strctProps(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    flnm='strctProps'
    N_itm = 5
    
    if (stageNo==1 .and. stageType==1) then
       
       if (strctNum.le.9) then
          write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
       elseif (strctNum.gt.9 .and. strctNum.le.99) then
          write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
       endif
    
    elseif (stageNo==4 .and. stageType==1) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T1.dat'
    elseif (stageNo==4 .and. stageType==2) then
       write(flnmbr,'(i1.1,a)') strctNum,'S4T2.dat'
    endif
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=21,file=trim(adjustl(full_flnm)))
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       if (strctNum==1 .and. i==4) then
          N_mvCoordnte1 = N_mvCoordnte
          jmax=N_mvCoordnte1
          allocate(coordntesXY_strct1(1:N_mvCoordnte1))
          
          NMCWL0A01 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct1(1:NMCWL0A01))
          
       elseif (strctNum==2 .and. i==4) then
          N_mvCoordnte2 = N_mvCoordnte
          jmax=N_mvCoordnte2
          allocate(coordntesXY_strct2(1:N_mvCoordnte2))
          
          NMCWL0A02 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct2(1:NMCWL0A02))
          
       elseif (strctNum==3 .and. i==4) then
          
          if (stageNo==4) then
             continue
          elseif (stageNo==1 .and. stageType==1) then
             N_mvCoordnte3 = N_mvCoordnte
             jmax=N_mvCoordnte3
             allocate(coordntesXY_strct3(1:N_mvCoordnte3))
             
             NMCWL0A03 = N_mvCoordnte_withl0A0
             allocate(coordntes_strct3(1:NMCWL0A03))
          endif
          
       elseif (strctNum==4 .and. i==4) then
          N_mvCoordnte4 = N_mvCoordnte
          jmax=N_mvCoordnte4
          allocate(coordntesXY_strct4(1:N_mvCoordnte4))
          
          NMCWL0A04 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct4(1:NMCWL0A04))
          
       elseif (strctNum==5 .and. i==4) then
          N_mvCoordnte5 = N_mvCoordnte
          jmax=N_mvCoordnte5
          allocate(coordntesXY_strct5(1:N_mvCoordnte5))
          
          NMCWL0A05 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct5(1:NMCWL0A05))
          
       elseif (strctNum==6 .and. i==4) then
          N_mvCoordnte6 = N_mvCoordnte
          jmax=N_mvCoordnte6
          allocate(coordntesXY_strct6(1:N_mvCoordnte6))
          
          NMCWL0A06 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct6(1:NMCWL0A06))
          
       elseif (strctNum==7 .and. i==4) then
          N_mvCoordnte7 = N_mvCoordnte
          jmax=N_mvCoordnte7
          allocate(coordntesXY_strct7(1:N_mvCoordnte7))
          
          NMCWL0A07 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct7(1:NMCWL0A07))
          
       elseif (strctNum==8 .and. i==4) then
          N_mvCoordnte8 = N_mvCoordnte
          jmax=N_mvCoordnte8
          allocate(coordntesXY_strct8(1:N_mvCoordnte8))
          
          NMCWL0A08 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct8(1:NMCWL0A08))
          
       elseif (strctNum==9 .and. i==4) then
          N_mvCoordnte9 = N_mvCoordnte
          jmax=N_mvCoordnte9
          allocate(coordntesXY_strct9(1:N_mvCoordnte9))
          
          NMCWL0A09 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct9(1:NMCWL0A09))
          
       elseif (strctNum==10 .and. i==4) then
          N_mvCoordnte10 = N_mvCoordnte
          jmax=N_mvCoordnte10
          allocate(coordntesXY_strct10(1:N_mvCoordnte10))
          
          NMCWL0A010 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct10(1:NMCWL0A010))
          
       elseif (strctNum==11 .and. i==4) then
          N_mvCoordnte11 = N_mvCoordnte
          jmax=N_mvCoordnte11
          allocate(coordntesXY_strct11(1:N_mvCoordnte11))
          
          NMCWL0A011 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct11(1:NMCWL0A011))
          
       elseif (strctNum==12 .and. i==4) then
          N_mvCoordnte12 = N_mvCoordnte
          jmax=N_mvCoordnte12
          allocate(coordntesXY_strct12(1:N_mvCoordnte12))
          
          NMCWL0A012 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct12(1:NMCWL0A012))
          
       elseif (strctNum==13 .and. i==4) then
          N_mvCoordnte13 = N_mvCoordnte
          jmax=N_mvCoordnte13
          allocate(coordntesXY_strct13(1:N_mvCoordnte13))
          
          NMCWL0A013 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct13(1:NMCWL0A013))
          
       elseif (strctNum==14 .and. i==4) then
          N_mvCoordnte14 = N_mvCoordnte
          jmax=N_mvCoordnte14
          allocate(coordntesXY_strct14(1:N_mvCoordnte14))
          
          NMCWL0A014 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct14(1:NMCWL0A014))
          
       elseif (strctNum==15 .and. i==4) then
          N_mvCoordnte15 = N_mvCoordnte
          jmax=N_mvCoordnte15
          allocate(coordntesXY_strct15(1:N_mvCoordnte15))
          
          NMCWL0A015 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct15(1:NMCWL0A015))
          
       elseif (strctNum==16 .and. i==4) then
          N_mvCoordnte16 = N_mvCoordnte
          jmax=N_mvCoordnte16
          allocate(coordntesXY_strct16(1:N_mvCoordnte16))
          
          NMCWL0A016 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct16(1:NMCWL0A016))
          
       elseif (strctNum==17 .and. i==4) then
          N_mvCoordnte17 = N_mvCoordnte
          jmax=N_mvCoordnte17
          allocate(coordntesXY_strct17(1:N_mvCoordnte17))
          
          NMCWL0A017 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct17(1:NMCWL0A017))
          
       elseif (strctNum==18 .and. i==4) then
          N_mvCoordnte18 = N_mvCoordnte
          jmax=N_mvCoordnte18
          allocate(coordntesXY_strct18(1:N_mvCoordnte18))
          
          NMCWL0A018 = N_mvCoordnte_withl0A0
          allocate(coordntes_strct18(1:NMCWL0A018))
          
       endif
       
       if (i==5) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) read(21,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j)
             if (i==2) read(21,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j)
             if (i==3) read(21,*) CgX_strctTN(strctNum,j)
             
          elseif (modelID==2) then
             write(*,*) "Not ready due to coordntesXY_strct1,2...18 arrays"
             write(*,*) "sb: read_strctProps inside fl: build_strctures"
             stop
          endif
          
          if (strctNum==1 .and. i==4) then
             read(21,*) coordntesXY_strct1(j)
          elseif (strctNum==2 .and. i==4) then
             read(21,*) coordntesXY_strct2(j)
          elseif (strctNum==3 .and. i==4) then
             
             if (stageNo==4) then
                continue
             elseif (stageNo==1 .and. stageType==1) then
                read(21,*) coordntesXY_strct3(j)
             endif
             
          elseif (strctNum==4 .and. i==4) then
             read(21,*) coordntesXY_strct4(j)
          elseif (strctNum==5 .and. i==4) then
             read(21,*) coordntesXY_strct5(j)
          elseif (strctNum==6 .and. i==4) then
             read(21,*) coordntesXY_strct6(j)
          elseif (strctNum==7 .and. i==4) then
             read(21,*) coordntesXY_strct7(j)
          elseif (strctNum==8 .and. i==4) then
             read(21,*) coordntesXY_strct8(j)
          elseif (strctNum==9 .and. i==4) then
             read(21,*) coordntesXY_strct9(j)
          elseif (strctNum==10 .and. i==4) then
             read(21,*) coordntesXY_strct10(j) 
          elseif (strctNum==11 .and. i==4) then
             read(21,*) coordntesXY_strct11(j)
          elseif (strctNum==12 .and. i==4) then
             read(21,*) coordntesXY_strct12(j)
          elseif (strctNum==13 .and. i==4) then
             read(21,*) coordntesXY_strct13(j)
          elseif (strctNum==14 .and. i==4) then
             read(21,*) coordntesXY_strct14(j)
          elseif (strctNum==15 .and. i==4) then
             read(21,*) coordntesXY_strct15(j)
          elseif (strctNum==16 .and. i==4) then
             read(21,*) coordntesXY_strct16(j)
          elseif (strctNum==17 .and. i==4) then
             read(21,*) coordntesXY_strct17(j)
          elseif (strctNum==18 .and. i==4) then
             read(21,*) coordntesXY_strct18(j)
             
          endif

          if (modelID==1) then
             if (i==5) read(21,*)nodeXY_strctTN(strctNum,j,1:N_dmnsn)
          elseif (modelID==2) then
             write(*,*) "Not ready due to coordntesXY_strct1,2...18 arrays"
             write(*,*) "sb: read_strctProps inside fl: build_strctures"
             stop
          endif
       enddo
       
    enddo
    
    close(21)
    
  end subroutine read_strctProps
  
  
  subroutine Strct_to_actual(struct_No)
    implicit none
    integer, intent(in) :: struct_No
    integer :: i,j
    
    if (modelID == 1) then
       l0(1:N_spr)      = l0_strctTN(struct_No,1:N_spr)
       k_spr(1:N_spr)   = ks_strctTN(struct_No,1:N_spr)
       A0(1:N_cell)     = A0_strctTN(struct_No,1:N_cell)
       k_area(1:N_cell) = ka_strctTN(struct_No,1:N_cell)
       CgXNode(1:N_node) = CgX_strctTN(struct_No,1:N_node)
       
    elseif (modelID==2) then
       l0(1:N_spr)      = l0_strctNI(struct_No,1:N_spr)
       k_spr(1:N_spr)   = ks_strctNI(struct_No,1:N_spr)
       A0(1:N_cell)     = A0_strctNI(struct_No,1:N_cell)
       k_area(1:N_cell) = ka_strctNI(struct_No,1:N_cell)
       CgXNode(1:N_node) = CgX_strctNI(struct_No,1:N_node)
       
    endif
    
    if (modelID==1) then
       
       if (struct_No==1) then
          coordntes_xy = coordntesXY_strct1
       elseif (struct_No==2) then
          coordntes_xy = coordntesXY_strct2
          
       elseif (struct_No==3) then
          if (stageNo==4) then
             continue
          elseif (stageNo==1 .and. stageType==1) then
             coordntes_xy = coordntesXY_strct3
          endif
          
       elseif (struct_No==4) then
          coordntes_xy = coordntesXY_strct4
       elseif (struct_No==5) then
          coordntes_xy = coordntesXY_strct5
       elseif (struct_No==6) then
          coordntes_xy = coordntesXY_strct6
       elseif (struct_No==7) then
          coordntes_xy = coordntesXY_strct7
       elseif (struct_No==8) then
          coordntes_xy = coordntesXY_strct8
       elseif (struct_No==9) then
          coordntes_xy = coordntesXY_strct9
       elseif (struct_No==10) then
          coordntes_xy = coordntesXY_strct10
       elseif (struct_No==11) then
          coordntes_xy = coordntesXY_strct11
       elseif (struct_No==12) then
          coordntes_xy = coordntesXY_strct12
       elseif (struct_No==13) then
          coordntes_xy = coordntesXY_strct13
       elseif (struct_No==14) then
          coordntes_xy = coordntesXY_strct14
       elseif (struct_No==15) then
          coordntes_xy = coordntesXY_strct15
       elseif (struct_No==16) then
          coordntes_xy = coordntesXY_strct16
       elseif (struct_No==17) then
          coordntes_xy = coordntesXY_strct17
       elseif (struct_No==18) then
          coordntes_xy = coordntesXY_strct18
       endif
       
    elseif(modelID==2) then
       write(*,*) "coordntesXY_strct1,2,...18 are not defined for NI model"
       write(*,*) "sb: Strct_to_actual, fl:build_strctures"
       write(*,*) "But its fine in this routine" 
    endif
    
    do i = 1,N_node
       do j = 1,N_dmnsn
          if (modelID==1) node_xy(i,j) = nodeXY_strctTN(struct_No,i,j)
          if (modelID==2) node_xy(i,j) = nodeXY_strctNI(struct_No,i,j)
       enddo
    enddo
    
  end subroutine Strct_to_actual
  
  
  subroutine adjustPress_VaryingConsec_A0orKa_of_NI_andSaveNI(seqNo,struct_No,A0V,kaV,l0V,ksV)
    implicit none
    integer, intent(in) :: seqNo,struct_No
    real*8 , intent(in) :: A0V(1:N_cell),kaV(1:N_cell)
    real*8 , intent(in) :: l0V(1:N_spr) ,ksV(1:N_spr)
    
    integer :: saveWhom
    real*8  :: E
    
    A0 = A0V ; k_area = kaV
    l0 = l0V ; k_spr  = ksV
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    
    call adjustPress_VaryingConsec_A0orKa_of_NI(seqNo)
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    
    call Find_Analytical_and_Numerical_Mismatch
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine adjustPress_VaryingConsec_A0orKa_of_NI_andSaveNI
  
  
  subroutine NI_switch_and_read_strctProps_withAddedCellNI(struct_no,A0V,kaV,l0V,ksV)
    implicit none
    integer, intent(in) :: struct_no
    real*8 , intent(in) :: A0V(1:N_cell),kaV(1:N_cell)
    real*8 , intent(in) :: l0V(1:N_spr) ,ksV(1:N_spr)
    
    real*8  :: E
    integer :: cntrlSt
    real*8, allocatable :: TnsnVal(:),PresVal(:)
    
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
    
    open(unit=78,file='NI_switch_and_readStrctProps.dat',position='append')
    
    A0 = A0V ; k_area = kaV
    l0 = l0V ; k_spr  = ksV
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    call read_strctProps_withAddedCell(struct_no)
    !call making_sameProp_for_segmntd_Membrne(struct_No)
    
    A0(1:N_cell)      = A0_strctNI(struct_no,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(struct_no,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(struct_no,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(struct_no,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(struct_no,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    
    write(78,*) (Frame_NI-1),"Frame_NI"
    
    allocate(TnsnVal(1:N_spr),PresVal(1:N_cell))
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    if (struct_no==13) cntrlSt = 1
    if (struct_no==14) cntrlSt = 2
    if (struct_no==15) cntrlSt = 3
    if (struct_no==16) cntrlSt = 4
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    call Find_Analytical_and_Numerical_Mismatch
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(78)
    
  end subroutine NI_switch_and_read_strctProps_withAddedCellNI
  
  
  subroutine NI_switch_and_borrow_strctPropsfrom_HlfCycl1(struct_no,A0V,kaV,l0V,ksV)
    implicit none
    integer, intent(in) :: struct_no
    real*8 , intent(in) :: A0V(1:N_cell),kaV(1:N_cell)
    real*8 , intent(in) :: l0V(1:N_spr) ,ksV(1:N_spr)
    
    real*8  :: E
    integer :: structToRead
    integer :: i,j,jmax
    integer :: cntrlSt,saveWhom
    real*8, allocatable :: TnsnVal(:),PresVal(:)
    
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
    
    open(unit=79,file='NI_switch_and_borrowStrctProp.dat',position='append')
    
    A0 = A0V ; k_area = kaV
    l0 = l0V ; k_spr  = ksV
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    structToRead = struct_no-1
    call read_strctProps_withAddedCell(structToRead)
    
    A0(1:N_cell)      = A0_strctNI(structToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(structToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(structToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(structToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(structToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    do i = 1,3
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if(i==1) write(*,*) k_spr(j),l0(j),j,"ks,l0,j"
          if(i==2) write(*,*) k_area(j),A0(j),j,"ka,A0,j"
          if(i==3) write(*,*) CgXNode(j),CgYNode(j),j,"CgX-YNode,j"
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    
    !after equilibrating the structure should match to go over pulley, now save this for strct3 properties
    
    A0_strctNI(struct_no,1:N_cell)  = A0_strctNI(structToRead,1:N_cell)
    ka_strctNI(struct_no,1:N_cell)  = ka_strctNI(structToRead,1:N_cell)
    l0_strctNI(struct_no,1:N_spr)   = l0_strctNI(structToRead,1:N_spr)
    ks_strctNI(struct_no,1:N_spr)   = ks_strctNI(structToRead,1:N_spr)
    CgX_strctNI(struct_no,1:N_node) = CgX_strctNI(structToRead,1:N_node)
    CgY_strctNI(struct_no,1:N_node) = 0.0d0
    
    saveWhom = 3 !for TN or NI or Both selection
    call saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    
    allocate(TnsnVal(1:N_spr),PresVal(1:N_cell))
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    if (struct_no==15)   cntrlSt = 3
    if (struct_no.ne.15) write(*,*) "cntrlst mismatch"
    if (struct_no.ne.15) stop
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    write(79,*) (Frame_NI-1),"Frame_NI"
    
    call Find_Analytical_and_Numerical_Mismatch
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(79)
    
  end subroutine NI_switch_and_borrow_strctPropsfrom_HlfCycl1
  
  
  subroutine print_trnsfrm_at_diffrntCellsMeet
    implicit none
    integer            :: i,j,jmax,N_itm
    character(len=100) :: flnm,flnmbr,full_flnm
    
    if (modelID==1) flnm   ='TrnsfrmsTNAtCellsMeetEquals'
    if (modelID==2) flnm   ='TrnsfrmsNIAtCellsMeetEquals'
    
    write(flnmbr,'(i2.2,a)') CellsMeet,'.dat'
    
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=87,file=trim(adjustl(full_flnm)))

    N_itm=2
    
    do i = 1,N_itm
       if (i==1) write(87,fmt=*) "node_spr:"
       if (i==2) write(87,fmt=*) "node_area:"
       
       do j = 1,N_node
          !write(87,fmt=*) "node_nm =",j
          
          if (i==1) then
             write(87,fmt=*) j,node_spr(j,0:max_spr_node)
          elseif (i==2) then
             write(87,fmt=*) j,node_area(j,0:max_area_node)
          endif
          
       enddo
       write(87,fmt=*) " "
       
    enddo
    
    
    do i = 1,N_itm
       if (i==1) write(87,fmt=*) "spr_node:"
       if (i==2) write(87,fmt=*) "spr_area:"
       
       do j = 1,N_spr
          !write(87,fmt=*) "spr_nm =",j
          
          if (i==1) then
             write(87,fmt=*) j,spr_node(j,0:max_node_spr) 
          elseif (i==2) then
             write(87,fmt=*) j,spr_area(j,0:max_area_spr)
          endif
          
       enddo
       write(87,fmt=*) " "
       
    enddo
    
    do i = 1,N_itm
       if (i==1) write(87,fmt=*) "area_node:"
       if (i==2) write(87,fmt=*) "area_spr:"
       
       do j = 1,N_cell
          !write(87,fmt=*) "area_nm =",j
          
          if (i==1) then
             write(87,fmt=*) j,area_node(j,0:max_node_area) 
          elseif (i==2) then
             write(87,fmt=*) j,area_spr(j,0:max_spr_area)
          endif
          
       enddo
       write(87,fmt=*) " "
       
    enddo
    
       
    do i = 1,N_node
       write(87,fmt=*) "Nlist for  node_nm =",i
       
       jmax = node_area(i,0)
       
       do j = 1,jmax
          write(87,fmt=*) Nlist(i,j,1:2)
       enddo
       
       write(87,fmt=*) " "
    enddo
    

    close(87)
    
  end subroutine print_trnsfrm_at_diffrntCellsMeet
  
  
  subroutine run_simltn_to_show_strctr_formt(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8, allocatable :: A0_fs(:),l0_fs(:),ka_fs(:),ks_fs(:) !fs=formed structure
    
    integer :: i,j,jmax
    integer :: elmntTyp,N_itrtn
    real*8  :: stpSize
    
    allocate(A0_fs(1:N_cell),l0_fs(1:N_spr))
    allocate(ka_fs(1:N_cell),ks_fs(1:N_spr))
    
    A0_fs(1:N_cell) = A0(1:N_cell) ; ka_fs(1:N_cell) = k_area(1:N_cell)  
    l0_fs(1:N_spr)  = l0(1:N_spr)  ; ks_fs(1:N_spr)  = k_spr(1:N_spr)
    
    A0(1:N_cell) = A(1:N_cell) ; l0(1:N_spr) = l(1:N_spr)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    write(*,*)
    N_itrtn  = 5 
    stpSize  = (1.0d0)/(N_itrtn-1)
    
    do elmntTyp = 1,2
       
       if (elmntTyp == 1) jmax = N_cell
       if (elmntTyp == 2) jmax = N_spr
       
       !do i = 1,N_itrtn
          
          !do j = 1,jmax
           !  if (elmntTyp == 1) A0(j) = + (i-1)*
             
          !enddo
          
       !enddo
       
    enddo
    
    
  end subroutine run_simltn_to_show_strctr_formt
  
  
  
  subroutine force_prtrbtn_applied_at_boundry_of_CS2(callLoctn,strctToRead,ExpNo,FrmNo)
    !improved routine of ForcePrtrbtnAnalysis_atCS2, which was in nodeInsrtd_cntrlStates.f08
    
    implicit none
    integer, intent(in)    :: callLoctn
    integer, intent(in)    :: strctToRead
    integer, intent(inout) :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8 :: E
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    real*8, allocatable :: kaV(:),A0V(:),Aval(:)
    real*8, allocatable :: ksV(:),l0V(:),lval(:)
    real*8, allocatable :: ks_stre(:)
    
    integer :: pullingORpushing,caseNo,ExpNoCase,FrmNoCase
    real*8  :: hmuchCgX
    integer :: NonUniorUni,TnsnOrExtrpltn
    integer :: flnmInit,flnmFinl
    integer :: incrORdcrs
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    allocate(kaV(1:N_cell),A0V(1:N_cell),Aval(1:N_cell))
    allocate(ksV(1:N_spr) ,l0V(1:N_spr) ,lval(1:N_spr))
    allocate(ks_stre(1:N_spr))
    
    if (callLoctn == 1) then
       
       call read_strctProps_withAddedCell(strctToRead)
       
       A0(1:N_cell)     = A0_strctNI(strctToRead,1:N_cell)
       k_area(1:N_cell) = ka_strctNI(strctToRead,1:N_cell)
       l0(1:N_spr)      = l0_strctNI(strctToRead,1:N_spr)
       k_spr(1:N_spr)   = ks_strctNI(strctToRead,1:N_spr)
       CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
       
    elseif (callLoctn == 2) then
       
       write(*,*) CyclVal,strctToRead,"CyclVal,strctToRead, IMPORTANT VALUE CHECK IT"
       
       A0(1:N_cell)     = A0_cycl(CyclVal,strctToRead,1:N_cell) !strctToRead must be 2
       k_area(1:N_cell) = ka_cycl(CyclVal,strctToRead,1:N_cell)
       l0(1:N_spr)      = l0_cycl(CyclVal,strctToRead,1:N_spr)
       k_spr(1:N_spr)   = ks_cycl(CyclVal,strctToRead,1:N_spr)
       write(*,*) CgXNode(1),"CgXNode(1) SHOULD BE ZERO HERE"
       !CgXNode does not need to be read, they are ZERO anyways
       
    endif
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    
    !**** go to subroutine ForcePrtrbtnAnalysis_atCS2 *****
    ! Adjustment Needed For Reading structural_info's
    
    call find_the_invagniating_and_invaginated_cell
    call find_the_frst_four_cellPair
    call find_the_last_three_cellPair
    
    !***********
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    call store_props_bfr_chng(kaV,A0V,Aval,ksV,l0V,lval)
    
    pullingOrpushing = 1
    caseNo           = 3
    call get_CgVal_ExpNo_and_FrmNo(pullingORpushing,caseNo,hmuchCgX,ExpNoCase,FrmNoCase)
    
    NonUniorUni    = 2 ! For NonUni=1 and Uni=2
    TnsnOrExtrpltn = 1 ! For Tnsn=1 and Extrplt=2
    
    if (NonUniOrUni==1) then
       
       call apply_force_and_adjust_anomalies_bfr_implmnting_meetingAlgs(pullingOrpushing,hmuchCgX,caseNo,ExpNoCase,FrmNoCase)
       
       if (TnsnOrExtrpltn==1) then
          write(*,*) (FrmNo-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=1 "
          incrORdcrs = 1
          call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          
       elseif (TnsnOrExtrpltn==2) then
          flnmInit = strctToRead-1 ; flnmFinl = strctToRead
          ExpNo = ExpNoCase ; FrmNo = FrmNoCase ! imp for Extrpltn
          
          call prep_vars_for_Extrplt(callLoctn,NonUniorUni,flnmInit,flnmFinl)
       endif
       
    elseif (NonUniOrUni==2) then
       continue
    endif
    
    if (TnsnOrExtrpltn==1) stop !stopping as I dont wanna go it to Extrpltn
    if (TnsnOrExtrpltn==2) continue
    
  end subroutine force_prtrbtn_applied_at_boundry_of_CS2
  
  subroutine prep_vars_for_Extrplt(callLoctn,NonUniorUni,flnmInit,flnmFinl)
    implicit none
    integer, intent(in) :: callLoctn
    integer, intent(in) :: NonUniorUni
    integer, intent(in) :: flnmInit,flnmFinl
    
    integer :: whatPropsToReadforExtrplt
    integer :: CSvalI,CSvalF
    
    
    allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
    allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
    
    A0_I(1:N_cell)=1.0d30 ; ka_I(1:N_cell)=1.0d30 ; l0_I(1:N_spr)=1.0d30 ; ks_I(1:N_spr)=1.0d30
    A0_F(1:N_cell)=1.0d30 ; ka_F(1:N_cell)=1.0d30 ; l0_F(1:N_spr)=1.0d30 ; ks_F(1:N_spr)=1.0d30
    
    if (NonUniorUni==1) whatPropsToReadforExtrplt = 1
    if (NonUniorUni==2) whatPropsToReadforExtrplt = 2   
    
    if (callLoctn==1) then
       
       if (whatPropsToReadforExtrplt == 1) then
          
          A0_F(1:N_cell) = A0_strctNI(flnmFinl,1:N_cell)
          A0_I(1:N_cell) = A0_strctNI(flnmInit,1:N_cell)
          ka_F(1:N_cell) = ka_strctNI(flnmFinl,1:N_cell)
          ka_I(1:N_cell) = ka_strctNI(flnmInit,1:N_cell)
          l0_F(1:N_spr)  = l0_strctNI(flnmFinl,1:N_spr)
          l0_I(1:N_spr)  = l0_strctNI(flnmInit,1:N_spr)
          ks_F(1:N_spr)  = ks_strctNI(flnmFinl,1:N_spr)
          ks_I(1:N_spr)  = ks_strctNI(flnmInit,1:N_spr)
          
       elseif (whatPropsToReadforExtrplt == 2) then
          
          call read_Init_Finl_props_from_propsFlnm(flnmInit,flnmFinl,A0_F,A0_I,ka_F,ka_I,l0_F,&
               l0_I,ks_F,ks_I)
       endif
       
    elseif (callLoctn==2) then ! callLoctn=2 is frm calling insd RunCyclicProcessForAddedCellNI
       
       if (whatPropsToReadforExtrplt == 1) then
          
          A0_F(1:N_cell) = A0_cycl(cyclVal,flnmFinl,1:N_cell)
          A0_I(1:N_cell) = A0_cycl(cyclVal,flnmInit,1:N_cell)
          ka_F(1:N_cell) = ka_cycl(cyclVal,flnmFinl,1:N_cell)
          ka_I(1:N_cell) = ka_cycl(cyclVal,flnmInit,1:N_cell)
          
          l0_F(1:N_spr)  = l0_cycl(cyclVal,flnmFinl,1:N_spr)
          l0_I(1:N_spr)  = l0_cycl(cyclVal,flnmInit,1:N_spr)
          ks_F(1:N_spr)  = ks_cycl(cyclVal,flnmFinl,1:N_spr)
          ks_I(1:N_spr)  = ks_cycl(cyclVal,flnmInit,1:N_spr)
          
       elseif (whatPropsToReadforExtrplt == 2) then
          write(*,*) "HAVE TO ADJUST FOR UNIFORM ONES"
          write(*,*) "should not enter here as of now, sb: prep_vars_for_Extrplt"
       endif
       
    endif
    
  end subroutine prep_vars_for_Extrplt
  
  
  subroutine readORwritSysBfrEditingInitiationPrjct(readORwrite)
    implicit none
    integer, intent(in) ::  readORwrite
    
    integer :: i,j,jmax
    integer :: N_itm,entityNum
    
    
    open(unit=534,file='bfr_editing_for_initiationProjecrt.dat')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (readORwrite == 1) then
             
             if (i==1) read(534,*) k_area(j),A0(j),entityNum
             if (i==2) read(534,*) k_spr(j),l0(j),entityNum
             if (i==3) read(534,*) CgXNode(j),entityNum
             
          elseif (readORwrite == 2) then
             
             if (i==1) write(534,*) k_area(j),A0(j),entityNum
             if (i==2) write(534,*) k_spr(j),l0(j),entityNum
             if (i==3) write(534,*) CgXNode(j),entityNum
             
          endif
          
       enddo
       
       if (readORwrite==2) write(534,*) " "
       
    enddo
    
    close(534)
    
  end subroutine readORwritSysBfrEditingInitiationPrjct
  
  
  subroutine save_ka_andCg_and_makethoseZero(kaSaveV,CgXSaveV,CgYSaveV)
    implicit none
    real*8, intent(out) :: kaSaveV(1:N_cell)
    real*8, intent(out) :: CgXSaveV(1:N_node), CgYSaveV(1:N_node)
    
    integer :: i,j,imax,jmax
    
    kaSaveV  = k_area
    CgXSaveV = CgXNode
    CgYSaveV = CgYNode
    
    k_area  = 0.0d0
    CgXNode = 0.0d0
    CgYNode = 0.0d0
    
    open(unit=911,file='kaCgXCgYSave.dat')
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_node
       
       do j = 1,jmax
          if (i==1) write(911,*) kaSaveV(j),j,"ka"
          if (i==2) write(911,*) CgXSaveV(j),CgYSaveV(j),j,"Cg"
       enddo
       
       write(911,*) " "
    enddo
    
    close(911)
    
  end subroutine save_ka_andCg_and_makethoseZero
  
  subroutine updt_nodeValsFrmPrvNIfiles
    implicit none
    character(len=100) :: flnm
    integer :: i,primNode
    integer :: nodeL,nodeR
    
    write(flnm,'(a,i3.3,a)')"NI_modelInitiatnNW",(Frame_NI-1),'.dat'
    write(*,*)trim(adjustl(flnm))
    
    open(unit=765,file=trim(adjustl(flnm)))
    
    primNode = (N_cell+1)*2
    
    do i = 1,N_node
       
       if (i.le.primNode) then
          read(765,*) node_xy(i,1:2)
          
       elseif (i.gt.primNode .and. i.le.(primNode+2)) then
          
          if ((i-primNode)==1) then
             nodeL=(Hlf_Ncell)*2+1 ; nodeR = (Hlf_Ncell+1)*2+nodeL
          endif
          
          if ((i-primNode)==1) then
             
             write(*,*) node_xy(i,1:2),i,"node chk11"
             node_xy(i,1:2)   = node_xy(nodeL,1:2) 
             node_xy(nodeL,2) = node_xy(i,2) - ydis_frmVitlnMem
             write(*,*) node_xy(i,1:2),i,"node chk12"
             write(*,*) node_xy(nodeL,1:2),"nodeL"
             
          elseif ((i-primNode)==2) then
             
             write(*,*) node_xy(i,1:2),i,"node chk21"
             node_xy(i,1:2)   = node_xy(nodeR,1:2)
             node_xy(nodeR,2) = node_xy(i,2) - ydis_frmVitlnMem
             write(*,*) node_xy(i,1:2),i,"node chk22"
             write(*,*) node_xy(nodeR,1:2),"nodeR"
             
          endif
          
       else
          read(765,*) node_xy(i,1:2)
       endif

       write(*,*) node_xy(i,1:2),i,"in updt"
       
    enddo
    
    close(765)
    
  end subroutine updt_nodeValsFrmPrvNIfiles
  
  
  
  subroutine taken_from_inside_of_tbsitsS1wc ! == trnsfrm_btn_strcts_inTwoStps_S1_with_cndn
    implicit none
    real*8  :: nodeXYpnt1(1:150,1:2)=-100.0d0,nodeXYpnt2(1:150,1:2)=-100.0d0
    integer :: forCellsMeet
    integer :: redcSimlnNum,cnt
    integer :: ApclSpr1(1:(NAEC_Apcl+1)) ,BsalSpr1(1:(NAEC_Bsal+1))
    integer :: LtrlSprB1(1:(NAEC_Ltrl+1)),LtrlSprB2(1:(NAEC_Ltrl+1))
    integer :: EndCellNum,cntMax
    integer :: cellNm,ApBsLt
    real*8  :: hmL0
    integer :: PCSorICS
    
    SystemTyp = 1
    call adding_node_inIC_ApclMem_and_redefine_system
    !call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt,Frame)
    Frame = Frame+1
    call switch_to_NI_model
    call updt_nodeValsFrmPrvNIfiles
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
    nodeXYpnt2(1:N_node,1:2) = node_xy(1:N_node,1:2)
    call save_only_NI_config
    
    forCellsMeet = 1           ! IMPORTANT LINE 
    redcSimlnNum = 1
    write(*,*) CellsMeet,"CMv"
    
    
    if (forCellsMeet==1) then
       
       write(*,*) Frame_NI-1,"Frm_NI bfr"
       call apply_Force_aft_ManuallyShifting_Nodes(Exprmnt,Frame)
       write(*,*) Frame_NI-1,"Frm_NI aft"
       
       CellsMeet = 1
       call combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM1
       call save_only_NI_config
       call Equilibrate_only_NI_model
       
    elseif (forCellsMeet==2) then
       
       call apply_Force_aft_ManuallyShifting_Nodes(Exprmnt,Frame)
       call Apcl_shorten_Adjacent_of_IC_and_incr_vrtcl_force_to_hold(Exprmnt,Frame)
       PCSorICs=2 ; call tilt_reduction_of_the_boundry(forCellsMeet,PCSorICS)
       
       CellsMeet = 2
       call combining_cntrl_and_neighApcl_sprs_and_redefine_system_aft_CM2
       call save_only_NI_config
       call Equilibrate_only_NI_model
       
    endif
    
    !Frame_NI=59
    !call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    !call Equilibrate_only_NI_model
    
  end subroutine taken_from_inside_of_tbsitsS1wc
  
  
  subroutine taken_from_inside_of_tbsitsS1wc_to_modelVF
    implicit none
    real*8  :: nodeXYpnt1(1:150,1:2)=-100.0d0,nodeXYpnt2(1:150,1:2)=-100.0d0
    integer :: forCellsMeet
    integer :: redcSimlnNum,cnt
    integer :: ApclSpr1(1:(NAEC_Apcl+1)) ,BsalSpr1(1:(NAEC_Bsal+1))
    integer :: LtrlSprB1(1:(NAEC_Ltrl+1)),LtrlSprB2(1:(NAEC_Ltrl+1))
    integer :: EndCellNum,cntMax
    integer :: cellNm,ApBsLt
    real*8  :: hmL0,Fctr
    integer :: PCSorICS,node_typ_chng_case
    integer :: N_diffForces
    integer :: ApclSide_PrtialSgmntn
    real*8  :: VF_area1,VF_area2,VF_area3
    
    !call Tst_for_pulley_point_force
    
    addedNCPair = 3
    
    node_typ_chng_case=4 ; call redefine_system_for_changing_only_node_typ(node_typ_chng_case)
    call Equilibrate_only_TN_model(Exprmnt,Frame)
    write(*,*) N_node,N_spr,N_cell,"VF0" ; write(*,*) (Frame-1),"Frame"
    
    NI_incldAS_woPP = 1 ; IN_apclSideCaseNo = 1 ; call switch_to_NI_model_withNodesInAllmembranes
    Frame_NI=1 ; call Equilibrate_only_NI_model   !'after switching all sides'
    
    call High_ksTst_onExistingCondition() ;
    !Frame_NI=8 ; call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI) ;call Equilibrate_only_NI_model
    cnt=1 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt)
    
    VF_regionModelled=1 ; call add_cellFor_VitellineFluid_region_and_redefine_system
    call WO_Equilibrate_onlysaveConfig_withVF_region ; call Equilibrate_only_NI_model_withVF_region
    cnt=6 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt)
    
    call SRyp_Test ; write(*,*) " " ; call print_prop_of_system_wtSRyp !; stop 'VF'
    !call Tst_LtrlSpr_IC_Shortn_wtVF_region ; stop 'Aft Ltrl Short modelVF'
    call print_all_trnsfrms_WOtxt
    
    NI_AS_wtCrtclSurfc=1 ; NCP_CrtclApSrfc=1 ; NAEC_ApclCrtcl = 5
    addApNodesPerApSpr = NAEC_ApclCrtcl-NAEC_Apcl ; addApNodesPerSgmntdApSpr = 1
    
    call change_num_of_InsrtdNodes_inApclSide_and_redefine_system()
    cnt=7 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt) !; stop 'CRTCL DESGND'
    
    call make_the_lateralTnsn_IC_equaltoDesirdval() ; stop 'Aft make_the_lateralTnsn_IC_equaltoDesirdVal' 
    
    !Frame_NIVF=1
    !N_diffForces=21; call apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion(N_diffForces)
    !Frame_NIVF=21 ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    !call Equilibrate_only_NI_model_withVF_region
    
    !Fctr=3.00d0 ; call Effct_of_VF_manpltng_ka(Fctr)
    
    !call manplts_for_JHT_articl_Img_1M ; stop 'AFT JHT'
    !call Diffrnrt_Level_Forces_AtFrstJointNode ; stop 'AFT JHT'
    call readFrm_FrceDownLtrlNC1_manpltn_withDiffrntPrtrbtn() ; stop 'AFT JHT1'
    
    !call manpltns_for_CF_WT_Nrt_P_2G5C_m()
    !call manpltns_for_CF_WT_Nrt_P_1G1B1_m()
    !call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa ; stop 'aft P_1G1B1_m'
    
    
  end subroutine taken_from_inside_of_tbsitsS1wc_to_modelVF
  
  !N_diffForces=4 ; call apply_force_at_Apcl_surfc_ofIC_bfr_adding_VFregion(N_diffForces)
  !VF_area1=-1.0d30 ; VF_area2=-1.0d30 ; VF_area3=-1.0d30    
  !VF_area1 = vitelline_fluid_region(N_cell,node_xy)
  !call vitelline_fluid_region_methd3(VF_area2,node_xy)
  !call vitelline_fluid_region_methd3_write(VF_area3,node_xy) 
  !write(*,*) VF_area1,VF_area2,VF_area3,"VF_area S"
  
  !cnt=3 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt) 
  
  subroutine Tst_LtrlSpr_IC_Shortn_wtVF_region
    implicit none
    integer :: cellNm
    
    CgXNode=0.00d0  ; CgYNode=0.00d0 ; call Equilibrate_only_NI_model_withVF_region
    cellNm=Hlf_Ncell; call Ltrl_Shorten_of_Cell_NIsys(cellNm)
    
  end subroutine Tst_LtrlSpr_IC_Shortn_wtVF_region
  
  subroutine Tst_for_pulley_point_force()
    implicit none
    character(len=200) :: flplcV
    integer            :: ExpNo,FrmNoToBeRead,sprNo,cnt
    logical            :: lgcl_rdfn
    real*8             :: tol_Rdfn,E
    integer            :: LftNeigh=0,RghtNeigh=0
    
    write(flplcV,*)&
         "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/NI_DataInitPhase/"
    ExpNo = Exprmnt ; FrmNoToBeRead = 7
    call read_config_and_start_simlnFrm_there_flpceLoc(flplcV,ExpNo,FrmNoToBeRead)
    call Equilibrate_only_TN_model(Exprmnt,Frame)
    
    lgcl_rdfn=.False. ; tol_Rdfn=0.05d0 ; call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
       
    if (lgcl_rdfn.eqv..True.) then
       write(*,*) " ";write(*,*) "THE CELLS MEET 1st Time !!!!!!!";write(*,*) " "
       CellsMeet = 1
    elseif (lgcl_rdfn.eqv..False.) then
       stop ':( Cells did not meet yet, stop here to adjust (1st)'
    endif
    
    call taking_lftrghtOriginNeigh_downwards(LftNeigh,RghtNeigh)
    sprNo = cntrlCell_Spring(1) ; write(*,*) sprNo,cntrlCell_Spring(1),"cntrlCell_Spring 1"
    dmlshDecsn = 1 ; call demolish_spring_and_redefine_system(sprNo)
    
    call Equilibrate_only_TN_model(Exprmnt,Frame)
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0) ; write(*,*) E,"E" ; call get_gradient(node_xy,l0,A0,gradient)
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1 ; cnt=5 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt)
    
    ExpNo = Exprmnt_NI ; FrmNoToBeRead = 11
    call read_config_and_start_simlnFrm_there_flpceLoc(flplcV,ExpNo,FrmNoToBeRead)
    call Equilibrate_only_NI_model
    call displace_pulley_virtually_to_check_force
    
    stop 'aft pulley_point_tst'
    
  end subroutine Tst_for_pulley_point_force
  
  subroutine SRyp_Test
    implicit none
    real*8  :: Energy_Val,Es,Ea,Eg,Eb,Er
    integer :: cnt
    
    call deallocate_all_gradient_variables
    
    SRyp_lgcl=.True.        ! Activating SRyp
    
    call get_positionalNodes()      ; call allocate_and_initialize_SRypVars
    call get_SRYP_props             ; call get_actvtn_fctrs() 
    call get_all_gradient_variables
    
    Energy_Val = Energy(node_xy,l0,A0) ; write(*,*) Energy_Val,"Ev"
    Es         = spr_E(node_xy,l0)     ; write(*,*) Es,"Es"
    Ea         = area_E(node_xy,A0)    ; write(*,*) Ea,"Ea"
    Eg         = grvtnl_E(node_xy)     ; write(*,*) Eg,"Eg"
    Eb         = bend_E(node_xy)       ; write(*,*) Eb,"Eb"
    Er         = SRyp_E(node_xy)       ; write(*,*) Er,"Er"
    
    write(*,*) Energy_Val,(Energy_Val-(Es+Ea+Eg+Eb)),"chk bfr stop"
    cnt=3 ; call Find_Analytical_and_Numerical_Mismatch_allInfo(cnt) 
    
  end subroutine SRyp_Test
  
  
  subroutine deformation_tests_using_one_or_multiple_cells
    implicit none
    integer :: node_typ_chng_case,i
    
    node_typ_chng_case=1 ; call redefine_system_for_changing_only_node_typ(node_typ_chng_case)
    call Equilibrate_bothTN_NImodel(Exprmnt,Frame)
    
  end subroutine deformation_tests_using_one_or_multiple_cells
  
  
  !write(*,*) "Frame NI checkpoint 011:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 012:",(Frame_NI-1) 
  !write(*,*) "Frame NI checkpoint 021:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 022:",(Frame_NI-1)  
  !write(*,*) "Frame NI checkpoint 031:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 032:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 041:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 042:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 051:",(Frame_NI-1)
  !write(*,*) "Frame NI checkpoint 052:",(Frame_NI-1)
  
end module build_struct
