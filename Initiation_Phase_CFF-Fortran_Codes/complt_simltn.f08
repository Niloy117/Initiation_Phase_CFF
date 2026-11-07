module Complt_Simltn
  
  use sys_building_info !input_frm_bash_what_to_build.f08
  use cell_info
  use end_info
  use system_parameters 
  
  use generating_the_shape !generating_system.f08
  use perturbing_the_shape !generating_system.f08
  
  use changing_parameters !manipulating_properties.f08
  
  use node_variables
  use spring_variables
  use area_variables
  
  use spr_to_other_transfrm
  use node_to_other_transfrm
  use area_to_other_transfrm
  use neighbour_info_and_trnsfrmInfo_dependent_info
  use transfrm_info
  
  use conversion_routines
  use Energy_module
  use gradient_module
  use PenF_module
  use PenF_minimizer_mod
  use Wfunc_and_its_derivative
  
  use mltiple_calls_together
  use moving_coordnte_info
  use moving_coordnte_variables
  
  use force_calc
  use storing_changing_restoring_routines
  use Manipulate_and_equilibrate
  
  use calls_for_tests
  use nodeInsrtd_cntrlStates
  use Adding_cells
  use build_struct
  
  use readjust_struct
  use mltpl_cycle
  use redefining_system_module
  
contains
  
  subroutine stage1_type1
    implicit none
    integer :: iMain,whichEnd
    integer :: strct_INPT,strct_OTPT
    real*8  :: hm
    integer :: i1,j,jmax,cnt
    integer :: readBoth=-1
    real*8  :: Es=0.0d0,Ea=0.0d0
    integer :: saveWhom
    integer :: readORwrite
    integer :: ExpNo,FrmNo
    integer :: strctBfrMeet,strctAftMeet
    integer :: printORread
    
    integer, allocatable :: dataUpdtd(:)
    
    N_struct = 19 !16
    N_seq    = 10
    SaveOrReadData = 3
    !nodeInsrtnStrts=1
    
    allocate(dataUpdtd(1:N_seq))
    dataUpdtd(1:N_seq) = 0
    dataUpdtd(1:8)     = 1
    !dataUpdtd(2)       = 0
    !dataUpdtd(5)       = 0
    
    if (SaveOrReadData==1) then
       Exprmnt = 21 ; Frame=1
       
       do iMain=1,N_struct
          !if (iMain==14) nodeInsrtnStrts=1
          
          strctNo=iMain
          call build_struct_And_saveProp_S1
       enddo
       
    elseif (SaveOrReadData==2 .or. SaveOrReadData==3) then
       
       !call checkAreaRoutineChk
       
       if (SaveOrReadData==2) Exprmnt = 23 ; Frame=1
       if (SaveOrReadData==3) Exprmnt = 27 ; Frame=1
       
       do iMain = 1,N_seq
          
          write(*,*) iMain,"h1"
          
          strct_INPT = (2*iMain-1)-((iMain-1)/2)
          strct_OTPT = (2*iMain-0)-((iMain-1)/2)
          
          if ((iMain==1).or.(iMain==N_seq)) then
             readBoth = 1
          else
             if (mod(iMain,2)==1) readBoth = 0
             if (mod(iMain,2)==0) readBoth = 1
          endif
          
          write(*,*) strct_INPT,strct_OTPT,"h2"
          
          if (iMain==N_seq-1)  hlfCycl = 1
          if (iMain==N_seq)    hlfCycl = 2
          
          call adjusting_routines_in_trnsfrm_strcts_S1(strct_INPT,strct_OTPT)
          
          write(*,*) SaveOrReadData,"h3"
          
          if (SaveOrReadData==2) then
             
             if (readBoth==1) then
                call read_strctProps(strct_INPT)
                call read_strctProps(strct_OTPT)
             elseif (readBoth==0) then
                call read_strctProps(strct_OTPT)
             endif
             
          elseif (SaveOrReadData==3) then
             
             call get_active_region_arrays(iMain)
             call get_active_region_S1
             
             if (readBoth==1) then
                
                if (dataUpdtd(iMain)==0) then
                   
                   if (iMain==1) then
                      call read_strctPropsA0l0Cg_withActvRgn(strct_INPT)
                   elseif (iMain==2 .and. dataUpdtd(iMain-1)==1) then
                      
                      ! SPECIAL CASE, [creates the strctPropsAddedTN3S1t1 and
                      ! strctPropsAddedCellNI from strctPropsAddedCell(TN/NI)2S1T1
                      ! that is from before the meet structures]
                      
                      strctBfrMeet = strct_INPT-1
                      strctAftMeet = strct_INPT
                      call read_strctPropsA0l0Cg_FrmStrcBfrMeet(strctBfrMeet,strctAftMeet)
                      !write(*,*) "STOPPED AFT CALLING THE ROUTINE SPEC"
                      !stop
                      
                   else
                      call read_strctPropsA0l0Cg_withActvRgn(strct_INPT)
                   endif
                   
                   call read_strctPropsA0l0Cg_withActvRgn(strct_OTPT)
                   
                elseif (dataUpdtd(iMain)==1) then
                   call read_strctPropsA0l0Cg_Updtd(strct_INPT)
                   call read_strctPropsA0l0Cg_Updtd(strct_OTPT)
                endif
                
             elseif (readBoth==0) then
                if (dataUpdtd(iMain)==0) call read_strctPropsA0l0Cg_withActvRgn(strct_OTPT)
                if (dataUpdtd(iMain)==1) call read_strctPropsA0l0Cg_Updtd(strct_OTPT)
             endif
             
             call deallocate_active_region_arrays
             
          endif
          
          saveWhom = 1
          !call saveA0l0_ofStrct_AddedCell_TNorNI(strct_INPT,saveWhom)
          !call incrAllA0val_andReduceAll_l0Val(Exprmnt,Frame)
          
          if (iMain.ge.9) then
             write(*,*) "ENTERING and STOPPING greater than iMain = 9"
             !stop
             
             nodeInsrtnStrts = 1 !make it 1 if it needed to
             call trnsfrm_btwn_strcts_S1_withPressAdjstmnt(strct_INPT,strct_OTPT,iMain)
             
             
          elseif (iMain.lt.9) then
             
             write(*,*) "Entering for iMain =",iMain
             
             nodeInsrtnStrts  = 1
             RUN_WTorWO_Force = 1 ! IMPORTANT 
             
             if (RUN_WTorWO_Force == 0) then
                
                if (mod(iMain,2)==1) then
                   
                   printORread=1 ; call printORread_NIorTN_trnsfrms_inFiles(printORread)
                   call switch_to_NI_model
                   printORread=1 ; call printORread_NIorTN_trnsfrms_inFiles(printORread)
                   call deallocate_repetitive_arrays
                   call switchback_to_TN_model
                   
                endif
                
                call trnsfrm_btwn_strcts_inTwoSteps_S1(strct_INPT,strct_OTPT) 
                
                if (iMain==8) then
                   write(*,*) "stopping after initiation phase without force case"
                   stop
                endif
                
             elseif (RUN_WTorWO_Force == 1) then
                call runningFor_WT_Force_InitPhase(iMain,whichEnd,dataUpdtd,strct_INPT,strct_OTPT,saveWhom)
             endif
             
             !if (iMain==8) then
                !readORwrite=2 ; call readORwritSysBfrEditingInitiationPrjct(readORwrite)
                !will uncomment the top line when we continue for the progression
                !stage with previously used system
                !stop
             !  continue
             !endif
             
          endif
          
          !saveWhom = 1
          !call saveA0l0_ofStrct_AddedCell_TNorNI(strct_OTPT,saveWhom)
          !if (dataUpdtd(1)==1) call saveA0l0_ofStrctAft_1stMeet_fromStrct_WOmeet()
          
          write(*,*) select_xy,"select_xy in iMain=",iMain
          
       enddo
       
       write(*,*) "Bfr Added"
       
       call switchto_ACmodel_and_mayOrmaynt_switchback
       call RunCyclicProcessForAddedCellNI
       
       stop
       Exprmnt=28 ; Frame=1 ; call Property_Trnsfr_ToNextCell
       !Exprmnt=29 ; Frame=1 ; call PrgrsnStgCyclicMovement
       !Exprmnt=31 ; Frame=1 ; call PrgrsnStgCyclicMovement_Separately
       
    endif
    
  end subroutine stage1_type1
  
  
  subroutine runningFor_WT_Force_InitPhase(iMain,whichEnd,dataUpdtd,strct_INPT,strct_OTPT,saveWhom)
    implicit none
    integer, intent(in)    :: iMain
    integer, intent(inout) :: whichEnd
    integer, intent(in)    :: dataUpdtd(1:N_seq)
    integer, intent(in)    :: strct_INPT,strct_OTPT
    integer, intent(inout) :: saveWhom
    
    write(*,*) "getting inside runningFor_WT"
    write(*,*) iMain,whichEnd,'[',dataUpdtd,']',strct_INPT,strct_OTPT,saveWhom,"chking inside runningFor_WT"
    !stop
    
    if (dataUpdtd(iMain)==0) then
       
       whichEnd = -1
       if ((iMain==1) .or. (iMain.ge.2)) then
          call change_shape_Initiation_wt_ExpFrmNo(iMain,whichEnd)
          saveWhom=1 ; call saveA0l0_ofStrct_AddedCell_TNorNI(strct_INPT,saveWhom)
       endif
       
    endif
    
    write(*,*) strct_INPT,strct_OTPT,iMain,dataUpdtd(iMain),"values bfr trnsfrm"
    
    if (iMain.ne.1) call trnsfrm_btwn_strcts_inTwoSteps_S1(strct_INPT,strct_OTPT) !regular cases
    if (iMain==1)   call trnsfrm_btwn_strcts_inTwoSteps_S1_with_condn(iMain,strct_INPT,strct_OTPT)!pulling bfr dimnsh spr
    
    if (dataUpdtd(iMain)==0) then
       
       whichEnd = +1
       if (iMain==2) then ! strct_INPT=3 (whichEnd=-1) ; strct_OTPT=4 (whichEnd=+1)
          call read_strctPropsA0l0Cg_Updtd(strct_INPT)!readFrm flnm:strctPrpsAddedCellTN [read(43,*) ks_strctTN]
          call get_strctProps_readFrmA0l0Cg_Updtd(strct_INPT)!insrt val,e.g:k_spr(1:N_spr)=ks_strctT(st,1:N_sp)
          call change_shape_Initiation_wt_ExpFrmNo(iMain,whichEnd)
          saveWhom=1 ; call saveA0l0_ofStrct_AddedCell_TNorNI(strct_OTPT,saveWhom)
       elseif (iMain.ge.3) then
          call change_shape_Initiation_wt_ExpFrmNo(iMain,whichEnd)
          saveWhom=1 ; call saveA0l0_ofStrct_AddedCell_TNorNI(strct_OTPT,saveWhom)
       endif
        
    endif
    
  end subroutine runningFor_WT_Force_InitPhase
  
  
  
  subroutine stage1_type2
    implicit none
    
    continue
    
  end subroutine stage1_type2
  
  subroutine stage2_type1
    implicit none
    
    continue
    
  end subroutine stage2_type1
  
  subroutine stage2_type2
    implicit none
    
    continue
    
  end subroutine stage2_type2
  
  subroutine stage3_type1
    implicit none
    
    continue
    
  end subroutine stage3_type1
  
  subroutine stage3_type2
    implicit none
    
    continue
    
  end subroutine stage3_type2
  
  subroutine stage4_type1
    implicit none
    integer :: iMain
    integer :: strct_INPT,strct_TRGT
    real*8  :: hm
    integer :: i1,j,jmax
    
    N_struct = 3
    
    do iMain = (N_struct-1),1,-1
       
       strctNo = iMain
       write(*,*) strctNo
       
       call build_struct_And_saveProp_S4
       !write(*,*) A0(1),iMain,"A01,iMain"
       
       if (iMain.ne.1) then
          call destroy_struct
          call reinitialize_system_params
       endif
       
    enddo
    
    !write(*,*) l0_strct(1,1),l0_strct(2,1),"l0_strct"  
    !write(*,*) l0_strct(1,13),l0_strct(2,10:11),"l0_strct"
    !stop
    
    
    !strct_INPT = 1 !
    !strct_TRGT = 2 !
    
    !call get_thirdStrctProp(strct_INPT,strct_TRGT) !
    !stop
    
    strctNo = 3 !keep this 3 line (this one + below 2) or top 4 line with stop
    call read_strctProps(strctNo)
    call Strct_to_actual(strctNo)
    
    call TwoStep_trnsfrmtn !This line and next if one cycle
    stop
    
    strctNo = 4
    call destroy_and_rebuilt_for_cyclic_trnsfrmtn
    !stop
    strctNo = 1 !!careful at this line
    
    EquilAlgrthm = 1
    
    do CyclNo = 1,3
       call TwoStep_trnsfrmtn
       
       open(unit=124,file='Varying_CgXNode.dat',position='append')
       
       do i1=1,N_node
          write(124,fmt=*) CgXNode(i1),i1
       enddo
       
       close(124)
       
       if (CyclNo==1) stop
       
       if (CyclNo.ne.3) then
          call readjust_strctProp
       endif
       
    enddo
    
    !strctNo = 4
    !call destroy_and_rebuilt_for_cyclic_trnsfrmtn
    
    
    !call TwoStep_trnsfrmtn
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call get_Springwise_Areawise_and_Nodewise_forces
    
  end subroutine stage4_type1

  
  subroutine stage4_type2
    implicit none
    integer :: iMain
    integer :: strct_INPT,strct_TRGT
    real*8  :: hm
    integer :: i1,j,jmax
    
    N_struct = 3
    
    do iMain = (N_struct-1),1,-1
       
       strctNo = iMain
       write(*,*) strctNo
       
       call build_struct_And_saveProp_S4
       !write(*,*) A0(1),iMain,"A01,iMain"
       
       if (iMain.ne.1) then
          call destroy_struct
          call reinitialize_system_params
       endif
       
    enddo
    
    !write(*,*) l0_strct(1,1),l0_strct(2,1),"l0_strct"  
    !write(*,*) l0_strct(1,13),l0_strct(2,10:11),"l0_strct"
    !stop
    
    
    strct_INPT = 1 !
    strct_TRGT = 2 !
    
    call get_thirdStrctProp(strct_INPT,strct_TRGT) !
    !stop
    
    strctNo = 3
    call read_strctProps(strctNo)
    call Strct_to_actual(strctNo)
    
    strctNo = 4
    call destroy_and_rebuilt_for_cyclic_trnsfrmtn
    stop
    strctNo = 1 !!careful at this line
    
    EquilAlgrthm = 1
    
    do CyclNo = 1,3
       call TwoStep_trnsfrmtn
       
       open(unit=124,file='Varying_CgXNode.dat',position='append')
       
       do i1=1,N_node
          write(124,fmt=*) CgXNode(i1),i1
       enddo
       
       close(124)
       
       if (CyclNo==1) stop
       
       if (CyclNo.ne.3) then
          call readjust_strctProp
       endif
       
    enddo
    
    !strctNo = 4
    !call destroy_and_rebuilt_for_cyclic_trnsfrmtn
    
    
    !call TwoStep_trnsfrmtn
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call get_Springwise_Areawise_and_Nodewise_forces
    
  end subroutine stage4_type2
  
  
  
  subroutine change_shape_Initiation_wt_ExpFrmNo(hlfCyclNo,whichEnd)
    implicit none
    integer, intent(in) :: hlfCyclNo,whichEnd
    integer             :: ExpNo,FrmNo
    
    ExpNo = Exprmnt ; FrmNo = Frame
    call change_shape_of_InitiatnPhase(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    Exprmnt = ExpNo ; Frame = FrmNo
    
    ! This routine is noting but to combine the previous three lines so
    ! that the codes do not take too much space and it will be easy to read
  end subroutine change_shape_Initiation_wt_ExpFrmNo
  
  
end module Complt_Simltn
