module MAE_initaitn
  use Wfunc_and_its_derivative
  use changing_parameters
  use storing_changing_restoring_routines
  use switch_and_unswitch_models
  use nodeInsrtd_cntrlStates
  use Adding_cells
  
  implicit none
  integer              :: N_manpltdAngls
  real*8               :: ka_tobeCopied,A0_tobeCopied
  real*8               :: lowrF=0.20d0,medmF=0.25d0,higrF=0.50d0,zeroF=0.00d0
  
  real*8,  allocatable :: ks_tobeCopied(:),l0_tobeCopied(:)
  integer, allocatable :: cellNmbrs(:),crnrVal(:)
  real*8 , allocatable :: criticalAngls(:),cmprDirctn(:)
  real*8 , allocatable :: crnrAngls(:)
  integer, allocatable :: switchingXcoor(:),bfrSwitchXcoor(:)
  
  real*8, allocatable  :: NrmLa_Sim(:),NrmLa_Exp(:),diffNrm_La(:)
  real*8, allocatable  :: NrmLb_Sim(:),NrmLb_Exp(:),diffNrm_Lb(:)
  real*8, allocatable  :: NrmLl_Sim(:),NrmLl_Exp(:),diffNrm_Ll(:)
  real*8, allocatable  :: NrmAv_Sim(:),NrmAv_Exp(:),diffNrm_Av(:)
  
  real*8 :: simL1a,simL1b,simL1l,simA1
  real*8 :: tol_NrmL=0.05d0
  real*8 :: ZERO=0.00000000000000000d0
  
contains
  
  subroutine change_shape_of_InitiatnPhase(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: hlfCyclNo
    integer, intent(in)    :: whichEnd
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    !write(*,*) Exprmnt,Frame,"DO I KNOW U?"
    
    if (hlfCyclNo==1)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo) 
    if (hlfCyclNo==2)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo) 
    if (hlfCyclNo==3)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo) 
    if (hlfCyclNo==4)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo) 
    if (hlfCyclNo==5)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    if (hlfCyclNo==6)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    if (hlfCyclNo==7)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    if (hlfCyclNo==8)call combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    
  end subroutine change_shape_of_InitiatnPhase
  
  
  subroutine combine_diffrnt_manpltn_and_variabls(hlfCyclNo,whichEnd,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: hlfCyclNo
    integer, intent(in)    :: whichEnd
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    logical :: lgcl_AnglManplt=.False.
    logical :: lgcl_strght=.False. 
    integer :: cnt_Manplt
    integer :: CC
    integer :: writORread,manpltNum
    integer :: cellNm,cellNmToCmpr,ApBsLt,TnsnDirc,maxCnt,l0ORks
    real*8  :: hmuch,hmuchA0
    integer :: upToWhchCell,idealCell
    integer :: cell_A,ApBsLt_A,cell_B,ApBsLt_B
    integer :: diffRteOrNot
    real*8  :: limTnsn
    integer :: cntInRtn,deallctnNeeds
    integer :: dirctn
    
    
    if (hlfCyclNo==1) then
       
       if (whichEnd==-1) then
          call adjust_the_firstFrm_of_EpithelialLayer(ExpNo,FrmNo)
          writORread=1 ; manpltNum=38 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       elseif (whichEnd==+1) then
          
          call get_variabls_for_diff_hlfCycl(hlfCyclNo)
          
          !  *** Manplt=1 (start)**** 
          !cnt_Manplt = 0
          !do
          !  call manplt_angls_apcl_btwn_ltrl(ExpNo,FrmNo,lgcl_AnglManplt)
          ! cnt_Manplt = cnt_Manplt+1
          ! if (lgcl_AnglManplt .eqv. .True.) exit
          ! if (cnt_Manplt==10) exit
          !enddo
          
          !write(*,*) cnt_Manplt,"--- > HOW MANY TIMES THE manplt_angls_apcl_btwn_ltrl LOOP RUNS"
          !stop
          
          !  *** Manplt=1 (End) ***
          
          !  *** Manplt=2 (start)****
          !CC=1;call apcl_bsal_l0_manplt_of_cells_make_systm_strght(ExpNo,FrmNo,CC,lgcl_strght)
          !  *** Manplt=2 (End)****
          
          !  *** Manplt=3 (start)****
          CC=2;call apcl_bsal_l0_manplt_of_cells_make_systm_strght(ExpNo,FrmNo,CC,lgcl_strght)
          !  *** Manplt=3 (End)****
          !stop
       endif
       
    elseif (hlfCyclNo==2) then
       
       if (whichEnd==-1) then
          
          !call incrPressAndMaintainShapeOFInitiatorCell(ExpNo,FrmNo)
          !writORread=2 ; manpltNum=1 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1 ; ApBsLt=2 ; TnsnDirc=+1 ; maxCnt=40
          !call change_l0ksSpr_gvnCellNmAndEquill(cellNm,ApBsLt,TnsnDirc,maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=2 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1 ; ApBsLt=1 ; TnsnDirc=-1 ; maxCnt=5
          !call change_l0ksSpr_gvnCellNmAndEquill(cellNm,ApBsLt,TnsnDirc,maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=3 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !** Not using now
          
          !TnsnDirc=+1 ; maxCnt=20
          !call initiatr_cell_bsal_l0ksIncr(TnsnDirc,maxCnt,ExpNo,FrmNo)
          !writORread=1 ; manpltNum=4 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !** Not using now
          
          !cellNm=Hlf_Ncell; ApBsLt=1; l0ORks=2; TnsnDirc=-1; hmuch=0.75d0; maxCnt=1
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo) ! reducing the apcal mem ks for cellNm=Hlf_Ncell
          !writORread=2 ; manpltNum=5 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !TnsnDirc=-1 ; maxCnt=10
          !call initiatr_cell_TwoLtrl_l0ksIncr(TnsnDirc,maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=6 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-2; ApBsLt=2 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=10
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=7 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-2; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=+1 ; hmuch=0.10d0 ; maxCnt=10
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=8 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-2; ApBsLt=2 ; l0ORks=1 ; TnsnDirc=+1 ; hmuch=0.10d0 ; maxCnt=3
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=9 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1; ApBsLt=2 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=3
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=10 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          
          !cellNm=Hlf_Ncell-1; ApBsLt=2 ; l0ORks=2 ; TnsnDirc=-1 ; hmuch=0.50d0 ; maxCnt=5
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=11 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !upToWhchCell=Hlf_Ncell-2 ; idealCell=1
          !call making_InitiatnPhase_IncmingCellsIdntcl(upToWhchCell,idealCell,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=12 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=2
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=13 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cell_A=Hlf_Ncell-1 ; ApBsLt_A=1 ; cell_B=Hlf_Ncell-1 ; ApBsLt_B=2
          !call copy_prp_frm_sprATosprB(cell_A,ApBsLt_A,cell_B,ApBsLt_B,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=14 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1; ApBsLt=3 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=1
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=15 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          
          !write(*,*) "ENTERED AS IT IS"
          
          !A0(N_cell) = 1.10d0*A0(N_cell)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !A0(N_cell) = 1.10d0*A0(N_cell)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !l0(33) = 0.90d0*l0(33) ; l0(66) = 0.90d0*l0(66) ; l0(67) = 0.90d0*l0(67)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !l0(32) = 0.90d0*l0(32)  ; l0(65) = 0.90d0*l0(65)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !stop
          
          !cellNm=Hlf_Ncell ; ApBsLt=3 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.075d0 ; maxCnt=2
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          writORread=2 ; manpltNum=16 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !CgXNode(23) = 0.10d0 ; CgXNode(47) = 0.10d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !CgXNode(23) = 0.15d0 ; CgXNode(47) = 0.15d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !CgXNode(23) = 0.20d0 ; CgXNode(47) = 0.20d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !CgXNode(23) = 0.25d0 ; CgXNode(47) = 0.25d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !CgXNode(23) = 0.30d0 ; CgXNode(47) = 0.30d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !write(*,*) "stopping in ManipltngEquilibrateForInitiatn.f08"
          !stop
          
       elseif (whichEnd==+1) then
          
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !call incrPressAndMaintainShapeOFInitiatorCell(ExpNo,FrmNo)
          !writORread=2 ; manpltNum=17 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell; ApBsLt=2 ; l0ORks=3 ; TnsnDirc=+1 ; hmuch=0.10d0 ; maxCnt=15
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=18 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1; ApBsLt=2 ; l0ORks=3 ; TnsnDirc=+1 ; hmuch=0.10d0 ; maxCnt=15
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=19 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell; ApBsLt=1 ; l0ORks=3 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=5
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          !writORread=1 ; manpltNum=20 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !CgXNode(2)  = 0.15d0 ; CgXNode(26) = -0.15d0
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=2 ; manpltNum=23 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !hmuchA0=1.25d0; diffRteOrNot=1
          !call incrRigidNessOfInitiatorCell(hmuchA0,diffRteOrNot,ExpNo,FrmNo)
          !writORread=1 ; manpltNum=24 ; call savePropFrmPrvManplt(writORread,manpltNum)
          
          !cellNm=Hlf_Ncell ; dirctn=-1
          !call apply_ThreeDiff_force_at_TopNode_of_LS_ofAcell_andEquill(cellNm,dirctn,ExpNo,FrmNo) 
          !write(*,*) "stopping at ManipltngEquilibrateForInitiatn.f08 f"
          !stop
          
       endif
       
    elseif (hlfCyclNo==3) then
       
       if (whichEnd==-1) then
          continue
          
       elseif (whichEnd==+1) then
          
          !cellNm=N_cell ; cellNmToCmpr=Hlf_Ncell
          !call incrPressAndCmprToACell(cellNm,cellNmToCmpr,ExpNo,FrmNo)
          !writORread=1 ; manpltNum=25 ; call savePropFrmPrvManplt(writORread,manpltNum)
          
          !A0(N_cell) = 1.85d0*A0(N_cell)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          !writORread=1 ; manpltNum=26 ; call savePropFrmPrvManplt(writORread,manpltNum)
          
          
          !cellNm=Hlf_Ncell-1;ApBsLt=3;l0ORks=1 ;TnsnDirc=+1;hmuch=0.05d0;limTnsn=0.00d0
          !call change_l0KsorBothOfAspr_TochangeTnsn(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     limTnsn,ExpNo,FrmNo)
          !writORread=2 ; manpltNum=27 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-1; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.05d0 ; maxCnt=6
          !call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     maxCnt,ExpNo,FrmNo)
          writORread=2 ; manpltNum=28 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          cellNm=N_cell ; cntInRtn=1 ; deallctnNeeds=0
          call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          
       endif
       
    elseif (hlfCyclNo==4) then

       if (whichEnd==-1) then
          writORread=2 ; manpltNum=28 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       elseif (whichEnd==+1) then
          
          !cellNm=N_cell ; cntInRtn=2 ; deallctnNeeds=0
          !call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=2 ; manpltNum=29 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       endif
       
    elseif (hlfCyclNo==5) then
       
       
       if (whichEnd==-1) then
          call very_stiff_areas(ExpNo,FrmNo)
          write(*,*) "Stopped at ManpltingEquilibrateForInitiatn, hlCycl=5,whichEnd=-1"
          stop
          continue
       elseif (whichEnd==+1) then
          
          !cellNm=N_cell ; cntInRtn=2 ; deallctnNeeds=0
          !call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=2 ; manpltNum=30 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !cellNm=Hlf_Ncell-2;ApBsLt=3;l0ORks=1 ;TnsnDirc=+1;hmuch=0.05d0;limTnsn=-0.05d0
          !call change_l0KsorBothOfAspr_TochangeTnsn(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
          !     limTnsn,ExpNo,FrmNo)
          writORread=2 ; manpltNum=31 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          
          cellNm=Hlf_Ncell; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=2
          call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
               maxCnt,ExpNo,FrmNo)
          writORread=1 ; manpltNum=32 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          cellNm=Hlf_Ncell-2; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=-1 ; hmuch=0.10d0 ; maxCnt=7
          call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
               maxCnt,ExpNo,FrmNo)
          writORread=1 ; manpltNum=33 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          !stop
          
       endif
       
    elseif (hlfCyclNo==6) then
       
       if (whichEnd==-1) then
          writORread=2 ; manpltNum=33 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       elseif (whichEnd==+1) then
          
          cellNm=N_cell ; cntInRtn=2 ; deallctnNeeds=0
          call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=1 ; manpltNum=34 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       endif
       
    elseif (hlfCyclNo==7) then
       
       if (whichEnd==-1) then
          continue
       elseif (whichEnd==+1) then
          
          cellNm=N_cell ; cntInRtn=2 ; deallctnNeeds=0
          call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=1 ; manpltNum=35 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       endif
          
    elseif (hlfCyclNo==8) then
       
       if (whichEnd==-1) then
          writORread=2 ; manpltNum=35 ; call savePropFrmPrvManplt(writORread,manpltNum)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       elseif (whichEnd==+1) then
          
          cellNm=N_cell ; cntInRtn=2 ; deallctnNeeds=0
          call propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
          call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          writORread=1 ; manpltNum=36 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
          cellNm=Hlf_Ncell-1; ApBsLt=1 ; l0ORks=1 ; TnsnDirc=+1 ; hmuch=0.10d0 ; maxCnt=2
          call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
               maxCnt,ExpNo,FrmNo)
          writORread=1 ; manpltNum=37 ; call savePropFrmPrvManplt(writORread,manpltNum)
          !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
          
       endif
          
    endif
    
  end subroutine combine_diffrnt_manpltn_and_variabls
  
  
  
  
  subroutine get_variabls_for_diff_hlfCycl(hlfCyclNo)
    implicit none
    integer, intent(in) :: hlfCyclNo
    integer :: i,j
    integer :: N_itm
    integer :: whichSwitch
    
    if (hlfCyclNo == 1) then
       
       N_manpltdAngls = 9
       
       allocate(cellNmbrs(1:N_manpltdAngls))
       allocate(crnrVal(1:N_manpltdAngls))
       allocate(criticalAngls(1:N_manpltdAngls))
       allocate(cmprDirctn(1:N_manpltdAngls))
       allocate(crnrAngls(1:N_manpltdAngls))
       allocate(switchingXcoor(1:N_manpltdAngls),bfrSwitchXcoor(1:N_manpltdAngls))

       
       do i = 1,N_manpltdAngls
          cellNmbrs(i) = i
       enddo
       
       crnrVal(1:N_manpltdAngls)       = 3   !set as 3 as we need one cornr angle to compre
       criticalAngls(1:N_manpltdAngls) = 90.00d0
       cmprDirctn(1:N_manpltdAngls)    = +1       ! Either can be +1 or -1
       crnrAngls(1:N_manpltdAngls)     = +1000.0d0

       whichSwitch=1 ; call get_SwitchXcoorVals(whichSwitch)
       
       open(unit=214,file='get_variables_for_diff_hlfCycl.dat',position='append')
       
       N_itm = 6
       
       do i = 1,N_itm
          
          do j = 1,N_manpltdAngls
             
             if (i==1) write(214,*) cellNmbrs(j),j,"cellNmbrs"
             if (i==2) write(214,*) crnrVal(j),j,"crnrVals"
             if (i==3) write(214,*) criticalAngls(j),j,"criticalAngls"
             if (i==4) write(214,*) cmprDirctn(j),j,"cmprDirctns"
             if (i==5) write(214,*) crnrAngls(j),j,"crnrAngls"
             if (i==6) write(214,*) bfrSwitchXcoor(j),j,"bfrSwitch"
             
          enddo
          
          write(214,*) " "
       enddo
       
       close(214)
       
    elseif (hlfCyclNo == 2) then
       
       continue
    elseif (hlfCyclNo == 3) then
       
       continue
    endif
    
  end subroutine get_variabls_for_diff_hlfCycl
  
  
  
  
  subroutine incrRigidNessOfInitiatorCell(hmuchA0,diffRteOrNot,ExpNo,FrmNo)
    implicit none
    real*8 , intent(in)    :: hmuchA0
    integer, intent(in)    :: diffRteOrNot
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: ICB !ICB=InitiatorCellBottom
    real*8  :: areaBfrChng,areaAftChng
    real*8  :: prcntChngArea,tolPrcnt
    real*8  :: hmuchl0,hmuchks,hmuchl0_incr,hmuchks_incr
    
    integer :: sprsInICB,sprNm
    integer :: i,j
    integer :: cnt
    
    ICB         = N_cell
    areaBfrChng = A(ICB) !area here is ACTUALLY NI area
    
    A0(ICB) = (hmuchA0)*(A0(ICB))
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    call switchto_NI_model_run_and_switchbackto_TN
    
    hmuchl0 = 1.0d0-(abs(hmuchA0-1.0d0)/4.0d0)
    hmuchks = 1.0d0+(abs(hmuchA0-1.0d0)/4.0d0)

    hmuchl0_incr = 1.0d0-(3.0d0*abs(hmuchA0-1.0d0)/4.0d0)
    hmuchks_incr = 1.0d0+(3.0d0*abs(hmuchA0-1.0d0)/4.0d0)

    write(*,*) hmuchl0,hmuchks,"hmuchl0-hmuchks"
    
    if (hmuchl0.gt.1.0d0) then
       write(*,*) 'hmuchl0 cant be greater than 1.0d0'
       stop
    endif
    
    tolPrcnt = 1.0d0
    cnt = 1
    
    do
       sprsInICB = area_spr(ICB,0)
       
       do i = 1,sprsInICB
          sprNm = area_spr(ICB,i)
          
          if (diffRteOrNot==1) then
             if (i.ne.(sprsInICB-1)) then
                l0(sprNm)    = (hmuchl0)*l0(sprNm)
                k_spr(sprNm) = (hmuchks)*k_spr(sprNm)
             elseif (i==(sprsInICB-1)) then
                l0(sprNm)    = (hmuchl0_incr)*l0(sprNm)
                k_spr(sprNm) = (hmuchks_incr)*k_spr(sprNm)
             endif
          elseif (diffRteOrNot==0) then
             l0(sprNm)    = (hmuchl0)*l0(sprNm)
             k_spr(sprNm) = (hmuchks)*k_spr(sprNm)
          endif
          
          
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       call switchto_NI_model_run_and_switchbackto_TN
       areaAftChng = A(ICB)
       
       prcntChngArea=(abs(areaAftChng-areaBfrChng)/(areaBfrChng))*(100.0d0)
       if (prcntChngArea.le.tolPrcnt) then
          write(*,*) "areaGOTreduced and count =",cnt
          exit
       endif
       
       cnt=cnt+1
    enddo
    
  end subroutine incrRigidNessOfInitiatorCell

  
  subroutine incrPressAndCmprToACell(cellNm,cellNmToCmpr,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm,cellNmToCmpr,ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: hmuchA0
    integer :: diffRteOrNot,countLp
    integer :: extOrNot
    real*8  :: PressVal(1:N_cell)
    real*8  :: areaBfrLp
    integer :: sgmntdOrNot,whrTo
    interface
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8  :: Pressure(1:N_cell)
         real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8  :: dumA0(1:N_cell)
       end function Pressure
    end interface
    
    if (cellNm==N_cell) then
       
       countLp = 1
       areaBfrLp = A(cellNm)
       write(*,*) areaBfrLp,"areaBfrLp"

       
       
       do
          
          hmuchA0 = 1.06d0 ; diffRteOrNot=0
          call incrRigidNessOfInitiatorCell(hmuchA0,diffRteOrNot,ExpNo,FrmNo)
          
          sgmntdOrNot=1 ; whrTo=1
          call switchto_NI_model_run_switchbackto_TN_seprtlyOrtogthrNOFlnm(sgmntdOrNot,whrTo)
          PressVal(1:N_cell) = Pressure(node_xy,A0)
          if (PressVal(cellNm).gt.PressVal(cellNmToCmpr)) then
             write(*,*) "countLp is =",countLp
             exit
          endif
          
          write(*,*) PressVal(N_cell),PressVal(Hlf_Ncell),"Press of N_cell,Hlf_Ncell"
          write(*,*) A(N_cell),areaBfrLp,"area of N_cell after-before"
          
          sgmntdOrNot=1 ; whrTo=2
          call switchto_NI_model_run_switchbackto_TN_seprtlyOrtogthrNOFlnm(sgmntdOrNot,whrTo)
          countLp=countLp+1
       enddo
       
       write(*,*) "FrmNo and countLp =",(FrmNo-1),countLp
       
    endif
    
  end subroutine incrPressAndCmprToACell
  
  
  subroutine change_l0KsorBothOfAspr_TochangeTnsn(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
       limTnsn,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm,ApBsLt,l0ORks
    integer, intent(in)    :: TnsnDirc
    real*8 , intent(in)    :: hmuch,limTnsn
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: maxCnt,LpCntExt,sprNm ! will not be needed, just for calling purpose
    real*8  :: TnsnVal(1:N_spr)
    
    interface
       
       function TnsnComprsn(dum_Nodes,duml0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: TnsnComprsn(1:N_spr)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: duml0(1:N_spr)
       end function TnsnComprsn
       
    end interface
    
    maxCnt   = 1
    LpCntExt = 0
    
    do
       
       call change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
            maxCnt,ExpNo,FrmNo)
       TnsnVal(1:N_spr) = TnsnComprsn(node_xy,l0)

       if (ApBsLt==1) sprNm=(cellNm-1)*3+1
       if (ApBsLt==2) sprNm=(cellNm-1)*3+2
       if (ApBsLt==3) sprNm=(cellNm-1)*3+3
       
       write(*,*) TnsnVal(sprNm),limTnsn,TnsnDirc,"TnsnDirc"
       
       if (TnsnDirc==+1) then
          
          if (TnsnVal(sprNm) .le. limTnsn) then
             if (LpCntExt==1) then
                exit
             else
                LpCntExt = LpCntExt+1
             endif
          endif
          
       elseif (TnsnDirc==-1) then
          
          if (TnsnVal(sprNm) .gt. limTnsn) then
             if (LpCntExt==1) then ! just to Equilibrate 1 more time after matching condition
                exit
             else
                LpCntExt = LpCntExt+1
             endif
          endif
       endif
       
    enddo

    write(*,*) TnsnVal(sprNm),limTnsn,"TnsnS"
    call sleep(1)
    
  end subroutine change_l0KsorBothOfAspr_TochangeTnsn
  
  subroutine manplt_angls_apcl_btwn_ltrl(ExpNo,FrmNo,lgcl_AnglManplt)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    logical, intent(out)   :: lgcl_AnglManplt
    
    integer :: area_nmL,crnr_nmL,node_nmL,cmpr_val
    integer :: cnt_AnglChng
    real*8  :: Phi=-1000.0d0
    integer :: i,j
    integer :: switchV,choiceSpr
    
    lgcl_AnglManplt = .False.
    cnt_AnglChng    = 0
    
    open(unit=215,file='manplt_angls_apcl_btwn_ltrl.dat',position='append')
    
    do i = 1,N_manpltdAngls
       
       area_nmL = cellNmbrs(i)
       crnr_nmL = crnrVal(i)
       node_nmL = (area_nmL-1)*2 + crnr_nmL 
       cmpr_val = cmprDirctn(i)
       
       call get_Angle_AtANode(node_nmL,area_nmL,Phi) ; CrnrAngls(i) = Phi
       call get_SingleSwitchXcoorVals(i,switchV) ; switchingXcoor(i) = switchV
       
       write(215,*) area_nmL,crnr_nmL,node_nmL,switchingXcoor(i),bfrSwitchXcoor(i),i,"var values"
       
       
       if (switchingXcoor(i) .ne. bfrSwitchXcoor(i)) then
          cnt_AnglChng = cnt_AnglChng+1   
          
       elseif (switchingXcoor(i) == bfrSwitchXcoor(i)) then
          choiceSpr = 0 ! choice=0 as we wanna alter both apical and basal
          call manplt_apclbsalprop_forAngls(area_nmL,cmpr_val,choiceSpr)
       endif
       
    enddo
    
    write(215,*) cnt_AnglChng,"cnt_AnglChng"
    
    if (cnt_AnglChng == (N_manpltdAngls-7)) then
       write(*,*) "Angle manipulation is done"
       lgcl_AnglManplt = .True.
    elseif (cnt_AnglChng .ne. N_manpltdAngls) then
       continue
    endif
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    call switchto_NI_model_run_and_switchbackto_TN
    
    close(215)
    !stop
    
  end subroutine manplt_angls_apcl_btwn_ltrl
  
  
  subroutine manplt_apclbsalprop_forAngls(areaNmL,cmprVal,choice)
    implicit none
    integer, intent(in)  :: areaNmL,cmprVal
    integer, intent(in)  :: choice
    
    integer :: i,j,jmax
    integer :: nsprsInACell
    integer :: nApclSpr,nBsalSpr,nLtrlSpr
    integer :: sprNmL,sprNmR
    real*8  :: hmuchA,hmuchB
    
    real*8  :: l0BfrL,l0BfrR,l0AftL,l0AftR
    
    integer, allocatable :: apclSprL(:),bsalSprL(:),ltrlSprL(:)
    integer, allocatable :: apclSprR(:),bsalSprR(:),ltrlSprR(:)
    
    open(unit=216,file='manplt_apclbsalprop_forAngls.dat',position='append')
    
    if (modelID==1) then
       nApclSpr = 1 ; nBsalSpr=1 ; nLtrlSpr=1
       nsprsInACell = (napclSpr+nbsalSpr+nltrlSpr)
       
    elseif (modelID==2) then
       nApclSpr=NAEC_Apcl+1 ; nBsalSpr=NAEC_Bsal+1 ; nLtrlSpr=NAEC_Ltrl+1
       nsprsInACell = ((NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)) 
    endif
    
    if (cmprVal == +1) then
       hmuchA = -0.05d0
       hmuchB = +0.05d0
    elseif (cmprVal == -1) then
       hmuchA = +0.05d0
       hmuchB = -0.05d0
    endif
    
    allocate(apclSprL(1:napclSpr),bsalSprL(1:nbsalSpr),ltrlSprL(1:nltrlSpr))
    allocate(apclSprR(1:napclSpr),bsalSprR(1:nbsalSpr),ltrlSprR(1:nltrlSpr))
    
    
    do i = 1,2
       
       if (i==1) jmax = nApclSpr
       if (i==2) jmax = nBsalSpr
       
       do j = 1,jmax
          
          if (i==1) then
             
             apclSprL(j) = (areaNmL-1)*(nsprsInACell) + j
             apclSprR(j) = apclSprL(j) + (Hlf_Ncell*nsprsInACell)
             
             sprNmL = apclSprL(j) ; sprNmR = apclSprR(j)
             
             l0BfrL = l0(sprNmL) ; l0BfrR = l0(sprNmR)
             
             if (choice==0 .or. choice==1) then ! choice=0 for both, choice=1 for Apical only
                l0(sprNmL) = (1.0d0-hmuchA) * l0(sprNmL)
                l0(sprNmR) = l0(sprNmL)
             endif
             
             l0AftL = l0(sprNmL) ; l0AftR = l0(sprNmR)
             
             write(216,*) l0BfrL,l0AftL,(l0AftL/l0BfrL),sprNmL,"l0chng ApclL"
             write(216,*) l0BfrR,l0AftR,(l0AftR/l0BfrR),sprNmR,"l0chng ApclR"
             
          elseif (i==2) then
             
             bsalSprL(j) = (areaNmL-1)*(nsprsInACell) + nApclSpr + j
             bsalSprR(j) = bsalSprL(j) + (Hlf_Ncell*nsprsInACell)
             
             sprNmL = bsalSprL(j) ; sprNmR = bsalSprR(j)
             
             l0BfrL = l0(sprNmL) ; l0BfrR = l0(sprNmR)
             
             if (choice==0 .or. choice==2) then ! choice=0 for both, choice=2 for Basal only
                l0(sprNmL) = (1.0d0-hmuchB) * l0(sprNmL)
                l0(sprNmR) = l0(sprNmL)
             endif
             
             l0AftL = l0(sprNmL) ; l0AftR = l0(sprNmR)
             
             write(216,*) l0BfrL,l0AftL,(l0AftL/l0BfrL),sprNmL,"l0chng BsalL"
             write(216,*) l0BfrR,l0AftR,(l0AftR/l0BfrR),sprNmR,"l0chng BsalR"
             
          endif
          
       enddo
       
       write(216,*) " "
       
    enddo
    
    write(216,*) "End OF the ROUTINE"
    
    close(216)
    
  end subroutine manplt_apclbsalprop_forAngls
  
  
  subroutine manpltns_for_CF_WT_Nrt_P_2G5C_m()
    implicit none
    integer :: cell_val
    real*8  :: hmL0Ap1,hmL0Ap2,hmL0Ap3
    
    integer, allocatable   :: ApclSp1(:), BsalSp1(:), LtrlSp1(:)
    integer, allocatable   :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    integer, allocatable   :: ApclSp3(:), BsalSp3(:), LtrlSp3(:)
    
    if (VF_regionModelled .ne. 1) stop 'for VF only'
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),  BsalSp1(1:(NAEC_Bsal+1)),  LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),  BsalSp2(1:(NAEC_Bsal+1)),  LtrlSp2(1:(NAEC_Ltrl+1)))
    allocate(ApclSp3(1:(NAEC_Apcl+1)),  BsalSp3(1:(NAEC_Bsal+1)),  LtrlSp3(1:(NAEC_Ltrl+1)))
    
    call get_ApclBsal_singl_cell_NIsys_with_CM0_and_VFarea(ApclSp1,BsalSp1)
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    cell_val=Hlf_Ncell-1  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    
    hmL0Ap1=0.05d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_Apcl,hmL0Ap1)
    hmL0Ap2=0.20d0 ; call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap2)
    hmL0Ap3=1.20d0 ; call manplt_spr_NIsys(ApclSp3,NAEC_Apcl,hmL0Ap3)
    
    call Equilibrate_only_NI_model_withVF_region
    
  end subroutine manpltns_for_CF_WT_Nrt_P_2G5C_m

  
  subroutine manpltns_for_CF_WT_Nrt_P_1G1B1_m()
    implicit none
    integer                :: cell_val
    real*8                 :: hmL0Ap1,hmL0Bs1,hmL0Ap2,hmL0Bs2
    integer, allocatable   :: ApclSp1(:), BsalSp1(:)
    integer, allocatable   :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    
    if (VF_regionModelled .ne. 1) stop 'for VF only'
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:(NAEC_Bsal+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:(NAEC_Bsal+1)),LtrlSp2(1:(NAEC_Ltrl+1)))
    
    call get_ApclBsal_singl_cell_NIsys_with_CM0_and_VFarea(ApclSp1,BsalSp1)
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    
    hmL0Ap1=0.50d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_Apcl,hmL0Ap1)
    hmL0Bs1=1.20d0 ; call manplt_spr_of_singl_cell_NIsys(BsalSp1,NAEC_Bsal,hmL0Bs1)
    hmL0Ap2=1.30d0 ; call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap2)
    hmL0Bs2=0.90d0 ; call manplt_spr_NIsys(BsalSp2,NAEC_Bsal,hmL0Bs2)
    
    call Equilibrate_only_NI_model_withVF_region
    
  end subroutine manpltns_for_CF_WT_Nrt_P_1G1B1_m
  
  
  subroutine manplts_for_JHT_articl_Img_1M
    implicit none
    real*8                 :: Fctr
    integer                :: i,j,imax,jmax
    integer                :: N_diffForces,numLoop
    integer                :: cellNm,NAEC_ApclV1,NAEC_ApclV2
    real*8                 :: hmL0Ap1,  hmL0Bs1,  hmL0Ap2,  hmL0Bs2,  hmL0Lt2
    real*8                 :: hmL0Ap1ES,hmL0Bs1ES,hmL0Ap2ES,hmL0Bs2ES,hmL0Lt2ES,hmA0ES ! ES = Each Step
    
    integer, allocatable   :: ApclSp1(:), BsalSp1(:)
    integer, allocatable   :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    
    if (VF_regionModelled  .ne. 1) stop 'VF region actvtn needed'
    if (NI_AS_wtCrtclSurfc .ne. 1) stop 'nodes In Apcl Side with Crtcl Surfc is needed'
    
    Frame_NIVF=1
    N_diffForces=3 ; call apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion(N_diffForces)
    Frame_NIVF=21  ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    call Equilibrate_only_NI_model_withVF_region
    
    !Fctr=3.00d0     ; call Effct_of_VF_manpltng_ka(Fctr) 
    Frame_NIVF=27 ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    call Equilibrate_only_NI_model_withVF_region
    
    call get_ApclBsal_singl_cell_NIsys_with_CM0_VFarea_and_CrtclSrfc(ApclSp1,BsalSp1,NAEC_ApclV1)
    cellNm=Hlf_Ncell; call get_ApclBsalLtrl_OfNIsystmCrtclApclSrfc(cellNm,ApclSp2,BsalSp2,LtrlSp2,NAEC_ApclV2)
    
    write(*,*) ApclSp1,"Ap1" ; write(*,*) BsalSp1,"Bs1"
    write(*,*) ApclSp2,"Ap2" ; write(*,*) BsalSp2,"Bs2" ; write(*,*) LtrlSp2,"Lt2"
    
    !numLoop=1 ; hmL0Lt2ES=-0.10d0 ; call manplt_spr_NIsys_inLoop(numLoop,LtrlSp2,NAEC_Ltrl,hmL0Lt2ES)
    
    !hmL0Ap1=0.60d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_ApclV1,hmL0Ap1)
    !call Equilibrate_only_NI_model_withVF_region
    !hmL0Lt2=0.95d0 ; call manplt_spr_NIsys_withCrtclApSrfc(LtrlSp2,NAEC_Ltrl,hmL0Lt2)
    !call Equilibrate_only_NI_model_withVF_region
    !hmL0Lt2=0.95d0 ; call manplt_spr_NIsys_withCrtclApSrfc(LtrlSp2,NAEC_Ltrl,hmL0Lt2)
    !call Equilibrate_only_NI_model_withVF_region
    
    Frame_NIVF=31 ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    call Equilibrate_only_NI_model_withVF_region
    
    !call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    Frame_NIVF=42 ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    call Equilibrate_only_NI_model_withVF_region
    
    !hmL0Ap1=0.60d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_ApclV1,hmL0Ap1)
    !call Equilibrate_only_NI_model_withVF_region
    
    Frame_NIVF=43 ; call read_config_and_start_simlnFrm_there(Exprmnt_NIVF,Frame_NIVF)
    call Equilibrate_only_NI_model_withVF_region
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
  end subroutine manplts_for_JHT_articl_Img_1M
  !numLoop=1 ; cellNm=N_cell-1   ; hmA0ES =-0.02d0 ; call manplt_press_withVF_region(numLoop,cellNm,hmA0ES)
  
  subroutine Diffrnrt_Level_Forces_AtFrstJointNode
    implicit none
    integer                :: N_diffForces,nodeL,cntLp,bfrORaftFrceApply
    integer                :: cellNm,NAEC_ApclV1,NAEC_ApclV2,Dirctnlity
    real*8                 :: hmL0Ap1,  hmL0Bs1,  hmL0Ap2,  hmL0Bs2,  hmL0Lt2
    real*8                 :: hmL0Ap1ES,hmL0Bs1ES,hmL0Ap2ES,hmL0Bs2ES,hmL0Lt2ES,hmA0ES ! ES = Each Step
    real*8                 :: bsalSrfcMvmnt,bsalSrfcBoundry,bsalSrfcMvmnt_aftForceApp
    integer, allocatable   :: ApclSp1(:), BsalSp1(:)
    integer, allocatable   :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    
    if (VF_regionModelled  .ne. 1) stop 'VF region actvtn needed'
    if (NI_AS_wtCrtclSurfc .ne. 1) stop 'nodes In Apcl Side with Crtcl Surfc is needed'
    
    call get_ApclBsal_singl_cell_NIsys_with_CM0_VFarea_and_CrtclSrfc(ApclSp1,BsalSp1,NAEC_ApclV1)
    cellNm=Hlf_Ncell; call get_ApclBsalLtrl_OfNIsystmCrtclApclSrfc(cellNm,ApclSp2,BsalSp2,LtrlSp2,NAEC_ApclV2)
    write(*,*) ApclSp1,"Ap1" ; write(*,*) BsalSp1,"Bs1"
    write(*,*) ApclSp2,"Ap2" ; write(*,*) BsalSp2,"Bs2" ; write(*,*) LtrlSp2,"Lt2"
    
    Frame_NIVF      = 1
    bsalSrfcBoundry = node_xy(2,2)
    
    nodeL = (2*Hlf_Ncell+1) ; CgYNode(nodeL)=0.0000d0
    N_diffForces=21 ; call apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion_frmCurrntVal(N_diffForces)
    
    call get_bsalSrfcMvmntInfo(bsalSrfcBoundry,bsalSrfcMvmnt,Dirctnlity)
    bsalSrfcMvmnt_aftForceApp = bsalSrfcMvmnt ;write(*,*)bsalSrfcMvmnt,Dirctnlity,"bsalSrfcMvmnt aft ForcApply"
    
    cntLp = 1
    do
       hmL0Lt2=0.95d0 ; call manplt_spr_NIsys_withCrtclApSrfc(LtrlSp2,NAEC_Ltrl,hmL0Lt2)
       call Equilibrate_only_NI_model_withVF_region
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call get_bsalSrfcMvmntInfo(bsalSrfcBoundry,bsalSrfcMvmnt,Dirctnlity)
       write(*,*) bsalSrfcMvmnt,bsalSrfcMvmnt_aftForceApp,cntLp,Dirctnlity,"bsalSrfc Movement"
       
       if (Dirctnlity==+1 .or. cntLp==10) then
          write(*,*) Dirctnlity,cntLp,"Dir-Cnt" ; exit
       endif
       
       cntLp = cntLp+1
    enddo
    
  end subroutine Diffrnrt_Level_Forces_AtFrstJointNode
  
  subroutine get_bsalSrfcMvmntInfo(bsalSrfcBoundry,bsalSrfcMvmnt,Dirctnlity)
    implicit none
    real*8,  intent(in)  :: bsalSrfcBoundry
    real*8,  intent(out) :: bsalSrfcMvmnt
    integer, intent(out) :: Dirctnlity
    integer              :: bsalICNodes(1:NAEC_Bsal),nodesBfr,i,j
    real*8               :: yNodebsalIC(1:NAEC_Bsal),yNodesSum,yNodesAvrg
    
    nodesBfr = numRegulrNode + (2*(Hlf_Ncell-NCP_CrtclApSrfc))*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + &
         (2*NCP_CrtclApSrfc)*(NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl) + (NAEC_ApclCrtcl)
    
    yNodesSum = 0.0000d0
    
    do j = 1,NAEC_Bsal
       bsalICNodes(j) = nodesBfr+j ; write(*,*) bsalICNodes(j),j,"BsalIC"
       yNodebsalIC(j) = node_xy(bsalICNodes(j),2)
       yNodesSum      = yNodesSum + yNodebsalIC(j) ; write(*,*) yNodesSum,yNodebsalIC(j),"yN,yS"
    enddo
    
    yNodesAvrg    = yNodesSum/(NAEC_Bsal*1.0000d0) 
    bsalSrfcMvmnt = abs(bsalSrfcBoundry-yNodesAvrg) ; write(*,*)bsalSrfcBoundry,yNodesAvrg,bsalSrfcMvmnt,"srfc"
    
    if (yNodesAvrg .le. bsalSrfcBoundry) Dirctnlity=-1
    if (yNodesAvrg .gt. bsalSrfcBoundry) Dirctnlity=+1
    
  end subroutine get_bsalSrfcMvmntInfo
  
  subroutine readFrm_FrceDownLtrlNC1_manpltn_withDiffrntPrtrbtn()
    implicit none
    character(len=200)    :: flplce1,flplce2,flplce
    integer               :: ExpNo,FrmNoToBeRead,cellNm,NAEC_ApclV1,NAEC_ApclV2
    integer, allocatable  :: ApclSp1(:), BsalSp1(:)
    integer, allocatable  :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    real*8                :: hmL0Ap1,hmL0Bs1,hmL0Ap2,hmL0Bs2,hmL0Lt2,Fctr,hmL0ES
    real*8                :: l0Str(1:N_spr),ksStr(1:N_spr),A0Str(1:N_cell),kaStr(1:N_cell)
    real*8                :: CgXnStr(1:N_node),CgYnStr(1:N_node)
    integer               :: i,j,numLoop
    
    if (VF_regionModelled  .ne. 1) stop 'VF region actvtn needed'
    if (NI_AS_wtCrtclSurfc .ne. 1) stop 'nodes In Apcl Side with Crtcl Surfc is needed'
    
    call get_ApclBsal_singl_cell_NIsys_with_CM0_VFarea_and_CrtclSrfc(ApclSp1,BsalSp1,NAEC_ApclV1)
    cellNm=Hlf_Ncell ; call get_ApclBsalLtrl_OfNIsystmCrtclApclSrfc(cellNm,ApclSp2,BsalSp2,LtrlSp2,NAEC_ApclV2)
    
    write(flplce1,*)"/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure"
    write(flplce2,*)"/Force_LtrlSprNC1_comboTst_CgY=1.00/OrgnlData_With_ForceDown_And_LtrlNC1Shortn/"
    flplce=trim(adjustl(flplce1))//trim(adjustl(flplce2))
    
    ExpNo = Exprmnt_NIVF ; Frame_NIVF = 21 ; FrmNoToBeRead = Frame_NIVF
    call read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
    call Equilibrate_only_NI_model_withVF_region
    !call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    l0Str=l0 ; ksStr=k_spr ; A0Str=A0 ; kaStr=k_area ; CgXnStr=CgXNode ; CgYnStr=CgYNode

    numLoop=4 ; hmL0ES = -0.1000d0
    call manplt_NC1Lt_withBndryPull_with_BndryTiltAdjstmnt(numLoop,LtrlSp2,NAEC_Ltrl,hmL0ES)
    stop 'aft Force Down + ltrl shorten + boundary pulling + always boundary force adjustmnt tilt'
    
    numLoop=6 ; hmL0ES = -0.0500d0
    call manplt_spr_NIsys_inLoop_with_BndryTiltAdjstmnt(numLoop,LtrlSp2,NAEC_Ltrl,hmL0ES)
    stop 'aft Force Down + ltrl shorten + always boundary force adjustmnt to reduce tilting after each change'
    
    numLoop=6 ; hmL0ES = -0.0500d0
    call manplt_spr_NIsys_inLoop(numLoop,LtrlSp2,NAEC_Ltrl,hmL0ES)
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    stop 'Aft Force Down Ltrl Shorten Then Bndry Adjstmnt'
    
    hmL0Ap1=0.60d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_ApclV1,hmL0Ap1)
    call Equilibrate_only_NI_model_withVF_region
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    l0=l0Str ; k_spr=ksStr ; A0=A0Str ; k_area=kaStr ; CgXNode=CgXnStr ; CgYNode=CgYnStr
    call Equilibrate_only_NI_model_withVF_region ; write(*,*) Frame_NIVF-1,"aft restore NIVF"
    Fctr=3.00d0    ; call Effct_of_VF_manpltng_ka(Fctr)
    hmL0Ap1=0.60d0 ; call manplt_spr_of_singl_cell_NIsys(ApclSp1,NAEC_ApclV1,hmL0Ap1)
    call Equilibrate_only_NI_model_withVF_region
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
  end subroutine ReadFrm_FrceDownLtrlNC1_Manpltn_WithDiffrntPrtrbtn
  
  
  subroutine Effct_of_VF_manpltng_ka(Fctr)
    implicit none
    real*8, intent(in) :: Fctr
    real*8             :: kaI,kaF,stepF
    integer            :: Nstp
    real*8             :: NstpReal
    
    kaI   = k_area(N_cell) ; kaF=(Fctr)*kaI ; Nstp=5 ; NstpReal=5.00d0
    stepF = (kaF-kaI)/(NstpReal)
    
    do i = 0,Nstp
       k_area(N_cell) = kaI+(real(i))*stepF ; write(*,*) k_area(N_cell),kaI,"ka"
       call Equilibrate_only_NI_model_withVF_region
    enddo
    
  end subroutine Effct_of_VF_manpltng_ka

  
  subroutine manplt_spr_NIsys_inLoop(numLoop,Sprs,NAEC_val,hmL0ES)
    implicit none
    integer, intent(in) :: numLoop
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0ES
    real*8              :: l0ValStr(1:(NAEC_val+1))
    integer             :: i,j,imax,jmax
    real*8              :: hmL0
    integer             :: strORrstore
    
    l0ValStr(1:(NAEC_val+1)) = -1.0d30
    strORrstore=1 ; call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore) ! storing
    
    do i = 1,numLoop
       hmL0 = 1.0000d0 + (i*1.00d0)*(hmL0ES)
       if (NI_AS_wtCrtclSurfc==0) call manplt_spr_NIsys(Sprs,NAEC_val,hmL0)
       if (NI_AS_wtCrtclSurfc==1) call manplt_spr_NIsys_withCrtclApSrfc(Sprs,NAEC_val,hmL0)
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
       strORrstore=2
       if (i.ne.numLoop) call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore)
    enddo
    
  end subroutine manplt_spr_NIsys_inLoop
  
  
  subroutine manplt_spr_NIsys_inLoop_with_BndryTiltAdjstmnt(numLoop,Sprs,NAEC_val,hmL0ES)
    implicit none
    integer, intent(in) :: numLoop
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0ES
    real*8              :: l0ValStr(1:(NAEC_val+1))
    integer             :: i,j,imax,jmax
    real*8              :: hmL0
    integer             :: strORrstore
    
    l0ValStr(1:(NAEC_val+1)) = -1.0d30
    strORrstore=1 ; call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore) ! storing
    
    do i = 1,numLoop
       hmL0 = 1.0000d0 + (i*1.00d0)*(hmL0ES)
       
       if (NI_AS_wtCrtclSurfc==0) call manplt_spr_NIsys(Sprs,NAEC_val,hmL0)
       if (NI_AS_wtCrtclSurfc==1) call manplt_spr_NIsys_withCrtclApSrfc(Sprs,NAEC_val,hmL0)
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       
       strORrstore=2
       if (i.ne.numLoop) call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore)
    enddo
    
  end subroutine manplt_spr_NIsys_inLoop_with_BndryTiltAdjstmnt
  
  subroutine manplt_NC1Lt_withBndryPull_with_BndryTiltAdjstmnt(numLoop,Sprs,NAEC_val,hmL0ES)
    implicit none
    integer, intent(in) :: numLoop
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0ES
    real*8              :: l0ValStr(1:(NAEC_val+1))
    integer             :: i,j,imax,jmax,pullingOrpushing
    real*8              :: hmL0,hmuch
    integer             :: strORrstore
    
    l0ValStr(1:(NAEC_val+1)) = -1.0d30
    strORrstore=1 ; call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore)!storing
    pullingOrpushing = 1
    
    do i = 1,numLoop
       hmL0 = 1.0000d0 + ((i-1)*1.00d0)*(hmL0ES)
       
       if (NI_AS_wtCrtclSurfc==0) call manplt_spr_NIsys(Sprs,NAEC_val,hmL0)
       if (NI_AS_wtCrtclSurfc==1) call manplt_spr_NIsys_withCrtclApSrfc(Sprs,NAEC_val,hmL0)
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region

       jmax = 5
       do j = 1,jmax
          
          if (j==1)   hmuch = 0.0000d0
          if (j.gt.1) hmuch = 0.2500d0
          if (j==1) write(*,*) Frame_NIVF,"at beg"
          
          pullingOrpushing = 1    ; call pulling_or_pushing_frmCurrntVal(pullingOrpushing,hmuch)
          if (VF_regionModelled==0) call Equilibrate_only_NI_model
          if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
          
          call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
          
          if ((j==jmax) .and. (i.ne.numLoop)) then
             pullingOrpushing = 2 ; hmuch = ((jmax-1)*0.2500d0)
             call pulling_or_pushing_frmCurrntVal(pullingOrpushing,hmuch)
             if (VF_regionModelled==0) call Equilibrate_only_NI_model
             if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
             call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
          endif
          
       enddo
       
       strORrstore=2
       if (i.ne.numLoop) call str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore)
    enddo
    
  end subroutine manplt_NC1Lt_withBndryPull_with_BndryTiltAdjstmnt
  
  subroutine manplt_press_withVF_region(numLoop,cellNm,hmA0ES)
    implicit none
    integer, intent(in) :: numLoop
    integer, intent(in) :: cellNm
    real*8 , intent(in) :: hmA0ES
    real*8              :: A0ValStr
    integer             :: i,j,imax,jmax
    real*8              :: hmA0
    integer             :: cellNmL,cellNmR
    
    if (VF_regionModelled==0) then
       
       if (cellNm.ne.N_cell) then
          cellNmL = cellNm
          cellNmR = cellNm + ((N_cell-1)/2)
       elseif (cellNm==N_cell) then
          cellNmL = cellNm
          cellNmR = cellNm
       endif
       
    elseif (VF_regionModelled==1) then
       
       if (cellNm.le.(N_cell-2)) then
          cellNmL = cellNm
          cellNmR = cellNm + ((N_cell-1)/2)
          
       elseif (cellNm.gt.(N_cell-2)) then
          cellNmL = cellNm
          cellNmR = cellNm
       endif
       
    endif
    
    A0ValStr = A0(cellNmL) ; write(*,*) A0ValStr,"A0ValStr" ; write(*,*) cellNmL,cellNmR,"cellNm L/R"
    write(*,*) cellNm,hmA0ES,N_cell,VF_regionModelled,"chk Prp"
    
    do i = 1,numLoop
       hmA0 = 1.0000d0 + (i*1.00d0)*(hmA0ES) ; write(*,*) hmA0,hmA0ES,"Lp1"
       A0(cellNmL) = (hmA0*A0(cellNmL)) ; A0(cellNmR) = A0(cellNmL) ; write(*,*) A0(cellNmL),A0(cellNmR),"Lp2"
       call Equilibrate_only_NI_model_withVF_region
       if (i.ne.numLoop) A0(cellNmL) = A0ValStr
    enddo
    
  end subroutine manplt_press_withVF_region
  
  subroutine get_the_apcls_bsals_ltrls(NcellToFind,cellsToFind,ApclsL,BsalsL,LtrlsL)
    implicit none
    integer, intent(in)  :: NcellToFind
    integer, intent(in)  :: cellsToFind(1:NcellToFind)
    integer, intent(out) :: ApclsL(1:NcellToFind,1:(NAEC_Apcl+1))
    integer, intent(out) :: BsalsL(1:NcellToFind,1:(NAEC_Bsal+1))
    integer, intent(out) :: LtrlsL(1:NcellToFind,1:(NAEC_Ltrl+1))
    
    integer :: i1,i2,i3,i3max
    integer :: cellNm,nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i1 = 1,NcellToFind
       
       cellNm = cellsToFind(i1) ; write(*,*) cellNm,"cellNm for i1=",i1
       
       do i2 = 1,3 ! Apcl/Basal/Ltrl
          
          if (i2==1) i3max=NAEC_Apcl+1
          if (i2==2) i3max=NAEC_Bsal+1
          if (i2==3) i3max=NAEC_Ltrl+1
          
          do i3 = 1,i3max
             
             if (i2==1) ApclsL(i1,i3) = (cellNm-1)*(nsprsInACell) + i3
             if (i2==2) BsalsL(i1,i3) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + i3
             if (i2==3) LtrlsL(i1,i3) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+NAEC_Bsal+2) + i3
             
          enddo
          
          if (i2==1) write(*,*) ApclsL(i1,1:(NAEC_Apcl+1)),"ApclsL"
          if (i2==2) write(*,*) BsalsL(i1,1:(NAEC_Bsal+1)),"BsalsL"
          if (i2==3) write(*,*) LtrlsL(i1,1:(NAEC_Ltrl+1)),"LtrlsL"
          
       enddo
       
    enddo
    
  end subroutine get_the_apcls_bsals_ltrls
  
  subroutine get_only_apcls_or_bsals_or_ltrls(N_cellsInvlvd,cellNms,whichOne,NAEC_val,SprNms)
    implicit none
    integer, intent(in)  :: N_cellsInvlvd
    integer, intent(in)  :: cellNms(1:N_cellsInvlvd)
    integer, intent(in)  :: whichOne,NAEC_val
    integer, intent(out) :: SprNms(1:N_cellsInvlvd,1:(NAEC_val+1)) 
    
    integer :: i,j
    integer :: cellNm
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i = 1,N_cellsInvlvd
       
       cellNm = cellNms(i)
       
       do j = 1,(NAEC_val+1)

          if (whichOne == 1) then
             SprNms(i,j) = (cellNm-1)*(nsprsInACell) + j
          elseif (whichOne == 2) then
             SprNms(i,j) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
          elseif (whichOne == 3) then
             SprNms(i,j) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+NAEC_Bsal+2) + j
          endif
          
       enddo
       
    enddo
    
  end subroutine get_only_apcls_or_bsals_or_ltrls
  
  
  subroutine get_the_apcls_bsals_and_both_ltrls(NcellToFind,cellsToFind,ApclsL,BsalsL,LtrlsBL,LtrlsAL)
    implicit none
    integer, intent(in)  :: NcellToFind
    integer, intent(in)  :: cellsToFind(1:NcellToFind)
    integer, intent(out) :: ApclsL(1:NcellToFind,1:(NAEC_Apcl+1))
    integer, intent(out) :: BsalsL(1:NcellToFind,1:(NAEC_Bsal+1))
    integer, intent(out) :: LtrlsBL(1:NcellToFind,1:(NAEC_Ltrl+1))
    integer, intent(out) :: LtrlsAL(1:NcellToFind,1:(NAEC_Ltrl+1))
    
    integer :: i1,i2,i3,i3max
    integer :: cellNm,nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i1 = 1,NcellToFind
       
       cellNm = cellsToFind(i1) ; write(*,*) cellNm,"cellNm for i1=",i1
       
       do i2 = 1,4! Apcl/Basal/Ltrl1/Ltrl2
          
          if (i2==1) i3max=NAEC_Apcl+1
          if (i2==2) i3max=NAEC_Bsal+1
          if (i2==3) i3max=NAEC_Ltrl+1
          if (i2==4) i3max=NAEC_Ltrl+1
          
          do i3 = 1,i3max
             
             if (i2==1) ApclsL(i1,i3)  = (cellNm-1)*(nsprsInACell) + i3
             if (i2==2) BsalsL(i1,i3)  = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + i3
             if (i2==3) LtrlsBL(i1,i3) = (cellNm-2)*(nsprsInACell) + (NAEC_Apcl+NAEC_Bsal+2) + i3
             if (i2==4) LtrlsAL(i1,i3) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+NAEC_Bsal+2) + i3
             
          enddo
          
          if (i2==1) write(*,*) ApclsL(i1,1:(NAEC_Apcl+1)), "ApclsL"
          if (i2==2) write(*,*) BsalsL(i1,1:(NAEC_Bsal+1)), "BsalsL"
          if (i2==3) write(*,*) LtrlsBL(i1,1:(NAEC_Ltrl+1)),"LtrlsBL"
          if (i2==4) write(*,*) LtrlsAL(i1,1:(NAEC_Ltrl+1)),"LtrlsAL"
          
       enddo
       
    enddo
    
  end subroutine get_the_apcls_bsals_and_both_ltrls
  
  
  subroutine get_SwitchXcoorVals(whichOne)
    implicit none
    integer, intent(in) :: whichOne
    
    integer :: i,j
    integer :: nodeAtTop,nodeAtBot
    real*8  :: xValAtTop,xvalAtBot

    open(unit=217,file='get_SwitchXcoor.dat',position='append')
    
    do i = 1,N_manpltdAngls
       
       nodeAtTop = crnrVal(i) + (i-1)*2
       nodeAtBot = nodeAtTop+1
       xvalAtTop = abs(node_xy(nodeAtTop,1))
       xvalAtBot = abs(node_xy(nodeAtBot,1))
       
       write(217,*) nodeAtTop,nodeAtBot,xvalAtTop,xvalAtBot,i,"node-x"
       
       if (whichOne==1) then
          
          if (xvalAtTop .gt. xvalAtBot) bfrSwitchXcoor(i) = +1
          if (xvalAtTop .le. xvalAtBot) bfrSwitchXcoor(i) = -1
          
          write(217,*) bfrSwitchXcoor(i),i,"bfrSwitch"
          
       elseif (whichOne==2) then
          
          if (xvalAtTop .gt. xvalAtBot) switchingXcoor(i) = +1 
          if (xvalAtTop .le. xvalAtBot) switchingXcoor(i) = -1
          
          write(217,*) switchingXcoor(i),i,"SwitchXcoor"
          
       endif
       
    enddo
    
    close(217)
    
  end subroutine get_SwitchXcoorVals
  
  subroutine get_SingleSwitchXcoorVals(seq,switchV)
    implicit none
    integer, intent(in)  :: seq
    integer, intent(out) :: switchV
    
    integer :: nodeAtTop,nodeAtBot
    real*8  :: xValAtTop,xvalAtBot

    open(unit=218,file='SingleSwitch1.dat',position='append')
    open(unit=219,file='SingleSwitch2.dat',position='append')
    
    nodeAtTop = crnrVal(seq) + (seq-1)*2
    nodeAtBot = nodeAtTop+1
    xvalAtTop = abs(node_xy(nodeAtTop,1))
    xvalAtBot = abs(node_xy(nodeAtBot,1))
    
    write(218,*) nodeAtTop,nodeAtBot,xvalAtTop,xvalAtBot,i,"node-x"
    
    if (xvalAtTop .gt. xvalAtBot) switchV = +1
    if (xvalAtTop .le. xvalAtBot) switchV = -1
    
    write(219,*) seq,switchV,"seq-switchV"
    
    close(218)
    close(219)
    
  end subroutine get_SingleSwitchXcoorVals
  
  
  subroutine apcl_bsal_l0_manplt_of_cells_make_systm_strght(ExpNo,FrmNo,CC,lgcl_strght)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    integer, intent(in)    :: CC    ! CC=Choice of Cells
    logical, intent(out)   :: lgcl_strght
    
    integer, allocatable :: cellNm(:)
    integer, allocatable :: l0_manpltDirctn(:)
    integer, allocatable :: choicSpr(:)
    
    integer :: choicSprV
    integer :: NcellMan
    integer :: cellNmV,l0_manpltDirctnV,choiceSprV
    
    if (CC==1) then
       
       NcellMan=2
       allocate(cellNm(1:NcellMan),l0_manpltDirctn(1:NcellMan),choicSpr(1:NcellMan))
       
       cellNm(1)                   = hlf_Ncell-1
       cellNm(2)                   = hlf_Ncell
       l0_manpltDirctn(1:NcellMan) = -1 ! (-1) --- Apcl contrct + Bsal expnsn
       choicSpr(1:NcellMan)        = 0  ! both apical and basal manplt
       
    elseif (CC==2) then
       
       NcellMan  = 1
       allocate(cellNm(1:NcellMan),l0_manpltDirctn(1:NcellMan),choicSpr(1:NcellMan))
       
       cellNm(1)                   = hlf_Ncell-2
       l0_manpltDirctn(1:NcellMan) = +1 ! (+1) --- Apcl expnsn + Bsal contrct
       choicSpr(1:NcellMan)        =  1 ! only apical   
       
    endif
    
    do
       
       do i = 1,NcellMan
          cellNmV          = cellNm(i)
          l0_manpltDirctnV = l0_manpltDirctn(i)
          choicSprV        = choicSpr(i) 
          call manplt_apclbsalprop_forAngls(cellNmV,l0_manpltDirctnV,choiceSprV)
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       call switchto_NI_model_run_and_switchbackto_TN
       
       call get_the_strght_lgcl(lgcl_strght)
       if (lgcl_strght .eqv. .True.) exit
    enddo
    
  end subroutine apcl_bsal_l0_manplt_of_cells_make_systm_strght
  
  
  subroutine get_the_strght_lgcl(lgcl_strght)
    implicit none
    logical, intent(out) :: lgcl_strght
    integer :: cellNmV
    integer :: first_TN,secnd_TN,third_TN,furth_TN
    real*8  :: tol_len,curr_dfrmtn
    
    lgcl_strght = .False. !strgtness is measured based on cell1
    
    cellNmV   = 1
    first_TN  = 2*(cellNmV-1)+1
    secnd_TN  = 2*(cellNmV-1)+2
    third_TN  = 2*(cellNmV-1)+3
    furth_TN  = 2*(cellNmV-1)+4
    
    tol_len     = 0.05d0*(node_xy(third_TN,1)-node_xy(first_TN,1))
    curr_dfrmtn = abs(node_xy(secnd_TN,1)-node_xy(first_TN,1))
    
    if (curr_dfrmtn .le. tol_len) lgcl_strght=.True.
    
  end subroutine get_the_strght_lgcl
  
  
  subroutine change_l0ksSpr_gvnCellNmAndEquill(cellNm,ApBsLt,TnsnDirc,maxCnt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm,ApBsLt,TnsnDirc,ExpNo
    integer, intent(inout) :: FrmNo 
    integer, intent(in)    :: maxCnt
    
    integer :: i,imax,countLp
    integer :: nsprInAp,nsprInBs,nsprInLt
    integer :: nsprsInACell,sprNmL,sprNmR
    real*8  :: hmuchL0orKS
    
    integer, allocatable :: apclSpr(:),bsalSpr(:),ltrlSpr(:)
    integer, allocatable :: sprsInTheSide(:)
    
    if (modelID==1) then
       nsprInAp=1 ; nsprInBs=1 ; nsprInLt=1 
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
       
    elseif (modelID==2) then
       nsprInAp=(NAEC_Apcl+1) ; nsprInBs=(NAEC_Bsal+1) ; nsprInLt=(NAEC_Ltrl+1)
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
    endif
    
    allocate(apclSpr(1:nsprInAp),bsalSpr(1:nsprInBs),ltrlSpr(1:nsprInLt))
    apclSpr=0 ; bsalSpr=0 ; ltrlSpr=0
    
    if (ApBsLt==1) imax=nsprInAp
    if (ApBsLt==2) imax=nsprInBs
    if (ApBsLt==3) imax=nsprInLt
    
    allocate(sprsInTheSide(1:imax))
    
    do i = 1,imax
       if (ApBsLt==1) apclSpr(i) = (cellNm-1)*(nsprsInACell)+i
       if (ApBsLt==2) bsalSpr(i) = (cellNm-1)*(nsprsInACell)+(nsprInAp)+i
       if (ApBsLt==3) ltrlSpr(i) = (cellNm-1)*(nsprsInACell)+(nsprInAp)+(nsprInBs)+i
    enddo
    
    if (ApBsLt==1) sprsInTheSide = apclSpr
    if (ApBsLt==2) sprsInTheSide = bsalSpr
    if (ApBsLt==3) sprsInTheSide = ltrlSpr
    
    if (TnsnDirc==+1) hmuchL0orKS = +0.05d0
    if (TnsnDirc==-1) hmuchL0orKS = -0.05d0
    
    countLp = 0
    
    do
       
       do i = 1,imax
          sprNmL        = sprsInTheSide(i)       ; sprNmR = sprNmL+(Hlf_Ncell)*(nsprsInACell)
          l0(sprNmL)    = (1.0d0-hmuchL0orKS)*l0(sprNmL)    ; l0(sprNmR)    = l0(sprNmL)
          k_spr(sprNmL) = (1.0d0+hmuchL0orKS)*k_spr(sprNmL) ; k_spr(sprNmR) = k_spr(sprNmL)
          write(*,*) sprNmL,sprNmR,i,"sprNmL-R"
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       call switchto_NI_model_run_and_switchbackto_TN
       countLp = countLp+1

       if (countLp==maxCnt) then
          write(*,*) "FrmNo bfr exiting in change_l0ksSpr routine",(FrmNo-1)
          exit
       endif
       
    enddo
    
  end subroutine change_l0ksSpr_gvnCellNmAndEquill
  
  !subroutine   
  
  
  subroutine change_l0_or_ks_gvnCellNmAndEquill(cellNm,ApBsLt,l0ORks,TnsnDirc,hmuch,&
       maxCnt,ExpNo,FrmNo)
    
    implicit none
    integer, intent(in)    :: cellNm,ApBsLt,l0ORks,TnsnDirc
    real*8 , intent(in)    :: hmuch
    integer, intent(in)    :: maxCnt,ExpNo
    integer, intent(inout) :: FrmNo 
    
    integer :: i,imax,countLp
    integer :: nsprInAp,nsprInBs,nsprInLt
    integer :: nsprsInACell,sprNmL,sprNmR
    real*8  :: hmuchL0orKS
    
    integer, allocatable :: apclSpr(:),bsalSpr(:),ltrlSpr(:)
    integer, allocatable :: sprsInTheSide(:)
    
    if (modelID==1) then
       nsprInAp=1 ; nsprInBs=1 ; nsprInLt=1 
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
       
    elseif (modelID==2) then
       nsprInAp=(NAEC_Apcl+1) ; nsprInBs=(NAEC_Bsal+1) ; nsprInLt=(NAEC_Ltrl+1)
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
    endif
    
    allocate(apclSpr(1:nsprInAp),bsalSpr(1:nsprInBs),ltrlSpr(1:nsprInLt))
    apclSpr=0 ; bsalSpr=0 ; ltrlSpr=0
    
    if (ApBsLt==1) imax=nsprInAp
    if (ApBsLt==2) imax=nsprInBs
    if (ApBsLt==3) imax=nsprInLt
    
    allocate(sprsInTheSide(1:imax))
    
    do i = 1,imax
       if (ApBsLt==1) apclSpr(i) = (cellNm-1)*(nsprsInACell)+i
       if (ApBsLt==2) bsalSpr(i) = (cellNm-1)*(nsprsInACell)+(nsprInAp)+i
       if (ApBsLt==3) ltrlSpr(i) = (cellNm-1)*(nsprsInACell)+(nsprInAp)+(nsprInBs)+i
    enddo
    
    if (ApBsLt==1) sprsInTheSide = apclSpr
    if (ApBsLt==2) sprsInTheSide = bsalSpr
    if (ApBsLt==3) sprsInTheSide = ltrlSpr
    
    if (TnsnDirc==+1) hmuchL0orKS = +abs(hmuch)
    if (TnsnDirc==-1) hmuchL0orKS = -abs(hmuch)
    
    countLp = 0
    
    do
       
       do i = 1,imax
          
          sprNmL        = sprsInTheSide(i)
          sprNmR = sprNmL+(Hlf_Ncell)*(nsprsInACell)
          
          if (l0ORks==1) then
             l0(sprNmL)    = (1.0d0-hmuchL0orKS)*l0(sprNmL)
             l0(sprNmR)    = l0(sprNmL)
             
          elseif (l0ORks==2) then
             k_spr(sprNmL) = (1.0d0+hmuchL0orKS)*k_spr(sprNmL)
             k_spr(sprNmR) = k_spr(sprNmL)
             
          elseif (l0ORks==3) then
             l0(sprNmL)    = (1.0d0-hmuchL0orKS)*l0(sprNmL)
             l0(sprNmR)    = l0(sprNmL)
             k_spr(sprNmL) = (1.0d0+hmuchL0orKS)*k_spr(sprNmL)
             k_spr(sprNmR) = k_spr(sprNmL)
          endif
          
          write(*,*) sprNmL,sprNmR,i,"sprNmL-R"
          
       enddo
       
       call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
       countLp = countLp+1
       
       if (countLp==maxCnt) then
          write(*,*) "FrmNo bfr exiting in change_l0ksSpr routine",(FrmNo-1)
          exit
       endif
       
    enddo
    
  end subroutine change_l0_or_ks_gvnCellNmAndEquill
  
  subroutine initiatr_cell_bsal_l0ksIncr(TnsnDirc,maxCnt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: TnsnDirc,maxCnt,ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,imax,countLp
    integer :: nsprInAp,nsprInBs,nsprInLt
    integer :: nsprsInACell,sprNm
    real*8  :: hmuchL0orKS
    
    
    if (modelID==1) then
       nsprInAp=1 ; nsprInBs=1 ; nsprInLt=1 
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
       
    elseif (modelID==2) then
       nsprInAp=(NAEC_Apcl+1) ; nsprInBs=(NAEC_Bsal+1) ; nsprInLt=(NAEC_Ltrl+1)
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
    endif
    
    if (TnsnDirc==+1) hmuchL0orKS = +0.05d0
    if (TnsnDirc==-1) hmuchL0orKS = -0.05d0
    
    countLp = 0
    imax    = nsprInBs
    
    do
       
       do i = 1,imax
          sprNm        = (N_cell-1)*(nsprsInACell) + i
          l0(sprNm)    = (1.0d0-hmuchL0orKS)*l0(sprNm)
          k_spr(sprNm) = (1.0d0+(5.0*hmuchL0orKS))*k_spr(sprNm)
          write(*,*) sprNm,i,"sprNm"
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       call switchto_NI_model_run_and_switchbackto_TN
       countLp = countLp+1
       
       if (countLp==maxCnt) then
          write(*,*) "FrmNo bfr exiting in change_l0ksSpr routine",(FrmNo-1)
          exit
       endif
       
    enddo
    
  end subroutine initiatr_cell_bsal_l0ksIncr
  
  
  
  subroutine initiatr_cell_TwoLtrl_l0ksIncr(TnsnDirc,maxCnt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: TnsnDirc,maxCnt,ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,imax,countLp
    integer :: nsprInAp,nsprInBs,nsprInLt
    integer :: nsprsInACell,sprNmL,sprNmR
    real*8  :: hmuchL0orKS
    
    
    if (modelID==1) then
       nsprInAp=1 ; nsprInBs=1 ; nsprInLt=1 
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
       
    elseif (modelID==2) then
       nsprInAp=(NAEC_Apcl+1) ; nsprInBs=(NAEC_Bsal+1) ; nsprInLt=(NAEC_Ltrl+1)
       nsprsInACell = nsprInAp+nsprInBs+nsprInLt
    endif
    
    if (TnsnDirc==+1) hmuchL0orKS = +0.05d0
    if (TnsnDirc==-1) hmuchL0orKS = -0.05d0
    
    countLp = 0
    imax    = nsprInBs
    
    do
       
       do i = 1,imax
          sprNmL   = (Hlf_Ncell*nsprsInACell)-(nsprInLt)+i; sprNmR = sprNmL+(Hlf_Ncell*nsprsInACell)
          l0(sprNmL)    = (1.0d0-hmuchL0orKS)*l0(sprNmL)    ; l0(sprNmR)    = l0(sprNmL)
          k_spr(sprNmL) = (1.0d0+hmuchL0orKS)*k_spr(sprNmL) ; k_spr(sprNmR) = k_spr(sprNmL)
          write(*,*) sprNmL,sprNmR,i,"sprNm"
       enddo
       
       call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
       countLp = countLp+1
       
       if (countLp==maxCnt) then
          write(*,*) "FrmNo bfr exiting in change_l0ksSpr routine",(FrmNo-1)
          exit
       endif
       
    enddo
    
  end subroutine initiatr_cell_TwoLtrl_l0ksIncr
  
  subroutine making_InitiatnPhase_IncmingCellsIdntcl(upToWhchCell,idealCell,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: upToWhchCell
    integer, intent(in)    :: idealCell
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: nsprsInACell
    integer :: apclSpr,bsalSpr,ltrlSpr
    
    real*8  :: ks_apcl,ks_bsal,ks_ltrl
    real*8  :: l0_apcl,l0_bsal,l0_ltrl
    real*8  :: ka_apcl,ka_bsal,ka_ltrl
    real*8  :: A0_apcl,A0_bsal,A0_ltrl
    real*8  :: ka_ideal,A0_ideal
    
    integer :: cellNmL,cellNmR,sprNmL,sprNmR
    integer :: i,j,jmax
    
    if (modelID==1) continue
    if (modelID==2) then
       write(*,*) "stop in making_InitiatnPhase_IncmingCells, not for modelID==2"
       stop
    endif
    
    nsprsInACell = 3
    
    ka_ideal = k_area(idealCell)
    A0_ideal = A0(idealCell)
    
    apclSpr = (idealCell-1)*(nsprsInACell)+1 ; bsalSpr = apclSpr+1 ; ltrlSpr = apclSpr+2
    
    ks_apcl = k_spr(apclSpr)  ; ks_bsal = k_spr(bsalSpr)  ; ks_ltrl = k_spr(ltrlSpr)
    l0_apcl = l0(apclSpr)     ; l0_bsal = l0(bsalSpr)     ; l0_ltrl = l0(ltrlSpr)
    ka_apcl = k_area(apclSpr) ; ka_bsal = k_area(bsalSpr) ; ka_ltrl = k_area(ltrlSpr)
    A0_apcl = A0(apclSpr)     ; A0_bsal = A0(bsalSpr)     ; A0_ltrl = A0(ltrlSpr)
    
    do i = 1,upToWhchCell
       
       cellNmL = i ; cellNmR = i+Hlf_Ncell
       
       k_area(cellNmL) = ka_ideal ; k_area(cellNmR) = k_area(cellNmL)
       A0(cellNmL)     = A0_ideal ; A0(cellNmR)     = A0(cellNmL)
       
       write(*,*) k_area(i),A0(i),i,"ka-A0-cell"
       write(*,*) " "
       
       do j = 1,nsprsInACell
          sprNmL = (cellNmL-1)*(nsprsInACell) + j
          sprNmR = sprNmL + (Hlf_Ncell*nsprsInACell)
          
          if (j==1) then
             k_spr(sprNmL) = ks_apcl ; l0(sprNmL) = l0_apcl
             k_spr(sprNmR) = ks_apcl ; l0(sprNmR) = l0_apcl
             
          elseif (j==2) then
             k_spr(sprNmL) = ks_bsal ; l0(sprNmL) = l0_bsal
             k_spr(sprNmR) = ks_bsal ; l0(sprNmR) = l0_bsal
             
          elseif (j==3) then
             k_spr(sprNmL) = ks_ltrl ; l0(sprNmL) = l0_ltrl
             k_spr(sprNmR) = ks_ltrl ; l0(sprNmR) = l0_ltrl
          endif
          
          write(*,*) k_spr(sprNmL),l0(sprNmL),sprNmL,"ks-l0-sprNmL"
          write(*,*) k_spr(sprNmR),l0(sprNmR),sprNmR,"ks-l0-sprNmR"
       enddo

       write(*,*) " "
    enddo
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
    
  end subroutine making_InitiatnPhase_IncmingCellsIdntcl

  
  subroutine copy_prp_frm_sprATosprB(cell_A,ApBsLt_A,cell_B,ApBsLt_B,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cell_A,ApBsLt_A
    integer, intent(in)    :: cell_B,ApBsLt_B
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,j
    integer :: spr_A,spr_B
    integer :: spr_AL,spr_BL,spr_AR,spr_BR
    integer :: nsprsInACell

    if (modelID==1) continue
    if (modelID==2) then
       write(*,*) "stopped as modelID==2"
       stop
    endif
    
    nsprsInACell = 3
    
    if (ApBsLt_A==1) spr_A = (cell_A-1)*(nsprsInACell) + 1
    if (ApBsLt_A==2) spr_A = (cell_A-1)*(nsprsInACell) + 2
    if (ApBsLt_A==3) spr_A = (cell_A-1)*(nsprsInACell) + 3
    
    if (ApBsLt_B==1) spr_B = (cell_B-1)*(nsprsInACell) + 1
    if (ApBsLt_B==2) spr_B = (cell_B-1)*(nsprsInACell) + 2
    if (ApBsLt_B==3) spr_B = (cell_B-1)*(nsprsInACell) + 3

    spr_BL = spr_B ; spr_BR = spr_BL+(Hlf_Ncell*nsprsInACell)
    spr_AL = spr_A ; spr_AR = spr_AL+(Hlf_Ncell*nsprsInACell)
    
    k_spr(spr_BL) = k_spr(spr_AL) ; k_spr(spr_BR) = k_spr(spr_AR)
    l0(spr_BL)    = l0(spr_AL)    ; l0(spr_BR)    = l0(spr_AR)
    
    write(*,*) spr_AL,spr_BL,"spr_AL-spr_BL"
    write(*,*) spr_AR,spr_BR,"spr_AR-spr_BR"
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
  end subroutine copy_prp_frm_sprATosprB
  
  
  subroutine savePropFrmPrvManplt(writORread,manpltNum)
    implicit none
    integer, intent(in) :: writORread
    integer, intent(in) :: manpltNum

    character(len=100) :: flnm,flnmbr,full_flnm
    
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 3

    flnm='InitiatnManplt'
    
    if ((manpltNum.gt.0) .and. (manpltNum.le.9)) then
       write(flnmbr,'(i1.1,a)') manpltNum,'.dat'
    elseif ((manpltNum.gt.9) .and. (manpltNum.le.99)) then
       write(flnmbr,'(i2.2,a)') manpltNum,'.dat'
    endif
    
    write(*,*) trim(adjustl(flnm))//trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=82,file=trim(adjustl(full_flnm)))
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (writORread==1) write(82,*) k_spr(j),l0(j),j
             if (writORread==2) read(82,*)  k_spr(j),l0(j)
             
          elseif (i==2) then 
             
             if (writORread==1) write(82,*) k_area(j),A0(j),j
             if (writORread==2) read(82,*)  k_area(j),A0(j)
             
          elseif (i==3) then
             
             if (writORread==1) write(82,*) CgXNode(j),j
             if (writORread==2) read(82,*)  CgXNode(j)
             
          endif
          
       enddo

       if (writORread==1) write(82,*) " "
    enddo
    
    close(82)
    
  end subroutine savePropFrmPrvManplt
  
  
  subroutine propsOfAllCompnts_writeORread(cellNm,cntInRtn,deallctnNeeds)
    implicit none
    integer, intent(in) :: cellNm
    integer, intent(in) :: cntInRtn
    integer, intent(in) :: deallctnNeeds
    !ka_tobeCopied,A0_tobeCopied,ks_tobeCopied(:),l0_tobeCopied(:)
    integer :: i,sprInArea,sprNm
    
    if (deallctnNeeds==1) then
       deallocate(ks_tobeCopied)
       deallocate(l0_tobeCopied)
    endif
    
    if (cntInRtn==1) then
       
       open(unit=232,file='propsOfAllCompntsSave.dat')
       write(232,*) ka_tobeCopied,A0_tobeCopied
       
       ka_tobeCopied = k_area(cellNm)
       A0_tobeCopied = A0(cellNm)
       
       
       sprInArea = area_spr(cellNm,0)
       allocate(ks_tobeCopied(1:sprInArea))
       allocate(l0_tobeCopied(1:sprInArea))
       
       do i = 1,sprInArea
          sprNm            = area_spr(cellNm,i)
          ks_tobeCopied(i) = k_spr(sprNm)
          l0_tobeCopied(i) = l0(sprNm)
          write(232,*) ks_tobeCopied(i),l0_tobeCopied(i)
       enddo
       
        
       close(232)
       
    elseif (cntInRtn==2) then

       open(unit=233,file='propsOfAllCompntsSave.dat')
       read(233,*) ka_tobeCopied,A0_tobeCopied
       
       k_area(cellNm) = ka_toBeCopied
       A0(cellNm)     = A0_tobeCopied
       
       sprInArea = area_spr(cellNm,0)
       allocate(ks_tobeCopied(1:sprInArea))
       allocate(l0_tobeCopied(1:sprInArea))
       
       do i = 1,sprInArea
          
          read(233,*) ks_tobeCopied(i),l0_tobeCopied(i)
          
          sprNm        = area_spr(cellNm,i)
          k_spr(sprNm) = ks_tobeCopied(i)
          l0(sprNm)    = l0_tobeCopied(i)
       enddo
       
       close(233)
       
    endif
    
  end subroutine propsOfAllCompnts_writeORread
  
  
  subroutine adjust_the_firstFrm_of_EpithelialLayer(ExpNo,FrmNo)
    implicit none        !we will take the strctpropsAddedCell_TN and copy the ideal cell properties
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: N_itm
    integer :: i,j,jmax
    real*8  :: ksV(1:N_spr),kaV(1:N_cell)
    real*8  :: l0V(1:N_spr),A0V(1:N_cell)
    real*8  :: CgXV(1:N_node)
    integer :: idealCell
    integer :: sprNmIdeal,nsprsInTN
    integer :: cellNmL,cellNmR
    integer :: sprNmL,sprNmR,sprNm
    
    real*8               :: kaIdeal,A0Ideal
    real*8, allocatable  :: ksIdeal(:),l0Ideal(:)
    
    open(unit=78,file='strctPropsAddedCellTN2S1T1.dat')
    
    N_itm = 3
    
    do i = 1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          if (i==1) read(78,*) ksV(j),l0V(j)
          if (i==2) read(78,*) kaV(j),A0V(j)
          if (i==3) read(78,*) CgXV(j)
       enddo
          
    enddo
    
    close(78)
    
    idealCell = 1
    
    kaIdeal = kaV(idealCell)
    A0Ideal = A0V(idealCell)
    
    nsprsInTN = 3
    allocate(ksIdeal(1:nsprsInTN))
    allocate(l0Ideal(1:nsprsInTN))
    
    do i = 1,nsprsInTN
       sprNmIdeal = (idealCell-1)*(nsprsInTN) + i 
       ksIdeal(i) = ksV(sprNmIdeal)
       l0Ideal(i) = l0V(sprNmIdeal)
    enddo
    
    
    do i = 1,Hlf_Ncell
       
       cellNmL=i                 ; cellNmR=i+Hlf_Ncell
       k_area(cellNmL) = kaIdeal ; A0(cellNmL) = A0Ideal
       k_area(cellNmR) = kaIdeal ; A0(cellNmR) = A0Ideal
       
       do j = 1,nsprsInTN
          sprNmL=(cellNmL-1)*(nsprsInTN)+j ; sprNmR=sprNmL+(Hlf_Ncell*nsprsInTN)
          k_spr(sprNmL) = ksIdeal(j) ; l0(sprNmL) = l0Ideal(j)
          k_spr(sprNmR) = ksIdeal(j) ; l0(sprNmR) = l0Ideal(j)
       enddo
          
    enddo
    
    k_area(N_cell) = kaIdeal ; A0(N_cell) = A0Ideal

    do i = 1,2
       sprNm        = (N_cell-1)*(nsprsInTN)+i
       k_spr(sprNm) = ksIdeal(i)
       l0(sprNm)    = l0Ideal(i)
    enddo
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
  end subroutine adjust_the_firstFrm_of_EpithelialLayer
  
  
  
  subroutine incrAllA0val_andReduceAll_l0Val(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: Aval(1:N_cell)
    integer :: i,j
    
    Aval(1:N_cell) = A(1:N_cell)
    
    do i = 1,N_cell
       A0(i) = 1.20d0*A0(i)
    enddo
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)

    do j = 1,5
       do i = 1,N_spr
          l0(i) = 0.95d0*l0(i)
       enddo
       call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    enddo
       
  end subroutine incrAllA0val_andReduceAll_l0Val
  
  
  subroutine apply_ThreeDiff_force_at_TopNode_of_LS_ofAcell_andEquill(cellNm,dirctn,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm
    integer, intent(in)    :: dirctn ! +1 for upward force, -1 for downward force
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: TopNodeOfLS1,TopNodeOfLS2 ! LS=LateralSide
    real*8  :: CgXSign
    
    TopNodeOfLS1 = (cellNm*2) + 1
    TopNodeOfLS2 = (Hlf_Ncell+1)*2 + (cellNm*2) + 1
    
    write(*,*) TopNodeOfLS1,TopNodeOfLS2,"TopNode L,R"
    
    if (dirctn==+1) CgXSign=-1.0d0
    if (dirctn==-1) CgXSign=+1.0d0
    
    write(*,*) ExpNo,FrmNo,"FrmNo in apply Force at Begin"
    
    CgXNode(TopNodeOfLS1) = (CgXSign*lowrF) ; CgXNode(TopNodeOfLS2) = (CgXSign*lowrF)
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    call check_grd_at_both_TNandNI
    
    CgXNode(TopNodeOfLS1) = (CgXSign*medmF) ; CgXNode(TopNodeOfLS2) = (CgXSign*medmF)
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    call check_grd_at_both_TNandNI
    
    CgXNode(TopNodeOfLS1) = (CgXSign*higrF) ; CgXNode(TopNodeOfLS2) = (CgXSign*higrF)
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    call check_grd_at_both_TNandNI
    
    CgXNode(TopNodeOfLS1) = zeroF ; CgXNode(TopNodeOfLS2) = zeroF
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
    write(*,*) CgXNode(1),CgXNode(TopNodeOfLS1),"CGNODE"
    write(*,*) ExpNo,(FrmNo-1),"FrmNo in apply force at End"
    
    
  end subroutine apply_ThreeDiff_force_at_TopNode_of_LS_ofAcell_andEquill
  
  
  subroutine apply_ThreeDiff_force_at_IntrmdNode(cellNm,ApBsLt,sgmntNo,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm,ApBsLt,sgmntNo,ExpNo
    integer, intent(inout) :: FrmNo
    
    !call switch_to_NI
    
    ! if (ApBsLt==1) then
       
    !    if (NAEC_Apcl.ne.0) then
    !       continue
    !    else
    !       write(*,*) "Not in Apical side"
    !    endif
          
    ! elseif (ApBsLt==2) then
       
    !    if (NAEC_Bsal.ne.0) then
          
    !    else
    !       continue
    !    endif
       
    ! elseif (ApBsLt==3) then
       
    !    if (NAEC_Ltrl.ne.0) then
          
    !    else
          
    !    endif
       
    ! endif
    
  end subroutine apply_ThreeDiff_force_at_IntrmdNode
  
  ! subroutine chk_forces_at_joint_node()
  !   implicit none
  !   real*8, allocatable :: forceAtTN(:,:),forceAtNI(:,:)
  !   real*8, allocatable :: CgXSave(:)
    
  !   integer :: i,j
    
  !   allocate(forceAtTN(1:N_node,1:N_dmnsn))
  !   allocate(CgXSave(1:N_node))
    
  !   open(unit=651,file='chk_forces_at_joint_node.dat')
    
  !   CgXSave(1:N_node) = CgXNode(1:N_node)

  !   write(651,*) select_xy,"slctXT" 
  !   write(651,*) "writing before changes"
  !   do i = 1,N_node
  !      write(651,*) CgXNode(i),i,"CgXVal"
  !   enddo
    
  !   close(651)
    
  !   call get_
  !   forceAtTN()
    
  ! end subroutine chk_forces_at_joint_node
  
  
  subroutine very_stiff_areas(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8, allocatable :: kaSave(:),A0Save(:)
    real*8  :: incrmntFac
    integer :: i,j
    
    allocate(kaSave(1:N_cell),A0Save(1:N_cell))
    
    kaSave(1:N_cell) = k_area(1:N_cell)
    A0Save(1:N_cell) = A0(1:N_cell)
    
    !now I will equate the A0 as A, and increase the k_area by say 10 times
    
    A0(1:N_cell) = A(1:N_cell)
    
    incrmntFac = 10.0d0
    
    do i = 1,N_cell
       k_area(i) = (incrmntFac)*(k_area(i))
    enddo
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
  end subroutine very_stiff_areas
  
  ! routines for applyling force before the lateral membranes actually meet
  
  subroutine read_get_the_conditions_force_bfr_apclSpr_demlish(demlishingSpr,initialL,tolrncL)
    implicit none
    integer, intent(out) :: demlishingSpr
    real*8 , intent(out) :: initialL,tolrncL
    
    integer  :: nsprsInACellTN,nsprsInACellNI 
    real*8   :: prcntTol
    
    nsprsInACellTN  = 1+1+1 
    nsprsInACellNI  = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1) 
    
    write(*,*) nsprsInACellTN, nsprsInACellNI,"nsrpsInACell TN,NI" 
    
    demlishingSpr  = (N_cell-1)*(nsprsInACellTN) + 1
    initialL       = l(demlishingSpr)
    prcntTol       = 0.85d0
    tolrncL        = (prcntTol)*(initialL)
    write(*,*) initialL,tolrncL,"initialL-tol L"
    
  end subroutine read_get_the_conditions_force_bfr_apclSpr_demlish
  
  
  subroutine chk_condn_bfrPutting_downwardForce(demlishingSpr,tolrncL,lgcl_condn)
    implicit none
    integer, intent(in)  :: demlishingSpr
    real*8 , intent(in)  :: tolrncL
    logical, intent(out) :: lgcl_condn
    
    real*8 :: currLen
    
    lgcl_condn = .False.
    currLen = l(demlishingSpr)
    
    if (currLen.lt.tolrncL) then
       lgcl_condn = .True.
    endif
    
  end subroutine chk_condn_bfrPutting_downwardForce
  
  subroutine chkCondtnAndPutForceBfrPulleyPntCreatn(demlishingSpr,tolrncL,ExpNo,FrmNo) 
    implicit none
    integer, intent(in)    :: demlishingSpr,ExpNo
    real*8 , intent(in)    :: tolrncL
    integer, intent(inout) :: FrmNo
    
    real*8   :: currLen 
    integer  :: nsprsInACellTN,nsprsInACellNI 
    integer  :: nodeL=-1,nodeR=-2
    
    call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
    
    currLen = l(demlishingSpr)
    write(*,*) currLen,tolrncL,"curr-tol L"
    
    if (currLen.lt.tolrncL) then 
       call change_NodeTypForThetobejoiningNodes(nodeL,nodeR)
       call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
       call applyForceAtToBeJoiningNodes(nodeL,nodeR,ExpNo,FrmNo) 
    endif
    
  end subroutine chkCondtnAndPutForceBfrPulleyPntCreatn
  
  
  
  subroutine applyForceAtToBeJoiningNodes(nodeL,nodeR,ExpNo,FrmNo) 
    implicit none 
    integer, intent(in)    :: nodeL,nodeR 
    integer, intent(in)    :: ExpNo 
    integer, intent(inout) :: FrmNo 
    
    integer  :: i, N_diffForces 
    real*8   :: incrCgX 
    
    
    N_diffForces = 3 
    
    do i = 1,N_diffForces 
       
       incrCgX = real(i-1)*0.05d0 
       
       CgYNode(nodeL) = 0.15d0 + incrCgX 
       CgYNode(nodeR) = CgYNode(nodeL) 
       
       call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
       
    enddo
    write(*,*) "stopping in applyForceAtToBeJoiningNodes"
    !stop
    
  end subroutine applyForceAtToBeJoiningNodes
  
  subroutine apply_Force_aft_ManuallyShifting_Nodes(ExpNo,FrmNo) 
    implicit none 
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer  :: nodeL,nodeR 
    integer  :: i,i1,j,jmax,N_diffForces 
    real*8   :: incrCgY
    integer  :: fileNo
    real*8   :: E
    integer  :: cell_val,crtclNode
    
    call get_the_shifting_nodes(nodeL,nodeR)
    N_diffForces = 21
    
    call find_the_cell_with_diagonal_spr(cell_val)
    crtclNode = (2*cell_val)+3
    write(*,*) cell_val,crtclNode,"cell_val and crtclNode"
    
    do i = 1,N_diffForces
       
       incrCgY        = real(i-1)*0.05d0   
       CgYNode(nodeL) = incrCgY 
       CgYNode(nodeR) = CgYNode(nodeL) 
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside apply_Force_aft_Man"
       
       ! do i1 = 1,4
          
       !    if (i1==1) jmax=N_spr
       !    if (i1==2) jmax=N_cell
       !    if (i1==3) jmax=N_node
       !    if (i1==4) jmax=N_node
       !   
       !    do j = 1,jmax
       !       if (i1==1) write(*,*) k_spr(j), l0(j),l(j),j,"sprprop  bfr calling"
       !       if (i1==2) write(*,*) k_area(j),A0(j),A(j),j,"cellprop bfr calling"
       !       if (i1==3) write(*,*) CgXNode(j),CgYNode(j),k_phi(j,1:4),j,"CgX/Y,kphi prop bfr calling"
       !       if (i1==4) write(*,*) node_xy(j,1),node_xy(j,2),j,"node_xy(j,1:2)"
       !    enddo
          
       !    write(*,*) " "
       ! enddo
       
       call Equilibrate_only_NI_model
       
       !call Equilibrate_bothTN_NImodel(ExpNo,FrmNo)
       !write(*,*) (FrmNo-1),i,"FrmNos in apply_force"
       
       !fileNo=i+1 ; call propCheckForDiffFile(fileNo)
       !write(*,*) A0(N_cell),"A0-N_cell"
       
       ! if (i==1) then
       !    call switch_to_NI_model
       !    E = Energy(node_xy,l0,A0)
       !    allocate(ksv2NI(1:N_spr),kav2NI(1:N_cell),l0v2NI(1:N_spr),A0v2NI(1:N_cell))
       !    allocate(CgXv2NI(1:N_node),CgYv2NI(1:N_node))
       !    ksv2NI=k_spr ; kav2NI=k_area ; l0v2NI=l0 ; A0v2NI=A0 ; CgXv2NI=CgXNode ; CgYv2NI=CgYNode
       !    call not_Equlibrate_backto_TN
       !    !exit
       ! endif
       
       restrctring_xval = 0.31d0
       
       if (abs(node_xy(crtclNode,1)) .le. restrctring_xval) then
          write(*,*) i,"loop val Val at exit"
          exit
       endif
       
    enddo
    
    write(*,*) "apply_Force_aft_ManuallyShifting_Nodes"
    !stop
    
  end subroutine apply_Force_aft_ManuallyShifting_Nodes
  
  subroutine apply_force_at_Apcl_surfc_ofIC_bfr_adding_VFregion(N_diffForces)
    implicit none
    integer, intent(in) :: N_diffForces
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax
    real*8              :: incrCgY
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. 11) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    do i = 1,N_diffForces
       
       incrCgY        = real(i-1)*(0.0200d0)   
       CgYNode(nodeL) = incrCgY 
       CgYNode(nodeR) = CgYNode(nodeL)
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside apply_Force"
       call Equilibrate_only_NI_model
       
    enddo
    
  end subroutine apply_force_at_Apcl_surfc_ofIC_bfr_adding_VFregion
  
  subroutine applyForceApclSurfcOfIC_WT_BndryAdjst_bfr_adding_VFregion(N_diffForces,IncrAmnt)
    implicit none
    integer, intent(in) :: N_diffForces
    real*8,  intent(in) :: IncrAmnt
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax
    real*8              :: incrCgY
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. 11) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    do i = 1,N_diffForces       
       incrCgY        = real(i)*IncrAmnt   
       CgYNode(nodeL) = incrCgY 
       CgYNode(nodeR) = CgYNode(nodeL)
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside apply_Force"
       call Equilibrate_only_NI_model ; write(*,*) (Frame_NI-1),"aft force appcln"
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*)Frame_NI-1,"aft pull push"
    enddo
    
  end subroutine applyForceApclSurfcOfIC_WT_BndryAdjst_bfr_adding_VFregion
  
  subroutine apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion(N_diffForces)
    implicit none
    integer, intent(in) :: N_diffForces
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax
    real*8              :: incrCgY
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. 11) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    do i = 1,N_diffForces
       
       incrCgY        = real(i-1)*(0.0500d0)   
       CgYNode(nodeL) = incrCgY 
       CgYNode(nodeR) = CgYNode(nodeL)
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside apply_Force VF"
       call Equilibrate_only_NI_model_withVF_region
       
    enddo
    
  end subroutine apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion
  
  subroutine apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion_frmCurrntVal(N_diffForces)
    implicit none
    integer, intent(in) :: N_diffForces
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax
    real*8              :: incrCgY,initCgYNode
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. 11) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    initCgYNode = CgYNode(nodeL)
    
    do i = 1,N_diffForces
       
       incrCgY        = real(i-1)*(0.0500d0)   
       CgYNode(nodeL) = initCgYnode + incrCgY 
       CgYNode(nodeR) = CgYNode(nodeL)
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside apply_Force VF"
       call Equilibrate_only_NI_model_withVF_region
       
    enddo
    
  end subroutine apply_force_at_Apcl_surfc_ofIC_aft_adding_VFregion_frmCurrntVal
  
  subroutine adjust_force_at_joint_ApclNode_ofIC_to_maintain_PP(CMv,ExpctdPPbfrVal,ExpctdPPaftVal)
    implicit none
    integer, intent(in) :: CMv
    real*8 , intent(in) :: ExpctdPPbfrVal, ExpctdPPaftVal
    real*8              :: CurrPPbfrVal,CurrPPaftVal
    integer             :: nodeL,nodeR,TN_bfr,TN_aft
    integer             :: i,j,jmax
    real*8              :: modfCgYperStp
    integer             :: movmntDirctn,movmntDirctnBegin,cntLp
    
    call get_the_shifting_nodes(nodeL,nodeR)
    call get_the_bfraft_TNof_PPnt(CMv,TN_bfr,TN_aft)
    
    cntLp        = 0
    CurrPPbfrVal = abs(node_xy(TN_bfr,1)) ; CurrPPaftVal = abs(node_xy(TN_aft,2))
    write(*,*) CurrPPbfrVal,CurrPPaftVal,cntLp,"CurrPPbfr~aft at cntLp=0"
    
    if (CurrPPbfrVal.gt.ExpctdPPbfrVal) then
       movmntDirctnBegin  = +1
       modfCgYperStp      = +0.0200d0
       
    elseif (CurrPPbfrVal.le.ExpctdPPbfrVal) then
       movmntDirctnBegin  = -1
       modfCgYperStp      = -0.0200d0
    endif
    
    do   
       CgYNode(nodeL) = CgYNode(nodeL) + modfCgYperStp  
       CgYNode(nodeR) = CgYNode(nodeL)
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,"CgY vals inside maintain PP"
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
       CurrPPbfrVal = abs(node_xy(TN_bfr,1)) ; CurrPPaftVal = abs(node_xy(TN_aft,2))
       write(*,*) CurrPPbfrVal,CurrPPaftVal,cntLp,"CurrPPbfr~aft,cntLp"
       
       if (CurrPPbfrVal.gt.ExpctdPPbfrVal) movmntDirctn = +1
       if (CurrPPbfrVal.le.ExpctdPPbfrVal) movmntDirctn = -1
       
       if (movmntDirctn.ne.movmntDirctnBegin) then
          CgYNode(nodeL) = CgYNode(nodeL) - 0.5000d0*modfCgYperStp  
          CgYNode(nodeR) = CgYNode(nodeL)
          call Equilibrate_only_NI_model
          write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,cntLp,"at Exit"
          exit
       endif
       
    enddo
    
  end subroutine adjust_force_at_joint_ApclNode_ofIC_to_maintain_PP
  
  subroutine pull_push_atBndry(pullpushIndctr,CgXChng)
    implicit none
    integer, intent(in) :: pullpushIndctr
    real*8,  intent(in) :: CgXChng
    
    integer             :: bndryLT,bndryLB
    integer             :: bndryRT,bndryRB
    
    bndryLT=1 ; bndryRT=bndryLT+(Hlf_Ncell+1)*2
    bndryLB=2 ; bndryRB=bndryLB+(Hlf_Ncell+1)*2 ; write(*,*) bndryLT,bndryRT,bndryLB,bndryRB,"LRTB"
    
    CgXNode(bndryLT) = CgXNode(bndryLT) + (pullpushIndctr)*(CgXChng)
    CgXNode(bndryLB) = CgXNode(bndryLB) + (pullpushIndctr)*(CgXChng)
    CgXNode(bndryRT) = CgXNode(bndryRT) - (pullpushIndctr)*(CgXChng)
    CgXNode(bndryRB) = CgXNode(bndryRB) - (pullpushIndctr)*(CgXChng)
    
    write(*,*)CgXNode(bndryLT),CgXNode(bndryLB),CgXNode(bndryRT),CgXNode(bndryRB),"CgXNode"
    
    call Equilibrate_only_NI_model
    
  end subroutine pull_push_atBndry
  
  
  subroutine remveForceFrmApclSrfcOfIC_frmCurrntVal()
    implicit none
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax,zeroIndctr
    real*8              :: incrCgY,initCgYNode
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. 11) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    initCgYNode = CgYNode(nodeL)

    zeroIndctr=0
    i = 1
    
    do   
       incrCgY        = real(i-1)*(0.0500d0)
       
       if (incrCgY.gt.initCgYNode) then
          CgYNode(nodeL) = 0.00d0
          CgYNode(nodeR) = 0.00d0
          zeroIndctr     = 1
       else
          CgYNode(nodeL) = initCgYnode - incrCgY 
          CgYNode(nodeR) = CgYNode(nodeL)
       endif
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,i,"CgY vals inside removal_Force"
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
       if (zeroIndctr==1) exit
       i=i+1
    enddo
    
  end subroutine remveForceFrmApclSrfcOfIC_frmCurrntVal
  
  
  subroutine remveOrAdd_ForceAtApclSrfcOfIC_ToMaintainSpcfcVal(TIIF_value)
    implicit none
    real*8, intent(in)  :: TIIF_value
    integer             :: nodeL,nodeR
    integer             :: i,j,jmax,extIndctr
    real*8              :: modfCgY,initCgYNode,initCgYDiff
    
    write(*,*) Hlf_Ncell,"Hlf_Ncell"
    if (Hlf_Ncell .ne. (N_cell-1)/2) stop 'wrong Hlf_Ncell'
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    initCgYNode = CgYNode(nodeL)
    initCgYDiff = abs(initCgYNode-TIIF_value)
    
    extIndctr=0
    i = 1
    
    do   
       modfCgY = real(i-1)*(0.0500d0)
       
       if (modfCgY.gt.initCgYDiff) then
          CgYNode(nodeL) = TIIF_value
          CgYNode(nodeR) = TIIF_value
          extIndctr     = 1
       else
          
          if (CgYNode(nodeL).gt.TIIF_value) then
             CgYNode(nodeL) = initCgYnode - modfCgY 
             CgYNode(nodeR) = CgYNode(nodeL)
          elseif (CgYNode(nodeL).lt.TIIF_value) then
             CgYNode(nodeL) = initCgYnode + modfCgY 
             CgYNode(nodeR) = CgYNode(nodeL)
          endif
          
       endif
       
       write(*,*) CgYNode(nodeL),CgYNode(nodeR),nodeL,nodeR,i,"CgY vals inside removal_Force"
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
       if (extIndctr==1) exit
       i=i+1
    enddo
    
  end subroutine remveOrAdd_ForceAtApclSrfcOfIC_ToMaintainSpcfcVal
  
  
  subroutine get_the_shifting_nodes(nodeL,nodeR)
    implicit none
    integer, intent(out) :: nodeL,nodeR
    
    nodeL = (Hlf_Ncell)*2   + 1
    nodeR = (Hlf_Ncell+1)*2 + nodeL
    
    write(*,*) nodeL,nodeR,"nodeL-R in the get the shifting nodes"
    
  end subroutine get_the_shifting_nodes
  
  subroutine get_the_bfraft_TNof_PPnt(CMv,TN_bfr,TN_aft)
    implicit none
    integer, intent(in)  :: CMv
    integer, intent(out) :: TN_bfr,TN_aft
    
    if (CMv==0) stop 'no PP in CMv=0'
    
    TN_aft = 2*(Hlf_Ncell+1)-1-(CMv-1)*2
    TN_bfr = 2*(Hlf_Ncell+1)-1-2-(CMv-1)*2
    
    write(*,*) CMv,TN_bfr,TN_aft,"CMv,TN_bfr,TN_aft"
    
  end subroutine get_the_bfraft_TNof_PPnt
  
  subroutine control_of_vertical_force_inICApcl_vs_VitellineFluid_strength
    implicit none
    integer :: nodeL,nodeR
    
    nodeL = (2*Hlf_Ncell) + 1       ; write(*,*) nodeL,Hlf_Ncell,"in control L"
    nodeR = 2*(Hlf_Ncell+1) + nodeL ; write(*,*) nodeR,"in control R"
    
  end subroutine control_of_vertical_force_inICApcl_vs_VitellineFluid_strength
  
  
  subroutine control_of_vertical_force_inICApcl_vs_boundary_force
    implicit none
    integer :: nodeL,nodeR
    real*8  :: forceCoeff,modfCgY
    real*8  :: quotientV,remaindrV
    integer :: numOfSteps
    integer :: i,j
    real*8  :: forceCoeffAtStp1
    
    if (modelID==1) stop 'Edit the Equilibrate portion for TN system exprmnt'
    
    nodeL = (2*Hlf_Ncell) + 1       ; write(*,*) nodeL,Hlf_Ncell,"in control L2"
    nodeR = 2*(Hlf_Ncell+1) + nodeL ; write(*,*) nodeR,"in control R2"
    
    forceCoeff = CgYNode(nodeL) ; write(*,*) forceCoeff,CgYNode(nodeL),CgYNode(nodeR),"forceCoeff"
    modfCgY    = 0.0500d0
    quotientV  = (forceCoeff/modfCgY) ; remaindrV = forceCoeff-(quotientV*modfCgY)
    write(*,*)    forceCoeff,quotientV,remaindrV,modfCgY,"fqrm"
    numOfSteps = int(quotientV)+1 ; write(*,*) numOfSteps,"nos"

    
    do i = 0,numOfSteps
       
       if (i==0) then
          continue
       elseif (i==1) then
          CgYNode(nodeL)   = CgYNode(nodeL) - remaindrV 
          CgYNode(nodeR)   = CgYNode(nodeL)
          forceCoeffAtStp1 = CgYNode(nodeL)
       else
          CgYNode(nodeL)   = forceCoeffAtStp1 - (i-1)*(modfCgY)
       endif
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
    enddo
    
  end subroutine control_of_vertical_force_inICApcl_vs_boundary_force
  
  subroutine change_NodeTypForThetobejoiningNodes(nodeL,nodeR)
    implicit none 
    integer, intent(out) :: nodeL,nodeR 
    
    nodeL=(Hlf_Ncell)*2+1 
    nodeR=(Hlf_Ncell+1)*2+nodeL
    
    write(*,*) nodeL,nodeR,"nodes" 
    write(*,*) node_typ(nodeL),node_typ(nodeR),"check node_typs" 
    
    if (node_typ(nodeL) .ne. 2) then
       write(*,*) "node_typ must be 2 here"
    endif
    
    node_typ(nodeL) = 1 
    node_typ(nodeR) = 1 ! make both the nodes free 
    
  end subroutine change_NodeTypForThetobejoiningNodes
  
  subroutine propCheckForDiffFile(fileNo)
    implicit none
    integer, intent(in) :: fileNo
    character(len=100)  :: fileNm

    integer :: i,j,imax,jmax
    
    write(fileNm,'(a,i2.2,a)')"propCheck",fileNo,'.dat'
    write(*,*)trim(adjustl(fileNm))
    
    open(unit=244,file=trim(adjustl(fileNm)))

    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node 
       
       do j = 1,jmax
          if (i==1) write(244,*) k_spr(j),l0(j),j
          if (i==2) write(244,*) k_area(j),A0(j),j
          if (i==3) write(244,*) CgXNode(j),CgYNode(j),j
       enddo
    
    enddo
    
    close(244)
    
  end subroutine propCheckForDiffFile
  
  
  subroutine Apcl_shorten_Adjacent_of_IC_and_incr_vrtcl_force_to_hold(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: hmL0Ap,hmL0Bs,hmL0Lt
    integer :: cell_val
    real*8  :: CgY_incr
    integer :: nodeL,nodeR
    real*8  :: y_nodeL,y_nodeR
    integer :: cntLp1,cntLp2
    
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    
    hmL0Ap   = 0.85d0 ; CgY_incr = 0.05d0
    cell_val = Hlf_Ncell
    
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)))
    call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp,BsalSp,LtrlSp)
    
    call get_the_shifting_nodes(nodeL,nodeR)
    y_nodeL = node_xy(nodeL,2) ; y_nodeR = node_xy(nodeR,2)
    write(*,*) nodeL,nodeR,y_nodeL,y_nodeR,"nodeL,nodeR,y_node L and R"
    
    write(*,*) CgYNode(nodeL),CgYNode(nodeR),"prv_CgYnodeL"
    !stop
    
    
    do cntLp1=1,6
       
       write(*,*) "cntLp1 value and Frame_NI",cntLp1,Frame_NI-1
       
       call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
       call Equilibrate_only_NI_model
       write(*,*) node_xy(nodeL,2),node_xy(nodeR,2),"ynodeL aft Apcl Incr"
       write(*,*) Frame_NI-1,"Frame_NI aft Apcl Incr"
       
       cntLp2 = 0
       
       do
          
          CgYNode(nodeL) = CgYNode(nodeL) + CgY_incr
          CgYNode(nodeR) = CgYNode(nodeL)
          call Equilibrate_only_NI_model
          write(*,*) node_xy(nodeL,2),node_xy(nodeR,2),"ynodeL aft Cg Incr"
          write(*,*) CgYNode(nodeL),CgYNode(nodeR),"CgYNode in Lp"
          
          write(*,*) Frame_NI-1,"Frame_NI aft Cg Incr"
          cntLp2 = cntLp2+1

          if (cntLp2.ge.5) then
             write(*,*) "Exiting for bounding"
             exit
          endif
          
          if (abs(node_xy(nodeL,2)) .gt. abs(y_nodeL)) then
             write(*,*) abs(node_xy(nodeL,2)),abs(y_nodeL),"abs val"
             exit
          endif
          
       enddo
       
    enddo
    
  end subroutine Apcl_shorten_Adjacent_of_IC_and_incr_vrtcl_force_to_hold
  
  
  subroutine Ltrl_Shorten_of_Cell_NIsys(cellNm)
    implicit none
    integer, intent(in)  :: cellNm
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    
    integer :: cntLp
    real*8  :: hmL0
    
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)))
    ApclSp = -1 ; BsalSp = -1 ; LtrlSp = -1
    
    call get_ApclBsalLtrl_OfNIsystm(cellNm,ApclSp,BsalSp,LtrlSp)
    write(*,*) ApclSp,"Apcl" ; write(*,*) BsalSp,"Bsal" ; write(*,*) LtrlSp,"Ltrl"
    
    hmL0 = 0.9500d0
    
    do cntLp = 1,5
       call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0)
       call Equilibrate_only_NI_model_withVF_region
    enddo
    
  end subroutine Ltrl_Shorten_of_Cell_NIsys
  
  
  subroutine Sprshorten_of_Cell_NIsys(cellNm,ApBsLt,hmL0,cntLpMax)
    implicit none
    integer, intent(in)  :: cellNm
    integer, intent(in)  :: ApBsLt
    real*8 , intent(in)  :: hmL0
    integer, intent(in)  :: cntLpMax
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    integer              :: cntLp
    
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)))
    ApclSp = -1 ; BsalSp = -1 ; LtrlSp = -1
    
    call get_ApclBsalLtrl_OfNIsystm(cellNm,ApclSp,BsalSp,LtrlSp)
    write(*,*) ApclSp,"Apcl" ; write(*,*) BsalSp,"Bsal" ; write(*,*) LtrlSp,"Ltrl"
    
    do cntLp = 1,cntLpMax
       
       if (ApBsLt==1) call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0)
       if (ApBsLt==2) call manplt_spr_NIsys(BsalSp,NAEC_Bsal,hmL0)
       if (ApBsLt==3) call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0)
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
    enddo
    
  end subroutine Sprshorten_of_Cell_NIsys
  
  
  subroutine tilt_reduction_of_the_boundry(forCellsMeet,PCSorICS)
    implicit none
    integer, intent(in) :: forCellsMeet,PCSorICS
    
    real*8  :: hmL0Ap0,hmL0Ap1,hmL0Ap2,hmL0Ap3,hmL0Ap4,hmL0Ap5,hmL0Ap6,hmL0Ap7,hmL0Ap8,hmL0Ap9
    real*8  :: hmL0Bs0,hmL0Bs1,hmL0Bs2,hmL0Bs3,hmL0Bs4,hmL0Bs5,hmL0Bs6,hmL0Bs7,hmL0Bs8,hmL0Bs9
    real*8  :: hmL0Lt1,hmL0Lt2,hmL0Lt3,hmL0Lt4,hmL0Lt5,hmL0Lt6,hmL0Lt7,hmL0Lt8,hmL0Lt9
    real*8  :: hmL0Ap10,hmL0Ap11,hmL0Bs10,hmL0Bs11,hmL0Lt10,hmL0Lt11
    
    integer, allocatable :: ApclSp1(:), BsalSp1(:), LtrlSp1(:)
    integer, allocatable :: ApclSp2(:), BsalSp2(:), LtrlSp2(:)
    integer, allocatable :: ApclSp3(:), BsalSp3(:), LtrlSp3(:)
    integer, allocatable :: ApclSp4(:), BsalSp4(:), LtrlSp4(:)
    integer, allocatable :: ApclSp5(:), BsalSp5(:), LtrlSp5(:)
    integer, allocatable :: ApclSp6(:), BsalSp6(:), LtrlSp6(:)
    integer, allocatable :: ApclSp7(:), BsalSp7(:), LtrlSp7(:)
    integer, allocatable :: ApclSp8(:), BsalSp8(:), LtrlSp8(:)
    integer, allocatable :: ApclSp9(:), BsalSp9(:), LtrlSp9(:)
    integer, allocatable :: ApclSp10(:),BsalSp10(:),LtrlSp10(:)
    integer, allocatable :: ApclSp11(:),BsalSp11(:),LtrlSp11(:)
    integer, allocatable :: ApclSp0(:), BsalSp0(:)
    
    integer :: readTheFile,FrmNoIncr
    integer :: cnt1,cnt2,cnt,cell_val
    integer :: caseVal
    
    integer :: i1,i1max,i2
    integer :: forceNodeL,forceNodeR
    integer :: BtmBndryNodeL,BtmBndryNodeR
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),  BsalSp1(1:(NAEC_Bsal+1)),  LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),  BsalSp2(1:(NAEC_Bsal+1)),  LtrlSp2(1:(NAEC_Ltrl+1)))
    allocate(ApclSp3(1:(NAEC_Apcl+1)),  BsalSp3(1:(NAEC_Bsal+1)),  LtrlSp3(1:(NAEC_Ltrl+1)))
    allocate(ApclSp4(1:(NAEC_Apcl+1)),  BsalSp4(1:(NAEC_Bsal+1)),  LtrlSp4(1:(NAEC_Ltrl+1)))
    allocate(ApclSp5(1:(NAEC_Apcl+1)),  BsalSp5(1:(NAEC_Bsal+1)),  LtrlSp5(1:(NAEC_Ltrl+1)))
    allocate(ApclSp6(1:(NAEC_Apcl+1)),  BsalSp6(1:(NAEC_Bsal+1)),  LtrlSp6(1:(NAEC_Ltrl+1)))
    allocate(ApclSp7(1:(NAEC_Apcl+1)),  BsalSp7(1:(NAEC_Bsal+1)),  LtrlSp7(1:(NAEC_Ltrl+1)))
    allocate(ApclSp8(1:(NAEC_Apcl+1)),  BsalSp8(1:(NAEC_Bsal+1)),  LtrlSp8(1:(NAEC_Ltrl+1)))
    allocate(ApclSp9(1:(NAEC_Apcl+1)),  BsalSp9(1:(NAEC_Bsal+1)),  LtrlSp9(1:(NAEC_Ltrl+1)))
    allocate(ApclSp10(1:(NAEC_Apcl+1)), BsalSp10(1:(NAEC_Bsal+1)), LtrlSp10(1:(NAEC_Ltrl+1)))
    allocate(ApclSp11(1:(NAEC_Apcl+1)), BsalSp11(1:(NAEC_Bsal+1)), LtrlSp11(1:(NAEC_Ltrl+1)))
    allocate(ApclSp0(1:(NAEC_Apcl+1)),  BsalSp0(1:(NAEC_Bsal+1)))
    
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp1,BsalSp1,LtrlSp1)
    cell_val=Hlf_Ncell-1  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    cell_val=Hlf_Ncell-2  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    cell_val=Hlf_Ncell-3  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp4,BsalSp4,LtrlSp4)
    cell_val=Hlf_Ncell-4  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp5,BsalSp5,LtrlSp5)
    cell_val=Hlf_Ncell-5  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp6,BsalSp6,LtrlSp6)
    cell_val=Hlf_Ncell-6  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp7,BsalSp7,LtrlSp7)
    cell_val=Hlf_Ncell-7  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp8,BsalSp8,LtrlSp8)
    cell_val=Hlf_Ncell-8  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp9,BsalSp9,LtrlSp9)
    cell_val=Hlf_Ncell-9  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp10,BsalSp10,LtrlSp10)
    cell_val=Hlf_Ncell-10 ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp11,BsalSp11,LtrlSp11)
    
    if (CellsMeet==0)   call get_ApclBsal_singl_cell_NIsys_with_CM0(ApclSp0,BsalSp0)
    if (CellsMeet.gt.0) call get_Bsal_singl_cell_NIsys_with_CM_gt_0(BsalSp0)
    
    write(*,*) ApclSp1, "Apcl1", BsalSp1, "Bsal1", LtrlSp1,"Ltrl1"
    write(*,*) ApclSp2, "Apcl2", BsalSp2, "Bsal2", LtrlSp2,"Ltrl2"
    write(*,*) ApclSp3, "Apcl3", BsalSp3, "Bsal3", LtrlSp3,"Ltrl3"
    write(*,*) ApclSp4, "Apcl4", BsalSp4, "Bsal4", LtrlSp4,"Ltrl4"
    write(*,*) ApclSp5, "Apcl5", BsalSp5, "Bsal5", LtrlSp5,"Ltrl5"
    write(*,*) ApclSp6, "Apcl6", BsalSp6, "Bsal6", LtrlSp6,"Ltrl6"
    write(*,*) ApclSp7, "Apcl7", BsalSp7, "Bsal7", LtrlSp7,"Ltrl7"
    write(*,*) ApclSp8, "Apcl8", BsalSp8, "Bsal8", LtrlSp8,"Ltrl8"
    write(*,*) ApclSp9, "Apcl9", BsalSp9, "Bsal9", LtrlSp9,"Ltrl9"
    write(*,*) ApclSp10,"Apcl10",BsalSp10,"Bsal10",LtrlSp10,"Ltrl10"
    write(*,*) ApclSp11,"Apcl11",BsalSp11,"Bsal11",LtrlSp11,"Ltrl11"
    
    if (CellsMeet==0)   write(*,*) ApclSp0,"Apcl0",BsalSp0,"Bsal0"
    if (CellsMeet.gt.0) write(*,*) BsalSp0,"Bsal0"
    
    readTheFile = 0    !!!!!!!!!!!!!!!!!!!!!!!
    
    if (forCellsMeet==1 .and. PCSorICS==1) then
       caseVal=1 ; cnt1=10  ; cnt2=4
       
    elseif (forCellsMeet==2 .and. PCSorICS==2) then !(cn1=10;cnt2=7)
       caseVal=2 ; cnt1=8   ; cnt2=0
       
    elseif (forCellsMeet==2 .and. PCSorICS==1) then
       caseVal=3 ; cnt1=5   ; cnt2=0
       
    elseif (forCellsMeet==3 .and. PCSorICS==2) then
       caseVal=4 ; cnt1=3  ; cnt2=0
       
    elseif (forCellsMeet==3 .and. PCSorICS==1) then
       caseVal=5 ; cnt1=5  ; cnt2=0
       
    elseif (forCellsMeet==4 .and. PCSorICS==2) then
       caseVal=6 ; cnt1=11  ; cnt2=0
       
    endif
    
    if (readTheFile==0) then
       
       if (cnt1==0) then
          continue
       elseif (cnt1.gt.0) then
          
          if (caseVal.le.2) then
             
             do cnt = 1,cnt1
                hmL0Bs1=0.95d0 ; call manplt_spr_NIsys(BsalSp1,NAEC_Bsal,hmL0Bs1)
                hmL0Bs2=0.95d0 ; call manplt_spr_NIsys(BsalSp2,NAEC_Bsal,hmL0Bs2)
                hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
                hmL0Bs0=0.95d0 ; call manplt_spr_of_singl_cell_NIsys(BsalSp0,NAEC_Bsal,hmL0Bs0)
                call Equilibrate_only_NI_model
             enddo
             
          elseif (caseVal==3) then
             
             do cnt = 1,cnt1
                hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
                hmL0Bs4=0.95d0 ; call manplt_spr_NIsys(BsalSp4,NAEC_Bsal,hmL0Bs4)
                hmL0Bs5=0.95d0 ; call manplt_spr_NIsys(BsalSp5,NAEC_Bsal,hmL0Bs5)
                hmL0Bs6=0.95d0 ; call manplt_spr_NIsys(BsalSp6,NAEC_Bsal,hmL0Bs6)
                call Equilibrate_only_NI_model
             enddo
             
          elseif (caseval==4) then
             
             do cnt = 1,cnt1   
                hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
                hmL0Bs4=0.95d0 ; call manplt_spr_NIsys(BsalSp4,NAEC_Bsal,hmL0Bs4)
                hmL0Bs5=0.95d0 ; call manplt_spr_NIsys(BsalSp5,NAEC_Bsal,hmL0Bs5)
                hmL0Bs6=0.95d0 ; call manplt_spr_NIsys(BsalSp6,NAEC_Bsal,hmL0Bs6)
                hmL0Bs7=0.95d0 ; call manplt_spr_NIsys(BsalSp7,NAEC_Bsal,hmL0Bs7)
                call Equilibrate_only_NI_model
             enddo
             
          elseif (caseval==5) then
             
             do cnt = 1,cnt1   
                hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
                hmL0Bs4=0.95d0 ; call manplt_spr_NIsys(BsalSp4,NAEC_Bsal,hmL0Bs4)
                hmL0Bs5=0.95d0 ; call manplt_spr_NIsys(BsalSp5,NAEC_Bsal,hmL0Bs5)
                hmL0Bs6=0.95d0 ; call manplt_spr_NIsys(BsalSp6,NAEC_Bsal,hmL0Bs6)
                hmL0Bs7=0.95d0 ; call manplt_spr_NIsys(BsalSp7,NAEC_Bsal,hmL0Bs7)
                call Equilibrate_only_NI_model
             enddo
             
          elseif (caseval==6) then
             
             do cnt = 1,cnt1   
                hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
                hmL0Bs4=0.95d0 ; call manplt_spr_NIsys(BsalSp4,NAEC_Bsal,hmL0Bs4)
                hmL0Bs5=0.95d0 ; call manplt_spr_NIsys(BsalSp5,NAEC_Bsal,hmL0Bs5)
                hmL0Bs6=0.95d0 ; call manplt_spr_NIsys(BsalSp6,NAEC_Bsal,hmL0Bs6)
                hmL0Bs7=0.95d0 ; call manplt_spr_NIsys(BsalSp7,NAEC_Bsal,hmL0Bs7)
                call Equilibrate_only_NI_model
             enddo
             
          endif
          
             
       endif
       
          
       if (cnt2==0) then
          continue
       elseif (cnt2.gt.0) then
          
          if (caseVal.le.2) then
             
             do cnt = 1,cnt2
                hmL0Ap2=1.05d0 ; call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap2)
                hmL0Ap3=1.05d0 ; call manplt_spr_NIsys(ApclSp3,NAEC_Apcl,hmL0Ap3)
                
                call Equilibrate_only_NI_model
             enddo
             
          elseif (caseVal==3) then
             
             do cnt = 1,cnt2
                hmL0Ap3=1.05d0 ; call manplt_spr_NIsys(ApclSp3,NAEC_Apcl,hmL0Ap3)
                hmL0Ap4=1.05d0 ; call manplt_spr_NIsys(ApclSp4,NAEC_Apcl,hmL0Ap4)
                hmL0Ap5=1.05d0 ; call manplt_spr_NIsys(ApclSp5,NAEC_Apcl,hmL0Ap5)
                hmL0Ap6=1.05d0 ; call manplt_spr_NIsys(ApclSp6,NAEC_Apcl,hmL0Ap6)
                
                call Equilibrate_only_NI_model
             enddo
             
          endif
          
       endif
       
    elseif (readTheFile==1) then
       
       write(*,*) Frame_NI,"Frame_NI bfr increasing"
       
       FrmNoIncr = cnt1+cnt2
       Frame_NI  = Frame_NI+FrmNoIncr-1
       
       write(*,*) Frame_NI,"Frame_NI aft increasing"
       
       call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
       call Equilibrate_only_NI_model
       
    endif
    
  end subroutine tilt_reduction_of_the_boundry
  
  
  subroutine tilt_reduction_general_purpose
    implicit none
    
  end subroutine tilt_reduction_general_purpose
  
  
  subroutine print_symm_prperties
    implicit none
    integer :: i,j,jmax
    integer :: nsprsInACell
    integer :: Hlf_NsprV
    integer :: sprL,sprR
    integer :: cellL,cellR
    real*8  :: dksV,dl0V
    real*8  :: dkaV,dA0V
    
    if (modelID.ne.2) stop 'only for modelID=2'
    
    open(unit=541,file='print_symm_prperties.dat')
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    Hlf_NsprV    = Hlf_Ncell*nsprsInACell
    
    write(541,*) Hlf_Ncell,Hlf_NsprV,"Hlf_Ncell-Hlf_Nspr"
    
    do i = 1,3
       
       if (i==1) jmax=Hlf_NsprV
       if (i==2) jmax=Hlf_Ncell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             sprL = j           ; sprR = sprL + Hlf_NsprV
             dksV = k_spr(sprL) - k_spr(sprR)
             dl0V = l0(sprL)    - l0(sprR)
             write(541,*) sprL,sprR,k_spr(sprL),k_spr(sprR),dksV,l0(sprL),l0(sprR),dl0V,"spr Prp" 
             
             if (j==jmax) write(541,*) (sprR+1),(sprR+2),(sprR+3),"N_spr"
             if (j==jmax) write(541,*) k_spr(sprR+1),k_spr(sprR+2),k_spr(sprR+3),l0(sprR+1),l0(sprR+2),l0(sprR+3)
             
          elseif (i==2) then
             
             cellL = j             ; cellR = cellL + Hlf_Ncell
             dkaV  = k_area(cellL) - k_area(cellR)
             dA0V  = A0(cellL)     - A0(cellR)
             write(541,*)cellL,cellR,k_area(cellL),k_area(cellR),dkaV,A0(cellL),A0(cellR),dA0V,"cell"
             
             if (j==jmax) write(541,*) (cellR+1),N_cell
             if (j==jmax) write(541,*) k_area(cellR+1),A0(cellR+1),"IC"
             
          elseif (i==3) then
             
             write(541,*) CgXNode(j),CgYNode(j),k_phi(j,1:max_Phi_node),j
             
          endif
          
       enddo
       
       write(541,*) " "
       
    enddo
    
    close(541)
    
  end subroutine print_symm_prperties
  
  subroutine CM1_match_curvatureWise
    implicit none
    integer, allocatable :: Apcls1(:),Apcls2(:) 
    integer, allocatable :: Bsals1(:),Bsals2(:)
    integer, allocatable :: Ltrlsb1(:),Ltrlsb2(:)
    integer, allocatable :: Ltrlsa1(:),Ltrlsa2(:)
    
    integer :: cell1L,cell2L,cell3L,cell4L
    integer :: cell1R,cell2R,cell3R,cell4R
    integer :: NcellToFind
    real*8  :: hmuchA0,hmuchL0
    real*8  :: A0Decr,A0Incr,L0Decr,L0Incr
    integer :: cnt,cnt1,cnt2
    integer :: incrORdecr,cntA
    integer :: cntmax,cnt1max,cnt2max
    
    integer :: nsprsInACell
    integer :: sprNmL,sprNmR
    
    integer, allocatable :: cellsToFind(:)
    integer, allocatable :: ApclsL(:,:),BsalsL(:,:),LtrlsBL(:,:),LtrlsAL(:,:)
    
    cell1L=11 ; cell2L=cell1L-1 ; cell3L=cell1L-2 ; cell4L=cell1L-3
    
    cell1R=cell1L+Hlf_Ncell     ; cell2R=cell2L+Hlf_Ncell
    cell3R=cell3L+Hlf_Ncell     ; cell4R=cell4L+Hlf_Ncell
    
    write(*,*) cell1L,cell2L,cell3L,cell4L,"cell1L to cell4L"
    write(*,*) cell1R,cell2R,cell3R,cell4R,"cell1R to cell4R"
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    write(*,*) nsprsInACell,"nsprsInACell"
    
    NcellToFind=4 ; allocate(cellsToFind(1:NcellToFind))
    cellsToFind(1)=cell1L ; cellsToFind(2)=cell2L
    cellsToFind(3)=cell3L ; cellsToFind(4)=cell4L
    
    allocate(ApclsL(1:NcellToFind,1:(NAEC_Apcl+1)), BsalsL(1:NcellToFind,1:(NAEC_Bsal+1)))
    allocate(LtrlsBL(1:NcellToFind,1:(NAEC_Ltrl+1)),LtrlsAL(1:NcellToFind,1:(NAEC_Ltrl+1)))
    
    call get_the_apcls_bsals_and_both_ltrls(NcellToFind,cellsToFind,ApclsL,BsalsL,LtrlsBL,LtrlsAL)
    incrORdecr=-1 ; cntA=4 ; call Press_modfctn_and_maintain_apprx_lngth(cell1L,incrORdecr,cntA)
    
    do cnt = 1,5
       A0(23)=0.98d0*A0(23)
       call Equilibrate_only_NI_model
    enddo
    
  end subroutine CM1_match_curvatureWise
  
  
  subroutine flatten_the_IC_and_maintain_Area_and_Length
    implicit none
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp1(:),LtrlSp2(:)
    integer :: nsprInACell
    integer :: i,j,imax,jmax
    
    if (modelID==1) then
       write(*,*) "not for modelID=1",modelID,"sb: flatten_the_IC_and_maintain_Area_and_Length"
       stop
    endif
    
    nsprInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)))
    allocate(LtrlSp1(1:(NAEC_Ltrl+1)),LtrlSp2(1:(NAEC_Ltrl+1)))
    
    !do i = 1,4
    !   if (i==1) jmax=NAEC_Apcl+1
    !   do j = 1,jmax
    !      ApclSp(j) = (N_cell-1)*(nsprInACell) + 1
    !   enddo
    !enddo
    
  end subroutine flatten_the_IC_and_maintain_Area_and_Length
  
  subroutine matching_cell_shape_CM1
    implicit none
    integer :: nsprsInACell
    integer :: N_cellsToAlter,cell1,cell2,cell3,cell4
    integer :: ApclSprs
    integer :: i,j
    integer :: whichOne
    
    integer, allocatable :: ApclsL(:,:),BsalsL(:,:),LtrlsL(:,:)
    integer, allocatable :: cellToAlter(:)
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    call allocate_NrmVars
    call read_the_normlzd_lngths_diff_in_SimExp
    
    N_cellsToAlter=4  ; allocate(cellToAlter(1:N_cellsToAlter))
    cellToAlter(1)=8  ; cellToAlter(2)=9
    cellToAlter(3)=10 ; cellToAlter(4)=11
    
    allocate(ApclsL(1:N_cellsToAlter,1:(NAEC_Apcl+1))) ; ApclsL=-1
    allocate(BsalsL(1:N_cellsToAlter,1:(NAEC_Bsal+1))) ; BsalsL=-1
    allocate(LtrlsL(1:N_cellsToAlter,1:(NAEC_Ltrl+1))) ; LtrlsL=-1
    
    call get_the_apcls_bsals_ltrls(N_cellsToAlter,cellToAlter,ApclsL,BsalsL,LtrlsL)
    
    whichOne=1;call alter_the_apcls_or_bsals_or_ltrls(N_cellsToAlter,cellToAlter,whichOne,NAEC_Apcl,ApclsL)
    whichOne=2;call alter_the_apcls_or_bsals_or_ltrls(N_cellsToAlter,cellToAlter,whichOne,NAEC_Bsal,BsalsL)
    whichOne=3;call alter_the_apcls_or_bsals_or_ltrls(N_cellsToAlter,cellToAlter,whichOne,NAEC_Ltrl,LtrlsL)
    
  end subroutine matching_cell_shape_CM1
  
  subroutine allocate_NrmVars
    implicit none
    
    allocate(NrmLa_Sim(1:Hlf_Ncell),NrmLa_Exp(1:Hlf_Ncell),diffNrm_La(1:Hlf_Ncell))
    allocate(NrmLb_Sim(1:Hlf_Ncell),NrmLb_Exp(1:Hlf_Ncell),diffNrm_Lb(1:Hlf_Ncell))
    allocate(NrmLl_Sim(1:Hlf_Ncell),NrmLl_Exp(1:Hlf_Ncell),diffNrm_Ll(1:Hlf_Ncell))
    allocate(NrmAv_Sim(1:Hlf_Ncell),NrmAv_Exp(1:Hlf_Ncell),diffNrm_Av(1:Hlf_Ncell))
    
    NrmLa_Sim=-10.0d20 ; NrmLa_Exp=-10.0d20 ; diffNrm_La=-10.0d20
    NrmLb_Sim=-10.0d20 ; NrmLb_Exp=-10.0d20 ; diffNrm_Lb=-10.0d20
    NrmLl_Sim=-10.0d20 ; NrmLl_Exp=-10.0d20 ; diffNrm_Ll=-10.0d20
    NrmAv_Sim=-10.0d20 ; NrmAv_Exp=-10.0d20 ; diffNrm_Av=-10.0d20
    
  end subroutine allocate_NrmVars
  
  subroutine read_the_normlzd_lngths_diff_in_SimExp
    implicit none
    integer :: i,j
    real*8  :: dum1,dum3,dum4
    integer :: cellNm
    
    write(*,*) CellsMeet,"CellsMeet"
    
    open(unit=32,file='CM1_Diff_SIM_EXP.dat')
    
    do i = 1,4
       
       do j = 1,Hlf_Ncell
          
          if (i==1) read(32,*) dum1,simL1a,NrmLa_Sim(j),dum3,dum4,NrmLa_Exp(j),diffNrm_La(j),cellNm
          if (i==2) read(32,*) dum1,simL1b,NrmLb_Sim(j),dum3,dum4,NrmLb_Exp(j),diffNrm_Lb(j),cellNm
          if (i==3) read(32,*) dum1,simL1l,NrmLl_Sim(j),dum3,dum4,NrmLl_Exp(j),diffNrm_Ll(j),cellNm
          if (i==4) read(32,*) dum1,simA1, NrmAv_Sim(j),dum3,dum4,NrmAv_Exp(j),diffNrm_Av(j),cellNm
          
       enddo
       
    enddo
    
    close(32)
    
  end subroutine read_the_normlzd_lngths_diff_in_SimExp
  
  
  subroutine alter_the_apcls_or_bsals_or_ltrls(N_cellsToAlter,cellNms,whichOne,NAEC_val,sprNms)
    implicit none
    integer, intent(in) :: N_cellsToAlter,whichOne,NAEC_val
    integer, intent(in) :: cellNms(1:N_cellsToAlter)
    integer, intent(in) :: sprNms(1:N_cellsToAlter,1:(NAEC_val+1))
    
    real*8  :: diffNrm 
    integer :: i,j,sprNm
    real*8  ::  hmuch(1:N_cellsToAlter)
    integer :: sprNmIncrForRight
    integer :: cntEndLp

    write(*,*) " "
    write(*,*) " "
    
    sprNmIncrForRight = (Hlf_Ncell)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    write(*,*) sprNmIncrForRight,"sprNmIncr"
    
    cntEndLp = 0
    
    do 
       cntEndLp = cntEndLp+1
       
       do i = 1,N_cellsToAlter
          if (whichOne==1) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"apcls before for i=",i
          if (whichOne==2) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"bsals before for i=",i
          if (whichOne==3) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"ltrls before for i=",i
       enddo
       
       
       do i = 1,N_cellsToAlter
          
          if (whichOne==1) then
             diffNrm = diffNrm_La(cellNms(i))
             write(*,*) diffNrm,i,"diffNrm for La"
          elseif (whichOne==2) then
             diffNrm = diffNrm_Lb(cellNms(i))
             write(*,*) diffNrm,i,"diffNrm for Lb"
          elseif (whichOne==3) then
             diffNrm = diffNrm_Ll(cellNms(i))
             write(*,*) diffNrm,i,"diffNrm for Ll"
          endif
          
          if (abs(diffNrm) .gt. tol_NrmL) then
             if (diffNrm .le. ZERO) hmuch(i) = 1.05d0
             if (diffNrm .gt. ZERO) hmuch(i) = 0.95d0
          elseif (abs(diffNrm) .le. tol_NrmL) then
             hmuch(i) = 1.00d0
          endif
          
          write(*,*) hmuch(i),i,"hmuch,i"
          
          do j = 1,(NAEC_val+1)
             sprNm                       = sprNms(i,j)
             l0(sprNm)                   = hmuch(i)*l0(sprNm)
             l0(sprNm+sprNmIncrForRight) = l0(sprNm)
          enddo
          
       enddo
       
       do i = 1,N_cellsToAlter
          if (whichOne==1) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"apcls after for i=",i
          if (whichOne==2) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"bsals after for i=",i
          if (whichOne==3) write(*,*) l0(sprNms(i,1:(NAEC_val+1))),"ltrls after for i=",i
       enddo
       
       call Equilibrate_only_NI_model
       
       if (whichOne==1) call calclt_the_sim_norm(N_cellsToAlter,cellNms,whichOne,NAEC_Apcl)
       if (whichOne==2) call calclt_the_sim_norm(N_cellsToAlter,cellNms,whichOne,NAEC_Bsal)
       if (whichOne==3) call calclt_the_sim_norm(N_cellsToAlter,cellNms,whichOne,NAEC_Ltrl)
       
       if (whichOne==1) then
          if (cntEndLp==3) exit
       elseif (whichOne==2) then
          if (cntEndLp==5) exit
       elseif (whichOne==3) then
          if (cntEndLp==5) exit
       endif
       
    enddo
    
  end subroutine alter_the_apcls_or_bsals_or_ltrls
  
  
  
  subroutine Press_modfctn_and_maintain_apprx_lngth(cellNm,incrORdecr,cntA)
    implicit none
    integer, intent(in) :: cellNm
    integer, intent(in) :: incrORdecr,cntA
    
    real*8  :: hmuchA0,hmuchL0
    real*8  :: A0mod,L0mod
    integer :: NcellToFind,nsprsInACell
    integer :: cnt,cnt1,cnt2,cnt2max
    integer :: sprNmL,sprNmR
    integer :: cntL
    
    integer, allocatable :: cellsToFind(:)
    integer, allocatable :: ApclsL(:,:),BsalsL(:,:),LtrlsBL(:,:),LtrlsAL(:,:)
    
    hmuchA0 = 0.04d0 ; hmuchL0 = hmuchA0/2.0d0
    
    A0mod = 1.00d0 + (incrORdecr*(+1.00d0)) * (hmuchA0)
    L0mod = 1.00d0 + (incrORdecr*(-1.00d0)) * (hmuchL0)
    
    write(*,*) A0mod,L0mod,"A0mod and L0mod"
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    NcellToFind  = 1 ; allocate(cellsToFind(1:NcellToFind))
    cellsToFind(1)=cellNm
    
    allocate(ApclsL(1:NcellToFind,1:(NAEC_Apcl+1)), BsalsL(1:NcellToFind,1:(NAEC_Bsal+1)))
    allocate(LtrlsBL(1:NcellToFind,1:(NAEC_Ltrl+1)),LtrlsAL(1:NcellToFind,1:(NAEC_Ltrl+1)))
    
    call get_the_apcls_bsals_and_both_ltrls(NcellToFind,cellsToFind,ApclsL,BsalsL,LtrlsBL,LtrlsAL)
    
    write(*,*) ApclsL(1,1:(NAEC_Apcl+1)), "apcls"   ,BsalsL(1,1:(NAEC_Bsal+1)), "bsals"
    write(*,*) LtrlsBL(1,1:(NAEC_Ltrl+1)),"ltrlsBfr",LtrlsAL(1,1:(NAEC_Ltrl+1)),"ltrlsAft"
    
    do cnt = 1,cntA
       A0(cellNm)=(A0mod)*A0(cellNm) ; A0(cellNm+Hlf_Ncell)=A0(cellNm) ! cellNm=11
       call Equilibrate_only_NI_model
    enddo
    
    cntL = cntA/2 ; write(*,*) cntL,cntA,"cntL-cntA"
    
    do cnt = 1,cntL
       
       do cnt1 = 1,4
          
          if (cnt1==1) cnt2max = NAEC_Apcl+1 
          if (cnt1==2) cnt2max = NAEC_Bsal+1
          if (cnt1==3) cnt2max = NAEC_Ltrl+1
          if (cnt1==4) cnt2max = NAEC_Ltrl+1
          
          do cnt2 = 1,cnt2max
             
             if (cnt1==1) sprNmL = ApclsL(1,cnt2)
             if (cnt1==2) sprNmL = BsalsL(1,cnt2) 
             if (cnt1==3) sprNmL = LtrlsBL(1,cnt2)
             if (cnt1==4) sprNmL = LtrlsAL(1,cnt2)
             
             sprNmR     = sprNmL+(Hlf_Ncell*nsprsInACell)
             write(*,*) sprNmL,sprNmR,"sprNmLR"
             
             l0(sprNmL) = L0mod*l0(sprNmL) ; l0(sprNmR) = l0(sprNmL)
             
          enddo
          
       enddo
       
       call Equilibrate_only_NI_model
       
    enddo
    
  end subroutine Press_modfctn_and_maintain_apprx_lngth
  
  
  subroutine make_cells_identical_to_frst_cell(EndcellNum)
    implicit none
    integer, intent(in) :: EndcellNum
    integer             :: nsprsInACell,sprNm,cell1,cellV
    integer             :: cnt1,cnt2,i,j,jmax,j1,j1max
    integer             :: sprStr,sprEnd,sprNmAdd
    
    real*8              :: A0_cell1
    real*8, allocatable :: l0_apclC1(:),l0_bsalC1(:),l0_ltrlC1(:)    
    
    cell1        = 1
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    allocate(l0_apclC1(1:(NAEC_Apcl+1)),l0_bsalC1(1:(NAEC_Bsal+1)),l0_ltrlC1(1:(NAEC_ltrl+1)))
    
    A0_cell1 = A0(cell1)
    write(*,*) A0_cell1,"cell1"
    
    do i = 1,3
       
       if (i==1) jmax=NAEC_Apcl+1
       if (i==2) jmax=NAEC_Bsal+1
       if (i==3) jmax=NAEC_Ltrl+1
       
       do j = 1,jmax
          
          if (i==1) then
             sprNm        = (cell1-1)*(nsprsInACell) + j
             l0_apclC1(j) = l0(sprNm)
             write(*,*) l0_apclC1(j),j,"l0-apcl"
             
          elseif (i==2) then
             sprNm        = (cell1-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
             l0_bsalC1(j) = l0(sprNm)
             write(*,*) l0_bsalC1(j),j,"l0-bsal"
             
          elseif (i==3) then
             sprNm        = (cell1-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             l0_ltrlC1(j) = l0(sprNm)
             write(*,*) l0_ltrlC1(j),j,"l0-ltrl"
          endif
          
       enddo
       
    enddo
    
    sprNmAdd = (Hlf_Ncell)*(nsprsInACell)
    
    do i = 2,EndcellNum
       
       cellV = i
       A0(i) = A0_cell1
       
       do j = 1,3
          
          if (j==1) j1max=NAEC_Apcl+1
          if (j==2) j1max=NAEC_Bsal+1
          if (j==3) j1max=NAEC_Ltrl+1
          
          do j1 = 1,j1max
             
             if (j==1) then
                sprNm              = (cellV-1)*(nsprsInACell) + j1
                l0(sprNm)          = l0_apclC1(j1)
                l0(sprNm+sprNmAdd) = l0(sprNm)
                
             elseif (j==2) then
                sprNm              = (cellV-1)*(nsprsInACell) + (NAEC_Apcl+1) + j1
                l0(sprNm)          = l0_bsalC1(j1)
                l0(sprNm+sprNmAdd) = l0(sprNm)
                
             elseif (j==3) then
                sprNm              = (cellV-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j1
                l0(sprNm)          = l0_ltrlC1(j1)
                l0(sprNm+sprNmAdd) = l0(sprNm)
             endif
             
          enddo
       enddo
       
       sprStr=(cellV-1)*(nsprsInACell)+1 ; sprEnd=(cellV)*(nsprsInACell)
       write(*,*) l0(sprStr:sprEnd),"l0 value"
       
    enddo
    
    call Equilibrate_only_NI_model
    
  end subroutine make_cells_identical_to_frst_cell
  
  
  subroutine reduce_force_stepbystep(cntMax)
    implicit none
    integer, intent(in) :: cntMax
    integer             :: force_ANL,force_ANR,cnt
    
    force_ANL=23 ; force_ANR=47
    
    do cnt = 1,cntMax
       CgYNode(force_ANL)=CgYNode(force_ANL)-0.05d0
       CgYNode(force_ANR)=CgYNode(force_ANL)
    
       call Equilibrate_only_NI_model
    enddo
    
  end subroutine reduce_force_stepbystep
  
  
  subroutine calclt_the_sim_norm(N_cellsToCalculte,cellNms,whichOne,NAEC_val)
    implicit none
    integer, intent(in) :: N_cellsToCalculte,whichOne,NAEC_val
    integer, intent(in) :: cellNms(1:N_cellsToCalculte)
    
    integer :: cellNm
    
    do i = 1,N_cellsToCalculte
       
       cellNm = cellNms(i)
       
       if (whichOne==1) then
          NrmLa_Sim(cellNm)  = len_add(cellNm,whichOne,NAEC_val)/simL1a 
          diffNrm_La(cellNm) = NrmLa_Sim(cellNm)-NrmLa_Exp(cellNm)
          write(*,*) NrmLa_Sim(cellNm),NrmLa_Exp(cellNm),diffNrm_La(cellNm),cellNm,"NrmLa_Sim"
          
       elseif (whichOne==2) then
          NrmLb_Sim(cellNm)  = len_add(cellNm,whichOne,NAEC_val)/simL1b
          diffNrm_Lb(cellNm) = NrmLb_Sim(cellNm)-NrmLb_Exp(cellNm)
          write(*,*) NrmLb_Sim(cellNm),NrmLb_Exp(cellNm),diffNrm_Lb(cellNm),cellNm,"NrmLb_Sim"
          
       elseif (whichOne==3) then
          NrmLl_Sim(cellNm)  = len_add(cellNm,whichOne,NAEC_val)/simL1l
          diffNrm_Ll(cellNm) = NrmLl_Sim(cellNm)-NrmLl_Exp(cellNm)
          write(*,*) NrmLl_Sim(cellNm),NrmLl_Exp(cellNm),diffNrm_Ll(cellNm),cellNm,"NrmLl_Sim"
       endif
       
    enddo
    
  end subroutine calclt_the_sim_norm
  
  
  real*8 function len_Add(cellNm,whichOne,NAEC_val)
    implicit none
    integer, intent(in) :: cellNm,whichOne,NAEC_val
    integer             :: i,nsprsInACell
    integer             :: sprNm
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    len_Add      = 0.00d0
    
    do i = 1,(NAEC_val+1)
       
       if (whichOne==1) sprNm = (cellNm-1)*(nsprsInACell) + i
       if (whichOne==2) sprNm = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + i
       if (whichOne==3) sprNm = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       
       write(*,*) sprNm,i,"sprNm in len_Add function"
       len_Add = len_Add+l(sprNm)
       
    enddo
    
    write(*,*) len_add,"len_add value"
    
  end function len_Add
  
  
  
  subroutine manpltns_neededtodo_for_cntrlStates(EXEorREAD,CMv,PCSorICS)
    implicit none
    integer, intent(in) :: EXEorREAD
    integer, intent(in) :: CMv,PCSorICS
    
    integer :: cntMax
    integer :: cellNm,EndCellNum,ApBsLt
    real*8  :: hmL0
    
    if (EXEorREAD == 1) then
       
       if ((CMv==1) .and. (PCSorICS==2)) then ! PCSorICS=2 comes first
          continue
       elseif ((CMv==1) .and. (PCSorICS==1)) then ! then PCSorICS=1 comes
          
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          call matching_cell_shape_CM1 !(Frm aft this = 40, look into /CM1PCS directory)
          
          call CM1_match_curvatureWise
          cntMax=8      ; call reduce_force_stepbystep(cntMax)
          EndCellNum=10 ; call make_cells_identical_to_frst_cell(EndCellNum)!(Frm No=59, /CM1PCS)
          
          EndCellNum=9 ; call make_cells_identical_to_frst_cell(EndCellNum)
          cellNm=11 ; ApBsLt=1 ; cntMax=4 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          ApBsLt=2 ; cntMax=5 ; hmL0=1.05d0
          call manplt_spr_NIsys_for_Ncell(ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          !%%%%%%%%%  upto this ; AFTER this perturbation; not NEEDED in actual PCS/ICS %%%%%%%%
          
          cellNm=9 ; ApBsLt=3 ; cntMax=10 ; hmL0=0.95d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
       elseif ((CMv==2) .and. (PCSorICS==2)) then
          
          cellNm=11 ; ApBsLt=3 ; cntMax=5 ; hmL0=0.95d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=10 ; ApBsLt=3 ; cntMax=5 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          
          cellNm=10 ; ApBsLt=1 ; cntMax=3 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
       elseif ((CMv==2) .and. (PCSorICS==1)) then   
          
          cellNm=10 ; ApBsLt=3 ; cntMax=5 ; hmL0=0.90d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          
       elseif ((CMv==3) .and. (PCSorICS==2)) then
          
          cellNm=11 ; ApBsLt=3 ; cntMax=5 ; hmL0=0.95d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=10 ; ApBsLt=3 ; cntMax=5 ; hmL0=0.95d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=9  ; ApBsLt=3 ; cntMax=5 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          cellNm=9  ; ApBsLt=1 ; cntMax=10 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          
       elseif ((CMv==3) .and. (PCSorICS==1)) then
          
          cellNm=9 ; ApBsLt=3 ; cntMax=3 ; hmL0=0.90d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          
       elseif ((CMv==4) .and. (PCSorICS==2)) then
          
          cellNm=11 ; ApBsLt=1 ; cntMax=7 ; hmL0=0.90d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          cellNm=9 ; ApBsLt=3 ; cntMax=5 ; hmL0=0.95d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=8 ; ApBsLt=3 ; cntMax=5 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          cellNm=10 ; ApBsLt=1 ; cntMax=7 ; hmL0=0.90d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=9 ; ApBsLt=1 ; cntMax=7 ; hmL0=0.90d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          cellNm=8 ; ApBsLt=1 ; cntMax=3 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          cellNm=7 ; ApBsLt=1 ; cntMax=3 ; hmL0=1.05d0
          call manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,Exprmnt_NI,Frame_NI)
          
          call tilt_reduction_of_the_boundry(CMv,PCSorICS)
          
       endif
       
    elseif (EXEorREAD == 2) then
       
       write(*,*) Frame_NI,"Frame_NI val 1"
       call read_props_from_CMdataFiles(CMv,PCSorICS)
       call Equilibrate_only_NI_model
       
    endif
    
  end subroutine manpltns_neededtodo_for_cntrlStates
  
  subroutine continue_simln_frm_prp
    implicit none
    
    Frame_NI=45 ; call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI,"Frm - NI"
    
  end subroutine continue_simln_frm_prp
  
  subroutine read_props_from_CMdataFiles(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    character(len=100)  :: prfS,prfA,prfN,suff
    character(len=100)  :: flnmS,flnmA,flnmN
    
    integer :: i,j,jmax
    integer :: sprNm,cellNm,nodeNm
    real*8  :: TnsnV,ksV,lV,l0V
    real*8  :: PresV,kaV,AV,A0V
    integer :: nodeTyp
    real*8  :: CgX,CgY
    
    prfS='sprsPrps_CM'
    prfA='cellPrps_CM'
    prfN='nodePrps_CM'
    
    if (PCSorICS==1) suff='PCS.dat'
    if (PCSorICS==2) suff='ICS.dat'
    
    write(flnmS,'(a,i1.1,a)')trim(adjustl(prfS)),CMv,trim(adjustl(suff))
    write(flnmA,'(a,i1.1,a)')trim(adjustl(prfA)),CMv,trim(adjustl(suff))
    write(flnmN,'(a,i1.1,a)')trim(adjustl(prfN)),CMv,trim(adjustl(suff))
    
    write(*,*)trim(adjustl(flnmS))
    write(*,*)trim(adjustl(flnmA))
    write(*,*)trim(adjustl(flnmN))
    
    open(unit=731,file=trim(adjustl(flnmS)))
    open(unit=732,file=trim(adjustl(flnmA)))
    open(unit=733,file=trim(adjustl(flnmN)))
    
    do i = 1,3
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) then
             read(731,*) TnsnV,ksV,lV,l0V,sprNm
             write(*,*)  TnsnV,ksV,lV,l0V,sprNm
             k_spr(j) = ksV ; l0(j) = l0V
             
          elseif (i==2) then
             read(732,*) cellNm,PresV,kaV,AV,A0V
             write(*,*)  cellNm,PresV,kaV,AV,A0V
             k_area(j) = kaV ; A0(j) = A0V
             
          elseif (i==3) then
             read(733,*) nodeNm,nodeTyp,CgX,CgY
             write(*,*)  nodeNm,nodeTyp,CgX,CgY
             CgXNode(j) = CgX ; CgYNode(j) = CgY
          endif
          
       enddo
    enddo
    
    close(731)
    close(732)
    close(733)
    
  end subroutine read_props_from_CMdataFiles
  
  
  ! subroutine tilt_reduction_of_the_boundry(forCellsMeet,PCSorICS)
  !   implicit none
  !   integer, intent(in) :: forCellsMeet,PCSorICS
    
  !   integer :: cell_val,cnt
  !   real*8  :: hmL0Ap1,hmL0Bs1,hmL0Ap2,hmL0Bs2,hmL0Ap3,hmL0Bs3,hmL0Ap4,hmL0Bs4
    
  !   integer, allocatable :: ApclSp1(:),BsalSp1(:),LtrlSp1(:)
  !   integer, allocatable :: ApclSp2(:),BsalSp2(:),LtrlSp2(:)
  !   integer, allocatable :: ApclSp3(:),BsalSp3(:),LtrlSp3(:)
  !   integer, allocatable :: ApclSp4(:),BsalSp4(:)
    
  !   integer :: readTheFile
  !   integer :: FrmNoIncr
  !   integer :: cnt1,cnt2
    
  !   allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:(NAEC_Bsal+1)),LtrlSp1(1:(NAEC_Ltrl+1)))
  !   allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:(NAEC_Bsal+1)),LtrlSp2(1:(NAEC_Ltrl+1)))
  !   allocate(ApclSp3(1:(NAEC_Apcl+1)),BsalSp3(1:(NAEC_Bsal+1)),LtrlSp3(1:(NAEC_Ltrl+1)))
  !   allocate(ApclSp4(1:(NAEC_Apcl+1)),BsalSp4(1:(NAEC_Bsal+1)))
    
  !   cell_val=Hlf_Ncell-1 ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp1,BsalSp1,LtrlSp1)
  !   cell_val=Hlf_Ncell-2 ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
  !   cell_val=Hlf_Ncell-0 ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    
  !   if (CellsMeet==0)   call get_ApclBsal_singl_cell_NIsys_with_CM0(ApclSp4,BsalSp4)
  !   if (CellsMeet.gt.0) call get_Bsal_singl_cell_NIsys_with_CM_gt_0(BsalSp4)
    
  !   write(*,*) ApclSp1,"Apcl1",BsalSp1,"Bsal1"
  !   write(*,*) ApclSp2,"Apcl2",BsalSp2,"Bsal2"
  !   write(*,*) ApclSp3,"Apcl3",BsalSp3,"Bsal3"
    
  !   if (CellsMeet==0)   write(*,*) ApclSp4,"Apcl4",BsalSp4,"Bsal4"
  !   if (CellsMeet.gt.0) write(*,*) BsalSp4,"Bsal4"
    
  !   readTheFile = 0 !!!!!!!!!!!!!!!!!!!!!!!
    
  !   if (forCellsMeet==1 .and. PCSorICS==1) then
  !      cnt1=10 ; cnt2=4
       
  !   elseif (forCellsMeet==2 .and. PCSorICS==2) then
  !      !cnt1=10 ; cnt2=7
  !      cnt1=8 ; cnt2=0  ! ICS
       
  !   elseif (forCellsMeet==2 .and. PCSorICS==1) then`
  !      cnt1=8 ; cnt2=0 !PCS 
  !   endif
    
  !   if (readTheFile==0) then
       
  !      if (cnt1==0) then
  !         continue
  !      elseif (cnt1.gt.0) then
  !         do cnt = 1,cnt1
  !            hmL0Bs1=0.95d0 ; call manplt_spr_NIsys(BsalSp1,NAEC_Bsal,hmL0Bs1)
  !            hmL0Bs2=0.95d0 ; call manplt_spr_NIsys(BsalSp2,NAEC_Bsal,hmL0Bs2)
  !            hmL0Bs3=0.95d0 ; call manplt_spr_NIsys(BsalSp3,NAEC_Bsal,hmL0Bs3)
  !            hmL0Bs4=0.95d0 ; call manplt_spr_of_singl_cell_NIsys(BsalSp4,NAEC_Bsal,hmL0Bs4)
  !            call Equilibrate_only_NI_model
  !         enddo
  !      endif
       
          
  !      if (cnt2==0) then
  !         continue
  !      elseif (cnt2.gt.0) then
  !         do cnt = 1,cnt2
  !            hmL0Ap1=1.05d0 ; call manplt_spr_NIsys(ApclSp1,NAEC_Apcl,hmL0Ap1)
  !            hmL0Ap2=1.05d0 ; call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap2)
  !            call Equilibrate_only_NI_model
  !         enddo
  !      endif
       
  !   elseif (readTheFile==1) then
       
  !      write(*,*) Frame_NI,"Frame_NI bfr increasing"
       
  !      FrmNoIncr = cnt1+cnt2
  !      Frame_NI  = Frame_NI+FrmNoIncr-1
       
  !      write(*,*) Frame_NI,"Frame_NI aft increasing"
       
  !      call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
  !      call Equilibrate_only_NI_model
       
  !   endif
    
  ! end subroutine tilt_reduction_of_the_boundry
  
  
  subroutine displace_pulley_virtually_to_check_force
    implicit none
    real*8  :: displcmntofPP
    integer :: pulley_pnt
    real*8  :: grdVlue(1:N_node,1:N_dmnsn),force(1:N_node,1:N_dmnsn)
    real*8  :: grdV1(1:N_node,1:N_dmnsn),grdSprV1(1:N_node,1:N_dmnsn),grdAreaV1(1:N_node,1:N_dmnsn)
    real*8  :: grdGrvtnlV1(1:N_node,1:N_dmnsn),grdBendV1(1:N_node,1:N_dmnsn),grdSRypV1(1:N_node,1:N_dmnsn)
    real*8  :: grdV2(1:N_node,1:N_dmnsn),grdSprV2(1:N_node,1:N_dmnsn) ,grdAreaV2(1:N_node,1:N_dmnsn)
    real*8  :: grdGrvtnlV2(1:N_node,1:N_dmnsn),grdBendV2(1:N_node,1:N_dmnsn),grdSRypV2(1:N_node,1:N_dmnsn)
    
    real*8  :: forceHorzntPP,forceVrtcalPP
    real*8  :: CgXStr(1:N_node),CgYStr(1:N_node)
    
    displcmntofPP         = 0.0700d0                       ; write(*,*) displcmntofPP,"displcmntPP"
    pulley_pnt            = (2*(Hlf_Ncell+1))*2 + 1        ; write(*,*) pulley_pnt,"pp"
    node_xy(pulley_pnt,1) = node_xy(pulley_pnt,1)          ; write(*,*) node_xy(pulley_pnt,1),"pp x"
    node_xy(pulley_pnt,1) = node_xy(pulley_pnt,1)-displcmntofPP ; write(*,*) node_xy(pulley_pnt,1),"pp x"
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call get_gradient_with_changing_nodeTyp_test_withCompnts(node_xy,l0,A0,grdV1,grdSprV1,grdAreaV1,&
         grdGrvtnlV1,grdBendV1,grdSRypV1)
    
    force         = -grdV1
    forceHorzntPP = force(pulley_pnt,1)
    forceVrtcalPP = force(pulley_pnt,2)
    
    write(*,*) forceHorzntPP,forceVrtcalPP,"force aft PP displcement"
    write(*,*) CgYNode(23),CgYNode(21),"CgY 23-21"
    
    CgXStr(1:N_node)  = CgXNode(1:N_node) ; CgYStr(1:N_node)  = CgYNode(1:N_node)
    CgXNode(1:N_node) = 0.0000d0          ; CgYNode(1:N_node) = 0.0000d0
    
    call get_gradient_with_changing_nodeTyp_test_withCompnts(node_xy,l0,A0,grdV2,grdSprV2,grdAreaV2,&
         grdGrvtnlV2,grdBendV2,grdSRypV2)
    
    force         = -grdV2
    forceHorzntPP = force(pulley_pnt,1)
    forceVrtcalPP = force(pulley_pnt,2)
    
    write(*,*) forceHorzntPP,forceVrtcalPP,"force aft PP displcement with Zero Force"
    write(*,*) CgYNode(23),CgYNode(21),"CgY 23-21"
    
    !call cmpreTwoGrdCompnts(grdV1,grdSprV1,grdAreaV1,grdGrvtnlV1,grdBendV1,&
    !     grdV2,grdSprV2,grdAreaV2,grdGrvtnlV2,grdBendV2)
         
  end subroutine displace_pulley_virtually_to_check_force
  
  
  subroutine turnOffbndryForceToChk
    implicit none
    integer :: BL1,BL2,BR1,BR2
    real*8  :: CgX_BL1,CgX_BL2,CgX_BR1,CgX_BR2
    
    open(unit=742,file='turnOffbndryForceToChk.dat',position='append')
    
    BL1 = 1 ; BR1 = (Hlf_Ncell+1)*2 + 1
    BL2 = 2 ; BR2 = (Hlf_Ncell+1)*2 + 2
    
    write(742,*) BL1,BL2,BR1,BR2,"BRs"
    write(742,*) Hlf_Ncell,"Hlf_Ncell"
    
    CgX_BL1 = CgXNode(BL1) ; CgX_BR1 = CgXNode(BR1) 
    CgX_BL2 = CgXNode(BL2) ; CgX_BR2 = CgXNode(BR2)
    
    write(742,*) CgX_BL1,CgX_BL2,CgX_BR1,CgX_BR2,"CgX of boundary"
    
    CgXNode(BL1) = 0.00d0 ; CgXNode(BR1) = 0.00d0
    CgXNode(BL2) = 0.00d0 ; CgXNode(BR2) = 0.00d0
    
    call Equilibrate_only_NI_model ; write(742,*) Frame_NI-1,"Frm - NI (1)"
    
    CgXNode(BL1) = CgX_BL1 ; CgXNode(BR1) = CgX_BR1
    CgXNode(BL2) = CgX_BL2 ; CgXNode(BR2) = CgX_BR2
    
    call Equilibrate_only_NI_model ; write(742,*) Frame_NI-1,"Frm - NI (2)"
    
    close(742)
    
  end subroutine turnOffbndryForceToChk
  
  
  subroutine make_the_lateralTnsn_IC_equaltoDesirdval()
    implicit none
    character(len=200)    :: flplce1,flplce2,flplce
    integer               :: ExpNo,FrmNoToBeRead,CgYCase
    real*8                :: TnsnB,TnsnA,TnsnNeeded,l0B,l0A,l0Needed
    integer               :: nsprsInACellT1,nsprsInACellT2,i,j
    integer               :: LtrlSprInitCellLft,LtrlSprInitCellRght
    
    CgYCase=3
    
    if (CgYCase==3) then 
       
       write(flplce1,*)"/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure"
       write(flplce2,*)"/Force_LtrlSprNC1_comboTst_CgY=0.75/ForceDownLtrlShortnALWAYS_TiltReduced/"
       flplce=trim(adjustl(flplce1))//trim(adjustl(flplce2))
       
       ExpNo = Exprmnt_NIVF ; FrmNoToBeRead = 29
       call read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
       Frame_NIVF = 1 ; call Equilibrate_only_NI_model_withVF_region
       
       TnsnNeeded = 0.18369396898325441
       TnsnB      = 0.18357391956781655 ; TnsnA = 0.26663969878163840
       l0B        = 2.250000000000000d0 ; l0A   = 2.125000000000000d0 ; l0Needed = 1.0d30
       
    elseif (CgYCase==2) then 
       
       write(flplce1,*)"/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure"
       write(flplce2,*)"/Force_LtrlSprNC1_comboTst_CgY=0.50/ForceDownLtrlShortnALWAYS_TiltReduced/"
       flplce=trim(adjustl(flplce1))//trim(adjustl(flplce2))
       
       ExpNo = Exprmnt_NIVF ; FrmNoToBeRead = 18
       call read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
       Frame_NIVF = 1 ; call Equilibrate_only_NI_model_withVF_region
       
       TnsnNeeded = 0.18369396898325441
       TnsnB      = 0.15039723392083815 ; TnsnA = 0.22709480657922221 
       l0B        = 2.375000000000000d0 ; l0A   = 2.250000000000000d0 ; l0Needed = 1.0d30
       
    elseif (CgYCase==1) then 
       
       write(flplce1,*)"/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure"
       write(flplce2,*)"/Force_LtrlSprNC1_comboTst_CgY=0.25/ForceDownLtrlShortnALWAYS_TiltReduced/"
       flplce=trim(adjustl(flplce1))//trim(adjustl(flplce2))
       
       ExpNo = Exprmnt_NIVF ; FrmNoToBeRead = 8
       call read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
       Frame_NIVF = 1 ; call Equilibrate_only_NI_model_withVF_region
       
       TnsnNeeded = 0.18369396898325441
       TnsnB      = 0.13234509312852083 ; TnsnA = 0.20331566027179448
       l0B        = 2.500000000000000d0 ; l0A   = 2.375000000000000d0 ; l0Needed = 1.0d30
       
    endif
    
    call get_the_l0val_to_achieveDesiredLtrlTnsn(TnsnB,TnsnA,TnsnNeeded,l0B,l0A,l0Needed)
    
    nsprsInACellT1 = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    nsprsInACellT2 = (NAEC_ApclCrtcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,(NAEC_Ltrl+1)
       
       LtrlSprInitCellLft  = (Hlf_Ncell-NCP_CrtclApSrfc)*(nsprsInACellT1) + &
            (NCP_CrtclApSrfc-1)*(nsprsInACellT2) + (NAEC_ApclCrtcl+1+NAEC_Bsal+1) + i 
       LtrlSprInitCellRght = LtrlSprInitCellLft +  (Hlf_Ncell-NCP_CrtclApSrfc)*(nsprsInACellT1) + &
            (NCP_CrtclApSrfc)*(nsprsInACellT2)
       write(*,*) LtrlSprInitCellLft,LtrlSprInitCellRght,i,"lft-rght"
       
       l0(LtrlSprInitCellLft)  = l0Needed
       l0(LtrlSprInitCellRght) = l0(LtrlSprInitCellLft)
       
    enddo
    
    call Equilibrate_only_NI_model_withVF_region
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
  end subroutine make_the_lateralTnsn_IC_equaltoDesirdval
  
  
  subroutine get_the_l0val_to_achieveDesiredLtrlTnsn(TnsnB,TnsnA,TnsnNeeded,l0B,l0A,l0Needed)
    implicit none
    real*8, intent(in)  :: TnsnB,TnsnA,TnsnNeeded,l0B,l0A
    real*8, intent(out) :: l0Needed
    integer             :: caseVal
    real*8              :: fctr
    
    if ((TnsnNeeded.gt.TnsnB).and.(TnsnNeeded.lt.TnsnA)) then
       caseVal=1
    elseif ((TnsnNeeded.gt.TnsnA).and.(TnsnNeeded.lt.TnsnB)) then
       caseVal=2
    else
       caseVal=0
       write(*,*) "Caseval is =",caseVal,"which should not be the case"
       stop 'inconsistent case value' 
    endif
    
    fctr     = (TnsnNeeded-TnsnB)/(TnsnA-TnsnB)
    l0Needed = l0B + fctr*(l0A-l0B)
    
    write(*,*) fctr,"fctr"
    write(*,*) TnsnNeeded,TnsnA,TnsnB,"TnsnNeeded-TnsnA-TnsnB"
    write(*,*) l0Needed,l0A,l0B,"l0Needed-l0A-l0B"
    
  end subroutine get_the_l0val_to_achieveDesiredLtrlTnsn
  
  
  subroutine manage_buckld_LtrlSpr_based_onTnsn()!ROUTINE IS FOR PulleyPnt SIM
    implicit none
    integer :: LtrlSp1(1:Hlf_Ncell), LtrlSp2(1:Hlf_Ncell), LtrlSp3(1:N_cell)
    
    call get_LtrlSps(LtrlSp1,LtrlSp2,LtrlSp3)
    call manplt_Tnsn_for_buckld_LtrlSprs(LtrlSp1,LtrlSp2,LtrlSp3)
    
  end subroutine manage_buckld_LtrlSpr_based_onTnsn
  
  subroutine manplt_Tnsn_for_buckld_LtrlSprs(LtrlSp1,LtrlSp2,LtrlSp3)
    implicit none
    integer,intent(in)    :: LtrlSp1(1:Hlf_Ncell),LtrlSp2(1:Hlf_Ncell),LtrlSp3(1:Hlf_Ncell)
    real*8                :: LtrlSpTnsn(1:Hlf_Ncell)
    real*8                :: TnsnLim,ZeroV
    integer               :: i,j,countNonBuckld
    integer               :: sprL1,sprL2,sprL3,sprR1,sprR2,sprR3
    integer               :: nsprsInACell
    
    LtrlSpTnsn   = -1.0d30
    TnsnLim      = 0.0200d0 ! TnsnLim
    ZeroV        = 0.0000d0 
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    call calculate_ltrlTnsn(LtrlSp1,LtrlSp2,LtrlSp3,LtrlSpTnsn)
    
    do 
       countNonBuckld = 0
       
       do i = 1,Hlf_Ncell
          
          if ((LtrlSpTnsn(i).gt.ZeroV) .or. (abs(LtrlSpTnsn(i)).lt.TnsnLim)) then
             ! both if tension is postive or tension is negative but very low value
             
             sprL1 = LtrlSp1(i) ; sprR1 = sprL1 + (nsprsInACell)*(Hlf_Ncell)
             sprL2 = LtrlSp2(i) ; sprR2 = sprL2 + (nsprsInACell)*(Hlf_Ncell)
             sprL3 = LtrlSp3(i) ; sprR3 = sprL3 + (nsprsInACell)*(Hlf_Ncell)
             
             l0(sprL1) = 0.9950d0*l0(sprL1) ; l0(sprR1) = l0(sprL1)
             l0(sprL2) = l0(sprL1)          ; l0(sprR2) = l0(sprL2)
             l0(sprL3) = l0(sprL1)          ; l0(sprR3) = l0(sprL3)
             
          else
             countNonBuckld = countNonBuckld+1 
          endif
       enddo
       
       call Equilibrate_only_NI_model
       call calculate_ltrlTnsn(LtrlSp1,LtrlSp2,LtrlSp3,LtrlSpTnsn)
       
       write(*,*) countNonBuckld,"countNonBuckld"
       if (countNonBuckld==Hlf_Ncell) exit
       
    enddo
    
  end subroutine manplt_Tnsn_for_buckld_LtrlSprs
  
  subroutine calculate_ltrlTnsn(LtrlSp1,LtrlSp2,LtrlSp3,LtrlSpTnsn)
    implicit none
    integer,intent(in)    :: LtrlSp1(1:Hlf_Ncell),LtrlSp2(1:Hlf_Ncell),LtrlSp3(1:Hlf_Ncell)
    real*8, intent(inout) :: LtrlSpTnsn(1:Hlf_Ncell)
    real*8                :: TnsnV(1:N_spr) 
    integer               :: i,j
    
    TnsnV(1:N_spr) = -1.0d30
    TnsnV(1:N_spr) = k_spr(1:N_spr)*(l0(1:N_spr)-l(1:N_spr))
    
    do i = 1,Hlf_Ncell
       LtrlSpTnsn(i) = (TnsnV(LtrlSp1(i))+TnsnV(LtrlSp2(i))+TnsnV(LtrlSp3(i)))/real(NAEC_Ltrl+1)
       write(*,*) LtrlSp1(i),LtrlSp2(i),LtrlSp3(i),LtrlSpTnsn(i),i,"LtrlSp"
    enddo
    
  end subroutine calculate_ltrlTnsn
  
  
  subroutine get_LtrlSps(LtrlSp1,LtrlSp2,LtrlSp3)
    implicit none
    integer, intent(out) :: LtrlSp1(1:Hlf_Ncell),LtrlSp2(1:Hlf_Ncell),LtrlSp3(1:Hlf_Ncell)
    integer              :: i,j,nsprsInACell

    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,Hlf_Ncell
       LtrlSp1(i)  = (i-1)*(nsprsInACell)+(NAEC_Apcl+1)+(NAEC_Bsal+1)+1
       LtrlSp2(i)  = LtrlSp1(i)+1
       LtrlSp3(i)  = LtrlSp1(i)+2
       write(*,*) LtrlSp1(i),LtrlSp2(i),LtrlSp3(i),i,"LtrlSp"
    enddo
    
  end subroutine get_LtrlSps
  
  
end module MAE_initaitn
