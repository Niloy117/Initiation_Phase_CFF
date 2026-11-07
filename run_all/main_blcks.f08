program main
  use Complt_Simltn
  
  implicit none
  
  stageNo   = 1 !1=1E , 2=1L, 3=2E, 4=2L
  stageType = 1 !1 for odd no of N_cell(symm), 2 for even no of N_cell(antisymm)
  
  call get_stage_info 
  call initialize_system_params
  
  call get_sys_lgcls(.True.,.True.,.True.,.True.,.False.)
  call get_sys_lgcls_cnstr(.False.,.True.,.True.,.False.,.False.)
  
  call get_l0A0_varitn_lgcls(.False.,.False.)
  
  call input_from_bash_about_grdnt_choice !gradient_calc.f08
  
  EquilAlgrthm = 1
  
  if (stageNo==1 .and. stageType==1) call stage1_type1
  if (stageNo==1 .and. stageType==2) call stage1_type2
  if (stageNo==2 .and. stageType==1) call stage2_type1
  if (stageNo==2 .and. stageType==2) call stage2_type2
  if (stageNo==3 .and. stageType==1) call stage3_type1
  if (stageNo==3 .and. stageType==2) call stage3_type2
  if (stageNo==4 .and. stageType==1) call stage4_type1
  if (stageNo==4 .and. stageType==2) call stage4_type2
  
end program main
