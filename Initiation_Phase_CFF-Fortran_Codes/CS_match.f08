module Control_State_Length_Match
  use MAE_initaitn
  
  implicit none
  integer :: NCMI  = 4 ! NCMI=N_cellsMeet_InitPhase
  integer :: NumCS = 7 ! Num of CntlStates To Work[should be 8,but we arent working with the CM0ICS]
  real*8  :: ksIncrFctr
  
  integer :: find_allorNot = 0         ! wont be needed when all Exprmnt CS is measured
  real*8  :: TOL_fctr      = 0.9800d0  
  real*8  :: yvalBotmLayr
  real*8  :: BndryPress
  real*8  :: TolPressDiffToMaintnOrientn
  real*8  :: TOL_PDwrtBP    = 1.5000d0    !Tolerance Pressure Difference wrt Boundary Pressure
  
  real*8, allocatable  :: Exprmntl_ApLn_PCS(:,:), Exprmntl_ApLn_ICS(:,:)
  real*8, allocatable  :: Exprmntl_BsLn_PCS(:,:), Exprmntl_BsLn_ICS(:,:)
  real*8, allocatable  :: Exprmntl_LtLn_PCS(:,:), Exprmntl_LtLn_ICS(:,:)
  
  real*8, allocatable  :: Simulted_ApLn_PCS(:,:), Simulted_ApLn_ICS(:,:)
  real*8, allocatable  :: Simulted_BsLn_PCS(:,:), Simulted_BsLn_ICS(:,:)
  real*8, allocatable  :: Simulted_LtLn_PCS(:,:), Simulted_LtLn_ICS(:,:)
  
  real*8, allocatable  :: fctrAR_BSwrtAP_PCS(:) !BSwrtAP = Basal side with_respect_to Apical side
  real*8, allocatable  :: fctrAR_LTwrtAP_PCS(:) !LTwrtAP = Lateral side with_respect_to Apical side
  real*8, allocatable  :: fctrAR_BSwrtAP_ICS(:) !AR=Aspect Ratio
  real*8, allocatable  :: fctrAR_LTwrtAP_ICS(:)
  
  real*8, allocatable  :: fctr_Apcl_PCS(:),fctr_Apcl_ICS(:)
  real*8, allocatable  :: fctr_Bsal_PCS(:),fctr_Bsal_ICS(:)
  real*8, allocatable  :: fctr_Ltrl_PCS(:),fctr_Ltrl_ICS(:)
  real*8, allocatable  :: Unit_Fctr_PCS(:),Unit_Fctr_ICS(:) 
  
  real*8, allocatable  :: Expected_ApLn_PCS(:,:), Expected_ApLn_ICS(:,:)
  real*8, allocatable  :: Expected_BsLn_PCS(:,:), Expected_BsLn_ICS(:,:)
  real*8, allocatable  :: Expected_LtLn_PCS(:,:), Expected_LtLn_ICS(:,:)
  
  real*8, allocatable  :: diff_A0andA(:),PressStr(:),TOL_Press(:)
  
  real*8               :: TolPrcnt1L=0.950d0,TolPrcnt2L=0.990d0
  real*8               :: TolPrcnt1G=1.050d0,TolPrcnt2G=1.010d0
  
  real*8,  allocatable :: Simulted_CHLtPCS(:,:),Simulted_CHBsPCS(:,:)
  real*8,  allocatable :: Exprmntl_CHLtPCS(:,:),Exprmntl_CHBsPCS(:,:)
  real*8,  allocatable :: Expected_CHLtPCS(:,:),Expected_CHBsPCS(:,:)
  integer, allocatable :: CurvtrSign_LtPCS(:,:),CurvtrSign_BsPCS(:,:) 
  
  real*8, allocatable  :: unitCurvCalcIC_BsPCS(:),unitCurvCalcIC_LtPCS(:)
  real*8, allocatable  :: ratioMeasuremnt(:)
  
  integer              :: writeForCurvtr=0
  real*8, allocatable  :: Exprmntl_PPbfr(:),Exprmntl_PPaft(:)
  real*8, allocatable  :: Expected_PPbfr(:),Expected_PPaft(:)
  
  real*8               :: disTnceBfrPull
  real*8               :: intended_TIIF=0.5000d0
  
contains
  
  subroutine parameter_values_for_module_CS_match
    implicit none
    
    BndryPress                  =  0.1000d0
    TolPressDiffToMaintnOrientn = (0.02500d0)*(BndryPress)
    
  end subroutine parameter_values_for_module_CS_match
  
  subroutine allocate_and_initialize_lenValues
    implicit none
    
    allocate(Exprmntl_ApLn_PCS(1:NCMI,1:Hlf_Ncell),    Exprmntl_ApLn_ICS(1:NCMI,1:Hlf_Ncell))
    allocate(Exprmntl_BsLn_PCS(1:NCMI,1:(Hlf_Ncell+1)),Exprmntl_BsLn_ICS(1:NCMI,1:(Hlf_Ncell+1)))
    allocate(Exprmntl_LtLn_PCS(1:NCMI,1:Hlf_Ncell),    Exprmntl_LtLn_ICS(1:NCMI,1:Hlf_Ncell))
    
    Exprmntl_ApLn_PCS=-1.0d30 ; Exprmntl_ApLn_ICS=-1.0d30
    Exprmntl_BsLn_PCS=-1.0d30 ; Exprmntl_BsLn_ICS=-1.0d30
    Exprmntl_LtLn_PCS=-1.0d30 ; Exprmntl_LtLn_ICS=-1.0d30
    
    allocate(Simulted_ApLn_PCS(1:NCMI,1:Hlf_Ncell),    Simulted_ApLn_ICS(1:NCMI,1:Hlf_Ncell))
    allocate(Simulted_BsLn_PCS(1:NCMI,1:(Hlf_Ncell+1)),Simulted_BsLn_ICS(1:NCMI,1:(Hlf_Ncell+1)))
    allocate(Simulted_LtLn_PCS(1:NCMI,1:Hlf_Ncell),    Simulted_LtLn_ICS(1:NCMI,1:Hlf_Ncell))
    
    Simulted_ApLn_PCS=-1.0d30 ; Simulted_ApLn_ICS=-1.0d30
    Simulted_BsLn_PCS=-1.0d30 ; Simulted_BsLn_ICS=-1.0d30
    Simulted_LtLn_PCS=-1.0d30 ; Simulted_LtLn_ICS=-1.0d30
    
    allocate(Expected_ApLn_PCS(1:NCMI,1:Hlf_Ncell),    Expected_ApLn_ICS(1:NCMI,1:Hlf_Ncell))
    allocate(Expected_BsLn_PCS(1:NCMI,1:(Hlf_Ncell+1)),Expected_BsLn_ICS(1:NCMI,1:(Hlf_Ncell+1)))
    allocate(Expected_LtLn_PCS(1:NCMI,1:Hlf_Ncell),    Expected_LtLn_ICS(1:NCMI,1:Hlf_Ncell))
    
    Expected_ApLn_PCS=-1.0d30 ; Expected_ApLn_ICS=-1.0d30
    Expected_BsLn_PCS=-1.0d30 ; Expected_BsLn_ICS=-1.0d30
    Expected_LtLn_PCS=-1.0d30 ; Expected_LtLn_ICS=-1.0d30
    
    allocate(diff_A0andA(1:Hlf_Ncell+1)) ; diff_A0andA = -1.0d30
    allocate(PressStr(1:Hlf_Ncell+1))    ; PressStr    = -1.0d30
    allocate(TOL_Press(1:(Hlf_Ncell+1))) ; TOL_Press   = -1.0d30
    
    allocate(Simulted_CHBsPCS(1:NCMI,1:Hlf_Ncell),Simulted_CHLtPCS(1:NCMI,1:Hlf_Ncell))
    allocate(Exprmntl_CHBsPCS(1:NCMI,1:Hlf_Ncell),Exprmntl_CHLtPCS(1:NCMI,1:Hlf_Ncell))
    allocate(Expected_CHBsPCS(1:NCMI,1:Hlf_Ncell),Expected_CHLtPCS(1:NCMI,1:Hlf_Ncell))
    allocate(CurvtrSign_BsPCS(1:NCMI,1:Hlf_Ncell),CurvtrSign_LtPCS(1:NCMI,1:Hlf_Ncell))
    
    Simulted_CHBsPCS=-1.0d30 ; Simulted_CHLtPCS=-1.0d30
    Exprmntl_CHBsPCS=-1.0d30 ; Exprmntl_CHLtPCS=-1.0d30
    Expected_CHBsPCS=-1.0d30 ; Expected_CHLtPCS=-1.0d30
    CurvtrSign_BsPCS=-10     ; CurvtrSign_LtPCS=-10
    
    allocate(unitCurvCalcIC_BsPCS(1:NCMI)) ; unitCurvCalcIC_BsPCS = -1.0d30
    allocate(unitCurvCalcIC_LtPCS(1:NCMI)) ; unitCurvCalcIC_LtPCS = -1.0d30
    allocate(ratioMeasuremnt(1:NCMI))      ; ratioMeasuremnt      = -1.0d30 
    
  end subroutine allocate_and_initialize_lenValues
  
  subroutine allocate_and_initialize_PPValues
    implicit none
    
    allocate(Exprmntl_PPbfr(1:NCMI),Exprmntl_PPaft(1:NCMI)) 
    allocate(Expected_PPbfr(1:NCMI),Expected_PPaft(1:NCMI)) 
    
    Exprmntl_PPbfr=-1.0d30 ; Exprmntl_PPbfr=-1.0d30
    Expected_PPbfr=-1.0d30 ; Expected_PPaft=-1.0d30
    
  end subroutine allocate_and_initialize_PPValues
  
  subroutine get_CL_frm_all_Sim_CS
    implicit none
    integer  :: CMv,PCSorICS
    integer  :: cntLp,cntLpMax
    integer  :: lpInfo
    integer  :: i,j,jmax
    
    cntLp = 1 ; cntLpMax=7
    
    do
       if (cntLp==1) lpInfo=11
       if (cntLp==2) lpInfo=22
       if (cntLp==3) lpInfo=21
       if (cntLp==4) lpInfo=32
       if (cntLp==5) lpInfo=31
       if (cntLp==6) lpInfo=42
       if (cntLp==7) lpInfo=41
       
       CMv = lpInfo/10 ; PCSorICS = mod(lpInfo,10)
       write(*,*) CMv,PCSorICS,"CMv,PCSorICS"
       
       call get_CL_frm_Sim_CS(CMv,PCSorICS)
       
       if (cntLp.ge.cntLpMax) exit
       cntLp = cntLp+1
       
    enddo
    
    ! call print_Simlted_Lengths
    
  end subroutine get_CL_frm_all_Sim_CS
  
  
  subroutine get_CL_frm_Sim_CS(CMv,PCSorICS) ! CL=Current_length, CS=Control_state
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    character(len=100)  :: sprFlnm,sprFlnmPrfx
    character(len=100)  :: convSprFlnm,convSprFlnmPrfx
    character(len=100)  :: PCSorICS_txt
    
    integer :: i,j,jmax
    integer :: sprNm,cellNm
    real*8  :: Ta,Tb,Tl,ksA,ksB,ksL,lA,lB,lL,l0A,l0B,l0L
    real*8  :: TnV,ksV,lV,l0V
    integer :: nsprsInACell,symmSpr
    
    sprFlnmPrfx='sprsPrps_CM' ; convSprFlnmPrfx='convSprs_CM'
    
    if (PCSorICS==1) PCSorICS_txt='PCS.dat'
    if (PCSorICS==2) PCSorICS_txt='ICS.dat'
    
    write(sprFlnm,'(a,i1.1,a)')trim(adjustl(sprFlnmPrfx)),CMv,trim(adjustl(PCSorICS_txt)) 
    write(convSprFlnm,'(a,i1.1,a)')trim(adjustl(convSprFlnmPrfx)),CMv,trim(adjustl(PCSorICS_txt))
    
    write(*,*) trim(adjustl(sprFlnm))
    write(*,*) trim(adjustl(convSprFlnm))
    
    open(unit=677,file=trim(adjustl(sprFlnm)))
    open(unit=678,file=trim(adjustl(convSprFlnm)))
    
    do i = 1,Hlf_Ncell
       read(678,*) cellNm,Ta,Tb,Tl,ksA,ksB,ksL,lA,lB,lL,l0A,l0B,l0L
       
       if (PCSorICS==1) then   
          Simulted_ApLn_PCS(CMv,i) = lA
          Simulted_BsLn_PCS(CMv,i) = lB
          Simulted_LtLn_PCS(CMv,i) = lL
       elseif (PCSorICS==2) then
          Simulted_ApLn_ICS(CMv,i) = lA
          Simulted_BsLn_ICS(CMv,i) = lB
          Simulted_LtLn_ICS(CMv,i) = lL
       endif
       
    enddo
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    symmSpr      = 2*(Hlf_Ncell)*(nsprsInACell) ; write(*,*) symmSpr,"symmSpr"

    if (PCSorICS==1) Simulted_BsLn_PCS(CMv,Hlf_Ncell+1) = 0.0d0
    if (PCSorICS==2) Simulted_BsLn_ICS(CMv,Hlf_Ncell+1) = 0.0d0
    
    do i = 1,N_spr
       read(677,*) TnV,ksV,lV,l0V,sprNm
       
       if (i.gt.symmSpr) then
          
          if(PCSorICS==1) Simulted_BsLn_PCS(CMv,Hlf_Ncell+1) = Simulted_BsLn_PCS(CMv,Hlf_Ncell+1)+lV
          if(PCSorICS==2) Simulted_BsLn_ICS(CMv,Hlf_Ncell+1) = Simulted_BsLn_ICS(CMv,Hlf_Ncell+1)+lV
          
          if (CMv==2.and.PCSorICS==1) write(*,*) TnV,ksV,lV,l0V,sprNm,"Tn-ks-l-lV-spr"
          if (CMv==3.and.PCSorICS==1) write(*,*) TnV,ksV,lV,l0V,sprNm,"Tn-ks-l-lV-spr"
          if (CMv==4.and.PCSorICS==1) write(*,*) TnV,ksV,lV,l0V,sprNm,"Tn-ks-l-lV-spr"
          
       endif
       
    enddo
    
    close(677)
    close(678)
    
  end subroutine get_CL_frm_Sim_CS
  
  
  subroutine print_Simlted_Lengths
    implicit none
    integer :: i,j,jmax
    
    do i = 1,4

       jmax = Hlf_Ncell
       
       do j = 1,jmax
          
          write(*,*) Simulted_ApLn_PCS(i,j),Simulted_ApLn_ICS(i,j),i,j,"Sim-Ap"
          write(*,*) Simulted_BsLn_PCS(i,j),Simulted_BsLn_ICS(i,j),i,j,"Sim-Bs"
          write(*,*) Simulted_LtLn_PCS(i,j),Simulted_LtLn_ICS(i,j),i,j,"Sim-Lt" 
          
          if (j==jmax) write(*,*) Simulted_BsLn_PCS(i,j+1),Simulted_BsLn_ICS(i,j+1),i,j+1,"Sim-Bs"
          
       enddo
       
    enddo
    
  end subroutine print_Simlted_Lengths
  
  subroutine get_CL_frm_Exp_CS()
    implicit none
    integer :: CMv
    real*8  :: ratioMeasuremnt_val=-1.0d30
    integer :: i,j,imax,jmax
    
    CMv = 1 ! (2G9b-m)
    
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val)
    ratioMeasuremnt(CMv) = ratioMeasuremnt_val ; write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt insd"
    
    
    Exprmntl_ApLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,1)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,2)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,3)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,4)  = 0.4300d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,5)  = 0.3400d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,6)  = 0.3200d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,7)  = 0.2100d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,8)  = 0.3800d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,9)  = 0.3100d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,10) = 0.4600d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,11) = 0.4900d0 * (ratioMeasuremnt(CMv))!27+22
    
    Exprmntl_BsLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,1)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,2)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,3)  = 0.3900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,4)  = 0.4000d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,5)  = 0.3600d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,6)  = 0.2900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,7)  = 0.2100d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,8)  = 0.3000d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,9)  = 0.3200d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,10) = 0.3000d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,11) = 0.2800d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,12) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,12) = 0.3700d0 * (ratioMeasuremnt(CMv))
    
    Exprmntl_LtLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,1)  = 1.8600d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,2)  = 1.8600d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,3)  = 1.8600d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,4)  = 1.8500d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,5)  = 1.8800d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,6)  = 1.8800d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,7)  = 1.8800d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,8)  = 1.9200d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,9)  = 1.8900d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,10) = 1.9300d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,11) = 1.7300d0 * (ratioMeasuremnt(CMv))
    
    ! the image that I had used here 2G9b-m, in that image [Exprmntl_LtLnPCS(CMv,1:2) = 1.76]; however I used
    ! the value as 1.86 to match same length.
    
    write(*,*) " "
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_ApLn_PCS(CMv,i),Exprmntl_BsLn_PCS(CMv,i),Exprmntl_LtLn_PCS(CMv,i),i,"Exprmnt CMv1"
    enddo
    write(*,*) Exprmntl_BsLn_PCS(CMv,Hlf_Ncell+1),"Exprmnt CMv1"
    
    
    CMv = 2 ! (3G6b-m)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !   ONE IMPORTANT THING TO NOTICE THAT THE FIRST THREE CELLS EXPERIMENTAL LENGTH
    !   IS NOT MODELLED BY COPYING FROM THE 4TH CELL WHICH IS THE LAST MEASURED LENGTH
    !   INSTEAD THE FIRST 3 CELLS ARE MODELLED AS FAR CELLS WHERE THE APICAL AND
    !   BASAL MEMBRANE LENGTHS ARE EQUAL AND THE LATERAL MEMBRANE LENGTH IS MEASURED
    !   AS THE VERICAL DISTANCE  OF APICAL MEMBRANE LAYER (FIXED HERE = 0.00D0 AND
    !   BASAL MEMBRANE LAYER (MEASURED FROM
    !   4TH CELL'S ACTUAL (POSTERIOR DIRECTION) LATERAL MEMBRANE'S Y-NODE VALUE)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !   FOR EXAMPLE UNIT DISTANCE = (4TH CELL'S SIM_APLN / 4TH CELL'S EXP APLN)
    !                               [HERE I USED THE 4TH CELL TO UNIT MEASUREMENT]
    !                             = (1.764980547982183 / 0.48) = 3.6770428082962145
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !   FOR CELL (1-3), APICAL AND BASAL ARE CONSIDERED SAME
    !   AND LATERAL LENGTH EXPERIMNTL VALUE = (ABSOLUTE_VAL(Y9-Y10)) / UNIT DISTANCE
    !                                       = (0-(-8.3382621590942083)) / 3.6770428082962145
    !                                       = 2.2676543608035407
    
    !   DAT FILE USED NI_modelInitiatnNW036.dat,
    !   AND FILE_STORE_LOCTN: ~/ restructuring_Data_Structure/NI_DataInitPhase/HighKTest/WT_AR_Fix
    
    
    Exprmntl_ApLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,1)  = 0.480d0
    Exprmntl_ApLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,2)  = 0.480d0
    Exprmntl_ApLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,3)  = 0.480d0
    Exprmntl_ApLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,4)  = 0.480d0
    Exprmntl_ApLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,5)  = 0.460d0
    Exprmntl_ApLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,6)  = 0.480d0
    Exprmntl_ApLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,7)  = 0.350d0
    Exprmntl_ApLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,8)  = 0.440d0
    Exprmntl_ApLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,9)  = 0.430d0
    Exprmntl_ApLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,10) = 0.780d0 !0.49 + 0.29
    Exprmntl_ApLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,11) = 0.470d0
    
    Exprmntl_BsLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,1)  = 0.480d0
    Exprmntl_BsLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,2)  = 0.480d0
    Exprmntl_BsLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,3)  = 0.480d0
    Exprmntl_BsLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,4)  = 0.310d0
    Exprmntl_BsLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,5)  = 0.350d0
    Exprmntl_BsLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,6)  = 0.300d0
    Exprmntl_BsLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,7)  = 0.470d0
    Exprmntl_BsLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,8)  = 0.200d0
    Exprmntl_BsLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,9)  = 0.300d0
    Exprmntl_BsLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,10) = 0.310d0
    Exprmntl_BsLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,11) = 0.370d0
    Exprmntl_BsLn_ICS(CMv,12) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,12) = 0.700d0
    
    Exprmntl_LtLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,1)  = 2.2676543608035407d0
    Exprmntl_LtLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,2)  = 2.2676543608035407d0
    Exprmntl_LtLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,3)  = 2.2676543608035407d0
    Exprmntl_LtLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,4)  = 2.310d0
    Exprmntl_LtLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,5)  = 2.310d0
    Exprmntl_LtLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,6)  = 2.320d0
    Exprmntl_LtLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,7)  = 2.310d0
    Exprmntl_LtLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,8)  = 2.330d0
    Exprmntl_LtLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,9)  = 2.360d0
    Exprmntl_LtLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,10) = 2.210d0
    Exprmntl_LtLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,11) = 1.730d0
    
    
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!
    
    CMv = 3 ! (2G8b-m)
    
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val)
    ratioMeasuremnt(CMv) = ratioMeasuremnt_val ; write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt insd"
    
    Exprmntl_ApLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,1)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,2)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,3)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,4)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,5)  = 0.520d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,6)  = 0.550d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,7)  = 0.250d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,8)  = 0.380d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,9)  = 0.760d0 * (ratioMeasuremnt(CMv))!34+42
    Exprmntl_ApLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,10) = 0.670d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,11) = 0.340d0 * (ratioMeasuremnt(CMv))
    
    Exprmntl_BsLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,1)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,2)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,3)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,4)  = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,5)  = 0.240d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,6)  = 0.280d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,7)  = 0.300d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,8)  = 0.150d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,9)  = 0.230d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,10) = 0.120d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,11) = 0.400d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,12) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,12) = 0.850d0 * (ratioMeasuremnt(CMv))
    
    Exprmntl_LtLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,1)  = 2.050d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,2)  = 2.050d0 * (ratioMeasuremnt(CMv)) 
    Exprmntl_LtLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,3)  = 2.050d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,4)  = 2.050d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,5)  = 2.070d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,6)  = 2.130d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,7)  = 2.060d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,8)  = 2.180d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,9)  = 2.010d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,10) = 1.490d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,11) = 1.190d0 * (ratioMeasuremnt(CMv))

    write(*,*) " "
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_ApLn_PCS(CMv,i),Exprmntl_BsLn_PCS(CMv,i),Exprmntl_LtLn_PCS(CMv,i),i,"Exprmnt CMv3"
    enddo
    write(*,*) Exprmntl_BsLn_PCS(CMv,Hlf_Ncell+1),"Exprmnt CMv3"
    
    
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!
    !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!
    
    CMv = 4 ! (7G12b-m)
    
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val)
    ratioMeasuremnt(CMv) = ratioMeasuremnt_val ; write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt insd"
    
    Exprmntl_ApLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,1)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,2)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,3)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,4)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,5)  = 0.370d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,6)  = 0.470d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,7)  = 0.640d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,8)  = 1.150d0 * (ratioMeasuremnt(CMv)) ! 0.94+0.21
    Exprmntl_ApLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,9)  = 0.560d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,10) = 0.470d0 * (ratioMeasuremnt(CMv))
    Exprmntl_ApLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_ApLn_PCS(CMv,11) = 0.130d0 * (ratioMeasuremnt(CMv))
    
    Exprmntl_BsLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,1)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,2)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,3)  = 0.340d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,4)  = 0.320d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,5)  = 0.370d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,6)  = 0.400d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,7)  = 0.370d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,8)  = 0.230d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,9)  = 0.250d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,10) = 0.360d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,11) = 0.410d0 * (ratioMeasuremnt(CMv))
    Exprmntl_BsLn_ICS(CMv,12) = -1.0d30 ;   Exprmntl_BsLn_PCS(CMv,12) = 0.700d0 * (ratioMeasuremnt(CMv))
    
    Exprmntl_LtLn_ICS(CMv,1)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,1)  = 2.050d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,2)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,2)  = 2.050d0 * (ratioMeasuremnt(CMv)) 
    Exprmntl_LtLn_ICS(CMv,3)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,3)  = 2.050d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,4)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,4)  = 2.070d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,5)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,5)  = 2.060d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,6)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,6)  = 2.030d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,7)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,7)  = 2.150d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,8)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,8)  = 2.430d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,9)  = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,9)  = 1.930d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,10) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,10) = 1.640d0 * (ratioMeasuremnt(CMv))
    Exprmntl_LtLn_ICS(CMv,11) = -1.0d30 ;   Exprmntl_LtLn_PCS(CMv,11) = 1.300d0 * (ratioMeasuremnt(CMv))
    
    write(*,*) " "
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_ApLn_PCS(CMv,i),Exprmntl_BsLn_PCS(CMv,i),Exprmntl_LtLn_PCS(CMv,i),i,"Exprmnt CMv4"
    enddo
    write(*,*) Exprmntl_BsLn_PCS(CMv,Hlf_Ncell+1),"Exprmnt CMv4"
    
    
  end subroutine get_CL_frm_Exp_CS
  
  
  subroutine get_EL_frm_taking_the_4thCell_asBase
    implicit none
    integer :: i,j,jmax
    
    call get_fctr_frm_taking_the_4thCell_asBase!4th as we've 8 cell in Exprmntl,but in Simln 11
    
    do i = 1,NCMI
       
       if (i==1 .or. i==2 .or. i==3 .or. i==4) then
          
          jmax = Hlf_Ncell
          
          do j = 1,jmax
             
             Expected_ApLn_PCS(i,j) = Unit_Fctr_PCS(i) * Exprmntl_ApLn_PCS(i,j) 
             Expected_BsLn_PCS(i,j) = Unit_Fctr_PCS(i) * Exprmntl_BsLn_PCS(i,j)
             Expected_LtLn_PCS(i,j) = Unit_Fctr_PCS(i) * Exprmntl_LtLn_PCS(i,j)
             
             if (j==jmax) Expected_BsLn_PCS(i,j+1) = Unit_Fctr_PCS(i) * Exprmntl_BsLn_PCS(i,j+1)
             
             write(*,*) " "
             write(*,*) Expected_ApLn_PCS(i,j),fctr_Apcl_PCS(i),Exprmntl_ApLn_PCS(i,j),i,j,"Ap-PCS-i-j"
             write(*,*) Expected_BsLn_PCS(i,j),fctr_Bsal_PCS(i),Exprmntl_BsLn_PCS(i,j),i,j,"Bs-PCS-i-j"
             write(*,*) Expected_LtLn_PCS(i,j),fctr_Ltrl_PCS(i),Exprmntl_LtLn_PCS(i,j),i,j,"Lt-PCS-i-j"
             
             if (j==jmax) then
                write(*,*) " "
                write(*,*) Expected_BsLn_PCS(i,j+1),fctr_Bsal_PCS(i),Exprmntl_BsLn_PCS(i,j+1),i,j+1,"Bs-PCS-i-j+1"
             endif
             
          enddo
          
       endif
             
    enddo
    
  end subroutine get_EL_frm_taking_the_4thCell_asBase
  
  
  subroutine adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val)
    implicit none
    integer, intent(in)  :: CMv
    real*8 , intent(out) :: ratioMeasuremnt_val 
    real*8               :: ExpLtLnforPCS_CMv2
    real*8               :: MeasuredExpLtLnforPCS_CMv1,MeasuredExpLtLnforPCS_CMv3,MeasuredExpLtLnforPCS_CMv4
    
    if (CMv == 1) then
       
       ExpLtLnforPCS_CMv2         = 2.2676543608035407d0
       MeasuredExpLtLnforPCS_CMv1 = 1.8600d0
       ratioMeasuremnt_val        = ExpLtLnforPCS_CMv2/MeasuredExpLtLnforPCS_CMv1
       write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt_val"
       
    elseif (CMv == 2) then
       
       ratioMeasuremnt_val        = 1.0000d0   
       write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt_val"
       
    elseif (CMv == 3) then
       
       ExpLtLnforPCS_CMv2         = 2.2676543608035407d0
       MeasuredExpLtLnforPCS_CMv3 = 2.0500d0
       ratioMeasuremnt_val        = ExpLtLnforPCS_CMv2/MeasuredExpLtLnforPCS_CMv3
       write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt_val"
       
    elseif (CMv==4) then
       
       ExpLtLnforPCS_CMv2         = 2.2676543608035407d0
       MeasuredExpLtLnforPCS_CMv4 = 2.1100d0
       ratioMeasuremnt_val        = ExpLtLnforPCS_CMv2/MeasuredExpLtLnforPCS_CMv4
       write(*,*) ratioMeasuremnt_val,"ratioMeasuremnt_val"
       
    endif
    
  end subroutine adjustmnt_Exprmntl_len_measuremnt
  
  subroutine allocate_and_initialize_fctrValues
    implicit none
    
    allocate(fctrAR_BSwrtAP_PCS(1:NCMI),fctrAR_LTwrtAP_PCS(1:NCMI))
    allocate(fctrAR_BSwrtAP_ICS(1:NCMI),fctrAR_LTwrtAP_ICS(1:NCMI))
    
    fctrAR_BSwrtAP_PCS = -1.00d0 ; fctrAR_LTwrtAP_PCS = -1.00d0
    fctrAR_BSwrtAP_ICS = -1.00d0 ; fctrAR_LTwrtAP_ICS = -1.00d0
    
    allocate(fctr_Apcl_ICS(1:NCMI)) ; allocate(fctr_Apcl_PCS(1:NCMI))
    allocate(fctr_Bsal_ICS(1:NCMI)) ; allocate(fctr_Bsal_PCS(1:NCMI))
    allocate(fctr_Ltrl_ICS(1:NCMI)) ; allocate(fctr_Ltrl_PCS(1:NCMI))
    
    fctr_Apcl_ICS = -1.0d0 ; fctr_Apcl_PCS = -1.0d0
    fctr_Bsal_ICS = -1.0d0 ; fctr_Bsal_PCS = -1.0d0
    fctr_Ltrl_ICS = -1.0d0 ; fctr_Ltrl_PCS = -1.0d0
    
    allocate(Unit_Fctr_PCS(1:NCMI),Unit_Fctr_ICS(1:NCMI))
    Unit_Fctr_PCS=-1.00d0 ; Unit_Fctr_ICS=-1.00d0
    
  end subroutine allocate_and_initialize_fctrValues
  
  subroutine get_fctr_frm_taking_the_4thCell_asBase 
    implicit none
    integer :: i
    
    call allocate_and_initialize_fctrValues
    
    do i = 1,NCMI
       
       if (find_allorNot==1) then
          
          fctr_Apcl_ICS(i) = Simulted_ApLn_ICS(i,4)/Exprmntl_ApLn_ICS(i,4)
          fctr_Bsal_ICS(i) = Simulted_BsLn_ICS(i,4)/Exprmntl_BsLn_ICS(i,4)
          fctr_Ltrl_ICS(i) = Simulted_LtLn_ICS(i,4)/Exprmntl_LtLn_ICS(i,4)
          
          fctrAR_BSwrtAP_ICS(i) = Exprmntl_BsLn_ICS(i,4)/Exprmntl_ApLn_ICS(i,4)
          fctrAR_LTwrtAP_ICS(i) = Exprmntl_LtLn_ICS(i,4)/Exprmntl_ApLn_ICS(i,4)
          
          Unit_Fctr_ICS(i) = Simulted_ApLn_ICS(i,4)/Exprmntl_ApLn_ICS(i,4)
          
          fctr_Apcl_PCS(i) = Simulted_ApLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
          fctr_Bsal_PCS(i) = Simulted_BsLn_PCS(i,4)/Exprmntl_BsLn_PCS(i,4)
          fctr_Ltrl_PCS(i) = Simulted_LtLn_PCS(i,4)/Exprmntl_LtLn_PCS(i,4)
          
          fctrAR_BSwrtAP_PCS(i) = Exprmntl_BsLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
          fctrAR_LTwrtAP_PCS(i) = Exprmntl_LtLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
          
          Unit_Fctr_PCS(i) = Simulted_ApLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
          
       elseif (find_allorNot==0) then
          
          if (i==1 .or. i==2 .or. i==3 .or. i==4) then
                              
             fctr_Apcl_PCS(i) = Simulted_ApLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
             fctr_Bsal_PCS(i) = Simulted_BsLn_PCS(i,4)/Exprmntl_BsLn_PCS(i,4)
             fctr_Ltrl_PCS(i) = Simulted_LtLn_PCS(i,4)/Exprmntl_LtLn_PCS(i,4)
             
             write(*,*) fctr_Apcl_PCS(i),fctr_Bsal_PCS(i),fctr_Ltrl_PCS(i),i,"fctr_Apcl-Bsal-Ltrl"
             
             fctrAR_BSwrtAP_PCS(i) = Exprmntl_BsLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
             fctrAR_LTwrtAP_PCS(i) = Exprmntl_LtLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
             
             write(*,*) fctrAR_BSwrtAP_PCS(i),fctrAR_LTwrtAP_PCS(i),i,"fctr AR"
             
             if (i==1) Unit_Fctr_PCS(i) = 3.5839837092863056d0
             if (i==2) Unit_Fctr_PCS(i) = Simulted_ApLn_PCS(i,4)/Exprmntl_ApLn_PCS(i,4)
             if (i==3) Unit_Fctr_PCS(i) = Unit_Fctr_PCS(2)
             if (i==4) Unit_Fctr_PCS(i) = Unit_Fctr_PCS(2)
             
             write(*,*) Unit_Fctr_PCS(i),Simulted_ApLn_PCS(i,4),Exprmntl_ApLn_PCS(i,4),"USE"
             
          endif
          
       endif
       
       write(*,*) fctr_Apcl_ICS(i),fctr_Apcl_PCS(i),i,"fctr_Ap_ICS/PCS"
       write(*,*) fctr_Bsal_ICS(i),fctr_Bsal_PCS(i),i,"fctr_Bs_ICS/PCS"
       write(*,*) fctr_Ltrl_ICS(i),fctr_Ltrl_PCS(i),i,"fctr_Lt_ICS/PCS"
       write(*,*) " "
       
    enddo
    
  end subroutine get_fctr_frm_taking_the_4thCell_asBase
  
  
  subroutine save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    integer :: i,j,imax,jmax
    integer :: sprApTN,sprBsTN,sprLtTN
    integer :: nspInACellTN,nspInACellNI
    integer :: sprT,NumOFSgmnt
    
    integer :: sprApNI(1:(NAEC_Apcl+1)),sprBsNI(1:(NAEC_Bsal+1)),sprLtNI(1:(NAEC_Ltrl+1))
    real*8  :: EL_Ap,EL_Bs,EL_Lt
    real*8  :: l0_storage(1:N_spr)
    
    diff_A0andA=-1.00d30 ; PressStr=-1.00d30 ; TOL_Press=-1.00d30
    call get_diffA0andA_and_PressStr ; TOL_Press=TOL_fctr*PressStr
    
    sprApNI=-1 ; sprBsNI=-1 ; sprLtNI=-1 
    
    l0_storage(1:N_spr) = l0(1:N_spr)
    nspInACellTN=3 ; nspInACellNI=(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i = 1,Hlf_Ncell
       
       sprApTN=(i-1)*(nspInACellTN)+1;sprBsTN=(i-1)*(nspInACellTN)+2;sprLtTN=(i-1)*(nspInACellTN)+3 
       
       sprT=1 ; NumOFSgmnt=NAEC_Apcl+1 ; call get_the_NIsprs(i,sprT,NumOFSgmnt,sprApNI)
       sprT=2 ; NumOFSgmnt=NAEC_Bsal+1 ; call get_the_NIsprs(i,sprT,NumOFSgmnt,sprBsNI)
       sprT=3 ; NumOFSgmnt=NAEC_Ltrl+1 ; call get_the_NIsprs(i,sprT,NumOFSgmnt,sprLtNI)
       
       if (PCSorICS==1) then
          EL_Ap = Expected_ApLn_PCS(CMv,i)
          EL_Bs = Expected_BsLn_PCS(CMv,i)
          EL_Lt = Expected_LtLn_PCS(CMv,i)
          
       elseif (PCSorICS==2) then
          EL_Ap = Expected_ApLn_ICS(CMv,i)
          EL_Bs = Expected_BsLn_ICS(CMv,i)
          EL_Lt = Expected_LtLn_ICS(CMv,i)
       endif
       
       call disburse_EL_to_the_Sprs(sprApNI,sprBsNI,sprLtNI,EL_Ap,EL_Bs,EL_Lt)
       
    enddo
    
    call disburse_EL_to_nonSymm_sprs(CMv,PCSorICS)
    
    do i = 1,N_spr
       write(*,*) l0_storage(i),l0(i),i,"l0-strg,l0"
    enddo
    
    do i = 1,(nspInACellNI*Hlf_Ncell)
       write(*,*)l0(i),l0(i+(nspInACellNI*Hlf_Ncell)),(l0(i)-l0(i+(nspInACellNI*Hlf_Ncell))),i,"l0-symmtr"
    enddo
    
    
  end subroutine save_DiffA0A_set_l0_as_Expected_lengths
  
  
  subroutine change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
    implicit none
    integer :: i,imax
    integer :: cellL,cellR
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       
       if (i.ne.imax) then
          
          cellL     = i  ; cellR = cellL+Hlf_Ncell
          write(*,*) cellL,cellR,A0(cellL),A0(cellR),"bfr chnging cellL-cellR-A0L,A0R"
          
          A0(cellL) = diff_A0andA(i) + A(cellL)
          A0(cellR) = A0(cellL)
          write(*,*) cellL,cellR,A0(cellL),A0(cellR),"aft chnging cellL-cellR-A0L,A0R"
          
       elseif (i==imax) then
          
          cellL = N_cell
          write(*,*) cellL,A0(cellL),"bfr chnging cellL-A0L (imax)"
          
          A0(cellL) = diff_A0andA(imax) + A(cellL)
          write(*,*) cellL,A0(cellL),"aft chnging cellL-A0L (imax)"
       endif
       
    enddo
    
  end subroutine change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
  
  
  subroutine get_diffA0andA_and_PressStr
    implicit none
    integer :: i,imax
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       if (i.ne.imax) then
          diff_A0andA(i) = A0(i)-A(i)
          PressStr(i)    = k_area(i)*diff_A0andA(i)
       elseif (i == imax) then
          diff_A0andA(i) = A0(N_cell)-A(N_cell)
          PressStr(i)    = k_area(N_cell)*diff_A0andA(i)
       endif
       
       if (i.ne.imax) write(*,*) PressStr(i),diff_A0andA(i),A0(i),A(i),i,"diff_A0andA"
       if (i == imax) write(*,*) PressStr(i),diff_A0andA(i),A0(N_cell),A(N_cell),i,"diff_A0andA im"
    enddo
    
  end subroutine get_diffA0andA_and_PressStr
  
  
  subroutine keep_increasing_A0_to_make_press_within_TOL(CMv)
    implicit none
    integer, intent(in) :: CMv
    
    integer :: i,imax,j,jmax
    integer :: cellL,cellR
    real*8  :: CurrPress(1:(Hlf_Ncell+1))
    integer :: cnt_Press(1:(Hlf_Ncell+1))
    integer :: cnt_PressSumVal
    real*8  :: diffBtwn_TOLandCurrPress
    real*8  :: ZeroVal=0.0000d0
    integer :: diffIndctr=0
    
    open(unit=542,file='keepIncrA0_ToMAtchPress.dat',position='append')
    
    CurrPress=-1.0d30 ; call get_the_currPress(CurrPress)
    cnt_Press=0
    imax=Hlf_Ncell+1
    
    do 
       
       cnt_PressSumVal = 0 ; jmax = Hlf_Ncell+1
       
       do j = 1,jmax
          cnt_PressSumVal = cnt_PressSumVal + cnt_Press(j)
          write(542,*) cnt_Press(j),cnt_PressSumVal,j,"cnt_Press and Sum"
       enddo
       
       write(542,*) " "
       
       if (cnt_PressSumVal == (Hlf_Ncell+1)) exit
       
       
       do i = 1,imax
          
          if (i.ne.imax) then
             cellL=i ; cellR=cellL+Hlf_Ncell
          elseif (i==imax) then
             cellL=N_cell ; cellR=0
          endif
          
          if (cnt_Press(i) == 0) then
             
             if (CurrPress(i) .lt. TOL_Press(i)) then
                
                write(542,*) cellL,CurrPress(i),TOL_Press(i),"cellL-currP-TolP L" 
                
                if (i.ne.imax) then
                   A0(cellL) = 1.02d0*A0(cellL) ; A0(cellR)=A0(cellL)
                elseif (i == imax) then
                   A0(cellL) = 1.02d0*A0(cellL)
                endif
                
             elseif (CurrPress(i) .gt. TOL_Press(i)) then
                
                diffBtwn_TOLandCurrPress = TOL_Press(i) - CurrPress(i)
                
                if (diffBtwn_TOLandCurrPress.gt.ZeroVal) diffIndctr = +1
                if (diffBtwn_TOLandCurrPress.le.ZeroVal) diffIndctr = -1
                
                write(542,*) cellL,CurrPress(i),TOL_Press(i),diffIndctr,"cellL-currP-TolP-diffI G"
                cnt_Press(i) = 1
                
             endif
             
          elseif (cnt_Press(i) == 1) then
             continue
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       CurrPress=-1.0d30 ; call get_the_currPress(CurrPress)
       
       write(542,*) " "
       
    enddo
    
    close(542)
    
    
  end subroutine keep_increasing_A0_to_make_press_within_TOL

  subroutine get_the_currPress(CurrPress)
    implicit none
    real*8, intent(inout) :: CurrPress(1:(Hlf_Ncell+1))
    integer               :: imax

    imax =Hlf_Ncell+1
    
    do i = 1,imax
       if (i.ne.imax) CurrPress(i) = k_area(i)      * (A0(i)      - A(i))
       if (i == imax) CurrPress(i) = k_area(N_cell) * (A0(N_cell) - A(N_cell))
    enddo
    
  end subroutine get_the_currPress
  
  subroutine disburse_EL_to_the_Sprs(sprApNI,sprBsNI,sprLtNI,EL_Ap,EL_Bs,EL_Lt)
    implicit none
    integer, intent(in) :: sprApNI(1:(NAEC_Apcl+1)),sprBsNI(1:(NAEC_Bsal+1)),sprLtNI(1:(NAEC_Ltrl+1))
    real*8 , intent(in) :: EL_Ap,EL_Bs,EL_Lt
    
    integer :: i,j,jmax
    integer :: sprNmL,sprNmR
    integer :: nsprsInACellNI,sprNmIncrForRghtSide
    real*8  :: ELv,divFctr
    
    nsprsInACellNI       = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    sprNmIncrForRghtSide = (Hlf_Ncell)*(nsprsInACellNI)
    
    
    do i = 1,3
       
       if (i==1) then
          jmax = NAEC_Apcl+1  
          ELv  = EL_Ap
          
       elseif (i==2) then
          jmax = NAEC_Bsal+1
          ELv  = EL_Bs
          
       elseif (i==3) then
          jmax = NAEC_Ltrl+1
          ELv  = EL_Lt
       endif
          
       if (jmax==1)   divFctr=1.0000d0
       if (jmax.gt.1) divFctr=1.0000d0*(jmax)
       
       do j = 1,jmax
          
          if (i==1) sprNmL = sprApNI(j)
          if (i==2) sprNmL = sprBsNI(j)
          if (i==3) sprNmL = sprLtNI(j)
          
          sprNmR = sprNmL+sprNmIncrForRghtSide
          
          if (j==1)   l0(sprNmL) = ELv/divFctr     
          if (j.gt.1) l0(sprNmL) = l0(sprNmL-(j-1))
          
          l0(sprNmR) = l0(sprNmL)
          
       enddo
       
    enddo
    
  end subroutine disburse_EL_to_the_Sprs
  
  subroutine disburse_EL_to_nonSymm_sprs(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    integer :: symmSpr
    real*8  :: EL_Bs,divFctr
    integer :: i,imax
    
    if (CellsMeet==0) stop 'not for CM=0'
    
    symmSpr = (N_cell-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    if (PCSorICS==1) EL_Bs = Expected_BsLn_PCS(CMv,Hlf_Ncell+1) 
    if (PCSorICS==2) EL_Bs = Expected_BsLn_ICS(CMv,Hlf_Ncell+1)
    
    imax    = NAEC_Bsal+1
    divFctr = 1.0000d0*(imax) 
    
    do i = 1,imax
       if (i==1)   l0(symmSpr+i) = EL_Bs/divFctr   
       if (i.gt.1) l0(symmSpr+i) = l0(symmSpr+1)
    enddo
    
  end subroutine disburse_EL_to_nonSymm_sprs
  
  
  subroutine Control_State_match_expected_length_calctns_calls_togthr
    implicit none
    integer :: CMv,PCSorICS
    
    call get_CL_frm_Exp_CS
    call get_CL_frm_all_Sim_CS
    call get_EL_frm_taking_the_4thCell_asBase
    
    call measured_curvtrFrm_ExprmntlImg_and_Expctd_curvatr_Height
    
  end subroutine Control_State_match_expected_length_calctns_calls_togthr
  
  subroutine higher_kspr_tst_with_same_tension_length_and_adjstd_restLength(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    integer             :: i,j,jmax
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr)
    integer             :: numOfCellsToAdd
    
    write(*,*) CMv,PCSorICS,"CMv-PCSorICS"
    
    if (CMv.lt.1 .or. CMv.gt.4)  stop 'CMv not equal to 1 or 2 or 3 or 4'
    if (PCSorICS.ne.1) stop 'PCSorICS not equal to 1' 
    
    call parameter_values_for_module_CS_match
    call allocate_and_initialize_lenValues
    
    if (CMv==1 .and. PCSorICS==1) call routines_to_fix_CM1_PCS(CMv,PCSorICS)
    if (CMv==2 .and. PCSorICS==1) call routines_to_fix_CM2_PCS(CMv,PCSorICS)
    if (CMv==3 .and. PCSorICS==1) call routines_to_fix_CM3_PCS(CMv,PCSorICS)
    if (CMv==4 .and. PCSorICS==1) call routines_to_fix_CM4_PCS(CMv,PCSorICS)
    
  end subroutine higher_kspr_tst_with_same_tension_length_and_adjstd_restLength
  
  
  subroutine alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    real*8  :: lowest_press,PL1,PL2,PL3 ! PL = Pressure Level
    integer :: NumCellPL1,NumCellPL2,NumCellPL3    
    real*8  :: A0mA(1:Hlf_Ncell+1)
    integer :: Dcsn_A0chng(1:(Hlf_Ncell+1))
    integer :: BsORLt,cellNum
    integer :: cellNum1,cellNum2,cntLp
    real*8  :: pres1,pres2
    
    if (CMv==2 .and. PCSorICS==-1) then
       
       A0mA(1:(Hlf_Ncell+1)) = -1.0d30                ;  Dcsn_A0chng(1:(Hlf_Ncell+1)) = -1
       lowest_press          = 0.089697071251372371d0 !  read frm NI_modelInitiatinAW037.dat
       
       PL1  = 0.100d0      ; NumCellPL1 = 5
       PL2  = 0.095d0      ; NumCellPL2 = 6
       PL3  = lowest_press ; NumCellPL3 = 8
       
       !call matchBiggstCurve_and_Achieving_HighLowRngePress(CMv,PCSorICS)
       call altrPresUsingA0mA_CM2PCS(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3,Dcsn_A0chng,A0mA)
       call Equilibrate_only_NI_model
       call altrPresUsingTolrnc_CM2PCS(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3)
       
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       
    elseif (CMv==2 .and. PCSorICS==1) then
       
       !call altrPresUsingTolrnc_generalized(CMv,PCSorICS)
       call altrCurvtrUsingTolrnc_generalized(CMv,PCSorICS)
       
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       call rscale_basedOn_bndryPress() ; write(*,*) "AFT REscaled at PULL-PUSH CM2PCS"
       
    elseif (CMv==3 .and. PCSorICS==1) then
       
       !call altrPresUsingTolrnc_generalized(CMv,PCSorICS)
       call altrCurvtrUsingTolrnc_generalized(CMv,PCSorICS)
       
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       call rscale_basedOn_bndryPress() ; write(*,*) "AFT REscaled at PULL-PUSH CM3PCS"
       
    elseif (CMv==1 .and. PCSorICS==1) then
       
       !call altrPresUsingTolrnc_generalized(CMv,PCSorICS)
       call altrCurvtrUsingTolrnc_generalized(CMv,PCSorICS)
       
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       call rscale_basedOn_bndryPress() ; write(*,*) "AFT REscaled at PULL-PUSH"
       
    elseif (CMv==4 .and. PCSorICS==1) then
       
       !call altrPresUsingTolrnc_generalized(CMv,PCSorICS)
       
       call altrCurvtrUsingTolrnc_generalized(CMv,PCSorICS)
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
       call print_relevant_info_Actual_and_ExpctdLen(CMv)
       call rscale_basedOn_bndryPress() ; write(*,*) "AFT REscaled at PULL-PUSH CMv4"
       
       write(*,*) "CM4 DONE"
       !stop 'CM4 DONE'
       
    endif
    
    BsORLt=2 ; call get_Simulted_Curvtr_Height(CMv,PCSorICS,BsORLt)
    
  end subroutine alter_pressure_to_match_the_exp_img
  
  subroutine altrPresUsingTolrnc_CM2PCS(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3)
    implicit none
    real*8,  intent(in)  :: PL1,PL2,PL3
    integer, intent(in)  :: NumCellPL1,NumCellPL2,NumCellPL3
    integer, allocatable :: DistntCells(:)
    
    integer :: i,imax,j,jmax
    integer :: NxtToDistntCell,cellNmToAltr
    integer :: cell1L,cell2L
    integer :: numofcellstoadjst
    integer :: Dirctn_dcsn,itrnNum
    
    allocate(DistntCells(1:NumCellPL1)) ; DistntCells = -1
    call get_diffPLcells_CM2PCS(NumCellPL1,DistntCells,NxtToDistntCell)
     
    !call altrPressMulipleCellsSamePressLevel_1by1(PL1,NumCellPL1,DistntCells)
    !call altrPressSingleCell(PL2,NxtToDistntCell)
    !call altrPressSingleCell(PL3,NxtToDistntCell+1)
    !call altrPressSingleCell(PL3,NxtToDistntCell+2)
    
    Frame_NI = 142
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    !call altrPressMulipleCellsSamePressLevel_1by1(PL1,NumCellPL1,DistntCells)
    !call altrPressSingleCell(PL2,NxtToDistntCell)
    !call altrPressSingleCell(PL3,NxtToDistntCell+1)
    !call altrPressSingleCell(PL3,NxtToDistntCell+2)
    
    Frame_NI = 206
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    !cellNmToAltr=1 ; call altrPressSingleCell(PL1,cellNmToAltr)
    !cellNmToAltr=2 ; call altrPressSingleCell(PL1,cellNmToAltr)
    !cellNmToAltr=3 ; call altrPressSingleCell(PL1,cellNmToAltr)
    !cellNmToAltr=4 ; call altrPressSingleCell(PL1,cellNmToAltr)
    !cellNmToAltr=5 ; call altrPressSingleCell(PL1,cellNmToAltr)
    
    Frame_NI = 240
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% THE FIRST SIX CELLS WITH EQUIVALENT PRESSURE
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! make_the_first_six_cells_Pressure_close
    Frame_NI = 337
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% THE FIRST NINE CELLS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !numofcellstoadjst = 9
    !call make_cellPress_inDefinedRnge(numofcellstoadjst)
    Frame_NI = 479
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
  end subroutine altrPresUsingTolrnc_CM2PCS
  
  
  
  subroutine altrPresUsingTolrnc_generalized(CMv,PCSorICS)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    
    integer :: i,imax,j,jmax
    integer :: cellNmToAltr
    integer :: cell1L,cell2L
    integer :: numofcellstoadjst
    real*8  :: TOL_prcnt
    integer :: Dirctn_dcsn,itrnNum
    real*8  :: currBndryPress,rsFctr
    integer :: PressChoice,CntChoice
    
    if (CMv==2 .and. PCSorICS==1) then
       numofcellstoadjst = 12
       PressChoice       = 3
       CntChoice         = 2
       TOL_prcnt         = 5.0000d0
       
    elseif (CMv==3 .and. PCSorICS==1) then
       numofcellstoadjst = 12
       PressChoice       = 3
       CntChoice         = 2
       TOL_prcnt         = 5.0000d0
       
    elseif (CMv==1 .and. PCSorICS==1) then
       numofcellstoadjst = 12
       PressChoice       = 3
       CntChoice         = 2
       TOL_prcnt         = 5.0000d0
       
    elseif (CMv==4 .and. PCSorICS==1) then
       numofcellstoadjst = 12
       PressChoice       = 3
       CntChoice         = 2
       TOL_prcnt         = 5.0000d0
       
    endif
    
    call make_cellPress_inDefinedRnge_generalized(CMv,PCSorICS,numofcellstoadjst,&
         PressChoice,CntChoice,TOL_prcnt)
    call rscale_basedOn_bndryPress() ; write(*,*) "AFT REscaled at pos 01"
    !call print_the_curvtr_params(CMv,LowRngeCurvtr,AchievingCurvtr,HighRngeCurvtr)
    
  end subroutine altrPresUsingTolrnc_generalized
  
  
  subroutine altrCurvtrUsingTolrnc_generalized(CMv,PCSorICS)
    implicit none
    integer,intent(in) ::CMv,PCSorICS
    real*8             ::rsFctr,currBndryPress
    real*8             ::AchievingCurvtr(1:Hlf_Ncell),LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    integer            ::numCurvtr
    integer,allocatable::Curvtrs(:)
    
    call make_curvtr_inDefinedRnge_generalized(CMv,PCSorICS);write(*,*) Frame_NI-1,CMv,"Frm-CMv at altrCurv"
    
    call get_Achieving_HighLowRnge_Curvtr(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,HighRngeCurvtr)
    call print_the_curvtr_params(CMv,LowRngeCurvtr,AchievingCurvtr,HighRngeCurvtr)
    call rscale_basedOn_bndryPress() ; write(*,*) "Rescaled in altrCurvtr"
    
  end subroutine altrCurvtrUsingTolrnc_generalized
  
  subroutine fix_one_curvtr(curvtrNo,curvtrVal)
    implicit none
    integer, intent(in) :: curvtrNo
    real*8 , intent(in) :: curvtrVal
    
  end subroutine fix_one_curvtr
  
  subroutine altrPressMulipleCellsSamePressLevel_1by1(PL,NumCellPL,cellNums)
    implicit none
    real*8,  intent(in) :: PL
    integer, intent(in) :: NumCellPL
    integer, intent(in) :: cellNums(1:NumCellPL)
    
    integer :: i,j,imax,jmax,cntLp
    real*8  :: TP1L,TP2L,TP1G,TP2G,CP
    integer :: cellL,cellR
    
    TP1L = (TolPrcnt1L)*PL ; TP2L = (TolPrcnt2L)*PL
    TP1G = (TolPrcnt1G)*PL ; TP2G = (TolPrcnt2G)*PL
    
    write(*,*) TP1L,TP2L,TP1G,TP2G,"TP vals"
    
    imax = NumCellPL
    if (imax .gt. Hlf_NCell) stop 'imax cant be greater than Hlf_Ncell here'
    
    
    do i = 1,NumCellPL
       
       cellL = cellNums(i)
       cntLp = 1
       
       CP = k_area(cellL) * (A0(cellL)-A(cellL))
       
       if (CP .le. PL) call doloop_for_IncrsPress(cellL,PL,TP1L,TP2L,cntLp)
       if (CP .gt. PL) call doloop_for_DecrsPress(cellL,PL,TP1G,TP2G,cntLp)
       
    enddo
    
  end subroutine altrPressMulipleCellsSamePressLevel_1by1
  
  
  subroutine altrPressSingleCell(PL,singlCellNum)
    implicit none
    real*8,  intent(in) :: PL
    integer, intent(in) :: singlCellNum
    
    integer :: i,j,imax,jmax,cntLp
    real*8  :: TP1L,TP2L,TP1G,TP2G,CP
    integer :: cellL
    
    TP1L = (TolPrcnt1L)*PL ; TP2L = (TolPrcnt2L)*PL
    TP1G = (TolPrcnt1G)*PL ; TP2G = (TolPrcnt2G)*PL
    
    write(*,*) TP1L,TP2L,TP1G,TP2G,"TP vals"
    
    cellL = singlCellNum
    cntLp = 1
       
    CP = k_area(cellL) * (A0(cellL)-A(cellL))
    
    if (CP .le. PL) call doloop_for_IncrsPress(cellL,PL,TP1L,TP2L,cntLp)
    if (CP .gt. PL) call doloop_for_DecrsPress(cellL,PL,TP1G,TP2G,cntLp)
    
  end subroutine altrPressSingleCell
  
  
  subroutine altrPresUsingA0mA_CM2PCS(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3,&
       Dcsn_A0chng,A0mA)
    implicit none
    real*8,  intent(in)    :: PL1,PL2,PL3                          ! PL = Pressure Level  
    integer, intent(in)    :: NumCellPL1,NumCellPL2,NumCellPL3  
    integer, intent(inout) :: Dcsn_A0chng(1:(Hlf_Ncell+1))
    real*8 , intent(inout) :: A0mA(1:(Hlf_Ncell+1))
    
    call get_A0mA_and_DcsnA0(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3,Dcsn_A0chng,A0mA)
    call get_A0_vals_using_A0mA_and_DcsnA0(Dcsn_A0chng,A0mA)
    call Equilibrate_only_NI_model
    
  end subroutine altrPresUsingA0mA_CM2PCS
  
  
  subroutine get_A0mA_and_DcsnA0(PL1,PL2,PL3,NumCellPL1,NumCellPL2,NumCellPL3,Dcsn_A0chng,A0mA)
    implicit none
    real*8,  intent(in)    :: PL1,PL2,PL3                          ! PL = Pressure Level  
    integer, intent(in)    :: NumCellPL1,NumCellPL2,NumCellPL3  
    integer, intent(inout) :: Dcsn_A0chng(1:(Hlf_Ncell+1))
    real*8 , intent(inout) :: A0mA(1:(Hlf_Ncell+1))
    
    integer :: i,imax,j,jmax
    integer :: lim1,lim2,lim3,lim4
    
    imax = Hlf_Ncell+1
    
    lim1 = NumCellPL1
    lim2 = NumCellPL2
    lim3 = NumCellPL3
    lim4 = imax
    
    do i = 1,imax
       
       if (i.le.lim1) then
          
          Dcsn_A0chng(i) = 1
          A0mA(i)        = PL1/k_area(i)
          write(*,*) Dcsn_A0chng(i),A0mA(i),PL1,k_area(i),i,"A0mA1"
          
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          
          Dcsn_A0chng(i) = 1
          A0mA(i)        = PL2/k_area(i)
          write(*,*) Dcsn_A0chng(i),A0mA(i),PL2,k_area(i),i,"A0mA2"
          
       elseif ((i.gt.lim2) .and. (i.le.lim3)) then
          
          if ((i-lim2)==1) then
             
             Dcsn_A0chng(i) = 0
             A0mA(i)        = A0(i)-A(i)
             write(*,*) Dcsn_A0chng(i),A0mA(i),PL3,k_area(i),i,"A0mA31"
             
          elseif ((i-lim2)==2) then
             
             Dcsn_A0chng(i) = 1
             A0mA(i)        = PL3/k_area(i)
             write(*,*) Dcsn_A0chng(i),A0mA(i),PL3,k_area(i),i,"A0mA32"
             
          endif
          
       elseif ((i.gt.lim3) .and. (i.le.lim4)) then
          
          Dcsn_A0chng(i) = 0
          
          if (i.ne.imax) A0mA(i) = A0(i)      - A(i)
          if (i == imax) A0mA(i) = A0(N_cell) - A(N_cell)
          
          write(*,*) Dcsn_A0chng(i),A0mA(i),i,"A0mA4"
          
       endif
       
    enddo
    
  end subroutine get_A0mA_and_DcsnA0
  
  subroutine get_A0_vals_using_A0mA_and_DcsnA0(Dcsn_A0chng,A0mA)
    implicit none
    integer, intent(in) :: Dcsn_A0chng(1:(Hlf_Ncell+1))
    real*8 , intent(in) :: A0mA(1:(Hlf_Ncell+1))
    
    real*8  :: A0_str(1:(Hlf_Ncell+1))
    real*8  :: curr_A0mA(1:(Hlf_Ncell+1))
    integer :: cellL,cellR
    integer :: i,imax
    
    A0_str(1:Hlf_Ncell) = A0(1:Hlf_Ncell)
    A0_str(Hlf_Ncell+1) = A0(N_cell)
    
    curr_A0mA(1:Hlf_Ncell) = A0(1:Hlf_Ncell) - A(1:Hlf_Ncell)
    curr_A0mA(Hlf_Ncell+1) = A0(N_cell)      - A(N_cell)
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       
       if (i.ne.imax) cellL = i 
       if (i == imax) cellL = N_cell
       if (i.ne.imax) cellR = cellL+Hlf_Ncell
       if (i == imax) cellR = 0
       
       if (Dcsn_A0chng(i) == 0) then
          continue
       elseif (Dcsn_A0Chng(i) == 1) then
          
          !A0(cellL) = A0(cellL) + A0mA(i)
          
          A0(cellL) = A0(cellL) + (A0mA(i)-curr_A0mA(i))
          
          if (i.ne.imax) A0(cellR) = A0(cellL)
          if (i == imax) continue
          
          write(*,*) i,A0_str(i),A0(cellL),(A0(cellL)-A0_str(i)),A0mA(i),"A0 bfr aft"
          
       endif
       
    enddo
    
  end subroutine get_A0_vals_using_A0mA_and_DcsnA0
  
  subroutine changeNodeY2_of_bndry_frm_4thCell_ltrl_mmbne_Y10
    implicit none
    integer :: cellNmbr,nodeNmbr
    integer :: bndryNodeL,bndryNodeR
    
    cellNmbr     = 4
    nodeNmbr     = (2*cellNmbr) + 2
    
    yvalBotmLayr = node_xy(nodeNmbr,2) ; write(*,*) yvalBotmLayr,"yvalBotmLayr"
    bndryNodeL   = 2                   ; bndryNodeR = bndryNodeL + (Hlf_Ncell+1)*2
    
    node_xy(bndryNodeL,2) = yvalBotmLayr
    node_xy(bndryNodeR,2) = node_xy(bndryNodeL,2)
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    call Equilibrate_only_NI_model
    
  end subroutine changeNodeY2_of_bndry_frm_4thCell_ltrl_mmbne_Y10
  
  subroutine get_bndryBotm_yNode_frm_CM2PCS
    implicit none
    integer :: lftBndryBotm,rghtBndryBotm
    
    lftBndryBotm=2 ; rghtBndryBotm=2*(Hlf_Ncell+1)+2
    
    node_xy(lftBndryBotm,2)  = -8.2101593017578125000d0    
    node_xy(rghtBndryBotm,2) = -8.2101593017578125000d0
    
    call nodes_to_coordntes(node_xy,coordntes_xy)
    
  end subroutine get_bndryBotm_yNode_frm_CM2PCS
  
  subroutine print_actual_and_expected_lngth
    implicit none
    integer :: nsprsInACellNI
    integer :: i,j,jmax
    integer :: remaindr,cellNm
    real*8  :: comp1,comp2,prcntV
    
    nsprsInACellNI = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    
    !open(unit=155,file='actual_expected_len.dat',position='append')
    !open(unit=156,file='Remaindr_chk.dat')
    
    do i = 1,(nsprsInACellNI*Hlf_Ncell)
       
       remaindr = mod(i,nsprsInACellNI)
       cellNm   = ((i-1)/nsprsInACellNI)+1
       !write(156,*) remaindr,cellNm,i
       
       if (remaindr==1) then
          comp1 = l(i)
          comp2 = Expected_ApLn_PCS(2,cellNm)
          
       elseif (remaindr==2.or.remaindr==3.or.remaindr==4) then
          comp1 = l(i) 
          comp2 = (Expected_BsLn_PCS(2,cellNm)/3.000d0)
          
       elseif (remaindr==5.or.remaindr==6.or.remaindr==0) then
          comp1 = l(i) 
          comp2 = (Expected_LtLn_PCS(2,cellNm)/3.000d0)
       endif
       
       prcntV = ((abs(comp1-comp2))/comp2)*100.0000d0
       
       !write(155,*) comp1,comp2,prcntV,i,"len vartn than Expected"
       
    enddo
    
    !write(155,*) " "
    !close(155)
    !close(156)
    
  end subroutine print_actual_and_expected_lngth
  
  
  subroutine print_relevant_info_Actual_and_ExpctdLen(CMv)
    implicit none
    integer, intent(in) :: CMv
    
    integer :: i,imax
    integer :: ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1))
    integer :: nsprsInACell
    real*8  :: ApLn,BsLn,LtLn
    real*8  :: prcntAp,prcntBs,prcntLt
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    imax         = Hlf_Ncell
    
    open(unit=332,file='diff_A0A_AreaPrpchk.dat',position='append')
    open(unit=334,file='diff_A0A_SprsPrpChk.dat',position='append')
    
    do i = 1,imax
       
       ApclSp(1) = (i-1)*(nsprsInACell) + 1  
       BsalSp(1) = (i-1)*(nsprsInACell) + 2 ; BsalSp(2) = BsalSp(1)+1 ; BsalSp(3) = BsalSp(1)+2 
       LtrlSp(1) = BsalSp(3)+1 ; LtrlSp(2) = BsalSp(3)+2 ; LtrlSp(3) = BsalSp(3)+3
       
       ApLn = l(ApclSp(1))
       BsLn = l(BsalSp(1)) + l(BsalSp(2)) + l(BsalSp(3))
       LtLn = l(LtrlSp(1)) + l(LtrlSp(2)) + l(LtrlSp(3)) 
       
       prcntAp = (abs(ApLn-Expected_ApLn_PCS(CMv,i))/Expected_ApLn_PCS(CMv,i))*100.00d0
       prcntBs = (abs(BsLn-Expected_BsLn_PCS(CMv,i))/Expected_BsLn_PCS(CMv,i))*100.00d0
       prcntLt = (abs(LtLn-Expected_LtLn_PCS(CMv,i))/Expected_LtLn_PCS(CMv,i))*100.00d0
       
       write(332,*) k_area(i),A0(i),A(i),(A0(i)-A(i)),(k_area(i)*(A0(i)-A(i))),i,"A0A-diff"
       write(334,*) ApLn,Expected_ApLn_PCS(CMv,i),prcntAp,"Ap"
       write(334,*) BsLn,Expected_BsLn_PCS(CMv,i),prcntBs,"Bs"
       write(334,*) LtLn,Expected_LtLn_PCS(CMv,i),prcntLt,"Lt"
       write(334,*) " "
       
    enddo
    
    BsalSp(1) = (N_cell-1)*(nsprsInACell)+1 ; BsalSp(2)=BsalSp(1)+1 ; BsalSp(3)=BsalSp(1)+2
    BsLn      = l(BsalSp(1)) + l(BsalSp(2)) + l(BsalSp(3))
    prcntBs   = (abs(BsLn-Expected_BsLn_PCS(CMv,i))/Expected_BsLn_PCS(CMv,i))*100.00d0
    
    write(332,*) k_area(N_cell),A0(N_cell),A(N_cell),(A0(N_cell)-A(N_cell)),&
         (k_area(N_cell)*(A0(N_cell)-A(N_cell))),imax+1,"A0A-diff"
    write(332,*) " "
    
    write(334,*) BsLn,Expected_BsLn_PCS(CMv,imax+1),prcntBs,"Bs"
    write(334,*) " "
    
    close(332)
    close(334)
    
  end subroutine print_relevant_info_Actual_and_ExpctdLen
  
  
  subroutine get_diffPLcells_CM2PCS(NumCellPL1,DistntCells,NxtToDistntCell)
    implicit none
    integer, intent(in)    :: NumCellPL1
    integer, intent(inout) :: DistntCells(1:NumCellPL1)
    integer, intent(inout) :: NxtToDistntCell
    integer                :: i,j
     
    do i = 1,NumCellPL1
       DistntCells(i) = i
    enddo
    
    NxtToDistntCell = DistntCells(NumCellPL1) + 1 
    
  end subroutine get_diffPLcells_CM2PCS
  
  
  subroutine doloop_for_IncrsPress(cellL,PL,TP1,TP2,cntLp)
    implicit none
    integer, intent(in)    :: cellL
    real*8 , intent(in)    :: PL,TP1,TP2
    integer, intent(inout) :: cntLp
    
    real*8  :: CP
    integer :: cellR
    
    cellR = cellL+Hlf_Ncell ; write(*,*) cellL,cellR,"cellL-cellR"
    
    do 
       
       CP = k_area(cellL)*(A0(cellL)-A(cellL)) ; write(*,*) CP,PL,A0(cellL),A(cellL),cntLp,"cntLp"
       
       if (CP .le. PL) then
          
          
          if (CP .le. TP1) then
             
             A0(cellL) = 1.01d0*A0(cellL)
             A0(cellR) = A0(cellL)
             
          elseif ((CP .gt. TP1) .and. (CP .le. TP2)) then
             
             A0(cellL) = 1.005d0*A0(cellL) 
             A0(cellR) = A0(cellL)
             
          elseif ((CP .gt. TP2) .and. (CP .le. PL)) then
             write(*,*) "bracketd press 1I",CP,PL,TP1,TP2
             exit
          endif
          
          call Equilibrate_only_NI_model
          
       elseif (CP .gt. PL) then
          write(*,*) "bracketd press 2I",CP,PL,TP1,TP2
          exit
       endif
       
       cntLp = cntLp+1
       
    enddo
    
    
  end subroutine doloop_for_IncrsPress
  
  
  subroutine doloop_for_DecrsPress(cellL,PL,TP1,TP2,cntLp)
    implicit none
    integer, intent(in)    :: cellL
    real*8 , intent(in)    :: PL,TP1,TP2
    integer, intent(inout) :: cntLp
    
    real*8  :: CP
    integer :: cellR
    
    cellR = cellL+Hlf_Ncell ; write(*,*) cellL,cellR,"cellL-cellR"
    
    
    do 
       
       CP = k_area(cellL)*(A0(cellL)-A(cellL)) ; write(*,*) CP,PL,A0(cellL),A(cellL),cntLp,"cntLp"
       
       if (CP .gt. PL) then
          
          
          if (CP .gt. TP1) then
             
             A0(cellL) = 0.99d0*A0(cellL)
             A0(cellR) = A0(cellL)
             
          elseif ((CP .lt. TP1) .and. (CP .gt. TP2)) then
             
             A0(cellL) = 0.995d0*A0(cellL) 
             A0(cellR) = A0(cellL)
             
          elseif ((CP .lt. TP2) .and. (CP .gt. PL)) then
             write(*,*) "bracketd press 1D",CP,PL,TP1,TP2
             exit
          endif
          
          call Equilibrate_only_NI_model
          
       elseif (CP .le. PL) then
          write(*,*) "bracketd press 2D",CP,PL,TP1,TP2
          exit
       endif
       
       cntLp = cntLp+1
       
    enddo
    
    
  end subroutine doloop_for_DecrsPress
  
  
  subroutine make_two_cells_equal_press(cell1L,cell2L)
    implicit none
    integer, intent(in) :: cell1L,cell2L
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    IDEA IS TO MAKE THE HIGH PRESSURE TO LOW PRESSURE    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    integer :: cell1R,cell2R
    integer :: LPC,HPC                    !LPC = Low Pressure Cell, HPC = High Pressure Cell
    real*8  :: LPV,HPV,LPV_init,HPV_init  !LPV = Low Pressure Value, HPV = High Pressure Value
    integer :: i,imax,j,jmax,cntLp
    real*8  :: Press1,Press2,A0_str
    real*8  :: TolPressLvl1,TolPressLvl2
    real*8  :: chk_diff
    
    cell1R = cell1L + Hlf_Ncell
    cell2R = cell2L + Hlf_Ncell
    
    Press1 = k_area(cell1L) * (A0(cell1L)-A(cell1L))
    write(*,*) Press1,k_area(cell1L),A0(cell1L),A(cell1L),"P1" 
    Press2 = k_area(cell2L) * (A0(cell2L)-A(cell2L))
    write(*,*) Press2,k_area(cell2L),A0(cell2L),A(cell2L),"P2"
    
    if (Press1 .le. Press2) then
       LPC      = cell1L ; HPC      = cell2L
       LPV      = Press1 ; HPV      = Press2
       LPV_init = LPV    ; HPV_init = HPV
       
    elseif (Press1 .gt. Press2) then
       LPC      = cell2L ; HPC      = cell1L
       LPV      = Press2 ; HPV      = Press1
       LPV_init = LPV    ; HPV_init = HPV
    endif
    
    TolPressLvl1 = (1.015)*(LPV)
    TolPressLvl2 = (1.005)*(LPV)
    
    write(*,*) " "
    
    cntLp = 1
    
    do 
       
       if (HPV .gt. LPV) then
          
          chk_diff = (abs(HPV-LPV)/LPV_init)*100.0d0 ; write(*,*) chk_diff,"chk_diff"
          
          if (chk_diff .le. 1.00d0) exit
          
          if (HPV.gt.LPV_init) then
             
             A0(HPC)           = 0.995d0*A0(HPC)
             A0(HPC+Hlf_Ncell) = A0(HPC)
             
          elseif (HPV.lt.LPV_init) then
             
             A0(LPC)           = 1.005d0*A0(LPC)
             A0(LPC+Hlf_Ncell) = A0(LPC)
             
          endif
          
       endif
              
       call Equilibrate_only_NI_model
       write(*,*) cntLp,Frame_NI-1,"FrmNI-cntLp"
       cntLp = cntLp+1
       
       LPV = k_area(LPC) * (A0(LPC)  - A(LPC)) 
       HPV = k_area(HPC) * (A0(HPC) - A(HPC))
       
    enddo
    
    call Equilibrate_only_NI_model
    write(*,*) cntLp,Frame_NI-1,"FrmNI-cntLp"
    
    LPV = k_area(LPC) * (A0(LPC)  - A(LPC)) 
    HPV = k_area(HPC) * (A0(HPC) - A(HPC))
    
  end subroutine make_two_cells_equal_press
    
  
  subroutine make_the_first_six_cells_Pressure_close
    implicit none
    integer :: cell1,cell2,cell3,cell4,cell5,cell6
    real*8  :: AchievingPress,LowRngePress,HighRngePress
    integer :: NcellsInsideLimitPress
    integer :: CellPressPostn(1:6)
    integer :: i,imax,j,jmax,cntLp
    integer :: cellL,cellR
    
    cell1=1 ; cell2=2 ; cell3=3 ; cell4=4 ; cell5=5 ; cell6=6
    
    AchievingPress = 0.1000d0
    LowRngePress   = AchievingPress - 0.005d0*AchievingPress
    HighRngePress  = AchievingPress + 0.005d0*AchievingPress
    
    cntLp = 1
    
    do
       
       NcellsInsideLimitPress=0 ; CellPressPostn=-10
       call get_cntCellInsideLimit(LowRngePress,HighRngePress,CellPressPostn,NcellsInsideLimitPress)
       
       if (NcellsInsideLimitPress == 6) exit
       
       imax = 6
       
       do i = 1,imax
          
          if (i==1) cellL = cell1
          if (i==2) cellL = cell2
          if (i==3) cellL = cell3
          if (i==4) cellL = cell4
          if (i==5) cellL = cell5
          if (i==6) cellL = cell6
          
          cellR = cellL+Hlf_Ncell
          
          if (CellPressPostn(i) == 0 ) then
             continue
          elseif (CellPressPostn(i) == +1) then
             
             A0(cellL) = 0.9975d0*A0(cellL) ! 0.995
             A0(cellR) = A0(cellL)
             
          elseif (CellPressPostn(i) == -1) then
             
             A0(cellL) = 1.0025d0*A0(cellL) ! 1.005
             A0(cellR) = A0(cellL)
             
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
    enddo
    
  end subroutine make_the_first_six_cells_Pressure_close
  
  
  subroutine get_cntCellInsideLimit(LowRngePress,HighRngePress,CellPressPostn,NcellsInsideLimitPress)
    implicit none
    real*8,  intent(in)    :: LowRngePress,HighRngePress
    integer, intent(inout) :: CellPressPostn(1:6)
    integer, intent(inout) :: NcellsInsideLimitPress
    
    real*8  :: Press1,Press2,Press3,Press4,Press5,Press6,PressVal
    integer :: i,j,imax,jmax
    
    Press1 = k_area(1)*(A0(1)-A(1)) ; Press2 = k_area(2)*(A0(2)-A(2))
    Press3 = k_area(3)*(A0(3)-A(3)) ; Press4 = k_area(4)*(A0(4)-A(4))
    Press5 = k_area(5)*(A0(5)-A(5)) ; Press6 = k_area(6)*(A0(6)-A(6))
    
    imax = 6
    
    do i = 1,imax
       
       if (i==1) PressVal = Press1
       if (i==2) PressVal = Press2
       if (i==3) PressVal = Press3
       if (i==4) PressVal = Press4
       if (i==5) PressVal = Press5
       if (i==6) PressVal = Press6
       
       if ((PressVal.gt.LowRngePress) .and. (PressVal.le.HighRngePress)) then
          CellPressPostn(i) =  0
       elseif (PressVal.gt.HighRngePress) then
          CellPressPostn(i) = +1
       elseif (PressVal.le.LowRngePress) then
          CellPressPostn(i) = -1
       endif
       
       if (CellPressPostn(i) == 0) NcellsInsideLimitPress = NcellsInsideLimitPress+1
       
    enddo
    
    write(*,*) CellPressPostn(1:imax),"CellPressPostn"
    write(*,*) NcellsInsideLimitPress,"NcellsInsideLimitPress"
    
  end subroutine get_cntCellInsideLimit
  

  subroutine get_signTrnsfrmtnInfo()
    implicit none
    integer :: currnt_CurvtrSign(1:Hlf_Ncell)
    
    call get_currnt_CurvtrSign(currnt_CurvtrSign)
    
    
  end subroutine get_signTrnsfrmtnInfo
  
  subroutine get_currnt_CurvtrSign(currnt_CurvtrSign)
    implicit none
    integer, intent(out) :: currnt_CurvtrSign(1:Hlf_Ncell)
    
    integer :: i,imax
    integer :: cellBfr,cellAft
    real*8  :: pressBfr,pressAft
    
    open(unit=874,file='CurrCurvtrSign.dat')
    
    imax = Hlf_Ncell ; currnt_CurvtrSign(1:Hlf_Ncell) = -10
    
    do i = 1,imax
       
       cellBfr = i
       
       if (i.ne.imax) cellAft = i+1
       if (i == imax) cellAft = N_cell
       write(874,*) 'cellBfr =',cellBfr,'cellAft =',cellAft,'for curvature =',i
       
       pressBfr = k_area(cellBfr) * (A0(cellBfr)-A(cellBfr))
       pressAft = k_area(cellAft) * (A0(cellAft)-A(cellAft))
       write(874,*) "press Bfr-Aft =",pressBfr,pressAft,"for curvature =",i
       
       if (pressBfr .gt. pressAft) currnt_CurvtrSign(i) = +1
       if (pressBfr .lt. pressAft) currnt_CurvtrSign(i) = -1
       
    enddo
    
    close(874)
    
  end subroutine get_currnt_CurvtrSign
  
  
  subroutine get_cntCellInsideLimit_Gnrl(numOfCellsToAdjst,LowRngePress,HighRngePress,&
       CellPressPostn,NcellsInsideLimitPress)
    implicit none
    integer, intent(in)    :: numOfCellsToAdjst
    real*8,  intent(in)    :: LowRngePress(1:numOfCellsToAdjst),HighRngePress(1:numOfCellsToAdjst)
    integer, intent(inout) :: CellPressPostn(1:numOfCellsToAdjst)
    integer, intent(inout) :: NcellsInsideLimitPress
    
    real*8  :: PressValue(1:numOfCellsToAdjst),PressVal
    real*8  :: LowRngePressVal,HighRngePressVal
    integer :: i,j,imax,jmax
    
    imax = numOfCellsToAdjst
    
    do i = 1,imax
       if (i.le.Hlf_Ncell) PressValue(i) = k_area(i)*(A0(i)-A(i))
       if (i.gt.Hlf_Ncell) PressValue(i) = k_area(N_cell)*(A0(N_cell)-A(N_cell))
    enddo
    
    
    do i = 1,imax
       
       PressVal         = PressValue(i)
       LowRngePressVal  = LowRngePress(i)
       HighRngePressVal = HighRngePress(i)
       
       if ((PressVal.gt.LowRngePressVal) .and. (PressVal.le.HighRngePressVal)) then
          CellPressPostn(i) =  0
       elseif (PressVal.gt.HighRngePressVal) then
          CellPressPostn(i) = +1
       elseif (PressVal.le.LowRngePressVal) then
          CellPressPostn(i) = -1
       endif
       
       if (CellPressPostn(i) == 0) NcellsInsideLimitPress = NcellsInsideLimitPress+1
       
    enddo
    
    write(*,*) CellPressPostn(1:imax),"CellPressPostn"
    write(*,*) NcellsInsideLimitPress,"NcellsInsideLimitPress"
    
  end subroutine get_cntCellInsideLimit_Gnrl
  
  
  subroutine get_cntCell_BasedOn_PressOrCurvtrOrientn(CMv,numOfCellsToAdjst,LowRngePress,HighRngePress,&
       CellPressPostn,NcellsInsideLimitPress,NcellsInsideLimitCurvtr)
    implicit none
    integer, intent(in)    :: CMv,numOfCellsToAdjst
    real*8,  intent(in)    :: LowRngePress(1:numOfCellsToAdjst),HighRngePress(1:numOfCellsToAdjst)
    integer, intent(inout) :: CellPressPostn(1:numOfCellsToAdjst)
    integer, intent(inout) :: NcellsInsideLimitPress,NcellsInsideLimitCurvtr
    
    real*8                 :: PressValue(1:numOfCellsToAdjst),PressVal,pressBfr,pressAft
    real*8                 :: LowRngePressVal,HighRngePressVal,CellPressPostnVal
    integer                :: i,j,imax,jmax
    integer                :: curvtrNum,cellBfr,cellAft
    integer                :: ExprmntlCurvtrSign,CurrntCurvtrSign,CntCurvtr
    
    imax = numOfCellsToAdjst
    
    do i = 1,imax
       if (i.le.Hlf_Ncell) PressValue(i) = k_area(i)*(A0(i)-A(i))
       if (i.gt.Hlf_Ncell) PressValue(i) = k_area(N_cell)*(A0(N_cell)-A(N_cell))
    enddo
    
    do i = 1,imax
       PressVal=PressValue(i) ; LowRngePressVal=LowRngePress(i) ; HighRngePressVal=HighRngePress(i)
       call find_ifPressureInthe_Rnge(PressVal,LowRngePressVal,HighRngePressVal,CellPressPostnVal)
       CellPressPostn(i) = CellPressPostnVal
       if (CellPressPostn(i) == 0) NcellsInsideLimitPress = NcellsInsideLimitPress+1
    enddo
    write(*,*)CellPressPostn(1:imax),"CellPressPostn";write(*,*)NcellsInsideLimitPress,"NcellsInsideLimitPress"
    
    
    !numCurvtr=Hlf_Ncell ; write(*,*) numCurvtr,"numCurv" ; imax=numCurvtr
    
    NcellsInsideLimitCurvtr = 0
    
    if (NcellsInsideLimitPress .ge. (Hlf_Ncell-1)) then
       write(*,*) "Inside Curvtr Cnt",(Frame_NI-1)
       
       do i = 1,Hlf_Ncell
          
          curvtrNum = i      ; cellBfr = i
          if (i.ne.Hlf_Ncell)  cellAft = i+1 
          if (i == Hlf_Ncell)  cellAft = N_cell
          
          pressBfr           = k_area(cellBfr)*(A0(cellBfr)-A(cellBfr))
          pressAft           = k_area(cellAft)*(A0(cellAft)-A(cellAft))
          ExprmntlCurvtrSign = CurvtrSign_LtPCS(CMv,curvtrNum)
          
          call find_ifCurvtrInthe_Rnge(curvtrNum,pressAft,pressBfr,&
               ExprmntlCurvtrSign,CurrntCurvtrSign,CntCurvtr)
          if (CntCurvtr==1) NcellsInsideLimitCurvtr = NcellsInsideLimitCurvtr+1
          
       enddo
    endif
    
  end subroutine get_cntCell_BasedOn_PressOrCurvtrOrientn
  
  
  subroutine find_ifPressureInthe_Rnge(PressVal,LowRngePressVal,HighRngePressVal,CellPressPostnVal)
    implicit none
    real*8, intent(in)  :: PressVal,LowRngePressVal,HighRngePressVal
    real*8, intent(out) :: CellPressPostnVal
    
    if ((PressVal.gt.LowRngePressVal) .and. (PressVal.le.HighRngePressVal)) then
       CellPressPostnVal =  0
    elseif (PressVal.gt.HighRngePressVal) then
       CellPressPostnVal = +1
    elseif (PressVal.le.LowRngePressVal) then
       CellPressPostnVal = -1
    endif
    
  end subroutine find_ifPressureInthe_Rnge
  
  subroutine find_ifCurvtrInthe_Rnge(curvtrNum,pressAft,pressBfr,&
       ExprmntlCurvtrSign,CurrntCurvtrSign,CntCurvtr)
    implicit none
    integer, intent(in)  :: curvtrNum
    real*8 , intent(in)  :: pressAft,pressBfr
    integer, intent(in)  :: ExprmntlCurvtrSign
    integer, intent(out) :: CurrntCurvtrSign,CntCurvtr
    real*8               :: absPDval,prcntOfBndryPress
    
    write(*,*) curvtrNum,pressBfr,pressAft,abs(pressBfr-pressAft),"input for find_ifCurvtr_Rnge"
    
    open(unit=439,file='find_ifCurvtrInthe_Rnge.dat',position='append')
    if (curvtrNum==1) then
       write(439,*) "PRINTING FOR THE CURVATURE IN RANGE"
       write(439,*) " " ; write(439,*) " "
    endif
    
    write(439,*) "CurvtrNum = ",curvtrNum,"which Exp Curvtr Sign = ",ExprmntlCurvtrSign
    
    absPDval = abs(pressAft-pressBfr)
    write(439,*) "pressBfr = ",pressBfr,"pressAft = ",pressAft,"absPD = ",absPDval
    
    if (ExprmntlCurvtrSign ==  0) then
       
       !does not matter who is higher
       
       prcntOfBndryPress = (absPDval/BndryPress)*(100.0000d0)
       
       write(439,*) "is pressBfr ~~ pressAft ?"
       write(439,*) "prcntOfBP = ",prcntOfBndryPress,"TOL = ",TOL_PDwrtBP
       
       if (prcntOfBndryPress .gt. TOL_PDwrtBP) then ! Pressure Diff with respect to Boundary Press
          CntCurvtr = 0
       elseif (prcntOfBndryPress .le. TOL_PDwrtBP) then
          CntCurvtr = 1
       endif
       
       write(439,*) "CntCurvtr = ", CntCurvtr,"and CurvtrNum = ",curvtrNum
       
       if (CntCurvtr == 0) then
          if (pressBfr.gt.pressAft) CurrntCurvtrSign = +1
          if (pressBfr.le.pressAft) CurrntCurvtrSign = -1
       elseif (CntCurvtr == 1) then
          CurrntCurvtrSign = 0
       endif
       
       write(439,*) "CurrntCurvtrSign = ", CurrntCurvtrSign," and CurvtrNum = ", curvtrNum
       
    elseif (ExprmntlCurvtrSign == +1) then ! --- ) ---
       
       
       if (pressBfr .gt. pressAft) then
          
          prcntOfBndryPress = (absPDval/BndryPress)*(100.0000d0)
          
          write(439,*) "pressBfr > pressAft, however should it be counted ?"
          write(439,*) "prcntOfBP = ",prcntOfBndryPress,"TOL = ",TOL_PDwrtBP
          
          if (prcntOfBndryPress .gt. TOL_PDwrtBP) then ! Pressure Diff with respect to Boundary Press
             CntCurvtr = 1 
          elseif (prcntOfBndryPress .le. TOL_PDwrtBP) then
             CntCurvtr = 0
          endif
          
          write(439,*) "CntCurvtr = ", CntCurvtr,"and CurvtrNum = ",curvtrNum
          
          CurrntCurvtrSign = +1
          
       elseif (pressBfr .le. pressAft) then
          write(439,*) "pressBfr < pressAft, so it should NOT be counted"
          CntCurvtr = -1
          write(439,*) "CntCurvtr = ", CntCurvtr,"and CurvtrNum = ",curvtrNum
          
          CurrntCurvtrSign = -1
       endif
       
       write(439,*) "CurrntCurvtrSign = ", CurrntCurvtrSign," and CurvtrNum = ", curvtrNum
       
    elseif (ExprmntlCurvtrSign == -1) then ! --- ( ---
       
       if (pressBfr .le. pressAft) then
          
          prcntOfBndryPress = (absPDval/BndryPress)*(100.0000d0)
          
          write(439,*) "pressBfr < pressAft, however should it be counted ?"
          write(439,*) "prcntOfBP = ",prcntOfBndryPress,"TOL = ",TOL_PDwrtBP
          
          if (prcntOfBndryPress .gt. TOL_PDwrtBP) then ! Pressure Diff with respect to Boundary Press
             CntCurvtr = 1 
          elseif (prcntOfBndryPress .le. TOL_PDwrtBP) then
             CntCurvtr = 0
          endif
          
          write(439,*) "CntCurvtr = ", CntCurvtr,"and CurvtrNum = ",curvtrNum
          
          CurrntCurvtrSign = -1
          
       elseif (pressBfr .gt. pressAft) then
          write(439,*) "pressBfr > pressAft, so it should NOT be counted"
          CntCurvtr = -1
          write(439,*) "CntCurvtr = ", CntCurvtr,"and CurvtrNum = ",curvtrNum
          
          CurrntCurvtrSign = +1
       endif
       
       write(439,*) "CurrntCurvtrSign = ", CurrntCurvtrSign," and CurvtrNum = ", curvtrNum
       
    endif
    
    close(439)
    
    write(*,*) ExprmntlCurvtrSign,CurrntCurvtrSign,CntCurvtr,"outputs for find_ifCurvtrInthe_Rnge"
    
  end subroutine find_ifCurvtrInthe_Rnge
  
  subroutine get_cntCurvtrInsideLimit_Gnrl(LowRngeCurvtr,HighRngeCurvtr,CellCurvtrPostn,NcellsInsideLimitPress)
    implicit none
    real*8,  intent(in)    :: LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    integer, intent(inout) :: CellCurvtrPostn(1:Hlf_Ncell)
    integer, intent(inout) :: NcellsInsideLimitPress
    
    real*8  :: CurvtrValue(1:Hlf_Ncell),CurvtrVal
    real*8  :: LowRngeCurvtrVal,HighRngeCurvtrVal
    integer :: i,j,imax,jmax
    integer :: BsORLt
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurvtrValue)
    
    imax = Hlf_Ncell
    
    do i = 1,imax
       
       CurvtrVal         = CurvtrValue(i)
       LowRngeCurvtrVal  = LowRngeCurvtr(i)
       HighRngeCurvtrVal = HighRngeCurvtr(i)
       
       if ((CurvtrVal.gt.LowRngeCurvtrVal) .and. (CurvtrVal.le.HighRngeCurvtrVal)) then
          CellCurvtrPostn(i) =  0
       elseif (CurvtrVal.gt.HighRngeCurvtrVal) then
          CellCurvtrPostn(i) = +1
       elseif (CurvtrVal.le.LowRngeCurvtrVal) then
          CellCurvtrPostn(i) = -1
       endif
       
       if (CellCurvtrPostn(i) == 0) NcellsInsideLimitPress = NcellsInsideLimitPress+1
       
    enddo
    
    write(*,*) CellCurvtrPostn(1:imax),"CellCurvtrPostn"
    write(*,*) NcellsInsideLimitPress,"NcellsInsideLimitPress"
    
  end subroutine get_cntCurvtrInsideLimit_Gnrl
  
  subroutine get_cntCurvtrInsideLimit_CurvtrSpcific(numCurvtr,Curvtrs,LowRngeCurvtr,HighRngeCurvtr,&
       spcfCellCurvtrPostn,NspcfCellsInsideLimit)
    implicit none
    integer, intent(in)    :: numCurvtr
    integer, intent(in)    :: Curvtrs(1:numCurvtr)
    real*8,  intent(in)    :: LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    integer, intent(inout) :: SpcfCellCurvtrPostn(1:numCurvtr),NspcfCellsInsideLimit
    
    real*8  :: CurvtrValue(1:Hlf_Ncell),CurvtrVal
    real*8  :: LowRngeCurvtrVal,HighRngeCurvtrVal
    integer :: i,j,imax,jmax,BsORLt,CurvtrNumber
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurvtrValue)
    
    imax = numCurvtr
    
    do i = 1,imax
       
       CurvtrNumber      = Curvtrs(i)
       CurvtrVal         = CurvtrValue(CurvtrNumber)
       LowRngeCurvtrVal  = LowRngeCurvtr(CurvtrNumber)
       HighRngeCurvtrVal = HighRngeCurvtr(CurvtrNumber)
       
       if ((CurvtrVal.gt.LowRngeCurvtrVal) .and. (CurvtrVal.le.HighRngeCurvtrVal)) then
          SpcfCellCurvtrPostn(i) =  0
       elseif (CurvtrVal.gt.HighRngeCurvtrVal) then
          SpcfCellCurvtrPostn(i) = +1
       elseif (CurvtrVal.le.LowRngeCurvtrVal) then
          SpcfCellCurvtrPostn(i) = -1
       endif
       
       if (SpcfCellCurvtrPostn(i) == 0) NspcfCellsInsideLimit = NspcfCellsInsideLimit+1
       
    enddo
    
    write(*,*) SpcfCellCurvtrPostn(1:imax),"SpcfCellCurvtrPostn"
    write(*,*) NspcfCellsInsideLimit,"NspcfCellsInsideLimit"
    
  end subroutine get_cntCurvtrInsideLimit_CurvtrSpcific
  
  
  subroutine make_cellPress_inDefinedRnge(numOfCellsToAdjst)
    implicit none
    integer, intent(in) :: numOfCellsToAdjst
    
    real*8  :: AchievingPress(1:numOfCellsToAdjst)
    real*8  :: LowRngePress(1:numOfCellsToAdjst)
    real*8  :: HighRngePress(1:numOfCellsToAdjst)
    integer :: CellPressPostn(1:numOfCellsToAdjst)
    
    integer :: NcellsInsideLimitPress
    integer :: i,imax,j,jmax,cntLp
    integer :: cellL,cellR
    
    
    call get_Achieving_HighRnge_LowRngePress(numOfCellsToAdjst,&
         AchievingPress,LowRngePress,HighRngePress)
    
    cntLp = 1
    
    do
       
       NcellsInsideLimitPress=0 ; CellPressPostn=-10
       call get_cntCellInsideLimit_Gnrl(numOfCellsToAdjst,LowRngePress,HighRngePress,&
            CellPressPostn,NcellsInsideLimitPress)
       
       if (NcellsInsideLimitPress == numOfCellsToAdjst) exit
       
       imax = numOfCellsToAdjst
       
       do i = 1,imax
          
          if (i.le.Hlf_Ncell) then
             cellL = i
             cellR = cellL+Hlf_Ncell
          elseif (i.gt.Hlf_Ncell) then
             cellL = (2*Hlf_Ncell)+1
          endif
          
          if (CellPressPostn(i) == 0 ) then
             continue
          elseif (CellPressPostn(i) == +1) then
             
             if (i.le.Hlf_Ncell) then
                A0(cellL) = 0.9975d0*A0(cellL) ! 0.995
                A0(cellR) = A0(cellL)
             elseif (i.gt.Hlf_Ncell) then
                A0(cellL) = 0.9975d0*A0(cellL)
             endif
                
          elseif (CellPressPostn(i) == -1) then
             
             if (i.le.Hlf_Ncell) then
                A0(cellL) = 1.0025d0*A0(cellL) ! 1.005
                A0(cellR) = A0(cellL)
             elseif (i.gt.Hlf_Ncell) then
                A0(cellL) = 1.0025d0*A0(cellL)
             endif
             
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
    enddo
    
  end subroutine make_cellPress_inDefinedRnge
  
  
  
  subroutine make_cellPress_inDefinedRnge_generalized(CMv,PCSorICS,numOfCellsToAdjst,PressChoice,CntChoice,TOL_prcnt)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    integer, intent(in) :: numOfCellsToAdjst
    integer, intent(in) :: PressChoice,CntChoice
    real*8 , intent(in) :: TOL_prcnt
    
    real*8  :: AchievingPress(1:numOfCellsToAdjst)
    real*8  :: LowRngePress(1:numOfCellsToAdjst)
    real*8  :: HighRngePress(1:numOfCellsToAdjst)
    integer :: CellPressPostn(1:numOfCellsToAdjst)
    
    integer :: NcellsInsideLimitPress,NcellsInsideLimitCurvtr
    integer :: i,imax,j,jmax,cntLp
    integer :: cellL,cellR
    real*8  :: curve_LtH(1:Hlf_Ncell)
    integer :: BsORLt
    integer :: FrameNumSaved1,FrameNumSaved2
    
    if (PressChoice == 1) then
       call get_Achieving_HighRnge_LowRngePress_generalized(CMv,PCSorICS,numOfCellsToAdjst,&
            AchievingPress,LowRngePress,HighRngePress)  
       
    elseif (PressChoice == 2) then
       
       FrameNumSaved1 = Frame_NI-1 ; write(*,*) FrameNumSaved1,"FrmNumSaved1"
       call matchBiggstCurve_and_Achieving_HighLowRngePress(CMv,PCSorICS,AchievingPress,&
            TOL_prcnt,LowRngePress,HighRngePress)
       FrameNumSaved2 = Frame_NI-1     ; write(*,*) FrameNumSaved2, "FrmNumSaved2"
       Frame_NI       = FrameNumSaved1
       call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
       call Equilibrate_only_NI_model
       Frame_NI = FrameNumSaved2
       
    elseif (PressChoice == 3) then
       call rscale_and_guessd_AchiLowHighPress_basedCurvtrOrientn(CMv,PCSorICS,TOL_prcnt,&
            AchievingPress,LowRngePress,HighRngePress)
    endif
    
    cntLp = 1
    
    do
       
       NcellsInsideLimitPress=0 ; CellPressPostn=-10
       
       if (CntChoice == 1) then
          call get_cntCellInsideLimit_Gnrl(numOfCellsToAdjst,LowRngePress,HighRngePress,&
               CellPressPostn,NcellsInsideLimitPress)
       elseif (CntChoice == 2) then
          call get_cntCell_BasedOn_PressOrCurvtrOrientn(CMv,numOfCellsToAdjst,LowRngePress,HighRngePress,&
               CellPressPostn,NcellsInsideLimitPress,NcellsInsideLimitCurvtr)
       endif
       
       if (NcellsInsideLimitPress == numOfCellsToAdjst) exit
       if (NcellsInsideLimitPress == (Hlf_Ncell))       exit
       
       imax = numOfCellsToAdjst
       
       do i = 1,imax
          
          if (i.le.Hlf_Ncell) then
             cellL = i
             cellR = cellL+Hlf_Ncell
          elseif (i.gt.Hlf_Ncell) then
             cellL = (2*Hlf_Ncell)+1
          endif
          
          if (CellPressPostn(i) == 0 ) then
             continue
          elseif (CellPressPostn(i) == +1) then
             
             if (i.le.Hlf_Ncell) then
                A0(cellL) = 0.9975d0*A0(cellL) ! 0.995,0.9975d0,0.99875
                A0(cellR) = A0(cellL)
             elseif (i.gt.Hlf_Ncell) then  
                A0(cellL) = 0.9975d0*A0(cellL)
             endif
             
          elseif (CellPressPostn(i) == -1) then
             
             if (i.le.Hlf_Ncell) then
                A0(cellL) = 1.0025d0*A0(cellL) ! 1.005,1.0025d0,1.00125
                A0(cellR) = A0(cellL)
             elseif (i.gt.Hlf_Ncell) then
                A0(cellL) = 1.0025d0*A0(cellL)
             endif
             
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
    enddo
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    do i = 1,Hlf_Ncell
       write(*,*) Expected_CHLtPCS(CMv,i),curve_LtH(i),&
            (Expected_CHLtPCS(CMv,i)-curve_LtH(i)),i,"curveCheck"
    enddo
    
    write(*,*) "curveCheck DONE"
    !stop 'curveCheck DONE'
    
  end subroutine make_cellPress_inDefinedRnge_generalized
  
  
  
  subroutine make_curvtr_inDefinedRnge_generalized(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    real*8  :: AchievingCurvtr(1:Hlf_Ncell)
    real*8  :: LowRngeCurvtr(1:Hlf_Ncell)
    real*8  :: HighRngeCurvtr(1:Hlf_Ncell)
    integer :: CellCurvtrPostn(1:Hlf_Ncell)
    
    integer :: NcellsInsideLimitPress
    integer :: i,imax,j,jmax,cntLp
    integer :: cellLb,cellRb,cellLa,cellRa
    real*8  :: curve_LtH(1:Hlf_Ncell)
    integer :: BsORLt
    real*8  :: pressBfrCell,pressAftCell
    real*8  :: f1,f2,f3,f4,fctrV
    
    
    call get_Achieving_HighLowRnge_Curvtr(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,HighRngeCurvtr)

    if (CMv==1) then
       
       !call get_first_bndryPress_and_curvtrOrientatn_right(CMv,PCSorICS)
       
       call altrPresUsingTolrnc_generalized(CMv,PCSorICS) ; write(*,*) "Pos2 FrmNI=",(Frame_NI-1)
       
    elseif (CMv==2) then
       
       !call get_first_bndryPress_and_curvtrOrientatn_right(CMv,PCSorICS)
       
       call altrPresUsingTolrnc_generalized(CMv,PCSorICS) ; write(*,*) "Pos2 FrmNI=",(Frame_NI-1)
       
    elseif (CMv==3) then
       
       !call get_first_bndryPress_and_curvtrOrientatn_right(CMv,PCSorICS)
       
       call altrPresUsingTolrnc_generalized(CMv,PCSorICS)
       
    elseif (CMv==4) then
       
       !call get_first_bndryPress_and_curvtrOrientatn_right(CMv,PCSorICS)
       
       call altrPresUsingTolrnc_generalized(CMv,PCSorICS) ; write(*,*) "Frame_NI = ",(Frame_NI-1)
       
       
    endif
    
    
    f1 = 0.005000d0 ; f2 = 0.002500d0 ; f3 = 0.001250d0 ; f4 = 0.000625d0
    
    cntLp = 1
    
    do
       
       open(unit=784,file='CurvtrInRnge.dat',position='append')
       
       write(784,*) "Logs For Cnt =",cntLp ; write(784,*) " "
       NcellsInsideLimitPress=0 ; CellCurvtrPostn=-10
       
       call get_cntCurvtrInsideLimit_Gnrl(LowRngeCurvtr,HighRngeCurvtr,CellCurvtrPostn,NcellsInsideLimitPress)
       BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
       
       do i = 1,Hlf_Ncell
          write(784,*) LowRngeCurvtr(i),AchievingCurvtr(i),HighRngeCurvtr(i),Curve_LtH(i),CellCurvtrPostn(i),&
               CurvtrSign_LtPCS(CMv,i),i,"L/H/C"
       enddo
       write(784,*) " " ; write(784,*) NcellsInsideLimitPress,"NcellsInsideLimitPress" ; write(784,*) " "
       
       if (NcellsInsideLimitPress  ==  (Hlf_Ncell-1)) fctrV = f4/4.0000000d0
       if (NcellsInsideLimitPress .lt. (Hlf_Ncell-1)) fctrV = f2
       if (NcellsInsideLimitPress  ==   Hlf_Ncell)    exit
       !if (NcellsInsideLimitPress ==  (Hlf_Ncell-2))    exit
       
       imax = Hlf_Ncell
       
       do i = 1,imax
          
          cellLb = i ; cellRb = cellLb+Hlf_Ncell
          
          if (i.ne.imax) then
             cellLa = i+1
             cellRa = cellLa+Hlf_Ncell
          elseif (i == imax) then
             cellLa = N_cell
             cellRa = N_cell
          endif
          
          if (CurvtrSign_LtPCS(CMv,i) == +1) then ! ---)---
             
             if (CellCurvtrPostn(i) == 0) then
                
                write(784,*) "Entering Block 11 for curveNum =",i
                continue
                
             elseif (CellCurvtrPostn(i) == +1) then
                
                write(784,*) "Entering Block 12 for curveNum =",i
                
                A0(cellLb) = (1.0000d0-fctrV) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
                
                !A0(cellLa) = (1.0000d0+fctrV) * A0(cellLa)
                
                !if (i.ne.imax) A0(cellRa) = A0(cellLa)
                !if (i == imax) continue ! to not changing A0(N_cell) twice
                
             elseif (CellCurvtrPostn(i) == -1) then
                
                write(784,*) "Entering Block 13 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
                
                !A0(cellLa) = (1.000000d0-fctrV) * A0(cellLa)
                
                !if (i.ne.imax) A0(cellRa) = A0(cellLa) 
                !if (i == imax) continue ! to not changing A0(N_cell) twice
                
             endif
             
          elseif (CurvtrSign_LtPCS(CMv,i) == -1) then ! ---(---
             
             if (CellCurvtrPostn(i) == 0) then
                
                write(784,*) "Entering Block 21 for curveNum =",i
                continue
                   
             elseif (CellCurvtrPostn(i) == +1) then
                
                write(784,*) "Entering Block 22 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV) * (A0(cellLb))         ! 1.005,1.0025d0,1.00125
                A0(cellRb) = A0(cellLb)
                
                !A0(cellLa) = (1.000000d0-fctrV) * (A0(cellLa))
                
                !if (i.ne.imax) A0(cellRa) = A0(cellLa) 
                !if (i == imax) continue ! to not changing A0(N_cell) twice
                
                
             elseif (CellCurvtrPostn(i) == -1) then
                
                write(784,*) "Entering Block 23 for curveNum =",i
                
                A0(cellLb) = (1.000000d0-fctrV) * A0(cellLb)         ! 0.995,0.9975d0,0.99875   
                A0(cellRb) = A0(cellLb)
                
                !A0(cellLa) = (1.000000d0+fctrV) * A0(cellLa)
                
                !if (i.ne.imax) A0(cellRa) = A0(cellLa) 
                !if (i == imax) continue ! to not changing A0(N_cell) twice
                
             endif
                
             
          elseif (CurvtrSign_LtPCS(CMv,i) == 0) then ! ---|---
             
             if (CellCurvtrPostn(i) == 0) then
                   
                write(784,*) "Entering Block 31 for curveNum =",i
                continue
                
             elseif (CellCurvtrPostn(i) == +1) then
                
                write(784,*) "Entering Block 32 for curveNum =",i
                   
                pressBfrCell = k_area(cellLb)*(A0(cellLb)-A(cellLb))
                pressAftCell = k_area(cellLa)*(A0(cellLa)-A(cellLa))
                
                if (pressBfrCell .gt. pressAftCell) then
                   
                   A0(cellLb) = (1.000000d0-fctrV)*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)
                   
                   !A0(cellLa) = (1.000000d0+fctrV) * A0(cellLa)
                   
                   !if (i.ne.imax) A0(cellRa) = A0(cellLa) 
                   !if (i == imax) continue ! to not changing A0(N_cell) twice
                   
                elseif (pressBfrCell .le. pressAftCell) then
                   
                   A0(cellLb) = (1.000000d0+fctrV)*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)
                   
                   !A0(cellLa) = (1.000000d0-fctrV)*A0(cellLa)
                   
                   !if (i.ne.imax) A0(cellRa) = A0(cellLa) 
                   !if (i == imax) continue ! to not changing A0(N_cell) twice
                   
                endif
                
             elseif (CellCurvtrPostn(i) == -1) then
                
                write(784,*) "Entering Block 33 for curveNum =",i
                write(*,*) "MUST NOT ENTER Block 33"
                stop
                
             endif
             
          endif
          
          !endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
       !if (cntLp==600) exit
       
       close(784)
       
    enddo
    
    !close(784)
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    do i = 1,Hlf_Ncell
       write(*,*) Expected_CHLtPCS(CMv,i),curve_LtH(i),(Expected_CHLtPCS(CMv,i)-curve_LtH(i)),i,"curveCheck"
    enddo
    
  end subroutine make_curvtr_inDefinedRnge_generalized
  
  
  subroutine make_curvtr_inDefinedRnge_gen_with_correct_curvtrOrientn(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    real*8  :: AchievingCurvtr(1:Hlf_Ncell)
    real*8  :: LowRngeCurvtr(1:Hlf_Ncell)
    real*8  :: HighRngeCurvtr(1:Hlf_Ncell)
    real*8  :: ExpctdCurvtr(1:Hlf_Ncell)
    integer :: CurvtrSignLt(1:Hlf_Ncell)
    integer :: CellCurvtrPostn(1:Hlf_Ncell)
    
    integer :: NcellsInsideLimitPress
    integer :: i,imax,j,jmax,cntLp
    integer :: cellLb,cellRb,cellLa,cellRa
    real*8  :: curve_LtH(1:Hlf_Ncell)
    integer :: BsORLt
    real*8  :: pressBfrCell,pressAftCell
    real*8  :: f1,f2,f3,f4,fctrV
    
    call get_Achieving_HighLowRnge_CurvtrRead(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,&
         HighRngeCurvtr,ExpctdCurvtr,CurvtrSignLt)
    
    f1 = 0.005000d0 ; f2 = 0.002500d0 ; f3 = 0.001250d0 ; f4 = 0.000625d0
    
    cntLp = 1
    
    do
       
       open(unit=784,file='CurvtrInRngeV2.dat',position='append')
       
       write(784,*) "Logs For Cnt =",cntLp ; write(784,*) " "
       NcellsInsideLimitPress=0 ; CellCurvtrPostn=-10
       
       call get_cntCurvtrInsideLimit_Gnrl(LowRngeCurvtr,HighRngeCurvtr,CellCurvtrPostn,&
            NcellsInsideLimitPress)
       BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
       
       do i = 1,Hlf_Ncell
          write(784,*) LowRngeCurvtr(i),AchievingCurvtr(i),HighRngeCurvtr(i),Curve_LtH(i),&
               CellCurvtrPostn(i),CurvtrSignLt(i),i,"L/H/C"
       enddo
       write(784,*)" ";write(784,*) NcellsInsideLimitPress,"NcellsInsideLimitPress";write(784,*)" "
       
       if (NcellsInsideLimitPress  ==  (Hlf_Ncell-1)) fctrV = f4/4.0000000d0
       if (NcellsInsideLimitPress .lt. (Hlf_Ncell-1)) fctrV = f2
       if (NcellsInsideLimitPress  ==   Hlf_Ncell)    exit
       
       imax = Hlf_Ncell
       
       do i = 1,imax
          
          cellLb = i ; cellRb = cellLb+Hlf_Ncell
          
          if (i.ne.imax) then
             cellLa = i+1
             cellRa = cellLa+Hlf_Ncell
          elseif (i == imax) then
             cellLa = N_cell
             cellRa = N_cell
          endif
          
          if (CurvtrSignLt(i) == +1) then ! ---)---
             
             if (CellCurvtrPostn(i) == 0) then
                write(784,*) "Entering Block 11 for curveNum =",i
                continue      
             elseif (CellCurvtrPostn(i) == +1) then
                write(784,*) "Entering Block 12 for curveNum =",i
                
                A0(cellLb) = (1.0000d0-fctrV) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
             elseif (CellCurvtrPostn(i) == -1) then
                write(784,*) "Entering Block 13 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
             endif
             
          elseif (CurvtrSignLt(i) == -1) then ! ---(---
             
             if (CellCurvtrPostn(i) == 0) then
                write(784,*) "Entering Block 21 for curveNum =",i
                continue
             elseif (CellCurvtrPostn(i) == +1) then
                write(784,*) "Entering Block 22 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV) * (A0(cellLb))         ! 1.005,1.0025d0,1.00125
                A0(cellRb) = A0(cellLb)
             elseif (CellCurvtrPostn(i) == -1) then
                write(784,*) "Entering Block 23 for curveNum =",i
                
                A0(cellLb) = (1.000000d0-fctrV) * A0(cellLb)         ! 0.995,0.9975d0,0.99875   
                A0(cellRb) = A0(cellLb)
             endif
                
             
          elseif (CurvtrSignLt(i) == 0) then ! ---|---
             
             if (CellCurvtrPostn(i) == 0) then
                write(784,*) "Entering Block 31 for curveNum =",i
                continue
             elseif (CellCurvtrPostn(i) == +1) then
                write(784,*) "Entering Block 32 for curveNum =",i
                
                pressBfrCell = k_area(cellLb)*(A0(cellLb)-A(cellLb))
                pressAftCell = k_area(cellLa)*(A0(cellLa)-A(cellLa))
                
                if (pressBfrCell .gt. pressAftCell) then
                   A0(cellLb) = (1.000000d0-fctrV)*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)   
                elseif (pressBfrCell .le. pressAftCell) then
                   A0(cellLb) = (1.000000d0+fctrV)*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)
                endif
                
             elseif (CellCurvtrPostn(i) == -1) then
                write(784,*) "Entering Block 33 for curveNum =",i
                write(*,*) "MUST NOT ENTER Block 33"
                stop
             endif
             
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
       close(784)
       
    enddo
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    do i = 1,Hlf_Ncell
       write(*,*) ExpctdCurvtr(i),curve_LtH(i),(ExpctdCurvtr(i)-curve_LtH(i)),i,"curveCheck"
    enddo
    
  end subroutine make_curvtr_inDefinedRnge_gen_with_correct_curvtrOrientn
  
  subroutine make_curvtr_inDefinedRnge_generalized_varSpeed(CMv,PCSorICS,numCurvtr,Curvtrs)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    integer, intent(in) :: numCurvtr
    integer, intent(in) :: Curvtrs(1:numCurvtr)
    
    real*8  :: AchievingCurvtr(1:Hlf_Ncell),LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    integer :: CellCurvtrPostn(1:Hlf_Ncell),SpcfCellCurvtrPostn(1:numCurvtr)
    
    integer :: NcellsInsideLimitPress,NspcfCellsInsideLimit
    integer :: i,imax,j,jmax,cntLp
    integer :: cellLb,cellRb,cellLa,cellRa
    real*8  :: curve_LtH(1:Hlf_Ncell)
    integer :: BsORLt,CurvtrNumber
    real*8  :: pressBfrCell,pressAftCell
    real*8  :: f1,f2,f3,f4
    real*8  :: fctrV(1:Hlf_Ncell)
    
    call get_Achieving_HighLowRnge_Curvtr(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,HighRngeCurvtr)
    
    if (CMv==2) continue
    if (CMv==3) continue
    
    f1 = 0.005000d0 ; f2 = 0.002500d0 ; f3 = 0.001250d0 ; f4 = 0.000625d0
    fctrV(1:Hlf_Ncell)   = 0.000000d0
    fctrV(Curvtrs(1))    = f2         ; fctrV(Curvtrs(2)) = f1 
    
    cntLp = 1
    
    do
       open(unit=787,file='CurvtrInRngeVarSpeed.dat',position='append')
       write(787,*) "Logs For Cnt =",cntLp ; write(787,*) " "
       
       NcellsInsideLimitPress     = 0 ; CellCurvtrPostn     = -10
       NspcfCellsInsideLimit = 0 ; SpcfCellCurvtrPostn = -10
       
       call get_cntCurvtrInsideLimit_Gnrl(LowRngeCurvtr,HighRngeCurvtr,CellCurvtrPostn,NcellsInsideLimitPress)
       call get_cntCurvtrInsideLimit_CurvtrSpcific(numCurvtr,Curvtrs,LowRngeCurvtr,HighRngeCurvtr,&
            spcfCellCurvtrPostn,NspcfCellsInsideLimit)
       
       BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
       
       write(787,*) " " ; write(787,*) NcellsInsideLimitPress,"NcellsInsideLimitPress"         ; write(787,*) " "
       write(787,*) " " ; write(787,*) NspcfCellsInsideLimit,"NspcfCellsInsideLimit" ; write(787,*) " "
       
       if (NspcfCellsInsideLimit  == numCurvtr) then
          write(*,*) NspcfCellsInsideLimit,"NspcfCellsInsideLimit at exit"
          exit
       endif
       
       imax = numCurvtr
       
       write(787,*) "Starting of cntLp =",cntLp
       
       do i = 1,imax
          
          CurvtrNumber = Curvtrs(i)   ; write(*,*) CurvtrNumber,i,"CurvtrNumber-i"
          cellLb       = CurvtrNumber ; cellRb = cellLb+Hlf_Ncell
          
          write(787,*) CurvtrNumber,i,"CurvtrNumber-i"
          write(787,*) cellLb,cellRb,i,"cellLb-cellRb"
          write(787,*) CurvtrSign_LtPCS(CMv,CurvtrNumber),"Curvtr_LtPCS"
          write(787,*) LowRngeCurvtr(CurvtrNumber),HighRngeCurvtr(CurvtrNumber),curve_LtH(CurvtrNumber),"LHC"
          write(787,*) CellCurvtrPostn(CurvtrNumber),"CellCurvtrPostn"
          write(787,*) A0(cellLb),A0(cellRb),i,"cellLb-cellRb"
          
          cellLa = cellLb+1
          cellRa = cellLa+Hlf_Ncell
          
          if (i == Hlf_Ncell) then
             cellLa = N_cell
             cellRa = N_cell
          endif
          
          if (CurvtrSign_LtPCS(CMv,CurvtrNumber) == +1) then ! ---)---
             
             if (CellCurvtrPostn(CurvtrNumber) == 0) then
                
                write(787,*) "Entering BlockNum 11 for curveNum =",i
                continue
                
             elseif (CellCurvtrPostn(CurvtrNumber) == +1) then
                
                write(787,*) "Entering BlockNum 12 for curveNum =",i
                
                A0(cellLb) = (1.0000d0-fctrV(CurvtrNumber)) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
                
             elseif (CellCurvtrPostn(CurvtrNumber) == -1) then
                   
                write(787,*) "Entering BlockNum 13 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV(CurvtrNumber)) * A0(cellLb)
                A0(cellRb) = A0(cellLb)
                
             endif
             
          elseif (CurvtrSign_LtPCS(CMv,CurvtrNumber) == -1) then ! ---(---
             
             if (CellCurvtrPostn(CurvtrNumber) == 0) then
                
                write(787,*) "Entering BlockNum 21 for curveNum =",i
                continue
                   
             elseif (CellCurvtrPostn(CurvtrNumber) == +1) then
                
                write(787,*) "Entering BlockNum 22 for curveNum =",i
                
                A0(cellLb) = (1.000000d0+fctrV(CurvtrNumber)) * (A0(cellLb))         ! 1.005,1.0025d0,1.00125
                A0(cellRb) = A0(cellLb)
                
             elseif (CellCurvtrPostn(CurvtrNumber) == -1) then
                
                write(787,*) "Entering BlockNum 23 for curveNum =",i
                
                A0(cellLb) = (1.000000d0-fctrV(CurvtrNumber)) * A0(cellLb)         ! 0.995,0.9975d0,0.99875   
                A0(cellRb) = A0(cellLb)
                
             endif
             
             
          elseif (CurvtrSign_LtPCS(CMv,CurvtrNumber) == 0) then ! ---|---
             
             if (CellCurvtrPostn(CurvtrNumber) == 0) then
                
                write(787,*) "Entering BlockNum 31 for curveNum =",i
                continue
                
             elseif (CellCurvtrPostn(CurvtrNumber) == +1) then
                
                write(787,*) "Entering BlockNum 32 for curveNum =",i
                
                pressBfrCell = k_area(cellLb)*(A0(cellLb)-A(cellLb))
                pressAftCell = k_area(cellLa)*(A0(cellLa)-A(cellLa))
                
                if (pressBfrCell .gt. pressAftCell) then
                   
                   A0(cellLb) = (1.000000d0-fctrV(CurvtrNumber))*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)
                   
                elseif (pressBfrCell .le. pressAftCell) then
                   
                   A0(cellLb) = (1.000000d0+fctrV(CurvtrNumber))*A0(cellLb)       
                   A0(cellRb) = A0(cellLb)
                   
                endif
                
             elseif (CellCurvtrPostn(CurvtrNumber) == -1) then
                
                write(787,*) "Entering BlockNum 33 for curveNum =",i
                write(*,*) "MUST NOT ENTER BlockNum 33"
                stop
                
             endif
                
          endif
          
       enddo
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
       if (cntLp==20) exit
       
       close(787)
       
    enddo
    
    !close(787)
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    do i = 1,Hlf_Ncell
       write(*,*) Expected_CHLtPCS(CMv,i),curve_LtH(i),(Expected_CHLtPCS(CMv,i)-curve_LtH(i)),i,"curveCheck"
    enddo
    
  end subroutine make_curvtr_inDefinedRnge_generalized_varSpeed
  
  
  ! subroutine hold_NcellsLimit(CMv,CellCurvtrPostn,fctrV)
  !   implicit none
  !   integer, intent(in) :: CMv,CellCurvtrPostn(1:Hlf_Ncell)
  !   real*8 , intent(in) :: fctrV
    
  !   integer             :: i,j,imax,jmax
  !   integer             :: cellL,cellR
  !   real*8              :: pressBfrCell,pressAftCell
    
  !   do i = 1,imax
       
  !      cellL = i
  !      cellR = cellL+Hlf_Ncell
       
  !      if (CurvtrSign_LtPCS(CMv,i) == +1) then ! ---)---
          
  !         if (CellCurvtrPostn(i) == 0) then
  !            continue
  !         elseif (CellCurvtrPostn(i) == +1) then
  !            A0(cellL) = (1.000000d0-fctrV) * A0(cellL) ; A0(cellR) = A0(cellL)
  !         elseif (CellCurvtrPostn(i) == -1) then
  !            A0(cellL) = (1.000000d0+fctrV) * A0(cellL) ; A0(cellR) = A0(cellL)
  !         endif
          
  !      elseif (CurvtrSign_LtPCS(CMv,i) == -1) then ! ---(---
          
  !         if (CellCurvtrPostn(i) == 0) then
  !            continue
  !         elseif (CellCurvtrPostn(i) == +1) then
  !            A0(cellL) = (1.000000d0+fctrV) * (A0(cellL)) ; A0(cellR) = A0(cellL)
  !         elseif (CellCurvtrPostn(i) == -1) then
  !            A0(cellL) = (1.000000d0-fctrV) * (A0(cellL)) ; A0(cellR) = A0(cellL)
  !         endif
          
          
  !      elseif (CurvtrSign_LtPCS(CMv,i) == 0) then ! ---|---
          
  !         if (CellCurvtrPostn(i) == 0) then
  !            continue
  !         elseif (CellCurvtrPostn(i) == +1) then
  !            pressBfrCell = k_area(i)*(A0(i)-A(i))
             
  !            if (i.ne.Hlf_Ncell) pressAftCell = k_area(i+1)    * (A0(i+1)    - A(i+1))
  !            if (i == Hlf_Ncell) pressAftCell = k_area(N_cell) * (A0(N_cell) - A(N_cell))
             
  !            if (pressBfrCell .gt. pressAftCell) then   
  !               A0(cellL) = (1.000000d0-fctrV)*A0(cellL) ; A0(cellR) = A0(cellL)
  !            elseif (pressBfrCell .le. pressAftCell) then
  !               A0(cellL) = (1.000000d0+fctrV)*A0(cellL) ; A0(cellR) = A0(cellL)
  !            endif
                
  !         elseif (CellCurvtrPostn(i) == -1) then
  !            stop 'MUST NOT ENTER Block 33'
  !         endif
          
  !      endif
       
  !   enddo
    
  !   call Equilibrate_only_NI_model

    
  !   call get_cntCurvtrInsideLimit_Gnrl(LowRngeCurvtr,HighRngeCurvtr,CellCurvtrPostn,NcellsInsideLimitPress)
  !   BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    
  ! end subroutine hold_NcellsLimit
  
  
  
  

  
  
  
  subroutine print_the_curvtr_params(CMv,LowRngeCurvtr,AchievingCurvtr,HighRngeCurvtr)
    implicit none
    integer, intent(in) :: CMv
    real*8,  intent(in) :: LowRngeCurvtr(1:Hlf_Ncell),AchievingCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    
    integer :: i,j,BsORLt
    real*8  :: Curve_LtH(1:Hlf_Ncell)
    
    BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,Curve_LtH) 
    
    do i = 1,Hlf_Ncell
       write(*,*) LowRngeCurvtr(i),AchievingCurvtr(i),HighRngeCurvtr(i),Curve_LtH(i),CurvtrSign_LtPCS(CMv,i),i,"CURVE"
    enddo
    
  end subroutine print_the_curvtr_params
  
  
  subroutine get_Achieving_HighRnge_LowRngePress(numOfCellsToAdjst,&
       AchievingPress,LowRngePress,HighRngePress) ! %%%%%%%%%%%%%%%%%%%%%%% CM2PCS %%%%%%%%%%%
    implicit none
    integer, intent(in)  :: numOfCellsToAdjst
    real*8,  intent(out) :: AchievingPress(1:numOfCellsToAdjst)
    real*8,  intent(out) :: LowRngePress(1:numOfCellsToAdjst)
    real*8,  intent(out) :: HighRngePress(1:numOfCellsToAdjst)
    
    integer :: i,imax,j,jmax
    integer :: lim1,lim2,lim3
    integer :: N_IncmingCellEquivPress,N_LowPressThnIncmCell
    
    
    imax = numOfCellsToAdjst
    
    N_IncmingCellEquivPress = 6 ; lim1 = N_IncmingCellEquivPress 
    N_LowPressThnIncmCell   = 9 ; lim2 = N_LowPressThnIncmCell 
    
    
    do i = 1,imax
       
       if (i.le.lim1) then   
          AchievingPress(i) = 0.1000d0
       elseif ((i.gt.lim1) .and. (i.le.lim2)) then
          AchievingPress(i) = 0.0930d0
       endif
       
       LowRngePress(i)   = AchievingPress(i) - 0.005d0*AchievingPress(i)
       HighRngePress(i)  = AchievingPress(i) + 0.005d0*AchievingPress(i)
       
    enddo
    
  end subroutine get_Achieving_HighRnge_LowRngePress
  
  
  subroutine get_Achieving_HighRnge_LowRngePress_generalized(CMv,PCSorICS,numOfCellsToAdjst,&
       AchievingPress,LowRngePress,HighRngePress)  
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    integer, intent(in)  :: numOfCellsToAdjst
    real*8,  intent(out) :: AchievingPress(1:numOfCellsToAdjst)
    real*8,  intent(out) :: LowRngePress(1:numOfCellsToAdjst)
    real*8,  intent(out) :: HighRngePress(1:numOfCellsToAdjst)
    
    integer :: i,imax,j,jmax
    integer :: lim1,lim2,lim3,lim4,lim5,lim6,lim7,lim8,lim9
    
    integer :: N_IncmingCellEquivPress,N_LowPressThnIncmCell
    integer :: N_LowrThanIncmngCellPress,N_MoreLowrThanIncmngCellPress,N_lowestCellPress
    integer :: N_highPressAftLowst_lvl1,N_highPressAftLowst_lvl2
    integer :: N_highPressAftLowst_lvl3,N_highPressAftLowst_lvl4
    
    imax = numOfCellsToAdjst
    
    if (CMv==2 .and. PCSorICS==1) then
       
       N_IncmingCellEquivPress = 6; lim1 = N_IncmingCellEquivPress   ; AchievingPress(1:lim1)        = 0.1000d0
       N_LowPressThnIncmCell   = 9; lim2 = N_LowPressThnIncmCell     ; AchievingPress((lim1+1):lim2) = 0.0930d0
       
    elseif (CMv==3 .and. PCSorICS==1) then
       
       N_IncmingCellEquivPress  = 5; lim1 = N_IncmingCellEquivPress  ; AchievingPress(1:lim1)        = 0.1000d0
       N_LowrThanIncmngCellPress= 6; lim2 = N_LowrThanIncmngCellPress; AchievingPress((lim1+1):lim2) = 0.0975d0
       N_lowestCellPress        = 7; lim3 = N_lowestCellPress        ; AchievingPress((lim2+1):lim3) = 0.0950d0
       N_highPressAftLowst_lvl1 = 8; lim4 = N_highPressAftLowst_lvl1 ; AchievingPress((lim3+1):lim4) = 0.1000d0
       N_highPressAftLowst_lvl2 = 9; lim5 = N_highPressAftLowst_lvl2 ; AchievingPress((lim4+1):lim5) = 0.1050d0
       N_highPressAftLowst_lvl3 = 10;lim6 = N_highPressAftLowst_lvl3 ; AchievingPress((lim5+1):lim6) = 0.1400d0
       N_highPressAftLowst_lvl4 = 11;lim7 = N_highPressAftLowst_lvl4 ; AchievingPress((lim6+1):lim7) = 0.1600d0
       
    elseif (CMv==1 .and. PCSorICS==1) then
       
       N_IncmingCellEquivPress  = 8; lim1 = N_IncmingCellEquivPress  ; AchievingPress(1:lim1)        = 0.1000d0
       N_highPressAftLowst_lvl1 = 9; lim2 = N_highPressAftLowst_lvl1 ; AchievingPress((lim1+1):lim2) = 0.1050d0
       N_highPressAftLowst_lvl2 = 10;lim3 = N_highPressAftLowst_lvl2 ; AchievingPress((lim2+1):lim3) = 0.1100d0
       N_highPressAftLowst_lvl3 = 11;lim4 = N_highPressAftLowst_lvl3 ; AchievingPress((lim3+1):lim4) = 0.1250d0
       
    elseif (CMv==4 .and. PCSorICS==1) then
       
       lim1 = 3  ; AchievingPress(1:lim1)        = 0.1000d0
       lim2 = 4  ; AchievingPress((lim1+1):lim2) = 0.1050d0
       lim3 = 5  ; AchievingPress((lim2+1):lim3) = 0.1025d0
       lim4 = 7  ; AchievingPress((lim3+1):lim4) = 0.1000d0
       lim5 = 8  ; AchievingPress((lim4+1):lim5) = 0.1075d0
       lim6 = 9  ; AchievingPress((lim5+1):lim6) = 0.1200d0
       lim7 = 10 ; AchievingPress((lim6+1):lim7) = 0.1350d0
       lim8 = 11 ; AchievingPress((lim7+1):lim8) = 0.1550d0
       
    endif
    
    do i = 1,imax
       LowRngePress(i)   = AchievingPress(i) - 0.005d0*AchievingPress(i)
       HighRngePress(i)  = AchievingPress(i) + 0.005d0*AchievingPress(i)
    enddo
    
  end subroutine get_Achieving_HighRnge_LowRngePress_generalized
  
    
  subroutine rscale_and_guessd_AchiLowHighPress_basedCurvtrOrientn(CMv,PCSorICS,TOL_prcnt&
       ,AchievingPress,LowRngePress,HighRngePress)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    real*8 , intent(in)  :: TOL_prcnt
    real*8 , intent(out) :: AchievingPress(1:(Hlf_Ncell+1))
    real*8 , intent(out) :: LowRngePress(1:(Hlf_Ncell+1))
    real*8 , intent(out) :: HighRngePress(1:(Hlf_Ncell+1))
    
    call rscale_basedOn_bndryPress() ; write(*,*) "Rescale at pos 00"
    
    AchievingPress(1) = 0.1000d0
    
    do i = 1,Hlf_Ncell
       
       if (CurvtrSign_LtPCS(CMv,i) ==  0) AchievingPress(i+1) = AchievingPress(i) 
       if (CurvtrSign_LtPCS(CMv,i) == +1) AchievingPress(i+1) = AchievingPress(i)-(0.100d0)*(BndryPress) 
       if (CurvtrSign_LtPCS(CMv,i) == -1) AchievingPress(i+1) = AchievingPress(i)+(0.200d0)*(BndryPress)
       
       if (i==1) then
          LowRngePress(i)   = AchievingPress(i) - (TOL_prcnt/100.0000d0)*AchievingPress(i)
          HighRngePress(i)  = AchievingPress(i) + (TOL_prcnt/100.0000d0)*AchievingPress(i)
          write(*,*) LowRngePress(i),AchievingPress(i),HighRngePress(i),i,"LAH 1"
       endif
       
       LowRngePress(i+1)   = AchievingPress(i+1) - (TOL_prcnt/100.0000d0)*AchievingPress(i+1)
       HighRngePress(i+1)  = AchievingPress(i+1) + (TOL_prcnt/100.0000d0)*AchievingPress(i+1)
       
       write(*,*) CurvtrSign_LtPCS(CMv,i),LowRngePress(i+1),AchievingPress(i+1),HighRngePress(i+1),i+1,"LAH 2"
       
    enddo
    
  end subroutine rscale_and_guessd_AchiLowHighPress_basedCurvtrOrientn
  
  
  subroutine rscale_basedOn_bndryPress()
    implicit none
    real*8 :: bndryPressCurr,rsFctr

    BndryPress     = 0.1000d0
    
    bndryPressCurr = k_area(1)*(A0(1)-A(1))    ; write(*,*)bndryPressCurr,BndryPress,"bndryPressCurr,BndryPress"
    rsFctr         = BndryPress/bndryPressCurr ; write(*,*)rsFctr,"rsFctr"
    call rescalingSystem_PresTnsnCg(rsFctr)
    write(*,*) "At rescaling based on BP, Frame_NI =", (Frame_NI-1)
    
  end subroutine rscale_basedOn_bndryPress
  
  subroutine get_Achieving_HighLowRnge_Curvtr(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,HighRngeCurvtr)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    real*8 , intent(out) :: AchievingCurvtr(1:Hlf_Ncell)
    real*8 , intent(out) :: LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    
    integer :: i,j,jmax
    real*8  :: ZERO=0.000000000d0,TINY=0.0000000001d0
    
    AchievingCurvtr(1:Hlf_Ncell) = abs(Expected_CHLtPCS(CMv,1:Hlf_Ncell))
    
    do i = 1,Hlf_Ncell
       
       if (CurvtrSign_LtPCS(CMv,i) == 0) then 
          LowRngeCurvtr(i)   = 0.0000d0 
          HighRngeCurvtr(i)  = 0.0200d0
          
       elseif (CurvtrSign_LtPCS(CMv,i) .ne. 0) then
          LowRngeCurvtr(i)   = AchievingCurvtr(i) - (0.02500d0)*AchievingCurvtr(i)
          HighRngeCurvtr(i)  = AchievingCurvtr(i) + (0.02500d0)*AchievingCurvtr(i)
       endif
       
       write(*,*) "LAH insd =",LowRngeCurvtr(i),AchievingCurvtr(i),HighRngeCurvtr(i),i
    enddo
    
    write(*,*) "WAITING FOR SOME SECONDS" ; call sleep(2)
    
  end subroutine get_Achieving_HighLowRnge_Curvtr
  
  subroutine get_Achieving_HighLowRnge_CurvtrRead(CMv,PCSorICS,AchievingCurvtr,LowRngeCurvtr,&
       HighRngeCurvtr,ExpctdCurvtr,CurvtrSignLt)
    ! same as the top routine, only diff it reads the curvature info when Expected_CHLtPCS and CurvtrSign_LtPCS is not allocated
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    real*8 , intent(out) :: AchievingCurvtr(1:Hlf_Ncell)
    real*8 , intent(out) :: LowRngeCurvtr(1:Hlf_Ncell),HighRngeCurvtr(1:Hlf_Ncell)
    real*8 , intent(out) :: ExpctdCurvtr(1:Hlf_Ncell)
    integer, intent(out) :: CurvtrSignLt(1:Hlf_Ncell)
    
    integer :: i,j,jmax
    real*8  :: ZERO=0.000000000d0,TINY=0.0000000001d0
    
    call read_expctd_curvtr_height_and_sign(CMv,ExpctdCurvtr,CurvtrSignLt)
    AchievingCurvtr(1:Hlf_Ncell) = abs(ExpctdCurvtr(1:Hlf_Ncell))
    
    do i = 1,Hlf_Ncell
       
       if (CurvtrSignLt(i) == 0) then 
          LowRngeCurvtr(i)   = 0.0000d0 
          HighRngeCurvtr(i)  = 0.0200d0
          
       elseif (CurvtrSignLt(i) .ne. 0) then
          LowRngeCurvtr(i)   = AchievingCurvtr(i) - (0.02500d0)*AchievingCurvtr(i)
          HighRngeCurvtr(i)  = AchievingCurvtr(i) + (0.02500d0)*AchievingCurvtr(i)
       endif
       
       write(*,*) "LAH insd =",LowRngeCurvtr(i),AchievingCurvtr(i),HighRngeCurvtr(i),i
    enddo
    
  end subroutine get_Achieving_HighLowRnge_CurvtrRead
  
  
  subroutine get_first_bndryPress_and_curvtrOrientatn_right(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    
    real*8  :: bndryPressCurr,rghtTolrncBP ! BP = Boundry Pressure
    integer :: i,j,jmax,cntLp,cnt_LpOrt     ! Ort=Orientation
    real*8  :: bracketPress
    integer :: cellBfr,cellAft
    real*8  :: pressBfr,pressAft,absPD
    integer :: signVal,ProceedVal
    real*8  :: rsFctr
    real*8  :: TINY=0.0500d0
    integer :: ManpltBfrAft,BsORLt
    real*8  :: curve_LtH(1:Hlf_Ncell)
    
    rghtTolrncBP   = BndryPress + 0.0050d0*BndryPress ; write(*,*) rghtTolrncBP,"rghtTolrncBP" 
    bracketPress   = 0.20d0*BndryPress
    
    cntLp = 1 ; write(*,*) Frame_NI-1,"Frame_NI before first_Bndry_Match"
    
    call rscale_basedOn_bndryPress() ; write(*,*) (Frame_NI-1),"Frame_NI aft rescling in orientn sbrtn"
    
    do
       
       bndryPressCurr = k_area(1) * (A0(1)-A(1))  ; write(*,*) bndryPressCurr,cntLp,"bndryPressCurr"
       
       if (abs(bndryPressCurr-BndryPress) .lt. TINY) then
          write(*,*) bndryPressCurr,BndryPress,"bfr exit"
          exit
          
       elseif (bndryPressCurr .lt. BndryPress) then 
          
          if ((BndryPress-bndryPressCurr).gt.bracketPress) A0(1) = 1.0050d0*A0(1)
          if ((BndryPress-bndryPressCurr).le.bracketPress) A0(1) = 1.0025d0*A0(1)
          
          A0(1+Hlf_Ncell) = A0(1)
          
       elseif(bndryPressCurr .gt. rghtTolrncBP) then
          
          if ((bndryPressCurr-rghtTolrncBP).gt.bracketPress) A0(1) = 0.9950d0*A0(1)
          if ((bndryPressCurr-rghtTolrncBP).le.bracketPress) A0(1) = 0.9975d0*A0(1)
          
          A0(1+Hlf_Ncell) = A0(1)
          
       elseif ((bndryPressCurr.gt.BndryPress) .and. (bndryPressCurr.le.rghtTolrncBP)) then
          write(*,*) "exiting for bndryPressCurr =",bndryPressCurr,"with cntLp =",cntLp
          exit
       endif
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       
    enddo
    
    write(*,*) Frame_NI-1,"Frame_NI after first_Bndry_Match"
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !ORIENTATION!
    
    ManpltBfrAft = 1 ! IMPORTANT --> DETERMINES WHICH SIDE OF THE CURVE'S CELL TO BE CHANGED
    
    do i = 1,Hlf_Ncell
       
       open(unit=542,file='cell_bfrAft.dat',position='append')
       write(542,*) " " ;  write(542,*) " " ; write(542,*) " " ;  write(542,*) " "
       
       cellBfr                     = i
       if (i.ne.Hlf_Ncell) cellAft = i+1 
       if (i == Hlf_Ncell) cellAft = N_cell
       
       write(542,*) cellBfr,cellAft,i,               "cell_BA_1"
       write(542,*) k_area(cellBfr),k_area(cellAft), "k_area_BA_1"
       write(542,*) A0(cellBfr),    A0(cellAft),     "A0_BA_1"
       write(542,*) A(cellBfr),     A(cellAft),      "A_BA_1"
       
       pressBfr  = k_area(cellBfr) * (A0(cellBfr)-A(cellBfr))
       pressAft  = k_area(cellAft) * (A0(cellAft)-A(cellAft))
       absPD     = abs(pressBfr-pressAft)
       
       write(542,*) pressBfr,pressAft,absPD,          "press_BA_diff"
       write(542,*) CMv,i," ",CurvtrSign_LtPCS(CMv,i),"CurvtrSign"
       
       if (CurvtrSign_LtPCS(CMv,i)==0) then ! curvtrDirctn ---|---
          
          cnt_LpOrt = 1
          do 
             if (absPD.gt.TolPressDiffToMaintnOrientn) then
                
                if (pressBfr .gt. pressAft) then
                   
                   write(542,*) " " ; write(542,*) " "
                   write(542,*) pressBfr,pressAft,absPD,TolPressDiffToMaintnOrientn,"PR... inside L00"
                   write(542,*) cellAft,A0(cellAft),"cellA,A0 lft inside L01"
                   if (cellAft == N_cell) continue
                   if (cellAft.ne.N_cell) write(542,*) (cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L00"
                   
                   if (ManpltBfrAft == 1) then
                      A0(cellBfr)                               = 0.9975d0*A0(cellBfr)
                      A0(cellBfr+Hlf_Ncell)                     = A0(cellBfr)
                      
                   elseif (ManpltBfrAft == 2) then
                      A0(cellAft)                               = 1.0025d0*A0(cellAft)
                      if (i.ne.Hlf_Ncell) A0(cellAft+Hlf_Ncell) = A0(cellAft)
                      if (i == Hlf_Ncell) continue
                   endif
                   
                elseif (pressBfr .lt. pressAft) then
                   
                   write(542,*) " " ; write(542,*) " "
                   write(542,*) pressBfr,pressAft,absPD,TolPressDiffToMaintnOrientn,"PR... inside L01"
                   write(542,*) cellAft,A0(cellAft),"cellA,A0 lft inside L01"
                   if (cellAft == N_cell) continue
                   if (cellAft.ne.N_cell) write(542,*) (cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L01"
                   
                   if (ManpltBfrAft == 1) then
                      A0(cellBfr)                               = 1.0025d0*A0(cellBfr)
                      A0(cellBfr+Hlf_Ncell)                     = A0(cellBfr)
                      
                   elseif (ManpltBfrAft == 2) then
                      A0(cellAft)                               = 0.9975d0*A0(cellAft)
                      if (i.ne.Hlf_Ncell) A0(cellAft+Hlf_Ncell) = A0(cellAft)
                      if (i == Hlf_Ncell) continue
                   endif
                   
                endif
                
                call Equilibrate_only_NI_model
                pressBfr  = k_area(cellBfr) * (A0(cellBfr)-A(cellBfr))
                pressAft  = k_area(cellAft) * (A0(cellAft)-A(cellAft))
                absPD     = abs(pressBfr-pressAft)
                cnt_LpOrt = cnt_LpOrt+1
                
             elseif (absPD.le.TolPressDiffToMaintnOrientn) then
                write(542,*) " "
                write(542,*) absPD,TolPressDiffToMaintnOrientn,cnt_LpOrt,"absPD,TolPressDiffBtwn"
                write(542,*)"less than TolPressDiff"
                exit
             endif
             
          enddo
          
       elseif (CurvtrSign_LtPCS(CMv,i)==+1) then ! curvtrDirctn ---)--- 
          
          cnt_LpOrt = 1
          do 
             if (pressBfr .gt. pressAft) then
                
                !%%%%%%%%%%%%%%% PRINTING VALUES START %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(542,*) " " ; write(542,*) " "
                write(542,*) pressBfr,pressAft,absPD,bracketPress,   "press and bracket inside L10"
                write(542,*) cellAft,A0(cellAft),                         "cellA,A0 lft inside L10"
                if (cellAft==N_cell)   continue
                if (cellAft.ne.N_cell) write(542,*) (cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L10"
                !%%%%%%%%%%%%%%% PRINTING VALUES FINISH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (absPD.le.TolPressDiffToMaintnOrientn) then
                   write(542,*) "+1 case, already (pressBfr>pressAft), but (absPD<TolPressToMaintainOrinetn)"
                   write(542,*) "To check this statement check the values:"
                   write(542,*) pressBfr,pressAft,absPD,TolPressDiffToMaintnOrientn,"pB,pA,absPD,TolP"
                   
                   if (ManpltBfrAft == 1) then
                      A0(cellBfr)                                  = 1.0025d0*A0(cellBfr)
                      A0(cellBfr+Hlf_Ncell)                        = A0(cellBfr)
                   elseif (ManpltBfrAft == 2) then
                      A0(cellAft)                                  = 0.9975d0*A0(cellAft)
                      if (cellAft.ne.N_cell) A0(cellAft+Hlf_Ncell) = A0(cellAft)
                      if (cellAft == N_cell) continue
                   endif
                   
                elseif (absPD.gt.TolPressDiffToMaintnOrientn) then
                   write(*,*) "(+1) case condn_passed for curveNo =",i,"Frame_NI =",(Frame_NI-1)
                   write(542,*) "(+1) case condn_passed for curveNo =",i,"Frame_NI =",(Frame_NI-1)
                   write(542,*) "and the values to pass the condition are:"
                   write(542,*) "pressBfr =",pressBfr,"pressAft =",pressAft
                   write(542,*) "absPD =",absPD,"TolPressDiffToMaintnOrientn",TolPressDiffToMaintnOrientn
                   write(542,*) " " ; write(542,*) " "
                   exit
                endif
                
             elseif (pressBfr .lt. pressAft) then
                
                !%%%%%%%%%%%%%%% PRINTING VALUES START %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(542,*) " " ; write(542,*) " "
                write(542,*) pressBfr,pressAft,absPD,bracketPress,   "press and bracket inside L11"
                write(542,*) cellAft,A0(cellAft),                         "cellA,A0 lft inside L11"
                if (cellAft == N_cell) continue
                if (cellAft.ne.N_cell) write(542,*) (cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L11"
                !%%%%%%%%%%%%%%% PRINTING VALUES FINISH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (ManpltBfrAft == 1) then
                   if (absPD .gt. bracketPress) A0(cellBfr)           = 1.0050d0*A0(cellBfr) 
                   if (absPD .lt. bracketPress) A0(cellBfr)           = 1.0025d0*A0(cellBfr)
                   A0(cellBfr+Hlf_Ncell)                              = A0(cellBfr)
                elseif (ManpltBfrAft == 2) then
                   if (absPD .gt. bracketPress) A0(cellAft)           = 0.9950d0*A0(cellAft) 
                   if (absPD .lt. bracketPress) A0(cellAft)           = 0.9975d0*A0(cellAft)
                   if (cellAft.ne.N_cell)       A0(cellAft+Hlf_Ncell) = A0(cellAft)
                   if (cellAft == N_cell)       continue
                endif
                
             endif
             
             call Equilibrate_only_NI_model
             pressBfr  = k_area(cellBfr) * (A0(cellBfr)-A(cellBfr))
             pressAft  = k_area(cellAft) * (A0(cellAft)-A(cellAft))
             absPD     = abs(pressBfr-pressAft)
             cnt_LpOrt = cnt_LpOrt+1
          enddo
          
       elseif (CurvtrSign_LtPCS(CMv,i)==-1) then ! curvtrDirctn ---(---
          
          cnt_LpOrt = 1
          do
             if (pressBfr .le. pressAft) then
                
                !%%%%%%%%%%%%%%% PRINTING VALUES START %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(542,*) " " ; write(542,*) " "
                write(542,*) pressBfr,pressAft,absPD,bracketPress,   "press and bracket inside L10"
                write(542,*) cellAft,A0(cellAft),                         "cellA,A0 lft inside L10"
                if (cellAft==N_cell)   continue
                if (cellAft.ne.N_cell) write(542,*)(cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L10"
                
                !%%%%%%%%%%%%%%% PRINTING VALUES FINISH %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (absPD.le.TolPressDiffToMaintnOrientn) then
                   write(542,*) "-1 case, already (pressBfr<pressAft), but (absPD<TolPressToMaintainOrinetn)"
                   write(542,*) "To check this statement check the values:"
                   write(542,*) pressBfr,pressAft,absPD,TolPressDiffToMaintnOrientn,"pB,pA,absPD,TolP"
                   
                   if (ManpltBfrAft == 1) then
                      A0(cellBfr)                                    = 0.9975d0*A0(cellBfr)
                      A0(cellBfr+Hlf_Ncell)                          = A0(cellBfr)
                   elseif (ManpltBfrAft == 2) then
                      A0(cellAft)                                    = 1.0025d0*A0(cellAft)
                      if (cellAft.ne.N_cell)   A0(cellAft+Hlf_Ncell) = A0(cellAft)
                      if (cellAft == N_cell)   continue
                   endif
                   
                elseif (absPD.gt.TolPressDiffToMaintnOrientn) then
                   write(*,*) "(-1) case condn_passed for curveNo =",i,"Frame_NI =",(Frame_NI-1)
                   write(542,*) "(-1) case condn_passed for curveNo =",i,"Frame_NI =",(Frame_NI-1)
                   write(542,*) "and the values to pass the condition are:"
                   write(542,*) "pressBfr =",pressBfr,"pressAft =",pressAft
                   write(542,*) "absPD =",absPD,"TolPressDiffToMaintnOrientn",TolPressDiffToMaintnOrientn
                   write(542,*) " " ; write(542,*) " "
                   exit
                endif
                
             elseif (pressBfr .gt. pressAft) then
                
                write(542,*) " " ; write(542,*) " "
                write(542,*) pressBfr,pressAft,absPD,bracketPress,"press and bracket inside L12"
                write(542,*) cellAft,A0(cellAft),"cellA,A0 lft inside L12"
                if (cellAft==N_cell)   continue
                if (cellAft.ne.N_cell) write(542,*)(cellAft+Hlf_Ncell),Hlf_Ncell,A0(cellAft+Hlf_Ncell),"cellA,A0 rgt L12"
                
                if (ManpltBfrAft == 1) then
                   if (absPD .gt. bracketPress) A0(cellBfr)           = 0.9950d0*A0(cellBfr)
                   if (absPD .lt. bracketPress) A0(cellBfr)           = 0.9975d0*A0(cellBfr)
                   A0(cellBfr+Hlf_Ncell)                              = A0(cellBfr)
                elseif (ManpltBfrAft == 2) then
                   if (absPD .gt. bracketPress) A0(cellAft)           = 1.0050d0*A0(cellAft)
                   if (absPD .lt. bracketPress) A0(cellAft)           = 1.0025d0*A0(cellAft)
                   if (cellAft.ne.N_cell)       A0(cellAft+Hlf_Ncell) = A0(cellAft)
                   if (cellAft == N_cell)       continue   
                endif
                
             endif
             
             call Equilibrate_only_NI_model
             pressBfr  = k_area(cellBfr) * (A0(cellBfr)-A(cellBfr))
             pressAft  = k_area(cellAft) * (A0(cellAft)-A(cellAft))
             absPD     = abs(pressBfr-pressAft)
             cnt_LpOrt = cnt_LpOrt+1
             
          enddo
             
       endif
       
       close(542)
       !stop 'VE'
    enddo
    
    !close(542)
    
    write(*,*) Frame_NI-1,"Frame_NI after Orientation_Match"
    
    
    do i = 1,Hlf_Ncell
       
       pressBfr = k_area(i)*(A0(i)-A(i))
       
       if (i.ne.Hlf_Ncell) pressAft = k_area(i+1)    * (A0(i+1)   - A(i+1))
       if (i == Hlf_Ncell) pressAft = k_area(N_cell) * (A0(N_cell)- A(N_cell))
       
       if (pressBfr.gt.pressAft) signVal = +1
       if (pressBfr.lt.pressAft) signVal = -1
       
       write(*,*) signVal,CurvtrSign_LtPCS(CMv,i),i,"currentSign,actualSign,i"
    enddo
    
    write(*,*) Frame_NI-1,'aft orientation 1'
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
    do i = 1,Hlf_Ncell
       write(*,*) curve_LtH(i),i,"curve_LtH-i"
    enddo
    
    write(*,*) Frame_NI-1,'aft orientation 2'
    
    !stop 'to check aft curvtrOrientn_right'
    
  end subroutine get_first_bndryPress_and_curvtrOrientatn_right
  
  
  subroutine matchBiggstCurve_and_Achieving_HighLowRngePress(CMv,PCSorICS,AchievingPress,&
       TOL_prcnt,LowRngePress,HighRngePress)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    real*8,  intent(out) :: AchievingPress(1:(Hlf_Ncell+1))
    real*8,  intent(in)  :: TOL_prcnt
    real*8,  intent(out) :: LowRngePress(1:(Hlf_Ncell+1))
    real*8,  intent(out) :: HighRngePress(1:(Hlf_Ncell+1))
    
    integer :: i,j,imax,jmax
    real*8  :: Expected_CurvtrPCS(1:Hlf_Ncell)
    real*8  :: max_curvtrPCS
    integer :: maxCurvtrNum
    real*8  :: TolL_maxCurvtrPCS,TolG_maxCurvtrPCS,currMaxCurvtrPCS
    real*8  :: PressDiffBtwnMaxCurvtrCells,cf_curvtr,PD_cells(1:Hlf_Ncell)
    
    
    Expected_CurvtrPCS(1:Hlf_Ncell) = Expected_CHLtPCS(CMv,1:Hlf_Ncell)
    call get_max_Curvtr(Expected_CurvtrPCS,max_curvtrPCS,maxCurvtrNum,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS)
    call match_the_Expected_curvtr_for_maxOne(maxCurvtrNum,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS,Expected_CurvtrPCS)
    
    call find_the_PD_for_max_curvtr(maxCurvtrNum,PressDiffBtwnMaxCurvtrCells)
    call find_the_coeff_ofPD_and_other_PDs(PressDiffBtwnMaxCurvtrCells,max_curvtrPCS,Expected_CurvtrPCS,&
         cf_curvtr,PD_cells)
    
    call find_the_achieving_low_high_press(PD_cells,TOL_prcnt,AchievingPress,LowRngePress,HighRngePress)
    call write_Vals_after_matching_BiggstCUrvtr(CMv,Expected_CurvtrPCS,max_curvtrPCS,maxCurvtrNum,&
       PressDiffBtwnMaxCurvtrCells,cf_curvtr,PD_cells,AchievingPress,LowRngePress,HighRngePress)
    
  end subroutine matchBiggstCurve_and_Achieving_HighLowRngePress
  
  
  subroutine get_max_Curvtr(curvtrVal,maxCurvtrVal,maxCurvtrNum,TolL_maxCurvtr,TolR_maxCurvtr)
    implicit none
    real*8,  intent(in)  :: curvtrVal(1:Hlf_Ncell) 
    real*8,  intent(out) :: maxCurvtrVal
    integer, intent(out) :: maxCurvtrNum
    real*8,  intent(out) :: TolL_maxCurvtr,TolR_maxCurvtr
    
    integer :: i,j,jmax
    real*8  :: TolFctr
    
    ! %%%%%%%%%%%%___EXAMPLE_OF_ONE_CALL___%%%%%%%%%%%%%
    ! it was called inside the ROUTINE named: matchBiggstCurve_and_Achieving_HighLowRNgePress
    ! get_max_Curvtr(Expected_CurvtrPCS,max_curvtrPCS,maxCurvtrNum,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS)
    ! curvtrVal(1:Hlf_Ncell) ----- >>>  Expected_CurvtrPCS(CMv,1:Hlf_Ncell)
    ! maxCurvtrVal ---- >>> max_CurvtrPCS
    ! maxCurvtrNum ---- >>> maxCurvtrNum
    ! TolL_maxCurvtr --- >>> TolL_maxCurvtrPCS
    ! TolR_maxCurvtr --- >>> TolR_maxCurvtrPCS
    ! %%%%%%%%%%%%___EXAMPLE_DONE___%%%%%%%%%%%%%
    
    maxCurvtrVal = 0.00d0 ; maxCurvtrNum = -1
    
    do i = 1,Hlf_Ncell   
       if (abs(curvtrVal(i)) .gt. maxCurvtrVal) then
          maxCurvtrVal = abs(curvtrVal(i))
          maxCurvtrNum = i
       endif
    enddo
    
    TolFctr        = 0.000100d0
    TolL_maxCurvtr = (1.0000d0-TolFctr)*(maxCurvtrVal)
    TolR_maxCurvtr = (1.0000d0+TolFctr)*(maxCurvtrVal)
    
    write(*,*) maxCurvtrVal,maxCurvtrNum,"maxCurvtrPCS,maxCurvtrNum"
    write(*,*) TolL_maxCurvtr,TolR_maxCurvtr,"Tol(L/R)_maxCurvtr"
    
  end subroutine get_max_Curvtr
  
  
  subroutine match_the_Expected_curvtr_for_maxOne(maxCurvtrNum,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS,Expected_CurvtrPCS)
    implicit none
    integer, intent(in) :: maxCurvtrNum
    real*8 , intent(in) :: TolL_maxCurvtrPCS,TolG_maxCurvtrPCS
    real*8 , intent(in) :: Expected_CurvtrPCS(1:Hlf_Ncell)
    
    integer :: BsORLt
    real*8  :: Curve_LtH(1:Hlf_Ncell),currMaxCurvtrPCS
    integer :: i,j,jmax,cntLp
    
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,Curve_LtH) ; currMaxCurvtrPCS = Curve_LtH(maxCurvtrNum)
    
    cntLp = 1
    
    do
       
       if (currMaxCurvtrPCS .gt. TolG_maxCurvtrPCS) then
          
          if (Expected_CurvtrPCS(maxCurvtrNum) .gt. ZERO) then ! cell Left High Pressure 
             A0(maxCurvtrNum) = (0.9975d0)*A0(maxCurvtrNum)
             A0(maxCurvtrNum+Hlf_Ncell) = A0(maxCurvtrNum)
             
          elseif (Expected_CurvtrPCS(maxCurvtrNum) .lt. ZERO) then ! cell Left Low Pressure
             A0(maxCurvtrNum)           = (1.0025d0)*A0(maxCurvtrNum)
             A0(maxCurvtrNum+Hlf_Ncell) = A0(maxCurvtrNum)
          endif
          
       elseif (currMaxCurvtrPCS .le. TolL_maxCurvtrPCS) then
          
          if (Expected_CurvtrPCS(maxCurvtrNum) .gt. ZERO) then ! cell Left High Pressure 
             A0(maxCurvtrNum)           = (1.0025d0)*A0(maxCurvtrNum)
             A0(maxCurvtrNum+Hlf_Ncell) = A0(maxCurvtrNum)
             
          elseif (Expected_CurvtrPCS(maxCurvtrNum) .lt. ZERO) then ! cell Left Low Pressure
             A0(maxCurvtrNum)           = (0.9975d0)*A0(maxCurvtrNum)
             A0(maxCurvtrNum+Hlf_Ncell) = A0(maxCurvtrNum)
          endif
          
       elseif ( (currMaxCurvtrPCS .gt.TolL_maxCurvtrPCS) .and. (currMaxCurvtrPCS.le.TolG_maxCurvtrPCS) ) then
          write(*,*) currMaxCurvtrPCS,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS,"currMax"
          write(*,*) cntLp,"cntLp during Exit"
          exit
       endif
       
       call Equilibrate_only_NI_model
       cntLp = cntLp+1
       call get_Curvtr_Height_general(BsORLt,Curve_LtH) ; currMaxCurvtrPCS = Curve_LtH(maxCurvtrNum)
       write(*,*) currMaxCurvtrPCS,TolL_maxCurvtrPCS,TolG_maxCurvtrPCS,"currMax-NonExit"
       
    enddo
    
    
  end subroutine match_the_Expected_curvtr_for_maxOne
  
  
  subroutine find_the_PD_for_max_curvtr(maxCurvtrNum,PressDiffBtwnMaxCurvtrCells)
    implicit none
    integer, intent(in)  :: maxCurvtrNum
    real*8,  intent(out) :: PressDiffBtwnMaxCurvtrCells
    
    integer :: lftCellmaxCurvtr,rgtCellmaxCurvtr
    real*8  :: PressLftCellMaxCurvtr,PressRgtCellMaxCurvtr
    
    lftCellmaxCurvtr = maxCurvtrNum
    if (maxCurvtrNum .ne. Hlf_Ncell) rgtCellmaxCurvtr = lftCellmaxCurvtr+1 
    if (maxCurvtrNum  ==  Hlf_Ncell) rgtCellmaxCurvtr = N_cell
    
    PressLftCellMaxCurvtr = k_area(lftCellmaxCurvtr) * (A0(lftCellmaxCurvtr)-A(lftCellmaxCurvtr))
    PressRgtCellMaxCurvtr = k_area(rgtCellmaxCurvtr) * (A0(rgtCellmaxCurvtr)-A(rgtCellmaxCurvtr))
    
    PressDiffBtwnMaxCurvtrCells = abs(PressLftCellMaxCurvtr - PressRgtCellMaxCurvtr) 
    write(*,*) PressLftCellMaxCurvtr,PressRgtCellMaxCurvtr,PressDiffBtwnMaxCurvtrCells,"PresDiff"
  
  end subroutine find_the_PD_for_max_curvtr
  
  
  subroutine find_the_coeff_ofPD_and_other_PDs(PressDiffBtwnMaxCurvtrCells,max_curvtrPCS,Expected_CurvtrPCS,&
       cf_curvtr,PD_cells)
    implicit none
    real*8, intent(in)  :: PressDiffBtwnMaxCurvtrCells, max_curvtrPCS, Expected_CurvtrPCS(1:Hlf_Ncell)
    real*8, intent(out) :: cf_curvtr,PD_cells(1:Hlf_Ncell)
    integer             :: i,j
    
    cf_curvtr = PressDiffBtwnMaxCurvtrCells/max_curvtrPCS
    write(*,*) cf_curvtr,PressDiffBtwnMaxCurvtrCells,max_curvtrPCS,"cf_curvtr,PD_maxCells,max_cur"
    
    do i = 1,Hlf_Ncell
       PD_cells(i) = cf_curvtr * Expected_CurvtrPCS(i)
       write(*,*) i,Expected_CurvtrPCS(i),PD_cells(i),"Expected_CurvtrPCS-PD_cells"
    enddo
    
  end subroutine find_the_coeff_ofPD_and_other_PDs
  
  subroutine find_the_achieving_low_high_press(PD_cells,TOL_prcnt,AchievingPress,LowRngePress,HighRngePress)
    implicit none
    real*8, intent(in)  :: PD_cells(1:Hlf_Ncell)
    real*8, intent(in)  :: TOL_prcnt
    real*8, intent(out) :: AchievingPress(1:(Hlf_Ncell+1)),LowRngePress(1:(Hlf_Ncell+1))
    real*8, intent(out) :: HighRngePress(1:(Hlf_Ncell+1)) 
    
    integer :: i,imax
    
    do i = 1,Hlf_Ncell
       
       if (i==1) AchievingPress(i) = 0.10d0
       AchievingPress(i+1) = AchievingPress(i) - PD_cells(i)
       
       write(*,*) AchievingPress(i),i,"Ach-i"
       if (i==Hlf_Ncell) write(*,*) AchievingPress(i+1),N_cell,"Ach-i"
    enddo
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       LowRngePress(i)  = AchievingPress(i) - (TOL_prcnt/100.0000000d0)*(AchievingPress(i)) !-0.005
       HighRngePress(i) = AchievingPress(i) + (TOL_prcnt/100.0000000d0)*(AchievingPress(i)) !+0.005
    enddo
    
    
  end subroutine find_the_achieving_low_high_press
  
  
  subroutine write_Vals_after_matching_BiggstCUrvtr(CMv,Expected_CurvtrPCS,max_curvtrPCS,maxCurvtrNum,&
       PressDiffBtwnMaxCurvtrCells,cf_curvtr,PD_cells,AchievingPress,LowRngePress,HighRngePress)
    implicit none
    integer , intent(in) :: CMv
    real*8  , intent(in) :: Expected_CurvtrPCS(1:Hlf_Ncell)
    real*8  , intent(in) :: max_curvtrPCS
    integer , intent(in) :: maxCurvtrNum
    real*8  , intent(in) :: PressDiffBtwnMaxCurvtrCells
    real*8  , intent(in) :: cf_curvtr
    real*8  , intent(in) :: PD_cells(1:Hlf_Ncell)
    real*8  , intent(in) :: AchievingPress(1:(Hlf_Ncell+1))
    real*8  , intent(in) :: LowRngePress(1:(Hlf_Ncell+1))
    real*8  , intent(in) :: HighRngePress(1:(Hlf_Ncell+1))
    
    integer              :: i,j,jmax
    character(len=100)   :: fileNamePrfx1,fileNamePrfx2,fileNameSufx,fileName1,fileName2
    
    fileNamePrfx1='InfoForBiggstCurvtrMatch_WTExplntn_CMv'
    fileNamePrfx2='InfoForBiggstCurvtrMatch_WOExplntn_CMv'
    
    write(fileNameSufx,'(i2.2,a)')CMv,'PCS.dat'
    
    fileName1=trim(adjustl(fileNamePrfx1))//trim(adjustl(fileNameSufx))
    fileName2=trim(adjustl(fileNamePrfx2))//trim(adjustl(fileNameSufx))
    
    write(*,*)trim(adjustl(fileName1)) ; write(*,*)trim(adjustl(fileName2))
    
    open(unit=562,file=trim(adjustl(fileName1)))
    open(unit=563,file=trim(adjustl(fileName2)))
    
    write(562,*) "This file is to print the data regarding the biggest curvature matching Data"
    write(562,*) "The CMv value is =",CMv
    write(563,*) CMv
    write(562,*) " "

    write(562,*) "The Expected Curvatures are:"
    
    do i = 1,Hlf_Ncell
       write(562,*) "Expected_Curvtr is =", Expected_CurvtrPCS(i),"for i =",i
       write(563,*) Expected_CurvtrPCS(i),i
    enddo

    write(562,*) " "
    write(562,*) "The maxCurvtr amount is =", max_CurvtrPCS, "for the curvature num =",maxCurvtrNum
    write(563,*) max_CurvtrPCS,maxCurvtrNum
    write(562,*) " "
    
    write(562,*) "The pressure difference between cells of max curvtr is =",PressDiffBtwnMaxCurvtrCells
    write(562,*) "The coeffiecient of curvature is =", cf_curvtr
    write(563,*) PressDiffBtwnMaxCurvtrCells
    write(563,*) cf_curvtr
    
    write(562,*) " "
    write(562,*) "The Pressure Difference between cells is ="
    
    do i = 1,Hlf_Ncell
       write(562,*) "Pressure diff is =",PD_cells(i),"for cuvature num =",i
       write(563,*) PD_cells(i),i
    enddo
    
    write(562,*) " "
    write(562,*) "The Achieving,High-Low Range press is = "
    
    do i = 1,(Hlf_Ncell+1)
       write(562,*) "Achieving-High-Low Range is =",AchievingPress(i),HighRngePress(i),LowRngePress(i),"for curvature =",i
       write(563,*) AchievingPress(i),HighRngePress(i),LowRngePress(i),i
    enddo
    
    write(562,*) " "
    
    close(562)
    close(563)
    
  end subroutine write_Vals_after_matching_BiggstCUrvtr
  
  subroutine get_Simulted_Curvtr_Height(CMv,PCSorICS,BsORLt)
    implicit none
    integer, intent(in) :: CMv,PCSorICS,BsORLt
    
    real*8  :: xVB(1:Hlf_Ncell,1:(NAEC_Ltrl+2)),yVB(1:Hlf_Ncell,1:(NAEC_Ltrl+2))
    real*8  :: xVL(1:Hlf_Ncell,1:(NAEC_Ltrl+2)),yVL(1:Hlf_Ncell,1:(NAEC_Ltrl+2))
    real*8  :: xML(1:Hlf_Ncell),yML(1:Hlf_Ncell),yML_WOCF(1:Hlf_Ncell)            !WOCF = WithOUT CoeFficient
    real*8  :: Curve_H(1:Hlf_Ncell)
    
    integer :: NAEC_val
    integer :: i,j,jmax
    
    real*8, allocatable :: CF_mat(:,:)
    
    
    if (BsORLt==1) NAEC_val=NAEC_Bsal
    if (BsORLt==2) NAEC_val=NAEC_Ltrl
    
    allocate(CF_mat(1:Hlf_Ncell,1:(NAEC_val+2))) ; CF_mat = -1.0d30
    
    call get_xVyV_for_curvtr(BsORLt,NAEC_val,xVL,yVL)   ! READ THE NODES COMMON FOR ALL APPROACH
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! PREDICT CURVATURE WITH FOUR NODES [DOES NOT WORK GOOD, WILL SEE WHAT HAPPENS]
    
    call get_xM_and_yMWOCF_for_curvtr(NAEC_val,xVL,yVL,xML,yML_WOCF)
    call use_xVyV_to_find_Coeff_mat_FortheCurve(BsORLt,NAEC_val,xVL,yVL,CF_mat)
    call get_curve_height_using_coeffs(CF_mat,NAEC_val,xML,yML)
    
    do i = 1,Hlf_Ncell
       write(*,*) i,xML(i),yML(i),yML_WOCF(i),(yML(i)-yML_WOCF(i)),"xmL-yML and Compr"
    enddo
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CURVE HEIGHT WITHOUT CURVATURE PREDICTION
    
    call get_curve_height_based_on_PLorMP_method(NAEC_val,xVL,yVL,Curve_H)
    
    if (BsOrLt==1 .and. PCSorICS==1) Simulted_CHBsPCS(CMv,1:Hlf_Ncell) = Curve_H(1:Hlf_Ncell)
    if (BsORLt==2 .and. PCSorICS==1) Simulted_CHLtPCS(CMv,1:Hlf_Ncell) = Curve_H(1:Hlf_Ncell)
    if (BsOrLt==1 .and. PCSorICS==2) write(*,*)"Curve_H not needed for ICS since no Exprmnt to cmpr1"
    if (BsOrLt==2 .and. PCSorICS==2) write(*,*)"Curve_H not needed for ICS since no Exprmnt to cmpr2"
    
  end subroutine get_Simulted_Curvtr_Height
  
  
  subroutine get_Curvtr_Height_general(BsORLt,Curve_H)
    implicit none
    integer, intent(in)  :: BsOrLt
    real*8,  intent(out) :: Curve_H(1:Hlf_Ncell)
    
    real*8  :: xVB(1:Hlf_Ncell,1:(NAEC_Ltrl+2)),yVB(1:Hlf_Ncell,1:(NAEC_Ltrl+2))
    real*8  :: xVL(1:Hlf_Ncell,1:(NAEC_Ltrl+2)),yVL(1:Hlf_Ncell,1:(NAEC_Ltrl+2))
    integer :: NAEC_val
    integer :: i,j,jmax
    
    if (BsORLt==1) NAEC_val=NAEC_Bsal
    if (BsORLt==2) NAEC_val=NAEC_Ltrl
    
    call get_xVyV_for_curvtr(BsORLt,NAEC_val,xVL,yVL)   ! READ THE NODES
    call get_curve_height_based_on_PLorMP_method(NAEC_val,xVL,yVL,Curve_H)
    
  end subroutine get_Curvtr_Height_general
  
  subroutine get_xVyV_for_curvtr(BsORLt,NAEC_val,xV,yV)
    implicit none
    integer, intent(in)  :: BsORLt,NAEC_val
    real*8 , intent(out) :: xV(1:Hlf_Ncell,1:(NAEC_val+2)),yV(1:Hlf_Ncell,1:(NAEC_val+2))
    
    integer  :: i,j,jmax
    integer  :: nsprsInACell
    integer  :: BsSpr1,BsSpr2,BsSpr3,LtSpr1,LtSpr2,LtSpr3
    integer  :: nodeNums(1:(NAEC_val+2))
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i = 1,Hlf_Ncell
       
       BsSpr1=(i-1)*(nsprsInACell)+2; BsSpr2=(i-1)*(nsprsInACell)+3; BsSpr3=(i-1)*(nsprsInACell)+4
       LtSpr1=(i-1)*(nsprsInACell)+5; LtSpr2=(i-1)*(nsprsInACell)+6; LtSpr3=(i-1)*(nsprsInACell)+7

       if (writeForCurvtr==1) then
          if (i==1 .or. i==2) then
             write(*,*) BsSpr1,BsSpr2,BsSpr3,"BsSpr 1-2-3"
             write(*,*) LtSpr1,LtSpr2,LtSpr3,"LtSpr 1-2-3"
          endif
       endif
       
       if (BsORLt==1) then
          nodeNums(1) = spr_node(BsSpr1,1) ; nodeNums(2) = spr_node(BsSpr1,2)
          nodeNums(3) = spr_node(BsSpr2,2) ; nodeNums(4) = spr_node(BsSpr3,2)
       elseif (BsORLt==2) then
          nodeNums(1) = spr_node(LtSpr1,1) ; nodeNums(2) = spr_node(LtSpr1,2)
          nodeNums(3) = spr_node(LtSpr2,2) ; nodeNums(4) = spr_node(LtSpr3,2)
       endif
       
       jmax=NAEC_val+2
       
       do j = 1,jmax
          xV(i,j) = node_xy(nodeNums(j),1)
          yV(i,j) = node_xy(nodeNums(j),2)
       enddo
       
       if (writeForCurvtr==1) then
          if (i==1 .or. i==2) then
             write(*,*) nodeNums(1),xV(i,1),yV(i,1),"xV-yV 1"
             write(*,*) nodeNums(2),xV(i,2),yV(i,2),"xV-yV 2"
             write(*,*) nodeNums(3),xV(i,3),yV(i,3),"xV-yV 3"
             write(*,*) nodeNums(4),xV(i,4),yV(i,4),"xV-yV 4"
          endif
       endif
       
    enddo
    
  end subroutine get_xVyV_for_curvtr
  
  subroutine get_xM_and_yMWOCF_for_curvtr(NAEC_val,xV,yV,xM,yM_WOCF)
    implicit none
    integer, intent(in)  :: NAEC_val
    real*8,  intent(in)  :: xV(1:Hlf_Ncell,1:(NAEC_val+2)),yV(1:Hlf_Ncell,1:(NAEC_val+2))
    real*8,  intent(out) :: xM(1:Hlf_Ncell),yM_WOCF(1:Hlf_Ncell) ! WOCF = WithOUT CoeFficient
    
    integer :: i,imax
    
    do i = 1,Hlf_Ncell
       xM(i)      = 0.5000d0*(xV(i,2)+xV(i,3))
       yM_WOCF(i) = 0.5000d0*(yV(i,2)+yV(i,3))
    enddo
    
  end subroutine get_xM_and_yMWOCF_for_curvtr
  
  
  subroutine use_xVyV_to_find_Coeff_mat_FortheCurve(BsORLt,NAEC_val,xV,yV,CF_mat)
    implicit none
    integer, intent(in)  :: BsORLt,NAEC_val
    real*8 , intent(in)  :: xV(1:Hlf_Ncell,1:(NAEC_val+2)),yV(1:Hlf_Ncell,1:(NAEC_val+2))
    real*8 , intent(out) :: CF_mat(1:Hlf_Ncell,1:(NAEC_val+2))
    
    integer :: i,j,jmax
    real*8  :: x1,y1,x2,y2,x3,y3,x4,y4
    real*8  :: xValue(1:(NAEC_val+2)),yValue(1:(NAEC_val+2))
    real*8  :: A_mat(1:(NAEC_val+2),1:(NAEC_val+2)),B_mat(1:(NAEC_val+2))
    real*8  :: A_matInv(1:(NAEC_val+2),1:(NAEC_val+2)),cf_matValue(1:(NAEC_val+2))
    
    interface
       function MAT_INV4b4(A1)
         implicit none
         real*8, intent(in)  :: A1(4,4)
         real*8              :: MAT_INV4b4(4,4) 
       end function MAT_INV4b4
    end interface
    
    
    do i = 1,Hlf_Ncell
       
       do j = 1,(NAEC_val+2)   
          xValue(j) = xV(i,j) ; yValue(j) = yV(i,j)
       enddo
       
       call get_A_and_B_matrix(xValue,yValue,NAEC_val,A_mat,B_mat)
       
       A_matInv                 = MAT_INV4b4(A_mat)
       cf_matValue              = matmul(A_matInv,B_mat)
       CF_mat(i,1:(NAEC_val+2)) = cf_matValue(1:(NAEC_val+2)) 
       
       if (i==7) then
          
          write(*,*) " "
          write(*,*) node_xy(15,1),node_xy(76,1),node_xy(77,1),node_xy(16,1),"nodeX"
          write(*,*) xValue(1:4),"xValue"
          write(*,*) node_xy(15,2),node_xy(76,2),node_xy(77,2),node_xy(16,2),"nodeX"
          write(*,*) yValue(1:4),"yValue"
          write(*,*) " "
          write(*,*) (xValue(1))**3,(xValue(1))**2,xValue(1),"x1"
          write(*,*) A_mat(1,1:4),"A row1"
          write(*,*) (xValue(2))**3,(xValue(2))**2,xValue(2),"X2"
          write(*,*) A_mat(2,1:4),"A row2"
          write(*,*) (xValue(3))**3,(xValue(3))**2,xValue(3),"X3"
          write(*,*) A_mat(3,1:4),"A row3"
          write(*,*) (xValue(4))**3,(xValue(4))**2,xValue(4),"X4"
          write(*,*) A_mat(4,1:4),"A row4"
          write(*,*) yValue,"yval"
          write(*,*) B_mat,"B_mat"
          write(*,*) " "
          write(*,*) A_matInV(1,1:4),"A inv Row1"
          write(*,*) A_matInV(2,1:4),"A inv Row2"
          write(*,*) A_matInV(3,1:4),"A inv Row3"
          write(*,*) A_matInV(4,1:4),"A inv Row4"
          write(*,*) " "
          write(*,*) cf_matValue(1:(NAEC_val+2)),"CF_matV"
          
       endif
       
    enddo
    
    
  end subroutine use_xVyV_to_find_Coeff_mat_FortheCurve
  
  subroutine get_curve_height_using_coeffs(CF_mat,NAEC_val,xM,yM)
    implicit none
    real*8,  intent(in)  :: CF_mat(1:Hlf_Ncell,1:(NAEC_val+2)),xM(1:Hlf_Ncell)
    integer, intent(in)  :: NAEC_val
    real*8,  intent(out) :: yM(1:Hlf_Ncell)
    
    integer :: i,j,jmax
    
    do i = 1,Hlf_Ncell
       yM(i) = CF_mat(i,1)*(xM(i))**3 + CF_mat(i,2)*(xM(i))**2 + CF_mat(i,3)*xM(i) + CF_mat(i,4)
       write(*,*) xM(i),yM(i),i,"xM-yM-i"
    enddo
    
  end subroutine get_curve_height_using_coeffs
  
  subroutine get_A_and_B_matrix(xValue,yValue,NAEC_val,A_mat,B_mat)
    implicit none
    real*8,  intent(in)  :: xValue(1:(NAEC_val+2)),yValue(1:(NAEC_val+2))
    integer, intent(in)  :: NAEC_val
    real*8,  intent(out) :: A_mat(1:(NAEC_val+2),1:(NAEC_val+2)),B_mat(1:(NAEC_val+2))
    integer              :: i,imax,j,jmax
    
    imax=NAEC_val+2 ; jmax=NAEC_val+2
    
    do i = 1,imax
       
       do j = 1,jmax   
          if (j==1) A_mat(i,j) = xValue(i)**3
          if (j==2) A_mat(i,j) = xValue(i)**2
          if (j==3) A_mat(i,j) = xValue(i)
          if (j==4) A_mat(i,j) = 1
       enddo
       
       B_mat(i) = yValue(i)
       
       write(*,*) xValue(i),yValue(i),i,"x-y-i"
       write(*,*) A_mat(i,1:(NAEC_val+2)),"A-mat"
       write(*,*) B_mat(i),"B-mat"
       
    enddo
    
  end subroutine get_A_and_B_matrix
  
  subroutine get_xM_by_differentiation
    implicit none
    continue
  end subroutine get_xM_by_differentiation
  
  subroutine get_curve_height_based_on_PLorMP_method(NAEC_val,xV,yV,Curve_H)
    implicit none
    integer, intent(in)  :: NAEC_val
    real*8,  intent(in)  :: xV(1:Hlf_Ncell,1:(NAEC_val+2)),yV(1:Hlf_Ncell,1:(NAEC_val+2))
    real*8,  intent(out) :: Curve_H(1:Hlf_Ncell)
    integer              :: i,j,imax,jmax
    real*8               :: curve_distnce
    real*8               :: Line1N(1:2,1:2),Line2N(1:2,1:2)
    
    do i = 1,Hlf_Ncell
       
       Line1N(1,1)=xV(i,1) ; Line1N(1,2)=yV(i,1) ; Line1N(2,1)=xV(i,4) ; Line1N(2,2)=yV(i,4)
       Line2N(1,1)=xV(i,2) ; Line2N(1,2)=yV(i,2) ; Line2N(2,1)=xV(i,3) ; Line2N(2,2)=yV(i,3)
       
       call curve_dist_measuremnt_using_PLorMP_methd(i,Line1N,Line2N,curve_distnce)
       Curve_H(i) = curve_distnce
       write(*,*) Curve_H(i),i,"Height"
       
    enddo
    
  end subroutine get_curve_height_based_on_PLorMP_method
  
  
  subroutine curve_dist_measuremnt_using_PLorMP_methd(serialNo,Line1N,Line2N,curve_distnce)
    implicit none
    integer, intent(in)  :: serialNo
    real*8,  intent(in)  :: Line1N(1:2,1:2),Line2N(1:2,1:2)
    real*8,  intent(out) :: curve_distnce
    
    real*8  :: x1L1,y1L1,x2L1,y2L1
    real*8  :: x1L2,y1L2,x2L2,y2L2
    real*8  :: mL1,mL2,diffInmL,cL1,cL2
    real*8  :: qnt1,qnt2,diffInqnt
    real*8  :: MidPL1(1:2),MidPL2(1:2)
    integer :: CLmaxV
    
    integer :: node1,node2,node3,node4
    integer :: LtSp1,LtSp2,LtSp3
    integer :: PLmethd_or_MPmethd
    real*8  :: distnce_usingPL,distnce_usingMP
    
    x1L1=Line1N(1,1) ; y1L1=Line1N(1,2) ; x2L1=Line1N(2,1) ; y2L1=Line1N(2,2) 
    x1L2=Line2N(1,1) ; y1L2=Line2N(1,2) ; x2L2=Line2N(2,1) ; y2L2=Line2N(2,2) 
    
    if (writeForCurvtr==1) write(*,*) x1L1,y1L1,x2L1,y2L1,"Line1 Two TNs"
    if (writeForCurvtr==1) write(*,*) x1L2,y1L2,x2L2,y2L2,"Line2 Two NIs"
    
    mL1 = (y2L1-y1L1)/(x2L1-x1L1)
    mL2 = (y2L2-y1L2)/(x2L2-x1L2)
    
    if (writeForCurvtr==1) write(*,*) mL1,mL2,"mL1-mL2"
    
    diffInmL = abs((mL1-mL2)/mL1)*100.0000d0
    if (writeForCurvtr==1) write(*,*) mL1,mL2,diffInmL,"mL difference"
    
    if (diffInmL .le. 1.0000d0) then
       PLmethd_or_MPmethd = 1
    elseif (diffInmL .gt. 1.0000d0) then
       write(*,*) 'not exact parallel lines',diffInmL
       PLmethd_or_MPmethd = 2   
    endif
    
    qnt1 = y1L1 - mL1*x1L1 ; qnt2 = y2L1 - mL1*x2L1
    diffInqnt = abs((qnt1-qnt2)/qnt1)*100.0000d0
    if (writeForCurvtr==1) write(*,*) qnt1,qnt2,diffInqnt,"qnt difference"
    if (diffInqnt .gt. 1.0000d0) stop 'not same line for cL1'
    
    cL1 = qnt1 ! %%%%%%%%% c value for Line 1
    
    qnt1 = y1L2 - (mL2)*x1L2 ; qnt2 = y2L2 - (mL2)*x2L2
    diffInqnt = abs((qnt1-qnt2)/qnt1)*100.0000d0
    if (writeForCurvtr==1) write(*,*) qnt1,qnt2,diffInqnt,"qnt difference"
    if (diffInqnt .gt. 1.0000d0) stop 'not same line for cL2'
    
    cL2 = qnt1 ! %%%%%%%%% c value for Line 2
    
    MidPL1(1) = 0.5000d0*(x1L1+x2L1) ; MidPL1(2) = 0.5000d0*(y1L1+y2L1) ! %%% MidPoint of Line1
    MidPL2(1) = 0.5000d0*(x1L2+x2L2) ; MidPL2(2) = 0.5000d0*(y1L2+y2L2) ! %%% MidPoint of Line2 
    
    distnce_usingPL = abs(cL1-cL2)/sqrt(1.0000d0+(mL1**2))
    distnce_usingMP = sqrt((MidPL1(1)-MidPL2(1))**2 + (MidPL1(2)-MidPL2(2))**2)
    
    if (writeForCurvtr==1) write(*,*) distnce_usingPL,distnce_usingMP,"distnce PL/MP"

    CLmaxV = 1000.0000d0
    if (PLmethd_or_MPmethd==1) then
       if (cL1 .gt. CLmaxV) PLmethd_or_MPmethd=2 !check if there is high C value, vertical line
       if (cL1 .le. CLmaxV) PLmethd_or_MPmethd=1
    endif
    
    if (PLmethd_or_MPmethd==1) curve_distnce = distnce_usingPL
    if (PLmethd_or_MPmethd==2) curve_distnce = distnce_usingMP
    
    if (writeForCurvtr==1) write(*,*) curve_distnce,"curve distnce"
    
    if (serialNo==1 .or. serialNo==2) then
       
       LtSp1 = (serialNo-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3) + (NAEC_Apcl+NAEC_Bsal+2) + 1
       LtSp2 = LtSp1+1 ; LtSp3 = LtSp2+1
       
       node1 = spr_node(LtSp1,1) ; node2 = spr_node(LtSp1,2)
       node3 = spr_node(LtSp2,2) ; node4 = spr_node(LtSp3,2)
       
       if (writeForCurvtr==1) then
          write(*,*) "CHEKING FOR SERIAL NO = 1 or 2"
          write(*,*) node1,node2,node3,node4,"nodes"
          write(*,*) node_xy(node1,1:2),node_xy(node4,1:2)
          write(*,*) node_xy(node2,1:2),node_xy(node3,1:2)
       endif
       
    endif
    
  end subroutine curve_dist_measuremnt_using_PLorMP_methd
  
  subroutine measured_curvtrFrm_ExprmntlImg_and_Expctd_curvatr_Height
    implicit none
    integer :: CMv
    real*8  :: ratioMeasuremnt_val
    integer :: i,imax
    
    !%%%%%%%%%% CMv=1 %%%%%%%%%%%%%%%%% CF-WT-Nrt-P-2G9b-m [Jeff Images]
    CMv = 1
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val) ; write(*,*)ratioMeasuremnt_val,"rMvcurvtr"
    
    Exprmntl_CHLtPCS(CMv,1:7)  = +0.0000d0 * ratioMeasuremnt_val; CurvtrSign_LtPCS(CMv,1:7)  =  0
    Exprmntl_CHLtPCS(CMv,8)    = -0.0800d0 * ratioMeasuremnt_val; CurvtrSign_LtPCS(CMv,8)    = -1
    Exprmntl_CHLtPCS(CMv,9)    = -0.0800d0 * ratioMeasuremnt_val; CurvtrSign_LtPCS(CMv,9)    = -1
    Exprmntl_CHLtPCS(CMv,10)   = -0.2200d0 * ratioMeasuremnt_val; CurvtrSign_LtPCS(CMv,10)   = -1
    Exprmntl_CHLtPCS(CMv,11)   = -0.1900d0 * ratioMeasuremnt_val; CurvtrSign_LtPCS(CMv,11)   = -1
    
    Expected_CHLtPCS(CMv,1:Hlf_Ncell) = Unit_Fctr_PCS(CMv) * Exprmntl_CHLtPCS(CMv,1:Hlf_Ncell)
    
    write(*,*) "Curvature Height Printing for CMv =",CMv
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_CHLtPCS(CMv,i)/ratioMeasuremnt_val,Exprmntl_CHLtPCS(CMv,i),Expected_CHLtPCS(CMv,i),i
    enddo
    write(*,*) " " 
    
    
    
    !%%%%%%%%%% CMv=2 %%%%%%%%%%%%%%%%% CF-WT-Nrt-P-2G1b-m [Jeff Images]
    CMv = 2
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val) ; write(*,*)ratioMeasuremnt_val,"rMvcurvtr"
    
    Exprmntl_CHLtPCS(CMv,1:3) = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,1:3) =  0  
    Exprmntl_CHLtPCS(CMv,4)   = -0.0500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,4)   = -1
    Exprmntl_CHLtPCS(CMv,5)   = -0.0500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,5)   = -1
    Exprmntl_CHLtPCS(CMv,6)   = +0.0700d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,6)   = +1
    Exprmntl_CHLtPCS(CMv,7:8) = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,7:8) =  0
    Exprmntl_CHLtPCS(CMv,9)   = -0.1500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,9)   = -1
    Exprmntl_CHLtPCS(CMv,10)  = -0.2700d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,10)  = -1
    Exprmntl_CHLtPCS(CMv,11)  = -0.1600d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,11)  = -1
    
    Expected_CHLtPCS(CMv,1:Hlf_Ncell) = Unit_Fctr_PCS(CMv) * Exprmntl_CHLtPCS(CMv,1:Hlf_Ncell)
    
    write(*,*) "Curvature Height Printing for CMv =",CMv
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_CHLtPCS(CMv,i)/ratioMeasuremnt_val,Exprmntl_CHLtPCS(CMv,i),Expected_CHLtPCS(CMv,i),i
    enddo
    write(*,*) " " 
    
    !%%%%%%%%%% CMv=3 %%%%%%%%%%%%%%%%% CF-WT-Nrt-P-2G8b-m [Jeff Images]
    CMv = 3
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val) ; write(*,*)ratioMeasuremnt_val,"rMvcurvtr"
    
    Exprmntl_CHLtPCS(CMv,1:4) = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,1:4) =   0
    Exprmntl_CHLtPCS(CMv,5)   = +0.0600d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,5)   =  +1
    Exprmntl_CHLtPCS(CMv,6)   = +0.0500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,6)   =  +1
    Exprmntl_CHLtPCS(CMv,7)   = -0.0900d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,7)   =  -1
    Exprmntl_CHLtPCS(CMv,8)   = -0.0800d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,8)   =  -1
    Exprmntl_CHLtPCS(CMv,9)   = -0.1700d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,9)   =  -1
    Exprmntl_CHLtPCS(CMv,10)  = -0.1700d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,10)  =  -1
    Exprmntl_CHLtPCS(CMv,11)  = -0.0800d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,11)  =  -1
    
    Expected_CHLtPCS(CMv,1:Hlf_Ncell) = Unit_Fctr_PCS(CMv) * Exprmntl_CHLtPCS(CMv,1:Hlf_Ncell)
    
    write(*,*) "Curvature Height Printing for CMv =",CMv
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_CHLtPCS(CMv,i)/ratioMeasuremnt_val,Exprmntl_CHLtPCS(CMv,i),Expected_CHLtPCS(CMv,i),i
    enddo
    write(*,*) " "
    
    
    !%%%%%%%%%%%%%%%%%%%% CMv=4 %%%%%%%%%%%%%%%%%%%%% CF-WT-Nrt-P-7G12b-m.tif [Jeff Images]
    CMv = 4
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val) ; write(*,*)ratioMeasuremnt_val,"rMvcurvtr"
    
    Exprmntl_CHLtPCS(CMv,1:2) = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,1:2) =   0
    Exprmntl_CHLtPCS(CMv,3)   = -0.0500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,3)   =  -1
    Exprmntl_CHLtPCS(CMv,4)   = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,4)   =   0
    Exprmntl_CHLtPCS(CMv,5)   = +0.0500d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,5)   =  +1
    Exprmntl_CHLtPCS(CMv,6)   = +0.0000d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,6)   =   0
    Exprmntl_CHLtPCS(CMv,7)   = -0.1400d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,7)   =  -1
    Exprmntl_CHLtPCS(CMv,8)   = -0.3600d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,8)   =  -1
    Exprmntl_CHLtPCS(CMv,9)   = -0.2700d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,9)   =  -1
    Exprmntl_CHLtPCS(CMv,10)  = -0.2900d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,10)  =  -1
    Exprmntl_CHLtPCS(CMv,11)  = -0.0400d0 * ratioMeasuremnt_val ; CurvtrSign_LtPCS(CMv,11)  =  -1
    
    Expected_CHLtPCS(CMv,1:Hlf_Ncell) = Unit_Fctr_PCS(CMv) * Exprmntl_CHLtPCS(CMv,1:Hlf_Ncell)
    
    write(*,*) "Curvature Height Printing for CMv =",CMv
    do i = 1,Hlf_Ncell
       write(*,*) Exprmntl_CHLtPCS(CMv,i)/ratioMeasuremnt_val,Exprmntl_CHLtPCS(CMv,i),Expected_CHLtPCS(CMv,i),i
    enddo
    write(*,*) " "
    
  end subroutine measured_curvtrFrm_ExprmntlImg_and_Expctd_curvatr_Height
  
  subroutine calculate_the_Expected_curve_height(CMv,BsORLt)
    implicit none
    integer, intent(in) :: CMv
    integer, intent(in) :: BsORLt
    integer             :: i,j,jmax
    integer             :: ICv          ! Initiator Cell Value
    
    ICv = Hlf_Ncell
    
    if (BsORLt==1) unitCurvCalcIC_BsPCS(CMv)=Simulted_CHBsPCS(CMv,ICv)/abs(Exprmntl_CHBsPCS(CMv,ICv))
    if (BsORLt==2) unitCurvCalcIC_LtPCS(CMv)=Simulted_CHLtPCS(CMv,ICv)/abs(Exprmntl_CHLtPCS(CMv,ICv))
    
    write(*,*) unitCurvCalcIC_BsPCS(CMv),Exprmntl_CHBsPCS(CMv,ICv),Simulted_CHBsPCS(CMv,ICv),"Bs"
    write(*,*) unitCurvCalcIC_LtPCS(CMv),Exprmntl_CHLtPCS(CMv,ICv),Simulted_CHLtPCS(CMv,ICv),"Lt"
    
    do i = 1,Hlf_Ncell
       
       if (BsORLt==1) Expected_CHBsPCS(CMv,i) = Exprmntl_CHBsPCS(CMv,ICv) * unitCurvCalcIC_BsPCS(CMv)
       if (BsORLt==2) Expected_CHLtPCS(CMv,i) = Exprmntl_CHLtPCS(CMv,ICv) * unitCurvCalcIC_LtPCS(CMv)
       
    enddo
    
  end subroutine calculate_the_Expected_curve_height
  
  
  
  
  subroutine routines_to_fix_CM2_PCS_apprch1(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr)
    integer             :: numOfCellsToAdd
    integer             :: i,j,jmax
    
    call Control_State_match_expected_length_calctns_calls_togthr
    call get_bndryBotm_yNode_frm_CM2PCS ; call Equilibrate_only_NI_model
    
    !%%%%%%%%%%%%%%% HERE I am testing on the CMv=2 and PCS case %%%%%%%%%%%%%%%%%%%
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS) ; call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS) ; call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS) ; call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS) ; call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model ; call print_relevant_info_Actual_and_ExpctdLen(CMv) 
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS) ; call Equilibrate_only_NI_model
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    !call changeNodeY2_of_bndry_frm_4thCell_ltrl_mmbne_Y10
    
    call change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    call Equilibrate_only_NI_model
    call print_relevant_info_Actual_and_ExpctdLen(CMv) !; stop 'tst_complt A0A'
    
    call keep_increasing_A0_to_make_press_within_TOL(CMv)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT TO NOTE THAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! upto this the Frame_NI increases 20 (18 to 37)
    
    call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    
    ! upto this the Frame_NI increases 452 (38 to 489)
    
  end subroutine routines_to_fix_CM2_PCS_apprch1


  subroutine routines_to_fix_CM2_PCS(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    integer             :: numOfCellsToAdd
    integer             :: i,j,jmax
    
    call Control_State_match_expected_length_calctns_calls_togthr
    call get_bndryBotm_yNode_frm_CM2PCS ; call Equilibrate_only_NI_model
    
    !%%%%%%%%%%%%%%% HERE I am testing on the CMv=2 and PCS case %%%%%%%%%%%%%%%%%%%
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    l_Val(1:N_spr)    = l(1:N_spr)
    
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos11"
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos12"
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos13"
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos14"
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos21"
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos22"
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos23"
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL 
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos24"
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos31"
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos32"
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos33"
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL 
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos34"
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos41"
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos42"
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model   ; write(*,*) "Frame_NI = ",(Frame_NI-1),"kIncrTst Pos43"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call change_A0_toCompensate_pressureChng_dueto_l0chngto_EL ;call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call Equilibrate_only_NI_model                                 ;call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call keep_increasing_A0_to_make_press_within_TOL(CMv)
    call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)        ;call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
  end subroutine routines_to_fix_CM2_PCS
  
  
  subroutine routines_to_fix_CM3_PCS(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    integer             :: numOfCellsToAdd,i,j,jmax
    
    call get_bndryBotm_yNode_frm_CM2PCS ; call Equilibrate_only_NI_model
    call Control_State_match_expected_length_calctns_calls_togthr
    
    ! HERE I am testing on the CMv=3 and PCS case.
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    l_Val(1:N_spr)    = l(1:N_spr)
    
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr)+(1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr)+(1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI-1,"Frame_NI aft ks=10 (CM3PCS)"
    
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    write(*,*) Frame_NI-1,"Frame_NI aft set l0"
    
    call change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call Equilibrate_only_NI_model
    
    write(*,*) Frame_NI-1,"Frame_NI aft set l0 2"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call keep_increasing_A0_to_make_press_within_TOL(CMv)
    
    call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
  end subroutine routines_to_fix_CM3_PCS
  
  
  subroutine routines_to_fix_CM1_PCS(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    integer             :: numOfCellsToAdd,i,j,jmax
    
    call get_bndryBotm_yNode_frm_CM2PCS ; call Equilibrate_only_NI_model
    call Control_State_match_expected_length_calctns_calls_togthr
    
    !%%%%%%%%%%%%%%% HERE I am testing on the CMv=1 and PCS case %%%%%%%%%%%%%%%%%%%
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    l_Val(1:N_spr)    = l(1:N_spr)
    
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
    write(*,*) Frame_NI-1,"Frame_NI aft ks=10 (CM1PCS)"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    write(*,*) Frame_NI-1,"Frame_NI aft set_l0_as_Expected_lengths"
    !stop 'test IS finished'
    
    call change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call Equilibrate_only_NI_model
    
    write(*,*) Frame_NI-1,"Frame_NI aft change_A0_toCompensate_pressureChng_dueto_l0chngto_EL"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call keep_increasing_A0_to_make_press_within_TOL(CMv)
    
    write(*,*) Frame_NI-1,"Frame_NI aft keep_increasing_A0_to_make_press_within_TOL"
    
    call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    call manplt_CM1PCS_to_reduce_buckling
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
  end subroutine routines_to_fix_CM1_PCS
  
  
  
  subroutine routines_to_fix_CM4_PCS(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    integer             :: numOfCellsToAdd,i,j,jmax
    
    call get_bndryBotm_yNode_frm_CM2PCS ; call Equilibrate_only_NI_model
    call Control_State_match_expected_length_calctns_calls_togthr
    
    ! %%%%%%%%%%%%%%% HERE I am testing on the CMv=4 and PCS case %%%%%%%%%%%%%%%%%%%
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    l_Val(1:N_spr)    = l(1:N_spr)
    
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    l_Val(1:N_spr) = l(1:N_spr)
    
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_Val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
    write(*,*) Frame_NI-1,"Frame_NI aft ks=10 (CM4PCS)"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)
    call Equilibrate_only_NI_model
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    write(*,*) Frame_NI-1,"Frame_NI aft set_l0_as_Expected_lengths"
    
    call change_A0_toCompensate_pressureChng_dueto_l0chngto_EL
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call Equilibrate_only_NI_model
    write(*,*)Frame_NI-1,"Frame_NI aft change_A0_toCompensate_pressureChng_dueto_l0chngto_EL"
    
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    call keep_increasing_A0_to_make_press_within_TOL(CMv)
    write(*,*) Frame_NI-1,"Frame_NI aft keep_increasing_A0_to_make_press_within_TOL"
    
    call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
  end subroutine routines_to_fix_CM4_PCS
  
  
  subroutine routines_to_fix_CM4_PCS_wt_help_of_CM3_PCS(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    integer             :: numOfCellsToAdd,i,j,jmax
    integer             :: CMvToRead
    
    write(*,*) Frame_NI,"Frame_NI pnt 1"
    call get_bndryBotm_yNode_frm_CM2PCS                               ; call Equilibrate_only_NI_model
    call Control_State_match_expected_length_calctns_calls_togthr
    
    CMvToRead=CMv-1
    call read_the_prps_of_High_ks_tst_wo_nodePrps(CMvToRead,PCSorICS) ; call Equilibrate_only_NI_model
    call save_DiffA0A_set_l0_as_Expected_lengths(CMv,PCSorICS)        ; call Equilibrate_only_NI_model
    
    stop 'routines_to_fix_CM4_PCS_wt_help_of_CM3_PCS' 
    
  end subroutine routines_to_fix_CM4_PCS_wt_help_of_CM3_PCS
  
  
  subroutine High_ksTst_onExistingCondition()
    implicit none
    real*8 :: k_sprVal(1:N_spr)
    real*8 :: l0_Val(1:N_spr)
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr) ! STORING CURRENT SPR VALUE 
    l0_Val(1:N_spr)   = l0(1:N_spr)
    
    ksIncrFctr     = 2.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    !l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    !l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    
    ksIncrFctr     = 4.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    !l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    !l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    
    ksIncrFctr     = 8.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    !l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
    k_spr(1:N_spr) = k_sprVal(1:N_spr) ! GOING BACK TO PRV VAL
    !l0(1:N_spr)    = l0_Val(1:N_spr)   ! GOING BACK TO PRV VAL
    call Equilibrate_only_NI_model
    
    ksIncrFctr     = 10.0d0
    k_spr(1:N_spr) = ksIncrFctr*k_sprVal(1:N_spr)
    !l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
  end subroutine High_ksTst_onExistingCondition
  
  subroutine storing_currnt_sprPrps
    implicit none
    
  end subroutine storing_currnt_sprPrps
  
  subroutine change_kspr_Equil_change_l0_Equil_and_goback(ksIncrFctr)
    implicit none
    real*8, intent(in) :: ksIncrFctr
    real*8              :: k_sprVal(1:N_spr),l0_Val(1:N_spr),l_Val(1:N_spr)
    
    k_sprVal(1:N_spr) = k_spr(1:N_spr)
    l0_Val(1:N_spr)   = l0(1:N_spr)
    l_Val(1:N_spr)    = l(1:N_spr)
    
    k_spr(1:N_spr) = (ksIncrFctr)*k_sprVal(1:N_spr)
    call Equilibrate_only_NI_model
    l0(1:N_spr)    = (1.0d0-(1.0d0/ksIncrFctr))*l_val(1:N_spr) + (1.0d0/ksIncrFctr)*(l0_Val(1:N_spr))
    call Equilibrate_only_NI_model
    
  end subroutine change_kspr_Equil_change_l0_Equil_and_goback
  
  
  subroutine save_the_prps_of_High_ks_tst(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    character(len=100)  :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer             :: i,imax,j,jmax
    real*8              :: TnsnV,ApLength,BsLength,LtLength,Apks,Bsks,Ltks,ApL0,BsL0,LtL0,TnsnAp,TnsnBs,TnsnLt
    integer             :: sprAp,sprBs1,sprBs2,sprBs3,sprLt1,sprLt2,sprLt3
    real*8              :: PressV
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=232,file=trim(adjustl(wflnmSprs)))
    open(unit=233,file=trim(adjustl(wflnmConv)))
    open(unit=234,file=trim(adjustl(wflnmCell)))
    open(unit=235,file=trim(adjustl(wflnmNode)))
    
    imax = 4
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = Hlf_Ncell
       if (i==3) jmax = N_cell
       if (i==4) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             TnsnV = k_spr(j)*(l0(j)-l(j))
             write(232,*) TnsnV,k_spr(j),l(j),l0(j),j
             
          elseif (i==2) then
             
             sprAp  = (j-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)+1
             sprBs1 = sprAp+1 ; sprBs2 = sprAp+2 ; sprBs3 = sprAp+3
             sprLt1 = sprAp+4 ; sprLt2 = sprAp+5 ; sprLt3 = sprAp+6
             
             ApLength = l(sprAp)
             BsLength = l(sprBs1)+l(sprBs2)+l(sprBs3)
             LtLength = l(sprLt1)+l(sprLt2)+l(sprLt3)
             
             Apks = k_spr(sprAp)
             Bsks = 1.0d0/((1.0d0/k_spr(sprBs1)) + (1.0d0/k_spr(sprBs2)) + (1.0d0/k_spr(sprBs3)))
             Ltks = 1.0d0/((1.0d0/k_spr(sprLt1)) + (1.0d0/k_spr(sprLt2)) + (1.0d0/k_spr(sprLt3)))
             
             ApL0 = l0(sprAp)
             BsL0 = l0(sprBs1)+l0(sprBs2)+l0(sprBs3)
             LtL0 = l0(sprLt1)+l0(sprLt2)+l0(sprLt3)

             TnsnAp = k_spr(sprAp) *(l0(sprAp) -l(sprAp))
             TnsnBs = k_spr(sprBs1)*(l0(sprBs1)-l(sprBs1))
             TnsnLt = k_spr(sprLt1)*(l0(sprLt1)-l(sprLt1))

             write(233,*) j,TnsnAp,TnsnBs,TnsnLt,Apks,Bsks,Ltks,ApLength,BsLength,LtLength,ApL0,BsL0,LtL0
             
          elseif (i==3) then
             
             PressV = k_area(j)*(A0(j)-A(j))
             write(234,*) j,PressV,k_area(j),A(j),A0(j)

          elseif (i==4) then 
             write(235,*) j,node_typ(j),CgXNode(j),CgYNode(j),node_xy(j,1),node_xy(j,2)
          endif
          
       enddo
       
    enddo
    
    
    close(232)
    close(233)
    close(234)
    close(235)
    
  end subroutine save_the_prps_of_High_ks_tst
  
  
  subroutine save_the_prps_of_High_ks_tst_buffr(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    character(len=100)  :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer             :: i,imax,j,jmax
    real*8              :: TnsnV,ApLength,BsLength,LtLength,Apks,Bsks,Ltks,ApL0,BsL0,LtL0,TnsnAp,TnsnBs,TnsnLt
    integer             :: sprAp,sprBs1,sprBs2,sprBs3,sprLt1,sprLt2,sprLt3
    real*8              :: PressV
    
    call flnm_management_for_High_ks_tst_buffr(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
  
    
    open(unit=232,file=trim(adjustl(wflnmSprs)))
    open(unit=233,file=trim(adjustl(wflnmConv)))
    open(unit=234,file=trim(adjustl(wflnmCell)))
    open(unit=235,file=trim(adjustl(wflnmNode)))
    
    imax = 4
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = Hlf_Ncell
       if (i==3) jmax = N_cell
       if (i==4) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             
             TnsnV = k_spr(j)*(l0(j)-l(j))
             write(232,*) TnsnV,k_spr(j),l(j),l0(j),j
             
          elseif (i==2) then
             
             sprAp  = (j-1)*(NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)+1
             sprBs1 = sprAp+1 ; sprBs2 = sprAp+2 ; sprBs3 = sprAp+3
             sprLt1 = sprAp+4 ; sprLt2 = sprAp+5 ; sprLt3 = sprAp+6
             
             ApLength = l(sprAp)
             BsLength = l(sprBs1)+l(sprBs2)+l(sprBs3)
             LtLength = l(sprLt1)+l(sprLt2)+l(sprLt3)
             
             Apks = k_spr(sprAp)
             Bsks = 1.0d0/((1.0d0/k_spr(sprBs1)) + (1.0d0/k_spr(sprBs2)) + (1.0d0/k_spr(sprBs3)))
             Ltks = 1.0d0/((1.0d0/k_spr(sprLt1)) + (1.0d0/k_spr(sprLt2)) + (1.0d0/k_spr(sprLt3)))
             
             ApL0 = l0(sprAp)
             BsL0 = l0(sprBs1)+l0(sprBs2)+l0(sprBs3)
             LtL0 = l0(sprLt1)+l0(sprLt2)+l0(sprLt3)

             TnsnAp = k_spr(sprAp) *(l0(sprAp) -l(sprAp))
             TnsnBs = k_spr(sprBs1)*(l0(sprBs1)-l(sprBs1))
             TnsnLt = k_spr(sprLt1)*(l0(sprLt1)-l(sprLt1))

             write(233,*) j,TnsnAp,TnsnBs,TnsnLt,Apks,Bsks,Ltks,ApLength,BsLength,LtLength,ApL0,BsL0,LtL0

          elseif (i==3) then
             
             PressV = k_area(j)*(A0(j)-A(j))
             write(234,*) j,PressV,k_area(j),A(j),A0(j)

          elseif (i==4) then 
             write(235,*) j,node_typ(j),CgXNode(j),CgYNode(j),node_xy(j,1),node_xy(j,2)
          endif
          
       enddo
       
    enddo
    
    
    close(232)
    close(233)
    close(234)
    close(235)
    
  end subroutine save_the_prps_of_High_ks_tst_buffr
  
  subroutine read_the_prps_of_High_ks_tst(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    character(len=100)  :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer             :: i,imax,j,jmax
    real*8              :: TnsnV,PressV
    integer             :: sprNm,cellNm,nodeNm
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=532,file=trim(adjustl(wflnmSprs)))
    open(unit=534,file=trim(adjustl(wflnmCell)))
    open(unit=535,file=trim(adjustl(wflnmNode)))
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) read(532,*) TnsnV,k_spr(j),l(j),l0(j),sprNm
          if (i==2) read(534,*) cellNm,PressV,k_area(j),A(j),A0(j)   
          if (i==3) read(535,*) nodeNm,node_typ(j),CgXNode(j),CgYNode(j),node_xy(j,1),node_xy(j,2)
       enddo
       
    enddo
    
    close(532)
    close(534)
    close(535)
    
  end subroutine read_the_prps_of_High_ks_tst
  
  
  subroutine read_the_prps_of_High_ks_tst_woNodeTYP(CMv,PCSorICS)
    implicit none
    integer, intent(in) :: CMv,PCSorICS
    character(len=100)  :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer             :: i,imax,j,jmax
    real*8              :: TnsnV,PressV
    integer             :: sprNm,cellNm,nodeNm
    real*8              :: nodeXv,nodeYv
    integer             :: nodeTypV
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=532,file=trim(adjustl(wflnmSprs)))
    open(unit=534,file=trim(adjustl(wflnmCell)))
    open(unit=535,file=trim(adjustl(wflnmNode)))
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) read(532,*) TnsnV,k_spr(j),l(j),l0(j),sprNm
          if (i==2) read(534,*) cellNm,PressV,k_area(j),A(j),A0(j)   
          if (i==3) read(535,*) nodeNm,nodeTypV,CgXNode(j),CgYNode(j),nodeXv,nodeYv
       enddo
       
    enddo
    
    close(532)
    close(534)
    close(535)
    
  end subroutine read_the_prps_of_High_ks_tst_woNodeTYP
  
  subroutine read_the_prps_of_High_ks_tst_wo_nodePrps(CMv,PCSorICS)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    character(len=100)   :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer              :: i,imax,j,jmax
    real*8               :: TnsnV,PressV
    integer              :: sprNm,cellNm,nodeNm
    integer              :: currCell,tobeCell,cell_bndry,cell_chnging
    integer              :: nodeTypV,nsprsInACell,sprAdd
    real*8               :: CgXval(1:N_node),CgYval(1:N_node),nodeXv,nodeYv
    real*8               :: lenV,areaV
    real*8               :: ksVal(1:N_spr) ,l0Val(1:N_spr)
    real*8               :: kaVal(1:N_cell),A0Val(1:N_cell)
    integer, allocatable :: currSprs(:),tobeSprs(:),intgrVal(:)
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=532,file=trim(adjustl(wflnmSprs)))
    open(unit=534,file=trim(adjustl(wflnmCell)))
    open(unit=535,file=trim(adjustl(wflnmNode)))
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          
          if (i==1) read(532,*) TnsnV,  ksVal(j), lenV,      l0Val(j),  sprNm
          if (i==2) read(534,*) cellNm, PressV,   kaVal(j),  areaV,     A0Val(j)   
          if (i==3) read(535,*) nodeNm, nodeTypV, CgXval(j), CgYval(j), nodeXv,nodeYv
          
       enddo
       
    enddo
    
    close(532)
    close(534)
    close(535)
    
    k_spr(1:N_spr)    = ksVal(1:N_spr)   ; l0(1:N_spr)       = l0Val(1:N_spr)
    k_area(1:N_cell)  = kaVal(1:N_cell)  ; A0(1:N_cell)      = A0Val(1:N_cell)
    CgXNode(1:N_node) = CgXval(1:N_node) ; CgYNode(1:N_node) = CgYval(1:N_node)
    
    cell_chnging  = Hlf_Ncell-1 ; write(*,*) cell_chnging,"cell_chnging"
    cell_bndry    = 3
    nsprsInACell  = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    sprAdd        = Hlf_Ncell*nsprsInACell
    
    allocate(currSprs(1:nsprsInACell),tobeSprs(1:nsprsInACell),intgrVal(1:nsprsInACell))
    currSprs = -1 ; tobeSprs = -1
    do i = 1,nsprsInACell
       intgrVal(i) = i
    enddo
    
    
    do i = (cell_bndry+1),cell_chnging
       
       currCell = i             ; tobeCell = i+1
       currSprs(1:nsprsInACell) = (currCell-1)*(nsprsInACell) + intgrVal(1:nsprsInACell)
       tobeSprs(1:nsprsInACell) = (tobeCell-1)*(nsprsInACell) + intgrVal(1:nsprsInACell)
       
       A0(currCell)                              = A0Val(tobeCell)
       k_area(currCell)                          = kaVal(tobeCell)
       l0(currSprs(1):currSprs(nsprsInACell))    = l0(tobeSprs(1):tobeSprs(nsprsInACell))
       k_spr(currSprs(1):currSprs(nsprsInACell)) = k_spr(tobeSprs(1):tobeSprs(nsprsInACell))
       
       A0(currCell+Hlf_Ncell)                                      = A0(currCell)
       k_area(currCell+Hlf_Ncell)                                  = k_area(currCell)
       l0((currSprs(1)+sprAdd):(currSprs(nsprsInACell)+sprAdd))    = l0(currSprs(1):currSprs(nsprsInACell))
       k_spr((currSprs(1)+sprAdd):(currSprs(nsprsInACell)+sprAdd)) = k_spr(currSprs(1):currSprs(nsprsInACell))
       
    enddo
    
    !open(unit=987,file='prprty_var.dat')
    
    !do i = (cell_bndry+1),cell_chnging
       
    !   currCell = i             ; tobeCell = i+1
    !   currSprs(1:nsprsInACell) = (currCell-1)*(nsprsInACell) + intgrVal(1:nsprsInACell)
    !   tobeSprs(1:nsprsInACell) = (tobeCell-1)*(nsprsInACell) + intgrVal(1:nsprsInACell)
       
    !   write(987,*) A0(currCell),A0Val(tobeCell),currCell,tobeCell,
    !enddo
    
    !close(987)
    
  end subroutine read_the_prps_of_High_ks_tst_wo_nodePrps
  
  
  subroutine read_the_prps_of_High_ks_tst_PCSorICS(CMv,PCSorICS,ksVal,l0Val,kaVal,A0Val,CgXval,CgYval)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    real*8 , intent(out) :: ksVal(1:N_spr)   , l0Val(1:N_spr)
    real*8 , intent(out) :: kaVal(1:N_cell)  , A0Val(1:N_cell)
    real*8 , intent(out) :: CgxVal(1:N_node) , CgyVal(1:N_node)
    
    character(len=100)   :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer              :: i,imax,j,jmax
    real*8               :: TnsnV,PressV
    integer              :: sprNm,cellNm,nodeNm
    integer              :: nodeTypV,nsprsInACell,sprAdd
    real*8               :: nodeXv,nodeYv
    real*8               :: lenV,areaV
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=632,file=trim(adjustl(wflnmSprs)))
    open(unit=634,file=trim(adjustl(wflnmCell)))
    open(unit=635,file=trim(adjustl(wflnmNode)))
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          
          if (i==1) read(632,*) TnsnV,  ksVal(j), lenV,      l0Val(j),  sprNm
          if (i==2) read(634,*) cellNm, PressV,   kaVal(j),  areaV,     A0Val(j)   
          if (i==3) read(635,*) nodeNm, nodeTypV, CgxVal(j), CgyVal(j), nodeXv,nodeYv
          
       enddo
       
    enddo
    
    close(632)
    close(634)
    close(635)
    
  end subroutine read_the_prps_of_High_ks_tst_PCSorICS
  
  subroutine read_the_nodePrps_frm_PCS_data(CMv,PCSorICS)
    implicit none
    integer, intent(in)  :: CMv,PCSorICS
    character(len=100)   :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer              :: i,imax,j,jmax
    real*8               :: TnsnV,PressV
    integer              :: sprNm,cellNm,nodeNm
    integer              :: nodeTypV,nsprsInACell,sprAdd
    real*8               :: nodeXv,nodeYv
    real*8               :: lenV,areaV
    
    call flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=837,file=trim(adjustl(wflnmNode)))
    do i = 1,N_node     
       read(837,*) nodeNm,node_typ(i),CgXNode(i),CgYNode(i),node_xy(i,1),node_xy(i,2)
    enddo
    close(837)
    
  end subroutine read_the_nodePrps_frm_PCS_data
  
  
  subroutine read_TIIFfile_props(CMv,TIIF_ValIndctr)
    implicit none
    integer, intent(in) :: CMv,TIIF_ValIndctr
    character(len=100)  :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    integer             :: i,imax,j,jmax
    real*8              :: TnsnV,PressV
    integer             :: sprNm,cellNm,nodeNm
    
    call flnm_management_forTIIFfile(CMv,TIIF_ValIndctr,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    
    open(unit=532,file=trim(adjustl(wflnmSprs)))
    open(unit=534,file=trim(adjustl(wflnmCell)))
    open(unit=535,file=trim(adjustl(wflnmNode)))
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) read(532,*) TnsnV,k_spr(j),l(j),l0(j),sprNm
          if (i==2) read(534,*) cellNm,PressV,k_area(j),A(j),A0(j)   
          if (i==3) read(535,*) nodeNm,node_typ(j),CgXNode(j),CgYNode(j),node_xy(j,1),node_xy(j,2)
       enddo
       
    enddo
    
    close(532)
    close(534)
    close(535)
    
  end subroutine read_TIIFfile_props
  
  
  subroutine flnm_management_for_High_ks_tst(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    implicit none
    integer, intent(in)             :: CMv,PCSorICS
    character(len=100), intent(out) :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    character(len=100)              :: wflnmSprsPrfx,wflnmConvPrfx,wflnmCellPrfx,wflnmNodePrfx,wflnmSufx
    
    write(wflnmSprsPrfx,'(a)') 'sprsPrps_CM'
    write(wflnmConvPrfx,'(a)') 'convSprs_CM'
    write(wflnmCellPrfx,'(a)') 'cellPrps_CM'
    write(wflnmNodePrfx,'(a)') 'nodePrps_CM'
    
    if (PCSorICS==1) write(wflnmSufx,'(i1.1,a)') CMv,'PCS_HighKs.dat'
    if (PCSorICS==2) write(wflnmSufx,'(i1.1,a)') CMv,'ICS_HighKs.dat'
    
    write(*,*) trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
    wflnmSprs=trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    wflnmConv=trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    wflnmCell=trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    wflnmNode=trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
  end subroutine flnm_management_for_High_ks_tst
  
  subroutine flnm_management_for_High_ks_tst_buffr(CMv,PCSorICS,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    implicit none
    integer, intent(in)             :: CMv,PCSorICS
    character(len=100), intent(out) :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    character(len=100)              :: wflnmSprsPrfx,wflnmConvPrfx,wflnmCellPrfx,wflnmNodePrfx,wflnmSufx
    
    write(wflnmSprsPrfx,'(a)') 'sprsPrps_CM'
    write(wflnmConvPrfx,'(a)') 'convSprs_CM'
    write(wflnmCellPrfx,'(a)') 'cellPrps_CM'
    write(wflnmNodePrfx,'(a)') 'nodePrps_CM'
    
    if (PCSorICS==1) write(wflnmSufx,'(i1.1,a)') CMv,'PCS_HighKs_buffr.dat'
    if (PCSorICS==2) write(wflnmSufx,'(i1.1,a)') CMv,'ICS_HighKs_buffr.dat'
    
    write(*,*) trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
    wflnmSprs=trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    wflnmConv=trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    wflnmCell=trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    wflnmNode=trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
  end subroutine flnm_management_for_High_ks_tst_buffr
  
  subroutine flnm_management_forTIIFfile(CMv,TIIF_ValIndctr,wflnmSprs,wflnmConv,wflnmCell,wflnmNode)
    implicit none
    integer, intent(in)             :: CMv,TIIF_ValIndctr
    character(len=100), intent(out) :: wflnmSprs,wflnmConv,wflnmCell,wflnmNode
    character(len=100)              :: wflnmSprsPrfx,wflnmConvPrfx,wflnmCellPrfx,wflnmNodePrfx,wflnmSufx
    
    write(wflnmSprsPrfx,'(a)') 'sprsPrps_CM'
    write(wflnmConvPrfx,'(a)') 'convSprs_CM'
    write(wflnmCellPrfx,'(a)') 'cellPrps_CM'
    write(wflnmNodePrfx,'(a)') 'nodePrps_CM'
    
    if (TIIF_ValIndctr==0) write(wflnmSufx,'(i1.1,a)') CMv,'PCS_WOPull_TIIF_000.dat'
    if (TIIF_ValIndctr==1) write(wflnmSufx,'(i1.1,a)') CMv,'PCS_WOPull_TIIF_025.dat'
    if (TIIF_ValIndctr==2) write(wflnmSufx,'(i1.1,a)') CMv,'PCS_WOPull_TIIF_050.dat'
    
    write(*,*) trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    write(*,*) trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
    wflnmSprs=trim(adjustl(wflnmSprsPrfx))//trim(adjustl(wflnmSufx))
    wflnmConv=trim(adjustl(wflnmConvPrfx))//trim(adjustl(wflnmSufx))
    wflnmCell=trim(adjustl(wflnmCellPrfx))//trim(adjustl(wflnmSufx))
    wflnmNode=trim(adjustl(wflnmNodePrfx))//trim(adjustl(wflnmSufx))
    
  end subroutine flnm_management_forTIIFfile
  
  subroutine Y_force_perturbation_to_PCS
    implicit none
    integer :: forceNodeL,forceNodeR
    integer :: i,j,jmax
    real*8  :: forceNodeYval
    
    
    forceNodeL    = (2*Hlf_Ncell)+1       ; forceNodeR = forceNodeL+((Hlf_Ncell+1)*2)
    forceNodeYval = node_xy(forceNodeL,2) 
    
    do i = 1,11
       
       CgYNode(forceNodeL) = CgYNode(forceNodeL) - 0.0500d0
       CgYNode(forceNodeR) = CgYNode(forceNodeL)
       
       if (abs(CgYNode(forceNodeL)) .lt. 1.00d-15) then
          CgYNode(forceNodeL) = 0.000000000000000d0
          CgYNode(forceNodeR) = CgYNode(forceNodeL)
       endif
       
       call Equilibrate_only_NI_model
       
    enddo
    
    
    do i = 1,15
       
       CgYNode(forceNodeL) = CgYNode(forceNodeL) + 0.0500d0
       CgYNode(forceNodeR) = CgYNode(forceNodeL)
       
       if (abs(CgYNode(forceNodeL)) .lt. 1.00d-15) then
          CgYNode(forceNodeL) = 0.000000000000000d0
          CgYNode(forceNodeR) = CgYNode(forceNodeL)
       endif
       
       call Equilibrate_only_NI_model
       
    enddo
    
  end subroutine Y_force_perturbation_to_PCS
  
  
  subroutine SensitvtyTst_vrtclForce_onTrmnlNode_ofIC_Apcl
    implicit none
    real*8  :: ksStr(1:N_spr),l0Str(1:N_spr),kaStr(1:N_cell),A0Str(1:N_cell)
    real*8  :: CgXStr(1:N_node),CgYStr(1:N_node)
    real*8  :: lvPos1(1:N_spr),lvPos2(1:N_spr),lvPos3(1:N_spr)
    real*8  :: PressPos1(1:N_cell),PressPos2(1:N_cell),PressPos3(1:N_cell)
    real*8  :: curveLtHPos1(1:Hlf_Ncell),curveLtHPos2(1:Hlf_Ncell),curveLtHPos3(1:Hlf_Ncell)
    real*8  :: bndryTiltPos1,ydisPos1,bndryTiltPos2,ydisPos2,bndryTiltPos3,ydisPos3
    integer :: i,j
    
    ksStr=k_spr ; l0Str=l0 ; kaStr=k_area ; A0Str=A0 ; CgXStr=CgXNode ; CgYStr=CgYNode
    
    call get_length_press_and_LtrlCurvtr(lvPos1,PressPos1,curveLtHPos1)
    call get_bndryTilt_and_yDisplcmnt_ofTrmnlApclNodeIC(bndryTiltPos1,ydisPos1)
    
    call  remveForceFrmApclSrfcOfIC_frmCurrntVal()
    
    call get_length_press_and_LtrlCurvtr(lvPos2,PressPos2,curveLtHPos2)
    call get_bndryTilt_and_yDisplcmnt_ofTrmnlApclNodeIC(bndryTiltPos2,ydisPos2)
    call write_chnges_SensitvtyTst_vrtclForce(lvPos1,lvPos2,PressPos1,PressPos2,curveLtHPos1,curveLtHPos2,&
         bndryTiltPos1,bndryTiltPos2,ydisPos1,ydisPos2)
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    call get_length_press_and_LtrlCurvtr(lvPos3,PressPos3,curveLtHPos3)
    call get_bndryTilt_and_yDisplcmnt_ofTrmnlApclNodeIC(bndryTiltPos3,ydisPos3)
    call write_chnges_SensitvtyTst_vrtclForce(lvPos2,lvPos3,PressPos2,PressPos3,curveLtHPos2,&
         curveLtHPos3,bndryTiltPos2,bndryTiltPos3,ydisPos2,ydisPos3)
    call write_chnges_SensitvtyTst_vrtclForce(lvPos1,lvPos3,PressPos1,PressPos3,curveLtHPos1,curveLtHPos3,&
         bndryTiltPos1,bndryTiltPos3,ydisPos1,ydisPos3)
    
    k_spr=ksStr ; l0=l0Str ; k_area=kaStr ; A0=A0Str ; CgXNode=CgXStr ; CgYNode=CgYStr
    
  end subroutine SensitvtyTst_vrtclForce_onTrmnlNode_ofIC_Apcl
  
  
  subroutine get_length_press_and_LtrlCurvtr(lengthV,PressV,curve_LtH)
    implicit none
    real*8, intent(out) :: lengthV(1:N_spr),PressV(1:N_cell)
    real*8 ,intent(out) :: curve_LtH(1:Hlf_Ncell)
    integer             :: i,j,BsOrLt

    interface
       function Pressure(dum_Nodes,dumA0)
         use system_parameters
         use transfrm_info
         implicit none
         real*8 :: Pressure(1:N_cell)
         real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8 :: dumA0(1:N_cell)
       end function Pressure
    end interface
    
    lengthV(1:N_spr) = l(1:N_spr)
    PressV(1:N_cell) = Pressure(node_xy,A0)
    BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,curve_LtH)
    
  end subroutine get_length_press_and_LtrlCurvtr
  
  subroutine get_bndryTilt_and_yDisplcmnt_ofTrmnlApclNodeIC(bndryTilt,ydis)
    implicit none
    real*8,  intent(out) :: bndryTilt,ydis
    integer              :: FrstJointNodeL,FrstJointNodeR
    
    FrstJointNodeL = (2*Hlf_Ncell)+1
    FrstJointNodeR = FrstJointNodeL + 2*(Hlf_Ncell+1); write(*,*) FrstJointNodeL,FrstJointNodeR,"L/R"
    
    bndryTilt = abs(node_xy(1,1)-node_xy(2,1))
    ydis      = 0.50d0 * abs(node_xy(FrstJointNodeL,2)+node_xy(FrstJointNodeR,2))
    
  end subroutine get_bndryTilt_and_yDisplcmnt_ofTrmnlApclNodeIC
  
  subroutine write_chnges_SensitvtyTst_vrtclForce(lvBfr,lvAft,PressVBfr,PressVAft,curveLtHBfr,curveLtHAft,&
       bndryTiltPosBfr,bndryTiltPosAft,ydisBfr,ydisAft)
    implicit none
    real*8, intent(in) :: lvBfr(1:N_spr),lvAft(1:N_spr)
    real*8, intent(in) :: PressVBfr(1:N_cell),PressVAft(1:N_cell)
    real*8, intent(in) :: curveLtHBfr(1:N_cell),curveLtHAft(1:N_cell)
    real*8, intent(in) :: bndryTiltPosBfr,bndryTiltPosAft,ydisBfr,ydisAft
    integer            :: i,j,jmax
    
    open(unit=672,file='Sensitvty_vrtclForce.dat',position='append')
    
    write(672,*) "CellsMeet value =",CellsMeet
    write(672,*)  bndryTiltPosBfr,bndryTiltPosAft,ydisBfr,ydisAft,"bndryTile B/A and ydis B/A"
    write(672,*) " "
    
    do i = 1,3
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = Hlf_Ncell
       
       do j = 1,jmax
          if (i==1) write(672,*) j,abs((lvAft(j)-lvBfr(j))/lvBfr(j))*100.00d0,"lc"
          if (i==2) write(672,*) j,abs((PressVAft(j)-PressVBfr(j))/PressVBfr(j))*100.00d0,"pc"
          if (i==3) write(672,*) j,abs((curveLtHAft(j)-curveLtHBfr(j))/curveLtHBfr(j))*100.00d0,"cc"
       enddo
       write(672,*) " "
    enddo
    write(672,*) " "
    close(672)
    
  end subroutine write_chnges_SensitvtyTst_vrtclForce
  
  subroutine manplt_CM1PCS_to_reduce_buckling
    implicit none
    integer              :: BsORLt,i,j,cell_val
    real*8               :: CurrCurvtr(1:Hlf_Ncell),ExpctdCurvtr(1:Hlf_Ncell)
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    real*8               :: hmL0Lt
    integer              :: PCSorICSval
    
    write(*,*) CellsMeet,"CellsMeet"
    if (CellsMeet.ne.1) stop 'CellsMeet not equal to 1'
    
    BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CellsMeet,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:NAEC_Bsal+1),LtrlSp(1:(NAEC_Ltrl+1)))
    ApclSp=-1 ; BsalSp=-1 ; LtrlSp=-1
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp,BsalSp,LtrlSp)
    
    hmL0Lt=0.9900d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
    call Equilibrate_only_NI_model
    hmL0Lt=0.9900d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
    call Equilibrate_only_NI_model
    
    writeForCurvtr=1 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine manplt_CM1PCS_to_reduce_buckling
  
  subroutine read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    implicit none
    integer, intent(in)  :: CMv
    real*8 , intent(out) :: ExpctdCurvtr(1:Hlf_Ncell)
    
    real*8  :: ratioV
    integer :: i,j,cellV
    integer :: CMvRead
    real*8  :: ExprmntlCurvtrOrg(1:Hlf_Ncell),ExprmntlCurvtrConv(1:Hlf_Ncell)
    
    open(unit=431,file='CurvatureInfoRead.dat')
    
    do 
       read(431,*) ratioV
       read(431,*) CMvRead
       
       do i = 1,Hlf_Ncell
          read(431,*) ExprmntlCurvtrOrg(i),ExprmntlCurvtrConv(i),ExpctdCurvtr(i),cellV
       enddo
       
       if (CMvRead==CMv) exit
    enddo
    
    close(431)
    
  end subroutine read_expctd_curvtr_height
  
  
  subroutine read_expctd_curvtr_height_and_sign(CMv,ExpctdCurvtr,CurvtrSignLt)
    implicit none
    integer, intent(in)  :: CMv
    real*8 , intent(out) :: ExpctdCurvtr(1:Hlf_Ncell)
    integer, intent(out) :: CurvtrSignLt(1:Hlf_Ncell)
    
    real*8  :: ratioV
    integer :: i,j,cellV
    integer :: CMvRead
    real*8  :: ExprmntlCurvtrOrg(1:Hlf_Ncell),ExprmntlCurvtrConv(1:Hlf_Ncell)
    
    
    open(unit=431,file='CurvatureInfoWithSignRead.dat')
    
    do 
       read(431,*) ratioV
       read(431,*) CMvRead
       
       do i = 1,Hlf_Ncell
          read(431,*) ExprmntlCurvtrOrg(i),ExprmntlCurvtrConv(i),ExpctdCurvtr(i),CurvtrSignLt(i),cellV
       enddo
       
       if (CMvRead==CMv) exit
    enddo
    
    close(431)
    
  end subroutine read_expctd_curvtr_height_and_sign
  
  
  subroutine compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    implicit none
    real*8, intent(in) :: CurrCurvtr(1:Hlf_Ncell),ExpctdCurvtr(1:Hlf_Ncell)
    integer            :: i,j
    real*8             :: CurvtrDeviatn(1:Hlf_Ncell)
    real*8             :: TINY_val
    
    TINY_val = 0.00000000001d0
    
    open(unit=434,file='CurvatureComprsn.dat',position='append')
    
    do i = 1,Hlf_Ncell
       if (abs(ExpctdCurvtr(i)).le.TINY_val) then
          CurvtrDeviatn(i) = 0.0000d0
       elseif (abs(ExpctdCurvtr(i)).gt.TINY_val) then
          CurvtrDeviatn(i) = abs((abs(CurrCurvtr(i))-abs(ExpctdCurvtr(i)))/ExpctdCurvtr(i))*100.00d0
       endif
       write(434,*) CurrCurvtr(i),ExpctdCurvtr(i),CurvtrDeviatn(i),i,CellsMeet
    enddo
    
    write(434,*) " " 
    close(434)
    
  end subroutine compare_curvtrs_writing
  
  subroutine read_TIIFfile_props_and_apply_force(CMv)
    implicit none
    integer, intent(in) :: CMv
    integer             :: TIIF_valIndctr,pullpushIndctr
    real*8              :: CgXChng
    
    TIIF_ValIndctr=0 ; call apply_force_on_diff_TIIF_PCS(CMv,TIIF_ValIndctr)
    TIIF_ValIndctr=1 ; call apply_force_on_diff_TIIF_PCS(CMv,TIIF_ValIndctr)
    TIIF_ValIndctr=2 ; call apply_force_on_diff_TIIF_PCS(CMv,TIIF_ValIndctr)
    
  end subroutine read_TIIFfile_props_and_apply_force
  
  subroutine apply_force_on_diff_TIIF_PCS(CMv,TIIF_ValIndctr)
    implicit none
    integer, intent(in) :: CMv,TIIF_ValIndctr
    integer             :: pullpushIndctr
    real*8              :: CgXChng
    
    call read_TIIFfile_props(CMv,TIIF_ValIndctr)
    call Equilibrate_only_NI_model
    
    pullpushIndctr = 1                   ! 1 for pull, -1 for push
    
    CgXChng=0.125d0 ; call pull_push_atBndry(pullpushIndctr,CgXChng) ; write(*,*) Frame_NI-1,"L11"
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()     ; write(*,*) Frame_NI-1,"L21"
    CgXChng=0.125d0 ; call pull_push_atBndry(pullpushIndctr,CgXChng) ; write(*,*) Frame_NI-1,"L12"
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()     ; write(*,*) Frame_NI-1,"L22"
    
    
  end subroutine apply_force_on_diff_TIIF_PCS
  
  
  subroutine make_equivalent_TIIF_for_all_CMPCS(CMv)
    implicit none
    integer, intent(in) :: CMv
    
    integer :: BsOrLt
    real*8  :: ExpctdCurvtr(1:Hlf_Ncell),CurrCurvtr(1:Hlf_Ncell)
    real*8  :: applied_TIIF
    real*8  :: d1v,d2v,d3v
    integer :: nodeL,nodeR
    real*8  :: diff_TIIF,TIIF_val,TINY_val
    
    write(*,*) "make_equivalent_TIIF CMv value =",CMv
    
    TINY_val = 0.00000001d0
    
    call get_the_shifting_nodes(nodeL,nodeR)
    d1v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce at beginning"
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    call read_TIIF_value(applied_TIIF)
    diff_TIIF = abs(applied_TIIF-intended_TIIF)
    TIIF_val  = intended_TIIF
    
    if (diff_TIIF.gt.TINY_val)call remveOrAdd_ForceAtApclSrfcOfIC_ToMaintainSpcfcVal(TIIF_val)
    
    d2v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft intnded"
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"loc2"
    d2v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft bndry adjst"
    
    writeForCurvtr=0 ; BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine make_equivalent_TIIF_for_all_CMPCS
  
  subroutine make_InptTIIF_for_all_CMPCS_frmZeroTIIF(CMv,TIIF_val)
    implicit none
    integer, intent(in) :: CMv
    real*8,  intent(in) :: TIIF_val
    
    real*8  :: TINY_val,TIIF_curr,IncrAmnt
    integer :: nodeL,nodeR,N_diffForces
    
    
    TINY_val = 0.00000001d0
    
    call get_the_shifting_nodes(nodeL,nodeR) ; write(*,*) nodeL,nodeR,"node L/R in sb:make_InptTIIF.....ZeroTIIF"
    TIIF_curr = CgYNode(nodeL) ; write(*,*) TIIF_curr,"TIIF_curr"
    
    if (TIIF_curr.gt.TINY_val) then
       write(*,*) "TIIF value is not ZERO and =",TIIF_curr
       stop 'inconsistent TIIF'
    endif
    
    
    N_diffForces = 10 ; IncrAmnt = 0.0500d0
    call applyForceApclSurfcOfIC_WT_BndryAdjst_bfr_adding_VFregion(N_diffForces,IncrAmnt)
    
    
  end subroutine make_InptTIIF_for_all_CMPCS_frmZeroTIIF
    
  
  subroutine read_TIIF_value(applied_TIIF)
    implicit none
    real*8, intent(out) :: applied_TIIF
    integer             :: nodeOfTIIF
    
    nodeOfTIIF   = 2*(Hlf_Ncell+1)-1        ; write(*,*) nodeOfTIIF,"nodeOfTIIF"
    applied_TIIF = abs(CgYNode(nodeOfTIIF)) ; write(*,*) applied_TIIF,"applied_TIIF" 
    
  end subroutine read_TIIF_value
  
  
  subroutine remove_vrtclForce_and_toChkIf_the_PCScanreturn(CMv)
    implicit none
    integer, intent(in)  :: CMv
    integer              :: nodeL,nodeR,cell_val
    real*8               :: d1v ! distnceBfrForceRmv
    real*8               :: d2v ! distnceAftForceRmv
    real*8               :: d3v ! distnce after adjsting tilt
    real*8               :: d4v ! aft manplt Ap-Lt
    real*8               :: d5v ! distnce after adjsting tilt
    real*8               :: distnceToBeRecovrd,RemDistnceToBeRecovrd,fctrForSprMan
    real*8               :: hmL0Ap,hmL0Lt
    
    integer              :: BsORLt
    real*8               :: ExpctdCurvtr(1:Hlf_Ncell),CurrCurvtr(1:Hlf_Ncell)
    
    integer, allocatable :: ApclSp1(:),BsalSp1(:),LtrlSp1(:)
    integer, allocatable :: ApclSp2(:),BsalSp2(:),LtrlSp2(:)
    integer, allocatable :: ApclSp3(:),BsalSp3(:),LtrlSp3(:)
    integer, allocatable :: ApclSp4(:),BsalSp4(:),LtrlSp4(:)
    integer, allocatable :: ApclSp(:), BsalSp(:), LtrlSp(:)
    
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    d1v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce bfr pull/push"
    
    call remveForceFrmApclSrfcOfIC_frmCurrntVal() ; write(*,*) Frame_NI-1,"loc1"
    d2v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft rmv"
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"loc2"
    d3v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft bndry adjst"
    
    distnceToBeRecovrd = abs(d1v-d3v) ; write(*,*) distnceToBeRecovrd,"dR" 
    
    fctrForSprMan = 0.9000d0
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:NAEC_Bsal+1),LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:NAEC_Bsal+1),LtrlSp2(1:(NAEC_Ltrl+1)))
    allocate(ApclSp3(1:(NAEC_Apcl+1)),BsalSp3(1:NAEC_Bsal+1),LtrlSp3(1:(NAEC_Ltrl+1)))
    allocate(ApclSp4(1:(NAEC_Apcl+1)),BsalSp4(1:NAEC_Bsal+1),LtrlSp4(1:(NAEC_Ltrl+1)))
    allocate(ApclSp(1:(NAEC_Apcl+1)), BsalSp(1:NAEC_Bsal+1), LtrlSp(1:(NAEC_Ltrl+1)))
    
    ApclSp1=-1 ; BsalSp1=-1 ; LtrlSp1=-1 ; ApclSp2=-1 ; BsalSp2=-1 ; LtrlSp2=-1 
    ApclSp3=-1 ; BsalSp3=-1 ; LtrlSp3=-1 ; ApclSp4=-1 ; BsalSp4=-1 ; LtrlSp4=-1 
    ApclSp=-1  ; BsalSp=-1  ; LtrlSp=-1
    
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp1,BsalSp1,LtrlSp1)
    cell_val=Hlf_Ncell-1  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    cell_val=Hlf_Ncell-2  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    cell_val=Hlf_Ncell-3  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp4,BsalSp4,LtrlSp4)
    
    d4v = 0.00d0
    
    do
       
       if (abs(d1v).gt.abs(d4v)) continue
       if (abs(d1v).le.abs(d4v)) exit
       
       do i = 1,CMv
          
          if (i==1) then
             ApclSp = ApclSp1
             BsalSp = BsalSp1
             LtrlSp = LtrlSp1
          elseif (i==2) then
             ApclSp = ApclSp2
             BsalSp = BsalSp2
             LtrlSp = LtrlSp2
          elseif (i==3) then
             ApclSp = ApclSp3
             BsalSp = BsalSp3
             LtrlSp = LtrlSp3
          elseif (i==4) then
             ApclSp = ApclSp4
             BsalSp = BsalSp4
             LtrlSp = LtrlSp4
          endif   
          
          hmL0Ap=1.0100d0 ; call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
          hmL0Lt=0.9900d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
          call Equilibrate_only_NI_model
          
       enddo
       
       d4v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft ApLt man"
       RemDistnceToBeRecovrd = abs(d1v-d4v) ; write(*,*) RemDistnceToBeRecovrd,"Rem dR"
       
       if (RemDistnceToBeRecovrd .lt. ((1.00d0-fctrForSprMan)*distnceToBeRecovrd)) then
          write(*,*) RemDistnceToBeRecovrd,distnceToBeRecovrd,Frame_NI-1,"exit condn"
          exit
       endif
       
    enddo
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"loc3"
    d5v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft bndry adjst"
    
    writeForCurvtr=0 ; BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine remove_vrtclForce_and_toChkIf_the_PCScanreturn
  
  
  
  subroutine pullAtBndry_and_toChkIf_the_SystmCanReturn(CMv,lpVal)
    implicit none
    integer, intent(in)  :: CMv,lpVal
    real*8               :: CgXChng
    integer              :: nodeL,nodeR,cell_val
    real*8               :: d1v                  ! distnceBfrForceRmv
    real*8               :: d2v                  ! distnceAftForceRmv
    real*8               :: d3v                  ! distnce after adjsting tilt
    real*8               :: d4v                  ! aft manplt Ap-Lt
    real*8               :: d5v                  ! distnce after adjsting tilt
    real*8               :: distnceToBeRecovrd
    real*8               :: RemDistnceToBeRecovrd,fctrForSprMan
    real*8               :: hmL0Ap,hmL0Lt
    
    integer              :: BsORLt
    real*8               :: ExpctdCurvtr(1:Hlf_Ncell),CurrCurvtr(1:Hlf_Ncell)
    
    integer, allocatable :: ApclSp1(:),BsalSp1(:),LtrlSp1(:)
    integer, allocatable :: ApclSp2(:),BsalSp2(:),LtrlSp2(:)
    integer, allocatable :: ApclSp3(:),BsalSp3(:),LtrlSp3(:)
    integer, allocatable :: ApclSp4(:),BsalSp4(:),LtrlSp4(:)
    integer, allocatable :: ApclSp(:), BsalSp(:), LtrlSp(:)
    
    integer :: bndryLT,bndryLB,bndryRT,bndryRB
    integer :: pullpushIndctr
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    call get_the_shifting_nodes(nodeL,nodeR)
    
    if (lpVal==1) then
       d1v            = abs(node_xy(nodeL,2))
       disTnceBfrPull = d1v
       write(*,*) node_xy(nodeL,2),nodeL,"distnce bfr rmv"
    elseif (lpVal.gt.1) then
       d1v = disTnceBfrPull
    endif
    
    pullpushIndctr = 1                   ! 1 for pull, -1 for push
    CgXChng        = 0.12500d0
    call pull_push_atBndry(pullpushIndctr,CgXChng) ; write(*,*) Frame_NI-1,"loc1"
    
    d2v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft pull/push"
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"loc2"
    d3v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft bndry adjst"
    
    distnceToBeRecovrd = abs(d1v-d3v) ; write(*,*) distnceToBeRecovrd,"dR" 
    
    fctrForSprMan = 0.9500d0 !0.9000d0
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:NAEC_Bsal+1),LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:NAEC_Bsal+1),LtrlSp2(1:(NAEC_Ltrl+1)))
    allocate(ApclSp3(1:(NAEC_Apcl+1)),BsalSp3(1:NAEC_Bsal+1),LtrlSp3(1:(NAEC_Ltrl+1)))
    allocate(ApclSp4(1:(NAEC_Apcl+1)),BsalSp4(1:NAEC_Bsal+1),LtrlSp4(1:(NAEC_Ltrl+1)))
    allocate(ApclSp(1:(NAEC_Apcl+1)), BsalSp(1:NAEC_Bsal+1), LtrlSp(1:(NAEC_Ltrl+1)))
    
    ApclSp1=-1 ; BsalSp1=-1 ; LtrlSp1=-1 ; ApclSp2=-1 ; BsalSp2=-1 ; LtrlSp2=-1 
    ApclSp3=-1 ; BsalSp3=-1 ; LtrlSp3=-1 ; ApclSp4=-1 ; BsalSp4=-1 ; LtrlSp4=-1 
    ApclSp =-1 ; BsalSp =-1 ; LtrlSp =-1
    
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp1,BsalSp1,LtrlSp1)
    cell_val=Hlf_Ncell-1  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    cell_val=Hlf_Ncell-2  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    cell_val=Hlf_Ncell-3  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp4,BsalSp4,LtrlSp4)
    
    d4v = 0.00d0
    
    do
       
       if (abs(d1v).gt.abs(d4v)) continue
       if (abs(d1v).le.abs(d4v)) exit
       
       do i = 1,CMv
          
          if (i==1) then
             ApclSp = ApclSp1
             BsalSp = BsalSp1
             LtrlSp = LtrlSp1
          elseif (i==2) then
             ApclSp = ApclSp2
             BsalSp = BsalSp2
             LtrlSp = LtrlSp2
          elseif (i==3) then
             ApclSp = ApclSp3
             BsalSp = BsalSp3
             LtrlSp = LtrlSp3
          elseif (i==4) then
             ApclSp = ApclSp4
             BsalSp = BsalSp4
             LtrlSp = LtrlSp4
          endif   
          
          hmL0Ap=1.0100d0 ; call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
          hmL0Lt=0.9900d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
          call Equilibrate_only_NI_model
          
       enddo
       
       d4v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft ApLt man"
       RemDistnceToBeRecovrd = abs(d1v-d4v) ; write(*,*) RemDistnceToBeRecovrd,"Rem dR"
       
       if (RemDistnceToBeRecovrd .lt. ((1.00d0-fctrForSprMan)*distnceToBeRecovrd)) then
          write(*,*) RemDistnceToBeRecovrd,distnceToBeRecovrd,Frame_NI-1,"exit condn"
          exit
       endif
       
    enddo
    
    write(*,*) Frame_NI-1,"Frame NI bfr buckle manage"
    call manage_buckld_LtrlSpr_based_onTnsn() ; write(*,*) Frame_NI-1,"Frame NI aft buckle manage"
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"loc3"
    d5v = abs(node_xy(nodeL,2)) ; write(*,*) node_xy(nodeL,2),nodeL,"distnce aft bndry adjst"
    
    writeForCurvtr=0 ; BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine pullAtBndry_and_toChkIf_the_SystmCanReturn
  
  subroutine edit_after_a_specific_frame
    implicit none
    integer :: cellNm
    integer :: ApBsLt,cntLpMax
    real*8  :: hmL0
    integer :: BsORLt,CMval
    real*8  :: CurrCurvtr(1:Hlf_Ncell),ExpctdCurvtr(1:Hlf_Ncell)
    
    Frame_NI = 85 
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    cellNm = Hlf_Ncell-0 ; ApBsLt=3 ; hmL0=1.0100d0 ; cntLpMax = 3
    call Sprshorten_of_Cell_NIsys(cellNm,ApBsLt,hmL0,cntLpMax)    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    CMval = 1
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMval,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine edit_after_a_specific_frame
  
  subroutine take_currnt_dats_and_edit_forces_to_match(CMv)
    implicit none
    integer, intent(in) :: CMv
    real*8              :: CurrCurvtr(1:Hlf_Ncell),ExpctdCurvtr(1:Hlf_Ncell)
    real*8              :: ExpctdPPbfrVal,ExpctdPPaftVal
    integer             :: BsORLt,PCSorICS
    
    call allocate_and_initialize_PPValues
    
    Frame_NI = 49 !67(CM3PCS) !84 (CM2PCS) !84 (CM1PCS)
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    !call get_pulley_prtns_BfrAft(CMv)
    !ExpctdPPbfrVal = Expected_PPbfr(CMv) ; ExpctdPPaftVal = Expected_PPaft(CMv) 
    !call adjust_force_at_joint_ApclNode_ofIC_to_maintain_PP(CMv,ExpctdPPbfrVal,ExpctdPPaftVal)
    
    !writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    !call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    !call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    Frame_NI = 411 !98(CM3PCS) !365(CM2PCS) !142 (CM1PCS)
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    !PCSorICS=1 ; call make_curvtr_inDefinedRnge_gen_with_correct_curvtrOrientn(CMv,PCSorICS)
    !write(*,*) Frame_NI-1,"FrmNIII" ; call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine take_currnt_dats_and_edit_forces_to_match
  
  subroutine get_pulley_prtns_BfrAft(CMv)
    implicit none
    integer, intent(in) :: CMv
    real*8              :: ratioMeasuremnt_val,Unit_Fctr_PCS_val
    integer             :: TN_bfr,TN_aft
    
    call adjustmnt_Exprmntl_len_measuremnt(CMv,ratioMeasuremnt_val)
    write(*,*) CMv,ratioMeasuremnt_val,"CMv,ratioMeasuremnt_val"
    
    if (CMv==1) then
       Exprmntl_PPbfr(CMv) = 0.2700d0*ratioMeasuremnt_val
       Exprmntl_PPaft(CMv) = 0.2200d0*ratioMeasuremnt_val
       
    elseif (CMv==2) then
       Exprmntl_PPbfr(CMv) = 0.4900d0*ratioMeasuremnt_val
       Exprmntl_PPaft(CMv) = 0.2900d0*ratioMeasuremnt_val
       
    elseif (CMv==3) then
       Exprmntl_PPbfr(CMv) = 0.3400d0*ratioMeasuremnt_val
       Exprmntl_PPaft(CMv) = 0.4200d0*ratioMeasuremnt_val
       
    elseif (CMv==4) then
       Exprmntl_PPbfr(CMv) = 0.9400d0*ratioMeasuremnt_val
       Exprmntl_PPaft(CMv) = 0.2100d0*ratioMeasuremnt_val
    endif
    
    Unit_Fctr_PCS_val  = 3.5839837092863056d0  ! put the value here because the defining subroutine
                                               ! is not used during this test
    
    Expected_PPbfr(CMv) = Unit_Fctr_PCS_val*Exprmntl_PPbfr(CMv)
    Expected_PPaft(CMv) = Unit_Fctr_PCS_val*Exprmntl_PPaft(CMv)
    
    write(*,*) Exprmntl_PPbfr(CMv),Exprmntl_PPaft(CMv),CMv,"Exprmntl-CMv"
    write(*,*) Expected_PPbfr(CMv),Expected_PPaft(CMv),CMv,"Expected-CMv"
    
  end subroutine get_pulley_prtns_BfrAft
  
  
  subroutine Edit_forces_to_zero_and_match(CMv)
    implicit none
    integer, intent(in) :: CMv
    real*8              :: CurrCurvtr(1:Hlf_Ncell),ExpctdCurvtr(1:Hlf_Ncell)
    integer             :: BsORLt,PCSorICS
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    call remveForceFrmApclSrfcOfIC_frmCurrntVal()
    
    PCSorICS=1 ; call make_curvtr_inDefinedRnge_gen_with_correct_curvtrOrientn(CMv,PCSorICS)
    write(*,*) Frame_NI-1,"FrmNIiI"; call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine Edit_forces_to_zero_and_match
  
  subroutine apply_forces_at_ApclSrfc_to_create_flatBase
    implicit none
    real*8  :: base_level=-1.0d30
    integer :: Nforce
    integer :: nodeEndBsal1,nodeEndBsal2
    real*8  :: nodeYEndBsal,CgYNodeV
    integer :: nodeL,nodeR
    
    Frame_NI = 70 ! 37
    call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
    call Equilibrate_only_NI_model
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
    call get_the_base_level(base_level)      ; write(*,*) base_level,"base lev aft call 0"
    call get_the_shifting_nodes(nodeL,nodeR) ; write(*,*) nodeL,nodeR,"nodeL-nodeR"
    
    Nforce       = 1
    nodeEndBsal1 = 2*(Hlf_Ncell+1) ; nodeEndBsal2 = 2*(Hlf_Ncell+1)+nodeEndBsal1
    write(*,*) nodeEndBsal1,nodeEndBsal2,"nodeEndBsal1,nodeEndBsal2"
    
    do
       CgYNodeV       = CgYNode(nodeL)
       CgYNode(nodeL) = CgYNode(nodeL) + 0.0200d0 ; CgYNode(nodeR) = CgYNode(nodeL)  
       call Equilibrate_only_NI_model
       
       call get_the_base_level(base_level) ; write(*,*) base_level,"base lev aft call 1"
       
       nodeYEndBsal = 0.5000d0*(node_xy(nodeEndBsal1,2)+node_xy(nodeEndBsal2,2))
       write(*,*) nodeYEndBsal,node_xy(nodeEndBsal1,2),node_xy(nodeEndBsal2,2),"nodeYEndBsal"
       
       if (abs(nodeYEndBsal).le.abs(base_level)) then
          continue
       elseif (abs(nodeYEndBsal).gt.abs(base_level)) then
          CgYNode(nodeL) = CgYNodeV + 0.0100d0 ; CgYNode(nodeR) = CgYNode(nodeL)
          write(*,*) CgYNode(nodeL),CgYNode(nodeR),"at exit value" ; exit
       endif
       
    enddo
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    
  end subroutine apply_forces_at_ApclSrfc_to_create_flatBase
  
  subroutine get_the_base_level(base_level)
    implicit none
    real*8, intent(out) :: base_level
    integer             :: node1,node2,node3
    real*8              :: nodeY1,nodeY2,nodeY3
    
    node1  = 2*Hlf_Ncell      ; node2  = node1-2          ; node3  = node2-2
    nodeY1 = node_xy(node1,2) ; nodeY2 = node_xy(node2,2) ; nodeY3 = node_xy(node3,2) 
    
    write(*,*) node1,node2,node3,"nodes" ; write(*,*) nodeY1,nodeY2,nodeY3,"node Ys"
    
    base_level = (nodeY1+nodeY2+nodeY3)/3.0000d0 ; write(*,*) base_level,"bl"
    
  end subroutine get_the_base_level
  
  subroutine take_equivlnt_TIIF_systm_corrct_IC_nodePstn_and_curvtr(CMv)
    implicit none
    integer, intent(in) :: CMv
    integer             :: nodeToRead=-100,PCSorICS,BsORLt
    real*8              :: yval=-1.0d30
    real*8              :: ExpctdCurvtr(1:Hlf_Ncell),CurrCurvtr(1:Hlf_Ncell)
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
    call read_the_ynode_of_ppEquivalnceFile(CMv,nodeToRead,yval); write(*,*)nodeToRead,yval,"aa"
    call manplt_invgnting_rgn_toRetain_y_OfICApclNode(CMv,nodeToRead,yval)   
    
    call parameter_values_for_module_CS_match
    call allocate_and_initialize_lenValues
    call Control_State_match_expected_length_calctns_calls_togthr
    
    PCSorICS=1 ; call alter_pressure_to_match_the_exp_img(CMv,PCSorICS)
    call print_relevant_info_Actual_and_ExpctdLen(CMv)
    
    writeForCurvtr=0 ; BsORLt=2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine take_equivlnt_TIIF_systm_corrct_IC_nodePstn_and_curvtr
  
  subroutine read_the_ynode_of_ppEquivalnceFile(CMv,nodeToRead,yval)
    implicit none
    integer, intent(in)  :: CMv
    integer, intent(out) :: nodeToRead
    real*8,  intent(out) :: yval
    
    integer :: nodeNmbr,nodeTyp
    real*8  :: CgXval,CgYval,xNval,yNval
    
    if (CMv==1) open(unit=754,file='nodePrps_CM1PCS_TIIF_PulleyPntEquivalnce.dat')
    if (CMv==2) open(unit=754,file='nodePrps_CM2PCS_TIIF_PulleyPntEquivalnce.dat')
    if (CMv==3) open(unit=754,file='nodePrps_CM3PCS_TIIF_PulleyPntEquivalnce.dat')
    if (CMv==4) open(unit=754,file='nodePrps_CM4PCS_TIIF_PulleyPntEquivalnce.dat')
    
    nodeToRead = 2*(Hlf_Ncell)+1 ; write(*,*) nodeToRead,"nodeToRead"
    
    do
       read(754,*) nodeNmbr,nodeTyp,CgXval,CgYval,xNval,yNval
       if (nodeNmbr==nodeToRead) exit
    enddo
    
    yval = abs(yNval) ; write(*,*) yval,"for CM =",CMv
    
    close(754)
    
  end subroutine read_the_ynode_of_ppEquivalnceFile
  
  subroutine manplt_invgnting_rgn_toRetain_y_OfICApclNode(CMv,nodeToRead,yToAchieve)
    implicit none
    integer, intent(in) :: CMv
    integer, intent(in) :: nodeToRead
    real*8,  intent(in) :: yToAchieve
    
    integer :: i,j,cell_val,lpCnt,chng_indctn=0
    real*8  :: LL_yToAchieve,HL_yToAchieve ! LL=LowLimit, HL=HighLimit
    real*8  :: y_curr,hmL0Ap,hmL0Lt
    
    integer, allocatable :: ApclSp1(:),BsalSp1(:),LtrlSp1(:)
    integer, allocatable :: ApclSp2(:),BsalSp2(:),LtrlSp2(:)
    integer, allocatable :: ApclSp3(:),BsalSp3(:),LtrlSp3(:)
    integer, allocatable :: ApclSp4(:),BsalSp4(:),LtrlSp4(:)
    integer, allocatable :: ApclSp(:), BsalSp(:), LtrlSp(:)
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:NAEC_Bsal+1),LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:NAEC_Bsal+1),LtrlSp2(1:(NAEC_Ltrl+1)))
    allocate(ApclSp3(1:(NAEC_Apcl+1)),BsalSp3(1:NAEC_Bsal+1),LtrlSp3(1:(NAEC_Ltrl+1)))
    allocate(ApclSp4(1:(NAEC_Apcl+1)),BsalSp4(1:NAEC_Bsal+1),LtrlSp4(1:(NAEC_Ltrl+1)))
    allocate(ApclSp(1:(NAEC_Apcl+1)), BsalSp(1:NAEC_Bsal+1), LtrlSp(1:(NAEC_Ltrl+1)))
    
    ApclSp1=-1 ; BsalSp1=-1 ; LtrlSp1=-1 ; ApclSp2=-1 ; BsalSp2=-1 ; LtrlSp2=-1 
    ApclSp3=-1 ; BsalSp3=-1 ; LtrlSp3=-1 ; ApclSp4=-1 ; BsalSp4=-1 ; LtrlSp4=-1 
    ApclSp =-1 ; BsalSp =-1 ; LtrlSp =-1
    
    cell_val=Hlf_Ncell-0  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp1,BsalSp1,LtrlSp1)
    cell_val=Hlf_Ncell-1  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp2,BsalSp2,LtrlSp2)
    cell_val=Hlf_Ncell-2  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp3,BsalSp3,LtrlSp3)
    cell_val=Hlf_Ncell-3  ; call get_ApclBsalLtrl_OfNIsystm(cell_val,ApclSp4,BsalSp4,LtrlSp4)
    
    LL_yToAchieve = yToAchieve-0.0500d0 
    HL_yToAchieve = yToAchieve+0.0500d0
    lpCnt         = 1
    
    do
       
       y_curr = abs(node_xy(nodeToRead,2)) ; write(*,*) "y_curr =",y_curr,"for lpCnt =",lpCnt
       
       if ((y_curr.gt.LL_yToachieve).and.(y_curr.lt.HL_yToachieve)) exit
       if ((y_curr.lt.LL_yToachieve)) chng_indctn = +1 
       if ((y_curr.gt.HL_yToachieve)) chng_indctn = -1
       
       do i = 1,CMv
          if (i==1) then
             ApclSp = ApclSp1
             BsalSp = BsalSp1
             LtrlSp = LtrlSp1
          elseif (i==2) then
             ApclSp = ApclSp2
             BsalSp = BsalSp2
             LtrlSp = LtrlSp2
          elseif (i==3) then
             ApclSp = ApclSp3
             BsalSp = BsalSp3
             LtrlSp = LtrlSp3
          elseif (i==4) then
             ApclSp = ApclSp4
             BsalSp = BsalSp4
             LtrlSp = LtrlSp4
          endif
          
          if (chng_indctn==+1) then
             hmL0Ap=1.0100d0 ; call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
             hmL0Lt=0.9900d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
          elseif (chng_indctn==-1) then
             hmL0Ap=0.9900d0 ; call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
             hmL0Lt=1.0100d0 ; call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
          endif
          
       enddo
       
       call Equilibrate_only_NI_model ; lpCnt = lpCnt+1
       call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"ppE"
       
    enddo
    
    write(*,*) "Aft exiting loop in sb:manplt_invgnting_rgn_toRetain_y_OfICApclNode"
    
  end subroutine manplt_invgnting_rgn_toRetain_y_OfICApclNode
  
  subroutine TIIF_adjstmnt_ToASpcfcVal(CMv,TIIF_val)
    implicit none
    integer,intent(in) :: CMv 
    real*8, intent(in) :: TIIF_val
    real*8             :: curr_TIIF,diff_TIIF,tol_TIIF
    integer            :: LnodeOfTIIF,RnodeOfTIIF
    integer            :: BsORLt
    real*8             :: ExpctdCurvtr(1:Hlf_Ncell),CurrCurvtr(1:Hlf_Ncell)
    
    
    LnodeOfTIIF = 2*(Hlf_Ncell+1)-1 ; RnodeOfTIIF = LnodeOfTIIF+2*(Hlf_Ncell+1) 
    write(*,*) LnodeOfTIIF,RnodeOfTIIF,"L/RnodeofTIIF"
    
    tol_TIIF   = 0.1700d0
    curr_TIIF  = CgYNode(LnodeOfTIIF)
    diff_TIIF  = abs(TIIF_val-curr_TIIF)
    
    write(*,*) curr_TIIF,TIIF_val,diff_TIIF,"TIIF vals"
    
    if (diff_TIIF.gt.tol_TIIF) then
       write(*,*) "diff_TIIF is higher, decide on what to do"
       write(*,*) diff_TIIF,tol_TIIF,"diff-tol"
    elseif (diff_TIIF.le.tol_TIIF) then
       
       CgYNode(LnodeOfTIIF) = TIIF_val
       CgYNode(RnodeOfTIIF) = CgYNode(LnodeOfTIIF)
       call Equilibrate_only_NI_model
       
    endif
    
    call pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa() ; write(*,*) Frame_NI-1,"in sb:TIIF_adjstmnt_ToASpcfcVal"
    
    writeForCurvtr=0 ; BsORLt = 2 ; call get_Curvtr_Height_general(BsORLt,CurrCurvtr)
    call read_expctd_curvtr_height(CMv,ExpctdCurvtr)
    call compare_curvtrs_writing(CurrCurvtr,ExpctdCurvtr)
    
  end subroutine TIIF_adjstmnt_ToASpcfcVal
  
  
end module Control_State_Length_Match
