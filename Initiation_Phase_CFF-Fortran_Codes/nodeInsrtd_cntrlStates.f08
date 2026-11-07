module nodeInsrtd_cntrlStates
  use switch_and_unswitch_models
  
  implicit none
  real*8, allocatable :: A0_prgrsnStrct1(:,:) ,l0_prgrsnStrct1(:,:)
  real*8, allocatable :: ka_prgrsnStrct1(:,:) ,ks_prgrsnStrct1(:,:)
  real*8, allocatable :: CgX_prgrsnStrct1(:,:),CgY_prgrsnStrct1(:,:)
  
  real*8, allocatable :: A0_prgrsnStrct2(:,:), l0_prgrsnStrct2(:,:)
  real*8, allocatable :: ka_prgrsnStrct2(:,:), ks_prgrsnStrct2(:,:)
  real*8, allocatable :: CgX_prgrsnStrct2(:,:),CgY_prgrsnStrct2(:,:)
  
  real*8, allocatable :: A0_prgrsnStrct3(:,:), l0_prgrsnStrct3(:,:)
  real*8, allocatable :: ka_prgrsnStrct3(:,:), ks_prgrsnStrct3(:,:)
  real*8, allocatable :: CgX_prgrsnStrct3(:,:),CgY_prgrsnStrct3(:,:)
  
  real*8, allocatable :: A0_prgrsnStrct4(:,:) ,l0_prgrsnStrct4(:,:)
  real*8, allocatable :: ka_prgrsnStrct4(:,:) ,ks_prgrsnStrct4(:,:)
  real*8, allocatable :: CgX_prgrsnStrct4(:,:),CgY_prgrsnStrct4(:,:)
  
  real*8, allocatable :: A0_prgrsnNIstrct1(:,:) ,l0_prgrsnNIstrct1(:,:)
  real*8, allocatable :: ka_prgrsnNIstrct1(:,:) ,ks_prgrsnNIstrct1(:,:)
  real*8, allocatable :: CgX_prgrsnNIstrct1(:,:),CgY_prgrsnNIstrct1(:,:)
  
  real*8, allocatable :: A0_prgrsnNIstrct2(:,:) ,l0_prgrsnNIstrct2(:,:)
  real*8, allocatable :: ka_prgrsnNIstrct2(:,:) ,ks_prgrsnNIstrct2(:,:)
  real*8, allocatable :: CgX_prgrsnNIstrct2(:,:),CgY_prgrsnNIstrct2(:,:)
  
  real*8, allocatable :: A0_prgrsnNIstrct3(:,:) ,l0_prgrsnNIstrct3(:,:)
  real*8, allocatable :: ka_prgrsnNIstrct3(:,:) ,ks_prgrsnNIstrct3(:,:)
  real*8, allocatable :: CgX_prgrsnNIstrct3(:,:),CgY_prgrsnNIstrct3(:,:)
  
  real*8, allocatable :: A0_prgrsnNIstrct4(:,:) ,l0_prgrsnNIstrct4(:,:)
  real*8, allocatable :: ka_prgrsnNIstrct4(:,:) ,ks_prgrsnNIstrct4(:,:)
  real*8, allocatable :: CgX_prgrsnNIstrct4(:,:),CgY_prgrsnNIstrct4(:,:)
  
  real*8, allocatable :: ks_strctIF_NI(:,:),l0_strctIF_NI(:,:) !IF =InterFace matching ...
  real*8, allocatable :: ka_strctIF_NI(:,:),A0_strctIF_NI(:,:)
  
  integer              :: N_rectCellPair,N_rectCells,N_chngCells
  integer, allocatable :: RectCells(:),ChngCells(:)
  
  real*8,  allocatable :: A0I_NI(:),kaI_NI(:),l0I_NI(:),ksI_NI(:)
  real*8,  allocatable :: A0F_NI(:),kaF_NI(:),l0F_NI(:),ksF_NI(:)
  
  integer :: invgnting_cells(1:2),invgnated_cells(1:2),nxt_invgnated_cells(1:2)
  integer :: frstCells(1:2),scndCells(1:2),thrdCells(1:2),frthCells(1:2)
  integer :: thrdlastCells(1:2),scndlastCells(1:2),frstlastCells(1:2)
  
  integer :: CyclVal !will be used for read_strctProps_forWTpropsOfMultiplCycl
  integer :: prtrbtnTestcallLoctn
  
  real*8, allocatable :: LenAp_Exp(:),LenBs_Exp(:),LenLt_Exp(:) !Length(Apical/Basal/Lateral)_Experimental
  real*8, allocatable :: RL_Ap(:),RL_Bs(:),RL_Lt(:) !RatioLength_(Apical/Basal/Lateral)
  
contains
  
  subroutine allocate_ProgrsnStgProp
    implicit none
    
    allocate(A0_prgrsnStrct1(1:N_progCycl,1:N_cell), l0_prgrsnStrct1(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnStrct1(1:N_progCycl,1:N_cell), ks_prgrsnStrct1(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnStrct1(1:N_progCycl,1:N_node),CgY_prgrsnStrct1(1:N_progCycl,1:N_node))
    
    allocate(A0_prgrsnStrct2(1:N_progCycl,1:N_cell), l0_prgrsnStrct2(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnStrct2(1:N_progCycl,1:N_cell), ks_prgrsnStrct2(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnStrct2(1:N_progCycl,1:N_node),CgY_prgrsnStrct2(1:N_progCycl,1:N_node))
    
    allocate(A0_prgrsnStrct3(1:N_progCycl,1:N_cell), l0_prgrsnStrct3(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnStrct3(1:N_progCycl,1:N_cell), ks_prgrsnStrct3(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnStrct3(1:N_progCycl,1:N_node),CgY_prgrsnStrct3(1:N_progCycl,1:N_node))
    
    allocate(A0_prgrsnStrct4(1:N_progCycl,1:N_cell), l0_prgrsnStrct4(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnStrct4(1:N_progCycl,1:N_cell), ks_prgrsnStrct4(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnStrct4(1:N_progCycl,1:N_node),CgY_prgrsnStrct4(1:N_progCycl,1:N_node))
    
  end subroutine allocate_ProgrsnStgProp
  
  
  subroutine allocate_ProgrsnStgNIProp
    implicit none
    real*8  :: E
    
    call switch_to_NI_model
    
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    allocate(A0_prgrsnNIstrct1(1:N_progCycl,1:N_cell), l0_prgrsnNIstrct1(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnNIstrct1(1:N_progCycl,1:N_cell), ks_prgrsnNIstrct1(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnNIstrct1(1:N_progCycl,1:N_node),CgY_prgrsnNIstrct1(1:N_progCycl,1:N_node))
     
    allocate(A0_prgrsnNIstrct2(1:N_progCycl,1:N_cell), l0_prgrsnNIstrct2(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnNIstrct2(1:N_progCycl,1:N_cell), ks_prgrsnNIstrct2(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnNIstrct2(1:N_progCycl,1:N_node),CgY_prgrsnNIstrct2(1:N_progCycl,1:N_node))
    
    allocate(A0_prgrsnNIstrct3(1:N_progCycl,1:N_cell), l0_prgrsnNIstrct3(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnNIstrct3(1:N_progCycl,1:N_cell), ks_prgrsnNIstrct3(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnNIstrct3(1:N_progCycl,1:N_node),CgY_prgrsnNIstrct3(1:N_progCycl,1:N_node))
     
    allocate(A0_prgrsnNIstrct4(1:N_progCycl,1:N_cell), l0_prgrsnNIstrct4(1:N_progCycl,1:N_spr))
    allocate(ka_prgrsnNIstrct4(1:N_progCycl,1:N_cell), ks_prgrsnNIstrct4(1:N_progCycl,1:N_spr))
    allocate(CgX_prgrsnNIstrct4(1:N_progCycl,1:N_node),CgY_prgrsnNIstrct4(1:N_progCycl,1:N_node))
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine allocate_ProgrsnStgNIProp
  
  
  subroutine switchto_NI_model_save_and_switchbackto_TN(ProgCycl,strctV1,strctV2,strghtOrNI,itrnNo)
    implicit none
    integer, intent(in)    :: ProgCycl
    integer, intent(in)    :: strctV1,strctV2
    integer, intent(inout) :: strghtOrNI
    integer, intent(in)    :: itrnNo
    
    real*8  :: E
    logical :: lgcl_meet
    integer :: strctV3
    integer :: i
    
    real*8, allocatable :: A0V(:),l0V(:)
    real*8, allocatable :: kaV(:),ksV(:)
    
    open(unit=95,file='switchNI_and_save.dat',position='append')
    write(95,*) ProgCycl,strctV1,strctV2,strghtOrNI,itrnNo,nodeInsrtnStrts
    write(95,*)"ProgCycl,strctV1,strctV2,strghtOrNI,itrnNo,nodeInsrtnStrts"
    close(95)
    
    if (strghtOrNI.ne.2) then
       write(*,*) "Error in strghtOrNI variable value which is =",strghtOrNI
       stop
    endif
    
    if (nodeInsrtnStrts==1) then
       
       call switch_to_NI_model
       
       E = Energy(node_xy,l0,A0)
       write(*,*) E,"E"
       call get_gradient(node_xy,l0,A0,gradient)
       
       if (NI_alrdySaved==0) then
          call Equilibrate_system
          call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          
          lgcl_meet = .False.
          call decision_of_meeting_NI(lgcl_meet)
          
       elseif (NI_alrdySaved==1) then
          lgcl_meet = .False.
          !no need for decision_of_meeting_NI as it has been done already
       endif
       
       allocate(A0V(1:N_cell),l0V(1:N_spr))
       allocate(kaV(1:N_cell),ksV(1:N_spr))
       
       if (itrnNo==0) then
          
          A0V(1:N_cell) = A0(1:N_cell)     ; l0V(1:N_spr) = l0(1:N_spr)
          kaV(1:N_cell) = k_area(1:N_cell) ; ksV(1:N_spr) = k_spr(1:N_spr)
          
          if (strctV1==1) then
             call save_ProgrsnStgProp(ProgCycl,strctV1,strghtOrNI,A0V,l0V,kaV,ksV)
             call sleep(1)
             
          elseif (strctV1==3) then
             
             open(unit=78,file='notSavingStrct3.dat')
             write(78,*) "strct3 is not being saved here"
             close(78)
             
          endif
          
       elseif (itrnNo.ne.0) then
          
          if (lgcl_meet.eqv..True.) then
             
             A0V(1:N_cell) = A0(1:N_cell)     ; l0V(1:N_spr) = l0(1:N_spr)
             kaV(1:N_cell) = k_area(1:N_cell) ; ksV(1:N_spr) = k_spr(1:N_spr)
             
             call save_ProgrsnStgProp(ProgCycl,strctV2,strghtOrNI,A0V,l0V,kaV,ksV)
             
             if (strctV2==2) then !adjustment for strct 3
                strctV3 = 3
                call save_ProgrsnStgProp(ProgCycl,strctV3,strghtOrNI,A0V,l0V,kaV,ksV)
             endif
             
             Frame_NI = Frame_NI+1
             NI_alrdySaved = 1
             
             open(unit=96,file='switchNI_and_saveInsd.dat',position='append')
             write(96,*) ProgCycl,strctV1,strctV2,strghtOrNI,"ProgCycl,strctV1-2,strghtOrNI"
             write(96,*) NI_alrdySaved,Frame_NI,"NI_alrdySaved,Frame_NI"
             close(96)
             
          elseif (lgcl_meet.eqv..False.) then
             
             if (itrnNo==Nstp_PT) then
                
                A0V(1:N_cell) = A0(1:N_cell)     ; l0V(1:N_spr) = l0(1:N_spr)
                kaV(1:N_cell) = k_area(1:N_cell) ; ksV(1:N_spr) = k_spr(1:N_spr)
                call save_ProgrsnStgProp(ProgCycl,strctV2,strghtOrNI,A0V,l0V,kaV,ksV)
                call sleep(1)
                
             endif
             
          endif
          
       endif
       
       if (NI_alrdySaved==0) then
          Frame_NI = Frame_NI+1
       elseif (NI_alrdySaved==1) then
          continue
       endif
       
       call Find_Analytical_and_Numerical_Mismatch
       call sleep(1)
       call deallocate_repetitive_arrays
       
       deallocate(A0V,l0V,kaV,ksV)
       call switchback_to_TN_model
       strghtOrNI = 1
       
    endif
    
    
  end subroutine switchto_NI_model_save_and_switchbackto_TN
  
  
  subroutine adjustPress_Varying_A0orKa_of_NI(seqNo)
    implicit none
    integer, intent(in) :: seqNo
    
    integer :: strtSeq
    integer :: i
    integer :: LastCellToCmpr
    integer :: L1Cell,L2Cell,R1Cell,R2Cell
    real*8  :: PressureVal(1:N_cell)
    real*8  :: E,hmuch
    integer :: no_of_cells
    integer :: A0orKA
    integer :: Lstrt,Lfnsh,Rstrt,Rfnsh
    integer :: storedOrNot
    
    integer, allocatable :: cell_nos(:)
    real*8 , allocatable :: sprlenBfr(:),sprlenVar(:)
    
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
    
    
    strtSeq = 9 !9 comes from iMain = 1,10 [as I am only concerned for last two cycle]
    A0orKA  = 1
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    
    allocate(sprlenBfr(1:N_spr),sprlenVar(1:N_spr))
    
    PressureVal(1:N_cell) = Pressure(node_xy,A0)
    
    open(unit=34,file='Frame_NI.dat',position='append')
    open(unit=37,file='AdjustPressChk.dat',position='append')
    
    write(34,*) Frame_NI,"start"
    Frame_NI = Frame_NI+1
    
    LastCellToCmpr = 6 - (seqNo-strtSeq) !6
    
    no_of_cells = 4 ! 4 due to 2 pair of cells, Left and Right
    allocate(cell_nos(1:no_of_cells))
    
    do i = 1,LastCellToCmpr
       
       L1Cell=LastCellToCmpr-(i-1)   ; R1Cell=(N_cell/2)+LastCellToCmpr-(i-1)
       L2Cell=LastCellToCmpr-(i-1)+1 ; R2Cell=(N_cell/2)+LastCellToCmpr-(i-1)+1
       
       write(37,*) L1Cell,L2Cell,R1Cell,R2Cell,"LR Cells"
       
       cell_nos(1) = L1Cell ; cell_nos(2) = R1Cell
       cell_nos(3) = L2Cell ; cell_nos(4) = R2Cell
       
       write(37,*) cell_nos(1),cell_nos(3),cell_nos(2),cell_nos(4),"cell_nos"
       
       
       storedOrNot = 0
       
       do 
          
          if (PressureVal(L1Cell) .le. PressureVal(L2Cell)) then
             
             if (storedOrNot==0) then
                sprlenBfr(1:N_spr) = l(1:N_spr)
                storedOrNot = 1
             endif
             
             hmuch = 1.00d0 + (0.03d0) ! 1 prcnt increase
             
             call increase_A0orkarea_of_Prtclr_cells(A0orKA,cell_nos,no_of_cells,hmuch)
             call Equilibrate_system
             call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
             write(34,*) Frame_NI,"in Do loop"
             Frame_NI = Frame_NI+1
             
             PressureVal(1:N_cell) = Pressure(node_xy,A0)
             
             if (PressureVal(L1Cell) .gt. PressureVal(L2Cell)) then
                write(34,*) Frame_NI,"Ending"
                
                sprlenVar(1:N_spr) = abs(l(1:N_spr)-sprlenBfr(1:N_spr))
                !close(34)
                call reduce_l0_of_Prtclr_cells_Spr_to_retainShape(cell_nos,no_of_cells,sprlenVar,sprlenBfr)
                !stop
                exit
             endif
             
          elseif (PressureVal(L1Cell) .gt. PressureVal(L2Cell)) then
             exit
          endif
          
       enddo
       
    enddo
    
    Lstrt = 1            ; Lfnsh = (LastCellToCmpr+1)
    Rstrt = (N_cell/2)+1 ; Rfnsh = (N_cell/2)+ (LastCellToCmpr+1)
    
    write(37,*) PressureVal(Lstrt:Lfnsh),"Lcell Pressures"
    write(37,*) PressureVal(Rstrt:Rfnsh),"Rcell Pressures"
    write(37,*) (Frame_NI-1),"Frame NI no"
    
    close(34)
    close(37)
    
    call Find_Analytical_and_Numerical_Mismatch
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine adjustPress_Varying_A0orKa_of_NI
  
  subroutine increase_A0orkarea_of_Prtclr_cells(A0orKA,cell_nos,no_of_cells,hmuch)
    implicit none
    integer, intent(in) :: A0orKA
    integer, intent(in) :: no_of_cells
    integer, intent(in) :: cell_nos(1:no_of_cells)
    real*8 , intent(in) :: hmuch
    
    integer :: i
    integer :: cell_no
    real*8  :: prev_A0,prev_ka
    integer :: hlf_NoOfCells
    
    hlf_NoOfCells = no_of_cells/2
    
    do i = 1,hlf_NoOfCells
       
       cell_no = cell_nos(i)
       !write(*,*) cell_no,"cell_no"
       
       if (A0orKA==1) then
          prev_A0 = A0(cell_no)
          A0(cell_no) = (hmuch)*(prev_A0)
       elseif (A0orKA==2) then
          prev_ka = k_area(cell_no)
          k_area(cell_no) = (hmuch)*(prev_ka) 
       endif
       
    enddo
    
    write(*,*) cell_nos(1:no_of_cells),"cell nos"
    
  end subroutine increase_A0orkarea_of_Prtclr_cells
  
  
  subroutine reduce_l0_of_Prtclr_cells_Spr_to_retainShape(cell_nos,no_of_cells,sprlenVarPrv,sprlenBfr)
    implicit none
    integer, intent(in) :: no_of_cells
    integer, intent(in) :: cell_nos(1:no_of_cells)
    real*8 , intent(in) :: sprlenVarPrv(1:N_spr),sprlenBfr(1:N_spr)
    
    integer :: i,j 
    integer :: hlf_NoOfCells,cell_no,spr_no
    integer :: nspr_in_cell,nspr_toAdjust
    real*8  :: currl0(1:N_spr)
    real*8  :: sprlenVarNow(1:N_spr)
    real*8  :: hmuch
    integer :: cnt_Adjstd
    
    integer :: spr_typ
    real*8  :: tol_prcnt
    logical :: lgcl_Adjstd
    
    integer, allocatable :: sprToBeAdjstd(:)
    real*8 , allocatable :: prcnt_chng(:)
    
    open(unit=27,file='reduce_l0_of_Prtclr_cells_Spr.dat',position='append')
    write(27,*) "start of the call"
    
    hlf_NoOfCells = no_of_cells/2
    nspr_toAdjust = 0
    
    do i = 1,hlf_NoOfCells
       
       cell_no = cell_nos(i)
       
       nspr_in_cell  = area_spr(cell_no,0)
       nspr_toAdjust = nspr_toAdjust + nspr_in_cell
       write(27,*) nspr_in_cell,nspr_toAdjust,"nspr_in_cell and toAdjust"
       
    enddo
    
    allocate(sprToBeAdjstd(1:nspr_toAdjust))
    cnt_Adjstd = 0
    
    do i = 1,hlf_NoOfCells
       
       cell_no = cell_nos(i)
       nspr_in_cell = area_spr(cell_no,0)
       
       do j = 1,nspr_in_cell
          
          cnt_Adjstd = cnt_Adjstd+1
          sprToBeAdjstd(cnt_Adjstd) = area_spr(cell_no,j)
          write(27,*) sprToBeAdjstd(cnt_Adjstd),typ_spr(sprToBeAdjstd(cnt_Adjstd)),cnt_Adjstd,"sprToBeAdjstd,cnt"
       enddo
       
    enddo
    
    currl0(1:N_spr) = l0(1:N_spr)
    hmuch           = 0.05d0 ! 5 percent change
    
    allocate(prcnt_chng(1:nspr_toAdjust))
    
    
    do 
       do i = 1,nspr_toAdjust
          spr_no = sprToBeAdjstd(i)
          l0(spr_no) = (1.00-hmuch)*l0(spr_no)
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       sprlenVarNow(1:N_spr) = abs(l(1:N_spr)-sprlenBfr(1:N_spr))
       
       
       do i = 1,nspr_toAdjust
          
          spr_no = sprToBeAdjstd(i)
          write(*,*) sprlenVarNow(spr_no),sprlenVarPrv(spr_no),spr_no,"sprlenVar now,prv,spr_no"
          !prcnt_chng(i) = ( abs((sprlenVarNow(spr_no) - sprlenVarPrv(spr_no)))/(sprlenVarPrv(spr_no)) ) * (100.00)
          
          prcnt_chng(i) = (abs(l(spr_no)-sprlenBfr(spr_no)) / sprlenBfr(spr_no)) * (100.00)
          write(*,*)  prcnt_chng(i),i,"prcnt"
          write(27,*) prcnt_chng(i),i,"prcnt"
          
       enddo
       
       tol_prcnt   = 5.00d0
       lgcl_Adjstd = .False.
       
       !%%%%% (CHEKING FOR TOLERATING CONDITIONS) %%%%%!
       
       do i = 1,nspr_toAdjust
          
          spr_no  = sprToBeAdjstd(i)
          spr_typ = typ_spr(spr_no)
          
          if (spr_typ==1 .or. spr_typ==2 .or. spr_typ==3 .or. spr_typ==4) then
             
             if (prcnt_chng(i) .le. tol_prcnt) then
                lgcl_Adjstd = .True.
             elseif (prcnt_chng(i) .gt. tol_prcnt) then
                lgcl_Adjstd = .False.
                exit
             endif
             
          endif
          
       enddo
       
       if (lgcl_Adjstd .eqv. .True.) then
          exit
       elseif (lgcl_Adjstd .eqv. .False.) then
          continue
       endif
       
       !if (Frame_NI==9) exit
       
    enddo
    
    close(27)
    
  end subroutine reduce_l0_of_Prtclr_cells_Spr_to_retainShape
  
  
  subroutine change_A0of_CellPairOrSinglCell(pairOrSingl,cellNm,hmuch)
    implicit none
    integer, intent(in) :: pairOrSingl
    integer, intent(in) :: cellNm
    real*8 , intent(in) :: hmuch
    
    integer :: i,cellNo
    real*8  :: prev_A0
    
    open(unit=84,file='chngA0_0fPairorSingl.dat',position='append')
    
    if (pairOrSingl==1) then
       cellNo     = cellNm
       prev_A0    = A0(cellNo)
       A0(cellNo) = (hmuch)*(prev_A0)
       
       write(84,*) cellNo,prev_A0,A0(cellNo),"singl"
       
    elseif (pairOrSingl==2) then
       
       do i = 1,2
          
          if (i==1) cellNo = cellNm
          if (i==2) cellNo = Hlf_Ncell + cellNm
          
          prev_A0    = A0(cellNo)
          A0(cellNo) = (hmuch)*prev_A0
          
          write(84,*) cellNo,prev_A0,A0(cellNo),i,"pair"
          
       enddo
       
    endif
    
    close(84)
    
  end subroutine change_A0of_CellPairOrSinglCell
  
  subroutine change_l0of_CellPairOrSinglCell(pairOrSingl,cellNm,hmuch)
    implicit none
    integer, intent(in) :: pairOrSingl
    integer, intent(in) :: cellNm
    real*8 , intent(in) :: hmuch
    
    integer :: i,j
    integer :: cellNo,sprNo
    real*8  :: prev_l0
    integer :: nsprs_inCell
    
    open(unit=85,file='chngl0_0fPairorSingl.dat',position='append')
    
    if (pairOrSingl==1) then
       cellNo     = cellNm
       
       nsprs_inCell = area_spr(cellNo,0)
       
       do j = 1,nsprs_inCell
          sprNo      = area_spr(cellNo,j)
          prev_l0    = l0(sprNo)
          l0(sprNo)  = (hmuch)*(prev_l0)
          
          write(85,*) sprNo,prev_l0,l0(sprNo),j,"singl"
       enddo
       
       write(85,*) " "
       
    elseif (pairOrSingl==2) then
       
       do i = 1,2
          
          if (i==1) cellNo = cellNm
          if (i==2) cellNo = Hlf_Ncell + cellNm
          
          nsprs_inCell = area_spr(cellNo,0)
          
          do j = 1,nsprs_inCell
             sprNo      = area_spr(cellNo,j)
             prev_l0    = l0(sprNo)
             l0(sprNo) = (hmuch)*(prev_l0)
             
             write(85,*) sprNo,prev_l0,l0(sprNo),i,j,"pair"
          enddo
          
       enddo
       
    endif
    
    write(85,*) " "
    
    close(85)
    
  end subroutine change_l0of_CellPairOrSinglCell
    
    
  subroutine adjustPress_VaryingConsec_A0orKa_of_NI(seqNo)
    implicit none
    integer, intent(in) :: seqNo
    
    integer :: strtSeq
    integer :: i,j
    integer :: LastCellToCmpr
    integer :: L1Cell,L2Cell,R1Cell,R2Cell
    real*8  :: PressureVal(1:N_cell)
    real*8  :: E,hmuch
    integer :: no_of_CC
    integer :: A0orKA
    integer :: Lstrt,Lfnsh,Rstrt,Rfnsh
    integer :: storedOrNot
    
    integer, allocatable :: consecCells(:)
    real*8 , allocatable :: sprlenBfr(:)  ,sprlenVar(:)
    real*8 , allocatable :: cellAreaBfr(:),cellAreaVar(:)
    
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
    
    
    strtSeq = 9 ! 9 is frm iMain=1,10[as I'm only concerned for last two cycls]
    A0orKA  = 1
    
    allocate(sprlenBfr(1:N_spr),sprlenVar(1:N_spr))
    sprlenBfr(1:N_spr) = 1.0d-30 ; sprlenVar(1:N_spr) = 1.0d-30
    
    allocate(cellAreaBfr(1:N_cell),cellAreaVar(1:N_cell))
    
    PressureVal(1:N_cell) = Pressure(node_xy,A0)
    
    open(unit=34,file='Frame_NI.dat',position='append')
    open(unit=37,file='AdjustPressChk.dat',position='append')
    
    write(34,*) Frame_NI,"start"
    Frame_NI = Frame_NI+1
    
    LastCellToCmpr = 6 - (seqNo-strtSeq) !6
    
    close(34)
    
    do i = 1,LastCellToCmpr
       
       open(unit=54,file='pressChkInLoop.dat',position='append')
       write(54,*)PressureVal(1:(LastCellToCmpr+1)),LastCellToCmpr,"PressureVal"
       write(54,*)seqNo,Frame_NI,"seqNo,Frame_NI"
       
       L1Cell = i ; L2Cell = (i+1)
       
       if (PressureVal(L1Cell) .le. PressureVal(L2Cell)) then
          no_of_CC = i
          allocate(consecCells(1:no_of_CC))
          
          do j = 1,no_of_CC
             consecCells(j) = j
          enddo
          
          call pressureIncreament_uptothe_problemeticCell(A0orKA,consecCells,no_of_CC,&
               LastCellToCmpr,sprlenBfr,sprlenVar,cellAreaBfr,cellAreaVar,PressureVal)
          
          deallocate(consecCells)
          
       elseif (PressureVal(L1Cell) .gt. PressureVal(L2Cell)) then
          continue
       endif
       
       if (i==LastCellToCmpr) then
          write(54,*) "i at End"
          
          no_of_CC = i
          allocate(consecCells(1:no_of_CC))
          
          do j = 1,no_of_CC
             consecCells(j) = j
          enddo
          
          call make_higherdP_inLastCell(A0orKA,consecCells,no_of_CC,LastCellToCmpr)
          
          deallocate(consecCells)
          
       endif
       
       close(54)
       
    enddo
    
    !stop
    
    !if (seqNo==10) stop

    call maintain_precise_dP_in_bfrInvaginatingCells(seqNo,LastCellToCmpr)
    
    Lstrt = 1            ; Lfnsh = (LastCellToCmpr+1)
    Rstrt = (N_cell/2)+1 ; Rfnsh = (N_cell/2)+ (LastCellToCmpr+1)
    
    write(37,*) PressureVal(Lstrt:Lfnsh),"Lcell Pressures"
    write(37,*) PressureVal(Rstrt:Rfnsh),"Rcell Pressures"
    write(37,*) (Frame_NI-1),"Frame NI no"
    
    close(37)
    
    !if (hlfCyclNo==2) then
    !   call reduce_buckling_of_LtrlMmbrne(seqNo,LastCellToCmpr)
    !endif
    
    !call Find_Analytical_and_Numerical_Mismatch
    !call deallocate_repetitive_arrays
    !call switchback_to_TN_model
    
    !if (seqNo==10) stop
    
  end subroutine adjustPress_VaryingConsec_A0orKa_of_NI
  
  
  subroutine pressureIncreament_uptothe_problemeticCell(A0orKA,consecCells,no_of_CC,LastCellToCmpr,&
       sprlenBfr,sprlenVar,cellAreaBfr,cellAreaVar,PressureVal)
    implicit none
    integer,intent(in)    :: A0orKA,no_of_CC
    integer,intent(in)    :: consecCells(1:no_of_CC)
    integer,intent(in)    :: LastCellToCmpr
    real*8, intent(inout) :: sprlenBfr(1:N_spr),sprlenVar(1:N_spr)
    real*8, intent(inout) :: cellAreaBfr(1:N_cell),cellAreaVar(1:N_cell)
    real*8, intent(inout) :: PressureVal(1:N_cell)
    
    real*8  :: hmuch
    integer :: L1Cell,L2Cell
    
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
    
    open(unit=34,file='Frame_NI.dat',position='append')
    
    sprlenBfr(1:N_spr)    = l(1:N_spr)
    cellAreaBfr(1:N_cell) = A(1:N_cell)
    
    L1Cell = no_of_CC ; L2Cell = (no_of_CC+1) 
    write(*,*) L1Cell,LastCellToCmpr,"L1Cell,LastCellToCmpr"
    
    do 
       
       if (PressureVal(L1Cell) .le. PressureVal(L2Cell)) then
          
          hmuch = 1.00d0 + (0.03d0) ! 3 prcnt increase
          
          call increase_A0_ofConsecutiveCells(A0orKA,consecCells,no_of_CC,hmuch)
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          write(34,*) Frame_NI,"in Do loop"
          Frame_NI = Frame_NI+1
          
          PressureVal(1:N_cell) = Pressure(node_xy,A0)
          
          if (PressureVal(L1Cell) .gt. PressureVal(L2Cell)) then
             write(34,*) Frame_NI,"Ending"
             
             !if (L1Cell==LastCellToCmpr) then
              !  call make_higherdP_inLastCell(A0orKA,consecCells,no_of_CC,LastCellToCmpr)
             !endif
             
             sprlenVar(1:N_spr)    = abs(l(1:N_spr)-sprlenBfr(1:N_spr))
             cellAreaVar(1:N_cell) = abs(A(1:N_cell)-cellAreaBfr(1:N_cell))
             
             call reduce_l0_of_Consec_cells_Spr_to_retainShape(consecCells,no_of_CC,sprlenVar,sprlenBfr,cellAreaVar,cellAreaBfr)
             
             PressureVal(1:N_cell) = Pressure(node_xy,A0)
             
             if (PressureVal(L1Cell).le.PressureVal(L2Cell)) then
                write(*,*) "Pressure altered after reducing l0, should not happen,fl:nodeInsrtd; sb:pressureIncr"
                write(*,*) PressureVal(L1Cell),PressureVal(L2Cell),L1Cell,L2Cell," PressureVals of L1,L2 Cells"
                call sleep(1)
                continue !stop
             endif
             
             exit
          endif
          
       elseif (PressureVal(L1Cell) .gt. PressureVal(L2Cell)) then
          exit
       endif
       
    enddo
    
    close(34)
    
  end subroutine pressureIncreament_uptothe_problemeticcell
  
  subroutine increase_A0_ofConsecutiveCells(A0orKA,consecCells,no_of_CC,hmuch)
    implicit none
    integer, intent(in) :: A0orKA
    integer, intent(in) :: no_of_CC
    integer, intent(in) :: consecCells(1:no_of_CC) !CC=Consec Cells
    real*8 , intent(in) :: hmuch
    
    integer :: i,LftCell,RghtCell
    real*8  :: prev_A0L,prev_A0R
    real*8  :: prev_KAL,prev_KAR
    
    do i = 1,no_of_CC
       LftCell = consecCells(i) ; RghtCell = (N_cell/2)+consecCells(i)
       
       if (A0orKA==1) then
          
          prev_A0L = A0(LftCell) ; prev_A0R = A0(RghtCell) 
          
          A0(LftCell)  = hmuch*prev_A0L
          A0(RghtCell) = hmuch*prev_A0R
          
       elseif (A0orKA==2) then
          
          prev_KAL = k_area(LftCell) ; prev_KAR = k_area(RghtCell) 
          
          k_area(LftCell)  = hmuch*prev_KAL
          k_area(RghtCell) = hmuch*prev_KAR
          
       endif
       
    enddo
    
  end subroutine increase_A0_ofConsecutiveCells
  
  subroutine make_higherdP_inLastCell(A0orKA,consecCells,no_of_CC,LastCellToCmpr)
    implicit none
    integer, intent(in) :: A0orKA, no_of_CC
    integer, intent(in) :: consecCells(1:no_of_CC)
    integer, intent(in) :: LastCellToCmpr
    real*8, allocatable :: dP_calc(:)
    
    integer :: i,cellA,cellB
    real*8  :: LastCelldP,maxdP
    real*8  :: diffFrmMaxdP,tol_dP
    real*8  :: hmuch
    real*8  :: PressureVal(1:N_cell)
    
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
    
    write(*,*) "EnteredInto makehigerdP"
    open(unit=56,file='make_higherdP.dat')
    
    allocate(dP_calc(1:LastCellToCmpr))
    dP_calc = 1.0d30
    
    
    do
       
       open(unit=56,file='make_higherdP.dat',position='append')
       PressureVal(1:N_cell) = Pressure(node_xy,A0)
       
       do i = 1,LastCellToCmpr
          cellA = i ; cellB = (i+1)
          dP_calc(i) = PressureVal(cellA) - PressureVal(cellB)
          write(56,*) cellA,cellB,"cell A-B"
          write(56,*) PressureVal(cellA),PressureVal(cellB),dP_calc(i),"PressureVals"
       enddo
       
       LastCelldP = dP_calc(LastCellToCmpr)
       
       maxdP = 0.00d0
       
       do i = 1,LastCellToCmpr
          
          if (dP_calc(i).gt.maxdP) then
             maxdP = dP_calc(i)
             write(56,*) maxdP,i,"maxdP and cell no"
          endif
          
       enddo
       
       diffFrmMaxdP = abs(LastCelldP-maxdP)
       tol_dP = 0.000001d0
       
       write(56,*) diffFrmMaxdP,"diffFrmMaxdP"
       
       if (diffFrmMaxdP .gt. tol_dP) then
          hmuch = 1.00d0 + (0.01)
       
          call increase_A0_ofConsecutiveCells(A0orKA,consecCells,no_of_CC,hmuch)
          call Equilibrate_system
          call save_config_and_generate_colordata(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          Frame_NI = Frame_NI+1
          
       elseif (diffFrmMaxdP .le. tol_dP) then
          exit
       endif
       
       close(56)
       
    enddo
       
  end subroutine make_higherdP_inLastCell
  
  subroutine reduce_l0_of_Consec_cells_Spr_to_retainShape(consecCells,no_of_CC,sprlenVarPrv,sprlenBfr,cellAreaVarPrv,cellAreaBfr)
    implicit none
    integer, intent(in) :: no_of_CC
    integer, intent(in) :: consecCells(1:no_of_CC)
    real*8 , intent(in) :: sprlenVarPrv(1:N_spr),sprlenBfr(1:N_spr)
    real*8 , intent(in) :: cellAreaVarPrv(1:N_cell),cellAreaBfr(1:N_cell)
    
    integer :: i,j,jmax,jmaxHlf
    integer :: Lcell,Rcell
    integer :: nspr_in_Lcell,nspr_in_Rcell
    integer :: nspr_toAdjust,ncell_toAdjust
    integer :: cnt_Adjstd
    real*8  :: hmuch
    integer :: sprLocL,sprLocR
    
    real*8  :: sprlenVarNow(1:N_spr),currl0(1:N_spr)
    real*8  :: cellAreaVarNow(1:N_cell)
    
    integer :: spr_no,spr_typ
    integer :: cell_no
    real*8  :: tol_prcnt
    logical :: lgcl_Adjstd
    
    integer, allocatable :: sprToBeAdjstd(:),cellToBeAdjstd(:)
    real*8 , allocatable :: prcnt_chng(:),prcnt_chngC(:)
    
    
    open(unit=28,file='reduce_l0_of_Consec_cells_Spr.dat',position='append')
    
    write(28,*)"start of the call"
    
    nspr_toAdjust  = 0
    ncell_toAdjust = 0
    
    do i = 1,no_of_CC
       Lcell = consecCells(i) ; Rcell = (N_cell/2)+consecCells(i)
       write(28,*) Lcell,Rcell,"L-R Cells"
       
       nspr_in_Lcell = area_spr(Lcell,0)
       nspr_in_Rcell = area_spr(Rcell,0)
       
       if (i==1) then
          nspr_toAdjust=nspr_toAdjust+(nspr_in_Lcell+nspr_in_Rcell)
       elseif (i.ne.1) then
          nspr_toAdjust=nspr_toAdjust+(nspr_in_Lcell+nspr_in_Rcell)-(2*(NAEC+1))
       endif
       
    enddo
    
    ncell_toAdjust = 2*no_of_CC
    
    allocate(sprToBeAdjstd(1:nspr_toAdjust))
    allocate(cellToBeAdjstd(1:ncell_toAdjust))
    
    cnt_Adjstd = 0
    
    do i = 1,no_of_CC
       Lcell = consecCells(i) ; Rcell = (N_cell/2)+consecCells(i)
       
       cellToBeAdjstd(2*i-1) = Lcell
       cellToBeAdjstd(2*i-0) = Rcell
       write(28,*) cellToBeAdjstd(2*i-1),cellToBeAdjstd(2*i-1),(2*i-1),(2*i),"cellToBeAdjstd"
       
       nspr_in_Lcell = area_spr(Lcell,0)
       nspr_in_Rcell = area_spr(Rcell,0)
       
       if (i==1) then
          jmax    = (nspr_in_Lcell+nspr_in_Rcell)
          jmaxHlf = (nspr_in_Lcell)
          sprLocL = 0
          sprLocR = 0
          
       elseif (i.ne.1) then
          jmax    = (nspr_in_Lcell+nspr_in_Rcell)-(2*(NAEC+1))
          jmaxHlf = (nspr_in_Lcell) - (NAEC+1)
          sprLocL = NAEC+1
          sprLocR = NAEC+1
          
       endif
       
       do j = 1,jmax
          
          if (j.le.jmaxHlf) then
             cnt_Adjstd = cnt_Adjstd+1
             sprToBeAdjstd(cnt_Adjstd) = area_spr(Lcell,(j+sprLocL))
             
          elseif (j.gt.jmaxHlf) then
             cnt_Adjstd = cnt_Adjstd+1
             sprToBeAdjstd(cnt_Adjstd) = area_spr(Rcell,(j-jmaxHlf+sprLocR))
             
          endif
          
          write(28,*) sprToBeAdjstd(cnt_Adjstd),typ_spr(sprToBeAdjstd(cnt_Adjstd)),cnt_Adjstd,"sprToBeAdjstd,cnt"
          
       enddo
       
    enddo
    
    tol_prcnt       = 10.00d0
    currl0(1:N_spr) = l0(1:N_spr)
    hmuch           = 0.05d0 ! 5 percent change
    
    allocate(prcnt_chng(1:nspr_toAdjust))
    allocate(prcnt_chngC(1:ncell_toAdjust))
    
    do 
       do i = 1,nspr_toAdjust
          spr_no = sprToBeAdjstd(i)
          l0(spr_no) = (1.00-hmuch)*l0(spr_no)
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       sprlenVarNow(1:N_spr)    = abs(l(1:N_spr)-sprlenBfr(1:N_spr))
       cellAreaVarNow(1:N_cell) = abs(A(1:N_cell)-cellAreaBfr(1:N_cell))
       
       open(unit=29,file='reduce_l0forprcntChnge.dat',position='append')
       write(29,*) (Frame_NI-1),"Frame_NI"
       
       do i = 1,nspr_toAdjust
          
          spr_no = sprToBeAdjstd(i)
          !write(*,*) sprlenVarNow(spr_no),sprlenVarPrv(spr_no),spr_no,"sprlenVar now,prv,spr_no"
          prcnt_chng(i) = (abs(l(spr_no)-sprlenBfr(spr_no)) / sprlenBfr(spr_no)) * (100.00)
          write(*,*)  prcnt_chng(i),i,"prcnt"
          write(29,*) prcnt_chng(i),l(spr_no),sprlenBfr(spr_no),i,spr_no,"spr"
          
       enddo
       
       write(29,*) " "
       
       do i= 1,ncell_toAdjust
          
          cell_no = cellToBeAdjstd(i)
          write(*,*) cellAreaVarNow(cell_no),cellAreaVarPrv(cell_no),cell_no,"cellAreaVar"
          prcnt_chngC(i) = (abs(A(cell_no)-cellAreaBfr(cell_no)) / cellAreaBfr(cell_no)) * (100.00)
          write(29,*) prcnt_chngC(i),A(cell_no),cellAreaBfr(cell_no),i,cell_no,"cell"
          
       enddo
       
       close(29)
       
       lgcl_Adjstd = .False.
       
       !%%%%% (CHEKING FOR TOLERATING CONDITIONS OVER SPRINGS) %%%%%!
       
       ! do i = 1,nspr_toAdjust
          
       !    spr_no  = sprToBeAdjstd(i)
       !    spr_typ = typ_spr(spr_no)
          
       !    if (spr_typ==1 .or. spr_typ==2 .or. spr_typ==3 .or. spr_typ==4) then
             
       !       if (prcnt_chng(i) .le. tol_prcnt) then
       !          lgcl_Adjstd = .True.
       !       elseif (prcnt_chng(i) .gt. tol_prcnt) then
       !          lgcl_Adjstd = .False.
       !          exit
       !       endif
             
       !    endif
          
       ! enddo
       
       ! if (lgcl_Adjstd .eqv. .True.) then
       !    exit
       ! elseif (lgcl_Adjstd .eqv. .False.) then
       !    continue
       ! endif
       
       !%%%%% (CHEKING FOR TOLERATING CONDITIONS OVER CELLS) %%%%%!
       
       do i = 1,ncell_toAdjust
          
          cell_no = cellToBeAdjstd(i)
          
          if (prcnt_chngC(i) .le. tol_prcnt) then
             lgcl_Adjstd = .True.
          elseif (prcnt_chngC(i) .gt. tol_prcnt) then
             lgcl_Adjstd = .False.
             exit
          endif
          
       enddo
       
       if (lgcl_Adjstd .eqv. .True.) then
          exit
       elseif (lgcl_Adjstd .eqv. .False.) then
          continue
       endif
       
       
    enddo
    
    
    close(28)
    
  end subroutine reduce_l0_of_Consec_cells_Spr_to_retainShape
  
  
  subroutine maintain_precise_dP_in_bfrInvaginatingCells(seqNo,LastCellToCmpr)
    implicit none
    integer, intent(in) :: seqNo,LastCellToCmpr
    
    integer :: N_dPs ! Number of deltaP's 
    real*8  :: PressureVal(1:N_cell)
    integer :: cell_A,cell_B ! Cell_A=Cell Ahead and Cell_B=Cell Behind
    
    real*8  :: dP_valDFTC,dP_valSBTC,dP_valDPTC
    integer :: PosOrNeg
    integer :: i,j
    real*8  :: tol_dP,tol_dPprcnt,maxdP
    real*8  :: ZERO
    
    real*8, allocatable :: dP_rltv(:)
    real*8, allocatable :: dP_valNow(:),dP_valShdBe(:)
    real*8, allocatable :: dP_valDiff(:),dP_valDiffPrcnt(:)
    
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
    
    PressureVal(1:N_cell) = Pressure(node_xy,A0)
    
    open(unit=43,file='maintain_dP.dat',position='append')
    
    maxdP = PressureVal(LastCellToCmpr) - PressureVal(LastCellToCmpr+1)
    N_dps = LastCellToCmpr
    ZERO  = 0.00d0
    
    allocate(dP_rltv(1:N_dPs))
    allocate(dP_valNow(1:N_dPs),dP_valShdBe(1:N_dPs))
    allocate(dP_valDiff(1:N_dPs),dP_valDiffPrcnt(1:N_dPs))
    
    
    if (seqNo==9) then
       
       dP_rltv(1:2)               = 0.00d0
       dP_rltv(3)                 = 0.25d0
       dP_rltv(4)                 = 0.50d0
       dP_rltv((N_dPs-1):(N_dPs)) = 1.00d0
       
    elseif (seqNo==10) then
       
       dP_rltv(1)                 = 0.00d0
       dP_rltv(2)                 = 0.25d0
       dP_rltv(3)                 = 0.50d0
       dP_rltv((N_dPs-1):(N_dPs)) = 1.00d0
       
    endif
    
    
    do i = 1,N_dPs
       cell_A = (i+1) ; cell_B = i
       write(43,*) cell_A,cell_B,"cell ahead and behind"
       
       dP_valNow(i)   = PressureVal(cell_B) - PressureVal(cell_A)
       dP_valShdBe(i) = (maxdP)*(dP_rltv(i)) 
       dP_valDiff(i)  = dP_valNow(i)-dP_valShdBe(i)
       
       write(*,*) dP_valNow(i),dP_valShdBe(i),dP_valDiff(i),i,"dP_vals"
       dP_valDiffPrcnt(i) = (abs(dP_valNow(i)-dP_valShdBe(i))/dP_valNow(i))*100.00
       
       write(43,*) maxdP,dP_rltv(i),PressureVal(cell_B),PressureVal(cell_A),"maxdP"
       write(43,*) dP_valNow(i),dP_valShdBe(i),dP_valDiff(i),dP_valDiffPrcnt(i),i,"dP_vals"
       
    enddo
    
    close(43)
    open(unit=43,file='maintain_dP.dat',position='append')
    
    tol_dP      = 1.0d-3
    tol_dPprcnt = 10.00d0
    
    do i = (N_dPs-1),1,-1
       
       dP_valDFTC = dP_valDiff(i)      ! DFTC=Difference For This Cell
       dP_valSBTC = dP_valShdBe(i)     ! SBTC=Should Be for This Cell
       dP_valDPTC = dP_valDiffPrcnt(i) ! DPTC=Difference Percent fot This Cell
       
       if ( (dP_valDiff(i).gt.ZERO) .and. (dP_valDiffPrcnt(i).gt.tol_dPprcnt) ) then
          
          PosOrNeg   = -1
          call pressureAdjstmnt_forPrecisedP(i,dP_valDFTC,dP_valSBTC,dP_valDPTC,PosOrNeg,tol_dPprcnt)
          
          ! %update vals 
          
          PressureVal(1:N_cell) = Pressure(node_xy,A0)
          
          do j = 1,N_dPs
             cell_A = (j+1) ; cell_B = j
             !write(43,*) cell_A,cell_B,"cell ahead and behind"
             
             dP_valNow(j)   = PressureVal(cell_B) - PressureVal(cell_A)
             dP_valDiff(j)  = dP_valNow(i)-dP_valShdBe(j)   
             !write(43,*) dP_valNow(i),dP_valShldBe(i),dP_valDiff(i),i,"dP_vals"
             
          enddo
          
       elseif ( (dP_valDiff(i).gt.ZERO) .and. (dP_valDiffPrcnt(i).gt.tol_dPprcnt) ) then
          continue
       elseif ( (dP_valDiff(i).le.ZERO) .and. (dP_valDiffPrcnt(i).gt.tol_dPprcnt) ) then
          
          PosOrNeg   = 1
          call pressureAdjstmnt_forPrecisedP(i,dP_valDFTC,dP_valSBTC,dP_valDPTC,PosOrNeg,tol_dPprcnt)
          
          ! %update vals
          
          PressureVal(1:N_cell) = Pressure(node_xy,A0)
          
          do j = 1,N_dPs
             cell_A = (j+1) ; cell_B = j
             !write(43,*) cell_A,cell_B,"cell ahead and behind"
             
             dP_valNow(j)   = PressureVal(cell_B) - PressureVal(cell_A)
             dP_valDiff(j)  = dP_valNow(j)-dP_valShdBe(j)   
             !write(43,*) dP_valNow(i),dP_valShldBe(i),dP_valDiff(i),i,"dP_vals"
          enddo
          
       elseif ( (dP_valDiff(i).le.ZERO) .and. (dP_valDiffPrcnt(i).gt.tol_dPprcnt) ) then
          continue
       endif
       
    enddo
    
    
    write(43,*) PressureVal(1:(LastCellToCmpr+1)),"PressureAft"
    write(43,*) dP_valDiff(1:N_dPs),"dP_valDiff"
    
    close(43)
    
  end subroutine maintain_precise_dP_in_bfrInvaginatingCells
  
  
  subroutine pressureAdjstmnt_forPrecisedP(dPsNo,dP_valDFTC,dP_valSBTC,dP_valDPTC,PosOrNeg,tol_dPprcnt)
    implicit none
    integer, intent(in)    :: dPsNo
    real*8 , intent(inout) :: dP_valDFTC,dP_valDPTC
    real*8 , intent(in)    :: dP_valSBTC
    integer, intent(in)    :: PosOrNeg
    real*8 , intent(in)    :: tol_dPprcnt
    
    real*8  :: StrtdP_valNFTC,StrtdP_valDFTC,prdct_of_DFTC,ZERO
    integer :: A0orKA,no_of_CC
    integer :: L1Cell,L2Cell,R1Cell,R2Cell
    
    integer, allocatable  :: consecCells(:)
    
    real*8  :: sprlenBfr(1:N_spr),sprlenVar(1:N_spr)
    real*8  :: cellAreaBfr(1:N_cell),cellAreaVar(1:N_cell)
    real*8  :: hmuch,dP_valNFTC 
    real*8  :: PressureVal(1:N_cell)
    
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
    
    open(unit=27,file='Frame_NI.dat',position='append')
    open(unit=28,file='prcnt_indP.dat',position='append')
    
    L1Cell = dPsNo ; L2Cell = dPsNo+1
    R1Cell = (N_cell/2) + L1Cell ; R2Cell = (N_cell/2) + L2Cell
    
    sprlenBfr(1:N_spr)    = l(1:N_spr)
    cellAreaBfr(1:N_cell) = A(1:N_cell)
    
    A0orKA   = 1
    no_of_CC = dPsNo
    ZERO     = 0.0d0
    
    allocate(consecCells(1:no_of_CC))
    
    do i = 1,no_of_CC
       consecCells(i) = i
    enddo
    
    
    if (PosOrNeg==1) then
       hmuch = 1.00d0 + (0.01d0) !1 prcnt
    elseif (PosOrNeg==-1) then
       hmuch = 1.00d0 - (0.01d0)
    endif
    
    StrtdP_valNFTC = PressureVal(dPsNo) - PressureVal(dPsNo+1)
    StrtdP_valDFTC = StrtdP_valNFTC - dP_valSBTC
    
    do
       
       call increase_A0_ofConsecutiveCells(A0orKA,consecCells,no_of_CC,hmuch)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       write(27,*) Frame_NI,"in Do loop Precise"
       Frame_NI = Frame_NI+1
       
       PressureVal(1:N_cell) = Pressure(node_xy,A0)
       
       dP_valNFTC = PressureVal(dPsNo) - PressureVal(dPsNo+1) !NFTC =Now for this cell
       dP_valDFTC = dP_valNFTC - dP_valSBTC
       dP_valDPTC = (abs(dP_valNFTC - dP_valSBTC)/dP_valNFTC)*100.00d0
       
       write(28,*) dP_valNFTC,dP_valDFTC,dP_valDPTC,"dP_NFTC,DFTC,DPTC"
       
       prdct_of_DFTC = (StrtdP_valDFTC)*(dP_valDFTC)
       
       if ((dP_valDPTC.lt.tol_dPprcnt) .or. (prdct_of_DFTC.lt.ZERO)) then
          write(27,*) Frame_NI,"Ending Precise"
          
          sprlenVar(1:N_spr)    = abs(l(1:N_spr)-sprlenBfr(1:N_spr))
          cellAreaVar(1:N_cell) = abs(A(1:N_cell)-cellAreaBfr(1:N_cell))
          
          !call reduce_l0_of_Consec_cells_Spr_to_retainShape(consecCells,no_of_CC,sprlenVar,sprlenBfr,cellAreaVar,cellAreaBfr)
          !PressureVal(1:N_cell) = Pressure(node_xy,A0)
          exit
          
          
       endif
       
    enddo
    
    close(27)
    close(28)
    
  end subroutine pressureAdjstmnt_forPrecisedP
  
  
  !!!*** ROUTINES FOR CHECKING IF PRESSURE DIFFERENCES ARE NEEDED (STARTING POINT) ***!!!
  
  
  subroutine ForcePrtrbtnAnalysis_atCS2(strctToRead,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(inout) :: ExpNo
    integer, intent(inout) :: FrmNo
    
    
    real*8,  allocatable :: PresVal(:),TnsnVal(:)
    real*8,  allocatable :: kaV(:),A0V(:),Aval(:)
    real*8,  allocatable :: ksV(:),l0V(:),lval(:)
    
    integer, allocatable :: sprsLft(:),sprsRgt(:)
    integer, allocatable :: sprsLftD(:),sprsRgtD(:)
    integer, allocatable :: sprsBtmLft(:),sprsBtmRgt(:)
    integer, allocatable :: sprsBtmLftD(:),sprsBtmRgtD(:)
    
    real*8  :: PresToBeMatched
    real*8  :: AvalsBeMatched(1:2)
    logical :: lgcl_matched,lgcl_Aval
    integer :: countIndvd(1:3)
    real*8  :: E
    real*8  :: tolrncePress
    integer :: nsprs_invgntingCells,nsprs_invgnatedCells
    integer :: nsprs_btmCell1,nsprs_btmCell2,nsprs_btmCell3
    integer :: nsprs_BTC !BTC=Both Type of (Invgnting) Cells
    integer :: nsprs_ETC !ETC=End Type of (last three bottom pair) Cells
    
    integer :: flnmGen
    real*8  :: hmuchEach,hmuchCgX
    integer :: cellChoice
    integer :: pullingOrpushing
    integer :: RedOrRlx
    real*8  :: hmuchl0
    integer :: howmanyT,ManpltTech
    integer :: flnmInit,flnmFinl
    integer :: whatPropsToReadforExtrplt
    
    integer :: i,j,jmax
    real*8  :: distnce
    integer :: count
    real*8  :: distance
    real*8  :: loopIncr
    integer :: incrORdcrs
    
    integer :: caseNo
    integer :: ExpNoCase,FrmNoCase
    integer :: NonUniOrUni,TnsnOrExtrpltn
    integer :: cellNm
    
    real*8, allocatable :: ks_stre(:)
    
    open(unit=168,file='parametersChk_ip_invagniating_invaginated_cell_rs_atCS2.dat')
    
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
    
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(*,*) (FrmNo-1),"bfr changing A0 of invgnting and invgnated cells"
    
    !!!!!!!!!  HERE  !!!!!!!!!
    
    call find_the_invagniating_and_invaginated_cell
    call find_the_frst_four_cellPair
    call find_the_last_three_cellPair
    
    nsprs_invgntingCells = area_spr(invgnting_cells(1),0)
    nsprs_invgnatedCells = area_spr(invgnated_cells(1),0)
    nsprs_BTC            = (nsprs_invgntingCells+nsprs_invgnatedCells)-(NAEC_Ltrl+1)
    
    allocate(sprsLft(1:nsprs_BTC),sprsRgt(1:nsprs_BTC))
    allocate(sprsLftD(1:nsprs_BTC),sprsRgtD(1:nsprs_BTC))
    sprsLftD(1:nsprs_BTC) = 1 ; sprsRgtD(1:nsprs_BTC) = 1
    call find_sprsList_InvgntingInvgnatedCells(nsprs_BTC,sprsLft,sprsRgt)
    
    nsprs_btmCell1 = area_spr(thrdlastCells(1),0)
    nsprs_btmCell2 = area_spr(scndlastCells(1),0)
    nsprs_btmCell3 = area_spr(frstlastCells(1),0)
    nsprs_ETC      = (nsprs_btmCell1)+(nsprs_btmCell2)+(nsprs_btmCell3) - 2*(NAEC_Ltrl+1)
    
    allocate(sprsBtmLft(1:nsprs_ETC),sprsBtmRgt(1:nsprs_ETC))
    allocate(sprsBtmLftD(1:nsprs_ETC),sprsBtmRgtD(1:nsprs_ETC))
    sprsBtmLftD(1:nsprs_ETC) = 1 ; sprsBtmRgtD(1:nsprs_ETC) = 1
    call find_sprsList_BottmCells(nsprs_ETC,sprsBtmLft,sprsBtmRgt)
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    call store_props_bfr_chng(kaV,A0V,Aval,ksV,l0V,lval)
    
    !! FROM HERE COMMENTING STARTS !!
    
    pullingOrpushing = 1
    caseNo           = 3
    call get_CgVal_ExpNo_and_FrmNo(pullingORpushing,caseNo,hmuchCgX,ExpNoCase,FrmNoCase)
    
    NonUniorUni    = 1
    TnsnOrExtrpltn = 1 !1 for DM tension, 2 for Extrpltn
    
    if (NonUniOrUni==1) then ! 1 for NonUnifrmPres, 2 for UnifrmPres
       
       write(*,*) "starting point (BEFORE APPLYING FORCE)"
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
       FrmNoCase = FrmNoCase+1
       
       !stop
       
       call pulling_or_pushing(pullingOrpushing,hmuchCgX)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
       FrmNoCase = FrmNoCase+1
       
       call find_distance_frm_pulley(distance)
       if (distance .le. 0.05d0) then
          write(*,*) "This decreasing MIDLINE APCL tension block is happening for CaseNo =",caseNo
          write(*,*) (FrmNoCase-1),"FrmNo bfr calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
          incrORdcrs = 2
         call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
       endif
       
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING DIAG TNSN (BEFORE)"
       call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNoCase,FrmNoCase)
       write(*,*) (FrmNoCase-1),"WE HAVE CHECKED DIAG TNSN (AFTER)"
       
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIGHTNESS (BEFORE)"
       call make_diagonal_spr_more_strght(ExpNoCase,FrmNoCase)
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIHTNESS (AFTER)"
       
       
       if (TnsnOrExtrpltn == 1) then
          write(*,*) (FrmNo-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=1 "
          incrORdcrs = 1
          call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=1"
          
       elseif (TnsnOrExtrpltn == 2) then
          flnmInit = strctToRead-1 ; flnmFinl = strctToRead
          ExpNo = ExpNoCase ; FrmNo = FrmNoCase ! imp for Extrpltn
          call prep_vars_for_Extrplt_loctn1(NonUniorUni,flnmInit,flnmFinl)
       endif
       
    elseif (NonUniorUni==2) then
       
       cellChoice = 1 ; call get_PresValMatched(PresVal,cellChoice,PresToBeMatched)
       call get_AvalsMatched(AvalsBeMatched)
       
       ! do 
       
       !    call increasingA0s_of_invgnAndinvgntd_Cells(PresVal,PresToBeMatched)
       !    call Equilibrate_system
       !    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       !    FrmNo = FrmNo+1
       
       !    RedOrRlx=1 ; call RedRlx_l0s_selected_sprs(RedOrRlx,sprsLft,sprsRgt,sprsLftD,sprsRgtD,nsprs_BTC)
       !    call Equilibrate_system
       !    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       !    FrmNo = FrmNo+1
       
       !    call find_Pressure_and_Tension(PresVal,TnsnVal)
       !    call find_if_PresVal_matched(PresVal,PresToBeMatched,lgcl_matched)
       !    call find_if_l0sNeededToBeReduced(AvalsBeMatched,sprsLftD,sprsRgtD,nsprs_BTC,lgcl_Aval)
       
       !    write(*,*) (FrmNo-1),"FrmNo"   
       
       !    if (lgcl_matched .eqv. .True.) then
       !       flnmGen=4 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       !       exit
       !    endif
       
       ! enddo
       
       flnmGen=4 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"POSITION VALUE 1"
       
       write(*,*) "START DECREASING A0s OF BOTTOM THREE"
       
       ! call find_Pressure_and_Tension(PresVal,TnsnVal)
       ! lgcl_matched = .False.
       
       ! do
       !    call decreasingA0s_of_lastthreePair_Cells(PresVal,PresToBeMatched)
       !    call Equilibrate_system
       !    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       !    FrmNo = FrmNo+1
       
       !    RedOrRlx=2 ; call RedRlx_l0s_selected_sprs(RedOrRlx,sprsBtmLft,sprsBtmRgt,sprsBtmLftD,sprsBtmRgtD,nsprs_ETC)
       !    call Equilibrate_system
       !    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       !    FrmNo = FrmNo+1
       
       !    call find_Pressure_and_Tension(PresVal,TnsnVal)
       !    call find_if_PresVal_matchedDecr(PresVal,PresToBeMatched,lgcl_matched,countIndvd)
       
       !    if (lgcl_matched .eqv. .True.) then
       !       flnmGen=14 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       !       exit
       !    endif
       
       ! enddo
       
       write(*,*) "START REDUCING THE REST LENGTHS OF APCL FROM INVGNTED CELLS"
       
       flnmGen=14 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"POSITION VALUE 2"
       
       !ManpltTech=1 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
       !flnmGen=17   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       flnmGen=17 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"Bfr Applying Force in CS2"
       write(*,*) (FrmNo-1),"POSITION VALUE 3"
       
       ks_stre(1:N_spr) = k_spr(1:N_spr)
       
       ! count = 1
       
       ! do 
       !    k_spr(50) = (1.25d0+(count-1)*0.05d0)*ks_stre(50)
       !    k_spr(127) = (1.25d0+(count-1)*0.05d0)*ks_stre(127)
       !    k_spr(57) = (1.25d0+(count-1)*0.05d0)*ks_stre(57)
       !    k_spr(134) = (1.25d0+(count-1)*0.05d0)*ks_stre(134)
       !    k_spr(64) = (1.25d0+(count-1)*0.05d0)*ks_stre(64)
       !    k_spr(141) = (1.25d0+(count-1)*0.05d0)*ks_stre(141)
       
       !    call Equilibrate_system
       !    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       !    FrmNo = FrmNo+1
       !    write(*,*) (FrmNo-1),"FrmNo in loop"
       
       !    call find_distance_frm_pulley(distance)
       
       !    if (distance .le. 0.05d0) then
       !       write(*,*) count,"count"
       !       exit
       !    endif
       
       !    count = count+1
       ! enddo
       
       !flnmGen=26   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       
       write(*,*) "Reading flnmGen=26"
       
       flnmGen=26 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"POSITION VALUE 4"
       
       !ManpltTech=1 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
       !flnmGen=27   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       
       flnmGen=27 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"POSITION VALUE 5"
       
       !if ((FrmNo-1) .gt. 14) stop
       
       write(*,*) "WE ARE ADJUSING DIAGONAL TENSION IF NEEDED, Frm Before Adjusting is =",(FrmNo-1)
       call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNo,FrmNo)
       write(*,*) "WE ARE ADJUSING DIAGONAL TENSION IF NEEDED, Frm After Adjusting is =",(FrmNo-1)
       flnmGen=28  ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       
       write(*,*) "WE ARE GETTING UNIFORM PRES AFT THIS (NO FORCE), Frame Num will be =",(FrmNo-1)
          
       ! write(*,*) "NOW WE WILL SEE IF WE RELEASE TSNN IN MIDLINE APCL, CELLS WILL GO FAR FROM MEETING",(FrmNo-1)
       
       ! call find_distance_frm_pulley(distance)
       ! if (distance .le. 0.05d0) then
       !    write(*,*) (FrmNo-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=2 "
       !    incrORdcrs = 2
       !    call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNo,FrmNo)
       !    write(*,*) (FrmNo-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=2 "
       ! endif
       
       ! write(*,*) "CELLS WENT FAR FROM MEETING AFT RELEASING MIDLINE APCL",(FrmNo-1)
       ! write(*,*) "WE WILL REGAIN STATE BFR MIDLINE RELEASING"
       
       ! flnmGen=28 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       ! call Equilibrate_system
       ! call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       ! FrmNo = FrmNo+1
       
       !write(*,*) "AFTER REGAINING STATE WITH MIDLINE RELEASING",(FrmNo-1)
       
       !write(*,*) (FrmNo-1),"Frame No before Pulling at CgX=0.15d0"
       
       write(*,*) "starting point (BEFORE APPLYING FORCE)"
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
       FrmNoCase = FrmNoCase+1
       
       !stop
       
       call pulling_or_pushing(pullingOrpushing,hmuchCgX)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
       FrmNoCase = FrmNoCase+1
       
       call find_distance_frm_pulley(distance)
       if (distance .le. 0.05d0) then
          write(*,*) "This decreasing MIDLINE APCL tension block is happening for CaseNo =",caseNo
          write(*,*)(FrmNoCase-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
          incrORdcrs = 2
          call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
       endif
       
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING DIAG TNSN (BEFORE)"
       call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNoCase,FrmNoCase)
       write(*,*) (FrmNoCase-1),"WE HAVE CHECKED DIAG TNSN (AFTER)"
       
       
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIGHTNESS (BEFORE)"
       call make_diagonal_spr_more_strght(ExpNoCase,FrmNoCase)
       write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIHTNESS (AFTER)"
       
       flnmGen=32 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       write(*,*) (FrmNoCase-1),"FrmNoCase"
       
       ! flnmGen=31   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       ! write(*,*) (FrmNo-1),"Frame No"
       ! stop
       
       if (TnsnOrExtrpltn == 1) then
          write(*,*) (FrmNo-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=1 "
          incrORdcrs = 1
          call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=1"
          
       elseif (TnsnOrExtrpltn == 2) then
          flnmInit = 25  ; flnmFinl = 32
          ExpNo = ExpNoCase ; FrmNo = FrmNoCase ! imp for Extrpltn
          call prep_vars_for_Extrplt_loctn1(NonUniorUni,flnmInit,flnmFinl)
       endif
       
    endif
    
    
    !!COMMENT START
    
    ! hmuchl0 = 0.05d0 ; howmanyT = 2
    ! call reducing_l0s_of_diagonal_spr(hmuchl0,howmanyT,ExpNo,FrmNo)
    ! flnmGen=9 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    ! write(*,*) (FrmNo-1),"FrmNo Bfr l0+ks"
    
    ! hmuchEach=0.10d0;call reducing_l0s_incrsing_ks_apcl_frm_invgnted_Cells(hmuchEach,ExpNo,FrmNo)
    ! flnmGen=12      ;call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    ! write(*,*) (FrmNo-1),"FrmNo Bfr Neg2"
    ! call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNo,FrmNo)
    ! write(*,*) (FrmNo-1),"FrmNo Aft Neg2"
    
    !!COMMENT END
    
    close(168)
    
    if (TnsnOrExtrpltn==1) then
       
       if (NonUniorUni==1) then
          call incoming_region_pressure_drop(ExpNoCase,FrmNoCase)
       elseif (NonUniorUni==2) then
          
          cellNm = Hlf_Ncell-1
          call tilted_rectngle_to_rectngleMatchedwithPage(cellNm,ExpNoCase,FrmNoCase)
          write(*,*) (FrmNo-1),"FrmNo before calling manplt_apcl_membrne_tnsn with incrORdcrs=2 "
          incrORdcrs = 2
          call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
          call find_distance_frm_pulley(distance)
          write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=1"
          
          stop
          
       endif
       
    elseif (TnsnOrExtrpltn==2) then
       continue
    endif
      
  end subroutine ForcePrtrbtnAnalysis_atCS2
  
  
  subroutine incr_press_invaginating_invaginated_cell_and_retainShape_atCS1(strctToRead,ExpNo,FrmNo)
    implicit none
    integer, intent(in)  :: strctToRead
    integer, intent(in)  :: ExpNo
    integer, intent(out) :: FrmNo
    
    real*8,  allocatable :: PresVal(:),TnsnVal(:)
    real*8,  allocatable :: kaV(:),A0V(:),Aval(:)
    real*8,  allocatable :: ksV(:),l0V(:),lval(:)
    
    real*8  :: E
    integer :: flnmGen
    real*8  :: hmuchEach,hmuchCgX
    integer :: cellChoice
    integer :: pullingOrpushing
    integer :: RedOrRlx
    real*8  :: hmuchl0
    integer :: howmanyT,ManpltTech
    
    open(unit=158,file='parametersChk_ip_invagniating_invaginated_cell_rs_atCS1.dat')
    
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
    
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(*,*) (FrmNo-1),"bfr chaging A0 of invgnting and invgnated cells"
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    call store_props_bfr_chng(kaV,A0V,Aval,ksV,l0V,lval)
    
    !ManpltTech=2 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=19   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=19 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(*,*) "After Reading" 
    
    
    !ManpltTech=2 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=20   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=20 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(*,*) (FrmNo-1),"After reading 20 (whatToFix=2)" 
    
    !ManpltTech=2 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=21   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=21 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    !ManpltTech=2 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=22   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=22 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    !ManpltTech=2 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=23   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=23 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(*,*) "After reading flnmGen=23", (FrmNo-1)
    
    !ManpltTech=1 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=24   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    flnmGen=24 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(*,*) "After reading flnmGen = 24", (FrmNo-1)
    
    !ManpltTech=1 ; call find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    !flnmGen=30   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    
    !write(*,*) "After reading flnmGen=30", (FrmNo-1)
    
    !flnmGen=30 ; call read_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    !call Equilibrate_system
    !call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    !FrmNo = FrmNo+1
    
    !write(*,*) "After getting flnmGen=30", (FrmNo-1)
    
    
    l0(68)  = l0(68)*0.94d0  ; l0(69)  = l0(69)*0.94d0  ; l0(70)  = l0(70)*0.94d0
    l0(145) = l0(145)*0.94d0 ; l0(146) = l0(146)*0.94d0 ; l0(147) = l0(147)*0.94d0
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    flnmGen=25   ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
    write(*,*) FrmNo-1,"FrmNo aft CS1"
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(158)
    
  end subroutine incr_press_invaginating_invaginated_cell_and_retainShape_atCS1
  
  
  subroutine loop_txt_for_incrPres_and_reducel0_separately(nsprs_BTC)
    implicit none
    integer, intent(in) :: nsprs_BTC
    
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    real*8  :: PresToBeMatched
    integer :: ExpNo,FrmNo
    logical :: lgcl_matched,lgcl_Aval
    
    integer, allocatable :: sprsLft(:),sprsRgt(:),sprsLftD(:),sprsRgtD(:)
    real*8  :: AvalsBeMatched(1:2)
    integer :: RedOrRlx
    integer :: flnmGen !! COPY ONLY THE LOOP : VARIABLES ONLY FOR REMOVING ERROR
    
    allocate(sprsLft(1:nsprs_BTC),sprsRgt(1:nsprs_BTC),sprsLftD(1:nsprs_BTC),sprsRgtD(1:nsprs_BTC))
    
    !!! COPY FROM BELOW !!!
    
    do 
       
       call increasingA0s_of_invgnAndinvgntd_Cells(PresVal,PresToBeMatched)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       call find_Pressure_and_Tension(PresVal,TnsnVal)
       call find_if_PresVal_matched(PresVal,PresToBeMatched,lgcl_matched)
       write(*,*) (FrmNo-1),"FrmNo"
       
       
       if (lgcl_matched .eqv. .True.) then
          flnmGen=1 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
       endif
       
       if (lgcl_matched .eqv. .True.) then
          
          call get_AvalsMatched(AvalsBeMatched)
          
          do
             RedOrRlx=1 ; call RedRlx_l0s_selected_sprs(RedOrRlx,sprsLft,sprsRgt,sprsLftD,sprsRgtD,nsprs_BTC)
             call Equilibrate_system
             call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
             FrmNo = FrmNo+1
             
             call find_if_l0sNeededToBeReduced(AvalsBeMatched,sprsLftD,sprsRgtD,nsprs_BTC,lgcl_Aval)
             if (lgcl_Aval.eqv..True.) then
                flnmGen=2 ; call save_props_aft_PresIncr(flnmGen,A0,k_area,l0,k_spr,CgXNode,CgYNode)
                exit
             endif
          enddo
          
          exit
          
       endif
       
    enddo
    
    write(*,*) A0(invgnting_cells(1)),A0(invgnting_cells(2)),"invgnting A0"
    write(*,*) A0(invgnated_cells(1)),A0(invgnated_cells(2)),"invgnated A0"
    write(*,*) PresVal(invgnting_cells(1)),PresVal(invgnting_cells(2)),"Pres ting cells"
    write(*,*) PresVal(invgnated_cells(1)),PresVal(invgnated_cells(2)),"Pres ted cells"
    
  end subroutine loop_txt_for_incrPres_and_reducel0_separately
  
  subroutine find_the_invagniating_and_invaginated_cell
    implicit none
    integer :: i
    integer :: node_cnnctd
    
    do i = 1,N_cell
       
       node_cnnctd = area_node(i,0)
       
       if (node_cnnctd == max_node_area) then
          invgnting_cells(1)     = i-1        ; invgnting_cells(2)     = i-1+Hlf_Ncell
          invgnated_cells(1)     = i+0        ; invgnated_cells(2)     = i+0+Hlf_Ncell
          nxt_invgnated_cells(1) = i+1        ; nxt_invgnated_cells(2) = i+1+Hlf_Ncell
          
          write(168,*) invgnting_cells(1:2),"invgnting_cells"
          write(168,*) invgnated_cells(1:2),"invgnated_cells"
          write(168,*) nxt_invgnated_cells(1:2),"nxt_of_invgnated_cells"
          exit
       endif
       
    enddo
    
  end subroutine find_the_invagniating_and_invaginated_cell
  
  subroutine find_the_frst_four_cellPair
    implicit none
    
    frstCells(1) = 1 ; frstCells(2) = 1+Hlf_Ncell
    scndCells(1) = 2 ; scndCells(2) = 2+Hlf_Ncell
    thrdCells(1) = 3 ; thrdCells(2) = 3+Hlf_Ncell
    frthCells(1) = 4 ; frthCells(2) = 4+Hlf_Ncell
    
    write(*,*) frstCells(1:2),"First  Cell Pair"
    write(*,*) scndCells(1:2),"Second Cell Pair"
    write(*,*) thrdCells(1:2),"Third  Cell Pair"
    
  end subroutine find_the_frst_four_cellPair

  subroutine find_the_last_three_cellPair
    implicit none
    
    thrdlastCells(1) = Hlf_Ncell-2 ; thrdlastCells(2) = Hlf_Ncell-2+Hlf_Ncell
    scndlastCells(1) = Hlf_Ncell-1 ; scndlastCells(2) = Hlf_Ncell-1+Hlf_Ncell
    frstlastCells(1) = Hlf_Ncell-0 ; frstlastCells(2) = Hlf_Ncell-0+Hlf_Ncell
    
  end subroutine find_the_last_three_cellPair
  
  subroutine find_Pressure_and_Tension(PresVal,TnsnVal)
    implicit none
    real*8, intent(out) :: PresVal(1:N_cell)
    real*8, intent(out) :: TnsnVal(1:N_spr)
    
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
         real*8  :: Pressure(1:N_cell)
         real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8  :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    
  end subroutine find_Pressure_and_Tension
  
  subroutine get_the_diagSprs_and_strTheDiagPrps(diagsprL,diagsprR,l0_diagsprLstr,l0_diagsprRstr,&
         ks_diagsprLstr,ks_diagsprRstr)
    implicit none
    integer, intent(out) :: diagsprL(1:(NAEC_Ltrl+1))      ,diagsprR(1:(NAEC_Ltrl+1))
    real*8 , intent(out) :: l0_diagsprLstr(1:(NAEC_Ltrl+1)),l0_diagsprRstr(1:(NAEC_Ltrl+1))
    real*8 , intent(out) :: ks_diagsprLstr(1:(NAEC_Ltrl+1)),ks_diagsprRstr(1:(NAEC_Ltrl+1))
    
    integer :: i,j
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    
    do i = 1,(NAEC_Ltrl+1)
       
       diagsprL(i) = (invgnting_cells(1)-1)*(nsprsInACell)    + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       diagsprR(i) = diagsprL(i) + (nsprsInACell)*(Hlf_Ncell) ;  write(*,*) diagsprL(i),diagsprR(i),i,"diag L-R"
       
       l0_diagsprLstr(i) = l0(diagsprL(i))     ; l0_diagsprRstr(i) = l0(diagsprR(i))
       ks_diagsprLstr(i) = k_spr(diagsprL(i))  ; ks_diagsprRstr(i) = k_spr(diagsprR(i))
       
       write(*,*) l0_diagsprLstr(i),l0_diagsprRstr(i),i,"l0 L/R"
       write(*,*) ks_diagsprLstr(i),ks_diagsprRstr(i),i,"ks L/R"
       
    enddo
    
  end subroutine get_the_diagSprs_and_strTheDiagPrps
  
  subroutine set_diagSprPrp_back(diagsprL,diagsprR,l0_diagsprLstr,l0_diagsprRstr,&
         ks_diagsprLstr,ks_diagsprRstr)
    implicit none
    integer, intent(in) :: diagsprL(1:(NAEC_Ltrl+1))      ,diagsprR(1:(NAEC_Ltrl+1))
    real*8 , intent(in) :: l0_diagsprLstr(1:(NAEC_Ltrl+1)),l0_diagsprRstr(1:(NAEC_Ltrl+1))
    real*8 , intent(in) :: ks_diagsprLstr(1:(NAEC_Ltrl+1)),ks_diagsprRstr(1:(NAEC_Ltrl+1))
    
    integer :: i,j
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    
    do i = 1,(NAEC_Ltrl+1)
       l0(diagsprL(i))    = l0_diagsprLstr(i) ; l0(diagSprR(i))    = l0_diagsprRstr(i)
       k_spr(diagsprL(i)) = ks_diagsprLstr(i) ; k_spr(diagSprR(i)) = ks_diagsprRstr(i)
       
       write(*,*) k_spr(diagsprL(i)),k_spr(diagsprR(i)),i,"k_spr in set"
       write(*,*) l0(diagsprL(i))   ,l0(diagsprR(i))   ,i,"l0    in set"
       
    enddo
    
    
  end subroutine set_diagSprPrp_back
  
  subroutine store_props_bfr_chng(kaV,A0V,Aval,ksV,l0V,lval)
    implicit none
    real*8, intent(inout) :: kaV(1:N_cell),A0V(1:N_cell),Aval(1:N_cell)
    real*8, intent(inout) :: ksV(1:N_spr) ,l0V(1:N_spr) ,lval(1:N_spr)
    
    kaV(1:N_cell) = k_area(1:N_cell) ; A0V(1:N_cell) = A0(1:N_cell) ; Aval(1:N_cell) = A(1:N_cell)
    ksV(1:N_spr)  = k_spr(1:N_spr)   ; l0V(1:N_spr)  = l0(1:N_spr)  ; lval(1:N_spr)  = l(1:N_spr)
    
  end subroutine store_props_bfr_chng
  
  subroutine apply_force_and_adjust_anomalies_bfr_implmnting_meetingAlgs(pullingOrpushing,hmuchCgX,caseNo,ExpNoCase,FrmNoCase)
    implicit none
    integer, intent(in)    :: pullingOrpushing
    real*8 , intent(in)    :: hmuchCgX
    integer, intent(in)    :: caseNo
    integer, intent(in)    :: ExpNoCase
    integer, intent(inout) :: FrmNoCase
    
    real*8  :: distance
    integer :: incrORdcrs
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
    FrmNoCase = FrmNoCase+1
    
    call pulling_or_pushing(pullingOrpushing,hmuchCgX)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNoCase,FrmNoCase) 
    FrmNoCase = FrmNoCase+1
    
    call find_distance_frm_pulley(distance)
    
    if (distance .le. 0.05d0) then
       write(*,*) "This decreasing MIDLINE APCL tension block is happening for CaseNo =",caseNo
       write(*,*) (FrmNoCase-1),"FrmNo bfr calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
       incrORdcrs = 2
       call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNoCase,FrmNoCase)
       write(*,*) (FrmNoCase-1),"FrmNo after calling manplt_apcl_membrne_tnsn with incrORdcrs=2"
    endif
    
    write(*,*) (FrmNoCase-1),"WE ARE CHECKING DIAG TNSN (BEFORE)"
    call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNoCase,FrmNoCase)
    write(*,*) (FrmNoCase-1),"WE HAVE CHECKED DIAG TNSN (AFTER)"
    
    write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIGHTNESS (BEFORE)"
    call make_diagonal_spr_more_strght(ExpNoCase,FrmNoCase)
    write(*,*) (FrmNoCase-1),"WE ARE CHECKING STRAIHTNESS (AFTER)"
    
  end subroutine apply_force_and_adjust_anomalies_bfr_implmnting_meetingAlgs
  
  
  subroutine increasingA0s_of_invgnAndinvgntd_Cells(PresVal,PresToBeMatched)
    implicit none
    real*8 , intent(in) :: PresVal(1:N_cell),PresToBeMatched
    
    real*8  :: hmuch,hmuchVal(1:2)
    integer :: i,j,cellNo
    
    call get_A0chnge_speed_for_bothInvgnInvgntdCellType(PresVal,PresToBeMatched,hmuchVal)
    
    do i = 1,2
       
       do j = 1,2
          
          if (i==1) cellNo = invgnting_cells(j)
          if (i==2) cellNo = invgnated_cells(j)
          
          A0(cellNo) = (hmuchVal(i))*A0(cellNo)
          
       enddo
       
    enddo
    
  end subroutine increasingA0s_of_invgnAndinvgntd_Cells
  
  subroutine get_A0chnge_speed_for_bothInvgnInvgntdCellType(PresVal,PresToBeMatched,hmuchVal)
    implicit none
    real*8 , intent(in)  :: PresVal(1:N_cell),PresToBeMatched
    real*8 , intent(out) :: hmuchVal(1:2)
    
    integer :: i,j,cellNo
    real*8  :: fastforwardVal,firstforwardPresMatched
    
    fastforwardVal          = 0.80d0
    firstforwardPresMatched = (fastforwardVal)*(PresToBeMatched) 
    
    do i = 1,2
   
       if (i==1) cellNo = invgnting_cells(1)
       if (i==2) cellNo = invgnated_cells(1)
       
       if (PresVal(cellNo) .le. firstforwardPresMatched) then
          hmuchVal(i) = 1.03d0
       elseif (PresVal(cellNo) .gt. firstforwardPresMatched) then
          hmuchVal(i) = 1.01d0
       endif
       
    enddo
    
  end subroutine get_A0chnge_speed_for_bothInvgnInvgntdCellType
  
  subroutine decreasingA0s_of_lastthreePair_Cells(PresVal,PresToBeMatched)
    implicit none
    real*8 , intent(in) :: PresVal(1:N_cell),PresToBeMatched
    
    real*8  :: hmuch,hmuchVal(1:3)
    integer :: countIndvd(1:3)
    integer :: i,j,cellNo
    logical :: lgcl_matched
    
    call find_if_PresVal_matchedDecr(PresVal,PresToBeMatched,lgcl_matched,countIndvd)
    call get_A0chnge_speed_for_lastThreeCellType(PresVal,PresToBeMatched,countIndvd,hmuchVal)
    
    do i = 1,3
       
       do j = 1,2
          
          if (i==1) cellNo = thrdlastCells(j)
          if (i==2) cellNo = scndlastCells(j)
          if (i==3) cellNo = frstlastCells(j)
          
          A0(cellNo) = (hmuchVal(i))*A0(cellNo)
          
       enddo
       
    enddo
    
  end subroutine decreasingA0s_of_lastthreePair_Cells
  
  subroutine get_A0chnge_speed_for_lastThreeCellType(PresVal,PresToBeMatched,countIndvd,hmuchVal)
    implicit none
    real*8 , intent(in)  :: PresVal(1:N_cell),PresToBeMatched
    integer, intent(in)  :: countIndvd(1:3)
    real*8 , intent(out) :: hmuchVal(1:3)
    
    integer :: i,j,cellNo
    real*8  :: fastforwardVal,firstforwardPresMatched
    
    fastforwardVal          = 1.15d0
    firstforwardPresMatched = (fastforwardVal)*(PresToBeMatched) 
    
    do i = 1,3
       
       if (i==1) cellNo = thrdlastCells(1)
       if (i==2) cellNo = scndlastCells(1)
       if (i==3) cellNo = frstlastCells(1)
       
       if (PresVal(cellNo) .gt. firstforwardPresMatched) then
          
          if (countIndvd(i) == 0) hmuchVal(i) = 0.97d0
          if (countIndvd(i) == 1) hmuchVal(i) = 1.00d0
          
       elseif (PresVal(cellNo) .le. firstforwardPresMatched) then
          
          if (countIndvd(i) == 0) hmuchVal(i) = 0.99d0
          if (countIndvd(i) == 1) hmuchVal(i) = 1.00d0
          
       endif
       
    enddo
    
  end subroutine get_A0chnge_speed_for_lastThreeCellType
  
  
  subroutine find_if_PresVal_matched(PresVal,PresToBeMatched,lgcl_matched)
    implicit none
    real*8, intent(in)  :: PresVal(1:N_cell),PresToBeMatched
    logical,intent(out) :: lgcl_matched
    
    integer :: count,cellNo,i,j
    
    count = 0
    
    do i = 1,2
       
       do j = 1,2
          
          if (i==1) cellNo = invgnting_cells(j)
          if (i==2) cellNo = invgnated_cells(j)
          
          if (PresVal(cellNo) .le. PresToBeMatched) then
             continue
          elseif (PresVal(cellNo) .gt. PresToBeMatched) then
             count = count+1
          endif
          
       enddo
       
    enddo
    
    lgcl_matched = .False.
    
    if (count==4) then
       lgcl_matched = .True.
    elseif (count .lt.4) then
       lgcl_matched = .False.
    else
       write(*,*) "count value not okay",count
       stop
    endif
    
    write(*,*) count,"countVal"
    
  end subroutine find_if_PresVal_matched
  
  subroutine find_if_PresVal_matchedDecr(PresVal,PresToBeMatched,lgcl_matched,countIndvd)
    implicit none
    real*8, intent(in)  :: PresVal(1:N_cell),PresToBeMatched
    logical,intent(out) :: lgcl_matched
    integer,intent(out) :: countIndvd(1:3)
    
    integer :: count,cellNo,i,j
    
    count = 0
    countIndvd(1:3) = 0
    
    do i = 1,3
       
       do j = 1,2
          
          if (i==1) cellNo = thrdlastCells(j)
          if (i==2) cellNo = scndlastCells(j)
          if (i==3) cellNo = frstlastCells(j)
          
          if (PresVal(cellNo) .gt. PresToBeMatched) then
             continue
          elseif (PresVal(cellNo) .le. PresToBeMatched) then
             if (j==1) countIndvd(i) = 1
             count = count+1
          endif
          
       enddo
       
    enddo
    
    lgcl_matched = .False.
    
    if (count==6) then
       lgcl_matched = .True.
    elseif (count .lt.6) then
       lgcl_matched = .False.
    else
       write(*,*) "count value not okay",count
       stop
    endif
    
    write(*,*) count,"countVal in Decr"
    
  end subroutine find_if_PresVal_matchedDecr
  
  subroutine find_a_more_fine_tuning_of_Pressure(ManpltTech,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ManpltTech
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: PresVal(1:N_cell)
    real*8  :: TnsnVal(1:N_spr)
    integer :: BufferCell
    real*8  :: BufferPres,DcsnPres
    integer :: PresManplt(1:Hlf_Ncell),ManpltdSpr(1:N_spr)
    real*8  :: l0_bfrMan(1:N_spr),A0_bfrMan(1:N_cell)
    
    integer :: countPresMan
    integer :: countLoop
    integer :: whatToFix
    
    DcsnPres = 0.16d0 ! This is a rather observational number
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    !call find_the_buffring_cell(PresVal,BufferCell,BufferPres,DcsnPres)
    
    BufferPres = 0.16932325055561295 !0.16228661918519233d0
    BufferCell = 4
    countLoop  = 0
    
    l0_bfrMan(1:N_spr)  = l0(1:N_spr) ; A0_bfrMan(1:N_cell) = A0(1:N_cell)
    
    do 
       !l0_bfrMan(1:N_spr)  = l0(1:N_spr) ; A0_bfrMan(1:N_cell) = A0(1:N_cell)
       
       call find_Pressure_and_Tension(PresVal,TnsnVal)
       call print_pressTnsn(countLoop,PresVal,TnsnVal)
       
       if (ManpltTech==1) then
          call get_the_PresManArry(BufferCell,BufferPres,PresVal,PresManplt)
       elseif (ManpltTech==2) then
          whatToFix=5
          call get_the_PresManArry_SegmntedWay(BufferCell,whatToFix,BufferPres,PresVal,PresManplt)
       endif
       
       call get_the_ManpltdSprArry(PresManplt,ManpltdSpr)
       call get_count_of_ZeroPresManplt(PresManplt,countPresMan)
       if (countPresMan==(Hlf_Ncell)) exit
       
       if (ManpltTech==1) then
          call manipulatePresWith_l0_manipltOnly(PresManplt,ManpltdSpr,countLoop,l0_bfrMan,A0_bfrMan,ExpNo,FrmNo)
       elseif (ManpltTech==2) then
          call manipulatePresWith_A0andl0_maniplt(PresManplt,ManpltdSpr,countLoop,BufferPres,l0_bfrMan,A0_bfrMan,ExpNo,FrmNo)
       endif
       
       countLoop = countLoop+1
       !if (countLoop==21) exit
       
    enddo
    
  end subroutine find_a_more_fine_tuning_of_Pressure
  
  
  subroutine find_the_buffring_cell(PresVal,BufferCell,BufferPres,DcsnPres)
    implicit none
    real*8,  intent(in)  :: PresVal(1:N_cell)
    integer, intent(out) :: BufferCell
    real*8 , intent(out) :: BufferPres
    real*8 , intent(in)  :: DcsnPres
    
    integer :: i
    
    do i = 1,Hlf_Ncell
       if (PresVal(i) .le. DcsnPres) then
          BufferCell = (i-1) ; write(*,*) BufferCell,"BufferCell"
          BufferPres = PresVal(i-1)
          exit
       else
          continue
       endif
    enddo
    
  end subroutine find_the_buffring_cell
  
  subroutine get_the_PresManArry(BufferCell,BufferPres,PresVal,PresManplt)
    implicit none
    integer, intent(in)  :: BufferCell
    real*8 , intent(in)  :: BufferPres
    real*8 , intent(in)  :: PresVal(1:N_cell)
    integer, intent(out) :: PresManplt(1:Hlf_Ncell)
    
    integer :: i
    real*8  :: TOL
    real*8  :: postvBufferPres,negtvBufferPres
    
    PresManplt(1:Hlf_Ncell) = 0 ; TOL = 0.01d0
    
    postvBufferPres = BufferPres + (TOL*BufferPres)
    negtvBufferPres = BufferPres - (TOL*BufferPres)
    
    write(*,*) postvBufferPres,negtvBufferPres,"postv-negtv Buffer Pres"
    
    do i = 1,Hlf_Ncell
       
       if (i==BufferCell) then
          continue
       else
          
          if (PresVal(i).le.BufferPres) then
             
             if (PresVal(i) .le. negtvBufferPres) PresManplt(i) = +1
             if (PresVal(i) .gt. negtvBufferPres) PresManplt(i) = 0
             
          elseif (PresVal(i).gt.BufferPres) then
             
             if (PresVal(i) .gt. postvBufferPres) PresManplt(i) = -1
             if (PresVal(i) .le. postvBufferPres) PresManplt(i) = 0
             
          endif
          
       endif
       
       write(*,*) PresManplt(i),i,"PresManplt"
       
    enddo
    
  end subroutine get_the_PresManArry
  
  
  subroutine get_the_PresManArry_SegmntedWay(BufferCell,whatToFix,BufferPres,PresVal,PresManplt)
    implicit none
    integer, intent(in)  :: BufferCell
    integer, intent(in)  :: whatToFix
    real*8 , intent(in)  :: BufferPres
    real*8 , intent(in)  :: PresVal(1:N_cell)
    integer, intent(out) :: PresManplt(1:Hlf_Ncell)
    
    integer :: i
    real*8  :: TOL
    real*8  :: postvBufferPres,negtvBufferPres
    integer :: CellsOfChoice(1:Hlf_Ncell)
    
    PresManplt(1:Hlf_Ncell) = 0 ; TOL = 0.01d0
    
    postvBufferPres = BufferPres + (TOL*BufferPres)
    negtvBufferPres = BufferPres - (TOL*BufferPres)
    
    write(*,*) postvBufferPres,negtvBufferPres,"postv-negtv Buffer Pres in Segmntd"
    
    call get_the_cells_to_fix(whatToFix,CellsOfChoice)
    
    do i = 1,Hlf_Ncell
       
       if (i==BufferCell) then
          continue
       else
          
          if (CellsOfChoice(i) == 1) then
             
             if (PresVal(i).le.BufferPres) then
                
                if (PresVal(i) .le. negtvBufferPres) PresManplt(i) = +1
                if (PresVal(i) .gt. negtvBufferPres) PresManplt(i) = 0
                
             elseif (PresVal(i).gt.BufferPres) then
                
                if (PresVal(i) .gt. postvBufferPres) PresManplt(i) = -1
                if (PresVal(i) .le. postvBufferPres) PresManplt(i) = 0
                
             endif
             
          elseif (CellsOfChoice(i) == 0) then
             continue
          endif
          
       endif
       
       write(*,*) PresManplt(i),i,"PresManplt"
       
    enddo
    
  end subroutine get_the_PresManArry_SegmntedWay
  
  subroutine get_the_cells_to_fix(whatToFix,CellsOfChoice)
    implicit none
    integer, intent(in)  :: whatToFix
    integer, intent(out) :: CellsOfChoice(1:Hlf_Ncell)
    
    CellsOfChoice(1:Hlf_Ncell)    = 0
    
    if (whatToFix==1) then
       CellsofChoice(Hlf_Ncell)   = 1
    elseif (whatToFix==2) then
       CellsofChoice(Hlf_Ncell-2) = 1
    elseif (whatToFix==3) then
       CellsofChoice(Hlf_Ncell-3) = 1
    elseif (whatToFix==4) then
       CellsofChoice(Hlf_Ncell-4) = 1
    elseif (whatToFix==5) then
       CellsofChoice(Hlf_Ncell-1) = 1
    endif
    
  end subroutine get_the_cells_to_fix
  
  
  subroutine find_prcntChng(countLoop,l0_bfrMan,A0_bfrMan)
    implicit none
    integer, intent(in) :: countLoop
    real*8,  intent(in) :: l0_bfrMan(1:N_spr)
    real*8,  intent(in) :: A0_bfrMan(1:N_cell)
    
    integer :: i,j,jmax
    real*8  :: l0_chngPrcnt(1:N_spr),A0_chngPrcnt(1:N_cell)
    
    open(unit=67,file='A0l0_prcntChng.dat',position='append')
    
    write(67,*) countLoop,"cloop"
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             A0_chngPrcnt(j) = (A0_bfrMan(j)-A0(j))/(A0_bfrMan(j))
             write(67,*) A0_chngPrcnt(j),A0_bfrMan(j),A0(j),j,"A0_prcnt"
          elseif (i==2) then
             l0_chngPrcnt(j) = (l0_bfrMan(j)-l0(j))/(l0_bfrMan(j))
             write(67,*) l0_chngPrcnt(j),l0_bfrMan(j),l0(j),j,"l0_prcnt"
          endif
          
          write(67,*) " "
          
       enddo
       
    enddo
       
    close(67)
    
  end subroutine find_prcntChng
  
  subroutine get_the_ManpltdSprArry(PresManplt,ManpltdSpr)
    implicit none
    integer, intent(in)  :: PresManplt(1:N_cell)
    integer, intent(out) :: ManpltdSpr(1:N_spr)
    
    integer :: i,j,jmax
    integer :: cellL,cellR
    integer :: sprL,sprR
    
    open(unit=65,file='ManpltdSprArray.dat',position='append')
    
    ManpltdSpr(1:N_spr) = 0
    
    do i = 1,Hlf_Ncell
       
       if (PresManplt(i) .ne. 0) then
          continue
       elseif (PresManplt(i) == 0) then
          
          cellL = i ; cellR = i+Hlf_Ncell
          jmax  = area_spr(cellL,0)
          
          do j = 1,jmax
             sprL = area_spr(cellL,j) ; sprR = area_spr(cellR,j)
             ManpltdSpr(sprL) = 1 ; ManpltdSpr(sprR) = 1 
          enddo
          
       endif
       
       write(65,*) ManpltdSpr(i),i,"ManpltdSpr"
       
    enddo

    write(65,*) " "
    close(65)
    
  end subroutine get_the_ManpltdSprArry
  
  subroutine manipulatePresWith_l0_manipltOnly(PresManplt,ManpltdSpr,countLoop,l0_bfrMan,A0_bfrMan,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: PresManplt(1:Hlf_Ncell)
    integer, intent(inout) :: ManpltdSpr(1:N_spr)
    integer, intent(in)    :: countLoop
    real*8 , intent(in)    :: l0_bfrMan(1:N_spr),A0_bfrMan(1:N_cell)
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,j,jmax
    integer :: cellL,cellR
    integer :: sprL,sprR
    real*8  :: hmuch
    
    do i = 1,Hlf_Ncell
       
       if (PresManplt(i) == -1) hmuch = +0.03d0
       if (PresManplt(i) == +1) hmuch = -0.03d0
       if (PresManplt(i) == 0)  hmuch = 0.000d0
       
       cellL = i ; cellR = (i+Hlf_Ncell)
       jmax  = area_spr(cellL,0) 
       
       do j = 1,jmax
          sprL = area_spr(cellL,j) ; sprR = area_spr(cellR,j)
             
          if ((ManpltdSpr(sprL)==0) .and. (ManpltdSpr(sprR)==0)) then
             l0(sprL)         = (1.00d0+hmuch)*l0(sprL) ; l0(sprR) = l0(sprL)
             ManpltdSpr(sprL) = 1                       ; ManpltdSpr(sprR) = 1
             
          elseif ((ManpltdSpr(sprL)==1) .and. (ManpltdSpr(sprR)==1)) then
             continue
          endif
             
       enddo
       
    enddo
    
    call find_prcntChng(countLoop,l0_bfrMan,A0_bfrMan)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
  end subroutine manipulatePresWith_l0_manipltOnly
  
  
  subroutine manipulatePresWith_A0andl0_maniplt(PresManplt,ManpltdSpr,countLoop,BufferPres,l0_bfrMan,A0_bfrMan,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: PresManplt(1:Hlf_Ncell)
    integer, intent(inout) :: ManpltdSpr(1:N_spr)
    integer, intent(in)    :: countLoop
    real*8 , intent(in)    :: BufferPres
    real*8 , intent(in)    :: l0_bfrMan(1:N_spr),A0_bfrMan(1:N_cell)
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: i,j,jmax
    integer :: cellL,cellR
    integer :: sprL,sprR
    real*8  :: hmuchL,hmuchA
    integer :: ManpltdSprStr(1:N_spr)
    real*8  :: AreaBfr(1:N_cell),AreaAft(1:N_cell),Area_changes(1:N_cell)
    real*8  :: TOL_prcnt
    real*8  :: AreaChngPrcnt,AreaToMatch
    integer :: PresManpltVal,prvPresMan
    integer :: cellNoL,ManPltSwitch
    real*8  :: A0_inLoop(1:N_cell),l0_inLoop(1:N_spr)
    
    ManpltdSprStr(1:N_spr) = ManpltdSpr(1:N_spr)
    AreaBfr(1:N_cell)      = A(1:N_cell)
    TOL_prcnt              = 1.00d0
    
    A0_inLoop(1:N_cell) = A0(1:N_cell)
    l0_inLoop(1:N_spr)  = l0(1:N_spr)
    
    
    ! do i = 1,Hlf_Ncell
       
    !    if (PresManplt(i) == -1) then
    !       hmuchL = +0.02d0 ; hmuchA = -0.04d0
    !    elseif (PresManplt(i) == +1) then
    !       hmuchL = -0.02d0 ; hmuchA = +0.04d0
    !    elseif (PresManplt(i) == 0)  then
    !       hmuchL = 0.000d0 ; hmuchA = 0.000d0
    !    endif
       
    !    cellL     = i                         ; cellR     = (i+Hlf_Ncell)
    !    A0(cellL) = (1.00d0+hmuchA)*A0(cellL) ; A0(cellR) = A0(cellL)
       
    !    jmax  = area_spr(cellL,0) 
       
    !    do j = 1,jmax
    !       sprL = area_spr(cellL,j) ; sprR = area_spr(cellR,j)
             
    !       if ((ManpltdSpr(sprL)==0) .and. (ManpltdSpr(sprR)==0)) then
    !          l0(sprL)         = (1.00d0+hmuchL)*l0(sprL) ; l0(sprR)     = l0(sprL)
    !          ManpltdSpr(sprL) = 1                        ; ManpltdSpr(sprR) = 1
             
    !       elseif ((ManpltdSpr(sprL)==1) .and. (ManpltdSpr(sprR)==1)) then
    !          continue
    !       endif
             
    !    enddo   
    ! enddo
    
    call get_the_singlPair_cells_frm_PresManplt(PresManplt,cellL,cellR)
    PresManpltVal = PresManplt(cellL)
    call get_the_hmuchVal(PresManpltVal,hmuchL,hmuchA)
    call with_hmucVal_toChngeA0l0(cellL,cellR,hmuchL,hmuchA)
    
    call find_prcntChng(countLoop,l0_bfrMan,A0_bfrMan)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    AreaAft(1:N_cell) = A(1:N_cell)
    call find_and_print_changes(AreaAft,AreaBfr,Area_changes)
    
    write(*,*) cellL,BufferPres,PresManpltVal,"bfr calling find_if_Pres 1"
    call find_if_PresManpltSwitches(cellL,BufferPres,PresManpltVal,ManpltSwitch)
    
    if (ManpltSwitch == +1) then
       
       do 
          l0(1:N_spr) = l0_inLoop(1:N_spr) ; A0(1:N_cell) = A0_inLoop(1:N_cell)
       
          hmuchL = hmuchL/2.0d0 ; hmuchA = hmuchA/2.0d0
          call with_hmucVal_toChngeA0l0(cellL,cellR,hmuchL,hmuchA)
          
          call find_prcntChng(countLoop,l0_bfrMan,A0_bfrMan)
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
          FrmNo = FrmNo+1
          
          prvPresMan = PresManplt(cellL)
          write(*,*) cellL,BufferPres,PresManpltVal,"bfr calling find_if_Pres 2"
          call find_if_PresManpltSwitches(cellL,BufferPres,prvPresMan,ManpltSwitch)
          
          if (ManpltSwitch == +1) then
             continue
             
          elseif (ManpltSwitch == -1) then
             AreaAft(1:N_cell) = A(1:N_cell)
             call find_and_print_changes(AreaAft,AreaBfr,Area_changes)
             exit
             
          elseif (ManpltSwitch ==  0) then
             write(*,*) "Exiting with MS=0, BEFORE match_area loop, Frame is =",(FrmNo-1)
             exit
          endif
          
       enddo
       
       if (ManpltSwitch == -1) then
          
          AreaChngPrcnt = Area_changes(cellL)*100.00d0
          AreaToMatch   = AreaBfr(cellL)
          
          write(*,*) cellL,AreaToMatch,prvPresMan,"match area 1"
          call match_area_manplt_l0(AreaChngPrcnt,cellL,AreaToMatch,prvPresMan,BufferPres,TOL_prcnt,ExpNo,FrmNo)
          
       endif
       
    elseif (ManpltSwitch == -1) then
       
       AreaChngPrcnt = Area_changes(cellL)*100.00d0
       AreaToMatch   = AreaBfr(cellL)

       write(*,*) cellL,AreaToMatch,PresManpltVal,"match area 2"
       call match_area_manplt_l0(AreaChngPrcnt,cellL,AreaToMatch,PresManpltVal,BufferPres,TOL_prcnt,ExpNo,FrmNo)
       
    elseif (ManpltSwitch == 0) then
       continue
    endif
    
    ! do i = 1,Hlf_Ncell
       
    !    if (PresManplt(i) .ne. 0) then
          
    !       AreaChngPrcnt = Area_changes(i)*100.00d0
    !       AreaToMatch   = AreaBfr(i)
    !       cellNoL       = i
          
    !       if(abs(AreaChngPrcnt) .gt. TOL_prcnt) then
    !          call match_area_manplt_l0(AreaChngPrcnt,cellNoL,AreaToMatch,TOL_prcnt,ExpNo,FrmNo)
    !       endif
          
    !    elseif (PresManplt(i) == 0) then
    !       continue
    !    endif
       
    ! enddo

    
    
  end subroutine manipulatePresWith_A0andl0_maniplt
  
  
  subroutine match_area_manplt_l0(AreaChngPrcnt,cellNoL,AreaToMatch,prvPresMan,BufferPres,TOL_prcnt,ExpNo,FrmNo)
    implicit none
    real*8, intent(inout) :: AreaChngPrcnt
    integer,intent(in)    :: cellNoL
    real*8, intent(in)    :: AreaToMatch
    integer,intent(in)    :: prvPresMan
    real*8, intent(in)    :: BufferPres
    real*8 ,intent(in)    :: TOL_prcnt
    integer,intent(in)    :: ExpNo
    integer,intent(inout) :: FrmNo
    
    real*8  :: ZERO,hmuchL
    integer :: cellNoR,sprL,sprR
    integer :: i,imax
    integer :: ManpltSwitch
    
    ZERO = 0.00d0
    
    write(*,*) AreaChngPrcnt,cellNoL,AreaToMatch,"AreaChngPrcnt At the Beginning"
    
    if (AreaChngPrcnt .gt. ZERO) hmuchL = -0.01d0
    if (AreaChngPrcnt .le. ZERO) hmuchL = +0.01d0
    
    cellNoR = cellNoL + Hlf_Ncell
    imax    = area_spr(cellNoL,0)

    do 
       
       do i = 1,imax
          sprL     = area_spr(cellNoL,i)      ; sprR = area_spr(cellNoR,i)
          l0(sprL) = (1.00d0+hmuchL)*l0(sprL) ; l0(sprR) = l0(sprL)
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       AreaChngPrcnt = ((A(cellNoL)-AreaToMatch)/AreaToMatch)*100.0d0
       write(*,*) AreaChngPrcnt,"AreaChngePrcnt inside loop"
       
       if (abs(AreaChngPrcnt).le.TOL_prcnt) then
          write(*,*) "Exiting with less area actually matched, Frame No =",(FrmNo-1)
          exit
       endif
       
       write(*,*) cellNoL,BufferPres,prvPresMan,"bfr calling find_if_Pres 3"
       call find_if_PresManpltSwitches(cellNoL,BufferPres,prvPresMan,ManpltSwitch)
       
       if (ManpltSwitch == +1) then
          write(*,*) "Exiting with MS=+1, insd match_area loop, Frame is =",(FrmNo-1)
          exit
       elseif (ManpltSwitch == -1) then
          continue
       elseif (ManpltSwitch ==  0) then
          write(*,*) "Exiting with MS=0, insd match_area loop, Frame is =",(FrmNo-1)
          exit
       endif
          
    enddo
    
  end subroutine match_area_manplt_l0
  
  


  
  subroutine get_count_of_ZeroPresManplt(PresManplt,countPresMan)
    implicit none
    integer, intent(in)  :: PresManplt(1:Hlf_Ncell)
    integer, intent(out) :: countPresMan
    
    integer :: i,j
    
    countPresMan = 0
    
    do i = 1,Hlf_Ncell
       
       if (PresManplt(i) == 0) then
          countPresMan = countPresMan+1
       else
          continue
       endif
       
    enddo
    
    write(*,*) PresManplt(1:Hlf_Ncell),"All PresManplt"
    write(*,*) countPresMan,"countPresMan"
    
  end subroutine get_count_of_ZeroPresManplt
  
  subroutine get_the_singlPair_cells_frm_PresManplt(PresManplt,cellL,cellR)
    implicit none
    integer, intent(in)  :: PresManplt(1:Hlf_Ncell)
    integer, intent(out) :: cellL,cellR
    integer :: i
    
    do i = 1,Hlf_Ncell
       if (PresManplt(i).ne.0) then
          cellL = i ; cellR = cellL+Hlf_Ncell
          exit
       else
          continue
       endif
    enddo
    
  end subroutine get_the_singlPair_cells_frm_PresManplt
  
  subroutine get_the_hmuchVal(PresManpltVal,hmuchL,hmuchA)
    implicit none
    integer, intent(in)  :: PresManpltVal
    real*8 , intent(out) :: hmuchL,hmuchA
    
    if (PresManpltVal == -1) then
       hmuchL = + 0.02d0 ; hmuchA = -0.04d0
    elseif (PresManpltVal == +1) then
       hmuchL = -0.02d0 ; hmuchA = +0.04d0
    elseif (PresManpltVal == 0) then
       hmuchL = 0.000d0 ; hmuchA = 0.000d0
    endif
    
  end subroutine get_the_hmuchVal
  
  subroutine with_hmucVal_toChngeA0l0(cellL,cellR,hmuchL,hmuchA)
    implicit none
    integer, intent(in) :: cellL,cellR
    real*8 , intent(in) :: hmuchL,hmuchA
    
    integer :: i,imax
    integer :: nsprsIncell
    integer :: sprL,sprR
    
    nsprsInCell = area_spr(cellL,0)
    A0(cellL)   = (1.00d0+hmuchA)*A0(cellL) ; A0(cellR) = A0(cellL)
    
    do i = 1,nsprsInCell
       sprL = area_spr(cellL,i) ; sprR = area_spr(cellR,i)
       l0(sprL) = (1.00d0+hmuchL)*l0(sprL) ; l0(sprR) = l0(sprL)
    enddo
    
  end subroutine with_hmucVal_toChngeA0l0
  
  subroutine find_if_PresManpltSwitches(cellL,BufferPres,prvPresMan,ManpltSwitch)
    implicit none
    integer, intent(in)  :: cellL
    real*8 , intent(in)  :: BufferPres
    integer, intent(in)  :: prvPresMan
    integer, intent(out) :: ManPltSwitch
    
    real*8  :: TOL
    real*8  :: postvBufferPres,negtvBufferPres
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    integer :: curPresMan
    
    TOL = 0.01d0
    postvBufferPres = BufferPres + (TOL*BufferPres)
    negtvBufferPres = BufferPres - (TOL*BufferPres)
    ManpltSwitch    = 10 ! unrealistic value
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    
    if (PresVal(cellL).le.BufferPres) then
       
       if (PresVal(cellL) .le. negtvBufferPres) curPresMan = +1
       if (PresVal(cellL) .gt. negtvBufferPres) curPresMan = 0
       
    elseif (PresVal(cellL).gt.BufferPres) then
       
       if (PresVal(cellL) .gt. postvBufferPres) curPresMan = -1
       if (PresVal(cellL) .le. postvBufferPres) curPresMan = 0
       
    endif

    write(*,*) PresVal(cellL),BufferPres,"PresVal+BufferPres"
    
    ManPltSwitch = -(curPresMan*prvPresMan)
    write(*,*) curPresMan,prvPresMan,ManPltSwitch,"ManPltSwitch"
    
  end subroutine find_if_PresManpltSwitches
  
  subroutine find_and_print_changes(AreaAft,AreaBfr,Area_changes)
    implicit none
    real*8, intent(in)  :: AreaAft(1:N_cell)
    real*8, intent(in)  :: AreaBfr(1:N_cell)
    real*8, intent(out) :: Area_changes(1:N_cell)
    
    integer :: i
    
    open(unit=76,file='Area_Changes.dat',position='append')
    
    do i = 1,N_cell
       Area_changes(i) = (AreaAft(i)-AreaBfr(i))/(AreaBfr(i))
       write(76,*) AreaAft(i),AreaBfr(i),Area_changes(i),i,"Area_changes"
    enddo
    
    close(76)
    
  end subroutine find_and_print_changes
  
  subroutine reduce_l0_oneLtrlMembrane_and_manpltApclBsalToMaintainShape(cellTopL,cellBotL)
    implicit none
    integer, intent(in) :: cellTopL,cellBotL
    
    integer :: cmnLtrlL(1:(NAEC_Ltrl+1)),cmnLtrlR(1:(NAEC_Ltrl+1))
    integer :: topApclL(1:(NAEC_Apcl+1)),topApclR(1:(NAEC_Apcl+1))
    integer :: botApclL(1:(NAEC_Apcl+1)),botApclR(1:(NAEC_Apcl+1))
    integer :: topBsalL(1:(NAEC_Bsal+1)),topBsalR(1:(NAEC_Bsal+1))
    integer :: botBsalL(1:(NAEC_Bsal+1)),botBsalR(1:(NAEC_Bsal+1))
    
    integer :: nsprsInACell
    integer :: i,j,jmax
    
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    
    do i = 1,3
       
       if (i==1) jmax = NAEC_Ltrl+1
       if (i==2) jmax = NAEC_Apcl+1
       if (i==3) jmax = NAEC_Bsal+1
       
       do j = 1,jmax
          
          if (i==1) then
             
             cmnLtrlL(j) = (cellTopL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
             cmnLtrlR(j) = cmnLtrlL(j) + (Hlf_Ncell)*(nsprsInACell)
             
             write(*,*) cmnLtrlL(j),cmnLtrlR(j),j,"cmnLtrl"
             
          elseif (i==2) then
             
             topApclL(j) = (cellTopL-1)*(nsprsInACell) + j
             botApclL(j) = (cellBotL-1)*(nsprsInACell) + j
             
             topApclR(j) = topApclL(j) + (Hlf_Ncell)*(nsprsInACell)
             botApclR(j) = botApclL(j) + (Hlf_Ncell)*(nsprsInACell)
             
             write(*,*) topApclL(j),topApclR(j),j,"topApcl L,R"
             write(*,*) botApclL(j),botApclR(j),j,"topApcl L,R"
             
          elseif (i==3) then
             topBsalL(j) = (cellTopL-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
             botBsalL(j) = (cellBotL-1)*(nsprsInACell) + (NAEC_Apcl+1) + j
             
             topBsalR(j) = topBsalL(j) + (Hlf_Ncell)*(nsprsInACell)
             botBsalR(j) = botBsalL(j) + (Hlf_Ncell)*(nsprsInACell)
             
             write(*,*) topBsalL(j),topBsalR(j),j,"topApcl L,R"
             write(*,*) botBsalL(j),botBsalR(j),j,"topApcl L,R"
             
          endif
          
       enddo
       
    enddo
    
    
  end subroutine reduce_l0_oneLtrlMembrane_and_manpltApclBsalToMaintainShape
  
  subroutine get_CgVal_ExpNo_and_FrmNo(pullORpush,caseNo,hmuchCgX,ExpNoCase,FrmNoCase)
    implicit none
    integer, intent(in)  :: pullORpush,caseNo
    real*8,  intent(out) :: hmuchCgX
    integer, intent(out) :: ExpNoCase,FrmNoCase
    
    if (caseNo==1) then
       if (pullORpush==1) hmuchCgX   = +0.15d0
       if (pullORpush==2) hmuchCgX   = -0.15d0
       ExpNoCase = 33 ;   FrmNoCase = 1
    elseif (caseNo==2) then
       if (pullORpush==1) hmuchCgX   = +0.20d0
       if (pullORpush==2) hmuchCgX   = -0.20d0
       ExpNoCase = 34 ;   FrmNoCase = 1
    elseif (caseNo==3) then
       if (pullORpush==1) hmuchCgX   = +0.25d0
       if (pullORpush==2) hmuchCgX   = -0.25d0
       ExpNoCase = 35 ;   FrmNoCase = 1
    elseif (caseNo==4) then
       if (pullORpush==1) hmuchCgX   = +0.30d0
       if (pullORpush==2) hmuchCgX   = -0.30d0
       ExpNoCase = 36 ;   FrmNoCase = 1
    endif
    
  end subroutine get_CgVal_ExpNo_and_FrmNo
  
  
  subroutine save_props_aft_PresIncr(flnmGen,A0V,kaV,l0V,ksV,CgXV,CgYV)
    implicit none
    integer,intent(in) :: flnmGen
    real*8, intent(in) :: A0V(1:N_cell),kaV(1:N_cell)
    real*8, intent(in) :: l0V(1:N_spr) ,ksV(1:N_spr)
    real*8, intent(in) :: CgXV(1:N_node),CgYV(1:N_node)
    
    integer :: i,j,jmax 
    
    if (flnmGen==1) open(unit=178,file='propsAftPresIncr.dat')
    if (flnmGen==2) open(unit=178,file='propsAftPresIncrAndAreaRed.dat')
    if (flnmGen==3) open(unit=178,file='propsAftPresIncrAndAreaRedAndksIncr.dat')
    if (flnmGen==4) open(unit=178,file='propsflGen4.dat')
    if (flnmGen==5) open(unit=178,file='propsflGen5.dat')
    if (flnmGen==6) open(unit=178,file='propsflGen6.dat')
    if (flnmGen==7) open(unit=178,file='propsflGen7.dat')
    if (flnmGen==8) open(unit=178,file='propsflGen8.dat')
    if (flnmGen==9) open(unit=178,file='propsflGen9.dat')
    
    if (flnmGen==10) open(unit=178,file='propsflGen10.dat')
    if (flnmGen==11) open(unit=178,file='propsflGen11.dat')
    if (flnmGen==12) open(unit=178,file='propsflGen12.dat')
    if (flnmGen==13) open(unit=178,file='propsflGen13.dat')
    
    if (flnmGen==14) open(unit=178,file='propsflGen14.dat')
    if (flnmGen==15) open(unit=178,file='propsflGen15.dat')
    if (flnmGen==16) open(unit=178,file='propsflGen16.dat')
    if (flnmGen==17) open(unit=178,file='propsflGen17.dat')
    
    if (flnmGen==18) open(unit=178,file='propsflGen18.dat')
    if (flnmGen==19) open(unit=178,file='propsflGen19.dat')
    if (flnmGen==20) open(unit=178,file='propsflGen20.dat')
    if (flnmGen==21) open(unit=178,file='propsflGen21.dat')
    if (flnmGen==22) open(unit=178,file='propsflGen22.dat')
    if (flnmGen==23) open(unit=178,file='propsflGen23.dat')
    if (flnmGen==24) open(unit=178,file='propsflGen24.dat')
    
    if (flnmGen==25) open(unit=178,file='propsflGen25.dat')
    if (flnmGen==26) open(unit=178,file='propsflGen26.dat')
    if (flnmGen==27) open(unit=178,file='propsflGen27.dat')
    if (flnmGen==28) open(unit=178,file='propsflGen28.dat')
    if (flnmGen==29) open(unit=178,file='propsflGen29.dat')

    if (flnmGen==30) open(unit=178,file='propsflGen30.dat')
    if (flnmGen==31) open(unit=178,file='propsflGen31.dat')
    if (flnmGen==32) open(unit=178,file='propsflGen32.dat')
    
    
    do i = 1,3
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) write(178,*) ksV(j),l0V(j),j
          if (i==2) write(178,*) kaV(j),A0V(j),j 
          if (i==3) write(178,*) CgXV(j),CgYV(j),j
       enddo
       
       write(178,*) " "
       
    enddo
    
    close(178)
    
  end subroutine save_props_aft_PresIncr
  
  subroutine read_props_aft_PresIncr(flnmGen,A0V,kaV,l0V,ksV,CgXV,CgYV)
    implicit none
    integer,intent(in)  :: flnmGen
    real*8, intent(out) :: A0V(1:N_cell),kaV(1:N_cell)
    real*8, intent(out) :: l0V(1:N_spr) ,ksV(1:N_spr)
    real*8, intent(out) :: CgXV(1:N_node),CgYV(1:N_node)
    
    integer :: i,j,jmax 
    
    if (flnmGen==1) open(unit=178,file='propsAftPresIncr.dat')
    if (flnmGen==2) open(unit=178,file='propsAftPresIncrAndAreaRed.dat')
    if (flnmGen==3) open(unit=178,file='propsAftPresIncrAndAreaRedAndksIncr.dat')
    if (flnmGen==4) open(unit=178,file='propsflGen4.dat')
    if (flnmGen==5) open(unit=178,file='propsflGen5.dat')
    if (flnmGen==6) open(unit=178,file='propsflGen6.dat')
    if (flnmGen==7) open(unit=178,file='propsflGen7.dat')
    if (flnmGen==8) open(unit=178,file='propsflGen8.dat')
    if (flnmGen==9) open(unit=178,file='propsflGen9.dat')
    
    if (flnmGen==10) open(unit=178,file='propsflGen10.dat')
    if (flnmGen==11) open(unit=178,file='propsflGen11.dat')
    if (flnmGen==12) open(unit=178,file='propsflGen12.dat')
    if (flnmGen==13) open(unit=178,file='propsflGen13.dat')
    
    if (flnmGen==14) open(unit=178,file='propsflGen14.dat')
    if (flnmGen==15) open(unit=178,file='propsflGen15.dat')
    if (flnmGen==16) open(unit=178,file='propsflGen16.dat')
    if (flnmGen==17) open(unit=178,file='propsflGen17.dat')
    
    if (flnmGen==18) open(unit=178,file='propsflGen18.dat')
    if (flnmGen==19) open(unit=178,file='propsflGen19.dat')
    if (flnmGen==20) open(unit=178,file='propsflGen20.dat')
    if (flnmGen==21) open(unit=178,file='propsflGen21.dat')
    if (flnmGen==22) open(unit=178,file='propsflGen22.dat')
    if (flnmGen==23) open(unit=178,file='propsflGen23.dat')
    if (flnmGen==24) open(unit=178,file='propsflGen24.dat')
    
    if (flnmGen==25) open(unit=178,file='propsflGen25.dat')
    if (flnmGen==26) open(unit=178,file='propsflGen26.dat')
    if (flnmGen==27) open(unit=178,file='propsflGen27.dat')
    if (flnmGen==28) open(unit=178,file='propsflGen28.dat')
    if (flnmGen==29) open(unit=178,file='propsflGen29.dat')
    
    if (flnmGen==30) open(unit=178,file='propsflGen30.dat')
    if (flnmGen==31) open(unit=178,file='propsflGen31.dat')
    if (flnmGen==32) open(unit=178,file='propsflGen32.dat')
    
    
    do i = 1,3
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          if (i==1) read(178,*) ksV(j),l0V(j)
          if (i==2) read(178,*) kaV(j),A0V(j) 
          if (i==3) read(178,*) CgXV(j)
       enddo
       
    enddo
    
    close(178)
    
    CgYV(1:N_Node) = 0.0d0
    
  end subroutine read_props_aft_PresIncr
  
  
  subroutine get_l0chnge_speed_for_bothCellType(AvalsBeMatched,hmuchVal)
    implicit none
    real*8 , intent(in)  :: AvalsBeMatched(1:2)
    real*8 , intent(out) :: hmuchVal(1:2)
    
    integer :: i,j,cellNo
    real*8  :: fastforwardVal
    real*8  :: firstforwardAvalMatched(1:2)
    
    
    fastforwardVal             = 0.80d0
    
    firstforwardAvalMatched(1) = (fastforwardVal)*AvalsBeMatched(1)
    firstforwardAvalMatched(2) = (fastforwardVal)*AvalsBeMatched(2)
    
    do i = 1,2
       
       if (i==1) cellNo = invgnting_cells(1)
       if (i==2) cellNo = invgnated_cells(1)
       
       if (A(cellNo) .le. firstforwardAvalMatched(i)) then
          hmuchVal(i) = 0.95d0
       elseif (A(cellNo) .gt. firstforwardAvalMatched(i)) then
          hmuchVal(i) = 0.98d0
       endif
       
    enddo
    
  end subroutine get_l0chnge_speed_for_bothCellType
  
  subroutine get_presValMatched(PresVal,cellChoice,PresToBeMatched)
    implicit none
    real*8,  intent(in)  :: PresVal(1:N_cell)
    integer, intent(in)  :: cellChoice
    real*8,  intent(out) :: PresToBeMatched
    
    real*8 :: AvgPresOfNxtInvgnatedCells
    real*8 :: tolrncePress
    
    if (cellChoice == 1) then
       AvgPresOfNxtInvgnatedCells = 0.50d0*(PresVal(nxt_invgnated_cells(1)) + PresVal(nxt_invgnated_cells(2)))
       tolrncePress    = 0.99d0 
       PresToBeMatched = (tolrncePress)*(AvgPresOfNxtInvgnatedCells)
       
    elseif (cellChoice == 2) then
       AvgPresOfNxtInvgnatedCells = 0.50d0*(PresVal(nxt_invgnated_cells(1)) + PresVal(nxt_invgnated_cells(2)))
       tolrncePress    = 0.75d0 
       PresToBeMatched = (tolrncePress)*(AvgPresOfNxtInvgnatedCells)
    endif
    
    write(*,*) PresToBeMatched,"PresToBeMatched"
    
  end subroutine get_presValMatched
  
  subroutine get_AvalsMatched(AvalsBeMatched)
    implicit none
    real*8,  intent(out) :: AvalsBeMatched(1:2)
    
    integer :: i,cell1,cell2
    real*8  :: AvalsAvrg(1:2)
    real*8  :: tolrnceArea
    
    tolrnceArea = 1.00d0
    
    do i = 1,2
       
       if (i==1) then
          cell1 = invgnated_cells(1) ; cell2 = invgnated_cells(2)  
       elseif (i==2) then
          cell1 = invgnated_cells(1) ; cell2 = invgnated_cells(2)
       endif
       
       AvalsAvrg(i)      = (0.50d0)*(A(cell1)+A(cell2))
       AvalsBeMatched(i) = (tolrnceArea)*(AvalsAvrg(i))
       
       write(*,*) AvalsBeMatched(i),i,A(cell1),A(cell2),"Avals"
    enddo
    
  end subroutine get_AvalsMatched
  
  
  subroutine find_sprsList_InvgntingInvgnatedCells(nsprs_BTC,sprsLft,sprsRgt)
    implicit none
    integer, intent(in)  :: nsprs_BTC
    integer, intent(out) :: sprsLft(1:nsprs_BTC),sprsRgt(1:nsprs_BTC)
    
    integer :: i,j,jmax,count,cellNo
    integer :: nsprsInACell
    
    count = 1
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    do i = 1,2
       
       if (i==1) cellNo = invgnting_cells(1)
       if (i==2) cellNo = invgnated_cells(2)
       
       jmax = area_spr(cellNo,0)
       
       do j = 1,jmax
          
          if (i==1) then
             
             sprsLft(count) = area_spr(cellNo,j)
             sprsRgt(count) = sprsLft(count) + (nsprsInACell*Hlf_Ncell)
             count = count+1
             
          elseif (i==2) then
             
             if (j.le.(NAEC_Ltrl+1)) then
                continue
             elseif (j.gt.(NAEC_Ltrl+1)) then
                sprsLft(count) = area_spr(cellNo,j)
                sprsRgt(count) = sprsLft(count) + (nsprsInACell*Hlf_Ncell)
                count = count+1
             endif
             
          endif
          
          write(*,*)i,j,area_spr(cellNo,j),nsprsInACell,Hlf_Ncell,"cellNo"
          write(*,*)i,j,(count-1),sprsLft(count-1),sprsRgt(count-1),"sprsLft,sprsRgt"
          
       enddo
       
    enddo
    
  end subroutine find_sprsList_InvgntingInvgnatedCells
  
  subroutine find_sprsList_BottmCells(nsprs_ETC,sprsBtmLft,sprsBtmRgt)
    implicit none
    integer, intent(in)  :: nsprs_ETC
    integer, intent(out) :: sprsBtmLft(1:nsprs_ETC),sprsBtmRgt(1:nsprs_ETC)
    
    integer :: i,j,jmax,count,cellNo
    integer :: nsprsInACell
    
    count = 1
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    do i = 1,3
       
       if (i==1) cellNo = thrdlastCells(1)
       if (i==2) cellNo = scndlastCells(1)
       if (i==3) cellNo = frstlastCells(1)
       
       jmax = area_spr(cellNo,0)
       
       do j = 1,jmax
          
          if (i==1) then
             sprsBtmLft(count) = area_spr(cellNo,j)
             sprsBtmRgt(count) = sprsBtmLft(count) + (nsprsInACell*Hlf_Ncell)
             count = count+1
             
          elseif (i.ne.1) then
             
             if (j.le.(NAEC_Ltrl+1)) then
                continue
             elseif (j.gt.(NAEC_Ltrl+1)) then
                sprsBtmLft(count) = area_spr(cellNo,j)
                sprsBtmRgt(count) = sprsBtmLft(count) + (nsprsInACell*Hlf_Ncell)
                count = count+1
             endif
             
          endif
          
          write(*,*)i,j,area_spr(cellNo,j),nsprsInACell,Hlf_Ncell,"cellNo"
          write(*,*)i,j,(count-1),sprsBtmLft(count-1),sprsBtmRgt(count-1),"sprsLft,sprsRgt"
          
       enddo
       
    enddo
    
    
  end subroutine find_sprsList_BottmCells
  
  subroutine RedRlx_l0s_selected_sprs(RedOrRlx,sprsL,sprsR,sprsLD,sprsRD,nsprs_val)
    implicit none
    integer, intent(in) :: RedOrRlx
    integer, intent(in) :: sprsL(1:nsprs_val),sprsR(1:nsprs_val)
    integer, intent(in) :: sprsLD(1:nsprs_val),sprsRD(1:nsprs_val)
    integer, intent(in) :: nsprs_val
    
    integer :: i
    real*8  :: hmuch
    integer :: sprNmL,sprNmR
    
    if (RedOrRlx == 1) hmuch = 0.97d0
    if (RedOrRlx == 2) hmuch = 1.03d0
    
    do i = 1,nsprs_val
       
       write(*,*) sprsLD(i),sprsRD(i),i,"Dcsn"
       if (sprsLD(i) .ne. sprsRD(i)) then
          write(*,*) "Lft And Rgt Decsion must be the same"
          stop
       endif
       
       sprNmL = sprsL(i) ; sprNmR = sprsR(i)
       
       if (sprsLD(i)==1) then
          l0(sprNmL) = (hmuch)*l0(sprNmL)
          l0(sprNmR) = l0(sprNmL)
       elseif (sprsLD(i) == 0) then
          continue
       endif
       
       write(*,*) l0(sprNmL),l0(sprNmR),sprNmL,sprNmR,i,"l0s lft rgt"
       
    enddo
    
  end subroutine RedRlx_l0s_selected_sprs !RedRlx=ReduceOrRelax
  
  subroutine find_if_l0sNeededToBeReduced(AvalsMatched,sprsLftD,sprsRgtD,nsprs_BTC,lgcl_Aval)
    implicit none
    real*8,  intent(in)  :: AvalsMatched(1:2)
    integer, intent(in)  :: nsprs_BTC
    integer, intent(out) :: sprsLftD(1:nsprs_BTC),sprsRgtD(1:nsprs_BTC)
    logical, intent(out) :: lgcl_Aval
    
    integer :: i,cellNo,strtNum,fnshNum
    integer :: otcm1,otcm2
    integer :: nsprsInACell
    
    lgcl_Aval = .False.
    nsprsInACell = (NAEC_Apcl+1+NAEC_Bsal+1+NAEC_Ltrl+1)
    
    do i = 1,2
       
       if (i==1) cellNo = invgnting_cells(1)
       if (i==2) cellNo = invgnated_cells(1)
       
       if (A(cellNo) .gt. AvalsMatched(i)) then
          
          if (i==1) otcm1 = 1
          if (i==2) otcm2 = 1
          
       elseif (A(cellNo) .le. AvalsMatched(i)) then
          
          if (i==1) otcm1 = 0
          if (i==2) otcm2 = 0
          
       endif
       
    enddo
    
    
    if ((otcm1==1).and.(otcm2==1)) then
       strtNum = 1
       fnshNum = nsprs_BTC
    elseif ((otcm1==1).and.(otcm2==0)) then
       strtNum = 1
       fnshNum = area_spr(cellNo,0)
    elseif ((otcm1==0).and.(otcm2==1)) then
       strtNum = area_spr(cellNo-1,0) - (NAEC_Ltrl+1) + 1
       fnshNum = nsprs_BTC
    elseif ((otcm1==0).and.(otcm2==0)) then
       strtNum = 0
    endif
    
    sprsLftD(1:nsprs_BTC) = 0
    sprsRgtD(1:nsprs_BTC) = 0
    
    
    if (strtNum.ne.0) then
       sprsLftD(strtNum:fnshNum) = 1
       sprsRgtD(strtNum:fnshNum) = 1

       sprsLftD((nsprsInACell+1):(nsprsInACell+NAEC_Ltrl+1)) = 0
       sprsRgtD((nsprsInACell+1):(nsprsInACell+NAEC_Ltrl+1)) = 0
       
    elseif (strtNum == 0) then
       write(*,*) "strtNum is ZERO"
       lgcl_Aval = .True.
    endif
    
    
  end subroutine find_if_l0sNeededToBeReduced
  
  
  subroutine reducing_l0s_incrsing_ks_apcl_frm_invgnted_Cells(hmuchEach,ExpNo,FrmNo)
    implicit none
    real*8, intent(in)    :: hmuchEach
    integer,intent(in)    :: ExpNo
    integer,intent(inout) :: FrmNo
    
    real*8  :: l0_prv(1:N_spr),ks_prv(1:N_spr)
    integer :: strting_cellL,fnshing_cellL
    integer :: strting_cellR,fnshing_cellR
    integer :: numCells
    integer :: count,count_max,i,j
    integer :: nsprsInACell
    integer :: lftCell,rgtCell,apclLft,apclRgt
    real*8  :: step
    real*8  :: distnce
    
    strting_cellL = invgnated_cells(1) ; strting_cellR = invgnated_cells(2)
    fnshing_cellL = Hlf_Ncell-1        ; fnshing_cellR = fnshing_cellL+Hlf_Ncell
    
    numCells = fnshing_cellL-strting_cellL+1
    write(*,*) numCells,"numCells"
    
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    write(*,*) nsprsInACell,"nsprsInACell"
    
    l0_prv(1:N_spr) = l0(1:N_spr)
    ks_prv(1:N_spr) = k_spr(1:N_spr)
    
    count_max = 10
    
    do count = 1,count_max
       
       step = (count-1)*(hmuchEach)
       write(*,*) step,count,"step"
       
       do i = 1,numCells 
          
          lftCell = strting_cellL+i-1 ; rgtCell = strting_cellR+i-1
          apclLft = (lftCell-1)*(nsprsInACell)+1 ; apclRgt = (rgtCell-1)*(nsprsInACell)+1
          write(*,*) apclLft,apclRgt,lftCell,rgtCell,"apcls"
          
          l0(apclLft) = l0_prv(apclLft)*(1.00d0-step) ; l0(apclRgt) = l0(apclLft)
          
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       if (count==5) then !!! VERY NONGENERAL
          write(*,*) (FrmNo-1),"FrmNo Bfr Insd"
          call make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNo,FrmNo)
          write(*,*) (FrmNo-1),"FrmNo Aft Insd"
       endif
          
       call find_distance_frm_pulley(distnce)
       if (distnce.le.0.05d0) then
          write(*,*) "MET IN L0 SHORTEN"
       endif
       
    enddo
    
    call find_distance_frm_pulley(distnce)
    write(*,*) (FrmNo-1),"FrmNo Bfr ks incr"
    
    if (distnce .le. 0.05d0) then
       continue
    else
       count_max = 3
       
       do count = 1,count_max
          
          step = (count-1)*(hmuchEach*10.00)
          write(*,*) step,count,"step"
          
          do i = 1,numCells
             
             lftCell = strting_cellL+i-1 ; rgtCell = strting_cellR+i-1
             apclLft = (lftCell-1)*(nsprsInACell)+1 ; apclRgt = (rgtCell-1)*(nsprsInACell)+1
             write(*,*) apclLft,apclRgt,lftCell,rgtCell,"apcls"
             
             k_spr(apclLft) = ks_prv(apclLft)*(1.00d0+step) ; k_spr(apclRgt) = k_spr(apclLft)
             
          enddo
          
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
          FrmNo = FrmNo+1
          
          call find_distance_frm_pulley(distnce)
          if (distnce .le. 0.05d0) then
             write(*,*) "MET in ks incr"
             exit
          endif
       enddo
       
    endif
    
  end subroutine reducing_l0s_incrsing_ks_apcl_frm_invgnted_Cells
  
  subroutine make_diagonal_spr_tnsn_closeToZeroAndNegtv(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: diagsprL(1:(NAEC_Ltrl+1)),diagsprR(1:(NAEC_Ltrl+1))
    real*8  :: l0V(1:N_spr),l0_Str(1:N_spr)
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    real*8  :: TolrnceTnsn,ZERO
    integer :: sprL,sprR
    real*8  :: TolrncTnsn,sumTnsn,avrgTnsn
    integer :: i,j,jmax
    integer :: cellL,cellR,nsprsInACell
    real*8  :: hmuch
    
    l0V(1:N_spr) = l0(1:N_spr)
    
    cellL = invgnting_cells(1) ; cellR = invgnting_cells(2)
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,(NAEC_Ltrl+1)
       diagsprL(i) = (cellL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       diagsprR(i) = (cellR-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       
       write(*,*) diagsprL(i),diagsprR(i),i,"diagspr"
    enddo
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    TolrncTnsn = 0.02d0 ; ZERO = 0.00d0
    call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
    
    if (avrgTnsn .le. ZERO) then
       
       if (abs(avrgTnsn) .gt. TolrncTnsn) then
          hmuch = 0.02d0
          
          do 
             
             l0_Str(1:N_spr) = l0(1:N_spr)
             
             do i = 1,(NAEC_Ltrl+1)
                sprL = diagsprL(i) ; sprR = diagsprR(i)
                l0(sprL) = (1.00d0+hmuch)*l0(sprL) 
                l0(sprR) = l0(sprL)
             enddo
             
             call Equilibrate_system
             call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
             FrmNo = FrmNo+1
             
             call find_Pressure_and_Tension(PresVal,TnsnVal)
             call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
             
             if (avrgTnsn .gt. ZERO) then
                l0(1:N_spr) = l0_Str(1:N_spr)
                hmuch = (hmuch/2.0d0)
             elseif (avrgTnsn .le. ZERO) then
                
                if (abs(avrgTnsn).le.TolrncTnsn) then
                   write(*,*) "FIRST NEGATIVE TENSION, EDIT AND TOLRNCE LIMIT"
                   exit
                elseif (abs(avrgTnsn).gt.TolrncTnsn) then
                   continue
                endif
                
             endif
                
          enddo
          
       elseif (abs(avrgTnsn) .le. TolrncTnsn) then
          write(*,*) "TOLRNCE LIMIT,NO LOOP NEEDED"
       endif
       
    elseif (avrgTnsn .gt. ZERO) then
       
       hmuch = 0.02d0
       
       do 
       
          l0_Str(1:N_spr) = l0(1:N_spr)
          
          do i = 1,(NAEC_Ltrl+1)
             sprL = diagsprL(i) ; sprR = diagsprR(i)
             l0(sprL) = (1.00d0-hmuch)*l0(sprL) 
             l0(sprR) = l0(sprL)
          enddo
          
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
          FrmNo = FrmNo+1
          
          call find_Pressure_and_Tension(PresVal,TnsnVal)
          call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
          
          if (avrgTnsn .le. ZERO) then
             
             if (abs(avrgTnsn) .gt. TolrncTnsn) then
                l0(1:N_spr) = l0_Str(1:N_spr)
                hmuch = (hmuch/2.0d0)
                
             elseif (abs(avrgTnsn) .le. TolrncTnsn) then
                write(*,*) "FIRST POSITIVE, THEN NEGATIVE IN TOLRNCE LIMIT"
                exit
             endif
             
          elseif (avrgTnsn .gt. ZERO) then
             continue
          endif
          
       enddo
       
    endif
    
  end subroutine make_diagonal_spr_tnsn_closeToZeroAndNegtv
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine make_diagonal_spr_tnsn_reduced_toFindICS
    implicit none
    integer :: diagsprL(1:(NAEC_Ltrl+1)),diagsprR(1:(NAEC_Ltrl+1))
    real*8  :: l0V(1:N_spr),l0_Str(1:N_spr)
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    real*8  :: TolrnceTnsn,ZERO
    integer :: sprL,sprR
    real*8  :: TolrncTnsn,sumTnsn,avrgTnsn
    integer :: i,j,jmax
    integer :: cellL,cellR,nsprsInACell
    real*8  :: hmuch,tol_Rdfn
    logical :: lgcl_Rdfn
    
    l0V(1:N_spr) = l0(1:N_spr)
    
    call find_the_invagniating_and_invaginated_cell
    cellL = invgnting_cells(1) ; cellR = invgnting_cells(2)
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,(NAEC_Ltrl+1)
       diagsprL(i) = (cellL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       diagsprR(i) = (cellR-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       write(*,*) diagsprL(i),diagsprR(i),i,"diagspr"
    enddo
    
    call find_Pressure_and_Tension(PresVal,TnsnVal)
    TolrncTnsn = 0.02d0 ; ZERO = 0.00d0
    call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
    
    if (avrgTnsn .le. ZERO) then
       
       if (abs(avrgTnsn) .gt. TolrncTnsn) then
          hmuch = 0.02d0
          
          do 
             
             l0_Str(1:N_spr) = l0(1:N_spr)
             
             do i = 1,(NAEC_Ltrl+1)
                sprL = diagsprL(i)                 ; sprR = diagsprR(i)
                l0(sprL) = (1.00d0+hmuch)*l0(sprL) ; l0(sprR) = l0(sprL)
             enddo
             
             call Equilibrate_only_NI_model
             
             call find_Pressure_and_Tension(PresVal,TnsnVal)
             call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
             
             if (avrgTnsn .gt. ZERO) then
                l0(1:N_spr) = l0_Str(1:N_spr)
                hmuch = (hmuch/2.0d0)
             elseif (avrgTnsn .le. ZERO) then
                if (abs(avrgTnsn).le.TolrncTnsn) then
                   write(*,*) "starts with NEGATIVE TENSION, aft EDIT,get inside TOLRNCE LIMIT"
                   exit
                elseif (abs(avrgTnsn).gt.TolrncTnsn) then
                   continue
                endif
             endif
             
             lgcl_Rdfn=.False. ; tol_Rdfn=0.07d0
             call get_decision_of_redefining_diffWay(tol_Rdfn,lgcl_Rdfn)
             if (lgcl_Rdfn .eqv. .True.) then
                write(*,*) "Exit with lgcl_Rdfn=True" ; exit
             endif
             
          enddo
          
       elseif (abs(avrgTnsn) .le. TolrncTnsn) then
          write(*,*) "TOLRNCE LIMIT,NO LOOP NEEDED"
       endif
       
    elseif (avrgTnsn .gt. ZERO) then
       
       hmuch = 0.02d0
       
       do 
       
          l0_Str(1:N_spr) = l0(1:N_spr)
          
          do i = 1,(NAEC_Ltrl+1)
             sprL = diagsprL(i) ; sprR = diagsprR(i)
             l0(sprL) = (1.00d0-hmuch)*l0(sprL) 
             l0(sprR) = l0(sprL)
          enddo
          
          call Equilibrate_only_NI_model
          
          call find_Pressure_and_Tension(PresVal,TnsnVal)
          call find_theAvrgDiagTnsn(TnsnVal,diagsprL,avrgTnsn)
          
          if (avrgTnsn .le. ZERO) then
             
             if (abs(avrgTnsn) .gt. TolrncTnsn) then
                l0(1:N_spr) = l0_Str(1:N_spr)
                hmuch = (hmuch/2.0d0)
                
             elseif (abs(avrgTnsn) .le. TolrncTnsn) then
                write(*,*) "starts with FIRST POSITIVE, then NEGATIVE & also IN TOLRNCE LIMIT"
                exit
             endif
             
          elseif (avrgTnsn .gt. ZERO) then
             continue
          endif
          
       enddo
       
    endif
    
  end subroutine make_diagonal_spr_tnsn_reduced_toFindICS
  
  subroutine make_diagonal_spr_more_strght(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: diagsprL(1:(NAEC_Ltrl+1)),diagsprR(1:(NAEC_Ltrl+1))
    real*8  :: ratio,incrmntPrcnt,admissibleRatio
    
    real*8  :: l0V(1:N_spr),l0_Str(1:N_spr)
    real*8  :: PresVal(1:N_cell),TnsnVal(1:N_spr)
    real*8  :: TolrnceTnsn,ZERO
    integer :: sprL,sprR
    real*8  :: TolrncTnsn,sumTnsn,avrgTnsn
    integer :: i,j,jmax
    integer :: cellL,cellR,nsprsInACell
    real*8  :: hmuch
    integer :: count
    
    l0V(1:N_spr) = l0(1:N_spr)
    
    cellL = invgnting_cells(1) ; cellR = invgnting_cells(2)
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,(NAEC_Ltrl+1)
       diagsprL(i) = (cellL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       diagsprR(i) = (cellR-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       
       write(*,*) diagsprL(i),diagsprR(i),i,"diagspr"
    enddo
    
    call find_the_ratio_btwn_strghtAndsgmntdLine(diagsprL,ratio,incrmntPrcnt)
    
    admissibleRatio = 1.50d0 
    hmuch = 0.01d0
    
    if (incrmntPrcnt .gt. admissibleRatio) then
       
       count = 1
       
       do
          l0_str(1:N_spr) = l0(1:N_spr)
          
          do i = 1,(NAEC_Ltrl+1)
             sprL = diagsprL(i) ; sprR = diagsprR(i)
             l0(sprL) = (1.00d0-hmuch)*l0(sprL)
             l0(sprR) = l0(sprL)
          enddo
         
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
          FrmNo = FrmNo+1
          count = count+1
          
          call find_the_ratio_btwn_strghtAndsgmntdLine(diagsprL,ratio,incrmntPrcnt)
          write(*,*) ratio,incrmntPrcnt,(FrmNo-1),"ratio-incrmntPrcnt-Frame"
          
          if (incrmntPrcnt .le. admissibleRatio) then
             write(*,*) FrmNo-1,"Frame No at exiting"
             exit
          endif

          write(*,*) count,"count"
          if (count==11) exit
          
       enddo
       
    elseif (ratio .le. admissibleRatio) then
       continue
    endif
    
  end subroutine make_diagonal_spr_more_strght
  
  
  subroutine manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: incrORdcrs
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: ks_stre(1:N_spr),l0_stre(1:N_spr)
    integer :: count
    real*8  :: strtPointK,loopIncrK,incrK
    real*8  :: strtPointL,loopIncrL,incrL
    real*8  :: distance
    
    integer, allocatable :: DoubleMembrns(:,:)
    integer :: num_of_DM
    
    ks_stre(1:N_spr) = k_spr(1:N_spr)
    l0_stre(1:N_spr) = l0(1:N_spr)
    
    count    = 1
    
    if (incrORdcrs == 1) then
       incrK = 1.25d0 ; incrL = 0.95d0
    elseif (incrORdcrs == 2) then
       incrK = 0.75d0 ; incrL = 1.05d0   
    endif
    
    
    call find_the_number_of_double_membrnes(num_of_DM)
    allocate(DoubleMembrns(1:num_of_DM,1:2))
    call find_the_double_membrnes(DoubleMembrns,num_of_DM)
    
    do   
       call change_kspr_and_l0_of_DMs(DoubleMembrns,num_of_DM,incrK,incrL)
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       write(*,*) (FrmNo-1),"FrmNo in loop"
       
       call find_distance_frm_pulley(distance)
       
       if (incrORdcrs == 1) then
          if (distance .le. 0.05d0) then
             write(*,*) count,"count"
             exit
          endif
       elseif (incrORdcrs == 2) then
          if (distance .gt. 0.10d0) then
             write(*,*) count,"count"
             exit
          endif
       endif
       count = count+1
    enddo
    
  end subroutine manplt_apcl_membrne_tnsn_invgnting_and_botmCells
  
  subroutine find_the_number_of_double_membrnes(num_of_DM)
    implicit none
    integer, intent(out) :: num_of_DM !number of Double membranes
    
    num_of_DM = Hlf_Ncell-(invgnated_cells(1))+1
    write(*,*) num_of_DM,"number of Double Membranes"
    
  end subroutine find_the_number_of_double_membrnes
  
  subroutine find_the_double_membrnes(DoubleMembrns,num_of_DM)
    implicit none
    integer, intent(in)  :: num_of_DM
    integer, intent(out) :: DoubleMembrns(1:num_of_DM,1:2)
    
    integer :: i,j
    integer :: nsprsInAcell,cellNum
    
    nsprsInAcell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    write(*,*) nsprsInAcell,"nsprsInAcell"

    write(*,*) "Check the Double Membranes"
    
    do i = 1,num_of_DM
       cellNum = invgnated_cells(1) + (i-1)
       
       DoubleMembrns(i,1) = (cellNum-1)*(nsprsInAcell) + 1
       DoubleMembrns(i,2) = DoubleMembrns(i,1) + (nsprsInAcell)*(Hlf_Ncell)
    enddo
    
  end subroutine find_the_double_membrnes
  
  subroutine change_kspr_and_l0_of_DMs(DoubleMembrns,num_of_DM,incrK,incrL)
    implicit none
    integer, intent(in) :: num_of_DM
    integer, intent(in) :: DoubleMembrns(1:num_of_DM,1:2)
    real*8 , intent(in) :: incrK,incrL
    
    integer :: i,j
    integer :: sprNmL,sprNmR
    
    do i = 1,num_of_DM
       
       sprNmL=DoubleMembrns(i,1) ; sprNmR=DoubleMembrns(i,2) 
       
       k_spr(sprNmL) = (incrK)*(k_spr(sprNmL)) ; k_spr(sprNmR) = k_spr(sprNmL)
       l0(sprNmL)    = (incrL)*(l0(sprNmL))    ; l0(sprNmR)    = l0(sprNmL)
       
       write(*,*) i,sprNmL,sprNmR,"i,sprNms"
       write(*,*) k_spr(sprNmL),k_spr(sprNmR),"k_spr of sprNmL,sprNmR"
       write(*,*) l0(sprNmL),l0(sprNmR),"l0 of sprNmL,sprNmR"
       write(*,*) " "
       
    enddo
    
  end subroutine change_kspr_and_l0_of_DMs
  
  subroutine find_the_ratio_btwn_strghtAndsgmntdLine(diagSpr,incrInsegmntLen,incrmntPrcnt)
    implicit none
    integer, intent(in)  :: diagSpr(1:(NAEC_Ltrl+1))
    real*8 , intent(out) :: incrInsegmntLen,incrmntPrcnt
    real*8, allocatable  :: x(:),y(:)
    
    integer :: i,j,jmax
    integer :: N_sgmnts,N_points
    real*8  :: sgmntLen,lenDiagSpr,strgtLineDis
    real*8  :: x_diffSqr,y_diffSqr
    
    integer :: nodesIntheSpr
    integer :: frstNode,scndNode
    integer :: sprNm
    
    N_sgmnts   = NAEC_Ltrl+1 ; N_points = N_sgmnts+1
    lenDiagSpr = 0.00d0
    
    allocate(x(1:N_points),y(1:N_points))
    
    do i = 1,N_sgmnts
       
       sprNm = diagSpr(i)

       nodesIntheSpr = spr_node(sprNm,0)
       if (nodesIntheSpr .ne. 2) then
          write(*,*) "nodes can't be higher than 2"
          stop
       endif
       
       frstNode = spr_node(sprNm,1)
       scndNode = spr_node(sprNm,2)
       
       if (i==1) then
          x(i)   = node_xy(frstNode,1)
          y(i)   = node_xy(frstNode,2)
          x(i+1) = node_xy(scndNode,1)
          y(i+1) = node_xy(scndNode,2)
          
       elseif (i.gt.1) then
          x(i+1) = node_xy(scndNode,1)
          y(i+1) = node_xy(scndNode,2)
       endif
       
    enddo
     
    do i = 1,N_sgmnts
       sgmntLen   = sqrt(((x(i+1)-x(i)))**2 + ((y(i+1)-y(i)))**2)
       lenDiagSpr = lenDiagSpr + sgmntLen
       
       write(*,*) lenDiagSpr,sgmntLen,i,"sgmnt length"
    enddo
    
    strgtLineDis = sqrt((x(N_points)-x(1))**2 + (y(N_points)-y(1))**2)
    write(*,*) strgtLineDis,"strgtLineDis"
    
    incrInsegmntLen = (lenDiagSpr/strgtLineDis)*100.0d0
    write(*,*) incrInsegmntLen,"ratio"

    incrmntPrcnt = incrInsegmntLen - 100.00d0
    write(*,*) incrmntPrcnt,"incrmnt Prcnt"
    
  end subroutine find_the_ratio_btwn_strghtAndsgmntdLine
  
  
  subroutine find_theAvrgDiagTnsn(TnsnVal,diagspr,avrgTnsn)
    implicit none
    real*8,  intent(in)  :: TnsnVal(1:N_spr)
    integer, intent(in)  :: diagspr(1:(NAEC_Ltrl+1))
    real*8,  intent(out) :: avrgTnsn
    
    integer :: i,sprNm
    real*8  :: sumTnsn
    
    sumTnsn = 0.00d0
    
    do i = 1,(NAEC_Ltrl+1)
       sprNm   = diagspr(i) 
       sumTnsn = sumTnsn + TnsnVal(sprNm)
    enddo
    
    avrgTnsn = sumTnsn/(NAEC_Ltrl+1)
    
  end subroutine find_theAvrgDiagTnsn
  
  subroutine reducing_l0s_of_diagonal_spr(hmuch,howmanyT,ExpNo,FrmNo)
    implicit none
    real*8,  intent(in) :: hmuch
    integer, intent(in) :: howmanyT
    integer, intent(in) :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: cellL,cellR
    integer :: diagsprL(1:(NAEC_Ltrl+1)),diagsprR(1:(NAEC_Ltrl+1))
    real*8  :: l0V(1:N_spr)
    integer :: i,j,nsprsInACell
    real*8  :: hmuchVal
    
    l0V(1:N_spr) = l0(1:N_spr)
    cellL = invgnting_cells(1) ; cellR = invgnting_cells(2)
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do i = 1,(NAEC_Ltrl+1)
       diagsprL(i) = (cellL-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       diagsprR(i) = (cellR-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       
       write(*,*) diagsprL(i),diagsprR(i),i,"diagspr"
    enddo
    
    do i = 1,howmanyT
       hmuchVal = (i-1)*hmuch

       do j = 1,(NAEC_Ltrl+1)
          l0(diagsprL(j)) = (1.00d0-hmuchVal)*l0V(diagsprL(j))
          l0(diagsprR(j)) = (1.00d0-hmuchVal)*l0V(diagsprR(j))
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
    enddo
    
  end subroutine reducing_l0s_of_diagonal_spr
  
  !!!*** ROUTINES FOR CHECKING IF PRESSURE DIFFERENCES ARE NEEDED (FINISHING POINT) ***!!!
  
  subroutine reduce_buckling_of_LtrlMmbrne(ExpNo,FrmNo,hmuch,LastCellToCmpr)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    integer, intent(in)    :: LastCellToCmpr
    real*8 , intent(in)    :: hmuch
    
    integer                :: cellL,cellR
    integer, allocatable   :: LtsprL(:),LtsprR(:)
    integer                :: sprCntEachCell
    integer                :: sprCntApcl,sprCntBsal,sprCntLtrl
    integer, allocatable   :: sP(:,:)!spring Pair
    real*8 , allocatable   :: Initl0(:,:)
    integer                :: Nstp,sprPair(1:2),i,j
    real*8                 :: hm,stepS

    if (modelID==1) then
       write(*,*) "routine: reduce_buckling_of_LtrlMmbrne, not for modelID=1"
       stop
    endif
    
    allocate(LtsprL(1:(NAEC+1)),LtsprR(1:(NAEC+1)))
    allocate(sP(1:(NAEC+1),1:2))
    allocate(Initl0(1:(NAEC+1),1:2))
    
    cellL = (LastCellToCmpr+2) ; cellR = (N_cell/2) + (LastCellToCmpr)
    
    sprCntApcl = 1      !apical is not divided
    sprCntBsal = NAEC+1
    sprCntLtrl = NAEC+1
    
    sprCntEachCell = sprCntApcl + sprCntBsal + sprCntLtrl
    
    do i = 1,(NAEC+1)
       LtsprL(i) = (sprCntEachCell*cellL) + (sprCntApcl+sprCntBsal) + i
       LtsprR(i) = (sprCntEachCell*cellR) + (sprCntApcl+sprCntBsal) + i
    enddo
    
    open(unit=16,file="reduce_buckling.dat")
    
    Nstp  = 5
    stepS = (hmuch)/(Nstp-1)
    
    
    do i = 1,(NAEC+1)
       sP(i,1)     = LtsprL(i)   ; sP(i,2)     = LtsprR(i)
       Initl0(i,1) = l0(sP(i,1)) ; Initl0(i,2) = l0(sP(i,2))
       write(16,*) sP(i,1:2),i,"sP and i"
       write(16,*) Initl0(i,1:2),i,"Initl0 and i"
    enddo
    
    do i = 1,Nstp
       
       hm = (i-1)*stepS
       
       do j = 1,(NAEC+1)
          sprPair(1:2) = sP(j,1:2)
          call shortening_sprPair(sprPair,hm)
          write(16,*) l0(sprPair(1)),l0(sprPair(2)),j,"l0 AftShrtn insd main sb"
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo + 1
       
       if (i.ne.Nstp) then
          
          do j = 1,(NAEC+1)
             l0(sP(j,1)) = Initl0(j,1)
             l0(sP(j,2)) = Initl0(j,2)
          enddo
          
       endif
       
    enddo
    
    write(16,*) ExpNo,FrmNo,"ExpNo,FrmNo"
    write(16,*) " "
    close(16)
    
    
  end subroutine reduce_buckling_of_LtrlMmbrne
  
  
  subroutine borrowPropsOf_InitCellAtBottom(strct1,strct4,modelWntd,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strct1,strct4
    integer, intent(in)    :: modelWntd
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: ICB, nsprs_ICB ! Initiator Cell at Bottom
    integer :: spr_nm
    real*8  :: E
    
    if (modelWntd==1) then
       continue
    elseif (modelWntd==2) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    endif
    
    ICB = N_cell
    nsprs_ICB = area_spr(ICB,0)
    
    if (modelID==1) then
       k_area(ICB) = ka_strctTN(strct1,ICB)
       A0(ICB)     = A0_strctTN(strct1,ICB)
       
       do i = 1,nsprs_ICB
          spr_nm = area_spr(ICB,i)
          k_spr(spr_nm) = ks_strctTN(strct1,spr_nm)
          l0(spr_nm)    = l0_strctTN(strct1,spr_nm)
       enddo
       
    elseif (modelID==2) then
       k_area(ICB) = ka_strctNI(strct1,ICB)
       A0(ICB)     = A0_strctNI(strct1,ICB)
       
       do i = 1,nsprs_ICB
          spr_nm = area_spr(ICB,i)
          k_spr(spr_nm) = ks_strctNI(strct1,spr_nm)
          l0(spr_nm)    = l0_strctNI(strct1,spr_nm)
       enddo
       
    endif
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    
    if (modelWntd==1) then
       continue
    elseif (modelWntd==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
  end subroutine borrowPropsOf_InitCellAtBottom
  
  subroutine borrowPropsOf_bottomCells(strct1,strct2,strctToRead,howManyPBC,IfBotmMostIncluded,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strct1,strct2,strctToRead
    integer, intent(in)    :: howManyPBC ! PBC=Pair of Bottom Cells
    integer, intent(in)    :: IfBotmMostIncluded
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: spr_nm
    integer :: i,j,imax
    integer :: cellL,cellR
    integer :: nsprs_cellL,nsprs_cellR
    integer :: saveWhom,cntrlSt
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
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
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    open(unit=46,file='borrowPropsOfBC.dat')
    
    cellL = 0 ; cellR = 0
    
    nsprs_cellL = 0
    nsprs_cellR = 0
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       
       if (ifBotmMostIncluded==0) then
          
          if (i.lt.(imax-(howManyPBC-1))) then ! if 3 botm then imax-2
             continue
          elseif (i.ge.(imax-(howManyPBC-1))) then
             cellL = i ; cellR = Hlf_Ncell+i
             nsprs_cellL = area_spr(cellL,0)
             nsprs_cellR = area_spr(cellR,0)
             
             write(46,*) i,cellL,cellR,nsprs_cellL,nsprs_cellR,"cells & sprs"
          endif    
          
       elseif (ifBotmMostIncluded==1) then
          
          if (i.lt.(imax-(howManyPBC-1))) then
             continue
             
          elseif (i.ge.(imax-(howManyPBC-1))) then
             
             if (i.ne.imax) then
                cellL = i ; cellR = Hlf_Ncell+i
                nsprs_cellL = area_spr(cellL,0)
                nsprs_cellR = area_spr(cellR,0)
                
             elseif (i==imax) then
                cellL = N_cell ; cellR = 0
                nsprs_cellL = area_spr(cellL,0)
                nsprs_cellR = 0
             endif
             
          endif
          
       endif
       
       if (cellL .ne. 0) then
          
          if (modelID==1) k_area(cellL) = ka_strctTN(strct1,cellL)
          if (modelID==2) k_area(cellL) = ka_strctNI(strct1,cellL)
          if (modelID==1) A0(cellL)     = A0_strctTN(strct1,cellL)
          if (modelID==2) A0(cellL)     = A0_strctNI(strct1,cellL)
          
          write(46,*) cellL,k_area(cellL),A0(cellL),"cellL and k_area,A0"
          
          do j = 1,nsprs_cellL
             
             spr_nm = area_spr(cellL,j)
             
             if (modelID==1) k_spr(spr_nm) = ks_strctTN(strct1,spr_nm)
             if (modelID==2) k_spr(spr_nm) = ks_strctNI(strct1,spr_nm)
             if (modelID==1) l0(spr_nm)    = l0_strctTN(strct1,spr_nm)
             if (modelID==2) l0(spr_nm)    = l0_strctNI(strct1,spr_nm)
             
             write(46,*) spr_nm,k_spr(spr_nm),l0(spr_nm),"sprLs' and k_spr,l0 frm L"
          enddo
          
       endif
       
       
       if (cellR .ne. 0) then
          
          if (modelID==1) k_area(cellR) = ka_strctTN(strct1,cellR)
          if (modelID==2) k_area(cellR) = ka_strctNI(strct1,cellR)
          if (modelID==1) A0(cellR)     = A0_strctTN(strct1,cellR)
          if (modelID==2) A0(cellR)     = A0_strctNI(strct1,cellR)
          
          write(46,*) cellR,k_area(cellR),A0(cellR),"cellR and k_area,A0"
          
          do j = 1,nsprs_cellR
             
             spr_nm = area_spr(cellR,j)
             
             if (modelID==1) k_spr(spr_nm) = ks_strctTN(strct1,spr_nm)
             if (modelID==2) k_spr(spr_nm) = ks_strctNI(strct1,spr_nm)
             if (modelID==1) l0(spr_nm)    = l0_strctTN(strct1,spr_nm)
             if (modelID==2) l0(spr_nm)    = l0_strctNI(strct1,spr_nm)
             
             write(46,*) spr_nm,k_spr(spr_nm),l0(spr_nm),"sprLs' and k_spr,l0 frm R"
          enddo
          
       endif
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(46,*) (FrmNo-1),"aftBorrowing Bot3 prop"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strct2,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    cntrlSt = 2
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(46)
    
  end subroutine borrowPropsOf_bottomCells
  
  
  subroutine borrowPropsOf_FirstCells(strctFrm,strctTo,howManyPFC,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctFrm,strctTo
    integer, intent(in)    :: howManyPFC ! PFC=Pair of First Cells
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: strctToRead
    real*8  :: E
    integer :: spr_nm
    integer :: i,j,m,imax
    integer :: cellL,cellR,cellNm
    integer :: nsprs_cellL,nsprs_cellR,nsprs_cell
    integer :: saveWhom
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
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
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    strctToRead = strctTo
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    open(unit=46,file='borrowPropsOfFC.dat')

    write(46,*) (FrmNo-1),"FrmNo bfrBorrowFC"
    cellL = 0 ; cellR = 0
    
    nsprs_cellL = 0
    nsprs_cellR = 0
    
    imax = howManyPFC
    
    do i = 1,imax
       
       cellL = i ; cellR = Hlf_Ncell+i
       nsprs_cellL = area_spr(cellL,0)
       nsprs_cellR = area_spr(cellR,0)
       
       write(46,*) i,cellL,cellR,nsprs_cellL,nsprs_cellR,"cells & sprs"
       
       do j = 1,2 ! 1 for left, 2 for right
          
          if (j==1) then
             cellNm = cellL
             nsprs_cell = nsprs_cellL
          elseif (j==2) then
             cellNm = cellR
             nsprs_cell = nsprs_cellR
          endif
          
          k_area(cellNm) = ka_strctNI(strctFrm,cellnm)
          A0(cellNm)     = A0_strctNI(strctFrm,cellnm)
          
          write(46,*) cellnm,k_area(cellnm),A0(cellnm),j,"cellnm and k_area,A0,j frm L/R"
          
          do m = 1,nsprs_cell
             
             spr_nm = area_spr(cellNm,m)
             
             k_spr(spr_nm) = ks_strctNI(strctFrm,spr_nm)
             l0(spr_nm)    = l0_strctNI(strctFrm,spr_nm)
             
             write(46,*) spr_nm,k_spr(spr_nm),l0(spr_nm),j,"sprLs' and k_spr,l0,j frm L/R"
          enddo
          
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(46,*) (FrmNo-1),"aftBorrowing FCell prop"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctTo,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(46)
    
  end subroutine borrowPropsOf_FirstCells
  
  subroutine copyShiftdPropsFrmCs1ToCs4_or_Cs4toCs1(copyDirectn,strctToRead,cntrlSt,CS1,CS4,TopOrBottom,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: copyDirectn
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: CS1,CS4
    integer, intent(in)    :: TopOrBottom
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: cellL,cellLToBe
    integer :: cellR,cellRToBe
    integer :: sprCntEachCell,sprCntLastCell
    integer :: sprNm,sprToBeNm
    integer :: i,j,m,istrt,imax,mmax
    integer :: saveWhom
    real*8  :: E
    integer :: LastCellAbove
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    real*8, allocatable :: A0InCS1(:),kaInCS1(:)
    real*8, allocatable :: l0InCS1(:),ksInCS1(:)
    real*8, allocatable :: CgXVInCS1(:)
    real*8, allocatable :: A0InCS4(:),kaInCS4(:)
    real*8, allocatable :: l0InCS4(:),ksInCS4(:)
    real*8, allocatable :: CgXVInCS4(:)

    integer, allocatable :: sprL(:),sprLToBe(:)
    integer, allocatable :: sprR(:),sprRToBe(:)
    
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

    open(unit=64,file='copyShiftedProps.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
       
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    
    allocate(A0InCS1(1:N_cell),kaInCS1(1:N_cell))
    allocate(l0InCS1(1:N_spr),ksInCS1(1:N_spr))
    allocate(CgXVInCS1(1:N_node))
    
    allocate(A0InCS4(1:N_cell),kaInCS4(1:N_cell))
    allocate(l0InCS4(1:N_spr),ksInCS4(1:N_spr))
    allocate(CgXVInCS4(1:N_node))
    
    call read_strctProps_withAddedCell(CS1)
    call read_strctProps_withAddedCell(CS4)
    
    A0InCS1(1:N_cell)  = A0_strctNI(CS1,1:N_cell)
    kaInCS1(1:N_cell)  = ka_strctNI(CS1,1:N_cell)
    l0InCS1(1:N_spr)   = l0_strctNI(CS1,1:N_spr)
    ksInCS1(1:N_spr)   = ks_strctNI(CS1,1:N_spr)
    CgXVInCS1(1:N_node) = CgX_strctNI(CS1,1:N_node)
    
    A0InCS4(1:N_cell)  = A0_strctNI(CS4,1:N_cell)
    kaInCS4(1:N_cell)  = ka_strctNI(CS4,1:N_cell)
    l0InCS4(1:N_spr)   = l0_strctNI(CS4,1:N_spr)
    ksInCS4(1:N_spr)   = ks_strctNI(CS4,1:N_spr)
    CgXVInCS4(1:N_node) = CgX_strctNI(CS4,1:N_node)
    
    if ((strctToRead.ne.CS1) .AND. (strctToRead.ne.CS4))  then
       write(*,*) strctToRead,CS1,CS4,"strctTRead"
       write(*,*) " strctToRead needs to equal to either CS1 or CS4"
       write(*,*) " chk in fl: nodeInstrdCntrlStates ; sb:copyShiftdPropsFrmCs1ToCs4_or_Cs4toCs1"
       stop
    endif
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0 
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(64,*) (FrmNo-1),"Frm aft reading CS1/CS4"
    
    
    if (TopOrBottom==1) then !TOTAL APPROACH
       call find_LastCellofTopLayer(LastCellAbove)
       istrt = 1
       imax  = Hlf_Ncell+1
       
    elseif (TopOrBottom==2) then !NON-TOTAL APPROACH
       write(*,*) "VERY_IMPRTNT,CHANGE THIS EVERYTIME BFR USING THIS ROUTINE"
       call sleep(5)
       
       istrt = Hlf_Ncell-1 !!VERY_IMPRTNT,CHANGE THIS EVERYTIME BFR USING THIS ROUTINE
       imax  = Hlf_Ncell-1
       
    endif
    
    if (TopOrBottom == 1) then
       write(64,*)"For TOP"
       
       do i = istrt,imax
          
          if (i.le.LastCellAbove) then
             
             if (copyDirectn == 1) then !CS1 to CS4
                
                cellL     = i       ; cellR     = Hlf_Ncell+i
                cellLToBe = cellL+1 ; cellRToBe = cellR+1
                
                write(64,*) cellL,cellLToBe,"cellL's TOP CoPYDIR 1"
                write(64,*) cellR,cellRToBe,"cellR's TOP CoPYDIR 1"
                
             elseif (copyDirectn == 2) then ! CS4 to CS1
                
                cellL     = i       ; cellR     = Hlf_Ncell+i
                cellLToBe = cellL-1 ; cellRToBe = cellR-1
                
                write(64,*) cellL,cellLToBe,"cellL's TOP CoPYDIR 2"
                write(64,*) cellR,cellRToBe,"cellR's TOP CoPYDIR 2"
                
             endif
             
          elseif ((i.gt.LastCellAbove) .and. (i.lt.imax)) then
             !cellL     = i     ; cellR     = Hlf_Ncell+i
             !cellLToBe = cellL ; cellRToBe = cellR
             
             cellL     = 0 ; cellR     = 0
             cellLToBe = 0 ; cellRToBe = 0
             
          elseif (i==imax) then
             !cellL     = N_cell ; cellR     = 0
             !cellLToBe = cellL  ; cellRToBe = 0
             
             cellL     = 0 ; cellR     = 0
             cellLToBe = 0 ; cellRToBe = 0
          endif
          
          write(64,*) cellL,cellLToBe,"cellL's"
          write(64,*) cellR,cellRToBe,"cellR's"
          
          if (i==istrt) then
             sprCntEachCell = 1 + (NAEC+1)*1 + (NAEC+1)*1 !apcl + bsal + ltrl
             allocate(sprL(1:sprCntEachCell),sprLToBe(1:sprCntEachCell))
             allocate(sprR(1:sprCntEachCell),sprRToBe(1:sprCntEachCell))
             write(64,*) sprCntEachCell,"sprCntEachCell"
             
             sprL = 0 ; sprLToBe = 0
             sprR = 0 ; sprRToBe = 0
             
          elseif (i==imax) then
             deallocate(sprL,sprLToBe)
             deallocate(sprR,sprRToBe)
             
             sprCntLastCell = (NAEC+1)*1 !ltrl
             allocate(sprL(1:sprCntLastCell),sprLToBe(1:sprCntLastCell))
             allocate(sprR(1:sprCntLastCell),sprRToBe(1:sprCntLastCell))
             write(64,*) sprCntLastCell,"sprCntLastCell"
             
             sprL = 0 ; sprLToBe = 0
             sprR = 0 ; sprRToBe = 0
          endif
          
          do j = 1,2 ! 1 for left, 2 for right
             
             if (i.ne.imax) continue
             if ((i==imax) .and. (j==2)) exit
             
             if (j==1) then
                
                if (cellL.ne.0) then
                   
                   if (copyDirectn == 1) then
                      k_area(cellL) = ka_strctNI(CS1,cellLToBe)
                      A0(cellL)     = A0_strctNI(CS1,cellLToBe)
                      
                   elseif (copyDirectn == 2) then
                      k_area(cellL) = ka_strctNI(CS4,cellLToBe)
                      A0(cellL)     = A0_strctNI(CS4,cellLToBe)
                   endif
                   
                elseif (cellL==0) then
                   continue
                endif
                
             elseif (j==2) then
                
                if (cellR.ne.0) then
                   
                   if (copyDirectn == 1) then
                      k_area(cellR) = ka_strctNI(CS1,cellRToBe)
                      A0(cellR)     = A0_strctNI(CS1,cellRToBe)

                   elseif (copyDirectn == 2) then
                      k_area(cellR) = ka_strctNI(CS4,cellRToBe)
                      A0(cellR)     = A0_strctNI(CS4,cellRToBe)
                   endif
                   
                elseif (cellR==0) then
                   continue
                endif
                
             endif
             
             if (i.ne.imax) mmax = sprCntEachCell
             if (i == imax) mmax = sprCntLastCell   
             
             do m = 1,mmax
                
                if (j==1) then
                   
                   if (cellL.ne.0) then
                      sprL(m)     = (sprCntEachCell)*(cellL-1) + m
                      sprLToBe(m) = (sprCntEachCell)*(cellLToBe-1) + m
                      
                      sprNm = sprL(m) ; sprToBeNm = sprLToBe(m)
                      write(64,*) sprNm,sprToBeNm,"spr's L"
                      
                      if (copyDirectn == 1) then
                         k_spr(sprNm) = ks_strctNI(CS1,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS1,sprToBeNm)
                         
                      elseif (copyDirectn == 2) then
                         k_spr(sprNm) = ks_strctNI(CS4,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS4,sprToBeNm)
                      endif
                      
                   elseif (cellL==0) then
                      continue
                   endif
                   
                elseif (j==2) then
                   
                   if (cellR.ne.0) then
                      sprR(m)     = (sprCntEachCell)*(cellR-1) + m
                      sprRToBe(m) = (sprCntEachCell)*(cellRToBe-1) + m
                      
                      sprNm = sprR(m) ; sprToBeNm = sprRToBe(m)
                      write(64,*) sprNm,sprToBeNm,"spr's R"
                      
                      if (copyDirectn == 1) then
                         k_spr(sprNm) = ks_strctNI(CS1,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS1,sprToBeNm)
                         
                      elseif (copyDirectn == 2) then
                         k_spr(sprNm) = ks_strctNI(CS4,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS4,sprToBeNm)
                      endif
                      
                   elseif (cellR==0) then
                      continue
                   endif
                   
                endif
                
             enddo
             
          enddo
          
       enddo
       
    elseif (TopOrBottom==2) then
       write(64,*)"For BOTTOM"
       
       do i = istrt,imax
          
          if (copyDirectn == 1) then !CS1 to CS4
             
             cellL     = i       ; cellR     = Hlf_Ncell+i
             cellLToBe = cellL+1 ; cellRToBe = cellR+1
             
             write(64,*) cellL,cellLToBe,"cellL's BOTTOM CoPYDIR 1"
             write(64,*) cellR,cellRToBe,"cellR's BOTTOM CoPYDIR 1"
             
          elseif (copyDirectn == 2) then ! CS4 to CS1
             
             cellL     = i       ; cellR     = Hlf_Ncell+i
             cellLToBe = cellL-1 ; cellRToBe = cellR-1
             
             write(64,*) cellL,cellLToBe,"cellL's BOTTOM CoPYDIR 2"
             write(64,*) cellR,cellRToBe,"cellR's BOTTOM CoPYDIR 2"
             
          endif
          
          
          if (i==istrt) then
             
             sprCntEachCell = 1 + (NAEC+1)*1 + (NAEC+1)*1 !apcl + bsal + ltrl
             allocate(sprL(1:sprCntEachCell),sprLToBe(1:sprCntEachCell))
             allocate(sprR(1:sprCntEachCell),sprRToBe(1:sprCntEachCell))
             write(64,*) sprCntEachCell,"sprCntEachCell"
             
             sprL = 0 ; sprLToBe = 0
             sprR = 0 ; sprRToBe = 0
             
          endif
          
          do j = 1,2 ! 1 for left, 2 for right
             
             if (j==1) then
                
                if (cellL.ne.0) then
                   
                   if (copyDirectn == 1) then
                      k_area(cellL) = ka_strctNI(CS1,cellLToBe)
                      A0(cellL)     = A0_strctNI(CS1,cellLToBe)
                      
                   elseif (copyDirectn == 2) then
                      k_area(cellL) = ka_strctNI(CS4,cellLToBe)
                      A0(cellL)     = A0_strctNI(CS4,cellLToBe)
                   endif
                   
                elseif (cellL==0) then
                   continue
                endif
                
             elseif (j==2) then
                
                if (cellR.ne.0) then
                   
                   if (copyDirectn == 1) then
                      k_area(cellR) = ka_strctNI(CS1,cellRToBe)
                      A0(cellR)     = A0_strctNI(CS1,cellRToBe)
                      
                   elseif (copyDirectn == 2) then
                      k_area(cellR) = ka_strctNI(CS4,cellRToBe)
                      A0(cellR)     = A0_strctNI(CS4,cellRToBe)
                   endif
                   
                elseif (cellR==0) then
                   continue
                endif
                
             endif
             
             mmax = sprCntEachCell
             
             do m = 1,mmax
                
                if (j==1) then
                   
                   if (cellL.ne.0) then
                      sprL(m)     = (sprCntEachCell)*(cellL-1) + m
                      sprLToBe(m) = (sprCntEachCell)*(cellLToBe-1) + m
                      
                      sprNm = sprL(m) ; sprToBeNm = sprLToBe(m)
                      write(64,*) sprNm,sprToBeNm,"spr's L"
                      
                      if (copyDirectn == 1) then
                         k_spr(sprNm) = ks_strctNI(CS1,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS1,sprToBeNm)
                         
                      elseif (copyDirectn == 2) then
                         k_spr(sprNm) = ks_strctNI(CS4,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS4,sprToBeNm)
                      endif
                      
                   elseif (cellL==0) then
                      continue
                   endif
                   
                elseif (j==2) then
                   
                   if (cellR.ne.0) then
                      sprR(m)     = (sprCntEachCell)*(cellR-1) + m
                      sprRToBe(m) = (sprCntEachCell)*(cellRToBe-1) + m
                      
                      sprNm = sprR(m) ; sprToBeNm = sprRToBe(m)
                      write(64,*) sprNm,sprToBeNm,"spr's R"
                      
                      if (copyDirectn == 1) then
                         k_spr(sprNm) = ks_strctNI(CS1,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS1,sprToBeNm)
                         
                      elseif (copyDirectn == 2) then
                         k_spr(sprNm) = ks_strctNI(CS4,sprToBeNm)
                         l0(sprNm)    = l0_strctNI(CS4,sprToBeNm)
                      endif
                      
                   elseif (cellR==0) then
                      continue
                   endif
                   
                endif
                
             enddo
             
          enddo
          
       enddo
       
    endif
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(64,*) (FrmNo-1),"aft copying frm CS1 prop"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(64)
    
  end subroutine copyShiftdPropsFrmCs1ToCs4_or_Cs4toCs1
  
  
  subroutine make_Invaginated_Rectangular_areas_Equivalent(strctToRead,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E,AreaOfFrstCell
    real*8  :: AreaDecr
    real*8  :: ApprxAreaRectCell
    real*8  :: firstChkPoint,TolArea
    integer :: cnt_TrueLgcls
    
    integer :: cellNm
    integer :: i,j
    real*8  :: Area_diff
    real*8  :: ZERO
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
    logical, allocatable :: lgcl_Rect(:)
    real*8 , allocatable :: Area_clsNess(:)
    real*8 , allocatable :: AreaIncrEachStep(:) 
    integer, allocatable :: countForFrstChkPoint(:)
    
    
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
    
    open(unit=64,file='Reduce_ActualAreas_ofRectangularCells.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(64,*) (FrmNo-1),"bfr chaging A0 of rectangular cells"
    
    ZERO = 0.00d0
    
    AreaOfFrstCell = A(1)
    AreaDecr       = 1.10d0
    ApprxAreaRectCell = A(1)/AreaDecr
    
    firstChkPoint = 0.10d0 ; TolArea = 0.01d0
    
    call Find_the_Rect_Cells(cntrlSt)
    
    allocate(lgcl_Rect(1:N_rectCells))
    allocate(Area_clsNess(1:N_rectCells))
    allocate(AreaIncrEachStep(1:N_rectCells))
    allocate(countForFrstChkPoint(1:N_rectCells))
    
    lgcl_Rect(1:N_rectCells)            = .False.
    AreaIncrEachStep(1:N_rectCells)     = 1.05d0 
    countForFrstChkPoint(1:N_rectCells) = 0
    cnt_TrueLgcls                       = 0
    
    do 
       
       do i = 1,N_rectCells
          
          cellNm = RectCells(i)
          
          if (lgcl_Rect(i) .eqv. .False.) then
             A0(cellNm) = A0(cellNm)*AreaIncrEachStep(i)
             
          elseif (lgcl_Rect(i) .eqv. .True.) then
             continue
          endif
          
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       
       do i = 1,N_rectCells
          
          cellNm = RectCells(i)
          
          Area_diff       = ApprxAreaRectCell-A(cellNm)
          Area_clsNess(i) = Area_diff/ApprxAreaRectCell
          
          if (Area_diff.lt.ZERO) then
             write(*,*) "Area increased too much for cellNm=",cellNm,A(cellNm),A0(cellNm)
          endif
          
          if ((Area_clsNess(i) .lt. firstChkPoint) .and. (Area_clsNess(i).gt.TolArea)) then
             
             if (countForFrstChkPoint(i) == 0) then
                AreaIncrEachStep(i) = 1.005d0
                countForFrstChkPoint(i) = 1
                
             elseif (countForFrstChkPoint(i) == 1) then
                continue
             endif
             
          elseif (Area_clsNess(i) .lt. TolArea) then
             AreaIncrEachStep(i) = 1.00d0
             lgcl_Rect(i) = .True.
          endif
          
       enddo
       
       
       do i = 1,N_rectCells
          
          if (lgcl_Rect(i) .eqv. .True.) then
             cnt_TrueLgcls = cnt_TrueLgcls+1
          elseif (lgcl_Rect(i) .eqv. .False.) then
             continue
          endif
          
       enddo
       
       if (cnt_TrueLgcls == N_rectCells) then
          exit
       elseif (cnt_TrueLgcls .lt. N_rectCells) then
          cnt_TrueLgcls = 0
       elseif (cnt_TrueLgcls .gt. N_rectCells) then
          write(*,*) "cnt_TrueLgcls cant be greater than N_rectCells"
       endif
       
    enddo
    
    close(64)
    
    write(*,*) "after increasing Areas of Rectangular Cells"
    
  end subroutine make_Invaginated_Rectangular_areas_Equivalent
  
  
  subroutine Find_the_Rect_Cells(cntrl_st)
    implicit none
    integer, intent(in)  :: cntrl_st
    
    
    if (cntrl_st==1 .or. cntrl_st==4) then
       N_rectCellPair = 1
       N_rectCells    = 2*N_rectCellPair
       
    elseif (cntrl_st==2) then
       N_rectCellPair = 1
       N_rectCells    = 2*N_rectCells
       
    !elseif (cntrl_st==4) then
     !  N_rectCellPair = 2
     !  N_rectCells    = 2*N_rectCellPair
    endif
    
    allocate(RectCells(1:N_rectCells))
    RectCells(1:N_rectCells) = 0
    
    do i = 1,N_rectCells
       
       if (mod(i,2) ==1) then
          RectCells(i) = Hlf_Ncell-1 - (i/2)
          
       elseif (mod(i,2)==0) then
          RectCells(i) = (N_cell-1) - 1 - ((i-1)/2)
       endif
       
    enddo
    
  end subroutine Find_the_Rect_Cells
  
  
  subroutine making_approaching_and_invaginating_cells_similar_area
    implicit none
    
    call Find_the_approaching_and_invaginating_cells_similar_area
    
  end subroutine making_approaching_and_invaginating_cells_similar_area
  
  
  subroutine Find_the_approaching_and_invaginating_cells_similar_area
    implicit none
    
  end subroutine Find_the_approaching_and_invaginating_cells_similar_area
  
  
  subroutine convert_oneside_tilted_trapezoid_to_recangle_with_same_area(strctToRead,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: cell_no
    real*8  :: area_Trapezoid 
    real*8  :: DBPS ! DBPS = Distance Between Perpendicular Sides
    real*8  :: APLS ! ALPS = Average Length of Perpendicular Sides 
    integer :: nsprs_InEachCell
    integer :: i,j,jmax
    integer :: N_itm,N_step
    real*8  :: sum_lngths
    real*8  :: ZERO,ONE
    real*8  :: firstChkPoint,TolLen
    
    real*8 , allocatable :: PresVal(:),TnsnVal(:)
    integer, allocatable :: Apcls(:),Bsals(:),Ltrls1(:),Ltrls2(:) 
    
    real*8 , allocatable :: chngOflenBsal(:),chngOflenLtrl1(:),chngOflenLtrl2(:)
    real*8 , allocatable :: lenChngEachStepBsal(:),lenChngEachStepLtrl1(:),lenChngEachStepLtrl2(:)
    real*8 , allocatable :: length_clsNessBsal(:),length_clsNessLtrl1(:),length_clsNessLtrl2(:)
    real*8 , allocatable :: l0stepBsal(:),l0stepLtrl1(:),l0stepLtrl2(:)
    
    logical, allocatable :: lgcl_Bsals(:),lgcl_Ltrls1(:),lgcl_Ltrls2(:)
    integer :: cnt_lgclsBsal,cnt_lgclsLtrl1,cnt_lgclsLtrl2,N_postvLgcls,N_memBsLt1Lt2
    real*8  :: avg_lngthBsal,avg_lngthLtrl,len_diff
    integer :: sprNm
    integer :: cnt_insdLoop
    
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
    
    
    open(unit=62,file='convertTrapToRectangularCells.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(62,*) (FrmNo-1),"bfr chaging A0 of rectangular cells"
    
    
    cell_no          = Hlf_Ncell-1
    nsprs_InEachCell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    N_itm            = 4
    cnt_insdLoop     = 0
    
    ZERO          = 0.00d0 ; ONE       = 1.00d0
    cnt_lgclsBsal = 0 ; cnt_lgclsLtrl1 = 0       ; cnt_lgclsLtrl2 = 0
    firstChkPoint = 0.10d0 ; TolLen    = 0.025d0
    
    N_memBsLt1Lt2 = (NAEC_Bsal+1) + (NAEC_Ltrl+1) + (NAEC_Ltrl+1)
    
    allocate(Apcls(1:(NAEC_Apcl+1)),Bsals(1:(NAEC_Bsal+1)),Ltrls1(1:(NAEC_Ltrl+1)),Ltrls2(1:(NAEC_Ltrl+1)))
    
    Apcls(1:(NAEC_Apcl+1))  = 0
    Bsals(1:(NAEC_Bsal+1))  = 0
    Ltrls1(1:(NAEC_Ltrl+1)) = 0
    Ltrls2(1:(NAEC_Ltrl+1)) = 0
    
    allocate(l0stepBsal(1:(NAEC_Bsal+1)),l0stepLtrl1(1:(NAEC_Ltrl+1)),l0stepLtrl2(1:(NAEC_Ltrl+1)))
    allocate(lgcl_Bsals(1:(NAEC_Bsal+1)),lgcl_Ltrls1(1:(NAEC_Ltrl+1)),lgcl_Ltrls2(1:(NAEC_Ltrl+1)))
    allocate(chngOflenBsal(1:(NAEC_Bsal+1)),chngOflenLtrl1(1:(NAEC_Ltrl+1)),chngOflenLtrl2(1:(NAEC_Ltrl+1)))
    allocate(lenChngEachStepBsal(1:(NAEC_Bsal+1)),lenChngEachStepLtrl1(1:(NAEC_Ltrl+1)),lenChngEachStepLtrl2(1:(NAEC_Ltrl+1)))
    allocate(length_clsNessBsal(1:(NAEC_Bsal+1)),length_clsNessLtrl1(1:(NAEC_Ltrl+1)),length_clsNessLtrl2(1:(NAEC_Ltrl+1)))
    
    lgcl_Bsals(1:(NAEC_Bsal+1))  = .False.
    lgcl_Ltrls1(1:(NAEC_Ltrl+1)) = .False.
    lgcl_Ltrls2(1:(NAEC_Ltrl+1)) = .False.
    
    call get_the_apcls_bsals_and_ltrls12(cell_no,Apcls,Bsals,Ltrls1,Ltrls2,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl)
    
    write(62,*) Apcls(1:(NAEC_Apcl+1)),"Apcls"
    write(62,*) Bsals(1:(NAEC_Bsal+1)),"Bsals"
    write(62,*) Ltrls1(1:(NAEC_Ltrl+1)),"Ltrls1"
    write(62,*) Ltrls2(1:(NAEC_Ltrl+1)),"Ltrls2"
    write(62,*) " "
    
    area_Trapezoid = A(cell_no)
    DBPS           = l(Apcls(1))
    sum_lngths     = 0.0d0
    
    write(62,*) area_Trapezoid,DBPS,"area_Trp,DBPS"
    
    do i = 1,(NAEC_Ltrl+1)
       sum_lngths = sum_lngths + l(Ltrls1(i)) + l(Ltrls2(i))   
    enddo
    
    avg_lngthBsal = DBPS/(NAEC_Bsal+1)
    avg_lngthLtrl = sum_lngths/((NAEC_Ltrl+1)+(NAEC_Ltrl+1))
    
    write(62,*) avg_lngthBsal,DBPS,"avg-DBPS Bsal"
    write(62,*) avg_lngthLtrl,sum_lngths,"avg-sum Ltrl"
    
    
    do i = 1,(N_itm-1)
       
       if (i==1) jmax = NAEC_Bsal+1
       if (i==2) jmax = NAEC_Ltrl+1
       if (i==3) jmax = NAEC_Ltrl+1
       
       do j = 1,jmax
          
          if (i==1) chngOflenBsal(j)  = avg_lngthBsal - l(Bsals(j))
          if (i==2) chngOflenLtrl1(j) = avg_lngthLtrl - l(Ltrls1(j))
          if (i==3) chngOflenLtrl2(j) = avg_lngthLtrl - l(Ltrls2(j))
          
       enddo
       
    enddo
    
    write(62,*) " "
    write(62,*) chngOflenBsal(1:(NAEC_Bsal+1)),"chngOflenBsal"
    write(62,*) chngOflenLtrl1(1:(NAEC_Ltrl+1)),"chngOflenLtrl1"
    write(62,*) chngOflenLtrl2(1:(NAEC_Ltrl+1)),"chngOflenLtrl2"
    
    
    do i = 1,(N_itm-1)
       
       if (i==1) jmax = NAEC_Bsal+1
       if (i==2) jmax = NAEC_Ltrl+1
       if (i==3) jmax = NAEC_Ltrl+1
       
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (chngOflenBsal(j) .lt. ZERO) lenChngEachStepBsal(j) = (0.10d0)*chngOflenBsal(j)!1.1 
             if (chngOflenBsal(j) .gt. ZERO) lenChngEachStepBsal(j) = (0.10d0)*chngOflenBsal(j)!0.9
             
          elseif (i==2) then
             
             if (chngOflenLtrl1(j) .lt. ZERO)lenChngEachStepLtrl1(j)=(0.10d0)*chngOflenLtrl1(j)!1.10
             if (chngOflenLtrl1(j) .gt. ZERO)lenChngEachStepLtrl1(j)=(0.10d0)*chngOflenLtrl1(j)!0.90
             
          elseif (i==3) then
             
             if (chngOflenLtrl2(j) .lt. ZERO)lenChngEachStepLtrl2(j)=(0.10d0)*chngOflenLtrl2(j)!1.10
             if (chngOflenLtrl2(j) .gt. ZERO)lenChngEachStepLtrl2(j)=(0.10d0)*chngOflenLtrl2(j)!0.90
             
          endif
          
       enddo
       
    enddo
    
    write(62,*) " "
    write(62,*) lenChngEachStepBsal(1:(NAEC_Bsal+1)),"step Bsal"
    write(62,*) lenChngEachStepLtrl1(1:(NAEC_Ltrl+1)),"step Ltrl1"
    write(62,*) lenChngEachStepLtrl2(1:(NAEC_Ltrl+1)),"step Ltrl2"
    
    write(62,*) " "
    write(62,*) "Info's before changing"
    
    write(62,*) l0(Apcls(1)),"Apcl l0"
    write(62,*) l(Apcls(1)) ,"Apcl length"
    write(62,*) " "
    write(62,*) l0(Bsals(1)) ,l0(Bsals(2)) ,l0(Bsals(3)),"Basal l0"
    write(62,*) l(Bsals(1))  ,l(Bsals(2))  ,l(Bsals(3)),"Basal length"
    write(62,*) " "
    write(62,*) l0(Ltrls1(1)),l0(Ltrls1(2)),l0(Ltrls1(3)),"Lateral1 l0"
    write(62,*) l(Ltrls1(1)) ,l(Ltrls1(2)) ,l(Ltrls1(3)) ,"Lateral1 length"
    write(62,*) " "
    write(62,*) l0(Ltrls2(1)),l0(Ltrls2(2)),l0(Ltrls2(3)),"Lateral2 l0"
    write(62,*) l(Ltrls2(1)) ,l(Ltrls2(2)) ,l(Ltrls2(3)) ,"Lateral length"
    write(62,*) " "
    
    
    do
       
       do i = 1,(N_itm-1)
          
          if (i==1) jmax = NAEC_Bsal+1
          if (i==2) jmax = NAEC_Ltrl+1
          if (i==3) jmax = NAEC_Ltrl+1
          
          do j = 1,jmax
             
             if (i==1) then
                
                sprNm     = Bsals(j)
                
                if (lgcl_Bsals(j) .eqv. .False.) l0(sprNm) = l0(sprNm) + lenChngEachStepBsal(j)
                if (lgcl_Bsals(j) .eqv. .True.)  continue
                
             elseif (i==2) then
                
                sprNm     = Ltrls1(j)
                
                if (lgcl_Ltrls1(j) .eqv. .False.) l0(sprNm) = l0(sprNm) + lenChngEachStepLtrl1(j)
                if (lgcl_Ltrls1(j) .eqv. .True.)  continue
                
             elseif (i==3) then
                
                sprNm     = Ltrls2(j)
                
                if (lgcl_Ltrls2(j) .eqv. .False.) l0(sprNm) = l0(sprNm) + lenChngEachStepLtrl2(j)
                if (lgcl_Ltrls2(j) .eqv. .True.)  continue !l0(sprNm)*lenChngEachStepLtrl2(j)
                
             endif
             
          enddo
          
       enddo
       
       call create_samel0Chnge_inRight(cell_no)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       
       write(62,*) l0(Apcls(1)),"Apcl l0"
       write(62,*) l(Apcls(1)) ,"Apcl length"
       write(62,*) " "
       write(62,*) l0(Bsals(1)) ,l0(Bsals(2)) ,l0(Bsals(3)),"Basal l0"
       write(62,*) l(Bsals(1))  ,l(Bsals(2))  ,l(Bsals(3)),"Basal length"
       write(62,*) " "
       write(62,*) l0(Ltrls1(1)),l0(Ltrls1(2)),l0(Ltrls1(3)),"Lateral1 l0"
       write(62,*) l(Ltrls1(1)) ,l(Ltrls1(2)) ,l(Ltrls1(3)) ,"Lateral1 length"
       write(62,*) " "
       write(62,*) l0(Ltrls2(1)),l0(Ltrls2(2)),l0(Ltrls2(3)),"Lateral2 l0"
       write(62,*) l(Ltrls2(1)) ,l(Ltrls2(2)) ,l(Ltrls2(3)) ,"Lateral length"
       write(62,*) " "
       
       
       do i = 1,(N_itm-1)
          
          if (i==1) jmax = NAEC_Bsal+1
          if (i==2) jmax = NAEC_Ltrl+1
          if (i==3) jmax = NAEC_Ltrl+1
          
          do j = 1,jmax
             
             if (i==1) then
                
                sprNm    = Bsals(j)
                len_diff = avg_lngthBsal - l(sprNm)
                length_clsNessBsal(j)  = abs(len_diff) / avg_lngthBsal
                
                if ((length_clsNessBsal(j).lt.firstChkPoint) .and. (length_clsNessBsal(j).gt.TolLen)) then
                   
                   if (lenChngEachStepBsal(j).gt.ONE)lenChngEachStepBsal(j) = (0.10d0)*chngOflenBsal(j)
                   if (lenChngEachStepBsal(j).lt.ONE)lenChngEachStepBsal(j) = (0.10d0)*chngOflenBsal(j)
                   
                elseif (length_clsNessBsal(j) .lt. TolLen) then
                   
                   if (cnt_lgclsBsal .ne. (NAEC_Bsal+1)) then
                      lgcl_Bsals(j) = .True.
                      cnt_lgclsBsal = cnt_lgclsBsal+1
                   elseif (cnt_lgclsBsal == (NAEC_Bsal+1)) then
                      continue
                   endif
                   
                endif
                
             elseif (i==2) then
                
                sprNm    = Ltrls1(j)
                len_diff = avg_lngthLtrl - l(sprNm)
                length_clsNessLtrl1(j) = abs(len_diff) / avg_lngthLtrl
                
                if ((length_clsNessLtrl1(j).lt.firstChkPoint) .and. (length_clsNessLtrl1(j).gt.TolLen)) then
                   
                   if(lenChngEachStepLtrl1(j).gt.ONE) lenChngEachStepLtrl1(j)=(0.10d0)*chngOflenLtrl1(j)
                   if(lenChngEachStepLtrl1(j).lt.ONE) lenChngEachStepLtrl1(j)=(0.10d0)*chngOflenLtrl1(j)
                   
                elseif (length_clsNessLtrl1(j) .lt. TolLen) then
                   
                   if (cnt_lgclsLtrl1 .ne. (NAEC_Ltrl+1)) then
                      lgcl_Ltrls1(j) = .True.
                      cnt_lgclsLtrl1 = cnt_lgclsLtrl1+1

                   elseif (cnt_lgclsLtrl1 == (NAEC_Ltrl+1)) then
                      continue
                   endif
                   
                endif
                
             elseif (i==3) then
                
                sprNm    = Ltrls2(j)
                len_diff = avg_lngthLtrl - l(sprNm)
                length_clsNessLtrl2(j) = abs(len_diff) / avg_lngthLtrl
                
                if ((length_clsNessLtrl2(j).lt.firstChkPoint) .and. (length_clsNessLtrl2(j).gt.TolLen)) then
                   
                   if(lenChngEachStepLtrl2(j).gt.ONE) lenChngEachStepLtrl2(j)=(0.10d0)*chngOflenLtrl2(j)
                   if(lenChngEachStepLtrl2(j).lt.ONE) lenChngEachStepLtrl2(j)=(0.10d0)*chngOflenLtrl2(j)
                   
                elseif (length_clsNessLtrl2(j) .lt. TolLen) then
                   
                   if (cnt_lgclsLtrl2 .ne. (NAEC_Ltrl+1)) then
                      lgcl_Ltrls2(j) = .True. 
                      cnt_lgclsLtrl2 = cnt_lgclsLtrl2+1

                   elseif (cnt_lgclsLtrl2 == (NAEC_Ltrl+1)) then
                      continue
                   endif
                   
                endif
                
             endif
             
          enddo
          
       enddo
       
       write(62,*) " "
       write(62,*) length_clsNessBsal(1:(NAEC_Bsal+1)),"clsNess for Bsal"
       write(62,*) length_clsNessLtrl1(1:(NAEC_Ltrl+1)),"clsNess for Ltrl1"
       write(62,*) length_clsNessLtrl2(1:(NAEC_Ltrl+1)),"clsNess for Ltrl2"
       
       write(62,*) " "
       write(62,*) lgcl_Bsals(1:(NAEC_Bsal+1)),"lgcl_Bsals bfr"
       write(62,*) lgcl_Ltrls1(1:(NAEC_Ltrl+1)),"lgcl_Ltrls1 bfr"
       write(62,*) lgcl_Ltrls2(1:(NAEC_Ltrl+1)),"lgcl_ltrls2 bfr"
       
       write(62,*) "cnt_lgcls before"
       write(62,*) cnt_lgclsBsal,cnt_lgclsLtrl1,cnt_lgclsLtrl2,"cnt_lgcls"
       
       if ((cnt_lgclsBsal.ne.0) .and. (cnt_lgclsBsal.ne.(NAEC_Bsal+1))) then
          write(62,*) "Every logicals (Basals) are not the same"
          write(62,*) lgcl_Bsals(1:(NAEC_Bsal+1)),"T?F"
          
          cnt_lgclsBsal = 0
          lgcl_Bsals(1:(NAEC_Bsal+1)) = .False.
       endif
       
       if ((cnt_lgclsLtrl1.ne.0) .and. (cnt_lgclsLtrl1.ne.(NAEC_Ltrl+1))) then
          write(62,*) "Every logicals (Laterals1) are not the same"
          write(62,*) lgcl_Ltrls1(1:(NAEC_Ltrl+1)),"T?F"
          
          cnt_lgclsLtrl1 = 0
          lgcl_Ltrls1(1:(NAEC_Ltrl+1)) = .False.
       endif
       
       if ((cnt_lgclsLtrl2.ne.0) .and. (cnt_lgclsLtrl2.ne.(NAEC_Ltrl+1))) then
          write(62,*) "Every logicals (Laterals2) are not the same"
          write(62,*) lgcl_Ltrls2(1:(NAEC_Ltrl+1)),"T?F"
          
          cnt_lgclsLtrl2 = 0
          lgcl_Ltrls2(1:(NAEC_Ltrl+1)) = .False.
       endif
       
       write(62,*) " "
       write(62,*) lgcl_Bsals(1:(NAEC_Bsal+1)),"lgcl_Bsals aft"
       write(62,*) lgcl_Ltrls1(1:(NAEC_Ltrl+1)),"lgcl_Ltrls1 aft"
       write(62,*) lgcl_Ltrls2(1:(NAEC_Ltrl+1)),"lgcl_ltrls2 aft"
       
       write(62,*) " "
       write(62,*) cnt_lgclsBsal,cnt_lgclsLtrl1,cnt_lgclsLtrl2,"cnt_lgcls after"
       
       N_postvLgcls = cnt_lgclsBsal + cnt_lgclsLtrl1 + cnt_lgclsLtrl2
       write(62,*) N_postvLgcls,"N_postvLgcls"
       
       if (N_postvLgcls==N_memBsLt1Lt2) then
          write(*,*) "All the logicals are True ..."
          write(*,*) lgcl_Bsals(1:(NAEC_Bsal+1)),"T/F Bsals"
          write(*,*) lgcl_Ltrls1(1:(NAEC_Ltrl+1)),"T/F Ltrl1"
          write(*,*) lgcl_Ltrls2(1:(NAEC_Ltrl+1)),"T/F Ltrl2"
          exit
       endif
       
       cnt_insdLoop = cnt_insdLoop+1
       if (cnt_insdLoop == 10) exit
       
    enddo
    
    close(62)
    stop
    
  end subroutine convert_oneside_tilted_trapezoid_to_recangle_with_same_area
  
  
  
  subroutine create_samel0Chnge_inRight(cell_no)
    implicit none
    integer, intent(in) :: cell_no
    
    integer :: nsprsInEachCell
    integer :: i,imax
    integer :: spr_lft,spr_rgt,sprNmToBeAdded
    
    if (cell_no .le. Hlf_Ncell) then
       write(*,*) "cell_no is less than N_cell",""
    endif
    
    open(unit=127,file='createSamel0Chng.dat',position='append')
    
    nsprsInEachCell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    sprNmToBeAdded  = (Hlf_Ncell)*(nsprsInEachCell)
    
    write(127,*) nsprsInEachCell,sprNmToBeAdded,"nsprsInEachCell,sprNmToBeAdded"
    
    imax = area_spr(cell_no,0)
    
    write(127,*) cell_no,imax,"cell_no,imax"
    
    do i = 1,imax
       
       spr_lft = area_spr(cell_no,i)
       spr_rgt = spr_lft + sprNmToBeAdded
       
       write(127,*) spr_lft,spr_rgt,"spr lft-rgt"
       l0(spr_rgt) = l0(spr_lft)
       
    enddo
    
    close(127)
    
  end subroutine create_samel0Chnge_inRight
  
  
  subroutine change_AR_of_Rectangle(cell_No)
    implicit none
    integer,intent(in) :: cell_no
    
    real*8  :: target_AR
    real*8  :: curr_AR,rect_len,rect_wid
    real*8  :: ONE,HMUCH
    integer :: AR_dcsn
    
    integer, allocatable :: Apcls(:),Bsals(:),Ltrls1(:),Ltrls2(:)

    allocate(Apcls(1:(NAEC_Apcl+1)),Bsals(1:(NAEC_Bsal+1)),Ltrls1(1:(NAEC_Ltrl+1)),Ltrls2(1:(NAEC+1)))
    Apcls = 0 ; Bsals = 0 ; Ltrls1 = 0 ; Ltrls2 = 0
    
    call get_the_apcls_bsals_and_ltrls12(cell_No,Apcls,Bsals,Ltrls1,Ltrls2,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl)
    
    ONE = 1.00d0
    call calc_current_AR(cell_no,curr_AR,rect_len,rect_wid)
    target_AR = AR_cntrl ! Important
    
    if ((curr_AR .gt. ONE) .and. (target_AR .gt. ONE)) then
       continue
    else
       write(*,*) "not equipped with less than ONE"
    endif
    
    if (curr_AR .gt. target_AR) then
       write(*,*) "have to make curr_AR small, means width small,"
       AR_dcsn = -1
    elseif (curr_AR .lt.target_AR) then
       write(*,*) "have to make curr AR big"
       AR_dcsn = +1
    endif
    
    l0_str(1:N_spr) = l0(1:N_spr) !storing l0's
    
    call take_A_good_guess_for_finding_HMUCH(rect_len,rect_wid,curr_AR,target_AR,AR_dcsn,HMUCH)
    call change_l0s_for_getting_AR(AR_dcsn,HMUCH,Apcls,Bsals,Ltrls1,Ltrls2,NAEC_Apcl,NAEC_Bsal,NAEC_Ltrl)
    
    
  end subroutine change_AR_of_Rectangle
  
  
  subroutine get_the_apcls_bsals_and_ltrls12(cell_no,Apcls,Bsals,Ltrls1,Ltrls2,N_Ap,N_Bs,N_Lt)
    implicit none
    integer, intent(in)  :: cell_no
    integer, intent(in)  :: N_Ap,N_Bs,N_Lt
    integer, intent(out) :: Apcls(1:(N_Ap+1)),Bsals(1:(N_Bs+1)),Ltrls1(1:(N_Lt+1)),Ltrls2(1:(N_Lt+1))
    
    integer :: i,j,jmax
    integer :: N_itm,nsprs_InEachCell
    
    N_itm            = 4
    nsprs_InEachCell = (N_Ap+1) + (N_Bs+1) + (N_Lt+1)
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_Ap+1
       if (i==2) jmax = N_Bs+1
       if (i==3) jmax = N_Lt+1
       if (i==4) jmax = N_Lt+1
       
       do j = 1,jmax
          
          if (i==1) Apcls(j)  = (nsprs_InEachCell)*(cell_no-1) + j
          if (i==2) Bsals(j)  = (nsprs_InEachCell)*(cell_no-1) + (N_Ap+1) + j
          if (i==3) Ltrls1(j) = (nsprs_InEachCell)*(cell_no-2) + (N_Ap+1) + (N_Bs+1) + j 
          if (i==4) Ltrls2(j) = (nsprs_InEachCell)*(cell_no-1) + (N_Ap+1) + (N_Bs+1) + j
          
       enddo
       
    enddo
    
  end subroutine get_the_apcls_bsals_and_ltrls12
  
  
  subroutine get_the_apcls_bsals_and_oneltrls(cell_no,Apcls,Bsals,Ltrls,N_Ap,N_Bs,N_Lt)
    implicit none
    integer, intent(in)  :: cell_no
    integer, intent(in)  :: N_Ap,N_Bs,N_Lt
    integer, intent(out) :: Apcls(1:(N_Ap+1)),Bsals(1:(N_Bs+1)),Ltrls(1:(N_Lt+1))
    
    integer :: i,j,jmax
    integer :: N_itm,nsprs_InEachCell
    
    N_itm            = 3
    nsprs_InEachCell = (N_Ap+1) + (N_Bs+1) + (N_Lt+1)
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_Ap+1
       if (i==2) jmax = N_Bs+1
       if (i==3) jmax = N_Lt+1
       
       do j = 1,jmax
          
          if (i==1) Apcls(j) = (nsprs_InEachCell)*(cell_no-1) + j
          if (i==2) Bsals(j) = (nsprs_InEachCell)*(cell_no-1) + (N_Ap+1) + j
          if (i==3) Ltrls(j) = (nsprs_InEachCell)*(cell_no-1) + (N_Ap+1) + (N_Bs+1) + j
          
       enddo
       
    enddo
    
  end subroutine get_the_apcls_bsals_and_oneltrls
  
  
  subroutine calc_current_AR(cell_no,curr_AR,rect_len,rect_wid)
    implicit none
    integer, intent(in)    :: cell_no
    real*8 , intent(inout) :: curr_AR
    real*8 , intent(inout) :: rect_len,rect_wid
    
    integer, allocatable :: apcls(:),ltrls1(:)
    integer :: i,j,jmax
    integer :: nsprsInACell,N_itm
    
    allocate(apcls(1:(NAEC_Apcl+1)),ltrls1(1:(NAEC_Ltrl+1)))
    
    N_itm        = 2
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    write(*,*) NAEC_Apcl,NAEC_Ltrl,nsprsInACell,"NAEC and nsprsInACell"
    
    do i = 1,N_itm
       
       if (i==1) jmax = NAEC_Apcl+1
       if (i==2) jmax = NAEC_Ltrl+1
       
       do j = 1,jmax
          
          if (i==1) apcls(j)  = (nsprsInACell)*(cell_no-1) + j
          if (i==2) ltrls1(j) = (nsprsInACell)*(cell_no-2) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + j
          
       enddo
       
    enddo

    write(*,*) apcls(1:(NAEC_Apcl+1)),"Apcls"
    write(*,*) ltrls1(1:(NAEC_Ltrl+1)),"Ltrls1"
    
    rect_len = 0.0d0 ; rect_wid = 0.0d0
    
    do i = 1,N_itm
       
       if (i==1) jmax = NAEC_Apcl+1
       if (i==2) jmax = NAEC_Ltrl+1
       
       do j = 1,jmax
          if (i==1) rect_wid = rect_wid + l(apcls(j))
          if (i==2) rect_len = rect_len + l(ltrls1(j))
       enddo
       
    enddo
    
    write(*,*) rect_len,rect_wid,"len-wid"
    
    curr_AR = rect_len/rect_wid

    write(*,*) curr_AR,"curr_AR"
    
  end subroutine calc_current_AR
  
  subroutine take_A_good_guess_for_finding_HMUCH(rect_len,rect_wid,curr_AR,target_AR,AR_dcsn,HMUCH)
    implicit none
    real*8, intent(in)  :: rect_len,rect_wid
    real*8, intent(in)  :: curr_AR,target_AR
    integer,intent(in)  :: AR_dcsn
    real*8, intent(out) :: HMUCH
    
    real*8  :: initial_guess
    real*8  :: imprv_rectLen,imprv_rectWid
    real*8  :: imprv_AR
    integer :: imprv_ARdcsn
    
    initial_guess = 0.05d0
    
    do 
       
       if (AR_dcsn == +1) then
          imprv_rectLen = rect_len + initial_guess*rect_len
          imprv_rectWid = rect_wid - initial_guess*rect_wid
          
       elseif (AR_dcsn == -1) then
          imprv_rectLen = rect_len - initial_guess*rect_len 
          imprv_rectWid = rect_wid + initial_guess*rect_wid
       endif
       
       imprv_AR = imprv_rectLen/imprv_rectWid
       
       if (imprv_AR .gt. target_AR) imprv_ARdcsn = -1
       if (imprv_AR .lt. target_AR) imprv_ARdcsn = +1
       
       if (AR_dcsn == imprv_ARdcsn) then
          write(*,*) "It means we can work with the initial guess"
          exit
       elseif (AR_dcsn .ne. imprv_ARdcsn) then
          write(*,*) "We cant work the initial guess"
          initial_guess = initial_guess/2.0d0
       endif
       
    enddo
    
    HMUCH = initial_guess
    
  end subroutine take_A_good_guess_for_finding_HMUCH
  
  
  subroutine change_l0s_for_getting_AR(AR_dcsn,HMUCH,Apcls,Bsals,Ltrls1,Ltrls2,N_Ap,N_Bs,N_Lt)
    implicit none
    integer, intent(in) :: AR_dcsn
    real*8 , intent(in) :: HMUCH
    integer, intent(in) :: N_Ap,N_Bs,N_Lt
    integer, intent(in) :: Apcls(1:N_Ap),Bsals(1:N_Bs),Ltrls1(1:N_Lt),Ltrls2(1:N_Lt)
    
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 4
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_Ap+1
       if (i==2) jmax = N_Bs+1
       if (i==3) jmax = N_Lt+1
       if (i==4) jmax = N_Lt+1
       
       do j = 1,jmax
          
          if (AR_dcsn == -1) then
             
             if (i==1) l0(Apcls(j))  = (1.0d0-HMUCH)*l0(Apcls(j))
             if (i==2) l0(Bsals(j))  = (1.0d0-HMUCH)*l0(Bsals(j))
             if (i==3) l0(Ltrls1(j)) = (1.0d0+HMUCH)*l0(Ltrls1(j))
             if (i==4) l0(Ltrls2(j)) = (1.0d0+HMUCH)*l0(Ltrls2(j))
             
          elseif (AR_dcsn == +1) then
             
             if (i==1) l0(Apcls(j))  = (1.0d0+HMUCH)*l0(Apcls(j))
             if (i==2) l0(Bsals(j))  = (1.0d0+HMUCH)*l0(Bsals(j))
             if (i==3) l0(Ltrls1(j)) = (1.0d0-HMUCH)*l0(Ltrls1(j))
             if (i==4) l0(Ltrls2(j)) = (1.0d0-HMUCH)*l0(Ltrls2(j))
             
          endif
          
       enddo
       
    enddo
    
    
  end subroutine change_l0s_for_getting_AR
  
  subroutine make_ManpltArea_invaginatedRegion_basedOnTopLayerCell(cell_no,readingDcsn,strctToRead,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cell_no
    integer, intent(in)    :: readingDcsn
    integer, intent(in)    :: strctToRead,cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E,AreaOfFrstCell
    real*8  :: AreaDecr
    real*8  :: ApprxAreaChngCell
    real*8  :: firstChkPoint,TolArea
    integer :: cnt_TrueLgcls
    
    integer :: i,j
    real*8  :: Area_diff
    real*8  :: ZERO
    integer :: cellNm
    real*8  :: AreaIncrPrcnt
    
    real*8,  allocatable :: PresVal(:),TnsnVal(:)
    logical, allocatable :: lgcl_Chng(:)
    real*8 , allocatable :: Area_clsNess(:)
    real*8 , allocatable :: AreaIncrEachStep(:) 
    integer, allocatable :: countForFrstChkPoint(:)
    
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
    
    
    open(unit=74,file='area_Manpltn.dat',position='append')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
       
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    
    if (readingDcsn == 1) then
       call read_strctProps_withAddedCell(strctToRead)
       
       A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
       k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
       l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
       k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
       CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
       CgYNode(1:N_node) = 0.0d0
       
    elseif (readingDcsn == 0) then
       continue
    endif
       
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(74,*) (FrmNo-1),"bfr chaging A0 of rectangular cells"
    
    AreaOfFrstCell = A(1)
    
    if (cell_no == (Hlf_Ncell-2)) AreaDecr = 1.65d0
    if (cell_no == (Hlf_Ncell-0)) AreaDecr = 1.65d0   
    if (cell_no == N_cell)        AreaDecr = 1.43d0
    
    if (cell_no .ne. N_cell) AreaIncrPrcnt = 1.05d0 !write A routine that will decide automatically guess this value
    if (cell_no == N_cell)   AreaIncrPrcnt = 1.01d0
    
    ApprxAreaChngCell = A(1)/AreaDecr
    firstChkPoint = 0.10d0 ; TolArea = 0.01d0
    
    if (cell_no .ne. N_cell) then
       N_chngCells = 2
       allocate(ChngCells(1:N_chngCells))
       ChngCells(1) = cell_no
       ChngCells(2) = Hlf_Ncell+cell_no
       
    elseif (cell_no == N_cell) then
       N_chngCells = 1
       allocate(ChngCells(1:N_chngCells))
       ChngCells(1) = N_cell
    endif
    
    allocate(lgcl_Chng(1:N_chngCells))
    allocate(Area_clsNess(1:N_chngCells))
    allocate(AreaIncrEachStep(1:N_chngCells))
    allocate(countForFrstChkPoint(1:N_chngCells))
    
    lgcl_Chng(1:N_chngCells)            = .False.
    AreaIncrEachStep(1:N_chngCells)     = AreaIncrPrcnt
    countForFrstChkPoint(1:N_chngCells) = 0
    cnt_TrueLgcls                       = 0
    
    do
       
       do i = 1,N_chngCells
          
          cellNm = ChngCells(i)         
          
          if (lgcl_Chng(i) .eqv. .False.) then
             A0(cellNm) = A0(cellNm)*AreaIncrEachStep(i)
             
          elseif (lgcl_Chng(i) .eqv. .True.) then
             continue
          endif
          
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       do i = 1,N_chngCells
          
          cellNm = ChngCells(i)
          
          Area_diff       = ApprxAreaChngCell-A(cellNm)
          Area_clsNess(i) = Area_diff/ApprxAreaChngCell
          
          if (Area_diff.lt.ZERO) then
             write(*,*) "Area increased too much for cellNm=",cellNm,A(cellNm),A0(cellNm)
          endif
          
          if ((Area_clsNess(i) .lt. firstChkPoint) .and. (Area_clsNess(i).gt.TolArea)) then
             
             if (countForFrstChkPoint(i) == 0) then
                AreaIncrEachStep(i) = 1.005d0
                countForFrstChkPoint(i) = 1
                
             elseif (countForFrstChkPoint(i) == 1) then
                continue
             endif
             
          elseif (Area_clsNess(i) .lt. TolArea) then
             AreaIncrEachStep(i) = 1.00d0
             lgcl_Chng(i) = .True.
          endif
          
       enddo
       
       
       do i = 1,N_chngCells
          
          if (lgcl_Chng(i) .eqv. .True.) then
             cnt_TrueLgcls = cnt_TrueLgcls+1
          elseif (lgcl_Chng(i) .eqv. .False.) then
             continue
          endif
          
       enddo
       
       if (cnt_TrueLgcls == N_chngCells) then
          exit
       elseif (cnt_TrueLgcls .lt. N_chngCells) then
          cnt_TrueLgcls = 0
       elseif (cnt_TrueLgcls .gt. N_chngCells) then
          write(*,*) "cnt_TrueLgcls cant be greater than N_chngCells"
       endif
       
    enddo
    
    deallocate(ChngCells)
    
    close(74)
    
  end subroutine make_ManpltArea_invaginatedRegion_basedOnTopLayerCell
  
  
  subroutine BorrowCellPropsFromOneCS_And_ManpltOtherCS(strctFrm,strctTo,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctFrm,strctTo
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: strctToRead
    real*8  :: E
    integer :: i,imax,j,jmaxL,jmaxR,jmax
    integer :: cellL,cellR
    integer :: sprNmL,sprNmR
    integer :: nsprsInACell
    integer :: NcellsToCopy
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
    real*8, allocatable :: A0_strctFrm(:),ka_strctFrm(:)
    real*8, allocatable :: l0_strctFrm(:),ks_strctFrm(:)
    real*8, allocatable :: CgX_strctFrm(:),CgY_strctFrm(:)
    
    real*8, allocatable :: A0_cs1(:),A0_cs2(:),A0_cs3(:),A0_cs4(:)
    real*8, allocatable :: ka_cs1(:),ka_cs2(:),ka_cs3(:),ka_cs4(:)
    real*8, allocatable :: l0_cs1(:),l0_cs2(:),l0_cs3(:),l0_cs4(:)
    real*8, allocatable :: ks_cs1(:),ks_cs2(:),ks_cs3(:),ks_cs4(:)
    !real*8, allocatable :: CgX_cs1(:),CgX_cs2(:),CgX_cs3(:),CgX_cs4(:)
    !real*8, allocatable :: nodeXY_cs1(:,:),nodeXY_cs2(:,:),nodeXY_cs3(:,:),nodeXY_cs4(:,:)
    
    integer, allocatable :: cellNos(:)
    
    integer :: CS1,CS2,CS3,CS4
    integer :: cellNoLInCS1,cellNoLInCS4,cellNoRInCS1,cellNoRInCS4
    integer :: jmaxFrmCS1,jmaxFrmCS4
    integer :: sprNoLInCS1,sprNoLInCS4,sprNoRInCS1,sprNoRInCS4
    integer :: N_cellsChngd
    integer :: cellNoL,cellNoR,sprNoL,sprNoR

    integer :: saveWhom
    real*8  :: weight,weight1,weight2
    
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
    
    open(unit=78,file='BorrowCellProps.dat')
    open(unit=68,file='BeforeBorrowingProps.dat')
    open(unit=58,file='CellSprNos.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
       
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    allocate(A0_strctFrm(1:N_cell),ka_strctFrm(1:N_cell))
    allocate(l0_strctFrm(1:N_spr) ,ks_strctFrm(1:N_spr) )
    allocate(CgX_strctFrm(1:N_node),CgY_strctFrm(1:N_node))
    
    strctToRead = strctTo
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(78,*) (FrmNo-1),"beginning of BorrowCells"
    
    strctToRead = strctFrm
    call read_strctProps_withAddedCell(strctToRead)
    
    A0_strctFrm(1:N_cell)  = A0_strctNI(strctToRead,1:N_cell)
    ka_strctFrm(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0_strctFrm(1:N_spr)   = l0_strctNI(strctToRead,1:N_spr)
    ks_strctFrm(1:N_spr)   = ks_strctNI(strctToRead,1:N_spr)
    CgX_strctFrm(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgY_strctFrm(1:N_node) = 0.0d0
    
    ! *** PRINTING VALUES BEFORE COPYING *** !
    
    
    NcellsToCopy = 4
    allocate(cellNos(1:NcellsToCopy))
    
    cellNos(1) = Hlf_Ncell - 2
    cellNos(2) = Hlf_Ncell - 1
    cellNos(3) = Hlf_Ncell - 0
    cellNos(4) = N_cell
    
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    do i = 1,NcellsToCopy
       
       cellL = cellNos(i)
       
       if (cellL .ne.  N_cell) then
          
          cellR = cellNos(i) + Hlf_Ncell
          
          write(68,*) A0(cellL),A0(cellR),cellL,cellR,"A0s Bfr"
          write(68,*) k_area(cellL),k_area(cellR),cellL,cellR,"kaS Bfr"
          write(68,*) " "
          
          jmaxL = area_spr(cellL,0)
          jmaxR = area_spr(cellR,0)
          
          if (jmaxL .ne. jmaxR) write(68,*) "jmaxL is not equal to jmaxR",jmaxL,jmaxR
          
          do j = 1,jmaxL
             
             sprNmL = area_spr(cellL,j) ; sprNmR = area_spr(cellR,j)
             
             write(68,*) l0(sprNmL),l0(sprNmR),sprNmL,sprNmR,"l0s Bfr"
             write(68,*) k_spr(sprNmL),k_spr(sprNmR),sprNmL,sprNmR,"ksS Bfr"
             write(68,*) " "
             
          enddo
          
          
       elseif (cellL == N_cell) then
          
          write(68,*) A0(cellL),cellL,"A0s Bfr at N_cell"
          write(68,*) k_area(cellL),cellL,"kaS Bfr at N_cell"
          write(68,*) " "
          
          jmaxL = area_spr(cellL,0)
          
          do j = 1,jmaxL
             
             sprNmL = area_spr(cellL,j)
             
             write(68,*) l0(sprNmL),sprNmL,"l0s Bfr at N_cell"
             write(68,*) k_spr(sprNmL),sprNmL,"ksS Bfr at N_cell"
             write(68,*) " "
             
          enddo
          
       endif
       
    enddo
    
    
    ! *** PRINTING FINISHED  *** !
    
    ! This routine is for copying properties of three bottom pair cells and one initiator cell from control state 4 to control state 1 and predict an average of those (l0 and A0) for control state 2 and control state 3 (This is a classic example of an extremely non-generic routine, I hate these routines !!! )
    
    
    CS1 = 13 ; call read_strctPropsIncludingNodeXY_withAddedCell(CS1)
    CS2 = 14 ; call read_strctPropsIncludingNodeXY_withAddedCell(CS2)
    CS3 = 15 ; call read_strctPropsIncludingNodeXY_withAddedCell(CS3) 
    CS4 = 16 ; call read_strctPropsIncludingNodeXY_withAddedCell(CS4)
    
    allocate(A0_cs1(1:N_cell),A0_cs2(1:N_cell),A0_cs3(1:N_cell),A0_cs4(1:N_cell))
    allocate(ka_cs1(1:N_cell),ka_cs2(1:N_cell),ka_cs3(1:N_cell),ka_cs4(1:N_cell))
    allocate(l0_cs1(1:N_spr) ,l0_cs2(1:N_spr) ,l0_cs3(1:N_spr) ,l0_cs4(1:N_spr) )
    allocate(ks_cs1(1:N_spr) ,ks_cs2(1:N_spr) ,ks_cs3(1:N_spr) ,ks_cs4(1:N_spr) )
    
    A0_CS1(1:N_cell) = A0_strctNI(CS1,1:N_cell) ; ka_CS1(1:N_cell) = ka_strctNI(CS1,1:N_cell)
    A0_CS2(1:N_cell) = A0_strctNI(CS2,1:N_cell) ; ka_CS2(1:N_cell) = ka_strctNI(CS2,1:N_cell)
    A0_CS3(1:N_cell) = A0_strctNI(CS3,1:N_cell) ; ka_CS3(1:N_cell) = ka_strctNI(CS3,1:N_cell)
    A0_CS4(1:N_cell) = A0_strctNI(CS4,1:N_cell) ; ka_CS4(1:N_cell) = ka_strctNI(CS4,1:N_cell)
    
    l0_CS1(1:N_spr) = l0_strctNI(CS1,1:N_spr) ; ks_CS1(1:N_spr) = ks_strctNI(CS1,1:N_spr)
    l0_CS2(1:N_spr) = l0_strctNI(CS2,1:N_spr) ; ks_CS2(1:N_spr) = ks_strctNI(CS2,1:N_spr)
    l0_CS3(1:N_spr) = l0_strctNI(CS3,1:N_spr) ; ks_CS3(1:N_spr) = ks_strctNI(CS3,1:N_spr)
    l0_CS4(1:N_spr) = l0_strctNI(CS4,1:N_spr) ; ks_CS4(1:N_spr) = ks_strctNI(CS4,1:N_spr)
    
    call print_Props_strctPropsAddedCellNI(CS1)
    call print_Props_strctPropsAddedCellNI(CS2)
    call print_Props_strctPropsAddedCellNI(CS3)
    call print_Props_strctPropsAddedCellNI(CS4)
    
    imax = NcellsToCopy-1
    
    do i = 1,imax
       
       if (i==1) then
          
          cellNoLInCS1 = Hlf_Ncell-1             ; cellNoLInCS4 = Hlf_Ncell-2
          cellNoRInCS1 = cellNoLInCS1+Hlf_Ncell  ; cellNoRInCS4 = cellNoLInCS4+Hlf_Ncell
          
          A0_CS1(cellNoLInCS1) = A0_CS4(cellNoLInCS4) ; A0_CS1(cellNoRInCS1) = A0_CS4(cellNoRInCS4)
          
          write(58,*) cellNoLInCS1,cellNoRInCS1,cellNoLInCS4,cellNoRInCS4,i,"cellNoL/RInCS1-4,i"
          write(58,*) A0_CS1(cellNoLInCS1),A0_CS4(cellNoLInCS4),A0_CS1(cellNoRInCS1),A0_CS4(cellNoRInCS4),"A0s"

          jmaxFrmCS1 = area_spr(cellNoLInCS1,0) ; jmaxFrmCS4 = area_spr(cellNoLInCS4,0)
          write(58,*) jmaxFrmCS1,jmaxFrmCS4,"jmax in i=1"
          
       elseif (i==2) then
          
          cellNoLInCS1 = Hlf_Ncell              ; cellNoLInCS4 = Hlf_Ncell
          cellNoRInCS1 = cellNoLInCS1+Hlf_Ncell ; cellNoRInCS4 = cellNoLInCS4+Hlf_Ncell
          
          A0_CS1(cellNoLInCS1) = A0_CS4(cellNoLInCS4) ; A0_CS1(cellNoRInCS1) = A0_CS4(cellNoRInCS4)
          
          write(58,*) cellNoLInCS1,cellNoRInCS1,cellNoLInCS4,cellNoRInCS4,i,"cellNoL/RInCS1-4,i"
          write(58,*) A0_CS1(cellNoLInCS1),A0_CS4(cellNoLInCS4),A0_CS1(cellNoRInCS1),A0_CS4(cellNoRInCS4),"A0s"
          
          jmaxFrmCS1 = nsprsInACell ; jmaxFrmCS4 = nsprsInACell
          write(58,*) jmaxFrmCS1,jmaxFrmCS4,"jmax in i=2"
          
       elseif (i==3) then
          
          cellNoLInCS1 = N_cell      ; cellNoLInCS4 = N_cell
          A0_CS1(cellNoLInCS1) = A0_CS4(cellNoLInCS4)
          
          write(58,*) cellNoLInCS1,cellNoLInCS4,i,"cellNoLInCS1-4,i"
          write(58,*) A0_CS1(cellNoLInCS1),A0_CS4(cellNoLInCS4),"A0s"

          jmaxFrmCS1 = NAEC_Ltrl+1 ; jmaxFrmCS4 = NAEC_Ltrl+1
          write(58,*) jmaxFrmCS1,jmaxFrmCS4,"jmax in i=3"
          
       endif
       
       if (jmaxFrmCS1 .ne. jmaxFrmCS4) then
          write(*,*) "jmaxS are not equal",jmaxFrmCS1,jmaxFrmCS4,i
          stop
       endif
       
       
       do j = 1,jmaxFrmCS1
          
          if (i==1) then
             
             sprNoLInCS1 = area_spr(cellNoLInCS1,j) ; sprNoLInCS4 = area_spr(cellNoLInCS4,j) 
             sprNoRInCS1 = area_spr(cellNoRInCS1,j) ; sprNoRInCS4 = area_spr(cellNoRInCS4,j)
             
             l0_CS1(sprNoLInCS1) = l0_CS4(sprNoLInCS4) ; l0_CS1(sprNoRInCS1) = l0_CS4(sprNoRInCS4) 
             
             write(58,*) sprNoLInCS1,sprNoRInCS1,sprNoLInCS4,sprNoRInCS4,j,"sprNoL/RInCS1-4,j"
             write(58,*) l0_CS1(sprNoLInCS1),l0_CS4(sprNoLInCS4),l0_CS1(sprNoRInCS1),l0_CS4(sprNoRInCS4),"l0s"
             
          elseif (i==2) then
             
             sprNoLInCS1 = (cellNoLInCS1-1)*(nsprsInACell) + j ; sprNoLInCS4 = (cellNoLInCS4-1)*(nsprsInACell) + j
             sprNoRInCS1 = (cellNoRInCS1-1)*(nsprsInACell) + j ; sprNoRInCS4 = (cellNoRInCS4-1)*(nsprsInACell) + j
             
             l0_CS1(sprNoLInCS1) = l0_CS4(sprNoLInCS4) ; l0_CS1(sprNoRInCS1) = l0_CS4(sprNoRInCS4)

             write(58,*) sprNoLInCS1,sprNoRInCS1,sprNoLInCS4,sprNoRInCS4,j,"sprNoL/RInCS1-4,j"
             write(58,*) l0_CS1(sprNoLInCS1),l0_CS4(sprNoLInCS4),l0_CS1(sprNoRInCS1),l0_CS4(sprNoRInCS4),"l0s"
             
          elseif (i==3) then
             
             sprNoLInCS1 = (cellNoLInCS1-1)*(nsprsInACell) + j ; sprNoLInCS4 = (cellNoLInCS4-1)*(nsprsInACell) + j
             l0_CS1(sprNoLInCS1) = l0_CS4(sprNoLInCS4)
             
             write(58,*) sprNoLInCS1,sprNoLInCS4,j,"sprNoLInCS1-4,j"
             write(58,*) l0_CS1(sprNoLInCS1),l0_CS4(sprNoLInCS4),"l0s"
             
          endif
          
       enddo
       
    enddo

    write(58,*) " "
    write(58,*) " "
    write(58,*) "CS2 and CS3"
    
    ! CS1 is edited upto this, now we have to edit the properties for CS2 and CS3
    
    N_cellsChngd = 4
    weight1 = 0.50d0 ; weight2 = 0.00d0
    
    do i = 1,N_cellsChngd
       
       if (i.le.2) weight = weight1
       if ((i.gt.2) .and. (i.le.N_cellsChngd)) weight = weight2
       
       if (i.ne.N_cellsChngd) then
          
          cellNoL = Hlf_Ncell-2 + (i-1)
          cellNoR = cellNoL + Hlf_Ncell
          
          A0_CS2(cellNoL) = A0_CS1(cellNoL) + (weight)*(A0_CS4(cellNoL) - A0_CS1(cellNoL)) ; A0_CS3(cellNoL) = A0_CS2(cellNoL)
          A0_CS2(cellNoR) = A0_CS1(cellNoR) + (weight)*(A0_CS4(cellNoR) - A0_CS1(cellNoR)) ; A0_CS3(cellNoR) = A0_CS2(cellNoR)
          
          write(58,*) " "
          write(58,*) cellNoL,cellNoR,i,"cellnos"
          
          write(58,*) A0_CS1(cellNoL),A0_CS4(cellNoL),cellNoL,"A0_CS1/CS4 to calc CS2/CS3 L"
          write(58,*) A0_CS1(cellNoR),A0_CS4(cellNoR),cellNoR,"A0_CS1/CS4 to calc CS2/CS3 R"
          write(58,*) A0_CS2(cellNoL),A0_CS3(cellNoL),i,"A0s in CS2/3 in L"
          write(58,*) A0_CS2(cellNoR),A0_CS3(cellNoR),i,"A0s in CS2/3 in R"
          
       elseif (i==N_cellsChngd) then
          
          cellNoL = N_cell
          A0_CS2(cellNoL) = A0_CS1(cellNoL) + (weight)*(A0_CS4(cellNoL) - A0_CS1(cellNoL)) ; A0_CS3(cellNoL) = A0_CS2(cellNoL)
          
          write(58,*) " "
          write(58,*) cellNoL,i,"cellnos"
          
          write(58,*) A0_CS1(cellNoL),A0_CS4(cellNoL),"A0_CS1/CS4 to calc CS2/CS3"
          write(58,*) A0_CS2(cellNoL),A0_CS3(cellNoL),i,"A0s in CS2/3 in L"
          
       endif
       
       if (i==1)               jmax = area_spr(cellNoL,0)
       if ((i==2) .or. (i==3)) jmax = nsprsInACell 
       if (i==N_cellsChngd)    jmax = NAEC_Ltrl+1
       
       do j = 1,jmax
          
          if (i==1) then
             
             sprNoL = area_spr(cellNoL,j) ; sprNoR = area_spr(cellNoR,j)
             
             l0_CS2(sprNoL) = l0_CS1(sprNoL) + (weight)*(l0_CS4(sprNoL) - l0_CS1(sprNoL)) ; l0_CS3(sprNoL) = l0_CS2(sprNoL)
             l0_CS2(sprNoR) = l0_CS1(sprNoR) + (weight)*(l0_CS4(sprNoR) - l0_CS1(sprNoR)) ; l0_CS3(sprNoR) = l0_CS2(sprNoR)
             
             write(58,*) sprNoL,sprNoR,j,"sprnos"
             write(58,*) l0_CS2(sprNoL),l0_CS3(sprNoL),j,"l0s in CS2/3 in L"
             write(58,*) l0_CS2(sprNoR),l0_CS3(sprNoR),j,"l0s in CS2/3 in R"
             
          elseif ((i==2) .or. (i==3)) then
             
             sprNoL = (cellNoL-1)*(nsprsInACell)+j ; sprNoR = (cellNoR-1)*(nsprsInACell)+j
             
             l0_CS2(sprNoL) = l0_CS1(sprNoL) + (weight)*(l0_CS4(sprNoL) - l0_CS1(sprNoL)) ; l0_CS3(sprNoL) = l0_CS2(sprNoL)
             l0_CS2(sprNoR) = l0_CS1(sprNoR) + (weight)*(l0_CS4(sprNoR) - l0_CS1(sprNoR)) ; l0_CS3(sprNoR) = l0_CS2(sprNoR)
             
             write(58,*) sprNoL,sprNoR,j,"sprnos"
             write(58,*) l0_CS2(sprNoL),l0_CS3(sprNoL),j,"l0s in CS2/3 in L"
             write(58,*) l0_CS2(sprNoR),l0_CS3(sprNoR),j,"l0s in CS2/3 in R"
             
          elseif (i==N_cellsChngd) then
             
             sprNoL = (cellNoL-1)*(nsprsInACell)+j
             l0_CS2(sprNoL) = l0_CS1(sprNoL) + (weight)*(l0_CS4(sprNoL) - l0_CS1(sprNoL)) ; l0_CS3(sprNoL) = l0_CS2(sprNoL)
             
             write(58,*) sprNoL,j,"sprnos"
             write(58,*) l0_CS2(sprNoL),l0_CS3(sprNoL),j,"l0s in CS2/3 in L"
             
          endif
          
       enddo
       
    enddo
    
    
    A0_strctNI(CS1,1:N_cell) = A0_CS1(1:N_cell)
    A0_strctNI(CS2,1:N_cell) = A0_CS2(1:N_cell)
    A0_strctNI(CS3,1:N_cell) = A0_CS3(1:N_cell)
    A0_strctNI(CS4,1:N_cell) = A0_CS4(1:N_cell)
    
    l0_strctNI(CS1,1:N_spr) = l0_CS1(1:N_spr)
    l0_strctNI(CS2,1:N_spr) = l0_CS2(1:N_spr)
    l0_strctNI(CS3,1:N_spr) = l0_CS3(1:N_spr)
    l0_strctNI(CS4,1:N_spr) = l0_CS4(1:N_spr)
    
    !saveWhom = 3 ; call saveA0l0_ofStrct_AddedCell_TNorNI(CS1,saveWhom)
    !saveWhom = 3 ; call saveA0l0_ofStrct_AddedCell_TNorNI(CS2,saveWhom)
    !saveWhom = 3 ; call saveA0l0_ofStrct_AddedCell_TNorNI(CS3,saveWhom)
    !saveWhom = 3 ; call saveA0l0_ofStrct_AddedCell_TNorNI(CS4,saveWhom)
    
    call saveA0l0_ofStrct_AddedCell_withknownStrct(CS1)
    call saveA0l0_ofStrct_AddedCell_withknownStrct(CS2)
    call saveA0l0_ofStrct_AddedCell_withknownStrct(CS3)
    call saveA0l0_ofStrct_AddedCell_withknownStrct(CS4)
    
    call Find_Analytical_and_Numerical_Mismatch
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(58)
    close(68)
    close(78)
    
  end subroutine BorrowCellPropsFromOneCS_And_ManpltOtherCS
  
  
  subroutine print_Props_strctPropsAddedCellNI(CS_No)
    implicit none
    integer, intent(in) :: CS_No
    
    character(len=100)  :: flnm
    character(len=100)  :: flnmbr
    character(len=100)  :: full_flnm

    integer :: i,j,jmax
    integer :: N_itm
    
    if (modelID.ne.2) then
       write(*,*) "This routine is compatible for modelID = 2"
       stop
    endif
    
    flnm='strctPropsAddedCellNIPrint'
    write(flnmbr,'(i2.2,a)') CS_No,'S1T1.dat' 
    
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr)) 
    write(*,*)trim(adjustl(full_flnm))
    
    open(unit=87,file=trim(adjustl(full_flnm)))
    
    N_itm = 4
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) write(87,*) ks_strctNI(CS_No,j),l0_strctNI(CS_No,j),j
          if (i==2) write(87,*) ka_strctNI(CS_No,j),A0_strctNI(CS_No,j),j
          if (i==3) write(87,*) CgX_strctNI(CS_No,j),CgY_strctNI(CS_No,j)
          if (i==4) write(87,*) nodeXY_strctNI(CS_No,j,1:N_dmnsn)
          
       enddo
       
       write(87,*) " "
       
    enddo
    
    close(87)
    
  end subroutine print_Props_strctPropsAddedCellNI
  
  
  subroutine making_incoming_cells_identical(strctToRead,howManyIC,cntrlSt,ExpNo,FrmNo)
    !IC = IncomingCells,not Initiator Cell 
    implicit none
    integer, intent(in)    :: howManyIC
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: IICL,IICR !Inner Incoming Cell left/Right
    integer :: i,j,m,imax,mmax
    real*8  :: A0_IICL,A0_IICR
    real*8  :: ka_IICL,ka_IICR
    integer :: sprCntEachCell
    integer :: sprNm,sprIICL,sprIICR
    integer :: cellL,cellR
    integer :: saveWhom
    
    real*8,  allocatable :: PresVal(:),TnsnVal(:)
    integer, allocatable :: spr_IICL(:),spr_IICR(:)
    
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
    
    open(unit=63,file='identicalCells.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(63,*) (FrmNo-1),"bfr Identical Incoming cell props"
    
    imax = howManyIC
    
    IICL = howManyIC       ; IICR = Hlf_Ncell + howManyIC
    A0_IICL = A0(IICL)     ; A0_IICR = A0(IICR)
    ka_IICL = k_area(IICL) ; ka_IICR = k_area(IICR)
    
    write(63,*) IICL,IICR,"IICL/R"
    write(63,*) A0_IICL,A0_IICR,ka_IICL,ka_IICR,"A0,ka of IICL/R"
    
    sprCntEachCell = 1 + (NAEC+1)*1 + (NAEC+1)*1 ! Apcl + Bsal + Ltrl
    allocate(spr_IICL(1:sprCntEachCell),spr_IICR(1:sprCntEachCell))
    
    do m = 1,sprCntEachCell
       spr_IICL(m) = sprCntEachCell*(IICL-1) + m
       spr_IICR(m) = sprCntEachCell*(IICR-1) + m
    enddo
    
    do i = 1,(imax-1)
       
       cellL = i ; cellR = Hlf_Ncell+i
       
       do j = 1,2
          
          if (j==1) then
             k_area(cellL) = ka_IICL
             A0(cellL)     = A0_IICL
          elseif (j==2) then
             k_area(cellR) = ka_IICR
             A0(cellR)     = A0_IICR
          endif
          
          mmax = sprCntEachCell
          
          do m = 1,mmax
             
             if (j==1) then
                sprNm   = sprCntEachCell*(cellL-1) + m
                sprIICL = spr_IICL(m)
                
                k_spr(sprNm) = k_spr(sprIICL)
                l0(sprNm)    = l0(sprIICL)
                
                write(63,*) sprNm,sprIICL,m,"spr in current cell,IICL,m"
                write(63,*) l0(sprNm),l0(sprIICL),k_spr(sprNm),k_spr(sprIICL),"l0,ks of spr/sprIICL"
                
             elseif (j==2) then
                sprNm   = sprCntEachCell*(cellR-1) + m
                sprIICR = spr_IICR(m)
                
                k_spr(sprNm) = k_spr(sprIICR)
                l0(sprNm)    = l0(sprIICR)
                
                write(63,*) sprNm,sprIICR,m,"spr in current cell,IICR,m"
                write(63,*) l0(sprNm),l0(sprIICR),k_spr(sprNm),k_spr(sprIICR),"l0,ks of spr/sprIICR"
             endif
             
          enddo
          
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(63,*) (FrmNo-1),"aft Identical Incoming cell props"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(63)
    
  end subroutine making_incoming_cells_identical
  
  
  subroutine making_intrmdteCS_anAverage(strctToRead,cntrlSt,CS1,TopOrBottom,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: CS1
    integer, intent(in)    :: TopOrBottom
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: i,j,m,istrt,imax,mmax
    integer :: cellL,cellR
    integer :: cellL_St,cellL_Nxt
    integer :: cellR_St,cellR_Nxt
    real*8  :: A0_St,A0_Nxt,ka_St,ka_Nxt
    
    integer :: sprNmL,sprNmL_St,sprNmL_Nxt
    integer :: sprNmR,sprNmR_St,sprNmR_Nxt
    real*8  :: l0_St,l0_Nxt,ks_St,ks_Nxt
    
    integer :: sprCntEachCell
    integer :: LastCellAbove,FrstCellAdjstBottm
    integer :: saveWhom
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    real*8, allocatable :: A0_In(:),ka_In(:)
    real*8, allocatable :: l0_In(:),ks_In(:)
    
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
    
    open(unit=43,file='making_intrmdCS.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(43,*) (FrmNo-1),"bfr avrging IntrmdStrct props"
    
    allocate(A0_In(1:N_cell),ka_In(1:N_cell))
    allocate(l0_In(1:N_spr),ks_In(1:N_spr))
    
    A0_In(1:N_cell) = A0_strctNI(CS1,1:N_cell)
    ka_In(1:N_cell) = ka_strctNI(CS1,1:N_cell)
    l0_In(1:N_spr)  = l0_strctNI(CS1,1:N_spr)
    ks_In(1:N_spr)  = ks_strctNI(CS1,1:N_spr)
    
    
    if (TopOrBottom == 1) then
       call find_LastCellofTopLayer(LastCellAbove)
       istrt = 1
       imax  = LastCellAbove !(Hlf_Ncell+1) - 5 !leaving last 5 pair incldin odd pair
       
    elseif (TopOrBottom == 2) then
       call find_FrstCellofBottomEndToAdjust(FrstCellAdjstBottm)
       istrt = FrstCellAdjstBottm
       imax  = Hlf_Ncell
       
    endif
    
    sprCntEachCell = 1 + (NAEC+1)*1 + (NAEC+1)*1 ! Apcl + Bsal + Ltrl
    
    do i = istrt,imax
       cellL = i ; cellR = Hlf_Ncell+i
       
       cellL_St = cellL ; cellL_Nxt = cellL+1
       cellR_St = cellR ; cellR_Nxt = cellR+1
       
       write(43,*) cellL,cellL_St,cellL_Nxt,"cellL's"
       write(43,*) cellR,cellR_St,cellR_Nxt,"cellR's"
       
       do j = 1,2
          
          if (j==1) then
             
             A0_St  = A0_In(cellL_St)  ; ka_St  = ka_In(cellL_St)
             A0_Nxt = A0_In(cellL_Nxt) ; ka_Nxt = ka_In(cellL_Nxt)
             
             A0(cellL)     = (A0_St+A0_Nxt)/2.0d0
             !k_area(cellL) = (ka_St+ka_Nxt)/2.0d0
             
             write(43,*) A0_St,A0_Nxt,A0(cellL),"A0's of L"
             !write(43,*) ka_St,ka_Nxt,k_area(cellL),"k_area's of L"
             
          elseif (j==2) then
             
             A0_St  = A0_In(cellR_St)  ; ka_St  = ka_In(cellR_St)
             A0_Nxt = A0_In(cellR_Nxt) ; ka_Nxt = ka_In(cellR_Nxt)
             
             A0(cellR)     = (A0_St+A0_Nxt)/2.0d0
             !k_area(cellR) = (ka_St+ka_Nxt)/2.0d0
             
             write(43,*) A0_St,A0_Nxt,A0(cellR),"A0's of R"
             !write(43,*) ka_St,ka_Nxt,k_area(cellR),"k_area's of R"
             
          endif
          
          write(43,*) " "
          
          mmax = sprCntEachCell
          
          do m = 1,mmax
             
             if (j==1) then
                
                sprNmL     = sprCntEachCell*(cellL-1) + m 
                sprNmL_St  = sprNmL
                sprNmL_Nxt = sprNmL + sprCntEachCell
                write(43,*) sprNmL,sprNmL_St,sprNmL_Nxt,"sprL's"
                
                l0_St  = l0_In(sprNmL_St)  ; ks_St  = ks_In(sprNmL_St)
                l0_Nxt = l0_In(sprNmL_Nxt) ; ks_Nxt = ks_In(sprNmL_Nxt)
                
                l0(sprNmL)    = (l0_St+l0_Nxt)/2.0d0
                !k_spr(sprNmL) = (ks_St+ks_Nxt)/2.0d0
                
                write(43,*) l0_St,l0_Nxt,l0(sprNmL),"l0's of L"
                !write(43,*) ks_St,ks_Nxt,k_spr(sprNmL),"k_spr's of L"
                
             elseif (j==2) then
                
                sprNmR     = sprCntEachCell*(cellR-1) + m 
                sprNmR_St  = sprNmR
                sprNmR_Nxt = sprNmR + sprCntEachCell
                write(43,*) sprNmR,sprNmR_St,sprNmR_Nxt,"sprR's"
                
                l0_St  = l0_In(sprNmR_St)  ; ks_St  = ks_In(sprNmR_St)
                l0_Nxt = l0_In(sprNmR_Nxt) ; ks_Nxt = ks_In(sprNmR_Nxt)
                
                l0(sprNmR)    = (l0_St+l0_Nxt)/2.0d0
                !k_spr(sprNmR) = (ks_St+ks_Nxt)/2.0d0
                
                write(43,*) l0_St,l0_Nxt,l0(sprNmR),"l0's of R"
                !write(43,*) ks_St,ks_Nxt,k_spr(sprNmR),"k_spr's of R"
                
             endif
             
          enddo
          
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(43,*) (FrmNo-1),"aft avrging IntrmdStrct props"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(43)
    
  end subroutine making_intrmdteCS_anAverage
  


  
  
  subroutine making_intrmdteCS_anOfAvrgCS1CS4(strctToRead,cntrlSt,CS1,CS4,TopOrBottom,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: CS1,CS4
    integer, intent(in)    :: TopOrBottom
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: i,j,m,istrt,imax,mmax

    integer :: cellL,cellR
    real*8  :: A0_St,A0_Fn,ka_St,ka_Fn !St=Start,Fn=Finish
    
    integer :: sprNmL,sprNmR
    real*8  :: l0_St,l0_Fn,ks_St,ks_Fn
    
    integer :: sprCntEachCell
    integer :: LastCellAbove,FrstCellAdjstBottm
    integer :: saveWhom
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    real*8, allocatable :: A0_In(:),ka_In(:) ! In=Input
    real*8, allocatable :: l0_In(:),ks_In(:)
    real*8, allocatable :: A0_Ot(:),ka_Ot(:) ! Ot=Output
    real*8, allocatable :: l0_Ot(:),ks_Ot(:)
    
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
    
    open(unit=43,file='making_intrmdCS_withCS1CS4.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(43,*) (FrmNo-1),"bfr avrging IntrmdStrct props"
    
    allocate(A0_In(1:N_cell),ka_In(1:N_cell))
    allocate(l0_In(1:N_spr) ,ks_In(1:N_spr))
    allocate(A0_Ot(1:N_cell),ka_Ot(1:N_cell))
    allocate(l0_Ot(1:N_spr) ,ks_Ot(1:N_spr))
    
    call read_strctProps_withAddedCell(CS1)
    call read_strctProps_withAddedCell(CS4)
    
    A0_In(1:N_cell) = A0_strctNI(CS1,1:N_cell)
    ka_In(1:N_cell) = ka_strctNI(CS1,1:N_cell)
    l0_In(1:N_spr)  = l0_strctNI(CS1,1:N_spr)
    ks_In(1:N_spr)  = ks_strctNI(CS1,1:N_spr)
    
    A0_Ot(1:N_cell) = A0_strctNI(CS4,1:N_cell)
    ka_Ot(1:N_cell) = ka_strctNI(CS4,1:N_cell)
    l0_Ot(1:N_spr)  = l0_strctNI(CS4,1:N_spr)
    ks_Ot(1:N_spr)  = ks_strctNI(CS4,1:N_spr)
    
    
    if (TopOrBottom == 1) then
       call find_LastCellofTopLayer(LastCellAbove)
       istrt = 1
       imax  = LastCellAbove !(Hlf_Ncell+1) - 5 !leaving last 5 pair incldin odd pair
       
    elseif (TopOrBottom == 2) then
       call find_FrstCellofBottomEndToAdjust(FrstCellAdjstBottm)
       istrt = FrstCellAdjstBottm
       imax  = Hlf_Ncell
    endif
    
    sprCntEachCell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    do i = istrt,imax
       cellL = i ; cellR = Hlf_Ncell+i
       write(43,*) cellL,cellR,"cellL,cellR"
       
       do j = 1,2
          
          if (j==1) then
             
             A0_St = A0_In(cellL) ; ka_St = ka_In(cellL)
             A0_Fn = A0_Ot(cellL) ; ka_Fn = ka_Ot(cellL)
             
             A0(cellL)     = (A0_St+A0_Fn)/2.0d0
             k_area(cellL) = (ka_St+ka_Fn)/2.0d0
             
             write(43,*) A0_St,A0_Fn,A0(cellL),"A0's of L"
             write(43,*) ka_St,ka_Fn,k_area(cellL),"k_area's of L"
             
          elseif (j==2) then
             
             A0_St = A0_In(cellR) ; ka_St = ka_In(cellR)
             A0_Fn = A0_Ot(cellR) ; ka_Fn = ka_Ot(cellR)
             
             A0(cellR)     = (A0_St+A0_Fn)/2.0d0
             k_area(cellR) = (ka_St+ka_Fn)/2.0d0
             
             write(43,*) A0_St,A0_Fn,A0(cellR),"A0's of R"
             write(43,*) ka_St,ka_Fn,k_area(cellR),"k_area's of R"
             
          endif
          
          write(43,*) " "
          
          if (i==istrt .and. TopOrBottom==2) then
             write(43,*) "istrt-Botm Started"
             
             mmax = (NAEC_Ltrl+1)
             
             do m = 1,mmax
                
                if (j==1) then
                   
                   sprNmL = sprCntEachCell*(cellL-2) + (NAEC_apcl+1) + (NAEC_bsal) + m 
                   write(43,*) sprNmL,"sprL's for istrt-Botm"
                   
                   l0_St = l0_In(sprNmL) ; ks_St = ks_In(sprNmL)
                   l0_Fn = l0_Ot(sprNmL) ; ks_Fn = ks_Ot(sprNmL)
                   
                   l0(sprNmL)    = (l0_St+l0_Fn)/2.0d0
                   k_spr(sprNmL) = (ks_St+ks_Fn)/2.0d0
                   
                   write(43,*) l0_St,l0_Fn,l0(sprNmL),"l0's of L for istrt-Botm"
                   write(43,*) ks_St,ks_Fn,k_spr(sprNmL),"k_spr's of L for istrt-Botm"
                   
                elseif (j==2) then
                   
                   sprNmR = sprCntEachCell*(cellR-2) + (NAEC_apcl+1) + (NAEC_bsal) + m 
                   write(43,*) sprNmR,"sprR's for istrt-Botm"
                   
                   l0_St = l0_In(sprNmR) ; ks_St = ks_In(sprNmR)
                   l0_Fn = l0_Ot(sprNmR) ; ks_Fn = ks_Ot(sprNmR)
                   
                   l0(sprNmR)    = (l0_St+l0_Fn)/2.0d0
                   k_spr(sprNmR) = (ks_St+ks_Fn)/2.0d0
                   
                   write(43,*) l0_St,l0_Fn,l0(sprNmR),"l0's of R for istrt-Botm"
                   write(43,*) ks_St,ks_Fn,k_spr(sprNmR),"k_spr's of R for istrt-Botm"
                   
                endif
                
             enddo
             
             write(43,*) "istrt-Botm Finished"
             
          endif
          
          mmax = sprCntEachCell
          
          do m = 1,mmax
             
             if (j==1) then
                
                sprNmL     = sprCntEachCell*(cellL-1) + m 
                write(43,*) sprNmL,"sprL's"
                
                l0_St = l0_In(sprNmL) ; ks_St = ks_In(sprNmL)
                l0_Fn = l0_Ot(sprNmL) ; ks_Fn = ks_Ot(sprNmL)
                
                l0(sprNmL)    = (l0_St+l0_Fn)/2.0d0
                k_spr(sprNmL) = (ks_St+ks_Fn)/2.0d0
                
                write(43,*) l0_St,l0_Fn,l0(sprNmL),"l0's of L"
                write(43,*) ks_St,ks_Fn,k_spr(sprNmL),"k_spr's of L"
                
             elseif (j==2) then
                
                sprNmR     = sprCntEachCell*(cellR-1) + m 
                write(43,*) sprNmR,"sprR's"
                
                l0_St = l0_In(sprNmR) ; ks_St = ks_In(sprNmR)
                l0_Fn = l0_Ot(sprNmR) ; ks_Fn = ks_Ot(sprNmR)
                
                l0(sprNmR)    = (l0_St+l0_Fn)/2.0d0
                k_spr(sprNmR) = (ks_St+ks_Fn)/2.0d0
                
                write(43,*) l0_St,l0_Fn,l0(sprNmR),"l0's of R"
                write(43,*) ks_St,ks_Fn,k_spr(sprNmR),"k_spr's of R"
                
             endif
             
          enddo
          
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(43,*) (FrmNo-1),"aft avrging IntrmdStrct props"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(43)
    
  end subroutine making_intrmdteCS_anOfAvrgCS1CS4
  
  
  
  subroutine getting_weightedAvgpropsFrmCS1andCS4(strctToRead,CS1,CS4,cntrlSt,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead,CS1,CS4,cntrlSt
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer, allocatable :: sprLgcl(:),cellLgcl(:),nodeLgcl(:)!not actual logicl but will perform same
    real*8 , allocatable :: weightspr(:),weightCell(:),weightNode(:)
    real*8 , allocatable :: PresVal(:),TnsnVal(:)
    real*8 , allocatable :: A0I(:),A0F(:),kaI(:),kaF(:)
    real*8 , allocatable :: l0I(:),l0F(:),ksI(:),ksF(:)
    real*8 , allocatable :: CgXI(:),CgXF(:)
    
    integer :: sprCntEachCell,saveWhom
    integer :: cellL,cellR
    integer :: apSprL,apSprR
    integer :: i,j,jmax
    real*8  :: A0_cs1,A0_cs4,l0_cs1,l0_cs4,E
    
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
    
    open(unit=41,file='weightdAvg.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    allocate(sprLgcl(1:N_spr),cellLgcl(1:N_cell),nodeLgcl(1:N_node))
    
    allocate(A0I(1:N_cell),kaI(1:N_cell))
    allocate(l0I(1:N_spr), ksI(1:N_spr))
    allocate(CgXI(1:N_node))
    
    allocate(A0F(1:N_cell),kaF(1:N_cell))
    allocate(l0F(1:N_spr), ksF(1:N_spr))
    allocate(CgXF(1:N_node))
    
    allocate(weightCell(1:N_cell),weightSpr(1:N_spr),weightNode(1:N_node))
    
    call read_strctProps_withAddedCell(CS1)
    
    A0I(1:N_cell)  = A0_strctNI(CS1,1:N_cell)
    kaI(1:N_cell)  = ka_strctNI(CS1,1:N_cell)
    l0I(1:N_spr)   = l0_strctNI(CS1,1:N_spr)
    ksI(1:N_spr)   = ks_strctNI(CS1,1:N_spr)
    CgXI(1:N_node) = CgX_strctNI(CS1,1:N_node)
    
    call read_strctProps_withAddedCell(CS4)
    
    A0F(1:N_cell)  = A0_strctNI(CS4,1:N_cell)
    kaF(1:N_cell)  = ka_strctNI(CS4,1:N_cell)
    l0F(1:N_spr)   = l0_strctNI(CS4,1:N_spr)
    ksF(1:N_spr)   = ks_strctNI(CS4,1:N_spr)
    CgXF(1:N_node) = CgX_strctNI(CS4,1:N_node)
    
    sprLgcl(1:N_spr)   = 0 ; weightSpr(1:N_spr)   = 0.0d0
    cellLgcl(1:N_cell) = 0 ; weightCell(1:N_cell) = 0.0d0
    nodeLgcl(1:N_node) = 0 ; weightNode(1:N_node) = 0.0d0
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(41,*) (FrmNo-1),"bfr weighted avrging IntrmdStrct props"
    
    
    sprCntEachCell = 1 + 1*(NAEC+1) + 1*(NAEC+1)
    
    !!! get lgcls and weighted !!!
    
    cellL = 10 ; cellR = Hlf_Ncell+10
    !apSprL = (cellL-1)*(sprCntEachCell) + 1
    !apSprR = (cellR-1)*(sprCntEachCell) + 1
    
    !sprLgcl(apSprL) = 1 ; sprLgcl(apSprR) = 1
    !weightSpr(apSprL) = 0.40d0 ; weightSpr(apSprR) = 0.40d0
    
    call lgcls_and_weighted_arrays(cellL,cellR,sprLgcl,cellLgcl,weightSpr,weightCell)
    !!! get lgcls and weighted !!!
    
    do i = 1,2 ! will be 1,3 if it is to do for CgXNodes
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          
          if (i==1) then
             
             if (cellLgcl(j)==1) then
                A0_cs1 = A0I(j)
                A0_cs4 = A0F(j)
                A0(j) = weightCell(j)*(A0_cs1+A0_cs4)
                write(41,*) A0(j),A0_cs1,A0_cs4,weightCell(j),j,"A0 and weights"
                
             elseif (cellLgcl(j)==0) then
                continue
             else
                write(*,*) cellLgcl(j),"cellLgcl has different value than 0,1"
                call sleep(10)
             endif
             
          elseif (i==2) then
             
             if (sprLgcl(j)==1) then
                l0_cs1 = l0I(j)
                l0_cs4 = l0F(j)
                l0(j) = weightSpr(j)*(l0_cs1+l0_cs4)
                write(41,*) l0(j),l0_cs1,l0_cs4,weightSpr(j),j,"l0 and weights"
                
             elseif (sprLgcl(j)==0) then
                continue
             else
                write(*,*) sprLgcl(j),"sprLgcl has different value than 0,1"
                call sleep(10)
             endif
             
          endif
          
       enddo
       
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(41,*) (FrmNo-1),"aft weighted avrging IntrmdStrct props"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(41)
    
  end subroutine getting_weightedAvgpropsFrmCS1andCS4
  
  
  subroutine ApprchnCellSprShortenWithSaving(strctToRead,cntrlSt,hmuch,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead
    integer, intent(in)    :: cntrlSt
    real*8 , intent(inout) :: hmuch
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: E
    integer :: LCAL,LCAR
    integer :: saveWhom
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
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

    open(unit=42,file='DiagSprShortenNI.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(42,*) (FrmNo-1),"Frm Bfr DiagShorten"
    
    call find_LastCellofTopLayer(LCAL)
    LCAR = Hlf_Ncell + LCAL
    write(42,*) LCAL,LCAR,"LastCellAbv"
    
    call ApprchnCellSprShortening_and_Equilibrate_UntilMetNI(LCAL,LCAR,hmuch,ExpNo,FrmNo)
    write(42,*) (FrmNo-1),"FrmNo aft DiagSprShorten"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(42)
    
  end subroutine ApprchnCellSprShortenWithSaving
  
  
  subroutine set_lower_bndry_to_AvgPostn_and_Equilibrate(struct_No,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: struct_No
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    integer :: i,j
    integer :: topNode,botNode
    real*8  :: xtop,ytop,xbot,ybot
    
    integer :: hlf_NnodeS
    real*8  :: ZERO,E
    integer :: cntBotNode
    real*8  :: tot_yNodeLft, tot_yNodeRgt
    real*8  :: yNodeBndryLft,yNodeBndryRgt
    integer :: botNodeLft,botNodeRgt
    integer :: BndryLft,BndryRgt
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "stageNo not= 1 and stageType not=1"
    endif
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    call read_strctProps_withAddedCell(struct_No)
    
    A0(1:N_cell)      = A0_strctNI(struct_No,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(struct_No,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(struct_No,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(struct_No,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(struct_No,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    open(unit=45,file='AvgPostnSet.dat')
    
    ZERO = 0.00d0
    hlf_NnodeS = (N_nodeS/2)
    write(45,*) hlf_NnodeS,"hlfNodeS"
    
    do i = 1,hlf_NnodeS
       topNode = (2*i-1)         ; botNode = (2*i)
       xtop = node_xy(topNode,1) ; ytop = node_xy(topNode,2)
       xbot = node_xy(botNode,1) ; ybot = node_xy(botNode,2)
       
       if ((ytop.lt.ZERO) .and. (ybot.lt.ZERO)) then
          cntBotNode = i-1
          write(45,*) cntBotNode,"cntBotNode"
          exit
       endif
       
    enddo
    
    tot_yNodeLft = 0.00d0 ; tot_yNodeRgt = 0.00d0
    
    do i = 2,cntBotNode !strts frm 2 due to avrging from 2nd to last
       botNodeLft = (2*i)
       botNodeRgt = hlf_NnodeS + (2*i)
       
       tot_yNodeLft = tot_yNodeLft + node_xy(botNodeLft,2)
       tot_yNodeRgt = tot_yNodeRgt + node_xy(botNodeRgt,2)
       write(45,*) node_xy(botNodeLft,2),node_xy(botNodeRgt,2),"yBot lft & rgt"
    enddo
    
    yNodeBndryLft = tot_yNodeLft/(cntBotNode-1)
    yNodeBndryRgt = tot_yNodeRgt/(cntBotNode-1)
    write(45,*) yNodeBndryLft,yNodeBndryRgt,"yNodeBndryLft & Rgt"
    
    BndryLft = 2 ; BndryRgt = hlf_NnodeS+2
    
    node_xy(BndryLft,2) = yNodeBndryLft
    node_xy(BndryRgt,2) = yNodeBndryRgt
    
    write(45,*) node_xy(BndryLft,1:2),"node_xyBndryLftNew"
    write(45,*) node_xy(BndryRgt,1:2),"node_xyBndryRgtNew"
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(45,*) N_nodeS,N_node,"N_nodes"
    
    node_xyStr(1:N_nodeS,1:N_dmnsn) = node_xy(1:N_nodeS,1:N_dmnsn) !update node_xyStr
    close(45)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
  end subroutine set_lower_bndry_to_AvgPostn_and_Equilibrate
  
  
  subroutine rescaling_the_system(struct_No,rsFctr,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: struct_No
    real*8 , intent(in)    :: rsFctr
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8, allocatable :: TnsnVal(:),TnsnRscld(:)
    real*8, allocatable :: PresVal(:),PresRscld(:)
    real*8, allocatable :: l0_Rscld(:),A0_Rscld(:)
    
    integer :: i,j,imax,jmax
    real*8  :: E
    integer :: cntrlSt
    integer :: saveWhom
    
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
    
    open(unit=54,file='rescale.dat',position='append')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    call read_strctProps_withAddedCell(struct_No)
    
    A0(1:N_cell)      = A0_strctNI(struct_No,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(struct_No,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(struct_No,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(struct_No,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(struct_No,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(54,*) FrmNo,"FrmNo Aft Read"
    write(54,*) " "
    
    allocate(TnsnVal(1:N_spr),PresVal(1:N_cell))
    allocate(TnsnRscld(1:N_spr),PresRscld(1:N_cell))
    allocate(l0_Rscld(1:N_spr),A0_Rscld(1:N_cell))
    
    TnsnVal   = -1.0d30 ; PresVal    = 1.0d30
    TnsnRscld = -1.0d30 ; PresRscld  = 1.0d30
    l0_Rscld  = -1.0d30 ; A0_Rscld   = -1.0d30
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    imax = 2
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
     
       do j = 1,jmax
          
          if (i==1) then
             TnsnRscld(j) = rsFctr * TnsnVal(j)
             l0_Rscld(j) = (TnsnRscld(j)/k_spr(j)) + l(j)
             write(54,*) TnsnRscld(j),TnsnVal(j),k_spr(j),j,"Tnsn's"
             write(54,*) l0_Rscld(j),l0(j),j,"l0's"
             
          elseif (i==2) then
             PresRscld(j) = rsFctr * PresVal(j)
             A0_Rscld(j) = (PresRscld(j)/k_area(j)) + A(j)
             write(54,*) PresRscld(j),PresVal(j),k_area(j),j,"Pres's"
             write(54,*) A0_Rscld(j),A0(j),j,"A0's"
          endif
          
       enddo
       write(54,*) " "
    enddo
    
    l0(1:N_spr)  = l0_Rscld(1:N_spr)
    A0(1:N_cell) = A0_Rscld(1:N_cell)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(54,*) FrmNo,"FrmNo"
    write(54,*) " "
    
    saveWhom = 3 !for TN or NI or Both selection
    call saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    if (struct_no == 13) then
       cntrlSt = 1
    elseif (struct_no == 16) then
       cntrlSt = 4
    endif
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(54)
    
  end subroutine rescaling_the_system
  
  
  subroutine rescalingSystem_PresTnsnCg(rsFctr)
    implicit none
    real*8, intent(in) :: rsFctr
    
    real*8  :: grdCalcltn(1:N_node,1:N_dmnsn)
    real*8  :: forceActv(1:N_node,1:N_dmnsn)
    
    real*8  :: TnsnVal(1:N_spr),    TnsnRscld(1:N_spr)
    real*8  :: PresVal(1:N_cell),   PresRscld(1:N_cell)
    real*8  :: l0_Rscld(1:N_spr),   A0_Rscld(1:N_cell)
    real*8  :: CgX_Rscld(1:N_node), CgY_Rscld(1:N_node)
    
    integer :: i,imax,j,jmax
    
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
    
    open(unit=514,file='rescaleOnlyTnsnPress.dat')
    
    grdCalcltn = -1.0d30 ; call get_gradient(node_xy,l0,A0,grdCalcltn)
    forceActv  = -grdCalcltn 
    
    TnsnVal   = -1.0d30 ; PresVal    = -1.0d30
    TnsnRscld = -1.0d30 ; PresRscld  = -1.0d30
    l0_Rscld  = -1.0d30 ; A0_Rscld   = -1.0d30
    CgX_Rscld = -1.0d30 ; CgY_Rscld  = -1.0d30
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    imax = 3
    
    do i = 1,imax
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             TnsnRscld(j) = rsFctr * TnsnVal(j)
             l0_Rscld(j) = (TnsnRscld(j)/k_spr(j)) + l(j)
             write(514,*) l0_Rscld(j),l0(j),j,"l0's"
          elseif (i==2) then
             PresRscld(j) = rsFctr * PresVal(j)
             A0_Rscld(j) = (PresRscld(j)/k_area(j)) + A(j)
             write(514,*) A0_Rscld(j),A0(j),j,"A0's"
          elseif (i==3) then
             CgX_Rscld(j) = rsFctr * CgXNode(j) ; CgY_Rscld(j) = rsFctr * CgYNode(j)
             write(514,*) CgX_Rscld(j),CgY_Rscld(j),j,"Cg's"
          endif
          
       enddo
       write(514,*) " "
    enddo
    
    l0(1:N_spr)       = l0_Rscld(1:N_spr)
    A0(1:N_cell)      = A0_Rscld(1:N_cell)
    CgXNode(1:N_node) = CgX_Rscld(1:N_node)
    CgYNode(1:N_node) = CgY_Rscld(1:N_node)
    
    call Equilibrate_only_NI_model
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    close(514)
    
  end subroutine rescalingSystem_PresTnsnCg
  
  
  subroutine pressure_incr_bottomPart(struct_no,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: struct_no
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8, allocatable :: TnsnVal(:),PresVal(:),PressStrt(:) 
    real*8, allocatable :: A0_incrsed(:)
    integer :: cellCnt,cellL,cellR
    real*8  :: PressRead,E
    integer :: cntrlSt,saveWhom
    integer :: i,imax
    
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
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    call read_strctProps_withAddedCell(struct_No)
    
    A0(1:N_cell)      = A0_strctNI(struct_No,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(struct_No,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(struct_No,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(struct_No,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(struct_No,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(54,*) " "
    
    allocate(TnsnVal(1:N_spr))
    allocate(PresVal(1:N_cell),PressStrt(1:N_cell))
    allocate(A0_incrsed(1:N_cell))
    
    open(unit=68,file='St1Press.dat')
    open(unit=67,file='PressIncrBottm3.dat')
    
    imax = Hlf_Ncell+1
    
    do i = 1,imax
       read(68,*) cellCnt,PressRead
       
       if (i.ne.imax) then
          cellL = i ; cellR = Hlf_Ncell+i
          PressStrt(cellL) = PressRead
          PressStrt(cellR) = PressRead
          
          write(67,*) i,cellL,cellR,PressStrt(cellL),PressStrt(cellR),"cell & PressStrt"
          
          if (i.lt.(imax-2)) then
             A0_incrsed(cellL) = A0(cellL)
             A0_incrsed(cellR) = A0(cellR)
             
          elseif (i.ge.(imax-2)) then
             A0_incrsed(cellL) = (PressStrt(cellL)/k_area(cellL)) + A(cellL)
             A0_incrsed(cellR) = (PressStrt(cellR)/k_area(cellR)) + A(cellR)
          endif

          write(67,*) A0(cellL),A0_incrsed(cellL),cellL,"A0s cmpr"
          write(67,*) A0(cellR),A0_incrsed(cellR),cellR,"A0s cmpr"
             
          
       elseif (i==imax) then
          cellL = N_cell
          PressStrt(cellL)  = PressRead
          write(67,*) i,cellL,PressStrt(cellL),"cell & PressStrt"
          A0_incrsed(cellL) = (PressStrt(cellL)/k_area(cellL)) + A(cellL)
          write(67,*) A0(cellL),A0_incrsed(cellL),cellL,"A0s cmpr"
          
       endif
       
    enddo
    
    
    A0(1:N_cell) = A0_incrsed(1:N_cell) ! ******Imprtnt
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(67,*) FrmNo,"FrmNo"
    
    saveWhom = 3 !for TN or NI or Both selection
    call saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    cntrlSt = 2
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(68)
    close(67)
    
  end subroutine pressure_incr_bottomPart
  
  
  subroutine making_sameProp_for_segmntd_Membrne(struct_No)
    implicit none
    integer, intent(in) :: struct_No
    
    integer :: i,j,m
    integer :: cellL,cellR
    integer :: sprCntEachCell
    integer :: sprStrt,sprFnsh,sprCnt
    integer :: NAEC_membrne,NsprInMembrne,midlSpr
    real*8  :: ksMidl,l0Midl
    integer :: saveWhom
    
    open(unit=47,file='makingSamePrpSegmntd.dat',position='append')
    
    sprCntEachCell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    write(47,*) NAEC_Apcl, NAEC_Bsal,NAEC_Ltrl,"NAEC ..."
    
    do i = 1,Hlf_Ncell
       
       cellL = i ; cellR = Hlf_Ncell+i
       write(47,*) cellL,cellR,"cell ..."
       
       do j = 1,2
          
          if (j==1) sprCnt = (cellL-1)*sprCntEachCell
          if (j==2) sprCnt = (cellR-1)*sprCntEachCell
          
          do m = 1,3
             
             if (m==1) NAEC_membrne = NAEC_Apcl
             if (m==2) NAEC_membrne = NAEC_Bsal
             if (m==3) NAEC_membrne = NAEC_Ltrl
             
             NsprInMembrne = NAEC_membrne+1
             write(47,*) m,NAEC_membrne,NsprInMembrne,"NsprInMem"
             
             if (NsprInMembrne.gt.1) then
                
                if (mod(NsprInMembrne,2) == 1) then
                   midlSpr = sprCnt + (NsprInMembrne/2) + 1
                elseif (mod(NsprInMembrne,2) == 0) then
                   midlSpr = sprCnt + (NsprInMembrne/2)
                endif
                
                write(47,*) midlSpr,"midlSpr"
                
                ksMidl = k_spr(midlSpr) ; l0Midl = l0(midlSpr)
                write(47,*) ksMidl,l0Midl,"midlSpr Props"
                
                sprStrt = sprCnt + 1
                sprFnsh = sprCnt + NsprInMembrne
                
                ks_strctNI(struct_No,sprStrt:sprFnsh) = ksMidl
                l0_strctNI(struct_no,sprStrt:sprFnsh) = l0Midl
                
             endif
             
             sprCnt = sprCnt+NsprInMembrne
             
          enddo
          
       enddo
       
    enddo
    
    write(47,*) " "
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    
    close(47)
    
  end subroutine making_sameProp_for_segmntd_Membrne
  
  subroutine CS1_and_CS4_smthTrnsition(CS1,CS2,CS3,CS4)
    implicit none
    integer, intent(in) :: CS1,CS2,CS3,CS4
    
    real*8,  allocatable :: A0I(:),A0F(:),A0P(:) !I=Initial;F=Final
    real*8,  allocatable :: kaI(:),kaF(:),kaP(:) !P=Pulley State
    real*8,  allocatable :: l0I(:),l0F(:),l0P(:)
    real*8,  allocatable :: ksI(:),ksF(:),ksP(:)
    integer, allocatable :: cellLgcl(:)
    integer, allocatable :: sprCurs(:),sprNxts(:)
    
    integer :: i,j,m,jmax,mmax
    real*8  :: E
    
    integer :: cellsToMatch
    integer :: cellCur,cellNxt,cellPrv
    real*8  :: A0Init,kaInit
    real*8  :: A0Finl,kaFinl
    real*8  :: A0InitN,kaInitN
    real*8  :: A0Pul1,kaPul1
    real*8  :: A0Pul2,kaPul2
    
    integer :: sprCntEachCell,sprCntCurCell,sprCntNxtCell
    integer :: sprCur,sprNxt
    real*8  :: l0Init,ksInit
    real*8  :: l0Finl,ksFinl
    real*8  :: l0InitN,ksInitN
    real*8  :: l0Pul1,l0Pul2
    real*8  :: ksPul1,ksPul2
    
    open(unit=37,file='CS1_CS4_smthTrnsitn.dat')
    open(unit=38,file='cellLgclChk.dat')
    
    call switch_to_NI_model
    E = Energy(node_xy,l0,A0)
    call get_gradient(node_xy,l0,A0,gradient)
    
    allocate(ks_strctIF_NI(1:2,1:N_spr) ,l0_strctIF_NI(1:2,1:N_spr))
    allocate(ka_strctIF_NI(1:2,1:N_cell),A0_strctIF_NI(1:2,1:N_cell))
    
    call read_strctProps_withAddedCell(CS1)
    call read_strctProps_withAddedCell(CS2)
    call read_strctProps_withAddedCell(CS3)
    call read_strctProps_withAddedCell(CS4)
    
    allocate(A0I(1:N_cell),A0F(1:N_cell),A0P(1:N_cell))
    allocate(kaI(1:N_cell),kaF(1:N_cell),kaP(1:N_cell))
    allocate(l0I(1:N_spr) ,l0F(1:N_spr) ,l0P(1:N_spr) )
    allocate(ksI(1:N_spr) ,ksF(1:N_spr) ,ksP(1:N_spr) )
    
    A0I(1:N_cell) = A0_strctNI(CS1,1:N_cell)
    A0F(1:N_cell) = A0_strctNI(CS4,1:N_cell)
    A0P(1:N_cell) = A0_strctNI(CS2,1:N_cell)
    
    kaI(1:N_cell) = ka_strctNI(CS1,1:N_cell)
    kaF(1:N_cell) = ka_strctNI(CS4,1:N_cell)
    kaP(1:N_cell) = ka_strctNI(CS2,1:N_cell)
    
    l0I(1:N_spr)  = l0_strctNI(CS1,1:N_spr)
    l0F(1:N_spr)  = l0_strctNI(CS4,1:N_spr)
    l0P(1:N_spr)  = l0_strctNI(CS2,1:N_spr)
    
    ksI(1:N_spr)  = ks_strctNI(CS1,1:N_spr)
    ksF(1:N_spr)  = ks_strctNI(CS4,1:N_spr)
    ksP(1:N_spr)  = ks_strctNI(CS2,1:N_spr)
    
    allocate(cellLgcl(1:N_cell))
    cellsToMatch = 6 !!! IMPORTANT !!!
    call get_lgclsForCellsMatch(cellLgcl)
    
    sprCntEachcell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    do i = 1,Hlf_Ncell
       
       do j = 1,2 !1 for Left, 2 for Right
          
          if (j==1) cellCur = i           ; cellNxt = cellCur+1 ; cellPrv = cellCur-1
          if (j==2) cellCur = Hlf_Ncell+i ; cellNxt = cellCur+1 ; cellPrv = cellCur-1
          
          write(37,*) cellCur,cellNxt,"cellCur & Nxt"
          
          if (cellLgcl(cellCur)==1) then
             write(38,*) cellCur,cellLgcl(cellCur),"cellCur & Lgcl"
          elseif (cellLgcl(cellCur)==0) then
             write(38,*) cellCur,cellLgcl(cellCur),"cellCur & Lgcl, cycle"
             cycle
          endif
             
          A0Init  = A0I(cellCur) ; kaInit  = kaI(cellCur)
          A0Finl  = A0F(cellCur) ; kaFinl  = kaF(cellCur)
          A0InitN = A0I(cellNxt) ; kaInitN = kaI(cellNxt)
          
          write(37,*) A0Init,A0Finl,A0InitN,"A0's bfr"
          write(37,*) kaInit,kaFinl,kaInitN,"ka's bfr"
          
          A0Pul1 = A0P(cellCur) ; kaPul1 = kaP(cellCur)
          A0Pul2 = A0P(cellNxt) ; kaPul2 = kaP(cellNxt)
          
          write(37,*) A0Pul1,A0Pul2,"A0 puls"
          write(37,*) kaPul1,kaPul2,"ka puls"
          
          A0Finl = (A0Pul1+A0Pul2)/2.0d0 !updtd A0Finl
          kaFinl = (kaPul1+kaPul2)/2.0d0 !updtd kaFinl
          
          write(37,*) A0Finl,kaFinl,"A0,ka updtd"
          write(37,*) A0F(cellCur),A0I(cellNxt),kaF(cellCur),kaI(cellNxt),"prop intrfce bfr"
          
          A0F(cellCur) = A0Finl
          kaF(cellCur) = kaFinl
          
          A0I(cellNxt) = A0Finl
          kaI(cellNxt) = kaFinl
          
          write(37,*) A0F(cellCur),A0I(cellNxt),kaF(cellCur),kaI(cellNxt),"prop intrfce aft"
          
          sprCntCurCell = area_spr(cellCur,0)
          sprCntNxtCell = area_spr(cellNxt,0)
          
          mmax = sprCntCurCell
          allocate(sprCurs(1:mmax),sprNxts(1:mmax))
          sprCurs = 0 ; sprNxts = 0
          
          write(37,*) sprCntCurCell,sprCntNxtCell,mmax,"sprCntCurCell-NxtCell,mmax"
          
          call get_sprCur_and_sprNxt(cellCur,cellNxt,mmax,sprCntEachCell,sprCurs,sprNxts)
          
          do m = 1,mmax
                
             !sprCur = (cellCur-1)*(sprCntEachCell) + m
             !sprNxt = (cellNxt-1)*(sprCntEachCell) + m

             sprCur = sprCurs(m) ; sprNxt = sprNxts(m)
             
             write(37,*) sprCur,sprNxt,"sprCur & Nxt"
             
             l0Init  = l0I(sprCur) ; ksInit  = ksI(sprCur)
             l0Finl  = l0F(sprCur) ; ksFinl  = ksF(sprCur)
             l0InitN = l0I(sprNxt) ; ksInitN = ksI(sprNxt)
             
             write(37,*) l0Init,l0Finl,l0InitN,"l0's bfr"
             write(37,*) ksInit,ksFinl,ksInitN,"ks's bfr"
             
             l0Pul1 = l0P(sprCur) ; ksPul1 = ksP(sprCur)
             l0Pul2 = l0P(sprNxt) ; ksPul2 = ksP(sprNxt)
             
             write(37,*) l0Pul1,l0Pul2,"l0 puls"
             write(37,*) ksPul1,ksPul2,"ks puls"
             
             l0Finl = (l0Pul1+l0Pul2)/2.0d0
             ksFinl = (ksPul1+ksPul2)/2.0d0
             
             write(37,*) l0Finl,ksFinl,"l0,ks updtd"
             write(37,*) l0F(sprCur),l0I(sprNxt),ksF(sprCur),ksI(sprNxt),"prop intrfce bfr"
             
             l0F(sprCur) = l0Finl
             ksF(sprCur) = ksFinl
             
             l0I(sprNxt) = l0Finl
             ksI(sprNxt) = ksFinl
             
             write(37,*) l0F(sprCur),l0I(sprNxt),ksF(sprCur),ksI(sprNxt),"prop intrfce aft"
             
          enddo

          deallocate(sprCurs,sprNxts)
          
       enddo
          
    enddo
    
    close(37)
    close(38)
    
    A0_strctIF_NI(1,1:N_cell) = A0I(1:N_cell)
    A0_strctIF_NI(2,1:N_cell) = A0F(1:N_cell)
    
    ka_strctIF_NI(1,1:N_cell) = kaI(1:N_cell)
    ka_strctIF_NI(2,1:N_cell) = kaF(1:N_cell)
    
    l0_strctIF_NI(1,1:N_spr) = l0I(1:N_spr)
    l0_strctIF_NI(2,1:N_spr) = l0F(1:N_spr)
    
    ks_strctIF_NI(1,1:N_spr) = ksI(1:N_spr)
    ks_strctIF_NI(2,1:N_spr) = ksF(1:N_spr)
    
    open(unit=38,file='pulleyAvrgPropCS1.dat')
    open(unit=39,file='pulleyAvrgPropCS4.dat')
    
    do i = 1,2
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             write(38,*) ks_strctIF_NI(1,j),l0_strctIF_NI(1,j)
             write(39,*) ks_strctIF_NI(2,j),l0_strctIF_NI(2,j)
             
          elseif (i==2) then
             write(38,*) ka_strctIF_NI(1,j),A0_strctIF_NI(1,j)
             write(39,*) ka_strctIF_NI(2,j),A0_strctIF_NI(2,j)
          endif
          
       enddo

       write(38,*) " "
       write(39,*) " "
    enddo
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(38)
    close(39)
    
    
  end subroutine CS1_and_CS4_smthTrnsition
  
  
  subroutine pressure_adjstmnt_at_bottomCells(strctToRead,cntrlSt,LastRglrCell,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: strctToRead,cntrlSt
    integer, intent(in)    :: LastRglrCell
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: PresReadLastRglrCell
    real*8  :: TOL_pres,diffPres,ZERO 
    integer :: pairOrSingl,cellNm
    real*8  :: E
    real*8  :: hmuch,hmuchChng
    integer :: hmuchShrtnr
    
    integer :: saveWhom
    real*8  :: init_dP,curr_dP
    integer :: IsigndP,CsigndP,sign_chng
    integer :: LftRghtCopy
    integer :: cellNmRead,cntrlStRead
    
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    real*8, allocatable :: A0_bfrChnging(:),ka_bfrChnging(:)
    real*8, allocatable :: l0_bfrChnging(:),ks_bfrChnging(:)
    
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
    
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    cellNm = LastRglrCell+2 !N_cell
    pairOrSingl = 2
    
    open(unit=24,file='PresReadLastRglrCell.dat')

    if (cntrlSt==3) then
       write(*,*) "Not for cntrlSt=3"
       stop
    endif
    
    do
       read(24,*) cellNmRead,cntrlStRead,PresReadLastRglrCell
       write(*,*) " "
       write(*,*) PresReadLastRglrCell,"PresReadLastRglrCell"
       write(*,*) " "
       
       if (cellNmRead==cellNm .and. cntrlStRead==cntrlSt) then
          write(*,*) cellNmRead,cellNm,cntrlStRead,cntrlSt,"cellNm,cntrlSt"
          exit
       endif
       
    enddo
    close(24)
    
    open(unit=426,file='pressure_adjstmnt_at_bottomCells.dat',position='append')
    open(unit=427,file='pressure_adjstmnt_at_bottomCells1.dat',position='append')
    open(unit=429,file='whyisnotopening.dat',position='append')
    
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    call read_strctProps_withAddedCell(strctToRead)
    
    A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
    k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
    l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
    k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
    CgXNode(1:N_node) = CgX_strctNI(strctToRead,1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(426,*) (FrmNo-1),"Frm Bfr pressure adjustment"
    write(427,*) (FrmNo-1),"Frm Bfr pressure adjustment"
    write(429,*) (FrmNo-1),"FrmNotopening"
    write(*,*)   (FrmNo-1),"Frm Bfr pressure adjustment"
    
    allocate(A0_bfrChnging(1:N_cell),ka_bfrChnging(1:N_cell))
    allocate(l0_bfrChnging(1:N_spr) ,ks_bfrChnging(1:N_spr))
    
    A0_bfrChnging(1:N_cell) = -1.0d30 ; ka_bfrChnging(1:N_cell) = -1.0d30
    l0_bfrChnging(1:N_spr)  = -1.0d30 ; ks_bfrChnging(1:N_spr)  = -1.0d30
    
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    
    TOL_Pres = 0.002d0 ; ZERO = 0.00d0 !! TOL !!
    hmuchShrtnr = 1
    
    init_dP = PresReadLastRglrCell - PresVal(cellNm)
    
    if (init_dP .gt. ZERO) IsigndP = +1 
    if (init_dP .le. ZERO) IsigndP = -1
    
    
    diffPres = PresReadLastRglrCell-PresVal(cellNm)
    write(426,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),cellNm,"pres"
    write(427,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),cellNm,"pres"
    write(429,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),cellNm,"pres"
    write(*,*)   diffPres,PresReadLastRglrCell,PresVal(cellNm),cellNm,"pres"
    
    
    if (abs(diffPres) .gt. TOL_pres) then
       hmuchChng = 0.010d0 ! 1 prcnt
       
       if (diffPres .gt. ZERO) then ! Pres(cellNm) < PresReadLastRgl
          hmuch = 1.00d0 + hmuchChng
       elseif (diffPres .lt. ZERO) then ! Pres(cellNm) > PresReadLastRgl
          hmuch = 1.00d0 - hmuchChng
       endif
       
       do
          
          A0_bfrChnging(1:N_cell) = A0(1:N_cell) ; ka_bfrChnging(1:N_cell) = k_area(1:N_cell)
          l0_bfrChnging(1:N_spr)  = l0(1:N_spr)  ; ks_bfrChnging(1:N_spr)  = k_spr(1:N_spr)
          
          call change_A0of_CellPairOrSinglCell(pairOrSingl,cellNm,hmuch)
          !call change_l0of_CellPairOrSinglCell(pairOrSingl,cellNm,-hmuch) !if A0 reduce,l0 increase (revrse)
          
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          Frame_NI = Frame_NI+1
          write(426,*) (FrmNo-1),"insd do loop"
          write(427,*) (FrmNo-1),"insd do loop"
          write(*,*)   (FrmNo-1),"insd do loop"
          
          PresVal(1:N_cell) = Pressure(node_xy,A0)
          diffPres = PresReadLastRglrCell-PresVal(cellNm)
          
          curr_dP  = diffPres 
          
          if (curr_dP .gt. ZERO) CsigndP = +1 
          if (curr_dP .le. ZERO) CsigndP = -1
          
          sign_chng = IsigndP*CsigndP
          
          write(426,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),"diff & PresRead,PresCell"
          write(426,*) sign_chng,"sign_chng"
          write(427,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),"diff & PresRead,PresCell"
          write(427,*) sign_chng,"sign_chng"
          
          write(*,*) diffPres,PresReadLastRglrCell,PresVal(cellNm),"diff & PresRead,PresCell"
          write(*,*) sign_chng,"sign_chng"
          
          
          if (abs(diffPres) .le. TOL_pres) then
             
             !now I have to wirte where to put the bar
             
             if (diffPres.ge.ZERO) then !
                write(426,*) (FrmNo-1),"Before Exit"
                write(427,*) (FrmNo-1),"Before Exit"
                write(429,*)(FrmNo-1),"Before Exit"
                write(*,*)  (FrmNo-1),"Before Exit"
                exit
                
             elseif (diffPres.lt.ZERO) then
                
                if (hmuchShrtnr == 1) then
                   hmuchChng   = hmuchChng/3.0d0
                   hmuchShrtnr = hmuchShrtnr+1
                   
                elseif (hmuchShrtnr == 2) then
                   continue
                endif
                
                write(426,*) hmuchChng,"hmuchChng with lesserhmuch"
                write(427,*) hmuchChng,"hmuchChng with lesserhmuch"
                write(*,*)  hmuchChng,"hmuchChng with lesserhmuch"
                
                if (diffPres .gt. ZERO) then
                   hmuch = 1.00d0 + hmuchChng
                elseif (diffPres .lt. ZERO) then
                   hmuch = 1.00d0 - hmuchChng
                endif
                
                write(426,*) hmuch,"hmuch with lesserhmuch"
                write(427,*) hmuch,"hmuch with lesserhmuch"
                write(*,*)  hmuch,"hmuch with lesserhmuch"
                
             endif
             
          elseif (sign_chng==(-1)) then
             A0(1:N_cell) = A0_bfrChnging(1:N_cell) ; k_area(1:N_cell) = ka_bfrChnging(1:N_cell)
             l0(1:N_spr)  = l0_bfrChnging(1:N_spr)  ; k_spr(1:N_spr)   = ks_bfrChnging(1:N_spr)
             
             call Equilibrate_system
             call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
             Frame_NI = Frame_NI+1
             write(426,*) (FrmNo-1),"going to restored values"
             write(427,*) (FrmNo-1),"going to restored values"
             write(*,*) (FrmNo-1),"going to restored values"
             call sleep(2)
             
             PresVal(1:N_cell) = Pressure(node_xy,A0)
             diffPres = PresReadLastRglrCell-PresVal(cellNm)
             
             hmuchChng = hmuchChng/2.0d0
             write(426,*) hmuchChng,"hmuchChng"
             write(427,*) hmuchChng,"hmuchChng"
             write(*,*) hmuchChng,"hmuchChng"
             
             
             if (diffPres .gt. ZERO) then
                hmuch = 1.00d0 + hmuchChng
             elseif (diffPres .lt. ZERO) then
                hmuch = 1.00d0 - hmuchChng
             endif
             
             write(426,*) hmuch,"hmuch"
             write(427,*) hmuch,"hmuch"
             write(*,*) hmuch,"hmuch"
          endif
          
          write(426,*) " "
          write(427,*) " "
          write(*,*) " "
          
       enddo
       
    elseif (abs(diffPres) .le. TOL_pres) then
       continue
    endif
    
    !if (cntrlSt==2) then
    !   LftRghtCopy = 1
    !   call symmetricize_node_for_NIsystem(LftRghtCopy)
    !   call nodes_to_coordntes(node_xy,coordntes_xy)
    !   call Equilibrate_system
    !   call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    !   Frame_NI = Frame_NI+1
    !   write(329,*) (Frame_NI-1),"Frm At Symm"
    !endif
        
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(strctToRead,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    call deallocate_repetitive_arrays
    call switchback_to_TN_model
    
    close(326)
    close(327)
    close(329)
    
    !if (cntrlSt==2) stop
    
  end subroutine pressure_adjstmnt_at_bottomCells
  
  
  subroutine making_PropsOneCellShftd_INCS1andCS4(cntrlSt,CS1,CS4,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cntrlSt
    integer, intent(in)    :: CS1,CS4
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8, allocatable :: A0InCS1(:),kaInCS1(:)
    real*8, allocatable :: l0InCS1(:),ksInCS1(:)
    real*8, allocatable :: CgXVInCS1(:)
    real*8, allocatable :: PresVal(:),TnsnVal(:)
    
    integer :: rgn,i,istrt,ifnsh,j
    integer :: cellCurr,cellToBe
    integer :: sprCurr,sprToBe
    
    integer :: saveWhom,nsprs_InEachCell
    real*8  :: E
    
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
    
    open(unit=84,file='making_PropsOneCellShftd_InCS1andCS4.dat')
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
       
    elseif (modelID==2) then
       continue
    endif
    
    allocate(kaInCS1(1:N_cell),A0InCS1(1:N_cell))
    allocate(ksInCS1(1:N_spr) ,l0InCS1(1:N_spr))
    allocate(CgXVInCS1(1:N_node))
    allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
    
    call read_strctProps_withAddedCell(CS1)
    
    kaInCS1(1:N_cell)   = ka_strctNI(CS1,1:N_cell)
    A0InCS1(1:N_cell)   = A0_strctNI(CS1,1:N_cell)
    ksInCS1(1:N_spr)    = ks_strctNI(CS1,1:N_spr)
    l0InCS1(1:N_spr)    = l0_strctNI(CS1,1:N_spr)
    CgXVInCS1(1:N_node) = CgX_strctNI(CS1,1:N_node)

    nsprs_InEachCell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    write(84,*) nsprs_InEachCell,"nsprs_InEachCell"
    
    do rgn = 1,2
       
       if (rgn==1) then
          istrt = 1
          ifnsh = Hlf_Ncell-2
          
       elseif (rgn==2) then
          istrt = Hlf_Ncell+1 
          ifnsh = Hlf_Ncell+Hlf_Ncell-2
       endif
       
       do i = 1,(Hlf_Ncell-2)
          
          cellCurr = i ; cellToBe = i+1
          
          k_area(cellCurr) = kaInCS1(cellToBe)
          A0(cellCurr)     = A0InCS1(cellToBe)
          
          write(84,*) cellCurr,cellToBe,"cellCuur,ToBe"
          
          do j = 1,nsprs_InEachCell
             
             sprCurr = (cellCurr-1)*(nsprs_InEachCell) + j
             sprToBe = (cellToBe-1)*(nsprs_InEachCell) + j
             
             k_spr(sprCurr) = ksInCS1(sprToBe)
             l0(sprCurr)    = l0InCS1(sprToBe)
             
             write(84,*) sprCurr,sprToBe,"sprCurr,ToBe"
          enddo
          
          write(84,*) " "
          
       enddo
       
       
    enddo
    
    CgXNode(1:N_node) = CgXVInCS1(1:N_node)
    CgYNode(1:N_node) = 0.0d0
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(84,*) (FrmNo-1),"Aft making Props Shifted"
    
    saveWhom = 3
    call saveA0l0_ofStrct_AddedCell_TNorNI(CS4,saveWhom)
    
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    
    call writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PresVal,TnsnVal)
    
    if (modelID==1) then
       continue
    elseif (modelID==2) then
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
    endif
    
    close(84)
    
  end subroutine making_PropsOneCellShftd_INCS1andCS4
  
  
  ! subroutine pulling_or_pushing_forces_at_CS(strctToRead,cntrlSt,ExpNo,FrmNo)
  !   implicit none
  !   integer, intent(in)    :: strctToRead,cntrlSt,ExpNo
  !   integer, intent(inout) :: FrmNo
    
  !   integer :: RghtstrtTN
  !   real*8  :: E
    
  !   real*8, allocatable :: PresVal(:),TnsnVal(:)
    
  !   interface
       
  !      function TnsnComprsn(dum_Nodes,duml0)
  !        use system_parameters
  !        use transfrm_info
  !        implicit none
  !        real*8 :: TnsnComprsn(1:N_spr)
  !        real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
  !        real*8 :: duml0(1:N_spr)
  !      end function TnsnComprsn
       
  !      function Pressure(dum_Nodes,dumA0)
  !        use system_parameters
  !        use transfrm_info
  !        implicit none
  !        real*8 :: Pressure(1:N_cell)
  !        real*8 :: dum_Nodes(1:N_node,1:N_dmnsn)
  !        real*8 :: dumA0(1:N_cell)
  !      end function Pressure
       
  !   end interface
    
  !   open(unit=64,file='PullingPushingAtNI.dat')
    
  !   if (modelID==1) then
  !      call switch_to_NI_model
  !      E = Energy(node_xy,l0,A0)
  !      call get_gradient(node_xy,l0,A0,gradient)
  !   elseif (modelID==2) then
  !      continue
  !   endif
    
  !   allocate(PresVal(1:N_cell),TnsnVal(1:N_spr))
  !   call read_strctProps_withAddedCell(strctToRead)
    
  !   A0(1:N_cell)      = A0_strctNI(strctToRead,1:N_cell)
  !   k_area(1:N_cell)  = ka_strctNI(strctToRead,1:N_cell)
  !   l0(1:N_spr)       = l0_strctNI(strctToRead,1:N_spr)
  !   k_spr(1:N_spr)    = ks_strctNI(strctToRead,1:N_spr)
  !   CgXNode(1:N_spr)  = CgX_strctNI(strctToRead,1:N_node)
  !   CgYNode(1:N_node) = 0.0d0
  
  !   call Equilibrate_system
  !   call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
  !   FrmNo = FrmNo+1
    
  !   RghtstrtTN = (Hlf_Ncell+1)*2 + 1
  !   CgXNode(1:2) = 0.15d0 ; CgXNode(RghtstrtTN:(RghtstrtTN+1)) = -0.15d0
    
  !   write(*,*) node_xy(1,1:2),node_xy(2,1:2),"node_xy 1-2"
  !   write(*,*) node_xy(RghtstrtTN,1:2),node_xy(RghtstrtTN+1,1:2),"node_xy R 1-2"
    
    
  !   if (cntrlSt == 1) then
       
  !      write(*,*) cntrlSt,(Frame_NI-1),"cntrlSt1,Frame_NI"
  !      call Equilibrate_system
  !      call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,Frame_NI)
  !      Frame_NI = Frame_NI+1
       
  !   elseif (cntrlSt == 2) then
       
  !      write(*,*) cntrlSt,(Frame_NI-1),"cntrlSt2,Frame_NI"
  !      call Equilibrate_system
  !      call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,Frame_NI)
  !      Frame_NI = Frame_NI+1
       
  !      allocate(A0_I(1:N_cell),ka_I(1:N_cell),l0_I(1:N_spr),ks_I(1:N_spr))
  !      allocate(A0_F(1:N_cell),ka_F(1:N_cell),l0_F(1:N_spr),ks_F(1:N_spr))
       
  !      A0_I(1:N_cell)=1.0d30 ; ka_I(1:N_cell)=1.0d30 ; l0_I(1:N_spr)=1.0d30 ; ks_I(1:N_spr)=1.0d30
  !      A0_F(1:N_cell)=1.0d30 ; ka_F(1:N_cell)=1.0d30 ; l0_F(1:N_spr)=1.0d30 ; ks_F(1:N_spr)=1.0d30
       
  !      allocate(A0I_NI(1:N_cell),kaI_NI(1:N_cell),l0I_NI(1:N_spr),ksI_NI(1:N_spr))
  !      allocate(A0F_NI(1:N_cell),kaF_NI(1:N_cell),l0F_NI(1:N_spr),ksF_NI(1:N_spr))
       
  !      A0I_NI(1:N_cell)=1.0d30 ; kaI_NI(1:N_cell)=1.0d30 ; l0I_NI(1:N_spr)=1.0d30 ; ksI_NI(1:N_spr)=1.0d30
  !      A0F_NI(1:N_cell)=1.0d30 ; kaF_NI(1:N_cell)=1.0d30 ; l0F_NI(1:N_spr)=1.0d30 ; ksF_NI(1:N_spr)=1.0d30
       
  !      A0_F(1:N_cell) = A0_strctNI(strctToRead,1:N_cell) ; A0_I(1:N_cell) = A0_strctNI((strctToRead-1),1:N_cell)
  !      ka_F(1:N_cell) = ka_strctNI(strctToRead,1:N_cell) ; ka_I(1:N_cell) = ka_strctNI((strctToRead-1),1:N_cell)
  !      l0_F(1:N_spr)  = l0_strctNI(strctToRead,1:N_spr)  ; l0_I(1:N_spr)  = l0_strctNI((strctToRead-1),1:N_spr)
  !      ks_F(1:N_spr)  = ks_strctNI(strctToRead,1:N_spr)  ; ks_I(1:N_spr)  = ks_strctNI((strctToRead-1),1:N_spr)
       
  !      write(*,*) A0_F(5),A0_I(5),ka_F(5),ka_I(5),"A0F+kaF in cell 5"
  !      write(*,*) l0_F(9),l0_I(9),ks_F(9),ks_I(9),"l0F+ksF in spr  9"
  
  !      call sleep(10)
  !      call Extrapolate_and_adjust_strctProp_ifNeeded_in_MAE_woSwitchBack(ExpNo,Frame_NI)
       
  !      A0I_NI(1:N_cell) = A0_I(1:N_cell) ; kaI_NI(1:N_cell) = ka_I(1:N_cell)
  !      l0I_NI(1:N_spr)  = l0_I(1:N_spr)  ; ksI_NI(1:N_spr)  = ks_I(1:N_spr)
       
  !      A0F_NI(1:N_cell) = A0_F(1:N_cell) ; kaF_NI(1:N_cell) = ka_F(1:N_cell)
  !      l0F_NI(1:N_spr)  = l0_F(1:N_spr)  ; ksF_NI(1:N_spr)  = ks_F(1:N_spr)

  !      call Equilibrate_system
  !      call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,Frame_NI)
  !      Frame_NI = Frame_NI+1
       
  !   elseif (cntrlSt==3) then
       
  !      write(*,*) cntrlSt,(Frame_NI-1),"cntrlSt3,Frame_NI"
       
  !      A0(1:N_cell) = A0F_NI(1:N_cell)   ; k_area(1:N_cell) = kaF_NI(1:N_cell)
  !      l0(1:N_spr)  = l0F_NI(1:N_spr)    ; k_spr(1:N_spr)   = ksF_NI(1:N_spr)
       
  !      call Equilibrate_system
  !      call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,Frame_NI)
  !      Frame_NI = Frame_NI+1
       
  !   elseif (cntrlSt==4) then
       
  !      write(*,*) cntrlSt,(Frame_NI-1),"cntrlSt4,Frame_NI"
       
  !      call Equilibrate_system
  !      call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,Frame_NI)
  !      Frame_NI = Frame_NI+1
       
  !   endif
    
  !   call switchto_NI_model_run_and_switchbackto_TN
    
    
  !   close(64)
    
  ! end subroutine pulling_or_pushing_forces_at_CS
  
  subroutine pulling_or_pushing(pullingOrpushing,hmuch)
    implicit none
    integer, intent(in) :: pullingOrpushing
    real*8 , intent(in) :: hmuch
    integer :: lftBndryNodes(1:2),rgtBndryNodes(1:2)
    integer :: lftNodes
    
    lftNodes = (Hlf_Ncell+1)*2
    
    lftBndryNodes(1) = 1 ; rgtBndryNodes(1) = lftNodes+1
    lftBndryNodes(2) = 2 ; rgtBndryNodes(2) = lftNodes+2
    
    write(*,*) lftBndryNodes(1:2),"lft bndry"
    write(*,*) rgtBndryNodes(1:2),"rgt bndry"
    
    if (pullingOrpushing == 1) then !pulling
       
       CgXNode(lftBndryNodes(1)) = +hmuch ; CgXNode(lftBndryNodes(2)) = +hmuch
       CgXNode(rgtBndryNodes(1)) = -hmuch ; CgXNode(rgtBndryNodes(2)) = -hmuch
       
    elseif (pullingOrpushing == 2) then !pushing
       
       CgXNode(lftBndryNodes(1)) = -hmuch ; CgXNode(lftBndryNodes(2)) = -hmuch
       CgXNode(rgtBndryNodes(1)) = +hmuch ; CgXNode(rgtBndryNodes(2)) = +hmuch
       
    endif
    
  end subroutine pulling_or_pushing
  

  subroutine pulling_or_pushing_frmCurrntVal(pullingOrpushing,hmuch)
    implicit none
    integer, intent(in) :: pullingOrpushing
    real*8 , intent(in) :: hmuch
    integer             :: lftBndryNodes(1:2),rgtBndryNodes(1:2)
    integer             :: lftNodes
    
    lftNodes = (Hlf_Ncell+1)*2
    
    lftBndryNodes(1) = 1 ; rgtBndryNodes(1) = lftNodes+1
    lftBndryNodes(2) = 2 ; rgtBndryNodes(2) = lftNodes+2
    
    write(*,*) lftBndryNodes(1:2),"lft bndry"
    write(*,*) rgtBndryNodes(1:2),"rgt bndry"
    
    if (pullingOrpushing == 1) then !pulling
       
       CgXNode(lftBndryNodes(1)) = CgXNode(lftBndryNodes(1))+hmuch
       CgXNode(lftBndryNodes(2)) = CgXNode(lftBndryNodes(2))+hmuch
       CgXNode(rgtBndryNodes(1)) = CgXNode(rgtBndryNodes(1))-hmuch
       CgXNode(rgtBndryNodes(2)) = CgXNode(rgtBndryNodes(2))-hmuch
       
    elseif (pullingOrpushing == 2) then !pushing
       
       CgXNode(lftBndryNodes(1)) = CgXNode(lftBndryNodes(1))-hmuch
       CgXNode(lftBndryNodes(2)) = CgXNode(lftBndryNodes(2))-hmuch
       CgXNode(rgtBndryNodes(1)) = CgXNode(rgtBndryNodes(1))+hmuch
       CgXNode(rgtBndryNodes(2)) = CgXNode(rgtBndryNodes(2))+hmuch
       
    endif
    
  end subroutine pulling_or_pushing_frmCurrntVal
  
  subroutine pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa()
    implicit none
    integer :: Dirctn_dcsnInit,Dirctn_dcsnCurr
    integer :: lftBndryNodes(1:2),rgtBndryNodes(1:2)
    integer :: lftNodes
    real*8  :: hmuch
    integer :: i,j,jmax,cntLp
    real*8  :: X_AbsValInitTop,X_AbsValInitBot
    real*8  :: X_AbsValCurrTop,X_AbsValCurrBot
    real*8  :: CgXNodeLT,CgXNodeLB,CgXNodeRT,CgXNodeRB
    
    lftNodes = (Hlf_Ncell+1)*2
    
    lftBndryNodes(1) = 1 ; rgtBndryNodes(1) = lftNodes+1
    lftBndryNodes(2) = 2 ; rgtBndryNodes(2) = lftNodes+2
    
    write(*,*) lftBndryNodes(1:2),"lft bndry" ; write(*,*) rgtBndryNodes(1:2),"rgt bndry"
    
    X_AbsValInitTop = abs(node_xy(lftBndryNodes(1),1))
    X_AbsValInitBot = abs(node_xy(lftBndryNodes(2),1))
    
    if (X_AbsValInitTop .gt. X_AbsValInitBot) Dirctn_dcsnCurr = -1  
    if (X_AbsValInitTop .lt. X_AbsValInitBot) Dirctn_dcsnCurr = +1
    
    Dirctn_dcsnInit = Dirctn_dcsnCurr 
    
    write(*,*) X_AbsValInitTop,X_AbsValInitBot,Dirctn_dcsnInit,"(xt,xb)-abs,Dirctn_dcsnInit"
    
    hmuch = 0.02d0
    cntLp = 0
    
    do
       
       CgXNodeLT = CgXNode(lftBndryNodes(1)) ; CgXNodeLB = CgXNode(lftBndryNodes(2))
       CgXNodeRT = CgXNode(rgtBndryNodes(1)) ; CgXNodeRB = CgXNode(rgtBndryNodes(2))
       
       if (Dirctn_dcsnCurr == +1) then ! TopBndryPULLING_and_BtmBndryPUSHING
          
          CgXNode(lftBndryNodes(1)) = CgXNode(lftBndryNodes(1)) + hmuch
          CgXNode(rgtBndryNodes(1)) = CgXNode(rgtBndryNodes(1)) - hmuch
          
          CgXNode(lftBndryNodes(2)) = CgXNode(lftBndryNodes(2)) - hmuch
          CgXNode(rgtBndryNodes(2)) = CgXNode(rgtBndryNodes(2)) + hmuch
          
       elseif (Dirctn_dcsnCurr == -1) then ! TopBndryPUSHING_and_BtmBndryPULLING
          
          CgXNode(lftBndryNodes(1)) = CgXNode(lftBndryNodes(1)) - hmuch
          CgXNode(rgtBndryNodes(1)) = CgXNode(rgtBndryNodes(1)) + hmuch
          
          CgXNode(lftBndryNodes(2)) = CgXNode(lftBndryNodes(2)) + hmuch
          CgXNode(rgtBndryNodes(2)) = CgXNode(rgtBndryNodes(2)) - hmuch
          
       endif
       
       if (VF_regionModelled==0) call Equilibrate_only_NI_model
       if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
       
       cntLp = cntLp+1
       
       X_AbsValCurrTop = abs(node_xy(lftBndryNodes(1),1))
       X_AbsValCurrBot = abs(node_xy(lftBndryNodes(2),1))
       
       if (X_AbsValCurrTop .gt. X_AbsValCurrBot) Dirctn_dcsnCurr = -1  
       if (X_AbsValCurrTop .lt. X_AbsValCurrBot) Dirctn_dcsnCurr = +1
       
       if (Dirctn_dcsnCurr == Dirctn_dcsnInit) then
          continue
       elseif (Dirctn_dcsnCurr .ne. Dirctn_dcsnInit) then
          
          CgXNode(lftBndryNodes(1)) = 0.5000d0 * (CgXNode(lftBndryNodes(1)) + CgXNodeLT)
          CgXNode(lftBndryNodes(2)) = 0.5000d0 * (CgXNode(lftBndryNodes(2)) + CgXNodeLB)
          CgXNode(rgtBndryNodes(1)) = 0.5000d0 * (CgXNode(rgtBndryNodes(1)) + CgXNodeRT)
          CgXNode(rgtBndryNodes(2)) = 0.5000d0 * (CgXNode(rgtBndryNodes(2)) + CgXNodeRB)
          
          if (VF_regionModelled==0) call Equilibrate_only_NI_model
          if (VF_regionModelled==1) call Equilibrate_only_NI_model_withVF_region
          
          cntLp = cntLp-1
          write(*,*) cntLp,"cntLp Val"
          exit
          
       endif
       
    enddo
    
    if (VF_regionModelled==0) write(*,*) "at time of exit from pullingToppushingBotm, Frame_NI =",  (Frame_NI-1)
    if (VF_regionModelled==1) write(*,*) "at time of exit from pullingToppushingBotm, Frame_NIVF =",(Frame_NIVF-1)
    
  end subroutine pullingTopBndry_or_pushingBotm_slowly_Or_viceVrsa
  
  
  
  subroutine read_Init_Finl_props_from_propsFlnm(flnmInit,flnmFinl,A0_F,A0_I,ka_F,ka_I,l0_F,l0_I,ks_F,ks_I)
    implicit none
    integer, intent(in)  :: flnmInit,flnmFinl
    real*8 , intent(out) :: A0_F(1:N_cell),A0_I(1:N_cell)
    real*8 , intent(out) :: ka_I(1:N_cell),ka_F(1:N_cell)
    real*8 , intent(out) :: l0_F(1:N_spr) ,l0_I(1:N_spr)
    real*8 , intent(out) :: ks_I(1:N_spr) ,ks_F(1:N_spr)
    
    character(len=100)  :: prfxName,flnm1,flnm2,full_flnm1,full_flnm2
    integer             :: i,j,jmax,sprNm,cellNm
    
    prfxName='propsflGen'
    
    write(flnm1,'(i2.2,a)') flnmInit,'.dat'
    write(flnm2,'(i2.2,a)') flnmFinl,'.dat'
    
    write(*,*) trim(adjustl(prfxName)),trim(adjustl(flnm1))
    write(*,*) trim(adjustl(prfxName)),trim(adjustl(flnm2))
    
    full_flnm1=trim(adjustl(prfxName))//trim(adjustl(flnm1))
    full_flnm2=trim(adjustl(prfxName))//trim(adjustl(flnm2))
    
    open(unit=411,file=trim(adjustl(full_flnm1))) 
    open(unit=412,file=trim(adjustl(full_flnm2)))
    
    do i = 1,2
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       
       do j = 1,jmax
          
          if (i==1) then
             read(411,*) ks_I(j),l0_I(j),sprNm
             read(412,*) ks_F(j),l0_F(j),sprNm
             
          elseif (i==2) then
             read(411,*) ka_I(j),A0_I(j),cellNm
             read(412,*) ka_F(j),A0_F(j),cellNm
          endif
             
       enddo
       
    enddo
    
    close(411)
    close(412)
    
  end subroutine read_Init_Finl_props_from_propsFlnm
  
  subroutine save_ProgrsnStgProp(ProgCycl,strctV,strghtOrNI,A0V,l0V,kaV,ksV)
    implicit none
    integer, intent(in) :: ProgCycl
    integer, intent(in) :: strctV
    integer, intent(in) :: strghtOrNI
    real*8 , intent(in) :: A0V(1:N_cell),l0V(1:N_spr)
    real*8 , intent(in) :: kaV(1:N_cell),ksV(1:N_spr)
    
    open(unit=94,file='save_ProgrsnStgProp.dat',position='append')
    write(94,*) ProgCycl,strctV,strghtOrNI,"ProgCycl,strctV,strghtOrNI"
    close(94)
    
    if (strctV==1) then
       
       if (strghtOrNI==1) then ! 1 is for straight spr model
          A0_prgrsnStrct1(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnStrct1(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnStrct1(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnStrct1(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnStrct1(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnStrct1(ProgCycl,1:N_node) = 0.0d0
          
       elseif (strghtOrNI==2) then ! 2 is for curvy spr model
          A0_prgrsnNIstrct1(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnNIstrct1(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnNIstrct1(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnNIstrct1(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnNIstrct1(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnNIstrct1(ProgCycl,1:N_node) = 0.0d0
          
       endif
          
    elseif (strctV==2) then
       
       if (strghtOrNI==1) then
          A0_prgrsnStrct2(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnStrct2(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnStrct2(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnStrct2(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnStrct2(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnStrct2(ProgCycl,1:N_node) = 0.0d0
          
       elseif (strghtOrNI==2) then
          A0_prgrsnNIstrct2(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnNIstrct2(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnNIstrct2(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnNIstrct2(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnNIstrct2(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnNIstrct2(ProgCycl,1:N_node) = 0.0d0
          
       endif
       
    elseif (strctV==3) then
       
       if (strghtOrNI==1) then
          A0_prgrsnStrct3(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnStrct3(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnStrct3(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnStrct3(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnStrct3(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnStrct3(ProgCycl,1:N_node) = 0.0d0
          
       elseif (strghtOrNI==2) then
          A0_prgrsnNIstrct3(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnNIstrct3(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnNIstrct3(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnNIstrct3(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnNIstrct3(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnNIstrct3(ProgCycl,1:N_node) = 0.0d0
          
       endif
       
    elseif (strctV==4) then
       
       if (strghtOrNI==1) then
          A0_prgrsnStrct4(ProgCycl,1:N_cell)  = A0V(1:N_cell)
          l0_prgrsnStrct4(ProgCycl,1:N_spr)   = l0V(1:N_spr)
          ka_prgrsnStrct4(ProgCycl,1:N_cell)  = kaV(1:N_cell)
          ks_prgrsnStrct4(ProgCycl,1:N_spr)   = ksV(1:N_spr)
          CgX_prgrsnStrct4(ProgCycl,1:N_node) = CgXNode(1:N_node)
          CgY_prgrsnStrct4(ProgCycl,1:N_node) = 0.0d0
          
       elseif (strghtOrNI==2) then
          A0_prgrsnNIstrct4(ProgCycl,1:N_cell)   = A0V(1:N_cell)
          l0_prgrsnNIstrct4(ProgCycl,1:N_spr)    = l0V(1:N_spr)
          ka_prgrsnNIstrct4(ProgCycl,1:N_cell)   = kaV(1:N_cell)
          ks_prgrsnNIstrct4(ProgCycl,1:N_spr)    = ksV(1:N_spr)
          CgX_prgrsnNIstrct4(ProgCycl,1:N_node)  = CgXNode(1:N_node)
          CgY_prgrsnNIstrct4(ProgCycl,1:N_node)  = 0.0d0
          
       endif
       
    endif
    
    call writeProgrsnProps(ProgCycl,strctV,strghtOrNI)
    
  end subroutine save_ProgrsnStgProp
  
  
  subroutine writeProgrsnProps(ProgCyl,strctV,strghtOrNI)
    implicit none
    integer, intent(in) :: ProgCyl
    integer, intent(in) :: strctV
    integer, intent(in) :: strghtOrNI
    
    character(len=100)  :: wflnm
    character(len=100)  :: wflnmbr
    character(len=100)  :: wfull_flnm
    
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 3
    
    if (strghtOrNI==1) then
       wflnm = 'PrgrsnCyclProp'
    elseif (strghtOrNI==2) then
       wflnm = 'PrgrsnCyclPropNI'
    endif
    
    write(wflnmbr,'(i1.1,a,i1.1,a)') ProgCyl,'Strct',strctV,'.dat'
    
    wfull_flnm=trim(adjustl(wflnm))//trim(adjustl(wflnmbr))
    
    open(unit=32,file=trim(adjustl(wfull_flnm)))
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (strctV==1) then
             
             if (i==1) then
                
                if (strghtOrNI==1) write(32,*) ks_prgrsnStrct1(ProgCyl,j)  , l0_prgrsnStrct1(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ks_prgrsnNIstrct1(ProgCyl,j), l0_prgrsnNIstrct1(ProgCyl,j)
                
             elseif (i==2) then
                
                if (strghtOrNI==1) write(32,*) ka_prgrsnStrct1(ProgCyl,j)  , A0_prgrsnStrct1(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ka_prgrsnNIstrct1(ProgCyl,j), A0_prgrsnNIstrct1(ProgCyl,j)
                
             elseif (i==3) then
                
                if (strghtOrNI==1) write(32,*) CgX_prgrsnStrct1(ProgCyl,j)   ,CgY_prgrsnStrct1(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) CgX_prgrsnNIstrct1(ProgCyl,j) ,CgY_prgrsnNIstrct1(ProgCyl,j)
                
             endif
                
          elseif (strctV==2) then
             
             if (i==1) then
                
                if (strghtOrNI==1) write(32,*) ks_prgrsnStrct2(ProgCyl,j)  ,l0_prgrsnStrct2(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ks_prgrsnNIstrct2(ProgCyl,j),l0_prgrsnNIstrct2(ProgCyl,j)
                
             elseif (i==2) then
                
                if (strghtOrNI==1) write(32,*) ka_prgrsnStrct2(ProgCyl,j)  ,A0_prgrsnStrct2(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ka_prgrsnNIstrct2(ProgCyl,j),A0_prgrsnNIstrct2(ProgCyl,j)
                
             elseif (i==3) then
                
                if (strghtOrNI==1) write(32,*) CgX_prgrsnStrct2(ProgCyl,j)   ,CgY_prgrsnStrct2(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) CgX_prgrsnNIstrct2(ProgCyl,j) ,CgY_prgrsnNIstrct2(ProgCyl,j)
                
             endif
             
          elseif (strctV==3) then
             
             if (i==1) then
                
                if (strghtOrNI==1) write(32,*) ks_prgrsnStrct3(ProgCyl,j)  ,l0_prgrsnStrct3(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ks_prgrsnNIstrct3(ProgCyl,j),l0_prgrsnNIstrct3(ProgCyl,j)
                
             elseif (i==2) then
                
                if (strghtOrNI==1) write(32,*) ka_prgrsnStrct3(ProgCyl,j)  ,A0_prgrsnStrct3(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ka_prgrsnNIstrct3(ProgCyl,j),A0_prgrsnNIstrct3(ProgCyl,j)
                
             elseif (i==3) then
                
                if (strghtOrNI==1) write(32,*) CgX_prgrsnStrct3(ProgCyl,j)   ,CgY_prgrsnStrct3(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) CgX_prgrsnNIstrct3(ProgCyl,j) ,CgY_prgrsnNIstrct3(ProgCyl,j)
                
             endif
             
          elseif (strctV==4) then
             
             if (i==1) then
                
                if (strghtOrNI==1) write(32,*) ks_prgrsnStrct4(ProgCyl,j)  ,l0_prgrsnStrct4(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ks_prgrsnNIstrct4(ProgCyl,j),l0_prgrsnNIstrct4(ProgCyl,j)
                
             elseif (i==2) then
                
                if (strghtOrNI==1) write(32,*) ka_prgrsnStrct4(ProgCyl,j)  ,A0_prgrsnStrct4(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) ka_prgrsnNIstrct4(ProgCyl,j),A0_prgrsnNIstrct4(ProgCyl,j)
                
             elseif (i==3) then
                
                if (strghtOrNI==1) write(32,*) CgX_prgrsnStrct4(ProgCyl,j)   ,CgY_prgrsnStrct4(ProgCyl,j)
                if (strghtOrNI==2) write(32,*) CgX_prgrsnNIstrct4(ProgCyl,j) ,CgY_prgrsnNIstrct4(ProgCyl,j)
                
             endif
             
          endif
          
       enddo
       
       write(32,*) " "
       
    enddo
    
    close(32)
    
  end subroutine writeProgrsnProps
  
  subroutine readProgrsnProps(ProgCyl,strctV,strghtOrNI,ksV,l0V,kaV,A0V,CgXV)
    implicit none
    integer, intent(in) :: ProgCyl
    integer, intent(in) :: strctV
    integer, intent(in) :: strghtOrNI
    
    real*8, intent(out) :: ksV(1:N_spr),kaV(1:N_cell)
    real*8, intent(out) :: l0V(1:N_spr),A0V(1:N_cell)
    real*8, intent(out) :: CgXV(1:N_node)
    
    character(len=100)  :: rflnm
    character(len=100)  :: rflnmbr
    character(len=100)  :: rfull_flnm
    
    integer :: i,j,jmax
    integer :: N_itm
    
    N_itm = 3
    
    if (strghtOrNI==1) rflnm = 'PrgrsnCyclProp'
    if (strghtOrNI==2) rflnm = 'PrgrsnCyclPropNI'
    
    write(rflnmbr,'(i1.1,a,i1.1,a)') ProgCyl,'Strct',strctV,'.dat'
    
    rfull_flnm=trim(adjustl(rflnm))//trim(adjustl(rflnmbr))
    write(*,*)trim(adjustl(rfull_flnm))
    
    open(unit=39,file=trim(adjustl(rfull_flnm)))
    
    do i = 1,N_itm
       
       if (i==1) jmax = N_spr
       if (i==2) jmax = N_cell
       if (i==3) jmax = N_node
       
       do j = 1,jmax
          
          if (i==1) then
             read(39,*) ksV(j),l0V(j)
          elseif (i==2) then
             read(39,*) kaV(j),A0V(j)
          elseif (i==3) then 
             read(39,*) CgXV(j)
          endif
          
       enddo
       
    enddo
    
    close(39)
    
    
  end subroutine readProgrsnProps
  
  subroutine read_strctProps_withAddedCell(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "Not for stage 1, type 1, routine: read_strctProps_withAddedCell"
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
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) read(21,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j)
             if (i==2) read(21,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j)
             if (i==3) read(21,*) CgX_strctTN(strctNum,j),CgY_strctTN(strctNum,j)
             
          elseif (modelID==2) then
             if (i==1) read(21,*) ks_strctNI(strctNum,j),l0_strctNI(strctNum,j)
             if (i==2) read(21,*) ka_strctNI(strctNum,j),A0_strctNI(strctNum,j)
             if (i==3) read(21,*) CgX_strctNI(strctNum,j),CgY_strctNI(strctNum,j)
          endif
          
       enddo
       
    enddo
    
    close(21)
    
  end subroutine read_strctProps_withAddedCell
  
  
  subroutine replaceProps_however_retain_only_Cg()
    implicit none
    integer :: strctNum
    real*8  :: CgXNodeVal(1:N_node),CgYNodeVal(1:N_node)
    integer :: choiceVal
    
    choiceVal = 4
    
    CgXNodeVal = CgXNode
    CgYNodeVal = CgYNode
    
    
    if (choiceVal==0) then
       strctNum = 13
       call replacePropsFrm_strctProps_withAddedCell(strctNum)
    else
       call readPrpsFrmDiffrntFiles(choiceVal)
    endif
    
    CgXNode = CgXNodeVal
    CgYNode = CgYNodeVal
    
  end subroutine replaceProps_however_retain_only_Cg
  
  subroutine modify_prp_frm_CM3PCS_wo_Highks_value_and_change_stpBystp
    implicit none
    real*8             :: kaV_Fl1(1:N_cell),A0V_Fl1(1:N_cell)
    real*8             :: ksV_Fl1(1:N_spr), l0V_Fl1(1:N_spr)
    real*8             :: CgXN_Fl1(1:N_node),CgYN_Fl1(1:N_node)
    
    real*8             :: kaV_Fl2(1:N_cell),A0V_Fl2(1:N_cell)
    real*8             :: ksV_Fl2(1:N_spr), l0V_Fl2(1:N_spr)
    real*8             :: CgXN_Fl2(1:N_node),CgYN_Fl2(1:N_node)
    
    integer            :: i,j,jmax
    integer            :: cellNm_1,sprNm_1,cellNm_2,sprNm_2
    
    real*8             :: PressV_1,kaV_1,Av_1,A0V_1
    real*8             :: TnsnV_1,ksV_1,lv_1,l0V_1
    real*8             :: PressV_2,kaV_2,Av_2,A0V_2
    real*8             :: TnsnV_2,ksV_2,lv_2,l0V_2
    real*8             :: nodeNm_1,nodeTyp_1,CgXN_1,CgYN_1
    real*8             :: nodeNm_2,nodeTyp_2,CgXN_2,CgYN_2
    
    character(len=100) :: cellFile1,sprFile1,nodeFile1
    character(len=100) :: cellFile2,sprFile2,nodeFile2
    
    cellFile1='cellPrps_CM3PCS.dat'
    sprFile1='sprsPrps_CM3PCS.dat'
    nodeFile1='nodePrps_CM3PCS.dat'
    
    cellFile2='cellPrps_CM3PCS_HighKs.dat'
    sprFile2='sprsPrps_CM3PCS_HighKs.dat'
    nodeFile2='nodePrps_CM3PCS_HighKs.dat'
    
    open(unit=343,file=trim(adjustl(cellFile1)))
    open(unit=443,file=trim(adjustl(sprFile1)))
    open(unit=743,file=trim(adjustl(nodeFile1)))
    
    open(unit=543,file=trim(adjustl(cellFile2)))
    open(unit=643,file=trim(adjustl(sprFile2)))
    open(unit=843,file=trim(adjustl(nodeFile2)))
    
    do i = 1,3
       
       if (i==1) jmax=N_cell
       if (i==2) jmax=N_spr
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) read(343,*) cellNm_1,PressV_1,kaV_1,Av_1, A0V_1
          if (i==1) read(543,*) cellNm_2,PressV_2,kaV_2,Av_2, A0V_2
          
          if (i==2) read(443,*) TnsnV_1, ksV_1,   lV_1, l0V_1,sprNm_1
          if (i==2) read(643,*) TnsnV_2, ksV_2,   lV_2, l0V_2,sprNm_2
          
          if (i==3) read(743,*) nodeNm_1,nodeTyp_1,CgXN_1,CgYN_1
          if (i==3) read(843,*) nodeNm_2,nodeTyp_2,CgXN_2,CgYN_2 
          
          if (i==1) then   
             kaV_Fl1(j)  = kaV_1  ; kaV_Fl2(j)  = kaV_2
             A0V_Fl1(j)  = A0V_1  ; A0V_Fl2(j)  = A0V_2
             
          elseif (i==2) then
             ksV_Fl1(j) = ksV_1 ; ksV_Fl2(j) = ksV_2
             l0V_Fl1(j) = l0V_1 ; l0V_Fl2(j) = l0V_2
             
          elseif (i==3) then
             CgXN_Fl1(j) = CgXN_1 ; CgXN_Fl2(j) = CgXN_2
             CgYN_Fl1(j) = CgYN_1 ; CgYN_Fl2(j) = CgYN_2
          endif
          
       enddo
       
    enddo
    
    close(343)
    close(443)
    close(743)
    close(543)
    close(643)
    close(843)
    
    k_area(1:N_cell)  = kaV_Fl1(1:N_cell)
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI-1,"Frame_NO aft kaV"
    
    k_spr(1:N_spr)    = ksV_Fl1(1:N_spr)
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI-1,"Frame_NO aft ksV"
    
    A0(1:N_cell)      = A0V_Fl2(1:N_cell)
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI-1,"Frame_NO aft A0V"

    CgXNode(1:N_node) = CgXN_Fl1(1:N_node) ; CgYNode(1:N_node) = CgYN_Fl1(1:N_node)
    call Equilibrate_only_NI_model ; write(*,*) Frame_NI-1,"Frame_NO aft CgX-CgY"
    
  end subroutine modify_prp_frm_CM3PCS_wo_Highks_value_and_change_stpBystp
    
  subroutine readPrpsFrmDiffrntFiles(choiceVal)
    implicit none
    integer, intent(in) :: choiceVal
    character(len=100)  :: flnm01,flnm02,flnm03,flnm04,flnm05,flnm06,flnm07,flnm08,flnm09,flnm10 ,flnm11
    integer             :: flTyp1,flTyp2,flTyp3,flTyp4,flTyp5,flTyp6,flTyp7,flTyp8,flTyp9,flTyp10,flTyp11
    
    integer             :: i,j
    
    flnm01='strctPropsAddedCellNI13Bfrpressure_adjstmnt_at_bottomCells.dat' ; flTyp1=0
    flnm02='strctPropsAddedCellNI13S1T1_BkUp.dat'                           ; flTyp2=1
    flnm03='strctPropsAddedCellNI13S1T1_TC.dat'                             ; flTyp3=1
    flnm04='strctPropsAddedCellNI13S1T1_BfrCgAppl.dat'                      ; flTyp4=0
    flnm05='strctPropsAddedCellNI13WithLowerThanSetValCell10.dat'           ; flTyp5=0
    flnm06='strctPropsAddedCellNI13WithLowerThanSetValCell1011.dat'         ; flTyp6=0
    
    if (choiceVal==1) call replaceProps_frmAFile(flnm01,flTyp1)
    if (choiceVal==2) call replaceProps_frmAFile(flnm02,flTyp2)
    if (choiceVal==3) call replaceProps_frmAFile(flnm03,flTyp3)
    if (choiceVal==4) call replaceProps_frmAFile(flnm04,flTyp4)
    if (choiceVal==5) call replaceProps_frmAFile(flnm05,flTyp5)
    if (choiceVal==6) call replaceProps_frmAFile(flnm06,flTyp6)
    
  end subroutine readPrpsFrmDiffrntFiles
  
  
  subroutine replacePropsFrm_strctProps_withAddedCell(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    N_itm = 3
    
    if (modelID==1) flnm='strctPropsAddedCellTN'
    if (modelID==2) flnm='strctPropsAddedCellNI'
    
    if (strctNum.le.9)                      write(flnmbr,'(i1.1,a)') strctNum,'S1T1.dat'
    if (strctNum.gt.9 .and. strctNum.le.99) write(flnmbr,'(i2.2,a)') strctNum,'S1T1.dat'
    
    write(*,*) trim(adjustl(flnm)),trim(adjustl(flnmbr))
    full_flnm=trim(adjustl(flnm))//trim(adjustl(flnmbr))
    
    open(unit=21,file=trim(adjustl(full_flnm)))
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) read(21,*) k_spr(j) ,l0(j)
          if (i==2) read(21,*) k_area(j),A0(j)
          if (i==3) read(21,*) CgXNode(j),CgYNode(j)
          
       enddo
       
    enddo
    
    close(21)
    
  end subroutine replacePropsFrm_strctProps_withAddedCell
  
  subroutine replaceProps_frmAFile(flNm,flTyp)
    implicit none
    character(len=100), intent(in) :: flNm
    integer, intent(in)            :: flTyp

    integer :: i,j,jmax
    integer :: N_itm

    N_itm = 3
    
    open(unit=21,file=trim(adjustl(flNm)))
    
    do i=1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       
       do j = 1,jmax
          
          if (i==1) read(21,*) k_spr(j) ,l0(j)
          if (i==2) read(21,*) k_area(j),A0(j)
          
          if (i==3) then
             
             if (flTyp==0) then
                
                if (i==3) read(21,*) CgXNode(j)
                CgYNode(j) = 0.0d0
                
             elseif (flTyp==1) then
                if (i==3) read(21,*) CgXNode(j),CgYNode(j)
             endif
             
          endif
          
       enddo
       
    enddo
    
    close(21)
    
  end subroutine replaceProps_frmAFile
  
  
  subroutine read_strctPropsIncludingNodeXY_withAddedCell(strctNum)
    implicit none
    integer, intent(in) :: strctNum
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm
    character(len=100) :: flnmbr
    character(len=100) :: full_flnm
    
    if (stageNo==1 .and. stageType==1) then
       continue
    else
       write(*,*) "Not for stg1,type1, routine:read_strctPropsIncludingNodeXY_withAddedCell"
    endif
    
    N_itm = 4
    
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
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) read(21,*) ks_strctTN(strctNum,j),l0_strctTN(strctNum,j)
             if (i==2) read(21,*) ka_strctTN(strctNum,j),A0_strctTN(strctNum,j)
             if (i==3) read(21,*) CgX_strctTN(strctNum,j),CgY_strctTN(strctNum,j)
             if (i==4) read(21,*) nodeXY_strctTN(strctNum,j,1:N_dmnsn)
             
          elseif (modelID==2) then
             if (i==1) read(21,*) ks_strctNI(strctNum,j),l0_strctNI(strctNum,j)
             if (i==2) read(21,*) ka_strctNI(strctNum,j),A0_strctNI(strctNum,j)
             if (i==3) read(21,*) CgX_strctNI(strctNum,j),CgY_strctNI(strctNum,j)
             if (i==4) read(21,*) nodeXY_strctNI(strctNum,j,1:N_dmnsn)
             
          endif
          
       enddo
       
    enddo
    
    close(21)
    
  end subroutine read_strctPropsIncludingNodeXY_withAddedCell
  
  subroutine read_strctProps_forWTpropsOfMultiplCycl
    
    !***READ THE COMMENT ABOUT THIS ROUTINE
    ! This routine is to read properties frm '(ks/l0/ks/A0)WT(1/2/3)Cycl(1/2/3/4/5/6/7)CS(1/2/3/4)'
    ! We need some values but that will created in later modules, so have created some dummy or
    ! passing variable to solve this problem instead of massively resuffling routines
    ! *** The variables are listed in the beginning of this file in module section
    
    implicit none
    
    
  end subroutine read_strctProps_forWTpropsOfMultiplCycl


  subroutine prep_vars_for_Extrplt_loctn1(NonUniorUni,flnmInit,flnmFinl)
    implicit none
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
    
    
  end subroutine prep_vars_for_Extrplt_loctn1
  
  
  
  subroutine saveA0l0_ofStrct_AddedCell_TNorNI(struct_No,saveWhom)
    implicit none
    integer, intent(in) :: struct_No
    integer, intent(in) :: saveWhom
    real*8 :: E
    
    if (saveWhom==1) then !save both TN and NI
       call saveA0l0_ofStrct_AddedCell(struct_No)
       
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
       call Equilibrate_system
       call save_config_and_generate_ColorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       !Frame_NI = Frame_NI+1
       
       call saveA0l0_ofStrct_AddedCell(struct_No)
       
       call deallocate_repetitive_arrays
       call switchback_to_TN_model
       
    elseif (saveWhom==2) then !save TN only
       call saveA0l0_ofStrct_AddedCell(struct_No)
    elseif (saveWhom==3) then ! save NI only
       call saveA0l0_ofStrct_AddedCell(struct_No)
    endif
    
  end subroutine saveA0l0_ofStrct_AddedCell_TNorNI
  
  
  subroutine saveA0l0_ofStrct_AddedCell(struct_No)
    implicit none
    integer, intent(in) :: struct_No
    
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm,flnmPrf,flnmSuf
    
    !This routine is for S1T1 only
    
    if (modelID == 1) flnmPrf='strctPropsAddedCellTN'
    if (modelID == 2) flnmPrf='strctPropsAddedCellNI'
    
    if (struct_No.gt.0 .and. struct_No.le.9) then
       write(flnmSuf,'(i1.1,a)') struct_No,'S1T1.dat'
    elseif (struct_No.gt.9 .and. struct_No.le.99) then
       write(flnmSuf,'(i2.2,a)') struct_No,'S1T1.dat'
    endif
    
    flnm=trim(adjustl(flnmPrf))//trim(adjustl(flnmSuf))
    write(*,*)trim(adjustl(flnmPrf))//trim(adjustl(flnmSuf))
    
    open(unit=311,file=trim(adjustl(flnm)))
    
    call coordntes_to_nodes(coordntes_xy,node_xy)
    call coordntesxyl0A0_to_coordntes(coordntes_xy,l0,A0,coordntes)
    
    if (modelID == 1) then
       
       l0_strctTN(struct_No,1:N_spr)   = l0(1:N_spr)
       ks_strctTN(struct_No,1:N_spr)   = k_spr(1:N_spr)
       A0_strctTN(struct_No,1:N_cell)  = A0(1:N_cell)
       ka_strctTN(struct_No,1:N_cell)  = k_area(1:N_cell)
       CgX_strctTN(struct_No,1:N_node) = CgXNode(1:N_node)
       CgY_strctTN(struct_No,1:N_node) = CgYNode(1:N_node)
       
       nodeXY_strctTN(struct_No,1:N_node,1:N_dmnsn) = node_xy(1:N_node,1:N_dmnsn)
       
    elseif (modelID == 2) then
       
       l0_strctNI(struct_No,1:N_spr)   = l0(1:N_spr)
       ks_strctNI(struct_No,1:N_spr)   = k_spr(1:N_spr)
       A0_strctNI(struct_No,1:N_cell)  = A0(1:N_cell)
       ka_strctNI(struct_No,1:N_cell)  = k_area(1:N_cell)
       CgX_strctNI(struct_No,1:N_node) = CgXNode(1:N_node)
       CgY_strctNI(struct_No,1:N_node) = CgYNode(1:N_node)
       
       nodeXY_strctNI(struct_No,1:N_node,1:N_dmnsn) = node_xy(1:N_node,1:N_dmnsn)
    endif
    
    N_itm = 4
    
    do i = 1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) write(311,*) ks_strctTN(struct_No,j)  ,l0_strctTN(struct_No,j)
             if (i==2) write(311,*) ka_strctTN(struct_No,j)  ,A0_strctTN(struct_No,j)
             if (i==3) write(311,*) CgX_strctTN(struct_No,j) ,CgY_strctTN(struct_No,j)
             if (i==4) write(311,*) nodeXY_strctTN(struct_No,j,1:N_dmnsn)
             
          elseif (modelID==2) then
             if (i==1) write(311,*) ks_strctNI(struct_No,j)  ,l0_strctNI(struct_No,j)
             if (i==2) write(311,*) ka_strctNI(struct_No,j)  ,A0_strctNI(struct_No,j)
             if (i==3) write(311,*) CgX_strctNI(struct_No,j) ,CgY_strctNI(struct_No,j)
             if (i==4) write(311,*) nodeXY_strctNI(struct_No,j,1:N_dmnsn)

             if (i==3 .and. j==50)  write(*,*) CgX_strctNI(struct_No,j) ,CgY_strctNI(struct_No,j),"CgXYchck"
          endif
          
       enddo
       
       write(311,*) " "
       
    enddo
    
    close(311)
    
  end subroutine saveA0l0_ofStrct_AddedCell
  
  
  
  subroutine saveA0l0_ofStrct_AddedCell_withknownStrct(struct_No)
    implicit none
    integer, intent(in) :: struct_No
    
    integer :: i,j,jmax
    integer :: N_itm
    
    character(len=100) :: flnm,flnmPrf,flnmSuf
    
    !This routine is for S1T1 only
    
    if (modelID == 1) flnmPrf='strctPropsAddedCellTN'
    if (modelID == 2) flnmPrf='strctPropsAddedCellNI'
       
    if (struct_No.gt.0 .and. struct_No.le.9) then
       write(flnmSuf,'(i1.1,a)') struct_No,'S1T1.dat'
    elseif (struct_No.gt.9 .and. struct_No.le.99) then
       write(flnmSuf,'(i2.2,a)') struct_No,'S1T1.dat'
    endif
    
    flnm=trim(adjustl(flnmPrf))//trim(adjustl(flnmSuf))
    write(*,*)trim(adjustl(flnmPrf))//trim(adjustl(flnmSuf))
    
    open(unit=311,file=trim(adjustl(flnm)))
      
    N_itm = 4
    
    do i = 1,N_itm
       
       if (i==1) jmax=N_spr
       if (i==2) jmax=N_cell
       if (i==3) jmax=N_node
       if (i==4) jmax=N_node
       
       do j = 1,jmax
          
          if (modelID==1) then
             if (i==1) write(311,*) ks_strctTN(struct_No,j),l0_strctTN(struct_No,j)
             if (i==2) write(311,*) ka_strctTN(struct_No,j),A0_strctTN(struct_No,j)
             if (i==3) write(311,*) CgX_strctTN(struct_No,j),CgY_strctTN(struct_No,j)
             if (i==4) write(311,*) nodeXY_strctTN(struct_No,j,1:N_dmnsn)
             
          elseif (modelID==2) then
             if (i==1) write(311,*) ks_strctNI(struct_No,j),l0_strctNI(struct_No,j)
             if (i==2) write(311,*) ka_strctNI(struct_No,j),A0_strctNI(struct_No,j)
             if (i==3) write(311,*) CgX_strctNI(struct_No,j),CgY_strctNI(struct_No,j)
             if (i==4) write(311,*) nodeXY_strctNI(struct_No,j,1:N_dmnsn)
             
          endif
          
       enddo
       
       write(311,*) " "
       
    enddo
    
    close(311)
    
  end subroutine saveA0l0_ofStrct_AddedCell_withknownStrct
  
  
  subroutine writeDatFileForPressure_And_SpecificTnsn(cntrlSt,PressV,TnsnV)
    implicit none
    integer, intent(in) :: cntrlSt
    real*8 , intent(in) :: PressV(1:N_cell)
    real*8 , intent(in) :: TnsnV(1:N_spr)
    
    integer :: i,j,imax,jmax
    integer :: middleLtrlSpr
    integer :: cell_no,nsprsCntInEachCell
    
    character(len=100) :: full_flnmT,flnmT,flnmbrT
    character(len=100) :: full_flnmP,flnmP,flnmbrP
    
    if (modelID==1) then
       write(*,*) "routine:is for modelID=2"
       stop
    endif
    
    write(flnmP,*)'St'
    write(flnmT,*)'St'
    
    write(flnmbrP,'(i1.1,a)') cntrlSt,'Press.dat'
    write(flnmbrT,'(i1.1,a)') cntrlSt,'Tnsn.dat'
    
    full_flnmP=trim(adjustl(flnmP))//trim(adjustl(flnmbrP))
    full_flnmT=trim(adjustl(flnmT))//trim(adjustl(flnmbrT))
    
    write(*,*) trim(adjustl(flnmP))//trim(adjustl(flnmbrP))
    write(*,*) trim(adjustl(flnmT))//trim(adjustl(flnmbrT))
    
    open(unit=52,file=trim(adjustl(flnmP))//trim(adjustl(flnmbrP)))
    open(unit=53,file=trim(adjustl(flnmT))//trim(adjustl(flnmbrT)))
    
    open(unit=54,file='writeDatFileForPressure.dat')
    
    nsprsCntInEachcell = 1+(NAEC+1)+(NAEC+1) 
    imax = Hlf_Ncell+1
    
    write(54,*) nsprsCntInEachCell,imax,"nspsCnt"
    
    do i = 1,imax
       
       if (i.ne.imax) cell_no = i
       if (i == imax) cell_no = N_cell
       
       if (i.ne.imax) then
          middleLtrlSpr = (cell_no-1)*nsprsCntInEachcell+ 1+(NAEC+1)+ ((NAEC+1)/2)+1
       elseif (i==imax) then
          middleLtrlSpr = (cell_no-1)*nsprsCntInEachcell+ ((NAEC+1)/2)+1
       endif
       
       write(52,*) i,PressV(cell_no)
       write(53,*) i,TnsnV(middleLtrlSpr)
       write(54,*) middleLtrlSpr,cell_no,i,"middleLtrl"
       
    enddo
    
    close(52)
    close(53)
    close(54)
    
  end subroutine writeDatFileForPressure_And_SpecificTnsn

  
  subroutine decision_of_meeting_NI(lgcl_meet)
    implicit none
    logical, intent(out) :: lgcl_meet
    
    real*8  :: tol_rdfn
    logical :: lgcl_rdfn 
    
    lgcl_meet = .False.
    
    tol_rdfn  = 0.14d0
    lgcl_rdfn = .False.
    
    call get_decision_of_redefining(tol_Rdfn,lgcl_rdfn)
    
    if (lgcl_rdfn .eqv. .True.) then
       lgcl_meet = .True.
    elseif (lgcl_rdfn .eqv. .False.) then
       lgcl_meet = .False.
    endif
    
  end subroutine decision_of_meeting_NI

  subroutine find_distance_frm_pulley(distnce)
    implicit none
    real*8, intent(out) :: distnce
    
    integer :: cellNmL,cellNmR
    integer :: nodeNmL,nodeNmR
    real*8  :: disFrmPulleyL,disFrmPulleyR
    
    cellNmL = invgnated_cells(1) ; cellNmR = invgnated_cells(2)
    
    nodeNmL = (2*cellNmL)-1
    nodeNmR = (Hlf_Ncell+1)*2 + (2*cellNmL)-1
    
    write(*,*) cellNmL,cellNmR,"cellNms"
    write(*,*) nodeNmL,nodeNmR,"nodeNms"

    write(*,*) node_xy(nodeNmL,1:2),"L nodes"
    write(*,*) node_xy(nodeNmR,1:2),"R nodes"
    
    disFrmPulleyL = abs(node_xy(nodeNmL,1) - 0.00d0)
    disFrmPulleyR = abs(node_xy(nodeNmR,1) - 0.00d0)
    
    write(*,*) disFrmPulleyL,disFrmPulleyR,"dis"
    distnce = disFrmPulleyL
    
  end subroutine find_distance_frm_pulley
  
  
  subroutine save_NI_ifNeeded(savingFlg)
    implicit none
    integer, intent(in) :: savingFlg
    
    
    if (savingFlg==1) then
       continue
    elseif (savingFlg==0) then
       continue
    else
       write(*,*) "savingFlg has unrealistic value, which is =",savingFlg
       stop
    endif
    
    
  end subroutine save_NI_ifNeeded

  subroutine print_pressTnsn(countLoop,PresVal,TnsnVal)
    implicit none
    integer, intent(in) :: countLoop
    real*8,  intent(in) :: PresVal(1:N_cell)
    real*8,  intent(in) :: TnsnVal(1:N_spr)
    
    integer :: i,j,jmax
    
    open(unit=87,file='print_presTnsn.dat',position='append')
    
    write(87,*) countLoop,"countLoop"
    write(87,*) " "
    
    do i = 1,2
       
       if (i==1) jmax = N_cell
       if (i==2) jmax = N_spr
       
       do j = 1,jmax
          if (i==1) write(87,*) PresVal(j),j
          if (i==2) write(87,*) TnsnVal(j),j
       enddo
       
       write(87,*) " "
       
    enddo
    
    close(87)
    
  end subroutine print_pressTnsn
  
  
  subroutine tilted_rectngle_to_rectngleMatchedwithPage(cellNm,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    !This is also a semi non-generalized_routine
    
    integer :: cellNmL,cellNmR
    integer :: nxtCellNmL,nxtCellNmR
    integer :: nsprsInACell,sprNmL,sprNmR
    real*8  :: E,hmuchA0,hmuchKa,hmuchL0,hmuchKs
    real*8  :: area_chngs,areaBfrChngs,areaAftChngs,area_chngsPrcnt
    real*8  :: area_recoversPrcnt,TOL_areaPrcnt,TOL_forReductn
    real*8  :: distance
    integer :: count,i
    integer :: whtSpr
    integer :: incrORdcrs
    
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
    
    
    if (modelID==1) then
       call switch_to_NI_model
       E = Energy(node_xy,l0,A0)
       call get_gradient(node_xy,l0,A0,gradient)
    elseif (modelID==2) then
       continue
    endif
    
    cellNmL = cellNm 
    cellNmR = (cellNmL)+(Hlf_Ncell)
    
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + 2*(NAEC_Ltrl+1)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    write(*,*) (FrmNo-1),"Before Changing"
    
    hmuchA0      = 1.20d0 ; hmuchL0 = 0.96d0 
    
    A0(cellNmL)  = (hmuchA0)*A0(cellNmL)
    A0(cellNmR)  = A0(cellNmL)
    areaBfrChngs = A(cellNmL)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    area_chngs      = abs(A(cellNmL)-areaBfrChngs)
    area_chngsPrcnt = (area_chngs/areaBfrChngs)*100.0d0
    count           = 0
    
    !area_chngs is difference, areaBfrChngs and areaAftChngs are the actual areas
    write(*,*) area_chngs,area_chngsPrcnt,"area_chngs and area_chngsPrcnt"
    
    TOL_areaPrcnt = 1.00d0 ; TOL_forReductn = 5.0d0
    
    !*** reducing l0 and k of DMs if the nodes are too close to Pulley Node ***
    write(*,*) (FrmNo-1),"FrmNo Bfr manplt in tilted"
    
    call find_distance_frm_pulley(distance)
    if (distance .le. 0.05d0) then
       incrORdcrs = 2
       call manplt_apcl_membrne_tnsn_invgnting_and_botmCells(incrORdcrs,ExpNo,FrmNo)
       write(*,*) (FrmNo-1),"FrmNo Aft manplt in tilted"
    endif
    !****
    
    do
       
       do i = 1,nsprsInACell   
          sprNmL     = area_spr(cellNmL,i)  ; sprNmR     = area_spr(cellNmR,i)
          l0(sprNmL) = (hmuchL0)*l0(sprNmL) ; l0(sprNmR) = l0(sprNmL)
       enddo
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       areaAftChngs       = A(cellNmL)
       area_recoversPrcnt = (abs(areaAftChngs-areaBfrChngs)/areaBfrChngs)*100.0d0
       
       if (area_recoversPrcnt .le. TOL_areaPrcnt) exit
       if (A(cellNm) .le. areaBfrChngs) exit !IMPORTANT EXIT CONDITION, IF I decrease A0 at first, THEN the (le) will be (gt) 
       
       if ((area_recoversPrcnt.gt.TOL_areaPrcnt).and.(area_recoversPrcnt.gt.TOL_forReductn)) then
          if (count==0) then
             hmuchL0 = 1.0d0-((1.0d0-hmuchL0)/2.0d0)
             count=1
             write(*,*) hmuchL0,"hmuchL0 Changes"
          endif
       endif
       
    enddo
    
    
    
    whtSpr = 4 ; hmuchKs = 5.0d0 ; hmuchL0 = 0.80d0
    call change_sprsProp_of_a_cell(cellNmL+1,whtSpr,hmuchKs,hmuchL0)
    call change_sprsProp_of_a_cell(cellNmR+1,whtSpr,hmuchKs,hmuchL0)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    write(*,*) (FrmNo-1),"FrmNo after the tilted routine"
    
    
    hmuchA0      = 1.20d0 
    nxtCellNmL   = cellNmL+1 ;  nxtCellNmR = cellNmR+1 
    
    A0(nxtcellNmL)  = (hmuchA0)*A0(nxtcellNmL)
    A0(nxtcellNmR)  = A0(nxtcellNmL)
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
    hmuchL0=0.33d0 ; hmuchKs=3.0d0
    call change_l0sOf_OnlyInitiatorCell_BasalSpr(hmuchL0,hmuchKs,ExpNo,FrmNo)
    
    write(*,*) "AFTER INITIATOR CELL CHANGE"
    write(*,*) "SLEEPING FOR 1 seconds"
    call sleep(1)
    
  end subroutine tilted_rectngle_to_rectngleMatchedwithPage
  
  
  subroutine change_sprsProp_of_a_cell(cellNm,whtSpr,hmuchKs,hmuchL0)
    implicit none
    integer, intent(in) :: cellNm
    integer, intent(in) :: whtSpr
    real*8 , intent(in) :: hmuchKs
    real*8 , intent(in) :: hmuchL0
    
    integer :: nsprsInACell,nsprsInThreeSide
    integer :: sprNm,i,imax
    integer, allocatable :: apcls(:),bsals(:),ltrls1(:),ltrls2(:)
    
    allocate(apcls(1:(NAEC_apcl+1)))
    allocate(bsals(1:(NAEC_Bsal+1)))
    allocate(ltrls1(1:(NAEC_Ltrl+1)))
    allocate(ltrls2(1:(NAEC_Ltrl+1)))
    
    apcls = 0 ; bsals = 0 ; ltrls1 = 0 ; ltrls2 = 0 
    
    nsprsInACell     = (NAEC_Apcl+1) + (NAEC_Bsal+1) + 2*(NAEC_Ltrl+1)
    nsprsInThreeSide = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    
    if (whtSpr==1) then
       imax = NAEC_Apcl+1
       do i = 1,imax
          apcls(i) = (cellNm-1)*(nsprsInThreeSide) + i
       enddo
    elseif (whtSpr==2) then
       imax = NAEC_Bsal+1
       do i = 1,imax
          bsals(i) = (cellNm-1)*(nsprsInThreeSide) + (NAEC_Apcl+1) + i
       enddo
    elseif (whtSpr==3) then
       imax = NAEC_Ltrl+1
       do i = 1,imax
          ltrls1(i) = (cellNm-2)*(nsprsInThreeSide) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       enddo
    elseif (whtSpr==4) then
       imax = NAEC_Ltrl+1
       do i = 1,imax
          ltrls2(i) = (cellNm-1)*(nsprsInThreeSide) + (NAEC_Apcl+1) + (NAEC_Bsal+1) + i
       enddo
    endif
    
    do i = 1,imax
       
       if (whtSpr==1) sprNm = apcls(i) 
       if (whtSpr==2) sprNm = bsals(i)
       if (whtSpr==3) sprNm = ltrls1(i)
       if (whtSpr==4) sprNm = ltrls2(i)
       
       l0(sprNm)    = (hmuchL0) * l0(sprNm)
       k_spr(sprNm) = (hmuchKs) * k_spr(sprNm)
       
    enddo
    
  end subroutine change_sprsProp_of_a_cell
  
  
  subroutine change_l0sOf_OnlyInitiatorCell_BasalSpr(hmuchL0,hmuchKs,ExpNo,FrmNo)
    implicit none
    real*8 , intent(in)    :: hmuchL0,hmuchKs
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: initiatorCell
    integer :: nsprsInACell
    integer :: i,imax
    integer :: sprNm
    
    initiatorCell = N_cell
    
    if (modelID==1) then
       nsprsInACell = 1+1+1
       imax         = 1
       
    elseif (modelID==2) then
       nsprsInACell  = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
       imax          = 3
    endif
    
    do i = 1,imax
       sprNm        = (initiatorCell-1)*(nsprsInACell) + i
       l0(sprNm)    = (hmuchL0)*l0(sprNm)
       k_spr(sprNm) = (hmuchKs)*k_spr(sprNm)
    enddo
    
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    
  end subroutine change_l0sOf_OnlyInitiatorCell_BasalSpr
  
  
  subroutine incrPressAndMaintainShapeOFInitiatorCell(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer :: ICB !ICB=InitiatorCellBottom
    real*8  :: areaBfrChng,areaAftChng
    real*8  :: prcntChngArea,tolPrcnt
    real*8  :: hmuchA0,hmuchl0,hmuchks
    integer :: sprsInICB,sprNm
    integer :: i,j
    integer :: cnt
    
    ICB         = N_cell
    areaBfrChng = 0.95d0*A(ICB)
    
    hmuchA0 = 1.20d0
    !A0(ICB) = (hmuchA0)*(A0(ICB))
    
    call Equilibrate_system
    call save_config_and_generate_ColorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
    FrmNo = FrmNo+1
    call switchto_NI_model_run_and_switchbackto_TN
    
    hmuchl0 = 1.0d0-(abs(hmuchA0-1.0d0)/4.0d0)
    hmuchks = 1.0d0+(abs(hmuchA0-1.0d0)/4.0d0)
    
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
          sprNm        = area_spr(ICB,i)
          l0(sprNm)    = (hmuchl0)*l0(sprNm)
          k_spr(sprNm) = (hmuchks)*k_spr(sprNm)
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
    
  end subroutine incrPressAndMaintainShapeOFInitiatorCell
  
  
  subroutine incoming_region_pressure_drop(ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    real*8  :: PresVal(1:N_cell)
    real*8  :: TnsnVal(1:N_spr)
    real*8  :: trgt_press
    integer :: trgt_cell
    
    real*8 , allocatable :: dP(:)
    integer, allocatable :: Lcells(:),Rcells(:),Lsprs(:),Rsprs(:)
    
    integer :: pressIncrOrDcrs
    integer :: N_incmingCell,N_incmingSprs
    logical :: lgcl_matched
    integer :: LoopCnt
    
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
         real*8  :: Pressure(1:N_cell)
         real*8  :: dum_Nodes(1:N_node,1:N_dmnsn)
         real*8  :: dumA0(1:N_cell)
       end function Pressure
       
    end interface
    
    
    PresVal(1:N_cell) = Pressure(node_xy,A0)
    TnsnVal(1:N_spr)  = TnsnComprsn(node_xy,l0)
    
    call  find_the_target_press(PresVal,trgt_cell,trgt_press)
    allocate(dP(1:trgt_cell)) ; dP(1:trgt_cell) = 0.0d0
    call get_the_diff_press(trgt_cell,trgt_press,PresVal,dP)
    
    call get_Nincmings(trgt_cell,N_incmingCell,N_incmingSprs)
    allocate(Lcells(1:N_incmingCell),Rcells(1:N_incmingCell))
    allocate(Lsprs(1:N_incmingSprs), Rsprs(1:N_incmingSprs))
    Lcells=0 ; Rcells=0 ; Lsprs=0 ; Rsprs=0
    call get_incming_cells_and_sprs(N_incmingCell,N_incmingSprs,Lcells,Rcells,Lsprs,Rsprs)
    
    LoopCnt = 1
    
    do 
       
       pressIncrOrDcrs = -1
       call get_incmingRgn_A0l0chnge(trgt_cell,pressIncrOrDcrs,N_incmingCell,N_incmingSprs,&
            Lcells,Rcells,Lsprs,Rsprs)
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,ExpNo,FrmNo)
       FrmNo = FrmNo+1
       
       call find_Pressure_and_Tension(PresVal,TnsnVal)
       call find_if_PresVal_matched_incmingRgn(PresVal,trgt_cell,trgt_press,N_incmingCell,&
            Lcells,Rcells,lgcl_matched)
       
       if (lgcl_matched .eqv. .True.) exit
       if (LoopCnt==5) exit
       
       LoopCnt = LoopCnt+1
       
    enddo
    stop
    
  end subroutine incoming_region_pressure_drop
  
  subroutine find_the_target_press(PresVal,trgt_cell,trgt_press)
    implicit none
    real*8,  intent(in)  :: PresVal(1:N_cell)
    integer, intent(out) :: trgt_cell
    real*8 , intent(out) :: trgt_press
    
    integer :: i,j
    integer :: nodeNm
    
    do i = 1,Hlf_Ncell
       nodeNm = area_node(i,0)
       if (nodeNm == max_node_area) then
          trgt_cell = (i-1)
          exit
       endif
    enddo
    
    trgt_press = PresVal(trgt_cell)
    
  end subroutine find_the_target_press
  
  subroutine get_the_diff_press(trgt_cell,trgt_press,PresVal,dP)
    implicit none
    integer, intent(in)  :: trgt_cell
    real*8 , intent(in)  :: trgt_press
    real*8 , intent(in)  :: PresVal(1:N_cell)
    real*8 , intent(out) :: dP(1:trgt_cell)
    
    integer :: i,j
    
    open(unit=432,file='diff_pressVal.dat')
    do i = 1,trgt_cell
       dP(i) = trgt_press - PresVal(i)
       write(432,*) i,trgt_press,PresVal(i),"cellNm,trgtP,PresCell"
    enddo
    close(432)
    
  end subroutine get_the_diff_press
  
  subroutine get_Nincmings(trgt_cell,N_incmingCell,N_incmingSprs)
    implicit none
    integer, intent(in)  :: trgt_cell
    integer, intent(out) :: N_incmingCell,N_incmingSprs
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    N_incmingCell = trgt_cell-1
    N_incmingSprs = (N_incmingCell)*(nsprsInACell)
    
  end subroutine get_Nincmings
  
  subroutine get_incmingRgn_A0l0chnge(trgt_cell,pressIncrOrDcrs,N_incmingCell,N_incmingSprs,&
       Lcells,Rcells,Lsprs,Rsprs)
    implicit none
    integer, intent(in) :: trgt_cell,pressIncrOrDcrs
    integer, intent(in) :: N_incmingCell,N_incmingSprs
    integer, intent(in) :: Lcells(1:N_incmingCell),Rcells(1:N_incmingCell)
    integer, intent(in) :: Lsprs(1:N_incmingSprs),Rsprs(1:N_incmingSprs)
    
    real*8  :: hmuchA0,hmuchL0
    integer :: nsprsInACell
    
    if (pressIncrOrDcrs==+1) then !pressure Increase
       hmuchA0 = +1.02d0 ; hmuchL0 = +0.98d0
    elseif (pressIncrOrDcrs==-1) then !pressure Decrease
       hmuchA0 = +0.98d0 ; hmuchL0 = +1.02d0
    endif
    
    !allocate(Lcells(1:N_incmingCell),Rcells(1:N_incmingCell))
    !allocate(Lsprs(1:N_incmingSprs), Rsprs(1:N_incmingSprs))
    !call get_incming_cells_and_sprs(N_incmingCell,N_incmingSprs,Lcells,Rcells,Lsprs,Rsprs)
    
    call alterA0_incmingCells(N_incmingCell,Lcells,Rcells,hmuchA0)
    call alterL0_incmingSprs(N_incmingSprs,Lsprs,Rsprs,hmuchL0)
    
  end subroutine get_incmingRgn_A0l0chnge
  
  
  subroutine get_incming_cells_and_sprs(N_incmingCell,N_incmingSprs,Lcells,Rcells,Lsprs,Rsprs)
    implicit none
    integer, intent(in)  :: N_incmingCell,N_incmingSprs
    integer, intent(out) :: Lcells(1:N_incmingCell),Rcells(1:N_incmingCell)
    integer, intent(out) :: Lsprs(1:N_incmingSprs),Rsprs(1:N_incmingSprs)
    
    integer :: i,j
    integer :: sprCnt
    integer :: nsprsInACell
    
    nsprsInACell = (NAEC_Apcl+1) + (NAEC_Bsal+1) + (NAEC_Ltrl+1)
    sprCnt       = 1
    
    open(unit=433,file='incmingCellAndSprList.dat')
    
    do i = 1,N_incmingCell
       Lcells(i) = i ; Rcells(i) = Hlf_Ncell+i
       
       write(433,*) Lcells(i),Rcells(i),i,"Lcells-Rcells"
       write(433,*) " "
       write(433,*) "List For The Springs :"
       
       do j = 1,nsprsInACell
          Lsprs(sprCnt) = (Lcells(i)-1)*(nsprsInACell)  + j
          Rsprs(sprCnt) = (Hlf_Ncell)*(nsprsInACell) + Lsprs(sprCnt)
          write(433,*) Lsprs(sprCnt),Rsprs(sprCnt),"Lsprs || Rsprs"
          sprCnt = sprCnt+1
       enddo
       
    enddo
    
    close(433)
    
  end subroutine get_incming_cells_and_sprs
  
  
  subroutine alterA0_incmingCells(N_incmingCell,Lcells,Rcells,hmuchA0)
    implicit none
    integer, intent(in) :: N_incmingCell
    integer, intent(in) :: Lcells(1:N_incmingCell),Rcells(1:N_incmingCell)
    real*8 , intent(in) :: hmuchA0
    
    integer :: i,j
    integer :: cellNmL,cellNmR
    real*8  :: prvA0
    
    open(unit=322,file='alterA0Chk.dat')
    
    do i = 1,N_incmingCell
       cellNmL = Lcells(i) ; cellNmR = Rcells(i)
       prvA0   = A0(cellNmL)
       
       A0(cellNmL) = (hmuchA0)*A0(cellNmL)
       A0(cellNmR) = A0(cellNmL)
       write(322,*) A0(cellNmL),A0(cellNmR),(A0(cellNmL)/prvA0),"A0 L,R and Ratio"
    enddo
    
    close(322)
    
  end subroutine alterA0_incmingCells
  
  
  subroutine alterL0_incmingSprs(N_incmingSprs,Lsprs,Rsprs,hmuchL0)
    implicit none
    integer, intent(in) :: N_incmingSprs
    integer, intent(in) :: Lsprs(1:N_incmingSprs),Rsprs(1:N_incmingSprs)
    real*8 , intent(in) :: hmuchL0
    
    integer :: i,j
    integer :: sprNmL,sprNmR
    real*8  :: prvL0
    
    open(unit=622,file='alterL0Chk.dat')
    
    do i = 1,N_incmingSprs
       sprNmL = Lsprs(i) ; sprNmR = Rsprs(i)
       prvL0  = l0(sprNmL)
       
       l0(sprNmL) = (hmuchL0)*l0(sprNmL)
       l0(sprNmR) = l0(sprNmL)
       write(622,*) l0(sprNmL),l0(sprNmR),(l0(sprNmL)/prvL0),"L0 L,R and Ratio"
    enddo
    
    close(622)
    
  end subroutine alterL0_incmingSprs
  
  
  subroutine find_if_PresVal_matched_incmingRgn(PresVal,trgt_cell,trgt_press,N_incmingCell,&
       Lcells,Rcells,lgcl_matched)
    implicit none
    real*8,  intent(in)  :: PresVal(1:N_cell)
    integer, intent(in)  :: trgt_cell
    real*8,  intent(in)  :: trgt_press
    integer, intent(in)  :: N_incmingCell
    integer, intent(in)  :: Lcells(1:N_incmingCell),Rcells(1:N_incmingCell)
    logical, intent(out) :: lgcl_matched
    
    integer :: count,cellNo,i,j
    integer :: prvOrAft,cntForlgcl
    real*8  :: trgt_pressForCmprsn
    
    prvOrAft = -1 ; call cnt_for_lgcl_match(trgt_cell,prvOrAft,cntForlgcl)
    write(*,*) cntForlgcl,"cntForlgcl"
    
    count = 0
    trgt_pressForCmprsn = (0.98d0)*(trgt_press)
    
    do i = 1,2
       
       do j = 1,2
          
          if (i==1) cellNo = Lcells(j)
          if (i==2) cellNo = Rcells(j)
          
          if (PresVal(cellNo) .le. trgt_pressForCmprsn) then
             continue
          elseif (PresVal(cellNo) .gt. trgt_pressForCmprsn) then
             count = count+1
          endif
          
       enddo
       
    enddo
    
    lgcl_matched = .False.
    
    if (count==cntForlgcl) then
       lgcl_matched = .True.
       
    elseif (count .lt.cntForlgcl) then
       lgcl_matched = .False.
    else
       write(*,*) "count value not okay",count
       stop
    endif
    
    write(*,*) count,"countVal"
    
  end subroutine find_if_PresVal_matched_incmingRgn
  
  
  subroutine cnt_for_lgcl_match(trgt_cell,prvOrAft,cntForlgcl)
    implicit none
    integer, intent(in)  :: trgt_cell,prvOrAft
    integer, intent(out) :: cntForlgcl
    
    integer :: prvCells,aftCells
    
    if (prvOrAft == (-1)) then
       prvCells   = trgt_cell-1
       cntForlgcl = 2*prvCells
       
    elseif (prvOrAft == (+1)) then
       aftCells   = Hlf_Ncell-trgt_cell
       cntForlgcl = 2*aftCells
    endif
    
  end subroutine cnt_for_lgcl_match
  
  
  subroutine make_control_state_perfect_for_supp_doc(cellNm) ! [[[INPUT IS ONLY LEFT SIDE CELL NUMBER]]]
    
    implicit none
    integer, intent(in) :: cellNm
    
    integer :: cellL,cellR,nsprsInACellNI,cntLp
    integer :: cell1,cell2
    real*8  :: hmL0Ap,hmL0Bs,hmL0Lt
    real*8  :: tolrncForStrghtNess,Ratio_SC
    integer :: cellVal,sprT
    integer :: FrmNoIncr,readTheFile
    
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    integer, allocatable :: ApclSp1(:),BsalSp1(:),LtrlSp1(:)
    integer, allocatable :: ApclSp2(:),BsalSp2(:),LtrlSp2(:)
    
    cellL=cellNm ; cellR=(cellL)+(Hlf_Ncell)
    if (modelID==1) call switch_to_NI_and_Equilibrate
    write(*,*) modelID,"modelID make_control_state_perfect_for_supp_doc"
    
    allocate(ApclSp(1:(NAEC_Apcl+1)), BsalSp(1:(NAEC_Bsal+1)), LtrlSp(1:(NAEC_Ltrl+1)))
    allocate(ApclSp1(1:(NAEC_Apcl+1)),BsalSp1(1:(NAEC_Bsal+1)),LtrlSp1(1:(NAEC_Ltrl+1)))
    allocate(ApclSp2(1:(NAEC_Apcl+1)),BsalSp2(1:(NAEC_Bsal+1)),LtrlSp2(1:(NAEC_Ltrl+1)))
    
    tolrncForStrghtNess = 0.99d0
    cntLp=1
    
    readTheFile = 1
    
    if (readTheFile == 0) then
       
       do
          
          if (cntLp.le.5) then
             
             if (cntLp==1) then
                call get_ApclBsalLtrl_OfNIsystm(cellL,ApclSp,BsalSp,LtrlSp)
                write(*,*) cellL,ApclSp,BsalSp,LtrlSp,"CellL and Sprs"
             endif
             
             hmL0Ap=1.05d0 ; hmL0Bs=0.95d0
             call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
             call manplt_spr_NIsys(BsalSp,NAEC_Bsal,hmL0Bs)
             
          elseif (cntLp.gt.5 .and. cntLp.le.10) then
             
             if (cntLp.gt.5 .and. cntLp.le.7) then
                A0(cellL)=0.95d0*A0(cellL) ; A0(cellR) = A0(cellL)
                
             elseif (cntLp.gt.7 .and. cntLp.le.10) then
                hmL0Ap=1.05d0 ; hmL0Bs=1.05d0
                call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
                call manplt_spr_NIsys(BsalSp,NAEC_Bsal,hmL0Bs)
             endif
             
          elseif (cntLp.gt.10 .and. cntLp.le.15) then
             A0(Hlf_Ncell-1)=0.95d0*A0(Hlf_Ncell-1) ; A0(N_cell-2) = A0(Hlf_Ncell-1)
             
          elseif (cntLp.gt.15 .and. cntLp.le.20) then
             
             if (cntLp==16) then
                call get_ApclBsalLtrl_OfNIsystm(cellL-1,ApclSp,BsalSp,LtrlSp)
                write(*,*) cellL-1,ApclSp,BsalSp,LtrlSp,"CellL-1 and Sprs"
             endif
             
             if (cntLp.gt.15 .and. cntLp.le.17) then 
                A0(cellL-1)=0.95d0*A0(cellL-1) ; A0(cellR-1) = A0(cellL-1)
                
             elseif (cntLp.gt.17 .and. cntLp.le.20) then
                hmL0Ap=1.05d0 ; hmL0Bs=1.05d0
                call manplt_spr_NIsys(ApclSp,NAEC_Apcl,hmL0Ap)
                call manplt_spr_NIsys(BsalSp,NAEC_Bsal,hmL0Bs)
             endif
             
          elseif (cntLp.gt.20 .and. cntLp.le.21) then
             
             call ApBsLt_DrawBoardMeasuredLengths()
             call ApBsLt_LengthRatios_BasedOnFirstCell()
             
             write(*,*) Frame_NI,"starts for 4"
             cellVal=4 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 5"
             cellVal=5 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 6"
             cellVal=6 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 7"
             cellVal=7 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 8"
             cellVal=8 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 9"
             cellVal=9 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             write(*,*) Frame_NI,"starts for 10"
             cellVal=10 ; sprT=2
             call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
             !write(*,*) Frame_NI,"starts for 11"
             !cellVal=11 ; sprT=2
             !call match_ratios_ApBsLt_BasedOnEXP_Images(cellVal,sprT)
             
          elseif (cntLp.gt.21 .and. cntLp.le.23) then
             cellVal=11
             A0(cellVal)=0.95d0*A0(cellVal) ; A0(cellVal+Hlf_Ncell) = A0(cellVal)
             
          elseif (cntLp.gt.23) then
             call shortn_Bsal_ofEndCell
          endif
          
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          Frame_NI = Frame_NI+1
          
          call get_Ratio_of_Strght_Curved_len_Of_LtrlSide(LtrlSp,Ratio_SC)
          
          !if (Ratio_SC .gt. tolrncForStrghtNess) then
          !   exit
          !endif
          
          
          if (cntLp==24) then
             write(*,*) "exits here for cntLp=24",cntLp
             exit
          endif
          
          cntLp = cntLp+1
          
       enddo
       
       
    elseif (readTheFile == 1) then
       
       write(*,*) Frame_NI,"Frame_NI bfr increasing"
       
       FrmNoIncr = 375
       Frame_NI  = Frame_NI+FrmNoIncr-1
       
       write(*,*) Frame_NI,"Frame_NI aft increasing"
       
       call read_config_and_start_simlnFrm_there(Exprmnt_NI,Frame_NI)
       call Equilibrate_only_NI_model

       call find_the_cell_with_diagonal_spr(cell1) ; cell2 = cell1+1
       call get_ApclBsalLtrl_OfNIsystm(cell1,ApclSp1,BsalSp1,LtrlSp1)
       call get_ApclBsalLtrl_OfNIsystm(cell2,ApclSp2,BsalSp2,LtrlSp2)
       
       do cntLp = 1,10 
          
          hmL0Ap=0.95d0
          call manplt_spr_NIsys(ApclSp1,NAEC_Apcl,hmL0Ap)
          call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap)
          
          call Equilibrate_system
          call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
          Frame_NI = Frame_NI+1
          
       enddo
       
    endif
    
    
  end subroutine make_control_state_perfect_for_supp_doc
  
  
  subroutine manplt_spr_NIsys_with_cellNm(cellNm,ApBsLt,cntMax,hmL0,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: cellNm,ApBsLt,cntMax
    real*8 , intent(in)    :: hmL0
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer, allocatable   :: Sprs(:)
    integer                :: cnt,NAEC_val
    integer                :: nsprsInACell
    
    
    if (ApBsLt==1) then
       allocate(Sprs(1:(NAEC_Apcl+1)))
       NAEC_val = NAEC_Apcl
       
    elseif (ApBsLt==2) then
       allocate(Sprs(1:(NAEC_Bsal+1)))
       NAEC_val = NAEC_Bsal
       
    elseif (ApBsLt==3) then
       allocate(Sprs(1:(NAEC_Ltrl+1)))
       NAEC_val = NAEC_Ltrl
    endif
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    do cnt = 1,(NAEC_val+1)
       if (ApBsLt==1) Sprs(cnt) = (cellNm-1)*(nsprsInACell) + cnt
       if (ApBsLt==2) Sprs(cnt) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + cnt
       if (ApBsLt==3) Sprs(cnt) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+NAEC_Bsal+2) + cnt
    enddo
    
    do cnt = 1,cntMax
       call manplt_spr_NIsys(Sprs,NAEC_val,hmL0)
       call Equilibrate_only_NI_model
    enddo
    
  end subroutine manplt_spr_NIsys_with_cellNm
  
  
  subroutine manplt_spr_NIsys_for_Ncell(ApBsLt,cntMax,hmL0,ExpNo,FrmNo)
    implicit none
    integer, intent(in)    :: ApBsLt,cntMax
    real*8 , intent(in)    :: hmL0
    integer, intent(in)    :: ExpNo
    integer, intent(inout) :: FrmNo
    
    integer, allocatable   :: Sprs(:)
    integer                :: cnt,NAEC_val
    integer                :: nsprsInACell
    integer                :: sprNm,numSprBfr
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    if (CellsMeet==0) then
       
       if (ApBsLt==1) then
          
          NAEC_val  = NAEC_Apcl
          numSprBfr = (N_cell-1)*(nsprsInACell)
          allocate(Sprs(1:(NAEC_val+1)))
          Sprs = 0
          
       elseif (ApBsLt==2) then
          
          NAEC_val = NAEC_Bsal
          numSprBfr = (N_cell-1)*(nsprsInACell) + (NAEC_Apcl+1)
          allocate(Sprs(1:(NAEC_val+1)))
          Sprs = 0
          
       elseif (ApBsLt==3) then
          write(*,*) "ApBsLt ca't be 3 here, sb: manplt_spr_NIsys_for_Ncell"
          stop
       endif
       
    elseif (CellsMeet.gt.0) then
       
       if (ApBsLt.ne.2) then
          write(*,*) "ApBsLt ca't be =",ApBsLt,"here, sb: manplt_spr_NIsys_for_Ncell"
       elseif (ApBsLt==2) then
          
          NAEC_val  = NAEC_Bsal
          numSprBfr = (N_cell-1)*(nsprsInACell)
          allocate(Sprs(1:(NAEC_val+1)))
          Sprs = 0
          
       endif
       
    endif
    
    do i = 1,(NAEC_val+1)
       Sprs(i) = numSprBfr + i
       write(*,*) i,Sprs(i),"Sprs in sb: manplt_spr_NIsys_for_Ncell"
    enddo
    
    do cnt = 1,cntMax
       
       do i = 1,(NAEC_val+1)
          sprNm     = Sprs(i)
          l0(sprNm) = (hmL0)*l0(sprNm)
          write(*,*) sprNm,l0(sprNm),"sprNm-L0s"
       enddo
       call Equilibrate_only_NI_model
       
    enddo
    
  end subroutine manplt_spr_NIsys_for_Ncell
  
  subroutine manplt_spr_NIsys(Sprs,NAEC_val,hmL0)
    implicit none
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0
    
    integer :: i,nsprsInACell
    integer :: sprL,sprR
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    do i = 1,(NAEC_val+1)
       sprL     = Sprs(i)          ; sprR = sprL+(Hlf_Ncell*nsprsInACell)
       l0(sprL) = (hmL0)*l0(sprL)  ; l0(sprR) = l0(sprL)
       write(*,*) sprL,l0(sprL),sprR,l0(sprR),"sprL-L0s"
    enddo
    
  end subroutine manplt_spr_NIsys
  

  subroutine manplt_spr_NIsys_withCrtclApSrfc(Sprs,NAEC_val,hmL0)
    implicit none
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0
    
    integer :: i,nsprsInACellT1,nsprsInACellT2
    integer :: sprL,sprR
    
    nsprsInACellT1 = (NAEC_Apcl     +NAEC_Bsal+NAEC_Ltrl+3)
    nsprsInACellT2 = (NAEC_ApclCrtcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    write(*,*) Sprs(1:(NAEC_val+1)),"hh"
    
    do i = 1,(NAEC_val+1)
       sprL     = Sprs(i)
       sprR     = sprL + (Hlf_Ncell-NCP_CrtclApSrfc)*(nsprsInACellT1) + (NCP_CrtclApSrfc*nsprsInACellT2)
       l0(sprL) = (hmL0)*l0(sprL)  ; l0(sprR) = l0(sprL)
       write(*,*) sprL,l0(sprL),sprR,l0(sprR),"sprL-L0s"
    enddo
    
  end subroutine manplt_spr_NIsys_withCrtclApSrfc
  
  subroutine manplt_spr_of_singl_cell_NIsys(Sprs,NAEC_val,hmL0)
    implicit none
    integer, intent(in) :: Sprs(1:(NAEC_val+1))
    integer, intent(in) :: NAEC_val
    real*8 , intent(in) :: hmL0
    
    integer :: i,nsprsInACell
    integer :: sprL
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl+3)
    
    do i = 1,(NAEC_val+1)
       sprL     = Sprs(i)         
       l0(sprL) = (hmL0)*l0(sprL)
       write(*,*) sprL,l0(sprL),"sprL-L0s"
    enddo
    
  end subroutine manplt_spr_of_singl_cell_NIsys
  
  
  subroutine str_or_rstore_l0_bfr_manplt(Sprs,NAEC_val,l0ValStr,strORrstore)
    implicit none
    integer, intent(in)    :: Sprs(1:(NAEC_val+1))
    integer, intent(in)    :: NAEC_val
    real*8,  intent(inout) :: l0ValStr(1:(NAEC_val+1))
    integer, intent(in)    :: strORrstore
    integer                :: j
    
    do j = 1,(NAEC_val+1)
       if (strORrstore==1) l0ValStr(j) = l0(Sprs(j))
       if (strORrstore==2) l0(Sprs(j)) = l0ValStr(j)
    enddo
  end subroutine str_or_rstore_l0_bfr_manplt
  
  subroutine get_Ratio_of_Strght_Curved_len_Of_LtrlSide(LtrlSp,Ratio_SC)
    implicit none
    integer, intent(in)  :: LtrlSp(1:(NAEC_Ltrl+1))
    real*8 , intent(out) :: Ratio_SC
    
    integer :: sprNm,node1,node2
    real*8  :: N1(1:N_dmnsn),N2(1:N_dmnsn)
    real*8  :: curvdLen,strghtLen
    integer :: RglrNode
    
    interface
       real*8 function Length(N1,N2,N_dmnsn)
         implicit none
         integer, intent(in) :: N_dmnsn
         real*8 , intent(in) :: N1(1:N_dmnsn),N2(1:N_dmnsn)
       end function Length
    end interface
    
    if (CellsMeet==0)   RglrNode=((Hlf_Ncell+1)*2)*2
    if (CellsMeet.gt.0) RglrNode=((Hlf_Ncell+1)*2)*2 + 1
    
    curvdLen = 0.0d0 ; strghtLen = 0.0d0
    
    do i = 1,(NAEC_Ltrl+1)
       
       sprNm = LtrlSp(i) 
       node1 = spr_node(sprNm,1) ; node2 = spr_node(sprNm,2)
       
       N1(1:N_dmnsn) = node_xy(node1,1:N_dmnsn)
       N2(1:N_dmnsn) = node_xy(node2,1:N_dmnsn)
       
       curvdLen = curvdLen + Length(N1,N2,N_dmnsn)
       write(*,*) curvdLen,i,"curvdLen"
       
    enddo
    
    node1 = spr_node(LtrlSp(1),1) ; node2 = spr_node(LtrlSp(NAEC_Ltrl+1),2)
    if ((node1.gt.RglrNode) .or. (node2.gt.RglrNode)) then
       write(*,*) 'node1 or node2 cant be greater than Rglr Node for strghtLen Calc'
       stop
    endif
    
    N1(1:N_dmnsn) = node_xy(node1,1:N_dmnsn)
    N2(1:N_dmnsn) = node_xy(node2,1:N_dmnsn)
    
    strghtLen = Length(N1,N2,N_dmnsn)
    write(*,*) strghtLen,"strghtLen"
    
    Ratio_SC = strghtLen/curvdLen 
    write(*,*) Ratio_SC,"Ratio_SC"
    
  end subroutine get_Ratio_of_Strght_Curved_len_Of_LtrlSide
  
  
  subroutine ApBsLt_DrawBoardMeasuredLengths()
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! ALTHOUGH FROM DRAWBOARD IT WAS 8 Cells, We are TAKING
    ! 11 Cells, keeping the FIRST 3 same as the 4th ONE
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    
    implicit none
    
    allocate(LenAp_Exp(1:Hlf_Ncell))
    allocate(LenBs_Exp(1:Hlf_Ncell))
    allocate(LenLt_Exp(1:Hlf_Ncell))
    
    LenAp_Exp(1) = 1.57d0 ; LenBs_Exp(1)  = 1.41d0 ; LenLt_Exp(1)  = 5.20d0
    LenAp_Exp(2) = 1.57d0 ; LenBs_Exp(2)  = 1.41d0 ; LenLt_Exp(2)  = 5.20d0
    LenAp_Exp(3) = 1.57d0 ; LenBs_Exp(3)  = 1.41d0 ; LenLt_Exp(3)  = 5.20d0
    LenAp_Exp(4) = 1.57d0 ; LenBs_Exp(4)  = 1.41d0 ; LenLt_Exp(4)  = 5.20d0
    LenAp_Exp(5) = 2.22d0 ; LenBs_Exp(5)  = 1.18d0 ; LenLt_Exp(5)  = 5.55d0
    LenAp_Exp(6) = 1.89d0 ; LenBs_Exp(6)  = 0.23d0 ; LenLt_Exp(6)  = 6.30d0
    LenAp_Exp(7) = 2.62d0 ; LenBs_Exp(7)  = 0.15d0 ; LenLt_Exp(7)  = 5.65d0
    LenAp_Exp(8) = 2.73d0 ; LenBs_Exp(8)  = 0.13d0 ; LenLt_Exp(8)  = 4.05d0
    LenAp_Exp(9) = 1.97d0 ; LenBs_Exp(9)  = 0.57d0 ; LenLt_Exp(9)  = 3.79d0
    LenAp_Exp(10)= 1.98d0 ; LenBs_Exp(10) = 1.81d0 ; LenLt_Exp(10) = 3.60d0
    LenAp_Exp(11)= 0.57d0 ; LenBs_Exp(11) = 3.87d0 ; LenLt_Exp(11) = 3.23d0
    
  end subroutine ApBsLt_DrawBoardMeasuredLengths
  
  subroutine ApBsLt_LengthRatios_BasedOnFirstCell()
    implicit none
    integer :: i
    
    open(unit=131,file='RL_val.dat')
    
    allocate(RL_Ap(1:Hlf_Ncell),RL_Bs(1:Hlf_Ncell),RL_Lt(1:Hlf_Ncell))
    
    do i = 1,Hlf_Ncell
       RL_Ap(i) = LenAp_Exp(i) / LenAp_Exp(1)
       RL_Bs(i) = LenBs_Exp(i) / LenBs_Exp(1)
       RL_Lt(i) = LenLt_Exp(i) / LenLt_Exp(1)
       write(131,*) RL_Ap(i),RL_Bs(i),RL_Lt(i),i,"RL"
    enddo
    
    close(131)
    
  end subroutine ApBsLt_LengthRatios_BasedOnFirstCell
  
  
  subroutine match_ratios_ApBsLt_BasedOnEXP_Images(cellNm,sprT)
    implicit none
    integer, intent(in)  :: cellNm,sprT
    integer, allocatable :: Sprs(:)
    
    integer :: nsprsInACell,imax
    real*8  :: CR,AR ! = CURRENT and ACTUAL RATIO
    real*8  :: hmL0
    real*8  :: closeNessRatio,extRatio
    integer :: cntTheEndless
    
    nsprsInACell = (NAEC_Apcl+1)+(NAEC_Bsal+1)+(NAEC_Ltrl+1)
    
    if (sprT==1) then
       allocate(Sprs(1:(NAEC_Apcl+1)))
       imax = NAEC_Apcl+1
       AR   = RL_Ap(cellNm)
    elseif (sprT==2) then
       allocate(Sprs(1:(NAEC_Bsal+1)))
       imax = NAEC_Bsal+1
       AR   = RL_Bs(cellNm)
    elseif (sprT==3) then
       allocate(Sprs(1:(NAEC_Ltrl+1)))
       imax = NAEC_Ltrl+1
       AR   = RL_Lt(cellNm)
    endif
    
    call get_the_NIsprs(cellNm,sprT,imax,Sprs)
    extRatio = 0.05d0 
    
    cntTheEndless = 1
    
    do 
       
       call get_curr_Ratio(Sprs,imax,sprT,CR)
       
       if (CR.le.AR) then ! lengths got smaller
          hmL0 = 1.05d0
          call manplt_spr_NIsys(Sprs,imax-1,hmL0)
          
       elseif (CR.gt.AR) then ! lengths got bigger
          hmL0 = 0.95d0
          call manplt_spr_NIsys(Sprs,imax-1,hmL0)
       endif
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       
       closeNessRATIO = (abs(CR-AR)/AR)
       cntTheEndless  = cntTheEndless+1
       write(*,*) cntTheEndLess,"cntTheEndLess"
       
       if (closeNessRATIO.lt.extRatio) exit
       if (cntTheEndless.gt.50) then
          write(*,*) cntTheEndless,"EXITING FOR cntTheEndless=51"
          exit
       endif
       
    enddo
    
  end subroutine match_ratios_ApBsLt_BasedOnEXP_Images
  
  subroutine get_the_NIsprs(cellNm,sprT,NumOFSgmnt,Sprs)
    implicit none
    integer, intent(in)  :: cellNm,sprT
    integer, intent(in)  :: NumOFSgmnt
    integer, intent(out) :: Sprs(1:NumOFSgmnt)
    
    integer :: nsprsInACell,i
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl)+3
    
    do i = 1,NumOFSgmnt
       if (sprT==1) Sprs(i) = (cellNm-1)*(nsprsInACell) + i
       if (sprT==2) Sprs(i) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + i
       if (sprT==3) Sprs(i) = (cellNm-1)*(nsprsInACell) + (NAEC_Apcl+1) + (NAEC_Bsal+1) +i
    enddo
    
    write(*,*) Sprs(1:NumOFSgmnt),"Sprs"
    
  end subroutine get_the_NIsprs
  
  
  subroutine get_curr_Ratio(Sprs,NumOFsgmnt,sprT,CR)
    implicit none
    integer, intent(in)  :: Sprs(1:NumOFsgmnt)
    integer, intent(in)  :: NumOFsgmnt
    integer, intent(in)  :: sprT
    real*8 , intent(out) :: CR
    
    integer :: SprsREF(1:NumOFSgmnt)
    integer :: REFcell
    real*8  :: REF_L,spr_L
    
    REFcell = 3                      ! Refernce cell has been set to 3 HERE, as cell 4 has been disrupted
    REF_L   = 0.0d0 ; spr_L = 0.0d0
    
    call get_the_NIsprs(REFcell,sprT,NumOFSgmnt,SprsREF)
    
    do i = 1,NumOFsgmnt
       REF_L = REF_L+l(SprsREF(i))
       Spr_L = Spr_L+l(Sprs(i))
    enddo
    
    CR = Spr_L/REF_L
    write(*,*) Spr_L,REF_L,CR,"Spr_L,REF_L,CR"
    
  end subroutine get_curr_Ratio
  
  subroutine shortn_Bsal_ofEndCell ! very very generalized routine
    implicit none
    integer :: nsprsInACell,i,cnt
    integer :: Sprs(1:(NAEC_Ltrl+1))
    real*8  :: initLen,currLen
    real*8  :: Ratio,hmL0
    real*8  :: RatioForExit
    
    if (modelID==1) then
       write(*,*) "Not for modelID=1",modelID
       stop
    endif
    
    nsprsInACell = (NAEC_Apcl+NAEC_Bsal+NAEC_Ltrl) + 3
    initLen      = 0.0d0
    currLen      = 0.0d0
    RatioForExit = 1.40183d0 ! COMPLETELY BASED ON IMAGE LENGTH
    
    do i = 1,(NAEC_Ltrl+1)
       Sprs(i) = (N_cell-1)*(nsprsInACell) + i
       initLen = initLen + l(SprS(i))
    enddo
    
    cnt = 1
    
    do
       hmL0 = 0.95d0
       call manplt_spr_of_singl_cell_NIsys(Sprs,NAEC_Ltrl,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       
       do i = 1,(NAEC_Ltrl+1)
          currLen = currLen + l(SprS(i))
       enddo
       
       Ratio   = initLen/currLen
       write(*,*) Ratio,currLen,initLen,cnt,Frame_NI-1,"params in loop"
       currLen = 0.0d0
       
       if (Ratio.gt.RatioForExit) then
          write(*,*) Ratio,"Ratio at time of EXIT"
          exit
       endif
       
       cnt = cnt+1
    enddo
    
  end subroutine shortn_Bsal_ofEndCell
  
  subroutine diagonal_tnsn_rls_test
    implicit none
    integer :: cellL,cntLp
    real*8  :: hmL0Ap,hmL0Bs,hmL0Lt
    integer :: crtclNode
    
    integer, allocatable :: ApclSp(:),BsalSp(:),LtrlSp(:)
    
    integer, allocatable :: ApclSp1(:),ApclSp2(:),ApclSp3(:)
    integer, allocatable :: BsalSp1(:),BsalSp2(:),BsalSp3(:)
    integer, allocatable :: LtrlSp1(:),LtrlSp2(:),LtrlSp3(:)
    
    if (modelID==1) call switch_to_NI_and_Equilibrate
    !if (modelID==1) call switch_to_NI_and_not_Equilibrate
    write(*,*) modelID,"modelID in diagonal_tnsn_rls_test"
    write(*,*) Frame_NI-1,"Frame_NI in diagonal_tnsn_rls_test"
    
    allocate(ApclSp(1:(NAEC_Apcl+1)),BsalSp(1:(NAEC_Bsal+1)),LtrlSp(1:(NAEC_Ltrl+1)))
    
    allocate(ApclSp1(1:(NAEC_Apcl+1)),ApclSp2(1:(NAEC_Apcl+1)),ApclSp3(1:(NAEC_Apcl+1)))
    allocate(BsalSp1(1:(NAEC_Bsal+1)),BsalSp2(1:(NAEC_Bsal+1)),BsalSp3(1:(NAEC_Bsal+1)))
    allocate(LtrlSp1(1:(NAEC_Ltrl+1)),LtrlSp2(1:(NAEC_Ltrl+1)),LtrlSp3(1:(NAEC_Ltrl+1)))
    
    call find_the_cell_with_diagonal_spr(cellL)
    
    cntLp=1
    
    do 
       
       if (cntLp.le.10) then
          
          if (cntLp==1) then
             call get_ApclBsalLtrl_OfNIsystm(cellL,ApclSp,BsalSp,LtrlSp)
             write(*,*) cellL,ApclSp,BsalSp,LtrlSp,"CellL and Sprs"
          endif
          
          hmL0Lt = 0.95d0
          call manplt_spr_NIsys(LtrlSp,NAEC_Ltrl,hmL0Lt)
          
       elseif (cntLp.gt.10) then
          
          if (cntLp==11) then
             
             call get_ApclBsalLtrl_OfNIsystm((cellL+1),ApclSp1,BsalSp1,LtrlSp1)
             call get_ApclBsalLtrl_OfNIsystm((cellL+2),ApclSp2,BsalSp2,LtrlSp2)
             call get_ApclBsalLtrl_OfNIsystm((cellL+3),ApclSp3,BsalSp3,LtrlSp3)
             
             write(*,*) cellL+1,ApclSp1,BsalSp1,LtrlSp1,"CellL and Sprs 1"
             write(*,*) cellL+2,ApclSp2,BsalSp2,LtrlSp2,"CellL and Sprs 2"
             write(*,*) cellL+3,ApclSp3,BsalSp3,LtrlSp3,"CellL and Sprs 3"
             
          endif
          
          hmL0Ap = 0.95d0
          call manplt_spr_NIsys(ApclSp1,NAEC_Apcl,hmL0Ap)
          call manplt_spr_NIsys(ApclSp2,NAEC_Apcl,hmL0Ap)
          call manplt_spr_NIsys(ApclSp3,NAEC_Apcl,hmL0Ap)
          
       endif
       
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       
       crtclNode = 2*cellL+1
       write(*,*) node_xy(crtclNode,1:N_dmnsn),"crtclNode"
       
       if (cntLp.gt.10) then
          if (abs(node_xy(crtclNode,1)) .le. 0.05d0) then
             write(*,*) cntLp,"cntLp Val at exit"
             exit
          endif
       endif
       
       cntLp = cntLp+1
       
    enddo
    
  end subroutine diagonal_tnsn_rls_test
  
  
  subroutine find_the_cell_with_diagonal_spr(cellL)
    implicit none
    integer, intent(out) :: cellL
    
    integer :: i,j
    integer :: nodeInCell_curr,nodeInCell_prv,nodeInCell_max
    integer :: nodeInCell
    
    nodeInCell_prv  = 0
    nodeInCell_curr = 0
    nodeInCell_max  = 0
    
    do i = 1,Hlf_Ncell
       
       nodeInCell_curr = area_node(i,0)
       write(*,*) nodeInCell_curr,i,"currNode"
       
       if (i.gt.1) nodeInCell_prv = area_node(i-1,0)
       write(*,*) nodeInCell_prv,i,"prvNode"
       
       if (nodeInCell_curr.gt.nodeInCell_prv) then
          nodeInCell_max = nodeInCell_curr
       else
          continue
       endif
       
    enddo
    
    
    do i = 1,Hlf_Ncell
       nodeInCell = area_node(i,0)
       
       if (nodeInCell==nodeInCell_max) then
          cellL = i-1 ! gives the cell before that pulley point holdng cell
          write(*,*) "cellNum = ",cellL
          exit
       endif
    enddo
    
    write(*,*) cellL,"cellL in find"
    
  end subroutine find_the_cell_with_diagonal_spr
  
  subroutine boundary_force_test
    implicit none
    character(len=200)   :: flplce
    character(len=100)   :: flnm
    integer              :: ExpNo,FrmNoToBeRead
    integer              :: NumCells,sprTyp
    integer, allocatable :: cellNums(:)
    real*8               :: hmL0
    integer              :: cntLp,FrmNoIncr
    integer              :: apcl_or_bsal
    
    if (modelID==1) call switch_to_NI_and_Equilibrate
    write(*,*) modelID,"modelID boundary_force_test"
    
    write(flplce,*)&
         "/home/rniloy/PROJECT_CELL/2D_codes/unit_blck/restructuring_Data_Structure/NonUniformPress_PullingCasesDoubleMembrneTnsn/"
    
    ExpNo = 34 ; FrmNoToBeRead = 6
    call read_config_and_start_simlnFrm_there_flpceLoc(flplce,ExpNo,FrmNoToBeRead)
    call Equilibrate_system
    call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
    Frame_NI = Frame_NI+1
    
    write(*,*) "stopping after two force case"
    
    apcl_or_bsal = 3
    
    if (apcl_or_bsal==1) call manplts_needed_for_fixing_CS_aft_mkingForce0_at_apcl_bndry(apcl_or_bsal)
    if (apcl_or_bsal==2) call manplts_needed_for_fixing_CS_aft_mkingForce0_at_bsal_bndry(apcl_or_bsal)
    if (apcl_or_bsal==3) call manplts_needed_for_fixing_CS_having_Force_at_both_bndry(apcl_or_bsal)
    
  end subroutine boundary_force_test
  
  subroutine release_Xforce_at_the_apcl_or_basal_boundary(apcl_or_bsal)
    implicit none
    integer, intent(in) :: apcl_or_bsal
    
    real*8  :: releasePerItr
    integer :: Apcl_bound(1:2),Bsal_bound(1:2),i,imax
    
    call get_the_boundary_node(Apcl_bound,Bsal_bound)
    
    releasePerItr = 0.050d0
    imax = int(CgXNode(Apcl_bound(1))*100)/5
    write(*,*) imax,"imax inside release_Xforce_at_the_apcl_or_bsal_boundary"
    
    do i = 1,imax
       
       if (apcl_or_bsal == 1) then
          CgXNode(Apcl_bound(1)) = CgXNode(Apcl_bound(1)) - releasePerItr 
          CgXNode(Apcl_bound(2)) = -CgXNode(Apcl_bound(1))

          if (i==imax) CgXNode(Apcl_bound(1:2)) = 0.00d0
          
       elseif (apcl_or_bsal == 2) then
          CgXNode(Bsal_bound(1)) = CgXNode(Bsal_bound(1)) - releasePerItr 
          CgXNode(Bsal_bound(2)) = -CgXNode(Bsal_bound(1))

          if (i==imax) CgXNode(Bsal_bound(1:2)) = 0.00d0
       endif
       
       
       write(*,*) CgXNode(Bsal_bound(1)),CgXNode(Bsal_bound(2)),"CgXNode at the Bsal bound"
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
       
    enddo
    
  end subroutine release_Xforce_at_the_apcl_or_basal_boundary
  
  subroutine get_the_boundary_node(Apcl_bound,Bsal_bound)
    implicit none
    integer, intent(out) :: Apcl_bound(1:2),Bsal_bound(1:2)
    
    Apcl_bound(1) = 1 ; Apcl_bound(2) = (Hlf_Ncell+1)*2+1
    Bsal_bound(1) = 2 ; Bsal_bound(2) = (Hlf_Ncell+1)*2+2
    
  end subroutine get_the_boundary_node
  
  subroutine manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
    implicit none
    integer, intent(in)  :: NumCells
    integer, intent(in)  :: cellNums(1:NumCells)
    integer, intent(in)  :: sprTyp
    real*8 , intent(in)  :: hmL0
    
    integer, allocatable :: ApclSp(:,:),BsalSp(:,:),LtrlSp(:,:)
    
    integer              :: i,j,cellV
    integer, allocatable :: ApclSpV(:),BsalSpV(:),LtrlSpV(:)
    
    if (modelID==1) call switch_to_NI_and_Equilibrate
    write(*,*) modelID,"modelID in manplt_multiple_cells_NIsys"
    
    allocate(ApclSp(1:NumCells,1:(NAEC_Apcl+1)),BsalSp(1:NumCells,1:(NAEC_Bsal+1)),LtrlSp(1:NumCells,1:(NAEC_Ltrl+1)))
    allocate(ApclSpV(1:(NAEC_Apcl+1)),BsalSpV(1:(NAEC_Bsal+1)),LtrlSpV(1:(NAEC_Ltrl+1)))
    
    do i = 1,NumCells
       
       cellV = cellNums(i)
       call get_ApclBsalLtrl_OfNIsystm(cellV,ApclSpV,BsalSpV,LtrlSpV)
       
       ApclSp(i,1:(NAEC_Apcl+1)) = ApclSpV(1:(NAEC_Apcl+1))
       BsalSp(i,1:(NAEC_Bsal+1)) = BsalSpV(1:(NAEC_Bsal+1))
       LtrlSp(i,1:(NAEC_Ltrl+1)) = LtrlSpV(1:(NAEC_Ltrl+1))
       
       write(*,*) ApclSp(i,1:(NAEC_Apcl+1)),"for Apcl where i =",i,"and cellV =",cellV
       write(*,*) BsalSp(i,1:(NAEC_Bsal+1)),"for Bsal where i =",i,"and cellV =",cellV
       write(*,*) LtrlSp(i,1:(NAEC_Ltrl+1)),"for Ltrl where i =",i,"and cellV =",cellV   
       
       if (sprTyp == 1) call manplt_spr_NIsys(ApclSpV,NAEC_Apcl,hmL0)
       if (sprTyp == 2) call manplt_spr_NIsys(BsalSpV,NAEC_Bsal,hmL0)
       if (sprTyp == 3) call manplt_spr_NIsys(LtrlSpV,NAEC_Ltrl,hmL0)
       
    enddo
    
  end subroutine manplt_multiple_cells_spr_NIsys
  
  
  subroutine manplts_needed_for_fixing_CS_aft_mkingForce0_at_bsal_bndry(apcl_or_bsal)
    implicit none
    integer, intent(in)  :: apcl_or_bsal
    integer              :: NumCells,sprTyp
    integer, allocatable :: cellNums(:)
    real*8               :: hmL0
    integer              :: cntLp,FrmNoIncr
    
    if (apcl_or_bsal .ne. 2) then
       write(*,*) "apcl_or_bsal val is wrong in Force0_at_bsal_bndry",apcl_or_bsal
       stop
    endif
    
    call release_Xforce_at_the_apcl_or_basal_boundary(apcl_or_bsal)
    
    NumCells=2 ; allocate(cellNums(1:NumCells))
    cellNums(1)=5 ; cellNums(2)=6 ; sprTyp=2 ; hmL0=1.05d0
    
    do cntLp = 1,10
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    NumCells=2 ; allocate(cellNums(1:NumCells))
    cellNums(1)=4 ; cellNums(2)=7 ; sprTyp=2 ; hmL0=1.10d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    NumCells=2 ; allocate(cellNums(1:NumCells))
    cellNums(1)=2 ; cellNums(2)=3 ; sprTyp=2 ; hmL0=1.10d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    NumCells=3 ; allocate(cellNums(1:NumCells))
    cellNums(1)=8 ; cellNums(2)=9 ; cellNums(3)=10 ; sprTyp=1 ; hmL0=0.90d0
    
    do cntLp = 1,5
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    
  end subroutine manplts_needed_for_fixing_CS_aft_mkingForce0_at_bsal_bndry
  
  
  subroutine manplts_needed_for_fixing_CS_aft_mkingForce0_at_apcl_bndry(apcl_or_bsal)
    implicit none
    integer, intent(in)  :: apcl_or_bsal
    integer              :: NumCells,sprTyp
    integer, allocatable :: cellNums(:)
    real*8               :: hmL0
    integer              :: cntLp,FrmNoIncr
    
    if (apcl_or_bsal .ne. 1) then
       write(*,*) "apcl_or_bsal val is wrong in Force0_at_bsal_bndry",apcl_or_bsal
       stop
    endif
    
    NumCells=5 ; allocate(cellNums(1:NumCells))
    cellNums(1)=2 ; cellNums(2)=3 ; cellNums(3)=4 ; cellNums(4)=5 ; cellNums(5)=6
    sprTyp=2 ; hmL0=0.90d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    !call release_Xforce_at_the_apcl_or_basal_boundary(apcl_or_bsal)
    
    NumCells=3 ; allocate(cellNums(1:NumCells))
    cellNums(1)=8 ; cellNums(2)=9 ; cellNums(3)=10 ; sprTyp=1 ; hmL0=1.10d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    call release_Xforce_at_the_apcl_or_basal_boundary(apcl_or_bsal)
    
    NumCells=4 ; allocate(cellNums(1:NumCells))
    cellNums(1)=3 ; cellNums(2)=4 ; cellNums(3)=5 ; cellNums(4)=6 ; sprTyp=2 ; hmL0=0.90d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    NumCells=1 ; allocate(cellNums(1:NumCells))
    cellNums(1)=7 ; sprTyp=1 ; hmL0=1.05d0
    
    do cntLp = 1,3
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    NumCells=1 ; allocate(cellNums(1:NumCells))
    cellNums(1)=1 ; sprTyp=2 ; hmL0=0.90d0
    
    do cntLp = 1,4
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
  end subroutine manplts_needed_for_fixing_CS_aft_mkingForce0_at_apcl_bndry
  
  
  subroutine manplts_needed_for_fixing_CS_having_Force_at_both_bndry(apcl_or_bsal)
    implicit none
    integer, intent(in)  :: apcl_or_bsal
    integer              :: NumCells,sprTyp
    integer, allocatable :: cellNums(:)
    real*8               :: hmL0
    integer              :: cntLp,FrmNoIncr
    
    if (apcl_or_bsal .ne. 3) then
       write(*,*) "apcl_or_bsal val is wrong in Force_at_both_bndry",apcl_or_bsal
       stop
    endif
    
    NumCells=3 ; allocate(cellNums(1:NumCells))
    cellNums(1)=8 ; cellNums(2)=9 ; cellNums(3)=10 ; sprTyp=1 ; hmL0=1.05d0
    
    do cntLp = 1,4
       call manplt_multiple_cells_spr_NIsys(NumCells,cellNums,sprTyp,hmL0)
       call Equilibrate_system
       call save_config_and_generate_colorData(coordntes_xy,l0,A0,Exprmnt_NI,Frame_NI)
       Frame_NI = Frame_NI+1
    enddo
    deallocate(cellNums)
    
    
  end subroutine manplts_needed_for_fixing_CS_having_Force_at_both_bndry
  
  
end module nodeInsrtd_cntrlStates

































